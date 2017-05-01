
abstract class Face extends Transformable {
	'constructor': { new (surface: Surface, contour: Edge[], holes?: Edge[][], name?: string, info?: any): Face }
	allEdges: Edge[]
	protected aabb: AABB
	constructor(readonly surface: Surface,
	            readonly contour: Edge[],
	            readonly holes: Edge[][] = [],
	            readonly name?: string,
				readonly info?: any) {
		super()
		//assert(name)
		Edge.assertLoop(contour)
		assert(contour.every(f => f instanceof Edge), 'contour.every(f => f instanceof Edge)' + contour.toSource())
		// contour.forEach(e => !surface.containsCurve(e.curve) &&
		// console.log('FAIL:'+surface.distanceToPoint(e.curve.anchor)))
		contour.forEach(e => assert(surface.containsCurve(e.curve), 'edge not in surface ' + e + surface))
		assert(surface.edgeLoopCCW(contour), surface.toString()+contour.join('\n'))
		holes && holes.forEach(hole => Edge.assertLoop(hole))
		holes && holes.forEach(hole => assert(!surface.edgeLoopCCW(hole)))
		assert(!holes || holes.constructor == Array, holes && holes.toString())
		this.allEdges = Array.prototype.concat.apply(this.contour, this.holes)
	}

	static assembleFacesFromLoops(loops: Edge[][], surface: Surface, faceConstructor): Face[] {
		type LoopInfo = {loop: Edge[], ccw: boolean, subloops: LoopInfo[]}
		function placeRecursively(newLoopInfo: LoopInfo, loopInfos: LoopInfo[]) {
			if (loopInfos.length == 0) {
				loopInfos.push(newLoopInfo)
			} else {
				const subLoopInfo = loopInfos.find(loopInfo => B2.loop1ContainsLoop2(loopInfo.loop, newLoopInfo.loop, surface))
				if (subLoopInfo) {
					placeRecursively(newLoopInfo, subLoopInfo.subloops)
				} else {
					// newLoopInfo isnt contained by any other subLoopInfo
					for (let i = loopInfos.length; --i >= 0;) {
						const subLoopInfo = loopInfos[i]
						//console.log('cheving subLoopInfo', surface.loopContainsPoint(newLoopInfo.edges, subLoopInfo.edges[0].a))
						if (B2.loop1ContainsLoop2(newLoopInfo.loop, subLoopInfo.loop, surface)) {
							newLoopInfo.subloops.push(subLoopInfo)
							loopInfos.splice(i, 1) // remove it
						}
					}
					loopInfos.push(newLoopInfo)
				}
			}
		}

		function newFacesRecursive(loopInfo: LoopInfo): void {
			newFaces.push(new faceConstructor(surface,
				loopInfo.ccw ? loopInfo.loop : Edge.reverseLoop(loopInfo.loop),
				loopInfo.subloops.map(sl => sl.ccw ? Edge.reverseLoop(sl.loop) : sl.loop)))
			loopInfo.subloops.forEach(sl => sl.subloops.forEach(sl2 => newFacesRecursive(sl2)))
		}

		const newFaces: Face[] = []
		const topLevelLoops: LoopInfo[] = []
		loops.forEach(loop => placeRecursively({loop: loop, ccw: surface.edgeLoopCCW(loop), subloops: []}, topLevelLoops))
		topLevelLoops.forEach(tll => newFacesRecursive(tll))
		return newFaces
	}

	fromLoops(loops: Edge[][], surface) {
		type LoopInfo = {loop: Edge[], ccw: boolean, subloops: LoopInfo[]}
		function placeRecursively(newLoopInfo: LoopInfo, loopInfos: LoopInfo[]) {
			if (loopInfos.length == 0) {
				loopInfos.push(newLoopInfo)
			} else {
				const subLoopInfo = loopInfos.find(loopInfo => B2.loop1ContainsLoop2(loopInfo.loop, newLoopInfo.loop, surface))
				if (subLoopInfo) {
					placeRecursively(newLoopInfo, subLoopInfo.subloops)
				} else {
					// newLoopInfo isnt contained by any other subLoopInfo
					for (let i = loopInfos.length; --i >= 0;) {
						const subLoopInfo = loopInfos[i]
						//console.log('cheving subLoopInfo', surface.loopContainsPoint(newLoopInfo.edges, subLoopInfo.edges[0].a))
						if (B2.loop1ContainsLoop2(newLoopInfo.loop, subLoopInfo.loop, surface)) {
							newLoopInfo.subloops.push(subLoopInfo)
							loopInfos.splice(i, 1) // remove it
						}
					}
					loopInfos.push(newLoopInfo)
				}
			}
		}

		function newFacesRecursive(loopInfo: LoopInfo): void {
			// CW loops can be top level, if they are holes in the original face not contained in the new face
			if (loopInfo.ccw) {
				if (loopInfo.subloops.every(sl => !sl.ccw)) {
					const newFace = new faceConstructor(surface, loopInfo.loop, loopInfo.subloops.map(sl => sl.loop))
					newFaces.push(newFace)
					loopInfo.subloops.forEach(sl => sl.subloops.forEach(slsl => slsl.ccw && newFacesRecursive(slsl)))
				} else {
					loopInfo.subloops.forEach(sl => sl.ccw && newFacesRecursive(sl))
				}
			}
		}

		const newFaces: Face[] = []
		const topLevelLoops:LoopInfo[] = []
		loops.forEach(loop => placeRecursively({loop: loop, ccw: surface.edgeLoopCCW(loop), subloops: []}, topLevelLoops))
		topLevelLoops.forEach(tll => newFacesRecursive(tll))
		return newFaces
	}


	abstract intersectFace(face2: Face,
	                       thisBrep: B2,
	                       face2Brep: B2,
	                       faceMap: Map<Face, Edge[]>,
	                       thisEdgePoints: NLA.CustomMap<Edge, IntersectionPointInfo[]>,
	                       otherEdgePoints: NLA.CustomMap<Edge, IntersectionPointInfo[]>,
	                       likeSurfaceFaces: Set<string>): void

	transform(m4: M4): this {
		const newEdges = this.contour.map(e => e.transform(m4))
		const newHoles = this.holes.map(hole => hole.map(e => e.transform(m4)))
		return new this.constructor(this.surface.transform(m4), newEdges, newHoles, this.name, this.info) as this
	}


	flipped() {
		const newEdges = this.contour.map(e => e.flipped()).reverse()
		const newHoles = this.holes.map(hole => hole.map(e => e.flipped()).reverse())
		return new this.constructor(this.surface.flipped(), newEdges, newHoles, this.name, this.info)
	}

	toString() {
		return 'new ' + this.constructor.name + '(' + this.surface + ', [' + this.contour.map(e => '\n\t' + e).join() + ']'
			+ this.holes.map(hole => '\n\t\thole: ' + hole.join()) + ')'
	}

	toSource() {
		return `new ${this.constructor.name}(${this.surface.toSource()}, [${this.contour.map(e => '\n\t' + e.toSource()).join(',')}], [${
			this.holes.map(hole => '[' + hole.map(e => '\n\t' + e.toSource()).join(',') + ']').join(',')}])`
	}

	equals(obj: any): boolean {
		function loopsEqual(a, b) {
			return a.length == b.length &&
				NLA.arrayRange(0, a.length, 1)
					.some(offset => a.every((edge, i) => edge.equals(b[(offset + i) % a.length])))

		}

		return this == obj ||
			Object.getPrototypeOf(this) == Object.getPrototypeOf(obj)
			&& this.holes.length == obj.holes.length
			&& loopsEqual(this.contour, obj.contour)
			&& this.holes.every(hole => obj.holes.some(hole2 => loopsEqual(hole, hole2)))
	}

	hashCode() {
		function arrayHashCode(array) {
			let hashCode = 0
			for (const val of array) {
				hashCode = hashCode * 31 + val | 0
			}
			return hashCode
		}
		function loopHashCode(loop) { return arrayHashCode(loop.map(edge => edge.hashCode()).sort()) }
		let hashCode = 0
		hashCode = hashCode * 31 + arrayHashCode(this.holes.map(loop => loopHashCode(loop)).sort()) | 0
		hashCode = hashCode * 31 + loopHashCode(this.contour) | 0
		hashCode = hashCode * 31 + this.surface.hashCode() | 0
		return hashCode
	}

	likeFace(face2) {
		function loopsLike(a, b) {
			return a.length == b.length &&
				NLA.arrayRange(0, a.length, 1)
					.some(offset => a.every((edge, i) => edge.like(b[(offset + i) % a.length])))

		}
		assertInst(Face, face2)
		return this.surface.like(face2.surface)
			&& this.holes.length == face2.holes.length
			&& loopsLike(this.contour, face2.contour)
			&& this.holes.every(hole => face2.holes.some(hole2 => loopsLike(hole, hole2)))
	}

	getAllEdges(): Edge[] {
		return this.allEdges
	}

	addEdgeLines(mesh) {
		assert(false, 'buggy, fix')
		const vertices = this.contour.map(edge => edge.getVerticesNo0()).concatenated(), mvl = mesh.vertices.length
		for (let i = 0; i < vertices.length; i++) {
			mesh.vertices.push(vertices[i])
			mesh.lines.push(mvl + i, mvl + (i + 1) % vertices.length)

		}
	}

	containsPoint(p: V3): boolean {
		assertVectors(p)
		return this.surface.loopContainsPoint(this.contour, p) != PointVsFace.OUTSIDE
			&& !this.holes.some(hole => this.surface.loopContainsPoint(hole, p) != PointVsFace.OUTSIDE)
	}

	containsPoint2(p: V3): PointVsFace {
		assertVectors(p)
		const contourContainsPoint = this.surface.loopContainsPoint(this.contour, p)
		if (contourContainsPoint != PointVsFace.INSIDE) return contourContainsPoint
		for (const hole of this.holes) {
			const loopContainsPoint = this.surface.loopContainsPoint(hole, p)
			if (loopContainsPoint != PointVsFace.OUTSIDE) {
				return loopContainsPoint == PointVsFace.ON_EDGE ? PointVsFace.ON_EDGE : PointVsFace.OUTSIDE
			}
		}
		return PointVsFace.INSIDE
	}

	/**
	 *
	 * @param line
	 * @returns t param of the line if there is an intersection, NaN otherwise
	 */
	intersectsLine(line: L3): number {
		assertInst(L3, line)
		if (!this.getAABB().intersectsLine(line)) return NaN
		const containedIntersectionsTs = this.surface.isTsForLine(line).filter(t => this.containsPoint(line.at(t)))
		const nearestPointT = containedIntersectionsTs.withMax(t => -t)

		return undefined != nearestPointT ? nearestPointT : NaN
	}

	toMesh() {
		const mesh = new GL.Mesh({triangles: true, normals: true, lines: true})
		this.addToMesh(mesh)
		//mesh.compile()
		return mesh
	}

	abstract addToMesh(mesh: GL.Mesh)

	zDirVolume(): {centroid: V3, volume: number} {
		return this.surface.zDirVolume(this.getAllEdges())
	}

	calcArea(): number {
		return this.surface.calculateArea(this.getAllEdges())
	}

	getLoops(): Edge[][] {
		return this.holes.concat(this.contour)
	}

	getAABB(): AABB {
		return this.aabb || (this.aabb = AABB.forAABBs(this.contour.map(e => e.getAABB())))
	}

	static create(surface: Surface, faceEdges: Edge[], holes?: Edge[][], faceName?: string, info?: any) {
		return surface instanceof PlaneSurface
			? new PlaneFace(surface, faceEdges, holes, faceName, info)
			: new RotationFace(surface, faceEdges, holes, faceName, info)
	}

	pointsToInside3(p: V3, curve: Curve, curveT: number, dir: -1 | 1): PointVsFace {
		const eps = 1e-6
		const normal = this.surface.normalAt(p)
		const curveTangent = curve.tangentAt(curveT).times(dir)
		const up = normal.cross(curveTangent)
		const ecd = curve.at(curveT).to(curve.at(curveT + dir * eps)).dot(up)
		let minValue = Infinity, result, advanced = false
		for (const edge of this.getAllEdges()) {
			const aEqP = edge.a.like(p), bEqP = edge.b.like(p)
			assert(aEqP == edge.a.like(p))
			assert(bEqP == edge.b.like(p))
			if (!aEqP && !bEqP) continue
			const edgeTangent = aEqP ? edge.aDir : edge.bDir.negated()
			const angle = curveTangent.angleRelativeNormal(edgeTangent, normal)
			if (eq0(angle)) {
				if (curve.isColinearTo(edge.curve)) {
					return PointVsFace.ON_EDGE
				}
				const edgeT = aEqP ? edge.aT : edge.bT
				const edgeDir = (aEqP ? 1 : -1) * sign(edge.deltaT())
				const iscd = edge.curve.diff(edgeT, edgeDir * eps).dot(up)
				//const iscd = edge.curve.at(edgeT).to(curve.at(edgeT + edgeDir * eps)).dot(up)
				const diff = iscd - ecd
				if (diff > 0 && (!advanced || diff < minValue)) {
					advanced = true
					minValue = diff
					result = aEqP ? PointVsFace.OUTSIDE : PointVsFace.INSIDE
				}
			} else if (!advanced) {
				const angle2 = (angle + TAU) % TAU
				if (angle2 < minValue) {
					minValue = angle2
					result = aEqP ? PointVsFace.OUTSIDE : PointVsFace.INSIDE
				}
			}
		}
		return result
	}

	pointsToInside2(p: V3, dir: V3): PointVsFace {
		return this.pointsToInside3(p, L3.anchorDirection(p, dir), 0, 1)
		const normal = this.surface.normalAt(p)
		let minAngle = Infinity, inOut = false
		function test(v, b) {
			const angle = (dir.angleRelativeNormal(v, normal) + TAU + NLA_PRECISION / 2) % TAU
			if (angle <= 2 * NLA_PRECISION) {
				return true
			}
			if (angle < minAngle) {
				minAngle = angle
				inOut = b
			}
		}
		for (const edge of this.getAllEdges()) {
			assert(edge.a.equals(p) || !edge.a.like(p))
			assert(edge.b.equals(p) || !edge.b.like(p))
			if (edge.a.equals(p) && test(edge.aDir, false)) return PointVsFace.ON_EDGE
			if (edge.b.equals(p) && test(edge.bDir.negated(), true)) return PointVsFace.ON_EDGE
		}
		return inOut ? PointVsFace.INSIDE : PointVsFace.OUTSIDE
	}
}

class PlaneFace extends Face {

	surface: PlaneSurface

	constructor(p: P3 | PlaneSurface, contour: Edge[], holes?: Edge[][], name?: string, info?: any) {
		assert(p instanceof P3 || p instanceof PlaneSurface)
		super(p instanceof P3 ? new PlaneSurface(p) : p, contour, holes, name, info)
	}


	zDirVolume(): {centroid: V3, volume: number} {
		let {centroid, area} = this.calculateArea()
		return {volume: this.surface.plane.normal1.z * centroid.z * area,
			centroid: new V3(centroid.x, centroid.y, centroid.z / 2) }

	}

	calculateArea(): {centroid: V3, area: number} {
		let centroid = V3.O, tcs = 0, tct = 0, totalArea = 0
		let r1 = this.surface.right, u1 = this.surface.up
		this.contour.forEach(edge => {
			let edgeCentroid, edgeArea: number, centroidS, centroidT
			if (edge instanceof StraightEdge) {
				const midPoint = edge.a.lerp(edge.b, 0.5)
				edgeCentroid = new V3(midPoint.x, centroid.y, centroid.z / 2)
				centroidS = midPoint.dot(r1) / 2
				centroidT = midPoint.dot(u1)
				const edgeLength = edge.a.distanceTo(edge.b)
				edgeArea = edgeLength * edge.curve.dir1.dot(r1)
				edgeArea = (edge.a.dot(u1) + edge.b.dot(u1)) / 2 * edge.b.to(edge.a).dot(r1)
			} else {
				let curve = edge.curve
				if (curve instanceof SemiEllipseCurve) {
					let info = curve.getAreaInDir(r1, u1, edge.aT, edge.bT)
					edgeArea = info.area
					let parametricCentroid = this.surface.pointToParameterFunction()(info.centroid)
					centroidS = parametricCentroid.x
					centroidT = parametricCentroid.y
				} else if (curve instanceof BezierCurve) {
					edgeArea = curve.getAreaInDirSurface(u1, this.surface, edge.aT, edge.bT).area
				} else {
					assertNever()
				}
			}


			tcs += edgeArea * centroidS
			tct += edgeArea * centroidT
			totalArea += edgeArea
		})
		centroid = r1.times(tcs).plus(u1.times(tct))
		assert(isFinite(totalArea))
		return {area: totalArea, centroid: centroid}
	}

	addToMesh(mesh) {
		const mvl = mesh.vertices.length
		const normal = this.surface.plane.normal1
		const vertices = this.contour.flatMap(edge => edge.getVerticesNo0())
		for (let i = 0; i < vertices.length; i++) { mesh.lines.push(mvl + i, mvl + (i + 1) % vertices.length) }
		const holeStarts = []
		this.holes.forEach(hole => {
			holeStarts.push(vertices.length)
			vertices.pushAll(hole.flatMap(edge => edge.getVerticesNo0()))
		})
		const triangles = triangulateVertices(normal, vertices, holeStarts).map(index => index + mvl)
		Array.prototype.push.apply(mesh.vertices, vertices)
		Array.prototype.push.apply(mesh.triangles, triangles)
		Array.prototype.push.apply(mesh.normals, NLA.arrayFromFunction(vertices.length, () => normal))
	}

	intersectsLine(line): number {
		assertInst(L3, line)
		const lambda = line.intersectWithPlaneLambda(this.surface.plane)
		if (!Number.isFinite(lambda)) {
			return NaN
		}
		const inside = this.containsPoint(line.at(lambda))
		return inside ? lambda : NaN
	}

	withHole(holeEdges) {
		return new PlaneFace(this.surface, this.contour, [holeEdges])
	}

	intersectFace(face2, thisBrep, face2Brep, faceMap, thisEdgePoints, otherEdgePoints, checkedPairs) {
		//const f = face2 instanceof PlaneFace ? this.intersectPlaneFace : RotationFace.intersectFace
		//f.call(this, face2, thisBrep, face2Brep, faceMap, thisEdgePoints, otherEdgePoints, checkedPairs)
		//if (face2 instanceof PlaneFace) {
		//	this.intersectPlaneFace(face2, thisBrep, face2Brep, faceMap, thisEdgePoints, otherEdgePoints, checkedPairs)
		//	return
		//}
		RotationFace.prototype.intersectFace.call(
			this, face2, thisBrep, face2Brep, faceMap, thisEdgePoints, otherEdgePoints, checkedPairs)
		//RotationFace.prototype.intersectFace.call(
		//	face2, this, face2Brep, thisBrep, faceMap, otherEdgePoints, thisEdgePoints, checkedPairs)
		return
		assert(false)
	}

	intersectPlaneFace(face2: PlaneFace,
	                   thisBrep: B2,
	                   face2Brep: B2,
	                   faceMap: Map<Face, Edge[]>,
	                   thisEdgePoints: NLA.CustomMap<Edge, { edge: Edge, edgeT: number, p: V3, passEdge?: Edge }[]>,
	                   otherEdgePoints: NLA.CustomMap<Edge, { edge: Edge, edgeT: number, p: V3, passEdge?: Edge }[]>,
	                   checkedPairs: NLA.CustomSet<NLA.Pair<Equalable, Equalable>>) {
		assertInst(NLA.CustomMap, thisEdgePoints, otherEdgePoints)

		function hasPair(a: NLA.Equalable, b: NLA.Equalable) {
			return checkedPairs.has(new NLA.Pair(a, b))
		}
		function addPair(a: NLA.Equalable, b: NLA.Equalable) {
			return checkedPairs.add(new NLA.Pair(a, b))
		}

		/**
		 * @param newEdge generated segment
		 * @param col1 if newEdge is colinear to an edge of this, the edge in question
		 * @param col2 same for face2
		 */
		function handleNewEdge(newEdge: StraightEdge, col1: Edge, col2: Edge) {
			if (!col1 && !col2) {
				NLA.mapPush(faceMap, face, newEdge)
				NLA.mapPush(faceMap, face2, newEdge.flipped())
				return true
			}
			function handleEdgeInFace(col1, col2, face, face2, thisBrep, face2Brep, coplanarSameIsInside: boolean, has, add) {
				if (col1 && !col2) {
					if (hasPair(col1.getCanon(), face2)) return

					//add(col1.getCanon(), face2)
					const face2Plane = face2.surface.plane

					// NB: a new edge is inserted even though it may be the same as an old one
					// however it indicates that it intersects the other volume here, i.e. the old edge cannot
					// be counted as 'inside' for purposes of reconstitution
					thisBrep.edgeFaces.get(col1.getCanon()).forEach(faceInfo => {
						//const dot = NLA.snap0(face2Plane.normal1.dot(faceInfo.inside))
						//if (dot == 0 ? !coplanarSameIsInside : dot < 0) {
						const pointsInsideFace = fff(faceInfo, face2.surface)
						const edgeInside = pointsInsideFace == INSIDE || !coplanarSameIsInside && pointsInsideFace == COPLANAR_SAME
						const pushEdge = (faceInfo.edge.aDir.like(newEdge.aDir)) ? newEdge : newEdge.flipped()
						assert(faceInfo.edge.aDir.like(pushEdge.aDir))
						edgeInside && NLA.mapPush(faceMap, faceInfo.face, pushEdge)
					})

					const newEdgeInside = face2Plane.normal1.cross(newEdge.aDir)
					const sVEF1 = splitsVolumeEnclosingFaces(thisBrep, col1.getCanon(), newEdgeInside, face2Plane.normal1)
					let addNewEdge, addNewEdgeFlipped
					if (addNewEdge = sVEF1 == INSIDE || coplanarSameIsInside && sVEF1 == COPLANAR_SAME) {
						NLA.mapPush(faceMap, face2, newEdge)
					}
					const sVEF2 = splitsVolumeEnclosingFaces(thisBrep, col1.getCanon(), newEdgeInside.negated(), face2Plane.normal1)
					if (addNewEdgeFlipped = sVEF2 == INSIDE || coplanarSameIsInside && sVEF2 == COPLANAR_SAME) {
						NLA.mapPush(faceMap, face2, newEdge.flipped())
					}
					if (addNewEdge || addNewEdgeFlipped || sVEF1 == COPLANAR_SAME && sVEF2 == INSIDE || sVEF2 == COPLANAR_SAME && sVEF1 == INSIDE) {
						return true
					}
				}
			}
			const c1 = handleEdgeInFace(col1, col2, face, face2, thisBrep, face2Brep, false, hasPair, addPair)
			const c2 = handleEdgeInFace(col2, col1, face2, face, face2Brep, thisBrep, true, (a, b) => hasPair(b, a), (a, b) => addPair(b, a))
			if (c1 || c2) return true

			if (col1 && col2) {
				if (hasPair(col1.getCanon(), col2.getCanon())) return

				addPair(col1.getCanon(), col2.getCanon())

				function handleColinearEdgeFaces(col1, col2, thisBrep, face2Brep, coplanarSameIsInside: boolean, thisEdgePoints, has, add) {
					// not entirely sure for what i had the dirInsides in?
					//const aDirNegatedInside = (newEdge.a.like(col2.a) || newEdge.a.like(col2.b)) && splitsVolumeEnclosingCone(face2Brep, newEdge.a, newEdge.aDir.negated()) == INSIDE
					//const bDirInside = (newEdge.b.like(col2.a) || newEdge.b.like(col2.b)) && splitsVolumeEnclosingCone(face2Brep, newEdge.b, newEdge.bDir) == INSIDE
					thisBrep.edgeFaces.get(col1.getCanon()).forEach(faceInfo => {
						const sVEF = splitsVolumeEnclosingFaces(face2Brep, col2.getCanon(), faceInfo.inside, faceInfo.normalAtCanonA)
						const edgeInside = sVEF == INSIDE || coplanarSameIsInside && sVEF == COPLANAR_SAME
						const pushEdge = (faceInfo.edge.aDir.like(newEdge.aDir)) ? newEdge : newEdge.flipped()
						edgeInside && NLA.mapPush(faceMap, faceInfo.face, pushEdge)
					})
				}
				handleColinearEdgeFaces(col1, col2, thisBrep, face2Brep, true, thisEdgePoints, hasPair, addPair)
				handleColinearEdgeFaces(col2, col1, face2Brep, thisBrep, false, otherEdgePoints, (a, b) => hasPair(b, a), (a, b) => addPair(b, a))
			}
		}


		// what needs to be generated: new edges on face
		// points on edges where they are cut by faces so that sub edges will be generated for loops
		// points on ends of edges where the edge will be an edge in the new volume where it goes from A to B
		//         you don't want thos to be marked as 'inside', otherwise invalid faces will be added
		// if a face cuts a corner, nothings needs to be done, as that alone does not limit what adjacent faces will be
		function handleEndPoint(a: IntersectionPointInfo, b: IntersectionPointInfo, newEdge: Edge) {
			// ends in the middle of b's face
			if (a && !b) {
				if (!a.colinear && a.edgeT != a.edge.aT && a.edgeT != a.edge.bT) {
					NLA.mapPush(thisEdgePoints, a.edge.getCanon(), a)
					assert(a.edge.isValidT(a.edgeT))
				}
				// else colinear segment ends in middle of other face, do nothing
			}
			// ends in the middle of a's face
			if (b && !a) {
				if (!b.colinear && b.edgeT != b.edge.aT && b.edgeT != b.edge.bT) {
					NLA.mapPush(otherEdgePoints, b.edge.getCanon(), b)
					assert(b.edge.isValidT(b.edgeT))
				}
				// else colinear segment ends in middle of other face, do nothing
			}
			if (a && b) {
				// if a or b is colinear the correct points will already have been added to the edge by handleNewEdge
				// segment starts/ends on edge/edge intersection
				function foo(a, b, face, face2, thisPlane, face2Plane, thisBrep, face2Brep, first, thisEdgePoints) {
					if (!a.colinear && a.edgeT != a.edge.aT && a.edgeT != a.edge.bT) {
						if (!hasPair(a.edge.getCanon(), b.edge.getCanon())) {
							addPair(a.edge.getCanon(), b.edge.getCanon())
							// ends on a, on colinear segment b bT != a.edge.bT &&
							// b can be colinear, so edgeT == aT is possible
							if (a.p.like(b.edge.a) || a.p.like(b.edge.b)) {
								const corner = a.p.like(b.edge.a) ? b.edge.a : b.edge.b
								// face2brep corner on edge
								const sVEC1 = splitsVolumeEnclosingCone(face2Brep, corner, a.edge.aDir)
								const sVEC2 = splitsVolumeEnclosingCone(face2Brep, corner, a.edge.aDir.negated())
								// if either of these return ALONG_EDGE_OR_PLANE, then the breps share a colinear edge

								if (INSIDE == sVEC1 || INSIDE == sVEC2) {
									NLA.mapPush(thisEdgePoints, a.edge.getCanon(), a)
									assert(a.edge.isValidT(a.edgeT))
								}
							} else {
								// edge / edge center intersection
								const aEdgeDir = a.edge.tangentAt(a.edgeT)
								const bEdgeDir = b.edge.tangentAt(b.edgeT)
								const testVector = aEdgeDir.rejectedFrom(bEdgeDir)
								assert(!testVector.isZero())
								const sVEF1 = splitsVolumeEnclosingFaces(face2Brep, b.edge.getCanon(), testVector, thisPlane.normal1)
								const sVEF2 = splitsVolumeEnclosingFaces(face2Brep, b.edge.getCanon(), testVector.negated(), thisPlane.normal1)
								if (INSIDE == sVEF1 || INSIDE == sVEF2) {
									NLA.mapPush(thisEdgePoints, a.edge.getCanon(), a)
									assert(a.edge.isValidT(a.edgeT))
								}
							}
						}
					}
				}

				foo(a, b, face, face2, thisPlane, face2Plane, thisBrep, face2Brep, true, thisEdgePoints)
				foo(b, a, face2, face, face2Plane, thisPlane, face2Brep, thisBrep, false, otherEdgePoints)

			}
		}


		assertInst(PlaneFace, face2)
		const face: PlaneFace = this
		// get intersection
		const thisPlane = this.surface.plane, face2Plane = face2.surface.plane
		if (thisPlane.isParallelToPlane(face2Plane)) {
			if (thisPlane.like(face2Plane)) {
				// normal1 same and same location in space
				// addLikeSurfaceFaces(likeSurfaceFaces, this, face2)
			}
			return
		}
		const isLine = L3.fromPlanes(thisPlane, face2Plane)
		// get intersections of newCurve with other edges of face and face2
		const ps1 = planeFaceEdgeISPsWithPlane(face, isLine, face2Plane)
		const ps2 = planeFaceEdgeISPsWithPlane(face2, isLine, thisPlane)
		if (ps1.length == 0 || ps2.length == 0) {
			// faces to not intersect
			return
		}

		let col1: IntersectionPointInfo, col2: IntersectionPointInfo
		let in1 = false, in2 = false
		let i = 0, j = 0, last
		let startP, startDir, startT, startA, startB
		while (i < ps1.length || j < ps2.length) {
			assert(i <= ps1.length)
			assert(j <= ps2.length)
			const a = ps1[i], b = ps2[j]
			assert(a || b)
			if (j == ps2.length || i < ps1.length && NLA.lt(a.t, b.t)) {
				last = a
				in1 = !in1
				a.used = true
				in1 && (col1 = a.colinear && a)
				i++
			} else if (i == ps1.length || NLA.gt(a.t, b.t)) {
				last = b
				in2 = !in2
				b.used = true
				in2 && (col2 = b.colinear && b)
				j++
			} else {
				// TODO: this will break if 3 points on the same t
				last = a
				in1 = !in1
				in2 = !in2
				//if (in1 == in2) {
				a.used = true
				b.used = true
				in1 && (col1 = a.colinear && a)
				in2 && (col2 = b.colinear && b)
				//}
				i++
				j++
			}
			if (startP && !(in1 && in2)) {
				// segment end
				const newEdge = new StraightEdge(isLine, startP, last.p, startT, last.t, null, 'genseg' + globalId++)
				startP = undefined
				last.used = true
				if (handleNewEdge(newEdge, col1 && col1.edge, col2 && col2.edge)) {
					handleEndPoint(startA || col1, startB || col2, newEdge)
					handleEndPoint(a && a.used && a || col1, b && b.used && b || col2, newEdge)
				}
			} else if (in1 && in2) {
				// new segment just started
				startP = last.p
				startDir = last.insideDir
				startT = last.t
				startA = a && a.used && a
				startB = b && b.used && b
			}
			if (!in1 && a && last == a && a.colinear) {
				checkedPairs.add(new NLA.Pair(a.edge.getCanon(), face2))
			}
			if (!in2 && b && (last == b || b.used) && b.colinear) {
				checkedPairs.add(new NLA.Pair(b.edge.getCanon(), face))
			}
		}
	}


	static forVertices(planeSurface, vs, ...holeVss): PlaneFace {
		if (planeSurface instanceof P3) {
			planeSurface = new PlaneSurface(planeSurface)
		}
		assert(isCCW(vs, planeSurface.plane.normal1), 'isCCW(vs, planeSurface.plane.normal1)')
		const edges = StraightEdge.chain(vs)
		holeVss.forEach(vs => assert(doubleSignedArea(vs, planeSurface.plane.normal1) >= 0, 'doubleSignedArea(vs, planeSurface.plane.normal1) >= 0'))
		const holes = holeVss.map(hvs => StraightEdge.chain(hvs))
		return new PlaneFace(planeSurface, edges, holes)
	}

	pointsToInside(p: V3, dir: V3): PointVsFace {
		return this.containsPoint2(p.plus(dir.times(NLA_PRECISION * 8)))
	}
}
NLA.registerClass(PlaneFace)



class RotationFace extends Face {
	constructor(rot: Surface, contour: Edge[], holes?: Edge[][], name?: string, info?: any) {
		super(rot, contour, holes, name, info)
	}

	static loopDoesNotCrossPlane(loop: Edge[], seamPlane: P3) {
		let side = 0
		// returns true if d is on the other side as previous calls
		function checkSide(d) {
			if (side == 0) {
				side = d
			} else {
				return !side || side * d < 0
			}
		}
		for (const edge of loop) {
			const ts = edge.edgeISTsWithPlane(seamPlane)
			if (ts.length == 0) {
				if (!(edge.curve instanceof L3) && checkSide(seamPlane.distanceToPointSigned(edge.a))) return false
			} else {
				for (const t of ts) {
					// TODO: this part probably should be in a separate function
					// check 'backwards' only if if aT != t
					if (edge.aT != t) {
						if (checkSide(dotCurve2(edge.curve, t, seamPlane.normal1, -sign(edge.bT - edge.aT)))) return false
					}
					if (edge.bT != t) {
						if (checkSide(dotCurve2(edge.curve, t, seamPlane.normal1, sign(edge.bT - edge.aT)))) return false
					}
				}
			}
		}
		return true
	}

	getAABB() {
		if (this.aabb) return this.aabb
		if (this.surface instanceof SemiEllipsoidSurface || this.surface instanceof EllipsoidSurface) {
			const aabb = AABB.forAABBs(this.contour.map(e => e.getAABB()))
			aabb.addPoints(this.surface.getExtremePoints().filter(p => this.containsPoint(p)))
			return aabb
		} else {
			return super.getAABB()
		}
	}

	getCanonSeamU(): number {
		const pointToParameterFunction = this.surface.pointToParameterFunction()
		for (const edge of this.contour) {
			// check edge.a
			let u = pointToParameterFunction(edge.a, PI).x
			// if u is not PI, or ~0, return its sign
			if (u != PI && !NLA.eq0(u)) {
				return sign(u) * PI
			}
			// check midpoint between edge.a and edge.b
			u = pointToParameterFunction(edge.curve.at((edge.aT + edge.bT) / 2), PI).x
			if (u != PI && !NLA.eq0(u)) {
				return sign(u) * PI
			}
		}
		const localEdge = this.contour[0].transform(this.surface.inverseMatrix)
		if (P3.ZX.containsCurve(localEdge.curve)) {
			const insideVector = localEdge.a.cross(localEdge.aDir)
			return sign(insideVector.dot(V3.Y)) * PI
		}
		assert(false, 'Couldn\'t find canon seam u')
	}


	unrollLoop(edgeLoop) {
		const vs = []
		const reverseFunc = this.surface.pointToParameterFunction()
		const verticesNo0s = edgeLoop.map(edge => edge.getVerticesNo0())
		const startEdgeIndex = verticesNo0s.findIndex(edgeVertices => !NLA.eq(reverseFunc(edgeVertices[0], Math.PI).x, Math.PI))
		assert(-1 != startEdgeIndex)
		// console.log(startEdgeIndex)
		let hint = Math.PI
		for (let i = 0; i < edgeLoop.length; i++) {
			let edgeIndex = (i + startEdgeIndex) % edgeLoop.length
			for (let j = 0; j < verticesNo0s[edgeIndex].length; j++) {
				let p = verticesNo0s[edgeIndex][j]
				let localP = reverseFunc(p, hint)
				if (Math.abs(localP.x) < Math.PI - NLA_PRECISION) {
					// update hint
					hint = localP.x
				}
				// console.log(hint, p.sce, localP.sce)
				vs.push(localP)
			}
		}
		// edgeLoop.forEach((edge, e) => {
		// 	var hint = edge.bDir
		// 	if (edge instanceof StraightEdge && edge.curve.dir1.isParallelTo(this.surface.dir || this.surface.dir1)) {
		// 		hint = this.surface.normalAt(edge.b).cross(edge.bDir)
		// 	}
		// 	edge.getVerticesNo0().forEach(p => {
		// 		vs.push(reverseFunc(p, hint))
		// 	})
		// })
		// console.log('vs\n', vs.join('\n'), vs.length)
		return vs
	}

	/**
	 * f1 cos t + f2 sin t
	 * tan(phi) = sin / cos
	 *          = (f1x cos t + f2x sin t) / (f1y cos t + f2y sin t)
	 *
	 *          = (-f1x sin t + f2x cos t) / (-f1y sin t + f2y cos t)
	 */
	unrollEllipsoidLoops(edgeLoops: Edge[][], uStep: number, vStep: number) {
		const verticesUV = [], vertices = [], normals = [], loopStarts = []
		const ellipsoid: SemiEllipsoidSurface = this.surface as SemiEllipsoidSurface
		const ptpf = ellipsoid.pointToParameterFunction()
		for (const edgeLoop of edgeLoops) {
			loopStarts.push(verticesUV.length)
			// console.log(startEdgeIndex)
			const hint = this.getCanonSeamU()
			for (let i = 0; i < edgeLoop.length; i++) {
				const ipp = (i + 1) % edgeLoop.length
				const verticesNo0 = edgeLoop[i].getVerticesNo0()
				vertices.pushAll(verticesNo0)
				normals.pushAll(verticesNo0.map(v => ellipsoid.normalAt(v)))
				verticesUV.pushAll(verticesNo0.map(v => { const uv = ptpf(v, hint); return new V3(uv.x / uStep, uv.y / vStep, 0) }))
				const nextStart = edgeLoop[ipp].a
				//console.log('BLAH', nextStart.str, ellipsoid.center.plus(ellipsoid.f3).str)
				if (nextStart.like(ellipsoid.center.plus(ellipsoid.f3)) || nextStart.like(ellipsoid.center.minus(ellipsoid.f3))) {
					console.log('FIXING')
					const localbDir = ellipsoid.inverseMatrix.transformVector(edgeLoop[i].bDir), localaDir = ellipsoid.inverseMatrix.transformVector(edgeLoop[ipp].aDir)
					let inAngle = Math.atan2(-localbDir.y, -localbDir.x)
					if (abs(inAngle) > Math.PI - NLA_PRECISION) {
						assert(hint == -PI || hint == PI)
						inAngle = hint
					}
					let outAngle = Math.atan2(localaDir.y, localaDir.x)
					if (abs(outAngle) > Math.PI - NLA_PRECISION) {
						assert(hint == -PI || hint == PI)
						outAngle = hint
					}

					const uvLast = verticesUV.pop()
					verticesUV.push(new V3(inAngle / uStep, uvLast.y, 0), new V3(outAngle / uStep, uvLast.y, 0))
					vertices.push(vertices.last())
					normals.push(normals.last())
				}
				verticesUV.forEach(({x: u, y: v}) => {
					assert(isFinite(u))
					assert(isFinite(v))
				})
			}
		}
		assert(vertices.length == vertices.length)
		//console.log(verticesUV.map(v => v.str).join('\n'))
		return {verticesUV: verticesUV, vertices: vertices, normals: normals, loopStarts: loopStarts}
	}

	unrollCylinderLoops(loops, uStep, vStep) {
		const vertexLoops = loops.map(loop => loop.map(edge => edge.getVerticesNo0()).concatenated())
		const vertices: V3[] = vertexLoops.concatenated()
		const normals: V3[] = vertices.map(v => this.surface.normalAt(v))
		// this.unrollLoop(loop).map(v => new V3(v.x / uStep, v.y / vStep, 0)))
		const loopStarts = vertexLoops.reduce((arr, loop) => (arr.push(arr.last() + loop.length), arr), [0])
		const pointToParameterFunction = this.surface.pointToParameterFunction()
		const verticesUV = vertices.map(v => { const uv = pointToParameterFunction(v); return new V3(uv.x / uStep, uv.y / vStep, 0) })
		return {verticesUV: verticesUV, vertices: vertices, normals: normals, loopStarts: loopStarts}
	}

	/**
	 * at(s, t) = new V3(s cos t, s sin t, t + )
	 *
	 * x = 0
	 *
	 * s cos t = 0
	 * ==> s = 0 || cos t = 0
	 * ==> L3.Z || V3(0, +-s, k * 2 pi)
	 *
	 * x = c
	 * s cos t = c
	 * ==> V3(c, c sin t / cos t = c tan t, t)
	 * ==> V3(c, c t, arctan t)
	 *
	 *
	 * x . n = w
	 *      s cos t nx + s sin t ny + t nz = w
	 *      s = (w - t nz) / (cos t nx + sub t ny)
	 * ==> V3(
	 *          cos t (w - t nz) / (cos t nx + sin t ny)
	 *          sin t (w - t nz) / (cos t nx + sin t ny)
	 *          t)
	 *
	 *  ==> V3(
	 *          (w - z arctan t) / (x + t y)
	 *          (w - z arctan t) / (y + x / t)
	 *          arctan t)
	 *
	 *
	 *
	 */

	addToMesh(mesh: GL.Mesh, uStep: number, vStep: number) {
		vStep = vStep || this.surface.vStep
		uStep = uStep || this.surface.uStep
		assertf(() => uStep > 0 && vStep > 0, uStep, vStep, 'Surface: ' + this.surface)
		const triangles = []
		const f = (i, j) => this.surface.parametricFunction()(i * uStep, j * vStep)
		const normalF = (i, j) => this.surface.parametricNormal()(i * uStep, j * vStep)
		const loops = [this.contour].concat(this.holes)
		const {vertices, verticesUV, normals, loopStarts} = this.surface instanceof SemiEllipsoidSurface
			? this.unrollEllipsoidLoops(loops, uStep, vStep)
			: this.unrollCylinderLoops(loops, uStep, vStep)
		loopStarts.push(vertices.length)

		for (let vertexLoopIndex = 0; vertexLoopIndex < loops.length; vertexLoopIndex++) {
			const vertexLoopStart = loopStarts[vertexLoopIndex]
			const vertexLoopLength = loopStarts[vertexLoopIndex + 1] - vertexLoopStart
			const base = mesh.vertices.length + loopStarts[vertexLoopIndex]
			for (let i = 0; i < vertexLoopLength; i++) {
				mesh.lines.push(base + i, base + (i + 1) % vertexLoopLength)
			}
		}

		disableConsole()
		let minU = Infinity, maxU = -Infinity, minV = Infinity, maxV = -Infinity
		//console.log('surface', this.surface.str)
		//console.log(verticesUV)
		//drPs.pushAll(verticesUV.map((v, i) => ({p: vertices[i], text: `${i} uv: ${v.toString(x => NLA.round10(x, -4))}`})))
		verticesUV.forEach(({x: u, y: v}) => {
			assert(isFinite(u))
			assert(isFinite(v))
			minU = min(minU, u)
			maxU = max(maxU, u)
			minV = min(minV, v)
			maxV = max(maxV, v)
		})
		const uOffset = floor(minU + NLA_PRECISION), vOffset = floor(minV + NLA_PRECISION)
		const uRes = ceil(maxU - NLA_PRECISION) - uOffset, vRes = ceil(maxV - NLA_PRECISION) - vOffset
		console.log(uStep, vStep, uRes, vRes)
		if (uRes == 1 && vRes == 1) {
			// triangulate this face as if it were a plane
			const polyTriangles = triangulateVertices(V3.Z, verticesUV, loopStarts.slice(1, 1 + this.holes.length))
			triangles.pushAll(polyTriangles)
		} else {
			const partss: int[][][] = new Array(uRes * vRes)

			function fixUpPart(part, baseU, baseV) {
				assert(baseU < uRes && baseV < vRes, `${baseU}, ${baseV}, ${uRes}, ${vRes}`)
				console.log('complete part', part, baseU, baseV)
				//console.trace()
				assert(part.length)
				const cellU = baseU + uOffset, cellV = baseV + vOffset
				for (const index of part) {
					assert(NLA.le(cellU, verticesUV[index].x) && NLA.le(verticesUV[index].x, cellU + 1), `${index} ${verticesUV[index].str} ${cellU} ${cellU}`)
					assert(NLA.le(cellV, verticesUV[index].y) && NLA.le(verticesUV[index].y, cellV + 1))
				}
				const pos = baseV * uRes + baseU
				;(partss[pos] || (partss[pos] = [])).push(part)
				//const outline = partss[pos] || (partss[pos] = [minU + baseU * uStep, minV + baseV * vStep, minU + (baseU + 1) * uStep, minV + (baseV + 1) * vStep])
			}

			// 'some' instead of forEach so we can return out of the entire function if this.edges crosses no borders and
			for (let vertexLoopIndex = 0; vertexLoopIndex < loops.length; vertexLoopIndex++) {
				let part: int[], firstPart, firstPartBaseU, firstPartBaseV
				let lastBaseV = -1, lastBaseU = -1
				let partCount = 0
				const vertexLoopStart = loopStarts[vertexLoopIndex]
				const vertexLoopLength = loopStarts[vertexLoopIndex + 1] - vertexLoopStart
				for (let vlvi = 0; vlvi < vertexLoopLength; vlvi++) {
					const vx0index = vertexLoopStart + vlvi, vx0 = verticesUV[vx0index]
					const vx1index = vertexLoopStart + (vlvi + 1) % vertexLoopLength, vx1 = verticesUV[vx1index]
					//console.log('dask', vx0index, vx1index)
					const vx01 = vx0.to(vx1)
					assert(vx0)
					const di = vx01.x, dj = vx01.y
					let vxIndex = vx0index, vx = vx0, currentT = 0
					let whileLimit = 400
					while (--whileLimit) {
						const vxu = vx.x, vxv = vx.y
						// points which are on a grid line are assigned to the cell into which they are going (+ NLA_PRECISION * sign(di))
						// if they are parallel to the gridline (eq0(di)), they belong the the cell for which they are a CCW boundary
						const baseU = floor(vxu + (!eq0(di) ? sign(di) : -sign(dj)) * NLA_PRECISION) - uOffset
						const baseV = floor(vxv + (!eq0(dj) ? sign(dj) : sign(di)) * NLA_PRECISION) - vOffset
						assert(baseU < uRes && baseV < vRes, `${baseU}, ${baseV}, ${uRes}, ${vRes}`)
						// figure out the next intersection with a gridline:
						// iNext is the positive horizontal distance to the next vertical gridline
						const iNext = ceil(sign(di) * vxu + NLA_PRECISION) - sign(di) * vxu
						const jNext = ceil(sign(dj) * vxv + NLA_PRECISION) - sign(dj) * vxv
						const iNextT = currentT + iNext / abs(di)
						const jNextT = currentT + jNext / abs(dj)
						//console.log(vxIndex, vx.str, 'vij', vxu, vxv, 'd', di, dj, 'ijNext', iNext, jNext, 'nextT', iNextT, jNextT)
						if (lastBaseU != baseU || lastBaseV != baseV) {
							if (part) {
								if (!firstPart) {
									firstPart = part
									firstPartBaseU = lastBaseU
									firstPartBaseV = lastBaseV
								} else {
									partCount++
									fixUpPart(part, lastBaseU, lastBaseV)
								}
							}
							part = [vxIndex]
						}
						lastBaseU = baseU
						lastBaseV = baseV
						currentT = min(iNextT, jNextT)
						if (NLA.ge(currentT, 1)) {
							//console.log('breaking ', vx1index)
							part.push(vx1index)
							break
						} else {
							const nextPoint = vx0.lerp(vx1, currentT)
							const nextPointIndex = addVertex(nextPoint.x, nextPoint.y)

							//console.log('pushing ', nextPointIndex)
							part.push(nextPointIndex)
							vx = nextPoint
							vxIndex = nextPointIndex
						}
					}
					assert(whileLimit, 'whileLimit')
				}
				if (0 == partCount) {
					// complete loop
					assert(false, 'found a hole, try increasing resolution')
				}
				// at this point, the firstPart hasn't been added, and the last part also hasn't been added
				// either they belong to the same cell, or not
				if (firstPartBaseU == lastBaseU && firstPartBaseV == lastBaseV) {
					part.pop()
					fixUpPart(part.concat(firstPart), lastBaseU, lastBaseV)
				} else {
					fixUpPart(firstPart, firstPartBaseU, firstPartBaseV)
					fixUpPart(part, lastBaseU, lastBaseV)
				}
				console.log('firstPart', firstPart)
			}
			console.log('calculated parts', partss)
			const fieldVertexIndices = new Array((uRes + 1) * (vRes + 1))

			function addVertex(u, v): int {
				verticesUV.push(new V3(u, v, 0))
				normals.push(normalF(u, v))
				return vertices.push(f(u, v)) - 1
			}

			function getGridVertexIndex(i, j): int {
				const index = j * (uRes + 1) + i
				return fieldVertexIndices[index] || (fieldVertexIndices[index] = addVertex(i + uOffset, j + vOffset))
			}

			for (let col = 0; col < uRes; col++) {
				let inside = false
				for (let row = 0; row < vRes; row++) {
					const pos = row * uRes + col
					const fieldU = uOffset + col, fieldV = vOffset + row
					const fieldCU = uOffset + col + 0.5, fieldCV = vOffset + row + 0.5
					const parts = partss[pos]
					if (!parts) {
						if (inside) {
							pushQuad(triangles, false,
								getGridVertexIndex(col, row), getGridVertexIndex(col + 1, row),
								getGridVertexIndex(col, row + 1), getGridVertexIndex(col + 1, row + 1))
						}
					} else {
						// assemble the field with segments in in
						function opos(index: int) {
							const p = verticesUV[index], u1 = p.x - fieldU, v1 = p.y - fieldV
							assert(-NLA_PRECISION < u1 && u1 < 1 + NLA_PRECISION && -NLA_PRECISION < v1 && v1 < 1 + NLA_PRECISION,
								'oob u1 v1 ' + u1 + ' ' + v1 + ' ' + index + ' ' + p.str + 'IF THIS FAILS check canonSeamU is correct')
							return v1 < u1 ? u1 + v1 : 4 - u1 - v1
						}

						while (parts.length) {
							const outline = [], outlineVertexIndices = []
							const startPart = parts[0]
							assert(startPart.length > 0)
							let currentPart = startPart
							do {
								outline.pushAll(currentPart)
								const currentPartEndOpos = opos(currentPart.last())
								const nextPartIndex = parts.indexWithMax(part => -NLA.mod(opos(part[0]) - currentPartEndOpos, 4))
								const nextPart = parts.removeIndex(nextPartIndex)
								let currentOpos = currentPartEndOpos
								const nextPartStartOpos = opos(nextPart[0]) > currentOpos ? opos(nextPart[0]) : opos(nextPart[0]) + 4
								let nextOpos = ceil(currentOpos + NLA_PRECISION)
								let flipping = eq0((currentOpos + NLA_PRECISION) % 1 - NLA_PRECISION)
								//inside = inside != (!eq0(currentOpos % 1) && currentOpos % 2 < 1)
								while (NLA.lt(nextOpos, nextPartStartOpos)) {
									switch (nextOpos % 4) {
										case 0:
											outline.push(getGridVertexIndex(col, row))
											break
										case 1:
											inside = inside != flipping
											outline.push(getGridVertexIndex(col + 1, row))
											break
										case 2:
											outline.push(getGridVertexIndex(col + 1, row + 1))
											break
										case 3:
											inside = inside != flipping
											outline.push(getGridVertexIndex(col, row + 1))
											break

									}
									flipping = true
									nextOpos++
								}
								// if the next loop would have completed a top or bottom segment
								inside = inside != (flipping && nextOpos % 2 == 1 && eq(nextOpos, nextPartStartOpos))
								currentOpos = nextOpos
								currentPart = nextPart
							} while (currentPart != startPart)
							// triangulate outline
							if (outline.length == 3) {
								// its just a triangle
								triangles.pushAll(outline)
							} else {
								const polyTriangles = triangulateVertices(V3.Z, outline.map(i => verticesUV[i]), []).map(i => outline[i])
								triangles.pushAll(polyTriangles)
							}
							//console.log('outline', col, row, outline)
						}
					}
				}
			}

		}
		//console.log('trinagle', triangles.max(), vertices.length, triangles.length, triangles.toSource(), triangles.map(col => vertices[col].$).toSource() )
		//assert(normals.every(n => n.hasLength(1)), normals.find(n => !n.hasLength(1)).length() +' '+normals.findIndex(n => !n.hasLength(1)))
		Array.prototype.push.apply(mesh.triangles, triangles.map(index => index + mesh.vertices.length))
		Array.prototype.push.apply(mesh.vertices, vertices)
		Array.prototype.push.apply(mesh.normals, normals)
		//this.addEdgeLines(mesh)
		enableConsole()
	}
	//addToMesh(mesh: GL.Mesh, uStep: number, vStep: number) {
	//    const closed = false
	//    let minU = Infinity, maxU = -Infinity, minV = Infinity, maxV = -Infinity
	//    const f = this.surface.parametricFunction()
	//    const normalF = this.surface.parametricNormal()
	//    const vertexLoops = this.holes.concat([this.edges]).map(loop => this.unrollLoop(loop))
	//    vertexLoops.forEach(vertexLoop => {
	//        vertexLoop.forEach(({x: u, y: v}) => {
	//            minU = min(minU, u)
	//            maxU = max(maxU, u)
	//            minV = min(minV, v)
	//            maxV = max(maxV, v)
	//        })
	//    })
	//const uRes = ceil(maxU / uStep) - floor(minU / uStep)
	//const vRes = ceil(maxV / vStep) - floor(minV / vStep)
	//   const fields = new Array(uRes * vRes)
	//const flippeds = {}
	//function addSegment(p0, p1, pos) {
	//	const arr = fields[pos] || (fields[pos] = [])
	//	arr.push(p0, p1)
	//}
	//function fixUpPart(part, row, col) {
	//	const pos = row * uRes + col
	//	;(fields[pos] || (fields[pos] = [])).push(part)
	//	//const outline = fields[pos] || (fields[pos] = [minU + row * uStep, minV + col * vStep, minU + (row + 1) * uStep, minV + (col + 1) * vStep])
	//}
	//   console.log('u', minU, maxU, 'v', minV, maxV, vertexLoops[0].toSource().replace(/\), /g, ',\n'))
	//    vertexLoops.forEach(vertexLoop => {
	//       let part, firstPart
	//       let lastCol = -1, lastRow = -1
	//       let blarhj = 0
	//        vertexLoop.forEach((v0, v0index, vs) => {
	//            const v1 = vs[(v0index + 1) % vs.length], v01 = v0.to(v1)
	//           const di = v01.x / uStep, dj = v01.y / vStep
	//		let v = v0, iNextT = 0, jNextT = 0
	//           while (true) {
	//            const vi = (v.x - minU) / uStep, vj = (v.y - minV) / vStep
	//            // figure out the next intersection with a gridline:
	//            const iNext = NLA.mod(-sign(di) * vi + 2 * NLA_PRECISION, 1)
	//            const jNext = NLA.mod(-sign(dj) * vj + 2 * NLA_PRECISION, 1)
	//            iNextT += iNext / di
	//            jNextT += jNext / dj
	//            const row = floor(vi + (!eq0(di) ? sign(di) : -sign(dj)) * NLA_PRECISION)
	//            const col = floor(vj + (!eq0(dj) ? sign(dj) : -sign(di)) + NLA_PRECISION)
	//            if (lastRow != row || lastCol != col) {
	//            	if (part) {
	//            	    if (!firstPart) {
	//            	    	firstPart = part
	//	                } else {
	//            	    	blarhj++
	//		                fixUpPart(part, lastRow, lastCol)
	//	                }
	//	            }
	//	            part = [v]
	//            }
	//            lastRow = row
	//            lastCol = col
	//            const nextT = min(iNextT, jNextT)
	//            if (NLA.ge(nextT, 1)) {
	//	            part.push(v1)
	//	            break
	//            } else {
	//	            const nextPoint = v0.lerp(v1, nextT)
	//
	//	            part.push(nextPoint)
	//	            v = nextPoint
	//	            part = []
	//            }
	//           }
	//        })
	//       fixUpPart(part.concat(firstPart), lastRow, lastCol)
	//    })
	//const fieldVertexIndices = new Array((uRes + 1) * (vRes + 1))
	//let vertices = [], triangles = [], normals = []
	//function getVertexIndex(i, j) {
	//    	const index = i * (uRes + 1) + j
	//	return fieldVertexIndices[index]
	//		|| (normals.push(normalF(i * uStep + minU, j * vStep + minV)),
	//			fieldVertexIndices[index] = vertices.push(f(i * uStep + minU, j * vStep + minV)) - 1)
	//}
	//for (let i = 0; i < uRes; i++) {
	//    	let inside = false
	//    	for (let j = 0; j < vRes; j++) {
	//        const pos = i * uRes + j
	//        const fieldCU = minU + (i + 1) * uStep, fieldCV = minV + (j + 1) * vStep
	//        const field = fields[pos]
	//		if (!field) {
	//			if (inside) {
	//				pushQuad(triangles, false,
	//					getVertexIndex(i, j), getVertexIndex(i + 1, j),
	//					getVertexIndex(i, j + 1), getVertexIndex(i + 1, j + 1))
	//			}
	//		} else {
	//			// assemble the field with segments in in
	//			function pos(p) {
	//				const fieldU =
	//				return (p.x < fieldCU ? p.x : uStep + (uStep - ))
	//			}
	//			const ends = []
	//			field.forEach(f => {
	//				ends.push({start: true, part: f, p: f[0], pos: pos(f[0])},
	//					{start: false, part: f, p: f.last(), pos: pos(f.last())})
	//			})
	//			ends.sort((a, b), a.pos - b.pos)
	//			if (flippeds[pos]) {
	//				inside = !inside
	//			}
	//		}
	//       }
	//}
	//    //console.log('trinagle', triangles.max(), vertices.length, triangles.length, triangles.toSource(), triangles.map(i => vertices[i].$).toSource() )
	//    triangles = triangles.map(index => index + mesh.vertices.length)
	//    //assert(normals.every(n => n.hasLength(1)), normals.find(n => !n.hasLength(1)).length() +' '+normals.findIndex(n => !n.hasLength(1)))
	//    Array.prototype.push.apply(mesh.vertices, vertices)
	//    Array.prototype.push.apply(mesh.triangles, triangles)
	//    Array.prototype.push.apply(mesh.normals, normals)
	//    //this.addEdgeLines(mesh)
	//
	//}
	addToMesh2(mesh) {
		const closed = false
		const hSplit = 12800, zSplit = 8
		const ribs = []
		let minZ = Infinity, maxZ = -Infinity
		let cmp = (a, b) => a.value - b.value
		const f = this.surface.parametricFunction()
		const normalF = this.surface.parametricNormal()
		const vertexLoops = this.holes.concat([this.contour]).map(loop => this.unrollLoop(loop))
		vertexLoops.forEach(vertexLoop => {
			vertexLoop.forEach(({x: d, y: z}) => {
				let index0 = ribs.binaryIndexOf(d, (a, b) => NLA.snap(a.value - b, 0))
				if (index0 < 0) {
					ribs.splice(-index0-1, 0, {value: d, left: [], right: []})
				}
				minZ = min(minZ, z)
				maxZ = max(maxZ, z)
			})
		})
		console.log('zzzs', minZ, maxZ, vertexLoops[0].toSource().replace(/\), /g, ',\n'))
		const correction = 1
		vertexLoops.forEach(vertexLoop => {
			vertexLoop.forEach((v0, i, vs) => {
				let v1 = vs[(i + 1) % vs.length], dDiff = v1.x - v0.x
				//console.log(v0.sce, v1.sce)
				if (NLA.eq0(dDiff)) {
					return
				}
				if (dDiff < 0) {
					[v0, v1] = [v1, v0]
					dDiff = -dDiff
				}
				const index0 = ribs.binaryIndexOf(v0.x, (a, b) => NLA.snap(a.value - b, 0))
				const index1 = ribs.binaryIndexOf(v1.x, (a, b) => NLA.snap(a.value - b, 0))
				ribs[index0].right.binaryInsert(v0.y)
				for (let j = (index0 + correction) % ribs.length; j != index1; j = (j + correction) % ribs.length) {
					const x = ribs[j].value
					const part = (x - v0.x) / dDiff
					const interpolated = v1.y * part + v0.y * (1 - part)
					ribs[j].left.binaryInsert(interpolated)
					ribs[j].right.binaryInsert(interpolated)
				}
				ribs[index1].left.binaryInsert(v1.y)
				// console.log(ribs.map(r=>r.toSource()).join('\n'))
			})
		})
		let vertices = [], triangles = [], normals = []
		for (let i = 0; i < ribs.length; i++) {
			let ribLeft = ribs[i], ribRight = ribs[(i + 1) % ribs.length]
			assert(ribLeft.right.length == ribRight.left.length)
			for (let j = 0; j < ribLeft.right.length; j++) {
				vertices.push(f(ribLeft.value, ribLeft.right[j]), f(ribRight.value, ribRight.left[j]))
				normals.push(normalF(ribLeft.value, ribLeft.right[j]), normalF(ribRight.value, ribRight.left[j]))
			}
		}
		//console.log(ribs.map(r=>r.toSource()).join('\n'))
		const vss = vertices.length, detailVerticesStart = vss
		const zInterval = maxZ - minZ, zStep = zInterval / zSplit
		const detailZs = NLA.arrayFromFunction(zSplit - 1, i => minZ + (1 + i) * zStep)
		console.log('detailsZs', detailZs)
		for (let i = 0; i < ribs.length; i++) {
			const d = ribs[i].value
			for (let j = 0; j < detailZs.length; j++) {
				vertices.push(f(d, detailZs[j]))
				normals.push(normalF(d, detailZs[j]))
			}
		}
		// console.log('detailVerticesStart', detailVerticesStart, 'vl', vertices.length, vertices.length - detailVerticesStart, ribs.length)
		// finally, fill in the ribs
		let vsStart = 0
		const flipped2 = true
		//for (var i = 0; i < 1; i++) {
		const end = closed ? ribs.length : ribs.length - 1
		for (let i = 0; i < end; i++) {
			const ipp = (i + 1) % ribs.length
			let inside = false, colPos = 0, ribLeft = ribs[i], ribRight = ribs[(i + 1) % ribs.length]
			for (let j = 0; j < detailZs.length + 1; j++) {
				const detailZ = detailZs[j] || 100000
				if (!inside) {
					if (ribLeft.right[colPos] < detailZ && ribRight.left[colPos] < detailZ) {
						if (ribLeft.right[colPos + 1] < detailZ || ribRight.left[colPos + 1] < detailZ) {
							pushQuad(triangles, flipped2,
								vsStart + colPos * 2,
								vsStart + (colPos + 1) * 2,
								vsStart + colPos * 2 + 1,
								vsStart + (colPos + 1) * 2 + 1)
							colPos += 2
							if (ribLeft.right[colPos] < detailZ || ribRight.left[colPos] < detailZ) {
								j--
							}
						} else {
							pushQuad(triangles, flipped2,
								vsStart + colPos * 2,
								vsStart + colPos * 2 + 1,
								detailVerticesStart + i * detailZs.length + j,
								detailVerticesStart + ipp * detailZs.length + j)
							inside = true
							colPos++
						}
					}
				} else {
					if (ribLeft.right[colPos] < detailZ || ribRight.left[colPos] < detailZ) {
						pushQuad(triangles, flipped2,
							detailVerticesStart + i * detailZs.length + j - 1,
							detailVerticesStart + ipp * detailZs.length + j - 1,
							vsStart + colPos * 2,
							vsStart + colPos * 2 + 1)
						inside = false
						colPos++
						if (ribLeft.right[colPos] < detailZ || ribRight.left[colPos] < detailZ) {
							j--
						}
					} else {
						pushQuad(triangles, flipped2,
							detailVerticesStart + i * detailZs.length + j,
							detailVerticesStart + i * detailZs.length + j - 1,
							detailVerticesStart + ipp * detailZs.length + j,
							detailVerticesStart + ipp * detailZs.length + j - 1)
					}
				}
			}
			vsStart += ribLeft.right.length * 2
		}
		//console.log('trinagle', triangles.max(), vertices.length, triangles.length, triangles.toSource(), triangles.map(i => vertices[i].$).toSource() )
		triangles = triangles.map(index => index + mesh.vertices.length)
		//assert(normals.every(n => n.hasLength(1)), normals.find(n => !n.hasLength(1)).length() +' '+normals.findIndex(n => !n.hasLength(1)))
		Array.prototype.push.apply(mesh.vertices, vertices)
		Array.prototype.push.apply(mesh.triangles, triangles)
		Array.prototype.push.apply(mesh.normals, normals)
		//this.addEdgeLines(mesh)

	}


	intersectFace(face2: PlaneFace,
	              thisBrep: B2,
	              face2Brep: B2,
	              faceMap: Map<Face, Edge[]>,
	              thisEdgePoints: Map<Edge, { edge: Edge, edgeT: number, p: V3, passEdge?: Edge }[]>,
	              otherEdgePoints: Map<Edge, { edge: Edge, edgeT: number, p: V3, passEdge?: Edge }[]>,
	              checkedPairs: Set<NLA.Pair<any, any>>) {

		const __this = this, _thisEdgePoints = thisEdgePoints
		//thisEdgePoints = {
		//   get(key) {
		//       return _thisEdgePoints.get(key)
		//    },
		//    set(key, value) {
		//       assert(thisBrep.edgeFaces.get(key))
		//        _thisEdgePoints.set(key, value)
		//    }
		//}
		function hasPair(a: NLA.Equalable, b: NLA.Equalable) {
			return checkedPairs.has(new NLA.Pair(a, b))
		}
		function addPair(a: NLA.Equalable, b: NLA.Equalable) {
			return checkedPairs.add(new NLA.Pair(a, b))
		}

		/**
		 * @param newEdge generated segment
		 * @param col1 if newEdge is colinear to an edge of this, the edge in question
		 * @param col2 same for face2
		 */
		function handleNewEdge(newEdge: Edge, col1: Edge, col2: Edge) {
			if (!col1 && !col2) {
				if (!(newEdge.aDir.cross(face.surface.normalAt(newEdge.a)).dot(face2.surface.normalAt(newEdge.a)) > 0)) {
					newEdge = newEdge.flipped()
				}
				NLA.mapPush(faceMap, face, newEdge)
				NLA.mapPush(faceMap, face2, newEdge.flipped())
				return true
			}
			function handleEdgeInFace(col1, col2, face, face2, thisBrep, face2Brep, coplanarSameIsInside: boolean, has, add) {
				if (col1 && !col2) {
					if (hasPair(col1.getCanon(), face2)) return

					//add(col1.getCanon(), face2)
					const surface2 = face2.surface

					// NB: a new edge is inserted even though it may be the same as an old one
					// however it indicates that it intersects the other volume here, i.e. the old edge cannot
					// be counted as 'inside' for purposes of reconstitution
					thisBrep.edgeFaces.get(col1.getCanon()).forEach(faceInfo => {
						//const dot = NLA.snap0(surface2.normal1.dot(faceInfo.inside))
						//if (dot == 0 ? !coplanarSameIsInside : dot < 0) {
						const pointsInsideFace = fff(faceInfo, face2.surface)
						const edgeInside = pointsInsideFace == INSIDE || !coplanarSameIsInside && pointsInsideFace == COPLANAR_SAME
						const pushEdge = faceInfo.edge.tangentAt(faceInfo.edge.curve.pointT(newEdge.a)).like(newEdge.aDir) ? newEdge : newEdge.flipped()
						assert(faceInfo.edge.tangentAt(faceInfo.edge.curve.pointT(pushEdge.a)).like(pushEdge.aDir))
						edgeInside && NLA.mapPush(faceMap, faceInfo.face, pushEdge)
					})

					const surface2NormalAtNewEdgeA = surface2.normalAt(newEdge.a)
					const newEdgeInside = surface2NormalAtNewEdgeA.cross(newEdge.aDir)
					const sVEF1 = splitsVolumeEnclosingFacesP(thisBrep, col1.getCanon(), newEdge.a, newEdgeInside, surface2NormalAtNewEdgeA)
					let addNewEdge, addNewEdgeFlipped
					if (addNewEdge = sVEF1 == INSIDE || coplanarSameIsInside && sVEF1 == COPLANAR_SAME) {
						NLA.mapPush(faceMap, face2, newEdge)
					}
					const sVEF2 = splitsVolumeEnclosingFacesP(thisBrep, col1.getCanon(), newEdge.a, newEdgeInside.negated(), surface2NormalAtNewEdgeA)
					if (addNewEdgeFlipped = sVEF2 == INSIDE || coplanarSameIsInside && sVEF2 == COPLANAR_SAME) {
						NLA.mapPush(faceMap, face2, newEdge.flipped())
					}
					if (addNewEdge || addNewEdgeFlipped || sVEF1 == COPLANAR_SAME && sVEF2 == INSIDE || sVEF2 == COPLANAR_SAME && sVEF1 == INSIDE) {
						return true
					}
				}
			}
			const c1 = handleEdgeInFace(col1, col2, face, face2, thisBrep, face2Brep, false, hasPair, addPair)
			const c2 = handleEdgeInFace(col2, col1, face2, face, face2Brep, thisBrep, true, (a, b) => hasPair(b, a), (a, b) => addPair(b, a))
			if (c1 || c2) return true

			if (col1 && col2) {
				if (hasPair(col1.getCanon(), col2.getCanon())) return

				addPair(col1.getCanon(), col2.getCanon())

				function handleColinearEdgeFaces(col1, col2, thisBrep, face2Brep, coplanarSameIsInside: boolean, thisEdgePoints, has, add) {
					// not entirely sure for what i had the dirInsides in?
					//const aDirNegatedInside = (newEdge.a.like(col2.a) || newEdge.a.like(col2.b)) && splitsVolumeEnclosingCone(face2Brep, newEdge.a, newEdge.aDir.negated()) == INSIDE
					//const bDirInside = (newEdge.b.like(col2.a) || newEdge.b.like(col2.b)) && splitsVolumeEnclosingCone(face2Brep, newEdge.b, newEdge.bDir) == INSIDE
					thisBrep.edgeFaces.get(col1.getCanon()).forEach(faceInfo => {
						const sVEF = splitsVolumeEnclosingFaces(face2Brep, col2.getCanon(), faceInfo.inside, faceInfo.normalAtCanonA)
						const edgeInside = sVEF == INSIDE || coplanarSameIsInside && sVEF == COPLANAR_SAME
						const pushEdge = (faceInfo.edge.aDir.like(newEdge.aDir)) ? newEdge : newEdge.flipped()
						edgeInside && NLA.mapPush(faceMap, faceInfo.face, pushEdge)
					})
				}
				handleColinearEdgeFaces(col1, col2, thisBrep, face2Brep, true, thisEdgePoints, hasPair, addPair)
				handleColinearEdgeFaces(col2, col1, face2Brep, thisBrep, false, otherEdgePoints, (a, b) => hasPair(b, a), (a, b) => addPair(b, a))
			}
		}


		// what needs to be generated: new edges on face
		// points on edges where they are cut by faces so that sub edges will be generated for loops
		// points on ends of edges where the edge will be an edge in the new volume where it goes from A to B
		//         you don't want thos to be marked as 'inside', otherwise invalid faces will be added
		// if a face cuts a corner, nothings needs to be done, as that alone does not limit what adjacent faces will be
		function handleEndPoint(a: IntersectionPointInfo, b: IntersectionPointInfo, newEdge: Edge) {
			// ends in the middle of b's face
			if (a && !b) {
				if (!a.colinear && a.edgeT != a.edge.aT && a.edgeT != a.edge.bT) {
					NLA.mapPush(thisEdgePoints, a.edge.getCanon(), a)
					assert(a.edge.isValidT(a.edgeT))
				}
				// else colinear segment ends in middle of other face, do nothing
			}
			// ends in the middle of a's face
			if (b && !a) {
				if (!b.colinear && b.edgeT != b.edge.aT && b.edgeT != b.edge.bT) {
					NLA.mapPush(otherEdgePoints, b.edge.getCanon(), b)
					assert(b.edge.isValidT(b.edgeT))
				}
				// else colinear segment ends in middle of other face, do nothing
			}
			if (a && b) {
				assert(a.colinear || b.colinear || eq(a.t, b.t))
				// if a or b is colinear the correct points will already have been added to the edge by handleNewEdge
				// segment starts/ends on edge/edge intersection
				function foo(a, b, face, face2, thisPlane, face2Plane, thisBrep, face2Brep, first, thisEdgePoints) {
					if (!a.colinear && a.edgeT != a.edge.aT && a.edgeT != a.edge.bT) {
						//if (!hasPair(a.edge.getCanon(), b.edge.getCanon())) {
						addPair(a.edge.getCanon(), b.edge.getCanon())
						// ends on a, on colinear segment b bT != a.edge.bT &&
						// b can be colinear, so edgeT == aT is possible
						if (a.p.like(b.edge.a) || a.p.like(b.edge.b)) {
							const corner = a.p.like(b.edge.a) ? b.edge.a : b.edge.b
							// face2brep corner on edge
							const sVEC1 = splitsVolumeEnclosingCone2(face2Brep, corner, a.edge.curve, a.edgeT, 1)
							const sVEC2 = splitsVolumeEnclosingCone2(face2Brep, corner, a.edge.curve, a.edgeT, -1)
							// if either of these return ALONG_EDGE_OR_PLANE, then the breps share a colinear edge

							if (INSIDE == sVEC1 || INSIDE == sVEC2) {
								NLA.mapPush(thisEdgePoints, a.edge.getCanon(), a)
								assert(a.edge.isValidT(a.edgeT))
							}
						} else {
							// edge / edge center intersection
							// todo: is this even necessary considering we add edges anyway? i think so...
							const testVector = a.edge.tangentAt(a.edgeT).rejectedFrom(b.edge.tangentAt(b.edge.curve.pointT(a.p)))
							assert(!testVector.isZero())
							const sVEF1 = splitsVolumeEnclosingFacesP2(face2Brep, b.edge.getCanon(), a.p, a.edge.curve, a.edgeT, 1, thisPlane.normalAt(a.p))
							const sVEF2 = splitsVolumeEnclosingFacesP2(face2Brep, b.edge.getCanon(), a.p, a.edge.curve, a.edgeT, -1, thisPlane.normalAt(a.p))
							if (INSIDE == sVEF1 || INSIDE == sVEF2) {
								NLA.mapPush(thisEdgePoints, a.edge.getCanon(), a)
								assert(a.edge.isValidT(a.edgeT))
							}
						}
						//}
					}
				}

				foo(a, b, face, face2, surface, surface2, thisBrep, face2Brep, true, thisEdgePoints)
				foo(b, a, face2, face, surface2, surface, face2Brep, thisBrep, false, otherEdgePoints)

			}
		}


		assertInst(Face, face2)


		const face = this
		const surface = face.surface, surface2 = face2.surface
		if (!this.getAABB().fuzzyTouchesAABB(face2.getAABB())) {
			return
		}
		if (surface.isCoplanarTo(surface2)) {
			return
		}
		const isCurves = surface.isCurvesWithSurface(surface2)
		if (0 == isCurves.length) {
			return
		}
		for (const isCurve of isCurves) {
			const t = (isCurve.tMin + isCurve.tMax) / 2, p = isCurve.at(t), dp = isCurve.tangentAt(t)
			const normal1 = surface.normalAt(p), normal2 = surface2.normalAt(p), dp2 = normal1.cross(normal2)
			assert(surface.containsCurve(isCurve))
			assert(surface2.containsCurve(isCurve))
			if (!dp2.isZero()) {
				//assert(dp2.dot(dp) > 0)
				assert(dp2.isParallelTo(dp))
			}
		}

		for (let isCurveIndex = 0; isCurveIndex < isCurves.length; isCurveIndex++ ) {
			// get intersections of newCurve with other edges of face and face2
			const isCurve = isCurves[isCurveIndex]
			const ps1 = faceEdgeISPsWithSurface(face, isCurve, face2.surface)
			const ps2 = faceEdgeISPsWithSurface(face2, isCurve, face.surface)
			// for non-endless curves, e.g. ellipses, the intersections of the faces can be non-zero, even if one of
			// the faces doesn't register any points on the curve. For example, if a cylinder is cut entirely by a
			// plane face (all its edges around the cylinder), then the face will contain the entire curve and
			// 'ps' for the plane face will be empty
			// TODO: behavior when curves touch face?
			// !! start in does depend on insidedir... TODO
			assertf(() => (0 == ps1.length) || !NLA.eq0(ps1[0].insideDir.dot(isCurve.tangentAt(ps1[0].t))), () => ps1[0].insideDir.dot(isCurve.tangentAt(ps1[0].t)))
			assertf(() => (0 == ps2.length) || !NLA.eq0(ps2[0].insideDir.dot(isCurve.tangentAt(ps2[0].t))), () => ps2[0].insideDir.dot(isCurve.tangentAt(ps2[0].t)))
			function startsInside(ps, face) {
				if (0 == ps.length) {
					return isFinite(isCurve.tMin) && face.containsPoint2(isCurve.at(isCurve.tMin)) == PointVsFace.INSIDE
				} else {
					return ps[0].insideDir.dot(isCurve.tangentAt(ps[0].t)) < 0
				}
			}
			// they can't both be empty currently
			// they can't both start 'inside'
			let in1 = startsInside(ps1, face)
			let in2 = startsInside(ps2, face2)
			if (0 == ps1.length && !in1 || 0 == ps2.length && !in2) {
				continue
			}
			//assert(!in1 || !in2)
			let col1: IntersectionPointInfo, col2: IntersectionPointInfo
			let i = 0, j = 0, last
			let startP = in1 && in2 && isCurve.at(isCurve.tMin), startDir, startT = isCurve.tMin, startA, startB
			while (i < ps1.length || j < ps2.length) {
				assert(i <= ps1.length)
				assert(j <= ps2.length)
				const a = ps1[i], b = ps2[j]
				assert(a || b)
				if (j == ps2.length || i < ps1.length && NLA.lt(a.t, b.t)) {
					last = a
					in1 = !in1
					a.used = true
					in1 && (col1 = a.colinear && a)
					i++
				} else if (i == ps1.length || NLA.gt(a.t, b.t)) {
					last = b
					b.used = true
					in2 = !in2
					in2 && (col2 = b.colinear && b)
					j++
				} else {
					last = a
					a.used = true
					b.used = true
					in1 = !in1
					in2 = !in2
					//if (in1 == in2) {
					in1 && (col1 = a.colinear && a)
					in2 && (col2 = b.colinear && b)
					//}
					i++
					j++
				}
				if (startP && !(in1 && in2)) {
					// segment end
					startDir = isCurve.tangentAt(startT)
					if (eq(startT, last.t)) {
						startP = undefined
						continue
					}
					assert(lt(startT, last.t))
					startT > last.t && (startDir = startDir.negated())
					let endDir = isCurve.tangentAt(last.t)
					startT > last.t && (endDir = endDir.negated())
					const newEdge = Edge.create(isCurve, startP, last.p, startT, last.t, null, startDir, endDir, 'genseg' + globalId++)
					startP = undefined
					if (handleNewEdge(newEdge, col1 && col1.edge, col2 && col2.edge)) {
						handleEndPoint(startA || col1, startB || col2, newEdge)
						handleEndPoint(a && a.used && a || col1, b && b.used && b || col2, newEdge)
					}
				} else if (in1 && in2) {
					// new segment just started
					startP = last.p
					startDir = last.insideDir
					startT = last.t
					startA = a && a.used && a
					startB = b && b.used && b
				}
			}
			if (in1 && in2 && startT !== isCurve.tMax) {
				const endT = isCurve.tMax
				startDir = isCurve.tangentAt(startT)
				startT > endT && (startDir = startDir.negated())
				let endDir = isCurve.tangentAt(endT)
				startT > endT && (endDir = endDir.negated())
				const newEdge = Edge.create(isCurve, startP, isCurve.at(endT), startT, endT, null, startDir, endDir, 'genseg' + globalId++)
				if (handleNewEdge(newEdge, col1 && col1.edge, col2 && col2.edge)) {
					handleEndPoint(startA || col1, startB || col2, newEdge)
				}
			}
		}
		face.getAllEdges().forEach(edge => {
			checkedPairs.add(new NLA.Pair(edge.getCanon(), face2))
		})
		face2.getAllEdges().forEach(edge => {
			checkedPairs.add(new NLA.Pair(edge.getCanon(), face))
		})
	}
}

NLA.registerClass(RotationFace)