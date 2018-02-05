import {Equalable, Pair} from 'javasetmap.ts'
import {
	AABB, arrayFromFunction, arrayRange, assert, assertf, assertInst, assertVectors, disableConsole, doubleSignedArea,
	enableConsole, eq, eq0, ge, GOLDEN_RATIO, gt, int, isCCW, le, lerp, lt, M4, mapPush, MINUS, mod, NLA_PRECISION,
	snap, TAU, Transformable, V3,
} from 'ts3dutils'
import {Mesh, pushQuad} from 'tsgl'

import {
	BRep, ConicSurface, COPLANAR_SAME, Curve, dotCurve, dotCurve2, Edge, EllipsoidSurface, fff, getGlobalId, INSIDE,
	IntersectionPointInfo, L3, P3, ParametricSurface, PlaneSurface, PointVsFace, SemiEllipsoidSurface,
	splitsVolumeEnclosingCone2, splitsVolumeEnclosingFaces, splitsVolumeEnclosingFacesP, splitsVolumeEnclosingFacesP2,
	StraightEdge, Surface, triangulateVertices, EPS,
} from './index'

const {PI,  min, max,  sign, ceil, floor, abs} = Math


export abstract class Face extends Transformable {
	'constructor': new (surface: Surface, contour: Edge[], holes?: Edge[][], name?: string, info?: any) => this
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
		assert(contour.every(f => f instanceof Edge), () => 'contour.every(f => f instanceof Edge)' + contour)
		// contour.forEach(e => !surface.containsCurve(e.curve) &&
		// console.log('FAIL:'+surface.distanceToPoint(e.curve.anchor)))
		contour.forEach(e => {
			assert(surface.containsCurve(e.curve), 'edge not in surface ' + e + surface)
		})
		assert(surface.edgeLoopCCW(contour), surface.toString() + contour.join('\n'))
		holes && holes.forEach(hole => Edge.assertLoop(hole))
		holes && holes.forEach(hole => assert(!surface.edgeLoopCCW(hole)))
		assert(!holes || holes.constructor == Array, holes && holes.toString())
		this.allEdges = Array.prototype.concat.apply(this.contour, this.holes)
	}

	static assembleFacesFromLoops(loops: Edge[][], surface: Surface, faceConstructor: typeof Face.prototype.constructor): Face[] {
		type LoopInfo = { loop: Edge[], ccw: boolean, subloops: LoopInfo[] }

		function placeRecursively(newLoopInfo: LoopInfo, loopInfos: LoopInfo[]) {
			if (loopInfos.length == 0) {
				loopInfos.push(newLoopInfo)
			} else {
				const subLoopInfo = loopInfos.find(
					loopInfo => BRep.loop1ContainsLoop2(loopInfo.loop, loopInfo.ccw, newLoopInfo.loop, newLoopInfo.ccw, surface))
				if (subLoopInfo) {
					placeRecursively(newLoopInfo, subLoopInfo.subloops)
				} else {
					// newLoopInfo isnt contained by any other subLoopInfo
					for (let i = loopInfos.length; --i >= 0;) {
						const subLoopInfo = loopInfos[i]
						//console.log('cheving subLoopInfo', surface.loopContainsPoint(newLoopInfo.edges,
						// subLoopInfo.edges[0].a))
						if (BRep.loop1ContainsLoop2(newLoopInfo.loop, newLoopInfo.ccw, subLoopInfo.loop, subLoopInfo.ccw, surface)) {
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
				loopInfo.ccw ? loopInfo.loop : Edge.reversePath(loopInfo.loop),
				loopInfo.subloops.map(sl => sl.ccw ? Edge.reversePath(sl.loop) : sl.loop)))
			loopInfo.subloops.forEach(sl => sl.subloops.forEach(sl2 => newFacesRecursive(sl2)))
		}

		const newFaces: Face[] = []
		const topLevelLoops: LoopInfo[] = []
		loops.forEach(loop => placeRecursively({
			loop: loop,
			ccw: surface.edgeLoopCCW(loop),
			subloops: [],
		}, topLevelLoops))
		topLevelLoops.forEach(tll => newFacesRecursive(tll))
		return newFaces
	}

	//fromLoops(loops: Edge[][], surface: Surface) {
	//	type LoopInfo = {loop: Edge[], ccw: boolean, subloops: LoopInfo[]}
	//	function placeRecursively(newLoopInfo: LoopInfo, loopInfos: LoopInfo[]) {
	//		if (loopInfos.length == 0) {
	//			loopInfos.push(newLoopInfo)
	//		} else {
	//			const subLoopInfo = loopInfos.find(loopInfo => BRep.loop1ContainsLoop2(loopInfo.loop, loopInfo.ccw,
	// newLoopInfo.loop, newLoopInfo.ccw, surface)) if (subLoopInfo) { placeRecursively(newLoopInfo,
	// subLoopInfo.subloops) } else { // newLoopInfo isnt contained by any other subLoopInfo for (let i =
	// loopInfos.length; --i >= 0;) { const subLoopInfo = loopInfos[i] //console.log('cheving subLoopInfo',
	// surface.loopContainsPoint(newLoopInfo.edges, subLoopInfo.edges[0].a)) if
	// (BRep.loop1ContainsLoop2(newLoopInfo.loop, subLoopInfo.loop, surface)) { newLoopInfo.subloops.push(subLoopInfo)
	// loopInfos.splice(i, 1) // remove it } } loopInfos.push(newLoopInfo) } } }  function newFacesRecursive(loopInfo:
	// LoopInfo): void { // CW loops can be top level, if they are holes in the original face not contained in the new
	// face if (loopInfo.ccw) { if (loopInfo.subloops.every(sl => !sl.ccw)) { const newFace = new
	// faceConstructor(surface, loopInfo.loop, loopInfo.subloops.map(sl => sl.loop)) newFaces.push(newFace)
	// loopInfo.subloops.forEach(sl => sl.subloops.forEach(slsl => slsl.ccw && newFacesRecursive(slsl))) } else {
	// loopInfo.subloops.forEach(sl => sl.ccw && newFacesRecursive(sl)) } } }  const newFaces: Face[] = [] const
	// topLevelLoops:LoopInfo[] = [] loops.forEach(loop => placeRecursively({loop: loop, ccw:
	// surface.edgeLoopCCW(loop), subloops: []}, topLevelLoops)) topLevelLoops.forEach(tll => newFacesRecursive(tll))
	// return newFaces }

	static create(surface: Surface, faceEdges: Edge[], holes?: Edge[][], faceName?: string, info?: any) {
		return surface instanceof PlaneSurface
			? new PlaneFace(surface, faceEdges, holes, faceName, info)
			: new RotationFace(surface, faceEdges, holes, faceName, info)
	}

	intersectFace(face2: Face,
				  thisBrep: BRep,
				  face2Brep: BRep,
				  faceMap: Map<Face, Edge[]>,
				  thisEdgePoints: Map<Edge, IntersectionPointInfo[]>,
				  otherEdgePoints: Map<Edge, IntersectionPointInfo[]>,
				  checkedPairs: Set<Pair<any, any>>) {

		//thisEdgePoints = {
		//   get(key) {
		//       return _thisEdgePoints.get(key)
		//    },
		//    set(key, value) {
		//       assert(thisBrep.edgeFaces.get(key))
		//        _thisEdgePoints.set(key, value)
		//    }
		//}
		function hasPair(a: Equalable, b: Equalable) {
			return checkedPairs.has(new Pair(a, b))
		}

		function addPair(a: Equalable, b: Equalable) {
			return checkedPairs.add(new Pair(a, b))
		}

		/**
		 * @param newEdge generated segment
		 * @param col1 if newEdge is colinear to an edge of this, the edge in question
		 * @param col2 same for face2
         * @return whether new edge was added.
		 */
		function handleNewEdge(newEdge: Edge, col1: Edge, col2: Edge): boolean {
			if (!col1 && !col2) {
				let correctDir = face.surface.normalP(newEdge.a).cross(face2.surface.normalP(newEdge.a))
				if (correctDir.likeO()) {
					const t = lerp(newEdge.aT, newEdge.bT, 1 / GOLDEN_RATIO), p = newEdge.curve.at(t)
					correctDir = face.surface.normalP(p).cross(face2.surface.normalP(p))
				}
				if (!correctDir.likeO()) {
					if (correctDir.dot(newEdge.aDir) < 0) {
						newEdge = newEdge.flipped()
					}
					mapPush(faceMap, face, newEdge)
					mapPush(faceMap, face2, newEdge.flipped())
				} else {
					const p = newEdge.a
					const plane = P3.normalOnAnchor(newEdge.aDir, p)
					const up = face.surface.normalP(p)
					const sameDir = up.dot(face2.surface.normalP(p)) > 0
					const canonDir = plane.normal1.cross(up)
					const curve = face.surface.isCurvesWithPlane(plane)[0], curveT = curve.pointT(p),
						curveDir = sign(canonDir.dot(curve.tangentAt(curveT)))
					const curve2 = face2.surface.isCurvesWithPlane(plane)[0], curve2T = curve2.pointT(p),
						curve2Dir = sign(canonDir.dot(curve.tangentAt(curve2T)))
					const foo = curve.diff(curveT, EPS * curveDir).dot(up)
					const foo2 = curve2.diff(curve2T, EPS * curve2Dir).dot(up)
					if (foo2 < foo) {
						mapPush(faceMap, face2, sameDir ? newEdge.flipped() : newEdge)
					}
					if (up.dot(face2.surface.normalP(p)) < 0 == foo2 < foo) {
						mapPush(faceMap, face, newEdge.flipped())
					}
					const bar = curve.diff(curveT, EPS * curveDir).dot(up)
					const bar2 = curve2.diff(curve2T, EPS * curve2Dir).dot(up)
					if (bar2 < bar) {
						mapPush(faceMap, face2, sameDir ? newEdge : newEdge.flipped())
					}
					if (sameDir != bar2 < bar) {
						mapPush(faceMap, face, newEdge)
					}
				}
				return true
			}

			function handleEdgeInFace(
			    col1: Edge | undefined, col2: Edge | undefined,
                face: Face, face2: Face,
                thisBrep: BRep, face2Brep: BRep,
                coplanarSameIsInside: boolean,
                has: typeof hasPair, add: typeof addPair
            ): boolean {
				if (col1 && !col2) {
					if (hasPair(col1.getCanon(), face2)) return false

					//add(col1.getCanon(), face2)
					const surface2 = face2.surface

					// NB: a new edge is inserted even though it may be the same as an old one
					// however it indicates that it intersects the other volume here, i.e. the old edge cannot
					// be counted as 'inside' for purposes of reconstitution
					thisBrep.edgeFaces.get(col1.getCanon()).forEach(faceInfo => {
						//const dot = snap0(surface2.normal1.dot(faceInfo.inside))
						//if (dot == 0 ? !coplanarSameIsInside : dot < 0) {
						const pointsInsideFace = fff(faceInfo, face2.surface)
						const edgeInside = pointsInsideFace == INSIDE || !coplanarSameIsInside && pointsInsideFace == COPLANAR_SAME
						const pushEdge = faceInfo.edge.tangentAt(faceInfo.edge.curve.pointT(newEdge.a)).like(newEdge.aDir)
							? newEdge
							: newEdge.flipped()
						assert(faceInfo.edge.tangentAt(faceInfo.edge.curve.pointT(pushEdge.a)).like(pushEdge.aDir))
						edgeInside && mapPush(faceMap, faceInfo.face, pushEdge)
					})

					const surface2NormalAtNewEdgeA = surface2.normalP(newEdge.a)
					const newEdgeInside = surface2NormalAtNewEdgeA.cross(newEdge.aDir)
					const sVEF1 = splitsVolumeEnclosingFacesP(thisBrep, col1.getCanon(), newEdge.a, newEdgeInside, surface2NormalAtNewEdgeA)
					let addNewEdge, addNewEdgeFlipped
					if (addNewEdge = sVEF1 == INSIDE || coplanarSameIsInside && sVEF1 == COPLANAR_SAME) {
						mapPush(faceMap, face2, newEdge)
					}
					const sVEF2 = splitsVolumeEnclosingFacesP(thisBrep, col1.getCanon(), newEdge.a, newEdgeInside.negated(), surface2NormalAtNewEdgeA)
					if (addNewEdgeFlipped = sVEF2 == INSIDE || coplanarSameIsInside && sVEF2 == COPLANAR_SAME) {
						mapPush(faceMap, face2, newEdge.flipped())
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
				if (hasPair(col1.getCanon(), col2.getCanon())) return false

				addPair(col1.getCanon(), col2.getCanon())
                let added = false
				function handleColinearEdgeFaces(col1: Edge, col2: Edge, thisBrep: BRep, face2Brep: BRep,
					coplanarSameIsInside: boolean, thisEdgePoints: Map<Edge, IntersectionPointInfo[]>,
					 has: typeof hasPair, add: typeof addPair) {
					// not entirely sure for what i had the dirInsides in?
					//const aDirNegatedInside = (newEdge.a.like(col2.a) || newEdge.a.like(col2.b)) &&
					// splitsVolumeEnclosingCone(face2Brep, newEdge.a, newEdge.aDir.negated()) == INSIDE const
					// bDirInside = (newEdge.b.like(col2.a) || newEdge.b.like(col2.b)) &&
					// splitsVolumeEnclosingCone(face2Brep, newEdge.b, newEdge.bDir) == INSIDE
					for (const faceInfo of thisBrep.edgeFaces.get(col1.getCanon())) {
						const sVEF = splitsVolumeEnclosingFaces(face2Brep, col2.getCanon(), faceInfo.inside, faceInfo.normalAtCanonA)
						const edgeInside = sVEF == INSIDE || coplanarSameIsInside && sVEF == COPLANAR_SAME
						const pushEdge = (faceInfo.edge.aDir.like(newEdge.aDir)) ? newEdge : newEdge.flipped()
						if (edgeInside) {
							 mapPush(faceMap, faceInfo.face, pushEdge)
							 const aT = col1.getCanon().curve.pointT(newEdge.a)
							 if (!eq(aT, col1.aT) && !eq(aT, col1.bT)) {
							 	// newEdge.a is in center of col1
							 	if (splitsVolumeEnclosingCone2(face2Brep, newEdge.a, newEdge.curve, newEdge.aT, -Math.sign(newEdge.deltaT()) as -1 | 1) == INSIDE) {
							 		mapPush(thisEdgePoints, col1.getCanon(), {
							 			p: newEdge.a,
							 			edgeT: aT
							 		})
							 	}
							 }
							 const bT = col1.getCanon().curve.pointT(newEdge.b)
							 if (!eq(bT, col1.aT) && !eq(bT, col1.bT)) {
							 	if (splitsVolumeEnclosingCone2(face2Brep, newEdge.b, newEdge.curve, newEdge.bT, Math.sign(newEdge.deltaT()) as -1 | 1) == INSIDE) {
							 		mapPush(thisEdgePoints, col1.getCanon(), {
							 			p: newEdge.b,
							 			edgeT: bT
							 		})
							 	}
							 }
                            added = true
                        }
					}
				}

				handleColinearEdgeFaces(col1, col2, thisBrep, face2Brep, true, thisEdgePoints, hasPair, addPair)
				handleColinearEdgeFaces(col2, col1, face2Brep, thisBrep, false, otherEdgePoints, (a, b) => hasPair(b, a), (a, b) => addPair(b, a))
                return false
			}
		}


		// what needs to be generated: new edges on face
		// points on edges where they are cut by faces so that sub edges will be generated for loops
		// points on ends of edges where the edge will be an edge in the new volume where it goes from A to B
		//         you don't want those to be marked as 'inside', otherwise invalid faces will be added
		// if a face cuts a corner, nothing needs to be done, as that alone does not limit what adjacent faces will be
		function handleEndPoint(a: IntersectionPointInfo | false, b: IntersectionPointInfo | false, newEdge: Edge) {
			// ends in the middle of b's face
			if (a && !b) {
				if (!a.colinear && a.edgeT != a.edge.aT && a.edgeT != a.edge.bT) {
					mapPush(thisEdgePoints, a.edge.getCanon(), a)
					assert(a.edge.isValidT(a.edgeT))
				}
				// else colinear segment ends in middle of other face, do nothing
			}
			// ends in the middle of a's face
			if (b && !a) {
				if (!b.colinear && b.edgeT != b.edge.aT && b.edgeT != b.edge.bT) {
					mapPush(otherEdgePoints, b.edge.getCanon(), b)
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
								mapPush(thisEdgePoints, a.edge.getCanon(), a)
								assert(a.edge.isValidT(a.edgeT))
							}
						} else {
							// edge / edge center intersection
							// todo: is this even necessary considering we add edges anyway? i think so...
							// const testVector =
							// a.edge.tangentAt(a.edgeT).rejectedFrom(b.edge.tangentAt(b.edge.curve.pointT(a.p)))
							// assert(!testVector.likeO())
							const sVEF1 = splitsVolumeEnclosingFacesP2(face2Brep, b.edge.getCanon(), a.p, a.edge.curve, a.edgeT, 1, thisPlane.normalP(a.p))
							const sVEF2 = splitsVolumeEnclosingFacesP2(face2Brep, b.edge.getCanon(), a.p, a.edge.curve, a.edgeT, -1, thisPlane.normalP(a.p))
							if (INSIDE == sVEF1 || INSIDE == sVEF2) {
								mapPush(thisEdgePoints, a.edge.getCanon(), a)
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
			const normal1 = surface.normalP(p), normal2 = surface2.normalP(p), dp2 = normal1.cross(normal2)
			assert(surface.containsCurve(isCurve))
			assert(surface2.containsCurve(isCurve))
			if (!dp2.likeO()) {
				//assert(dp2.dot(dp) > 0)
				// TODO assert(dp2.isParallelTo(dp))
			}
		}

		for (let isCurveIndex = 0; isCurveIndex < isCurves.length; isCurveIndex++) {
			// get intersections of newCurve with other edges of face and face2
			const isCurve = isCurves[isCurveIndex]
			const ps1 = face.edgeISPsWithSurface(isCurve, face2.surface)
			const ps2 = face2.edgeISPsWithSurface(isCurve, face.surface)
			// for non-endless curves, e.g. ellipses, the intersections of the faces can be non-zero, even if one of
			// the faces doesn't register any points on the curve. For example, if a cylinder is cut entirely by a
			// plane face (all its edges around the cylinder), then the face will contain the entire curve and
			// 'ps' for the plane face will be empty
			// TODO: behavior when curves touch face?
			// !! start in does depend on insideDir... TODO
			assertf(() => (0 == ps1.length) || !eq0(ps1[0].insideDir.dot(isCurve.tangentAt(ps1[0].t))), () => ps1[0].insideDir.dot(isCurve.tangentAt(ps1[0].t)))
			assertf(() => (0 == ps2.length) || !eq0(ps2[0].insideDir.dot(isCurve.tangentAt(ps2[0].t))), () => ps2[0].insideDir.dot(isCurve.tangentAt(ps2[0].t)))

			function startsInside(ps: IntersectionPointInfo[], face: Face) {
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
			let col1: IntersectionPointInfo | false, col2: IntersectionPointInfo | false
			let i = 0, j = 0, last
			let startP = in1 && in2 && isCurve.at(isCurve.tMin), startDir, startT = isCurve.tMin, startA, startB
			while (i < ps1.length || j < ps2.length) {
				assert(i <= ps1.length)
				assert(j <= ps2.length)
				const a = ps1[i], b = ps2[j]
				assert(a || b)
				if (j == ps2.length || i < ps1.length && lt(a.t, b.t)) {
					last = a
					in1 = !in1
					a.used = true
					in1 && (col1 = a.colinear && a)
					i++
				} else if (i == ps1.length || gt(a.t, b.t)) {
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
					const newEdge = Edge.create(isCurve, startP, last.p, startT, last.t, undefined, startDir, endDir, 'genseg' + getGlobalId())
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
				const newEdge = Edge.create(isCurve, startP, isCurve.at(endT), startT, endT, undefined, startDir, endDir, 'genseg' + getGlobalId())
				if (handleNewEdge(newEdge, col1 && col1.edge, col2 && col2.edge)) {
					handleEndPoint(startA || col1, startB || col2, newEdge)
				}
			}
		}
		face.getAllEdges().forEach(edge => {
			checkedPairs.add(new Pair(edge.getCanon(), face2))
		})
		face2.getAllEdges().forEach(edge => {
			checkedPairs.add(new Pair(edge.getCanon(), face))
		})
	}

	edgeISPsWithSurface(isCurve: Curve, surface2: Surface): IntersectionPointInfo[] {
		const face = this
		const surface = face.surface
		const loops = face.holes.concat([face.contour])
		const ps = []
		for (const loop of loops) {
			const colinearEdges: boolean[] = loop.map(edge => edge.curve.isColinearTo(isCurve))
			//const colinearSides = loop.map((edge, edgeIndex) => -1 != colinearEdges[edgeIndex]
			//            && -sign(isCurves[colinearEdges[edgeIndex]].tangentAt(edge.aT).dot(edge.aDir)))
			for (let edgeIndex = 0; edgeIndex < loop.length; edgeIndex++) {
				const edge = loop[edgeIndex]
				const nextEdgeIndex = (edgeIndex + 1) % loop.length, nextEdge = loop[nextEdgeIndex]
				//console.log(edge.toSource()) {p:V3(2, -2.102, 0),
				if (colinearEdges[edgeIndex]) {
					if (isCurve.containsPoint(edge.a)) {
						const prevEdgeIndex = (edgeIndex - 1 + loop.length) % loop.length,
							prevEdge = loop[prevEdgeIndex]
						const curveAT = isCurve.pointT(edge.a)
						const colinearOutA = edge.aDir.cross(surface.normalP(edge.a))
						if (!colinearEdges[prevEdgeIndex] && dotCurve2(prevEdge.curve, prevEdge.bT, colinearOutA, -sign(prevEdge.deltaT())) > 0) {
							ps.push({
								p: prevEdge.b,
								insideDir: edge.aDir.negated(),
								t: curveAT,
								edge: prevEdge,
								edgeT: prevEdge.bT,
								colinear: false,
							})
						}
						ps.push({
							p: edge.a,
							insideDir: edge.aDir,
							t: curveAT,
							edge: edge,
							edgeT: edge.aT,
							colinear: true,
						})
					}
					if (isCurve.containsPoint(edge.b)) {
						const curveBT = isCurve.pointT(edge.b)
						const colinearOutB = edge.bDir.cross(surface.normalP(edge.b))
						if (!colinearEdges[nextEdgeIndex] && dotCurve2(nextEdge.curve, nextEdge.aT, colinearOutB, sign(nextEdge.deltaT())) > 0) {
							ps.push({
								p: edge.b,
								insideDir: edge.bDir,
								t: curveBT,
								edge: nextEdge,
								edgeT: nextEdge.aT,
								colinear: false,
							})
						}
						ps.push({
							p: edge.b,
							insideDir: edge.bDir.negated(),
							t: curveBT,
							edge: edge,
							edgeT: edge.bT,
							colinear: true,
						})
					}

				} else {
					const edgeTs = edge.edgeISTsWithSurface(surface2)
					for (const edgeT of edgeTs) {
						const p = edge.curve.at(edgeT)
						if (!isCurve.containsPoint(p)) continue
						const curveT = isCurve.pointT(p)
						assert(!isNaN(curveT))
						const insideDir = edge.tangentAt(edgeT).cross(surface.normalP(p)).negated()

						const isTangent = isCurve.tangentAt(curveT)
						const dirFactor = sign(isTangent.dot(edge.curve.tangentAt(edgeT)))
						const eps = 1e-4
						const normVector = surface2.normalP(p)
						//if(!eq0(insideDir.dot(isTangent))) {
						// Edge.edgeISTsWithSurface returns snapped values, so comparison with == is ok:
						if (edgeT == edge.bT) {
							// endpoint lies on intersection line
							if (!colinearEdges[nextEdgeIndex]) {
								if (!eq(curveT, isCurve.tMax)) {
									const pointsToInside = this.pointsToInside3(edge.b, isCurve, curveT, 1)
									assert(pointsToInside != PointVsFace.ON_EDGE)
									if (PointVsFace.INSIDE == pointsToInside) {
										ps.push({
											p: edge.b,
											insideDir: isTangent,
											t: curveT,
											edge: edge,
											edgeT: edge.bT,
											colinear: false,
										})
									}
								}
								if (!eq(curveT, isCurve.tMin)) {
									const pointsToInside = this.pointsToInside3(edge.b, isCurve, curveT, -1)
									assert(pointsToInside != PointVsFace.ON_EDGE)
									if (PointVsFace.INSIDE == pointsToInside) {
										ps.push({
											p: edge.b,
											insideDir: isTangent.negated(),
											t: curveT,
											edge: edge,
											edgeT: edge.bT,
											colinear: false,
										})
									}
								}
								//let thisSide = -normVector.dot(edge.bDir)
								//if (eq0(thisSide)) {
								//    // advanced test
								//    const dir = -sign(edge.deltaT())
								//    const iscd = isCurve.at(curveT).to(isCurve.at(curveT + dir * dirFactor *
								// eps)).dot(normVector) const ecd = edge.curve.at(edgeT).to(edge.curve.at(edgeT + dir
								// * eps)).dot(normVector) thisSide = sign(ecd - iscd) } let nextSide =
								// normVector.dot(nextEdge.aDir) if (eq0(nextSide)) { // advanced test const dirFactor
								// = sign(snap0(isTangent.dot(nextEdge.curve.tangentAt(nextEdge.aT)))) assert(dirFactor
								// !== 0) const dir = sign(nextEdge.deltaT()) const iscd =
								// isCurve.at(curveT).to(isCurve.at(curveT + dir * dirFactor * eps)).dot(normVector)
								// const ecd = nextEdge.curve.at(nextEdge.aT).to(nextEdge.curve.at(nextEdge.aT + dir *
								// eps)).dot(normVector) nextSide = sign(ecd - iscd) } if (nextSide < 0 || thisSide <
								// 0) { assert(!eq0(insideDir.dot(isTangent))) // next segment is not colinear and ends
								// on different side ps.push({ p: edge.b, insideDir: insideDir, t: curveT, edge: edge,
								// edgeT: edge.bT, colinear: false}) }
							}
						} else if (edgeT != edge.aT) {
							// edge crosses/touches an intersection curve, neither starts nor ends on it
							if (eq0(insideDir.dot(isTangent))) {
								const dirFactor = sign(isTangent.dot(edge.curve.tangentAt(edgeT)))
								const eps = 1e-4
								for (const dir of [-1, 1]) {
									if (-1 == dir * dirFactor && edgeT == edge.minT ||
										1 == dir * dirFactor && edgeT == edge.maxT ||
										-1 == dir && curveT == isCurve.tMin ||
										1 == dir && curveT == isCurve.tMax) continue
									const iscd = isCurve.at(curveT).to(isCurve.at(curveT + dir * eps)).dot(insideDir)
									const ecd = edge.curve.at(edgeT).to(edge.curve.at(edgeT + dir * dirFactor * eps)).dot(insideDir)
									if (iscd - ecd > 0) {
										ps.push({
											p,
											insideDir: isTangent.times(dir),
											t: curveT,
											edge: edge,
											edgeT: edgeT,
											colinear: false,
										})
									}
								}
							} else {
								ps.push({
									p: p,
									insideDir: insideDir,
									t: curveT,
									edge: edge,
									edgeT: edgeT,
									colinear: false,
								})
							}
						}
						//} else {
						//
						//	const dirFactor = sign(isTangent.dot(edge.curve.tangentAt(edgeT)))
						//	const eps = 1e-4
						//	const normVector = surface2.normalP(p)
						//	for (const dir of [-1, 1]) {
						//		if (-1 == dir * dirFactor && edgeT == edge.minT ||
						//			1 == dir * dirFactor && edgeT == edge.maxT ||
						//			-1 == dir && curveT == isCurve.tMin ||
						//			1 == dir && curveT == isCurve.tMax) continue
						//		const iscd = isCurve.at(curveT).to(isCurve.at(curveT + dir * eps)).dot(normVector)
						//		const ecd = edge.curve.at(edgeT).to(edge.curve.at(edgeT + dir * dirFactor *
						// eps)).dot(normVector) if (iscd > ecd) { ps.push({p, insideDir: isTangent.times(dir *
						// dirFactor), t: curveT, edge: edge, edgeT: edgeT, colinear: false}) } }
						// curveVsSurface(isCurve, curveT, p, surface2) }
					}
				}
			}
		}
		// duplicate 't's are ok, as sometimes a segment needs to stop and start again
		// should be sorted so that back facing ones are first
		ps.sort((a, b) => a.t - b.t || a.insideDir.dot(isCurve.tangentAt(a.t)))
		return ps
	}

	transform(m4: M4): Face {
		const mirroring = m4.isMirroring()
		const newEdges = Edge.reversePath(this.contour.map(e => e.transform(m4)), mirroring)
		const newHoles = this.holes.map(hole => Edge.reversePath(hole.map(e => e.transform(m4)), mirroring))
		return new this.constructor(this.surface.transform(m4), newEdges, newHoles, this.name, this.info)
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
		function loopsEqual(a: Edge[], b: Edge[]) {
			return a.length == b.length &&
				arrayRange(0, a.length, 1)
					.some(offset => a.every((edge, i) => edge.equals(b[(offset + i) % a.length])))

		}

		return this == obj ||
			Object.getPrototypeOf(this) == Object.getPrototypeOf(obj)
			&& this.holes.length == obj.holes.length
			&& loopsEqual(this.contour, obj.contour)
			&& this.holes.every(hole => (obj as Face).holes.some(hole2 => loopsEqual(hole, hole2)))
	}

	hashCode() {
		function arrayHashCode(array: number[]) {
			let hashCode = 0
			for (const val of array) {
				hashCode = hashCode * 31 + val | 0
			}
			return hashCode
		}

		function loopHashCode(loop: Edge[]) { return arrayHashCode(loop.map(edge => edge.hashCode()).sort(MINUS)) }

		let hashCode = 0
		hashCode = hashCode * 31 + arrayHashCode(this.holes.map(loop => loopHashCode(loop)).sort(MINUS)) | 0
		hashCode = hashCode * 31 + loopHashCode(this.contour) | 0
		hashCode = hashCode * 31 + this.surface.hashCode() | 0
		return hashCode
	}

	likeFace(face2: Face) {
		function loopsLike(a: Edge[], b: Edge[]) {
			return a.length == b.length &&
				arrayRange(0, a.length, 1)
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

	addEdgeLines(mesh: Mesh) {
		assert(false, 'buggy, fix')
		const vertices = this.contour.flatMap(edge => edge.getVerticesNo0()), mvl = mesh.vertices!.length
		for (let i = 0; i < vertices.length; i++) {
			mesh.vertices!.push(vertices[i])
			mesh.LINES!.push(mvl + i, mvl + (i + 1) % vertices.length)

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
		const mesh = new Mesh()
			.addIndexBuffer('TRIANGLES')
			.addIndexBuffer('LINES')
			.addVertexBuffer('normals', 'LGL_Normal')
		this.addToMesh(mesh)
		//mesh.compile()
		return mesh
	}

	abstract addToMesh(mesh: Mesh & { TRIANGLES: int[], normals: V3[] }): void

	zDirVolume(): { centroid: V3, volume: number } {
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

	pointsToInside3(p: V3, curve: Curve, curveT: number, dir: -1 | 1): PointVsFace {
		const eps = 1e-6
		const normal = this.surface.normalP(p)
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
		if (result == undefined) throw new Error()
		return result
	}

	pointsToInside2(p: V3, dir: V3): PointVsFace {
		return this.pointsToInside3(p, L3.anchorDirection(p, dir), 0, 1)
		//const normal = this.surface.normalP(p)
		//let minAngle = Infinity, inOut = false
		//function test(v, b) {
		//	const angle = (dir.angleRelativeNormal(v, normal) + TAU + NLA_PRECISION / 2) % TAU
		//	if (angle <= 2 * NLA_PRECISION) {
		//		return true
		//	}
		//	if (angle < minAngle) {
		//		minAngle = angle
		//		inOut = b
		//	}
		//}
		//for (const edge of this.getAllEdges()) {
		//	assert(edge.a.equals(p) || !edge.a.like(p))
		//	assert(edge.b.equals(p) || !edge.b.like(p))
		//	if (edge.a.equals(p) && test(edge.aDir, false)) return PointVsFace.ON_EDGE
		//	if (edge.b.equals(p) && test(edge.bDir.negated(), true)) return PointVsFace.ON_EDGE
		//}
		//return inOut ? PointVsFace.INSIDE : PointVsFace.OUTSIDE
	}
}

export class PlaneFace extends Face {

	surface: PlaneSurface

	constructor(p: P3 | PlaneSurface, contour: Edge[], holes?: Edge[][], name?: string, info?: any) {
		assert(p instanceof P3 || p instanceof PlaneSurface)
		super(p instanceof P3 ? new PlaneSurface(p) : p, contour, holes, name, info)
	}

	static forVertices(planeSurface: PlaneSurface | P3, vs: V3[], ...holeVss: V3[][]): PlaneFace {
		const _planeSurface = planeSurface instanceof P3 ? new PlaneSurface(planeSurface) : planeSurface
		assert(isCCW(vs, _planeSurface.plane.normal1), 'isCCW(vs, planeSurface.plane.normal1)')
		const edges = StraightEdge.chain(vs)
		holeVss.forEach(vs => assert(doubleSignedArea(vs, _planeSurface.plane.normal1) >= 0, 'doubleSignedArea(vs, planeSurface.plane.normal1) >= 0'))
		const holes = holeVss.map(hvs => StraightEdge.chain(hvs))
		return new PlaneFace(planeSurface, edges, holes)
	}

	addToMesh(mesh: Mesh & { TRIANGLES: int[], normals: V3[] }) {
		const mvl = mesh.vertices!.length
		const normal = this.surface.plane.normal1
		const vertices = this.contour.flatMap(edge => edge.getVerticesNo0())
		for (let i = 0; i < vertices.length; i++) { mesh.LINES!.push(mvl + i, mvl + (i + 1) % vertices.length) }
		const holeStarts: number[] = []
		this.holes.forEach(hole => {
			holeStarts.push(vertices.length)
			vertices.push(...hole.flatMap(edge => edge.getVerticesNo0()))
		})
		const triangles = triangulateVertices(normal, vertices, holeStarts).map(index => index + mvl)
		Array.prototype.push.apply(mesh.vertices, vertices)
		Array.prototype.push.apply(mesh.TRIANGLES, triangles)
		Array.prototype.push.apply(mesh.normals, arrayFromFunction(vertices.length, () => normal))
	}

	intersectsLine(line: L3): number {
		assertInst(L3, line)
		const lambda = line.isTWithPlane(this.surface.plane)
		if (!Number.isFinite(lambda)) {
			return NaN
		}
		const inside = this.containsPoint(line.at(lambda))
		return inside ? lambda : NaN
	}

	//intersectPlaneFace(face2: PlaneFace,
	//                   thisBrep: BRep,
	//                   face2Brep: BRep,
	//                   faceMap: Map<Face, Edge[]>,
	//                   thisEdgePoints: CustomMap<Edge, { edge: Edge, edgeT: number, p: V3, passEdge?: Edge }[]>,
	//                   otherEdgePoints: CustomMap<Edge, { edge: Edge, edgeT: number, p: V3, passEdge?: Edge }[]>,
	//                   checkedPairs: CustomSet<Pair<Equalable, Equalable>>) {
	//	assertInst(CustomMap, thisEdgePoints, otherEdgePoints)
	//
	//	function hasPair(a: Equalable, b: Equalable) {
	//		return checkedPairs.has(new Pair(a, b))
	//	}
	//	function addPair(a: Equalable, b: Equalable) {
	//		return checkedPairs.add(new Pair(a, b))
	//	}
	//
	//	/**
	//	 * @param newEdge generated segment
	//	 * @param col1 if newEdge is colinear to an edge of this, the edge in question
	//	 * @param col2 same for face2
	//	 */
	//	function handleNewEdge(newEdge: StraightEdge, col1: Edge, col2: Edge) {
	//		if (!col1 && !col2) {
	//			mapPush(faceMap, face, newEdge)
	//			mapPush(faceMap, face2, newEdge.flipped())
	//			return true
	//		}
	//		function handleEdgeInFace(col1, col2, face, face2, thisBrep, face2Brep, coplanarSameIsInside: boolean,
	// has, add) { if (col1 && !col2) { if (hasPair(col1.getCanon(), face2)) return  //add(col1.getCanon(), face2)
	// const face2Plane = face2.surface.plane  // NB: a new edge is inserted even though it may be the same as an old
	// one // however it indicates that it intersects the other volume here, i.e. the old edge cannot // be counted as
	// 'inside' for purposes of reconstitution thisBrep.edgeFaces.get(col1.getCanon()).forEach(faceInfo => { //const
	// dot = snap0(face2Plane.normal1.dot(faceInfo.inside)) //if (dot == 0 ? !coplanarSameIsInside : dot < 0) { const
	// pointsInsideFace = fff(faceInfo, face2.surface) const edgeInside = pointsInsideFace == INSIDE ||
	// !coplanarSameIsInside && pointsInsideFace == COPLANAR_SAME const pushEdge =
	// (faceInfo.edge.aDir.like(newEdge.aDir)) ? newEdge : newEdge.flipped()
	// assert(faceInfo.edge.aDir.like(pushEdge.aDir)) edgeInside && mapPush(faceMap, faceInfo.face, pushEdge) })  const
	// newEdgeInside = face2Plane.normal1.cross(newEdge.aDir) const sVEF1 = splitsVolumeEnclosingFaces(thisBrep,
	// col1.getCanon(), newEdgeInside, face2Plane.normal1) let addNewEdge, addNewEdgeFlipped if (addNewEdge = sVEF1 ==
	// INSIDE || coplanarSameIsInside && sVEF1 == COPLANAR_SAME) { mapPush(faceMap, face2, newEdge) } const sVEF2 =
	// splitsVolumeEnclosingFaces(thisBrep, col1.getCanon(), newEdgeInside.negated(), face2Plane.normal1) if
	// (addNewEdgeFlipped = sVEF2 == INSIDE || coplanarSameIsInside && sVEF2 == COPLANAR_SAME) { mapPush(faceMap,
	// face2, newEdge.flipped()) } if (addNewEdge || addNewEdgeFlipped || sVEF1 == COPLANAR_SAME && sVEF2 == INSIDE ||
	// sVEF2 == COPLANAR_SAME && sVEF1 == INSIDE) { return true } } } const c1 = handleEdgeInFace(col1, col2, face,
	// face2, thisBrep, face2Brep, false, hasPair, addPair) const c2 = handleEdgeInFace(col2, col1, face2, face,
	// face2Brep, thisBrep, true, (a, b) => hasPair(b, a), (a, b) => addPair(b, a)) if (c1 || c2) return true  if (col1
	// && col2) { if (hasPair(col1.getCanon(), col2.getCanon())) return  addPair(col1.getCanon(), col2.getCanon())
	// function handleColinearEdgeFaces(col1, col2, thisBrep, face2Brep, coplanarSameIsInside: boolean, thisEdgePoints,
	// has, add) { // not entirely sure for what i had the dirInsides in? //const aDirNegatedInside =
	// (newEdge.a.like(col2.a) || newEdge.a.like(col2.b)) && splitsVolumeEnclosingCone(face2Brep, newEdge.a,
	// newEdge.aDir.negated()) == INSIDE //const bDirInside = (newEdge.b.like(col2.a) || newEdge.b.like(col2.b)) &&
	// splitsVolumeEnclosingCone(face2Brep, newEdge.b, newEdge.bDir) == INSIDE
	// thisBrep.edgeFaces.get(col1.getCanon()).forEach(faceInfo => { const sVEF = splitsVolumeEnclosingFaces(face2Brep,
	// col2.getCanon(), faceInfo.inside, faceInfo.normalAtCanonA) const edgeInside = sVEF == INSIDE ||
	// coplanarSameIsInside && sVEF == COPLANAR_SAME const pushEdge = (faceInfo.edge.aDir.like(newEdge.aDir)) ? newEdge
	// : newEdge.flipped() edgeInside && mapPush(faceMap, faceInfo.face, pushEdge) }) } handleColinearEdgeFaces(col1,
	// col2, thisBrep, face2Brep, true, thisEdgePoints, hasPair, addPair) handleColinearEdgeFaces(col2, col1,
	// face2Brep, thisBrep, false, otherEdgePoints, (a, b) => hasPair(b, a), (a, b) => addPair(b, a)) } }   // what
	// needs to be generated: new edges on face // points on edges where they are cut by faces so that sub edges will
	// be generated for loops // points on ends of edges where the edge will be an edge in the new volume where it goes
	// from A to B //         you don't want thos to be marked as 'inside', otherwise invalid faces will be added // if
	// a face cuts a corner, nothings needs to be done, as that alone does not limit what adjacent faces will be
	// function handleEndPoint(a: IntersectionPointInfo, b: IntersectionPointInfo, newEdge: Edge) { // ends in the
	// middle of b's face if (a && !b) { if (!a.colinear && a.edgeT != a.edge.aT && a.edgeT != a.edge.bT) {
	// mapPush(thisEdgePoints, a.edge.getCanon(), a) assert(a.edge.isValidT(a.edgeT)) } // else colinear segment ends
	// in middle of other face, do nothing } // ends in the middle of a's face if (b && !a) { if (!b.colinear &&
	// b.edgeT != b.edge.aT && b.edgeT != b.edge.bT) { mapPush(otherEdgePoints, b.edge.getCanon(), b)
	// assert(b.edge.isValidT(b.edgeT)) } // else colinear segment ends in middle of other face, do nothing } if (a &&
	// b) { // if a or b is colinear the correct points will already have been added to the edge by handleNewEdge //
	// segment starts/ends on edge/edge intersection function foo(a, b, face, face2, thisPlane, face2Plane, thisBrep,
	// face2Brep, first, thisEdgePoints) { if (!a.colinear && a.edgeT != a.edge.aT && a.edgeT != a.edge.bT) { if
	// (!hasPair(a.edge.getCanon(), b.edge.getCanon())) { addPair(a.edge.getCanon(), b.edge.getCanon()) // ends on a,
	// on colinear segment b bT != a.edge.bT && // b can be colinear, so edgeT == aT is possible if (a.p.like(b.edge.a)
	// || a.p.like(b.edge.b)) { const corner = a.p.like(b.edge.a) ? b.edge.a : b.edge.b // face2brep corner on edge
	// const sVEC1 = splitsVolumeEnclosingCone(face2Brep, corner, a.edge.aDir) const sVEC2 =
	// splitsVolumeEnclosingCone(face2Brep, corner, a.edge.aDir.negated()) // if either of these return
	// ALONG_EDGE_OR_PLANE, then the breps share a colinear edge  if (INSIDE == sVEC1 || INSIDE == sVEC2) {
	// mapPush(thisEdgePoints, a.edge.getCanon(), a) assert(a.edge.isValidT(a.edgeT)) } } else { // edge / edge center
	// intersection const aEdgeDir = a.edge.tangentAt(a.edgeT) const bEdgeDir = b.edge.tangentAt(b.edgeT) const
	// testVector = aEdgeDir.rejectedFrom(bEdgeDir) assert(!testVector.likeO()) const sVEF1 =
	// splitsVolumeEnclosingFaces(face2Brep, b.edge.getCanon(), testVector, thisPlane.normal1) const sVEF2 =
	// splitsVolumeEnclosingFaces(face2Brep, b.edge.getCanon(), testVector.negated(), thisPlane.normal1) if (INSIDE ==
	// sVEF1 || INSIDE == sVEF2) { mapPush(thisEdgePoints, a.edge.getCanon(), a) assert(a.edge.isValidT(a.edgeT)) } } }
	// } }  foo(a, b, face, face2, thisPlane, face2Plane, thisBrep, face2Brep, true, thisEdgePoints) foo(b, a, face2,
	// face, face2Plane, thisPlane, face2Brep, thisBrep, false, otherEdgePoints)  } }   assertInst(PlaneFace, face2)
	// const face: PlaneFace = this // get intersection const thisPlane = this.surface.plane, face2Plane =
	// face2.surface.plane if (thisPlane.isParallelToPlane(face2Plane)) { if (thisPlane.like(face2Plane)) { // normal1
	// same and same location in space // addLikeSurfaceFaces(likeSurfaceFaces, this, face2) } return } const isLine =
	// L3.fromPlanes(thisPlane, face2Plane) // get intersections of newCurve with other edges of face and face2 const
	// ps1 = planeFaceEdgeISPsWithPlane(face, isLine, face2Plane) const ps2 = planeFaceEdgeISPsWithPlane(face2, isLine,
	// thisPlane) if (ps1.length == 0 || ps2.length == 0) { // faces to not intersect return }  let col1:
	// IntersectionPointInfo, col2: IntersectionPointInfo let in1 = false, in2 = false let i = 0, j = 0, last let
	// startP, startDir, startT, startA, startB while (i < ps1.length || j < ps2.length) { assert(i <= ps1.length)
	// assert(j <= ps2.length) const a = ps1[i], b = ps2[j] assert(a || b) if (j == ps2.length || i < ps1.length &&
	// lt(a.t, b.t)) { last = a in1 = !in1 a.used = true in1 && (col1 = a.colinear && a) i++ } else if (i == ps1.length
	// || gt(a.t, b.t)) { last = b in2 = !in2 b.used = true in2 && (col2 = b.colinear && b) j++ } else { // TODO: this
	// will break if 3 points on the same t last = a in1 = !in1 in2 = !in2 //if (in1 == in2) { a.used = true b.used =
	// true in1 && (col1 = a.colinear && a) in2 && (col2 = b.colinear && b) //} i++ j++ } if (startP && !(in1 && in2))
	// { // segment end const newEdge = new StraightEdge(isLine, startP, last.p, startT, last.t, undefined, 'genseg' +
	// getGlobalId()) startP = undefined last.used = true if (handleNewEdge(newEdge, col1 && col1.edge, col2 &&
	// col2.edge)) { handleEndPoint(startA || col1, startB || col2, newEdge) handleEndPoint(a && a.used && a || col1, b
	// && b.used && b || col2, newEdge) } } else if (in1 && in2) { // new segment just started startP = last.p startDir
	// = last.insideDir startT = last.t startA = a && a.used && a startB = b && b.used && b } if (!in1 && a && last ==
	// a && a.colinear) { checkedPairs.add(new Pair(a.edge.getCanon(), face2)) } if (!in2 && b && (last == b || b.used)
	// && b.colinear) { checkedPairs.add(new Pair(b.edge.getCanon(), face)) } } }

	withHole(holeEdges: Edge[]) {
		return new PlaneFace(this.surface, this.contour, [holeEdges])
	}

	pointsToInside(p: V3, dir: V3): PointVsFace {
		return this.containsPoint2(p.plus(dir.times(NLA_PRECISION * 8)))
	}

	edgeISPsWithPlane(isLine: L3, plane2: P3): IntersectionPointInfo[] {
		const face = this
		assert(face.surface.plane.containsLine(isLine))
		assert(plane2.containsLine(isLine))
		const plane = face.surface.plane
		const ps: IntersectionPointInfo[] = []
		const loops = [face.contour].concat(face.holes)
		loops.forEach(loop => {
			const colinearEdges = loop.map((edge) => edge.colinearToLine(isLine) && -sign(edge.aDir.dot(isLine.dir1)))
			const isLineOut = isLine.dir1.cross(plane.normal1)

			loop.forEach((edge, edgeIndex, edges) => {
				const nextEdgeIndex = (edgeIndex + 1) % edges.length, nextEdge = edges[nextEdgeIndex]
				//console.log(edge.toSource()) {p:V3(2, -2.102, 0),
				if (colinearEdges[edgeIndex]) {
					// edge colinear to intersection line
					const curveAT = isLine.pointT(edge.a), curveBT = isLine.pointT(edge.b)
					// add interval for colinear segment
					ps.push(
						{p: edge.a, insideDir: edge.aDir, t: curveAT, edge: edge, edgeT: edge.aT, colinear: true},
						{
							p: edge.b,
							insideDir: edge.bDir.negated(),
							t: curveBT,
							edge: edge,
							edgeT: edge.bT,
							colinear: true,
						})
					// open next interval if necessary
					const nextSide = colinearEdges[nextEdgeIndex] || dotCurve(isLineOut, nextEdge.aDir, nextEdge.aDDT)
					if (colinearEdges[edgeIndex] * nextSide < 0) {
						// side changes
						ps.push({
							p: nextEdge.a,
							insideDir: edge.bDir,
							t: curveBT,
							edge: nextEdge,
							edgeT: nextEdge.aT,
							colinear: false,
						})
					}
				} else {
					// not necessarily a straight edge, so multiple intersections are possible
					const edgeTs = edge.edgeISTsWithPlane(plane2)
					assert(edgeTs.every(t => plane2.containsPoint(edge.curve.at(t))), edgeTs)
					for (const edgeT of edgeTs) {
						if (edgeT == edge.bT) {
							// endpoint lies on intersection line
							const side = -dotCurve(isLineOut, edge.bDir, edge.bDDT)
							const nextSide = colinearEdges[nextEdgeIndex] || dotCurve(isLineOut, nextEdge.aDir, nextEdge.aDDT)
							if (side * nextSide < 0) {
								// next segment is not colinear and ends on different side
								ps.push({
									p: edge.b,
									insideDir: plane2.normal1.negated(),
									t: isLine.pointT(edge.b),
									edge: edge,
									edgeT: edge.bT,
									colinear: false,
								})
							}
						} else if (edgeT != edge.aT) {
							// edge crosses intersection line, neither starts nor ends on it
							const p = edge.curve.at(edgeT)
							assert(plane2.containsPoint(p), edge.toString(), p, edgeT, plane2.distanceToPoint(p))
							assert(isLine.containsPoint(p), edge.toString(), p, edgeT, isLine.distanceToPoint(p))
							const insideDir = plane2.normal1.negated()
							ps.push({
								p: p,
								insideDir: insideDir,
								t: isLine.pointT(p),
								edge: edge,
								edgeT: edgeT,
								colinear: false,
							})
						}
					}
				}
			})
		})
		// duplicate 't's are ok, as sometimes a segment needs to stop and start again
		// should be sorted so that back facing ones are first
		ps.sort((a, b) => a.t - b.t || a.insideDir.dot(isLine.dir1))
		return ps
	}
}


export class RotationFace extends Face {
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
			this.aabb = AABB.forAABBs(this.contour.map(e => e.getAABB()))
			this.aabb.addPoints(this.surface.getExtremePoints().filter(p => this.containsPoint(p)))
			return this.aabb
		} else {
			return super.getAABB()
		}
	}

	getCanonSeamU(): number {
		const stPFunc = this.surface.stPFunc()
		for (const edge of this.contour) {
			// check edge.a
			let u = stPFunc(edge.a, PI).x
			// if u is not PI, or ~0, return its sign
			if (u != PI && !eq0(u)) {
				return sign(u) * PI
			}
			// check midpoint between edge.a and edge.b
			u = stPFunc(edge.curve.at((edge.aT + edge.bT) / 2), PI).x
			if (u != PI && !eq0(u)) {
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


	unrollLoop(edgeLoop: Edge[]) {
		const vs: V3[] = []
		const reverseFunc = this.surface.stPFunc()
		const verticesNo0s = edgeLoop.map(edge => edge.getVerticesNo0())
		const startEdgeIndex = verticesNo0s.findIndex(edgeVertices => !eq(reverseFunc(edgeVertices[0], Math.PI).x, Math.PI))
		assert(-1 != startEdgeIndex)
		// console.log(startEdgeIndex)
		let hint = Math.PI
		for (let i = 0; i < edgeLoop.length; i++) {
			const edgeIndex = (i + startEdgeIndex) % edgeLoop.length
			for (let j = 0; j < verticesNo0s[edgeIndex].length; j++) {
				const p = verticesNo0s[edgeIndex][j]
				const localP = reverseFunc(p, hint)
				if (Math.abs(localP.x) < Math.PI - NLA_PRECISION) {
					// update hint
					hint = localP.x
				}
				// console.log(hint, p.sce, localP.sce)
				vs.push(localP)
			}
		}
		edgeLoop.forEach((edge, e) => {
			let hint = edge.bDir
			if (edge instanceof StraightEdge && edge.curve.dir1.isParallelTo(this.surface.dir || this.surface.dir1)) {
				hint = this.surface.normalP(edge.b).cross(edge.bDir)
			}
			edge.getVerticesNo0().forEach(p => {
				vs.push(reverseFunc(p, hint))
			})
		})
		console.log('vs\n', vs.join('\n'), vs.length)
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
		const verticesST: V3[] = [], vertices: V3[] = [], loopStarts = []
		const ellipsoid: SemiEllipsoidSurface = this.surface as SemiEllipsoidSurface
		const ptpf = ellipsoid.stPFunc()
		const testDegeneratePoint = ellipsoid instanceof SemiEllipsoidSurface
			? (nextStart: V3) => nextStart.like(ellipsoid.center.plus(ellipsoid.f3)) || nextStart.like(ellipsoid.center.minus(ellipsoid.f3))
			: (nextStart: V3) => nextStart.like((this.surface as ConicSurface).center)
		for (const edgeLoop of edgeLoops) {
			loopStarts.push(verticesST.length)
			// console.log(startEdgeIndex)
			const hint = this.getCanonSeamU()
			for (let i = 0; i < edgeLoop.length; i++) {
				const ipp = (i + 1) % edgeLoop.length
				const verticesNo0 = edgeLoop[i].getVerticesNo0()
				vertices.push(...verticesNo0)
				verticesST.push(...verticesNo0.map(v => ptpf(v)))
				const nextStart = edgeLoop[ipp].a
				//console.log('BLAH', nextStart.str, ellipsoid.center.plus(ellipsoid.f3).str)

				if (testDegeneratePoint(nextStart)) {
					const bDirLC = ellipsoid.inverseMatrix.transformVector(edgeLoop[i].bDir),
						aDirLC = ellipsoid.inverseMatrix.transformVector(edgeLoop[ipp].aDir)
					let inAngle = Math.atan2(-bDirLC.y, -bDirLC.x)
					if (abs(inAngle) > Math.PI - NLA_PRECISION) {
						assert(hint == -PI || hint == PI)
						inAngle = hint
					}
					let outAngle = Math.atan2(aDirLC.y, aDirLC.x)
					if (abs(outAngle) > Math.PI - NLA_PRECISION) {
						assert(hint == -PI || hint == PI)
						outAngle = hint
					}

					const stLast = verticesST.pop()!
					verticesST.push(new V3(inAngle, stLast.y, 0), new V3(outAngle, stLast.y, 0))
					vertices.push(vertices.last)
				}
				verticesST.forEach(({x: u, y: v}) => {
					assert(isFinite(u))
					assert(isFinite(v))
				})
			}
		}
		let normals
		if (this.surface instanceof EllipsoidSurface) {
			normals = vertices.map(v => ellipsoid.normalP(v))
		} else {
			const pN = ellipsoid.normalSTFunc()
			normals = verticesST.map(({x, y}) => pN(x, y))
		}
		assert(vertices.length == vertices.length)
		//console.log(verticesST.map(v => v.str).join('\n'))
		return {
			verticesUV: verticesST.map(vST => new V3(vST.x / uStep, vST.y / vStep, 0)),
			vertices: vertices,
			normals: normals,
			loopStarts: loopStarts,
		}
	}

	unrollCylinderLoops(loops: Edge[][], uStep: number, vStep: number) {
		const vertexLoops = loops.map(loop => loop.flatMap(edge => edge.getVerticesNo0()))
		const surface = this.surface as ParametricSurface
		const vertices: V3[] = vertexLoops.concatenated()
		// this.unrollLoop(loop).map(v => new V3(v.x / uStep, v.y / vStep, 0)))
		const loopStarts = vertexLoops.reduce((arr, loop) => (arr.push(arr.last + loop.length), arr), [0])
		const stPFunc = surface.stPFunc()
		const verticesST = vertices.map(v => stPFunc(v))
		const verticesUV = verticesST.map(st => new V3(st.x / uStep, st.y / vStep, 0))
		const normalST = surface.normalSTFunc()
		const normals: V3[] = verticesST.map(({x, y}) => normalST(x, y))
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

	addToMesh(this: this & { surface: ParametricSurface }, mesh: Mesh, uStep: number = this.surface.uStep, vStep: number = this.surface.vStep) {
		assertf(() => uStep > 0 && vStep > 0, uStep, vStep, 'Surface: ' + this.surface)
		const triangles: int[] = []
		const pIJFunc = (i: number, j: number) => this.surface.pSTFunc()(i * uStep, j * vStep)
		const normalIJFunc = (i: number, j: number) => this.surface.normalSTFunc()(i * uStep, j * vStep)
		const loops = [this.contour].concat(this.holes)
		const {vertices, verticesUV, normals, loopStarts} = this.surface instanceof SemiEllipsoidSurface || this.surface instanceof ConicSurface
			? this.unrollEllipsoidLoops(loops, uStep, vStep)
			: this.unrollCylinderLoops(loops, uStep, vStep)
		loopStarts.push(vertices.length)

		for (let vertexLoopIndex = 0; vertexLoopIndex < loops.length; vertexLoopIndex++) {
			const vertexLoopStart = loopStarts[vertexLoopIndex]
			const vertexLoopLength = loopStarts[vertexLoopIndex + 1] - vertexLoopStart
			const base = mesh.vertices!.length + loopStarts[vertexLoopIndex]
			for (let i = 0; i < vertexLoopLength; i++) {
				mesh.LINES!.push(base + i, base + (i + 1) % vertexLoopLength)
			}
		}

		disableConsole()
		let minU = Infinity, maxU = -Infinity, minV = Infinity, maxV = -Infinity
		//console.log('surface', this.surface.str)
		//console.log(verticesUV)
		//drPs.push(...verticesUV.map((v, i) => ({p: vertices[i], text: `${i} uv: ${v.toString(x => round10(x,
		// -4))}`})))
		verticesUV.forEach(({x: u, y: v}) => {
			assert(isFinite(u))
			assert(isFinite(v))
			minU = min(minU, u)
			maxU = max(maxU, u)
			minV = min(minV, v)
			maxV = max(maxV, v)
		})
		if (ParametricSurface.is(this.surface)) {
			assert(this.surface.boundsSigned(minU * uStep, minV * vStep) > -NLA_PRECISION)
			assert(this.surface.boundsSigned(maxU * uStep, maxV * vStep) > -NLA_PRECISION)
		}
		const uOffset = floor(minU + NLA_PRECISION), vOffset = floor(minV + NLA_PRECISION)
		const uRes = ceil(maxU - NLA_PRECISION) - uOffset, vRes = ceil(maxV - NLA_PRECISION) - vOffset
		console.log(uStep, vStep, uRes, vRes)
		if (uRes == 1 && vRes == 1) {
			// triangulate this face as if it were a plane
			const polyTriangles = triangulateVertices(V3.Z, verticesUV, loopStarts.slice(1, 1 + this.holes.length))
			triangles.push(...polyTriangles)
		} else {
			const partss: int[][][] = new Array(uRes * vRes)

			function fixUpPart(part: number[], baseU: int, baseV: int) {
				assert(baseU < uRes && baseV < vRes, `${baseU}, ${baseV}, ${uRes}, ${vRes}`)
				console.log('complete part', part, baseU, baseV)
				//console.trace()
				assert(part.length)
				const cellU = baseU + uOffset, cellV = baseV + vOffset
				for (const index of part) {
					assert(le(cellU, verticesUV[index].x) && le(verticesUV[index].x, cellU + 1), `${index} ${verticesUV[index].str} ${cellU} ${cellU}`)
					assert(le(cellV, verticesUV[index].y) && le(verticesUV[index].y, cellV + 1))
				}
				const pos = baseV * uRes + baseU
				;(partss[pos] || (partss[pos] = [])).push(part)
				//const outline = partss[pos] || (partss[pos] = [minU + baseU * uStep, minV + baseV * vStep, minU +
				// (baseU + 1) * uStep, minV + (baseV + 1) * vStep])
			}

			// 'some' instead of forEach so we can return out of the entire function if this.edges crosses no borders
			// and
			for (let vertexLoopIndex = 0; vertexLoopIndex < loops.length; vertexLoopIndex++) {
				let part: int[] | undefined = undefined, firstPart, firstPartBaseU: int = -1, firstPartBaseV: int = -1
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
						// points which are on a grid line are assigned to the cell into which they are going (+
						// NLA_PRECISION * sign(di)) if they are parallel to the gridline (eq0(di)), they belong the
						// the cell for which they are a CCW boundary
						const baseU = floor(vxu + (!eq0(di) ? sign(di) : -sign(dj)) * NLA_PRECISION) - uOffset
						const baseV = floor(vxv + (!eq0(dj) ? sign(dj) : sign(di)) * NLA_PRECISION) - vOffset
						assert(baseU < uRes && baseV < vRes, `${baseU}, ${baseV}, ${uRes}, ${vRes}`)
						// figure out the next intersection with a gridline:
						// iNext is the positive horizontal distance to the next vertical gridline
						const iNext = ceil(sign(di) * vxu + NLA_PRECISION) - sign(di) * vxu
						const jNext = ceil(sign(dj) * vxv + NLA_PRECISION) - sign(dj) * vxv
						const iNextT = currentT + iNext / abs(di)
						const jNextT = currentT + jNext / abs(dj)
						//console.log(vxIndex, vx.str, 'vij', vxu, vxv, 'd', di, dj, 'ijNext', iNext, jNext, 'nextT',
						// iNextT, jNextT)
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
						if (ge(currentT, 1)) {
							//console.log('breaking ', vx1index)
							part!.push(vx1index)
							break
						} else {
							const nextPoint = vx0.lerp(vx1, currentT)
							const nextPointIndex = addVertex(nextPoint.x, nextPoint.y)

							//console.log('pushing ', nextPointIndex)
							part!.push(nextPointIndex)
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
					part!.pop()
					fixUpPart(part!.concat(firstPart!), lastBaseU, lastBaseV)
				} else {
					fixUpPart(firstPart!, firstPartBaseU!, firstPartBaseV!)
					fixUpPart(part!, lastBaseU, lastBaseV)
				}
				console.log('firstPart', firstPart)
			}
			console.log('calculated parts', partss)
			const fieldVertexIndices = new Array((uRes + 1) * (vRes + 1))

			function addVertex(u: number, v: number): int {
				verticesUV.push(new V3(u, v, 0))
				normals.push(normalIJFunc(u, v))
				return vertices.push(pIJFunc(u, v)) - 1
			}

			function getGridVertexIndex(i: int, j: int): int {
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
							const outline: int[] = [], outlineVertexIndices = []
							const startPart = parts[0]
							assert(startPart.length > 0)
							let currentPart = startPart
							do {
								outline.push(...currentPart)
								const currentPartEndOpos = opos(currentPart.last)
								const nextPartIndex = parts.indexWithMax(part => -mod(opos(part[0]) - currentPartEndOpos, 4))
								const nextPart = parts.removeIndex(nextPartIndex)
								let currentOpos = currentPartEndOpos
								const nextPartStartOpos = opos(nextPart[0]) > currentOpos
									? opos(nextPart[0])
									: opos(nextPart[0]) + 4
								let nextOpos = ceil(currentOpos + NLA_PRECISION)
								let flipping = eq0((currentOpos + NLA_PRECISION) % 1 - NLA_PRECISION)
								//inside = inside != (!eq0(currentOpos % 1) && currentOpos % 2 < 1)
								while (lt(nextOpos, nextPartStartOpos)) {
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
								triangles.push(...outline)
							} else {
								const polyTriangles = triangulateVertices(V3.Z, outline.map(i => verticesUV[i]), []).map(i => outline[i])
								triangles.push(...polyTriangles)
							}
							//console.log('outline', col, row, outline)
						}
					}
				}
			}

		}
		//console.log('trinagle', triangles.max(), vertices.length, triangles.length, triangles.toSource(),
		// triangles.map(col => vertices[col].$).toSource() ) assert(normals.every(n => n.hasLength(1)), normals.find(n
		// => !n.hasLength(1)).length() +' '+normals.findIndex(n => !n.hasLength(1)))
		Array.prototype.push.apply(mesh.TRIANGLES, triangles.map(index => index + mesh.vertices!.length))
		Array.prototype.push.apply(mesh.vertices, vertices)
		Array.prototype.push.apply(mesh.normals, normals)
		//this.addEdgeLines(mesh)
		enableConsole()
	}

	addToMesh2(this: this & { surface: ParametricSurface }, mesh: Mesh) {
		const closed = false
		const hSplit = 12800, zSplit = 8
		const ribs: { value: number, left: number[], right: number[] }[] = []
		let minZ = Infinity, maxZ = -Infinity
		//let cmp = (a, b) => a.value - b.value
		const f = this.surface.pSTFunc()
		const normalF = this.surface.normalSTFunc()
		const vertexLoops = this.holes.concat([this.contour]).map(loop => this.unrollLoop(loop))
		vertexLoops.forEach(vertexLoop => {
			vertexLoop.forEach(({x: d, y: z}) => {
				const index0 = ribs.binaryIndexOf(d, (a, b) => snap(a.value - b, 0))
				if (index0 < 0) {
					ribs.splice(-index0 - 1, 0, {value: d, left: [], right: []})
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
				if (eq0(dDiff)) {
					return
				}
				if (dDiff < 0) {
					[v0, v1] = [v1, v0]
					dDiff = -dDiff
				}
				const index0 = ribs.binaryIndexOf(v0.x, (a, b) => snap(a.value - b, 0))
				const index1 = ribs.binaryIndexOf(v1.x, (a, b) => snap(a.value - b, 0))
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
		const vertices = [], triangles0: int[] = [], normals = []
		for (let i = 0; i < ribs.length; i++) {
			const ribLeft = ribs[i], ribRight = ribs[(i + 1) % ribs.length]
			assert(ribLeft.right.length == ribRight.left.length)
			for (let j = 0; j < ribLeft.right.length; j++) {
				vertices.push(f(ribLeft.value, ribLeft.right[j]), f(ribRight.value, ribRight.left[j]))
				normals.push(normalF(ribLeft.value, ribLeft.right[j]), normalF(ribRight.value, ribRight.left[j]))
			}
		}
		//console.log(ribs.map(r=>r.toSource()).join('\n'))
		const vss = vertices.length, detailVerticesStart = vss
		const zInterval = maxZ - minZ, zStep = zInterval / zSplit
		const detailZs = arrayFromFunction(zSplit - 1, i => minZ + (1 + i) * zStep)
		console.log('detailsZs', detailZs)
		for (let i = 0; i < ribs.length; i++) {
			const d = ribs[i].value
			for (let j = 0; j < detailZs.length; j++) {
				vertices.push(f(d, detailZs[j]))
				normals.push(normalF(d, detailZs[j]))
			}
		}
		// console.log('detailVerticesStart', detailVerticesStart, 'vl', vertices.length, vertices.length -
		// detailVerticesStart, ribs.length) finally, fill in the ribs
		let vsStart = 0
		const flipped2 = true
		//for (var i = 0; i < 1; i++) {
		const end = closed ? ribs.length : ribs.length - 1
		for (let i = 0; i < end; i++) {
			const ipp = (i + 1) % ribs.length
			let inside = false, colPos = 0
			const ribLeft = ribs[i], ribRight = ribs[(i + 1) % ribs.length]
			for (let j = 0; j < detailZs.length + 1; j++) {
				const detailZ = detailZs[j] || 100000
				if (!inside) {
					if (ribLeft.right[colPos] < detailZ && ribRight.left[colPos] < detailZ) {
						if (ribLeft.right[colPos + 1] < detailZ || ribRight.left[colPos + 1] < detailZ) {
							pushQuad(triangles0, flipped2,
								vsStart + colPos * 2,
								vsStart + (colPos + 1) * 2,
								vsStart + colPos * 2 + 1,
								vsStart + (colPos + 1) * 2 + 1)
							colPos += 2
							if (ribLeft.right[colPos] < detailZ || ribRight.left[colPos] < detailZ) {
								j--
							}
						} else {
							pushQuad(triangles0, flipped2,
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
						pushQuad(triangles0, flipped2,
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
						pushQuad(triangles0, flipped2,
							detailVerticesStart + i * detailZs.length + j,
							detailVerticesStart + i * detailZs.length + j - 1,
							detailVerticesStart + ipp * detailZs.length + j,
							detailVerticesStart + ipp * detailZs.length + j - 1)
					}
				}
			}
			vsStart += ribLeft.right.length * 2
		}
		//console.log('trinagle', triangles0.max(), vertices.length, triangles0.length, triangles0.toSource(),
		// triangles0.map(i => vertices[i].$).toSource() )
		const triangles = triangles0.map(index => index + mesh.vertices!.length)
		//assert(normals.every(n => n.hasLength(1)), normals.find(n => !n.hasLength(1)).length() +'
		// '+normals.findIndex(n => !n.hasLength(1)))
		Array.prototype.push.apply(mesh.vertices, vertices)
		Array.prototype.push.apply(mesh.TRIANGLES, triangles)
		Array.prototype.push.apply(mesh.normals, normals)
		//this.addEdgeLines(mesh)

	}


}

