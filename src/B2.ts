"use strict";
import registerClass = NLA.registerClass;

var eps = 1e-5


class Face extends Transformable {
    surface: Surface
    edges: Edge[]
    holes: Edge[][]
    id: int
    name: string

	constructor(surface:Surface, contour:Edge[], holes:Edge[][], name?:string) {
		super()
		this.assertChain(contour)
		//assert(surface.edgeLoopCCW(contour), surface.toString()+contour.join("\n"))
		assert(contour.every(f => f instanceof Edge), 'contour.every(f => f instanceof Edge)' + contour.toSource())
		// contour.forEach(e => !surface.containsCurve(e.curve) && console.log("FAIL:"+surface.distanceToPoint(e.curve.anchor)))
		contour.forEach(e => assert(surface.containsCurve(e.curve), e + surface))
		holes && holes.forEach(hole => this.assertChain(hole))
		//holes && holes.forEach(hole => assert(!surface.edgeLoopCCW(hole)))
		assert (!holes || holes.constructor == Array, holes && holes.toString())
		this.surface = surface
		this.edges = contour // TODO refactor to contour
		this.holes = holes || []
		this.id = globalId++
		this.name = name
	}

	transform(m4) {
		var newEdges = this.edges.map(e => e.transform(m4))
		var newHoles = this.holes.map(hole => hole.map(e => e.transform(m4)))
		return new this.constructor(this.surface.transform(m4), newEdges, newHoles)
	}

	assertChain(edges) {
		edges.forEach((edge, i) => {
			var j = (i + 1) % edges.length
			assert(edge.b.like(edges[j].a), `edges[${i}].b != edges[${j}].a (${edges[i].b.sce} != ${edges[j].a.sce})`)
		})
	}

	flipped() {
		var newEdges = this.edges.map(e => e.flipped()).reverse()
		var newHoles = this.holes.map(hole => hole.map(e => e.flipped()).reverse())
		return new this.constructor(this.surface.flipped(), newEdges, newHoles, this.name)
	}

	toString() {
		return `new ${this.constructor.name}(${this.surface}, [${this.edges.map(e => '\n\t' + e).join()}]`
			+`${this.holes.map(hole => '\n\t\thole: ' + hole.join())})`
	}

	toSource() {
		return `new ${this.constructor.name}(${this.surface.toSource()}, [${this.edges.map(e => '\n\t' + e.toSource()).join(',')}], [${
			this.holes.map(hole => '['+hole.map(e => '\n\t' + e.toSource()).join(',')+']').join(',')}])`
	}

	equals(face) {
		//TODO		assert(false)
		var edgeCount = this.edges.length

		return this.surface.equalsSurface(face.surface) &&
			this.edges.length == face.edges.length &&
			NLA.arrayRange(0, edgeCount, 1)
				.some(offset => this.edges.every((edge, i) => edge.equals(face.edges[(offset + i) % edgeCount])))
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
			&& loopsLike(this.edges, face2.edges)
			&& this.holes.every(hole => face2.holes.some(hole2 => loopsLike(hole, hole2)))
	}

	getAllEdges() {
		return Array.prototype.concat.apply(this.edges, this.holes)
	}

	addEdgeLines(mesh) {
		assert(false, "buggy, fix")
		var vertices = this.edges.map(edge => edge.getVerticesNo0()).concatenated(), mvl = mesh.vertices.length
		for (var i = 0; i < vertices.length; i++) {
			mesh.vertices.push(vertices[i])
			mesh.lines.push(mvl + i, mvl + (i + 1) % vertices.length)

		}
	}

    inB2(): B2 {
        return new B2([this])
    }

    containsPoint(p: V3): boolean {
        assertVectors(p)
        return this.surface.edgeLoopContainsPoint(this.edges, p)
            && !this.holes.some(hole => this.surface.edgeLoopContainsPoint(hole, p))
    }

	/**
	 *
	 * @param line
	 * @returns t param of the line if there is an intersection, NaN otherwise
	 */
	intersectsLine(line) {
		assertInst(L3, line)
		let containedIntersectionsTs = this.surface.isTsForLine(line).filter(t => this.containsPoint(line.at(t)))
		let nearestPointT = containedIntersectionsTs.withMax(t => -t)

		return nearestPointT || NaN
	}

	toMesh() {
		let mesh = new GL.Mesh({triangles: true, normals: true, lines: true})
		mesh.faceIndexes = new Map()
		this.addToMesh(mesh)
		mesh.compile()
		return mesh
	}
}

class PlaneFace extends Face {
    surface: PlaneSurface

	constructor(planeSurface, contour, holes, name) {
		assertInst(PlaneSurface, planeSurface)
		super(planeSurface, contour, holes, name)
	}

	calculateArea() {
		assert(false)
	}

	addToMesh(mesh) {
		var mvl = mesh.vertices.length
		var normal = this.surface.plane.normal
		var vertices = this.edges.flatMap(edge => edge.getVerticesNo0())
		for (var i = 0; i < vertices.length; i++) { mesh.lines.push(mvl + i, mvl + (i + 1) % vertices.length) }
		var holeStarts = []
		this.holes.forEach(hole => {
			holeStarts.push(vertices.length)
			vertices.pushAll(hole.flatMap(edge => edge.getVerticesNo0()))
		})
		var triangles = triangulateVertices(normal, vertices, holeStarts).map(index => index + mvl)
		mesh.faceIndexes && mesh.faceIndexes.set(this, {start: mesh.triangles.length, count: triangles.length})
		Array.prototype.push.apply(mesh.vertices, vertices)
		Array.prototype.push.apply(mesh.triangles, triangles)
		Array.prototype.push.apply(mesh.normals, NLA.arrayFromFunction(vertices.length, i => normal))
	}

	containsPoint(p) {
		assertVectors (p)
		return this.surface.edgeLoopContainsPoint(this.edges, p)
			&& !this.holes.some(hole => this.surface.edgeLoopContainsPoint(hole, p))
	}

	intersectsLine(line) {
		assertInst(L3, line)
		var lambda = line.intersectWithPlaneLambda(this.surface.plane)
		if (!Number.isFinite(lambda)) {
			return NaN
		}
		var inside = this.containsPoint(line.at(lambda))
		return inside ? lambda : NaN
	}

	withHole(holeEdges) {
		return new PlaneFace(this.surface, this.edges, [holeEdges])
	}

	intersectFace(face2, thisBrep, face2Brep, faceMap, edgeMap, likeSurfaceFaces) {
		if (face2 instanceof PlaneFace) {
			this.intersectPlaneFace(face2, thisBrep, face2Brep, faceMap, edgeMap, likeSurfaceFaces)
			return
		}
		if (face2 instanceof RotationFace) {
			face2.intersectFace(this, face2Brep, thisBrep, faceMap, edgeMap, likeSurfaceFaces)
			return
		}
		assert(false)
	}

	intersectPlaneFace(face2: PlaneFace,
                       thisBrep: B2,
                       face2Brep: B2,
                       faceMap: Map<Face, Edge[]>,
                       edgeMap: Map<Edge, {edgeT: number, p: V3}[]>,
                       likeSurfaceFaces) {
        assertInst(PlaneFace, face2)
		let face:PlaneFace = this
		// get intersection
		let thisPlane = this.surface.plane, face2Plane = face2.surface.plane
		if (thisPlane.isParallelToPlane(face2Plane)) {
			if (thisPlane.like(face2Plane)) {
				// normal same and same location in space
				addLikeSurfaceFaces(likeSurfaceFaces, this, face2)
			}
			return
		}
		let intersectionLine = L3.fromPlanes(thisPlane, face2Plane)
		let thisDir = true
		// get intersections of newCurve with other edges of face and face2
		let ps1 = planeFaceEdgeISPsWithPlane(thisBrep, this, intersectionLine, face2Plane)
		let ps2 = planeFaceEdgeISPsWithPlane(face2Brep, face2, intersectionLine, thisPlane)
		console.log('ps1\n', ps1.map(m => m.toSource()).join('\n'), '\nps2\n', ps2.map(m => m.toSource()).join('\n'))
		if (ps1.length == 0 || ps2.length == 0) {
			// faces to not intersect
			return
		}
		console.log(''+thisPlane+face2Plane)

		/**
		 *
		 * @param seg generated segment
		 * @param col1 if seg is colinear to an edge of this, the edge in question
		 * @param col2 same for face2
		 * @param in1
		 * @param in2
		 * @param a
		 * @param b
		 */
		function handleGeneratedSegment(seg:StraightEdge, col1, col2, in1, in2, a, b) {
			console.log("handling seg", col1 && col1.toSource(), col2 && col2.toSource(), seg.toSource())
			if (!col1 && !col2) {
				console.log("adding")
				NLA.mapAdd(faceMap, face, seg)
				NLA.mapAdd(faceMap, face2, seg.flipped())
			}
			if (!col1 && col2) {
				console.log('!col1 && col2')
				if (col2.aDir.cross(face2Plane.normal).dot(thisPlane.normal) > 0) {
					// NB: a new edge is inserted even though it may be the same as an old one
					// however it indicates that it intersects the other volume here, i.e. the old edge cannot
					// be counted as "inside" for purposes of reconstitution
					NLA.mapAdd(faceMap, face2, seg.flipped())

					NLA.mapAdd(edgeMap, col2, {edgeT: col2.curve.pointLambda(seg.a), p: seg.a})
					NLA.mapAdd(edgeMap, col2, {edgeT: col2.curve.pointLambda(seg.b), p: seg.b})
					let testVector = face2Plane.normal.negated().rejectedFrom(thisPlane.normal)
					console.log("testVector", testVector.sce, face2Brep.sce, splitsVolumeEnclosingFaces(face2Brep, col2, testVector, thisPlane.normal) == INSIDE)
					if (splitsVolumeEnclosingFaces(face2Brep, col2, testVector, thisPlane.normal) == INSIDE) {
						NLA.mapAdd(faceMap, face, seg)
					}
				}
			}
			if (!col2 && col1) {
				if (col1.aDir.cross(thisPlane.normal).dot(face2Plane.normal) > 0) {
					NLA.mapAdd(faceMap, face, seg)

					NLA.mapAdd(edgeMap, col1, {edgeT: col1.curve.pointLambda(seg.a), p: seg.a})
					NLA.mapAdd(edgeMap, col1, {edgeT: col1.curve.pointLambda(seg.b), p: seg.b})
					let testVector = thisPlane.normal.negated().rejectedFrom(face2Plane.normal)
					if (splitsVolumeEnclosingFaces(thisBrep, col1, testVector, face2Plane.normal) == INSIDE) {
						NLA.mapAdd(faceMap, face2, seg.flipped())
					}
				}
			}
			if (col1 && col2) {
				// faces are touching edge-edge
				let testVector1 = col1.aDir.cross(thisPlane.normal).negated()
				let testVector2 = col2.aDir.cross(face2Plane.normal).negated()
				console.log("testVector", testVector1.sce, face2Brep.sce, splitsVolumeEnclosingFaces(face2Brep, col2, testVector1, thisPlane.normal) == INSIDE)
				if (splitsVolumeEnclosingFaces(face2Brep, col2, testVector1, thisPlane.normal) == INSIDE
					&& splitsVolumeEnclosingFaces(thisBrep, col1, testVector2, face2Plane.normal) == INSIDE) {
					//NLA.mapAdd(faceMap, face2, seg.flipped())
					NLA.mapAdd(edgeMap, col1, {edgeT: col1.curve.pointLambda(seg.a), p: seg.a})
					NLA.mapAdd(edgeMap, col1, {edgeT: col1.curve.pointLambda(seg.b), p: seg.b})

					NLA.mapAdd(edgeMap, col2, {edgeT: col2.curve.pointLambda(seg.a), p: seg.a})
					NLA.mapAdd(edgeMap, col2, {edgeT: col2.curve.pointLambda(seg.b), p: seg.b})
				}
			}
		}

		/**
		 * What does this do??
		 * @param a
		 * @param b
		 * @param useA
		 * @param useB
		 */
		function handlePoint(a, b, useA, useB) {
			console.log("logging point", useA, a.toSource())
			console.log("logging point", useB, b.toSource())
			if (useA && !useB) {
				if (!a.colinear) {
					NLA.mapAdd(edgeMap, a.edge, a)
				}
				// else colinear segment ends in middle of other face, do nothing
			}
			if (useB && !useA) {
				if (!b.colinear) {
					NLA.mapAdd(edgeMap, b.edge, b)
				}
				// else colinear segment ends in middle of other face, do nothing
			}
			if (useA && useB) {
				// segment start/ends on edge/edge intersection
				if (!a.colinear && !NLA.eq(a.edgeT, a.edge.aT) && !NLA.eq(a.edgeT, a.edge.bT)) {
					//						assert(a.edgeT != a.edge.aT && a.edgeT != a.edge.bT, "implement point intersections")
					// ends on a, on colinear segment b bT != a.edge.bT &&
					let testVector = a.edge.aDir.rejectedFrom(b.edge.aDir)
					let sVEF1 = splitsVolumeEnclosingFaces(face2Brep, b.edge, testVector, thisPlane.normal)
					let sVEF2 = splitsVolumeEnclosingFaces(face2Brep, b.edge, testVector.negated(), thisPlane.normal)
					if (!(INSIDE == sVEF1 && INSIDE == sVEF2
						|| (OUTSIDE == sVEF1 || COPLANAR_OPPOSITE == sVEF1) && (OUTSIDE == sVEF2 || COPLANAR_OPPOSITE == sVEF2))) {
						NLA.mapAdd(edgeMap, a.edge, a)
					}
				}
				if (!b.colinear && !NLA.eq(b.edgeT, b.edge.aT) && !NLA.eq(b.edgeT, b.edge.bT)) {
					//						assert(b.edgeT != b.edge.aT && b.edgeT != b.edge.bT, "implement point intersections")
					// ends on b, on colinear segment a
					let testVector = b.edge.aDir.rejectedFrom(a.edge.aDir)
					let sVEF1 = splitsVolumeEnclosingFaces(thisBrep, a.edge, testVector, face2Plane.normal)
					let sVEF2 = splitsVolumeEnclosingFaces(thisBrep, a.edge, testVector.negated(), face2Plane.normal)
					if (!(INSIDE == sVEF1 && INSIDE == sVEF2
						|| (OUTSIDE == sVEF1 || COPLANAR_OPPOSITE == sVEF1) && (OUTSIDE == sVEF2 || COPLANAR_OPPOSITE == sVEF2))) {
						NLA.mapAdd(edgeMap, b.edge, b)
					}
				}
			}
		}

		let in1 = false, in2 = false, colinear1 = false, colinear2 = false
		let i = 0, j = 0, last, segments = []
		let startP, startDir, startT
		while (i < ps1.length || j < ps2.length) {
			assert(i <= ps1.length)
			assert(j <= ps2.length)
			let a = ps1[i], b = ps2[j]
			assert(a || b)
			if (j == ps2.length || i < ps1.length && NLA.lt(a.t, b.t)) {
				last = a
				in1 = !in1
				// ": colinear1" to remember the colinear segment, as segments are only handled once it ends (in1 false)
				colinear1 = in1 ? a.colinear && a.edge : colinear1
				i++
			} else if (i == ps1.length || NLA.gt(a.t, b.t)) {
				last = b
				in2 = !in2
				colinear2 = in2 ? b.colinear && b.edge : colinear2
				j++
			} else {
				last = a
				in1 = !in1
				in2 = !in2
				//if (in1 == in2) {
				a.used = true
				b.used = true
				colinear1 = in1 ? a.colinear && a.edge : colinear1
				colinear2 = in2 ? b.colinear && b.edge : colinear2
				//}
				i++
				j++
			}
			console.log("as", a, b, in1, in2)
			if (startP && !(in1 && in2)) {
				// segment end
				let newEdge = new StraightEdge(intersectionLine, startP, last.p, startT, last.t, null, 'genseg' + globalId++)
				handleGeneratedSegment(newEdge, colinear1, colinear2, in1, in2, a, b)
				startP = undefined
				last.used = true
				handlePoint(a, b, colinear1 || !!a.used, colinear2 || !!b.used)
			} else if (in1 && in2) {
				// new segment just started
				startP = last.p
				startDir = last.insideDir
				startT = last.t
				last.used = true
				handlePoint(a, b, colinear1 || !!a.used, colinear2 || !!b.used)
			}
		}
	}


    static forVertices(planeSurface, vs, ...holeVs) {
        if (planeSurface instanceof P3) {
            planeSurface = new PlaneSurface(planeSurface)
        }
        assert(isCCW(vs, planeSurface.plane.normal), 'isCCW(vs, planeSurface.plane.normal)')
        var edges = vs.map((a, i, vs) => {
            var b = vs[(i + 1) % vs.length]
            return StraightEdge.throughPoints(a, b)
        })
        holeVs.forEach(vs => assert(doubleSignedArea(vs, planeSurface.plane.normal) >= 0, 'doubleSignedArea(vs, planeSurface.plane.normal) >= 0'))
        var holes = holeVs.map(hvs => hvs.map((a, i, vs) => {
            var b = vs[(i + 1) % vs.length]
            return StraightEdge.throughPoints(a, b)
        }))
        return new PlaneFace(planeSurface, edges, holes)
    }
}
NLA.registerClass(PlaneFace)



class RotationFace extends Face {
	constructor(rot, contour, holes, name) {
		super(rot, contour, holes, name)
	}

	unrollLoop(edgeLoop) {
		var vs = []
		var reverseFunc = this.surface.pointToParameterFunction()
		let verticesNo0s = edgeLoop.map(edge => edge.getVerticesNo0())
		let startEdgeIndex = verticesNo0s.findIndex(edgeVertices => !NLA.eq(reverseFunc(edgeVertices[0], Math.PI).x, Math.PI))
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
		// console.log("vs\n", vs.join("\n"), vs.length)
		return vs
	}

	addToMesh(mesh) {
		var closed = false
		var hSplit = 32, zSplit = 1
		var ribs = []
		var minZ = Infinity, maxZ = -Infinity
		var cmp = (a, b) => a.value - b.value
		var f = this.surface.parametricFunction()
		var normalF = this.surface.parametricNormal()
		var vertexLoops = this.holes.concat([this.edges]).map(loop => this.unrollLoop(loop))
		vertexLoops.forEach(vertexLoop => {
			vertexLoop.forEach(({x: d, y: z}) => {
				var index0 = ribs.binaryIndexOf(d, (a, b) => NLA.snapTo(a.value - b, 0))
				if (index0 < 0) {
					ribs.splice(-index0-1, 0, {value: d, left: [], right: []})
				}
				minZ = min(minZ, z)
				maxZ = max(maxZ, z)
			})
		})
		var correction = 1
		vertexLoops.forEach(vertexLoop => {
			vertexLoop.forEach((v0, i, vs) => {
				var v1 = vs[(i + 1) % vs.length], dDiff = v1.x - v0.x
				//console.log(v0.sce, v1.sce)
				if (NLA.eq0(dDiff)) {
					return
				}
				if (dDiff < 0) {
					[v0, v1] = [v1, v0]
					dDiff = -dDiff
				}
				var index0 = ribs.binaryIndexOf(v0.x, (a, b) => NLA.snapTo(a.value - b, 0))
				var index1 = ribs.binaryIndexOf(v1.x, (a, b) => NLA.snapTo(a.value - b, 0))
				ribs[index0].right.binaryInsert(v0.y)
				for (var j = (index0 + correction) % ribs.length; j != index1; j = (j + correction) % ribs.length) {
					var x = ribs[j].value
					var part = (x - v0.x) / dDiff
					var interpolated = v1.y * part + v0.y * (1 - part)
					ribs[j].left.binaryInsert(interpolated)
					ribs[j].right.binaryInsert(interpolated)
				}
				ribs[index1].left.binaryInsert(v1.y)
				// console.log(ribs.map(r=>r.toSource()).join('\n'))
			})
		})
		var vertices = [], triangles = [], normals = []
		for (let i = 0; i < ribs.length; i++) {
			let ribLeft = ribs[i], ribRight = ribs[(i + 1) % ribs.length]
			assert(ribLeft.right.length == ribRight.left.length)
			for (let j = 0; j < ribLeft.right.length; j++) {
				vertices.push(f(ribLeft.value, ribLeft.right[j]), f(ribRight.value, ribRight.left[j]))
				normals.push(normalF(ribLeft.value, ribLeft.right[j]), normalF(ribRight.value, ribRight.left[j]))
			}
		}
		//console.log(ribs.map(r=>r.toSource()).join('\n'))
		var vss = vertices.length, detailVerticesStart = vss
		var zInterval = maxZ - minZ, zStep = zInterval / zSplit
		var detailZs = NLA.arrayFromFunction(zSplit - 1, i => minZ + (1 + i) * zStep)
		for (let i = 0; i < ribs.length; i++) {
			var d = ribs[i].value
			for (let j = 0; j < detailZs.length; j++) {
				vertices.push(f(d, detailZs[j]))
				normals.push(normalF(d, detailZs[j]))
			}
		}
		// console.log('detailVerticesStart', detailVerticesStart, 'vl', vertices.length, vertices.length - detailVerticesStart, ribs.length)
		// finally, fill in the ribs
		var vsStart = 0
		var flipped2 = true
		//for (var i = 0; i < 1; i++) {
		var end = closed ? ribs.length : ribs.length - 1
		for (let i = 0; i < end; i++) {
			var ipp = (i + 1) % ribs.length
			let inside = false, colPos = 0, ribLeft = ribs[i], ribRight = ribs[(i + 1) % ribs.length]
			for (let j = 0; j < detailZs.length + 1; j++) {
				var detailZ = detailZs[j] || 100000
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
		//assert(normals.every(n => n.hasLength(1)), normals.find(n => !n.hasLength(1)).length() +" "+normals.findIndex(n => !n.hasLength(1)))
		mesh.faceIndexes.set(this, {start: mesh.triangles.length, count: triangles.length})
		Array.prototype.push.apply(mesh.vertices, vertices)
		Array.prototype.push.apply(mesh.triangles, triangles)
		Array.prototype.push.apply(mesh.normals, normals)
		//this.addEdgeLines(mesh)

	}

	intersectFace(face2:Face, thisBrep:B2, face2Brep:B2, faceMap:Map<Face, Edge[]>, edgeMap:Map<Edge, any[]>, likeSurfaceFaces) {

		/**
		 * @param seg generated segment
		 * @param col1 if seg is colinear to an edge of this, the edge in question
		 * @param col2 same for face2
		 * @param in1
		 * @param in2
		 * @param a
		 * @param b
		 */
		function handleGeneratedSegment(seg, col1, col2, in1, in2, a, b) {
			console.log("handling seg", col1 && col1.toSource(), col2 && col2.toSource(), seg.toSource())
			if (!col1 && !col2) {
				console.log("adding")
				NLA.mapAdd(faceMap, face, seg)
				NLA.mapAdd(faceMap, face2, seg.flipped())
			}
			if (!col1 && col2) {
				console.log('!col1 && col2')
				if (col2.aDir.cross(surface2.normal).dot(surface.normal) > 0) {
					// NB: a new edge is inserted even though it may be the same as an old one
					// however it indicates that it intersects the other volume here, i.e. the old edge cannot
					// be counted as "inside" for purposes of reconstitution
					NLA.mapAdd(faceMap, face2, seg.flipped())

					NLA.mapAdd(edgeMap, col2, {edgeT: col2.curve.pointLambda(seg.a), p: seg.a})
					NLA.mapAdd(edgeMap, col2, {edgeT: col2.curve.pointLambda(seg.b), p: seg.b})
					let testVector = surface2.normalAt(col2.a).negated().rejectedFrom(surface.normalAt(col2.a))
					console.log("testVector", testVector.sce, face2Brep.sce, splitsVolumeEnclosingFaces(face2Brep, col2, testVector, surface.normal) == INSIDE)
					if (splitsVolumeEnclosingFaces(face2Brep, col2, testVector, surface.normal) == INSIDE) {
						NLA.mapAdd(faceMap, face, seg)
					}
				}
			}
			if (!col2 && col1) {
				if (col1.aDir.cross(surface.normal).dot(surface2.normal) > 0) {
					NLA.mapAdd(faceMap, face, seg)

					NLA.mapAdd(edgeMap, col1, {edgeT: col1.curve.pointLambda(seg.a), p: seg.a})
					NLA.mapAdd(edgeMap, col1, {edgeT: col1.curve.pointLambda(seg.b), p: seg.b})
					let testVector = surface.normal.negated().rejectedFrom(surface2.normal)
					if (splitsVolumeEnclosingFaces(thisBrep, col1, testVector, surface2.normal) == INSIDE) {
						NLA.mapAdd(faceMap, face2, seg.flipped())
					}
				}
			}
			if (col1 && col2) {
				// faces are touching edge-edge
				let testVector1 = col1.aDir.cross(surface.normal).negated()
				let testVector2 = col2.aDir.cross(surface2.normal).negated()
				console.log("testVector", testVector1.sce, face2Brep.sce, splitsVolumeEnclosingFaces(face2Brep, col2, testVector1, surface.normal) == INSIDE)
				if (splitsVolumeEnclosingFaces(face2Brep, col2, testVector1, surface.normal) == INSIDE
					&& splitsVolumeEnclosingFaces(thisBrep, col1, testVector2, surface2.normal) == INSIDE) {
					//NLA.mapAdd(faceMap, face2, seg.flipped())
					NLA.mapAdd(edgeMap, col1, {edgeT: col1.curve.pointLambda(seg.a), p: seg.a})
					NLA.mapAdd(edgeMap, col1, {edgeT: col1.curve.pointLambda(seg.b), p: seg.b})

					NLA.mapAdd(edgeMap, col2, {edgeT: col2.curve.pointLambda(seg.a), p: seg.a})
					NLA.mapAdd(edgeMap, col2, {edgeT: col2.curve.pointLambda(seg.b), p: seg.b})
				}
			}
		}

		/**
		 * What does this do?? todo
		 * @param a
		 * @param b
		 * @param useA
		 * @param useB
		 */
		function handlePoint(a, b, useA, useB) {
			console.log("logging point", useA, a && a.toSource())
			console.log("logging point", useB, b && b.toSource())
			if (useA && !useB) {
				if (!a.colinear) {
					NLA.mapAdd(edgeMap, a.edge, a)
				}
				// else colinear segment ends in middle of other face, do nothing
			}
			if (useB && !useA) {
				if (!b.colinear) {
					NLA.mapAdd(edgeMap, b.edge, b)
				}
				// else colinear segment ends in middle of other face, do nothing
			}
			if (useA && useB) {
				// segment start/ends on edge/edge intersection
				if (!a.colinear && !NLA.eq(a.edgeT, a.edge.aT) && !NLA.eq(a.edgeT, a.edge.bT)) {
					//						assert(a.edgeT != a.edge.aT && a.edgeT != a.edge.bT, "implement point intersections")
					// ends on a, on colinear segment b bT != a.edge.bT &&
					let testVector = a.edge.aDir.rejectedFrom(b.edge.aDir)
					let sVEF1 = splitsVolumeEnclosingFaces(face2Brep, b.edge, testVector, surface.normalAt(b.edge.a))
					let sVEF2 = splitsVolumeEnclosingFaces(face2Brep, b.edge, testVector.negated(), surface.normalAt(b.edge.a))
					if (!(INSIDE == sVEF1 && INSIDE == sVEF2
						|| (OUTSIDE == sVEF1 || COPLANAR_OPPOSITE == sVEF1) && (OUTSIDE == sVEF2 || COPLANAR_OPPOSITE == sVEF2))) {
						NLA.mapAdd(edgeMap, a.edge, a)
					}
				}
				if (!b.colinear && !NLA.eq(b.edgeT, b.edge.aT) && !NLA.eq(b.edgeT, b.edge.bT)) {
					//						assert(b.edgeT != b.edge.aT && b.edgeT != b.edge.bT, "implement point intersections")
					// ends on b, on colinear segment a
					let testVector = b.edge.aDir.rejectedFrom(a.edge.aDir)
					let sVEF1 = splitsVolumeEnclosingFaces(thisBrep, a.edge, testVector, surface2.normal)
					let sVEF2 = splitsVolumeEnclosingFaces(thisBrep, a.edge, testVector.negated(), surface2.normal)
					if (!(INSIDE == sVEF1 && INSIDE == sVEF2
						|| (OUTSIDE == sVEF1 || COPLANAR_OPPOSITE == sVEF1) && (OUTSIDE == sVEF2 || COPLANAR_OPPOSITE == sVEF2))) {
						NLA.mapAdd(edgeMap, b.edge, b)
					}
				}
			}
		}


		assertInst(Face, face2)

		const face = this
		const surface = face.surface, surface2 = face2.surface
		if (surface.isCoplanarTo(surface2)) {
			addLikeSurfaceFaces(likeSurfaceFaces, face, face2)
			return
		}
		const isCurves = surface.isCurvesWithSurface(surface2)
		// get intersections of newCurve with other edges of face and face2
		const pss1 = faceEdgeISPsWithSurface(thisBrep, face, isCurves, face2.surface)
		const pss2 = faceEdgeISPsWithSurface(face2Brep, face2, isCurves, face.surface)
		console.log('pss1\n', pss1.map(m => m.toSource()).join('\n'), '\npss2\n', pss2.map(m => m.toSource()).join('\n'))
		console.log(''+surface+surface2)

		isCurves.forEach((isCurve, isCurveIndex) => {
			const ps1 = pss1[isCurveIndex], ps2 = pss2[isCurveIndex]
			// for non-endless curves, e.g. ellipses, the intersections of the faces can be non-zero, even if one of
			// the faces doesn't register any points on the curve. For example, if a cylinder is cut entirely by a
			// plane face (all its edges around the cylinder), then the face will contain the entire curve and
			// "ps" for the plane face will be empty
			const curvePoint = isCurve.at(0)
			// TODO: behavior when curves touch face?
			const faceContainsCurvePoint = face.containsPoint(curvePoint)
			const face2ContainsCurvePoint = face2.containsPoint(curvePoint)
			if (ps1.length == 0 && !faceContainsCurvePoint || ps2.length == 0 && !face2ContainsCurvePoint) {
				return
			}
			// !! start in does depend on insidedir... TODO
			assertf(() => (0 == ps1.length) || !NLA.eq0(ps1[0].insideDir.dot(isCurve.tangentAt(ps1[0].t))), () => ps1[0].insideDir.dot(isCurve.tangentAt(ps1[0].t)))
			assertf(() => (0 == ps2.length) || !NLA.eq0(ps2[0].insideDir.dot(isCurve.tangentAt(ps2[0].t))), () => ps2[0].insideDir.dot(isCurve.tangentAt(ps2[0].t)))
			// if ps is empty, in must be true or we wouldn't have gotten this far
			let in1 = (0 == ps1.length) || ps1[0].insideDir.dot(isCurve.tangentAt(ps1[0].t)) < 0
			let in2 = (0 == ps2.length) || ps2[0].insideDir.dot(isCurve.tangentAt(ps2[0].t)) < 0
			let colinear1:boolean|Edge = false, colinear2:boolean|Edge = false
			let i = 0, j = 0, last, segments = []
			let startP, startDir, startT
			while (i < ps1.length || j < ps2.length) {
				let a = ps1[i], b = ps2[j]
				if (j >= ps2.length || i < ps1.length && NLA.lt(a.t, b.t)) {
					last = a
					in1 = !in1
					// ": colinear1" to remember the colinear segment, as segments are only handled once it ends (in1 false)
					colinear1 = in1 ? a.colinear && a.edge : colinear1
					i++
				} else if (i >= ps1.length || NLA.gt(a.t, b.t)) {
					last = b
					in2 = !in2
					colinear2 = in2 ? b.colinear && b.edge : colinear2
					j++
				} else {
					last = a
					in1 = !in1
					in2 = !in2
					//if (in1 == in2) {
					a.used = true
					b.used = true
					colinear1 = in1 ? a.colinear && a.edge : colinear1
					colinear2 = in2 ? b.colinear && b.edge : colinear2
					//}
					i++
					j++
				}
				console.log("as", a, b, in1, in2)
				if (startP && !(in1 && in2)) {
					// segment end
					handleGeneratedSegment(Edge.create(isCurve, startP, last.p, startT, last.t, null, isCurve.tangentAt(startT), isCurve.tangentAt(last.t)), colinear1, colinear2, in1, in2, a, b)
					startP = undefined
					last.used = true
					handlePoint(a, b, colinear1 || a && !!a.used, colinear2 || b && !!b.used)
				} else if (in1 && in2) {
					// new segment just started
					startP = last.p
					startDir = last.insideDir
					startT = last.t
					last.used = true
					handlePoint(a, b, colinear1 || a && !!a.used, colinear2 || b && !!b.used)
				}
			}
		})
	}


}

function addLikeSurfaceFaces(likeSurfaceFaces, face1, face2) {
	// There cannot be two subgroups which will later be connected, as the "graph" of like surface faces is fully connected
	for (let i = 0; i < likeSurfaceFaces.length; i++) {
		let faceGroup = likeSurfaceFaces[i]
		let foundFace1 = false, foundFace2 = false
		for (let j = 0; j < faceGroup.length; j++) {
			let face = faceGroup[j]
			if (face == face1) {
				foundFace1 = true
			}
			if (face == face2) {
				foundFace2 = true
			}
		}
		if (foundFace1 != foundFace2) {
			faceGroup.push(foundFace1 ? face2 : face1)
			return
		} else if (foundFace1) {
			// found both
			return
		}
	}
	// nothing found, add a new group
	likeSurfaceFaces.push([face1, face2])
}

/**
 *
 */
class B2 extends Transformable {
	faces: Face[]
	infiniteVolume: boolean
	generator: string
	vertexNames: Map<V3, string>

	constructor(faces, infiniteVolume, generator, vertexNames) {
		super()
		this.faces = faces
		assertInst.apply(undefined, [Face].concat(faces))
		this.infiniteVolume = !!infiniteVolume
		this.generator = generator
		this.vertexNames = vertexNames
	}

	calculateVolume(): number {
		return this.faces.map(face => face.zDirVolume()).sum()
	}

	toMesh(): GL.Mesh {
		var mesh = new GL.Mesh({triangles: true, normals: true, lines: true})
		mesh.faceIndexes = new Map()
		this.faces.forEach((face, i) => {
			face.addToMesh(mesh)
		})
		mesh.compile()
		return mesh
	}

	minus(brep2: B2): B2 {
		return this.intersection(brep2.flipped(), true, true)
	}

	plus(brep2: B2): B2 {
		let result = this.flipped().intersection(brep2.flipped(), true, true).flipped()
		result.generator = `${this.generator}.plus(${brep2.generator})`
		return result
	}

	xor(brep2: B2): B2 {
		assertNever()
	}

	equals(brep): boolean {
		return this.faces.length == brep.faces.length &&
			this.faces.every((face) => brep.faces.some((face2) => face.equals(face2)))
	}

	toString(): string {
		return `new B2([\n${this.faces.join(',\n').replace(/^/gm, '\t')}])`
	}

	toSource(): string {
		return this.generator || `new B2([\n${this.faces.map(face => face.toSource()).join(',\n').replace(/^/gm, '\t')}])`
	}

	assembleFacesFromLoops(newFaces: Face[], loops, oldFace): boolean {
		function placeRecursively(newLoop, arr) {
			if (arr.length == 0) {
				arr.push(newLoop)
			} else {
				var sl = arr.find(subloop => {
					var contains = surface.edgeLoopContainsPoint(subloop.edges, newLoop.edges[0].a)
					console.log(newLoop.edges[0].a.sce, newLoop.ccw ? contains : !contains)
					return contains
				})
				console.log("here", sl)
				if (sl) {
					placeRecursively(newLoop, sl.subloops)
				} else {
					// newLoop isnt contained by any other subloop
					for (var i = arr.length; --i >= 0;) {
						var subloop = arr[i]
						//console.log("cheving subloop", surface.edgeLoopContainsPoint(newLoop.edges, subloop.edges[0].a))
						if (surface.edgeLoopContainsPoint(newLoop.edges, subloop.edges[0].a)) {
							newLoop.subloops.push(subloop)
							arr.splice(i, 1) // remove it
						}
					}
					arr.push(newLoop)
				}
			}
		}

		function newFacesRecursive(loop) {
			var face = new oldFace.constructor(oldFace.surface, loop.edges, loop.subloops.map(sl => sl.edges))
			loop.subloops.forEach(sl => sl.subloops.forEach(tlLoop => newFacesRecursive(tlLoop)))
			console.log(face.toString(), face.holes.toString())
			newFaces.push(face)
		}

		var surface = oldFace.surface
		console.log(loops.map(loop=> loop.join('\n')).join('\n\n'))
		var topLevelLoops = []
		loops.forEach(loop => {
			var ccw = surface.edgeLoopCCW(loop)
			console.log('CCW', ccw)
			placeRecursively({edges: loop, ccw: ccw, subloops: []}, topLevelLoops)
			console.log(topLevelLoops)
		})
		if (topLevelLoops[0].ccw) {
			let firstSubLoop = topLevelLoops[0].subloops[0]
			if (firstSubLoop && firstSubLoop.ccw) {
				newFacesRecursive(firstSubLoop)
			} else {
				topLevelLoops.forEach(tlLoop => newFacesRecursive(tlLoop))
			}
			return false
		} else {
			newFacesRecursive({edges: oldFace.edges, ccw: true, subloops: topLevelLoops})
			return true
		}
	}

	/**
	 * Rightmost next segment doesn't work, as the correct next segment isn't obvious from the current corner
	 * alone.
	 * (at least, not without extensive pre-analysos on the face edges, which shouldn't be necessary, as the
	 * correct new faces are defined by the new edges already.) Leftmost edge should work. Holes which touch the
	 * edge of the face will be added to the face contour.
	 *
	 * New segments will always be part left-er than exisiting ones, so no special check is required.
	 *
	 */
	reconstituteFaces(oldFaces, edgeLooseSegments, faceMap, newFaces) {
		// reconstitute faces
		var insideEdges = []
		oldFaces.forEach(face => {
			var allOldFaceEdges = face.getAllEdges()
			console.log('reconstituting face', face.toString())
			var faceLooseSegments = face.getAllEdges().map(edge => edgeLooseSegments.get(edge)).filter(x=>x).concatenated()
			var els = face.edges.map(edge => edgeLooseSegments.get(edge) && edgeLooseSegments.get(edge).join('\n')).map(s => '\n\n' + s).join()
			console.log('edgeLooseSegments', els)
			var newSegments = faceMap.get(face)
			if (!newSegments) {
				face.insideOutside = 'undecided'
			} else {
				face.insideOutside = 'part'
				console.log('newSegments\n', newSegments.map(e=>e.toString()).join('\n'))
				var loops = []
				var currentEdge
				var edgeCond = face instanceof PlaneFace
					? (edge => edge.a.like(currentEdge.b))
					: (edge => (edge.curve == currentEdge.curve // TODO: ??
					? NLA.eq(edge.aT, currentEdge.bT)
					: edge.a.like(currentEdge.b)))
				var getNextStart = function () {
					return newSegments.find(edge => !edge.visited)
				}
				allOldFaceEdges.forEach(edge => edge.visited = false)
				while (currentEdge = getNextStart()) {
					var cancelLoop = false
					var startEdge = currentEdge, edges = [], i = 0
					// wether only new edges are used (can include looseSegments)
					var looseLoop = true
					do {
						currentEdge.visited = true
						console.log('currentEdge', currentEdge.b.sce, currentEdge.toSource())
						edges.push(currentEdge)
						// find next edge
						var possibleLooseEdges = newSegments.filter(edge => edgeCond(edge)), looseSegments
						var possibleNewEdgesCount = possibleLooseEdges.length
						possibleLooseEdges.pushAll(faceLooseSegments.filter(edgeCond))
						var noSegments = false
						if (possibleNewEdgesCount == possibleLooseEdges.length) {
							possibleLooseEdges.pushAll(allOldFaceEdges.filter(edgeCond))
							noSegments = true
						}
						var index = possibleLooseEdges.indexWithMax((edge, index) => currentEdge.bDir.angleRelativeNormal(edge.aDir, face.surface.normalAt(currentEdge.b)) + (index < possibleNewEdgesCount) * NLA_PRECISION)
						//console.log('possibleLooseEdges\n', possibleLooseEdges.map(e=>e.toString()).join('\n'), allOldFaceEdges.find(edge => edgeCond(edge)), index, possibleLooseEdges[index])
						// TODO assert(possibleLooseEdges.length < 2)
						assert(0 < possibleLooseEdges.length)
						currentEdge = possibleLooseEdges[index]
						if (currentEdge.visited) {
							console.log("breaking")
							break
						}
						if (index < possibleNewEdgesCount) {
							possibleLooseEdges.forEach(possibleLooseEdge => possibleLooseEdge.isCoEdge(currentEdge) && (possibleLooseEdge.visited = true))
						} else {
							looseLoop = false
							noSegments && insideEdges.push(currentEdge)
						}
						assert(currentEdge)
						assert(currentEdge != startEdge)
					} while (++i < 200)
					if (200 == i) {
						assert(false, "too many")
					}
					if (currentEdge == startEdge) {
						console.log('finished loop')
						loops.push(edges)
					} else {
						insideEdges.removeAll(edges)
					}
				}

				if (this.assembleFacesFromLoops(newFaces, loops, face)) {
					insideEdges.pushAll(face.edges)
				}
			}
		})
		while (insideEdges.length != 0) {
			var edge = insideEdges.pop()
			var adjoiningFaces = facesWithEdge(edge, oldFaces)
			adjoiningFaces.forEach(info => {
				if (info.face.insideOutside == 'undecided') {
					info.face.insideOutside = 'inside'
					insideEdges.push.apply(insideEdges, info.face.edges)
				}
			})
		}
		oldFaces.forEach(face => {
			if (face.insideOutside == 'inside') {
				newFaces.push(face)
			}
		})
	}

	reconstituteCoplanarFaces(likeSurfacePlanes, edgeLooseSegments, faceMap, newFaces) {
		likeSurfacePlanes.forEach(faceGroup => {
			// calculate total contours
			let surface = faceGroup[0].surface, bag = []
			faceGroup.forEach(face => {
				Array.prototype.push.apply(bag, faceMap(face))
				face.getAllEdges().forEach(edge => {
					let edgeSubSegments
					if (edgeSubSegments = edgeLooseSegments.get(edge)) {
						Array.prototype.push.apply(bag, edgeSubSegments)
					} else {
						bag.push(edge)
					}
				})
			})
			let currentEdge, loops = []
			while (currentEdge = bag.find(edge => !edge.visited)) {
				let path = []
				do {
					currentEdge.visited = true
					path.push(currentEdge)
					let possibleNextEdges = bag.filter(edge => currentEdge.b.like(edge.a))
					// lowest angle, i.e. the right-most next edge
					let nextEdgeIndex = possibleNextEdges.indexWithMax((edge, index) => -currentEdge.bDir.angleRelativeNormal(edge.aDir, surface.normalAt(currentEdge.b)))
					currentEdge = possibleNextEdges[nextEdgeIndex]
				} while (!currentEdge.visited)
				let startIndex = path.find(currentEdge)
				if (-1 != startIndex) {
					loops.push(path.slice(startIndex))
				}
			}
		})
	}

	getLooseEdgeSegmentsInsideDirs(edgeMap: Map<Edge, {edgeT: number, p: V3}[]>): Map<Edge, Edge> {
		var edgeLooseSegments = new Map()
		console.log("edgeMap", edgeMap)
		edgeMap.forEach((pointInfos, baseEdge) => {
			// TODO: make sure it works with loop
			// TODO: dont generate unnecessarry segments
			var looseSegments = []
			if (!baseEdge.reversed) {
				pointInfos.sort((a, b) => NLA.snapTo(a.edgeT - b.edgeT, 0) || (b.insideDir.dot(baseEdge.curve.dir1) - a.insideDir.dot(baseEdge.curve.dir1)))
			} else {
				pointInfos.sort((b, a) => NLA.snapTo(a.edgeT - b.edgeT, 0) || (b.insideDir.dot(baseEdge.curve.dir1) - a.insideDir.dot(baseEdge.curve.dir1)))
			}
			var startP = baseEdge.a, startDir = baseEdge.aDir, startT = baseEdge.aT, startInfo
			for (var i = 0; i < pointInfos.length; i++) {
				var info = pointInfos[i]
				assert(info.insideDir, info.toSource())
				console.log("info.insideDir.dot(baseEdge.curve.dir1)", info.insideDir.dot(baseEdge.curve.dir1))
				// ignore start, end and repeating points
				if (NLA.eq(info.edgeT, baseEdge.bT) || NLA.eq(info.edgeT, baseEdge.aT) || NLA.eq(info.edgeT, startT)) {
					continue
				}
				var pDir = baseEdge.tangentAt(info.edgeT)
				// add segment only if insideDir points backwards
				if (info.insideDir.dot(baseEdge.curve.dir1) < 0) {
					looseSegments.push(new baseEdge.constructor(baseEdge.curve, startP, info.p, startT, info.edgeT, null, startDir, pDir))
				}
				startP = info.p
				startT = info.edgeT
				startInfo = info
				startDir = pDir
			}
			if (startInfo && startInfo.insideDir.dot(baseEdge.curve.dir1) > 0) {
				looseSegments.push(new baseEdge.constructor(baseEdge.curve, startP, baseEdge.b, startT, baseEdge.bT, null, startDir, baseEdge.bDir))
			}
			edgeLooseSegments.set(baseEdge, looseSegments)
		})
		return edgeLooseSegments
	}

	getLooseEdgeSegments(edgeMap: Map<Edge, {edgeT: number, p: V3}[]>): Map<Edge, Edge[]> {
		var edgeLooseSegments = new Map()
		console.log("edgeMap", edgeMap)
		edgeMap.forEach((pointInfos, baseEdge) => {
			// TODO: make sure it works with loop
			// TODO: dont generate unnecessarry segments
			var looseSegments = []
			if (!baseEdge.reversed) {
				pointInfos.sort((a, b) => NLA.snapTo(a.edgeT - b.edgeT, 0) ||
				(b.insideDir.dot(baseEdge.curve.tangentAt(a.edgeT)) - a.insideDir.dot(baseEdge.curve.tangentAt(a.edgeT))))
			} else {
				pointInfos.sort((b, a) => NLA.snapTo(a.edgeT - b.edgeT, 0) ||
				(b.insideDir.dot(baseEdge.curve.tangentAt(a.edgeT)) - a.insideDir.dot(baseEdge.curve.tangentAt(a.edgeT))))
			}
			var startP = baseEdge.a, startDir = baseEdge.aDir, startT = baseEdge.aT, startInfo
			for (var i = 0; i < pointInfos.length; i++) {
				var info = pointInfos[i]
				// ignore start, end and repeating points
				if (NLA.eq(info.edgeT, baseEdge.bT) || NLA.eq(info.edgeT, baseEdge.aT) || NLA.eq(info.edgeT, startT)) {
					continue
				}
				var pDir = baseEdge.tangentAt(info.edgeT)
				looseSegments.push(Edge.create(baseEdge.curve, startP, info.p, startT, info.edgeT, null, startDir, pDir, 'looseSegment' + globalId++))
				startP = info.p
				startT = info.edgeT
				startInfo = info
				startDir = pDir
			}
			looseSegments.push(Edge.create(baseEdge.curve, startP, baseEdge.b, startT, baseEdge.bT, null, startDir, baseEdge.bDir, 'looseSegment' + globalId++))
			edgeLooseSegments.set(baseEdge, looseSegments)
		})
		return edgeLooseSegments;
	}

	getIntersectionEdges(brep2) {
		var faceMap = new Map(), edgeMap = new Map()

		let likeSurfaceFaces = []

		this.faces.forEach(face => {
			//console.log('face', face.toString())
			brep2.faces.forEach(face2 => {
				//console.log('face2', face2.toString())
				face.intersectFace(face2, this, brep2, faceMap, edgeMap, likeSurfaceFaces)
			})
		})

		return Array.from(faceMap.values()).concatenated()

	}

	intersection(brep2: B2, buildThis: boolean, buildBREP2: boolean, buildCoplanar): B2 {
		var faceMap = new Map(), edgeMap = new Map()

		let likeSurfaceFaces = []

		this.faces.forEach(face => {
			//console.log('face', face.toString())
			brep2.faces.forEach(face2 => {
				//console.log('face2', face2.toString())
				face.intersectFace(face2, this, brep2, faceMap, edgeMap, likeSurfaceFaces)
			})
		})
		var newFaces = []

		/*
		 TODO:
		 faceMap.forEach((faceLooses, face) => {
		 faceLooses.forEach(edge => {
		 face.edges.forEach(faceEdge => {
		 var edgeT = faceEdge.getEdgeT(edge.a)
		 if (undefined !== edgeT) {
		 console.log("WAARGH", edge.a.$, faceEdge.toString(), edgeT)
		 NLA.mapAdd(edgeMap, faceEdge, {edgeT: edgeT, p: edge.a})
		 }
		 })
		 })
		 })
		 */
		var edgeLooseSegments = this.getLooseEdgeSegments(edgeMap);

		buildThis && this.reconstituteFaces(this.faces, edgeLooseSegments, faceMap, newFaces, brep2.infiniteVolume)
		buildBREP2 && this.reconstituteFaces(brep2.faces, edgeLooseSegments, faceMap, newFaces, this.infiniteVolume)
		//buildCoplanar && this.reconstituteCoplanarFaces(likeSurfaceFaces, edgeLooseSegments, faceMap, newFaces, this.infiniteVolume, brep2.infiniteVolume)

		return new B2(newFaces, this.infiniteVolume && brep2.infiniteVolume)
	}

	transform(m4: M4, desc?: string): this {
		desc = desc || ''

		let vertexNames: Map<V3,string>
		if (this.vertexNames) {
			vertexNames = new Map()
			this.vertexNames.forEach((name, vertex) => vertexNames.set(m4.transformPoint(vertex), name + desc))
		}
		return new B2(
			this.faces.map(f => f.transform(m4)),
			this.infiniteVolume,
			this.generator && this.generator + desc,
			vertexNames
		) as this
	}

	flipped(): B2 {
		return new B2(
			this.faces.map(f => f.flipped()),
			!this.infiniteVolume,
			this.generator && this.generator + '.flipped()',
			this.vertexNames)
	}


	static EMPTY = new B2([], false, 'B2.EMPTY', new Map())
	static R3 = new B2([], true, 'B2.EMPTY', new Map())
}
namespace B2 {}

function facesWithEdge(edge, faces) {
	return faces.mapFilter((face) => {
		var matchingEdge = face.getAllEdges().find(e => e.isCoEdge(edge))
		if (matchingEdge) {
			return {face: face, reversed: !edge.a.like(matchingEdge.a), angle: NaN, normalAtEdgeA: null, edge: matchingEdge}
		}
	})
}
/**
 *
 * @param brep
 * @param brepFace
 * @param isLine
 * @param plane2
 * @returns {Array}
 * p: intersection point
 * insideDir: currently not needed
 *
 * t: param on intersection line
 * edge: face edge doing the intersection
 * edgeT: !!
 * colinear: whether edge is colinear to intersection line
 */
function planeFaceEdgeISPsWithPlane(brep, brepFace, isLine, plane2):PSSInfo[] {
	assert(brepFace.surface.plane.containsLine(isLine))
	assert(plane2.containsLine(isLine))
	var facePlane = brepFace.surface.plane
	var ps = []
	var loops = brepFace.holes.concat([brepFace.edges])
	loops.forEach(loop => {
		var colinearSegments = loop.map((edge) => edge.colinearToLine(isLine))
		var intersectionLinePerpendicular = isLine.dir1.cross(facePlane.normal)

		loop.forEach((edge, i, edges) => {
			var j = (i + 1) % edges.length, nextEdge = edges[j]
			//console.log(edge.toSource()) {p:V3(2, -2.102, 0),
			if (colinearSegments[i]) {
				// edge colinear to intersection line
				var h = (i + edges.length - 1) % edges.length, prevEdge = edges[h]
				var colinearOutside = edge.aDir.cross(facePlane.normal)
				if (prevEdge.bDir.dot(colinearOutside) < 0) {
					ps.push({p: prevEdge.b, insideDir: edge.aDir.negated(), t: isLine.pointLambda(prevEdge.b), edge: prevEdge, edgeT: prevEdge.bT,
						colinear: false})
				}
				ps.push(
					{p: edge.a, insideDir: edge.aDir, t: isLine.pointLambda(edge.a), edge: edge, edgeT: edge.aT,
						colinear: true},
					{p: edge.b, insideDir: edge.bDir.negated(), t: isLine.pointLambda(edge.b), edge: edge, edgeT: edge.bT,
						colinear: true})
				if (nextEdge.aDir.dot(colinearOutside) > 0) {
					ps.push({p: edge.b, insideDir: edge.bDir, t: isLine.pointLambda(edge.b), edge: edge, edgeT: prevEdge.bT,
						colinear: false})
				}
			} else {
				// not necessarily a straight edge, so multiple intersections are possible
				var edgeTs = edge.edgeISTsWithPlane(plane2)
				assert(edgeTs.every(t => plane2.containsPoint(edge.curve.at(t))), edgeTs)
				for (var k = 0; k < edgeTs.length; k++) {
					var edgeT = edgeTs[k]
					if (edgeT == edge.bT) {
						// endpoint lies on intersection line
						console.log('endpoint lies on intersection line',
							intersectionLinePerpendicular.dot(edge.bDir) , intersectionLinePerpendicular.dot(nextEdge.aDir),
							intersectionLinePerpendicular.dot(edge.bDir) * intersectionLinePerpendicular.dot(nextEdge.aDir), intersectionLinePerpendicular.sce,
							edge.bDir.sce, nextEdge.aDir.sce)
						if (!colinearSegments[j] && intersectionLinePerpendicular.dot(edge.bDir) * intersectionLinePerpendicular.dot(nextEdge.aDir) > 0) {
							// next segment is not colinear and ends on different side
							console.log("adding")
							ps.push({p: edge.b, insideDir: plane2.normal.negated(), t: isLine.pointLambda(edge.b), edge: edge, edgeT: edge.bT, colinear: false})
							//console.log('end on line, next other side')
						}
					} else if (edgeT != edge.aT) {
						// edge crosses intersection line, neither starts nor ends on it
						var p = edge.curve.at(edgeT)
						assert(plane2.containsPoint(p), edge.toString(), p, edgeT, plane2.distanceToPoint(p))
						assert(isLine.containsPoint(p), edge.toString(), p, edgeT, isLine.distanceToPoint(p))
						var insideDir = plane2.normal.negated()
						ps.push({p: p, insideDir: insideDir, t: isLine.pointLambda(p), edge: edge, edgeT: edgeT, colinear: false})
						console.log('middle')
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
type PSSInfo = {p: V3, insideDir: V3, t: number, edge: Edge, edgeT: number, colinear: boolean, used?: boolean}
function faceEdgeISPsWithSurface(brep:B2, brepFace:Face, isCurves:Curve[], surface2:Surface)
:PSSInfo[][] {
	let faceSurface = brepFace.surface
	let pss = NLA.arrayFromFunction(isCurves.length, i => [])

	let loops = brepFace.holes.concat([brepFace.edges])
	loops.forEach(loop => {
		let colinearEdges:boolean[] = loop.map(edge => surface2.containsCurve(edge.curve))

		loop.forEach((edge, i, edges) => {
			let j = (i + 1) % edges.length, nextEdge = edges[j]
			//console.log(edge.toSource()) {p:V3(2, -2.102, 0),
			if (colinearEdges[i]) {
				// edge colinear to intersection line
				var h = (i + edges.length - 1) % edges.length, prevEdge = edges[h]
				var colinearOutside = edge.aDir.cross(faceSurface.normalAt(edge.a))
				if (prevEdge.bDir.dot(colinearOutside) < 0) {
					ps.push({p: prevEdge.b, insideDir: edge.aDir.negated(), t: line.pointLambda(prevEdge.b), edge: prevEdge, edgeT: prevEdge.bT,
						colinear: false})
				}
				ps.push(
					{p: edge.a, insideDir: edge.aDir, t: line.pointLambda(edge.a), edge: edge, edgeT: edge.aT,
						colinear: true},
					{p: edge.b, insideDir: edge.bDir.negated(), t: line.pointLambda(edge.b), edge: edge, edgeT: edge.bT,
						colinear: true})
				if (nextEdge.aDir.dot(colinearOutside) > 0) {
					ps.push({p: edge.b, insideDir: edge.bDir, t: line.pointLambda(edge.b), edge: edge, edgeT: prevEdge.bT,
						colinear: false})
				}
			} else {
				// not necessarily a straight edge, so multiple intersections are possible
				const edgeTs = edge.edgeISTsWithSurface(surface2)
				for (var k = 0; k < edgeTs.length; k++) {
					const edgeT = edgeTs[k]
					if (edgeT == edge.bT) {
						assert(false)
						// endpoint lies on intersection line
						console.log('endpoint lies on intersection line',
							intersectionLinePerpendicular.dot(edge.bDir) , intersectionLinePerpendicular.dot(nextEdge.aDir),
							intersectionLinePerpendicular.dot(edge.bDir) * intersectionLinePerpendicular.dot(nextEdge.aDir), intersectionLinePerpendicular.sce,
							edge.bDir.sce, nextEdge.aDir.sce)
						if (!colinearEdges[j] && intersectionLinePerpendicular.dot(edge.bDir) * intersectionLinePerpendicular.dot(nextEdge.aDir) > 0) {
							// next segment is not colinear and ends on different side
							console.log("adding")
							ps.push({p: edge.b, insideDir: plane2.normal.negated(), t: line.pointLambda(edge.b), edge: edge, edgeT: edge.bT, colinear: false})
							//console.log('end on line, next other side')
						}
					} else if (edgeT != edge.aT) {
						// edge crosses an intersection curve, neither starts nor ends on it
						let p = edge.curve.at(edgeT)
						if (!isCurves.some(isCurve => isCurve.containsPoint(p), edge.toString() + p+edgeT)) {
							console.log(isCurves, isCurves[0].sce, p.sce, edge.sce, surface2.sce, brepFace.surface.sce)
							isCurves.forEach(isCurve => (isCurve.debugToMesh(mesh1, 'curve1')))
							drPs.push({p:p, text:'This is the problematic point'})
							assert(false, isCurves.map(isc => isc.distanceToPoint(p)).join(' '))
						}
						let isCurveIndex
						if (1 == isCurves.length) {
							isCurveIndex = 0
						} else {
							isCurveIndex = isCurves.findIndex(isCurve => isCurve.containsPoint(p))
							assert(isCurves.slice(isCurveIndex + 1).every(isCurve => !isCurve.containsPoint(p)))
						}
						const isCurve = isCurves[isCurveIndex]
						let curveT = isCurve.pointLambda(p)
						if (isNaN(curveT)) {
							if (isCurve instanceof EllipseCurve) {
								let hint = edge.curve.tangentAt(edgeT).cross(isCurve.f1).dot(isCurve.f2)
								curveT = isCurve.pointLambda(p, hint)
							} else {
								assert(false)
							}
						}
						assert(!isNaN(curveT))
						const insideDir = edge.tangentAt(edgeT).cross(faceSurface.normalAt(p)).negated()
						assert(!NLA.eq0(insideDir.dot(isCurve.tangentAt(curveT))))
						pss[isCurveIndex].push({p: p, insideDir: insideDir, t: curveT, edge: edge, edgeT: edgeT, colinear: false})
						console.log('middle')
					}
				}
			}
		})
	})
	// duplicate 't's are ok, as sometimes a segment needs to stop and start again
	// should be sorted so that back facing ones are first
	pss.forEach((ps, isCurveIndex) => ps.sort((a, b) => a.t - b.t || a.insideDir.dot(isCurves[isCurveIndex].tangentAt(a.t))))
	if (isCurves[0] instanceof EllipseCurve) {
		drPs.pushAll(pss[0].map(info => ({p:info.p, text:info.toSource()})))
	}
	return pss
}


var INSIDE = 0, OUTSIDE = 1, COPLANAR_SAME = 2, COPLANAR_OPPOSITE= 3, ALONG_EDGE_OR_PLANE = 4
/**
 *
 * @param {B2} brep BREP to check
 * @param {Edge} edge edge to check
 * @param {V3} dirAtEdgeA the direction vector to check
 * @param {V3} faceNormal If dirAtEdgeA doesn't split a volume, but is along a face, the returned value depends on
 *     wether that faces normal points in the same direction as faceNormal
 * @returns {number} INSIDE, OUTSIDE, COPLANAR_SAME or COPLANAR_OPPOSITE
 */
function splitsVolumeEnclosingFaces(brep, edge, dirAtEdgeA, faceNormal) {
	assert(arguments.length == 4)
	//assert(p.equals(edge.a))
	var ab1 = edge.aDir.normalized()
	var relFaces = facesWithEdge(edge, brep.faces)
	relFaces.forEach(faceInfo => {
		faceInfo.normalAtEdgeA = faceInfo.face.surface.normalAt(edge.a)
		faceInfo.edgeDirAtEdgeA = !faceInfo.reversed
				? faceInfo.edge.aDir
				: faceInfo.edge.bDir
		faceInfo.outsideVector = faceInfo.edgeDirAtEdgeA.cross(faceInfo.normalAtEdgeA)
		faceInfo.angle = (dirAtEdgeA.angleRelativeNormal(faceInfo.outsideVector.negated(), ab1) + 2 * Math.PI + NLA_PRECISION / 2) % (2 * Math.PI)
	})
	assert(relFaces.length != 0, edge.toSource())
	relFaces.sort((a, b) => a.angle - b.angle)
	// assert(relFaces.length % 2 == 0, edge.toSource()) // even number of touching faces

	if (NLA.eq0(relFaces[0].angle)) {
		var coplanarSame = relFaces[0].normalAtEdgeA.dot(faceNormal) > 0
		return coplanarSame ? COPLANAR_SAME : COPLANAR_OPPOSITE
	} else {
		return !relFaces[0].reversed ? INSIDE : OUTSIDE
	}
}
function splitsVolumeEnclosingCone(brep, point, dir) {
	var testPlane = P3.forAnchorAndPlaneVectors(point, dir, dir.getPerpendicular())
	var rays = []
	for (var k = 0; k < brep.faces.length; k++) {
		var face = brep.faces[k]
		if (face.getAllEdges().some(edge => edge.a.like(point))) {
			if (testPlane.isParallelToPlane(face.surface.plane)) {
				return ALONG_EDGE_OR_PLANE
			}
			var intersectionLine = L3.fromPlanes(testPlane, face.surface.plane)
			var ps = planeFaceEdgeISPsWithPlane(null, face, intersectionLine, testPlane)
			var i = 0
			while (i < ps.length) {
				var a = ps[i++], b = ps[i++]
				var out = a.p.like(point)
				if (out || b.p.like(point)) {
					var dir2 = out ? intersectionLine.dir1 : intersectionLine.dir1.negated()
					var angle = (dir.angleRelativeNormal(dir2, testPlane.normal) + 2 * Math.PI + NLA_PRECISION / 2) % (2 * Math.PI)
					rays.push({angle: angle, out: out})
				}
			}
		}
	}
	rays.sort((a, b) => a.angle - b.angle)
	//console.log("testPlane", testPlane.toSource(), "rays", rays.toSource())

	if (NLA.eq0(rays[0].angle)) {
		return ALONG_EDGE_OR_PLANE
	} else {
		return rays[0].out ? OUTSIDE : INSIDE
	}
}


/**
 * Solves a quadratic system of equations of the form
 *      a * x + b * y = c
 *      x^2 + y^2 = 1
 * This can be understood as the intersection of the unit circle with a line.
 *
 * @param {number} a double
 * @param {number} b double
 * @param {number} c double
 * @returns {{x1: number, y1: number, x2: number, y2: number}} with x1 >= x2 and y1 <= y2
 * a * b + (b -c) * (b + c)
 */
function intersectionUnitCircleLine(a, b, c) {
	assertNumbers(a, b, c)
	// TODO: disambiguate on a < b
	var term = sqrt(a * a + b * b - c * c)
	return {
		x1: (a * c + b * term) / (a * a + b * b),
		x2: (a * c - b * term) / (a * a + b * b),
		y1: (b * c - a * term) / (a * a + b * b),
		y2: (b * c + a * term) / (a * a + b * b)
	}
}
/**
 * Solves a quadratic system of equations of the form
 *      a * x + b * y = c
 *      x^2 - y^2 = 1
 * This can be understood as the intersection of the unit hyperbola with a line.
 *
 * @param {number} a
 * @param {number} b
 * @param {number} c
 * @returns {{x1:number, y1:number, x2:number, y2:number}} with x1 >= x2 and y1 <= y2
 * a * b + (b -c) * (b + c)
 */
function intersectionUnitHyperbolaLine(a, b, c) {
	assertNumbers(a, b, c)
	var aa = a*a, bb = b*b, cc = c*c
	// TODO: disambiguate on a < b
	//var xTerm = sqrt(4*cc*aa-4*(bb-aa)*(-cc-bb))
	var xTerm = 2 * sqrt(bb*cc + bb*bb - aa*bb)
	var yTerm = sqrt(4*cc*bb-4*(bb-aa)*(cc-aa))
	return {
		x1: (-2 * a * c + xTerm) / 2 / (bb - aa),
		x2: (-2 * a * c - xTerm) / 2 / (bb - aa),
		y1: (2 * b * c - yTerm) / 2 / (bb - aa),
		y2: (2 * b * c + yTerm) / 2 / (bb - aa)
	}
}
/**
 *
 * @param {number} a
 * @param {number} b
 * @param {number} c
 * @param {number} r
 * @returns {{x1: number, x2: number, y1: number, y2: number}}
 */
function intersectionCircleLine(a, b, c, r) {
	assertNumbers(a, b, c, r)
	var term = sqrt(r * r * (a * a + b * b) - c * c)
	return {
		x1: (a * c + b * term) / (a * a + b * b),
		x2: (a * c - b * term) / (a * a + b * b),
		y1: (b * c - a * term) / (a * a + b * b),
		y2: (b * c + a * term) / (a * a + b * b)
	}
}

/**
 *
 * @param {Curve} curve
 * @param {number} startT
 * @param {number} endT
 * @param {number} steps integer
 * @returns {number}
 */
function integrateCurve(curve, startT, endT, steps) {
	var step = (endT - startT) / steps
	var length = 0
	var p = curve.at(startT)
	for (var i = 0, t = startT + step; i < steps; i++, t += step) {
		var next = curve.at(t)
		length += p.distanceTo(next)
		p = next
	}
	return length
}

var gaussLegendreXs = [
	-0.0640568928626056260850430826247450385909,
	0.0640568928626056260850430826247450385909,
	-0.1911188674736163091586398207570696318404,
	0.1911188674736163091586398207570696318404,
	-0.3150426796961633743867932913198102407864,
	0.3150426796961633743867932913198102407864,
	-0.4337935076260451384870842319133497124524,
	0.4337935076260451384870842319133497124524,
	-0.5454214713888395356583756172183723700107,
	0.5454214713888395356583756172183723700107,
	-0.6480936519369755692524957869107476266696,
	0.6480936519369755692524957869107476266696,
	-0.7401241915785543642438281030999784255232,
	0.7401241915785543642438281030999784255232,
	-0.8200019859739029219539498726697452080761,
	0.8200019859739029219539498726697452080761,
	-0.8864155270044010342131543419821967550873,
	0.8864155270044010342131543419821967550873,
	-0.9382745520027327585236490017087214496548,
	0.9382745520027327585236490017087214496548,
	-0.9747285559713094981983919930081690617411,
	0.9747285559713094981983919930081690617411,
	-0.9951872199970213601799974097007368118745,
	0.9951872199970213601799974097007368118745
]
var gaussLegendreWeights = [
	0.1279381953467521569740561652246953718517,
	0.1279381953467521569740561652246953718517,
	0.1258374563468282961213753825111836887264,
	0.1258374563468282961213753825111836887264,
	0.1216704729278033912044631534762624256070,
	0.1216704729278033912044631534762624256070,
	0.1155056680537256013533444839067835598622,
	0.1155056680537256013533444839067835598622,
	0.1074442701159656347825773424466062227946,
	0.1074442701159656347825773424466062227946,
	0.0976186521041138882698806644642471544279,
	0.0976186521041138882698806644642471544279,
	0.0861901615319532759171852029837426671850,
	0.0861901615319532759171852029837426671850,
	0.0733464814110803057340336152531165181193,
	0.0733464814110803057340336152531165181193,
	0.0592985849154367807463677585001085845412,
	0.0592985849154367807463677585001085845412,
	0.0442774388174198061686027482113382288593,
	0.0442774388174198061686027482113382288593,
	0.0285313886289336631813078159518782864491,
	0.0285313886289336631813078159518782864491,
	0.0123412297999871995468056670700372915759,
	0.0123412297999871995468056670700372915759
]

function gaussLegendreQuadrature24(fn, startT, endT) {
	let result = 0
	for (let i = 0; i < gaussLegendreXs.length; i++) {
		// gauss-legendre goes from -1 to 1, so we need to scale
		let t = startT + (gaussLegendreXs[i] + 1) / 2 * (endT - startT)
		result += gaussLegendreWeights[i] * fn(t)
	}
	// again, [-1,1], so div by 2
	return result / 2 * (endT - startT)
}

function curveLengthByDerivative(df, startT, endT, steps) {
	var dt = (endT - startT) / steps
	var length = 0
	for (var i = 0, t = startT + dt / 2; i < steps; i++, t += dt) {
		length += df(t) * dt
	}
	return length
}



function curvePoint(implicitCurve, startPoint) {
	var eps = 1e-5
	var p = startPoint
	for (var i = 0; i < 4; i++) {
		var fp = implicitCurve(p.x, p.y)
		var dfpdx = (implicitCurve(p.x + eps, p.y) - fp) / eps,
			dfpdy = (implicitCurve(p.x, p.y + eps) - fp) / eps
		var scale = fp / (dfpdx * dfpdx + dfpdy * dfpdy)
		//console.log(p.$)
		p = p.minus(V3(scale * dfpdx, scale * dfpdy))
	}
	return p
}
function followAlgorithm (implicitCurve, startPoint, endPoint, stepLength, startDir, tangentEndPoints, boundsFunction) {
	NLA.assertNumbers(stepLength, implicitCurve(0, 0))
	NLA.assertVectors(startPoint, endPoint)
	assert (!startDir || startDir instanceof V3)
	var points = []
	tangentEndPoints = tangentEndPoints || []
	assert (NLA.eq0(implicitCurve(startPoint.x, startPoint.y)), 'NLA.isZero(implicitCurve(startPoint.x, startPoint.y))')
	stepLength = stepLength || 0.5
	var eps = 1e-5
	var p = startPoint, prevp = startDir ? p.minus(startDir) : p
	var i = 0
	do {
		var fp = implicitCurve(p.x, p.y)
		var dfpdx = (implicitCurve(p.x + eps, p.y) - fp) / eps,
			dfpdy = (implicitCurve(p.x, p.y + eps) - fp) / eps
		var tangent = V3.create(-dfpdy, dfpdx, 0)
		var reversedDir = p.minus(prevp).dot(tangent) < 0
		tangent = tangent.toLength(reversedDir ? -stepLength : stepLength)
		var tangentEndPoint = p.plus(tangent)
		points.push(p)
		tangentEndPoints.push(tangentEndPoint)
		prevp = p
		p = curvePoint(implicitCurve, tangentEndPoint)
	} while (i++ < 100 && (i < 4 || prevp.distanceTo(endPoint) > 1.1 * stepLength) && boundsFunction(p.x, p.x))
	// TODO gleichmäßige Verteilung der Punkte
	return points
}
// both curves must be in the same s-t coordinates for this to make sense
function intersectionICurveICurve(pCurve1, startParams1, endParams1, startDir, stepLength, pCurve2) {
	NLA.assertNumbers(stepLength, pCurve1(0, 0), pCurve2(0, 0))
	NLA.assertVectors(startParams1, endParams1)
	assert (!startDir || startDir instanceof V3)
	var vertices = []
	assert (NLA.eq0(pCurve1(startParams1.x, startParams1.y)))
	stepLength = stepLength || 0.5
	var eps = 1e-5
	var p = startParams1, prevp = p // startDir ? p.minus(startDir) : p
	var i = 0
	while (i++ < 1000 && (i < 4 || p.distanceTo(endParams1) > 1.1 * stepLength)) {
		var fp = pCurve1(p.x, p.y)
		var dfpdx = (pCurve1(p.x + eps, p.y) - fp) / eps,
			dfpdy = (pCurve1(p.x, p.y + eps) - fp) / eps
		var tangent = V3(-dfpdy, dfpdx, 0).toLength(stepLength)
		if (p.minus(prevp).dot(tangent) < 0) tangent = tangent.negated()
		prevp = p
		p = curvePoint(pCurve1, p.plus(tangent))
		vertices.push(p)
	}
	// TODO gleichmäßige Verteilung der Punkte
	return vertices

}
function intersectionICurveICurve(iCurve1, loopPoints1, iCurve2) {
	var p = loopPoints1[0], val = iCurve2(p.x, p.y), lastVal
	var iss = []
	for (var i = 0; i < loopPoints1.length; i++) {
		lastVal = val
		p = loopPoints1[i]
		val = iCurve2(p)
		if (val * lastVal <= 0) { // TODO < ?
			iss.push(newtonIterate2d(iCurve1, iCurve2, p.x, p.y))
		}
	}
	return iss
}

function intersectionPCurveISurface(parametricCurve, searchStart, searchEnd, searchStep, implicitSurface) {
	assertNumbers(searchStart, searchEnd, searchStep)
	var iss = []
	var val = implicitSurface(parametricCurve(searchStart)), lastVal
	for (var t = searchStart + searchStep; t <= searchEnd; t += searchStep) {
		lastVal = val
		val = implicitSurface(parametricCurve(t))
		if (val * lastVal <= 0) {
			iss.push(newtonIterate1d(t => implicitSurface(parametricCurve(t)), t))
		}
	}
	return iss
}
function intersectionICurvePSurface(f0, f1, parametricSurface) {

}
function cassini(a, c) {
	return (x,y) => (x*x+y*y) * (x*x+y*y) - 2 * c * c * (x * x - y * y) - (a * a * a * a - c * c * c * c)
}



