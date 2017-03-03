let eps = 1e-5

//function transformCached(what: Transformable, cache: Map<Transformable, Transformable>) {
//	let cached = cache.get(what)
//	if (!cached) {
//		cached =
//		cache.set(what, cache)
//	}
//	return cached
//}

abstract class Face extends Transformable {
    surface: Surface
    edges: Edge[]
    holes: Edge[][]
    id: int
    name: string
    "constructor": { new (surface: Surface, contour: Edge[], holes?: Edge[][], name?: string): Face }

    constructor(surface: Surface, contour: Edge[], holes?: Edge[][], name?: string) {
        super()
        Edge.assertLoop(contour)
        //assert(surface.edgeLoopCCW(contour), surface.toString()+contour.join("\n"))
        assert(contour.every(f => f instanceof Edge), 'contour.every(f => f instanceof Edge)' + contour.toSource())
        // contour.forEach(e => !surface.containsCurve(e.curve) &&
        // console.log("FAIL:"+surface.distanceToPoint(e.curve.anchor)))
        contour.forEach(e => assert(surface.containsCurve(e.curve), 'edge not in surface ' + e + surface))
        holes && holes.forEach(hole => Edge.assertLoop(hole))
        //holes && holes.forEach(hole => assert(!surface.edgeLoopCCW(hole)))
        assert(!holes || holes.constructor == Array, holes && holes.toString())
        this.surface = surface
        this.edges = contour // TODO refactor to contour
        this.holes = holes || []
        this.id = globalId++
        this.name = name
    }

    abstract intersectFace(face2: Face, thisBrep: B2, face2Brep: B2, faceMap: Map<Face, Edge[]>, edgeMap: Map<Edge, any[]>, likeSurfaceFaces: Set<string>)

    transform(m4: M4): this {
        const newEdges = this.edges.map(e => e.transform(m4))
        const newHoles = this.holes.map(hole => hole.map(e => e.transform(m4)))
        return new this.constructor(this.surface.transform(m4), newEdges, newHoles) as this
    }


    flipped() {
        const newEdges = this.edges.map(e => e.flipped()).reverse()
        const newHoles = this.holes.map(hole => hole.map(e => e.flipped()).reverse())
        return new this.constructor(this.surface.flipped(), newEdges, newHoles, this.name)
    }

    toString() {
        return 'new ' + this.constructor.name + '(' + this.surface + ', [' + this.edges.map(e => '\n\t' + e).join() + ']'
            + this.holes.map(hole => '\n\t\thole: ' + hole.join()) + ')'
    }

    toSource() {
        return `new ${this.constructor.name}(${this.surface.toSource()}, [${this.edges.map(e => '\n\t' + e.toSource()).join(',')}], [${
            this.holes.map(hole => '[' + hole.map(e => '\n\t' + e.toSource()).join(',') + ']').join(',')}])`
    }

    equals(obj: any): boolean {
	    function loopsEqual(a, b) {
		    return a.length == b.length &&
			    NLA.arrayRange(0, a.length, 1)
				    .some(offset => a.every((edge, i) => edge.equals(b[(offset + i) % a.length])))

	    }

	    return this == obj ||
		    Object.getPrototypeOf(this) == Object.getPrototypeOf(obj) // TODO?
		    && this.holes.length == obj.holes.length
		    && loopsEqual(this.edges, obj.edges)
		    && this.holes.every(hole => obj.holes.some(hole2 => loopsEqual(hole, hole2)))
    }

    hashCode() {
	    function arrayHashCode(array) {
		    let hashCode = 0
		    for (const val of array) {
			    hashCode = hashCode * 31 + val
		    }
		    return hashCode
	    }
	    function loopHashCode(loop) { return arrayHashCode(loop.map(edge => edge.hashCode()).sort()) }
	    return (arrayHashCode(this.holes.map(loop => loopHashCode(loop)).sort()) * 31 + loopHashCode(this.edges)) * 31 + this.surface.hashCode()
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

    getAllEdges():Edge[] {
        return Array.prototype.concat.apply(this.edges, this.holes)
    }

    addEdgeLines(mesh) {
        assert(false, "buggy, fix")
        const vertices = this.edges.map(edge => edge.getVerticesNo0()).concatenated(), mvl = mesh.vertices.length
        for (let i = 0; i < vertices.length; i++) {
            mesh.vertices.push(vertices[i])
            mesh.lines.push(mvl + i, mvl + (i + 1) % vertices.length)

        }
    }

    containsPoint(p: V3): boolean {
        assertVectors(p)
        return this.surface.loopContainsPoint(this.edges, p) != PointVsFace.OUTSIDE
            && !this.holes.some(hole => this.surface.loopContainsPoint(hole, p) != PointVsFace.OUTSIDE)
    }

    containsPoint2(p: V3): PointVsFace {
        assertVectors(p)
        const contourContainsPoint = this.surface.loopContainsPoint(this.edges, p)
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
        const containedIntersectionsTs = this.surface.isTsForLine(line).filter(t => this.containsPoint(line.at(t)))
        const nearestPointT = containedIntersectionsTs.withMax(t => -t)

        return undefined != nearestPointT ? nearestPointT : NaN
    }

    toMesh() {
        let mesh = new GL.Mesh({triangles: true, normals: true, lines: true})
        this.addToMesh(mesh)
        //mesh.compile()
	    console.log(mesh)
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
        return this.holes.concat(this.edges)
    }

    getAABB(): AABB {
        return AABB.forAAABBs(this.edges.map(e => e.getAABB()))
    }

	static create(surface: Surface, faceEdges: Edge[], holes?: Edge[][], faceName?: string) {
		return surface instanceof PlaneSurface
			? new PlaneFace(surface, faceEdges, holes, faceName)
			: new RotationFace(surface, faceEdges, holes, faceName)
	}
}

class PlaneFace extends Face {

    surface: PlaneSurface

    constructor(p: P3 | PlaneSurface, contour: Edge[], holes?: Edge[][], name?: string) {
        super(p instanceof P3 ? new PlaneSurface(p) : p, contour, holes, name)
    }


    zDirVolume(): {centroid: V3, volume: number} {
        let {centroid, area} = this.calculateArea()
        return {volume: this.surface.plane.normal.z * centroid.z * area,
            centroid: new V3(centroid.x, centroid.y, centroid.z / 2) }

    }

    calculateArea(): {centroid: V3, area: number} {
        let centroid = V3.ZERO, tcs = 0, tct = 0, totalArea = 0
        let r1 = this.surface.right, u1 = this.surface.up
        this.edges.forEach(edge => {
            let edgeCentroid, edgeArea, centroidS, centroidT
            if (edge instanceof StraightEdge) {
                let midPoint = edge.a.lerp(edge.b, 0.5)
                edgeCentroid = new V3(midPoint.x, centroid.y, centroid.z / 2)
                centroidS = midPoint.dot(r1) / 2
                centroidT = midPoint.dot(u1)
                let edgeLength = edge.a.distanceTo(edge.b)
                edgeArea = edgeLength * edge.curve.dir1.dot(r1)
            } else {
                let curve = edge.curve
                if (curve instanceof EllipseCurve) {
                    let info = curve.getAreaInDir(u1, r1, edge.aT, edge.bT)
                    edgeArea = info.area
                    let parametricCentroid = this.surface.pointToParameterFunction()(info.centroid)
                    centroidS = parametricCentroid.x
                    centroidT = parametricCentroid.y
                } else if (curve instanceof BezierCurve) {
                    edgeArea = curve.getAreaInDirSurface(u1, this.surface, edge.aT, edge.bT)
                } else {
                    assertNever()
                }
            }


            tcs += edgeArea * centroidS
            tct += edgeArea * centroidT
            totalArea += edgeArea
        })
        centroid = r1.times(tcs).plus(u1.times(tct))
	    return {area: totalArea, centroid: centroid}
    }

    addToMesh(mesh) {
        const mvl = mesh.vertices.length
        const normal = this.surface.plane.normal
        const vertices = this.edges.flatMap(edge => edge.getVerticesNo0())
        for (let i = 0; i < vertices.length; i++) { mesh.lines.push(mvl + i, mvl + (i + 1) % vertices.length) }
        const holeStarts = []
        this.holes.forEach(hole => {
            holeStarts.push(vertices.length)
            vertices.pushAll(hole.flatMap(edge => edge.getVerticesNo0()))
        })
        const triangles = triangulateVertices(normal, vertices, holeStarts).map(index => index + mvl)
        Array.prototype.push.apply(mesh.vertices, vertices)
        Array.prototype.push.apply(mesh.triangles, triangles)
        Array.prototype.push.apply(mesh.normals, NLA.arrayFromFunction(vertices.length, i => normal))
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
        return new PlaneFace(this.surface, this.edges, [holeEdges])
    }

    intersectFace(face2, thisBrep, face2Brep, faceMap, edgeMap, edgePairs) {
        if (face2 instanceof PlaneFace) {
            this.intersectPlaneFace(face2, thisBrep, face2Brep, faceMap, edgeMap, edgePairs)
            return
        }
        if (face2 instanceof RotationFace) {
            face2.intersectFace(this, face2Brep, thisBrep, faceMap, edgeMap, edgePairs)
            return
        }
        assert(false)
    }

    intersectPlaneFace(face2: PlaneFace,
                       thisBrep: B2,
                       face2Brep: B2,
                       faceMap: Map<Face, Edge[]>,
                       faceEdgePoints: Map<Face, { edge: Edge, edgeT: number, p: V3, passEdge?: Edge }[]>,
                       colinearEdgePairs: Set<NLA.Pair<any, any>>) {

        function hasPair(a: NLA.Equalable, b: NLA.Equalable) {
            return colinearEdgePairs.has(new NLA.Pair(a, b))
        }
        function addPair(a: NLA.Equalable, b: NLA.Equalable) {
            return colinearEdgePairs.add(new NLA.Pair(a, b))
        }

        /**
         *
         * @param newEdge generated segment
         * @param col1 if newEdge is colinear to an edge of this, the edge in question
         * @param col2 same for face2
         * @param in1
         * @param in2
         * @param a
         * @param b
         */
        function handleNewEdge(newEdge: StraightEdge, col1: Edge, col2: Edge, in1, in2, a, b) {
            // console.log("handling newEdge", col1 && col1.toSource(), col2 && col2.toSource(), newEdge.toSource())
            if (!col1 && !col2) {
                NLA.mapPush(faceMap, face, newEdge)
                NLA.mapPush(faceMap, face2, newEdge.flipped())
                return true
            }
            function handleEdgeInFace(col1, col2, face, face2, thisBrep, face2Brep, coplanarSameIsInside: boolean) {
                if (col1 && !col2) {
                    if (hasPair(col1.getCanon(), face2)) return

                    addPair(col1.getCanon(), face2)
                    const face2Plane = face2.surface.plane

                    //         // NB: a new edge is inserted even though it may be the same as an old one
                    //         // however it indicates that it intersects the other volume here, i.e. the old edge cannot
                    //         // be counted as "inside" for purposes of reconstitution
                    thisBrep.edgeFaces.get(col1.getCanon()).forEach(faceInfo => {
                        //const dot = NLA.snap0(face2Plane.normal.dot(faceInfo.inside))
                        //if (dot == 0 ? !coplanarSameIsInside : dot < 0) {
                        const pointsInsideFace = fff(faceInfo, face2.surface)
                        if (pointsInsideFace == INSIDE || !coplanarSameIsInside && pointsInsideFace == COPLANAR_SAME) {
                            // we need to pass the edge through so that the curve is correct in the final B2
                            const faceInfoEdgeAT = faceInfo.edge.curve.pointLambda(newEdge.a)
                            const faceInfoEdgeBT = faceInfo.edge.curve.pointLambda(newEdge.b)
                            const flipEdge = faceInfo.edge.tangentAt(faceInfoEdgeAT).dot(newEdge.aDir) < 0
                            NLA.mapPush(faceEdgePoints, faceInfo.face,
                                {edge: faceInfo.edge, edgeT: faceInfoEdgeAT, p: newEdge.a, passEdge: flipEdge ? newEdge.flipped() : undefined})
                            NLA.mapPush(faceEdgePoints, faceInfo.face,
                                {edge: faceInfo.edge, edgeT: faceInfoEdgeBT, p: newEdge.b, passEdge: !flipEdge ? newEdge : undefined})


                            for (const [face, ps] of faceEdgePoints.entries()) {
                                assert(ps.every(p => p.edge.isValidT(p.edgeT)), face.sce, ps.map(p => p.str).sce)
                            }
                        }
                    })
                    const newEdgeInside = face2Plane.normal.cross(newEdge.aDir)
                    const sVEF1 = splitsVolumeEnclosingFaces(thisBrep, col1.getCanon(), newEdgeInside, face2Plane.normal)
                    let addNewEdge, addNewEdgeFlipped
                    if (addNewEdge = sVEF1 == INSIDE || coplanarSameIsInside && sVEF1 == COPLANAR_SAME) {
                        NLA.mapPush(faceMap, face2, newEdge)
                    }
                    const sVEF2 = splitsVolumeEnclosingFaces(thisBrep, col1.getCanon(), newEdgeInside.negated(), face2Plane.normal)
                    if (addNewEdgeFlipped = sVEF2 == INSIDE || coplanarSameIsInside && sVEF2 == COPLANAR_SAME) {
                        NLA.mapPush(faceMap, face2, newEdge.flipped())
                    }
                    if (addNewEdge || addNewEdgeFlipped || sVEF1 == COPLANAR_SAME && sVEF2 == INSIDE || sVEF2 == COPLANAR_SAME && sVEF1 == INSIDE) {
                        return true
                    }
                }
            }
            const c1 = handleEdgeInFace(col1, col2, face, face2, thisBrep, face2Brep, false)
            const c2 = handleEdgeInFace(col2, col1, face2, face, face2Brep, thisBrep, true)
            if (c1 || c2) return true

            if (col1 && col2) {
                if (hasPair(col1.getCanon(), col2.getCanon())) return

                addPair(col1.getCanon(), col2.getCanon())

                function handleColinearEdgeFaces(col1, col2, thisBrep, face2Brep, coplanarSameIsInside: boolean) {
                    // not entirely sure for what i had the dirInsides in?
                    //const aDirNegatedInside = (newEdge.a.like(col2.a) || newEdge.a.like(col2.b)) && splitsVolumeEnclosingCone(face2Brep, newEdge.a, newEdge.aDir.negated()) == INSIDE
                    //const bDirInside = (newEdge.b.like(col2.a) || newEdge.b.like(col2.b)) && splitsVolumeEnclosingCone(face2Brep, newEdge.b, newEdge.bDir) == INSIDE
                    for (const faceInfo of thisBrep.edgeFaces.get(col1.getCanon())) {
                        const sVEF = splitsVolumeEnclosingFaces(face2Brep, col2.getCanon(), faceInfo.inside, faceInfo.normalAtCanonA)
                        const edgeInside = sVEF == INSIDE || coplanarSameIsInside && sVEF == COPLANAR_SAME
                        //if (edgeInside || aDirNegatedInside) {
                        //   NLA.mapPush(faceEdgePoints, face, {edge: faceInfo.edge, edgeT: col1.curve.pointLambda(newEdge.a), p: newEdge.a})
                        //}
                        //if (edgeInside || bDirInside) {
                        //    NLA.mapPush(faceEdgePoints, face, {edge: faceInfo.edge, edgeT: col1.curve.pointLambda(newEdge.b), p: newEdge.b})
                        //}
                        if (edgeInside) {
                            // we need to pass the edge through so that the curve is correct in the final B2
                            const faceInfoEdgeAT = faceInfo.edge.curve.pointLambda(newEdge.a)
                            const faceInfoEdgeBT = faceInfo.edge.curve.pointLambda(newEdge.b)
                            const flipEdge = faceInfo.edge.tangentAt(faceInfoEdgeAT).dot(newEdge.aDir) < 0
                            NLA.mapPush(faceEdgePoints, faceInfo.face,
                                {edge: faceInfo.edge, edgeT: faceInfoEdgeAT, p: newEdge.a, passEdge: flipEdge ? newEdge.flipped() : undefined})
                            NLA.mapPush(faceEdgePoints, faceInfo.face,
                                {edge: faceInfo.edge, edgeT: faceInfoEdgeBT, p: newEdge.b, passEdge: !flipEdge ? newEdge : undefined})
                        }
                    }
                }
                handleColinearEdgeFaces(col1, col2, thisBrep, face2Brep, true)
                handleColinearEdgeFaces(col2, col1, face2Brep, thisBrep, false)
            }
        }


        // what needs to be generated: new edges on face
        // points on edges where they are cut by faces so that sub edges will be generated for loops
        // points on ends of edges where the edge will be an edge in the new volume where it goes from A to B
        //         you don't want thos to be marked as "inside", otherwise invalid faces will be added
        // if a face cuts a corner, nothings needs to be done, as that alone does not limit what adjacent faces will be
        function handleEndPoint(a, b, useA, useB, newEdge) {
            console.log("logging point", useA, a.toSource())
            console.log("logging point", useB, b.toSource())
            if (useA && !useB) {
                if (!a.colinear) {
                    NLA.mapPush(faceEdgePoints, face, a)
                    assert(a.edge.isValidT(a.edgeT))
                }
                // else colinear segment ends in middle of other face, do nothing
            }
            if (useB && !useA) {
                if (!b.colinear) {
                    NLA.mapPush(faceEdgePoints, face2, b)
                    assert(b.edge.isValidT(b.edgeT))
                }
                // else colinear segment ends in middle of other face, do nothing
            }
            if (useA && useB) {
                // if a or b is colinear the correct points will already have been added to the edge by handleNewEdge
                // segment starts/ends on edge/edge intersection
                function foo(a, b, face, face2, thisPlane, face2Plane, thisBrep, face2Brep, first) {
                    if (!a.colinear && a.edgeT != a.edge.aT && a.edgeT != a.edge.bT) {
                        let addedPoint = false
                        //						assert(a.edgeT != a.edge.aT && a.edgeT != a.edge.bT, "implement point intersections") todo
                        // ends on a, on colinear segment b bT != a.edge.bT &&
                        // b can be colinear, so edgeT == aT is possible
                        if (a.p.like(b.edge.a) || a.p.like(b.edge.b)) {
                            // face2brep corner on edge
                            const sVEC1 = splitsVolumeEnclosingCone(face2Brep, a.p, a.edge.aDir)
                            const sVEC2 = splitsVolumeEnclosingCone(face2Brep, a.p, a.edge.aDir.negated())
                            // if either of these return ALONG_EDGE_OR_PLANE, then the breps share a colinear edge

                            // sVEC1 != ALONG_EDGE_OR_PLANE && sVEC2 != ALONG_EDGE_OR_PLANE && sVEC1 != sVEC2
                            // while ALONG_EDGE_OR_PLANE implies that there are some colinear edges, these wont necessarily be added
                            // and the point is still necessary for a
                            // i.e. when cutting corner from box, on a box face which gets changed, at the point where the new edge
                            // touches the existing edge

                            // sVEC1 == ALONG_EDGE_OR_PLANE || sVEC2 == ALONG_EDGE_OR_PLANE || sVEC1 != sVEC2
                            // doesn't work because ALONG_EDGE_OR_PLANE means there is a colinear part, but that shouldn't necessarily be added




                            if ((sVEC1 != INSIDE) != (sVEC2 != INSIDE)) {
                                NLA.mapPush(faceEdgePoints, face, a)
                                assert(a.edge.isValidT(a.edgeT))
                                addedPoint = true
                            }
                        } else {
                            const testVector = a.edge.aDir.rejectedFrom(b.edge.aDir)
                            assert(!testVector.isZero())
                            const sVEF1 = splitsVolumeEnclosingFaces(face2Brep, b.edge.getCanon(), testVector, thisPlane.normal)
                            const sVEF2 = splitsVolumeEnclosingFaces(face2Brep, b.edge.getCanon(), testVector.negated(), thisPlane.normal)
                            if (!(INSIDE == sVEF1 && INSIDE == sVEF2
                                || (OUTSIDE == sVEF1 || COPLANAR_OPPOSITE == sVEF1) && (OUTSIDE == sVEF2 || COPLANAR_OPPOSITE == sVEF2))) {
                                NLA.mapPush(faceEdgePoints, face, a)
                                assert(a.edge.isValidT(a.edgeT))
                                addedPoint = true
                            }
                        }

                        //if (addedPoint && !hasPair(a, b)) {
                        //	addPair(a, b)
                        //
                        //}
                    }
                }
                const canonA = a, canonB = b
                foo(a, b, face, face2, thisPlane, face2Plane, thisBrep, face2Brep, true)
                foo(b, a, face2, face, face2Plane, thisPlane, face2Brep, thisBrep, false)
                //if (!b.colinear && b.edgeT != b.edge.aT && b.edgeT != b.edge.bT) {
                //	//						assert(b.edgeT != b.edge.aT && b.edgeT != b.edge.bT, "implement point intersections")
                //	// ends on b, on colinear segment a
                //	const testVector = b.edge.aDir.rejectedFrom(a.edge.aDir)
                //	const sVEF1 = splitsVolumeEnclosingFaces(thisBrep, a.edge, testVector, face2Plane.normal)
                //	const sVEF2 = splitsVolumeEnclosingFaces(thisBrep, a.edge, testVector.negated(), face2Plane.normal)
                //	if (!(INSIDE == sVEF1 && INSIDE == sVEF2
                //		|| (OUTSIDE == sVEF1 || COPLANAR_OPPOSITE == sVEF1) && (OUTSIDE == sVEF2 || COPLANAR_OPPOSITE == sVEF2))) {
                //		NLA.mapPush(faceEdgePoints, b.edge, b)
                //	}
                //}
            }
        }


        assertInst(PlaneFace, face2)
        const face: PlaneFace = this
        // get intersection
        const thisPlane = this.surface.plane, face2Plane = face2.surface.plane
        if (thisPlane.isParallelToPlane(face2Plane)) {
            if (thisPlane.like(face2Plane)) {
                // normal same and same location in space
                // addLikeSurfaceFaces(likeSurfaceFaces, this, face2)
            }
            return
        }
        const isLine = L3.fromPlanes(thisPlane, face2Plane)
        // get intersections of newCurve with other edges of face and face2
        const ps1 = planeFaceEdgeISPsWithPlane(face, isLine, face2Plane)
        const ps2 = planeFaceEdgeISPsWithPlane(face2, isLine, thisPlane)
        console.log('ps1\n', ps1.map(m => m.toSource()).join('\n'), '\nps2\n', ps2.map(m => m.toSource()).join('\n'))
        if (ps1.length == 0 || ps2.length == 0) {
            // faces to not intersect
            return
        }
        console.log(''+thisPlane+face2Plane)

        let in1 = false, in2 = false, col1: Edge = undefined, col2: Edge = undefined
        let i = 0, j = 0, last, segments = []
        let startP, startDir, startT, startA, startB, startCol1, startCol2
        while (i < ps1.length || j < ps2.length) {
            assert(i <= ps1.length)
            assert(j <= ps2.length)
            const a = ps1[i], b = ps2[j]
            assert(a || b)
            if (j == ps2.length || i < ps1.length && NLA.lt(a.t, b.t)) {
                last = a
                in1 = !in1
                // ": col1" to remember the colinear segment, as segments are only handled once it ends (in1 false)
                col1 = in1 ? a.colinear && a.edge : col1
                i++
            } else if (i == ps1.length || NLA.gt(a.t, b.t)) {
                last = b
                in2 = !in2
                col2 = in2 ? b.colinear && b.edge : col2
                j++
            } else {
                last = a
                in1 = !in1
                in2 = !in2
                //if (in1 == in2) {
                a.used = true
                b.used = true
                col1 = in1 ? a.colinear && a.edge : col1
                col2 = in2 ? b.colinear && b.edge : col2
                //}
                i++
                j++
            }
            console.log("as", a, b, in1, in2)
            if (startP && !(in1 && in2)) {
                // segment end
                let newEdge = new StraightEdge(isLine, startP, last.p, startT, last.t, null, 'genseg' + globalId++)
                startP = undefined
                last.used = true
                if (handleNewEdge(newEdge, col1, col2, in1, in2, a, b)) {
                    handleEndPoint(a, b, col1 || !!a.used, col2 || !!b.used, newEdge)
                    handleEndPoint(startA, startB, startCol1, startCol2, newEdge)
                }
            } else if (in1 && in2) {
                // new segment just started
                startP = last.p
                startDir = last.insideDir
                startT = last.t
                last.used = true
                startA = a
                startB = b
                startCol1 = col1 || !!a.used
                startCol2 = col2 || !!b.used
            }
        }

        for (const [face, ps] of faceEdgePoints.entries()) {
            assert(ps.every(p => p.edge.isValidT(p.edgeT)), face.sce, ps.map(p => p.str).sce)
        }
    }


    static forVertices(planeSurface, vs, ...holeVss): PlaneFace {
        if (planeSurface instanceof P3) {
            planeSurface = new PlaneSurface(planeSurface)
        }
        assert(isCCW(vs, planeSurface.plane.normal), 'isCCW(vs, planeSurface.plane.normal)')
        const edges = StraightEdge.chain(vs)
        holeVss.forEach(vs => assert(doubleSignedArea(vs, planeSurface.plane.normal) >= 0, 'doubleSignedArea(vs, planeSurface.plane.normal) >= 0'))
        const holes = holeVss.map(hvs => StraightEdge.chain(hvs))
        return new PlaneFace(planeSurface, edges, holes)
    }

    pointsToInside(p: V3, dir: V3): PointVsFace {
        return this.containsPoint2(p.plus(dir.times(NLA_PRECISION * 8)))
    }
}
NLA.registerClass(PlaneFace)



class RotationFace extends Face {
    constructor(rot: Surface, contour: Edge[], holes?: Edge[][], name?) {
        super(rot, contour, holes, name)
        if (rot instanceof CylinderSurface || rot instanceof ConicSurface || rot instanceof EllipsoidSurface) {
            assertf(() => RotationFace.loopDoesNotCrossPlane(this.edges, this.surface.getSeamPlane()))
        }
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
				    // check "backwards" only if if aT != t
				    if (edge.aT != t) {
					    if (checkSide(dotCurve2(edge.curve, t, seamPlane.normal, -sign(edge.bT - edge.aT)))) return false
	    			}
	    			if (edge.bT != t) {
					    if (checkSide(dotCurve2(edge.curve, t, seamPlane.normal, sign(edge.bT - edge.aT)))) return false
				    }
			    }
		    }
	    }
	    return true
    }

	getCanonSeamU(): number {
		const pointToParameterFunction = this.surface.pointToParameterFunction()
		for (const edge of this.edges) {
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
        const localEdge = this.edges[0].transform(this.surface.inverseMatrix)
		if (P3.ZX.containsCurve(localEdge.curve)) {
		    const insideVector = localEdge.a.cross(localEdge.aDir)
            return sign(insideVector.dot(V3.Y)) * PI
        }
		assert(false, "Couldn't find canon seam u")
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
        // console.log("vs\n", vs.join("\n"), vs.length)
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
		const ellipsoid: EllipsoidSurface = this.surface as EllipsoidSurface
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
		        verticesUV.pushAll(verticesNo0.map(v => { const uv = ptpf(v, hint) return new V3(uv.x / uStep, uv.y / vStep, 0) }))
		        const nextStart = edgeLoop[ipp].a
		        console.log("BLAH", nextStart.str, ellipsoid.center.plus(ellipsoid.f3).str)
		        if (nextStart.like(ellipsoid.center.plus(ellipsoid.f3)) || nextStart.like(ellipsoid.center.minus(ellipsoid.f3))) {
			        console.log("FIXING")
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
	        }
        }
        assert(vertices.length == vertices.length)
		console.log(verticesUV.map(v => v.str).join("\n"))
		return {verticesUV: verticesUV, vertices: vertices, normals: normals, loopStarts: loopStarts}
    }

    unrollCylinderLoops(loops, uStep, vStep) {
	    const vertexLoops = loops.map(loop => loop.map(edge => edge.getVerticesNo0()).concatenated())
	    const vertices: V3[] = vertexLoops.concatenated()
	    const normals: V3[] = vertices.map(v => this.surface.normalAt(v))
	    // this.unrollLoop(loop).map(v => new V3(v.x / uStep, v.y / vStep, 0)))
	    const loopStarts = vertexLoops.reduce((arr, loop) => (arr.push(arr.last() + loop.length), arr), [0])
	    const pointToParameterFunction = this.surface.pointToParameterFunction()
	    const hint = this.getCanonSeamU()
	    const verticesUV = vertices.map(v => { const uv = pointToParameterFunction(v, hint) return new V3(uv.x / uStep, uv.y / vStep, 0) })
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
		assertf(() => uStep > 0 && vStep > 0, uStep, vStep, "Surface: " + this.surface)
		const triangles = []
		const f = (i, j) => this.surface.parametricFunction()(i * uStep, j * vStep)
		const normalF = (i, j) => this.surface.parametricNormal()(i * uStep, j * vStep)
		const loops = [this.edges].concat(this.holes)
		const {vertices, verticesUV, normals, loopStarts} = this.surface instanceof EllipsoidSurface
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
		//console.log("surface", this.surface.str)
		//console.log(verticesUV)
		//drPs.pushAll(verticesUV.map((v, i) => ({p: vertices[i], text: `${i} uv: ${v.toString(x => NLA.round10(x, -4))}`})))
		verticesUV.forEach(({x: u, y: v}) => {
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
			PlaneFace.prototype.addToMesh.call(this, mesh)
		}
		const partss: int[][][] = new Array(uRes * vRes)
		function fixUpPart(part, baseU, baseV) {
		    assert(baseU < uRes && baseV < vRes, `${baseU}, ${baseV}, ${uRes}, ${vRes}`)
			console.log("complete part", part, baseU, baseV)
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
		// "some" instead of forEach so we can return out of the entire function if this.edges crosses no borders and
		for (let vertexLoopIndex = 0; vertexLoopIndex < loops.length; vertexLoopIndex++) {
			let part: int[], firstPart, firstPartBaseU, firstPartBaseV
			let lastBaseV = -1, lastBaseU = -1
			let partCount = 0
			const vertexLoopStart = loopStarts[vertexLoopIndex]
			const vertexLoopLength = loopStarts[vertexLoopIndex + 1] - vertexLoopStart
			for (let vlvi = 0; vlvi < vertexLoopLength; vlvi++) {
				const vx0index = vertexLoopStart + vlvi, vx0 = verticesUV[vx0index]
				const vx1index = vertexLoopStart + (vlvi + 1) % vertexLoopLength, vx1 = verticesUV[vx1index]
				//console.log("dask", vx0index, vx1index)
				const vx01 = vx0.to(vx1)
				assert(vx0)
				const di = vx01.x, dj = vx01.y
				let vxIndex = vx0index, vx = vx0, currentT = 0
				let whileLimit = 40
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
					//console.log(vxIndex, vx.str, "vij", vxu, vxv, "d", di, dj, "ijNext", iNext, jNext, "nextT", iNextT, jNextT)
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
						//console.log("breaking ", vx1index)
						part.push(vx1index)
						break
					} else {
						const nextPoint = vx0.lerp(vx1, currentT)
						const nextPointIndex = addVertex(nextPoint.x, nextPoint.y)

						//console.log("pushing ", nextPointIndex)
						part.push(nextPointIndex)
						vx = nextPoint
						vxIndex = nextPointIndex
					}
				}
				assert(whileLimit, "whileLimit")
			}
			if (0 == partCount) {
				// complete loop
				assert(false, "found a hole, try increasing resolution")
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
			console.log("firstPart", firstPart)
		}
		console.log("calculated parts", partss)
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
                            "oob u1 v1 " + u1 + " " + v1 + " " + index + " " + p.str + "IF THIS FAILS check canonSeamU is correct")
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
							let flipping = eq0(currentOpos % 1)
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
						//console.log("outline", col, row, outline)
					}
				}
			}
		}
		//console.log('trinagle', triangles.max(), vertices.length, triangles.length, triangles.toSource(), triangles.map(col => vertices[col].$).toSource() )
		//assert(normals.every(n => n.hasLength(1)), normals.find(n => !n.hasLength(1)).length() +" "+normals.findIndex(n => !n.hasLength(1)))
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
	 //   console.log("u", minU, maxU, "v", minV, maxV, vertexLoops[0].toSource().replace(/\), /g, ",\n"))
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
    //    //assert(normals.every(n => n.hasLength(1)), normals.find(n => !n.hasLength(1)).length() +" "+normals.findIndex(n => !n.hasLength(1)))
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
        const vertexLoops = this.holes.concat([this.edges]).map(loop => this.unrollLoop(loop))
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
	    console.log("zzzs", minZ, maxZ, vertexLoops[0].toSource().replace(/\), /g, ",\n"))
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
	    console.log("detailsZs", detailZs)
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
        //assert(normals.every(n => n.hasLength(1)), normals.find(n => !n.hasLength(1)).length() +" "+normals.findIndex(n => !n.hasLength(1)))
        Array.prototype.push.apply(mesh.vertices, vertices)
        Array.prototype.push.apply(mesh.triangles, triangles)
        Array.prototype.push.apply(mesh.normals, normals)
        //this.addEdgeLines(mesh)

    }


    intersectFace(face2: PlaneFace,
                  thisBrep: B2,
                  face2Brep: B2,
                  faceMap: Map<Face, Edge[]>,
                  faceEdgePoints: Map<Face, { edge: Edge, edgeT: number, p: V3, passEdge?: Edge }[]>,
                  colinearEdgePairs: Set<NLA.Pair<any, any>>) {

        function hasPair(a: NLA.Equalable, b: NLA.Equalable) {
            return colinearEdgePairs.has(new NLA.Pair(a, b))
        }
        function addPair(a: NLA.Equalable, b: NLA.Equalable) {
            return colinearEdgePairs.add(new NLA.Pair(a, b))
        }


        /**
         *
         * @param newEdge generated segment
         * @param col1 if newEdge is colinear to an edge of this, the edge in question
         * @param col2 same for face2
         * @param in1
         * @param in2
         * @param a
         * @param b
         */
        function handleNewEdge(newEdge: Edge, col1: Edge, col2: Edge, in1, in2, a, b) {
            // console.log("handling newEdge", col1 && col1.toSource(), col2 && col2.toSource(), newEdge.toSource())
            if (!col1 && !col2) {
                NLA.mapPush(faceMap, face, newEdge)
                NLA.mapPush(faceMap, face2, newEdge.flipped())
                return true
            }
            function handleEdgeInFace(col1, col2, face, face2, thisBrep, face2Brep, coplanarSameIsInside: boolean) {
                if (col1 && !col2) {
                    if (hasPair(col1.getCanon(), face2)) return

                    addPair(col1.getCanon(), face2)
                    const face2Plane = face2.surface.plane

                    // NB: a new edge is inserted even though it may be the same as an old one
                    // however it indicates that it intersects the other volume here, i.e. the old edge cannot
                    // be counted as "inside" for purposes of reconstitution
                    thisBrep.edgeFaces.get(col1.getCanon()).forEach(faceInfo => {
                        const dot = NLA.snap0(face2Plane.normal.dot(faceInfo.inside))
                        const pointsInsideFace = fff(faceInfo, face2.surface)
                        if (pointsInsideFace == INSIDE || !coplanarSameIsInside && pointsInsideFace == COPLANAR_SAME) {
                            // we need to pass the edge through so that the curve is correct in the final B2
                            const faceInfoEdgeAT = faceInfo.edge.curve.pointLambda(newEdge.a)
                            const faceInfoEdgeBT = faceInfo.edge.curve.pointLambda(newEdge.b)
                            const flipEdge = faceInfo.edge.tangentAt(faceInfoEdgeAT).dot(newEdge.aDir) < 0
                            NLA.mapPush(faceEdgePoints, faceInfo.face,
                                {edge: faceInfo.edge, edgeT: faceInfoEdgeAT, p: newEdge.a, passEdge: flipEdge ? newEdge.flipped() : undefined})
                            NLA.mapPush(faceEdgePoints, faceInfo.face,
                                {edge: faceInfo.edge, edgeT: faceInfoEdgeBT, p: newEdge.b, passEdge: !flipEdge ? newEdge : undefined})


                            for (const [face, ps] of faceEdgePoints.entries()) {
                                assert(ps.every(p => p.edge.isValidT(p.edgeT)), face.sce, ps.map(p => p.str).sce)
                            }
                        }
                    })
                    const newEdgeInside = face2Plane.normal.cross(newEdge.aDir)
                    const sVEF1 = splitsVolumeEnclosingFaces(thisBrep, col1.getCanon(), newEdgeInside, face2Plane.normal)
                    let addNewEdge, addNewEdgeFlipped
                    if (addNewEdge = sVEF1 == INSIDE || coplanarSameIsInside && sVEF1 == COPLANAR_SAME) {
                        NLA.mapPush(faceMap, face2, newEdge)
                    }
                    const sVEF2 = splitsVolumeEnclosingFaces(thisBrep, col1.getCanon(), newEdgeInside.negated(), face2Plane.normal)
                    if (addNewEdgeFlipped = sVEF2 == INSIDE || coplanarSameIsInside && sVEF2 == COPLANAR_SAME) {
                        NLA.mapPush(faceMap, face2, newEdge.flipped())
                    }
                    if (addNewEdge || addNewEdgeFlipped || sVEF1 == COPLANAR_SAME && sVEF2 == INSIDE || sVEF2 == COPLANAR_SAME && sVEF1 == INSIDE) {
                        return true
                    }
                }
            }
            const c1 = handleEdgeInFace(col1, col2, face, face2, thisBrep, face2Brep, false)
            const c2 = handleEdgeInFace(col2, col1, face2, face, face2Brep, thisBrep, true)
            if (c1 || c2) return true

            if (col1 && col2) {
                if (hasPair(col1.getCanon(), col2.getCanon())) return

                addPair(col1.getCanon(), col2.getCanon())

                function handleColinearEdgeFaces(col1, col2, thisBrep, face2Brep, coplanarSameIsInside: boolean) {
                    // not entirely sure for what i had the dirInsides in?
                    //const aDirNegatedInside = (newEdge.a.like(col2.a) || newEdge.a.like(col2.b)) && splitsVolumeEnclosingCone(face2Brep, newEdge.a, newEdge.aDir.negated()) == INSIDE
                    //const bDirInside = (newEdge.b.like(col2.a) || newEdge.b.like(col2.b)) && splitsVolumeEnclosingCone(face2Brep, newEdge.b, newEdge.bDir) == INSIDE
                    for (const faceInfo of thisBrep.edgeFaces.get(col1.getCanon())) {
                        const sVEF = splitsVolumeEnclosingFaces(face2Brep, col2.getCanon(), faceInfo.inside, faceInfo.normalAtCanonA)
                        const edgeInside = sVEF == INSIDE || coplanarSameIsInside && sVEF == COPLANAR_SAME
                        //if (edgeInside || aDirNegatedInside) {
                        //   NLA.mapPush(faceEdgePoints, face, {edge: faceInfo.edge, edgeT: col1.curve.pointLambda(newEdge.a), p: newEdge.a})
                        //}
                        //if (edgeInside || bDirInside) {
                        //    NLA.mapPush(faceEdgePoints, face, {edge: faceInfo.edge, edgeT: col1.curve.pointLambda(newEdge.b), p: newEdge.b})
                        //}
                        if (edgeInside) {
                            // we need to pass the edge through so that the curve is correct in the final B2
                            const faceInfoEdgeAT = faceInfo.edge.curve.pointLambda(newEdge.a)
                            const faceInfoEdgeBT = faceInfo.edge.curve.pointLambda(newEdge.b)
                            const flipEdge = faceInfo.edge.tangentAt(faceInfoEdgeAT).dot(newEdge.aDir) < 0
                            NLA.mapPush(faceEdgePoints, faceInfo.face,
                                {edge: faceInfo.edge, edgeT: faceInfoEdgeAT, p: newEdge.a, passEdge: flipEdge ? newEdge.flipped() : undefined})
                            NLA.mapPush(faceEdgePoints, faceInfo.face,
                                {edge: faceInfo.edge, edgeT: faceInfoEdgeBT, p: newEdge.b, passEdge: !flipEdge ? newEdge : undefined})
                        }
                    }
                }
                handleColinearEdgeFaces(col1, col2, thisBrep, face2Brep, true)
                handleColinearEdgeFaces(col2, col1, face2Brep, thisBrep, false)
            }
        }


        // what needs to be generated: new edges on face
        // points on edges where they are cut by faces so that sub edges will be generated for loops
        // points on ends of edges where the edge will be an edge in the new volume where it goes from A to B
        //         you don't want thos to be marked as "inside", otherwise invalid faces will be added
        // if a face cuts a corner, nothings needs to be done, as that alone does not limit what adjacent faces will be
        function handleEndPoint(a, b, useA, useB, newEdge) {
            console.log("logging point", useA, a.toSource())
            console.log("logging point", useB, b.toSource())
            if (useA && !useB) {
                if (!a.colinear) {
                    NLA.mapPush(faceEdgePoints, face, a)
                    assert(a.edge.isValidT(a.edgeT))
                }
                // else colinear segment ends in middle of other face, do nothing
            }
            if (useB && !useA) {
                if (!b.colinear) {
                    NLA.mapPush(faceEdgePoints, face2, b)
                    assert(b.edge.isValidT(b.edgeT))
                }
                // else colinear segment ends in middle of other face, do nothing
            }
            if (useA && useB) {
                // if a or b is colinear the correct points will already have been added to the edge by handleNewEdge
                // segment starts/ends on edge/edge intersection
                function foo(a, b, face, face2, thisPlane, face2Plane, thisBrep, face2Brep, first) {
                    if (!a.colinear && a.edgeT != a.edge.aT && a.edgeT != a.edge.bT) {
                        let addedPoint = false
                        //						assert(a.edgeT != a.edge.aT && a.edgeT != a.edge.bT, "implement point intersections") todo
                        // ends on a, on colinear segment b bT != a.edge.bT &&
                        // b can be colinear, so edgeT == aT is possible
                        if (a.p.like(b.edge.a) || a.p.like(b.edge.b)) {
                            // face2brep corner on edge
                            const sVEC1 = splitsVolumeEnclosingCone(face2Brep, a.p, a.edge.aDir)
                            const sVEC2 = splitsVolumeEnclosingCone(face2Brep, a.p, a.edge.aDir.negated())
                            // if either of these return ALONG_EDGE_OR_PLANE, then the breps share a colinear edge

                            // sVEC1 != ALONG_EDGE_OR_PLANE && sVEC2 != ALONG_EDGE_OR_PLANE && sVEC1 != sVEC2
                            // while ALONG_EDGE_OR_PLANE implies that there are some colinear edges, these wont necessarily be added
                            // and the point is still necessary for a
                            // i.e. when cutting corner from box, on a box face which gets changed, at the point where the new edge
                            // touches the existing edge

                            // sVEC1 == ALONG_EDGE_OR_PLANE || sVEC2 == ALONG_EDGE_OR_PLANE || sVEC1 != sVEC2
                            // doesn't work because ALONG_EDGE_OR_PLANE means there is a colinear part, but that shouldn't necessarily be added




                            if ((sVEC1 != INSIDE) != (sVEC2 != INSIDE)) {
                                NLA.mapPush(faceEdgePoints, face, a)
                                assert(a.edge.isValidT(a.edgeT))
                                addedPoint = true
                            }
                        } else {
                            const testVector = a.edge.aDir.rejectedFrom(b.edge.aDir)
                            assert(!testVector.isZero())
                            const sVEF1 = splitsVolumeEnclosingFaces(face2Brep, b.edge.getCanon(), testVector, thisPlane.normal)
                            const sVEF2 = splitsVolumeEnclosingFaces(face2Brep, b.edge.getCanon(), testVector.negated(), thisPlane.normal)
                            if (!(INSIDE == sVEF1 && INSIDE == sVEF2
                                || (OUTSIDE == sVEF1 || COPLANAR_OPPOSITE == sVEF1) && (OUTSIDE == sVEF2 || COPLANAR_OPPOSITE == sVEF2))) {
                                NLA.mapPush(faceEdgePoints, face, a)
                                assert(a.edge.isValidT(a.edgeT))
                                addedPoint = true
                            }
                        }

                        //if (addedPoint && !hasPair(a, b)) {
                        //	addPair(a, b)
                        //
                        //}
                    }
                }
                const canonA = a, canonB = b
                foo(a, b, face, face2, surface, surface2, thisBrep, face2Brep, true)
                foo(b, a, face2, face, surface2, surface, face2Brep, thisBrep, false)
                //if (!b.colinear && b.edgeT != b.edge.aT && b.edgeT != b.edge.bT) {
                //	//						assert(b.edgeT != b.edge.aT && b.edgeT != b.edge.bT, "implement point intersections")
                //	// ends on b, on colinear segment a
                //	const testVector = b.edge.aDir.rejectedFrom(a.edge.aDir)
                //	const sVEF1 = splitsVolumeEnclosingFaces(thisBrep, a.edge, testVector, face2Plane.normal)
                //	const sVEF2 = splitsVolumeEnclosingFaces(thisBrep, a.edge, testVector.negated(), face2Plane.normal)
                //	if (!(INSIDE == sVEF1 && INSIDE == sVEF2
                //		|| (OUTSIDE == sVEF1 || COPLANAR_OPPOSITE == sVEF1) && (OUTSIDE == sVEF2 || COPLANAR_OPPOSITE == sVEF2))) {
                //		NLA.mapPush(faceEdgePoints, b.edge, b)
                //	}
                //}
            }
        }


        assertInst(Face, face2)

        const face = this
        const surface = face.surface, surface2 = face2.surface
        if (surface.isCoplanarTo(surface2)) {
            return
        }
        const isCurves = surface.isCurvesWithSurface(surface2)
        // get intersections of newCurve with other edges of face and face2
        const pss1 = faceEdgeISPsWithSurface(face, isCurves, face2.surface)
        const pss2 = faceEdgeISPsWithSurface(face2, isCurves, face.surface)
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
            function startsInside(ps, otherFace) {
                if (0 == ps.length) {
                    return undefined !== isCurve.minT && otherFace.containsPoint(isCurve.at(isCurve.minT))
                } else {
                    return ps[0].insideDir.dot(isCurve.tangentAt(ps[0].t)) < 0
                }
            }
            // they can't both be empty currently
            // they can't both start "inside"
            let in1 = startsInside(ps1, face2)
            let in2 = startsInside(ps2, face)
            let col1: Edge, col2: Edge
            let i = 0, j = 0, last, segments = []
            let startP, startDir, startT = isCurve.minT, startA, startB, startCol1, startCol2
            while (i < ps1.length || j < ps2.length) {
                assert(i <= ps1.length)
                assert(j <= ps2.length)
                const a = ps1[i], b = ps2[j]
                assert(a || b)
                if (j == ps2.length || i < ps1.length && NLA.lt(a.t, b.t)) {
                    last = a
                    in1 = !in1
                    // ": col1" to remember the colinear segment, as segments are only handled once it ends (in1 false)
                    col1 = in1 ? a.colinear && a.edge : col1
                    i++
                } else if (i == ps1.length || NLA.gt(a.t, b.t)) {
                    last = b
                    in2 = !in2
                    col2 = in2 ? b.colinear && b.edge : col2
                    j++
                } else {
                    last = a
                    in1 = !in1
                    in2 = !in2
                    //if (in1 == in2) {
                    a.used = true
                    b.used = true
                    col1 = in1 ? a.colinear && a.edge : col1
                    col2 = in2 ? b.colinear && b.edge : col2
                    //}
                    i++
                    j++
                }
                console.log("as", a, b, in1, in2)
                if (startP && !(in1 && in2)) {
                    // segment end
                    startDir = isCurve.tangentAt(startT)
                    startT > last.t && (startDir = startDir.negated())
                    let endDir = isCurve.tangentAt(last.t)
                    startT > last.t && (endDir = endDir.negated())
                    const newEdge = Edge.create(isCurve, startP, last.p, startT, last.t, null, startDir, endDir, 'genseg' + globalId++)
                    startP = undefined
                    last.used = true
                    if (handleNewEdge(newEdge, col1, col2, in1, in2, a, b)) {
                        handleEndPoint(a || {colinear: false}, b, col1 || a && !!a.used, col2 || !!b.used, newEdge)
                        handleEndPoint(startA || {colinear: false}, startB, startCol1, startCol2, newEdge)
                    }
                } else if (in1 && in2) {
                    // new segment just started
                    startP = last.p
                    startDir = last.insideDir
                    startT = last.t
                    last.used = true
                    startA = a
                    startB = b
                    startCol1 = col1 || (a && !!a.used)
                    startCol2 = col2 || !!b.used
                }
            }
        })
    }


}

NLA.registerClass(RotationFace)

function addLikeSurfaceFaces(likeSurfaceFaces: Face[][], face1: Face, face2: Face) {
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

function assembleFaceFromLooseEdges(edges: Edge[], surface: Surface, faceConstructor: typeof Face.constructor): Face {
    const visited = new Set()
    function nextStart() { return edges.find(edge => !visited.has(edge)) }
    const loops = []
    let startEdge, currentEdge
    while (startEdge = nextStart()) {
        currentEdge = startEdge
        const loop = []
        let total = 0
        do {
            visited.add(currentEdge)
            loop.push(currentEdge)
            const possibleEdges = edges.filter(edge => currentEdge.b.like(edge.a))
            const normalAtCurrentB = surface.normalAt(currentEdge.b)
            const nextEdgeIndex = possibleEdges.indexWithMax(
                (edge, index) => currentEdge.bDir.angleRelativeNormal(edge.aDir, normalAtCurrentB))
            currentEdge = possibleEdges[nextEdgeIndex]
        } while (startEdge != currentEdge && total++ < 200)
        console.log(edges)
        assert(total != 201)
        loops.push(loop)
    }


    console.log(loops)
    const assembledFaces = B2.assembleFacesFromLoops(loops, surface, faceConstructor)
    assertf(() => 1 == assembledFaces.length)
    return assembledFaces[0]
}

class B2 extends Transformable {
    faces: Face[]
    infiniteVolume: boolean
    generator: string
    vertexNames: Map<V3, string>
    edgeFaces: NLA.CustomMap<Edge, {face: Face, edge: Edge, normalAtCanonA: V3, inside: V3, reversed: boolean, angle: number}[]>

    constructor(faces: Face[], infiniteVolume: boolean, generator?: string, vertexNames?) {
        super()
        this.faces = faces
        assertInst.apply(undefined, [Face].concat(faces))
        this.infiniteVolume = !!infiniteVolume
        this.generator = generator
        this.vertexNames = vertexNames
        this.edgeFaces = undefined
        this.assertSanity()
    }

    containsPoint(p: V3): boolean {
        const testLine = new L3(p, V3.randomUnit())
        for (const face of this.faces) {

        }
    }

    withMergedFaces():B2 {
        const likeSurfaceFaces = []
        for (let i = 0; i < this.faces.length; i++) {
            let addedToGroup = false
            for (let j = 0; j < i; j++) {
                if (this.faces[i].surface.isCoplanarTo(this.faces[j].surface)) {
                    const faceGroup = likeSurfaceFaces.find(faceGroup => faceGroup.contains(this.faces[j]))
                    if (faceGroup) {
                        faceGroup.push(this.faces[i])
                        addedToGroup = true
                    }
                }
            }
            !addedToGroup && likeSurfaceFaces.push([this.faces[i]])
        }

        console.log("likeSurfaceFaces", likeSurfaceFaces)
        if (likeSurfaceFaces.every(group => group.length == 1)) return this

        const newFaces = []
        let total = 0
        for (const faceGroup of likeSurfaceFaces) {
            console.log(faceGroup)
            if (faceGroup.length == 1) {
                newFaces.push(faceGroup[0])
            } else {
                const allEdges = faceGroup.map(face => face.getAllEdges()).concatenated()
                for (let i = allEdges.length; i-- > 0;) {
                    for (let j = 0; j < i; j++) {
                        console.log("blugh", total)
                        assert(i >= 0 && j >= 0 && total++ < 500, i + ' '+ j+' '+ total)
                        if (allEdges[i].isCoEdge(allEdges[j])) {
                            // remove both
                            allEdges.splice(i, 1)
                            allEdges.splice(j, 1)
                            i--
                            break
                        }
                    }
                }
                const newFace = assembleFaceFromLooseEdges(allEdges, faceGroup[0].surface, faceGroup[0].constructor)
                newFaces.push(newFace)
            }
        }

        return new B2(newFaces, this.infiniteVolume, this.generator && this.generator + '.withMergedFaces()', this.vertexNames)
    }

    calculateVolume(): number {
        return this.faces.map(face => face.zDirVolume().volume).sum()
    }

    toMesh(): GL.Mesh & {faceIndexes: Map<Face, {start: int, count: int}>} {
        const mesh = new GL.Mesh({triangles: true, normals: true, lines: true}) as any
        mesh.faceIndexes = new Map()
        this.faces.forEach((face, i) => {
            let triangleStart = mesh.triangles.length
            face.addToMesh(mesh)
            mesh.faceIndexes.set(face, {start: triangleStart, count: mesh.triangles.length - triangleStart})
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

    and(other: B2): B2 {
        return this.intersection(other, true, true,
            this.generator && other.generator && this.generator + '.and(' + other.generator + ')')
    }

    xor(other: B2): B2 {
        return new B2(
            this.minus(other).faces.concat(other.minus(this).faces),
            this.infiniteVolume != other.infiniteVolume,
            this.generator && other.generator && this.generator + '.xor(' + other.generator + ')')
    }

    equals(brep): boolean {
        return this.faces.length == brep.faces.length &&
            this.faces.every((face) => brep.faces.some((face2) => face.equals(face2)))
    }

    like(brep): boolean {
        return this.faces.length == brep.faces.length &&
            this.faces.every((face) => brep.faces.some((face2) => face.likeFace(face2)))
    }

    toString(): string {
        return `new B2([\n${this.faces.join(',\n').replace(/^/gm, '\t')}], ${this.infiniteVolume})`
    }

    toSource(): string {
        return this.generator ||
            `new B2([\n${this.faces.map(SCE).join(',\n').replace(/^/gm, '\t')}], ${this.infiniteVolume})`
    }

    static loop1ContainsLoop2(loop1: Edge[], loop2: Edge[], surface: Surface): boolean {
        for (const edge of loop2) {
            const loop1ContainsPoint = surface.loopContainsPoint(loop1, edge.a)
            if (PointVsFace.ON_EDGE != loop1ContainsPoint) return PointVsFace.INSIDE == loop1ContainsPoint
        }
        for (const edge of loop2) {
            const edgePoint = edge.curve.at(edge.aT * 0.2 + edge.bT * 0.8)
            const loop1ContainsPoint = surface.loopContainsPoint(loop1, edgePoint)
            if (PointVsFace.ON_EDGE != loop1ContainsPoint) return PointVsFace.INSIDE == loop1ContainsPoint
        }
        throw new Error(loop1.sce+ loop2.sce)
    }

    static assembleFacesFromLoops(loops: Edge[][], surface: Surface, faceConstructor): Face[] {
        type LoopInfo = {loop: Edge[], ccw: boolean, subloops: LoopInfo[]}
        function placeRecursively(newLoopInfo: LoopInfo, loopInfos: LoopInfo[]) {
            if (loopInfos.length == 0) {
                loopInfos.push(newLoopInfo)
            } else {
                const subLoopInfo = loopInfos.find(loopInfo => B2.loop1ContainsLoop2(loopInfo.loop, newLoopInfo.loop, surface))
                console.log("here", subLoopInfo)
                if (subLoopInfo) {
                    placeRecursively(newLoopInfo, subLoopInfo.subloops)
                } else {
                    // newLoopInfo isnt contained by any other subLoopInfo
                    for (let i = loopInfos.length; --i >= 0;) {
                        const subLoopInfo = loopInfos[i]
                        //console.log("cheving subLoopInfo", surface.loopContainsPoint(newLoopInfo.edges, subLoopInfo.edges[0].a))
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
            assert(loopInfo.ccw)
            if (loopInfo.subloops.every(sl => !sl.ccw)) {
                const newFace = new faceConstructor(surface, loopInfo.loop, loopInfo.subloops.map(sl => sl.loop))
                newFaces.push(newFace)
                loopInfo.subloops.forEach(sl => sl.subloops.forEach(slsl => newFacesRecursive(slsl)))
            } else {
                loopInfo.subloops.forEach(sl => sl.ccw && newFacesRecursive(sl))
            }
        }

        const newFaces: Face[] = []
        console.log(loops.map(loop=> loop.join('\n')).join('\n\n'))
        const topLevelLoops:LoopInfo[] = []
        loops.forEach(loop => placeRecursively({loop: loop, ccw: surface.edgeLoopCCW(loop), subloops: []}, topLevelLoops))
        topLevelLoops.forEach(tll => newFacesRecursive(tll))
        return newFaces
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
    static reconstituteFaces(oldFaces: Face[], edgeSubEdges, faceMap: Map<Face, Edge[]>, newFaces: Face[]): void {
        const oldFaceStatuses: Map<Face, string> = new Map()
        // reconstitute faces
        const insideEdges = []
        oldFaces.forEach((face, faceIndex) => {
            console.log('reconstituting face', face.toString())
            const usableOldEdges = face.getAllEdges().filter(edge => !edgeSubEdges.get(edge))
            const subEdges = face.getAllEdges().mapFilter(edge => edgeSubEdges.get(edge)).concatenated()
            const els = face.edges.map(edge => edgeSubEdges.get(edge) && edgeSubEdges.get(edge).join('\n')).map(s => '\n\n' + s).join()
            console.log('edgeSubEdges', els)
            const newEdges = faceMap.get(face) || []
            if (newEdges.length || subEdges.length) {
                oldFaceStatuses.set(face, 'partial')
                console.log('newEdges\n', newEdges.map(e=>e.toString()).join('\n'))
                const loops = []
                let edgeCond = face instanceof PlaneFace
                    ? (edge => edge.a.like(currentEdge.b))
                    : (edge => (edge.curve == currentEdge.curve // TODO: ??
                    ? NLA.eq(edge.aT, currentEdge.bT)
                    : edge.a.like(currentEdge.b)))
                // new edges are definitely part of a resulting loop
                // old edges (both contour and holes) can either be part of a new loop, in which case they will already
                // have been visited when starting a loop search with a new edge, OR they can be stranded, OR they can
                // remain in their old loop
                function getNextStart() {
                    return newEdges.find(edge => !edge.visited)
                        || subEdges.find(edge => !edge.visited)
                        || usableOldEdges.find(edge => !edge.visited)
                }
                usableOldEdges.forEach(edge => edge.visited = false)

                // search for a loop:
                let currentEdge
                while (currentEdge = getNextStart()) {
                    let cancelLoop = false
                    let startEdge = currentEdge, edges = [], i = 0
                    // wether only new edges are used (can include looseSegments)
                    do {
                        currentEdge.visited = true
                        console.log('currentEdge', currentEdge.b.sce, currentEdge.toSource())
                        edges.push(currentEdge)
                        // find next edge
                        const possibleOldEdges = usableOldEdges.filter(edge => currentEdge.b.like(edge.a))
                        const possibleSubEdges = subEdges.filter(edge => currentEdge.b.like(edge.a))
                        const possibleNewEdges = newEdges.filter(edge => currentEdge.b.like(edge.a))
                        const possibleEdges = possibleOldEdges.concat(possibleSubEdges, possibleNewEdges)
                        assertf(() => 0 < possibleEdges.length, face.sce)
                        const faceNormalAtCurrentB = face.surface.normalAt(currentEdge.b)
                        const nextEdgeIndex = possibleEdges.indexWithMax(
                            (edge, index) => (currentEdge.bDir.angleRelativeNormal(edge.aDir, faceNormalAtCurrentB) + NLA_PRECISION + PI) % TAU)
                        currentEdge = possibleEdges[nextEdgeIndex]
                        if (currentEdge.visited) {
                            console.log("breaking")
                            break
                        }
                        assert(currentEdge)
                        assert(currentEdge != startEdge)
                    } while (++i < 200)
                    if (200 == i) {
                        assert(false, "too many")
                    }
                    // check if we found a loop
                    if (currentEdge == startEdge) {
                        console.log('finished loop')
                        loops.push(edges)
                    }
                }

                const faceNewFaces = B2.assembleFacesFromLoops(loops, face.surface, face.constructor)
                newFaces.pushAll(faceNewFaces)
                const faceNewFacesEdges = faceNewFaces.map(face => face.getAllEdges()).concatenated()
                insideEdges.pushAll(usableOldEdges.filter(edge => faceNewFacesEdges.includes(edge)))
            }
        })
        console.log("INSIDE EDGES", insideEdges, oldFaceStatuses)
        while (insideEdges.length != 0) {
            const insideEdge = insideEdges.pop()
            const adjacentFaces = facesWithEdge(insideEdge, oldFaces)
            adjacentFaces.forEach(info => {
                if (!oldFaceStatuses.has(info.face)) {
                    oldFaceStatuses.set(info.face, 'inside')
                    insideEdges.push.apply(insideEdges, info.face.edges)
                }
            })
        }
        newFaces.pushAll(oldFaces.filter(face => oldFaceStatuses.get(face) == 'inside'))
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

    //getLooseEdgeSegmentsInsideDirs(edgeMap: Map<Edge, IntersectionPointInfo[]>): Map<Edge, Edge> {
    //	var edgeLooseSegments = new Map()
    //	console.log("edgeMap", edgeMap)
    //	edgeMap.forEach((pointInfos, baseEdge) => {
    //		// TODO: make sure it works with loop
    //		// TODO: dont generate unnecessarry segments
    //		var looseSegments = []
    //		if (!baseEdge.reversed) {
    //			pointInfos.sort((a, b) => NLA.snap(a.edgeT - b.edgeT, 0) || (b.insideDir.dot(baseEdge.curve.dir1) - a.insideDir.dot(baseEdge.curve.dir1)))
    //		} else {
    //			pointInfos.sort((b, a) => NLA.snap(a.edgeT - b.edgeT, 0) || (b.insideDir.dot(baseEdge.curve.dir1) - a.insideDir.dot(baseEdge.curve.dir1)))
    //		}
    //		var startP = baseEdge.a, startDir = baseEdge.aDir, startT = baseEdge.aT, startInfo
    //		for (var i = 0; i < pointInfos.length; i++) {
    //			var info = pointInfos[i]
    //			assert(info.insideDir, info.toSource())
    //			console.log("info.insideDir.dot(baseEdge.curve.dir1)", info.insideDir.dot(baseEdge.curve.dir1))
    //			// ignore start, end and repeating points
    //			if (NLA.eq(info.edgeT, baseEdge.bT) || NLA.eq(info.edgeT, baseEdge.aT) || NLA.eq(info.edgeT, startT)) {
    //				continue
    //			}
    //			var pDir = baseEdge.tangentAt(info.edgeT)
    //			// add segment only if insideDir points backwards
    //			if (info.insideDir.dot(baseEdge.curve.dir1) < 0) {
    //				looseSegments.push(Edge.create(baseEdge.curve, startP, info.p, startT, info.edgeT, null, startDir, pDir))
    //			}
    //			startT = info.edgeT
    //			startInfo = info
    //			startDir = pDir
    //		}
    //		if (startInfo && startInfo.insideDir.dot(baseEdge.curve.dir1) > 0) {
    //			looseSegments.push(Edge.create(baseEdge.curve, startP, baseEdge.b, startT, baseEdge.bT, null, startDir, baseEdge.bDir))
    //		}
    //		edgeLooseSegments.set(baseEdge, looseSegments)
    //	})
    //	return edgeLooseSegments
    //}

    getLooseEdgeSegments(edgePointInfoss: Map<Face, { edge: Edge, edgeT: number, p: V3, passEdge?: Edge }[]>): Map<Face, Map<Edge, Edge[]>> {
        const result = new Map()
        console.log("edgePointInfoss", edgePointInfoss)
        edgePointInfoss.forEach((pointInfos, face) => {
            // TODO: dont generate unnecessarry segments
            // if there are no point info, the original edge will be kept, so we should return nothing
            // otherwise, something will be returned, even if it a new edge identical to the base edge
	        const edgeMap = new Map<Edge, { edge: Edge, edgeT: number, p: V3, passEdge?: Edge }[]>()
	        for (const info of pointInfos) {
	        	NLA.mapPush(edgeMap, info.edge, info)
	        }
	        for (const [baseEdge, pointInfos] of edgeMap) {
		        if (0 == pointInfos.length) return []
		        const subEdges = []
		        if (!baseEdge.reversed) {
			        pointInfos.sort((a, b) => a.edgeT - b.edgeT)
		        } else {
			        pointInfos.sort((b, a) => a.edgeT - b.edgeT)
		        }
		        let startP = baseEdge.a, startDir = baseEdge.aDir, startT = baseEdge.aT, startInfo
		        for (let i = 0; i < pointInfos.length; i++) {
                    const info = pointInfos[i]
		            if (info.passEdge) {
                        //assert(startP.equals(info.passEdge.a), startP.str, info.passEdge.sce, info.sce)
		                subEdges.push(info.passEdge)
                        startP = info.p
                        startT = info.edgeT
                        startInfo = info
                        startDir = baseEdge.tangentAt(info.edgeT)
                        continue
                    }
			        // ignore start, end and repeating points
			        if (NLA.eq(info.edgeT, baseEdge.bT) || NLA.eq(info.edgeT, baseEdge.aT) || NLA.eq(info.edgeT, startT)) {
				        continue
			        }
			        const pDir = baseEdge.tangentAt(info.edgeT)
			        subEdges.push(Edge.create(baseEdge.curve, startP, info.p, startT, info.edgeT, null, startDir, pDir, 'looseSegment' + globalId++))
			        startP = info.p
			        startT = info.edgeT
			        startInfo = info
			        startDir = pDir
		        }
		        if (!NLA.eq(baseEdge.bT, startT)) {
                    subEdges.push(Edge.create(baseEdge.curve, startP, baseEdge.b, startT, baseEdge.bT, null, startDir, baseEdge.bDir, 'looseSegment' + globalId++))
                }
		        result.set(baseEdge, subEdges)
	        }
        })
        return result
    }

    getIntersectionEdges(brep2) {
        const faceMap = new Map(), edgeMap = new Map()

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

	assertSanity(): void {
		const allFaceEdges = this.faces.map(face => face.getAllEdges()).concatenated()
		for (const {i, j} of NLA.combinations(allFaceEdges.length)) {
			const a = allFaceEdges[i], b = allFaceEdges[j]
			//assert(i == j || !a.isCoEdge(b) || a == b || a.flippedOf == b, 'coedges not linked properly', a, b)

			assert(i == j
				|| !a.curve.isColinearTo(b.curve)
				|| (a.curve.equals(b.curve) && a.isCoEdge(b))
                || !a.overlaps(b), 'colinear edges overlap', a, b)
		}
	}


    buildAdjacencies(): void {
        if (this.edgeFaces) return

        this.edgeFaces = new NLA.CustomMap()
        for (const face of this.faces) {
            for (const edge of face.getAllEdges()) {
                const canon = edge.getCanon()
                const normalAtCanonA = face.surface.normalAt(canon.a)
                const inside = normalAtCanonA.cross(canon == edge ? edge.aDir : edge.bDir)
                NLA.mapPush(this.edgeFaces, canon,
                    {face: face, edge: edge, normalAtCanonA: normalAtCanonA, reversed: canon != edge, inside: inside, angle: 0})
            }
        }

        for (const [canonEdge, edgeFaceInfos] of this.edgeFaces) {
	        // TODO handle curved faces
            assert(edgeFaceInfos.length % 2 == 0, canonEdge, edgeFaceInfos.sce)
            const faceInfo0 = edgeFaceInfos.find(faceInfo => faceInfo.reversed)
	        if (!faceInfo0) {
		        console.warn("invalid brep")
		        continue
	        }
            edgeFaceInfos.forEach(faceInfo => {
	            if (faceInfo != faceInfo0) {
		            faceInfo.angle = faceInfo0.inside.angleRelativeNormal(faceInfo.inside, canonEdge.aDir.normalized())
		            if (faceInfo.angle < 0) faceInfo.angle += 2 * Math.PI
	            }
            })
            edgeFaceInfos.sort((a, b) => NLA.snap(a.angle - b.angle, 0) || assertNever())
        }
    }

    /**
     * Cases for volumes A and B
     *
     *          1.  Volumes do not touch.
     *          2.  face/face Face surfaces intersect each other.
     *              implies edges going through faces.
     *              e.g. box(5, 5, 5) - box(5, 5, 5).translate(1, 1, 1)
     *          3.  face/edge Edge of A lies in a face of B
     *              implies vertices of A lying in face of B
     *              e.g. box(5, 5, 5) - box(3, 3, 3).rotateZ([0, 1, 2] * PI / 2).translate(0, 1, 1)
     *          4.  edge/edge Two edges are colinear.
     *              implies vertex of A lying in edge of B
     *           5.  vertex/edge Vertex of A lies on edge of B (but no edge/edge)
     *          6.  vertex/vertex with/without edge/edge, edge/face and face/face intersections
     *          7.  vertex lies in face
     *
     *
     *
     */
    intersection(other: B2, buildThis: boolean, buildOther: boolean, name: string): B2 {
        this.assertSanity()
        other.assertSanity()
        this.buildAdjacencies()
        other.buildAdjacencies()

        const faceMap = new Map(), edgeMap = new Map()

        const likeSurfaceFaces = new NLA.CustomSet()

        for (const thisFace of this.faces) {
            for (const otherFace of other.faces) {
                thisFace.intersectFace(otherFace, this, other, faceMap, edgeMap, likeSurfaceFaces)
            }
        }
        const newFaces = []
        const edgeLooseSegments = this.getLooseEdgeSegments(edgeMap)

        if (0 == faceMap.size && 0 == edgeMap.size) {
            const thisInOther = other.containsPoint(this.faces[0].edges[0].a)
            const otherInThis = !thisInOther && contians(other, this)
        } else {
            buildThis && B2.reconstituteFaces(this.faces, edgeLooseSegments, faceMap, newFaces)
            buildOther && B2.reconstituteFaces(other.faces, edgeLooseSegments, faceMap, newFaces)
        }
        //buildCoplanar && this.reconstituteCoplanarFaces(likeSurfaceFaces, edgeLooseSegments, faceMap, newFaces, this.infiniteVolume, other.infiniteVolume)

        return new B2(newFaces, this.infiniteVolume && other.infiniteVolume)

    }

    //intersection(other: B2, buildThis: boolean, buildOther: boolean, buildCoplanar): B2 {
    //    const faceMap = new Map(), edgeMap = new Map(), pairs = new NLA.CustomSet()
    //
    //    this.faces.forEach(face => {
    //        //console.log('face', face.toString())
    //        other.faces.forEach(face2 => {
    //            //console.log('face2', face2.toString())
    //            face.intersectFace(face2, this, other, faceMap, edgeMap, pairs)
    //        })
    //    })
    //    const newFaces = [];
    //
    //    /*
    //     TODO:
    //     faceMap.forEach((faceLooses, face) => {
    //     faceLooses.forEach(edge => {
    //     face.edges.forEach(faceEdge => {
    //     var edgeT = faceEdge.getEdgeT(edge.a)
    //     if (undefined !== edgeT) {
    //     console.log("WAARGH", edge.a.$, faceEdge.toString(), edgeT)
    //     NLA.mapPush(edgeMap, faceEdge, {edgeT: edgeT, p: edge.a})
    //     }
    //     })
    //     })
    //     })
    //     */
    //    const edgeLooseSegments = this.getLooseEdgeSegments(edgeMap);
    //
    //    buildThis && B2.reconstituteFaces(this.faces, edgeLooseSegments, faceMap, newFaces, other.infiniteVolume)
    //    buildOther && B2.reconstituteFaces(other.faces, edgeLooseSegments, faceMap, newFaces, this.infiniteVolume)
    //    //buildCoplanar && this.reconstituteCoplanarFaces(likeSurfaceFaces, edgeLooseSegments, faceMap, newFaces, this.infiniteVolume, other.infiniteVolume)
    //
    //    return new B2(newFaces, this.infiniteVolume && other.infiniteVolume)
    //}

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
            this.generator && desc && this.generator + desc, // if desc isn't set, the generator will be invalid
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
    static R3 = new B2([], true, 'B2.R3', new Map())
}
namespace B2 {
	export const asldk = 0
}

function facesWithEdge(edge: Edge, faces: Face[]): { face: Face, reversed: boolean, angle: number, normalAtEdgeA: V3, edge: Edge }[] {
    return faces.mapFilter((face) => {
        const matchingEdge = face.getAllEdges().find(e => e.isCoEdge(edge))
        if (matchingEdge) {
            return {face: face, reversed: !edge.a.like(matchingEdge.a), angle: NaN, normalAtEdgeA: null, edge: matchingEdge}
        }
    })
}
type IntersectionPointInfo = {
    p: V3, // intersection point
    insideDir: V3, // currently not needed todo
    t: number, // param on intersection curve
    edge: Edge, // face edge doing the intersection
    edgeT: number,
    colinear: boolean, // whether edge is colinear to intersection line
    used?: boolean }

function planeFaceEdgeISPsWithPlane(face: PlaneFace, isLine: L3, plane2: P3): IntersectionPointInfo[] {
    assert(face.surface.plane.containsLine(isLine))
    assert(plane2.containsLine(isLine))
    const plane = face.surface.plane
    const ps = []
    const loops = [face.edges].concat(face.holes)
    loops.forEach(loop => {
        const colinearEdges = loop.map((edge) => edge.colinearToLine(isLine) && -sign(edge.aDir.dot(isLine.dir1)))
        const isLineOut = isLine.dir1.cross(plane.normal)

        loop.forEach((edge, edgeIndex, edges) => {
            const nextEdgeIndex = (edgeIndex + 1) % edges.length, nextEdge = edges[nextEdgeIndex]
            //console.log(edge.toSource()) {p:V3(2, -2.102, 0),
            if (colinearEdges[edgeIndex]) {
                // edge colinear to intersection line
                const curveAT = isLine.pointLambda(edge.a), curveBT = isLine.pointLambda(edge.b)
                // add interval for colinear segment
                ps.push(
                    {p: edge.a, insideDir: edge.aDir,           t: curveAT, edge: edge, edgeT: edge.aT, colinear: true},
                    {p: edge.b, insideDir: edge.bDir.negated(), t: curveBT, edge: edge, edgeT: edge.bT, colinear: true})
                // open next interval if necessary
                const nextSide = colinearEdges[nextEdgeIndex] || dotCurve(isLineOut, nextEdge.aDir, nextEdge.aDDT)
                if (colinearEdges[edgeIndex] * nextSide < 0) {
                    // side changes
                    ps.push({p: nextEdge.a, insideDir: edge.bDir, t: curveBT, edge: nextEdge, edgeT: nextEdge.aT, colinear: false})
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
                            ps.push({p: edge.b, insideDir: plane2.normal.negated(), t: isLine.pointLambda(edge.b), edge: edge, edgeT: edge.bT, colinear: false})
                        }
                    } else if (edgeT != edge.aT) {
                        // edge crosses intersection line, neither starts nor ends on it
                        const p = edge.curve.at(edgeT)
                        assert(plane2.containsPoint(p), edge.toString(), p, edgeT, plane2.distanceToPoint(p))
                        assert(isLine.containsPoint(p), edge.toString(), p, edgeT, isLine.distanceToPoint(p))
                        const insideDir = plane2.normal.negated()
                        ps.push({p: p, insideDir: insideDir, t: isLine.pointLambda(p), edge: edge, edgeT: edgeT, colinear: false})
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
function planeFaceEdgeISPsWithPlane2(face: PlaneFace, isLine: L3, plane2: P3): IntersectionPointInfo[] {
    assert(face.surface.plane.containsLine(isLine))
    assert(plane2.containsLine(isLine))
    let plane = face.surface.plane
    let ps = []
    let loops = [face.edges].concat(face.holes)
    loops.forEach(loop => {
        let colinearEdges = loop.map((edge) => edge.colinearToLine(isLine))
        let isLineOut = isLine.dir1.cross(plane.normal)

        loop.forEach((edge, edgeIndex, edges) => {
            let nextEdgeIndex = (edgeIndex + 1) % edges.length, nextEdge = edges[nextEdgeIndex]
            //console.log(edge.toSource()) {p:V3(2, -2.102, 0),
            if (colinearEdges[edgeIndex]) {
                // edge colinear to intersection line
                let prevEdgeIndex = (edgeIndex - 1 + edges.length) % edges.length, prevEdge = edges[prevEdgeIndex]
                // if colinear, edge.curve must be a line, so colinearOut is constant
                let colinearOut = edge.aDir.cross(plane.normal)
                let curveAT = isLine.pointLambda(edge.a), curveBT = isLine.pointLambda(edge.b)
                // close previous interval if necessary
                if (!colinearEdges[prevEdgeIndex] && dotCurve(colinearOut, prevEdge.bDir, prevEdge.bDDT) < 0) {
                    ps.push({p: prevEdge.b, insideDir: edge.aDir.negated(), t: curveAT, edge: prevEdge, edgeT: prevEdge.bT,
                        colinear: false})
                }
                // add interval for colinear segment
                ps.push(
                    {p: edge.a, insideDir: edge.aDir, t: curveAT, edge: edge, edgeT: edge.aT,
                        colinear: true},
                    {p: edge.b, insideDir: edge.bDir.negated(), t: curveBT, edge: edge, edgeT: edge.bT,
                        colinear: true})
                // open next interval if necessary
                if (!colinearEdges[nextEdgeIndex] && dotCurve(colinearOut, nextEdge.aDir, nextEdge.aDDT) > 0) {
                    ps.push({p: nextEdge.a, insideDir: edge.bDir, t: curveBT, edge: nextEdge, edgeT: nextEdge.aT, colinear: false})
                }
            } else {
                // not necessarily a straight edge, so multiple intersections are possible
                let edgeTs = edge.edgeISTsWithPlane(plane2)
                assert(edgeTs.every(t => plane2.containsPoint(edge.curve.at(t))), edgeTs)
                for (let k = 0; k < edgeTs.length; k++) {
                    let edgeT = edgeTs[k]
                    if (edgeT == edge.bT) {
                        // endpoint lies on intersection line
                        if (!colinearEdges[nextEdgeIndex]
                            && edgesDifferentSidesOfDir(isLineOut, edge, nextEdge)) {
                            // next segment is not colinear and ends on different side
                            ps.push({p: edge.b, insideDir: plane2.normal.negated(), t: isLine.pointLambda(edge.b), edge: edge, edgeT: edge.bT, colinear: false})
                        }
                    } else if (edgeT != edge.aT) {
                        // edge crosses intersection line, neither starts nor ends on it
                        let p = edge.curve.at(edgeT)
                        assert(plane2.containsPoint(p), edge.toString(), p, edgeT, plane2.distanceToPoint(p))
                        assert(isLine.containsPoint(p), edge.toString(), p, edgeT, isLine.distanceToPoint(p))
                        let insideDir = plane2.normal.negated()
                        ps.push({p: p, insideDir: insideDir, t: isLine.pointLambda(p), edge: edge, edgeT: edgeT, colinear: false})
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
function faceEdgeISPsWithSurface(face: Face, isCurves: Curve[], surface2: Surface): IntersectionPointInfo[][] {
    const surface = face.surface
    const pss = NLA.arrayFromFunction(isCurves.length, i => [])

    const loops = face.holes.concat([face.edges])
    loops.forEach(loop => {
        const colinearEdges: int[] = loop.map(edge => isCurves.findIndex(curve => edge.curve.isColinearTo(curve)))
        const colinearSides = loop.map((edge, edgeIndex) => -1 != colinearEdges[edgeIndex]
                    && -sign(isCurves[colinearEdges[edgeIndex]].tangentAt(edge.aT).dot(edge.aDir)))
// todo: this assumes that isCurves do not touch
        loop.forEach((edge, edgeIndex, edges) => {
            let nextEdgeIndex = (edgeIndex + 1) % edges.length, nextEdge = edges[nextEdgeIndex]
            //console.log(edge.toSource()) {p:V3(2, -2.102, 0),
            let colinearEdgeCurveIndex = colinearEdges[edgeIndex]
            if (-1 != colinearEdgeCurveIndex) {
                let isCurve = isCurves[colinearEdgeCurveIndex]
                // edge colinear to an intersection curve
                let prevEdgeIndex = (edgeIndex - 1 + edges.length) % edges.length, prevEdge = edges[prevEdgeIndex]
                const bHint = isCurve instanceof EllipseCurve
                    && sign(isCurve.tangentAt(isCurve.pointLambda(edge.a, PI)).dot(edge.aDir)) * PI
                const curveAT = isCurve.pointLambda(edge.a, -bHint), curveBT = isCurve.pointLambda(edge.b, bHint)
                let colinearOutA = edge.aDir.cross(surface.normalAt(edge.a))
                let ps = pss[colinearEdgeCurveIndex]
                if (-1 != colinearEdges[prevEdgeIndex] && dotCurve2(prevEdge.curve, prevEdge.bT, colinearOutA, -1) > 0) {
                //if (-1 != colinearEdges[prevEdgeIndex] && dotCurve(colinearOutA, prevEdge.bDir, prevEdge.bDDT) < 0) {
                    ps.push({p: prevEdge.b, insideDir: edge.aDir.negated(), t: curveAT, edge: prevEdge, edgeT: prevEdge.bT, colinear: false})
                }
                ps.push(
                    {p: edge.a, insideDir: edge.aDir, t: curveAT, edge: edge, edgeT: edge.aT, colinear: true},
                    {p: edge.b, insideDir: edge.bDir.negated(), t: curveBT, edge: edge, edgeT: edge.bT, colinear: true})
                let colinearOutB = edge.bDir.cross(surface.normalAt(edge.b))
                if (-1 != colinearEdges[nextEdgeIndex] && dotCurve2(nextEdge.curve, nextEdge.aT, colinearOutB, 1) > 0) {
                //if (-1 != colinearEdges[nextEdgeIndex] && dotCurve(colinearOutB, nextEdge.aDir, nextEdge.aDDT) > 0) {
                    ps.push({p: edge.b, insideDir: edge.bDir, t: curveBT, edge: edge, edgeT: prevEdge.bT, colinear: false})
                }
            } else {
                const edgeTs = edge.edgeISTsWithSurface(surface2)
                for (let k = 0; k < edgeTs.length; k++) {
                    const edgeT = edgeTs[k]
                    let p = edge.curve.at(edgeT)
                    if (!isCurves.some(isCurve => isCurve.containsPoint(p), edge.toString() + p+edgeT)) {
                        console.log(isCurves, isCurves[0].sce, p.sce, edge.sce, surface2.sce, face.surface.sce)
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
                    assert(!(isCurve instanceof EllipseCurve) || PI != isCurve.pointLambda(p, PI) || !eq0(edge.curve.tangentAt(edgeT).cross(isCurve.f1).dot(isCurve.f2)))
                    const hint = isCurve instanceof EllipseCurve &&
                        sign(edge.curve.tangentAt(edgeT).cross(isCurve.f1).dot(isCurve.f2)) * PI
                    const curveT = isCurve.pointLambda(p, hint)
                    assert(!isNaN(curveT))
                    const insideDir = edge.tangentAt(edgeT).cross(surface.normalAt(p)).negated()
                    assert(!NLA.eq0(insideDir.dot(isCurve.tangentAt(curveT))))
                    // Edge.edgeISTsWithSurface returns snapped values, so comparison with == is ok:
                    if (edgeT == edge.bT) {
                        // endpoint lies on intersection line
                        if (-1 != colinearEdges[nextEdgeIndex]) {
                            if (edgesDifferentSidesOfDir(isCurve.tangentAt(curveT), edge, nextEdge)) {
                                // next segment is not colinear and ends on different side
                                console.log("adding")
                                pss[isCurveIndex].push({ p: edge.b, insideDir: insideDir, t: curveT, edge: edge, edgeT: edge.bT, colinear: false})
                            }
                        }
                    } else if (edgeT != edge.aT) {
                        // edge crosses an intersection curve, neither starts nor ends on it
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
    return pss
}
function dotCurve(v: V3, cDir: V3, cDDT: V3): number {
	let dot = v.dot(cDir)
	if (NLA.eq0(dot)) { dot = v.dot(cDDT) }
	assert(!NLA.eq0(dot))
	return dot
}
function dotCurve2(curve: Curve, t: number, normal: V3, sign: number): number {
	assert(sign == 1 || sign == -1, sign)
	const tangentDot = curve.tangentAt(t).dot(normal)
	// if tangentDot != 0 the curve simply crosses the plane
	if (!NLA.eq0(tangentDot)) { return sign * tangentDot }
	const ddtDot = curve.ddt(t).dot(normal)
	// tangentDot == 0 ==> critical point at t, if ddtDot != 0, then it is a turning point, otherwise we can't be sure and must do a numeric test
	if (!NLA.eq0(ddtDot)) { return ddtDot }
	const numericDot = curve.at(t).to(curve.at(t + sign * 4 * NLA_PRECISION)).dot(normal)
	assert(!(curve instanceof L3))
	return numericDot
}
function edgesDifferentSidesOfDir(dir: V3, e1: Edge, e2: Edge): boolean {
    let factor1 = dir.dot(e1.bDir)
    if (NLA.eq0(factor1)) { factor1 = dir.dot(e1.bDDT) }
    assert(!NLA.eq0(factor1))
    let factor2 = dir.dot(e2.aDir)
    if (NLA.eq0(factor2)) { factor2 = dir.dot(e2.aDDT) }
    assert(!NLA.eq0(factor2))
    return factor1 * factor2 > 0
}


const INSIDE = 0, OUTSIDE = 1, COPLANAR_SAME = 2, COPLANAR_OPPOSITE= 3, ALONG_EDGE_OR_PLANE = 4
/**
 *
 * @param brep BREP to check
 * @param edge edge to check
 * @param dirAtEdgeA the direction vector to check
 * @param faceNormal If dirAtEdgeA doesn't split a volume, but is along a face, the returned value depends on
 *     wether that faces normal points in the same direction as faceNormal
 * @returns INSIDE, OUTSIDE, COPLANAR_SAME or COPLANAR_OPPOSITE
 */
//function splitsVolumeEnclosingFaces(brep: B2, edge: Edge, dirAtEdgeA: V3, faceNormal: V3): int {
//    assert(arguments.length == 4)
//    //assert(p.equals(edge.a))
//    const ab1 = edge.aDir.normalized()
//    const relFaces = facesWithEdge(edge, brep.faces) as any[]
//    relFaces.forEach(faceInfo => {
//        faceInfo.normalAtEdgeA = faceInfo.face.surface.normalAt(edge.a)
//        faceInfo.edgeDirAtEdgeA = !faceInfo.reversed
//            ? faceInfo.edge.aDir
//            : faceInfo.edge.bDir
//        faceInfo.outsideVector = faceInfo.edgeDirAtEdgeA.cross(faceInfo.normalAtEdgeA)
//        faceInfo.angle = (dirAtEdgeA.angleRelativeNormal(faceInfo.outsideVector.negated(), ab1) + 2 * Math.PI + NLA_PRECISION / 2) % (2 * Math.PI)
//    })
//    assert(relFaces.length != 0, edge.toSource())
//    relFaces.sort((a, b) => a.angle - b.angle)
//    // assert(relFaces.length % 2 == 0, edge.toSource()) // even number of touching faces
//
//    if (NLA.eq0(relFaces[0].angle)) {
//        //assert(false) todo
//        const coplanarSame = relFaces[0].normalAtEdgeA.dot(faceNormal) > 0;
//        return coplanarSame ? COPLANAR_SAME : COPLANAR_OPPOSITE
//    } else {
//        return !relFaces[0].reversed ? INSIDE : OUTSIDE
//    }
//}
function splitsVolumeEnclosingFaces(brep: B2, canonEdge: Edge, dirAtEdgeA: V3, faceNormal: V3): int {
    assert(arguments.length == 4)
    assert(canonEdge == canonEdge.getCanon())
    //assert(p.equals(canonEdge.a))
    const edgeFaceInfos = brep.edgeFaces.get(canonEdge) as any[]
    assertf(() => edgeFaceInfos.length % 2 == 0)
    assertf(() => brep.edgeFaces)
    const faceInfo0 = edgeFaceInfos[0]
    const aDir1 = canonEdge.aDir.normalized()
    const angleToCanon = (faceInfo0.inside.angleRelativeNormal(dirAtEdgeA, aDir1) + 2 * Math.PI + NLA_PRECISION) % (2 * Math.PI) - NLA_PRECISION
    const nearestFaceInfoIndex = edgeFaceInfos.findIndex(faceInfo => NLA.lt(angleToCanon, faceInfo.angle))
	const nearestFaceInfo = edgeFaceInfos[nearestFaceInfoIndex == -1 ? edgeFaceInfos.length - 1 : nearestFaceInfoIndex - 1]
    if (NLA.eq(nearestFaceInfo.angle, angleToCanon)) {
        //assert(false) todo
        const coplanarSame = nearestFaceInfo.normalAtCanonA.dot(faceNormal) > 0
        return coplanarSame ? COPLANAR_SAME : COPLANAR_OPPOSITE
    } else {
        return nearestFaceInfo.reversed ? INSIDE : OUTSIDE
    }
}
function splitsVolumeEnclosingCone(brep: B2, p: V3, dir: V3) {
    const testPlane = P3.forAnchorAndPlaneVectors(p, dir, dir.getPerpendicular())
    const rays = []
    for (let k = 0; k < brep.faces.length; k++) {
        const face = brep.faces[k] as PlaneFace
        assertf(() => face instanceof PlaneFace)
        if (face.getAllEdges().some(edge => edge.a.like(p))) {
            if (testPlane.isParallelToPlane(face.surface.plane)) {
                if (face.pointsToInside(p, dir) != PointVsFace.OUTSIDE) {
                    return ALONG_EDGE_OR_PLANE
                }
            } else {
                const isLine = L3.fromPlanes(testPlane, face.surface.plane)
                const ps = planeFaceEdgeISPsWithPlane(face, isLine, testPlane)
                let i = 0
                while (i < ps.length) {
                    const a = ps[i++], b = ps[i++]
                    const out = a.p.like(p)
                    if (out || b.p.like(p)) {
                        const dir2 = out ? isLine.dir1 : isLine.dir1.negated()
                        const angle = (dir.angleRelativeNormal(dir2, testPlane.normal) + 2 * Math.PI + NLA_PRECISION / 2) % (2 * Math.PI)
                        rays.push({angle: angle, out: out})
                    }
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
function fff(info: {face: Face, edge: Edge, normalAtCanonA: V3, inside: V3, reversed: boolean, angle: number}, surface: Surface): int {
    const canonA = info.edge.reversed ? info.edge.b : info.edge.a
    const surfaceNormalAtCanonA = surface.normalAt(canonA)
    const dot = NLA.snap0(info.inside.dot(surfaceNormalAtCanonA))
    if (0 !== dot) {
        return 0 < dot ? OUTSIDE : INSIDE
    }
    if (surface.isCoplanarTo(info.face.surface)) {
        return 0 < info.normalAtCanonA.dot(surfaceNormalAtCanonA) ? COPLANAR_SAME : COPLANAR_OPPOSITE
    }
    assert(false)
}

/**
 * Solves a quadratic system of equations of the form
 *      a * x + b * y = c
 *      x^2 + y^2 = 1
 * This can be understood as the intersection of the unit circle with a line.
 *
 * a * b + (b -c) * (b + c)
 */
function intersectionUnitCircleLine(a: number, b: number, c: number): {x1: number, y1: number, x2: number, y2: number} {
    assertNumbers(a, b, c)
    // TODO: disambiguate on a < b
    const term = sqrt(a * a + b * b - c * c)
    return {
        x1: (a * c + b * term) / (a * a + b * b),
        x2: (a * c - b * term) / (a * a + b * b),
        y1: (b * c - a * term) / (a * a + b * b),
        y2: (b * c + a * term) / (a * a + b * b)
    }
}
function intersectionCircleLine(a: number, b: number, c: number, r: number): {x1: number, x2: number, y1: number, y2: number} {
    assertNumbers(a, b, c, r)
    const term = sqrt(r * r * (a * a + b * b) - c * c)
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
 * @returns with x1 >= x2 and y1 <= y2
 * a * b + (b -c) * (b + c)
 */
function intersectionUnitHyperbolaLine(a: number, b: number, c: number): { x1: number, y1: number, x2: number, y2: number } {
    assertNumbers(a, b, c)
    const aa = a * a, bb = b * b, cc = c * c
    // TODO: disambiguate on a < b
    //var xTerm = sqrt(4*cc*aa-4*(bb-aa)*(-cc-bb))
    const xTerm = 2 * sqrt(bb * cc + bb * bb - aa * bb)
    const yTerm = sqrt(4 * cc * bb - 4 * (bb - aa) * (cc - aa))
    return {
        x1: (-2 * a * c + xTerm) / 2 / (bb - aa),
        x2: (-2 * a * c - xTerm) / 2 / (bb - aa),
        y1: (2 * b * c - yTerm) / 2 / (bb - aa),
        y2: (2 * b * c + yTerm) / 2 / (bb - aa)
    }
}




function curvePoint(implicitCurve, startPoint) {
    const eps = 1e-5
    let p = startPoint
    for (let i = 0; i < 4; i++) {
        const fp = implicitCurve(p.x, p.y)
        const dfpdx = (implicitCurve(p.x + eps, p.y) - fp) / eps,
            dfpdy = (implicitCurve(p.x, p.y + eps) - fp) / eps
        const scale = fp / (dfpdx * dfpdx + dfpdy * dfpdy)
        //console.log(p.$)
        p = p.minus(new V3(scale * dfpdx, scale * dfpdy, 0))
    }
    return p
}
function followAlgorithm (implicitCurve, startPoint, endPoint, stepLength, startDir, tangentEndPoints, boundsFunction) {
    assertNumbers(stepLength, implicitCurve(0, 0))
    assertVectors(startPoint, endPoint)
    assert (!startDir || startDir instanceof V3)
    const points = []
    tangentEndPoints = tangentEndPoints || []
    assert (NLA.eq0(implicitCurve(startPoint.x, startPoint.y)), 'NLA.isZero(implicitCurve(startPoint.x, startPoint.y))')
    stepLength = stepLength || 0.5
    const eps = 1e-5
    let p = startPoint, prevp = startDir ? p.minus(startDir) : p
    let i = 0
    do {
        const fp = implicitCurve(p.x, p.y)
        const dfpdx = (implicitCurve(p.x + eps, p.y) - fp) / eps,
            dfpdy = (implicitCurve(p.x, p.y + eps) - fp) / eps
        let tangent = new V3(-dfpdy, dfpdx, 0)
        const reversedDir = p.minus(prevp).dot(tangent) < 0
        tangent = tangent.toLength(reversedDir ? -stepLength : stepLength)
        const tangentEndPoint = p.plus(tangent)
        points.push(p)
        tangentEndPoints.push(tangentEndPoint)
        prevp = p
        p = curvePoint(implicitCurve, tangentEndPoint)
    } while (i++ < 100 && (i < 4 || prevp.distanceTo(endPoint) > 1.1 * stepLength) && boundsFunction(p.x, p.x))
    // TODO gleichmäßige Verteilung der Punkte
    return points
}
// both curves must be in the same s-t coordinates for this to make sense
function intersectionICurveICurve(
    iCurve1: (s: number, t: number) => number,
    startParams1: V3,
    endParams1: V3,
    startDir,
    stepLength: number,
    iCurve2: (s: number, t: number) => number) {

    assertNumbers(stepLength, iCurve1(0, 0), iCurve2(0, 0))
    assertVectors(startParams1, endParams1)
    assert (!startDir || startDir instanceof V3)
    const vertices = []
    assert (NLA.eq0(iCurve1(startParams1.x, startParams1.y)))
    stepLength = stepLength || 0.5
    const eps = 1e-5
    let p = startParams1, prevp = p // startDir ? p.minus(startDir) : p
    let i = 0
    while (i++ < 1000 && (i < 4 || p.distanceTo(endParams1) > 1.1 * stepLength)) {
        const fp = iCurve1(p.x, p.y)
        const dfpdx = (iCurve1(p.x + eps, p.y) - fp) / eps,
            dfpdy = (iCurve1(p.x, p.y + eps) - fp) / eps
        let tangent = new V3(-dfpdy, dfpdx, 0).toLength(stepLength)
        if (p.minus(prevp).dot(tangent) < 0) tangent = tangent.negated()
        prevp = p
        p = curvePoint(iCurve1, p.plus(tangent))
        vertices.push(p)
    }
    // TODO gleichmäßige Verteilung der Punkte
    return vertices

}
function intersectionICurveICurve2(iCurve1, loopPoints1, iCurve2) {
    let p = loopPoints1[0], val = iCurve2(p.x, p.y), lastVal
    const iss = []
    for (let i = 0; i < loopPoints1.length; i++) {
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
    const iss = []
    let val = implicitSurface(parametricCurve(searchStart)), lastVal
    for (let t = searchStart + searchStep; t <= searchEnd; t += searchStep) {
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
    return (x, y) => (x*x + y*y)*(x*x + y*y) - 2*c*c*(x*x - y*y) - (a*a*a*a - c*c*c*c)
}



