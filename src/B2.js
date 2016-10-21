var __extends = (this && this.__extends) || function (d, b) {
    for (var p in b) if (b.hasOwnProperty(p)) d[p] = b[p];
    function __() { this.constructor = d; }
    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
};
var eps = 1e-5;
var Face = (function (_super) {
    __extends(Face, _super);
    function Face(surface, contour, holes, name) {
        _super.call(this);
        Edge.assertLoop(contour);
        //assert(surface.edgeLoopCCW(contour), surface.toString()+contour.join("\n"))
        assert(contour.every(function (f) { return f instanceof Edge; }), 'contour.every(f => f instanceof Edge)' + contour.toSource());
        // contour.forEach(e => !surface.containsCurve(e.curve) &&
        // console.log("FAIL:"+surface.distanceToPoint(e.curve.anchor)))
        contour.forEach(function (e) { return assert(surface.containsCurve(e.curve), '' + e + surface); });
        holes && holes.forEach(function (hole) { return Edge.assertLoop(hole); });
        //holes && holes.forEach(hole => assert(!surface.edgeLoopCCW(hole)))
        assert(!holes || holes.constructor == Array, holes && holes.toString());
        this.surface = surface;
        this.edges = contour; // TODO refactor to contour
        this.holes = holes || [];
        this.id = globalId++;
        this.name = name;
    }
    Face.prototype.transform = function (m4) {
        var newEdges = this.edges.map(function (e) { return e.transform(m4); });
        var newHoles = this.holes.map(function (hole) { return hole.map(function (e) { return e.transform(m4); }); });
        return new this.constructor(this.surface.transform(m4), newEdges, newHoles);
    };
    Face.prototype.flipped = function () {
        var newEdges = this.edges.map(function (e) { return e.flipped(); }).reverse();
        var newHoles = this.holes.map(function (hole) { return hole.map(function (e) { return e.flipped(); }).reverse(); });
        return new this.constructor(this.surface.flipped(), newEdges, newHoles, this.name);
    };
    Face.prototype.toString = function () {
        return 'new ' + this.constructor.name + '(' + this.surface + ', ' + this.edges.map(function (e) { return '\n\t' + e; }).join() + ']'
            + this.holes.map(function (hole) { return '\n\t\thole: ' + hole.join(); }) + ')';
    };
    Face.prototype.toSource = function () {
        return "new " + this.constructor.name + "(" + this.surface.toSource() + ", [" + this.edges.map(function (e) { return '\n\t' + e.toSource(); }).join(',') + "], [" + this.holes.map(function (hole) { return '[' + hole.map(function (e) { return '\n\t' + e.toSource(); }).join(',') + ']'; }).join(',') + "])";
    };
    Face.prototype.equals = function (face) {
        var _this = this;
        //TODO		assert(false)
        var edgeCount = this.edges.length;
        return this.surface.equals(face.surface) &&
            this.edges.length == face.edges.length &&
            NLA.arrayRange(0, edgeCount, 1)
                .some(function (offset) { return _this.edges.every(function (edge, i) { return edge.equals(face.edges[(offset + i) % edgeCount]); }); });
    };
    Face.prototype.likeFace = function (face2) {
        function loopsLike(a, b) {
            return a.length == b.length &&
                NLA.arrayRange(0, a.length, 1)
                    .some(function (offset) { return a.every(function (edge, i) { return edge.like(b[(offset + i) % a.length]); }); });
        }
        assertInst(Face, face2);
        return this.surface.like(face2.surface)
            && this.holes.length == face2.holes.length
            && loopsLike(this.edges, face2.edges)
            && this.holes.every(function (hole) { return face2.holes.some(function (hole2) { return loopsLike(hole, hole2); }); });
    };
    Face.prototype.getAllEdges = function () {
        return Array.prototype.concat.apply(this.edges, this.holes);
    };
    Face.prototype.addEdgeLines = function (mesh) {
        assert(false, "buggy, fix");
        var vertices = this.edges.map(function (edge) { return edge.getVerticesNo0(); }).concatenated(), mvl = mesh.vertices.length;
        for (var i = 0; i < vertices.length; i++) {
            mesh.vertices.push(vertices[i]);
            mesh.lines.push(mvl + i, mvl + (i + 1) % vertices.length);
        }
    };
    Face.prototype.inB2 = function () {
        return new B2([this]);
    };
    Face.prototype.containsPoint = function (p) {
        var _this = this;
        assertVectors(p);
        return this.surface.edgeLoopContainsPoint(this.edges, p) != PointVsFace.OUTSIDE
            && !this.holes.some(function (hole) { return _this.surface.edgeLoopContainsPoint(hole, p) != PointVsFace.OUTSIDE; });
    };
    Face.prototype.containsPoint2 = function (p) {
        assertVectors(p);
        var contourContainsPoint = this.surface.edgeLoopContainsPoint(this.edges, p);
        if (contourContainsPoint != PointVsFace.INSIDE)
            return contourContainsPoint;
        for (var _i = 0, _a = this.holes; _i < _a.length; _i++) {
            var hole = _a[_i];
            var loopContainsPoint = this.surface.edgeLoopContainsPoint(hole, p);
            if (loopContainsPoint != PointVsFace.OUTSIDE) {
                return loopContainsPoint == PointVsFace.ON_EDGE ? PointVsFace.ON_EDGE : PointVsFace.OUTSIDE;
            }
        }
        return PointVsFace.INSIDE;
    };
    /**
     *
     * @param line
     * @returns t param of the line if there is an intersection, NaN otherwise
     */
    Face.prototype.intersectsLine = function (line) {
        var _this = this;
        assertInst(L3, line);
        var containedIntersectionsTs = this.surface.isTsForLine(line).filter(function (t) { return _this.containsPoint(line.at(t)); });
        var nearestPointT = containedIntersectionsTs.withMax(function (t) { return -t; });
        return undefined != nearestPointT ? nearestPointT : NaN;
    };
    Face.prototype.toMesh = function () {
        var mesh = new GL.Mesh({ triangles: true, normals: true, lines: true });
        this.addToMesh(mesh);
        mesh.compile();
        return mesh;
    };
    Face.prototype.getLoops = function () {
        return this.holes.concat(this.edges);
    };
    return Face;
}(Transformable));
var PlaneFace = (function (_super) {
    __extends(PlaneFace, _super);
    function PlaneFace(planeSurface, contour, holes, name) {
        assertInst(PlaneSurface, planeSurface);
        _super.call(this, planeSurface, contour, holes, name);
    }
    PlaneFace.prototype.zDirVolume = function () {
        var _a = this.calculateArea(), centroid = _a.centroid, area = _a.area;
        return { volume: this.surface.plane.normal.z * centroid.z * area,
            centroid: new V3(centroid.x, centroid.y, centroid.z / 2) };
    };
    PlaneFace.prototype.calculateArea = function () {
        var _this = this;
        var centroid = V3.ZERO, tcs = 0, tct = 0, totalArea = 0;
        var r1 = this.surface.right, u1 = this.surface.up;
        this.edges.forEach(function (edge) {
            var edgeCentroid, edgeArea, centroidS, centroidT;
            if (edge instanceof StraightEdge) {
                var midPoint = edge.a.lerp(edge.b, 0.5);
                edgeCentroid = new V3(midPoint.x, centroid.y, centroid.z / 2);
                centroidS = midPoint.dot(r1) / 2;
                centroidT = midPoint.dot(u1);
                var edgeLength = edge.a.distanceTo(edge.b);
                edgeArea = edgeLength * edge.curve.dir1.dot(r1);
            }
            else {
                var curve = edge.curve;
                if (curve instanceof EllipseCurve) {
                    var info = curve.getAreaInDir(u1, r1, edge.aT, edge.bT);
                    edgeArea = info.area;
                    var parametricCentroid = _this.surface.pointToParameterFunction()(info.centroid);
                    centroidS = parametricCentroid.x;
                    centroidT = parametricCentroid.y;
                }
                else if (curve instanceof BezierCurve) {
                    edgeArea = curve.getAreaInDir(u1, r1, edge.aT, edge.bT);
                }
                else {
                    assertNever();
                }
            }
            tcs += edgeArea * centroidS;
            tct += edgeArea * centroidT;
            totalArea += edgeArea;
        });
        centroid = r1.times(tcs).plus(u1.times(tct));
        assertNever();
    };
    PlaneFace.prototype.addToMesh = function (mesh) {
        var mvl = mesh.vertices.length;
        var normal = this.surface.plane.normal;
        var vertices = this.edges.flatMap(function (edge) { return edge.getVerticesNo0(); });
        for (var i = 0; i < vertices.length; i++) {
            mesh.lines.push(mvl + i, mvl + (i + 1) % vertices.length);
        }
        var holeStarts = [];
        this.holes.forEach(function (hole) {
            holeStarts.push(vertices.length);
            vertices.pushAll(hole.flatMap(function (edge) { return edge.getVerticesNo0(); }));
        });
        var triangles = triangulateVertices(normal, vertices, holeStarts).map(function (index) { return index + mvl; });
        Array.prototype.push.apply(mesh.vertices, vertices);
        Array.prototype.push.apply(mesh.triangles, triangles);
        Array.prototype.push.apply(mesh.normals, NLA.arrayFromFunction(vertices.length, function (i) { return normal; }));
    };
    PlaneFace.prototype.intersectsLine = function (line) {
        assertInst(L3, line);
        var lambda = line.intersectWithPlaneLambda(this.surface.plane);
        if (!Number.isFinite(lambda)) {
            return NaN;
        }
        var inside = this.containsPoint(line.at(lambda));
        return inside ? lambda : NaN;
    };
    PlaneFace.prototype.withHole = function (holeEdges) {
        return new PlaneFace(this.surface, this.edges, [holeEdges]);
    };
    PlaneFace.prototype.intersectFace = function (face2, thisBrep, face2Brep, faceMap, edgeMap, likeSurfaceFaces) {
        if (face2 instanceof PlaneFace) {
            this.intersectPlaneFace(face2, thisBrep, face2Brep, faceMap, edgeMap, likeSurfaceFaces);
            return;
        }
        if (face2 instanceof RotationFace) {
            face2.intersectFace(this, face2Brep, thisBrep, faceMap, edgeMap, likeSurfaceFaces);
            return;
        }
        assert(false);
    };
    PlaneFace.prototype.intersectPlaneFace = function (face2, thisBrep, face2Brep, faceMap, edgeMap, likeSurfaceFaces) {
        assertInst(PlaneFace, face2);
        var face = this;
        // get intersection
        var thisPlane = this.surface.plane, face2Plane = face2.surface.plane;
        if (thisPlane.isParallelToPlane(face2Plane)) {
            if (thisPlane.like(face2Plane)) {
                // normal same and same location in space
                addLikeSurfaceFaces(likeSurfaceFaces, this, face2);
            }
            return;
        }
        var intersectionLine = L3.fromPlanes(thisPlane, face2Plane);
        var thisDir = true;
        // get intersections of newCurve with other edges of face and face2
        var ps1 = planeFaceEdgeISPsWithPlane(thisBrep, face, intersectionLine, face2Plane);
        var ps2 = planeFaceEdgeISPsWithPlane(face2Brep, face2, intersectionLine, thisPlane);
        console.log('ps1\n', ps1.map(function (m) { return m.toSource(); }).join('\n'), '\nps2\n', ps2.map(function (m) { return m.toSource(); }).join('\n'));
        if (ps1.length == 0 || ps2.length == 0) {
            // faces to not intersect
            return;
        }
        console.log('' + thisPlane + face2Plane);
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
        function handleNewEdge(seg, col1, col2, in1, in2, a, b) {
            console.log("handling seg", col1 && col1.toSource(), col2 && col2.toSource(), seg.toSource());
            if (!col1 && !col2) {
                console.log("adding");
                NLA.mapAdd(faceMap, face, seg);
                NLA.mapAdd(faceMap, face2, seg.flipped());
            }
            if (!col1 && col2) {
                if (col2.aDir.cross(face2Plane.normal).dot(thisPlane.normal) > 0) {
                    // NB: a new edge is inserted even though it may be the same as an old one
                    // however it indicates that it intersects the other volume here, i.e. the old edge cannot
                    // be counted as "inside" for purposes of reconstitution
                    //NLA.mapAdd(faceMap, face2, seg.flipped())
                    NLA.mapAdd(edgeMap, col2, { edgeT: col2.curve.pointLambda(seg.a), p: seg.a });
                    NLA.mapAdd(edgeMap, col2, { edgeT: col2.curve.pointLambda(seg.b), p: seg.b });
                    var testVector = face2Plane.normal.negated().rejectedFrom(thisPlane.normal);
                    var sVEF = splitsVolumeEnclosingFaces(face2Brep, col2, testVector, thisPlane.normal);
                    if (sVEF == INSIDE) {
                        NLA.mapAdd(faceMap, face, seg);
                    }
                }
            }
            if (col1 && !col2) {
                if (col1.aDir.cross(thisPlane.normal).dot(face2Plane.normal) > 0) {
                    NLA.mapAdd(edgeMap, col1, { edgeT: col1.curve.pointLambda(seg.a), p: seg.a });
                    NLA.mapAdd(edgeMap, col1, { edgeT: col1.curve.pointLambda(seg.b), p: seg.b });
                    var testVector = thisPlane.normal.negated().rejectedFrom(face2Plane.normal);
                    var sVEF = splitsVolumeEnclosingFaces(thisBrep, col1, testVector, face2Plane.normal);
                    if (sVEF == INSIDE) {
                        NLA.mapAdd(faceMap, face2, seg.flipped());
                    }
                }
            }
            if (col1 && col2) {
                // faces are touching edge-edge
                var testVector1 = col1.aDir.cross(thisPlane.normal).negated();
                var testVector2 = col2.aDir.cross(face2Plane.normal).negated();
                console.log("testVector", testVector1.sce, face2Brep.sce, splitsVolumeEnclosingFaces(face2Brep, col2, testVector1, thisPlane.normal) == INSIDE);
                var sVEF1 = splitsVolumeEnclosingFaces(face2Brep, col2, testVector1, thisPlane.normal);
                var sVEF2 = splitsVolumeEnclosingFaces(thisBrep, col1, testVector2, face2Plane.normal);
                // todo this doesn't work when the only faces which split VEFs of the other volume are coplanar to those faces, so never actually get paired
                //if ((sVEF1 == INSIDE || sVEF1 == COPLANAR_OPPOSITE)
                //	&& (sVEF2 == INSIDE || sVEF2 == COPLANAR_OPPOSITE)) {
                //	//NLA.mapAdd(faceMap, face2, seg.flipped())
                //	NLA.mapAdd(edgeMap, col1, {edgeT: col1.curve.pointLambda(seg.a), p: seg.a})
                //	NLA.mapAdd(edgeMap, col1, {edgeT: col1.curve.pointLambda(seg.b), p: seg.b})
                //
                //	NLA.mapAdd(edgeMap, col2, {edgeT: col2.curve.pointLambda(seg.a), p: seg.a})
                //	NLA.mapAdd(edgeMap, col2, {edgeT: col2.curve.pointLambda(seg.b), p: seg.b})
                //}
                if (sVEF1 == INSIDE) {
                    NLA.mapAdd(edgeMap, col1, { edgeT: col1.curve.pointLambda(seg.a), p: seg.a });
                    NLA.mapAdd(edgeMap, col1, { edgeT: col1.curve.pointLambda(seg.b), p: seg.b });
                }
                if (sVEF2 == INSIDE) {
                    NLA.mapAdd(edgeMap, col2, { edgeT: col2.curve.pointLambda(seg.a), p: seg.a });
                    NLA.mapAdd(edgeMap, col2, { edgeT: col2.curve.pointLambda(seg.b), p: seg.b });
                }
            }
        }
        // what needs to be generated: new edges on face
        // points on edges where they are cut by faces so that sub edges will be generated for loops
        // points on ends of edges where the edge will be an edge in the new volume where it goes from A to B
        //         you don't want thos to be marked as "inside", otherwise invalid faces will be added
        // if a face cuts a corner, nothings needs to be done, as that alone does not limit what adjacent faces will be
        function handleEndPoint(a, b, useA, useB) {
            console.log("logging point", useA, a.toSource());
            console.log("logging point", useB, b.toSource());
            if (useA && !useB) {
                if (!a.colinear) {
                    NLA.mapAdd(edgeMap, a.edge, a);
                }
            }
            if (useB && !useA) {
                if (!b.colinear) {
                    NLA.mapAdd(edgeMap, b.edge, b);
                }
            }
            if (useA && useB) {
                // if a or b is colinear the correct points will already have been added to the edge by handleNewEdge
                // segment starts/ends on edge/edge intersection
                function foo(a, b, face, face2, thisPlane, face2Plane, thisBrep, face2Brep, first) {
                    if (!a.colinear && a.edgeT != a.edge.aT && a.edgeT != a.edge.bT) {
                        //						assert(a.edgeT != a.edge.aT && a.edgeT != a.edge.bT, "implement point intersections") todo
                        // ends on a, on colinear segment b bT != a.edge.bT &&
                        // b can be colinear, so edgeT == aT is possible
                        if (a.p.like(b.edge.a) || a.p.like(b.edge.b)) {
                            // face2brep corner on edge
                            var sVEC1 = splitsVolumeEnclosingCone(face2Brep, a.p, a.edge.aDir);
                            var sVEC2 = splitsVolumeEnclosingCone(face2Brep, a.p, a.edge.aDir.negated());
                            // if either of these return ALONG_EDGE_OR_PLANE, then the breps share a colinear edge
                            if ((sVEC1 != OUTSIDE) != (sVEC2 != OUTSIDE)) {
                                NLA.mapAdd(edgeMap, a.edge, a);
                            }
                        }
                        else {
                            var testVector = a.edge.aDir.rejectedFrom(b.edge.aDir);
                            assert(!testVector.isZero());
                            var sVEF1 = splitsVolumeEnclosingFaces(face2Brep, b.edge, testVector, thisPlane.normal);
                            var sVEF2 = splitsVolumeEnclosingFaces(face2Brep, b.edge, testVector.negated(), thisPlane.normal);
                            if (!(INSIDE == sVEF1 && INSIDE == sVEF2
                                || (OUTSIDE == sVEF1 || COPLANAR_OPPOSITE == sVEF1) && (OUTSIDE == sVEF2 || COPLANAR_OPPOSITE == sVEF2))) {
                                NLA.mapAdd(edgeMap, a.edge, a);
                            }
                        }
                    }
                }
                foo(a, b, face, face2, thisPlane, face2Plane, thisBrep, face2Brep, true);
                foo(b, a, face2, face, face2Plane, thisPlane, face2Brep, thisBrep, false);
            }
        }
        var in1 = false, in2 = false, colinear1 = false, colinear2 = false;
        var i = 0, j = 0, last, segments = [];
        var startP, startDir, startT;
        while (i < ps1.length || j < ps2.length) {
            assert(i <= ps1.length);
            assert(j <= ps2.length);
            var a = ps1[i], b = ps2[j];
            assert(a || b);
            if (j == ps2.length || i < ps1.length && NLA.lt(a.t, b.t)) {
                last = a;
                in1 = !in1;
                // ": colinear1" to remember the colinear segment, as segments are only handled once it ends (in1 false)
                colinear1 = in1 ? a.colinear && a.edge : colinear1;
                i++;
            }
            else if (i == ps1.length || NLA.gt(a.t, b.t)) {
                last = b;
                in2 = !in2;
                colinear2 = in2 ? b.colinear && b.edge : colinear2;
                j++;
            }
            else {
                last = a;
                in1 = !in1;
                in2 = !in2;
                //if (in1 == in2) {
                a.used = true;
                b.used = true;
                colinear1 = in1 ? a.colinear && a.edge : colinear1;
                colinear2 = in2 ? b.colinear && b.edge : colinear2;
                //}
                i++;
                j++;
            }
            console.log("as", a, b, in1, in2);
            if (startP && !(in1 && in2)) {
                // segment end
                var newEdge = new StraightEdge(intersectionLine, startP, last.p, startT, last.t, null, 'genseg' + globalId++);
                handleNewEdge(newEdge, colinear1, colinear2, in1, in2, a, b);
                startP = undefined;
                last.used = true;
                handleEndPoint(a, b, colinear1 || !!a.used, colinear2 || !!b.used);
            }
            else if (in1 && in2) {
                // new segment just started
                startP = last.p;
                startDir = last.insideDir;
                startT = last.t;
                last.used = true;
                handleEndPoint(a, b, colinear1 || !!a.used, colinear2 || !!b.used);
            }
        }
        var _loop_1 = function(edge, ps) {
            assert(ps.every(function (p) { return edge.isValidT(p.edgeT); }), edge.sce, ps.sce);
        };
        for (var _i = 0, _a = edgeMap.entries(); _i < _a.length; _i++) {
            var _b = _a[_i], edge = _b[0], ps = _b[1];
            _loop_1(edge, ps);
        }
    };
    PlaneFace.forVertices = function (planeSurface, vs) {
        var holeVs = [];
        for (var _i = 2; _i < arguments.length; _i++) {
            holeVs[_i - 2] = arguments[_i];
        }
        if (planeSurface instanceof P3) {
            planeSurface = new PlaneSurface(planeSurface);
        }
        assert(isCCW(vs, planeSurface.plane.normal), 'isCCW(vs, planeSurface.plane.normal)');
        var edges = vs.map(function (a, i, vs) {
            var b = vs[(i + 1) % vs.length];
            return StraightEdge.throughPoints(a, b);
        });
        holeVs.forEach(function (vs) { return assert(doubleSignedArea(vs, planeSurface.plane.normal) >= 0, 'doubleSignedArea(vs, planeSurface.plane.normal) >= 0'); });
        var holes = holeVs.map(function (hvs) { return hvs.map(function (a, i, vs) {
            var b = vs[(i + 1) % vs.length];
            return StraightEdge.throughPoints(a, b);
        }); });
        return new PlaneFace(planeSurface, edges, holes);
    };
    PlaneFace.prototype.pointsToInside = function (p, dir) {
        return this.containsPoint2(p.plus(dir.times(NLA_PRECISION * 8)));
    };
    return PlaneFace;
}(Face));
NLA.registerClass(PlaneFace);
var RotationFace = (function (_super) {
    __extends(RotationFace, _super);
    function RotationFace(rot, contour, holes, name) {
        _super.call(this, rot, contour, holes, name);
    }
    RotationFace.prototype.unrollLoop = function (edgeLoop) {
        var vs = [];
        var reverseFunc = this.surface.pointToParameterFunction();
        var verticesNo0s = edgeLoop.map(function (edge) { return edge.getVerticesNo0(); });
        var startEdgeIndex = verticesNo0s.findIndex(function (edgeVertices) { return !NLA.eq(reverseFunc(edgeVertices[0], Math.PI).x, Math.PI); });
        assert(-1 != startEdgeIndex);
        // console.log(startEdgeIndex)
        var hint = Math.PI;
        for (var i = 0; i < edgeLoop.length; i++) {
            var edgeIndex = (i + startEdgeIndex) % edgeLoop.length;
            for (var j = 0; j < verticesNo0s[edgeIndex].length; j++) {
                var p = verticesNo0s[edgeIndex][j];
                var localP = reverseFunc(p, hint);
                if (Math.abs(localP.x) < Math.PI - NLA_PRECISION) {
                    // update hint
                    hint = localP.x;
                }
                // console.log(hint, p.sce, localP.sce)
                vs.push(localP);
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
        return vs;
    };
    RotationFace.prototype.addToMesh = function (mesh) {
        var _this = this;
        var closed = false;
        var hSplit = 32, zSplit = 1;
        var ribs = [];
        var minZ = Infinity, maxZ = -Infinity;
        var cmp = function (a, b) { return a.value - b.value; };
        var f = this.surface.parametricFunction();
        var normalF = this.surface.parametricNormal();
        var vertexLoops = this.holes.concat([this.edges]).map(function (loop) { return _this.unrollLoop(loop); });
        vertexLoops.forEach(function (vertexLoop) {
            vertexLoop.forEach(function (_a) {
                var d = _a.x, z = _a.y;
                var index0 = ribs.binaryIndexOf(d, function (a, b) { return NLA.snapTo(a.value - b, 0); });
                if (index0 < 0) {
                    ribs.splice(-index0 - 1, 0, { value: d, left: [], right: [] });
                }
                minZ = min(minZ, z);
                maxZ = max(maxZ, z);
            });
        });
        var correction = 1;
        vertexLoops.forEach(function (vertexLoop) {
            vertexLoop.forEach(function (v0, i, vs) {
                var v1 = vs[(i + 1) % vs.length], dDiff = v1.x - v0.x;
                //console.log(v0.sce, v1.sce)
                if (NLA.eq0(dDiff)) {
                    return;
                }
                if (dDiff < 0) {
                    _a = [v1, v0], v0 = _a[0], v1 = _a[1];
                    dDiff = -dDiff;
                }
                var index0 = ribs.binaryIndexOf(v0.x, function (a, b) { return NLA.snapTo(a.value - b, 0); });
                var index1 = ribs.binaryIndexOf(v1.x, function (a, b) { return NLA.snapTo(a.value - b, 0); });
                ribs[index0].right.binaryInsert(v0.y);
                for (var j = (index0 + correction) % ribs.length; j != index1; j = (j + correction) % ribs.length) {
                    var x = ribs[j].value;
                    var part = (x - v0.x) / dDiff;
                    var interpolated = v1.y * part + v0.y * (1 - part);
                    ribs[j].left.binaryInsert(interpolated);
                    ribs[j].right.binaryInsert(interpolated);
                }
                ribs[index1].left.binaryInsert(v1.y);
                var _a;
                // console.log(ribs.map(r=>r.toSource()).join('\n'))
            });
        });
        var vertices = [], triangles = [], normals = [];
        for (var i = 0; i < ribs.length; i++) {
            var ribLeft = ribs[i], ribRight = ribs[(i + 1) % ribs.length];
            assert(ribLeft.right.length == ribRight.left.length);
            for (var j = 0; j < ribLeft.right.length; j++) {
                vertices.push(f(ribLeft.value, ribLeft.right[j]), f(ribRight.value, ribRight.left[j]));
                normals.push(normalF(ribLeft.value, ribLeft.right[j]), normalF(ribRight.value, ribRight.left[j]));
            }
        }
        //console.log(ribs.map(r=>r.toSource()).join('\n'))
        var vss = vertices.length, detailVerticesStart = vss;
        var zInterval = maxZ - minZ, zStep = zInterval / zSplit;
        var detailZs = NLA.arrayFromFunction(zSplit - 1, function (i) { return minZ + (1 + i) * zStep; });
        for (var i = 0; i < ribs.length; i++) {
            var d = ribs[i].value;
            for (var j = 0; j < detailZs.length; j++) {
                vertices.push(f(d, detailZs[j]));
                normals.push(normalF(d, detailZs[j]));
            }
        }
        // console.log('detailVerticesStart', detailVerticesStart, 'vl', vertices.length, vertices.length - detailVerticesStart, ribs.length)
        // finally, fill in the ribs
        var vsStart = 0;
        var flipped2 = true;
        //for (var i = 0; i < 1; i++) {
        var end = closed ? ribs.length : ribs.length - 1;
        for (var i = 0; i < end; i++) {
            var ipp = (i + 1) % ribs.length;
            var inside = false, colPos = 0, ribLeft = ribs[i], ribRight = ribs[(i + 1) % ribs.length];
            for (var j = 0; j < detailZs.length + 1; j++) {
                var detailZ = detailZs[j] || 100000;
                if (!inside) {
                    if (ribLeft.right[colPos] < detailZ && ribRight.left[colPos] < detailZ) {
                        if (ribLeft.right[colPos + 1] < detailZ || ribRight.left[colPos + 1] < detailZ) {
                            pushQuad(triangles, flipped2, vsStart + colPos * 2, vsStart + (colPos + 1) * 2, vsStart + colPos * 2 + 1, vsStart + (colPos + 1) * 2 + 1);
                            colPos += 2;
                            if (ribLeft.right[colPos] < detailZ || ribRight.left[colPos] < detailZ) {
                                j--;
                            }
                        }
                        else {
                            pushQuad(triangles, flipped2, vsStart + colPos * 2, vsStart + colPos * 2 + 1, detailVerticesStart + i * detailZs.length + j, detailVerticesStart + ipp * detailZs.length + j);
                            inside = true;
                            colPos++;
                        }
                    }
                }
                else {
                    if (ribLeft.right[colPos] < detailZ || ribRight.left[colPos] < detailZ) {
                        pushQuad(triangles, flipped2, detailVerticesStart + i * detailZs.length + j - 1, detailVerticesStart + ipp * detailZs.length + j - 1, vsStart + colPos * 2, vsStart + colPos * 2 + 1);
                        inside = false;
                        colPos++;
                        if (ribLeft.right[colPos] < detailZ || ribRight.left[colPos] < detailZ) {
                            j--;
                        }
                    }
                    else {
                        pushQuad(triangles, flipped2, detailVerticesStart + i * detailZs.length + j, detailVerticesStart + i * detailZs.length + j - 1, detailVerticesStart + ipp * detailZs.length + j, detailVerticesStart + ipp * detailZs.length + j - 1);
                    }
                }
            }
            vsStart += ribLeft.right.length * 2;
        }
        //console.log('trinagle', triangles.max(), vertices.length, triangles.length, triangles.toSource(), triangles.map(i => vertices[i].$).toSource() )
        triangles = triangles.map(function (index) { return index + mesh.vertices.length; });
        //assert(normals.every(n => n.hasLength(1)), normals.find(n => !n.hasLength(1)).length() +" "+normals.findIndex(n => !n.hasLength(1)))
        Array.prototype.push.apply(mesh.vertices, vertices);
        Array.prototype.push.apply(mesh.triangles, triangles);
        Array.prototype.push.apply(mesh.normals, normals);
        //this.addEdgeLines(mesh)
    };
    RotationFace.prototype.intersectFace = function (face2, thisBrep, face2Brep, faceMap, edgeMap, likeSurfaceFaces) {
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
            console.log("handling seg", col1 && col1.toSource(), col2 && col2.toSource(), seg.toSource());
            if (!col1 && !col2) {
                console.log("adding");
                NLA.mapAdd(faceMap, face, seg);
                NLA.mapAdd(faceMap, face2, seg.flipped());
            }
            var segASurfaceNormal = surface.normalAt(seg.a);
            var segASurface2Normal = surface2.normalAt(seg.a);
            if (!col1 && col2) {
                console.log('!col1 && col2');
                // TODO: replace with curve.tangent DOT edge.tangent
                if (col2.aDir.cross(surface2.normalAt(col2.a)).dot(surface.normalAt(col2.a)) > 0) {
                    // if (seg.aDir.cross(segASurface2Normal).dot(segASurfaceNormal) > 0) {
                    // NB: a new edge is inserted even though it may be the same as an old one
                    // however it indicates that it intersects the other volume here, i.e. the old edge cannot
                    // be counted as "inside" for purposes of reconstitution
                    //NLA.mapAdd(faceMap, face2, seg.flipped()) todo
                    NLA.mapAdd(edgeMap, col2, { edgeT: col2.curve.pointLambda(seg.a), p: seg.a });
                    NLA.mapAdd(edgeMap, col2, { edgeT: col2.curve.pointLambda(seg.b), p: seg.b });
                    var testVector = surface2.normalAt(col2.a).negated().rejectedFrom(surface.normalAt(col2.a));
                    console.log("testVector", testVector.sce, face2Brep.sce, splitsVolumeEnclosingFaces(face2Brep, col2, testVector, surface.normal) == INSIDE);
                    if (splitsVolumeEnclosingFaces(face2Brep, col2, testVector, surface.normalAt(col2.a)) == INSIDE) {
                        NLA.mapAdd(faceMap, face, seg);
                    }
                }
            }
            if (!col2 && col1) {
                if (col1.aDir.cross(surface.normalAt(col1.a)).dot(surface2.normalAt(col1.a)) > 0) {
                    //NLA.mapAdd(faceMap, face, seg) todo
                    NLA.mapAdd(edgeMap, col1, { edgeT: col1.curve.pointLambda(seg.a), p: seg.a });
                    NLA.mapAdd(edgeMap, col1, { edgeT: col1.curve.pointLambda(seg.b), p: seg.b });
                    var testVector = surface.normalAt(col1.a).negated().rejectedFrom(surface2.normalAt(col1.a));
                    if (splitsVolumeEnclosingFaces(thisBrep, col1, testVector, surface2.normalAt(col1.a)) == INSIDE) {
                        NLA.mapAdd(faceMap, face2, seg.flipped());
                    }
                }
            }
            if (col1 && col2) {
                // faces are touching edge-edge
                var testVector1 = col1.aDir.cross(surface.normalAt(col1.a)).negated();
                var testVector2 = col2.aDir.cross(surface2.normalAt(col2.a)).negated();
                console.log("testVector", testVector1.sce, face2Brep.sce, splitsVolumeEnclosingFaces(face2Brep, col2, testVector1, surface.normal) == INSIDE);
                if (splitsVolumeEnclosingFaces(face2Brep, col2, testVector1, surface.normalAt(col2.a)) == INSIDE
                    && splitsVolumeEnclosingFaces(thisBrep, col1, testVector2, surface2.normalAt(col1.a)) == INSIDE) {
                    //NLA.mapAdd(faceMap, face2, seg.flipped()) todo
                    NLA.mapAdd(edgeMap, col1, { edgeT: col1.curve.pointLambda(seg.a), p: seg.a });
                    NLA.mapAdd(edgeMap, col1, { edgeT: col1.curve.pointLambda(seg.b), p: seg.b });
                    NLA.mapAdd(edgeMap, col2, { edgeT: col2.curve.pointLambda(seg.a), p: seg.a });
                    NLA.mapAdd(edgeMap, col2, { edgeT: col2.curve.pointLambda(seg.b), p: seg.b });
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
            console.log("logging point", useA, a && a.toSource());
            console.log("logging point", useB, b && b.toSource());
            if (useA && !useB) {
                if (!a.colinear) {
                    NLA.mapAdd(edgeMap, a.edge, a);
                }
            }
            if (useB && !useA) {
                if (!b.colinear) {
                    NLA.mapAdd(edgeMap, b.edge, b);
                }
            }
            if (useA && useB) {
                // segment start/ends on edge/edge intersection
                if (!a.colinear && !NLA.eq(a.edgeT, a.edge.aT) && !NLA.eq(a.edgeT, a.edge.bT)) {
                    //						assert(a.edgeT != a.edge.aT && a.edgeT != a.edge.bT, "implement point intersections")
                    // ends on a, on colinear segment b bT != a.edge.bT &&
                    var testVector = a.edge.aDir.rejectedFrom(b.edge.aDir);
                    var sVEF1 = splitsVolumeEnclosingFaces(face2Brep, b.edge, testVector, surface.normalAt(b.edge.a));
                    var sVEF2 = splitsVolumeEnclosingFaces(face2Brep, b.edge, testVector.negated(), surface.normalAt(b.edge.a));
                    if (!(INSIDE == sVEF1 && INSIDE == sVEF2
                        || (OUTSIDE == sVEF1 || COPLANAR_OPPOSITE == sVEF1) && (OUTSIDE == sVEF2 || COPLANAR_OPPOSITE == sVEF2))) {
                        NLA.mapAdd(edgeMap, a.edge, a);
                    }
                }
                if (!b.colinear && !NLA.eq(b.edgeT, b.edge.aT) && !NLA.eq(b.edgeT, b.edge.bT)) {
                    //						assert(b.edgeT != b.edge.aT && b.edgeT != b.edge.bT, "implement point intersections")
                    // ends on b, on colinear segment a
                    var testVector = b.edge.aDir.rejectedFrom(a.edge.aDir);
                    var sVEF1 = splitsVolumeEnclosingFaces(thisBrep, a.edge, testVector, surface2.normal);
                    var sVEF2 = splitsVolumeEnclosingFaces(thisBrep, a.edge, testVector.negated(), surface2.normal);
                    if (!(INSIDE == sVEF1 && INSIDE == sVEF2
                        || (OUTSIDE == sVEF1 || COPLANAR_OPPOSITE == sVEF1) && (OUTSIDE == sVEF2 || COPLANAR_OPPOSITE == sVEF2))) {
                        NLA.mapAdd(edgeMap, b.edge, b);
                    }
                }
            }
        }
        assertInst(Face, face2);
        var face = this;
        var surface = face.surface, surface2 = face2.surface;
        if (surface.isCoplanarTo(surface2)) {
            addLikeSurfaceFaces(likeSurfaceFaces, face, face2);
            return;
        }
        var isCurves = surface.isCurvesWithSurface(surface2);
        // get intersections of newCurve with other edges of face and face2
        var pss1 = faceEdgeISPsWithSurface(thisBrep, face, isCurves, face2.surface);
        var pss2 = faceEdgeISPsWithSurface(face2Brep, face2, isCurves, face.surface);
        console.log('pss1\n', pss1.map(function (m) { return m.toSource(); }).join('\n'), '\npss2\n', pss2.map(function (m) { return m.toSource(); }).join('\n'));
        console.log('' + surface + surface2);
        isCurves.forEach(function (isCurve, isCurveIndex) {
            var ps1 = pss1[isCurveIndex], ps2 = pss2[isCurveIndex];
            // for non-endless curves, e.g. ellipses, the intersections of the faces can be non-zero, even if one of
            // the faces doesn't register any points on the curve. For example, if a cylinder is cut entirely by a
            // plane face (all its edges around the cylinder), then the face will contain the entire curve and
            // "ps" for the plane face will be empty
            var curvePoint = isCurve.at(0);
            // TODO: behavior when curves touch face?
            var faceContainsCurvePoint = face.containsPoint(curvePoint);
            var face2ContainsCurvePoint = face2.containsPoint(curvePoint);
            if (ps1.length == 0 && !faceContainsCurvePoint || ps2.length == 0 && !face2ContainsCurvePoint) {
                return;
            }
            // !! start in does depend on insidedir... TODO
            assertf(function () { return (0 == ps1.length) || !NLA.eq0(ps1[0].insideDir.dot(isCurve.tangentAt(ps1[0].t))); }, function () { return ps1[0].insideDir.dot(isCurve.tangentAt(ps1[0].t)); });
            assertf(function () { return (0 == ps2.length) || !NLA.eq0(ps2[0].insideDir.dot(isCurve.tangentAt(ps2[0].t))); }, function () { return ps2[0].insideDir.dot(isCurve.tangentAt(ps2[0].t)); });
            // if ps is empty, in must be true or we wouldn't have gotten this far
            var in1 = (0 == ps1.length) || ps1[0].insideDir.dot(isCurve.tangentAt(ps1[0].t)) < 0;
            var in2 = (0 == ps2.length) || ps2[0].insideDir.dot(isCurve.tangentAt(ps2[0].t)) < 0;
            var colinear1 = false, colinear2 = false;
            var i = 0, j = 0, last, segments = [];
            var startP, startDir, startT;
            while (i < ps1.length || j < ps2.length) {
                var a = ps1[i], b = ps2[j];
                if (j >= ps2.length || i < ps1.length && NLA.lt(a.t, b.t)) {
                    last = a;
                    in1 = !in1;
                    // ": colinear1" to remember the colinear segment, as segments are only handled once it ends (in1 false)
                    colinear1 = in1 ? a.colinear && a.edge : colinear1;
                    i++;
                }
                else if (i >= ps1.length || NLA.gt(a.t, b.t)) {
                    last = b;
                    in2 = !in2;
                    colinear2 = in2 ? b.colinear && b.edge : colinear2;
                    j++;
                }
                else {
                    last = a;
                    in1 = !in1;
                    in2 = !in2;
                    //if (in1 == in2) {
                    a.used = true;
                    b.used = true;
                    colinear1 = in1 ? a.colinear && a.edge : colinear1;
                    colinear2 = in2 ? b.colinear && b.edge : colinear2;
                    //}
                    i++;
                    j++;
                }
                console.log("as", a, b, in1, in2);
                if (startP && !(in1 && in2)) {
                    // segment end
                    handleGeneratedSegment(Edge.create(isCurve, startP, last.p, startT, last.t, null, isCurve.tangentAt(startT), isCurve.tangentAt(last.t)), colinear1, colinear2, in1, in2, a, b);
                    startP = undefined;
                    last.used = true;
                    handlePoint(a, b, colinear1 || a && !!a.used, colinear2 || b && !!b.used);
                }
                else if (in1 && in2) {
                    // new segment just started
                    startP = last.p;
                    startDir = last.insideDir;
                    startT = last.t;
                    last.used = true;
                    handlePoint(a, b, colinear1 || a && !!a.used, colinear2 || b && !!b.used);
                }
            }
        });
    };
    return RotationFace;
}(Face));
NLA.registerClass(RotationFace);
function addLikeSurfaceFaces(likeSurfaceFaces, face1, face2) {
    // There cannot be two subgroups which will later be connected, as the "graph" of like surface faces is fully connected
    for (var i = 0; i < likeSurfaceFaces.length; i++) {
        var faceGroup = likeSurfaceFaces[i];
        var foundFace1 = false, foundFace2 = false;
        for (var j = 0; j < faceGroup.length; j++) {
            var face = faceGroup[j];
            if (face == face1) {
                foundFace1 = true;
            }
            if (face == face2) {
                foundFace2 = true;
            }
        }
        if (foundFace1 != foundFace2) {
            faceGroup.push(foundFace1 ? face2 : face1);
            return;
        }
        else if (foundFace1) {
            // found both
            return;
        }
    }
    // nothing found, add a new group
    likeSurfaceFaces.push([face1, face2]);
}
var B2 = (function (_super) {
    __extends(B2, _super);
    function B2(faces, infiniteVolume, generator, vertexNames) {
        _super.call(this);
        this.faces = faces;
        assertInst.apply(undefined, [Face].concat(faces));
        this.infiniteVolume = !!infiniteVolume;
        this.generator = generator;
        this.vertexNames = vertexNames;
    }
    B2.prototype.calculateVolume = function () {
        return this.faces.map(function (face) { return face.zDirVolume(); }).sum();
    };
    B2.prototype.toMesh = function () {
        var mesh = new GL.Mesh({ triangles: true, normals: true, lines: true });
        mesh.faceIndexes = new Map();
        this.faces.forEach(function (face, i) {
            var triangleStart = mesh.triangles.length;
            face.addToMesh(mesh);
            mesh.faceIndexes.set(face, { start: triangleStart, count: mesh.triangles.length - triangleStart });
        });
        mesh.compile();
        return mesh;
    };
    B2.prototype.minus = function (brep2) {
        return this.intersection(brep2.flipped(), true, true);
    };
    B2.prototype.plus = function (brep2) {
        var result = this.flipped().intersection(brep2.flipped(), true, true).flipped();
        result.generator = this.generator + ".plus(" + brep2.generator + ")";
        return result;
    };
    B2.prototype.and = function (other) {
        return this.intersection(other, true, true, this.generator && other.generator && this.generator + '.and(' + other.generator + ')');
    };
    B2.prototype.xor = function (other) {
        return new B2(this.minus(other).faces.concat(other.minus(this).faces), this.infiniteVolume != other.infiniteVolume, this.generator && other.generator && this.generator + '.xor(' + other.generator + ')');
    };
    B2.prototype.equals = function (brep) {
        return this.faces.length == brep.faces.length &&
            this.faces.every(function (face) { return brep.faces.some(function (face2) { return face.equals(face2); }); });
    };
    B2.prototype.toString = function () {
        return "new B2([\n" + this.faces.join(',\n').replace(/^/gm, '\t') + "])";
    };
    B2.prototype.toSource = function () {
        return this.generator ||
            'new B2([\n' + this.faces.map(function (face) { return face.toSource(); }).join(',\n').replace(/^/gm, '\t') + '])';
    };
    B2.loop1ContainsLoop2 = function (loop1, loop2, surface) {
        for (var _i = 0, loop2_1 = loop2; _i < loop2_1.length; _i++) {
            var edge = loop2_1[_i];
            var loop1ContainsPoint = surface.loopContainsPoint(loop1, edge.a);
            if (PointVsFace.ON_EDGE != loop1ContainsPoint)
                return PointVsFace.INSIDE == loop1ContainsPoint;
        }
        for (var _a = 0, loop2_2 = loop2; _a < loop2_2.length; _a++) {
            var edge = loop2_2[_a];
            var edgePoint = edge.curve.at(edge.aT * 0.2 + edge.bT * 0.8);
            var loop1ContainsPoint = surface.loopContainsPoint(loop1, edgePoint);
            if (PointVsFace.ON_EDGE != loop1ContainsPoint)
                return PointVsFace.INSIDE == loop1ContainsPoint;
        }
        throw new Error();
    };
    B2.assembleFacesFromLoops = function (loops, surface, faceConstructor) {
        function placeRecursively(newLoopInfo, loopInfos) {
            if (loopInfos.length == 0) {
                loopInfos.push(newLoopInfo);
            }
            else {
                var subLoopInfo = loopInfos.find(function (loopInfo) { return B2.loop1ContainsLoop2(loopInfo.loop, newLoopInfo.loop, surface); });
                console.log("here", subLoopInfo);
                if (subLoopInfo) {
                    placeRecursively(newLoopInfo, subLoopInfo.subloops);
                }
                else {
                    // newLoopInfo isnt contained by any other subLoopInfo
                    for (var i = loopInfos.length; --i >= 0;) {
                        var subLoopInfo_1 = loopInfos[i];
                        //console.log("cheving subLoopInfo", surface.loopContainsPoint(newLoopInfo.edges, subLoopInfo.edges[0].a))
                        if (B2.loop1ContainsLoop2(newLoopInfo.loop, subLoopInfo_1.loop, surface)) {
                            newLoopInfo.subloops.push(subLoopInfo_1);
                            loopInfos.splice(i, 1); // remove it
                        }
                    }
                    loopInfos.push(newLoopInfo);
                }
            }
        }
        function newFacesRecursive(loopInfo) {
            assert(loopInfo.ccw);
            if (loopInfo.subloops.every(function (sl) { return !sl.ccw; })) {
                var newFace = new faceConstructor(surface, loopInfo.loop, loopInfo.subloops.map(function (sl) { return sl.loop; }));
                newFaces.push(newFace);
                loopInfo.subloops.forEach(function (sl) { return sl.subloops.forEach(function (slsl) { return newFacesRecursive(slsl); }); });
            }
            else {
                loopInfo.subloops.forEach(function (sl) { return sl.ccw && newFacesRecursive(sl); });
            }
        }
        var newFaces = [];
        console.log(loops.map(function (loop) { return loop.join('\n'); }).join('\n\n'));
        var topLevelLoops = [];
        loops.forEach(function (loop) { return placeRecursively({ loop: loop, ccw: surface.edgeLoopCCW(loop), subloops: [] }, topLevelLoops); });
        topLevelLoops.forEach(function (tll) { return newFacesRecursive(tll); });
        return newFaces;
    };
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
    B2.reconstituteFaces = function (oldFaces, edgeSubEdges, faceMap, newFaces) {
        var oldFaceStatuses = new Map();
        // reconstitute faces
        var insideEdges = [];
        oldFaces.forEach(function (face, faceIndex) {
            console.log('reconstituting face', face.toString());
            var usableOldEdges = face.getAllEdges().filter(function (edge) { return !edgeSubEdges.get(edge); });
            var subEdges = face.getAllEdges().mapFilter(function (edge) { return edgeSubEdges.get(edge); }).concatenated();
            var els = face.edges.map(function (edge) { return edgeSubEdges.get(edge) && edgeSubEdges.get(edge).join('\n'); }).map(function (s) { return '\n\n' + s; }).join();
            console.log('edgeSubEdges', els);
            var newEdges = faceMap.get(face) || [];
            if (newEdges.length || subEdges.length) {
                oldFaceStatuses.set(face, 'partial');
                console.log('newEdges\n', newEdges.map(function (e) { return e.toString(); }).join('\n'));
                var loops = [];
                var edgeCond = face instanceof PlaneFace
                    ? (function (edge) { return edge.a.like(currentEdge_1.b); })
                    : (function (edge) { return (edge.curve == currentEdge_1.curve // TODO: ??
                        ? NLA.eq(edge.aT, currentEdge_1.bT)
                        : edge.a.like(currentEdge_1.b)); });
                // new edges are definitely part of a resulting loop
                // old edges (both contour and holes) can either be part of a new loop, in which case they will already
                // have been visited when starting a loop search with a new edge, OR they can be stranded, OR they can
                // remain in their old loop
                function getNextStart() {
                    return newEdges.find(function (edge) { return !edge.visited; })
                        || subEdges.find(function (edge) { return !edge.visited; })
                        || usableOldEdges.find(function (edge) { return !edge.visited; });
                }
                usableOldEdges.forEach(function (edge) { return edge.visited = false; });
                // search for a loop:
                var currentEdge_1;
                while (currentEdge_1 = getNextStart()) {
                    var cancelLoop = false;
                    var startEdge = currentEdge_1, edges = [], i = 0;
                    // wether only new edges are used (can include looseSegments)
                    var _loop_2 = function() {
                        currentEdge_1.visited = true;
                        console.log('currentEdge', currentEdge_1.b.sce, currentEdge_1.toSource());
                        edges.push(currentEdge_1);
                        // find next edge
                        var possibleOldEdges = usableOldEdges.filter(function (edge) { return currentEdge_1.b.like(edge.a); });
                        var possibleSubEdges = subEdges.filter(function (edge) { return currentEdge_1.b.like(edge.a); });
                        var possibleNewEdges = newEdges.filter(function (edge) { return currentEdge_1.b.like(edge.a); });
                        var possibleEdges = possibleOldEdges.concat(possibleSubEdges, possibleNewEdges);
                        assert(0 < possibleEdges.length);
                        var faceNormalAtCurrentB = face.surface.normalAt(currentEdge_1.b);
                        var nextEdgeIndex = possibleEdges.indexWithMax(function (edge, index) { return currentEdge_1.bDir.angleRelativeNormal(edge.aDir, faceNormalAtCurrentB); });
                        currentEdge_1 = possibleEdges[nextEdgeIndex];
                        if (currentEdge_1.visited) {
                            // this start edge doesn't lead to a valid loop
                            console.log("breaking");
                            return "break";
                        }
                        assert(currentEdge_1);
                        assert(currentEdge_1 != startEdge);
                    };
                    do {
                        var state_2 = _loop_2();
                        if (state_2 === "break") break;
                    } while (++i < 200);
                    if (200 == i) {
                        assert(false, "too many");
                    }
                    if (currentEdge_1 == startEdge) {
                        console.log('finished loop');
                        loops.push(edges);
                    }
                }
                var faceNewFaces = B2.assembleFacesFromLoops(loops, face.surface, face.constructor);
                newFaces.pushAll(faceNewFaces);
                var faceNewFacesEdges_1 = faceNewFaces.map(function (face) { return face.getAllEdges(); }).concatenated();
                insideEdges.pushAll(usableOldEdges.filter(function (edge) { return faceNewFacesEdges_1.includes(edge); }));
            }
        });
        console.log("INSIDE EDGES", insideEdges, oldFaceStatuses);
        while (insideEdges.length != 0) {
            var insideEdge = insideEdges.pop();
            var adjacentFaces = facesWithEdge(insideEdge, oldFaces);
            adjacentFaces.forEach(function (info) {
                if (!oldFaceStatuses.has(info.face)) {
                    oldFaceStatuses.set(info.face, 'inside');
                    insideEdges.push.apply(insideEdges, info.face.edges);
                }
            });
        }
        newFaces.pushAll(oldFaces.filter(function (face) { return oldFaceStatuses.get(face) == 'inside'; }));
    };
    B2.prototype.reconstituteCoplanarFaces = function (likeSurfacePlanes, edgeLooseSegments, faceMap, newFaces) {
        likeSurfacePlanes.forEach(function (faceGroup) {
            // calculate total contours
            var surface = faceGroup[0].surface, bag = [];
            faceGroup.forEach(function (face) {
                Array.prototype.push.apply(bag, faceMap(face));
                face.getAllEdges().forEach(function (edge) {
                    var edgeSubSegments;
                    if (edgeSubSegments = edgeLooseSegments.get(edge)) {
                        Array.prototype.push.apply(bag, edgeSubSegments);
                    }
                    else {
                        bag.push(edge);
                    }
                });
            });
            var currentEdge, loops = [];
            while (currentEdge = bag.find(function (edge) { return !edge.visited; })) {
                var path = [];
                do {
                    currentEdge.visited = true;
                    path.push(currentEdge);
                    var possibleNextEdges = bag.filter(function (edge) { return currentEdge.b.like(edge.a); });
                    // lowest angle, i.e. the right-most next edge
                    var nextEdgeIndex = possibleNextEdges.indexWithMax(function (edge, index) { return -currentEdge.bDir.angleRelativeNormal(edge.aDir, surface.normalAt(currentEdge.b)); });
                    currentEdge = possibleNextEdges[nextEdgeIndex];
                } while (!currentEdge.visited);
                var startIndex = path.find(currentEdge);
                if (-1 != startIndex) {
                    loops.push(path.slice(startIndex));
                }
            }
        });
    };
    //getLooseEdgeSegmentsInsideDirs(edgeMap: Map<Edge, IntersectionPointInfo[]>): Map<Edge, Edge> {
    //	var edgeLooseSegments = new Map()
    //	console.log("edgeMap", edgeMap)
    //	edgeMap.forEach((pointInfos, baseEdge) => {
    //		// TODO: make sure it works with loop
    //		// TODO: dont generate unnecessarry segments
    //		var looseSegments = []
    //		if (!baseEdge.reversed) {
    //			pointInfos.sort((a, b) => NLA.snapTo(a.edgeT - b.edgeT, 0) || (b.insideDir.dot(baseEdge.curve.dir1) - a.insideDir.dot(baseEdge.curve.dir1)))
    //		} else {
    //			pointInfos.sort((b, a) => NLA.snapTo(a.edgeT - b.edgeT, 0) || (b.insideDir.dot(baseEdge.curve.dir1) - a.insideDir.dot(baseEdge.curve.dir1)))
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
    B2.prototype.getLooseEdgeSegments = function (edgePointInfoss) {
        var result = new Map();
        console.log("edgePointInfoss", edgePointInfoss);
        edgePointInfoss.forEach(function (pointInfos, baseEdge) {
            // TODO: dont generate unnecessarry segments
            // if there are no point info, the original edge will be kept, so we should return nothing
            // otherwise, something will be returned, even if it a new edge identical to the base edge
            if (0 == pointInfos.length)
                return [];
            var subEdges = [];
            if (!baseEdge.reversed) {
                pointInfos.sort(function (a, b) { return a.edgeT - b.edgeT; });
            }
            else {
                pointInfos.sort(function (b, a) { return a.edgeT - b.edgeT; });
            }
            var startP = baseEdge.a, startDir = baseEdge.aDir, startT = baseEdge.aT, startInfo;
            for (var i = 0; i < pointInfos.length; i++) {
                var info = pointInfos[i];
                // ignore start, end and repeating points
                if (NLA.eq(info.edgeT, baseEdge.bT) || NLA.eq(info.edgeT, baseEdge.aT) || NLA.eq(info.edgeT, startT)) {
                    continue;
                }
                var pDir = baseEdge.tangentAt(info.edgeT);
                subEdges.push(Edge.create(baseEdge.curve, startP, info.p, startT, info.edgeT, null, startDir, pDir, 'looseSegment' + globalId++));
                startP = info.p;
                startT = info.edgeT;
                startInfo = info;
                startDir = pDir;
            }
            subEdges.push(Edge.create(baseEdge.curve, startP, baseEdge.b, startT, baseEdge.bT, null, startDir, baseEdge.bDir, 'looseSegment' + globalId++));
            result.set(baseEdge, subEdges);
        });
        return result;
    };
    B2.prototype.getIntersectionEdges = function (brep2) {
        var _this = this;
        var faceMap = new Map(), edgeMap = new Map();
        var likeSurfaceFaces = [];
        this.faces.forEach(function (face) {
            //console.log('face', face.toString())
            brep2.faces.forEach(function (face2) {
                //console.log('face2', face2.toString())
                face.intersectFace(face2, _this, brep2, faceMap, edgeMap, likeSurfaceFaces);
            });
        });
        return Array.from(faceMap.values()).concatenated();
    };
    B2.prototype.intersection = function (other, buildThis, buildOther, buildCoplanar) {
        var _this = this;
        var faceMap = new Map(), edgeMap = new Map();
        var likeSurfaceFaces = [];
        this.faces.forEach(function (face) {
            //console.log('face', face.toString())
            other.faces.forEach(function (face2) {
                //console.log('face2', face2.toString())
                face.intersectFace(face2, _this, other, faceMap, edgeMap, likeSurfaceFaces);
            });
        });
        var newFaces = [];
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
        buildThis && B2.reconstituteFaces(this.faces, edgeLooseSegments, faceMap, newFaces, other.infiniteVolume);
        buildOther && B2.reconstituteFaces(other.faces, edgeLooseSegments, faceMap, newFaces, this.infiniteVolume);
        //buildCoplanar && this.reconstituteCoplanarFaces(likeSurfaceFaces, edgeLooseSegments, faceMap, newFaces, this.infiniteVolume, other.infiniteVolume)
        return new B2(newFaces, this.infiniteVolume && other.infiniteVolume);
    };
    B2.prototype.transform = function (m4, desc) {
        desc = desc || '';
        var vertexNames;
        if (this.vertexNames) {
            vertexNames = new Map();
            this.vertexNames.forEach(function (name, vertex) { return vertexNames.set(m4.transformPoint(vertex), name + desc); });
        }
        return new B2(this.faces.map(function (f) { return f.transform(m4); }), this.infiniteVolume, this.generator && this.generator + desc, vertexNames);
    };
    B2.prototype.flipped = function () {
        return new B2(this.faces.map(function (f) { return f.flipped(); }), !this.infiniteVolume, this.generator && this.generator + '.flipped()', this.vertexNames);
    };
    B2.EMPTY = new B2([], false, 'B2.EMPTY', new Map());
    B2.R3 = new B2([], true, 'B2.EMPTY', new Map());
    return B2;
}(Transformable));
function facesWithEdge(edge, faces) {
    return faces.mapFilter(function (face) {
        var matchingEdge = face.getAllEdges().find(function (e) { return e.isCoEdge(edge); });
        if (matchingEdge) {
            return { face: face, reversed: !edge.a.like(matchingEdge.a), angle: NaN, normalAtEdgeA: null, edge: matchingEdge };
        }
    });
}
function planeFaceEdgeISPsWithPlane(brep, face, isLine, plane2) {
    // todo remove param brep
    assert(face.surface.plane.containsLine(isLine));
    assert(plane2.containsLine(isLine));
    var plane = face.surface.plane;
    var ps = [];
    var loops = [face.edges].concat(face.holes);
    loops.forEach(function (loop) {
        var colinearEdges = loop.map(function (edge) { return edge.colinearToLine(isLine); });
        var isLineOut = isLine.dir1.cross(plane.normal);
        loop.forEach(function (edge, edgeIndex, edges) {
            var nextEdgeIndex = (edgeIndex + 1) % edges.length, nextEdge = edges[nextEdgeIndex];
            //console.log(edge.toSource()) {p:V3(2, -2.102, 0),
            if (colinearEdges[edgeIndex]) {
                // edge colinear to intersection line
                var prevEdgeIndex = (edgeIndex - 1 + edges.length) % edges.length, prevEdge = edges[prevEdgeIndex];
                // if colinear, edge.curve must be a line, so colinearOut is constant
                var colinearOut = edge.aDir.cross(plane.normal);
                var curveAT = isLine.pointLambda(edge.a), curveBT = isLine.pointLambda(edge.b);
                // close previous interval if necessary
                if (!colinearEdges[prevEdgeIndex] && dotCurve(colinearOut, prevEdge.bDir, prevEdge.bDDT) < 0) {
                    ps.push({ p: prevEdge.b, insideDir: edge.aDir.negated(), t: curveAT, edge: prevEdge, edgeT: prevEdge.bT,
                        colinear: false });
                }
                // add interval for colinear segment
                ps.push({ p: edge.a, insideDir: edge.aDir, t: curveAT, edge: edge, edgeT: edge.aT,
                    colinear: true }, { p: edge.b, insideDir: edge.bDir.negated(), t: curveBT, edge: edge, edgeT: edge.bT,
                    colinear: true });
                // open next interval if necessary
                if (!colinearEdges[nextEdgeIndex] && dotCurve(colinearOut, nextEdge.aDir, nextEdge.aDDT) > 0) {
                    ps.push({ p: nextEdge.a, insideDir: edge.bDir, t: curveBT, edge: nextEdge, edgeT: nextEdge.aT, colinear: false });
                }
            }
            else {
                // not necessarily a straight edge, so multiple intersections are possible
                var edgeTs = edge.edgeISTsWithPlane(plane2);
                assert(edgeTs.every(function (t) { return plane2.containsPoint(edge.curve.at(t)); }), edgeTs);
                for (var k = 0; k < edgeTs.length; k++) {
                    var edgeT = edgeTs[k];
                    if (edgeT == edge.bT) {
                        // endpoint lies on intersection line
                        if (!colinearEdges[nextEdgeIndex]
                            && edgesDifferentSidesOfDir(isLineOut, edge, nextEdge)) {
                            // next segment is not colinear and ends on different side
                            ps.push({ p: edge.b, insideDir: plane2.normal.negated(), t: isLine.pointLambda(edge.b), edge: edge, edgeT: edge.bT, colinear: false });
                        }
                    }
                    else if (edgeT != edge.aT) {
                        // edge crosses intersection line, neither starts nor ends on it
                        var p = edge.curve.at(edgeT);
                        assert(plane2.containsPoint(p), edge.toString(), p, edgeT, plane2.distanceToPoint(p));
                        assert(isLine.containsPoint(p), edge.toString(), p, edgeT, isLine.distanceToPoint(p));
                        var insideDir = plane2.normal.negated();
                        ps.push({ p: p, insideDir: insideDir, t: isLine.pointLambda(p), edge: edge, edgeT: edgeT, colinear: false });
                    }
                }
            }
        });
    });
    // duplicate 't's are ok, as sometimes a segment needs to stop and start again
    // should be sorted so that back facing ones are first
    ps.sort(function (a, b) { return a.t - b.t || a.insideDir.dot(isLine.dir1); });
    return ps;
}
function faceEdgeISPsWithSurface(brep, face, isCurves, surface2) {
    var surface = face.surface;
    var pss = NLA.arrayFromFunction(isCurves.length, function (i) { return []; });
    var loops = face.holes.concat([face.edges]);
    loops.forEach(function (loop) {
        var colinearEdges = loop.map(function (edge) { return isCurves.findIndex(function (curve) { return edge.curve.isColinearTo(curve); }); });
        // todo: this assumes that isCurves do not touch
        loop.forEach(function (edge, edgeIndex, edges) {
            var nextEdgeIndex = (edgeIndex + 1) % edges.length, nextEdge = edges[nextEdgeIndex];
            //console.log(edge.toSource()) {p:V3(2, -2.102, 0),
            var colinearEdgeCurveIndex = colinearEdges[edgeIndex];
            if (-1 != colinearEdgeCurveIndex) {
                var isCurve = isCurves[colinearEdgeCurveIndex];
                // edge colinear to an intersection curve
                var prevEdgeIndex = (edgeIndex - 1 + edges.length) % edges.length, prevEdge = edges[prevEdgeIndex];
                var curveAT = isCurve.pointLambda(edge.a), curveBT = isCurve.pointLambda(edge.b);
                var colinearOutA = edge.aDir.cross(surface.normalAt(edge.a));
                var ps = pss[colinearEdgeCurveIndex];
                if (-1 != colinearEdges[prevEdgeIndex] && dotCurve(colinearOutA, prevEdge.bDir, prevEdge.bDDT) < 0) {
                    ps.push({ p: prevEdge.b, insideDir: edge.aDir.negated(), t: curveAT, edge: prevEdge, edgeT: prevEdge.bT, colinear: false });
                }
                ps.push({ p: edge.a, insideDir: edge.aDir, t: curveAT, edge: edge, edgeT: edge.aT, colinear: true }, { p: edge.b, insideDir: edge.bDir.negated(), t: curveBT, edge: edge, edgeT: edge.bT, colinear: true });
                var colinearOutB = edge.bDir.cross(surface.normalAt(edge.b));
                if (-1 != colinearEdges[nextEdgeIndex] && dotCurve(colinearOutB, nextEdge.aDir, nextEdge.aDDT) > 0) {
                    ps.push({ p: edge.b, insideDir: edge.bDir, t: curveBT, edge: edge, edgeT: prevEdge.bT, colinear: false });
                }
            }
            else {
                var edgeTs = edge.edgeISTsWithSurface(surface2);
                var _loop_3 = function(k) {
                    var edgeT = edgeTs[k];
                    var p = edge.curve.at(edgeT);
                    if (!isCurves.some(function (isCurve) { return isCurve.containsPoint(p); }, edge.toString() + p + edgeT)) {
                        console.log(isCurves, isCurves[0].sce, p.sce, edge.sce, surface2.sce, face.surface.sce);
                        isCurves.forEach(function (isCurve) { return (isCurve.debugToMesh(mesh1, 'curve1')); });
                        drPs.push({ p: p, text: 'This is the problematic point' });
                        assert(false, isCurves.map(function (isc) { return isc.distanceToPoint(p); }).join(' '));
                    }
                    var isCurveIndex = void 0;
                    if (1 == isCurves.length) {
                        isCurveIndex = 0;
                    }
                    else {
                        isCurveIndex = isCurves.findIndex(function (isCurve) { return isCurve.containsPoint(p); });
                        assert(isCurves.slice(isCurveIndex + 1).every(function (isCurve) { return !isCurve.containsPoint(p); }));
                    }
                    var isCurve = isCurves[isCurveIndex];
                    var curveT = isCurve.pointLambda(p);
                    if (isNaN(curveT)) {
                        if (isCurve instanceof EllipseCurve) {
                            var hint = edge.curve.tangentAt(edgeT).cross(isCurve.f1).dot(isCurve.f2);
                            curveT = isCurve.pointLambda(p, hint);
                        }
                        else {
                            assert(false);
                        }
                    }
                    assert(!isNaN(curveT));
                    var insideDir = edge.tangentAt(edgeT).cross(surface.normalAt(p)).negated();
                    assert(!NLA.eq0(insideDir.dot(isCurve.tangentAt(curveT))));
                    // Edge.edgeISTsWithSurface returns snapped values, so comparison with == is ok:
                    if (edgeT == edge.bT) {
                        // endpoint lies on intersection line
                        if (-1 != colinearEdges[nextEdgeIndex]) {
                            if (edgesDifferentSidesOfDir(isCurve.tangentAt(curveT), edge, nextEdge)) {
                                // next segment is not colinear and ends on different side
                                console.log("adding");
                                pss[isCurveIndex].push({ p: edge.b, insideDir: insideDir, t: curveT, edge: edge, edgeT: edge.bT, colinear: false });
                            }
                        }
                    }
                    else if (edgeT != edge.aT) {
                        // edge crosses an intersection curve, neither starts nor ends on it
                        pss[isCurveIndex].push({ p: p, insideDir: insideDir, t: curveT, edge: edge, edgeT: edgeT, colinear: false });
                        console.log('middle');
                    }
                };
                for (var k = 0; k < edgeTs.length; k++) {
                    _loop_3(k);
                }
            }
        });
    });
    // duplicate 't's are ok, as sometimes a segment needs to stop and start again
    // should be sorted so that back facing ones are first
    pss.forEach(function (ps, isCurveIndex) { return ps.sort(function (a, b) { return a.t - b.t || a.insideDir.dot(isCurves[isCurveIndex].tangentAt(a.t)); }); });
    return pss;
}
function dotCurve(v, cDir, cDDT) {
    var dot = v.dot(cDir);
    if (NLA.eq0(dot)) {
        dot = v.dot(cDDT);
    }
    assert(!NLA.eq0(dot));
    return dot;
}
function edgesDifferentSidesOfDir(dir, e1, e2) {
    var factor1 = dir.dot(e1.bDir);
    if (NLA.eq0(factor1)) {
        factor1 = dir.dot(e1.bDDT);
    }
    assert(!NLA.eq0(factor1));
    var factor2 = dir.dot(e2.aDir);
    if (NLA.eq0(factor2)) {
        factor2 = dir.dot(e2.aDDT);
    }
    assert(!NLA.eq0(factor2));
    return factor1 * factor2 > 0;
}
var INSIDE = 0, OUTSIDE = 1, COPLANAR_SAME = 2, COPLANAR_OPPOSITE = 3, ALONG_EDGE_OR_PLANE = 4;
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
    assert(arguments.length == 4);
    //assert(p.equals(edge.a))
    var ab1 = edge.aDir.normalized();
    var relFaces = facesWithEdge(edge, brep.faces);
    relFaces.forEach(function (faceInfo) {
        faceInfo.normalAtEdgeA = faceInfo.face.surface.normalAt(edge.a);
        faceInfo.edgeDirAtEdgeA = !faceInfo.reversed
            ? faceInfo.edge.aDir
            : faceInfo.edge.bDir;
        faceInfo.outsideVector = faceInfo.edgeDirAtEdgeA.cross(faceInfo.normalAtEdgeA);
        faceInfo.angle = (dirAtEdgeA.angleRelativeNormal(faceInfo.outsideVector.negated(), ab1) + 2 * Math.PI + NLA_PRECISION / 2) % (2 * Math.PI);
    });
    assert(relFaces.length != 0, edge.toSource());
    relFaces.sort(function (a, b) { return a.angle - b.angle; });
    // assert(relFaces.length % 2 == 0, edge.toSource()) // even number of touching faces
    if (NLA.eq0(relFaces[0].angle)) {
        //assert(false) todo
        var coplanarSame = relFaces[0].normalAtEdgeA.dot(faceNormal) > 0;
        return coplanarSame ? COPLANAR_SAME : COPLANAR_OPPOSITE;
    }
    else {
        return !relFaces[0].reversed ? INSIDE : OUTSIDE;
    }
}
function splitsVolumeEnclosingCone(brep, p, dir) {
    var testPlane = P3.forAnchorAndPlaneVectors(p, dir, dir.getPerpendicular());
    var rays = [];
    var _loop_4 = function(k) {
        var face = brep.faces[k];
        assertf(function () { return face instanceof PlaneFace; });
        if (face.getAllEdges().some(function (edge) { return edge.a.like(p); })) {
            if (testPlane.isParallelToPlane(face.surface.plane)) {
                if (face.pointsToInside(p, dir)) {
                    return { value: ALONG_EDGE_OR_PLANE };
                }
            }
            else {
                var isLine = L3.fromPlanes(testPlane, face.surface.plane);
                var ps = planeFaceEdgeISPsWithPlane(null, face, isLine, testPlane);
                var i = 0;
                while (i < ps.length) {
                    var a = ps[i++], b = ps[i++];
                    var out = a.p.like(p);
                    if (out || b.p.like(p)) {
                        var dir2 = out ? isLine.dir1 : isLine.dir1.negated();
                        var angle = (dir.angleRelativeNormal(dir2, testPlane.normal) + 2 * Math.PI + NLA_PRECISION / 2) % (2 * Math.PI);
                        rays.push({ angle: angle, out: out });
                    }
                }
            }
        }
    };
    for (var k = 0; k < brep.faces.length; k++) {
        var state_4 = _loop_4(k);
        if (typeof state_4 === "object") return state_4.value;
    }
    rays.sort(function (a, b) { return a.angle - b.angle; });
    //console.log("testPlane", testPlane.toSource(), "rays", rays.toSource())
    if (NLA.eq0(rays[0].angle)) {
        return ALONG_EDGE_OR_PLANE;
    }
    else {
        return rays[0].out ? OUTSIDE : INSIDE;
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
    assertNumbers(a, b, c);
    // TODO: disambiguate on a < b
    var term = sqrt(a * a + b * b - c * c);
    return {
        x1: (a * c + b * term) / (a * a + b * b),
        x2: (a * c - b * term) / (a * a + b * b),
        y1: (b * c - a * term) / (a * a + b * b),
        y2: (b * c + a * term) / (a * a + b * b)
    };
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
    assertNumbers(a, b, c);
    var aa = a * a, bb = b * b, cc = c * c;
    // TODO: disambiguate on a < b
    //var xTerm = sqrt(4*cc*aa-4*(bb-aa)*(-cc-bb))
    var xTerm = 2 * sqrt(bb * cc + bb * bb - aa * bb);
    var yTerm = sqrt(4 * cc * bb - 4 * (bb - aa) * (cc - aa));
    return {
        x1: (-2 * a * c + xTerm) / 2 / (bb - aa),
        x2: (-2 * a * c - xTerm) / 2 / (bb - aa),
        y1: (2 * b * c - yTerm) / 2 / (bb - aa),
        y2: (2 * b * c + yTerm) / 2 / (bb - aa)
    };
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
    assertNumbers(a, b, c, r);
    var term = sqrt(r * r * (a * a + b * b) - c * c);
    return {
        x1: (a * c + b * term) / (a * a + b * b),
        x2: (a * c - b * term) / (a * a + b * b),
        y1: (b * c - a * term) / (a * a + b * b),
        y2: (b * c + a * term) / (a * a + b * b)
    };
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
    var step = (endT - startT) / steps;
    var length = 0;
    var p = curve.at(startT);
    var i = 0, t = startT + step;
    for (; i < steps; i++, t += step) {
        var next = curve.at(t);
        length += p.distanceTo(next);
        p = next;
    }
    return length;
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
];
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
];
function gaussLegendreQuadrature24(fn, startT, endT) {
    var result = 0;
    for (var i = 0; i < gaussLegendreXs.length; i++) {
        // gauss-legendre goes from -1 to 1, so we need to scale
        var t = startT + (gaussLegendreXs[i] + 1) / 2 * (endT - startT);
        result += gaussLegendreWeights[i] * fn(t);
    }
    // again, [-1,1], so div by 2
    return result / 2 * (endT - startT);
}
function curveLengthByDerivative(df, startT, endT, steps) {
    var dt = (endT - startT) / steps;
    var length = 0;
    var i = 0, t = startT + dt / 2;
    for (; i < steps; i++, t += dt) {
        length += df(t) * dt;
    }
    return length;
}
function curvePoint(implicitCurve, startPoint) {
    var eps = 1e-5;
    var p = startPoint;
    for (var i = 0; i < 4; i++) {
        var fp = implicitCurve(p.x, p.y);
        var dfpdx = (implicitCurve(p.x + eps, p.y) - fp) / eps, dfpdy = (implicitCurve(p.x, p.y + eps) - fp) / eps;
        var scale = fp / (dfpdx * dfpdx + dfpdy * dfpdy);
        //console.log(p.$)
        p = p.minus(new V3(scale * dfpdx, scale * dfpdy, 0));
    }
    return p;
}
function followAlgorithm(implicitCurve, startPoint, endPoint, stepLength, startDir, tangentEndPoints, boundsFunction) {
    assertNumbers(stepLength, implicitCurve(0, 0));
    assertVectors(startPoint, endPoint);
    assert(!startDir || startDir instanceof V3);
    var points = [];
    tangentEndPoints = tangentEndPoints || [];
    assert(NLA.eq0(implicitCurve(startPoint.x, startPoint.y)), 'NLA.isZero(implicitCurve(startPoint.x, startPoint.y))');
    stepLength = stepLength || 0.5;
    var eps = 1e-5;
    var p = startPoint, prevp = startDir ? p.minus(startDir) : p;
    var i = 0;
    do {
        var fp = implicitCurve(p.x, p.y);
        var dfpdx = (implicitCurve(p.x + eps, p.y) - fp) / eps, dfpdy = (implicitCurve(p.x, p.y + eps) - fp) / eps;
        var tangent = new V3(-dfpdy, dfpdx, 0);
        var reversedDir = p.minus(prevp).dot(tangent) < 0;
        tangent = tangent.toLength(reversedDir ? -stepLength : stepLength);
        var tangentEndPoint = p.plus(tangent);
        points.push(p);
        tangentEndPoints.push(tangentEndPoint);
        prevp = p;
        p = curvePoint(implicitCurve, tangentEndPoint);
    } while (i++ < 100 && (i < 4 || prevp.distanceTo(endPoint) > 1.1 * stepLength) && boundsFunction(p.x, p.x));
    // TODO gleichmäßige Verteilung der Punkte
    return points;
}
// both curves must be in the same s-t coordinates for this to make sense
function intersectionICurveICurve(pCurve1, startParams1, endParams1, startDir, stepLength, pCurve2) {
    assertNumbers(stepLength, pCurve1(0, 0), pCurve2(0, 0));
    assertVectors(startParams1, endParams1);
    assert(!startDir || startDir instanceof V3);
    var vertices = [];
    assert(NLA.eq0(pCurve1(startParams1.x, startParams1.y)));
    stepLength = stepLength || 0.5;
    var eps = 1e-5;
    var p = startParams1, prevp = p; // startDir ? p.minus(startDir) : p
    var i = 0;
    while (i++ < 1000 && (i < 4 || p.distanceTo(endParams1) > 1.1 * stepLength)) {
        var fp = pCurve1(p.x, p.y);
        var dfpdx = (pCurve1(p.x + eps, p.y) - fp) / eps, dfpdy = (pCurve1(p.x, p.y + eps) - fp) / eps;
        var tangent = new V3(-dfpdy, dfpdx, 0).toLength(stepLength);
        if (p.minus(prevp).dot(tangent) < 0)
            tangent = tangent.negated();
        prevp = p;
        p = curvePoint(pCurve1, p.plus(tangent));
        vertices.push(p);
    }
    // TODO gleichmäßige Verteilung der Punkte
    return vertices;
}
function intersectionICurveICurve2(iCurve1, loopPoints1, iCurve2) {
    var p = loopPoints1[0], val = iCurve2(p.x, p.y), lastVal;
    var iss = [];
    for (var i = 0; i < loopPoints1.length; i++) {
        lastVal = val;
        p = loopPoints1[i];
        val = iCurve2(p);
        if (val * lastVal <= 0) {
            iss.push(newtonIterate2d(iCurve1, iCurve2, p.x, p.y));
        }
    }
    return iss;
}
function intersectionPCurveISurface(parametricCurve, searchStart, searchEnd, searchStep, implicitSurface) {
    assertNumbers(searchStart, searchEnd, searchStep);
    var iss = [];
    var val = implicitSurface(parametricCurve(searchStart)), lastVal;
    for (var t = searchStart + searchStep; t <= searchEnd; t += searchStep) {
        lastVal = val;
        val = implicitSurface(parametricCurve(t));
        if (val * lastVal <= 0) {
            iss.push(newtonIterate1d(function (t) { return implicitSurface(parametricCurve(t)); }, t));
        }
    }
    return iss;
}
function intersectionICurvePSurface(f0, f1, parametricSurface) {
}
function cassini(a, c) {
    return function (x, y) { return (x * x + y * y) * (x * x + y * y) - 2 * c * c * (x * x - y * y) - (a * a * a * a - c * c * c * c); };
}
