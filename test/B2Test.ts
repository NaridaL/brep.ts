QUnit.module('B2Test')

registerTests({

    'BREP.isCCW'(assert) {
        const vertices = [V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)]
        assert.ok(isCCW(vertices, V(0, 0, 1)))
        assert.notOk(isCCW(vertices, V(0, 0, -1)))
    },
    'Face.equals'(assert) {
        const a = PlaneFace.forVertices(P3.XY, [V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)])
        const b = PlaneFace.forVertices(P3.XY, [V(0, 10, 0), V(0, 0, 0), V(10, 0, 0), V(10, 10, 0)])
        const c = PlaneFace.forVertices(new P3(V(0, 0, -1), 0), [V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)].slice().reverse())
        assert.ok(a.equals(a))

        assert.ok(a.equals(b))
        assert.ok(b.equals(a))

        assert.notOk(a.equals(c))
        assert.notOk(c.equals(a))

        assert.notOk(b.equals(c))
        assert.notOk(c.equals(b))
    },
    'B2.like'(assert) {
        const a = B2T.tetrahedron(V(5, 5, 5), V(5, 5, -5), V(10, 12, 1), V(0, 12, 1))
        const b = B2T.tetrahedron(V(5, 5, 5), V(10, 12, 1), V(5, 5, -5), V(0, 12, 1))
        const c = B2T.tetrahedron(V(5, 5, 5), V(12, 12, 1), V(5, 5, -5), V(0, 12, 1))

        assert.ok(a.like(a))

        assert.ok(a.like(b))
        assert.ok(b.like(a))

        assert.notOk(a.like(c))
        assert.notOk(c.like(a))

        assert.notOk(b.like(c))
        assert.notOk(c.like(b))
    },

    'B2T.tetrahedron'(assert) {
        const a = B2T.tetrahedron(V(5, 5, 5), V(5, 5, -5), V(9, 11, 1), V(1, 9, 0))
        const result = new B2([
            PlaneFace.forVertices(new P3(V(0.8320502943378437, -0.5547001962252291, 0), 1.386750490563073), [V(5, 5, 5), V(5, 5, -5), V(9, 11, 1)]),
            PlaneFace.forVertices(new P3(V(-0.7071067811865475, -0.7071067811865475, 0), -7.071067811865475), [V(5, 5, 5), V(1, 9, 0), V(5, 5, -5)]),
            PlaneFace.forVertices(new P3(V(-0.1003911722115382, 0.7362019295512802, -0.6692744814102547), 6.525426193749983), [V(5, 5, -5), V(1, 9, 0), V(9, 11, 1)]),
            PlaneFace.forVertices(new P3(V(-0.25177250044296223, 0.6474150011390457, 0.7193500012656064), 5.57496250980845), [V(9, 11, 1), V(1, 9, 0), V(5, 5, 5)])], false)
        assert.B2equals(a, result)
    },


    'Face.prototype.containsPoint'(assert) {
        let a = PlaneFace.forVertices(P3.XY, [V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)])
        assert.ok(a.containsPoint(V(5, 5, 0)))
        assert.notOk(a.containsPoint(V(11, 5, 0)))


        let b = PlaneFace.forVertices(P3.XY, [V(0, 10, 0), V(0, 0, 0), V(5, 0, 0), V(6, 5, 0)])
        assert.ok(b.containsPoint(V(2, 5, 0)))

        let c = PlaneFace.forVertices(P3.XY, [V(0, 10, 0), V(0, 0, 0), V(5, 0, 0), V(6, 5, 0), V(10, 5, 0)])
        assert.ok(c.containsPoint(V(2, 5, 0)))


        a = a.rotateZ(30 * DEG)
        const m = M4.rotationZ(30 * DEG)
        console.log(a.toString())
        assert.ok(a.containsPoint(m.transformPoint(V(5, 5, 0))))
        assert.notOk(a.containsPoint(m.transformPoint(V(-5, 5, 0))))

        b = b.rotateZ(30 * DEG)
        assert.ok(b.containsPoint(m.transformPoint(V(2, 5, 0))))

        c = c.rotateZ(30 * DEG)
        assert.ok(c.containsPoint(m.transformPoint(V(2, 5, 0))))
    },

    'Face.prototype.containsPoint 2'(assert) {
        const a = PlaneFace.forVertices(P3.XY, [V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)])
        assert.notOk(a.containsPoint(V(-0.00000001, 11, 0)))

        assert.ok(PlaneFace.forVertices(new P3(V(0, 0, -1), 0), [V(-1, -10, 0), V(0, 25, 0), V(25, 0, 0)]).containsPoint(V(0, 0, 0)))
    },
    'Face.withHole'(assert) {
        const a = PlaneFace.forVertices(P3.XY, [V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)])
        let holeVertices = [V(2, 3, 0), V(8, 7, 0), V(7, 2, 0)]


        assert.notOk(a.containsPoint(V(-0.00000001, 11, 0)))

    },

    'splitsVolumeEnclosingFaces'(assert) {
        const brep = B2T.tetrahedron(V(0, 0, 0), V(10, 0, 0), V(0, 10, 0), V(0, 0, 10))
        brep.buildAdjacencies()
        // pointing into tetrahedon
        const edge = (a, b) => brep.faces.map(face => face.getAllEdges()).concatenated().find(edge => edge.a.like(a) && edge.b.like(b)).getCanon()
        assert.equal(splitsVolumeEnclosingFaces(brep, edge(V(0, 0, 0), V(10, 0, 0)), V(0, 1, 1), V(0, -1, 1)), INSIDE)
        assert.equal(splitsVolumeEnclosingFaces(brep, edge(V(0, 0, 0), V(10, 0, 0)), V(0, 1, 1), V(0, 1, -1)), INSIDE)

        // pointing out of tetrahedon
        assert.equal(splitsVolumeEnclosingFaces(brep, edge(V(0, 0, 0), V(10, 0, 0)), V(0, -1, 0), V(0, 1, 1)), OUTSIDE)
        assert.equal(splitsVolumeEnclosingFaces(brep, edge(V(0, 0, 0), V(10, 0, 0)), V(0, -1, -1), V(0, -1, 1)), OUTSIDE)
        assert.equal(splitsVolumeEnclosingFaces(brep, edge(V(0, 0, 0), V(10, 0, 0)), V(0, -1, -1), V(0, 1, -1)), OUTSIDE)
    },


    'splitsVolumeEnclosingFaces 2'(assert) {
        const brep = B2T.tetrahedron(V(0, 0, 0), V(10, 0, 0), V(0, 10, 0), V(0, 0, 10))
        brep.buildAdjacencies()
        // pointing out of tetrahedon
        const edge = (a, b) => brep.faces.map(face => face.getAllEdges()).concatenated().find(edge => edge.a.like(a) && edge.b.like(b)).getCanon()
        assert.equal(splitsVolumeEnclosingFaces(brep, edge(V(0, 0, 0), V(10, 0, 0)), V(0, 1, 0), V(0, 0, 1)), COPLANAR_OPPOSITE)
        assert.equal(splitsVolumeEnclosingFaces(brep, edge(V(0, 0, 0), V(10, 0, 0)), V(0, 1, 0), V(0, 0, -1)), COPLANAR_SAME)
        assert.equal(splitsVolumeEnclosingFaces(brep, edge(V(0, 0, 0), V(10, 0, 0)), V(0, 0, 1), V(0, 1, 0)), COPLANAR_OPPOSITE)
    },
    'splitsVolumeEnclosingFaces 3'(assert) {
        const brep = B2T.box(5, 5, 5).flipped()
        brep.buildAdjacencies()
        console.log(brep.sce)

        const edge = (a, b) => brep.faces.map(face => face.getAllEdges()).concatenated().find(edge => edge.a.like(a) && edge.b.like(b)).getCanon()
        assert.equal(splitsVolumeEnclosingFaces(brep, edge(V(0, 5, 0), V(0, 0, 0)), V(0, 0, -1), V(1, 0, 0)), INSIDE)
        assert.equal(splitsVolumeEnclosingFaces(brep, edge(V(0, 0, 0), V(0, 5, 0)), V(0, 0, -1), V(1, 0, 0)), INSIDE)

        assert.equal(splitsVolumeEnclosingFaces(brep, edge(V(0, 0, 0), V(0, 5, 0)), V(0, 0, 1), V(-1, 0, 0)), COPLANAR_OPPOSITE)
        assert.equal(splitsVolumeEnclosingFaces(brep, edge(V(0, 0, 0), V(0, 5, 0)), V(0, 0, 1), V(1, 0, 0)), COPLANAR_SAME)
    },
    'splitsVolumeEnclosingCone'(assert) {
        const brep = B2T.box(5, 5, 5)

        assert.equal(splitsVolumeEnclosingCone(brep, V3.ZERO, V(1, 1, 1)), INSIDE)
        assert.equal(splitsVolumeEnclosingCone(brep, V3.ZERO, V(-1, 1, 1)), OUTSIDE)
        assert.equal(splitsVolumeEnclosingCone(brep, V3.ZERO, V3.X), ALONG_EDGE_OR_PLANE)
        assert.equal(splitsVolumeEnclosingCone(brep, V3.ZERO, V3.X.negated()), OUTSIDE)
        assert.equal(splitsVolumeEnclosingCone(brep, V(0, 5, 5), V3.Y), OUTSIDE)
        assert.equal(splitsVolumeEnclosingCone(brep, V3.ZERO, V(1, 1, 0)), ALONG_EDGE_OR_PLANE)
    },

    'pointsToInside'(assert) {
        const face = PlaneFace.forVertices(P3.XY, [V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)])

        assert.equal(face.pointsToInside(V3.ZERO, V3.X.negated()), PointVsFace.OUTSIDE)

        const v = V(10, 10, 0)
        assert.equal(face.pointsToInside(v, V(-1, -1, 0)), PointVsFace.INSIDE)
        assert.equal(face.pointsToInside(v, V(1, 1, 0)), PointVsFace.OUTSIDE)
        assert.equal(face.pointsToInside(v, V(-1, 0, 0)), PointVsFace.ON_EDGE)
        assert.equal(face.pointsToInside(v, V(0, -1, 0)), PointVsFace.ON_EDGE)


        const face2 = new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 0)), [
            StraightEdge.throughPoints(V(0, 0, 0), V(0, 5, 0)),
            StraightEdge.throughPoints(V(0, 5, 0), V(5, 5, 0)),
            StraightEdge.throughPoints(V(5, 5, 0), V(5, 0, 0)),
            StraightEdge.throughPoints(V(5, 0, 0), V(0, 0, 0))], [])


    },
    'planeFaceEdgeISPsWithPlane'(assert) {
        const brep = B2T.extrudeVertices([
                V(3.318793683290006e-10, 9.12934465653568e-10, 3.8212177478055603e-10), // 0 0 0
                V(9.999999999857609, 4.561232401125957e-13, 5.033029171552123e-10), // 10 0 0
                V(2.984222284459465e-10, 9.999999999852161, -1.0461618780772852e-9)], // 0 10 0
            P3.XY, V(-1.8192047548394726e-10, -3.0150747681012967e-10, -5.0000000009123795), 'ex0')
        /*var result = planeFaceEdgeISPsWithPlane(brep,
         new BREP.Face([
         V(9.999999999675689, -3.010513535700171e-10, -5.000000000409076), // 10 0 -5
         V(1.4995889284505332e-10, 6.114269888434384e-10, -5.000000000530258), // 0 0 -5
         V(3.318793683290006e-10, 9.12934465653568e-10, 3.8212177478055603e-10), // 0 0 0
         V(9.999999999857609, 4.561232401125957e-13, 5.033029171552123e-10)], // 10 0 0
         new P3(V(9.12478342464039e-11, 1, -6.03014953543423e-11), 9.129344656608087e-10)), // 0 1 0
         L3(V(-1.3833878355530264e-10, 6.114269894465992e-10, -4.999999990964091), V(-1, 9.12478342480723e-11, 2.7667756772219476e-11)),
         new P3(V(2.766775686256173e-11, 9.90075577448337e-10, 1), -4.999999990964091),
         true, true)
         assert.deepEqual(result, [])*/
        console.log(brep.faces[2].toSource())
        let line = L3.X.translate(0, 0, -1)
        let result = planeFaceEdgeISPsWithPlane(brep, brep.faces[2], L3.Y.translate(0, 0, -1), P3.XY.translate(0, 0, -1), true, true, new NLA.CustomSet()).map(is => is.p)
        assert.compareV3arraysLike(result, [V(0, 0, -1), V(0, 10, -1)])
        result = planeFaceEdgeISPsWithPlane(brep, brep.faces[2], L3.Y, P3.XY, true, true, new NLA.CustomSet()).map(is => is.p)
        assert.compareV3arraysLike(result, [V(0, 0, 0), V(0, 10, 0)])
        result = planeFaceEdgeISPsWithPlane(brep.translate(0, 0, 10), brep.translate(0, 0, 10).faces[2], L3.Y.translate(0, 0, 6), P3.XY.translate(0, 0, 6), true, true, new NLA.CustomSet()).map(is => is.p)
        assert.compareV3arraysLike(result, [V(0, 0, 6), V(0, 10, 6)])
    },



    'planeFaceEdgeISPsWithPlane 2'(assert) {
        const brep = BREP.extrude([
                V(3.318793683290006e-10, 9.12934465653568e-10, 3.8212177478055603e-10), // 0 0 0
                V(9.999999999857609, 4.561232401125957e-13, 5.033029171552123e-10), // 10 0 0
                V(2.984222284459465e-10, 9.999999999852161, -1.0461618780772852e-9)], // 0 10 0
            P3.XY, V(-1.8192047548394726e-10, -3.0150747681012967e-10, -5.0000000009123795), 'ex0').flipped()
        const result = planeFaceEdgeISPsWithPlane(brep,
            new BREP.Face([
                    V(9.999999999675689, -3.010513535700171e-10, -5.000000000409076), // 10 0 -5
                    V(1.4995889284505332e-10, 6.114269888434384e-10, -5.000000000530258), // 0 0 -5
                    V(3.318793683290006e-10, 9.12934465653568e-10, 3.8212177478055603e-10), // 0 0 0
                    V(9.999999999857609, 4.561232401125957e-13, 5.033029171552123e-10)], // 10 0 0
                new P3(V(9.12478342464039e-11, 1, -6.03014953543423e-11), 9.129344656608087e-10)), // 0 1 0
            L3(V(-1.3833878355530264e-10, 6.114269894465992e-10, -4.999999990964091), V(-1, 9.12478342480723e-11, 2.7667756772219476e-11)),
            new P3(V(2.766775686256173e-11, 9.90075577448337e-10, 1), -4.999999990964091),
            true, true, new NLA.CustomSet()).map(is => is.p)
        assert.compareV3arraysLike(result, [])
        console.log(brep.faces[2])
        const result2 = planeFaceEdgeISPsWithPlane(brep, brep.faces[2], L3.X, P3.XY, true, true, new NLA.CustomSet())
        assert.compareV3arraysLike(result2, [])
        console.log(result2)
    },

    'B2.prototype.minus B2T.box(5, 5, 5).minus(B2T.box(1, 1, 6))'(assert) {
        const a = B2T.box(5, 5, 5, 'a')
        const b = B2T.box(1, 1, 6, 'b')
        const result = new B2([
            PlaneFace.forVertices(new P3(V(0, 0, -1), 0), [V(0, 1, 0), V(0, 5, 0), V(5, 5, 0), V(5, 0, 0), V(1, 0, 0), V(1, 1, 0)]),
            PlaneFace.forVertices(new P3(V(0, 0, 1), 5), [V(0, 1, 5), V(1, 1, 5), V(1, 0, 5), V(5, 0, 5), V(5, 5, 5), V(0, 5, 5)]),
            PlaneFace.forVertices(new P3(V(-1, 0, 0), 0), [V(0, 1, 0), V(0, 1, 5), V(0, 5, 5), V(0, 5, 0)]),
            PlaneFace.forVertices(new P3(V(0, 1, 0), 5), [V(5, 5, 0), V(0, 5, 0), V(0, 5, 5), V(5, 5, 5)]),
            PlaneFace.forVertices(new P3(V(1, 0, 0), 5), [V(5, 0, 0), V(5, 5, 0), V(5, 5, 5), V(5, 0, 5)]),
            PlaneFace.forVertices(new P3(V(0, -1, 0), 0), [V(1, 0, 0), V(5, 0, 0), V(5, 0, 5), V(1, 0, 5)]),
            PlaneFace.forVertices(new P3(V(0, -1, 0), -1), [V(0, 1, 0), V(1, 1, 0), V(1, 1, 5), V(0, 1, 5)]),
            PlaneFace.forVertices(new P3(V(-1, 0, 0), -1), [V(1, 1, 0), V(1, 0, 0), V(1, 0, 5), V(1, 1, 5)])], false)
        assert.b2Equal(a, b, a.minus(b), result)
    },
    'B2.prototype.minus B2T.box(5, 5, 5).minus(B2T.box(1, 1, 5))'(assert) {
        const a = B2T.box(5, 5, 5)
        const b = B2T.box(1, 1, 5)
        const result = new B2([
            PlaneFace.forVertices(new P3(V(0, 0, -1), 0), [V(0, 1, 0), V(0, 5, 0), V(5, 5, 0), V(5, 0, 0), V(1, 0, 0), V(1, 1, 0)]),
            PlaneFace.forVertices(new P3(V(0, 0, 1), 5), [V(0, 1, 5), V(1, 1, 5), V(1, 0, 5), V(5, 0, 5), V(5, 5, 5), V(0, 5, 5)]),
            PlaneFace.forVertices(new P3(V(-1, 0, 0), 0), [V(0, 1, 0), V(0, 1, 5), V(0, 5, 5), V(0, 5, 0)]),
            PlaneFace.forVertices(new P3(V(0, 1, 0), 5), [V(5, 5, 0), V(0, 5, 0), V(0, 5, 5), V(5, 5, 5)]),
            PlaneFace.forVertices(new P3(V(1, 0, 0), 5), [V(5, 0, 0), V(5, 5, 0), V(5, 5, 5), V(5, 0, 5)]),
            PlaneFace.forVertices(new P3(V(0, -1, 0), 0), [V(1, 0, 0), V(5, 0, 0), V(5, 0, 5), V(1, 0, 5)]),
            PlaneFace.forVertices(new P3(V(0, -1, 0), -1), [V(0, 1, 0), V(1, 1, 0), V(1, 1, 5), V(0, 1, 5)]),
            PlaneFace.forVertices(new P3(V(-1, 0, 0), -1), [V(1, 1, 0), V(1, 0, 0), V(1, 0, 5), V(1, 1, 5)])], false)
        assert.b2Equal(a, b, a.minus(b), result)
    },
    'B2 edge/face intersection'(assert) {
        const wideBox = B2T.box(10, 10, 5)
        const extrusion = B2T.extrudeVertices([V(1, 0), V(0, 3), V(-2, 5), V(5, 5)], P3.XY.flipped(), V(0, 0, 10), 'lol')
            .translate(0, 1, 1)
        const result = new B2([
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 5)), [
                StraightEdge.throughPoints(V(0, 4, 5), V(0, 0, 5)),
                StraightEdge.throughPoints(V(0, 0, 5), V(10, 0, 5)),
                StraightEdge.throughPoints(V(10, 0, 5), V(10, 10, 5)),
                StraightEdge.throughPoints(V(10, 10, 5), V(0, 10, 5)),
                StraightEdge.throughPoints(V(0, 10, 5), V(0, 6, 5)),
                StraightEdge.throughPoints(V(0, 6, 5), V(5, 6, 5)),
                StraightEdge.throughPoints(V(5, 6, 5), V(1, 1, 5)),
                StraightEdge.throughPoints(V(1, 1, 5), V(0, 4, 5))]),
            new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), 0)), [
                StraightEdge.throughPoints(V(0, 6, 1), V(0, 6, 5)),
                StraightEdge.throughPoints(V(0, 6, 5), V(0, 10, 5)),
                StraightEdge.throughPoints(V(0, 10, 5), V(0, 10, 0)),
                StraightEdge.throughPoints(V(0, 10, 0), V(0, 0, 0)),
                StraightEdge.throughPoints(V(0, 0, 0), V(0, 0, 5)),
                StraightEdge.throughPoints(V(0, 0, 5), V(0, 4.000000000000002, 5)),
                StraightEdge.throughPoints(V(0, 4, 5), V(0, 4, 1)),
                StraightEdge.throughPoints(V(0, 4, 1), V(0, 6, 1))]),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 0)), [
                StraightEdge.throughPoints(V(0, 0, 0), V(0, 10, 0)),
                StraightEdge.throughPoints(V(0, 10, 0), V(10, 10, 0)),
                StraightEdge.throughPoints(V(10, 10, 0), V(10, 0, 0)),
                StraightEdge.throughPoints(V(10, 0, 0), V(0, 0, 0))]),
            new PlaneFace(new PlaneSurface(new P3(V(0, 1, 0), 10)), [
                StraightEdge.throughPoints(V(10, 10, 0), V(0, 10, 0)),
                StraightEdge.throughPoints(V(0, 10, 0), V(0, 10, 5)),
                StraightEdge.throughPoints(V(0, 10, 5), V(10, 10, 5)),
                StraightEdge.throughPoints(V(10, 10, 5), V(10, 10, 0))]),
            new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 10)), [
                StraightEdge.throughPoints(V(10, 0, 0), V(10, 10, 0)),
                StraightEdge.throughPoints(V(10, 10, 0), V(10, 10, 5)),
                StraightEdge.throughPoints(V(10, 10, 5), V(10, 0, 5)),
                StraightEdge.throughPoints(V(10, 0, 5), V(10, 0, 0))]),
            new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0)), [
                StraightEdge.throughPoints(V(0, 0, 0), V(10, 0, 0)),
                StraightEdge.throughPoints(V(10, 0, 0), V(10, 0, 5)),
                StraightEdge.throughPoints(V(10, 0, 5), V(0, 0, 5)),
                StraightEdge.throughPoints(V(0, 0, 5), V(0, 0, 0))]),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), -1)), [
                StraightEdge.throughPoints(V(0, 4, 1), V(-2, 6, 1)),
                StraightEdge.throughPoints(V(-2, 6, 1), V(0, 6, 1)),
                StraightEdge.throughPoints(V(0, 6, 1), V(0, 4, 1))]),
            new PlaneFace(new PlaneSurface(new P3(V(-0.9486832980505139, -0.316227766016838, 0), -1.2649110640673518)), [
                StraightEdge.throughPoints(V(1, 1, 5), V(1, 1, 11)),
                StraightEdge.throughPoints(V(1, 1, 11), V(0, 4, 11)),
                StraightEdge.throughPoints(V(0, 4, 11), V(0, 4, 5)),
                StraightEdge.throughPoints(V(0, 4, 5), V(1, 1, 5))]),
            new PlaneFace(new PlaneSurface(new P3(V(-0.7071067811865477, -0.7071067811865472, 0), -2.8284271247461894)), [
                StraightEdge.throughPoints(V(0, 4, 5), V(0, 4, 11)),
                StraightEdge.throughPoints(V(0, 4, 11), V(-2, 6, 11)),
                StraightEdge.throughPoints(V(-2, 6, 11), V(-2, 6, 1)),
                StraightEdge.throughPoints(V(-2, 6, 1), V(0, 4, 1)),
                StraightEdge.throughPoints(V(0, 4, 1), V(0, 4, 5))]),
            new PlaneFace(new PlaneSurface(new P3(V(0, 1, 0), 6)), [
                StraightEdge.throughPoints(V(0, 6, 5), V(0, 6, 1)),
                StraightEdge.throughPoints(V(0, 6, 1), V(-2, 6, 1)),
                StraightEdge.throughPoints(V(-2, 6, 1), V(-2, 6, 11)),
                StraightEdge.throughPoints(V(-2, 6, 11), V(5, 6, 11)),
                StraightEdge.throughPoints(V(5, 6, 11), V(5, 6, 5)),
                StraightEdge.throughPoints(V(5, 6, 5), V(0, 6, 5))]),
            new PlaneFace(new PlaneSurface(new P3(V(0.7808688094430303, -0.6246950475544244, 0), 0.15617376188860604)), [
                StraightEdge.throughPoints(V(5, 6, 5), V(5, 6, 11)),
                StraightEdge.throughPoints(V(5, 6, 11), V(1, 1, 11)),
                StraightEdge.throughPoints(V(1, 1, 11), V(1, 1, 5)),
                StraightEdge.throughPoints(V(1, 1, 5), V(5, 6, 5))]),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 11)), [
                StraightEdge.throughPoints(V(5, 6, 11), V(-2, 6, 11)),
                StraightEdge.throughPoints(V(-2, 6, 11), V(0, 4, 11)),
                StraightEdge.throughPoints(V(0, 4, 11), V(1, 1, 11)),
                StraightEdge.throughPoints(V(1, 1, 11), V(5, 6, 11))])], true)
        assert.b2Equal(wideBox, extrusion, wideBox.plus(extrusion), result)
    },
    'B2T.box(10, 10, 5) + B2T.box(10, 10, 5).translate(0, 0, 5) touching coplanar faces result contains both'(assert) {
        const a = B2T.box(10, 10, 5), b = a.translate(0, 0, 5)

        // the expected result is the union of a and b's faces, minus the ones where they touch, which we remove with an AABB test
        const testAABB = AABB.forXYZ(12, 12, 2).translate(-1, -1, 4)
        const expected = new B2(a.faces.concat(b.faces).filter(face => !testAABB.containsAABB(face.getAABB())), false)
        assert.b2Equal(a, b, b.plus(a), expected)
    },
    'B2T.box(10, 10, 5) - B2T.box(10, 5, 5).translate(0, 0, 5) - B2T.box(5, 5, 5).translate(5, 5, 5)'(assert) {
        const box = B2T.box(10, 10, 10), box2 = B2T.box(10, 5, 5).translate(0, 0, 5), box3 = B2T.box(5, 5, 5).translate(5, 5, 5)
        const actual = box.minus(box2).minus(box3)
        const expected = new B2([
            new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 10)), [
                StraightEdge.throughPoints(V(10, 10, 5), V(10, 5, 5)),
                StraightEdge.throughPoints(V(10, 5, 5), V(10, 0, 5)),
                StraightEdge.throughPoints(V(10, 0, 5), V(10, 0, 0)),
                StraightEdge.throughPoints(V(10, 0, 0), V(10, 10, 0)),
                StraightEdge.throughPoints(V(10, 10, 0), V(10, 10, 5))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 10)), [
                StraightEdge.throughPoints(V(5, 5, 10), V(5, 10, 10)),
                StraightEdge.throughPoints(V(5, 10, 10), V(0, 10, 10)),
                StraightEdge.throughPoints(V(0, 10, 10), V(0, 5, 10)),
                StraightEdge.throughPoints(V(0, 5, 10), V(5, 5, 10))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 1, 0), 10)), [
                StraightEdge.throughPoints(V(5, 10, 10), V(5, 10, 5)),
                StraightEdge.throughPoints(V(5, 10, 5), V(10, 10, 5)),
                StraightEdge.throughPoints(V(10, 10, 5), V(10, 10, 0)),
                StraightEdge.throughPoints(V(10, 10, 0), V(0, 10, 0)),
                StraightEdge.throughPoints(V(0, 10, 0), V(0, 10, 10)),
                StraightEdge.throughPoints(V(0, 10, 10), V(5, 10, 10))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), -5)), [
                StraightEdge.throughPoints(V(5, 5, 5), V(5, 5, 10)),
                StraightEdge.throughPoints(V(5, 5, 10), V(0, 5, 10)),
                StraightEdge.throughPoints(V(0, 5, 10), V(0, 5, 5)),
                StraightEdge.throughPoints(V(0, 5, 5), V(5, 5, 5))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 5)), [
                StraightEdge.throughPoints(V(10, 5, 5), V(5, 5, 5)),
                StraightEdge.throughPoints(V(5, 5, 5), V(0, 5, 5)),
                StraightEdge.throughPoints(V(0, 5, 5), V(0, 0, 5)),
                StraightEdge.throughPoints(V(0, 0, 5), V(10, 0, 5)),
                StraightEdge.throughPoints(V(10, 0, 5), V(10, 5, 5))], []),
            new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), 0)), [
                StraightEdge.throughPoints(V(0, 5, 5), V(0, 5, 10)),
                StraightEdge.throughPoints(V(0, 5, 10), V(0, 10, 10)),
                StraightEdge.throughPoints(V(0, 10, 10), V(0, 10, 0)),
                StraightEdge.throughPoints(V(0, 10, 0), V(0, 0, 0)),
                StraightEdge.throughPoints(V(0, 0, 0), V(0, 0, 5)),
                StraightEdge.throughPoints(V(0, 0, 5), V(0, 5, 5))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0)), [
                StraightEdge.throughPoints(V(10, 0, 5), V(0, 0, 5)),
                StraightEdge.throughPoints(V(0, 0, 5), V(0, 0, 0)),
                StraightEdge.throughPoints(V(0, 0, 0), V(10, 0, 0)),
                StraightEdge.throughPoints(V(10, 0, 0), V(10, 0, 5))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 0)), [
                StraightEdge.throughPoints(V(0, 0, 0), V(0, 10, 0)),
                StraightEdge.throughPoints(V(0, 10, 0), V(10, 10, 0)),
                StraightEdge.throughPoints(V(10, 10, 0), V(10, 0, 0)),
                StraightEdge.throughPoints(V(10, 0, 0), V(0, 0, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 5)), [
                StraightEdge.throughPoints(V(5, 10, 5), V(5, 10, 10)),
                StraightEdge.throughPoints(V(5, 10, 10), V(5, 5, 10)),
                StraightEdge.throughPoints(V(5, 5, 10), V(5, 5, 5)),
                StraightEdge.throughPoints(V(5, 5, 5), V(5, 10, 5))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 5)), [
                StraightEdge.throughPoints(V(5, 5, 5), V(10, 5, 5)),
                StraightEdge.throughPoints(V(10, 5, 5), V(10, 10, 5)),
                StraightEdge.throughPoints(V(10, 10, 5), V(5, 10, 5)),
                StraightEdge.throughPoints(V(5, 10, 5), V(5, 5, 5))], [])], false)
        assert.b2Equal(box.minus(box2), box3, actual, expected)
    },

    'B2T.box(10, 10, 5) && B2T.box(10, 10, 5).translate(3, 3, 0) overlapping faces result contains intersection of faces'(assert) {
        let box = B2T.box(10, 10, 5, 'a')
        let box2 = B2T.box(10, 10, 4, 'b').translate(3, 3, 0)
        assert.b2Equal(box, box2, box.intersection(box2, true, true), B2T.box(7, 7, 4).translate(3, 3, 0))
    },

    'B2T.box(10, 10, 5) + B2T.box(10, 10, 5).translate(3, 3, 0) overlapping faces result contains union of faces'(assert) {
        const boxA = B2T.box(10, 10, 5, 'boxA').flipped(), boxB = boxA.translate(3, 3, 0)
        //const result = B2T.extrudeVertices([V3.ZERO, V(0, 10), V(3, 10), V(3, 13), V(13, 13), V(13, 3), V(10, 3), V(10, 0)], P3.XY.flipped(), V(0, 0, 5), 'result').flipped()
        const result = new B2([
            new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), -10)), [
                StraightEdge.throughPoints(V(3, 10, 0), V(3, 10, 5)),
                StraightEdge.throughPoints(V(3, 10, 5), V(0, 10, 5)),
                StraightEdge.throughPoints(V(0, 10, 5), V(0, 10, 0)),
                StraightEdge.throughPoints(V(0, 10, 0), V(3, 10, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), -10)), [
                StraightEdge.throughPoints(V(10, 3, 5), V(10, 3, 0)),
                StraightEdge.throughPoints(V(10, 3, 0), V(10, 0, 0)),
                StraightEdge.throughPoints(V(10, 0, 0), V(10, 0, 5)),
                StraightEdge.throughPoints(V(10, 0, 5), V(10, 3, 5))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 0)), [
                StraightEdge.throughPoints(V(3, 3, 0), V(3, 10, 0)),
                StraightEdge.throughPoints(V(3, 10, 0), V(0, 10, 0)),
                StraightEdge.throughPoints(V(0, 10, 0), V(0, 0, 0)),
                StraightEdge.throughPoints(V(0, 0, 0), V(10, 0, 0)),
                StraightEdge.throughPoints(V(10, 0, 0), V(10, 3, 0)),
                StraightEdge.throughPoints(V(10, 3, 0), V(3, 3, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), -5)), [
                StraightEdge.throughPoints(V(3, 10, 5), V(3, 3, 5)),
                StraightEdge.throughPoints(V(3, 3, 5), V(10, 3, 5)),
                StraightEdge.throughPoints(V(10, 3, 5), V(10, 0, 5)),
                StraightEdge.throughPoints(V(10, 0, 5), V(0, 0, 5)),
                StraightEdge.throughPoints(V(0, 0, 5), V(0, 10, 5)),
                StraightEdge.throughPoints(V(0, 10, 5), V(3, 10, 5))], []),
            new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 0)), [
                StraightEdge.throughPoints(V(0, 10, 0), V(0, 10, 5)),
                StraightEdge.throughPoints(V(0, 10, 5), V(0, 0, 5)),
                StraightEdge.throughPoints(V(0, 0, 5), V(0, 0, 0)),
                StraightEdge.throughPoints(V(0, 0, 0), V(0, 10, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 1, 0), 0)), [
                StraightEdge.throughPoints(V(0, 0, 0), V(0, 0, 5)),
                StraightEdge.throughPoints(V(0, 0, 5), V(10, 0, 5)),
                StraightEdge.throughPoints(V(10, 0, 5), V(10, 0, 0)),
                StraightEdge.throughPoints(V(10, 0, 0), V(0, 0, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 3)), [
                StraightEdge.throughPoints(V(3, 10, 5), V(3, 10, 0)),
                StraightEdge.throughPoints(V(3, 10, 0), V(3, 13, 0)),
                StraightEdge.throughPoints(V(3, 13, 0), V(3, 13, 5)),
                StraightEdge.throughPoints(V(3, 13, 5), V(3, 10, 5))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 1, 0), 3)), [
                StraightEdge.throughPoints(V(10, 3, 0), V(10, 3, 5)),
                StraightEdge.throughPoints(V(10, 3, 5), V(13, 3, 5)),
                StraightEdge.throughPoints(V(13, 3, 5), V(13, 3, 0)),
                StraightEdge.throughPoints(V(13, 3, 0), V(10, 3, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 0)), [
                StraightEdge.throughPoints(V(10, 10, 0), V(3, 10, 0)),
                StraightEdge.throughPoints(V(3, 10, 0), V(3, 3, 0)),
                StraightEdge.throughPoints(V(3, 3, 0), V(10, 3, 0)),
                StraightEdge.throughPoints(V(10, 3, 0), V(10, 10, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 0)), [
                StraightEdge.throughPoints(V(3, 10, 0), V(10, 10, 0)),
                StraightEdge.throughPoints(V(10, 10, 0), V(10, 3, 0)),
                StraightEdge.throughPoints(V(10, 3, 0), V(13, 3, 0)),
                StraightEdge.throughPoints(V(13, 3, 0), V(13, 13, 0)),
                StraightEdge.throughPoints(V(13, 13, 0), V(3, 13, 0)),
                StraightEdge.throughPoints(V(3, 13, 0), V(3, 10, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), -5)), [
                StraightEdge.throughPoints(V(3, 10, 5), V(10, 10, 5)),
                StraightEdge.throughPoints(V(10, 10, 5), V(10, 3, 5)),
                StraightEdge.throughPoints(V(10, 3, 5), V(3, 3, 5)),
                StraightEdge.throughPoints(V(3, 3, 5), V(3, 10, 5))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), -5)), [
                StraightEdge.throughPoints(V(10, 10, 5), V(3, 10, 5)),
                StraightEdge.throughPoints(V(3, 10, 5), V(3, 13, 5)),
                StraightEdge.throughPoints(V(3, 13, 5), V(13, 13, 5)),
                StraightEdge.throughPoints(V(13, 13, 5), V(13, 3, 5)),
                StraightEdge.throughPoints(V(13, 3, 5), V(10, 3, 5)),
                StraightEdge.throughPoints(V(10, 3, 5), V(10, 10, 5))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), -13)), [
                StraightEdge.throughPoints(V(13, 13, 0), V(13, 13, 5)),
                StraightEdge.throughPoints(V(13, 13, 5), V(3, 13, 5)),
                StraightEdge.throughPoints(V(3, 13, 5), V(3, 13, 0)),
                StraightEdge.throughPoints(V(3, 13, 0), V(13, 13, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), -13)), [
                StraightEdge.throughPoints(V(13, 3, 0), V(13, 3, 5)),
                StraightEdge.throughPoints(V(13, 3, 5), V(13, 13, 5)),
                StraightEdge.throughPoints(V(13, 13, 5), V(13, 13, 0)),
                StraightEdge.throughPoints(V(13, 13, 0), V(13, 3, 0))], [])], true)
        assert.b2Equal(boxA, boxB, boxA.intersection(boxB, true, true), result)
    },


    'B2T.box(10, 10, 5) + B2T.box(4, 10, 2).translate(2, 0, 3) overlapping faces result contains union of faces 2'(assert) {
        const box = B2T.box(10, 10, 5), box2 = B2T.box(4, 10, 2).translate(2, 0, 3)
        const expected = new B2([
            new PlaneFace(new PlaneSurface(new P3(V(0, 1, 0), 10)), [
                StraightEdge.throughPoints(V(2, 10, 3), V(6, 10, 3)),
                StraightEdge.throughPoints(V(6, 10, 3), V(6, 10, 5)),
                StraightEdge.throughPoints(V(6, 10, 5), V(10, 10, 5)),
                StraightEdge.throughPoints(V(10, 10, 5), V(10, 10, 0)),
                StraightEdge.throughPoints(V(10, 10, 0), V(0, 10, 0)),
                StraightEdge.throughPoints(V(0, 10, 0), V(0, 10, 5)),
                StraightEdge.throughPoints(V(0, 10, 5), V(2, 10, 5)),
                StraightEdge.throughPoints(V(2, 10, 5), V(2, 10, 3))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 1, 0), 10)), [
                StraightEdge.throughPoints(V(2, 10, 5), V(6, 10, 5)),
                StraightEdge.throughPoints(V(6, 10, 5), V(6, 10, 3)),
                StraightEdge.throughPoints(V(6, 10, 3), V(2, 10, 3)),
                StraightEdge.throughPoints(V(2, 10, 3), V(2, 10, 5))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0)), [
                StraightEdge.throughPoints(V(2, 0, 5), V(0, 0, 5)),
                StraightEdge.throughPoints(V(0, 0, 5), V(0, 0, 0)),
                StraightEdge.throughPoints(V(0, 0, 0), V(10, 0, 0)),
                StraightEdge.throughPoints(V(10, 0, 0), V(10, 0, 5)),
                StraightEdge.throughPoints(V(10, 0, 5), V(6, 0, 5)),
                StraightEdge.throughPoints(V(6, 0, 5), V(6, 0, 3)),
                StraightEdge.throughPoints(V(6, 0, 3), V(2, 0, 3)),
                StraightEdge.throughPoints(V(2, 0, 3), V(2, 0, 5))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0)), [
                StraightEdge.throughPoints(V(2, 0, 3), V(6, 0, 3)),
                StraightEdge.throughPoints(V(6, 0, 3), V(6, 0, 5)),
                StraightEdge.throughPoints(V(6, 0, 5), V(2, 0, 5)),
                StraightEdge.throughPoints(V(2, 0, 5), V(2, 0, 3))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 5)), [
                StraightEdge.throughPoints(V(2, 10, 5), V(0, 10, 5)),
                StraightEdge.throughPoints(V(0, 10, 5), V(0, 0, 5)),
                StraightEdge.throughPoints(V(0, 0, 5), V(2, 0, 5)),
                StraightEdge.throughPoints(V(2, 0, 5), V(2, 10, 5))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 5)), [
                StraightEdge.throughPoints(V(2, 0, 5), V(6, 0, 5)),
                StraightEdge.throughPoints(V(6, 0, 5), V(6, 10, 5)),
                StraightEdge.throughPoints(V(6, 10, 5), V(2, 10, 5)),
                StraightEdge.throughPoints(V(2, 10, 5), V(2, 0, 5))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 5)), [
                StraightEdge.throughPoints(V(6, 0, 5), V(10, 0, 5)),
                StraightEdge.throughPoints(V(10, 0, 5), V(10, 10, 5)),
                StraightEdge.throughPoints(V(10, 10, 5), V(6, 10, 5)),
                StraightEdge.throughPoints(V(6, 10, 5), V(6, 0, 5))], []),
            new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), 0)), [
                StraightEdge.throughPoints(V(0, 10, 0), V(0, 0, 0)),
                StraightEdge.throughPoints(V(0, 0, 0), V(0, 0, 5)),
                StraightEdge.throughPoints(V(0, 0, 5), V(0, 10, 5)),
                StraightEdge.throughPoints(V(0, 10, 5), V(0, 10, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 10)), [
                StraightEdge.throughPoints(V(10, 0, 0), V(10, 10, 0)),
                StraightEdge.throughPoints(V(10, 10, 0), V(10, 10, 5)),
                StraightEdge.throughPoints(V(10, 10, 5), V(10, 0, 5)),
                StraightEdge.throughPoints(V(10, 0, 5), V(10, 0, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 0)), [
                StraightEdge.throughPoints(V(0, 0, 0), V(0, 10, 0)),
                StraightEdge.throughPoints(V(0, 10, 0), V(10, 10, 0)),
                StraightEdge.throughPoints(V(10, 10, 0), V(10, 0, 0)),
                StraightEdge.throughPoints(V(10, 0, 0), V(0, 0, 0))], [])], false)
        assert.b2Equal(box, box2, box.plus(box2), expected)
    },

    'B2T.box(10, 10, 10) + B2T.box(10, 12, 12).translate(0, 0, -2) overlapping faces result contains union of faces 3'(assert) {
        const box = B2T.box(10, 10, 10, 'box0'), box2 = B2T.box(10, 12, 12, 'box').translate(0, 0, -2)
        const expected = new B2([
            new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), 0)), [
                StraightEdge.throughPoints(V(0, 10, 10), V(0, 10, 0)),
                StraightEdge.throughPoints(V(0, 10, 0), V(0, 0, 0)),
                StraightEdge.throughPoints(V(0, 0, 0), V(0, 0, 10)),
                StraightEdge.throughPoints(V(0, 0, 10), V(0, 10, 10))], []),
            new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 10)), [
                StraightEdge.throughPoints(V(10, 0, 0), V(10, 10, 0)),
                StraightEdge.throughPoints(V(10, 10, 0), V(10, 10, 10)),
                StraightEdge.throughPoints(V(10, 10, 10), V(10, 0, 10)),
                StraightEdge.throughPoints(V(10, 0, 10), V(10, 0, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0)), [
                StraightEdge.throughPoints(V(0, 0, 0), V(10, 0, 0)),
                StraightEdge.throughPoints(V(10, 0, 0), V(10, 0, 10)),
                StraightEdge.throughPoints(V(10, 0, 10), V(0, 0, 10)),
                StraightEdge.throughPoints(V(0, 0, 10), V(0, 0, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 10)), [
                StraightEdge.throughPoints(V(0, 0, 10), V(10, 0, 10)),
                StraightEdge.throughPoints(V(10, 0, 10), V(10, 10, 10)),
                StraightEdge.throughPoints(V(10, 10, 10), V(0, 10, 10)),
                StraightEdge.throughPoints(V(0, 10, 10), V(0, 0, 10))], []),
            new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), 0)), [
                StraightEdge.throughPoints(V(0, 10, 10), V(0, 12, 10)),
                StraightEdge.throughPoints(V(0, 12, 10), V(0, 12, -2)),
                StraightEdge.throughPoints(V(0, 12, -2), V(0, 0, -2)),
                StraightEdge.throughPoints(V(0, 0, -2), V(0, 0, 0)),
                StraightEdge.throughPoints(V(0, 0, 0), V(0, 10, 0)),
                StraightEdge.throughPoints(V(0, 10, 0), V(0, 10, 10))], []),
            new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 10)), [
                StraightEdge.throughPoints(V(10, 10, 0), V(10, 0, 0)),
                StraightEdge.throughPoints(V(10, 0, 0), V(10, 0, -2)),
                StraightEdge.throughPoints(V(10, 0, -2), V(10, 12, -2)),
                StraightEdge.throughPoints(V(10, 12, -2), V(10, 12, 10)),
                StraightEdge.throughPoints(V(10, 12, 10), V(10, 10, 10)),
                StraightEdge.throughPoints(V(10, 10, 10), V(10, 10, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0)), [
                StraightEdge.throughPoints(V(0, 0, 0), V(0, 0, -2)),
                StraightEdge.throughPoints(V(0, 0, -2), V(10, 0, -2)),
                StraightEdge.throughPoints(V(10, 0, -2), V(10, 0, 0)),
                StraightEdge.throughPoints(V(10, 0, 0), V(0, 0, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 10)), [
                StraightEdge.throughPoints(V(10, 10, 10), V(10, 12, 10)),
                StraightEdge.throughPoints(V(10, 12, 10), V(0, 12, 10)),
                StraightEdge.throughPoints(V(0, 12, 10), V(0, 10, 10)),
                StraightEdge.throughPoints(V(0, 10, 10), V(10, 10, 10))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 1, 0), 12)), [
                StraightEdge.throughPoints(V(10, 12, -2), V(0, 12, -2)),
                StraightEdge.throughPoints(V(0, 12, -2), V(0, 12, 10)),
                StraightEdge.throughPoints(V(0, 12, 10), V(10, 12, 10)),
                StraightEdge.throughPoints(V(10, 12, 10), V(10, 12, -2))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 2)), [
                StraightEdge.throughPoints(V(0, 0, -2), V(0, 12, -2)),
                StraightEdge.throughPoints(V(0, 12, -2), V(10, 12, -2)),
                StraightEdge.throughPoints(V(10, 12, -2), V(10, 0, -2)),
                StraightEdge.throughPoints(V(10, 0, -2), V(0, 0, -2))], [])], false)
        assert.b2Equal(box, box2, box.plus(box2), expected)
    },

    'B2T.box(10, 10, 10) - B2T.box(5, 5, 5)'(assert) {
        const box = B2T.box(10, 10, 10, 'box0'), box2 = B2T.box(1,1,1).scale(8 / 3 ** 0.5).transform(M4.rotationAB(V3.ONES, V3.X))
        const expected = new B2([
            new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), 0)), [
                StraightEdge.throughPoints(V(0, 10, 10), V(0, 10, 0)),
                StraightEdge.throughPoints(V(0, 10, 0), V(0, 0, 0)),
                StraightEdge.throughPoints(V(0, 0, 0), V(0, 0, 10)),
                StraightEdge.throughPoints(V(0, 0, 10), V(0, 10, 10))], []),
            new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 10)), [
                StraightEdge.throughPoints(V(10, 0, 0), V(10, 10, 0)),
                StraightEdge.throughPoints(V(10, 10, 0), V(10, 10, 10)),
                StraightEdge.throughPoints(V(10, 10, 10), V(10, 0, 10)),
                StraightEdge.throughPoints(V(10, 0, 10), V(10, 0, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0)), [
                StraightEdge.throughPoints(V(0, 0, 0), V(10, 0, 0)),
                StraightEdge.throughPoints(V(10, 0, 0), V(10, 0, 10)),
                StraightEdge.throughPoints(V(10, 0, 10), V(0, 0, 10)),
                StraightEdge.throughPoints(V(0, 0, 10), V(0, 0, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 10)), [
                StraightEdge.throughPoints(V(0, 0, 10), V(10, 0, 10)),
                StraightEdge.throughPoints(V(10, 0, 10), V(10, 10, 10)),
                StraightEdge.throughPoints(V(10, 10, 10), V(0, 10, 10)),
                StraightEdge.throughPoints(V(0, 10, 10), V(0, 0, 10))], []),
            new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), 0)), [
                StraightEdge.throughPoints(V(0, 10, 10), V(0, 12, 10)),
                StraightEdge.throughPoints(V(0, 12, 10), V(0, 12, -2)),
                StraightEdge.throughPoints(V(0, 12, -2), V(0, 0, -2)),
                StraightEdge.throughPoints(V(0, 0, -2), V(0, 0, 0)),
                StraightEdge.throughPoints(V(0, 0, 0), V(0, 10, 0)),
                StraightEdge.throughPoints(V(0, 10, 0), V(0, 10, 10))], []),
            new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 10)), [
                StraightEdge.throughPoints(V(10, 10, 0), V(10, 0, 0)),
                StraightEdge.throughPoints(V(10, 0, 0), V(10, 0, -2)),
                StraightEdge.throughPoints(V(10, 0, -2), V(10, 12, -2)),
                StraightEdge.throughPoints(V(10, 12, -2), V(10, 12, 10)),
                StraightEdge.throughPoints(V(10, 12, 10), V(10, 10, 10)),
                StraightEdge.throughPoints(V(10, 10, 10), V(10, 10, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0)), [
                StraightEdge.throughPoints(V(0, 0, 0), V(0, 0, -2)),
                StraightEdge.throughPoints(V(0, 0, -2), V(10, 0, -2)),
                StraightEdge.throughPoints(V(10, 0, -2), V(10, 0, 0)),
                StraightEdge.throughPoints(V(10, 0, 0), V(0, 0, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 10)), [
                StraightEdge.throughPoints(V(10, 10, 10), V(10, 12, 10)),
                StraightEdge.throughPoints(V(10, 12, 10), V(0, 12, 10)),
                StraightEdge.throughPoints(V(0, 12, 10), V(0, 10, 10)),
                StraightEdge.throughPoints(V(0, 10, 10), V(10, 10, 10))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 1, 0), 12)), [
                StraightEdge.throughPoints(V(10, 12, -2), V(0, 12, -2)),
                StraightEdge.throughPoints(V(0, 12, -2), V(0, 12, 10)),
                StraightEdge.throughPoints(V(0, 12, 10), V(10, 12, 10)),
                StraightEdge.throughPoints(V(10, 12, 10), V(10, 12, -2))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 2)), [
                StraightEdge.throughPoints(V(0, 0, -2), V(0, 12, -2)),
                StraightEdge.throughPoints(V(0, 12, -2), V(10, 12, -2)),
                StraightEdge.throughPoints(V(10, 12, -2), V(10, 0, -2)),
                StraightEdge.throughPoints(V(10, 0, -2), V(0, 0, -2))], [])], false)
        assert.b2Equal(box, box2, box.minus(box2), expected)
    },

    'B2.withMergedFaces'(assert) {
        const box = B2T.box(5, 5, 5)
        const boxToMerge = new B2(box.faces.filter((face:PlaneFace) => face.surface.plane.normal.x != 1).concat(
            box.translate(5, 0, 0).faces.filter((face:PlaneFace) => face.surface.plane.normal.x != -1)
        ), false)

        assert.equal(boxToMerge.faces.length, 10)
        const boxMerged = boxToMerge.withMergedFaces()
        assert.b2Equal(boxToMerge, B2.EMPTY, boxMerged, B2T.box(10, 5, 5))
    }
})