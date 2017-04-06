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
            new PlaneFace(new P3(V(0.8320502943378437, -0.5547001962252291, 0), 1.386750490563073), [
                new StraightEdge(new L3(V(5, 5, 5), V(0, 0, -1)), V(5, 5, 5), V(5, 5, -5), 0, 10),
                new StraightEdge(new L3(V(5, 5, -5), V(0.42640143271122083, 0.6396021490668313, 0.6396021490668313)), V(5, 5, -5), V(9, 11, 1), 0, 9.38083151964686),
                new StraightEdge(new L3(V(5, 5, 5), V(0.48507125007266594, 0.7276068751089989, -0.48507125007266594)), V(9, 11, 1), V(5, 5, 5), 8.246211251235321, 0)]),
            new PlaneFace(new P3(V(-0.7071067811865475, -0.7071067811865475, 0), -7.071067811865475), [
                new StraightEdge(new L3(V(5, 5, 5), V(-0.5298129428260175, 0.5298129428260175, -0.6622661785325219)), V(5, 5, 5), V(1, 9, 0), 0, 7.54983443527075),
                new StraightEdge(new L3(V(5, 5, -5), V(-0.5298129428260175, 0.5298129428260175, 0.6622661785325219)), V(1, 9, 0), V(5, 5, -5), 7.54983443527075, 0),
                new StraightEdge(new L3(V(5, 5, 5), V(0, 0, -1)), V(5, 5, -5), V(5, 5, 5), 10, 0)]),
            new PlaneFace(new P3(V(-0.1003911722115382, 0.7362019295512802, -0.6692744814102547), 6.525426193749983), [
                new StraightEdge(new L3(V(5, 5, -5), V(-0.5298129428260175, 0.5298129428260175, 0.6622661785325219)), V(5, 5, -5), V(1, 9, 0), 0, 7.54983443527075),
                new StraightEdge(new L3(V(9, 11, 1), V(-0.9630868246861536, -0.2407717061715384, -0.1203858530857692)), V(1, 9, 0), V(9, 11, 1), 8.306623862918075, 0),
                new StraightEdge(new L3(V(5, 5, -5), V(0.42640143271122083, 0.6396021490668313, 0.6396021490668313)), V(9, 11, 1), V(5, 5, -5), 9.38083151964686, 0)]),
            new PlaneFace(new P3(V(-0.25177250044296223, 0.6474150011390457, 0.7193500012656064), 5.57496250980845), [
                new StraightEdge(new L3(V(9, 11, 1), V(-0.9630868246861536, -0.2407717061715384, -0.1203858530857692)), V(9, 11, 1), V(1, 9, 0), 0, 8.306623862918075),
                new StraightEdge(new L3(V(5, 5, 5), V(-0.5298129428260175, 0.5298129428260175, -0.6622661785325219)), V(1, 9, 0), V(5, 5, 5), 7.54983443527075, 0),
                new StraightEdge(new L3(V(5, 5, 5), V(0.48507125007266594, 0.7276068751089989, -0.48507125007266594)), V(5, 5, 5), V(9, 11, 1), 0, 8.246211251235321)])], false)
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

    'splitsVolumeEnclosingFacesP '(assert) {
        const brep = B2T.box().flipped()
        brep.buildAdjacencies()
        const edge = (a, b) => brep.faces.map(face => face.getAllEdges()).concatenated().find(edge => edge.a.like(a) && edge.b.like(b)).getCanon()
        assert.equal(splitsVolumeEnclosingFacesP(brep, edge(V(0, 0, 0), V(1, 0, 0)), V(0.5, 0, 0), V(0, 0, -1), V3.O), INSIDE)
        assert.equal(splitsVolumeEnclosingFacesP(brep, edge(V(0, 0, 1), V(1, 0, 1)), V(0.5, 0, 1), V(0, 0, 1), V3.O), INSIDE)


        const brep2 = B2T.box(12, 2, 3, "").translate(-6, 0, 1).flipped()
        brep2.buildAdjacencies()
        assert.equal(splitsVolumeEnclosingFacesP(
            brep2,
            new StraightEdge(new L3(V(6, 0, 1), V(-1, 0, 0)), V(6, 0, 1), V(-6, 0, 1), 0, 12),
            V(-4.898979485566356, 5.999519546087386e-16, 1.000000000000001),
            V(0, 1.2246467991473547e-16, -4.898979485566356),
            V(-0.9797958971132712, 1.199903909217477e-16, 0.20000000000000023)), INSIDE)
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

        assert.equal(splitsVolumeEnclosingCone(brep, V3.O, V(1, 1, 1)), INSIDE)
        assert.equal(splitsVolumeEnclosingCone(brep, V3.O, V(-1, 1, 1)), OUTSIDE)
        assert.equal(splitsVolumeEnclosingCone(brep, V3.O, V3.X), ALONG_EDGE_OR_PLANE)
        assert.equal(splitsVolumeEnclosingCone(brep, V3.O, V3.X.negated()), OUTSIDE)
        assert.equal(splitsVolumeEnclosingCone(brep, V(0, 5, 5), V3.Y), OUTSIDE)
        assert.equal(splitsVolumeEnclosingCone(brep, V3.O, V(1, 1, 0)), ALONG_EDGE_OR_PLANE)
    },
    'splitsVolumeEnclosingCone2'(assert) {
        const brep = B2T.box(5, 5, 5)

        assert.equal(splitsVolumeEnclosingCone2(brep, V3.O, new L3(V3.O, V(1, 1, 1).unit()), 0, 1), INSIDE)
        assert.equal(splitsVolumeEnclosingCone2(brep, V3.O, new L3(V3.O, V(-1, 1, 1).unit()), 0, 1), OUTSIDE)
        assert.equal(splitsVolumeEnclosingCone2(brep, V3.O, new L3(V3.O, V3.X), 0, -1), OUTSIDE)
        assert.equal(splitsVolumeEnclosingCone2(brep, V(0, 5, 5), new L3(V3.O, V3.Y), 0, 1), OUTSIDE)
        assert.equal(splitsVolumeEnclosingCone2(brep, V3.O, new L3(V3.O, V3.X), 0, 1), ALONG_EDGE_OR_PLANE)
        assert.equal(splitsVolumeEnclosingCone2(brep, V3.O, new L3(V3.O, V(1, 1, 0).unit()), 0, 1), ALONG_EDGE_OR_PLANE)
    },

    'pointsToInside'(assert) {
        const face = PlaneFace.forVertices(P3.XY, [V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)])

        assert.equal(face.pointsToInside(V3.O, V3.X.negated()), PointVsFace.OUTSIDE)

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
    'pointsToInside2'(assert) {
        const face = PlaneFace.forVertices(P3.XY, [V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)])

        assert.equal(face.pointsToInside2(V3.O, V3.X.negated()), PointVsFace.OUTSIDE)

        const v = V(10, 10, 0)
        assert.equal(face.pointsToInside2(v, V(-1, -1, 0)), PointVsFace.INSIDE)
        assert.equal(face.pointsToInside2(v, V(1, 1, 0)), PointVsFace.OUTSIDE)
        assert.equal(face.pointsToInside2(v, V(-1, 0, 0)), PointVsFace.ON_EDGE)
        assert.equal(face.pointsToInside2(v, V(0, -1, 0)), PointVsFace.ON_EDGE)

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
        let result = planeFaceEdgeISPsWithPlane(brep.faces[2], L3.Y.translate(0, 0, -1), P3.XY.translate(0, 0, -1)).map(is => is.p)
        assert.V3ArraysLike(result, [V(0, 0, -1), V(0, 10, -1)])
        result = planeFaceEdgeISPsWithPlane(brep.faces[2], L3.Y, P3.XY).map(is => is.p)
        assert.V3ArraysLike(result, [V(0, 0, 0), V(0, 10, 0)])
        result = planeFaceEdgeISPsWithPlane(brep.translate(0, 0, 10).faces[2], L3.Y.translate(0, 0, 6), P3.XY.translate(0, 0, 6)).map(is => is.p)
        assert.V3ArraysLike(result, [V(0, 0, 6), V(0, 10, 6)])
    },

    'B2.prototype.minus B2T.box(5, 5, 5).minus(B2T.box(1, 1, 6))'(assert) {
        const a = B2T.box(5, 5, 5, 'a')
        const b = B2T.box(1, 1, 6, 'b')
        const result = new B2([
            new PlaneFace(new P3(V(-1, 0, 0), 0), [
                new StraightEdge(new L3(V(0, 1, 0), V(0, 0, 1)), V(0, 1, 0), V(0, 1, 5), 0, 5),
                new StraightEdge(new L3(V(0, 0, 5), V(0, 1, 0)), V(0, 1, 5), V(0, 5, 5), 1, 5),
                new StraightEdge(new L3(V(0, 5, 0), V(0, 0, 1)), V(0, 5, 5), V(0, 5, 0), 5, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 5, 0), V(0, 1, 0), 5, 1)]),
            new PlaneFace(new P3(V(0, -1, 0), 0), [
                new StraightEdge(new L3(V(1, 0, 0), V(0, 0, -1)), V(1, 0, 5), V(1, 0, 0), -5, 0),
                new StraightEdge(new L3(V(5, 0, 0), V(-1, 0, 0)), V(1, 0, 0), V(5, 0, 0), 4, 0),
                new StraightEdge(new L3(V(5, 0, 0), V(0, 0, 1)), V(5, 0, 0), V(5, 0, 5), 0, 5),
                new StraightEdge(new L3(V(5, 0, 5), V(-1, 0, 0)), V(5, 0, 5), V(1, 0, 5), 0, 4)]),
            new PlaneFace(new P3(V(0, 0, -1), 0), [
                new StraightEdge(new L3(V(0, 1, 0), V(-1, 0, 0)), V(1, 1, 0), V(0, 1, 0), -1, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 1, 0), V(0, 5, 0), 1, 5),
                new StraightEdge(new L3(V(0, 5, 0), V(1, 0, 0)), V(0, 5, 0), V(5, 5, 0), 0, 5),
                new StraightEdge(new L3(V(5, 5, 0), V(0, -1, 0)), V(5, 5, 0), V(5, 0, 0), 0, 5),
                new StraightEdge(new L3(V(5, 0, 0), V(-1, 0, 0)), V(5, 0, 0), V(1, 0, 0), 0, 4),
                new StraightEdge(new L3(V(1, 0, 0), V(0, 1, 0)), V(1, 0, 0), V(1, 1, 0), 0, 1)]),
            new PlaneFace(new P3(V(0, 0, 1), 5), [
                new StraightEdge(new L3(V(0, 1, 5), V(1, 0, 0)), V(0, 1, 5), V(1, 1, 5), 0, 1),
                new StraightEdge(new L3(V(1, 0, 5), V(0, -1, 0)), V(1, 1, 5), V(1, 0, 5), -1, 0),
                new StraightEdge(new L3(V(5, 0, 5), V(-1, 0, 0)), V(1, 0, 5), V(5, 0, 5), 4, 0),
                new StraightEdge(new L3(V(5, 5, 5), V(0, -1, 0)), V(5, 0, 5), V(5, 5, 5), 5, 0),
                new StraightEdge(new L3(V(0, 5, 5), V(1, 0, 0)), V(5, 5, 5), V(0, 5, 5), 5, 0),
                new StraightEdge(new L3(V(0, 0, 5), V(0, 1, 0)), V(0, 5, 5), V(0, 1, 5), 5, 1)]),
            new PlaneFace(new P3(V(0, 1, 0), 5), [
                new StraightEdge(new L3(V(0, 5, 0), V(1, 0, 0)), V(5, 5, 0), V(0, 5, 0), 5, 0),
                new StraightEdge(new L3(V(0, 5, 0), V(0, 0, 1)), V(0, 5, 0), V(0, 5, 5), 0, 5),
                new StraightEdge(new L3(V(0, 5, 5), V(1, 0, 0)), V(0, 5, 5), V(5, 5, 5), 0, 5),
                new StraightEdge(new L3(V(5, 5, 0), V(0, 0, 1)), V(5, 5, 5), V(5, 5, 0), 5, 0)]),
            new PlaneFace(new P3(V(1, 0, 0), 5), [
                new StraightEdge(new L3(V(5, 5, 0), V(0, -1, 0)), V(5, 0, 0), V(5, 5, 0), 5, 0),
                new StraightEdge(new L3(V(5, 5, 0), V(0, 0, 1)), V(5, 5, 0), V(5, 5, 5), 0, 5),
                new StraightEdge(new L3(V(5, 5, 5), V(0, -1, 0)), V(5, 5, 5), V(5, 0, 5), 0, 5),
                new StraightEdge(new L3(V(5, 0, 0), V(0, 0, 1)), V(5, 0, 5), V(5, 0, 0), 5, 0)]),
            new PlaneFace(new P3(V(0, -1, 0), -1), [
                new StraightEdge(new L3(V(0, 1, 5), V(1, 0, 0)), V(1, 1, 5), V(0, 1, 5), 1, 0),
                new StraightEdge(new L3(V(0, 1, 0), V(0, 0, 1)), V(0, 1, 5), V(0, 1, 0), 5, 0),
                new StraightEdge(new L3(V(0, 1, 0), V(-1, 0, 0)), V(0, 1, 0), V(1, 1, 0), 0, -1),
                new StraightEdge(new L3(V(1, 1, 0), V(0, 0, 1)), V(1, 1, 0), V(1, 1, 5), 0, 5)]),
            new PlaneFace(new P3(V(-1, 0, 0), -1), [
                new StraightEdge(new L3(V(1, 0, 5), V(0, -1, 0)), V(1, 0, 5), V(1, 1, 5), 0, -1),
                new StraightEdge(new L3(V(1, 1, 0), V(0, 0, 1)), V(1, 1, 5), V(1, 1, 0), 5, 0),
                new StraightEdge(new L3(V(1, 0, 0), V(0, 1, 0)), V(1, 1, 0), V(1, 0, 0), 1, 0),
                new StraightEdge(new L3(V(1, 0, 0), V(0, 0, -1)), V(1, 0, 0), V(1, 0, 5), 0, -5)])], false)
        b2Equal(assert, a, b, a.minus(b), result)
    },
    'B2.prototype.minus B2T.box(5, 5, 5).minus(B2T.box(1, 1, 5))'(assert) {
        const a = B2T.box(5, 5, 5)
        const b = B2T.box(1, 1, 5)
        const result = new B2([
            new PlaneFace(new P3(V(-1, 0, 0), 0), [
                new StraightEdge(new L3(V(0, 1, 0), V(0, 0, 1)), V(0, 1, 0), V(0, 1, 5), 0, 5),
                new StraightEdge(new L3(V(0, 0, 5), V(0, 1, 0)), V(0, 1, 5), V(0, 5, 5), 1, 5),
                new StraightEdge(new L3(V(0, 5, 0), V(0, 0, 1)), V(0, 5, 5), V(0, 5, 0), 5, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 5, 0), V(0, 1, 0), 5, 1)]),
            new PlaneFace(new P3(V(0, -1, 0), 0), [
                new StraightEdge(new L3(V(1, 0, 0), V(0, 0, -1)), V(1, 0, 5), V(1, 0, 0), -5, 0),
                new StraightEdge(new L3(V(5, 0, 0), V(-1, 0, 0)), V(1, 0, 0), V(5, 0, 0), 4, 0),
                new StraightEdge(new L3(V(5, 0, 0), V(0, 0, 1)), V(5, 0, 0), V(5, 0, 5), 0, 5),
                new StraightEdge(new L3(V(5, 0, 5), V(-1, 0, 0)), V(5, 0, 5), V(1, 0, 5), 0, 4)]),
            new PlaneFace(new P3(V(0, 0, -1), 0), [
                new StraightEdge(new L3(V(0, 1, 0), V(-1, 0, 0)), V(1, 1, 0), V(0, 1, 0), -1, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 1, 0), V(0, 5, 0), 1, 5),
                new StraightEdge(new L3(V(0, 5, 0), V(1, 0, 0)), V(0, 5, 0), V(5, 5, 0), 0, 5),
                new StraightEdge(new L3(V(5, 5, 0), V(0, -1, 0)), V(5, 5, 0), V(5, 0, 0), 0, 5),
                new StraightEdge(new L3(V(5, 0, 0), V(-1, 0, 0)), V(5, 0, 0), V(1, 0, 0), 0, 4),
                new StraightEdge(new L3(V(1, 0, 0), V(0, 1, 0)), V(1, 0, 0), V(1, 1, 0), 0, 1)]),
            new PlaneFace(new P3(V(0, 0, 1), 5), [
                new StraightEdge(new L3(V(0, 1, 5), V(1, 0, 0)), V(0, 1, 5), V(1, 1, 5), 0, 1),
                new StraightEdge(new L3(V(1, 0, 5), V(0, -1, 0)), V(1, 1, 5), V(1, 0, 5), -1, 0),
                new StraightEdge(new L3(V(5, 0, 5), V(-1, 0, 0)), V(1, 0, 5), V(5, 0, 5), 4, 0),
                new StraightEdge(new L3(V(5, 5, 5), V(0, -1, 0)), V(5, 0, 5), V(5, 5, 5), 5, 0),
                new StraightEdge(new L3(V(0, 5, 5), V(1, 0, 0)), V(5, 5, 5), V(0, 5, 5), 5, 0),
                new StraightEdge(new L3(V(0, 0, 5), V(0, 1, 0)), V(0, 5, 5), V(0, 1, 5), 5, 1)]),
            new PlaneFace(new P3(V(0, 1, 0), 5), [
                new StraightEdge(new L3(V(0, 5, 0), V(1, 0, 0)), V(5, 5, 0), V(0, 5, 0), 5, 0),
                new StraightEdge(new L3(V(0, 5, 0), V(0, 0, 1)), V(0, 5, 0), V(0, 5, 5), 0, 5),
                new StraightEdge(new L3(V(0, 5, 5), V(1, 0, 0)), V(0, 5, 5), V(5, 5, 5), 0, 5),
                new StraightEdge(new L3(V(5, 5, 0), V(0, 0, 1)), V(5, 5, 5), V(5, 5, 0), 5, 0)]),
            new PlaneFace(new P3(V(1, 0, 0), 5), [
                new StraightEdge(new L3(V(5, 5, 0), V(0, -1, 0)), V(5, 0, 0), V(5, 5, 0), 5, 0),
                new StraightEdge(new L3(V(5, 5, 0), V(0, 0, 1)), V(5, 5, 0), V(5, 5, 5), 0, 5),
                new StraightEdge(new L3(V(5, 5, 5), V(0, -1, 0)), V(5, 5, 5), V(5, 0, 5), 0, 5),
                new StraightEdge(new L3(V(5, 0, 0), V(0, 0, 1)), V(5, 0, 5), V(5, 0, 0), 5, 0)]),
            new PlaneFace(new P3(V(0, -1, 0), -1), [
                new StraightEdge(new L3(V(0, 1, 5), V(1, 0, 0)), V(1, 1, 5), V(0, 1, 5), 1, 0),
                new StraightEdge(new L3(V(0, 1, 0), V(0, 0, 1)), V(0, 1, 5), V(0, 1, 0), 5, 0),
                new StraightEdge(new L3(V(0, 1, 0), V(-1, 0, 0)), V(0, 1, 0), V(1, 1, 0), 0, -1),
                new StraightEdge(new L3(V(1, 1, 0), V(0, 0, 1)), V(1, 1, 0), V(1, 1, 5), 0, 5)]),
            new PlaneFace(new P3(V(-1, 0, 0), -1), [
                new StraightEdge(new L3(V(1, 0, 5), V(0, -1, 0)), V(1, 0, 5), V(1, 1, 5), 0, -1),
                new StraightEdge(new L3(V(1, 1, 0), V(0, 0, 1)), V(1, 1, 5), V(1, 1, 0), 5, 0),
                new StraightEdge(new L3(V(1, 0, 0), V(0, 1, 0)), V(1, 1, 0), V(1, 0, 0), 1, 0),
                new StraightEdge(new L3(V(1, 0, 0), V(0, 0, -1)), V(1, 0, 0), V(1, 0, 5), 0, -5)])], false)
        b2Equal(assert, a, b, a.minus(b), result)
    },
    'B2 edge/face intersection'(assert) {
        const wideBox = B2T.box(10, 10, 5)
        const extrusion = B2T.extrudeVertices([V(1, 0), V(0, 3), V(-2, 5), V(5, 5)], P3.XY.flipped(), V(0, 0, 10), 'lol')
            .translate(0, 1, 1)
        const result = new B2([
            new PlaneFace(new P3(V(-1, 0, 0), 0), [
                new StraightEdge(new L3(V(0, 0, 1), V(0, -1, 0)), V(0, 4, 1), V(0, 6, 1), -4, -6),
                new StraightEdge(new L3(V(0, 6, 0), V(0, 0, -1)), V(0, 6, 1), V(0, 6, 5), -1, -5),
                new StraightEdge(new L3(V(0, 0, 5), V(0, 1, 0)), V(0, 6, 5), V(0, 10, 5), 6, 10),
                new StraightEdge(new L3(V(0, 10, 0), V(0, 0, 1)), V(0, 10, 5), V(0, 10, 0), 5, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 10, 0), V(0, 0, 0), 10, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 0), V(0, 0, 5), 0, 5),
                new StraightEdge(new L3(V(0, 0, 5), V(0, 1, 0)), V(0, 0, 5), V(0, 4, 5), 0, 4),
                new StraightEdge(new L3(V(0, 4, 0), V(0, 0, 1)), V(0, 4, 5), V(0, 4, 1), 5, 1)]),
            new PlaneFace(new P3(V(0, 0, 1), 5), [
                new StraightEdge(new L3(V(0, 0, 5), V(0, 1, 0)), V(0, 4, 5), V(0, 0, 5), 4, 0),
                new StraightEdge(new L3(V(10, 0, 5), V(-1, 0, 0)), V(0, 0, 5), V(10, 0, 5), 10, 0),
                new StraightEdge(new L3(V(10, 10, 5), V(0, -1, 0)), V(10, 0, 5), V(10, 10, 5), 10, 0),
                new StraightEdge(new L3(V(0, 10, 5), V(1, 0, 0)), V(10, 10, 5), V(0, 10, 5), 10, 0),
                new StraightEdge(new L3(V(0, 0, 5), V(0, 1, 0)), V(0, 10, 5), V(0, 6, 5), 10, 6),
                new StraightEdge(new L3(V(0, 6, 5), V(-1, 0, 0)), V(0, 6, 5), V(5, 6, 5), 0, -5),
                new StraightEdge(new L3(V(0.12195121951219434, -0.09756097560975575, 5), V(0.6246950475544243, 0.7808688094430303, 0)), V(5, 6, 5), V(1, 1, 5), 7.808688094430304, 1.4055638569974547),
                new StraightEdge(new L3(V(1.2, 0.4, 5), V(0.316227766016838, -0.9486832980505138, 0)), V(1, 1, 5), V(0, 4, 5), -0.6324555320336758, -3.7947331922020537)]),
            new PlaneFace(new P3(V(0, 1, 0), 10), [
                new StraightEdge(new L3(V(0, 10, 0), V(1, 0, 0)), V(10, 10, 0), V(0, 10, 0), 10, 0),
                new StraightEdge(new L3(V(0, 10, 0), V(0, 0, 1)), V(0, 10, 0), V(0, 10, 5), 0, 5),
                new StraightEdge(new L3(V(0, 10, 5), V(1, 0, 0)), V(0, 10, 5), V(10, 10, 5), 0, 10),
                new StraightEdge(new L3(V(10, 10, 0), V(0, 0, 1)), V(10, 10, 5), V(10, 10, 0), 5, 0)]),
            new PlaneFace(new P3(V(1, 0, 0), 10), [
                new StraightEdge(new L3(V(10, 10, 0), V(0, -1, 0)), V(10, 0, 0), V(10, 10, 0), 10, 0),
                new StraightEdge(new L3(V(10, 10, 0), V(0, 0, 1)), V(10, 10, 0), V(10, 10, 5), 0, 5),
                new StraightEdge(new L3(V(10, 10, 5), V(0, -1, 0)), V(10, 10, 5), V(10, 0, 5), 0, 10),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 0, 1)), V(10, 0, 5), V(10, 0, 0), 5, 0)]),
            new PlaneFace(new P3(V(0, -1, 0), 0), [
                new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(0, 0, 0), V(10, 0, 0), 10, 0),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 0, 1)), V(10, 0, 0), V(10, 0, 5), 0, 5),
                new StraightEdge(new L3(V(10, 0, 5), V(-1, 0, 0)), V(10, 0, 5), V(0, 0, 5), 0, 10),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 5), V(0, 0, 0), 5, 0)]),
            new PlaneFace(new P3(V(0, 0, -1), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 0, 0), V(0, 10, 0), 0, 10),
                new StraightEdge(new L3(V(0, 10, 0), V(1, 0, 0)), V(0, 10, 0), V(10, 10, 0), 0, 10),
                new StraightEdge(new L3(V(10, 10, 0), V(0, -1, 0)), V(10, 10, 0), V(10, 0, 0), 0, 10),
                new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(10, 0, 0), V(0, 0, 0), 0, 10)]),
            new PlaneFace(new P3(V(-0.9486832980505142, -0.3162277660168372, 0), -1.2649110640673595), [
                new StraightEdge(new L3(V(1, 1, 1), V(0, 0, 1)), V(1, 1, 5), V(1, 1, 11), 4, 10),
                new StraightEdge(new L3(V(1, 1, 11), V(-0.31622776601683794, 0.9486832980505138, 0)), V(1, 1, 11), V(0, 4, 11), 0, 3.1622776601683795),
                new StraightEdge(new L3(V(0, 4, 1), V(0, 0, 1)), V(0, 4, 11), V(0, 4, 5), 10, 4),
                new StraightEdge(new L3(V(1.2, 0.4, 5), V(0.316227766016838, -0.9486832980505138, 0)), V(0, 4, 5), V(1, 1, 5), -3.7947331922020537, -0.6324555320336758)]),
            new PlaneFace(new P3(V(-0.7071067811865485, -0.7071067811865467, 0), -2.8284271247461907), [
                new StraightEdge(new L3(V(0, 4, 11), V(-0.7071067811865475, 0.7071067811865475, 0)), V(0, 4, 11), V(-2, 6, 11), 0, 2.8284271247461903),
                new StraightEdge(new L3(V(-2, 6, 1), V(0, 0, 1)), V(-2, 6, 11), V(-2, 6, 1), 10, 0),
                new StraightEdge(new L3(V(0, 4, 1), V(-0.7071067811865475, 0.7071067811865475, 0)), V(-2, 6, 1), V(0, 4, 1), 2.8284271247461903, 0),
                new StraightEdge(new L3(V(0, 4, 0), V(0, 0, 1)), V(0, 4, 1), V(0, 4, 5), 1, 5),
                new StraightEdge(new L3(V(0, 4, 1), V(0, 0, 1)), V(0, 4, 5), V(0, 4, 11), 4, 10)]),
            new PlaneFace(new P3(V(0, 1, 0), 6), [
                new StraightEdge(new L3(V(-2, 6, 1), V(1, 0, 0)), V(0, 6, 1), V(-2, 6, 1), 2, 0),
                new StraightEdge(new L3(V(-2, 6, 1), V(0, 0, 1)), V(-2, 6, 1), V(-2, 6, 11), 0, 10),
                new StraightEdge(new L3(V(-2, 6, 11), V(1, 0, 0)), V(-2, 6, 11), V(5, 6, 11), 0, 7),
                new StraightEdge(new L3(V(5, 6, 1), V(0, 0, 1)), V(5, 6, 11), V(5, 6, 5), 10, 4),
                new StraightEdge(new L3(V(0, 6, 5), V(-1, 0, 0)), V(5, 6, 5), V(0, 6, 5), -5, 0),
                new StraightEdge(new L3(V(0, 6, 0), V(0, 0, -1)), V(0, 6, 5), V(0, 6, 1), -5, -1)]),
            new PlaneFace(new P3(V(0.7808688094430296, -0.6246950475544252, 0), 0.15617376188861254), [
                new StraightEdge(new L3(V(5, 6, 1), V(0, 0, 1)), V(5, 6, 5), V(5, 6, 11), 4, 10),
                new StraightEdge(new L3(V(5, 6, 11), V(-0.6246950475544243, -0.7808688094430303, 0)), V(5, 6, 11), V(1, 1, 11), 0, 6.403124237432849),
                new StraightEdge(new L3(V(1, 1, 1), V(0, 0, 1)), V(1, 1, 11), V(1, 1, 5), 10, 4),
                new StraightEdge(new L3(V(0.12195121951219434, -0.09756097560975575, 5), V(0.6246950475544243, 0.7808688094430303, 0)), V(1, 1, 5), V(5, 6, 5), 1.4055638569974547, 7.808688094430304)]),
            new PlaneFace(new P3(V(0, 0, -1), -1), [
                new StraightEdge(new L3(V(0, 4, 1), V(-0.7071067811865475, 0.7071067811865475, 0)), V(0, 4, 1), V(-2, 6, 1), 0, 2.8284271247461903),
                new StraightEdge(new L3(V(-2, 6, 1), V(1, 0, 0)), V(-2, 6, 1), V(0, 6, 1), 0, 2),
                new StraightEdge(new L3(V(0, 0, 1), V(0, -1, 0)), V(0, 6, 1), V(0, 4, 1), -6, -4)]),
            new PlaneFace(new P3(V(0, 0, 1), 11), [
                new StraightEdge(new L3(V(5, 6, 11), V(-0.6246950475544243, -0.7808688094430303, 0)), V(1, 1, 11), V(5, 6, 11), 6.403124237432849, 0),
                new StraightEdge(new L3(V(-2, 6, 11), V(1, 0, 0)), V(5, 6, 11), V(-2, 6, 11), 7, 0),
                new StraightEdge(new L3(V(0, 4, 11), V(-0.7071067811865475, 0.7071067811865475, 0)), V(-2, 6, 11), V(0, 4, 11), 2.8284271247461903, 0),
                new StraightEdge(new L3(V(1, 1, 11), V(-0.31622776601683794, 0.9486832980505138, 0)), V(0, 4, 11), V(1, 1, 11), 3.1622776601683795, 0)])], false)
        b2Equal(assert, wideBox, extrusion, wideBox.plus(extrusion), result)
    },
    'B2T.box(10, 10, 5) + B2T.box(10, 10, 5).translate(0, 0, 5) touching coplanar faces result contains both'(assert) {
        const a = B2T.box(10, 10, 5), b = a.translate(0, 0, 5)

        // the expected result is the union of a and b's faces, minus the ones where they touch, which we remove with an AABB test
        const testAABB = AABB.forXYZ(12, 12, 2).translate(-1, -1, 4)
        const result = new B2(a.faces.concat(b.faces).filter(face => !testAABB.containsAABB(face.getAABB())), false)
        b2Equal(assert, a, b, b.plus(a), result)
    },
    'B2T.box(10, 10, 5) - B2T.box(10, 5, 5).translate(0, 0, 5) - B2T.box(5, 5, 5).translate(5, 5, 5)'(assert) {
        const box = B2T.box(10, 10, 10), box2 = B2T.box(10, 5, 5).translate(0, 0, 5), box3 = B2T.box(5, 5, 5).translate(5, 5, 5)
        const actual = box.minus(box2).minus(box3)
        const result = new B2([
            new PlaneFace(new P3(V(1, 0, 0), 10), [
                new StraightEdge(new L3(V(10, 0, 5), V(0, -1, 0)), V(10, 10, 5), V(10, 5, 5), -10, -5),
                new StraightEdge(new L3(V(10, 0, 5), V(0, -1, 0)), V(10, 5, 5), V(10, 0, 5), -5, 0),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 0, 1)), V(10, 0, 5), V(10, 0, 0), 5, 0),
                new StraightEdge(new L3(V(10, 10, 0), V(0, -1, 0)), V(10, 0, 0), V(10, 10, 0), 10, 0),
                new StraightEdge(new L3(V(10, 10, 0), V(0, 0, 1)), V(10, 10, 0), V(10, 10, 5), 0, 5)]),
            new PlaneFace(new P3(V(0, 0, 1), 10), [
                new StraightEdge(new L3(V(5, 0, 10), V(0, 1, 0)), V(5, 5, 10), V(5, 10, 10), 5, 10),
                new StraightEdge(new L3(V(0, 10, 10), V(1, 0, 0)), V(5, 10, 10), V(0, 10, 10), 5, 0),
                new StraightEdge(new L3(V(0, 0, 10), V(0, 1, 0)), V(0, 10, 10), V(0, 5, 10), 10, 5),
                new StraightEdge(new L3(V(0, 5, 10), V(1, 0, 0)), V(0, 5, 10), V(5, 5, 10), 0, 5)]),
            new PlaneFace(new P3(V(0, 1, 0), 10), [
                new StraightEdge(new L3(V(5, 10, 0), V(0, 0, -1)), V(5, 10, 10), V(5, 10, 5), -10, -5),
                new StraightEdge(new L3(V(0, 10, 5), V(1, 0, 0)), V(5, 10, 5), V(10, 10, 5), 5, 10),
                new StraightEdge(new L3(V(10, 10, 0), V(0, 0, 1)), V(10, 10, 5), V(10, 10, 0), 5, 0),
                new StraightEdge(new L3(V(0, 10, 0), V(1, 0, 0)), V(10, 10, 0), V(0, 10, 0), 10, 0),
                new StraightEdge(new L3(V(0, 10, 0), V(0, 0, 1)), V(0, 10, 0), V(0, 10, 10), 0, 10),
                new StraightEdge(new L3(V(0, 10, 10), V(1, 0, 0)), V(0, 10, 10), V(5, 10, 10), 0, 5)]),
            new PlaneFace(new P3(V(0, -1, 0), -5), [
                new StraightEdge(new L3(V(5, 5, 0), V(0, 0, 1)), V(5, 5, 5), V(5, 5, 10), 5, 10),
                new StraightEdge(new L3(V(0, 5, 10), V(1, 0, 0)), V(5, 5, 10), V(0, 5, 10), 5, 0),
                new StraightEdge(new L3(V(0, 5, 0), V(0, 0, 1)), V(0, 5, 10), V(0, 5, 5), 10, 5),
                new StraightEdge(new L3(V(0, 5, 5), V(1, 0, 0)), V(0, 5, 5), V(5, 5, 5), 0, 5)]),
            new PlaneFace(new P3(V(0, 0, 1), 5), [
                new StraightEdge(new L3(V(0, 5, 5), V(-1, 0, 0)), V(10, 5, 5), V(5, 5, 5), -10, -5),
                new StraightEdge(new L3(V(0, 5, 5), V(1, 0, 0)), V(5, 5, 5), V(0, 5, 5), 5, 0),
                new StraightEdge(new L3(V(0, 0, 5), V(0, 1, 0)), V(0, 5, 5), V(0, 0, 5), 5, 0),
                new StraightEdge(new L3(V(0, 0, 5), V(-1, 0, 0)), V(0, 0, 5), V(10, 0, 5), 0, -10),
                new StraightEdge(new L3(V(10, 0, 5), V(0, -1, 0)), V(10, 0, 5), V(10, 5, 5), 0, -5)]),
            new PlaneFace(new P3(V(-1, 0, 0), 0), [
                new StraightEdge(new L3(V(0, 5, 0), V(0, 0, 1)), V(0, 5, 5), V(0, 5, 10), 5, 10),
                new StraightEdge(new L3(V(0, 0, 10), V(0, 1, 0)), V(0, 5, 10), V(0, 10, 10), 5, 10),
                new StraightEdge(new L3(V(0, 10, 0), V(0, 0, 1)), V(0, 10, 10), V(0, 10, 0), 10, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 10, 0), V(0, 0, 0), 10, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 0), V(0, 0, 5), 0, 5),
                new StraightEdge(new L3(V(0, 0, 5), V(0, 1, 0)), V(0, 0, 5), V(0, 5, 5), 0, 5)]),
            new PlaneFace(new P3(V(0, -1, 0), 0), [
                new StraightEdge(new L3(V(0, 0, 5), V(-1, 0, 0)), V(10, 0, 5), V(0, 0, 5), -10, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 5), V(0, 0, 0), 5, 0),
                new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(0, 0, 0), V(10, 0, 0), 10, 0),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 0, 1)), V(10, 0, 0), V(10, 0, 5), 0, 5)]),
            new PlaneFace(new P3(V(0, 0, -1), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 0, 0), V(0, 10, 0), 0, 10),
                new StraightEdge(new L3(V(0, 10, 0), V(1, 0, 0)), V(0, 10, 0), V(10, 10, 0), 0, 10),
                new StraightEdge(new L3(V(10, 10, 0), V(0, -1, 0)), V(10, 10, 0), V(10, 0, 0), 0, 10),
                new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(10, 0, 0), V(0, 0, 0), 0, 10)]),
            new PlaneFace(new P3(V(1, 0, 0), 5), [
                new StraightEdge(new L3(V(5, 10, 0), V(0, 0, -1)), V(5, 10, 5), V(5, 10, 10), -5, -10),
                new StraightEdge(new L3(V(5, 0, 10), V(0, 1, 0)), V(5, 10, 10), V(5, 5, 10), 10, 5),
                new StraightEdge(new L3(V(5, 5, 0), V(0, 0, 1)), V(5, 5, 10), V(5, 5, 5), 10, 5),
                new StraightEdge(new L3(V(5, 5, 5), V(0, 1, 0)), V(5, 5, 5), V(5, 10, 5), 0, 5)]),
            new PlaneFace(new P3(V(0, 0, 1), 5), [
                new StraightEdge(new L3(V(0, 5, 5), V(-1, 0, 0)), V(5, 5, 5), V(10, 5, 5), -5, -10),
                new StraightEdge(new L3(V(10, 0, 5), V(0, -1, 0)), V(10, 5, 5), V(10, 10, 5), -5, -10),
                new StraightEdge(new L3(V(0, 10, 5), V(1, 0, 0)), V(10, 10, 5), V(5, 10, 5), 10, 5),
                new StraightEdge(new L3(V(5, 5, 5), V(0, 1, 0)), V(5, 10, 5), V(5, 5, 5), 5, 0)])], false)
        b2Equal(assert, box.minus(box2), box3, actual, result)
    },

    'B2T.box(10, 10, 5) && B2T.box(10, 10, 5).translate(3, 3, 0) overlapping faces result contains intersection of faces'(assert) {
        const box = B2T.box(10, 10, 5, 'a')
        const box2 = B2T.box(10, 10, 4, 'b').translate(3, 3, 0)
        b2Equal(assert, box, box2, box.intersection(box2, true, true), B2T.box(7, 7, 4).translate(3, 3, 0))
    },

    'B2T.box(10, 10, 5) + B2T.box(10, 10, 5).translate(3, 3, 0) overlapping faces result contains union of faces'(assert) {
        const boxA = B2T.box(10, 10, 5, 'boxA').flipped(), boxB = boxA.translate(3, 3, 0)
        //const result = B2T.extrudeVertices([V3.O, V(0, 10), V(3, 10), V(3, 13), V(13, 13), V(13, 3), V(10, 3), V(10, 0)], P3.XY.flipped(), V(0, 0, 5), 'result').flipped()
        const result = new B2([
            new PlaneFace(new P3(V(0, -1, 0), -10), [
                new StraightEdge(new L3(V(3, 10, 0), V(0, 0, 1)), V(3, 10, 0), V(3, 10, 5), 0, 5),
                new StraightEdge(new L3(V(0, 10, 5), V(1, 0, 0)), V(3, 10, 5), V(0, 10, 5), 3, 0),
                new StraightEdge(new L3(V(0, 10, 0), V(0, 0, 1)), V(0, 10, 5), V(0, 10, 0), 5, 0),
                new StraightEdge(new L3(V(0, 10, 0), V(1, 0, 0)), V(0, 10, 0), V(3, 10, 0), 0, 3)]),
            new PlaneFace(new P3(V(-1, 0, 0), -10), [
                new StraightEdge(new L3(V(10, 3, 0), V(0, 0, -1)), V(10, 3, 5), V(10, 3, 0), -5, 0),
                new StraightEdge(new L3(V(10, 10, 0), V(0, -1, 0)), V(10, 3, 0), V(10, 0, 0), 7, 10),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 0, 1)), V(10, 0, 0), V(10, 0, 5), 0, 5),
                new StraightEdge(new L3(V(10, 10, 5), V(0, -1, 0)), V(10, 0, 5), V(10, 3, 5), 10, 7)]),
            new PlaneFace(new P3(V(0, 0, 1), 0), [
                new StraightEdge(new L3(V(3, 0, 0), V(0, 1, 0)), V(3, 3, 0), V(3, 10, 0), 3, 10),
                new StraightEdge(new L3(V(0, 10, 0), V(1, 0, 0)), V(3, 10, 0), V(0, 10, 0), 3, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 10, 0), V(0, 0, 0), 10, 0),
                new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(0, 0, 0), V(10, 0, 0), 10, 0),
                new StraightEdge(new L3(V(10, 10, 0), V(0, -1, 0)), V(10, 0, 0), V(10, 3, 0), 10, 7),
                new StraightEdge(new L3(V(0, 3, 0), V(-1, 0, 0)), V(10, 3, 0), V(3, 3, 0), -10, -3)]),
            new PlaneFace(new P3(V(0, 0, 1), 0), [
                new StraightEdge(new L3(V(3, 0, 0), V(0, 1, 0)), V(3, 10, 0), V(3, 3, 0), 10, 3),
                new StraightEdge(new L3(V(0, 3, 0), V(-1, 0, 0)), V(3, 3, 0), V(10, 3, 0), -3, -10),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 1, 0)), V(10, 3, 0), V(10, 10, 0), 3, 10),
                new StraightEdge(new L3(V(0, 10, 0), V(-1, 0, 0)), V(10, 10, 0), V(3, 10, 0), -10, -3)]),
            new PlaneFace(new P3(V(0, 0, -1), -5), [
                new StraightEdge(new L3(V(3, 0, 5), V(0, -1, 0)), V(3, 10, 5), V(3, 3, 5), -10, -3),
                new StraightEdge(new L3(V(0, 3, 5), V(1, 0, 0)), V(3, 3, 5), V(10, 3, 5), 3, 10),
                new StraightEdge(new L3(V(10, 10, 5), V(0, -1, 0)), V(10, 3, 5), V(10, 0, 5), 7, 10),
                new StraightEdge(new L3(V(10, 0, 5), V(-1, 0, 0)), V(10, 0, 5), V(0, 0, 5), 0, 10),
                new StraightEdge(new L3(V(0, 0, 5), V(0, 1, 0)), V(0, 0, 5), V(0, 10, 5), 0, 10),
                new StraightEdge(new L3(V(0, 10, 5), V(1, 0, 0)), V(0, 10, 5), V(3, 10, 5), 0, 3)]),
            new PlaneFace(new P3(V(0, 0, -1), -5), [
                new StraightEdge(new L3(V(3, 0, 5), V(0, -1, 0)), V(3, 3, 5), V(3, 10, 5), -3, -10),
                new StraightEdge(new L3(V(0, 10, 5), V(1, 0, 0)), V(3, 10, 5), V(10, 10, 5), 3, 10),
                new StraightEdge(new L3(V(10, 0, 5), V(0, -1, 0)), V(10, 10, 5), V(10, 3, 5), -10, -3),
                new StraightEdge(new L3(V(0, 3, 5), V(1, 0, 0)), V(10, 3, 5), V(3, 3, 5), 10, 3)]),
            new PlaneFace(new P3(V(1, 0, 0), 0), [
                new StraightEdge(new L3(V(0, 10, 0), V(0, 0, 1)), V(0, 10, 0), V(0, 10, 5), 0, 5),
                new StraightEdge(new L3(V(0, 0, 5), V(0, 1, 0)), V(0, 10, 5), V(0, 0, 5), 10, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 5), V(0, 0, 0), 5, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 0, 0), V(0, 10, 0), 0, 10)]),
            new PlaneFace(new P3(V(0, 1, 0), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 0), V(0, 0, 5), 0, 5),
                new StraightEdge(new L3(V(10, 0, 5), V(-1, 0, 0)), V(0, 0, 5), V(10, 0, 5), 10, 0),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 0, 1)), V(10, 0, 5), V(10, 0, 0), 5, 0),
                new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(10, 0, 0), V(0, 0, 0), 0, 10)]),
            new PlaneFace(new P3(V(1, 0, 0), 3), [
                new StraightEdge(new L3(V(3, 10, 0), V(0, 0, 1)), V(3, 10, 5), V(3, 10, 0), 5, 0),
                new StraightEdge(new L3(V(3, 3, 0), V(0, 1, 0)), V(3, 10, 0), V(3, 13, 0), 7, 10),
                new StraightEdge(new L3(V(3, 13, 0), V(0, 0, 1)), V(3, 13, 0), V(3, 13, 5), 0, 5),
                new StraightEdge(new L3(V(3, 3, 5), V(0, 1, 0)), V(3, 13, 5), V(3, 10, 5), 10, 7)]),
            new PlaneFace(new P3(V(0, 1, 0), 3), [
                new StraightEdge(new L3(V(10, 3, 0), V(0, 0, -1)), V(10, 3, 0), V(10, 3, 5), 0, -5),
                new StraightEdge(new L3(V(13, 3, 5), V(-1, 0, 0)), V(10, 3, 5), V(13, 3, 5), 3, 0),
                new StraightEdge(new L3(V(13, 3, 0), V(0, 0, 1)), V(13, 3, 5), V(13, 3, 0), 5, 0),
                new StraightEdge(new L3(V(13, 3, 0), V(-1, 0, 0)), V(13, 3, 0), V(10, 3, 0), 0, 3)]),
            new PlaneFace(new P3(V(0, 0, 1), 0), [
                new StraightEdge(new L3(V(0, 10, 0), V(-1, 0, 0)), V(3, 10, 0), V(10, 10, 0), -3, -10),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 1, 0)), V(10, 10, 0), V(10, 3, 0), 10, 3),
                new StraightEdge(new L3(V(13, 3, 0), V(-1, 0, 0)), V(10, 3, 0), V(13, 3, 0), 3, 0),
                new StraightEdge(new L3(V(13, 13, 0), V(0, -1, 0)), V(13, 3, 0), V(13, 13, 0), 10, 0),
                new StraightEdge(new L3(V(3, 13, 0), V(1, 0, 0)), V(13, 13, 0), V(3, 13, 0), 10, 0),
                new StraightEdge(new L3(V(3, 3, 0), V(0, 1, 0)), V(3, 13, 0), V(3, 10, 0), 10, 7)]),
            new PlaneFace(new P3(V(0, 0, -1), -5), [
                new StraightEdge(new L3(V(0, 10, 5), V(1, 0, 0)), V(10, 10, 5), V(3, 10, 5), 10, 3),
                new StraightEdge(new L3(V(3, 3, 5), V(0, 1, 0)), V(3, 10, 5), V(3, 13, 5), 7, 10),
                new StraightEdge(new L3(V(3, 13, 5), V(1, 0, 0)), V(3, 13, 5), V(13, 13, 5), 0, 10),
                new StraightEdge(new L3(V(13, 13, 5), V(0, -1, 0)), V(13, 13, 5), V(13, 3, 5), 0, 10),
                new StraightEdge(new L3(V(13, 3, 5), V(-1, 0, 0)), V(13, 3, 5), V(10, 3, 5), 0, 3),
                new StraightEdge(new L3(V(10, 0, 5), V(0, -1, 0)), V(10, 3, 5), V(10, 10, 5), -3, -10)]),
            new PlaneFace(new P3(V(0, -1, 0), -13), [
                new StraightEdge(new L3(V(13, 13, 0), V(0, 0, 1)), V(13, 13, 0), V(13, 13, 5), 0, 5),
                new StraightEdge(new L3(V(3, 13, 5), V(1, 0, 0)), V(13, 13, 5), V(3, 13, 5), 10, 0),
                new StraightEdge(new L3(V(3, 13, 0), V(0, 0, 1)), V(3, 13, 5), V(3, 13, 0), 5, 0),
                new StraightEdge(new L3(V(3, 13, 0), V(1, 0, 0)), V(3, 13, 0), V(13, 13, 0), 0, 10)]),
            new PlaneFace(new P3(V(-1, 0, 0), -13), [
                new StraightEdge(new L3(V(13, 3, 0), V(0, 0, 1)), V(13, 3, 0), V(13, 3, 5), 0, 5),
                new StraightEdge(new L3(V(13, 13, 5), V(0, -1, 0)), V(13, 3, 5), V(13, 13, 5), 10, 0),
                new StraightEdge(new L3(V(13, 13, 0), V(0, 0, 1)), V(13, 13, 5), V(13, 13, 0), 5, 0),
                new StraightEdge(new L3(V(13, 13, 0), V(0, -1, 0)), V(13, 13, 0), V(13, 3, 0), 0, 10)])], true)
        b2Equal(assert, boxA, boxB, boxA.intersection(boxB, true, true), result)
    },


    'B2T.box(10,10,5) + B2T.box(4,10,2).translate(2, 0, 3) overlapping faces result contains union of faces 2'(assert) {
        const box = B2T.box(10, 10, 5), box2 = B2T.box(4, 10, 2).translate(2, 0, 3)
        const result = new B2([
            new PlaneFace(new P3(V(0, 1, 0), 10), [
                new StraightEdge(new L3(V(0, 10, 3), V(-1, 0, 0)), V(2, 10, 3), V(6, 10, 3), -2, -6),
                new StraightEdge(new L3(V(6, 10, 0), V(0, 0, -1)), V(6, 10, 3), V(6, 10, 5), -3, -5),
                new StraightEdge(new L3(V(0, 10, 5), V(1, 0, 0)), V(6, 10, 5), V(10, 10, 5), 6, 10),
                new StraightEdge(new L3(V(10, 10, 0), V(0, 0, 1)), V(10, 10, 5), V(10, 10, 0), 5, 0),
                new StraightEdge(new L3(V(0, 10, 0), V(1, 0, 0)), V(10, 10, 0), V(0, 10, 0), 10, 0),
                new StraightEdge(new L3(V(0, 10, 0), V(0, 0, 1)), V(0, 10, 0), V(0, 10, 5), 0, 5),
                new StraightEdge(new L3(V(0, 10, 5), V(1, 0, 0)), V(0, 10, 5), V(2, 10, 5), 0, 2),
                new StraightEdge(new L3(V(2, 10, 0), V(0, 0, 1)), V(2, 10, 5), V(2, 10, 3), 5, 3)]),
            new PlaneFace(new P3(V(0, 1, 0), 10), [
                new StraightEdge(new L3(V(0, 10, 5), V(1, 0, 0)), V(2, 10, 5), V(6, 10, 5), 2, 6),
                new StraightEdge(new L3(V(6, 10, 0), V(0, 0, -1)), V(6, 10, 5), V(6, 10, 3), -5, -3),
                new StraightEdge(new L3(V(0, 10, 3), V(-1, 0, 0)), V(6, 10, 3), V(2, 10, 3), -6, -2),
                new StraightEdge(new L3(V(2, 10, 0), V(0, 0, 1)), V(2, 10, 3), V(2, 10, 5), 3, 5)]),
            new PlaneFace(new P3(V(0, -1, 0), 0), [
                new StraightEdge(new L3(V(10, 0, 5), V(-1, 0, 0)), V(2, 0, 5), V(0, 0, 5), 8, 10),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 5), V(0, 0, 0), 5, 0),
                new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(0, 0, 0), V(10, 0, 0), 10, 0),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 0, 1)), V(10, 0, 0), V(10, 0, 5), 0, 5),
                new StraightEdge(new L3(V(10, 0, 5), V(-1, 0, 0)), V(10, 0, 5), V(6, 0, 5), 0, 4),
                new StraightEdge(new L3(V(6, 0, 0), V(0, 0, 1)), V(6, 0, 5), V(6, 0, 3), 5, 3),
                new StraightEdge(new L3(V(0, 0, 3), V(1, 0, 0)), V(6, 0, 3), V(2, 0, 3), 6, 2),
                new StraightEdge(new L3(V(2, 0, 0), V(0, 0, -1)), V(2, 0, 3), V(2, 0, 5), -3, -5)]),
            new PlaneFace(new P3(V(0, -1, 0), 0), [
                new StraightEdge(new L3(V(0, 0, 3), V(1, 0, 0)), V(2, 0, 3), V(6, 0, 3), 2, 6),
                new StraightEdge(new L3(V(6, 0, 0), V(0, 0, 1)), V(6, 0, 3), V(6, 0, 5), 3, 5),
                new StraightEdge(new L3(V(0, 0, 5), V(-1, 0, 0)), V(6, 0, 5), V(2, 0, 5), -6, -2),
                new StraightEdge(new L3(V(2, 0, 0), V(0, 0, -1)), V(2, 0, 5), V(2, 0, 3), -5, -3)]),
            new PlaneFace(new P3(V(0, 0, 1), 5), [
                new StraightEdge(new L3(V(0, 10, 5), V(1, 0, 0)), V(2, 10, 5), V(0, 10, 5), 2, 0),
                new StraightEdge(new L3(V(0, 0, 5), V(0, 1, 0)), V(0, 10, 5), V(0, 0, 5), 10, 0),
                new StraightEdge(new L3(V(10, 0, 5), V(-1, 0, 0)), V(0, 0, 5), V(2, 0, 5), 10, 8),
                new StraightEdge(new L3(V(2, 0, 5), V(0, -1, 0)), V(2, 0, 5), V(2, 10, 5), 0, -10)]),
            new PlaneFace(new P3(V(0, 0, 1), 5), [
                new StraightEdge(new L3(V(0, 0, 5), V(-1, 0, 0)), V(2, 0, 5), V(6, 0, 5), -2, -6),
                new StraightEdge(new L3(V(6, 0, 5), V(0, 1, 0)), V(6, 0, 5), V(6, 10, 5), 0, 10),
                new StraightEdge(new L3(V(0, 10, 5), V(1, 0, 0)), V(6, 10, 5), V(2, 10, 5), 6, 2),
                new StraightEdge(new L3(V(2, 0, 5), V(0, -1, 0)), V(2, 10, 5), V(2, 0, 5), -10, 0)]),
            new PlaneFace(new P3(V(0, 0, 1), 5), [
                new StraightEdge(new L3(V(10, 0, 5), V(-1, 0, 0)), V(6, 0, 5), V(10, 0, 5), 4, 0),
                new StraightEdge(new L3(V(10, 10, 5), V(0, -1, 0)), V(10, 0, 5), V(10, 10, 5), 10, 0),
                new StraightEdge(new L3(V(0, 10, 5), V(1, 0, 0)), V(10, 10, 5), V(6, 10, 5), 10, 6),
                new StraightEdge(new L3(V(6, 0, 5), V(0, 1, 0)), V(6, 10, 5), V(6, 0, 5), 10, 0)]),
            new PlaneFace(new P3(V(-1, 0, 0), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 10, 0), V(0, 0, 0), 10, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 0), V(0, 0, 5), 0, 5),
                new StraightEdge(new L3(V(0, 0, 5), V(0, 1, 0)), V(0, 0, 5), V(0, 10, 5), 0, 10),
                new StraightEdge(new L3(V(0, 10, 0), V(0, 0, 1)), V(0, 10, 5), V(0, 10, 0), 5, 0)]),
            new PlaneFace(new P3(V(1, 0, 0), 10), [
                new StraightEdge(new L3(V(10, 10, 0), V(0, -1, 0)), V(10, 0, 0), V(10, 10, 0), 10, 0),
                new StraightEdge(new L3(V(10, 10, 0), V(0, 0, 1)), V(10, 10, 0), V(10, 10, 5), 0, 5),
                new StraightEdge(new L3(V(10, 10, 5), V(0, -1, 0)), V(10, 10, 5), V(10, 0, 5), 0, 10),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 0, 1)), V(10, 0, 5), V(10, 0, 0), 5, 0)]),
            new PlaneFace(new P3(V(0, 0, -1), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 0, 0), V(0, 10, 0), 0, 10),
                new StraightEdge(new L3(V(0, 10, 0), V(1, 0, 0)), V(0, 10, 0), V(10, 10, 0), 0, 10),
                new StraightEdge(new L3(V(10, 10, 0), V(0, -1, 0)), V(10, 10, 0), V(10, 0, 0), 0, 10),
                new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(10, 0, 0), V(0, 0, 0), 0, 10)])], false)
        b2Equal(assert, box, box2, box.plus(box2), result)
    },

    'B2T.box(10,10,10) + B2T.box(10,12,12).translate(0, 0, -2) overlapping faces result contains union of faces 3'(assert) {
        const box = B2T.box(10, 10, 10, 'box0'), box2 = B2T.box(10, 12, 12, 'box').translate(0, 0, -2)
        const result = new B2([
            new PlaneFace(new P3(V(-1, 0, 0), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 10, 0), V(0, 0, 0), 10, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 0), V(0, 0, 10), 0, 10),
                new StraightEdge(new L3(V(0, 0, 10), V(0, 1, 0)), V(0, 0, 10), V(0, 10, 10), 0, 10),
                new StraightEdge(new L3(V(0, 10, 0), V(0, 0, 1)), V(0, 10, 10), V(0, 10, 0), 10, 0)]),
            new PlaneFace(new P3(V(1, 0, 0), 10), [
                new StraightEdge(new L3(V(10, 0, 0), V(0, -1, 0)), V(10, 0, 0), V(10, 10, 0), 0, -10),
                new StraightEdge(new L3(V(10, 10, 0), V(0, 0, -1)), V(10, 10, 0), V(10, 10, 10), 0, -10),
                new StraightEdge(new L3(V(10, 0, 10), V(0, -1, 0)), V(10, 10, 10), V(10, 0, 10), -10, 0),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 0, -1)), V(10, 0, 10), V(10, 0, 0), -10, 0)]),
            new PlaneFace(new P3(V(0, -1, 0), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(-1, 0, 0)), V(0, 0, 0), V(10, 0, 0), 0, -10),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 0, -1)), V(10, 0, 0), V(10, 0, 10), 0, -10),
                new StraightEdge(new L3(V(0, 0, 10), V(-1, 0, 0)), V(10, 0, 10), V(0, 0, 10), -10, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 10), V(0, 0, 0), 10, 0)]),
            new PlaneFace(new P3(V(0, 0, 1), 10), [
                new StraightEdge(new L3(V(0, 0, 10), V(-1, 0, 0)), V(0, 0, 10), V(10, 0, 10), 0, -10),
                new StraightEdge(new L3(V(10, 0, 10), V(0, -1, 0)), V(10, 0, 10), V(10, 10, 10), 0, -10),
                new StraightEdge(new L3(V(0, 10, 10), V(1, 0, 0)), V(10, 10, 10), V(0, 10, 10), 10, 0),
                new StraightEdge(new L3(V(0, 0, 10), V(0, 1, 0)), V(0, 10, 10), V(0, 0, 10), 10, 0)]),
            new PlaneFace(new P3(V(-1, 0, 0), 0), [
                new StraightEdge(new L3(V(0, 0, 10), V(0, 1, 0)), V(0, 10, 10), V(0, 12, 10), 10, 12),
                new StraightEdge(new L3(V(0, 12, -2), V(0, 0, 1)), V(0, 12, 10), V(0, 12, -2), 12, 0),
                new StraightEdge(new L3(V(0, 0, -2), V(0, 1, 0)), V(0, 12, -2), V(0, 0, -2), 12, 0),
                new StraightEdge(new L3(V(0, 0, -2), V(0, 0, 1)), V(0, 0, -2), V(0, 0, 0), 0, 2),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 0, 0), V(0, 10, 0), 0, 10),
                new StraightEdge(new L3(V(0, 10, 0), V(0, 0, 1)), V(0, 10, 0), V(0, 10, 10), 0, 10)]),
            new PlaneFace(new P3(V(1, 0, 0), 10), [
                new StraightEdge(new L3(V(10, 0, 0), V(0, -1, 0)), V(10, 10, 0), V(10, 0, 0), -10, 0),
                new StraightEdge(new L3(V(10, 0, -2), V(0, 0, 1)), V(10, 0, 0), V(10, 0, -2), 2, 0),
                new StraightEdge(new L3(V(10, 12, -2), V(0, -1, 0)), V(10, 0, -2), V(10, 12, -2), 12, 0),
                new StraightEdge(new L3(V(10, 12, -2), V(0, 0, 1)), V(10, 12, -2), V(10, 12, 10), 0, 12),
                new StraightEdge(new L3(V(10, 12, 10), V(0, -1, 0)), V(10, 12, 10), V(10, 10, 10), 0, 2),
                new StraightEdge(new L3(V(10, 10, 0), V(0, 0, -1)), V(10, 10, 10), V(10, 10, 0), -10, 0)]),
            new PlaneFace(new P3(V(0, -1, 0), 0), [
                new StraightEdge(new L3(V(0, 0, -2), V(0, 0, 1)), V(0, 0, 0), V(0, 0, -2), 2, 0),
                new StraightEdge(new L3(V(10, 0, -2), V(-1, 0, 0)), V(0, 0, -2), V(10, 0, -2), 10, 0),
                new StraightEdge(new L3(V(10, 0, -2), V(0, 0, 1)), V(10, 0, -2), V(10, 0, 0), 0, 2),
                new StraightEdge(new L3(V(0, 0, 0), V(-1, 0, 0)), V(10, 0, 0), V(0, 0, 0), -10, 0)]),
            new PlaneFace(new P3(V(0, 0, 1), 10), [
                new StraightEdge(new L3(V(10, 12, 10), V(0, -1, 0)), V(10, 10, 10), V(10, 12, 10), 2, 0),
                new StraightEdge(new L3(V(0, 12, 10), V(1, 0, 0)), V(10, 12, 10), V(0, 12, 10), 10, 0),
                new StraightEdge(new L3(V(0, 0, 10), V(0, 1, 0)), V(0, 12, 10), V(0, 10, 10), 12, 10),
                new StraightEdge(new L3(V(0, 10, 10), V(1, 0, 0)), V(0, 10, 10), V(10, 10, 10), 0, 10)]),
            new PlaneFace(new P3(V(0, 1, 0), 12), [
                new StraightEdge(new L3(V(0, 12, -2), V(1, 0, 0)), V(10, 12, -2), V(0, 12, -2), 10, 0),
                new StraightEdge(new L3(V(0, 12, -2), V(0, 0, 1)), V(0, 12, -2), V(0, 12, 10), 0, 12),
                new StraightEdge(new L3(V(0, 12, 10), V(1, 0, 0)), V(0, 12, 10), V(10, 12, 10), 0, 10),
                new StraightEdge(new L3(V(10, 12, -2), V(0, 0, 1)), V(10, 12, 10), V(10, 12, -2), 12, 0)]),
            new PlaneFace(new P3(V(0, 0, -1), 2), [
                new StraightEdge(new L3(V(0, 0, -2), V(0, 1, 0)), V(0, 0, -2), V(0, 12, -2), 0, 12),
                new StraightEdge(new L3(V(0, 12, -2), V(1, 0, 0)), V(0, 12, -2), V(10, 12, -2), 0, 10),
                new StraightEdge(new L3(V(10, 12, -2), V(0, -1, 0)), V(10, 12, -2), V(10, 0, -2), 0, 12),
                new StraightEdge(new L3(V(10, 0, -2), V(-1, 0, 0)), V(10, 0, -2), V(0, 0, -2), 0, 10)])], false)
        b2Equal(assert, box, box2, box.plus(box2), result)
    },

    'B2T.box(10, 10, 10) - B2T.box().scale(8 / 3 ** 0.5).rotateAB(V3.ONES, V3.X)'(assert) {
        const box = B2T.box(10, 10, 10, 'box0'), box2 = B2T.box(1,1,1).scale(8 / 3 ** 0.5).rotateAB(V3.ONES, V3.X)
        const result = new B2([
            new PlaneFace(new P3(V(0, -1, 0), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(0.7071067811865476, 0, 0.7071067811865475)), V(0, 0, 0), V(3.381197846482994, 0, 3.3811978464829933), 0, 4.781735851562952),
                new StraightEdge(new L3(V(2.791322082992153, 0, 3.813016875511774), V(0.8068982213550735, 0, -0.5906904945688721)), V(3.381197846482994, 0, 3.3811978464829933), V(8, 0, 0), 0.7310411002024879, 6.45518577084059),
                new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(8, 0, 0), V(10, 0, 0), 2, 0),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 0, 1)), V(10, 0, 0), V(10, 0, 10), 0, 10),
                new StraightEdge(new L3(V(10, 0, 10), V(-1, 0, 0)), V(10, 0, 10), V(0, 0, 10), 0, 10),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 10), V(0, 0, 0), 10, 0)]),
            new PlaneFace(new P3(V(0, 0, -1), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(-0.7071067811865476, -0.7071067811865475, 0)), V(3.381197846482994, 3.3811978464829933, 0), V(0, 0, 0), -4.781735851562952, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 0, 0), V(0, 10, 0), 0, 10),
                new StraightEdge(new L3(V(0, 10, 0), V(1, 0, 0)), V(0, 10, 0), V(10, 10, 0), 0, 10),
                new StraightEdge(new L3(V(10, 10, 0), V(0, -1, 0)), V(10, 10, 0), V(10, 0, 0), 0, 10),
                new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(10, 0, 0), V(8, 0, 0), 0, 2),
                new StraightEdge(new L3(V(2.7913220829921492, 3.813016875511774, 0), V(-0.8068982213550737, 0.590690494568872, 0)), V(8, 0, 0), V(3.381197846482994, 3.3811978464829933, 0), -6.455185770840593, -0.7310411002024895)]),
            new PlaneFace(new P3(V(-1, 0, 0), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 10, 0), V(0, 0, 0), 10, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 0), V(0, 0, 10), 0, 10),
                new StraightEdge(new L3(V(0, 0, 10), V(0, 1, 0)), V(0, 0, 10), V(0, 10, 10), 0, 10),
                new StraightEdge(new L3(V(0, 10, 0), V(0, 0, 1)), V(0, 10, 10), V(0, 10, 0), 10, 0)]),
            new PlaneFace(new P3(V(0, 1, 0), 10), [
                new StraightEdge(new L3(V(0, 10, 0), V(1, 0, 0)), V(10, 10, 0), V(0, 10, 0), 10, 0),
                new StraightEdge(new L3(V(0, 10, 0), V(0, 0, 1)), V(0, 10, 0), V(0, 10, 10), 0, 10),
                new StraightEdge(new L3(V(0, 10, 10), V(1, 0, 0)), V(0, 10, 10), V(10, 10, 10), 0, 10),
                new StraightEdge(new L3(V(10, 10, 0), V(0, 0, 1)), V(10, 10, 10), V(10, 10, 0), 10, 0)]),
            new PlaneFace(new P3(V(1, 0, 0), 10), [
                new StraightEdge(new L3(V(10, 10, 0), V(0, -1, 0)), V(10, 0, 0), V(10, 10, 0), 10, 0),
                new StraightEdge(new L3(V(10, 10, 0), V(0, 0, 1)), V(10, 10, 0), V(10, 10, 10), 0, 10),
                new StraightEdge(new L3(V(10, 10, 10), V(0, -1, 0)), V(10, 10, 10), V(10, 0, 10), 0, 10),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 0, 1)), V(10, 0, 10), V(10, 0, 0), 10, 0)]),
            new PlaneFace(new P3(V(0, 0, 1), 10), [
                new StraightEdge(new L3(V(10, 0, 10), V(-1, 0, 0)), V(0, 0, 10), V(10, 0, 10), 10, 0),
                new StraightEdge(new L3(V(10, 10, 10), V(0, -1, 0)), V(10, 0, 10), V(10, 10, 10), 10, 0),
                new StraightEdge(new L3(V(0, 10, 10), V(1, 0, 0)), V(10, 10, 10), V(0, 10, 10), 10, 0),
                new StraightEdge(new L3(V(0, 0, 10), V(0, 1, 0)), V(0, 10, 10), V(0, 0, 10), 10, 0)]),
            new PlaneFace(new P3(V(0.5773502691896265, -0.5773502691896265, -0.5773502691896243), -1.509903313490213e-14), [
                new StraightEdge(new L3(V(0, 0, 0), V(0.7071067811865476, 0, 0.7071067811865475)), V(3.381197846482994, 0, 3.3811978464829933), V(0, 0, 0), 4.781735851562952, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(-0.7071067811865476, -0.7071067811865475, 0)), V(0, 0, 0), V(3.381197846482994, 3.3811978464829933, 0), 0, -4.781735851562952),
                new StraightEdge(new L3(V(2.666666666666668, 3.642734410091836, -0.9760677434251697), V(0.5773502691896258, -0.21132486540518713, 0.7886751345948129)), V(3.381197846482994, 3.3811978464829933, 0), V(5.333333333333336, 2.666666666666666, 2.666666666666666), 1.2376043070340121, 4.618802153517006),
                new StraightEdge(new L3(V(2.666666666666668, -0.9760677434251697, 3.642734410091836), V(0.5773502691896258, 0.7886751345948129, -0.21132486540518713)), V(5.333333333333336, 2.666666666666666, 2.666666666666666), V(3.381197846482994, 0, 3.3811978464829933), 4.618802153517006, 1.2376043070340121)]),
            new PlaneFace(new P3(V(-0.5773502691896246, -0.7886751345948133, 0.2113248654051885), -4.618802153517008), [
                new StraightEdge(new L3(V(2.7913220829921492, 3.813016875511774, 0), V(-0.8068982213550737, 0.590690494568872, 0)), V(3.381197846482994, 3.3811978464829933, 0), V(8, 0, 0), -0.7310411002024895, -6.455185770840593),
                new StraightEdge(new L3(V(5.333333333333336, 2.666666666666666, 2.666666666666666), V(0.5773502691896257, -0.5773502691896258, -0.5773502691896258)), V(8, -6.661338147750939e-16, -8.881784197001252e-16), V(5.333333333333336, 2.666666666666666, 2.666666666666666), 4.618802153517006, 0),
                new StraightEdge(new L3(V(2.666666666666668, 3.642734410091836, -0.9760677434251697), V(0.5773502691896258, -0.21132486540518713, 0.7886751345948129)), V(5.333333333333336, 2.666666666666666, 2.666666666666666), V(3.381197846482994, 3.3811978464829933, 0), 4.618802153517006, 1.2376043070340121)]),
            new PlaneFace(new P3(V(-0.5773502691896261, 0.2113248654051888, -0.7886751345948121), -4.618802153517), [
                new StraightEdge(new L3(V(2.791322082992153, 0, 3.813016875511774), V(0.8068982213550735, 0, -0.5906904945688721)), V(8, 0, 0), V(3.381197846482994, 0, 3.3811978464829933), 6.45518577084059, 0.7310411002024879),
                new StraightEdge(new L3(V(2.666666666666668, -0.9760677434251697, 3.642734410091836), V(0.5773502691896258, 0.7886751345948129, -0.21132486540518713)), V(3.381197846482994, 0, 3.3811978464829933), V(5.333333333333336, 2.666666666666666, 2.666666666666666), 1.2376043070340121, 4.618802153517006),
                new StraightEdge(new L3(V(5.333333333333336, 2.666666666666666, 2.666666666666666), V(0.5773502691896257, -0.5773502691896258, -0.5773502691896258)), V(5.333333333333336, 2.666666666666666, 2.666666666666666), V(8, -6.661338147750939e-16, -8.881784197001252e-16), 0, 4.618802153517006)])], false)
        b2Equal(assert, box, box2, box.minus(box2), result)
    },

    'B2T.box(10, 10, 10) - B2T.tetrahedron(V(5,0,10),V(5,3,10),V(7,1,9),V(7,1,11))'(assert) {
        const box = B2T.box(10, 10, 10, 'box0'), box2 = B2T.tetrahedron(V(5,0,10),V(5,3,10),V(7,1,9),V(7,1,11)).flipped()
        const result = new B2([
            new PlaneFace(new P3(V(0, -1, 0), 0), [
                new StraightEdge(new L3(V(10, 0, 10), V(-1, 0, 0)), V(10, 0, 10), V(5, 0, 10), 0, 5),
                new StraightEdge(new L3(V(10, 0, 10), V(-1, 0, 0)), V(5, 0, 10), V(0, 0, 10), 5, 10),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 10), V(0, 0, 0), 10, 0),
                new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(0, 0, 0), V(10, 0, 0), 10, 0),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 0, 1)), V(10, 0, 0), V(10, 0, 10), 0, 10)]),
            new PlaneFace(new P3(V(0, 0, 1), 10), [
                new StraightEdge(new L3(V(5, 0, 10), V(0, 1, 0)), V(5, 0, 10), V(5, 3, 10), 0, 3),
                new StraightEdge(new L3(V(4, 4, 10), V(0.7071067811865475, -0.7071067811865475, 0)), V(5, 3, 10), V(7, 1, 10), 1.414213562373095, 4.242640687119285),
                new StraightEdge(new L3(V(1, -2, 10), V(-0.8944271909999159, -0.4472135954999579, 0)), V(7, 1, 10), V(5, 0, 10), -6.708203932499369, -4.47213595499958),
                new StraightEdge(new L3(V(10, 0, 10), V(-1, 0, 0)), V(5, 0, 10), V(10, 0, 10), 5, 0),
                new StraightEdge(new L3(V(10, 10, 10), V(0, -1, 0)), V(10, 0, 10), V(10, 10, 10), 10, 0),
                new StraightEdge(new L3(V(0, 10, 10), V(1, 0, 0)), V(10, 10, 10), V(0, 10, 10), 10, 0),
                new StraightEdge(new L3(V(0, 0, 10), V(0, 1, 0)), V(0, 10, 10), V(0, 0, 10), 10, 0),
                new StraightEdge(new L3(V(10, 0, 10), V(-1, 0, 0)), V(0, 0, 10), V(5, 0, 10), 10, 5)]),
            new PlaneFace(new P3(V(-1, 0, 0), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 10, 0), V(0, 0, 0), 10, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 0), V(0, 0, 10), 0, 10),
                new StraightEdge(new L3(V(0, 0, 10), V(0, 1, 0)), V(0, 0, 10), V(0, 10, 10), 0, 10),
                new StraightEdge(new L3(V(0, 10, 0), V(0, 0, 1)), V(0, 10, 10), V(0, 10, 0), 10, 0)]),
            new PlaneFace(new P3(V(0, 1, 0), 10), [
                new StraightEdge(new L3(V(0, 10, 0), V(1, 0, 0)), V(10, 10, 0), V(0, 10, 0), 10, 0),
                new StraightEdge(new L3(V(0, 10, 0), V(0, 0, 1)), V(0, 10, 0), V(0, 10, 10), 0, 10),
                new StraightEdge(new L3(V(0, 10, 10), V(1, 0, 0)), V(0, 10, 10), V(10, 10, 10), 0, 10),
                new StraightEdge(new L3(V(10, 10, 0), V(0, 0, 1)), V(10, 10, 10), V(10, 10, 0), 10, 0)]),
            new PlaneFace(new P3(V(1, 0, 0), 10), [
                new StraightEdge(new L3(V(10, 10, 0), V(0, -1, 0)), V(10, 0, 0), V(10, 10, 0), 10, 0),
                new StraightEdge(new L3(V(10, 10, 0), V(0, 0, 1)), V(10, 10, 0), V(10, 10, 10), 0, 10),
                new StraightEdge(new L3(V(10, 10, 10), V(0, -1, 0)), V(10, 10, 10), V(10, 0, 10), 0, 10),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 0, 1)), V(10, 0, 10), V(10, 0, 0), 10, 0)]),
            new PlaneFace(new P3(V(0, 0, -1), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 0, 0), V(0, 10, 0), 0, 10),
                new StraightEdge(new L3(V(0, 10, 0), V(1, 0, 0)), V(0, 10, 0), V(10, 10, 0), 0, 10),
                new StraightEdge(new L3(V(10, 10, 0), V(0, -1, 0)), V(10, 10, 0), V(10, 0, 0), 0, 10),
                new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(10, 0, 0), V(0, 0, 0), 0, 10)]),
            new PlaneFace(new P3(V(0.4472135954999579, 0, 0.8944271909999159), 11.180339887498949), [
                new StraightEdge(new L3(V(5, 0, 10), V(0, 1, 0)), V(5, 3, 10), V(5, 0, 10), 3, 0),
                new StraightEdge(new L3(V(5, 0, 10), V(0.8164965809277261, 0.4082482904638631, -0.4082482904638631)), V(5, 0, 10), V(7, 1, 9), 0, 2.449489742783178),
                new StraightEdge(new L3(V(5, 3, 10), V(0.6666666666666666, -0.6666666666666666, -0.3333333333333333)), V(7, 1, 9), V(5, 3, 10), 3, 0)]),
            new PlaneFace(new P3(V(-0.7071067811865475, -0.7071067811865475, 0), -5.65685424949238), [
                new StraightEdge(new L3(V(4, 4, 10), V(0.7071067811865475, -0.7071067811865475, 0)), V(7, 1, 10), V(5, 3, 10), 4.242640687119285, 1.414213562373095),
                new StraightEdge(new L3(V(5, 3, 10), V(0.6666666666666666, -0.6666666666666666, -0.3333333333333333)), V(5, 3, 10), V(7, 1, 9), 0, 3),
                new StraightEdge(new L3(V(7, 1, 9), V(0, 0, 1)), V(7, 1, 9), V(7, 1, 10), 0, 1)]),
            new PlaneFace(new P3(V(-0.4472135954999579, 0.8944271909999159, 0), -2.23606797749979), [
                new StraightEdge(new L3(V(1, -2, 10), V(-0.8944271909999159, -0.4472135954999579, 0)), V(5, 0, 10), V(7, 1, 10), -4.47213595499958, -6.708203932499369),
                new StraightEdge(new L3(V(7, 1, 9), V(0, 0, 1)), V(7, 1, 10), V(7, 1, 9), 1, 0),
                new StraightEdge(new L3(V(5, 0, 10), V(0.8164965809277261, 0.4082482904638631, -0.4082482904638631)), V(7, 1, 9), V(5, 0, 10), 2.449489742783178, 0)])], false)
        b2Equal(assert, box, box2, box.and(box2), result)
    },

    'B2T.box(10, 10, 10) - B2T.box().rotateAB(V3.ONES, V3.Y).translate(5,0,10)'(assert) {
        const box = B2T.box(10, 10, 10, 'box0'), box2 = B2T.box().rotateAB(V3.ONES, V3.Y).translate(5,0,10).flipped()
        const result = new B2([
            new PlaneFace(new P3(V(0, -1, 0), 0), [
                new StraightEdge(new L3(V(10, 0, 10), V(-1, 0, 0)), V(10, 0, 10), V(5, 0, 10), 0, 5),
                new StraightEdge(new L3(V(10, 0, 10), V(-1, 0, 0)), V(5, 0, 10), V(0, 0, 10), 5, 10),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 10), V(0, 0, 0), 10, 0),
                new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(0, 0, 0), V(10, 0, 0), 10, 0),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 0, 1)), V(10, 0, 0), V(10, 0, 10), 0, 10)]),
            new PlaneFace(new P3(V(0, 0, 1), 10), [
                new StraightEdge(new L3(V(3.2554236981299045, 2.3831355471948585, 10), V(-0.5906904945688723, 0.8068982213550733, 0)), V(5, 0, 10), V(4.267949192431122, 1.0000000000000002, 10), -2.953452472844362, -1.7141387979168856),
                new StraightEdge(new L3(V(1.6339745962155634, -1.6339745962155598, 10), V(0.7071067811865471, 0.707106781186548, 0)), V(4.267949192431122, 1.0000000000000002, 10), V(5, 1.7320508075688772, 10), 3.725002596914241, 4.7602787773243245),
                new StraightEdge(new L3(V(4.080966067942895, 2.987474505698783, 10), V(0.5906904945688725, -0.8068982213550733, 0)), V(5, 1.7320508075688772, 10), V(5.732050807568878, 0.7320508075688772, 10), 1.5558637569204163, 2.7951774318478924),
                new StraightEdge(new L3(V(2.499999999999999, -2.5, 10), V(-0.7071067811865477, -0.7071067811865474, 0)), V(5.732050807568878, 0.7320508075688772, 10), V(5, 0, 10), -4.570810086342822, -3.535533905932738),
                new StraightEdge(new L3(V(10, 0, 10), V(-1, 0, 0)), V(5, 0, 10), V(10, 0, 10), 5, 0),
                new StraightEdge(new L3(V(10, 10, 10), V(0, -1, 0)), V(10, 0, 10), V(10, 10, 10), 10, 0),
                new StraightEdge(new L3(V(0, 10, 10), V(1, 0, 0)), V(10, 10, 10), V(0, 10, 10), 10, 0),
                new StraightEdge(new L3(V(0, 0, 10), V(0, 1, 0)), V(0, 10, 10), V(0, 0, 10), 10, 0),
                new StraightEdge(new L3(V(10, 0, 10), V(-1, 0, 0)), V(0, 0, 10), V(5, 0, 10), 10, 5)]),
            new PlaneFace(new P3(V(-1, 0, 0), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 10, 0), V(0, 0, 0), 10, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 0), V(0, 0, 10), 0, 10),
                new StraightEdge(new L3(V(0, 0, 10), V(0, 1, 0)), V(0, 0, 10), V(0, 10, 10), 0, 10),
                new StraightEdge(new L3(V(0, 10, 0), V(0, 0, 1)), V(0, 10, 10), V(0, 10, 0), 10, 0)]),
            new PlaneFace(new P3(V(0, 1, 0), 10), [
                new StraightEdge(new L3(V(0, 10, 0), V(1, 0, 0)), V(10, 10, 0), V(0, 10, 0), 10, 0),
                new StraightEdge(new L3(V(0, 10, 0), V(0, 0, 1)), V(0, 10, 0), V(0, 10, 10), 0, 10),
                new StraightEdge(new L3(V(0, 10, 10), V(1, 0, 0)), V(0, 10, 10), V(10, 10, 10), 0, 10),
                new StraightEdge(new L3(V(10, 10, 0), V(0, 0, 1)), V(10, 10, 10), V(10, 10, 0), 10, 0)]),
            new PlaneFace(new P3(V(1, 0, 0), 10), [
                new StraightEdge(new L3(V(10, 10, 0), V(0, -1, 0)), V(10, 0, 0), V(10, 10, 0), 10, 0),
                new StraightEdge(new L3(V(10, 10, 0), V(0, 0, 1)), V(10, 10, 0), V(10, 10, 10), 0, 10),
                new StraightEdge(new L3(V(10, 10, 10), V(0, -1, 0)), V(10, 10, 10), V(10, 0, 10), 0, 10),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 0, 1)), V(10, 0, 10), V(10, 0, 0), 10, 0)]),
            new PlaneFace(new P3(V(0, 0, -1), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 0, 0), V(0, 10, 0), 0, 10),
                new StraightEdge(new L3(V(0, 10, 0), V(1, 0, 0)), V(0, 10, 0), V(10, 10, 0), 0, 10),
                new StraightEdge(new L3(V(10, 10, 0), V(0, -1, 0)), V(10, 10, 0), V(10, 0, 0), 0, 10),
                new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(10, 0, 0), V(0, 0, 0), 0, 10)]),
            new PlaneFace(new P3(V(0.7886751345948128, 0.577350269189626, -0.21132486540518747), 1.8301270189221892), [
                new StraightEdge(new L3(V(3.2554236981299045, 2.3831355471948585, 10), V(-0.5906904945688723, 0.8068982213550733, 0)), V(4.267949192431122, 1.0000000000000002, 10), V(5, 0, 10), -1.7141387979168856, -2.953452472844362),
                new StraightEdge(new L3(V(5, 0, 10), V(-0.5773502691896258, 0.5773502691896257, -0.5773502691896258)), V(5, 0, 10), V(4.422649730810374, 0.5773502691896256, 9.422649730810374), 0, 0.9999999999999999),
                new StraightEdge(new L3(V(4.422649730810374, 0.5773502691896256, 9.422649730810374), V(-0.21132486540518713, 0.5773502691896258, 0.7886751345948129)), V(4.422649730810374, 0.5773502691896256, 9.422649730810374), V(4.267949192431122, 1.0000000000000002, 10), 0, 0.7320508075688779)]),
            new PlaneFace(new P3(V(0.5773502691896257, -0.5773502691896251, 0.5773502691896264), 7.660254037844392), [
                new StraightEdge(new L3(V(1.6339745962155634, -1.6339745962155598, 10), V(0.7071067811865471, 0.707106781186548, 0)), V(5, 1.7320508075688772, 10), V(4.267949192431122, 1.0000000000000002, 10), 4.7602787773243245, 3.725002596914241),
                new StraightEdge(new L3(V(4.422649730810374, 0.5773502691896256, 9.422649730810374), V(-0.21132486540518713, 0.5773502691896258, 0.7886751345948129)), V(4.267949192431122, 1.0000000000000002, 10), V(4.422649730810374, 0.5773502691896256, 9.422649730810374), 0.7320508075688779, 0),
                new StraightEdge(new L3(V(4.422649730810374, 0.5773502691896256, 9.422649730810374), V(0.7886751345948129, 0.5773502691896258, -0.21132486540518713)), V(4.422649730810374, 0.5773502691896256, 9.422649730810374), V(5.211324865405187, 1.1547005383792515, 9.211324865405187), 0, 0.9999999999999998),
                new StraightEdge(new L3(V(5.211324865405187, 1.1547005383792515, 9.211324865405187), V(-0.21132486540518713, 0.5773502691896258, 0.7886751345948129)), V(5.211324865405187, 1.1547005383792515, 9.211324865405187), V(5, 1.7320508075688772, 10), 0, 0.9999999999999998)]),
            new PlaneFace(new P3(V(-0.7886751345948126, -0.577350269189626, 0.21132486540518705), -2.8301270189221928), [
                new StraightEdge(new L3(V(4.080966067942895, 2.987474505698783, 10), V(0.5906904945688725, -0.8068982213550733, 0)), V(5.732050807568878, 0.7320508075688772, 10), V(5, 1.7320508075688772, 10), 2.7951774318478924, 1.5558637569204163),
                new StraightEdge(new L3(V(5.211324865405187, 1.1547005383792515, 9.211324865405187), V(-0.21132486540518713, 0.5773502691896258, 0.7886751345948129)), V(5, 1.7320508075688772, 10), V(5.211324865405187, 1.1547005383792515, 9.211324865405187), 0.9999999999999998, 0),
                new StraightEdge(new L3(V(5.211324865405187, 1.1547005383792515, 9.211324865405187), V(0.5773502691896258, -0.5773502691896257, 0.5773502691896258)), V(5.211324865405187, 1.1547005383792515, 9.211324865405187), V(5.788675134594813, 0.5773502691896257, 9.788675134594813), 0, 0.9999999999999999),
                new StraightEdge(new L3(V(5.788675134594813, 0.5773502691896257, 9.788675134594813), V(-0.21132486540518713, 0.5773502691896258, 0.7886751345948129)), V(5.788675134594813, 0.5773502691896257, 9.788675134594813), V(5.732050807568878, 0.7320508075688772, 10), 0, 0.2679491924311225)]),
            new PlaneFace(new P3(V(-0.5773502691896255, 0.5773502691896258, -0.577350269189626), -8.660254037844387), [
                new StraightEdge(new L3(V(2.499999999999999, -2.5, 10), V(-0.7071067811865477, -0.7071067811865474, 0)), V(5, 0, 10), V(5.732050807568878, 0.7320508075688772, 10), -3.535533905932738, -4.570810086342822),
                new StraightEdge(new L3(V(5.788675134594813, 0.5773502691896257, 9.788675134594813), V(-0.21132486540518713, 0.5773502691896258, 0.7886751345948129)), V(5.732050807568878, 0.7320508075688772, 10), V(5.788675134594813, 0.5773502691896257, 9.788675134594813), 0.2679491924311225, 0),
                new StraightEdge(new L3(V(5.788675134594813, 0.5773502691896257, 9.788675134594813), V(-0.7886751345948129, -0.5773502691896258, 0.21132486540518713)), V(5.788675134594813, 0.5773502691896257, 9.788675134594813), V(5, 0, 10), 0, 0.9999999999999998)]),
            new PlaneFace(new P3(V(-0.2113248654051868, 0.5773502691896251, 0.7886751345948136), 6.830127018922202), [
                new StraightEdge(new L3(V(5.788675134594813, 0.5773502691896257, 9.788675134594813), V(-0.7886751345948129, -0.5773502691896258, 0.21132486540518713)), V(5, 0, 10), V(5.788675134594813, 0.5773502691896257, 9.788675134594813), 0.9999999999999998, 0),
                new StraightEdge(new L3(V(5.211324865405187, 1.1547005383792515, 9.211324865405187), V(0.5773502691896258, -0.5773502691896257, 0.5773502691896258)), V(5.788675134594813, 0.5773502691896257, 9.788675134594813), V(5.211324865405187, 1.1547005383792515, 9.211324865405187), 0.9999999999999999, 0),
                new StraightEdge(new L3(V(4.422649730810374, 0.5773502691896256, 9.422649730810374), V(0.7886751345948129, 0.5773502691896258, -0.21132486540518713)), V(5.211324865405187, 1.1547005383792515, 9.211324865405187), V(4.422649730810374, 0.5773502691896256, 9.422649730810374), 0.9999999999999998, 0),
                new StraightEdge(new L3(V(5, 0, 10), V(-0.5773502691896258, 0.5773502691896257, -0.5773502691896258)), V(4.422649730810374, 0.5773502691896256, 9.422649730810374), V(5, 0, 10), 0.9999999999999999, 0)])], false)
        b2Equal(assert, box, box2, box.and(box2), result)
    },

    'menger(1) - B2T.box(2,2,1).rotateZ(-45*DEG)'(assert) {
        const box = B2T.box(1,1,1, 'box0').minus(B2T.box(1/3,1/3,1).translate(1/3,1/3,0))
        const box2 = B2T.box(2,2,1).rotateZ(-45*DEG).flipped()
        const result = new B2([
            new PlaneFace(new P3(V(0, 0, -1), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(-0.7071067811865475, -0.7071067811865476, 0)), V(1, 1, 0), V(0.6666666666666666, 0.6666666666666666, 0), -1.414213562373095, -0.9428090415820634),
                new StraightEdge(new L3(V(0, 0.6666666666666666, 0), V(-1, 0, 0)), V(0.6666666666666666, 0.6666666666666666, 0), V(0.3333333333333333, 0.6666666666666666, 0), -0.6666666666666666, -0.3333333333333333),
                new StraightEdge(new L3(V(0.3333333333333333, 0, 0), V(0, -1, 0)), V(0.3333333333333333, 0.6666666666666666, 0), V(0.3333333333333333, 0.3333333333333333, 0), -0.6666666666666666, -0.3333333333333333),
                new StraightEdge(new L3(V(0, 0, 0), V(-0.7071067811865475, -0.7071067811865476, 0)), V(0.3333333333333333, 0.3333333333333333, 0), V(0, 0, 0), -0.4714045207910317, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 0, 0), V(0, 1, 0), 0, 1),
                new StraightEdge(new L3(V(0, 1, 0), V(1, 0, 0)), V(0, 1, 0), V(1, 1, 0), 0, 1)]),
            new PlaneFace(new P3(V(0, 0, 1), 1), [
                new StraightEdge(new L3(V(0, 0, 1), V(0.7071067811865475, 0.7071067811865476, 0)), V(0, 0, 1), V(0.3333333333333333, 0.3333333333333333, 1), 0, 0.4714045207910317),
                new StraightEdge(new L3(V(0.3333333333333333, 0, 1), V(0, 1, 0)), V(0.3333333333333333, 0.3333333333333333, 1), V(0.3333333333333333, 0.6666666666666666, 1), 0.3333333333333333, 0.6666666666666666),
                new StraightEdge(new L3(V(0, 0.6666666666666666, 1), V(1, 0, 0)), V(0.3333333333333333, 0.6666666666666666, 1), V(0.6666666666666666, 0.6666666666666666, 1), 0.3333333333333333, 0.6666666666666666),
                new StraightEdge(new L3(V(0, 0, 1), V(0.7071067811865475, 0.7071067811865476, 0)), V(0.6666666666666666, 0.6666666666666666, 1), V(1, 1, 1), 0.9428090415820634, 1.414213562373095),
                new StraightEdge(new L3(V(0, 1, 1), V(1, 0, 0)), V(1, 1, 1), V(0, 1, 1), 1, 0),
                new StraightEdge(new L3(V(0, 0, 1), V(0, 1, 0)), V(0, 1, 1), V(0, 0, 1), 1, 0)]),
            new PlaneFace(new P3(V(-1, 0, 0), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 0), V(0, 0, 1), 0, 1),
                new StraightEdge(new L3(V(0, 0, 1), V(0, 1, 0)), V(0, 0, 1), V(0, 1, 1), 0, 1),
                new StraightEdge(new L3(V(0, 1, 0), V(0, 0, 1)), V(0, 1, 1), V(0, 1, 0), 1, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 1, 0), V(0, 0, 0), 1, 0)]),
            new PlaneFace(new P3(V(0, 1, 0), 1), [
                new StraightEdge(new L3(V(0.9999999999999999, 1, 0), V(0, 0, -1)), V(1, 1, 1), V(1, 1, 0), -1, 0),
                new StraightEdge(new L3(V(0, 1, 0), V(1, 0, 0)), V(1, 1, 0), V(0, 1, 0), 1, 0),
                new StraightEdge(new L3(V(0, 1, 0), V(0, 0, 1)), V(0, 1, 0), V(0, 1, 1), 0, 1),
                new StraightEdge(new L3(V(0, 1, 1), V(1, 0, 0)), V(0, 1, 1), V(1, 1, 1), 0, 1)]),
            new PlaneFace(new P3(V(1, 0, 0), 0.3333333333333333), [
                new StraightEdge(new L3(V(0.3333333333333333, 0.33333333333333337, 0), V(0, 0, -1)), V(0.3333333333333333, 0.3333333333333333, 1), V(0.3333333333333333, 0.3333333333333333, 0), -1, 0),
                new StraightEdge(new L3(V(0.3333333333333333, 0, 0), V(0, -1, 0)), V(0.3333333333333333, 0.3333333333333333, 0), V(0.3333333333333333, 0.6666666666666666, 0), -0.3333333333333333, -0.6666666666666666),
                new StraightEdge(new L3(V(0.3333333333333333, 0.6666666666666666, 0), V(0, 0, 1)), V(0.3333333333333333, 0.6666666666666666, 0), V(0.3333333333333333, 0.6666666666666666, 1), 0, 1),
                new StraightEdge(new L3(V(0.3333333333333333, 0, 1), V(0, 1, 0)), V(0.3333333333333333, 0.6666666666666666, 1), V(0.3333333333333333, 0.3333333333333333, 1), 0.6666666666666666, 0.3333333333333333)]),
            new PlaneFace(new P3(V(0, -1, 0), -0.6666666666666666), [
                new StraightEdge(new L3(V(0.6666666666666665, 0.6666666666666666, 0), V(0, 0, 1)), V(0.6666666666666666, 0.6666666666666666, 0), V(0.6666666666666666, 0.6666666666666666, 1), 0, 1),
                new StraightEdge(new L3(V(0, 0.6666666666666666, 1), V(1, 0, 0)), V(0.6666666666666666, 0.6666666666666666, 1), V(0.3333333333333333, 0.6666666666666666, 1), 0.6666666666666666, 0.3333333333333333),
                new StraightEdge(new L3(V(0.3333333333333333, 0.6666666666666666, 0), V(0, 0, 1)), V(0.3333333333333333, 0.6666666666666666, 1), V(0.3333333333333333, 0.6666666666666666, 0), 1, 0),
                new StraightEdge(new L3(V(0, 0.6666666666666666, 0), V(-1, 0, 0)), V(0.3333333333333333, 0.6666666666666666, 0), V(0.6666666666666666, 0.6666666666666666, 0), -0.3333333333333333, -0.6666666666666666)]),
            new PlaneFace(new P3(V(0.7071067811865476, -0.7071067811865475, 0), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(-0.7071067811865475, -0.7071067811865476, 0)), V(0.6666666666666666, 0.6666666666666666, 0), V(1, 1, 0), -0.9428090415820634, -1.414213562373095),
                new StraightEdge(new L3(V(0.9999999999999999, 1, 0), V(0, 0, -1)), V(1, 1, 0), V(1, 1, 1), 0, -1),
                new StraightEdge(new L3(V(0, 0, 1), V(0.7071067811865475, 0.7071067811865476, 0)), V(1, 1, 1), V(0.6666666666666666, 0.6666666666666666, 1), 1.414213562373095, 0.9428090415820634),
                new StraightEdge(new L3(V(0.6666666666666665, 0.6666666666666666, 0), V(0, 0, 1)), V(0.6666666666666666, 0.6666666666666666, 1), V(0.6666666666666666, 0.6666666666666666, 0), 1, 0)]),
            new PlaneFace(new P3(V(0.7071067811865476, -0.7071067811865475, 0), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(-0.7071067811865475, -0.7071067811865476, 0)), V(0, 0, 0), V(0.3333333333333333, 0.3333333333333333, 0), 0, -0.4714045207910317),
                new StraightEdge(new L3(V(0.3333333333333333, 0.33333333333333337, 0), V(0, 0, -1)), V(0.3333333333333333, 0.3333333333333333, 0), V(0.3333333333333333, 0.3333333333333333, 1), 0, -1),
                new StraightEdge(new L3(V(0, 0, 1), V(0.7071067811865475, 0.7071067811865476, 0)), V(0.3333333333333333, 0.3333333333333333, 1), V(0, 0, 1), 0.4714045207910317, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 1), V(0, 0, 0), 1, 0)])], false)
        b2Equal(assert, box, box2, box.and(box2), result)
    },

    'B2T.sphere(5) - B2T.box(12,2,3).translate(-6, 0,1)'(assert) {
        const a = B2T.sphere(5)
        const b = B2T.box(12,2,3).translate(-6, 0,1).flipped()
        const result = new B2([
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(5, 0, 0), V(0, 5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 2, 0), V(0, 0, 4.58257569495584), V(4.58257569495584, 0, 0), 0, 3.141592653589793), V(2.2360679774997863, 2, 4), V(4.472135954999578, 2, 1), 0.5097396788315063, 1.3508083493994372, null, V(4.000000000000001, 0, -2.236067977499787), V(1, 0, -4.47213595499958)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 1), V(-4.898979485566356, 0, 0), V(0, 4.898979485566356, 0), 0, 3.141592653589793), V(4.472135954999578, 2, 1), V(4.898979485566356, 0, 1.000000000000001), 2.721058318305828, 3.141592653589793, null, V(1.9999999999999998, -4.47213595499958, 0), V(5.999519546087385e-16, -4.898979485566356, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0)), V(4.898979485566356, 0, 1.000000000000001), V(0, 0, -5), 1.7721542475852277, 0, null, V(1.000000000000001, 0, -4.898979485566356), V(-5, 0, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(-5, 6.123233995736766e-16, 0), 0, 3.141592653589793), V(0, 0, -5), V(-4.898979485566356, 5.999519546087386e-16, 1.000000000000001), 0, 1.7721542475852277, null, V(-5, 6.123233995736766e-16, 0), V(1.000000000000001, -1.2246467991473547e-16, 4.898979485566356)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 1), V(-4.898979485566356, 0, 0), V(0, 4.898979485566356, 0), 0, 3.141592653589793), V(-4.898979485566356, 5.999519546087386e-16, 1.000000000000001), V(-4.472135954999578, 2, 1), 1.2246467991473535e-16, 0.4205343352839653, null, V(5.999519546087386e-16, 4.898979485566356, 0), V(2.0000000000000004, 4.472135954999579, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 2, 0), V(0, 0, -4.58257569495584), V(-4.58257569495584, 0, 0), 0, 3.141592653589793), V(-4.472135954999578, 2, 1), V(-2.236067977499787, 2, 4), 1.7907843041903562, 2.631852974758287, null, V(1.0000000000000004, 0, 4.47213595499958), V(4.000000000000001, 0, 2.2360679774997876)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 4), V(2.999999999999999, 0, 0), V(0, 2.999999999999999, 0), 0, 3.141592653589793), V(-2.236067977499787, 2, 4), V(-3, 3.67394039744206e-16, 4), 2.411864997362826, 3.141592653589793, null, V(-2.000000000000001, -2.2360679774997876, 0), V(-3.6739403974420584e-16, -2.999999999999999, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(-5, 6.123233995736766e-16, 0), 0, 3.141592653589793), V(-3, 3.67394039744206e-16, 4), V(0, 0, 5), 2.498091544796509, 3.141592653589793, null, V(4, -4.898587196589413e-16, 3), V(5, -6.123233995736766e-16, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0)), V(0, 0, 5), V(3, 0, 4), 3.141592653589793, 2.498091544796509, null, V(5, 0, 0), V(4, 0, -3)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 4), V(2.999999999999999, 0, 0), V(0, 2.999999999999999, 0), 0, 3.141592653589793), V(3, 0, 4), V(2.2360679774997863, 2, 4), 0, 0.7297276562269671, null, V(0, 2.999999999999999, 0), V(-2.000000000000001, 2.2360679774997876, 0))]),
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(-5, 6.123233995736766e-16, 0), V(-6.123233995736766e-16, -5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0), 0, 3.141592653589793), V(4.898979485566356, 0, 1), V(2.9999999999999987, 0, 4), 1.7721542475852277, 2.498091544796509, null, V(-1.000000000000001, 0, 4.898979485566356), V(-4, 0, 3)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0)), V(3, 0, 4), V(0, 0, 5), 2.498091544796509, 3.141592653589793, null, V(-4, 0, 3), V(-5, 0, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(-5, 6.123233995736766e-16, 0), 0, 3.141592653589793), V(0, 0, 5), V(-3, 3.67394039744206e-16, 4), 3.141592653589793, 2.498091544796509, null, V(-5, 6.123233995736766e-16, 0), V(-4, 4.898587196589413e-16, -3)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, 5), V(-5, 0, 0), 0, 3.141592653589793), V(-2.9999999999999982, 0, 4), V(-4.898979485566356, 0, 1), 0.643501108793284, 1.3694384060045657, null, V(-4.000000000000001, 0, -2.9999999999999987), V(-1.0000000000000007, 0, -4.898979485566356)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(-5, 6.123233995736766e-16, 0), 0, 3.141592653589793), V(-4.898979485566356, 5.999519546087386e-16, 1.000000000000001), V(0, 0, -5), 1.7721542475852277, 0, null, V(-1.000000000000001, 1.2246467991473547e-16, -4.898979485566356), V(5, -6.123233995736766e-16, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0)), V(0, 0, -5), V(4.898979485566356, 0, 1.000000000000001), 0, 1.7721542475852277, null, V(5, 0, 0), V(-1.000000000000001, 0, 4.898979485566356))]),
            new PlaneFace(new P3(V(0, -1, 0), -2), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 2, 0), V(0, 0, 4.58257569495584), V(4.58257569495584, 0, 0), 0, 3.141592653589793), V(4.472135954999578, 2, 1), V(2.2360679774997863, 2, 4), 1.3508083493994372, 0.5097396788315063, null, V(-1, 0, 4.47213595499958), V(-4.000000000000001, 0, 2.236067977499787)),
                new StraightEdge(new L3(V(-6, 2, 4), V(1, 0, 0)), V(2.2360679774997863, 2, 4), V(-2.236067977499787, 2, 4), 8.236067977499786, 3.763932022500213),
                new PCurveEdge(new SemiEllipseCurve(V(0, 2, 0), V(0, 0, -4.58257569495584), V(-4.58257569495584, 0, 0), 0, 3.141592653589793), V(-2.236067977499787, 2, 4), V(-4.472135954999578, 2, 1), 2.631852974758287, 1.7907843041903562, null, V(-4.000000000000001, 0, -2.2360679774997876), V(-1.0000000000000004, 0, -4.47213595499958)),
                new StraightEdge(new L3(V(-6, 2, 1), V(1, 0, 0)), V(-4.472135954999578, 2, 1), V(4.472135954999578, 2, 1), 1.5278640450004222, 10.472135954999578)]),
            new PlaneFace(new P3(V(0, 1, 0), 0), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0), 0, 3.141592653589793), V(2.9999999999999987, 0, 4), V(4.898979485566356, 0, 1), 2.498091544796509, 1.7721542475852277, null, V(4, 0, -3), V(1.000000000000001, 0, -4.898979485566356)),
                new StraightEdge(new L3(V(6, 0, 1), V(-1, 0, 0)), V(4.898979485566356, 0, 1), V(-4.898979485566356, 0, 1), 1.1010205144336442, 10.898979485566356),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, 5), V(-5, 0, 0), 0, 3.141592653589793), V(-4.898979485566356, 0, 1), V(-2.9999999999999982, 0, 4), 1.3694384060045657, 0.643501108793284, null, V(1.0000000000000007, 0, 4.898979485566356), V(4.000000000000001, 0, 2.9999999999999987)),
                new StraightEdge(new L3(V(6, 0, 4), V(-1, 0, 0)), V(-2.9999999999999982, 0, 4), V(2.9999999999999987, 0, 4), 8.999999999999998, 3.0000000000000013)]),
            new PlaneFace(new P3(V(0, 0, 1), 1), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 1), V(-4.898979485566356, 0, 0), V(0, 4.898979485566356, 0), 0, 3.141592653589793), V(-4.472135954999578, 2, 1), V(-4.898979485566356, 5.999519546087386e-16, 1.000000000000001), 0.4205343352839653, 1.2246467991473535e-16, null, V(-2.0000000000000004, -4.472135954999579, 0), V(-5.999519546087386e-16, -4.898979485566356, 0)),
                new StraightEdge(new L3(V(6, 0, 1), V(-1, 0, 0)), V(-4.898979485566356, 0, 1), V(4.898979485566356, 0, 1), 10.898979485566356, 1.1010205144336442),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 1), V(-4.898979485566356, 0, 0), V(0, 4.898979485566356, 0), 0, 3.141592653589793), V(4.898979485566356, 0, 1.000000000000001), V(4.472135954999578, 2, 1), 3.141592653589793, 2.721058318305828, null, V(-5.999519546087385e-16, 4.898979485566356, 0), V(-1.9999999999999998, 4.47213595499958, 0)),
                new StraightEdge(new L3(V(-6, 2, 1), V(1, 0, 0)), V(4.472135954999578, 2, 1), V(-4.472135954999578, 2, 1), 10.472135954999578, 1.5278640450004222)]),
            new PlaneFace(new P3(V(0, 0, -1), -4), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 4), V(2.999999999999999, 0, 0), V(0, 2.999999999999999, 0), 0, 3.141592653589793), V(2.2360679774997863, 2, 4), V(3, 0, 4), 0.7297276562269671, 0, null, V(2.000000000000001, -2.2360679774997876, 0), V(0, -2.999999999999999, 0)),
                new StraightEdge(new L3(V(6, 0, 4), V(-1, 0, 0)), V(2.9999999999999987, 0, 4), V(-2.9999999999999982, 0, 4), 3.0000000000000013, 8.999999999999998),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 4), V(2.999999999999999, 0, 0), V(0, 2.999999999999999, 0), 0, 3.141592653589793), V(-3, 3.67394039744206e-16, 4), V(-2.236067977499787, 2, 4), 3.141592653589793, 2.411864997362826, null, V(3.6739403974420584e-16, 2.999999999999999, 0), V(2.000000000000001, 2.2360679774997876, 0)),
                new StraightEdge(new L3(V(-6, 2, 4), V(1, 0, 0)), V(-2.236067977499787, 2, 4), V(2.2360679774997863, 2, 4), 3.763932022500213, 8.236067977499786)])], undefined)
        b2Equal(assert, a, b, a.and(b), result)
    },

    'B2T.sphere(5) AND B2T.octahedron().scale(6)'(assert) {
        const a = B2T.sphere(5)
        const b = B2T.octahedron().scale(6)
        const result = new B2([
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(5, 0, 0), V(0, 5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(2, 2.0000000000000018, 2), V(1.4719601443879746, -2.943920288775947, 1.4719601443879746), V(-2.5495097567963914, 0, 2.5495097567963914), 0.8238977337258584, 3.141592653589793), V(1.1291713066130307, 0, 4.87082869338697), V(0, 1.1291713066130296, 4.87082869338697), 0.8238977337258595, 1.270497368667335, null, V(-2.8121742573035204, 2.1602468994692843, 0.6519273578342348), V(-2.1602468994692874, 2.8121742573035187, -0.6519273578342336)),
                new PCurveEdge(new SemiEllipseCurve(V(-2.0000000000000018, 2, 2), V(1.4719601443879742, 2.9439202887759484, -1.4719601443879733), V(2.549509756796391, 0, 2.5495097567963927), 0, 2.3176949198639334), V(0, 1.1291713066130296, 4.87082869338697), V(-1.1291713066130307, 1.3828360263326826e-16, 4.87082869338697), 1.8710952849224574, 2.317694919863933, null, V(-2.1602468994692847, -2.8121742573035204, 0.6519273578342343), V(-2.8121742573035196, -2.160246899469287, -0.6519273578342344)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(-5, 6.123233995736766e-16, 0), 0, 3.141592653589793), V(-1.1291713066130307, 1.3828360263326826e-16, 4.87082869338697), V(0, 0, 5), 2.9137933168918813, 3.141592653589793, null, V(4.87082869338697, -5.965044768551438e-16, 1.1291713066130307), V(5, -6.123233995736766e-16, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0)), V(0, 0, 5), V(1.1291713066130307, 0, 4.87082869338697), 3.141592653589793, 2.9137933168918813, null, V(5, 0, 0), V(4.87082869338697, 0, -1.1291713066130307))], []),
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(5, 0, 0), V(0, 5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(2, 2.0000000000000018, 2), V(1.4719601443879746, -2.943920288775947, 1.4719601443879746), V(-2.5495097567963914, 0, 2.5495097567963914), 0.8238977337258584, 3.141592653589793), V(0, 4.8708286933869696, 1.1291713066130304), V(0.528039855612025, 4.943920288775949, 0.5280398556120257), 2.918292836119055, 3.141592653589793, null, V(2.1602468994692847, 0.6519273578342353, -2.8121742573035204), V(2.5495097567963914, 3.605262558594415e-16, -2.5495097567963914)),
                new PCurveEdge(new SemiEllipseCurve(V(2, 2.0000000000000018, 2), V(-1.4719601443879746, 2.943920288775947, -1.4719601443879746), V(2.5495097567963914, 0, -2.5495097567963914), 0, 2.3176949198639347), V(0.5280398556120254, 4.943920288775949, 0.5280398556120254), V(1.1291713066130296, 4.87082869338697, 0), 0, 0.22329981747073774, null, V(2.5495097567963914, 0, -2.5495097567963914), V(2.8121742573035196, -0.6519273578342337, -2.1602468994692856)),
                new PCurveEdge(new SemiEllipseCurve(V(1.9999999999999998, 2.0000000000000013, -2.000000000000001), V(1.471960144387974, -2.9439202887759475, -1.4719601443879748), V(2.5495097567963922, 0, 2.549509756796391), 0.8238977337258588, 3.141592653589793), V(1.1291713066130296, 4.87082869338697, 0), V(0.5280398556120262, 4.943920288775949, -0.5280398556120257), 2.9182928361190554, 3.141592653589793, null, V(-2.8121742573035204, 0.6519273578342341, -2.1602468994692847), V(-2.5495097567963922, 3.6052625585944157e-16, -2.549509756796391)),
                new PCurveEdge(new SemiEllipseCurve(V(1.9999999999999998, 2.0000000000000013, -2.000000000000001), V(-1.471960144387974, 2.9439202887759475, 1.4719601443879748), V(-2.5495097567963922, 0, -2.549509756796391), 0, 2.3176949198639343), V(0.5280398556120258, 4.943920288775949, -0.528039855612026), V(0, 4.8708286933869696, -1.1291713066130304), 0, 0.22329981747073796, null, V(-2.5495097567963922, 0, -2.549509756796391), V(-2.160246899469286, -0.6519273578342345, -2.8121742573035196)),
                new PCurveEdge(new SemiEllipseCurve(V(-2.0000000000000013, 1.9999999999999998, -2.000000000000001), V(-1.4719601443879737, -2.9439202887759484, -1.4719601443879728), V(2.549509756796391, 0, -2.5495097567963922), 0.8238977337258602, 3.141592653589793), V(0, 4.8708286933869696, -1.1291713066130304), V(-0.5280398556120273, 4.943920288775948, -0.5280398556120284), 2.9182928361190554, 3.141592653589793, null, V(-2.160246899469285, 0.6519273578342343, 2.8121742573035204), V(-2.549509756796391, 3.6052625585944167e-16, 2.5495097567963922)),
                new PCurveEdge(new SemiEllipseCurve(V(-2.0000000000000013, 1.9999999999999998, -2.000000000000001), V(1.4719601443879737, 2.9439202887759484, 1.4719601443879728), V(-2.549509756796391, 0, 2.5495097567963922), 0, 2.317694919863933), V(-0.5280398556120276, 4.943920288775948, -0.528039855612028), V(-1.1291713066130296, 4.87082869338697, 0), 0, 0.22329981747073774, null, V(-2.549509756796391, 0, 2.5495097567963922), V(-2.812174257303519, -0.6519273578342341, 2.1602468994692865)),
                new PCurveEdge(new SemiEllipseCurve(V(-2.0000000000000018, 2, 2), V(-1.4719601443879742, -2.9439202887759484, 1.4719601443879733), V(-2.549509756796391, 0, -2.5495097567963927), 0.8238977337258597, 3.141592653589793), V(-1.1291713066130296, 4.87082869338697, 0), V(-0.5280398556120279, 4.943920288775948, 0.5280398556120264), 2.918292836119056, 3.141592653589793, null, V(2.812174257303519, 0.651927357834233, 2.1602468994692883), V(2.549509756796391, 3.6052625585944167e-16, 2.5495097567963927)),
                new PCurveEdge(new SemiEllipseCurve(V(-2.0000000000000018, 2, 2), V(1.4719601443879742, 2.9439202887759484, -1.4719601443879733), V(2.549509756796391, 0, 2.5495097567963927), 0, 2.3176949198639334), V(-0.5280398556120276, 4.943920288775948, 0.5280398556120267), V(0, 4.8708286933869696, 1.1291713066130304), 0, 0.22329981747073813, null, V(2.549509756796391, 0, 2.5495097567963927), V(2.1602468994692847, -0.6519273578342352, 2.8121742573035213))], []),
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(5, 0, 0), V(0, 5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(2, 2.0000000000000018, 2), V(-1.4719601443879746, 2.943920288775947, -1.4719601443879746), V(2.5495097567963914, 0, -2.5495097567963914), 0, 2.3176949198639347), V(4.8708286933869696, 1.1291713066130304, 0), V(4.8708286933869696, 0, 1.129171306613032), 1.871095284922458, 2.3176949198639343, null, V(0.6519273578342338, -2.8121742573035187, 2.160246899469287), V(-0.6519273578342364, -2.1602468994692834, 2.8121742573035213)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0)), V(4.8708286933869696, 0, 1.129171306613032), V(4.87082869338697, 0, -1.1291713066130316), 1.7985956634928089, 1.3429969900969845, null, V(1.129171306613032, 0, -4.8708286933869696), V(-1.1291713066130316, 0, -4.87082869338697)),
                new PCurveEdge(new SemiEllipseCurve(V(1.9999999999999998, 2.0000000000000013, -2.000000000000001), V(1.471960144387974, -2.9439202887759475, -1.4719601443879748), V(2.5495097567963922, 0, 2.549509756796391), 0.8238977337258588, 3.141592653589793), V(4.87082869338697, 0, -1.1291713066130316), V(4.8708286933869696, 1.1291713066130304, 0), 0.8238977337258596, 1.2704973686673355, null, V(0.6519273578342355, 2.1602468994692847, 2.8121742573035204), V(-0.6519273578342341, 2.8121742573035196, 2.1602468994692865))], []),
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(5, 0, 0), V(0, 5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(-2.0000000000000018, 2, 2), V(-1.4719601443879742, -2.9439202887759484, 1.4719601443879733), V(-2.549509756796391, 0, -2.5495097567963927), 0.8238977337258597, 3.141592653589793), V(-4.87082869338697, 5.965044768551438e-16, 1.129171306613031), V(-4.8708286933869696, 1.1291713066130304, 0), 0.8238977337258594, 1.2704973686673355, null, V(-0.6519273578342353, 2.1602468994692847, -2.8121742573035204), V(0.6519273578342346, 2.8121742573035204, -2.1602468994692856)),
                new PCurveEdge(new SemiEllipseCurve(V(-2.0000000000000013, 1.9999999999999998, -2.000000000000001), V(1.4719601443879737, 2.9439202887759484, 1.4719601443879728), V(-2.549509756796391, 0, 2.5495097567963922), 0, 2.317694919863933), V(-4.8708286933869696, 1.1291713066130304, 0), V(-4.87082869338697, 5.965044768551438e-16, -1.1291713066130304), 1.8710952849224576, 2.3176949198639334, null, V(-0.6519273578342345, -2.8121742573035204, -2.1602468994692847), V(0.6519273578342342, -2.1602468994692856, -2.8121742573035196)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(-5, 6.123233995736766e-16, 0), 0, 3.141592653589793), V(-4.87082869338697, 5.965044768551438e-16, -1.1291713066130316), V(-4.87082869338697, 5.965044768551438e-16, 1.129171306613031), 1.3429969900969845, 1.7985956634928086, null, V(-1.1291713066130316, 1.3828360263326838e-16, 4.87082869338697), V(1.129171306613031, -1.382836026332683e-16, 4.87082869338697))], []),
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(5, 0, 0), V(0, 5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(1.9999999999999998, 2.0000000000000013, -2.000000000000001), V(-1.471960144387974, 2.9439202887759475, 1.4719601443879748), V(-2.5495097567963922, 0, -2.549509756796391), 0, 2.3176949198639343), V(0, 1.1291713066130296, -4.87082869338697), V(1.1291713066130313, 0, -4.87082869338697), 1.8710952849224582, 2.317694919863934, null, V(2.160246899469287, -2.8121742573035187, -0.6519273578342336), V(2.812174257303521, -2.1602468994692843, 0.6519273578342348)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0)), V(1.1291713066130313, 0, -4.87082869338697), V(0, 0, -5), 0.22779933669791208, 0, null, V(-4.87082869338697, 0, -1.1291713066130313), V(-5, 0, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(-5, 6.123233995736766e-16, 0), 0, 3.141592653589793), V(0, 0, -5), V(-1.1291713066130313, 1.3828360263326835e-16, -4.87082869338697), 0, 0.2277993366979121, null, V(-5, 6.123233995736766e-16, 0), V(-4.87082869338697, 5.965044768551438e-16, 1.1291713066130313)),
                new PCurveEdge(new SemiEllipseCurve(V(-2.0000000000000013, 1.9999999999999998, -2.000000000000001), V(-1.4719601443879737, -2.9439202887759484, -1.4719601443879728), V(2.549509756796391, 0, -2.5495097567963922), 0.8238977337258602, 3.141592653589793), V(-1.1291713066130313, 1.3828360263326835e-16, -4.87082869338697), V(0, 1.1291713066130296, -4.87082869338697), 0.8238977337258602, 1.270497368667336, null, V(2.812174257303519, 2.1602468994692865, -0.6519273578342348), V(2.1602468994692843, 2.8121742573035204, 0.6519273578342343))], []),
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(-5, 6.123233995736766e-16, 0), V(-6.123233995736766e-16, -5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(-1.9999999999999996, -2.0000000000000018, 2), V(-1.4719601443879744, 2.9439202887759466, 1.4719601443879757), V(2.5495097567963922, -3.122248963055649e-16, 2.549509756796391), 0.823897733725858, 3.141592653589793), V(-1.1291713066130307, 1.3828360263326826e-16, 4.87082869338697), V(0, -1.1291713066130296, 4.87082869338697), 0.8238977337258592, 1.2704973686673346, null, V(2.8121742573035213, -2.1602468994692834, 0.6519273578342344), V(2.160246899469288, -2.8121742573035178, -0.6519273578342334)),
                new PCurveEdge(new SemiEllipseCurve(V(2.0000000000000018, -1.9999999999999996, 2), V(-1.4719601443879744, -2.9439202887759492, -1.4719601443879724), V(-2.5495097567963905, 3.1222489630556463e-16, 2.549509756796393), 0, 2.317694919863933), V(0, -1.1291713066130296, 4.87082869338697), V(1.1291713066130307, 0, 4.87082869338697), 1.8710952849224571, 2.317694919863933, null, V(2.1602468994692847, 2.8121742573035213, 0.6519273578342338), V(2.812174257303519, 2.1602468994692874, -0.6519273578342353)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0)), V(1.1291713066130307, 0, 4.87082869338697), V(0, 0, 5), 2.9137933168918813, 3.141592653589793, null, V(-4.87082869338697, 0, 1.1291713066130307), V(-5, 0, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(-5, 6.123233995736766e-16, 0), 0, 3.141592653589793), V(0, 0, 5), V(-1.1291713066130307, 1.3828360263326826e-16, 4.87082869338697), 3.141592653589793, 2.9137933168918813, null, V(-5, 6.123233995736766e-16, 0), V(-4.87082869338697, 5.965044768551438e-16, -1.1291713066130307))], []),
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(-5, 6.123233995736766e-16, 0), V(-6.123233995736766e-16, -5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(-1.9999999999999996, -2.0000000000000018, 2), V(-1.4719601443879744, 2.9439202887759466, 1.4719601443879757), V(2.5495097567963922, -3.122248963055649e-16, 2.549509756796391), 0.823897733725858, 3.141592653589793), V(0, -4.8708286933869696, 1.1291713066130304), V(-0.5280398556120248, -4.943920288775948, 0.5280398556120246), 2.918292836119055, 3.141592653589793, null, V(-2.1602468994692856, -0.6519273578342348, -2.8121742573035204), V(-2.5495097567963922, -4.8301359553876594e-17, -2.549509756796391)),
                new PCurveEdge(new SemiEllipseCurve(V(-1.9999999999999996, -2.0000000000000018, 2), V(1.4719601443879744, -2.9439202887759466, -1.4719601443879757), V(-2.5495097567963922, 3.122248963055649e-16, -2.549509756796391), 0, 2.317694919863935), V(-0.5280398556120252, -4.943920288775948, 0.5280398556120243), V(-1.1291713066130296, -4.87082869338697, 0), 0, 0.22329981747073757, null, V(-2.5495097567963922, 3.122248963055649e-16, -2.549509756796391), V(-2.8121742573035204, 0.6519273578342335, -2.160246899469285)),
                new PCurveEdge(new SemiEllipseCurve(V(-1.9999999999999996, -2.0000000000000018, -2.000000000000001), V(-1.4719601443879742, 2.943920288775947, -1.4719601443879755), V(-2.5495097567963922, 3.122248963055649e-16, 2.549509756796391), 0.8238977337258584, 3.141592653589793), V(-1.1291713066130296, -4.87082869338697, 0), V(-0.5280398556120257, -4.943920288775949, -0.528039855612025), 2.9182928361190554, 3.141592653589793, null, V(2.8121742573035204, -0.6519273578342343, -2.1602468994692847), V(2.5495097567963922, -6.727511521650064e-16, -2.549509756796391)),
                new PCurveEdge(new SemiEllipseCurve(V(-1.9999999999999996, -2.0000000000000018, -2.000000000000001), V(1.4719601443879742, -2.943920288775947, 1.4719601443879755), V(2.5495097567963922, -3.122248963055649e-16, -2.549509756796391), 0, 2.3176949198639347), V(-0.5280398556120254, -4.943920288775949, -0.5280398556120254), V(0, -4.8708286933869696, -1.1291713066130304), 0, 0.22329981747073802, null, V(2.5495097567963922, -3.122248963055649e-16, -2.549509756796391), V(2.160246899469286, 0.6519273578342343, -2.81217425730352)),
                new PCurveEdge(new SemiEllipseCurve(V(2.0000000000000018, -1.9999999999999996, -2.000000000000001), V(1.4719601443879744, 2.9439202887759497, -1.4719601443879726), V(-2.549509756796391, 3.1222489630556473e-16, -2.549509756796393), 0.8238977337258606, 3.141592653589793), V(0, -4.8708286933869696, -1.1291713066130304), V(0.528039855612027, -4.943920288775949, -0.5280398556120286), 2.9182928361190554, 3.141592653589793, null, V(2.160246899469285, -0.6519273578342348, 2.8121742573035213), V(2.549509756796391, -6.727511521650066e-16, 2.549509756796393)),
                new PCurveEdge(new SemiEllipseCurve(V(2.0000000000000018, -1.9999999999999996, -2.000000000000001), V(-1.4719601443879744, -2.9439202887759497, 1.4719601443879726), V(2.549509756796391, -3.1222489630556473e-16, 2.549509756796393), 0, 2.3176949198639325), V(0.5280398556120274, -4.943920288775949, -0.5280398556120283), V(1.1291713066130296, -4.87082869338697, 0), 0, 0.2232998174707378, null, V(2.549509756796391, -3.1222489630556473e-16, 2.549509756796393), V(2.812174257303519, 0.6519273578342342, 2.1602468994692874)),
                new PCurveEdge(new SemiEllipseCurve(V(2.0000000000000018, -1.9999999999999996, 2), V(1.4719601443879744, 2.9439202887759492, 1.4719601443879724), V(2.5495097567963905, -3.1222489630556463e-16, -2.549509756796393), 0.8238977337258602, 3.141592653589793), V(1.1291713066130296, -4.87082869338697, 0), V(0.5280398556120277, -4.943920288775949, 0.5280398556120273), 2.9182928361190554, 3.141592653589793, null, V(-2.8121742573035187, -0.6519273578342341, 2.160246899469288), V(-2.5495097567963905, -4.8301359553877186e-17, 2.549509756796393)),
                new PCurveEdge(new SemiEllipseCurve(V(2.0000000000000018, -1.9999999999999996, 2), V(-1.4719601443879744, -2.9439202887759492, -1.4719601443879724), V(-2.5495097567963905, 3.1222489630556463e-16, 2.549509756796393), 0, 2.317694919863933), V(0.5280398556120274, -4.943920288775949, 0.5280398556120276), V(0, -4.8708286933869696, 1.1291713066130304), 0, 0.2232998174707378, null, V(-2.5495097567963905, 3.1222489630556463e-16, 2.549509756796393), V(-2.1602468994692847, 0.6519273578342347, 2.812174257303521))], []),
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(-5, 6.123233995736766e-16, 0), V(-6.123233995736766e-16, -5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(-1.9999999999999996, -2.0000000000000018, 2), V(1.4719601443879744, -2.9439202887759466, -1.4719601443879757), V(-2.5495097567963922, 3.122248963055649e-16, -2.549509756796391), 0, 2.317694919863935), V(-4.8708286933869696, -1.1291713066130304, 0), V(-4.8708286933869696, 5.965044768551438e-16, 1.129171306613032), 1.8710952849224578, 2.317694919863934, null, V(-0.6519273578342342, 2.8121742573035187, 2.1602468994692874), V(0.6519273578342359, 2.1602468994692834, 2.8121742573035213)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(-5, 6.123233995736766e-16, 0), 0, 3.141592653589793), V(-4.87082869338697, 5.965044768551438e-16, 1.129171306613031), V(-4.87082869338697, 5.965044768551438e-16, -1.1291713066130316), 1.7985956634928086, 1.3429969900969845, null, V(-1.129171306613031, 1.382836026332683e-16, -4.87082869338697), V(1.1291713066130316, -1.3828360263326838e-16, -4.87082869338697)),
                new PCurveEdge(new SemiEllipseCurve(V(-1.9999999999999996, -2.0000000000000018, -2.000000000000001), V(-1.4719601443879742, 2.943920288775947, -1.4719601443879755), V(-2.5495097567963922, 3.122248963055649e-16, 2.549509756796391), 0.8238977337258584, 3.141592653589793), V(-4.87082869338697, 5.965044768551438e-16, -1.1291713066130316), V(-4.8708286933869696, -1.1291713066130304, 0), 0.8238977337258597, 1.2704973686673355, null, V(-0.6519273578342353, -2.1602468994692847, 2.8121742573035204), V(0.6519273578342343, -2.812174257303519, 2.160246899469287))], []),
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(-5, 6.123233995736766e-16, 0), V(-6.123233995736766e-16, -5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(2.0000000000000018, -1.9999999999999996, 2), V(1.4719601443879744, 2.9439202887759492, 1.4719601443879724), V(2.5495097567963905, -3.1222489630556463e-16, -2.549509756796393), 0.8238977337258602, 3.141592653589793), V(4.8708286933869696, 0, 1.129171306613032), V(4.8708286933869696, -1.1291713066130304, 0), 0.8238977337258591, 1.2704973686673355, null, V(0.6519273578342355, -2.1602468994692847, -2.8121742573035204), V(-0.6519273578342349, -2.8121742573035213, -2.1602468994692847)),
                new PCurveEdge(new SemiEllipseCurve(V(2.0000000000000018, -1.9999999999999996, -2.000000000000001), V(-1.4719601443879744, -2.9439202887759497, 1.4719601443879726), V(2.549509756796391, -3.1222489630556473e-16, 2.549509756796393), 0, 2.3176949198639325), V(4.8708286933869696, -1.1291713066130304, 0), V(4.87082869338697, 0, -1.1291713066130316), 1.8710952849224576, 2.317694919863934, null, V(0.6519273578342352, 2.8121742573035218, -2.1602468994692847), V(-0.651927357834235, 2.1602468994692856, -2.8121742573035204)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0)), V(4.87082869338697, 0, -1.1291713066130316), V(4.8708286933869696, 0, 1.129171306613032), 1.3429969900969845, 1.7985956634928089, null, V(1.1291713066130316, 0, 4.87082869338697), V(-1.129171306613032, 0, 4.8708286933869696))], []),
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(-5, 6.123233995736766e-16, 0), V(-6.123233995736766e-16, -5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(-1.9999999999999996, -2.0000000000000018, -2.000000000000001), V(1.4719601443879742, -2.943920288775947, 1.4719601443879755), V(2.5495097567963922, -3.122248963055649e-16, -2.549509756796391), 0, 2.3176949198639347), V(0, -1.1291713066130296, -4.87082869338697), V(-1.1291713066130307, 1.3828360263326826e-16, -4.87082869338697), 1.8710952849224587, 2.3176949198639343, null, V(-2.160246899469288, 2.812174257303518, -0.6519273578342328), V(-2.8121742573035213, 2.1602468994692834, 0.6519273578342355)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(-5, 6.123233995736766e-16, 0), 0, 3.141592653589793), V(-1.1291713066130313, 1.3828360263326835e-16, -4.87082869338697), V(0, 0, -5), 0.2277993366979121, 0, null, V(4.87082869338697, -5.965044768551438e-16, -1.1291713066130313), V(5, -6.123233995736766e-16, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0)), V(0, 0, -5), V(1.1291713066130313, 0, -4.87082869338697), 0, 0.22779933669791208, null, V(5, 0, 0), V(4.87082869338697, 0, 1.1291713066130313)),
                new PCurveEdge(new SemiEllipseCurve(V(2.0000000000000018, -1.9999999999999996, -2.000000000000001), V(1.4719601443879744, 2.9439202887759497, -1.4719601443879726), V(-2.549509756796391, 3.1222489630556473e-16, -2.549509756796393), 0.8238977337258606, 3.141592653589793), V(1.1291713066130313, 0, -4.87082869338697), V(0, -1.1291713066130296, -4.87082869338697), 0.8238977337258604, 1.2704973686673362, null, V(-2.812174257303519, -2.160246899469288, -0.6519273578342348), V(-2.1602468994692847, -2.812174257303522, 0.6519273578342342))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0.5773502691896257, 0.5773502691896257, 0.5773502691896257), 3.4641016151377553)), [
                new PCurveEdge(new SemiEllipseCurve(V(2, 2.0000000000000018, 2), V(1.4719601443879746, -2.943920288775947, 1.4719601443879746), V(-2.5495097567963914, 0, 2.5495097567963914), 0.8238977337258584, 3.141592653589793), V(0, 1.1291713066130296, 4.87082869338697), V(1.1291713066130307, 0, 4.87082869338697), 1.270497368667335, 0.8238977337258595, null, V(2.1602468994692874, -2.8121742573035187, 0.6519273578342336), V(2.8121742573035204, -2.1602468994692843, -0.6519273578342348)),
                new StraightEdge(new L3(V(6, 0, 0), V(-0.7071067811865475, 0, 0.7071067811865475)), V(1.1291713066130296, 0, 4.87082869338697), V(4.8708286933869696, 0, 1.1291713066130304), 6.8883919981838755, 1.5968893760546963),
                new PCurveEdge(new SemiEllipseCurve(V(2, 2.0000000000000018, 2), V(-1.4719601443879746, 2.943920288775947, -1.4719601443879746), V(2.5495097567963914, 0, -2.5495097567963914), 0, 2.3176949198639347), V(4.8708286933869696, 0, 1.129171306613032), V(4.8708286933869696, 1.1291713066130304, 0), 2.3176949198639343, 1.871095284922458, null, V(0.6519273578342364, 2.1602468994692834, -2.8121742573035213), V(-0.6519273578342338, 2.8121742573035187, -2.160246899469287)),
                new StraightEdge(new L3(V(6, 0, 0), V(-0.7071067811865475, 0.7071067811865475, 0)), V(4.8708286933869696, 1.1291713066130304, 0), V(1.1291713066130296, 4.87082869338697, 0), 1.5968893760546963, 6.8883919981838755),
                new PCurveEdge(new SemiEllipseCurve(V(2, 2.0000000000000018, 2), V(-1.4719601443879746, 2.943920288775947, -1.4719601443879746), V(2.5495097567963914, 0, -2.5495097567963914), 0, 2.3176949198639347), V(1.1291713066130296, 4.87082869338697, 0), V(0.5280398556120254, 4.943920288775949, 0.5280398556120254), 0.22329981747073774, 0, null, V(-2.8121742573035196, 0.6519273578342337, 2.1602468994692856), V(-2.5495097567963914, 0, 2.5495097567963914)),
                new PCurveEdge(new SemiEllipseCurve(V(2, 2.0000000000000018, 2), V(1.4719601443879746, -2.943920288775947, 1.4719601443879746), V(-2.5495097567963914, 0, 2.5495097567963914), 0.8238977337258584, 3.141592653589793), V(0.528039855612025, 4.943920288775949, 0.5280398556120257), V(0, 4.8708286933869696, 1.1291713066130304), 3.141592653589793, 2.918292836119055, null, V(-2.5495097567963914, -3.605262558594415e-16, 2.5495097567963914), V(-2.1602468994692847, -0.6519273578342353, 2.8121742573035204)),
                new StraightEdge(new L3(V(0, 6, 0), V(0, -0.7071067811865475, 0.7071067811865475)), V(0, 4.8708286933869696, 1.1291713066130304), V(0, 1.1291713066130296, 4.87082869338697), 1.5968893760546963, 6.8883919981838755)], []),
            new PlaneFace(new PlaneSurface(new P3(V(-0.5773502691896257, 0.5773502691896257, 0.5773502691896257), 3.4641016151377553)), [
                new PCurveEdge(new SemiEllipseCurve(V(-2.0000000000000018, 2, 2), V(-1.4719601443879742, -2.9439202887759484, 1.4719601443879733), V(-2.549509756796391, 0, -2.5495097567963927), 0.8238977337258597, 3.141592653589793), V(-4.8708286933869696, 1.1291713066130304, 0), V(-4.87082869338697, 5.965044768551438e-16, 1.129171306613031), 1.2704973686673355, 0.8238977337258594, null, V(-0.6519273578342346, -2.8121742573035204, 2.1602468994692856), V(0.6519273578342353, -2.1602468994692847, 2.8121742573035204)),
                new StraightEdge(new L3(V(-6, 0, 0), V(0.7071067811865475, 0, 0.7071067811865475)), V(-4.8708286933869696, 0, 1.1291713066130304), V(-1.1291713066130296, 0, 4.87082869338697), 1.5968893760546963, 6.8883919981838755),
                new PCurveEdge(new SemiEllipseCurve(V(-2.0000000000000018, 2, 2), V(1.4719601443879742, 2.9439202887759484, -1.4719601443879733), V(2.549509756796391, 0, 2.5495097567963927), 0, 2.3176949198639334), V(-1.1291713066130307, 1.3828360263326826e-16, 4.87082869338697), V(0, 1.1291713066130296, 4.87082869338697), 2.317694919863933, 1.8710952849224574, null, V(2.8121742573035196, 2.160246899469287, 0.6519273578342344), V(2.1602468994692847, 2.8121742573035204, -0.6519273578342343)),
                new StraightEdge(new L3(V(0, 6, 0), V(0, -0.7071067811865475, 0.7071067811865475)), V(0, 1.1291713066130296, 4.87082869338697), V(0, 4.8708286933869696, 1.1291713066130304), 6.8883919981838755, 1.5968893760546963),
                new PCurveEdge(new SemiEllipseCurve(V(-2.0000000000000018, 2, 2), V(1.4719601443879742, 2.9439202887759484, -1.4719601443879733), V(2.549509756796391, 0, 2.5495097567963927), 0, 2.3176949198639334), V(0, 4.8708286933869696, 1.1291713066130304), V(-0.5280398556120276, 4.943920288775948, 0.5280398556120267), 0.22329981747073813, 0, null, V(-2.1602468994692847, 0.6519273578342352, -2.8121742573035213), V(-2.549509756796391, 0, -2.5495097567963927)),
                new PCurveEdge(new SemiEllipseCurve(V(-2.0000000000000018, 2, 2), V(-1.4719601443879742, -2.9439202887759484, 1.4719601443879733), V(-2.549509756796391, 0, -2.5495097567963927), 0.8238977337258597, 3.141592653589793), V(-0.5280398556120279, 4.943920288775948, 0.5280398556120264), V(-1.1291713066130296, 4.87082869338697, 0), 3.141592653589793, 2.918292836119056, null, V(-2.549509756796391, -3.6052625585944167e-16, -2.5495097567963927), V(-2.812174257303519, -0.651927357834233, -2.1602468994692883)),
                new StraightEdge(new L3(V(-6, 0, 0), V(0.7071067811865475, 0.7071067811865475, 0)), V(-1.1291713066130296, 4.87082869338697, 0), V(-4.8708286933869696, 1.1291713066130304, 0), 6.8883919981838755, 1.5968893760546963)], []),
            new PlaneFace(new PlaneSurface(new P3(V(-0.5773502691896257, -0.5773502691896257, 0.5773502691896257), 3.4641016151377553)), [
                new PCurveEdge(new SemiEllipseCurve(V(-1.9999999999999996, -2.0000000000000018, 2), V(-1.4719601443879744, 2.9439202887759466, 1.4719601443879757), V(2.5495097567963922, -3.122248963055649e-16, 2.549509756796391), 0.823897733725858, 3.141592653589793), V(0, -1.1291713066130296, 4.87082869338697), V(-1.1291713066130307, 1.3828360263326826e-16, 4.87082869338697), 1.2704973686673346, 0.8238977337258592, null, V(-2.160246899469288, 2.8121742573035178, 0.6519273578342334), V(-2.8121742573035213, 2.1602468994692834, -0.6519273578342344)),
                new StraightEdge(new L3(V(-6, 0, 0), V(0.7071067811865475, 0, 0.7071067811865475)), V(-1.1291713066130296, 0, 4.87082869338697), V(-4.8708286933869696, 0, 1.1291713066130304), 6.8883919981838755, 1.5968893760546963),
                new PCurveEdge(new SemiEllipseCurve(V(-1.9999999999999996, -2.0000000000000018, 2), V(1.4719601443879744, -2.9439202887759466, -1.4719601443879757), V(-2.5495097567963922, 3.122248963055649e-16, -2.549509756796391), 0, 2.317694919863935), V(-4.8708286933869696, 5.965044768551438e-16, 1.129171306613032), V(-4.8708286933869696, -1.1291713066130304, 0), 2.317694919863934, 1.8710952849224578, null, V(-0.6519273578342359, -2.1602468994692834, -2.8121742573035213), V(0.6519273578342342, -2.8121742573035187, -2.1602468994692874)),
                new StraightEdge(new L3(V(-6, 0, 0), V(0.7071067811865475, -0.7071067811865475, 0)), V(-4.8708286933869696, -1.1291713066130304, 0), V(-1.1291713066130296, -4.87082869338697, 0), 1.5968893760546963, 6.8883919981838755),
                new PCurveEdge(new SemiEllipseCurve(V(-1.9999999999999996, -2.0000000000000018, 2), V(1.4719601443879744, -2.9439202887759466, -1.4719601443879757), V(-2.5495097567963922, 3.122248963055649e-16, -2.549509756796391), 0, 2.317694919863935), V(-1.1291713066130296, -4.87082869338697, 0), V(-0.5280398556120252, -4.943920288775948, 0.5280398556120243), 0.22329981747073757, 0, null, V(2.8121742573035204, -0.6519273578342335, 2.160246899469285), V(2.5495097567963922, -3.122248963055649e-16, 2.549509756796391)),
                new PCurveEdge(new SemiEllipseCurve(V(-1.9999999999999996, -2.0000000000000018, 2), V(-1.4719601443879744, 2.9439202887759466, 1.4719601443879757), V(2.5495097567963922, -3.122248963055649e-16, 2.549509756796391), 0.823897733725858, 3.141592653589793), V(-0.5280398556120248, -4.943920288775948, 0.5280398556120246), V(0, -4.8708286933869696, 1.1291713066130304), 3.141592653589793, 2.918292836119055, null, V(2.5495097567963922, 4.8301359553876594e-17, 2.549509756796391), V(2.1602468994692856, 0.6519273578342348, 2.8121742573035204)),
                new StraightEdge(new L3(V(0, -6, 0), V(0, 0.7071067811865475, 0.7071067811865475)), V(0, -4.8708286933869696, 1.1291713066130304), V(0, -1.1291713066130296, 4.87082869338697), 1.5968893760546963, 6.8883919981838755)], []),
            new PlaneFace(new PlaneSurface(new P3(V(0.5773502691896257, -0.5773502691896257, 0.5773502691896257), 3.4641016151377553)), [
                new PCurveEdge(new SemiEllipseCurve(V(2.0000000000000018, -1.9999999999999996, 2), V(1.4719601443879744, 2.9439202887759492, 1.4719601443879724), V(2.5495097567963905, -3.1222489630556463e-16, -2.549509756796393), 0.8238977337258602, 3.141592653589793), V(4.8708286933869696, -1.1291713066130304, 0), V(4.8708286933869696, 0, 1.129171306613032), 1.2704973686673355, 0.8238977337258591, null, V(0.6519273578342349, 2.8121742573035213, 2.1602468994692847), V(-0.6519273578342355, 2.1602468994692847, 2.8121742573035204)),
                new StraightEdge(new L3(V(6, 0, 0), V(-0.7071067811865475, 0, 0.7071067811865475)), V(4.8708286933869696, 0, 1.1291713066130304), V(1.1291713066130296, 0, 4.87082869338697), 1.5968893760546963, 6.8883919981838755),
                new PCurveEdge(new SemiEllipseCurve(V(2.0000000000000018, -1.9999999999999996, 2), V(-1.4719601443879744, -2.9439202887759492, -1.4719601443879724), V(-2.5495097567963905, 3.1222489630556463e-16, 2.549509756796393), 0, 2.317694919863933), V(1.1291713066130307, 0, 4.87082869338697), V(0, -1.1291713066130296, 4.87082869338697), 2.317694919863933, 1.8710952849224571, null, V(-2.812174257303519, -2.1602468994692874, 0.6519273578342353), V(-2.1602468994692847, -2.8121742573035213, -0.6519273578342338)),
                new StraightEdge(new L3(V(0, -6, 0), V(0, 0.7071067811865475, 0.7071067811865475)), V(0, -1.1291713066130296, 4.87082869338697), V(0, -4.8708286933869696, 1.1291713066130304), 6.8883919981838755, 1.5968893760546963),
                new PCurveEdge(new SemiEllipseCurve(V(2.0000000000000018, -1.9999999999999996, 2), V(-1.4719601443879744, -2.9439202887759492, -1.4719601443879724), V(-2.5495097567963905, 3.1222489630556463e-16, 2.549509756796393), 0, 2.317694919863933), V(0, -4.8708286933869696, 1.1291713066130304), V(0.5280398556120274, -4.943920288775949, 0.5280398556120276), 0.2232998174707378, 0, null, V(2.1602468994692847, -0.6519273578342347, -2.812174257303521), V(2.5495097567963905, -3.1222489630556463e-16, -2.549509756796393)),
                new PCurveEdge(new SemiEllipseCurve(V(2.0000000000000018, -1.9999999999999996, 2), V(1.4719601443879744, 2.9439202887759492, 1.4719601443879724), V(2.5495097567963905, -3.1222489630556463e-16, -2.549509756796393), 0.8238977337258602, 3.141592653589793), V(0.5280398556120277, -4.943920288775949, 0.5280398556120273), V(1.1291713066130296, -4.87082869338697, 0), 3.141592653589793, 2.9182928361190554, null, V(2.5495097567963905, 4.8301359553877186e-17, -2.549509756796393), V(2.8121742573035187, 0.6519273578342341, -2.160246899469288)),
                new StraightEdge(new L3(V(6, 0, 0), V(-0.7071067811865475, -0.7071067811865475, 0)), V(1.1291713066130296, -4.87082869338697, 0), V(4.8708286933869696, -1.1291713066130304, 0), 6.8883919981838755, 1.5968893760546963)], []),
            new PlaneFace(new PlaneSurface(new P3(V(0.5773502691896257, 0.5773502691896257, -0.5773502691896257), 3.4641016151377553)), [
                new PCurveEdge(new SemiEllipseCurve(V(1.9999999999999998, 2.0000000000000013, -2.000000000000001), V(1.471960144387974, -2.9439202887759475, -1.4719601443879748), V(2.5495097567963922, 0, 2.549509756796391), 0.8238977337258588, 3.141592653589793), V(4.8708286933869696, 1.1291713066130304, 0), V(4.87082869338697, 0, -1.1291713066130316), 1.2704973686673355, 0.8238977337258596, null, V(0.6519273578342341, -2.8121742573035196, -2.1602468994692865), V(-0.6519273578342355, -2.1602468994692847, -2.8121742573035204)),
                new StraightEdge(new L3(V(6, 0, 0), V(-0.7071067811865475, 0, -0.7071067811865475)), V(4.8708286933869696, 0, -1.1291713066130304), V(1.1291713066130296, 0, -4.87082869338697), 1.5968893760546963, 6.8883919981838755),
                new PCurveEdge(new SemiEllipseCurve(V(1.9999999999999998, 2.0000000000000013, -2.000000000000001), V(-1.471960144387974, 2.9439202887759475, 1.4719601443879748), V(-2.5495097567963922, 0, -2.549509756796391), 0, 2.3176949198639343), V(1.1291713066130313, 0, -4.87082869338697), V(0, 1.1291713066130296, -4.87082869338697), 2.317694919863934, 1.8710952849224582, null, V(-2.812174257303521, 2.1602468994692843, -0.6519273578342348), V(-2.160246899469287, 2.8121742573035187, 0.6519273578342336)),
                new StraightEdge(new L3(V(0, 6, 0), V(0, -0.7071067811865475, -0.7071067811865475)), V(0, 1.1291713066130296, -4.87082869338697), V(0, 4.8708286933869696, -1.1291713066130304), 6.8883919981838755, 1.5968893760546963),
                new PCurveEdge(new SemiEllipseCurve(V(1.9999999999999998, 2.0000000000000013, -2.000000000000001), V(-1.471960144387974, 2.9439202887759475, 1.4719601443879748), V(-2.5495097567963922, 0, -2.549509756796391), 0, 2.3176949198639343), V(0, 4.8708286933869696, -1.1291713066130304), V(0.5280398556120258, 4.943920288775949, -0.528039855612026), 0.22329981747073796, 0, null, V(2.160246899469286, 0.6519273578342345, 2.8121742573035196), V(2.5495097567963922, 0, 2.549509756796391)),
                new PCurveEdge(new SemiEllipseCurve(V(1.9999999999999998, 2.0000000000000013, -2.000000000000001), V(1.471960144387974, -2.9439202887759475, -1.4719601443879748), V(2.5495097567963922, 0, 2.549509756796391), 0.8238977337258588, 3.141592653589793), V(0.5280398556120262, 4.943920288775949, -0.5280398556120257), V(1.1291713066130296, 4.87082869338697, 0), 3.141592653589793, 2.9182928361190554, null, V(2.5495097567963922, -3.6052625585944157e-16, 2.549509756796391), V(2.8121742573035204, -0.6519273578342341, 2.1602468994692847)),
                new StraightEdge(new L3(V(6, 0, 0), V(-0.7071067811865475, 0.7071067811865475, 0)), V(1.1291713066130296, 4.87082869338697, 0), V(4.8708286933869696, 1.1291713066130304, 0), 6.8883919981838755, 1.5968893760546963)], []),
            new PlaneFace(new PlaneSurface(new P3(V(-0.5773502691896257, 0.5773502691896257, -0.5773502691896257), 3.4641016151377553)), [
                new PCurveEdge(new SemiEllipseCurve(V(-2.0000000000000013, 1.9999999999999998, -2.000000000000001), V(-1.4719601443879737, -2.9439202887759484, -1.4719601443879728), V(2.549509756796391, 0, -2.5495097567963922), 0.8238977337258602, 3.141592653589793), V(0, 1.1291713066130296, -4.87082869338697), V(-1.1291713066130313, 1.3828360263326835e-16, -4.87082869338697), 1.270497368667336, 0.8238977337258602, null, V(-2.1602468994692843, -2.8121742573035204, -0.6519273578342343), V(-2.812174257303519, -2.1602468994692865, 0.6519273578342348)),
                new StraightEdge(new L3(V(-6, 0, 0), V(0.7071067811865475, 0, -0.7071067811865475)), V(-1.1291713066130296, 0, -4.87082869338697), V(-4.8708286933869696, 0, -1.1291713066130304), 6.8883919981838755, 1.5968893760546963),
                new PCurveEdge(new SemiEllipseCurve(V(-2.0000000000000013, 1.9999999999999998, -2.000000000000001), V(1.4719601443879737, 2.9439202887759484, 1.4719601443879728), V(-2.549509756796391, 0, 2.5495097567963922), 0, 2.317694919863933), V(-4.87082869338697, 5.965044768551438e-16, -1.1291713066130304), V(-4.8708286933869696, 1.1291713066130304, 0), 2.3176949198639334, 1.8710952849224576, null, V(-0.6519273578342342, 2.1602468994692856, 2.8121742573035196), V(0.6519273578342345, 2.8121742573035204, 2.1602468994692847)),
                new StraightEdge(new L3(V(-6, 0, 0), V(0.7071067811865475, 0.7071067811865475, 0)), V(-4.8708286933869696, 1.1291713066130304, 0), V(-1.1291713066130296, 4.87082869338697, 0), 1.5968893760546963, 6.8883919981838755),
                new PCurveEdge(new SemiEllipseCurve(V(-2.0000000000000013, 1.9999999999999998, -2.000000000000001), V(1.4719601443879737, 2.9439202887759484, 1.4719601443879728), V(-2.549509756796391, 0, 2.5495097567963922), 0, 2.317694919863933), V(-1.1291713066130296, 4.87082869338697, 0), V(-0.5280398556120276, 4.943920288775948, -0.528039855612028), 0.22329981747073774, 0, null, V(2.812174257303519, 0.6519273578342341, -2.1602468994692865), V(2.549509756796391, 0, -2.5495097567963922)),
                new PCurveEdge(new SemiEllipseCurve(V(-2.0000000000000013, 1.9999999999999998, -2.000000000000001), V(-1.4719601443879737, -2.9439202887759484, -1.4719601443879728), V(2.549509756796391, 0, -2.5495097567963922), 0.8238977337258602, 3.141592653589793), V(-0.5280398556120273, 4.943920288775948, -0.5280398556120284), V(0, 4.8708286933869696, -1.1291713066130304), 3.141592653589793, 2.9182928361190554, null, V(2.549509756796391, -3.6052625585944167e-16, -2.5495097567963922), V(2.160246899469285, -0.6519273578342343, -2.8121742573035204)),
                new StraightEdge(new L3(V(0, 6, 0), V(0, -0.7071067811865475, -0.7071067811865475)), V(0, 4.8708286933869696, -1.1291713066130304), V(0, 1.1291713066130296, -4.87082869338697), 1.5968893760546963, 6.8883919981838755)], []),
            new PlaneFace(new PlaneSurface(new P3(V(-0.5773502691896257, -0.5773502691896257, -0.5773502691896257), 3.4641016151377553)), [
                new PCurveEdge(new SemiEllipseCurve(V(-1.9999999999999996, -2.0000000000000018, -2.000000000000001), V(-1.4719601443879742, 2.943920288775947, -1.4719601443879755), V(-2.5495097567963922, 3.122248963055649e-16, 2.549509756796391), 0.8238977337258584, 3.141592653589793), V(-4.8708286933869696, -1.1291713066130304, 0), V(-4.87082869338697, 5.965044768551438e-16, -1.1291713066130316), 1.2704973686673355, 0.8238977337258597, null, V(-0.6519273578342343, 2.812174257303519, -2.160246899469287), V(0.6519273578342353, 2.1602468994692847, -2.8121742573035204)),
                new StraightEdge(new L3(V(-6, 0, 0), V(0.7071067811865475, 0, -0.7071067811865475)), V(-4.8708286933869696, 0, -1.1291713066130304), V(-1.1291713066130296, 0, -4.87082869338697), 1.5968893760546963, 6.8883919981838755),
                new PCurveEdge(new SemiEllipseCurve(V(-1.9999999999999996, -2.0000000000000018, -2.000000000000001), V(1.4719601443879742, -2.943920288775947, 1.4719601443879755), V(2.5495097567963922, -3.122248963055649e-16, -2.549509756796391), 0, 2.3176949198639347), V(-1.1291713066130307, 1.3828360263326826e-16, -4.87082869338697), V(0, -1.1291713066130296, -4.87082869338697), 2.3176949198639343, 1.8710952849224587, null, V(2.8121742573035213, -2.1602468994692834, -0.6519273578342355), V(2.160246899469288, -2.812174257303518, 0.6519273578342328)),
                new StraightEdge(new L3(V(0, -6, 0), V(0, 0.7071067811865475, -0.7071067811865475)), V(0, -1.1291713066130296, -4.87082869338697), V(0, -4.8708286933869696, -1.1291713066130304), 6.8883919981838755, 1.5968893760546963),
                new PCurveEdge(new SemiEllipseCurve(V(-1.9999999999999996, -2.0000000000000018, -2.000000000000001), V(1.4719601443879742, -2.943920288775947, 1.4719601443879755), V(2.5495097567963922, -3.122248963055649e-16, -2.549509756796391), 0, 2.3176949198639347), V(0, -4.8708286933869696, -1.1291713066130304), V(-0.5280398556120254, -4.943920288775949, -0.5280398556120254), 0.22329981747073802, 0, null, V(-2.160246899469286, -0.6519273578342343, 2.81217425730352), V(-2.5495097567963922, 3.122248963055649e-16, 2.549509756796391)),
                new PCurveEdge(new SemiEllipseCurve(V(-1.9999999999999996, -2.0000000000000018, -2.000000000000001), V(-1.4719601443879742, 2.943920288775947, -1.4719601443879755), V(-2.5495097567963922, 3.122248963055649e-16, 2.549509756796391), 0.8238977337258584, 3.141592653589793), V(-0.5280398556120257, -4.943920288775949, -0.528039855612025), V(-1.1291713066130296, -4.87082869338697, 0), 3.141592653589793, 2.9182928361190554, null, V(-2.5495097567963922, 6.727511521650064e-16, 2.549509756796391), V(-2.8121742573035204, 0.6519273578342343, 2.1602468994692847)),
                new StraightEdge(new L3(V(-6, 0, 0), V(0.7071067811865475, -0.7071067811865475, 0)), V(-1.1291713066130296, -4.87082869338697, 0), V(-4.8708286933869696, -1.1291713066130304, 0), 6.8883919981838755, 1.5968893760546963)], []),
            new PlaneFace(new PlaneSurface(new P3(V(0.5773502691896257, -0.5773502691896257, -0.5773502691896257), 3.4641016151377553)), [
                new PCurveEdge(new SemiEllipseCurve(V(2.0000000000000018, -1.9999999999999996, -2.000000000000001), V(1.4719601443879744, 2.9439202887759497, -1.4719601443879726), V(-2.549509756796391, 3.1222489630556473e-16, -2.549509756796393), 0.8238977337258606, 3.141592653589793), V(0, -1.1291713066130296, -4.87082869338697), V(1.1291713066130313, 0, -4.87082869338697), 1.2704973686673362, 0.8238977337258604, null, V(2.1602468994692847, 2.812174257303522, -0.6519273578342342), V(2.812174257303519, 2.160246899469288, 0.6519273578342348)),
                new StraightEdge(new L3(V(6, 0, 0), V(-0.7071067811865475, 0, -0.7071067811865475)), V(1.1291713066130296, 0, -4.87082869338697), V(4.8708286933869696, 0, -1.1291713066130304), 6.8883919981838755, 1.5968893760546963),
                new PCurveEdge(new SemiEllipseCurve(V(2.0000000000000018, -1.9999999999999996, -2.000000000000001), V(-1.4719601443879744, -2.9439202887759497, 1.4719601443879726), V(2.549509756796391, -3.1222489630556473e-16, 2.549509756796393), 0, 2.3176949198639325), V(4.87082869338697, 0, -1.1291713066130316), V(4.8708286933869696, -1.1291713066130304, 0), 2.317694919863934, 1.8710952849224576, null, V(0.651927357834235, -2.1602468994692856, 2.8121742573035204), V(-0.6519273578342352, -2.8121742573035218, 2.1602468994692847)),
                new StraightEdge(new L3(V(6, 0, 0), V(-0.7071067811865475, -0.7071067811865475, 0)), V(4.8708286933869696, -1.1291713066130304, 0), V(1.1291713066130296, -4.87082869338697, 0), 1.5968893760546963, 6.8883919981838755),
                new PCurveEdge(new SemiEllipseCurve(V(2.0000000000000018, -1.9999999999999996, -2.000000000000001), V(-1.4719601443879744, -2.9439202887759497, 1.4719601443879726), V(2.549509756796391, -3.1222489630556473e-16, 2.549509756796393), 0, 2.3176949198639325), V(1.1291713066130296, -4.87082869338697, 0), V(0.5280398556120274, -4.943920288775949, -0.5280398556120283), 0.2232998174707378, 0, null, V(-2.812174257303519, -0.6519273578342342, -2.1602468994692874), V(-2.549509756796391, 3.1222489630556473e-16, -2.549509756796393)),
                new PCurveEdge(new SemiEllipseCurve(V(2.0000000000000018, -1.9999999999999996, -2.000000000000001), V(1.4719601443879744, 2.9439202887759497, -1.4719601443879726), V(-2.549509756796391, 3.1222489630556473e-16, -2.549509756796393), 0.8238977337258606, 3.141592653589793), V(0.528039855612027, -4.943920288775949, -0.5280398556120286), V(0, -4.8708286933869696, -1.1291713066130304), 3.141592653589793, 2.9182928361190554, null, V(-2.549509756796391, 6.727511521650066e-16, -2.549509756796393), V(-2.160246899469285, 0.6519273578342348, -2.8121742573035213)),
                new StraightEdge(new L3(V(0, -6, 0), V(0, 0.7071067811865475, -0.7071067811865475)), V(0, -4.8708286933869696, -1.1291713066130304), V(0, -1.1291713066130296, -4.87082869338697), 1.5968893760546963, 6.8883919981838755)], [])], undefined)
        b2Equal(assert, a, b, a.and(b), result)
    },

    'B2T.sphere(5) - B2T.box(12,2,3).translate(-6, 0,1) 2'(assert) {
        const a = B2T.sphere(5)
        const b = B2T.box(12,2,2).translate(-6, -1,1).rotateX(-10*DEG).flipped()
        const result = new B2([
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(5, 0, 0), V(0, 5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0.984807753012208, -0.1736481776669303), V(0, -0.8506988600962733, -4.824552979233506), V(-4.898979485566356, 0, 0), 0, 3.141592653589793), V(-4.795831523312719, 1.1584559306791384, 0.8111595753452777), V(-3.8729833462074166, 1.505752286012999, 2.7807750813696934), 1.7763652579560707, 2.229854362621306, null, V(1.000000000000001, 0.832787404420872, 4.722972066298714), V(3.0000000000000004, 0.6725365002032875, 3.8141440266322277)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0.520944533000792, 2.954423259036625), V(0, 3.939231012048832, -0.6945927106677224), V(-3.9999999999999996, 0, 0), 0, 1.703430096372499), V(-3.8729833462074166, 1.505752286012999, 2.7807750813696934), V(-3.964868114183388, 4.855563045076089e-16, 3.0462798356572343), 1.3181160716528182, 1.703430096372499, null, V(-0.999999999999999, -3.814144026632229, 0.6725365002032889), V(0.5289809421253958, -3.904632858518692, 0.6884921227176649)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(-5, 6.123233995736766e-16, 0), 0, 3.141592653589793), V(-3.964868114183388, 4.855563045076089e-16, 3.0462798356572343), V(0, 0, 5), 2.2259182971449905, 3.141592653589793, null, V(3.0462798356572343, -3.730616850044757e-16, 3.964868114183388), V(5, -6.123233995736766e-16, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0)), V(0, 0, 5), V(3.964868114183388, 0, 3.0462798356572343), 3.141592653589793, 2.2259182971449905, null, V(5, 0, 0), V(3.0462798356572343, 0, -3.964868114183388)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0.520944533000792, 2.954423259036625), V(0, -3.939231012048832, 0.6945927106677224), V(3.9999999999999996, 0, 0), 1.438162557217294, 3.141592653589793), V(3.964868114183388, 0, 3.0462798356572343), V(3.872983346207416, 1.505752286012999, 2.7807750813696934), 1.4381625572172942, 1.8234765819369751, null, V(0.5289809421253954, 3.904632858518692, -0.6884921227176649), V(-0.9999999999999994, 3.814144026632229, -0.6725365002032889)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0.984807753012208, -0.1736481776669303), V(0, 0.8506988600962733, 4.824552979233506), V(4.898979485566356, 0, 0), 0, 3.141592653589793), V(3.872983346207416, 1.505752286012999, 2.7807750813696934), V(4.795831523312719, 1.1584559306791384, 0.8111595753452777), 0.9117382909684875, 1.3652273956337226, null, V(3.0000000000000004, -0.6725365002032876, -3.814144026632228), V(1.0000000000000004, -0.832787404420872, -4.722972066298714)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0.17364817766693008, 0.9848077530122084), V(0, 4.8245529792335065, -0.8506988600962722), V(4.898979485566356, 0, 0), 0, 1.606796696749177), V(4.795831523312719, 1.1584559306791384, 0.8111595753452777), V(4.895805224462492, 0, 1.015426611885745), 1.3652273956337226, 1.6067966967491771, null, V(1.0000000000000004, -4.7229720662987145, 0.8327874044208708), V(-0.1763269807084656, -4.821426942288336, 0.8501476554401471)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0)), V(4.895805224462492, 0, 1.015426611885745), V(0, 0, -5), 1.775304209444905, 0, null, V(1.015426611885745, 0, -4.895805224462492), V(-5, 0, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(-5, 6.123233995736766e-16, 0), 0, 3.141592653589793), V(0, 0, -5), V(-4.895805224462492, 5.995632197386881e-16, 1.015426611885745), 0, 1.775304209444905, null, V(-5, 6.123233995736766e-16, 0), V(1.015426611885745, -1.243538950014919e-16, 4.895805224462492)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0.17364817766693008, 0.9848077530122084), V(0, -4.8245529792335065, 0.8506988600962722), V(-4.898979485566356, 0, 0), 1.5347959568406162, 3.141592653589793), V(-4.895805224462492, 5.995632197386881e-16, 1.015426611885745), V(-4.795831523312719, 1.1584559306791384, 0.8111595753452777), 1.5347959568406164, 1.7763652579560707, null, V(-0.176326980708464, 4.821426942288336, -0.8501476554401471), V(1.000000000000001, 4.7229720662987145, -0.8327874044208708))]),
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(-5, 6.123233995736766e-16, 0), V(-6.123233995736766e-16, -5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(1.8173711827954756e-17, -0.984807753012208, 0.17364817766693041), V(3.959914808152817e-15, 0.8506988600962738, 4.824552979233506), V(-4.898979485566356, 5.999519546087386e-16, 3.915215160272503e-15), 0, 3.141592653589793), V(-3.872983346207416, -0.463863220011417, 3.1280714367035545), V(-4.795831523312719, -0.8111595753452777, 1.1584559306791384), 0.9117382909684884, 1.3652273956337235, null, V(-3, -0.6725365002032881, -3.8141440266322286), V(-1, -0.8327874044208725, -4.722972066298714)),
                new PCurveEdge(new SemiEllipseCurve(V(-1.3428701023778413e-17, 0.17364817766693008, 0.9848077530122084), V(4.898979485566356, -5.999519546087386e-16, 1.7258951700871463e-16), V(-6.20807191446493e-16, -4.8245529792335065, 0.8506988600962722), 0.03600036995428033, 3.105592283635513), V(-4.795831523312719, -0.8111595753452777, 1.1584559306791384), V(-4.895805224462492, 5.995632197386881e-16, 1.015426611885745), 2.936023722428619, 3.105592283635513, null, V(-1.0000000000000002, 4.7229720662987145, -0.8327874044208708), V(-0.1763269807084637, 4.821426942288336, -0.8501476554401471)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(-5, 6.123233995736766e-16, 0), 0, 3.141592653589793), V(-4.895805224462492, 5.995632197386881e-16, 1.015426611885745), V(0, 0, -5), 1.775304209444905, 0, null, V(-1.015426611885745, 1.243538950014919e-16, -4.895805224462492), V(5, -6.123233995736766e-16, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0)), V(0, 0, -5), V(4.895805224462492, 0, 1.015426611885745), 0, 1.775304209444905, null, V(5, 0, 0), V(-1.015426611885745, 0, 4.895805224462492)),
                new PCurveEdge(new SemiEllipseCurve(V(-1.3428701023778413e-17, 0.17364817766693008, 0.9848077530122084), V(4.898979485566356, -5.999519546087386e-16, 1.7258951700871463e-16), V(-6.20807191446493e-16, -4.8245529792335065, 0.8506988600962722), 0.03600036995428033, 3.105592283635513), V(4.895805224462492, 0, 1.015426611885745), V(4.795831523312719, -0.8111595753452777, 1.1584559306791384), 0.036000369954280194, 0.2055689311611737, null, V(-0.17632698070846461, -4.821426942288336, 0.8501476554401471), V(-0.9999999999999997, -4.7229720662987145, 0.8327874044208708)),
                new PCurveEdge(new SemiEllipseCurve(V(1.8173711827954756e-17, -0.984807753012208, 0.17364817766693041), V(-3.959914808152817e-15, -0.8506988600962738, -4.824552979233506), V(4.898979485566356, -5.999519546087386e-16, -3.915215160272503e-15), 0, 3.141592653589793), V(4.795831523312719, -0.8111595753452777, 1.1584559306791384), V(3.8729833462074166, -0.463863220011417, 3.1280714367035545), 1.7763652579560714, 2.2298543626213068, null, V(-1.0000000000000004, 0.8327874044208724, 4.722972066298714), V(-3.0000000000000013, 0.6725365002032878, 3.8141440266322273)),
                new PCurveEdge(new SemiEllipseCurve(V(6.379730548727326e-17, 0.520944533000792, 2.954423259036625), V(-3.9999999999999996, 4.898587196589413e-16, 0), V(-4.824166650007591e-16, -3.939231012048832, 0.6945927106677224), 0.13263376957760245, 3.008958884012191), V(3.8729833462074166, -0.463863220011417, 3.1280714367035545), V(3.964868114183388, 0, 3.0462798356572343), 2.8889123984477143, 3.008958884012191, null, V(1.0000000000000013, 3.8141440266322286, -0.6725365002032889), V(0.528980942125396, 3.904632858518692, -0.6884921227176649)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0)), V(3.964868114183388, 0, 3.0462798356572343), V(0, 0, 5), 2.2259182971449905, 3.141592653589793, null, V(-3.0462798356572343, 0, 3.964868114183388), V(-5, 0, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(-5, 6.123233995736766e-16, 0), 0, 3.141592653589793), V(0, 0, 5), V(-3.964868114183388, 4.855563045076089e-16, 3.0462798356572343), 3.141592653589793, 2.2259182971449905, null, V(-5, 6.123233995736766e-16, 0), V(-3.0462798356572343, 3.730616850044757e-16, -3.964868114183388)),
                new PCurveEdge(new SemiEllipseCurve(V(6.379730548727326e-17, 0.520944533000792, 2.954423259036625), V(-3.9999999999999996, 4.898587196589413e-16, 0), V(-4.824166650007591e-16, -3.939231012048832, 0.6945927106677224), 0.13263376957760245, 3.008958884012191), V(-3.964868114183388, 4.855563045076089e-16, 3.0462798356572343), V(-3.872983346207416, -0.463863220011417, 3.1280714367035545), 0.13263376957760237, 0.25268025514207904, null, V(0.5289809421253949, -3.904632858518692, 0.6884921227176649), V(1.0000000000000009, -3.8141440266322286, 0.6725365002032889))]),
            new PlaneFace(new P3(V(0, -0.9848077530122081, 0.1736481776669303), -1), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0.984807753012208, -0.1736481776669303), V(0, -0.8506988600962733, -4.824552979233506), V(-4.898979485566356, 0, 0), 0, 3.141592653589793), V(-3.8729833462074166, 1.505752286012999, 2.7807750813696934), V(-4.795831523312719, 1.1584559306791384, 0.8111595753452777), 2.229854362621306, 1.7763652579560707, null, V(-3.0000000000000004, -0.6725365002032875, -3.8141440266322277), V(-1.000000000000001, -0.832787404420872, -4.722972066298714)),
                new StraightEdge(new L3(V(-6, 1.1584559306791384, 0.8111595753452777), V(1, 0, 0)), V(-4.795831523312719, 1.1584559306791384, 0.8111595753452777), V(4.795831523312719, 1.1584559306791384, 0.8111595753452777), 1.2041684766872809, 10.79583152331272),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0.984807753012208, -0.1736481776669303), V(0, 0.8506988600962733, 4.824552979233506), V(4.898979485566356, 0, 0), 0, 3.141592653589793), V(4.795831523312719, 1.1584559306791384, 0.8111595753452777), V(3.872983346207416, 1.505752286012999, 2.7807750813696934), 1.3652273956337226, 0.9117382909684875, null, V(-1.0000000000000004, 0.832787404420872, 4.722972066298714), V(-3.0000000000000004, 0.6725365002032876, 3.814144026632228)),
                new StraightEdge(new L3(V(-6, 1.505752286012999, 2.7807750813696934), V(1, 0, 0)), V(3.872983346207416, 1.505752286012999, 2.7807750813696934), V(-3.8729833462074166, 1.505752286012999, 2.7807750813696934), 9.872983346207416, 2.1270166537925834)]),
            new PlaneFace(new P3(V(0, 0.9848077530122081, -0.17364817766693044), -1), [
                new PCurveEdge(new SemiEllipseCurve(V(1.8173711827954756e-17, -0.984807753012208, 0.17364817766693041), V(3.959914808152817e-15, 0.8506988600962738, 4.824552979233506), V(-4.898979485566356, 5.999519546087386e-16, 3.915215160272503e-15), 0, 3.141592653589793), V(-4.795831523312719, -0.8111595753452777, 1.1584559306791384), V(-3.872983346207416, -0.463863220011417, 3.1280714367035545), 1.3652273956337235, 0.9117382909684884, null, V(1, 0.8327874044208725, 4.722972066298714), V(3, 0.6725365002032881, 3.8141440266322286)),
                new StraightEdge(new L3(V(6, -0.463863220011417, 3.1280714367035545), V(-1, 0, 0)), V(-3.872983346207416, -0.463863220011417, 3.1280714367035545), V(3.8729833462074166, -0.463863220011417, 3.1280714367035545), 9.872983346207416, 2.1270166537925834),
                new PCurveEdge(new SemiEllipseCurve(V(1.8173711827954756e-17, -0.984807753012208, 0.17364817766693041), V(-3.959914808152817e-15, -0.8506988600962738, -4.824552979233506), V(4.898979485566356, -5.999519546087386e-16, -3.915215160272503e-15), 0, 3.141592653589793), V(3.8729833462074166, -0.463863220011417, 3.1280714367035545), V(4.795831523312719, -0.8111595753452777, 1.1584559306791384), 2.2298543626213068, 1.7763652579560714, null, V(3.0000000000000013, -0.6725365002032878, -3.8141440266322273), V(1.0000000000000004, -0.8327874044208724, -4.722972066298714)),
                new StraightEdge(new L3(V(6, -0.8111595753452777, 1.1584559306791384), V(-1, 0, 0)), V(4.795831523312719, -0.8111595753452777, 1.1584559306791384), V(-4.795831523312719, -0.8111595753452777, 1.1584559306791384), 1.2041684766872809, 10.79583152331272)]),
            new PlaneFace(new P3(V(0, 0.1736481776669303, 0.9848077530122081), 1), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0.17364817766693008, 0.9848077530122084), V(0, -4.8245529792335065, 0.8506988600962722), V(-4.898979485566356, 0, 0), 1.5347959568406162, 3.141592653589793), V(-4.795831523312719, 1.1584559306791384, 0.8111595753452777), V(-4.895805224462492, 5.995632197386881e-16, 1.015426611885745), 1.7763652579560707, 1.5347959568406164, null, V(-1.000000000000001, -4.7229720662987145, 0.8327874044208708), V(0.176326980708464, -4.821426942288336, 0.8501476554401471)),
                new PCurveEdge(new SemiEllipseCurve(V(-1.3428701023778413e-17, 0.17364817766693008, 0.9848077530122084), V(4.898979485566356, -5.999519546087386e-16, 1.7258951700871463e-16), V(-6.20807191446493e-16, -4.8245529792335065, 0.8506988600962722), 0.03600036995428033, 3.105592283635513), V(-4.895805224462492, 5.995632197386881e-16, 1.015426611885745), V(-4.795831523312719, -0.8111595753452777, 1.1584559306791384), 3.105592283635513, 2.936023722428619, null, V(0.1763269807084637, -4.821426942288336, 0.8501476554401471), V(1.0000000000000002, -4.7229720662987145, 0.8327874044208708)),
                new StraightEdge(new L3(V(6, -0.8111595753452777, 1.1584559306791384), V(-1, 0, 0)), V(-4.795831523312719, -0.8111595753452777, 1.1584559306791384), V(4.795831523312719, -0.8111595753452777, 1.1584559306791384), 10.79583152331272, 1.2041684766872809),
                new PCurveEdge(new SemiEllipseCurve(V(-1.3428701023778413e-17, 0.17364817766693008, 0.9848077530122084), V(4.898979485566356, -5.999519546087386e-16, 1.7258951700871463e-16), V(-6.20807191446493e-16, -4.8245529792335065, 0.8506988600962722), 0.03600036995428033, 3.105592283635513), V(4.795831523312719, -0.8111595753452777, 1.1584559306791384), V(4.895805224462492, 0, 1.015426611885745), 0.2055689311611737, 0.036000369954280194, null, V(0.9999999999999997, 4.7229720662987145, -0.8327874044208708), V(0.17632698070846461, 4.821426942288336, -0.8501476554401471)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0.17364817766693008, 0.9848077530122084), V(0, 4.8245529792335065, -0.8506988600962722), V(4.898979485566356, 0, 0), 0, 1.606796696749177), V(4.895805224462492, 0, 1.015426611885745), V(4.795831523312719, 1.1584559306791384, 0.8111595753452777), 1.6067966967491771, 1.3652273956337226, null, V(0.1763269807084656, 4.821426942288336, -0.8501476554401471), V(-1.0000000000000004, 4.7229720662987145, -0.8327874044208708)),
                new StraightEdge(new L3(V(-6, 1.1584559306791384, 0.8111595753452777), V(1, 0, 0)), V(4.795831523312719, 1.1584559306791384, 0.8111595753452777), V(-4.795831523312719, 1.1584559306791384, 0.8111595753452777), 10.79583152331272, 1.2041684766872809)]),
            new PlaneFace(new P3(V(0, -0.17364817766693055, -0.9848077530122081), -3), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0.520944533000792, 2.954423259036625), V(0, -3.939231012048832, 0.6945927106677224), V(3.9999999999999996, 0, 0), 1.438162557217294, 3.141592653589793), V(3.872983346207416, 1.505752286012999, 2.7807750813696934), V(3.964868114183388, 0, 3.0462798356572343), 1.8234765819369751, 1.4381625572172942, null, V(0.9999999999999994, -3.814144026632229, 0.6725365002032889), V(-0.5289809421253954, -3.904632858518692, 0.6884921227176649)),
                new PCurveEdge(new SemiEllipseCurve(V(6.379730548727326e-17, 0.520944533000792, 2.954423259036625), V(-3.9999999999999996, 4.898587196589413e-16, 0), V(-4.824166650007591e-16, -3.939231012048832, 0.6945927106677224), 0.13263376957760245, 3.008958884012191), V(3.964868114183388, 0, 3.0462798356572343), V(3.8729833462074166, -0.463863220011417, 3.1280714367035545), 3.008958884012191, 2.8889123984477143, null, V(-0.528980942125396, -3.904632858518692, 0.6884921227176649), V(-1.0000000000000013, -3.8141440266322286, 0.6725365002032889)),
                new StraightEdge(new L3(V(6, -0.463863220011417, 3.1280714367035545), V(-1, 0, 0)), V(3.8729833462074166, -0.463863220011417, 3.1280714367035545), V(-3.872983346207416, -0.463863220011417, 3.1280714367035545), 2.1270166537925834, 9.872983346207416),
                new PCurveEdge(new SemiEllipseCurve(V(6.379730548727326e-17, 0.520944533000792, 2.954423259036625), V(-3.9999999999999996, 4.898587196589413e-16, 0), V(-4.824166650007591e-16, -3.939231012048832, 0.6945927106677224), 0.13263376957760245, 3.008958884012191), V(-3.872983346207416, -0.463863220011417, 3.1280714367035545), V(-3.964868114183388, 4.855563045076089e-16, 3.0462798356572343), 0.25268025514207904, 0.13263376957760237, null, V(-1.0000000000000009, 3.8141440266322286, -0.6725365002032889), V(-0.5289809421253949, 3.904632858518692, -0.6884921227176649)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0.520944533000792, 2.954423259036625), V(0, 3.939231012048832, -0.6945927106677224), V(-3.9999999999999996, 0, 0), 0, 1.703430096372499), V(-3.964868114183388, 4.855563045076089e-16, 3.0462798356572343), V(-3.8729833462074166, 1.505752286012999, 2.7807750813696934), 1.703430096372499, 1.3181160716528182, null, V(-0.5289809421253958, 3.904632858518692, -0.6884921227176649), V(0.999999999999999, 3.814144026632229, -0.6725365002032889)),
                new StraightEdge(new L3(V(-6, 1.505752286012999, 2.7807750813696934), V(1, 0, 0)), V(-3.8729833462074166, 1.505752286012999, 2.7807750813696934), V(3.872983346207416, 1.505752286012999, 2.7807750813696934), 2.1270166537925834, 9.872983346207416)])], undefined)
        b2Equal(assert, a, b, a.and(b), result)
    },

    'B2T.sphere(5) AND B2T.tetrahedron(V3.ZERO, V3.Y, V3.Z, V(7,0,0))'(assert) {
        const a = B2T.sphere(5)
        const b = B2T.tetrahedron(V3.O, V3.Y, V3.Z, V(7,0,0))
        const result = new B2([
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(5, 0, 0), V(0, 5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(5, 0, 0), V(0, 5, 0), 0, 3.141592653589793), V(5, 0, -3.061616997868383e-16), V(4.991762566325767, 0.2868910619534618, 0), 0, 0.057409743118180596, null, V(0, 5, 0), V(-0.2868910619534618, 4.991762566325767, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0.07070707070707062, 0.4949494949494949, 0.4949494949494949), V(-0.492519286346973, 3.5179949024783825, -3.447635004428815), V(4.900505024479568, 0, -0.7000721463542231), 0, 1.7119554405905957), V(4.991762566325767, 0.2868910619534618, 0), V(4.991762566325767, 0, 0.2868910619534626), 1.6299720553223171, 1.711955440590596, null, V(0.20183545628617083, -3.5118371007669364, 3.4830034641546264), V(-0.2018354562861731, -3.483003464154626, 3.5118371007669364)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, 5), V(5, 0, 0), 0, 3.141592653589793), V(4.991762566325767, 0, 0.28689106195346187), V(5, 0, 0), 1.513386583676716, 1.5707963267948966, null, V(0.2868910619534621, 0, -4.991762566325767), V(3.061616997868383e-16, 0, -5))]),
            new PlaneFace(new P3(V(0, 0, -1), 0), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(5, 0, 0), V(0, 5, 0), 0, 3.141592653589793), V(4.991762566325767, 0.2868910619534618, 0), V(5, 0, -3.061616997868383e-16), 0.057409743118180596, 0, null, V(0.2868910619534618, -4.991762566325767, 0), V(0, -5, 0)),
                new StraightEdge(new L3(V(0, 0, 0), V(1, 0, 0)), V(5, 0, 0), V(0, 0, 0), 5, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 0, 0), V(0, 1, 0), 0, 1),
                new StraightEdge(new L3(V(0, 1, 0), V(0.9899494936611665, -0.1414213562373095, 0)), V(0, 1, 0), V(4.991762566325767, 0.2868910619534618, 0), 0, 5.042441658174448)]),
            new PlaneFace(new P3(V(0.10050378152592121, 0.7035264706814485, 0.7035264706814485), 0.7035264706814485), [
                new PCurveEdge(new SemiEllipseCurve(V(0.07070707070707062, 0.4949494949494949, 0.4949494949494949), V(-0.492519286346973, 3.5179949024783825, -3.447635004428815), V(4.900505024479568, 0, -0.7000721463542231), 0, 1.7119554405905957), V(4.991762566325767, 0, 0.2868910619534626), V(4.991762566325767, 0.2868910619534618, 0), 1.711955440590596, 1.6299720553223171, null, V(0.2018354562861731, 3.483003464154626, -3.5118371007669364), V(-0.20183545628617083, 3.5118371007669364, -3.4830034641546264)),
                new StraightEdge(new L3(V(0, 1, 0), V(0.9899494936611665, -0.1414213562373095, 0)), V(4.991762566325767, 0.2868910619534618, 0), V(0, 1, 0), 5.042441658174448, 0),
                new StraightEdge(new L3(V(0, 1, 0), V(0, -0.7071067811865475, 0.7071067811865475)), V(0, 1, 0), V(0, 0, 1), 0, 1.4142135623730951),
                new StraightEdge(new L3(V(7, 0, 0), V(-0.9899494936611665, 0, 0.1414213562373095)), V(0, 0, 1), V(4.991762566325767, 0, 0.28689106195346187), 7.0710678118654755, 2.028626153691028)]),
            new PlaneFace(new P3(V(0, -1, 0), 0), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, 5), V(5, 0, 0), 0, 3.141592653589793), V(5, 0, 0), V(4.991762566325767, 0, 0.28689106195346187), 1.5707963267948966, 1.513386583676716, null, V(-3.061616997868383e-16, 0, 5), V(-0.2868910619534621, 0, 4.991762566325767)),
                new StraightEdge(new L3(V(7, 0, 0), V(-0.9899494936611665, 0, 0.1414213562373095)), V(4.991762566325767, 0, 0.28689106195346187), V(0, 0, 1), 2.028626153691028, 7.0710678118654755),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 1), V(0, 0, 0), 1, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(1, 0, 0)), V(0, 0, 0), V(5, 0, 0), 0, 5)]),
            new PlaneFace(new P3(V(-1, 0, 0), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 0), V(0, 0, 1), 0, 1),
                new StraightEdge(new L3(V(0, 1, 0), V(0, -0.7071067811865475, 0.7071067811865475)), V(0, 0, 1), V(0, 1, 0), 1.4142135623730951, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 1, 0), V(0, 0, 0), 1, 0)])], undefined)
        b2Equal(assert, a, b, a.and(b), result)
    },

    'B2T.sphere(5) AND B2T.tetrahedron(V3.ZERO, V3.X, V3.Y, V3.Z).scale(7)'(assert) {
        const a = B2T.sphere(5)
        const b = B2T.tetrahedron(V3.O, V3.X, V3.Y, V3.Z).scale(7)
        const result = new B2([
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(5, 0, 0), V(0, 5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, 5), V(5, 0, 0), 0, 3.141592653589793), V(0, 0, 5), V(3.0000000000000058, 0, 3.9999999999999942), 0, 0.643501108793286, null, V(5, 0, 0), V(3.999999999999995, 0, -3.0000000000000067)),
                new PCurveEdge(new SemiEllipseCurve(V(2.3333333333333326, 2.3333333333333313, 2.333333333333336), V(1.201850425154661, -2.4037008503093267, 1.2018504251546627), V(-2.081665999466134, 0, 2.081665999466131), 0.24256387409548896, 3.141592653589793), V(3.0000000000000107, 0, 3.999999999999992), V(0, 3.0000000000000053, 3.9999999999999947), 0.24256387409547972, 1.8518312282977136, null, V(-2.3094010767585003, 0.5773502691896129, 1.7320508075688847), V(-0.5773502691896154, 2.309401076758501, -1.732050807568882)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(0, 5, 0), 0, 3.141592653589793), V(0, 3.0000000000000053, 3.9999999999999947), V(0, 0, 5), 2.4980915447965075, 3.141592653589793, null, V(0, -3.9999999999999956, 3.0000000000000053), V(0, -5, 6.123233995736766e-16))]),
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(5, 0, 0), V(0, 5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, 5), V(5, 0, 0), 0, 3.141592653589793), V(3.9999999999999947, 0, 3.0000000000000053), V(5, 0, 0), 0.9272952180016107, 1.5707963267948966, null, V(3.000000000000006, 0, -3.9999999999999956), V(3.061616997868383e-16, 0, -5)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(5, 0, 0), V(0, 5, 0), 0, 3.141592653589793), V(5, 0, -3.061616997868383e-16), V(3.9999999999999947, 3.0000000000000053, 0), 0, 0.6435011087932858, null, V(0, 5, 0), V(-3.0000000000000053, 3.9999999999999956, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(2.3333333333333326, 2.3333333333333313, 2.333333333333336), V(-1.201850425154661, 2.4037008503093267, -1.2018504251546627), V(2.081665999466134, 0, -2.081665999466131), 0, 2.899028779494304), V(3.9999999999999947, 3.0000000000000053, 0), V(3.9999999999999925, 0, 3.0000000000000098), 1.2897614252920795, 2.899028779494311, null, V(1.7320508075688816, -2.309401076758501, 0.5773502691896176), V(-1.7320508075688834, -0.5773502691896192, 2.3094010767584994))]),
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(5, 0, 0), V(0, 5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(5, 0, 0), V(0, 5, 0), 0, 3.141592653589793), V(3.0000000000000058, 3.9999999999999942, 0), V(0, 5, 0), 0.9272952180016106, 1.5707963267948966, null, V(-3.999999999999995, 3.0000000000000067, 0), V(-5, 3.061616997868383e-16, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(0, 5, 0), 0, 3.141592653589793), V(0, 5, 0), V(0, 3.9999999999999942, 3.0000000000000058), 1.5707963267948966, 2.2142974355881826, null, V(0, 3.061616997868383e-16, 5), V(0, -3.0000000000000067, 3.999999999999995)),
                new PCurveEdge(new SemiEllipseCurve(V(2.3333333333333326, 2.3333333333333313, 2.333333333333336), V(1.201850425154661, -2.4037008503093267, 1.2018504251546627), V(-2.081665999466134, 0, 2.081665999466131), 0.24256387409548896, 3.141592653589793), V(0, 3.9999999999999942, 3.0000000000000058), V(1.1314829081786715, 4.737034183642658, 1.1314829081786737), 2.336958976488679, 3.141592653589793, null, V(0.5773502691896236, 1.7320508075688807, -2.3094010767585003), V(2.081665999466134, 2.943684552439088e-16, -2.081665999466131)),
                new PCurveEdge(new SemiEllipseCurve(V(2.3333333333333326, 2.3333333333333313, 2.333333333333336), V(-1.201850425154661, 2.4037008503093267, -1.2018504251546627), V(2.081665999466134, 0, -2.081665999466131), 0, 2.899028779494304), V(1.1314829081786717, 4.737034183642658, 1.1314829081786735), V(3.0000000000000058, 3.9999999999999942, 0), 0, 0.8046336771011146, null, V(2.081665999466134, 0, -2.081665999466131), V(2.309401076758501, -1.7320508075688812, -0.5773502691896201))]),
            new PlaneFace(new P3(V(0, -1, 0), 0), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, 5), V(5, 0, 0), 0, 3.141592653589793), V(3.0000000000000058, 0, 3.9999999999999942), V(0, 0, 5), 0.643501108793286, 0, null, V(-3.999999999999995, 0, 3.0000000000000067), V(-5, 0, 0)),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 5), V(0, 0, 0), 5, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(1, 0, 0)), V(0, 0, 0), V(5, 0, 0), 0, 5),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, 5), V(5, 0, 0), 0, 3.141592653589793), V(5, 0, 0), V(3.9999999999999947, 0, 3.0000000000000053), 1.5707963267948966, 0.9272952180016107, null, V(-3.061616997868383e-16, 0, 5), V(-3.000000000000006, 0, 3.9999999999999956)),
                new StraightEdge(new L3(V(7, 0, 0), V(-0.7071067811865475, 0, 0.7071067811865475)), V(3.9999999999999947, 0, 3.0000000000000053), V(3.0000000000000058, 0, 3.9999999999999942), 4.242640687119293, 5.656854249492373)]),
            new PlaneFace(new P3(V(0, 0, -1), 0), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(5, 0, 0), V(0, 5, 0), 0, 3.141592653589793), V(3.9999999999999947, 3.0000000000000053, 0), V(5, 0, -3.061616997868383e-16), 0.6435011087932858, 0, null, V(3.0000000000000053, -3.9999999999999956, 0), V(0, -5, 0)),
                new StraightEdge(new L3(V(0, 0, 0), V(1, 0, 0)), V(5, 0, 0), V(0, 0, 0), 5, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 0, 0), V(0, 5, 0), 0, 5),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(5, 0, 0), V(0, 5, 0), 0, 3.141592653589793), V(0, 5, 0), V(3.0000000000000058, 3.9999999999999942, 0), 1.5707963267948966, 0.9272952180016106, null, V(5, -3.061616997868383e-16, 0), V(3.999999999999995, -3.0000000000000067, 0)),
                new StraightEdge(new L3(V(7, 0, 0), V(-0.7071067811865475, 0.7071067811865475, 0)), V(3.0000000000000058, 3.9999999999999942, 0), V(3.9999999999999947, 3.0000000000000053, 0), 5.656854249492373, 4.242640687119293)]),
            new PlaneFace(new P3(V(0.5773502691896257, 0.5773502691896257, 0.5773502691896257), 4.0414518843273814), [
                new PCurveEdge(new SemiEllipseCurve(V(2.3333333333333326, 2.3333333333333313, 2.333333333333336), V(1.201850425154661, -2.4037008503093267, 1.2018504251546627), V(-2.081665999466134, 0, 2.081665999466131), 0.24256387409548896, 3.141592653589793), V(0, 3.0000000000000053, 3.9999999999999947), V(3.0000000000000107, 0, 3.999999999999992), 1.8518312282977136, 0.24256387409547972, null, V(0.5773502691896154, -2.309401076758501, 1.732050807568882), V(2.3094010767585003, -0.5773502691896129, -1.7320508075688847)),
                new StraightEdge(new L3(V(7, 0, 0), V(-0.7071067811865475, 0, 0.7071067811865475)), V(3.0000000000000058, 0, 3.9999999999999942), V(3.9999999999999947, 0, 3.0000000000000053), 5.656854249492373, 4.242640687119293),
                new PCurveEdge(new SemiEllipseCurve(V(2.3333333333333326, 2.3333333333333313, 2.333333333333336), V(-1.201850425154661, 2.4037008503093267, -1.2018504251546627), V(2.081665999466134, 0, -2.081665999466131), 0, 2.899028779494304), V(3.9999999999999925, 0, 3.0000000000000098), V(3.9999999999999947, 3.0000000000000053, 0), 2.899028779494311, 1.2897614252920795, null, V(1.7320508075688834, 0.5773502691896192, -2.3094010767584994), V(-1.7320508075688816, 2.309401076758501, -0.5773502691896176)),
                new StraightEdge(new L3(V(7, 0, 0), V(-0.7071067811865475, 0.7071067811865475, 0)), V(3.9999999999999947, 3.0000000000000053, 0), V(3.0000000000000058, 3.9999999999999942, 0), 4.242640687119293, 5.656854249492373),
                new PCurveEdge(new SemiEllipseCurve(V(2.3333333333333326, 2.3333333333333313, 2.333333333333336), V(-1.201850425154661, 2.4037008503093267, -1.2018504251546627), V(2.081665999466134, 0, -2.081665999466131), 0, 2.899028779494304), V(3.0000000000000058, 3.9999999999999942, 0), V(1.1314829081786717, 4.737034183642658, 1.1314829081786735), 0.8046336771011146, 0, null, V(-2.309401076758501, 1.7320508075688812, 0.5773502691896201), V(-2.081665999466134, 0, 2.081665999466131)),
                new PCurveEdge(new SemiEllipseCurve(V(2.3333333333333326, 2.3333333333333313, 2.333333333333336), V(1.201850425154661, -2.4037008503093267, 1.2018504251546627), V(-2.081665999466134, 0, 2.081665999466131), 0.24256387409548896, 3.141592653589793), V(1.1314829081786715, 4.737034183642658, 1.1314829081786737), V(0, 3.9999999999999942, 3.0000000000000058), 3.141592653589793, 2.336958976488679, null, V(-2.081665999466134, -2.943684552439088e-16, 2.081665999466131), V(-0.5773502691896236, -1.7320508075688807, 2.3094010767585003)),
                new StraightEdge(new L3(V(0, 0, 7), V(0, 0.7071067811865475, -0.7071067811865475)), V(0, 3.9999999999999942, 3.0000000000000058), V(0, 3.0000000000000053, 3.9999999999999947), 5.656854249492373, 4.242640687119293)]),
            new PlaneFace(new P3(V(-1, 0, 0), 0), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(0, 5, 0), 0, 3.141592653589793), V(0, 3.9999999999999942, 3.0000000000000058), V(0, 5, 0), 2.2142974355881826, 1.5707963267948966, null, V(0, 3.0000000000000067, -3.999999999999995), V(0, -3.061616997868383e-16, -5)),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 5, 0), V(0, 0, 0), 5, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 0), V(0, 0, 5), 0, 5),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(0, 5, 0), 0, 3.141592653589793), V(0, 0, 5), V(0, 3.0000000000000053, 3.9999999999999947), 3.141592653589793, 2.4980915447965075, null, V(0, 5, -6.123233995736766e-16), V(0, 3.9999999999999956, -3.0000000000000053)),
                new StraightEdge(new L3(V(0, 0, 7), V(0, 0.7071067811865475, -0.7071067811865475)), V(0, 3.0000000000000053, 3.9999999999999947), V(0, 3.9999999999999942, 3.0000000000000058), 4.242640687119293, 5.656854249492373)])], undefined)
        b2Equal(assert, a, b, a.and(b), result)
    },

    'B2T.sphere(5) AND B2T.tetrahedron(V3.ZERO, V3.X, V3.Y, V3.Z.negated()).scale(7)'(assert) {
        const a = B2T.sphere(5)
        const b = B2T.tetrahedron(V(-6,0,0), V(6,0,0), V(0,-4,0), V(0,0,-6))
        const result = new B2([
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(-5, 6.123233995736766e-16, 0), V(-6.123233995736766e-16, -5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0), 0, 3.141592653589793), V(0, 0, -5), V(1.1291713066130296, 0, -4.87082869338697), 0, 0.22779933669791175, null, V(5, 0, 0), V(4.87082869338697, 0, 1.1291713066130296)),
                new PCurveEdge(new SemiEllipseCurve(V(1.411764705882352, -2.1176470588235303, -1.4117647058823533), V(2.0917534572581817, 2.7890046096775745, -2.091753457258183), V(-2.874840149008801, 3.5206637865439285e-16, -2.874840149008799), 0.7085839061491113, 3.141592653589793), V(1.1291713066130291, 0, -4.87082869338697), V(0, -0.6994424542963525, -4.950836318555471), 0.7085839061491117, 1.0373562345961393, null, V(-3.5440484447865424, -1.8149704259460573, -0.8215928058674522), V(-3.262983117260863, -2.401508361946856, 0.3392794256594245)),
                new PCurveEdge(new SemiEllipseCurve(V(-1.4117647058823535, -2.117647058823529, -1.4117647058823533), V(2.0917534572581813, -2.789004609677577, 2.0917534572581813), V(2.8748401490088, -3.520663786543927e-16, -2.8748401490088007), 0, 2.43300874744068), V(0, -0.6994424542963525, -4.950836318555471), V(-1.1291713066130296, 1.382836026332681e-16, -4.87082869338697), 2.1042364189936533, 2.4330087474406805, null, V(-3.2629831172608608, 2.401508361946859, -0.33927942565942426), V(-3.5440484447865406, 1.8149704259460617, 0.8215928058674513)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, 5), V(-5, 0, 0), 0, 3.141592653589793), V(-1.1291713066130296, 0, -4.87082869338697), V(0, 0, -5), 2.9137933168918817, 3.141592653589793, null, V(4.870828693386971, 0, -1.1291713066130287), V(5, 0, -6.123233995736766e-16))], []),
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(-5, 6.123233995736766e-16, 0), V(-6.123233995736766e-16, -5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0), 0, 3.141592653589793), V(4.8708286933869696, 0, -1.1291713066130304), V(5, 0, 0), 1.3429969900969847, 1.5707963267948966, null, V(1.1291713066130304, 0, 4.87082869338697), V(3.061616997868383e-16, 0, 5)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(5, -6.123233995736766e-16, 0), V(-6.123233995736766e-16, -5, 0), 0, 3.141592653589793), V(5, 0, -3.061616997868383e-16), V(4.950836318555471, -0.6994424542963525, 0), 1.2246467991473532e-16, 0.1403487973449601, null, V(-1.2246467991473533e-15, -5, 0), V(-0.6994424542963528, -4.950836318555472, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(1.411764705882352, -2.1176470588235303, -1.4117647058823533), V(-2.0917534572581817, -2.7890046096775745, 2.091753457258183), V(2.874840149008801, -3.5206637865439285e-16, 2.874840149008799), 0, 2.433008747440682), V(4.950836318555471, -0.6994424542963525, 0), V(4.87082869338697, 0, -1.1291713066130293), 2.1042364189936533, 2.4330087474406805, null, V(0.3392794256594245, 2.4015083619468567, -3.262983117260862), V(-0.8215928058674515, 1.81497042594606, -3.544048444786541))], []),
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(-5, 6.123233995736766e-16, 0), V(-6.123233995736766e-16, -5, 0), V(0, 0, 5)), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, 5), V(-5, 0, 0), 0, 3.141592653589793), V(-4.999999999999999, 0, 0), V(-4.8708286933869696, 0, -1.1291713066130304), 1.5707963267948966, 1.7985956634928086, null, V(-3.061616997868383e-16, 0, -5), V(1.129171306613031, 0, -4.87082869338697)),
                new PCurveEdge(new SemiEllipseCurve(V(-1.4117647058823535, -2.117647058823529, -1.4117647058823533), V(-2.0917534572581813, 2.789004609677577, -2.0917534572581813), V(-2.8748401490088, 3.520663786543927e-16, 2.8748401490088007), 0.708583906149113, 3.141592653589793), V(-4.87082869338697, 5.965044768551438e-16, -1.1291713066130304), V(-4.950836318555471, -0.6994424542963525, 0), 0.7085839061491122, 1.03735623459614, null, V(-0.8215928058674524, -1.8149704259460604, 3.544048444786542), V(0.33927942565942515, -2.401508361946859, 3.2629831172608608)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(5, -6.123233995736766e-16, 0), V(-6.123233995736766e-16, -5, 0), 0, 3.141592653589793), V(-4.950836318555471, -0.6994424542963525, 0), V(-5, 6.123233995736766e-16, -3.061616997868383e-16), 3.001243856244833, 3.141592653589793, null, V(-0.6994424542963525, 4.950836318555472, 0), V(0, 5, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 1, 0), 0)), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0), 0, 3.141592653589793), V(1.1291713066130296, 0, -4.87082869338697), V(0, 0, -5), 0.22779933669791175, 0, null, V(-4.87082869338697, 0, -1.1291713066130296), V(-5, 0, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, 5), V(-5, 0, 0), 0, 3.141592653589793), V(0, 0, -5), V(-1.1291713066130296, 0, -4.87082869338697), 3.141592653589793, 2.9137933168918817, null, V(-5, 0, 6.123233995736766e-16), V(-4.870828693386971, 0, 1.1291713066130287)),
                new StraightEdge(new L3(V(-6, 0, 0), V(0.7071067811865475, 0, -0.7071067811865475)), V(-1.1291713066130296, 0, -4.87082869338697), V(-4.8708286933869696, 0, -1.1291713066130304), 6.8883919981838755, 1.5968893760546963),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, 5), V(-5, 0, 0), 0, 3.141592653589793), V(-4.8708286933869696, 0, -1.1291713066130304), V(-4.999999999999999, 0, 0), 1.7985956634928086, 1.5707963267948966, null, V(-1.129171306613031, 0, 4.87082869338697), V(3.061616997868383e-16, 0, 5)),
                new StraightEdge(new L3(V(-6, 0, 0), V(1, 0, 0)), V(-4.999999999999999, 0, 0), V(5, 0, 0), 1.0000000000000009, 11),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0), 0, 3.141592653589793), V(5, 0, 0), V(4.8708286933869696, 0, -1.1291713066130304), 1.5707963267948966, 1.3429969900969847, null, V(-3.061616997868383e-16, 0, -5), V(-1.1291713066130304, 0, -4.87082869338697)),
                new StraightEdge(new L3(V(6, 0, 0), V(-0.7071067811865475, 0, -0.7071067811865475)), V(4.8708286933869696, 0, -1.1291713066130304), V(1.1291713066130296, 0, -4.87082869338697), 1.5968893760546963, 6.8883919981838755)], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 0)), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(5, -6.123233995736766e-16, 0), V(-6.123233995736766e-16, -5, 0), 0, 3.141592653589793), V(4.950836318555471, -0.6994424542963525, 0), V(5, 0, -3.061616997868383e-16), 0.1403487973449601, 1.2246467991473532e-16, null, V(0.6994424542963528, 4.950836318555472, 0), V(1.2246467991473533e-15, 5, 0)),
                new StraightEdge(new L3(V(-6, 0, 0), V(1, 0, 0)), V(5, 0, 0), V(-4.999999999999999, 0, 0), 11, 1.0000000000000009),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(5, -6.123233995736766e-16, 0), V(-6.123233995736766e-16, -5, 0), 0, 3.141592653589793), V(-5, 6.123233995736766e-16, -3.061616997868383e-16), V(-4.950836318555471, -0.6994424542963525, 0), 3.141592653589793, 3.001243856244833, null, V(0, -5, 0), V(0.6994424542963525, -4.950836318555472, 0)),
                new StraightEdge(new L3(V(-6, 0, 0), V(0.8320502943378436, -0.554700196225229, 0)), V(-4.950836318555471, -0.6994424542963525, 0), V(0, -4, 0), 1.2609378166009386, 7.211102550927979),
                new StraightEdge(new L3(V(6, 0, 0), V(-0.8320502943378436, -0.554700196225229, 0)), V(0, -4, 0), V(4.950836318555471, -0.6994424542963525, 0), 7.211102550927979, 1.2609378166009386)], []),
            new PlaneFace(new PlaneSurface(new P3(V(0.485071250072666, -0.727606875108999, -0.485071250072666), 2.910427500435996)), [
                new PCurveEdge(new SemiEllipseCurve(V(1.411764705882352, -2.1176470588235303, -1.4117647058823533), V(2.0917534572581817, 2.7890046096775745, -2.091753457258183), V(-2.874840149008801, 3.5206637865439285e-16, -2.874840149008799), 0.7085839061491113, 3.141592653589793), V(0, -0.6994424542963525, -4.950836318555471), V(1.1291713066130291, 0, -4.87082869338697), 1.0373562345961393, 0.7085839061491117, null, V(3.262983117260863, 2.401508361946856, -0.3392794256594245), V(3.5440484447865424, 1.8149704259460573, 0.8215928058674522)),
                new StraightEdge(new L3(V(6, 0, 0), V(-0.7071067811865475, 0, -0.7071067811865475)), V(1.1291713066130296, 0, -4.87082869338697), V(4.8708286933869696, 0, -1.1291713066130304), 6.8883919981838755, 1.5968893760546963),
                new PCurveEdge(new SemiEllipseCurve(V(1.411764705882352, -2.1176470588235303, -1.4117647058823533), V(-2.0917534572581817, -2.7890046096775745, 2.091753457258183), V(2.874840149008801, -3.5206637865439285e-16, 2.874840149008799), 0, 2.433008747440682), V(4.87082869338697, 0, -1.1291713066130293), V(4.950836318555471, -0.6994424542963525, 0), 2.4330087474406805, 2.1042364189936533, null, V(0.8215928058674515, -1.81497042594606, 3.544048444786541), V(-0.3392794256594245, -2.4015083619468567, 3.262983117260862)),
                new StraightEdge(new L3(V(6, 0, 0), V(-0.8320502943378436, -0.554700196225229, 0)), V(4.950836318555471, -0.6994424542963525, 0), V(0, -4, 0), 1.2609378166009386, 7.211102550927979),
                new StraightEdge(new L3(V(0, 0, -6), V(0, -0.554700196225229, 0.8320502943378436)), V(0, -4, 0), V(0, -0.6994424542963525, -4.950836318555471), 7.211102550927979, 1.2609378166009386)], []),
            new PlaneFace(new PlaneSurface(new P3(V(-0.485071250072666, -0.727606875108999, -0.485071250072666), 2.910427500435996)), [
                new PCurveEdge(new SemiEllipseCurve(V(-1.4117647058823535, -2.117647058823529, -1.4117647058823533), V(-2.0917534572581813, 2.789004609677577, -2.0917534572581813), V(-2.8748401490088, 3.520663786543927e-16, 2.8748401490088007), 0.708583906149113, 3.141592653589793), V(-4.950836318555471, -0.6994424542963525, 0), V(-4.87082869338697, 5.965044768551438e-16, -1.1291713066130304), 1.03735623459614, 0.7085839061491122, null, V(-0.33927942565942515, 2.401508361946859, -3.2629831172608608), V(0.8215928058674524, 1.8149704259460604, -3.544048444786542)),
                new StraightEdge(new L3(V(-6, 0, 0), V(0.7071067811865475, 0, -0.7071067811865475)), V(-4.8708286933869696, 0, -1.1291713066130304), V(-1.1291713066130296, 0, -4.87082869338697), 1.5968893760546963, 6.8883919981838755),
                new PCurveEdge(new SemiEllipseCurve(V(-1.4117647058823535, -2.117647058823529, -1.4117647058823533), V(2.0917534572581813, -2.789004609677577, 2.0917534572581813), V(2.8748401490088, -3.520663786543927e-16, -2.8748401490088007), 0, 2.43300874744068), V(-1.1291713066130296, 1.382836026332681e-16, -4.87082869338697), V(0, -0.6994424542963525, -4.950836318555471), 2.4330087474406805, 2.1042364189936533, null, V(3.5440484447865406, -1.8149704259460617, -0.8215928058674513), V(3.2629831172608608, -2.401508361946859, 0.33927942565942426)),
                new StraightEdge(new L3(V(0, 0, -6), V(0, -0.554700196225229, 0.8320502943378436)), V(0, -0.6994424542963525, -4.950836318555471), V(0, -4, 0), 1.2609378166009386, 7.211102550927979),
                new StraightEdge(new L3(V(-6, 0, 0), V(0.8320502943378436, -0.554700196225229, 0)), V(0, -4, 0), V(-4.950836318555471, -0.6994424542963525, 0), 7.211102550927979, 1.2609378166009386)], [])], undefined)
        b2Equal(assert, a, b, a.and(b), result)
    },

    'B2T.sphere() - boxx'(assert) {
        const a = B2T.sphere()
        const b = B2T.box(0.2, 0.2,0.2, "").translate(0, 0.95, 0).rotateAB(V3.Y, V(0, -0.9341723589627158, 0.35682208977308993)).flipped()
        const result = new B2([
            new RotationFace(new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(0, 0, 1), V(0, 0, -1), 3.141592653589793, 0, null, V(-1, 1.2246467991473532e-16, 0), V(1, -1.2246467991473532e-16, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0, 0, -1), V(0, 0, 1), 0, 3.141592653589793, null, V(1, 0, 0), V(-1, 0, 0))], [[
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(0, -1, 0), 0, 3.141592653589793), V(0, -0.9866626624629129, 0.16277834776651348), V(0, -0.9341723589627158, 0.3568220897730897), 1.734302234118049, 1.9356601549083796, null, V(0, 0.16277834776651368, 0.9866626624629129), V(0, 0.35682208977308977, 0.9341723589627158)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(-1, 1.2246467991473532e-16, -4.6777345306052285e-17), V(-1.310943090317052e-16, -0.9341723589627158, 0.3568220897730897), 0, 3.141592653589793), V(0, -0.9341723589627158, 0.3568220897730897), V(0.2, -0.9152982445082949, 0.34961281955905665), 1.5707963267948968, 1.7721542475852277, null, V(1, 2.776169300024433e-17, -1.0604023140853287e-17), V(0.9797958971132712, 0.18683447179254326, -0.07136441795461798)),
                new PCurveEdge(new SemiEllipseCurve(V(0.2, -2.2884754904439333e-18, 0), V(0, 0, 0.9797958971132712), V(-1.1211194480906224e-17, -0.9797958971132712, 0), 2.2662332591841976e-17, 3.141592653589793), V(0.2, -0.9152982445082949, 0.34961281955905665), V(0.2, -0.9673910674187775, 0.1554172534770778), 1.2059324986814137, 1.4115014298425876, null, V(-4.000401843529488e-18, -0.3496128195590566, -0.9152982445082949), V(-1.7783428768720205e-18, -0.1554172534770777, -0.9673910674187775)),
                new PCurveEdge(new SemiEllipseCurve(V(-4.129479174250033e-19, -0.07136441795461794, -0.18683447179254323), V(1.276729399545399e-16, 0.9152982445082949, -0.3496128195590565), V(0.9797958971132712, -1.1999039092174769e-16, 4.366667272259081e-17), 1.4927486281828333, 3.141592653589793), V(0.2, -0.9673910674187775, 0.1554172534770778), V(0, -0.9866626624629129, 0.16277834776651348), 2.936023722428619, 3.141592653589793, null, V(-0.9591663046625438, -0.1868344717925432, 0.07136441795461794), V(-0.9797958971132712, 7.898684381520212e-18, -8.5145068120284055e-19))]]),
            new RotationFace(new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0, 0, 1), V(0, 0, -1), 3.141592653589793, 0, null, V(1, 0, 0), V(-1, 0, 0)),
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(0, 0, -1), V(0, 0, 1), 0, 3.141592653589793, null, V(-1, 1.2246467991473532e-16, 0), V(1, -1.2246467991473532e-16, 0))], []),
            new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 0)), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(0, -1, 0), 0, 3.141592653589793), V(0, -0.9341723589627158, 0.3568220897730897), V(0, -0.9866626624629129, 0.16277834776651348), 1.9356601549083796, 1.734302234118049, null, V(0, -0.35682208977308977, -0.9341723589627158), V(0, -0.16277834776651368, -0.9866626624629129)),
                new StraightEdge(new L3(V(0, -0.9588281589691978, 0.15214651349189204), V(0, -0.9341723589627158, 0.3568220897730897)), V(0, -0.9866626624629129, 0.16277834776651348), V(0, -0.9588281589691978, 0.15214651349189204), 0.029795897113271352, 0),
                new StraightEdge(new L3(V(0, -0.8874637410145799, 0.3389809852844352), V(0, -0.3568220897730897, -0.9341723589627158)), V(0, -0.9588281589691978, 0.15214651349189204), V(0, -0.8874637410145799, 0.3389809852844352), 0.2, 0),
                new StraightEdge(new L3(V(0, -0.8874637410145799, 0.3389809852844352), V(0, -0.9341723589627158, 0.3568220897730897)), V(0, -0.8874637410145799, 0.3389809852844352), V(0, -0.9341723589627158, 0.3568220897730897), 0, 0.050000000000000044)], []),
            new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), -0.2)), [
                new PCurveEdge(new SemiEllipseCurve(V(0.2, -2.2884754904439333e-18, 0), V(0, 0, 0.9797958971132712), V(-1.1211194480906224e-17, -0.9797958971132712, 0), 2.2662332591841976e-17, 3.141592653589793), V(0.2, -0.9673910674187775, 0.1554172534770778), V(0.2, -0.9152982445082949, 0.34961281955905665), 1.4115014298425876, 1.2059324986814137, null, V(1.7783428768720205e-18, 0.1554172534770777, 0.9673910674187775), V(4.000401843529488e-18, 0.3496128195590566, 0.9152982445082949)),
                new StraightEdge(new L3(V(0.2, -1.074298212807123, 0.41034540323905316), V(0, 0.9341723589627158, -0.3568220897730897)), V(0.2, -0.9152982445082949, 0.34961281955905665), V(0.2, -0.8874637410145799, 0.3389809852844352), 0.17020410288672871, 0.2),
                new StraightEdge(new L3(V(0.2, -0.8874637410145799, 0.3389809852844352), V(0, -0.3568220897730897, -0.9341723589627158)), V(0.2, -0.8874637410145799, 0.3389809852844352), V(0.2, -0.9588281589691978, 0.15214651349189204), 0, 0.2),
                new StraightEdge(new L3(V(0.2, -1.145662630761741, 0.22351093144651), V(0, 0.9341723589627158, -0.3568220897730897)), V(0.2, -0.9588281589691978, 0.15214651349189204), V(0.2, -0.9673910674187775, 0.1554172534770778), 0.2, 0.1908336953374561)], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, -0.3568220897730897, -0.9341723589627158), 0)), [
                new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(-1, 1.2246467991473532e-16, -4.6777345306052285e-17), V(-1.310943090317052e-16, -0.9341723589627158, 0.3568220897730897), 0, 3.141592653589793), V(0.2, -0.9152982445082949, 0.34961281955905665), V(0, -0.9341723589627158, 0.3568220897730897), 1.7721542475852277, 1.5707963267948968, null, V(-0.9797958971132712, -0.18683447179254326, 0.07136441795461798), V(-1, -2.776169300024433e-17, 1.0604023140853287e-17)),
                new StraightEdge(new L3(V(0, -0.8874637410145799, 0.3389809852844352), V(0, -0.9341723589627158, 0.3568220897730897)), V(0, -0.9341723589627158, 0.3568220897730897), V(0, -0.8874637410145799, 0.3389809852844352), 0.050000000000000044, 0),
                new StraightEdge(new L3(V(0.2, -0.8874637410145799, 0.3389809852844352), V(-1, 0, 0)), V(0, -0.8874637410145799, 0.3389809852844352), V(0.2, -0.8874637410145799, 0.3389809852844352), 0.2, 0),
                new StraightEdge(new L3(V(0.2, -1.074298212807123, 0.41034540323905316), V(0, 0.9341723589627158, -0.3568220897730897)), V(0.2, -0.8874637410145799, 0.3389809852844352), V(0.2, -0.9152982445082949, 0.34961281955905665), 0.2, 0.17020410288672871)], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0.3568220897730897, 0.9341723589627159), -0.20000000000000004)), [
                new PCurveEdge(new SemiEllipseCurve(V(-4.129479174250033e-19, -0.07136441795461794, -0.18683447179254323), V(1.276729399545399e-16, 0.9152982445082949, -0.3496128195590565), V(0.9797958971132712, -1.1999039092174769e-16, 4.366667272259081e-17), 1.4927486281828333, 3.141592653589793), V(0, -0.9866626624629129, 0.16277834776651348), V(0.2, -0.9673910674187775, 0.1554172534770778), 3.141592653589793, 2.936023722428619, null, V(0.9797958971132712, -7.898684381520212e-18, 8.5145068120284055e-19), V(0.9591663046625438, 0.1868344717925432, -0.07136441795461794)),
                new StraightEdge(new L3(V(0.2, -1.145662630761741, 0.22351093144651), V(0, 0.9341723589627158, -0.3568220897730897)), V(0.2, -0.9673910674187775, 0.1554172534770778), V(0.2, -0.9588281589691978, 0.15214651349189204), 0.1908336953374561, 0.2),
                new StraightEdge(new L3(V(0.2, -0.9588281589691978, 0.15214651349189204), V(-1, 0, 0)), V(0.2, -0.9588281589691978, 0.15214651349189204), V(0, -0.9588281589691978, 0.15214651349189204), 0, 0.2),
                new StraightEdge(new L3(V(0, -0.9588281589691978, 0.15214651349189204), V(0, -0.9341723589627158, 0.3568220897730897)), V(0, -0.9588281589691978, 0.15214651349189204), V(0, -0.9866626624629129, 0.16277834776651348), 0, 0.029795897113271352)], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, -0.9341723589627158, 0.3568220897730897), 0.95)), [
                new StraightEdge(new L3(V(0, -0.8874637410145799, 0.3389809852844352), V(0, -0.3568220897730897, -0.9341723589627158)), V(0, -0.8874637410145799, 0.3389809852844352), V(0, -0.9588281589691978, 0.15214651349189204), 0, 0.2),
                new StraightEdge(new L3(V(0.2, -0.9588281589691978, 0.15214651349189204), V(-1, 0, 0)), V(0, -0.9588281589691978, 0.15214651349189204), V(0.2, -0.9588281589691978, 0.15214651349189204), 0.2, 0),
                new StraightEdge(new L3(V(0.2, -0.8874637410145799, 0.3389809852844352), V(0, -0.3568220897730897, -0.9341723589627158)), V(0.2, -0.9588281589691978, 0.15214651349189204), V(0.2, -0.8874637410145799, 0.3389809852844352), 0.2, 0),
                new StraightEdge(new L3(V(0.2, -0.8874637410145799, 0.3389809852844352), V(-1, 0, 0)), V(0.2, -0.8874637410145799, 0.3389809852844352), V(0, -0.8874637410145799, 0.3389809852844352), 0, 0.2)], [])], undefined)
        b2Equal(assert, a, b, a.and(b), result)
    },

    'splitsVolumeEnclosingCone 2'(assert) {
        const a = B2T.box(1,1,1,'box').minus(B2T.box(1 / 3, 1 / 3, 1,'cut'))
        assert.equal(splitsVolumeEnclosingCone(a, V(1/3,1/3,0), V3.X), ALONG_EDGE_OR_PLANE)
    },

    'Face.containsPoint2 2'(assert) {
        const face = new PlaneFace(new P3(V(0, 0, -1), 0), [
            new StraightEdge(new L3(V(0, 0.3333333333333333, 0), V(-1, 0, 0)), V(0.3333333333333333, 0.3333333333333333, 0), V(0, 0.3333333333333333, 0), -0.3333333333333333, 0),
            new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 0.3333333333333333, 0), V(0, 1, 0), 0.3333333333333333, 1),
            new StraightEdge(new L3(V(0, 1, 0), V(1, 0, 0)), V(0, 1, 0), V(1, 1, 0), 0, 1),
            new StraightEdge(new L3(V(1, 1, 0), V(0, -1, 0)), V(1, 1, 0), V(1, 0, 0), 0, 1),
            new StraightEdge(new L3(V(1, 0, 0), V(-1, 0, 0)), V(1, 0, 0), V(0.33333333333333326, 0, 0), 0, 0.6666666666666667),
            new StraightEdge(new L3(V(0.3333333333333333, 0, 0), V(0, 1, 0)), V(0.33333333333333326, 0, 0), V(0.3333333333333333, 0.3333333333333333, 0), 0, 0.3333333333333333)])
        assert.equal(face.containsPoint2(V(0.33333345254262287, 0.3333333333333333, 0)), PointVsFace.INSIDE)
    },

    'fuz test'(assert) {
        const a = B2T.box(1,1,1,'box').minus(B2T.box(1 / 3, 1 / 3, 1,'cut'))
        const b = B2T.box(3, 1 / 3, 1 / 3).translate(-1).flipped()
        const result = new B2([
            new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), 0)), [
                new StraightEdge(new L3(V(0, 0.3333333333333333, 0), V(0, 0, 1)), V(0, 0.3333333333333333, 0), V(0, 0.3333333333333333, 0.3333333333333333), 0, 0.3333333333333333),
                new StraightEdge(new L3(V(0, 0.3333333333333333, 0), V(0, 0, 1)), V(0, 0.3333333333333333, 0.3333333333333333), V(0, 0.3333333333333333, 1), 0.3333333333333333, 1),
                new StraightEdge(new L3(V(0, 0, 1), V(0, 1, 0)), V(0, 0.3333333333333333, 1), V(0, 1, 1), 0.3333333333333333, 1),
                new StraightEdge(new L3(V(0, 1, 0), V(0, 0, 1)), V(0, 1, 1), V(0, 1, 0), 1, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 1, 0), V(0, 0.3333333333333333, 0), 1, 0.3333333333333333)], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0)), [
                new StraightEdge(new L3(V(0, 0, 0.3333333333333333), V(1, 0, 0)), V(0.3333333333333333, 0, 0.3333333333333333), V(1, 0, 0.3333333333333333), 0.3333333333333333, 1),
                new StraightEdge(new L3(V(1, 0, 0), V(0, 0, 1)), V(1, 0, 0.3333333333333333), V(1, 0, 1), 0.3333333333333333, 1),
                new StraightEdge(new L3(V(1, 0, 1), V(-1, 0, 0)), V(1, 0, 1), V(0.33333333333333326, 0, 1), 0, 0.6666666666666667),
                new StraightEdge(new L3(V(0.3333333333333333, 0, 0), V(0, 0, -1)), V(0.33333333333333326, 0, 1), V(0.3333333333333333, 0, 0.3333333333333333), -1, -0.3333333333333333)], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 0)), [
                new StraightEdge(new L3(V(0, 0.3333333333333333, 0), V(-1, 0, 0)), V(1, 0.33333333333333326, 0), V(0.3333333333333333, 0.3333333333333333, 0), -1, -0.3333333333333333),
                new StraightEdge(new L3(V(0, 0.3333333333333333, 0), V(-1, 0, 0)), V(0.3333333333333333, 0.3333333333333333, 0), V(0, 0.3333333333333333, 0), -0.3333333333333333, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 0.3333333333333333, 0), V(0, 1, 0), 0.3333333333333333, 1),
                new StraightEdge(new L3(V(0, 1, 0), V(1, 0, 0)), V(0, 1, 0), V(1, 1, 0), 0, 1),
                new StraightEdge(new L3(V(1, 1, 0), V(0, -1, 0)), V(1, 1, 0), V(1, 0.33333333333333326, 0), 0, 0.6666666666666667)], []),
            new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 1)), [
                new StraightEdge(new L3(V(1, 0.3333333333333333, 0), V(0, 0, -1)), V(1, 0.3333333333333333, 0.3333333333333333), V(1, 0.33333333333333326, 0), -0.3333333333333333, 0),
                new StraightEdge(new L3(V(1, 1, 0), V(0, -1, 0)), V(1, 0.33333333333333326, 0), V(1, 1, 0), 0.6666666666666667, 0),
                new StraightEdge(new L3(V(1, 1, 0), V(0, 0, 1)), V(1, 1, 0), V(1, 1, 1), 0, 1),
                new StraightEdge(new L3(V(1, 1, 1), V(0, -1, 0)), V(1, 1, 1), V(1, 0, 1), 0, 1),
                new StraightEdge(new L3(V(1, 0, 0), V(0, 0, 1)), V(1, 0, 1), V(1, 0, 0.3333333333333333), 1, 0.3333333333333333),
                new StraightEdge(new L3(V(1, 0, 0.3333333333333333), V(0, 1, 0)), V(1, 0, 0.3333333333333333), V(1, 0.3333333333333333, 0.3333333333333333), 0, 0.3333333333333333)], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), -0.3333333333333333)), [
                new StraightEdge(new L3(V(0, 0.3333333333333333, 0), V(0, 0, 1)), V(0, 0.3333333333333333, 0.3333333333333333), V(0, 0.3333333333333333, 0), 0.3333333333333333, 0),
                new StraightEdge(new L3(V(0, 0.3333333333333333, 0), V(-1, 0, 0)), V(0, 0.3333333333333333, 0), V(0.3333333333333333, 0.3333333333333333, 0), 0, -0.3333333333333333),
                new StraightEdge(new L3(V(0.3333333333333333, 0.3333333333333333, 0), V(0, 0, 1)), V(0.3333333333333333, 0.3333333333333333, 0), V(0.33333333333333326, 0.3333333333333333, 0.3333333333333333), 0, 0.3333333333333333),
                new StraightEdge(new L3(V(0, 0.3333333333333333, 0.3333333333333333), V(1, 0, 0)), V(0.3333333333333333, 0.3333333333333333, 0.3333333333333333), V(0, 0.3333333333333333, 0.3333333333333333), 0.3333333333333333, 0)], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), -0.3333333333333333)), [
                new StraightEdge(new L3(V(0, 0.3333333333333333, 0.3333333333333333), V(1, 0, 0)), V(0, 0.3333333333333333, 0.3333333333333333), V(0.3333333333333333, 0.3333333333333333, 0.3333333333333333), 0, 0.3333333333333333),
                new StraightEdge(new L3(V(0.3333333333333333, 0.3333333333333333, 0), V(0, 0, 1)), V(0.3333333333333333, 0.3333333333333333, 0.3333333333333333), V(0.3333333333333333, 0.3333333333333333, 1), 0.3333333333333333, 1),
                new StraightEdge(new L3(V(0, 0.3333333333333333, 1), V(1, 0, 0)), V(0.3333333333333333, 0.3333333333333333, 1), V(0, 0.3333333333333333, 1), 0.3333333333333333, 0),
                new StraightEdge(new L3(V(0, 0.3333333333333333, 0), V(0, 0, 1)), V(0, 0.3333333333333333, 1), V(0, 0.3333333333333333, 0.3333333333333333), 1, 0.3333333333333333)], []),
            new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), -0.3333333333333333)), [
                new StraightEdge(new L3(V(0.3333333333333333, 0, 0.3333333333333333), V(0, -1, 0)), V(0.3333333333333333, 0.3333333333333333, 0.3333333333333333), V(0.3333333333333333, 0, 0.3333333333333333), -0.3333333333333333, 0),
                new StraightEdge(new L3(V(0.3333333333333333, 0, 0), V(0, 0, -1)), V(0.3333333333333333, 0, 0.3333333333333333), V(0.33333333333333326, 0, 1), -0.3333333333333333, -1),
                new StraightEdge(new L3(V(0.3333333333333333, 0, 1), V(0, -1, 0)), V(0.33333333333333326, 0, 1), V(0.3333333333333333, 0.3333333333333333, 1), 0, -0.3333333333333333),
                new StraightEdge(new L3(V(0.3333333333333333, 0.3333333333333333, 0), V(0, 0, 1)), V(0.3333333333333333, 0.3333333333333333, 1), V(0.3333333333333333, 0.3333333333333333, 0.3333333333333333), 1, 0.3333333333333333)], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 1)), [
                new StraightEdge(new L3(V(0, 0.3333333333333333, 1), V(1, 0, 0)), V(0, 0.3333333333333333, 1), V(0.3333333333333333, 0.3333333333333333, 1), 0, 0.3333333333333333),
                new StraightEdge(new L3(V(0.3333333333333333, 0, 1), V(0, -1, 0)), V(0.3333333333333333, 0.3333333333333333, 1), V(0.33333333333333326, 0, 1), -0.3333333333333333, 0),
                new StraightEdge(new L3(V(1, 0, 1), V(-1, 0, 0)), V(0.33333333333333326, 0, 1), V(1, 0, 1), 0.6666666666666667, 0),
                new StraightEdge(new L3(V(1, 1, 1), V(0, -1, 0)), V(1, 0, 1), V(1, 1, 1), 1, 0),
                new StraightEdge(new L3(V(0, 1, 1), V(1, 0, 0)), V(1, 1, 1), V(0, 1, 1), 1, 0),
                new StraightEdge(new L3(V(0, 0, 1), V(0, 1, 0)), V(0, 1, 1), V(0, 0.3333333333333333, 1), 1, 0.3333333333333333)], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 1, 0), 1)), [
                new StraightEdge(new L3(V(0, 1, 0), V(1, 0, 0)), V(1, 1, 0), V(0, 1, 0), 1, 0),
                new StraightEdge(new L3(V(0, 1, 0), V(0, 0, 1)), V(0, 1, 0), V(0, 1, 1), 0, 1),
                new StraightEdge(new L3(V(0, 1, 1), V(1, 0, 0)), V(0, 1, 1), V(1, 1, 1), 0, 1),
                new StraightEdge(new L3(V(1, 1, 0), V(0, 0, 1)), V(1, 1, 1), V(1, 1, 0), 1, 0)], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), -0.3333333333333333)), [
                new StraightEdge(new L3(V(0, 0.3333333333333333, 0), V(-1, 0, 0)), V(0.3333333333333333, 0.3333333333333333, 0), V(1, 0.33333333333333326, 0), -0.3333333333333333, -1),
                new StraightEdge(new L3(V(1, 0.3333333333333333, 0), V(0, 0, -1)), V(1, 0.33333333333333326, 0), V(1, 0.3333333333333333, 0.3333333333333333), 0, -0.3333333333333333),
                new StraightEdge(new L3(V(-1, 0.3333333333333333, 0.3333333333333333), V(1, 0, 0)), V(1, 0.3333333333333333, 0.3333333333333333), V(0.33333333333333326, 0.3333333333333333, 0.3333333333333333), 2, 1.3333333333333333),
                new StraightEdge(new L3(V(0.3333333333333333, 0.3333333333333333, 0), V(0, 0, 1)), V(0.33333333333333326, 0.3333333333333333, 0.3333333333333333), V(0.3333333333333333, 0.3333333333333333, 0), 0.3333333333333333, 0)], []),
            new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), -0.3333333333333333)), [
                new StraightEdge(new L3(V(0, 0, 0.3333333333333333), V(1, 0, 0)), V(1, 0, 0.3333333333333333), V(0.3333333333333333, 0, 0.3333333333333333), 1, 0.3333333333333333),
                new StraightEdge(new L3(V(0.3333333333333333, 0, 0.3333333333333333), V(0, -1, 0)), V(0.3333333333333333, 0, 0.3333333333333333), V(0.3333333333333333, 0.3333333333333333, 0.3333333333333333), 0, -0.3333333333333333),
                new StraightEdge(new L3(V(-1, 0.3333333333333333, 0.3333333333333333), V(1, 0, 0)), V(0.33333333333333326, 0.3333333333333333, 0.3333333333333333), V(1, 0.3333333333333333, 0.3333333333333333), 1.3333333333333333, 2),
                new StraightEdge(new L3(V(1, 0, 0.3333333333333333), V(0, 1, 0)), V(1, 0.3333333333333333, 0.3333333333333333), V(1, 0, 0.3333333333333333), 0.3333333333333333, 0)], [])], false)
        b2Equal(assert, a, b, a.and(b), result)
    },

    'B2.withMergedFaces'(assert) {
        const box = B2T.box(5, 5, 5)
        const boxToMerge = new B2(box.faces.filter((face:PlaneFace) => face.surface.plane.normal.x != 1).concat(
            box.translate(5, 0, 0).faces.filter((face:PlaneFace) => face.surface.plane.normal.x != -1)
        ), false)

        assert.equal(boxToMerge.faces.length, 10)
        const boxMerged = boxToMerge.withMergedFaces()
        b2Equal(assert, boxToMerge, B2.EMPTY, boxMerged, B2T.box(10, 5, 5))
    }
})