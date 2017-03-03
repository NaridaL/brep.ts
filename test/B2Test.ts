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
        let box = B2T.box(10, 10, 5, 'a')
        let box2 = B2T.box(10, 10, 4, 'b').translate(3, 3, 0)
        b2Equal(assert, box, box2, box.intersection(box2, true, true), B2T.box(7, 7, 4).translate(3, 3, 0))
    },

    'B2T.box(10, 10, 5) + B2T.box(10, 10, 5).translate(3, 3, 0) overlapping faces result contains union of faces'(assert) {
        const boxA = B2T.box(10, 10, 5, 'boxA').flipped(), boxB = boxA.translate(3, 3, 0)
        //const result = B2T.extrudeVertices([V3.ZERO, V(0, 10), V(3, 10), V(3, 13), V(13, 13), V(13, 3), V(10, 3), V(10, 0)], P3.XY.flipped(), V(0, 0, 5), 'result').flipped()
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

    'B2T.box(10, 10, 10) - B2T.box(5, 5, 5)'(assert) {
        const box = B2T.box(10, 10, 10, 'box0'), box2 = B2T.box(1,1,1).scale(8 / 3 ** 0.5).transform(M4.rotationAB(V3.ONES, V3.X))
        const result = new B2([
            new PlaneFace(new P3(V(0, -1, 0), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(0.7071067811865476, 0, 0.7071067811865475)), V(0, 0, 0), V(3.381197846482994, 0, 3.3811978464829933), 0, 4.781735851562952),
                new StraightEdge(new L3(V(2.791322082992153, 0, 3.813016875511774), V(0.8068982213550735, 0, -0.5906904945688721)), V(3.381197846482994, 0, 3.3811978464829933), V(8, 0, 0), 0.7310411002024879, 6.45518577084059),
                new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(8, 0, 0), V(10, 0, 0), 1.9999999999999984, 0),
                new StraightEdge(new L3(V(10, 0, 0), V(0, 0, 1)), V(10, 0, 0), V(10, 0, 10), 0, 10),
                new StraightEdge(new L3(V(10, 0, 10), V(-1, 0, 0)), V(10, 0, 10), V(0, 0, 10), 0, 10),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 10), V(0, 0, 0), 10, 0)]),
            new PlaneFace(new P3(V(0, 0, -1), 0), [
                new StraightEdge(new L3(V(0, 0, 0), V(-0.7071067811865476, -0.7071067811865475, 0)), V(3.381197846482994, 3.3811978464829933, 0), V(0, 0, 0), -4.781735851562952, 0),
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 0, 0), V(0, 10, 0), 0, 10),
                new StraightEdge(new L3(V(0, 10, 0), V(1, 0, 0)), V(0, 10, 0), V(10, 10, 0), 0, 10),
                new StraightEdge(new L3(V(10, 10, 0), V(0, -1, 0)), V(10, 10, 0), V(10, 0, 0), 0, 10),
                new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(10, 0, 0), V(8.000000000000004, 0, 0), 0, 1.9999999999999962),
                new StraightEdge(new L3(V(2.7913220829921492, 3.813016875511774, 0), V(-0.8068982213550737, 0.590690494568872, 0)), V(8.000000000000004, 0, 0), V(3.381197846482994, 3.3811978464829933, 0), -6.455185770840593, -0.7310411002024895)]),
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
                new StraightEdge(new L3(V(2.7913220829921492, 3.813016875511774, 0), V(-0.8068982213550737, 0.590690494568872, 0)), V(3.381197846482994, 3.3811978464829933, 0), V(8.000000000000004, 0, 0), -0.7310411002024895, -6.455185770840593),
                new StraightEdge(new L3(V(5.333333333333336, 2.666666666666666, 2.666666666666666), V(0.5773502691896257, -0.5773502691896258, -0.5773502691896258)), V(8, -6.661338147750939e-16, -8.881784197001252e-16), V(5.333333333333336, 2.666666666666666, 2.666666666666666), 4.618802153517006, 0),
                new StraightEdge(new L3(V(2.666666666666668, 3.642734410091836, -0.9760677434251697), V(0.5773502691896258, -0.21132486540518713, 0.7886751345948129)), V(5.333333333333336, 2.666666666666666, 2.666666666666666), V(3.381197846482994, 3.3811978464829933, 0), 4.618802153517006, 1.2376043070340121)]),
            new PlaneFace(new P3(V(-0.5773502691896261, 0.2113248654051888, -0.7886751345948121), -4.618802153517), [
                new StraightEdge(new L3(V(2.791322082992153, 0, 3.813016875511774), V(0.8068982213550735, 0, -0.5906904945688721)), V(8, 0, 0), V(3.381197846482994, 0, 3.3811978464829933), 6.45518577084059, 0.7310411002024879),
                new StraightEdge(new L3(V(2.666666666666668, -0.9760677434251697, 3.642734410091836), V(0.5773502691896258, 0.7886751345948129, -0.21132486540518713)), V(3.381197846482994, 0, 3.3811978464829933), V(5.333333333333336, 2.666666666666666, 2.666666666666666), 1.2376043070340121, 4.618802153517006),
                new StraightEdge(new L3(V(5.333333333333336, 2.666666666666666, 2.666666666666666), V(0.5773502691896257, -0.5773502691896258, -0.5773502691896258)), V(5.333333333333336, 2.666666666666666, 2.666666666666666), V(8, -6.661338147750939e-16, -8.881784197001252e-16), 0, 4.618802153517006)])], false)
        b2Equal(assert, box, box2, box.minus(box2), result)
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