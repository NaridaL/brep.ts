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
	'Face.trasnform'(assert) {
    	const loop = Edge.pathFromSVG('m 2 0, 1 1, -2 1 Z')
		const face = new PlaneFace(P3.XY, loop)
		const faceMirrored = face.mirroredX()
	},
	'rotateEdges'(assert){
		const chain = Edge.pathFromSVG('m 2 0, 1 1, -2 1 Z')
			.map(e => e.transform(M4.rotationX(90 * DEG)))
		linkB3(assert, {edges: chain})
		const a = B2T.rotateEdges(chain,  TAU / 3, 'rot').flipped()
		a.assertSanity()
		linkB3(assert, {a})
		b2equals(assert, a, B2.EMPTY)
		const b = B2T.box(4,4,1.4).translate(0, 0, -0.2).rotateX(5*DEG).rotateY(-10*DEG)
		b2EqualAnd(assert, a, b, B2.EMPTY)
	},
	'rotateEdges 2'(assert){
		const chain = [
			new StraightEdge(new L3(V(198.46477746372744, 0, 244.94352900661897), V(-0.16064282877602398, 0, -0.9870126045612776)), V(173.9128996557253, 0, 94.093266677553), V(198.46477746372744, 0, 244.94352900661897), 152.83519342300417, 0),
			new StraightEdge(new L3(V(131.35224103228387, 0, 180.2100595549249), V(0.7197488536413841, 0, 0.6942345336281635)), V(198.46477746372744, 0, 244.94352900661897), V(131.35224103228387, 0, 180.2100595549249), 93.24438113642698, 0),
			new StraightEdge(new L3(V(173.9128996557253, 0, 94.093266677553), V(-0.44306356566594673, 0, 0.8964901989310186)), V(131.35224103228387, 0, 180.2100595549249), V(173.9128996557253, 0, 94.093266677553), 96.05993794472955, 0)
		]
		Edge.assertLoop(chain)
		linkB3(assert, {edges: chain})
		const a = B2T.rotateEdges(chain,  TAU / 3, 'rot')
		a.assertSanity()
		//const a = B2T.rotateEdges(Edge.reverseLoop(chain),  TAU / 3, 'rot')
		linkB3(assert, {a})
		b2equals(assert, a, B2.EMPTY)
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

    'tetrahedron'(assert) {
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
        b2equals(assert, a, result)
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
    	function tst(b2, p, curve, t0, dir, result) {
    		linkB3(assert, {a: b2, edges: [Edge.forCurveAndTs(curve, t0, t0 + dir)]})
		    assert.equal(splitsVolumeEnclosingCone2(b2, p, curve, t0, dir), result)
	    }
        const brep = B2T.box(5, 5, 5)

        tst(brep, V3.O, new L3(V3.O, V(1, 1, 1).unit()), 0, 1, INSIDE)
        tst(brep, V3.O, new L3(V3.O, V(-1, 1, 1).unit()), 0, 1, OUTSIDE)
        tst(brep, V3.O, new L3(V3.O, V3.X), 0, -1, OUTSIDE)
        tst(brep, V(0, 5, 5), new L3(V(0, 5, 5), V3.Y), 0, 1, OUTSIDE)
        tst(brep, V3.O, new L3(V3.O, V3.X), 0, 1, ALONG_EDGE_OR_PLANE)
	    tst(brep, V3.O, new L3(V3.O, V(1, 1, 0).unit()), 0, 1, ALONG_EDGE_OR_PLANE)

	    tst(brep, V(5,0,0), SemiEllipseCurve.UNIT.translate(4), 0, 1, ALONG_EDGE_OR_PLANE)
	    tst(brep, V(5,0,0), SemiEllipseCurve.UNIT.translate(4).rotateX(5 * DEG), 0, 1, INSIDE)
	    tst(brep, V(0,0,0), SemiEllipseCurve.UNIT.translate(-1), 0, 1, OUTSIDE)

	    const semicyl = B2T.cylinder(1,1,TAU /2)
	    tst(semicyl, V(-1,0,0), L3.X, -1, -1, OUTSIDE)
	    tst(semicyl, V(-1,0,0), L3.X, -1, 1, ALONG_EDGE_OR_PLANE)
	    tst(semicyl.flipped(), V(-1,0,0), L3.X, -1, -1, INSIDE)
	    tst(semicyl.flipped(), V(-1,0,0), L3.X, -1, 1, ALONG_EDGE_OR_PLANE)
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
        let line = L3.X.translate(0, 0, -1)
        let result = planeFaceEdgeISPsWithPlane(brep.faces[2], L3.Y.translate(0, 0, -1), P3.XY.translate(0, 0, -1)).map(is => is.p)
        assert.V3ArraysLike(result, [V(0, 0, -1), V(0, 10, -1)])
        result = planeFaceEdgeISPsWithPlane(brep.faces[2], L3.Y, P3.XY).map(is => is.p)
        assert.V3ArraysLike(result, [V(0, 0, 0), V(0, 10, 0)])
        result = planeFaceEdgeISPsWithPlane(brep.translate(0, 0, 10).faces[2], L3.Y.translate(0, 0, 6), P3.XY.translate(0, 0, 6)).map(is => is.p)
        assert.V3ArraysLike(result, [V(0, 0, 6), V(0, 10, 6)])
    },

    'box(5, 5, 5) - box(1, 1, 6)'(assert) {
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

    'box(5, 5, 5) - box(1, 1, 5)'(assert) {
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

    'box(10, 10, 5) + box(10, 10, 5).translate(0, 0, 5) touching coplanar faces result contains both'(assert) {
        const a = B2T.box(10, 10, 5), b = a.translate(0, 0, 5)

        // the expected result is the union of a and b's faces, minus the ones where they touch, which we remove with an AABB test
        const testAABB = AABB.forXYZ(12, 12, 2).translate(-1, -1, 4)
        const result = new B2(a.faces.concat(b.faces).filter(face => !testAABB.containsAABB(face.getAABB())), false)
        b2Equal(assert, a, b, b.plus(a), result)
    },

    'box(10, 10, 5) - box(10, 5, 5).translate(0, 0, 5) - box(5, 5, 5).translate(5, 5, 5)'(assert) {
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

    'box(10, 10, 5) && box(10, 10, 5).translate(3, 3, 0) overlapping faces result contains intersection of faces'(assert) {
        const box = B2T.box(10, 10, 5, 'a')
        const box2 = B2T.box(10, 10, 4, 'b').translate(3, 3, 0)
        b2Equal(assert, box, box2, box.intersection(box2, true, true), B2T.box(7, 7, 4).translate(3, 3, 0))
    },

    'box(10, 10, 5) + box(10, 10, 5).translate(3, 3, 0) overlapping faces result contains union of faces'(assert) {
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


    'box(10,10,5) + box(4,10,2).translate(2, 0, 3) overlapping faces result contains union of faces 2'(assert) {
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

    'box(10,10,10) + box(10,12,12).translate(0, 0, -2) overlapping faces result contains union of faces 3'(assert) {
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

    'box(10, 10, 10) - box().scale(8 / 3 ** 0.5).rotateAB(V3.ONES, V3.X)'(assert) {
        const box = B2T.box(10, 10, 10, 'box0'), box2 = B2T.box(1,1,1).scale(8 / 3 ** 0.5).rotateAB(V3.XYZ, V3.X)
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

    'box(10, 10, 10) - tetrahedron(V(5,0,10),V(5,3,10),V(7,1,9),V(7,1,11))'(assert) {
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

    'box(10, 10, 10) - box().rotateAB(V3.ONES, V3.Y).translate(5,0,10)'(assert) {
        const box = B2T.box(10, 10, 10, 'box0'), box2 = B2T.box().rotateAB(V3.XYZ, V3.Y).translate(5,0,10).flipped()
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

    'menger(1) - box(2,2,1).rotateZ(-45*DEG)'(assert) {
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

    'sphere(5) - box(12,2,3).translate(-6, 0,1)'(assert) {
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
                new StraightEdge(new L3(V(-6, 2, 4), V(1, 0, 0)), V(-2.236067977499787, 2, 4), V(2.2360679774997863, 2, 4), 3.763932022500213, 8.236067977499786)])], false)
        b2EqualAnd(assert, a, b, result)
    },

    'sphere(5) AND octahedron().scale(6)'(assert) {
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
                new StraightEdge(new L3(V(0, -6, 0), V(0, 0.7071067811865475, -0.7071067811865475)), V(0, -4.8708286933869696, -1.1291713066130304), V(0, -1.1291713066130296, -4.87082869338697), 1.5968893760546963, 6.8883919981838755)], [])], false)
        b2EqualAnd(assert, a, b, result)
    },

    'sphere(5) - box(12,2,3).translate(-6, 0,1) 2'(assert) {
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
                new StraightEdge(new L3(V(-6, 1.505752286012999, 2.7807750813696934), V(1, 0, 0)), V(-3.8729833462074166, 1.505752286012999, 2.7807750813696934), V(3.872983346207416, 1.505752286012999, 2.7807750813696934), 2.1270166537925834, 9.872983346207416)])], false)
        b2EqualAnd(assert, a, b, result)
    },

    'sphere(5) AND tetrahedron(V3.ZERO, V3.Y, V3.Z, V(7,0,0))'(assert) {
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
                new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 1, 0), V(0, 0, 0), 1, 0)])], false)
        b2EqualAnd(assert, a, b, result)
    },

    'sphere(5) AND tetrahedron(V3.ZERO, V3.X, V3.Y, V3.Z).scale(7)'(assert) {
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
                new StraightEdge(new L3(V(0, 0, 7), V(0, 0.7071067811865475, -0.7071067811865475)), V(0, 3.0000000000000053, 3.9999999999999947), V(0, 3.9999999999999942, 3.0000000000000058), 4.242640687119293, 5.656854249492373)])], false)
        b2EqualAnd(assert, a, b, result)
    },

	'sphere(5) AND tetrahedron(V3.ZERO, V3.X, V3.Y, V3.Z.negated()).scale(7)'(assert) {
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
				new StraightEdge(new L3(V(-6, 0, 0), V(0.8320502943378436, -0.554700196225229, 0)), V(0, -4, 0), V(-4.950836318555471, -0.6994424542963525, 0), 7.211102550927979, 1.2609378166009386)], [])], false)
		b2EqualAnd(assert, a, b, result)
	},

	'box - doublebarrel'(assert) {
		const a = B2T.extrudeEdges([
			new StraightEdge(new L3(V(-409.7544757659938, 183.97231686845578, 0), V(0.9961411421023264, -0.08776573939227537, 0)), V(-409.7544757659938, 183.97231686845578, 0), V(-42.014475765994064, 151.57231686845583, 0), 0, 369.16455355301895),
			new StraightEdge(new L3(V(-42.014475765994064, 151.57231686845583, 0), V(-0.35632458233389064, 0.9343622381199803, 0)), V(-42.014475765994064, 151.57231686845583, 0), V(-114.91447576599401, 342.73231686845577, 0), 0, 204.5887474911559),
			new StraightEdge(new L3(V(-114.91447576599401, 342.73231686845577, 0), V(-0.9951217646853636, 0.09865431287829038, 0)), V(-114.91447576599401, 342.73231686845577, 0), V(-490.7544757659938, 379.99231686845576, 0), 0, 377.68242373719204),
			new StraightEdge(new L3(V(-490.7544757659938, 379.99231686845576, 0), V(0.3819019948313259, -0.9242028274918087, 0)), V(-490.7544757659938, 379.99231686845576, 0), V(-409.7544757659938, 183.97231686845578, 0), 0, 212.09629982628172)
		], new P3(V3.Z, 0), V(0, 0, -100), "box").minus(B2T.extrudeEdges([
			new PCurveEdge(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, 0),V(21.599999999999966, -23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(-358.31485781208175, 223.03687021651945, 0),V(-401.5148578120817, 269.83687021651934, 0),0,3.141592653589793,null,V(23.399999999999945, 21.59999999999997, 0),V(-23.39999999999995, -21.599999999999966, 0)),
			new PCurveEdge(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, 0),V(-21.599999999999966, 23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(-401.5148578120817, 269.83687021651934, 0),V(-358.31485781208175, 223.03687021651945, 0),0,3.141592653589793,null,V(-23.39999999999995, -21.599999999999966, 0),V(23.399999999999952, 21.599999999999962, 0))
		], new P3(V3.Z, 0), V(0, 0, -50), "cyl1"))
		const b = B2T.extrudeEdges([
			new PCurveEdge(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, 0),V(21.599999999999966, -23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(-312.10471643749133, 242.13196995808573, 0),V(-355.30471643749127, 288.9319699580857, 0),0,3.141592653589793,null,V(23.399999999999945, 21.59999999999997, 0),V(-23.39999999999995, -21.599999999999966, 0)),
			new PCurveEdge(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, 0),V(-21.599999999999966, 23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(-355.30471643749127, 288.9319699580857, 0),V(-312.10471643749133, 242.13196995808573, 0),0,3.141592653589793,null,V(-23.39999999999995, -21.599999999999966, 0),V(23.399999999999952, 21.599999999999962, 0))
		], new P3(V3.Z, 0), V(0, 0, -50), "cyl2").flipped()
		const result = new B2([
			new PlaneFace(new PlaneSurface(new P3(V3.Z, 0),V3.Y,V(-1, 0, 0)), [
				new StraightEdge(new L3(V(-409.7544757659938, 183.97231686845578, 0), V(0.9961411421023264, -0.08776573939227537, 0)), V(-409.7544757659938, 183.97231686845578, 0), V(-42.014475765994064, 151.57231686845583, 0), 0, 369.16455355301895),
				new StraightEdge(new L3(V(-42.014475765994064, 151.57231686845583, 0), V(-0.35632458233389064, 0.9343622381199803, 0)), V(-42.014475765994064, 151.57231686845583, 0), V(-114.91447576599401, 342.73231686845577, 0), 0, 204.5887474911559),
				new StraightEdge(new L3(V(-114.91447576599401, 342.73231686845577, 0), V(-0.9951217646853636, 0.09865431287829038, 0)), V(-114.91447576599401, 342.73231686845577, 0), V(-490.7544757659938, 379.99231686845576, 0), 0, 377.68242373719204),
				new StraightEdge(new L3(V(-490.7544757659938, 379.99231686845576, 0), V(0.3819019948313259, -0.9242028274918087, 0)), V(-490.7544757659938, 379.99231686845576, 0), V(-409.7544757659938, 183.97231686845578, 0), 0, 212.09629982628172)], [[
				new PCurveEdge(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, 0),V(-21.599999999999966, 23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(-355.30471643749127, 288.9319699580857, 0),V(-312.10471643749133, 242.13196995808573, 0),1.7763568394002505e-15,3.141592653589793,null,V(23.399999999999988, 21.599999999999923, 0),V(-23.399999999999945, -21.59999999999997, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, 0),V(21.599999999999966, -23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(-312.10471643749133, 242.13196995808573, 0),V(-349.27634070974653, 237.75347976290573, 0),0,1.256337123208411,null,V(-23.39999999999995, -21.599999999999966, 0),V(-27.77849019517994, 15.571624272255224, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, 0),V(-21.599999999999966, 23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(-349.27634070974653, 237.75347976290573, 0),V(-358.31485781208175, 223.03687021651945, 0),2.5923876672862365,3.141592653589793,null,V(-8.683390453613667, -30.63851710233521, 0),V(-23.399999999999945, -21.59999999999997, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, 0),V(21.599999999999966, -23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(-358.31485781208175, 223.03687021651945, 0),V(-401.5148578120817, 269.83687021651934, 0),0,3.1415926535897922,null,V(-23.39999999999995, -21.599999999999966, 0),V(23.399999999999928, 21.59999999999999, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, 0),V(-21.599999999999966, 23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(-401.5148578120817, 269.83687021651934, 0),V(-364.3432335398265, 274.2153604116993, 0),8.881784197001268e-16,1.2563371232084128,null,V(23.399999999999967, 21.599999999999945, 0),V(27.778490195179913, -15.571624272255272, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, 0),V(21.599999999999966, -23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(-364.3432335398265, 274.2153604116993, 0),V(-355.30471643749127, 288.9319699580857, 0),2.5923876672862347,3.1415926535897913,null,V(8.683390453613615, 30.63851710233523, 0),V(23.399999999999906, 21.600000000000012, 0))]]),
			new RotationFace(new SemiCylinderSurface(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, 0),V(21.599999999999966, -23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(0, 0, -1),-Infinity,Infinity), [
				new StraightEdge(new L3(V(-364.3432335398265, 274.21536041169935, 0), V3.Z), V(-364.3432335398265, 274.2153604116993, -50), V(-364.3432335398265, 274.2153604116993, 0), -50, 0),
				new PCurveEdge(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, 0),V(-21.599999999999966, 23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(-364.3432335398265, 274.2153604116993, 0),V(-401.5148578120817, 269.83687021651934, 0),1.2563371232084128,8.881784197001268e-16,null,V(-27.778490195179913, 15.571624272255272, 0),V(-23.399999999999967, -21.599999999999945, 0)),
				new StraightEdge(new L3(V(-401.5148578120817, 269.83687021651934, 0), V(0, 0, -1)), V(-401.5148578120817, 269.83687021651934, 0), V(-401.5148578120817, 269.83687021651934, -50), 0, 50),
				new PCurveEdge(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, -50),V(21.599999999999966, -23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(-401.5148578120817, 269.83687021651934, -50),V(-364.3432335398265, 274.2153604116993, -50),3.141592653589793,1.8852555303813805,null,V(23.39999999999995, 21.599999999999966, 0),V(27.778490195179916, -15.57162427225527, 0))], []),
			new RotationFace(new SemiCylinderSurface(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, 0),V(21.599999999999966, -23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(0, 0, -1),-Infinity,Infinity), [
				new StraightEdge(new L3(V(-349.27634070974653, 237.75347976290573, 0), V(0, 0, -1)), V(-349.27634070974653, 237.75347976290573, 0), V(-349.27634070974653, 237.75347976290573, -50), 0, 50),
				new PCurveEdge(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, -50),V(21.599999999999966, -23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(-349.27634070974653, 237.75347976290573, -50),V(-358.31485781208175, 223.03687021651945, -50),0.5492049863035567,0,null,V(-8.683390453613667, -30.63851710233521, 0),V(-23.399999999999945, -21.59999999999997, 0)),
				new StraightEdge(new L3(V(-358.31485781208175, 223.03687021651945, 0), V(0, 0, -1)), V(-358.31485781208175, 223.03687021651945, -50), V(-358.31485781208175, 223.03687021651945, 0), 50, 0),
				new PCurveEdge(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, 0),V(-21.599999999999966, 23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(-358.31485781208175, 223.03687021651945, 0),V(-349.27634070974653, 237.75347976290573, 0),3.141592653589793,2.5923876672862365,null,V(23.399999999999945, 21.59999999999997, 0),V(8.683390453613667, 30.63851710233521, 0))], []),
			new PlaneFace(new PlaneSurface(new P3(V3.Z, -50),V3.Y,V(-1, 0, 0)), [
				new PCurveEdge(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, -50),V(21.599999999999966, -23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(-349.2763407097464, 237.75347976290567, -50),V(-364.34323353982654, 274.2153604116993, -50),0.5492049863035546,1.8852555303813814,null,V(8.683390453613734, 30.638517102335193, 0),V(-27.778490195179927, 15.571624272255246, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, -50),V(21.599999999999966, -23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(-364.3432335398265, 274.2153604116993, -50),V(-349.27634070974653, 237.75347976290573, -50),2.5923876672862347,1.256337123208411,null,V(-8.683390453613615, -30.63851710233523, 0),V(27.77849019517994, -15.571624272255224, 0))], []),
			new PlaneFace(new PlaneSurface(new P3(V3.Z, -50),V3.Y,V(-1, 0, 0)), [
				new PCurveEdge(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, -50),V(21.599999999999966, -23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(-349.27634070974653, 237.75347976290573, -50),V(-364.3432335398265, 274.2153604116993, -50),1.256337123208411,2.5923876672862347,null,V(-27.77849019517994, 15.571624272255224, 0),V(8.683390453613615, 30.63851710233523, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, -50),V(21.599999999999966, -23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(-364.3432335398265, 274.2153604116993, -50),V(-401.5148578120817, 269.83687021651934, -50),1.8852555303813805,3.141592653589793,null,V(-27.778490195179916, 15.57162427225527, 0),V(-23.39999999999995, -21.599999999999966, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, -50),V(-21.599999999999966, 23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(-401.5148578120817, 269.83687021651934, -50),V(-358.31485781208175, 223.03687021651945, -50),0,3.141592653589793,null,V(-23.39999999999995, -21.599999999999966, 0),V(23.399999999999952, 21.599999999999962, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, -50),V(21.599999999999966, -23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(-358.31485781208175, 223.03687021651945, -50),V(-349.27634070974653, 237.75347976290573, -50),0,0.5492049863035567,null,V(23.399999999999945, 21.59999999999997, 0),V(8.683390453613667, 30.63851710233521, 0))], []),
			new PlaneFace(new PlaneSurface(new P3(V(-0.08776573939227537, -0.9961411421023264, 0), -147.29998930565802),V(0.9961411421023264, -0.08776573939227537, 0),V3.Z), [
				new StraightEdge(new L3(V(-409.7544757659938, 183.97231686845578, 0), V(0.9961411421023264, -0.08776573939227537, 0)), V(-42.014475765994064, 151.57231686845583, 0), V(-409.7544757659938, 183.97231686845578, 0), 369.16455355301895, 0),
				new StraightEdge(new L3(V(-409.7544757659938, 183.97231686845578, 0), V(0, 0, -1)), V(-409.7544757659938, 183.97231686845578, 0), V(-409.7544757659938, 183.97231686845578, -100), 0, 100),
				new StraightEdge(new L3(V(-409.7544757659938, 183.97231686845578, -100), V(0.9961411421023264, -0.08776573939227537, 0)), V(-409.7544757659938, 183.97231686845578, -100), V(-42.014475765994064, 151.57231686845583, -100), 0, 369.16455355301895),
				new StraightEdge(new L3(V(-42.014475765994064, 151.57231686845583, 0), V(0, 0, -1)), V(-42.014475765994064, 151.57231686845583, -100), V(-42.014475765994064, 151.57231686845583, 0), 100, 0)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0.9343622381199802, 0.3563245823338906, 0), 14.752202891380762),V(-0.3563245823338906, 0.9343622381199802, 0),V3.Z), [
				new StraightEdge(new L3(V(-42.014475765994064, 151.57231686845583, 0), V(-0.35632458233389064, 0.9343622381199803, 0)), V(-114.91447576599401, 342.73231686845577, 0), V(-42.014475765994064, 151.57231686845583, 0), 204.5887474911559, 0),
				new StraightEdge(new L3(V(-42.014475765994064, 151.57231686845583, 0), V(0, 0, -1)), V(-42.014475765994064, 151.57231686845583, 0), V(-42.014475765994064, 151.57231686845583, -100), 0, 100),
				new StraightEdge(new L3(V(-42.014475765994064, 151.57231686845583, -100), V(-0.35632458233389064, 0.9343622381199803, 0)), V(-42.014475765994064, 151.57231686845583, -100), V(-114.91447576599401, 342.73231686845577, -100), 0, 204.5887474911559),
				new StraightEdge(new L3(V(-114.91447576599401, 342.73231686845577, 0), V(0, 0, -1)), V(-114.91447576599401, 342.73231686845577, -100), V(-114.91447576599401, 342.73231686845577, 0), 100, 0)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0.09865431287829038, 0.9951217646853637, 0), 329.72357933037785),V(-0.9951217646853637, 0.09865431287829038, 0),V3.Z), [
				new StraightEdge(new L3(V(-114.91447576599401, 342.73231686845577, 0), V(-0.9951217646853636, 0.09865431287829038, 0)), V(-490.7544757659938, 379.99231686845576, 0), V(-114.91447576599401, 342.73231686845577, 0), 377.68242373719204, 0),
				new StraightEdge(new L3(V(-114.91447576599401, 342.73231686845577, 0), V(0, 0, -1)), V(-114.91447576599401, 342.73231686845577, 0), V(-114.91447576599401, 342.73231686845577, -100), 0, 100),
				new StraightEdge(new L3(V(-114.91447576599401, 342.73231686845577, -100), V(-0.9951217646853636, 0.09865431287829038, 0)), V(-114.91447576599401, 342.73231686845577, -100), V(-490.7544757659938, 379.99231686845576, -100), 0, 377.68242373719204),
				new StraightEdge(new L3(V(-490.7544757659938, 379.99231686845576, 0), V(0, 0, -1)), V(-490.7544757659938, 379.99231686845576, -100), V(-490.7544757659938, 379.99231686845576, 0), 100, 0)], []),
			new PlaneFace(new PlaneSurface(new P3(V(-0.9242028274918087, -0.381901994831326, 0), 308.4368502745512),V(0.381901994831326, -0.9242028274918087, 0),V3.Z), [
				new StraightEdge(new L3(V(-490.7544757659938, 379.99231686845576, 0), V(0.3819019948313259, -0.9242028274918087, 0)), V(-409.7544757659938, 183.97231686845578, 0), V(-490.7544757659938, 379.99231686845576, 0), 212.09629982628172, 0),
				new StraightEdge(new L3(V(-490.7544757659938, 379.99231686845576, 0), V(0, 0, -1)), V(-490.7544757659938, 379.99231686845576, 0), V(-490.7544757659938, 379.99231686845576, -100), 0, 100),
				new StraightEdge(new L3(V(-490.7544757659938, 379.99231686845576, -100), V(0.3819019948313259, -0.9242028274918087, 0)), V(-490.7544757659938, 379.99231686845576, -100), V(-409.7544757659938, 183.97231686845578, -100), 0, 212.09629982628172),
				new StraightEdge(new L3(V(-409.7544757659938, 183.97231686845578, 0), V(0, 0, -1)), V(-409.7544757659938, 183.97231686845578, -100), V(-409.7544757659938, 183.97231686845578, 0), 100, 0)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 100),V3.Y,V3.X), [
				new StraightEdge(new L3(V(-490.7544757659938, 379.99231686845576, -100), V(0.3819019948313259, -0.9242028274918087, 0)), V(-409.7544757659938, 183.97231686845578, -100), V(-490.7544757659938, 379.99231686845576, -100), 212.09629982628172, 0),
				new StraightEdge(new L3(V(-114.91447576599401, 342.73231686845577, -100), V(-0.9951217646853636, 0.09865431287829038, 0)), V(-490.7544757659938, 379.99231686845576, -100), V(-114.91447576599401, 342.73231686845577, -100), 377.68242373719204, 0),
				new StraightEdge(new L3(V(-42.014475765994064, 151.57231686845583, -100), V(-0.35632458233389064, 0.9343622381199803, 0)), V(-114.91447576599401, 342.73231686845577, -100), V(-42.014475765994064, 151.57231686845583, -100), 204.5887474911559, 0),
				new StraightEdge(new L3(V(-409.7544757659938, 183.97231686845578, -100), V(0.9961411421023264, -0.08776573939227537, 0)), V(-42.014475765994064, 151.57231686845583, -100), V(-409.7544757659938, 183.97231686845578, -100), 369.16455355301895, 0)], []),
			new RotationFace(new SemiCylinderSurface(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, 0),V(-21.599999999999966, 23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(0, 0, -1),-Infinity,Infinity), [
				new PCurveEdge(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, 0),V(21.599999999999966, -23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(-401.5148578120817, 269.83687021651934, 0),V(-358.31485781208175, 223.03687021651945, 0),3.1415926535897922,0,null,V(-23.399999999999928, -21.59999999999999, 0),V(23.39999999999995, 21.599999999999966, 0)),
				new StraightEdge(new L3(V(-358.31485781208175, 223.03687021651945, 0), V(0, 0, -1)), V(-358.31485781208175, 223.03687021651945, 0), V(-358.31485781208175, 223.03687021651945, -50), 0, 50),
				new PCurveEdge(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, -50),V(-21.599999999999966, 23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(-358.31485781208175, 223.03687021651945, -50),V(-401.5148578120817, 269.83687021651934, -50),3.141592653589793,0,null,V(-23.399999999999952, -21.599999999999962, 0),V(23.39999999999995, 21.599999999999966, 0)),
				new StraightEdge(new L3(V(-401.5148578120817, 269.83687021651934, 0), V(0, 0, -1)), V(-401.5148578120817, 269.83687021651934, -50), V(-401.5148578120817, 269.83687021651934, 0), 50, 0)], []),
			new RotationFace(new SemiCylinderSurface(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, 0),V(21.599999999999966, -23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(0, 0, -1),-Infinity,Infinity), [
				new PCurveEdge(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, 0),V(-21.599999999999966, 23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(-312.10471643749133, 242.13196995808573, 0),V(-355.30471643749127, 288.9319699580857, 0),3.141592653589793,1.7763568394002505e-15,null,V(23.399999999999945, 21.59999999999997, 0),V(-23.399999999999988, -21.599999999999923, 0)),
				new StraightEdge(new L3(V(-355.30471643749127, 288.9319699580857, 0), V(0, 0, -1)), V(-355.30471643749127, 288.9319699580857, 0), V(-355.30471643749127, 288.9319699580857, -50), 0, 50),
				new PCurveEdge(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, -50),V(21.599999999999966, -23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(-355.30471643749127, 288.9319699580857, -50),V(-312.10471643749133, 242.13196995808573, -50),3.141592653589793,0,null,V(23.39999999999995, 21.599999999999966, 0),V(-23.399999999999945, -21.59999999999997, 0)),
				new StraightEdge(new L3(V(-312.10471643749133, 242.13196995808573, 0), V(0, 0, -1)), V(-312.10471643749133, 242.13196995808573, -50), V(-312.10471643749133, 242.13196995808573, 0), 50, 0)], []),
			new RotationFace(new SemiCylinderSurface(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, 0),V(-21.599999999999966, 23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(0, 0, -1),-Infinity,Infinity), [
				new PCurveEdge(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, 0),V(21.599999999999966, -23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(-349.27634070974653, 237.75347976290573, 0),V(-312.10471643749133, 242.13196995808573, 0),1.256337123208411,0,null,V(27.77849019517994, -15.571624272255224, 0),V(23.39999999999995, 21.599999999999966, 0)),
				new StraightEdge(new L3(V(-312.10471643749133, 242.13196995808573, 0), V(0, 0, -1)), V(-312.10471643749133, 242.13196995808573, 0), V(-312.10471643749133, 242.13196995808573, -50), 0, 50),
				new PCurveEdge(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, -50),V(-21.599999999999966, 23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(-312.10471643749133, 242.13196995808573, -50),V(-349.2763407097464, 237.75347976290567, -50),3.141592653589793,1.885255530381386,null,V(-23.399999999999952, -21.599999999999962, 0),V(-27.778490195180005, 15.571624272255118, 0)),
				new StraightEdge(new L3(V(-349.27634070974653, 237.75347976290573, 0), V(0, 0, -1)), V(-349.27634070974653, 237.75347976290573, -50), V(-349.27634070974653, 237.75347976290573, 0), 50, 0)], []),
			new RotationFace(new SemiCylinderSurface(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, 0),V(-21.599999999999966, 23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(0, 0, -1),-Infinity,Infinity), [
				new PCurveEdge(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, 0),V(21.599999999999966, -23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(-355.30471643749127, 288.9319699580857, 0),V(-364.3432335398265, 274.2153604116993, 0),3.1415926535897913,2.5923876672862347,null,V(-23.399999999999906, -21.600000000000012, 0),V(-8.683390453613615, -30.63851710233523, 0)),
				new StraightEdge(new L3(V(-364.3432335398265, 274.21536041169935, 0), V3.Z), V(-364.3432335398265, 274.2153604116993, 0), V(-364.3432335398265, 274.2153604116993, -50), 0, -50),
				new PCurveEdge(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, -50),V(-21.599999999999966, 23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(-364.34323353982654, 274.2153604116993, -50),V(-355.30471643749127, 288.9319699580857, -50),0.5492049863035583,0,null,V(8.68339045361362, 30.638517102335225, 0),V(23.39999999999995, 21.599999999999966, 0)),
				new StraightEdge(new L3(V(-355.30471643749127, 288.9319699580857, 0), V(0, 0, -1)), V(-355.30471643749127, 288.9319699580857, -50), V(-355.30471643749127, 288.9319699580857, 0), 50, 0)], []),
			new PlaneFace(new PlaneSurface(new P3(V3.Z, -50),V3.Y,V(-1, 0, 0)), [
				new PCurveEdge(new SemiEllipseCurve(V(-379.9148578120817, 246.4368702165194, -50),V(21.599999999999966, -23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(-364.34323353982654, 274.2153604116993, -50),V(-349.2763407097464, 237.75347976290567, -50),1.8852555303813814,0.5492049863035546,null,V(27.778490195179927, -15.571624272255246, 0),V(-8.683390453613734, -30.638517102335193, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, -50),V(-21.599999999999966, 23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(-349.2763407097464, 237.75347976290567, -50),V(-312.10471643749133, 242.13196995808573, -50),1.885255530381386,3.141592653589793,null,V(27.778490195180005, -15.571624272255118, 0),V(23.399999999999952, 21.599999999999962, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, -50),V(21.599999999999966, -23.39999999999995, 0),V(23.39999999999995, 21.599999999999966, 0),0,3.141592653589793),V(-312.10471643749133, 242.13196995808573, -50),V(-355.30471643749127, 288.9319699580857, -50),0,3.141592653589793,null,V(23.399999999999945, 21.59999999999997, 0),V(-23.39999999999995, -21.599999999999966, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(-333.7047164374913, 265.5319699580857, -50),V(-21.599999999999966, 23.39999999999995, 0),V(-23.39999999999995, -21.599999999999966, 0),0,3.141592653589793),V(-355.30471643749127, 288.9319699580857, -50),V(-364.34323353982654, 274.2153604116993, -50),0,0.5492049863035583,null,V(-23.39999999999995, -21.599999999999966, 0),V(-8.68339045361362, -30.638517102335225, 0))], [])], false)
		b2EqualAnd(assert, a, b, result)
	},

    'sphere() - cubelet near -V3.Y'(assert) {
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
                new StraightEdge(new L3(V(0.2, -0.8874637410145799, 0.3389809852844352), V(-1, 0, 0)), V(0.2, -0.8874637410145799, 0.3389809852844352), V(0, -0.8874637410145799, 0.3389809852844352), 0, 0.2)], [])], false)
        b2EqualAnd(assert, a, b, result)
    },

    'sphere() - B2 w/ PCS'(assert) {
        const a = B2T.sphere()
        const b = B2T.extrudeEdges([PCurveEdge.forCurveAndTs(BezierCurve.QUARTER_CIRCLE, 0, 1), StraightEdge.throughPoints(V3.Y, V3.X)], P3.XY, V3.Z.negated())
            .scale(0.2, 0.2, 2)
            .rotateX(85 * DEG)
            .translate(0.1, 0.1, 0.4)
            .flipped()
        const result = new B2([
	        new RotationFace(new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), [
		        new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0, 0, 1), V(0, 0, -1), 3.141592653589793, 0, null, V(1, 0, 0), V(-1, 0, 0)),
		        new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(0, 0, -1), V(0, 0, 1), 0, 3.141592653589793, null, V(-1, 1.2246467991473532e-16, 0), V(1, -1.2246467991473532e-16, 0))], [[
		        new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, 0.1, 0.4), V(0.30000000000000004, 0.10962695751594148, 0.5100366279236725), V(0.2104569499661587, 0.11743114854953163, 0.5992389396183492), V(0.10000000000000002, 0.11743114854953163, 0.5992389396183492), 0, 1), V(0, 0.9961946980917455, -0.08715574274765814), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.30000000000000004, 0.8948730420942135, 0.3304576198743918), V(0.10000000000000002, 0.8381569952434025, 0.5361835985225125), -1), V(0.30000000000000004, 0.8948730420942135, 0.33045761987439176), V(0.1, 0.8381569952434027, 0.5361835985225124), 0, 100, null, V(-1.2097686775990997e-18, -0.001254118535048534, 0.0033961294916794536), V(-0.003309386532465258, 0.00041824948412341, -0.0000365920883863643)),
		        new PCurveEdge(new SemiEllipseCurve(V(0.35359672675573195, 0.030817985353536543, 0.3522511844566559), V(0.03780987199589104, -0.8643439107956781, 0.037665994017834215), V(-0.6111789986634738, 0, 0.613513603148274), 1.5351339957098988, 3.141592653589793), V(0.1, 0.8381569952434027, 0.5361835985225124), V(0.30000000000000004, 0.8948730420942135, 0.33045761987439176), 2.7763814662109167, 3.1157389336980312, null, V(0.5573670738448863, 0.3086974558359978, -0.586503653363929), V(0.6099973325422855, 0.022344015987121178, -0.6142822713509873))]]),
	        new RotationFace(new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), [
		        new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(0, 0, 1), V(0, 0, -1), 3.141592653589793, 0, null, V(-1, 1.2246467991473532e-16, 0), V(1, -1.2246467991473532e-16, 0)),
		        new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0, 0, -1), V(0, 0, 1), 0, 3.141592653589793, null, V(1, 0, 0), V(-1, 0, 0))], []),
	        new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, 0.1, 0.4), V(0.30000000000000004, 0.10962695751594148, 0.5100366279236725), V(0.2104569499661587, 0.11743114854953163, 0.5992389396183492), V(0.10000000000000002, 0.11743114854953163, 0.5992389396183492), 0, 1), V(0, 0.9961946980917455, -0.08715574274765814), 0, 1, -Infinity, Infinity), [
		        new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, 0.1, 0.4), V(0.30000000000000004, 0.10962695751594148, 0.5100366279236725), V(0.2104569499661587, 0.11743114854953163, 0.5992389396183492), V(0.10000000000000002, 0.11743114854953163, 0.5992389396183492), 0, 1), V(0, 0.9961946980917455, -0.08715574274765814), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.30000000000000004, 0.8948730420942135, 0.3304576198743918), V(0.10000000000000002, 0.8381569952434025, 0.5361835985225125), -1), V(0.1, 0.8381569952434027, 0.5361835985225124), V(0.30000000000000004, 0.8948730420942135, 0.33045761987439176), 100, 0, null, V(0.003309386532465258, -0.00041824948412341, 0.0000365920883863643), V(1.2097686775990997e-18, 0.001254118535048534, -0.0033961294916794536)),
		        new StraightEdge(new L3(V(0.30000000000000004, 0.1, 0.4), V(0, 0.9961946980917455, -0.08715574274765814)), V(0.30000000000000004, 0.8948730420942135, 0.33045761987439176), V(0.30000000000000004, 0.1, 0.4), 0.7979093279825997, 0),
		        new PCurveEdge(new BezierCurve(V(0.30000000000000004, 0.1, 0.4), V(0.30000000000000004, 0.10962695751594148, 0.5100366279236725), V(0.2104569499661587, 0.11743114854953163, 0.5992389396183492), V(0.10000000000000002, 0.11743114854953163, 0.5992389396183492), 0, 1), V(0.30000000000000004, 0.1, 0.4), V(0.10000000000000002, 0.11743114854953163, 0.5992389396183492), 0, 1, null, V(0, 0.028880872547824416, 0.33010988377101746), V(-0.33137084989847604, 0, 0)),
		        new StraightEdge(new L3(V(0.1, 0.11743114854953163, 0.5992389396183492), V(0, 0.9961946980917455, -0.08715574274765814)), V(0.1, 0.11743114854953163, 0.5992389396183492), V(0.1, 0.8381569952434027, 0.5361835985225124), 0, 0.7234789023415331)], []),
	        new PlaneFace(new PlaneSurface(new P3(V(0.7071067811865475, 0.06162841671621931, 0.7044160264027587), 0.5000612865886896), V(0.701755757082144, 0.06116204423593941, -0.7097873355780223), V(-0.08682659386424756, 0.9962234400966147, 0)), [
		        new PCurveEdge(new SemiEllipseCurve(V(0.35359672675573195, 0.030817985353536543, 0.3522511844566559), V(0.03780987199589104, -0.8643439107956781, 0.037665994017834215), V(-0.6111789986634738, 0, 0.613513603148274), 1.5351339957098988, 3.141592653589793), V(0.30000000000000004, 0.8948730420942135, 0.33045761987439176), V(0.1, 0.8381569952434027, 0.5361835985225124), 3.1157389336980312, 2.7763814662109167, null, V(-0.6099973325422855, -0.022344015987121178, 0.6142822713509873), V(-0.5573670738448863, -0.3086974558359978, 0.586503653363929)),
		        new StraightEdge(new L3(V(0.1, 0.11743114854953163, 0.5992389396183492), V(0, 0.9961946980917455, -0.08715574274765814)), V(0.1, 0.8381569952434027, 0.5361835985225124), V(0.1, 0.11743114854953163, 0.5992389396183492), 0.7234789023415331, 0),
		        new StraightEdge(new L3(V(0.1, 0.11743114854953163, 0.5992389396183492), V(0.7071067811865476, -0.06162841671621934, -0.7044160264027587)), V(0.1, 0.11743114854953163, 0.5992389396183492), V(0.30000000000000004, 0.1, 0.4), 0, 0.282842712474619),
		        new StraightEdge(new L3(V(0.30000000000000004, 0.1, 0.4), V(0, 0.9961946980917455, -0.08715574274765814)), V(0.30000000000000004, 0.1, 0.4), V(0.30000000000000004, 0.8948730420942135, 0.33045761987439176), 0, 0.7979093279825997)], []),
	        new PlaneFace(new PlaneSurface(new P3(V(0, 0.9961946980917455, -0.08715574274765814), 0.0647571727101113), V(0, -0.08715574274765814, -0.9961946980917455), V(-1, 0, 0)), [
		        new StraightEdge(new L3(V(0.1, 0.11743114854953163, 0.5992389396183492), V(0.7071067811865476, -0.06162841671621934, -0.7044160264027587)), V(0.30000000000000004, 0.1, 0.4), V(0.1, 0.11743114854953163, 0.5992389396183492), 0.282842712474619, 0),
		        new PCurveEdge(new BezierCurve(V(0.30000000000000004, 0.1, 0.4), V(0.30000000000000004, 0.10962695751594148, 0.5100366279236725), V(0.2104569499661587, 0.11743114854953163, 0.5992389396183492), V(0.10000000000000002, 0.11743114854953163, 0.5992389396183492), 0, 1), V(0.10000000000000002, 0.11743114854953163, 0.5992389396183492), V(0.30000000000000004, 0.1, 0.4), 1, 0, null, V(0.33137084989847604, 0, 0), V(0, -0.028880872547824416, -0.33010988377101746))], [])], false)
        b2EqualAnd(assert, a, b, result)
    },

    'sphere() - B2 w/ PCS 2'(assert) {
        const a = B2T.sphere()
        const b = B2T.extrudeEdges([PCurveEdge.forCurveAndTs(BezierCurve.QUARTER_CIRCLE, 0, 1), StraightEdge.throughPoints(V3.Y, V3.X)], P3.XY, V3.Z.negated())
            .scale(0.2, 0.2, 2)
            .translate(0.1, -0.1, 1.2)
            .flipped()
        const result = new B2([
	        new RotationFace(new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), [
		        new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.10000000000000002, 0.1, 0.9899494936611666), V(0.2732435609477567, -6.071532165918825e-18, 0.9619448822051031), -1), V(0.1, 0.1, 0.9899494936611666), V(0.2732435609477568, 0, 0.961944882205103), 0, 67, null, V(0.003311853592333425, 0, -0.00033454773334799895), V(0.0015464506679683772, -0.0026733520656420703, -0.0004392743234696266)),
		        new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0.2732435609477568, 0, 0.961944882205103), V(0, 0, -1), 2.8648293508136553, 0, null, V(0.961944882205103, 0, -0.2732435609477568), V(-1, 0, 0)),
		        new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(0, 0, -1), V(0, 0, 1), 0, 3.141592653589793, null, V(-1, 1.2246467991473532e-16, 0), V(1, -1.2246467991473532e-16, 0)),
		        new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0, 0, 1), V(0.19999999999999976, 0, 0.9797958971132713), 3.141592653589793, 2.9402347327994627, null, V(1, 0, 0), V(0.9797958971132713, 0, -0.19999999999999976)),
		        new PCurveEdge(new SemiEllipseCurve(V(0.09999999999999994, 0.09999999999999995, 0), V(0.7000000000000001, -0.7, 0), V(0, 0, 0.9899494936611666), 1.4274487578895312, 3.141592653589793), V(0.19999999999999976, 0, 0.9797958971132713), V(0.1, 0.1, 0.9899494936611666), 1.4274487578895314, 1.5707963267948966, null, V(-0.6928203230275509, 0.6928203230275509, 0.14142135623730936), V(-0.7000000000000001, 0.7, 6.061692393648453e-17))], []),
	        new RotationFace(new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), [
		        new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.27324356094775665, -3.2161642623246884e-17, 0.9619448822051031), V(0.30000000000000004, -0.1, 0.9486832980505138), -1), V(0.2732435609477568, 0, 0.961944882205103), V(0.30000000000000004, -0.1, 0.9486832980505138), 0, 33, null, V(0.0015464443998091284, -0.002673341229863844, -0.00043927254297868013), V(0, -0.003311691350426493, -0.00034908291916088506)),
		        new PCurveEdge(new SemiEllipseCurve(V(0.09999999999999994, 0.09999999999999995, 0), V(0, 0, -0.9899494936611666), V(0.7000000000000001, -0.7, 0), 0.14334756890536535, 2.9982450846844277), V(0.30000000000000004, -0.1, 0.9486832980505138), V(0.19999999999999976, 0, 0.9797958971132713), 2.851840952153746, 2.9982450846844277, null, V(-0.670820393249937, 0.6708203932499369, 0.282842712474619), V(-0.6928203230275509, 0.6928203230275509, 0.1414213562373092)),
		        new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0.19999999999999976, 0, 0.9797958971132713), V(0, 0, 1), 2.9402347327994627, 3.141592653589793, null, V(-0.9797958971132713, 0, 0.19999999999999976), V(-1, 0, 0)),
		        new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(0, 0, 1), V(0, 0, -1), 3.141592653589793, 0, null, V(-1, 1.2246467991473532e-16, 0), V(1, -1.2246467991473532e-16, 0)),
		        new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0, 0, -1), V(0.2732435609477568, 0, 0.961944882205103), 0, 2.8648293508136553, null, V(1, 0, 0), V(-0.961944882205103, 0, 0.2732435609477568))], []),
	        new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
		        new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.10000000000000002, 0.1, 0.9899494936611666), V(0.2732435609477567, -6.071532165918825e-18, 0.9619448822051031), -1), V(0.2732435609477568, 0, 0.961944882205103), V(0.1, 0.1, 0.9899494936611666), 67, 0, null, V(-0.0015464506679683772, 0.0026733520656420703, 0.0004392743234696266), V(-0.003311853592333425, 0, 0.00033454773334799895)),
		        new StraightEdge(new L3(V(0.1, 0.1, 1.2), V(0, 0, -1)), V(0.1, 0.1, 0.9899494936611666), V(0.1, 0.1, -0.8), 0.21005050633883338, 2),
		        new PCurveEdge(new BezierCurve(V(0.30000000000000004, -0.1, -0.8), V(0.30000000000000004, 0.010456949966158674, -0.8), V(0.2104569499661587, 0.1, -0.8), V(0.10000000000000002, 0.1, -0.8), 0, 1), V(0.10000000000000002, 0.1, -0.8), V(0.30000000000000004, -0.1, -0.8), 1, 0, null, V(0.33137084989847604, 0, 0), V(0, -0.33137084989847604, 0)),
		        new StraightEdge(new L3(V(0.30000000000000004, -0.1, 1.2), V(0, 0, -1)), V(0.30000000000000004, -0.1, -0.8), V(0.30000000000000004, -0.1, 0.9486832980505138), 2, 0.2513167019494862),
		        new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.27324356094775665, -3.2161642623246884e-17, 0.9619448822051031), V(0.30000000000000004, -0.1, 0.9486832980505138), -1), V(0.30000000000000004, -0.1, 0.9486832980505138), V(0.2732435609477568, 0, 0.961944882205103), 33, 0, null, V(0, 0.003311691350426493, 0.00034908291916088506), V(-0.0015464443998091284, 0.002673341229863844, 0.00043927254297868013))], []),
	        new PlaneFace(new PlaneSurface(new P3(V(0.7071067811865475, 0.7071067811865476, 0), 0.14142135623730942), V(0, 0, -1), V(-0.7071067811865476, 0.7071067811865475, 0)), [
		        new PCurveEdge(new SemiEllipseCurve(V(0.09999999999999994, 0.09999999999999995, 0), V(0.7000000000000001, -0.7, 0), V(0, 0, 0.9899494936611666), 1.4274487578895312, 3.141592653589793), V(0.1, 0.1, 0.9899494936611666), V(0.19999999999999976, 0, 0.9797958971132713), 1.5707963267948966, 1.4274487578895314, null, V(0.7000000000000001, -0.7, -6.061692393648453e-17), V(0.6928203230275509, -0.6928203230275509, -0.14142135623730936)),
		        new PCurveEdge(new SemiEllipseCurve(V(0.09999999999999994, 0.09999999999999995, 0), V(0, 0, -0.9899494936611666), V(0.7000000000000001, -0.7, 0), 0.14334756890536535, 2.9982450846844277), V(0.19999999999999976, 0, 0.9797958971132713), V(0.30000000000000004, -0.1, 0.9486832980505138), 2.9982450846844277, 2.851840952153746, null, V(0.6928203230275509, -0.6928203230275509, -0.1414213562373092), V(0.670820393249937, -0.6708203932499369, -0.282842712474619)),
		        new StraightEdge(new L3(V(0.30000000000000004, -0.1, 1.2), V(0, 0, -1)), V(0.30000000000000004, -0.1, 0.9486832980505138), V(0.30000000000000004, -0.1, -0.8), 0.2513167019494862, 2),
		        new StraightEdge(new L3(V(0.1, 0.1, -0.8), V(0.7071067811865475, -0.7071067811865475, 0)), V(0.30000000000000004, -0.1, -0.8), V(0.1, 0.1, -0.8), 0.28284271247461906, 0),
		        new StraightEdge(new L3(V(0.1, 0.1, 1.2), V(0, 0, -1)), V(0.1, 0.1, -0.8), V(0.1, 0.1, 0.9899494936611666), 2, 0.21005050633883338)], []),
	        new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), -0.8), V(-1, 0, 0), V(0, -1, 0)), [
		        new PCurveEdge(new BezierCurve(V(0.30000000000000004, -0.1, -0.8), V(0.30000000000000004, 0.010456949966158674, -0.8), V(0.2104569499661587, 0.1, -0.8), V(0.10000000000000002, 0.1, -0.8), 0, 1), V(0.30000000000000004, -0.1, -0.8), V(0.10000000000000002, 0.1, -0.8), 0, 1, null, V(0, 0.33137084989847604, 0), V(-0.33137084989847604, 0, 0)),
		        new StraightEdge(new L3(V(0.1, 0.1, -0.8), V(0.7071067811865475, -0.7071067811865475, 0)), V(0.1, 0.1, -0.8), V(0.30000000000000004, -0.1, -0.8), 0, 0.28284271247461906)], [])], false)
        b2EqualAnd(assert, a, b, result)
    },

    'box() - B2 w/ PCS 2'(assert) {
        const a = B2T.box()
        const b = B2T.extrudeEdges([PCurveEdge.forCurveAndTs(BezierCurve.QUARTER_CIRCLE, 0, 1), StraightEdge.throughPoints(V3.Y, V3.X)], P3.XY, V3.Z.negated())
            .scale(0.2, 0.2, 4)
	        .translate(-0.1,0.4,2)
	        .rotateX(10*DEG)
	        .flipped()
        const result = new B2([
	        new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), 0), V(0, 0, -1), V(0, -1, 0)), [
		        new StraightEdge(new L3(V(3.3306690738754696e-16, 0.21723834785181606, 2.0691582057422955), V(0, 0.17364817766693036, -0.9848077530122082)), V(0, 0.582086766878481, 0), V(0, 0.40575978617006886, 1), 2.1010783063124836, 1.0856516944267476),
		        new StraightEdge(new L3(V(0, 0, 1), V(0, 1, 0)), V(0, 0.40575978617006886, 1), V(0, 1, 1), 0.40575978617006886, 1),
		        new StraightEdge(new L3(V(0, 1, 0), V(0, 0, 1)), V(0, 1, 1), V(0, 1, 0), 1, 0),
		        new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 1, 0), V(0, 0.582086766878481, 0), 1, 0.582086766878481)], []),
	        new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), 0), V(0, 0, -1), V(0, -1, 0)), [
		        new StraightEdge(new L3(V(0, 0.49240387650610384, 0.08682408883346515), V(0, 0.1736481776669304, -0.9848077530122082)), V(0, 0.33138632523440736, 1), V(0, 0.5077133059428723, 0), -0.9272631215315127, 0.08816349035423249),
		        new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 0.5077133059428723, 0), V(0, 0, 0), 0.5077133059428723, 0),
		        new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 0), V(0, 0, 1), 0, 1),
		        new StraightEdge(new L3(V(0, 0, 1), V(0, 1, 0)), V(0, 0, 1), V(0, 0.33138632523440736, 1), 0, 0.33138632523440736)], []),
	        new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 0), V(-1, 0, 0), V(0, 1, 0)), [
		        new PCurveEdge(new BezierCurve(V(-0.09999999999999999, 0.6092559671314469, 0), V(0.010456949966158688, 0.6092559671314469, 0), V(0.1, 0.5183315712176677, 0), V(0.1, 0.4061706447542979, 0), 0, 1), V(0.1, 0.406170644754298, 0), V(0, 0.582086766878481, 0), 0.9999999999999998, 0.3298001108864301, null, V(-1.1929530701127963e-16, 0.3364827793901093, 0), V(-0.26759242130902194, 0.15718180225447817, 0)),
		        new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 0.582086766878481, 0), V(0, 1, 0), 0.582086766878481, 1),
		        new StraightEdge(new L3(V(0, 1, 0), V(1, 0, 0)), V(0, 1, 0), V(1, 1, 0), 0, 1),
		        new StraightEdge(new L3(V(1, 1, 0), V(0, -1, 0)), V(1, 1, 0), V(1, 0, 0), 0, 1),
		        new StraightEdge(new L3(V(1, 0, 0), V(-1, 0, 0)), V(1, 0, 0), V(0, 0, 0), 0, 1),
		        new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 0, 0), V(0, 0.5077133059428723, 0), 0, 0.5077133059428723),
		        new StraightEdge(new L3(V(0.2538269089126336, 0.24997070782028524, 0), V(0.7016738431598636, -0.7124982932086696, 0)), V(0, 0.5077133059428723, 0), V(0.1, 0.406170644754298, 0), -0.3617448639236276, -0.2192285068229158)], []),
	        new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 1), V(1, 0, 0), V(0, 1, 0)), [
		        new PCurveEdge(new BezierCurve(V(0.1, 0.22984366404583298, 1), V(0.1, 0.34200459050920273, 1), V(0.010456949966158688, 0.43292898642298194, 1), V(-0.09999999999999999, 0.43292898642298194, 1), 0, 1), V(0, 0.40575978617006886, 1), V(0.1, 0.22984366404583303, 1), 0.6701998891136561, 1.6497471678008118e-16, null, V(0.2675924213090445, -0.15718180225444303, 0), V(8.863403591374567e-17, -0.3364827793901093, 0)),
		        new StraightEdge(new L3(V(0.16567374856161604, 0.16315679205407474, 1), V(-0.7016738431598636, 0.7124982932086696, 0)), V(0.1, 0.22984366404583303, 1), V(0, 0.33138632523440736, 1), 0.09359583402149652, 0.2361121911222083),
		        new StraightEdge(new L3(V(0, 0, 1), V(0, 1, 0)), V(0, 0.33138632523440736, 1), V(0, 0, 1), 0.33138632523440736, 0),
		        new StraightEdge(new L3(V(1, 0, 1), V(-1, 0, 0)), V(0, 0, 1), V(1, 0, 1), 1, 0),
		        new StraightEdge(new L3(V(1, 1, 1), V(0, -1, 0)), V(1, 0, 1), V(1, 1, 1), 1, 0),
		        new StraightEdge(new L3(V(0, 1, 1), V(1, 0, 0)), V(1, 1, 1), V(0, 1, 1), 1, 0),
		        new StraightEdge(new L3(V(0, 0, 1), V(0, 1, 0)), V(0, 1, 1), V(0, 0.40575978617006886, 1), 1, 0.40575978617006886)], []),
	        new PlaneFace(new PlaneSurface(new P3(V(0, 1, 0), 1), V(0, 0, -1), V(-1, 0, 0)), [
		        new StraightEdge(new L3(V(0, 1, 0), V(1, 0, 0)), V(1, 1, 0), V(0, 1, 0), 1, 0),
		        new StraightEdge(new L3(V(0, 1, 0), V(0, 0, 1)), V(0, 1, 0), V(0, 1, 1), 0, 1),
		        new StraightEdge(new L3(V(0, 1, 1), V(1, 0, 0)), V(0, 1, 1), V(1, 1, 1), 0, 1),
		        new StraightEdge(new L3(V(1, 1, 0), V(0, 0, 1)), V(1, 1, 1), V(1, 1, 0), 1, 0)], []),
	        new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 1), V(0, 0, -1), V(0, 1, 0)), [
		        new StraightEdge(new L3(V(1, 1, 0), V(0, -1, 0)), V(1, 0, 0), V(1, 1, 0), 1, 0),
		        new StraightEdge(new L3(V(1, 1, 0), V(0, 0, 1)), V(1, 1, 0), V(1, 1, 1), 0, 1),
		        new StraightEdge(new L3(V(1, 1, 1), V(0, -1, 0)), V(1, 1, 1), V(1, 0, 1), 0, 1),
		        new StraightEdge(new L3(V(1, 0, 0), V(0, 0, 1)), V(1, 0, 1), V(1, 0, 0), 1, 0)], []),
	        new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0), V(0, 0, -1), V(1, 0, 0)), [
		        new StraightEdge(new L3(V(1, 0, 0), V(-1, 0, 0)), V(0, 0, 0), V(1, 0, 0), 1, 0),
		        new StraightEdge(new L3(V(1, 0, 0), V(0, 0, 1)), V(1, 0, 0), V(1, 0, 1), 0, 1),
		        new StraightEdge(new L3(V(1, 0, 1), V(-1, 0, 0)), V(1, 0, 1), V(0, 0, 1), 0, 1),
		        new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 1), V(0, 0, 0), 1, 0)], []),
	        new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.1, 0.04662674587102256, 2.039074777091188), V(0.1, 0.1554056065717771, 2.058255425163459), V(0.010456949966158688, 0.2435882964734642, 2.073804412624574), V(-0.09999999999999999, 0.2435882964734642, 2.073804412624574), 0, 1), V(0, 0.17364817766693036, -0.9848077530122082), 0, 1, -Infinity, Infinity), [
		        new StraightEdge(new L3(V(3.3306690738754696e-16, 0.21723834785181606, 2.0691582057422955), V(0, 0.17364817766693036, -0.9848077530122082)), V(0, 0.40575978617006886, 1), V(0, 0.582086766878481, 0), 1.0856516944267476, 2.1010783063124836),
		        new PCurveEdge(new BezierCurve(V(-0.09999999999999999, 0.6092559671314469, 0), V(0.010456949966158688, 0.6092559671314469, 0), V(0.1, 0.5183315712176677, 0), V(0.1, 0.4061706447542979, 0), 0, 1), V(0, 0.582086766878481, 0), V(0.1, 0.406170644754298, 0), 0.3298001108864301, 0.9999999999999998, null, V(0.26759242130902194, -0.15718180225447817, 0), V(1.1929530701127963e-16, -0.3364827793901093, 0)),
		        new StraightEdge(new L3(V(0.1, 0.04662674587102256, 2.039074777091188), V(0, 0.17364817766693036, -0.9848077530122081)), V(0.1, 0.406170644754298, 0), V(0.1, 0.22984366404583303, 1), 2.0705307922833858, 1.0551041803976409),
		        new PCurveEdge(new BezierCurve(V(0.1, 0.22984366404583298, 1), V(0.1, 0.34200459050920273, 1), V(0.010456949966158688, 0.43292898642298194, 1), V(-0.09999999999999999, 0.43292898642298194, 1), 0, 1), V(0.1, 0.22984366404583303, 1), V(0, 0.40575978617006886, 1), 1.6497471678008118e-16, 0.6701998891136561, null, V(-8.863403591374567e-17, 0.3364827793901093, 0), V(-0.2675924213090445, 0.15718180225444303, 0))], []),
	        new PlaneFace(new PlaneSurface(new P3(V(0.7071067811865475, 0.6963642403200191, 0.12278780396897289), 0.35355339059327373), V(0.08748610075473388, 0.08615699030406915, -0.9924329474561377), V(-0.7016738431598636, 0.7124982932086696, 0)), [
		        new StraightEdge(new L3(V(0, 0.49240387650610384, 0.08682408883346515), V(0, 0.1736481776669304, -0.9848077530122082)), V(0, 0.5077133059428723, 0), V(0, 0.33138632523440736, 1), 0.08816349035423249, -0.9272631215315127),
		        new StraightEdge(new L3(V(0.16567374856161604, 0.16315679205407474, 1), V(-0.7016738431598636, 0.7124982932086696, 0)), V(0, 0.33138632523440736, 1), V(0.1, 0.22984366404583303, 1), 0.2361121911222083, 0.09359583402149652),
		        new StraightEdge(new L3(V(0.1, 0.04662674587102256, 2.039074777091188), V(0, 0.17364817766693036, -0.9848077530122081)), V(0.1, 0.22984366404583303, 1), V(0.1, 0.406170644754298, 0), 1.0551041803976409, 2.0705307922833858),
		        new StraightEdge(new L3(V(0.2538269089126336, 0.24997070782028524, 0), V(0.7016738431598636, -0.7124982932086696, 0)), V(0.1, 0.406170644754298, 0), V(0, 0.5077133059428723, 0), -0.2192285068229158, -0.3617448639236276)], [])], false)
	    b2EqualAnd(assert, a, b, result)
    },

	'B2 w/ PCS - sphere()'(assert) {
		const a = B2T.sphere().flipped()
		const b = B2T.extrudeEdges([PCurveEdge.forCurveAndTs(BezierCurve.QUARTER_CIRCLE, 0, 1), StraightEdge.throughPoints(V3.Y, V3.X)], P3.XY, V3.Z.negated())
			.scale(0.2, 0.2, 2)
			.translate(0.1, -0.1, 1.2)
		const result = new B2([
			new RotationFace(new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V(0, 0, -1)), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, 1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V(0, 0, -1)), V(0.10000000000000002, 0.1, 0.9899494936611666), V(0.2732435609477567, -6.071532165918825e-18, 0.9619448822051031), -1), V(0.1, 0.1, 0.9899494936611666), V(0.2732435609477568, 0, 0.961944882205103), 0, 67, null, V(0.003311853592333425, 0, -0.00033454773334799895), V(0.0015464506679683768, -0.0026733520656420694, -0.00043927432346962643)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0.2732435609477568, 0, 0.961944882205103), V(0.19999999999999976, 0, 0.9797958971132713), 2.8648293508136553, 2.9402347327994627, null, V(-0.961944882205103, 0, 0.2732435609477568), V(-0.9797958971132713, 0, 0.19999999999999976)),
				new PCurveEdge(new SemiEllipseCurve(V(0.09999999999999994, 0.09999999999999995, 0), V(0.7000000000000001, -0.7, 0), V(0, 0, 0.9899494936611666), 1.4274487578895312, 3.141592653589793), V(0.19999999999999976, 0, 0.9797958971132713), V(0.1, 0.1, 0.9899494936611666), 1.4274487578895314, 1.5707963267948966, null, V(-0.6928203230275509, 0.6928203230275509, 0.14142135623730936), V(-0.7000000000000001, 0.7, 6.061692393648453e-17))], []),
			new RotationFace(new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V(0, 0, -1)), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, 1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V(0, 0, -1)), V(0.27324356094775665, -3.2161642623246884e-17, 0.9619448822051031), V(0.30000000000000004, -0.1, 0.9486832980505138), -1), V(0.2732435609477568, 0, 0.961944882205103), V(0.30000000000000004, -0.1, 0.9486832980505138), 0, 33, null, V(0.0015464443998091284, -0.002673341229863844, -0.00043927254297868013), V(0, -0.0033116913504264997, -0.00034908291916088577)),
				new PCurveEdge(new SemiEllipseCurve(V(0.09999999999999994, 0.09999999999999995, 0), V(0, 0, -0.9899494936611666), V(0.7000000000000001, -0.7, 0), 0.14334756890536535, 2.9982450846844277), V(0.30000000000000004, -0.1, 0.9486832980505138), V(0.19999999999999976, 0, 0.9797958971132713), 2.851840952153746, 2.9982450846844277, null, V(-0.670820393249937, 0.6708203932499369, 0.282842712474619), V(-0.6928203230275509, 0.6928203230275509, 0.1414213562373092)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0.19999999999999976, 0, 0.9797958971132713), V(0.2732435609477568, 0, 0.961944882205103), 2.9402347327994627, 2.8648293508136553, null, V(0.9797958971132713, 0, -0.19999999999999976), V(0.961944882205103, 0, -0.2732435609477568))], []),
			new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, 1), 0, 1, -Infinity, Infinity), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, 1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V(0, 0, -1)), V(0.10000000000000002, 0.1, 0.9899494936611666), V(0.2732435609477567, -6.071532165918825e-18, 0.9619448822051031), -1), V(0.2732435609477568, 0, 0.961944882205103), V(0.1, 0.1, 0.9899494936611666), 67, 0, null, V(-0.0015464506679683768, 0.0026733520656420694, 0.00043927432346962643), V(-0.003311853592333425, 0, 0.00033454773334799895)),
				new StraightEdge(new L3(V(0.1, 0.1, 1.2), V(0, 0, -1)), V(0.1, 0.1, 0.9899494936611666), V(0.1, 0.1, 1.2), 0.21005050633883338, 0),
				new PCurveEdge(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0.10000000000000002, 0.1, 1.2), V(0.30000000000000004, -0.1, 1.2), 1, 0, null, V(0.33137084989847604, 0, 0), V(0, -0.33137084989847604, 0)),
				new StraightEdge(new L3(V(0.30000000000000004, -0.1, 1.2), V(0, 0, -1)), V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, -0.1, 0.9486832980505138), 0, 0.2513167019494862),
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, 1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V(0, 0, -1)), V(0.27324356094775665, -3.2161642623246884e-17, 0.9619448822051031), V(0.30000000000000004, -0.1, 0.9486832980505138), -1), V(0.30000000000000004, -0.1, 0.9486832980505138), V(0.2732435609477568, 0, 0.961944882205103), 33, 0, null, V(0, 0.0033116913504264997, 0.00034908291916088577), V(-0.0015464443998091284, 0.002673341229863844, 0.00043927254297868013))], []),
			new PlaneFace(new PlaneSurface(new P3(V(-0.7071067811865475, -0.7071067811865476, 0), -0.14142135623730942), V(0, 0, -1), V(0.7071067811865476, -0.7071067811865475, 0)), [
				new PCurveEdge(new SemiEllipseCurve(V(0.09999999999999994, 0.09999999999999995, 0), V(0.7000000000000001, -0.7, 0), V(0, 0, 0.9899494936611666), 1.4274487578895312, 3.141592653589793), V(0.1, 0.1, 0.9899494936611666), V(0.19999999999999976, 0, 0.9797958971132713), 1.5707963267948966, 1.4274487578895314, null, V(0.7000000000000001, -0.7, -6.061692393648453e-17), V(0.6928203230275509, -0.6928203230275509, -0.14142135623730936)),
				new PCurveEdge(new SemiEllipseCurve(V(0.09999999999999994, 0.09999999999999995, 0), V(0, 0, -0.9899494936611666), V(0.7000000000000001, -0.7, 0), 0.14334756890536535, 2.9982450846844277), V(0.19999999999999976, 0, 0.9797958971132713), V(0.30000000000000004, -0.1, 0.9486832980505138), 2.9982450846844277, 2.851840952153746, null, V(0.6928203230275509, -0.6928203230275509, -0.1414213562373092), V(0.670820393249937, -0.6708203932499369, -0.282842712474619)),
				new StraightEdge(new L3(V(0.30000000000000004, -0.1, 1.2), V(0, 0, -1)), V(0.30000000000000004, -0.1, 0.9486832980505138), V(0.30000000000000004, -0.1, 1.2), 0.2513167019494862, 0),
				new StraightEdge(new L3(V(0.1, 0.1, 1.2), V(0.7071067811865475, -0.7071067811865475, 0)), V(0.30000000000000004, -0.1, 1.2), V(0.1, 0.1, 1.2), 0.28284271247461906, 0),
				new StraightEdge(new L3(V(0.1, 0.1, 1.2), V(0, 0, -1)), V(0.1, 0.1, 1.2), V(0.1, 0.1, 0.9899494936611666), 0, 0.21005050633883338)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 1.2), V(1, 0, 0), V(0, 1, 0)), [
				new PCurveEdge(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0.30000000000000004, -0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1, null, V(0, 0.33137084989847604, 0), V(-0.33137084989847604, 0, 0)),
				new StraightEdge(new L3(V(0.1, 0.1, 1.2), V(0.7071067811865475, -0.7071067811865475, 0)), V(0.1, 0.1, 1.2), V(0.30000000000000004, -0.1, 1.2), 0, 0.28284271247461906)], [])], false)
		b2EqualAnd(assert, a, b, result)
	},

	'sphere() - B2 w/ PCS - sphere(0.9)'(assert) {
		const a = B2T.sphere(0.9).flipped()
		const b = B2T.extrudeEdges([PCurveEdge.forCurveAndTs(BezierCurve.QUARTER_CIRCLE, 0, 1), StraightEdge.throughPoints(V3.Y, V3.X)], P3.XY, V3.Z.negated())
			.scale(0.2, 0.2, 2)
			.translate(0.1, -0.1, 1.2).flipped()
		const c = B2T.sphere()
		const d = a.and(b)
		const result = new B2([
			new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.2732435609477567, -6.071532165918825e-18, 0.9619448822051031), V(0.10000000000000002, 0.1, 0.9899494936611666), 1), V(0.2732435609477568, 0, 0.961944882205103), V(0.1, 0.1, 0.9899494936611666), 0, 67, null, V(-0.0021067547553675317, 0.0036419507545378837, 0.0005984305151461936), V(-0.003311899411325918, 0, 0.0003345523617651845)),
				new StraightEdge(new L3(V(0.1, 0.1, 1.2), V(0, 0, -1)), V(0.1, 0.1, 0.9899494936611666), V(0.1, 0.1, 0.8888194417315591), 0.21005050633883338, 0.3111805582684408),
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(0.27324356094775665, 4.6837533851373795e-18, 0.8575184874978438), V(0.1, 0.09999999999999999, 0.8888194417315589), -1), V(0.1, 0.1, 0.8888194417315591), V(0.273243560947757, 0, 0.8575184874978437), 67, 0, null, V(0.003311454530207254, 0, -0.0003725677426402853), V(0.0015460590308877795, -0.002672675041912773, -0.000492643227165731)),
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(0.3, -0.1, 0.8426149773176359), V(0.2732435609477566, -3.970768974374599e-17, 0.8575184874978438), -1), V(0.273243560947757, 0, 0.8575184874978437), V(0.30000000000000004, -0.1, 0.8426149773176361), 33, 0, null, V(0.0015460577093711644, -0.0026726727574046118, -0.0004926428060717128), V(0, -0.0033111490223668035, -0.0003929610927291431)),
				new StraightEdge(new L3(V(0.30000000000000004, -0.1, 1.2), V(0, 0, -1)), V(0.30000000000000004, -0.1, 0.8426149773176361), V(0.30000000000000004, -0.1, 0.9486832980505138), 0.35738502268236383, 0.2513167019494862),
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.30000000000000004, -0.1, 0.9486832980505138), V(0.27324356094775665, -3.2161642623246884e-17, 0.9619448822051031), 1), V(0.30000000000000004, -0.1, 0.9486832980505138), V(0.2732435609477568, 0, 0.961944882205103), 0, 33, null, V(0, 0.0033116888624497654, 0.0003490826569051109), V(-0.0015464438590823478, 0.002673340295108481, 0.0004392723893834859))], []),
			new PlaneFace(new PlaneSurface(new P3(V(0.7071067811865475, 0.7071067811865476, 0), 0.14142135623730942), V(0, 0, -1), V(-0.7071067811865476, 0.7071067811865475, 0)), [
				new PCurveEdge(new SemiEllipseCurve(V(0.09999999999999994, 0.09999999999999995, 0), V(0.7000000000000001, -0.7, 0), V(0, 0, 0.9899494936611666), 1.4274487578895312, 3.141592653589793), V(0.1, 0.1, 0.9899494936611666), V(0.19999999999999976, 0, 0.9797958971132713), 1.5707963267948966, 1.4274487578895314, null, V(0.7000000000000001, -0.7, -6.061692393648453e-17), V(0.6928203230275509, -0.6928203230275509, -0.14142135623730936)),
				new PCurveEdge(new SemiEllipseCurve(V(0.09999999999999994, 0.09999999999999995, 0), V(0, 0, -0.9899494936611666), V(0.7000000000000001, -0.7, 0), 0.14334756890536535, 2.9982450846844277), V(0.19999999999999976, 0, 0.9797958971132713), V(0.30000000000000004, -0.1, 0.9486832980505138), 2.9982450846844277, 2.851840952153746, null, V(0.6928203230275509, -0.6928203230275509, -0.1414213562373092), V(0.670820393249937, -0.6708203932499369, -0.282842712474619)),
				new StraightEdge(new L3(V(0.30000000000000004, -0.1, 1.2), V(0, 0, -1)), V(0.30000000000000004, -0.1, 0.9486832980505138), V(0.30000000000000004, -0.1, 0.8426149773176361), 0.2513167019494862, 0.35738502268236383),
				new PCurveEdge(new SemiEllipseCurve(V(0.0999999999999999, 0.09999999999999994, 0), V(0, 0, 0.8888194417315588), V(0.6284902544988268, -0.6284902544988266, 0), 0.15979057883583087, 2.981802074753962), V(0.30000000000000004, -0.1, 0.8426149773176361), V(0.20000000000000012, 0, 0.8774964387392122), 0.32385436390177574, 0.15979057883583106, null, V(-0.5958187643906493, 0.5958187643906491, 0.282842712474619), V(-0.6204836822995429, 0.6204836822995426, 0.14142135623730964)),
				new PCurveEdge(new SemiEllipseCurve(V(0.09999999999999992, 0.09999999999999994, 0), V(-0.6284902544988268, 0.6284902544988267, 0), V(0, 0, 0.8888194417315588), 0, 1.7305869056307275), V(0.20000000000000012, 0, 0.8774964387392122), V(0.1, 0.1, 0.8888194417315591), 1.7305869056307275, 1.5707963267948968, null, V(-0.6204836822995429, 0.6204836822995428, 0.14142135623730967), V(-0.6284902544988268, 0.6284902544988267, 1.4293306757214634e-16)),
				new StraightEdge(new L3(V(0.1, 0.1, 1.2), V(0, 0, -1)), V(0.1, 0.1, 0.8888194417315591), V(0.1, 0.1, 0.9899494936611666), 0.3111805582684408, 0.21005050633883338)], []),
			new RotationFace(new SemiEllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(0.27324356094775665, 4.6837533851373795e-18, 0.8575184874978438), V(0.1, 0.09999999999999999, 0.8888194417315589), -1), V(0.273243560947757, 0, 0.8575184874978437), V(0.1, 0.1, 0.8888194417315591), 0, 67, null, V(-0.0015460590308877795, 0.002672675041912773, 0.000492643227165731), V(-0.003311454530207254, 0, 0.0003725677426402853)),
				new PCurveEdge(new SemiEllipseCurve(V(0.09999999999999992, 0.09999999999999994, 0), V(-0.6284902544988268, 0.6284902544988267, 0), V(0, 0, 0.8888194417315588), 0, 1.7305869056307275), V(0.1, 0.1, 0.8888194417315591), V(0.20000000000000012, 0, 0.8774964387392122), 1.5707963267948968, 1.7305869056307275, null, V(0.6284902544988268, -0.6284902544988267, -1.4293306757214634e-16), V(0.6204836822995429, -0.6204836822995428, -0.14142135623730967)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(0.9, 0, 0)), V(0.20000000000000012, 0, 0.8774964387392122), V(0, 0, 0.9), 2.9174995612884222, 3.141592653589793, null, V(-0.8774964387392122, 0, 0.20000000000000012), V(-0.9, 0, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(-0.9, 1.1021821192326179e-16, 0), 0, 3.141592653589793), V(0, 0, 0.9), V(0, 0, -0.9), 3.141592653589793, 0, null, V(-0.9, 1.1021821192326179e-16, 0), V(0.9, -1.1021821192326179e-16, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(0.9, 0, 0)), V(0, 0, -0.9), V(0.273243560947757, 0, 0.8575184874978437), 0, 2.833119770489437, null, V(0.9, 0, 0), V(-0.8575184874978437, 0, 0.273243560947757))], []),
			new RotationFace(new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(0.3, -0.1, 0.8426149773176359), V(0.2732435609477566, -3.970768974374599e-17, 0.8575184874978438), -1), V(0.30000000000000004, -0.1, 0.8426149773176361), V(0.273243560947757, 0, 0.8575184874978437), 0, 33, null, V(0, 0.0033111490223668035, 0.0003929610927291431), V(-0.0015460577093711644, 0.0026726727574046118, 0.0004926428060717128)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(0.9, 0, 0)), V(0.273243560947757, 0, 0.8575184874978437), V(0, 0, -0.9), 2.833119770489437, 0, null, V(0.8575184874978437, 0, -0.273243560947757), V(-0.9, 0, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(-0.9, 1.1021821192326179e-16, 0), 0, 3.141592653589793), V(0, 0, -0.9), V(0, 0, 0.9), 0, 3.141592653589793, null, V(-0.9, 1.1021821192326179e-16, 0), V(0.9, -1.1021821192326179e-16, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(0.9, 0, 0)), V(0, 0, 0.9), V(0.20000000000000012, 0, 0.8774964387392122), 3.141592653589793, 2.9174995612884222, null, V(0.9, 0, 0), V(0.8774964387392122, 0, -0.20000000000000012)),
				new PCurveEdge(new SemiEllipseCurve(V(0.0999999999999999, 0.09999999999999994, 0), V(0, 0, 0.8888194417315588), V(0.6284902544988268, -0.6284902544988266, 0), 0.15979057883583087, 2.981802074753962), V(0.20000000000000012, 0, 0.8774964387392122), V(0.30000000000000004, -0.1, 0.8426149773176361), 0.15979057883583106, 0.32385436390177574, null, V(0.6204836822995429, -0.6204836822995426, -0.14142135623730964), V(0.5958187643906493, -0.5958187643906491, -0.282842712474619))], []),
			new RotationFace(new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.2732435609477567, -6.071532165918825e-18, 0.9619448822051031), V(0.10000000000000002, 0.1, 0.9899494936611666), 1), V(0.1, 0.1, 0.9899494936611666), V(0.2732435609477568, 0, 0.961944882205103), 67, 0, null, V(0.003311899411325918, 0, -0.0003345523617651845), V(0.0021067547553675317, -0.0036419507545378837, -0.0005984305151461936)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0.2732435609477568, 0, 0.961944882205103), V(0, 0, -1), 2.8648293508136553, 0, null, V(0.961944882205103, 0, -0.2732435609477568), V(-1, 0, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(0, 0, -1), V(0, 0, 1), 0, 3.141592653589793, null, V(-1, 1.2246467991473532e-16, 0), V(1, -1.2246467991473532e-16, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0, 0, 1), V(0.19999999999999976, 0, 0.9797958971132713), 3.141592653589793, 2.9402347327994627, null, V(1, 0, 0), V(0.9797958971132713, 0, -0.19999999999999976)),
				new PCurveEdge(new SemiEllipseCurve(V(0.09999999999999994, 0.09999999999999995, 0), V(0.7000000000000001, -0.7, 0), V(0, 0, 0.9899494936611666), 1.4274487578895312, 3.141592653589793), V(0.19999999999999976, 0, 0.9797958971132713), V(0.1, 0.1, 0.9899494936611666), 1.4274487578895314, 1.5707963267948966, null, V(-0.6928203230275509, 0.6928203230275509, 0.14142135623730936), V(-0.7000000000000001, 0.7, 6.061692393648453e-17))], []),
			new RotationFace(new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004, 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.30000000000000004, -0.1, 0.9486832980505138), V(0.27324356094775665, -3.2161642623246884e-17, 0.9619448822051031), 1), V(0.2732435609477568, 0, 0.961944882205103), V(0.30000000000000004, -0.1, 0.9486832980505138), 33, 0, null, V(0.0015464438590823478, -0.002673340295108481, -0.0004392723893834859), V(0, -0.0033116888624497654, -0.0003490826569051109)),
				new PCurveEdge(new SemiEllipseCurve(V(0.09999999999999994, 0.09999999999999995, 0), V(0, 0, -0.9899494936611666), V(0.7000000000000001, -0.7, 0), 0.14334756890536535, 2.9982450846844277), V(0.30000000000000004, -0.1, 0.9486832980505138), V(0.19999999999999976, 0, 0.9797958971132713), 2.851840952153746, 2.9982450846844277, null, V(-0.670820393249937, 0.6708203932499369, 0.282842712474619), V(-0.6928203230275509, 0.6928203230275509, 0.1414213562373092)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0.19999999999999976, 0, 0.9797958971132713), V(0, 0, 1), 2.9402347327994627, 3.141592653589793, null, V(-0.9797958971132713, 0, 0.19999999999999976), V(-1, 0, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(0, 0, 1), V(0, 0, -1), 3.141592653589793, 0, null, V(-1, 1.2246467991473532e-16, 0), V(1, -1.2246467991473532e-16, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0, 0, -1), V(0.2732435609477568, 0, 0.961944882205103), 0, 2.8648293508136553, null, V(1, 0, 0), V(-0.961944882205103, 0, 0.2732435609477568))], [])], false)
		b2Equal(assert, d, c, d.and(c), result)
	},

	'sphere() - cylinder(0.05, 4).scale(10,1,1).translate(0.5,0,-2)'(assert) {
		const a = B2T.sphere()
		const b = B2T.cylinder(0.05, 4).scale(10,1,1).translate(0.5,0,-2).flipped()
		const result = B2.EMPTY
		b2EqualAnd(assert, a, b, result)
	},

	'box() - cylinder(0.2,2).translate(0.5,0.2)'(assert) {
		const a = B2T.box()
		const b = B2T.cylinder(0.2,2).translate(0.5,0.2).flipped()
		const result = B2.EMPTY
		b2EqualAnd(assert, a, b, result)
	},

	//'cylinder() - triangle(1.5)'(assert) {
		// file:///C:/Users/aval/Desktop/cs/viewer.html?a=B2T.cylinder()&b=B2T.extrudeEdges(Edge.ngon(3,1.5))&c=a.and(b).translate(3)
		//const a = B2T.cylinder()
		//const b = B2T.extrudeEdges(Edge.round(Edge.ngon(3,1.6),0.5))
		//const result = B2.EMPTY
		//b2EqualAnd(assert, a, b, result)
	//},
	'box - snug sphere'(assert) {
		// file:///C:/Users/aval/Desktop/cs/viewer.html?a=B2T.cylinder()&b=B2T.extrudeEdges(Edge.ngon(3,1.5))&c=a.and(b).translate(3)
		const a = B2T.box(2,2,3)
		const b = B2T.sphere().translate(1,1).flipped()
		const result = new B2([
			new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), 0),V(0, -1, 0),V3.Z), [
				new StraightEdge(new L3(V3.O, V3.Y), V3.Y, V3.O, 1, 0),
				new StraightEdge(new L3(V3.O, V3.Z), V3.O, V(0, 0, 3), 0, 3),
				new StraightEdge(new L3(V(0, 0, 3), V3.Y), V(0, 0, 3), V(0, 2, 3), 0, 2),
				new StraightEdge(new L3(V(0, 2, 0), V3.Z), V(0, 2, 3), V(0, 2, 0), 3, 0),
				new StraightEdge(new L3(V3.O, V3.Y), V(0, 2, 0), V3.Y, 2, 1)], []),
			new PlaneFace(new PlaneSurface(new P3(V3.Y, 2),V(-1, 0, 0),V3.Z), [
				new StraightEdge(new L3(V(0, 2, 0), V3.X), V(1, 2, 0), V(0, 2, 0), 1, 0),
				new StraightEdge(new L3(V(0, 2, 0), V3.Z), V(0, 2, 0), V(0, 2, 3), 0, 3),
				new StraightEdge(new L3(V(0, 2, 3), V3.X), V(0, 2, 3), V(2, 2, 3), 0, 2),
				new StraightEdge(new L3(V(2, 2, 0), V3.Z), V(2, 2, 3), V(2, 2, 0), 3, 0),
				new StraightEdge(new L3(V(0, 2, 0), V3.X), V(2, 2, 0), V(1, 2, 0), 2, 1)], []),
			new PlaneFace(new PlaneSurface(new P3(V3.X, 2),V3.Y,V3.Z), [
				new StraightEdge(new L3(V(2, 2, 0), V(0, -1, 0)), V(2, 1, 0), V(2, 2, 0), 1, 0),
				new StraightEdge(new L3(V(2, 2, 0), V3.Z), V(2, 2, 0), V(2, 2, 3), 0, 3),
				new StraightEdge(new L3(V(2, 2, 3), V(0, -1, 0)), V(2, 2, 3), V(2, 0, 3), 0, 2),
				new StraightEdge(new L3(V(2, 0, 0), V3.Z), V(2, 0, 3), V(2, 0, 0), 3, 0),
				new StraightEdge(new L3(V(2, 2, 0), V(0, -1, 0)), V(2, 0, 0), V(2, 1, 0), 2, 1)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0),V3.X,V3.Z), [
				new StraightEdge(new L3(V(2, 0, 0), V(-1, 0, 0)), V3.X, V(2, 0, 0), 0.9999999999999999, 0),
				new StraightEdge(new L3(V(2, 0, 0), V3.Z), V(2, 0, 0), V(2, 0, 3), 0, 3),
				new StraightEdge(new L3(V(2, 0, 3), V(-1, 0, 0)), V(2, 0, 3), V(0, 0, 3), 0, 2),
				new StraightEdge(new L3(V3.O, V3.Z), V(0, 0, 3), V3.O, 3, 0),
				new StraightEdge(new L3(V(2, 0, 0), V(-1, 0, 0)), V3.O, V3.X, 2, 0.9999999999999999)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 0),V3.Y,V3.X), [
				new PCurveEdge(new SemiEllipseCurve(V(1, 1, 0),V(-1, 0, 0),V3.Y,0,3.141592653589793),V(1, 2, 0),V(0, 1.0000000000000002, -6.123233995736766e-17),1.5707963267948966,2.220446049250313e-16,null,V(-1, -6.123233995736766e-17, 0),V(-2.220446049250313e-16, -1, 0)),
				new StraightEdge(new L3(V3.O, V3.Y), V3.Y, V(0, 2, 0), 1, 2),
				new StraightEdge(new L3(V(0, 2, 0), V3.X), V(0, 2, 0), V(1, 2, 0), 0, 1)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 0),V3.Y,V3.X), [
				new PCurveEdge(new SemiEllipseCurve(V(1, 1, 0),V(-1, 0, 0),V3.Y,0,3.141592653589793),V(2, 1, -6.123233995736766e-17),V(1, 2, 0),3.141592653589793,1.5707963267948966,null,V(-1.2246467991473532e-16, 1, 0),V(-1, -6.123233995736766e-17, 0)),
				new StraightEdge(new L3(V(0, 2, 0), V3.X), V(1, 2, 0), V(2, 2, 0), 1, 2),
				new StraightEdge(new L3(V(2, 2, 0), V(0, -1, 0)), V(2, 2, 0), V(2, 1, 0), 0, 1)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 0),V3.Y,V3.X), [
				new PCurveEdge(new SemiEllipseCurve(V(1, 1, 0),V(1, -1.2246467991473532e-16, 0),V(-1.2246467991473532e-16, -1, 0),0,3.141592653589793),V3.X,V(2, 1, -6.123233995736766e-17),1.5707963267948966,0,null,V(1, -6.123233995736766e-17, 0),V(1.2246467991473532e-16, 1, 0)),
				new StraightEdge(new L3(V(2, 2, 0), V(0, -1, 0)), V(2, 1, 0), V(2, 0, 0), 1, 2),
				new StraightEdge(new L3(V(2, 0, 0), V(-1, 0, 0)), V(2, 0, 0), V3.X, 0, 0.9999999999999999)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 0),V3.Y,V3.X), [
				new PCurveEdge(new SemiEllipseCurve(V(1, 1, 0),V(1, -1.2246467991473532e-16, 0),V(-1.2246467991473532e-16, -1, 0),0,3.141592653589793),V(0, 1.0000000000000002, -6.123233995736766e-17),V3.X,3.141592653589793,1.5707963267948966,null,V(0, -1, 0),V(1, -6.123233995736766e-17, 0)),
				new StraightEdge(new L3(V(2, 0, 0), V(-1, 0, 0)), V3.X, V3.O, 0.9999999999999999, 2),
				new StraightEdge(new L3(V3.O, V3.Y), V3.O, V3.Y, 0, 1)], []),
			new PlaneFace(new PlaneSurface(new P3(V3.Z, 3),V3.Y,V(-1, 0, 0)), [
				new StraightEdge(new L3(V(2, 0, 3), V(-1, 0, 0)), V(0, 0, 3), V(2, 0, 3), 2, 0),
				new StraightEdge(new L3(V(2, 2, 3), V(0, -1, 0)), V(2, 0, 3), V(2, 2, 3), 2, 0),
				new StraightEdge(new L3(V(0, 2, 3), V3.X), V(2, 2, 3), V(0, 2, 3), 2, 0),
				new StraightEdge(new L3(V(0, 0, 3), V3.Y), V(0, 2, 3), V(0, 0, 3), 2, 0)], []),
			new RotationFace(new SemiEllipsoidSurface(V(1, 1, 0), V3.X, V3.Y, V(0, 0, -1)), [
				new PCurveEdge(new SemiEllipseCurve(V(1, 1, 0),V(-1, 0, 0),V3.Y,0,3.141592653589793),V(0, 1.0000000000000002, -6.123233995736766e-17),V(1, 2, 0),2.220446049250313e-16,1.5707963267948966,null,V(2.220446049250313e-16, 1, 0),V(1, 6.123233995736766e-17, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(1, 1, 0),V(-1, 0, 0),V3.Y,0,3.141592653589793),V(1, 2, 0),V(2, 1, -6.123233995736766e-17),1.5707963267948966,3.141592653589793,null,V(1, 6.123233995736766e-17, 0),V(1.2246467991473532e-16, -1, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(1, 1, 0),V(0, 0, -1),V3.X,0,3.141592653589793),V(2, 1, -6.123233995736766e-17),V3.XYZ,1.5707963267948966,3.141592653589793,null,V(6.123233995736766e-17, 0, 1),V(-1, 0, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(1, 1, 0),V(0, 0, -1),V(-1, 1.2246467991473532e-16, 0),0,3.141592653589793),V3.XYZ,V(0, 1.0000000000000002, -6.123233995736766e-17),3.141592653589793,1.5707963267948966,null,V(-1, 1.2246467991473532e-16, 0),V(6.123233995736766e-17, -7.498798913309288e-33, -1))], []),
			new RotationFace(new SemiEllipsoidSurface(V(1, 1, 0), V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V(0, 0, -1)), [
				new PCurveEdge(new SemiEllipseCurve(V(1, 1, 0),V(1, -1.2246467991473532e-16, 0),V(-1.2246467991473532e-16, -1, 0),0,3.141592653589793),V(2, 1, -6.123233995736766e-17),V3.X,0,1.5707963267948966,null,V(-1.2246467991473532e-16, -1, 0),V(-1, 6.123233995736766e-17, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(1, 1, 0),V(1, -1.2246467991473532e-16, 0),V(-1.2246467991473532e-16, -1, 0),0,3.141592653589793),V3.X,V(0, 1.0000000000000002, -6.123233995736766e-17),1.5707963267948966,3.141592653589793,null,V(-1, 6.123233995736766e-17, 0),V3.Y),
				new PCurveEdge(new SemiEllipseCurve(V(1, 1, 0),V(0, 0, -1),V(-1, 1.2246467991473532e-16, 0),0,3.141592653589793),V(0, 1.0000000000000002, -6.123233995736766e-17),V3.XYZ,1.5707963267948966,3.141592653589793,null,V(-6.123233995736766e-17, 7.498798913309288e-33, 1),V(1, -1.2246467991473532e-16, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(1, 1, 0),V(0, 0, -1),V3.X,0,3.141592653589793),V3.XYZ,V(2, 1, -6.123233995736766e-17),3.141592653589793,1.5707963267948966,null,V3.X,V(-6.123233995736766e-17, 0, -1))], [])], false)
		b2EqualAnd(assert, a, b, result)
	},
	'pcs - triangle'(assert) {
		// file:///C:/Users/aval/Desktop/cs/viewer.html?a=B2T.cylinder()&b=B2T.extrudeEdges(Edge.ngon(3,1.5))&c=a.and(b).translate(3)
		const e1 = Edge.forCurveAndTs(BezierCurve.EX3D,0.5).projXY()
		const edges = [e1, StraightEdge.throughPoints(e1.b, e1.a)]
		const a = B2T.extrudeEdges(edges, P3.XY.flipped(), V3.Z)
		const p = e1.curve.at(0.95).plus(V(-NLA_PRECISION*2,0,0))
		const b = B2T.extrudeVertices([p, V(0,-0.2), V(0,0.2)], P3.XY.flipped(), V3.Z).flipped()
		const result = B2.EMPTY
		b2EqualAnd(assert, a, b, result)
	},
	'box() - cylinder(0.2,3).translate(0.5,0.2,-1).rotateX(10*DEG)'(assert) {
		const a = B2T.box()
		const b = B2T.cylinder(0.2,1).translate(0.5,0.2,-0.5).rotateX(-10*DEG).flipped()
		const result = new B2([
			new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0),V3.X,V3.Z), [
				new StraightEdge(new L3(V3.X, V(-1, 0, 0)), V(0.5, 0, 0), V3.X, 0.5, 0),
				new StraightEdge(new L3(V3.X, V3.Z), V3.X, V(1, 0, 1), 0, 1),
				new StraightEdge(new L3(V(1, 0, 1), V(-1, 0, 0)), V(1, 0, 1), V3.Z, 0, 1),
				new StraightEdge(new L3(V3.O, V3.Z), V3.Z, V3.O, 1, 0),
				new StraightEdge(new L3(V3.X, V(-1, 0, 0)), V3.O, V(0.5, 0, 0), 1, 0.5)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 0),V3.Y,V3.X), [
				new PCurveEdge(new SemiEllipseCurve(V(0.5, 0.203085322377149, 0),V(-0.2, 0, 0),V(0, 0.203085322377149, 0),0,3.141592653589793),V(0.7, 0.203085322377149, 0),V(0.3, 0.20308532237714905, 0),3.141592653589793,2.220446049250313e-16,null,V(-2.4492935982947065e-17, 0.203085322377149, 0),V(-4.4408920985006264e-17, -0.203085322377149, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0.5, 0.203085322377149, 0),V(0.2, -2.4870779000298388e-17, 0),V(-2.4492935982947065e-17, -0.203085322377149, 0),0,3.141592653589793),V(0.3, 0.20308532237714905, 0),V(0.5, 0, 0),3.141592653589793,1.5707963267948966,null,V(0, -0.203085322377149, 0),V(0.2, -1.2435389500149195e-17, 0)),
				new StraightEdge(new L3(V3.X, V(-1, 0, 0)), V(0.5, 0, 0), V3.O, 0.5, 1),
				new StraightEdge(new L3(V3.O, V3.Y), V3.O, V3.Y, 0, 1),
				new StraightEdge(new L3(V3.Y, V3.X), V3.Y, V(1, 1, 0), 0, 1),
				new StraightEdge(new L3(V(1, 1, 0), V(0, -1, 0)), V(1, 1, 0), V3.X, 0, 1),
				new StraightEdge(new L3(V3.X, V(-1, 0, 0)), V3.X, V(0.5, 0, 0), 0, 0.5),
				new PCurveEdge(new SemiEllipseCurve(V(0.5, 0.203085322377149, 0),V(0.2, -2.4870779000298388e-17, 0),V(-2.4492935982947065e-17, -0.203085322377149, 0),0,3.141592653589793),V(0.5, 0, 0),V(0.7, 0.203085322377149, 0),1.5707963267948966,0,null,V(0.2, -1.2435389500149195e-17, 0),V(2.4492935982947065e-17, 0.203085322377149, 0))], []),
			new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), 0),V(0, -1, 0),V3.Z), [
				new StraightEdge(new L3(V3.O, V3.Y), V3.Y, V3.O, 1, 0),
				new StraightEdge(new L3(V3.O, V3.Z), V3.O, V3.Z, 0, 1),
				new StraightEdge(new L3(V3.Z, V3.Y), V3.Z, V(0, 1, 1), 0, 1),
				new StraightEdge(new L3(V3.Y, V3.Z), V(0, 1, 1), V3.Y, 1, 0)], []),
			new PlaneFace(new PlaneSurface(new P3(V3.Y, 1),V(-1, 0, 0),V3.Z), [
				new StraightEdge(new L3(V3.Y, V3.X), V(1, 1, 0), V3.Y, 1, 0),
				new StraightEdge(new L3(V3.Y, V3.Z), V3.Y, V(0, 1, 1), 0, 1),
				new StraightEdge(new L3(V(0, 1, 1), V3.X), V(0, 1, 1), V3.XYZ, 0, 1),
				new StraightEdge(new L3(V(1, 1, 0), V3.Z), V3.XYZ, V(1, 1, 0), 1, 0)], []),
			new PlaneFace(new PlaneSurface(new P3(V3.X, 1),V3.Y,V3.Z), [
				new StraightEdge(new L3(V(1, 1, 0), V(0, -1, 0)), V3.X, V(1, 1, 0), 1, 0),
				new StraightEdge(new L3(V(1, 1, 0), V3.Z), V(1, 1, 0), V3.XYZ, 0, 1),
				new StraightEdge(new L3(V3.XYZ, V(0, -1, 0)), V3.XYZ, V(1, 0, 1), 0, 1),
				new StraightEdge(new L3(V3.X, V3.Z), V(1, 0, 1), V3.X, 1, 0)], []),
			new PlaneFace(new PlaneSurface(new P3(V3.Z, 1),V3.Y,V(-1, 0, 0)), [
				new StraightEdge(new L3(V(1, 0, 1), V(-1, 0, 0)), V3.Z, V(1, 0, 1), 1, 0),
				new StraightEdge(new L3(V3.XYZ, V(0, -1, 0)), V(1, 0, 1), V3.XYZ, 1, 0),
				new StraightEdge(new L3(V(0, 1, 1), V3.X), V3.XYZ, V(0, 1, 1), 1, 0),
				new StraightEdge(new L3(V3.Z, V3.Y), V(0, 1, 1), V3.Z, 1, 0)], []),
			new RotationFace(new SemiCylinderSurface(new SemiEllipseCurve(V(0.5, 0.11013746176897644, -0.5271335120394901),V(0.2, 0, 0),V(0, 0.1969615506024416, -0.034729635533386066),0,3.141592653589793),V(0, -0.17364817766693036, -0.9848077530122082),-Infinity,Infinity), [
				new PCurveEdge(new SemiEllipseCurve(V(0.5, 0.203085322377149, 0),V(-0.2, 0, 0),V(0, 0.203085322377149, 0),0,3.141592653589793),V(0.3, 0.20308532237714905, 0),V(0.7, 0.203085322377149, 0),2.220446049250313e-16,3.141592653589793,null,V(4.4408920985006264e-17, 0.203085322377149, 0),V(2.4492935982947065e-17, -0.203085322377149, 0)),
				new StraightEdge(new L3(V(0.7, 0.11013746176897644, -0.5271335120394901), V(0, 0.17364817766693036, 0.9848077530122081)), V(0.7, 0.203085322377149, 0), V(0.7, 0.2837856394359068, 0.4576742409727179), 0.535265396141693, 0.9999999999999999),
				new PCurveEdge(new SemiEllipseCurve(V(0.5, 0.2837856394359068, 0.4576742409727179),V(0.2, 0, 0),V(0, 0.1969615506024416, -0.034729635533386066),0,3.141592653589793),V(0.7, 0.2837856394359068, 0.4576742409727179),V(0.3, 0.2837856394359068, 0.4576742409727179),0,3.141592653589793,null,V(0, 0.1969615506024416, -0.034729635533386066),V(-2.4492935982947065e-17, -0.1969615506024416, 0.034729635533386066)),
				new StraightEdge(new L3(V(0.3, 0.11013746176897647, -0.5271335120394901), V(0, 0.17364817766693036, 0.9848077530122081)), V(0.3, 0.2837856394359068, 0.4576742409727179), V(0.3, 0.20308532237714905, 0), 0.9999999999999999, 0.535265396141693)], []),
			new RotationFace(new SemiCylinderSurface(new SemiEllipseCurve(V(0.5, 0.11013746176897644, -0.5271335120394901),V(-0.2, 2.4120833250037956e-17, -4.2531536991515426e-18),V(-2.4492935982947065e-17, -0.1969615506024416, 0.034729635533386066),0,3.141592653589793),V(0, -0.17364817766693036, -0.9848077530122082),-Infinity,Infinity), [
				new PCurveEdge(new SemiEllipseCurve(V(0.5, 0.203085322377149, 0),V(0.2, -2.4870779000298388e-17, 0),V(-2.4492935982947065e-17, -0.203085322377149, 0),0,3.141592653589793),V(0.7, 0.203085322377149, 0),V(0.5, 0, 0),0,1.5707963267948966,null,V(-2.4492935982947065e-17, -0.203085322377149, 0),V(-0.2, 1.2435389500149195e-17, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0.5, 0.203085322377149, 0),V(0.2, -2.4870779000298388e-17, 0),V(-2.4492935982947065e-17, -0.203085322377149, 0),0,3.141592653589793),V(0.5, 0, 0),V(0.3, 0.20308532237714905, 0),1.5707963267948966,3.141592653589793,null,V(-0.2, 1.2435389500149195e-17, 0),V(0, 0.203085322377149, 0)),
				new StraightEdge(new L3(V(0.3, 0.11013746176897647, -0.5271335120394901), V(0, 0.17364817766693036, 0.9848077530122081)), V(0.3, 0.20308532237714905, 0), V(0.3, 0.2837856394359068, 0.4576742409727179), 0.535265396141693, 0.9999999999999999),
				new PCurveEdge(new SemiEllipseCurve(V(0.5, 0.2837856394359068, 0.4576742409727179),V(-0.2, 2.4120833250037956e-17, -4.2531536991515426e-18),V(-2.4492935982947065e-17, -0.1969615506024416, 0.034729635533386066),0,3.141592653589793),V(0.3, 0.2837856394359068, 0.4576742409727179),V(0.7, 0.2837856394359068, 0.4576742409727179),0,3.141592653589793,null,V(-2.4492935982947065e-17, -0.1969615506024416, 0.034729635533386066),V(4.898587196589413e-17, 0.1969615506024416, -0.034729635533386066)),
				new StraightEdge(new L3(V(0.7, 0.11013746176897644, -0.5271335120394901), V(0, 0.17364817766693036, 0.9848077530122081)), V(0.7, 0.2837856394359068, 0.4576742409727179), V(0.7, 0.203085322377149, 0), 0.9999999999999999, 0.535265396141693)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, -0.1736481776669304, -0.9848077530122081), -0.5),V(-1, 0, 0),V(0, 0.9848077530122081, -0.1736481776669304)), [
				new StraightEdge(new L3(V(0.3, 0.2837856394359068, 0.4576742409727179), V(1, -1.2060416625018976e-16, 2.1265768495757713e-17)), V(0.5, 0.2837856394359068, 0.4576742409727179), V(0.3, 0.2837856394359068, 0.4576742409727179), 0.2, 0),
				new PCurveEdge(new SemiEllipseCurve(V(0.5, 0.2837856394359068, 0.4576742409727179),V(0.2, 0, 0),V(0, 0.1969615506024416, -0.034729635533386066),0,3.141592653589793),V(0.3, 0.2837856394359068, 0.4576742409727179),V(0.7, 0.2837856394359068, 0.4576742409727179),3.141592653589793,0,null,V(2.4492935982947065e-17, 0.1969615506024416, -0.034729635533386066),V(0, -0.1969615506024416, 0.034729635533386066)),
				new StraightEdge(new L3(V(0.7, 0.2837856394359068, 0.4576742409727179), V(-1, 0, 0)), V(0.7, 0.2837856394359068, 0.4576742409727179), V(0.5, 0.2837856394359068, 0.4576742409727179), 0, 0.2)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, -0.1736481776669304, -0.9848077530122081), -0.5),V(-1, 0, 0),V(0, 0.9848077530122081, -0.1736481776669304)), [
				new StraightEdge(new L3(V(0.7, 0.2837856394359068, 0.4576742409727179), V(-1, 0, 0)), V(0.5, 0.2837856394359068, 0.4576742409727179), V(0.7, 0.2837856394359068, 0.4576742409727179), 0.2, 0),
				new PCurveEdge(new SemiEllipseCurve(V(0.5, 0.2837856394359068, 0.4576742409727179),V(-0.2, 2.4120833250037956e-17, -4.2531536991515426e-18),V(-2.4492935982947065e-17, -0.1969615506024416, 0.034729635533386066),0,3.141592653589793),V(0.7, 0.2837856394359068, 0.4576742409727179),V(0.3, 0.2837856394359068, 0.4576742409727179),3.141592653589793,0,null,V(-4.898587196589413e-17, -0.1969615506024416, 0.034729635533386066),V(2.4492935982947065e-17, 0.1969615506024416, -0.034729635533386066)),
				new StraightEdge(new L3(V(0.3, 0.2837856394359068, 0.4576742409727179), V(1, -1.2060416625018976e-16, 2.1265768495757713e-17)), V(0.3, 0.2837856394359068, 0.4576742409727179), V(0.5, 0.2837856394359068, 0.4576742409727179), 0, 0.2)], [])], false)
		b2EqualAnd(assert, a, b, result)
	},
	'star - ball'(assert) {
		const a = B2T.extrudeEdges(Edge.star(4,2,1),P3.XY,V3.Z.negated())
		const b = B2T.sphere().flipped()
		const result = B2.EMPTY
		b2EqualAnd(assert, a, b, result)
	},
	//'fail'(assert) {
	//	const a = B2T.extrudeEdges([
	//		new PCurveEdge(new BezierCurve(V(-165.04412089048037, 233.67721659018514, 0), V(162.0276742021058, 92.61628384754499, 0), V(-26.526878671628538, 103.9576526746689, 0), V(106.49691429238996, -38.442642927294784, 0), -0.1, 1.1),V(106.49691429238996, -38.442642927294784, 0),V(-165.04412089048037, 233.67721659018514, 0),1,0,null,V(-399.07137889205546, 427.2008868058911, 0),V(-981.2153852777585, 423.18279822792044, 0)),
	//		new PCurveEdge(new BezierCurve(V(-234.98404318638632, -41.414503714148154, 0), V(-263.95799566517644, 2.0493326379292967, 0), V(-96.30326822511165, 115.1793952369976, 0), V(-165.04412089048037, 233.67721659018514, 0), -0.1, 1.1),V(-165.04412089048037, 233.67721659018514, 0),V(-234.98404318638632, -41.414503714148154, 0),1,0,null,V(206.22255799610616, -355.4934640595626, 0),V(86.92185743637037, -130.39150905623234, 0)),
	//		new PCurveEdge(new BezierCurve(V(106.49691429238996, -38.442642927294784, 0), V(6.920773903871436, -168.34765938596584, 0), V(-96.36814809386107, 9.19324831183017, 0), V(-234.98404318638632, -41.414503714148154, 0), -0.1, 1.1),V(-234.98404318638632, -41.414503714148154, 0),V(106.49691429238996, -38.442642927294784, 0),1,0,null,V(415.84768527757575, 151.82325607793496, 0),V(298.7284211655556, 389.7150493760132, 0))
	//	], new P3(V3.Z, 0), V(0, 0, -100), "extrude42")
	//	const b = B2T.extrudeEdges([
	//		new StraightEdge(new L3(V(-113.15294177340922, 90.2593377922355, 0), V(0.9138115486202569, -0.40613846605344794, 0)), V(-13.152941773409395, 45.814893347791084, 0), V(-113.15294177340922, 90.2593377922355, 0), 109.43175335328989, 0),
	//		new StraightEdge(new L3(V(-150.09473256378868, -35.10025916844387, 0), V(0.2826685651896956, 0.9592176407122622, 0)), V(-113.15294177340922, 90.2593377922355, 0), V(-150.09473256378868, -35.10025916844387, 0), 130.68941983551744, 0),
	//		new StraightEdge(new L3(V(-13.152941773409395, 45.814893347791084, 0), V(-0.8609402861542653, -0.5087060287401869, 0)), V(-150.09473256378868, -35.10025916844387, 0), V(-13.152941773409395, 45.814893347791084, 0), 159.06073045098708, 0)
	//	], new P3(V3.Z, 0), V(0, 0, -100), "extrude90").flipped()
	//	const result = B2.EMPTY
	//	b2EqualAnd(assert, a, b, result)
	//},

	async 'sphere() - "a"'(assert) {
		const a = B2T.sphere()
		const b = B2T.text('a', 64, 64, await B2T.loadFont('fonts/FiraSansMedium.woff')).scale(0.5/32).translate(-0.25,-0.25,1.2).flipped()
		const result = new B2([
			new RotationFace(new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.20296874999999998, 0.11703124999999998, 1.2), V(0.20296874999999998, 0.2240625, 1.2), V(0.14500000000000002, 0.2890625, 1.2), V(0.010937499999999989, 0.2890625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.010937499999999989, 0.2890625, 0.9572477433702835), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), -1), V(0.010937499999999989, 0.2890625, 0.9572477433702835), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), 0, 100, null, V(0.004021832534498103, 0, -0.00004595340511453907), V(0, -0.0032085426335165305, 0.0003862505247495081)),
				new PCurveEdge(new SemiEllipseCurve(V(0.20296874999999998, 0, 0), V(0, 0, -0.9791852156376941), V(0, 0.9791852156376941, 0), 0, 3.141592653589793), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), V(0.2029687500000001, 0, 0.9791852156376941), 3.0217872455230754, 3.141592653589793, null, V(0, -0.9721663299286162, 0.1170312499999999), V(0, -0.9791852156376941, 1.199156040103113e-16)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0.2029687500000001, 0, 0.9791852156376941), V(0, 0, -1), 2.937203822794261, 0, null, V(0.9791852156376941, 0, -0.2029687500000001), V(-1, 0, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(0, 0, -1), V(-0.171046939653212, 2.0947208715025655e-17, 0.9852628808776216), 0, 2.969700482947388, null, V(-1, 1.2246467991473532e-16, 0), V(0.9852628808776216, -1.2065990333854792e-16, 0.171046939653212)),
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(-0.1710469396532132, 6.938893903907228e-18, 0.9852628808776215), V(0.04296875, 0.07093749999999999, 0.9965548442595558), -1), V(-0.171046939653212, 2.0947208715025655e-17, 0.9852628808776216), V(0.04296875, 0.07093749999999999, 0.9965548442595558), 0, 63, null, V(0.0018504236016956058, 0.0021729551024794553, 0.0003212434978268739), V(0.0047984440668687395, 0, -0.00020689593220678292)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0.07093749999999999, 0), V(0, 0, 0.9974807622674986), V(0.9974807622674986, 0, 0), 0, 3.141592653589793), V(0.04296875, 0.07093749999999999, 0.9965548442595558), V(0.11093750000000002, 0.07093749999999999, 0.9912924604714292), 0.04309060575909555, 0.11144825166905915, null, V(0.9965548442595558, 0, -0.04296875), V(0.9912924604714292, 0, -0.11093750000000004)),
				new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000002, 0, 0), V(0, 0, 0.9938273849586506), V(0, 0.9938273849586506, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.07093749999999999, 0.9912924604714292), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), 0.07143883874192354, 0.11296719097084001, null, V(0, 0.9912924604714292, -0.07093749999999999), V(0, 0.9874927190198353, -0.11203125000000001)),
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.0029687499999999922, 0.2140625, 1.2), V(0.07406249999999998, 0.2140625, 1.2), V(0.11093750000000002, 0.18593749999999998, 1.2), V(0.11093750000000002, 0.11203125000000003, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), -1), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), 0, 100, null, V(0, 0.0022164863970481866, -0.00025146083296267513), V(-0.0023109368785162285, 0, -0.000007023428862442937)),
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1409375, 0.18703124999999998, 1.2), V(-0.08703125, 0.2059375, 1.2), V(-0.037968749999999996, 0.2140625, 1.2), V(-0.0029687499999999922, 0.2140625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), V(-0.1409375, 0.18703124999999998, 0.9721913045369145), -1), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), V(-0.1409375, 0.18703124999999998, 0.9721913045369144), 0, 100, null, V(-0.001049999946536169, 0, -0.0000031911732417376288), V(-0.0016170602564199374, -0.0005671428725414859, -0.00012531586009935784)),
				new PCurveEdge(new SemiEllipseCurve(V(-0.06640184370318536, -0.0235077791500932, 0), V(0, 0, 0.997516004619599), V(-0.33289781807779234, 0.9403282523625955, 0), 0.02500215049153284, 3.1165905030982604), V(-0.1409375, 0.18703124999999998, 0.9721913045369144), V(-0.16499999999999998, 0.255, 0.9527591510974849), 0.22581372635962282, 0.3006922248899664, null, V(-0.32444628711291806, 0.9164554213903857, -0.22334333850612303), V(-0.31796125684715715, 0.8981373164189183, -0.29544573016418446)),
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.010937499999999989, 0.2890625, 1.2), V(-0.04093749999999999, 0.2890625, 1.2), V(-0.1, 0.2790625, 1.2), V(-0.16499999999999998, 0.255, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(-0.16499999999999998, 0.255, 0.9527591510974849), V(0.010937499999999989, 0.2890625, 0.9572477433702835), -1), V(-0.16499999999999998, 0.255, 0.9527591510974849), V(0.010937499999999989, 0.2890625, 0.9572477433702835), 0, 100, null, V(0.0019497964548261053, 0.0007217996491423563, 0.0001444830042896346), V(0.0015562475240402148, 0, -0.000017781663537028133))], []),
			new RotationFace(new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), V(-0.03286362438479655, -2.6020852139652106e-18, 0.9994598452125503), -1), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), 0, 29, null, V(-0.003388613668575339, 0, 0.00016274305244215148), V(-0.002184873829710333, -0.0006707444983087241, -0.00007184167849434948)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), V(0, 0, 1), 3.108723110778215, 3.141592653589793, null, V(0.9994598452125503, -1.2239853003158589e-16, 0.032863624384797126), V(1, -1.2246467991473532e-16, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0, 0, 1), V(0.11093749999999993, 0, 0.9938273849586506), 3.141592653589793, 3.030426330354509, null, V(1, 0, 0), V(0.9938273849586506, 0, -0.11093749999999993)),
				new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000002, 0, 0), V(0, 0, 0.9938273849586506), V(0, 0.9938273849586506, 0), 0, 3.141592653589793), V(0.11093749999999993, 0, 0.9938273849586506), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), 0, 0.010062279327831384, null, V(0, 0.9938273849586506, 0), V(0, 0.993777073137507, -0.010000000000000007)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0.010000000000000009, 0), V(0, 0, -0.9999499987499375), V(0.9999499987499375, 0, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), 3.0304207486077566, 3.0936030871101416, null, V(-0.9937770731375071, 0, 0.11093749999999983), V(-0.9987987780446257, 0, 0.0479687500000002))], []),
			new RotationFace(new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), V(-0.17104693965321316, 2.0947208715025797e-17, 0.9852628808776215), -1), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), V(-0.171046939653212, 2.0947208715025655e-17, 0.9852628808776216), 0, 36, null, V(0, 0.003270063928473244, 0.00033267819828226424), V(0.0018504199701244795, 0.0021729508379202474, 0.00032124286736657804)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(-0.171046939653212, 2.0947208715025655e-17, 0.9852628808776216), V(0, 0, -1), 2.969700482947388, 0, null, V(-0.9852628808776216, 1.2065990333854792e-16, -0.171046939653212), V(1, -1.2246467991473532e-16, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0, 0, -1), V(0.2029687500000001, 0, 0.9791852156376941), 0, 2.937203822794261, null, V(1, 0, 0), V(-0.9791852156376941, 0, 0.2029687500000001)),
				new PCurveEdge(new SemiEllipseCurve(V(0.20296874999999998, -2.3224450485052117e-18, 0), V(0, 0, 0.9791852156376941), V(-1.1204206832959602e-17, -0.9791852156376941, 0), 2.3013070043406876e-17, 3.141592653589793), V(0.2029687500000001, 0, 0.9791852156376941), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), 2.3013070043406876e-17, 0.13412262887088133, null, V(-1.1204206832959602e-17, -0.9791852156376941, -2.3224450485052132e-18), V(-1.1103582248632314e-17, -0.9703911879325716, -0.13093749999999996)),
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.24203124999999998, -0.19796875, 1.2), V(0.21500000000000002, -0.18796875, 1.2), V(0.20296874999999998, -0.17296875, 1.2), V(0.20296874999999998, -0.13093749999999998, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), V(0.24203124999999998, -0.19796875, 0.9498574882827818), -1), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), V(0.24203124999999998, -0.19796875, 0.9498574882827818), 0, 100, null, V(0, -0.0012607550302782095, -0.00017011707631925008), V(0.0008106296264206332, -0.0002998861045717958, -0.0002690569713737753)),
				new PCurveEdge(new SemiEllipseCurve(V(0.27731540188902115, -0.09536944308866366, 0), V(0.31091053038645644, 0.9040660812655796, 0), V(0, 0, 0.9560339100680944), 1.465110231866663, 3.141592653589793), V(0.24203124999999998, -0.19796875, 0.9498574882827818), V(0.21999999999999997, -0.26203125, 0.9396486705271484), 1.684527864719763, 1.7562036900306723, null, V(-0.30890190438173115, -0.8982253957199245, -0.10849695458036647), V(-0.30558190818745745, -0.8885715060770011, -0.17624191662783179)),
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.12296875000000002, -0.18796875, 1.2), V(0.13906249999999998, -0.235, 1.2), V(0.17296875, -0.2559375, 1.2), V(0.21999999999999997, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.21999999999999997, -0.26203125, 0.9396486705271484), V(0.12296875000000002, -0.18796875, 0.9744467330474638), -1), V(0.21999999999999997, -0.26203125, 0.9396486705271484), V(0.12296875000000002, -0.18796875, 0.9744467330474638), 0, 100, null, V(-0.001409912819409077, 0.0001826797340762595, 0.0003810452039896991), V(-0.00048254782297658694, 0.0014101640263684755, 0.0003329120627492527)),
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.034062499999999996, -0.26203125, 1.2), V(0.030000000000000027, -0.26203125, 1.2), V(0.0859375, -0.23796875, 1.2), V(0.12296875000000002, -0.18796875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.12296875000000002, -0.18796875, 0.9744467330474638), V(-0.034062499999999996, -0.26203125, 0.9644580706890205), -1), V(0.12296875000000002, -0.18796875, 0.9744467330474638), V(-0.034062499999999996, -0.26203125, 0.9644580706890203), 0, 100, null, V(-0.0011108139461180212, -0.0014998331761931086, -0.0001491373104062419), V(-0.001921830714021067, 0, -0.0000678747585673429)),
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.20500000000000002, -0.0990625, 1.2), V(-0.20500000000000002, -0.19703125, 1.2), V(-0.14, -0.26203125, 1.2), V(-0.034062499999999996, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.034062499999999996, -0.26203125, 0.9644580706890205), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), -1), V(-0.034062499999999996, -0.26203125, 0.9644580706890203), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), 0, 100, null, V(-0.0031779248165223913, 0, -0.0001122371903482131), V(0, 0.0029377479944999226, 0.00029887021513792583))], []),
			new RotationFace(new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.009062500000000001, -0.19296875, 1.2), V(-0.07500000000000001, -0.19296875, 1.2), V(-0.1040625, -0.16, 1.2), V(-0.1040625, -0.095, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.1040625, -0.095, 0.9900232300778351), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), -1), V(-0.1040625, -0.095, 0.9900232300778351), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), 0, 100, null, V(0, -0.001949658715744448, -0.00018708407274560705), V(0.0019781216920816535, 0, 0.000018270895766511516)),
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.11093750000000002, -0.12, 1.2), V(0.0809375, -0.16593750000000002, 1.2), V(0.040937500000000016, -0.19296875, 1.2), V(-0.009062500000000001, -0.19296875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), V(0.11093750000000002, -0.12, 0.9865560658643532), -1), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), V(0.11093750000000002, -0.12, 0.9865560658643532), 0, 100, null, V(0.0014999985603503193, 0, 0.000013854717561505174), V(0.0008999801702201017, 0.00137809463564953, 0.00006642278975471942)),
				new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000002, 2.7017833632379675e-19, 0), V(0, 0, -0.9938273849586506), V(2.4203775050019737e-18, -0.9938273849586506, 0), 1.3942163371701867e-17, 3.141592653589793), V(0.11093750000000002, -0.12, 0.9865560658643532), V(0.11093749999999993, 0, 0.9938273849586506), 3.02055199779192, 3.141592653589793, null, V(-2.4026688591808874e-18, 0.9865560658643532, 0.11999999999999994), V(-2.4203775050019737e-18, 0.9938273849586506, 1.217087525894596e-16)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0.11093749999999993, 0, 0.9938273849586506), V(0, 0, 1), 3.030426330354509, 3.141592653589793, null, V(-0.9938273849586506, 0, 0.11093749999999993), V(-1, 0, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(0, 0, 1), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), 3.141592653589793, 3.108723110778215, null, V(-1, 1.2246467991473532e-16, 0), V(-0.9994598452125503, 1.2239853003158589e-16, -0.032863624384797126)),
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.03286362438479653, 2.7235906341395923e-18, 0.9994598452125503), V(-0.1040625, -0.095, 0.9900232300778351), -1), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), V(-0.1040625, -0.095, 0.9900232300778351), 0, 70, null, V(-0.0021848918014005826, -0.0006707500155208376, -0.00007184226942842322), V(0, -0.0019213481561159574, -0.00018436746662668277))], []),
			new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.20296874999999998, 0.11703124999999998, 1.2), V(0.20296874999999998, 0.2240625, 1.2), V(0.14500000000000002, 0.2890625, 1.2), V(0.010937499999999989, 0.2890625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.20296874999999998, 0.11703124999999998, 1.2), V(0.20296874999999998, 0.2240625, 1.2), V(0.14500000000000002, 0.2890625, 1.2), V(0.010937499999999989, 0.2890625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.010937499999999989, 0.2890625, 0.9572477433702835), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), -1), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), V(0.010937499999999989, 0.2890625, 0.9572477433702835), 100, 0, null, V(0, 0.0032085426335165305, -0.0003862505247495081), V(-0.004021832534498103, 0, 0.00004595340511453907)),
				new StraightEdge(new L3(V(0.010937499999999989, 0.2890625, 1.2), V(0, 0, -1)), V(0.010937499999999989, 0.2890625, 0.9572477433702835), V(0.010937499999999989, 0.2890625, 0.19999999999999996), 0.24275225662971645, 1),
				new PCurveEdge(new BezierCurve(V(0.20296874999999998, 0.11703124999999998, 0.19999999999999996), V(0.20296874999999998, 0.2240625, 0.19999999999999996), V(0.14500000000000002, 0.2890625, 0.19999999999999996), V(0.010937499999999989, 0.2890625, 0.19999999999999996), 0, 1), V(0.010937499999999989, 0.2890625, 0.19999999999999996), V(0.20296874999999998, 0.11703124999999998, 0.19999999999999996), 1, 0, null, V(0.4021875000000001, 0, 0), V(0, -0.32109375000000007, 0)),
				new StraightEdge(new L3(V(0.20296874999999998, 0.11703124999999998, 1.2), V(0, 0, -1)), V(0.20296874999999998, 0.11703124999999998, 0.19999999999999996), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), 1, 0.22783367007138378)], []),
			new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.010937499999999989, 0.2890625, 1.2), V(-0.04093749999999999, 0.2890625, 1.2), V(-0.1, 0.2790625, 1.2), V(-0.16499999999999998, 0.255, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.010937499999999989, 0.2890625, 1.2), V(-0.04093749999999999, 0.2890625, 1.2), V(-0.1, 0.2790625, 1.2), V(-0.16499999999999998, 0.255, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(-0.16499999999999998, 0.255, 0.9527591510974849), V(0.010937499999999989, 0.2890625, 0.9572477433702835), -1), V(0.010937499999999989, 0.2890625, 0.9572477433702835), V(-0.16499999999999998, 0.255, 0.9527591510974849), 100, 0, null, V(-0.0015562475240402148, 0, 0.000017781663537028133), V(-0.0019497964548261053, -0.0007217996491423563, -0.0001444830042896346)),
				new StraightEdge(new L3(V(-0.16499999999999998, 0.255, 1.2), V(0, 0, -1)), V(-0.16499999999999998, 0.255, 0.9527591510974849), V(-0.16499999999999998, 0.255, 0.19999999999999996), 0.24724084890251508, 1),
				new PCurveEdge(new BezierCurve(V(0.010937499999999989, 0.2890625, 0.19999999999999996), V(-0.04093749999999999, 0.2890625, 0.19999999999999996), V(-0.1, 0.2790625, 0.19999999999999996), V(-0.16499999999999998, 0.255, 0.19999999999999996), 0, 1), V(-0.16499999999999998, 0.255, 0.19999999999999996), V(0.010937499999999989, 0.2890625, 0.19999999999999996), 1, 0, null, V(0.19499999999999995, 0.07218749999999996, 0), V(0.15562499999999993, 0, 0)),
				new StraightEdge(new L3(V(0.010937499999999989, 0.2890625, 1.2), V(0, 0, -1)), V(0.010937499999999989, 0.2890625, 0.19999999999999996), V(0.010937499999999989, 0.2890625, 0.9572477433702835), 1, 0.24275225662971645)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0.942669839890126, 0.33372679389213644, 0), -0.0704401911393759), V(0, 0, -1), V(-0.33372679389213644, 0.942669839890126, 0)), [
				new PCurveEdge(new SemiEllipseCurve(V(-0.06640184370318536, -0.0235077791500932, 0), V(0, 0, 0.997516004619599), V(-0.33289781807779234, 0.9403282523625955, 0), 0.02500215049153284, 3.1165905030982604), V(-0.16499999999999998, 0.255, 0.9527591510974849), V(-0.1409375, 0.18703124999999998, 0.9721913045369144), 0.3006922248899664, 0.22581372635962282, null, V(0.31796125684715715, -0.8981373164189183, 0.29544573016418446), V(0.32444628711291806, -0.9164554213903857, 0.22334333850612303)),
				new StraightEdge(new L3(V(-0.1409375, 0.18703124999999998, 1.2), V(0, 0, -1)), V(-0.1409375, 0.18703124999999998, 0.9721913045369144), V(-0.1409375, 0.18703124999999998, 0.19999999999999996), 0.22780869546308558, 1),
				new StraightEdge(new L3(V(-0.16499999999999998, 0.255, 0.19999999999999996), V(0.33372679389213644, -0.942669839890126, 0)), V(-0.1409375, 0.18703124999999998, 0.19999999999999996), V(-0.16499999999999998, 0.255, 0.19999999999999996), 0.07210239165806154, 0),
				new StraightEdge(new L3(V(-0.16499999999999998, 0.255, 1.2), V(0, 0, -1)), V(-0.16499999999999998, 0.255, 0.19999999999999996), V(-0.16499999999999998, 0.255, 0.9527591510974849), 1, 0.24724084890251508)], []),
			new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(-0.1409375, 0.18703124999999998, 1.2), V(-0.08703125, 0.2059375, 1.2), V(-0.037968749999999996, 0.2140625, 1.2), V(-0.0029687499999999922, 0.2140625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1409375, 0.18703124999999998, 1.2), V(-0.08703125, 0.2059375, 1.2), V(-0.037968749999999996, 0.2140625, 1.2), V(-0.0029687499999999922, 0.2140625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), V(-0.1409375, 0.18703124999999998, 0.9721913045369145), -1), V(-0.1409375, 0.18703124999999998, 0.9721913045369144), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), 100, 0, null, V(0.0016170602564199374, 0.0005671428725414859, 0.00012531586009935784), V(0.001049999946536169, 0, 0.0000031911732417376288)),
				new StraightEdge(new L3(V(-0.0029687499999999922, 0.2140625, 1.2), V(0, 0, -1)), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), V(-0.0029687499999999922, 0.2140625, 0.19999999999999996), 0.22318454526088316, 1),
				new PCurveEdge(new BezierCurve(V(-0.1409375, 0.18703124999999998, 0.19999999999999996), V(-0.08703125, 0.2059375, 0.19999999999999996), V(-0.037968749999999996, 0.2140625, 0.19999999999999996), V(-0.0029687499999999922, 0.2140625, 0.19999999999999996), 0, 1), V(-0.0029687499999999922, 0.2140625, 0.19999999999999996), V(-0.1409375, 0.18703124999999998, 0.19999999999999996), 1, 0, null, V(-0.10500000000000001, 0, 0), V(-0.16171874999999997, -0.05671875000000004, 0)),
				new StraightEdge(new L3(V(-0.1409375, 0.18703124999999998, 1.2), V(0, 0, -1)), V(-0.1409375, 0.18703124999999998, 0.19999999999999996), V(-0.1409375, 0.18703124999999998, 0.9721913045369144), 1, 0.22780869546308558)], []),
			new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(-0.0029687499999999922, 0.2140625, 1.2), V(0.07406249999999998, 0.2140625, 1.2), V(0.11093750000000002, 0.18593749999999998, 1.2), V(0.11093750000000002, 0.11203125000000003, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.0029687499999999922, 0.2140625, 1.2), V(0.07406249999999998, 0.2140625, 1.2), V(0.11093750000000002, 0.18593749999999998, 1.2), V(0.11093750000000002, 0.11203125000000003, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), -1), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), 100, 0, null, V(0.0023109368785162285, 0, 0.000007023428862442937), V(0, -0.0022164863970481866, 0.00025146083296267513)),
				new StraightEdge(new L3(V(0.11093750000000002, 0.11203125000000003, 1.2), V(0, 0, -1)), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), V(0.11093750000000002, 0.11203125000000003, 0.19999999999999996), 0.21250728098016458, 1),
				new PCurveEdge(new BezierCurve(V(-0.0029687499999999922, 0.2140625, 0.19999999999999996), V(0.07406249999999998, 0.2140625, 0.19999999999999996), V(0.11093750000000002, 0.18593749999999998, 0.19999999999999996), V(0.11093750000000002, 0.11203125000000003, 0.19999999999999996), 0, 1), V(0.11093750000000002, 0.11203125000000003, 0.19999999999999996), V(-0.0029687499999999922, 0.2140625, 0.19999999999999996), 1, 0, null, V(0, 0.22171874999999985, 0), V(-0.2310937499999999, 0, 0)),
				new StraightEdge(new L3(V(-0.0029687499999999922, 0.2140625, 1.2), V(0, 0, -1)), V(-0.0029687499999999922, 0.2140625, 0.19999999999999996), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), 1, 0.22318454526088316)], []),
			new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 0.11093750000000002), V(0, 0, -1), V(0, 1, 0)), [
				new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000002, 0, 0), V(0, 0, 0.9938273849586506), V(0, 0.9938273849586506, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), V(0.11093750000000002, 0.07093749999999999, 0.9912924604714292), 0.11296719097084001, 0.07143883874192354, null, V(0, -0.9874927190198353, 0.11203125000000001), V(0, -0.9912924604714292, 0.07093749999999999)),
				new StraightEdge(new L3(V(0.11093750000000002, 0.07093749999999999, 1.2), V(0, 0, -1)), V(0.11093750000000002, 0.07093749999999999, 0.9912924604714292), V(0.11093750000000002, 0.07093749999999999, 0.19999999999999996), 0.2087075395285708, 1),
				new StraightEdge(new L3(V(0.11093750000000002, 0.11203125000000003, 0.19999999999999996), V(0, -1, 0)), V(0.11093750000000002, 0.07093749999999999, 0.19999999999999996), V(0.11093750000000002, 0.11203125000000003, 0.19999999999999996), 0.04109375000000004, 0),
				new StraightEdge(new L3(V(0.11093750000000002, 0.11203125000000003, 1.2), V(0, 0, -1)), V(0.11093750000000002, 0.11203125000000003, 0.19999999999999996), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), 1, 0.21250728098016458)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), -0.07093749999999999), V(0, 0, -1), V(1, 0, 0)), [
				new PCurveEdge(new SemiEllipseCurve(V(0, 0.07093749999999999, 0), V(0, 0, 0.9974807622674986), V(0.9974807622674986, 0, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.07093749999999999, 0.9912924604714292), V(0.04296875, 0.07093749999999999, 0.9965548442595558), 0.11144825166905915, 0.04309060575909555, null, V(-0.9912924604714292, 0, 0.11093750000000004), V(-0.9965548442595558, 0, 0.04296875)),
				new StraightEdge(new L3(V(0.04296875, 0.07093749999999999, 1.2), V(0, 0, -1)), V(0.04296875, 0.07093749999999999, 0.9965548442595558), V(0.04296875, 0.07093749999999999, 0.19999999999999996), 0.2034451557404442, 1),
				new StraightEdge(new L3(V(0.11093750000000002, 0.07093749999999999, 0.19999999999999996), V(-1, 0, 0)), V(0.04296875, 0.07093749999999999, 0.19999999999999996), V(0.11093750000000002, 0.07093749999999999, 0.19999999999999996), 0.06796875000000002, 0),
				new StraightEdge(new L3(V(0.11093750000000002, 0.07093749999999999, 1.2), V(0, 0, -1)), V(0.11093750000000002, 0.07093749999999999, 0.19999999999999996), V(0.11093750000000002, 0.07093749999999999, 0.9912924604714292), 1, 0.2087075395285708)], []),
			new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(-0.1710469396532132, 6.938893903907228e-18, 0.9852628808776215), V(0.04296875, 0.07093749999999999, 0.9965548442595558), -1), V(0.04296875, 0.07093749999999999, 0.9965548442595558), V(-0.171046939653212, 2.0947208715025655e-17, 0.9852628808776216), 63, 0, null, V(-0.0047984440668687395, 0, 0.00020689593220678292), V(-0.0018504236016956058, -0.0021729551024794553, -0.0003212434978268739)),
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), V(-0.17104693965321316, 2.0947208715025797e-17, 0.9852628808776215), -1), V(-0.171046939653212, 2.0947208715025655e-17, 0.9852628808776216), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), 36, 0, null, V(-0.0018504199701244795, -0.0021729508379202474, -0.00032124286736657804), V(0, -0.003270063928473244, -0.00033267819828226424)),
				new StraightEdge(new L3(V(-0.20500000000000002, -0.0990625, 1.2), V(0, 0, -1)), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), V(-0.20500000000000002, -0.0990625, 0.19999999999999996), 0.22626409068282272, 1),
				new PCurveEdge(new BezierCurve(V(0.04296875, 0.07093749999999999, 0.19999999999999996), V(-0.11703125, 0.07093749999999999, 0.19999999999999996), V(-0.20500000000000002, 0.010000000000000009, 0.19999999999999996), V(-0.20500000000000002, -0.0990625, 0.19999999999999996), 0, 1), V(-0.20500000000000002, -0.0990625, 0.19999999999999996), V(0.04296875, 0.07093749999999999, 0.19999999999999996), 1, 0, null, V(0, 0.3271875, 0), V(0.48, 0, 0)),
				new StraightEdge(new L3(V(0.04296875, 0.07093749999999999, 1.2), V(0, 0, -1)), V(0.04296875, 0.07093749999999999, 0.19999999999999996), V(0.04296875, 0.07093749999999999, 0.9965548442595558), 1, 0.2034451557404442)], []),
			new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(-0.20500000000000002, -0.0990625, 1.2), V(-0.20500000000000002, -0.19703125, 1.2), V(-0.14, -0.26203125, 1.2), V(-0.034062499999999996, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.20500000000000002, -0.0990625, 1.2), V(-0.20500000000000002, -0.19703125, 1.2), V(-0.14, -0.26203125, 1.2), V(-0.034062499999999996, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.034062499999999996, -0.26203125, 0.9644580706890205), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), -1), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), V(-0.034062499999999996, -0.26203125, 0.9644580706890203), 100, 0, null, V(0, -0.0029377479944999226, -0.00029887021513792583), V(0.0031779248165223913, 0, 0.0001122371903482131)),
				new StraightEdge(new L3(V(-0.034062499999999996, -0.26203125, 1.2), V(0, 0, -1)), V(-0.034062499999999996, -0.26203125, 0.9644580706890203), V(-0.034062499999999996, -0.26203125, 0.19999999999999996), 0.23554192931097961, 1),
				new PCurveEdge(new BezierCurve(V(-0.20500000000000002, -0.0990625, 0.19999999999999996), V(-0.20500000000000002, -0.19703125, 0.19999999999999996), V(-0.14, -0.26203125, 0.19999999999999996), V(-0.034062499999999996, -0.26203125, 0.19999999999999996), 0, 1), V(-0.034062499999999996, -0.26203125, 0.19999999999999996), V(-0.20500000000000002, -0.0990625, 0.19999999999999996), 1, 0, null, V(-0.3178125, 0, 0), V(0, 0.29390625, 0)),
				new StraightEdge(new L3(V(-0.20500000000000002, -0.0990625, 1.2), V(0, 0, -1)), V(-0.20500000000000002, -0.0990625, 0.19999999999999996), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), 1, 0.22626409068282272)], []),
			new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(-0.034062499999999996, -0.26203125, 1.2), V(0.030000000000000027, -0.26203125, 1.2), V(0.0859375, -0.23796875, 1.2), V(0.12296875000000002, -0.18796875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.034062499999999996, -0.26203125, 1.2), V(0.030000000000000027, -0.26203125, 1.2), V(0.0859375, -0.23796875, 1.2), V(0.12296875000000002, -0.18796875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.12296875000000002, -0.18796875, 0.9744467330474638), V(-0.034062499999999996, -0.26203125, 0.9644580706890205), -1), V(-0.034062499999999996, -0.26203125, 0.9644580706890203), V(0.12296875000000002, -0.18796875, 0.9744467330474638), 100, 0, null, V(0.001921830714021067, 0, 0.0000678747585673429), V(0.0011108139461180212, 0.0014998331761931086, 0.0001491373104062419)),
				new StraightEdge(new L3(V(0.12296875000000002, -0.18796875, 1.2), V(0, 0, -1)), V(0.12296875000000002, -0.18796875, 0.9744467330474638), V(0.12296875000000002, -0.18796875, 0.19999999999999996), 0.2255532669525362, 1),
				new PCurveEdge(new BezierCurve(V(-0.034062499999999996, -0.26203125, 0.19999999999999996), V(0.030000000000000027, -0.26203125, 0.19999999999999996), V(0.0859375, -0.23796875, 0.19999999999999996), V(0.12296875000000002, -0.18796875, 0.19999999999999996), 0, 1), V(0.12296875000000002, -0.18796875, 0.19999999999999996), V(-0.034062499999999996, -0.26203125, 0.19999999999999996), 1, 0, null, V(-0.11109375000000005, -0.15000000000000002, 0), V(-0.19218750000000007, 0, 0)),
				new StraightEdge(new L3(V(-0.034062499999999996, -0.26203125, 1.2), V(0, 0, -1)), V(-0.034062499999999996, -0.26203125, 0.19999999999999996), V(-0.034062499999999996, -0.26203125, 0.9644580706890203), 1, 0.23554192931097961)], []),
			new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.12296875000000002, -0.18796875, 1.2), V(0.13906249999999998, -0.235, 1.2), V(0.17296875, -0.2559375, 1.2), V(0.21999999999999997, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.12296875000000002, -0.18796875, 1.2), V(0.13906249999999998, -0.235, 1.2), V(0.17296875, -0.2559375, 1.2), V(0.21999999999999997, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.21999999999999997, -0.26203125, 0.9396486705271484), V(0.12296875000000002, -0.18796875, 0.9744467330474638), -1), V(0.12296875000000002, -0.18796875, 0.9744467330474638), V(0.21999999999999997, -0.26203125, 0.9396486705271484), 100, 0, null, V(0.00048254782297658694, -0.0014101640263684755, -0.0003329120627492527), V(0.001409912819409077, -0.0001826797340762595, -0.0003810452039896991)),
				new StraightEdge(new L3(V(0.21999999999999997, -0.26203125, 1.2), V(0, 0, -1)), V(0.21999999999999997, -0.26203125, 0.9396486705271484), V(0.21999999999999997, -0.26203125, 0.19999999999999996), 0.26035132947285156, 1),
				new PCurveEdge(new BezierCurve(V(0.12296875000000002, -0.18796875, 0.19999999999999996), V(0.13906249999999998, -0.235, 0.19999999999999996), V(0.17296875, -0.2559375, 0.19999999999999996), V(0.21999999999999997, -0.26203125, 0.19999999999999996), 0, 1), V(0.21999999999999997, -0.26203125, 0.19999999999999996), V(0.12296875000000002, -0.18796875, 0.19999999999999996), 1, 0, null, V(-0.1410937499999999, 0.01828125, 0), V(-0.04828124999999989, 0.14109375000000002, 0)),
				new StraightEdge(new L3(V(0.12296875000000002, -0.18796875, 1.2), V(0, 0, -1)), V(0.12296875000000002, -0.18796875, 0.19999999999999996), V(0.12296875000000002, -0.18796875, 0.9744467330474638), 1, 0.2255532669525362)], []),
			new PlaneFace(new PlaneSurface(new P3(V(-0.9456422745519422, 0.325208684662986, 0), -0.2932561385545256), V(0, 0, -1), V(-0.325208684662986, -0.9456422745519422, 0)), [
				new PCurveEdge(new SemiEllipseCurve(V(0.27731540188902115, -0.09536944308866366, 0), V(0.31091053038645644, 0.9040660812655796, 0), V(0, 0, 0.9560339100680944), 1.465110231866663, 3.141592653589793), V(0.21999999999999997, -0.26203125, 0.9396486705271484), V(0.24203124999999998, -0.19796875, 0.9498574882827818), 1.7562036900306723, 1.684527864719763, null, V(0.30558190818745745, 0.8885715060770011, 0.17624191662783179), V(0.30890190438173115, 0.8982253957199245, 0.10849695458036647)),
				new StraightEdge(new L3(V(0.24203124999999998, -0.19796875, 1.2), V(0, 0, -1)), V(0.24203124999999998, -0.19796875, 0.9498574882827818), V(0.24203124999999998, -0.19796875, 0.19999999999999996), 0.25014251171721813, 1),
				new StraightEdge(new L3(V(0.21999999999999997, -0.26203125, 0.19999999999999996), V(0.32520868466298514, 0.9456422745519424, 0)), V(0.24203124999999998, -0.19796875, 0.19999999999999996), V(0.21999999999999997, -0.26203125, 0.19999999999999996), 0.06774496204746519, 0),
				new StraightEdge(new L3(V(0.21999999999999997, -0.26203125, 1.2), V(0, 0, -1)), V(0.21999999999999997, -0.26203125, 0.19999999999999996), V(0.21999999999999997, -0.26203125, 0.9396486705271484), 1, 0.26035132947285156)], []),
			new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.24203124999999998, -0.19796875, 1.2), V(0.21500000000000002, -0.18796875, 1.2), V(0.20296874999999998, -0.17296875, 1.2), V(0.20296874999999998, -0.13093749999999998, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.24203124999999998, -0.19796875, 1.2), V(0.21500000000000002, -0.18796875, 1.2), V(0.20296874999999998, -0.17296875, 1.2), V(0.20296874999999998, -0.13093749999999998, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), V(0.24203124999999998, -0.19796875, 0.9498574882827818), -1), V(0.24203124999999998, -0.19796875, 0.9498574882827818), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), 100, 0, null, V(-0.0008106296264206332, 0.0002998861045717958, 0.0002690569713737753), V(0, 0.0012607550302782095, 0.00017011707631925008)),
				new StraightEdge(new L3(V(0.20296874999999998, -0.13093749999999998, 1.2), V(0, 0, -1)), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), V(0.20296874999999998, -0.13093749999999998, 0.19999999999999996), 0.22960881206742834, 1),
				new PCurveEdge(new BezierCurve(V(0.24203124999999998, -0.19796875, 0.19999999999999996), V(0.21500000000000002, -0.18796875, 0.19999999999999996), V(0.20296874999999998, -0.17296875, 0.19999999999999996), V(0.20296874999999998, -0.13093749999999998, 0.19999999999999996), 0, 1), V(0.20296874999999998, -0.13093749999999998, 0.19999999999999996), V(0.24203124999999998, -0.19796875, 0.19999999999999996), 1, 0, null, V(0, -0.12609375, 0), V(0.08109374999999985, -0.030000000000000006, 0)),
				new StraightEdge(new L3(V(0.24203124999999998, -0.19796875, 1.2), V(0, 0, -1)), V(0.24203124999999998, -0.19796875, 0.19999999999999996), V(0.24203124999999998, -0.19796875, 0.9498574882827818), 1, 0.25014251171721813)], []),
			new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), -0.20296874999999998), V(0, 0, -1), V(0, -1, 0)), [
				new PCurveEdge(new SemiEllipseCurve(V(0.20296874999999998, 0, 0), V(0, 0, -0.9791852156376941), V(0, 0.9791852156376941, 0), 0, 3.141592653589793), V(0.2029687500000001, 0, 0.9791852156376941), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), 3.141592653589793, 3.0217872455230754, null, V(0, 0.9791852156376941, -1.199156040103113e-16), V(0, 0.9721663299286162, -0.1170312499999999)),
				new StraightEdge(new L3(V(0.20296874999999998, 0.11703124999999998, 1.2), V(0, 0, -1)), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), V(0.20296874999999998, 0.11703124999999998, 0.19999999999999996), 0.22783367007138378, 1),
				new StraightEdge(new L3(V(0.20296874999999998, -0.13093749999999998, 0.19999999999999996), V(0, 1, 0)), V(0.20296874999999998, 0.11703124999999998, 0.19999999999999996), V(0.20296874999999998, -0.13093749999999998, 0.19999999999999996), 0.24796874999999996, 0),
				new StraightEdge(new L3(V(0.20296874999999998, -0.13093749999999998, 1.2), V(0, 0, -1)), V(0.20296874999999998, -0.13093749999999998, 0.19999999999999996), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), 1, 0.22960881206742834),
				new PCurveEdge(new SemiEllipseCurve(V(0.20296874999999998, -2.3224450485052117e-18, 0), V(0, 0, 0.9791852156376941), V(-1.1204206832959602e-17, -0.9791852156376941, 0), 2.3013070043406876e-17, 3.141592653589793), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), V(0.2029687500000001, 0, 0.9791852156376941), 0.13412262887088133, 2.3013070043406876e-17, null, V(1.1103582248632314e-17, 0.9703911879325716, 0.13093749999999996), V(1.1204206832959602e-17, 0.9791852156376941, 2.3224450485052132e-18))], []),
			new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(-0.009062500000000001, -0.19296875, 1.2), V(-0.07500000000000001, -0.19296875, 1.2), V(-0.1040625, -0.16, 1.2), V(-0.1040625, -0.095, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.009062500000000001, -0.19296875, 1.2), V(-0.07500000000000001, -0.19296875, 1.2), V(-0.1040625, -0.16, 1.2), V(-0.1040625, -0.095, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.1040625, -0.095, 0.9900232300778351), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), -1), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), V(-0.1040625, -0.095, 0.9900232300778351), 100, 0, null, V(-0.0019781216920816535, 0, -0.000018270895766511516), V(0, 0.001949658715744448, 0.00018708407274560705)),
				new StraightEdge(new L3(V(-0.1040625, -0.095, 1.2), V(0, 0, -1)), V(-0.1040625, -0.095, 0.9900232300778351), V(-0.1040625, -0.095, 0.19999999999999996), 0.20997676992216485, 1),
				new PCurveEdge(new BezierCurve(V(-0.009062500000000001, -0.19296875, 0.19999999999999996), V(-0.07500000000000001, -0.19296875, 0.19999999999999996), V(-0.1040625, -0.16, 0.19999999999999996), V(-0.1040625, -0.095, 0.19999999999999996), 0, 1), V(-0.1040625, -0.095, 0.19999999999999996), V(-0.009062500000000001, -0.19296875, 0.19999999999999996), 1, 0, null, V(0, -0.195, 0), V(0.19781250000000003, 0, 0)),
				new StraightEdge(new L3(V(-0.009062500000000001, -0.19296875, 1.2), V(0, 0, -1)), V(-0.009062500000000001, -0.19296875, 0.19999999999999996), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), 1, 0.21883694901551276)], []),
			new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), V(-0.03286362438479655, -2.6020852139652106e-18, 0.9994598452125503), -1), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), 29, 0, null, V(0.002184873829710333, 0.0006707444983087241, 0.00007184167849434948), V(0.003388613668575339, 0, -0.00016274305244215148)),
				new StraightEdge(new L3(V(0.047968750000000004, 0.010000000000000009, 1.2), V(0, 0, -1)), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), V(0.047968750000000004, 0.010000000000000009, 0.19999999999999996), 0.20120122195537427, 1),
				new PCurveEdge(new BezierCurve(V(-0.1040625, -0.095, 0.19999999999999996), V(-0.1040625, -0.030937500000000007, 0.19999999999999996), V(-0.065, 0.010000000000000009, 0.19999999999999996), V(0.047968750000000004, 0.010000000000000009, 0.19999999999999996), 0, 1), V(0.047968750000000004, 0.010000000000000009, 0.19999999999999996), V(-0.1040625, -0.095, 0.19999999999999996), 1, 0, null, V(-0.33890625, 0, 0), V(0, -0.19218749999999998, 0)),
				new StraightEdge(new L3(V(-0.1040625, -0.095, 1.2), V(0, 0, -1)), V(-0.1040625, -0.095, 0.19999999999999996), V(-0.1040625, -0.095, 0.9900232300778351), 1, 0.20997676992216485),
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.03286362438479653, 2.7235906341395923e-18, 0.9994598452125503), V(-0.1040625, -0.095, 0.9900232300778351), -1), V(-0.1040625, -0.095, 0.9900232300778351), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), 70, 0, null, V(0, 0.0019213481561159574, 0.00018436746662668277), V(0.0021848918014005826, 0.0006707500155208376, 0.00007184226942842322))], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 1, 0), 0.010000000000000009), V(0, 0, -1), V(-1, 0, 0)), [
				new PCurveEdge(new SemiEllipseCurve(V(0, 0.010000000000000009, 0), V(0, 0, -0.9999499987499375), V(0.9999499987499375, 0, 0), 0, 3.141592653589793), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), 3.0936030871101416, 3.0304207486077566, null, V(0.9987987780446257, 0, -0.0479687500000002), V(0.9937770731375071, 0, -0.11093749999999983)),
				new StraightEdge(new L3(V(0.11093750000000002, 0.010000000000000009, 1.2), V(0, 0, -1)), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), V(0.11093750000000002, 0.010000000000000009, 0.19999999999999996), 0.20622292686249288, 1),
				new StraightEdge(new L3(V(0.047968750000000004, 0.010000000000000009, 0.19999999999999996), V(1, 0, 0)), V(0.11093750000000002, 0.010000000000000009, 0.19999999999999996), V(0.047968750000000004, 0.010000000000000009, 0.19999999999999996), 0.06296875000000002, 0),
				new StraightEdge(new L3(V(0.047968750000000004, 0.010000000000000009, 1.2), V(0, 0, -1)), V(0.047968750000000004, 0.010000000000000009, 0.19999999999999996), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), 1, 0.20120122195537427)], []),
			new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 0.11093750000000002), V(0, 0, -1), V(0, 1, 0)), [
				new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000002, 0, 0), V(0, 0, 0.9938273849586506), V(0, 0.9938273849586506, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), V(0.11093749999999993, 0, 0.9938273849586506), 0.010062279327831384, 0, null, V(0, -0.993777073137507, 0.010000000000000007), V(0, -0.9938273849586506, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000002, 2.7017833632379675e-19, 0), V(0, 0, -0.9938273849586506), V(2.4203775050019737e-18, -0.9938273849586506, 0), 1.3942163371701867e-17, 3.141592653589793), V(0.11093749999999993, 0, 0.9938273849586506), V(0.11093750000000002, -0.12, 0.9865560658643532), 3.141592653589793, 3.02055199779192, null, V(2.4203775050019737e-18, -0.9938273849586506, -1.217087525894596e-16), V(2.4026688591808874e-18, -0.9865560658643532, -0.11999999999999994)),
				new StraightEdge(new L3(V(0.11093750000000002, -0.12, 1.2), V(0, 0, -1)), V(0.11093750000000002, -0.12, 0.9865560658643532), V(0.11093750000000002, -0.12, 0.19999999999999996), 0.21344393413564677, 1),
				new StraightEdge(new L3(V(0.11093750000000002, 0.010000000000000009, 0.19999999999999996), V(0, -1, 0)), V(0.11093750000000002, -0.12, 0.19999999999999996), V(0.11093750000000002, 0.010000000000000009, 0.19999999999999996), 0.13, 0),
				new StraightEdge(new L3(V(0.11093750000000002, 0.010000000000000009, 1.2), V(0, 0, -1)), V(0.11093750000000002, 0.010000000000000009, 0.19999999999999996), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), 1, 0.20622292686249288)], []),
			new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.11093750000000002, -0.12, 1.2), V(0.0809375, -0.16593750000000002, 1.2), V(0.040937500000000016, -0.19296875, 1.2), V(-0.009062500000000001, -0.19296875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.11093750000000002, -0.12, 1.2), V(0.0809375, -0.16593750000000002, 1.2), V(0.040937500000000016, -0.19296875, 1.2), V(-0.009062500000000001, -0.19296875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), V(0.11093750000000002, -0.12, 0.9865560658643532), -1), V(0.11093750000000002, -0.12, 0.9865560658643532), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), 100, 0, null, V(-0.0008999801702201017, -0.00137809463564953, -0.00006642278975471942), V(-0.0014999985603503193, 0, -0.000013854717561505174)),
				new StraightEdge(new L3(V(-0.009062500000000001, -0.19296875, 1.2), V(0, 0, -1)), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), V(-0.009062500000000001, -0.19296875, 0.19999999999999996), 0.21883694901551276, 1),
				new PCurveEdge(new BezierCurve(V(0.11093750000000002, -0.12, 0.19999999999999996), V(0.0809375, -0.16593750000000002, 0.19999999999999996), V(0.040937500000000016, -0.19296875, 0.19999999999999996), V(-0.009062500000000001, -0.19296875, 0.19999999999999996), 0, 1), V(-0.009062500000000001, -0.19296875, 0.19999999999999996), V(0.11093750000000002, -0.12, 0.19999999999999996), 1, 0, null, V(0.15000000000000005, 0, 0), V(0.09000000000000008, 0.1378125, 0)),
				new StraightEdge(new L3(V(0.11093750000000002, -0.12, 1.2), V(0, 0, -1)), V(0.11093750000000002, -0.12, 0.19999999999999996), V(0.11093750000000002, -0.12, 0.9865560658643532), 1, 0.21344393413564677)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 0.19999999999999996), V(-1, 0, 0), V(0, -1, 0)), [
				new PCurveEdge(new BezierCurve(V(0.20296874999999998, 0.11703124999999998, 0.19999999999999996), V(0.20296874999999998, 0.2240625, 0.19999999999999996), V(0.14500000000000002, 0.2890625, 0.19999999999999996), V(0.010937499999999989, 0.2890625, 0.19999999999999996), 0, 1), V(0.20296874999999998, 0.11703124999999998, 0.19999999999999996), V(0.010937499999999989, 0.2890625, 0.19999999999999996), 0, 1, null, V(0, 0.32109375000000007, 0), V(-0.4021875000000001, 0, 0)),
				new PCurveEdge(new BezierCurve(V(0.010937499999999989, 0.2890625, 0.19999999999999996), V(-0.04093749999999999, 0.2890625, 0.19999999999999996), V(-0.1, 0.2790625, 0.19999999999999996), V(-0.16499999999999998, 0.255, 0.19999999999999996), 0, 1), V(0.010937499999999989, 0.2890625, 0.19999999999999996), V(-0.16499999999999998, 0.255, 0.19999999999999996), 0, 1, null, V(-0.15562499999999993, 0, 0), V(-0.19499999999999995, -0.07218749999999996, 0)),
				new StraightEdge(new L3(V(-0.16499999999999998, 0.255, 0.19999999999999996), V(0.33372679389213644, -0.942669839890126, 0)), V(-0.16499999999999998, 0.255, 0.19999999999999996), V(-0.1409375, 0.18703124999999998, 0.19999999999999996), 0, 0.07210239165806154),
				new PCurveEdge(new BezierCurve(V(-0.1409375, 0.18703124999999998, 0.19999999999999996), V(-0.08703125, 0.2059375, 0.19999999999999996), V(-0.037968749999999996, 0.2140625, 0.19999999999999996), V(-0.0029687499999999922, 0.2140625, 0.19999999999999996), 0, 1), V(-0.1409375, 0.18703124999999998, 0.19999999999999996), V(-0.0029687499999999922, 0.2140625, 0.19999999999999996), 0, 1, null, V(0.16171874999999997, 0.05671875000000004, 0), V(0.10500000000000001, 0, 0)),
				new PCurveEdge(new BezierCurve(V(-0.0029687499999999922, 0.2140625, 0.19999999999999996), V(0.07406249999999998, 0.2140625, 0.19999999999999996), V(0.11093750000000002, 0.18593749999999998, 0.19999999999999996), V(0.11093750000000002, 0.11203125000000003, 0.19999999999999996), 0, 1), V(-0.0029687499999999922, 0.2140625, 0.19999999999999996), V(0.11093750000000002, 0.11203125000000003, 0.19999999999999996), 0, 1, null, V(0.2310937499999999, 0, 0), V(0, -0.22171874999999985, 0)),
				new StraightEdge(new L3(V(0.11093750000000002, 0.11203125000000003, 0.19999999999999996), V(0, -1, 0)), V(0.11093750000000002, 0.11203125000000003, 0.19999999999999996), V(0.11093750000000002, 0.07093749999999999, 0.19999999999999996), 0, 0.04109375000000004),
				new StraightEdge(new L3(V(0.11093750000000002, 0.07093749999999999, 0.19999999999999996), V(-1, 0, 0)), V(0.11093750000000002, 0.07093749999999999, 0.19999999999999996), V(0.04296875, 0.07093749999999999, 0.19999999999999996), 0, 0.06796875000000002),
				new PCurveEdge(new BezierCurve(V(0.04296875, 0.07093749999999999, 0.19999999999999996), V(-0.11703125, 0.07093749999999999, 0.19999999999999996), V(-0.20500000000000002, 0.010000000000000009, 0.19999999999999996), V(-0.20500000000000002, -0.0990625, 0.19999999999999996), 0, 1), V(0.04296875, 0.07093749999999999, 0.19999999999999996), V(-0.20500000000000002, -0.0990625, 0.19999999999999996), 0, 1, null, V(-0.48, 0, 0), V(0, -0.3271875, 0)),
				new PCurveEdge(new BezierCurve(V(-0.20500000000000002, -0.0990625, 0.19999999999999996), V(-0.20500000000000002, -0.19703125, 0.19999999999999996), V(-0.14, -0.26203125, 0.19999999999999996), V(-0.034062499999999996, -0.26203125, 0.19999999999999996), 0, 1), V(-0.20500000000000002, -0.0990625, 0.19999999999999996), V(-0.034062499999999996, -0.26203125, 0.19999999999999996), 0, 1, null, V(0, -0.29390625, 0), V(0.3178125, 0, 0)),
				new PCurveEdge(new BezierCurve(V(-0.034062499999999996, -0.26203125, 0.19999999999999996), V(0.030000000000000027, -0.26203125, 0.19999999999999996), V(0.0859375, -0.23796875, 0.19999999999999996), V(0.12296875000000002, -0.18796875, 0.19999999999999996), 0, 1), V(-0.034062499999999996, -0.26203125, 0.19999999999999996), V(0.12296875000000002, -0.18796875, 0.19999999999999996), 0, 1, null, V(0.19218750000000007, 0, 0), V(0.11109375000000005, 0.15000000000000002, 0)),
				new PCurveEdge(new BezierCurve(V(0.12296875000000002, -0.18796875, 0.19999999999999996), V(0.13906249999999998, -0.235, 0.19999999999999996), V(0.17296875, -0.2559375, 0.19999999999999996), V(0.21999999999999997, -0.26203125, 0.19999999999999996), 0, 1), V(0.12296875000000002, -0.18796875, 0.19999999999999996), V(0.21999999999999997, -0.26203125, 0.19999999999999996), 0, 1, null, V(0.04828124999999989, -0.14109375000000002, 0), V(0.1410937499999999, -0.01828125, 0)),
				new StraightEdge(new L3(V(0.21999999999999997, -0.26203125, 0.19999999999999996), V(0.32520868466298514, 0.9456422745519424, 0)), V(0.21999999999999997, -0.26203125, 0.19999999999999996), V(0.24203124999999998, -0.19796875, 0.19999999999999996), 0, 0.06774496204746519),
				new PCurveEdge(new BezierCurve(V(0.24203124999999998, -0.19796875, 0.19999999999999996), V(0.21500000000000002, -0.18796875, 0.19999999999999996), V(0.20296874999999998, -0.17296875, 0.19999999999999996), V(0.20296874999999998, -0.13093749999999998, 0.19999999999999996), 0, 1), V(0.24203124999999998, -0.19796875, 0.19999999999999996), V(0.20296874999999998, -0.13093749999999998, 0.19999999999999996), 0, 1, null, V(-0.08109374999999985, 0.030000000000000006, 0), V(0, 0.12609375, 0)),
				new StraightEdge(new L3(V(0.20296874999999998, -0.13093749999999998, 0.19999999999999996), V(0, 1, 0)), V(0.20296874999999998, -0.13093749999999998, 0.19999999999999996), V(0.20296874999999998, 0.11703124999999998, 0.19999999999999996), 0, 0.24796874999999996)], [[
				new PCurveEdge(new BezierCurve(V(-0.009062500000000001, -0.19296875, 0.19999999999999996), V(-0.07500000000000001, -0.19296875, 0.19999999999999996), V(-0.1040625, -0.16, 0.19999999999999996), V(-0.1040625, -0.095, 0.19999999999999996), 0, 1), V(-0.009062500000000001, -0.19296875, 0.19999999999999996), V(-0.1040625, -0.095, 0.19999999999999996), 0, 1, null, V(-0.19781250000000003, 0, 0), V(0, 0.195, 0)),
				new PCurveEdge(new BezierCurve(V(-0.1040625, -0.095, 0.19999999999999996), V(-0.1040625, -0.030937500000000007, 0.19999999999999996), V(-0.065, 0.010000000000000009, 0.19999999999999996), V(0.047968750000000004, 0.010000000000000009, 0.19999999999999996), 0, 1), V(-0.1040625, -0.095, 0.19999999999999996), V(0.047968750000000004, 0.010000000000000009, 0.19999999999999996), 0, 1, null, V(0, 0.19218749999999998, 0), V(0.33890625, 0, 0)),
				new StraightEdge(new L3(V(0.047968750000000004, 0.010000000000000009, 0.19999999999999996), V(1, 0, 0)), V(0.047968750000000004, 0.010000000000000009, 0.19999999999999996), V(0.11093750000000002, 0.010000000000000009, 0.19999999999999996), 0, 0.06296875000000002),
				new StraightEdge(new L3(V(0.11093750000000002, 0.010000000000000009, 0.19999999999999996), V(0, -1, 0)), V(0.11093750000000002, 0.010000000000000009, 0.19999999999999996), V(0.11093750000000002, -0.12, 0.19999999999999996), 0, 0.13),
				new PCurveEdge(new BezierCurve(V(0.11093750000000002, -0.12, 0.19999999999999996), V(0.0809375, -0.16593750000000002, 0.19999999999999996), V(0.040937500000000016, -0.19296875, 0.19999999999999996), V(-0.009062500000000001, -0.19296875, 0.19999999999999996), 0, 1), V(0.11093750000000002, -0.12, 0.19999999999999996), V(-0.009062500000000001, -0.19296875, 0.19999999999999996), 0, 1, null, V(-0.09000000000000008, -0.1378125, 0), V(-0.15000000000000005, 0, 0))]])], false)
		b2EqualAnd(assert, a, b, result)
	},

	'cylinder(1,2) AND cylinder(1,2).rotateZ(PI/2).translate(0,0,1)'(assert) {
		const a = B2T.cylinder(1,2)
		const b = B2T.cylinder(1,2).rotateZ(PI/2).translate(0,0,1)
		const result = new B2([
			new RotationFace(new SemiCylinderSurface(new SemiEllipseCurve(V(0, 0, 0), V(1, 0, 0), V(0, 1, 0), 0, PI), V3.Z), [
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 1), V(1, 0, 0), V(0, 1, 0), 0, PI), V(0, 1, 1), V(-1, 0, 1), 1.5707963267948966, PI, null, V(-1, 0, 0), V(0, -1, 0)),
				new StraightEdge(new L3(V(-1, 0, 0), V(0, 0, 1)), V(-1, 0, 1), V(-1, 0, 2), 1, 2),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 2), V(1, 0, 0), V(0, 1, 0), 0, PI), V(-1, 0, 2), V(1, 0, 2), PI, 0, null, V(0, 1, 0), V(0, -1, 0)),
				new StraightEdge(new L3(V(1, 0, 0), V(0, 0, 1)), V(1, 0, 2), V(1, 0, 1), 2, 1),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 1), V(1, 0, 0), V(0, 1, 0), 0, PI), V(1, 0, 1), V(0, 1, 1), 0, 1.5707963267948966, null, V(0, 1, 0), V(-1, 0, 0))], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 2), V(1, 0, 0), V(0, 1, 0)), [
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 2), V(0, -1, 0), V(-1, 0, 0), 0, PI), V(0, 1, 2), V(-1, 0, 2), PI, 1.5707963267948968, null, V(-1, 0, 0), V(-2.220446049250313e-16, -1, 0)),
				new StraightEdge(new L3(V(-1, 0, 2), V(1, 0, 0)), V(-1, 0, 2), V(0, 0, 2), 0, 1),
				new StraightEdge(new L3(V(1, 0, 2), V(-1, 0, 0)), V(0, 0, 2), V(1, 0, 2), 1, 0),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 2), V(0, 1, 0), V(1, 0, 0), 0, PI), V(1, 0, 2), V(0, 1, 2), 1.5707963267948963, 1.224646799147353e-16, null, V(-9.957992501029601e-17, 1, 0), V(-1, 3.0616169978683826e-16, 0))], []),
			new RotationFace(new SemiCylinderSurface(new SemiEllipseCurve(V(0, 0, 0), V(-1, 0, 0), V(0, -1, 0), 0, PI), V3.Z), [
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 1), V(-1, 0, 0), V(0, -1, 0), 0, PI), V(-1, 0, 1), V(0, -1, 1), 0, 1.5707963267948966, null, V(0, -1, 0), V(1, 0, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 1), V(-1, 0, 0), V(0, -1, 0), 0, PI), V(0, -1, 1), V(1, 0, 1), 1.5707963267948966, PI, null, V(1, 0, 0), V(2.4492935982947064e-16, 1, 0)),
				new StraightEdge(new L3(V(1, 0, 0), V(0, 0, 1)), V(1, 0, 1), V(1, 0, 2), 1, 2),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 2), V(-1, 0, 0), V(0, -1, 0), 0, PI), V(1, 0, 2), V(-1, 0, 2), PI, 0, null, V(-2.4492935982947064e-16, -1, 0), V(0, 1, 0)),
				new StraightEdge(new L3(V(-1, 0, 0), V(0, 0, 1)), V(-1, 0, 2), V(-1, 0, 1), 2, 1)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 2), V(1, 0, 0), V(0, 1, 0)), [
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 2), V(0, -1, 0), V(-1, 0, 0), 0, PI), V(-1, 0, 2), V(0, -1, 2), 1.5707963267948968, 1.224646799147353e-16, null, V(-2.220446049250313e-16, -1, 0), V(1, 0, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 2), V(0, 1, 0), V(1, 0, 0), 0, PI), V(0, -1, 2), V(1, 0, 2), PI, 1.5707963267948963, null, V(1, -6.123233995736765e-17, 0), V(-9.957992501029601e-17, 1, 0)),
				new StraightEdge(new L3(V(1, 0, 2), V(-1, 0, 0)), V(1, 0, 2), V(0, 0, 2), 0, 1),
				new StraightEdge(new L3(V(-1, 0, 2), V(1, 0, 0)), V(0, 0, 2), V(-1, 0, 2), 1, 0)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), -1), V(-1, 0, 0), V(0, 1, 0)), [
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 1), V(1, 0, 0), V(0, 1, 0), 0, PI), V(-1, 0, 1), V(0, 1, 1), PI, 1.5707963267948966, null, V(0, 1, 0), V(1, 0, 0)),
				new StraightEdge(new L3(V(0, 0, 1), V(0, 1, 0)), V(0, 1, 1), V(0, 0, 1), 1, 0),
				new StraightEdge(new L3(V(0, 0, 1), V(0, -1, 0)), V(0, 0, 1), V(0, -1, 1), 0, 1),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 1), V(-1, 0, 0), V(0, -1, 0), 0, PI), V(0, -1, 1), V(-1, 0, 1), 1.5707963267948966, 0, null, V(-1, 0, 0), V(0, 1, 0))], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), -1), V(-1, 0, 0), V(0, 1, 0)), [
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 1), V(1, 0, 0), V(0, 1, 0), 0, PI), V(0, 1, 1), V(1, 0, 1), 1.5707963267948966, 0, null, V(1, 0, 0), V(0, -1, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 1), V(-1, 0, 0), V(0, -1, 0), 0, PI), V(1, 0, 1), V(0, -1, 1), PI, 1.5707963267948966, null, V(-2.4492935982947064e-16, -1, 0), V(-1, 0, 0)),
				new StraightEdge(new L3(V(0, 0, 1), V(0, -1, 0)), V(0, -1, 1), V(0, 0, 1), 1, 0),
				new StraightEdge(new L3(V(0, 0, 1), V(0, 1, 0)), V(0, 0, 1), V(0, 1, 1), 0, 1)], [])], false)
		b2EqualAnd(assert, a, b, result)
	},

	'box - semicylinder'(assert) {
    	const box = B2T.box(4, 2, 2)
		const cyl = B2T.cylinder(1,2,180 * DEG).translate(2).flipped()
		const result = new B2([
			new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0),V3.X,V3.Z), [
				new StraightEdge(new L3(V(1, 1.2246467991473532e-16, 0), V(0, 0, -1)), V(1, 1.2246467991473532e-16, 0), V(1, 1.2246467991473532e-16, 2), 0, -2),
				new StraightEdge(new L3(V(4, 0, 2), V(-1, 0, 0)), V(1, 0, 2), V(0, 0, 2), 3, 4),
				new StraightEdge(new L3(V3.O, V3.Z), V(0, 0, 2), V3.O, 2, 0),
				new StraightEdge(new L3(V(4, 0, 0), V(-1, 0, 0)), V3.O, V3.X, 4, 3)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0),V3.X,V3.Z), [
				new StraightEdge(new L3(V(3, 0, 0), V3.Z), V(3, 0, 2), V(3, 0, 0), 2, 0),
				new StraightEdge(new L3(V(4, 0, 0), V(-1, 0, 0)), V(3, 0, 0), V(4, 0, 0), 1, 0),
				new StraightEdge(new L3(V(4, 0, 0), V3.Z), V(4, 0, 0), V(4, 0, 2), 0, 2),
				new StraightEdge(new L3(V(4, 0, 2), V(-1, 0, 0)), V(4, 0, 2), V(3, 0, 2), 0, 1)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 0),V3.Y,V3.X), [
				new PCurveEdge(new SemiEllipseCurve(V(2, 0, 0),V(-1, 0, 0),V3.Y,0,3.141592653589793),V(3, 0, 0),V(1, 1.2246467991473532e-16, 0),3.141592653589793,1.2246467991473532e-16,null,V(-1.2246467991473532e-16, 1, 0),V(-1.2246467991473532e-16, -1, 0)),
				new StraightEdge(new L3(V(4, 0, 0), V(-1, 0, 0)), V3.X, V3.O, 3, 4),
				new StraightEdge(new L3(V3.O, V3.Y), V3.O, V(0, 2, 0), 0, 2),
				new StraightEdge(new L3(V(0, 2, 0), V3.X), V(0, 2, 0), V(4, 2, 0), 0, 4),
				new StraightEdge(new L3(V(4, 2, 0), V(0, -1, 0)), V(4, 2, 0), V(4, 0, 0), 0, 2),
				new StraightEdge(new L3(V(4, 0, 0), V(-1, 0, 0)), V(4, 0, 0), V(3, 0, 0), 0, 1)], []),
			new PlaneFace(new PlaneSurface(new P3(V3.Z, 2),V3.Y,V(-1, 0, 0)), [
				new PCurveEdge(new SemiEllipseCurve(V(2, 0, 2),V3.X,V3.Y,0,3.141592653589793),V(1, 1.2246467991473532e-16, 2),V(3, 0, 2),3.141592653589793,0,null,V(1.2246467991473532e-16, 1, 0),V(0, -1, 0)),
				new StraightEdge(new L3(V(4, 0, 2), V(-1, 0, 0)), V(3, 0, 2), V(4, 0, 2), 1, 0),
				new StraightEdge(new L3(V(4, 2, 2), V(0, -1, 0)), V(4, 0, 2), V(4, 2, 2), 2, 0),
				new StraightEdge(new L3(V(0, 2, 2), V3.X), V(4, 2, 2), V(0, 2, 2), 4, 0),
				new StraightEdge(new L3(V(0, 0, 2), V3.Y), V(0, 2, 2), V(0, 0, 2), 2, 0),
				new StraightEdge(new L3(V(4, 0, 2), V(-1, 0, 0)), V(0, 0, 2), V(1, 0, 2), 4, 3)], []),
			new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), 0),V(0, -1, 0),V3.Z), [
				new StraightEdge(new L3(V3.O, V3.Y), V(0, 2, 0), V3.O, 2, 0),
				new StraightEdge(new L3(V3.O, V3.Z), V3.O, V(0, 0, 2), 0, 2),
				new StraightEdge(new L3(V(0, 0, 2), V3.Y), V(0, 0, 2), V(0, 2, 2), 0, 2),
				new StraightEdge(new L3(V(0, 2, 0), V3.Z), V(0, 2, 2), V(0, 2, 0), 2, 0)], []),
			new PlaneFace(new PlaneSurface(new P3(V3.Y, 2),V(-1, 0, 0),V3.Z), [
				new StraightEdge(new L3(V(0, 2, 0), V3.X), V(4, 2, 0), V(0, 2, 0), 4, 0),
				new StraightEdge(new L3(V(0, 2, 0), V3.Z), V(0, 2, 0), V(0, 2, 2), 0, 2),
				new StraightEdge(new L3(V(0, 2, 2), V3.X), V(0, 2, 2), V(4, 2, 2), 0, 4),
				new StraightEdge(new L3(V(4, 2, 0), V3.Z), V(4, 2, 2), V(4, 2, 0), 2, 0)], []),
			new PlaneFace(new PlaneSurface(new P3(V3.X, 4),V3.Y,V3.Z), [
				new StraightEdge(new L3(V(4, 2, 0), V(0, -1, 0)), V(4, 0, 0), V(4, 2, 0), 2, 0),
				new StraightEdge(new L3(V(4, 2, 0), V3.Z), V(4, 2, 0), V(4, 2, 2), 0, 2),
				new StraightEdge(new L3(V(4, 2, 2), V(0, -1, 0)), V(4, 2, 2), V(4, 0, 2), 0, 2),
				new StraightEdge(new L3(V(4, 0, 0), V3.Z), V(4, 0, 2), V(4, 0, 0), 2, 0)], []),
			new RotationFace(new SemiCylinderSurface(new SemiEllipseCurve(V(2, 0, 0),V3.X,V3.Y,0,3.141592653589793),V(0, 0, -1),-Infinity,Infinity), [
				new StraightEdge(new L3(V(1, 1.2246467991473532e-16, 0), V(0, 0, -1)), V(1, 1.2246467991473532e-16, 2), V(1, 1.2246467991473532e-16, 0), -2, 0),
				new PCurveEdge(new SemiEllipseCurve(V(2, 0, 0),V(-1, 0, 0),V3.Y,0,3.141592653589793),V(1, 1.2246467991473532e-16, 0),V(3, 0, 0),1.2246467991473532e-16,3.141592653589793,null,V(1.2246467991473532e-16, 1, 0),V(1.2246467991473532e-16, -1, 0)),
				new StraightEdge(new L3(V(3, 0, 0), V3.Z), V(3, 0, 0), V(3, 0, 2), 0, 2),
				new PCurveEdge(new SemiEllipseCurve(V(2, 0, 2),V3.X,V3.Y,0,3.141592653589793),V(3, 0, 2),V(1, 1.2246467991473532e-16, 2),0,3.141592653589793,null,V3.Y,V(-1.2246467991473532e-16, -1, 0))], [])], false)
		b2EqualAnd(assert, box, cyl, result)
	},

	//async 'B2T.sphere() - "a" - B2T.sphere(0.9)'(assert) {
	//	const a = B2T.sphere(0.9).flipped()
	//	const b = B2T.text('a', await B2T.loadFont('fonts/FiraSansMedium.woff'), 64, 64).scale(0.5/32).translate(-0.25,-0.25,1.2).flipped()
	//	const c = B2T.sphere()
	//	const d = a.and(b)
	//	const result = new B2([
	//		new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.20296874999999998, 0.11703124999999998, 1.2), V(0.20296874999999998, 0.2240625, 1.2), V(0.14500000000000002, 0.2890625, 1.2), V(0.010937499999999989, 0.2890625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.20296874999999998, 0.11703124999999998, 1.2), V(0.20296874999999998, 0.2240625, 1.2), V(0.14500000000000002, 0.2890625, 1.2), V(0.010937499999999989, 0.2890625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), V(0.010937499999999989, 0.2890625, 0.9572477433702835), 1), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), V(0.010937499999999989, 0.2890625, 0.9572477433702835), 0, 100, null, V(0, 0.003208541417521995, -0.0003862503783657504), V(-0.004021833637101672, 0, 0.00004595341771287279)),
	//			new StraightEdge(new L3(V(0.010937499999999989, 0.2890625, 1.2), V(0, 0, -1)), V(0.010937499999999989, 0.2890625, 0.9572477433702835), V(0.010937499999999989, 0.2890625, 0.852245998633904), 0.24275225662971645, 0.3477540013660959),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.20296874999999998, 0.11703124999999998, 1.2), V(0.20296874999999998, 0.2240625, 1.2), V(0.14500000000000002, 0.2890625, 1.2), V(0.010937499999999989, 0.2890625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(0.20296874999999995, 0.11703124999999996, 0.8689691438980299), V(0.010937499999999987, 0.2890625, 0.8522459986339039), -1), V(0.010937499999999989, 0.2890625, 0.852245998633904), V(0.20296874999999998, 0.11703124999999998, 0.8689691438980301), 100, 0, null, V(0.004021823173538629, 0, -0.0000516150161233843), V(0, -0.0032079393598033606, 0.0004320396826956172)),
	//			new StraightEdge(new L3(V(0.20296874999999998, 0.11703124999999998, 1.2), V(0, 0, -1)), V(0.20296874999999998, 0.11703124999999998, 0.8689691438980301), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), 0.3310308561019698, 0.22783367007138378)], []),
	//		new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.010937499999999989, 0.2890625, 1.2), V(-0.04093749999999999, 0.2890625, 1.2), V(-0.1, 0.2790625, 1.2), V(-0.16499999999999998, 0.255, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.010937499999999989, 0.2890625, 1.2), V(-0.04093749999999999, 0.2890625, 1.2), V(-0.1, 0.2790625, 1.2), V(-0.16499999999999998, 0.255, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.010937499999999989, 0.2890625, 0.9572477433702835), V(-0.16499999999999998, 0.255, 0.9527591510974849), 1), V(0.010937499999999989, 0.2890625, 0.9572477433702835), V(-0.16499999999999998, 0.255, 0.9527591510974849), 0, 100, null, V(-0.0015562475396635686, 0, 0.000017781663715540364), V(-0.0019497960757584164, -0.0007217995088144136, -0.0001444829762001187)),
	//			new StraightEdge(new L3(V(-0.16499999999999998, 0.255, 1.2), V(0, 0, -1)), V(-0.16499999999999998, 0.255, 0.9527591510974849), V(-0.16499999999999998, 0.255, 0.8472012747865765), 0.24724084890251508, 0.3527987252134235),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.010937499999999989, 0.2890625, 1.2), V(-0.04093749999999999, 0.2890625, 1.2), V(-0.1, 0.2790625, 1.2), V(-0.16499999999999998, 0.255, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(0.010937499999999987, 0.2890625, 0.8522459986339039), V(-0.16499999999999998, 0.25499999999999995, 0.8472012747865764), -1), V(-0.16499999999999998, 0.255, 0.8472012747865765), V(0.010937499999999989, 0.2890625, 0.852245998633904), 100, 0, null, V(0.001949741977734284, 0.0007217794821420186, 0.0001624804665392168), V(0.0015562468960623505, 6.688077170406378e-21, -0.00001997246153454078)),
	//			new StraightEdge(new L3(V(0.010937499999999989, 0.2890625, 1.2), V(0, 0, -1)), V(0.010937499999999989, 0.2890625, 0.852245998633904), V(0.010937499999999989, 0.2890625, 0.9572477433702835), 0.3477540013660959, 0.24275225662971645)], []),
	//		new PlaneFace(new PlaneSurface(new P3(V(0.942669839890126, 0.33372679389213644, 0), -0.0704401911393759), V(0, 0, -1), V(-0.33372679389213644, 0.942669839890126, 0)), [
	//			new PCurveEdge(new SemiEllipseCurve(V(-0.06640184370318536, -0.0235077791500932, 0), V(0, 0, 0.997516004619599), V(-0.33289781807779234, 0.9403282523625955, 0), 0.02500215049153284, 3.1165905030982604), V(-0.16499999999999998, 0.255, 0.9527591510974849), V(-0.1409375, 0.18703124999999998, 0.9721913045369144), 0.3006922248899664, 0.22581372635962282, null, V(0.31796125684715715, -0.8981373164189183, 0.29544573016418446), V(0.32444628711291806, -0.9164554213903857, 0.22334333850612303)),
	//			new StraightEdge(new L3(V(-0.1409375, 0.18703124999999998, 1.2), V(0, 0, -1)), V(-0.1409375, 0.18703124999999998, 0.9721913045369144), V(-0.1409375, 0.18703124999999998, 0.8689970843548256), 0.22780869546308558, 0.33100291564517437),
	//			new PCurveEdge(new SemiEllipseCurve(V(-0.06640184370318533, -0.02350777915009319, 0), V(0, 0, -0.8972391985820994), V(-0.299432761097154, 0.8458003316705326, 0), 0.027797112240687767, 3.113795541349105), V(-0.1409375, 0.18703124999999998, 0.8689970843548256), V(-0.16499999999999998, 0.255, 0.8472012747865765), 2.8900247146126086, 2.8060483910328458, null, V(-0.2900076108633503, 0.8191773423737494, -0.2233433385061232), V(-0.282733765215855, 0.7986310900577724, -0.2954457301641845)),
	//			new StraightEdge(new L3(V(-0.16499999999999998, 0.255, 1.2), V(0, 0, -1)), V(-0.16499999999999998, 0.255, 0.8472012747865765), V(-0.16499999999999998, 0.255, 0.9527591510974849), 0.3527987252134235, 0.24724084890251508)], []),
	//		new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(-0.1409375, 0.18703124999999998, 1.2), V(-0.08703125, 0.2059375, 1.2), V(-0.037968749999999996, 0.2140625, 1.2), V(-0.0029687499999999922, 0.2140625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1409375, 0.18703124999999998, 1.2), V(-0.08703125, 0.2059375, 1.2), V(-0.037968749999999996, 0.2140625, 1.2), V(-0.0029687499999999922, 0.2140625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(-0.1409375, 0.18703124999999998, 0.9721913045369145), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), 1), V(-0.1409375, 0.18703124999999998, 0.9721913045369144), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), 0, 100, null, V(0.001617060512887581, 0.0005671429624910072, 0.00012531587997459948), V(0.0010499999469858546, 0, 0.0000031911732431043194)),
	//			new StraightEdge(new L3(V(-0.0029687499999999922, 0.2140625, 1.2), V(0, 0, -1)), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), V(-0.0029687499999999922, 0.2140625, 0.874167279539327), 0.22318454526088316, 0.3258327204606729),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1409375, 0.18703124999999998, 1.2), V(-0.08703125, 0.2059375, 1.2), V(-0.037968749999999996, 0.2140625, 1.2), V(-0.0029687499999999922, 0.2140625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(-0.14093749999999997, 0.18703124999999995, 0.8689970843548255), V(-0.002968749999999992, 0.21406249999999996, 0.8741672795393268), -1), V(-0.0029687499999999922, 0.2140625, 0.874167279539327), V(-0.1409375, 0.18703124999999998, 0.8689970843548256), 100, 0, null, V(-0.001049999933945326, 0, -0.0000035658933671629643), V(-0.0016170285671746015, -0.0005671317583423969, -0.00014019448879870969)),
	//			new StraightEdge(new L3(V(-0.1409375, 0.18703124999999998, 1.2), V(0, 0, -1)), V(-0.1409375, 0.18703124999999998, 0.8689970843548256), V(-0.1409375, 0.18703124999999998, 0.9721913045369144), 0.33100291564517437, 0.22780869546308558)], []),
	//		new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(-0.0029687499999999922, 0.2140625, 1.2), V(0.07406249999999998, 0.2140625, 1.2), V(0.11093750000000002, 0.18593749999999998, 1.2), V(0.11093750000000002, 0.11203125000000003, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.0029687499999999922, 0.2140625, 1.2), V(0.07406249999999998, 0.2140625, 1.2), V(0.11093750000000002, 0.18593749999999998, 1.2), V(0.11093750000000002, 0.11203125000000003, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), 1), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), 0, 100, null, V(0.0023109369300238886, 0, 0.000007023429018985676), V(0, -0.002216485816645523, 0.00025146076711587485)),
	//			new StraightEdge(new L3(V(0.11093750000000002, 0.11203125000000003, 1.2), V(0, 0, -1)), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), V(0.11093750000000002, 0.11203125000000003, 0.8860823156553729), 0.21250728098016458, 0.3139176843446271),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.0029687499999999922, 0.2140625, 1.2), V(0.07406249999999998, 0.2140625, 1.2), V(0.11093750000000002, 0.18593749999999998, 1.2), V(0.11093750000000002, 0.11203125000000003, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(-0.002968749999999992, 0.21406249999999996, 0.8741672795393268), V(0.11093750000000001, 0.11203125000000001, 0.8860823156553727), -1), V(0.11093750000000002, 0.11203125000000003, 0.8860823156553729), V(-0.0029687499999999922, 0.2140625, 0.874167279539327), 100, 0, null, V(0, 0.0022163159641569587, -0.0002802184892673451), V(-0.0023109367883072116, 3.1668370133166326e-21, -0.00000784814731787084)),
	//			new StraightEdge(new L3(V(-0.0029687499999999922, 0.2140625, 1.2), V(0, 0, -1)), V(-0.0029687499999999922, 0.2140625, 0.874167279539327), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), 0.3258327204606729, 0.22318454526088316)], []),
	//		new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 0.11093750000000002), V(0, 0, -1), V(0, 1, 0)), [
	//			new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000002, 0, 0), V(0, 0, 0.9938273849586506), V(0, 0.9938273849586506, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), V(0.11093750000000002, 0.07093749999999999, 0.9912924604714292), 0.11296719097084001, 0.07143883874192354, null, V(0, -0.9874927190198353, 0.11203125000000001), V(0, -0.9912924604714292, 0.07093749999999999)),
	//			new StraightEdge(new L3(V(0.11093750000000002, 0.07093749999999999, 1.2), V(0, 0, -1)), V(0.11093750000000002, 0.07093749999999999, 0.9912924604714292), V(0.11093750000000002, 0.07093749999999999, 0.8903149679677975), 0.2087075395285708, 0.30968503203220243),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000001, 0, 0), V(0, 0, -0.8931365355273235), V(0, 0.8931365355273235, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.07093749999999999, 0.8903149679677975), V(0.11093750000000002, 0.11203125000000003, 0.8860823156553729), 3.0620837623516257, 3.015825616877026, null, V(0, 0.8903149679677973, -0.0709375), V(0, 0.8860823156553727, -0.11203125000000015)),
	//			new StraightEdge(new L3(V(0.11093750000000002, 0.11203125000000003, 1.2), V(0, 0, -1)), V(0.11093750000000002, 0.11203125000000003, 0.8860823156553729), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), 0.3139176843446271, 0.21250728098016458)], []),
	//		new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), -0.07093749999999999), V(0, 0, -1), V(1, 0, 0)), [
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0.07093749999999999, 0), V(0, 0, 0.9974807622674986), V(0.9974807622674986, 0, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.07093749999999999, 0.9912924604714292), V(0.04296875, 0.07093749999999999, 0.9965548442595558), 0.11144825166905915, 0.04309060575909555, null, V(-0.9912924604714292, 0, 0.11093750000000004), V(-0.9965548442595558, 0, 0.04296875)),
	//			new StraightEdge(new L3(V(0.04296875, 0.07093749999999999, 1.2), V(0, 0, -1)), V(0.04296875, 0.07093749999999999, 0.9965548442595558), V(0.04296875, 0.07093749999999999, 0.8961704958417165), 0.2034451557404442, 0.30382950415828347),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0.07093749999999997, 0), V(0, 0, -0.8972000173282154), V(0.8972000173282154, 0, 0), 0, 3.141592653589793), V(0.04296875, 0.07093749999999999, 0.8961704958417165), V(0.11093750000000002, 0.07093749999999999, 0.8903149679677975), 3.093682274615096, 3.0176268184353994, null, V(0.8961704958417163, 0, -0.04296875000000003), V(0.8903149679677973, 0, -0.11093750000000019)),
	//			new StraightEdge(new L3(V(0.11093750000000002, 0.07093749999999999, 1.2), V(0, 0, -1)), V(0.11093750000000002, 0.07093749999999999, 0.8903149679677975), V(0.11093750000000002, 0.07093749999999999, 0.9912924604714292), 0.30968503203220243, 0.2087075395285708)], []),
	//		new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.04296875, 0.07093749999999999, 0.9965548442595558), V(-0.1710469396532132, 6.938893903907228e-18, 0.9852628808776215), 1), V(0.04296875, 0.07093749999999999, 0.9965548442595558), V(-0.171046939653212, 2.0947208715025655e-17, 0.9852628808776216), 0, 64, null, V(-0.0047989723215347765, 0, 0.0002069187091194757), V(-0.002384621933525827, -0.0028002649734856517, -0.00041398320374781613)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.17104693965321316, 2.0947208715025797e-17, 0.9852628808776215), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), 1), V(-0.171046939653212, 2.0947208715025655e-17, 0.9852628808776216), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), 0, 36, null, V(-0.0018504236016956069, -0.002172955102479456, -0.000321243497826874), V(0, -0.003269989649016917, -0.0003326706415016505)),
	//			new StraightEdge(new L3(V(-0.20500000000000002, -0.0990625, 1.2), V(0, 0, -1)), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), V(-0.20500000000000002, -0.0990625, 0.8707247677043248), 0.22626409068282272, 0.3292752322956751),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.17104693965321313, 2.0947208715025794e-17, 0.883596595984429), V(-0.205, -0.09906249999999998, 0.8707247677043246), -1), V(-0.20500000000000002, -0.0990625, 0.8707247677043248), V(-0.17104693965321202, 2.0947208715025658e-17, 0.8835965959844292), 36, 0, null, V(0, 0.0032695078711222625, 0.0003719724481215541), V(0.00215411320093391, 0.002529578236571627, 0.00041699399065229685)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(0.04296874999999999, 0.07093749999999997, 0.8961704958417162), V(-0.17104693965321316, 3.122502256758253e-18, 0.8835965959844289), -1), V(-0.17104693965321202, 2.0947208715025658e-17, 0.8835965959844292), V(0.04296875, 0.07093749999999999, 0.8961704958417165), 64, 0, null, V(0.0023970206921961742, 0.0028148248536625227, 0.00046401610819787214), V(0.004798729293096285, 5.284086361068957e-20, -0.000230085012025602)),
	//			new StraightEdge(new L3(V(0.04296875, 0.07093749999999999, 1.2), V(0, 0, -1)), V(0.04296875, 0.07093749999999999, 0.8961704958417165), V(0.04296875, 0.07093749999999999, 0.9965548442595558), 0.30382950415828347, 0.2034451557404442)], []),
	//		new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(-0.20500000000000002, -0.0990625, 1.2), V(-0.20500000000000002, -0.19703125, 1.2), V(-0.14, -0.26203125, 1.2), V(-0.034062499999999996, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.20500000000000002, -0.0990625, 1.2), V(-0.20500000000000002, -0.19703125, 1.2), V(-0.14, -0.26203125, 1.2), V(-0.034062499999999996, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), V(-0.034062499999999996, -0.26203125, 0.9644580706890205), 1), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), V(-0.034062499999999996, -0.26203125, 0.9644580706890203), 0, 100, null, V(0, -0.002937749569980493, -0.00029887037541859586), V(0.003177924394426102, 0, 0.00011223717544071708)),
	//			new StraightEdge(new L3(V(-0.034062499999999996, -0.26203125, 1.2), V(0, 0, -1)), V(-0.034062499999999996, -0.26203125, 0.9644580706890203), V(-0.034062499999999996, -0.26203125, 0.8603367771501972), 0.23554192931097961, 0.33966322284980277),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.20500000000000002, -0.0990625, 1.2), V(-0.20500000000000002, -0.19703125, 1.2), V(-0.14, -0.26203125, 1.2), V(-0.034062499999999996, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.205, -0.09906249999999998, 0.8707247677043246), V(-0.034062499999999996, -0.26203124999999994, 0.8603367771501969), -1), V(-0.034062499999999996, -0.26203125, 0.8603367771501972), V(-0.20500000000000002, -0.0990625, 0.8707247677043248), 100, 0, null, V(-0.003177872773799387, 0, -0.0001258185099515327), V(0, 0.0029374208172160575, 0.0003341908493916618)),
	//			new StraightEdge(new L3(V(-0.20500000000000002, -0.0990625, 1.2), V(0, 0, -1)), V(-0.20500000000000002, -0.0990625, 0.8707247677043248), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), 0.3292752322956751, 0.22626409068282272)], []),
	//		new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(-0.034062499999999996, -0.26203125, 1.2), V(0.030000000000000027, -0.26203125, 1.2), V(0.0859375, -0.23796875, 1.2), V(0.12296875000000002, -0.18796875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.034062499999999996, -0.26203125, 1.2), V(0.030000000000000027, -0.26203125, 1.2), V(0.0859375, -0.23796875, 1.2), V(0.12296875000000002, -0.18796875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.034062499999999996, -0.26203125, 0.9644580706890205), V(0.12296875000000002, -0.18796875, 0.9744467330474638), 1), V(-0.034062499999999996, -0.26203125, 0.9644580706890203), V(0.12296875000000002, -0.18796875, 0.9744467330474638), 0, 100, null, V(0.0019218307292635764, 0, 0.00006787475910567421), V(0.0011108137429678576, 0.0014998329018975289, 0.00014913728313140515)),
	//			new StraightEdge(new L3(V(0.12296875000000002, -0.18796875, 1.2), V(0, 0, -1)), V(0.12296875000000002, -0.18796875, 0.9744467330474638), V(0.12296875000000002, -0.18796875, 0.871519612829726), 0.2255532669525362, 0.32848038717027395),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.034062499999999996, -0.26203125, 1.2), V(0.030000000000000027, -0.26203125, 1.2), V(0.0859375, -0.23796875, 1.2), V(0.12296875000000002, -0.18796875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.034062499999999996, -0.26203124999999994, 0.8603367771501969), V(0.12296874999999999, -0.18796874999999996, 0.8715196128297258), -1), V(0.12296875000000002, -0.18796875, 0.871519612829726), V(-0.034062499999999996, -0.26203125, 0.8603367771501972), 100, 0, null, V(-0.0011107827301533816, -0.0014997910280551989, -0.0001667458526657328), V(-0.0019218193657008743, 0, -0.00007608877579431633)),
	//			new StraightEdge(new L3(V(-0.034062499999999996, -0.26203125, 1.2), V(0, 0, -1)), V(-0.034062499999999996, -0.26203125, 0.8603367771501972), V(-0.034062499999999996, -0.26203125, 0.9644580706890203), 0.33966322284980277, 0.23554192931097961)], []),
	//		new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.12296875000000002, -0.18796875, 1.2), V(0.13906249999999998, -0.235, 1.2), V(0.17296875, -0.2559375, 1.2), V(0.21999999999999997, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.12296875000000002, -0.18796875, 1.2), V(0.13906249999999998, -0.235, 1.2), V(0.17296875, -0.2559375, 1.2), V(0.21999999999999997, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.12296875000000002, -0.18796875, 0.9744467330474638), V(0.21999999999999997, -0.26203125, 0.9396486705271484), 1), V(0.12296875000000002, -0.18796875, 0.9744467330474638), V(0.21999999999999997, -0.26203125, 0.9396486705271484), 0, 100, null, V(0.0004825448774820869, -0.0014101554186612473, -0.00033291003063843776), V(0.0014098883430614425, -0.00018267656272224708, -0.0003810385889743954)),
	//			new StraightEdge(new L3(V(0.21999999999999997, -0.26203125, 1.2), V(0, 0, -1)), V(0.21999999999999997, -0.26203125, 0.9396486705271484), V(0.21999999999999997, -0.26203125, 0.8324299514214022), 0.26035132947285156, 0.36757004857859776),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.12296875000000002, -0.18796875, 1.2), V(0.13906249999999998, -0.235, 1.2), V(0.17296875, -0.2559375, 1.2), V(0.21999999999999997, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(0.12296874999999999, -0.18796874999999996, 0.8715196128297258), V(0.21999999999999995, -0.2620312499999999, 0.832429951421402), -1), V(0.21999999999999997, -0.26203125, 0.8324299514214022), V(0.12296875000000002, -0.18796875, 0.871519612829726), 100, 0, null, V(-0.0014095926581339242, 0.00018263825138612356, 0.00043002695120080944), V(-0.0004824780014564421, 0.0014099599848387352, 0.0003721753680202314)),
	//			new StraightEdge(new L3(V(0.12296875000000002, -0.18796875, 1.2), V(0, 0, -1)), V(0.12296875000000002, -0.18796875, 0.871519612829726), V(0.12296875000000002, -0.18796875, 0.9744467330474638), 0.32848038717027395, 0.2255532669525362)], []),
	//		new PlaneFace(new PlaneSurface(new P3(V(-0.9456422745519422, 0.325208684662986, 0), -0.2932561385545256), V(0, 0, -1), V(-0.325208684662986, -0.9456422745519422, 0)), [
	//			new PCurveEdge(new SemiEllipseCurve(V(0.27731540188902115, -0.09536944308866366, 0), V(0.31091053038645644, 0.9040660812655796, 0), V(0, 0, 0.9560339100680944), 1.465110231866663, 3.141592653589793), V(0.21999999999999997, -0.26203125, 0.9396486705271484), V(0.24203124999999998, -0.19796875, 0.9498574882827818), 1.7562036900306723, 1.684527864719763, null, V(0.30558190818745745, 0.8885715060770011, 0.17624191662783179), V(0.30890190438173115, 0.8982253957199245, 0.10849695458036647)),
	//			new StraightEdge(new L3(V(0.24203124999999998, -0.19796875, 1.2), V(0, 0, -1)), V(0.24203124999999998, -0.19796875, 0.9498574882827818), V(0.24203124999999998, -0.19796875, 0.8439367559520533), 0.25014251171721813, 0.35606324404794665),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.2773154018890211, -0.09536944308866366, 0), V(-0.27671434201165657, -0.8046303562041052, 0), V(0, 0, 0.8508823874073836), 0, 1.6896013946143438), V(0.24203124999999998, -0.19796875, 0.8439367559520533), V(0.21999999999999997, -0.26203125, 0.8324299514214022), 1.4429371319818958, 1.3621575275257576, null, V(-0.2744555623419146, -0.7980622734764866, -0.10849695458036643), V(-0.27071344957582744, -0.7871809526672972, -0.17624191662783154)),
	//			new StraightEdge(new L3(V(0.21999999999999997, -0.26203125, 1.2), V(0, 0, -1)), V(0.21999999999999997, -0.26203125, 0.8324299514214022), V(0.21999999999999997, -0.26203125, 0.9396486705271484), 0.36757004857859776, 0.26035132947285156)], []),
	//		new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.24203124999999998, -0.19796875, 1.2), V(0.21500000000000002, -0.18796875, 1.2), V(0.20296874999999998, -0.17296875, 1.2), V(0.20296874999999998, -0.13093749999999998, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.24203124999999998, -0.19796875, 1.2), V(0.21500000000000002, -0.18796875, 1.2), V(0.20296874999999998, -0.17296875, 1.2), V(0.20296874999999998, -0.13093749999999998, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.24203124999999998, -0.19796875, 0.9498574882827818), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), 1), V(0.24203124999999998, -0.19796875, 0.9498574882827818), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), 0, 100, null, V(-0.0008106439109762247, 0.00029989138903166777, 0.0002690617125763413), V(0, 0.0012607577897218358, 0.00017011744865842043)),
	//			new StraightEdge(new L3(V(0.20296874999999998, -0.13093749999999998, 1.2), V(0, 0, -1)), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), V(0.20296874999999998, -0.13093749999999998, 0.8669827320178805), 0.22960881206742834, 0.3330172679821195),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.24203124999999998, -0.19796875, 1.2), V(0.21500000000000002, -0.18796875, 1.2), V(0.20296874999999998, -0.17296875, 1.2), V(0.20296874999999998, -0.13093749999999998, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(0.2420312499999999, -0.19796874999999997, 0.843936755952053), V(0.20296874999999995, -0.13093749999999996, 0.8669827320178802), -1), V(0.20296874999999998, -0.13093749999999998, 0.8669827320178805), V(0.24203124999999998, -0.19796875, 0.8439367559520533), 100, 0, null, V(-1.6025994766705188e-19, -0.0012607132877938707, -0.00019040130792028905), V(0.0008105656446711983, -0.000299862435022873, -0.00030280184601271074)),
	//			new StraightEdge(new L3(V(0.24203124999999998, -0.19796875, 1.2), V(0, 0, -1)), V(0.24203124999999998, -0.19796875, 0.8439367559520533), V(0.24203124999999998, -0.19796875, 0.9498574882827818), 0.35606324404794665, 0.25014251171721813)], []),
	//		new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), -0.20296874999999998), V(0, 0, -1), V(0, -1, 0)), [
	//			new PCurveEdge(new SemiEllipseCurve(V(0.20296874999999998, 0, 0), V(0, 0, -0.9791852156376941), V(0, 0.9791852156376941, 0), 0, 3.141592653589793), V(0.2029687500000001, 0, 0.9791852156376941), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), 3.141592653589793, 3.0217872455230754, null, V(0, 0.9791852156376941, -1.199156040103113e-16), V(0, 0.9721663299286162, -0.1170312499999999)),
	//			new StraightEdge(new L3(V(0.20296874999999998, 0.11703124999999998, 1.2), V(0, 0, -1)), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), V(0.20296874999999998, 0.11703124999999998, 0.8689691438980301), 0.22783367007138378, 0.3310308561019698),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.20296874999999995, 0, 0), V(0, 0, 0.8768145108992196), V(0, 0.8768145108992196, 0), 0, 3.141592653589793), V(0.20296874999999998, 0.11703124999999998, 0.8689691438980301), V(0.20296875000000006, 0, 0.8768145108992196), 0.13387273064879698, 0, null, V(0, -0.86896914389803, 0.11703124999999995), V(0, -0.8768145108992196, 0)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.20296874999999995, 4.943121956121308e-19, 0), V(0, 0, -0.8768145108992196), V(2.1354031397797386e-18, -0.8768145108992196, 0), 2.891240380027181e-17, 3.141592653589793), V(0.20296875000000006, 0, 0.8768145108992196), V(0.20296874999999998, -0.13093749999999998, 0.8669827320178805), 3.141592653589793, 2.9916987961562977, null, V(2.1354031397797386e-18, -0.8768145108992196, -1.0737880842186813e-16), V(2.111458723678207e-18, -0.8669827320178803, -0.13093749999999996)),
	//			new StraightEdge(new L3(V(0.20296874999999998, -0.13093749999999998, 1.2), V(0, 0, -1)), V(0.20296874999999998, -0.13093749999999998, 0.8669827320178805), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), 0.3330172679821195, 0.22960881206742834),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.20296874999999998, -2.3224450485052117e-18, 0), V(0, 0, 0.9791852156376941), V(-1.1204206832959602e-17, -0.9791852156376941, 0), 2.3013070043406876e-17, 3.141592653589793), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), V(0.2029687500000001, 0, 0.9791852156376941), 0.13412262887088133, 2.3013070043406876e-17, null, V(1.1103582248632314e-17, 0.9703911879325716, 0.13093749999999996), V(1.1204206832959602e-17, 0.9791852156376941, 2.3224450485052132e-18))], []),
	//		new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(-0.009062500000000001, -0.19296875, 1.2), V(-0.07500000000000001, -0.19296875, 1.2), V(-0.1040625, -0.16, 1.2), V(-0.1040625, -0.095, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.009062500000000001, -0.19296875, 1.2), V(-0.07500000000000001, -0.19296875, 1.2), V(-0.1040625, -0.16, 1.2), V(-0.1040625, -0.095, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), V(-0.1040625, -0.095, 0.9900232300778351), 1), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), V(-0.1040625, -0.095, 0.9900232300778351), 0, 100, null, V(-0.001978121698253082, 0, -0.000018270895823513832), V(0, 0.0019496588340396214, 0.00018708408409689773)),
	//			new StraightEdge(new L3(V(-0.1040625, -0.095, 1.2), V(0, 0, -1)), V(-0.1040625, -0.095, 0.9900232300778351), V(-0.1040625, -0.095, 0.8889015671567637), 0.20997676992216485, 0.31109843284323624),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.009062500000000001, -0.19296875, 1.2), V(-0.07500000000000001, -0.19296875, 1.2), V(-0.1040625, -0.16, 1.2), V(-0.1040625, -0.095, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.0090625, -0.19296874999999997, 0.879022714505824), V(-0.10406249999999999, -0.09499999999999999, 0.8889015671567634), -1), V(-0.1040625, -0.095, 0.8889015671567637), V(-0.009062500000000001, -0.19296875, 0.8790227145058244), 100, 0, null, V(0, -0.0019495768662389622, -0.0002083580557576388), V(0.001978120886365625, 1.734719868511089e-20, 0.000020393921837124158)),
	//			new StraightEdge(new L3(V(-0.009062500000000001, -0.19296875, 1.2), V(0, 0, -1)), V(-0.009062500000000001, -0.19296875, 0.8790227145058244), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), 0.3209772854941756, 0.21883694901551276)], []),
	//		new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(-0.03286362438479655, -2.6020852139652106e-18, 0.9994598452125503), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), 1), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), 0, 29, null, V(0.0022900449065656465, 0.0007030314524671974, 0.00007529984920773373), V(0.003388522845144273, 0, -0.0001627386905160511)),
	//			new StraightEdge(new L3(V(0.047968750000000004, 0.010000000000000009, 1.2), V(0, 0, -1)), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), V(0.047968750000000004, 0.010000000000000009, 0.8986651206224918), 0.20120122195537427, 0.30133487937750814),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(-0.03286362438479655, -2.3418766925686897e-18, 0.8993997899667839), V(0.04796875, 0.010000000000000007, 0.8986651206224917), -1), V(0.047968750000000004, 0.010000000000000009, 0.8986651206224918), V(-0.03286362438479702, 4.024633241122258e-18, 0.8993997899667838), 29, 0, null, V(-0.003388395809325897, 0, 0.00018086504944802423), V(-0.002292032896076846, -0.0007036417545401136, -0.00008374975068322923)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.10406249999999999, -0.09499999999999999, 0.8889015671567634), V(-0.03286362438479653, 4.0246332411221976e-18, 0.8993997899667839), -1), V(-0.03286362438479702, 4.024633241122258e-18, 0.8993997899667838), V(-0.1040625, -0.095, 0.8889015671567637), 71, 0, null, V(-0.0022509978133461728, -0.0006910442051507778, -0.00008225034901502365), V(0, -0.0019214697267915689, -0.00020535414807408821)),
	//			new StraightEdge(new L3(V(-0.1040625, -0.095, 1.2), V(0, 0, -1)), V(-0.1040625, -0.095, 0.8889015671567637), V(-0.1040625, -0.095, 0.9900232300778351), 0.31109843284323624, 0.20997676992216485),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.1040625, -0.095, 0.9900232300778351), V(-0.03286362438479653, 2.7235906341395923e-18, 0.9994598452125503), 1), V(-0.1040625, -0.095, 0.9900232300778351), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), 0, 71, null, V(0, 0.0019215482684322, 0.00018438666887312053), V(0.0022530279161436587, 0.0006916674357757679, 0.00007408267927848543))], []),
	//		new PlaneFace(new PlaneSurface(new P3(V(0, 1, 0), 0.010000000000000009), V(0, 0, -1), V(-1, 0, 0)), [
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0.010000000000000009, 0), V(0, 0, -0.9999499987499375), V(0.9999499987499375, 0, 0), 0, 3.141592653589793), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), 3.0936030871101416, 3.0304207486077566, null, V(0.9987987780446257, 0, -0.0479687500000002), V(0.9937770731375071, 0, -0.11093749999999983)),
	//			new StraightEdge(new L3(V(0.11093750000000002, 0.010000000000000009, 1.2), V(0, 0, -1)), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), V(0.11093750000000002, 0.010000000000000009, 0.8930805512907278), 0.20622292686249288, 0.3069194487092721),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0.010000000000000007, 0), V(0, 0, 0.899944442729661), V(0.899944442729661, 0, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.010000000000000009, 0.8930805512907278), V(0.047968750000000004, 0.010000000000000009, 0.8986651206224918), 0.12358585566041683, 0.053327173050216115, null, V(-0.8930805512907277, 0, 0.11093750000000001), V(-0.8986651206224918, 0, 0.04796875)),
	//			new StraightEdge(new L3(V(0.047968750000000004, 0.010000000000000009, 1.2), V(0, 0, -1)), V(0.047968750000000004, 0.010000000000000009, 0.8986651206224918), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), 0.30133487937750814, 0.20120122195537427)], []),
	//		new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 0.11093750000000002), V(0, 0, -1), V(0, 1, 0)), [
	//			new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000002, 0, 0), V(0, 0, 0.9938273849586506), V(0, 0.9938273849586506, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), V(0.11093749999999993, 0, 0.9938273849586506), 0.010062279327831384, 0, null, V(0, -0.993777073137507, 0.010000000000000007), V(0, -0.9938273849586506, 0)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000002, 2.7017833632379675e-19, 0), V(0, 0, -0.9938273849586506), V(2.4203775050019737e-18, -0.9938273849586506, 0), 1.3942163371701867e-17, 3.141592653589793), V(0.11093749999999993, 0, 0.9938273849586506), V(0.11093750000000002, -0.12, 0.9865560658643532), 3.141592653589793, 3.02055199779192, null, V(2.4203775050019737e-18, -0.9938273849586506, -1.217087525894596e-16), V(2.4026688591808874e-18, -0.9865560658643532, -0.11999999999999994)),
	//			new StraightEdge(new L3(V(0.11093750000000002, -0.12, 1.2), V(0, 0, -1)), V(0.11093750000000002, -0.12, 0.9865560658643532), V(0.11093750000000002, -0.12, 0.8850383444200315), 0.21344393413564677, 0.31496165557996847),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000001, 2.701783363237952e-19, 0), V(0, 0, 0.8931365355273235), V(2.175153967583285e-18, -0.8931365355273235, 0), 1.5513981584219775e-17, 3.141592653589793), V(0.11093750000000002, -0.12, 0.8850383444200315), V(0.11093750000000006, 0, 0.8931365355273235), 0.13476551593405545, 1.5513981584219775e-17, null, V(-2.155431549099e-18, 0.8850383444200314, 0.11999999999999998), V(-2.175153967583285e-18, 0.8931365355273235, 2.7017833632379536e-19)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000001, 0, 0), V(0, 0, -0.8931365355273235), V(0, 0.8931365355273235, 0), 0, 3.141592653589793), V(0.11093750000000006, 0, 0.8931365355273235), V(0.11093750000000002, 0.010000000000000009, 0.8930805512907278), 3.141592653589793, 3.130395923246913, null, V(0, 0.8931365355273235, -1.093776799435093e-16), V(0, 0.8930805512907277, -0.01000000000000018)),
	//			new StraightEdge(new L3(V(0.11093750000000002, 0.010000000000000009, 1.2), V(0, 0, -1)), V(0.11093750000000002, 0.010000000000000009, 0.8930805512907278), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), 0.3069194487092721, 0.20622292686249288)], []),
	//		new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.11093750000000002, -0.12, 1.2), V(0.0809375, -0.16593750000000002, 1.2), V(0.040937500000000016, -0.19296875, 1.2), V(-0.009062500000000001, -0.19296875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.11093750000000002, -0.12, 1.2), V(0.0809375, -0.16593750000000002, 1.2), V(0.040937500000000016, -0.19296875, 1.2), V(-0.009062500000000001, -0.19296875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.11093750000000002, -0.12, 0.9865560658643532), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), 1), V(0.11093750000000002, -0.12, 0.9865560658643532), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), 0, 100, null, V(-0.0008999801458405889, -0.001378094598318401, -0.00006642278795539617), V(-0.001499998572758444, 0, -0.000013854717676112659)),
	//			new StraightEdge(new L3(V(-0.009062500000000001, -0.19296875, 1.2), V(0, 0, -1)), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), V(-0.009062500000000001, -0.19296875, 0.8790227145058244), 0.21883694901551276, 0.3209772854941756),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.11093750000000002, -0.12, 1.2), V(0.0809375, -0.16593750000000002, 1.2), V(0.040937500000000016, -0.19296875, 1.2), V(-0.009062500000000001, -0.19296875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(0.11093750000000001, -0.11999999999999997, 0.8850383444200313), V(-0.0090625, -0.19296874999999997, 0.879022714505824), -1), V(-0.009062500000000001, -0.19296875, 0.8790227145058244), V(0.11093750000000002, -0.12, 0.8850383444200315), 100, 0, null, V(0.0014999982255756404, 0, 0.000015464599145110234), V(0.0008999753301001092, 0.0013780872242157918, 0.00007404137248522919)),
	//			new StraightEdge(new L3(V(0.11093750000000002, -0.12, 1.2), V(0, 0, -1)), V(0.11093750000000002, -0.12, 0.8850383444200315), V(0.11093750000000002, -0.12, 0.9865560658643532), 0.31496165557996847, 0.21344393413564677)], []),
	//		new RotationFace(new SemiEllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.20296874999999998, 0.11703124999999998, 1.2), V(0.20296874999999998, 0.2240625, 1.2), V(0.14500000000000002, 0.2890625, 1.2), V(0.010937499999999989, 0.2890625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(0.20296874999999995, 0.11703124999999996, 0.8689691438980299), V(0.010937499999999987, 0.2890625, 0.8522459986339039), -1), V(0.20296874999999998, 0.11703124999999998, 0.8689691438980301), V(0.010937499999999989, 0.2890625, 0.852245998633904), 0, 100, null, V(0, 0.0032079393598033606, -0.0004320396826956172), V(-0.004021823173538629, 0, 0.0000516150161233843)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.010937499999999989, 0.2890625, 1.2), V(-0.04093749999999999, 0.2890625, 1.2), V(-0.1, 0.2790625, 1.2), V(-0.16499999999999998, 0.255, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(0.010937499999999987, 0.2890625, 0.8522459986339039), V(-0.16499999999999998, 0.25499999999999995, 0.8472012747865764), -1), V(0.010937499999999989, 0.2890625, 0.852245998633904), V(-0.16499999999999998, 0.255, 0.8472012747865765), 0, 100, null, V(-0.0015562468960623505, -6.688077170406378e-21, 0.00001997246153454078), V(-0.001949741977734284, -0.0007217794821420186, -0.0001624804665392168)),
	//			new PCurveEdge(new SemiEllipseCurve(V(-0.06640184370318533, -0.02350777915009319, 0), V(0, 0, -0.8972391985820994), V(-0.299432761097154, 0.8458003316705326, 0), 0.027797112240687767, 3.113795541349105), V(-0.16499999999999998, 0.255, 0.8472012747865765), V(-0.1409375, 0.18703124999999998, 0.8689970843548256), 2.8060483910328458, 2.8900247146126086, null, V(0.282733765215855, -0.7986310900577724, 0.2954457301641845), V(0.2900076108633503, -0.8191773423737494, 0.2233433385061232)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1409375, 0.18703124999999998, 1.2), V(-0.08703125, 0.2059375, 1.2), V(-0.037968749999999996, 0.2140625, 1.2), V(-0.0029687499999999922, 0.2140625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(-0.14093749999999997, 0.18703124999999995, 0.8689970843548255), V(-0.002968749999999992, 0.21406249999999996, 0.8741672795393268), -1), V(-0.1409375, 0.18703124999999998, 0.8689970843548256), V(-0.0029687499999999922, 0.2140625, 0.874167279539327), 0, 100, null, V(0.0016170285671746015, 0.0005671317583423969, 0.00014019448879870969), V(0.001049999933945326, 0, 0.0000035658933671629643)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.0029687499999999922, 0.2140625, 1.2), V(0.07406249999999998, 0.2140625, 1.2), V(0.11093750000000002, 0.18593749999999998, 1.2), V(0.11093750000000002, 0.11203125000000003, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(-0.002968749999999992, 0.21406249999999996, 0.8741672795393268), V(0.11093750000000001, 0.11203125000000001, 0.8860823156553727), -1), V(-0.0029687499999999922, 0.2140625, 0.874167279539327), V(0.11093750000000002, 0.11203125000000003, 0.8860823156553729), 0, 100, null, V(0.0023109367883072116, -3.1668370133166326e-21, 0.00000784814731787084), V(0, -0.0022163159641569587, 0.0002802184892673451)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000001, 0, 0), V(0, 0, -0.8931365355273235), V(0, 0.8931365355273235, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.11203125000000003, 0.8860823156553729), V(0.11093750000000002, 0.07093749999999999, 0.8903149679677975), 3.015825616877026, 3.0620837623516257, null, V(0, -0.8860823156553727, 0.11203125000000015), V(0, -0.8903149679677973, 0.0709375)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0.07093749999999997, 0), V(0, 0, -0.8972000173282154), V(0.8972000173282154, 0, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.07093749999999999, 0.8903149679677975), V(0.04296875, 0.07093749999999999, 0.8961704958417165), 3.0176268184353994, 3.093682274615096, null, V(-0.8903149679677973, 0, 0.11093750000000019), V(-0.8961704958417163, 0, 0.04296875000000003)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(0.04296874999999999, 0.07093749999999997, 0.8961704958417162), V(-0.17104693965321316, 3.122502256758253e-18, 0.8835965959844289), -1), V(0.04296875, 0.07093749999999999, 0.8961704958417165), V(-0.17104693965321202, 2.0947208715025658e-17, 0.8835965959844292), 0, 64, null, V(-0.004798729293096285, -5.284086361068957e-20, 0.000230085012025602), V(-0.0023970206921961742, -0.0028148248536625227, -0.00046401610819787214)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(-0.9, 1.1021821192326179e-16, 0), 0, 3.141592653589793), V(-0.17104693965321202, 2.0947208715025658e-17, 0.8835965959844292), V(0, 0, -0.9), 2.9503773839343412, 0, null, V(-0.8835965959844292, 1.0820937430098282e-16, -0.17104693965321202), V(0.9, -1.1021821192326179e-16, 0)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(0.9, 0, 0)), V(0, 0, -0.9), V(0.20296875000000006, 0, 0.8768145108992196), 0, 2.9141150448015316, null, V(0.9, 0, 0), V(-0.8768145108992196, 0, 0.20296875000000006)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.20296874999999995, 0, 0), V(0, 0, 0.8768145108992196), V(0, 0.8768145108992196, 0), 0, 3.141592653589793), V(0.20296875000000006, 0, 0.8768145108992196), V(0.20296874999999998, 0.11703124999999998, 0.8689691438980301), 0, 0.13387273064879698, null, V(0, 0.8768145108992196, 0), V(0, 0.86896914389803, -0.11703124999999995))], []),
	//		new RotationFace(new SemiEllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(-0.03286362438479655, -2.3418766925686897e-18, 0.8993997899667839), V(0.04796875, 0.010000000000000007, 0.8986651206224917), -1), V(-0.03286362438479702, 4.024633241122258e-18, 0.8993997899667838), V(0.047968750000000004, 0.010000000000000009, 0.8986651206224918), 0, 29, null, V(0.002292032896076846, 0.0007036417545401136, 0.00008374975068322923), V(0.003388395809325897, 0, -0.00018086504944802423)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0.010000000000000007, 0), V(0, 0, 0.899944442729661), V(0.899944442729661, 0, 0), 0, 3.141592653589793), V(0.047968750000000004, 0.010000000000000009, 0.8986651206224918), V(0.11093750000000002, 0.010000000000000009, 0.8930805512907278), 0.053327173050216115, 0.12358585566041683, null, V(0.8986651206224918, 0, -0.04796875), V(0.8930805512907277, 0, -0.11093750000000001)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000001, 0, 0), V(0, 0, -0.8931365355273235), V(0, 0.8931365355273235, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.010000000000000009, 0.8930805512907278), V(0.11093750000000006, 0, 0.8931365355273235), 3.130395923246913, 3.141592653589793, null, V(0, -0.8930805512907277, 0.01000000000000018), V(0, -0.8931365355273235, 1.093776799435093e-16)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(0.9, 0, 0)), V(0.11093750000000006, 0, 0.8931365355273235), V(0, 0, 0.9), 3.018014465996861, 3.141592653589793, null, V(-0.8931365355273235, 0, 0.11093750000000006), V(-0.9, 0, 0)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(-0.9, 1.1021821192326179e-16, 0), 0, 3.141592653589793), V(0, 0, 0.9), V(-0.03286362438479702, 4.024633241122258e-18, 0.8993997899667838), 3.141592653589793, 3.1050693959027966, null, V(-0.9, 1.1021821192326179e-16, 0), V(-0.8993997899667838, 1.1014470739366237e-16, -0.03286362438479702))], []),
	//		new RotationFace(new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.17104693965321313, 2.0947208715025794e-17, 0.883596595984429), V(-0.205, -0.09906249999999998, 0.8707247677043246), -1), V(-0.17104693965321202, 2.0947208715025658e-17, 0.8835965959844292), V(-0.20500000000000002, -0.0990625, 0.8707247677043248), 0, 36, null, V(-0.00215411320093391, -0.002529578236571627, -0.00041699399065229685), V(0, -0.0032695078711222625, -0.0003719724481215541)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.20500000000000002, -0.0990625, 1.2), V(-0.20500000000000002, -0.19703125, 1.2), V(-0.14, -0.26203125, 1.2), V(-0.034062499999999996, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.205, -0.09906249999999998, 0.8707247677043246), V(-0.034062499999999996, -0.26203124999999994, 0.8603367771501969), -1), V(-0.20500000000000002, -0.0990625, 0.8707247677043248), V(-0.034062499999999996, -0.26203125, 0.8603367771501972), 0, 100, null, V(0, -0.0029374208172160575, -0.0003341908493916618), V(0.003177872773799387, 0, 0.0001258185099515327)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.034062499999999996, -0.26203125, 1.2), V(0.030000000000000027, -0.26203125, 1.2), V(0.0859375, -0.23796875, 1.2), V(0.12296875000000002, -0.18796875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.034062499999999996, -0.26203124999999994, 0.8603367771501969), V(0.12296874999999999, -0.18796874999999996, 0.8715196128297258), -1), V(-0.034062499999999996, -0.26203125, 0.8603367771501972), V(0.12296875000000002, -0.18796875, 0.871519612829726), 0, 100, null, V(0.0019218193657008743, 0, 0.00007608877579431633), V(0.0011107827301533816, 0.0014997910280551989, 0.0001667458526657328)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.12296875000000002, -0.18796875, 1.2), V(0.13906249999999998, -0.235, 1.2), V(0.17296875, -0.2559375, 1.2), V(0.21999999999999997, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(0.12296874999999999, -0.18796874999999996, 0.8715196128297258), V(0.21999999999999995, -0.2620312499999999, 0.832429951421402), -1), V(0.12296875000000002, -0.18796875, 0.871519612829726), V(0.21999999999999997, -0.26203125, 0.8324299514214022), 0, 100, null, V(0.0004824780014564421, -0.0014099599848387352, -0.0003721753680202314), V(0.0014095926581339242, -0.00018263825138612356, -0.00043002695120080944)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.2773154018890211, -0.09536944308866366, 0), V(-0.27671434201165657, -0.8046303562041052, 0), V(0, 0, 0.8508823874073836), 0, 1.6896013946143438), V(0.21999999999999997, -0.26203125, 0.8324299514214022), V(0.24203124999999998, -0.19796875, 0.8439367559520533), 1.3621575275257576, 1.4429371319818958, null, V(0.27071344957582744, 0.7871809526672972, 0.17624191662783154), V(0.2744555623419146, 0.7980622734764866, 0.10849695458036643)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.24203124999999998, -0.19796875, 1.2), V(0.21500000000000002, -0.18796875, 1.2), V(0.20296874999999998, -0.17296875, 1.2), V(0.20296874999999998, -0.13093749999999998, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(0.2420312499999999, -0.19796874999999997, 0.843936755952053), V(0.20296874999999995, -0.13093749999999996, 0.8669827320178802), -1), V(0.24203124999999998, -0.19796875, 0.8439367559520533), V(0.20296874999999998, -0.13093749999999998, 0.8669827320178805), 0, 100, null, V(-0.0008105656446711983, 0.000299862435022873, 0.00030280184601271074), V(1.6025994766705188e-19, 0.0012607132877938707, 0.00019040130792028905)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.20296874999999995, 4.943121956121308e-19, 0), V(0, 0, -0.8768145108992196), V(2.1354031397797386e-18, -0.8768145108992196, 0), 2.891240380027181e-17, 3.141592653589793), V(0.20296874999999998, -0.13093749999999998, 0.8669827320178805), V(0.20296875000000006, 0, 0.8768145108992196), 2.9916987961562977, 3.141592653589793, null, V(-2.111458723678207e-18, 0.8669827320178803, 0.13093749999999996), V(-2.1354031397797386e-18, 0.8768145108992196, 1.0737880842186813e-16)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(0.9, 0, 0)), V(0.20296875000000006, 0, 0.8768145108992196), V(0, 0, -0.9), 2.9141150448015316, 0, null, V(0.8768145108992196, 0, -0.20296875000000006), V(-0.9, 0, 0)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(-0.9, 1.1021821192326179e-16, 0), 0, 3.141592653589793), V(0, 0, -0.9), V(-0.17104693965321202, 2.0947208715025658e-17, 0.8835965959844292), 0, 2.9503773839343412, null, V(-0.9, 1.1021821192326179e-16, 0), V(0.8835965959844292, -1.0820937430098282e-16, 0.17104693965321202))], []),
	//		new RotationFace(new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.009062500000000001, -0.19296875, 1.2), V(-0.07500000000000001, -0.19296875, 1.2), V(-0.1040625, -0.16, 1.2), V(-0.1040625, -0.095, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.0090625, -0.19296874999999997, 0.879022714505824), V(-0.10406249999999999, -0.09499999999999999, 0.8889015671567634), -1), V(-0.009062500000000001, -0.19296875, 0.8790227145058244), V(-0.1040625, -0.095, 0.8889015671567637), 0, 100, null, V(-0.001978120886365625, -1.734719868511089e-20, -0.000020393921837124158), V(0, 0.0019495768662389622, 0.0002083580557576388)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.10406249999999999, -0.09499999999999999, 0.8889015671567634), V(-0.03286362438479653, 4.0246332411221976e-18, 0.8993997899667839), -1), V(-0.1040625, -0.095, 0.8889015671567637), V(-0.03286362438479702, 4.024633241122258e-18, 0.8993997899667838), 0, 71, null, V(0, 0.0019214697267915689, 0.00020535414807408821), V(0.0022509978133461728, 0.0006910442051507778, 0.00008225034901502365)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(-0.9, 1.1021821192326179e-16, 0), 0, 3.141592653589793), V(-0.03286362438479702, 4.024633241122258e-18, 0.8993997899667838), V(0, 0, 0.9), 3.1050693959027966, 3.141592653589793, null, V(0.8993997899667838, -1.1014470739366237e-16, 0.03286362438479702), V(0.9, -1.1021821192326179e-16, 0)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(0.9, 0, 0)), V(0, 0, 0.9), V(0.11093750000000006, 0, 0.8931365355273235), 3.141592653589793, 3.018014465996861, null, V(0.9, 0, 0), V(0.8931365355273235, 0, -0.11093750000000006)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000001, 2.701783363237952e-19, 0), V(0, 0, 0.8931365355273235), V(2.175153967583285e-18, -0.8931365355273235, 0), 1.5513981584219775e-17, 3.141592653589793), V(0.11093750000000006, 0, 0.8931365355273235), V(0.11093750000000002, -0.12, 0.8850383444200315), 1.5513981584219775e-17, 0.13476551593405545, null, V(2.175153967583285e-18, -0.8931365355273235, -2.7017833632379536e-19), V(2.155431549099e-18, -0.8850383444200314, -0.11999999999999998)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.11093750000000002, -0.12, 1.2), V(0.0809375, -0.16593750000000002, 1.2), V(0.040937500000000016, -0.19296875, 1.2), V(-0.009062500000000001, -0.19296875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(0.11093750000000001, -0.11999999999999997, 0.8850383444200313), V(-0.0090625, -0.19296874999999997, 0.879022714505824), -1), V(0.11093750000000002, -0.12, 0.8850383444200315), V(-0.009062500000000001, -0.19296875, 0.8790227145058244), 0, 100, null, V(-0.0008999753301001092, -0.0013780872242157918, -0.00007404137248522919), V(-0.0014999982255756404, 0, -0.000015464599145110234))], []),
	//		new RotationFace(new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.20296874999999998, 0.11703124999999998, 1.2), V(0.20296874999999998, 0.2240625, 1.2), V(0.14500000000000002, 0.2890625, 1.2), V(0.010937499999999989, 0.2890625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), V(0.010937499999999989, 0.2890625, 0.9572477433702835), 1), V(0.010937499999999989, 0.2890625, 0.9572477433702835), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), 100, 0, null, V(0.004021833637101672, 0, -0.00004595341771287279), V(0, -0.003208541417521995, 0.0003862503783657504)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.20296874999999998, 0, 0), V(0, 0, -0.9791852156376941), V(0, 0.9791852156376941, 0), 0, 3.141592653589793), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), V(0.2029687500000001, 0, 0.9791852156376941), 3.0217872455230754, 3.141592653589793, null, V(0, -0.9721663299286162, 0.1170312499999999), V(0, -0.9791852156376941, 1.199156040103113e-16)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0.2029687500000001, 0, 0.9791852156376941), V(0, 0, -1), 2.937203822794261, 0, null, V(0.9791852156376941, 0, -0.2029687500000001), V(-1, 0, 0)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(0, 0, -1), V(-0.171046939653212, 2.0947208715025655e-17, 0.9852628808776216), 0, 2.969700482947388, null, V(-1, 1.2246467991473532e-16, 0), V(0.9852628808776216, -1.2065990333854792e-16, 0.171046939653212)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.04296875, 0.07093749999999999, 0.9965548442595558), V(-0.1710469396532132, 6.938893903907228e-18, 0.9852628808776215), 1), V(-0.171046939653212, 2.0947208715025655e-17, 0.9852628808776216), V(0.04296875, 0.07093749999999999, 0.9965548442595558), 64, 0, null, V(0.002384621933525827, 0.0028002649734856517, 0.00041398320374781613), V(0.0047989723215347765, 0, -0.0002069187091194757)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0.07093749999999999, 0), V(0, 0, 0.9974807622674986), V(0.9974807622674986, 0, 0), 0, 3.141592653589793), V(0.04296875, 0.07093749999999999, 0.9965548442595558), V(0.11093750000000002, 0.07093749999999999, 0.9912924604714292), 0.04309060575909555, 0.11144825166905915, null, V(0.9965548442595558, 0, -0.04296875), V(0.9912924604714292, 0, -0.11093750000000004)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000002, 0, 0), V(0, 0, 0.9938273849586506), V(0, 0.9938273849586506, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.07093749999999999, 0.9912924604714292), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), 0.07143883874192354, 0.11296719097084001, null, V(0, 0.9912924604714292, -0.07093749999999999), V(0, 0.9874927190198353, -0.11203125000000001)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.0029687499999999922, 0.2140625, 1.2), V(0.07406249999999998, 0.2140625, 1.2), V(0.11093750000000002, 0.18593749999999998, 1.2), V(0.11093750000000002, 0.11203125000000003, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), 1), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), 100, 0, null, V(0, 0.002216485816645523, -0.00025146076711587485), V(-0.0023109369300238886, 0, -0.000007023429018985676)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1409375, 0.18703124999999998, 1.2), V(-0.08703125, 0.2059375, 1.2), V(-0.037968749999999996, 0.2140625, 1.2), V(-0.0029687499999999922, 0.2140625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(-0.1409375, 0.18703124999999998, 0.9721913045369145), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), 1), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), V(-0.1409375, 0.18703124999999998, 0.9721913045369144), 100, 0, null, V(-0.0010499999469858546, 0, -0.0000031911732431043194), V(-0.001617060512887581, -0.0005671429624910072, -0.00012531587997459948)),
	//			new PCurveEdge(new SemiEllipseCurve(V(-0.06640184370318536, -0.0235077791500932, 0), V(0, 0, 0.997516004619599), V(-0.33289781807779234, 0.9403282523625955, 0), 0.02500215049153284, 3.1165905030982604), V(-0.1409375, 0.18703124999999998, 0.9721913045369144), V(-0.16499999999999998, 0.255, 0.9527591510974849), 0.22581372635962282, 0.3006922248899664, null, V(-0.32444628711291806, 0.9164554213903857, -0.22334333850612303), V(-0.31796125684715715, 0.8981373164189183, -0.29544573016418446)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.010937499999999989, 0.2890625, 1.2), V(-0.04093749999999999, 0.2890625, 1.2), V(-0.1, 0.2790625, 1.2), V(-0.16499999999999998, 0.255, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.010937499999999989, 0.2890625, 0.9572477433702835), V(-0.16499999999999998, 0.255, 0.9527591510974849), 1), V(-0.16499999999999998, 0.255, 0.9527591510974849), V(0.010937499999999989, 0.2890625, 0.9572477433702835), 100, 0, null, V(0.0019497960757584164, 0.0007217995088144136, 0.0001444829762001187), V(0.0015562475396635686, 0, -0.000017781663715540364))], []),
	//		new RotationFace(new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(-0.03286362438479655, -2.6020852139652106e-18, 0.9994598452125503), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), 1), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), 29, 0, null, V(-0.003388522845144273, 0, 0.0001627386905160511), V(-0.0022900449065656465, -0.0007030314524671974, -0.00007529984920773373)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), V(0, 0, 1), 3.108723110778215, 3.141592653589793, null, V(0.9994598452125503, -1.2239853003158589e-16, 0.032863624384797126), V(1, -1.2246467991473532e-16, 0)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0, 0, 1), V(0.11093749999999993, 0, 0.9938273849586506), 3.141592653589793, 3.030426330354509, null, V(1, 0, 0), V(0.9938273849586506, 0, -0.11093749999999993)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000002, 0, 0), V(0, 0, 0.9938273849586506), V(0, 0.9938273849586506, 0), 0, 3.141592653589793), V(0.11093749999999993, 0, 0.9938273849586506), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), 0, 0.010062279327831384, null, V(0, 0.9938273849586506, 0), V(0, 0.993777073137507, -0.010000000000000007)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0.010000000000000009, 0), V(0, 0, -0.9999499987499375), V(0.9999499987499375, 0, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), 3.0304207486077566, 3.0936030871101416, null, V(-0.9937770731375071, 0, 0.11093749999999983), V(-0.9987987780446257, 0, 0.0479687500000002))], []),
	//		new RotationFace(new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.17104693965321316, 2.0947208715025797e-17, 0.9852628808776215), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), 1), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), V(-0.171046939653212, 2.0947208715025655e-17, 0.9852628808776216), 36, 0, null, V(0, 0.003269989649016917, 0.0003326706415016505), V(0.0018504236016956069, 0.002172955102479456, 0.000321243497826874)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(-0.171046939653212, 2.0947208715025655e-17, 0.9852628808776216), V(0, 0, -1), 2.969700482947388, 0, null, V(-0.9852628808776216, 1.2065990333854792e-16, -0.171046939653212), V(1, -1.2246467991473532e-16, 0)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0, 0, -1), V(0.2029687500000001, 0, 0.9791852156376941), 0, 2.937203822794261, null, V(1, 0, 0), V(-0.9791852156376941, 0, 0.2029687500000001)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.20296874999999998, -2.3224450485052117e-18, 0), V(0, 0, 0.9791852156376941), V(-1.1204206832959602e-17, -0.9791852156376941, 0), 2.3013070043406876e-17, 3.141592653589793), V(0.2029687500000001, 0, 0.9791852156376941), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), 2.3013070043406876e-17, 0.13412262887088133, null, V(-1.1204206832959602e-17, -0.9791852156376941, -2.3224450485052132e-18), V(-1.1103582248632314e-17, -0.9703911879325716, -0.13093749999999996)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.24203124999999998, -0.19796875, 1.2), V(0.21500000000000002, -0.18796875, 1.2), V(0.20296874999999998, -0.17296875, 1.2), V(0.20296874999999998, -0.13093749999999998, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.24203124999999998, -0.19796875, 0.9498574882827818), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), 1), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), V(0.24203124999999998, -0.19796875, 0.9498574882827818), 100, 0, null, V(0, -0.0012607577897218358, -0.00017011744865842043), V(0.0008106439109762247, -0.00029989138903166777, -0.0002690617125763413)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.27731540188902115, -0.09536944308866366, 0), V(0.31091053038645644, 0.9040660812655796, 0), V(0, 0, 0.9560339100680944), 1.465110231866663, 3.141592653589793), V(0.24203124999999998, -0.19796875, 0.9498574882827818), V(0.21999999999999997, -0.26203125, 0.9396486705271484), 1.684527864719763, 1.7562036900306723, null, V(-0.30890190438173115, -0.8982253957199245, -0.10849695458036647), V(-0.30558190818745745, -0.8885715060770011, -0.17624191662783179)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.12296875000000002, -0.18796875, 1.2), V(0.13906249999999998, -0.235, 1.2), V(0.17296875, -0.2559375, 1.2), V(0.21999999999999997, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.12296875000000002, -0.18796875, 0.9744467330474638), V(0.21999999999999997, -0.26203125, 0.9396486705271484), 1), V(0.21999999999999997, -0.26203125, 0.9396486705271484), V(0.12296875000000002, -0.18796875, 0.9744467330474638), 100, 0, null, V(-0.0014098883430614425, 0.00018267656272224708, 0.0003810385889743954), V(-0.0004825448774820869, 0.0014101554186612473, 0.00033291003063843776)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.034062499999999996, -0.26203125, 1.2), V(0.030000000000000027, -0.26203125, 1.2), V(0.0859375, -0.23796875, 1.2), V(0.12296875000000002, -0.18796875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.034062499999999996, -0.26203125, 0.9644580706890205), V(0.12296875000000002, -0.18796875, 0.9744467330474638), 1), V(0.12296875000000002, -0.18796875, 0.9744467330474638), V(-0.034062499999999996, -0.26203125, 0.9644580706890203), 100, 0, null, V(-0.0011108137429678576, -0.0014998329018975289, -0.00014913728313140515), V(-0.0019218307292635764, 0, -0.00006787475910567421)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.20500000000000002, -0.0990625, 1.2), V(-0.20500000000000002, -0.19703125, 1.2), V(-0.14, -0.26203125, 1.2), V(-0.034062499999999996, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), V(-0.034062499999999996, -0.26203125, 0.9644580706890205), 1), V(-0.034062499999999996, -0.26203125, 0.9644580706890203), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), 100, 0, null, V(-0.003177924394426102, 0, -0.00011223717544071708), V(0, 0.002937749569980493, 0.00029887037541859586))], []),
	//		new RotationFace(new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), [
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.009062500000000001, -0.19296875, 1.2), V(-0.07500000000000001, -0.19296875, 1.2), V(-0.1040625, -0.16, 1.2), V(-0.1040625, -0.095, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), V(-0.1040625, -0.095, 0.9900232300778351), 1), V(-0.1040625, -0.095, 0.9900232300778351), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), 100, 0, null, V(0, -0.0019496588340396214, -0.00018708408409689773), V(0.001978121698253082, 0, 0.000018270895823513832)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.11093750000000002, -0.12, 1.2), V(0.0809375, -0.16593750000000002, 1.2), V(0.040937500000000016, -0.19296875, 1.2), V(-0.009062500000000001, -0.19296875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.11093750000000002, -0.12, 0.9865560658643532), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), 1), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), V(0.11093750000000002, -0.12, 0.9865560658643532), 100, 0, null, V(0.001499998572758444, 0, 0.000013854717676112659), V(0.0008999801458405889, 0.001378094598318401, 0.00006642278795539617)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000002, 2.7017833632379675e-19, 0), V(0, 0, -0.9938273849586506), V(2.4203775050019737e-18, -0.9938273849586506, 0), 1.3942163371701867e-17, 3.141592653589793), V(0.11093750000000002, -0.12, 0.9865560658643532), V(0.11093749999999993, 0, 0.9938273849586506), 3.02055199779192, 3.141592653589793, null, V(-2.4026688591808874e-18, 0.9865560658643532, 0.11999999999999994), V(-2.4203775050019737e-18, 0.9938273849586506, 1.217087525894596e-16)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0.11093749999999993, 0, 0.9938273849586506), V(0, 0, 1), 3.030426330354509, 3.141592653589793, null, V(-0.9938273849586506, 0, 0.11093749999999993), V(-1, 0, 0)),
	//			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(0, 0, 1), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), 3.141592653589793, 3.108723110778215, null, V(-1, 1.2246467991473532e-16, 0), V(-0.9994598452125503, 1.2239853003158589e-16, -0.032863624384797126)),
	//			new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.1040625, -0.095, 0.9900232300778351), V(-0.03286362438479653, 2.7235906341395923e-18, 0.9994598452125503), 1), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), V(-0.1040625, -0.095, 0.9900232300778351), 71, 0, null, V(-0.0022530279161436587, -0.0006916674357757679, -0.00007408267927848543), V(0, -0.0019215482684322, -0.00018438666887312053))], [])], false)
	//	b2Equal(assert, d, c, d.and(c), result)
	//},

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
        b2EqualAnd(assert, a, b, result)
    },

     //'B2.withMergedFaces'(assert) {
    //    const box = B2T.box(5, 5, 5)
    //    const boxToMerge = new B2(box.faces.filter((face:PlaneFace) => face.surface.plane.normal1.x != 1).concat(
    //        box.translate(5, 0, 0).faces.filter((face:PlaneFace) => face.surface.plane.normal1.x != -1)
    //    ), false)
    //
    //    assert.equal(boxToMerge.faces.length, 10)
    //    const boxMerged = boxToMerge.withMergedFaces()
    //    b2Equal(assert, boxToMerge, B2.EMPTY, boxMerged, B2T.box(10, 5, 5))
    //}
})