{
	QUnit.module('SemiCylinderSurface')
	registerTests({
		'testSurface'(assert) {
			const ps = new SemiCylinderSurface(SemiEllipseCurve.UNIT, V3.Z, 0, 1)
			testParametricSurface(assert, ps)
			testParametricSurface(assert, ps.rotateZ(PI))
		},
		'is curves w/ SemiCylinderSurface'(assert) {
			const cyl = SemiCylinderSurface.UNIT.scale(5,5,1)
			const ell = new SemiCylinderSurface(new SemiEllipseCurve(V(6, 1, 4), V(3, 1, 4), V(4, 0, 0)), V3.Z)
			testISCurves(assert, cyl, ell, 1)
		},
		'is curves w/ SemiCylinderSurface 2'(assert) {
			const cyl = SemiCylinderSurface.UNIT.scale(5,5,1)
			const scs = new SemiCylinderSurface(new SemiEllipseCurve(V(-1.5, 2.5980762113533173, 0), V(1.250, 2.1650635094610964, 0), V(-2.165063509461096, 1.25, 0), 0, 3.141592653589793), V(0, 0, -1), -2, 2)
			testISCurves(assert, cyl, scs, 2)
		},
		'is curves w/ SemiEllipsoidSurface'(assert) {
			const cyl = SemiCylinderSurface.UNIT.scale(0.05,0.05, 4).translate(0.5,0.5,-2)
			const sphere = SemiEllipsoidSurface.UNIT
			testISCurves(assert, sphere, cyl, 2)
		},
		'is curves w/ SemiEllipsoidSurface 2'(assert) {
			const cyl = SemiCylinderSurface.UNIT.scale(0.5,0.05, 4).translate(0,0,-2)
			const sphere = SemiEllipsoidSurface.UNIT
			testISCurves(assert, sphere, cyl, 2)
		},
		'is curves w/ SemiEllipsoidSurface 3'(assert) {
			const cyl = SemiCylinderSurface.UNIT.rotateZ(PI).scale(0.5,0.05, 4).translate(0.25,0,-2)
			const sphere = SemiEllipsoidSurface.UNIT
			testISCurves(assert, sphere, cyl, 0)
		},
		'is curves w/ plane'(assert) {
			testISCurves(assert, SemiCylinderSurface.UNIT, P3.XY, 1)
			testISCurves(assert, SemiCylinderSurface.UNIT, P3.XY.flipped(), 1)
			testISCurves(assert, SemiCylinderSurface.UNIT.flipped(), P3.XY, 1)
			testISCurves(assert, SemiCylinderSurface.UNIT.flipped(), P3.XY.flipped(), 1)


			testISCurves(assert,
				new SemiCylinderSurface(new SemiEllipseCurve(V(0.5, 0.2, 0), V(-0.2, 2.4492935982947065e-17, 0), V(-2.4492935982947065e-17, -0.2, 0), 0, 3.141592653589793), V(0, 0, -1), -Infinity, Infinity),
				new PlaneSurface(new P3(V(0, -1, 0), 0), V(0, 0, -1), V(1, 0, 0)),
				1)


		},
		'intersectionLine 2'(assert) {
			const cylSurface = new SemiCylinderSurface(new SemiEllipseCurve(V3.O, V(8, 0, 0), V(0, 5, 0)), V3.Z)
			const line = L3.throughPoints(V(10, 0, 0), V(-10, 2, 10))
			testISTs(assert, line, cylSurface, 2)
		},
		'intersectionLine'(assert) {
			const cylSurface = SemiCylinderSurface.semicylinder(5)
			const line = L3.throughPoints(V(10, 0, 0), V(-10, 2, 10))
			const isPoints = cylSurface.isTsForLine(line).map(line.at, line)

			assert.equal(isPoints.length, 2, 'no of points')
			assert.notOk(isPoints[0].like(isPoints[1]))

			assert.ok(cylSurface.containsPoint(isPoints[0]))
			assert.ok(cylSurface.containsPoint(isPoints[1]))

			assert.ok(line.containsPoint(isPoints[0]), '' + line.distanceToPoint(isPoints[0]))
			assert.ok(line.containsPoint(isPoints[1]), '' + line.distanceToPoint(isPoints[1]))
		},
		'zDirVolume'(assert) {
			const face = B2T.extrudeEdges([Edge.forCurveAndTs(SemiEllipseCurve.forAB(-1, 1), 0, PI), StraightEdge.throughPoints(V3.X, V3.X.negated())], P3.XY.flipped(), V3.Z, 'cyl').faces.find(face => face.surface instanceof SemiCylinderSurface)
			const face2 = B2T.extrudeEdges([Edge.forCurveAndTs(SemiEllipseCurve.UNIT, PI, 0), StraightEdge.throughPoints(V3.X, V3.X.negated())], P3.XY.flipped(), V3.Z, 'cyl').faces.find(face => face.surface instanceof SemiCylinderSurface)
			const face3 = B2T.extrudeEdges([Edge.forCurveAndTs(SemiEllipseCurve.UNIT, PI, 0).rotateY(-80 * DEG), StraightEdge.throughPoints(V3.X, V3.X.negated()).rotateY(-80 * DEG)], P3.XY.flipped().rotateY(-80 * DEG), new V3(-10, -1, 0).unit(), 'cyl').faces.find(face => face.surface instanceof SemiCylinderSurface)
			const modface = face.rotateY(-45 * DEG).translate(1, 0, 2)
			const e0 = modface.contour[0].project(new P3(modface.surface.dir1, 0))
			const face4 = Face.create(modface.surface, [e0, StraightEdge.throughPoints(e0.b, modface.contour[2].a), modface.contour[2], StraightEdge.throughPoints(modface.contour[2].b, e0.a)])

			testZDirVolume(assert, face)
			testZDirVolume(assert, face.rotateY(-45 * DEG).translate(1, 0, 2))
			testZDirVolume(assert, face.rotateY(90 * DEG).translate(1, 0, 2))

			testZDirVolume(assert, face2)
			testZDirVolume(assert, face2.rotateY(-45 * DEG).translate(1, 0, 2))
			testZDirVolume(assert, face2.rotateY(90 * DEG).translate(1, 0, 2))

			testZDirVolume(assert, face3)
			testZDirVolume(assert, face3.translate(1, 0, 2))

			testZDirVolume(assert, face4)
			testZDirVolume(assert, face4.translate(1, 0, 2))
		},
		'loopContainsPoint'(assert) {
			const surface = new SemiCylinderSurface(new SemiEllipseCurve(V(0, 0, 0), V(8, 0, 0), V(0, 8, 0)), V(0, 0, -1))
			const loop = [
				StraightEdge.throughPoints(V(1, 7.937253933193773, 4), V(1, 7.937253933193773, 1)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 1), V(8, 0, 0), V(0, 8, 0)), V(1, 7.937253933193773, 1), V(6, 5.291502622129181, 1), 1.4454684956268313, 0.7227342478134156, null, V(7.937253933193772, -0.9999999999999991, 0), V(5.2915026221291805, -6, 0)),
				StraightEdge.throughPoints(V(6, 5.291502622129181, 1), V(6, 5.291502622129181, 4)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 4), V(8, 0, 0), V(0, 8, 0)), V(6, 5.291502622129181, 4), V(1, 7.937253933193773, 4), 0.7227342478134156, 1.4454684956268313, null, V(-5.2915026221291805, 6, 0), V(-7.937253933193772, 0.9999999999999991, 0))
			]
			testLoopContainsPoint(assert, surface, loop, V(8, 0, 0), PointVsFace.OUTSIDE)
			testLoopContainsPoint(assert, surface, loop, V(1, 7.937253933193773, 3), PointVsFace.ON_EDGE)
		},
	})
}