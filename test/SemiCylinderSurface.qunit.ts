{
	QUnit.module('SemiCylidnderSurface')
	QUnit.test('testSurface', function (assert) {
		const ps = new SemiCylinderSurface(SemiEllipseCurve.UNIT, V3.Z, 0, 1)
		testParametricSurface(assert, ps)
		testParametricSurface(assert, ps.rotateZ(PI))
	})
	registerTests({
		'is curves w/ SemiCylinderSurface'(assert) {
			const cyl = SemiCylinderSurface.semicylinder(5)
			const ell = new SemiCylinderSurface(new SemiEllipseCurve(V(6, 1, 4), V(3, 1, 4), V(4, 0, 0)), V3.Z)
			testISCurves(assert, cyl, ell, 1)
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
		'SemiCylinderSurface.intersectionLine 2'(assert) {
			const cylSurface = new SemiCylinderSurface(new SemiEllipseCurve(V3.O, V(8, 0, 0), V(0, 5, 0)), V3.Z)
			const line = L3.throughPoints(V(10, 0, 0), V(-10, 2, 10))
			testISTs(assert, line, cylSurface, 2)
		},
	})
}