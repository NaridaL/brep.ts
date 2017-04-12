{
	QUnit.module('SemiCylidnderSurface')
	QUnit.test('testSurface', function (assert) {
		const ps = new SemiCylinderSurface(SemiEllipseCurve.UNIT, V3.Z, 0, 1)
		testParametricSurface(assert, ps)
		testParametricSurface(assert, ps.rotateZ(PI))
	})
	registerTests({
		'SemiCylinderSurface.isCurvesWithSurface'(assert) {
			const cyl = SemiCylinderSurface.semicylinder(5)
			const ell = new SemiCylinderSurface(new SemiEllipseCurve(V(6, 1, 4), V(3, 1, 4), V(4, 0, 0)), V3.Z)
			testISCurves(assert, cyl, ell, 1)
		},
		'SemiCylinderSurface.intersectionLine 2'(assert) {
			const cylSurface = new SemiCylinderSurface(new SemiEllipseCurve(V3.O, V(8, 0, 0), V(0, 5, 0)), V3.Z)
			const line = L3.throughPoints(V(10, 0, 0), V(-10, 2, 10))
			testISTs(assert, line, cylSurface, 2)
		},
	})
}