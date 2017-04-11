{
	QUnit.module('SemiEllipsoidSurface')
	const ses2 = SemiEllipsoidSurface.UNIT.scale(2)
	QUnit.test('testSurface', function (assert) {
		testParametricSurface(assert, ses2)
	})
	QUnit.test('testSurface', function (assert) {
		testISCurves(assert, SemiEllipsoidSurface.UNIT, new PlaneSurface(new P3(V(-1.249000902703301e-16, 1, 0), 0.11006944444444443)), 2)
	})
	QUnit.test('getAABB', function (assert) {
		assert.V3ArraysLike(ses2.getExtremePoints(), [V(0,2,0)])
		assert.V3ArraysLike(ses2.rotateZ(30 * DEG).getExtremePoints(), [V(-2,0,0),V(0,2,0)])
	})
}