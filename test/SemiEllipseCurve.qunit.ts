QUnit.module('SemiEllipseCurve')
{
	const curve = SemiEllipseCurve.UNIT.shearX(2, 1)
	registerTests({
		'withBounds'(assert) {
			const newCurve = curve.withBounds(1, 2)
			assert.equal(newCurve.tMin, 1)
			assert.equal(newCurve.tMax, 2)
		},
		'testCurve'(assert) {
			testCurve(assert, SemiEllipseCurve.UNIT)
			testCurve(assert, curve)
			testCurve(assert, curve.reversed())
		},
		'isTsWithPlane'(assert) {
			const plane = new P3(V(2, 7, 1).unit(), 2)
			testISTs(assert, curve.scale(1, 3, 1), plane, 2)
		},
		//'rightAngled'(assert) {
		//	const curveRA = curve.rightAngled()
		//	assert.ok(curveRA.f1.isPerpendicularTo(curveRA.f2))
		//	assert.ok(curveRA.isColinearTo(curve))
		//	arrayRange(-10, 10, 1).forEach(t => assert.ok(curveRA.containsPoint(curve.at(t))))
		//},
		'isTsWithSurface(SemiEllipsoidSurface)'(assert) {
			const s = SemiEllipsoidSurface.sphere(5)
			const c = new SemiEllipseCurve(V(5, 2), V3.Z.negated(), V(-1, 1.2246467991473532e-16, 0), 0, PI)
			testISTs(assert, c, s, 2)
		},
		'isTsWithSurface(PlaneSurface)'(assert) {
			const c = SemiEllipseCurve.UNIT.translate(1.2, -1)
			const s = new PlaneSurface(P3.ZX)
			testISTs(assert, c, s, 1)
		},
		'distanceToPoint'(assert) {
			const curve = SemiEllipseCurve.forAB(10, 15)
			const p = V(12, 12)
			const closestT = curve.closestTToPoint(p)
			const pDist = curve.at(closestT).distanceTo(p)
			const EPS = 0.001
			assert.push(pDist < curve.at(closestT - EPS).distanceTo(p), curve.at(closestT - EPS).distanceTo(p), '> ' + pDist,
				'' + (pDist - curve.at(closestT - EPS).distanceTo(p)) + 'larger')
			assert.push(pDist < curve.at(closestT + EPS).distanceTo(p), curve.at(closestT + EPS).distanceTo(p), '> ' + pDist)
		},
		'isColinearTo'(assert) {
			assert.ok(SemiEllipseCurve.forAB(1, 2).isColinearTo(SemiEllipseCurve.forAB(1, -2)))
		},
		'isInfosWithEllipse'(assert) {
			const c1 = SemiEllipseCurve.semicircle(5), c2 = SemiEllipseCurve.semicircle(5, V(3, 0))
			testCurveISInfos(assert, c1, c2, 1)

			const verticalEllipse = new SemiEllipseCurve(V(2, 0), V(1, 1), V(1, 10))
			testCurveISInfos(assert, c1, verticalEllipse, 2)

			const verticalEllipse2 = new SemiEllipseCurve(V(10, 2), V(1, 1), V(1, 10))
			testCurveISInfos(assert, c1, verticalEllipse2, 0)

			const smallEllipse = SemiEllipseCurve.forAB(2, 3)
			testCurveISInfos(assert, c1, smallEllipse, 0)

			const test = new SemiEllipseCurve(V(6, 1, 0), V(3, 1, 0), V(4, 0, 0))
			testCurveISInfos(assert, c1, test, 1)
		},
		'EllipseCurve.isInfosWithEllipse'(assert) {
			const c1 = EllipseCurve.circle(5), c2 = EllipseCurve.circle(5, V(3, 0))
			testCurveISInfos(assert, c1, c2, 2)

			const verticalEllipse = new EllipseCurve(V(2, 0), V(1, 1), V(1, 10))
			testCurveISInfos(assert, c1, verticalEllipse, 4)

			const verticalEllipse2 = new EllipseCurve(V(10, 2), V(1, 1), V(1, 10))
			testCurveISInfos(assert, c1, verticalEllipse2, 0)

			const smallEllipse = EllipseCurve.forAB(2, 3)
			testCurveISInfos(assert, c1, smallEllipse, 0)

			const test = new EllipseCurve(V(6, 1, 0), V(3, 1, 0), V(4, 0, 0))
			testCurveISInfos(assert, c1, test, 2)

			const e1 = SemiEllipseCurve.UNIT
			const e2 = SemiEllipseCurve.UNIT.scale(0.5, 0.1, 1).translate(0.5)
			testCurveISInfos(assert, e1, e2, 1)
		},
		'SemiEllipseCurve.isInfosWithBezier2D()'(assert) {
			const ell = SemiEllipseCurve.forAB(3, 1)
			const bez = BezierCurve.graphXY(2, -3, -3, 2, -2, 3)
			testCurveISInfos(assert, ell, bez, 3)
		},
	})
}