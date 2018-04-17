import { suite, test, testCurve, testCurveISInfos, testISTs } from './manager'

import { V, V3 } from 'ts3dutils'
import { BezierCurve, P3, PlaneSurface, EllipseCurve, EllipsoidSurface } from '..'
import { PI, sin } from '../src/math'

suite('EllipseCurve', () => {
	const curve = EllipseCurve.UNIT.shearX(2, 1)

	test('withBounds', assert => {
		const newCurve = curve.withBounds(1, 2)
		assert.equal(newCurve.tMin, 1)
		assert.equal(newCurve.tMax, 2)
	})
	test('testCurve', assert => {
		testCurve(assert, EllipseCurve.UNIT)
		testCurve(assert, curve)
	})
	test('UNIT.shearX(2, 3)', assert => testCurve(assert, EllipseCurve.UNIT.shearX(2, 2)))
	test('isTsWithPlane', assert => {
		const plane = new P3(V(2, 7, 1).unit(), 2)
		testISTs(assert, curve.scale(1, 3, 1), plane, 2)
	})
	//test('rightAngled', assert => {
	//	const curveRA = curve.rightAngled()
	//	assert.ok(curveRA.f1.isPerpendicularTo(curveRA.f2))
	//	assert.ok(curveRA.isColinearTo(curve))
	//	arrayRange(-10, 10, 1).forEach(t => assert.ok(curveRA.containsPoint(curve.at(t))))
	//},
	test('isTsWithSurface(EllipsoidSurface)', assert => {
		const s = EllipsoidSurface.sphere(5)
		const c = new EllipseCurve(V(5, 2), V3.Z.negated(), V(-1, 1.2246467991473532e-16, 0), 0, PI)
		testISTs(assert, c, s, 2)
	})
	test('isTsWithSurface(PlaneSurface)', assert => {
		const c = EllipseCurve.UNIT.translate(1.2, -1)
		const s = new PlaneSurface(P3.ZX)
		testISTs(assert, c, s, 1)
	})
	test('isTsWithSurface(PlaneSurface) 2', assert => {
		const c = EllipseCurve.UNIT
		const s = P3.YZ.translate(0.5, 0)
		testISTs(assert, c, s, 1)
	})
	test('distanceToPoint', assert => {
		const curve = EllipseCurve.forAB(10, 15)
		const p = V(12, 12)
		const closestT = curve.closestTToPoint(p)
		const pDist = curve.at(closestT).distanceTo(p)
		const EPS = 0.001
		assert.push(
			pDist < curve.at(closestT - EPS).distanceTo(p),
			curve.at(closestT - EPS).distanceTo(p),
			'> ' + pDist,
			'' + (pDist - curve.at(closestT - EPS).distanceTo(p)) + 'larger',
		)
		assert.push(
			pDist < curve.at(closestT + EPS).distanceTo(p),
			curve.at(closestT + EPS).distanceTo(p),
			'> ' + pDist,
		)
	})
	test('isColinearTo', assert => {
		assert.ok(EllipseCurve.forAB(1, 2).isColinearTo(EllipseCurve.forAB(1, -2)))
	})
	const c1 = EllipseCurve.semicircle(5)
	test('isInfosWithEllipse', assert => {
		const c1 = EllipseCurve.semicircle(5),
			c2 = EllipseCurve.semicircle(5, V(3, 0))
		testCurveISInfos(assert, c1, c2, 1, 'c1 c2')

		const verticalEllipse = new EllipseCurve(V(2, 0), V(1, 1), V(1, 10))
		testCurveISInfos(assert, c1, verticalEllipse, 2, 'c1 verticalEllipse')

		const verticalEllipse2 = new EllipseCurve(V(10, 2), V(1, 1), V(1, 10))
		testCurveISInfos(assert, c1, verticalEllipse2, 0, 'c1 verticalEllipse2')

		const smallEllipse = EllipseCurve.forAB(2, 3)
		testCurveISInfos(assert, c1, smallEllipse, 0, 'c1 smallEllipse')
	})
	test('c1 test', assert => {
		const test = new EllipseCurve(V(6, 1, 0), V(3, 1, 0), V(4, 0, 0))
		testCurveISInfos(assert, c1, test, 1, 'c1 test')
	})
	test('isInfosWithBezier2D', assert => {
		const ell = EllipseCurve.forAB(3, 1)
		const bez = BezierCurve.graphXY(2, -3, -3, 2, -2, 3)
		testCurveISInfos(assert, ell, bez, 3)
	})
})
