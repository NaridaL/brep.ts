import {suite, test, inDifferentSystems} from './manager'

suite('BezierCurve', () => {
	test('isTsWithSurface(SemiCylinderSurface)', inDifferentSystems((assert, m4) => {
		const bez = BezierCurve.graphXY(2, -3, -3, 2, -2, 3).rotateX(15 * DEG).translate(0, 0, 100).transform(m4)
		const cyl = new SemiCylinderSurface(SemiEllipseCurve.forAB(4, 1).rotateY(10 * DEG), V3.Z, undefined, undefined).transform(m4)
		testISTs(assert, bez, cyl, 3)
	}))

	test('testCurve', assert => {
		const curve = BezierCurve.graphXY(2, -3, -3, 2, -2, 3)//.rotateZ(PI/3)
		testCurve(assert, curve)
	})
	test('pointT', assert => {
		const curve = new BezierCurve(
			V(92.48132002394416, 253.35277539335377, 0),
			V(99.18055157018783, 225.4322156490681, 0),
			V(63.52151476563836, 168.59279980361327, 0),
			V(151.89049176954802, 231.21343792479922, 0))
		const p = V(90.8280915025532, 214.7764313721318, 0)
		assert.ok(eq0(curve.distanceToPoint(p)))
		assert.ok(isFinite(curve.pointT(p)))
	})
	test('pointT 2', assert => {
		const curve = new BezierCurve(V(67.44, 3.02, 0), V(67.42, 2.8200000000000003, 0), V(67.39333333333333, 2.64, 0), V(67.36, 2.48, 0), 0, 1)
		const t = 56.58829486216517
		const p = curve.at(t)
		assert.push(eq(t, curve.pointT2(p)), curve.pointT(p), t)
	})
	test('pointT 3', assert => {
		const curve = new BezierCurve(V(75.07, 17.86, 0), V(75.07, 28.16, 0), V(70.27, 34.5, 0), V(61.44, 34.5, 0), 0, 1)
		const p = curve.p3
		assert.push(eq(1, curve.pointT(p)), curve.pointT(p), 1)
	})
	test('pointT 4', assert => {
		const curve = new BezierCurve(V(11.74, 11.49, 0), V(14.18, 12.74, 0), V(15.39, 14.34, 0), V(15.39, 16.29, 0), 0, 1)
		const p = curve.p3
		assert.push(eq(1, curve.pointT(p)), curve.pointT(p), 1)
	})
	test('distanceToPoint', assert => {
		const curve = BezierCurve.graphXY(0, 0, 0, 1)//.rotateZ(PI/3)
		//        assert.ok(eq2(curve.distanceToPoint(V(0.5, 0)), 1, NLA_PRECISION))

		const curve2 = BezierCurve.graphXY(2, -3, -3, 2)
		const p = V(0.5, 0.2)
		const closestT = curve2.closestTToPoint(p)
		const pDist = curve2.at(closestT).distanceTo(p)
		const EPS = NLA_PRECISION
		assert.push(pDist < curve2.at(closestT - EPS).distanceTo(p), curve2.at(closestT - EPS).distanceTo(p), '> ' + pDist, '' + (pDist - curve2.at(closestT - EPS).distanceTo(p)) + 'larger')
		assert.push(pDist < curve2.at(closestT + EPS).distanceTo(p), curve2.at(closestT + EPS).distanceTo(p), '> ' + pDist)

		const curve3 = BezierCurve.graphXY(2, -3, -3, 2).scale(100, 100, 100)
		const p3 = V(71, -65, 0)
		const closestT3 = curve3.closestTToPoint(p3)
		const pDist3 = curve3.at(closestT3).distanceTo(p3)
		assert.push(pDist3 < curve3.at(closestT3 - EPS).distanceTo(p3), curve3.at(closestT3 - EPS).distanceTo(p3), '> ' + pDist3, '' + (pDist3 - curve3.at(closestT3 - EPS).distanceTo(p3)) + 'larger')
		assert.push(pDist3 < curve3.at(closestT3 + EPS).distanceTo(p3), curve3.at(closestT3 + EPS).distanceTo(p3), '> ' + pDist3)

	})
	test('isPointsWithBezier()', assert => {
		const curve1 = BezierCurve.graphXY(2, -3, -3, 2, -2, 3)
		const curve2 = curve1.transform(M4.rotateLine(V(0.5, 0), V3.Z, PI / 2))
		testCurveISInfos(assert, curve1, curve2, 9)

			// test self-intersections
			;[
				new BezierCurve(V(133, 205, 0), V(33, 240, 0), V(63, 168, 0), V(151, 231, 0)),
				new BezierCurve(V3.O, V(10, 10), V(-9, 10), V3.X)
			].forEach(curve => {
				assert.push(true, undefined, undefined, 'Testing ' + curve.sce)
				const isInfos = curve.isInfosWithBezier(curve, 0, 1, 0, 1)
				assert.equal(isInfos.length, 1)
				isInfos.forEach(info => {
					const p = info.p
					assert.ok(eq0(curve.distanceToPoint(p)), `curve.distanceToPoint(${p}) = ${curve.distanceToPoint(p, -2, 3)}`)
				})
			})
	})

	test('isInfosWithLine', assert => {
		const curve = BezierCurve.graphXY(2, -3, -3, 2, -3, 4)
		const line = new L3(V3.Y, V3.X)
		testCurveISInfos(assert, curve, line, 3)

		const line2 = new L3(V(0, 2, 1), V3.Z)
		testCurveISInfos(assert, curve, line2, 1)

		const line3 = new L3(V3.Z, V3.X)
		testCurveISInfos(assert, curve, line3, 0)
	})
	test('isInfosWithLine 2', assert => {
		const curve = new BezierCurve(V(0.1, 0.39392310120488316, 0.06945927106677233), V(0.1, 0.5027019619056378, 0.08863991913904307), V(0.010456949966158688, 0.5908846518073247, 0.10418890660015839), V(-0.09999999999999999, 0.5908846518073247, 0.10418890660015839), 0, 1)
		const line = new L3(V(0, 0, 0), V(0, 0.9848077530122081, 0.17364817766693033))
		testCurveISInfos(assert, curve, line, 1)
	})
	test('isTsWithPlane', assert => {
		const curve = BezierCurve.graphXY(2, -3, -3, 2, -2, 3)
		const plane = new P3(V(0, 1, 1).unit(), 1)

		testISTs(assert, curve, plane, 3)
	})
	test('isTsWithEllipsoidSurface', assert => {
		//const curve = ParabolaCurve.XY.asBezier().scale(5).translate(0, 1)
		const curve = BezierCurve.graphXY(2, -3, -3, 2, -2, 3)
		const s = EllipsoidSurface.UNIT


		testISTs(assert, curve.translate(0.2), s, 4)
		testISTs(assert, curve, s, 4)
		// TODO
		//testISTs(assert, curve.translate(-0.00635), s, 3)
		//console.log(arrayRange(-0.00640, -0.00630, 0.000005).map(i =>
		// curve.translate(i).isTsWithSurface(s).length))
	})
	test('isInfosWithEllipseCurve', assert => {
		const bc = new BezierCurve(V(2, 0, 1.2185119342628936), V(1, 0, 0.8800363969676455), V(-1, 0, 0.5415608596723974), V(-2, 0, 0.2030853223771491), -0.10000000000000009, 1.1)
		const sec = new SemiEllipseCurve(V3.O, V3.X, V3.Z)
		//const curve = ParabolaCurve.XY.asBezier().scale(5).translate(0, 1)
		testCurveISInfos(assert, bc, sec, 2)
	})
})