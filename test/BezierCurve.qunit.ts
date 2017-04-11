{
	QUnit.module('BezierCurve')
	QUnit.test('testCurve', function(assert) {
		const curve = BezierCurve.graphXY(2, -3, -3, 2, -2, 3)//.rotateZ(PI/3)
		testCurve(assert, curve)
	})
	QUnit.test('BezierCurve.pointT', function(assert) {
		let curve = new BezierCurve(
			V(92.48132002394416, 253.35277539335377, 0),
			V(99.18055157018783, 225.4322156490681, 0),
			V(63.52151476563836, 168.59279980361327, 0),
			V(151.89049176954802, 231.21343792479922, 0))
		let p = V(90.8280915025532, 214.7764313721318, 0)
		assert.ok(NLA.eq0(curve.distanceToPoint(p)))
		assert.ok(isFinite(curve.pointT(p)))
	})
	QUnit.test('BezierCurve.pointT 2', function(assert) {
		const curve = new BezierCurve(V(67.44, 3.02, 0), V(67.42, 2.8200000000000003, 0), V(67.39333333333333, 2.64, 0), V(67.36, 2.48, 0), 0, 1)
		const t = 56.58829486216517
		const p = curve.at(t)
		assert.push(eq(t, curve.pointT2(p)), curve.pointT(p), t)
	})
	QUnit.test('BezierCurve.pointT 3', function(assert) {
		const curve = new BezierCurve(V(75.07, 17.86, 0), V(75.07, 28.16, 0), V(70.27, 34.5, 0), V(61.44, 34.5, 0), 0, 1)
		const p = curve.p3
		assert.push(eq(1, curve.pointT(p)), curve.pointT(p), 1)
	})
	QUnit.test('BezierCurve.pointT 4', function(assert) {
		const curve = new BezierCurve(V(11.74, 11.49, 0), V(14.18, 12.74, 0), V(15.39, 14.34, 0), V(15.39, 16.29, 0), 0, 1)
		const p = curve.p3
		assert.push(eq(1, curve.pointT(p)), curve.pointT(p), 1)
	})
	QUnit.test('BezierCurve.distanceToPoint', function(assert) {
		let curve = BezierCurve.graphXY(0, 0, 0, 1)//.rotateZ(PI/3)
//        assert.ok(NLA.eq2(curve.distanceToPoint(V(0.5, 0)), 1, NLA_PRECISION))

		let curve2 = BezierCurve.graphXY(2, -3, -3, 2)
		let p = V(0.5, 0.2)
		let closestT = curve2.closestTToPoint(p)
		let pDist = curve2.at(closestT).distanceTo(p)
		const EPS = NLA_PRECISION
		assert.push(pDist < curve2.at(closestT - EPS).distanceTo(p), curve2.at(closestT - EPS).distanceTo(p), '> ' + pDist, '' + (pDist - curve2.at(closestT - EPS).distanceTo(p)) + 'larger')
		assert.push(pDist < curve2.at(closestT + EPS).distanceTo(p), curve2.at(closestT + EPS).distanceTo(p), '> ' + pDist)

		let curve3 = BezierCurve.graphXY(2, -3, -3, 2).scale(100, 100, 100)
		let p3 = V(71, -65, 0)
		let closestT3 = curve3.closestTToPoint(p3)
		let pDist3 = curve3.at(closestT3).distanceTo(p3)
		assert.push(pDist3 < curve3.at(closestT3 - EPS).distanceTo(p3), curve3.at(closestT3 - EPS).distanceTo(p3), '> ' + pDist3, '' + (pDist3 - curve3.at(closestT3 - EPS).distanceTo(p3)) + 'larger')
		assert.push(pDist3 < curve3.at(closestT3 + EPS).distanceTo(p3), curve3.at(closestT3 + EPS).distanceTo(p3), '> ' + pDist3)

	})
	QUnit.test('BezierCurve.isPointsWithBezier()', function(assert) {
		let curve1 = BezierCurve.graphXY(2, -3, -3, 2, -2, 3)
		let curve2 = curve1.transform(M4.rotationLine(V(0.5, 0), V3.Z, PI / 2))
		let isInfos = curve1.isInfosWithBezier(curve2)
		assert.equal(isInfos.length, 9)
		console.log(isInfos.map(SCE))
		isInfos.forEach(info => {
			let p = info.p
			assert.ok(curve1.containsPoint(p), `curve1.distanceToPoint(${p}) = ${curve1.distanceToPoint(p, -2, 3)}`)
			assert.ok(curve2.containsPoint(p), `curve2.distanceToPoint(${p}) = ${curve2.distanceToPoint(p, -2, 3)}`)
		})

		// test self-intersections
		;[  new BezierCurve(V(133, 205, 0), V(33, 240, 0), V(63, 168, 0), V(151, 231, 0)),
			new BezierCurve(V3.O, V(10, 10), V(-9, 10), V3.X)].forEach(curve => {
			assert.push(true, undefined, undefined, 'Testing ' + curve.sce)
			const isInfos = curve.isInfosWithBezier(curve, 0, 1, 0, 1)
			assert.equal(isInfos.length, 1)
			console.log(isInfos.map(SCE))
			isInfos.forEach(info => {
				let p = info.p
				assert.ok(NLA.eq0(curve.distanceToPoint(p)), `curve.distanceToPoint(${p}) = ${curve.distanceToPoint(p, -2, 3)}`)
			})
		})
	})
	registerTests({
		'isInfosWithLine'(assert) {
			console.log(solveCubicReal2(1, 0, 0, 0))
			let curve = BezierCurve.graphXY(2, -3, -3, 2, -3, 4)
			let line = new L3(V3.Y, V3.X)
			let isInfos = curve.isInfosWithLine(line.anchor, line.dir1)
			assert.equal(isInfos.length, 3)
			isInfos.forEach(info => {
				let p = info.p
				assert.ok(line.containsPoint(p))
				assert.ok(curve.containsPoint(p))
			})


			let line2 = new L3(V(0, 2, 1), V3.Z)
			let isInfos2 = curve.isInfosWithLine(line2.anchor, line2.dir1)
			assert.equal(isInfos2.length, 1)
			assert.ok(V(0, 2, 0).like(isInfos2[0].p))


			let line3 = new L3(V3.Z, V3.X)
			assert.equal(curve.isInfosWithLine(line3.anchor, line3.dir1).length, 0)

		},
		'isTsWithPlane'(assert) {
			const curve = BezierCurve.graphXY(2, -3, -3, 2, -2, 3)
			const plane = new P3(V(0, 1, 1).unit(), 1)

			testISTs(assert, curve, plane, 3)
		},
		'isTsWithEllipsoidSurface'(assert) {
			//const curve = ParabolaCurve.XY.asBezier().scale(5).translate(0, 1)
			const curve = BezierCurve.graphXY(2, -3, -3, 2, -2, 3)
			const s = EllipsoidSurface.UNIT


			testISTs(assert, curve.translate(0.2), s, 4)
			testISTs(assert, curve, s, 4)
			// TODO
			//testISTs(assert, curve.translate(-0.00635), s, 3)
			//console.log(NLA.arrayRange(-0.00640, -0.00630, 0.000005).map(i => curve.translate(i).isTsWithSurface(s).length))
		},
	})
}