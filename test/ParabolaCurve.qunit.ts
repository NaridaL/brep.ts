QUnit.module('ParabolaCurve')
{
	const curve = new ParabolaCurve(V(1, 1), V(4, 1, -2), V(1, 10, 2))
	registerTests({
		'testCurve'(assert) {
			testCurve(assert, curve)
			testCurve(assert, ParabolaCurve.XY)
		},
		'isTsWithPlane'(assert) {
			const plane = new P3(V(2, 7, 1).unit(), 10)
			testISTs(assert, curve, plane, 2)
		},
		'rightAngled'(assert) {
			const curveRA = curve.rightAngled()
			assert.ok(curveRA.isColinearTo(curve))
			NLA.arrayRange(-10, 10, 1).forEach(t => assert.ok(curveRA.containsPoint(curve.at(t))))
		}
	})
}