QUnit.module('ParabolaCurve')
{
	function testCurvesColinear(test: Assert, curve1: Curve, curve2: Curve): void {
		test.ok(curve1.isColinearTo(curve2))
		const t = (curve1.tMin + curve1.tMax) / 2
		test.notOk(curve1.translate(curve1.tangentAt(t).getPerpendicular().unit()).isColinearTo(curve2))
		linkB3(test, {edges: [curve1, curve2].map(c => Edge.forCurveAndTs(c))})
		for (let i = 0; i < 10; i++) {
			const t = lerp(curve1.tMin, curve1.tMax, i / 9)
			if (!curve2.containsPoint(curve1.at(t))) {
				test.push(false)
			}
		}
	}
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
			testCurvesColinear(assert, curve, curveRA)
		},
		'isColinearTo'(assert) {
			testCurvesColinear(assert, ParabolaCurve.XY, ParabolaCurve.XY.scale(2, 4, 1))
		},
	})
}