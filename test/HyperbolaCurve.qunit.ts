QUnit.module('HyperbolaCurve')

registerTests({
	'HyperbolaCurve'(assert) {
		const hb = HyperbolaCurve.UNIT
		testCurve(assert, hb)

		const hbSheared = hb.shearedX(2, 3)
		assert.notOk(hbSheared.isOrthogonal())
		const hbScaledRA = hbSheared.rightAngled()
		assert.ok(hbScaledRA.isOrthogonal(), 'hbScaledRA.isOrthogonal()')
		//TODO:assert.ok(hbSheared.isColinearTo(hbScaledRA))
		testCurve(assert, hbScaledRA)

		assert.deepEqual(intersectionUnitHyperbolaLine(1, 0, 2), {x1: 2, y1: sqrt(3), x2: 2, y2: -sqrt(3)})
	},
})