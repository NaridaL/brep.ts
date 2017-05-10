///<reference path="helperfunctions.ts"/>
QUnit.module('HyperbolaCurve')

registerTests({
	'testCurve'(assert) {
		testCurve(assert, HyperbolaCurve.XY)
	},
	'HyperbolaCurve'(assert) {
		const hbSheared = HyperbolaCurve.XY.shearedX(2, 3)
		assert.notOk(hbSheared.isOrthogonal())
		const hbShearedRA = hbSheared.rightAngled()
		assert.ok(hbShearedRA.isOrthogonal(), 'hbShearedRA.isOrthogonal()')
		assert.ok(hbSheared.isColinearTo(hbShearedRA))
		testCurve(assert, hbShearedRA)

		assert.deepEqual(intersectionUnitHyperbolaLine(1, 0, 2), {x1: 2, y1: sqrt(3), x2: 2, y2: -sqrt(3)})
	},
	'isTsWithPlane'(assert) {
		testISTs(assert, HyperbolaCurve.XY, P3.YZ, 0)
		testISTs(assert, HyperbolaCurve.XY, P3.YZ.translate(1), 1)
		testISTs(assert, HyperbolaCurve.XY, P3.YZ.translate(2), 2)
		testISTs(assert, HyperbolaCurve.XY, new P3(V(1, 2).unit(), 2), 1)
		testISTs(assert, HyperbolaCurve.XY, new P3(V(1, 2).unit(), 2).flipped(), 1)
		testISTs(assert, HyperbolaCurve.XY, new P3(V(1, -2).unit(), 2), 1)

		testISTs(assert, HyperbolaCurve.XY, new P3(V(1, 1).unit(), 2), 1)
		testISTs(assert, HyperbolaCurve.XY, new P3(V(1, 1).unit(), 2).flipped(), 1)

		testISTs(assert, HyperbolaCurve.XY, new P3(V(2, 1).unit(), 2), 2)
		testISTs(assert, HyperbolaCurve.XY, new P3(V(2, 1).unit(), 2).flipped(), 2)
		testISTs(assert, HyperbolaCurve.XY, new P3(V(2, 1).unit(), 0.85), 2)
		testISTs(assert, HyperbolaCurve.XY, new P3(V(2, 1).unit(), 0.5), 0)
	},
	'isTsWithPlane no IS with planes X < 0'(assert) {
		testISTs(assert, HyperbolaCurve.XY, P3.YZ, 0)
		testISTs(assert, HyperbolaCurve.XY, P3.YZ.translate(-2), 0)

		testISTs(assert, HyperbolaCurve.XY, new P3(V(1, 2).unit(), -2).flipped(), 1)
		testISTs(assert, HyperbolaCurve.XY, new P3(V(1, -2).unit(), -2), 1)

		testISTs(assert, HyperbolaCurve.XY, new P3(V(1, 1).unit(), -2), 0)
		testISTs(assert, HyperbolaCurve.XY, new P3(V(1, 1).unit(), 0), 0)

		testISTs(assert, HyperbolaCurve.XY, new P3(V(2, 1).unit(), -2), 0)
		testISTs(assert, HyperbolaCurve.XY, new P3(V(2, -1).unit(), -2), 0)
	},
})