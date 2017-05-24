{
	QUnit.module('RotationREqFOfZ')
	const sinSurface = new RotationREqFOfZ(M4.IDENTITY, z => sin(z) + 2 + 0.2 * z, 0, 2 * TAU, 1, z => cos(z) + 0.2)
	registerTests({
		'testSurface'(assert) {
			testParametricSurface(assert, sinSurface)
		},
	})
}