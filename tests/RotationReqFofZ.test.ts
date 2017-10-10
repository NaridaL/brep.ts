import {suite, test, testParametricSurface} from './manager'

import {M4, TAU} from 'ts3dutils'
import {RotationREqFOfZ} from '..'

const {sin, cos} = Math

suite('RotationREqFOfZ', () => {
	const sinSurface = new RotationREqFOfZ(M4.IDENTITY, z => sin(z) + 2 + 0.2 * z, 0, 2 * TAU, 1, z => cos(z) + 0.2)

	test('testSurface', assert => {
		testParametricSurface(assert, sinSurface)
	})
})