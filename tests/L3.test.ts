import { suite, test } from './manager'

import { V, V3 } from 'ts3dutils'
import { L3 } from '..'

suite('L3', () => {
	test('isInfosWithLine', assert => {
		const l1 = new L3(V3.Y, V3.X)
		const res = l1.isInfosWithLine(V3.X, V(0, 2))[0]
		assert.equal(res.tThis, 1)
		assert.equal(res.tOther, 0.5)
		assert.v3like(res.p, V(1, 1))
	})
})
