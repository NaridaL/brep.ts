import {lerp, V} from 'ts3dutils'
import {Curve, Edge, P3, ParabolaCurve} from '..'
import {Assert, outputLink, suite, test, testCurve, testISTs, testCurvesColinear} from './manager'

suite('ParabolaCurve', () => {
	const curve = new ParabolaCurve(V(1, 1), V(4, 1, -2), V(1, 10, 2))

	test('testCurve', assert => {
		testCurve(assert, curve)
		testCurve(assert, ParabolaCurve.XY)
	})
	test('isTsWithPlane', assert => {
		const plane = new P3(V(2, 7, 1).unit(), 10)
		testISTs(assert, curve, plane, 2)
	})
	test('rightAngled', assert => {
		const curveRA = curve.rightAngled()
		testCurvesColinear(assert, curve, curveRA)
	})
	test('isColinearTo', assert => {
		testCurvesColinear(assert, ParabolaCurve.XY, ParabolaCurve.XY.scale(2, 4, 1))
	})
})
