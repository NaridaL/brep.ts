import {
	outputLink,
	rotateEdge,
	suite,
	suiteSurface,
	surfaceVolumeAndAreaTests,
	test,
	testISCurves,
	testISTs,
	testLoopContainsPoint,
    skip
} from './manager'

import { DEG, M4, raddd, TAU, V, V3 } from 'ts3dutils'
import {
	B2T,
	Edge,
	EllipseCurve,
	Face,
	HyperbolaCurve,
	L3,
	P3,
	PlaneSurface,
	rotateCurve,
	RotatedCurveSurface,
	StraightEdge,
} from '..'
import { cos, PI, sin } from '../src/math'

suite('NURBSSurface', () => {
	const baseCurve = EllipseCurve.forAB(2, 2)
		.rotateZ(-20 * DEG)
		.translate(4, 3)
		.rotateX(90 * DEG)
	const torusSurface = rotateCurve(baseCurve, undefined, undefined, 100 * DEG, false) as RotatedCurveSurface
	const s = torusSurface.asNURBSSurface().transform4(M4.perspective(45, 1, 1, 5).times(M4.rotateX(20*DEG).translate(0, 0, -7)))
	test('b', assert => {
		console.log(s.pUV(0.1, 0.1))
		console.log(s.sce)
	})
	suite('a', () => suiteSurface(s))
})
// 122.96083177519078
