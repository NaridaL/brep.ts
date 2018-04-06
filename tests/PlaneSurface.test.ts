import { suite, surfaceVolumeAndAreaTests, test, testLoopContainsPoint } from './manager'

import { Edge, Face, L3, P3, PlaneSurface, PointVsFace, SemiEllipseCurve, StraightEdge } from '..'

import { DEG, V, V3 } from 'ts3dutils'

suite('PlaneSurface', () => {
	test('loopContainsPoint', assert => {
		const loop = StraightEdge.chain([V(0, 0), V(10, 0), V(10, 10), V(0, 10)], true)
		assert.equal(new PlaneSurface(P3.XY).loopContainsPoint(loop, V(8, 10)), PointVsFace.ON_EDGE)
	})
	test('loopContainsPoint 2', assert => {
		const loop = [
			new StraightEdge(new L3(V(2, 10, 0), V3.Z), V(2, 10, 3), V(2, 10, 5), 3, 5),
			new StraightEdge(new L3(V(0, 10, 5), V3.X), V(2, 10, 5), V(0, 10, 5), 2, 0),
			new StraightEdge(new L3(V(0, 10, 0), V3.Z), V(0, 10, 5), V(0, 10, 0), 5, 0),
			new StraightEdge(new L3(V(0, 10, 0), V3.X), V(0, 10, 0), V(10, 10, 0), 0, 10),
			new StraightEdge(new L3(V(10, 10, 0), V3.Z), V(10, 10, 0), V(10, 10, 5), 0, 5),
			new StraightEdge(new L3(V(0, 10, 5), V3.X), V(10, 10, 5), V(6, 10, 5), 10, 6),
			new StraightEdge(new L3(V(6, 10, 0), V(0, 0, -1)), V(6, 10, 5), V(6, 10, 3), -5, -3),
			new StraightEdge(new L3(V(0, 10, 3), V(-1, 0, 0)), V(6, 10, 3), V(2, 10, 3), -6, -2),
		]
		const p = V(6, 10, 3)
		testLoopContainsPoint(assert, new PlaneSurface(new P3(V(0, -1, 0), -10)), loop, p, PointVsFace.ON_EDGE)
	})

	const triangleFace = Face.create(new PlaneSurface(P3.XY), StraightEdge.chain([V(1, 1), V(3, 2), V(2, 3)])).rotateX(
		10 * DEG,
	)
	suite('triangleFace', () => surfaceVolumeAndAreaTests(triangleFace))
	suite('triangleFace.translate(2, 2, 2)', () => surfaceVolumeAndAreaTests(triangleFace.translate(2, 2, 2)))
	suite('triangleFace.shearX(2, 2)', () => surfaceVolumeAndAreaTests(triangleFace.shearX(2, 2)))
	suite('triangleFace.foo()', () => surfaceVolumeAndAreaTests(triangleFace.foo()))

	const faceWithEllipses = Face.create(new PlaneSurface(P3.XY), [
		Edge.forCurveAndTs(SemiEllipseCurve.UNIT),
		Edge.forCurveAndTs(SemiEllipseCurve.circleThroughPoints(V3.X.negated(), V(0, -0.5), V3.X)),
	])
	suite('faceWithEllipses', () => surfaceVolumeAndAreaTests(faceWithEllipses))
	suite('faceWithEllipses.foo()', () => surfaceVolumeAndAreaTests(faceWithEllipses.foo()))
})
