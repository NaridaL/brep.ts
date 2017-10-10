import {DEG, M4, V, V3} from 'ts3dutils'
import {BezierCurve, L3, P3, PCurveEdge, PlaneSurface, RotationFace, SemiEllipsoidSurface, StraightEdge,ProjectedCurveSurface} from '..'
import {inDifferentSystems, suite, test, testISTs, testParametricSurface} from './manager'

suite('ProjectedCurveSurface', () => {
	test('ProjectedCurveSurface', inDifferentSystems((assert, m4) => {
		const baseCurve = BezierCurve.graphXY(2, -3, -3, 2)
		const pcs = new ProjectedCurveSurface(baseCurve, V3.Z, undefined, undefined, -100, 100).transform(m4)
		testParametricSurface(assert, pcs)
	}))
	test('Face line intersection test', inDifferentSystems((assert, m4) => {
		const curve = BezierCurve.graphXY(2, -3, -3, 2)
		const edge = PCurveEdge.forCurveAndTs(curve, 0, 1)
		const edges = [
			edge,
			StraightEdge.throughPoints(curve.at(1), curve.at(1).plus(V(0, 0, 10))),
			edge.flipped().transform(M4.translate(0, 0, 10)),
			StraightEdge.throughPoints(curve.at(0).plus(V(0, 0, 10)), curve.at(0))]
		const surface = new ProjectedCurveSurface(curve, V3.Z)
		const face = new RotationFace(surface, edges).transform(m4)
		const line = new L3(V3.Z, V3.X).transform(m4)
		const d = face.intersectsLine(line)
		assert.ok(d)
	}))
	const pcs = new ProjectedCurveSurface(BezierCurve.EX2D, V3.Z, undefined, undefined, -1, 1)
	test('testSurface', assert => {
		testParametricSurface(assert, pcs)
	})
	test('isTsWithSurface', assert => {
		const pcs = new ProjectedCurveSurface(BezierCurve.EX2D, V3.Z, undefined, undefined, -2, 2).scale(0.2, 0.2, 1).rotateX(-90 * DEG)
		const ses = SemiEllipsoidSurface.UNIT
		const pic = ses.isCurvesWithSurface(pcs)[0]
		testISTs(assert, pic, new PlaneSurface(P3.XY), 1)
	})


	test('Face containsPoint', assert => {
		const face = new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(142.87578921496748, -191.46078243076332, 0), V(161.78547089700214, -252.13248349581008, 0), V(284.63214994898954, -163.59789158697575, 0), V(372.40411211189405, -210.3992206435476, 0), -3, 4), V(0, 0, 1), -3, 4, -100, 100), [
			PCurveEdge.forCurveAndTs(
				new BezierCurve(V(142.87578921496748, -191.46078243076332, 0), V(161.78547089700214, -252.13248349581008, 0), V(284.63214994898954, -163.59789158697575, 0), V(372.40411211189405, -210.3992206435476, 0)), 1, 0),
			StraightEdge.throughPoints(V(142.87578921496748, -191.46078243076332, 0), V(142.87578921496748, -191.46078243076332, -100)),
			PCurveEdge.forCurveAndTs(new BezierCurve(V(142.87578921496748, -191.46078243076332, -100), V(161.78547089700214, -252.13248349581008, -100), V(284.63214994898954, -163.59789158697575, -100), V(372.40411211189405, -210.3992206435476, -100)), 0, 1),
			StraightEdge.throughPoints(V(372.40411211189405, -210.3992206435476, -100), V(372.40411211189405, -210.3992206435476, 0))], [])
		const line = new L3(V(1241.5987, -1214.1894, 38.9886), V(-0.6705, 0.7386, -0.0696).unit())
		testISTs(assert, line, face.surface, 3)
	})
})
}