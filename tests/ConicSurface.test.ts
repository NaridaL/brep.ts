import {
	suite,
	suiteSurface,
	surfaceVolumeAndAreaTests,
	test,
	testContainsCurve,
	testISCurves,
	testLoopCCW,
	testLoopContainsPoint,
} from './manager'

import { DEG, V, V3 } from 'ts3dutils'
import {
	B2T,
	ConicSurface,
	Edge,
	HyperbolaCurve,
	L3,
	P3,
	ParabolaCurve,
	PCurveEdge,
	PointVsFace,
	RotationFace,
	SemiEllipseCurve,
	SemiEllipsoidSurface,
	StraightEdge,
} from '..'
import { PI } from '../src/math'

suite('ConicSurface', () => {
	const UCS = ConicSurface.UNIT

	suite('UNIT', () => suiteSurface(ConicSurface.UNIT))
	suite('UNIT.scale(2, 2, 1)', () => suiteSurface(ConicSurface.UNIT.scale(2, 2, 1)))
	suite('weird', () => suiteSurface(new ConicSurface(V(2, 0.2, 1.1), V(0, 0.6, 0), V(0, 0, -2.4), V(-12, 0, 0))))

	const testFace = new RotationFace(
		new ConicSurface(V3.Z, V(-1, 0, 0), V3.Y, V(0, 0, -1)),
		[
			new PCurveEdge(
				new HyperbolaCurve(
					V(0.10792279653395696, 0.0905579787672639, 1),
					V(0, 0, -0.1408832052805518),
					V(0.0905579787672639, -0.10792279653395696, 0),
					-7,
					7,
				),
				V(-0.4634542598189221, 0.7714986384017117, 0.1),
				V(0.1, 0.1, 0.8585786437626904),
				-2.5414277085137025,
				-0.08737743596203365,
				undefined,
				V(0.5785088487178852, -0.6894399988070802, 0.8889049006895381),
				V(0.09090389553440875, -0.10833504408394042, 0.012325683343243887),
				'genseg17',
			),
			new PCurveEdge(
				new HyperbolaCurve(
					V(-0.00792279653395693, 0.009442021232736126, 1),
					V(0, 0, -0.012325683343243885),
					V(0.009442021232736126, 0.00792279653395693, 0),
					-7,
					7,
				),
				V(0.1, 0.1, 0.8585786437626904),
				V(0.6814525440390486, 0.5878966152502569, 0.1),
				3.131301331471644,
				4.983809888872043,
				undefined,
				V(0.10833504408394039, 0.09090389553440875, -0.14088320528055173),
				V(0.6894399988070802, 0.5785088487178854, -0.8999155946699233),
				'genseg18',
			),
			new PCurveEdge(
				new SemiEllipseCurve(V(0, 0, 0.09999999999999998), V(0.9, 0, 0), V(0, 0.9, 0), 0, 3.141592653589793),
				V(0.6814525440390486, 0.5878966152502569, 0.1),
				V(-0.4634542598189221, 0.7714986384017117, 0.1),
				0.7118273326574678,
				2.1117446875459924,
				undefined,
				V(-0.5878966152502568, 0.6814525440390486, 0),
				V(-0.7714986384017115, -0.4634542598189224, 0),
				'genseg19',
			),
		],
		[],
	)
	suite('testFace surface', () => suiteSurface(testFace.surface as ConicSurface))
	suite('testFace', () => surfaceVolumeAndAreaTests(testFace))
	suite('testFace.scale(2)', () => surfaceVolumeAndAreaTests(testFace.scale(2)))
	suite('testFace.shearX(2, 2)', () => surfaceVolumeAndAreaTests(testFace.shearX(2, 2)))
	suite('testFace.foo()', () => surfaceVolumeAndAreaTests(testFace.foo()))
	test('testLoopCCW', assert => {
		const surface = new ConicSurface(
			V(0, 0, 53.51411369448604),
			V(198.46477746372744, 0, 0),
			V(0, 198.46477746372744, 0),
			V(0, 0, 191.42941531213293),
		).scale(1 / 200)
		const loop = [
			new StraightEdge(
				new L3(V(131.35224103228387, 0, 180.2100595549249), V(0.7197488536413841, 0, 0.6942345336281635)),
				V(131.35224103228387, 0, 180.2100595549249),
				V(198.46477746372744, 0, 244.94352900661897),
				0,
				93.24438113642698,
			),
			new PCurveEdge(
				new SemiEllipseCurve(
					V(0, 0, 244.94352900661897),
					V(198.46477746372744, 0, 0),
					V(0, 198.46477746372744, 0),
					0,
					3.141592653589793,
				),
				V(198.46477746372744, 0, 244.94352900661897),
				V(-198.46477746372744, 2.4304925446444556e-14, 244.94352900661897),
				0,
				3.141592653589793,
				null,
				V(0, 198.46477746372744, 0),
				V(-2.4304925446444556e-14, -198.46477746372744, 0),
			),
			new StraightEdge(
				new L3(
					V(-131.35224103228387, 1.6086010154101807e-14, 180.2100595549249),
					V(-0.7197488536413841, 8.814381298018978e-17, 0.6942345336281635),
				),
				V(-198.46477746372744, 2.4304925446444556e-14, 244.94352900661897),
				V(-131.35224103228387, 1.6086010154101807e-14, 180.2100595549249),
				93.24438113642698,
				0,
			),
			new PCurveEdge(
				new SemiEllipseCurve(
					V(0, 0, 180.2100595549249),
					V(131.35224103228387, 0, 0),
					V(0, 131.35224103228387, 0),
					0,
					3.141592653589793,
				),
				V(-131.35224103228387, 1.6086010154101807e-14, 180.2100595549249),
				V(131.35224103228387, 0, 180.2100595549249),
				3.141592653589793,
				0,
				null,
				V(1.6086010154101807e-14, 131.35224103228387, 0),
				V(0, -131.35224103228387, 0),
			),
		].map(e => e.scale(1 / 200))
		testLoopCCW(assert, surface, Edge.reversePath(loop))
		testLoopCCW(assert, surface.flipped(), loop)
	})
	test('isCoplanarTo', assert => {
		const unitCone = ConicSurface.UNIT
		assert.ok(unitCone.matrix.isIdentity(), 'UCS.matrix.isIdentity()')
		assert.v3like(unitCone.pSTFunc()(0, 3), V(3, 0, 3))
		const ellipseAtZ3 = SemiEllipseCurve.UNIT.scale(3, 3, 3).translate(0, 0, 3)
		const planeAtZ3 = P3.XY.translate(0, 0, 3)
		const issAtZ3 = unitCone.isCurvesWithPlane(planeAtZ3)
		assert.equal(issAtZ3.length, 1)
		assert.ok(ellipseAtZ3.isColinearTo(issAtZ3[0]))
		assert.ok(unitCone.containsEllipse(ellipseAtZ3))

		const scaledUnit = ConicSurface.UNIT.scale(2, 2, 1)
		assert.notOk(scaledUnit.isCoplanarTo(unitCone))
		assert.notOk(unitCone.isCoplanarTo(scaledUnit))
		const ell1 = unitCone.isCurvesWithPlane(new P3(V(2, 3, 10).unit(), 10))[0] as SemiEllipseCurve
		assert.ok(unitCone.containsEllipse(ell1), 'UCS.containsEllipse(ell1)')
		const ell2 = unitCone.isCurvesWithPlane(new P3(V(1, 1, 2).unit(), 4))[0] as SemiEllipseCurve
		const ell1Cone = ConicSurface.atApexThroughEllipse(V3.O, ell1)
		const ell2Cone = ConicSurface.atApexThroughEllipse(V3.O, ell2)
		console.log(ell1Cone)
		assert.ok(unitCone.isCoplanarTo(ell1Cone))
		assert.ok(unitCone.isCoplanarTo(ell2Cone))
		assert.ok(ell1Cone.isCoplanarTo(ell2Cone))
		assert.ok(ell2Cone.isCoplanarTo(ell1Cone))
		assert.ok(ell1Cone.foo().isCoplanarTo(ell2Cone.foo()))
	})
	test('isCurvesWithPlane', assert => {
		testISCurves(assert, UCS, new P3(V(1, 1, 2).unit(), 4), 2)
		testISCurves(assert, UCS, P3.XY.translate(0, 0, 3), 1)
		testISCurves(assert, UCS, P3.XY.translate(0, 0, 3).flipped(), 1)
	})
	test('isCurvesWithPlane 2', assert => {
		testISCurves(assert, ConicSurface.UNIT, P3.ZX, 2)
		testISCurves(assert, ConicSurface.UNIT, P3.ZX.flipped(), 2)
		testISCurves(assert, ConicSurface.UNIT, P3.YZ, 2)
		testISCurves(assert, ConicSurface.UNIT, P3.YZ.flipped(), 2)
		testISCurves(assert, ConicSurface.UNIT, new P3(V(1, 0, 1).unit(), 4), 1)
	})
	test('isCurvesWithPlane hyperbolas', assert => {
		const plane = new P3(V(2, 0, -1).unit(), 1)
		testISCurves(assert, UCS, plane, 1)
		testISCurves(assert, UCS, plane.flipped(), 1)
	})
	test('isCurvesWithEllipsoid', assert => {
		const a = ConicSurface.UNIT.scale(0.05, 0.2)
			.rotateZ(90 * DEG)
			.rotateY(-90 * DEG)
			.translate(2, 0.2, 1.1)
			.flipped()
		const cone = new ConicSurface(V(2, 0.2, 1.1), V(0, 0.6, 0), V(0, 0, -2.4), V(-12, 0, 0))
		const sphere = new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V(0, 0, -1))
		const b = SemiEllipsoidSurface.UNIT
		testISCurves(assert, a, b, 2)
		testISCurves(assert, cone, sphere, 2)
	})
	test('isCurvesWithEllipsoid 2', assert => {
		const cone = new ConicSurface(
			V(2, 0.2, 0.7),
			V(2.2496396739927868e-33, 0.6000000000000001, 3.67394039744206e-17),
			V(-1.469576158976824e-16, 1.469576158976824e-16, -2.4000000000000004),
			V(-12, 0, 7.347880794884119e-16),
		)
		const sphere = new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V(0, 0, -1))

		testISCurves(assert, cone, sphere, 2)
	})
	test('containsParabola', assert => {
		const pb = UCS.isCurvesWithPlane(new P3(V(1, 0, 1).unit(), 4))[0] as ParabolaCurve
		assert.ok(pb instanceof ParabolaCurve)
		testContainsCurve(assert, UCS, pb)

		const c2 = UCS.shearX(2, 3)
		const pb2 = c2.isCurvesWithPlane(new P3(V(1, 0, 1).unit(), 4).shearX(2, 3))[0] as ParabolaCurve
		assert.ok(pb2 instanceof ParabolaCurve)
		testContainsCurve(assert, c2, pb2)
	})
	test('containsHyperbola', assert => {
		const s = new ConicSurface(
			V(-242.1625189124994, 38.960257711878945, 0),
			V(197.87979681325515, -15.226749714620981, 2.4304925446444556e-14),
			V(2.4233285978328154e-14, -1.8647390299428456e-15, -198.46477746372744),
			V(14.686977871964286, 190.86517159433123, 0),
		)
		const c = new HyperbolaCurve(
			V(-242.16251891249937, 38.960257711878945, -100.00000000000003),
			V(7.400294429901329, 96.17080372320217, 0),
			V(-99.70524711843181, 7.672268051394617, -1.8369701987210304e-14),
		)

		testContainsCurve(assert, s, c)
	})
	test('containsPoint', assert => {
		const face = B2T.cone(1, 1, PI).faces.find(face => face.surface instanceof ConicSurface)
		testLoopContainsPoint(assert, face.surface, face.contour, V(-0.2, 0, 0.2), PointVsFace.ON_EDGE)
	})
})
