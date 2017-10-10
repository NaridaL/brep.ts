import {DEG, V, V3} from 'ts3dutils'
import {
	B2T, ConicSurface, Edge, HyperbolaCurve, L3, P3, PCurveEdge, PointVsFace, SemiEllipseCurve, SemiEllipsoidSurface,
	StraightEdge,
} from '..'
import {linkB2, suite, test, testISCurves, testLoopCCW, testLoopContainsPoint, testParametricSurface} from './manager'

suite('ConicSurface', () => {
	const UCS = ConicSurface.UNIT

	test('testSurface', assert => {
		testParametricSurface(assert, UCS)
		testParametricSurface(assert, ConicSurface.UNIT.scale(2, 2, 1))
		testParametricSurface(assert, new ConicSurface(
			V(2, 0.2, 1.1),
			V(0, 0.6, 0),
			V(0, 0, -2.4),
			V(-12, 0, 0)))
	})
	test('testLoopCCW', assert => {
		const surface = new ConicSurface(V(0, 0, 53.51411369448604), V(198.46477746372744, 0, 0), V(0, 198.46477746372744, 0), V(0, 0, 191.42941531213293)).scale(1 / 200)
		const loop = [
			new StraightEdge(new L3(V(131.35224103228387, 0, 180.2100595549249), V(0.7197488536413841, 0, 0.6942345336281635)), V(131.35224103228387, 0, 180.2100595549249), V(198.46477746372744, 0, 244.94352900661897), 0, 93.24438113642698),
			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 244.94352900661897), V(198.46477746372744, 0, 0), V(0, 198.46477746372744, 0), 0, 3.141592653589793), V(198.46477746372744, 0, 244.94352900661897), V(-198.46477746372744, 2.4304925446444556e-14, 244.94352900661897), 0, 3.141592653589793, null, V(0, 198.46477746372744, 0), V(-2.4304925446444556e-14, -198.46477746372744, 0)),
			new StraightEdge(new L3(V(-131.35224103228387, 1.6086010154101807e-14, 180.2100595549249), V(-0.7197488536413841, 8.814381298018978e-17, 0.6942345336281635)), V(-198.46477746372744, 2.4304925446444556e-14, 244.94352900661897), V(-131.35224103228387, 1.6086010154101807e-14, 180.2100595549249), 93.24438113642698, 0),
			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 180.2100595549249), V(131.35224103228387, 0, 0), V(0, 131.35224103228387, 0), 0, 3.141592653589793), V(-131.35224103228387, 1.6086010154101807e-14, 180.2100595549249), V(131.35224103228387, 0, 180.2100595549249), 3.141592653589793, 0, null, V(1.6086010154101807e-14, 131.35224103228387, 0), V(0, -131.35224103228387, 0)),
		].map(e => e.scale(1 / 200))
		testLoopCCW(assert, surface, Edge.reversePath(loop))
		testLoopCCW(assert, surface.flipped(), loop)
	})
	test('isCoplanarTo', assert => {
		assert.ok(UCS.matrix.isIdentity(), 'UCS.matrix.isIdentity()')
		assert.v3like(UCS.pSTFunc()(0, 3), V(3, 0, 3))
		const ellipseAtZ3 = SemiEllipseCurve.UNIT.scale(3, 3, 3).translate(0, 0, 3)
		const planeAtZ3 = P3.XY.translate(0, 0, 3)


		const CS2 = ConicSurface.UNIT.scale(2, 2, 1)
		assert.notOk(CS2.isCoplanarTo(UCS))
		assert.notOk(UCS.isCoplanarTo(CS2))
		const ell1 = UCS.isCurvesWithPlane(new P3(V(2, 3, 10).unit(), 10))[0]
		assert.ok(UCS.containsEllipse(ell1), 'UCS.containsEllipse(ell1)')
		const ell2 = UCS.isCurvesWithPlane(new P3(V(1, 1, 2).unit(), 4))[0]
		const ell1Cone = ConicSurface.atApexThroughEllipse(V3.O, ell1, 1)
		const ell2Cone = ConicSurface.atApexThroughEllipse(V3.O, ell2, 1)
		assert.ok(UCS.isCoplanarTo(ell1Cone))
		assert.ok(UCS.isCoplanarTo(ell2Cone))
		assert.ok(ell1Cone.isCoplanarTo(ell2Cone))
		assert.ok(ell2Cone.isCoplanarTo(ell1Cone))
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
		const a = ConicSurface.UNIT
			.scale(0.05, 0.2)
			.rotateZ(90 * DEG)
			.rotateY(-90 * DEG)
			.translate(2, 0.2, 1.1).flipped()
		const cone = new ConicSurface(
			V(2, 0.2, 1.1),
			V(0, 0.6, 0),
			V(0, 0, -2.4),
			V(-12, 0, 0))
		const sphere = new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V(0, 0, -1))
		const b = SemiEllipsoidSurface.UNIT
		testISCurves(assert, a, b, 2)
		testISCurves(assert, cone, sphere, 2)
	})
	test('isCurvesWithEllipsoid 2', assert => {

		const cone = new ConicSurface(V(2, 0.2, 0.7), V(2.2496396739927868e-33, 0.6000000000000001, 3.67394039744206e-17), V(-1.469576158976824e-16, 1.469576158976824e-16, -2.4000000000000004), V(-12, 0, 7.347880794884119e-16))
		const sphere = new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V(0, 0, -1))


		testISCurves(assert, cone, sphere, 2)
	})
	test('containsParabola', assert => {
		const pb = UCS.isCurvesWithPlane(new P3(V(1, 0, 1).unit(), 4))[0]
		assert.ok(UCS.containsParabola(pb))

		const c2 = UCS.shearX(2, 3)
		const pb2 = c2.isCurvesWithPlane(new P3(V(1, 0, 1).unit(), 4).shearX(2, 3))[0]
		assert.ok(c2.containsParabola(pb2))
	})
	test('containsHyperbola', assert => {
		const s = new ConicSurface(V(-242.1625189124994, 38.960257711878945, 0), V(197.87979681325515, -15.226749714620981, 2.4304925446444556e-14), V(2.4233285978328154e-14, -1.8647390299428456e-15, -198.46477746372744), V(14.686977871964286, 190.86517159433123, 0))
		const c = new HyperbolaCurve(V(-242.16251891249937, 38.960257711878945, -100.00000000000003), V(7.400294429901329, 96.17080372320217, 0), V(-99.70524711843181, 7.672268051394617, -1.8369701987210304e-14))

		linkB2(assert, `viewer.html?mesh=[${s.sce}.toMesh()]&edges=[Edge.forCurveAndTs(${c.sce})]`)
		assert.ok(s.containsCurve(c))
	})
	test('containsPoint', assert => {
		const face = B2T.cone(1, 1, PI).faces.find(face => face.surface instanceof ConicSurface)
		testLoopContainsPoint(assert, face.surface, face.contour, V(-0.2, 0, 0.2), PointVsFace.ON_EDGE)
	})
})