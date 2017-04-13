{
	QUnit.module('SemiEllipsoidSurface')
	const ses2 = SemiEllipsoidSurface.UNIT.scale(2)
	QUnit.test('testSurface', function (assert) {
		testParametricSurface(assert, ses2)
	})
	QUnit.test('testSurface', function (assert) {
		testISCurves(assert, SemiEllipsoidSurface.UNIT, new PlaneSurface(new P3(V(-1.249000902703301e-16, 1, 0), 0.11006944444444443)), 2)
	})
	QUnit.test('getAABB', function (assert) {
		assert.V3ArraysLike(ses2.getExtremePoints(), [V(0,2,0)])
		assert.V3ArraysLike(ses2.rotateZ(30 * DEG).getExtremePoints(), [V(-2,0,0),V(0,2,0)])
	})
	registerTests({

		'SemiEllipsoidSurface.loopContainsPoint'(assert) {
			const s = new SemiEllipsoidSurface(V3.O, V(5, 0, 0), V(0, 5, 0), V(0, 0, 5))
			const loop = [
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0)), V(0, 0, 5), V(0, 0, -5), PI, 0, null, V(5, 0, 0), V(-5, 0, 0)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(-5, 0, 0)), V(0, 0, -5), V(0, 0, 5), 0, PI, null, V(-5, 0, 0), V(5, 0, 0))]
			const p = V(5, 0, 0)
			testLoopContainsPoint(assert, s, loop, p, PointVsFace.ON_EDGE)
		},
		'SemiEllipsoidSurface.loopContainsPoint 3'(assert) {
			const s = new SemiEllipsoidSurface(V3.O, V(-5, 6.123233995736766e-16, 0), V(-6.123233995736766e-16, -5, 0), V(0, 0, 5))
			const loop = [
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0), 0, 3.141592653589793), V(0, 0, -5), V(1.1291713066130296, 0, -4.87082869338697), 0, 0.22779933669791175, null, V(5, 0, 0), V(4.87082869338697, 0, 1.1291713066130296)),
				new PCurveEdge(new SemiEllipseCurve(V(1.411764705882352, -2.1176470588235303, -1.4117647058823533), V(2.0917534572581817, 2.7890046096775745, -2.091753457258183), V(-2.874840149008801, 3.5206637865439285e-16, -2.874840149008799), 0.7085839061491113, 3.141592653589793), V(1.1291713066130291, 0, -4.87082869338697), V(0, -0.6994424542963525, -4.950836318555471), 0.7085839061491117, 1.0373562345961393, null, V(-3.544048444786543, -1.8149704259460577, -0.821592805867452), V(-3.262983117260863, -2.401508361946856, 0.3392794256594245)),
				new PCurveEdge(new SemiEllipseCurve(V(-1.4117647058823535, -2.117647058823529, -1.4117647058823533), V(2.0917534572581813, -2.789004609677577, 2.0917534572581813), V(2.8748401490088, -3.520663786543927e-16, -2.8748401490088007), 0, 2.43300874744068), V(0, -0.6994424542963525, -4.950836318555471), V(-1.1291713066130296, 1.382836026332681e-16, -4.87082869338697), 2.1042364189936533, 2.4330087474406805, null, V(-3.2629831172608608, 2.401508361946859, -0.33927942565942426), V(-3.5440484447865406, 1.8149704259460617, 0.8215928058674513)),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, 5), V(-5, 0, 0), 0, 3.141592653589793), V(-1.1291713066130296, 0, -4.87082869338697), V(0, 0, -5), 2.9137933168918817, 3.141592653589793, null, V(4.870828693386971, 0, -1.1291713066130287), V(5, 0, -6.123233995736766e-16))]
			const p = V(-4.999999999999999, 0, 0)
			testLoopContainsPoint(assert, s, loop, p, PointVsFace.OUTSIDE)
		},
		'SemiEllipsoidSurface.intersect SES'(assert) {
			const a = SemiEllipsoidSurface.sphere(5)
			const b = SemiEllipsoidSurface.sphere(1).rotateAB(V3.Y, V3.X.negated()).translate(5,2)
			testISCurves(assert, a, b, 2)
			testISCurves(assert, a, b.flipped(), 2)
			testISCurves(assert, a.flipped(), b, 2)
			testISCurves(assert, a.flipped(), b.flipped(), 2)

			testISCurves(assert, b, a, 2)
			testISCurves(assert, SemiEllipsoidSurface.sphere(1), SemiEllipsoidSurface.sphere(2).translate(-1.2), 1)
		},
		'SemiEllipsoidSurface.isCurvesWithProjectedCurveSurface is curves both y > 0'(assert) {
			const s1 = SemiEllipsoidSurface.UNIT
			const s2 = new ProjectedCurveSurface(BezierCurve.EX2D, V3.Z, undefined, undefined, -2, 2)
			testISCurves(assert, s1, s2.translate(0.2).rotateZ(90 * DEG).rotateX(10 * DEG), 2)
		},
		'SemiEllipsoidSurface.isCurvesWithProjectedCurveSurface is curves both cross y = 0'(assert) {
			const s1 = SemiEllipsoidSurface.UNIT
			const s2 = new ProjectedCurveSurface(BezierCurve.EX2D, V3.Z, undefined, undefined, -2, 2)
			testISCurves(assert, s1, s2.translate(0.2), 2)
		},
		'SemiEllipsoidSurface.isCurvesWithProjectedCurveSurface is curves both y < 0'(assert) {
			const s1 = SemiEllipsoidSurface.UNIT
			const s2 = new ProjectedCurveSurface(BezierCurve.EX2D, V3.Z, undefined, undefined, -2, 2)
			testISCurves(assert, s1, s2.translate(0.2).rotateZ(-90 * DEG), 0)
		},
		'SemiEllipsoidSurface.isCurvesWithProjectedCurveSurface one isCurve cross y = 0 twice, one isCurve y > 0'(assert) {
			const s1 = SemiEllipsoidSurface.UNIT
			const s2 = new ProjectedCurveSurface(BezierCurve.EX2D, V3.Z, undefined, undefined, -2, 2)
			testISCurves(assert, s1, s2.translate(0.2).rotateZ(90 * DEG).rotateX(80 * DEG), 3)
		},
		'loopCCW'(assert) {
			const s1 = SemiEllipsoidSurface.UNIT
			const loop =[
				new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), V(-0.03286362438479655, -2.6020852139652106e-18, 0.9994598452125503), -1), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), 0, 29, null, V(-0.003388613668575339, 0, 0.00016274305244215148), V(-0.002184873829710333, -0.0006707444983087241, -0.00007184167849434948)),new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), V(0, 0, 1), 3.108723110778215, 3.141592653589793, null, V(0.9994598452125503, -1.2239853003158589e-16, 0.032863624384797126), V(1, -1.2246467991473532e-16, 0)),new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0, 0, 1), V(0.11093749999999993, 0, 0.9938273849586506), 3.141592653589793, 3.030426330354509, null, V(1, 0, 0), V(0.9938273849586506, 0, -0.11093749999999993)),new PCurveEdge(new SemiEllipseCurve(V(0.11093750000000002, 0, 0), V(0, 0, 0.9938273849586506), V(0, 0.9938273849586506, 0), 0, 3.141592653589793), V(0.11093749999999993, 0, 0.9938273849586506), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), 0, 0.010062279327831384, null, V(0, 0.9938273849586506, 0), V(0, 0.993777073137507, -0.010000000000000007)),new PCurveEdge(new SemiEllipseCurve(V(0, 0.010000000000000009, 0), V(0, 0, -0.9999499987499375), V(0.9999499987499375, 0, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), 3.0304207486077566, 3.0936030871101416, null, V(-0.9937770731375071, 0, 0.11093749999999983), V(-0.9987987780446257, 0, 0.0479687500000002))]
			linkB3(assert, {edges: loop})
			assert.ok(s1.edgeLoopCCW(loop))
		},
	})
}