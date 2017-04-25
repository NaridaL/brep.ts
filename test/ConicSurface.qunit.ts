{
	QUnit.module('ConicSurface')
	const UCS = ConicSurface.UNIT
	registerTests({
		'testSurface'(assert) {
			testParametricSurface(assert, UCS)
			testParametricSurface(assert, ConicSurface.UNIT.scale(2, 2, 1))
		},
		'isCoplanarTo'(assert) {
			assert.ok(UCS.matrix.isIdentity(), 'UCS.matrix.isIdentity()')
			assert.V3like(UCS.parametricFunction()(0, 3), V(3, 0, 3))
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
		},
		'isCurvesWithPlane'(assert) {
			testISCurves(assert, UCS, new P3(V(1, 1, 2).unit(), 4), 1)
			testISCurves(assert, UCS, P3.XY.translate(0, 0, 3), 1)
			testISCurves(assert, UCS, P3.XY.translate(0, 0, 3).flipped(), 1)
		},
		'containsParabola'(assert) {
			const pb = UCS.isCurvesWithPlane(new P3(V(1, 0, 1).unit(), 4))[0]
			assert.ok(UCS.containsParabola(pb))

			const c2 = UCS.shearedX(2, 3)
			const pb2 = c2.isCurvesWithPlane(new P3(V(1, 0, 1).unit(), 4).shearedX(2, 3))[0]
			assert.ok(c2.containsParabola(pb2))
		},
		'containsPoint'(assert) {
			const face = B2T.cone(1,1,PI).faces.find(face => face.surface instanceof ConicSurface)
			testLoopContainsPoint(assert, face.surface, face.contour, V(-0.2, 0, 0.2), PointVsFace.ON_EDGE)
		}
	})
}