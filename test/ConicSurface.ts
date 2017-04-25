/// <reference path="../src/NLA.ts" />
/// <reference path="../src/Vector3.ts" />
/// <reference path="../src/surface/ConicSurface.ts" />
/// <reference path="../src/curve/SemiEllipseCurve.ts" />
/// <reference path="../src/Plane3.ts" />

import 'mocha'
import {assert} from 'chai'


suite('ConicSurface', function () {
    test('ConicSurface.isCoplanarTo', function () {
        const unitCone = ConicSurface.UNIT
        assert.ok(unitCone.matrix.isIdentity(), 'UCS.matrix.isIdentity()')
        assert.V3like(unitCone.parametricFunction()(0, 3), V(3, 0, 3))
        const ellipseAtZ3 = SemiEllipseCurve.UNIT.scale(3, 3, 3).translate(0, 0, 3)
        const planeAtZ3 = P3.XY.translate(0, 0, 3)
        const issAtZ3 = unitCone.isCurvesWithPlane(planeAtZ3)
        assert.equal(issAtZ3.length, 1)
        assert.ok(ellipseAtZ3.isColinearTo(issAtZ3[0]))
        assert.ok(unitCone.containsEllipse(ellipseAtZ3))


        const scaledUnit = ConicSurface.UNIT.scale(2, 2, 1)
        assert.notOk(scaledUnit.isCoplanarTo(unitCone))
        assert.notOk(unitCone.isCoplanarTo(scaledUnit))
        const ell1 = unitCone.isCurvesWithPlane(new P3(V(2, 3, 10).unit(), 10))[0]
        assert.ok(unitCone.containsEllipse(ell1), 'UCS.containsEllipse(ell1)')
        const ell2 = unitCone.isCurvesWithPlane(new P3(V(1, 1, 2).unit(), 4))[0]
        const ell1Cone = ConicSurface.atApexThroughEllipse(V3.O, ell1, 1)
        const ell2Cone = ConicSurface.atApexThroughEllipse(V3.O, ell2, 1)
        console.log(ell1Cone)
        assert.ok(unitCone.isCoplanarTo(ell1Cone))
        assert.ok(unitCone.isCoplanarTo(ell2Cone))
        assert.ok(ell1Cone.isCoplanarTo(ell2Cone))
        assert.ok(ell2Cone.isCoplanarTo(ell1Cone))
        assert.ok(ell1Cone.foo().isCoplanarTo(ell2Cone.foo()))
    })
    test('ConicSurface.containsParabola', function () {
        const unitCone = ConicSurface.UNIT
        const pb = unitCone.isCurvesWithPlane(new P3(V(1, 0, 1).unit(), 4))[0]
        assert.ok(unitCone.containsParabola(pb))

        const c2 = unitCone.shearedX(2, 3)
        const pb2 = c2.isCurvesWithPlane(new P3(V(1, 0, 1).unit(), 4).shearedX(2, 3))[0]
        assert.ok(c2.containsParabola(pb2))
    })
})