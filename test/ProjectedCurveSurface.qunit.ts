QUnit.module('PCS')
const pcs = new ProjectedCurveSurface(BezierCurve.EX2D, V3.Z, undefined, undefined, -1, 1)
QUnit.test('testSurface', function (assert) {
    testParametricSurface(assert, pcs)
})
QUnit.test('isTsWithSurface', function (assert) {
    const pcs = new ProjectedCurveSurface(BezierCurve.EX2D, V3.Z, undefined, undefined, -2, 2).scale(0.2,0.2,1).rotateX(-90*DEG)
    const ses = SemiEllipsoidSurface.UNIT
    const pic = ses.isCurvesWithSurface(pcs)[0]
    assert.ok(true, `<html><a href='brep2.html?edges=[${Edge.forCurveAndTs(pic).sce}]'>view</a>`)
    testISTs(assert, pic, new PlaneSurface(P3.XY), 1)
})