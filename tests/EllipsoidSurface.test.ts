import {
  outputLink,
  suite,
  suiteSurface,
  surfaceVolumeAndAreaTests,
  test,
  testCurve,
  testImplicitSurface,
  testISCurves,
  testLoopContainsPoint,
} from "./manager"

import { assertf, DEG, M4, V, V3, lt } from "ts3dutils"
import {
  B2T,
  BezierCurve,
  Edge,
  Face,
  P3,
  PCurveEdge,
  PICurve,
  PlaneSurface,
  PointVsFace,
  ProjectedCurveSurface,
  EllipseCurve,
  EllipsoidSurface,
  StraightEdge,
} from ".."

import { PI, sin, sqrt } from "../src/math"
suite("EllipsoidSurface", () => {
  const ses2 = EllipsoidSurface.UNIT.scale(2)
  const ses3 = EllipsoidSurface.UNIT.scale(2).translate(1, 2, 3).shearX(2, 3)
  suite("UNIT", () => suiteSurface(EllipsoidSurface.UNIT))
  suite("UNIT.scale(2)", () => suiteSurface(ses2))
  suite("UNIT.shearX(2, 3)", () =>
    suiteSurface(EllipsoidSurface.UNIT.shearX(2, 3)),
  )
  suite("UNIT.foo()", () => suiteSurface(EllipsoidSurface.UNIT.foo()))
  suite("UNIT.foo().flipped()", () =>
    suiteSurface(EllipsoidSurface.UNIT.foo().flipped()),
  )
  test("testSurface", (assert) => {
    testISCurves(
      assert,
      EllipsoidSurface.UNIT,
      new PlaneSurface(
        new P3(V(-1.249000902703301e-16, 1, 0), 0.11006944444444443),
      ),
      2,
    )
  })
  test("getAABB", (assert) => {
    assert.v3ArraysLike(ses2.getExtremePoints(), [V(0, 2, 0)])
    assert.v3ArraysLike(ses2.rotateZ(30 * DEG).getExtremePoints(), [
      V(-2, 0, 0),
      V(0, 2, 0),
    ])
  })
  test("loopContainsPoint", (assert) => {
    const s = new EllipsoidSurface(V3.O, V(5, 0, 0), V(0, 5, 0), V(0, 0, 5))
    const loop = [
      new PCurveEdge(
        new EllipseCurve(V(0, 0, 0), V(0, 0, -5), V(5, 0, 0)),
        V(0, 0, 5),
        V(0, 0, -5),
        PI,
        0,
        null,
        V(5, 0, 0),
        V(-5, 0, 0),
      ),
      new PCurveEdge(
        new EllipseCurve(V(0, 0, 0), V(0, 0, -5), V(-5, 0, 0)),
        V(0, 0, -5),
        V(0, 0, 5),
        0,
        PI,
        null,
        V(-5, 0, 0),
        V(5, 0, 0),
      ),
    ]
    const p = V(5, 0, 0)
    testLoopContainsPoint(assert, s, loop, p, PointVsFace.ON_EDGE)
  })
  test("loopContainsPoint 3", (assert) => {
    const s = new EllipsoidSurface(
      V3.O,
      V(-5, 6.123233995736766e-16, 0),
      V(-6.123233995736766e-16, -5, 0),
      V(0, 0, 5),
    )
    const loop = [
      new PCurveEdge(
        new EllipseCurve(
          V(0, 0, 0),
          V(0, 0, -5),
          V(5, 0, 0),
          0,
          3.141592653589793,
        ),
        V(0, 0, -5),
        V(1.1291713066130296, 0, -4.87082869338697),
        0,
        0.22779933669791175,
        null,
        V(5, 0, 0),
        V(4.87082869338697, 0, 1.1291713066130296),
      ),
      new PCurveEdge(
        new EllipseCurve(
          V(1.411764705882352, -2.1176470588235303, -1.4117647058823533),
          V(2.0917534572581817, 2.7890046096775745, -2.091753457258183),
          V(-2.874840149008801, 3.5206637865439285e-16, -2.874840149008799),
          0.7085839061491113,
          3.141592653589793,
        ),
        V(1.1291713066130291, 0, -4.87082869338697),
        V(0, -0.6994424542963525, -4.950836318555471),
        0.7085839061491117,
        1.0373562345961393,
        null,
        V(-3.544048444786543, -1.8149704259460577, -0.821592805867452),
        V(-3.262983117260863, -2.401508361946856, 0.3392794256594245),
      ),
      new PCurveEdge(
        new EllipseCurve(
          V(-1.4117647058823535, -2.117647058823529, -1.4117647058823533),
          V(2.0917534572581813, -2.789004609677577, 2.0917534572581813),
          V(2.8748401490088, -3.520663786543927e-16, -2.8748401490088007),
          0,
          2.43300874744068,
        ),
        V(0, -0.6994424542963525, -4.950836318555471),
        V(-1.1291713066130296, 1.382836026332681e-16, -4.87082869338697),
        2.1042364189936533,
        2.4330087474406805,
        null,
        V(-3.2629831172608608, 2.401508361946859, -0.33927942565942426),
        V(-3.5440484447865406, 1.8149704259460617, 0.8215928058674513),
      ),
      new PCurveEdge(
        new EllipseCurve(
          V(0, 0, 0),
          V(0, 0, 5),
          V(-5, 0, 0),
          0,
          3.141592653589793,
        ),
        V(-1.1291713066130296, 0, -4.87082869338697),
        V(0, 0, -5),
        2.9137933168918817,
        3.141592653589793,
        null,
        V(4.870828693386971, 0, -1.1291713066130287),
        V(5, 0, -6.123233995736766e-16),
      ),
    ]
    const p = V(-4.999999999999999, 0, 0)
    testLoopContainsPoint(assert, s, loop, p, PointVsFace.OUTSIDE)
  })
  test("intersect SES", (assert) => {
    const a = EllipsoidSurface.sphere(5)
    const b = EllipsoidSurface.sphere(1)
      .rotateAB(V3.Y, V3.X.negated())
      .translate(5, 2)
    testISCurves(assert, a, b, 2)
    testISCurves(assert, a, b.flipped(), 2)
    testISCurves(assert, a.flipped(), b, 2)
    testISCurves(assert, a.flipped(), b.flipped(), 2)

    testISCurves(assert, b, a, 2)
    testISCurves(
      assert,
      EllipsoidSurface.sphere(1),
      EllipsoidSurface.sphere(2).translate(-1.2),
      1,
    )
  })
  test("isCurvesWithProjectedCurveSurface is curves both y > 0", (assert) => {
    const s1 = EllipsoidSurface.UNIT
    const s2 = new ProjectedCurveSurface(
      BezierCurve.EX2D,
      V3.Z,
      undefined,
      undefined,
      -2,
      2,
    )
    testISCurves(
      assert,
      s1,
      s2
        .translate(0.2)
        .rotateZ(90 * DEG)
        .rotateX(10 * DEG),
      2,
    )
  })
  test("isCurvesWithProjectedCurveSurface is curves both cross y = 0", (assert) => {
    const s1 = EllipsoidSurface.UNIT
    const s2 = new ProjectedCurveSurface(
      BezierCurve.EX2D,
      V3.Z,
      undefined,
      undefined,
      -2,
      2,
    )
    testISCurves(assert, s1, s2.translate(0.2), 2)
  })
  test("isCurvesWithProjectedCurveSurface is curves both y < 0", (assert) => {
    const s1 = EllipsoidSurface.UNIT
    const s2 = new ProjectedCurveSurface(
      BezierCurve.EX2D,
      V3.Z,
      undefined,
      undefined,
      -2,
      2,
    )
    testISCurves(assert, s1, s2.translate(0.2).rotateZ(-90 * DEG), 0)
  })
  test("isCurvesWithProjectedCurveSurface one isCurve cross y = 0 twice, one isCurve y > 0", (assert) => {
    const s1 = EllipsoidSurface.UNIT
    const s2 = new ProjectedCurveSurface(
      BezierCurve.EX2D,
      V3.Z,
      undefined,
      undefined,
      -2,
      2,
    )
    testISCurves(
      assert,
      s1,
      s2
        .translate(0.2)
        .rotateZ(90 * DEG)
        .rotateX(80 * DEG),
      2,
    )
  })
  test("isCurvesWithProjectedCurveSurface", (assert) => {
    const s1 = EllipsoidSurface.UNIT.rotateZ(PI).flipped()
    const s2 = new ProjectedCurveSurface(
      BezierCurve.QUARTER_CIRCLE,
      V3.Z,
      undefined,
      undefined,
      -2,
      2,
    )
      .scale(0.2, 0.2, 2)
      .translate(0.1, -0.1, 1.2)
    testISCurves(assert, s1, s2, 2)
  })
  test("isCurvesWithProjectedCurveSurface 2", (assert) => {
    const s1 = EllipsoidSurface.UNIT
    const s2 = new ProjectedCurveSurface(
      new BezierCurve(
        V(0.30000000000000004, -0.1, 1.2),
        V(0.30000000000000004, 0.010456949966158674, 1.2),
        V(0.2104569499661587, 0.1, 1.2),
        V(0.10000000000000002, 0.1, 1.2),
        0,
        1,
      ),
      V(0, 0, 2),
      0,
      1,
      -1,
      0,
    )
    testISCurves(assert, s1, s2, 1)
  })
  test("loopCCW", (assert) => {
    const s1 = EllipsoidSurface.UNIT
    const loop = [
      new PCurveEdge(
        PICurve.forParametricStartEnd(
          new ProjectedCurveSurface(
            new BezierCurve(
              V(-0.1040625, -0.095, 1.2),
              V(-0.1040625, -0.030937500000000007, 1.2),
              V(-0.065, 0.010000000000000009, 1.2),
              V(0.047968750000000004, 0.010000000000000009, 1.2),
              0,
              1,
            ),
            V(0, 0, -1),
            0,
            1,
            -Infinity,
            Infinity,
          ),
          new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z),
          V(0.7084172378801449, 0.2005401547874497, 0),
          V(1, 0.20120122195537427, 0),
          0.01,
        ),
        V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257),
        V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503),
        29,
        0,
        null,
        V(-0.003388613668575339, 0, 0.00016274305244215148),
        V(
          -0.002184873829710333,
          -0.0006707444983087241,
          -0.00007184167849434948,
        ),
      ),
      new PCurveEdge(
        new EllipseCurve(
          V(0, 0, 0),
          V(0, 0, -1),
          V(-1, 1.2246467991473532e-16, 0),
          0,
          3.141592653589793,
        ),
        V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503),
        V(0, 0, 1),
        3.108723110778215,
        3.141592653589793,
        null,
        V(0.9994598452125503, -1.2239853003158589e-16, 0.032863624384797126),
        V(1, -1.2246467991473532e-16, 0),
      ),
      new PCurveEdge(
        new EllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)),
        V(0, 0, 1),
        V(0.11093749999999993, 0, 0.9938273849586506),
        3.141592653589793,
        3.030426330354509,
        null,
        V(1, 0, 0),
        V(0.9938273849586506, 0, -0.11093749999999993),
      ),
      new PCurveEdge(
        new EllipseCurve(
          V(0.11093750000000002, 0, 0),
          V(0, 0, 0.9938273849586506),
          V(0, 0.9938273849586506, 0),
          0,
          3.141592653589793,
        ),
        V(0.11093749999999993, 0, 0.9938273849586506),
        V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071),
        0,
        0.010062279327831384,
        null,
        V(0, 0.9938273849586506, 0),
        V(0, 0.993777073137507, -0.010000000000000007),
      ),
      new PCurveEdge(
        new EllipseCurve(
          V(0, 0.010000000000000009, 0),
          V(0, 0, -0.9999499987499375),
          V(0.9999499987499375, 0, 0),
          0,
          3.141592653589793,
        ),
        V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071),
        V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257),
        3.0304207486077566,
        3.0936030871101416,
        null,
        V(-0.9937770731375071, 0, 0.11093749999999983),
        V(-0.9987987780446257, 0, 0.0479687500000002),
      ),
    ]
    assertf(() => Edge.isLoop(loop))
    outputLink(assert, { edges: loop })
    assert.ok(s1.edgeLoopCCW(loop))
  })

  suite("triangular face from P3.XY to V3.Z", () => {
    const surface = EllipsoidSurface.UNIT
    const loop = [
      Edge.forCurveAndTs(EllipseCurve.UNIT, 10 * DEG, 40 * DEG),
      Edge.forCurveAndTs(
        new EllipseCurve(V3.O, V3.sphere(40 * DEG, 0), V3.Z),
        0,
        PI / 2,
      ),
      Edge.forCurveAndTs(
        new EllipseCurve(V3.O, V3.sphere(10 * DEG, 0), V3.Z),
        PI / 2,
        0,
      ),
    ]
    const face = Face.create(surface, loop)

    surfaceVolumeAndAreaTests(face, undefined, ((4 / 3) * PI * (30 / 360)) / 2)
    surfaceVolumeAndAreaTests(
      face.scale(1, 1, 2),
      ".scale(1, 1, 2)",
      (((4 / 3) * PI * (30 / 360)) / 2) * 2,
    )
    surfaceVolumeAndAreaTests(face.shearX(2, 2), ".shearX(2, 2)")
    surfaceVolumeAndAreaTests(face.foo(), ".foo()")
  })
  //suite('triangular face from P3.XY to V3.Z', () => {
  //
  //    const surface = EllipsoidSurface.UNIT
  //    const loop = [
  //        StraightEdge.throughPoints(V(1, 0, 1), V(1, 0, 0)),
  //        Edge.forCurveAndTs(EllipseCurve.forAB(1, 1), 0, PI / 2),
  //        Edge.forCurveAndTs(new EllipseCurve(V3.O, V(1, 0, 1), V(0, 1, 0)), PI / 2, 0),
  //    ]
  //    const face = Face.create(surface, loop)
  //
  //    surfaceVolumeAndAreaTests(face)
  //    surfaceVolumeAndAreaTests(face.transform(M4.FOO), '.transform(M4.foo)')
  //})

  test("mainAxes", (assert) => {
    const es = new EllipsoidSurface(V3.O, V(5, 0, -1), V(5, 1, 1), V(5, -1, 1))
    outputLink(assert, {
      mesh: `${es.sce}.toMesh()`,
      drPs: [V(-5, 1, -1)],
      edges: [
        StraightEdge.throughPoints(
          V(-5, 1, -1),
          V(-5, 1, -1).plus(
            V(
              -9.660064978873681e-7,
              0.999999999962889,
              0.000008560880502697739,
            ).times(100),
          ),
        ),
        StraightEdge.throughPoints(
          V(-5, 1, -1),
          V(-5, 1, -1).plus(
            V(
              0.14427746420619014,
              -0.6925318281897137,
              -1.413919149220665,
            ).times(100),
          ),
        ),
      ],
    })
    const esn = es.mainAxes()
    assert.ok(esn.f1.isPerpendicularTo(esn.f2))
    assert.ok(esn.f2.isPerpendicularTo(esn.f3))
    assert.ok(esn.f3.isPerpendicularTo(esn.f1))
    assert.ok(es.matrixInverse.times(esn.matrix).isOrthogonal())
  })

  test("loopContainsPoint", (assert) => {
    const testFace = B2T.rotateEdges(
      [
        Edge.forCurveAndTs(EllipseCurve.UNIT, 0, 90 * DEG).rotateX(90 * DEG),
        StraightEdge.throughPoints(V3.Z, V3.X),
      ],
      45 * DEG,
      "blah",
    ).faces.find((face) => face.surface instanceof EllipsoidSurface)

    const p1 = V3.sphere(10 * DEG, 10 * DEG)
    testLoopContainsPoint(
      assert,
      testFace.surface,
      testFace.contour,
      p1,
      PointVsFace.INSIDE,
    )
    const p2 = V3.sphere(10 * DEG, -10 * DEG)
    testLoopContainsPoint(
      assert,
      testFace.surface,
      testFace.contour,
      p2,
      PointVsFace.OUTSIDE,
    )
    testLoopContainsPoint(
      assert,
      testFace.surface.foo(),
      testFace.contour.map((edge) => edge.foo()),
      M4.FOO.transformPoint(p1),
      PointVsFace.INSIDE,
    )
  })
  test("loopContainsPoint 2", (assert) => {
    const testFace = B2T.sphere(1).faces[0]

    const p1 = new V3(0, 1, 0)
    testLoopContainsPoint(
      assert,
      testFace.surface,
      testFace.contour,
      p1,
      PointVsFace.INSIDE,
    )
  })
  test("isCurvesWithPlane", (assert) => {
    const es = EllipsoidSurface.sphere(5)
    testISCurves(
      assert,
      es,
      new PlaneSurface(new P3(V(0, -1, 0.1).unit(), 4)),
      0,
    )
    testISCurves(
      assert,
      es,
      new PlaneSurface(new P3(V(0, -1, 0.1).unit(), 4)).flipped(),
      0,
    )
    testISCurves(
      assert,
      es,
      new PlaneSurface(new P3(V3.sphere(-PI / 2, sin(3 / 5)), 4)),
      0,
    )
    testISCurves(
      assert,
      es,
      new PlaneSurface(new P3(V(0, -0.1, 1).unit(), 4)),
      1,
    )
    testISCurves(assert, es, new PlaneSurface(new P3(V(0, 0, 1), 4)), 1)
    testISCurves(
      assert,
      es,
      new PlaneSurface(new P3(V(0, 0.1, 1).unit(), 4)),
      2,
    )
    testISCurves(
      assert,
      es,
      new PlaneSurface(new P3(V(0, 1, 0.1).unit(), 4)),
      2,
    )
    testISCurves(assert, es, new PlaneSurface(P3.XY), 1)
    testISCurves(assert, es, new PlaneSurface(P3.XY.flipped()), 1)

    // slices perpendicular to V3.Y
    testISCurves(assert, es, new PlaneSurface(P3.ZX), 2)
    testISCurves(assert, es, new PlaneSurface(P3.ZX.flipped()), 2)
    testISCurves(assert, es, new PlaneSurface(P3.ZX).translate(0, 2, 0), 2)
    testISCurves(assert, es, new PlaneSurface(P3.ZX).translate(0, -2, 0), 0)

    // slices perpendicular to V3.X
    testISCurves(assert, es, new PlaneSurface(P3.YZ), 1)
  })

  test("transform4", (assert) => {
    const s = new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z, -PI, PI).flipped()
    // prettier-ignore
    const m4 = new M4(
            1,0,0,0,
            0,1,0,0,
            0,0,1,0,
            1,2,3,4,
        ).translate(2, 3)
    const m5 = EllipsoidSurface.unitTransform4(m4)
    const s2 = s.transform(m5)
    outputLink(assert, {
      mesh:
        "[" +
        s.toSource() +
        ".toMesh(), " +
        s2.toSource() +
        ".toMesh(), " +
        s.toSource() +
        ".toMesh().transform(" +
        m4.toSource() +
        ")]",
      drPs: [m5.getTranslation()],
      boxes: [m5],
    })
  })
})
