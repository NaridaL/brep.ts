import {
  Assert,
  outputLink,
  suite,
  test,
  testCurve,
  testCurveISInfos,
  testISTs,
} from "./manager"

import {
  arrayFromFunction,
  arraySamples,
  checkDerivate,
  DEG,
  getLast,
  lerp,
  M4,
  V,
  V3,
  Vector,
  VV,
} from "ts3dutils"
import {
  BezierCurve,
  Edge,
  edgeForCurveAndTs,
  EllipsoidSurface,
  HyperbolaCurve,
  L3,
  NURBS,
  P3,
  XiEtaCurve,
} from ".."
import { PI, SQRT1_2, SQRT2 } from "../src/math"

suite("NURBS", () => {
  test("openUniformKnots", (assert) => {
    assert.deepEqual(NURBS.openUniformKnots(7, 3, 0, 4), [
      0,
      0,
      0,
      0,
      1,
      2,
      3,
      4,
      4,
      4,
      4,
    ])
  })

  test("linear", (assert) => {
    const nurbs = new NURBS([VV(0, 0, 0, 1), VV(1, 1, 0, 1)], 1, [0, 0, 1, 1])
    testCurve(assert, nurbs)
  })

  test("quadratic", (assert) => {
    const nurbs = new NURBS(
      [
        VV(-3, 0, 0, 1),
        VV(-2, 0, 0, 1),
        VV(-1, 0, 0, 1),
        VV(0, 6, 0, 1),
        VV(1, 0, 0, 1),
        VV(2, 0, 0, 1),
        VV(3, 0, 0, 1),
      ],
      3,
      [-2, -2, -2, -2, -1, 0, 1, 2, 2, 2, 2],
    )
    assert.equal(nurbs.tMin, -2)
    assert.equal(nurbs.tMax, 2)
    testCurve(assert, nurbs)
  })

  test("circle", (assert) => {
    const nurbs = new NURBS(
      [
        VV(1, 0, 0, 1),
        VV(1, 1, 0, 1).times(SQRT1_2),
        VV(0, 1, 0, 1),
        VV(-1, 1, 0, 1).times(SQRT1_2),
        VV(-1, 0, 0, 1),
        VV(-1, -1, 0, 1).times(SQRT1_2),
        VV(0, -1, 0, 1),
      ],
      2,
      [0, 0, 0, PI / 2, PI / 2, PI, PI, (3 * PI) / 2, (3 * PI) / 2, 2 * PI],
    )
    testCurve(assert, nurbs)
  })

  test("quartic", (assert) => {
    const nurbs = new NURBS(
      [
        VV(-4, 0, 0, 1),
        VV(-3, 0, 0, 1),
        VV(-2, 0, 0, 1),
        VV(-1, 0, 0, 1),
        VV(0, 1, 0, 1),
        VV(1, 0, 0, 1),
        VV(2, 0, 0, 1),
        VV(3, 0, 0, 1),
        VV(4, 0, 0, 1),
      ],
      4,
      [0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5, 5],
    )
    testCurve(assert, nurbs)
  })

  test("EX2D", (assert) => testCurve(assert, NURBS.EX2D))
  test("EX2D transformed", (assert) =>
    testCurve(
      assert,
      NURBS.EX2D.transform(
        M4.forSys(
          V(0, 0.05, 0),
          V(-0.05, 0, 0),
          V(0, 0, 0.05),
          V(1.3860347309027767, -5.6602156597222235, 0),
        ),
      ),
    ))
  test("EX3D", (assert) => testCurve(assert, NURBS.EX3D))

  test("fromBezier", (assert) => {
    const nurbs = NURBS.fromBezier(BezierCurve.EX3D)
    const xx = arrayFromFunction(100, (i) => {
      const t = lerp(nurbs.tMin, nurbs.tMax, i / 99)
      return [t, nurbs.at(t), nurbs.tangentAt(t)] as [number, V3, V3]
    })
    outputLink(assert, {
      drPs: xx.map((x) => x[1]),
      drVs: xx
        .filter((x) => isFinite(x[2].x) && !x[2].likeO())
        .map((x) => ({ anchor: x[1], v: x[2] })),
      edges: [edgeForCurveAndTs(BezierCurve.EX3D), edgeForCurveAndTs(nurbs)],
    })
  })

  test("bezier", (assert) => {
    testCurve(assert, NURBS.fromBezier(BezierCurve.EX2D))
  })

  test("UnitCircle", (assert) => {
    console.log(NURBS.UnitCircle().str)
    testCurve(
      assert,
      NURBS.UnitCircle(),
      undefined,
      "half circle in 2 segments",
    )
    testCurve(
      assert,
      NURBS.UnitCircle(5, -PI, PI),
      undefined,
      "full circle in 5 segments",
    )
    testCurve(
      assert,
      NURBS.UnitCircle(3, -PI, PI),
      undefined,
      "full circle in 3 segments",
    )
  })

  test("getSegments 1", (assert) => {
    const nurbs = NURBS.fromV3s(
      [
        V(178, 256),
        V(29, 245),
        V(12, 159),
        V(24, 50),
        V(121, 15),
        V(247, 11),
        V(249, 252),
      ],
      4,
    )
    const segments = nurbs.getSegments()
    outputLink(assert, {
      edges: [
        edgeForCurveAndTs(nurbs).translate(V3.Z),
        ...segments.map((s) => edgeForCurveAndTs(s)),
      ],
    })
    assert.equal(
      segments.length,
      nurbs.points.length - nurbs.degree,
      "number of segments",
    )
    assert.v3like(segments[0].points[0].p3(), nurbs.points[0].p3())
    assert.v3like(
      getLast(getLast(segments).points).p3(),
      getLast(nurbs.points).p3(),
    )
    segments.forEach((s) => {
      assert.ok(s instanceof NURBS)
      assert.equal(s.degree, nurbs.degree)
      assert.ok(s.isBezier(), "s.isBezier()")
    })
  })

  test("getSegments on NURBS with duplicate knots", (assert) => {
    const nurbs = NURBS.UnitCircle()
    const segments = nurbs.getSegments()
    outputLink(assert, {
      edges: [
        edgeForCurveAndTs(nurbs).translate(V3.Z),
        ...segments.map((s) => edgeForCurveAndTs(s)),
      ],
    })
    assert.equal(
      segments.length,
      nurbs.points.length - nurbs.degree - 1,
      "number of segments",
    )
    assert.v3like(segments[0].points[0].p3(), nurbs.points[0].p3())
    assert.v3like(
      getLast(getLast(segments).points).p3(),
      getLast(nurbs.points).p3(),
    )
    segments.forEach((s) => {
      assert.ok(s instanceof NURBS)
      assert.equal(s.degree, nurbs.degree)
      assert.ok(s.isBezier(), "s.isBezier()")
    })
  })

  test("getSegments degree=4, multi=3", (assert) => {
    const nurbs = NURBS.EX2D.withKnot(NURBS.EX2D.knots[5]).withKnot(
      NURBS.EX2D.knots[5],
    )
    const segments = nurbs.getSegments()
    outputLink(assert, {
      //
      edges: [
        edgeForCurveAndTs(nurbs).translate(V3.Z),
        ...segments.map((s) => edgeForCurveAndTs(s)),
      ],
    })
    assert.equal(
      segments.length,
      nurbs.points.length - nurbs.degree - 2,
      "number of segments",
    )
    assert.v3like(segments[0].points[0].p3(), nurbs.points[0].p3())
    assert.v3like(
      getLast(getLast(segments).points).p3(),
      getLast(nurbs.points).p3(),
    )
    segments.forEach((s) => {
      assert.ok(s instanceof NURBS)
      assert.equal(s.degree, nurbs.degree)
      assert.ok(s.isBezier(), "s.isBezier()")
    })
  })

  test("split on new knot", (assert) => {
    const nurbs = NURBS.EX2D
    const t0 = 0.5
    assert.notOk(nurbs.knots.includes(t0))
    const [n0, n1] = nurbs.split(t0)
    console.log(n0.points.length, n1.points.length)
    outputLink(assert, {
      //drPs: arrayFromFunction(12, (i, l) => n0.at(lerp(n0.tMin, n0.tMax, i/(l - 1)))),
      edges: [
        edgeForCurveAndTs(nurbs).translate(V3.Z),
        ...[n0, n1].map((s) => edgeForCurveAndTs(s)),
      ],
    })
    assert.equal(n0.tMax, t0)
    assert.equal(n0.tMin, nurbs.tMin)
    assert.v3like(n0.at(t0), nurbs.at(t0))
    assert.equal(n1.tMin, t0)
    assert.equal(n1.tMax, nurbs.tMax)
    assert.v3like(n1.at(t0), nurbs.at(t0))
  })

  test("split on existing knot with multi=1", (assert) => {
    const nurbs = NURBS.EX2D
    const t0 = nurbs.knots[7]
    const [n0, n1] = nurbs.split(t0)
    console.log(n0.points.length, n1.points.length)
    outputLink(assert, {
      drPs: arrayFromFunction(12, (i, l) =>
        n0.at(lerp(n0.tMin, n0.tMax, i / (l - 1))),
      ),
      edges: [
        edgeForCurveAndTs(nurbs).translate(V3.Z),
        ...[n0, n1].map((s) => edgeForCurveAndTs(s)),
      ],
    })
    assert.equal(n0.tMax, t0)
    assert.equal(n0.tMin, nurbs.tMin)
    assert.v3like(n0.at(t0), nurbs.at(t0))
    assert.equal(n1.tMin, t0)
    assert.equal(n1.tMax, nurbs.tMax)
    assert.v3like(n1.at(t0), nurbs.at(t0))
  })

  test("elevateDegreeBezier", (assert) => {
    const nurbsBezier = NURBS.EX2D.getSegments()[0]
    const nurbsBezier2 = nurbsBezier.elevateDegreeBezier()
    const curves = [nurbsBezier, nurbsBezier2.translate(0, 0, 10)]
    outputLink(assert, {
      edges: curves.map((n) => edgeForCurveAndTs(n)),
    })
  })

  test("elevateDegree", (assert) => {
    const nurbs = NURBS.EX2D
    const nurbs2 = nurbs.elevateDegree()
    const curves = [nurbs, nurbs2.translate(0, 0, 10)]
    const edges = curves.map((n) => edgeForCurveAndTs(n))
    outputLink(assert, {
      drPs: edges.flatMap((e) => e.getVerticesNo0()),
      edges: edges,
    })
    assert.equal(nurbs2.degree, nurbs.degree + 1)
    assert.v3like(nurbs2.points[0].p3(), nurbs.points[0].p3())
    assert.v3like(getLast(nurbs2.points).p3(), getLast(nurbs.points).p3())
  })

  test("elevateDegree 2", (assert) => {
    // test taken from https://pages.mtu.edu/~shene/COURSES/cs3621/LAB/curve/elevation.html
    // I got a better result (less control points)!
    const nurbs = NURBS.fromV3s(
      [
        V(50, 229),
        V(13, 186),
        V(31, 64),
        V(72, 189),
        V(133, 191),
        V(91, 42),
        V(221, 56),
        V(203, 227),
      ],
      5,
    )
    const nurbs2 = nurbs.elevateDegree()
    const curves = [nurbs, nurbs2.translate(0, 0, 10)]
    const edges = curves.map((n) => edgeForCurveAndTs(n))
    outputLink(assert, {
      edges: edges,
    })
    assert.equal(nurbs.points.length, 8)
    assert.equal(nurbs2.points.length, 11)
    assert.equal(nurbs2.degree, nurbs.degree + 1)
    assert.v3like(nurbs2.points[0].p3(), nurbs.points[0].p3())
    assert.v3like(getLast(nurbs2.points).p3(), getLast(nurbs.points).p3())
  })

  test("removeKnots multi=1", (assert) => {
    const nurbs = NURBS.EX2D.getSegments()[0].withKnot(0.1)
    const nurbs2 = nurbs.removeKnot(0.1)!
    assert.equal(nurbs2.degree, nurbs.degree)
    assert.equal(nurbs2.knots.length, nurbs.knots.length - 1)
    outputLink(assert, {
      edges: [nurbs, nurbs2.translate(0, 0, 10)].map((n) =>
        edgeForCurveAndTs(n),
      ),
    })
  })

  test("removeKnots 2 multi>1", (assert) => {
    const t = NURBS.EX2D.knots[NURBS.EX2D.degree + 2]
    const target = NURBS.EX2D.withKnot(t)
    const nurbs = NURBS.EX2D.withKnot(t).withKnot(t)
    const nurbs2 = nurbs.removeKnot(t)!
    assert.equal(nurbs2.degree, nurbs.degree)
    assert.equal(nurbs2.knots.length, nurbs.knots.length - 1)
    outputLink(assert, {
      edges: [nurbs, nurbs2.translate(0, 0, 10)].map((n) =>
        edgeForCurveAndTs(n),
      ),
    })
  })

  test("isTsWithPlane 1", (assert) =>
    testISTs(assert, NURBS.EX2D, new P3(V3.XY.unit(), 120), 4))
  test("isTsWithPlane 2", (assert) =>
    testISTs(
      assert,
      NURBS.EX2D,
      P3.normalOnAnchor(
        NURBS.EX2D.tangentAt(0.65).cross(V3.Z),
        NURBS.EX2D.at(0.65),
      ),
      2,
    ))
  test("isTsWithPlane nurbs contained by plane", (assert) =>
    testISTs(assert, NURBS.EX2D, P3.XY, 0))
  test("isTsWithSurface EllipsoidSurface which does not intersect control poly", (assert) => {
    const surface = EllipsoidSurface.UNIT.rotateZ(-90 * DEG)
      .scale(20)
      .translate(NURBS.EX2D.at(0.65))
    testISTs(assert, NURBS.EX2D, surface, 2)
  })
  test("isTsWithPlane EX3D P3.XY", (assert) =>
    testISTs(assert, NURBS.EX3D, P3.XY, 3))
  test("isTsWithPlane EX3D P3.YZ", (assert) =>
    testISTs(assert, NURBS.EX3D, P3.YZ, 4))

  suite("isInfosWithLine", () => {
    test("line in curve plane", (assert) =>
      testCurveISInfos(
        assert,
        NURBS.EX2D,
        new L3(NURBS.EX2D.at(0.65), NURBS.EX2D.tangentAt(0.65).unit()),
        2,
      ))
    test("parallel to curve plane", (assert) =>
      testCurveISInfos(
        assert,
        NURBS.EX2D,
        new L3(
          NURBS.EX2D.at(0.65).plus(V3.Z),
          NURBS.EX2D.tangentAt(0.65).unit(),
        ),
        0,
      ))
    test("crosses curve plane, is", (assert) =>
      testCurveISInfos(
        assert,
        NURBS.EX2D,
        new L3(
          NURBS.EX2D.at(0.65),
          NURBS.EX2D.tangentAt(0.65).plus(V3.Z).unit(),
        ),
        1,
      ))
    test("crosses curve plane, no is", (assert) =>
      testCurveISInfos(
        assert,
        NURBS.EX2D,
        new L3(
          NURBS.EX2D.at(0.65).plus(V3.Z),
          NURBS.EX2D.tangentAt(0.65).plus(V3.Z).unit(),
        ),
        0,
      ))
    test("EX3D, no is", (assert) =>
      testCurveISInfos(assert, NURBS.EX3D, new L3(V3.O, V(0, -1, 1).unit()), 0))
    test("EX3D, 1 is", (assert) =>
      testCurveISInfos(
        assert,
        NURBS.EX3D,
        new L3(NURBS.EX3D.at(5), V(0, -1, 1).unit()),
        1,
      ))
    test("EX3D, 2 is", (assert) =>
      testCurveISInfos(
        assert,
        NURBS.EX3D,
        L3.throughPoints(NURBS.EX3D.at(-5), NURBS.EX3D.at(10), -100, 200),
        2,
      ))
  })

  function testSimplify(assert: Assert, c: NURBS) {
    const simple = c.simplify()
    assert.ok(simple instanceof XiEtaCurve, "simple instanceof XiEtaCurve")
    if (c.at(c.tMin).like(simple.at(simple.tMin))) {
      assert.v3like(c.at(c.tMin), simple.at(simple.tMin))
      assert.v3like(c.at(c.tMax), simple.at(simple.tMax))
    } else {
      assert.v3like(c.at(c.tMin), simple.at(simple.tMax))
      assert.v3like(c.at(c.tMax), simple.at(simple.tMin))
    }
    arraySamples(0, 1, 8).forEach((t) => {
      assert.ok(simple.containsPoint(c.at(t)))
    })
    outputLink(assert, {
      edges: [edgeForCurveAndTs(simple), edgeForCurveAndTs(c.translate(V3.Z))],
      drPs: arraySamples(0, 1, 8).map((t) => c.at(t)),
    })
  }

  test("simplify ellipse", (assert) => {
    const c = NURBS.fromV3s(
      [V(1, 0, 0), V(1, 1, 0), V(0, 1, 0)],
      2,
      undefined,
      [1.1, SQRT2 / 3, 1],
    )
    testSimplify(assert, c)
  })

  test("simplify hyperbola", (assert) => {
    const c = NURBS.fromV3s(
      [V(1, 0, 0), V(1, 1, 0), V(0, 1, 0)],
      2,
      undefined,
      [1, 2, 1.3],
    )
    testSimplify(assert, c)
  })
  test("simplify hyperbola with w linear", (assert) => {
    const c = NURBS.fromV3s(
      [V(1, 0, 0), V(1, 1, 0), V(0, 1, 0)],
      2,
      undefined,
      [1, 2, 3],
    )
    testSimplify(assert, c)
  })
  test("simplify circle", (assert) => {
    const c = NURBS.UnitCircle(1, 0, PI / 2)
    testSimplify(assert, c)
  })
  test("simplify parabola", (assert) => {
    const c = NURBS.fromV3s([V3.X, V3.O, V3.Y], 2, undefined, [
      3,
      3,
      3,
    ]).translate(V3.Z)
    testSimplify(assert, c)
  })

  test("simplifyUnit2", (assert) => {
    //const [w0, w1, w2] = [1, SQRT2/3, 1]
    //const [w0, w1, w2] = [1, SQRT2/3, 1].map(x => 2*x)
    //const [w0, w1, w2] = [1/4, 2, 4]
    const [w0, w1, w2] = [1, 1, 1]
    const c1 = NURBS.fromV3s([V3.X, V3.O, V3.Y], 2, undefined, [
      w0,
      w1,
      w2,
    ]).translate(V3.Z)
    const c2 = NURBS.simplifyUnit2(w0, w1, w2)
    outputLink(assert, {
      edges: [c1, c2].map((c) => edgeForCurveAndTs(c)),
      drPs: arraySamples(0, 1, 10).map((t) => c1.at(t)),
    })
  })

  test("fromHyperbola", (assert) => {
    const hb = HyperbolaCurve.XY.withBounds(-1, 2)
    const nurbs = NURBS.fromHyperbola(hb)
    outputLink(assert, {
      edges: [hb, nurbs].map((c) => edgeForCurveAndTs(c)),
      drPs: arraySamples(0, 1, 33).map((t) => nurbs.at(t)),
    })
  })
  test("fromHyperbola FOO", (assert) => {
    const hb = HyperbolaCurve.XY.withBounds(-1, 2).transform(M4.FOO)
    const nurbs = NURBS.fromHyperbola(hb)
    outputLink(assert, {
      edges: [hb, nurbs].map((c) => edgeForCurveAndTs(c)),
      drPs: arraySamples(0, 1, 33).map((t) => nurbs.at(t)),
    })
  })
})
