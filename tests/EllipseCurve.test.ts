import {
  outputLink,
  testCurve,
  testCurveISInfos,
  testCurveTransform,
  testISTs,
} from "./manager"

import {
  arrayFromFunction,
  arraySamples,
  DEG,
  lerp,
  M4,
  V,
  V3,
  Vector,
} from "ts3dutils"
import {
  BezierCurve,
  EllipseCurve,
  EllipsoidSurface,
  P3,
  ParabolaCurve,
  PlaneSurface,
  parabola4Projection,
  edgeForCurveAndTs,
} from "../src"

const { PI } = Math

describe("EllipseCurve", () => {
  const curve = EllipseCurve.UNIT.shearX(2, 1)

  test("withBounds", () => {
    const newCurve = curve.withBounds(1, 2)
    expect(newCurve.tMin).toBe(1)
    expect(newCurve.tMax).toBe(2)
  })
  test("testCurve", () => {
    testCurve(EllipseCurve.UNIT)
    testCurve(curve)
  })
  test("UNIT.shearX(2, 3)", () => testCurve(EllipseCurve.UNIT.shearX(2, 2)))
  test("isTsWithPlane", () => {
    const plane = new P3(V(2, 7, 1).unit(), 2)
    testISTs(curve.scale(1, 3, 1), plane, 2)
  })
  //test('rightAngled', () => {
  //	const curveRA = curve.rightAngled()
  //	assert.ok(curveRA.f1.isPerpendicularTo(curveRA.f2))
  //	assert.ok(curveRA.isColinearTo(curve))
  //	arrayRange(-10, 10, 1).forEach(t => assert.ok(curveRA.containsPoint(curve.at(t))))
  //},
  test("isTsWithSurface(EllipsoidSurface)", () => {
    const s = EllipsoidSurface.sphere(5)
    const c = new EllipseCurve(
      V(5, 2),
      V3.Z.negated(),
      V(-1, 1.2246467991473532e-16, 0),
      0,
      PI,
    )
    testISTs(c, s, 2)
  })

  test("isTsWithSurface(PlaneSurface)", () => {
    const c = EllipseCurve.UNIT.translate(1.2, -1)
    const s = new PlaneSurface(P3.ZX)
    testISTs(c, s, 1)
  })

  test("isTsWithSurface(PlaneSurface) 2", () => {
    const c = EllipseCurve.UNIT
    const s = P3.YZ.translate(0.5, 0)
    testISTs(c, s, 1)
  })

  test("distanceToPoint", () => {
    const curve = EllipseCurve.forAB(10, 15)
    const p = V(12, 12)
    const closestT = curve.closestTToPoint(p)
    const pDist = curve.at(closestT).distanceTo(p)
    const EPS = 0.001
    expect(pDist).toBeLessThan(curve.at(closestT - EPS).distanceTo(p))
    expect(pDist).toBeLessThan(curve.at(closestT + EPS).distanceTo(p))
  })

  test("isColinearTo", () => {
    expect(
      EllipseCurve.forAB(1, 2).isColinearTo(EllipseCurve.forAB(1, -2)),
    ).toBeTruthy()
  })

  const c1 = EllipseCurve.semicircle(5)
  test("isInfosWithEllipse", () => {
    const c1 = EllipseCurve.semicircle(5)
    const c2 = EllipseCurve.semicircle(5, V(3, 0))
    testCurveISInfos(c1, c2, 1, "c1 c2")

    const verticalEllipse = new EllipseCurve(V(2, 0), V(1, 1), V(1, 10))
    testCurveISInfos(c1, verticalEllipse, 2, "c1 verticalEllipse")

    const verticalEllipse2 = new EllipseCurve(V(10, 2), V(1, 1), V(1, 10))
    testCurveISInfos(c1, verticalEllipse2, 0, "c1 verticalEllipse2")

    const smallEllipse = EllipseCurve.forAB(2, 3)
    testCurveISInfos(c1, smallEllipse, 0, "c1 smallEllipse")
  })
  test("c1 test", () => {
    const test = new EllipseCurve(V(6, 1, 0), V(3, 1, 0), V(4, 0, 0))
    testCurveISInfos(c1, test, 1, "c1 test")
  })
  test("isInfosWithBezier2D", () => {
    const ell = EllipseCurve.forAB(3, 1)
    const bez = BezierCurve.graphXY(2, -3, -3, 2, -2, 3)
    testCurveISInfos(ell, bez, 3)
  })
  test("transform4", () => {
    const c = EllipseCurve.UNIT.withBounds(0, 3)
    const ps = arrayFromFunction(128, (i, l) =>
      c.at(lerp(c.tMin, c.tMax, i / (l - 1))),
    )
    const p3 = new P3(V3.X, 2)
    const proj1 = M4.projectPlanePoint(V(-2, 0, -1), p3)
    testCurveTransform(c, proj1)
  })

  test("transform4 2", () => {
    const c = EllipseCurve.UNIT
    const ps = arrayFromFunction(128, (i, l) =>
      c.at(lerp(c.tMin, c.tMax, i / (l - 1))),
    )
    const mv = M4.scale(0.5, 1, 1)
      .rotateZ(20 * DEG)
      .translate(0, -2)
      .rotateX(90 * DEG)
    const perspective = M4.perspective(45, 1, 1, 2).times(mv)
    //const m4 = mv
    //const m4 = M4.product(M4.rotateX(90* DEG), perspective,M4.rotateX(-90* DEG))
    const m4 = perspective

    testCurveTransform(c, m4)
  })

  test("transform4 3", () => {
    const c = new EllipseCurve(V3.O, V3.X, V3.Y, -3, 2).translate(1, -4, 0)
    const m = M4.product(
      M4.rotateX(90 * DEG),
      M4.perspective(45, 1, 2, 5),
      M4.rotateX(-90 * DEG),
    )
    testCurveTransform(c, m)
  })

  test("transform4 4", () => {
    const c = new EllipseCurve(
      V(1, 0, -3),
      V3.X,
      V(0, 6.123233995736766e-17, 1),
      0,
      3.141592653589793,
    )
    const m = M4.perspective(45, 1, 2, 4)
    testCurveTransform(c, m)
  })

  test("transform4 map parabola onto itself", () => {
    //const c = new EllipseCurve(V3.O, V3.X, V3.Y, -3, 3)
    const c = ParabolaCurve.XY.withBounds(0.2, 2)
    const m = M4.product(
      M4.translate(0, -1),
      // prettier-ignore
      new M4(
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 1, 0, -1,
      ),
      M4.translate(0, 1),
    )
    console.log(m.toString())
    testCurveTransform(c, m)
  })

  test("transform4 map parabola onto line", () => {
    //const c = new EllipseCurve(V3.O, V3.X, V3.Y, -3, 3)
    const c = ParabolaCurve.XY.withBounds(0 - 2, 2)
    const m = new M4(0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0)
    const m2 = M4.translate(m.m[3], m.m[7], m.m[11]).transform(m)
    const P = M4.product(m2.as3x3().inversed())
    console.log(m.rank())
    console.log("P")
    console.log(P.toString())
    console.log("m2")
    console.log(m2.toString())
    console.log("P * m")
    console.log(P.times(m).toString())
    console.log()

    // prettier-ignore
    //const m = M4.product(M4.translate(0, -1),
    //   new M4(
    //   1, 0, 0, 0,
    //    0, 1, 0, 0,
    //    0, 0, 1, 0,
    //    0, 1, 0, -1
    //),M4.translate(0, 1))
    //console.log(m.toString())
    testCurveTransform(c, m)
  })
  test("transform4 at(0) on vanishing plane", () => {
    const c = new EllipseCurve(V3.O, V3.X, V3.Y, 0.2, 3)
    const m = M4.perspectivePlane(P3.YZ.translate(1, 0, 0))
    testCurveTransform(c, m)
  })

  function p4(p: Vector) {
    const [xw, yw, zw, w] = p.v
    return new V3(xw, yw, w)
  }

  function p42(p: Vector) {
    const [xw, yw, zw, w] = p.v
    return new V3(xw / w, yw / w, 1)
  }

  function p43(p: V3) {
    const [x, y, z] = p
    return new V3(x / z, y / z, 1)
  }

  test("transform4 fooooo", () => {
    const c = new ParabolaCurve(V3.O, V3.X, V3.Y, -2, 2)
    //const m = M4.translate(0, 1, 1);
    const m = M4.IDENTITY
    //const cm = c.transform(m);
    const p = M4.permutation4(2, 3) as M4
    const pm = m.times(p)
    const ss = arraySamples(c.tMin, c.tMax, 16).flatMap((t) =>
      ((p) => [p, p.div(p.z)])(pm.transformPoint(c.at(t))),
    )
    console.log(pm.toString())
    outputLink({
      edges: [c, parabola4Projection(pm, c.tMin, c.tMax)].map((c) =>
        edgeForCurveAndTs(c),
      ),
      drPs: ss,
      drLines: ss,
    })
  })
})
