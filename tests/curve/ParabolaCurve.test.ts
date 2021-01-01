import { DEG, lerp, M4, V } from "ts3dutils"

import {
  testCurve,
  testCurvesColinear,
  testCurveTransform,
  testISTs,
} from "../manager"

import { Curve, Edge, P3, ParabolaCurve } from "brepts"

describe("ParabolaCurve", () => {
  const curve = new ParabolaCurve(V(1, 1), V(4, 1, -2), V(1, 10, 2))

  test("testCurve", () => {
    testCurve(curve)
    testCurve(ParabolaCurve.XY)
  })
  test("isTsWithPlane", () => {
    const plane = new P3(V(2, 7, 1).unit(), 10)
    testISTs(curve, plane, 2)
  })
  test("rightAngled", () => {
    const curveRA = curve.rightAngled()
    testCurvesColinear(curve, curveRA)
  })
  test("isColinearTo", () => {
    testCurvesColinear(ParabolaCurve.XY, ParabolaCurve.XY.scale(2, 4, 1))
  })

  test("transform4", () => {
    const c = ParabolaCurve.XY.withBounds(-1, 1).translate(1, -4, 0)
    const m = M4.product(
      M4.rotateX(90 * DEG),
      M4.perspective(45, 1, 2, 5),
      M4.rotateX(-90 * DEG),
    )
    testCurveTransform(c, m)
  })
})
