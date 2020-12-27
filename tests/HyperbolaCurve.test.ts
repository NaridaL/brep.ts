import { testCurve, testCurveTransform, testISTs } from "./manager"

import { DEG, M4, V } from "ts3dutils"
import { HyperbolaCurve, intersectionUnitHyperbolaLine, P3 } from "../src"

import { sqrt } from "../src/math"

describe("HyperbolaCurve", () => {
  test("testCurve", () => {
    testCurve(HyperbolaCurve.XY)
  })
  test("testCurve 2", () => {
    const hbSheared = HyperbolaCurve.XY.shearX(2, 3)
    expect(hbSheared.isOrthogonal()).toBeFalsy()
    const hbShearedRA = hbSheared.rightAngled()
    expect(hbShearedRA.isOrthogonal()).toBeTruthy()
    expect(hbSheared.isColinearTo(hbShearedRA)).toBeTruthy()
    testCurve(hbShearedRA)

    expect(intersectionUnitHyperbolaLine(1, 0, 2)).toEqual({
      x1: 2,
      y1: sqrt(3),
      x2: 2,
      y2: -sqrt(3),
    })
  })
  test("isTsWithPlane", () => {
    testISTs(HyperbolaCurve.XY, P3.YZ, 0)
    testISTs(HyperbolaCurve.XY, P3.YZ.translate(1), 1)
    testISTs(HyperbolaCurve.XY, P3.YZ.translate(2), 2)
    testISTs(HyperbolaCurve.XY, new P3(V(1, 2).unit(), 2), 1)
    testISTs(HyperbolaCurve.XY, new P3(V(1, 2).unit(), 2).flipped(), 1)
    testISTs(HyperbolaCurve.XY, new P3(V(1, -2).unit(), 2), 1)

    testISTs(HyperbolaCurve.XY, new P3(V(1, 1).unit(), 2), 1)
    testISTs(HyperbolaCurve.XY, new P3(V(1, 1).unit(), 2).flipped(), 1)

    testISTs(HyperbolaCurve.XY, new P3(V(2, 1).unit(), 2), 2)
    testISTs(HyperbolaCurve.XY, new P3(V(2, 1).unit(), 2).flipped(), 2)
    testISTs(HyperbolaCurve.XY, new P3(V(2, 1).unit(), 0.85), 2)
    testISTs(HyperbolaCurve.XY, new P3(V(2, 1).unit(), 0.5), 0)
  })
  test("isTsWithPlane no IS with planes X < 0", () => {
    testISTs(HyperbolaCurve.XY, P3.YZ, 0)
    testISTs(HyperbolaCurve.XY, P3.YZ.translate(-2), 0)

    testISTs(HyperbolaCurve.XY, new P3(V(1, 2).unit(), -2).flipped(), 1)
    testISTs(HyperbolaCurve.XY, new P3(V(1, -2).unit(), -2), 1)

    testISTs(HyperbolaCurve.XY, new P3(V(1, 1).unit(), -2), 0)
    testISTs(HyperbolaCurve.XY, new P3(V(1, 1).unit(), 0), 0)

    testISTs(HyperbolaCurve.XY, new P3(V(2, 1).unit(), -2), 0)
    testISTs(HyperbolaCurve.XY, new P3(V(2, -1).unit(), -2), 0)
  })

  test("transform4", () => {
    const c = HyperbolaCurve.XY.withBounds(-1, 1).translate(1, -4, 0)
    const m = M4.product(
      M4.rotateX(90 * DEG),
      M4.perspective(45, 1, 2, 5),
      M4.rotateX(-90 * DEG),
    )
    testCurveTransform(c, m)
  })
})
