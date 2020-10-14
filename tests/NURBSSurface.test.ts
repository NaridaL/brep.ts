import { suite, suiteSurface, test } from "./manager"

import { DEG, M4, VV } from "ts3dutils"
import {
  EllipseCurve,
  NURBSSurface,
  rotateCurve,
  RotatedCurveSurface,
} from ".."

suite("NURBSSurface", () => {
  const baseCurve = EllipseCurve.forAB(2, 2)
    .rotateZ(-20 * DEG)
    .translate(4, 3)
    .rotateX(90 * DEG)
  const torusSurface = rotateCurve(
    baseCurve,
    undefined,
    undefined,
    100 * DEG,
    false,
  ) as RotatedCurveSurface
  const s = torusSurface
    .asNURBSSurface()
    .transform4(
      M4.perspective(45, 1, 1, 5).times(
        M4.rotateX(20 * DEG).translate(0, 0, -7),
      ),
    )
  const zs = [
    [1, 2, 3, 4, 5],
    [2, 3, 3, 5, 6],
    [1, 4, 7, 0, 1],
    [1, 0, 1, 0, 1],
    [1, 1, 1, 1, 1],
  ]
  const ps = [-2, -1, 0, 1, 2].flatMap((x, xi) =>
    [-2, -1, 0, 1, 2].map((y, yi) => VV(y, x, zs[yi][xi], 1)),
  )
  const s2 = new NURBSSurface(
    ps,
    [0, 1, 2, 3, 4, 5, 6, 7],
    [0, 1, 2, 3, 4, 5, 6, 7],
    2,
    2,
  )
  test("b", (assert) => {
    console.log(s.pUV(0.1, 0.1))
    console.log(s.sce)
  })
  suite("a", () => suiteSurface(s2))
})
