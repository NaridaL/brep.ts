import { suite, test, testParametricSurface } from "../../tests/manager"

import { M4, TAU } from "../../../ts3dutils/index"
import { RotationREqFOfZ } from "../index"

import { sin, cos } from "../math"

suite("RotationREqFOfZ", () => {
  const sinSurface = new RotationREqFOfZ(
    M4.IDENTITY,
    (z) => Math.sin(z) + 2 + 0.2 * z,
    0,
    2 * TAU,
    1,
    (z) => Math.cos(z) + 0.2,
  )

  test("testSurface", (assert) => {
    testParametricSurface(assert, sinSurface)
  })
})
