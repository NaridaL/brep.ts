import { testParametricSurface } from "./manager"

import { M4, TAU } from "ts3dutils"
import { RotationREqFOfZ } from ".."

import { sin, cos } from "../src/math"

describe.skip("RotationREqFOfZ", () => {
  const sinSurface = new RotationREqFOfZ(
    M4.IDENTITY,
    (z) => Math.sin(z) + 2 + 0.2 * z,
    0,
    2 * TAU,
    1,
    (z) => Math.cos(z) + 0.2,
  )

  test("testSurface", () => {
    testParametricSurface(sinSurface)
  })
})
