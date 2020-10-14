import {
  suite,
  test,
  testCurve,
  testCurveISInfos,
  testISCurves,
  testISTs,
  testPointT,
} from "./manager"

import { DEG, M4, V3 } from "ts3dutils"

import { Quaternion } from "../src"
import { PI } from "../src/math"

suite("Quaternion", () => {
  test("transformPoint", (assert) => {
    const q = Quaternion.axis(V3.Z, PI / 2)
    assert.v3like(q.rotatePoint(V3.X), V3.Y)
  })
  test("toM4", (assert) => {
    const v = new V3(3, 4, 5).unit()
    const q = Quaternion.axis(v, 2)
    const m = M4.rotate(2, v)
    console.log(m.str)
    assert.mlike(q.toM4(), m)
  })
  test("times", (assert) => {
    assert.mlike(
      Quaternion.axis(V3.X, 1).times(Quaternion.axis(V3.Y, 2)).toM4(),
      M4.rotateX(1).times(M4.rotateY(2)),
    )
  })
  test("unit", (assert) => {
    assert.fuzzyEqual(Quaternion.of(1, 2, 3, 4).unit().norm(), 1)
  })
  test("plus", (assert) => {
    const a = Quaternion.of(1, 2, 3, 4)
    const b = Quaternion.of(5, 6, 7, 8)
    const actual = a.plus(b)
    const expected = Quaternion.of(6, 8, 10, 12)
    assert.push(actual.like(expected), actual, expected)
  })
  test("inversed", (assert) => {
    const q = Quaternion.of(1, 2, 3, 4)
    const actual = q.times(q.inverse())
    const expected = Quaternion.O
    assert.push(actual.like(expected), actual, expected)
  })
  test("slerp", (assert) => {
    const actual = Quaternion.O.slerp(Quaternion.axis(V3.Z, 100 * DEG), 0.3)
    const expected = Quaternion.axis(V3.Z, 30 * DEG)
    assert.push(actual.like(expected), actual, expected)
  })
})
