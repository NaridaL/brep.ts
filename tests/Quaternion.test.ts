import "./manager"

import { DEG, M4, V3 } from "ts3dutils"

import { Quaternion } from "../src"

const { PI } = Math

function expectQLike(received: Quaternion, expected: Quaternion) {
  expect(received.toArray()).toBeLike(expected.toArray())
}

describe("Quaternion", () => {
  test("transformPoint", () => {
    const q = Quaternion.axis(V3.Z, PI / 2)
    expect(q.rotatePoint(V3.X)).toBeLike(V3.Y)
  })
  test("toM4", () => {
    const v = new V3(3, 4, 5).unit()
    const q = Quaternion.axis(v, 2)
    const m = M4.rotate(2, v)
    expect(q.toM4()).toBeLike(m)
  })
  test("times", () => {
    expect(
      Quaternion.axis(V3.X, 1).times(Quaternion.axis(V3.Y, 2)).toM4(),
    ).toBeLike(M4.rotateX(1).times(M4.rotateY(2)))
  })
  test("unit", () => {
    expect(Quaternion.of(1, 2, 3, 4).unit().norm()).toFuzzyEqual(1)
  })
  test("plus", () => {
    const a = Quaternion.of(1, 2, 3, 4)
    const b = Quaternion.of(5, 6, 7, 8)
    const actual = a.plus(b)
    const expected = Quaternion.of(6, 8, 10, 12)
    expectQLike(actual, expected)
  })
  test("inversed", () => {
    const q = Quaternion.of(1, 2, 3, 4)
    const actual = q.times(q.inverse())
    const expected = Quaternion.O
    expectQLike(actual, expected)
  })
  test("slerp", () => {
    const actual = Quaternion.O.slerp(Quaternion.axis(V3.Z, 100 * DEG), 0.3)
    const expected = Quaternion.axis(V3.Z, 30 * DEG)
    expectQLike(actual, expected)
  })
})
