import { surfaceVolumeAndAreaTests, testLoopContainsPoint } from "../manager"

import {
  Edge,
  Face,
  L3,
  P3,
  PlaneSurface,
  PointVsFace,
  EllipseCurve,
  StraightEdge,
  edgeForCurveAndTs,
} from "brepts"

import { DEG, V, V3 } from "ts3dutils"

describe("PlaneSurface", () => {
  test("loopContainsPoint", () => {
    const loop = StraightEdge.chain(
      [V(0, 0), V(10, 0), V(10, 10), V(0, 10)],
      true,
    )
    expect(new PlaneSurface(P3.XY).loopContainsPoint(loop, V(8, 10))).toBe(
      PointVsFace.ON_EDGE,
    )
  })
  test("loopContainsPoint 2", () => {
    const loop = [
      new StraightEdge(
        new L3(V(2, 10, 0), V3.Z),
        V(2, 10, 3),
        V(2, 10, 5),
        3,
        5,
      ),
      new StraightEdge(
        new L3(V(0, 10, 5), V3.X),
        V(2, 10, 5),
        V(0, 10, 5),
        2,
        0,
      ),
      new StraightEdge(
        new L3(V(0, 10, 0), V3.Z),
        V(0, 10, 5),
        V(0, 10, 0),
        5,
        0,
      ),
      new StraightEdge(
        new L3(V(0, 10, 0), V3.X),
        V(0, 10, 0),
        V(10, 10, 0),
        0,
        10,
      ),
      new StraightEdge(
        new L3(V(10, 10, 0), V3.Z),
        V(10, 10, 0),
        V(10, 10, 5),
        0,
        5,
      ),
      new StraightEdge(
        new L3(V(0, 10, 5), V3.X),
        V(10, 10, 5),
        V(6, 10, 5),
        10,
        6,
      ),
      new StraightEdge(
        new L3(V(6, 10, 0), V(0, 0, -1)),
        V(6, 10, 5),
        V(6, 10, 3),
        -5,
        -3,
      ),
      new StraightEdge(
        new L3(V(0, 10, 3), V(-1, 0, 0)),
        V(6, 10, 3),
        V(2, 10, 3),
        -6,
        -2,
      ),
    ]
    const p = V(6, 10, 3)
    testLoopContainsPoint(
      new PlaneSurface(new P3(V(0, -1, 0), -10)),
      loop,
      p,
      PointVsFace.ON_EDGE,
    )
  })

  const triangleFace = Face.create(
    new PlaneSurface(P3.XY),
    StraightEdge.chain([V(1, 1), V(3, 2), V(2, 3)]),
  ).rotateX(10 * DEG)
  describe("triangleFace", () => surfaceVolumeAndAreaTests(triangleFace))
  describe("triangleFace.translate(2, 2, 2)", () =>
    surfaceVolumeAndAreaTests(triangleFace.translate(2, 2, 2)))
  describe("triangleFace.shearX(2, 2)", () =>
    surfaceVolumeAndAreaTests(triangleFace.shearX(2, 2)))
  describe("triangleFace.foo()", () =>
    surfaceVolumeAndAreaTests(triangleFace.foo()))

  const faceWithEllipses = Face.create(new PlaneSurface(P3.XY), [
    edgeForCurveAndTs(EllipseCurve.UNIT),
    edgeForCurveAndTs(
      EllipseCurve.circleThroughPoints(V3.X.negated(), V(0, -0.5), V3.X),
    ),
  ])
  describe("faceWithEllipses", () =>
    surfaceVolumeAndAreaTests(faceWithEllipses))
  describe("faceWithEllipses.foo()", () =>
    surfaceVolumeAndAreaTests(faceWithEllipses.foo()))
})
