import { outputLink, suiteSurface } from "./manager"

import {
  addOwnProperties,
  arrayFromFunction,
  DEG,
  M4,
  Transformable,
  V,
  V3,
} from "ts3dutils"
import { Mesh } from "tsgl"
import {
  BezierCurve,
  Edge,
  EllipsoidSurface,
  P3,
  PointProjectedSurface,
  ProjectedCurveSurface,
} from ".."

test.skip("PointProjectedSurface", () => {
  const baseCurve = BezierCurve.graphXY(2, -3, -3, 2, -1, 2)
  const testSurface = new PointProjectedSurface(
    baseCurve,
    V(0, 0, 3),
    P3.XY,
    1,
    undefined,
    undefined,
    0,
    2,
  )

  // const edge = edgeForCurveAndTs(baseCurve, 0, 2)
  // const edges = [
  // 	edge,
  // 	StraightEdge.throughPoints(baseCurve.at(2), baseCurve.at(2).plus(V(0, 0, 2))),
  // 	edge.flipped().translate(0, 0, 2),
  // 	StraightEdge.throughPoints(baseCurve.at(0).plus(V(0, 0, 2)), baseCurve.at(0)),
  // ]
  // const testFace = new RotationFace(testSurface, edges)
  describe("projectedBezierSurface", () => suiteSurface(testSurface))
  describe("projectedBezierSurface.shearX(2, 2)", () =>
    suiteSurface(testSurface.shearX(2, 2)))
  describe("projectedBezierSurface.foo()", () =>
    suiteSurface(testSurface.foo()))
})
