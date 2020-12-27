import {
  outputLink,
  rotateEdge,
  suiteSurface,
  surfaceVolumeAndAreaTests,
  testISCurves,
  testISTs,
  testLoopContainsPoint,
} from "./manager"

import { DEG, M4, raddd, TAU, V, V3 } from "ts3dutils"
import {
  B2T,
  Edge,
  Face,
  HyperbolaCurve,
  L3,
  P3,
  PlaneSurface,
  rotateCurve,
  RotatedCurveSurface,
  EllipseCurve,
  StraightEdge,
  edgeForCurveAndTs,
} from "../src"
import { cos, PI, sin } from "../src/math"

describe("RotatedCurveSurface", () => {
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

  const torusFace = rotateEdge(edgeForCurveAndTs(baseCurve).flipped(), 40 * DEG)

  describe("torusSurface", () => suiteSurface(torusSurface))
  describe("torusSurface.scale(2)", () => suiteSurface(torusSurface.scale(2)))
  describe("torusSurface.shearX(2, 2)", () =>
    suiteSurface(torusSurface.shearX(2, 2)))

  test("getExtremePoints", () => {
    outputLink({
      mesh: torusSurface.sce + ".toMesh()",
      drPs: torusSurface.getExtremePoints(),
    })
  })
  test("getExtremePoints2", () => {
    outputLink({
      mesh: torusSurface.shearX(2, 2).sce + ".toMesh()",
      drPs: torusSurface.shearX(2, 2).getExtremePoints(),
    })
  })

  test("coplanar with curve translated up, surface down", () => {
    const torusSurface2 = rotateCurve(
      baseCurve.translate(0, 0, 2),
      undefined,
      undefined,
      100 * DEG,
      false,
    ).translate(0, 0, -2) as RotatedCurveSurface
    expect(torusSurface.isCoplanarTo(torusSurface2)).toBeTruthy()
  })
  test("is curves with plane perpendicular to rotation axis", () => {
    const cs = testISCurves(torusSurface, new P3(V3.Z, 4), 2)
    expect(cs[0]).toBeInstanceOf(EllipseCurve)
    expect(cs[1]).toBeInstanceOf(EllipseCurve)
  })
  test("is curves with plane through rotation axis", () => {
    const cs = testISCurves(torusSurface, P3.ZX.rotateZ(20 * DEG), 1)
    expect(cs[0]).toBeInstanceOf(EllipseCurve)
  })
  test("is curves with plane ", () => {
    testISCurves(torusSurface, new P3(V3.XYZ.unit(), 4), 1)
  })
  test("is curves with plane 2", () => {
    const torus = new RotatedCurveSurface(
      new EllipseCurve(
        V(2, 0, 0),
        V3.X,
        V(0, 6.123233995736766e-17, 1),
        0,
        3.141592653589793,
      ),
      M4.scale(-1, 1, 1),
      0,
      3.141592653589793,
    )
    const plane = new PlaneSurface(new P3(V3.Y, 2.5), V(-1, 0, 0), V3.Z)
    testISCurves(torus, plane, 1)
  })

  describe("torusFace", () => surfaceVolumeAndAreaTests(torusFace))
  console.log(
    "shearX22 centroid",
    torusFace.toMesh().shearX(2, 2).calcVolume().centroid,
  )
  describe("torusFace.shearX(2, 2)", () =>
    surfaceVolumeAndAreaTests(torusFace.shearX(2, 2)))

  describe("isTsWithLine", () => {
    test("vertical", () => testISTs(new L3(V(4, 4, 0), V3.Z), torusSurface, 1))
    test("vertical, 0 is", () =>
      testISTs(new L3(V(4, -4, 0), V3.Z), torusSurface, 0))
    test("in ZX plane", () =>
      testISTs(new L3(V(0, 0, 3.8), V(1, 0, -0.1).unit()), torusSurface, 3))
    test("through Z axis", () =>
      testISTs(new L3(V(4, 4, 3.8), V(1, 1, -0.1).unit()), torusSurface, 2))
    test("parallel to XY plane", () =>
      testISTs(new L3(V(4, 4, 3.8), V(1, -0.2, 0).unit()), torusSurface, 2))
    test("parallel to XY plane, 0 is", () =>
      testISTs(new L3(V(4, 4, 0), V(1, -0.2, 0).unit()), torusSurface, 0))
    test("general", () => {
      const line = new L3(V3.X, V(0, 1, 1).unit())
      testISTs(line, torusSurface, 1)
    })

    test("general 2", () => {
      const line = L3.throughPoints(V(4, 2, 3.8), V(-2, 0.2, 4))
      testISTs(line, torusSurface, 3)
    })
  })

  const baseEdge = edgeForCurveAndTs(baseCurve, 0, (2 / 3) * PI)
  const face = B2T.rotateEdges(
    [baseEdge, StraightEdge.throughPoints(baseEdge.b, baseEdge.a)],
    100 * DEG,
  ).faces.filter((x) => x.surface instanceof RotatedCurveSurface)
})
// 122.96083177519078
