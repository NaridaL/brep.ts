import { outputLink, testCurveCentralProjection } from "./manager"

import { M4, V, V3 } from "ts3dutils"
import { Mesh } from "tsgl"
import { Curve, L3, edgeForCurveAndTs } from "../src"

describe("L3", () => {
  test("isInfosWithLine", () => {
    const l1 = new L3(V3.Y, V3.X)
    const res = l1.isInfosWithLine(V3.X, V(0, 2))[0]
    expect(res.tThis).toBe(1)
    expect(res.tOther).toBe(0.5)
    expect(res.p).toBeLike(V(1, 1))
  })
  test("distanceToLine", () => {
    expect(L3.X.distanceToLine(new L3(V3.Z, V3.Y))).toBe(1)
    expect(L3.X.distanceToLine(new L3(V3.Z, V3.X))).toBe(1)
  })
  test("isTsForLine", () => {
    console.log(L3.X.isInfoWithLine(new L3(V(1, 1, 0), V3.Y)).sce)
    expect(L3.X.isInfoWithLine(new L3(V(1, 1, 0), V3.Y))).toEqualV(V3.X)
  })

  function testCurvePerspectiveTransform(assert: Assert, curve: Curve, m4: M4) {
    Mesh.prototype.compile = function () {
      return this
    }
    const cubeMesh = Mesh.cube().transform(M4.scale(2).translate(-1, -1, 2))
    const cubeLines = cubeMesh.LINES.map((i) => cubeMesh.vertices[i])
    const minv = m4.inversed()
    const curveTransformed = (curve as any).transform4(m4) as Curve
    outputLink({
      edges: [edgeForCurveAndTs(curve), edgeForCurveAndTs(curveTransformed)],
      drLines: [...cubeLines, ...m4.transformedPoints(cubeLines)],
    })
  }

  test("transform4", () => {
    const m4 = M4.perspective(45, 1, 1, 5)
    const line = new L3(V3.X, V3.Z, -10, -1)
    const lineT = line.transform4(m4)
    expect(lineT.containsPoint(m4.transformPoint(line.at(-1)))).toBeTruthy()
    expect(lineT.containsPoint(m4.transformPoint(line.at(-10)))).toBeTruthy()
  })

  test("transform4 central projection", () => {
    const l = new L3(V(1, 1, 0), V(1, 1, 0.3).unit(), 0.1, 10)
    testCurveCentralProjection(l)
  })
})
