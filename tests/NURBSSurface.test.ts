import { outputLink, suiteSurface, testISTs } from "./manager"

import { arrayFromFunction, DEG, M4, V, V3, VV } from "ts3dutils"
import {
  EllipseCurve,
  L3,
  NURBSSurface,
  rotateCurve,
  RotatedCurveSurface,
  StraightEdge,
} from "../src"

describe("NURBSSurface", () => {
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
  test("b", () => {
    console.log(s.pUV(0.1, 0.1))
    console.log(s.sce)
  })

  test("s2 guessUVForMeshPos", () => {
    //  let s2 = s
    const debugInfo = s2.debugInfo()
    const drPs: V3[] = []
    const uPointCount = s2.knotsU.length - s2.degreeU - 1
    for (let x = 0; x < uPointCount; x++) {
      for (let y = 0; y < s2.knotsV.length - s2.degreeV - 1; y++) {
        const uv = s2.pointFoot(s2.points[y * uPointCount + x].p3())
        console.log("x", x, "y", y, "uv", uv, s2.guessUVForMeshPos(x, y))
        uv && drPs.push(s2.points[y * uPointCount + x].p3(), s2.pUV(uv.x, uv.y))
      }
    }

    outputLink(
      {
        mesh: `[${s2}.toMesh()]`,
        edges: [
          ...(debugInfo && debugInfo.lines
            ? arrayFromFunction(debugInfo.lines.length / 2, (i) => {
                const a = debugInfo.lines[i * 2],
                  b = debugInfo.lines[i * 2 + 1]
                return !a.like(b) && StraightEdge.throughPoints(a, b)
              }).filter((x) => x)
            : []),
        ],
        drPs: [...(debugInfo?.points ?? [])],
        drLines: drPs,
      },
      "view",
    )
  })

  test("s2 ISTs line", () =>
    testISTs(L3.throughPoints(V(0, -2, 1.2), V(-1, -1, 3.2)), s2, 1))

  describe("a", () => suiteSurface(s2))
})
