import { b2equals, outputLink } from "./manager"

import * as chroma from "chroma.ts"
import * as fs from "fs"
import { AABB, DEG, M4, TAU, V } from "ts3dutils"
import { Mesh } from "tsgl"

import {
  B2T,
  BRep,
  EllipseCurve,
  StraightEdge,
  edgeForCurveAndTs,
  PlaneFace,
  P3,
  L3,
} from ".."

describe("BRep generators", () => {
  test("rotStep w/ straight edges", () => {
    const actual = B2T.rotStep(
      StraightEdge.chain([V(2, 0, 2), V(4, 0, 2), V(4, 0, 4)]),
      TAU,
      5,
    )
    expect(actual).toMatchBRepSnapshot()
  })

  test("rotStep w/ circle", () => {
    const actual = B2T.rotStep(
      [
        edgeForCurveAndTs(
          EllipseCurve.semicircle(2, V(3, 0)).rotateX(90 * DEG),
        ),
        StraightEdge.throughPoints(V(1, 0, 0), V(5, 0, 0)),
      ],
      [0.1, 1, 2],
    )
    expect(actual).toMatchBRepSnapshot()
  })

  test("rotateEdges w/ straight edges", () => {
    const actual = B2T.rotateEdges(
      StraightEdge.chain([V(2, 0, 2), V(4, 0, 2), V(4, 0, 4)]),
      TAU,
    )
    expect(actual).toMatchBRepSnapshot()
  })

  test("pyramidEdges", () => {
    const actual = B2T.pyramidEdges(
      [
        edgeForCurveAndTs(EllipseCurve.semicircle(2, V(3, 0))),
        StraightEdge.throughPoints(V(1, 0, 0), V(5, 0, 0)),
      ],
      V(0, 0, 4),
    )
    expect(actual).toMatchBRepSnapshot()
  })

  test.only("torus", () => {
    const actual = B2T.torus(1, 2)
    expect(actual).toMatchBRepSnapshot()
  })

  test("fromBPT", () => {
    const bpt = fs.readFileSync(__dirname + "/fixtures/teapotrim.bpt", "utf8")
    const actual = B2T.fromBPT(bpt)
    const expected = BRep.EMPTY
    actual[0].toMesh()
    expect(actual).toMatchBRepSnapshot()
  })

  test("tetrahedron", () => {
    const actual = B2T.tetrahedron(
      V(5, 5, 5),
      V(5, 5, -5),
      V(9, 11, 1),
      V(1, 9, 0),
    )
    expect(actual).toMatchBRepSnapshot()
  })

  test("chroma lab rgb", () => {
    const r = 8 * 2
    const drPs = Mesh.box(r, r, r).vertices.map((p) => {
      const c = chroma.gl(...p.toArray(), 1)
      return { p: V(c.lab()), color: c.hex() }
    })
    const aabb = new AABB(V(0, -87, -108), V(100, 99, 95))

    const obb = M4.forSys(
      V(-0.433, -0.665, -0.609).toLength(124.525),
      V(-0.354, 0.746, -0.564).toLength(250.71),
      V(0.829, -0.028, -0.558).toLength(86.463),
      V(99.707, -28.17, 152.469),
    )
    console.log("aabb vol", aabb.volume())
    console.log("obb vol", obb.determinant())
    //; cX: (45.848, -110.990, 76.664); cY: (10.883, 158.919, 11.178); cZ: (171.393, -30.598, 104.189)

    outputLink({ drPs, boxes: [obb, aabb.getM4()] })

    const wedgeVol = (3602639 / 0x1_00_00_00) * aabb.volume()
    console.log("wedgeVol " + wedgeVol)
    console.log("wedgeVol/aabb vol " + wedgeVol / aabb.volume())
    console.log("wedgeVol/obb vol " + wedgeVol / obb.determinant())
    //let c = 0,m=aabb.getM4()
    //for(let i=0;i<0x1_00_00_00;i++){
    //    const [x,y,z]=chroma.num(i).gl()
    //    const [l,a,b]=m.transformPoint(new V3(x,y,z))
    //    if (!chroma.lab(l,a,b).clipped())c++
    //
    //}
    //console.log("c",c)
  })
})
