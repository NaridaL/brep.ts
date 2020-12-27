import { inDifferentSystems } from "./manager"

import { M4, P3XY, snap0, V, V3 } from "ts3dutils"
import { L3, P3 } from "../src"
import { sign } from "../src/math"

describe("P3", () => {
  test("projectedVector", () => {
    const p = new P3(V3.Z, 2)
    expect(V(1, 1, 0).like(p.projectedVector(V(1, 1, 1)))).toBeTruthy()
  })
  test("transform IDENTITY", () => {
    const p = new P3(V3.Z, 2)
    expect(P3.XY.like(P3.XY.transform(M4.IDENTITY))).toBeTruthy()
  })
  test("intersectionWithPlane", () => {
    expect(P3.XY.intersectionWithPlane(P3.ZX).isColinearTo(L3.X)).toBeTruthy()
    expect(P3.ZX.intersectionWithPlane(P3.XY).isColinearTo(L3.X)).toBeTruthy()
    expect(P3.ZX.intersectionWithPlane(P3.XY).isColinearTo(L3.Y)).toBeFalsy()
  })

  describe("transform", () =>
    inDifferentSystems(
      () => {
        const plane = new P3(V3.X, 2)
        const m4 = M4.mirror(P3.YZ)
        const planet = plane.transform(m4)
        console.log(planet.sce)
        const ps = [plane.anchor, V(2, 2, 2), V(3, 0, 0)]
        ps.forEach((p) => {
          const pt = m4.transformPoint(p)
          expect(sign(snap0(planet.distanceToPointSigned(pt)))).toBe(
            sign(snap0(plane.distanceToPointSigned(p))),
          )
        })
      },
      M4.perspective(45, 1, 1, 5),
      M4.mirror(P3.YZ),
      M4.translate(1, 0, 0),
      M4.FOO,
      M4.IDENTITY,
    ))

  test("M4.projection", () => {
    const plane = new P3(V(1, 2, 3).unit(), 5)
    const proj = M4.project(plane)
    console.log(proj.transformPoint(V(2, 4, 6)))
    expect(proj.transformPoint(V(2, 4, 6))).toBeLike(plane.anchor)
    expect(proj.transformVector(V(2, 4, 6))).toBeLike(V3.O)
    const p2 = V(3, 5, 22)
    expect(
      proj.transformPoint(p2).minus(p2).isParallelTo(plane.normal1),
    ).toBeTruthy()
    expect(plane.containsPoint(proj.transformPoint(p2))).toBeTruthy()
    expect(
      proj.transformVector(p2).minus(p2).isParallelTo(plane.normal1),
    ).toBeTruthy()
    expect(
      proj.transformVector(p2).isPerpendicularTo(plane.normal1),
    ).toBeTruthy()
  })
  test("M4.projection 2", () => {
    ;[V(1, 1, 1), V(0, 0, -1)].forEach((dir) => {
      const plane = new P3(V(1, 2, 3).unit(), 5)
      const proj = M4.project(plane, dir)
      const p2 = V(3, 5, 22)
      expect(proj.transformPoint(p2).minus(p2).isParallelTo(dir)).toBeTruthy()
      expect(plane.containsPoint(proj.transformPoint(p2))).toBeTruthy()
      expect(proj.transformVector(p2).minus(p2).isParallelTo(dir)).toBeTruthy()
      expect(
        proj.transformVector(p2).isPerpendicularTo(plane.normal1),
      ).toBeTruthy()
      console.log(proj.transformPoint(p2).sce)
      console.log(proj.toString())
    })
  })
  test("M4.projectPlanePoint()", () => {
    const m4 = M4.projectPlanePoint(V3.Z.negated(), P3XY)
    expect(m4.transformPoint(V(4, 0, 1))).toBeLike(V(2, 0, 0))
    expect(m4.transformPoint(V(4, 8, 1))).toBeLike(V(2, 4, 0))
    expect(m4.transformPoint(V(4, 8, 2))).toBeLike(V(4 / 3, 8 / 3, 0))
    expect(
      M4.projectPlanePoint(
        M4.FOO.transformPoint(V3.Z.negated()),
        P3.XY.transform(M4.FOO),
      ).normalized2(),
    ).toBeLike(M4.product(M4.FOO, m4, M4.FOO_INV).normalized2())
  })
})
