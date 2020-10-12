import { inDifferentSystems, suite, test } from "./manager"

import { M4, P3XY, snap0, V, V3 } from "ts3dutils"
import { L3, P3 } from ".."
import { sign } from "../src/math"

suite("P3", () => {
  test("projectedVector", (assert) => {
    const p = new P3(V3.Z, 2)
    assert.ok(V(1, 1, 0).like(p.projectedVector(V(1, 1, 1))))
  })
  test("transform IDENTITY", (assert) => {
    const p = new P3(V3.Z, 2)
    assert.ok(
      P3.XY.like(P3.XY.transform(M4.IDENTITY)),
      P3.XY.transform(M4.IDENTITY).sce,
    )
  })
  test("intersectionWithPlane", (assert) => {
    assert.ok(P3.XY.intersectionWithPlane(P3.ZX).isColinearTo(L3.X))
    assert.ok(P3.ZX.intersectionWithPlane(P3.XY).isColinearTo(L3.X))
    assert.notOk(P3.ZX.intersectionWithPlane(P3.XY).isColinearTo(L3.Y))
  })

  suite("transform", () =>
    inDifferentSystems(
      (assert) => {
        const plane = new P3(V3.X, 2)
        const m4 = M4.mirror(P3.YZ)
        const planet = plane.transform(m4)
        console.log(planet.sce)
        const ps = [plane.anchor, V(2, 2, 2), V(3, 0, 0)]
        ps.forEach((p) => {
          const pt = m4.transformPoint(p)
          assert.equal(
            sign(snap0(planet.distanceToPointSigned(pt))),
            sign(snap0(plane.distanceToPointSigned(p))),
          )
        })
      },
      M4.perspective(45, 1, 1, 5),
      M4.mirror(P3.YZ),
      M4.translate(1, 0, 0),
      M4.FOO,
      M4.IDENTITY,
    ),
  )

  test("M4.projection", (assert) => {
    const plane = new P3(V(1, 2, 3).unit(), 5)
    const proj = M4.project(plane)
    console.log(proj.transformPoint(V(2, 4, 6)))
    assert.v3like(proj.transformPoint(V(2, 4, 6)), plane.anchor)
    assert.v3like(proj.transformVector(V(2, 4, 6)), V3.O)
    const p2 = V(3, 5, 22)
    assert.ok(proj.transformPoint(p2).minus(p2).isParallelTo(plane.normal1))
    assert.ok(plane.containsPoint(proj.transformPoint(p2)))
    assert.ok(proj.transformVector(p2).minus(p2).isParallelTo(plane.normal1))
    assert.ok(proj.transformVector(p2).isPerpendicularTo(plane.normal1))
  })
  test("M4.projection 2", (assert) => {
    ;[V(1, 1, 1), V(0, 0, -1)].forEach((dir) => {
      const plane = new P3(V(1, 2, 3).unit(), 5)
      const proj = M4.project(plane, dir)
      const p2 = V(3, 5, 22)
      assert.ok(proj.transformPoint(p2).minus(p2).isParallelTo(dir))
      assert.ok(plane.containsPoint(proj.transformPoint(p2)))
      assert.ok(proj.transformVector(p2).minus(p2).isParallelTo(dir))
      assert.ok(proj.transformVector(p2).isPerpendicularTo(plane.normal1))
      console.log(proj.transformPoint(p2).sce)
      console.log(proj.str)
    })
  })
  test("M4.projectPlanePoint()", (assert) => {
    const m4 = M4.projectPlanePoint(V3.Z.negated(), P3XY)
    assert.v3like(m4.transformPoint(V(4, 0, 1)), V(2, 0, 0))
    assert.v3like(m4.transformPoint(V(4, 8, 1)), V(2, 4, 0))
    assert.v3like(m4.transformPoint(V(4, 8, 2)), V(4 / 3, 8 / 3, 0))
    assert.m4equiv(
      M4.projectPlanePoint(
        M4.FOO.transformPoint(V3.Z.negated()),
        P3.XY.transform(M4.FOO),
      ),
      M4.product(M4.FOO, m4, M4.BAR),
    )
  })
})
