import {
  b2equals,
  inDifferentSystems,
  suite,
  test,
  testBRepOp,
  testLoopContainsPoint,
} from "./manager"

import { JavaSet as CustomSet } from "javasetmap.ts"
import { DEG, eq, M4, P3XY, V, V3 } from "ts3dutils"
import * as ts3dutils from "ts3dutils"
import * as brepts from ".."
import {
  B2T,
  BezierCurve,
  ClassSerializer,
  CylinderSurface,
  Edge,
  EllipseCurve,
  EllipsoidSurface,
  Face,
  intersectionCircleLine,
  intersectionUnitCircleLine,
  L3,
  P3,
  PCurveEdge,
  PlaneSurface,
  PointVsFace,
  StraightEdge,
} from ".."

import { cos, PI, sin, sqrt } from "../src/math"

suite("NLA", () => {
  suite("isPointsWithBezier()", () =>
    inDifferentSystems((assert, m4) => {
      const ell = new EllipseCurve(
        V(-223.34900663163222, -176.63214006755936, 0),
        V(-169.5891804980124, -35.54247345835796, 0),
        V(35.54247345835796, -169.5891804980124, 0),
      )
      const bez = new BezierCurve(
        V(-267.6481190901419, -368.37017217006473, 0),
        V(563.959389388763, 94.96018577817034, 0),
        V(-1110.7787051488917, -95.8394860073627, 0),
        V(-59.14331799274822, -299.7830665459221, 0),
      )
      const isPoints = ell.isInfosWithBezier(bez).map((info) => info.p)
      assert.equal(isPoints.length, 3)
      isPoints.forEach((p) => {
        assert.ok(
          ell.containsPoint(p),
          `ell.distanceToPoint(${p}) = ${ell.distanceToPoint(p)}`,
        )
        assert.ok(
          bez.containsPoint(p),
          `bez.distanceToPoint(${p}) = ${bez.distanceToPoint(p)}`,
        )
      })
    }),
  )

  suite("EllipseCurve.getAreaInDir", () =>
    inDifferentSystems(
      (assert, m4) => {
        const k = 1
        ;[
          {
            right: V3.X,
            up: V3.Y,
            s: 0,
            t: PI,
            result: PI / 2,
            c: V(0, 4 / 3 / PI),
          },
          {
            right: V3.X,
            up: V3.Y,
            s: PI,
            t: 0,
            result: -PI / 2,
            c: V(0, 4 / 3 / PI),
          },
          {
            right: V3.X,
            up: V3.Y,
            s: -PI / 2,
            t: PI / 2,
            result: PI / 2,
            c: V(4 / 3 / PI, 0),
          },
          {
            right: V3.X,
            up: V3.Y,
            s: -PI,
            t: 0,
            result: PI / 2,
            c: V(0, -4 / 3 / PI),
          },
          // let 'down' be X
          {
            right: V3.Y,
            up: V3.X.negated(),
            s: 0,
            t: PI,
            result: PI / 2,
            c: V(0, 4 / 3 / PI),
          },
          {
            right: V3.Y,
            up: V3.X.negated(),
            s: -PI,
            t: 0,
            result: PI / 2,
            c: V(0, -4 / 3 / PI),
          },
        ].forEach((test) => {
          ;[0, 4].forEach((yDiff) => {
            const r = m4.transformVector(test.right)
            const areaFactor = m4
              .transformVector(V3.X)
              .cross(m4.transformVector(V3.Y))
              .length()
            console.log(areaFactor)
            const ell = EllipseCurve.UNIT.translate(0, yDiff, 0).transform(m4)
            const up = m4.transformVector(test.up).unit()
            const offsetArea =
              yDiff * (1 - cos(test.t) - (1 - cos(test.s))) * test.up.dot(V3.Y)
            const totalArea = test.result + offsetArea
            const expectedArea = totalArea * areaFactor
            const result = ell.getAreaInDir(r, up, test.s, test.t)
            const offsetCentroid = V((cos(test.t) + cos(test.s)) / 2, yDiff / 2)
            const movedCentroid = test.c.plus(V(0, yDiff))
            const expectedCentroid = m4.transformPoint(
              movedCentroid
                .times(test.result)
                .plus(offsetCentroid.times(offsetArea))
                .div(totalArea),
            )
            console.log(test.t, test.s, 1 - cos(test.t), 1 - cos(test.s))
            console.log(
              test.c.times(test.result).str,
              offsetCentroid.str,
              offsetArea,
              offsetCentroid.times(offsetArea).str,
              test.c.times(test.result).plus(offsetCentroid.times(offsetArea))
                .str,
              totalArea,
              expectedCentroid.str,
            )
            assert.fuzzyEqual(
              result.area,
              expectedArea,
              `yDiff: ${yDiff}, ${
                test.sce
              }, offsetArea: ${offsetArea}, expected: ${expectedArea}, ${
                areaFactor * offsetArea
              }`,
            )

            assert.fuzzyEqual(
              result.centroid.x,
              expectedCentroid.x,
              "cx " + result.centroid.x,
            )
            assert.fuzzyEqual(
              result.centroid.y,
              expectedCentroid.y,
              "cy " + result.centroid.y,
            )
            // if (!k--) throw new Error()
          })
        })
      },
      M4.IDENTITY,
      M4.rotateZ(45 * DEG),
      M4.forRows(V(1, 2, 3), V3.Y, V3.Z),
      M4.FOO.as3x3(),
    ),
  )

  test("Edge.edgesIntersects", (assert) => {
    const curve1 = BezierCurve.graphXY(2, -3, -3, 2, -2, 3)
    const curve2 = curve1.transform(M4.rotateLine(V(0.5, 0), V3.Z, PI / 2))
    const edge1 = PCurveEdge.forCurveAndTs(curve1, 0, 1)
    const edge2 = PCurveEdge.forCurveAndTs(curve2, 0, 1)
    assert.ok(Edge.edgesIntersect(edge1, edge2))
    assert.notOk(Edge.edgesIntersect(edge1, edge1.translate(10, 0, 0)))
    assert.notOk(Edge.edgesIntersect(edge1, edge2.translate(10, 0, 0)))
  })
  test("V3.areDisjoint2", (assert) => {
    const s = new CustomSet()
    const a = V(0, 2.7499999999999996, -5),
      b = V(0, 2.749999999999999, -5)
    s.canonicalizeLike(a)
    assert.ok(s.canonicalizeLike(b) == a)
  })
  test("intersectionUnitCircleLine", (assert) => {
    // y = -x + 1 => x + y = 1
    assert.deepEqual(intersectionUnitCircleLine(1, 1, 1), {
      x1: 1,
      x2: 0,
      y1: 0,
      y2: 1,
    })
  })
  test("intersectionCircleLine", (assert) => {
    // y = -x + 2 => x + y = 2
    assert.deepEqual(intersectionCircleLine(1, 1, 2, 2), {
      x1: 2,
      x2: 0,
      y1: 0,
      y2: 2,
    })
  })
  //test('EllipsoidSurface.splitOnPlaneLoop', assert => {
  //    //const es = EllipsoidSurface.UNIT
  //    const a = V3.sphere(30 * DEG, 70 * DEG), z = a.z, xy = a.lengthXY(), center = V(0, 0, z), f1 = V(a.x, a.y,
  // 0), f2 = V(-a.y, a.x) const curve = new EllipseCurve(center, f1, f2) const seamCurve =
  // EllipseCurve.UNIT.rotateX(-PI / 2) const edge = Edge.forCurveAndTs(curve, -PI, PI) assert.ok(true, `<html><a
  // style='color: #0000ff text-decoration: underline' target='blank'
  // href='viewer.html?mesh=${es.sce}.toMesh()&points=[V(-5, 1, -1)]&edges=[${edge.str}]'>view</a>`) const [front,
  // back] = EllipsoidSurface.splitOnPlaneLoop([edge], true)  assert.ok(true, `<html><a style='color: #0000ff
  // text-decoration: underline' target='blank' href='viewer.html?mesh=${es.sce}.toMesh()&points=[V(-5, 1,
  // -1)]&edges=${back.sce}'>view</a>`) console.log(front, back) const expectedFront = [] const expectedBack =
  // [Edge.forCurveAndTs(curve, -120 * DEG, 60 * DEG), Edge.forCurveAndTs(seamCurve)] },

  //test('EllipseCurve.getVolZAnd', assert => {
  //
  //	assert.equal(EllipseCurve.UNIT.getVolZAnd(V3.Z, -PI, PI).volume, 0)
  //	assert.equal(EllipseCurve.UNIT.rotateY(90 * DEG).translate(1, 0, 0).getVolZAnd(V3.Z, -PI, PI).volume, PI)
  //},
})

suite("serialization", () => {
  test("serialization", (assert) => {
    const sz = new ClassSerializer(),
      unserialize = (x) => sz.unserialize(x),
      serialize = (x) => sz.serialize(x)
    let a: any = { a: 2, b: 3 }
    assert.equal(unserialize(serialize(a)).toString(), a.toString())

    a.c = a

    const a2 = unserialize(serialize(a))
    assert.equal(a2.a, 2)
    assert.equal(a2.b, 3)
    assert.equal(a2.c, a2)

    a = [1, 2, 3]
    assert.equal(unserialize(serialize(a)).toString(), a.toString())
  })
  test("brep", (assert) => {
    const cs = new ClassSerializer()
    cs.addNamespace(ts3dutils, "ts3dutils")
    cs.addNamespace(brepts, "brepts")
    const input = B2T.box().rotateX(20 * DEG)
    const str = cs.serialize(input)
    const output = cs.unserialize(str)
    b2equals(assert, output, input)
  })
})

suite("tsgl", () => {
  test("centroid of tetrahedron O X Y Z", (assert) => {
    const centroidMesh = B2T.tetrahedron(V3.O, V3.X, V3.Y, V3.Z).toMesh()
    const centroid = centroidMesh.calcVolume().centroid
    assert.v3like(centroid, V(0.25, 0.25, 0.25))

    const centroid2 = centroidMesh.translate(2, 2).calcVolume().centroid
    assert.v3like(centroid2, V(2.25, 2.25, 0.25))
  })
})
