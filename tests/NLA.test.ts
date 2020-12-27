import { b2equals, inDifferentSystems } from "./manager"

import { JavaSet as CustomSet } from "javasetmap.ts"
import { DEG, M4, V, V3 } from "ts3dutils"
import * as ts3dutils from "ts3dutils"
import * as brepts from "../src"
import {
  B2T,
  BezierCurve,
  ClassSerializer,
  Edge,
  edgeForCurveAndTs,
  EllipseCurve,
  intersectionCircleLine,
  intersectionUnitCircleLine,
} from "../src"

import { cos, PI } from "../src/math"

describe("NLA", () => {
  describe("isPointsWithBezier()", () =>
    inDifferentSystems((m4) => {
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
      expect(isPoints.length).toBe(3)
      isPoints.forEach((p) => {
        expect(ell.containsPoint(p)).toBeTruthy()
        expect(bez.containsPoint(p)).toBeTruthy()
      })
    }))

  describe("EllipseCurve.getAreaInDir", () =>
    inDifferentSystems(
      (m4) => {
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
              test.c.times(test.result).toString(),
              offsetCentroid.toString(),
              offsetArea,
              offsetCentroid.times(offsetArea).toString(),
              test.c
                .times(test.result)
                .plus(offsetCentroid.times(offsetArea))
                .toString(),
              totalArea,
              expectedCentroid.toString(),
            )
            expect(result.area).toFuzzyEqual(expectedArea)

            expect(result.centroid.x).toFuzzyEqual(expectedCentroid.x)
            expect(result.centroid.y).toFuzzyEqual(expectedCentroid.y)
            // if (!k--) throw new Error()
          })
        })
      },
      M4.IDENTITY,
      M4.rotateZ(45 * DEG),
      M4.forRows(V(1, 2, 3), V3.Y, V3.Z),
      M4.FOO.as3x3(),
    ))

  test("Edge.edgesIntersects", () => {
    const curve1 = BezierCurve.graphXY(2, -3, -3, 2, -2, 3)
    const curve2 = curve1.transform(M4.rotateLine(V(0.5, 0), V3.Z, PI / 2))
    const edge1 = edgeForCurveAndTs(curve1, 0, 1)
    const edge2 = edgeForCurveAndTs(curve2, 0, 1)
    expect(Edge.edgesIntersect(edge1, edge2)).toBeTruthy()
    expect(Edge.edgesIntersect(edge1, edge1.translate(10, 0, 0))).toBeFalsy()
    expect(Edge.edgesIntersect(edge1, edge2.translate(10, 0, 0))).toBeFalsy()
  })
  test("V3.areDisjoint2", () => {
    const s = new CustomSet()
    const a = V(0, 2.7499999999999996, -5),
      b = V(0, 2.749999999999999, -5)
    s.canonicalizeLike(a)
    expect(s.canonicalizeLike(b) == a).toBeTruthy()
  })
  test("intersectionUnitCircleLine", () => {
    // y = -x + 1 => x + y = 1
    expect(intersectionUnitCircleLine(1, 1, 1)).toEqual({
      x1: 1,
      x2: 0,
      y1: 0,
      y2: 1,
    })
  })
  test("intersectionCircleLine", () => {
    // y = -x + 2 => x + y = 2
    expect(intersectionCircleLine(1, 1, 2, 2)).toEqual({
      x1: 2,
      x2: 0,
      y1: 0,
      y2: 2,
    })
  })
  //test('EllipsoidSurface.splitOnPlaneLoop', () => {
  //    //const es = EllipsoidSurface.UNIT
  //    const a = V3.sphere(30 * DEG, 70 * DEG), z = a.z, xy = a.lengthXY(), center = V(0, 0, z), f1 = V(a.x, a.y,
  // 0), f2 = V(-a.y, a.x) const curve = new EllipseCurve(center, f1, f2) const seamCurve =
  // EllipseCurve.UNIT.rotateX(-PI / 2) const edge = edgeForCurveAndTs(curve, -PI, PI) assert.ok(true, `<html><a
  // style='color: #0000ff text-decoration: underline' target='blank'
  // href='viewer.html?mesh=${es.sce}.toMesh()&points=[V(-5, 1, -1)]&edges=[${edge.toString()}]'>view</a>`) const [front,
  // back] = EllipsoidSurface.splitOnPlaneLoop([edge], true)  assert.ok(true, `<html><a style='color: #0000ff
  // text-decoration: underline' target='blank' href='viewer.html?mesh=${es.sce}.toMesh()&points=[V(-5, 1,
  // -1)]&edges=${back.sce}'>view</a>`) console.log(front, back) const expectedFront = [] const expectedBack =
  // [edgeForCurveAndTs(curve, -120 * DEG, 60 * DEG), edgeForCurveAndTs(seamCurve)] },

  //test('EllipseCurve.getVolZAnd', () => {
  //
  //	assert.equal(EllipseCurve.UNIT.getVolZAnd(V3.Z, -PI, PI).volume, 0)
  //	assert.equal(EllipseCurve.UNIT.rotateY(90 * DEG).translate(1, 0, 0).getVolZAnd(V3.Z, -PI, PI).volume, PI)
  //},
})

describe("serialization", () => {
  test("serialization", () => {
    const sz = new ClassSerializer(),
      unserialize = (x) => sz.unserialize(x),
      serialize = (x) => sz.serialize(x)
    let a: any = { a: 2, b: 3 }
    expect(unserialize(serialize(a)).toString()).toBe(a.toString())

    a.c = a

    const a2 = unserialize(serialize(a))
    expect(a2.a).toBe(2)
    expect(a2.b).toBe(3)
    expect(a2.c).toBe(a2)

    a = [1, 2, 3]
    expect(unserialize(serialize(a)).toString()).toBe(a.toString())
  })
  test("brep", () => {
    const cs = new ClassSerializer()
    cs.addNamespace(ts3dutils, "ts3dutils")
    cs.addNamespace(brepts, "brepts")
    const input = B2T.box().rotateX(20 * DEG)
    const str = cs.serialize(input)
    const output = cs.unserialize(str)
    b2equals(output, input)
  })
})

describe("tsgl", () => {
  test("centroid of tetrahedron O X Y Z", () => {
    const centroidMesh = B2T.tetrahedron(V3.O, V3.X, V3.Y, V3.Z).toMesh()
    const centroid = centroidMesh.calcVolume().centroid
    expect(centroid).toBeLike(V(0.25, 0.25, 0.25))

    const centroid2 = centroidMesh.translate(2, 2).calcVolume().centroid
    expect(centroid2).toBeLike(V(2.25, 2.25, 0.25))
  })
})
