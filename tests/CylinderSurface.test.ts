import {
  suiteSurface,
  surfaceVolumeAndAreaTests,
  testContainsCurve,
  testISCurves,
  testISTs,
  testLoopContainsPoint,
} from "./manager"

import { DEG, V, V3 } from "ts3dutils"
import {
  B2T,
  BezierCurve,
  Edge,
  Face,
  L3,
  P3,
  PCurveEdge,
  PICurve,
  PlaneFace,
  PlaneSurface,
  PointVsFace,
  RotationFace,
  CylinderSurface,
  EllipseCurve,
  EllipsoidSurface,
  StraightEdge,
  edgeForCurveAndTs,
} from ".."
import { PI, sin } from "../src/math"

describe("CylinderSurface", () => {
  test("loopContainsPoint with PICurve", () => {
    const surface = new CylinderSurface(
      new EllipseCurve(
        V(0.5, 0, -2),
        V(0.05, 0, 0),
        V(0, 0.1, 0),
        0,
        3.141592653589793,
      ),
      V(0, 0, -4),
      0,
      3.141592653589793,
      -1,
      0,
    )

    const loop = [
      new PCurveEdge(
        PICurve.forParametricStartEnd(
          new CylinderSurface(
            new EllipseCurve(
              V(0.5, 0, -2),
              V(0.05, 0, 0),
              V(0, 0.1, 0),
              0,
              3.141592653589793,
            ),
            V(0, 0, -4),
            0,
            3.141592653589793,
            -1,
            0,
          ),
          new EllipsoidSurface(
            V3.O,
            V3.X,
            V3.Y,
            V3.Z,
            0,
            3.141592653589793,
            -1.5707963267948966,
            1.5707963267948966,
          ),
          V(0, -0.7087911636061258, 0),
          V(3.141592653589793, -0.7232571387436469, 0),
          0.05,
          V(0.05, 0, 0),
          0,
          63,
        ),
        V(0.5500000000000002, 0, 0.8351646544245032),
        V(0.4499999999999999, 0, 0.8930285549745877),
        0,
        63,
        undefined,
        V(-3.42237988979182e-11, 0.0049999958882951, -3.007402222491174e-12),
        V(4.834393140838787e-9, -0.004999126576363996, -5.523525999743552e-9),
        "genseg2",
      ),
      new StraightEdge(
        new L3(V(0.45, 1.2246467991473533e-17, -2), V3.Z),
        V(0.4499999999999999, 0, 0.8930285549745877),
        V(0.44999999999999996, 0, -0.8930285549745877),
        2.8930285549745878,
        1.1069714450254122,
      ),
      new PCurveEdge(
        PICurve.forParametricStartEnd(
          new CylinderSurface(
            new EllipseCurve(
              V(0.5, 0, -2),
              V(0.05, 0, 0),
              V(0, 0.1, 0),
              0,
              3.141592653589793,
            ),
            V(0, 0, -4),
            0,
            3.141592653589793,
            -1,
            0,
          ),
          new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z),
          V(3.1445936466433873, -0.27674290222555437, 0),
          V(-0.005364365408638514, -0.29120876102138094, 0),
          0.05,
          V(-0.04999999998136267, 0.000001365186530351914, 0),
          0.060019866926763324,
          62.89271268849603,
        ),
        V(0.44999999999999996, 0, -0.8930285549745877),
        V(0.55, 0, -0.8351646544245033),
        0.060019866926763324,
        62.89271268849603,
        undefined,
        V(-2.5860397529306077e-9, 0.004999647422699275, -2.954598107221308e-9),
        V(-3.917919863240466e-9, -0.004999401467841771, 3.444920250754088e-10),
        "genseg3",
      ),
      new StraightEdge(
        new L3(V(0.55, 0, -2), V(0, 0, -1)),
        V(0.55, 0, -0.8351646544245033),
        V(0.5500000000000002, 0, 0.8351646544245032),
        -1.1648353455754967,
        -2.8351646544245033,
      ),
    ]
    const p = V(0.55, 0, 2)
    testLoopContainsPoint(surface, loop, p, PointVsFace.OUTSIDE)
  })
  const simpleCyl = new CylinderSurface(
    EllipseCurve.UNIT,
    V3.Z,
    undefined,
    undefined,
    0,
    1,
  )
  describe("simpleCyl", () => suiteSurface(simpleCyl))
  describe("simpleCyl.rotateZ(PI)", () => suiteSurface(simpleCyl.rotateZ(PI)))
  test("is curves w/ CylinderSurface", () => {
    const cyl = CylinderSurface.UNIT.scale(5, 5, 1)
    const ell = new CylinderSurface(
      new EllipseCurve(V(6, 1, 4), V(3, 1, 4), V(4, 0, 0)),
      V3.Z,
      undefined,
      undefined,
    )
    testISCurves(cyl, ell, 1)
  })
  test("is curves w/ CylinderSurface 2", () => {
    const cyl = CylinderSurface.UNIT.scale(5, 5, 1)
    const scs = new CylinderSurface(
      new EllipseCurve(
        V(-1.5, 2.5980762113533173, 0),
        V(1.25, 2.1650635094610964, 0),
        V(-2.165063509461096, 1.25, 0),
        0,
        3.141592653589793,
      ),
      V(0, 0, -1),
      undefined,
      undefined,
      -2,
      2,
    )
    testISCurves(cyl, scs, 2)
  })
  test("is curves w/ EllipsoidSurface", () => {
    const cyl = CylinderSurface.UNIT.scale(0.05, 0.05, 4).translate(
      0.5,
      0.5,
      -2,
    )
    const sphere = EllipsoidSurface.UNIT
    testISCurves(sphere, cyl, 2)
  })
  test("is curves w/ EllipsoidSurface 2", () => {
    const cyl = CylinderSurface.UNIT.scale(0.5, 0.05, 4).translate(0, 0, -2)
    const sphere = EllipsoidSurface.UNIT
    testISCurves(sphere, cyl, 2)
  })
  test("is curves w/ EllipsoidSurface 3", () => {
    const cyl = CylinderSurface.UNIT.rotateZ(PI)
      .scale(0.5, 0.05, 4)
      .translate(0.25, 0, -2)
    const sphere = EllipsoidSurface.UNIT
    testISCurves(sphere, cyl, 0)
  })
  test("is curves w/ plane", () => {
    testISCurves(CylinderSurface.UNIT, P3.XY, 1)
    testISCurves(CylinderSurface.UNIT, P3.XY.flipped(), 1)
    testISCurves(CylinderSurface.UNIT.flipped(), P3.XY, 1)
    testISCurves(CylinderSurface.UNIT.flipped(), P3.XY.flipped(), 1)

    testISCurves(
      new CylinderSurface(
        new EllipseCurve(
          V(0.5, 0.2, 0),
          V(-0.2, 2.4492935982947065e-17, 0),
          V(-2.4492935982947065e-17, -0.2, 0),
          0,
          3.141592653589793,
        ),
        V(0, 0, -1),
        undefined,
        undefined,
        -Infinity,
        Infinity,
      ),
      new PlaneSurface(new P3(V(0, -1, 0), 0), V(0, 0, -1), V(1, 0, 0)),
      1,
    )
  })
  test("intersectionLine 2", () => {
    const cylSurface = new CylinderSurface(
      new EllipseCurve(V3.O, V(8, 0, 0), V(0, 5, 0)),
      V3.Z,
      undefined,
      undefined,
    )
    const line = L3.throughPoints(V(10, 0, 0), V(-10, 2, 10))
    testISTs(line, cylSurface, 2)
  })
  test("intersectionLine", () => {
    const cylSurface = CylinderSurface.semicylinder(5)
    const line = L3.throughPoints(V(10, 0, 0), V(-10, 2, 10))
    const isPoints = cylSurface.isTsForLine(line).map(line.at, line)

    expect(isPoints.length).toBe(2)
    expect(isPoints[0].like(isPoints[1])).toBeFalsy()

    expect(cylSurface.containsPoint(isPoints[0])).toBeTruthy()
    expect(cylSurface.containsPoint(isPoints[1])).toBeTruthy()

    expect(
      line.containsPoint(isPoints[0]),
      "" + line.distanceToPoint(isPoints[0]),
    ).toBeTruthy()
    expect(
      line.containsPoint(isPoints[1]),
      "" + line.distanceToPoint(isPoints[1]),
    ).toBeTruthy()
  })
  describe("zDirVolume", () => {
    const face = B2T.extrudeEdges(
      [
        edgeForCurveAndTs(EllipseCurve.forAB(-1, 1), 0, PI),
        StraightEdge.throughPoints(V3.X, V3.X.negated()),
      ],
      P3.XY.flipped(),
      V3.Z,
      "cyl",
    ).faces.find((face) => face.surface instanceof CylinderSurface)
    const face2 = B2T.extrudeEdges(
      [
        edgeForCurveAndTs(EllipseCurve.UNIT, PI, 0),
        StraightEdge.throughPoints(V3.X, V3.X.negated()),
      ],
      P3.XY.flipped(),
      V3.Z,
      "cyl",
    ).faces.find((face) => face.surface instanceof CylinderSurface)
    const face3 = B2T.extrudeEdges(
      [
        edgeForCurveAndTs(EllipseCurve.UNIT, PI, 0).rotateY(-80 * DEG),
        StraightEdge.throughPoints(V3.X, V3.X.negated()).rotateY(-80 * DEG),
      ],
      P3.XY.flipped().rotateY(-80 * DEG),
      new V3(-10, -1, 0).unit(),
      "cyl",
    ).faces.find((face) => face.surface instanceof CylinderSurface)
    const modface = face.rotateY(-45 * DEG).translate(1, 0, 2)
    const e0 = modface.contour[0].project(new P3(modface.surface.dir, 0))
    const face4 = Face.create(modface.surface, [
      e0,
      StraightEdge.throughPoints(e0.b, modface.contour[2].a),
      modface.contour[2],
      StraightEdge.throughPoints(modface.contour[2].b, e0.a),
    ])
    surfaceVolumeAndAreaTests(face)
    surfaceVolumeAndAreaTests(face.rotateY(-45 * DEG).translate(1, 0, 2), "")
    surfaceVolumeAndAreaTests(face.rotateY(90 * DEG).translate(1, 0, 2), "")
    surfaceVolumeAndAreaTests(face2, "")
    surfaceVolumeAndAreaTests(face2.rotateY(-45 * DEG).translate(1, 0, 2), "")
    surfaceVolumeAndAreaTests(face2.rotateY(90 * DEG).translate(1, 0, 2), "")
    surfaceVolumeAndAreaTests(face3, "")
    surfaceVolumeAndAreaTests(face3.translate(1, 0, 2), "")
    surfaceVolumeAndAreaTests(face4, "")
    surfaceVolumeAndAreaTests(face4.translate(1, 0, 2), "")
  })
  test("loopContainsPoint", () => {
    const surface = new CylinderSurface(
      new EllipseCurve(V(0, 0, 0), V(8, 0, 0), V(0, 8, 0)),
      V(0, 0, -1),
      undefined,
      undefined,
    )
    const loop = [
      StraightEdge.throughPoints(
        V(1, 7.937253933193773, 4),
        V(1, 7.937253933193773, 1),
      ),
      new PCurveEdge(
        new EllipseCurve(V(0, 0, 1), V(8, 0, 0), V(0, 8, 0)),
        V(1, 7.937253933193773, 1),
        V(6, 5.291502622129181, 1),
        1.4454684956268313,
        0.7227342478134156,
        null,
        V(7.937253933193772, -0.9999999999999991, 0),
        V(5.2915026221291805, -6, 0),
      ),
      StraightEdge.throughPoints(
        V(6, 5.291502622129181, 1),
        V(6, 5.291502622129181, 4),
      ),
      new PCurveEdge(
        new EllipseCurve(V(0, 0, 4), V(8, 0, 0), V(0, 8, 0)),
        V(6, 5.291502622129181, 4),
        V(1, 7.937253933193773, 4),
        0.7227342478134156,
        1.4454684956268313,
        null,
        V(-5.2915026221291805, 6, 0),
        V(-7.937253933193772, 0.9999999999999991, 0),
      ),
    ]
    testLoopContainsPoint(surface, loop, V(8, 0, 0), PointVsFace.OUTSIDE)
    testLoopContainsPoint(
      surface,
      loop,
      V(1, 7.937253933193773, 3),
      PointVsFace.ON_EDGE,
    )
  })

  test("containsCurve", () => {
    const surface = new CylinderSurface(
      new EllipseCurve(V3.O, V3.X, V3.Y, 0, 3.141592653589793),
      V3.Z,
      0,
      3.141592653589793,
      0,
      2,
    )
    const curve = new EllipseCurve(
      V(0, 0, 2),
      V(-6.123233995736766e-17, -1, 0),
      V(-1, 6.123233995736766e-17, 0),
      0,
      3.141592653589793,
    )
    testContainsCurve(assert, surface, curve)
  })

  const surface = CylinderSurface.UNIT
  // loop which is 1 high and goes around a quarter of the cylinder
  const loop = [
    StraightEdge.throughPoints(V(1, 0, 1), V(1, 0, 0)),
    edgeForCurveAndTs(EllipseCurve.UNIT, 0, PI / 2),
    StraightEdge.throughPoints(V(0, 1, 0), V(0, 1, 1)),
    edgeForCurveAndTs(EllipseCurve.UNIT.translate(0, 0, 1), PI / 2, 0),
  ]
  const cylinderFace2 = Face.create(surface, loop)
  describe("cylinderFace2", () => surfaceVolumeAndAreaTests(cylinderFace2))
  describe("cylinderFace2.foo()", () =>
    surfaceVolumeAndAreaTests(cylinderFace2.foo()))
})
