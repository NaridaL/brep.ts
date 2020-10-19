import {
  b2equals,
  FONT_PATH,
  outputLink,
  testBRepAnd,
  testBRepOp,
} from "./manager"

import { AABB, DEG, isCCW, M4, NLA_PRECISION, TAU, V, V3 } from "ts3dutils"
import {
  ALONG_EDGE_OR_PLANE,
  B2T,
  BezierCurve,
  BRep,
  ConicSurface,
  COPLANAR_OPPOSITE,
  COPLANAR_SAME,
  Edge,
  edgeForCurveAndTs,
  edgeNgon,
  edgePathFromSVG,
  edgeStar,
  EllipseCurve,
  INSIDE,
  L3,
  OUTSIDE,
  P3,
  PCurveEdge,
  PlaneFace,
  PlaneSurface,
  PointVsFace,
  RotatedCurveSurface,
  RotationFace,
  splitsVolumeEnclosingCone,
  splitsVolumeEnclosingCone2,
  splitsVolumeEnclosingFaces,
  splitsVolumeEnclosingFacesP,
  StraightEdge,
} from ".."

import { PI } from "../src/math"

describe("Edge", () => {
  test("like", () => {
    const edge1 = new PCurveEdge(
      new EllipseCurve(
        V(2, 0, 0),
        V3.X,
        V(0, 6.123233995736766e-17, 1),
        0,
        3.141592653589793,
      ),
      V(1, 7.498798913309288e-33, 1.2246467991473532e-16),
      V(3, 0, 0),
      3.141592653589793,
      0,
      undefined,
      V(1.2246467991473532e-16, 6.123233995736766e-17, 1),
      V(0, -6.123233995736766e-17, -1),
      "undefined.rotateX(1.5707963267948966)",
    )
    const edge2 = new PCurveEdge(
      new EllipseCurve(
        V(2, 0, 0),
        V(-1, 0, 0),
        V(0, -6.123233995736766e-17, -1),
        0,
        3.141592653589793,
      ),
      V3.X,
      V(3, -7.498798913309288e-33, -1.2246467991473532e-16),
      0,
      3.141592653589793,
      undefined,
      V(0, -6.123233995736766e-17, -1),
      V(1.2246467991473532e-16, 6.123233995736766e-17, 1),
      "undefined.rotateX(1.5707963267948966)",
    )
    outputLink({
      edges: [edge1, edge2],
    })
    expect(edge1.like(edge2)).toBeFalsy()
  })
})

describe("Face", () => {
  test("Face not like", () => {
    const face1 = new RotationFace( // region face1...
      new RotatedCurveSurface(
        new EllipseCurve(
          V(2, 0, 0),
          V(-1, 0, 0),
          V(0, -6.123233995736766e-17, -1),
          0,
          3.141592653589793,
        ),
        M4.forSys(
          V(-1, 1.2246467991473532e-16, 0),
          V(-1.2246467991473532e-16, -1, 0),
          V3.Z,
        ),
        0,
        3.141592653589793,
        0,
        3.141592653589793,
      ),
      [
        new PCurveEdge(
          new EllipseCurve(
            V(-2, 2.4492935982947064e-16, 0),
            V(1, -1.2246467991473532e-16, 0),
            V(7.498798913309288e-33, 6.123233995736766e-17, -1),
            0,
            3.141592653589793,
          ),
          V(-3, 3.6739403974420594e-16, -1.2246467991473532e-16),
          V(-1, 1.2246467991473532e-16, 0),
          3.141592653589793,
          0,
          undefined,
          V(1.2246467991473532e-16, 6.123233995736765e-17, -1),
          V(-7.498798913309288e-33, -6.123233995736766e-17, 1),
          "undefined.rotateX(1.5707963267948966)undefined",
        ),
        new PCurveEdge(
          new EllipseCurve(
            V3.O,
            V(-1, 1.2246467991473532e-16, 0),
            V(-1.2246467991473532e-16, -1, 0),
            0,
            3.141592653589793,
          ),
          V(-1, 1.2246467991473532e-16, 0),
          V3.X,
          0,
          3.141592653589793,
          undefined,
          V(-1.2246467991473532e-16, -1, 0),
          V(2.4492935982947064e-16, 1, 0),
          "torus2rib1",
        ),
        new PCurveEdge(
          new EllipseCurve(
            V(2, 0, 0),
            V(-1, 0, 0),
            V(0, -6.123233995736766e-17, -1),
            0,
            3.141592653589793,
          ),
          V3.X,
          V(3, -7.498798913309288e-33, -1.2246467991473532e-16),
          0,
          3.141592653589793,
          undefined,
          V(0, -6.123233995736766e-17, -1),
          V(1.2246467991473532e-16, 6.123233995736766e-17, 1),
          "undefined.rotateX(1.5707963267948966)",
        ),
        new PCurveEdge(
          new EllipseCurve(
            V3.O,
            V(-3, 3.6739403974420594e-16, 0),
            V(-3.6739403974420594e-16, -3, 0),
            0,
            3.141592653589793,
          ),
          V(3, 0, 0),
          V(-3, 3.6739403974420594e-16, 0),
          3.141592653589793,
          0,
          undefined,
          V(-7.347880794884119e-16, -3, 0),
          V(3.6739403974420594e-16, 3, 0),
          "torus2rib0",
        ),
      ],
      [],
    )
    // endregion
    const face2 = new RotationFace( //region face2...
      new RotatedCurveSurface(
        new EllipseCurve(
          V(2, 0, 0),
          V3.X,
          V(0, 6.123233995736766e-17, 1),
          0,
          3.141592653589793,
        ),
        M4.IDENTITY,
        0,
        3.141592653589793,
        0,
        3.141592653589793,
      ),
      [
        new PCurveEdge(
          new EllipseCurve(
            V(2, 0, 0),
            V3.X,
            V(0, 6.123233995736766e-17, 1),
            0,
            3.141592653589793,
          ),
          V(1, 7.498798913309288e-33, 1.2246467991473532e-16),
          V(3, 0, 0),
          3.141592653589793,
          0,
          undefined,
          V(1.2246467991473532e-16, 6.123233995736766e-17, 1),
          V(0, -6.123233995736766e-17, -1),
          "undefined.rotateX(1.5707963267948966)",
        ),
        new PCurveEdge(
          new EllipseCurve(V3.O, V(3, 0, 0), V(0, 3, 0), 0, 3.141592653589793),
          V(3, 0, 0),
          V(-3, 3.6739403974420594e-16, 0),
          0,
          3.141592653589793,
          undefined,
          V(0, 3, 0),
          V(-3.6739403974420594e-16, -3, 0),
          "torus2rib0",
        ),
        new PCurveEdge(
          new EllipseCurve(
            V(-2, 2.4492935982947064e-16, 0),
            V(-1, 1.2246467991473532e-16, 0),
            V(-7.498798913309288e-33, -6.123233995736766e-17, 1),
            0,
            3.141592653589793,
          ),
          V(-3, 3.6739403974420594e-16, 0),
          V(-1, 1.2246467991473532e-16, 1.2246467991473532e-16),
          0,
          3.141592653589793,
          undefined,
          V(-7.498798913309288e-33, -6.123233995736766e-17, 1),
          V(1.2246467991473532e-16, 6.123233995736765e-17, -1),
          "undefined.rotateX(1.5707963267948966)undefined",
        ),
        new PCurveEdge(
          new EllipseCurve(V3.O, V3.X, V3.Y, 0, 3.141592653589793),
          V(-1, 1.2246467991473532e-16, 0),
          V3.X,
          3.141592653589793,
          0,
          undefined,
          V(1.2246467991473532e-16, 1, 0),
          V(0, -1, 0),
          "torus2rib1",
        ),
      ],
      [],
    )
    //endregion
    outputLink({
      face: [face1, face2],
    })
    expect(face1.likeFace(face2)).toBeFalsy()
  })

  test("equals", () => {
    const a = PlaneFace.forVertices(P3.XY, [
      V(0, 0, 0),
      V(10, 0, 0),
      V(10, 10, 0),
      V(0, 10, 0),
    ])
    const b = PlaneFace.forVertices(P3.XY, [
      V(0, 10, 0),
      V(0, 0, 0),
      V(10, 0, 0),
      V(10, 10, 0),
    ])
    const c = PlaneFace.forVertices(
      new P3(V(0, 0, -1), 0),
      [V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)].slice().reverse(),
    )
    expect(a.equals(a)).toBeTruthy()

    expect(a.equals(b)).toBeTruthy()
    expect(b.equals(a)).toBeTruthy()

    expect(a.equals(c)).toBeFalsy()
    expect(c.equals(a)).toBeFalsy()

    expect(b.equals(c)).toBeFalsy()
    expect(c.equals(b)).toBeFalsy()
  })

  test("transform", () => {
    const loop = edgePathFromSVG("m 2 0, 1 1, -2 1 Z")
    const face = new PlaneFace(P3.XY, loop)
    const faceMirrored = face.mirroredX()
    expect(faceMirrored.mirroredX().likeFace(face)).toBeTruthy()
  })

  test("containsPoint", () => {
    let a = PlaneFace.forVertices(P3.XY, [
      V(0, 0, 0),
      V(10, 0, 0),
      V(10, 10, 0),
      V(0, 10, 0),
    ])
    expect(a.containsPoint(V(5, 5, 0))).toBeTruthy()
    expect(a.containsPoint(V(11, 5, 0))).toBeFalsy()

    let b = PlaneFace.forVertices(P3.XY, [
      V(0, 10, 0),
      V(0, 0, 0),
      V(5, 0, 0),
      V(6, 5, 0),
    ])
    expect(b.containsPoint(V(2, 5, 0))).toBeTruthy()

    let c = PlaneFace.forVertices(P3.XY, [
      V(0, 10, 0),
      V(0, 0, 0),
      V(5, 0, 0),
      V(6, 5, 0),
      V(10, 5, 0),
    ])
    expect(c.containsPoint(V(2, 5, 0))).toBeTruthy()

    a = a.rotateZ(30 * DEG)
    const m = M4.rotateZ(30 * DEG)
    expect(a.containsPoint(m.transformPoint(V(5, 5, 0)))).toBeTruthy()
    expect(a.containsPoint(m.transformPoint(V(-5, 5, 0)))).toBeFalsy()

    b = b.rotateZ(30 * DEG)
    expect(b.containsPoint(m.transformPoint(V(2, 5, 0)))).toBeTruthy()

    c = c.rotateZ(30 * DEG)
    expect(c.containsPoint(m.transformPoint(V(2, 5, 0)))).toBeTruthy()
  })

  test("containsPoint 2", () => {
    const a = PlaneFace.forVertices(P3.XY, [
      V(0, 0, 0),
      V(10, 0, 0),
      V(10, 10, 0),
      V(0, 10, 0),
    ])
    expect(a.containsPoint(V(-0.00000001, 11, 0))).toBeFalsy()

    expect(
      PlaneFace.forVertices(new P3(V(0, 0, -1), 0), [
        V(-1, -10, 0),
        V(0, 25, 0),
        V(25, 0, 0),
      ]).containsPoint(V(0, 0, 0)),
    ).toBeTruthy()
  })

  test("withHole", () => {
    const a = PlaneFace.forVertices(P3.XY, [
      V(0, 0, 0),
      V(10, 0, 0),
      V(10, 10, 0),
      V(0, 10, 0),
    ])
    const holeVertices = [V(2, 3, 0), V(8, 7, 0), V(7, 2, 0)]

    expect(a.containsPoint(V(-0.00000001, 11, 0))).toBeFalsy()
  })

  test("containsPoint2 3", () => {
    const face = new PlaneFace(new P3(V(0, 0, -1), 0), [
      new StraightEdge(
        new L3(V(0, 0.3333333333333333, 0), V(-1, 0, 0)),
        V(0.3333333333333333, 0.3333333333333333, 0),
        V(0, 0.3333333333333333, 0),
        -0.3333333333333333,
        0,
      ),
      new StraightEdge(
        new L3(V(0, 0, 0), V(0, 1, 0)),
        V(0, 0.3333333333333333, 0),
        V(0, 1, 0),
        0.3333333333333333,
        1,
      ),
      new StraightEdge(
        new L3(V(0, 1, 0), V(1, 0, 0)),
        V(0, 1, 0),
        V(1, 1, 0),
        0,
        1,
      ),
      new StraightEdge(
        new L3(V(1, 1, 0), V(0, -1, 0)),
        V(1, 1, 0),
        V(1, 0, 0),
        0,
        1,
      ),
      new StraightEdge(
        new L3(V(1, 0, 0), V(-1, 0, 0)),
        V(1, 0, 0),
        V(0.33333333333333326, 0, 0),
        0,
        0.6666666666666667,
      ),
      new StraightEdge(
        new L3(V(0.3333333333333333, 0, 0), V(0, 1, 0)),
        V(0.33333333333333326, 0, 0),
        V(0.3333333333333333, 0.3333333333333333, 0),
        0,
        0.3333333333333333,
      ),
    ])
    expect(
      face.containsPoint2(V(0.33333345254262287, 0.3333333333333333, 0)),
    ).toBe(PointVsFace.INSIDE)
  })

  test("pointsToInside", () => {
    const face = PlaneFace.forVertices(P3.XY, [
      V(0, 0, 0),
      V(10, 0, 0),
      V(10, 10, 0),
      V(0, 10, 0),
    ])

    expect(face.pointsToInside(V3.O, V3.X.negated())).toBe(PointVsFace.OUTSIDE)

    const v = V(10, 10, 0)
    expect(face.pointsToInside(v, V(-1, -1, 0))).toBe(PointVsFace.INSIDE)
    expect(face.pointsToInside(v, V(1, 1, 0))).toBe(PointVsFace.OUTSIDE)
    expect(face.pointsToInside(v, V(-1, 0, 0))).toBe(PointVsFace.ON_EDGE)
    expect(face.pointsToInside(v, V(0, -1, 0))).toBe(PointVsFace.ON_EDGE)

    const face2 = new PlaneFace(
      new PlaneSurface(new P3(V(0, 0, -1), 0)),
      [
        StraightEdge.throughPoints(V(0, 0, 0), V(0, 5, 0)),
        StraightEdge.throughPoints(V(0, 5, 0), V(5, 5, 0)),
        StraightEdge.throughPoints(V(5, 5, 0), V(5, 0, 0)),
        StraightEdge.throughPoints(V(5, 0, 0), V(0, 0, 0)),
      ],
      [],
    )
  })
  test("pointsToInside2", () => {
    const face = PlaneFace.forVertices(P3.XY, [
      V(0, 0, 0),
      V(10, 0, 0),
      V(10, 10, 0),
      V(0, 10, 0),
    ])

    expect(face.pointsToInside2(V3.O, V3.X.negated())).toBe(PointVsFace.OUTSIDE)

    const v = V(10, 10, 0)
    expect(face.pointsToInside2(v, V(-1, -1, 0))).toBe(PointVsFace.INSIDE)
    expect(face.pointsToInside2(v, V(1, 1, 0))).toBe(PointVsFace.OUTSIDE)
    expect(face.pointsToInside2(v, V(-1, 0, 0))).toBe(PointVsFace.ON_EDGE)
    expect(face.pointsToInside2(v, V(0, -1, 0))).toBe(PointVsFace.ON_EDGE)
  })
})

describe("BREP", () => {
  describe("assembleFacesFromLoops", () => {
    function ccwSquare(x, y, width, height) {
      return StraightEdge.chain(
        [V(x, y), V(x + width, y), V(x + width, y + height), V(x, y + height)],
        true,
      )
    }

    function reversedLoop(loop) {
      return loop.map((edge) => edge.flipped()).reverse()
    }

    const surface = new PlaneSurface(P3.XY)
    // 1(3(2), 4)
    const CCW1 = ccwSquare(0, 0, 10, 10) // outer loop
    const CCW2 = ccwSquare(5, 5, 1, 1),
      CW2 = reversedLoop(CCW2) // small in middle
    const CCW3 = ccwSquare(4, 4, 3, 3),
      CW3 = reversedLoop(CCW3) // around small in middle
    const CCW4 = ccwSquare(1, 1, 1, 1),
      CW4 = reversedLoop(CCW4)
    const CCW5 = ccwSquare(10, 10, 1, 1)
    const originalFace = new PlaneFace(surface, CCW1)

    function doLoopTest(loops, expected) {
      const actual = BRep.assembleFacesFromLoops(loops, surface, originalFace)
      expected = expected.length ? expected : [expected]
      expect(actual.length).toBe(expected.length)
      for (let i = 0; i < expected.length; i++) {
        expect(expected[i].likeFace(actual[i])).toBeTruthy()
      }
    }

    test("CCW: CCW", () => doLoopTest([CCW1], new PlaneFace(surface, CCW1, [])))

    test("CCW1(CW2): CCW1(CW2)", () =>
      doLoopTest([CCW1, CW2], new PlaneFace(surface, CCW1, [CW2])))

    test("CCW1(CCW4, CW2): CCW4", () =>
      doLoopTest([CCW1, CCW4, CW2], new PlaneFace(surface, CCW4, [])))

    test("CCW1(CCW3(CW2)): CCW3(CW2)", () =>
      doLoopTest([CCW1, CCW3, CW2], new PlaneFace(surface, CCW3, [CW2])))

    test("CCW1, CCW5: CCW1, CCW5", () =>
      doLoopTest(
        [CCW1, CCW5],
        [new PlaneFace(surface, CCW1, []), new PlaneFace(surface, CCW5, [])],
      ))

    test("CCW1(CW4, CW2): CCW1(CW4, CW2)", () =>
      doLoopTest([CCW1, CW4, CW2], new PlaneFace(surface, CCW1, [CW2, CW4])))
  })

  test("BREP.isCCW", () => {
    const vertices = [V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)]
    expect(isCCW(vertices, V(0, 0, 1))).toBeTruthy()
    expect(isCCW(vertices, V(0, 0, -1))).toBeFalsy()
  })

  test("edgePathFromSVG", () => {
    const path = `
      M80 80
      A 45 45, 0, 0, 0, 125 125
      L 125 80 Z
      M230 80
      A 45 45, 0, 1, 0, 275 125
      L 275 80 Z
      M80 230
      A 45 45, 0, 0, 1, 125 275
      L 125 230 Z
      M230 230
      A 45 45, 0, 1, 1, 275 275
      L 275 230 Z`
    const edges = edgePathFromSVG(path)
    outputLink({ edges })
    expect(edges.length).toBe(14)
  })

  test("rotateEdges", () => {
    const chain = edgePathFromSVG("m 2 0, 1 1, -2 1 Z").map((e) =>
      e.transform(M4.rotateX(90 * DEG)),
    )
    outputLink({ edges: chain })
    const a = B2T.rotateEdges(chain, TAU / 3, "rot").flipped()
    a.assertSanity()
    outputLink({ a })
    const rotatedEdges = new BRep(
      [
        new RotationFace(
          new ConicSurface(V(0, 0, -2), V(-3, 0, 0), V(0, 3, 0), V(0, 0, 3)),
          [
            new PCurveEdge(
              new EllipseCurve(
                V3.Z,
                V(3, 0, 0),
                V(0, 3, 0),
                0,
                3.141592653589793,
              ),
              V(3, 6.123233995736766e-17, 1),
              V(-1.4999999999999993, 2.598076211353316, 1),
              0,
              2.0943951023931953,
              undefined,
              V(0, 3, 0),
              V(-2.598076211353316, -1.4999999999999993, 0),
              "rotrib1",
            ),
            new StraightEdge(
              new L3(
                V(-0.9999999999999996, 1.7320508075688774, 0),
                V(-0.3535533905932737, 0.6123724356957946, 0.7071067811865476),
              ),
              V(-1.4999999999999993, 2.598076211353316, 1),
              V(-0.9999999999999996, 1.7320508075688774, 0),
              1.414213562373095,
              0,
            ),
            new PCurveEdge(
              new EllipseCurve(
                V3.O,
                V(2, 0, 0),
                V(0, 2, 0),
                0,
                3.141592653589793,
              ),
              V(-0.9999999999999996, 1.7320508075688774, 0),
              V(2, 0, 0),
              2.0943951023931953,
              0,
              undefined,
              V(1.7320508075688774, 0.9999999999999996, 0),
              V(0, -2, 0),
              "rotrib0",
            ),
            new StraightEdge(
              new L3(
                V(2, 0, 0),
                V(
                  0.7071067811865475,
                  4.329780281177466e-17,
                  0.7071067811865475,
                ),
              ),
              V(2, 0, 0),
              V(3, 6.123233995736766e-17, 1),
              0,
              1.4142135623730951,
            ),
          ],
          [],
        ),
        new RotationFace(
          new ConicSurface(V(0, 0, 2.5), V(3, 0, 0), V(0, 3, 0), V(0, 0, -1.5)),
          [
            new PCurveEdge(
              new EllipseCurve(V(0, 0, 2), V3.X, V3.Y, 0, 3.141592653589793),
              V(1, 1.2246467991473532e-16, 2),
              V(-0.4999999999999999, 0.8660254037844386, 2),
              0,
              2.0943951023931953,
              undefined,
              V3.Y,
              V(-0.8660254037844387, -0.4999999999999998, 0),
              "rotrib2",
            ),
            new StraightEdge(
              new L3(
                V(-1.4999999999999993, 2.598076211353316, 1),
                V(0.44721359549995776, -0.7745966692414835, 0.447213595499958),
              ),
              V(-0.4999999999999999, 0.8660254037844386, 2),
              V(-1.4999999999999993, 2.598076211353316, 1),
              2.2360679774997894,
              0,
            ),
            new PCurveEdge(
              new EllipseCurve(
                V3.Z,
                V(3, 0, 0),
                V(0, 3, 0),
                0,
                3.141592653589793,
              ),
              V(-1.4999999999999993, 2.598076211353316, 1),
              V(3, 6.123233995736766e-17, 1),
              2.0943951023931953,
              0,
              undefined,
              V(2.598076211353316, 1.4999999999999993, 0),
              V(0, -3, 0),
              "rotrib1",
            ),
            new StraightEdge(
              new L3(
                V(3, 6.123233995736766e-17, 1),
                V(
                  -0.8944271909999159,
                  2.738393491321013e-17,
                  0.4472135954999579,
                ),
              ),
              V(3, 6.123233995736766e-17, 1),
              V(1, 1.2246467991473532e-16, 2),
              0,
              2.23606797749979,
            ),
          ],
          [],
        ),
        new RotationFace(
          new ConicSurface(V(0, 0, 4), V(-2, 0, 0), V(0, 2, 0), V(0, 0, -4)),
          [
            new PCurveEdge(
              new EllipseCurve(
                V3.O,
                V(2, 0, 0),
                V(0, 2, 0),
                0,
                3.141592653589793,
              ),
              V(2, 0, 0),
              V(-0.9999999999999996, 1.7320508075688774, 0),
              0,
              2.0943951023931953,
              undefined,
              V(0, 2, 0),
              V(-1.7320508075688774, -0.9999999999999996, 0),
              "rotrib0",
            ),
            new StraightEdge(
              new L3(
                V(-0.4999999999999999, 0.8660254037844386, 2),
                V(-0.2236067977499788, 0.3872983346207417, -0.8944271909999159),
              ),
              V(-0.9999999999999996, 1.7320508075688774, 0),
              V(-0.4999999999999999, 0.8660254037844386, 2),
              2.23606797749979,
              0,
            ),
            new PCurveEdge(
              new EllipseCurve(V(0, 0, 2), V3.X, V3.Y, 0, 3.141592653589793),
              V(-0.4999999999999999, 0.8660254037844386, 2),
              V(1, 1.2246467991473532e-16, 2),
              2.0943951023931953,
              0,
              undefined,
              V(0.8660254037844387, 0.4999999999999998, 0),
              V(0, -1, 0),
              "rotrib2",
            ),
            new StraightEdge(
              new L3(
                V(1, 1.2246467991473532e-16, 2),
                V(
                  0.4472135954999579,
                  -5.476786982642026e-17,
                  -0.8944271909999159,
                ),
              ),
              V(1, 1.2246467991473532e-16, 2),
              V(2, 0, 0),
              0,
              2.23606797749979,
            ),
          ],
          [],
        ),
        new PlaneFace(
          new PlaneSurface(new P3(V3.Y, 0), V3.X, V(0, 0, -1)),
          [
            new StraightEdge(
              new L3(
                V(1, 1.2246467991473532e-16, 2),
                V(
                  0.4472135954999579,
                  -5.476786982642026e-17,
                  -0.8944271909999159,
                ),
              ),
              V(2, 0, 0),
              V(1, 1.2246467991473532e-16, 2),
              2.23606797749979,
              0,
            ),
            new StraightEdge(
              new L3(
                V(3, 6.123233995736766e-17, 1),
                V(
                  -0.8944271909999159,
                  2.738393491321013e-17,
                  0.4472135954999579,
                ),
              ),
              V(1, 1.2246467991473532e-16, 2),
              V(3, 6.123233995736766e-17, 1),
              2.23606797749979,
              0,
            ),
            new StraightEdge(
              new L3(
                V(2, 0, 0),
                V(
                  0.7071067811865475,
                  4.329780281177466e-17,
                  0.7071067811865475,
                ),
              ),
              V(3, 6.123233995736766e-17, 1),
              V(2, 0, 0),
              1.4142135623730951,
              0,
            ),
          ],
          [],
        ),
        new PlaneFace(
          new PlaneSurface(
            new P3(V(0.8660254037844387, 0.4999999999999998, 0), 0),
            V(0.4999999999999998, -0.8660254037844387, 0),
            V(0, 0, -1),
          ),
          [
            new StraightEdge(
              new L3(
                V(-0.9999999999999996, 1.7320508075688774, 0),
                V(-0.3535533905932737, 0.6123724356957946, 0.7071067811865476),
              ),
              V(-0.9999999999999996, 1.7320508075688774, 0),
              V(-1.4999999999999993, 2.598076211353316, 1),
              0,
              1.414213562373095,
            ),
            new StraightEdge(
              new L3(
                V(-1.4999999999999993, 2.598076211353316, 1),
                V(0.44721359549995776, -0.7745966692414835, 0.447213595499958),
              ),
              V(-1.4999999999999993, 2.598076211353316, 1),
              V(-0.4999999999999999, 0.8660254037844386, 2),
              0,
              2.2360679774997894,
            ),
            new StraightEdge(
              new L3(
                V(-0.4999999999999999, 0.8660254037844386, 2),
                V(-0.2236067977499788, 0.3872983346207417, -0.8944271909999159),
              ),
              V(-0.4999999999999999, 0.8660254037844386, 2),
              V(-0.9999999999999996, 1.7320508075688774, 0),
              0,
              2.23606797749979,
            ),
          ],
          [],
        ),
      ],
      true,
    )
    b2equals(a, rotatedEdges)
    const b = B2T.box(4, 4, 1.4)
      .translate(0, 0, -0.2)
      .rotateX(5 * DEG)
      .rotateY(-10 * DEG)
    testBRepAnd(a, b)
  })

  test("rotateEdges 2", () => {
    const chain = [
      new StraightEdge(
        new L3(
          V(198.46477746372744, 0, 244.94352900661897),
          V(-0.16064282877602398, 0, -0.9870126045612776),
        ),
        V(173.9128996557253, 0, 94.093266677553),
        V(198.46477746372744, 0, 244.94352900661897),
        152.83519342300417,
        0,
      ),
      new StraightEdge(
        new L3(
          V(131.35224103228387, 0, 180.2100595549249),
          V(0.7197488536413841, 0, 0.6942345336281635),
        ),
        V(198.46477746372744, 0, 244.94352900661897),
        V(131.35224103228387, 0, 180.2100595549249),
        93.24438113642698,
        0,
      ),
      new StraightEdge(
        new L3(
          V(173.9128996557253, 0, 94.093266677553),
          V(-0.44306356566594673, 0, 0.8964901989310186),
        ),
        V(131.35224103228387, 0, 180.2100595549249),
        V(173.9128996557253, 0, 94.093266677553),
        96.05993794472955,
        0,
      ),
    ]
    Edge.assertLoop(chain)
    const a = B2T.rotateEdges(chain, TAU / 3, "rot")
    expect(a).toMatchBRepSnapshot({ edges: chain })
  })

  test("like", () => {
    const a = B2T.tetrahedron(
      V(5, 5, 5),
      V(5, 5, -5),
      V(10, 12, 1),
      V(0, 12, 1),
    )
    const b = B2T.tetrahedron(
      V(5, 5, 5),
      V(10, 12, 1),
      V(5, 5, -5),
      V(0, 12, 1),
    )
    const c = B2T.tetrahedron(
      V(5, 5, 5),
      V(12, 12, 1),
      V(5, 5, -5),
      V(0, 12, 1),
    )

    expect(a.like(a)).toBeTruthy()

    expect(a.like(b)).toBeTruthy()
    expect(b.like(a)).toBeTruthy()

    expect(a.like(c)).toBeFalsy()
    expect(c.like(a)).toBeFalsy()

    expect(b.like(c)).toBeFalsy()
    expect(c.like(b)).toBeFalsy()
  })

  test("extrudeEdges", () => {
    const loop = StraightEdge.chain([V3.O, V3.X, V3.Y], true)
    const a = B2T.extrudeEdges(loop, P3.XY, V3.Z)
    expect(a.infiniteVolume).toBeTruthy()
    const b = B2T.extrudeEdges(loop, P3.XY, V3.Z.negated())
    expect(b.infiniteVolume).toBeFalsy()
  })

  test("splitsVolumeEnclosingFaces", () => {
    const brep = B2T.tetrahedron(
      V(0, 0, 0),
      V(10, 0, 0),
      V(0, 10, 0),
      V(0, 0, 10),
    )
    brep.buildAdjacencies()
    // pointing into tetrahedon
    const edge = (a, b) =>
      brep.faces
        .flatMap((face) => face.getAllEdges())
        .find((edge) => edge.a.like(a) && edge.b.like(b))
        .getCanon()
    expect(
      splitsVolumeEnclosingFaces(
        brep,
        edge(V(0, 0, 0), V(10, 0, 0)),
        V(0, 1, 1),
        V(0, -1, 1),
      ),
    ).toBe(INSIDE)
    expect(
      splitsVolumeEnclosingFaces(
        brep,
        edge(V(0, 0, 0), V(10, 0, 0)),
        V(0, 1, 1),
        V(0, 1, -1),
      ),
    ).toBe(INSIDE)

    // pointing out of tetrahedon
    expect(
      splitsVolumeEnclosingFaces(
        brep,
        edge(V(0, 0, 0), V(10, 0, 0)),
        V(0, -1, 0),
        V(0, 1, 1),
      ),
    ).toBe(OUTSIDE)
    expect(
      splitsVolumeEnclosingFaces(
        brep,
        edge(V(0, 0, 0), V(10, 0, 0)),
        V(0, -1, -1),
        V(0, -1, 1),
      ),
    ).toBe(OUTSIDE)
    expect(
      splitsVolumeEnclosingFaces(
        brep,
        edge(V(0, 0, 0), V(10, 0, 0)),
        V(0, -1, -1),
        V(0, 1, -1),
      ),
    ).toBe(OUTSIDE)
  })

  test("splitsVolumeEnclosingFacesP ", () => {
    const brep = B2T.box().flipped()
    brep.buildAdjacencies()
    const edge = (a, b) =>
      brep.faces
        .flatMap((face) => face.getAllEdges())
        .find((edge) => edge.a.like(a) && edge.b.like(b))
        .getCanon()
    expect(
      splitsVolumeEnclosingFacesP(
        brep,
        edge(V(0, 0, 0), V(1, 0, 0)),
        V(0.5, 0, 0),
        V(0, 0, -1),
        V3.O,
      ),
    ).toBe(INSIDE)
    expect(
      splitsVolumeEnclosingFacesP(
        brep,
        edge(V(0, 0, 1), V(1, 0, 1)),
        V(0.5, 0, 1),
        V(0, 0, 1),
        V3.O,
      ),
    ).toBe(INSIDE)

    const brep2 = B2T.box(12, 2, 3, "").translate(-6, 0, 1).flipped()
    brep2.buildAdjacencies()
    expect(
      splitsVolumeEnclosingFacesP(
        brep2,
        new StraightEdge(
          new L3(V(6, 0, 1), V(-1, 0, 0)),
          V(6, 0, 1),
          V(-6, 0, 1),
          0,
          12,
        ),
        V(-4.898979485566356, 5.999519546087386e-16, 1.000000000000001),
        V(0, 1.2246467991473547e-16, -4.898979485566356),
        V(-0.9797958971132712, 1.199903909217477e-16, 0.20000000000000023),
      ),
    ).toBe(INSIDE)
  })

  test("splitsVolumeEnclosingFaces 2", () => {
    const brep = B2T.tetrahedron(
      V(0, 0, 0),
      V(10, 0, 0),
      V(0, 10, 0),
      V(0, 0, 10),
    )
    brep.buildAdjacencies()
    // pointing out of tetrahedon
    const edge = (a, b) =>
      brep.faces
        .flatMap((face) => face.getAllEdges())
        .find((edge) => edge.a.like(a) && edge.b.like(b))
        .getCanon()
    expect(
      splitsVolumeEnclosingFaces(
        brep,
        edge(V(0, 0, 0), V(10, 0, 0)),
        V(0, 1, 0),
        V(0, 0, 1),
      ),
    ).toBe(COPLANAR_OPPOSITE)
    expect(
      splitsVolumeEnclosingFaces(
        brep,
        edge(V(0, 0, 0), V(10, 0, 0)),
        V(0, 1, 0),
        V(0, 0, -1),
      ),
    ).toBe(COPLANAR_SAME)
    expect(
      splitsVolumeEnclosingFaces(
        brep,
        edge(V(0, 0, 0), V(10, 0, 0)),
        V(0, 0, 1),
        V(0, 1, 0),
      ),
    ).toBe(COPLANAR_OPPOSITE)
  })

  test("splitsVolumeEnclosingFaces 3", () => {
    const brep = B2T.box(5, 5, 5).flipped()
    brep.buildAdjacencies()

    const edge = (a, b) =>
      brep.faces
        .flatMap((face) => face.getAllEdges())
        .find((edge) => edge.a.like(a) && edge.b.like(b))
        .getCanon()
    expect(
      splitsVolumeEnclosingFaces(
        brep,
        edge(V(0, 5, 0), V(0, 0, 0)),
        V(0, 0, -1),
        V(1, 0, 0),
      ),
    ).toBe(INSIDE)
    expect(
      splitsVolumeEnclosingFaces(
        brep,
        edge(V(0, 0, 0), V(0, 5, 0)),
        V(0, 0, -1),
        V(1, 0, 0),
      ),
    ).toBe(INSIDE)

    expect(
      splitsVolumeEnclosingFaces(
        brep,
        edge(V(0, 0, 0), V(0, 5, 0)),
        V(0, 0, 1),
        V(-1, 0, 0),
      ),
    ).toBe(COPLANAR_OPPOSITE)
    expect(
      splitsVolumeEnclosingFaces(
        brep,
        edge(V(0, 0, 0), V(0, 5, 0)),
        V(0, 0, 1),
        V(1, 0, 0),
      ),
    ).toBe(COPLANAR_SAME)
  })
  test("splitsVolumeEnclosingCone", () => {
    const brep = B2T.box(5, 5, 5)

    expect(splitsVolumeEnclosingCone(brep, V3.O, V(1, 1, 1))).toBe(INSIDE)
    expect(splitsVolumeEnclosingCone(brep, V3.O, V(-1, 1, 1))).toBe(OUTSIDE)
    expect(splitsVolumeEnclosingCone(brep, V3.O, V3.X)).toBe(
      ALONG_EDGE_OR_PLANE,
    )
    expect(splitsVolumeEnclosingCone(brep, V3.O, V3.X.negated())).toBe(OUTSIDE)
    expect(splitsVolumeEnclosingCone(brep, V(0, 5, 5), V3.Y)).toBe(OUTSIDE)
    expect(splitsVolumeEnclosingCone(brep, V3.O, V(1, 1, 0))).toBe(
      ALONG_EDGE_OR_PLANE,
    )
  })
  test("splitsVolumeEnclosingCone2", () => {
    function tst(b2, p, curve, t0, dir, result) {
      outputLink({ a: b2, edges: [edgeForCurveAndTs(curve, t0, t0 + dir)] })
      expect(splitsVolumeEnclosingCone2(b2, p, curve, t0, dir)).toBe(result)
    }

    const brep = B2T.box(5, 5, 5)

    tst(brep, V3.O, new L3(V3.O, V(1, 1, 1).unit()), 0, 1, INSIDE)
    tst(brep, V3.O, new L3(V3.O, V(-1, 1, 1).unit()), 0, 1, OUTSIDE)
    tst(brep, V3.O, new L3(V3.O, V3.X), 0, -1, OUTSIDE)
    tst(brep, V(0, 5, 5), new L3(V(0, 5, 5), V3.Y), 0, 1, OUTSIDE)
    tst(brep, V3.O, new L3(V3.O, V3.X), 0, 1, ALONG_EDGE_OR_PLANE)
    tst(brep, V3.O, new L3(V3.O, V(1, 1, 0).unit()), 0, 1, ALONG_EDGE_OR_PLANE)

    tst(
      brep,
      V(5, 0, 0),
      EllipseCurve.UNIT.translate(4),
      0,
      1,
      ALONG_EDGE_OR_PLANE,
    )
    tst(
      brep,
      V(5, 0, 0),
      EllipseCurve.UNIT.translate(4).rotateX(5 * DEG),
      0,
      1,
      INSIDE,
    )
    tst(brep, V(0, 0, 0), EllipseCurve.UNIT.translate(-1), 0, 1, OUTSIDE)

    const semicyl = B2T.cylinder(1, 1, TAU / 2)
    tst(semicyl, V(-1, 0, 0), L3.X, -1, -1, OUTSIDE)
    tst(semicyl, V(-1, 0, 0), L3.X, -1, 1, ALONG_EDGE_OR_PLANE)
    tst(semicyl.flipped(), V(-1, 0, 0), L3.X, -1, -1, INSIDE)
    tst(semicyl.flipped(), V(-1, 0, 0), L3.X, -1, 1, ALONG_EDGE_OR_PLANE)
  })

  test("splitsVolumeEnclosingCone 2", () => {
    const a = B2T.box(1, 1, 1, "box").minus(B2T.box(1 / 3, 1 / 3, 1, "cut"))
    expect(splitsVolumeEnclosingCone(a, V(1 / 3, 1 / 3, 0), V3.X)).toBe(
      ALONG_EDGE_OR_PLANE,
    )
  })
  test("planeFaceEdgeISPsWithPlane", () => {
    const brep = B2T.extrudeVertices(
      [
        V(3.318793683290006e-10, 9.12934465653568e-10, 3.8212177478055603e-10), // 0 0 0
        V(9.999999999857609, 4.561232401125957e-13, 5.033029171552123e-10), // 10 0 0
        V(2.984222284459465e-10, 9.999999999852161, -1.0461618780772852e-9),
      ], // 0 10 0
      P3.XY,
      V(-1.8192047548394726e-10, -3.0150747681012967e-10, -5.0000000009123795),
      "ex0",
    )
    /*var result = planeFaceEdgeISPsWithPlane(brep,
         new BREP.Face([
         V(9.999999999675689, -3.010513535700171e-10, -5.000000000409076), // 10 0 -5
         V(1.4995889284505332e-10, 6.114269888434384e-10, -5.000000000530258), // 0 0 -5
         V(3.318793683290006e-10, 9.12934465653568e-10, 3.8212177478055603e-10), // 0 0 0
         V(9.999999999857609, 4.561232401125957e-13, 5.033029171552123e-10)], // 10 0 0
         new P3(V(9.12478342464039e-11, 1, -6.03014953543423e-11), 9.129344656608087e-10)), // 0 1 0
         L3(V(-1.3833878355530264e-10, 6.114269894465992e-10, -4.999999990964091), V(-1, 9.12478342480723e-11, 2.7667756772219476e-11)),
         new P3(V(2.766775686256173e-11, 9.90075577448337e-10, 1), -4.999999990964091),
         true, true)
         assert.deepEqual(result, [])*/
    const face = brep.faces[2] as PlaneFace
    let result = face
      .edgeISPsWithPlane(L3.Y.translate(0, 0, -1), P3.XY.translate(0, 0, -1))
      .map((is) => is.p)
    expect(result).toBeLike([V(0, 0, -1), V(0, 10, -1)])
    result = face.edgeISPsWithPlane(L3.Y, P3.XY).map((is) => is.p)
    expect(result).toBeLike([V(0, 0, 0), V(0, 10, 0)])
    result = face
      .translate(0, 0, 10)
      .edgeISPsWithPlane(L3.Y.translate(0, 0, 6), P3.XY.translate(0, 0, 6))
      .map((is) => is.p)
    expect(result).toBeLike([V(0, 0, 6), V(0, 10, 6)])
  })

  describe("intersection", () => {
    test("box(5, 5, 5) - box(1, 1, 6)", () => {
      const a = B2T.box(5, 5, 5, "a")
      const b = B2T.box(1, 1, 6, "b").flipped()
      testBRepAnd(a, b)
    })

    test("sphere() - sphere(0.2).translate(0, 1.05)", () => {
      const a = B2T.sphere()
      const v = V(-0.5773502691896257, 0.5773502691896257, 0.5773502691896257)
      // const m4 = M4.forSys(v, v.getPerpendicular
      const b = B2T.sphere(0.2).translate(0, 1.05).rotateAB(V3.Y, v).flipped()
      testBRepAnd(a, b)
    })

    test.skip("sphere() - cone()", () => {
      const a = B2T.sphere()
      const b = B2T.cone()
        .flipped()
        //.scale(0.5, 0.6, 1)
        .translate(0.4, 0.6)
      const result = BRep.EMPTY
      testBRepAnd(a, b, result)
    })

    test("box(5, 5, 5) - box(1, 1, 5)", () => {
      const a = B2T.box(5, 5, 5)
      const b = B2T.box(1, 1, 5).flipped()
      testBRepAnd(a, b)
    })

    test("edge/face intersection", () => {
      const wideBox = B2T.box(10, 10, 5)
      const extrusion = B2T.extrudeVertices(
        [V(1, 0), V(0, 3), V(-2, 5), V(5, 5)],
        P3.XY.flipped(),
        V(0, 0, 10),
        "lol",
      ).translate(0, 1, 1)
      testBRepOp(wideBox, extrusion, () => wideBox.plus(extrusion))
    })

    test("edge/edge intersection", () => {
      const wideBox = B2T.box(10, 10, 10)
      const t = B2T.tetrahedron(
        V(0, 0, 5),
        V(0, 0, 15),
        V(8, 2, 13),
        V(-1, 7, 12),
      ).flipped()
      testBRepAnd(wideBox, t)
    })
    test("edge/edge intersection 2", () => {
      const box = B2T.box(10, 10, 10)
      const tetra = B2T.tetrahedron(
        V(0, 0, 5),
        V(0, 0, 15),
        V(8, 2, 13),
        V(1, 7, 12),
      ).flipped()
      testBRepAnd(box, tetra)
    })
    test.skip("edge/point intersection", () => {
      const wideBox = B2T.box(10, 10, 10)
      const t = B2T.tetrahedron(
        V(0, 0, 5),
        V(1, 1, 15),
        V(8, 2, 13),
        V(1, 7, 12),
      ).flipped()
      const result = BRep.EMPTY
      testBRepAnd(wideBox, t, result)
    })

    test("box(10, 10, 5) + box(10, 10, 5).translate(0, 0, 5) touching coplanar faces result contains both", () => {
      const a = B2T.box(10, 10, 5),
        b = a.translate(0, 0, 5)

      // the expected result is the union of a and b's faces, minus the ones where they touch, which we remove
      // with an AABB test
      const testAABB = AABB.forXYZ(12, 12, 2).translate(-1, -1, 4)
      const result = new BRep(
        a.faces
          .concat(b.faces)
          .filter((face) => !testAABB.containsAABB(face.getAABB())),
        false,
      )
      testBRepOp(a, b, () => b.plus(a), result)
    })

    test("box(10, 10, 5) - box(10, 5, 5).translate(0, 0, 5) - box(5, 5, 5).translate(5, 5, 5)", () => {
      const box = B2T.box(10, 10, 10)
      const box2 = B2T.box(10, 5, 5).translate(0, 0, 5)
      const box3 = B2T.box(5, 5, 5).translate(5, 5, 5).flipped()
      testBRepAnd(box.minus(box2), box3)
    })

    test("box(10, 10, 5) && box(10, 10, 5).translate(3, 3, 0) overlapping faces result contains intersection of faces", () => {
      const box = B2T.box(10, 10, 5, "a")
      const box2 = B2T.box(10, 10, 4, "b").translate(3, 3, 0)
      testBRepAnd(box, box2, B2T.box(7, 7, 4).translate(3, 3, 0))
    })

    test("box(10, 10, 5) + box(10, 10, 5).translate(3, 3, 0) overlapping faces result contains union of faces", () => {
      const boxA = B2T.box(10, 10, 5, "boxA").flipped(),
        boxB = boxA.translate(3, 3, 0)
      //const result = B2T.extrudeVertices([V3.O, V(0, 10), V(3, 10), V(3, 13), V(13, 13), V(13, 3), V(10, 3),
      // V(10, 0)], P3.XY.flipped(), V(0, 0, 5), 'result').flipped()
      testBRepAnd(boxA, boxB)
    })

    test("box(10,10,5) + box(4,10,2).translate(2, 0, 3) overlapping faces result contains union of faces 2", () => {
      const box = B2T.box(10, 10, 5),
        box2 = B2T.box(4, 10, 2).translate(2, 0, 3)
      testBRepOp(box, box2, () => box.plus(box2))
    })

    test("box(10,10,10) + box(10,12,12).translate(0, 0, -2) overlapping faces result contains union of faces 3", () => {
      const box = B2T.box(10, 10, 10, "box0"),
        box2 = B2T.box(10, 12, 12, "box").translate(0, 0, -2)
      testBRepOp(box, box2, () => box.plus(box2))
    })

    test("box(10, 10, 10) - box().scale(8 / 3 ** 0.5).rotateAB(V3.ONES, V3.X)", () => {
      const box = B2T.box(10, 10, 10, "box0")
      const box2 = B2T.box(1, 1, 1)
        .scale(8 / 3 ** 0.5)
        .rotateAB(V3.XYZ, V3.X)
        .flipped()
      testBRepAnd(box, box2)
    })

    test("box(10, 10, 10) - tetrahedron(V(5,0,10),V(5,3,10),V(7,1,9),V(7,1,11))", () => {
      const box = B2T.box(10, 10, 10, "box0")
      const box2 = B2T.tetrahedron(
        V(5, 0, 10),
        V(5, 3, 10),
        V(7, 1, 9),
        V(7, 1, 11),
      ).flipped()
      testBRepAnd(box, box2)
    })

    test("box(10, 10, 10) - box().rotateAB(V3.ONES, V3.Y).translate(5,0,10)", () => {
      const box = B2T.box(10, 10, 10, "box0"),
        box2 = B2T.box().rotateAB(V3.XYZ, V3.Y).translate(5, 0, 10).flipped()
      testBRepAnd(box, box2)
    })

    test("menger(1) - box(2,2,1).rotateZ(-45*DEG)", () => {
      const box = B2T.box(1, 1, 1, "box0").minus(
        B2T.box(1 / 3, 1 / 3, 1).translate(1 / 3, 1 / 3, 0),
      )
      const box2 = B2T.box(2, 2, 1)
        .rotateZ(-45 * DEG)
        .flipped()
      testBRepAnd(box, box2)
    })

    test("sphere(5) - box(12,2,3).translate(-6, 0,1)", () => {
      const a = B2T.sphere(5)
      const b = B2T.box(12, 2, 3).translate(-6, 0, 1).flipped()
      testBRepAnd(a, b)
    })

    test("sphere(5) AND octahedron().scale(6)", () => {
      const a = B2T.sphere(5)
      const b = B2T.octahedron().scale(6)
      testBRepAnd(a, b)
    })

    test("sphere(5) - box(12,2,3).translate(-6, 0,1) 2", () => {
      const a = B2T.sphere(5)
      const b = B2T.box(12, 2, 2)
        .translate(-6, -1, 1)
        .rotateX(-10 * DEG)
        .flipped()
      testBRepAnd(a, b)
    })

    test("sphere(5) AND tetrahedron(V3.ZERO, V3.Y, V3.Z, V(7,0,0))", () => {
      const a = B2T.sphere(5)
      const b = B2T.tetrahedron(V3.O, V3.Y, V3.Z, V(7, 0, 0))
      testBRepAnd(a, b)
    })

    test("sphere(5) AND tetrahedron(V3.ZERO, V3.X, V3.Y, V3.Z).scale(7)", () => {
      const a = B2T.sphere(5)
      const b = B2T.tetrahedron(V3.O, V3.X, V3.Y, V3.Z).scale(7)
      testBRepAnd(a, b)
    })

    test("sphere(5) AND tetrahedron(V3.ZERO, V3.X, V3.Y, V3.Z.negated()).scale(7)", () => {
      const a = B2T.sphere(5)
      const b = B2T.tetrahedron(
        V(-6, 0, 0),
        V(6, 0, 0),
        V(0, -4, 0),
        V(0, 0, -6),
      )
      testBRepAnd(a, b)
    })

    test("box - doublebarrel", () => {
      const a = B2T.extrudeEdges(
        [
          new StraightEdge(
            new L3(
              V(-409.7544757659938, 183.97231686845578, 0),
              V(0.9961411421023264, -0.08776573939227537, 0),
            ),
            V(-409.7544757659938, 183.97231686845578, 0),
            V(-42.014475765994064, 151.57231686845583, 0),
            0,
            369.16455355301895,
          ),
          new StraightEdge(
            new L3(
              V(-42.014475765994064, 151.57231686845583, 0),
              V(-0.35632458233389064, 0.9343622381199803, 0),
            ),
            V(-42.014475765994064, 151.57231686845583, 0),
            V(-114.91447576599401, 342.73231686845577, 0),
            0,
            204.5887474911559,
          ),
          new StraightEdge(
            new L3(
              V(-114.91447576599401, 342.73231686845577, 0),
              V(-0.9951217646853636, 0.09865431287829038, 0),
            ),
            V(-114.91447576599401, 342.73231686845577, 0),
            V(-490.7544757659938, 379.99231686845576, 0),
            0,
            377.68242373719204,
          ),
          new StraightEdge(
            new L3(
              V(-490.7544757659938, 379.99231686845576, 0),
              V(0.3819019948313259, -0.9242028274918087, 0),
            ),
            V(-490.7544757659938, 379.99231686845576, 0),
            V(-409.7544757659938, 183.97231686845578, 0),
            0,
            212.09629982628172,
          ),
        ],
        new P3(V3.Z, 0),
        V(0, 0, -100),
        "box",
      ).minus(
        B2T.extrudeEdges(
          [
            new PCurveEdge(
              new EllipseCurve(
                V(-379.9148578120817, 246.4368702165194, 0),
                V(21.599999999999966, -23.39999999999995, 0),
                V(23.39999999999995, 21.599999999999966, 0),
                0,
                3.141592653589793,
              ),
              V(-358.31485781208175, 223.03687021651945, 0),
              V(-401.5148578120817, 269.83687021651934, 0),
              0,
              3.141592653589793,
              null,
              V(23.399999999999945, 21.59999999999997, 0),
              V(-23.39999999999995, -21.599999999999966, 0),
            ),
            new PCurveEdge(
              new EllipseCurve(
                V(-379.9148578120817, 246.4368702165194, 0),
                V(-21.599999999999966, 23.39999999999995, 0),
                V(-23.39999999999995, -21.599999999999966, 0),
                0,
                3.141592653589793,
              ),
              V(-401.5148578120817, 269.83687021651934, 0),
              V(-358.31485781208175, 223.03687021651945, 0),
              0,
              3.141592653589793,
              null,
              V(-23.39999999999995, -21.599999999999966, 0),
              V(23.399999999999952, 21.599999999999962, 0),
            ),
          ],
          new P3(V3.Z, 0),
          V(0, 0, -50),
          "cyl1",
        ),
      )
      const b = B2T.extrudeEdges(
        [
          new PCurveEdge(
            new EllipseCurve(
              V(-333.7047164374913, 265.5319699580857, 0),
              V(21.599999999999966, -23.39999999999995, 0),
              V(23.39999999999995, 21.599999999999966, 0),
              0,
              3.141592653589793,
            ),
            V(-312.10471643749133, 242.13196995808573, 0),
            V(-355.30471643749127, 288.9319699580857, 0),
            0,
            3.141592653589793,
            null,
            V(23.399999999999945, 21.59999999999997, 0),
            V(-23.39999999999995, -21.599999999999966, 0),
          ),
          new PCurveEdge(
            new EllipseCurve(
              V(-333.7047164374913, 265.5319699580857, 0),
              V(-21.599999999999966, 23.39999999999995, 0),
              V(-23.39999999999995, -21.599999999999966, 0),
              0,
              3.141592653589793,
            ),
            V(-355.30471643749127, 288.9319699580857, 0),
            V(-312.10471643749133, 242.13196995808573, 0),
            0,
            3.141592653589793,
            null,
            V(-23.39999999999995, -21.599999999999966, 0),
            V(23.399999999999952, 21.599999999999962, 0),
          ),
        ],
        new P3(V3.Z, 0),
        V(0, 0, -50),
        "cyl2",
      ).flipped()
      testBRepAnd(a, b)
    })

    test("sphere() - cubelet near -V3.Y", () => {
      const a = B2T.sphere()
      const b = B2T.box(0.2, 0.2, 0.2, "")
        .translate(0, 0.95, 0)
        .rotateAB(V3.Y, V(0, -0.9341723589627158, 0.35682208977308993))
        .flipped()
      testBRepAnd(a, b)
    })

    test("sphere() - BRep w/ PCS", () => {
      const a = B2T.sphere()
      const b = B2T.extrudeEdges(
        [
          edgeForCurveAndTs(BezierCurve.QUARTER_CIRCLE, 0, 1),
          StraightEdge.throughPoints(V3.Y, V3.X),
        ],
        P3.XY,
        V3.Z.negated(),
      )
        .scale(0.2, 0.2, 2)
        .rotateX(85 * DEG)
        .translate(0.1, 0.1, 0.4)
        .flipped()
      testBRepAnd(a, b)
    })

    test("sphere() - BRep w/ PCS 2", () => {
      const a = B2T.sphere()
      const b = B2T.extrudeEdges(
        [
          edgeForCurveAndTs(BezierCurve.QUARTER_CIRCLE, 0, 1),
          StraightEdge.throughPoints(V3.Y, V3.X),
        ],
        P3.XY,
        V3.Z.negated(),
      )
        .scale(0.2, 0.2, 2)
        .translate(0.1, -0.1, 1.2)
        .flipped()
      testBRepAnd(a, b)
    })

    test("box() - BRep w/ PCS 2", () => {
      const a = B2T.box()
      const b = B2T.extrudeEdges(
        [
          edgeForCurveAndTs(BezierCurve.QUARTER_CIRCLE, 0, 1),
          StraightEdge.throughPoints(V3.Y, V3.X),
        ],
        P3.XY,
        V3.Z.negated(),
      )
        .scale(0.2, 0.2, 4)
        .translate(-0.1, 0.4, 2)
        .rotateX(10 * DEG)
        .flipped()
      testBRepAnd(a, b)
    })

    test("BRep w/ PCS - sphere()", () => {
      const a = B2T.sphere().flipped()
      const b = B2T.extrudeEdges(
        [
          edgeForCurveAndTs(BezierCurve.QUARTER_CIRCLE, 0, 1),
          StraightEdge.throughPoints(V3.Y, V3.X),
        ],
        P3.XY,
        V3.Z.negated(),
      )
        .scale(0.2, 0.2, 2)
        .translate(0.1, -0.1, 1.2)
      testBRepAnd(a, b)
    })

    test.skip("BRep w/ PCS and cylinder", () => {
      const bezierEdge = edgeForCurveAndTs(BezierCurve.EX2D, 0, 1)
      const a = B2T.extrudeEdges(
        [
          bezierEdge,
          ...StraightEdge.chain(
            [bezierEdge.b, V3.X.negated(), bezierEdge.a],
            false,
          ),
        ],
        P3.XY.flipped(),
      )
      const b = B2T.cylinder(0.2, 2)
        .rotateY(90 * DEG)
        .translate(0, 0, 0.5)
      const result = BRep.EMPTY
      testBRepAnd(a, b, result)
    })

    test("sphere() - BRep w/ PCS - sphere(0.9)", () => {
      const a = B2T.sphere(0.9).flipped()
      const b = B2T.extrudeEdges(
        [
          edgeForCurveAndTs(BezierCurve.QUARTER_CIRCLE, 0, 1),
          StraightEdge.throughPoints(V3.Y, V3.X),
        ],
        P3.XY,
        V3.Z.negated(),
      )
        .scale(0.2, 0.2, 2)
        .translate(0.1, -0.1, 1.2)
        .flipped()
      const c = B2T.sphere()
      const d = a.and(b)
      testBRepAnd(d, c)
    })

    test("sphere() - cylinder().scale(0.5,0.1,4).translate(0.5,0,-2)", () => {
      const a = B2T.sphere(1, undefined, PI)
      const b = B2T.cylinder()
        .scale(0.5, 0.1, 4)
        .translate(0.5, 0, -2)
        .flipped()
      testBRepAnd(a, b)
    })

    test("remove half of a half-pie", () => {
      const pie = B2T.puckman(8, 180 * DEG, 5, "pie/2")
      const boxKnife = B2T.box(11, 10, 7, "knife")
        .translate(-10, -1, -1)
        .flipped()

      const resultTopPoint = V(1, 8 * Math.sin(Math.acos(1 / 8)), 0)
      const result = B2T.extrudeEdges(
        [
          StraightEdge.throughPoints(V(8, 0, 0), V(1, 0, 0)),
          StraightEdge.throughPoints(V(1, 0, 0), resultTopPoint),
          edgeForCurveAndTs(EllipseCurve.semicircle(8), Math.acos(1 / 8), 0),
        ],
        P3.XY.flipped(),
        V(0, 0, 5),
        "pie/4",
      )
      testBRepAnd(pie, boxKnife, result)
    })

    test("box() - cylinder(0.2,2).translate(0.5,0.2)", () => {
      const a = B2T.box()
      const b = B2T.cylinder(0.2, 2).translate(0.5, 0.2).flipped()
      testBRepAnd(a, b)
    })

    //test('cylinder() - triangle(1.5)', () => {
    // file:///C:/Users/aval/Desktop/cs/viewer.html?a=B2T.cylinder()&b=B2T.extrudeEdges(Edge.ngon(3,1.5))&c=a.and(b).translate(3)
    // const a = B2T.cylinder() const b = B2T.extrudeEdges(Edge.round(Edge.ngon(3,1.6),0.5)) const result =
    // BRep.EMPTY testBRepAnd(assert, a, b, result) },
    test("box - snug sphere", () => {
      // file:///C:/Users/aval/Desktop/cs/viewer.html?a=B2T.cylinder()&b=B2T.extrudeEdges(Edge.ngon(3,1.5))&c=a.and(b).translate(3)
      const a = B2T.box(2, 2, 3)
      const b = B2T.sphere().translate(1, 1).flipped()
      testBRepAnd(a, b)
    })
    test("pcs - triangle", () => {
      // file:///C:/Users/aval/Desktop/cs/viewer.html?a=B2T.cylinder()&b=B2T.extrudeEdges(Edge.ngon(3,1.5))&c=a.and(b).translate(3)
      const e1 = edgeForCurveAndTs(BezierCurve.EX3D, 0.5).projectXY()
      const edges = [e1, StraightEdge.throughPoints(e1.b, e1.a)]
      const a = B2T.extrudeEdges(edges, P3.XY.flipped(), V3.Z)
      const p = e1.curve.at(0.95).plus(V(-NLA_PRECISION * 2, 0, 0))
      const b = B2T.extrudeVertices(
        [p, V(0, -0.2), V(0, 0.2)],
        P3.XY.flipped(),
        V3.Z,
      ).flipped()
      testBRepAnd(a, b)
    })
    test("box() - cylinder(0.2,3).translate(0.5,0.2,-1).rotateX(10*DEG)", () => {
      const a = B2T.box()
      const b = B2T.cylinder(0.2, 1)
        .translate(0.5, 0.2, -0.5)
        .rotateX(-10 * DEG)
        .flipped()
      testBRepAnd(a, b)
    })
    test("star - ball", () => {
      const a = B2T.extrudeEdges(edgeStar(4, 2, 1), P3.XY, V3.Z.negated())
      const b = B2T.sphere().flipped()
      testBRepAnd(a, b)
    })
    test("cylinder - cube", () => {
      const a = B2T.cylinder(4, 4)
      const b = B2T.extrudeEdges(edgeNgon(4, 4), P3.XY, V(0, 0, 4)).translate(
        0,
        0,
        2,
      )
      testBRepAnd(a, b)
    })
    test.skip("fail", () => {
      const a = B2T.extrudeEdges(
        [
          new PCurveEdge(
            new BezierCurve(
              V(-165.04412089048037, 233.67721659018514, 0),
              V(162.0276742021058, 92.61628384754499, 0),
              V(-26.526878671628538, 103.9576526746689, 0),
              V(106.49691429238996, -38.442642927294784, 0),
              -0.1,
              1.1,
            ),
            V(106.49691429238996, -38.442642927294784, 0),
            V(-165.04412089048037, 233.67721659018514, 0),
            1,
            0,
            null,
            V(-399.07137889205546, 427.2008868058911, 0),
            V(-981.2153852777585, 423.18279822792044, 0),
          ),
          new PCurveEdge(
            new BezierCurve(
              V(-234.98404318638632, -41.414503714148154, 0),
              V(-263.95799566517644, 2.0493326379292967, 0),
              V(-96.30326822511165, 115.1793952369976, 0),
              V(-165.04412089048037, 233.67721659018514, 0),
              -0.1,
              1.1,
            ),
            V(-165.04412089048037, 233.67721659018514, 0),
            V(-234.98404318638632, -41.414503714148154, 0),
            1,
            0,
            null,
            V(206.22255799610616, -355.4934640595626, 0),
            V(86.92185743637037, -130.39150905623234, 0),
          ),
          new PCurveEdge(
            new BezierCurve(
              V(106.49691429238996, -38.442642927294784, 0),
              V(6.920773903871436, -168.34765938596584, 0),
              V(-96.36814809386107, 9.19324831183017, 0),
              V(-234.98404318638632, -41.414503714148154, 0),
              -0.1,
              1.1,
            ),
            V(-234.98404318638632, -41.414503714148154, 0),
            V(106.49691429238996, -38.442642927294784, 0),
            1,
            0,
            null,
            V(415.84768527757575, 151.82325607793496, 0),
            V(298.7284211655556, 389.7150493760132, 0),
          ),
        ],
        new P3(V3.Z, 0),
        V(0, 0, -100),
        "extrude42",
      )
      const b = B2T.extrudeEdges(
        [
          new StraightEdge(
            new L3(
              V(-113.15294177340922, 90.2593377922355, 0),
              V(0.9138115486202569, -0.40613846605344794, 0),
            ),
            V(-13.152941773409395, 45.814893347791084, 0),
            V(-113.15294177340922, 90.2593377922355, 0),
            109.43175335328989,
            0,
          ),
          new StraightEdge(
            new L3(
              V(-150.09473256378868, -35.10025916844387, 0),
              V(0.2826685651896956, 0.9592176407122622, 0),
            ),
            V(-113.15294177340922, 90.2593377922355, 0),
            V(-150.09473256378868, -35.10025916844387, 0),
            130.68941983551744,
            0,
          ),
          new StraightEdge(
            new L3(
              V(-13.152941773409395, 45.814893347791084, 0),
              V(-0.8609402861542653, -0.5087060287401869, 0),
            ),
            V(-150.09473256378868, -35.10025916844387, 0),
            V(-13.152941773409395, 45.814893347791084, 0),
            159.06073045098708,
            0,
          ),
        ],
        new P3(V3.Z, 0),
        V(0, 0, -100),
        "extrude90",
      ).flipped()
      testBRepAnd(a, b)
    })

    test("box - torus", () => {
      const a = B2T.box(5, 5, 3).translate(-2.5, -2.5, -1.5)
      const b = B2T.torus(1, 3)
        .rotateX(10 * DEG)
        .flipped()
      testBRepAnd(a, b)
    })

    test('sphere() - "a"', async () => {
      const a = B2T.sphere()
      const b = B2T.text(
        "a",
        64,
        64,
        await B2T.loadFont(FONT_PATH + "/FiraSansMedium.woff"),
      )
        .scale(0.5 / 32)
        .translate(-0.25, -0.25, 1.2)
        .flipped()
      testBRepAnd(a, b)
    })

    test("cylinder(1,2) AND cylinder(1,2).rotateZ(PI/2).translate(0,0,1)", () => {
      const a = B2T.cylinder(1, 2)
      const b = B2T.cylinder(1, 2)
        .rotateZ(PI / 2)
        .translate(0, 0, 1)
      testBRepAnd(a, b)
    })

    test("box - semicylinder", () => {
      const box = B2T.box(4, 2, 2)
      const cyl = B2T.cylinder(1, 2, 180 * DEG)
        .translate(2)
        .flipped()
      testBRepAnd(box, cyl)
    })

    //async 'B2T.sphere() - 'atest(' - B2T.sphere(0.9)', () => {
    //	const a = B2T.sphere(0.9).flipped()
    //	const b = B2T.text('a', await B2T.loadFont('fonts/FiraSansMedium.woff'), 64,
    // 64).scale(0.5/32).translate(-0.25,-0.25,1.2).flipped() const c = B2T.sphere() const d = a.and(b) const
    // result = new BRep([ new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.20296874999999998,
    // 0.11703124999999998, 1.2), V(0.20296874999999998, 0.2240625, 1.2), V(0.14500000000000002, 0.2890625, 1.2),
    // V(0.010937499999999989, 0.2890625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [ new
    // PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.20296874999999998, 0.11703124999999998,
    // 1.2), V(0.20296874999999998, 0.2240625, 1.2), V(0.14500000000000002, 0.2890625, 1.2),
    // V(0.010937499999999989, 0.2890625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new
    // EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.20296874999999998, 0.11703124999999998,
    // 0.9721663299286162), V(0.010937499999999989, 0.2890625, 0.9572477433702835), 1), V(0.20296874999999998,
    // 0.11703124999999998, 0.9721663299286162), V(0.010937499999999989, 0.2890625, 0.9572477433702835), 0, 100,
    // null, V(0, 0.003208541417521995, -0.0003862503783657504), V(-0.004021833637101672, 0, 0.00004595341771287279)), new StraightEdge(new L3(V(0.010937499999999989, 0.2890625, 1.2), V(0, 0, -1)), V(0.010937499999999989, 0.2890625, 0.9572477433702835), V(0.010937499999999989, 0.2890625, 0.852245998633904), 0.24275225662971645, 0.3477540013660959), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.20296874999999998, 0.11703124999999998, 1.2), V(0.20296874999999998, 0.2240625, 1.2), V(0.14500000000000002, 0.2890625, 1.2), V(0.010937499999999989, 0.2890625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(0.20296874999999995, 0.11703124999999996, 0.8689691438980299), V(0.010937499999999987, 0.2890625, 0.8522459986339039), -1), V(0.010937499999999989, 0.2890625, 0.852245998633904), V(0.20296874999999998, 0.11703124999999998, 0.8689691438980301), 100, 0, null, V(0.004021823173538629, 0, -0.0000516150161233843), V(0, -0.0032079393598033606, 0.0004320396826956172)), new StraightEdge(new L3(V(0.20296874999999998, 0.11703124999999998, 1.2), V(0, 0, -1)), V(0.20296874999999998, 0.11703124999999998, 0.8689691438980301), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), 0.3310308561019698, 0.22783367007138378)], []), new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.010937499999999989, 0.2890625, 1.2), V(-0.04093749999999999, 0.2890625, 1.2), V(-0.1, 0.2790625, 1.2), V(-0.16499999999999998, 0.255, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.010937499999999989, 0.2890625, 1.2), V(-0.04093749999999999, 0.2890625, 1.2), V(-0.1, 0.2790625, 1.2), V(-0.16499999999999998, 0.255, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.010937499999999989, 0.2890625, 0.9572477433702835), V(-0.16499999999999998, 0.255, 0.9527591510974849), 1), V(0.010937499999999989, 0.2890625, 0.9572477433702835), V(-0.16499999999999998, 0.255, 0.9527591510974849), 0, 100, null, V(-0.0015562475396635686, 0, 0.000017781663715540364), V(-0.0019497960757584164, -0.0007217995088144136, -0.0001444829762001187)), new StraightEdge(new L3(V(-0.16499999999999998, 0.255, 1.2), V(0, 0, -1)), V(-0.16499999999999998, 0.255, 0.9527591510974849), V(-0.16499999999999998, 0.255, 0.8472012747865765), 0.24724084890251508, 0.3527987252134235), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.010937499999999989, 0.2890625, 1.2), V(-0.04093749999999999, 0.2890625, 1.2), V(-0.1, 0.2790625, 1.2), V(-0.16499999999999998, 0.255, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(0.010937499999999987, 0.2890625, 0.8522459986339039), V(-0.16499999999999998, 0.25499999999999995, 0.8472012747865764), -1), V(-0.16499999999999998, 0.255, 0.8472012747865765), V(0.010937499999999989, 0.2890625, 0.852245998633904), 100, 0, null, V(0.001949741977734284, 0.0007217794821420186, 0.0001624804665392168), V(0.0015562468960623505, 6.688077170406378e-21, -0.00001997246153454078)), new StraightEdge(new L3(V(0.010937499999999989, 0.2890625, 1.2), V(0, 0, -1)), V(0.010937499999999989, 0.2890625, 0.852245998633904), V(0.010937499999999989, 0.2890625, 0.9572477433702835), 0.3477540013660959, 0.24275225662971645)], []), new PlaneFace(new PlaneSurface(new P3(V(0.942669839890126, 0.33372679389213644, 0), -0.0704401911393759), V(0, 0, -1), V(-0.33372679389213644, 0.942669839890126, 0)), [ new PCurveEdge(new EllipseCurve(V(-0.06640184370318536, -0.0235077791500932, 0), V(0, 0, 0.997516004619599), V(-0.33289781807779234, 0.9403282523625955, 0), 0.02500215049153284, 3.1165905030982604), V(-0.16499999999999998, 0.255, 0.9527591510974849), V(-0.1409375, 0.18703124999999998, 0.9721913045369144), 0.3006922248899664, 0.22581372635962282, null, V(0.31796125684715715, -0.8981373164189183, 0.29544573016418446), V(0.32444628711291806, -0.9164554213903857, 0.22334333850612303)), new StraightEdge(new L3(V(-0.1409375, 0.18703124999999998, 1.2), V(0, 0, -1)), V(-0.1409375, 0.18703124999999998, 0.9721913045369144), V(-0.1409375, 0.18703124999999998, 0.8689970843548256), 0.22780869546308558, 0.33100291564517437), new PCurveEdge(new EllipseCurve(V(-0.06640184370318533, -0.02350777915009319, 0), V(0, 0, -0.8972391985820994), V(-0.299432761097154, 0.8458003316705326, 0), 0.027797112240687767, 3.113795541349105), V(-0.1409375, 0.18703124999999998, 0.8689970843548256), V(-0.16499999999999998, 0.255, 0.8472012747865765), 2.8900247146126086, 2.8060483910328458, null, V(-0.2900076108633503, 0.8191773423737494, -0.2233433385061232), V(-0.282733765215855, 0.7986310900577724, -0.2954457301641845)), new StraightEdge(new L3(V(-0.16499999999999998, 0.255, 1.2), V(0, 0, -1)), V(-0.16499999999999998, 0.255, 0.8472012747865765), V(-0.16499999999999998, 0.255, 0.9527591510974849), 0.3527987252134235, 0.24724084890251508)], []), new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(-0.1409375, 0.18703124999999998, 1.2), V(-0.08703125, 0.2059375, 1.2), V(-0.037968749999999996, 0.2140625, 1.2), V(-0.0029687499999999922, 0.2140625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1409375, 0.18703124999999998, 1.2), V(-0.08703125, 0.2059375, 1.2), V(-0.037968749999999996, 0.2140625, 1.2), V(-0.0029687499999999922, 0.2140625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(-0.1409375, 0.18703124999999998, 0.9721913045369145), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), 1), V(-0.1409375, 0.18703124999999998, 0.9721913045369144), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), 0, 100, null, V(0.001617060512887581, 0.0005671429624910072, 0.00012531587997459948), V(0.0010499999469858546, 0, 0.0000031911732431043194)), new StraightEdge(new L3(V(-0.0029687499999999922, 0.2140625, 1.2), V(0, 0, -1)), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), V(-0.0029687499999999922, 0.2140625, 0.874167279539327), 0.22318454526088316, 0.3258327204606729), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1409375, 0.18703124999999998, 1.2), V(-0.08703125, 0.2059375, 1.2), V(-0.037968749999999996, 0.2140625, 1.2), V(-0.0029687499999999922, 0.2140625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(-0.14093749999999997, 0.18703124999999995, 0.8689970843548255), V(-0.002968749999999992, 0.21406249999999996, 0.8741672795393268), -1), V(-0.0029687499999999922, 0.2140625, 0.874167279539327), V(-0.1409375, 0.18703124999999998, 0.8689970843548256), 100, 0, null, V(-0.001049999933945326, 0, -0.0000035658933671629643), V(-0.0016170285671746015, -0.0005671317583423969, -0.00014019448879870969)), new StraightEdge(new L3(V(-0.1409375, 0.18703124999999998, 1.2), V(0, 0, -1)), V(-0.1409375, 0.18703124999999998, 0.8689970843548256), V(-0.1409375, 0.18703124999999998, 0.9721913045369144), 0.33100291564517437, 0.22780869546308558)], []), new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(-0.0029687499999999922, 0.2140625, 1.2), V(0.07406249999999998, 0.2140625, 1.2), V(0.11093750000000002, 0.18593749999999998, 1.2), V(0.11093750000000002, 0.11203125000000003, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.0029687499999999922, 0.2140625, 1.2), V(0.07406249999999998, 0.2140625, 1.2), V(0.11093750000000002, 0.18593749999999998, 1.2), V(0.11093750000000002, 0.11203125000000003, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), 1), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), 0, 100, null, V(0.0023109369300238886, 0, 0.000007023429018985676), V(0, -0.002216485816645523, 0.00025146076711587485)), new StraightEdge(new L3(V(0.11093750000000002, 0.11203125000000003, 1.2), V(0, 0, -1)), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), V(0.11093750000000002, 0.11203125000000003, 0.8860823156553729), 0.21250728098016458, 0.3139176843446271), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.0029687499999999922, 0.2140625, 1.2), V(0.07406249999999998, 0.2140625, 1.2), V(0.11093750000000002, 0.18593749999999998, 1.2), V(0.11093750000000002, 0.11203125000000003, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(-0.002968749999999992, 0.21406249999999996, 0.8741672795393268), V(0.11093750000000001, 0.11203125000000001, 0.8860823156553727), -1), V(0.11093750000000002, 0.11203125000000003, 0.8860823156553729), V(-0.0029687499999999922, 0.2140625, 0.874167279539327), 100, 0, null, V(0, 0.0022163159641569587, -0.0002802184892673451), V(-0.0023109367883072116, 3.1668370133166326e-21, -0.00000784814731787084)), new StraightEdge(new L3(V(-0.0029687499999999922, 0.2140625, 1.2), V(0, 0, -1)), V(-0.0029687499999999922, 0.2140625, 0.874167279539327), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), 0.3258327204606729, 0.22318454526088316)], []), new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 0.11093750000000002), V(0, 0, -1), V(0, 1, 0)), [ new PCurveEdge(new EllipseCurve(V(0.11093750000000002, 0, 0), V(0, 0, 0.9938273849586506), V(0, 0.9938273849586506, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), V(0.11093750000000002, 0.07093749999999999, 0.9912924604714292), 0.11296719097084001, 0.07143883874192354, null, V(0, -0.9874927190198353, 0.11203125000000001), V(0, -0.9912924604714292, 0.07093749999999999)), new StraightEdge(new L3(V(0.11093750000000002, 0.07093749999999999, 1.2), V(0, 0, -1)), V(0.11093750000000002, 0.07093749999999999, 0.9912924604714292), V(0.11093750000000002, 0.07093749999999999, 0.8903149679677975), 0.2087075395285708, 0.30968503203220243), new PCurveEdge(new EllipseCurve(V(0.11093750000000001, 0, 0), V(0, 0, -0.8931365355273235), V(0, 0.8931365355273235, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.07093749999999999, 0.8903149679677975), V(0.11093750000000002, 0.11203125000000003, 0.8860823156553729), 3.0620837623516257, 3.015825616877026, null, V(0, 0.8903149679677973, -0.0709375), V(0, 0.8860823156553727, -0.11203125000000015)), new StraightEdge(new L3(V(0.11093750000000002, 0.11203125000000003, 1.2), V(0, 0, -1)), V(0.11093750000000002, 0.11203125000000003, 0.8860823156553729), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), 0.3139176843446271, 0.21250728098016458)], []), new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), -0.07093749999999999), V(0, 0, -1), V(1, 0, 0)), [ new PCurveEdge(new EllipseCurve(V(0, 0.07093749999999999, 0), V(0, 0, 0.9974807622674986), V(0.9974807622674986, 0, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.07093749999999999, 0.9912924604714292), V(0.04296875, 0.07093749999999999, 0.9965548442595558), 0.11144825166905915, 0.04309060575909555, null, V(-0.9912924604714292, 0, 0.11093750000000004), V(-0.9965548442595558, 0, 0.04296875)), new StraightEdge(new L3(V(0.04296875, 0.07093749999999999, 1.2), V(0, 0, -1)), V(0.04296875, 0.07093749999999999, 0.9965548442595558), V(0.04296875, 0.07093749999999999, 0.8961704958417165), 0.2034451557404442, 0.30382950415828347), new PCurveEdge(new EllipseCurve(V(0, 0.07093749999999997, 0), V(0, 0, -0.8972000173282154), V(0.8972000173282154, 0, 0), 0, 3.141592653589793), V(0.04296875, 0.07093749999999999, 0.8961704958417165), V(0.11093750000000002, 0.07093749999999999, 0.8903149679677975), 3.093682274615096, 3.0176268184353994, null, V(0.8961704958417163, 0, -0.04296875000000003), V(0.8903149679677973, 0, -0.11093750000000019)), new StraightEdge(new L3(V(0.11093750000000002, 0.07093749999999999, 1.2), V(0, 0, -1)), V(0.11093750000000002, 0.07093749999999999, 0.8903149679677975), V(0.11093750000000002, 0.07093749999999999, 0.9912924604714292), 0.30968503203220243, 0.2087075395285708)], []), new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.04296875, 0.07093749999999999, 0.9965548442595558), V(-0.1710469396532132, 6.938893903907228e-18, 0.9852628808776215), 1), V(0.04296875, 0.07093749999999999, 0.9965548442595558), V(-0.171046939653212, 2.0947208715025655e-17, 0.9852628808776216), 0, 64, null, V(-0.0047989723215347765, 0, 0.0002069187091194757), V(-0.002384621933525827, -0.0028002649734856517, -0.00041398320374781613)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.17104693965321316, 2.0947208715025797e-17, 0.9852628808776215), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), 1), V(-0.171046939653212, 2.0947208715025655e-17, 0.9852628808776216), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), 0, 36, null, V(-0.0018504236016956069, -0.002172955102479456, -0.000321243497826874), V(0, -0.003269989649016917, -0.0003326706415016505)), new StraightEdge(new L3(V(-0.20500000000000002, -0.0990625, 1.2), V(0, 0, -1)), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), V(-0.20500000000000002, -0.0990625, 0.8707247677043248), 0.22626409068282272, 0.3292752322956751), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.17104693965321313, 2.0947208715025794e-17, 0.883596595984429), V(-0.205, -0.09906249999999998, 0.8707247677043246), -1), V(-0.20500000000000002, -0.0990625, 0.8707247677043248), V(-0.17104693965321202, 2.0947208715025658e-17, 0.8835965959844292), 36, 0, null, V(0, 0.0032695078711222625, 0.0003719724481215541), V(0.00215411320093391, 0.002529578236571627, 0.00041699399065229685)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(0.04296874999999999, 0.07093749999999997, 0.8961704958417162), V(-0.17104693965321316, 3.122502256758253e-18, 0.8835965959844289), -1), V(-0.17104693965321202, 2.0947208715025658e-17, 0.8835965959844292), V(0.04296875, 0.07093749999999999, 0.8961704958417165), 64, 0, null, V(0.0023970206921961742, 0.0028148248536625227, 0.00046401610819787214), V(0.004798729293096285, 5.284086361068957e-20, -0.000230085012025602)), new StraightEdge(new L3(V(0.04296875, 0.07093749999999999, 1.2), V(0, 0, -1)), V(0.04296875, 0.07093749999999999, 0.8961704958417165), V(0.04296875, 0.07093749999999999, 0.9965548442595558), 0.30382950415828347, 0.2034451557404442)], []), new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(-0.20500000000000002, -0.0990625, 1.2), V(-0.20500000000000002, -0.19703125, 1.2), V(-0.14, -0.26203125, 1.2), V(-0.034062499999999996, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.20500000000000002, -0.0990625, 1.2), V(-0.20500000000000002, -0.19703125, 1.2), V(-0.14, -0.26203125, 1.2), V(-0.034062499999999996, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), V(-0.034062499999999996, -0.26203125, 0.9644580706890205), 1), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), V(-0.034062499999999996, -0.26203125, 0.9644580706890203), 0, 100, null, V(0, -0.002937749569980493, -0.00029887037541859586), V(0.003177924394426102, 0, 0.00011223717544071708)), new StraightEdge(new L3(V(-0.034062499999999996, -0.26203125, 1.2), V(0, 0, -1)), V(-0.034062499999999996, -0.26203125, 0.9644580706890203), V(-0.034062499999999996, -0.26203125, 0.8603367771501972), 0.23554192931097961, 0.33966322284980277), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.20500000000000002, -0.0990625, 1.2), V(-0.20500000000000002, -0.19703125, 1.2), V(-0.14, -0.26203125, 1.2), V(-0.034062499999999996, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.205, -0.09906249999999998, 0.8707247677043246), V(-0.034062499999999996, -0.26203124999999994, 0.8603367771501969), -1), V(-0.034062499999999996, -0.26203125, 0.8603367771501972), V(-0.20500000000000002, -0.0990625, 0.8707247677043248), 100, 0, null, V(-0.003177872773799387, 0, -0.0001258185099515327), V(0, 0.0029374208172160575, 0.0003341908493916618)), new StraightEdge(new L3(V(-0.20500000000000002, -0.0990625, 1.2), V(0, 0, -1)), V(-0.20500000000000002, -0.0990625, 0.8707247677043248), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), 0.3292752322956751, 0.22626409068282272)], []), new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(-0.034062499999999996, -0.26203125, 1.2), V(0.030000000000000027, -0.26203125, 1.2), V(0.0859375, -0.23796875, 1.2), V(0.12296875000000002, -0.18796875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.034062499999999996, -0.26203125, 1.2), V(0.030000000000000027, -0.26203125, 1.2), V(0.0859375, -0.23796875, 1.2), V(0.12296875000000002, -0.18796875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.034062499999999996, -0.26203125, 0.9644580706890205), V(0.12296875000000002, -0.18796875, 0.9744467330474638), 1), V(-0.034062499999999996, -0.26203125, 0.9644580706890203), V(0.12296875000000002, -0.18796875, 0.9744467330474638), 0, 100, null, V(0.0019218307292635764, 0, 0.00006787475910567421), V(0.0011108137429678576, 0.0014998329018975289, 0.00014913728313140515)), new StraightEdge(new L3(V(0.12296875000000002, -0.18796875, 1.2), V(0, 0, -1)), V(0.12296875000000002, -0.18796875, 0.9744467330474638), V(0.12296875000000002, -0.18796875, 0.871519612829726), 0.2255532669525362, 0.32848038717027395), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.034062499999999996, -0.26203125, 1.2), V(0.030000000000000027, -0.26203125, 1.2), V(0.0859375, -0.23796875, 1.2), V(0.12296875000000002, -0.18796875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.034062499999999996, -0.26203124999999994, 0.8603367771501969), V(0.12296874999999999, -0.18796874999999996, 0.8715196128297258), -1), V(0.12296875000000002, -0.18796875, 0.871519612829726), V(-0.034062499999999996, -0.26203125, 0.8603367771501972), 100, 0, null, V(-0.0011107827301533816, -0.0014997910280551989, -0.0001667458526657328), V(-0.0019218193657008743, 0, -0.00007608877579431633)), new StraightEdge(new L3(V(-0.034062499999999996, -0.26203125, 1.2), V(0, 0, -1)), V(-0.034062499999999996, -0.26203125, 0.8603367771501972), V(-0.034062499999999996, -0.26203125, 0.9644580706890203), 0.33966322284980277, 0.23554192931097961)], []), new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.12296875000000002, -0.18796875, 1.2), V(0.13906249999999998, -0.235, 1.2), V(0.17296875, -0.2559375, 1.2), V(0.21999999999999997, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.12296875000000002, -0.18796875, 1.2), V(0.13906249999999998, -0.235, 1.2), V(0.17296875, -0.2559375, 1.2), V(0.21999999999999997, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.12296875000000002, -0.18796875, 0.9744467330474638), V(0.21999999999999997, -0.26203125, 0.9396486705271484), 1), V(0.12296875000000002, -0.18796875, 0.9744467330474638), V(0.21999999999999997, -0.26203125, 0.9396486705271484), 0, 100, null, V(0.0004825448774820869, -0.0014101554186612473, -0.00033291003063843776), V(0.0014098883430614425, -0.00018267656272224708, -0.0003810385889743954)), new StraightEdge(new L3(V(0.21999999999999997, -0.26203125, 1.2), V(0, 0, -1)), V(0.21999999999999997, -0.26203125, 0.9396486705271484), V(0.21999999999999997, -0.26203125, 0.8324299514214022), 0.26035132947285156, 0.36757004857859776), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.12296875000000002, -0.18796875, 1.2), V(0.13906249999999998, -0.235, 1.2), V(0.17296875, -0.2559375, 1.2), V(0.21999999999999997, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(0.12296874999999999, -0.18796874999999996, 0.8715196128297258), V(0.21999999999999995, -0.2620312499999999, 0.832429951421402), -1), V(0.21999999999999997, -0.26203125, 0.8324299514214022), V(0.12296875000000002, -0.18796875, 0.871519612829726), 100, 0, null, V(-0.0014095926581339242, 0.00018263825138612356, 0.00043002695120080944), V(-0.0004824780014564421, 0.0014099599848387352, 0.0003721753680202314)), new StraightEdge(new L3(V(0.12296875000000002, -0.18796875, 1.2), V(0, 0, -1)), V(0.12296875000000002, -0.18796875, 0.871519612829726), V(0.12296875000000002, -0.18796875, 0.9744467330474638), 0.32848038717027395, 0.2255532669525362)], []), new PlaneFace(new PlaneSurface(new P3(V(-0.9456422745519422, 0.325208684662986, 0), -0.2932561385545256), V(0, 0, -1), V(-0.325208684662986, -0.9456422745519422, 0)), [ new PCurveEdge(new EllipseCurve(V(0.27731540188902115, -0.09536944308866366, 0), V(0.31091053038645644, 0.9040660812655796, 0), V(0, 0, 0.9560339100680944), 1.465110231866663, 3.141592653589793), V(0.21999999999999997, -0.26203125, 0.9396486705271484), V(0.24203124999999998, -0.19796875, 0.9498574882827818), 1.7562036900306723, 1.684527864719763, null, V(0.30558190818745745, 0.8885715060770011, 0.17624191662783179), V(0.30890190438173115, 0.8982253957199245, 0.10849695458036647)), new StraightEdge(new L3(V(0.24203124999999998, -0.19796875, 1.2), V(0, 0, -1)), V(0.24203124999999998, -0.19796875, 0.9498574882827818), V(0.24203124999999998, -0.19796875, 0.8439367559520533), 0.25014251171721813, 0.35606324404794665), new PCurveEdge(new EllipseCurve(V(0.2773154018890211, -0.09536944308866366, 0), V(-0.27671434201165657, -0.8046303562041052, 0), V(0, 0, 0.8508823874073836), 0, 1.6896013946143438), V(0.24203124999999998, -0.19796875, 0.8439367559520533), V(0.21999999999999997, -0.26203125, 0.8324299514214022), 1.4429371319818958, 1.3621575275257576, null, V(-0.2744555623419146, -0.7980622734764866, -0.10849695458036643), V(-0.27071344957582744, -0.7871809526672972, -0.17624191662783154)), new StraightEdge(new L3(V(0.21999999999999997, -0.26203125, 1.2), V(0, 0, -1)), V(0.21999999999999997, -0.26203125, 0.8324299514214022), V(0.21999999999999997, -0.26203125, 0.9396486705271484), 0.36757004857859776, 0.26035132947285156)], []), new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.24203124999999998, -0.19796875, 1.2), V(0.21500000000000002, -0.18796875, 1.2), V(0.20296874999999998, -0.17296875, 1.2), V(0.20296874999999998, -0.13093749999999998, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.24203124999999998, -0.19796875, 1.2), V(0.21500000000000002, -0.18796875, 1.2), V(0.20296874999999998, -0.17296875, 1.2), V(0.20296874999999998, -0.13093749999999998, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.24203124999999998, -0.19796875, 0.9498574882827818), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), 1), V(0.24203124999999998, -0.19796875, 0.9498574882827818), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), 0, 100, null, V(-0.0008106439109762247, 0.00029989138903166777, 0.0002690617125763413), V(0, 0.0012607577897218358, 0.00017011744865842043)), new StraightEdge(new L3(V(0.20296874999999998, -0.13093749999999998, 1.2), V(0, 0, -1)), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), V(0.20296874999999998, -0.13093749999999998, 0.8669827320178805), 0.22960881206742834, 0.3330172679821195), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.24203124999999998, -0.19796875, 1.2), V(0.21500000000000002, -0.18796875, 1.2), V(0.20296874999999998, -0.17296875, 1.2), V(0.20296874999999998, -0.13093749999999998, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(0.2420312499999999, -0.19796874999999997, 0.843936755952053), V(0.20296874999999995, -0.13093749999999996, 0.8669827320178802), -1), V(0.20296874999999998, -0.13093749999999998, 0.8669827320178805), V(0.24203124999999998, -0.19796875, 0.8439367559520533), 100, 0, null, V(-1.6025994766705188e-19, -0.0012607132877938707, -0.00019040130792028905), V(0.0008105656446711983, -0.000299862435022873, -0.00030280184601271074)), new StraightEdge(new L3(V(0.24203124999999998, -0.19796875, 1.2), V(0, 0, -1)), V(0.24203124999999998, -0.19796875, 0.8439367559520533), V(0.24203124999999998, -0.19796875, 0.9498574882827818), 0.35606324404794665, 0.25014251171721813)], []), new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), -0.20296874999999998), V(0, 0, -1), V(0, -1, 0)), [ new PCurveEdge(new EllipseCurve(V(0.20296874999999998, 0, 0), V(0, 0, -0.9791852156376941), V(0, 0.9791852156376941, 0), 0, 3.141592653589793), V(0.2029687500000001, 0, 0.9791852156376941), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), 3.141592653589793, 3.0217872455230754, null, V(0, 0.9791852156376941, -1.199156040103113e-16), V(0, 0.9721663299286162, -0.1170312499999999)), new StraightEdge(new L3(V(0.20296874999999998, 0.11703124999999998, 1.2), V(0, 0, -1)), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), V(0.20296874999999998, 0.11703124999999998, 0.8689691438980301), 0.22783367007138378, 0.3310308561019698), new PCurveEdge(new EllipseCurve(V(0.20296874999999995, 0, 0), V(0, 0, 0.8768145108992196), V(0, 0.8768145108992196, 0), 0, 3.141592653589793), V(0.20296874999999998, 0.11703124999999998, 0.8689691438980301), V(0.20296875000000006, 0, 0.8768145108992196), 0.13387273064879698, 0, null, V(0, -0.86896914389803, 0.11703124999999995), V(0, -0.8768145108992196, 0)), new PCurveEdge(new EllipseCurve(V(0.20296874999999995, 4.943121956121308e-19, 0), V(0, 0, -0.8768145108992196), V(2.1354031397797386e-18, -0.8768145108992196, 0), 2.891240380027181e-17, 3.141592653589793), V(0.20296875000000006, 0, 0.8768145108992196), V(0.20296874999999998, -0.13093749999999998, 0.8669827320178805), 3.141592653589793, 2.9916987961562977, null, V(2.1354031397797386e-18, -0.8768145108992196, -1.0737880842186813e-16), V(2.111458723678207e-18, -0.8669827320178803, -0.13093749999999996)), new StraightEdge(new L3(V(0.20296874999999998, -0.13093749999999998, 1.2), V(0, 0, -1)), V(0.20296874999999998, -0.13093749999999998, 0.8669827320178805), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), 0.3330172679821195, 0.22960881206742834), new PCurveEdge(new EllipseCurve(V(0.20296874999999998, -2.3224450485052117e-18, 0), V(0, 0, 0.9791852156376941), V(-1.1204206832959602e-17, -0.9791852156376941, 0), 2.3013070043406876e-17, 3.141592653589793), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), V(0.2029687500000001, 0, 0.9791852156376941), 0.13412262887088133, 2.3013070043406876e-17, null, V(1.1103582248632314e-17, 0.9703911879325716, 0.13093749999999996), V(1.1204206832959602e-17, 0.9791852156376941, 2.3224450485052132e-18))], []), new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(-0.009062500000000001, -0.19296875, 1.2), V(-0.07500000000000001, -0.19296875, 1.2), V(-0.1040625, -0.16, 1.2), V(-0.1040625, -0.095, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.009062500000000001, -0.19296875, 1.2), V(-0.07500000000000001, -0.19296875, 1.2), V(-0.1040625, -0.16, 1.2), V(-0.1040625, -0.095, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), V(-0.1040625, -0.095, 0.9900232300778351), 1), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), V(-0.1040625, -0.095, 0.9900232300778351), 0, 100, null, V(-0.001978121698253082, 0, -0.000018270895823513832), V(0, 0.0019496588340396214, 0.00018708408409689773)), new StraightEdge(new L3(V(-0.1040625, -0.095, 1.2), V(0, 0, -1)), V(-0.1040625, -0.095, 0.9900232300778351), V(-0.1040625, -0.095, 0.8889015671567637), 0.20997676992216485, 0.31109843284323624), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.009062500000000001, -0.19296875, 1.2), V(-0.07500000000000001, -0.19296875, 1.2), V(-0.1040625, -0.16, 1.2), V(-0.1040625, -0.095, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.0090625, -0.19296874999999997, 0.879022714505824), V(-0.10406249999999999, -0.09499999999999999, 0.8889015671567634), -1), V(-0.1040625, -0.095, 0.8889015671567637), V(-0.009062500000000001, -0.19296875, 0.8790227145058244), 100, 0, null, V(0, -0.0019495768662389622, -0.0002083580557576388), V(0.001978120886365625, 1.734719868511089e-20, 0.000020393921837124158)), new StraightEdge(new L3(V(-0.009062500000000001, -0.19296875, 1.2), V(0, 0, -1)), V(-0.009062500000000001, -0.19296875, 0.8790227145058244), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), 0.3209772854941756, 0.21883694901551276)], []), new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(-0.03286362438479655, -2.6020852139652106e-18, 0.9994598452125503), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), 1), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), 0, 29, null, V(0.0022900449065656465, 0.0007030314524671974, 0.00007529984920773373), V(0.003388522845144273, 0, -0.0001627386905160511)), new StraightEdge(new L3(V(0.047968750000000004, 0.010000000000000009, 1.2), V(0, 0, -1)), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), V(0.047968750000000004, 0.010000000000000009, 0.8986651206224918), 0.20120122195537427, 0.30133487937750814), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(-0.03286362438479655, -2.3418766925686897e-18, 0.8993997899667839), V(0.04796875, 0.010000000000000007, 0.8986651206224917), -1), V(0.047968750000000004, 0.010000000000000009, 0.8986651206224918), V(-0.03286362438479702, 4.024633241122258e-18, 0.8993997899667838), 29, 0, null, V(-0.003388395809325897, 0, 0.00018086504944802423), V(-0.002292032896076846, -0.0007036417545401136, -0.00008374975068322923)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.10406249999999999, -0.09499999999999999, 0.8889015671567634), V(-0.03286362438479653, 4.0246332411221976e-18, 0.8993997899667839), -1), V(-0.03286362438479702, 4.024633241122258e-18, 0.8993997899667838), V(-0.1040625, -0.095, 0.8889015671567637), 71, 0, null, V(-0.0022509978133461728, -0.0006910442051507778, -0.00008225034901502365), V(0, -0.0019214697267915689, -0.00020535414807408821)), new StraightEdge(new L3(V(-0.1040625, -0.095, 1.2), V(0, 0, -1)), V(-0.1040625, -0.095, 0.8889015671567637), V(-0.1040625, -0.095, 0.9900232300778351), 0.31109843284323624, 0.20997676992216485), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.1040625, -0.095, 0.9900232300778351), V(-0.03286362438479653, 2.7235906341395923e-18, 0.9994598452125503), 1), V(-0.1040625, -0.095, 0.9900232300778351), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), 0, 71, null, V(0, 0.0019215482684322, 0.00018438666887312053), V(0.0022530279161436587, 0.0006916674357757679, 0.00007408267927848543))], []), new PlaneFace(new PlaneSurface(new P3(V(0, 1, 0), 0.010000000000000009), V(0, 0, -1), V(-1, 0, 0)), [ new PCurveEdge(new EllipseCurve(V(0, 0.010000000000000009, 0), V(0, 0, -0.9999499987499375), V(0.9999499987499375, 0, 0), 0, 3.141592653589793), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), 3.0936030871101416, 3.0304207486077566, null, V(0.9987987780446257, 0, -0.0479687500000002), V(0.9937770731375071, 0, -0.11093749999999983)), new StraightEdge(new L3(V(0.11093750000000002, 0.010000000000000009, 1.2), V(0, 0, -1)), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), V(0.11093750000000002, 0.010000000000000009, 0.8930805512907278), 0.20622292686249288, 0.3069194487092721), new PCurveEdge(new EllipseCurve(V(0, 0.010000000000000007, 0), V(0, 0, 0.899944442729661), V(0.899944442729661, 0, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.010000000000000009, 0.8930805512907278), V(0.047968750000000004, 0.010000000000000009, 0.8986651206224918), 0.12358585566041683, 0.053327173050216115, null, V(-0.8930805512907277, 0, 0.11093750000000001), V(-0.8986651206224918, 0, 0.04796875)), new StraightEdge(new L3(V(0.047968750000000004, 0.010000000000000009, 1.2), V(0, 0, -1)), V(0.047968750000000004, 0.010000000000000009, 0.8986651206224918), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), 0.30133487937750814, 0.20120122195537427)], []), new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 0.11093750000000002), V(0, 0, -1), V(0, 1, 0)), [ new PCurveEdge(new EllipseCurve(V(0.11093750000000002, 0, 0), V(0, 0, 0.9938273849586506), V(0, 0.9938273849586506, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), V(0.11093749999999993, 0, 0.9938273849586506), 0.010062279327831384, 0, null, V(0, -0.993777073137507, 0.010000000000000007), V(0, -0.9938273849586506, 0)), new PCurveEdge(new EllipseCurve(V(0.11093750000000002, 2.7017833632379675e-19, 0), V(0, 0, -0.9938273849586506), V(2.4203775050019737e-18, -0.9938273849586506, 0), 1.3942163371701867e-17, 3.141592653589793), V(0.11093749999999993, 0, 0.9938273849586506), V(0.11093750000000002, -0.12, 0.9865560658643532), 3.141592653589793, 3.02055199779192, null, V(2.4203775050019737e-18, -0.9938273849586506, -1.217087525894596e-16), V(2.4026688591808874e-18, -0.9865560658643532, -0.11999999999999994)), new StraightEdge(new L3(V(0.11093750000000002, -0.12, 1.2), V(0, 0, -1)), V(0.11093750000000002, -0.12, 0.9865560658643532), V(0.11093750000000002, -0.12, 0.8850383444200315), 0.21344393413564677, 0.31496165557996847), new PCurveEdge(new EllipseCurve(V(0.11093750000000001, 2.701783363237952e-19, 0), V(0, 0, 0.8931365355273235), V(2.175153967583285e-18, -0.8931365355273235, 0), 1.5513981584219775e-17, 3.141592653589793), V(0.11093750000000002, -0.12, 0.8850383444200315), V(0.11093750000000006, 0, 0.8931365355273235), 0.13476551593405545, 1.5513981584219775e-17, null, V(-2.155431549099e-18, 0.8850383444200314, 0.11999999999999998), V(-2.175153967583285e-18, 0.8931365355273235, 2.7017833632379536e-19)), new PCurveEdge(new EllipseCurve(V(0.11093750000000001, 0, 0), V(0, 0, -0.8931365355273235), V(0, 0.8931365355273235, 0), 0, 3.141592653589793), V(0.11093750000000006, 0, 0.8931365355273235), V(0.11093750000000002, 0.010000000000000009, 0.8930805512907278), 3.141592653589793, 3.130395923246913, null, V(0, 0.8931365355273235, -1.093776799435093e-16), V(0, 0.8930805512907277, -0.01000000000000018)), new StraightEdge(new L3(V(0.11093750000000002, 0.010000000000000009, 1.2), V(0, 0, -1)), V(0.11093750000000002, 0.010000000000000009, 0.8930805512907278), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), 0.3069194487092721, 0.20622292686249288)], []), new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0.11093750000000002, -0.12, 1.2), V(0.0809375, -0.16593750000000002, 1.2), V(0.040937500000000016, -0.19296875, 1.2), V(-0.009062500000000001, -0.19296875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.11093750000000002, -0.12, 1.2), V(0.0809375, -0.16593750000000002, 1.2), V(0.040937500000000016, -0.19296875, 1.2), V(-0.009062500000000001, -0.19296875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.11093750000000002, -0.12, 0.9865560658643532), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), 1), V(0.11093750000000002, -0.12, 0.9865560658643532), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), 0, 100, null, V(-0.0008999801458405889, -0.001378094598318401, -0.00006642278795539617), V(-0.001499998572758444, 0, -0.000013854717676112659)), new StraightEdge(new L3(V(-0.009062500000000001, -0.19296875, 1.2), V(0, 0, -1)), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), V(-0.009062500000000001, -0.19296875, 0.8790227145058244), 0.21883694901551276, 0.3209772854941756), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.11093750000000002, -0.12, 1.2), V(0.0809375, -0.16593750000000002, 1.2), V(0.040937500000000016, -0.19296875, 1.2), V(-0.009062500000000001, -0.19296875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(0.11093750000000001, -0.11999999999999997, 0.8850383444200313), V(-0.0090625, -0.19296874999999997, 0.879022714505824), -1), V(-0.009062500000000001, -0.19296875, 0.8790227145058244), V(0.11093750000000002, -0.12, 0.8850383444200315), 100, 0, null, V(0.0014999982255756404, 0, 0.000015464599145110234), V(0.0008999753301001092, 0.0013780872242157918, 0.00007404137248522919)), new StraightEdge(new L3(V(0.11093750000000002, -0.12, 1.2), V(0, 0, -1)), V(0.11093750000000002, -0.12, 0.8850383444200315), V(0.11093750000000002, -0.12, 0.9865560658643532), 0.31496165557996847, 0.21344393413564677)], []), new RotationFace(new EllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.20296874999999998, 0.11703124999999998, 1.2), V(0.20296874999999998, 0.2240625, 1.2), V(0.14500000000000002, 0.2890625, 1.2), V(0.010937499999999989, 0.2890625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(0.20296874999999995, 0.11703124999999996, 0.8689691438980299), V(0.010937499999999987, 0.2890625, 0.8522459986339039), -1), V(0.20296874999999998, 0.11703124999999998, 0.8689691438980301), V(0.010937499999999989, 0.2890625, 0.852245998633904), 0, 100, null, V(0, 0.0032079393598033606, -0.0004320396826956172), V(-0.004021823173538629, 0, 0.0000516150161233843)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.010937499999999989, 0.2890625, 1.2), V(-0.04093749999999999, 0.2890625, 1.2), V(-0.1, 0.2790625, 1.2), V(-0.16499999999999998, 0.255, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(0.010937499999999987, 0.2890625, 0.8522459986339039), V(-0.16499999999999998, 0.25499999999999995, 0.8472012747865764), -1), V(0.010937499999999989, 0.2890625, 0.852245998633904), V(-0.16499999999999998, 0.255, 0.8472012747865765), 0, 100, null, V(-0.0015562468960623505, -6.688077170406378e-21, 0.00001997246153454078), V(-0.001949741977734284, -0.0007217794821420186, -0.0001624804665392168)), new PCurveEdge(new EllipseCurve(V(-0.06640184370318533, -0.02350777915009319, 0), V(0, 0, -0.8972391985820994), V(-0.299432761097154, 0.8458003316705326, 0), 0.027797112240687767, 3.113795541349105), V(-0.16499999999999998, 0.255, 0.8472012747865765), V(-0.1409375, 0.18703124999999998, 0.8689970843548256), 2.8060483910328458, 2.8900247146126086, null, V(0.282733765215855, -0.7986310900577724, 0.2954457301641845), V(0.2900076108633503, -0.8191773423737494, 0.2233433385061232)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1409375, 0.18703124999999998, 1.2), V(-0.08703125, 0.2059375, 1.2), V(-0.037968749999999996, 0.2140625, 1.2), V(-0.0029687499999999922, 0.2140625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(-0.14093749999999997, 0.18703124999999995, 0.8689970843548255), V(-0.002968749999999992, 0.21406249999999996, 0.8741672795393268), -1), V(-0.1409375, 0.18703124999999998, 0.8689970843548256), V(-0.0029687499999999922, 0.2140625, 0.874167279539327), 0, 100, null, V(0.0016170285671746015, 0.0005671317583423969, 0.00014019448879870969), V(0.001049999933945326, 0, 0.0000035658933671629643)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.0029687499999999922, 0.2140625, 1.2), V(0.07406249999999998, 0.2140625, 1.2), V(0.11093750000000002, 0.18593749999999998, 1.2), V(0.11093750000000002, 0.11203125000000003, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(-0.002968749999999992, 0.21406249999999996, 0.8741672795393268), V(0.11093750000000001, 0.11203125000000001, 0.8860823156553727), -1), V(-0.0029687499999999922, 0.2140625, 0.874167279539327), V(0.11093750000000002, 0.11203125000000003, 0.8860823156553729), 0, 100, null, V(0.0023109367883072116, -3.1668370133166326e-21, 0.00000784814731787084), V(0, -0.0022163159641569587, 0.0002802184892673451)), new PCurveEdge(new EllipseCurve(V(0.11093750000000001, 0, 0), V(0, 0, -0.8931365355273235), V(0, 0.8931365355273235, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.11203125000000003, 0.8860823156553729), V(0.11093750000000002, 0.07093749999999999, 0.8903149679677975), 3.015825616877026, 3.0620837623516257, null, V(0, -0.8860823156553727, 0.11203125000000015), V(0, -0.8903149679677973, 0.0709375)), new PCurveEdge(new EllipseCurve(V(0, 0.07093749999999997, 0), V(0, 0, -0.8972000173282154), V(0.8972000173282154, 0, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.07093749999999999, 0.8903149679677975), V(0.04296875, 0.07093749999999999, 0.8961704958417165), 3.0176268184353994, 3.093682274615096, null, V(-0.8903149679677973, 0, 0.11093750000000019), V(-0.8961704958417163, 0, 0.04296875000000003)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(0.04296874999999999, 0.07093749999999997, 0.8961704958417162), V(-0.17104693965321316, 3.122502256758253e-18, 0.8835965959844289), -1), V(0.04296875, 0.07093749999999999, 0.8961704958417165), V(-0.17104693965321202, 2.0947208715025658e-17, 0.8835965959844292), 0, 64, null, V(-0.004798729293096285, -5.284086361068957e-20, 0.000230085012025602), V(-0.0023970206921961742, -0.0028148248536625227, -0.00046401610819787214)), new PCurveEdge(new EllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(-0.9, 1.1021821192326179e-16, 0), 0, 3.141592653589793), V(-0.17104693965321202, 2.0947208715025658e-17, 0.8835965959844292), V(0, 0, -0.9), 2.9503773839343412, 0, null, V(-0.8835965959844292, 1.0820937430098282e-16, -0.17104693965321202), V(0.9, -1.1021821192326179e-16, 0)), new PCurveEdge(new EllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(0.9, 0, 0)), V(0, 0, -0.9), V(0.20296875000000006, 0, 0.8768145108992196), 0, 2.9141150448015316, null, V(0.9, 0, 0), V(-0.8768145108992196, 0, 0.20296875000000006)), new PCurveEdge(new EllipseCurve(V(0.20296874999999995, 0, 0), V(0, 0, 0.8768145108992196), V(0, 0.8768145108992196, 0), 0, 3.141592653589793), V(0.20296875000000006, 0, 0.8768145108992196), V(0.20296874999999998, 0.11703124999999998, 0.8689691438980301), 0, 0.13387273064879698, null, V(0, 0.8768145108992196, 0), V(0, 0.86896914389803, -0.11703124999999995))], []), new RotationFace(new EllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(0.9, 0, 0), V(0, 0.9, 0), V(0, 0, -0.9)), V(-0.03286362438479655, -2.3418766925686897e-18, 0.8993997899667839), V(0.04796875, 0.010000000000000007, 0.8986651206224917), -1), V(-0.03286362438479702, 4.024633241122258e-18, 0.8993997899667838), V(0.047968750000000004, 0.010000000000000009, 0.8986651206224918), 0, 29, null, V(0.002292032896076846, 0.0007036417545401136, 0.00008374975068322923), V(0.003388395809325897, 0, -0.00018086504944802423)), new PCurveEdge(new EllipseCurve(V(0, 0.010000000000000007, 0), V(0, 0, 0.899944442729661), V(0.899944442729661, 0, 0), 0, 3.141592653589793), V(0.047968750000000004, 0.010000000000000009, 0.8986651206224918), V(0.11093750000000002, 0.010000000000000009, 0.8930805512907278), 0.053327173050216115, 0.12358585566041683, null, V(0.8986651206224918, 0, -0.04796875), V(0.8930805512907277, 0, -0.11093750000000001)), new PCurveEdge(new EllipseCurve(V(0.11093750000000001, 0, 0), V(0, 0, -0.8931365355273235), V(0, 0.8931365355273235, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.010000000000000009, 0.8930805512907278), V(0.11093750000000006, 0, 0.8931365355273235), 3.130395923246913, 3.141592653589793, null, V(0, -0.8930805512907277, 0.01000000000000018), V(0, -0.8931365355273235, 1.093776799435093e-16)), new PCurveEdge(new EllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(0.9, 0, 0)), V(0.11093750000000006, 0, 0.8931365355273235), V(0, 0, 0.9), 3.018014465996861, 3.141592653589793, null, V(-0.8931365355273235, 0, 0.11093750000000006), V(-0.9, 0, 0)), new PCurveEdge(new EllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(-0.9, 1.1021821192326179e-16, 0), 0, 3.141592653589793), V(0, 0, 0.9), V(-0.03286362438479702, 4.024633241122258e-18, 0.8993997899667838), 3.141592653589793, 3.1050693959027966, null, V(-0.9, 1.1021821192326179e-16, 0), V(-0.8993997899667838, 1.1014470739366237e-16, -0.03286362438479702))], []), new RotationFace(new EllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.17104693965321313, 2.0947208715025794e-17, 0.883596595984429), V(-0.205, -0.09906249999999998, 0.8707247677043246), -1), V(-0.17104693965321202, 2.0947208715025658e-17, 0.8835965959844292), V(-0.20500000000000002, -0.0990625, 0.8707247677043248), 0, 36, null, V(-0.00215411320093391, -0.002529578236571627, -0.00041699399065229685), V(0, -0.0032695078711222625, -0.0003719724481215541)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.20500000000000002, -0.0990625, 1.2), V(-0.20500000000000002, -0.19703125, 1.2), V(-0.14, -0.26203125, 1.2), V(-0.034062499999999996, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.205, -0.09906249999999998, 0.8707247677043246), V(-0.034062499999999996, -0.26203124999999994, 0.8603367771501969), -1), V(-0.20500000000000002, -0.0990625, 0.8707247677043248), V(-0.034062499999999996, -0.26203125, 0.8603367771501972), 0, 100, null, V(0, -0.0029374208172160575, -0.0003341908493916618), V(0.003177872773799387, 0, 0.0001258185099515327)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.034062499999999996, -0.26203125, 1.2), V(0.030000000000000027, -0.26203125, 1.2), V(0.0859375, -0.23796875, 1.2), V(0.12296875000000002, -0.18796875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.034062499999999996, -0.26203124999999994, 0.8603367771501969), V(0.12296874999999999, -0.18796874999999996, 0.8715196128297258), -1), V(-0.034062499999999996, -0.26203125, 0.8603367771501972), V(0.12296875000000002, -0.18796875, 0.871519612829726), 0, 100, null, V(0.0019218193657008743, 0, 0.00007608877579431633), V(0.0011107827301533816, 0.0014997910280551989, 0.0001667458526657328)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.12296875000000002, -0.18796875, 1.2), V(0.13906249999999998, -0.235, 1.2), V(0.17296875, -0.2559375, 1.2), V(0.21999999999999997, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(0.12296874999999999, -0.18796874999999996, 0.8715196128297258), V(0.21999999999999995, -0.2620312499999999, 0.832429951421402), -1), V(0.12296875000000002, -0.18796875, 0.871519612829726), V(0.21999999999999997, -0.26203125, 0.8324299514214022), 0, 100, null, V(0.0004824780014564421, -0.0014099599848387352, -0.0003721753680202314), V(0.0014095926581339242, -0.00018263825138612356, -0.00043002695120080944)), new PCurveEdge(new EllipseCurve(V(0.2773154018890211, -0.09536944308866366, 0), V(-0.27671434201165657, -0.8046303562041052, 0), V(0, 0, 0.8508823874073836), 0, 1.6896013946143438), V(0.21999999999999997, -0.26203125, 0.8324299514214022), V(0.24203124999999998, -0.19796875, 0.8439367559520533), 1.3621575275257576, 1.4429371319818958, null, V(0.27071344957582744, 0.7871809526672972, 0.17624191662783154), V(0.2744555623419146, 0.7980622734764866, 0.10849695458036643)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.24203124999999998, -0.19796875, 1.2), V(0.21500000000000002, -0.18796875, 1.2), V(0.20296874999999998, -0.17296875, 1.2), V(0.20296874999999998, -0.13093749999999998, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(0.2420312499999999, -0.19796874999999997, 0.843936755952053), V(0.20296874999999995, -0.13093749999999996, 0.8669827320178802), -1), V(0.24203124999999998, -0.19796875, 0.8439367559520533), V(0.20296874999999998, -0.13093749999999998, 0.8669827320178805), 0, 100, null, V(-0.0008105656446711983, 0.000299862435022873, 0.00030280184601271074), V(1.6025994766705188e-19, 0.0012607132877938707, 0.00019040130792028905)), new PCurveEdge(new EllipseCurve(V(0.20296874999999995, 4.943121956121308e-19, 0), V(0, 0, -0.8768145108992196), V(2.1354031397797386e-18, -0.8768145108992196, 0), 2.891240380027181e-17, 3.141592653589793), V(0.20296874999999998, -0.13093749999999998, 0.8669827320178805), V(0.20296875000000006, 0, 0.8768145108992196), 2.9916987961562977, 3.141592653589793, null, V(-2.111458723678207e-18, 0.8669827320178803, 0.13093749999999996), V(-2.1354031397797386e-18, 0.8768145108992196, 1.0737880842186813e-16)), new PCurveEdge(new EllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(0.9, 0, 0)), V(0.20296875000000006, 0, 0.8768145108992196), V(0, 0, -0.9), 2.9141150448015316, 0, null, V(0.8768145108992196, 0, -0.20296875000000006), V(-0.9, 0, 0)), new PCurveEdge(new EllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(-0.9, 1.1021821192326179e-16, 0), 0, 3.141592653589793), V(0, 0, -0.9), V(-0.17104693965321202, 2.0947208715025658e-17, 0.8835965959844292), 0, 2.9503773839343412, null, V(-0.9, 1.1021821192326179e-16, 0), V(0.8835965959844292, -1.0820937430098282e-16, 0.17104693965321202))], []), new RotationFace(new EllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.009062500000000001, -0.19296875, 1.2), V(-0.07500000000000001, -0.19296875, 1.2), V(-0.1040625, -0.16, 1.2), V(-0.1040625, -0.095, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.0090625, -0.19296874999999997, 0.879022714505824), V(-0.10406249999999999, -0.09499999999999999, 0.8889015671567634), -1), V(-0.009062500000000001, -0.19296875, 0.8790227145058244), V(-0.1040625, -0.095, 0.8889015671567637), 0, 100, null, V(-0.001978120886365625, -1.734719868511089e-20, -0.000020393921837124158), V(0, 0.0019495768662389622, 0.0002083580557576388)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(-0.10406249999999999, -0.09499999999999999, 0.8889015671567634), V(-0.03286362438479653, 4.0246332411221976e-18, 0.8993997899667839), -1), V(-0.1040625, -0.095, 0.8889015671567637), V(-0.03286362438479702, 4.024633241122258e-18, 0.8993997899667838), 0, 71, null, V(0, 0.0019214697267915689, 0.00020535414807408821), V(0.0022509978133461728, 0.0006910442051507778, 0.00008225034901502365)), new PCurveEdge(new EllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(-0.9, 1.1021821192326179e-16, 0), 0, 3.141592653589793), V(-0.03286362438479702, 4.024633241122258e-18, 0.8993997899667838), V(0, 0, 0.9), 3.1050693959027966, 3.141592653589793, null, V(0.8993997899667838, -1.1014470739366237e-16, 0.03286362438479702), V(0.9, -1.1021821192326179e-16, 0)), new PCurveEdge(new EllipseCurve(V(0, 0, 0), V(0, 0, -0.9), V(0.9, 0, 0)), V(0, 0, 0.9), V(0.11093750000000006, 0, 0.8931365355273235), 3.141592653589793, 3.018014465996861, null, V(0.9, 0, 0), V(0.8931365355273235, 0, -0.11093750000000006)), new PCurveEdge(new EllipseCurve(V(0.11093750000000001, 2.701783363237952e-19, 0), V(0, 0, 0.8931365355273235), V(2.175153967583285e-18, -0.8931365355273235, 0), 1.5513981584219775e-17, 3.141592653589793), V(0.11093750000000006, 0, 0.8931365355273235), V(0.11093750000000002, -0.12, 0.8850383444200315), 1.5513981584219775e-17, 0.13476551593405545, null, V(2.175153967583285e-18, -0.8931365355273235, -2.7017833632379536e-19), V(2.155431549099e-18, -0.8850383444200314, -0.11999999999999998)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.11093750000000002, -0.12, 1.2), V(0.0809375, -0.16593750000000002, 1.2), V(0.040937500000000016, -0.19296875, 1.2), V(-0.009062500000000001, -0.19296875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-0.9, 1.1021821192326179e-16, 0), V(-1.1021821192326179e-16, -0.9, 0), V(0, 0, -0.9)), V(0.11093750000000001, -0.11999999999999997, 0.8850383444200313), V(-0.0090625, -0.19296874999999997, 0.879022714505824), -1), V(0.11093750000000002, -0.12, 0.8850383444200315), V(-0.009062500000000001, -0.19296875, 0.8790227145058244), 0, 100, null, V(-0.0008999753301001092, -0.0013780872242157918, -0.00007404137248522919), V(-0.0014999982255756404, 0, -0.000015464599145110234))], []), new RotationFace(new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.20296874999999998, 0.11703124999999998, 1.2), V(0.20296874999999998, 0.2240625, 1.2), V(0.14500000000000002, 0.2890625, 1.2), V(0.010937499999999989, 0.2890625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), V(0.010937499999999989, 0.2890625, 0.9572477433702835), 1), V(0.010937499999999989, 0.2890625, 0.9572477433702835), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), 100, 0, null, V(0.004021833637101672, 0, -0.00004595341771287279), V(0, -0.003208541417521995, 0.0003862503783657504)), new PCurveEdge(new EllipseCurve(V(0.20296874999999998, 0, 0), V(0, 0, -0.9791852156376941), V(0, 0.9791852156376941, 0), 0, 3.141592653589793), V(0.20296874999999998, 0.11703124999999998, 0.9721663299286162), V(0.2029687500000001, 0, 0.9791852156376941), 3.0217872455230754, 3.141592653589793, null, V(0, -0.9721663299286162, 0.1170312499999999), V(0, -0.9791852156376941, 1.199156040103113e-16)), new PCurveEdge(new EllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0.2029687500000001, 0, 0.9791852156376941), V(0, 0, -1), 2.937203822794261, 0, null, V(0.9791852156376941, 0, -0.2029687500000001), V(-1, 0, 0)), new PCurveEdge(new EllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(0, 0, -1), V(-0.171046939653212, 2.0947208715025655e-17, 0.9852628808776216), 0, 2.969700482947388, null, V(-1, 1.2246467991473532e-16, 0), V(0.9852628808776216, -1.2065990333854792e-16, 0.171046939653212)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.04296875, 0.07093749999999999, 0.9965548442595558), V(-0.1710469396532132, 6.938893903907228e-18, 0.9852628808776215), 1), V(-0.171046939653212, 2.0947208715025655e-17, 0.9852628808776216), V(0.04296875, 0.07093749999999999, 0.9965548442595558), 64, 0, null, V(0.002384621933525827, 0.0028002649734856517, 0.00041398320374781613), V(0.0047989723215347765, 0, -0.0002069187091194757)), new PCurveEdge(new EllipseCurve(V(0, 0.07093749999999999, 0), V(0, 0, 0.9974807622674986), V(0.9974807622674986, 0, 0), 0, 3.141592653589793), V(0.04296875, 0.07093749999999999, 0.9965548442595558), V(0.11093750000000002, 0.07093749999999999, 0.9912924604714292), 0.04309060575909555, 0.11144825166905915, null, V(0.9965548442595558, 0, -0.04296875), V(0.9912924604714292, 0, -0.11093750000000004)), new PCurveEdge(new EllipseCurve(V(0.11093750000000002, 0, 0), V(0, 0, 0.9938273849586506), V(0, 0.9938273849586506, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.07093749999999999, 0.9912924604714292), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), 0.07143883874192354, 0.11296719097084001, null, V(0, 0.9912924604714292, -0.07093749999999999), V(0, 0.9874927190198353, -0.11203125000000001)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.0029687499999999922, 0.2140625, 1.2), V(0.07406249999999998, 0.2140625, 1.2), V(0.11093750000000002, 0.18593749999999998, 1.2), V(0.11093750000000002, 0.11203125000000003, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), 1), V(0.11093750000000002, 0.11203125000000003, 0.9874927190198354), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), 100, 0, null, V(0, 0.002216485816645523, -0.00025146076711587485), V(-0.0023109369300238886, 0, -0.000007023429018985676)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1409375, 0.18703124999999998, 1.2), V(-0.08703125, 0.2059375, 1.2), V(-0.037968749999999996, 0.2140625, 1.2), V(-0.0029687499999999922, 0.2140625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(-0.1409375, 0.18703124999999998, 0.9721913045369145), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), 1), V(-0.0029687499999999922, 0.2140625, 0.9768154547391168), V(-0.1409375, 0.18703124999999998, 0.9721913045369144), 100, 0, null, V(-0.0010499999469858546, 0, -0.0000031911732431043194), V(-0.001617060512887581, -0.0005671429624910072, -0.00012531587997459948)), new PCurveEdge(new EllipseCurve(V(-0.06640184370318536, -0.0235077791500932, 0), V(0, 0, 0.997516004619599), V(-0.33289781807779234, 0.9403282523625955, 0), 0.02500215049153284, 3.1165905030982604), V(-0.1409375, 0.18703124999999998, 0.9721913045369144), V(-0.16499999999999998, 0.255, 0.9527591510974849), 0.22581372635962282, 0.3006922248899664, null, V(-0.32444628711291806, 0.9164554213903857, -0.22334333850612303), V(-0.31796125684715715, 0.8981373164189183, -0.29544573016418446)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.010937499999999989, 0.2890625, 1.2), V(-0.04093749999999999, 0.2890625, 1.2), V(-0.1, 0.2790625, 1.2), V(-0.16499999999999998, 0.255, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(0.010937499999999989, 0.2890625, 0.9572477433702835), V(-0.16499999999999998, 0.255, 0.9527591510974849), 1), V(-0.16499999999999998, 0.255, 0.9527591510974849), V(0.010937499999999989, 0.2890625, 0.9572477433702835), 100, 0, null, V(0.0019497960757584164, 0.0007217995088144136, 0.0001444829762001187), V(0.0015562475396635686, 0, -0.000017781663715540364))], []), new RotationFace(new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z), V(-0.03286362438479655, -2.6020852139652106e-18, 0.9994598452125503), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), 1), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), 29, 0, null, V(-0.003388522845144273, 0, 0.0001627386905160511), V(-0.0022900449065656465, -0.0007030314524671974, -0.00007529984920773373)), new PCurveEdge(new EllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), V(0, 0, 1), 3.108723110778215, 3.141592653589793, null, V(0.9994598452125503, -1.2239853003158589e-16, 0.032863624384797126), V(1, -1.2246467991473532e-16, 0)), new PCurveEdge(new EllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0, 0, 1), V(0.11093749999999993, 0, 0.9938273849586506), 3.141592653589793, 3.030426330354509, null, V(1, 0, 0), V(0.9938273849586506, 0, -0.11093749999999993)), new PCurveEdge(new EllipseCurve(V(0.11093750000000002, 0, 0), V(0, 0, 0.9938273849586506), V(0, 0.9938273849586506, 0), 0, 3.141592653589793), V(0.11093749999999993, 0, 0.9938273849586506), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), 0, 0.010062279327831384, null, V(0, 0.9938273849586506, 0), V(0, 0.993777073137507, -0.010000000000000007)), new PCurveEdge(new EllipseCurve(V(0, 0.010000000000000009, 0), V(0, 0, -0.9999499987499375), V(0.9999499987499375, 0, 0), 0, 3.141592653589793), V(0.11093750000000002, 0.010000000000000009, 0.9937770731375071), V(0.047968750000000004, 0.010000000000000009, 0.9987987780446257), 3.0304207486077566, 3.0936030871101416, null, V(-0.9937770731375071, 0, 0.11093749999999983), V(-0.9987987780446257, 0, 0.0479687500000002))], []), new RotationFace(new EllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.04296875, 0.07093749999999999, 1.2), V(-0.11703125, 0.07093749999999999, 1.2), V(-0.20500000000000002, 0.010000000000000009, 1.2), V(-0.20500000000000002, -0.0990625, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.17104693965321316, 2.0947208715025797e-17, 0.9852628808776215), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), 1), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), V(-0.171046939653212, 2.0947208715025655e-17, 0.9852628808776216), 36, 0, null, V(0, 0.003269989649016917, 0.0003326706415016505), V(0.0018504236016956069, 0.002172955102479456, 0.000321243497826874)), new PCurveEdge(new EllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(-0.171046939653212, 2.0947208715025655e-17, 0.9852628808776216), V(0, 0, -1), 2.969700482947388, 0, null, V(-0.9852628808776216, 1.2065990333854792e-16, -0.171046939653212), V(1, -1.2246467991473532e-16, 0)), new PCurveEdge(new EllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0, 0, -1), V(0.2029687500000001, 0, 0.9791852156376941), 0, 2.937203822794261, null, V(1, 0, 0), V(-0.9791852156376941, 0, 0.2029687500000001)), new PCurveEdge(new EllipseCurve(V(0.20296874999999998, -2.3224450485052117e-18, 0), V(0, 0, 0.9791852156376941), V(-1.1204206832959602e-17, -0.9791852156376941, 0), 2.3013070043406876e-17, 3.141592653589793), V(0.2029687500000001, 0, 0.9791852156376941), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), 2.3013070043406876e-17, 0.13412262887088133, null, V(-1.1204206832959602e-17, -0.9791852156376941, -2.3224450485052132e-18), V(-1.1103582248632314e-17, -0.9703911879325716, -0.13093749999999996)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.24203124999999998, -0.19796875, 1.2), V(0.21500000000000002, -0.18796875, 1.2), V(0.20296874999999998, -0.17296875, 1.2), V(0.20296874999999998, -0.13093749999999998, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.24203124999999998, -0.19796875, 0.9498574882827818), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), 1), V(0.20296874999999998, -0.13093749999999998, 0.9703911879325716), V(0.24203124999999998, -0.19796875, 0.9498574882827818), 100, 0, null, V(0, -0.0012607577897218358, -0.00017011744865842043), V(0.0008106439109762247, -0.00029989138903166777, -0.0002690617125763413)), new PCurveEdge(new EllipseCurve(V(0.27731540188902115, -0.09536944308866366, 0), V(0.31091053038645644, 0.9040660812655796, 0), V(0, 0, 0.9560339100680944), 1.465110231866663, 3.141592653589793), V(0.24203124999999998, -0.19796875, 0.9498574882827818), V(0.21999999999999997, -0.26203125, 0.9396486705271484), 1.684527864719763, 1.7562036900306723, null, V(-0.30890190438173115, -0.8982253957199245, -0.10849695458036647), V(-0.30558190818745745, -0.8885715060770011, -0.17624191662783179)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.12296875000000002, -0.18796875, 1.2), V(0.13906249999999998, -0.235, 1.2), V(0.17296875, -0.2559375, 1.2), V(0.21999999999999997, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.12296875000000002, -0.18796875, 0.9744467330474638), V(0.21999999999999997, -0.26203125, 0.9396486705271484), 1), V(0.21999999999999997, -0.26203125, 0.9396486705271484), V(0.12296875000000002, -0.18796875, 0.9744467330474638), 100, 0, null, V(-0.0014098883430614425, 0.00018267656272224708, 0.0003810385889743954), V(-0.0004825448774820869, 0.0014101554186612473, 0.00033291003063843776)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.034062499999999996, -0.26203125, 1.2), V(0.030000000000000027, -0.26203125, 1.2), V(0.0859375, -0.23796875, 1.2), V(0.12296875000000002, -0.18796875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.034062499999999996, -0.26203125, 0.9644580706890205), V(0.12296875000000002, -0.18796875, 0.9744467330474638), 1), V(0.12296875000000002, -0.18796875, 0.9744467330474638), V(-0.034062499999999996, -0.26203125, 0.9644580706890203), 100, 0, null, V(-0.0011108137429678576, -0.0014998329018975289, -0.00014913728313140515), V(-0.0019218307292635764, 0, -0.00006787475910567421)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.20500000000000002, -0.0990625, 1.2), V(-0.20500000000000002, -0.19703125, 1.2), V(-0.14, -0.26203125, 1.2), V(-0.034062499999999996, -0.26203125, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), V(-0.034062499999999996, -0.26203125, 0.9644580706890205), 1), V(-0.034062499999999996, -0.26203125, 0.9644580706890203), V(-0.20500000000000002, -0.0990625, 0.9737359093171772), 100, 0, null, V(-0.003177924394426102, 0, -0.00011223717544071708), V(0, 0.002937749569980493, 0.00029887037541859586))], []), new RotationFace(new EllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), [ new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.009062500000000001, -0.19296875, 1.2), V(-0.07500000000000001, -0.19296875, 1.2), V(-0.1040625, -0.16, 1.2), V(-0.1040625, -0.095, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), V(-0.1040625, -0.095, 0.9900232300778351), 1), V(-0.1040625, -0.095, 0.9900232300778351), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), 100, 0, null, V(0, -0.0019496588340396214, -0.00018708408409689773), V(0.001978121698253082, 0, 0.000018270895823513832)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(0.11093750000000002, -0.12, 1.2), V(0.0809375, -0.16593750000000002, 1.2), V(0.040937500000000016, -0.19296875, 1.2), V(-0.009062500000000001, -0.19296875, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(0.11093750000000002, -0.12, 0.9865560658643532), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), 1), V(-0.009062500000000001, -0.19296875, 0.9811630509844872), V(0.11093750000000002, -0.12, 0.9865560658643532), 100, 0, null, V(0.001499998572758444, 0, 0.000013854717676112659), V(0.0008999801458405889, 0.001378094598318401, 0.00006642278795539617)), new PCurveEdge(new EllipseCurve(V(0.11093750000000002, 2.7017833632379675e-19, 0), V(0, 0, -0.9938273849586506), V(2.4203775050019737e-18, -0.9938273849586506, 0), 1.3942163371701867e-17, 3.141592653589793), V(0.11093750000000002, -0.12, 0.9865560658643532), V(0.11093749999999993, 0, 0.9938273849586506), 3.02055199779192, 3.141592653589793, null, V(-2.4026688591808874e-18, 0.9865560658643532, 0.11999999999999994), V(-2.4203775050019737e-18, 0.9938273849586506, 1.217087525894596e-16)), new PCurveEdge(new EllipseCurve(V(0, 0, 0), V(0, 0, -1), V(1, 0, 0)), V(0.11093749999999993, 0, 0.9938273849586506), V(0, 0, 1), 3.030426330354509, 3.141592653589793, null, V(-0.9938273849586506, 0, 0.11093749999999993), V(-1, 0, 0)), new PCurveEdge(new EllipseCurve(V(0, 0, 0), V(0, 0, -1), V(-1, 1.2246467991473532e-16, 0), 0, 3.141592653589793), V(0, 0, 1), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), 3.141592653589793, 3.108723110778215, null, V(-1, 1.2246467991473532e-16, 0), V(-0.9994598452125503, 1.2239853003158589e-16, -0.032863624384797126)), new PCurveEdge(new PICurve(new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -Infinity, Infinity), new EllipsoidSurface(V3.O, V(-1, 1.2246467991473532e-16, 0), V(-1.2246467991473532e-16, -1, 0), V3.Z), V(-0.1040625, -0.095, 0.9900232300778351), V(-0.03286362438479653, 2.7235906341395923e-18, 0.9994598452125503), 1), V(-0.032863624384797126, 4.024633241122271e-18, 0.9994598452125503), V(-0.1040625, -0.095, 0.9900232300778351), 71, 0, null, V(-0.0022530279161436587, -0.0006916674357757679, -0.00007408267927848543), V(0, -0.0019215482684322, -0.00018438666887312053))], [])], false) testBRepOp(assert, d, c, d.and(c), result) }('fuz test', () => { const a = B2T.box(1, 1, 1, 'box').minus(B2T.box(1 / 3, 1 / 3, 1, 'cut')) const b = B2T.box(3, 1 / 3, 1 / 3).translate(-1).flipped() const result = new BRep([ new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), 0)), [ new StraightEdge(new L3(V(0, 0.3333333333333333, 0), V(0, 0, 1)), V(0, 0.3333333333333333, 0), V(0, 0.3333333333333333, 0.3333333333333333), 0, 0.3333333333333333), new StraightEdge(new L3(V(0, 0.3333333333333333, 0), V(0, 0, 1)), V(0, 0.3333333333333333, 0.3333333333333333), V(0, 0.3333333333333333, 1), 0.3333333333333333, 1), new StraightEdge(new L3(V(0, 0, 1), V(0, 1, 0)), V(0, 0.3333333333333333, 1), V(0, 1, 1), 0.3333333333333333, 1), new StraightEdge(new L3(V(0, 1, 0), V(0, 0, 1)), V(0, 1, 1), V(0, 1, 0), 1, 0), new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 1, 0), V(0, 0.3333333333333333, 0), 1, 0.3333333333333333)], []), new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0)), [ new StraightEdge(new L3(V(0, 0, 0.3333333333333333), V(1, 0, 0)), V(0.3333333333333333, 0, 0.3333333333333333), V(1, 0, 0.3333333333333333), 0.3333333333333333, 1), new StraightEdge(new L3(V(1, 0, 0), V(0, 0, 1)), V(1, 0, 0.3333333333333333), V(1, 0, 1), 0.3333333333333333, 1), new StraightEdge(new L3(V(1, 0, 1), V(-1, 0, 0)), V(1, 0, 1), V(0.33333333333333326, 0, 1), 0, 0.6666666666666667), new StraightEdge(new L3(V(0.3333333333333333, 0, 0), V(0, 0, -1)), V(0.33333333333333326, 0, 1), V(0.3333333333333333, 0, 0.3333333333333333), -1, -0.3333333333333333)], []), new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 0)), [ new StraightEdge(new L3(V(0, 0.3333333333333333, 0), V(-1, 0, 0)), V(1, 0.33333333333333326, 0), V(0.3333333333333333, 0.3333333333333333, 0), -1, -0.3333333333333333), new StraightEdge(new L3(V(0, 0.3333333333333333, 0), V(-1, 0, 0)), V(0.3333333333333333, 0.3333333333333333, 0), V(0, 0.3333333333333333, 0), -0.3333333333333333, 0), new StraightEdge(new L3(V(0, 0, 0), V(0, 1, 0)), V(0, 0.3333333333333333, 0), V(0, 1, 0), 0.3333333333333333, 1), new StraightEdge(new L3(V(0, 1, 0), V(1, 0, 0)), V(0, 1, 0), V(1, 1, 0), 0, 1), new StraightEdge(new L3(V(1, 1, 0), V(0, -1, 0)), V(1, 1, 0), V(1, 0.33333333333333326, 0), 0, 0.6666666666666667)], []), new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 1)), [ new StraightEdge(new L3(V(1, 0.3333333333333333, 0), V(0, 0, -1)), V(1, 0.3333333333333333, 0.3333333333333333), V(1, 0.33333333333333326, 0), -0.3333333333333333, 0), new StraightEdge(new L3(V(1, 1, 0), V(0, -1, 0)), V(1, 0.33333333333333326, 0), V(1, 1, 0), 0.6666666666666667, 0), new StraightEdge(new L3(V(1, 1, 0), V(0, 0, 1)), V(1, 1, 0), V(1, 1, 1), 0, 1), new StraightEdge(new L3(V(1, 1, 1), V(0, -1, 0)), V(1, 1, 1), V(1, 0, 1), 0, 1), new StraightEdge(new L3(V(1, 0, 0), V(0, 0, 1)), V(1, 0, 1), V(1, 0, 0.3333333333333333), 1, 0.3333333333333333), new StraightEdge(new L3(V(1, 0, 0.3333333333333333), V(0, 1, 0)), V(1, 0, 0.3333333333333333), V(1, 0.3333333333333333, 0.3333333333333333), 0, 0.3333333333333333)], []), new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), -0.3333333333333333)), [ new StraightEdge(new L3(V(0, 0.3333333333333333, 0), V(0, 0, 1)), V(0, 0.3333333333333333, 0.3333333333333333), V(0, 0.3333333333333333, 0), 0.3333333333333333, 0), new StraightEdge(new L3(V(0, 0.3333333333333333, 0), V(-1, 0, 0)), V(0, 0.3333333333333333, 0), V(0.3333333333333333, 0.3333333333333333, 0), 0, -0.3333333333333333), new StraightEdge(new L3(V(0.3333333333333333, 0.3333333333333333, 0), V(0, 0, 1)), V(0.3333333333333333, 0.3333333333333333, 0), V(0.33333333333333326, 0.3333333333333333, 0.3333333333333333), 0, 0.3333333333333333), new StraightEdge(new L3(V(0, 0.3333333333333333, 0.3333333333333333), V(1, 0, 0)), V(0.3333333333333333, 0.3333333333333333, 0.3333333333333333), V(0, 0.3333333333333333, 0.3333333333333333), 0.3333333333333333, 0)], []), new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), -0.3333333333333333)), [ new StraightEdge(new L3(V(0, 0.3333333333333333, 0.3333333333333333), V(1, 0, 0)), V(0, 0.3333333333333333, 0.3333333333333333), V(0.3333333333333333, 0.3333333333333333, 0.3333333333333333), 0, 0.3333333333333333), new StraightEdge(new L3(V(0.3333333333333333, 0.3333333333333333, 0), V(0, 0, 1)), V(0.3333333333333333, 0.3333333333333333, 0.3333333333333333), V(0.3333333333333333, 0.3333333333333333, 1), 0.3333333333333333, 1), new StraightEdge(new L3(V(0, 0.3333333333333333, 1), V(1, 0, 0)), V(0.3333333333333333, 0.3333333333333333, 1), V(0, 0.3333333333333333, 1), 0.3333333333333333, 0), new StraightEdge(new L3(V(0, 0.3333333333333333, 0), V(0, 0, 1)), V(0, 0.3333333333333333, 1), V(0, 0.3333333333333333, 0.3333333333333333), 1, 0.3333333333333333)], []), new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), -0.3333333333333333)), [ new StraightEdge(new L3(V(0.3333333333333333, 0, 0.3333333333333333), V(0, -1, 0)), V(0.3333333333333333, 0.3333333333333333, 0.3333333333333333), V(0.3333333333333333, 0, 0.3333333333333333), -0.3333333333333333, 0), new StraightEdge(new L3(V(0.3333333333333333, 0, 0), V(0, 0, -1)), V(0.3333333333333333, 0, 0.3333333333333333), V(0.33333333333333326, 0, 1), -0.3333333333333333, -1), new StraightEdge(new L3(V(0.3333333333333333, 0, 1), V(0, -1, 0)), V(0.33333333333333326, 0, 1), V(0.3333333333333333, 0.3333333333333333, 1), 0, -0.3333333333333333), new StraightEdge(new L3(V(0.3333333333333333, 0.3333333333333333, 0), V(0, 0, 1)), V(0.3333333333333333, 0.3333333333333333, 1), V(0.3333333333333333, 0.3333333333333333, 0.3333333333333333), 1, 0.3333333333333333)], []), new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 1)), [ new StraightEdge(new L3(V(0, 0.3333333333333333, 1), V(1, 0, 0)), V(0, 0.3333333333333333, 1), V(0.3333333333333333, 0.3333333333333333, 1), 0, 0.3333333333333333), new StraightEdge(new L3(V(0.3333333333333333, 0, 1), V(0, -1, 0)), V(0.3333333333333333, 0.3333333333333333, 1), V(0.33333333333333326, 0, 1), -0.3333333333333333, 0), new StraightEdge(new L3(V(1, 0, 1), V(-1, 0, 0)), V(0.33333333333333326, 0, 1), V(1, 0, 1), 0.6666666666666667, 0), new StraightEdge(new L3(V(1, 1, 1), V(0, -1, 0)), V(1, 0, 1), V(1, 1, 1), 1, 0), new StraightEdge(new L3(V(0, 1, 1), V(1, 0, 0)), V(1, 1, 1), V(0, 1, 1), 1, 0), new StraightEdge(new L3(V(0, 0, 1), V(0, 1, 0)), V(0, 1, 1), V(0, 0.3333333333333333, 1), 1, 0.3333333333333333)], []), new PlaneFace(new PlaneSurface(new P3(V(0, 1, 0), 1)), [ new StraightEdge(new L3(V(0, 1, 0), V(1, 0, 0)), V(1, 1, 0), V(0, 1, 0), 1, 0), new StraightEdge(new L3(V(0, 1, 0), V(0, 0, 1)), V(0, 1, 0), V(0, 1, 1), 0, 1), new StraightEdge(new L3(V(0, 1, 1), V(1, 0, 0)), V(0, 1, 1), V(1, 1, 1), 0, 1), new StraightEdge(new L3(V(1, 1, 0), V(0, 0, 1)), V(1, 1, 1), V(1, 1, 0), 1, 0)], []), new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), -0.3333333333333333)), [ new StraightEdge(new L3(V(0, 0.3333333333333333, 0), V(-1, 0, 0)), V(0.3333333333333333, 0.3333333333333333, 0), V(1, 0.33333333333333326, 0), -0.3333333333333333, -1), new StraightEdge(new L3(V(1, 0.3333333333333333, 0), V(0, 0, -1)), V(1, 0.33333333333333326, 0), V(1, 0.3333333333333333, 0.3333333333333333), 0, -0.3333333333333333), new StraightEdge(new L3(V(-1, 0.3333333333333333, 0.3333333333333333), V(1, 0, 0)), V(1, 0.3333333333333333, 0.3333333333333333), V(0.33333333333333326, 0.3333333333333333, 0.3333333333333333), 2, 1.3333333333333333), new StraightEdge(new L3(V(0.3333333333333333, 0.3333333333333333, 0), V(0, 0, 1)), V(0.33333333333333326, 0.3333333333333333, 0.3333333333333333), V(0.3333333333333333, 0.3333333333333333, 0), 0.3333333333333333, 0)], []), new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), -0.3333333333333333)), [ new StraightEdge(new L3(V(0, 0, 0.3333333333333333), V(1, 0, 0)), V(1, 0, 0.3333333333333333), V(0.3333333333333333, 0, 0.3333333333333333), 1, 0.3333333333333333), new StraightEdge(new L3(V(0.3333333333333333, 0, 0.3333333333333333), V(0, -1, 0)), V(0.3333333333333333, 0, 0.3333333333333333), V(0.3333333333333333, 0.3333333333333333, 0.3333333333333333), 0, -0.3333333333333333), new StraightEdge(new L3(V(-1, 0.3333333333333333, 0.3333333333333333), V(1, 0, 0)), V(0.33333333333333326, 0.3333333333333333, 0.3333333333333333), V(1, 0.3333333333333333, 0.3333333333333333), 1.3333333333333333, 2), new StraightEdge(new L3(V(1, 0, 0.3333333333333333), V(0, 1, 0)), V(1, 0.3333333333333333, 0.3333333333333333), V(1, 0, 0.3333333333333333), 0.3333333333333333, 0)], [])], false) testBRepAnd(assert, a, b, result) })
  })

  //test('BRep.withMergedFaces', () => {
  //    const box = B2T.box(5, 5, 5)
  //    const boxToMerge = new BRep(box.faces.filter((face:PlaneFace) => face.surface.plane.normal1.x != 1).concat(
  //        box.translate(5, 0, 0).faces.filter((face:PlaneFace) => face.surface.plane.normal1.x != -1)
  //    ), false)
  //
  //    assert.equal(boxToMerge.faces.length, 10)
  //    const boxMerged = boxToMerge.withMergedFaces()
  //    testBRepOp(assert, boxToMerge, BRep.EMPTY, boxMerged, B2T.box(10, 5, 5))
  //}
})
