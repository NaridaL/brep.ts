import {
  outputLink,
  suiteSurface,
  surfaceVolumeAndAreaTests,
  testContainsCurve,
  testISCurves,
  testLoopCCW,
  testLoopContainsPoint,
  testSurfaceTransform,
} from "./manager"

import { DEG, M4, V, V3 } from "ts3dutils"
import {
  B2T,
  ConicSurface,
  Edge,
  EllipseCurve,
  EllipsoidSurface,
  HyperbolaCurve,
  L3,
  P3,
  ParabolaCurve,
  PCurveEdge,
  PointVsFace,
  RotationFace,
  StraightEdge,
} from "../src"
import { PI } from "../src/math"

describe("ConicSurface", () => {
  const UCS = ConicSurface.UNIT

  describe("UNIT", () => suiteSurface(ConicSurface.UNIT))
  describe("UNIT.scale(2, 2, 1)", () =>
    suiteSurface(ConicSurface.UNIT.scale(2, 2, 1)))
  describe("weird", () =>
    suiteSurface(
      new ConicSurface(
        V(2, 0.2, 1.1),
        V(0, 0.6, 0),
        V(0, 0, -2.4),
        V(-12, 0, 0),
      ),
    ))

  const testFace = new RotationFace(
    new ConicSurface(V3.Z, V(-1, 0, 0), V3.Y, V(0, 0, -1)),
    [
      new PCurveEdge(
        new HyperbolaCurve(
          V(0.10792279653395696, 0.0905579787672639, 1),
          V(0, 0, -0.1408832052805518),
          V(0.0905579787672639, -0.10792279653395696, 0),
          -7,
          7,
        ),
        V(-0.4634542598189221, 0.7714986384017117, 0.1),
        V(0.1, 0.1, 0.8585786437626904),
        -2.5414277085137025,
        -0.08737743596203365,
        undefined,
        V(0.5785088487178852, -0.6894399988070802, 0.8889049006895381),
        V(0.09090389553440875, -0.10833504408394042, 0.012325683343243887),
        "genseg17",
      ),
      new PCurveEdge(
        new HyperbolaCurve(
          V(-0.00792279653395693, 0.009442021232736126, 1),
          V(0, 0, -0.012325683343243885),
          V(0.009442021232736126, 0.00792279653395693, 0),
          -7,
          7,
        ),
        V(0.1, 0.1, 0.8585786437626904),
        V(0.6814525440390486, 0.5878966152502569, 0.1),
        3.131301331471644,
        4.983809888872043,
        undefined,
        V(0.10833504408394039, 0.09090389553440875, -0.14088320528055173),
        V(0.6894399988070802, 0.5785088487178854, -0.8999155946699233),
        "genseg18",
      ),
      new PCurveEdge(
        new EllipseCurve(
          V(0, 0, 0.09999999999999998),
          V(0.9, 0, 0),
          V(0, 0.9, 0),
          0,
          3.141592653589793,
        ),
        V(0.6814525440390486, 0.5878966152502569, 0.1),
        V(-0.4634542598189221, 0.7714986384017117, 0.1),
        0.7118273326574678,
        2.1117446875459924,
        undefined,
        V(-0.5878966152502568, 0.6814525440390486, 0),
        V(-0.7714986384017115, -0.4634542598189224, 0),
        "genseg19",
      ),
    ],
    [],
  )
  describe("testFace surface", () =>
    suiteSurface(testFace.surface as ConicSurface))
  describe("testFace", () => surfaceVolumeAndAreaTests(testFace))
  describe("testFace.scale(2)", () =>
    surfaceVolumeAndAreaTests(testFace.scale(2)))
  describe("testFace.shearX(2, 2)", () =>
    surfaceVolumeAndAreaTests(testFace.shearX(2, 2)))
  describe("testFace.foo()", () => surfaceVolumeAndAreaTests(testFace.foo()))
  test("testLoopCCW", () => {
    const surface = new ConicSurface(
      V(0, 0, 53.51411369448604),
      V(198.46477746372744, 0, 0),
      V(0, 198.46477746372744, 0),
      V(0, 0, 191.42941531213293),
    ).scale(1 / 200)
    const loop = [
      new StraightEdge(
        new L3(
          V(131.35224103228387, 0, 180.2100595549249),
          V(0.7197488536413841, 0, 0.6942345336281635),
        ),
        V(131.35224103228387, 0, 180.2100595549249),
        V(198.46477746372744, 0, 244.94352900661897),
        0,
        93.24438113642698,
      ),
      new PCurveEdge(
        new EllipseCurve(
          V(0, 0, 244.94352900661897),
          V(198.46477746372744, 0, 0),
          V(0, 198.46477746372744, 0),
          0,
          3.141592653589793,
        ),
        V(198.46477746372744, 0, 244.94352900661897),
        V(-198.46477746372744, 2.4304925446444556e-14, 244.94352900661897),
        0,
        3.141592653589793,
        null,
        V(0, 198.46477746372744, 0),
        V(-2.4304925446444556e-14, -198.46477746372744, 0),
      ),
      new StraightEdge(
        new L3(
          V(-131.35224103228387, 1.6086010154101807e-14, 180.2100595549249),
          V(-0.7197488536413841, 8.814381298018978e-17, 0.6942345336281635),
        ),
        V(-198.46477746372744, 2.4304925446444556e-14, 244.94352900661897),
        V(-131.35224103228387, 1.6086010154101807e-14, 180.2100595549249),
        93.24438113642698,
        0,
      ),
      new PCurveEdge(
        new EllipseCurve(
          V(0, 0, 180.2100595549249),
          V(131.35224103228387, 0, 0),
          V(0, 131.35224103228387, 0),
          0,
          3.141592653589793,
        ),
        V(-131.35224103228387, 1.6086010154101807e-14, 180.2100595549249),
        V(131.35224103228387, 0, 180.2100595549249),
        3.141592653589793,
        0,
        null,
        V(1.6086010154101807e-14, 131.35224103228387, 0),
        V(0, -131.35224103228387, 0),
      ),
    ].map((e) => e.scale(1 / 200))
    testLoopCCW(surface, Edge.reversePath(loop))
    testLoopCCW(surface.flipped(), loop)
  })
  test("isCoplanarTo", () => {
    const unitCone = ConicSurface.UNIT
    expect(unitCone.matrix.isIdentity()).toBeTruthy()
    expect(unitCone.pUVFunc()(0, 3)).toBeLike(V(3, 0, 3))
    const ellipseAtZ3 = EllipseCurve.UNIT.scale(3, 3, 3).translate(0, 0, 3)
    const planeAtZ3 = P3.XY.translate(0, 0, 3)
    const issAtZ3 = unitCone.isCurvesWithPlane(planeAtZ3)
    expect(issAtZ3.length).toBe(1)
    expect(ellipseAtZ3.isColinearTo(issAtZ3[0])).toBeTruthy()
    expect(unitCone.containsEllipse(ellipseAtZ3)).toBeTruthy()

    const scaledUnit = ConicSurface.UNIT.scale(2, 2, 1)
    expect(scaledUnit.isCoplanarTo(unitCone)).toBeFalsy()
    expect(unitCone.isCoplanarTo(scaledUnit)).toBeFalsy()
    const ell1 = unitCone.isCurvesWithPlane(
      new P3(V(2, 3, 10).unit(), 10),
    )[0] as EllipseCurve
    expect(unitCone.containsEllipse(ell1)).toBeTruthy()
    const ell2 = unitCone.isCurvesWithPlane(
      new P3(V(1, 1, 2).unit(), 4),
    )[0] as EllipseCurve
    const ell1Cone = ConicSurface.atApexThroughEllipse(V3.O, ell1)
    const ell2Cone = ConicSurface.atApexThroughEllipse(V3.O, ell2)
    console.log(ell1Cone)
    expect(unitCone.isCoplanarTo(ell1Cone)).toBeTruthy()
    expect(unitCone.isCoplanarTo(ell2Cone)).toBeTruthy()
    expect(ell1Cone.isCoplanarTo(ell2Cone)).toBeTruthy()
    expect(ell2Cone.isCoplanarTo(ell1Cone)).toBeTruthy()
    expect(ell1Cone.foo().isCoplanarTo(ell2Cone.foo())).toBeTruthy()
  })
  test("isCurvesWithPlane", () => {
    testISCurves(UCS, new P3(V(1, 1, 2).unit(), 4), 2)
    testISCurves(UCS, P3.XY.translate(0, 0, 3), 1)
    testISCurves(UCS, P3.XY.translate(0, 0, 3).flipped(), 1)
  })
  test("isCurvesWithPlane 2", () => {
    testISCurves(ConicSurface.UNIT, P3.ZX, 2)
    testISCurves(ConicSurface.UNIT, P3.ZX.flipped(), 2)
    testISCurves(ConicSurface.UNIT, P3.YZ, 2)
    testISCurves(ConicSurface.UNIT, P3.YZ.flipped(), 2)
    testISCurves(ConicSurface.UNIT, new P3(V(1, 0, 1).unit(), 4), 1)
  })
  test("isCurvesWithPlane hyperbolas", () => {
    const plane = new P3(V(2, 0, -1).unit(), 1)
    testISCurves(UCS, plane, 1)
    testISCurves(UCS, plane.flipped(), 1)
  })
  test("isCurvesWithEllipsoid", () => {
    const a = ConicSurface.UNIT.scale(0.05, 0.2)
      .rotateZ(90 * DEG)
      .rotateY(-90 * DEG)
      .translate(2, 0.2, 1.1)
      .flipped()
    const cone = new ConicSurface(
      V(2, 0.2, 1.1),
      V(0, 0.6, 0),
      V(0, 0, -2.4),
      V(-12, 0, 0),
    )
    const sphere = new EllipsoidSurface(V3.O, V3.X, V3.Y, V(0, 0, -1))
    const b = EllipsoidSurface.UNIT
    testISCurves(a, b, 2)
    testISCurves(cone, sphere, 2)
  })
  test("isCurvesWithEllipsoid 2", () => {
    const cone = new ConicSurface(
      V(2, 0.2, 0.7),
      V(2.2496396739927868e-33, 0.6000000000000001, 3.67394039744206e-17),
      V(-1.469576158976824e-16, 1.469576158976824e-16, -2.4000000000000004),
      V(-12, 0, 7.347880794884119e-16),
    )
    const sphere = new EllipsoidSurface(V3.O, V3.X, V3.Y, V(0, 0, -1))

    testISCurves(cone, sphere, 2)
  })
  test("containsParabola", () => {
    const pb = UCS.isCurvesWithPlane(
      new P3(V(1, 0, 1).unit(), 4),
    )[0] as ParabolaCurve
    expect(pb instanceof ParabolaCurve).toBeTruthy()
    testContainsCurve(UCS, pb)
  })
  test("containsParabola 2", () => {
    const c2 = UCS.shearX(2, 3)
    const pb2 = c2.isCurvesWithPlane(
      new P3(V(1, 0, 1).unit(), 4).shearX(2, 3),
    )[0] as ParabolaCurve
    expect(pb2 instanceof ParabolaCurve).toBeTruthy()
    testContainsCurve(c2, pb2)
  })
  test("containsHyperbola", () => {
    const s = new ConicSurface(
      V(-242.1625189124994, 38.960257711878945, 0),
      V(197.87979681325515, -15.226749714620981, 2.4304925446444556e-14),
      V(2.4233285978328154e-14, -1.8647390299428456e-15, -198.46477746372744),
      V(14.686977871964286, 190.86517159433123, 0),
    )
    const c = new HyperbolaCurve(
      V(-242.16251891249937, 38.960257711878945, -100.00000000000003),
      V(7.400294429901329, 96.17080372320217, 0),
      V(-99.70524711843181, 7.672268051394617, -1.8369701987210304e-14),
    )

    testContainsCurve(s, c)
  })
  test("containsHyperbola 2", () => {
    const pb = UCS.isCurvesWithPlane(new P3(V3.Y, 2))[0] as HyperbolaCurve
    testContainsCurve(UCS, pb)
    testContainsCurve(UCS, pb.translate(0, 2, 0), false, ".translate(0, 2, 0)")
    testContainsCurve(UCS, pb.scale(1, 2, 1), false, ".scale(1, 2, 1)")
    testContainsCurve(
      UCS,
      pb.rotate(pb.center.plus(pb.f1), pb.f2, -20 * DEG),
      false,
      ".rotate(pb.center.plus(pb.f1), pb.f2, -20 * DEG)",
    )
  })
  test("containsHyperbola 3", () => {
    const pb = UCS.isCurvesWithPlane(
      new P3(V(0, 1, -0.1).unit(), 2),
    )[0] as HyperbolaCurve
    console.log(pb)
    testContainsCurve(UCS, pb)
    testContainsCurve(UCS, pb.translate(0, 2, 0), false, ".translate(0, 2, 0)")
    testContainsCurve(UCS, pb.scale(1, 2, 1), false, ".scale(1, 2, 1)")
    testContainsCurve(
      UCS,
      pb.rotate(pb.center.plus(pb.f1), pb.f2, -20 * DEG),
      false,
      ".rotate(pb.center.plus(pb.f1), pb.f2, -20 * DEG)",
    )
  })
  test("containsPoint", () => {
    const face = B2T.cone(1, 1, PI).faces.find(
      (face) => face.surface instanceof ConicSurface,
    )
    testLoopContainsPoint(
      face.surface,
      face.contour,
      V(-0.2, 0, 0.2),
      PointVsFace.ON_EDGE,
    )
  })

  test("atApexThroughEllipse", () => {
    const apex = V(0, 0, 3)
    const ellipse = new EllipseCurve(
      V(0.8047378541243649, 0, 0.3333333333333333),
      V(2.414213562373095, 0, 0),
      V(0, 2.414213562373095, 0),
      0,
      3.141592653589793,
    )
    const cone = ConicSurface.atApexThroughEllipse(apex, ellipse)
    expect(cone.containsPoint(ellipse.at(0))).toBeTruthy()
    outputLink({ mesh: cone.sce + ".toMesh()", drPs: [ellipse.at(0)] })
  })

  test("transform4", () => {
    const c = new ConicSurface(V3.O, V3.X, V3.Y, V3.Z, 0, PI, 1, 2)
      .rotateZ(20 * DEG)
      .scale(1, 0.8, -2)
      .rotateY(10 * DEG)
      .translate(1, 0, -1)
    const m = M4.perspective(45, 1, 2, 6)
    testSurfaceTransform(c, m)
  })

  test("transform4 2", () => {
    const c = new ConicSurface(V3.O, V3.X, V3.Y, V3.Z, 0, PI, 1, 2)
      .scale(1, 1, -2)
      .rotateY(10 * DEG)
      .translate(1, 0, 0)
    console.log(c.sce)
    const m = M4.perspective(45, 1, 2, 6)
    testSurfaceTransform(c, m)
  })
})
