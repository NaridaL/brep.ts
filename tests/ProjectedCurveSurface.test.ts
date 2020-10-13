import {
  inDifferentSystems,
  outputLink,
  suite,
  suiteSurface,
  surfaceVolumeAndAreaTests,
  test,
  testISCurves,
  testISTs,
  testLoopCCW,
  testSurfaceTransform,
} from "./manager"

import { DEG, M4, V, V3 } from "ts3dutils"
import {
  B2T,
  BezierCurve,
  CylinderSurface,
  Edge,
  EllipseCurve,
  EllipsoidSurface,
  Face,
  L3,
  P3,
  PCurveEdge,
  PICurve,
  PlaneSurface,
  ProjectedCurveSurface,
  RotationFace,
  StraightEdge,
} from ".."

suite("ProjectedCurveSurface", () => {
  const baseCurve = BezierCurve.graphXY(2, -3, -3, 2, 0, 2)
  const testSurface = new ProjectedCurveSurface(
    baseCurve,
    V3.Z,
    undefined,
    undefined,
    0,
    2,
  )

  const edge = PCurveedgeForCurveAndTs(baseCurve, 0, 2)
  const edges = [
    edge,
    StraightEdge.throughPoints(
      baseCurve.at(2),
      baseCurve.at(2).plus(V(0, 0, 2)),
    ),
    edge.flipped().translate(0, 0, 2),
    StraightEdge.throughPoints(
      baseCurve.at(0).plus(V(0, 0, 2)),
      baseCurve.at(0),
    ),
  ]
  const testFace = new RotationFace(testSurface, edges)
  suite("projectedBezierSurface", () => suiteSurface(testSurface))
  suite("projectedBezierSurface.shearX(2, 2)", () =>
    suiteSurface(testSurface.shearX(2, 2)),
  )
  suite("projectedBezierSurface.foo()", () => suiteSurface(testSurface.foo()))
  const curve = BezierCurve.graphXY(2, -3, -3, 2)
  const surface = new ProjectedCurveSurface(curve, V3.Z)

  //test('(face w/ dir=V3.Z).rotateY(20 * DEG) ps test', assert => testSurface(assert, pcsFace.rotateY(20 *
  // DEG).surface as ParametricSurface))
  suite("face w/ dir=V3.Z", () => surfaceVolumeAndAreaTests(testFace))
  suite("(face w/ dir=V3.Z).rotateY(90 * DEG).translate(0, 0, 1)", () =>
    surfaceVolumeAndAreaTests(testFace.rotateY(90 * DEG).translate(0, 0, 1)),
  )
  suite("(face w/ dir=V3.Z).shearX(2, 2)", () =>
    surfaceVolumeAndAreaTests(testFace.shearX(2, 2)),
  )
  suite("(face w/ dir=V3.Z).foo()", () => {
    console.log("pUV(0,0)", testSurface.foo().pUV(0, 0))
    const { u, v } = testSurface.foo().uvP(V(-0.59566, 12.37799, 2.98851))
    console.log(u, v)
    console.log(testSurface.foo().pUV(u, v))
    surfaceVolumeAndAreaTests(testFace.foo())
  })
  suite("Face line intersection test", () =>
    inDifferentSystems((assert, m4) => {
      const line = new L3(V3.Z, V3.X).transform(m4)
      const d = testFace.transform(m4).intersectsLine(line)
      assert.ok(d)
    }),
  )
  const pcs1 = new ProjectedCurveSurface(
    BezierCurve.EX2D,
    V3.Z,
    undefined,
    undefined,
    -1,
    4,
  )
  suite("pcs1", () => suiteSurface(pcs1))
  test("isTsWithSurface", (assert) => {
    const pcs = new ProjectedCurveSurface(
      BezierCurve.EX2D,
      V3.Z,
      undefined,
      undefined,
      -2,
      2,
    )
      .scale(0.2, 0.2, 1)
      .rotateX(-90 * DEG)
    const ses = EllipsoidSurface.UNIT
    const pic = ses.isCurvesWithSurface(pcs)[0]
    testISTs(assert, pic, new PlaneSurface(P3.XY), 1)
  })

  // create a pcs face which includes a PICurve
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
  const loop = [
    new PCurveEdge(
      PICurve.forParametricStartEnd(
        new ProjectedCurveSurface(
          new BezierCurve(
            V(0, 2, 0),
            V(0.3333333333333333, 1, 0),
            V(0.6666666666666666, -1, 0),
            V(1, -2, 0),
            -0.1,
            1.1,
          ),
          V(0, 0, -1),
          -0.1,
          1.1,
          -1,
          0,
        ),
        new CylinderSurface(
          new EllipseCurve(
            V(0, 0, 0.5),
            V(1.2246467991473533e-17, 0, -0.2),
            V(0, 0.2, 0),
            0,
            3.141592653589793,
          ),
          V(1, 0, 6.123233995736766e-17),
          0,
          3.141592653589793,
          0,
          2,
        ),
        V(0.5020636899504077, -0.30021571809025155, 0),
        V(0.50440179447819, -0.6990167000940748, 0),
        0.05,
        V(-0.048940803585290286, 0.010237076947353528, 0),
        0.21324043802451342,
        10.421498564770445,
      ),
      V(0.49999999999999994, 0, 0.30000000000000004),
      V(0.4579408460333606, 0.1891173898820845, 0.4349260970573627),
      0.21324043802451342,
      3.9862958940211684,
      undefined,
      V(-0.04683824539832595, 0.21076813481893145, -0.0014356938902703453),
      V(-0.003821083373315509, 0.017154318833024252, 0.04985377941394385),
      "genseg2",
    ),
    new PCurveEdge(
      PICurve.forParametricStartEnd(
        new ProjectedCurveSurface(
          new BezierCurve(
            V(0, 2, 0),
            V(0.3333333333333333, 1, 0),
            V(0.6666666666666666, -1, 0),
            V(1, -2, 0),
            -0.1,
            1.1,
          ),
          V(0, 0, -1),
          0,
          1,
          -1,
          0,
        ),
        new CylinderSurface(
          new EllipseCurve(
            V(0, 0, 0.5),
            V(1.2246467991473533e-17, 0, -0.2),
            V(0, 0.2, 0),
            0,
            3.141592653589793,
          ),
          V(1, 0, 6.123233995736766e-17),
          0,
          3.141592653589793,
          0,
          2,
        ),
        V(0.45794084603336055, -0.4349260970573627, 0),
        V(0.543312790740714, -0.5455452074045042, 0),
        0.05,
        V(-0.003821083373315509, -0.04985377941394385, 0),
        0,
        6.3461388249955215,
      ),
      V(0.4579408460333606, 0.1891173898820845, 0.4349260970573627),
      V(0.49999999999999994, 2.4492935982947065e-17, 0.7),
      0,
      6.3461388249955215,
      undefined,
      V(-0.003821083373315509, 0.017154318833024252, 0.04985377941394385),
      V(0.04658715677301388, -0.209638049867369, 0.0007477447876345841),
      "genseg3",
    ),
    new PCurveEdge(
      PICurve.forParametricStartEnd(
        new ProjectedCurveSurface(
          new BezierCurve(
            V(0, 2, 0),
            V(0.3333333333333333, 1, 0),
            V(0.6666666666666666, -1, 0),
            V(1, -2, 0),
            -0.1,
            1.1,
          ),
          V(0, 0, -1),
          0,
          1,
          -1,
          0,
        ),
        new CylinderSurface(
          new EllipseCurve(
            V(0, 0, 0.5),
            V(-1.2246467991473533e-17, 2.4492935982947065e-17, 0.2),
            V(-1.4997597826618578e-33, -0.2, 2.4492935982947065e-17),
            0,
            3.141592653589793,
          ),
          V(1, 0, 6.123233995736766e-17),
          0,
          3.141592653589793,
          0,
          2,
        ),
        V(0.4839777423387698, -0.686554960163734, 0),
        V(0.4944895751568088, -0.3015431287112728, 0),
        0.05,
        V(0.024931684232969407, -0.04334064052719462, 0),
        1.1569694674108177,
        11.336604127893224,
      ),
      V(0.49999999999999994, 2.4492935982947065e-17, 0.7),
      V(0.543312790740714, -0.19474504892931377, 0.5455452074045042),
      1.1569694674108177,
      5.303950261557475,
      undefined,
      V(0.04658715677301388, -0.209638049867369, 0.0007477447876345841),
      V(0.002601553486113264, -0.011677707635158832, -0.04993227332556462),
      "genseg5",
    ),
    new PCurveEdge(
      PICurve.forParametricStartEnd(
        new ProjectedCurveSurface(
          new BezierCurve(
            V(0, 2, 0),
            V(0.3333333333333333, 1, 0),
            V(0.6666666666666666, -1, 0),
            V(1, -2, 0),
            -0.1,
            1.1,
          ),
          V(0, 0, -1),
          0,
          1,
          -1,
          0,
        ),
        new CylinderSurface(
          new EllipseCurve(
            V(0, 0, 0.5),
            V(-1.2246467991473533e-17, 2.4492935982947065e-17, 0.2),
            V(-1.4997597826618578e-33, -0.2, 2.4492935982947065e-17),
            0,
            3.141592653589793,
          ),
          V(1, 0, 6.123233995736766e-17),
          0,
          3.141592653589793,
          0,
          2,
        ),
        V(0.543312790740714, -0.5455452074045042, 0),
        V(0.45794084603336055, -0.43492609705736274, 0),
        0.05,
        V(0.002601553486113264, 0.04993227332556462, 0),
        0,
        6.224025149750129,
      ),
      V(0.543312790740714, -0.19474504892931377, 0.5455452074045042),
      V(0.49999999999999994, 0, 0.30000000000000004),
      0,
      6.224025149750129,
      undefined,
      V(0.002601553486113264, -0.011677707635158832, -0.04993227332556462),
      V(-0.04683824539832592, 0.2107681348189313, -0.0014356938902701979),
      "genseg4",
    ),
  ]
  const piCurveFace = Face.create(
    a.faces.find((f) => f.surface instanceof ProjectedCurveSurface).surface,
    loop,
  )
  test("ccw for round face", (assert) => {
    testLoopCCW(
      assert,
      a.faces.find((f) => f.surface instanceof ProjectedCurveSurface).surface,
      loop,
    )
  })
  test("isCurvesWithCylinderSurface", (assert) => {
    testISCurves(
      assert,
      a.faces.find((f) => f.surface instanceof ProjectedCurveSurface).surface,
      b.faces.find((f) => f.surface instanceof CylinderSurface).surface,
      1,
    )
  })
  test("isCurvesWithSurface w/ PCS", (assert) => {
    const pcs = new ProjectedCurveSurface(
      BezierCurve.EX2D,
      V3.Z,
      -1.1,
      2.1,
      0,
      4,
    )
    const pcs2 = pcs.mirrorY().rotateX(80 * DEG)
    testISCurves(assert, pcs, pcs2, 2)
  })

  // TODO
  //suite('face w/ PICurve', () => surfaceVolumeAndAreaTests(piCurveFace))
  //suite('(face w/ PICurve).shearX(2, 2)', () => surfaceVolumeAndAreaTests(piCurveFace.shearX(2, 2)))
  //suite('(face w/ PICurve).foo()', () => surfaceVolumeAndAreaTests(piCurveFace.foo()))

  test("Face containsPoint", (assert) => {
    const face = new RotationFace(
      new ProjectedCurveSurface(
        new BezierCurve(
          V(142.87578921496748, -191.46078243076332, 0),
          V(161.78547089700214, -252.13248349581008, 0),
          V(284.63214994898954, -163.59789158697575, 0),
          V(372.40411211189405, -210.3992206435476, 0),
          -3,
          4,
        ),
        V(0, 0, 1),
        -3,
        4,
        -100,
        100,
      ),
      [
        PCurveedgeForCurveAndTs(
          new BezierCurve(
            V(142.87578921496748, -191.46078243076332, 0),
            V(161.78547089700214, -252.13248349581008, 0),
            V(284.63214994898954, -163.59789158697575, 0),
            V(372.40411211189405, -210.3992206435476, 0),
          ),
          1,
          0,
        ),
        StraightEdge.throughPoints(
          V(142.87578921496748, -191.46078243076332, 0),
          V(142.87578921496748, -191.46078243076332, -100),
        ),
        PCurveedgeForCurveAndTs(
          new BezierCurve(
            V(142.87578921496748, -191.46078243076332, -100),
            V(161.78547089700214, -252.13248349581008, -100),
            V(284.63214994898954, -163.59789158697575, -100),
            V(372.40411211189405, -210.3992206435476, -100),
          ),
          0,
          1,
        ),
        StraightEdge.throughPoints(
          V(372.40411211189405, -210.3992206435476, -100),
          V(372.40411211189405, -210.3992206435476, 0),
        ),
      ],
      [],
    )
    const line = new L3(
      V(1241.5987, -1214.1894, 38.9886),
      V(-0.6705, 0.7386, -0.0696).unit(),
    )
    testISTs(assert, line, face.surface, 3)
  })

  test("transform4", (assert) => {
    const s = CylinderSurface.UNIT.translate(1, 0, -3)
    const m = M4.perspective(45, 1, 2, 4)
    testSurfaceTransform(assert, s, m)
  })

  test("transform4 2", (assert) => {
    const s = CylinderSurface.UNIT.translate(1, 0, -3).flipped()
    const m = M4.perspective(45, 1, 2, 4)
    testSurfaceTransform(assert, s, m)
  })

  test("transform4 3", (assert) => {
    const s = CylinderSurface.UNIT.rotateX(90 * DEG).translate(1, 0, -3)
    const m = M4.perspective(45, 1, 2, 4)
    testSurfaceTransform(assert, s, m)
  })

  test("boxprespected", (assert) => {
    const punch = B2T.cylinder(0.5, 2).translate(1, 1, 0)
    const box = B2T.box(2, 2, 2)
      .minus(punch)
      .rotateY(10 * DEG)
      .translate(-0, -0, -5)
    const m = M4.perspective(45, 1, 2, 5)
    outputLink(assert, {
      a: box,
      b: box.transform4(m).translate(10),
    })
  })
})
