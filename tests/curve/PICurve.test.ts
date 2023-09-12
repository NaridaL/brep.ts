import { testCurve, testISCurves, testISTs, testPointT } from "../manager"
import dedent from "dedent"
import { DEG, M4, V, V3 } from "ts3dutils"
import {
  BezierCurve,
  CylinderSurface,
  EllipseCurve,
  EllipsoidSurface,
  P3,
  PCurveEdge,
  PICurve,
  PlaneSurface,
  ProjectedCurveSurface,
  RotatedCurveSurface,
} from "brepts"
import { doc } from "prettier"

describe("PICurve", () => {
  describe("pointT", () => {
    test("point not on curve => undefined", () => {
      const piCurve = PICurve.forStartEnd(
        new ProjectedCurveSurface(
          new BezierCurve(
            V(0.20296874999999998, 0.11703124999999998, 1.2),
            V(0.20296874999999998, 0.2240625, 1.2),
            V(0.14500000000000002, 0.2890625, 1.2),
            V(0.010937499999999989, 0.2890625, 1.2),
            0,
            1,
          ),
          V(0, 0, -1),
          0,
          1,
          -100,
          100,
        ),
        new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z),
        V(0.010937499999999989, 0.2890625, -0.9572477433702835),
        V(0.20296874999999998, 0.11703124999999998, -0.9721663299286162),
        0.02,
      )
      const p = V(0.010937499999999989, 0.2890625, 0.9572477433702835)
      testPointT(piCurve, p, NaN)
    })

    test("piCurve.points[8]", () => {
      const piCurve = PICurve.forParametricStartEnd(
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
      )
      const p = piCurve.points[8]
      testPointT(piCurve, p, 8)
    })

    test("3", () => {
      const piCurve = PICurve.forParametricStartEnd(
        new ProjectedCurveSurface(
          new BezierCurve(
            V(0.01100000000000001, 0.28900000000000003, 1.2),
            V(-0.04100000000000001, 0.28900000000000003, 1.2),
            V(-0.1, 0.279, 1.2),
            V(-0.16499999999999998, 0.255, 1.2),
            0,
            1,
          ),
          V(0, 0, -1),
          0,
          1,
          0,
          1,
        ),
        new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z),
        V(0, 0.2427341017251267, 0),
        V(1, 0.24724084890251513, 0),
        0.05,
        V(0.04999991966434234, -0.00008963012502444188, 0),
        0,
        20,
      )
      const p = V(-0.16499999999999998, 0.255, 0.9527591510974849)
      testPointT(piCurve, p, 20)
    })

    test("4", () => {
      const piCurve = PICurve.forParametricStartEnd(
        new ProjectedCurveSurface(
          new BezierCurve(
            V(4, 13.8, 15.7),
            V(3, 13.5, 13),
            V(1, 12.8, 6.9),
            V(0, 12.5, 4.200000000000001),
            -0.1,
            1.1,
          ),
          V(-1, -0.8, -5.5),
          0,
          1,
          -1,
          0,
        ),
        new CylinderSurface(
          new EllipseCurve(
            V(2.5, 13.4, 11.65),
            V(-0.2, -0.16000000000000003, -1.1),
            V(0.2, 0.08000000000000002, 0.68),
            0,
            3.141592653589793,
          ),
          V(6.123233995736766e-17, 0.30000000000000004, 2.1000000000000005),
          0,
          3.141592653589793,
          0,
          2,
        ),
        V(0.5002392685936982, -0.30000289827476356, 0),
        V(0.501457710232049, -0.6998923972267981, 0),
        0.05,
        V(-0.04998533362291978, 0.0012109593739698047, 0),
        0.018324617645703256,
        9.887672300683334,
      )
      const p = V(2.700000000000252, 13.710000000000084, 13.80000000000074)
      testPointT(piCurve, p)
    })

    test("5", () => {
      const piCurve = PICurve.forParametricStartEnd(
        new ProjectedCurveSurface(
          new BezierCurve(
            V(0.3333333333333333, -0.1111111111111111, -1.333333333333333),
            V(0.3333333333333333, 0.011618833295731858, -1.333333333333333),
            V(0.2338410555179541, 0.1111111111111111, -1.333333333333333),
            V(0.11111111111111112, 0.1111111111111111, -1.333333333333333),
            0,
            1,
          ),
          V(0, 0, -2.222222222222222),
          0,
          1,
          -1,
          0,
        ),
        new EllipsoidSurface(
          V3.O,
          V(0.9999999999999999, 0, 0),
          V(0, 0.9999999999999999, 0),
          V(0, 0, -0.9999999999999999),
        ),
        V(1.002075784554904, -0.15555163884043577, 0),
        V(-0.047646525196279994, -0.17960385445412552, 2.710505431213761e-20),
        -0.05,
        V(-0.04999136396194676, -0.000929262731508745, 2.710505431213761e-20),
        0,
        21,
      )
      const p = V(
        0.30360395660901607,
        -6.873697213268433e-13,
        -0.9527983194419218,
      )
      testPointT(piCurve, p, 13.449200584902428)
    })

    test("6", () => {
      const piCurve = PICurve.forParametricStartEnd(
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
        V(0.5103488594210569, -0.3054969051428064, 0),
        V(0.49441919231892134, -0.6984170426639139, 0),
        0.05,
        V(-0.03401735422540411, 0.03664450315536262, 0),
        0,
        11,
      )
      const p = V(0.49999999999999967, 0, 0.30000000000000004)
      testPointT(piCurve, p)
    })
  })

  test("is curves", () => {
    const pcs = new ProjectedCurveSurface(
      new BezierCurve(
        V(4, 13.8, 15.7),
        V(3, 13.5, 13),
        V(1, 12.8, 6.9),
        V(0, 12.5, 4.200000000000001),
        -0.1,
        1.1,
      ),
      V(-1, -0.8, -5.5),
      0,
      1,
      -1,
      0,
    )
    const ses = new CylinderSurface(
      new EllipseCurve(
        V(2.5, 13.4, 11.65),
        V(-0.2, -0.16000000000000003, -1.1),
        V(0.2, 0.08000000000000002, 0.68),
        0,
        3.141592653589793,
      ),
      V(6.123233995736766e-17, 0.30000000000000004, 2.1000000000000005),
      0,
      3.141592653589793,
      0,
      2,
    )
    console.log(pcs.isCurvesWithSurface(ses).sce)
    testISCurves(pcs, ses, 1)
    const edges = [
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
        V(0.49999999999915834, 3.7878311598404935e-12, 0.3),
        V(0.5000000000011606, -5.222877685895355e-12, 0.7000000000000001),
        0.21324043802451342,
        10.421498564770445,
        undefined,
        V(-0.05, 0.22500000000000006, 4.261306993203557e-12),
        V(0.05, -0.225, -5.8757343350152756e-12),
        undefined,
      ),
    ]
    const piCurve = PICurve.forParametricStartEnd(
      new ProjectedCurveSurface(
        new BezierCurve(
          V(0.01100000000000001, 0.28900000000000003, 1.2),
          V(-0.04100000000000001, 0.28900000000000003, 1.2),
          V(-0.1, 0.279, 1.2),
          V(-0.16499999999999998, 0.255, 1.2),
          0,
          1,
        ),
        V(0, 0, -1),
        0,
        1,
        0,
        1,
      ),
      new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z),
      V(0, 0.2427341017251267, 0),
      V(1, 0.24724084890251513, 0),
      0.05,
      V(0.04999991966434234, -0.00008963012502444188, 0),
      0,
      20,
    )
    const pic2 = new PCurveEdge(
      PICurve.forParametricStartEnd(
        new ProjectedCurveSurface(
          new BezierCurve(
            V(0.01100000000000001, 0.28900000000000003, 1.2),
            V(-0.04100000000000001, 0.28900000000000003, 1.2),
            V(-0.1, 0.279, 1.2),
            V(-0.16499999999999998, 0.255, 1.2),
            0,
            1,
          ),
          V(0, 0, -1),
          0,
          1,
          0,
          1,
        ),
        new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z),
        V(0, 0.2427341017251267, 0),
        V(1, 0.24724084890251513, 0),
        0.05,
        V(0.04999991966434234, -0.00008963012502444188, 0),
        0,
        20,
      ),
      V(0.01100000000000001, 0.28900000000000003, 0.9572658982748733),
      V(-0.16499999999999998, 0.255, 0.9527591510974849),
      0,
      20,
      undefined,
      V(-0.007799987467637408, 0, 0.00008963012502444188),
      V(-0.009748975193987298, -0.0035996216100876227, -0.0007249233929057179),
      undefined,
    )
    ;[
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
        V(0.49441919231892134, -0.6984170426639139, 0),
        V(0.5293754631948687, -0.34986900885706007, 0),
        0.05,
        V(0.04344758580313105, -0.02474484366245956, 0),
        0,
        0.9878894241992384,
      ),
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
        V(0.5293754631948687, -0.34986900885706007, 0),
        V(0.49441919231892134, -0.6984170426639139, 0),
        0.05,
        V(-0.0122531615274917, 0.048475354898971056, 0),
        2.011113594053313,
        10.99999999976717,
      ),
    ]
  })

  test("isTsWithSurface", () => {
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
    const piCurve = ses.isCurvesWithSurface(pcs)[0]
    console.log(piCurve)
    testISTs(piCurve, new PlaneSurface(P3.XY), 1)
  })

  test("isTsWithPlane", () => {
    const piCurve = PICurve.forParametricStartEnd(
      new RotatedCurveSurface(
        new EllipseCurve(
          V(3, 0, 0),
          V3.X,
          V(0, 6.123233995736766e-17, 1),
          0,
          3.141592653589793,
        ),
        M4.forSys(
          V(1, -1.2060416625018976e-16, -2.1265768495757713e-17),
          V(-1.2246467991473532e-16, -0.984807753012208, -0.17364817766693033),
          V(0, -0.17364817766693033, 0.984807753012208),
        ),
        0,
        3.141592653589793,
      ),
      new PlaneSurface(new P3(V(-1, 0, 0), 2.5), V(0, -1, 0), V3.Z),
      V(2.2459278597319283, 0, 0),
      V(3.141592653589793, 2.0943951023931957, 0),
      0.02,
      V(0, 0.02, 0),
      0,
      122,
    )
    const plane = new P3(
      V(0, -0.1220799925322713, -0.9925202644900105),
      0.8281304558850031,
    )
    testISTs(piCurve, plane, 1)
  })

  test("testCurve", () => {
    const curve = PICurve.forParametricStartEnd(
      new ProjectedCurveSurface(
        new BezierCurve(
          V(0.30000000000000004, -0.1, -1.2),
          V(0.30000000000000004, 0.010456949966158674, -1.2),
          V(0.2104569499661587, 0.1, -1.2),
          V(0.10000000000000002, 0.1, -1.2),
          0,
          1,
        ),
        V(0, 0, 2),
        0,
        1,
        0,
        1,
      ),
      new EllipsoidSurface(V3.O, V3.X, V3.Y, V(0, 0, -1)),
      V(1.0013147208257906, 0.10500326852305793, 0),
      V(-0.04846385366656537, 0.12648109268323876, 0),
      0.05,
      V(-0.049993023321055666, 0.0008352360267516577, 0),
      0,
      13.433240755090269,
    )
    testCurve(curve, true, "curve")
    const curveFoo = curve.transform(M4.FOO)
    testCurve(curveFoo, true, "curve.foo()", {
      drPs: [V(1.111182481682829, 12.285190464999221, 4.054554240538709)],
    })
    const roots = curveFoo.roots()
    expect(roots).toEqual([
      [1.4160324062686414],
      [11.253090623533353],
      [9.631006746785715],
    ])

    for (const dim of [0, 1, 2]) {
      for (const root of roots[dim]) {
        const tangent = curveFoo.tangentAt(root)
        expect(
          tangent.e(dim),
          dedent`
            dim:${dim}
            root:${root}
            p:${curveFoo.at(root)}
            tangent:${tangent}`,
        ).toFuzzyEqual(0)
      }
    }
  })

  test("testCurve 2", () => {
    const curve = PICurve.forParametricStartEnd(
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
      V(0.5103488594210569, -0.3054969051428064, 0),
      V(0.49441919231892134, -0.6984170426639139, 0),
      0.02,
      V(-0.03401735422540411, 0.03664450315536262, 0),
      2,
      10,
    )
    testCurve(curve, false, "curve")
    testCurve(curve.transform(M4.FOO), false, "curve.foo()")
  })
})
