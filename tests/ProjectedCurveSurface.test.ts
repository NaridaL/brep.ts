import {
    inDifferentSystems,
    suite,
    skip,
    test,
    testISTs,
    testParametricSurface,
    testZDirVolumeAndAreaAndArea,
    testCalculateArea,
    testZDirVolumeAndArea, surfaceVolumeAndArea, linkB3, testISCurves, testLoopCCW
} from './manager'

import {DEG, M4, V, V3} from 'ts3dutils'
import {
    BezierCurve,
    L3,
    P3,
    PCurveEdge,
    PlaneSurface,
    RotationFace,
    SemiEllipsoidSurface,
    StraightEdge,
    ProjectedCurveSurface,
    Edge,
    B2T,
    PICurve,
    SemiEllipseCurve,
    Face,
} from '..'
import {SemiCylinderSurface} from '../src/surface/SemiCylinderSurface'

suite('ProjectedCurveSurface', () => {
    suite('ProjectedCurveSurface', inDifferentSystems((assert, m4) => {
        const baseCurve = BezierCurve.graphXY(2, -3, -3, 2)
        const pcs = new ProjectedCurveSurface(baseCurve, V3.Z, undefined, undefined, -100, 100).transform(m4)
        testParametricSurface(assert, pcs)
    }))
    const curve = BezierCurve.graphXY(2, -3, -3, 2)
    const edge = PCurveEdge.forCurveAndTs(curve, 0, 1)
    const edges = [
        edge,
        StraightEdge.throughPoints(curve.at(1), curve.at(1).plus(V(0, 0, 10))),
        edge.flipped().transform(M4.translate(0, 0, 10)),
        StraightEdge.throughPoints(curve.at(0).plus(V(0, 0, 10)), curve.at(0))]
    const surface = new ProjectedCurveSurface(curve, V3.Z)
    const pcsFace = new RotationFace(surface, edges)

    suite('Face line intersection test', inDifferentSystems((assert, m4) => {
        const line = new L3(V3.Z, V3.X).transform(m4)
        const d = pcsFace.transform(m4).intersectsLine(line)
        assert.ok(d)
    }))
    const pcs = new ProjectedCurveSurface(BezierCurve.EX2D, V3.Z, undefined, undefined, -1, 4)
    test('testSurface', assert => {
        testParametricSurface(assert, pcs)
    })
    test('isTsWithSurface', assert => {
        const pcs = new ProjectedCurveSurface(BezierCurve.EX2D, V3.Z, undefined, undefined, -2, 2).scale(0.2, 0.2, 1).rotateX(-90 * DEG)
        const ses = SemiEllipsoidSurface.UNIT
        const pic = ses.isCurvesWithSurface(pcs)[0]
        testISTs(assert, pic, new PlaneSurface(P3.XY), 1)
    })
    suite('area and volume', surfaceVolumeAndArea(pcsFace))
    suite('area and volume 2', surfaceVolumeAndArea(pcsFace.transform(M4.FOO)))

    // create a pcs face which includes a PICurve
    const bezierEdge = Edge.forCurveAndTs(BezierCurve.EX2D, 0, 1)


    const a = B2T.extrudeEdges([bezierEdge, ...StraightEdge.chain([bezierEdge.b, V3.X.negated(), bezierEdge.a], false)],
        P3.XY.flipped())
    const b = B2T.cylinder(0.2, 2).rotateY(90 * DEG).translate(0, 0, 0.5)
    const loop = [
        new PCurveEdge(PICurve.forParametricStartEnd(new ProjectedCurveSurface(new BezierCurve(V(0, 2, 0),V(0.3333333333333333, 1, 0),V(0.6666666666666666, -1, 0),V(1, -2, 0),-0.1,1.1),V(0, 0, -1),0,1,-1,0),new SemiCylinderSurface(new SemiEllipseCurve(V(0, 0, 0.5),V(1.2246467991473533e-17, 0, -0.2),V(0, 0.2, 0),0,3.141592653589793),V(1, 0, 6.123233995736766e-17),0,3.141592653589793,0,2),V(0.543312790740714, -0.5455452074045042, 0),V(0.45794084603336055, -0.4349260970573627, 0),0.05,V(0.0026015534861132676, 0.04993227332556462, 0),6.22402514974994,10),V(0.49999999999999994, 0, 0.30000000000000004),V(0.4579408460333606, 0.1891173898820845, 0.4349260970573627),6.224025149750126,10,undefined,V(-0.04683824539832595, 0.21076813481893145, -0.0014356938902703453),V(-0.003821083373315509, 0.017154318833024252, 0.04985377941394385),'genseg2'),
            new PCurveEdge(PICurve.forParametricStartEnd(new ProjectedCurveSurface(new BezierCurve(V(0, 2, 0),V(0.3333333333333333, 1, 0),V(0.6666666666666666, -1, 0),V(1, -2, 0),-0.1,1.1),V(0, 0, -1),0,1,-1,0),new SemiCylinderSurface(new SemiEllipseCurve(V(0, 0, 0.5),V(1.2246467991473533e-17, 0, -0.2),V(0, 0.2, 0),0,3.141592653589793),V(1, 0, 6.123233995736766e-17),0,3.141592653589793,0,2),V(0.45794084603336055, -0.4349260970573627, 0),V(0.543312790740714, -0.5455452074045042, 0),0.05,V(-0.003821083373315509, -0.04985377941394385, 0),0,6.3461388249955215),V(0.4579408460333606, 0.1891173898820845, 0.4349260970573627),V(0.49999999999999994, 2.4492935982947065e-17, 0.7),0,6.3461388249955215,undefined,V(-0.003821083373315509, 0.017154318833024252, 0.04985377941394385),V(0.04658715677301388, -0.209638049867369, 0.0007477447876345841),'genseg3'),
            new PCurveEdge(PICurve.forParametricStartEnd(new ProjectedCurveSurface(new BezierCurve(V(0, 2, 0),V(0.3333333333333333, 1, 0),V(0.6666666666666666, -1, 0),V(1, -2, 0),-0.1,1.1),V(0, 0, -1),0,1,-1,0),new SemiCylinderSurface(new SemiEllipseCurve(V(0, 0, 0.5),V(-1.2246467991473533e-17, 2.4492935982947065e-17, 0.2),V(-1.4997597826618578e-33, -0.2, 2.4492935982947065e-17),0,3.141592653589793),V(1, 0, 6.123233995736766e-17),0,3.141592653589793,0,2),V(0.45794084603336055, -0.43492609705736274, 0),V(0.543312790740714, -0.5455452074045042, 0),0.05,V(-0.0038210833733155047, -0.04985377941394385, 0),6.3461388249955215,11),V(0.49999999999999994, 2.4492935982947065e-17, 0.7),V(0.543312790740714, -0.19474504892931377, 0.5455452074045042),6.346138825001836,11,undefined,V(0.04658715677301388, -0.209638049867369, 0.0007477447876345841),V(0.002601553486113264, -0.011677707635158832, -0.04993227332556462),'genseg5'),
            new PCurveEdge(PICurve.forParametricStartEnd(new ProjectedCurveSurface(new BezierCurve(V(0, 2, 0),V(0.3333333333333333, 1, 0),V(0.6666666666666666, -1, 0),V(1, -2, 0),-0.1,1.1),V(0, 0, -1),0,1,-1,0),new SemiCylinderSurface(new SemiEllipseCurve(V(0, 0, 0.5),V(-1.2246467991473533e-17, 2.4492935982947065e-17, 0.2),V(-1.4997597826618578e-33, -0.2, 2.4492935982947065e-17),0,3.141592653589793),V(1, 0, 6.123233995736766e-17),0,3.141592653589793,0,2),V(0.543312790740714, -0.5455452074045042, 0),V(0.45794084603336055, -0.43492609705736274, 0),0.05,V(0.002601553486113264, 0.04993227332556462, 0),0,6.224025149750129),V(0.543312790740714, -0.19474504892931377, 0.5455452074045042),V(0.49999999999999994, 0, 0.30000000000000004),0,6.224025149750129,undefined,V(0.002601553486113264, -0.011677707635158832, -0.04993227332556462),V(-0.04683824539832592, 0.2107681348189313, -0.0014356938902701979),'genseg4')
        ]
    const piCurveFace = Face.create(a.faces.find(f => f.surface instanceof ProjectedCurveSurface).surface, loop)
    test('ccw for round face', assert => {
        testLoopCCW(assert, a.faces.find(f => f.surface instanceof ProjectedCurveSurface).surface, loop)
    })
    test('isCurvesWithCylinderSurface', assert => {
        testISCurves(assert,
            a.faces.find(f => f.surface instanceof ProjectedCurveSurface).surface,
            b.faces.find(f => f.surface instanceof SemiCylinderSurface).surface,
            2)
        const c = a.and(b)
        linkB3(assert, {a, b, c: c.translate(10)})
        console.log(c.toSource())
    })

    suite('area and volume with PICurves', surfaceVolumeAndArea(piCurveFace))
    suite('area and volume with PICurves 2', surfaceVolumeAndArea(piCurveFace.transform(M4.FOO)))

    test('Face containsPoint', assert => {
        const face = new RotationFace(
            new ProjectedCurveSurface(
                new BezierCurve(
                    V(142.87578921496748, -191.46078243076332, 0),
                    V(161.78547089700214, -252.13248349581008, 0),
                    V(284.63214994898954, -163.59789158697575, 0),
                    V(372.40411211189405, -210.3992206435476, 0), -3, 4), V(0, 0, 1), -3, 4, -100, 100),
            [
                PCurveEdge.forCurveAndTs(
                    new BezierCurve(
                        V(142.87578921496748, -191.46078243076332, 0),
                        V(161.78547089700214, -252.13248349581008, 0),
                        V(284.63214994898954, -163.59789158697575, 0),
                        V(372.40411211189405, -210.3992206435476, 0)), 1, 0),
                StraightEdge.throughPoints(
                    V(142.87578921496748, -191.46078243076332, 0),
                    V(142.87578921496748, -191.46078243076332, -100)),
                PCurveEdge.forCurveAndTs(
                    new BezierCurve(
                        V(142.87578921496748, -191.46078243076332, -100),
                        V(161.78547089700214, -252.13248349581008, -100),
                        V(284.63214994898954, -163.59789158697575, -100),
                        V(372.40411211189405, -210.3992206435476, -100)), 0, 1),
                StraightEdge.throughPoints(
                    V(372.40411211189405, -210.3992206435476, -100),
                    V(372.40411211189405, -210.3992206435476, 0))], [])
        const line = new L3(V(1241.5987, -1214.1894, 38.9886), V(-0.6705, 0.7386, -0.0696).unit())
        testISTs(assert, line, face.surface, 3)
    })
})