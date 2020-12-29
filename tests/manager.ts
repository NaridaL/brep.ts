import * as path from "path"
import prettier from "prettier"
import diff from "jest-diff"

try {
  ;(global as any).WebGLRenderingContext = {}
} catch (e) {}

import {
  arrayFromFunction,
  arraySamples,
  DEG,
  eq0,
  glqInSteps,
  int,
  lerp,
  M4,
  raddd,
  toSource,
  V,
  V3,
} from "ts3dutils"

export * from "ts3dutils/lib/manager"
import * as fs from "fs"
import slug from "slug"

import { RenderObjects } from "../src/viewer"

import {
  BRep,
  Curve,
  CustomPlane,
  Edge,
  edgeForCurveAndTs,
  EllipseCurve,
  Face,
  ImplicitSurface,
  L3,
  P3,
  ParametricSurface,
  PlaneSurface,
  PointVsFace,
  rotateCurve,
  StraightEdge,
  Surface,
} from "../src"
import * as brepts from "../src"
import * as ts3dutils from "ts3dutils"

declare global {
  namespace jest {
    // noinspection JSUnusedGlobalSymbols
    interface Matchers<R> {
      toMatchBRepSnapshot(additionalStuff?): R

      toEqualBRep(expected: BRep)
    }
  }
}
const brepToString = (x: BRep) =>
  prettier.format(x.toSource(false), {
    parser: "typescript",
    semi: false,
    trailingComma: "all",
  })

function toEqualBRep(actual: BRep, expected: BRep, additionalStuff: {}) {
  const x = expected.getAABB().max.x - actual.getAABB().min.x
  let actualTranslated = actual.translate(x === -Infinity ? 0 : x + 1)
  outputLink(
    Object.assign(
      {
        a: expected,
        b: actualTranslated,
      },
      additionalStuff,
    ),
  )
  let message = ""
  if (actual.faces.length !== expected.faces.length) {
    message = `actual.faces.length !== expected.faces.length ${actual.faces.length} !== ${expected.faces.length}`
  }
  actual.faces.forEach((face) => {
    if (!expected.faces.some((expectedFace) => expectedFace.likeFace(face))) {
      message += "Unexpected face in result:" + face.toSource() + "\n"
    }
  })
  return {
    pass: !message,
    message: () =>
      message + "\n" + diff(brepToString(expected), brepToString(actual)),
  }
}

// noinspection JSUnusedGlobalSymbols
expect.extend({
  toMatchBRepSnapshot(actual: BRep, additionalStuff = {}) {
    const fileName = path.join(
      __dirname,
      "snapshots",
      sanitizeFilename(expect.getState().currentTestName) + ".ts",
    )
    console.log("snapshot", fileName)
    if (!fs.existsSync(fileName)) {
      // file does not exist
      fs.writeFileSync(fileName, brepToString(actual), "utf8")
      return { pass: true, message: () => "" }
    }
    const destructure = (x) => "{" + Object.keys(x).join(",") + "}"
    const expected: BRep = new Function(
      destructure(ts3dutils),
      destructure(brepts),
      "return " + fs.readFileSync(fileName, "utf8"),
    )(ts3dutils, brepts)

    return toEqualBRep(actual, expected, additionalStuff)
  },

  toEqualBRep(actual: BRep, expected: BRep, additionalStuff = {}) {
    return toEqualBRep(actual, expected, additionalStuff)
  },
})

const sanitizeFilenameOptions = {
  multicharmap: Object.defineProperties(Object.assign({}, slug.multicharmap), {
    str: {
      value: "str",
    },
    sce: {
      value: "sce",
    },
  }),
  charmap: Object.assign({}, slug.charmap, {
    "-": "minus",
    "+": "plus",
  }),
  replacement: "_",
}

function sanitizeFilename(s: string) {
  return slug(s, sanitizeFilenameOptions)
}

export function testCurveCentralProjection(curve: Curve) {
  const m4 = M4.projectPlanePoint(V3.O, new P3(V3.Z, 1))
  const pls = arraySamples(curve.tMin, curve.tMax, 16).flatMap((t) => {
    const p = curve.at(t)
    return [p, m4.transformPoint(p)]
  })
  let curveTransformed
  try {
    curveTransformed = (curve as any).transform4(m4) as Curve
  } catch (e) {}
  outputLink({
    edges: [
      edgeForCurveAndTs(curve),
      ...(curveTransformed ? [edgeForCurveAndTs(curveTransformed)] : []),
    ],
    drPs: pls,
    drLines: pls,
    planes: [
      new CustomPlane(
        V3.Z,
        V3.X,
        V3.Y,
        "Z=1",
        [0, 0, 0, 1],
        -100,
        100,
        -100,
        100,
      ),
    ],
  })
  if (curveTransformed) {
    arraySamples(curve.tMin, curve.tMax, 4).forEach((t) => {
      expect(
        curveTransformed.containsPoint(m4.transformPoint(curve.at(t))),
      ).toBeTruthy()
    })
  }
}

export function testBRepAnd(
  a: BRep,
  b: BRep,
  expected?: BRep,
  message?: string,
) {
  return testBRepOp(a, b, () => a.and(b), expected, message)
}

export function testBRepOp(
  a: BRep,
  b: BRep,
  calculateActual: () => BRep,
  expected?: BRep,
  message?: string,
) {
  let actual
  try {
    actual = calculateActual()
  } finally {
    const abWidth = a.getAABB().addAABB(b.getAABB()).size().x
    if (actual) {
      expect(actual).toMatchBRepSnapshot({
        c: a.translate(-abWidth - 1),
        d: b.translate(-abWidth - 1),
      })
    } else {
      outputLink({
        a,
        b,
      })
    }
  }
}

export function makeLink(values: any) {
  return Object.getOwnPropertyNames(values)
    .map((name) => {
      const val = values[name]
      return name + "=" + (typeof val == "string" ? val : val.toSource())
    })
    .join(";")
}

export function outputLink(
  values: { [K in keyof RenderObjects]?: RenderObjects[K] | string },
  msg = "view",
) {
  const script =
    "TEST_NAME = " +
    expect.getState().currentTestName.toSource() +
    "\n" +
    Object.getOwnPropertyNames(values)
      .map((name) => {
        const val = values[name]
        return (
          "const " +
          name +
          " = " +
          (typeof val == "string" ? val : toSource(val))
        )
      })
      .join("\n") +
    "\n" +
    "return {" +
    Object.keys(values).join(",") +
    "}"
  const o =
    sanitizeFilename(expect.getState().currentTestName + "_" + msg) + ".html"
  fs.mkdirSync(__dirname + "/results", { recursive: true })
  fs.writeFileSync(
    __dirname + "/results/" + o,
    demoFile.replace("/*INSERT*/", script),
    "utf8",
  )
  link("http://localhost:10001/tests/results/" + o, msg)
}

function link(url: string, msg?: string) {
  console.log(url, msg)
}

const demoFile = fs.readFileSync(__dirname + "/../viewer.html", "utf8")

export function testISCurves(
  surface1: Surface | P3,
  surface2: Surface | P3,
  expectedCurveCount: int,
): Curve[] | undefined {
  surface1 instanceof P3 && (surface1 = new PlaneSurface(surface1))
  surface2 instanceof P3 && (surface2 = new PlaneSurface(surface2))
  let isCurves: Curve[] | undefined
  try {
    isCurves = surface1.isCurvesWithSurface(surface2)
  } finally {
    if (isCurves) {
      outputLink({
        mesh: `[${surface1}.toMesh(), ${surface2}.toMesh()]`,
        edges: isCurves.map((c) => edgeForCurveAndTs(c)),
      })

      expect(isCurves.length).toBe(expectedCurveCount)
      for (const curve of isCurves) {
        expect(
          surface1.containsCurve(curve),
          "surface1.containsCurve(curve) " +
            surface1.toString() +
            " " +
            curve.toString(),
        ).toBeTruthy()
        expect(
          surface2.containsCurve(curve),
          "surface2.containsCurve(curve) " +
            surface2.toString() +
            " " +
            curve.toString(),
        ).toBeTruthy()
        const t = lerp(curve.tMin, curve.tMax, 0.5),
          p = curve.at(t),
          dp = curve.tangentAt(t)
        expect(
          surface1.containsPoint(p),
          "surface1.containsPoint(curve.at(curve.tMin))",
        ).toBeTruthy()
        expect(
          surface2.containsPoint(p),
          "surface2.containsPoint(curve.at(curve.tMin))",
        ).toBeTruthy()

        const pN1 = surface1.normalP(p)
        const pN2 = surface2.normalP(p)
        const expectedTangent = pN1.cross(pN2)
        // expectedTangent can be zero if the surfaces just touch and dont
        // cross each other !expectedTangent.likeZero() &&
        // assert.ok(expectedTangent.isParallelTo(dp),
        // 'pN1.cross(pN2).isParallelTo(dp)') !expectedTangent.likeZero() &&
        // assert.ok(expectedTangent.dot(dp) > 0, 'pN1.cross(pN2).dot(dp) > 0')
      }
    } else {
      outputLink({
        mesh: `[${surface1}.toMesh(), ${surface2}.toMesh()]`,
      })
      throw new Error("no isCurves returned: " + isCurves)
    }
  }
  return isCurves
}

export function testPointT(
  curve: Curve,
  p: V3,
  expectedT?: number,
  precision?: number,
) {
  outputLink({
    edges: [edgeForCurveAndTs(curve)],
    drPs: [p],
  })
  const actualT = curve.pointT(p)
  if (undefined !== expectedT) {
    if (isNaN(expectedT)) {
      expect(actualT).toBeNaN()
    } else {
      expect(actualT).toFuzzyEqual(expectedT, precision)
    }
  }
  if (!isNaN(actualT)) {
    expect(curve.at(actualT)).toBeLike(p)
  }
}

/**
 * Tests that the passed loop is CCW on the passed surface, and that the
 * reversed loop is CW.
 * @param assert
 * @param surface
 * @param loop
 */
export function testLoopCCW(surface: Surface, loop: Edge[]) {
  const points = [loop[0].a, loop[0].atAvgT()]
  outputLink(
    {
      mesh: `[${surface.toSource()}.toMesh()]`,
      edges: loop,
      drPs: points,
    },
    "testLoopCCW",
  )
  expect(surface.edgeLoopCCW(loop)).toBeTruthy()
  expect(!surface.edgeLoopCCW(Edge.reversePath(loop))).toBeTruthy()
}

export function surfaceVolumeAndAreaTests(
  face: Face,
  msg = "face",
  expectedVolume?: number,
) {
  const flippedFace = face.flipped()
  const faceMesh = face.toMesh()
  const faceMeshVol = faceMesh.calcVolume()

  test(msg + " area", () => {
    outputLink({ mesh: face.toSource() + ".toMesh()", face: [face] })
    const actualArea = face.calcArea()
    const expectedArea = faceMeshVol.area
    expect(actualArea).toFuzzyEqual(expectedArea, 0.05)
    console.log("OK! actual = " + actualArea + ", expected = " + expectedArea)
  })
  test(msg + " flipped() area", () => {
    outputLink({
      mesh: `[${face.flipped().toSource()}.toMesh()]`,
      face: [face.flipped()],
    })
    const actualArea = flippedFace.calcArea()
    const expectedArea = faceMeshVol.area
    expect(actualArea).toFuzzyEqual(expectedArea, 0.05)
    console.log("OK! actual = " + actualArea + ", expected = " + expectedArea)
  })
  test(msg + " volume", () => {
    outputLink({
      mesh: face.toSource() + ".toMesh()",
      edges: face.allEdges,
    })
    faceMesh.calcVolume()
    const actual = face.zDirVolume().volume,
      expected =
        undefined === expectedVolume ? faceMeshVol.volume : expectedVolume
    expect(actual).toFuzzyEqual(expected, 0.05)
    console.log(
      "OK! actual = " +
        actual +
        ", expected = " +
        expected +
        ", |dv| = " +
        (actual - expected),
    )
  })
  test(msg + " flipped() volume", () => {
    outputLink({
      mesh: flippedFace.toSource() + ".toMesh()",
      edges: flippedFace.allEdges,
    })
    const actual = flippedFace.zDirVolume().volume
    const expected =
      undefined === expectedVolume ? -faceMeshVol.volume : -expectedVolume
    expect(actual).toFuzzyEqual(expected, 0.05)
    console.log(
      "OK! actual = " +
        actual +
        ", expected = " +
        expected +
        ", |dv| = " +
        (actual - expected),
    )
  })
  test(msg + " centroid", () => {
    const actual = flippedFace.zDirVolume()
    const expected = faceMeshVol.centroid
    outputLink({
      mesh: face.toSource() + ".toMesh()",
      drPs: [expected, actual.centroid],
      edges: face.getAllEdges(),
    })
    if (!eq0(actual.volume)) {
      // centroid doesn't make sense when volume is 0
      expect(actual.centroid).toBeLike(expected, expected.length() / 100)
      console.log(
        "OK! actual = " +
          actual.centroid +
          ", expected = " +
          expected +
          ", |dv| = " +
          actual.centroid.distanceTo(expected),
      )
    }
  })
  test(msg + " flipped() centroid", () => {
    const actual = flippedFace.zDirVolume()
    const expected = faceMeshVol.centroid
    outputLink({
      mesh: face.toSource() + ".toMesh()",
      drPs: [expected, actual.centroid],
    })
    if (!eq0(actual.volume)) {
      // centroid doesn't make sense when volume is 0
      expect(actual.centroid).toBeLike(expected, expected.length() / 100)
    }
  })
  test(msg + " aabb", () => {
    const expected = faceMesh.getAABB()
    const actual = face.getAABB()
    expect(actual.min).toBeLike(expected.min, 0.1)
    expect(actual.max).toBeLike(expected.max, 0.1)
  })
}

export function testCurve(curve: Curve, checkTangents = true, msg?: string) {
  const edge = edgeForCurveAndTs(curve)
  const aabb = curve.getAABB()
  outputLink(
    {
      edges: [edge],
      drPs: [edge.a, edge.b],
      boxes: aabb ? [aabb.getM4()] : [],
    },
    msg,
  )
  const STEPS = 12
  for (let i = 0; i < STEPS; i++) {
    const t = lerp(curve.tMin, curve.tMax, i / (STEPS - 1))
    const p = curve.at(t)
    // check that pointT and containsPoint behave as expected
    expect(t).toFuzzyEqual(curve.pointT(p))
    expect(
      curve.containsPoint(p),
      `containsPoint(at(t == ${t}) == ${p})`,
    ).toBeTruthy()

    expect(aabb.containsPoint(p)).toBeTruthy()

    // check that tangentAt() behaves correctly
    if (checkTangents) {
      const eps = t != curve.tMax ? 1e-8 : -1e-8
      const expectedTangent = curve
        .at(t + eps)
        .minus(p)
        .div(eps)
      const actualTangent = curve.tangentAt(t)
      expect(actualTangent).toBeLike(expectedTangent, 1e-3)
    }
  }

  // test curve length
  if (curve.arcLength !== Curve.prototype.arcLength) {
    const expectedCurveLength = glqInSteps(
      (t) => curve.tangentAt(t).length(),
      curve.tMin,
      curve.tMax,
      4,
    )
    const actualCurveLength = curve.arcLength(curve.tMin, curve.tMax)
    expect(actualCurveLength).toFuzzyEqual(expectedCurveLength, 1e-6)
  }
}

export function suiteSurface(surface: Surface) {
  if (ParametricSurface.is(surface)) {
    test("ParametricSurface", () => testParametricSurface(surface))
  }
  if (ImplicitSurface.is(surface)) {
    test("ImplicitSurface", () => testImplicitSurface(surface))
  }
}

export function testParametricSurface(surf: ParametricSurface) {
  const debugInfo = surf.debugInfo && surf.debugInfo()
  outputLink(
    {
      mesh: `[${surf}.toMesh()]`,
      edges: [
        ...(debugInfo && debugInfo.lines
          ? arrayFromFunction(debugInfo.lines.length / 2, (i) => {
              const a = debugInfo.lines[i * 2],
                b = debugInfo.lines[i * 2 + 1]
              return !a.like(b) && StraightEdge.throughPoints(a, b)
            }).filter((x) => x)
          : []),
      ],
      drPs: [...((debugInfo && debugInfo.points) || [])],
    },
    "view",
  )

  expect(ParametricSurface.is(surf)).toBeTruthy()

  // test equals
  const clone = new surf.constructor(...surf.getConstructorParameters())
  expect(clone.equals(surf)).toBeTruthy()
  expect(clone.isCoplanarTo(surf)).toBeTruthy()
  expect(clone.like(surf)).toBeTruthy()

  const params = [
    V(0, 0),
    V(0, 1),
    V(1, 0),
    V(1, 1),
    V(0.25, 0.25),
    V(0.6, 0.25),
    V(0.25, 0.6),
    V(0.6, 0.7),
  ].map(
    (pm) =>
      new V3(
        lerp(surf.uMin, surf.uMax, pm.x),
        lerp(surf.vMin, surf.vMax, pm.y),
        0,
      ),
  )
  const points = params.map(({ x, y }) => surf.pUV(x, y))
  const psFlipped = surf.flipped()
  for (let i = 0; i < points.length; i++) {
    const p = points[i]
    const pm = params[i]
    expect(surf.containsPoint(p)).toBeTruthy()

    // test dpdu and dpdv
    const eps = 0.0001
    const dpdu = surf.dpdu()(pm.x, pm.y)
    const dpdv = surf.dpdv()(pm.x, pm.y)
    const dpduNumeric = p.to(surf.pUV(pm.x + eps, pm.y)).div(eps)
    const dpdvNumeric = p.to(surf.pUV(pm.x, pm.y + eps)).div(eps)
    expect(dpduNumeric).toBeLike(dpdu, 0.01)
    expect(dpdvNumeric).toBeLike(dpdv, 0.01)
    const pmNormal = surf.normalUV(pm.x, pm.y)
    expect(pmNormal.length()).toFuzzyEqual(1)
    const dpduXdpdv = dpdu.cross(dpdv)
    if (ImplicitSurface.is(surf)) {
      expect(eq0(surf.implicitFunction()(p))).toBeTruthy()
    }
    const pm2 = surf.uvP(p)
    const pNormal = surf.normalP(p)
    const psFlippedUV = psFlipped.uvP(p)
    const psFlippedNormal = psFlipped.normalP(p)
    if (!dpdu.likeO() && !dpdv.likeO()) {
      expect(pm2).toBeLike(pm) // pm == uvP(pUV(pm))
      expect(psFlipped.pUV(psFlippedUV.x, psFlippedUV.y)).toBeLike(p)

      expect(pmNormal).toBeLike(pNormal)

      const computedNormal = dpduXdpdv.unit()
      expect(computedNormal.angleTo(pNormal) < 5 * DEG).toBeTruthy()

      expect(psFlippedNormal).toBeLike(pNormal.negated())
    }

    // test pointFoot:
    const offsetPoint = p.plus(pNormal.toLength(0.5))
    const actualFoot = surf.pointFoot(offsetPoint)
    // the original params may not actually be the closest.
    //if (!actualFoot.like(pm) && surf.pUV(actualFoot.x,
    // actualFoot.y).distanceTo(p) > 0.5) { assert.v3like(actualFoot, pm,
    // 'distance to foot ' + surf.pUV(actualFoot.x,
    // actualFoot.y).distanceTo(p)) }
    const offsetPoint2 = p.plus(pNormal.toLength(-0.5))
    const actualFoot2 = surf.pointFoot(offsetPoint2)
    // the original params may not actually be the closest.
    //if (!actualFoot2.like(pm) && surf.pUV(actualFoot2.x,
    // actualFoot2.y).distanceTo(p) > 0.5) { assert.v3like(actualFoot2, pm,
    // 'distance to foot ' + surf.pUV(actualFoot2.x,
    // actualFoot2.y).distanceTo(p)) }
  }
  const matrices = [M4.mirror(P3.XY), M4.mirror(P3.YZ), M4.mirror(P3.ZX)]
  for (let mI = 0; mI < matrices.length; mI++) {
    const m = matrices[mI]
    for (let i = 0; i < points.length; i++) {
      const dpdu = surf.dpdu()(params[i].x, params[i].y)
      const dpdv = surf.dpdv()(params[i].x, params[i].y)
      const p = points[i],
        pNormal = surf.normalP(p)
      const normalMatrix = m.as3x3().inversed().transposed()
      const mNormal = normalMatrix.transformVector(pNormal)
      const mP = m.transformPoint(p)
      const mSurface = surf.transform(m)

      if (!dpdu.likeO() && !dpdv.likeO()) {
        expect(mSurface.normalP(mP)).toBeLike(mNormal)
      }

      expect(mSurface.containsPoint(mP)).toBeTruthy()

      //const mPSFlipped = mSurface.flipped()
      //assert.ok(mPSFlipped.normalP(mP).negated().like(mNormal))
      //assert(mPSFlipped.normalP(mP).negated().like(mNormal))
    }
  }
}

export function testContainsCurve(
  surface: Surface,
  curve: Curve,
  expected = true,
  msg?: string,
) {
  outputLink(
    {
      mesh: `[${surface.toSource()}.toMesh()]`,
      edges: [edgeForCurveAndTs(curve)],
    },
    msg,
  )
  expect(surface.containsCurve(curve)).toBe(expected)
}

export function rotateEdge(edge: Edge, angle: raddd) {
  const surface = rotateCurve(
    edge.curve,
    undefined,
    undefined,
    angle,
    edge.deltaT() > 0,
  )
  const edges = [
    edge,
    edgeForCurveAndTs(
      EllipseCurve.semicircle(edge.b.lengthXY(), V(0, 0, edge.b.z), 0, angle),
    ),
    edge.rotateZ(angle).flipped(),
    edgeForCurveAndTs(
      EllipseCurve.semicircle(edge.a.lengthXY(), V(0, 0, edge.a.z), 0, angle),
    ).flipped(),
  ]
  return Face.create(surface, edges)
}

function testImplicitSurface(surface: ImplicitSurface) {
  const EPS = 1e-8
  const testPoints = [
    V3.O.plus(V(0.2, 0, 0)), // V3.O fails on ellipsoidSurface
    V3.Y,
    V3.X,
    V3.Z.plus(V(0.2, 0, 0)),
    V3.XY,
    V3.XYZ,
    new V3(10, 10, 10),
    new V3(5, 6, 7),
  ]
  for (const testPoint of testPoints) {
    const i = surface.implicitFunction()(testPoint)
    const didpGuess = testPoint.map((el, dim) => {
      const i2 = surface.implicitFunction()(
        testPoint.plus(V3.O.withElement(dim, EPS)),
      )
      return (i2 - i) / EPS
    })
    const didp = surface.didp(testPoint)
    expect(didp).toBeLike(didpGuess, 1e-5)
  }
}

export function testCurveISInfos(
  c1: Curve,
  c2: Curve,
  expectedCount: int,
  msg: string = "view",
  f = "isInfosWithCurve",
) {
  let intersections
  try {
    intersections = c1[f](c2).map((info) => info.p)
    outputLink(
      {
        edges: [c1, c2].map((c) => edgeForCurveAndTs(c)),
        drPs: intersections,
      },
      msg,
    )
    expect(intersections.length).toBe(expectedCount)
    intersections.forEach((is, i) => {
      expect(
        intersections.every((is2, j) => j == i || !is.like(is2)),
        is.toSource() + " is not unique " + intersections,
      ).toBeTruthy()
      expect(
        c1.containsPoint(is),
        `e1.containsPoint(is): ${c1.toSource()}.containsPoint(${is.toSource()},`,
      ).toBeTruthy()
      expect(
        c2.containsPoint(is),
        `e2.containsPoint(is): ${c1.toSource()}.containsPoint(${is.toSource()},`,
      ).toBeTruthy()
    })
  } finally {
    !intersections &&
      outputLink({ edges: [c1, c2].map((c) => edgeForCurveAndTs(c)) }, msg)
  }
}

/**
 * Test intersections of a curve with a surface.
 * @param assert
 * @param curve
 * @param surface
 * @param expectedTCount
 */
export function testISTs(
  curve: Curve,
  surface: Surface | P3,
  expectedTCount: int,
) {
  surface instanceof P3 && (surface = new PlaneSurface(surface))
  let ists
  try {
    ists =
      curve instanceof L3
        ? surface.isTsForLine(curve)
        : curve.isTsWithSurface(surface)
    expect(ists.length).toBe(expectedTCount)
    for (const t of ists) {
      const p = curve.at(t)
      expect(surface.containsPoint(p)).toBeTruthy()
      // assert.ok(
      //   surface.containsPoint(p),
      //   "surface.containsPoint(p) " +
      //     surface.toString() +
      //     " " +
      //     p.toString() +
      //     " t: " +
      //     t +
      //     (ImplicitSurface.is(surface)
      //       ? " dist: " + surface.implicitFunction()(p)
      //       : ""),
      // )
    }
  } finally {
    outputLink(
      {
        mesh: `[${surface}.toMesh()]`,
        edges: [edgeForCurveAndTs(curve)],
        drPs: ists ? ists.map((t) => curve.at(t)) : [],
      },
      "view",
    )
  }
}

export function linkBRep(hash: string, message = "view") {
  const escapedHash = encodeURIComponent(
    hash.replace(/, /g, ",").replace(/([\n\t])+/g, ""),
  )
    .replace(/\(/g, "%28")
    .replace(/\)/g, "%29")
  link("http://localhost:10001/viewer.html#" + escapedHash, message)
}

export function testLoopContainsPoint(
  surface: Surface,
  loop: Edge[],
  p: V3,
  result: PointVsFace,
) {
  outputLink({
    mesh: `[${Face.create(surface, loop).toSource()}.toMesh()]`,
    drPs: [p],
  })
  expect(surface.loopContainsPoint(loop, p)).toBe(result)
}

export const FONT_PATH = __dirname + "/../fonts"

export function testCurvesColinear(curve1: Curve, curve2: Curve): void {
  expect(curve1.isColinearTo(curve2)).toBeTruthy()
  const t = (curve1.tMin + curve1.tMax) / 2
  expect(
    curve1
      .translate(curve1.tangentAt(t).getPerpendicular().unit())
      .isColinearTo(curve2),
  ).toBeFalsy()
  outputLink({
    edges: [curve1, curve2].map((c) => edgeForCurveAndTs(c)),
  })
  for (let i = 0; i < 10; i++) {
    const t = lerp(curve1.tMin, curve1.tMax, i / 9)
    expect(curve2.containsPoint(curve1.at(t))).toBeTruthy()
  }
}

export function testCurveTransform(
  curve: Curve,
  m4: M4,
  allowFlipped = false,
  msg?: string,
) {
  let curveT: Curve
  try {
    curveT = (curve.transform4 || curve.transform).call(curve, m4)
    console.log(curveT.sce)
  } finally {
    const ss = arraySamples(curve.tMin, curve.tMax, 16).flatMap((t) => [
      { p: curve.at(t), color: "red" },
      { p: m4.transformPoint(curve.at(t)), color: "green" },
    ])
    outputLink(
      {
        edges: [curve, curveT]
          .filter((x) => x)
          .map((c) => edgeForCurveAndTs(c)),
        drPs: ss,
        drLines: ss.map((x) => x.p),
        boxes: [curve.getAABB().getM4(), m4.times(curve.getAABB().getM4())],
        planes: ((x) => (x ? [CustomPlane.forPlane(x)] : []))(
          P3.vanishingPlane(m4),
        ),
      },
      msg,
    )
  }
  const cPTMin = m4.transformPoint(curve.at(curve.tMin))
  const cPTMax = m4.transformPoint(curve.at(curve.tMax))
  const cTPMin = curveT.at(curveT.tMin)
  const cTPMax = curveT.at(curveT.tMax)
  if (!allowFlipped || cPTMin.like(cTPMin)) {
    expect(cTPMin).toBeLike(cPTMin)
    expect(cTPMax).toBeLike(cPTMax)
  } else {
    expect(cTPMin).toBeLike(cPTMax)
    expect(cTPMax).toBeLike(cPTMin)
  }
  arraySamples(curve.tMin, curve.tMax, 8).forEach((t) => {
    const pT = m4.transformPoint(curve.at(t))
    expect(
      curveT.containsPoint(pT),
      pT.toSource() + curveT.distanceToPoint(pT),
    ).toBeTruthy()
  })
  return curveT
}

export function testSurfaceTransform(
  surface: ParametricSurface,
  transform: M4,
  allowFlipped = false,
  msg?: string,
) {
  const uvs = arraySamples(surface.uMin, surface.uMax, 4).flatMap((u) =>
    arraySamples(surface.vMin, surface.vMax, 4).map((v) => new V3(u, v, 0)),
  )
  let surfaceT
  try {
    surfaceT = (surface.transform4 || surface.transform).call(
      surface,
      transform,
    )
  } finally {
    const ss = uvs.flatMap((uv) => [
      { p: surface.pUV(uv.x, uv.y), color: "orange" },
      {
        p: transform.transformPoint(surface.pUV(uv.x, uv.y)),
        color: "lightgrey",
      },
    ])
    const aabbM4 = surface.getApproxAABB().getM4()
    outputLink(
      {
        mesh:
          "[" +
          [surface, surfaceT]
            .filter((x) => x)
            .map((s) => s.toSource() + ".toMesh()")
            .join(",") +
          "]",
        drPs: ss,
        drVs: [],
        drLines: ss.map((x) => x.p),
        boxes: [aabbM4, transform.times(aabbM4)],
        edges: [],
      },
      msg,
    )
  }
  uvs.forEach((uv) => {
    const pUV = surface.pUV(uv.x, uv.y)
    const pUVT = transform.transformPoint(pUV)
    const surfaceNT = P3.normalOnAnchor(
      surface.normalUV(uv.x, uv.y),
      pUV,
    ).transform(transform).normal1
    const surfaceTN = surfaceT.normalP(pUVT).unit()
    expect(
      surfaceT.containsPoint(transform.transformPoint(surface.pUV(uv.x, uv.y))),
    ).toBeTruthy()
    expect(
      surfaceNT.isParallelTo(surfaceTN),
      "uv " +
        uv.toSource() +
        " " +
        surfaceNT.toSource() +
        " " +
        surfaceTN.toSource(),
    ).toBeTruthy()
    expect(surfaceNT.dot(surfaceTN) > 0, "normal same dir").toBeTruthy()
  })
  return surfaceT
}
