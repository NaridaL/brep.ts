import {
  assert,
  assertVectors,
  callSource,
  getLast,
  M4,
  newtonIterate,
  Tuple3,
  V3,
} from "ts3dutils"

import {
  Curve,
  curvePointPP,
  followAlgorithmPP,
  ImplicitCurve,
  ImplicitSurface,
  ISInfo,
  P3,
  ParametricSurface,
  PlaneSurface,
  Surface,
  surfaceIsICurveIsInfosWithLine,
} from "../index"
import { abs, ceil, floor } from "../math"

export class PPCurve extends ImplicitCurve {
  constructor(
    points: ReadonlyArray<V3>,
    tangents: ReadonlyArray<V3>,
    readonly parametricSurface1: ParametricSurface,
    readonly parametricSurface2: ParametricSurface,
    readonly st1s: ReadonlyArray<V3>,
    readonly pmTangents: ReadonlyArray<V3> | undefined,
    readonly stepSize: number,
    dir: number = 1,
    generator?: string,
    tMin?: number,
    tMax?: number,
  ) {
    super(points, tangents, dir, generator, tMin, tMax)
    assert(ParametricSurface.is(parametricSurface1))
    assert(ParametricSurface.is(parametricSurface2))
    assert(Array.isArray(st1s))
    assert(dir == 1)
    assert(stepSize <= 1)
  }

  at(t: number) {
    assert(!isNaN(t))
    if (0 === t % 1) return this.points[t]
    const startPoint = V3.lerp(
      this.points[floor(t)],
      this.points[ceil(t)],
      t % 1,
    )
    return curvePointPP(
      this.parametricSurface1,
      this.parametricSurface2,
      startPoint,
    )!.p
  }

  isColinearTo(curve: Curve) {
    if (curve instanceof PPCurve) {
      if (this.equals(curve)) {
        return true
      }
      if (
        this.parametricSurface1.isCoplanarTo(curve.parametricSurface1) &&
        this.parametricSurface1.isCoplanarTo(curve.parametricSurface2)
      ) {
        // TODO
      }
      return false
    } else {
      return false
    }
  }

  containsPoint(p: V3) {
    assertVectors(p)
    // TODO: wrong, as there could be another curve
    return (
      this.parametricSurface1.containsPoint(p) &&
      this.parametricSurface2.containsPoint(p) &&
      !isNaN(this.pointT(p))
    )
  }

  rootPoints() {
    const pF1 = this.parametricSurface1.pUVFunc()
    const pF2 = this.parametricSurface2.pUVFunc()
    const pN1 = this.parametricSurface1.normalUVFunc()
    const pN2 = this.parametricSurface2.normalUVFunc()

    const rootsApprox = this.rootsApprox()
    const results: Tuple3<V3[]> = [[], [], []]
    for (let dim = 0; dim < 3; dim++) {
      for (let i = 0; i < rootsApprox[dim].length; i++) {
        const lambda = rootsApprox[dim][i]
        const p = this.at(lambda)
        assert(this.parametricSurface1.containsPoint(p))
        const pp1 = this.parametricSurface1.uvP(p)
        const { x: u, y: v } = this.parametricSurface2.uvP(p)
        const startValues = [pp1.x, pp1.y, u, v]

        function f(vals: number[]) {
          const [u1, v1, u2, v2] = vals
          const diff = pF1(u1, v1).minus(pF2(u2, v2))
          const n1 = pN1(u1, v1)
          const n2 = pN2(u2, v2)
          const tangent = n1.cross(n2)
          return [diff.x, diff.y, diff.z, tangent.e(dim)]
        }

        const pps = newtonIterate(f, startValues, 8)
        // assert(pF1(pps[0], pps[1]).like(pF2(pps[2], pps[3])),
        // 	pF1(pps[0], pps[1]).sce + pF2(pps[2], pps[3]).sce)
        const result = pF1(pps[0], pps[1])
        results[dim].push(result)
      }
    }
    return results
  }

  roots(): [number[], number[], number[]] {
    return this.rootPoints().map((ps) => ps.map((p) => this.pointT(p)))
  }

  pointTangent(pWC: V3): V3 {
    assertVectors(pWC)
    assert(this.containsPoint(pWC), "this.containsPoint(pWC)")
    const n1 = this.parametricSurface1.normalP(pWC)
    const n2 = this.parametricSurface2.normalP(pWC)
    return n1.cross(n2)
  }

  transform(m4: M4) {
    return new PPCurve(
      m4.transformedPoints(this.points),
      m4.transformedVectors(this.tangents),
      this.parametricSurface1.transform(m4),
      this.parametricSurface2.transform(m4),
      this.st1s,
      undefined,
      this.stepSize,
      this.dir,
      undefined,
    ) as this
  }

  toSource(): string {
    return callSource(
      "PPCurve.forStartEnd",
      this.parametricSurface1,
      this.parametricSurface2,
      this.points[0],
      getLast(this.points),
      this.stepSize,
    )
  }

  static forStartEnd(
    ps1: ParametricSurface,
    ps2: ParametricSurface,
    startPoint: V3,
    end: V3,
    stepSize = 0.02,
  ) {
    const { points, tangents, st1s } = followAlgorithmPP(
      ps1,
      ps2,
      startPoint,
      stepSize,
    )
    return new PPCurve(points, tangents, ps1, ps2, st1s, undefined, stepSize, 1)
  }

  isInfosWithLine(
    anchorWC: V3,
    dirWC: V3,
    tMin?: number | undefined,
    tMax?: number | undefined,
    lineMin?: number | undefined,
    lineMax?: number | undefined,
  ): ISInfo[] {
    return surfaceIsICurveIsInfosWithLine.call(
      this,
      this.parametricSurface1,
      this.parametricSurface2,
      anchorWC,
      dirWC,
      tMin,
      tMax,
      lineMin,
      lineMax,
    )
  }

  isTsWithSurface(surface: Surface): number[] {
    if (ImplicitSurface.is(surface)) {
      const result: number[] = []
      const iF = surface.implicitFunction()
      const pUV1 = this.parametricSurface1.pUVFunc()
      const pUV2 = this.parametricSurface2.pUVFunc()
      let prevSignedDistance = iF(this.points[0])
      for (let i = 1; i < this.points.length; i++) {
        const point = this.points[i]
        const signedDistance = iF(point)
        if (prevSignedDistance * signedDistance <= 0) {
          const startIndex =
            abs(prevSignedDistance) < abs(signedDistance) ? i - 1 : i
          const startPoint = this.points[startIndex]
          const startUV1 = this.st1s[startIndex]
          const startUV2 = this.parametricSurface2.uvP(startPoint)
          const isSTUV = newtonIterate(
            ([u1, v1, u2, v2]: number[]) => {
              const ps1p = pUV1(u1, v1)
              const ps2p = pUV2(u2, v2)
              return [...ps1p.to(ps2p), iF(ps1p)]
            },
            [startUV1.x, startUV1.y, startUV2.x, startUV2.y],
          )
          result.push(
            this.pointT(this.parametricSurface1.pUV(isSTUV[0], isSTUV[1])),
          )
        }
        prevSignedDistance = signedDistance
      }
      return result
    }
    throw new Error("Method not implemented.")
  }

  isTsWithPlane(planeWC: P3): number[] {
    return this.isTsWithSurface(new PlaneSurface(planeWC))
  }
}
