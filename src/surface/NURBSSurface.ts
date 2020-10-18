import {
  arrayFromFunction,
  assert,
  clamp,
  firstUnsorted,
  indexWithMax,
  int,
  lerp,
  M4,
  MINUS,
  newtonIterate,
  sliceStep,
  toSource,
  V,
  V3,
  Vector,
} from "ts3dutils"
import { Curve, L3, NURBS, ParametricSurface, Surface } from ".."
import { Edge } from "../Edge"
import { P3 } from "../P3"
import { PointVsFace } from "./Surface"
import { ceil, floor } from "../math"

export class NURBSSurface extends ParametricSurface {
  constructor(
    /**
     * Control points in u-major order. I.e. the first pointCountU points are a NURBS.
     */
    readonly points: ReadonlyArray<Vector>,
    readonly knotsU: ReadonlyArray<number>,
    readonly knotsV: ReadonlyArray<number>,
    readonly degreeU: int,
    readonly degreeV: int,
    uMin: number = knotsU[degreeU],
    uMax: number = knotsU[knotsU.length - degreeU - 1],
    vMin: number = knotsV[degreeV],
    vMax: number = knotsV[knotsV.length - degreeV - 1],
  ) {
    super(uMin, uMax, vMin, vMax)
    const pointCountU = knotsU.length - 1 - degreeU
    const pointCountV = knotsV.length - 1 - degreeV
    assert(pointCountU * pointCountV == points.length)
    assert(degreeU <= degreeV, "degreeU <= degreeV")
    assert(
      -1 === firstUnsorted(knotsU, MINUS),
      "knot values must be in ascending order",
    )
    assert(
      -1 === firstUnsorted(knotsV, MINUS),
      "knot values must be in ascending order",
    )
  }

  getConstructorParameters() {
    return [
      this.points,
      this.knotsU,
      this.knotsV,
      this.degreeU,
      this.degreeV,
      this.uMin,
      this.uMax,
      this.vMin,
      this.vMax,
    ]
  }

  transform(m4: M4) {
    return this.transform4(m4) as this
  }

  transform4(m4: M4) {
    return new NURBSSurface(
      this.points.map((p) => m4.timesVector(p)),
      this.knotsU,
      this.knotsV,
      this.degreeU,
      this.degreeV,
      this.uMin,
      this.uMax,
      this.vMin,
      this.vMax,
    )
  }

  pUV(u: number, v: number) {
    return this.isoparametricU(u).at(v)
  }

  dpdu() {
    return (u: number, v: number) => this.isoparametricV(v).tangentAt(u)
  }

  dpdv() {
    return (u: number, v: number) => this.isoparametricU(u).tangentAt(v)
  }

  normalUV(u: number, v: number) {
    const normal = this.dpdu()(u, v).cross(this.dpdv()(u, v))
    return normal.likeO() ? V3.X : normal.unit()
  }

  isoparametricU(u: number) {
    const pointCountU = this.knotsU.length - 1 - this.degreeU
    const pointCountV = this.knotsV.length - 1 - this.degreeV
    return new NURBS(
      arrayFromFunction(pointCountV, (i) => {
        return deBoor(
          this.points.slice(i * pointCountU, (i + 1) * pointCountU),
          this.degreeU,
          this.knotsU,
          u,
        )
      }),
      this.degreeV,
      this.knotsV,
      this.vMin,
      this.vMax,
    )
  }

  isoparametricV(v: number): NURBS {
    const pointCountU = this.knotsU.length - 1 - this.degreeU
    return new NURBS(
      arrayFromFunction(pointCountU, (i) => {
        return deBoor(
          sliceStep(this.points, i, this.points.length, pointCountU, 1),
          this.degreeV,
          this.knotsV,
          v,
        )
      }),
      this.degreeU,
      this.knotsU,
      this.uMin,
      this.uMax,
    )
  }

  debugInfo() {
    const pointCountU = this.knotsU.length - 1 - this.degreeU
    const pointCountV = this.knotsV.length - 1 - this.degreeV
    const grid = []
    for (let u = 0; u < pointCountU; u++) {
      for (let v = 0; v < pointCountV; v++) {
        const i = v * pointCountU + u
        if (u < pointCountU - 1) {
          const j = v * pointCountU + u + 1
          grid.push(this.points[i].p3(), this.points[j].p3())
        }
        if (v < pointCountV - 1) {
          const j = (v + 1) * pointCountU + u
          grid.push(this.points[i].p3(), this.points[j].p3())
        }
      }
    }
    return { points: this.points.map((p) => p.p3()), lines: grid }
  }

  flipped() {
    const pointCountU = this.knotsU.length - 1 - this.degreeU
    return new NURBSSurface(
      arrayFromFunction(this.points.length, (i) => {
        const u = i % pointCountU
        return this.points[i - u + (pointCountU - u - 1)]
      }),
      this.knotsU.map((x) => -x).reverse(),
      this.knotsV,
      this.degreeU,
      this.degreeV,
      -this.uMax,
      -this.uMin,
      this.vMin,
      this.vMax,
    ) as this
  }

  isCoplanarTo(surface: Surface): boolean {
    throw new Error("not implemented")
  }

  isTsForLine(line: L3): number[] {
    // intersect line with
    let t = 1

    throw new Error("not implemented")
  }

  pointFoot(pWC: V3, startU?: number, startV?: number): V3 {
    const closestPointIndex = indexWithMax(
      this.points,
      (p) => -p.p3().distanceTo(pWC),
    )
    const pointCountU = this.knotsU.length - this.degreeU - 1
    const closestPointPos = V(
      closestPointIndex % pointCountU,
      (closestPointIndex / pointCountU) | 0,
    )
    const start = this.guessUVForMeshPos(closestPointPos.x, closestPointPos.y)
    const dpdu = this.dpdu()
    const dpdv = this.dpdv()

    const [u, v] = newtonIterate(
      ([u, v]) => {
        const pUV = this.pUV(u, v)
        const pUVToPWC = pUV.to(pWC)
        return [pUVToPWC.dot(dpdu(u, v)), pUVToPWC.dot(dpdv(u, v))]
      },
      [start.x, start.y],
      8,
      undefined,
      0.1,
    )
    return new V3(u, v, 0)
    /**
     *
     */
    throw new Error("Method not implemented.")
  }

  isCurvesWithPlane(plane: P3): Curve[] {
    throw new Error("Method not implemented.")
  }

  containsPoint(pWC: V3): boolean {
    throw new Error("Method not implemented.")
  }

  loopContainsPoint(contour: Edge[], point: V3): PointVsFace {
    throw new Error("Method not implemented.")
  }

  guessUVForMeshPos(x: int, y: int): V3 {
    function eLerp<T>(
      arr: ArrayLike<T>,
      t: number,
      lerp: (a: T, b: T, t: number) => T,
    ): T {
      if (0 === t % 1) return arr[t]
      return lerp(arr[floor(t)], arr[ceil(t)], t % 1)
    }

    return new V3(
      clamp(
        eLerp(this.knotsU, x + (this.degreeU + 1) / 2, lerp),
        this.uMin,
        this.uMax,
      ),
      clamp(
        eLerp(this.knotsV, y + (this.degreeV + 1) / 2, lerp),
        this.vMin,
        this.vMax,
      ),
      0,
    )
  }
}

NURBSSurface.prototype.uStep = 1 / 8
NURBSSurface.prototype.vStep = 1 / 8

function getInterval(
  degree: int,
  knots: ReadonlyArray<number>,
  t: number,
): int {
  for (let s = degree; s < knots.length - 1 - degree; s++) {
    if (t >= knots[s] && t <= knots[s + 1]) {
      return s
    }
  }
  throw new Error(t + " " + knots)
}

function deBoor(
  points: ReadonlyArray<Vector>,
  degree: int,
  knots: ReadonlyArray<number>,
  t: number,
): Vector {
  // find s (the spline segment) for the [t] value provided
  const s = getInterval(degree, knots, t)

  const v = Vector.pack(points, new Float64Array(points.length * 4))

  // l (level) goes from 1 to the curve degree + 1
  for (let l = 1; l <= degree; l++) {
    // build level l of the pyramid
    for (let i = s; i > s - degree - 1 + l; i--) {
      const alpha = (t - knots[i]) / (knots[i + degree + 1 - l] - knots[i])

      // interpolate each component
      for (let d = 0; d < 4; d++) {
        v[i * 4 + d] = (1 - alpha) * v[(i - 1) * 4 + d] + alpha * v[i * 4 + d]
      }
    }
  }

  return new Vector(v.slice(s * 4, s * 4 + 4))
}
