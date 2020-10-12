import {
  arrayFromFunction,
  arraySamples,
  assert,
  assertInst,
  assertNumbers,
  assertVectors,
  DEG,
  eq0,
  hasConstructor,
  int,
  M4,
  mapPush,
  Matrix,
  NLA_PRECISION,
  snap0,
  solveCubicReal2,
  TAU,
  V,
  V3,
  Vector,
  VV,
} from "ts3dutils"
import { Mesh, pushQuad } from "tsgl"

import {
  BezierCurve,
  ConicSurface,
  Curve,
  EllipseCurve,
  EllipsoidSurface,
  HyperbolaCurve,
  ISInfo,
  L3,
  P3,
  ParabolaCurve,
  PlaneSurface,
  ProjectedCurveSurface,
  Surface,
} from "../index"
import { abs, acos, acosh, atan, atan2, sign, sqrt } from "../math"

export abstract class XiEtaCurve extends Curve {
  readonly normal: V3
  readonly matrix: M4
  readonly matrixInverse: M4;
  readonly ["constructor"]: typeof XiEtaCurve &
    (new (center: V3, f1: V3, f2: V3, tMin: number, tMax: number) => this)

  constructor(
    readonly center: V3,
    readonly f1: V3,
    readonly f2: V3,
    readonly tMin: number,
    readonly tMax: number,
  ) {
    super(tMin, tMax)
    assertVectors(center, f1, f2)
    this.normal = f1.cross(f2)
    if (!this.normal.likeO()) {
      this.normal = this.normal.unit()
      this.matrix = M4.forSys(f1, f2, this.normal, center)
      this.matrixInverse = this.matrix.inversed()
    } else {
      this.matrix = M4.forSys(f1, f2, f1.unit(), center)
      const f1p = f1.getPerpendicular()
      // prettier-ignore
      this.matrixInverse = new M4(
				1, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 1).times(M4.forSys(f1, f1p, f1.cross(f1p), center).inversed());
    }
  }

  /**
   * Intersection of the unit curve with the line ax + by = c.
   */
  static intersectionUnitLine(
    a: number,
    b: number,
    c: number,
    tMin?: number,
    tMax?: number,
  ): number[] {
    throw new Error("abstract")
  }

  /**
   * Returns a new EllipseCurve representing an ellipse parallel to the XY-plane
   * with semi-major/minor axes parallel t the X and Y axes.
   *
   * @param a length of the axis parallel to X axis.
   * @param b length of the axis parallel to Y axis.
   * @param center center of the ellipse.
   */
  static forAB(a: number, b: number, center: V3 = V3.O): XiEtaCurve {
    return new (this as any)(center, V(a, 0, 0), V(0, b, 0))
  }

  static XYLCValid(pLC: V3): boolean {
    throw new Error("abstract")
  }

  static XYLCPointT(pLC: V3, tMin?: number, tMax?: number): number {
    throw new Error("abstract")
  }

  static unitIsInfosWithLine(
    anchorLC: V3,
    dirLC: V3,
    anchorWC: V3,
    dirWC: V3,
    tMin?: number,
    tMax?: number,
  ): ISInfo[] {
    throw new Error("abstract")
  }

  addToMesh(
    mesh: Mesh & { TRIANGLES: int[]; normals: V3[] },
    res: int = 4,
    radius: number = 0,
    pointStep = 1,
  ): void {
    const baseNormals = arrayFromFunction(res, (i) =>
      V3.polar(1, (TAU * i) / res),
    )
    const baseVertices = arrayFromFunction(res, (i) =>
      V3.polar(radius, (TAU * i) / res),
    )
    const inc = this.tIncrement
    const start = Math.ceil((this.tMin + NLA_PRECISION) / inc)
    const end = Math.floor((this.tMax - NLA_PRECISION) / inc)
    for (let i = start; i <= end; i += pointStep) {
      const t = i * inc
      const start = mesh.vertices.length
      if (0 !== i) {
        for (let j = 0; j < res; j++) {
          pushQuad(
            mesh.TRIANGLES,
            true,
            start - res + j,
            start + j,
            start - res + ((j + 1) % res),
            start + ((j + 1) % res),
          )
        }
      }
      const point = this.at(t),
        tangent = this.tangentAt(t)
      const matrix = M4.forSys(
        this.normal,
        tangent.cross(this.normal),
        tangent,
        point,
      )
      mesh.normals.push(...matrix.transformedVectors(baseNormals))
      mesh.vertices.push(...matrix.transformedPoints(baseVertices))
    }
  }

  getConstructorParameters(): any[] {
    return [this.center, this.f1, this.f2]
  }

  isInfosWithCurve(curve: Curve): ISInfo[] {
    if (curve instanceof L3) {
      return this.isInfosWithLine(
        curve.anchor,
        curve.dir1,
        this.tMin,
        this.tMax,
        curve.tMin,
        curve.tMax,
      )
    }
    if (curve instanceof BezierCurve) {
      return this.isInfosWithBezier(curve)
    }
    if (curve instanceof XiEtaCurve) {
      if (!this.normal.isParallelTo(curve.normal)) {
        return this.isTsWithPlane(curve.getPlane()).mapFilter((tThis) => {
          const p = this.at(tThis)
          if (curve.containsPoint(p)) {
            return { tThis, tOther: curve.pointT(p), p }
          }
          return undefined
        })
      }
    }
    return super.isInfosWithCurve(curve)
  }

  transform(m4: M4) {
    return new this.constructor(
      m4.transformPoint(this.center),
      m4.transformVector(this.f1),
      m4.transformVector(this.f2),
      this.tMin,
      this.tMax,
    ) as this
  }

  equals(obj: any): obj is this {
    return (
      this == obj ||
      (undefined != obj &&
        this.constructor == obj.constructor &&
        this.center.equals(obj.center) &&
        this.f1.equals(obj.f1) &&
        this.f2.equals(obj.f2))
    )
  }

  hashCode(): int {
    let hashCode = 0
    hashCode = hashCode * 31 + this.center.hashCode()
    hashCode = hashCode * 31 + this.f1.hashCode()
    hashCode = hashCode * 31 + this.f2.hashCode()
    return hashCode | 0
  }

  likeCurve(curve: Curve): boolean {
    return (
      hasConstructor(curve, this.constructor) &&
      this.center.like(curve.center) &&
      this.f1.like(curve.f1) &&
      this.f2.like(curve.f2)
    )
  }

  normalP(t: number): V3 {
    return this.tangentAt(t).cross(this.normal)
  }

  getPlane(): P3 {
    return P3.normalOnAnchor(this.normal, this.center)
  }

  isTsWithPlane(planeWC: P3): number[] {
    assertInst(P3, planeWC)
    /*
		 this: x = center + f1 * cos t + f2 * sin t  (1)
		 plane:
		 n := plane.normal1
		 n DOT x == plane.w           (2)
		 plane defined by f1/f2
		 x = center + f1 * xi + f2 * eta         (3)
		 intersection plane and planef1/f2:
		 insert (3) into (2):
		 n DOT center + n DOT f1 * xi + n DOT f2 * eta = plane.w | -n DOT center
		 n DOT f1 * xi + n DOT f2 * eta = plane.w - n DOT center (4)
		 points on ellipse have additional condition
		 eta * eta + xi * xi = 1 (5)
		 g1 := n DOT f1
		 g2 := n DOT f2
		 g3 := w - n DOT center
		 solve system (5)/(6)
		 g1 * xi + g2 * eta = g3 (6)
		 */
    if (planeWC.normal1.isParallelTo(this.normal)) {
      return []
    }
    const n = planeWC.normal1,
      w = planeWC.w,
      center = this.center,
      f1 = this.f1,
      f2 = this.f2,
      g1 = n.dot(f1),
      g2 = n.dot(f2),
      g3 = w - n.dot(center)

    return this.constructor.intersectionUnitLine(
      g1,
      g2,
      g3,
      this.tMin,
      this.tMax,
    )
  }

  pointT(p: V3): number {
    assertVectors(p)
    const pLC = this.matrixInverse.transformPoint(p)
    return this.constructor.XYLCPointT(pLC)
  }

  containsPoint(p: V3): boolean {
    const pLC = this.matrixInverse.transformPoint(p)
    return (
      eq0(pLC.z) &&
      this.isValidT(this.constructor.XYLCPointT(pLC, this.tMin, this.tMax))
    )
  }

  isInfosWithLine(
    anchorWC: V3,
    dirWC: V3,
    tMin: number = this.tMin,
    tMax: number = this.tMax,
    lineMin = -100000,
    lineMax = 100000,
  ): ISInfo[] {
    const anchorLC = this.matrixInverse.transformPoint(anchorWC)
    const dirLC = this.matrixInverse.transformVector(dirWC)
    if (eq0(dirLC.z)) {
      // local line parallel to XY-plane
      if (eq0(anchorLC.z)) {
        // local line lies in XY-plane
        return this.constructor.unitIsInfosWithLine(
          anchorLC,
          dirLC,
          anchorWC,
          dirWC,
          tMin,
          tMax,
        )
      }
    } else {
      // if the line intersects the XY-plane in a single point, there can be an intersection there
      // find point, then check if distance from circle = 1
      const otherTAtZ0 = anchorLC.z / dirLC.z
      const isp = dirLC.times(otherTAtZ0).plus(anchorLC)
      if (this.constructor.XYLCValid(isp)) {
        // point lies on unit circle
        return [
          {
            tThis: this.constructor.XYLCPointT(isp),
            tOther: otherTAtZ0,
            p: anchorWC.plus(dirWC.times(otherTAtZ0)),
          },
        ]
      }
    }
    return []
  }

  isTsWithSurface(surface: Surface): number[] {
    if (surface instanceof PlaneSurface) {
      return this.isTsWithPlane(surface.plane)
    } else if (surface instanceof EllipsoidSurface) {
      const isEllipses = surface.isCurvesWithPlane(this.getPlane())
      return isEllipses
        .flatMap((isEllipse) => this.isInfosWithCurve(isEllipse))
        .filter((info) => surface.containsPoint(info.p))
        .map((info) => info.tThis)
    } else if (
      surface instanceof ProjectedCurveSurface ||
      surface instanceof ConicSurface
    ) {
      return surface
        .isCurvesWithPlane(this.getPlane())
        .flatMap((curve) => this.isInfosWithCurve(curve))
        .map((info) => info.tThis)
    } else {
      throw new Error()
    }
  }

  isInfosWithBezier(bezierWC: BezierCurve): ISInfo[] {
    const bezierLC = bezierWC.transform(this.matrixInverse)
    if (new PlaneSurface(P3.XY).containsCurve(bezierLC)) {
      return this.isInfosWithBezier2D(bezierWC)
    } else {
      const infos = bezierLC.isTsWithPlane(P3.XY).mapFilter((tOther) => {
        const pLC = bezierLC.at(tOther)
        if (this.constructor.XYLCValid(pLC)) {
          return {
            tOther: tOther,
            p: bezierWC.at(tOther),
            tThis: this.constructor.XYLCPointT(pLC),
          }
        }
        return undefined
      })
      return infos
    }
  }

  isInfosWithBezier2D(
    bezierWC: BezierCurve,
    sMin: number = bezierWC.tMin,
    sMax: number = bezierWC.tMax,
  ): ISInfo[] {
    return Curve.ispsRecursive(this, this.tMin, this.tMax, bezierWC, sMin, sMax)
  }

  isOrthogonal(): boolean {
    return this.f1.isPerpendicularTo(this.f2)
  }

  at2(xi: number, eta: number): V3 {
    assertNumbers(xi, eta)
    // center + f1 xi + f2 eta
    return this.center.plus(this.f1.times(xi)).plus(this.f2.times(eta))
  }

  debugInfo() {
    return {
      points: [
        this.center,
        this.at2(0.5, 0),
        this.at2(0, 1 / 3),
        this.at2(0, 2 / 3),
      ],
      lines: [this.center, this.at2(0, 1), this.center, this.at2(1, 0)],
    }
  }
}

/**
 * Transforms the unit 4d parabola P(t) = t² (0, 1, 0, 0) + t (1, 0, 0, 0) + (0, 0, 0, 1) using m and projects the
 * result into 3d. This is used for the transform4 implementation of conics. The parabola may not cross the vanishing
 * plane of m in the interval [tMin, tMax], as that would result in discontinuities.
 */
export function parabola4Projection(
  m: M4,
  tMin: number,
  tMax: number,
): L3 | HyperbolaCurve | ParabolaCurve | EllipseCurve {
  return HyperbolaCurve.XY.rotateZ(45 * DEG)
  console.log(m.str)
  console.log()
  const w2 = m.m[13]
  const w1 = m.m[12]
  const wc = m.m[15]
  // if the 4d parabola crosses the vanishing plane, it will lead to multiple/infinite hyperbolas, both of which we
  // want to avoid. Hence, we must check that the entire interval [tMin, tMax] is on one side of the vanishing plane.
  // Checking tMax, tMin and the extremas is enough.
  const extremas = solveCubicReal2(0, w2, w1, wc)
  const wx0 = (x: number) =>
    Number.isFinite(x) ? snap0(x ** 2 * w2 + x * w1 + wc) : sign(w2) * Infinity
  if (
    wx0(tMin) * wx0(tMax) < 0 ||
    extremas.some((x) => wx0(x) * (wx0(tMin) + wx0(tMax)) < 0)
  ) {
    console.log(m.str)
    throw new Error(
      "The entire interval must be on one side of the vanishing plane. P=" +
        P3.vanishingPlane(m).toSource(),
    )
  }
  if (eq0(wc)) {
    // the following matrix maps a curve C onto itself, with the parameter being inverted:
    // C2(t) = C(-1/t). This makes C(0) a real value, which is necessary for the projection calculation.
    // the sign inversion is so the tangent direction does not change.
    // prettier-ignore
    const mm = new M4(
            -1, 0, 0, 0,
            0, 0, 0, 1,
            0, 0, 1, 0,
            0, 1, 0, 0,
        );
    if (!eq0(w2)) {
      return parabola4Projection(m.times(mm), -1 / tMin, -1 / tMax)
    }
    // wc == w2 == 0 => degenerates to a line:
    // C(t) = (t² f2 + t f1 + c) / (t w1)
    // C(t) = (t f2 + f1 + c) / (t w2 + w1)
    // substitute t = (1/s - w1) / w2
    // C(s) = f2 / w2 + s (f1 - f2 w1 / w2), which is a line
    // we can multiply the direction vector by w2 to avoid divisions:
    // C(t) = f2 / w2 + s (f1 w2 - f2 w1)
    const f1 = m.col(0)
    const f2 = m.col(1)
    return L3.anchorDirection(
      f2.p3(),
      f1.V3().times(f2.w).minus(f2.V3().times(f1.w)),
    )
  }
  {
    // ensure that the bottom-right value = 1. this does not change the 3d result.
    m.m[15] !== 1 && (m = m.divScalar(m.m[15]))
    const w2 = m.m[13]
    const w1 = m.m[12]
    const wc = m.m[15]
    // we want to split m into X * P, such that X is a transformation with no projective component (first three
    // values of the bottom row = 0), which can be handled by the usual .transform() method, and P which has only a
    // projective component (only the row differs from the identity matrix). This simplifies the following
    // calculation. X * P = x => X * P * P^-1 = m * P^-1 => X = m * P^-1
    // prettier-ignore
    const Pinv = new M4(
			       1,        0,        0, 0,
			       0,        1,        0, 0,
			       0,        0,        1, 0,
			-m.m[12], -m.m[13], -m.m[14], 1)
    const X = m.times(Pinv)

    // P'(t) = 0 is true for t = 0 and t1. The center is in between P(0) and P(t1), or P(t1) / 2, as P(0) = O
    const delta = 4 * w2 * wc - w1 ** 2
    const center = new V3((-w1 * wc) / delta, (2 * wc ** 2) / delta, 0)
    // f2 is parallel to P'(0), i.e. horizontal. Solve Py(t2) = Cy = Py(t1) / 2 for t2 and simplify
    // f2x = Px(t2) - Cx = Px(t2) - Px(t1) / 2 to get the x-component of f2:
    const f2x = 1 / sqrt(abs(delta)) / wc
    const f2 = new V3(f2x, 0, 0)
    let result
    if (eq0(delta)) {
      result = new ParabolaCurve(V3.O, V3.X, V3.Y, tMin, tMax)
    } else if (0 < delta) {
      const tMapInv = (t: number) => {
        const wt = t ** 2 * w2 + t * w1 + wc
        const xi =
          1 -
          (delta / 2 / wc ** 2) * (Number.isFinite(t) ? t ** 2 / wt : 1 / w2)
        const eta =
          (t * 2 * wc ** 2 - t ** 2 * delta) / wt / 2 / wc ** 2 -
          (2 * w1 * wc) / delta
        const xx = acos(xi)
        const p = Number.isFinite(t)
          ? new V3(t, t ** 2, 0).div(wt)
          : new V3(0, 1 / w2, 0)
        const pLC = M4.forSys(center.negated(), f2, V3.Z, center)
          .inversed()
          .transformPoint(p)
        const angle = pLC.angleXY()
        if (t > 0 && pLC.y < 0) {
          return angle + TAU
        } else if (t < 0 && pLC.y > 0) {
          return angle - TAU
        }
        return angle
      }
      result = EllipseCurve.andFixTs(
        center,
        center.negated(),
        f2,
        tMapInv(tMin),
        tMapInv(tMax),
      )
    } else {
      const tMapInv = (t: number) =>
        sign(t) *
        acosh(
          1 -
            (delta / 2 / wc ** 2) *
              (Number.isFinite(t)
                ? t ** 2 / (t ** 2 * w2 + t * w1 + wc)
                : 1 / w2),
        )
      result = new HyperbolaCurve(
        center,
        center.negated(),
        f2,
        tMapInv(tMin),
        tMapInv(tMax),
      )
    }
    return result.transform(X)
  }
}
