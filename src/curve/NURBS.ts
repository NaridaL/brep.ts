import {
  arrayFromFunction,
  arraySamples,
  assert,
  between,
  DEG,
  eq,
  eq0,
  hasConstructor,
  int,
  le,
  lerp,
  M4,
  MINUS,
  newtonIterateWithDerivative2,
  NLA_DEBUG,
  snap,
  snap0,
  Tuple2,
  Tuple3,
  V,
  V3,
  vArrGet,
  Vector,
  VV,
} from "ts3dutils"

import {
  BezierCurve,
  Curve,
  EllipseCurve,
  HyperbolaCurve,
  L3,
  P3,
  ParabolaCurve,
} from "../index"
import { abs, cos, cosh, PI, sin, sinh, sqrt, SQRT1_2 } from "../math"
import { ISInfo } from "./Curve"

/**
 * Non-Uniform Rational B-Spline implementation.
 *
 * See https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/ for a good reference.
 *
 *
 */
export class NURBS extends Curve {
  constructor(
    /**
     * The control points of the NURBS curve, as 4D homogeneous coordinates.
     */
    readonly points: ReadonlyArray<Vector>,
    /**
     * The degree of the NURBS curve. Must be at least 1 (linear).
     */
    readonly degree: int,
    readonly knots: ReadonlyArray<number> = NURBS.openUniformKnots(
      points.length,
      degree,
    ),
    tMin: number = knots[degree],
    tMax: number = knots[knots.length - degree - 1],
  ) {
    super(tMin, tMax)
    const knotsLength = points.length + degree + 1
    NLA_DEBUG && Object.freeze(points)
    NLA_DEBUG && Object.freeze(knots)

    assert(
      knots.length === knotsLength,
      "bad knot vector length: expected " +
        knotsLength +
        " (degree = " +
        degree +
        " pcount = " +
        points.length +
        "), but was " +
        knots.length,
    )
    assert(knots[degree] <= tMin)
    assert(tMax <= knots[knots.length - degree - 1])
    for (let i = 0; i < points.length; i++) {
      assert(points[i].dim() == 4)
    }
    assert(degree >= 1, "degree must be at least 1 (linear)")
    assert(degree % 1 == 0)
    assert(
      -1 == knots.firstUnsorted(MINUS),
      "knot values must be in ascending order",
    )
  }

  getConstructorParameters() {
    return [this.points, this.degree, this.knots]
  }

  at4(t: number) {
    NLA_DEBUG && assert(between(t, this.tMin, this.tMax), t)
    const { points, degree, knots } = this

    // find s (the spline segment) for the [t] value provided
    const s = this.tInterval(t)

    const v = Vector.pack(
      points,
      new Float64Array((degree + 1) * 4),
      s - degree,
      0,
      degree + 1,
    )

    for (let level = 0; level < degree; level++) {
      // build level l of the pyramid
      for (let i = degree; i > level; i--) {
        const alpha =
          (t - knots[i + s - degree]) /
          (knots[i + s - level] - knots[i + s - degree])

        // interpolate each component
        for (let dim = 0; dim < 4; dim++) {
          v[i * 4 + dim] =
            (1 - alpha) * v[(i - 1) * 4 + dim] + alpha * v[i * 4 + dim]
        }
      }
    }

    return new Vector(v.slice(degree * 4, (degree + 1) * 4))
  }

  at(t: number) {
    return this.at4(t).p3()
  }

  /*
	d(k, i, t) = a(i, k, t) * d(k - 1, i, t) + (1 - a(i, k, t)) * d(k - 1, i - 1, t)
	a(i, k, t) = (t - knots[i]) / (knots[i + 1 + n - k] - knots[i])
	a'(i, k, t) = 1 / (knots[i + 1 + n - k] - knots[i])

	d/dt =  a(i, k, t) * d'(k - 1, i, t) + a'(i, k, t) * d(k - 1, i, t)
	+ (1 - a(i, k, t)) * d'(k - 1, i - 1, t) + a'(i, k, t) * d(k - 1, i - 1, t)
*/
  ptDtDdt4(t: number) {
    const { points, degree, knots } = this

    // find s (the spline segment) for the [t] value provided
    const s = this.tInterval(t)

    const v = Vector.pack(
      points,
      new Float64Array((degree + 1) * 4),
      s - degree,
      0,
      degree + 1,
    )

    let ddt: Vector = Vector.Zero(4),
      derivative: Vector

    for (let level = 0; level < degree; level++) {
      if (level == degree - 2) {
        // see https://www.globalspec.com/reference/61012/203279/10-8-derivatives
        const a = new Vector(v.slice(degree * 4, (degree + 1) * 4))
        const b = new Vector(v.slice((degree - 1) * 4, degree * 4))
        const c = new Vector(v.slice((degree - 2) * 4, (degree - 1) * 4))
        function step(
          k: int,
          i: int,
          dkMinus1iMinus1: Vector,
          dkMinus1i: Vector,
        ) {
          return dkMinus1i
            .minus(dkMinus1iMinus1)
            .times(k / (knots[i + degree - k] - knots[i - 1]))
        }
        ddt = step(
          degree,
          s + 1,
          step(degree - 1, s + 1, a, b),
          step(degree - 1, s, b, c),
        )
      }
      if (level == degree - 1) {
        const a = new Vector(v.slice(degree * 4, (degree + 1) * 4))
        const b = new Vector(v.slice((degree - 1) * 4, degree * 4))
        derivative = b!.minus(a!).times(degree / (knots[s] - knots[s + 1]))
      }
      for (let i = degree; i > level; i--) {
        const alpha =
          (t - knots[i + s - degree]) /
          (knots[i + s - level] - knots[i + s - degree])

        // interpolate each component
        for (let dim = 0; dim < 4; dim++) {
          v[i * 4 + dim] =
            (1 - alpha) * v[(i - 1) * 4 + dim] + alpha * v[i * 4 + dim]
        }
      }
    }

    const p = new Vector(v.slice(degree * 4, degree * 4 + 4))

    return [p, derivative!, ddt!]
  }

  tangentAt(t: number) {
    // x(t) = xw(t) / w(t)
    // quotient rule

    const [p, derivative] = this.ptDtDdt4(t)

    const expected = derivative
      .times(p.w)
      .minus(p.times(derivative.w))
      .div(p.w ** 2)
      .V3()

    return expected
  }

  ddt(t: number) {
    const [p, dt, ddt] = this.ptDtDdt4(t)
    // =(-w(t) x(t) w''(t) - 2 w(t) w'(t) x'(t) + 2 x(t) w'(t)^2 + w(t)^2 x''(t))/w(t)^3
    // =(x(t) ((-w(t)) w''(t) + 2 w'(t)^2) - x'(t) 2 w(t) w'(t) + x''(t) w(t)^2 )/w(t)^3
    // prettier-ignore
    return Vector.add(
		p.times(-p.w * ddt.w + 2 * dt.w ** 2),
		dt.times(-2 * p.w * dt.w),
		ddt.times(p.w ** 2)
	).div(p.w ** 3).V3();
  }

  ptDtDdt(t: number) {
    const [pt, dt4, ddt4] = this.ptDtDdt4(t)
    return [
      pt.p3(),
      dt4
        .times(pt.w)
        .minus(pt.times(dt4.w))
        .div(pt.w ** 2)
        .V3(),
      Vector.add(
        pt.times(-pt.w * ddt4.w + 2 * dt4.w ** 2), //
        dt4.times(-2 * pt.w * dt4.w),
        ddt4.times(pt.w ** 2),
      )
        .div(pt.w ** 3)
        .V3(),
    ]
  }

  pointT(pWC: V3) {
    return this.closestTToPoint(pWC)
  }

  closestTToPoint(
    p: V3,
    tStart?: number,
    tMin = this.tMin,
    tMax = this.tMax,
  ): number {
    // this.at(t) has minimal distance to p when this.tangentAt(t) is perpendicular to
    // the vector between this.at(t) and p. This is the case iff the dot product of the two is 0.
    // f = (this.at(t) - p) . (this.tangentAt(t)
    // df = this.tangentAt(t) . this.tangentAt(t) + (this.at(t) - p) . this.ddt(t)
    //    = this.tangentAt(t)² + (this.at(t) - p) . this.ddt(t)
    const f = (t: number): Tuple2<number> => {
      const [pt, dt, ddt] = this.ptDtDdt(t)
      return [pt.minus(p).dot(dt), dt.squared() + pt.minus(p).dot(ddt)]
    }
    //checkDerivate(f, df, tMin, tMax)

    const STEPS = 32
    if (undefined === tStart) {
      tStart = arraySamples(tMin, tMax, STEPS).withMax(
        (t) => -this.at(t).distanceTo(p),
      )
    }

    const result = newtonIterateWithDerivative2(
      f,
      tStart,
      8,
      this.tMin,
      this.tMax,
    )
    //assert(undefined !== result)
    return result!
  }

  containsPoint(pWC: V3) {
    const tGuess = this.closestTToPoint(pWC)
    return undefined === tGuess ? false : this.at(tGuess).like(pWC)
  }

  derivate(): NURBS {
    const k = this.degree
    const ps = arrayFromFunction(this.points.length - 1, (i) =>
      this.points[i]
        .to(this.points[i + 1])
        .times(k / (this.knots[i + k + 1] - this.knots[i + 1])),
    )
    return new NURBS(
      ps,
      this.degree - 1,
      this.knots.slice(1, -1),
      this.tMin,
      this.tMax,
    )
  }

  /**
   * Create a new NURBS of equal degree with the added knot [newKnot]. New NURBS will have one additional control
   * point.
   */
  withKnot(newKnot: number) {
    assert(between(newKnot, this.tMin, this.tMax))
    const k = this.tInterval(newKnot)
    const { knots, points, degree } = this
    const insertPoints = arrayFromFunction(this.degree, (j) => {
      const i = k - degree + 1 + j
      const aiNumerator = newKnot - knots[i]
      // 0/0 defined as 0:
      const ai =
        aiNumerator == 0 ? 0 : aiNumerator / (knots[i + degree] - knots[i])
      assert(between(ai, 0, 1))
      return Vector.lerp(points[i - 1], points[i], ai)
    })
    const newPoints = points.slice()
    newPoints.splice(k - degree + 1, degree - 1, ...insertPoints)
    const newKnots = knots.slice()
    newKnots.splice(k + 1, 0, newKnot)
    return new NURBS(newPoints, degree, newKnots, this.tMin, this.tMax)
  }

  removeKnot(t: number) {
    const { knots, points, degree } = this
    let k = this.tInterval(t),
      s = 0 // s = multiplicity of the knot
    while (knots[k + 1] == t) {
      k++
      s++
    }
    if (s == 0) throw new Error("There is no knot " + t + "!")
    // the points which were relevant when inserting were (k - p - 1) to (k - 1). (- 1) because the current k has
    // been increased by one due to the insertion.
    // p - 1 points were replaced by p points, hence we need to generate the original p - 1 point, + 1 to check if
    // this transformation is valid.
    const insertPoints = [points[k - degree - 1]]
    const oldKnots = knots.slice()
    oldKnots.splice(k, 1)
    for (let i = k - degree; i <= k - s; i++) {
      const alphaInv = (oldKnots[i + degree] - oldKnots[i]) / (t - oldKnots[i])
      const oldPoint = Vector.lerp(insertPoints.last, points[i], alphaInv)
      insertPoints.push(oldPoint)
    }
    if (insertPoints.last.like(points[k + 1 - s])) {
      const oldPoints = points.slice()
      oldPoints.splice(k - degree - 1, degree - s + 3, ...insertPoints)
      return new NURBS(oldPoints, degree, oldKnots)
    }
    return undefined
  }

  static openUniformKnots(pointCount: int, degree: int, tMin = 0, tMax = 1) {
    const knotsLength = pointCount + degree + 1
    return arrayFromFunction(knotsLength, (i) => {
      if (i <= degree) {
        return tMin
      } else if (i >= knotsLength - degree - 1) {
        return tMax
      } else {
        return lerp(tMin, tMax, (i - degree) / (knotsLength - degree * 2 - 1))
      }
    })
  }

  static bezierKnots(
    degree: int,
    tMin: number = 0,
    tMax: number = 1,
  ): number[] {
    const result = new Array((degree + 1) * 2)
    for (let i = 0; i < degree + 1; i++) {
      result[i] = tMin
      result[degree + 1 + i] = tMax
    }
    return result
  }

  static fromBezier(bezier: BezierCurve) {
    const bezier01 = bezier.selectPart(bezier.tMin, bezier.tMax)
    return NURBS.Bezier(bezier01.points)
  }

  static Bezier(points: (V3 | Vector)[], tMin = 0, tMax = 1) {
    return new NURBS(
      points.map((p) =>
        p instanceof V3 ? new Vector(new Float64Array([p.x, p.y, p.z, 1])) : p,
      ),
      points.length - 1,
      arrayFromFunction(points.length * 2, (i) => (i < points.length ? 0 : 1)),
      tMin,
      tMax,
    )
  }

  static fromHyperbola(
    hyperbola: HyperbolaCurve,
    tMin = hyperbola.tMin,
    tMax = hyperbola.tMax,
  ) {
    const p0 = HyperbolaCurve.XY.at(tMin)
    const p2 = HyperbolaCurve.XY.at(tMax)
    const p1 = new V3(
      (sinh(tMin) - sinh(tMax)) / sinh(tMin - tMax),
      (cosh(tMin) - cosh(tMax)) / sinh(tMin - tMax),
      0,
    )
    // M: midpoint between p0 and p2
    // X: intersection of line through p1 and M and unit hyperbola
    // result.at(1/2) = X
    // result.at(1/2) = (1/4 p0 + 1/2 p1 w + 1/4 p2) / (1/4 + 1/ 2 w + 1/4)
    // result.at(1/2) = (1/2 p0 + p1 w + 1/2 p2) / (1 + w)
    // result.at(1/2) = (M + p1 w) / (1 + w) = X
    // => w * (p1 - X) = (X - M)
    // as p1, X and M are all on the same line, we can solve this equation with only the x
    const M = p0.lerp(p2, 0.5)
    const Xx = 1 / sqrt(1 - (M.y / M.x) ** 2)
    const w = (Xx - M.x) / (p1.x - Xx)
    return NURBS.fromV3s([p0, p1, p2], 2, undefined, [1, w, 1]).transform(
      hyperbola.matrix,
    )
  }

  static fromParabola(parabola: ParabolaCurve) {
    return NURBS.fromBezier(parabola.asBezier())
  }

  static fromEllipse(ellipse: EllipseCurve) {
    const unitSemiEllipse = new NURBS(
      [
        VV(1, 0, 0, 1),
        VV(1, 1, 0, 1).times(SQRT1_2),
        VV(0, 1, 0, 1),
        VV(-1, 1, 0, 1).times(SQRT1_2),
        VV(-1, 0, 0, 1),
        VV(-1, -1, 0, 1).times(SQRT1_2),
        VV(0, -1, 0, 1),
      ],
      2,
      [0, 0, 0, PI / 2, PI / 2, PI, PI, (3 * PI) / 2, (3 * PI) / 2, 2 * PI],
    )
    return unitSemiEllipse.transform(ellipse.matrix)
  }

  /**
   * Create a new NURBS from V3s, with optional weights.
   * @param points
   * @param degree
   * @param knots
   * @param weights
   */
  static fromV3s(
    points: V3[],
    degree: int,
    knots?: number[],
    weights: number[] = arrayFromFunction(points.length, () => 1),
  ) {
    assert(points.length == weights.length)
    return new NURBS(
      points.map((p, i) => Vector.fromV3AndWeight(p, weights[i])),
      degree,
      knots,
    )
  }

  isUniform(precision = 0) {
    const intervals = arrayFromFunction(
      this.knots.length - 1,
      (i) => this.knots[i + 1] - this.knots[i],
    )
    const [min, max] = minAndMax(intervals)
    return eq(min, max, precision)
  }

  /**
   * NURBS is a B spline if control points all have the same weight.
   */
  isBSpline(precision = 0) {
    const [minWeight, maxWeight] = minAndMax(this.points.map((p) => p.w))
    return eq(minWeight, maxWeight, precision)
  }

  /**
   * Whether this is a (rational) bezier curve.
   */
  isBezier(precision = 0) {
    if (this.degree + 1 != this.points.length) return false
    const [min0, max0] = minAndMax(this.knots, 0, this.degree + 1)
    if (!eq(min0, max0, precision)) return false
    const [min1, max1] = minAndMax(this.knots, this.degree + 1)
    if (!eq(min1, max1, precision)) return false
    return true
  }

  /**
   * Splits NURBS curve into rational bezier curves.
   * See https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/subdivision.html
   */
  getSegments() {
    const { knots, points, degree } = this
    const result: NURBS[] = []
    const v = Vector.pack(points, new Float64Array(points.length * 4))
    const vectorFromV = (i: int) => new Vector(v.slice(i * 4, (i + 1) * 4))

    let k = degree + 1 // k = knot index we are duplicating
    while (k < knots.length - degree - 1) {
      const t = knots[k]
      const prevKnot = knots[k - 1]
      let s = 1 // s = multiplicity of the knot
      while (knots[k + 1] == t) {
        k++
        s++
      }
      const newNURBSPoints = new Array(degree + 1)
      // the first s + 1 points are identical to the current curve
      for (let i = 0; i < s + 1; i++) {
        newNURBSPoints[i] = vectorFromV(k - degree - s + i)
      }
      // we need to have multiplicity degree, so insert (degree - s) times
      for (let level = 1; level <= degree - s; level++) {
        for (let i = k - degree; i <= k - s - level; i++) {
          const alpha = (t - prevKnot) / (knots[i + degree + 1] - prevKnot)
          for (let dim = 0; dim < 4; dim++) {
            v[i * 4 + dim] =
              (1 - alpha) * v[i * 4 + dim] + alpha * v[(i + 1) * 4 + dim]
          }
        }
        newNURBSPoints[s + level] = vectorFromV(k - degree)
      }
      const newNURBSKnots = arrayFromFunction((degree + 1) * 2, (i) =>
        i < degree + 1 ? knots[k - s] : t,
      )
      result.push(new NURBS(newNURBSPoints, degree, newNURBSKnots))
      k++
    }

    // last curve
    const newNURBSPoints = arrayFromFunction(degree + 1, (i) =>
      vectorFromV(points.length - degree - 1 + i),
    )
    const newNURBSKnots = arrayFromFunction((degree + 1) * 2, (i) =>
      i < degree + 1 ? knots[k - 1] : knots[k],
    )
    result.push(new NURBS(newNURBSPoints, degree, newNURBSKnots))
    return result
  }

  split(t: number) {
    const { knots, points, degree } = this
    assert(le(this.tMin, t) && le(t, this.tMax))
    let k = this.tInterval(t),
      s = 0 // s = multiplicity of the knot
    while (knots[k + 1] == t) {
      k++
      s++
    }
    const vectorFromV = (i: int) => new Vector(v.slice(i * 4, (i + 1) * 4))

    const leftPoints = new Array(k + 1 - s)
    // the first k + s + 1 points are identical to the current curve
    for (let i = 0; i < k + s - degree + 1; i++) {
      leftPoints[i] = this.points[i]
    }
    const rightPointsLength = points.length - (k - degree)
    const v = Vector.pack(
      points,
      new Float64Array(rightPointsLength * 4),
      k - degree,
    )
    // we need to have multiplicity degree, so insert (degree - s) times
    for (let level = 1; level <= degree - s; level++) {
      for (let i = k - degree; i <= k - s - level; i++) {
        const alpha =
          (t - knots[i + level]) / (knots[i + degree + 1] - knots[i + level])
        const j = i - (k - degree)
        for (let dim = 0; dim < 4; dim++) {
          v[j * 4 + dim] =
            (1 - alpha) * v[j * 4 + dim] + alpha * v[(j + 1) * 4 + dim]
        }
      }
      leftPoints[k - degree + level] = vectorFromV(0)
    }
    const leftKnots = knots.slice(0, k + degree + 2 - s)
    for (let i = 0; i < degree - s + 1; i++) {
      leftKnots[k - s + 1 + i] = t
    }
    const rightKnots = knots.slice(k - degree)
    for (let i = 0; i < degree + 1; i++) {
      rightKnots[i] = t
    }

    const rightPoints = arrayFromFunction(rightPointsLength, (i) =>
      vArrGet(v, 4, i),
    )
    return [
      new NURBS(leftPoints, degree, leftKnots),
      new NURBS(rightPoints, degree, rightKnots),
    ]
  }

  simplify() {
    assert(this.isBezier())
    if (3 == this.degree && this.isBSpline()) {
      return new BezierCurve(
        this.points[0].p3(),
        this.points[1].p3(),
        this.points[2].p3(),
        this.points[3].p3(),
        this.tMin,
        this.tMax,
      )
    } else if (2 == this.degree) {
      const [P0, P1, P2] = this.points
      const [p0, p1, p2] = this.points.map((p) => p.p3())
      const c = NURBS.simplifyUnit2(P0.w, P1.w, P2.w).transform(
        M4.forSys(p1.to(p0), p1.to(p2), undefined, p1),
      )
      const [tMin, tMax] = [c.pointT(p0), c.pointT(p2)].sort()
      return c.withBounds(snap(tMin, c.tMin), snap(tMax, c.tMax))
    } else if (1 == this.degree) {
      return L3.throughPoints(this.points[0].p3(), this.points[1].p3())
    } else {
      return this
    }
  }

  static simplifyUnit2(w0: number, w1: number, w2: number) {
    // see https://math.stackexchange.com/a/2794874/230980
    const delta = w0 * w2 - w1 ** 2
    const cxy = (w0 * w2) / 2 / delta
    const center = new V3(cxy, cxy, 0)
    const k = (w1 ** 2 + delta - 2 * w1 * sqrt(abs(delta))) / 2 / delta
    const p = V3.X
    const q = new V3(k, cxy, 0)
    // const q = new V3(cxy, k, 0)
    if (eq0(delta)) {
      return new ParabolaCurve(
        new V3(1 / 4, 1 / 4, 0),
        new V3(1, -1, 0),
        new V3(1, 1, 0),
        -0.5,
        0.5,
      )
    } else if (delta < 0) {
      // hyperbola
      return new HyperbolaCurve(center, center.to(p), center.to(q))
    } else {
      // ellipse
      return new EllipseCurve(center, center.to(p), center.to(q), 0)
    }
  }

  elevateDegreeBezier() {
    assert(this.isBezier())
    const newPoints = new Array(this.points.length + 1)
    newPoints[0] = this.points[0]
    newPoints[this.points.length] = this.points[this.points.length - 1]
    for (let i = 1; i < this.points.length; i++) {
      newPoints[i] = Vector.lerp(
        this.points[i],
        this.points[i - 1],
        i / (this.degree + 1),
      )
    }
    const newKnots = NURBS.bezierKnots(
      this.degree + 1,
      this.knots[0],
      this.knots[this.degree + 1],
    )
    return new NURBS(
      newPoints,
      this.degree + 1,
      newKnots as number[],
      this.tMin,
      this.tMax,
    )
  }

  elevateDegree() {
    const segmentsElevated = this.getSegments().map((b) =>
      b.elevateDegreeBezier(),
    )
    // stitch together the segments
    const newPoints = new Array(2 + segmentsElevated.length * this.degree)
    newPoints[0] = segmentsElevated[0].points[0]
    newPoints.last = segmentsElevated.last.points.last
    for (let i = 0; i < segmentsElevated.length; i++) {
      for (let pi = 1; pi < segmentsElevated[i].points.length - 1; pi++) {
        newPoints[i * (segmentsElevated[0].points.length - 2) + pi] =
          segmentsElevated[i].points[pi]
      }
    }
    const newKnots = new Array(newPoints.length + this.degree + 2)
    for (let i = 0; i < this.degree + 2; i++) {
      newKnots[i] = this.knots[0]
    }
    for (let i = 0; i < segmentsElevated.length; i++) {
      for (let pi = 1; pi < segmentsElevated[i].points.length - 1; pi++) {
        newKnots[
          i * (segmentsElevated[0].points.length - 2) + pi + this.degree + 1
        ] = segmentsElevated[i].knots.last
      }
    }
    newKnots[newKnots.length - 1] = this.knots.last
    newKnots[newKnots.length - 2] = this.knots.last

    let result = new NURBS(
      newPoints,
      this.degree + 1,
      newKnots as number[],
      this.tMin,
      this.tMax,
    )
    for (let i = 0; i < segmentsElevated.length - 1; i++) {
      let optimization
      while (
        (optimization = result.removeKnot(segmentsElevated[i].knots.last))
      ) {
        result = optimization
      }
    }
    return result
  }

  transform(m4: M4) {
    return this.transform4(m4) as this
  }

  transform4(m4: M4) {
    return new NURBS(
      this.points.map((p) => m4.timesVector(p)),
      this.degree,
      this.knots,
      this.tMin,
      this.tMax,
    )
  }

  /**
   * Returns the index of the interval which contains the value t.
   */
  tInterval(t: number) {
    const { degree, knots } = this
    for (let s = degree; s < knots.length - 1 - degree; s++) {
      if (t >= knots[s] && t <= knots[s + 1]) {
        return s
      }
    }
    throw new Error(t + " " + knots)
  }

  static EX2D = NURBS.fromV3s(
    [
      V(51, 141),
      V(11, 76),
      V(29, 32),
      V(46, 102),
      V(74, 148),
      V(189, 107),
      V(56, 10),
      V(206, 10),
      V(211, 98),
      V(195, 141),
      V(139, 148),
    ],
    4,
  )
  static EX3D = new NURBS(
    [
      VV(94, 0, -34, 1),
      VV(69, 57, 45, 0.5),
      VV(-20, 44, 91, 1),
      VV(-89, -13, 47, 0.5),
      VV(-56, -97, -7, 1),
      VV(34, -83, -54, 0.5),
      VV(112, -53, 16, 1),
      VV(79, 30, 70, 0.5),
      VV(-2, -9, 141, 1),
      VV(-80, -40, 72, 0.5),
      VV(-38, -150, 43, 1),
      VV(43, -110, -29, 0.5),
      VV(130, -106, 65, 1),
    ],
    2,
    [-12, -12, -12, -8, -8, -4, -4, 0, 0, 4, 4, 8, 8, 12, 12, 12],
  )

  static UnitCircle(sections: int = 2, tMin: number = 0, tMax: number = PI) {
    const dt = tMax - tMin
    const tStep = dt / sections
    const w = sin(PI / 2 - tStep / 2)
    console.log(tStep / 2 / DEG)
    // cos
    const r = 1 / cos(tStep / 2)
    const points = arrayFromFunction(sections * 2 + 1, (i) => {
      const t = lerp(tMin, tMax, i / 2 / sections)
      if (i % 2 == 0) {
        // control point on circle
        return VV(cos(t), sin(t), 0, 1)
      } else {
        return VV(r * w * cos(t), r * w * sin(t), 0, w)
      }
    })
    const knots = []
    knots.push(tMin, tMin, tMin)
    for (let i = 0; i < sections - 1; i++) {
      const knot = lerp(tMin, tMax, (i + 1) / sections)
      knots.push(knot, knot)
    }
    knots.push(tMax, tMax, tMax)
    return new NURBS(points, 2, knots)
  }

  debugInfo() {
    return {
      points: [
        ...this.knots.slice(this.degree, -this.degree).map((t) => this.at(t)),
        ...this.points.map((p) => p.p3()),
      ],
      lines: this.points.flatMap((p, i, ps) =>
        ps[i + 1] ? [p.p3(), ps[i + 1].p3()] : [],
      ),
    }
  }

  isTsWithPlane(planeWC: P3) {
    const { knots, degree, points } = this
    const controlPointTs = [
      knots[degree],
      ...points
        .slice(1, -1)
        .map((p, i) =>
          this.closestTToPoint(
            p.p3(),
            undefined,
            knots[i + 3],
            knots[i + degree],
          ),
        ),
      knots[knots.length - degree - 1],
    ]
    const result: number[] = []
    for (let i = 0; i < this.points.length - 1; i++) {
      const findClosest = (startT: number) => {
        console.log("startT", startT)
        // try {
        const f = (t: number): Tuple2<number> => {
          const [p, dt] = this.ptDtDdt(t)
          return [planeWC.distanceToPointSigned(p), planeWC.normal1.dot(dt)]
        }
        let t = newtonIterateWithDerivative2(f, startT, 8, this.tMin, this.tMax)
        let [distanceAtT, distanceDtAtT] =
          undefined === t ? [undefined, undefined] : f(t)
        if (t === undefined || !eq0(distanceAtT) || eq0(distanceDtAtT)) {
          t = newtonIterateWithDerivative2(
            (t) => {
              const [, dt, ddt] = this.ptDtDdt(t)
              return [planeWC.normal1.dot(dt), planeWC.normal1.dot(ddt)]
            },
            startT,
            8,
            this.tMin,
            this.tMax,
          )
        }
        ;[distanceAtT, distanceDtAtT] = undefined === t ? [] : f(t)
        if (
          undefined !== t &&
          eq0(distanceAtT) &&
          !result.some((r) => eq(r, t))
        ) {
          result.push(t)
        }
      }
      const a = this.points[i].p3()
      const b = this.points[i + 1].p3()

      const ad = snap0(planeWC.distanceToPointSigned(a))
      const bd = snap0(planeWC.distanceToPointSigned(b))
      if (ad * bd < 0) {
        const startT = lerp(
          controlPointTs[i],
          controlPointTs[i + 1],
          ad / (ad - bd),
        )
        findClosest(startT)
      } else if (0 == bd) {
        findClosest(this.closestTToPoint(b, controlPointTs[i + 1]))
      }
    }
    return result
  }

  isInfosWithCurve(curveWC: Curve) {
    if (curveWC instanceof L3) {
      return this.isInfosWithLine(curveWC.anchor, curveWC.dir1)
    }
    return super.isInfosWithCurve(curveWC)
  }

  isInfosWithLine(anchor: V3, dir: V3): ISInfo[] {
    const thisPlane = P3.fromPoints(this.points.map((p) => p.p3()))!
    const l = L3.anchorDirection(anchor, dir)
    const maxDistanceToPlane = this.points
      .map((p) => thisPlane.distanceToPoint(p.p3()))
      .max()
    const thisIsPlanar = eq0(maxDistanceToPlane)
    if (thisIsPlanar && !thisPlane.containsLine(l)) {
      const [t] = l.isTsWithPlane(thisPlane)
      if (undefined === t) return []

      const p = l.at(t)
      return this.containsPoint(p)
        ? [{ tThis: this.pointT(p), tOther: L3.pointT(anchor, dir, p), p }]
        : []
    } else {
      const thisTs = this.isTsWithPlane(
        P3.normalOnAnchor(thisPlane.normal1.cross(dir), anchor),
      )
      const infos = thisTs.map((tThis) => {
        const p = this.at(tThis)
        return { tThis, tOther: L3.pointT(anchor, dir, p), p }
      })
      return thisIsPlanar
        ? infos
        : infos.filter((info) => L3.containsPoint(anchor, dir, info.p))
    }
  }

  roots() {
    console.log(this.tMin, this.tMax)
    arraySamples(this.tMin, this.tMax, 30).forEach((t) => {
      console.log(t + "," + this.tangentAt(t).z)
    })

    const result: Tuple3<number[]> = [[], [], []]
    for (let i = 0; i < this.points.length - 1; i++) {
      const findClosest = (startT: number, d: int) => {
        console.log("d", d, "startT", startT)
        // try {
        const root = newtonIterateWithDerivative2(
          (t) => {
            const [, dt, ddt] = this.ptDtDdt(t)
            return [dt.e(d), ddt.e(d)]
          },
          startT,
          8,
          this.tMin,
          this.tMax,
        )
        if (undefined !== root) {
          result[d].push(root)
        }
        console.log("d", d, "startT", startT, "root", root)
      }
      const a = this.points[i].p3()
      const b = this.points[i + 1].p3()
      const ab = a.to(b)
      for (let d = 0; d < 3; d++) {
        if (0 !== i && eq0(ab.e(d))) {
          const startT = lerp(
            this.knots[i],
            this.knots[i + this.degree + 2],
            0.5,
          )
          findClosest(startT, d)
        } else if (i < this.points.length - 2) {
          const bc = b.to(this.points[i + 2].p3())
          if (!eq0(bc.e(d)) && ab.e(d) * bc.e(d) < 0) {
            findClosest(
              this.closestTToPoint(b, this.guessTClosestToControlPoint(i + 1)),
              d,
            )
          }
        }
      }
    }
    console.log(result)
    return result
  }
  //getAABB() {
  //	return new AABB().addPoints(this.points.map(p => p.p3()))
  //}

  /**
   * Rough approximation of t param for points closest to control point.
   */
  guessTClosestToControlPoint(pointIndex: int) {
    return lerp(
      this.knots[pointIndex],
      this.knots[pointIndex + this.degree + 1],
      0.5,
    )
  }

  likeCurve(curve: Curve): boolean {
    return (
      this == curve ||
      (hasConstructor(curve, NURBS) &&
        this.degree === curve.degree &&
        this.points.every((p, i) => p.like(curve.points[i])) &&
        this.knots.every((k, i) => eq(k, curve.knots[i])))
    )
  }
  isColinearTo(curve: Curve): boolean {
    throw new Error("This doesn't even make sense.")
  }
}

NURBS.prototype.tIncrement = 1 / 128

function minAndMax(
  arr: ReadonlyArray<number>,
  start: int = 0,
  end: int = arr.length,
) {
  let min = Infinity,
    max = -Infinity
  for (let i = start; i < end; i++) {
    if (min > arr[i]) min = arr[i]
    if (max < arr[i]) max = arr[i]
  }
  return [min, max]
}
