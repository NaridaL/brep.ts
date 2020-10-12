import {
	AABB,
	arrayFromFunction,
	assert,
	assertf,
	assertInst,
	assertNever,
	assertNumbers,
	assertVectors,
	between,
	combinations,
	eq,
	eq0,
	fuzzyUniques,
	hasConstructor,
	int,
	lerp,
	M4,
	MINUS,
	newtonIterate1d,
	newtonIterate2dWithDerivatives,
	newtonIterateWithDerivative,
	NLA_PRECISION,
	solveCubicReal2,
	Tuple3,
	V,
	V3,
} from 'ts3dutils'
import { Mesh } from 'tsgl'

import {
	Curve,
	CylinderSurface,
	EllipseCurve,
	EllipsoidSurface,
	ISInfo,
	L3,
	P3,
	PlaneSurface,
	ProjectedCurveSurface,
	R2_R,
	Surface,
} from '../index'

import { abs, cos, PI, sin } from '../math'

/**
 * Bezier curve with degree 3.
 */
export class BezierCurve extends Curve {
	/**
	 * https://en.wikipedia.org/wiki/Cubic_function#/media/File:Graph_of_cubic_polynomial.svg
	 */
	static readonly EX2D = BezierCurve.graphXY(2, -3, -3, 2)
	static readonly EX3D = new BezierCurve(V3.O, V(-0.1, -1, 1), V(1.1, 1, 1), V3.X)
	static readonly QUARTER_CIRCLE = BezierCurve.approximateUnitArc(PI / 2)
	readonly p0: V3
	readonly p1: V3
	readonly p2: V3
	readonly p3: V3

	constructor(p0: V3, p1: V3, p2: V3, p3: V3, tMin: number = -0.1, tMax: number = 1.1) {
		super(tMin, tMax)
		assertVectors(p0, p1, p2, p3)
		assert(isFinite(tMin) && isFinite(tMax))
		//assert(!L3.throughPoints(p0, p3).containsPoint(p1) || !L3.throughPoints(p0, p3).containsPoint(p2))
		this.p0 = p0
		this.p1 = p1
		this.p2 = p2
		this.p3 = p3
	}

	get points(): V3[] {
		return [this.p0, this.p1, this.p2, this.p3]
	}

	/**
	 * Returns a curve with curve.at(x) == V(x, ax³ + bx² + cx + d, 0)
	 */
	static graphXY(a: number, b: number, c: number, d: number, tMin?: number, tMax?: number): BezierCurve {
		// d = p0y
		// c = -3 p0y + 3 p1y => p1y = c/3 + p0y
		// b = 3 p0y - 6 p1y + 3 p2y => p2y = b/3 - p0y + 2 p1y
		// a = -p0y + 3 p1y -3 p2y + p3y => p3y = a + p0y - 3 p1y + 3 p2y
		const p0y = d
		const p1y = c / 3 + p0y
		const p2y = b / 3 - p0y + 2 * p1y
		const p3y = a + p0y - 3 * p1y + 3 * p2y
		return new BezierCurve(V(0, p0y), V(1 / 3, p1y), V(2 / 3, p2y), V(1, p3y), tMin, tMax)
	}

	static quadratic(a: V3, b: V3, c: V3, tMin: number = 0, tMax: number = 1): BezierCurve | L3 {
		const line = L3.throughPoints(a, c)
		if (line.containsPoint(b)) {
			return line
		} else {
			// p1 = 1/3 a + 2/3 b
			// p2 = 1/3 c + 2/3 b
			return new BezierCurve(
				a,
				b
					.times(2)
					.plus(a)
					.div(3),
				b
					.times(2)
					.plus(c)
					.div(3),
				c,
				tMin,
				tMax,
			)
		}
	}

	/**
	 * Returns a bezier curve which approximates a CCW unit circle arc starting at V3.X of angle phi
	 * phi <= PI / 2 is recommended
	 *
	 * Formula from here: https://pomax.github.io/bezierinfo/#circles_cubic
	 */
	static approximateUnitArc(phi: number): BezierCurve {
		const f = (4 / 3) * Math.tan(phi / 4)
		return new BezierCurve(
			V3.X,
			new V3(1, f, 0),
			new V3(cos(phi) + f * sin(phi), sin(phi) - f * cos(phi), 0),
			V3.sphere(phi, 0),
			0,
			1,
		)
	}

	getConstructorParameters(): any[] {
		return [this.p0, this.p1, this.p2, this.p3]
	}

	at(t: number): V3 {
		// = s^3 p0 + 3 s^2 t p1 + 3 s t^2 p2 + t^3 p3
		assertNumbers(t)
		const p0 = this.p0,
			p1 = this.p1,
			p2 = this.p2,
			p3 = this.p3
		const s = 1 - t,
			c0 = s * s * s,
			c1 = 3 * s * s * t,
			c2 = 3 * s * t * t,
			c3 = t * t * t
		return new V3(
			p0.x * c0 + p1.x * c1 + p2.x * c2 + p3.x * c3,
			p0.y * c0 + p1.y * c1 + p2.y * c2 + p3.y * c3,
			p0.z * c0 + p1.z * c1 + p2.z * c2 + p3.z * c3,
		)
	}

	/**
	 * s := (1 - t)
	 * at(t) := s³ p0 + 3 s² t p1 + 3 s t² p2 + t³ p3
	 * tangent(t) := 3 s² (p1 - p0) + 6 s t (p2 - p1) + 3 t² (p3 - p2)
	 *            := 3 (1 - t)² (p1 - p0) + 6 (1 - t) t (p2 - p1) + 3 t² (p3 - p2)
	 *            := 3 (1 - 2 t + t²) (p1 - p0) + 6 (t - t²) (p2 - p1) + 3 t² (p3 - p2)
	 *            := (3 (p3 - p2) - 6 (p2 - p1) + 3 (p1 - p0)) t²*
	 *                + (-6 (p1 - p0) + (p2 - p1)) t
	 *                + 3 (p1 - p0)
	 */
	tangentAt(t: number): V3 {
		assertNumbers(t)
		const p0 = this.p0,
			p1 = this.p1,
			p2 = this.p2,
			p3 = this.p3
		const s = 1 - t,
			c01 = 3 * s * s,
			c12 = 6 * s * t,
			c23 = 3 * t * t
		return new V3(
			(p1.x - p0.x) * c01 + (p2.x - p1.x) * c12 + (p3.x - p2.x) * c23,
			(p1.y - p0.y) * c01 + (p2.y - p1.y) * c12 + (p3.y - p2.y) * c23,
			(p1.z - p0.z) * c01 + (p2.z - p1.z) * c12 + (p3.z - p2.z) * c23,
		)
	}

	ddt(t: number): V3 {
		assertNumbers(t)
		const p0 = this.p0,
			p1 = this.p1,
			p2 = this.p2,
			p3 = this.p3
		const c012 = 6 * (1 - t),
			c123 = 6 * t
		return new V3(
			(p2.x - 2 * p1.x + p0.x) * c012 + (p3.x - 2 * p2.x + p1.x) * c123,
			(p2.y - 2 * p1.y + p0.y) * c012 + (p3.y - 2 * p2.y + p1.y) * c123,
			(p2.z - 2 * p1.z + p0.z) * c012 + (p3.z - 2 * p2.z + p1.z) * c123,
		)
	}

	normalP(t: number): V3 {
		const tangent = this.tangentAt(t)
		const rot = tangent.cross(this.ddt(t))
		return rot.cross(tangent)
	}

	isTsWithPlane(planeWC: P3) {
		assertInst(P3, planeWC)
		/*
		 We are solving for t:
		 n := plane.normal1
		 this.at(t) DOT n == plane.w // according to plane definition
		 (a t³ + b t² + c t + d) DOT n == plane.w // bezier curve as cubic equation
		 (a DOT n) t³ + (b DOT n) t³ + (c DOT n) t + d DOT n - plane.w == 0 // multiply out DOT n, minus plane.w
		 */

		const { p0, p1, p2, p3 } = this
		const n = planeWC.normal1
		const a = p1
			.minus(p2)
			.times(3)
			.minus(p0)
			.plus(p3)
		const b = p0
			.plus(p2)
			.times(3)
			.minus(p1.times(6))
		const c = p1.minus(p0).times(3)
		const d = p0

		return solveCubicReal2(a.dot(n), b.dot(n), c.dot(n), d.dot(n) - planeWC.w).filter(t =>
			between(t, this.tMin, this.tMax),
		)
	}

	isTsWithSurface(surfaceWC: Surface): number[] {
		if (surfaceWC instanceof CylinderSurface) {
			const projPlane = new P3(surfaceWC.dir.unit(), 0)
			const projThis = this.project(projPlane)
			const projEllipse = surfaceWC.baseCurve.project(projPlane)
			return projEllipse.isInfosWithBezier2D(projThis).map(info => info.tOther)
		}
		return super.isTsWithSurface(surfaceWC)
	}

	likeCurve(curve: Curve): boolean {
		return (
			this == curve ||
			(hasConstructor(curve, BezierCurve) &&
				this.p0.like(curve.p0) &&
				this.p1.like(curve.p1) &&
				this.p2.like(curve.p2) &&
				this.p3.like(curve.p3))
		)
	}

	equals(obj: any): boolean {
		return (
			this == obj ||
			(hasConstructor(obj, BezierCurve) &&
				this.p0.equals(obj.p0) &&
				this.p1.equals(obj.p1) &&
				this.p2.equals(obj.p2) &&
				this.p3.equals(obj.p3))
		)
	}

	hashCode(): int {
		let hashCode = 0
		hashCode = hashCode * 31 + this.p0.hashCode()
		hashCode = hashCode * 31 + this.p1.hashCode()
		hashCode = hashCode * 31 + this.p2.hashCode()
		hashCode = hashCode * 31 + this.p3.hashCode()
		return hashCode | 0
	}

	/**
	 * Checks if this curve is colinear to the passed curve, i.e.
	 * for every t:number there exists a s:number with this.at(t) = curve.at(s)
	 */
	isColinearTo(curve: BezierCurve): boolean {
		if (this === curve || this.likeCurve(curve)) return true
		if (!(curve instanceof BezierCurve)) return false
		// first, find out where/if curve.p0 and curve.p3 are on this
		// then split this at curve.p0 --> curve.p3 to compare points p1 and p2
		let curveP0T, curveP3T
		// assign in if condition to exploit short-circuit
		if (isNaN((curveP0T = this.pointT(curve.p0))) || isNaN((curveP3T = this.pointT(curve.p3)))) {
			return false
		}
		let thisSplit
		if (eq(1, curveP0T)) {
			// this.split(curveP0T).right is degenerate in this case, so we need to handle it separately

			// this.split(curveP3T): 0 --> curveP3T --> 1
			// .right: curveP3T --> 1
			// .reversed(): 1 --> curveP3T
			thisSplit = this.split(curveP3T)[1].reversed()
		} else {
			// curveP3T describes the point on this
			// adjust it so it describes the same point on this.split(curveP0T).right
			// this:                       0           p0t        p3t      1
			//                             |            |          |       |
			// this.split(curveP0T).right:              0        p3tad     1
			const curveP3Tadjusted = (curveP3T - curveP0T) / (1 - curveP0T)
			thisSplit = this.split(curveP0T)[1].split(curveP3Tadjusted)[0]
		}

		return curve.likeCurve(thisSplit)
	}

	selectPart(t0: number, t1: number) {
		const t1Adjusted = (t1 - t0) / (1 - t0)
		return this.split(t0)[1].split(t1Adjusted)[0]
	}

	reversed(): BezierCurve {
		return new BezierCurve(this.p3, this.p2, this.p1, this.p0, 1 - this.tMax, 1 - this.tMin)
	}

	getCoefficients() {
		const { p0, p1, p2, p3 } = this
		// calculate cubic equation coefficients
		// a t³ + b t² + c t + d = 0
		// multiplying out the cubic Bézier curve equation gives:
		// a = -p0 + 3 p1 - 3 p2 + p3
		// b = 3 p0 - 6 p1 + 3 p2
		// c = -3 p0 + 3 p1
		// d = p0 - p
		const a = p1
			.minus(p2)
			.times(3)
			.minus(p0)
			.plus(p3)
		const b = p0
			.plus(p2)
			.times(3)
			.minus(p1.times(6))
		const c = p1.minus(p0).times(3)
		const d = p0
		return [a, b, c, d]
	}

	tangentCoefficients() {
		const { p0, p1, p2, p3 } = this
		const p01 = p1.minus(p0),
			p12 = p2.minus(p1),
			p23 = p3.minus(p2)
		const a = p01
			.plus(p23)
			.times(3)
			.minus(p12.times(6))
		const b = p12.minus(p01).times(6)
		const c = p01.times(3)
		return [V3.O, a, b, c]
	}

	pointT2(p: V3, tMin = this.tMin, tMax = this.tMax): number {
		const t = this.closestTToPoint(p, undefined, tMin, tMax)
		assert(this.at(t).like(p))
		return t
	}

	pointT(p: V3) {
		const { p0, p1, p2, p3 } = this
		// calculate cubic equation coefficients
		// a t³ + b t² + c t + d = 0
		// multiplying out the cubic Bézier curve equation gives:
		// a = -p0 + 3 p1 - 3 p2 + p3
		// b = 3 p0 - 6 p1 + 3 p2
		// c = -3 p0 + 3 p1
		// d = p0 - p
		const a = p1
			.minus(p2)
			.times(3)
			.minus(p0)
			.plus(p3)
		const b = p0
			.plus(p2)
			.times(3)
			.minus(p1.times(6))
		const c = p1.minus(p0).times(3)
		const d = p0.minus(p)

		// a t³ + b t² + c t + d = 0 is 3 cubic equations, some of which can be degenerate
		const maxDim =
			NLA_PRECISION < a.maxAbsElement()
				? a.maxAbsDim()
				: NLA_PRECISION < b.maxAbsElement()
				? b.maxAbsDim()
				: NLA_PRECISION < c.maxAbsElement()
				? c.maxAbsDim()
				: assertNever()

		const results = solveCubicReal2(a.e(maxDim), b.e(maxDim), c.e(maxDim), d.e(maxDim)).filter(t =>
			this.at(t).like(p),
		)
		if (0 == results.length) return NaN
		if (1 == results.length) return results[0]
		throw new Error('multiple intersection ' + this.toString() + p.sce)
	}

	pointT3(p: V3) {
		const { p0, p1, p2, p3 } = this
		// calculate cubic equation coefficients
		// a t³ + b t² + c t + d = 0
		// multiplying out the cubic Bézier curve equation gives:
		// a = -p0 + 3 p1 - 3 p2 + p3
		// b = 3 p0 - 6 p1 + 3 p2
		// c = -3 p0 + 3 p1
		// d = p0 - p
		const a = p1
			.minus(p2)
			.times(3)
			.minus(p0)
			.plus(p3)
			.els()
		const b = p0
			.plus(p2)
			.times(3)
			.minus(p1.times(6))
			.els()
		const c = p1
			.minus(p0)
			.times(3)
			.els()
		const d = p0.minus(p).els()
		let results = undefined

		// assume passed point is on curve and that curve does not self-intersect,
		// i.e. there is exactly one correct result for t
		// try to find a single result in the x-dimension, if multiple are found,
		// filter them by checking the other dimensions
		for (let dim = 0; dim < 3; dim++) {
			if (eq0(a[dim]) && eq0(b[dim]) && eq0(c[dim])) {
				// for case x:
				// ax == bx == cx == 0 => x(t) = dx
				// x value is constant
				// if x == 0 for all t, this does not limit the result, otherwise, there is no result, i.e
				// the passed point is not on the curve
				if (!eq0(d[dim])) return NaN
			} else {
				const newResults = solveCubicReal2(a[dim], b[dim], c[dim], d[dim])
				if (0 == newResults.length) return NaN
				if (1 == newResults.length) return newResults[0]
				if (results) {
					results = results.filter(t => newResults.some(t2 => eq(t, t2)))
					if (0 == results.length) return NaN
					if (1 == results.length) return results[0]
				} else {
					results = newResults
				}
			}
		}
		throw new Error('multiple intersection ' + results + this.toString() + p.sce)
	}

	transform(m4: M4) {
		// perspective projection turn bezier curve into rational spline
		assert(m4.isNoProj(), m4.str)
		return new BezierCurve(
			m4.transformPoint(this.p0),
			m4.transformPoint(this.p1),
			m4.transformPoint(this.p2),
			m4.transformPoint(this.p3),
			this.tMin,
			this.tMax,
		) as this
	}

	isClosed(): boolean {
		return this.p0.like(this.p3)
	}

	isQuadratic(): boolean {
		return this.p0.lerp(this.p1, 1.5).like(this.p3.lerp(this.p2, 1.5))
	}

	debugInfo() {
		return {
			lines: [0, 1, 1, 2, 2, 3].map(i => this.points[i]),
			points: this.points,
		}
	}

	split(t: number): [BezierCurve, BezierCurve] {
		// do de Casteljau's algorithm at t, the resulting points are the points needed to create 2 new curves
		const s = 1 - t
		const { p0, p1, p2, p3 } = this
		/*
		p3 // n3
		b01 = s p0 + t p1
		b11 = s p1 + t p2
		b21 = s p2 + t p3 // n2
		b02 = s b01 + t b11
		b12 = s b11 + t b21 // n1
		b03 = s b02 + t b12 // n0

		c01 =
		*/
		const b01 = p0.times(s).plus(p1.times(t)),
			b11 = p1.times(s).plus(p2.times(t)),
			b21 = p2.times(s).plus(p3.times(t))
		const b02 = b01.times(s).plus(b11.times(t)),
			b12 = b11.times(s).plus(b21.times(t))
		const b03 = b02.times(s).plus(b12.times(t))
		return [new BezierCurve(p0, b01, b02, b03), new BezierCurve(b03, b12, b21, p3)]
	}

	containsPoint(p: V3): boolean {
		return isFinite(this.pointT(p))
	}

	roots(): Tuple3<number[]> {
		/**
		 *            := (3 (p3 - p2) - 6 (p2 - p1) + 3 (p1 - p0)) t²*
		 *                + (-6 (p1 - p0) + 6 (p2 - p1)) t
		 *                + 3 (p1 - p0)
		 *                */
		const { p0, p1, p2, p3 } = this
		const p01 = p1.minus(p0),
			p12 = p2.minus(p1),
			p23 = p3.minus(p2)
		const a = p01
			.plus(p23)
			.times(3)
			.minus(p12.times(6))
		const b = p12.minus(p01).times(6)
		const c = p01.times(3)

		return arrayFromFunction(3, dim => solveCubicReal2(0, a.e(dim), b.e(dim), c.e(dim)))
	}

	isInfosWithLine(
		anchorWC: V3,
		dirWC: V3,
		tMin?: number,
		tMax?: number,
		lineMin = -100000,
		lineMax = 100000,
	): ISInfo[] {
		// const dirLength = dirWC.length()
		// // TODO: no:
		// let result = Curve.ispsRecursive(this, this.tMin, this.tMax, new L3(anchorWC, dirWC.unit()), lineMin, lineMax)
		// result = fuzzyUniquesF(result, info => info.tOther)
		// result.forEach(info => (info.tOther /= dirLength))
		// return result
		// looking for this.at(t) == line.at(s)
		// this.at(t).x == anchorWC.x + dirWC.x * s
		// (this.at(t).x - anchorWC.x) / dirWC.x == s (analogue for y and z) (1x, 1y, 1z)
		// (1x) - (1y):
		// (this.at(t).x - anchorWC.x) / dirWC.x - (this.at(t).y - anchorWC.y) / dirWC.y == 0
		// (this.at(t).x - anchorWC.x) * dirWC.y - (this.at(t).y - anchorWC.y) * dirWC.x == 0 (2)

		// cubic equation params (see #pointT):
		const { p0, p1, p2, p3 } = this
		const a = p1
			.minus(p2)
			.times(3)
			.minus(p0)
			.plus(p3)

		const v1 = V3.UNITS[a.minAbsDim()]
		const testPlane = P3.forAnchorAndPlaneVectors(anchorWC, dirWC, v1.isParallelTo(dirWC) ? a : v1)

		return this.isTsWithPlane(testPlane)
			.map(tThis => {
				const p = this.at(tThis)
				return { tThis, tOther: L3.pointT(anchorWC, dirWC, p), p }
			})
			.filter(info => L3.containsPoint(anchorWC, dirWC, info.p))
	}

	closestPointToLine(line: L3, tMin: number, tMax: number) {
		// (this(t)-line(s)) * line.dir == 0 (1)
		// (this(t)-line(s)) * this.tangentAt(t) == 0 (2)
		// this(t) * line.dir - line(s) * line.dir == 0
		// this(t) * line.dir - line.anchor * line.dir - s line.dir * line.dir == 0
		// this(t) * line.dir - line.anchor * line.dir == s (3)
		// insert (3) in (2)
		// (this(t)-line(this(t) * line.dir - line.anchor * line.dir)) * this.tangentAt(t) == 0 (4)
		// (4) is a 5th degree polynomial, solve numerically

		tMin = isFinite(tMin) ? tMin : this.tMin
		tMax = isFinite(tMax) ? tMax : this.tMax

		const anchorDotDir1 = line.anchor.dot(line.dir1)
		const f = (t: number) => {
			const atT = this.at(t)
			return atT.minus(line.at(atT.dot(line.dir1) - anchorDotDir1)).dot(this.tangentAt(t))
		}

		const STEPS = 32
		const startT = arrayFromFunction(STEPS, i => tMin + ((tMax - tMin) * i) / STEPS).withMax(t => -f(t))

		return newtonIterate1d(f, startT, 8)
	}

	/**
	 *
	 * @param bezier
	 * @param tMin
	 * @param tMax
	 * @param sMin
	 * @param {number=} sMax
	 * @returns
	 */
	isInfosWithBezier3(bezier: BezierCurve, tMin?: number, tMax?: number, sMin?: number, sMax?: number) {
		const handleStartTS = (startT: number, startS: number) => {
			if (!result.some(info => eq(info.tThis, startT) && eq(info.tOther, startS))) {
				const f1: R2_R = (t, s) => this.tangentAt(t).dot(this.at(t).minus(bezier.at(s)))
				const f2: R2_R = (t, s) => bezier.tangentAt(s).dot(this.at(t).minus(bezier.at(s)))
				// f = (b1, b2, t1, t2) = b1.tangentAt(t1).dot(b1.at(t1).minus(b2.at(t2)))
				const fdt1 = (b1: BezierCurve, b2: BezierCurve, t1: number, t2: number) =>
					b1.ddt(t1).dot(b1.at(t1).minus(b2.at(t2))) + b1.tangentAt(t1).squared()
				const fdt2 = (b1: BezierCurve, b2: BezierCurve, t1: number, t2: number) =>
					-b1.tangentAt(t1).dot(b2.tangentAt(t2))
				const ni = newtonIterate2dWithDerivatives(
					f1,
					f2,
					startT,
					startS,
					16,
					fdt1.bind(undefined, this, bezier),
					fdt2.bind(undefined, this, bezier),
					(t, s) => -fdt2(bezier, this, s, t),
					(t, s) => -fdt1(bezier, this, s, t),
				)!
				result.push({ tThis: ni.x, tOther: ni.y, p: this.at(ni.x) })
			}
		}

		tMin = undefined !== tMin ? tMin : this.tMin
		tMax = undefined !== tMax ? tMax : this.tMax
		sMin = undefined !== sMin ? sMin : bezier.tMin
		sMax = undefined !== sMax ? sMax : bezier.tMax

		// stack of indices:
		const indices = [tMin, tMax, sMin, sMax]
		const result: ISInfo[] = []
		while (indices.length) {
			const i = indices.length - 4
			const tMin = indices[i],
				tMax = indices[i + 1],
				sMin = indices[i + 2],
				sMax = indices[i + 3]
			indices.length -= 4
			const thisAABB = this.getAABB(tMin, tMax)
			const otherAABB = bezier.getAABB(sMin, sMax)
			// console.log(tMin, tMax, sMin, sMax, thisAABB.sce, otherAABB.sce)
			if (thisAABB && otherAABB && thisAABB.intersectsAABB2d(otherAABB)) {
				const tMid = (tMin + tMax) / 2
				const sMid = (sMin + sMax) / 2
				const EPS = 0.00001
				if (tMax - tMin < EPS || sMax - sMin < EPS) {
					console.log(tMin, tMax, sMin, sMax)
					console.log(thisAABB.sce)
					console.log(otherAABB.sce)
					console.log(tMid, sMid)
					handleStartTS(tMid, sMid)
				} else {
					indices.push(
						tMin,
						tMid,
						sMin,
						sMid,
						tMin,
						tMid,
						sMid,
						sMax,
						tMid,
						tMax,
						sMin,
						sMid,
						tMid,
						tMax,
						sMid,
						sMax,
					)
				}
			}
		}

		return result
	}

	isInfosWithBezier(bezier: BezierCurve, tMin?: number, tMax?: number, sMin?: number, sMax?: number): ISInfo[] {
		tMin = undefined !== tMin ? tMin : this.tMin
		tMax = undefined !== tMax ? tMax : this.tMax
		sMin = undefined !== sMin ? sMin : bezier.tMin
		sMax = undefined !== sMax ? sMax : bezier.tMax

		assertf(() => tMin! < tMax!)
		assertf(() => sMin! < sMax!)
		const result: ISInfo[] = []

		const likeCurves = this.likeCurve(bezier),
			colinearCurves = this.isColinearTo(bezier)
		if (likeCurves || colinearCurves) {
			if (!likeCurves) {
				// only colinear
				// recalculate sMin and sMax so they are valid on this, from then on we can ignore bezier
				sMin = this.pointT(bezier.at(sMin))
				sMax = this.pointT(bezier.at(sMax))
			}
			tMin = Math.min(tMin, sMin)
			tMax = Math.max(tMax, sMax)
			const splits = fuzzyUniques(
				this.roots()
					.concatenated()
					.filter(isFinite)
					.concat([tMin, tMax]),
			).sort(MINUS)
			//const aabbs = arrayFromFunction(splits.length - 1, i => this.getAABB(splits[i], splits[i + 1]))
			Array.from(combinations(splits.length - 1)).forEach(({ i, j }) => {
				// adjacent curves can't intersect
				if (Math.abs(i - j) > 2) {
					// console.log(splits[i], splits[i + 1], splits[j], splits[j + 1], aabbs[i], aabbs[j])
					//findRecursive(splits[i], splits[i + 1], splits[j], splits[j + 1], aabbs[i], aabbs[j])
					result.push(
						...Curve.ispsRecursive(this, splits[i], splits[i + 1], bezier, splits[j], splits[j + 1]),
					)
				}
			})
		} else {
			return Curve.ispsRecursive(this, tMin, tMax, bezier, sMin, sMax)
		}

		return result
	}

	selfIntersectionsInfo() {
		return this.isInfosWithBezier(this)
	}

	isInfosWithCurve(curve: Curve): ISInfo[] {
		if (curve instanceof L3) {
			return this.isInfosWithLine(curve.anchor, curve.dir1, curve.tMin, curve.tMax)
		}
		if (curve instanceof BezierCurve) {
			return this.isInfosWithBezier(curve)
		}
		return curve.isInfosWithCurve(this).map(({ tThis, tOther, p }) => ({ tThis: tOther, tOther: tThis, p }))
	}

	/**
	 * Approximate this bezier curve with a number of circular segments. This curve is recursively split in half until
	 * segments are close enough (relative error < REL_ERR in two test points) to an arc which goes through the start,
	 * end and mid points of the segment.
	 * @returns each EllipseCurve is circular and their tMin and tMax respectively define their start and end points.
	 * @param t0 Start parameter of segment which should be approximated.
	 * @param t1 End parameter of segment which should be approximated.
	 * @param REL_ERROR max allowable relative error.
	 * @param result Resulting circle arcs are stored in this array. Mainly used by the recursion.
	 */
	circleApprox(
		t0: number = this.tMin,
		t1: number = this.tMax,
		REL_ERROR = 1 / 1024,
		result: EllipseCurve[] = [],
	): EllipseCurve[] {
		const a = this.at(t0),
			b = this.at(t1),
			tMid = (t0 + t1) / 2,
			pMid = this.at(tMid),
			abLine = L3.throughPoints(a, b)
		if (!abLine.containsPoint(pMid) && between(abLine.pointT(pMid), 0, abLine.pointT(b))) {
			const arc = EllipseCurve.circleThroughPoints(a, pMid, b),
				arcRadius = arc.f1.length(),
				pTest1 = this.at(lerp(t0, t1, 0.25)),
				pTest2 = this.at(lerp(t0, t1, 0.75))
			if (
				abs(arc.center.distanceTo(pTest1) / arcRadius - 1) <= REL_ERROR &&
				abs(arc.center.distanceTo(pTest2) / arcRadius - 1) <= REL_ERROR
			) {
				result.push(arc)
				return result
			}
		}
		this.circleApprox(t0, tMid, REL_ERROR, result)
		this.circleApprox(tMid, t1, REL_ERROR, result)
		return result
	}
}

BezierCurve.prototype.hlol = Curve.hlol++
BezierCurve.prototype.tIncrement = 1 / 80
