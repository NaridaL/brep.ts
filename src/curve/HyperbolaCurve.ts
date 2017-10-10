import {arrayFromFunction, assertNumbers, eq, eq0, hasConstructor, le, snap0, V3} from 'ts3dutils'

import {Curve, XiEtaCurve, intersectionUnitHyperbolaLine} from '../index'

const {PI, cos, sin, min, max, tan, sign, ceil, floor, abs, sqrt, pow, atan2, round} = Math

/**
 * x² - y² = 1
 *
 */
export class HyperbolaCurve extends XiEtaCurve {
	static XY = new HyperbolaCurve(V3.O, V3.X, V3.Y)

	constructor(center: V3, f1: V3, f2: V3, tMin: number = -7, tMax: number = 7) {
		super(center, f1, f2, tMin, tMax)
	}

	static XYLCValid(pLC: V3): boolean {
		return pLC.x > 0 && eq(1, pLC.x * pLC.x - pLC.y * pLC.y)
	}

	static XYLCPointT(pLC: V3): number {
		return Math.asinh(pLC.y)
	}

	/**
	 * http://www.wolframalpha.com/input/?i=x%C2%B2-y%C2%B2%3D1,ax%2Bby%3Dc
	 * Minor empiric test shows asinh(eta) consistently gets more accurate results than atanh(eta/xi)
	 */
	static magic(a: number, b: number, c: number): number[] {
		if (eq0(b)) {
			const sqrtVal = snap0(c ** 2 / a ** 2 - 1)
			if (sqrtVal < 0 || c * a < 0) {
				return []
			} else if (sqrtVal == 0) {
				return [0]
			}
			const eta1 = Math.sqrt(sqrtVal)
			return [-Math.asinh(eta1), Math.asinh(eta1)]
		} else if (eq(abs(a), abs(b))) {
			if (le(c * a, 0)) {
				return []
			}
			const eta = sign(a * b) * (c ** 2 - a ** 2) / 2 / a / c
			return [Math.asinh(eta)]
		} else {
			const sqrtVal = snap0(b ** 2 * (-(a ** 2) + b ** 2 + c ** 2))
			if (sqrtVal < 0) {
				return []
			}
			const xi1 = (a * c - Math.sqrt(sqrtVal)) / (a ** 2 - b ** 2)
			const xi2 = (a * c + Math.sqrt(sqrtVal)) / (a ** 2 - b ** 2)
			const eta1 = (b ** 2 * c - a * Math.sqrt(sqrtVal)) / (b * (b ** 2 - a ** 2))
			const eta2 = (b ** 2 * c + a * Math.sqrt(sqrtVal)) / (b * (b ** 2 - a ** 2))
			const foo: number = 20
			const bar = foo > 0 && foo
			return [xi1 > 0 && Math.asinh(eta1), xi2 > 0 && Math.asinh(eta2)].filter((x: any) => x !== false)
		}

	}

	at(t: number): V3 {
		assertNumbers(t)
		// = center + f1 cosh t + f2 sinh t
		return this.center.plus(this.f1.times(Math.cosh(t))).plus(this.f2.times(Math.sinh(t)))
	}

	tangentAt(t: number): V3 {
		assertNumbers(t)
		// = f1 sinh t + f2 cosh t
		return this.f1.times(Math.sinh(t)).plus(this.f2.times(Math.cosh(t)))
	}

	tangentAt2(xi: number, eta: number): V3 {
		assertNumbers(xi, eta)
		// = f1 eta + f2 xi
		return this.f1.times(eta).plus(this.f2.times(xi))
	}

	ddt(t: number): V3 {
		assertNumbers(t)
		return this.f1.times(Math.cosh(t)).plus(this.f2.times(Math.sinh(t)))
	}

	isColinearTo(curve: Curve): boolean {
		if (!hasConstructor(curve, HyperbolaCurve)) return false
		if (!curve.center || !this.center.like(curve.center)) {
			return false
		}
		if (this === curve) {
			return true
		}
		const {f1: f1, f2: f2} = this.rightAngled(), {f1: c1, f2: c2} = curve.rightAngled()
		return eq(f1.squared(), Math.abs(f1.dot(c1)))
			&& eq(f2.squared(), Math.abs(f2.dot(c2)))
	}

	reversed() {
		return new HyperbolaCurve(this.center, this.f1, this.f2.negated(), -this.tMax, -this.tMin)
	}

	rightAngled(): HyperbolaCurve {
		const f1 = this.f1, f2 = this.f2, a = f1.dot(f2), b = f2.squared() + f1.squared()
		if (eq0(a)) {
			return this
		}
		const g1 = 2 * a, g2 = b + Math.sqrt(b * b - 4 * a * a)
		const {x1: xi, y1: eta} = intersectionUnitHyperbolaLine(g1, g2, 0)
		return new HyperbolaCurve(this.center, f1.times(xi).plus(f2.times(eta)), f1.times(eta).plus(f2.times(xi)))
	}

	eccentricity(): number {
		const mainAxes = this.rightAngled()
		const f1length = mainAxes.f1.length(), f2length = mainAxes.f1.length()
		const [a, b] = f1length > f2length ? [f1length, f2length] : [f2length, f1length]
		return Math.sqrt(1 + b * b / a / a)
	}

	roots(): [number[], number[], number[]] {
		// tangent(t) = f1 sinh t + f2 cosh t = 0
		// tangentAt2(xi, eta) = f1 eta + f2 xi = V3.O
		// xi² - eta² = 1 (by def for hyperbola)

		return arrayFromFunction(3, dim => {
			const a = this.f2.e(dim), b = this.f1.e(dim)
			return HyperbolaCurve.magic(a, b, 0)
		})
	}
}

HyperbolaCurve.prototype.tIncrement = PI / 16