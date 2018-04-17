import {
	arrayFromFunction,
	assert,
	assertf,
	assertNumbers,
	assertVectors,
	between,
	checkDerivate,
	eq,
	eq0,
	fuzzyBetween,
	hasConstructor,
	int,
	le,
	lerp,
	lt,
	M4,
	newtonIterate1d,
	newtonIterateSmart,
	pqFormula,
	TAU,
	V3,
} from 'ts3dutils'

import { Curve, intersectionUnitCircleLine, intersectionUnitCircleLine2, ISInfo, L3, P3, XiEtaCurve } from '../index'

import { atan2, max, min, PI } from '../math'

export class EllipseCurve extends XiEtaCurve {
	static readonly UNIT = new EllipseCurve(V3.O, V3.X, V3.Y)

	constructor(center: V3, f1: V3, f2: V3, tMin: number = 0, tMax: number = PI) {
		super(center, f1, f2, tMin, tMax)
		assert(-PI <= this.tMin && this.tMin < PI)
		assert(-PI < this.tMax && this.tMax <= PI)
	}

	static XYLCValid(pLC: V3): boolean {
		const { x, y } = pLC
		return eq0(x ** 2 + y ** 2 - 1)
	}

	static XYLCPointT(pLC: V3, tMin: number, tMax: number): number {
		assertNumbers(tMin, tMax)
		const t = atan2(pLC.y, pLC.x)
		const lowSplitter = lerp(tMin, tMax - TAU, 0.5)
		if (t < lowSplitter) {
			return t + TAU
		}
		const highSplitter = lerp(tMax, tMin + TAU, 0.5)
		if (t > highSplitter) {
			return t - TAU
		}
		return t
	}

	static intersectionUnitLine(a: number, b: number, c: number, tMin: number, tMax: number): number[] {
		const isLC = intersectionUnitCircleLine2(a, b, c)
		const result = []
		for (const [xi, eta] of isLC) {
			const t = EllipseCurve.XYLCPointT(new V3(xi, eta, 0), tMin, tMax)
			fuzzyBetween(t, tMin, tMax) && result.push(t)
		}
		return result
	}

	static unitIsInfosWithLine(anchorLC: V3, dirLC: V3, anchorWC: V3, dirWC: V3, tMin: number, tMax: number): ISInfo[] {
		// ell: x² + y² = 1 = p²
		// line(t) = anchor + t dir
		// anchor² - 1 + 2 t dir anchor + t² dir² = 0
		const pqDiv = dirLC.squared()
		const lineTs = pqFormula(2 * dirLC.dot(anchorLC) / pqDiv, (anchorLC.squared() - 1) / pqDiv)
		return lineTs.filter(tOther => le(0, anchorLC.y + tOther * dirLC.y)).map(tOther => ({
			tThis: EllipseCurve.XYLCPointT(dirLC.times(tOther).plus(anchorLC), tMin, tMax),
			tOther: tOther,
			p: L3.at(anchorWC, dirWC, tOther),
		}))
	}

	/**
	 * Returns a new EllipseCurve representing a circle parallel to the XY-plane.`
	 */
	static semicircle(radius: number, center: V3 = V3.O, tMin?: number, tMax?: number): EllipseCurve {
		return new EllipseCurve(center, new V3(radius, 0, 0), new V3(0, radius, 0), tMin, tMax)
	}

	static circleForCenter2P(center: V3, a: V3, b: V3, radius: number, tMin?: number, tMax?: number) {
		const f1 = center.to(a)
		const normal = f1.cross(center.to(b))
		const f2 = normal.cross(f1).toLength(f1.length())
		return new EllipseCurve(
			center,
			f1,
			f2,
			undefined !== tMin ? tMin : 0,
			undefined !== tMax ? tMax : f1.angleTo(center.to(b)),
		)
	}

	split(tMin = this.tMin, tMax = this.tMax): EllipseCurve[] {
		const result: EllipseCurve[] = []
		tMin < 0 &&
			result.push(
				new EllipseCurve(this.center, this.f1.negated(), this.f2.negated(), tMin + PI, min(0, tMax) + PI),
			)
		tMax > 0 && result.push(new EllipseCurve(this.center, this.f1, this.f2, max(0, tMin), tMax))
		return result
	}

	static forAB(a: number, b: number, center: V3 = V3.O): EllipseCurve {
		return super.forAB(a, b, center) as EllipseCurve
	}

	/**
	 * Create a circle curve which has a, b and c on it. a, b, c can't be on a straight line.
	 * tMin defaults to 0, tMax defaults to the value for c
	 */
	static circleThroughPoints(a: V3, b: V3, c: V3, tMin = 0, tMax?: number) {
		assertf(() => !L3.throughPoints(a, c).containsPoint(b))
		const normal = a.to(b).cross(b.to(c))
		const center = new L3(a.lerp(b, 0.5), normal.cross(a.to(b)).unit()).isInfoWithLine(
			new L3(b.lerp(c, 0.5), normal.cross(b.to(c)).unit()),
		)!
		const f1 = center.to(a).negated()
		return new EllipseCurve(
			center,
			f1,
			normal.unit().cross(f1),
			-PI,
			undefined === tMax ? f1.angleRelativeNormal(center.to(c), normal.unit()) : tMax,
		)
	}

	getAreaInDir(right: V3, up: V3, tStart: number, tEnd: number): { area: number; centroid: V3 } {
		//assertf(() => tStart < tEnd)
		assertf(() => right.isPerpendicularTo(this.normal))
		assertf(() => up.isPerpendicularTo(this.normal))
		//assertf(() => EllipseCurve.isValidT(tStart), tStart)
		//assertf(() => EllipseCurve.isValidT(tEnd), tEnd)

		const upLC = this.matrixInverse.transformVector(up)
		const rightLC = upLC.cross(V3.Z)
		const normTStart = tStart - rightLC.angleXY()
		const normTEnd = tEnd - rightLC.angleXY()
		const transformedOriginY = this.matrixInverse.getTranslation().dot(upLC.unit())
		// integral of sqrt(1 - x²) from 0 to cos(t)
		// Basically, we want
		// INTEGRAL[cos(t); PI/2] sqrt(1 - x²) dx
		// INTEGRAL[PI/2: cos(t)] -sqrt(1 - x²) dx
		// = INTEGRAL[cos(0); cos(t)] -sqrt(1 - x²) dx
		// = INTEGRAL[0; t] -sqrt(1 - cos²(t)) * -sin(t) dt
		// = INTEGRAL[0; t] -sin(t) * -sin(t) dt
		// = INTEGRAL[0; t] sin²(t) dt (partial integration / wolfram alpha)
		// = (1/2 * (t - sin(t) * cos(t)))[0; t] (this form has the distinct advantage of being defined everywhere)
		function fArea(t: number) {
			return (t - Math.sin(t) * Math.cos(t)) / 2
		}

		// for the centroid, we want
		// cx = 1 / area * INTEGRAL[cos(t); PI/2] x * f(x) dx
		// cx = 1 / area * INTEGRAL[cos(t); PI/2] x * sqrt(1 - x²) dx
		// cx = 1 / area * INTEGRAL[cos(0); cos(t)] x * -sqrt(1 - x²) dx
		// ...
		// cx = 1 / area * INTEGRAL[0; t] cos(t) * sin²(t) dt // WA
		// cx = 1 / area * (sin^3(t) / 3)[0; t]
		function cxTimesArea(t: number) {
			return Math.pow(Math.sin(t), 3) / 3
		}

		// cy = 1 / area * INTEGRAL[cos(t); PI/2] f²(x) / 2 dx
		// cy = 1 / area * INTEGRAL[cos(0); cos(t)] -(1 - x²) / 2 dx
		// cy = 1 / area * INTEGRAL[0; t] (cos²(t) - 1) * -sin(t) / 2 dt
		// cy = 1 / area * (cos (3 * t) - 9 * cos(t)) / 24 )[0; t]
		function cyTimesArea(t: number) {
			return (Math.cos(3 * t) - 9 * Math.cos(t)) / 24
		}

		const restArea = -transformedOriginY * (-Math.cos(normTEnd) + Math.cos(normTStart))
		const area = fArea(normTEnd) - fArea(normTStart) + restArea
		const cxt =
			(cxTimesArea(normTEnd) -
				cxTimesArea(normTStart) +
				-transformedOriginY * (-Math.cos(normTEnd) - Math.cos(normTStart)) / 2 * restArea) /
			area
		const cyt = (cyTimesArea(normTEnd) - cyTimesArea(normTStart) - -transformedOriginY / 2 * restArea) / area
		const factor = this.matrix.xyAreaFactor() // * upLC.length()
		//console.log('fctor', factor, 'area', area, 'resultarea', area* factor)
		assert(!eq0(factor))
		return {
			area: area * factor,
			centroid: this.matrix.transformPoint(M4.rotateZ(rightLC.angleXY()).transformPoint(new V3(cxt, cyt, 0))),
		}
	}

	at(t: number): V3 {
		assertNumbers(t)
		//assert(this.isValidT(t))
		// = center + f1 cos t + f2 sin t
		return this.center.plus(this.f1.times(Math.cos(t))).plus(this.f2.times(Math.sin(t)))
	}

	tangentAt(t: number): V3 {
		assertNumbers(t)
		//assert(this.isValidT(t))
		// ) f2 cos(t) - f1 sin(t)
		return this.f2.times(Math.cos(t)).minus(this.f1.times(Math.sin(t)))
	}

	ddt(t: number): V3 {
		assertNumbers(t)
		assert(this.isValidT(t))
		return this.f2.times(-Math.sin(t)).minus(this.f1.times(Math.cos(t)))
	}

	tangentAt2(xi: number, eta: number): V3 {
		return this.f2.times(xi).minus(this.f1.times(eta))
	}

	isCircular(): boolean {
		return eq(this.f1.length(), this.f2.length()) && this.f1.isPerpendicularTo(this.f2)
	}

	isColinearTo(curve: Curve): boolean {
		if (!hasConstructor(curve, EllipseCurve)) return false
		if (!this.center.like(curve.center)) {
			return false
		}
		if (this == curve) {
			return true
		}
		if (this.isCircular()) {
			return (
				curve.isCircular() && eq(this.f1.length(), curve.f1.length()) && this.normal.isParallelTo(curve.normal)
			)
		} else {
			let { f1: f1, f2: f2 } = this.rightAngled(),
				{ f1: c1, f2: c2 } = curve.rightAngled()
			if (f1.length() > f2.length()) {
				;[f1, f2] = [f2, f1]
			}
			if (c1.length() > c2.length()) {
				;[c1, c2] = [c2, c1]
			}
			return eq(f1.squared(), Math.abs(f1.dot(c1))) && eq(f2.squared(), Math.abs(f2.dot(c2)))
		}
	}

	pointT(pWC: V3) {
		assertVectors(pWC)
		assert(this.containsPoint(pWC))
		const pLC = this.matrixInverse.transformPoint(pWC)
		const t = EllipseCurve.XYLCPointT(pLC, this.tMin, this.tMax)
		assert(this.isValidT(t))
		return t
	}

	reversed(): EllipseCurve {
		return new EllipseCurve(this.center, this.f1.negated(), this.f2, PI - this.tMax, PI - this.tMin)
	}

	eccentricity() {
		const mainAxes = this.rightAngled()
		const f1length = mainAxes.f1.length(),
			f2length = mainAxes.f1.length()
		const [a, b] = f1length > f2length ? [f1length, f2length] : [f2length, f1length]
		return Math.sqrt(1 - b * b / a / a)
	}

	circumference(): number {
		return this.arcLength(-Math.PI, Math.PI)
	}

	arcLength(tStart: number = this.tMin, tEnd: number = this.tMax, steps: int = 2): number {
		assert(tStart < tEnd, 'startT < endT')
		const f1Length = this.f1.length()
		if (eq(f1Length, this.f2.length())) {
			return f1Length * (tEnd - tStart)
		}
		return super.arcLength(tStart, tEnd, steps)
	}

	circumferenceApproximate(): number {
		// approximate circumference by Ramanujan
		// https://en.wikipedia.org/wiki/Ellipse#Circumference
		const { f1, f2 } = this.rightAngled(),
			a = f1.length(),
			b = f2.length()
		const h = (a - b) ** 2 / (a + b) ** 2
		return Math.PI * (a + b) * (1 + 3 * h / (10 + Math.sqrt(4 - 3 * h)))
	}

	/**
	 * Radii of the ellipse are described by
	 * q(phi) = f1 * cos(phi) + f2 * sin(phi)
	 * or q(xi, eta) = f1 * xi + f2 * eta (1) with the added condition
	 * xi² + eta² = 1 (2)
	 * we want to find the radius where the corresponding tangent is perpendicular
	 * tangent: q'(phi) = f1 * -sin(phi) + f2 * cos(phi)
	 * tangent: q'(xi, eta) = f1 * -eta + f2 * xi
	 * perpendicular when: q'(xi, eta) DOT q(xi, eta) = 0
	 * (f1 * -eta + f2 * xi) DOT (f1 * xi + f2 * eta) = 0
	 * DOT is distributive:
	 * f1² * (-eta * xi) + f1 * f2 * (-eta² + xi²) + f2² * (xi * eta) = 0
	 * (f2² - f1²) * (eta * xi) + f1 * f2 * (-eta² + xi²) = 0
	 * a * (xi² - eta²) + b * xi * eta = 0 (2)
	 * with a = f1 * f2, b = f2² - f1²
	 * => (xi/eta)² + xi/eta * b/a + 1 = 0 (divide by a * eta²)
	 * xi/eta = b/a/2 +- sqrt(b²/a²/4 - 1) | * 2*a*eta
	 * 2 * a * xi = eta * (b +- sqrt(b² - 4 * a²))
	 * g1 * xi - g2 * eta = 0 (3)
	 * with g1 = 2 * a, g2 = b +- sqrt(b² - 4 * a²)
	 * Solve (3), (2) with intersectionUnitCircleLine
	 */
	rightAngled(): EllipseCurve {
		const f1 = this.f1,
			f2 = this.f2,
			a = f1.dot(f2),
			b = f2.squared() - f1.squared()
		if (eq0(a)) {
			return this
		}
		const g1 = 2 * a,
			g2 = b + Math.sqrt(b * b + 4 * a * a)
		const { x1: xi, y1: eta } = intersectionUnitCircleLine(g1, g2, 0)
		const f1RA = f1.times(xi).plus(f2.times(eta))
		const f2RA = f1.times(-eta).plus(f2.times(xi))
		return new EllipseCurve(this.center, f1RA, f2RA, -PI, PI)
	}

	isInfosWithEllipse(ellipse: EllipseCurve): ISInfo[] {
		if (this.normal.isParallelTo(ellipse.normal) && eq0(this.center.minus(ellipse.center).dot(ellipse.normal))) {
			// ellipses are coplanar
			const ellipseLCRA = ellipse.transform(this.matrixInverse).rightAngled()

			const r1 = ellipseLCRA.f1.lengthXY(),
				r2 = ellipseLCRA.f2.lengthXY(),
				centerDist = ellipseLCRA.center.lengthXY()
			const rMin = min(r1, r2),
				rMax = max(r1, r2)
			if (
				lt(centerDist + rMax, 1) || // entirely inside unit circle
				lt(1, centerDist - rMax) || // entirely outside unit circle
				lt(1, rMin - centerDist) || // contains unit circle
				(eq(1, r1) && eq(1, r2) && eq0(centerDist)) // also unit circle, return no IS
			) {
				return []
			}

			const f = (t: number) => ellipseLCRA.at(t).lengthXY() - 1
			const df = (t: number) =>
				ellipseLCRA
					.at(t)
					.xy()
					.dot(ellipseLCRA.tangentAt(t)) / ellipseLCRA.at(t).lengthXY()
			checkDerivate(f, df, -PI, PI, 1)
			const ellipseLCRATs: number[] = []
			for (let startT = -4 / 5 * PI; startT < PI; startT += PI / 4) {
				let t = newtonIterateSmart(f, startT, 16, df, 1e-4)
				le(t, -PI) && (t += TAU)
				assert(!isNaN(t))
				if (between(t, -PI, PI) && eq0(f(t)) && !ellipseLCRATs.some(r => eq(t, r))) {
					ellipseLCRATs.push(t)
				}
			}
			const result: ISInfo[] = []
			for (const ellipseLCRAT of ellipseLCRATs) {
				const p = this.matrix.transformPoint(ellipseLCRA.at(ellipseLCRAT))
				if (this.containsPoint(p) && ellipse.containsPoint(p)) {
					result.push({ tThis: this.pointT(p), tOther: ellipse.pointT(p), p })
				}
			}
			return result

			//const angle = ellipseLCRA.f1.angleXY()
			//const aSqr = ellipseLCRA.f1.squared(), bSqr = ellipseLCRA.f2.squared()
			//const a = Math.sqrt(aSqr), b = Math.sqrt(bSqr)
			//const {x: centerX, y: centerY} = ellipseLCRA.center
			//const rotCenterX = centerX * Math.cos(-angle) + centerY * -Math.sin(-angle)
			//const rotCenterY = centerX * Math.sin(-angle) + centerY * Math.cos(-angle)
			//const rotCenter = V(rotCenterX, rotCenterY)
			//const f = t => {
			//	const lex = Math.cos(t) - rotCenterX, ley = Math.sin(t) - rotCenterY
			//	return lex * lex / aSqr + ley * ley / bSqr - 1
			//}
			//const f2 = (x, y) => (x * x + y * y - 1)
			//const f3 = (x, y) => ((x - rotCenterX) * (x - rotCenterX) / aSqr + (y - rotCenterY) * (y - rotCenterY) /
			// bSqr - 1) const results = [] const resetMatrix = this.matrix.times(M4.rotateZ(angle)) for (let startT =
			// Math.PI / 4; startT < 2 * Math.PI; startT += Math.PI / 2) { const startP = EllipseCurve.XY.at(startT)
			// const p = newtonIterate2d(f3, f2, startP.x, startP.y, 10) if (p && !results.some(r => r.like(p))) {
			// results.push(p) } } const rotEl = new EllipseCurve(rotCenter, V(a, 0, 0), V(0, b, 0)) return
			// results.map(pLC => { const p = resetMatrix.transformPoint(pLC) return {tThis: this.pointT(p, PI),
			// tOther: ellipse.pointT(p, PI), p} })
		} else {
			return this.isTsWithPlane(P3.normalOnAnchor(ellipse.normal.unit(), ellipse.center)).mapFilter(t => {
				const p = this.at(t)
				if (ellipse.containsPoint(p)) {
					return { tThis: t, tOther: ellipse.pointT(p), p }
				}
				return undefined
			})
		}
	}

	isInfosWithCurve(curve: Curve): ISInfo[] {
		if (curve instanceof EllipseCurve) {
			return this.isInfosWithEllipse(curve)
		}
		return super.isInfosWithCurve(curve)
	}

	roots(): [number[], number[], number[]] {
		// tangent(t) = f2 cos t - f1 sin t
		// solve for each dimension separately
		// tangent(eta, xi) = f2 eta - f1 xi

		return arrayFromFunction(3, dim => {
			const a = this.f2.e(dim),
				b = -this.f1.e(dim)
			return intersectionUnitCircleLine2(a, b, 0)
				.map(([xi, eta]) => Math.atan2(eta, xi))
				.filter(t => this.isValidT(t))
		})
	}

	closestTToPoint(p: V3, tStart?: number): number {
		// (at(t) - p) * tangentAt(t) = 0
		// (xi f1 + eta f2 + q) * (xi f2 - eta f1) = 0
		// xi eta (f2^2-f1^2) + xi f2 q - eta² f1 f2 + xi² f1 f2 - eta f1 q = 0
		//  (xi² - eta²) f1 f2 + xi eta (f2^2-f1^2) + xi f2 q - eta f1 q = 0

		// atan2 of p is a good first approximation for the searched t
		tStart = tStart || this.matrixInverse.transformPoint(p).angleXY()
		const pRelCenter = p.minus(this.center)
		const f = (t: number) =>
			this.tangentAt(t).dot(
				this.f1
					.times(Math.cos(t))
					.plus(this.f2.times(Math.sin(t)))
					.minus(pRelCenter),
			)
		return newtonIterate1d(f, tStart, 8)
	}

	area(): number {
		// see
		// https://upload.wikimedia.org/wikipedia/commons/thumb/4/4e/Cross_product_parallelogram.svg/220px-Cross_product_parallelogram.svg.png
		return Math.PI * this.f1.cross(this.f2).length()
	}

	angleToT(phi: number): number {
		// atan2(y, x) = phi
		const phiDir = this.f1
			.unit()
			.times(Math.cos(phi))
			.plus(
				this.f2
					.rejectedFrom(this.f1)
					.unit()
					.times(Math.sin(phi)),
			)
		const dirLC = this.matrixInverse.transformVector(phiDir)
		return dirLC.angleXY()
	}
}

EllipseCurve.prototype.hlol = Curve.hlol++
EllipseCurve.prototype.tIncrement = 2 * Math.PI / (4 * 32)
