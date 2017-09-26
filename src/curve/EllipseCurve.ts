import { Equalable } from 'javasetmap.ts'
import { V3, assertNumbers, assert, Transformable, le, ge, arrayFromFunction, newtonIterateWithDerivative, NLA_PRECISION, int, callsce, eq, fuzzyUniquesF, clamp, AABB, glqInSteps, M4, newtonIterate2dWithDerivatives, V, eq0, getIntervals, assertf } from 'ts3dutils'
import { followAlgorithm2d } from '../B2'
import { P3 } from '../P3'
import { Surface } from '../surface/Surface'
import {XiEtaCurve} from './XiEtaCurve'

export class EllipseCurve extends XiEtaCurve {
	constructor(center: V3, f1: V3, f2: V3, tMin: number = -PI, tMax: number = PI) {
		super(center, f1, f2, tMin, tMax)
		assert(EllipseCurve.isValidT(tMin))
		assert(EllipseCurve.isValidT(tMax))
	}

	// TODO: there'S alsoa commented out test
	getVolZAnd(dir1: V3, tStart: number, tEnd: number): {volume: number, centroid: V3} {
		// let p = at(t)
		// integrate area [p -> plane.projectPoint(p)] to x axis...
		// INTEGRATE[tStart, tEnd] fp(this.at(t)) dt
		const ft = t => this.center.plus(this.f1.times(Math.cos(t))).plus(this.f2.times(Math.sin(t)))
		// f(t) = c + f1 cos + f2 sin
		// p dot d1 = (cx + f1x cos + f2x sin) dx + (cy + f1y cos + f2y sin) dy + (cz + f1z cos + f2z sin) dz
		function fp(p: V3) {
			const p0ToP = dir1.times(dir1.dot(p))
			const area = p0ToP.lengthXY() * (p.z - p0ToP.z / 2)
			return area
		}
		const f = (t: number) => fp(this.at(t)) * this.tangentAt(t).cross(this.normal).unit().z
		return {volume: glqInSteps(f, tStart, tEnd, 4), centroid: undefined}
	}

	getAreaInDir(right: V3, up: V3, tStart: number, tEnd: number): {area: number, centroid: V3} {
		//assertf(() => tStart < tEnd)
		assertf(() => right.isPerpendicularTo(this.normal))
		assertf(() => up.isPerpendicularTo(this.normal))
		//assertf(() => EllipseCurve.isValidT(tStart), tStart)
		//assertf(() => EllipseCurve.isValidT(tEnd), tEnd)

		const upLC = this.inverseMatrix.transformVector(up)
		const rightLC = upLC.cross(V3.Z)
		const normTStart = tStart - rightLC.angleXY()
		const normTEnd = tEnd - rightLC.angleXY()
		const transformedOriginY = this.inverseMatrix.getTranslation().dot(upLC.unit())
		//console.log(upLC.str, rightLC.str, normTStart, normTEnd, 'upLC.length()', upLC.length())
		//console.log('transformedOriginY', transformedOriginY)
		//assertf(() => upLC.hasLength(1), upLC.length())
		const fPi = Math.PI / 4
		// integral of sqrt(1 - x²) from 0 to cos(t)
		// Basically, we want
		// INTEGRAL[cos(t); PI/2] sqrt(1 - x²) dx
		// INTEGRAL[PI/2: cos(t)] -sqrt(1 - x²) dx
		// = INTEGRAL[cos(0); cos(t)] -sqrt(1 - x²) dx
		// = INTEGRAL[0; t] -sqrt(1 - cos²(t)) * -sin(t) dt
		// = INTEGRAL[0; t] -sin(t) * -sin(t) dt
		// = INTEGRAL[0; t] sin²(t) dt (partial integration / wolfram alpha)
		// = (1/2 * (t - sin(t) * cos(t)))[0; t] (this form has the distinct advantage of being defined everywhere)
		function fArea(t) { return (t - Math.sin(t) * Math.cos(t)) / 2 }

		// for the centroid, we want
		// cx = 1 / area * INTEGRAL[cos(t); PI/2] x * f(x) dx
		// cx = 1 / area * INTEGRAL[cos(t); PI/2] x * sqrt(1 - x²) dx
		// cx = 1 / area * INTEGRAL[cos(0); cos(t)] x * -sqrt(1 - x²) dx
		// ...
		// cx = 1 / area * INTEGRAL[0; t] cos(t) * sin²(t) dt // WA
		// cx = 1 / area * (sin^3(t) / 3)[0; t]
		function cxTimesArea(t) { return Math.pow(Math.sin(t), 3) / 3 }

		// cy = 1 / area * INTEGRAL[cos(t); PI/2] f²(x) / 2 dx
		// cy = 1 / area * INTEGRAL[cos(0); cos(t)] -(1 - x²) / 2 dx
		// cy = 1 / area * INTEGRAL[0; t] (cos²(t) - 1) * -sin(t) / 2 dt
		// cy = 1 / area * (cos (3 * t) - 9 * cos(t)) / 24 )[0; t]
		function cyTimesArea(t) { return (Math.cos(3 * t) - 9 * Math.cos(t)) / 24 }

		const restArea = -transformedOriginY * (-Math.cos(normTEnd) + Math.cos(normTStart) )
		const area = fArea(normTEnd) - fArea(normTStart) + restArea
		const cxt = (cxTimesArea(normTEnd) - cxTimesArea(normTStart) + -transformedOriginY * (-Math.cos(normTEnd) - Math.cos(normTStart)) / 2 * restArea) / area
		const cyt = (cyTimesArea(normTEnd) - cyTimesArea(normTStart) - -transformedOriginY / 2 * restArea) / area
		const factor = this.matrix.xyAreaFactor() // * upLC.length()
		//console.log('fctor', factor, 'area', area, 'resultarea', area* factor)
		assert(!eq0(factor))
		return {area: area * factor, centroid: this.matrix.transformPoint(M4.rotateZ(rightLC.angleXY()).transformPoint(new V3(cxt, cyt, 0)))}

	}

	static isValidT(t) {
		return -Math.PI <= t && t <= Math.PI
	}

	at(t: number): V3 {
		// = center + f1 cos t + f2 sin t
		return this.center.plus(this.f1.times(Math.cos(t))).plus(this.f2.times(Math.sin(t)))
	}

	tangentAt(t: number) {
		assertNumbers(t)
		// f2 cos(t) - f1 sin(t)
		return this.f2.times(Math.cos(t)).minus(this.f1.times(Math.sin(t)))
	}

	ddt(t: number): V3 {
		assertNumbers(t)
		return this.f2.times(-Math.sin(t)).minus(this.f1.times(Math.cos(t)))
	}

	tangentAt2(xi: number, eta: number): V3 {
		return this.f2.times(xi).minus(this.f1.times(eta))
	}

	isCircular() {
		return eq(this.f1.length(), this.f2.length()) && this.f1.isPerpendicularTo(this.f2)
	}

	reversed(): this {
		return new this.constructor(this.center, this.f1, this.f2.negated(), -this.tMax, -this.tMin) as this
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
			return curve.isCircular() && eq(this.f1.length(), curve.f1.length()) && this.normal.isParallelTo(curve.normal)
		} else {
			let {f1: f1, f2: f2} = this.rightAngled(), {f1: c1, f2: c2} = curve.rightAngled()
			if (f1.length() > f2.length()) {[f1, f2] = [f2, f1]}
			if (c1.length() > c2.length()) {[c1, c2] = [c2, c1]}
			return eq(f1.squared(), Math.abs(f1.dot(c1)))
				&& eq(f2.squared(), Math.abs(f2.dot(c2)))
		}
	}

	static XYLCValid(pLC: V3): boolean {
		return eq(1, pLC.lengthXY())
	}

	/**
	 * @param hint +-PI, whichever is correct
	 */
	static XYLCPointT(pLC: V3, hint?: number): number {
		const angle = pLC.angleXY()
		if (angle < -Math.PI + NLA_PRECISION || angle > Math.PI - NLA_PRECISION) {
			assert(isFinite(hint))
			return Math.sign(hint) * Math.PI
		}
		return angle
	}

	eccentricity() {
		const mainAxes = this.rightAngled()
		const f1length = mainAxes.f1.length(), f2length = mainAxes.f1.length()
		const [a, b] = f1length > f2length ? [f1length, f2length] : [f2length, f1length]
		return Math.sqrt(1 - b * b / a / a)
	}

	circumference(): number {
		return this.arcLength(-Math.PI, Math.PI)
	}

	arcLength(startT: number, endT: number, steps?: int): number {
		assert(startT < endT, 'startT < endT')
		if (this.isCircular()) {
			return this.f1.length() * (endT - startT)
		}
		return super.arcLength(startT, endT, steps)
	}

	circumferenceApproximate(): number {
		// approximate circumference by Ramanujan
		// https://en.wikipedia.org/wiki/Ellipse#Circumference
		const {f1, f2} = this.rightAngled(), a = f1.length(), b = f2.length()
		const h = (a - b) ** 2 / (a + b) ** 2
		return Math.PI * (a + b) * (1 + 3 * h / (10 + Math.sqrt(4 - 3 * h)))
	}

	rightAngled(): EllipseCurve {
		const f1 = this.f1, f2 = this.f2, a = f1.dot(f2), b = f2.squared() - f1.squared()
		if (eq0(a)) {
			return this
		}
		const g1 = 2 * a, g2 = b + Math.sqrt(b * b + 4 * a * a)
		const {x1: xi, y1: eta} = intersectionUnitCircleLine(g1, g2, 0)
		return new EllipseCurve(this.center,
			f1.times(xi).plus(f2.times(eta)),
			f1.times(-eta).plus(f2.times(xi)))
	}

	static magic(a: number, b: number, c: number): number[] {
		const isLC = intersectionUnitCircleLine2(a, b, c)
		return isLC.map(([xi, eta]) => Math.atan2(eta, xi))
	}

	isInfosWithEllipse(ellipse: EllipseCurve): ISInfo[] {
		if (this.normal.isParallelTo(ellipse.normal) && eq0(this.center.minus(ellipse.center).dot(ellipse.normal))) {

			// ellipses are coplanar
			const ellipseLCRA = ellipse.transform(this.inverseMatrix).rightAngled()

			const r1 = ellipseLCRA.f1.lengthXY(), r2 = ellipseLCRA.f2.lengthXY(), centerDist = ellipseLCRA.center.lengthXY()
			const rMin = min(r1, r2), rMax = max(r1, r2)
			if (lt(centerDist + rMax, 1) || // entirely inside unit circle
				lt(1, centerDist - rMax) || // entirely outside unit circle
				lt(1, rMin - centerDist) || // contains unit circle
				eq(1, r1) && eq(1, r2) && eq0(centerDist) // also unit circle, return no IS
			) {
				return []
			}

			const f = (t: number) => ellipseLCRA.at(t).lengthXY() - 1
			const df =(t: number) => ellipseLCRA.at(t).xy().dot(ellipseLCRA.tangentAt(t)) / ellipseLCRA.at(t).lengthXY()
			checkDerivate(f, df, -PI, PI, 1)
			const ts: number[] = []
			const tsvs = arrayRange(-4/5 * PI, PI, PI/4).map(startT => [startT, df(startT), newtonIterateSmart(f, startT, 16, df, 1e-4), f(newtonIterateSmart(f, startT, 16, df, 1e-4))])
			for (let startT = -4/5 * PI; startT < PI; startT += PI / 4) {
				let t = newtonIterateSmart(f, startT, 16, df, 1e-4)
				le(t, -PI) && (t += TAU)
				assert(!isNaN(t))
				if (ellipseLCRA.isValidT(t) && eq0(f(t)) && !ts.some(r => eq(t, r))) {
					ts.push(t)
				}
			}
			return ts.map(raT => {
				const p = this.matrix.transformPoint(ellipseLCRA.at(raT))
				return {tThis: this.pointT(p), tOther: ellipse.pointT(p, PI), p}
			})

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
			//const f3 = (x, y) => ((x - rotCenterX) * (x - rotCenterX) / aSqr + (y - rotCenterY) * (y - rotCenterY) / bSqr - 1)
			//const results = []
			//const resetMatrix = this.matrix.times(M4.rotateZ(angle))
			//for (let startT = Math.PI / 4; startT < 2 * Math.PI; startT += Math.PI / 2) {
			//	const startP = EllipseCurve.XY.at(startT)
			//	const p = newtonIterate2d(f3, f2, startP.x, startP.y, 10)
			//	if (p && !results.some(r => r.like(p))) {
			//		results.push(p)
			//	}
			//}
			//const rotEl = new EllipseCurve(rotCenter, V(a, 0, 0), V(0, b, 0))
			//return results.map(pLC => {
			//	const p = resetMatrix.transformPoint(pLC)
			//	return {tThis: this.pointT(p, PI), tOther: ellipse.pointT(p, PI), p}
			//})
		} else {
			return this.isTsWithPlane(ellipse.getPlane()).mapFilter(t => {
				const p = this.at(t)
				if (ellipse.containsPoint(p)) {
					return {tThis: t, tOther: ellipse.pointT(p), p}
				}
			})
		}
	}

	static unitIsInfosWithLine(anchorLC: V3, dirLC: V3, anchorWC: V3, dirWC: V3): ISInfo[] {
		// ell: x² + y² = 1 = p²
		// line(t) = anchor + t dir
		// anchor² - 1 + 2 t dir anchor + t² dir² = 0
		const pqDiv = dirLC.dot(dirLC)
		const lineTs = pqFormula(2 * dirLC.dot(anchorLC) / pqDiv, (anchorLC.dot(anchorLC) - 1) / pqDiv)
		return lineTs.map(tOther => ({
			tThis: Math.atan2(anchorLC.y + tOther * dirLC.y, anchorLC.x + tOther * dirLC.x),
			tOther: tOther,
			p: L3.at(anchorWC, dirWC, tOther)}))
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
			const a = this.f2.e(dim), b = -this.f1.e(dim)
			const {x1,y1,x2,y2} = intersectionUnitCircleLine(a, b, 0)
			return [Math.atan2(y1, x1), Math.atan2(y2, x2)]
		})
	}

	closestTToPoint(p: V3, tStart?: number): number {
		// (at(t) - p) * tangentAt(t) = 0
		// (xi f1 + eta f2 + q) * (xi f2 - eta f1) = 0
		// xi eta (f2^2-f1^2) + xi f2 q - eta² f1 f2 + xi² f1 f2 - eta f1 q = 0
		//  (xi² - eta²) f1 f2 + xi eta (f2^2-f1^2) + xi f2 q - eta f1 q = 0

		// atan2 of p is a good first approximation for the searched t
		const startT = this.inverseMatrix.transformPoint(p).angleXY()
		const pRelCenter = p.minus(this.center)
		const f = (t: number) => this.tangentAt(t).dot(this.f1.times(Math.cos(t)).plus(this.f2.times(Math.sin(t))).minus(pRelCenter))
		return newtonIterate1d(f, startT)
	}

	/**
	 * Returns a new EllipseCurve representing a circle parallel to the XY-plane.`
	 */
	static circle(radius: number, center: V3 = V3.O): EllipseCurve {
		return new EllipseCurve(center, new V3(radius, 0, 0), new V3(0, radius, 0))
	}

	area(): number {
		// see
		// https://upload.wikimedia.org/wikipedia/commons/thumb/4/4e/Cross_product_parallelogram.svg/220px-Cross_product_parallelogram.svg.png
		return Math.PI * this.f1.cross(this.f2).length()
	}

	angleToT(phi: number): number {
		// atan2(y, x) = phi
		const phiDir = this.f1.unit().times(Math.cos(phi)).plus(this.f2.rejectedFrom(this.f1).unit().times(Math.sin(phi)))
		const dirLC = this.inverseMatrix.transformVector(phiDir)
		return dirLC.angleXY()
	}

    static readonly XY = new EllipseCurve(V3.O, V3.X, V3.Y)

	static circleForCenter2P(center: V3, a: V3, b: V3, radius: number) {
		const f1 = center.to(a)
		const normal = f1.cross(center.to(b))
		const f2 = normal.cross(f1).toLength(f1.length())
		const tMax = f1.angleTo(center.to(b))
		return new EllipseCurve(center, f1, f2, 0, tMax)
	}
}
EllipseCurve.prototype.hlol = Curve.hlol++
EllipseCurve.prototype.tIncrement = 2 * Math.PI / (4 * 800)