import {
	arrayRange,
	assert,
	assertVectors,
	bisect,
	callsce,
	clamp,
	fuzzyUniques,
	int,
	M4,
	newtonIterate2dWithDerivatives,
	V3,
} from 'ts3dutils'

import {
	Curve,
	curvePoint,
	EllipsoidSurface,
	followAlgorithm2d,
	ImplicitCurve,
	ImplicitSurface,
	ISInfo,
	MathFunctionR2R,
	P3,
	ParametricSurface,
	PlaneSurface,
	ProjectedCurveSurface,
	R2_R,
	Surface,
	surfaceIsICurveIsInfosWithLine,
} from '../index'

import { abs, ceil, floor, max, min } from '../math'

export class PICurve extends ImplicitCurve {
	didu: (u: number, v: number) => number
	didv: (u: number, v: number) => number

	constructor(
		points: ReadonlyArray<V3>,
		tangents: ReadonlyArray<V3>,
		readonly parametricSurface: ParametricSurface,
		readonly implicitSurface: ImplicitSurface,
		readonly pmPoints: ReadonlyArray<V3>,
		readonly pmTangents: ReadonlyArray<V3>,
		readonly stepSize: number,
		dir: number = 1,
		generator?: string,
		tMin?: number,
		tMax?: number,
	) {
		super(points, tangents, dir, generator, tMin, tMax)
		assert(Array.isArray(pmPoints))
		assert(dir == 1)
		assert(stepSize <= 1)
		const pf = parametricSurface.pUVFunc()
		const dpdu = parametricSurface.dpdu()
		const dpdv = parametricSurface.dpdv()
		const didp = implicitSurface.didp.bind(implicitSurface)
		this.didu = (u, v) => didp(pf(u, v)).dot(dpdu(u, v))
		this.didv = (u, v) => didp(pf(u, v)).dot(dpdv(u, v))
		for (let i = 0; i < points.length - 1; i++) {
			assert(!points[i].equals(points[i + 1]))
			//assert(parametricSurface.pUV(pmPoints[i].x, pmPoints[i].y).equals(points[i]))
		}
		{
			const ps = this.parametricSurface
			const is = implicitSurface
			const pFunc = ps.pUVFunc(),
				iFunc = is.implicitFunction()
			const dpdu = ps.dpdu()
			const dpdv = ps.dpdv()
			const didp = is.didp.bind(is)
			const mf = MathFunctionR2R.forFFxFy(
				(x, y) => iFunc(pFunc(x, y)),
				(u, v) => didp(pFunc(u, v)).dot(dpdu(u, v)),
				(u, v) => didp(pFunc(u, v)).dot(dpdv(u, v)),
			)
			const { points } = followAlgorithm2d(
				mf,
				this.pmPoints[0],
				stepSize,
				ps,
				(u, v) => is.containsPoint(pFunc(u, v)),
				this.pmPoints.last,
				this.pmTangents[0],
			)
			if (points.length !== this.points.length) {
				followAlgorithm2d(
					mf,
					this.pmPoints[0],
					stepSize,
					ps,
					(u, v) => is.containsPoint(pFunc(u, v)),
					this.pmPoints.last,
					this.pmTangents[0],
				)
			}
			assert(points.length == this.points.length, points.length, this.points.length)
		}
	}

	static forParametricStartEnd(
		ps: ParametricSurface,
		is: ImplicitSurface,
		pmStart: V3,
		pmEnd: V3,
		stepSize: number = 0.02,
		startPMTangent?: V3,
		tMin?: number,
		tMax?: number,
	): PICurve {
		const pFunc = ps.pUVFunc(),
			iFunc = is.implicitFunction()
		const dpdu = ps.dpdu()
		const dpdv = ps.dpdv()
		const didp = is.didp.bind(is)
		const mf = MathFunctionR2R.forFFxFy(
			(x, y) => iFunc(pFunc(x, y)),
			(u, v) => didp(pFunc(u, v)).dot(dpdu(u, v)),
			(u, v) => didp(pFunc(u, v)).dot(dpdv(u, v)),
		)
		const { points, tangents } = followAlgorithm2d(
			mf,
			pmStart,
			stepSize,
			ps,
			(u, v) => is.containsPoint(pFunc(u, v)),
			pmEnd,
			startPMTangent,
		)
		return PICurve.forParametricPointsTangents(ps, is, points, tangents, stepSize, 1, tMin, tMax)
	}

	static forStartEnd(
		ps: ParametricSurface,
		is: ImplicitSurface,
		start: V3,
		end: V3,
		stepSize: number = 0.02,
		startTangent?: V3,
		min?: V3,
		max?: V3,
	): PICurve {
		const startPM = ps.uvP(start)
		const dpdu = ps.dpdu()(startPM.x, startPM.y),
			dpdv = ps.dpdv()(startPM.x, startPM.y)
		const startPMTangent =
			startTangent &&
			M4.forSys(dpdu, dpdv)
				.inversed()
				.transformVector(startTangent)
		// assert(dpdu.times(startPMTangent.x).plus(dpdv.times(startPMTangent.y)).like(startTangent))
		const curve = PICurve.forParametricStartEnd(ps, is, startPM, ps.uvP(end), stepSize, startPMTangent)

		return curve.withBounds(min && curve.pointT(min), max && curve.pointT(max))
	}

	static forParametricPointsTangents(
		ps: ParametricSurface,
		is: ImplicitSurface,
		pmPoints: V3[],
		pmTangents: V3[],
		stepSize: number,
		dir: number = 1,
		tMin?: number,
		tMax?: number,
	): PICurve {
		const pFunc = ps.pUVFunc(),
			dpdu = ps.dpdu()
		const dpdv = ps.dpdv()
		const points = pmPoints.map(({ x, y }) => pFunc(x, y))
		const tangents = pmPoints.map(({ x: u, y: v }, i) => {
			const ds = dpdu(u, v)
			const dt = dpdv(u, v)
			return ds.times(pmTangents[i].x).plus(dt.times(pmTangents[i].y))
			//const p = points[i]
			//return cs.normalP(p).cross(ses.normalP(p))
			//	.toLength(ds.times(pmTangents[i].x).plus(dt.times(pmTangents[i].y)).length())
		})
		return new PICurve(points, tangents, ps, is, pmPoints, pmTangents, stepSize, dir, undefined, tMin, tMax)
	}

	getConstructorParameters(): any[] {
		return [
			this.points,
			this.tangents,
			this.parametricSurface,
			this.implicitSurface,
			this.pmPoints,
			this.pmTangents,
			this.stepSize,
			this.dir,
			this.generator,
		]
	}

	implicitCurve(): R2_R {
		const pF = this.parametricSurface.pUVFunc()
		const iF = this.implicitSurface.implicitFunction()
		return (u, v) => iF(pF(u, v))
	}

	isColinearTo(curve: Curve) {
		if (curve instanceof PICurve) {
			if (this.equals(curve)) {
				return true
			}
			if (
				this.parametricSurface.isCoplanarTo(curve.parametricSurface) &&
				this.implicitSurface.isCoplanarTo(curve.implicitSurface)
			) {
				// TODO
			}
			return false
		} else {
			return false
		}
	}

	containsPoint(p: V3): boolean {
		assertVectors(p)
		const t = this.pointT(p)
		return !isNaN(t) && this.isValidT(t)
	}

	equals(obj: any): boolean {
		return (
			Object.getPrototypeOf(obj) == PICurve.prototype &&
			this.parametricSurface.equals(obj.parametricSurface) &&
			this.implicitSurface.equals(obj.implicitSurface) &&
			this.points[0].equals(obj.points[0]) &&
			this.tangents[0].equals(obj.tangents[0]) &&
			this.dir === obj.dir
		)
	}

	hashCode(): int {
		let hashCode = 0
		hashCode = hashCode * 31 + this.parametricSurface.hashCode()
		hashCode = hashCode * 31 + this.implicitSurface.hashCode()
		hashCode = hashCode * 31 + this.points[0].hashCode()
		hashCode = hashCode * 31 + this.tangents[0].hashCode()
		return hashCode | 0
	}

	tangentP(point: V3): V3 {
		assertVectors(point)
		assert(this.containsPoint(point), 'this.containsPoint(point)')
		const t = this.pointT(point)
		return this.tangentAt(t)
	}

	tangentAt(t: number): V3 {
		assert(!isNaN(t))
		if (0 === t % 1) return this.tangents[t]
		const uv = this.uvT(t)
		const uvTangent = new V3(-this.didv(uv.x, uv.y), this.didu(uv.x, uv.y), 0).toLength(this.stepSize)
		const du = this.parametricSurface.dpdu()(uv.x, uv.y)
		const dv = this.parametricSurface.dpdv()(uv.x, uv.y)
		return du.times(uvTangent.x).plus(dv.times(uvTangent.y))
	}

	at(t: number): V3 {
		assert(!isNaN(t))
		if (0 === t % 1) return this.points[t]
		const startParams = V3.lerp(this.pmPoints[floor(t)], this.pmPoints[ceil(t)], t % 1)
		return this.closestPointToParams(startParams)
	}

	uvT(t: number): V3 {
		assert(!isNaN(t))
		if (0 === t % 1) return this.pmPoints[t]
		const startParams = V3.lerp(this.pmPoints[floor(t)], this.pmPoints[ceil(t)], t % 1)
		return curvePoint(this.implicitCurve(), startParams, this.didu, this.didv)
	}

	closestTToPoint(p: V3, tStart?: number): number {
		// TODO
		return 0
	}

	closestPointToParams(startUV: V3): V3 {
		const pointParams = curvePoint(this.implicitCurve(), startUV, this.didu, this.didv)
		return this.parametricSurface.pUVFunc()(pointParams.x, pointParams.y)
	}

	isTsWithSurface(surface: Surface): number[] {
		if (surface instanceof EllipsoidSurface) {
			const pS = this.parametricSurface,
				iS = this.implicitSurface
			if (pS instanceof ProjectedCurveSurface && iS instanceof EllipsoidSurface) {
				const iscs = iS.isCurvesWithSurface(surface)
				const points = iscs.flatMap(isc => isc.isTsWithSurface(pS).map(t => isc.at(t)))
				const ts = fuzzyUniques(points.map(p => this.pointT(p)))
				return ts.filter(t => !isNaN(t) && this.isValidT(t))
			}
		} else if (ImplicitSurface.is(surface)) {
			const result: number[] = []
			const iF = surface.implicitFunction()
			let prevSignedDistance = iF(this.points[0])
			for (let i = 1; i < this.points.length; i++) {
				const point = this.points[i]
				const signedDistance = iF(point)
				if (prevSignedDistance * signedDistance <= 0) {
					const pF = this.parametricSurface.pUVFunc()
					const dpdu = this.parametricSurface.dpdu()
					const dpdv = this.parametricSurface.dpdv()
					const startUV = this.pmPoints[abs(prevSignedDistance) < abs(signedDistance) ? i - 1 : i]
					const isUV = newtonIterate2dWithDerivatives(
						this.implicitCurve(),
						(u, v) => iF(pF(u, v)),
						startUV.x,
						startUV.y,
						4,
						this.didu,
						this.didv,
						(u, v) => dpdu(u, v).dot(surface.didp(pF(u, v))),
						(u, v) => dpdv(u, v).dot(surface.didp(pF(u, v))),
					)!
					result.push(this.pointT(this.parametricSurface.pUV(isUV.x, isUV.y)))
				}
				prevSignedDistance = signedDistance
			}
			return result
		}
		throw new Error()
	}

	isTsWithPlane(planeWC: P3): number[] {
		return this.isTsWithSurface(new PlaneSurface(planeWC))
		// version which intersects the plane with the defining surfaces of this PICurve, but this causes
		// issues when they are PICurves too:
		// assertInst(P3, planeWC)
		// const ps = this.parametricSurface,
		// 	is = this.implicitSurface
		// const pscs = ps.isCurvesWithPlane(planeWC)
		// const iscs = is.isCurvesWithPlane(planeWC)
		// const infos = iscs.flatMap(isc => pscs.flatMap(psc => isc.isInfosWithCurve(psc)))
		// const ts = fuzzyUniques(infos.map(info => this.pointT(info.p)))
		// return ts.filter(t => !isNaN(t) && this.isValidT(t))
	}

	pointT(p: V3): number {
		assertVectors(p)
		if (!this.parametricSurface.containsPoint(p) || !this.implicitSurface.containsPoint(p)) {
			return NaN
		}
		const pmPoint = this.parametricSurface.uvPFunc()(p)
		const ps = this.points,
			pmps = this.pmPoints
		let t = 0,
			pmDistance = pmPoint.distanceTo(pmps[0])
		while (pmDistance > abs(this.stepSize) && t < ps.length - 1) {
			// TODO -1?
			//console.log(t, pmps[t].$, pmDistance)
			t = min(pmps.length - 1, t + max(1, Math.round(pmDistance / abs(this.stepSize) / 2 / 2)))
			pmDistance = pmPoint.distanceTo(pmps[t])
		}
		// if (t < this.pmPoints.length - 1 && pmDistance > pmPoint.distanceTo(pmps[t + 1])) {
		//     t++
		// }
		if (pmDistance > abs(this.stepSize) * 1.1) {
			// p is not on this curve
			return NaN
		}
		if (t == ps.length - 1) {
			t--
		}
		if (ps[t].like(p)) return t
		if (ps[t + 1].like(p)) return t + 1
		const startT = arrayRange(floor(this.tMin), ceil(this.tMax), 1).withMax(t => -pmPoint.distanceTo(pmps[t]))
		if (undefined === startT) throw new Error()
		if (ps[startT].like(p)) return startT
		//const [a, b] = 0 === startT
		//    ? [0, 1]
		//    : this.points.length - 1 === startT
		//        ? [startT - 1, startT]
		//        : pmPoint.distanceTo(pmps[startT - 1]) < pmPoint.distanceTo(pmps[startT + 1])
		//            ? [startT - 1, startT]
		//            : [startT, startT + 1]
		const a = max(0, startT - 1),
			b = min(this.points.length - 1, startT + 1)
		const tangent = this.tangentAt(startT)
		const f = (t: number) =>
			this.at(clamp(t, 0, this.points.length - 1))
				.to(p)
				.dot(tangent)
		// const df = (t: number) => -this.tangentAt(clamp(t, 0, this.points.length - 1)).dot(tangent)
		//checkDerivate(f, df, 0, this.points.length - 2, 3)
		// 8 steps necessary because df can currently be way off
		t = bisect(f, a, b, 32)
		if (!isFinite(t) || this.at(t).distanceTo(p) > abs(this.stepSize)) {
			return NaN
		}
		return t
	}

	transform(m4: M4): this {
		const dirFactor = m4.isMirroring() ? -1 : 1
		return PICurve.forStartEnd(
			this.parametricSurface.transform(m4),
			this.implicitSurface.transform(m4),
			m4.transformPoint(this.points[0]),
			m4.transformPoint(this.points.last),
			this.stepSize * dirFactor,
			m4.transformVector(this.tangents[0]),
			m4.transformPoint(this.at(this.tMin)),
			m4.transformPoint(this.at(this.tMax)),
		) as this
		//return PICurve.forParametricStartEnd(
		//	this.parametricSurface.transform(m4),
		//	this.implicitSurface.transform(m4),
		//	this.pmPoints[0],
		//	this.pmPoints.last,
		//	this.stepSize,
		//	this.dir,
		//	this.tMin,
		//	this.tMax)
		// TODO: pass transformed points?
		//return new PICurve(
		//	m4.transformedPoints(this.points),
		//	m4.transformedVectors(this.tangents),
		//    this.parametricSurface.transform(m4),
		//   this.implicitSurface.transform(m4),
		//   this.pmPoints,
		//   this.pmTangents,
		//this.stepSize,
		//   this.dir,
		//this.generator,
		//this.tMin, this.tMax)
	}

	roots(): [number[], number[], number[]] {
		const allTs = arrayRange(0, this.points.length)
		return [allTs, allTs, allTs]
	}

	isInfosWithLine(
		anchorWC: V3,
		dirWC: V3,
		tMin?: number | undefined,
		tMax?: number | undefined,
		lineMin?: number | undefined,
		lineMax?: number | undefined,
	): ISInfo[] {
		return surfaceIsICurveIsInfosWithLine.call(this, anchorWC, dirWC, tMin, tMax, lineMin, lineMax)
	}

	toSource(rounder: (x: number) => number = x => x): string {
		const result = callsce(
			'PICurve.forParametricStartEnd',
			this.parametricSurface,
			this.implicitSurface,
			this.pmPoints[0],
			this.pmPoints.last,
			this.stepSize,
			this.pmTangents[0],
			this.tMin,
			this.tMax,
		)
		return result
	}
}

PICurve.prototype.tIncrement = 1
