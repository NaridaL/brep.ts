import { Equalable } from 'javasetmap.ts'
import {
	AABB,
	arrayFromFunction,
	assert,
	assertNumbers,
	callsce,
	clamp,
	eq,
	eq0,
	fuzzyUniquesF,
	getIntervals,
	glqInSteps,
	hasConstructor,
	int,
	le,
	M4,
	newtonIterate1d,
	newtonIterate2dWithDerivatives,
	newtonIterateWithDerivative,
	NLA_PRECISION,
	Transformable,
	V,
	V3,
} from 'ts3dutils'

import {
	curvePointPP,
	CylinderSurface,
	EllipsoidSurface,
	followAlgorithm2d,
	followAlgorithmPP,
	ISInfo,
	MathFunctionR2R,
	P3,
	ParametricSurface,
	PlaneSurface,
	PPCurve,
	ProjectedCurveSurface,
	Surface,
} from '../index'

import { abs, ceil, floor } from '../math'

/**
 * Information about the intersection of two curves.
 */
export interface ISInfo {
	/** curve parameter of the "this" curve */
	tThis: number
	/** curve parameter of the "other" curve */
	tOther: number
	/** intersection point */
	p: V3
}

let insideIsInfosWithCurve = false

export interface Curve {
	/**
	 * Derivative of tangentAt for parameter t at t.
	 */
	ddt?(t: number): V3
}
export abstract class Curve extends Transformable implements Equalable {
	static hlol = 0
	tIncrement!: number
	hlol!: number;
	readonly ['constructor']: new (...args: any[]) => this

	constructor(readonly tMin: number, readonly tMax: number) {
		super()
		assertNumbers(tMin, tMax)
		assert('number' == typeof tMin && !isNaN(tMin))
		assert('number' == typeof tMax && !isNaN(tMax))
		assert(tMin < tMax, 'tMin < tMax ' + tMin + ' < ' + tMax)
	}

	static integrate(curve: Curve, startT: number, endT: number, steps: int): number {
		const step = (endT - startT) / steps
		let length = 0
		let p = curve.at(startT)
		let i = 0,
			t = startT + step
		for (; i < steps; i++, t += step) {
			const next = curve.at(t)
			length += p.distanceTo(next)
			p = next
		}
		return length
	}

	static ispsRecursive(
		curve1: Curve,
		tMin: number,
		tMax: number,
		curve2: Curve,
		sMin: number,
		sMax: number,
	): ISInfo[] {
		// the recursive function finds good approximates for the intersection points
		// curve1 function uses newton iteration to improve the result as much as possible
		function handleStartTS(startT: number, startS: number) {
			if (!result.some(info => eq(info.tThis, startT) && eq(info.tOther, startS))) {
				const f1 = (t: number, s: number) => curve1.tangentAt(t).dot(curve1.at(t).minus(curve2.at(s)))
				const f2 = (t: number, s: number) => curve2.tangentAt(s).dot(curve1.at(t).minus(curve2.at(s)))
				// f = (b1, b2, t1, t2) = b1.tangentAt(t1).dot(b1.at(t1).minus(b2.at(t2)))
				const dfdt1 = (b1: Curve, b2: Curve, t1: number, t2: number) =>
					b1.ddt!(t1).dot(b1.at(t1).minus(b2.at(t2))) + b1.tangentAt(t1).squared()
				const dfdt2 = (b1: Curve, b2: Curve, t1: number, t2: number) => -b1.tangentAt(t1).dot(b2.tangentAt(t2))
				const ni = newtonIterate2dWithDerivatives(
					f1,
					f2,
					startT,
					startS,
					16,
					dfdt1.bind(undefined, curve1, curve2),
					dfdt2.bind(undefined, curve1, curve2),
					(t, s) => -dfdt2(curve2, curve1, s, t),
					(t, s) => -dfdt1(curve2, curve1, s, t),
				)!
				assert(isFinite(ni.x))
				assert(isFinite(ni.y))
				if (ni == undefined) console.log(startT, startS, curve1.sce, curve2.sce)
				result.push({ tThis: ni.x, tOther: ni.y, p: curve1.at(ni.x) })
			}
		}

		// returns whether an intersection was immediately found (i.e. without further recursion)
		function findRecursive(
			tMin: number,
			tMax: number,
			sMin: number,
			sMax: number,
			curve1AABB: AABB,
			curve2AABB: AABB,
			depth = 0,
		) {
			const EPS = NLA_PRECISION
			if (curve1AABB.touchesAABBfuzzy(curve2AABB)) {
				const tMid = (tMin + tMax) / 2
				const sMid = (sMin + sMax) / 2
				if (Math.abs(tMax - tMin) < EPS || Math.abs(sMax - sMin) < EPS) {
					handleStartTS(tMid, sMid)
					return true
				} else {
					const curve1AABBleft = curve1.getAABB(tMin, tMid)
					const curve2AABBleft = curve2.getAABB(sMin, sMid)
					let curve1AABBright, curve2AABBright
					// if one of the following calls immediately finds an intersection, we don't want to call the others
					// as that will lead to the same intersection being output multiple times
					findRecursive(tMin, tMid, sMin, sMid, curve1AABBleft, curve2AABBleft, depth + 1) ||
						findRecursive(
							tMin,
							tMid,
							sMid,
							sMax,
							curve1AABBleft,
							(curve2AABBright = curve2.getAABB(sMid, sMax)),
							depth + 1,
						) ||
						findRecursive(
							tMid,
							tMax,
							sMin,
							sMid,
							(curve1AABBright = curve1.getAABB(tMid, tMax)),
							curve2AABBleft,
							depth + 1,
						) ||
						findRecursive(tMid, tMax, sMid, sMax, curve1AABBright, curve2AABBright, depth + 1)
				}
			}
			return false
		}

		const result: ISInfo[] = []
		findRecursive(tMin, tMax, sMin, sMax, curve1.getAABB(tMin, tMax), curve2.getAABB(sMin, sMax))
		return fuzzyUniquesF(result, info => info.tThis)
	}

	/**
	 * Searches a 2d area for (an) implicit curve(s).
	 * @param implicitCurve
	 * @param bounds Defines area to search.
	 * @param uStep Granularity of search in s-direction.
	 * @param vStep Granularity of search in t-direction.
	 * @param stepSize step size to take along the curve
	 * @return
	 */
	static breakDownIC(
		implicitCurve: MathFunctionR2R,
		bounds: AABB2,
		uStep: number,
		vStep: number,
		stepSize: number,
		validUV: R2<boolean>,
	): { points: V3[]; tangents: V3[] }[] {
		//undefined == didu && (didu = (u, v) => (implicitCurve(u + EPS, v) - implicitCurve(u, v)) / EPS)
		//undefined == didv && (didv = (u, v) => (implicitCurve(u, v + EPS) - implicitCurve(u, v)) / EPS)

		const { uMin, uMax, vMin, vMax } = bounds
		const deltaS = uMax - uMin,
			deltaT = vMax - vMin
		const sRes = ceil(deltaS / uStep),
			tRes = ceil(deltaT / vStep)
		const grid = new Array(sRes * tRes).fill(0)
		// const printGrid = () =>
		// 	console.log(
		// 		arrayFromFunction(tRes, i =>
		// 			grid
		// 				.slice(sRes * i, sRes * (i + 1))
		// 				.map(v => (v ? 'X' : '_'))
		// 				.join(''),
		// 		).join('\n'),
		// 	)
		const get = (i: int, j: int) => grid[j * sRes + i]
		const set = (i: int, j: int) => 0 <= i && i < sRes && 0 <= j && j < tRes && (grid[j * sRes + i] = 1)
		const result: { points: V3[]; tangents: V3[] }[] = []
		const logTable = []
		for (let i = 0; i < sRes; i++) {
			search: for (let j = 0; j < tRes; j++) {
				if (get(i, j)) continue
				set(i, j)
				let u = uMin + (i + 0.5) * uStep,
					v = vMin + (j + 0.5) * vStep
				const startS = u,
					startT = v
				// basically curvePoint
				for (let k = 0; k < 8; k++) {
					const fp = implicitCurve(u, v)
					const dfpdx = implicitCurve.x(u, v),
						dfpdy = implicitCurve.y(u, v)
					if (0 === dfpdx ** 2 + dfpdy ** 2) {
						// top of a hill, keep looking
						continue search
					}
					const scale = fp / (dfpdx ** 2 + dfpdy ** 2)
					u -= scale * dfpdx
					v -= scale * dfpdy
				}
				const li = floor((u - uMin) / uStep),
					lj = floor((v - vMin) / vStep)
				logTable.push({
					i,
					j,
					li,
					lj,
					startS,
					startT,
					u,
					v,
					'bounds(u, v)': uvInAABB2(bounds, u, v),
					'ic(s,t)': implicitCurve(u, v),
				})
				if (!(i == li && j == lj) && get(li, lj)) {
					continue search
				}
				set(li, lj)
				// u, v are now good starting coordinates to use follow algorithm
				if (uvInAABB2(bounds, u, v) && validUV(u, v) && eq0(implicitCurve(u, v))) {
					const subResult = mkcurves(implicitCurve, u, v, stepSize, bounds, validUV)
					for (const curveData of subResult) {
						assert(curveData.points.length > 2)
						for (const { x, y } of curveData.points) {
							const lif = (x - uMin) / uStep,
								ljf = (y - vMin) / vStep
							set((lif - 0.5) | 0, (ljf - 0.5) | 0)
							set((lif - 0.5) | 0, (ljf + 0.5) | 0)
							set((lif + 0.5) | 0, (ljf - 0.5) | 0)
							set((lif + 0.5) | 0, (ljf + 0.5) | 0)
						}
					}
					//printGrid()
					result.push(...subResult)
				}
			}
		}
		// console.table(logTable)
		for (const { points } of result) {
			for (let i = 0; i < points.length - 1; i++) {
				assert(!points[i].equals(points[i + 1]))
			}
		}
		return result
	}

	toString() {
		return this.toSource()
	}

	toSource(rounder: (x: number) => number = x => x): string {
		return callsce.call(
			undefined,
			'new ' + this.constructor.name,
			...this.getConstructorParameters(),
			this.tMin,
			this.tMax,
		)
	}

	/**
	 * Excludes tMin, tMax.
	 */
	abstract getConstructorParameters(): any[]

	withBounds<T extends Curve>(this: T, tMin = this.tMin, tMax = this.tMax): T {
		//assert(this.tMin <= tMin && tMin <= this.tMax)
		//assert(this.tMin <= tMax && tMax <= this.tMax)
		return new this.constructor(...this.getConstructorParameters(), tMin, tMax)
	}

	/**
	 * Curve parameter t for point p on curve.
	 */
	abstract pointT(p: V3, tMin?: number, tMax?: number): number

	/**
	 * The point on the line that is closest to the given point.
	 */
	closestPointToPoint(p: V3): V3 {
		return this.at(this.closestTToPoint(p))
	}

	isValidT(t: number): boolean {
		return le(this.tMin, t) && le(t, this.tMax)
	}

	diff(t: number, eps: number): V3 {
		return this.at(t).to(this.at(t + eps))
	}

	// TODO: tmin/tmax first
	closestTToPoint(p: V3, tStart?: number, tMin = this.tMin, tMax = this.tMax): number {
		// this.at(t) has minimal distance to p when this.tangentAt(t) is perpendicular to
		// the vector between this.at(t) and p. This is the case iff the dot product of the two is 0.
		// f = (this.at(t) - p) . (this.tangentAt(t)
		// df = this.tangentAt(t) . this.tangentAt(t) + (this.at(t) - p) . this.ddt(t)
		//    = this.tangentAt(t)² + (this.at(t) - p) . this.ddt(t)
		const f = (t: number) =>
			this.at(t)
				.minus(p)
				.dot(this.tangentAt(t)) // 5th degree polynomial
		const df = (t: number) =>
			this.tangentAt(t).squared() +
			this.at(t)
				.minus(p)
				.dot(this.ddt!(t))
		//checkDerivate(f, df, tMin, tMax)

		const STEPS = 32
		if (undefined === tStart) {
			tStart = arrayFromFunction(STEPS, i => tMin + ((tMax - tMin) * i) / (STEPS - 1)).withMax(
				t => -this.at(t).distanceTo(p),
			)
		}

		return newtonIterateWithDerivative(f, tStart, 16, df)
	}

	/**
	 * So different edges on the same curve do not have different vertices, they are always generated
	 * on fixed points this.at(k * this.tIncrement), with k taking integer values
	 *
	 */
	calcSegmentPoints(aT: number, bT: number, a: V3, b: V3, reversed: boolean, includeFirst: boolean): V3[] {
		assert(this.tIncrement, 'tIncrement not defined on ' + this)
		const inc = this.tIncrement
		const result = []
		if (includeFirst) result.push(a)
		assert(reversed != aT < bT)
		if (aT < bT) {
			const start = Math.ceil((aT + NLA_PRECISION) / inc)
			const end = Math.floor((bT - NLA_PRECISION) / inc)
			for (let i = start; i <= end; i++) {
				result.push(this.at(i * inc))
			}
		} else {
			const start = Math.floor((aT - NLA_PRECISION) / inc)
			const end = Math.ceil((bT + NLA_PRECISION) / inc)
			for (let i = start; i >= end; i--) {
				result.push(this.at(i * inc))
			}
		}
		result.push(b)
		return result
	}

	calcSegmentTs(aT: number, bT: number, reversed: boolean, includeFirst: boolean): number[] {
		assert(this.tIncrement, 'tIncrement not defined on ' + this)
		const inc = this.tIncrement
		const result = []
		if (includeFirst) result.push(aT)
		assert(reversed != aT < bT)
		if (aT < bT) {
			const start = Math.ceil((aT + NLA_PRECISION) / inc)
			const end = Math.floor((bT - NLA_PRECISION) / inc)
			for (let i = start; i <= end; i++) {
				result.push(i * inc)
			}
		} else {
			const start = Math.floor((aT - NLA_PRECISION) / inc)
			const end = Math.ceil((bT + NLA_PRECISION) / inc)
			for (let i = start; i >= end; i--) {
				result.push(i * inc)
			}
		}
		result.push(bT)
		return result
	}

	/**
	 *
	 * @param p
	 * @param tStart Defines interval with tEnd in which a start value for t will be searched.
	 * Result is not necessarily in this interval.
	 * @param tEnd
	 */
	distanceToPoint(p: V3, tStart?: number, tEnd?: number) {
		const closestT = this.closestTToPoint(p, tStart, tEnd)
		return this.at(closestT).distanceTo(p)
	}

	asSegmentDistanceToPoint(p: V3, tStart: number, tEnd: number) {
		let t = this.closestTToPoint(p, tStart, tEnd)
		t = clamp(t, tStart, tEnd)
		return this.at(t).distanceTo(p)
	}

	/**
	 * Behavior when curves are colinear: self intersections
	 */
	isInfosWithCurve(curve: Curve): ISInfo[] {
		if (insideIsInfosWithCurve) {
			return Curve.ispsRecursive(this, this.tMin, this.tMax, curve, curve.tMin, curve.tMax)
		} else {
			try {
				insideIsInfosWithCurve = true
				const infos = curve.isInfosWithCurve(this)
				return infos.map(info => {
					assert(info)
					const { tThis, tOther, p } = info
					return { tOther: tThis, tThis: tOther, p }
				})
			} finally {
				insideIsInfosWithCurve = false
			}
		}
	}

	abstract transform4(m4: M4): Curve

	/**
	 * Curve point at parameter t.
	 */
	abstract at(t: number): V3

	/**
	 * Tangent of curve at parameter t. This is also the first derivative of {@see at}
	 */
	abstract tangentAt(t: number): V3

	abstract containsPoint(p: V3): boolean

	abstract isInfosWithLine(
		anchorWC: V3,
		dirWC: V3,
		tMin?: number,
		tMax?: number,
		lineMin?: number,
		lineMax?: number,
	): ISInfo[]

	abstract transform(m4: M4, desc?: string): this

	isTsWithSurface(surface: Surface): number[] {
		if (surface instanceof PlaneSurface) {
			return this.isTsWithPlane(surface.plane)
		}
		if (surface instanceof ProjectedCurveSurface) {
			const projPlane = new P3(surface.dir.unit(), 0)
			const projThis = this.project(projPlane)
			const projEllipse = surface.baseCurve.project(projPlane)
			return projEllipse.isInfosWithCurve(projThis).map(info => info.tOther)
		}
		if (surface instanceof EllipsoidSurface) {
			const thisOC = this.transform(surface.matrixInverse)
			if (!thisOC.getAABB().touchesAABBfuzzy(new AABB(V3.XYZ.negated(), V3.XYZ))) {
				return []
			}
			const f = (t: number) => thisOC.at(t).length() - 1
			const df = (t: number) =>
				thisOC
					.at(t)
					.unit()
					.dot(thisOC.tangentAt(t))

			const stepSize = 1 / (1 << 11)
			const result: number[] = []
			for (let startT = this.tMin; startT <= this.tMax; startT += stepSize) {
				const dt = stepSize * thisOC.tangentAt(startT).length()
				if (abs(f(startT)) <= dt) {
					//const t = newtonIterate1d(f, startT, 16)
					let t = newtonIterateWithDerivative(f, startT, 16, df)
					if (!eq0(f(t)) || eq0(df(t))) {
						t = newtonIterate1d(df, startT, 16)
						//if (f(a) * f(b) < 0) {
						//    t = bisect(f, a, b, 16)
						//} else if (df(a) * df(b) < 0) {
						//    t = bisect(df, a, b, 16)
						//}
					}
					if (eq0(f(t)) && !result.some(r => eq(r, t))) {
						result.push(t)
					}
				}
			}
			return result.filter(t => surface.containsPoint(this.at(t)))
		}
		throw new Error()
	}

	abstract isTsWithPlane(planeWC: P3): number[]

	arcLength(startT: number, endT: number, steps: int = 1): number {
		assert(startT < endT, 'startT < endT')
		return glqInSteps(t => this.tangentAt(t).length(), startT, endT, steps)
	}

	/**
	 * iff for any t, this.at(t) == curve.at(t)
	 */
	abstract likeCurve(curve: Curve): boolean

	equals(obj: any): boolean {
		if (this === obj) return true
		return (
			hasConstructor(obj, this.constructor) &&
			this.getConstructorParameters().equals(obj.getConstructorParameters())
		)
	}

	hashCode(): int {
		return this.getConstructorParameters().hashCode()
	}

	/**
	 * Return whether the curves occupy the same points in space. They do
	 * not necessarily need to share the same parameter values.
	 *
	 *
	 * iff for every t, there is an s so that this.at(t) == curve.at(s)
	 * and for every s, there is a t so that curve.at(s) == this.a(t)
	 */
	abstract isColinearTo(curve: Curve): boolean

	getAABB(tMin = this.tMin, tMax = this.tMax): AABB {
		tMin = isFinite(tMin) ? tMin : this.tMin
		tMax = isFinite(tMax) ? tMax : this.tMax
		const tMinAt = this.at(tMin),
			tMaxAt = this.at(tMax)
		const roots = this.roots()
		const mins = new Array(3),
			maxs = new Array(3)
		for (let dim = 0; dim < 3; dim++) {
			const tRoots = roots[dim]
			mins[dim] = Math.min(tMinAt.e(dim), tMaxAt.e(dim))
			maxs[dim] = Math.max(tMinAt.e(dim), tMaxAt.e(dim))
			for (const tRoot of tRoots) {
				if (tMin < tRoot && tRoot < tMax) {
					mins[dim] = Math.min(mins[dim], this.at(tRoot).e(dim))
					maxs[dim] = Math.max(maxs[dim], this.at(tRoot).e(dim))
				}
			}
		}
		return new AABB(V3.fromArray(mins), V3.fromArray(maxs))
	}

	/**
	 * Calculates, for each dimension individually, the parameters at which the curve has
	 * local extrema. This is mainly used for calculating the AABB of the curve.
	 * @returns {number[][]}
	 */
	abstract roots(): number[][]

	reversed(): Curve {
		throw new Error()
	}

	clipPlane(plane: P3): Curve[] {
		const ists = this.isTsWithPlane(plane).filter(ist => this.tMin <= ist && ist <= this.tMax)
		return getIntervals(ists, this.tMin, this.tMax).mapFilter(([a, b]) => {
			const midT = (a + b) / 2
			return !eq(a, b) && plane.distanceToPointSigned(this.at(midT)) < 0 && this.withBounds(a, b)
		})
	}
}

function mkcurves(
	implicitCurve: MathFunctionR2R,
	sStart: number,
	tStart: number,
	stepSize: number,
	bounds: AABB2,
	validUV: R2<boolean>,
): { points: V3[]; tangents: V3[] }[] {
	const start = V(sStart, tStart)
	assert(stepSize > 0)
	// checkDerivate(s => implicitCurve(s, 0), s => didu(s, 0), -1, 1, 0)
	// checkDerivate(t => implicitCurve(0, t), t => didv(0, t), -1, 1, 0)
	const { points, tangents } = followAlgorithm2d(implicitCurve, start, stepSize, bounds, validUV)
	if (points.length > 4 && points[0].distanceTo(points.last) <= abs(stepSize)) {
		// this is a loop: split it
		for (let i = 0; i < points.length - 1; i++) {
			assert(!points[i].equals(points[i + 1]))
		}
		const half = floor(points.length / 2)
		const points1 = points.slice(0, half),
			points2 = points.slice(half - 1, points.length)
		const tangents1 = tangents.slice(0, half),
			tangents2 = tangents.slice(half - 1, tangents.length)
		//tangents2[tangents2.length - 1] = tangents1[0]
		//points2[tangents2.length - 1] = points1[0]
		for (let i = 0; i < points1.length - 1; i++) {
			assert(!points1[i].equals(points1[i + 1]))
		}
		for (let i = 0; i < points2.length - 1; i++) {
			assert(!points2[i].equals(points2[i + 1]))
		}
		return [{ points: points1, tangents: tangents1 }, { points: points2, tangents: tangents2 }]
	} else {
		// not a loop: check in the other direction
		const { points: reversePoints, tangents: reverseTangents } = followAlgorithm2d(
			implicitCurve,
			start,
			-stepSize,
			bounds,
			validUV,
		)
		const result = followAlgorithm2d(
			implicitCurve,
			reversePoints.last,
			stepSize,
			bounds,
			validUV,
			undefined,
			reverseTangents.last.negated(),
		)
		assert(result.points.length > 2)
		return [result]
	}
}

export function breakDownPPCurves(
	ps1: ParametricSurface,
	ps2: ParametricSurface,
	uStep: number,
	vStep: number,
	stepSize: number,
): Curve[] {
	const { uMin, uMax, vMin, vMax } = ps1
	const bounds = uvInAABB2.bind(undefined, ps1)
	const bounds2 = uvInAABB2.bind(undefined, ps2)
	const deltaU = uMax - uMin,
		deltaV = vMax - vMin
	const sRes = ceil(deltaU / uStep),
		tRes = ceil(deltaV / vStep)
	const grid = new Array(sRes * tRes).fill(0)
	//const printGrid = () => console.log(arrayFromFunction(tRes, i => grid.slice(sRes * i, sRes * (i + 1)).map(v => v ? 'X' : '_').join('')).join('\n'))
	const at = (i: int, j: int) => grid[j * sRes + i]
	const set = (i: int, j: int) => 0 <= i && i < sRes && 0 <= j && j < tRes && (grid[j * sRes + i] = 1)
	const result: { points: V3[]; tangents: V3[]; st1s: V3[] }[] = []
	const logTable = []
	for (let i = 0; i < sRes; i++) {
		search: for (let j = 0; j < tRes; j++) {
			if (at(i, j)) continue
			set(i, j)
			const startU = uMin + (i + 0.5) * uStep,
				startV = vMin + (j + 0.5) * vStep
			// assume point is valid, currently (TODO)
			const curvePointPPResult = curvePointPP(ps1, ps2, ps1.pUV(startU, startV))
			if (undefined === curvePointPPResult) {
				continue search
			}
			const {
				p: startP,
				st1: { x: u, y: v },
				st2: { x: u2, y: v2 },
			} = curvePointPPResult
			const li = floor((u - uMin) / uStep),
				lj = floor((v - vMin) / vStep)
			logTable.push({
				i,
				j,
				li,
				lj,
				startU,
				startV,
				u,
				v,
				'bounds(u, v)': bounds(u, v),
			})
			if (!(i == li && j == lj) && at(li, lj)) {
				continue search
			}
			set(li, lj)
			// u, v are now good starting coordinates to use follow algorithm
			if (bounds(u, v) && bounds2(u2, v2)) {
				console.log(V(u, v).sce)
				const subResult = mkPPCurves(ps1, ps2, startP, stepSize, bounds, bounds2)
				for (const curveData of subResult) {
					assert(curveData.st1s.length > 2)
					for (const { x, y } of curveData.st1s) {
						const lif = (x - uMin) / uStep,
							ljf = (y - vMin) / vStep
						set((lif - 0.5) | 0, (ljf - 0.5) | 0)
						set((lif - 0.5) | 0, (ljf + 0.5) | 0)
						set((lif + 0.5) | 0, (ljf - 0.5) | 0)
						set((lif + 0.5) | 0, (ljf + 0.5) | 0)
					}
				}
				//printGrid()
				result.push(...subResult)
			}
		}
	}
	console.table(logTable)
	for (const { points } of result) {
		for (let i = 0; i < points.length - 1; i++) {
			assert(!points[i].equals(points[i + 1]))
		}
	}
	return result.map(({ points, tangents, st1s }) => {
		return new PPCurve(points, tangents, ps1, ps2, st1s, undefined, stepSize, 1)
	})
}

function mkPPCurves(
	ps1: ParametricSurface,
	ps2: ParametricSurface,
	startPoint: V3,
	stepSize: number,
	bounds1: (u: number, v: number) => boolean,
	bounds2: (u: number, v: number) => boolean,
): { points: V3[]; tangents: V3[]; st1s: V3[] }[] {
	// checkDerivate(s => implicitCurve(s, 0), s => didu(s, 0), -1, 1, 0)
	// checkDerivate(t => implicitCurve(0, t), t => didv(0, t), -1, 1, 0)
	const { points, tangents, st1s } = followAlgorithmPP(ps1, ps2, startPoint, stepSize, bounds1, bounds2)
	if (points[0].distanceTo(points.last) < stepSize && points.length > 2) {
		// this is a loop: split it
		for (let i = 0; i < points.length - 1; i++) {
			assert(!points[i].equals(points[i + 1]))
		}
		const half = floor(points.length / 2)
		const points1 = points.slice(0, half),
			points2 = points.slice(half - 1, points.length)
		const tangents1 = tangents.slice(0, half),
			tangents2 = tangents.slice(half - 1, tangents.length)
		const st1s1 = st1s.slice(0, half),
			st1s2 = st1s.slice(half - 1, tangents.length)
		tangents2[tangents2.length - 1] = tangents1[0]
		points2[tangents2.length - 1] = points1[0]
		st1s2[tangents2.length - 1] = st1s1[0]
		for (let i = 0; i < points1.length - 1; i++) {
			assert(!points1[i].equals(points1[i + 1]))
		}
		for (let i = 0; i < points2.length - 1; i++) {
			assert(!points2[i].equals(points2[i + 1]))
		}
		return [
			{ points: points1, tangents: tangents1, st1s: st1s1 },
			{ points: points2, tangents: tangents2, st1s: st1s2 },
		]
	} else {
		// not a loop: check in the other direction
		const { points: reversePoints } = followAlgorithmPP(ps1, ps2, startPoint, -stepSize, bounds1, bounds2)
		const result = followAlgorithmPP(ps1, ps2, reversePoints.last, stepSize, bounds1, bounds2)
		assert(result.points.length > 2)
		return [result]
	}
}

export type R2_R = (u: number, v: number) => number
export type R2<R> = (u: number, v: number) => R

export interface AABB2 {
	uMin: number
	uMax: number
	vMin: number
	vMax: number
}

export function AABB2(uMin: number, uMax: number, vMin: number, vMax: number): AABB2 {
	return { uMin, uMax, vMin, vMax }
}

export function uvInAABB2(aabb2: AABB2, u: number, v: number) {
	return aabb2.uMin <= u && u <= aabb2.uMax && aabb2.vMin <= v && v <= aabb2.vMax
}

export function curvePoint(implicitCurve: R2_R, startPoint: V3, didu: R2_R, didv: R2_R) {
	let p = startPoint
	for (let i = 0; i < 8; i++) {
		const fp = implicitCurve(p.x, p.y)
		const dfpdx = didu(p.x, p.y),
			dfpdy = didv(p.x, p.y)
		const scale = fp / (dfpdx * dfpdx + dfpdy * dfpdy)
		p = p.minus(new V3(scale * dfpdx, scale * dfpdy, 0))
	}
	return p
}

export function curvePointMF(mf: MathFunctionR2R, startPoint: V3, steps: int = 8, eps: number = 1 / (1 << 30)) {
	let p = startPoint
	for (let i = 0; i < steps; i++) {
		const fp = mf(p.x, p.y)
		const dfpdx = mf.x(p.x, p.y),
			dfpdy = mf.y(p.x, p.y)
		const scale = fp / (dfpdx * dfpdx + dfpdy * dfpdy)
		p = p.minus(new V3(scale * dfpdx, scale * dfpdy, 0))
		if (abs(fp) <= eps) break
	}
	return p
}
