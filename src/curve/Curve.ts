type ISInfo = {tThis: number, tOther: number, p: V3}
abstract class Curve extends Transformable implements Equalable {
	tIncrement: number
	hlol: number

	constructor(readonly tMin: number, readonly tMax: number) {
		super()
		assertNumbers(tMin, tMax)
		assert(!isNaN(tMin))
		assert(!isNaN(tMax))
		assert(tMin < tMax)
	}

	toSource(rounder: (x: number) => number = x => x): string {
		return callsce.call(undefined, 'new ' + this.constructor.name, ...this.getConstructorParameters())
	}

	abstract getConstructorParameters(): any[]

    /**
     * Returns curve parameter t for point p on curve.
     */
    abstract pointT(p: V3, hint?): number

	/**
	 * Returns the point on the line that is closest to the given point.
	 */
	closestPointToPoint(p: V3): V3 {
		return this.at(this.closestTToPoint(p))
	}

	isValidT(t): boolean {
		return le(this.tMin, t) && le(t, this.tMax)
	}

	diff(t: number, eps: number): V3 {
		return this.at(t).to(this.at(t + eps))
	}


	closestTToPoint(p: V3, tStart?: number): number {
		// this.at(t) has minimal distance to p when this.tangentAt(t) is perpendicular to
		// the vector between this.at(t) and p. This is the case iff the dot product of the two is 0.
		// f = (this.at(t) - p) . (this.tangentAt(t)
		// df = this.tangentAt(t) . this.tangentAt(t) + (this.at(t) - p) . this.ddt(t)
		//    = this.tangentAt(t)² + (this.at(t) - p) . this.ddt(t)
		const f = t => this.at(t).minus(p).dot(this.tangentAt(t)) // 5th degree polynomial
		const df = t => this.tangentAt(t).squared() + (this.at(t).minus(p).dot(this.ddt(t)))

		const STEPS = 32
		const startT = undefined !== tStart ? tStart :
			arrayFromFunction(STEPS, i => this.tMin + (this.tMax - this.tMin) * i / STEPS)
				.withMax(t => -this.at(t).distanceTo(p))

		return newtonIterateWithDerivative(f, startT, 16, df)
	}

	/**
	 * So different edges on the same curve do not have different vertices, they are always generated
	 * on fixed points this.at(k * this.tIncrement), with k taking integer values
	 *
	 */
	calcSegmentPoints(aT: number, bT: number, a: V3, b: V3, reversed: boolean, includeFirst: boolean): V3[] {
		assert(this.tIncrement, 'tIncrement not defined on ' + this)
		const inc = this.tIncrement
		const points = []
		if (includeFirst) points.push(a)
		assert(reversed != aT < bT)
		if (aT < bT) {
			const start = Math.ceil((aT + NLA_PRECISION) / inc)
			const end = Math.floor((bT - NLA_PRECISION) / inc)
			for (let i = start; i <= end; i++) {
				points.push(this.at(i * inc))
			}
		} else {
			const start = Math.floor((aT - NLA_PRECISION) / inc)
			const end = Math.ceil((bT + NLA_PRECISION) / inc)
			for (let i = start; i >= end; i--) {
				points.push(this.at(i * inc))
			}
		}
		points.push(b)
		return points
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

	asSegmentDistanceToPoint(p, tStart, tEnd) {
		let t = this.closestTToPoint(p, tStart, tEnd)
		t = clamp(t, tStart, tEnd)
		return this.at(t).distanceTo(p)
	}

	/**
	 * Behavior when curves are colinear: self intersections
	 */
	isInfosWithCurve(curve: Curve): ISInfo[] {
		return Curve.ispsRecursive(this, this.tMin, this.tMax, curve, curve.tMin, curve.tMax)
	}


	static integrate(curve: Curve, startT: number, endT: number, steps: int): number {
		const step = (endT - startT) / steps
		let length = 0
		let p = curve.at(startT)
		let i = 0, t = startT + step
		for (; i < steps; i++, t += step) {
			const next = curve.at(t)
			length += p.distanceTo(next)
			p = next
		}
		return length
	}

	abstract at(t: number): V3

	abstract tangentAt(t: number): V3

	/**
	 * Derivative of tangentAt for parameter t
	 */
	abstract ddt(t: number): V3

	abstract containsPoint(p: V3): boolean

	abstract transform(m4: M4, desc?: string): Curve

	abstract isTsWithSurface(surface: Surface): number[]

	abstract isTsWithPlane(plane: P3): number[]

	arcLength(startT: number, endT: number, steps?: int): number {
		assert(startT < endT, 'startT < endT')
		return glqInSteps(t => this.tangentAt(t).length(), startT, endT, steps)
	}

	/**
	 * iff for any t, this.at(t) == curve.at(t)
	 */
	abstract likeCurve(curve: Curve): boolean

	abstract equals(obj: any): boolean

    abstract hashCode(): int

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
		const tMinAt = this.at(tMin), tMaxAt = this.at(tMax)
		const roots = this.roots()
		const mins = new Array(3), maxs = new Array(3)
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

	static hlol = 0

    abstract roots(): number[][]

    static ispsRecursive(curve1: Curve, tMin: number, tMax: number, curve2: Curve, sMin: number, sMax: number): ISInfo[] {
        // the recursive function finds good approximates for the intersection points
        // curve1 function uses newton iteration to improve the result as much as possible
        const handleStartTS = (startT, startS) => {
            if (!result.some(info => eq(info.tcurve1, startT) && eq(info.tOther, startS))) {
                const f1 = (t, s) => curve1.tangentAt(t).dot(curve1.at(t).minus(curve2.at(s)))
                const f2 = (t, s) => curve2.tangentAt(s).dot(curve1.at(t).minus(curve2.at(s)))
                // f = (b1, b2, t1, t2) = b1.tangentAt(t1).dot(b1.at(t1).minus(b2.at(t2)))
                const dfdt1 = (b1, b2, t1, t2) => b1.ddt(t1).dot(b1.at(t1).minus(b2.at(t2))) + (b1.tangentAt(t1).squared())
                const dfdt2 = (b1, b2, t1, t2) => -b1.tangentAt(t1).dot(b2.tangentAt(t2))
                const ni = newtonIterate2dWithDerivatives(f1, f2, startT, startS, 16,
                    dfdt1.bind(undefined, curve1, curve2), dfdt2.bind(undefined, curve1, curve2),
                    (t, s) => -dfdt2(curve2, curve1, s, t), (t, s) => -dfdt1(curve2, curve1, s, t))
                assert(isFinite(ni.x))
                assert(isFinite(ni.y))
                if (ni == null) console.log(startT, startS, curve1.sce, curve2.sce)
                result.push({tThis: ni.x, tOther: ni.y, p: curve1.at(ni.x)})
            }
        }

        // returns whether an intersection was immediately found (i.e. without further recursion)
        function findRecursive(tMin, tMax, sMin, sMax, curve1AABB: AABB, curve2AABB: AABB, depth = 0) {
            const EPS = NLA_PRECISION
            if (curve1AABB.fuzzyTouchesAABB(curve2AABB)) {
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
                    findRecursive(tMin, tMid, sMin, sMid, curve1AABBleft, curve2AABBleft, depth + 1)
                    || findRecursive(tMin, tMid, sMid, sMax, curve1AABBleft, curve2AABBright = curve2.getAABB(sMid, sMax), depth + 1)
                    || findRecursive(tMid, tMax, sMin, sMid, curve1AABBright = curve1.getAABB(tMid, tMax), curve2AABBleft, depth + 1)
                    || findRecursive(tMid, tMax, sMid, sMax, curve1AABBright, curve2AABBright, depth + 1)
                }
            }
            return false
        }

        const result = []
        findRecursive(tMin, tMax, sMin, sMax, curve1.getAABB(tMin, tMax), curve2.getAABB(sMin, sMax))
        return fuzzyUniquesF(result, info => info.tThis)
    }

	reversed(): Curve {
		throw new Error()
	}

	static breakDownIC(implicitCurve: (s, t) => number,
	                   sMin: number, sMax: number,
	                   tMin: number, tMax: number,
	                   sStep: number, tStep: number,
	                   stepSize: number,
	                   dids?: (s, t) => number,
	                   didt?: (s, t) => number): {points: V3[], tangents: V3[]}[] {
		const EPS = 1 / (1 << 20)
		undefined == dids && (dids = (s, t) => (implicitCurve(s + EPS, t) - implicitCurve(s, t)) / EPS)
		undefined == didt && (didt = (s, t) => (implicitCurve(s, t + EPS) - implicitCurve(s, t)) / EPS)

		const bounds = (s, t) => sMin <= s && s <= sMax && tMin <= t && t <= tMax
		const deltaS = sMax - sMin, deltaT = tMax - tMin
		const sRes = ceil(deltaS / sStep), tRes = ceil(deltaT / tStep)
		const grid = new Array(sRes * tRes)
		const at = (i, j) => grid[j * sRes + i]
		const set = (i, j) => 0 <= i && i < sRes && 0 <= j && j < tRes && (grid[j * sRes + i] = 1)
		const result = []
		for (let i = 0; i < sRes; i++) {
			search: for (let j = 0; j < tRes; j++) {
				if (at(i, j)) continue
				set(i, j)
				let s = sMin + i * sStep, t = tMin + j * tStep

				// basically curvePoint
				for (let i = 0; i < 8; i++) {
					const fp = implicitCurve(s, t)
					const dfpdx = dids(s, t), dfpdy = didt(s, t)
					const scale = fp / (dfpdx * dfpdx + dfpdy * dfpdy)
					s -= scale * dfpdx
					t -= scale * dfpdy
				}
				const li = floor((s - sMin) / sStep), lj = floor((t - tMin) / tStep)
				if (at(li, lj)) {
					continue search
				}
				if (0 <= li && li < sRes && 0 <= lj && lj < tRes) {
					set(li, lj)
				}
				// s, t are now good starting coordinates to use follow algo
				if (bounds(s, t)
					&& eq02(implicitCurve(s, t), 0.01)
				) {
					console.log(V(s, t).sce)
					const subresult = mkcurves(implicitCurve, s, t, stepSize, dids, didt, bounds)
					for (const curvedata of subresult) {
						for (const {x, y} of curvedata.points) {
							const lif = (x - sMin) / sStep, ljf = (y - tMin) / tStep
							set((lif - 0.5) | 0, (ljf - 0.5) | 0)
							set((lif - 0.5) | 0, (ljf + 0.5) | 0)
							set((lif + 0.5) | 0, (ljf - 0.5) | 0)
							set((lif + 0.5) | 0, (ljf + 0.5) | 0)
						}
					}
					result.pushAll(subresult)
				}

			}
		}
		return result
	}

	static test2() {
		const ic = (x, y) => sin(x+y)-cos(x*y)+1
		const dids = (x, y) => y * sin(x * y) + cos(x + y)
		const didt = (x, y) => x * sin(x * y) + cos(x + y)
		const ic2 = (x, y) => (3 * x ** 2 - y ** 2) ** 2 * y ** 2 - (x ** 2 + y ** 2) ** 4
		const di2ds = (x, y) => 4* x* (9* x**2* y**2 - 3* y**4 - 2* (x**2 + y**2)**3)
		const di2dt = (x, y) => 2 * y * (-4 * (x ** 2 + y ** 2) ** 3 + (3 * x ** 2 - y ** 2) ** 2 + 2 * y ** 2 * (y ** 2 - 3 * x ** 2))
		const start = V(-3.6339970071165784, 3.5625834844534974, 0) // curvePoint(ic, V(-4, 4))
		assert(eq02(ic(start.x, start.y), 0.1))
		const bounds = (s, t) => -5 <= s && s <= 5 && -5 <= t && t <= 5
		//const curves =  Curve.breakDownIC(ic, -5, 5, -5, 5, 0.1, 0.1, 0.05, dids, didt)
		const curves =  Curve.breakDownIC(ic2, -5, 5, -5, 5, 0.1, 0.1, 0.02, di2ds, di2dt)
		//const curves =  Curve.breakDownIC(cassini(1, 1.02), -5, 5, -5, 5, 0.1, 0.1, 0.02)
		//const curves = mkcurves(ic, start.x, start.y, 0.05, dids, didt, bounds)
			.map(({points, tangents}, i) => {
				const curve = new ImplicitCurve(ic, points, tangents)
				return Edge.forCurveAndTs(curve.translate(5, 0, 0.1 * i))
			})
		//checkDerivate(s => ic(s, 0), s => dids(s, 0), -5, 5, 0)
		//checkDerivate(t => ic(0, t), t => dids(0, t), -5, 5, 0)
		console.log(curves.length)
		return curves

	}
	static test() {
		const ses = SemiEllipsoidSurface.UNIT
		const cs = ConicSurface.UNIT.scale(0.05,0.2).translate(-0.1,0,-2).rotateX(5*DEG)
		const pf = cs.parametricFunction(), icc = ses.implicitFunction()
		const ic = (x, y) => icc(pf(x, y))
		//const start = V(-3.6339970071165784, 3.5625834844534974, 0) // curvePoint(ic, V(-4, 4))
		//assert(eq02(ic(start.x, start.y), 0.1))
		const sMin = 0, sMax = PI, tMin = 0, tMax = 5
		const bounds = (s, t) => sMin <= s && s <= sMax && tMin <= t && t <= tMax
		//const curves =  Curve.breakDownIC(ic, -5, 5, -5, 5, 0.1, 0.1, 0.05, dids, didt)
		const curves =  Curve.breakDownIC(ic, sMin, sMax, tMin, tMax, 0.1, 0.1, 0.02)
		//const curves =  Curve.breakDownIC(cassini(1, 1.02), -5, 5, -5, 5, 0.1, 0.1, 0.02)
		//const curves = mkcurves(ic, start.x, start.y, 0.05, dids, didt, bounds)
			.map(({points: pmPoints, tangents: pmTangents}, i) => {
				console.log(icc)
				const points = pmPoints.map(({x, y}) => pf(x, y))
				const curve = new ImplicitCurve(ic, points, pmTangents)
				return Edge.forCurveAndTs(curve)
				//return Edge.forCurveAndTs(curve.translate(5, 0, 0.1 * i))
			})
		//checkDerivate(s => ic(s, 0), s => dids(s, 0), -5, 5, 0)
		//checkDerivate(t => ic(0, t), t => dids(0, t), -5, 5, 0)
		console.log(curves.length)
		return curves

	}
}

function mkcurves(implicitCurve: (s: number, t: number) => number,
                  sStart: number, tStart: number,
                  stepSize: number,
                  dids: (s, t) => number,
                  didt: (s, t) => number,
                  bounds: (s, t) => boolean): {points: V3[], tangents: V3[]}[] {
	const start = V(sStart, tStart)
	checkDerivate(s => implicitCurve(s, 0), s => dids(s, 0), -1, 1, 0)
	checkDerivate(t => implicitCurve(0, t), t => didt(0, t), -1, 1, 0)
	const {points, tangents} = followAlgorithm2d(implicitCurve, start, stepSize, dids, didt, bounds)
	if (points[0].distanceTo(points.last()) < stepSize) {
		// this is a loop: split it
		const half = floor(points.length / 2)
		const points1 = points.slice(0, half), points2 = points.slice(half - 1, points.length)
		const tangents1 = tangents.slice(0, half), tangents2 = tangents.slice(half - 1, tangents.length)
		tangents2[tangents2.length - 1] = tangents1[0]
		points2[tangents2.length - 1] = points1[0]
		return [{points: points1, tangents: tangents1}, {points: points2, tangents: tangents2}]
	} else {
		// not a loop: check in the other direction
		const {points: points2, tangents: tangents2} = followAlgorithm2d(implicitCurve, start, -stepSize, dids, didt, bounds)
		const allPoints = points2.reverse().concat(points)
		const allTangents = tangents2.reverse().concat(tangents)
		return [{points: allPoints, tangents: allTangents}]
	}
}


function curvePoint(implicitCurve: (s: number, t: number) => number, startPoint: V3) {
	const eps = 1 / (1 << 20)
	let p = startPoint
	for (let i = 0; i < 4; i++) {
		const fp = implicitCurve(p.x, p.y)
		const dfpdx = (implicitCurve(p.x + eps, p.y) - fp) / eps,
			dfpdy = (implicitCurve(p.x, p.y + eps) - fp) / eps
		const scale = fp / (dfpdx * dfpdx + dfpdy * dfpdy)
		//console.log(p.$)
		p = p.minus(new V3(scale * dfpdx, scale * dfpdy, 0))
	}
	return p
}