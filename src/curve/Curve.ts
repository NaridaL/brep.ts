abstract class Curve extends Transformable {
	tIncrement: number
	hlol: number

	constructor(readonly tMin: number, readonly tMax: number) {
		super()
		assertNumbers(tMin, tMax)
		assert(!isNaN(tMin))
		assert(!isNaN(tMax))
		assert(tMin < tMax)
	}

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

	isValidT(t) {
		return le(this.tMin, t) && le(t, this.tMax)
	}


	/**
	 *
	 * @param p
	 * @param tMin Defines interval with tMax in which a start value for t will be searched.
	 * Result is not necessarily in this interval.
	 * @param tMax
	 * @param tStart
	 */
	closestTToPoint(p: V3, tStart?: number): number {


		// this.at(t) has minimal distance to p when this.tangentAt(t) is perpendicular to
		// the vector between this.at(t) and p. This is the case iff the dot product of the two is 0.
		// f = (this.at(t) - p) . (this.tangentAt(t)
		// df = this.tangentAt(t) . this.tangentAt(t) + (this.at(t) - p) . this.ddt(t)
		//    = this.tangentAt(t)² + (this.at(t) - p) . this.ddt(t)
		const f = t => this.at(t).minus(p).dot(this.tangentAt(t)) // 5th degree polynomial
		const df = t => this.tangentAt(t).squared() + (this.at(t).minus(p).dot(this.ddt(t)))

		const STEPS = 32
		const startT = 'undefined' === typeof tStart
			? NLA.arrayFromFunction(STEPS, i => this.tMin + (this.tMax - this.tMin) * i / STEPS)
				.withMax(t => -this.at(t).distanceTo(p))
			: tStart

		return newtonIterateWithDerivative(f, startT, 16, df)
	}

	/**
	 * So different edges on the same curve do not have different vertices, they are always generated
	 * on fixed points this.at(k * this.tIncrement), with k taking integer values
	 *
	 */
	calcSegmentPoints(aT: number, bT: number, a: V3, b: V3, reversed: boolean, includeFirst: boolean): V3[] {
		assert(this.tIncrement, "tIncrement not defined on " + this)
		const split = 4 * 62, inc = this.tIncrement
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

	distanceToPoint(p: V3): number {
		return this.at(this.closestTToPoint(p)).distanceTo(p)
	}

	/**
	 * Behavior when curves are colinear: self intersections
	 */
	abstract isInfosWithCurve(curve): {tThis: number, tOther: number, p: V3}[]

	abstract at(t: number): V3

	abstract tangentAt(t: number): V3

	/**
	 * Derivative of tangentAt for parameter t
	 */
	abstract ddt(t: number): V3

	abstract containsPoint(p: V3): boolean

	abstract isTsWithSurface(surface: Surface): number[]

	abstract isTsWithPlane(plane: P3): number[]

	/**
	 * Should really be abstract, but it works for all the conic is curves, so it's here.
	 */
	debugToMesh(mesh, bufferName) {
		mesh[bufferName] || mesh.addVertexBuffer(bufferName, bufferName)
		for (let t = 0; t < Math.PI; t += 0.1) {
			const p = this.at(t)
			mesh[bufferName].push(p, p.plus(this.tangentAt(t).toLength(1)))
			mesh[bufferName].push(p, p.plus(this.normalAt(t).toLength(1)))
		}
		mesh[bufferName].push(this.center, this.center.plus(this.f1.times(1.2)))
		mesh[bufferName].push(this.center, this.center.plus(this.f2))
		mesh[bufferName].push(this.center, this.center.plus(this.normal))
	}

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
		let tMinAt = this.at(tMin), tMaxAt = this.at(tMax)
		let roots = this.roots()
		let mins = new Array(3), maxs = new Array(3)
		for (let dim = 0; dim < 3; dim++) {
			let tRoots = roots[dim]
			mins[dim] = Math.min(tMinAt.e(dim), tMaxAt.e(dim))
			maxs[dim] = Math.max(tMinAt.e(dim), tMaxAt.e(dim))
			for (let j = 0; j < tRoots.length; j++) {
				let tRoot = tRoots[j]
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

    static ispsRecursive(curve1: Curve, tMin: number, tMax: number, curve2: Curve, sMin: number, sMax: number): { tThis: number, tOther: number, p: V3 }[] {
        // the recursive function finds good approximates for the intersection points
        // curve1 function uses newton iteration to improve the result as much as possible
        const handleStartTS = (startT, startS) => {
            if (!result.some(info => NLA.eq(info.tcurve1, startT) && NLA.eq(info.tOther, startS))) {
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
        return NLA.fuzzyUniquesF(result, info => info.tThis)
    }

	reversed(): Curve {
		throw new Error()
	}
}