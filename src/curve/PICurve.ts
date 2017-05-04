import fuzzyUniques = NLA.fuzzyUniques
const STEP_SIZE = 0.01
class PICurve extends Curve {
	parametricSurface: Surface
	implicitSurface: Surface
	startPoint: V3
	endPoint: V3
	isLoop: boolean
    points: V3[]
    pmPoints: V3[]
    tangents: V3[]
    pmTangents: V3[]
    dir: number // 1 | -1

	constructor(parametricSurface: Surface, implicitSurface: Surface, startPoint: V3, endPoint?: V3, dir: number = 1) {
		super(0, 1)
		assert(!startPoint.like(endPoint))
		assert(parametricSurface.parametricFunction, 'parametricSurface.parametricFunction')
		assert(implicitSurface.implicitFunction, 'implicitSurface.implicitFunction')
		this.parametricSurface = parametricSurface
		this.implicitSurface = implicitSurface
		if (!startPoint) {
			const pmPoint = curvePoint(this.implicitCurve(), V(1, 1, 0))
			this.startPoint = this.parametricSurface.parametricFunction()(pmPoint.x, pmPoint.y)
		} else {
			this.startPoint = startPoint
		}
		this.endPoint = endPoint
        this.dir = dir
		this.isLoop = false
        try {
            this.calcPoints(startPoint, endPoint)
            this.startPoint = startPoint
            this.endPoint = endPoint
        } catch (e) {
            this.calcPoints(this.endPoint, this.startPoint)
            this.startPoint = endPoint
            this.endPoint = startPoint
        }
        this.tMin = 0
        this.tMax = this.points.length - 1
	}

	reversed() {
	    return new PICurve(this.parametricSurface, this.implicitSurface, this.endPoint, this.startPoint, -this.dir)
    }

	toString() {
	    return `new PICurve(${this.parametricSurface.sce}, ${this.implicitSurface.sce}, ${this.startPoint.sce}, ${this.endPoint.sce}, ${this.dir})`
    }

	implicitCurve() {
		const pF = this.parametricSurface.parametricFunction()
		const iF = this.implicitSurface.implicitFunction()
		return function (s, t) {
			return iF(pF(s, t))
		}
	}

	isColinearTo(curve: Curve) {
		if (curve instanceof PICurve) {
			if (this.equals(curve)) {
				return true
			}
			if (this.parametricSurface.isCoplanarTo(curve.parametricSurface) && this.implicitSurface.isCoplanarTo(curve.implicitSurface)) {

			}
			return false
			assertNever()
		} else {
			return false
		}
	}

	containsPoint(p: V3): boolean {
		assertVectors(p)
		return !isNaN(this.pointT(p))
	}

	equals(obj: any): boolean {
		return this == obj ||
			Object.getPrototypeOf(obj) == PICurve.prototype
			&& this.parametricSurface.equals(obj.parametricSurface)
			&& this.implicitSurface.equals(obj.implicitSurface)
			&& this.startPoint.equals(obj.startPoint)
			&& this.endPoint.equals(obj.endPoint)
			&& this.dir === obj.dir
	}

	hashCode(): int {
		let hashCode = 0
		hashCode = hashCode * 31 + this.parametricSurface.hashCode()
		hashCode = hashCode * 31 + this.implicitSurface.hashCode()
		hashCode = hashCode * 31 + this.startPoint.hashCode()
		hashCode = hashCode * 31 + this.endPoint.hashCode()
		return hashCode | 0
	}


	//getVerticesNo0() {
    //
	//	// TODO
	//	let start, end, arr
	//	if (!this.canon) {
	//		start = Math.floor(this.aT + 1)
	//		end = ceil(this.bT)
	//		arr = sliceCyclic(this.curve.points, start, end)
	//	} else {
	//		start = Math.floor(this.bT + 1)
	//		end = ceil(this.aT)
	//		arr = sliceCyclic(this.curve.points, start, end)
	//		console.log("this.canon", !!this.canon, arr.length, start, end, this.aT)
	//		arr.reverse()
	//	}
	//	arr.push(this.b)
	//	return arr
	//}

	calcPoints(start, end) {
        if (!this.points) {
            const pF = this.parametricSurface.parametricFunction()
            const iF = this.implicitSurface.implicitFunction()
            const iBounds = this.implicitSurface.boundsFunction()
            const pBounds = this.parametricSurface.parametersValid.bind(this.parametricSurface)
            const curveFunction = (s, t) => iF(pF(s, t))
            const pTPF = this.parametricSurface.pointToParameterFunction()
            const startParams = pTPF(start), endParams = pTPF(end)
            //checkDerivate(x => iF(V(x, 0.1, 0.1)), x => this.implicitSurface.didp(V(x, 0.1, 0.1)).x, -2, 2)
            //checkDerivate(z => iF(V(0.1, 0.1, z)), z => this.implicitSurface.didp(V(0.1, 0.1, z)).z, -2, 2)
            //checkDerivate(y => iF(V(0.1, y, 0.1)), y => this.implicitSurface.didp(V(0.1, y, 0.1)).y, -2, 2)
            this.pmTangents = []
            this.pmPoints = followAlgorithm2(curveFunction, startParams, endParams, STEP_SIZE, null,
                this.pmTangents, pBounds, this.dir,
                (s, t) => this.implicitSurface.didp(pF(s, t)).dot(this.parametricSurface.dpds(s)),
                (s, t) => this.implicitSurface.didp(pF(s, t)).dot(this.parametricSurface.dpdt(t)),)
            assert(this.pmPoints.last().distanceTo(endParams) < STEP_SIZE)
            this.pmPoints[this.pmPoints.length - 1] = endParams
            this.isLoop = this.pmPoints[0].distanceTo(this.pmPoints[this.pmPoints.length - 1]) < STEP_SIZE * 1.1
            //if (!this.isLoop) {
            //	// the curve starting at curveStartPoint is not closed, so we need to find curve points in the other
            //	// direction until out of bounds
            //	const pmTangent0 = this.pmTangentEndPoints[0].minus(this.pmPoints[0])
            //	let pmTangentEndPoints2 = []
            //	const pmPoints2 = followAlgorithm(curveFunction, startParams, startParams, STEP_SIZE, pmTangent0.negated(),
            //		pmTangentEndPoints2, (s, t) => iBounds(pF(s, t)))
            //	pmTangentEndPoints2 = pmTangentEndPoints2.map((ep, i) => pmPoints2[i].times(2).minus(ep))
            //	this.startT = pmPoints2.length
            //	pmPoints2.reverse()
            //	pmPoints2.pop()
            //	this.pmPoints = pmPoints2.concat(this.pmPoints)
            //	pmTangentEndPoints2.reverse()
            //	pmTangentEndPoints2.pop()
            //	this.pmTangentEndPoints = pmTangentEndPoints2.concat(this.pmTangentEndPoints)
            //}
            this.points = this.pmPoints.map(({x: d, y: z}) => pF(d, z))
            this.tangents = this.points.map((p, i) => {
                const ds = this.parametricSurface.dpds(this.pmPoints[i].x)
                const dt = this.parametricSurface.dpdt(this.pmPoints[i].y)
                ds.times(this.pmTangents[i].x).plus(dt.times(this.pmTangents[i].y))
                return this.parametricSurface.normalAt(p).cross(this.implicitSurface.normalAt(p))
                    .toLength(this.dir * ds.times(this.pmTangents[i].x).plus(dt.times(this.pmTangents[i].y)).length())
            })
            //this.tangents = NLA.arrayFromFunction(this.points.length, i => {
            //    const ds = this.parametricSurface.dpds(this.pmPoints[i].x)
            //    const dt = this.parametricSurface.dpdt(this.pmPoints[i].y)
            //    return ds.times(this.pmTangents[i].x).plus(dt.times(this.pmTangents[i].y))
            //})

            for (let i = 0; i < this.points.length; i++) {
                assert(this.parametricSurface.containsPoint(this.points[i]))
                assert(this.implicitSurface.containsPoint(this.points[i]))
                const t = this.tangents[i], p = this.points[i]
                const t2 = this.parametricSurface.normalAt(p).cross(
                    this.implicitSurface.normalAt(p)).times(this.dir)
                assert(t.isParallelTo(t2))
            }
        }
	}

	tangentP(point: V3): V3 {
		assertVectors(point)
		assert(this.containsPoint(point), 'this.containsPoint(point)' + this.containsPoint(point))
		this.calcPoints(point)
		const t = this.pointT(point)
		return this.tangentAt(t)
	}

	tangentAt(t: number): V3 {
		return V3.lerp(this.tangents[floor(t)], this.tangents[ceil(t)], t % 1)
	}

	at(t: number): V3 {
	    assert(!isNaN(t))
	    if (t % 1 == 0) return this.points[t]
        const startParams = V3.lerp(this.pmPoints[floor(t)], this.pmPoints[ceil(t)], t % 1)
        return this.closestPointToParams(startParams)
    }

    closestPointToParams(startParams: V3): V3 {
        const pointParams = curvePoint(this.implicitCurve(), startParams)
        return this.parametricSurface.parametricFunction()(pointParams.x, pointParams.y)
    }

    isTsWithSurface(surface: Surface): number[] {
        if (surface instanceof PlaneSurface) {
            return this.isTsWithPlane(surface.plane)
        } else if (surface instanceof EllipsoidSurface || surface instanceof SemiEllipsoidSurface) {
	        const ps = this.parametricSurface, is = this.implicitSurface
	        if (ps instanceof ProjectedCurveSurface && is instanceof SemiEllipsoidSurface) {
		        const iscs = is.isCurvesWithSurface(surface)
		        const points = iscs.map(isc => isc.isTsWithSurface(ps).map(t => isc.at(t))).concatenated()
		        const ts = fuzzyUniques(points.map(p => this.pointT(p)))
		        return ts.filter(t => !isNaN(t))
	        }

        }
        assertNever()
    }

    ddt() {
		return V3.O
    }

    isTsWithPlane(plane: P3): number[] {
        assertInst(P3, plane)
        const ps = this.parametricSurface, is = this.implicitSurface
        if (ps instanceof ProjectedCurveSurface && is instanceof SemiEllipsoidSurface) {
            const pscs = ps.isCurvesWithPlane(plane)
            const iscs = is.isCurvesWithPlane(plane)
            const infos = iscs.map(isc => pscs.map(psc => isc.isInfosWithCurve(psc)).concatenated()).concatenated()
            const ts = fuzzyUniques(infos.map(info => this.pointT(info.p)))
            return ts.filter(t => !isNaN(t))
        }
        assertNever()
    }






	pointT(p: V3): number {
		assertVectors(p)
		if (!this.parametricSurface.containsPoint(p) || !this.implicitSurface.containsPoint(p)) {
		    return NaN
        }
		const pmPoint = this.parametricSurface.pointToParameterFunction()(p)
		let ps = this.points, pmps = this.pmPoints, t = 0, prevDistance, pmDistance = pmPoint.distanceTo(pmps[0])
		while (pmDistance > STEP_SIZE && t < ps.length - 1) { // TODO -1?
			//console.log(t, pmps[t].$, pmDistance)
			t = min(pmps.length - 1, t + max(1, Math.round(pmDistance / STEP_SIZE / 2)))
			pmDistance = pmPoint.distanceTo(pmps[t])
		}
        if (t < this.pmPoints.length - 1 && pmDistance > pmPoint.distanceTo(pmps[t + 1])) {
            t++
        }
		if (pmDistance > STEP_SIZE * 1.1) {
			// p is not on this curve
			return NaN
		}
		if (ps[t].like(p)) return t
        const startT = t + V3.inverseLerp(ps[t], ps[t + 1], p)
        return newtonIterate1d(t => this.at(t).distanceTo(p), startT, 2)
	}

	transform(m4: M4): PICurve {
	    // TODO: pass transformed points?
	    return new PICurve(
	        this.parametricSurface.transform(m4),
            this.implicitSurface.transform(m4),
            m4.transformPoint(this.startPoint),
            m4.transformPoint(this.endPoint),
            this.dir * (m4.isMirroring() ? -1 : 1))
    }

    roots(): number[][] {
		const allTs = NLA.arrayRange(0, this.points.length)
		return [allTs, allTs, allTs]
    }
}
PICurve.prototype.tIncrement = 1