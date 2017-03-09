class SemiEllipseCurve extends Curve {
    center: V3
    f1: V3
    f2: V3
    normal: V3
    matrix: M4
    inverseMatrix: M4
    tMin: number
    tMax: number

    constructor(center, f1, f2, minT?, maxT?) {
        super()
        assertVectors(center, f1, f2)
        this.center = center
		this.f1 = f1
		this.f2 = f2
        undefined != minT && (this.tMin = minT)
        undefined != maxT && (this.tMax = maxT)
        assert(undefined == minT || undefined == maxT || this.tMin < this.tMax)
        assert(undefined == minT || undefined == maxT || 0 <= this.tMin && this.tMin < PI)
        assert(undefined == minT || undefined == maxT || 0 < this.tMax && this.tMax <= PI)
		this.normal = f1.cross(f2)
		if (!this.normal.isZero()) {
			this.normal = this.normal.normalized()
			this.matrix = M4.forSys(f1, f2, this.normal, center)
			this.inverseMatrix = this.matrix.inversed()
		} else {
			this.matrix = M4.forSys(f1, f2, f1.normalized(), center)
			let f1p = f1.getPerpendicular()
			this.inverseMatrix = new M4(
				1, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 1).times(M4.forSys(f1, f1p, f1.cross(f1p), center).inversed())
		}
	}

	getVolZAnd(dir1: V3, tStart: number, tEnd: number): {volume: number, centroid: V3} {
		// let p = at(t)
		// integrate area [p -> plane.projectPoint(p)] to x axis...
		// INTEGRATE[tStart, tEnd] fp(this.at(t)) dt
		const ft = t => this.center.plus(this.f1.times(Math.cos(t))).plus(this.f2.times(Math.sin(t)))
		// f(t) = c + f1 cos + f2 sin
		// p dot d1 = (cx + f1x cos + f2x sin) dx + (cy + f1y cos + f2y sin) dy + (cz + f1z cos + f2z sin) dz
		function fp(p) {
			const p0ToP = dir1.times(dir1.dot(p))
			const area = p0ToP.lengthXY() * (p.z - p0ToP.z / 2)
			return area
		}
		const f = t => fp(this.at(t)) * this.tangentAt(t).cross(this.normal).normalized().z
		return {volume: glqInSteps(f, tStart, tEnd, 4), centroid: undefined}
	}

	getAreaInDir(right: V3, up: V3, tStart: number, tEnd: number): {area: number, centroid: V3} {
		//assertf(() => tStart < tEnd)
		assertf(() => right.isPerpendicularTo(this.normal))
		assertf(() => up.isPerpendicularTo(this.normal))
		//assertf(() => SemiEllipseCurve.isValidT(tStart), tStart)
		//assertf(() => SemiEllipseCurve.isValidT(tEnd), tEnd)

		let localUp = this.inverseMatrix.transformVector(up)
		let localRight = localUp.cross(V3.Z)
		let normTStart = tStart - localRight.angleXY()
		let normTEnd = tEnd - localRight.angleXY()
		let transformedOriginY = this.inverseMatrix.getTranslation().dot(localUp.normalized())
		//console.log(localUp.str, localRight.str, normTStart, normTEnd, 'localUp.length()', localUp.length())
		//console.log('transformedOriginY', transformedOriginY)
		//assertf(() => localUp.hasLength(1), localUp.length())
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

		let restArea = -transformedOriginY * (-Math.cos(normTEnd) + Math.cos(normTStart) )
		let area = fArea(normTEnd) - fArea(normTStart) + restArea
		let cxt = (cxTimesArea(normTEnd) - cxTimesArea(normTStart) + -transformedOriginY * (-Math.cos(normTEnd) - Math.cos(normTStart)) / 2 * restArea) / area
		let cyt = (cyTimesArea(normTEnd) - cyTimesArea(normTStart) - -transformedOriginY / 2 * restArea) / area
		let factor = this.matrix.xyAreaFactor() // * localUp.length()
		//console.log('fctor', factor, 'area', area, 'resultarea', area* factor)
		assert(!NLA.eq0(factor))
		return {area: area * factor, centroid: this.matrix.transformPoint(M4.rotationZ(localRight.angleXY()).transformPoint(new V3(cxt, cyt, 0)))}

	}

	toString(f?) {
		return `new SemiEllipseCurve(${this.center}, ${this.f1}, ${this.f2}${
		    this.hasOwnProperty('tMin') || this.hasOwnProperty('tMax') ? `, ${this.tMin}, ${this.tMax}` : ''})`
	}

	isValidT(t) {
		return le(this.tMin, t) && le(t, this.tMax)
	}

	at(t) {
        assertNumbers(t)
        assert(this.isValidT(t))
		return this.center.plus(this.f1.times(Math.cos(t))).plus(this.f2.times(Math.sin(t)))
	}

	at2(xi, eta) {
		// center + f1 xi + f2 eta
		return this.center.plus(this.f1.times(xi)).plus(this.f2.times(eta))
	}

	tangentAt(t) {
		assertNumbers(t)
        assert(this.isValidT(t))
		// f2 cos(t) - f1 sin(t)
		return this.f2.times(Math.cos(t)).minus(this.f1.times(Math.sin(t)))
	}

	ddt(t) {
		assertNumbers(t)
        assert(this.isValidT(t))
		return this.f2.times(-Math.sin(t)).minus(this.f1.times(Math.cos(t)))
	}

	tangentAt2(xi, eta) {
		return this.f2.times(xi).minus(this.f1.times(eta))
	}

	/**
	 *
	 */
	isCircular() {
		return NLA.eq(this.f1.length(), this.f2.length()) && this.f1.isPerpendicularTo(this.f2)
	}


	equals(obj: any): boolean {
		return this == obj ||
			Object.getPrototypeOf(obj) == SemiEllipseCurve.prototype
			&& this.center.equals(obj.center)
			&& this.f1.equals(obj.f1)
			&& this.f2.equals(obj.f2)
	}

    hashCode(): int {
	    let hashCode = 0
        hashCode = hashCode * 31 + this.center.hashCode()
        hashCode = hashCode * 31 + this.f1.hashCode()
        hashCode = hashCode * 31 + this.f2.hashCode()
        return hashCode | 0
    }

	/**
	 * @inheritDoc
	 */
	likeCurve(curve) {
		return curve.constructor == SemiEllipseCurve
			&& this.center.like(curve.center)
			&& this.f1.like(curve.f1)
			&& this.f2.like(curve.f2)
	}

	isColinearTo(curve) {
		if (curve.constructor != SemiEllipseCurve) {
			return false
		}
		const ell = curve
		if (!this.center.like(ell.center)) {
			return false
		}
		if (this == ell) {
			return true
		}
		if (this.isCircular()) {
			return ell.isCircular() && NLA.eq(this.f1.length(), ell.f1.length()) && this.normal.isParallelTo(ell.normal)
		} else {
			let {f1: f1, f2: f2} = this.mainAxes(), {f1: c1, f2: c2} = ell.mainAxes()
			if (f1.length() > f2.length()) {[f1, f2] = [f2, f1]}
			if (c1.length() > c2.length()) {[c1, c2] = [c2, c1]}
			console.log(f1.squared(), Math.abs(f1.dot(c1)), f2.squared(), Math.abs(f2.dot(c2)))
			return NLA.eq(f1.squared(), Math.abs(f1.dot(c1)))
				&& NLA.eq(f2.squared(), Math.abs(f2.dot(c2)))
		}
	}

	normalAt(t) {
		return this.tangentAt(t).cross(this.normal)
	}

    static validPlanePoint(x: number, y: number): boolean {
        return le(0, y) && eq0(x ** 2 + y ** 2 - 1)
    }
    static unitT(p): number {
        assert(le(0, p.y))
        const angle = Math.atan2(p.y, p.x)
        return angle < 0 ? (assert(eq0(angle) || eq(PI, abs(angle))) && abs(angle)) : angle
    }

	pointT(p: V3) {
		assertVectors(p)
        assert(this.containsPoint(p))
		const pLC = this.inverseMatrix.transformPoint(p)
        const t = SemiEllipseCurve.unitT(pLC)
        assert(this.isValidT(t))
        return t
	}

	isOrthogonal(p) {
		return this.f1.isPerpendicularTo(this.f2)
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
	mainAxes(): {f1: V3, f2: V3} {
		const f1 = this.f1, f2 = this.f2, a = f1.dot(f2), b = f2.squared() - f1.squared()
		if (NLA.eq0(a)) {
			return this
		}
		const g1 = 2 * a, g2 = b + Math.sqrt(b * b + 4 * a * a)
		const {x1: xi, y1: eta} = intersectionUnitCircleLine(g1, g2, 0)
		return {f1: f1.times(xi).plus(f2.times(eta)), f2: f1.times(-eta).plus(f2.times(xi))}
	}

	eccentricity() {
		const mainAxes = this.mainAxes()
		const f1length = mainAxes.f1.length(), f2length = mainAxes.f1.length()
		const [a, b] = f1length > f2length ? [f1length, f2length] : [f2length, f1length]
		return Math.sqrt(1 - b * b / a / a)
	}

	circumference(): number {
		return this.arcLength(-Math.PI, Math.PI)
	}

	arcLength(startT: number, endT: number, steps: int = 2): number {
		assert(startT < endT, 'startT < endT')
		const f1Length = this.f1.length()
		if (eq(f1Length, this.f2.length())) {
			return f1Length * (endT - startT)
		}
		return super.arcLength(startT, endT, steps)
	}

	circumferenceApproximate(): number {
		// approximate circumference by Ramanujan
		// https://en.wikipedia.org/wiki/Ellipse#Circumference
		const {f1, f2} = this.mainAxes(), a = f1.length(), b = f2.length()
		const h = (a - b) * (a - b) / (a + b) / (a + b) // (a - b)² / (a + b)²
		return Math.PI * (a + b) * (1 + 3 * h / (10 + Math.sqrt(4 - 3 * h)))
	}

	transform(m4): this {
		return new SemiEllipseCurve(
		    m4.transformPoint(this.center),
            m4.transformVector(this.f1),
            m4.transformVector(this.f2),
            this.tMin, this.tMax) as this
    }

    rightAngled(): SemiEllipseCurve {
        const {f1, f2} = this.mainAxes()
        return new SemiEllipseCurve(this.center, f1, f2)
    }

    isTsWithSurface(surface) {
        if (surface instanceof PlaneSurface) {
            return this.isTsWithPlane(surface.plane)
        } else if (surface instanceof SemiCylinderSurface) {
            const ellipseProjected = surface.baseEllipse.transform(M4.projection(this.getPlane(), surface.dir1))
            return this.isInfosWithEllipse(ellipseProjected).map(info => info.tThis)
		} else {
			assert(false)
		}
	}

	isTsWithPlane(plane) {
		assertInst(P3, plane)
		/*
		 this: x = center + f1 * cos t + f2 * sin t  (1)
		 n DOT center + n DOT f1 cos t + n DOT f2 sin t - w = 0
		 plane:
		 n := plane.normal
		 n DOT x == plane.w           (2)
		 plane defined by f1/f2
		 x = center + f1 * xi + f2 * eta         (3)
		 intersection plane and planef1/f2:
		 insert (3) into (2):
		 n DOT center + n DOT f1 * xi + n DOT f2 * eta = plane.w | -n DOT center
		 n DOT f1 * xi + n DOT f2 * eta = plane.w - n DOT center (4)
		 points on ellipse have additional condition
		 eta * eta + xi * xi = 1 (5)
		 g1 := n DOT f1
		 g2 := n DOT f2
		 g3 := w - n DOT center
		 solve system (5)/(6)
		 g1 * eta + g2 * eta = g3 (6)
		 */
		if (plane.normal.isParallelTo(this.normal)) {
			return []
		}
		const n = plane.normal, w = plane.w,
			center = this.center, f1 = this.f1, f2 = this.f2,
			g1 = n.dot(f1), g2 = n.dot(f2), g3 = w - n.dot(center)

		const {x1: xi1, y1: eta1, x2: xi2, y2: eta2} = intersectionUnitCircleLine(g1, g2, g3)
        const result = []
        le(0, eta1) && result.push(SemiEllipseCurve.unitT(new V3(xi1, eta1, 0)))
        le(0, eta2) && result.push(SemiEllipseCurve.unitT(new V3(xi2, eta2, 0)))
        for (const t of result) {
		    assert(plane.containsPoint(this.at(t)))
        }
		return result

	}

    getPlane(): P3 {
        return P3.normalOnAnchor(this.normal, this.center)
    }

	containsPoint(p) {
		const pLC = this.inverseMatrix.transformPoint(p)
		return NLA.eq0(pLC.z) && SemiEllipseCurve.validPlanePoint(pLC.x, pLC.y)
	}

	isInfosWithEllipse(ellipse: SemiEllipseCurve): {tThis: number, tOther: number, p: V3}[] {
		if (this.normal.isParallelTo(ellipse.normal) && NLA.eq0(this.center.minus(ellipse.center).dot(ellipse.normal))) {

			// ellipses are coplanar
			const ellipseLC = ellipse.transform(this.inverseMatrix).rightAngled()

			// check if colinear
			if (ellipseLC.f1.hasLength(1) && ellipseLC.f2.hasLength(1) && ellipseLC.center.isZero()) {
				return []
			}
            return Curve.ispsRecursive(this, this.tMin, this.tMax, ellipse, ellipse.tMin, ellipse.tMax)

			//new SemiEllipseCurve(V3.ZERO, V3.X, V3.Y).debugToMesh(dMesh, 'curve4')
			console.log(ellipseLC, ellipseLC.sce)
			//ellipseLC.debugToMesh(dMesh, 'curve3')
			const angle = ellipseLC.f1.angleXY()
			console.log('angle', angle)
			const aSqr = ellipseLC.f1.squared(), bSqr = ellipseLC.f2.squared()
			const a = Math.sqrt(aSqr), b = Math.sqrt(bSqr)
			const {x: centerX, y: centerY} = ellipseLC.center
			const rotCenterX = centerX * Math.cos(-angle) + centerY * -Math.sin(-angle)
			const rotCenterY = centerX * Math.sin(-angle) + centerY * Math.cos(-angle)
			const rotCenter = V(rotCenterX, rotCenterY)
			let f = t => {
				const lex = Math.cos(t) - rotCenterX, ley = Math.sin(t) - rotCenterY
				return lex * lex / aSqr + ley * ley / bSqr - 1
			}
			//uc.debugToMesh(dMesh, 'curve4')
			const f2 = (x, y) => 200 * (x * x + y * y - 1)
			const f3 = (x, y) => 200 * ((x - rotCenterX) * (x - rotCenterX) / aSqr + (y - rotCenterY) * (y - rotCenterY) / bSqr - 1)
			const results = []
			const resetMatrix = this.matrix.times(M4.rotationZ(angle))
			for (let da = Math.PI / 4; da < Math.PI; da += Math.PI / 2) {
				const startP = SemiEllipseCurve.UNIT.at(da)
				const p = newtonIterate2d(f3, f2, startP.x, startP.y, 10)
				if (p && !results.some(r => r.like(p))) {
					results.push(p)
					drPs.push(p)
				}
			}
			const rotEl = new SemiEllipseCurve(rotCenter, V(a, 0, 0), V(0, b, 0))
			//var rotEl = ellipseLC.transform(resetMatrix)
			console.log(rotEl, rotEl.sce)
			//rotEl.debugToMesh(dMesh, 'curve2')
			return results.map(localP => ({tThis: undefined, tOther: undefined, p: resetMatrix.transformPoint(localP)}))
			/*
			 // new rel center
			 var mat = M4.forSys(ellipseLC.f1.normalized(), ellipseLC.f2.normalized(), V3.Z, ellipseLC.center).inversed()
			 console.log(mat.toString())
			 var newCenter = mat.transformPoint(V3.ZERO)
			 var x0 = newCenter.x, y0 = newCenter.y
			 var c = (1 - bSqr / aSqr) / 2/ y0, d = -x0 / y0, e = (bSqr + x0 * x0 + y0 * y0) / 2 / y0

			 var ff = x => c*c*x*x*x*x+ 2*c*d*x*x*x+ (2*c*e+d*d+bSqr/aSqr)*x*x+2*d*e*x+e*e-bSqr
			 var newx1 = newtonIterate1d(ff, x0 - 1)
			 var newx2 = newtonIterate1d(ff, x0)
			 var f1 = (x, y) => 2 * x - y
			 var f2 = (x, y) => 2 * ((x - x0) * (x - x0) + (y - y0) * (y - y0) - 1)
			 var f3 = (x, y) => 2 * (x * x / aSqr + y * y / bSqr - 1)
			 ellipseLC = ellipseLC.transform(mat)
			 for (var a = PI / 4; a < 2 * PI; a+= PI / 2) {
			 var startP = circle.at(a)
			 //drPs.push(startP)
			 var p = newtonIterate2d(f3, f2, startP.x, startP.y, 10)
			 p && drPs.push(p)
			 p && console.log(p.$, p.minus(V(x0, y0)).length())
			 }
			 circle.debugToMesh(dMesh, 'curve1')
			 ellipseLC.debugToMesh(dMesh, 'curve2')
			 */
		} else {
			return this.isTsWithPlane(P3.normalOnAnchor(ellipse.normal.normalized(), ellipse.center)).mapFilter(t => {
				const p = this.at(t)
				if (ellipse.containsPoint(p)) {
					return {tThis: t, tOther: ellipse.pointT(p), p}
				}
			})
		}
	}

	isInfosWithLine(line) {
		const anchorLC = this.inverseMatrix.transformPoint(line.anchor)
		const dirLC = this.inverseMatrix.transformVector(line.dir1)
		if (NLA.eq0(dirLC.z)) {
			// local line parallel to XY-plane
			if (NLA.eq0(anchorLC.z)) {
				// local line lies in XY-plane
				// ell: x² + y² = 1 = p²
				// line(t) = anchor + t dir
				// anchor² - 1 + 2 t dir anchor + t² dir² = 0
				const pqDiv = dirLC.squared()
				const lineTs = pqFormula(2 * dirLC.dot(anchorLC) / pqDiv, (anchorLC.squared() - 1) / pqDiv)
				return lineTs.filter(tOther => le(0, anchorLC.y + tOther * dirLC.y))
                    .map(tOther => ({
					tThis: SemiEllipseCurve.unitT(dirLC.times(tOther).plus(anchorLC)),
					tOther: tOther,
					p: line.at(tOther)}))
			}
		} else {
			// if the line intersects the XY-plane in a single point, there can be an intersection there
			// find point, then check if distance from circle = 1
			const otherTAtZ0 = anchorLC.z / dirLC.z
			const isp = dirLC.times(otherTAtZ0).plus(anchorLC)
			if (SemiEllipseCurve.validPlanePoint(isp.x, isp.y)) {
				// point lies on unit circle
				return [{
					tThis: SemiEllipseCurve.unitT(isp),
					tOther: otherTAtZ0,
					p: line.at(otherTAtZ0)}]
			}
		}
		return []
	}

	isInfosWithCurve(curve): { tThis: number, tOther: number, p: V3 }[] {
        if (curve instanceof L3) {
            return this.isInfosWithLine(curve)
        }
        if (curve instanceof BezierCurve) {
            return this.isInfosWithBezier(curve)
        }
        if (curve instanceof SemiEllipseCurve) {
            return this.isInfosWithEllipse(curve)
        }
        assert(false)
    }

	/**
	 *
	 * @param bezier
	 */
    isPointsWithBezier(bezier: BezierCurve) {
        const bezierLC = bezier.transform(this.inverseMatrix)
        if (new PlaneSurface(P3.XY).containsCurve(bezier)) {
            // up to 6 solutions possible
            let f = t => bezierLC.at(t).squaredXY() - 1
            // f is polynome degree six, no explicit solution is possble
            const possibleOtherTs = NLA.arrayFromFunction(16, i => newtonIterate1d(f, i / 15, 8))
                .filter(t => SemiEllipseCurve.validPlanePoint(bezierLC.at(t).x, bezierLC.at(t).y))
            return NLA.fuzzyUniques(possibleOtherTs).map(t => bezier.at(t))
        } else {
            return bezierLC.isTsWithPlane(P3.XY)
                .filter(t => SemiEllipseCurve.validPlanePoint(bezierLC.at(t).x, bezierLC.at(t).y))
                .map(t => bezier.at(t))
        }
    }

	isInfosWithBezier(bezier: BezierCurve): {tThis: number, tOther: number, p: V3}[] {
		const bezierLC = bezier.transform(this.inverseMatrix)
		if (new PlaneSurface(P3.XY).containsCurve(bezier)) {
			return this.isInfosWithBezier2D(bezier)
		} else {
            return bezierLC.isTsWithPlane(P3.XY).mapFilter(tOther => {
				const pLC = bezierLC.at(tOther)
				if (SemiEllipseCurve.validPlanePoint(pLC.x, pLC.y)) {
					return {tOther: tOther, p: bezier.at(tOther), tThis: SemiEllipseCurve.unitT(pLC)}
				}
			})
		}
	}

    isInfosWithBezier2D(bezier: BezierCurve, sMin?: number, sMax?: number): { tThis: number, tOther: number, p: V3 }[] {
        sMin = isFinite(sMin) ? sMin : bezier.tMin
		sMax = isFinite(sMax) ? sMax : bezier.tMax
		assertf(() => 0 < Math.PI)
		assertf(() => sMin < sMax)
        return Curve.ispsRecursive(this, this.tMin, this.tMax, bezier, sMin, sMax)
	}


	roots(): number[][] {
		// tangent(t) = f2 cos t - f1 sin t
		// solve for each dimension separately
		// tangent(eta, xi) = f2 eta - f1 xi

		return NLA.arrayFromFunction(3, dim => {
			let a = this.f2.e(dim), b = -this.f1.e(dim)
			let {x1,y1,x2,y2} = intersectionUnitCircleLine(a, b, 0)
			return [Math.atan2(y1, x1), Math.atan2(y2, x2)]
		})
	}

	closestTToPoint(p) {
		// (at(t) - p) * tangentAt(t) = 0
		// (xi f1 + eta f2 + q) * (xi f2 - eta f1) = 0
		// xi eta (f2^2-f1^2) + xi f2 q - eta² f1 f2 + xi² f1 f2 - eta f1 q = 0
		//  (xi² - eta²) f1 f2 + xi eta (f2^2-f1^2) + xi f2 q - eta f1 q = 0

		// atan2 of p is a good first approximation for the searched t
		const startT = this.inverseMatrix.transformPoint(p).angleXY()
		const pRelCenter = p.minus(this.center)
		const f = t => this.tangentAt(t).dot(this.f1.times(Math.cos(t)).plus(this.f2.times(Math.sin(t))).minus(pRelCenter))
		return newtonIterate1d(f, startT)
	}

	/**
	 * Returns a new SemiEllipseCurve representing an ellipse parallel to the XY-plane
	 * with semi-major/minor axes parallel t the X and Y axes and of length a and b.
	 *
	 * @param a length of the axis parallel to X axis
	 * @param b length of the axis parallel to Y axis
	 * @param center Defaults to V3.ZERO
	 */
	static forAB(a: number, b: number, center: V3 = V3.ZERO): SemiEllipseCurve {
		return new SemiEllipseCurve(center, new V3(a, 0, 0), new V3(0, b, 0))
	}

	/**
	 * Returns a new SemiEllipseCurve representing a circle parallel to the XY-plane.`
	 */
	static semicircle(radius: number, center: V3 = V3.ZERO): SemiEllipseCurve {
		return new SemiEllipseCurve(center, new V3(radius, 0, 0), new V3(0, radius, 0))
	}

	area(): number {
		// see
		// https://upload.wikimedia.org/wikipedia/commons/thumb/4/4e/Cross_product_parallelogram.svg/220px-Cross_product_parallelogram.svg.png
		return Math.PI * this.f1.cross(this.f2).length()
	}

	static readonly UNIT = new SemiEllipseCurve(V3.ZERO, V3.X, V3.Y)

	angleToT(phi: number): number {
		// atan2(y, x) = phi
		const phiDir = this.f1.normalized().times(Math.cos(phi)).plus(this.f2.rejectedFrom(this.f1).normalized().times(Math.sin(phi)))
		const localDir = this.inverseMatrix.transformVector(phiDir)
		return localDir.angleXY()
	}
}
SemiEllipseCurve.prototype.hlol = Curve.hlol++
SemiEllipseCurve.prototype.tIncrement = 2 * Math.PI / (4 * 80)
SemiEllipseCurve.prototype.tMin = 0
SemiEllipseCurve.prototype.tMax = PI