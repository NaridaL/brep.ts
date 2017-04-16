class SemiEllipseCurve extends Curve {
	readonly normal: V3
	readonly matrix: M4
	readonly inverseMatrix: M4

	constructor(readonly center: V3,
	            readonly f1: V3,
	            readonly f2: V3,
	            readonly tMin: number = 0,
	            readonly tMax: number = PI) {
        super(tMin, tMax)
        assertVectors(center, f1, f2)
        assert(0 <= this.tMin && this.tMin < PI)
        assert(0 < this.tMax && this.tMax <= PI)
		this.normal = f1.cross(f2)
		if (!this.normal.isZero()) {
			this.normal = this.normal.unit()
			this.matrix = M4.forSys(f1, f2, this.normal, center)
			this.inverseMatrix = this.matrix.inversed()
		} else {
			this.matrix = M4.forSys(f1, f2, f1.unit(), center)
			const f1p = f1.getPerpendicular()
			this.inverseMatrix = new M4(
				1, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 1).times(M4.forSys(f1, f1p, f1.cross(f1p), center).inversed())
		}
	}

	getAreaInDir(right: V3, up: V3, tStart: number, tEnd: number): {area: number, centroid: V3} {
		return EllipseCurve.prototype.getAreaInDir.call(this, right, up, tStart, tEnd)
	}

	toSource() {
		return makeGen('new SemiEllipseCurve', this.center, this.f1, this.f2, this.tMin, this.tMax)
	}

	at(t: number): V3 {
        assertNumbers(t)
        //assert(this.isValidT(t))
        // center + f1 cos t + f2 sin t
		return this.center.plus(this.f1.times(Math.cos(t))).plus(this.f2.times(Math.sin(t)))
	}

	at2(xi: number, eta: number): V3 {
		// center + f1 xi + f2 eta
		return this.center.plus(this.f1.times(xi)).plus(this.f2.times(eta))
	}

	tangentAt(t: number): V3 {
		assertNumbers(t)
        //assert(this.isValidT(t))
		// f2 cos(t) - f1 sin(t)
		return this.f2.times(Math.cos(t)).minus(this.f1.times(Math.sin(t)))
	}

	ddt(t: number): V3 {
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

	likeCurve(curve: Curve): boolean {
		return curve.constructor == SemiEllipseCurve
			&& this.center.like(curve.center)
			&& this.f1.like(curve.f1)
			&& this.f2.like(curve.f2)
	}

	isColinearTo(curve: Curve): boolean {
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
			return NLA.eq(f1.squared(), Math.abs(f1.dot(c1)))
				&& NLA.eq(f2.squared(), Math.abs(f2.dot(c2)))
		}
	}

	normalAt(t: number): V3 {
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

	reversed(): SemiEllipseCurve {
		return new SemiEllipseCurve(this.center, this.f1.negated(), this.f2, PI - this.tMax, PI - this.tMin)
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

	transform(m4: M4): this {
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

    isTsWithSurface(surface): number[] {
        if (surface instanceof PlaneSurface) {
            return this.isTsWithPlane(surface.plane)
        } else if (surface instanceof SemiCylinderSurface) {
            const infos = surface.isCurvesWithPlane(this.getPlane())
	            .map(isc => this.isInfosWithCurve(isc))
	            .concatenated()
            return infos.map(info => info.tThis)
        } else if (surface instanceof EllipsoidSurface) {
            const isEllipse = surface.isCurvesWithPlane(this.getPlane())
            if (isEllipse.length < 1) return
            const infos = this.isInfosWithEllipse(isEllipse[0] as EllipseCurve)
            return infos.map(info => info.tThis)
        } else if (surface instanceof SemiEllipsoidSurface) {
            const isEllipse = surface.asEllipsoidSurface().isCurvesWithSurface(new PlaneSurface(this.getPlane()))
            if (isEllipse.length < 1) return []
            const possibleInfos = this.isInfosWithEllipse(isEllipse[0] as EllipseCurve)
            return possibleInfos.filter(info => surface.containsPoint(info.p)).map(info => info.tThis)
        } else if (surface instanceof ProjectedCurveSurface) {
            return surface.isCurvesWithPlane(this.getPlane())
                .map(curve => this.isInfosWithCurve(curve))
                .concatenated()
                .map(info => info.tThis)
        } else {
			assert(false)
		}
	}

	isTsWithPlane(plane: P3) {
		assertInst(P3, plane)
		/*
		 this: x = center + f1 * cos t + f2 * sin t  (1)
		 n DOT center + n DOT f1 cos t + n DOT f2 sin t - w = 0
		 plane:
		 n := plane.normal1
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
		if (plane.normal1.isParallelTo(this.normal)) {
			return []
		}
		const n = plane.normal1, w = plane.w,
			center = this.center, f1 = this.f1, f2 = this.f2,
			g1 = n.dot(f1), g2 = n.dot(f2), g3 = w - n.dot(center)

		const isLC = intersectionUnitCircleLine2(g1, g2, g3)
        const result = []
		let t
		for (const [xi, eta] of isLC) {
			le(0, eta) && this.isValidT(t = SemiEllipseCurve.unitT(new V3(xi, eta, 0))) && result.push(t)
		}
        for (const t of result) {
		    assert(plane.containsPoint(this.at(t)))
        }
		return result

	}

    getPlane(): P3 {
        return P3.normalOnAnchor(this.normal, this.center)
    }

	containsPoint(p: V3): boolean {
		const pLC = this.inverseMatrix.transformPoint(p)
		return NLA.eq0(pLC.z) && SemiEllipseCurve.validPlanePoint(pLC.x, pLC.y)
	}

	asEllipse(): EllipseCurve {
		return new EllipseCurve(this.center, this.f1, this.f2, this.tMin, this.tMax)
	}

	isInfosWithEllipse(ellipse: EllipseCurve | SemiEllipseCurve): {tThis: number, tOther: number, p: V3}[] {
		if (this.normal.isParallelTo(ellipse.normal) && NLA.eq0(this.center.minus(ellipse.center).dot(ellipse.normal))) {
			ellipse instanceof SemiEllipseCurve && (ellipse = ellipse.asEllipse())
			return this.asEllipse().isInfosWithCurve(ellipse).filter(info => this.isValidT(info.tThis) && ellipse.isValidT(info.tOther))
		} else {
			return this.isTsWithPlane(P3.normalOnAnchor(ellipse.normal.unit(), ellipse.center)).mapFilter(t => {
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

	isInfosWithCurve(curve: Curve): { tThis: number, tOther: number, p: V3 }[] {
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

    isPointsWithBezier(bezier: BezierCurve): V3[] {
        const bezierLC = bezier.transform(this.inverseMatrix)
        if (new PlaneSurface(P3.XY).containsCurve(bezier)) {
            // up to 6 solutions possible
            const f = t => bezierLC.at(t).squaredXY() - 1
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
			const a = this.f2.e(dim), b = -this.f1.e(dim)
			const {x1,y1,x2,y2} = intersectionUnitCircleLine(a, b, 0)
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
	 * @param center Defaults to V3.O
	 */
	static forAB(a: number, b: number, center: V3 = V3.O): SemiEllipseCurve {
		return new SemiEllipseCurve(center, new V3(a, 0, 0), new V3(0, b, 0))
	}

	/**
	 * Returns a new SemiEllipseCurve representing a circle parallel to the XY-plane.`
	 */
	static semicircle(radius: number, center: V3 = V3.O): SemiEllipseCurve {
		return new SemiEllipseCurve(center, new V3(radius, 0, 0), new V3(0, radius, 0))
	}

	area(): number {
		// see
		// https://upload.wikimedia.org/wikipedia/commons/thumb/4/4e/Cross_product_parallelogram.svg/220px-Cross_product_parallelogram.svg.png
		return Math.PI * this.f1.cross(this.f2).length()
	}

	static readonly UNIT = new SemiEllipseCurve(V3.O, V3.X, V3.Y)

	angleToT(phi: number): number {
		// atan2(y, x) = phi
		const phiDir = this.f1.unit().times(Math.cos(phi)).plus(this.f2.rejectedFrom(this.f1).unit().times(Math.sin(phi)))
		const localDir = this.inverseMatrix.transformVector(phiDir)
		return localDir.angleXY()
	}
}
SemiEllipseCurve.prototype.hlol = Curve.hlol++
SemiEllipseCurve.prototype.tIncrement = 2 * Math.PI / (4 * 32)
SemiEllipseCurve.prototype.tMin = 0
SemiEllipseCurve.prototype.tMax = PI