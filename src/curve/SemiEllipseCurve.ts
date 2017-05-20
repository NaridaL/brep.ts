///<reference path="../../node_modules/ts3dutils/out/complete.d.ts"/>


class SemiEllipseCurve extends XiEtaCurve {
	constructor(center: V3, f1: V3, f2: V3, tMin: number = 0, tMax: number = PI) {
		super(center, f1, f2, tMin, tMax)
		assert(0 <= this.tMin && this.tMin < PI)
		assert(0 < this.tMax && this.tMax <= PI)
	}

	getAreaInDir(right: V3, up: V3, tStart: number, tEnd: number): {area: number, centroid: V3} {
		return EllipseCurve.prototype.getAreaInDir.call(this, right, up, tStart, tEnd)
	}

	at(t: number): V3 {
        assertNumbers(t)
        //assert(this.isValidT(t))
        // center + f1 cos t + f2 sin t
		return this.center.plus(this.f1.times(Math.cos(t))).plus(this.f2.times(Math.sin(t)))
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

	isCircular(): boolean {
		return eq(this.f1.length(), this.f2.length()) && this.f1.isPerpendicularTo(this.f2)
	}

	isColinearTo(curve: Curve): boolean {
		if (!((x): x is SemiEllipseCurve => x.constructor == this.constructor)(curve)) {
			return false
		}
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
		const {x, y} = pLC
		return le(0, y) && eq0(x ** 2 + y ** 2 - 1)
	}

	static XYLCPointT(pLC: V3): number {
		assert(le(0, pLC.y))
		const angle = Math.atan2(pLC.y, pLC.x)
		return angle < 0 ? (assert(eq0(angle) || eq(PI, abs(angle))) && abs(angle)) : angle
	}

	isValidT(t: number) {
		return le(0, t) && le(t, PI)
	}

	pointT(p: V3) {
		assertVectors(p)
        assert(this.containsPoint(p))
		const pLC = this.inverseMatrix.transformPoint(p)
		const t = SemiEllipseCurve.XYLCPointT(pLC)
		assert(this.isValidT(t))
        return t
	}

	reversed(): SemiEllipseCurve {
		return new SemiEllipseCurve(this.center, this.f1.negated(), this.f2, PI - this.tMax, PI - this.tMin)
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
		const {f1, f2} = this.rightAngled(), a = f1.length(), b = f2.length()
		const h = (a - b) * (a - b) / (a + b) / (a + b) // (a - b)² / (a + b)²
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
    rightAngled(): SemiEllipseCurve {
	    const f1 = this.f1, f2 = this.f2, a = f1.dot(f2), b = f2.squared() - f1.squared()
	    if (eq0(a)) {
		    return this
	    }
	    const g1 = 2 * a, g2 = b + Math.sqrt(b * b + 4 * a * a)
	    const {x1: xi, y1: eta} = intersectionUnitCircleLine(g1, g2, 0)
	    const f1RA = f1.times(xi).plus(f2.times(eta))
	    const f2RA = f1.times(-eta).plus(f2.times(xi))
        return new SemiEllipseCurve(this.center, f1RA, f2RA)
    }

	static magic(a: number, b: number, c: number): number[] {
		const isLC = intersectionUnitCircleLine2(a, b, c)
		const result = []
		let t
		for (const [xi, eta] of isLC) {
			le(0, eta) && result.push(SemiEllipseCurve.XYLCPointT(new V3(xi, eta, 0)))
		}
		return result
	}

	asEllipse(): EllipseCurve {
		return new EllipseCurve(this.center, this.f1, this.f2, this.tMin, this.tMax)
	}

	isInfosWithEllipse(ellipse: EllipseCurve | SemiEllipseCurve): ISInfo[] {
		if (this.normal.isParallelTo(ellipse.normal) && eq0(this.center.minus(ellipse.center).dot(ellipse.normal))) {
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

	static unitIsInfosWithLine(anchorLC: V3, dirLC: V3, anchorWC: V3, dirWC: V3): ISInfo[] {
		// ell: x² + y² = 1 = p²
		// line(t) = anchor + t dir
		// anchor² - 1 + 2 t dir anchor + t² dir² = 0
		const pqDiv = dirLC.squared()
		const lineTs = pqFormula(2 * dirLC.dot(anchorLC) / pqDiv, (anchorLC.squared() - 1) / pqDiv)
		return lineTs.filter(tOther => le(0, anchorLC.y + tOther * dirLC.y))
			.map(tOther => ({
				tThis: SemiEllipseCurve.XYLCPointT(dirLC.times(tOther).plus(anchorLC)),
				tOther: tOther,
				p: L3.at(anchorWC, dirWC, tOther)}))
	}

	isInfosWithCurve(curve: Curve): ISInfo[] {
        if (curve instanceof SemiEllipseCurve || curve instanceof EllipseCurve) {
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

	closestTToPoint(p: V3) {
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

	static fromEllipse(curve: EllipseCurve, tMin: number, tMax: number): SemiEllipseCurve[] {
		return [
			tMin < 0 && new SemiEllipseCurve(curve.center, curve.f1.negated(), curve.f2.negated(), tMin + PI, min(0, tMax) + PI),
			tMax > 0 && new SemiEllipseCurve(curve.center, curve.f1, curve.f2, max(0, tMin), tMax)
		].filter(x => x)
	}
}
SemiEllipseCurve.prototype.hlol = Curve.hlol++
SemiEllipseCurve.prototype.tIncrement = 2 * Math.PI / (4 * 32)