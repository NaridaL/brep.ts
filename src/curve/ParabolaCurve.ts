/**
 * eta = xi²
 */
class ParabolaCurve extends Curve {
	readonly normal: V3
    readonly center: V3
    readonly f1: V3
    readonly f2: V3
    readonly matrix: M4
    readonly inverseMatrix: M4

	constructor(center: V3, f1: V3, f2: V3, tMin: number = -10, tMax: number = 10) {
		super(tMin, tMax)
		assertVectors(center, f1, f2)
		this.center = center
		this.f1 = f1
		this.f2 = f2
		this.normal = f1.cross(f2).unit()
		this.matrix = M4.forSys(f1, f2, this.normal, center)
		this.inverseMatrix = this.matrix.inversed()
	}

	toString() {
		return makeGen('new ParabolaCurve', this.center, this.f1, this.f2)
	}

	at(t) {
	    // center + f1 t + f2 t²
		return this.center.plus(this.f1.times(t)).plus(this.f2.times(t * t))
	}

	at2(xi, eta) {
		return this.center.plus(this.f1.times(xi)).plus(this.f2.times(eta))
	}

	tangentAt(t: number): V3 {
		assertNumbers(t)
        // f1 + f2 2 t
		return this.f1.plus(this.f2.times(2 * t))
	}

	ddt(t: number): V3 {
		assertNumbers(t)
		return this.f2.times(2)
	}

	tangentAt2(xi, eta) {
		assertNumbers(xi, eta)
		return this.f1.plus(this.f2.times(2 * eta))
	}

    /**
     * tangent: f1 + 2 * t * f2 = 0
     * t = -f1 / 2 / f2 (for individual dimensions)
     */
    roots(): number[][] {
        return NLA.arrayFromFunction(3, dim => eq0(this.f2.e(dim)) ? [] : [-this.f1.e(dim) / 2 / this.f2.e(dim)])
    }

    equals(obj: any): boolean {
        return this == obj ||
            Object.getPrototypeOf(obj) == ParabolaCurve.prototype
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

	isColinearTo(curve: Curve): boolean {
		if (curve.constructor != ParabolaCurve) {
			return false
		}
		const mainAxes = this.rightAngled(), curveMainAxes = curve.rightAngled()
		return mainAxes.center.like(curveMainAxes.center)
			&& mainAxes.f2.like(curveMainAxes.f2)
			&& mainAxes.f1.likeOrReversed(curveMainAxes.f1)
	}

	normalAt(t: number): V3 {
		return this.tangentAt(t).cross(this.normal)
	}

	pointT(p: V3): number {
		assertVectors(p)
		return this.inverseMatrix.transformPoint(p).x
	}

	isOrthogonal(p: V3): boolean {
		return this.f1.isPerpendicularTo(this.f2)
	}

	rightAngled() {
		// looking for vertex of parabola
		// this is the point where the tangent is perpendicular to the main axis (f2)
		// tangent = f1 + f2 * 2 * t0
		// f2 DOT (f1 + f2 * 2 * t0) == 0
		// f1 DOT f2 + f2 DOT f2 * 2 * t0 == 0
		// t0 == -(f1 DOT f2) / (f2 DOT f2 * 2)
		const f1 = this.f1, f2 = this.f2
        const f1DOTf2 = f1.dot(f2)
		if (NLA.eq0(f1DOTf2)) {
			return this
		}
		const t0 = -f1DOTf2 / f2.squared() / 2
		// we need to rearange tMin/tMax
		// tMin' = pointT(at(tMin)) =
		const raCenter = this.at(t0)
		const raF1 = this.tangentAt(t0)
		const repos = t => this.at(t).minus(raCenter).dot(raF1) / raF1.squared()
		return new ParabolaCurve(raCenter, raF1, f2, repos(this.tMin), repos(this.tMax))
	}

	arcLength(startT: number, endT: number): number {
		let f1 = this.f1
		const f2 = this.f2
		let f1DOTf2 = f1.dot(f2), t0 = 0
		if (!NLA.eq0(f1DOTf2)) {
			t0 = -f1DOTf2 / f2.squared() / 2
			f1 = f1.plus(f2.times(2 * t0))
		}
		const f1Length = f1.length()
		const a = f2.length() / f1Length

		function F(x) {
			return Math.asinh(a * 2 * x) / 4 / a + x * Math.sqrt(1 + a * a * 4 * x * x) / 2
		}

		return f1Length * (F(endT - t0) - F(startT - t0))
	}

	transform(m4) {
		return new ParabolaCurve(
			m4.transformPoint(this.center),
			m4.transformVector(this.f1),
			m4.transformVector(this.f2),
			this.tMin, this.tMax) as this
	}

	static eccentricity() {
		return 1
	}

	isTsWithSurface(surface: Surface) {
		if (surface instanceof PlaneSurface) {
			return this.isTsWithPlane(surface.plane)
		} else if (surface instanceof ConicSurface) {
			const ParabolaProjected = surface.baseParabola.transform(M4.projection(this.getPlane(), surface.dir))
			return this.intersectWithParabola(ParabolaProjected).map(p => this.pointT(p))
		} else {
			assert(false)
		}
	}

	isInfosWithLine(line: L3) {
		const anchorLC = this.inverseMatrix.transformPoint(line.anchor)
		const dirLC = this.inverseMatrix.transformVector(line.dir1)
		if (NLA.eq0(dirLC.z)) {
			// local line parallel to XY-plane
			if (NLA.eq0(anchorLC.z)) {
				// local line lies in XY-plane
				// para: x² = y
				// line(t) = anchor + t dir
				// (ax + t dx)² = ay + t dy
				// ax² + t ax dx + t² dx² = ay + t dy
				// t² dx² + t (ax dx + dy) + ay² + ay = 0
				const pqDiv = dirLC.x ** 2
				const lineTs = pqFormula((anchorLC.x * dirLC.x + dirLC.y) / pqDiv, (anchorLC.x ** 2 + anchorLC.y) / pqDiv)
				return lineTs.filter(tOther => le(0, anchorLC.y + tOther * dirLC.y))
					.map(tOther => ({
						tThis: dirLC.x * tOther + anchorLC.x,
						tOther: tOther,
						p: line.at(tOther)}))
			}
		} else {
			// if the line intersects the XY-plane in a single point, there can be an intersection there
			// find point, then check if distance from circle = 1
			const otherTAtZ0 = anchorLC.z / dirLC.z
			const isp = dirLC.times(otherTAtZ0).plus(anchorLC)
			if (fuzzyBetween(isp.x, this.tMin, this.tMax)) {
				// point lies on unit circle
				return [{
					tThis: isp.x,
					tOther: otherTAtZ0,
					p: line.at(otherTAtZ0)}]
			}
		}
		return []
	}

	isTsWithPlane(plane: P3) {
		assertInst(P3, plane)
		/*
		 this: x = center + f1 * cos t + f2 * sin t  (1)
		 plane:
		 n := plane.normal1
		 n DOT x == plane.w           (2)
		 plane defined by f1/f2
		 x = center + f1 * xi + f2 * eta         (3)
		 intersection plane and planef1/f2:
		 insert (3) into (2):
		 n DOT center + n DOT f1 * xi + n DOT f2 * eta = plane.w | -n DOT center
		 n DOT f1 * xi + n DOT f2 * eta = plane.w - n DOT center (4)
		 points on Parabola have additional condition
		 eta = xi * xi (5)
		 g1 := n DOT f1
		 g2 := n DOT f2
		 g3 := w - n DOT center
		 solve system (5)/(6)
		 g1 * xi + g2 * eta = g3 (6)
		 g1 * xi + g2 * xi * xi = g3
		 xi² + xi * g1/g2 - g3/g2 = 0
		 */
		if (plane.normal1.isParallelTo(this.normal)) {
			return []
		}
		// funnily enough, changing the order of the operations changes nothing...
		const n = plane.normal1, w = plane.w,
			g1 = n.dot(this.f1), g2 = n.dot(this.f2)
		let g3 = w - n.dot(this.center)
		// g2 not zero (!plane.normal1.isParallelTo(this.normal1))
		let p = g1 / g2
		const q = -g3 / g2
		const discriminant4 = p * p / 4 - q
		if (discriminant4 < -NLA_PRECISION) {
			return []
		} else if (discriminant4 <= NLA_PRECISION) {
			return [-p / 2]
		} else {
			const root = Math.sqrt(discriminant4)
			return [-p / 2 - root, -p / 2 + root]
		}
	}

	getPlane() {
		return P3.normalOnAnchor(this.normal, this.center)
	}

	containsPoint(p) {
		const pLC = this.inverseMatrix.transformPoint(p)
		return NLA.eq(pLC.x ** 2, pLC.y)
	}

    static forAB(a: V3, b: V3, center: V3 = V3.O): ParabolaCurve {
        return new ParabolaCurve(center, V(a, 0, 0), V(0, b, 0))
    }

	static quadratic(a: V3, b: V3, c: V3): ParabolaCurve {
        // (1 - t)² a + 2 * t * (1 - t) b + t² c
        // (1 -2t +t²)a + (2t -2t²) b + t² c
        // = t²(a - 2b + c) + t (-2a + 2b) + a
        // (2t - 2) a + (1 - 2t) b + 2t c = t(2a + 2b - 2c) - 2a + b
        // 2 a + -2 b + 2 c
        const f2 = a.plus(c).minus(b.times(2))
        const f1 = b.minus(a).times(2)
        const center = a
        return new ParabolaCurve(center, f1, f2, 0, 1)
    }

    asBezier() {
	    return BezierCurve.quadratic(
	        this.at(-1),
            new L3(this.at(-1), this.tangentAt(-1).unit()).isInfoWithLine(new L3(this.at(1), this.tangentAt(1).unit())),
            this.at(1))
    }

	static readonly XY = new ParabolaCurve(V3.O, V3.X, V3.Y)
	static readonly YZ = new ParabolaCurve(V3.O, V3.Y, V3.Z)
	static readonly ZX = new ParabolaCurve(V3.O, V3.Z, V3.X)
}
ParabolaCurve.prototype.tIncrement = 1/32
