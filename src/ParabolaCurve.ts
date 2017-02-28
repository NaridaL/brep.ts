/**
 * eta = xi²
 */
class ParabolaCurve extends Curve {
	normal: V3
	center: V3
	f1: V3
	f2: V3
	matrix: M4
	inverseMatrix: M4

	constructor(center, f1, f2) {
		super()
		assertVectors(center, f1, f2)
		this.center = center
		this.f1 = f1
		this.f2 = f2
		this.normal = f1.cross(f2).normalized()
		this.matrix = M4.forSys(f1, f2, this.normal, center)
		this.inverseMatrix = this.matrix.inversed()
	}

	toString() {
		return makeGen('new ParabolaCurve', this.center, this.f1, this.f2)
	}

	at(t) {
		return this.center.plus(this.f1.times(t)).plus(this.f2.times(t * t))
	}

	at2(xi, eta) {
		return this.center.plus(this.f1.times(xi)).plus(this.f2.times(eta))
	}

	tangentAt(t) {
		assertNumbers(t)
		return this.f1.plus(this.f2.times(2 * t))
	}

	ddt(t) {
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
            Object.getPrototypeOf(obj) == EllipseCurve.prototype
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

	isColinearTo(curve) {
		if (curve.constructor != ParabolaCurve) {
			return false
		}
		const mainAxes = this.rightAngled(), curveMainAxes = curve.rightAngled()
		return mainAxes.center.like(curveMainAxes.center)
			&& mainAxes.f2.like(curveMainAxes.f2)
			&& mainAxes.f1.likeOrReversed(curveMainAxes.f1)
	}

	normalAt(t) {
		return this.tangentAt(t).cross(this.normal)
	}

	pointLambda(p) {
		assertVectors(p)
		return this.inverseMatrix.transformPoint(p).x
	}

	isOrthogonal(p) {
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
		let f1DOTf2 = f1.dot(f2)
		if (NLA.eq0(f1DOTf2)) {
			return this
		}
		const t0 = -f1DOTf2 / f2.squared() / 2
		// can't use .tangentAt as that gets normalized
		return new ParabolaCurve(this.at(t0), f1.plus(f2.times(2 * t0)), f2)
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
		return new ParabolaCurve(m4.transformPoint(this.center), m4.transformVector(this.f1), m4.transformVector(this.f2))
	}

	static eccentricity() {
		return 1
	}

	/**
	 * @inheritDoc
	 */
	isTsWithSurface(surface) {
		if (surface instanceof PlaneSurface) {
			return this.isTsWithPlane(surface.plane)
		} else if (surface instanceof ConicSurface) {
			const ParabolaProjected = surface.baseParabola.transform(M4.projection(this.getPlane(), surface.dir))
			return this.intersectWithParabola(ParabolaProjected).map(p => this.pointLambda(p))
		} else {
			assert(false)
		}
	}

	/**
	 * @inheritDoc
	 */
	isTsWithPlane(plane) {
		assertInst(P3, plane)
		/*
		 this: x = center + f1 * cos t + f2 * sin t  (1)
		 plane:
		 n := plane.normal
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
		if (plane.normal.isParallelTo(this.normal)) {
			return []
		}
		// funnily enough, changing the order of the operations changes nothing...
		const n = plane.normal, w = plane.w,
			g1 = n.dot(this.f1), g2 = n.dot(this.f2)
		let g3 = w - n.dot(this.center)
		// g2 not zero (!plane.normal.isParallelTo(this.normal))
		let p = g1 / g2
		const q = -g3 / g2
		const discriminant4 = p * p / 4 - q
		console.log('pq', p, q, discriminant4)
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
		const localP = this.inverseMatrix.transformPoint(p)
		return NLA.eq(localP.x * localP.x, localP.y)
	}

	static forAB(a, b, center):ParabolaCurve {
		return new ParabolaCurve(center || V3.ZERO, V(a, 0, 0), V(0, b, 0))
	}

	static XY = new ParabolaCurve(V3.ZERO, V3.X, V3.Y)
	static YZ = new ParabolaCurve(V3.ZERO, V3.Y, V3.Z)
	static ZX = new ParabolaCurve(V3.ZERO, V3.Z, V3.X)
}
