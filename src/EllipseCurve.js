/**
 * Created by aval on 14/02/2016.
 */
class EllipseCurve extends Curve {
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

	toString(f) {
		return `new EllipseCurve(${this.center}, ${this.f1}, ${this.f2})`
	}

	static isValidT(t) {
		return -Math.PI <= t && t <= Math.PI
	}

	at(t) {
		return this.center.plus(this.f1.times(Math.cos(t))).plus(this.f2.times(Math.sin(t)))
	}

	at2(xi, eta) {
		return this.center.plus(this.f1.times(xi)).plus(this.f2.times(eta))
	}

	tangentAt(t) {
		assertNumbers(t)
		return this.f2.times(Math.cos(t)).minus(this.f1.times(Math.sin(t))).normalized()
	}

	tangentAt2(xi, eta) {
		return this.f2.times(xi).minus(this.f1.times(eta)).normalized()
	}

	isCircular() {
		return NLA.equals(this.f1.length(), this.f2.length())
	}


	/**
	 * @inheritDoc
	 */
	likeCurve(curve) {
		return curve.constructor == EllipseCurve
			&& this.center.like(curve.center)
			&& this.f1.like(curve.f1)
			&& this.f2.like(curve.f2)
	}

	isColinearTo(curve) {
		if (curve.constructor != EllipseCurve) {
			return false
		}
		if (!this.center.like(curve.center)) {
			return false
		}
		if (this == curve) {
			return true
		}
		if (this.isCircular()) {
			return curve.isCircular() && NLA.equals(this.f1.length(), curve.f1.length())
		} else {
			var {f1: f1, f2: f2} = this.mainAxes(), {f1: c1, f2: c2} = curve.mainAxes()
			return NLA.equals(f1.lengthSquared(), abs(f1.dot(c1)))
				&& NLA.equals(f2.lengthSquared(), abs(f2.dot(c2)))
		}
	}

	normalAt(t) {
		return this.tangentAt(t).cross(this.normal)
	}

	/**
	 *
	 * @param {V3} p
	 * @param {V3} hint TODO document hint
	 * @returns {*}
	 */
	pointLambda(p, hint) {
		assertVectors(p)
		var p2 = this.inverseMatrix.transformPoint(p)
		var angle = p2.angleXY()
		if (angle < -Math.PI + NLA.PRECISION || angle > Math.PI - NLA.PRECISION) {
			return hint.dot(this.f2) < 0
				? Math.PI
				: -Math.PI
		}
		return angle
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
	 * @returns {{f1: V3, f2: V3}}
	 */
	mainAxes() {
		var f1 = this.f1, f2 = this.f2, a = f1.dot(f2), b = f2.lengthSquared() - f1.lengthSquared()
		if (NLA.isZero(a)) {
			return this
		}
		var g1 = 2 * a, g2 = b + Math.sqrt(b * b + 4 * a * a)
		var {x1: xi, y1: eta} = intersectionUnitCircleLine(g1, g2, 0)
		return {f1: f1.times(xi).plus(f2.times(eta)), f2: f1.times(-eta).plus(f2.times(xi))}

	}

	eccentricity() {
		var mainAxes = this.mainAxes()
		var f1length = mainAxes.f1.length(), f2length = mainAxes.f1.length()
		var [a, b] = f1length > f2length ? [f1length, f2length] : [f2length, f1length]
		return Math.sqrt(1 - b * b / a / a)
	}

	circumference(p) {
		assert(false, "implementation?")
	}

	circumferenceApproximate(p) {
		// approximate circumference by Ramanujan
		// https://en.wikipedia.org/wiki/Ellipse#Circumference
		var {f1, f2} = this.mainAxes(), a = f1.length(), b = f2.length()
		var h = (a - b) * (a - b) / (a + b) / (a + b) // (a - b)² / (a + b)²
		return Math.PI * (a + b) * (1 + 3 * h / (10 + Math.sqrt(4 - 3 * h)))
	}

	transform(m4) {
		return new EllipseCurve(m4.transformPoint(this.center), m4.transformVector(this.f1), m4.transformVector(this.f2))
	}

	rightAngled() {
		var {f1, f2} = this.mainAxes()
		return new EllipseCurve(this.center, f1, f2)
	}

	/**
	 * @inheritDoc
	 */
	isTsWithSurface(surface) {
		if (surface instanceof PlaneSurface) {
			return this.isTsWithPlane(surface.plane)
		} else if (surface instanceof CylinderSurface) {
			var ellipseProjected = surface.baseEllipse.transform(M4.projection(this.getPlane(), surface.dir))
			return this.intersectWithEllipse(ellipseProjected).map(p => this.pointLambda(p))
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
		var
			n = plane.normal, w = plane.w,
			center = this.center, f1 = this.f1, f2 = this.f2,
			g1 = n.dot(f1), g2 = n.dot(f2), g3 = w - n.dot(center)

		var {x1: xi1, y1: eta1, x2: xi2, y2: eta2} = intersectionUnitCircleLine(g1, g2, g3)
		if (isNaN(xi1)) {
			return []
		}
		return [atan2(eta1, xi1), atan2(eta2, xi2)]

	}

	getPlane() {
		return P3.normalOnAnchor(this.normal, this.center)
	}

	containsPoint(p) {
		var localP = this.inverseMatrix.transformPoint(p)
		return NLA.equals(1, localP.lengthXY()) && NLA.isZero(localP.z)
	}

	intersectWithEllipse(ellipse) {
		if (this.normal.isParallelTo(ellipse.normal) && NLA.isZero(this.center.minus(ellipse.center).dot(ellipse.normal))) {

			// ellipses are coplanar
			var localEllipse = ellipse.transform(this.inverseMatrix).rightAngled()
			//new EllipseCurve(V3.ZERO, V3.X, V3.Y).debugToMesh(dMesh, 'curve4')
			console.log(localEllipse, localEllipse.ss)
			//localEllipse.debugToMesh(dMesh, 'curve3')
			var angle = localEllipse.f1.angleXY()
			console.log('angle', angle)
			var aSqr = localEllipse.f1.lengthSquared(), bSqr = localEllipse.f2.lengthSquared()
			var a = sqrt(aSqr), b = sqrt(bSqr)
			var {x: centerX, y: centerY} = localEllipse.center
			var rotCenterX = centerX * cos(-angle) + centerY * -sin(-angle)
			var rotCenterY = centerX * sin(-angle) + centerY * cos(-angle)
			var rotCenter = V3(rotCenterX, rotCenterY)
			var f = t => {
				var lex = cos(t) - rotCenterX, ley = sin(t) - rotCenterY;
				return lex * lex / aSqr + ley * ley / bSqr - 1
			}
			var uc = new EllipseCurve(V3.ZERO, V3.X, V3.Y)
			//uc.debugToMesh(dMesh, 'curve4')
			var f2 = (x, y) => 200 * (x * x + y * y - 1)
			var f3 = (x, y) => 200 * ((x - rotCenterX) * (x - rotCenterX) / aSqr + (y - rotCenterY) * (y - rotCenterY) / bSqr - 1)
			var results = []
			var resetMatrix = this.matrix.times(M4.rotationZ(angle))
			for (var da = PI / 4; da < 2 * PI; da += PI / 2) {
				var startP = uc.at(da)
				var p = newtonIterate2d(f3, f2, startP.x, startP.y, 10)
				if (p && !results.some(r => r.like(p))) {
					results.push(p)
					drPs.push(p)
				}
			}
			var rotEl = new EllipseCurve(rotCenter, V3(a, 0, 0), V3(0, b, 0))
			//var rotEl = localEllipse.transform(resetMatrix)
			console.log(rotEl, rotEl.ss)
			//rotEl.debugToMesh(dMesh, 'curve2')
			return results.map(p => resetMatrix.transformPoint(p))
			/*
			 // new rel center
			 var mat = M4.forSys(localEllipse.f1.normalized(), localEllipse.f2.normalized(), V3.Z, localEllipse.center).inversed()
			 console.log(mat.toString())
			 var newCenter = mat.transformPoint(V3.ZERO)
			 var x0 = newCenter.x, y0 = newCenter.y
			 var c = (1 - bSqr / aSqr) / 2/ y0, d = -x0 / y0, e = (bSqr + x0 * x0 + y0 * y0) / 2 / y0

			 var ff = x => c*c*x*x*x*x+ 2*c*d*x*x*x+ (2*c*e+d*d+bSqr/aSqr)*x*x+2*d*e*x+e*e-bSqr
			 var newx1 = newtonIterate(ff, x0 - 1)
			 var newx2 = newtonIterate(ff, x0)
			 var f1 = (x, y) => 2 * x - y
			 var f2 = (x, y) => 2 * ((x - x0) * (x - x0) + (y - y0) * (y - y0) - 1)
			 var f3 = (x, y) => 2 * (x * x / aSqr + y * y / bSqr - 1)
			 localEllipse = localEllipse.transform(mat)
			 for (var a = PI / 4; a < 2 * PI; a+= PI / 2) {
			 var startP = circle.at(a)
			 //drPs.push(startP)
			 var p = newtonIterate2d(f3, f2, startP.x, startP.y, 10)
			 p && drPs.push(p)
			 p && console.log(p.$, p.minus(V3(x0, y0)).length())
			 }
			 circle.debugToMesh(dMesh, 'curve1')
			 localEllipse.debugToMesh(dMesh, 'curve2')
			 */
		} else {
			assert(false)
		}
	}

	static forAB(a, b, center) {
		return new EllipseCurve(center || V3.ZERO, V3(a, 0, 0), V3(0, b, 0))
	}

	/**
	 *
	 * @param {number} radius
	 * @param {V3=} center
	 * @returns {EllipseCurve}
	 */
	static circle(radius, center) {
		return new EllipseCurve(center || V3.ZERO, V3(radius, 0, 0), V3(0, radius, 0))
	}
}
EllipseCurve.UNIT = new EllipseCurve(V3.ZERO, V3.X, V3.Y)
EllipseCurve.prototype.tIncrement = 2 * Math.PI / (4 * 8)