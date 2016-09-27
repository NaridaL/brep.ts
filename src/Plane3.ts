"use strict"

/**
 * Created by aval on 20/11/2015.
 */

class P3 extends Transformable {
	w: number
	normal: V3

	/**
	 * Oriented plane, i.e. splits R^3 in half, with one half being "in front" of the plane.
	 * Leads to multiple comparisons: isCoplanarToPlane returns if the plane occupies the same space,
	 * like returns if the plane occupies the same space and has the same orientation
	 *
	 * Points x on the plane fulfill the equation: normal DOT x = w
	 *
	 * @param normal1 normalized plane normal
	 * @param w signed (rel to normal1) distance from the origin
	 * @param prototype object to set as prototype of the new Plane3 object. Defaults to Plane3.prototype
	 * @constructor P3
	 * @augments {Transformable}
	 * @property {V3} normal
	 * @property {number} w
	 */
	constructor(normal1:V3, w:number, prototype?) {

		assertVectors(normal1)
		assertNumbers(w)
		assert(normal1.hasLength(1), "normal1.hasLength(1)" + normal1)
		var p = Object.create(prototype || P3.prototype)
		p.w = w
		p.normal = normal1
		return p
	}

	/**
	 *
	 * @returns {V3}
	 */
	axisIntercepts() {
		let w = this.w, n = this.normal
		return new V3(w / n.x, w / n.y, w / n.z)
	}

	/**
	 * 	 * @returns {V3}
	 */
	get anchor():V3 {
		return this.normal.times(this.w)
	}

	isCoplanarToPlane(plane) {
		assertInst(P3, plane)
		return this.like(plane) || this.likeFlipped(plane)
	}

	like(plane) {
		assertInst(P3, plane)
		return NLA.equals(this.w, plane.w) && this.normal.like(plane.normal)
	}

	likeFlipped(plane) {
		assertInst(P3, plane)
		return NLA.equals(this.w, -plane.w) && this.normal.like(plane.normal.negated())
	}

	/**
	 * True iff plane.normal is equal to this.normal or it's negation.
	 *
	 * @param plane
	 */
	isParallelToPlane(plane:P3) {
		assertInst(P3, plane)
		return NLA.equals(1, Math.abs(this.normal.dot(plane.normal)))
	}

	isParallelToLine(line) {
		assertInst(L3, line)
		return NLA.isZero(this.normal.dot(line.dir1))
	}

	isPerpendicularToLine(line) {
		assertInst(L3, line)
		// this.normal || line.dir1
		return NLA.equals(1, Math.abs(this.normal.dot(line.normal)))
	}

	isPerpendicularToPlane(plane) {
		assertInst(P3, plane)
		return NLA.isZero(this.normal.dot(plane.normal))
	}

	toString(roundFunction) {
		roundFunction = roundFunction || (v => v) //((v) => +v.toFixed(3))
		return "new P3("+this.normal.toString(roundFunction) + ", " + roundFunction(this.w) +")"
	}

	translated(offset) {
		return new P3(this.normal, this.w + offset.dot(this.normal))
	}

	transform(m4) {
		var mirror = m4.isMirroring()
		// get two vectors in the plane:
		var u = this.normal.getPerpendicular()
		var v = u.cross(this.normal)
		// get 3 points in the plane:
		var p1 = m4.transformPoint(this.anchor),
			p2 = m4.transformPoint(this.anchor.plus(v)),
			p3 = m4.transformPoint(this.anchor.plus(u))
		// and create a new plane from the transformed points:
		return P3.throughPoints(p1, !mirror ? p2 : p3, !mirror ? p3 : p2)
	}

	distanceToLine(line) {
		assertInst(L3, line)
		if (!this.isParallelToLine(line)) {
			return this.distanceToPoint(line.anchor)
		} else {
			return 0
		}
	}

	containsPoint(x) {
		assertVectors(x)
		return NLA.equals(this.w, this.normal.dot(x))
	}

	containsLine(line) {
		assertInst(L3, line)
		return this.containsPoint(line.anchor) && this.isParallelToLine(line)
	}

	distanceToPointSigned(point) {
		assertInst(V3, point)
		return this.normal.dot(point) - this.w
	}

	distanceToPoint(point) {
		assertInst(V3, point)
		return Math.abs(this.normal.dot(point) - this.w)
	}

	intersectionWithLine(line) {
		line.intersectionWithPlane(this)
	}

	/**
	 *
	 * @param plane
	 * @returns {L3}
	 */
	intersectionWithPlane(plane) {
		/*

		 this: n0 * x = w0
		 plane: n1 * x = w1
		 plane perpendicular to both which goes through origin:
		 n2 := n0 X x1
		 n2 * x = 0
		 */
		assertInst(P3, plane)
		assert(!this.isParallelToPlane(plane), "!this.isParallelToPlane(plane)")
		/*
		 var n0 = this.normal, n1 = plane.normal, n2 = n0.cross(n1).normalized(), m = M4.forSys(n0, n1, n2)
		 var x0 = this.anchor, x1 = plane.anchor, x2 = V3.ZERO
		 var p = n2.times(x2.dot(n2))
		 .plus(n1.cross(n2).times(x0.dot(n0)))
		 .plus(n2.cross(n0).times(x1.dot(n1)))
		 .div(m.determinant())
		 */
		var n0 = this.normal, n1 = plane.normal, n2 = n0.cross(n1).normalized()
		var p = M4.forRows(n0, n1, n2).inversed().transformVector(new V3(this.w, plane.w, 0))
		return new L3(p, n2)
	}

	/**
	 * Returns the point in the plane closest to the given point
	 *
	 * @param x
	 * @returns {V3}
	 */
	projectedPoint(x) {
		// See http://math.stackexchange.com/questions/444968/project-a-point-in-3d-on-a-given-plane
		// p = x - ((x - planeAnchor) * normal) * normal
		return x.minus(this.normal.times(x.minus(this.anchor).dot(this.normal)))
	}

	/**
	 *
	 * @param x
	 * @returns {V3}
	 */
	projectedVector(x) {
		// See V3.rejectedFrom. Simplified, as this.normal.length() == 1
		return x.minus(this.normal.times(x.dot(this.normal)))
	}

	/**
	 *
	 * @returns {P3}
	 */
	flipped() {
		return new P3(this.normal.negated(), -this.w)
	}
	
	static throughPoints(a:V3, b:V3, c:V3, prototype?) {
		assertVectors(a, b, c)
		var n1 = b.minus(a).cross(c.minus(a)).normalized();
		return new P3(n1, n1.dot(a), prototype)
	}

	/**
	 * @alias P3.normalOnAnchor
	 * @param normal
	 * @param anchor
	 * @param prototype
	 * @returns {P3}
	 */
	static normalOnAnchor(normal:V3, anchor:V3, prototype?):P3 {
		assertVectors(normal, anchor)
		var n1 = normal.normalized()
		return new P3(n1, n1.dot(anchor), prototype)
	}

	/**
	 * x/x0 + y/y0 + y/y0 = 1
	 *
	 * @param x0
	 * @param y0
	 * @param z0
	 */
	static forAxisIntercepts(x0:number, y0:number, z0:number):P3 {
		assertNumbers(x0, y0, z0)
		let normal = new V3(1 / x0, 1 / y0, 1 / z0)
		return new P3(normal.normalized(), normal.length())
	}

	/**
	 * @alias P3.forAnchorAndPlaneVectors
	 * @param anchor
	 * @param v0
	 * @param v1
	 * @param prototype
	 * @returns {P3}
	 */
	static forAnchorAndPlaneVectors(anchor:V3, v0:V3, v1:V3, prototype):P3 {
		assertVectors(anchor, v0, v1)
		return P3.normalOnAnchor(v0.cross(v1), anchor, prototype)
	}


// X-Y-Z planes
	static YZ = new P3(V3.X, 0)
	static ZX = new P3(V3.Y, 0)
	static XY = new P3(V3.Z, 0)
}
