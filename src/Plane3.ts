import between = NLA.between
import fuzzyBetween = NLA.fuzzyBetween
import clamp = NLA.clamp
import eq0 = NLA.eq0
import eq = NLA.eq
import lt = NLA.lt
import le = NLA.le
class P3 extends Transformable {

	/**
	 * Oriented plane, i.e. splits R^3 in half, with one half being "in front" of the plane.
	 * Leads to multiple comparisons: isCoplanarToPlane returns if the plane occupies the same space,
	 * like returns if the plane occupies the same space and has the same orientation
	 *
	 * Points x on the plane fulfill the equation: normal1 DOT x = w
	 *
	 * @param normal1 unit plane normal1
	 * @param w signed (rel to normal1) distance from the origin
	 */
	constructor(readonly normal1: V3,
	            readonly w: number = 0) {
		super()
		assertVectors(normal1)
		assertNumbers(w)
		assert(normal1.hasLength(1), 'normal1.hasLength(1)' + normal1)
	}

	axisIntercepts(): V3 {
		const w = this.w, n = this.normal1
		return new V3(w / n.x, w / n.y, w / n.z)
	}

	get anchor(): V3 {
		return this.normal1.times(this.w)
	}

	isCoplanarToPlane(plane: P3): boolean {
		assertInst(P3, plane)
		return this.like(plane) || this.likeFlipped(plane)
	}

	like(plane): boolean {
		assertInst(P3, plane)
		return eq(this.w, plane.w) && this.normal1.like(plane.normal1)
	}

	likeFlipped(plane): boolean {
		assertInst(P3, plane)
		return eq(this.w, -plane.w) && this.normal1.like(plane.normal1.negated())
	}

	/**
	 * True iff plane.normal1 is equal to this.normal1 or it's negation.
	 *
	 */
	isParallelToPlane(plane: P3): boolean {
		assertInst(P3, plane)
		return eq(1, Math.abs(this.normal1.dot(plane.normal1)))
	}

	isParallelToLine(line: L3): boolean {
		assertInst(L3, line)
		return eq0(this.normal1.dot(line.dir1))
	}

	isPerpendicularToLine(line: L3): boolean {
		assertInst(L3, line)
		// this.normal1 || line.dir1
		return eq(1, Math.abs(this.normal1.dot(line.dir1)))
	}

	isPerpendicularToPlane(plane: P3): boolean {
		assertInst(P3, plane)
		return eq0(this.normal1.dot(plane.normal1))
	}

	toString(roundFunction?) {
		roundFunction = roundFunction || (v => v) //((v) => +v.toFixed(3))
		return 'new P3(' + this.normal1.toString(roundFunction) + ', ' + roundFunction(this.w) + ')'
	}

	translated(offset: V3): P3 {
		return new P3(this.normal1, this.w + offset.dot(this.normal1))
	}

	transform(m4: M4): this {
		const mirror = m4.isMirroring()
		// get two vectors in the plane:
		const u = this.normal1.getPerpendicular()
		const v = u.cross(this.normal1)
		// get 3 points in the plane:
		const p1 = m4.transformPoint(this.anchor),
			p2 = m4.transformPoint(this.anchor.plus(v)),
			p3 = m4.transformPoint(this.anchor.plus(u))
		// and create a new plane from the transformed points:
		return P3.throughPoints(p1, !mirror ? p2 : p3, !mirror ? p3 : p2) as this
	}

	distanceToLine(line): number {
		assertInst(L3, line)
		if (!this.isParallelToLine(line)) {
			return this.distanceToPoint(line.anchor)
		} else {
			return 0
		}
	}

	containsPoint(x: V3): boolean {
		assertVectors(x)
		return eq(this.w, this.normal1.dot(x))
	}

	containsLine(line: L3): boolean {
		assertInst(L3, line)
		return this.containsPoint(line.anchor) && this.isParallelToLine(line)
	}

	distanceToPointSigned(point: V3): number {
		assertInst(V3, point)
		return this.normal1.dot(point) - this.w
	}

	distanceToPoint(point: V3): number {
		assertInst(V3, point)
		return Math.abs(this.normal1.dot(point) - this.w)
	}

	intersectionWithLine(line: L3): V3 {
		return line.intersectionWithPlane(this)
	}

	intersectionWithPlane(plane: P3): L3 | null {
		/*

		 this: n0 * x = w0
		 plane: n1 * x = w1
		 plane perpendicular to both which goes through origin:
		 n2 := n0 X x1
		 n2 * x = 0
		 */
		assertInst(P3, plane)
		assert(!this.isParallelToPlane(plane), '!this.isParallelToPlane(plane)')
		/*
		 var n0 = this.normal1, n1 = plane.normal1, n2 = n0.cross(n1).unit(), m = M4.forSys(n0, n1, n2)
		 var x0 = this.anchor, x1 = plane.anchor, x2 = V3.O
		 var p = n2.times(x2.dot(n2))
		 .plus(n1.cross(n2).times(x0.dot(n0)))
		 .plus(n2.cross(n0).times(x1.dot(n1)))
		 .div(m.determinant())
		 */
		const n0 = this.normal1, n1 = plane.normal1, n2 = n0.cross(n1).unit()
		const p = M4.forRows(n0, n1, n2).inversed().transformVector(new V3(this.w, plane.w, 0))
		return new L3(p, n2)
	}

	/**
	 * Returns the point in the plane closest to the given point
	 *
	 */
	projectedPoint(x: V3): V3 {
		// See http://math.stackexchange.com/questions/444968/project-a-point-in-3d-on-a-given-plane
		// p = x - ((x - planeAnchor) * normal1) * normal1
		return x.minus(this.normal1.times(x.minus(this.anchor).dot(this.normal1)))
	}

	projectedVector(x: V3): V3 {
		// See V3.rejectedFrom. Simplified, as this.normal1.length() == 1
		return x.minus(this.normal1.times(x.dot(this.normal1)))
	}

	flipped(): P3 {
		return new P3(this.normal1.negated(), -this.w)
	}

	static throughPoints(a: V3, b: V3, c: V3): P3 {
		assertVectors(a, b, c)
		const n1 = b.minus(a).cross(c.minus(a)).unit()
		return new P3(n1, n1.dot(a))
	}

	static normalOnAnchor(normal: V3, anchor: V3): P3 {
		assertVectors(normal, anchor)
		const n1 = normal.unit()
		return new P3(n1, n1.dot(anchor))
	}

	/**
	 * x/x0 + y/y0 + y/y0 = 1
	 *
	 */
	static forAxisIntercepts(x0: number, y0: number, z0: number): P3 {
		assertNumbers(x0, y0, z0)
		const normal = new V3(1 / x0, 1 / y0, 1 / z0)
		return new P3(normal.unit(), normal.length())
	}

	static forAnchorAndPlaneVectors(anchor: V3, v0: V3, v1: V3): P3 {
		assertVectors(anchor, v0, v1)
		return P3.normalOnAnchor(v0.cross(v1), anchor)
	}


	static readonly YZ = new P3(V3.X, 0)
	static readonly ZX = new P3(V3.Y, 0)
	static readonly XY = new P3(V3.Z, 0)

    containsCurve(curve: Curve) {
        if (curve instanceof L3) {
            return this.containsLine(curve)
        } else if (curve instanceof SemiEllipseCurve ||
            curve instanceof EllipseCurve ||
            curve instanceof HyperbolaCurve ||
            curve instanceof ParabolaCurve) {
            return this.containsPoint(curve.center) && this.normal1.isParallelTo(curve.normal)
        } else if (curve instanceof BezierCurve) {
            return curve.points.every(p => this.containsPoint(p))
        } else {
            throw new Error('' + curve)
        }

    }
}
NLA.registerClass(P3)