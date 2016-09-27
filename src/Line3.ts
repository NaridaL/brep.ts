/**
 * Created by aval on 20/11/2015.
 */
class L3 extends Curve {
	anchor:V3
	dir1:V3

	/**
	 * 3-dimensional line
	 *
	 * @param anchor
	 * @param dir1
	 * @class L3
	 * @constructor
	 * @property {V3} dir1 Normalized direction of the line.
	 * @property {V3} anchor Anchor of the line.
	 * @returns L3
	 * @augments Curve
	 */
	constructor(anchor, dir1) {
		super()
		assertVectors(anchor, dir1)
		assert(dir1.hasLength(1), "dir must be normalized" + dir1)
		assertf(() => !Number.isNaN(anchor.x))
		var l = Object.create(L3.prototype)
		l.dir1 = dir1
		l.anchor = anchor
		return l
	}

	containsPoint(point) {
		assertVectors(point)
		var dist = this.distanceToPoint(point);
		assertNumbers(dist)
		return NLA.isZero(dist)
	}

	likeCurve(curve) {
		return this.anchor.like(curve.anchor) && this.dir1.like(curve.dir1)
	}

	isColinearTo(obj) {
		return obj instanceof L3
			&& this.containsPoint(obj.anchor)
			&& NLA.equals(1, Math.abs(this.dir1.dot(obj.dir1)))
	}

	/**
	 *
	 * @param line
	 * @returns {number}
	 */
	distanceToLine(line: L3) {
		assertInst(L3, line)
		if (this.isParallelToLine(line)) {
			return this.distanceToPoint(line.anchor)
		}
		var dirCross1 = this.dir1.cross(line.dir1).normalized()
		var anchorDiff = this.anchor.minus(line.anchor)
		return Math.abs(anchorDiff.dot(dirCross1))
	}

	/**
	 *
	 * @param x
	 * @returns {number}
	 */
	distanceToPoint(x: V3) {
		assertVectors(x)
		// See http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
		var t = x.minus(this.anchor).dot(this.dir1)
		return this.at(t).distanceTo(x)

		//return x.minus(this.anchor).cross(x.minus(this.anchor.plus(this.dir1))).length()
	}

	asSegmentDistanceToPoint(x, sStart, sEnd) {
		var t = x.minus(this.anchor).dot(this.dir1)
		t = NLA.clamp(t, sStart, sEnd)
		return this.at(t).minus(x).length()
	}

	asSegmentDistanceToLine(line, sStart, sEnd) {
		assertInst(L3, line)
		var dirCross = this.dir1.cross(line.dir1)
		var div = dirCross.lengthSquared()
		if (NLA.isZero(div)) {
			return null
		} // lines parallel
		var anchorDiff = line.anchor.minus(this.anchor)
		// check if distance is zero (see also L3.distanceToLine)
		if (!NLA.isZero(anchorDiff.dot(dirCross.normalized()))) {
			return null
		}
		var t = this.infoClosestToLine(line).t
		t = NLA.clamp(t, sStart, sEnd)
		return this.at(NLA.clamp(t, sStart, sEnd))
	}

	/**
	 *
	 * @param lambda
	 * @returns {V3}
	 */
	at(lambda: number) {
		assertNumbers(lambda)
		return this.anchor.plus(this.dir1.times(lambda))
	}

	/**
	 * This function returns lambda for a given point x
	 *
	 * Every point x on this line is described by the equation
	 *      x = this.anchor + lambda * this.dir1 | - this.anchor
	 *      x - this.anchor = lambda * this.dir1 | DOT this.dir1
	 *      (x - this.anchor) DOT this.dir1 = lambda (dir1² is 1 as |dir1| == 1)
	 *
	 *  @param x
	 *  @returns {number}
	 */
	pointLambda(x: V3) {
		assertVectors(x)
		var t = x.minus(this.anchor).dot(this.dir1)
		return t
	}

	/**
	 * Returns true if the line is parallel (this.dir = line.dir || this.dir = -line.dir) to the argument.
	 *
	 * @param line
	 * @returns {boolean}
	 */
	isParallelToLine(line: L3) {
		assertInst(L3, line)
		// we know that 1 == this.dir1.length() == line.dir1.length(), we can check for parallelity simpler than isParallelTo()
		return NLA.equals(1, Math.abs(this.dir1.dot(line.dir1)))
	}

	/**
	 *
	 * @param line
	 * @returns {number}
	 */
	angleToLine(line: L3) {
		assertInst(L3, line)
		return this.dir1.angleTo(line.dir1)
	}

	/**
	 *
	 * @param line
	 * @returns {boolean} If the distance between the lines is zero
	 */
    intersectsLine(line) {
        return this.distanceToLine(line) <= NLA_PRECISION;
    }

	isInfosWithCurve(line) {
		assertInst(L3, line)

		let dirCross = this.dir1.cross(line.dir1)
		let div = dirCross.lengthSquared()
		if (NLA.isZero(div)) {
			// lines are parallel
			return []
		}
		let anchorDiff = line.anchor.minus(this.anchor)
		if (NLA.isZero(anchorDiff.dot(dirCross))) {
			let tThis = anchorDiff.cross(line.dir1).dot(dirCross) / div
			let tOther = anchorDiff.cross(this.dir1).dot(dirCross) / div
			let p = this.at(tThis)
			return [{tThis: tThis, tOther: tOther, p: p}]
		}
		return []
	}

	/**
	 *
	 * @param line
	 * @returns {V3}
	 */
	isInfoWithLine(line: L3) {
		assertInst(L3, line)
		var dirCross = this.dir1.cross(line.dir1)
		var div = dirCross.lengthSquared()
		if (NLA.isZero(div)) {
			return null
		} // lines parallel
		var anchorDiff = line.anchor.minus(this.anchor)
		// check if distance is zero (see also L3.distanceToLine)
		if (!NLA.isZero(anchorDiff.dot(dirCross.normalized()))) {
			return null
		}
		var t = anchorDiff.cross(line.dir1).dot(dirCross) / div
		return this.at(t)
	}

	/**
	 * returns s and t with this.at(s) == line.at(t)
	 *
	 * @param line
	 * @returns {{s: number, t: number}}
	 */
	intersectionLineST(line: L3) {
		// the two points on two lines the closest two each other are the ones whose
		// connecting
		// TODO Where does this come from?
		// TODO: return value when no IS?
		assertInst(L3, line)
		var dirCross = this.dir1.cross(line.dir1)
		var div = dirCross.lengthSquared()
		var anchorDiff = line.anchor.minus(this.anchor)
		var s = anchorDiff.cross(this.dir1).dot(dirCross) / div
		var t = anchorDiff.cross(line.dir1).dot(dirCross) / div
		return {s: s, t: t}
		//console.log(segmentIntersectsRay, a, b, "ab", ab, "p", p, "dir", dir, s > 0 && t / div >= 0 && t / div <= 1, "s", s, "t", t, "div", div)
	}

	ddt(t) {
		return V3.ZERO
	}

	toString(roundFunction?) {
		roundFunction = roundFunction || ((v) => +v.toFixed(4))
		return "new L3(" + this.anchor.toString(roundFunction) + ", " + this.dir1.toString(roundFunction) + ")"
	}

	closestTToPoint(p) {
		// similar logic as pointLambda; we project the vector (anchor -> p) onto dir1, then add anchor back to it
		let nearestT = p.minus(this.anchor).dot(this.dir1)
		return nearestT
	}

	infoClosestToLine(obj: L3):{t: number, s?: number, closest?: V3, closest2?: V3, distance: number} {
		/*
		 line = a + s*b
		 this = c + t*d

		 (this - line) * b = 0
		 (this - line) * d = 0

		 (a + s*b - c - t*d) * b = 0
		 (a + s*b - c - t*d) * d = 0

		 (a - c + s*b - t*d) * b = 0
		 (a - c + s*b - t*d) * d = 0

		 (a - c)*b + (s*b - t*d)*b = 0
		 (a - c)*d + (s*b - t*d)*d = 0

		 (a - c)*b + s*(b*b) - t*(d*b) = 0
		 (a - c)*d + s*(b*d) - t*(d*d) = 0

		 s = (t*(d*b) - (a - c)*b) / (b*b)
		 =>
		 (a - c)*d + (t*(d*b) - (a - c)*b) / (b*b)*(b*d) - t*(d*d) = 0 | * (b*b)
		 (a - c)*d * (b*b) + (t*(d*b) - (a - c)*b)*(b*d) - t*(d*d) * (b*b) = 0
		 (a - c)*d * (b*b) + t*(d*b)*(b*d) - (a - c)*b*(b*d) - t*(d*d) * (b*b) = 0
		 t = ((a - c)*b*(b*d) - (a - c)*d * (b*b)) / ((d*b)*(b*d) - (d*d) * (b*b))
		 */
		if (this.isParallelToLine(obj)) {
			return {t: NaN, s: NaN, distance: this.distanceToLine(obj)};
		}
		var a = obj.anchor, b = obj.dir1, c = this.anchor, d = this.dir1;
		var bd = b.dot(d), bb = b.lengthSquared(), dd = d.lengthSquared(), amc = a.minus(c), divisor = bd * bd - dd * bb;
		var t = (amc.dot(b) * bd - amc.dot(d) * bb) / divisor;
		var s = (amc.dot(b) * dd - amc.dot(d) * bd) / divisor;
		return {
			t: t,
			s: s,
			closest: this.at(t),
			closest2: obj.at(s),
			distance: this.at(t).distanceTo(obj.at(s))
		};
	}

	/**
	 *
	 * @param plane
	 * @returns {V3}
	 */
	intersectionWithPlane(plane: P3) {
		// plane: plane.normal * p = plane.w
		// line: p=line.point + lambda * line.dir1
		var lambda = (plane.w - plane.normal.dot(this.anchor)) / plane.normal.dot(this.dir1);
		var point = this.anchor.plus(this.dir1.times(lambda));
		return point;
	}

	tangentAt(t) {
		return this.dir1
	}

	intersectWithPlaneLambda(plane) {
		// plane: plane.normal * p = plane.w
		// line: p=line.point + lambda * line.dir1
		var div = plane.normal.dot(this.dir1)
		if (NLA.isZero(div)) return NaN
		var lambda = (plane.w - plane.normal.dot(this.anchor)) / div
		return lambda
	}

	isTsWithPlane(plane) {
		return [this.intersectWithPlaneLambda(plane)]
	}

	flipped() {
		return new L3(this.anchor, this.dir1.negated())
	}

	transform(m4: M4):this {
		var newAnchor = m4.transformPoint(this.anchor)
		var newDir = m4.transformVector(this.dir1)
		return new L3(newAnchor, newDir.normalized())
	}

	projectedOnPlane(plane) {
		assertInst(P3, plane)
		return new L3(plane.projectedPoint(this.anchor), plane.projectedVector(this.dir1).normalized())
	}

	debugToMesh(mesh, bufferName) {
		mesh[bufferName] || mesh.addVertexBuffer(bufferName, bufferName)
		mesh[bufferName].push(this.at(-1000), this.at(2000))
	}

	static throughPoints = (anchor: V3, b: V3):L3 => new L3(anchor, b.minus(anchor).normalized())
	static anchorDirection = (anchor: V3, dir: V3):L3 => new L3(anchor, dir.normalized())


	static pointLambdaNotNormalized(anchor, dir, x) {
		assertVectors(anchor, dir, x)
		return x.minus(anchor).dot(dir) / dir.lengthSquared()
	}

	/**
	 *
	 * @param p1
	 * @param p2
	 * @returns {V3}
	 */
	static fromPlanes(p1: P3, p2: P3) {
		assertInst(P3, p1, p2)
		var direction = p1.normal.cross(p2.normal);
		var l = direction.length();
		if (l < 1e-10) {
			throw new Error("Parallel planes");
		}

		return p1.intersectionWithPlane(p2)
	}

	static X:L3 = new L3(V3.ZERO, V3.X)
	static Y:L3 = new L3(V3.ZERO, V3.Y)
	static Z:L3 = new L3(V3.ZERO, V3.Z)

}
L3.prototype.hlol = Curve.hlol++
NLA.registerClass(L3)