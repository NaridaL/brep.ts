/**
 * Created by aval on 23/02/2016.
 */

/**
 * @extends {Curve.<BezierCurve>}
 */
class BezierCurve extends Curve {
	/**
	 *
	 * @param {V3} p0
	 * @param {V3} p1
	 * @param {V3} p2
	 * @param {V3} p3
	 */
	constructor(p0, p1, p2, p3) {
		super()
		assertVectors(p0, p1, p2, p3)
		this.p0 = p0
		this.p1 = p1
		this.p2 = p2
		this.p3 = p3
	}

	/**
	 *
	 * @returns {V3[]}
	 */
	get points() {
		return [this.p0, this.p1, this.p2, this.p3]
	}

	toString(f) {
		return `new BezierCurve(${this.p0}, ${this.p1}, ${this.p2}, ${this.p3})`
	}

	/**
	 * @inheritDoc
	 */
	at(t) {
		assertNumbers(t)
		var p0 = this.p0, p1 = this.p1, p2 = this.p2, p3 = this.p3
		var s = 1 - t, c0 = s * s * s, c1 = 3 * s * s * t, c2 = 3 * s * t * t, c3 = t * t * t
		return V3.create(
			p0.x * c0 + p1.x * c1 + p2.x * c2 + p3.x * c3,
			p0.y * c0 + p1.y * c1 + p2.y * c2 + p3.y * c3,
			p0.z * c0 + p1.z * c1 + p2.z * c2 + p3.z * c3
		)
	}

	tangentAt(t) {
		assertNumbers(t)
		var p0 = this.p0, p1 = this.p1, p2 = this.p2, p3 = this.p3
		var s = 1 - t, c01 = 3 * s * s, c12 = 6 * s * t, c23 = 3 * t * t
		return V3.create(
			(p1.x - p0.x) * c01 + (p2.x - p1.x) * c12 + (p3.x - p2.x) * c23,
			(p1.y - p0.y) * c01 + (p2.y - p1.y) * c12 + (p3.y - p2.y) * c23,
			(p1.z - p0.z) * c01 + (p2.z - p1.z) * c12 + (p3.z - p2.z) * c23
		)
	}

	ddt(t) {
		assertNumbers(t)
		var p0 = this.p0, p1 = this.p1, p2 = this.p2, p3 = this.p3
		var c012 = 6 * (1 - t), c123 = 6 * t
		return V3.create(
			(p2.x - 2 * p1.x + p0.x) * c012 + (p3.x - 2 * p2.x + p1.x) * c123,
			(p2.y - 2 * p1.y + p0.y) * c012 + (p3.y - 2 * p2.y + p1.y) * c123,
			(p2.z - 2 * p1.z + p0.z) * c012 + (p3.z - 2 * p2.z + p1.z) * c123
		)
	}

	normalAt(t) {
		var tangent = this.tangentAt(t)
		var rot = tangent.cross(this.ddt(t))
		return M4.rotation(Math.PI / 2, rot).transformVector(tangent)
	}


	/**
	 * @inheritDoc
	 */
	isTsWithPlane(plane) {
		assertInst(P3, plane)
		/*
		We are solving for t:
		 n := plane.normal
		 this.at(t) DOT n == plane.w // according to plane definition
		 (a t³ + b t² + c t + d) DOT n == plane.w // bezier curve as cubic equation
		 (a DOT n) t³ + (b DOT n) t³ + (c DOT n) t + d DOT n - plane.w == 0 // multiply out DOT n, minus plane.w
		 */

		let {p0, p1, p2, p3} = this
		let n = plane.normal
		let a = p1.minus(p2).times(3).minus(p0).plus(p3)
		let b = p0.plus(p2).times(3).minus(p1.times(6))
		let c = p1.minus(p0).times(3)
		let d = p0

		return solveCubicReal2(a.dot(n), b.dot(n), c.dot(n), d.dot(n) - plane.w)
	}

	/**
	 * @inheritDoc
	 */
	likeCurve(curve) {
		return curve.constructor == BezierCurve
			&& this.p0.like(curve.p0)
			&& this.p1.like(curve.p1)
			&& this.p2.like(curve.p2)
			&& this.p3.like(curve.p3)
	}

	/**
	 * Checks if this curve is colinear to the passed curve, i.e.
	 * for every t:number there exists a s:number with this.at(t) = curve.at(s)
	 * @param {BezierCurve} curve
	 * @returns {boolean}
	 */
	isColinearTo(curve) {
		// first, find out where/if curve.p0 and curve.p3 are on this
		// then split this at curve.p0 --> curve.p3 to compare points p1 and p2
		var curveP0T, curveP3T
		// assign in if condition to exploit short-circuit
		if (isNaN(curveP0T = this.pointLambda(curve.p0)) || isNaN(curveP3T = this.pointLambda(curve.p3))) {
			return false
		}
		var thisSplit
		if (NLA.equals(1, curveP0T)) {
			// this.split(curveP0T).right is degenerate in this case, so we need to handle it separately

			// this.split(curveP3T): 0 --> curveP3T --> 1
			// .right: curveP3T --> 1
			// .reversed(): 1 --> curveP3T
			thisSplit = this.split(curveP3T).right.reversed()
		} else {
			// curveP3T describes the point on this
			// adjust it so it describes the same point on this.split(curveP0T).right
			// this:                       0           p0t        p3t      1
			//                             |            |          |       |
			// this.split(curveP0T).right:              0        p3tad     1
			var curveP3Tadjusted = (curveP3T - curveP0T) / (1 - curveP0T)
			thisSplit = this.split(curveP0T).right.split(curveP3Tadjusted).left
		}

		return curve.likeCurve(thisSplit)
	}

	/**
	 * @returns {BezierCurve}
	 */
	reversed() {
		return new BezierCurve(this.p3, this.p2, this.p1, this.p0)
	}

	pointLambda(p) {
		var {p0, p1, p2, p3} = this
		// calculate cubic equation coefficients
		// a t³ + b t² + c t + d = 0
		// multiplying out the cubic Bézier curve equation gives:
		// a = -p0 + 3 p1 - 3 p2 + p3
		// b = 3 p0 - 6 p1 + 3 p2
		// c = -3 p0 + 3 p1
		// d = p0 - p
		var a = p1.minus(p2).times(3).minus(p0).plus(p3).els()
		var b = p0.plus(p2).times(3).minus(p1.times(6)).els()
		var c = p1.minus(p0).times(3).els()
		var d = p0.minus(p).els()
		var results = null

		// assume passed point is on curve and that curve does not self-intersect,
		// i.e. there is exactly one correct result for t
		// try to find a single result in the x-dimension, if multiple are found,
		// filter them by checking the other dimensions
		for (var dim = 0; dim < 3; dim++) {
			if (NLA.isZero(a[dim]) && NLA.isZero(b[dim]) && NLA.isZero(c[dim])) {
				// for case x:
				// ax == bx == cx == 0 => x(t) = dx
				// x value is constant
				// if x == 0 for all t, this does not limit the result, otherwise, there is no result, i.e
				// the passed point is not on the curve
				if (NLA.isZero(d[dim])) return NaN
			} else {
				var newResults = solveCubicReal2(a[dim], b[dim], c[dim], d[dim])
				if (0 == newResults.length) return NaN
				if (1 == newResults.length) return newResults[0]
				if (results) {
					results = results.filter(t => newResults.some(t2 => NLA.equals(t, t2)))
					if (results.length == 1) return results[0]
				} else {
					results = newResults
				}
			}
		}
		assert(false, 'multiple intersection ' + this.toString() + p.ss)
	}

	/**
	 * @inheritDoc
	 */
	transform(m4) {
		return new BezierCurve(
			m4.transformPoint(this.p0),
			m4.transformPoint(this.p1),
			m4.transformPoint(this.p2),
			m4.transformPoint(this.p3))
	}

	isClosed() {
		return this.p0.like(this.p3)
	}

	debugToMesh(mesh, bufferName) {
		mesh.addVertexBuffer(bufferName, bufferName)
		for (var t = -2; t <= 2; t += 0.01) {
			var p = this.at(t);
			mesh[bufferName].push(p, p.plus(this.tangentAt(t).toLength(1)))
			mesh[bufferName].push(p, p.plus(this.normalAt(t).toLength(1)))
		}
		mesh[bufferName].push(this.p0, this.p1)
		mesh[bufferName].push(this.p1, this.p2)
	}

	/**
	 *
	 * @param {number} t
	 * @returns {{left: BezierCurve, right: BezierCurve}}
	 */
	split(t) {
		var s = (1 - t)
		var {p0, p1, p2, p3} = this
		var b01 = p0.times(s).plus(p1.times(t)), b11 = p1.times(s).plus(p2.times(t)), b21 = p2.times(s).plus(p3.times(t))
		var b02 = b01.times(s).plus(b11.times(t)), b12 = b11.times(s).plus(b21.times(t))
		var b03 = b02.times(s).plus(b12.times(t))
		return {left: new BezierCurve(p0, b01, b02, b03), right: new BezierCurve(b03, b12, b21, p3)}
	}

	containsPoint(p) {
		return isFinite(this.pointLambda(p))
	}

	/**
	 *
	 * @param {V3} p
	 * @param {number=} tStart Defines interval with tEnd in which a start value for t will be searched.
	 * Result is not necessarily in this interval.
	 * @param {number=} tEnd
	 * @returns {number}
	 */
	closestTToPoint(p, tStart=0, tEnd=1) {
		let startT = NLA.arrayFromFunction(16, i => tStart + (tEnd - tStart) * i / 15)
			.withMax(t => -this.at(t).distanceTo(p))

		// this.at(t) has minimal distance to p when this.tangentAt(t) is perpendicular to
		// the vector between this.at(t) and p. This is the case iff the dot product of the two is 0.
		let dotTangentAndAtTToP = t => this.at(t).minus(p).dot(this.tangentAt(t))

		return newtonIterate(dotTangentAndAtTToP, startT, 8)
	}

	distanceToPoint(p, tStart, tEnd) {
		let t = this.closestTToPoint(p, tStart, tEnd)
		return this.at(t).distanceTo(p)
	}

	asSegmentDistanceToPoint(p, tStart, tEnd) {
		var t = this.closestTToPoint(p, tStart, tEnd)
		t = NLA.clamp(t, tStart, tEnd)
		return this.at(t).distanceTo(p)
	}

	/**
	 *
	 * @returns {number[]}
	 */
	isTsWithLine(/** V3 */ anchor, /** V3 */ dir1) {
		// looking for this.at(t) == at(s)
		// this.at(t).x == anchor.x + dir1.x * s
		// (this.at(t).x - anchor.x) / dir1.x == s (analogue for y and z)
		// (this.at(t).x - anchor.x) / dir1.x - (this.at(t).y - anchor.y) / dir1.y == 0
		// (this.at(t).x - anchor.x) * dir1.y - (this.at(t).y - anchor.y) * dir1.x == 0 (1)

		// cubic equation params (see #pointLambda):
		let {p0, p1, p2, p3} = this
		let a = p1.minus(p2).times(3).minus(p0).plus(p3)
		let b = p0.plus(p2).times(3).minus(p1.times(6))
		let c = p1.minus(p0).times(3)
		let d = p0

		// modifier cubic equation parameters to get (1)
		// let w = a.x * dir1.y - a.y * dir1.x
		// let x = b.x * dir1.y - b.y * dir1.x
		// let y = c.x * dir1.y - c.y * dir1.x
		// let z = (d.x - anchor.x) * dir1.y - (d.y - anchor.y) * dir1.x

		// the above version doesn't work for dir1.x == dir1.y == 0, so:
		let absMinDim = dir1.absMinDim()
		let [coord0, coord1] = [[1, 2], [2, 0], [0, 1]][absMinDim]

		let w = a.e(coord0) * dir1.e(coord1) - a.e(coord1) * dir1.e(coord0)
		let x = b.e(coord0) * dir1.e(coord1) - b.e(coord1) * dir1.e(coord0)
		let y = c.e(coord0) * dir1.e(coord1) - c.e(coord1) * dir1.e(coord0)
		let z = (d.e(coord0) - anchor.e(coord0)) * dir1.e(coord1) - (d.e(coord1) - anchor.e(coord1)) * dir1.e(coord0)

		return solveCubicReal2(w, x, y, z)
	}

	/**
	 *
	 * @param {L3} line
	 * @returns {V3[]}
	 */
	isPointsWithLine(line) {
		assertInst(L3, line)

		return this.isTsWithLine(line.anchor, line.dir1).map(t => this.at(t)).filter(p => line.containsPoint(p))
}

	/**
	 * Returns a curve with curve.at(x) == V3(x, ax³ + bx² + cx + d, 0)
	 * 
	 * @param a
	 * @param b
	 * @param c
	 * @param d
	 * @returns {BezierCurve}
	 */
	static graphXY(/** number */ a, /** number */ b, /** number */ c, /** number */ d) {
		// d = p0y
		// c = -3 p0y + 3 p1y => p1y = c/3 + p0y
		// b = 3 p0y - 6 p1y + 3 p2y => p2y = b/3 - p0y + 2 p1y
		// a = -p0y + 3 p1y -3 p2y + p3y => p3y = a + p0y - 3 p1y + 3 p2y
		let p0y = d
		let p1y = c / 3 + p0y
		let p2y = b / 3 - p0y + 2 * p1y
		let p3y = a + p0y - 3 * p1y + 3 * p2y
		return new BezierCurve(V3(0, p0y), V3(1 / 3, p1y), V3(2 / 3, p2y), V3(1, p3y))
	}
}
BezierCurve.prototype.tIncrement = 1 / (16)