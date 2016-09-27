/**
 * Created by aval on 23/02/2016.
 */

/**
 * Blah
 *
 * @extends {Curve.<BezierCurve>}
 */
class BezierCurve extends Curve {
	/**
	 *
	 * @param {V3} p0
	 * @param {V3} p1
	 * @param {V3} p2
	 * @param {V3} p3
	 * @param {number=} tMin
	 * @param {number=} tMax
	 */
	constructor(p0, p1, p2, p3, tMin, tMax) {
		super()
		assertVectors(p0, p1, p2, p3)
		this.p0 = p0
		this.p1 = p1
		this.p2 = p2
		this.p3 = p3
		this.tMin = isFinite(tMin) ? tMin : -0.1
		this.tMax = isFinite(tMax) ? tMax : 1.1
	}

	/**
	 *
	 * @returns {V3[]}
	 */
	get points() {
		return [this.p0, this.p1, this.p2, this.p3]
	}

	toString(f) {
		return `new BezierCurve(${this.p0}, ${this.p1}, ${this.p2}, ${this.p3}, ${this.tMin}, ${this.tMax})`
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

	/**
	 * s := (1 - t)
	 * at(t) := s³ p0 + 3 s² t p1 + 3 s t² p2 + t³ p3
	 * tangent(t) := 3 s² (p1 - p0) + 6 s t (p2 - p1) + 3 t² (p3 - p2)
	 *            := 3 (1 - t)² (p1 - p0) + 6 (1 - t) t (p2 - p1) + 3 t² (p3 - p2)
	 *            := 3 (1 - 2 t + t²) (p1 - p0) + 6 (t - t²) (p2 - p1) + 3 t² (p3 - p2)
	 *            := (3 (p3 - p2) - 6 (p2 - p1) + 3 (p1 - p0)) t²*
	 *                + (-6 (p1 - p0) + (p2 - p1)) t
	 *                + 3 (p1 - p0)
	 *
	 *
	 * @param t
	 * @returns {V3}
	 */
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
		return rot.cross(tangent)
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

	isTsWithSurface(surface) {
		if (surface instanceof PlaneSurface) {
			return this.isTsWithPlane(surface.plane)
		}
		if (surface instanceof CylinderSurface) {
			let projPlane = P3(surface.dir.normalized(), 0)
			let projThis = this.project(projPlane)
			let projEllipse = surface.baseEllipse.project(projPlane)
			return projEllipse.isInfosWithBezier2D(projThis).map(info => info.tOther)
		}
		assert(false)
	}

	/**
	 * @inheritDoc
	 */
	likeCurve(curve) {
		return this == curve
			||curve.constructor == BezierCurve
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
		if (curve == this) return true
		if (!(curve instanceof BezierCurve)) return false
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
			thisSplit = this.split(curveP3T).right.reversed4()
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
		return new BezierCurve(this.p3, this.p2, this.p1, this.p0, -this.tMax, -this.tMin)
	}
	
	getCoefficients() {
		let {p0, p1, p2, p3} = this
		// calculate cubic equation coefficients
		// a t³ + b t² + c t + d = 0
		// multiplying out the cubic Bézier curve equation gives:
		// a = -p0 + 3 p1 - 3 p2 + p3
		// b = 3 p0 - 6 p1 + 3 p2
		// c = -3 p0 + 3 p1
		// d = p0 - p
		let a = p1.minus(p2).times(3).minus(p0).plus(p3)
		let b = p0.plus(p2).times(3).minus(p1.times(6))
		let c = p1.minus(p0).times(3)
		let d = p0
		return [a, b, c, d]
	}
	
	tangentCoefficients() {
		let {p0, p1, p2, p3} = this
		let p01 = p1.minus(p0), p12 = p2.minus(p1), p23 = p3.minus(p2)
		let a = p01.plus(p23).times(3).minus(p12.times(6))
		let b = p12.minus(p01).times(6)
		let c = p01.times(3)
		return [V3.ZERO, a, b, c]
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
		var a = p1.minus(p2).times(3).minus(p0).plus(p3)
		var b = p0.plus(p2).times(3).minus(p1.times(6))
		var c = p1.minus(p0).times(3)
		var d = p0.minus(p)

		let absMinDim = a.absMinDim()
		let [coord0, coord1] = [[1, 2], [2, 0], [0, 1]][absMinDim]


		let matrix = M4.forSys(a, b, c, d) // B(t) = matrix * V3(1, t, t², t³)
		const g = matrix.gauss().U.m
		let pqp = g[6] / g[5], q = g[7] / g[5]
		let results = solveCubicReal2(0, g[5], g[6], g[7]).filter(t => NLA.isZero(t * t * g[1] + t * g[2] + g[3] + g[0] * t * t * t))
		if (0 == results.length) return NaN
		if (1 == results.length) return results[0]
		assert(false, 'multiple intersection ' + this.toString() + p.sce)
	}

	pointLambda2(p) {
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
				if (!NLA.isZero(d[dim])) return NaN
			} else {

				var newResults = solveCubicReal2(a[dim], b[dim], c[dim], d[dim])
				console.log([a[dim], b[dim], c[dim], d[dim]], dim, newResults)
				if (0 == newResults.length) return NaN
				if (1 == newResults.length) return newResults[0]
				if (results) {
					results = results.filter(t => newResults.some(t2 => NLA.equals(t, t2)))
					if (0 == results.length) return NaN
					if (1 == results.length) return results[0]
				} else {
					results = newResults
				}
			}
		}
		assert(false, 'multiple intersection ' + results + this.toString() + p.sce)
	}

	/**
	 * @inheritDoc
	 * @returns BezierCurve
	 */
	transform(m4) {
		return new BezierCurve(
			m4.transformPoint(this.p0),
			m4.transformPoint(this.p1),
			m4.transformPoint(this.p2),
			m4.transformPoint(this.p3),
			this.tMin, this.tMax)
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
		mesh[bufferName].push(this.p2, this.p3)
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

	roots() {
		/**
		 *            := (3 (p3 - p2) - 6 (p2 - p1) + 3 (p1 - p0)) t²*
		 *                + (-6 (p1 - p0) + 6 (p2 - p1)) t
		 *                + 3 (p1 - p0)
		 *                */
		let {p0, p1, p2, p3} = this
		let p01 = p1.minus(p0), p12 = p2.minus(p1), p23 = p3.minus(p2)
		var a = p01.plus(p23).times(3).minus(p12.times(6))
		var b = p12.minus(p01).times(6)
		var c = p01.times(3)

		return NLA.arrayFromFunction(3, dim => solveCubicReal2(0, a.e(dim), b.e(dim), c.e(dim)))
	}


	/**
	 *
	 * @param {V3} p
	 * @param {number=} tMin Defines interval with tMax in which a start value for t will be searched.
	 * Result is not necessarily in this interval.
	 * @param {number=} tMax
	 * @returns {number}
	 */
	closestTToPoint(p, tMin, tMax, tStart) {
		tMin = isFinite(tMin) ? tMin : this.tMin
		tMax = isFinite(tMax) ? tMax : this.tMax


		// this.at(t) has minimal distance to p when this.tangentAt(t) is perpendicular to
		// the vector between this.at(t) and p. This is the case iff the dot product of the two is 0.
		// f = (this.at(t) - p) . (this.tangentAt(t)
		// df = this.tangentAt(t) . this.tangentAt(t) + (this.at(t) - p) . this.ddt(t)
		//    = this.tangentAt(t)² + (this.at(t) - p) . this.ddt(t)
		let f = t => this.at(t).minus(p).dot(this.tangentAt(t)) // 5th degree polynomial
		let df = t => this.tangentAt(t).lengthSquared() + (this.at(t).minus(p).dot(this.ddt(t)))

		const STEPS = 32
		let startT = 'undefined' === typeof tStart
			? NLA.arrayFromFunction(STEPS, i => tMin + (tMax - tMin) * i / STEPS)
					.withMax(t => -this.at(t).distanceTo(p))
			: tStart

		return newtonIterateWithDerivative(f, startT, 16, df)
	}

	/**
	 *
	 * @param {V3} p
	 * @param {number=} tStart Defines interval with tEnd in which a start value for t will be searched.
	 * Result is not necessarily in this interval.
	 * @param {number=} tEnd
	 * @returns {number}
	 */
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
	 * @param {V3} anchor
	 * @param {V3} dir
	 * @returns {{tThis: number, tOther: number, p: V3}[]}
	 */
	isInfosWithLine(anchor, dir, tMin, tMax) {
		// looking for this.at(t) == line.at(s)
		// this.at(t).x == anchor.x + dir.x * s
		// (this.at(t).x - anchor.x) / dir.x == s (analogue for y and z) (1x, 1y, 1z)
		// (1x) - (1y):
		// (this.at(t).x - anchor.x) / dir.x - (this.at(t).y - anchor.y) / dir.y == 0
		// (this.at(t).x - anchor.x) * dir.y - (this.at(t).y - anchor.y) * dir.x == 0 (2)

		// cubic equation params (see #pointLambda):
		let {p0, p1, p2, p3} = this
		let a = p1.minus(p2).times(3).minus(p0).plus(p3)
		let b = p0.plus(p2).times(3).minus(p1.times(6))
		let c = p1.minus(p0).times(3)
		let d = p0

		// modifier cubic equation parameters to get (1)
		// let w = a.x * dir.y - a.y * dir.x
		// let x = b.x * dir.y - b.y * dir.x
		// let y = c.x * dir.y - c.y * dir.x
		// let z = (d.x - anchor.x) * dir.y - (d.y - anchor.y) * dir.x

		// the above version doesn't work for dir.x == dir.y == 0, so:
		let absMinDim = dir.absMinDim()
		let [coord0, coord1] = [[1, 2], [2, 0], [0, 1]][absMinDim]

		let w = a.e(coord0) * dir.e(coord1) - a.e(coord1) * dir.e(coord0)
		let x = b.e(coord0) * dir.e(coord1) - b.e(coord1) * dir.e(coord0)
		let y = c.e(coord0) * dir.e(coord1) - c.e(coord1) * dir.e(coord0)
		let z = (d.e(coord0) - anchor.e(coord0)) * dir.e(coord1) - (d.e(coord1) - anchor.e(coord1)) * dir.e(coord0)

		tMin = isFinite(tMin) ? tMin : this.tMin
		tMax = isFinite(tMax) ? tMax : this.tMax

		// we ignored a dimension in the previous step, so we need to check it too
		return solveCubicReal2(w, x, y, z).mapFilter(tThis => {
			if (tMin <= tThis && tThis <= tMax) {
				let p = this.at(tThis)
				// console.log(t*t*t*w+t*t*x+t*y+z, dir.length())
				let s = p.minus(anchor).dot(dir) / dir.dot(dir)
				let lineAtS = dir.times(s).plus(anchor)
				if (lineAtS.like(p)) return {tThis: tThis, tOther: s, p: p}
			}
		})
	}

	closestPointToLine(line, tMin, tMax) {
		// (this(t)-line(s)) * line.dir == 0 (1)
		// (this(t)-line(s)) * this.tangentAt(t) == 0 (2)
		// this(t) * line.dir - line(s) * line.dir == 0
		// this(t) * line.dir - line.anchor * line.dir - s line.dir * line.dir == 0
		// this(t) * line.dir - line.anchor * line.dir == s (3)
		// insert (3) in (2)
		// (this(t)-line(this(t) * line.dir - line.anchor * line.dir)) * this.tangentAt(t) == 0 (4)
		// (4) is a 5th degree polynomial, solve numerically

		tMin = isFinite(tMin) ? tMin : this.tMin
		tMax = isFinite(tMax) ? tMax : this.tMax

		let anchorDotDir1 = line.anchor.dot(line.dir1)
		let f = t => {
			let atT = this.at(t);
			return (atT.minus(line.at(atT.dot(line.dir1) - anchorDotDir1))).dot(this.tangentAt(t)) }

		const STEPS = 32
		let startT = NLA.arrayFromFunction(STEPS, i => tMin + (tMax - tMin) * i / STEPS).withMax(t => -f(t))

		return newtonIterateWithDerivative(f, startT, 8, df)


	}

	/**
	 *
	 * @param {BezierCurve} bezier
	 * @param {number=} tMin
	 * @param {number=} tMax
	 * @param {number=} sMin
	 * @param {number=} sMax
	 * @returns {{tThis: number, tOther: number, p: V3}[]}
	 */
	isInfosWithBezie3(bezier, tMin, tMax, sMin, sMax) {
		const handleStartTS = (startT, startS) => {
			if (!result.some(info => NLA.equals(info.tThis, startT) && NLA.equals(info.tOther, startS))) {
				let f1 = (t, s) => this.tangentAt(t).dot(this.at(t).minus(bezier.at(s)))
				let f2 = (t, s) => bezier.tangentAt(s).dot(this.at(t).minus(bezier.at(s)))
				// f = (b1, b2, t1, t2) = b1.tangentAt(t1).dot(b1.at(t1).minus(b2.at(t2)))
				let fdt1 = (b1, b2, t1, t2) => b1.ddt(t1).dot(b1.at(t1).minus(b2.at(t2))) + (b1.tangentAt(t1).lengthSquared())
				let fdt2 = (b1, b2, t1, t2) => -b1.tangentAt(t1).dot(b2.tangentAt(t2))
				let ni = newtonIterate2dWithDerivatives(f1, f2, startT, startS, 16,
					fdt1.bind(undefined, this, bezier), fdt2.bind(undefined, this, bezier),
					(t, s) => -fdt2(bezier, this, s, t), (t, s) => -fdt1(bezier, this, s, t))
				result.push({tThis: ni.x, tOther: ni.y, p: this.at(ni.x)})
			}
		}

		tMin = isFinite(tMin) ? tMin : this.tMin
		tMax = isFinite(tMax) ? tMax : this.tMax
		sMin = isFinite(sMin) ? sMin : bezier.tMin
		sMax = isFinite(sMax) ? sMax : bezier.tMax

		// stack of indices:
		let indices = [tMin, tMax, sMin, sMax]
		let tMid = (tMin + tMax) / 2
		let sMid = (sMin + sMax) / 2
		let aabbs = [this.getAABB(tMin, tMdid), this.getAABB(tMid, tMax), bezier.getAABB(sMin, sMin), bezier.getAABB(sMid, sMax)]
		let result = []
		while (indices.length) {
			let i = indices.length - 4
			let tMin = indices[i], tMax = indices[i + 1], sMin = indices[i + 2], sMax = indices[i + 3]
			indices.length -= 4
			let thisAABB = this.getAABB(tMin, tMax)
			let otherAABB = bezier.getAABB(sMin, sMax)
			// console.log(tMin, tMax, sMin, sMax, thisAABB.sce, otherAABB.sce)
			if (thisAABB && otherAABB && thisAABB.intersectsAABB2d(otherAABB)) {
				let tMid = (tMin + tMax) / 2
				let sMid = (sMin + sMax) / 2
				const EPS = 0.00001
				if (tMax - tMin < EPS || sMax - sMin < EPS) {
					console.log(tMin, tMax, sMin, sMax)
					console.log(thisAABB.sce)
					console.log(otherAABB.sce)
					console.log(tMid, sMid)
					handleStartTS(tMid, sMid)
				} else {
					Array.prototype.push.call(indices,
						tMin, tMid, sMin, sMid,
						tMin, tMid, sMid, sMax,
						tMid, tMax, sMin, sMid,
						tMid, tMax, sMid, sMax)
				}
			}
		}

		return result
	}
	/**
	 *
	 * @param {BezierCurve} bezier
	 * @param {number=} tMin
	 * @param {number=} tMax
	 * @param {number=} sMin
	 * @param {number=} sMax
	 * @returns {{tThis: number, tOther: number, p: V3}[]}
	 */
	isInfosWithBezier(bezier, tMin, tMax, sMin, sMax) {
		// the recursive function finds good approximates for the intersection points
		// this function uses newton iteration to improve the result as much as possible
		// is declared as an arrow function so this will be bound correctly
		const handleStartTS = (startT, startS) => {
			if (!result.some(info => NLA.equals(info.tThis, startT) && NLA.equals(info.tOther, startS))) {
				let f1 = (t, s) => this.tangentAt(t).dot(this.at(t).minus(bezier.at(s)))
				let f2 = (t, s) => bezier.tangentAt(s).dot(this.at(t).minus(bezier.at(s)))
				// f = (b1, b2, t1, t2) = b1.tangentAt(t1).dot(b1.at(t1).minus(b2.at(t2)))
				let dfdt1 = (b1, b2, t1, t2) => b1.ddt(t1).dot(b1.at(t1).minus(b2.at(t2))) + (b1.tangentAt(t1).lengthSquared())
				let dfdt2 = (b1, b2, t1, t2) => -b1.tangentAt(t1).dot(b2.tangentAt(t2))
				let ni = newtonIterate2dWithDerivatives(f1, f2, startT, startS, 16,
					dfdt1.bind(undefined, this, bezier), dfdt2.bind(undefined, this, bezier),
					(t, s) => -dfdt2(bezier, this, s, t), (t, s) => -dfdt1(bezier, this, s, t))
				if (ni == null) console.log(startT, startS, this.sce, bezier.sce)
				result.push({tThis: ni.x, tOther: ni.y, p: this.at(ni.x)})
			}
		}

		// is declared as an arrow function so this will be bound correctly
		// returns whether an intersection was immediately found (i.e. without further recursion)
		// is declared as an arrow function so this will be bound correctly
		const findRecursive = (tMin, tMax, sMin, sMax, thisAABB, otherAABB) => {
			const EPS = NLA_PRECISION
			if (thisAABB.touchesAABB(otherAABB)) {
				let tMid = (tMin + tMax) / 2
				let sMid = (sMin + sMax) / 2
				if (Math.abs(tMax - tMin) < EPS && Math.abs(sMax - sMin) < EPS) {
					// console.log("AAAAAABBSSS", thisAABB, otherAABB, tMax, tMin, sMax, sMin)
					handleStartTS(tMid, sMid)
					return true
				} else {
					let thisAABBleft = this.getAABB(tMin, tMid), thisAABBright
					let bezierAABBleft = bezier.getAABB(sMin, sMid), bezierAABBright
					// if one of the following calls immediately finds an intersection, we don't want to call the others
					// as that will lead to the same intersection being output multiple times
					findRecursive(tMin, tMid, sMin, sMid, thisAABBleft, bezierAABBleft)
						|| findRecursive(tMin, tMid, sMid, sMax, thisAABBleft, bezierAABBright = bezier.getAABB(sMid, sMax))
						|| findRecursive(tMid, tMax, sMin, sMid, thisAABBright = this.getAABB(tMid, tMax), bezierAABBleft)
						|| findRecursive(tMid, tMax, sMid, sMax, thisAABBright, bezierAABBright)
					return false
				}
			}
			return false
		}

		tMin = isFinite(tMin) ? tMin : this.tMin
		tMax = isFinite(tMax) ? tMax : this.tMax
		sMin = isFinite(sMin) ? sMin : bezier.tMin
		sMax = isFinite(sMax) ? sMax : bezier.tMax
		assertf(() => tMin < tMax)
		assertf(() => sMin < sMax)
		let result = []

		let likeCurves = this.likeCurve(bezier), colinearCurves = this.isColinearTo(bezier)
		if (likeCurves || colinearCurves) {
			if (!likeCurves) {
				// only colinear
				// recalculate sMin and sMax so they are valid on this, from then on we can ignore bezier
				sMin = this.pointLambda(bezier.at(sMin))
				sMax = this.pointLambda(bezier.at(sMax))
			}
			tMin = Math.min(tMin, sMin)
			tMax = Math.max(tMax, sMax)
			let splits = NLA.fuzzyUniques(this.roots().concatenated().filter(isFinite).concat([tMin, tMax])).sort(NLA.minus)
			console.log('splits', splits, this.roots().concatenated())
			let aabbs = NLA.arrayFromFunction(splits.length - 1, i => this.getAABB(splits[i], splits[i + 1]))
			Array.from(NLA.combinations(splits.length - 1)).forEach(({i, j}) => {
				// adjacent curves can't intersect
				if (Math.abs(i - j) > 2) {
					// console.log(splits[i], splits[i + 1], splits[j], splits[j + 1], aabbs[i], aabbs[j])
					findRecursive(splits[i], splits[i + 1], splits[j], splits[j + 1], aabbs[i], aabbs[j])
				}
			})
		} else {
			findRecursive(tMin, tMax, sMin, sMax, this.getAABB(tMin, tMax), bezier.getAABB(sMin, sMax))
		}

		return result
	}

	selfIntersectionsInfo() {
		return this.isInfosWithBezier(this)
	}

	/** @inheritDoc */
	isInfosWithCurve(curve) {
		if (curve instanceof L3) {
			return this.isInfosWithLine(curve.anchor, curve.dir1)
		}
		if (curve instanceof BezierCurve) {
			return this.isInfosWithBezier(curve)
		}
		assert(false)
	}

	/**
	 * Returns a curve with curve.at(x) == V3(x, ax³ + bx² + cx + d, 0)
	 * 
	 * @param a
	 * @param b
	 * @param c
	 * @param d
	 * @param tMin
	 * @param tMax
	 * @returns {BezierCurve}
	 */
	static graphXY(/** number */ a, /** number */ b, /** number */ c, /** number */ d, tMin, tMax) {
		// d = p0y
		// c = -3 p0y + 3 p1y => p1y = c/3 + p0y
		// b = 3 p0y - 6 p1y + 3 p2y => p2y = b/3 - p0y + 2 p1y
		// a = -p0y + 3 p1y -3 p2y + p3y => p3y = a + p0y - 3 p1y + 3 p2y
		let p0y = d
		let p1y = c / 3 + p0y
		let p2y = b / 3 - p0y + 2 * p1y
		let p3y = a + p0y - 3 * p1y + 3 * p2y
		return new BezierCurve(V3(0, p0y), V3(1 / 3, p1y), V3(2 / 3, p2y), V3(1, p3y), tMin, tMax)
	}
}
/**
 * https://en.wikipedia.org/wiki/Cubic_function#/media/File:Graph_of_cubic_polynomial.svg
 * @type {BezierCurve}
 */
BezierCurve.EX2D = BezierCurve.graphXY(2,-3,-3,2)
BezierCurve.EX3D = new BezierCurve(V3.ZERO, V3(-0.1, -1, 1), V3(1.1, 1, 1), V3.X)
BezierCurve.prototype.hlol = Curve.hlol++
BezierCurve.prototype.tIncrement = 1 / (64)