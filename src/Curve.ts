abstract class Curve extends Transformable {
	tIncrement: number
	hlol: number

	/**
	 * Returns curve parameter t for point p on curve.
	 */
	abstract pointLambda(p:V3, hint?): number

	/**
	 * Returns the point on the line that is closest to the given point.
	 */
	closestPointToPoint(p: V3): V3 {
		return this.at(this.closestTToPoint(p))
	}

	abstract closestTToPoint(p: V3): number

	/**
	 * So different edges on the same curve do not have different vertices, they are always generated
	 * on fixed points this.at(k * this.tIncrement), with k taking integer values
	 *
	 */
	calcSegmentPoints(aT: number, bT: number, a: V3, b: V3, reversed: boolean, includeFirst: boolean): V3[] {
		assert(this.tIncrement, "tIncrement not defined on " + this)
		var split = 4 * 62, inc = this.tIncrement
		var verts = []
		if (includeFirst) verts.push(a)
		if (!reversed) {
			assert(aT < bT)
			let start = Math.ceil((aT + NLA_PRECISION) / inc)
			let end = Math.floor((bT - NLA_PRECISION) / inc)
			for (let i = start; i <= end; i++) {
				verts.push(this.at(i * inc))
			}
		} else {
			assert(bT < aT)
			let start = Math.floor((aT - NLA_PRECISION) / inc)
			let end = Math.ceil((bT + NLA_PRECISION) / inc)
			for (let i = start; i >= end; i--) {
				verts.push(this.at(i * inc))
			}
		}
		verts.push(b)
		return verts
	}

	distanceToPoint(p: V3): number {
		return this.at(this.closestTToPoint(p)).distanceTo(p)
	}

	/**
	 * Behavior when curves are colinear: self intersections
	 */
	abstract isInfosWithCurve(curve): {tThis: number, tOther: number, p: V3}[]

	abstract at(t: number): V3

	abstract tangentAt(t: number): V3

	/**
	 * Derivative of tangentAt for parameter t
	 */
	abstract ddt(t: number): V3

	abstract containsPoint(p: V3): boolean

	abstract isTsWithSurface(surface: Surface): number[]

	abstract isTsWithPlane(plane: P3): number[]

	/**
	 * Should really be abstract, but it works for all the conic is curves, so it's here.
	 * @param mesh
	 * @param bufferName
	 */
	debugToMesh(mesh, bufferName) {
		mesh[bufferName] || mesh.addVertexBuffer(bufferName, bufferName)
		for (var t = -Math.PI; t < Math.PI; t += 0.1) {
			var p = this.at(t);
			mesh[bufferName].push(p, p.plus(this.tangentAt(t).toLength(1)))
			mesh[bufferName].push(p, p.plus(this.normalAt(t).toLength(1)))
		}
		mesh[bufferName].push(this.center, this.center.plus(this.f1.times(1.2)))
		mesh[bufferName].push(this.center, this.center.plus(this.f2))
		mesh[bufferName].push(this.center, this.center.plus(this.normal))
	}

	arcLength(startT, endT, steps) {
		assert(startT < endT, 'startT < endT')
		return gaussLegendreQuadrature24(t => this.tangentAt(t).length(), startT, endT)
		//return integrateCurve(this, startT, endT, steps || 1024)
	}

	/**
	 * iff for any t, this.at(t) == curve.at(t)
	 */
	abstract likeCurve(curve: Curve): boolean

	/**
	 * Return whether the curves occupy the same points in space. They do
	 * not necessarily need to share the same parameter values.
	 *
	 *
	 * iff for every t, there is an s so that this.at(t) == curve.at(s)
	 * and for every s, there is a t so that curve.at(s) == this.a(t)
	 */
	abstract isColinearTo(curve: Curve): boolean


	getAABB(tMin, tMax) {
		tMin = isFinite(tMin) ? tMin : this.tMin
		tMax = isFinite(tMax) ? tMax : this.tMax
		let tMinAt = this.at(tMin), tMaxAt = this.at(tMax)
		let roots = this.roots()
		let mins = new Array(3), maxs = new Array(3)
		for (let dim = 0; dim < 3; dim++) {
			let tRoots = roots[dim]
			mins[dim] = Math.min(tMinAt.e(dim), tMaxAt.e(dim))
			maxs[dim] = Math.max(tMinAt.e(dim), tMaxAt.e(dim))
			for (let j = 0; j < tRoots.length; j++) {
				let tRoot = tRoots[j]
				if (tMin < tRoot && tRoot < tMax) {
					mins[dim] = Math.min(mins[dim], this.at(tRoot).e(dim))
					maxs[dim] = Math.max(maxs[dim], this.at(tRoot).e(dim))
				}
			}
		}
		return new AABB(V(mins), V(maxs))
	}

	static hlol = 0
}