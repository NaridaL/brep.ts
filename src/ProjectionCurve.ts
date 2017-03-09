class ProjectionCurve extends Curve {
	constructor(originalCurve, projDir, targetSurface, isIndex) {
		super()
		this.originalCurve = originalCurve
		this.projDir = projDir
		this.targetSurface = targetSurface
		this.isIndex = isIndex
	}

	/**
	 * Returns curve parameter t for point p on curve.
	 * @param p This is a test
	 */
	pointT(p) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 * Returns the point on the line that is closest to the given point.
	 *
	 * @param p
	 * @returns {V3}
	 */
	closestPointToPoint(p: V3) {
		return this.at(this.closestTToPoint(p))
	}

	/**
	 * @abstract
	 * @param p
	 */
	closestTToPoint(p: V3) {
	}

	/**
	 *
	 * @param p
	 * @returns {number}
	 */
	distanceToPoint(p: V3) {
		return this.at(this.closestTToPoint(p)).distanceTo(p)
	}

	/**
	 * Behavior when curves are colinear: self intersections
	 * @abstract
	 * @param curve
	 * @returns {{tThis:number, tOther:number, p:V3}[]}
	 */
	isInfosWithCurve(curve) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 *
	 * @abstract
	 * @param t
	 * @returns {V3}
	 */
	at(t: number) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 *
	 * @abstract
	 * @param t
	 * @returns {V3}
	 */
	tangentAt(t: number) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 * Derivative of tangentAt(t) for parameter t
	 * @abstract
	 * @param t
	 * @returns {V3}
	 */
	ddt(t: number) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 * @abstract
	 * @param p
	 * @returns {boolean}
	 */
	containsPoint(p: V3) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 *
	 * @abstract
	 * @param surface
	 * @returns {number[]}
	 */
	isTsWithSurface(surface: Surface) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 *
	 * @abstract
	 * @param plane
	 * @returns {number[]}
	 */
	isTsWithPlane(plane: P3) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 * Should really be abstract, but it works for all teh conic is curves, so it's here.
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
	 *
	 * iff for any t, this.at(t) == curve.at(t)
	 *
	 * @abstract
	 * @param curve
	 * @returns {boolean}
	 */
	likeCurve(curve: Curve) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 * Return whether the curves occupy the same points in space. They do
	 * not necessarily need to share the same parameter values.
	 *
	 *
	 * iff for every t, there is an s so that this.at(t) == curve.at(s)
	 * and for every s, there is a t so that curve.at(s) == this.a(t)
	 *
	 * @abstract
	 * @param curve
	 * @returns {boolean}
	 */
	isColinearTo(curve: Curve) {
		assert(false, "Not implemented on " + this.constructor.name)
	}


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
}