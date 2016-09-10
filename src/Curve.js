/**
 *
 */
class Curve extends NLA.Transformable {

	/**
	 * Returns curve parameter t for point p on curve.
	 * @abstract
	 * @param {V3} p
	 * @returns {number}
	 */
	pointLambda(p) {
		assert(false)
	}

	/**
	 * @abstract
	 * @param {M4} m4
	 * @returns {Curve}
	 */
	transform(m4) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 * So different edges on the same curve do not have different vertices, they are always generated
	 * on fixed points this.at(k * this.tIncrement), with k taking integer values
	 *
	 * @param {number} aT
	 * @param {number} bT
	 * @param {V3} a
	 * @param {V3} b
	 * @param {boolean} reversed
	 * @param {boolean} includeFirst
	 * @returns {Array.<V3>}
	 */
	calcSegmentPoints(aT, bT, a, b, reversed, includeFirst) {
		assert(this.tIncrement, "tIncrement not defined on " + this)
		var split = 4 * 62, inc = this.tIncrement
		var verts = []
		if (includeFirst) verts.push(a)
		if (!reversed) {
			assert(aT < bT)
			let start = Math.ceil((aT + NLA.PRECISION) / inc)
			let end = Math.floor((bT - NLA.PRECISION) / inc)
			console.log(aT, bT, start, end, inc)
			for (let i = start; i <= end; i++) {
				verts.push(this.at(i * inc))
			}
		} else {
			assert(bT < aT)
			let start = Math.floor((aT - NLA.PRECISION) / inc)
			let end = ceil((bT + NLA.PRECISION) / inc)
			for (let i = start; i >= end; i--) {
				verts.push(this.at(i * inc))
			}
		}
		verts.push(b)
		return verts
	}

	/**
	 *
	 * @abstract
	 * @param {number} t
	 * @returns {V3}
	 */
	at(t) {
		assert(false)
	}

	/**
	 *
	 * @abstract
	 * @param {number} t
	 * @returns {V3}
	 */
	tangentAt(t) {
		assert(false)
	}

	/**
	 * @abstract
	 * @param {V3} p
	 * @returns {boolean}
	 */
	containsPoint(p) {
		assert(false)
	}

	/**
	 * @abstract
	 * @param {Curve} curve
	 * @returns {boolean}
	 */
	likeCurve(curve) {
		assert(false)
	}

	/**
	 *
	 * @abstract
	 * @param {Surface} surface
	 * @returns {number[]}
	 */
	isTsWithSurface(surface) {
		assert(false)
	}

	/**
	 *
	 * @abstract
	 * @param {P3} plane
	 * @returns {number[]}
	 */
	isTsWithPlane(plane) {
		assert(false)
	}

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
	 * @abstract
	 * @param {Curve} curve
	 * @returns {boolean}
	 */
	isColinearTo(curve) {
		assert(false)
	}
}