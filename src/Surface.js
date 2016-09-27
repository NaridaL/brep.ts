/**
 * @class Surface
 * @template T
 * @augments Transformable<T>
 */
class Surface extends Transformable {
	toString() {
		return this.toSource()
	}
	/**
	 * @abstract
	 */
	toSource() {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 *
	 * @abstract
	 * @param {L3} line
	 * @returns {number[]}
	 */
	isTsForLine(line) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 *
	 * @abstract
	 * @param {P3} plane
	 * @returns {Curve[]}
	 */
	isCurvesWithPlane(plane) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 * IMPORTANT: The tangents of the resulting curves need to be equal to the cross product of this and surface in the point.
	 * I.e.: for every point p p on a returned curve: curve.tangentAt(curve.pointLambda(p)) == this.normalAt(p) X surface.normalAt(p)
	 *
	 * Cross product is not commutative, so curve.tangentAt(curve.pointLambda(p)) == surface.normalAt(p) X this.normalAt(p)
	 * is not valid.
	 *
	 * @abstract
	 * @param {Surface} surface
	 * @returns {Curve[]}
	 */
	isCurvesWithSurface(surface) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 *
	 * @abstract
	 * @param {Curve} curve
	 * @returns {boolean}
	 */
	containsCurve(curve) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 *
	 * @abstract
	 * @param {V3} p
	 * @returns {boolean}
	 */
	containsPoint(p) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 * @return {T}
	 */
	flipped() {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 *
	 * @param {V3} p
	 * @returns {V3}
	 */
	normalAt(p) {
		var pmPoint = this.pointToParameterFunction()(p)
		return this.parametricNormal()(pmPoint.x, pmPoint.y)
	}

	/**
	 * @abstract
	 * @param {Edge[]} contour
	 * @param {V3} point
	 * @returns {boolean}
	 */
	edgeLoopContainsPoint(contour, point) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 * Returns true iff the surfac occupies the same space as the argument (not necessarily same normal)
	 * @abstract
	 * @param {Surface} surface
	 * @returns {boolean}
	 */
	isCoplanarTo(surface) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 * @abstract
	 * @returns {boolean}
	 */
	like(object) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 * @abstract
	 * @param {M4} m4
	 * @returns {T}
	 */
	transform(m4) {
		assert(false, "Not implemented on " + this.constructor.name)
	}

	/**
	 * @abstract
	 * @returns {function (V3):V3}
	 */
	pointToParameterFunction() {
		assert(false, "Not implemented on " + this.constructor.name)
	}
}