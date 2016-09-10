/**
 * @class Surface
 * @template T
 * @augments NLA.Transformable<T>
 */
class Surface extends NLA.Transformable {
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
	 * @returns {Array.<V3>}
	 */
	isPointsWithLine(line) {
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
	 * @param {Curve|L3} curve
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
	 * @abstract
	 * @param {Surface} surface
	 * @returns {boolean}
	 */
	isCoplanarTo(surface) {
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
}