abstract class Surface extends Transformable {
	toString():string {
		return this.toSource()
	}

	abstract toSource():string

	abstract isTsForLine(line:L3): number[]

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
	abstract isCurvesWithPlane(plane:P3): Curve[]

	abstract isCurvesWithSurface(surface: Surface):Curve[]

	abstract containsCurve(curve: Curve): boolean

	abstract containsPoint(p:V3):boolean

	abstract flipped():Surface

	/**
	 *
	 * @param p
	 * @returns {V3}
	 */
	normalAt(p:V3) {
		var pmPoint = this.pointToParameterFunction()(p)
		return this.parametricNormal()(pmPoint.x, pmPoint.y)
	}

	abstract edgeLoopContainsPoint(contour:Edge[], point:V3):boolean

	/**
	 * Returns true iff the surface occupies the same space as the argument (not necessarily same normal)
	 */
	abstract isCoplanarTo(surface:Surface): boolean

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