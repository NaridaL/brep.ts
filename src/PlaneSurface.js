var PlaneSurface = class PlaneSurface extends Surface {
	/**
	 *
	 * @param {P3} plane
	 * @param {V3=} right
	 * @param {V3=} up
	 * @constructor
	 * @extends {Surface}
	 */
	constructor(plane, right, up) {
		super()
		assertInst(P3, plane)
		this.plane = plane
		this.up = up || plane.normal.getPerpendicular().normalized()
		this.right = right || this.up.cross(this.plane.normal).normalized()
		assert(this.right.cross(this.up).like(this.plane.normal))
	}

	isCoplanarTo(surface) {
		return surface instanceof PlaneSurface && this.plane.isCoplanarToPlane(surface.plane)
	}

	like(surface) {
		return surface instanceof PlaneSurface && this.plane.like(surface.plane)
	}

	parametricFunction() {
		var matrix = M4.forSys(this.right, this.up, this.normal, this.plane.anchor)
		return function (s, t) {
			return matrix.transformPoint(V3.create(s, t, 0))
		}
	}

	implicitFunction() {
		return p => this.plane.distanceToPointSigned(p)
	}

	isCurvesWithISurface(implicitSurface) {
		assert (implicitSurface.implicitFunction, 'implicitSurface.implicitFunction')
		return new CurvePI(this, implicitSurface)
	}

	/**
	 * @inheritDoc
	 */
	isCurvesWithSurface(surface2) {
		assert(false)
	}

	/**
	 *
	 * @param {Array.<Edge>} contour
	 * @returns {boolean}
	 */
	edgeLoopCCW(contour) {
		var totalAngle = 0
		for (var i = 0; i < contour.length; i++) {
			var ipp = (i + 1) % contour.length
			var edge = contour[i], nextEdge = contour[ipp]
			assert(edge.b.like(nextEdge.a), "edges dont form a loop")
			if (edge.curve instanceof EllipseCurve) {
				totalAngle += edge.rotViaPlane(this.plane.normal)
				// console.log(edge.toString(), edge.rotViaPlane(this.plane.normal))
			}
			totalAngle += edge.bDir.angleRelativeNormal(nextEdge.aDir, this.plane.normal)
		}
		return totalAngle > 0
	}

	/**
	 *
	 * @param {Edge[]} contour
	 * @param {V3} p
	 * @returns {boolean}
	 */
	edgeLoopContainsPoint(contour, p) {
		assert(contour)
		assertVectors (p)
		var dir = this.right.plus(this.up.times(0.123)).normalized()
		var line = L3(p, dir)
		var plane = this.plane
		var intersectionLinePerpendicular = dir.cross(plane.normal)
		var plane2 = P3.normalOnAnchor(intersectionLinePerpendicular, p)
		var colinearSegments = contour.map((edge) => edge.colinearToLine(line))
		var colinearSegmentsInside = contour.map((edge, i) => edge.aDir.dot(dir) > 0)
		var inside = false
		function logIS(p) {
			if (line.pointLambda(p) > 0) {
				inside = !inside
			}
		}
		contour.forEach((/** Edge= */ edge, i, edges) => {
			var j = (i + 1) % edges.length, nextEdge = edges[j]
			//console.log(edge.toSource()) {p:V3(2, -2.102, 0),
			if (colinearSegments[i]) {
				// edge colinear to intersection
				var outVector = edge.bDir.cross(plane.normal)
				var insideNext = outVector.dot(nextEdge.aDir) > 0
				if (colinearSegmentsInside[i] != insideNext) {
					logIS(edge.b)
				}
			} else {
				var edgeTs = edge.edgeISTsWithPlane(plane2)
				for (var k = 0; k < edgeTs.length; k++) {
					var edgeT = edgeTs[k]
					if (edgeT == edge.bT) {
						// endpoint lies on intersection line
						if (colinearSegments[j]) {
							// next segment is colinear
							// we need to calculate if the section of the plane intersection line BEFORE the colinear segment is
							// inside or outside the face. It is inside when the colinear segment out vector and the current segment vector
							// point in the same direction (dot > 0)
							var colinearSegmentOutsideVector = nextEdge.aDir.cross(plane.normal)
							var insideFaceBeforeColinear = colinearSegmentOutsideVector.dot(edge.bDir) < 0
							// if the "inside-ness" changes, add intersection point
							//console.log("segment end on line followed by colinear", insideFaceBeforeColinear != colinearSegmentInsideFace, nextSegmentOutsideVector)
							if (colinearSegmentsInside[j] != insideFaceBeforeColinear) {
								logIS(edge.b)
							}
						} else if (intersectionLinePerpendicular.dot(edge.bDir) * intersectionLinePerpendicular.dot(nextEdge.aDir) > 0) {
							logIS(edge.b)
						}
					} else if (edgeT != edge.aT) {
						// edge crosses line, neither starts nor ends on it
						logIS(edge.curve.at(edgeT))
					}
				}
			}
		})
		return inside

	}

	pointToParameterFunction() {
		var matrix = M4.forSys(this.right, this.up, this.normal, this.plane.anchor)
		var matrixInverse = matrix.inversed()
		return function (pWC) {
			return matrixInverse.transformPoint(pWC)
		}
	}

	normalAt(p) {
		return this.plane.normal
	}

	containsPoint(p) { return this.plane.containsPoint(p) }

	containsCurve(curve) {
		if (curve instanceof L3) {
			return this.plane.containsLine(curve)
		} else if (curve instanceof EllipseCurve) {
			return this.plane.containsPoint(curve.center) && this.plane.normal.isParallelTo(curve.normal)
		} else if (curve instanceof BezierCurve) {
			return curve.points.every(p => this.plane.containsPoint(p))
		} else {
			assert(false, edge.toString())
		}
	}

	transform(m4) {
		return new PlaneSurface(this.plane.transform(m4))
	}

	flipped() {
		return new PlaneSurface(this.plane.flipped(), this.right, this.up.negated())
	}

	toString() {
		return this.plane.toString()
	}

	toSource() {
		return `new PlaneSurface(${this.plane})`
	}
}
PlaneSurface.throughPoints = function (a, b, c) {
	return new PlaneSurface(P3.throughPoints(a, b, c))
}
