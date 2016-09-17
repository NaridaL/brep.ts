/**
 * Created by aval on 30/08/2016.
 *
 * Surface normal is (t, z) => this.baseCurve.tangentAt(t) X this.dir1
 * Choose dir1 appropriately to select surface orientation.
 */
class ProjectedCurveSurface extends Surface {

	constructor(baseCurve, dir1, tMin=0, tMax=1) {
		super()
		assertInst(Curve, baseCurve)
		assertInst(V3, dir1)
		assertNumbers(tMin, tMax)
		assert(dir1.hasLength(1), "dir must be normalized" + dir1)
		this.baseCurve = baseCurve
		this.dir1 = dir1
		this.tMin = tMin
		this.tMax = tMax
	}

	toSource() {
		return `new ProjectedCurveSurface(${this.baseCurve}, ${this.dir1}, ${this.tMin}, ${this.tMax})`
	}

	toString() {
		return this.toSource()
	}

	toMesh(tStart, tEnd, zStart, zEnd, count) {
		tStart = tStart || this.tMin
		tEnd = tEnd || this.tMax
		zStart = zStart || -30
		zEnd = zEnd || 30
		count = count || 16

		let tInterval = tEnd - tStart, tStep = tInterval / (count - 1)
		let baseVertices =
			NLA.arrayFromFunction(count, i => this.baseCurve.at(tStart + i * tStep).plus(this.dir1.times(zStart)))
		let normalFunc = this.parametricNormal()
		let normals = NLA.arrayFromFunction(count, i => normalFunc(tStart + i * tStep, 0))

		return GL.Mesh.offsetVertices(baseVertices, this.dir1.times(zEnd - zStart), false, normals)
	}

	parametricFunction() {
		let baseCurve = this.baseCurve, dir1 = this.dir1
		return (
			/**
			 *
			 * @param {number} t curve Parameter
			 * @param {number} z projection factor
			 * @returns {V3}
			 */
			function (t, z) {
				return baseCurve.at(t).plus(dir1.times(z))
			})
	}

	parametricNormal() {
		let baseCurve = this.baseCurve, dir1 = this.dir1
		return (t, z) => {
			return baseCurve.tangentAt(t).cross(dir1).normalized()
		}
	}

	implicitFunction() {
		var z = this.l3Axis.dir1, x = z.getPerpendicular().normalized(), y = z.cross(x)
		var matrix = M4.forSys(x, y, z, this.l3Axis.anchor)
		var matrixInverse = matrix.inversed()
		var f = this.FofR
		return function (pWC) {
			var p = matrixInverse.transformPoint(pWC)
			var radiusLC = Math.sqrt(p.x * p.x + p.y * p.y)
			return f(p.z) - radiusLC
		}
	}

	boundsFunction() {
		let ptpf = this.pointToParameterFunction()
		let tMin = this.tMin, tMax = this.tMax
		return function (pWC) {
			let pointParameterT = ptpf(pWC).x
			return tMin <= pointParameterT && pointParameterT <= tMax
		}
	}

	pointToParameterFunction() {
		let projPlane = P3(this.dir1, 0)
		let dir1 = this.dir1, baseCurve = this.baseCurve
		let projBaseCurve = baseCurve.project(projPlane)
		let _this = this
		return function (pWC) {
			let projPoint = projPlane.projectedPoint(pWC)
			let t = projBaseCurve.pointLambda(projPoint)
			let z = pWC.minus(baseCurve.at(t)).dot(dir1)
			return V3.create(t, z, 0)
		}

		// let baseCurve = this.baseCurve, dir1 = this.dir1
		// let tMin = this.tMin, tMax = this.tMax
		// return function (pWC) {
		// 	let iterationFunc = t => {let d = baseCurve.at(t).minus(pWC); return d.dot(dir1) + d.length() }
		// 	let startT = NLA.arrayFromFunction(16, i => tMin + (tMax - tMin) * i / 15)
		// 		.withMax(t => -abs(iterationFunc(t)))
		// 	let t = newtonIterate(iterationFunc, startT, 16)
		// 	let z = pWC.minus(baseCurve.at(t)).dot(dir1)
		// 	return V3.create(t, z, 0)
		// }
	}

	/**
	 * @inheritDoc
	 */
	isCurvesWithSurface(surface) {
		// prefer other surface to be the paramteric one
		if (surface.parametricFunction) {
			return new PICurve(surface, this)
		} else if (surface.implicitFunction) {
			return new PICurve(this, surface)
		}
	}

	/**
	 * @inheritDoc
	 */
	containsPoint(p) {
		let pp = this.pointToParameterFunction()(p)
		return this.parametricFunction()(pp.x, pp.y).like(p)
	}

	/**
	 * @inheritDoc
	 */
	containsCurve(curve) {
		if (curve instanceof BezierCurve) {
			// project baseCurve and test curve onto a common plane and check if the curves are alike
			let projectionPlane = P3(this.dir1, 0)
			let baseCurveProjection = this.baseCurve.project(projectionPlane)
			if (curve instanceof L3 && curve.dir1.isParallelTo(this.dir1)) {
				// projection onto basePlane would be single point
				return baseCurveProjection.containsPoint(projectionPlane.projectedPoint(curve.anchor))
			}
			let curveProjection = curve.project(projectionPlane)
			return baseCurveProjection.likeCurve(curveProjection)
		}
		if (curve instanceof L3) {
			return this.dir1.isParallelTo(curve.dir1) && this.containsPoint(curve.anchor)
		}
		return false
	}

	isCoplanarTo(surface) {
		return ProjectedCurveSurface == surface.constructor
				&& this.dir1.isParallelTo(surface.dir1)
				&& this.containsCurve(surface.curve)
	}

	edgeLoopContainsPoint(contour, p) {
		// TODO: this is copied from CylinderSurface, mb create common super class?
		assertVectors(p)
		assert(isFinite(p.x), p.y, p.z)
		var testLine = L3(p, this.dir1)
		let ptpf = this.pointToParameterFunction()
		let pp = ptpf(p)
		if (isNaN(pp.x)) {
			console.log(this.sce, p.sce)
			assert(false)
		}
		var intersectionLinePerpendicular = this.baseCurve.tangentAt(pp.x).rejectedFrom(this.dir1)
		var plane2 = P3.normalOnAnchor(intersectionLinePerpendicular, p)
		var colinearSegments = contour.map((edge) => edge.colinearToLine(testLine))
		var colinearSegmentsInside = contour.map((edge, i) => edge.aDir.dot(this.dir1) > 0)
		var inside = false

		function logIS(p) {
			if (testLine.pointLambda(p) > 0) {
				inside = !inside
			}
		}

		contour.forEach(function (/** @type {B2.Edge} */ edge, /** @type number */ i, /** @type Edge[] */ edges) {
			var j = (i + 1) % edges.length, nextEdge = edges[j]
			//console.log(edge.toSource()) {p:V3(2, -2.102, 0),
			if (colinearSegments[i]) {
				// edge colinear to intersection
				var outVector = edge.bDir.cross(this.normalAt(edge.b))
				var insideNext = outVector.dot(nextEdge.aDir) > 0
				if (colinearSegmentsInside[i] != insideNext) {
					logIS(edge.b)
				}
			} else {
				var edgeTs = edge.edgeISTsWithPlane(plane2)
				for (var k = 0; k < edgeTs.length; k++) {
					var edgeT = edgeTs[k]
					let isp = edge.curve.at(edgeT)
					if (!NLA.equals(pp.x, ptpf(isp).x)) {
						// point is on plane, but not on line
						continue
					}
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


	/**
	 *
	 * @param {M4} m4
	 * @returns {ProjectedCurveSurface}
	 */
	transform(m4) {
		return new ProjectedCurveSurface(this.baseCurve.transform(m4), m4.transformVector(this.dir1).normalized(), this.tMin, this.tMax)
	}

	/**
	 * @inheritDoc
	 */
	isPointsWithLine(line) {
		assertInst(L3, line)
		let projectionPlane = P3(this.dir1, 0)
		let projDir = projectionPlane.projectedVector(line.dir1)
		if (projDir.isZero()) {
			// line is parallel to this.dir
			return []
		}
		let projAnchor = projectionPlane.projectedPoint(line.anchor)
		let projBaseCurve = this.baseCurve.project(projectionPlane)
		return projBaseCurve
			.isInfosWithLine(projAnchor, projDir)
			.map(info => line.at(info.tOther))
	}


	flipped() {
		return new ProjectedCurveSurface(this.baseCurve, this.dir1.negated(), this.tMin, this.tMax)
	}
}
NLA.registerClass(ProjectedCurveSurface)