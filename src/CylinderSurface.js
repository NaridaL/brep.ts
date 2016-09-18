class CylinderSurface extends Surface {
	constructor(baseEllipse, dir) {
		super()
		assert(2 == arguments.length)
		assertVectors(dir)
		assertInst(EllipseCurve, baseEllipse)
		//assert(!baseEllipse.normal.isPerpendicularTo(dir), !baseEllipse.normal.isPerpendicularTo(dir))
		assert(dir.hasLength(1))
		this.baseEllipse = baseEllipse
		this.dir = dir
		this.matrix = M4.forSys(baseEllipse.f1, baseEllipse.f2, dir, baseEllipse.center)
		this.inverseMatrix = this.matrix.inversed()
	}

	toSource() {
		return `new CylinderSurface(${this.baseEllipse.toSource()}, ${this.dir.toSource()})`
	}

	edgeLoopContainsPoint(contour, p) {
		assertVectors(p)
		var line = L3(p, this.dir)
		// create plane that goes through cylinder seam
		var seamBase = this.baseEllipse.at(PI)
		var intersectionLinePerpendicular = this.dir.cross(p.minus(seamBase))
		var plane2 = P3.normalOnAnchor(intersectionLinePerpendicular, p)
		var colinearSegments = contour.map((edge) => edge.colinearToLine(line))
		var colinearSegmentsInside = contour.map((edge, i) => edge.aDir.dot(this.dir) > 0)
		var inside = false

		function logIS(p) {
			if (line.pointLambda(p) > 0) {
				inside = !inside
			}
		}

		contour.forEach((edge, i, edges) => {
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
					if (edgeT == edge.bT) {
						// endpoint lies on intersection line
						if (colinearSegments[j]) {
							// next segment is colinear
							// we need to calculate if the section of the plane intersection line BEFORE the colinear
							// segment is inside or outside the face. It is inside when the colinear segment out vector
							// and the current segment vector point in the same direction (dot > 0)
							var colinearSegmentOutsideVector = nextEdge.aDir.cross(plane.normal)
							var insideFaceBeforeColinear = colinearSegmentOutsideVector.dot(edge.bDir) < 0
							// if the "inside-ness" changes, add intersection point
							//console.log("segment end on line followed by colinear", insideFaceBeforeColinear !=
							// colinearSegmentInsideFace, nextSegmentOutsideVector)
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
	 * @param {L3} line
	 * @returns {Array.<V3>}
	 */
	isPointsWithLine(line) {
		assertInst(L3, line)
		// transforming line manually has advantage that dir1 will not be renormalized,
		// meaning that calculated values t for localLine are directly transferable to line
		var localDir = this.inverseMatrix.transformVector(line.dir1)
		if (localDir.isParallelTo(V3.Z)) {
			// line is parallel to this.dir
			return []
		}
		var localAnchor = this.inverseMatrix.transformPoint(line.anchor)
		assert(!CylinderSurface.unitISLineTs(localAnchor, localDir).length || !isNaN(CylinderSurface.unitISLineTs(localAnchor, localDir)[0]), 'sad ' +localDir)
		return CylinderSurface.unitISLineTs(localAnchor, localDir).map(t => line.at(t))
	}

	isCoplanarTo(surface) {
		return this == surface ||
			surface instanceof CylinderSurface
			&& this.dir.isParallelTo(surface.dir)
			&& this.containsEllipse(surface.baseEllipse)
	}

	containsEllipse(ellipse) {
		var ellipseProjected = ellipse.transform(M4.projection(this.baseEllipse.getPlane(), this.dir))
		return this == ellipse || this.baseEllipse.isColinearTo(ellipseProjected)
	}

	containsLine(line) {
		return this.dir.isParallelTo(line.dir1) && this.containsPoint(line.anchor)
	}

	containsCurve(curve) {
		if (curve instanceof EllipseCurve) {
			return this.containsEllipse(curve)
		} else if (curve instanceof L3) {
			return this.containsLine(curve)
		} else {
			assert(false)
		}
	}

	/**
	 * @inheritDoc
	 */
	transform(m4) {
		return new CylinderSurface(
			this.baseEllipse.transform(m4),
			m4.transformVector(this.dir))
	}

	flipped() {
		return new CylinderSurface(
			this.baseEllipse,
			this.dir.negated())
	}

	toMesh(zStart, zEnd) {
		zStart = zStart || -30
		zEnd = zEnd || 30
		var mesh = new GL.Mesh({triangles: true, lines: false, normals: true})
		var pF = this.parametricFunction(), pN = this.parametricNormal()
		var split = 4 * 10, inc = 2 * PI / split
		var c = split * 2
		for (var i = 0; i < split; i++) {
			var v = pF(i * inc, zStart)
			mesh.vertices.push(pF(i * inc, zStart), pF(i * inc, zEnd))
			pushQuad(mesh.triangles, false, 2 * i, (2 * i + 2) % c, (2 * i + 1), (2 * i + 3) % c)
			var normal = pN(i * inc, 0)
			mesh.normals.push(normal, normal)
		}
		console.log(mesh)
		//mesh.computeNormalLi00nes()
		mesh.compile()
		return mesh
	}

	parametricNormal() {
		return (d, z) => {
			return this.baseEllipse.tangentAt(d).cross(this.dir).normalized()
		}
	}

	normalAt(p) {
		var localP = this.inverseMatrix.transformPoint(p)
		return this.parametricNormal()(localP.angleXY(), localP.z)
	}

	parametricFunction() {
		return (d, z) => {
			return this.baseEllipse.at(d).plus(this.dir.times(z))
		}
	}

	implicitFunction() {
		return (pWC) => {
			var p = this.inverseMatrix.transformPoint(pWC)
			var radiusLC = p.lengthXY()
			const normalDir = Math.sign(this.baseEllipse.normal.dot(this.dir))
			return normalDir * (1 - radiusLC)
		}
	}

	containsPoint(p) {
		return NLA.isZero(this.implicitFunction()(p))
	}

	boundsFunction() {
		assert(false)
	}

	pointToParameterFunction() {
		return (pWC, hint) => {
			var pLC = this.inverseMatrix.transformPoint(pWC)
			var angle = pLC.angleXY()
			if (angle < -Math.PI + NLA.PRECISION || angle > Math.PI - NLA.PRECISION) {
				angle = Math.sign(hint) * Math.PI
			}
			return V3.create(angle, pLC.z, 0)
		}
	}

	isCurvesWithSurface(surface2) {
		if (surface2 instanceof PlaneSurface) {
			return this.isTsWithPlane(surface2.plane)
		} else if (surface2 instanceof CylinderSurface) {
			if (surface2.dir.isParallelTo(this.dir)) {
				var ellipseProjected = surface2.baseEllipse.transform(M4.projection(this.baseEllipse.getPlane(), this.dir))
				return this.baseEllipse.isInfosWithEllipse(ellipseProjected).map(info => L3(info.p, this.dir))
			} else if (NLA.isZero(this.getCenterLine().distanceToLine(surface2.getCenterLine()))) {
				assert(false)
			} else {
				assert(false)
			}
		}
	}

	getCenterLine() {
		return L3(this.baseEllipse.center, this.dir)
	}

	/**
	 * @inheritDoc
	 */
	isCurvesWithPlane(plane) {
		assertInst(P3, plane)
		if (this.dir.isPerpendicularTo(plane.normal)) {
			var eTs = this.baseEllipse.isTsWithPlane(plane)
			return eTs.map(t => L3(this.baseEllipse.at(t), this.dir))
		} else {
			return [this.baseEllipse.transform(M4.projection(plane, this.dir))]
		}
	}

	edgeLoopCCW(contour) {
		if (contour.length < 56) {
			var totalAngle = 0
			for (var i = 0; i < contour.length; i++) {
				var ipp = (i + 1) % contour.length
				var edge = contour[i], nextEdge = contour[ipp]
				totalAngle += edge.bDir.angleRelativeNormal(nextEdge.aDir, this.normalAt(edge.b))
			}
			return totalAngle > 0
		} else {
			var ptpF = this.pointToParameterFunction()
			return isCCW(contour.map(e => ptpF(e.a)), V3.Z)
		}
	}

	/**
	 *
	 * @param {number} radius
	 * @returns {CylinderSurface}
	 */
	static cyl(radius) {
		return new CylinderSurface(new EllipseCurve(V3.ZERO, V3(radius, 0, 0), V3(0, radius, 0)), V3.Z)
	}

	/**
	 *
	 * @param {V3} anchor
	 * @param {V3} dir not necessarily normalized
	 * @returns {Array.<number>}
	 */
	static unitISLineTs(anchor, dir) {
		var {x: ax, y: ay, z: az} = anchor
		var {x: dx, y: dy, z: dz} = dir

		// this cylinder: x² + y² = 1
		// line: p = anchor + t * dir
		// split line equation into 3 component equations, insert into cylinder equation
		// x = ax + t * dx
		// y = ay + t * dy
		// (ax² + 2 ax t dx + t²dx²) + (ay² + 2 ay t dy + t²dy²) = 1
		// transform to form (a t² + b t + c = 0) and solve with pqFormula
		var a = dx * dx + dy * dy
		var b = 2 * (ax * dx + ay * dy)
		var c = ax * ax + ay * ay - 1
		return pqFormula(b / a, c / a)
	}
}
CylinderSurface.UNIT = new CylinderSurface(EllipseCurve.UNIT, V3.Z)
