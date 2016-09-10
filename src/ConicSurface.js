class ConicSurface extends Surface {
	constructor(baseEllipse, dir, normalDir) {
		super()
		console.log(baseEllipse.toString(), dir.ss, normalDir)
		assertVectors(dir)
		assert(1 == normalDir || -1 == normalDir, "normalDir == 1 || normalDir == -1" + normalDir)
		assertInst(EllipseCurve, baseEllipse)
		assert(!baseEllipse.normal.isPerpendicularTo(dir), !baseEllipse.normal.isPerpendicularTo(dir))
		this.baseEllipse = baseEllipse
		this.dir = dir
		this.normalDir = -1 * normalDir
		this.matrix = M4.forSys(baseEllipse.f1, baseEllipse.f2, dir, baseEllipse.center)
		this.inverseMatrix = this.matrix.inversed()
	}

	get apex() {
		return this.baseEllipse.center
	}

	/**
	 * @inheritDoc
	 */
	edgeLoopContainsPoint(contour, p) {
		assertVectors(p)
		var line = L3.throughPoints(this.apex, p)
		// create plane that goes through cylinder seam
		var seamBase = this.baseEllipse.at(Math.PI)
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
				var edgeTs = edge.isTsWithPlane(plane2)
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

	toString() {
		return "ConicSurface"
	}

	/**
	 *
	 * @param {L3} line
	 * @returns {Array.<V3>}
	 */
	intersectionWithLine(line) {
		// transforming line manually has advantage that dir1 will not be renormalized,
		// meaning that calculated values t for localLine are directly transferable to line
		var localAnchor = this.inverseMatrix.transformPoint(line.anchor)
		var localDir = this.inverseMatrix.transformVector(line.dir1)
		return ConicSurface.unitISLineTs(localAnchor, localDir).map(t => line.at(t))
	}

	/**
	 * Interestingly, two cones don't need to have parallel dirs to be coplanar.
	 * @param surface
	 * @returns {boolean}
	 */
	isCoplanarTo(surface) {
		if (this == surface) return true
		if (!(surface instanceof ConicSurface) || !this.apex.like(surface.apex)) return false
		// at this point apexes are equal
		var ell = surface.baseEllipse
		return this.containsEllipse(
			new EllipseCurve(ell.center.plus(surface.dir), ell.f1, ell.f2))
	}

	/**
	 *
	 * @param {EllipseCurve} ellipse
	 * @returns {boolean}
	 */
	containsEllipse(ellipse) {
		console.log(ellipse.toString())
		var localEllipse = ellipse.transform(this.inverseMatrix)
		if (localEllipse.center.z < 0) {
			return false
		}
		var {f1, f2} = localEllipse.mainAxes()
		console.log(f1.ss, f2.ss)
		var p1 = localEllipse.center.plus(f1), p2 = localEllipse.center.plus(f2)
		// check if both endpoints are on the cone's surface
		// and that one main axis is perpendicular to the Z-axis
		return NLA.equals(p1.x * p1.x + p1.y * p1.y, p1.z * p1.z)
			&& NLA.equals(p2.x * p2.x + p2.y * p2.y, p2.z * p2.z)
			&& (NLA.isZero(f1.z) || NLA.isZero(f2.z))
	}

	/**
	 *
	 * @param {L3} line
	 * @returns {boolean}
	 */
	containsLine(line) {
		var localLine = line.transform(this.inverseMatrix)
		var d = localLine.dir1
		return localLine.contains(V3.ZERO) && NLA.equals(d.x * d.x + d.y * d.y, d.z * d.z)
	}

	/**
	 *
	 * @param {ParabolaCurve} parabola
	 * @returns {boolean}
	 */
	containsParabola(parabola) {
		assertInst(ParabolaCurve, parabola)
		var local = parabola.transform(this.inverseMatrix)
		if (local.center.z < 0 || local.f2.z < 0) {
			return false
		}
		var {center, f1, f2} = local.rightAngled()
		console.log(f1.ss, f2.ss)
		// check if center is on the surface,
		// that tangent is perpendicular to the Z-axis
		// and that "y" axis is parallel to surface
		return NLA.equals(center.x * center.x + center.y * center.y, center.z * center.z)
			&& NLA.isZero(f1.z)
			&& NLA.equals(f2.x * f2.x + f2.y * f2.y, f2.z * f2.z)

	}

	/**
	 * @inheritDoc
	 */
	containsCurve(curve) {
		if (curve instanceof EllipseCurve) {
			return this.containsEllipse(/** @type EllipseCurve */ curve)
		} else if (curve instanceof L3) {
			return this.containsLine(curve)
		} else {
			assert(false)
		}
	}

	transform(m4) {
		return new ConicSurface(
			this.baseEllipse.transform(m4),
			m4.transformVector(this.dir),
			this.normalDir)
	}

	flipped() {
		return new ConicSurface(
			this.baseEllipse,
			this.dir,
			-this.normalDir)
	}

	toMesh(zStart, zEnd) {
		zStart = zStart || 0
		zEnd = zEnd || 30
		var mesh = new GL.Mesh({triangles: true, lines: false, normals: true})
		var pF = this.parametricFunction(), pN = this.parametricNormal()
		var split = 4 * 32, inc = 2 * Math.PI / split
		var c = split * 2
		for (var i = 0; i < split; i++) {
			var v = pF(i * inc, zStart)
			mesh.vertices.push(pF(i * inc, zStart), pF(i * inc, zEnd))
			pushQuad(mesh.triangles, -1 != this.normalDir, 2 * i, (2 * i + 2) % c, (2 * i + 1), (2 * i + 3) % c)
			var normal = pN(i * inc, 0)
			mesh.normals.push(normal, normal)
		}
		mesh.compile()
		return mesh
	}

	parametricNormal() {
		var {f1, f2} = this.baseEllipse, f3 = this.dir
		return (d, z) => {
			return f2.cross(f1).plus(f2.cross(f3.times(cos(d)))).plus(f3.cross(f1.times(sin(d))))
		}
	}

	normalAt(p) {
		var localP = this.inverseMatrix.transformPoint(p)
		return this.parametricNormal()(localP.angleXY(), localP.z)
	}

	parametricFunction() {
		return (d, z) => {
			var {center, f1, f2} = this.baseEllipse
			return center.plus(f1.times(Math.cos(d) * z)).plus(f2.times(Math.sin(d) * z)).plus(this.dir.times(z))
		}
	}

	implicitFunction() {
		return (pWC) => {
			var p = this.inverseMatrix.transformPoint(pWC)
			var radiusLC = p.lengthXY()
			return this.normalDir * (1 - radiusLC)
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
			var p2 = this.inverseMatrix.transformPoint(pWC)
			var angle = p2.angleXY()
			if (angle < -Math.PI + NLA.PRECISION || angle > Math.PI - NLA.PRECISION) {
				angle = hint.dot(this.baseEllipse.f2) < 0
					? Math.PI
					: -Math.PI
			}
			return V3.create(angle, p2.z, 0)
		}
	}

	/**
	 * @inheritDoc
	 */
	isCurvesWithSurface(surface2) {
		if (surface2 instanceof PlaneSurface) {
			return this.isCurvesWithPlane(surface2.plane)
		} else if (surface2 instanceof ConicSurface) {
			if (surface2.dir.isParallelTo(this.dir)) {
				var ellipseProjected = surface2.baseEllipse.transform(M4.projection(this.baseEllipse.getPlane(), this.dir))
				return this.baseEllipse.intersectWithEllipse(ellipseProjected).map(is => L3(is, this.dir))
			} else if (NLA.isZero(this.getCenterLine().distanceToLine(surface2.getCenterLine()))) {

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
		var localPlane = plane.transform(this.inverseMatrix)
		var planeNormal = localPlane.normal
		var c = planeNormal.z
		/** "rotate" plane normal when passing to {@link ConicSurface.unitISPlane} so that
		 *  y-component of normal is 0 */
		var a = planeNormal.lengthXY()
		var d = localPlane.w
		// generated curves need to be rotated back before transforming to world coordinates
		var wcMatrix = this.matrix.times(M4.rotationZ(planeNormal.angleXY()))
		return ConicSurface.unitISPlane(a, c, d).map(curve => curve.transform(wcMatrix))
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
			return isCCW(contour.map(e => ptpF(e.a)), V3.create(0, 0, this.normalDir))
		}
	}

	/**
	 * returns new cone C = {apex + f1 * z * cos(d) + f2 * z * sin(d) + f3 * z | -PI <= d <= PI, 0 <= z}
	 * @param apex Apex of cone.
	 * @param f1
	 * @param f2
	 * @param f3 Direction in which the cone opens. The ellipse spanned by f1, f2 is contained at (apex + f1).
	 * @param normalDir 1 or -1. 1 if surface normals point outwards, -1 if not.
	 * @returns {ConicSurface}
	 */
	static forApexF123(apex, f1, f2, f3, normalDir) {
		return new ConicSurface(new EllipseCurve(apex, f1, f2), f3, normalDir)
	}

	/**
	 *
	 * @param {V3} apex
	 * @param {EllipseCurve} ellipse
	 * @param {number} normalDir
	 * @returns {ConicSurface}
	 */
	static atApexThroughEllipse(apex, ellipse, normalDir) {
		assertVectors(apex)
		assertInst(EllipseCurve, ellipse)
		return new ConicSurface(new EllipseCurve(apex, ellipse.f1, ellipse.f2), ellipse.center.minus(apex), normalDir)
	}

	/**
	 *
	 * @param {V3} anchor
	 * @param {V3} dir
	 * @returns {Array.<number>}
	 */
	static unitISLineTs(anchor, dir) {
		var {x: ax, y: ay, z: az} = anchor
		var {x: dx, y: dy, z: dz} = dir

		// this cone: x² + y² = z²
		// line: p = anchor + t * dir1
		// split line equation into 3 component equations, insert into cone equation
		// transform to form (a t² + b t + c = 0) and solve with pqFormula
		var a = dx * dx + dy * dy - dz * dz
		var b = 2 * (ax * dx + ay * dy - az * dz)
		var c = ax * ax + ay * ay - az * az
		// cone only defined for 0 <= z, so filter invalid values
		return pqFormula(b / a, c / a).filter(t => 0 < az + t * dz)
	}

	/**
	 *
	 * @param {number} a
	 * @param {number} c
	 * @param {number} d
	 * @returns {Array.<L3|HyperbolaCurve|ParabolaCurve|EllipseCurve>}
	 */
	static unitISPlane(a, c, d) {
		// plane: ax + cz = d
		// cone: x² + y² = z²
		if (NLA.isZero(c)) {
			// plane is "vertical", i.e. parallel to Y and Z axes
			assert(!NLA.isZero(a))
			// z² - y² = d²/a²
			if (NLA.isZero(d)) {
				// d = 0 => z² - y² = 0 => z² = y² => z = y
				// plane goes through origin/V3.ZERO
				return [
					L3(V3.ZERO, V3.create(0, Math.sqrt(2) / 2, Math.sqrt(2) / 2)),
					V3(V3.ZERO, V3.create(0, Math.sqrt(2) / 2, -Math.sqrt(2) / 2))]
			} else {
				// hyperbola
				var center = V3.create(d / a, 0, 0)
				var f1 = V3.create(0, 0, Math.abs(d / a)) // abs, because we always want the hyperbola to be pointing up
				var f2 = V3.create(0, d / a, 0)
				return [new HyperbolaCurve(center, f1, f2)]
			}

		} else {
			// c != 0
			var aa = a * a, cc = c * c
			if (NLA.isZero(d)) {
				if (NLA.equals(aa, cc)) {
					return [L3(V3.ZERO, V3.create(c, 0, -a).normalized())]
				} else if (aa < cc) {
					assert(false, 'intersection is single point V3.ZERO')
				} else if (aa > cc) {
					return [
						L3(V3.ZERO, V3.create(c, Math.sqrt(aa - cc), -a).normalized()),
						V3(V3.ZERO, V3.create(c, -Math.sqrt(aa - cc), -a).normalized())]
				}
			} else {
				if (NLA.equals(aa, cc)) {
					console.log('acd', a, c, d)
					// parabola
					let parabolaVertex = V3.create(d / 2 / a, 0, d / 2 / c)
					let parabolaVertexTangentPoint = V3.create(d / 2 / a, d / c, d / 2 / c)
					let p2 = V3.create(0, 0, d / c)
					return [new ParabolaCurve(parabolaVertex, parabolaVertexTangentPoint.minus(parabolaVertex), p2.minus(parabolaVertex))]
				} else if (aa < cc) {
					// ellipse
					let center = V3.create(-a * d / (cc - aa), 0, d * c / (cc - aa))
					if (center.z < 0) {
						return []
					}
					let p1 = V3.create(d / (a - c), 0, -d / (a - c))
					let p2 = V3.create(-a * d / (cc - aa), d / Math.sqrt(cc - aa), d * c / (cc - aa))
					return [new EllipseCurve(center, p1.minus(center), p2.minus(center))]
				} else if (aa > cc) {
					// hyperbola
					let center = V3.create(-a * d / (cc - aa), 0, -d * c / (cc - aa))
					let f1 = V3.create(d / (a - c) - center.x, 0, d / (a - c) - center.z)
					if (f1.z < 0) {
						f1 = f1.negated()
					}
					// sqrt(cc - aa) flipped relative to ellipse case:
					let p2 = V3.create(-a * d / (cc - aa), d / Math.sqrt(aa - cc), -d * c / (cc - aa))
					console.log(center.ss, f1.ss, p2.ss)
					return [new HyperbolaCurve(center, f1, p2.minus(center))]
				}
			}
		}

	}

}
/**
 * Unit cone. x² + y² = z², 0 <= z
 * @type {ConicSurface}
 */
ConicSurface.UNIT = new ConicSurface(EllipseCurve.UNIT, V3.Z, 1)
