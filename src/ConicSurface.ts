class ConicSurface extends Surface {
	baseEllipse: SemiEllipseCurve
	dir: V3
	matrix: M4
	inverseMatrix: M4

	constructor(baseEllipse, dir) {
		super()
		assertVectors(dir)
		assertInst(SemiEllipseCurve, baseEllipse)
		assert(!baseEllipse.normal.isPerpendicularTo(dir), !baseEllipse.normal.isPerpendicularTo(dir))
		this.baseEllipse = baseEllipse
		this.dir = dir
		this.matrix = M4.forSys(baseEllipse.f1, baseEllipse.f2, dir, baseEllipse.center)
		this.inverseMatrix = this.matrix.inversed()
	}

	get apex() {
		return this.baseEllipse.center
	}

	getSeamPlane(): P3 {
		return P3.forAnchorAndPlaneVectors(this.baseEllipse.center, this.baseEllipse.f1, this.dir)
	}

	loopContainsPoint(contour: Edge[], p: V3): PointVsFace {
		assertVectors(p)
		const line = L3.throughPoints(this.apex, p)
		// create plane that goes through cylinder seam
		const seamBase = this.baseEllipse.at(Math.PI)
		const intersectionLinePerpendicular = this.dir.cross(p.minus(seamBase))
		const plane2 = P3.normalOnAnchor(intersectionLinePerpendicular, p)
		const colinearSegments = contour.map((edge) => edge.colinearToLine(line))
		const colinearSegmentsInside = contour.map((edge, i) => edge.aDir.dot(this.dir) > 0)
		let inside = false

		function logIS(p) {
			if (line.pointT(p) > 0) {
				inside = !inside
			}
		}

		contour.forEach((edge, i, edges) => {
			const j = (i + 1) % edges.length, nextEdge = edges[j]
			//console.log(edge.toSource()) {p:V(2, -2.102, 0),
			if (colinearSegments[i]) {
				// edge colinear to intersection
				const outVector = edge.bDir.cross(this.normalAt(edge.b))
				const insideNext = outVector.dot(nextEdge.aDir) > 0
				if (colinearSegmentsInside[i] != insideNext) {
					logIS(edge.b)
				}
			} else {
				const edgeTs = edge.isTsWithPlane(plane2)
				for (let k = 0; k < edgeTs.length; k++) {
					const edgeT = edgeTs[k]
					if (edgeT == edge.bT) {
						// endpoint lies on intersection line
						if (colinearSegments[j]) {
							// next segment is colinear
							// we need to calculate if the section of the plane intersection line BEFORE the colinear
							// segment is inside or outside the face. It is inside when the colinear segment out vector
							// and the current segment vector point in the same direction (dot > 0)
							const colinearSegmentOutsideVector = nextEdge.aDir.cross(plane.normal)
							const insideFaceBeforeColinear = colinearSegmentOutsideVector.dot(edge.bDir) < 0
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

	toString() {
		return makeGen('new ConicSurface', this.baseEllipse, this.dir)
	}

	isTsForLine(line: L3): number[] {
		// transforming line manually has advantage that dir1 will not be renormalized,
		// meaning that calculated values t for localLine are directly transferable to line
		const localAnchor = this.inverseMatrix.transformPoint(line.anchor)
		const localDir = this.inverseMatrix.transformVector(line.dir1)
		return ConicSurface.unitISLineTs(localAnchor, localDir)
	}

	/**
	 * Interestingly, two cones don't need to have parallel dirs to be coplanar.
	 */
	isCoplanarTo(surface: Surface): boolean {
		if (this == surface) return true
		if (!(surface instanceof ConicSurface) || !this.apex.like(surface.apex)) return false
		// at this point apexes are equal
		const ell = surface.baseEllipse
		return this.containsEllipse(
			new SemiEllipseCurve(ell.center.plus(surface.dir), ell.f1, ell.f2))
	}

	containsEllipse(ellipse: SemiEllipseCurve): boolean {
		console.log(ellipse.toString())
		const localEllipse = ellipse.transform(this.inverseMatrix)
		if (localEllipse.center.z < 0) {
			return false
		}
		const {f1, f2} = localEllipse.mainAxes()
		console.log(f1.sce, f2.sce)
		const p1 = localEllipse.center.plus(f1), p2 = localEllipse.center.plus(f2)
		// check if both endpoints are on the cone's surface
		// and that one main axis is perpendicular to the Z-axis
		return NLA.eq(p1.x ** 2 + p1.y ** 2, p1.z ** 2)
			&& NLA.eq(p2.x ** 2 + p2.y ** 2, p2.z ** 2)
			&& (NLA.eq0(f1.z) || NLA.eq0(f2.z))
	}

	containsLine(line: L3): boolean {
		const localLine = line.transform(this.inverseMatrix)
		const d = localLine.dir1
		return localLine.containsPoint(V3.ZERO) && NLA.eq(d.x * d.x + d.y * d.y, d.z * d.z)
	}

	containsParabola(parabola: ParabolaCurve): boolean {
		assertInst(ParabolaCurve, parabola)
		const local = parabola.transform(this.inverseMatrix)
		if (local.center.z < 0 || local.f2.z < 0) {
			return false
		}
		const {center, f1, f2} = local.rightAngled()
		console.log(f1.sce, f2.sce)
		// check if center is on the surface,
		// that tangent is perpendicular to the Z-axis
		// and that "y" axis is parallel to surface
		return NLA.eq(center.x * center.x + center.y * center.y, center.z * center.z)
			&& NLA.eq0(f1.z)
			&& NLA.eq(f2.x * f2.x + f2.y * f2.y, f2.z * f2.z)

	}

	containsCurve(curve) {
		if (curve instanceof SemiEllipseCurve) {
			return this.containsEllipse(curve)
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
			this.normalDir) as this
	}

	flipped() {
		return new ConicSurface(
			this.baseEllipse,
			this.dir.negated())
	}

	toMesh(zStart = 0, zEnd = 30) {
		const mesh = new GL.Mesh({triangles: true, lines: false, normals: true})
		const pF = this.parametricFunction(), pN = this.parametricNormal()
		const split = 4 * 128, inc = 2 * Math.PI / split
		const c = split * 2
		for (let i = 0; i < split; i++) {
			let v = pF(i * inc, zStart)
			mesh.vertices.push(pF(i * inc, zStart), pF(i * inc, zEnd))
			pushQuad(mesh.triangles, -1 != this.normalDir, 2 * i, (2 * i + 2) % c, (2 * i + 1), (2 * i + 3) % c)
			const normal = pN(i * inc, 0)
			mesh.normals.push(normal, normal)
		}
		//mesh.compile()
		return mesh
	}

	parametricNormal() {
		const {f1, f2} = this.baseEllipse, f3 = this.dir
		return (d, z) => {
			return f2.cross(f1).plus(f2.cross(f3.times(Math.cos(d)))).plus(f3.cross(f1.times(Math.sin(d))))
		}
	}

	normalAt(p) {
		const localP = this.inverseMatrix.transformPoint(p)
		return this.parametricNormal()(localP.angleXY(), localP.z)
	}

	parametricFunction() {
		return (d, z) => {
			// center + f1 z cos d + f2 z sin d + z dir
			return this.matrix.transformPoint(new V3(z * cos(d), z * sin(d), z))
		}
	}

	implicitFunction() {
		return (pWC) => {
			const pLC = this.inverseMatrix.transformPoint(pWC)
			const radiusLC = pLC.lengthXY()
			return this.normalDir * (1 - radiusLC)
		}
	}

	containsPoint(p) {
		return NLA.eq0(this.implicitFunction()(p))
	}

	boundsFunction() {
		assert(false)
	}

	pointToParameterFunction() {
		return (pWC, hint) => {
			const pLC = this.inverseMatrix.transformPoint(pWC)
			let angle = pLC.angleXY()
			if (abs(angle) > Math.PI - NLA_PRECISION) {
				assert(hint == -PI || hint == PI)
				angle = hint
			}
			return new V3(angle, pLC.z, 0)
		}
	}

	isCurvesWithSurface(surface2) {
		if (surface2 instanceof PlaneSurface) {
			return this.isCurvesWithPlane(surface2.plane)
		} else if (surface2 instanceof ConicSurface) {
			if (surface2.dir.isParallelTo(this.dir)) {
				const ellipseProjected = surface2.baseEllipse.transform(M4.projection(this.baseEllipse.getPlane(), this.dir))
				return this.baseEllipse.isPointsWithEllipse(ellipseProjected).map(is => new L3(is, this.dir))
			} else if (NLA.eq0(this.getCenterLine().distanceToLine(surface2.getCenterLine()))) {

			} else {
				assert(false)
			}
		}
	}

	getCenterLine() {
		return new L3(this.baseEllipse.center, this.dir)
	}

	isCurvesWithPlane(plane) {
		assertInst(P3, plane)
		const localPlane = plane.transform(this.inverseMatrix)
		const planeNormal = localPlane.normal
		const c = planeNormal.z
		/** "rotate" plane normal when passing to {@link ConicSurface.unitISPlane} so that
		 *  y-component of normal is 0 */
		const a = planeNormal.lengthXY()
		const d = localPlane.w
		// generated curves need to be rotated back before transforming to world coordinates
		const wcMatrix = this.matrix.times(M4.rotationZ(planeNormal.angleXY()))
		return ConicSurface.unitISPlane(a, c, d).map(curve => curve.transform(wcMatrix))
	}

	edgeLoopCCW(contour) {
		if (contour.length < 56) {
			let totalAngle = 0
			for (let i = 0; i < contour.length; i++) {
				const ipp = (i + 1) % contour.length
				const edge = contour[i], nextEdge = contour[ipp]
				totalAngle += edge.bDir.angleRelativeNormal(nextEdge.aDir, this.normalAt(edge.b))
			}
			return totalAngle > 0
		} else {
			const ptpF = this.pointToParameterFunction()
			return isCCW(contour.map(e => ptpF(e.a)), new V3(0, 0, this.normalDir))
		}
	}



	calculateArea(edges: Edge[]): number {
		// calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
		const totalArea = edges.map(edge => {
			if (edge.curve instanceof SemiEllipseCurve || edge.curve instanceof HyperbolaCurve || edge.curve instanceof ParabolaCurve) {
				const f = (t) => {
					const at = edge.curve.at(t), tangent = edge.tangentAt(t)
					console.log(t, at.distanceTo(this.baseEllipse.center), at.rejectedLength(this.dir), tangent.rejectedLength(this.dir))
					return at.minus(this.baseEllipse.center).cross(tangent.rejectedFrom(this.dir)).length() / 2
				}
				// ellipse with normal parallel to dir1 need to be counted negatively so CCW faces result in a positive area
				// hyperbola normal can be perpendicular to
				const sign = edge.curve instanceof SemiEllipseCurve
					? -Math.sign(edge.curve.normal.dot(this.dir))
					: -Math.sign(this.baseEllipse.center.to(edge.curve.center).cross(edge.curve.f1).dot(this.dir))
				return glqInSteps(f, edge.aT, edge.bT, 4) * sign
			} else if (edge.curve instanceof L3) {
				return 0
			} else {
				assertNever()
			}
		}).sum()
		// if the cylinder faces inwards, CCW faces will have been CW, so we need to reverse that here
		// Math.abs is not an option as "holes" may also be passed
		return totalArea * Math.sign(this.baseEllipse.normal.dot(this.dir))
	}

	/**
	 * at(t)
	 * |\                                    ^
	 * | \ at(t).projectedOn(dir)             \  dir
	 * |  \                                    \
	 * |   \ at(t).rejectedFrom(dir)
	 * |   |
	 * |___|
	 *        z = 0
	 *
	 *
	 * A = ((at(t) + at(t).rejectedFrom(dir)) / 2).z * at(t).projectedOn(dir).lengthXY()
	 * scaling = tangentAt(t) DOT dir.cross(V3.Z).normalized()
	 */
	zDirVolume(edges: Edge[]): {volume: number} {
		// INT[edge.at; edge.bT] (at(t) DOT dir) * (at(t) - at(t).projectedOn(dir) / 2).z
		const totalVolume = edges.map(edge => {
			if (edge.curve instanceof SemiEllipseCurve || edge.curve instanceof HyperbolaCurve || edge.curve instanceof ParabolaCurve) {
				const f = (t) => {
					const at = edge.curve.at(t), tangent = edge.tangentAt(t)
					console.log("subarea", t, (at.z + at.rejectedFrom(this.dir).z) / 2 * at.projectedOn(this.dir).lengthXY(), tangent.dot(V3.Z.cross(this.dir).normalized()))
					return (at.z + at.rejectedFrom(this.dir).z) / 2 * at.projectedOn(this.dir).lengthXY() *
						tangent.dot(V3.Z.cross(this.dir).normalized())
				}
				// ellipse with normal parallel to dir need to be counted negatively so CCW faces result in a positive area
				const sign = edge.curve instanceof SemiEllipseCurve
					? -Math.sign(edge.curve.normal.dot(this.dir))
					: -Math.sign(this.baseEllipse.center.to(edge.curve.center).cross(edge.curve.f1).dot(this.dir))
				const val = glqInSteps(f, edge.aT, edge.bT, 1)
				console.log("edge", edge, val, sign)
				return val * sign
			} else if (edge.curve instanceof L3) {
				return 0
			} else {
				assertNever()
			}
		}).sum()

		return {volume: totalVolume * Math.sign(this.baseEllipse.normal.dot(this.dir))}
	}



	/**
	 * returns new cone C = {apex + f1 * z * cos(d) + f2 * z * sin(d) + f3 * z | -PI <= d <= PI, 0 <= z}
	 * @param apex Apex of cone.
	 * @param f1
	 * @param f2
	 * @param f3 Direction in which the cone opens. The ellipse spanned by f1, f2 is contained at (apex + f1).
	 * @param normalDir 1 or -1. 1 if surface normals point outwards, -1 if not.
	 */
	static forApexF123(apex, f1, f2, f3, normalDir): ConicSurface {
		return new ConicSurface(new SemiEllipseCurve(apex, f1, f2), f3, normalDir)
	}

	static atApexThroughEllipse(apex: V3, ellipse: SemiEllipseCurve, normalDir: number): ConicSurface {
		assertVectors(apex)
		assertInst(SemiEllipseCurve, ellipse)
		return new ConicSurface(new SemiEllipseCurve(apex, ellipse.f1, ellipse.f2), ellipse.center.minus(apex), normalDir)
	}

	static unitISLineTs(anchor: V3, dir: V3): number[] {
		const {x: ax, y: ay, z: az} = anchor
		const {x: dx, y: dy, z: dz} = dir

		// this cone: x² + y² = z²
		// line: p = anchor + t * dir1
		// split line equation into 3 component equations, insert into cone equation
		// transform to form (a t² + b t + c = 0) and solve with pqFormula
		const a = dx * dx + dy * dy - dz * dz
		const b = 2 * (ax * dx + ay * dy - az * dz)
		const c = ax * ax + ay * ay - az * az
		// cone only defined for 0 <= z, so filter invalid values
		return pqFormula(b / a, c / a).filter(t => 0 < az + t * dz)
	}

	static unitISPlane(a: number, c: number, d: number): Curve[] {
		// plane: ax + cz = d
		// cone: x² + y² = z²
		if (NLA.eq0(c)) {
			// plane is "vertical", i.e. parallel to Y and Z axes
			assert(!NLA.eq0(a))
			// z² - y² = d²/a²
			if (NLA.eq0(d)) {
				// d = 0 => z² - y² = 0 => z² = y² => z = y
				// plane goes through origin/V3.ZERO
				return [new L3(V3.ZERO, new V3(0, Math.sqrt(2) / 2, Math.sqrt(2) / 2)),
					new L3(V3.ZERO, new V3(0, Math.sqrt(2) / 2, -Math.sqrt(2) / 2))]
			} else {
				// hyperbola
				const center = new V3(d / a, 0, 0)
				const f1 = new V3(0, 0, Math.abs(d / a)); // abs, because we always want the hyperbola to be pointing up
				const f2 = new V3(0, d / a, 0)
				return [new HyperbolaCurve(center, f1, f2)]
			}

		} else {
			// c != 0
			const aa = a * a, cc = c * c
			if (NLA.eq0(d)) {
				if (NLA.eq(aa, cc)) {
					return [new L3(V3.ZERO, new V3(c, 0, -a).normalized())]
				} else if (aa < cc) {
					assert(false, 'intersection is single point V3.ZERO')
				} else if (aa > cc) {
					return [new L3(V3.ZERO, new V3(c, Math.sqrt(aa - cc), -a).normalized()),
						new L3(V3.ZERO, new V3(c, -Math.sqrt(aa - cc), -a).normalized())]
				}
			} else {
				if (NLA.eq(aa, cc)) {
					console.log('acd', a, c, d)
					// parabola
					let parabolaVertex = new V3(d / 2 / a, 0, d / 2 / c)
					let parabolaVertexTangentPoint = new V3(d / 2 / a, d / c, d / 2 / c)
					let p2 = new V3(0, 0, d / c)
					return [new ParabolaCurve(parabolaVertex, parabolaVertexTangentPoint.minus(parabolaVertex), p2.minus(parabolaVertex)) as Curve]
				} else if (aa < cc) {
					// ellipse
					let center = new V3(-a * d / (cc - aa), 0, d * c / (cc - aa))
					if (center.z < 0) {
						return []
					}
					let p1 = new V3(d / (a - c), 0, -d / (a - c))
					let p2 = new V3(-a * d / (cc - aa), d / Math.sqrt(cc - aa), d * c / (cc - aa))
					return [new SemiEllipseCurve(center, p1.minus(center), p2.minus(center))]
				} else if (aa > cc) {
					// hyperbola
					let center = new V3(-a * d / (cc - aa), 0, -d * c / (cc - aa))
					let f1 = new V3(d / (a - c) - center.x, 0, d / (a - c) - center.z)
					if (f1.z < 0) {
						f1 = f1.negated()
					}
					// sqrt(cc - aa) flipped relative to ellipse case:
					let p2 = new V3(-a * d / (cc - aa), d / Math.sqrt(aa - cc), -d * c / (cc - aa))
					console.log(center.sce, f1.sce, p2.sce)
					return [new HyperbolaCurve(center, f1, p2.minus(center))]
				}
			}
		}

	}

	/**
	 * Unit cone. x² + y² = z², 0 <= z
	 */
	static readonly UNIT = new ConicSurface(SemiEllipseCurve.UNIT, V3.Z, 1)
}
