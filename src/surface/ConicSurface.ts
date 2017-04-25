class ConicSurface extends Surface {
	readonly matrix: M4
	readonly inverseMatrix: M4
	readonly normalMatrix: M4
	readonly normalDir: number // -1 | 1

	/**
	 * returns new cone C = {apex + f1 * z * cos(d) + f2 * z * sin(d) + f3 * z | -PI <= d <= PI, 0 <= z}
	 * @param f1
	 * @param f2
	 * @param f3 Direction in which the cone opens. The ellipse spanned by f1, f2 is contained at (apex + f1).
	 * @param normalDir 1 or -1. 1 if surface normals point outwards, -1 if not.
	 */
	constructor(readonly center: V3,
	            readonly f1: V3,
	            readonly f2: V3,
	            readonly dir: V3) {
		super()
		assertVectors(center, f1, f2, dir)
		this.matrix = M4.forSys(f1, f2, dir, center)
		this.inverseMatrix = this.matrix.inversed()
		this.normalDir = sign(this.f1.cross(this.f2).dot(this.dir))
		this.normalMatrix = this.matrix.as3x3().inversed().transposed().timesScalar(this.normalDir)
	}

	get apex() {
		return this.center
	}

	like(object) {
		assert(false)
		if (!this.isCoplanarTo(object)) return false
		// normals need to point in the same direction (outwards or inwards) for both
		return this.normalDir == object.normalDir
	}
	getVectors() {
		return [    {anchor: this.center, dir1: this.dir},
					{anchor: this.center.plus(this.dir), dir1: this.f1},
					{anchor: this.center.plus(this.dir), dir1: this.f2}]
	}

	getSeamPlane(): P3 {
		return P3.forAnchorAndPlaneVectors(this.center, this.f1, this.dir)
	}

	loopContainsPoint(contour: Edge[], p: V3): PointVsFace {
		assertVectors(p)
		const line = this.center.like(p)
			? new L3(p, this.matrix.transformVector(new V3(0, 1, 1).unit()))
			: L3.throughPoints(p, this.apex)
		const lineOut = line.dir1.cross(this.dir)

		return Surface.loopContainsPointGeneral(contour, p, line, lineOut)
	}

	toSource(): string {
		return makeGen('new ConicSurface', this.center, this.f1, this.f2, this.dir)
	}

	isTsForLine(line: L3): number[] {
		// transforming line manually has advantage that dir1 will not be renormalized,
		// meaning that calculated values t for lineLC are directly transferable to line
		const anchorLC = this.inverseMatrix.transformPoint(line.anchor)
		const dirLC = this.inverseMatrix.transformVector(line.dir1)
		return ConicSurface.unitISLineTs(anchorLC, dirLC)
	}

	/**
	 * Interestingly, two cones don't need to have parallel dirs to be coplanar.
	 */
	isCoplanarTo(surface: Surface): boolean {
		if (this === surface) return true
		if (!(surface instanceof ConicSurface) || !this.apex.like(surface.apex)) return false
		// at this point apexes are equal
		return this.containsEllipse(
			new SemiEllipseCurve(surface.center.plus(surface.dir), surface.f1, surface.f2))
	}

	containsEllipse(ellipse: SemiEllipseCurve): boolean {
		const ellipseLC = ellipse.transform(this.inverseMatrix)
		if (ellipseLC.center.z < 0) {
			return false
		}
		const {f1, f2} = ellipseLC.mainAxes()
		const p1 = ellipseLC.center.plus(f1), p2 = ellipseLC.center.plus(f2)
		// check if both endpoints are on the cone's surface
		// and that one main axis is perpendicular to the Z-axis
		return NLA.eq(p1.x ** 2 + p1.y ** 2, p1.z ** 2)
			&& NLA.eq(p2.x ** 2 + p2.y ** 2, p2.z ** 2)
			&& (NLA.eq0(f1.z) || NLA.eq0(f2.z))
	}

	containsLine(line: L3): boolean {
		const lineLC = line.transform(this.inverseMatrix)
		const d = lineLC.dir1
		return lineLC.containsPoint(V3.O) && NLA.eq(d.x * d.x + d.y * d.y, d.z * d.z)
	}

	containsParabola(parabola: ParabolaCurve): boolean {
		assertInst(ParabolaCurve, parabola)
		const parabolaLC = parabola.transform(this.inverseMatrix)
		if (parabolaLC.center.z < 0 || parabolaLC.f2.z < 0) {
			return false
		}
		const {center, f1, f2} = parabolaLC.rightAngled()
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

	transform(m4: M4): ConicSurface {
		return new ConicSurface(
			m4.transformPoint(this.center),
			m4.transformVector(this.f1).times(m4.isMirroring() ? -1 : 1),
			m4.transformVector(this.f2),
			m4.transformVector(this.dir))
	}

	rightAngled() {

	}

	flipped(): ConicSurface {
		return new ConicSurface(this.center, this.f1.negated(), this.f2, this.dir)
	}

	toMesh(zStart = 0, zEnd = 30) {
		return GL.Mesh.parametric(this.parametricFunction(), this.parametricNormal(), 0, PI, this.tMin, this.tMax, 16, 1)
	}

	parametricNormal() {
		const {f1, f2} = this, f3 = this.dir
		return (d, z) => {
			return f2.cross(f1).plus(f2.cross(f3.times(Math.cos(d)))).plus(f3.cross(f1.times(Math.sin(d)))).unit()
		}
	}

	normalAt(p) {
		const pLC = this.inverseMatrix.transformPoint(p)
		return this.parametricNormal()(pLC.angleXY(), pLC.z)
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
			return this.normalDir * (radiusLC - pLC.z)
		}
	}

	containsPoint(p: V3) {
		return NLA.eq0(this.implicitFunction()(p))
	}

	boundsFunction() {
		assert(false)
	}

	pointToParameterFunction() {
		return (pWC: V3, hint = PI) => {
			const pLC = this.inverseMatrix.transformPoint(pWC)
			let angle = pLC.angleXY()
			if (abs(angle) > Math.PI - NLA_PRECISION) {
				assert(hint == -PI || hint == PI)
				angle = hint
			}
			return new V3(angle, pLC.z, 0)
		}
	}

	isCurvesWithSurface(surface: Surface): Curve[] {
		if (surface instanceof PlaneSurface) {
			return this.isCurvesWithPlane(surface.plane)
		} else if (surface instanceof ConicSurface) {
			if (surface.dir.isParallelTo(this.dir)) {
				const ellipseProjected = surface.transform(M4.projection(this.getPlane(), this.dir))
				return this.isPointsWithEllipse(ellipseProjected).map(is => new L3(is, this.dir))
			} else if (NLA.eq0(this.getCenterLine().distanceToLine(surface.getCenterLine()))) {

			} else {
				assert(false)
			}
		}
	}

	getCenterLine() {
		return new L3(this.center, this.dir)
	}

	isCurvesWithPlane(plane: P3): Curve[] {
		assertInst(P3, plane)
		const planeLC = plane.transform(this.inverseMatrix)
		const planeNormal = planeLC.normal1
		const c = planeNormal.z
		/** "rotate" plane normal1 when passing to {@link ConicSurface.unitISPlane} so that
		 *  y-component of normal1 is 0 */
		const a = planeNormal.lengthXY()
		const d = planeLC.w
		// generated curves need to be rotated back before transforming to world coordinates
		const wcMatrix = eq0(planeNormal.lengthXY())
			? this.matrix
			: this.matrix.times(M4.rotationZ(planeNormal.angleXY()))
		return ConicSurface.unitISPlane(a, c, d).map(curve => curve.transform(wcMatrix)) as Curve[]
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
					return at.minus(this.center).cross(tangent.rejectedFrom(this.dir)).length() / 2
				}
				// ellipse with normal1 parallel to dir1 need to be counted negatively so CCW faces result in a positive area
				// hyperbola normal1 can be perpendicular to
				const sign = edge.curve instanceof SemiEllipseCurve
					? -Math.sign(edge.curve.normal.dot(this.dir))
					: -Math.sign(this.center.to(edge.curve.center).cross(edge.curve.f1).dot(this.dir))
				return glqInSteps(f, edge.aT, edge.bT, 4) * sign
			} else if (edge.curve instanceof L3) {
				return 0
			} else {
				assertNever()
			}
		}).sum()
		// if the cylinder faces inwards, CCW faces will have been CW, so we need to reverse that here
		// Math.abs is not an option as "holes" may also be passed
		return totalArea * Math.sign(this.normal.dot(this.dir))
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
	 * scaling = tangentAt(t) DOT dir.cross(V3.Z).unit()
	 */
	zDirVolume(edges: Edge[]): {volume: number, centroid: V3} {
		// INT[edge.at; edge.bT] (at(t) DOT dir) * (at(t) - at(t).projectedOn(dir) / 2).z
		const totalVolume = edges.map(edge => {
			if (edge.curve instanceof SemiEllipseCurve || edge.curve instanceof HyperbolaCurve || edge.curve instanceof ParabolaCurve) {
				const f = (t) => {
					const at = edge.curve.at(t), tangent = edge.tangentAt(t)
					return (at.z + at.rejectedFrom(this.dir).z) / 2 * at.projectedOn(this.dir).lengthXY() *
						tangent.dot(V3.Z.cross(this.dir).unit())
				}
				// ellipse with normal1 parallel to dir need to be counted negatively so CCW faces result in a positive area
				const sign = edge.curve instanceof SemiEllipseCurve
					? -Math.sign(edge.curve.normal.dot(this.dir))
					: -Math.sign(this.center.to(edge.curve.center).cross(edge.curve.f1).dot(this.dir))
				const val = glqInSteps(f, edge.aT, edge.bT, 1)
				return val * sign
			} else if (edge.curve instanceof L3) {
				return 0
			} else {
				assertNever()
			}
		}).sum()

		return {volume: totalVolume * Math.sign(this.normal.dot(this.dir))}
	}




	static atApexThroughEllipse(apex: V3, ellipse: SemiEllipseCurve): ConicSurface {
		assertVectors(apex)
		assertInst(SemiEllipseCurve, ellipse)
		return new ConicSurface(apex, ellipse.f1, ellipse.f2, apex.to(ellipse.center))
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
				// plane goes through origin/V3.O
				return [new L3(V3.O, new V3(0, Math.sqrt(2) / 2, Math.sqrt(2) / 2)),
					new L3(V3.O, new V3(0, Math.sqrt(2) / 2, -Math.sqrt(2) / 2))]
			} else {
				// hyperbola
				const center = new V3(d / a, 0, 0)
				const f1 = new V3(0, 0, Math.abs(d / a)) // abs, because we always want the hyperbola to be pointing up
				const f2 = new V3(0, d / a, 0)
				return [new HyperbolaCurve(center, f1, f2)]
			}

		} else {
			// c != 0
			const aa = a * a, cc = c * c
			if (NLA.eq0(d)) {
				if (NLA.eq(aa, cc)) {
					return [new L3(V3.O, new V3(c, 0, -a).unit())]
				} else if (aa < cc) {
					assert(false, 'intersection is single point V3.O')
				} else if (aa > cc) {
					return [new L3(V3.O, new V3(c, Math.sqrt(aa - cc), -a).unit()),
						new L3(V3.O, new V3(c, -Math.sqrt(aa - cc), -a).unit())]
				}
			} else {
				if (NLA.eq(aa, cc)) {
					console.log('acd', a, c, d)
					// parabola
					const parabolaVertex = new V3(d / 2 / a, 0, d / 2 / c)
					const parabolaVertexTangentPoint = new V3(d / 2 / a, d / c, d / 2 / c)
					const p2 = new V3(0, 0, d / c)
					return [new ParabolaCurve(parabolaVertex, parabolaVertexTangentPoint.minus(parabolaVertex), p2.minus(parabolaVertex))]
				} else if (aa < cc) {
					// ellipse
					const center = new V3(-a * d / (cc - aa), 0, d * c / (cc - aa))
					if (center.z < 0) {
						return []
					}
					const p1 = new V3(d / (a - c), 0, -d / (a - c))
					const p2 = new V3(-a * d / (cc - aa), d / Math.sqrt(cc - aa), d * c / (cc - aa))
					return [new SemiEllipseCurve(center, center.to(p1), center.to(p2))]
				} else if (aa > cc) {
					// hyperbola
					const center = new V3(-a * d / (cc - aa), 0, -d * c / (cc - aa))
					let f1 = new V3(d / (a - c) - center.x, 0, d / (a - c) - center.z)
					if (f1.z < 0) {
						f1 = f1.negated()
					}
					// sqrt(cc - aa) flipped relative to ellipse case:
					const p2 = new V3(-a * d / (cc - aa), d / Math.sqrt(aa - cc), -d * c / (cc - aa))
					return [new HyperbolaCurve(center, f1, p2.minus(center))]
				}
			}
		}

	}

	/**
	 * Unit cone. x² + y² = z², 0 <= z
	 */
	static readonly UNIT = new ConicSurface(V3.O, V3.X, V3.Y, V3.Z)
}
ConicSurface.prototype.uStep = PI / 16
ConicSurface.prototype.vStep = 256
ConicSurface.prototype.sMin = 0
ConicSurface.prototype.sMax = PI
ConicSurface.prototype.tMin = 0
ConicSurface.prototype.tMax = 16
