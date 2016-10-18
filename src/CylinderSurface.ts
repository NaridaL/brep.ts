class CylinderSurface extends Surface {
	baseEllipse: EllipseCurve
	dir: V3
	matrix:M4
	inverseMatrix:M4
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


	edgeLoopContainsPoint(loop: Edge[], p: V3): PointVsFace {
		assertVectors(p)

		// create plane that goes through cylinder seam
		const line = new L3(p, this.dir)
		const seamBase = this.baseEllipse.at(PI)
		const lineOut = this.dir.cross(p.minus(seamBase))
		return Surface.loopContainsPointGeneral(loop, p, line, lineOut)
	}

	/**
	 * @inheritDoc
	 */
	isTsForLine(line) {
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
		return CylinderSurface.unitISLineTs(localAnchor, localDir)
	}

	isCoplanarTo(surface) {
		return this == surface ||
			surface instanceof CylinderSurface
			&& this.dir.isParallelTo(surface.dir)
			&& this.containsEllipse(surface.baseEllipse)
	}

	/**
	 * @inheritDoc
	 */
	like(object) {
		if (!this.isCoplanarTo(object)) return false
		// normals need to point in the same direction (outwards or inwards) for both
		let thisFacesOut = 0 < this.baseEllipse.normal.dot(this.dir)
		let objectFacesOut = 0 < object.baseEllipse.normal.dot(object.dir)
		return thisFacesOut == objectFacesOut
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
			m4.transformVector(this.dir).normalized())
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
		return NLA.eq0(this.implicitFunction()(p))
	}

	boundsFunction() {
		assert(false)
	}

	pointToParameterFunction() {
		return (pWC, hint) => {
			var pLC = this.inverseMatrix.transformPoint(pWC)
			var angle = pLC.angleXY()
			if (angle < -Math.PI + NLA_PRECISION || angle > Math.PI - NLA_PRECISION) {
				angle = Math.sign(hint) * Math.PI
			}
			return new V3(angle, pLC.z, 0)
		}
	}

	isCurvesWithSurface(surface2) {
		if (surface2 instanceof PlaneSurface) {
			return this.isCurvesWithPlane(surface2.plane)
		} else if (surface2 instanceof CylinderSurface) {
			if (surface2.dir.isParallelTo(this.dir)) {
				var ellipseProjected = surface2.baseEllipse.transform(M4.projection(this.baseEllipse.getPlane(), this.dir))
				return this.baseEllipse.isInfosWithEllipse(ellipseProjected).map(info => new L3(info.p, this.dir))
			} else if (NLA.eq0(this.getCenterLine().distanceToLine(surface2.getCenterLine()))) {
				assert(false)
			} else {
				assert(false)
			}
		}
	}

	getCenterLine() {
		return new L3(this.baseEllipse.center, this.dir)
	}

	isCurvesWithPlane(plane):Curve[] {
		assertInst(P3, plane)
		if (this.dir.isPerpendicularTo(plane.normal)) {
			var ellipseTs = this.baseEllipse.isTsWithPlane(plane)
			return ellipseTs.map(t => {
				let l3dir = 0 < this.baseEllipse.tangentAt(t).dot(plane.normal)
					? this.dir
					: this.dir.negated()
				return new L3(this.baseEllipse.at(t), l3dir)
			})
		} else {
			let projEllipse = this.baseEllipse.transform(M4.projection(plane, this.dir))
			if (this.dir.dot(plane.normal) > 0) {
				// we need to flip the ellipse so the tangent is correct
				projEllipse = new EllipseCurve(projEllipse.center, projEllipse.f1, projEllipse.f2.negated())
		}
			return [projEllipse]
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

	static cylinder(radius:number):CylinderSurface {
		return new CylinderSurface(new EllipseCurve(V3.ZERO, V(radius, 0, 0), V(0, radius, 0)), V3.Z)
	}

	/**
	 *
	 * @param anchor
	 * @param dir not necessarily normalized
	 * @returns {Array.<number>}
	 */
	static unitISLineTs(anchor:V3, loop:V3):number[] {
		var {x: ax, y: ay, z: az} = anchor
		var {x: dx, y: dy, z: dz} = loop

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

	static UNIT = new CylinderSurface(EllipseCurve.UNIT, V3.Z)
}
