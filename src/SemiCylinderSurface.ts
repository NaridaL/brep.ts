class SemiCylinderSurface extends ProjectedCurveSurface {
	readonly matrix: M4
	readonly inverseMatrix: M4
	readonly normalMatrix: M4
	readonly normalDir: number
	readonly baseCurve: SemiEllipseCurve

	constructor(baseCurve: SemiEllipseCurve, dir1: V3, zMin = -Infinity, zMax = Infinity) {
		super(baseCurve, dir1, undefined, undefined, zMin, zMax)
		assertVectors(dir1)
		assertInst(SemiEllipseCurve, baseCurve)
		//assert(!baseCurve.normal.isPerpendicularTo(dir1), !baseCurve.normal.isPerpendicularTo(dir1))
		assert(dir1.hasLength(1))
		this.matrix = M4.forSys(baseCurve.f1, baseCurve.f2, dir1, baseCurve.center)
		this.inverseMatrix = this.matrix.inversed()
		this.normalDir = sign(this.baseCurve.normal.dot(this.dir1))
		this.normalMatrix = this.matrix.as3x3().inversed().transposed().timesScalar(this.normalDir)
	}

	toSource() {
		return `new SemiCylinderSurface(${this.baseCurve.toSource()}, ${this.dir1.toSource()}, ${this.tMin}, ${this.tMax})`
	}

	normalAt(p: V3): V3 {
		return this.normalMatrix.transformVector(this.inverseMatrix.transformPoint(p).xy()).unit()
	}

	loopContainsPoint(loop: Edge[], p: V3): PointVsFace {
		assertVectors(p)
        if (!this.containsPoint(p)) return OUTSIDE
		// create plane that goes through cylinder seam
		const line = new L3(p, this.dir1)
		const seamBase = this.baseCurve.at(PI)
		const lineOut = this.dir1.cross(this.normalAt(p))
		return Surface.loopContainsPointGeneral(loop, p, line, lineOut)
	}

	isTsForLine(line) {
		assertInst(L3, line)
		// transforming line manually has advantage that dir1 will not be renormalized,
		// meaning that calculated values t for localLine are directly transferable to line
		const dirLC = this.inverseMatrix.transformVector(line.dir1)
		if (dirLC.isParallelTo(V3.Z)) {
			// line is parallel to this.dir
			return []
		}
		const anchorLC = this.inverseMatrix.transformPoint(line.anchor)
		assert(!SemiCylinderSurface.unitISLineTs(anchorLC, dirLC).length || !isNaN(SemiCylinderSurface.unitISLineTs(anchorLC, dirLC)[0]), 'sad ' + dirLC)
		return SemiCylinderSurface.unitISLineTs(anchorLC, dirLC)
	}

	isCoplanarTo(surface) {
		return this == surface ||
			surface instanceof SemiCylinderSurface
			&& this.dir1.isParallelTo(surface.dir1)
			&& this.containsSemiEllipse(surface.baseCurve, false)
	}

	like(object) {
		if (!this.isCoplanarTo(object)) return false
		// normals need to point in the same direction (outwards or inwards) for both
		const thisFacesOut = 0 < this.baseCurve.normal.dot(this.dir1)
		const objectFacesOut = 0 < object.baseCurve.normal.dot(object.dir1)
		return thisFacesOut == objectFacesOut
	}

	containsSemiEllipse(ellipse, checkAABB = true) {
		const projEllipse = ellipse.transform(M4.projection(this.baseCurve.getPlane(), this.dir1))
		return this == ellipse || this.baseCurve.isColinearTo(projEllipse) &&
			(!checkAABB || le(0, ellipse.transform(this.inverseMatrix).getAABB().min.y))
	}

	containsCurve(curve) {
		if (curve instanceof L3) {
			return this.containsLine(curve)
		} else if (curve instanceof SemiEllipseCurve) {
			return this.containsSemiEllipse(curve)
		} else if (curve instanceof PICurve) {
			return curve.points.every(p => this.containsPoint(p))
		} else {
			assert(false)
		}
	}

	transform(m4) {
		const newDir = m4.transformVector(this.dir1)
		return new SemiCylinderSurface(
			this.baseCurve.transform(m4),
			newDir.toLength(m4.isMirroring() ? -1 : 1),
			this.tMin * newDir.length(), this.tMax * newDir.length()) as this
	}

	flipped() {
		return new SemiCylinderSurface(
			this.baseCurve,
			this.dir1.negated())
	}

	implicitFunction() {
		return (pWC) => {
			const pLC = this.inverseMatrix.transformPoint(pWC)
			const radiusLC = pLC.lengthXY()
			const normalDir = Math.sign(this.baseCurve.normal.dot(this.dir1))
			return normalDir * (1 - radiusLC)
		}
	}

	containsPoint(p) {
	    const pLC = this.inverseMatrix.transformPoint(p)
		return SemiEllipseCurve.validPlanePoint(pLC.x, pLC.y)
	}

	parameters(pWC) {
        assert(arguments.length == 1)
        const pLC = this.inverseMatrix.transformPoint(pWC)
        const u = SemiEllipseCurve.unitT(pLC)
        return new V3(u, pLC.z, 0)
	}

	isCurvesWithSurface(surface2) {
		if (surface2 instanceof PlaneSurface) {
			return this.isCurvesWithPlane(surface2.plane)
		} else if (surface2 instanceof SemiCylinderSurface) {
			if (surface2.dir1.isParallelTo(this.dir1)) {
				const ellipseProjected = surface2.baseCurve.transform(M4.projection(this.baseCurve.getPlane(), this.dir1))
				return this.baseCurve.isInfosWithEllipse(ellipseProjected).map(info => new L3(info.p, this.dir1))
			} else if (NLA.eq0(this.getCenterLine().distanceToLine(surface2.getCenterLine()))) {
				assert(false)
			} else {
				assert(false)
			}
		}
	}

	getCenterLine() {
		return new L3(this.baseCurve.center, this.dir1)
	}

	static semicylinder(radius: number): SemiCylinderSurface {
		return new SemiCylinderSurface(new SemiEllipseCurve(V3.O, new V3(radius, 0, 0), new V3(0, radius, 0)), V3.Z)
	}

	/**
	 *
	 * @param anchor
	 * @param dir not necessarily unit
	 */
	static unitISLineTs(anchor: V3, dir: V3): number[] {
		const {x: ax, y: ay} = anchor
		const {x: dx, y: dy} = dir

		// this cylinder: x² + y² = 1
		// line: p = anchor + t * dir
		// split line equation into 3 component equations, insert into cylinder equation
		// x = ax + t * dx
		// y = ay + t * dy
		// (ax² + 2 ax t dx + t²dx²) + (ay² + 2 ay t dy + t²dy²) = 1
		// transform to form (a t² + b t + c = 0) and solve with pqFormula
		const a = dx ** 2 + dy ** 2
		const b = 2 * (ax * dx + ay * dy)
		const c = ax ** 2 + ay ** 2 - 1
		return pqFormula(b / a, c / a).filter(t => SemiEllipseCurve.validPlanePoint(ax + dx * t, ay + dy * t))
	}

	/**
	 * Calculating the surface area of a projected ellipse is analogous to the circumference of the ellipse
	 * ==> Elliptic integrals/numeric calculation is necessary
	 */
	calculateArea(edges: Edge[]): number {
		// calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
		const totalArea = edges.map(edge => {
			if (edge.curve instanceof SemiEllipseCurve) {
				const f = (t) => {
					const at = edge.curve.at(t), tangent = edge.tangentAt(t)
					return at.dot(this.dir1) * tangent.rejected1Length(this.dir1)
				}
				// ellipse with normal parallel to dir1 need to be counted negatively so CCW faces result in a positive area
				const sign = -Math.sign(edge.curve.normal.dot(this.dir1))
				const val = glqInSteps(f, edge.aT, edge.bT, 4)
				console.log("edge", edge, val)
				return val * sign
			} else if (edge.curve instanceof L3) {
				return 0
			} else {
				assertNever()
			}
		}).sum()
		// if the cylinder faces inwards, CCW faces will have been CW, so we need to reverse that here
		// Math.abs is not an option as "holes" may also be passed
		return totalArea * Math.sign(this.baseCurve.normal.dot(this.dir1))
	}

	/**
	 * at(t)
	 * |\                                    ^
	 * | \ at(t).projectedOn(dir1)            \  dir1
	 * |  \                                    \
	 * |   \ at(t).rejectedFrom(dir1)
	 * |   |
	 * |___|
	 *        z = 0
	 *
	 *
	 * A = ((at(t) + at(t).rejectedFrom(dir1)) / 2).z * at(t).projectedOn(dir1).lengthXY()
	 * scaling = tangentAt(t) DOT dir1.cross(V3.Z).unit()
	 */
	zDirVolume(edges: Edge[]): {volume: number, centroid: any} {
		if (V3.Z.cross(this.dir1).isZero()) return {volume: 0}
		// the tangent needs to be projected onto a vector which is perpendicular to the volume-slices
		const scalingVector = this.dir1.cross(V3.Z).unit()
		// the length of the base of the trapezoid is calculated by dotting with the baseVector
		const baseVector = this.dir1.rejectedFrom(V3.Z).unit()
		// INT[edge.at; edge.bT] (at(t) DOT dir1) * (at(t) - at(t).projectedOn(dir) / 2).z
		console.log("scalingVector", scalingVector.sce)
		const totalArea = edges.map(edge => {
			if (edge.curve instanceof SemiEllipseCurve) {
				const f = (t) => {
					// use curve.tangent not edge.tangent, reverse edges are handled by the integration boundaries
					const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t)
					const area = (at.z + at.rejectedFrom(this.dir1).z) / 2 * at.projectedOn(this.dir1).dot(baseVector)
					const scale = tangent.dot(scalingVector)
					//assert(Math.sign(scale) == Math.sign(this.normalAt(at).dot(V3.Z)), this.normalAt(at).dot(V3.Z))
					//console.log(
					//	"", t,
					//	",", area,
					//	",", scale,
					//	"atz", at.z)
					return area * scale
				}
				// ellipse with normal parallel to dir1 need to be counted negatively so CCW faces result in a positive area
				const sign = -Math.sign(edge.curve.normal.dot(this.dir1))
				const val = glqInSteps(f, edge.aT, edge.bT, 1)
				console.log("edge", edge, val, sign)
				return val * sign
			} else if (edge.curve instanceof L3) {
				return 0
			} else {
				assertNever()
			}
		}).sum()

		return {volume: totalArea * Math.sign(this.baseCurve.normal.dot(this.dir1))}
	}

	facesOutwards(): boolean {
		return this.baseCurve.normal.dot(this.dir1) > 0
	}

	getSeamPlane(): P3 {
		return P3.forAnchorAndPlaneVectors(this.baseCurve.center, this.baseCurve.f1, this.dir1)
	}

	static readonly UNIT = new SemiCylinderSurface(SemiEllipseCurve.UNIT, V3.Z, 0, 1)
}
SemiCylinderSurface.prototype.uStep = TAU  / 32
SemiCylinderSurface.prototype.vStep = 256
SemiCylinderSurface.prototype.sMin = 0
SemiCylinderSurface.prototype.sMax = PI
SemiCylinderSurface.prototype.tMin = 0
SemiCylinderSurface.prototype.tMax = 1










