class PlaneSurface extends Surface {
	readonly matrix: M4
	constructor(readonly plane: P3,
	            readonly right: V3 = plane.normal1.getPerpendicular().unit(),
	            readonly up: V3 = plane.normal1.cross(right).unit()) {
		super()
		assertInst(P3, plane)
		assert(this.right.cross(this.up).like(this.plane.normal1))
		this.matrix = M4.forSys(right, up, plane.normal1, plane.anchor)
	}

	isCoplanarTo(surface) {
		return surface instanceof PlaneSurface && this.plane.isCoplanarToPlane(surface.plane)
	}


	isTsForLine(line: L3): number[] {
		return line.isTsWithPlane(this.plane)
	}

	like(surface) {
		return surface instanceof PlaneSurface && this.plane.like(surface.plane)
	}

	parametricFunction() {
		return (s, t) => this.matrix.transformPoint(new V3(s, t, 0))
	}

	implicitFunction() {
		return p => this.plane.distanceToPointSigned(p)
	}

	isCurvesWithISurface(implicitSurface) {
		assert(implicitSurface.implicitFunction, 'implicitSurface.implicitFunction')
		return new CurvePI(this, implicitSurface)
	}

	isCurvesWithSurface(surface2: Surface): Curve[] {
		if (surface2 instanceof PlaneSurface) {
			return this.isCurvesWithPlane(surface2.plane)
		}
		return super.isCurvesWithSurface(surface2)
	}

	isCurvesWithPlane(plane: P3): L3[] {
		if (this.plane.isParallelToPlane(plane)) {
			return []
		}
		return [this.plane.intersectionWithPlane(plane)]
	}

	edgeLoopCCW(contour: Edge[]): boolean {
		return isCCW(contour.map(edge => edge.points()).concatenated(), this.plane.normal1)
		let totalAngle = 0
		for (let i = 0; i < contour.length; i++) {
			const ipp = (i + 1) % contour.length
			const edge = contour[i], nextEdge = contour[ipp]
			assert(edge.b.like(nextEdge.a), 'edges dont form a loop')
			if (edge.curve instanceof SemiEllipseCurve) {
				totalAngle += edge.rotViaPlane(this.plane.normal1)
				// console.log(edge.toString(), edge.rotViaPlane(this.plane.normal1))
			}
			totalAngle += edge.bDir.angleRelativeNormal(nextEdge.aDir, this.plane.normal1)
		}
		const result = totalAngle > 0
		const result2 = PlaneFace.prototype.calculateArea.apply({surface: this, contour: contour}).area > 0
		//assert (result == result2)
		return result2
	}

	loopContainsPoint(loop: Edge[], p: V3): PointVsFace {
		const dir = this.right.plus(this.up.times(0.123)).unit()
		const line = new L3(p, dir)
		const lineOut = dir.cross(this.plane.normal1)
		return Surface.loopContainsPointGeneral(loop, p, line, lineOut)
	}

	pointToParameterFunction() {
		const matrix = M4.forSys(this.right, this.up, this.normal, this.plane.anchor)
		const matrixInverse = matrix.inversed()
		return function (pWC) {
			return matrixInverse.transformPoint(pWC)
		}
	}

	normalAt(p) {
		return this.plane.normal1
	}

	containsPoint(p) {
		return this.plane.containsPoint(p)
	}

	containsCurve(curve: Curve): boolean {
		return this.plane.containsCurve(curve)
	}

	transform(m4: M4): this {
		return new PlaneSurface(this.plane.transform(m4)) as this
	}

	flipped() {
		return new PlaneSurface(this.plane.flipped(), this.right, this.up.negated()) as this
	}

	toString() {
		return callsce('new PlaneSurface', this.plane, this.right, this.up)
	}

	toMesh(xMin: number = -10, xMax: number = 10, yMin: number = -10, yMax: number = 10) {
		const mesh = new GL.Mesh({triangles: true, lines: false, normals: true})
		const matrix = M4.forSys(this.right, this.up, this.plane.normal1, this.plane.anchor)
		mesh.vertices = [V(xMin, yMin), V(xMax, yMin), V(xMin, yMax), V(xMax, yMax)].map(p => matrix.transformPoint(p))
		mesh.normals = NLA.arrayFromFunction(4, i => this.plane.normal1)
		pushQuad(mesh.triangles, false, 0, 1, 2, 3)
		mesh.compile()
		return mesh
	}


	static throughPoints(a, b, c): PlaneSurface {
		return new PlaneSurface(P3.throughPoints(a, b, c))
	}
}
NLA.registerClass(PlaneSurface)