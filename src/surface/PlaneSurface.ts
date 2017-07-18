class PlaneSurface extends ParametricSurface implements ImplicitSurface {
	readonly matrix: M4
	constructor(readonly plane: P3,
	            readonly right: V3 = plane.normal1.getPerpendicular().unit(),
	            readonly up: V3 = plane.normal1.cross(right).unit(),
                readonly sMin: number = -100,
                readonly sMax: number = 100,
                readonly tMin: number = -100,
                readonly tMax: number = 100) {
		super()
		assertInst(P3, plane)
		assert(this.right.cross(this.up).like(this.plane.normal1))
		this.matrix = M4.forSys(right, up, plane.normal1, plane.anchor)
	}

	isCoplanarTo(surface: Surface): boolean {
		return surface instanceof PlaneSurface && this.plane.isCoplanarToPlane(surface.plane)
	}

	isTsForLine(line: L3): number[] {
		return line.isTsWithPlane(this.plane)
	}

	like(surface: Surface): boolean {
		return surface instanceof PlaneSurface && this.plane.like(surface.plane)
	}

    pST(s: number, t: number): V3 {
		return this.matrix.transformPoint(new V3(s, t, 0))
	}

    implicitFunction(): (pWC: V3) => number {
		return p => this.plane.distanceToPointSigned(p)
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
		return isCCW(contour.flatMap(edge => edge.points()), this.plane.normal1)
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

	stPFunc() {
		const matrixInverse = this.matrix.inversed()
		return function (pWC: V3) {
			return matrixInverse.transformPoint(pWC)
		}
	}

    pointFoot(pWC: V3): V3 {
        return this.stP(pWC)
    }

    normalP(pWC: V3): V3 {
		return this.plane.normal1
	}

	containsPoint(p) {
		return this.plane.containsPoint(p)
	}

	containsCurve(curve: Curve): boolean {
		return this.plane.containsCurve(curve)
	}

	transform(m4: M4): PlaneSurface {
		return new PlaneSurface(this.plane.transform(m4))
	}

	flipped() {
		return new PlaneSurface(this.plane.flipped(), this.right, this.up.negated())
	}

	getConstructorParameters(): any[] {
		return [this.plane, this.right, this.up]
	}

	toMesh(xMin: number = -10, xMax: number = 10, yMin: number = -10, yMax: number = 10) {
		const mesh = new Mesh({triangles: true, lines: false, normals: true})
		const matrix = M4.forSys(this.right, this.up, this.plane.normal1, this.plane.anchor)
		mesh.vertices = [V(xMin, yMin), V(xMax, yMin), V(xMin, yMax), V(xMax, yMax)].map(p => matrix.transformPoint(p))
		mesh.normals = arrayFromFunction(4, i => this.plane.normal1)
		pushQuad(mesh.triangles, false, 0, 1, 2, 3)
		mesh.compile()
		return mesh
	}


    dpds(): (s: number, t: number) => V3 {
        return () => this.right
    }

    dpdt(): (s: number, t: number) => V3 {
        return () => this.up
    }

    equals(obj: any): boolean {
        return null
    }

    didp(pWC: V3): V3 {
        return this.plane.normal1
    }

    static throughPoints(a: V3, b: V3, c: V3): PlaneSurface {
		return new PlaneSurface(P3.throughPoints(a, b, c))
	}
}