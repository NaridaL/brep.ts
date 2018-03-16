import {arrayFromFunction, assert, assertInst, isCCW, M4, V, V3, callsce} from 'ts3dutils'
import {Mesh, pushQuad} from 'tsgl'

import {Curve, Edge, ImplicitSurface, L3, P3, ParametricSurface, PointVsFace, Surface, ImplicitCurve,} from '../index'

export class PlaneSurface extends ParametricSurface implements ImplicitSurface {
	readonly matrix: M4

	constructor(readonly plane: P3,
		readonly right: V3 = plane.normal1.getPerpendicular().unit(),
		readonly up: V3 = plane.normal1.cross(right).unit(),
		sMin: number = -100,
		sMax: number = 100,
		tMin: number = -100,
		tMax: number = 100) {
		super(sMin, sMax, tMin, tMax)
		assertInst(P3, plane)
		assert(this.right.cross(this.up).like(this.plane.normal1))
		this.matrix = M4.forSys(right, up, plane.normal1, plane.anchor)
	}

	toSource(rounder: (x: number) => number = x => x): string {
		return callsce.call(undefined, 'new PlaneSurface', ...this.getConstructorParameters())
	}

	static throughPoints(a: V3, b: V3, c: V3): PlaneSurface {
		return new PlaneSurface(P3.throughPoints(a, b, c))
	}

	static forAnchorAndPlaneVectors(anchor: V3, v0: V3, v1: V3, sMin?: number, sMax?: number, tMin?: number, tMax?: number): PlaneSurface {
		return new PlaneSurface(
			P3.forAnchorAndPlaneVectors(anchor, v0, v1),
			v0, v1,
			sMin, sMax,
			tMin, tMax)
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
	    assert(Edge.isLoop(contour), 'isLoop')
		return isCCW(contour.flatMap(edge => edge.points()), this.plane.normal1)
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

	containsPoint(p: V3) {
		return this.plane.containsPoint(p)
	}

	containsCurve(curve: Curve): boolean {
		return curve instanceof ImplicitCurve ? super.containsCurve(curve) : this.plane.containsCurve(curve)
	}

	transform(m4: M4) {
		return new PlaneSurface(this.plane.transform(m4)) as this
	}

	flipped() {
		return new PlaneSurface(this.plane.flipped(), this.right, this.up.negated()) as this
	}

	getConstructorParameters(): any[] {
		return [this.plane, this.right, this.up]
	}

	// toMesh(xMin: number = -10, xMax: number = 10, yMin: number = -10, yMax: number = 10) {
	// 	const mesh = new Mesh()
	// 		.addIndexBuffer('TRIANGLES')
	// 		.addVertexBuffer('normals', 'ts_Normal')
	// 	const matrix = M4.forSys(this.right, this.up, this.plane.normal1, this.plane.anchor)
	// 	mesh.vertices = [V(xMin, yMin), V(xMax, yMin), V(xMin, yMax), V(xMax, yMax)].map(p => matrix.transformPoint(p))
	// 	mesh.normals = arrayFromFunction(4, i => this.plane.normal1)
	// 	pushQuad(mesh.TRIANGLES, false, 0, 1, 2, 3)
	// 	mesh.compile()
	// 	return mesh
	// }

	dpds(): (s: number, t: number) => V3 {
		return () => this.right
	}

	dpdt(): (s: number, t: number) => V3 {
		return () => this.up
	}

	equals(obj: any): boolean {
		return undefined
	}

	didp(pWC: V3): V3 {
		return this.plane.normal1
	}

	normalST() {
	    return this.plane.normal1
    }
}

PlaneSurface.prototype.uStep = 1e6
PlaneSurface.prototype.vStep = 1e6
