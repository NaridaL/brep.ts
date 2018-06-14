import { assert, assertInst, callsce, hasConstructor, isCCW, M4, V3 } from 'ts3dutils'

import { Curve, Edge, ImplicitCurve, ImplicitSurface, L3, P3, ParametricSurface, PointVsFace, Surface } from '../index'

export class PlaneSurface extends ParametricSurface implements ImplicitSurface {
	readonly matrix: M4

	constructor(
		readonly plane: P3,
		readonly right: V3 = plane.normal1.getPerpendicular().unit(),
		readonly up: V3 = plane.normal1.cross(right).unit(),
		uMin: number = -100,
		uMax: number = 100,
		vMin: number = -100,
		vMax: number = 100,
	) {
		super(uMin, uMax, vMin, vMax)
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

	static forAnchorAndPlaneVectors(
		anchor: V3,
		v0: V3,
		v1: V3,
		uMin?: number,
		uMax?: number,
		vMin?: number,
		vMax?: number,
	): PlaneSurface {
		return new PlaneSurface(P3.forAnchorAndPlaneVectors(anchor, v0, v1), v0, v1, uMin, uMax, vMin, vMax)
	}
	isCoplanarTo(surface: Surface): boolean {
		return hasConstructor(surface, PlaneSurface) && this.plane.isCoplanarToPlane(surface.plane)
	}

	isTsForLine(line: L3): number[] {
		return line.isTsWithPlane(this.plane)
	}

	like(surface: Surface): boolean {
		return hasConstructor(surface, PlaneSurface) && this.plane.like(surface.plane)
	}

	pUV(u: number, v: number): V3 {
		return this.matrix.transformPoint(new V3(u, v, 0))
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
		const result = this.plane.intersectionWithPlane(plane)
		return result ? [result] : []
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

	uvPFunc() {
		const matrixInverse = this.matrix.inversed()
		return function(pWC: V3) {
			return matrixInverse.transformPoint(pWC)
		}
	}

	pointFoot(pWC: V3): V3 {
		return this.uvP(pWC)
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
		return [this.plane, this.right, this.up, this.uMin, this.uMax, this.vMin, this.vMax]
	}

	dpdu(): (u: number, v: number) => V3 {
		return () => this.right
	}

	dpdv(): (u: number, v: number) => V3 {
		return () => this.up
	}

	didp(pWC: V3): V3 {
		return this.plane.normal1
	}

	normalUV() {
		return this.plane.normal1
	}
}

PlaneSurface.prototype.uStep = 1e6
PlaneSurface.prototype.vStep = 1e6
