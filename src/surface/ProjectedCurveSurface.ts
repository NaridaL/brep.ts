import { assert, assertInst, assertNumbers, assertVectors, hasConstructor, int, M4, V3 } from 'ts3dutils'

import {
	Curve,
	Edge,
	EllipsoidSurface,
	ImplicitCurve,
	ImplicitSurface,
	L3,
	P3,
	ParametricSurface,
	PlaneSurface,
	PointVsFace,
	Surface,
} from '../index'

import { sign } from '../math'

/**
 * Surface normal1 is (t, z) => this.baseCurve.tangentAt(t) X this.dir
 * Choose dir appropriately to select surface orientation.
 */
export class ProjectedCurveSurface extends ParametricSurface {
	readonly ['constructor']: typeof ProjectedCurveSurface & {
		new <T extends ProjectedCurveSurface>(
			baseCurve: Curve,
			dir: V3,
			uMin: number,
			uMax: number,
			vMin: number,
			vMax: number,
		): T
	}

	constructor(
		readonly baseCurve: Curve,
		readonly dir: V3,
		uMin: number = baseCurve.tMin,
		uMax: number = baseCurve.tMax,
		vMin: number = -100,
		vMax: number = 100,
	) {
		super(uMin, uMax, vMin, vMax)
		assertInst(Curve, baseCurve)
		assertInst(V3, dir)
		assert(uMin < uMax)
		assert(vMin < vMax)
	}

	getConstructorParameters() {
		return [this.baseCurve, this.dir, this.uMin, this.uMax, this.vMin, this.vMax]
	}

	equals(obj: any): boolean {
		return (
			this == obj ||
			(Object.getPrototypeOf(this) == Object.getPrototypeOf(obj) &&
				this.dir.equals(obj.dir) &&
				this.baseCurve.equals(obj.baseCurve))
		)
	}

	hashCode(): int {
		return [this.dir, this.baseCurve].hashCode()
	}

	containsLine(line: L3): boolean {
		return this.dir.isParallelTo(line.dir1) && this.containsPoint(line.anchor)
	}

	dpdu(): (u: number, v: number) => V3 {
		return (u, v) => this.baseCurve.tangentAt(u)
	}

	dpdv(): (u: number, v: number) => V3 {
		return (u, v) => this.dir
	}

	normalUV(u: number, v: number): V3 {
		return this.baseCurve
			.tangentAt(u)
			.cross(this.dir)
			.unit()
	}

	pUV(u: number, v: number): V3 {
		return this.baseCurve.at(u).plus(this.dir.times(v))
	}

	pointFoot(pWC: V3, ss?: number): V3 {
		const basePlane = new P3(this.dir.unit(), 0)
		const projCurve = this.baseCurve.project(basePlane)
		const projPoint = basePlane.projectedPoint(pWC)
		const t = projCurve.closestTToPoint(projPoint, ss, this.uMin, this.uMax)
		const z = L3.pointT(this.baseCurve.at(t), this.dir, pWC)
		return new V3(t, z, 0)
	}

	uvPFunc(): (pWC: V3) => V3 {
		const projPlane = new P3(this.dir.unit(), 0)
		const projBaseCurve = this.baseCurve.project(projPlane)
		return pWC => {
			const projPoint = projPlane.projectedPoint(pWC)
			assertNumbers(this.uMin)
			const t = projBaseCurve.pointT(projPoint, this.uMin, this.uMax)
			const z = L3.pointT(this.baseCurve.at(t), this.dir, pWC)
			return new V3(t, z, 0)
		}
	}

	isCurvesWithPlane(plane: P3): Curve[] {
		assertInst(P3, plane)
		if (this.dir.isPerpendicularTo(plane.normal1)) {
			const ts = this.baseCurve.isTsWithPlane(plane)
			return ts.map(t => {
				const l3dir = 0 < this.baseCurve.tangentAt(t).dot(plane.normal1) ? this.dir : this.dir.negated()
				return new L3(this.baseCurve.at(t), l3dir.unit())
			})
		} else {
			let projCurve = this.baseCurve.transform(M4.project(plane, this.dir))
			if (this.dir.dot(plane.normal1) > 0) {
				// we need to flip the ellipse so the tangent is correct
				projCurve = projCurve.reversed()
			}
			return [projCurve]
		}
	}

	isCurvesWithSurface(surface: Surface): Curve[] {
		if (surface instanceof PlaneSurface) {
			return this.isCurvesWithPlane(surface.plane)
		}
		if (surface instanceof ProjectedCurveSurface) {
			const dir1 = surface.dir
			if (this.dir.isParallelTo(dir1)) {
				const ts = surface.baseCurve.isTsWithSurface(this)
				return ts.map(t => {
					const p = surface.baseCurve.at(t)
					const correctDir = this.normalP(p).cross(surface.normalP(p))
					return new L3(p, dir1.times(sign(correctDir.dot(dir1))))
				})
			} else if (ImplicitSurface.is(surface)) {
				let curves2 = ParametricSurface.isCurvesParametricImplicitSurface(
					this,
					surface,
					0.1,
					0.1 / surface.dir.length(),
					0.05,
				)
				curves2 = surface.clipCurves(curves2)
				return curves2
			} else {
				let curves2 = ParametricSurface.isCurvesParametricParametricSurface(
					this,
					surface,
					0.05,
					0.1 / surface.dir.length(),
					0.05,
				)
				curves2 = this.clipCurves(curves2)
				curves2 = surface.clipCurves(curves2)
				return curves2
			}
		}
		if (surface instanceof EllipsoidSurface) {
			return surface.isCurvesWithSurface(this)
		}
		return super.isCurvesWithSurface(surface)
	}

	containsPoint(pWC: V3): boolean {
		const uv = this.uvPFunc()(pWC)
		return this.pUVFunc()(uv.x, uv.y).like(pWC)
	}

	containsCurve(curve: Curve): boolean {
		if (curve instanceof L3) {
			return this.dir.isParallelTo(curve.dir1) && this.containsPoint(curve.anchor)
		}
		if (curve instanceof ImplicitCurve) {
			return super.containsCurve(curve)
		}
		// project baseCurve and test curve onto a common plane and check if the curves are alike
		const projPlane = new P3(this.dir.unit(), 0)
		const projBaseCurve = this.baseCurve.project(projPlane)
		const projCurve = curve.project(projPlane)

		return projBaseCurve.isColinearTo(projCurve)
	}

	isCoplanarTo(surface: Surface): boolean {
		return (
			this == surface ||
			(hasConstructor(surface, ProjectedCurveSurface) &&
				this.dir.isParallelTo(surface.dir) &&
				this.containsCurve(surface.baseCurve))
		)
	}

	like(object: any): boolean {
		if (!this.isCoplanarTo(object)) return false
		// normals need to point in the same direction (outwards or inwards) for both
		const p00 = this.pUVFunc()(0, 0)
		const thisNormal = this.normalUVFunc()(0, 0)
		const otherNormal = object.normalP(p00)
		return 0 < thisNormal.dot(otherNormal)
	}

	loopContainsPoint(loop: Edge[], p: V3): PointVsFace {
		assertVectors(p)
		assert(isFinite(p.x), p.y, p.z)
		const line = new L3(p, this.dir.unit())
		const ptpf = this.uvPFunc()
		const pp = ptpf(p)
		if (isNaN(pp.x)) {
			console.log(this.sce, p.sce)
			assert(false)
		}
		const lineOut = this.baseCurve.tangentAt(pp.x).rejectedFrom(this.dir)

		return Surface.loopContainsPointGeneral(loop, p, line, lineOut)
	}

	transform(m4: M4): this {
		const f = m4.isMirroring() ? -1 : 1
		return new this.constructor<this>(
			this.baseCurve.transform(m4),
			m4.transformVector(this.dir).times(f),
			this.uMin,
			this.uMax,
			1 == f ? this.vMin : -this.vMax,
			1 == f ? this.vMax : -this.vMin,
		)
	}

	isTsForLine(line: L3): number[] {
		assertInst(L3, line)
		const projPlane = new P3(this.dir.unit(), 0)
		const projDir = projPlane.projectedVector(line.dir1)
		if (projDir.likeO()) {
			// line is parallel to this.dir
			return []
		}
		const projAnchor = projPlane.projectedPoint(line.anchor)
		const projBaseCurve = this.baseCurve.project(projPlane)
		return projBaseCurve
			.isInfosWithLine(projAnchor, projDir, this.uMin, this.uMax, line.tMin, line.tMax)
			.map(info => info.tOther)
	}

	flipped(): this {
		return new this.constructor<this>(
			this.baseCurve,
			this.dir.negated(),
			this.uMin,
			this.uMax,
			-this.vMax,
			-this.vMin,
		)
	}
}
ProjectedCurveSurface.prototype.uStep = 1 / 128
ProjectedCurveSurface.prototype.vStep = 256
