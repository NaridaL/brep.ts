import { assert, assertInst, assertVectors, hasConstructor, M4, newtonIterate, pqFormula, V3 } from 'ts3dutils'

import {
	Curve,
	Edge,
	EllipseCurve,
	ImplicitCurve,
	ImplicitSurface,
	L3,
	P3,
	ParametricSurface,
	PlaneSurface,
	PointVsFace,
	Surface,
} from '../index'

export class PointProjectedSurface extends ParametricSurface {
	pointFoot(pWC: V3, ss?: number, st?: number): V3 {
		if (undefined === ss || undefined === st) {
			// similar to stP
			if (undefined === ss) {
				ss = pWC.like(this.apex)
					? 0
					: this.curve.closestTToPoint(this.planeProjectionMatrix.transformPoint(pWC)) * this.normalDir
			}
			if (undefined === st) {
				st = V3.inverseLerp(this.apex, this.curve.at(ss), pWC)
			}
		}
		const f = ([s, t]: number[]) => {
			const pSTToPWC = this.pST(s, t).to(pWC)
			return [this.dpds()(s, t).dot(pSTToPWC), this.dpdt()(s).dot(pSTToPWC)]
		}
		const { 0: x, 1: y } = newtonIterate(f, [ss, st])
		return new V3(x, y, 0)
	}

	readonly planeProjectionMatrix: M4

	constructor(
		readonly curve: Curve,
		readonly apex: V3,
		readonly curvePlane: P3,
		readonly normalDir = 1,
		sMin: number = curve.tMin,
		sMax: number = curve.tMax,
		tMin: number = 0,
		tMax: number = 16,
	) {
		super(sMin, sMax, tMin, tMax)
		assertInst(Curve, curve)
		assert(!(curve instanceof L3), 'use PlaneSurface instead')
		assert(!(curve instanceof EllipseCurve), 'use ConicSurface instead')
		assert(!(curve instanceof ImplicitCurve), 'this just seems like a terrible idea')
		assert(new PlaneSurface(curvePlane).containsCurve(curve))
		assertVectors(apex)
		assert(0 <= tMin)
		this.planeProjectionMatrix = M4.projectPlanePoint(apex, curvePlane)
		this.uStep = curve.tIncrement
	}

	getConstructorParameters(): any[] {
		return [this.curve, this.apex, this.curvePlane, this.normalDir, this.sMin, this.sMax, this.tMin, this.tMax]
	}

	static unitISLineTs(anchor: V3, dir: V3): number[] {
		const { x: ax, y: ay, z: az } = anchor
		const { x: dx, y: dy, z: dz } = dir

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

	equals(obj: any): boolean {
		return (
			this == obj ||
			(hasConstructor(obj, PointProjectedSurface) && this.curve.equals(obj.curve) && this.apex.equals(this.apex))
		)
	}

	like(object: any): boolean {
		if (!this.isCoplanarTo(object)) return false
		// normals need to point in the same direction (outwards or inwards) for both
		return this.normalDir == object.normalDir
	}

	loopContainsPoint(contour: Edge[], p: V3): PointVsFace {
		assertVectors(p)
		const line = this.apex.like(p)
			? new L3(p, this.apex.to(this.curve.at(this.curve.tMin)).unit())
			: L3.throughPoints(p, this.apex)
		const lineOut = line.dir1.cross(this.curvePlane.normal1)

		return Surface.loopContainsPointGeneral(contour, p, line, lineOut)
	}

	isTsForLine(line: L3): number[] {
		// transforming line manually has advantage that dir1 will not be renormalized,
		// meaning that calculated values t for lineLC are directly transferable to line
		const anchorPlane = this.planeProjectionMatrix.transformPoint(line.anchor)
		const anchor2Plane = this.planeProjectionMatrix.transformPoint(line.anchor.plus(line.dir1))
		if (anchorPlane.like(anchor2Plane)) {
			// line projects onto a point in plane.
			// there are either no or infinite intersection points
			return []
		}
		return this.curve
			.isInfosWithLine(anchorPlane, anchorPlane.to(anchor2Plane), undefined, undefined, line.tMin, line.tMax)
			.map(info => info.tOther)
	}

	/**
	 * Interestingly, two cones don't need to have parallel dirs to be coplanar.
	 */
	isCoplanarTo(surface: Surface): boolean {
		if (this === surface) return true
		if (!(surface instanceof PointProjectedSurface) || !this.apex.like(surface.apex)) return false
		// at this point apexes are equal
		return this.containsCurve(surface.curve)
	}

	containsLine(line: L3): boolean {
		if (this.curvePlane.isParallelToLine(line)) {
			return false
		}
		if (!line.containsPoint(this.apex)) {
			return false
		}
		const p = this.curvePlane.intersectionWithLine(line)
		return this.curve.containsPoint(p)
	}

	containsCurve(curve: Curve): boolean {
		if (curve instanceof L3) {
			return this.containsLine(curve)
		} else if (!(curve instanceof ImplicitCurve)) {
			const otherCurveOnThisPlane = curve.transform(this.planeProjectionMatrix)
			return this.curve.isColinearTo(otherCurveOnThisPlane)
		} else {
			return super.containsCurve(curve)
		}
	}

	transform(m4: M4): this {
		return new PointProjectedSurface(
			this.curve.transform(m4),
			m4.transformPoint(this.apex),
			this.curvePlane.transform(m4),
			(m4.isMirroring() ? -1 : 1) * this.normalDir,
			this.sMin,
			this.sMax,
			this.tMin,
			this.tMax,
		) as this
	}

	flipped(): this {
		return new PointProjectedSurface(
			this.curve,
			this.apex,
			this.curvePlane,
			-this.normalDir,
			-this.sMax,
			-this.sMin,
			this.tMin,
			this.tMax,
		) as this
	}

	normalSTFunc(): (s: number, t: number) => V3 {
		const dpdt = this.dpdt()
		return (s, t) =>
			this.curve
				.tangentAt(s * this.normalDir)
				.times(this.normalDir)
				.cross(dpdt(s))
				.unit()
	}

	pSTFunc(): (s: number, t: number) => V3 {
		return (s, t) => {
			return this.apex.lerp(this.curve.at(s * this.normalDir), t)
		}
	}

	dpds(): (s: number, t: number) => V3 {
		return (s, t) => {
			return this.curve.tangentAt(s * this.normalDir).times(t * this.normalDir)
		}
	}

	dpdt(): (s: number) => V3 {
		return s => {
			return this.apex.to(this.curve.at(s * this.normalDir))
		}
	}

	containsPoint(pWC: V3) {
		return this.apex.like(pWC) || this.curve.containsPoint(this.planeProjectionMatrix.transformPoint(pWC))
	}

	stP(pWC: V3) {
		const s = pWC.like(this.apex) ? 0 : this.curve.pointT(this.planeProjectionMatrix.transformPoint(pWC))
		const t = V3.inverseLerp(this.apex, this.curve.at(s), pWC)
		return new V3(s * this.normalDir, t, 0)
	}

	isCurvesWithSurface(surface: Surface): Curve[] {
		if (surface instanceof PlaneSurface) {
			return this.isCurvesWithPlane(surface.plane)
		} else if (ImplicitSurface.is(surface)) {
			return ParametricSurface.isCurvesParametricImplicitSurface(
				this,
				surface,
				0.1,
				0.1 / this.curvePlane.distanceToPoint(this.apex),
				0.02,
			)
		}
		return super.isCurvesWithSurface(surface)
	}

	isCurvesWithPlane(plane: P3): Curve[] {
		if (plane.containsPoint(this.apex)) {
			if (plane.isParallelToPlane(this.curvePlane)) {
				return []
			}
			return this.curve.isTsWithPlane(plane).map(t => L3.throughPoints(this.apex, this.curve.at(t)))
		}
		return [this.curve.transform(M4.projectPlanePoint(this.apex, plane))]
	}
}

PointProjectedSurface.prototype.vStep = 256
