/**
 * @prettier
 */
import { assert, assertInst, eq0, M4, V3, lerp, hasConstructor, le, fuzzyBetween, DEG } from 'ts3dutils'

import {
	ParametricSurface,
	Curve,
	L3,
	P3,
	Edge,
	PointVsFace,
	Surface,
	SemiEllipseCurve,
	PlaneSurface,
	HyperbolaCurve,
} from '../index'

import { PI, cos, sin, abs } from '../math'

/**
 * Rotation surface with r = f(z)
 */
export class RotatedCurveSurface extends ParametricSurface {
	readonly matrixInverse: M4

	constructor(
		readonly curve: Curve,
		readonly tMin: number = curve.tMin,
		readonly tMax: number = curve.tMax,
		readonly matrix: M4 = M4.IDENTITY,
	) {
		// d/dz (r(z))
		super(0, PI, tMin, tMax)
		assertInst(M4, matrix)
		assert(matrix.isNoProj())
		assert(eq0(curve.at(tMin).y))
		this.matrixInverse = matrix.inversed()
	}

	getConstructorParameters(): any[] {
		return [this.curve, this.tMin, this.tMax, this.matrix]
	}

	flipped(): this {
		return new RotatedCurveSurface(this.curve, this.tMin, this.tMax, this.matrix.times(M4.mirror(P3.YZ))) as this
	}

	transform(m4: M4): this {
		return new RotatedCurveSurface(
			this.curve,
			this.tMin,
			this.tMax,
			m4.isMirroring() ? m4.times(this.matrix).times(M4.mirror(P3.YZ)) : m4.times(this.matrix),
		) as this
	}

	containsPoint(pWC: V3): boolean {
		const pLC = this.matrixInverse.transformPoint(pWC)
		const radius = pLC.lengthXY()
		return this.curve.containsPoint(new V3(radius, 0, pLC.z))
	}

	pSTFunc(): (s: number, t: number) => V3 {
		return (s, t) => {
			const { x: radius, z: z } = this.curve.at(t)
			return this.matrix.transformPoint(V3.polar(radius, s, z))
		}
	}

	dpds(): (s: number, t: number) => V3 {
		return (s, t) => {
			const radius = this.curve.at(t).x
			const resultLC = new V3(radius * -sin(s), radius * cos(s), 0)
			return this.matrix.transformVector(resultLC)
		}
	}

	dpdt(): (s: number, t: number) => V3 {
		return (s, t) => {
			const { x: drdt, z: dzdt } = this.curve.tangentAt(t)
			return this.matrix.transformVector(V3.polar(drdt, s, dzdt))
		}
	}

	normalSTFunc(): (s: number, t: number) => V3 {
		const matrix = this.matrix
			.inversed()
			.transposed()
			.as3x3()
		const normalLength = this.matrix.isMirroring() ? -1 : 1
		return (s, t) => {
			const { x: drdt, z: dzdt } = this.curve.tangentAt(t)
			return matrix.transformVector(V3.polar(dzdt, s, -drdt)).toLength(normalLength)
		}
	}

	stPFunc(): (pWC: V3) => V3 {
		return pWC => {
			const pLC = this.matrixInverse.transformPoint(pWC)
			const angle = abs(pLC.angleXY())
			const radius = pLC.lengthXY()
			return new V3(angle, this.curve.pointT(new V3(radius, 0, pLC.z)), 0)
		}
	}

	pointFoot(pWC: V3, startS?: number, startT?: number): V3 {
		const pLC = this.matrixInverse.transformPoint(pWC)
		const angle = abs(pLC.angleXY())
		const radius = pLC.lengthXY()
		return new V3(angle, this.curve.closestTToPoint(new V3(radius, 0, pLC.z)), 0)
	}

	isTsForLine(line: L3): number[] {
		const anchorLC = this.matrixInverse.transformPoint(line.anchor)
		const dirLC = this.matrixInverse.transformPoint(line.dir1)
		if (dirLC.isParallelTo(V3.Z)) {
			if (!fuzzyBetween(anchorLC.angleXY(), this.sMin, this.sMax)) return []
			return this.curve
				.isInfosWithLine(new V3(anchorLC.lengthXY(), 0, anchorLC.z), dirLC)
				.map(info => info.tOther)
		} else if (L3.containsPoint(anchorLC.xy(), dirLC.xy(), V3.O)) {
			// line goes through Z axis
			const dotter = dirLC.xy().unit()
			return [
				...this.curve.isInfosWithLine(
					new V3(dotter.dot(anchorLC), 0, anchorLC.z),
					new V3(dotter.dot(dirLC), 0, dirLC.z),
				),
				...this.curve.isInfosWithLine(
					new V3(-dotter.dot(anchorLC), 0, anchorLC.z),
					new V3(-dotter.dot(dirLC), 0, dirLC.z),
				),
			]
				.map(info => info.tOther)
				.filter(t => fuzzyBetween(L3.at(anchorLC, dirLC, t).angleXY(), this.sMin, this.sMax))
		} else if (dirLC.isPerpendicularTo(V3.Z)) {
			const sec = this.isCurvesWithPlane(new P3(V3.Z, anchorLC.z))[0] as SemiEllipseCurve | undefined
			if (!sec) return []
			assertInst(SemiEllipseCurve, sec)
			return sec.isInfosWithLine(anchorLC, dirLC).map(info => info.tOther)
		} else {
			// transform into hyperbola
			// f(t) = V(((ax + t dx)² + (ay + t dy)²) ** 1/2, 0, az + t dz)
			// f(t) = V((ax² + 2 ax t dx + t² dx² + ay² + 2 ay t dy + t² dy²) ** 1/2, 0, az + t dz)
			// f(t) = V((t² (dx² + dy²) + 2 t (ax dx + ay dy) + ax² + ay²) ** 1/2, 0, az + t * dz)

			// (anchorLC.xy + t * dirLC.xy) * dir.xy = 0
			// t * dirLC.xy² = -anchorLC.xy * dirLC.xy
			const closestTToZ = -anchorLC.xy().dot(dirLC.xy()) / dirLC.xy().squared()
			const closestPointToZ = L3.at(anchorLC, dirLC, closestTToZ)
			const scaleX = closestPointToZ.lengthXY()
			const lineGradientWC = dirLC.z / dirLC.lengthXY()
			const scaleZ = scaleX * lineGradientWC
			const hc = HyperbolaCurve.XY.transform(
				M4.rotateX(90 * DEG)
					.scale(scaleX, 0, scaleZ)
					.translate(0, 0, closestPointToZ.z),
			)
			console.log(hc.sce, closestPointToZ)

			const infos = hc.isInfosWithCurve(this.curve)
			return infos
				.map(info => (info.p.z - anchorLC.z) / dirLC.z)
				.filter(t => fuzzyBetween(L3.at(anchorLC, dirLC, t).angleXY(), this.sMin, this.sMax))
		}
	}

	isCurvesWithPlane(plane: P3): Curve[] {
		const planeLC = plane.transform(this.matrixInverse)
		if (planeLC.normal1.isParallelTo(V3.Z)) {
			return this.curve.isTsWithPlane(planeLC).map(t => {
				const { x: radius } = this.curve.at(t)
				return new SemiEllipseCurve(
					new V3(0, 0, planeLC.w),
					new V3(radius, 0, 0),
					new V3(0, radius, 0),
					this.sMin,
					this.sMax,
				).transform(this.matrix)
			})
		} else if (planeLC.normal1.isPerpendicularTo(V3.Z) && planeLC.containsPoint(V3.O)) {
			return [this.curve.rotateZ(V3.Y.angleRelativeNormal(planeLC.normal1, V3.Z)).transform(this.matrix)]
		} else {
			return ParametricSurface.isCurvesParametricImplicitSurface(this, new PlaneSurface(plane), 0.05, 0.05, 0.02)
		}
	}

	loopContainsPoint(contour: Edge[], point: V3): PointVsFace {
		throw new Error('Method not implemented.')
	}

	isCoplanarTo(surface: Surface): boolean {
		return (
			this == surface ||
			(hasConstructor(surface, RotatedCurveSurface) &&
				surface.curve.transform(this.matrixInverse).isColinearTo(this.curve) &&
				this.matrixInverse.times(surface.matrix).isZRotation())
		)
	}

	like(object: any): boolean {
		if (!this.isCoplanarTo(object)) return false
		// normals need to point in the same direction (outwards or inwards) for both
		const pSMinTMin = this.pSTFunc()(this.sMin, this.tMin)
		const thisNormal = this.normalSTFunc()(this.sMin, this.tMin)
		const otherNormal = object.normalP(pSMinTMin)
		return 0 < thisNormal.dot(otherNormal)
	}

	equals(obj: any): boolean {
		return (
			this == obj ||
			(hasConstructor(obj, RotatedCurveSurface) && this.curve.equals(obj.curve) && this.matrix.equals(obj.matrix))
		)
	}

	isCurvesWithSurface(surface: Surface): Curve[] {
		if (surface instanceof PlaneSurface) {
			return this.isCurvesWithPlane(surface.plane)
		}
		return super.isCurvesWithSurface(surface)
	}

	containsCurve(curve: Curve): boolean {
		if (curve.constructor == this.curve.constructor) {
			const curveLC = curve.transform(this.matrixInverse)
			// find a point on curveLC which isn't on the Z-axis
			const t = [0, 0.5, 1].map(x => lerp(curveLC.tMin, curveLC.tMax, x)).withMax(t => curveLC.at(t).lengthXY())
			const angle = curveLC.at(t).angleXY()
			const curveLCRotated = curveLC.rotateZ(-angle)
			if (this.curve.isColinearTo(curveLCRotated)) {
				return true
			}
		}
		if (curve instanceof SemiEllipseCurve) {
			const curveLC = curve.transform(this.matrixInverse)
			if (curveLC.normal.isParallelTo(V3.Z)) {
				return (
					curveLC.isCircular() && this.curve.containsPoint(new V3(curveLC.f1.length(), 0, curveLC.center.z))
				)
			}
			return false
		}
		return super.containsCurve(curve)
	}

	hashCode() {
		return [this.curve, this.matrix].hashCode()
	}
}

RotatedCurveSurface.prototype.uStep = PI / 16
RotatedCurveSurface.prototype.vStep = 1 / 4
