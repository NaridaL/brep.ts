import { assert, assertInst, DEG, eq0, fuzzyBetween, hasConstructor, lerp, M4, V3 } from 'ts3dutils'

import {
	Curve,
	Edge,
	HyperbolaCurve,
	intersectionUnitCircleLine2,
	L3,
	P3,
	ParametricSurface,
	PlaneSurface,
	PointVsFace,
	SemiEllipseCurve,
	Surface,
} from '../index'

import { abs, cos, PI, sin } from '../math'

/**
 * Rotation surface with r = f(z)
 */
export class RotatedCurveSurface extends ParametricSurface {
	readonly matrixInverse: M4

	constructor(
		readonly curve: Curve,
		readonly matrix: M4 = M4.IDENTITY,
		sMin: number = 0,
		sMax: number = PI,
		tMin: number = curve.tMin,
		tMax: number = curve.tMax,
	) {
		// d/dz (r(z))
		super(sMin, sMax, tMin, tMax)
		assertInst(M4, matrix)
		assert(matrix.isNoProj())
		assert(eq0(curve.at(tMin).y))
		this.matrixInverse = matrix.inversed()
	}

	getConstructorParameters(): any[] {
		return [this.curve, this.matrix, this.sMin, this.sMax, this.tMin, this.tMax]
	}

	flipped(): this {
		return new RotatedCurveSurface(
			this.curve,
			this.matrix.times(M4.mirror(P3.YZ)),
			this.sMin,
			this.sMax,
			this.tMin,
			this.tMax,
		) as this
	}

	transform(m4: M4): this {
		return new RotatedCurveSurface(
			this.curve,
			m4.isMirroring() ? m4.times(this.matrix).times(M4.mirror(P3.YZ)) : m4.times(this.matrix),
			this.sMin,
			this.sMax,
			this.tMin,
			this.tMax,
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
			const angle = SemiEllipseCurve.XYLCPointT(pLC, this.sMin, this.sMax)
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
		const dirLC = this.matrixInverse.transformVector(line.dir1)
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
			const secs = this.isCurvesWithPlaneLC(new P3(V3.Z, anchorLC.z))
			if (!secs) return []
			return secs.flatMap(sec => sec.isInfosWithLine(anchorLC, dirLC).map(info => info.tOther))
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

			const infos = hc.isInfosWithCurve(this.curve)
			return infos
				.map(info => (info.p.z - anchorLC.z) / dirLC.z)
				.filter(t => fuzzyBetween(L3.at(anchorLC, dirLC, t).angleXY(), this.sMin, this.sMax))
		}
	}

	private isCurvesWithPlaneLC(planeLC: P3): Curve[] | undefined {
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
		}
		return undefined
	}

	isCurvesWithPlane(plane: P3): Curve[] {
		const planeLC = plane.transform(this.matrixInverse)
		const planeLCCurves = this.isCurvesWithPlaneLC(planeLC)
		if (planeLCCurves) {
			return planeLCCurves.map(curve => curve.transform(this.matrix))
		} else {
			return ParametricSurface.isCurvesParametricImplicitSurface(this, new PlaneSurface(plane), 0.05, 0.05, 0.02)
		}
	}

	loopContainsPoint(loop: Edge[], pWC: V3): PointVsFace {
		const pLC = this.matrixInverse.transformPoint(pWC)
		const angle = SemiEllipseCurve.XYLCPointT(pLC, this.sMin, this.sMax)
		const testCurveLC = SemiEllipseCurve.semicircle(pLC.lengthXY(), new V3(0, 0, pLC.z))
		const testCurveWC = testCurveLC.transform(this.matrix)
		return Surface.loopContainsPointEllipse(loop, pWC, testCurveWC, angle)
		throw new Error('Method not implemented.')
	}

	isCoplanarTo(surface: Surface): boolean {
		if (this === surface) return true
		if (!hasConstructor(surface, RotatedCurveSurface)) return false
		const surfaceLCToThisLC = this.matrixInverse.times(surface.matrix)
		assert(!surfaceLCToThisLC.X.xy().likeO())
		const zRotation = surfaceLCToThisLC.X.angleXY()
		return surface.curve.transform(M4.rotateZ(-zRotation).times(surfaceLCToThisLC)).isColinearTo(this.curve)
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

	getExtremePoints(): V3[] {
		// this logic comes from EllipseCurve.roots
		const f1 = this.matrix.X
		const f2 = this.matrix.Y
		return [0, 1, 2].flatMap(dim => {
			const a = f2.e(dim),
				b = -f1.e(dim)
			const xiEtas = eq0(a) && eq0(b) ? [[1, 0]] : intersectionUnitCircleLine2(a, b, 0)
			return xiEtas.flatMap(([xi, eta]) => {
				const s = Math.atan2(eta, xi)
				if (!fuzzyBetween(s, this.sMin, this.sMax)) return []
				const testCurve = this.curve.transform(this.matrix.times(M4.rotateZ(s)))
				return testCurve.roots()[dim].map(t => this.pST(s, t))
			})
		})
	}
}

RotatedCurveSurface.prototype.uStep = PI / 16
RotatedCurveSurface.prototype.vStep = 1 / 4
