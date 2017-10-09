import {
	arrayFromFunction,
	assertf,
	assertInst,
	assertNumbers,
	assertVectors,
	eq0,
	hasConstructor,
	int,
	M4,
	NLA_PRECISION,
	TAU,
	V,
	V3,
} from 'ts3dutils'
import {Mesh, pushQuad} from 'tsgl'

import {BezierCurve, ConicSurface, Curve, EllipseCurve, ISInfo, L3, P3, ProjectedCurveSurface, Surface, EllipsoidSurface,
	SemiEllipsoidSurface,
	PlaneSurface,} from '../index'

const {PI} = Math

export abstract class XiEtaCurve extends Curve {
	readonly normal: V3
	readonly matrix: M4
	readonly inverseMatrix: M4
	'constructor': typeof XiEtaCurve & ( new(center: V3, f1: V3, f2: V3, tMin: number, tMax: number) => this )

	constructor(readonly center: V3,
				readonly f1: V3,
				readonly f2: V3,
				readonly tMin: number = -PI,
				readonly tMax: number = PI) {
		super(tMin, tMax)
		assertVectors(center, f1, f2)
		this.normal = f1.cross(f2)
		if (!this.normal.likeO()) {
			this.normal = this.normal.unit()
			this.matrix = M4.forSys(f1, f2, this.normal, center)
			this.inverseMatrix = this.matrix.inversed()
		} else {
			this.matrix = M4.forSys(f1, f2, f1.unit(), center)
			const f1p = f1.getPerpendicular()
			this.inverseMatrix = new M4(
				1, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 1).times(M4.forSys(f1, f1p, f1.cross(f1p), center).inversed())
		}
	}

	static magic(a: number, b: number, c: number): number[] {
		throw new Error('abstract')
	}

	/**
	 * Returns a new EllipseCurve representing an ellipse parallel to the XY-plane
	 * with semi-major/minor axes parallel t the X and Y axes and of length a and b.
	 *
	 * @param a length of the axis parallel to X axis
	 * @param b length of the axis parallel to Y axis
	 * @param center Defaults to V3.O
	 */
	static forAB(a: number, b: number, center: V3 = V3.O): XiEtaCurve {
		return new (this as any)(center, V(a, 0, 0), V(0, b, 0))
	}

	static XYLCValid(pLC: V3): boolean {
		throw new Error('abstract')
	}

	static XYLCPointT(pLC: V3): number {
		throw new Error('abstract')
	}

	static unitIsInfosWithLine(anchorLC: V3, dirLC: V3, anchorWC: V3, dirWC: V3): ISInfo[] {
		throw new Error('abstract')
	}

	addToMesh(mesh: Mesh & { TRIANGLES: int[], normals: V3[] }, res: int = 4, radius: number = 0, pointStep = 1): void {
		const baseNormals = arrayFromFunction(res, i => V3.polar(1, TAU * i / res))
		const baseVertices = arrayFromFunction(res, i => V3.polar(radius, TAU * i / res))
		const inc = this.tIncrement
		const start = Math.ceil((this.tMin + NLA_PRECISION) / inc)
		const end = Math.floor((this.tMax - NLA_PRECISION) / inc)
		for (let i = start; i <= end; i += pointStep) {
			const t = i * inc
			const start = mesh.vertices.length
			if (0 !== i) {
				for (let j = 0; j < res; j++) {
					pushQuad(mesh.TRIANGLES, true,
						start - res + j, start + j,
						start - res + (j + 1) % res, start + (j + 1) % res)
				}
			}
			const point = this.at(t), tangent = this.tangentAt(t)
			const matrix = M4.forSys(this.normal, tangent.cross(this.normal), tangent, point)
			mesh.normals.push(...matrix.transformedVectors(baseNormals))
			mesh.vertices.push(...matrix.transformedPoints(baseVertices))
		}
	}

	getConstructorParameters(): any[] {
		return [this.center, this.f1, this.f2, this.tMin, this.tMax]
	}

	isInfosWithCurve(curve: Curve): ISInfo[] {
		if (curve instanceof L3) {
			return this.isInfosWithLine(curve.anchor, curve.dir1, this.tMin, this.tMax, curve.tMin, curve.tMax)
		}
		if (curve instanceof BezierCurve) {
			return this.isInfosWithBezier(curve)
		}
		if (curve instanceof XiEtaCurve) {
			if (!this.normal.isParallelTo(curve.normal)) {
				return this.isTsWithPlane(curve.getPlane()).mapFilter(tThis => {
					const p = this.at(tThis)
					if (curve.containsPoint(p)) {
						return {tThis, tOther: curve.pointT(p), p}
					}
				})
			}
		}
		return super.isInfosWithCurve(curve)
	}

	transform(m4: M4) {
		return new this.constructor(
			m4.transformPoint(this.center),
			m4.transformVector(this.f1),
			m4.transformVector(this.f2),
			this.tMin, this.tMax) as this
	}

	equals(obj: any): boolean {
		return this == obj ||
			obj.constructor == this.constructor
			&& this.center.equals(obj.center)
			&& this.f1.equals(obj.f1)
			&& this.f2.equals(obj.f2)
	}

	hashCode(): int {
		let hashCode = 0
		hashCode = hashCode * 31 + this.center.hashCode()
		hashCode = hashCode * 31 + this.f1.hashCode()
		hashCode = hashCode * 31 + this.f2.hashCode()
		return hashCode | 0
	}

	likeCurve(curve: Curve): boolean {
		return hasConstructor(curve, this.constructor)
			&& this.center.like(curve.center)
			&& this.f1.like(curve.f1)
			&& this.f2.like(curve.f2)
	}

	normalP(t: number): V3 {
		return this.tangentAt(t).cross(this.normal)
	}

	getPlane(): P3 {
		return P3.normalOnAnchor(this.normal, this.center)
	}

	isTsWithPlane(plane: P3): number[] {
		assertInst(P3, plane)
		/*
		 this: x = center + f1 * cos t + f2 * sin t  (1)
		 plane:
		 n := plane.normal1
		 n DOT x == plane.w           (2)
		 plane defined by f1/f2
		 x = center + f1 * xi + f2 * eta         (3)
		 intersection plane and planef1/f2:
		 insert (3) into (2):
		 n DOT center + n DOT f1 * xi + n DOT f2 * eta = plane.w | -n DOT center
		 n DOT f1 * xi + n DOT f2 * eta = plane.w - n DOT center (4)
		 points on ellipse have additional condition
		 eta * eta + xi * xi = 1 (5)
		 g1 := n DOT f1
		 g2 := n DOT f2
		 g3 := w - n DOT center
		 solve system (5)/(6)
		 g1 * xi + g2 * eta = g3 (6)
		 */
		if (plane.normal1.isParallelTo(this.normal)) {
			return []
		}
		const n = plane.normal1, w = plane.w,
			center = this.center, f1 = this.f1, f2 = this.f2,
			g1 = n.dot(f1), g2 = n.dot(f2), g3 = w - n.dot(center)

		return this.constructor.magic(g1, g2, g3)
	}

	pointT(p: V3): number {
		assertVectors(p)
		const pLC = this.inverseMatrix.transformPoint(p)
		return this.constructor.XYLCPointT(pLC)
	}

	containsPoint(p: V3): boolean {
		const pLC = this.inverseMatrix.transformPoint(p)
		return eq0(pLC.z) && this.constructor.XYLCValid(pLC)
	}

	isInfosWithLine(anchorWC: V3, dirWC: V3, tMin?: number, tMax?: number, lineMin = -100000, lineMax = 100000): ISInfo[] {
		const anchorLC = this.inverseMatrix.transformPoint(anchorWC)
		const dirLC = this.inverseMatrix.transformVector(dirWC)
		if (eq0(dirLC.z)) {
			// local line parallel to XY-plane
			if (eq0(anchorLC.z)) {
				// local line lies in XY-plane
				return this.constructor.unitIsInfosWithLine(anchorLC, dirLC, anchorWC, dirWC)
			}
		} else {
			// if the line intersects the XY-plane in a single point, there can be an intersection there
			// find point, then check if distance from circle = 1
			const otherTAtZ0 = anchorLC.z / dirLC.z
			const isp = dirLC.times(otherTAtZ0).plus(anchorLC)
			if (this.constructor.XYLCValid(isp)) {
				// point lies on unit circle
				return [{
					tThis: this.constructor.XYLCPointT(isp),
					tOther: otherTAtZ0,
					p: anchorWC.plus(dirWC.times(otherTAtZ0)),
				}]
			}
		}
		return []
	}

	isTsWithSurface(surface: Surface): number[] {
		if (surface instanceof PlaneSurface) {
			return this.isTsWithPlane(surface.plane)
		} else if (surface instanceof SemiEllipsoidSurface) {
			const isEllipse = surface.asEllipsoidSurface().isCurvesWithSurface(new PlaneSurface(this.getPlane()))
			if (isEllipse.length < 1) return []
			const possibleInfos = this.isInfosWithCurve(isEllipse[0] as EllipseCurve)
			return possibleInfos.filter(info => surface.containsPoint(info.p)).map(info => info.tThis)
		} else if (surface instanceof ProjectedCurveSurface ||
			surface instanceof EllipsoidSurface ||
			surface instanceof ConicSurface) {
			return surface.isCurvesWithPlane(this.getPlane())
				.flatMap(curve => this.isInfosWithCurve(curve))
				.map(info => info.tThis)
		} else {
			throw new Error()
		}
	}

	isInfosWithBezier(bezierWC: BezierCurve): ISInfo[] {
		const bezierLC = bezierWC.transform(this.inverseMatrix)
		if (new PlaneSurface(P3.XY).containsCurve(bezierLC)) {
			return this.isInfosWithBezier2D(bezierWC)
		} else {
			const infos = bezierLC.isTsWithPlane(P3.XY).mapFilter(tOther => {
				const pLC = bezierLC.at(tOther)
				if (this.constructor.XYLCValid(pLC)) {
					return {tOther: tOther, p: bezierWC.at(tOther), tThis: this.constructor.XYLCPointT(pLC)}
				}
			})
			return infos
		}
	}

	isInfosWithBezier2D(bezierWC: BezierCurve, sMin?: number, sMax?: number): ISInfo[] {
		sMin = isFinite(sMin) ? sMin : bezierWC.tMin
		sMax = isFinite(sMax) ? sMax : bezierWC.tMax
		assertf(() => 0 < Math.PI)
		assertf(() => sMin < sMax)
		return Curve.ispsRecursive(this, this.tMin, this.tMax, bezierWC, sMin, sMax)
	}

	isOrthogonal(): boolean {
		return this.f1.isPerpendicularTo(this.f2)
	}

	at2(xi: number, eta: number): V3 {
		assertNumbers(xi, eta)
		// center + f1 xi + f2 eta
		return this.center.plus(this.f1.times(xi)).plus(this.f2.times(eta))
	}

	debugToMesh(mesh: Mesh, bufferName: string) {
		mesh[bufferName] || mesh.addVertexBuffer(bufferName, bufferName)
		for (let t = 0; t < Math.PI; t += 0.1) {
			const p = this.at(t)
			mesh[bufferName].push(p, p.plus(this.tangentAt(t).toLength(1)))
			mesh[bufferName].push(p, p.plus(this.normalP(t).toLength(1)))
		}
		mesh[bufferName].push(this.center, this.center.plus(this.f1.times(1.2)))
		mesh[bufferName].push(this.center, this.center.plus(this.f2))
		mesh[bufferName].push(this.center, this.center.plus(this.normal))
	}

}