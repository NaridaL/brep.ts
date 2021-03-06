import {
	arrayFromFunction,
	arrayRange,
	assert,
	assertVectors,
	bisect,
	clamp,
	eq,
	eq0,
	int,
	M4,
	TAU,
	Tuple3,
	V3,
} from 'ts3dutils'
import { MeshWith, pushQuad } from 'tsgl'

import { Curve, L3, PICurve, Surface } from '../index'
import { ceil, floor, max, min } from '../math'

export abstract class ImplicitCurve extends Curve {
	constructor(
		readonly points: ReadonlyArray<V3>,
		readonly tangents: ReadonlyArray<V3>,
		readonly dir: number = 1,
		readonly generator?: string,
		tMin: number = 1 == dir ? 0 : -(points.length - 1),
		tMax: number = 1 == dir ? points.length - 1 : 0,
	) {
		super(tMin, tMax)
		assert(points.length > 2)
		assert(0 <= tMin && tMin <= points.length - 1, tMin, points.length)
		assert(0 <= tMax && tMax <= points.length - 1, tMax, points.length)
	}

	likeCurve(curve: Curve): boolean {
		throw new Error('Method not implemented.')
	}

	toSource(rounder: (x: number) => number = x => x): string {
		return this.generator || super.toSource(rounder)
	}

	containsPoint(p: V3): boolean {
		assertVectors(p)
		return !isNaN(this.pointT(p))
	}

	equals(obj: any): boolean {
		return (
			this == obj ||
			(Object.getPrototypeOf(obj) == PICurve.prototype &&
				this.points[0].equals(obj.points[0]) &&
				this.tangents[0].equals(obj.tangents[0]))
		)
	}

	hashCode(): int {
		return [this.points[0], this.tangents[0]].hashCode()
	}

	tangentP(pWC: V3): V3 {
		assertVectors(pWC)
		assert(this.containsPoint(pWC), 'this.containsPoint(pWC)' + this.containsPoint(pWC))
		const t = this.pointT(pWC)
		return this.tangentAt(t)
	}

	tangentAt(t: number): V3 {
		t = clamp(t, this.tMin, this.tMax)
		return V3.lerp(this.tangents[floor(t)], this.tangents[ceil(t)], t % 1)
	}

	at(t: number): V3 {
		assert(isFinite(t))
		return V3.lerp(this.points[floor(t)], this.points[ceil(t)], t % 1)
	}

	getConstructorParameters(): any[] {
		throw new Error()
	}

	roots(): [number[], number[], number[]] {
		const allTs = arrayRange(0, this.points.length)
		return [allTs, allTs, allTs]
	}

	/**
	 * @param mesh
	 * @param res
	 * @param radius default to 0. Use the shader to achieve dynamic scaling.
	 * @param pointStep
	 */
	addToMesh(mesh: MeshWith<'normals' | 'TRIANGLES'>, res: int = 4, radius: number = 0, pointStep = 1): void {
		const baseNormals = arrayFromFunction(res, i => V3.polar(1, TAU * i / res))
		const baseVertices = arrayFromFunction(res, i => V3.polar(radius, TAU * i / res))
		let prevTangent = V3.Z,
			prevMatrix = M4.IDENTITY
		for (let i = 0; i < this.points.length; i += pointStep) {
			const start = mesh.vertices.length
			if (0 !== i) {
				for (let j = 0; j < res; j++) {
					pushQuad(
						mesh.TRIANGLES,
						true,
						start - res + j,
						start + j,
						start - res + (j + 1) % res,
						start + (j + 1) % res,
					)
				}
			}
			const point = this.points[i],
				tangent = this.tangents[i]
			const tangentMatrix = M4.rotateAB(prevTangent, tangent).times(prevMatrix)
			mesh.normals.push(...tangentMatrix.transformedVectors(baseNormals))
			const baseMatrix = M4.translate(point).times(tangentMatrix)
			mesh.vertices.push(...baseMatrix.transformedPoints(baseVertices))
			prevTangent = tangent
			prevMatrix = tangentMatrix
		}
	}

	rootsApprox() {
		const roots: Tuple3<number[]> = [[], [], []]
		const points = this.points
		let lastDiff = points[1].minus(points[0])
		for (let i = 2; i < points.length; i++) {
			const diff = points[i].minus(points[i - 1])
			for (let dim = 0; dim < 3; dim++) {
				if (Math.sign(lastDiff.e(dim)) != Math.sign(diff.e(dim))) {
					roots[dim].push(i)
				}
			}
			lastDiff = diff
		}
		return roots
	}

	pointT(pWC: V3): number {
		const startT = arrayRange(floor(this.tMin), ceil(this.tMax), 1).withMax(t => -pWC.distanceTo(this.points[t]))
		if (undefined === startT) throw new Error()
		if (this.points[startT].like(pWC)) return startT
		const a = max(0, startT - 1),
			b = min(this.points.length - 1, startT + 1)
		const tangent = this.tangentAt(startT)
		const f = (t: number) =>
			this.at(t)
				.to(pWC)
				.dot(tangent)
		// const df = (t: number) => -this.tangentAt(clamp(t, 0, this.points.length - 1)).dot(tangent)
		//checkDerivate(f, df, 0, this.points.length - 2, 3)
		const t = bisect(f, a, b, 32)
		if (!isFinite(t) || !eq0(this.at(t).distanceTo(pWC))) {
			return NaN
		}
		return t
	}
}

ImplicitCurve.prototype.tIncrement = 1

/**
 * isInfosWithLine for an ImplicitCurve defined as the intersection of two surfaces.
 */
export function surfaceIsICurveIsInfosWithLine(
	this: ImplicitCurve,
	surface1: Surface,
	surface2: Surface,
	anchorWC: V3,
	dirWC: V3,
	tMin?: number | undefined,
	tMax?: number | undefined,
	lineMin?: number | undefined,
	lineMax?: number | undefined,
) {
	const line = new L3(anchorWC, dirWC.unit())
	const psTs = surface1.isTsForLine(line)
	const isTs = surface2.isTsForLine(line)
	const commonTs = psTs.filter(psT => isTs.some(isT => eq(psT, isT)))
	const commonTInfos = commonTs.map(t => ({ tThis: 0, tOther: t / dirWC.length(), p: line.at(t) }))
	const result = commonTInfos.filter(info => this.containsPoint(info.p))
	result.forEach(info => (info.tThis = this.pointT(info.p)))
}
