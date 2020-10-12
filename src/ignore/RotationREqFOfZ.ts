import { arrayFromFunction, assert, assertInst, eq0, M4, newtonIterate1d, V, V3 } from '../../../ts3dutils/index'

import { EPS, ImplicitSurface, ParametricSurface } from '../index'

import { atan2, cos, PI, sin } from '../math'

/**
 * Rotation surface with r = f(z)
 */
export class RotationREqFOfZ extends ParametricSurface implements ImplicitSurface {
	readonly matrixInverse: M4

	constructor(
		readonly matrix: M4,
		readonly rt: (z: number) => number, // r(z)
		readonly vMin: number,
		readonly vMax: number,
		readonly normalDir: number,
		readonly drdz: (z: number) => number = z => (rt(z + EPS) - rt(z)) / EPS,
	) {
		// d/dz (r(z))
		super(0, PI, vMin, vMax)
		assertInst(M4, matrix)
		assert(matrix.isNoProj())
		assert(1 == normalDir || -1 == normalDir)
		this.matrixInverse = matrix.inversed()
	}

	getConstructorParameters(): any[] {
		return [this.matrix, this.rt, this.vMin, this.vMax, this.normalDir, this.drdz]
	}

	flipped(): this {
		return new RotationREqFOfZ(this.matrix, this.rt, this.vMin, this.vMax, -this.normalDir, this.drdz) as this
	}

	transform(m4: M4): this {
		return new RotationREqFOfZ(
			m4.times(this.matrix),
			this.rt,
			this.vMin,
			this.vMax,
			this.normalDir,
			this.drdz,
		) as this
	}

	containsPoint(p: V3): boolean {
		return eq0(this.implicitFunction()(p))
	}

	pUVFunc(): (u: number, v: number) => V3 {
		return (d, z) => {
			const radius = this.rt(z)
			return this.matrix.transformPoint(V3.polar(radius, d, z))
		}
	}

	dpdu(): (u: number, v: number) => V3 {
		return (u, v) => {
			const radius = this.rt(t)
			return this.matrix.transformVector(new V3(radius * -sin(s), radius * cos(s), 0))
		}
	}

	/**
	 * new V3(f(z) * cos d, f(z) * sin d, z)
	 */

	dpdv(): (u: number, v: number) => V3 {
		return (u, v) => {
			const drdt = this.drdz(t)
			return this.matrix.transformVector(new V3(drdt * cos(s), drdt * sin(s), 1))
		}
	}

	normalUVFunc(): (u: number, v: number) => V3 {
		/**
		 * (radius * -sin(s), radius * cos(s), 0) X (drds * cos(s), drds * sin(s), 1)
		 * =(radius * cos(s)*1,
		 * -radius * -sin(s)*1,
		 * radius * -sin(s)* drds * sin(s)- radius * cos(s)*drds * cos(s))
		 * div by radius
		 * => (cos s, sin s, -drds * (sin² + cos²))
		 */
		const matrix = this.matrix.inversed().transposed()
		return (d, z) => {
			const drdz = this.drdz(z)
			return matrix.transformVector(V3.polar(1, d, -drdz)).toLength(this.normalDir)
		}
	}

	implicitFunction(): (pWC: V3) => number {
		return pWC => {
			const pLC = this.matrixInverse.transformPoint(pWC)
			const radiusLC = pLC.lengthXY()
			return this.rt(pLC.z) - radiusLC
		}
	}

	uvPFunc(): (pWC: V3) => V3 {
		return pWC => {
			const pLC = this.matrixInverse.transformPoint(pWC)
			return new V3(atan2(pLC.y, pLC.x), pLC.z, 0)
		}
	}

	pointFoot(pWC: V3, startU?: number, startV?: number): V3 {
		const pLC = this.matrixInverse.transformPoint(pWC)
		return new V3(
			pLC.angleXY(),
			closestXToP(this.rt, this.drdz, this.vMin, this.vMax, new V3(pLC.z, pLC.lengthXY(), 0), startV),
			0,
		)
	}
}

RotationREqFOfZ.prototype.uStep = PI / 16
RotationREqFOfZ.prototype.vStep = 1 / 4

Object.assign(RotationREqFOfZ.prototype, ImplicitSurface.prototype)

/**
 * For a function f(x) = y, returns the parameter x for which (x, f(x)) is closest to a point (px, py)
 * @param f
 * @param df
 * @param xMin
 * @param xMax
 * @param p
 * @param startX
 */
function closestXToP(
	f: (x: number) => number,
	df: (x: number) => number,
	xMin: number,
	xMax: number,
	p: V3,
	startX?: number,
) {
	const STEPS = 32
	if (undefined === startX) {
		startX = arrayFromFunction(STEPS, i => xMin + ((xMax - xMin) * i) / STEPS).withMax(
			x => -Math.hypot(x - p.x, f(x) - p.y),
		)
	}

	return newtonIterate1d(
		(t: number) =>
			V(t, f(t))
				.minus(p)
				.dot(V(t, df(t))),
		startX,
		4,
	)
}
