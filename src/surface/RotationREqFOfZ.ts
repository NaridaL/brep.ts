/**
 * Rotation surface with r = f(z)
 */
class RotationREqFOfZ extends ParametricSurface implements ImplicitSurface {
	matrixInverse: M4
	constructor(readonly matrix: M4,
	            readonly rt: (z: number) => number, // r(z)
	            readonly tMin: number,
	            readonly tMax: number,
				readonly normalDir: number,
				readonly drdz: (z: number) => number = z => (rt(z + eps) - rt(z)) / eps) { // d/dz (r(z))
		super()
		assertInst(M4, matrix)
		assert(matrix.isNoProj())
		assert(1 == normalDir || -1 == normalDir)
		this.matrixInverse = matrix.inversed()
	}

	getConstructorParameters(): any[] {
		return [this.matrix, this.rt, this.tMin, this.tMax, this.normalDir, this.drdz]
	}

	flipped(): RotationREqFOfZ {
		return new RotationREqFOfZ(this.matrix, this.rt, this.tMin, this.tMax, -this.normalDir, this.drdz)
	}

	transform(m4: M4): RotationREqFOfZ {
		return new RotationREqFOfZ(m4.times(this.matrix), this.rt, this.tMin, this.tMax, this.normalDir, this.drdz)
	}

	containsPoint(p: V3): boolean {
		return eq0(this.implicitFunction()(p))
	}

	pSTFunc(): (s: number, t: number) => V3 {
		return (d, z) => {
			const radius = this.rt(z)
			return this.matrix.transformPoint(V3.polar(radius, d, z))
		}
	}

	dpds(): (s: number, t: number) => V3 {
		return (s, t) => {
			const radius = this.rt(t)
			return this.matrix.transformVector(new V3(radius * -sin(s), radius * cos(s), 0))
		}
	}

	/**
	 * new V3(f(z) * cos d, f(z) * sin d, z)
	 */

	dpdt(): (s: number, t: number) => V3 {
		return (s, t) => {
			const drdt = this.drdz(t)
			return this.matrix.transformVector(new V3(drdt * cos(s), drdt * sin(s), 1))
		}
	}

	normalSTFunc(): (s: number, t: number) => V3 {
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
		return (pWC) => {
			const pLC = this.matrixInverse.transformPoint(pWC)
			const radiusLC = pLC.lengthXY()
			return this.rt(pLC.z) - radiusLC
		}
	}

	stPFunc(): (pWC: V3) => V3 {
		return (pWC) => {
			const pLC = this.matrixInverse.transformPoint(pWC)
			return new V3(atan2(pLC.y, pLC.x), pLC.z, 0)
		}
	}
}
Object.assign(RotationREqFOfZ.prototype, ImplicitSurface.prototype)
RotationREqFOfZ.prototype.sMin = 0
RotationREqFOfZ.prototype.sMax = PI