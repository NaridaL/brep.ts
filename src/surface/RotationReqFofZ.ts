/**
 * Rotation surface with r = f(z)
 */
class RotationReqFofZ extends Surface {
	constructor(readonly axis: L3,
	            readonly FofR: (r: number) => number,
	            readonly tMin: number,
	            readonly tMax: number,
				readonly normalDir: number,
				readonly drdz: (r: number) => number = z => (FofR(z + eps) - FofR(z)) / eps) {
		super()
		assertInst(L3, axis)
		assert(1 == normalDir || -1 == normalDir)
	}


	getConstructorParameters(): any[] {
		return [this.axis, this.FofR, this.tMin, this.tMax, this.normalDir, this.drdz]
	}

	toMesh(): Mesh {
		return Mesh.parametric(this.parametricFunction(), this.parametricNormal(), this.sMin, this.sMax, this.tMin, this.tMax, 32, 40)
	}

	flipped(): RotationReqFofZ {
		return new RotationReqFofZ(this.axis, this.FofR, this.tMin, this.tMax, -this.normalDir, this.drdz)
	}


	transform(m4: M4): RotationReqFofZ {
		return new RotationReqFofZ(axis.transform(m4), )
	}

	containsPoint(p: V3): boolean {
		return eq0(this.implicitFunction()(p))
	}

	parametricFunction() {
		const z = this.axis.dir1, x = z.getPerpendicular().unit(), y = z.cross(x)
		const matrix = M4.forSys(x, y, z, this.axis.anchor)
		const f = this.FofR
		return function (d, z) {
			const radius = f(z)
			return matrix.transformPoint(new V3(radius * cos(d), radius * sin(d), z))
		}
	}

	parametricNormal() {
		const z = this.axis.dir1, x = z.getPerpendicular().unit(), y = z.cross(x)
		const matrix = M4.forSys(x, y, z, this.axis.anchor).inversed().transposed()
		return (d, z) => {
			const fz = this.FofR(z)
			const drdz = this.drdz(z)
			return matrix.transformVector(new V3(cos(d), sin(d), -drdz)).toLength(this.normalDir)
		}
	}

	implicitFunction() {
		const z = this.axis.dir1, x = z.getPerpendicular().unit(), y = z.cross(x)
		const matrix = M4.forSys(x, y, z, this.axis.anchor)
		const matrixInverse = matrix.inversed()
		const f = this.FofR
		return function (pWC) {
			const p = matrixInverse.transformPoint(pWC)
			const radiusLC = Math.sqrt(p.x * p.x + p.y * p.y)
			return f(p.z) - radiusLC
		}
	}

	boundsFunction() {
		const z = this.axis.dir1, x = z.getPerpendicular().unit(), y = z.cross(x)
		const matrix = M4.forSys(x, y, z, this.axis.anchor)
		const matrixInverse = matrix.inversed()
		const f = this.FofR, minZ = this.minZ, maxZ = this.maxZ
		return function (pWC) {
			const z = matrixInverse.transformPoint(pWC).z
			return minZ <= z && z <= maxZ
		}
	}

	pointToParameterFunction(p) {
		const z = this.axis.dir1, x = z.getPerpendicular().unit(), y = z.cross(x)
		const matrix = M4.forSys(x, y, z, this.axis.anchor)
		const matrixInverse = matrix.inversed()
		const f = this.FofR
		return function (pWC) {
			const p = matrixInverse.transformPoint(pWC)
			return new V3(atan2(p.y, p.x), p.z, 0)
		}
	}
}
RotationReqFofZ.prototype.sMin = 0
RotationReqFofZ.prototype.sMax = PI