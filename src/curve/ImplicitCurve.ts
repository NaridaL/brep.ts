class ImplicitCurve extends Curve {
	constructor(readonly points: V3[],
	            readonly tangents: V3[],
	            readonly dir: number = 1,
				readonly generator?: string,
				tMin: number = (1 == dir ? 0 : -(points.length - 1)),
				tMax: number = (1 == dir ? points.length - 1 : 0)) {
		super(tMin, tMax)
		assert(points.length > 2)
		assert(0 <= tMin && tMin <= points.length - 1)
		assert(0 <= tMax && tMax <= points.length - 1)
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
		return this == obj ||
			Object.getPrototypeOf(obj) == PICurve.prototype
			&& this.points[0].equals(obj.points[0])
			&& this.tangents[0].equals(obj.tangents[0])
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
		assert(!isNaN(t))
		return V3.lerp(this.points[floor(t)], this.points[ceil(t)], t % 1)
	}

	getConstructorParameters(): any[] {
		return []
	}

	transform(m4: M4): ImplicitCurve {
	    return new ImplicitCurve(
		    m4.transformedPoints(this.points),
		    m4.transformedVectors(this.tangents))
    }

    roots(): [number[], number[], number[]] {
		const allTs = arrayRange(0, this.points.length)
		return [allTs, allTs, allTs]
    }

	addToMesh(mesh: Mesh, res: int = 4, radius: number = 0, pointStep = 1): void {
		const baseNormals = arrayFromFunction(res, i => V3.polar(1, TAU * i / res))
		const baseVertices = arrayFromFunction(res, i => V3.polar(radius, TAU * i / res))
		let prevTangent = V3.Z, prevMatrix = M4.IDENTITY
		for (let i = ceil(this.tMin); i < floor(this.tMax); i += pointStep) {

			const start = mesh.vertices.length
			if (ceil(this.tMin) !== i) {
				for (let j = 0; j < res; j++) {
					pushQuad(mesh.triangles, true,
						start - res + j, start + j,
						start - res + (j + 1) % res, start + (j + 1) % res)
				}
			}
			const point = this.points[i], tangent = this.tangents[i]
			const tangentMatrix = M4.rotationAB(prevTangent, tangent).times(prevMatrix)
			mesh.normals.pushAll(tangentMatrix.transformedVectors(baseNormals))
			const baseMatrix = M4.translation(point).times(tangentMatrix)
			mesh.vertices.pushAll(baseMatrix.transformedPoints(baseVertices))
			prevTangent = tangent
			prevMatrix = tangentMatrix
		}
	}
}
ImplicitCurve.prototype.tIncrement = 1