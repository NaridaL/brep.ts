import { assert, assertNumbers, between, isCCW, V3 } from 'ts3dutils'
import { Mesh } from 'tsgl'

import { breakDownPPCurves, Curve, Edge, ImplicitSurface, MathFunctionR2R, PICurve, Surface } from '../index'

import { ceil, min } from '../math'

export abstract class ParametricSurface extends Surface {
	uStep!: number
	vStep!: number

	constructor(readonly uMin: number, readonly uMax: number, readonly vMin: number, readonly vMax: number) {
		super()
		assertNumbers(uMin, uMax, vMin, vMax)
		assert(uMin < uMax)
		assert(vMin < vMax)
		assert(
			(x => x[x.length - 4])(this.getConstructorParameters()) == this.uMin,
			this.getConstructorParameters(),
			this.uMin,
		)
	}

	static isCurvesParametricImplicitSurface(
		ps: ParametricSurface,
		is: ImplicitSurface,
		uStep: number,
		vStep: number = uStep,
		curveStepSize: number,
	): Curve[] {
		const pf = ps.pUVFunc(),
			icc = is.implicitFunction()
		const dpdu = ps.dpdu()
		const dpdv = ps.dpdv()
		const didp = is.didp.bind(is)
		const ist = (x: number, y: number) => icc(pf(x, y))
		const didu = (u: number, v: number) => didp(pf(u, v)).dot(dpdu(u, v))
		const didv = (u: number, v: number) => didp(pf(u, v)).dot(dpdv(u, v))
		const mf = MathFunctionR2R.forFFxFy(ist, didu, didv)
		const curves = Curve.breakDownIC(mf, ps, uStep, vStep, curveStepSize, (u, v) => is.containsPoint(pf(u, v))).map(
			({ points, tangents }, i) => PICurve.forParametricPointsTangents(ps, is, points, tangents, curveStepSize),
		)
		return curves
	}

	static isCurvesParametricParametricSurface(
		ps1: ParametricSurface,
		ps2: ParametricSurface,
		s1Step: number,
		t1Step: number = s1Step,
		curveStepSize: number,
	): Curve[] {
		return breakDownPPCurves(ps1, ps2, s1Step, t1Step, curveStepSize)
	}

	static is(obj: any): obj is ParametricSurface {
		return obj.pUVFunc
	}

	pUV(u: number, v: number): V3 {
		return this.pUVFunc()(u, v)
	}

	pUVFunc(): (u: number, v: number) => V3 {
		return this.pUV.bind(this)
	}

	abstract dpdu(): (u: number, v: number) => V3

	abstract dpdv(): (u: number, v: number) => V3

	uvP(pWC: V3): V3 {
		return this.uvPFunc()(pWC)
	}

	uvPFunc(): (pWC: V3) => V3 {
		return this.uvP.bind(this)
	}

	bounds(u: number, v: number): boolean {
		return this.uMin <= u && u <= this.uMax && this.vMin <= v && v <= this.vMax
	}

	/**
	 * Positive values are inside bounds.
	 */
	boundsSigned(u: number, v: number): number {
		return min(u - this.uMin, this.uMax - u, v - this.vMin, this.vMax - v)
	}

	normalP(p: V3): V3 {
		const pmPoint = this.uvPFunc()(p)
		return this.normalUV(pmPoint.x, pmPoint.y)
	}

	normalUVFunc(): (u: number, v: number) => V3 {
		return this.normalUV.bind(this)
	}

	normalUV(u: number, v: number): V3 {
		return this.normalUVFunc()(u, v)
	}

	parametersValid(u: number, v: number): boolean {
		return between(u, this.uMin, this.uMax) && between(v, this.vMin, this.vMax)
	}

	abstract pointFoot(pWC: V3, startU?: number, startV?: number): V3

	toMesh() {
		assert(isFinite(this.vMin) && isFinite(this.vMax) && isFinite(this.uMin) && isFinite(this.uMax))
		assert(isFinite(this.uStep) && isFinite(this.vStep))
		return Mesh.parametric(
			this.pUVFunc(),
			this.normalUVFunc(),
			this.uMin,
			this.uMax,
			this.vMin,
			this.vMax,
			ceil((this.uMax - this.uMin) / this.uStep),
			ceil((this.vMax - this.vMin) / this.vStep),
		)
	}

	isCurvesWithImplicitSurface(is: ImplicitSurface, uStep: number, vStep: number, stepSize: number): Curve[] {
		return ParametricSurface.isCurvesParametricImplicitSurface(this, is, uStep, vStep, stepSize)
	}

	edgeLoopCCW(contour: Edge[]): boolean {
		const ptpF = this.uvPFunc()
		return isCCW(contour.flatMap(e => e.getVerticesNo0()).map(v => ptpF(v)), V3.Z)
	}

	like(object: any): boolean {
		if (!this.isCoplanarTo(object)) return false
		// normals need to point in the same direction (outwards or inwards) for both
		const pSMinTMin = this.pUVFunc()(this.uMin, this.vMin)
		const thisNormal = this.normalUVFunc()(this.uMin, this.vMin)
		const otherNormal = object.normalP(pSMinTMin)
		return 0 < thisNormal.dot(otherNormal)
	}
}
