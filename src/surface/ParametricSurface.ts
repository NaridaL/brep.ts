import { assert, assertNumbers, between, isCCW, V3 } from 'ts3dutils'
import { Mesh } from 'tsgl'

import { breakDownPPCurves, Curve, Edge, ImplicitSurface, MathFunctionR2R, PICurve, Surface } from '../index'

import { ceil, min } from '../math'

export abstract class ParametricSurface extends Surface {
	uStep: number
	vStep: number

	constructor(readonly sMin: number, readonly sMax: number, readonly tMin: number, readonly tMax: number) {
		super()
		assertNumbers(sMin, sMax, tMin, tMax)
		assert(sMin < sMax)
		assert(tMin < tMax)
	}

	static isCurvesParametricImplicitSurface(
		ps: ParametricSurface,
		is: ImplicitSurface,
		sStep: number,
		tStep: number = sStep,
		curveStepSize: number,
	): Curve[] {
		const pf = ps.pSTFunc(),
			icc = is.implicitFunction()
		const dpds = ps.dpds()
		const dpdt = ps.dpdt()
		const didp = is.didp.bind(is)
		const ist = (x: number, y: number) => icc(pf(x, y))
		const dids = (s: number, t: number) => didp(pf(s, t)).dot(dpds(s, t))
		const didt = (s: number, t: number) => didp(pf(s, t)).dot(dpdt(s, t))
		const mf = MathFunctionR2R.forFFxFy(ist, dids, didt)
		const curves = Curve.breakDownIC(mf, ps, sStep, tStep, curveStepSize).map(({ points, tangents }, i) =>
			PICurve.forParametricPointsTangents(ps, is, points, tangents, curveStepSize),
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
		return obj.pSTFunc
	}

	pST(s: number, t: number): V3 {
		return this.pSTFunc()(s, t)
	}

	pSTFunc(): (s: number, t: number) => V3 {
		return this.pST.bind(this)
	}

	abstract dpds(): (s: number, t: number) => V3

	abstract dpdt(): (s: number, t: number) => V3

	stP(pWC: V3): V3 {
		return this.stPFunc()(pWC)
	}

	stPFunc(): (pWC: V3) => V3 {
		return this.stP.bind(this)
	}

	bounds(s: number, t: number): boolean {
		return this.sMin <= s && s <= this.sMax && this.tMin <= t && t <= this.tMax
	}

	/**
	 * Positive values are inside bounds.
	 */
	boundsSigned(s: number, t: number): number {
		return min(s - this.sMin, this.sMax - s, t - this.tMin, this.tMax - t)
	}

	normalP(p: V3): V3 {
		const pmPoint = this.stPFunc()(p)
		return this.normalST(pmPoint.x, pmPoint.y)
	}

	normalSTFunc(): (s: number, t: number) => V3 {
		return this.normalST.bind(this)
	}

	normalST(s: number, t: number): V3 {
		return this.normalSTFunc()(s, t)
	}

	parametersValid(s: number, t: number): boolean {
		return between(s, this.sMin, this.sMax) && between(t, this.tMin, this.tMax)
	}

	abstract pointFoot(pWC: V3, ss?: number, st?: number): V3

	toMesh() {
		assert(isFinite(this.tMin) && isFinite(this.tMax) && isFinite(this.sMin) && isFinite(this.sMax))
		assert(isFinite(this.uStep) && isFinite(this.vStep))
		return Mesh.parametric(
			this.pSTFunc(),
			this.normalSTFunc(),
			this.sMin,
			this.sMax,
			this.tMin,
			this.tMax,
			ceil((this.sMax - this.sMin) / this.uStep),
			ceil((this.tMax - this.tMin) / this.vStep),
		)
	}

	isCurvesWithImplicitSurface(is: ImplicitSurface, sStep: number, tStep: number, stepSize: number): Curve[] {
		return ParametricSurface.isCurvesParametricImplicitSurface(this, is, sStep, tStep, stepSize)
	}

	edgeLoopCCW(contour: Edge[]): boolean {
		const ptpF = this.stPFunc()
		return isCCW(contour.flatMap(e => e.getVerticesNo0()).map(v => ptpF(v)), V3.Z)
	}
}
