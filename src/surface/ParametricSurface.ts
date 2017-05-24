///<reference path="Surface.ts"/>

abstract class ParametricSurface extends Surface {
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

	footParameters(pWC: V3, ss: number, st: number): V3 {
		throw new Error()
	}

	toMesh(): Mesh {
		assert(isFinite(this.tMin) && isFinite(this.tMax)&&isFinite(this.sMin)&&isFinite(this.sMax))
		return Mesh.parametric(this.pSTFunc(), this.normalSTFunc(),
			this.sMin, this.sMax, this.tMin, this.tMax,
			ceil((this.sMax - this.sMin) / this.uStep),
			ceil((this.tMax - this.tMin) / this.vStep))
	}

	sMin: number
	sMax: number
	tMin: number
	tMax: number
	uStep: number
	vStep: number

	isCurvesWithImplicitSurface(is: ImplicitSurface, sStep: number, tStep: number, stepSize: number): Curve[] {
		return ParametricSurface.isCurvesParametricImplicitSurface(this, is, sStep, tStep, stepSize)
	}
	static isCurvesParametricImplicitSurface(ps: ParametricSurface,
	                                         is: ImplicitSurface,
	                                         sStep: number,
	                                         tStep: number = sStep,
	                                         curveStepSize: number): Curve[] {
		const pf = ps.pSTFunc(), icc = is.implicitFunction()
		const dpds = ps.dpds()
		const dpdt = ps.dpdt()
		const didp = is.didp.bind(is)
		const ist = (x: number, y: number) => icc(pf(x, y))
		const dids = (s: number, t: number) => didp(pf(s, t)).dot(dpds(s, t))
		const didt = (s: number, t: number) => didp(pf(s, t)).dot(dpdt(s, t))
		const mf = MathFunctionR2R.forFFxFy(ist, dids, didt)
		const curves
			= Curve.breakDownIC(mf, ps, sStep, tStep, curveStepSize, dids, didt)
				.map(({points, tangents}, i) => PICurve.forParametricPointsTangents(ps, is, points, tangents, curveStepSize))
		return curves
	}

	static is(surface: any): surface is ParametricSurface {
		return surface.pSTFunc
	}
}
abstract class ImplicitSurface extends Surface {
	abstract implicitFunction(): (pWC: V3) => number

	abstract didp(pWC: V3): V3

	static is(obj: any): obj is ImplicitSurface {
		return obj.implicitFunction
	}
}