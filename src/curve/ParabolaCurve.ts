import {int} from 'ts3dutils'
import {XiEtaCurve} from './XiEtaCurve'
import {BezierCurve} from './BezierCurve'
import {L3} from './Line3'

/**
 * eta = xi²
 */
class ParabolaCurve extends XiEtaCurve {
	constructor(center: V3, f1: V3, f2: V3, tMin: number = -10, tMax: number = 10) {
		super(center, f1, f2, tMin, tMax)
	}

	at(t: number): V3 {
	    // center + f1 t + f2 t²
		return this.center.plus(this.f1.times(t)).plus(this.f2.times(t * t))
	}

	tangentAt(t: number): V3 {
		assertNumbers(t)
        // f1 + f2 2 t
		return this.f1.plus(this.f2.times(2 * t))
	}

	ddt(t: number): V3 {
		assertNumbers(t)
		return this.f2.times(2)
	}

	tangentAt2(xi: number, eta: number): V3 {
		assertNumbers(xi, eta)
		return this.f1.plus(this.f2.times(2 * eta))
	}

	reversed() {
		return new this.constructor(this.center, this.f1.negated(), this.f2, -this.tMax, -this.tMin)
	}

    /**
     * tangent: f1 + 2 * t * f2 = 0
     * t = -f1 / 2 / f2 (for individual dimensions)
     */
    roots(): [number[], number[], number[]] {
	    const dimRoots = (dim: int) => eq0(this.f2.e(dim)) ? [] : [-this.f1.e(dim) / 2 / this.f2.e(dim)]
	    return arrayFromFunction(3, dimRoots) as [number[], number[], number[]]
    }

	isColinearTo(curve: Curve): boolean {
		if (!hasConstructor(curve, ParabolaCurve)) return false
		const thisRA = this.rightAngled(), curveRA = curve.rightAngled()
		return thisRA.center.like(curveRA.center)
			&& thisRA.f2.like(curveRA.f2)
			&& thisRA.f1.likeOrReversed(curveRA.f1)
	}

	rightAngled() {
		// looking for vertex of parabola
		// this is the point where the tangent is perpendicular to the main axis (f2)
		// tangent = f1 + f2 * 2 * t0
		// f2 DOT (f1 + f2 * 2 * t0) == 0
		// f1 DOT f2 + f2 DOT f2 * 2 * t0 == 0
		// t0 == -(f1 DOT f2) / (f2 DOT f2 * 2)
		const f1 = this.f1, f2 = this.f2
        const f1DOTf2 = f1.dot(f2)
		if (eq0(f1DOTf2) && f1.hasLength(1)) {
			return this
		}
		const t0 = -f1DOTf2 / f2.squared() / 2
		// we need to rearange tMin/tMax
		// tMin' = pointT(at(tMin)) =
		const raCenter = this.at(t0)
		const raF1 = this.tangentAt(t0), raF1Length = raF1.length(), raF11 = raF1.unit()
		const repos = (t: number) => this.at(t).minus(raCenter).dot(raF11)
		return new ParabolaCurve(raCenter, raF11, f2.div(raF1Length ** 2), repos(this.tMin), repos(this.tMax))
	}

	arcLength(startT: number, endT: number): number {
		let f1 = this.f1
		const f2 = this.f2
		const f1DOTf2 = f1.dot(f2)
		let t0 = 0
		if (!eq0(f1DOTf2)) {
			t0 = -f1DOTf2 / f2.squared() / 2
			f1 = f1.plus(f2.times(2 * t0))
		}
		const f1Length = f1.length()
		const a = f2.length() / f1Length

		function F(x: number) {
			return Math.asinh(a * 2 * x) / 4 / a + x * Math.sqrt(1 + a * a * 4 * x * x) / 2
		}

		return f1Length * (F(endT - t0) - F(startT - t0))
	}

	static eccentricity() {
		return 1
	}

	static unitIsInfosWithLine(anchorLC: V3, dirLC: V3, anchorWC: V3, dirWC: V3): ISInfo[] {
		// para: x² = y
		// line(t) = anchor + t dir
		// (ax + t dx)² = ay + t dy
		// ax² + t ax dx + t² dx² = ay + t dy
		// t² dx² + t (ax dx + dy) + ay² + ay = 0
		const pqDiv = dirLC.x ** 2
		const lineTs = pqFormula((anchorLC.x * dirLC.x + dirLC.y) / pqDiv, (anchorLC.x ** 2 + anchorLC.y) / pqDiv)
		return lineTs.filter(tOther => le(0, anchorLC.y + tOther * dirLC.y))
			.map(tOther => ({
				tThis: dirLC.x * tOther + anchorLC.x,
				tOther: tOther,
				p: L3.at(anchorWC, dirWC, tOther)}))
	}

	static magic(a: number, b: number, c: number): number[] {
		/*
		 solve system (5)/(6)
		 g1 * xi + g2 * eta = g3 (6)
		 g1 * xi + g2 * xi * xi = g3
		 xi² + xi * g1/g2 - g3/g2 = 0
		 */
    	return pqFormula(a / b, -c / b)
	}

	static XYLCValid(pLC: V3): boolean {
		return eq(pLC.x ** 2, pLC.y)
	}

	static XYLCPointT(pLC: V3): number {
		return pLC.x
	}

	static quadratic(a: V3, b: V3, c: V3): ParabolaCurve {
        // (1 - t)² a + 2 * t * (1 - t) b + t² c
        // (1 -2t +t²)a + (2t -2t²) b + t² c
        // = t²(a - 2b + c) + t (-2a + 2b) + a
        // (2t - 2) a + (1 - 2t) b + 2t c = t(2a + 2b - 2c) - 2a + b
        // 2 a + -2 b + 2 c
        const f2 = a.plus(c).minus(b.times(2))
        const f1 = b.minus(a).times(2)
        const center = a
        return new ParabolaCurve(center, f1, f2, 0, 1)
    }

    asBezier() {
	    return BezierCurve.quadratic(
	        this.at(-1),
            new L3(this.at(-1), this.tangentAt(-1).unit()).isInfoWithLine(new L3(this.at(1), this.tangentAt(1).unit())),
            this.at(1))
    }

	static readonly XY = new ParabolaCurve(V3.O, V3.X, V3.Y)
	static readonly YZ = new ParabolaCurve(V3.O, V3.Y, V3.Z)
	static readonly ZX = new ParabolaCurve(V3.O, V3.Z, V3.X)
}
ParabolaCurve.prototype.tIncrement = 1/32
