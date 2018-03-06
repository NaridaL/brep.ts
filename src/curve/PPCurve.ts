import {ParametricSurface, ImplicitCurve, curvePointPP, Curve, followAlgorithmPP} from '../index'
import {assert, V3, Tuple3, assertVectors, newtonIterate, Tuple4, assertNever, M4, callsce} from 'ts3dutils'

import {floor, ceil} from '../math'

export class PPCurve extends ImplicitCurve {
    constructor(points: ReadonlyArray<V3>,
                tangents: ReadonlyArray<V3>,
                readonly parametricSurface1: ParametricSurface,
                readonly parametricSurface2: ParametricSurface,
                readonly st1s: ReadonlyArray<V3>,
                readonly pmTangents: ReadonlyArray<V3>,
                readonly stepSize: number,
                dir: number = 1,
                generator?: string,
                tMin?: number, tMax?: number) {
        super(points, tangents, dir, generator, tMin, tMax)
		assert(ParametricSurface.is(parametricSurface1))
		assert(ParametricSurface.is(parametricSurface2))
        assert(Array.isArray(st1s))
        assert(dir == 1)
        assert(stepSize <= 1)
	}

	at(t: number) {
        assert(!isNaN(t))
        if (0 === t % 1) return this.points[t]
        const startPoint = V3.lerp(this.points[floor(t)], this.points[ceil(t)], t % 1)
        return curvePointPP(this.parametricSurface1, this.parametricSurface2, startPoint)!.p
    }

    isColinearTo(curve: Curve) {
        if (curve instanceof PPCurve) {
            if (this.equals(curve)) {
                return true
            }
            if (this.parametricSurface1.isCoplanarTo(curve.parametricSurface1) && this.parametricSurface1.isCoplanarTo(curve.parametricSurface2)) {
                // TODO
            }
            return false
        } else {
            return false
        }
    }

	containsPoint(p: V3) {
		assertVectors(p)
		// TODO: wrong, as there could be another curve
		return this.parametricSurface1.containsPoint(p) && this.parametricSurface2.containsPoint(p) && !isNaN(this.pointT(p))
	}

	rootsAprox() {
		const roots: Tuple3<number[]> = [[], [], []]
		const ps = this.points
		let lastDiff = ps[1].minus(ps[0])
		for (let i = 2; i < ps.length; i++) {
			const diff = ps[i].minus(ps[i - 1])
			for (let dim = 0; dim < 3; dim++) {
				if (Math.sign(lastDiff.e(dim)) != Math.sign(diff.e(dim))) {
					roots[dim].push(i)
				}
			}
			lastDiff = diff
		}
		return roots
	}

	rootPoints() {
		const pF1 = this.parametricSurface1.pSTFunc()
		const pF2 = this.parametricSurface2.pSTFunc()
		const pN1 = this.parametricSurface1.normalSTFunc()
		const pN2 = this.parametricSurface2.normalSTFunc()

		const rootsAprox = this.rootsAprox()
		const results: Tuple3<V3[]> = [[], [], []]
		for (let dim = 0; dim < 3; dim++) {
			for (let i = 0; i < rootsAprox[dim].length; i++) {
				const lambda = rootsAprox[dim][i]
				const p = this.at(lambda)
				assert(this.parametricSurface1.containsPoint(p))
				const pp1 = this.parametricSurface1.stP(p)
				const {x: u, y: v} = this.parametricSurface2.stP(p)
				const startValues = [pp1.x, pp1.y, u, v]

				function f(vals: Tuple4<number>) {
					const [s, t, u, v] = vals
					const diff = pF1(s, t).minus(pF2(u, v))
					const n1 = pN1(s, t)
					const n2 = pN2(u, v)
					const tangent = n1.cross(n2)
					return [diff.x, diff.y, diff.z, tangent.e(dim)]
				}

				const pps = newtonIterate(f, startValues, 8)
				// assert(pF1(pps[0], pps[1]).like(pF2(pps[2], pps[3])),
				// 	pF1(pps[0], pps[1]).sce + pF2(pps[2], pps[3]).sce)
				const result = pF1(pps[0], pps[1])
				results[dim].push(result)
			}
		}
		return results
	}

	roots(): Tuple3<number[]> {
        return this.rootPoints().map(ps => ps.map(p => this.pointT(p)))
    }

	pointTangent(pWC: V3): V3 {
		assertVectors(pWC)
		assert(this.containsPoint(pWC), 'this.containsPoint(pWC)')
		const n1 = this.parametricSurface1.normalP(pWC)
		const n2 = this.parametricSurface2.normalP(pWC)
		return n1.cross(n2)
	}

	pointT(point: V3) {
	}

    transform(m4: M4, desc?: string) {
        const dirFactor = m4.isMirroring() ? -1 : 1
        return new PPCurve(
            m4.transformedPoints(this.points),
            m4.transformedVectors(this.tangents),
            this.parametricSurface1.transform(m4),
            this.parametricSurface2.transform(m4),
            this.st1s,
            undefined,
            this.stepSize,
            this.dir,
            undefined) as this
    }

    toSource(rounder: (x: number) => number = x => x): string {
        const result = callsce('PPCurve.forStartEnd',
            this.parametricSurface1, this.parametricSurface2,
            this.points[0], this.points.last,
            this.stepSize)
        return result
    }

    static forStartEnd(ps1: ParametricSurface, ps2: ParametricSurface, startPoint: V3, end: V3, stepSize = 0.02) {
        const {points, tangents, st1s} = followAlgorithmPP(ps1, ps2, startPoint, stepSize)
        return new PPCurve(points, tangents, ps1, ps2, st1s, undefined, stepSize, 1)
    }

}
