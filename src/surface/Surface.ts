import { Equalable } from 'javasetmap.ts'
import { callsce, eq, eq0, int, le, NLA_PRECISION, Transformable, V3 } from 'ts3dutils'

import {
	CalculateAreaVisitor,
	Curve,
	dotCurve2,
	Edge,
	EllipseCurve,
	ImplicitCurve,
	L3,
	P3,
	PICurve,
	PPCurve,
	ZDirVolumeVisitor,
} from '../index'

import { ceil, floor, PI, sign } from '../math'

export abstract class Surface extends Transformable implements Equalable {
	readonly ['constructor']: new (...args: any[]) => this
	static loopContainsPointGeneral(loop: Edge[], pWC: V3, testLine: L3, lineOut: V3): PointVsFace {
		const testPlane = P3.normalOnAnchor(lineOut, pWC)
		// edges colinear to the testing line; these will always be counted as "inside" relative to the testing line
		const colinearEdges = loop.map(edge => edge.colinearToLine(testLine))
		let inside = false

		function logIS(isP: V3) {
			const isT = testLine.pointT(isP)
			if (eq0(isT)) {
				return true
			} else if (isT > 0) {
				inside = !inside
			}
			return false
		}

		for (let edgeIndex = 0; edgeIndex < loop.length; edgeIndex++) {
			const edge = loop[edgeIndex]
			const nextEdgeIndex = (edgeIndex + 1) % loop.length,
				nextEdge = loop[nextEdgeIndex]
			//console.log(edge.toSource()) {p:V(2, -2.102, 0),
			if (colinearEdges[edgeIndex]) {
				const lineAT = testLine.pointT(edge.a),
					lineBT = testLine.pointT(edge.b)
				if (Math.min(lineAT, lineBT) <= NLA_PRECISION && -NLA_PRECISION <= Math.max(lineAT, lineBT)) {
					return PointVsFace.ON_EDGE
				}
				// edge colinear to intersection
				const nextInside =
					colinearEdges[nextEdgeIndex] ||
					dotCurve2(nextEdge.curve, nextEdge.aT, lineOut, sign(nextEdge.deltaT())) < 0
				if (!nextInside) {
					if (logIS(edge.b)) return PointVsFace.ON_EDGE
				}
			} else {
				for (const edgeT of edge.edgeISTsWithPlane(testPlane)) {
					if (edgeT == edge.bT) {
						if (!testLine.containsPoint(edge.b)) continue
						// endpoint lies on intersection line
						if (edge.b.like(pWC)) {
							// TODO: refactor, dont check for different sides, just logIs everything
							return PointVsFace.ON_EDGE
						}
						const edgeInside = dotCurve2(edge.curve, edge.bT, lineOut, -sign(edge.deltaT())) < 0
						const nextInside =
							colinearEdges[nextEdgeIndex] ||
							dotCurve2(nextEdge.curve, nextEdge.aT, lineOut, sign(nextEdge.deltaT())) < 0
						if (edgeInside != nextInside) {
							if (logIS(edge.b)) return PointVsFace.ON_EDGE
						}
					} else if (edgeT != edge.aT) {
						const p = edge.curve.at(edgeT)
						if (!testLine.containsPoint(p)) continue
						// edge crosses line, neither starts nor ends on it
						if (logIS(p)) return PointVsFace.ON_EDGE
						// TODO: tangents?
					}
				}
			}
		}
		return inside ? PointVsFace.INSIDE : PointVsFace.OUTSIDE
	}

	static loopContainsPointEllipse(loop: Edge[], pWC: V3, testLine: EllipseCurve, pWCT?: number): PointVsFace {
		const lineOut = testLine.normal
		const testPlane = P3.normalOnAnchor(testLine.normal, pWC)
		const colinearEdges = loop.map(edge => testLine.isColinearTo(edge.curve))
		let inside = false
		if (undefined === pWCT) {
			pWCT = testLine.pointT(pWC)
		}
		const pT = pWCT!

		function logIS(isP: V3) {
			const isT = testLine.pointT(isP)
			if (eq(pT, isT)) {
				return true
			} else if (pT < isT && le(isT, PI)) {
				inside = !inside
			}
			return false
		}

		for (let edgeIndex = 0; edgeIndex < loop.length; edgeIndex++) {
			const edge = loop[edgeIndex]
			const nextEdgeIndex = (edgeIndex + 1) % loop.length,
				nextEdge = loop[nextEdgeIndex]
			//console.log(edge.toSource()) {p:V(2, -2.102, 0),
			if (colinearEdges[edgeIndex]) {
				let edgeT
				if (
					edge.curve.containsPoint(pWC) &&
					le(edge.minT, (edgeT = edge.curve.pointT(pWC))) &&
					le(edgeT, edge.maxT)
				) {
					return PointVsFace.ON_EDGE
				}
				// edge colinear to intersection
				const nextInside =
					colinearEdges[nextEdgeIndex] ||
					dotCurve2(nextEdge.curve, nextEdge.aT, lineOut, sign(nextEdge.deltaT())) < 0
				if (!nextInside && testLine.containsPoint(edge.b)) {
					if (logIS(edge.b)) return PointVsFace.ON_EDGE
				}
			} else {
				for (const edgeT of edge.edgeISTsWithPlane(testPlane)) {
					if (edgeT == edge.bT) {
						if (!testLine.containsPoint(edge.b)) continue
						// endpoint lies on intersection testLine
						const edgeInside = dotCurve2(edge.curve, edge.bT, lineOut, -sign(edge.deltaT())) < 0
						const nextInside =
							colinearEdges[nextEdgeIndex] ||
							dotCurve2(nextEdge.curve, nextEdge.aT, lineOut, sign(nextEdge.deltaT())) < 0
						if (edgeInside != nextInside) {
							if (logIS(edge.b)) return PointVsFace.ON_EDGE
						}
					} else if (edgeT != edge.aT) {
						const p = edge.curve.at(edgeT)
						if (!testLine.containsPoint(p)) continue
						// edge crosses testLine, neither starts nor ends on it
						if (logIS(p)) return PointVsFace.ON_EDGE
						// TODO: tangents?
					}
				}
			}
		}
		return inside ? PointVsFace.INSIDE : PointVsFace.OUTSIDE
	}

	toString(): string {
		return this.toSource()
	}

	toSource(rounder: (x: number) => number = x => x): string {
		return callsce.call(undefined, 'new ' + this.constructor.name, ...this.getConstructorParameters())
	}

	/**
	 * Return points which would touch AABB. Doesnt include borders due to paramtetric bounds, for example.
	 */
	getExtremePoints(): V3[] {
		return []
	}

	abstract normalP(p: V3): V3

	abstract getConstructorParameters(): any[]

	abstract isTsForLine(line: L3): number[]

	/**
	 * TODO remove constraint
	 * IMPORTANT: The tangents of the resulting curves need to be equal to the cross product of this and surface in the
	 * point. I.e.: for every point p p on a returned curve: curve.tangentAt(curve.pointT(p)) == this.normalP(p)
	 * X surface.normalP(p)
	 *
	 * Cross product is not commutative, so curve.tangentAt(curve.pointT(p)) == surface.normalP(p) X
	 * this.normalP(p) is not valid.
	 */
	abstract isCurvesWithPlane(plane: P3): Curve[]

	isCurvesWithSurface(surface: Surface): Curve[] {
		return surface.isCurvesWithSurface(this) //.map(curve => curve.reversed())
	}

	containsCurve(curve: Curve): boolean {
		if (curve instanceof PICurve) {
			// if (this.equals(curve.parametricSurface) || this.equals(curve.implicitSurface)) {
			// 	return true
			// }
		}
		if (curve instanceof PPCurve) {
			if (this.equals(curve.parametricSurface1) || this.equals(curve.parametricSurface2)) {
				return true
			}
		}
		if (curve instanceof ImplicitCurve) {
			for (let i = ceil(curve.tMin) + 1; i <= floor(curve.tMax) - 1; i++) {
				if (!this.containsPoint(curve.points[i])) {
					return false
				}
			}
			return true
		}
		return false
	}

	abstract containsPoint(pWC: V3): boolean

	abstract flipped(): this

	flipped2(doFlip: boolean): this {
		return doFlip ? this.flipped() : this
	}

	abstract loopContainsPoint(contour: Edge[], point: V3): PointVsFace

	/**
	 * Returns true iff the surface occupies the same space as the argument (not necessarily same normal1)
	 */
	abstract isCoplanarTo(surface: Surface): boolean

	/**
	 * coplanar and same normals
	 */
	abstract like(object: any): boolean

	abstract edgeLoopCCW(loop: Edge[]): boolean

	clipCurves(curves: Curve[]): Curve[] {
		return curves
	}

	equals(obj: any): boolean {
		return (
			this === obj ||
			(this.constructor === obj.constructor &&
				this.getConstructorParameters().equals(obj.getConstructorParameters()))
		)
	}

	hashCode(): int {
		return this.getConstructorParameters().hashCode()
	}

	zDirVolume(allEdges: Edge[]): { centroid: V3; volume: number } {
		return this.visit(ZDirVolumeVisitor, allEdges)
	}

	calculateArea(allEdges: Edge[]): number {
		return this.visit(CalculateAreaVisitor as any, allEdges)
	}
}

export enum PointVsFace {
	INSIDE,
	OUTSIDE,
	ON_EDGE,
}

export abstract class ImplicitSurface extends Surface {
	static is(obj: any): obj is ImplicitSurface {
		return obj.implicitFunction && obj.didp
	}

	abstract implicitFunction(): (pWC: V3) => number

	/**
	 * partial derivatives of this.implicitFunction in point pWC
	 * @param pWC
	 * @return
	 */
	abstract didp(pWC: V3): V3
}
