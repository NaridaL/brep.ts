import {Equalable} from 'javasetmap.ts'
import {int, Tuple3, Tuple2, Tuple4,AABB,DEG,EllipticE,EllipticF,GOLDEN_RATIO,M4,MINUS,Matrix,NLA_DEBUG,NLA_PRECISION,P3XY,P3YZ,P3ZX,SCE,STR,TAU,Transformable,V,V3,Vector,addOwnProperties,arithmeticGeometricMean,arrayCopy,arrayCopyBlocks,arrayCopyStep,arrayFromFunction,arrayRange,arraySwap,assert,assertInst,assertNever,assertNumbers,assertVectors,assertf,between,bisect,callsce,canonAngle,ceil10,checkDerivate,clamp,combinations,defaultRoundFunction,disableConsole,doubleSignedArea,enableConsole,eq,eq0,eq02,eq2,eqAngle,floatHashCode,floor10,forceFinite,fuzzyBetween,fuzzyUniques,fuzzyUniquesF,gaussLegendre24Weights,gaussLegendre24Xs,gaussLegendreQuadrature24,ge,getIntervals,getRoots,glq24_11,glqInSteps,gt,hasConstructor,isCCW,le,lerp,lt,mapPush,midpointRuleQuadrature,mod,newtonIterate,newtonIterate1d,newtonIterate2d,newtonIterate2dWithDerivatives,newtonIterateSmart,newtonIterateWithDerivative,numberToStr,pqFormula,rad2deg,randomColor,repeatString,round10,snap,snap0,snap2,snapEPS,solveCubicReal2,time,zeroAngle} from 'ts3dutils'
import {Mesh} from 'tsgl'

import {Curve, P3, SemiCylinderSurface, ProjectedCurveSurface, L3, ISInfo, ParametricSurface, ImplicitSurface,
    Edge,
    HyperbolaCurve,
    SemiEllipseCurve,
    ParabolaCurve,
    EllipseCurve} from '../index'

const {PI, cos, sin, min, max, tan, sign, ceil, floor, abs, sqrt, pow, atan2, round} = Math

export abstract class Surface extends Transformable implements Equalable {
	toString(): string {
		return this.toSource()
	}

	toSource(rounder: (x: number) => number = x => x): string {
		return callsce.call(undefined, 'new ' + this.constructor.name, ...this.getConstructorParameters())
	}

	abstract normalP(p: V3): V3

	abstract getConstructorParameters(): any[]

	abstract isTsForLine(line: L3): number[]

	/**
	 * IMPORTANT: The tangents of the resulting curves need to be equal to the cross product of this and surface in the
	 * point. I.e.: for every point p p on a returned curve: curve.tangentAt(curve.pointT(p)) == this.normalP(p)
	 * X surface.normalP(p)
	 *
	 * Cross product is not commutative, so curve.tangentAt(curve.pointT(p)) == surface.normalP(p) X
	 * this.normalP(p) is not valid.
	 */
	abstract isCurvesWithPlane(plane: P3): Curve[]

	isCurvesWithSurface(surface: Surface): Curve[] {
		return surface.isCurvesWithSurface(this).map(curve => curve.reversed())
	}

	containsCurve(curve: Curve): boolean {
		if (curve instanceof ImplicitCurve) {
			for (let i = ceil(curve.tMin); i <= floor(curve.tMax); i++) {
				if (!this.containsPoint(curve.points[i])) {
					return false
				}
			}
			return true
		} else {
			return false
		}
	}

	abstract containsPoint(pWC: V3): boolean

	abstract flipped(): this

	flipped2<T extends Surface>(this: T, doFlip: boolean): T {
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


	static loopContainsPointGeneral(loop: Edge[], p: V3, testLine: L3, lineOut: V3): PointVsFace {
		const testPlane = P3.normalOnAnchor(lineOut, p)
		// edges colinear to the testing line; these will always be counted as "inside" relative to the testing line
		const colinearEdges = loop.map((edge) => edge.colinearToLine(testLine))
		let inside = false

		function logIS(isP: V3) {
			const isT = testLine.pointT(isP)
			if (eq0(isT)) {
				return true
			} else if (isT > 0) {
				inside = !inside
			}
		}

		for (let edgeIndex = 0; edgeIndex < loop.length; edgeIndex++) {
			const edge = loop[edgeIndex]
			const nextEdgeIndex = (edgeIndex + 1) % loop.length, nextEdge = loop[nextEdgeIndex]
			//console.log(edge.toSource()) {p:V(2, -2.102, 0),
			if (colinearEdges[edgeIndex]) {
				const lineAT = testLine.pointT(edge.a), lineBT = testLine.pointT(edge.b)
				if (Math.min(lineAT, lineBT) <= NLA_PRECISION && -NLA_PRECISION <= Math.max(lineAT, lineBT)) {
					return PointVsFace.ON_EDGE
				}
				// edge colinear to intersection
				const nextInside = colinearEdges[nextEdgeIndex] || dotCurve(lineOut, nextEdge.aDir, nextEdge.aDDT) < 0
				if (!nextInside) {
					if (logIS(edge.b)) return PointVsFace.ON_EDGE
				}
			} else {
				for (const edgeT of edge.edgeISTsWithPlane(testPlane)) {
					if (edgeT == edge.bT) {
						if (!testLine.containsPoint(edge.b)) continue
						// endpoint lies on intersection line
						if (edge.b.like(p)) {
							// TODO: refactor, dont check for different sides, just logIs everything
							return PointVsFace.ON_EDGE
						}
						const edgeInside = dotCurve(lineOut, edge.bDir, edge.bDDT) > 0
						const nextInside = colinearEdges[nextEdgeIndex] || dotCurve(lineOut, nextEdge.aDir, nextEdge.aDDT) < 0
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

	clipCurves(curves: Curve[]): Curve[] {
		return curves
	}

	abstract equals(obj: any): boolean

	hashCode(): int {
		return this.getConstructorParameters().hashCode()
	}

	zDirVolume(allEdges: Edge[]): {centroid: V3, volume: number} {
		return this.visit(ZDirVolumeVisitor, allEdges)
	}

	calculateArea(allEdges: Edge[]): number {
		return this.visit(CalculateAreaVisitor, allEdges)
	}
}
export enum PointVsFace {INSIDE, OUTSIDE, ON_EDGE}
