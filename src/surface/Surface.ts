abstract class Surface extends Transformable implements Equalable {
	toString(): string {
		return this.toSource()
	}

	toSource(rounder: (x: number) => number = x => x): string {
		return callsce.call(undefined, 'new ' + this.constructor.name, ...this.getConstructorParameters())
	}

	abstract getConstructorParameters(): any[]

	abstract transform(m4: M4, desc?: string): Surface

	abstract isTsForLine(line: L3): number[]

	/**
	 * IMPORTANT: The tangents of the resulting curves need to be equal to the cross product of this and surface in the
	 * point. I.e.: for every point p p on a returned curve: curve.tangentAt(curve.pointT(p)) == this.normalAt(p)
	 * X surface.normalAt(p)
	 *
	 * Cross product is not commutative, so curve.tangentAt(curve.pointT(p)) == surface.normalAt(p) X
	 * this.normalAt(p) is not valid.
	 *
	 */
	abstract isCurvesWithPlane(plane: P3): Curve[]

	isCurvesWithSurface(surface: Surface): Curve[] {
		return surface.isCurvesWithSurface(this).map(curve => curve.reversed())
	}

	abstract containsCurve(curve: Curve): boolean

	abstract containsPoint(p: V3): boolean

	abstract flipped(): Surface

	flipped2(doFlip: boolean): Surface {
		return doFlip ? this.flipped() : this
	}

	normalAt(p: V3): V3 {
		const pmPoint = this.pointToParameterFunction()(p)
		return this.parametricNormal()(pmPoint.x, pmPoint.y)
	}

	abstract loopContainsPoint(contour: Edge[], point: V3): PointVsFace

	/**
	 * Returns true iff the surface occupies the same space as the argument (not necessarily same normal1)
	 */
	abstract isCoplanarTo(surface: Surface): boolean

	/**
	 * coplanar and same normals
	 */
	abstract like(object): boolean



	abstract edgeLoopCCW(loop: Edge[]): boolean




	static loopContainsPointGeneral(loop: Edge[], p: V3, testLine: L3, lineOut: V3): PointVsFace {
		const testPlane = P3.normalOnAnchor(lineOut, p)
		// edges colinear to the testing line; these will always be counted as "inside" relative to the testing line
		const colinearEdges = loop.map((edge) => edge.colinearToLine(testLine))
		let inside = false

		function logIS(isP) {
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
				if (nextInside) {
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

	abstract equals(obj: any): boolean
	abstract hashCode(): int {
		return 43
	}

	abstract zDirVolume(allEdges: Edge[]): {centroid: V3, volume: number}

	abstract calculateArea(allEdges: Edge[]): number

	static genericCurvesIP(ps: ParametricSurface, is: ImplicitSurface): Curve[] {
		const psf = ps.parametricFunction()
		const isf = is.implicitFunction()
		const ic = (s, t) => isf(psf(s, t))
		return Curve.breakDownIC(ic, ps.sMin, ps.sMax, ps.tMin, ps.tMax, 0.1, 0.1, 0.0001)
	}
}
enum PointVsFace {INSIDE, OUTSIDE, ON_EDGE}