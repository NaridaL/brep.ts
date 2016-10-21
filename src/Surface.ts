abstract class Surface extends Transformable {
	toString(): string {
		return this.toSource()
	}

	abstract toSource(): string

	abstract isTsForLine(line: L3): number[]

	/**
	 * IMPORTANT: The tangents of the resulting curves need to be equal to the cross product of this and surface in the
	 * point. I.e.: for every point p p on a returned curve: curve.tangentAt(curve.pointLambda(p)) == this.normalAt(p)
	 * X surface.normalAt(p)
	 *
	 * Cross product is not commutative, so curve.tangentAt(curve.pointLambda(p)) == surface.normalAt(p) X
	 * this.normalAt(p) is not valid.
	 *
	 */
	abstract isCurvesWithPlane(plane: P3): Curve[]

	abstract isCurvesWithSurface(surface: Surface): Curve[]

	abstract containsCurve(curve: Curve): boolean

	abstract containsPoint(p: V3): boolean

	abstract flipped(): Surface

	normalAt(p: V3): V3 {
		const pmPoint = this.pointToParameterFunction()(p)
		return this.parametricNormal()(pmPoint.x, pmPoint.y)
	}

	abstract loopContainsPoint(contour: Edge[], point: V3): PointVsFace

	/**
	 * Returns true iff the surface occupies the same space as the argument (not necessarily same normal)
	 */
	abstract isCoplanarTo(surface: Surface): boolean

	abstract equals(object): boolean

	abstract like(object): boolean

	abstract pointToParameterFunction?(): (pWC: V3, hint?: number)=>V3

	abstract parametricFunction?(): (s: number, t: number)=>V3

	abstract parametricNormal?(): (s: number, t: number)=>V3

	abstract edgeLoopCCW(loop: Edge[]): boolean




	static loopContainsPointGeneral(loop: Edge[], p: V3, line: L3, lineOut: V3):PointVsFace {
		const plane2 = P3.normalOnAnchor(lineOut, p)
		const colinearEdges = loop.map((edge) => edge.colinearToLine(line))
		const colinearEdgeInside = loop.map((edge, i) => colinearEdges[i] && edge.aDir.dot(line.dir1) > 0)
		let inside = false

		function logIS(p) {
			const t = line.pointLambda(p)
			if (NLA.eq0(t)) {
				return true
			} else if (t > 0) {
				inside = !inside
			}
		}

		for (let edgeIndex = 0; edgeIndex < loop.length; edgeIndex++) {
			const edge = loop[edgeIndex]
			const nextEdgeIndex = (edgeIndex + 1) % loop.length, nextEdge = loop[nextEdgeIndex]
			//console.log(edge.toSource()) {p:V(2, -2.102, 0),
			if (colinearEdges[edgeIndex]) {
				const lineAT = line.pointLambda(edge.a), lineBT = line.pointLambda(edge.b)
				if (Math.min(lineAT, lineBT) <= NLA_PRECISION && -NLA_PRECISION <= Math.max(lineAT, lineBT)) {
					return PointVsFace.ON_EDGE
				}
				// edge colinear to intersection
				const nextInside = colinearEdges[nextEdgeIndex]
					? colinearEdgeInside[nextEdgeIndex]
					: dotCurve(lineOut, nextEdge.aDir, nextEdge.aDDT)
				if (colinearEdgeInside[edgeIndex] != nextInside) {
					if (logIS(edge.b)) return PointVsFace.ON_EDGE
				}
			} else {
				for (const edgeT of edge.edgeISTsWithPlane(plane2)) {
					if (edgeT == edge.bT) {
						if (!line.containsPoint(edge.b)) continue
						// endpoint lies on intersection line
						const edgeInside = dotCurve(lineOut, edge.bDir, edge.bDDT)
						const nextInside = colinearEdges[nextEdgeIndex]
							? colinearEdgeInside[nextEdgeIndex]
							: dotCurve(lineOut, nextEdge.aDir, nextEdge.aDDT)
						if (edgeInside != nextInside) {
							if (logIS(edge.b)) return PointVsFace.ON_EDGE
						}
					} else if (edgeT != edge.aT) {
						const p = edge.curve.at(edgeT)
						if (!line.containsPoint(p)) continue
						// edge crosses line, neither starts nor ends on it
						if (logIS(p)) return PointVsFace.ON_EDGE
						 // TODO: tangents?
					}
				}
			}
		}
		return inside ? PointVsFace.INSIDE : PointVsFace.OUTSIDE

	}
}

enum PointVsFace {INSIDE, OUTSIDE, ON_EDGE}