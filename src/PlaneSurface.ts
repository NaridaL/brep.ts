class PlaneSurface extends Surface {
    plane: P3
    right: V3
    up: V3

    constructor(plane: P3, right?: V3, up?: V3) {
        super()
        assertInst(P3, plane)
        this.plane = plane
        this.up = up || plane.normal.getPerpendicular().normalized()
        this.right = right || this.up.cross(this.plane.normal).normalized()
        assert(this.right.cross(this.up).like(this.plane.normal))
    }

    isCoplanarTo(surface) {
        return surface instanceof PlaneSurface && this.plane.isCoplanarToPlane(surface.plane)
    }

    like(surface) {
        return surface instanceof PlaneSurface && this.plane.like(surface.plane)
    }

    parametricFunction() {
        var matrix = M4.forSys(this.right, this.up, this.plane.normal, this.plane.anchor)
        return function (s, t) {
            return matrix.transformPoint(new V3(s, t, 0))
        }
    }

    implicitFunction() {
        return p => this.plane.distanceToPointSigned(p)
    }

    isCurvesWithISurface(implicitSurface) {
        assert(implicitSurface.implicitFunction, 'implicitSurface.implicitFunction')
        return new CurvePI(this, implicitSurface)
    }

    isCurvesWithSurface(surface2: Surface): Curve[] {
        assert(false)
        return null
    }

    edgeLoopCCW(contour: Edge[]):boolean {
        var totalAngle = 0
        for (var i = 0; i < contour.length; i++) {
            var ipp = (i + 1) % contour.length
            var edge = contour[i], nextEdge = contour[ipp]
            assert(edge.b.like(nextEdge.a), "edges dont form a loop")
            if (edge.curve instanceof EllipseCurve) {
                totalAngle += edge.rotViaPlane(this.plane.normal)
                // console.log(edge.toString(), edge.rotViaPlane(this.plane.normal))
            }
            totalAngle += edge.bDir.angleRelativeNormal(nextEdge.aDir, this.plane.normal)
        }
        return totalAngle > 0
    }

    loopContainsPoint(loop: Edge[], p: V3): PointVsFace {
        assert(loop)
        assertVectors(p)
        const dir = this.right.plus(this.up.times(0.123)).normalized()
        const line = new L3(p, dir)
	    const lineOut = dir.cross(this.plane.normal)
        const plane2 = P3.normalOnAnchor(lineOut, p)
        const colinearEdges = loop.map((edge) => edge.colinearToLine(line))
        const colinearEdgeInside = loop.map((edge, i) => colinearEdges[i] && edge.aDir.dot(dir) > 0)
        let inside = false

        function logIS(p): boolean {
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
                        // endpoint lies on intersection line
	                    const edgeInside = dotCurve(lineOut, edge.bDir, edge.bDDT)
	                    const nextInside = colinearEdges[nextEdgeIndex]
		                    ? colinearEdgeInside[nextEdgeIndex]
		                    : dotCurve(lineOut, nextEdge.aDir, nextEdge.aDDT)
	                    if (edgeInside != nextInside) {
		                    if (logIS(edge.b)) return PointVsFace.ON_EDGE
	                    }
                    } else if (edgeT != edge.aT) {
                        // edge crosses line, neither starts nor ends on it
	                    if (logIS(edge.curve.at(edgeT))) return PointVsFace.ON_EDGE// TODO: tangents?
                    }
                }
            }
        }
        return inside ? PointVsFace.INSIDE : PointVsFace.OUTSIDE

    }

    pointToParameterFunction() {
        var matrix = M4.forSys(this.right, this.up, this.normal, this.plane.anchor)
        var matrixInverse = matrix.inversed()
        return function (pWC) {
            return matrixInverse.transformPoint(pWC)
        }
    }

    normalAt(p) {
        return this.plane.normal
    }

    containsPoint(p) {
        return this.plane.containsPoint(p)
    }

    containsCurve(curve: Curve): boolean {
        if (curve instanceof L3) {
            return this.plane.containsLine(curve)
        } else if (curve instanceof EllipseCurve || curve instanceof HyperbolaCurve || curve instanceof ParabolaCurve) {
            return this.plane.containsPoint(curve.center) && this.plane.normal.isParallelTo(curve.normal)
        } else if (curve instanceof BezierCurve) {
            return curve.points.every(p => this.plane.containsPoint(p))
        } else {
            throw new Error(curve)
        }
    }

    transform(m4): this {
        return new PlaneSurface(this.plane.transform(m4)) as this
    }

    flipped() {
        return new PlaneSurface(this.plane.flipped(), this.right, this.up.negated())
    }

    toString() {
        return this.plane.toString()
    }

    toSource() {
        return `new PlaneSurface(${this.plane})`
    }

    static throughPoints(a, b, c): PlaneSurface {
        return new PlaneSurface(P3.throughPoints(a, b, c))
    }
}
NLA.registerClass(PlaneSurface)