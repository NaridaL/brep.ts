class PlaneSurface extends Surface {
    plane: P3
    right: V3
    up: V3

    constructor(plane: P3, right?: V3, up?: V3) {
        super()
        assertInst(P3, plane)
        this.plane = plane
        this.up = up || plane.normal.getPerpendicular().unit()
        this.right = right || this.up.cross(this.plane.normal).unit()
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

    edgeLoopCCW(contour: Edge[]): boolean {
        return isCCW(contour.map(edge => edge.points()).concatenated(), this.plane.normal)
	    let totalAngle = 0
	    for (let i = 0; i < contour.length; i++) {
		    const ipp = (i + 1) % contour.length
		    const edge = contour[i], nextEdge = contour[ipp]
		    assert(edge.b.like(nextEdge.a), "edges dont form a loop")
		    if (edge.curve instanceof SemiEllipseCurve) {
			    totalAngle += edge.rotViaPlane(this.plane.normal)
			    // console.log(edge.toString(), edge.rotViaPlane(this.plane.normal))
		    }
		    totalAngle += edge.bDir.angleRelativeNormal(nextEdge.aDir, this.plane.normal)
	    }
	    const result = totalAngle > 0
        const result2 = PlaneFace.prototype.calculateArea.apply({surface: this, contour: contour}).area > 0
        //assert (result == result2)
	    return result2
    }

    loopContainsPoint(loop: Edge[], p: V3): PointVsFace {
        const dir = this.right.plus(this.up.times(0.123)).unit()
        const line = new L3(p, dir)
        const lineOut = dir.cross(this.plane.normal)
        return Surface.loopContainsPointGeneral(loop, p, line, lineOut)
    }

    pointToParameterFunction() {
        const matrix = M4.forSys(this.right, this.up, this.normal, this.plane.anchor)
        const matrixInverse = matrix.inversed()
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
        return this.plane.containsCurve(curve)
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