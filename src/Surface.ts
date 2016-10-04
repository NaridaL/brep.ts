abstract class Surface extends Transformable {
    toString(): string {
        return this.toSource()
    }

    abstract toSource(): string

    abstract isTsForLine(line: L3): number[]

    /**
     * IMPORTANT: The tangents of the resulting curves need to be equal to the cross product of this and surface in the point.
     * I.e.: for every point p p on a returned curve: curve.tangentAt(curve.pointLambda(p)) == this.normalAt(p) X surface.normalAt(p)
     *
     * Cross product is not commutative, so curve.tangentAt(curve.pointLambda(p)) == surface.normalAt(p) X this.normalAt(p)
     * is not valid.
     *
     */
    abstract isCurvesWithPlane(plane: P3): Curve[]

    abstract isCurvesWithSurface(surface: Surface): Curve[]

    abstract containsCurve(curve: Curve): boolean

    abstract containsPoint(p: V3): boolean

    abstract flipped(): Surface

    normalAt(p: V3) : V3 {
        var pmPoint = this.pointToParameterFunction()(p)
        return this.parametricNormal()(pmPoint.x, pmPoint.y)
    }

    abstract edgeLoopContainsPoint(contour: Edge[], point: V3): boolean

    /**
     * Returns true iff the surface occupies the same space as the argument (not necessarily same normal)
     */
    abstract isCoplanarTo(surface: Surface): boolean

    abstract like(object) :boolean

    abstract pointToParameterFunction?(): (pWC: V3, hint?: number)=>V3
    abstract parametricFunction?(): (s:number, t:number)=>V3
    abstract parametricNormal?(): (s:number, t:number)=>V3
}