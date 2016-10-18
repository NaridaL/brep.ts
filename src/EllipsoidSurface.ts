class EllipsoidSurface extends Surface {
    center: V3
    f1: V3
    f2: V3
    f3: V3
    matrix: M4
    inverseMatrix: M4
    matrixTransposedInverse: M4


    constructor(center: V3, f1: V3, f2: V3, f3: V3) {
        super()
        assertVectors(center, f1, f2, f3)
        this.center = center
        this.f1 = f1
        this.f2 = f2
        this.f3 = f3
        this.matrix = M4.forSys(f1, f2, f3, center)
        this.inverseMatrix = this.matrix.inversed()
        this.matrixTransposedInverse = this.matrix.inversed().transposed()
    }

    toSource() {
        return `new SphereSurface(${this.center.toSource()}, ${this.f1.toSource()}, ${this.f2.toSource()}, ${this.f3.toSource()})`
    }

    isTsForLine(line) {
        assertInst(L3, line)
        // transforming line manually has advantage that dir1 will not be renormalized,
        // meaning that calculated values t for localLine are directly transferable to line
        var localAnchor = this.inverseMatrix.transformPoint(line.anchor)
        var localDir = this.inverseMatrix.transformVector(line.dir1)
        return EllipsoidSurface.unitISTsWithLine(localAnchor, localDir)
    }

    isCoplanarTo(surface) {
        if (this == surface) return true
        if (surface.constructor != EllipsoidSurface) return false
        if (!this.center.like(surface.center)) return false
        if (this.isSphere()) return surface.isSphere() && NLA.eq(this.f1.length(), this.f2.length())

        let thisMA = this.mainAxes(), surfaceMA = surface.mainAxes()
        return thisMA.every(tma => surfaceMA.some(sma => tma.like(sma)))
    }

    containsEllipse(ellipse: EllipseCurve):boolean {
        let localEllipse = ellipse.transform(this.inverseMatrix)
        let distLocalEllipseCenter = localEllipse.center.length()
        let correctRadius = Math.sqrt(1 - distLocalEllipseCenter * distLocalEllipseCenter)
        return NLA.lt(distLocalEllipseCenter, 1) && localEllipse.isCircular() && localEllipse.f1.hasLength(correctRadius)
    }

    containsCurve(curve) {
        if (curve instanceof EllipseCurve) {
            return this.containsEllipse(curve)
        } else {
            return false
        }
    }

    transform(m4) {
        return new EllipsoidSurface(
            m4.transformPoint(this.center),
            m4.transformVector(this.f1),
            m4.transformVector(this.f2),
            m4.transformVector(this.f3)) as this
    }

    isInsideOut():boolean {
        return this.f1.cross(this.f2).dot(this.f3) < 0
    }

    /**
     *
     * @returns {EllipsoidSurface}
     */
    flipped() {
        return new EllipsoidSurface(
            this.center,
            this.f1,
            this.f2,
            this.f3.negated())
    }


    toMesh(subdivisions) {
        return GL.Mesh.sphere(subdivisions).transform(this.matrix)
        // let mesh = new GL.Mesh({triangles: true, lines: false, normals: true})
        // let pf = this.parametricFunction()
        // let pn = this.parametricNormal()
        // let aCount = 32, bCount = 16, vTotal = aCount * bCount
        // for (let i = 0, a = -PI; i < aCount; i++, a += 2 * PI / aCount) {
        // 	for (let j = 0, b = -Math.PI / 2; j < bCount; j++, b += Math.PI / (bCount - 1)) {
        // 		mesh.vertices.push(pf(a, b))
        // 		mesh.normals.push(pn(a, b))
        // 		j != (bCount - 1) && pushQuad(mesh.triangles, true,
        // 			i * bCount + j, i * bCount + j + 1,
        // 			((i + 1) * bCount + j) % vTotal, ((i + 1) * bCount + j + 1) % vTotal)
        // 	}
        // }
        // mesh.compile()
        // return mesh
    }

    parametricNormal() {
        // ugh
        // paramtric ellipsoid point q(a, b)
        // normal == (dq(a, b) / da) X (dq(a, b) / db) (Cross product of partial derivatives
        // normal == cos b * (f2 X f3 * cos b * cos a + f3 X f1 * cos b * sin a + f1 X f2 * sin b)
        return (a, b) => {
            let {f1, f2, f3} = this
            let normal = f2.cross(f3).times(Math.cos(b) * Math.cos(a))
                .plus(f3.cross(f1).times(Math.cos(b) * Math.sin(a)))
                .plus(f1.cross(f2).times(Math.sin(b)))
                //.times(Math.cos(b))
                .normalized()
            return normal
        }
    }

    parametricFunction() {
        // this(a, b) = f1 cos a cos b + f2 sin a cos b + f2 sin b
        return (alpha, beta) => {
            return this.matrix.transformPoint(V3.sphere(alpha, beta))
        }
    }

    pointToParameterFunction() {
        return (pWC:V3) => {
            const pLC = this.inverseMatrix.transformPoint(pWC)
            let alpha = pLC.angleXY()
            let beta = Math.asin(pLC.z)
            return new V3(alpha, beta, 0)
        }
    }

    isSphere():boolean {
        return NLA.eq(this.f1.length(), this.f2.length())
            && NLA.eq(this.f2.length(), this.f3.length())
            && NLA.eq(this.f3.length(), this.f1.length())
            && this.f1.isPerpendicularTo(this.f2)
            && this.f2.isPerpendicularTo(this.f3)
            && this.f3.isPerpendicularTo(this.f1)
    }

    implicitFunction() {
        return (pWC) => {
            let pLC = this.inverseMatrix.transformPoint(pWC)
            return pLC.length() - 1
        }
    }

    mainAxes(): {f1: V3, f2: V3, f3: V3} {
        // q(a, b) = f1 cos a cos b + f2 sin a cos b + f2 sin b
        // q(s, t, u) = s * f1 + t * f2 + u * f3 with s² + t² + u² = 1
        // (del q(a, b) / del a) = f1 (-sin a) cos b  + f2 cos a cos b
        // (del q(a, b) / del b) = f1 cos a (-sin b) + f2 sin a (-sin b) + f2 cos b
        // del q(s, t, u) / del a = -t f1 + s f2
        // (del q(a, b) / del a) DOT q(a, b) == 0
        // (f1 (-sin a) cos b  + f2 cos a cos b) DOT (f1 cos a cos b + f2 sin a cos b + f2 sin b) == 0
        // (del q(a, b) / del b) DOT q(a, b) == 0
        // (f1 cos a (-sin b) + f2 sin a (-sin b) + f2 cos b) DOT (f1 cos a cos b + f2 sin a cos b + f2 sin b) == 0

        // Solve[
        // (f1 (-sin a) cos b  + f2 cos a cos b) * (f1 cos a cos b + f2 sin a cos b + f2 sin b) = 0,
        // (f1 cos a (-sin b) + f2 sin a (-sin b) + f2 cos b) * (f1 cos a cos b + f2 sin a cos b + f2 sin b) = 0}, a, b]
        const {f1, f2, f3} = this
        const f = (s, t) => f1
        assert(false)
    }

    containsPoint(p) {
        return NLA.eq0(this.implicitFunction()(p))
    }

    boundsFunction() {
        assert(false)
    }

    /**
     * unit sphere: x² + y² + z² = 1
     * line: p = anchor + t * dir |^2
     * p² = (anchor + t * dir)^2
     * 1 == (anchor + t * dir)^2
     * 1 == anchor DOT anchor + 2 * anchor * t * dir + t² * dir DOT dir
     */
    static unitISTsWithLine(anchor: V3, dir: V3):number[] {
        // for 0 = a t² + b t + c
        let a = dir.dot(dir)
        let b = 2 * anchor.dot(dir)
        let c = anchor.dot(anchor) - 1
        return pqFormula(b / a, c / a)
    }

    /**
     * unit sphere: x² + y² + z² = 1
     * plane: normal DOT p = w
     */
    static unitISCurvesWithPlane(plane:P3):EllipseCurve[] {
        let distPlaneCenter = Math.abs(plane.w)
        if (NLA.lt(distPlaneCenter, 1)) {
            // result is a circle
            // radius of circle: imagine right angled triangle (origin -> center of intersection circle -> point on intersection circle)
            // pythagoras: 1² == distPlaneCenter² + isCircleRadius² => isCircleRadius == sqrt(1 - distPlaneCenter²)
            let isCircleRadius = Math.sqrt(1 - distPlaneCenter * distPlaneCenter)
            let center = plane.anchor
            let f1 = plane.normal.getPerpendicular().toLength(isCircleRadius)
            let f2 = plane.normal.cross(f1)
            return [new EllipseCurve(plane.anchor, f1, f2)]
        } else {
            return []
        }
    }

    static sphere(radius: number, center?: V3): EllipsoidSurface {
        assertNumbers(radius)
        center && assertVectors(center)
        return new EllipsoidSurface(center || V3.ZERO, new V3(radius, 0, 0), new V3(0, radius, 0), new V3(0, 0, radius))
    }

    /**
     * x²/a² + y²/b² + z²/c² = 1
     */
    static forABC(a: number, b: number, c: number, center?: V3): EllipsoidSurface {
        return new EllipsoidSurface(center || V3.ZERO, new V3(a, 0, 0), new V3(0, b, 0), new V3(0, 0, c))
    }

    volume(): number {
        return 4 / 4 * Math.PI * this.f1.dot(this.f2.cross(this.f3))
    }
}
