/**
 * Surface normal is (t, z) => this.baseCurve.tangentAt(t) X this.dir1
 * Choose dir1 appropriately to select surface orientation.
 */
class ProjectedCurveSurface extends Surface {
    baseCurve: BezierCurve
    dir1: V3
    tMin: number
    tMax: number
    zMin: number
    zMax: number


    constructor(baseCurve, dir1, tMin = baseCurve.tMin, tMax = baseCurve.tMax, zMin = -Infinity, zMax = Infinity) {
        super()
        assertInst(Curve, baseCurve)
        assertInst(V3, dir1)
        assertNumbers(tMin, tMax)
        assert(dir1.hasLength(1), 'dir must be unit' + dir1)
        assert(tMin < tMax)
        assert(zMin < zMax)
        this.baseCurve = baseCurve
        this.dir1 = dir1
        this.tMin = tMin
        this.tMax = tMax
        this.zMin = zMin
        this.zMax = zMax
    }

    boundsFunction() {
        return (t, z) => {
            return this.tMin <= t && t <= this.tMax
        }
    }

    zDirVolume(allEdges: Edge[]): {centroid: V3, volume: number} {
        assertNever()
    }

    toSource(): string {
        return `new ProjectedCurveSurface(${this.baseCurve}, ${this.dir1}, ${this.tMin}, ${this.tMax}, ${this.zMin}, ${this.zMax})`
    }

    toString(): string {
        return this.toSource()
    }


    toMesh(tStart = this.tMin, tEnd = this.tMax, zStart = this.zMin, zEnd = this.zMax, count = 128) {

        const tInterval = tEnd - tStart, tStep = tInterval / (count - 1)
        const baseVertices =
            NLA.arrayFromFunction(count, i => this.baseCurve.at(tStart + i * tStep).plus(this.dir1.times(zStart)))
        const normalFunc = this.parametricNormal()
        const normals = NLA.arrayFromFunction(count, i => normalFunc(tStart + i * tStep, 0))

        return GL.Mesh.offsetVertices(baseVertices, this.dir1.times(zEnd - zStart), false, normals)
    }

    parametricFunction() {
        return (t: number, z: number) => this.baseCurve.at(t).plus(this.dir1.times(z))
    }

    parametricNormal() {
        return (t: number, z: number) => this.baseCurve.tangentAt(t).cross(this.dir1).unit()
    }

    implicitFunction() {

    }

    parametersValid(t: number, z: number) {
        return NLA.between(t, this.tMin, this.tMax) && NLA.between(z, this.zMin, this.zMax)
    }

    footParameters(pWC, ss, st) {
        const basePlane = new P3(this.dir1, 0)
        const projCurve = this.baseCurve.project(basePlane)
        const projPoint = basePlane.projectedPoint(pWC)
        const t = projCurve.closestTToPoint(projPoint, undefined, undefined, ss)
        const z = pWC.minus(this.baseCurve.at(t)).dot(this.dir1)
        return new V3(t, z, 0)
    }

    pointToParameterFunction() {
        const projPlane = new P3(this.dir1, 0)
        const dir1 = this.dir1, baseCurve = this.baseCurve
        const projBaseCurve = baseCurve.project(projPlane)
        const _this = this
        return function (pWC) {
            const projPoint = projPlane.projectedPoint(pWC)
            const t = projBaseCurve.pointT(projPoint)
            const z = pWC.minus(baseCurve.at(t)).dot(dir1)
            return new V3(t, z, 0)
        }

        // const baseCurve = this.baseCurve, dir1 = this.dir1
        // const tMin = this.tMin, tMax = this.tMax
        // return function (pWC) {
        // 	const iterationFunc = t => {const d = baseCurve.at(t).minus(pWC); return d.dot(dir1) + d.length() }
        // 	const startT = NLA.arrayFromFunction(16, i => tMin + (tMax - tMin) * i / 15)
        // 		.withMax(t => -abs(iterationFunc(t)))
        // 	const t = newtonIterate1d(iterationFunc, startT, 16)
        // 	const z = pWC.minus(baseCurve.at(t)).dot(dir1)
        // 	return V(t, z, 0)
        // }
    }

    isCurvesWithPlane(plane): Curve[] {
        assertInst(P3, plane)
        if (this.dir1.isPerpendicularTo(plane.normal)) {

            const ts = this.baseCurve.isTsWithPlane(plane)
            return ts.map(t => {
                const l3dir = 0 < this.baseCurve.tangentAt(t).dot(plane.normal)
                    ? this.dir1
                    : this.dir1.negated()
                return new L3(this.baseCurve.at(t), l3dir)
            })
        } else {
            const projCurve = this.baseCurve.transform(M4.projection(plane, this.dir1))
            if (this.dir1.dot(plane.normal) > 0) {
                // we need to flip the ellipse so the tangent is correct
                console.log('FLIPPING')
                projCurve = projCurve.reversed()
            }
            return [projCurve]
        }
    }

    isCurvesWithSurface(surface) {
        if (surface instanceof PlaneSurface) {
            return this.isCurvesWithPlane(surface.plane)
        }
        if (surface instanceof ProjectedCurveSurface || surface instanceof SemiCylinderSurface) {
            const dir1 = surface instanceof ProjectedCurveSurface ? surface.dir1 : surface.dir.unit()
            if (this.dir1.isParallelTo(dir1)) {
                const otherCurve = surface instanceof ProjectedCurveSurface ? surface.baseCurve : surface.baseEllipse
                const infos = this.baseCurve.isInfosWithCurve(otherCurve)
                return infos.map(info => new L3(info.p, dir1))
            }
            if (surface instanceof ProjectedCurveSurface) {
                const line = new L3(this.baseCurve.at(0.5), this.dir1)
                const startPoint = line.at(surface.isTsForLine(line)[0])
                drVs.push({anchor: this.baseCurve.at(0.5), dir: this.dir1})
                console.log(startPoint)
                return [new PPCurve(this, surface, startPoint)]
                // const testVector = this.dir1.cross(surface.dir1).unit()
                // // look for points on surface.baseCurve where tangent DOT testVector == 0
                // const abcd1 = surface.baseCurve.tangentCoefficients().map(c => c.dot(testVector))
                // const ts1 = solveCubicReal2.apply(undefined, abcd1).concat(surface.tMin, surface.tMax)
                // const abcd2 = this.baseCurve.tangentCoefficients().map(c => c.dot(testVector))
                // const ts2 = solveCubicReal2.apply(undefined, abcd2)
                // const tt1 = ts1.map(t => surface.baseCurve.at(t).dot(testVector))
                // const tt2 = ts1.map(t => surface.baseCurve.at(t).dot(testVector))
                // console.log(ts1, ts2, tt1, tt2)
                // ts1.forEach(t => drPs.push(surface.baseCurve.at(t)))
                // ts2.forEach(t => drPs.push(this.baseCurve.at(t)))
                // return
            }
        }
        if (surface.parametricFunction) {
            return new PICurve(surface, this)
        } else if (surface.implicitFunction) {
            return new PICurve(this, surface)
        }
    }

    /**
     * @inheritDoc
     */
    containsPoint(p) {
        const uv = this.pointToParameterFunction()(p)
        return this.parametricFunction()(uv.x, uv.y).like(p)
    }

    /**
     * @inheritDoc
     */
    containsCurve(curve) {
        if (curve instanceof BezierCurve) {
            // project baseCurve and test curve onto a common plane and check if the curves are alike
            const projPlane = new P3(this.dir1, 0)
            const projBaseCurve = this.baseCurve.project(projPlane)
            const projCurve = curve.project(projPlane)

            return projBaseCurve.isColinearTo(projCurve)
        }
        if (curve instanceof L3) {
            return this.dir1.isParallelTo(curve.dir1) && this.containsPoint(curve.anchor)
        }
        return false
    }

    isCoplanarTo(surface) {
        return this == surface ||
            ProjectedCurveSurface == surface.constructor
            && this.dir1.isParallelTo(surface.dir1)
            && this.containsCurve(surface.baseCurve)
    }


    like(object) {
        if (!this.isCoplanarTo(object)) return false
        // normals need to point in the same direction (outwards or inwards) for both
        const p00 = this.parametricFunction()(0, 0)
        const thisNormal = this.parametricNormal()(0, 0)
        const otherNormal = object.normalAt(p00)
        return 0 < thisNormal.dot(otherNormal)
    }

    loopContainsPoint(loop: Edge[], p: V3): PointVsFace {
        assertVectors(p)
        assert(isFinite(p.x), p.y, p.z)
        const line = new L3(p, this.dir1)
        const ptpf = this.pointToParameterFunction()
        const pp = ptpf(p)
        if (isNaN(pp.x)) {
            console.log(this.sce, p.sce)
            assert(false)
        }
        const lineOut = this.baseCurve.tangentAt(pp.x).rejectedFrom(this.dir1)

        return Surface.loopContainsPointGeneral(loop, p, line, lineOut)
    }


    edgeLoopCCW(contour) {
        if (contour.length < 56) {
            var totalAngle = 0
            for (const i = 0; i < contour.length; i++) {
                const ipp = (i + 1) % contour.length
                const edge = contour[i], nextEdge = contour[ipp]
                totalAngle += edge.bDir.angleRelativeNormal(nextEdge.aDir, this.normalAt(edge.b))
            }
            return totalAngle > 0
        } else {
            const ptpF = this.pointToParameterFunction()
            return isCCW(contour.map(e => ptpF(e.a)), V3.Z)
        }
    }

    transform(m4: M4) {
        return new ProjectedCurveSurface(
            this.baseCurve.transform(m4),
            m4.transformVector(this.dir1).unit(),
            this.tMin, this.tMax, this.zMin, this.zMax) as this
    }

    isTsForLine(line) {
        assertInst(L3, line)
        let projectionPlane = new P3(this.dir1, 0)
        let projDir = projectionPlane.projectedVector(line.dir1)
        if (projDir.isZero()) {
            // line is parallel to this.dir
            return []
        }
        let projAnchor = projectionPlane.projectedPoint(line.anchor)
        let projBaseCurve = this.baseCurve.project(projectionPlane)
        return projBaseCurve
            .isInfosWithLine(projAnchor, projDir, this.tMin, this.tMax)
            .map(info => info.tOther)
    }


    flipped() {
        return new ProjectedCurveSurface(this.baseCurve, this.dir1.negated(), this.tMin, this.tMax)
    }
}
ProjectedCurveSurface.prototype.uStep = 1 / 4
ProjectedCurveSurface.prototype.vStep = 256
NLA.registerClass(ProjectedCurveSurface)