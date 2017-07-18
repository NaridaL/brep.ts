/**
 * Surface normal1 is (t, z) => this.baseCurve.tangentAt(t) X this.dir
 * Choose dir appropriately to select surface orientation.
 */
class ProjectedCurveSurface extends ParametricSurface {
    constructor(readonly baseCurve: Curve,
                readonly dir: V3,
                readonly sMin: number = baseCurve.tMin,
                readonly sMax: number = baseCurve.tMax,
                readonly tMin: number = -100,
                readonly tMax: number = 100) {
        super()
        assertInst(Curve, baseCurve)
        assertInst(V3, dir)
        assertNumbers(sMin, sMax, tMin, tMax)
        assert(sMin < sMax)
        assert(tMin < tMax)
    }

	'constructor': typeof ProjectedCurveSurface &
		{ new <T extends ProjectedCurveSurface>(baseCurve: Curve, dir: V3,
		                                        sMin: number, sMax: number, tMin: number, tMax: number): T }

	getConstructorParameters() {
		return [this.baseCurve, this.dir, this.sMin, this.sMax, this.tMin, this.tMax]
	}

	equals(obj: any): boolean {
		return this == obj ||
			Object.getPrototypeOf(this) == Object.getPrototypeOf(obj)
			&& this.dir.equals(obj.dir)
			&& this.baseCurve.equals(obj.baseCurve)
	}

	hashCode(): int {
		return [this.dir, this.baseCurve].hashCode()
	}

	containsLine(line: L3): boolean {
		return this.dir.isParallelTo(line.dir1) && this.containsPoint(line.anchor)
	}

	dpds(): (s: number, t: number) => V3 {
		return (s, t) => this.baseCurve.tangentAt(s)
	}

	dpdt(): (s: number, t: number) => V3 {
		return (s, t) => this.dir
	}

	normalST(s: number, t: number): V3 {
		return this.baseCurve.tangentAt(s).cross(this.dir).unit()
	}

	pST(s: number, t: number): V3 {
        return this.baseCurve.at(s).plus(this.dir.times(t))
    }

    pointFoot(pWC: V3, ss: number, st: number): V3 {
        const basePlane = new P3(this.dir, 0)
        const projCurve = this.baseCurve.project(basePlane)
        const projPoint = basePlane.projectedPoint(pWC)
        const t = projCurve.closestTToPoint(projPoint, ss)
        const z = pWC.minus(this.baseCurve.at(t)).dot(this.dir)
        return new V3(t, z, 0)
    }

	stPFunc(): (pWC: V3) => V3 {
        const projPlane = new P3(this.dir.unit(), 0)
	    const projBaseCurve = this.baseCurve.project(projPlane)
        return (pWC) => {
            const projPoint = projPlane.projectedPoint(pWC)
            const t = projBaseCurve.pointT(projPoint)
            const z = L3.pointT(this.baseCurve.at(t), this.dir, pWC)
            return new V3(t, z, 0)
        }
    }

    isCurvesWithPlane(plane: P3): Curve[] {
        assertInst(P3, plane)
        if (this.dir.isPerpendicularTo(plane.normal1)) {

            const ts = this.baseCurve.isTsWithPlane(plane)
            return ts.map(t => {
                const l3dir = 0 < this.baseCurve.tangentAt(t).dot(plane.normal1)
                    ? this.dir
                    : this.dir.negated()
                return new L3(this.baseCurve.at(t), l3dir.unit())
            })
        } else {
            let projCurve = this.baseCurve.transform(M4.project(plane, this.dir))
            if (this.dir.dot(plane.normal1) > 0) {
                // we need to flip the ellipse so the tangent is correct
                projCurve = projCurve.reversed()
            }
            return [projCurve]
        }
    }

    isCurvesWithSurface(surface: Surface): Curve[] {
        if (surface instanceof PlaneSurface) {
            return this.isCurvesWithPlane(surface.plane)
        }
        if (surface instanceof ProjectedCurveSurface || surface instanceof SemiCylinderSurface) {
            const dir1 = surface.dir
            if (this.dir.isParallelTo(dir1)) {
                const otherCurve = surface.baseCurve
                const infos = this.baseCurve.isInfosWithCurve(otherCurve)
                return infos.map(info => {
	                const correctDir = this.normalP(info.p).cross(surface.normalP(info.p))
	                return new L3(info.p, dir1.times(sign(correctDir.dot(dir1))))
                })
            }
            if (surface instanceof ProjectedCurveSurface) {
                const line = new L3(this.baseCurve.at(0.5), this.dir)
                const startPoint = line.at(surface.isTsForLine(line)[0])
                console.log(startPoint)
                return [new PPCurve(this, surface, startPoint)]
                // const testVector = this.dir.cross(surface.dir).unit()
                // // look for points on surface.baseCurve where tangent DOT testVector == 0
                // const abcd1 = surface.baseCurve.tangentCoefficients().map(c => c.dot(testVector))
                // const ts1 = solveCubicReal2.apply(undefined, abcd1).concat(surface.sMin, surface.sMax)
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
	    if (surface instanceof SemiEllipsoidSurface) {
		    return surface.isCurvesWithSurface(this)
	    }
	    assertNever()
    }

    containsPoint(pWC: V3): boolean {
        const uv = this.stPFunc()(pWC)
        return this.pSTFunc()(uv.x, uv.y).like(pWC)
    }

    containsCurve(curve: Curve): boolean {
	    if (curve instanceof L3) {
		    return this.dir.isParallelTo(curve.dir1) && this.containsPoint(curve.anchor)
	    }
	    if (curve instanceof PICurve) {
		    return super.containsCurve(curve)
	    }
        // project baseCurve and test curve onto a common plane and check if the curves are alike
        const projPlane = new P3(this.dir.unit(), 0)
        const projBaseCurve = this.baseCurve.project(projPlane)
        const projCurve = curve.project(projPlane)

        return projBaseCurve.isColinearTo(projCurve)
    }

    isCoplanarTo(surface: Surface): boolean {
        return this == surface ||
	        ((x): x is ProjectedCurveSurface => x.constructor == ProjectedCurveSurface)(surface)
	        //&& ProjectedCurveSurface == surface.constructor
	        //ProjectedCurveSurface.prototype == Object.getPrototypeOf(surface)
            && this.dir.isParallelTo(surface.dir)
            && this.containsCurve(surface.baseCurve)
    }

    like(object: any): boolean {
        if (!this.isCoplanarTo(object)) return false
        // normals need to point in the same direction (outwards or inwards) for both
        const p00 = this.pSTFunc()(0, 0)
        const thisNormal = this.normalSTFunc()(0, 0)
        const otherNormal = object.normalP(p00)
        return 0 < thisNormal.dot(otherNormal)
    }

    loopContainsPoint(loop: Edge[], p: V3): PointVsFace {
        assertVectors(p)
        assert(isFinite(p.x), p.y, p.z)
        const line = new L3(p, this.dir.unit())
        const ptpf = this.stPFunc()
        const pp = ptpf(p)
        if (isNaN(pp.x)) {
            console.log(this.sce, p.sce)
            assert(false)
        }
        const lineOut = this.baseCurve.tangentAt(pp.x).rejectedFrom(this.dir)

        return Surface.loopContainsPointGeneral(loop, p, line, lineOut)
    }


	edgeLoopCCW(loop: Edge[]): boolean {
        if (loop.length < 56) {
            let totalAngle = 0
            for (let i = 0; i < loop.length; i++) {
                const ipp = (i + 1) % loop.length
                const edge = loop[i], nextEdge = loop[ipp]
                totalAngle += edge.bDir.angleRelativeNormal(nextEdge.aDir, this.normalP(edge.b))
            }
            return totalAngle > 0
        } else {
            const ptpF = this.stPFunc()
            return isCCW(loop.map(e => ptpF(e.a)), V3.Z)
        }
    }

    transform<T extends ProjectedCurveSurface>(this: T, m4: M4): T {
    	const f = m4.isMirroring() ? -1 : 1
        return new this.constructor<T>(
            this.baseCurve.transform(m4),
	        m4.transformVector(this.dir).times(f),
            this.sMin, this.sMax, 1 == f ? this.tMin : -this.tMax, 1 == f ? this.tMax : -this.tMin)
    }

    isTsForLine(line: L3): number[] {
        assertInst(L3, line)
        const projPlane = new P3(this.dir.unit(), 0)
        const projDir = projPlane.projectedVector(line.dir1)
        if (projDir.likeO()) {
            // line is parallel to this.dir
            return []
        }
        const projAnchor = projPlane.projectedPoint(line.anchor)
        const projBaseCurve = this.baseCurve.project(projPlane)
        return projBaseCurve
            .isInfosWithLine(projAnchor, projDir, this.sMin, this.sMax, line.tMin, line.tMax)
            .map(info => info.tOther)
    }

    flipped<T extends ProjectedCurveSurface>(this: T): T {
        return new this.constructor<T>(this.baseCurve, this.dir.negated(), this.sMin, this.sMax, -this.tMax, -this.tMin)
    }
}
ProjectedCurveSurface.prototype.uStep = 1 / 40
ProjectedCurveSurface.prototype.vStep = 256