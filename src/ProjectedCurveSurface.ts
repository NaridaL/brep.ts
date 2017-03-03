/**
 * Surface normal is (t, z) => this.baseCurve.tangentAt(t) X this.dir1
 * Choose dir1 appropriately to select surface orientation.
 */
class ProjectedCurveSurface extends Surface {
	baseCurve: BezierCurve
	dir1: V3
	tMin: number
	tMax: number


	constructor(baseCurve, dir1, tMin = baseCurve.tMin, tMax = baseCurve.tMax) {
		super()
		assertInst(Curve, baseCurve)
		assertInst(V3, dir1)
		assertNumbers(tMin, tMax)
		assert(dir1.hasLength(1), "dir must be normalized" + dir1)
		this.baseCurve = baseCurve
		this.dir1 = dir1
		this.tMin = tMin
		this.tMax = tMax
	}

	boundsFunction() {
		return (t, z) => {
			return this.tMin <= t && t <= this.tMax
		}
	}

	toSource():string {
		return `new ProjectedCurveSurface(${this.baseCurve}, ${this.dir1}, ${this.tMin}, ${this.tMax})`
	}

	toString():string {
		return this.toSource()
	}


	toMesh(tStart, tEnd, zStart, zEnd, count) {
		tStart = tStart || this.tMin
		tEnd = tEnd || this.tMax
		zStart = zStart || -400
		zEnd = zEnd || 400
		count = count || 128

		let tInterval = tEnd - tStart, tStep = tInterval / (count - 1)
		let baseVertices =
			NLA.arrayFromFunction(count, i => this.baseCurve.at(tStart + i * tStep).plus(this.dir1.times(zStart)))
		let normalFunc = this.parametricNormal()
		let normals = NLA.arrayFromFunction(count, i => normalFunc(tStart + i * tStep, 0))

		return GL.Mesh.offsetVertices(baseVertices, this.dir1.times(zEnd - zStart), false, normals)
	}

	parametricFunction() {
		let baseCurve = this.baseCurve, dir1 = this.dir1
		return (
			/**
			 *
			 * @param t curve Parameter
			 * @param z projection factor
			 * @returns {V3}
			 */
			function (t:number, z:number) {
				return baseCurve.at(t).plus(dir1.times(z))
			})
	}

	parametricNormal() {
		let baseCurve = this.baseCurve, dir1 = this.dir1
		return (t, z) => {
			return baseCurve.tangentAt(t).cross(dir1).normalized()
		}
	}

	implicitFunction() {

	}

	boundsFunction() {
		let ptpf = this.pointToParameterFunction()
		let tMin = this.tMin, tMax = this.tMax
		return function (pWC) {
			let pointParameterT = ptpf(pWC).x
			return tMin <= pointParameterT && pointParameterT <= tMax
		}
	}

	footParameters(pWC, ss, st) {
		let basePlane = new P3(this.dir1, 0)
		let projCurve = this.baseCurve.project(basePlane)
		let projPoint = basePlane.projectedPoint(pWC)
		let t = projCurve.closestTToPoint(projPoint, undefined, undefined, ss)
		let z = pWC.minus(this.baseCurve.at(t)).dot(this.dir1)
		return new V3(t, z, 0)
	}

	pointToParameterFunction() {
		let projPlane = new P3(this.dir1, 0)
		let dir1 = this.dir1, baseCurve = this.baseCurve
		let projBaseCurve = baseCurve.project(projPlane)
		let _this = this
		return function (pWC) {
			let projPoint = projPlane.projectedPoint(pWC)
			let t = projBaseCurve.pointLambda(projPoint)
			let z = pWC.minus(baseCurve.at(t)).dot(dir1)
			return new V3(t, z, 0)
		}

		// let baseCurve = this.baseCurve, dir1 = this.dir1
		// let tMin = this.tMin, tMax = this.tMax
		// return function (pWC) {
		// 	let iterationFunc = t => {let d = baseCurve.at(t).minus(pWC); return d.dot(dir1) + d.length() }
		// 	let startT = NLA.arrayFromFunction(16, i => tMin + (tMax - tMin) * i / 15)
		// 		.withMax(t => -abs(iterationFunc(t)))
		// 	let t = newtonIterate1d(iterationFunc, startT, 16)
		// 	let z = pWC.minus(baseCurve.at(t)).dot(dir1)
		// 	return V(t, z, 0)
		// }
	}
	/**
	 * @inheritDoc
	 */
	isCurvesWithPlane(plane):Curve[] {
		assertInst(P3, plane)
			if (this.dir1.isPerpendicularTo(plane.normal)) {

				let ts = this.baseCurve.isTsWithPlane(plane)
			return ts.map(t => {
				let l3dir = 0 < this.baseCurve.tangentAt(t).dot(plane.normal)
					? this.dir1
					: this.dir1.negated()
				return new L3(this.baseCurve.at(t), l3dir)
			})
			} else {
			let projCurve = this.baseCurve.transform(M4.projection(plane, this.dir1))
			if (this.dir1.dot(plane.normal) > 0) {
				// we need to flip the ellipse so the tangent is correct
				console.log("FLIPPING")
				projCurve = projCurve.reversed()
			}
			return [projCurve]
		}
	}

	/**
	 * @inheritDoc
	 */
	isCurvesWithSurface(surface) {
		if (surface instanceof PlaneSurface) {
			return this.isCurvesWithPlane(surface.plane)
		}
		if (surface instanceof ProjectedCurveSurface || surface instanceof CylinderSurface) {
			let dir1 = surface instanceof ProjectedCurveSurface ? surface.dir1 : surface.dir.normalized()
			if (this.dir1.isParallelTo(dir1)) {
				let otherCurve = surface instanceof ProjectedCurveSurface ? surface.baseCurve : surface.baseEllipse
				let infos = this.baseCurve.isInfosWithCurve(otherCurve)
				return infos.map(info => new L3(info.p, dir1))
			}
			if (surface instanceof ProjectedCurveSurface) {
				const line = new L3(this.baseCurve.at(0.5), this.dir1)
				const startPoint = line.at(surface.isTsForLine(line)[0])
				drVs.push({anchor: this.baseCurve.at(0.5), dir: this.dir1})
				console.log(startPoint)
				return [new PPCurve(this, surface, startPoint)]
				// let testVector = this.dir1.cross(surface.dir1).normalized()
				// // look for points on surface.baseCurve where tangent DOT testVector == 0
				// let abcd1 = surface.baseCurve.tangentCoefficients().map(c => c.dot(testVector))
				// let ts1 = solveCubicReal2.apply(undefined, abcd1).concat(surface.tMin, surface.tMax)
				// let abcd2 = this.baseCurve.tangentCoefficients().map(c => c.dot(testVector))
				// let ts2 = solveCubicReal2.apply(undefined, abcd2)
				// let tt1 = ts1.map(t => surface.baseCurve.at(t).dot(testVector))
				// let tt2 = ts1.map(t => surface.baseCurve.at(t).dot(testVector))
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
		let pp = this.pointToParameterFunction()(p)
		return this.parametricFunction()(pp.x, pp.y).like(p)
	}

	/**
	 * @inheritDoc
	 */
	containsCurve(curve) {
		if (curve instanceof BezierCurve) {
			// project baseCurve and test curve onto a common plane and check if the curves are alike
			let projectionPlane = new P3(this.dir1, 0)
			let baseCurveProjection = this.baseCurve.project(projectionPlane)
			if (curve instanceof L3 && curve.dir1.isParallelTo(this.dir1)) {
				// projection onto basePlane would be single point
				return baseCurveProjection.containsPoint(projectionPlane.projectedPoint(curve.anchor))
			}
			let curveProjection = curve.project(projectionPlane)

			return baseCurveProjection.likeCurve(curveProjection)
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


	/**
	 * @inheritDoc
	 */
	like(object) {
		if (!this.isCoplanarTo(object)) return false
		// normals need to point in the same direction (outwards or inwards) for both
		let p00 = this.parametricFunction()(0, 0)
		let thisNormal = this.parametricNormal()(0, 0)
		let otherNormal = object.normalAt(p00)
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
			for (var i = 0; i < contour.length; i++) {
				var ipp = (i + 1) % contour.length
				var edge = contour[i], nextEdge = contour[ipp]
				totalAngle += edge.bDir.angleRelativeNormal(nextEdge.aDir, this.normalAt(edge.b))
			}
			return totalAngle > 0
		} else {
			var ptpF = this.pointToParameterFunction()
			return isCCW(contour.map(e => ptpF(e.a)), V3.Z)
		}
	}

	/**
	 *
	 * @param m4
	 * @returns {ProjectedCurveSurface}
	 */
	transform(m4:M4) {
		return new ProjectedCurveSurface(this.baseCurve.transform(m4), m4.transformVector(this.dir1).normalized(), this.tMin, this.tMax)
	}

	/**
	 * @inheritDoc
	 */
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
ProjectedCurveSurface.prototype.uStep = 1 / 32
ProjectedCurveSurface.prototype.vStep = 256
NLA.registerClass(ProjectedCurveSurface)