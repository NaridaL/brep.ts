import { assert, assertNever, eq, glqInSteps, V3, NLA_PRECISION } from 'ts3dutils'

import {
	BezierCurve, ConicSurface, CylinderSurface, dotCurve, Edge, EllipseCurve, EllipsoidSurface, HyperbolaCurve, L3,
	ParabolaCurve, PlaneSurface, ProjectedCurveSurface, SemiCylinderSurface, SemiEllipseCurve, StraightEdge, PICurve,
    ImplicitCurve,
} from '../index'

import { PI } from '../math'


export const CalculateAreaVisitor = {
	[ConicSurface.name](this: ConicSurface, edges: Edge[]): number {
		// calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
		const totalArea = edges.map(edge => {
			if (edge.curve instanceof SemiEllipseCurve || edge.curve instanceof HyperbolaCurve || edge.curve instanceof ParabolaCurve) {
				const f = (t: number) => {
					const at = edge.curve.at(t), tangent = edge.tangentAt(t)
					return at.minus(this.center).cross(tangent.rejectedFrom(this.dir)).length() / 2
				}
				// ellipse with normal1 parallel to dir1 need to be counted negatively so CCW faces result in a
				// positive area hyperbola normal1 can be perpendicular to
				const sign = edge.curve instanceof SemiEllipseCurve
					? -Math.sign(edge.curve.normal.dot(this.dir))
					: -Math.sign(this.center.to(edge.curve.center).cross(edge.curve.f1).dot(this.dir))
				return glqInSteps(f, edge.aT, edge.bT, 4) * sign
			} else if (edge.curve instanceof L3) {
				return 0
			} else {
				assertNever()
			}
		}).sum()
		// if the cylinder faces inwards, CCW faces will have been CW, so we need to reverse that here
		// Math.abs is not an option as "holes" may also be passed
		return totalArea * Math.sign(this.normal.dot(this.dir))
	},

	[PlaneSurface.name](this: PlaneSurface, edges: Edge[]) {
		let centroid = V3.O, tcs = 0, tct = 0, totalArea = 0
		let r1 = this.surface.right, u1 = this.surface.up
		for (const edge of edges) {
			let edgeCentroid, edgeArea: number, centroidS, centroidT
			if (edge instanceof StraightEdge) {
				const midPoint = edge.a.lerp(edge.b, 0.5)
				edgeCentroid = new V3(midPoint.x, centroid.y, centroid.z / 2)
				centroidS = midPoint.dot(r1) / 2
				centroidT = midPoint.dot(u1)
				const edgeLength = edge.a.distanceTo(edge.b)
				edgeArea = edgeLength * edge.curve.dir1.dot(r1)
				edgeArea = (edge.a.dot(u1) + edge.b.dot(u1)) / 2 * edge.b.to(edge.a).dot(r1)
			} else {
				const curve = edge.curve
				if (curve instanceof SemiEllipseCurve) {
					const info = curve.getAreaInDir(r1, u1, edge.aT, edge.bT)
					edgeArea = info.area
					const parametricCentroid = this.surface.stPFunc()(info.centroid)
					centroidS = parametricCentroid.x
					centroidT = parametricCentroid.y
				} else if (curve instanceof BezierCurve) {
					edgeArea = curve.getAreaInDirSurface(u1, this.surface, edge.aT, edge.bT).area
				} else {
					assertNever()
				}
			}


			tcs += edgeArea * centroidS
			tct += edgeArea * centroidT
			totalArea += edgeArea
		}
		centroid = r1.times(tcs).plus(u1.times(tct))
		assert(isFinite(totalArea))
		return { area: totalArea, centroid: centroid }
	},

	/**
	 * Calculating the surface area of a projected ellipse is analogous to the circumference of the ellipse
	 * ==> Elliptic integrals/numeric calculation is necessary
	 */
	[CylinderSurface.name](this: CylinderSurface, edges: Edge[]): number {
		// calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
		const totalArea = edges.map(edge => {
			if (edge.curve instanceof EllipseCurve) {
				const f = (t: number) => {
					const at = edge.curve.at(t), tangent = edge.tangentAt(t)
					return at.dot(this.dir) * tangent.rejected1Length(this.dir)
				}
				// ellipse with normal parallel to dir1 need to be counted negatively so CCW faces result in a positive
				// area
				const sign = -Math.sign(edge.curve.normal.dot(this.dir))
				const val = glqInSteps(f, edge.aT, edge.bT, 4)
				console.log('edge', edge, val)
				return val * sign
			} else if (edge.curve instanceof L3) {
				return 0
			} else {
				assertNever()
			}
		}).sum()
		// if the cylinder faces inwards, CCW faces will have been CW, so we need to reverse that here
		// Math.abs is not an option as "holes" may also be passed
		return totalArea * Math.sign(this.baseCurve.normal.dot(this.dir))
	},

	[EllipsoidSurface.name](this: EllipsoidSurface, edges: Edge[], canApproximate = true): number {
		assert(this.isVerticalSpheroid())
		const { f1, f2, f3 } = this
		// calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
		const circleRadius = f1.length()
		const f31 = f3.unit()
		const totalArea = edges.map(edge => {
			if (edge.curve instanceof EllipseCurve) {
				const f = (t: number) => {
					const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t)
					const localAt = this.inverseMatrix.transformPoint(at)
					let angleXY = localAt.angleXY()
					if (eq(Math.abs(angleXY), PI)) {
						if (edge.curve.normal.isParallelTo(this.f2)) {
							angleXY = PI * -Math.sign((edge.bT - edge.aT) * edge.curve.normal.dot(this.f2))
						} else {
							angleXY = PI * dotCurve(this.f2, tangent, edge.curve.ddt(t))
						}
						console.log(angleXY)
					}
					const arcLength = angleXY * circleRadius * Math.sqrt(1 - localAt.z ** 2)
					const dotter = this.matrix.transformVector(new V3(-localAt.z * localAt.x / localAt.lengthXY(), -localAt.z * localAt.y / localAt.lengthXY(), localAt.lengthXY())).unit()
					const df3 = tangent.dot(f31)
					//const scaling = df3 / localAt.lengthXY()
					const scaling = dotter.dot(tangent)
					//console.log(t, at.str, arcLength, scaling)
					return arcLength * scaling
				}
				const val = glqInSteps(f, edge.aT, edge.bT, 1)
				console.log('edge', edge, val)
				return val
			} else {
				assertNever()
			}
		}).sum()


		return totalArea * Math.sign(this.f1.cross(this.f2).dot(this.f3))
	},

	//[SemiCylinderSurface.name](this: SemiCylinderSurface, edges: Edge[], canApproximate = true): number {
    //	assert(this.isVerticalSpheroid())
    //	const { f1, f2, f3 } = this
    //	// calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
    //	const circleRadius = f1.length()
    //	const f31 = f3.unit()
    //	const totalArea = edges.map(edge => {
    //		if (edge.curve instanceof SemiEllipseCurve) {
    //			const f = (t: number) => {
    //				const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t)
    //				const localAt = this.inverseMatrix.transformPoint(at)
    //				let angleXY = localAt.angleXY()
    //				if (eq(Math.abs(angleXY), PI)) {
    //					if (edge.curve.normal.isParallelTo(this.f2)) {
    //						angleXY = PI * -Math.sign((edge.bT - edge.aT) * edge.curve.normal.dot(this.f2))
    //					} else {
    //						angleXY = PI * dotCurve(this.f2, tangent, edge.curve.ddt(t))
    //					}
    //				}
    //				const arcLength = angleXY * circleRadius * Math.sqrt(1 - localAt.z ** 2)
    //				const dotter = this.matrix.transformVector(new V3(-localAt.z * localAt.x / localAt.lengthXY(),
    // -localAt.z * localAt.y / localAt.lengthXY(), localAt.lengthXY())).unit() const df3 = tangent.dot(f31) //const
    // scaling = df3 / localAt.lengthXY() const scaling = dotter.dot(tangent) //console.log(t, at.str, arcLength,
    // scaling) return arcLength * scaling } const val = glqInSteps(f, edge.aT, edge.bT, 1) return val } else {
    // assertNever() } }).sum()

    //
	//	return totalArea * Math.sign(this.f1.cross(this.f2).dot(this.f3))
	//},

	[ProjectedCurveSurface.name](this: ProjectedCurveSurface, edges: Edge[]): number {
		// calculation cannot be done in local coordinate system, as the area doesn't scale proportionally
        const thisDir1 = this.dir.unit()
        const totalArea = edges.map(edge => {
            if (edge.curve instanceof L3) {
                return 0
            } else if (edge.curve instanceof ImplicitCurve) {
                const {points, tangents} = edge.curve
                const minT = edge.minT, maxT = edge.maxT
                let sum = 0
                const start = Math.ceil(minT + NLA_PRECISION)
                const end = Math.floor(maxT - NLA_PRECISION)
                for (let i = start; i <= end; i++) {
                    const at = points[i], tangent = tangents[i].toLength(edge.curve.stepSize)
                    console.log("at", at.sce, "tangent", tangent.sce, 'tangent.length()', tangent.length(), this.normalP(at).cross(thisDir1).unit().sce)
                    const scaling = this.normalP(at).cross(thisDir1).unit().dot(tangent)
                    console.log('partsum', at.dot(thisDir1) * scaling)
                    sum += at.dot(thisDir1) * scaling
                }
                const f = (t: number) => {
                    const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t)
                    const scaling = this.normalP(at).cross(thisDir1).unit().dot(tangent)
                    return at.dot(thisDir1) * scaling
                }
                console.log('foo', start - minT, maxT - end)
                console.log('start', start, 'end', end, 'minT', minT, 'maxT', maxT)
                //sum += f(minT) * (start - minT - 0.5) * 0.5
                //sum += f(maxT) * (maxT - end - 0.5) * 0.5
                console.log(sum)
                console.log('f(minT) * (start - minT - 0.5) * 0.5', f(minT) * (start - minT - 0.5) * 0.5)
                console.log('f(maxT) * (maxT - end - 0.5) * 0.5', f(maxT) * (maxT - end - 0.5) * 0.5)
                return sum * Math.sign(edge.deltaT())
            } else {
                const f = (t: number) => {
                    const at = edge.curve.at(t), tangent = edge.tangentAt(t)
                    return at.dot(thisDir1) * tangent.rejected1Length(thisDir1)
                }
                const val = glqInSteps(f, edge.aT, edge.bT, 1)
                const sign = Math.sign(this.normalP(edge.a).cross(this.dir).dot(edge.curve.tangentAt(edge.aT)))
                assert(0 !== sign)
                return val * sign
            }
        }).sum()
        console.log("totalArea", totalArea)
		return totalArea
	},

	//[SemiCylinderSurface.name](this: SemiCylinderSurface, edges: Edge[]): number {
    //	// calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
    //	const totalArea = edges.map(edge => {
    //		if (edge.curve instanceof SemiEllipseCurve) {
    //			const f = (t: number) => {
    //				const at = edge.curve.at(t), tangent = edge.tangentAt(t)
    //				return at.dot(this.dir) * tangent.rejectedLength(this.dir)
    //			}
    //			// ellipse with normal parallel to dir1 need to be counted negatively so CCW faces result in a
    // positive // area const sign = -Math.sign(edge.curve.normal.dot(this.dir)) const val = glqInSteps(f, edge.aT,
    // edge.bT, 4) return val * sign } else if (edge.curve instanceof L3) { return 0 } else { assertNever() } }).sum()
    // // if the cylinder faces inwards, CCW faces will have been CW, so we need to reverse that here // Math.abs is
    // not an option as "holes" may also be passed return totalArea * Math.sign(this.baseCurve.normal.dot(this.dir)) },
}
