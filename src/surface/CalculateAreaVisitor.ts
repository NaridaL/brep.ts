import {V3,assert,assertNever,eq,glqInSteps} from 'ts3dutils'

import {SemiCylinderSurface, Surface, ProjectedCurveSurface, L3, Edge,
    HyperbolaCurve,
    SemiEllipseCurve,
    ParabolaCurve,
    EllipseCurve, ConicSurface} from '../index'

const {PI} = Math




export const CalculateAreaVisitor: {[className: string]: <T extends Surface>(this: T, allEdges: Edge[]) => number } = {
	[ConicSurface.name](this: ConicSurface, edges: Edge[]): number {
		// calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
		const totalArea = edges.map(edge => {
			if (edge.curve instanceof SemiEllipseCurve || edge.curve instanceof HyperbolaCurve || edge.curve instanceof ParabolaCurve) {
				const f = (t) => {
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
            this.contour.forEach(edge => {
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
                    let curve = edge.curve
                    if (curve instanceof SemiEllipseCurve) {
                        let info = curve.getAreaInDir(r1, u1, edge.aT, edge.bT)
                        edgeArea = info.area
                        let parametricCentroid = this.surface.stPFunc()(info.centroid)
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
            })
            centroid = r1.times(tcs).plus(u1.times(tct))
            assert(isFinite(totalArea))
            return {area: totalArea, centroid: centroid}
    },

	/**
	 * Calculating the surface area of a projected ellipse is analogous to the circumference of the ellipse
	 * ==> Elliptic integrals/numeric calculation is necessary
	 */
	[CylinderSurface.name](this: CylinderSurface, edges: Edge[]): number {
		// calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
		const totalArea = edges.map(edge => {
			if (edge.curve instanceof EllipseCurve) {
				const f = (t) => {
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
		const {f1, f2, f3} = this
		// calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
		const circleRadius = f1.length()
		const f31 = f3.unit()
		const totalArea = edges.map(edge => {
			if (edge.curve instanceof EllipseCurve) {
				const f = (t) => {
					const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t)
					const localAt = this.inverseMatrix.transformPoint(at)
					let angleXY = localAt.angleXY()
					if(eq(Math.abs(angleXY), PI)) {
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

	[SemiCylinderSurface.name](this: SemiCylinderSurface, edges: Edge[], canApproximate = true): number {
		assert(this.isVerticalSpheroid())
		const {f1, f2, f3} = this
		// calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
		const circleRadius = f1.length()
		const f31 = f3.unit()
		const totalArea = edges.map(edge => {
			if (edge.curve instanceof SemiEllipseCurve) {
				const f = (t) => {
					const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t)
					const localAt = this.inverseMatrix.transformPoint(at)
					let angleXY = localAt.angleXY()
					if(eq(Math.abs(angleXY), PI)) {
						if (edge.curve.normal.isParallelTo(this.f2)) {
							angleXY = PI * -Math.sign((edge.bT - edge.aT) * edge.curve.normal.dot(this.f2))
						} else {
							angleXY = PI * dotCurve(this.f2, tangent, edge.curve.ddt(t))
						}
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
				return val
			} else {
				assertNever()
			}
		}).sum()



		return totalArea * Math.sign(this.f1.cross(this.f2).dot(this.f3))
	},

	[ProjectedCurveSurface.name](this: ProjectedCurveSurface, edges: Edge[]): number {
		// calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
		const totalArea = edges.map(edge => {
			if (edge.curve instanceof SemiEllipseCurve) {
				const f = (t) => {
					const at = edge.curve.at(t), tangent = edge.tangentAt(t)
					return at.dot(this.dir) * tangent.rejected1Length(this.dir)
				}
				// ellipse with normal1 parallel to dir1 need to be counted negatively so CCW faces result in a
                // positive area
				const sign = -Math.sign(edge.curve.normal.dot(this.dir))
				const val = glqInSteps(f, edge.aT, edge.bT, 4)
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

	/**
	 * Calculating the surface area of a projected ellipse is analogous to the circumference of the ellipse
	 * ==> Elliptic integrals/numeric calculation is necessary
	 */
	[SemiCylinderSurface.name](this: SemiCylinderSurface, edges: Edge[]): number {
		// calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
		const totalArea = edges.map(edge => {
			if (edge.curve instanceof SemiEllipseCurve) {
				const f = (t) => {
					const at = edge.curve.at(t), tangent = edge.tangentAt(t)
					return at.dot(this.dir) * tangent.rejected1Length(this.dir)
				}
				// ellipse with normal parallel to dir1 need to be counted negatively so CCW faces result in a positive
                // area
				const sign = -Math.sign(edge.curve.normal.dot(this.dir))
				const val = glqInSteps(f, edge.aT, edge.bT, 4)
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
}