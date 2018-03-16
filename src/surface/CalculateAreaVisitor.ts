/**
 * @prettier
 */
import { assert, glqInSteps, V3, NLA_PRECISION, M4, assertf, eq } from 'ts3dutils'

import {
	BezierCurve,
	ConicSurface,
	Edge,
	EllipseCurve,
	HyperbolaCurve,
	L3,
	ParabolaCurve,
	PlaneSurface,
	ProjectedCurveSurface,
	SemiEllipseCurve,
	StraightEdge,
	ImplicitCurve,
	SemiEllipsoidSurface,
	glqV3,
} from '../index'

export const CalculateAreaVisitor = {
	[ConicSurface.name](this: ConicSurface, edges: Edge[]): number {
		const dpds = this.dpds()
		const dpdt = this.dpdt()
		// calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
		const totalArea = edges
			.map(edge => {
				if (
					edge.curve instanceof SemiEllipseCurve ||
					edge.curve instanceof HyperbolaCurve ||
					edge.curve instanceof ParabolaCurve
				) {
					const f = (t: number) => {
						const at = edge.curve.at(t),
							tangentWC = edge.tangentAt(t)
						const stOfPWC = this.stP(at)
						// INTEGRATE [0; atST.y]
						//   dpds(atST.x, t) X dpdt(atST.x, t)
						// dt
						// dpdt is constant with respect to t
						// => dpdt(atST.x, 0) X (INTEGRATE [0; atST.y] dpds(atST.x, t) dt)
						// dpds(s, t) === t * dpds(s, 1)
						// => dpdt(atST.x, 0) X (1/2 t² dpds(atST.x, 1))[0; atST.y]
						// => dpdt(atST.x, 0) X dpds(atST.x, atST.y² / 2)

						const ds = M4.forSys(dpds(stOfPWC.x, stOfPWC.y), dpdt(stOfPWC.x))
							.inversed()
							.transformVector(tangentWC).x

						return (
							dpds(stOfPWC.x, stOfPWC.y ** 2 / 2)
								.cross(dpdt(stOfPWC.x))
								.length() * ds
						)
					}
					return glqInSteps(f, edge.aT, edge.bT, 1)
				} else if (edge.curve instanceof L3) {
					return 0
				} else {
					throw new Error()
				}
			})
			.sum()
		// if the cylinder faces inwards, CCW faces will have been CW, so we need to reverse that here
		// Math.abs is not an option as "holes" may also be passed
		return totalArea * this.normalDir
	},

	[PlaneSurface.name](this: PlaneSurface, edges: Edge[]) {
		return planeSurfaceAreaAndCentroid(this, edges)
	},

	[SemiEllipsoidSurface.name](this: SemiEllipsoidSurface, edges: Edge[], canApproximate = true): number {
		const areaParts = edges.map((edgeWC, ei) => {
			console.log('edge', ei, edgeWC.sce)
			const curveWC = edgeWC.curve
			if (edgeWC.curve instanceof ImplicitCurve) {
				const { points, tangents } = edgeWC.curve
				const minT = edgeWC.minT,
					maxT = edgeWC.maxT
				let sum = 0
				const start = Math.ceil(minT + NLA_PRECISION)
				const end = Math.floor(maxT - NLA_PRECISION)
				for (let i = start; i <= end; i++) {
					const at = points[i],
						tangent = tangents[i].toLength(edgeWC.curve.stepSize)
					// console.log(
					// 	'at',
					// 	at.sce,
					// 	'tangent',
					// 	tangent.sce,
					// 	'tangent.length()',
					// 	tangent.length(),
					// 	this.normalP(at)
					// 		.cross(thisDir1)
					// 		.unit().sce,
					// )
					const scaling = this.normalP(at)
						.cross(thisDir1)
						.unit()
						.dot(tangent)
					console.log('partsum', at.dot(thisDir1) * scaling)
					sum += at.dot(thisDir1) * scaling
				}
				return sum * Math.sign(edgeWC.deltaT())
			} else if (curveWC instanceof EllipseCurve || curveWC instanceof SemiEllipseCurve) {
				if (this.isVerticalSpheroid()) {
					const circleRadius = this.f1.length()
					const f = (t: number) => {
						const pWC = curveWC.at(t),
							tangent = curveWC.tangentAt(t)
						const pLC = this.inverseMatrix.transformPoint(pWC)
						const angleXY = pLC.angleXY()
						const arcLength = angleXY * circleRadius * Math.sqrt(1 - pLC.z ** 2)
						const dotter = this.matrix
							.transformVector(
								new V3(
									-pLC.z * pLC.x / pLC.lengthXY(),
									-pLC.z * pLC.y / pLC.lengthXY(),
									pLC.lengthXY(),
								),
							)
							.unit()
						const scaling = dotter.dot(tangent)
						return arcLength * scaling
					}
					return glqInSteps(f, edgeWC.aT, edgeWC.bT, 1)
				} else {
					const dpds = this.dpds(),
						dpdt = this.dpdt()
					const f2 = (curveT: number) => {
						const pWC = curveWC.at(curveT),
							tangentWC = curveWC.tangentAt(curveT)
						const stPWC = this.stP(pWC)
						const slice = (phi: number) => {
							//return this.dpds()(phi, st.y).length() * this.dpdt()(phi, st.y).length()
							return dpds(phi, stPWC.y)
								.cross(dpdt(phi, stPWC.y))
								.length()
						}
						// we need to do a coordinate transform from curveT to dt, as that is what we are integrating
						const dt = M4.forSys(dpds(stPWC.x, stPWC.y), dpdt(stPWC.x, stPWC.y))
							.inversed()
							.transformVector(tangentWC).y
						return glqInSteps(slice, 0, stPWC.x, 1) * dt
					}
					return glqInSteps(f2, edgeWC.aT, edgeWC.bT, 1)
				}
			} else {
				throw new Error()
			}
		})
		return areaParts.sum()
	},

	[ProjectedCurveSurface.name](this: ProjectedCurveSurface, edges: Edge[]): number {
		// calculation cannot be done in local coordinate system, as the area doesn't scale proportionally
		const thisDir1 = this.dir.unit()
		const totalArea = edges
			.map(edge => {
				if (edge.curve instanceof L3) {
					return 0
				} else if (edge.curve instanceof ImplicitCurve) {
					const { points, tangents } = edge.curve
					const minT = edge.minT,
						maxT = edge.maxT
					let sum = 0
					const start = Math.ceil(minT + NLA_PRECISION)
					const end = Math.floor(maxT - NLA_PRECISION)
					for (let i = start; i <= end; i++) {
						const at = points[i],
							tangent = tangents[i].toLength(edge.curve.stepSize)
						console.log(
							'at',
							at.sce,
							'tangent',
							tangent.sce,
							'tangent.length()',
							tangent.length(),
							this.normalP(at)
								.cross(thisDir1)
								.unit().sce,
						)
						const scaling = this.normalP(at)
							.cross(thisDir1)
							.unit()
							.dot(tangent)
						console.log('partsum', at.dot(thisDir1) * scaling)
						sum += at.dot(thisDir1) * scaling
					}
					const f = (t: number) => {
						const at = edge.curve.at(t),
							tangent = edge.curve.tangentAt(t)
						const scaling = this.normalP(at)
							.cross(thisDir1)
							.unit()
							.dot(tangent)
						return at.dot(thisDir1) * scaling
					}
					// sum += f(minT) * (start - minT)
					// sum += f(maxT) * (maxT - end)
					console.log('foo', start - minT, maxT - end)
					console.log('start', start, 'end', end, 'minT', minT, 'maxT', maxT)
					sum += f(minT) * (start - minT - 0.5) * 0.5
					sum += f(maxT) * (maxT - end - 0.5) * 0.5
					console.log(sum)
					console.log('f(minT) * (start - minT - 0.5) * 0.5', f(minT) * (start - minT - 0.5) * 0.5)
					console.log('f(maxT) * (maxT - end - 0.5) * 0.5', f(maxT) * (maxT - end - 0.5) * 0.5)
					return sum * Math.sign(edge.deltaT())
				} else {
					const f = (t: number) => {
						const at = edge.curve.at(t),
							tangent = edge.tangentAt(t)
						return at.dot(thisDir1) * tangent.rejected1Length(thisDir1)
					}
					const val = glqInSteps(f, edge.aT, edge.bT, 1)
					const sign = Math.sign(
						this.normalP(edge.a)
							.cross(this.dir)
							.dot(edge.curve.tangentAt(edge.aT)),
					)
					assert(0 !== sign)
					return val * sign
				}
			})
			.sum()
		console.log('totalArea', totalArea)
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

export function planeSurfaceAreaAndCentroid(surface: PlaneSurface, edges: Edge[]) {
	let centroid = V3.O,
		tcs = 0,
		tct = 0,
		totalArea = 0
	const r1 = surface.right,
		u1 = surface.up
	for (const edge of edges) {
		let edgeArea: number, centroidS, centroidT
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
				const parametricCentroid = surface.stPFunc()(info.centroid)
				centroidS = parametricCentroid.x
				centroidT = parametricCentroid.y
			} else if (curve instanceof BezierCurve) {
				const { aT, bT } = edge
				const dir1 = u1
				assertf(() => dir1.hasLength(1))
				// INT[aT; bT] at(t) * dir1 * tangentAt(t).rejectedFrom(dir1) dt
				const f = (t: number) => {
					const tangent = curve.tangentAt(t)
					const at = curve.at(t)
					const outsideVector = tangent.cross(surface.normalP(at))
					const sign = Math.sign(outsideVector.dot(dir1))
					return at.dot(dir1) * tangent.rejected1Length(dir1) * sign
					//return curve.at(t).dot(dir1) * tangent.minus(dir1.times(tangent.dot(dir1))).length()
				}
				const cx = (t: number) => {
					const height = curve.at(t).dot(dir1)
					//console.log(t, curve.at(t).minus(dir1.times(height / 2)).sce, f(t))
					return curve.at(t).minus(dir1.times(height / 2))
				}

				const area = gaussLegendreQuadrature24(f, aT, bT)
				const x = glqV3(cx, aT, bT).div(2 * (bT - aT) * area)
				return { area: area, centroid: x }
				edgeArea = curve.getAreaInDirSurface(u1, surface, edge.aT, edge.bT).area
			} else {
				throw new Error()
			}
		}

		tcs += edgeArea * centroidS
		tct += edgeArea * centroidT
		totalArea += edgeArea
	}
	centroid = r1.times(tcs).plus(u1.times(tct))
	assert(isFinite(totalArea))
	return { area: totalArea, centroid: centroid }
}
