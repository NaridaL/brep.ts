import { assert, assertf, eq, glqInSteps, M4, NLA_PRECISION } from 'ts3dutils'

import {
	ConicSurface,
	Edge,
	EllipseCurve,
	EllipsoidSurface,
	HyperbolaCurve,
	ImplicitCurve,
	L3,
	ParabolaCurve,
	PlaneSurface,
	ProjectedCurveSurface,
	RotatedCurveSurface,
} from '../index'
import { ceil, cos, floor, sign, sin } from '../math'

export const CalculateAreaVisitor = {
	[ConicSurface.name](this: ConicSurface, edges: Edge[]): number {
		const dpdu = this.dpdu()
		const dpdv = this.dpdv()
		// calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
		const totalArea = edges
			.map(edge => {
				if (
					edge.curve instanceof EllipseCurve ||
					edge.curve instanceof HyperbolaCurve ||
					edge.curve instanceof ParabolaCurve
				) {
					const f = (t: number) => {
						const at = edge.curve.at(t),
							tangentWC = edge.tangentAt(t)
						const uvOfPWC = this.uvP(at)
						// INTEGRATE [0; atUV.y]
						//   dpdu(atUV.x, t) X dpdv(atUV.x, t)
						// dt
						// dpdv is constant with respect to t
						// => dpdv(atUV.x, 0) X (INTEGRATE [0; atUV.y] dpdu(atUV.x, t) dt)
						// dpdu(u, v) === v * dpdu(u, 1)
						// => dpdv(atUV.x, 0) X (1/2 t² dpdu(atUV.x, 1))[0; atUV.y]
						// => dpdv(atUV.x, 0) X dpdu(atUV.x, atUV.y² / 2)

						const du = -M4.forSys(dpdu(uvOfPWC.x, uvOfPWC.y), dpdv(uvOfPWC.x))
							.inversed()
							.transformVector(tangentWC).x

						return (
							dpdu(uvOfPWC.x, uvOfPWC.y ** 2 / 2)
								.cross(dpdv(uvOfPWC.x))
								.length() * du
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
		return totalArea * this.normalDir
	},

	[PlaneSurface.name](this: PlaneSurface, edges: Edge[]) {
		let totalArea = 0
		const r1 = this.right,
			u1 = this.up
		for (const edge of edges) {
			let edgeArea: number
			const curve = edge.curve
			if (curve instanceof L3) {
				edgeArea = (edge.a.dot(u1) + edge.b.dot(u1)) / 2 * edge.b.to(edge.a).dot(r1)
			} else if (curve instanceof EllipseCurve) {
				// INTEGRATE[aT; bT] (curve.at(t) * u1) * (tangent(t) * r1) dt
				// INTEGRATE[aT; bT] (u1 f1 cos t + u1 f2 sin t + u1 c) * (r1 f1 (-sin t) + r1 f2 cos t) dt
				const { f1, f2, center } = curve
				const a = u1.dot(f1),
					b = u1.dot(f2),
					c = u1.dot(center),
					d = r1.dot(f1),
					e = r1.dot(f2)
				function fArea(t: number) {
					return (
						0.25 *
						(2 * (-b * d + a * e) * t +
							4 * c * d * cos(t) +
							4 * c * e * sin(t) +
							(a * d - b * e) * cos(2 * t) +
							(b * d + a * e) * sin(2 * t))
					)
				}
				edgeArea = -(fArea(edge.bT) - fArea(edge.aT))
			} else if (curve instanceof ImplicitCurve) {
				throw new Error('implement for implicitCurve')
			} else {
				const dir1 = u1
				assertf(() => dir1.hasLength(1))
				// INT[aT; bT] at(t) * dir1 * tangentAt(t).rejectedFrom(dir1) dt
				const f = (curveT: number) => {
					const at = curve.at(curveT)
					const tangent = curve.tangentAt(curveT)
					const ds = r1.dot(tangent)
					const t = u1.dot(at)
					return ds * t
				}
				edgeArea = glqInSteps(f, edge.aT, edge.bT, 3)
			}

			totalArea += edgeArea
		}
		assert(isFinite(totalArea))
		return totalArea
	},

	[RotatedCurveSurface.name](this: RotatedCurveSurface, edges: Edge[], canApproximate = true): number {
		const f1 = this.matrix.X,
			f2 = this.matrix.Y,
			f3 = this.matrix.Z
		const likeVerticalSpheroid =
			eq(f1.length(), f2.length()) &&
			f1.isPerpendicularTo(f2) &&
			f2.isPerpendicularTo(f3) &&
			f3.isPerpendicularTo(f1)

		const areaParts = edges.map((edgeWC, ei) => {
			console.log('edge', ei, edgeWC.sce)
			const curveWC = edgeWC.curve
			if (edgeWC.curve instanceof ImplicitCurve) {
				throw new Error()
			} else {
				if (likeVerticalSpheroid) {
					const f = (curveT: number) => {
						const pWC = curveWC.at(curveT),
							tangent = curveWC.tangentAt(curveT)
						const pLC = this.matrixInverse.transformPoint(pWC)
						const { x: angleXY, y: t } = this.uvP(pWC)
						const arcRadius = this.matrix.transformVector(pLC.xy()).length()
						const arcLength = angleXY * arcRadius
						const dpdv = this.dpdv()(angleXY, t).unit()
						const scaling = dpdv.dot(tangent)
						return arcLength * scaling
					}
					return glqInSteps(f, edgeWC.aT, edgeWC.bT, 1)
				} else {
					const dpdu = this.dpdu(),
						dpdv = this.dpdv()
					const f2 = (curveT: number) => {
						const pWC = curveWC.at(curveT),
							tangentWC = curveWC.tangentAt(curveT)
						const uvPWC = this.uvP(pWC)
						const slice = (phi: number) => {
							//return this.dpdu()(phi, st.y).length() * this.dpdv()(phi, st.y).length()
							return dpdu(phi, uvPWC.y)
								.cross(dpdv(phi, uvPWC.y))
								.length()
						}
						// we need to do a coordinate transform from curveT to dt, as that is what we are integrating
						const dt = M4.forSys(dpdu(uvPWC.x, uvPWC.y), dpdv(uvPWC.x, uvPWC.y))
							.inversed()
							.transformVector(tangentWC).y
						return glqInSteps(slice, 0, uvPWC.x, 1) * dt
					}
					return glqInSteps(f2, edgeWC.aT, edgeWC.bT, 1)
				}
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
					const start = ceil(minT + NLA_PRECISION)
					const end = floor(maxT - NLA_PRECISION)
					for (let i = start; i <= end; i++) {
						const at = points[i],
							tangent = tangents[i] //.toLength(edge.curve.stepSize)
						const scaling = this.normalP(at)
							.cross(thisDir1)
							.unit()
							.dot(tangent)
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
					sum += f(minT) * (start - minT - 0.5)
					sum += f(maxT) * (maxT - end - 0.5)
					return sum * sign(edge.deltaT())
				} else {
					const f = (t: number) => {
						const at = edge.curve.at(t)
						const tangent = edge.tangentAt(t)
						const scaling = tangent.rejected1Length(thisDir1)
						return at.dot(thisDir1) * scaling
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

	//[CylinderSurface.name](this: CylinderSurface, edges: Edge[]): number {
	//	// calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
	//	const totalArea = edges.map(edge => {
	//		if (edge.curve instanceof EllipseCurve) {
	//			const f = (t: number) => {
	//				const at = edge.curve.at(t), tangent = edge.tangentAt(t)
	//				return at.dot(this.dir) * tangent.rejectedLength(this.dir)
	//			}
	//			// ellipse with normal parallel to dir1 need to be counted negatively so CCW faces result in a
	// positive // area const sign = -Math.sign(edge.curve.normal.dot(this.dir)) const val = glqInSteps(f, edge.aT,
	// edge.bT, 4) return val * sign } else if (edge.curve instanceof L3) { return 0 } else { assertNever() } }).sum()
	// // if the cylinder faces inwards, CCW faces will have been CW, so we need to reverse that here //abs is
	// not an option as "holes" may also be passed return totalArea *sign(this.baseCurve.normal.dot(this.dir)) },
}
CalculateAreaVisitor[EllipsoidSurface.name] = CalculateAreaVisitor[RotatedCurveSurface.name]
