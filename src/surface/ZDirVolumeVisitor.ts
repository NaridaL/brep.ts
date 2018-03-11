/**
 * @prettier
 */
import {
	V3,
	assertNever,
	TAU,
	NLA_PRECISION,
	M4,
	gaussLegendreQuadrature24,
	glqInSteps,
	V,
	assert,
	eq,
	gaussLegendre24Xs,
	gaussLegendre24Weights,
} from 'ts3dutils'
import {
	ConicSurface,
	Edge,
	PlaneSurface,
	SemiCylinderSurface,
	SemiEllipseCurve,
	SemiEllipsoidSurface,
	HyperbolaCurve,
	ParabolaCurve,
	CylinderSurface,
	EllipseCurve,
	EllipsoidSurface,
	L3,
	ProjectedCurveSurface,
	ImplicitCurve,
	P3,
} from '../index'

import { PI, sin, cos, sqrt } from '../math'

export const ZDirVolumeVisitor: { [className: string]: (edges: Edge[]) => { volume: number; centroid: any } } = {
	/**
	 * at(t)
	 * |\                                    ^
	 * | \ at(t).projectedOn(dir)             \  dir
	 * |  \                                    \
	 * |   \ at(t).rejectedFrom(dir)
	 * |   |
	 * |___|
	 *        z = 0
	 *
	 *
	 * A = ((at(t) + at(t).rejectedFrom(dir)) / 2).z * at(t).projectedOn(dir).lengthXY()
	 * scaling = tangentAt(t) DOT dir.cross(V3.Z).unit()
	 */
	[ConicSurface.name](this: ConicSurface, edges: Edge[]): { volume: number; centroid: any } {
		// INT[edge.at; edge.bT] (at(t) DOT dir) * (at(t) - at(t).projectedOn(dir) / 2).z dt
		const totalVolume = edges
			.map(edge => {
				if (
					edge.curve instanceof SemiEllipseCurve ||
					edge.curve instanceof HyperbolaCurve ||
					edge.curve instanceof ParabolaCurve
				) {
					const f = t => {
						const at = edge.curve.at(t),
							tangent = edge.tangentAt(t)
						return (
							(at.z + at.rejectedFrom(this.dir).z) /
							2 *
							at.projectedOn(this.dir).lengthXY() *
							tangent.dot(V3.Z.cross(this.dir).unit())
						)
					}
					// ellipse with normal1 parallel to dir need to be counted negatively so CCW faces result in a positive
					// area
					const sign =
						edge.curve instanceof SemiEllipseCurve
							? -Math.sign(edge.curve.normal.dot(this.dir))
							: -Math.sign(
									this.center
										.to(edge.curve.center)
										.cross(edge.curve.f1)
										.dot(this.dir),
							  )
					const val = glqInSteps(f, edge.aT, edge.bT, 1)
					return val * sign
				} else if (edge.curve instanceof L3) {
					return 0
				} else {
					assertNever()
				}
			})
			.sum()

		return { volume: totalVolume * Math.sign(this.normal.dot(this.dir)) }
	},

	[PlaneSurface.name](this: PlaneSurface, edges: Edge[]): { centroid: V3; volume: number } {
		const { centroid, area } = this.calculateArea()
		return {
			volume: this.surface.plane.normal1.z * centroid.z * area,
			centroid: new V3(centroid.x, centroid.y, centroid.z / 2),
		}
	},
	/**
	 * at(t)
	 * |\                                    ^
	 * | \ at(t).projectedOn(dir1)            \  dir1
	 * |  \                                    \
	 * |   \ at(t).rejectedFrom(dir1)
	 * |   |
	 * |___|
	 *        z = 0
	 *
	 *
	 * A = ((at(t) + at(t).rejectedFrom(dir1)) / 2).z * at(t).projectedOn(dir1).lengthXY()
	 * scaling = tangentAt(t) DOT dir1.cross(V3.Z).unit()
	 */
	[SemiCylinderSurface.name](this: SemiCylinderSurface, edges: Edge[]): { volume: number; centroid: any } {
		if (V3.Z.cross(this.dir).likeO()) return { volume: 0 }
		// the tangent needs to be projected onto a vector which is perpendicular to the volume-slices
		const scalingVector = this.dir.cross(V3.Z).unit()
		// the length of the base of the trapezoid is calculated by dotting with the baseVector
		const baseVector = this.dir.rejectedFrom(V3.Z).unit()
		// INT[edge.at; edge.bT] (at(t) DOT dir1) * (at(t) - at(t).projectedOn(dir) / 2).z
		const totalArea = edges
			.map(edge => {
				if (edge.curve instanceof SemiEllipseCurve) {
					const f = t => {
						// use curve.tangent not edge.tangent, reverse edges are handled by the integration boundaries
						const at = edge.curve.at(t),
							tangent = edge.curve.tangentAt(t)
						const area = (at.z + at.rejectedFrom(this.dir).z) / 2 * at.projectedOn(this.dir).dot(baseVector)
						const scale = tangent.dot(scalingVector)
						//assert(Math.sign(scale) == Math.sign(this.normalP(at).dot(V3.Z)), this.normalP(at).dot(V3.Z))
						//console.log(
						//	"", t,
						//	",", area,
						//	",", scale,
						//	"atz", at.z)
						return area * scale
					}
					// ellipse with normal parallel to dir1 need to be counted negatively so CCW faces result in a positive
					// area
					const sign = -Math.sign(edge.curve.normal.dot(this.dir))
					const val = glqInSteps(f, edge.aT, edge.bT, 1)
					return val * sign
				} else if (edge.curve instanceof L3) {
					return 0
				} else {
					assertNever()
				}
			})
			.sum()

		return { volume: totalArea * Math.sign(this.baseCurve.normal.dot(this.dir)) }
	},

	/**
	 * at(t)
	 * |\                                    ^
	 * | \ at(t).projectedOn(dir1)            \  dir1
	 * |  \                                    \
	 * |   \ at(t).rejectedFrom(dir1)
	 * |   |
	 * |___|
	 *        z = 0
	 *
	 *
	 * A = ((at(t) + at(t).rejectedFrom(dir1)) / 2).z * at(t).projectedOn(dir1).lengthXY()
	 * scaling = tangentAt(t) DOT dir1.cross(V3.Z).unit()
	 */
	[ProjectedCurveSurface.name](this: ProjectedCurveSurface, edges: Edge[]) {
		if (V3.Z.cross(this.dir).likeO()) return { volume: 0, centroid: V3.O }
		// NB: the surface normal doesn't influence the resulting area:
		// normalize this.dir so it always points up
		const upDir1 = this.dir.toLength(Math.sign(this.dir.z) || 1)
		const scalingVector = V3.Z.cross(upDir1).unit()
		// the length of the base of the trapezoid is calculated by dotting with the baseVector
		const baseVector = upDir1.rejectedFrom(V3.Z).unit()
		// INT[edge.at; edge.bT] (at(t) DOT dir1) * (at(t) - at(t).projectedOn(dir) / 2).z
		const totalVolume = edges
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
							tangent = tangents[i]
						const area = (at.z + at.rejectedFrom1(upDir1).z) / 2 * at.projectedOn(this.dir).dot(baseVector)
						const scale = tangent.dot(scalingVector)
						sum += area * scale
					}
					const f = (t: number) => {
						const at = edge.curve.at(t),
							tangent = edge.curve.tangentAt(t)
						const area = (at.z + at.rejectedFrom1(upDir1).z) / 2 * at.projectedOn(this.dir).dot(baseVector)
						const scale = tangent.dot(scalingVector)
						return area * scale
					}
					sum += f(minT) * -(start - minT)
					sum += f(maxT) * -(maxT - end)
					return sum * Math.sign(edge.deltaT())
				} else {
					const f = (t: number) => {
						// use curve.tangent not edge.tangent, reverse edges are handled by the integration boundaries
						const at = edge.curve.at(t),
							tangent = edge.curve.tangentAt(t)
						const area = (at.z + at.rejectedFrom1(upDir1).z) / 2 * at.projectedOn(this.dir).dot(baseVector)
						const scale = tangent.dot(scalingVector)
						return area * scale
					}
					const val = glqInSteps(f, edge.aT, edge.bT, 1)
					return val
				}
			})
			.sum()

		return { volume: totalVolume }
	},

	zDirVolumeForLoop2(loop: Edge[]): number {
		const angles = this.inverseMatrix.getZ().toAngles()
		const T = M4.rotateY(-angles.theta)
			.times(M4.rotateZ(-angles.phi))
			.times(this.inverseMatrix)
		const rot90x = M4.rotateX(PI / 2)
		let totalVolume = 0
		assert(V3.X.isParallelTo(T.transformVector(V3.Z)))
		//const zDistanceFactor = toT.transformVector(V3.Z).length()
		loop.map(edge => edge.transform(T)).forEach((edge, edgeIndex, edges) => {
			const nextEdgeIndex = (edgeIndex + 1) % edges.length,
				nextEdge = edges[nextEdgeIndex]

			function f(t) {
				const at2d = edge.curve.at(t).withElement('x', 0)
				const result =
					1 /
					3 *
					(1 - (at2d.y ** 2 + at2d.z ** 2)) *
					edge.tangentAt(t).dot(rot90x.transformVector(at2d.unit()))
				console.log('at2d', at2d.sce, 'result', result)
				return result
			}

			//if (edge.)
			if (edge.b.like(V3.X)) {
				const angleDiff = (edge.bDir.angleRelativeNormal(nextEdge.aDir, V3.X) + 2 * PI) % (2 * PI)
				totalVolume += 2 / 3 * angleDiff
				console.log('xaa')
			}
			if (edge.b.like(V3.X.negated())) {
				const angleDiff = (edge.bDir.angleRelativeNormal(nextEdge.aDir, V3.X) + 2 * PI) % (2 * PI)
				totalVolume += 2 / 3 * angleDiff
				console.log('xbb')
			}
			const volume = gaussLegendreQuadrature24(f, edge.aT, edge.bT)
			console.log('edge', edge, 'volume', volume)
			totalVolume += volume
		})

		return totalVolume * this.f1.dot(this.f2.cross(this.f3))
	},

	// volume does scale linearly, so this could be done in the local coordinate system
	// however, shear matrices lead to point-to-plane distances having to calculated along a vector other than
	// the plane normal
	[SemiEllipsoidSurface.name](this: SemiEllipsoidSurface, loop: Edge[]): { volume: number; centroid: V3 } {
		const dpds = this.dpds()
		const dpdt = this.dpdt()
		const totalVolume = loop
			.map((edgeWC, edgeIndex, edges) => {
				const curveWC = edgeWC.curve

				const f = (curveT: number) => {
					const pWC = curveWC.at(curveT),
						tangentWC = curveWC.tangentAt(curveT)
					const stOfPWC = this.stP(pWC)
					// const slice = (phi: number) => {
					// 	const p = this.pST(phi, stOfPWC.y)
					// 	const normal = dpds(phi, stOfPWC.y).cross(dpdt(phi, stOfPWC.y))
					// 	return p.z * normal.z
					// }
					// const pz =
					// 	this.f1.z * cos(t) * cos(s) +
					// 	this.f2.z * cos(t) * sin(s) +
					// 	this.f3.z * sin(t) +
					// 	this.center.z
					// const dpdsy = this.f1.y * cos(t) * -sin(s) + this.f2.y * cos(t) * cos(s)
					// const dpdsx = this.f1.x * cos(t) * -sin(s) + this.f2.x * cos(t) * cos(s)
					// const dpdty = this.f1.y * -sin(t) * cos(s) + this.f2.y * -sin(t) * sin(s) + this.f3.y * cos(t)
					// const dpdtx = this.f1.x * -sin(t) * cos(s) + this.f2.x * -sin(t) * sin(s) + this.f3.x * cos(t)
					// const normalz = dpdsx * dpdty - dpdsy * dpdtx
					// result = pz * normalz
					const c = cos(stOfPWC.y),
						s = sin(stOfPWC.y)
					const { x: x_1, y: y_1, z: z_1 } = this.f1
					const { x: x_2, y: y_2, z: z_2 } = this.f2
					const { x: x_3, y: y_3, z: z_3 } = this.f3
					const { z: z_4 } = this.center
					// the following was generated by expanding the above and integrating with Wolfram-Alpha, hence the
					// ugliness
					const sliceIntegral = (p: number) => {
						return (
							0.25 *
							c *
							(x_1 *
								(4 * s * y_2 * (c * sin(p) * z_1 - c * cos(p) * z_2 + p * (s * z_3 + z_4)) +
									c *
										y_3 *
										(c * cos(2 * p) * z_1 +
											c * (-2 * p + sin(2 * p)) * z_2 +
											4 * cos(p) * (s * z_3 + z_4))) -
								c *
									x_3 *
									(y_1 *
										(c * cos(2 * p) * z_1 +
											c * (-2 * p + sin(2 * p)) * z_2 +
											4 * cos(p) * (s * z_3 + z_4)) +
										y_2 *
											(c * (2 * p + sin(2 * p)) * z_1 -
												c * cos(2 * p) * z_2 +
												4 * sin(p) * (s * z_3 + z_4))) +
								x_2 *
									(-4 * s * y_1 * (c * sin(p) * z_1 - c * cos(p) * z_2 + p * (s * z_3 + z_4)) +
										c *
											y_3 *
											(c * (2 * p + sin(2 * p)) * z_1 -
												c * cos(2 * p) * z_2 +
												4 * sin(p) * (s * z_3 + z_4))))
						)
					}
					const dt = M4.forSys(dpds(stOfPWC.x, stOfPWC.y), dpdt(stOfPWC.x, stOfPWC.y))
						.inversed()
						.transformVector(tangentWC).y
					const sliceIntegral0ToPWCS = sliceIntegral(stOfPWC.x) // - sliceIntegral(0) (always 0)
					const result = sliceIntegral0ToPWCS * dt
					return result
				}

				return gaussLegendreQuadrature24(f, edgeWC.aT, edgeWC.bT)
			})
			.sum()

		// calc centroid:
		const centroidZX2Parts = loop.map(edgeWC => {
			const f = (curveT: number) => {
				const curveWC = edgeWC.curve
				const pWC = curveWC.at(curveT),
					tangentWC = curveWC.tangentAt(curveT)
				const stOfPWC = this.stP(pWC)
				const slice = (phi: number) => {
					const p = this.pST(phi, stOfPWC.y)
					const normal = dpds(phi, stOfPWC.y).cross(dpdt(phi, stOfPWC.y))
					return p.times(p.z * normal.z)
				}
				const sliceIntegral0ToPWCS = glqV3(slice, 0, stOfPWC.x)
				const dt = M4.forSys(dpds(stOfPWC.x, stOfPWC.y), dpdt(stOfPWC.x, stOfPWC.y))
					.inversed()
					.transformVector(tangentWC).y
				const result = sliceIntegral0ToPWCS.times(dt)
				return result
			}

			return glqV3(f, edgeWC.aT, edgeWC.bT)
		})
		const centroid = V3.add(...centroidZX2Parts)
			.schur(new V3(1, 1, 0.5))
			.div(totalVolume)
		return { volume: totalVolume, centroid: centroid }
		const localBasePlane = P3.XY.transform(this.inverseMatrix)
		console.log(this)
		console.log(localBasePlane)
		//const zDistanceFactor = toT.transformVector(V3.Z).length()
		const localVolume = loop
			.map((edgeWC, edgeIndex, edges) => {
				const edgeLC = edgeWC.transform(this.inverseMatrix)
				// const nextEdgeIndex = (edgeIndex + 1) % edges.length, nextEdge = edges[nextEdgeIndex]

				// at(t) = edge.curve.at(t)
				// INT[edge.aT; edge.bT] (
				// 	INT[0; at(t).angleXY()] (
				// 		(V3.polar(phi, at(t).lengthXY()) + ()) dot (surface.normalAt(edge.curve.at(t)))
				// 	) dphi
				// ) dt

				// a horizontal slice of the surface at z:
				// Slice(phi) => {
				// 	const p = (r * cos phi, r * sin phi, z)hi, radiusAtZ) + V(0, 0, z) // also the surface normal in
				// LC return = (localBasePlane.normal1.dot(p) - localBasePlane.w) * localBasePlane.normal1.dot(p)
				// return = localBasePlane.normal1.dot(p) ** 2 - localBasePlane.w * localBasePlane.normal1.dot(p) }
				// integral(((r cos(ϕ)) n_1 + r sin(ϕ) n_2 + z n_3)^2 - w ((r cos(ϕ)) n_1 + r sin(ϕ) n_2 + z n_3)) dϕ =
				// 1/4 (n_1^2 r^2 (2 ϕ + sin(2 ϕ)) + 2 n_2^2 r^2 (ϕ - sin(ϕ) cos(ϕ)) + 4 n_2 r cos(ϕ) (w - 2 n_3 z) - 2
				// n_1 r (n_2 r cos(2 ϕ) + 2 sin(ϕ) (w - 2 n_3 z)) + 4 n_3 z ϕ (n_3 z - w)) + constant
				function f(t: number) {
					const atLC = edgeLC.curve.at(t),
						tangent = edgeLC.curve.tangentAt(t)
					const atLCST = SemiEllipsoidSurface.UNIT.stP(atLC)
					const { x: nx, y: ny, z: nz } = localBasePlane.normal1,
						w = localBasePlane.w
					const r = atLC.lengthXY(),
						z = atLC.z
					function sliceIntegral(phi: number) {
						return (
							1 /
							4 *
							(nx ** 2 * r ** 2 * (2 * phi + sin(2 * phi)) +
								2 * ny ** 2 * r ** 2 * (phi - sin(phi) * cos(phi)) +
								4 * ny * r * cos(phi) * (w - 2 * nz * z) -
								2 * nx * r * (ny * r * cos(2 * phi) + 2 * sin(phi) * (w - 2 * nz * z)) +
								4 * nz * z * phi * (nz * z - w))
						)
					}
					// slice = normal(s, t).length() * (plane.normal.dot(at(s, t)) - w)
					function slice(phi: number) {
						const p = V(r * cos(phi), r * sin(phi), z)
						const dpdt = V(r * -sin(phi), r * cos(phi), 0)
						// return (localBasePlane.normal1.dot(p) - localBasePlane.w) *localBasePlane.normal1.dot(p) * dpdt.cross(localBasePlane.normal1).length()
						return (
							SemiEllipsoidSurface.UNIT.dpds()(phi, atLCST.y)
								.cross(SemiEllipsoidSurface.UNIT.dpdt()(phi, atLCST.y))
								.length() *
							(localBasePlane.normal1.dot(p) - localBasePlane.w) * // height
							localBasePlane.normal1.dot(p)
						)
					}
					const testValue = gaussLegendreQuadrature24(slice, 0, atLCST.x)
					// scaling = sqrt(1 - tangent.z ** 2)
					const scaling = SemiEllipsoidSurface.UNIT.dpdt()(atLCST.x, atLCST.y).dot(tangent) // atLC.getPerpendicular().cross(atLC).unit().dot(tangent)
					console.log('scaling', scaling)
					const result = testValue * scaling
					// edge.tangentAt(t).dot(at2d.unit())', edge.tangentAt(t).dot(at2d.unit()))
					return result
				}

				return gaussLegendreQuadrature24(f, edgeLC.aT, edgeLC.bT)
			})
			.sum()
		const volumeScalingFactor = this.f1.dot(this.f2.cross(this.f3))
		console.log('localVolume', localVolume, 'volumeScalingFactor', volumeScalingFactor)
		return { volume: localVolume * volumeScalingFactor, centroid: undefined }
	},
}

function glqV3(f: (x: number) => V3, startT: number, endT: number) {
	return gaussLegendre24Xs
		.reduce((val, currVal, index) => {
			const x = startT + (currVal + 1) / 2 * (endT - startT)
			return val.plus(f(x).times(gaussLegendre24Weights[index]))
		}, V3.O)
		.times((endT - startT) / 2)
}
