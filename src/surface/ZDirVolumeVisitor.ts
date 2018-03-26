import {
	assert,
	eq,
	gaussLegendre24Weights,
	gaussLegendre24Xs,
	gaussLegendreQuadrature24,
	glqInSteps,
	M4,
	NLA_PRECISION,
	V,
	V3,
} from 'ts3dutils'
import {
	ConicSurface,
	Edge,
	HyperbolaCurve,
	ImplicitCurve,
	L3,
	P3,
	ParabolaCurve,
	ParametricSurface,
	PlaneSurface,
	planeSurfaceAreaAndCentroid,
	ProjectedCurveSurface,
	RotatedCurveSurface,
	SemiEllipseCurve,
	SemiEllipsoidSurface,
} from '../index'

import { cos, sin } from '../math'

/**
 * In general: the z-dir shadow volume of a face is the integral: SURFACE_INTEGRAL[p in face] (normal(p).z * p.z) dp
 * In general: the centroid of the z-dir shadow volume of a face is the integral:
 *     SURFACE_INTEGRAL[p in face] ((p schur (1, 1, 0.5)) * normal(p).z * p.z) dp
 *     dividing the z component by 2 is usually done at the very end
 */

export const ZDirVolumeVisitor: { [className: string]: (edges: Edge[]) => { volume: number; centroid: V3 } } = {
	[ConicSurface.name](this: ConicSurface, edges: Edge[]): { volume: number; centroid: V3 } {
		console.log(this)
		const dpds = this.dpds()
		const dpdt = this.dpdt()
		// INT[edge.at; edge.bT] (at(t) DOT dir) * (at(t) - at(t).projectedOn(dir) / 2).z dt
		const totalVolume = edges
			.map(edgeWC => {
				const curveWC = edgeWC.curve
				if (
					curveWC instanceof SemiEllipseCurve ||
					curveWC instanceof HyperbolaCurve ||
					curveWC instanceof ParabolaCurve
				) {
					const f = (curveT: number) => {
						const at = curveWC.at(curveT),
							tangentWC = curveWC.tangentAt(curveT)
						const stOfPWC = this.stP(at)
						// INTEGRATE [0; atST.y] (dpds(atST.x, t) X dpdt(atST.x)).z * pST(atST.x, t).z dt
						// dpds(s, t) === t * dpds(s, 1)
						// => INTEGRATE [0; atST.y] (t * dpds(atST.x, 1) X dpdt(atST.x)).z * pST(atST.x, t).z dt
						// => (dpds(atST.x, 1) X dpdt(atST.x)).z * INTEGRATE [0; atST.y] t * pST(atST.x, t).z dt
						// pST(s, t) === t * (pST(s, 1) - center) + center
						// => (dpds(atST.x, 1) X dpdt(atST.x)).z
						//      * INTEGRATE [0; atST.y] t² * (pST(atST.x, t) - center).z + t * center.z dt
						// => (dpds(atST.x, 1) X dpdt(atST.x)).z
						//      * INTEGRATE [0; atST.y] t² * (pST(atST.x, t) - center).z + t * center.z dt
						// => (dpds(atST.x, 1) X dpdt(atST.x)).z
						//      * (1/3 t³ pST(atST.x, 1).z + 1/2 t² center.z)[0; atST.y]

						const ds = -M4.forSys(dpds(stOfPWC.x, stOfPWC.y), dpdt(stOfPWC.x))
							.inversed()
							.transformVector(tangentWC).x
						const factor =
							stOfPWC.y ** 3 / 3 * (this.pST(stOfPWC.x, 1).z - this.center.z) +
							stOfPWC.y ** 2 / 2 * this.center.z
						const actual = dpds(stOfPWC.x, factor).cross(dpdt(stOfPWC.x)).z
						return actual * ds
					}
					const val = glqInSteps(f, edgeWC.aT, edgeWC.bT, 1)
					return val
				} else if (curveWC instanceof L3) {
					return 0
				} else {
					throw new Error()
				}
			})
			.sum()
		const centroidZX2Parts = edges.map(edgeWC => {
			const curveWC = edgeWC.curve
			if (
				curveWC instanceof SemiEllipseCurve ||
				curveWC instanceof HyperbolaCurve ||
				curveWC instanceof ParabolaCurve
			) {
				const f = (curveT: number) => {
					const at = curveWC.at(curveT),
						tangentWC = curveWC.tangentAt(curveT)
					const stOfPWC = this.stP(at)
					// INTEGRATE [0; atST.y] dpds(atST.x, t) X dpdt(atST.x, t) * pST(atST.x, t).z dt
					// dpdt is constant with respect to t
					// => (dpds(atST.x, t) X dpdt(atST.x, t)).z
					//      * (INTEGRATE [0; atST.y] t * pST(atST.x, t) * pST(atST.x, t).z dt)
					// dpds(s, t) === t * dpds(s, 1)
					// pST(s, t) === t * (pST(s, 1) - center) + center
					// INTEGRATE [0; atST.y] t * pST(atST.x, t) * pST(atST.x, t).z dt
					// = INTEGRATE [0; atST.y] t *
					//                         (t * (pST(s, 1) - center) + center) *
					//                         (t (pST(s, 1) - center).z + center.z) dt
					// = INTEGRATE [0; atST.y] t³ (pST(s, 1) - center) * (pST(s, 1) - center).z
					//                       + t² ((pST(s, 1) - center) * center.z + (pST(s, 1) - center).z * center)
					//                       + t center center.z dt
					// = (1/4 t^4 (pST(s, 1) - center) * (pST(s, 1) - center).z
					//   (1/3 t³ ((pST(s, 1) - center) * center.z + (pST(s, 1) - center).z * center)
					//   (1/2 t² center center.z dt)[0; atST.y]
					const pSTS1V = this.pST(stOfPWC.x, 1).minus(this.center)
					const factor = V3.add(
						pSTS1V.times(1 / 4 * stOfPWC.y ** 4 * pSTS1V.z + 1 / 3 * stOfPWC.y ** 3 * this.center.z),
						this.center.times(1 / 3 * stOfPWC.y ** 3 * pSTS1V.z + 1 / 2 * stOfPWC.y ** 2 * this.center.z),
					)
					const partialCentroid = factor.times(dpds(stOfPWC.x, 1).cross(dpdt(stOfPWC.x)).z)

					const ds = -M4.forSys(dpds(stOfPWC.x, stOfPWC.y), dpdt(stOfPWC.x))
						.inversed()
						.transformVector(tangentWC).x

					return partialCentroid.times(ds)
				}
				return glqV3(f, edgeWC.aT, edgeWC.bT)
			} else if (curveWC instanceof L3) {
				return V3.O
			} else {
				throw new Error()
			}
		})

		const centroid = V3.add(...centroidZX2Parts)
			.schur(new V3(1, 1, 0.5))
			.div(totalVolume)
		return { volume: totalVolume, centroid: centroid }
	},

	[PlaneSurface.name](this: PlaneSurface, edges: Edge[]): { centroid: V3; volume: number } {
		const { centroid, area } = planeSurfaceAreaAndCentroid(this, edges)
		return {
			volume: this.plane.normal1.z * centroid.z * area,
			centroid: new V3(centroid.x, centroid.y, centroid.z / 2),
		}
	},

	/**
	 * Generic implementation.
	 */
	[ParametricSurface.name](this: ParametricSurface, edges: Edge[]): { centroid: V3; volume: number } {
		const dpds = this.dpds()
		const dpdt = this.dpdt()
		const volume = edges
			.map(edgeWC => {
				const curveWC = edgeWC.curve
				if (curveWC instanceof ImplicitCurve) {
					throw new Error()
				} else {
					const f = (curveT: number) => {
						// use curve.tangent not edge.tangent, reverse edges are handled by the integration boundaries
						const pWC = curveWC.at(curveT),
							tangentWC = curveWC.tangentAt(curveT)
						const stOfPWC = this.stP(pWC)
						const slice = (t: number) => {
							const p = this.pST(stOfPWC.x, t)
							const normal = dpds(stOfPWC.x, t).cross(dpdt(stOfPWC.x, t))
							return p.z * normal.z
						}
						const sliceIntegral0ToPWCT = glqInSteps(slice, 0, stOfPWC.y, 1)
						// const dt = tangentWC.dot(scalingVector)
						const dt = -M4.forSys(dpds(stOfPWC.x, stOfPWC.y), dpdt(stOfPWC.x, stOfPWC.y))
							.inversed()
							.transformVector(tangentWC).x
						const result = sliceIntegral0ToPWCT * dt
						return result
					}
					const val = glqInSteps(f, edgeWC.aT, edgeWC.bT, 1)
					return val
				}
			})
			.sum()
		const centroidParts = edges.map(edgeWC => {
			const curveWC = edgeWC.curve
			console.log(edgeWC.sce)
			const f = (curveT: number) => {
				const pWC = curveWC.at(curveT),
					tangentWC = curveWC.tangentAt(curveT)
				const stOfPWC = this.stP(pWC)
				const slice = (t: number) => {
					const p = this.pST(stOfPWC.x, t)
					const normal = dpds(stOfPWC.x, t).cross(dpdt(stOfPWC.x, t))
					return p.times(p.z * normal.z)
				}
				const sliceIntegral0ToPWCT = glqV3(slice, 0, stOfPWC.y)
				// const dt = tangentWC.dot(scalingVector)
				const dt = -M4.forSys(dpds(stOfPWC.x, stOfPWC.y), dpdt(stOfPWC.x, stOfPWC.y))
					.inversed()
					.transformVector(tangentWC).x
				const result = sliceIntegral0ToPWCT.times(dt)
				return result
			}

			return glqV3(f, edgeWC.aT, edgeWC.bT)
		})
		const centroid = V3.add(...centroidParts).div(volume)
		return { volume, centroid }
	},

	/**
	 * at(t)
	 * |\                                    ^
	 * | \ at(t).projectedOn(dir1)            \  dir1
	 * |  \                                    \
	 * |   \ at(t).rejectedFrom(dir1) = b
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
		// normalize this.dir so it always points up
		const upDir1 = this.dir.toLength(Math.sign(this.dir.z) || 1)
		const scalingVector = V3.Z.cross(upDir1).unit()
		// the length of the base of the trapezoid is calculated by dotting with the baseVector
		const baseVector = upDir1.rejectedFrom(V3.Z).unit()
		// INT[edge.at; edge.bT] (at(t) DOT dir1) * (at(t) - at(t).projectedOn(dir) / 2).z
		const volume = edges
			.map(edgeWC => {
				if (edgeWC.curve instanceof L3) {
					return 0
				} else if (edgeWC.curve instanceof ImplicitCurve) {
					const { points, tangents } = edgeWC.curve
					const minT = edgeWC.minT,
						maxT = edgeWC.maxT
					let sum = 0
					const start = Math.ceil(minT + NLA_PRECISION)
					const end = Math.floor(maxT - NLA_PRECISION)
					for (let i = start; i <= end; i++) {
						const at = points[i],
							tangent = tangents[i]
						const area = (at.z + at.rejectedFrom1(upDir1).z) / 2 * at.projectedOn(upDir1).dot(baseVector)
						const scale = tangent.dot(scalingVector)
						sum += area * scale
					}
					const f = (t: number) => {
						const at = edgeWC.curve.at(t),
							tangent = edgeWC.curve.tangentAt(t)
						const area = (at.z + at.rejectedFrom1(upDir1).z) / 2 * at.projectedOn(upDir1).dot(baseVector)
						const scale = tangent.dot(scalingVector)
						return area * scale
					}
					sum += f(minT) * (start - minT - 0.5)
					sum += f(maxT) * (maxT - end - 0.5)
					return sum * Math.sign(edgeWC.deltaT())
				} else {
					const f = (curveT: number) => {
						// use curve.tangent not edge.tangent, reverse edges are handled by the integration boundaries
						const at = edgeWC.curve.at(curveT),
							tangent = edgeWC.curve.tangentAt(curveT)
						const b = at.rejectedFrom1(upDir1)
						const area = at.z * b.to(at).dot(baseVector) / 2 + b.z * b.to(at).dot(baseVector) / 2
						const area2 = (at.z + at.rejectedFrom1(upDir1).z) / 2 * at.projectedOn(upDir1).dot(baseVector)
						assert(eq(area, area2), area, area2)
						const scale = tangent.dot(scalingVector)
						return area * scale
					}
					const val = glqInSteps(f, edgeWC.aT, edgeWC.bT, 4)
					return val
				}
			})
			.sum()
		// calc centroid:
		const centroidParts = edges.map(edgeWC => {
			const fCentroid = (curveT: number) => {
				// use curve.tangent not edge.tangent, reverse edges are handled by the integration boundaries
				const at = edgeWC.curve.at(curveT),
					tangent = edgeWC.curve.tangentAt(curveT)
				const b = at.rejectedFrom1(upDir1)
				const areaCentroidA = V3.add(at.xy(), b, at).times(at.z * b.to(at).dot(baseVector) / 2 / 3)
				const areaCentroidB = V3.add(at.xy(), b, b.xy()).times(b.z * b.to(at).dot(baseVector) / 2 / 3)
				const scale = tangent.dot(scalingVector)

				return areaCentroidA.plus(areaCentroidB).times(scale)
			}

			return glqV3(fCentroid, edgeWC.aT, edgeWC.bT)
		})
		const centroid = V3.add(...centroidParts).div(volume)
		return { volume, centroid }
	},

	// volume does scale linearly, so this could be done in the local coordinate system
	// however, shear matrices lead to point-to-plane distances having to be calculated along a vector other than
	// the plane normal
	[RotatedCurveSurface.name](this: RotatedCurveSurface, edges: Edge[]): { volume: number; centroid: V3 } {
		const dpds = this.dpds()
		const dpdt = this.dpdt()
		const totalVolume = edges
			.map((edgeWC, edgeIndex, edges) => {
				const curveWC = edgeWC.curve

				const f = (curveT: number) => {
					const pWC = curveWC.at(curveT),
						tangentWC = curveWC.tangentAt(curveT)
					const stOfPWC = this.stP(pWC)
					const pLC = this.matrixInverse.transformPoint(pWC)
					const dpdtAtS0 =
						this instanceof RotatedCurveSurface
							? this.curve.tangentAt(stOfPWC.y)
							: V(-pLC.z, 0, pLC.lengthXY())
					// const slice = (phi: number) => {
					// 	const p = this.pST(phi, stOfPWC.y)
					// 	const normal = dpds(phi, stOfPWC.y).cross(dpdt(phi, stOfPWC.y))
					// 	return p.z * normal.z
					// }
					// const z = this.curve.at(stOfPWC.y).z
					// const r = this.curve.at(stOfPWC.y).lengthXY()
					// const pz =
					// 	this.f1.z * r * cos(s) +
					// 	this.f2.z * r * sin(s) +
					// 	this.f3.z * z +
					// 	this.center.z
					// const dpdsx = this.f1.x * r * -sin(s) + this.f2.x * r * cos(s)
					// const dpdsy = this.f1.y * r * -sin(s) + this.f2.y * r * cos(s)
					// const dpdtx = this.f1.x * dr * cos(s) + this.f2.x * dr * sin(s) + this.f3.x * dz
					// const dpdty = this.f1.y * dr * cos(s) + this.f2.y * dr * sin(s) + this.f3.y * dz
					// const normalz = dpdsx * dpdty - dpdsy * dpdtx
					// result = pz * normalz
					const r = pLC.lengthXY(),
						z = pLC.z
					const dr = dpdtAtS0.x
					const dz = dpdtAtS0.z
					const a = this.matrix.X.z * r,
						b = this.matrix.Y.z * r,
						c = this.matrix.Z.z * z + this.matrix.O.z
					const t0 = (this.matrix.X.x * this.matrix.Y.y - this.matrix.X.y * this.matrix.Y.x) * r * dr
					const t1 = (this.matrix.Y.x * this.matrix.X.y - this.matrix.Y.y * this.matrix.X.x) * r * dr
					const t2 = (this.matrix.X.x * this.matrix.X.y - this.matrix.X.y * this.matrix.X.x) * r * dr
					const t3 = (this.matrix.Y.x * this.matrix.Y.y - this.matrix.Y.y * this.matrix.Y.x) * r * dr
					const t4 = (this.matrix.Y.x * this.matrix.Z.y - this.matrix.Y.y * this.matrix.Z.x) * r * dz
					const t5 = (this.matrix.X.x * this.matrix.Z.y - this.matrix.X.y * this.matrix.Z.x) * r * dz
					const sliceIntegral = (p: number) => {
						return (
							(6 * (c * (-t0 + t1) + a * t4 - b * t5) * p +
								3 * (3 * b * t0 - b * t1 + a * (t2 - t3) + 4 * c * t5) * cos(p) +
								3 * (3 * a * t1 - a * t0 - b * (t2 - t3) + 4 * c * t4) * sin(p) +
								3 * (a * t5 - b * t4 + c * (t2 - t3)) * cos(2 * p) +
								3 * (a * t4 + b * t5 + c * (t0 + t1)) * sin(2 * p) +
								(a * (t2 - t3) - b * (t0 + t1)) * cos(3 * p) +
								(a * (t0 + t1) + b * (t2 - t3)) * sin(3 * p)) /
							12
						)
					}
					const dt = M4.forSys(dpds(stOfPWC.x, stOfPWC.y), dpdt(stOfPWC.x, stOfPWC.y))
						.inversed()
						.transformVector(tangentWC).y
					const sliceIntegral0ToPWCS = sliceIntegral(stOfPWC.x) //- sliceIntegral(0) //(always 0)
					const result = sliceIntegral0ToPWCS * dt
					return result
				}

				return gaussLegendreQuadrature24(f, edgeWC.aT, edgeWC.bT)
			})
			.sum()

		// calc centroid:
		const centroidZX2Parts = edges.map(edgeWC => {
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
		const localBasePlane = P3.XY.transform(this.matrixInverse)
		console.log(this)
		console.log(localBasePlane)
		//const zDistanceFactor = toT.transformVector(V3.Z).length()
		const localVolume = edges
			.map((edgeWC, edgeIndex, edges) => {
				const edgeLC = edgeWC.transform(this.matrixInverse)
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
					// slice = normal(s, t).length() * (plane.normal.dot(at(s, t)) - w)
					function slice(phi: number) {
						const p = V(r * cos(phi), r * sin(phi), z)
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
ZDirVolumeVisitor[SemiEllipsoidSurface.name] = ZDirVolumeVisitor[RotatedCurveSurface.name]

export function glqV3(f: (x: number) => V3, startT: number, endT: number) {
	return gaussLegendre24Xs
		.reduce((val, currVal, index) => {
			const x = startT + (currVal + 1) / 2 * (endT - startT)
			return val.plus(f(x).times(gaussLegendre24Weights[index]))
		}, V3.O)
		.times((endT - startT) / 2)
}
