import {
	assert,
	gaussLegendre24Weights,
	gaussLegendre24Xs,
	gaussLegendreQuadrature24,
	glqInSteps,
	M4,
	V,
	V3,
} from 'ts3dutils'
import {
	ConicSurface,
	Edge,
	EllipseCurve,
	EllipsoidSurface,
	HyperbolaCurve,
	ImplicitCurve,
	L3,
	ParabolaCurve,
	ParametricSurface,
	PlaneSurface,
	ProjectedCurveSurface,
	RotatedCurveSurface,
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
		const dpdu = this.dpdu()
		const dpdv = this.dpdv()
		// INT[edge.at; edge.bT] (at(t) DOT dir) * (at(t) - at(t).projectedOn(dir) / 2).z dt
		const totalVolume = edges
			.map(edgeWC => {
				const curveWC = edgeWC.curve
				if (
					curveWC instanceof EllipseCurve ||
					curveWC instanceof HyperbolaCurve ||
					curveWC instanceof ParabolaCurve
				) {
					const f = (curveT: number) => {
						const at = curveWC.at(curveT),
							tangentWC = curveWC.tangentAt(curveT)
						const uvOfPWC = this.uvP(at)
						// INTEGRATE [0; atUV.y] (dpdu(atUV.x, t) X dpdv(atUV.x)).z * pUV(atUV.x, t).z dt
						// dpdu(u, v) === t * dpdu(s, 1)
						// => INTEGRATE [0; atUV.y] (t * dpdu(atUV.x, 1) X dpdv(atUV.x)).z * pUV(atUV.x, t).z dt
						// => (dpdu(atUV.x, 1) X dpdv(atUV.x)).z * INTEGRATE [0; atUV.y] t * pUV(atUV.x, t).z dt
						// pUV(u, v) === t * (pUV(s, 1) - center) + center
						// => (dpdu(atUV.x, 1) X dpdv(atUV.x)).z
						//      * INTEGRATE [0; atUV.y] t² * (pUV(atUV.x, t) - center).z + t * center.z dt
						// => (dpdu(atUV.x, 1) X dpdv(atUV.x)).z
						//      * INTEGRATE [0; atUV.y] t² * (pUV(atUV.x, t) - center).z + t * center.z dt
						// => (dpdu(atUV.x, 1) X dpdv(atUV.x)).z
						//      * (1/3 t³ pUV(atUV.x, 1).z + 1/2 t² center.z)[0; atUV.y]

						const du = -M4.forSys(dpdu(uvOfPWC.x, uvOfPWC.y), dpdv(uvOfPWC.x))
							.inversed()
							.transformVector(tangentWC).x
						const factor =
							uvOfPWC.y ** 3 / 3 * (this.pUV(uvOfPWC.x, 1).z - this.center.z) +
							uvOfPWC.y ** 2 / 2 * this.center.z
						const actual = dpdu(uvOfPWC.x, factor).cross(dpdv(uvOfPWC.x)).z
						return actual * du
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
				curveWC instanceof EllipseCurve ||
				curveWC instanceof HyperbolaCurve ||
				curveWC instanceof ParabolaCurve
			) {
				const f = (curveT: number) => {
					const at = curveWC.at(curveT),
						tangentWC = curveWC.tangentAt(curveT)
					const uvOfPWC = this.uvP(at)
					// INTEGRATE [0; atUV.y] dpdu(atUV.x, t) X dpdv(atUV.x, t) * pUV(atUV.x, t).z dt
					// dpdv is constant with respect to t
					// => (dpdu(atUV.x, t) X dpdv(atUV.x, t)).z
					//      * (INTEGRATE [0; atUV.y] t * pUV(atUV.x, t) * pUV(atUV.x, t).z dt)
					// dpdu(u, v) === t * dpdu(s, 1)
					// pUV(u, v) === t * (pUV(s, 1) - center) + center
					// INTEGRATE [0; atUV.y] t * pUV(atUV.x, t) * pUV(atUV.x, t).z dt
					// = INTEGRATE [0; atUV.y] t *
					//                         (t * (pUV(s, 1) - center) + center) *
					//                         (t (pUV(s, 1) - center).z + center.z) dt
					// = INTEGRATE [0; atUV.y] t³ (pUV(s, 1) - center) * (pUV(s, 1) - center).z
					//                       + t² ((pUV(s, 1) - center) * center.z + (pUV(s, 1) - center).z * center)
					//                       + t center center.z dt
					// = (1/4 t^4 (pUV(s, 1) - center) * (pUV(s, 1) - center).z
					//   (1/3 t³ ((pUV(s, 1) - center) * center.z + (pUV(s, 1) - center).z * center)
					//   (1/2 t² center center.z dt)[0; atUV.y]
					const pUVS1V = this.pUV(uvOfPWC.x, 1).minus(this.center)
					const factor = V3.add(
						pUVS1V.times(1 / 4 * uvOfPWC.y ** 4 * pUVS1V.z + 1 / 3 * uvOfPWC.y ** 3 * this.center.z),
						this.center.times(1 / 3 * uvOfPWC.y ** 3 * pUVS1V.z + 1 / 2 * uvOfPWC.y ** 2 * this.center.z),
					)
					const partialCentroid = factor.times(dpdu(uvOfPWC.x, 1).cross(dpdv(uvOfPWC.x)).z)

					const ds = -M4.forSys(dpdu(uvOfPWC.x, uvOfPWC.y), dpdv(uvOfPWC.x))
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
		const r1 = this.right
		const u1 = this.up
		const c = this.plane.anchor
		assert(r1.hasLength(1))
		assert(u1.hasLength(1))
		assert(r1.isPerpendicularTo(u1))
		const volumeAndCentroidZX2Parts = edges.map((edgeWC): [number, V3] => {
			const curveWC = edgeWC.curve
			if (curveWC instanceof L3) {
				// split shadow volume into two triangle shadow volumes and use the same logic as for mesh triangles:
				function triangleShadowVolumeAndCentroid(a: V3, b: V3, c: V3): [number, V3] {
					const ab = b.minus(a),
						ac = c.minus(a)
					const normal = ab.cross(ac)
					const faceCentroid = V3.add(a, b, c).div(3)
					return [
						faceCentroid.z * normal.z / 2,
						V3.add(
							a.times(2 * a.z + b.z + c.z),
							b.times(a.z + 2 * b.z + c.z),
							c.times(a.z + b.z + 2 * c.z),
						).times(normal.z), // 1/24 factor is done at very end
					]
				}
				const a = edgeWC.a,
					b = edgeWC.b
				const as = a.dot(r1)
				const bs = b.dot(r1)
				const aBase = this.pUV(as, 0)
				const bBase = this.pUV(bs, 0)
				const [v1, c1] = triangleShadowVolumeAndCentroid(a, b, aBase)
				const [v2, c2] = triangleShadowVolumeAndCentroid(bBase, aBase, b)
				return [v1 + v2, c1.plus(c2).div(24)]
			} else if (curveWC instanceof ImplicitCurve) {
				throw new Error()
			} else {
				const sliceAreaAndCentroidZX2TimesDs = (curveT: number) => {
					const p = curveWC.at(curveT)
					const s = p.dot(r1)
					const t = p.dot(u1)
					const area = t * c.z + s * t * r1.z + 1 / 2 * t ** 2 * u1.z
					const ds = -curveWC.tangentAt(curveT).dot(r1)
					return [
						area * ds,
						...V3.add(
							c.times(area),
							r1.times(c.z * s * t + r1.z * s ** 2 * t + 1 / 2 * s * t ** 2 * u1.z),
							u1.times(1 / 2 * c.z * t ** 2 + 1 / 2 * r1.z * s * t ** 2 + 1 / 3 * t ** 3 * u1.z),
						).times(ds),
					]
				}
				const [vol, cx, cy, cz] = glqArray(sliceAreaAndCentroidZX2TimesDs, edgeWC.aT, edgeWC.bT, 4)
				return [vol * this.plane.normal1.z, new V3(cx, cy, cz).times(this.plane.normal1.z)]
			}
		})
		return mergeVolumeAndCentroidZX2Parts(volumeAndCentroidZX2Parts)
	},

	/**
	 * Generic implementation.
	 */
	[ParametricSurface.name](this: ParametricSurface, edges: Edge[]): { centroid: V3; volume: number } {
		const dpdu = this.dpdu()
		const dpdv = this.dpdv()
		const volume = edges.map((edgeWC): [number, V3] => {
			const curveWC = edgeWC.curve
			if (curveWC instanceof ImplicitCurve) {
				throw new Error()
			} else {
				const sliceAreaAndCentroidZX2TimesDs = (curveT: number) => {
					// use curve.tangent not edge.tangent, reverse edges are handled by the integration boundaries
					const pWC = curveWC.at(curveT),
						tangentWC = curveWC.tangentAt(curveT)
					const uvOfPWC = this.uvP(pWC)
					const slice = (t: number) => {
						const p = this.pUV(uvOfPWC.x, t)
						const normal = dpdu(uvOfPWC.x, t).cross(dpdv(uvOfPWC.x, t))
						return p.z * normal.z
					}
					const sliceIntegral0ToPWCT = glqInSteps(slice, 0, uvOfPWC.y, 1)
					// const dt = tangentWC.dot(scalingVector)
					const dt = -M4.forSys(dpdu(uvOfPWC.x, uvOfPWC.y), dpdv(uvOfPWC.x, uvOfPWC.y))
						.inversed()
						.transformVector(tangentWC).x
					const sliceAreaTimesDs = sliceIntegral0ToPWCT * dt
					const slice2 = (t: number) => {
						const p = this.pUV(uvOfPWC.x, t)
						const normal = dpdu(uvOfPWC.x, t).cross(dpdv(uvOfPWC.x, t))
						return p.times(p.z * normal.z)
					}
					const sliceIntegral0ToPWCT2 = glqV3(slice2, 0, uvOfPWC.y)
					// const dt = tangentWC.dot(scalingVector)
					const sliceCentroidZX2TimesDs = sliceIntegral0ToPWCT2.times(dt)
					return [sliceAreaTimesDs, ...sliceCentroidZX2TimesDs.toArray()]
				}
				const [vol, cx, cy, cz] = glqArray(sliceAreaAndCentroidZX2TimesDs, edgeWC.aT, edgeWC.bT, 4)
				return [vol, new V3(cx, cy, cz)]
			}
		})
		return mergeVolumeAndCentroidZX2Parts(volume)
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
		const volume = edges.map((edgeWC): [number, V3] => {
			if (edgeWC.curve instanceof L3) {
				return [0, V3.O]
			} else if (edgeWC.curve instanceof ImplicitCurve) {
				return [0, V3.O]
				// 	const { points, tangents } = edgeWC.curve
				// 	const minT = edgeWC.minT,
				// 		maxT = edgeWC.maxT
				// 	let sum = 0
				// 	const start = Math.ceil(minT + NLA_PRECISION)
				// 	const end = Math.floor(maxT - NLA_PRECISION)
				// 	for (let i = start; i <= end; i++) {
				// 		const at = points[i],
				// 			tangent = tangents[i]
				// 		const area = (at.z + at.rejectedFrom1(upDir1).z) / 2 * at.projectedOn(upDir1).dot(baseVector)
				// 		const scale = tangent.dot(scalingVector)
				// 		sum += area * scale
				// 	}
				// 	const f = (t: number) => {
				// 		const at = edgeWC.curve.at(t),
				// 			tangent = edgeWC.curve.tangentAt(t)
				// 		const area = (at.z + at.rejectedFrom1(upDir1).z) / 2 * at.projectedOn(upDir1).dot(baseVector)
				// 		const scale = tangent.dot(scalingVector)
				// 		return area * scale
				// 	}
				// 	sum += f(minT) * (start - minT - 0.5)
				// 	sum += f(maxT) * (maxT - end - 0.5)
				// 	return sum * Math.sign(edgeWC.deltaT())
			} else {
				const f = (curveT: number) => {
					// use curve.tangent not edge.tangent, reverse edges are handled by the integration boundaries
					const at = edgeWC.curve.at(curveT),
						tangent = edgeWC.curve.tangentAt(curveT)
					const b = at.rejectedFrom1(upDir1)
					const area = at.z * b.to(at).dot(baseVector) / 2 + b.z * b.to(at).dot(baseVector) / 2
					const areaCentroidA = V3.add(at.xy(), b, at).times(at.z * b.to(at).dot(baseVector) / 2 / 3)
					const areaCentroidB = V3.add(at.xy(), b, b.xy()).times(b.z * b.to(at).dot(baseVector) / 2 / 3)
					const scale = tangent.dot(scalingVector)
					return [
						area * scale,
						...areaCentroidA
							.plus(areaCentroidB)
							.times(scale)
							.schur(V(1, 1, 2)),
					]
				}
				const [vol, cx, cy, cz] = glqArray(f, edgeWC.aT, edgeWC.bT, 4)
				return [vol, new V3(cx, cy, cz)]
			}
		})
		return mergeVolumeAndCentroidZX2Parts(volume)
	},

	// volume does scale linearly, so this could be done in the local coordinate system
	// however, shear matrices lead to point-to-plane distances having to be calculated along a vector other than
	// the plane normal
	[RotatedCurveSurface.name](this: RotatedCurveSurface, edges: Edge[]): { volume: number; centroid: V3 } {
		const dpdu = this.dpdu()
		const dpdv = this.dpdv()
		const totalVolume = edges
			.map(edgeWC => {
				const curveWC = edgeWC.curve

				const f = (curveT: number) => {
					const pWC = curveWC.at(curveT),
						tangentWC = curveWC.tangentAt(curveT)
					const uvOfPWC = this.uvP(pWC)
					const pLC = this.matrixInverse.transformPoint(pWC)
					const dpdvAtS0 =
						this instanceof RotatedCurveSurface
							? this.curve.tangentAt(uvOfPWC.y)
							: V(-pLC.z, 0, pLC.lengthXY())
					// const slice = (phi: number) => {
					// 	const p = this.pUV(phi, uvOfPWC.y)
					// 	const normal = dpdu(phi, uvOfPWC.y).cross(dpdv(phi, uvOfPWC.y))
					// 	return p.z * normal.z
					// }
					// const z = this.curve.at(uvOfPWC.y).z
					// const r = this.curve.at(uvOfPWC.y).lengthXY()
					// const pz =
					// 	this.f1.z * r * cos(s) +
					// 	this.f2.z * r * sin(s) +
					// 	this.f3.z * z +
					// 	this.center.z
					// const dpdux = this.f1.x * r * -sin(s) + this.f2.x * r * cos(s)
					// const dpduy = this.f1.y * r * -sin(s) + this.f2.y * r * cos(s)
					// const dpdvx = this.f1.x * dr * cos(s) + this.f2.x * dr * sin(s) + this.f3.x * dz
					// const dpdvy = this.f1.y * dr * cos(s) + this.f2.y * dr * sin(s) + this.f3.y * dz
					// const normalz = dpdux * dpdvy - dpduy * dpdvx
					// result = pz * normalz
					const r = pLC.lengthXY(),
						z = pLC.z
					const dr = dpdvAtS0.x
					const dz = dpdvAtS0.z
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
					const dt = M4.forSys(dpdu(uvOfPWC.x, uvOfPWC.y), dpdv(uvOfPWC.x, uvOfPWC.y))
						.inversed()
						.transformVector(tangentWC).y
					const sliceIntegral0ToPWCS = sliceIntegral(uvOfPWC.x) //- sliceIntegral(0) //(always 0)
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
				const uvOfPWC = this.uvP(pWC)
				const slice = (phi: number) => {
					const p = this.pUV(phi, uvOfPWC.y)
					const normal = dpdu(phi, uvOfPWC.y).cross(dpdv(phi, uvOfPWC.y))
					return p.times(p.z * normal.z)
				}
				const sliceIntegral0ToPWCS = glqV3(slice, 0, uvOfPWC.x)
				const dt = M4.forSys(dpdu(uvOfPWC.x, uvOfPWC.y), dpdv(uvOfPWC.x, uvOfPWC.y))
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
	},
}
ZDirVolumeVisitor[EllipsoidSurface.name] = ZDirVolumeVisitor[RotatedCurveSurface.name]

export function glqV3(f: (x: number) => V3, startT: number, endT: number) {
	return gaussLegendre24Xs
		.reduce((val, currVal, index) => {
			const x = startT + (currVal + 1) / 2 * (endT - startT)
			return val.plus(f(x).times(gaussLegendre24Weights[index]))
		}, V3.O)
		.times((endT - startT) / 2)
}
export function glqArray(f: (x: number) => number[], startT: number, endT: number, numEls = 3) {
	const result = new Array(numEls).fill(0)
	for (let i = 0; i < 24; i++) {
		const x = startT + (gaussLegendre24Xs[i] + 1) / 2 * (endT - startT)
		const fx = f(x)
		for (let j = 0; j < numEls; j++) {
			result[j] += fx[j] * gaussLegendre24Weights[i]
		}
	}
	for (let j = 0; j < numEls; j++) {
		result[j] *= (endT - startT) / 2
	}
	return result
}

function mergeVolumeAndCentroidZX2Parts(volumeAndCentroidZX2Parts: [number, V3][]) {
	const volume = volumeAndCentroidZX2Parts.reduce((result, [volume]) => result + volume, 0)
	const weightedCentroid = V3.add(...volumeAndCentroidZX2Parts.map(([, centroidZX2]) => centroidZX2)).schur(
		new V3(1, 1, 0.5),
	)
	return { volume, centroid: weightedCentroid.div(volume) }
}
