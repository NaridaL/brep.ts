///<reference path="ConicSurface.ts"/>
const ZDirVolumeVisitor: {[className: string]: <T extends Surface>(this: T, allEdges: Edge[]) => {volume: number, centroid: any} } = {
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
	[ConicSurface.name](this: ConicSurface, allEdges: Edge[]): {volume: number, centroid: any} {
		// INT[edge.at; edge.bT] (at(t) DOT dir) * (at(t) - at(t).projectedOn(dir) / 2).z
		const totalVolume = allEdges.map(edge => {
			if (edge.curve instanceof SemiEllipseCurve || edge.curve instanceof HyperbolaCurve || edge.curve instanceof ParabolaCurve) {
				const f = (t) => {
					const at = edge.curve.at(t), tangent = edge.tangentAt(t)
					return (at.z + at.rejectedFrom(this.dir).z) / 2 * at.projectedOn(this.dir).lengthXY() *
						tangent.dot(V3.Z.cross(this.dir).unit())
				}
				// ellipse with normal1 parallel to dir need to be counted negatively so CCW faces result in a positive
				// area
				const sign = edge.curve instanceof SemiEllipseCurve
					? -Math.sign(edge.curve.normal.dot(this.dir))
					: -Math.sign(this.center.to(edge.curve.center).cross(edge.curve.f1).dot(this.dir))
				const val = glqInSteps(f, edge.aT, edge.bT, 1)
				return val * sign
			} else if (edge.curve instanceof L3) {
				return 0
			} else {
				assertNever()
			}
		}).sum()

		return {volume: totalVolume * Math.sign(this.normal.dot(this.dir))}
	},

    [PlaneSurface.name](): {centroid: V3, volume: number} {
        const {centroid, area} = this.calculateArea()
        return {volume: this.surface.plane.normal1.z * centroid.z * area,
            centroid: new V3(centroid.x, centroid.y, centroid.z / 2) }

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
	[SemiCylinderSurface.name](this: SemiCylinderSurface, allEdges: Edge[]): {volume: number, centroid: any} {
		if (V3.Z.cross(this.dir).likeO()) return {volume: 0}
		// the tangent needs to be projected onto a vector which is perpendicular to the volume-slices
		const scalingVector = this.dir.cross(V3.Z).unit()
		// the length of the base of the trapezoid is calculated by dotting with the baseVector
		const baseVector = this.dir.rejectedFrom(V3.Z).unit()
		// INT[edge.at; edge.bT] (at(t) DOT dir1) * (at(t) - at(t).projectedOn(dir) / 2).z
		const totalArea = edges.map(edge => {
			if (edge.curve instanceof SemiEllipseCurve) {
				const f = (t) => {
					// use curve.tangent not edge.tangent, reverse edges are handled by the integration boundaries
					const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t)
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
				// ellipse with normal parallel to dir1 need to be counted negatively so CCW faces result in a positive area
				const sign = -Math.sign(edge.curve.normal.dot(this.dir))
				const val = glqInSteps(f, edge.aT, edge.bT, 1)
				return val * sign
			} else if (edge.curve instanceof L3) {
				return 0
			} else {
				assertNever()
			}
		}).sum()

		return {volume: totalArea * Math.sign(this.baseCurve.normal.dot(this.dir))}
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
	zDirVolume(edges: Edge[]): {volume: number, centroid: any} {
		if (V3.Z.cross(this.dir).likeO()) return {volume: 0}
		// the tangent needs to be projected onto a vector which is perpendicular to the volume-slices
		const scalingVector = this.dir.cross(V3.Z).unit()
		// the length of the base of the trapezoid is calculated by dotting with the baseVector
		const baseVector = this.dir.rejectedFrom(V3.Z).unit()
		// INT[edge.at; edge.bT] (at(t) DOT dir1) * (at(t) - at(t).projectedOn(dir) / 2).z
		const totalArea = edges.map(edge => {
			if (edge.curve instanceof SemiEllipseCurve) {
				const f = (t) => {
					// use curve.tangent not edge.tangent, reverse edges are handled by the integration boundaries
					const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t)
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
				// ellipse with normal parallel to dir1 need to be counted negatively so CCW faces result in a positive area
				const sign = -Math.sign(edge.curve.normal.dot(this.dir))
				const val = glqInSteps(f, edge.aT, edge.bT, 1)
				return val * sign
			} else if (edge.curve instanceof L3) {
				return 0
			} else {
				assertNever()
			}
		}).sum()

		return {volume: totalArea * Math.sign(this.baseCurve.normal.dot(this.dir))}
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
	[CylinderSurface.name](this: CylinderSurface, edges: Edge[]): {volume: number} {
		if (V3.Z.cross(this.dir).likeO()) return {volume: 0}
		// the tangent needs to be projected onto a vector which is perpendicular to the volume-slices
		const scalingVector = this.dir.cross(V3.Z).unit()
		// the length of the base of the trapezoid is calculated by dotting with the baseVector
		const baseVector = this.dir.rejectedFrom(V3.Z).unit()
		// INT[edge.at; edge.bT] (at(t) DOT dir1) * (at(t) - at(t).projectedOn(dir) / 2).z
		console.log("scalingVector", scalingVector.sce)
		const totalArea = edges.map(edge => {
			if (edge.curve instanceof EllipseCurve) {
				const f = (t) => {
					// use curve.tangent not edge.tangent, reverse edges are handled by the integration boundaries
					const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t)
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
				// ellipse with normal parallel to dir1 need to be counted negatively so CCW faces result in a positive area
				const sign = -Math.sign(edge.curve.normal.dot(this.dir))
				const val = glqInSteps(f, edge.aT, edge.bT, 1)
				console.log("edge", edge, val, sign)
				return val * sign
			} else if (edge.curve instanceof L3) {
				return 0
			} else {
				assertNever()
			}
		}).sum()

		return {volume: totalArea * Math.sign(this.baseCurve.normal.dot(this.dir))}
	},

	// volume does scale linearly, so this can be done in the local coordinate system
	// first transform edges with inverse matrix
	// then rotate everything edges so the original world Z dir again points in Z dir
	// now we have a problem because edges which originally  did not cross the seam plane can now be anywhere
	// we need to split the transformed loop along the local seam plane
	// and then sum the zDir volumes of the resulting loops
	[EllipsoidSurface.name](this: EllipsoidSurface, loop: Edge[]): {centroid: V3, volume: number} {
		const angles = this.inverseMatrix.transformVector(V3.Z).toAngles()
		const T = M4.rotateAB(this.inverseMatrix.transformVector(V3.Z), V3.Z).times(M4.rotateZ(-angles.phi)).times(this.inverseMatrix)
		function calc(loop) {
			let totalVolume = 0
			assert(V3.Z.isParallelTo(T.transformVector(V3.Z)))
			//const zDistanceFactor = toT.transformVector(V3.Z).length()
			loop.map(edge => edge.transform(T)).forEach((edge, edgeIndex, edges) => {
				const nextEdgeIndex = (edgeIndex + 1) % edges.length, nextEdge = edges[nextEdgeIndex]

				function f(t) {
					const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t)
					const r = at.lengthXY()
					const at2d = at.withElement('z', 0)
					const angleAdjusted = (at.angleXY() + TAU - NLA_PRECISION) % TAU + NLA_PRECISION
					const result = angleAdjusted * Math.sqrt(1 - r * r) * r * Math.abs(tangent.dot(at2d.unit())) * Math.sign(tangent.z)
					//console.log("at2d", at2d.sce, "result", result, 'angle', angleAdjusted, ' edge.tangentAt(t).dot(at2d.unit())', edge.tangentAt(t).dot(at2d.unit()))
					return result
				}

				const volume = gaussLegendreQuadrature24(f, edge.aT, edge.bT)
				console.log("edge", edge, "volume", volume)
				totalVolume += volume
			})
			return totalVolume
		}
		const [front, back] = EllipsoidSurface.splitOnPlaneLoop(loop.map(edge => edge.transform(T)), ccw)
		const localVolume = calc(front, PI) + calc(back, -PI)

		return {area: localVolume * this.f1.dot(this.f2.cross(this.f3)), centroid: undefined}
	},
	zDirVolumeForLoop2(loop: Edge[]): number {
		const angles = this.inverseMatrix.getZ().toAngles()
		const T = M4.rotateY(-angles.theta).times(M4.rotateZ(-angles.phi)).times(this.inverseMatrix)
		const rot90x = M4.rotateX(PI / 2)
		let totalVolume = 0
		assert(V3.X.isParallelTo(T.transformVector(V3.Z)))
		//const zDistanceFactor = toT.transformVector(V3.Z).length()
		loop.map(edge => edge.transform(T)).forEach((edge, edgeIndex, edges) => {
			const nextEdgeIndex = (edgeIndex + 1) % edges.length, nextEdge = edges[nextEdgeIndex]
			function f (t) {
				const at2d = edge.curve.at(t).withElement('x', 0)
				const result = 1 / 3 * (1 - (at2d.y ** 2 + at2d.z ** 2)) * edge.tangentAt(t).dot(rot90x.transformVector(at2d.unit()))
				console.log("at2d", at2d.sce, "result", result)
				return result
			}
			//if (edge.)
			if (edge.b.like(V3.X)) {
				const angleDiff = (edge.bDir.angleRelativeNormal(nextEdge.aDir, V3.X) + 2 * PI) % (2 * PI)
				totalVolume += 2 / 3 * angleDiff
				console.log("xaa")
			}
			if (edge.b.like(V3.X.negated())) {
				const angleDiff = (edge.bDir.angleRelativeNormal(nextEdge.aDir, V3.X) + 2 * PI) % (2 * PI)
				totalVolume += 2 / 3 * angleDiff
				console.log("xbb")
			}
			const volume = gaussLegendreQuadrature24(f, edge.aT, edge.bT)
			console.log("edge", edge, "volume", volume)
			totalVolume += volume
		})

		return totalVolume * this.f1.dot(this.f2.cross(this.f3))
	},

	// volume does scale linearly, so this can be done in the local coordinate system
	// first transform edges with inverse matrix
	// then rotate everything edges so the original world Z dir again points in Z dir
	// now we have a problem because edges which originally  did not cross the seam plane can now be anywhere
	// we need to split the transformed loop along the local seam plane
	// and then sum the zDir volumes of the resulting loops
	[SemiEllipsoidSurface.name](loop: Edge[]): {volume: number, centroid: V3} {
		const angles = this.inverseMatrix.transformVector(V3.Z).toAngles()
		const T = M4.rotateAB(this.inverseMatrix.transformVector(V3.Z), V3.Z).times(M4.rotateZ(-angles.phi)).times(this.inverseMatrix)
		function calc(loop) {
			let totalVolume = 0
			assert(V3.Z.isParallelTo(T.transformVector(V3.Z)))
			//const zDistanceFactor = toT.transformVector(V3.Z).length()
			loop.map(edge => edge.transform(T)).forEach((edge, edgeIndex, edges) => {
				const nextEdgeIndex = (edgeIndex + 1) % edges.length, nextEdge = edges[nextEdgeIndex]

				function f(t) {
					const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t)
					const r = at.lengthXY()
					const at2d = at.withElement('z', 0)
					const angleAdjusted = (at.angleXY() + TAU - NLA_PRECISION) % TAU + NLA_PRECISION
					const result = angleAdjusted * Math.sqrt(1 - r * r) * r * Math.abs(tangent.dot(at2d.unit())) * Math.sign(tangent.z)
					//console.log("at2d", at2d.sce, "result", result, 'angle', angleAdjusted, ' edge.tangentAt(t).dot(at2d.unit())', edge.tangentAt(t).dot(at2d.unit()))
					return result
				}

				const volume = gaussLegendreQuadrature24(f, edge.aT, edge.bT)
				totalVolume += volume
			})
			return totalVolume
		}
		const [front, back] = SemiEllipsoidSurface.splitOnPlaneLoop(loop.map(edge => edge.transform(T)), ccw)
		const localVolume = calc(front, PI) + calc(back, -PI)

		return {volume: localVolume * this.f1.dot(this.f2.cross(this.f3)), centroid: undefined}
	},
}