///<reference path="ConicSurface.ts"/>
///<reference path="CylinderSurface.ts"/>
///<reference path="EllipsoidSurface.ts"/>
///<reference path="SemiCylinderSurface.ts"/>
///<reference path="../ignore/ProjectedCurveSurface.ts"/>
const CalculateAreaVisitor: {[className: string]: <T extends Surface>(this: T, allEdges: Edge[]) => number } = {
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
				// ellipse with normal parallel to dir1 need to be counted negatively so CCW faces result in a positive area
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
				console.log("edge", edge, val)
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
				// ellipse with normal1 parallel to dir1 need to be counted negatively so CCW faces result in a positive area
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
				// ellipse with normal parallel to dir1 need to be counted negatively so CCW faces result in a positive area
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
	}
}