import { int, V3 } from 'ts3dutils'

import chroma from 'chroma-js'
import { GL_COLOR } from 'tsgl'
import { getGlobalId, L3, P3, PlaneSurface } from './index'

export class CustomPlane extends P3 {
	readonly up: V3
	readonly right: V3
	readonly tMin: number
	readonly tMax: number
	readonly sMin: number
	readonly sMax: number
	readonly color: GL_COLOR
	readonly name: string

	constructor(
		anchor: V3,
		right: V3,
		up: V3,
		name: string,
		color: GL_COLOR = chroma.random().gl(),
		rightStart: number = -500,
		rightEnd: number = 500,
		upStart: number = -500,
		upEnd: number = 500,
	) {
		const { normal1, w } = P3.forAnchorAndPlaneVectors(anchor, right, up)
		super(normal1, w)
		this.up = up
		this.right = right
		this.sMin = rightStart
		this.sMax = rightEnd
		this.tMin = upStart
		this.tMax = upEnd
		this.name = name
		this.color = color
	}

	get plane() {
		return this
	}

	toPlaneSurface() {
		return new PlaneSurface(this, this.right, this.up)
	}

	static forPlane(plane: P3, color: GL_COLOR, name?: string) {
		//assert(!name)
		const up = plane.normal1.getPerpendicular().unit(),
			right = up.cross(plane.normal1)
		return new CustomPlane(plane.anchor, right, up, name, color)
	}

	static fromPlaneSurface(surface: PlaneSurface) {
		return new CustomPlane(surface.plane.anchor, surface.right, surface.up, 'genCustomPlane' + getGlobalId())
	}

	distanceTo(line: L3, mindist: number) {
		return [
			new L3(this.anchor.plus(this.right.times(this.sMin)), this.up),
			new L3(this.anchor.plus(this.right.times(this.sMax)), this.up),
			new L3(this.anchor.plus(this.up.times(this.tMin)), this.right),
			new L3(this.anchor.plus(this.up.times(this.tMax)), this.right),
		]
			.map((line2, line2Index) => {
				const info = line2.infoClosestToLine(line)
				if (
					(isNaN(info.t) || // parallel LINES
						(line2Index < 2 && this.tMin <= info.t && info.t <= this.tMax) ||
						(line2Index >= 2 && this.sMin <= info.t && info.t <= this.sMax)) &&
					info.distance <= mindist
				) {
					return info.s
				} else {
					return Infinity
				}
			})
			.min()
	}

	distanceTo2(line: L3, mindist: number) {
		return [
			new L3(this.anchor.plus(this.right.times(this.sMin)), this.up),
			new L3(this.anchor.plus(this.right.times(this.sMax)), this.up),
			new L3(this.anchor.plus(this.up.times(this.tMin)), this.right),
			new L3(this.anchor.plus(this.up.times(this.tMax)), this.right),
		]
			.map((line2, line2Index) => {
				const info = line2.infoClosestToLine(line)
				if (
					(isNaN(info.t) || // parallel LINES
						(line2Index < 2 && this.tMin <= info.t && info.t <= this.tMax) ||
						(line2Index >= 2 && this.sMin <= info.t && info.t <= this.sMax)) &&
					info.distance <= mindist
				) {
					return info.distance
				} else {
					return Infinity
				}
			})
			.min()
	}
}
