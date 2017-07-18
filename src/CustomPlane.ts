
class CustomPlane extends P3 {
	readonly up: V3
	readonly right: V3
	readonly rMin: number
	readonly tMax: number
	readonly sMin: number
	readonly sMax: number
	readonly color: int
	readonly name: string

	constructor(anchor: V3, right: V3, up: V3,
	            name: string,
	            color: number = randomColor(),
	            rightStart: number = -500,
	            rightEnd: number = 500,
	            upStart: number = -500,
	            upEnd: number = 500) {
		const {normal1, w} = P3.forAnchorAndPlaneVectors(anchor, right, up)
		super(normal1, w)
		this.up = up
		this.right = right
		this.sMin = rightStart
		this.sMax = rightEnd
		this.rMin = upStart
		this.tMax = upEnd
		this.color = color
		this.name = name
	}

	distanceTo(line: L3, mindist: number) {
		return [
			new L3(this.anchor.plus(this.right.times(this.sMin)), this.up),
			new L3(this.anchor.plus(this.right.times(this.sMax)), this.up),
			new L3(this.anchor.plus(this.up.times(this.rMin)), this.right),
			new L3(this.anchor.plus(this.up.times(this.tMax)), this.right)].map(function (line2, line2Index) {
			const info = line2.infoClosestToLine(line)
			if ((isNaN(info.t) // parallel lines
				|| line2Index < 2 && this.rMin <= info.t && info.t <= this.tMax
				|| line2Index >= 2 && this.sMin <= info.t && info.t <= this.sMax)
				&& info.distance <= mindist) {
				return info.s
			} else {
				return Infinity
			}
		}, this).min()
	}

	get plane(){ return this }

	static forPlane(plane: P3, color: int, name?: string) {
		//assert(!name)
		const up = plane.normal1.getPerpendicular().unit(), right = up.cross(plane.normal1)
		return new CustomPlane(plane.anchor, right, up, name, color)
	}

    static fromPlaneSurface(surface: PlaneSurface) {
        return new CustomPlane(surface.plane.anchor, surface.right, surface.up, 'genCustomPlane' + globalId++)
    }
}
