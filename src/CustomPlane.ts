import { callsce, min, V3 } from "ts3dutils"

import * as chroma from "chroma.ts"
import { GL_COLOR, GL_COLOR_BLACK } from "tsgl"
import { getGlobalId, L3, P3, PlaneSurface } from "."

export class CustomPlane extends P3 {
  readonly up: V3
  readonly right: V3
  readonly vMin: number
  readonly vMax: number
  readonly uMin: number
  readonly uMax: number
  readonly color: GL_COLOR
  readonly name: string

  constructor(
    anchor: V3,
    right: V3,
    up: V3,
    name: string = "CustomPlane" + getGlobalId(),
    color: GL_COLOR = chroma.random().gl(),
    uMin: number = -500,
    uMax: number = 500,
    vMin: number = -500,
    vMax: number = 500,
  ) {
    const { normal1, w } = P3.forAnchorAndPlaneVectors(anchor, right, up)
    super(normal1, w)
    this.up = up
    this.right = right
    this.uMin = uMin
    this.uMax = uMax
    this.vMin = vMin
    this.vMax = vMax
    this.name = name
    this.color = color
  }

  get plane() {
    return this
  }

  toPlaneSurface() {
    return new PlaneSurface(this, this.right, this.up)
  }

  toSource() {
    return callsce(
      "new CustomPlane",
      this.anchor,
      this.right,
      this.up,
      this.name,
      this.color,
      this.uMin,
      this.uMax,
      this.vMin,
      this.vMax,
    )
  }

  static forPlane(plane: P3, color: GL_COLOR = GL_COLOR_BLACK, name?: string) {
    //assert(!name)
    const up = plane.normal1.getPerpendicular().unit(),
      right = up.cross(plane.normal1)
    return new CustomPlane(plane.anchor, right, up, name, color)
  }

  static fromPlaneSurface(surface: PlaneSurface) {
    return new CustomPlane(
      surface.plane.anchor,
      surface.right,
      surface.up,
      "genCustomPlane" + getGlobalId(),
    )
  }

  distanceTo(line: L3, mindist: number) {
    return min(
      [
        new L3(this.anchor.plus(this.right.times(this.uMin)), this.up),
        new L3(this.anchor.plus(this.right.times(this.uMax)), this.up),
        new L3(this.anchor.plus(this.up.times(this.vMin)), this.right),
        new L3(this.anchor.plus(this.up.times(this.vMax)), this.right),
      ].map((line2, line2Index): number => {
        const info = line2.infoClosestToLine(line)
        if (
          (isNaN(info.t) || // parallel LINES
            (line2Index < 2 && this.vMin <= info.t && info.t <= this.vMax) ||
            (line2Index >= 2 && this.uMin <= info.t && info.t <= this.uMax)) &&
          info.distance <= mindist
        ) {
          return info.s
        } else {
          return Infinity
        }
      }),
    )
  }

  distanceTo2(line: L3, mindist: number) {
    return min(
      [
        new L3(this.anchor.plus(this.right.times(this.uMin)), this.up),
        new L3(this.anchor.plus(this.right.times(this.uMax)), this.up),
        new L3(this.anchor.plus(this.up.times(this.vMin)), this.right),
        new L3(this.anchor.plus(this.up.times(this.vMax)), this.right),
      ].map((line2, line2Index) => {
        const info = line2.infoClosestToLine(line)
        if (
          (isNaN(info.t) || // parallel LINES
            (line2Index < 2 && this.vMin <= info.t && info.t <= this.vMax) ||
            (line2Index >= 2 && this.uMin <= info.t && info.t <= this.uMax)) &&
          info.distance <= mindist
        ) {
          return info.distance
        } else {
          return Infinity
        }
      }),
    )
  }
}
