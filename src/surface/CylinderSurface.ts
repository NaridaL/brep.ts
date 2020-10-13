import {
  assert,
  assertInst,
  assertVectors,
  eq0,
  hasConstructor,
  M4,
  pqFormula,
  TAU,
  V3,
} from "ts3dutils"

import {
  BezierCurve,
  Curve,
  Edge,
  EllipseCurve,
  ImplicitSurface,
  L3,
  OUTSIDE,
  P3,
  PointVsFace,
  ProjectedCurveSurface,
  Surface,
} from ".."

import { sign } from "../math"

export class CylinderSurface
  extends ProjectedCurveSurface
  implements ImplicitSurface {
  static readonly UNIT = new CylinderSurface(
    EllipseCurve.UNIT,
    V3.Z,
    undefined,
    undefined,
    0,
    1,
  )
  readonly matrix: M4
  readonly matrixInverse: M4
  readonly pLCNormalWCMatrix: M4
  readonly pWCNormalWCMatrix: M4
  readonly normalDir: number
  // @ts-ignore
  // readonly baseCurve: EllipseCurve

  constructor(
    readonly baseCurve: EllipseCurve,
    dir1: V3,
    uMin: number = baseCurve.tMin,
    uMax: number = baseCurve.tMax,
    zMin = -Infinity,
    zMax = Infinity,
  ) {
    super(baseCurve, dir1, uMin, uMax, zMin, zMax)
    assertInst(EllipseCurve, baseCurve)
    //assert(!baseCurve.normal1.isPerpendicularTo(dir1), !baseCurve.normal1.isPerpendicularTo(dir1))
    this.matrix = M4.forSys(baseCurve.f1, baseCurve.f2, dir1, baseCurve.center)
    this.matrixInverse = this.matrix.inversed()
    this.normalDir = sign(this.baseCurve.normal.dot(this.dir))
    this.pLCNormalWCMatrix = this.matrix
      .as3x3()
      .inversed()
      .transposed()
      .scale(this.normalDir)
    this.pWCNormalWCMatrix = this.pLCNormalWCMatrix.times(this.matrixInverse)
  }

  static semicylinder(
    radius: number,
    sMin?: number,
    sMax?: number,
    tMin?: number,
    tMax?: number,
  ): CylinderSurface {
    return new CylinderSurface(
      new EllipseCurve(V3.O, new V3(radius, 0, 0), new V3(0, radius, 0)),
      V3.Z,
      sMin,
      sMax,
      tMin,
      tMax,
    )
  }

  /**
   *
   * @param anchorLC
   * @param dirLC not necessarily unit
   */
  static unitISLineTs(anchorLC: V3, dirLC: V3): number[] {
    const { x: ax, y: ay } = anchorLC
    const { x: dx, y: dy } = dirLC

    // this cylinder: x² + y² = 1
    // line: p = anchorLC + t * dirLC
    // split line equation into 3 component equations, insert into cylinder equation
    // x = ax + t * dx
    // y = ay + t * dy
    // (ax² + 2 ax t dx + t²dx²) + (ay² + 2 ay t dy + t²dy²) = 1
    // transform to form (a t² + b t + c = 0) and solve with pqFormula
    const a = dx ** 2 + dy ** 2
    const b = 2 * (ax * dx + ay * dy)
    const c = ax ** 2 + ay ** 2 - 1
    return pqFormula(b / a, c / a).filter((t) =>
      EllipseCurve.XYLCValid(new V3(ax + dx * t, ay + dy * t, 0)),
    )
  }

  normalP(p: V3): V3 {
    return this.pLCNormalWCMatrix
      .transformVector(this.matrixInverse.transformPoint(p).xy())
      .unit()
  }

  loopContainsPoint(loop: Edge[], p: V3): PointVsFace {
    assertVectors(p)
    if (!this.containsPoint(p)) return OUTSIDE
    const line = new L3(p, this.dir.unit())
    const lineOut = this.dir.cross(this.normalP(p))
    return Surface.loopContainsPointGeneral(loop, p, line, lineOut)
  }

  isTsForLine(line: L3) {
    assertInst(L3, line)
    // transforming line manually has advantage that dir1 will not be renormalized,
    // meaning that calculated values t for localLine are directly transferable to line
    const dirLC = this.matrixInverse.transformVector(line.dir1)
    if (dirLC.isParallelTo(V3.Z)) {
      // line is parallel to this.dir
      return []
    }
    const anchorLC = this.matrixInverse.transformPoint(line.anchor)
    assert(
      !CylinderSurface.unitISLineTs(anchorLC, dirLC).length ||
        !isNaN(CylinderSurface.unitISLineTs(anchorLC, dirLC)[0]),
      "sad " + dirLC,
    )
    return CylinderSurface.unitISLineTs(anchorLC, dirLC)
  }

  isCoplanarTo(surface: Surface): surface is CylinderSurface {
    return (
      this == surface ||
      (hasConstructor(surface, CylinderSurface) &&
        this.dir.isParallelTo(surface.dir) &&
        this.containsEllipse(surface.baseCurve, false))
    )
  }

  like(surface: Surface): boolean {
    if (!this.isCoplanarTo(surface)) return false
    // normals need to point in the same direction (outwards or inwards) for both
    const thisFacesOut = 0 < this.baseCurve.normal.dot(this.dir)
    const objectFacesOut = 0 < surface.baseCurve.normal.dot(surface.dir)
    return thisFacesOut == objectFacesOut
  }

  containsEllipse(ellipse: EllipseCurve, checkAABB: boolean = true) {
    const projEllipse = ellipse.transform(
      M4.project(this.baseCurve.getPlane(), this.dir),
    )
    return this.baseCurve == ellipse || this.baseCurve.isColinearTo(projEllipse)
  }

  containsCurve(curve: Curve) {
    if (curve instanceof L3) {
      return this.containsLine(curve)
    } else if (curve instanceof EllipseCurve) {
      return this.containsEllipse(curve)
    } else if (curve instanceof BezierCurve) {
      return false
    } else {
      return super.containsCurve(curve)
    }
  }

  implicitFunction() {
    return (pWC: V3) => {
      const pLC = this.matrixInverse.transformPoint(pWC)
      return (pLC.lengthXY() - 1) * this.normalDir
    }
  }

  didp(pWC: V3) {
    const pLC = this.matrixInverse.transformPoint(pWC)
    const pLCLengthXY = pLC.lengthXY()
    const didpLC = new V3(pLC.x / pLCLengthXY, pLC.y / pLCLengthXY, 0)
    return this.pLCNormalWCMatrix.transformVector(didpLC)
  }

  containsPoint(pWC: V3): boolean {
    const pLC = this.matrixInverse.transformPoint(pWC)
    return this.baseCurve.isValidT(
      EllipseCurve.XYLCPointT(pLC, this.uMin, this.uMax),
    )
  }

  uvP(pWC: V3): V3 {
    assert(arguments.length == 1)
    const pLC = this.matrixInverse.transformPoint(pWC)
    const u = EllipseCurve.XYLCPointT(pLC, this.vMin, this.vMax)
    return new V3(u, pLC.z, 0)
  }

  isCurvesWithSurface(surface2: Surface): Curve[] {
    if (surface2 instanceof ProjectedCurveSurface) {
      if (surface2.dir.isParallelTo(this.dir)) {
        const projectedCurve = surface2.baseCurve.transform(
          M4.project(this.baseCurve.getPlane(), this.dir),
        )
        return this.baseCurve.isInfosWithCurve(projectedCurve).map((info) => {
          const lineDir =
            sign(
              this.normalP(info.p)
                .cross(surface2.normalP(info.p))
                .dot(this.dir),
            ) || 1
          return new L3(info.p, this.dir.times(lineDir))
        })
      }
    }
    if (surface2 instanceof CylinderSurface) {
      if (eq0(this.getCenterLine().distanceToLine(surface2.getCenterLine()))) {
        throw new Error()
      }
    }
    return super.isCurvesWithSurface(surface2)
  }

  getCenterLine(): L3 {
    return new L3(this.baseCurve.center, this.dir)
  }

  facesOutwards(): boolean {
    return this.baseCurve.normal.dot(this.dir) > 0
  }

  getSeamPlane(): P3 {
    let normal = this.baseCurve.f1.cross(this.dir)
    normal = normal.times(-sign(normal.dot(this.baseCurve.f2)))
    return P3.normalOnAnchor(normal, this.baseCurve.center)
  }

  clipCurves(curves: Curve[]): Curve[] {
    return curves.flatMap((curve) => curve.clipPlane(this.getSeamPlane()))
  }
}

CylinderSurface.prototype.uStep = TAU / 32
CylinderSurface.prototype.vStep = 256
