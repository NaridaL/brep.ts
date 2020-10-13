import {
  arrayFromFunction,
  assert,
  assertf,
  assertInst,
  assertNumbers,
  assertVectors,
  between,
  checkDerivate,
  clamp,
  eq,
  eq0,
  fuzzyBetween,
  gaussLegendreQuadrature24,
  getIntervals,
  getRoots,
  glqInSteps,
  hasConstructor,
  le,
  lt,
  M4,
  MINUS,
  newtonIterate,
  NLA_PRECISION,
  pqFormula,
  snap,
  sum,
  toSource,
  V,
  V3,
} from "ts3dutils"

import {
  Curve,
  CylinderSurface,
  Edge,
  EllipseCurve,
  getExtremePointsHelper,
  ImplicitSurface,
  L3,
  P3,
  ParametricSurface,
  PICurve,
  PlaneSurface,
  PointVsFace,
  ProjectedCurveSurface,
  Surface,
} from ".."

import { abs, cos, max, min, PI, sign, sin, sqrt } from "../math"

class ArrayExt {}
declare global {
  interface Array<T> extends ArrayExt {}
}
export class EllipsoidSurface
  extends ParametricSurface
  implements ImplicitSurface {
  static readonly UNIT = new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z)
  readonly matrix: M4
  readonly matrixInverse: M4
  readonly pLCNormalWCMatrix: M4
  readonly pWCNormalWCMatrix: M4
  readonly normalDir: number // -1 | 1

  constructor(
    readonly center: V3,
    readonly f1: V3,
    readonly f2: V3,
    readonly f3: V3,
    uMin: number = 0,
    uMax: number = PI,
    vMin: number = -PI / 2,
    vMax: number = PI / 2,
  ) {
    super(uMin, uMax, vMin, vMax)
    assert(0 <= uMin && uMin <= PI)
    assert(0 <= uMax && uMax <= PI)
    assert(-PI / 2 <= vMin && vMin <= PI / 2)
    assert(-PI / 2 <= vMax && vMax <= PI / 2)
    assertVectors(center, f1, f2, f3)
    this.matrix = M4.forSys(f1, f2, f3, center)
    this.matrixInverse = this.matrix.inversed()
    this.normalDir = sign(this.f1.cross(this.f2).dot(this.f3))
    this.pLCNormalWCMatrix = this.matrix
      .as3x3()
      .inversed()
      .transposed()
      .scale(this.normalDir)
    this.pWCNormalWCMatrix = this.pLCNormalWCMatrix.times(this.matrixInverse)
  }

  static unitArea(contour: Edge[]) {
    const totalArea = sum(
      contour.map((edge) => {
        if (edge.curve instanceof PICurve) {
          const points = edge.curve.calcSegmentPoints(
            edge.aT,
            edge.bT,
            edge.a,
            edge.b,
            edge.aT > edge.bT,
            true,
          )
          let sum = 0
          for (let i = 0; i < points.length - 1; i++) {
            const p = points[i],
              ppp = points[i + 1]
            sum += ((abs(p.angleXY()) + abs(ppp.angleXY())) / 2) * (ppp.z - p.z)
          }
          return sum
        } else if (edge.curve instanceof EllipseCurve) {
          const f = (t: number) => {
            const at = edge.curve.at(t),
              tangent = edge.curve.tangentAt(t)
            const angleXY = abs(at.angleXY())
            //const arcLength = angleXY * Math.sqrt(1 - at.z ** 2) ( == at.lengthXY())
            //const scaling = tangent.z / at.lengthXY()
            return angleXY * tangent.z
          }
          const val = glqInSteps(f, edge.aT, edge.bT, 1)
          return val
        } else {
          throw new Error()
        }
      }),
    )
    return totalArea
  }

  /**
   * unit sphere: x² + y² + z² = 1
   * line: p = anchor + t * dir |^2
   * p² = (anchor + t * dir)^2
   * 1 == (anchor + t * dir)^2
   * 1 == anchor DOT anchor + 2 * anchor * t * dir + t² * dir DOT dir
   */
  static unitISTsWithLine(anchor: V3, dir: V3): number[] {
    // for 0 = a t² + b t + c
    const a = dir.dot(dir)
    const b = 2 * anchor.dot(dir)
    const c = anchor.dot(anchor) - 1
    return pqFormula(b / a, c / a).filter((t) => le(0, anchor.y + t * dir.y))
  }

  /**
   * unit sphere: x² + y² + z² = 1
   * plane: normal1 DOT p = w
   */
  static unitISCurvesWithPlane(plane: P3): EllipseCurve[] {
    const distPlaneCenter = Math.abs(plane.w)
    if (lt(distPlaneCenter, 1)) {
      // result is a circle
      // radius of circle: imagine right angled triangle (origin -> center of intersection circle -> point on
      // intersection circle) pythagoras: 1² == distPlaneCenter² + isCircleRadius² => isCircleRadius == sqrt(1 -
      // distPlaneCenter²)
      const isCircleRadius = Math.sqrt(1 - distPlaneCenter ** 2)
      const anchorY = plane.normal1.y * plane.w
      const d = abs(distPlaneCenter * isCircleRadius)
      if (le(anchorY, -d) && !eq0(distPlaneCenter)) {
        return []
      } else if (le(anchorY, 0) && !plane.normal1.isParallelTo(V3.Y)) {
        const f1 = plane.normal1.isParallelTo(V3.Y)
          ? V3.Z
          : plane.normal1.cross(V3.Y).toLength(isCircleRadius)
        const f2 = f1.cross(plane.normal1)
        const minEta = -anchorY / f2.y,
          minT = max(0, Math.asin(minEta))
        return [new EllipseCurve(plane.anchor, f1, f2, minT, PI - minT)]
      } else {
        const f2 = (plane.normal1.isParallelTo(V3.Y)
          ? V3.X
          : plane.normal1.cross(V3.Y)
        ).toLength(isCircleRadius)
        const f1 = f2.cross(plane.normal1)
        const minXi = eq0(f1.y) ? -1 : -anchorY / f1.y,
          maxT = Math.acos(max(-1, minXi - NLA_PRECISION))
        return [
          new EllipseCurve(plane.anchor, f1.negated(), f2, PI - maxT, PI),
          new EllipseCurve(plane.anchor, f1, f2.negated(), 0, maxT),
        ]
      }
    } else {
      return []
    }
  }

  static unitISCurvesWithEllipsoidSurface(surface: EllipsoidSurface): Curve[] {
    if (surface.isSphere()) {
      const surfaceRadius = surface.f1.length()
      const surfaceCenterDist = surface.center.length()
      if (
        le(1, surfaceCenterDist - surfaceRadius) ||
        le(surfaceCenterDist + surfaceRadius, 1) ||
        le(surfaceCenterDist - surfaceRadius, -1)
      ) {
        return []
      } else {
        // origin, surface.center and points on the intersection curves form a triangle.
        // the height on the segment origin - surface.center is the radius of the is curves
        // the distance from the origin to the lot point is the distance to the intersection plane
        function heron(a: number, b: number, c: number) {
          const p = (a + b + c) / 2
          return sqrt(p * (p - a) * (p - b) * (p - c))
        }

        const triangleArea = heron(1, surfaceRadius, surfaceCenterDist)
        const radius = (triangleArea * 2) / surfaceCenterDist
        const isCurvesCenterDist =
          sign(1 + surfaceCenterDist ** 2 - surfaceRadius ** 2) *
          sqrt(1 - radius ** 2)
        const plane = new P3(surface.center.unit(), isCurvesCenterDist)
        return EllipsoidSurface.unitISCurvesWithPlane(plane.flipped())
      }
    }
    throw new Error()
  }

  static unitISCurvesWithCylinderSurface(
    surface: CylinderSurface,
  ): EllipseCurve[] {
    if (new L3(surface.baseCurve.center, surface.dir).containsPoint(V3.O)) {
      const projEllipse = surface.baseCurve.transform(
        M4.project(new P3(surface.dir, 0)),
      )
      const f1Length = projEllipse.f1.length(),
        f2Length = projEllipse.f2.length()
      if (lt(1, min(f1Length, f2Length))) return []
      if (projEllipse.isCircular()) {
        const distISCurveCenter = Math.sqrt(1 - min(1, f1Length) ** 2)
        const isCurveCenter = (surface.dir.y < 0
          ? surface.dir.negated()
          : surface.dir
        ).times(distISCurveCenter)
        // isCurve.at(t).y = isCurveCenter.y + projEllipse.f1.y * cos(t) + projEllipse.f2.y * sin(t) = 0
        return [new EllipseCurve(isCurveCenter, projEllipse.f1, projEllipse.f2)]
      }
    }
    throw new Error()
  }

  static sphere(radius: number, center: V3 = V3.O): EllipsoidSurface {
    assertNumbers(radius)
    return new EllipsoidSurface(
      center,
      new V3(radius, 0, 0),
      new V3(0, radius, 0),
      new V3(0, 0, radius),
    )
  }

  /**
   * x²/a² + y²/b² + z²/c² = 1
   */
  static forABC(
    a: number,
    b: number,
    c: number,
    center: V3 = V3.O,
  ): EllipsoidSurface {
    return new EllipsoidSurface(
      center,
      new V3(a, 0, 0),
      new V3(0, b, 0),
      new V3(0, 0, c),
    )
  }

  static calculateAreaSpheroid(a: V3, b: V3, c: V3, edges: Edge[]): number {
    assertf(() => a.isPerpendicularTo(b))
    assertf(() => b.isPerpendicularTo(c))
    assertf(() => c.isPerpendicularTo(a))

    // handling discontinuities:
    // option 1: check for intersections with baseline, if there are any integrate parts separetely
    // "rotate" the edge so that there are no overlaps
    const matrix = M4.forSys(a, b, c),
      matrixInverse = matrix.inversed()
    const circleRadius = a.length()
    const c1 = c.unit()
    const totalArea = sum(
      edges.map((edge) => {
        if (edge.curve instanceof EllipseCurve) {
          const f = (t: number) => {
            const at = edge.curve.at(t),
              tangent = edge.tangentAt(t)
            const localAt = matrixInverse.transformPoint(at)
            const angleXY = localAt.angleXY()
            const arcLength =
              angleXY * circleRadius * Math.sqrt(1 + localAt.z ** 2)
            const scaling = Math.sqrt(1 + c1.dot(tangent) ** 2)
            return arcLength * scaling
          }
          const val = glqInSteps(f, edge.aT, edge.bT, 1)
          return val
        } else {
          throw new Error()
        }
      }),
    )

    return totalArea
  }

  getConstructorParameters(): any[] {
    return [
      this.center,
      this.f1,
      this.f2,
      this.f3,
      this.uMin,
      this.uMax,
      this.vMin,
      this.vMax,
    ]
  }

  equals(obj: any): boolean {
    return (
      this == obj ||
      (Object.getPrototypeOf(obj) == this.constructor.prototype &&
        this.matrix.equals(obj.matrix))
    )
  }

  edgeLoopCCW(loop: Edge[]): boolean {
    return (
      EllipsoidSurface.unitArea(
        loop.map((edge) => edge.transform(this.matrixInverse)),
      ) > 0
    )
    //let totalAngle = 0
    //for (let i = 0; i < contour.length; i++) {
    //    const ipp = (i + 1) % contour.length
    //    const edge = contour[i], nextEdge = contour[ipp]
    //    totalAngle += edge.bDir.angleRelativeNormal(nextEdge.aDir, this.normalP(edge.b))
    //}
    //return le(0, totalAngle)
  }

  like(object: any) {
    if (!this.isCoplanarTo(object)) return false
    // normals need to point in the same direction (outwards or inwards) for both
    return this.matrix.determinant3() * object.matrix.determinant3() > 0
  }

  rootPoints() {}

  toMesh() {
    return ParametricSurface.prototype.toMesh.call(this)
  }

  clipCurves(curves: Curve[]): Curve[] {
    return curves.flatMap((curve) => curve.clipPlane(this.getSeamPlane()))
  }
  dpdu(): (u: number, v: number) => V3 {
    // dp(u, v) = new V3(cos(t) * cos(s), cos(t) * sin(s), sin(t)
    return (u: number, v: number) =>
      this.matrix.transformVector(new V3(cos(v) * -sin(u), cos(v) * cos(u), 0))
  }

  dpdv(): (u: number, v: number) => V3 {
    return (u: number, v: number) =>
      this.matrix.transformVector(
        new V3(-sin(v) * cos(u), -sin(v) * sin(u), cos(v)),
      )
  }

  isCurvesWithPCS(surface: ProjectedCurveSurface): Curve[] {
    let curves2 = ParametricSurface.isCurvesParametricImplicitSurface(
      surface,
      this,
      0.1,
      0.1 / surface.dir.length(),
      0.05,
    )
    curves2 = this.clipCurves(curves2)
    return curves2
  }

  isCurvesWithPCSSmart(surface: ProjectedCurveSurface): Curve[] {
    const surfaceLC = surface.transform(this.matrixInverse)
    //const lcMinZ0RelO =
    const baseCurveLC = surfaceLC.baseCurve.project(new P3(surfaceLC.dir, 0))
    const ists = baseCurveLC.isTsWithSurface(EllipsoidSurface.UNIT)
    const insideIntervals = getIntervals(
      ists,
      baseCurveLC.tMin,
      baseCurveLC.tMax,
    ).filter(([a, b]) => baseCurveLC.at((a + b) / 2).length() < 1)
    const projectedCurves = [0, 1].map((id) => {
      return (t: number) => {
        const atSqr = snap(baseCurveLC.at(t).squared(), 1)
        const lineISTs = /* +- */ sqrt(1 - atSqr)
        //assert(!isNaN(lineISTs))
        return eq0(lineISTs)
          ? baseCurveLC.at(t)
          : baseCurveLC
              .at(t)
              .plus(surfaceLC.dir.times(sign(id - 0.5) * lineISTs))
      }
    })
    const dProjectedCurves = [0, 1].map((id) => {
      return (t: number) => {
        // d/dt sqrt(1 - baseCurveLC.at(t).squared())
        // = -1/2 * 1/sqrt(1 - baseCurveLC.at(t).squared()) * -2*baseCurveLC.at(t) * baseCurveLC.tangentAt(t)
        const atSqr = snap(baseCurveLC.at(t).squared(), 1)
        const lineISTs = /* +- */ baseCurveLC
          .at(t)
          .times(-1 / sqrt(1 - atSqr))
          .dot(baseCurveLC.tangentAt(t))
        //assert(!isNaN(lineISTs))
        return baseCurveLC
          .tangentAt(t)
          .plus(surfaceLC.dir.times(sign(id - 0.5) * lineISTs))
      }
    })
    //const f2 = t => sqrt(1 - baseCurveLC.at(t).squared())
    //const df2 = t => baseCurveLC.at(t).times(-1 / sqrt(1 -
    // baseCurveLC.at(t).squared())).dot(baseCurveLC.tangentAt(t)) checkDerivate(f2, df2, 0.31, 0.60)
    const curves = []
    for (const [aT, bT] of insideIntervals) {
      //const aLine = new L3(baseCurveLC.at(aT), surfaceLC.dir1)
      //const a = EllipsoidSurface.UNIT.isTsForLine(aLine).map(t => aLine.at(t))
      //const bLine = new L3(baseCurveLC.at(bT), surfaceLC.dir1)
      //const b = EllipsoidSurface.UNIT.isTsForLine(bLine).map(t => bLine.at(t))
      for (const i of [0, 1]) {
        const f = (t: number) => projectedCurves[i](t).y
        const df = (t: number) => dProjectedCurves[i](t).y
        checkDerivate(f, df, aT + 0.1, bT - 0.1)
        const tsAtY0 = getRoots(
          f,
          aT + NLA_PRECISION,
          bT - NLA_PRECISION,
          1 / (1 << 11),
          df,
        )
        const ii2 = getIntervals(tsAtY0, aT, bT).filter(
          ([a, b]) => f((a + b) / 2) > 0,
        )
        for (const [aT2, bT2] of ii2) {
          let aP = projectedCurves[i](aT2),
            bP = projectedCurves[i](bT2)
          0 === i && ([aP, bP] = [bP, aP])
          assert(EllipsoidSurface.UNIT.containsPoint(aP))
          assert(EllipsoidSurface.UNIT.containsPoint(bP))
          curves.push(
            PICurve.forStartEnd(
              surface,
              this,
              this.matrix.transformPoint(bP),
              this.matrix.transformPoint(aP),
              undefined,
            ),
          )
        }
      }
    }

    return surface.clipCurves(curves)
  }

  isCurvesWithSurface(surface: Surface): Curve[] {
    if (surface instanceof PlaneSurface) {
      return this.isCurvesWithPlane(surface.plane)
    } else if (surface instanceof CylinderSurface) {
      return this.isCurvesWithCylinderSurface(surface)
    } else if (surface instanceof EllipsoidSurface) {
      const surfaceLC = surface.transform(this.matrixInverse)
      const curves = EllipsoidSurface.unitISCurvesWithEllipsoidSurface(
        surfaceLC,
      ).map((c) => c.transform(this.matrix))
      return surface.clipCurves(curves)
    } else if (surface instanceof ProjectedCurveSurface) {
      return this.isCurvesWithPCS(surface)
    } else if (surface instanceof ParametricSurface) {
      let curves2 = ParametricSurface.isCurvesParametricImplicitSurface(
        surface,
        this,
        0.1,
        0.1,
        0.05,
      )
      curves2 = this.clipCurves(curves2)
      curves2 = surface.clipCurves(curves2)
      return curves2
    } else {
      throw new Error()
    }
  }

  isCurvesWithPlane(plane: P3) {
    const planeLC = plane.transform(this.matrixInverse)
    return EllipsoidSurface.unitISCurvesWithPlane(planeLC).map((c) =>
      c.transform(this.matrix),
    )
  }

  isCurvesWithCylinderSurface(surface: CylinderSurface): Curve[] {
    if (L3.containsPoint(surface.baseCurve.center, surface.dir, this.center)) {
      assert(this.isSphere())
      const ellipseProjected = surface.baseCurve.transform(
        M4.project(surface.baseCurve.getPlane(), surface.dir),
      )
      if (ellipseProjected.isCircular()) {
        const thisRadius = this.f1.length()
        const surfaceRadius = ellipseProjected.f1.length()
        // sphereRadius² = distanceISFromCenter² + isRadius²
        if (eq(thisRadius, surfaceRadius)) {
          // return
        } else if (surfaceRadius < thisRadius) {
        }
        assert(false)
      }
    }
    return this.isCurvesWithPCS(surface)
  }

  isTsForLine(line: L3) {
    assertInst(L3, line)
    // transforming line manually has advantage that dir1 will not be renormalized,
    // meaning that calculated values t for localLine are directly transferable to line
    const anchorLC = this.matrixInverse.transformPoint(line.anchor)
    const dirLC = this.matrixInverse.transformVector(line.dir1)
    return EllipsoidSurface.unitISTsWithLine(anchorLC, dirLC)
  }

  isCoplanarTo(surface: Surface): boolean {
    if (this === surface) return true
    if (!hasConstructor(surface, EllipsoidSurface)) return false
    if (!this.center.like(surface.center)) return false
    if (this.isSphere())
      return surface.isSphere() && eq(this.f1.length(), this.f2.length())

    const otherMatrixLC = this.matrixInverse.times(surface.matrix)
    // Ellipsoid with matrix otherMatrixLC is unit sphere iff otherMatrixLC is orthogonal
    return otherMatrixLC.like3x3() && otherMatrixLC.isOrthogonal()
  }

  containsEllipse(ellipse: EllipseCurve): boolean {
    const ellipseLC = ellipse.transform(this.matrixInverse)
    const distEllipseLCCenter = ellipseLC.center.length()
    const correctRadius = Math.sqrt(1 - distEllipseLCCenter ** 2)
    return (
      lt(distEllipseLCCenter, 1) &&
      ellipseLC.isCircular() &&
      ellipseLC.f1.hasLength(correctRadius)
    )
    //&& le(0, ellipseLC.getAABB().min.y)
  }

  containsCurve(curve: Curve): boolean {
    if (curve instanceof EllipseCurve) {
      return this.containsEllipse(curve)
    } else {
      return super.containsCurve(curve)
    }
  }

  transform(m4: M4): this {
    assert(m4.isNoProj(), () => m4.sce)
    return new EllipsoidSurface(
      m4.transformPoint(this.center),
      m4.transformVector(this.f1),
      m4.transformVector(this.f2),
      m4.transformVector(this.f3).times(m4.isMirroring() ? -1 : 1),
    ) as this
  }

  transform4(m4: M4): this {
    console.log("transform4")
    const resultMatrix = m4.times(this.matrix)
    console.log(resultMatrix.toString())
    const scaleDir = V(
      resultMatrix.m[12],
      resultMatrix.m[13],
      resultMatrix.m[14],
    )
    // need to find parameters where scaleDir is parallel to the normal
    const pLC = this.pLCNormalWCMatrix.inversed().transformPoint(scaleDir)
    const s = pLC.angleXY()
    const t = Math.asin(clamp(pLC.z, -1, 1))
    const fa = resultMatrix.transformPoint(scaleDir.unit())
    const fb = resultMatrix.transformPoint(scaleDir.unit().negated())
    const newCenter = V3.lerp(fa, fb, 0.5)
    console.log(scaleDir.sce, s, t, fa, fb, "newCenter", newCenter.sce)
    return new EllipsoidSurface(
      newCenter,
      m4.transformVector2(this.f1, this.center),
      m4.transformVector2(this.f2, this.center),
      m4
        .transformVector2(this.f3, this.center)
        .times(m4.isMirroring() ? -1 : 1),
    ) as this
  }

  isInsideOut(): boolean {
    return this.f1.cross(this.f2).dot(this.f3) < 0
  }

  flipped(): this {
    return new EllipsoidSurface(
      this.center,
      this.f1,
      this.f2,
      this.f3.negated(),
      this.uMin,
      this.uMax,
      -this.vMax,
      -this.vMin,
    ) as this
  }

  normalUVFunc(): (u: number, v: number) => V3 {
    // ugh
    // paramtric ellipsoid point q(a, b)
    // normal1 == (dq(a, b) / da) X (dq(a, b) / db) (cross product of partial derivatives)
    // normal1 == cos b * (f2 X f3 * cos b * cos a + f3 X f1 * cos b * sin a + f1 X f2 * sin b)
    return (a, b) => {
      const { f1, f2, f3 } = this
      const normal = f2
        .cross(f3)
        .times(Math.cos(b) * Math.cos(a))
        .plus(f3.cross(f1).times(Math.cos(b) * Math.sin(a)))
        .plus(f1.cross(f2).times(Math.sin(b)))
        //.times(Math.cos(b))
        .unit()
      return normal
    }
  }

  normalP(p: V3): V3 {
    return this.pLCNormalWCMatrix
      .transformVector(this.matrixInverse.transformPoint(p))
      .unit()
  }

  normalUV(u: number, v: number): V3 {
    return this.pLCNormalWCMatrix.transformVector(V3.sphere(u, v)).unit()
  }

  uvPFunc() {
    return (pWC: V3) => {
      const pLC = this.matrixInverse.transformPoint(pWC)
      const alpha = abs(pLC.angleXY())
      const beta = Math.asin(clamp(pLC.z, -1, 1))
      assert(isFinite(alpha))
      assert(isFinite(beta))
      return new V3(alpha, beta, 0)
    }
  }

  pUVFunc() {
    // this(a, b) = f1 cos a cos b + f2 sin a cos b + f2 sin b
    return (alpha: number, beta: number) => {
      return this.matrix.transformPoint(V3.sphere(alpha, beta))
    }
  }

  isSphere(): boolean {
    return (
      eq(this.f1.length(), this.f2.length()) &&
      eq(this.f2.length(), this.f3.length()) &&
      eq(this.f3.length(), this.f1.length()) &&
      this.f1.isPerpendicularTo(this.f2) &&
      this.f2.isPerpendicularTo(this.f3) &&
      this.f3.isPerpendicularTo(this.f1)
    )
  }

  isVerticalSpheroid(): boolean {
    return (
      eq(this.f1.length(), this.f2.length()) &&
      this.f1.isPerpendicularTo(this.f2) &&
      this.f2.isPerpendicularTo(this.f3) &&
      this.f3.isPerpendicularTo(this.f1)
    )
  }

  mainAxes(): EllipsoidSurface {
    // q(a, b) = f1 cos a cos b + f2 sin a cos b + f3 sin b
    // q(s, t, u) = s * f1 + t * f2 + u * f3 with s² + t² + u² = 1
    // (del q(a, b) / del a) = f1 (-sin a) cos b  + f2 cos a cos b
    // (del q(a, b) / del b) = f1 cos a (-sin b) + f2 sin a (-sin b) + f2 cos b
    // del q(s, t, u) / del a = -t f1 + s f2
    // (del q(a, b) / del a) DOT q(a, b) == 0
    // (f1 (-sin a) cos b  + f2 cos a cos b) DOT (f1 cos a cos b + f2 sin a cos b + f2 sin b) == 0
    // (del q(a, b) / del b) DOT q(a, b) == 0
    // (f1 cos a (-sin b) + f2 sin a (-sin b) + f2 cos b) DOT (f1 cos a cos b + f2 sin a cos b + f2 sin b) == 0

    // Solve[
    // (f1 (-sin a) cos b  + f2 cos a cos b) * (f1 cos a cos b + f2 sin a cos b + f2 sin b) = 0,
    // (f1 cos a (-sin b) + f2 sin a (-sin b) + f2 cos b) * (f1 cos a cos b + f2 sin a cos b + f2 sin b) = 0}, a, b]
    const { f1, f2, f3 } = this

    if (eq0(f1.dot(f2)) && eq0(f2.dot(f3)) && eq0(f3.dot(f1))) {
      return this
    }

    //const f = ([a, b], x?) => {
    //    const sinA = Math.sin(a), cosA = Math.cos(a), sinB = Math.sin(b), cosB = Math.cos(b)
    //    const centerToP = V3.add(f1.times(cosA * cosB), f2.times(sinA * cosB), f3.times(sinB))
    //    const centerToPdelA = f1.times(-sinA * cosB).plus(f2.times(cosA * cosB))
    //    const centerToPdelB = V3.add(f1.times(cosA * -sinB), f2.times(sinA * -sinB), f3.times(cosB))
    //    x && console.log(centerToP.sce, centerToPdelA.sce, centerToPdelB.sce)
    //    return [centerToP.dot(centerToPdelA), centerToP.dot(centerToPdelB)]
    //}
    //const mainF1Params = newtonIterate(f, [0, 0], 8), mainF1 = this.pUVFunc()(mainF1Params[0], mainF1Params[1])
    //console.log(f(mainF1Params, 1).sce)
    //const mainF2Params = newtonIterate(f, this.uvPFunc()(f2.rejectedFrom(mainF1)).toArray(2), 8),
    //   mainF2 = this.pUVFunc()(mainF2Params[0], mainF2Params[1])
    //console.log(this.normalUVFunc()(mainF2Params[0], mainF2Params[1]).sce)
    //assert(mainF1.isPerpendicularTo(mainF2), mainF1, mainF2, mainF1.dot(mainF2), mainF1Params)
    //const mainF3Params = this.uvPFunc()(mainF1.cross(mainF2)), mainF3 = this.pUVFunc()(mainF3Params[0],
    // mainF3Params[1]) return new EllipsoidSurface(this.center, mainF1, mainF2, mainF3)

    const { U, SIGMA } = this.matrix.svd3()
    assert(SIGMA.isDiagonal())
    assert(U.isOrthogonal())
    const U_SIGMA = U.times(SIGMA)
    // column vectors of U_SIGMA
    const [mainF1, mainF2, mainF3] = arrayFromFunction(
      3,
      (i) => new V3(U_SIGMA.m[i], U_SIGMA.m[i + 4], U_SIGMA.m[i + 8]),
    )
    return new EllipsoidSurface(this.center, mainF1, mainF2, mainF3)
  }

  containsPoint(p: V3) {
    return eq0(this.implicitFunction()(p))
  }

  boundsFunction() {
    return (a: number, b: number) => between(a, 0, PI) && between(b, -PI, PI)
  }

  volume(): number {
    return (4 / 3) * Math.PI * this.f1.dot(this.f2.cross(this.f3))
  }

  loopContainsPoint(loop: Edge[], pWC: V3): PointVsFace {
    if (!this.containsPoint(pWC)) return PointVsFace.OUTSIDE
    assertVectors(pWC)
    assert(Edge.isLoop(loop))
    const pLCXY = this.matrixInverse.transformPoint(pWC).xy()
    const testLine = new EllipseCurve(
      this.center,
      this.f3,
      pLCXY.likeO() ? this.f2 : this.matrix.transformVector(pLCXY.unit()),
    )

    if (P3.normalOnAnchor(this.f2.unit(), this.center).containsPoint(pWC)) {
      return loop.some(
        (edge) =>
          edge.curve.containsPoint(pWC) &&
          fuzzyBetween(edge.curve.pointT(pWC), edge.minT, edge.maxT),
      )
        ? PointVsFace.ON_EDGE
        : PointVsFace.OUTSIDE
    }

    return Surface.loopContainsPointEllipse(loop, pWC, testLine)
  }

  surfaceAreaApprox(): number {
    // See https://en.wikipedia.org/wiki/Ellipsoid#Surface_area
    const mainAxes = this.mainAxes(),
      a = mainAxes.f1.length(),
      b = mainAxes.f2.length(),
      c = mainAxes.f3.length()
    const p = 1.6075
    return (
      4 *
      PI *
      Math.pow(
        (Math.pow(a * b, p) + Math.pow(b * c, p) + Math.pow(c * a, p)) / 3,
        1 / p,
      )
    )
  }

  surfaceArea(): number {
    // See https://en.wikipedia.org/wiki/Ellipsoid#Surface_area
    const mainAxes = this.mainAxes(),
      f1l = mainAxes.f1.length(),
      f2l = mainAxes.f2.length(),
      f3l = mainAxes.f3.length(),
      [c, b, a] = [f1l, f2l, f3l].sort(MINUS)

    // https://en.wikipedia.org/w/index.php?title=Spheroid&oldid=761246800#Area
    function spheroidArea(a: number, c: number) {
      if (c < a) {
        const eccentricity2 = 1 - c ** 2 / a ** 2
        const eccentricity = Math.sqrt(eccentricity2)
        return (
          2 *
          PI *
          a ** 2 *
          (1 +
            ((1 - eccentricity2) / Math.sqrt(eccentricity)) *
              Math.atanh(eccentricity))
        )
      } else {
        const eccentricity = Math.sqrt(1 - a ** 2 / c ** 2)
        return (
          2 *
          PI *
          a ** 2 *
          (1 + (c / a / eccentricity) * Math.asin(eccentricity))
        )
      }
    }

    if (eq(a, b)) {
      return spheroidArea(a, c)
    } else if (eq(b, c)) {
      return spheroidArea(b, a)
    } else if (eq(c, a)) {
      return spheroidArea(c, b)
    }

    const phi = Math.acos(c / a)
    const kk = (a ** 2 * (b ** 2 - c ** 2)) / (b ** 2 * (a ** 2 - c ** 2))
    const incompleteEllipticInt1 = gaussLegendreQuadrature24(
      (phi) => Math.pow(1 - kk * Math.sin(phi) ** 2, -0.5),
      0,
      phi,
    )
    const incompleteEllipticInt2 = gaussLegendreQuadrature24(
      (phi) => Math.pow(1 - kk * Math.sin(phi) ** 2, 0.5),
      0,
      phi,
    )
    return (
      (2 * PI * c ** 2 + (2 * PI * a * b) / Math.sin(phi)) *
      (incompleteEllipticInt2 * Math.sin(phi) ** 2 +
        incompleteEllipticInt1 * Math.cos(phi) ** 2)
    )
  }

  getSeamPlane(): P3 {
    const plane = P3.forAnchorAndPlaneVectors(this.center, this.f1, this.f3)
    return plane.normal1.dot(this.f2) < 0 ? plane : plane.flipped()
  }

  getExtremePoints(): V3[] {
    return getExtremePointsHelper.call(
      this,
      new EllipseCurve(V3.O, V3.X, V3.Z, -PI / 2, PI / 2),
    )
  }

  pointFoot(pWC: V3, startS?: number, startT?: number): V3 {
    console.log(pWC.sce)
    if (undefined === startS || undefined === startT) {
      let pLC1 = this.matrixInverse.transformPoint(pWC).unit()
      if (pLC1.y < 0) pLC1 = pLC1.negated()
      ;({ x: startS, y: startT } = EllipsoidSurface.UNIT.uvP(pLC1))
    }
    const dpdu = this.dpdu()
    const dpdv = this.dpdv()
    const [u, v] = newtonIterate(
      ([u, v]) => {
        const p = this.pUV(u, v)
        console.log(
          [p, p.plus(dpdu(u, v)), p, p.plus(dpdv(u, v))].map(toSource).join() +
            ",",
        )
        const pUVToPWC = this.pUV(u, v).to(pWC)
        return [pUVToPWC.dot(dpdu(u, v)), pUVToPWC.dot(dpdv(u, v))]
      },
      [startS, startT],
      8,
      undefined,
      0.1,
    )
    return new V3(u, v, 0)
  }

  implicitFunction() {
    return (pWC: V3) => {
      const pLC = this.matrixInverse.transformPoint(pWC)
      return (pLC.length() - 1) * this.normalDir
    }
  }

  // = this.inverseMatrix.transformPoint(this.inverseMatrix.transformPoint(pWC).unit())
  didp(pWC: V3) {
    // i(pWC) = this.inverseMatrix.transformPoint(pWC).length() - 1
    // chain diff rule
    const pLC = this.matrixInverse.transformPoint(pWC)
    return this.pLCNormalWCMatrix.transformVector(pLC.unit()) //.times(this.normalDir)
  }

  /*+
   * An ellipsoid remains an ellipsoid after a perspective transform (as long as it does not intersect the vanishing
   * plane. This transforms a matrix with a perspective component into one which would return an identical ellipsoid,
   * but with no perspective component.
   */
  static unitTransform4(m: M4): M4 {
    m.m[15] !== 1 && (m = m.divScalar(m.m[15]))
    // X * P = m => X = m * P^-1
    // prettier-ignore
    const Pinv = new M4(
            1,        0,        0, 0,
            0,        1,        0, 0,
            0,        0,        1, 0,
            -m.m[12], -m.m[13], -m.m[14], 1,
        )
    const pn = new V3(m.m[12], m.m[13], m.m[14]),
      pw = m.m[15]
    const pwSqrMinusPnSqr = pw ** 2 - pn.squared()
    if (lt(pwSqrMinusPnSqr, 0)) {
      throw new Error("vanishing plane intersects unit sphere")
    }
    const c = pn.div(-pwSqrMinusPnSqr)
    const scale = pn.times(
      (pw * pn.length()) / (pn.squared() * -pwSqrMinusPnSqr),
    )
    const scale1 = pw / -pwSqrMinusPnSqr
    const scale2 = 1 / sqrt(pwSqrMinusPnSqr)
    const rotNX = M4.forSys(pn.unit(), pn.getPerpendicular().unit())
    return M4.product(
      m,
      Pinv,
      M4.translate(c),
      rotNX,
      M4.scale(scale1, scale2, scale2),
      rotNX.transposed(),
    )
  }
}
EllipsoidSurface.prototype.uStep = PI / 32
EllipsoidSurface.prototype.vStep = PI / 32
