import { SVGPathData } from "svg-pathdata"
import {
  arrayFromFunction,
  arrayRange,
  assert,
  DEG,
  eq0,
  getIntervals,
  int,
  M4,
  MINUS,
  mod,
  newtonIterate,
  NLA_PRECISION,
  PI,
  TAU,
  V,
  V3,
} from "ts3dutils"

import {
  BezierCurve,
  Curve,
  Edge,
  EllipseCurve,
  L3,
  ParabolaCurve,
  PCurveEdge,
  StraightEdge,
} from "./index"
import { floor, ceil, abs, sign } from "./math"

export function edgePathFromSVG(pathString: string): Edge[] {
  let currentPos: V3 = undefined!
  const parsed: any[] = new SVGPathData(pathString)
    .toAbs()
    .normalizeHVZ()
    .sanitize(NLA_PRECISION)
    .annotateArcs().commands
  const path: Edge[] = []
  for (const c of parsed) {
    assert("x" in c && "y" in c)
    const endPos = new V3(c.x, c.y, 0)
    switch (c.type) {
      case SVGPathData.LINE_TO:
        path.push(StraightEdge.throughPoints(currentPos, endPos))
        break
      case SVGPathData.CURVE_TO: {
        const c1 = new V3(c.x1, c.y1, 0)
        const c2 = new V3(c.x2, c.y2, 0)
        const curve = new BezierCurve(currentPos, c1, c2, endPos, 0, 1)
        const edge = new PCurveEdge(
          curve,
          currentPos,
          endPos,
          0,
          1,
          undefined,
          curve.tangentAt(0),
          curve.tangentAt(1),
        )
        path.push(edge)
        break
      }
      case SVGPathData.QUAD_TO: {
        const c1 = new V3(c.x1, c.y1, 0)
        const curve = ParabolaCurve.quadratic(
          currentPos,
          c1,
          endPos,
        ).rightAngled()
        const edge = new PCurveEdge(
          curve,
          currentPos,
          endPos,
          curve.tMin,
          curve.tMax,
          undefined,
          curve.tangentAt(curve.tMin),
          curve.tangentAt(curve.tMax),
        )
        path.push(edge)
        break
      }
      case SVGPathData.ARC: {
        const phi1 = c.phi1 * DEG,
          phi2 = c.phi2 * DEG,
          [phiMin, phiMax] = [phi1, phi2].sort(MINUS)
        const stops = arrayRange(-3, 4, 1)
          .map((n) => n * PI)
          .filter((stop) => phiMin <= stop && stop <= phiMax)
        const center = V(c.cX, c.cY)
        const f1 = V3.polar(c.rX, c.xRot * DEG)
        const f2 = V3.polar(c.rY, c.xRot * DEG + Math.PI / 2)
        const edges = getIntervals(stops, phiMin, phiMax).map(([t1, t2]) => {
          const deltaT = t2 - t1
          const t1_ = mod(t1, TAU)
          const t2_ = t1_ + deltaT
          assert(t1_ >= 0 == t2_ >= 0)
          const gtPI = t1_ > PI || t2_ > PI
          const aT = gtPI ? t1_ - PI : t1_
          const bT = gtPI ? t2_ - PI : t2_
          const curve = new EllipseCurve(
            center,
            gtPI ? f1.negated() : f1,
            gtPI ? f2.negated() : f2,
          )
          const a = phi1 == t1 ? currentPos : phi2 == t1 ? endPos : curve.at(aT)
          const b = phi1 == t2 ? currentPos : phi2 == t2 ? endPos : curve.at(bT)
          return new PCurveEdge(
            curve,
            a,
            b,
            aT,
            bT,
            undefined,
            curve.tangentAt(aT),
            curve.tangentAt(bT),
          )
        })
        path.push(...(c.phiDelta > 0 ? edges : Edge.reversePath(edges)))
        break
      }
    }
    currentPos = endPos
  }
  return path
}

/**
 * Create an axis-aligned rectangle of edges on the XY-plane with the bottom-left corner on the origin.
 * @param width
 * @param height
 */
export function edgeRect(width: number = 1, height: number = width): Edge[] {
  const vertices = [
    new V3(0, 0, 0),
    new V3(width, 0, 0),
    new V3(width, height, 0),
    new V3(0, height, 0),
  ]
  return StraightEdge.chain(vertices)
}

export function ngon(n: int = 3, radius: number = 1): Edge[] {
  return StraightEdge.chain(
    arrayFromFunction(n, (i) => V3.polar(radius, (TAU * i) / n)),
  )
}

export function star(
  pointCount: int = 5,
  r0: number = 1,
  r1: number = 0.5,
): Edge[] {
  const vertices = arrayFromFunction(pointCount * 2, (i) =>
    V3.polar(0 == i % 2 ? r0 : r1, (TAU * i) / pointCount / 2),
  )
  return StraightEdge.chain(vertices)
}

export function createEdge(
  curve: Curve,
  a: V3,
  b: V3,
  aT: number,
  bT: number,
  flippedOf: Edge | undefined,
  aDir: V3,
  bDir: V3,
  name?: string,
): Edge {
  if (curve instanceof L3) {
    return new StraightEdge(
      curve,
      a,
      b,
      aT,
      bT,
      flippedOf as StraightEdge,
      name,
    )
  } else {
    return new PCurveEdge(
      curve,
      a,
      b,
      aT,
      bT,
      flippedOf as PCurveEdge,
      aDir,
      bDir,
      name,
    )
  }
}

export function edgeForCurveAndTs(
  curve: Curve,
  aT: number = curve.tMin,
  bT: number = curve.tMax,
): Edge {
  return createEdge(
    curve,
    curve.at(aT),
    curve.at(bT),
    aT,
    bT,
    undefined,
    aT < bT ? curve.tangentAt(aT) : curve.tangentAt(aT).negated(),
    aT < bT ? curve.tangentAt(bT) : curve.tangentAt(bT).negated(),
  )
}

export function reuleaux(n: int = 3, radius: number = 1): Edge[] {
  assert(3 <= n)
  assert(1 == n % 2)
  const corners = arrayFromFunction(n, (i) => V3.polar(radius, (TAU * i) / n))
  return arrayFromFunction(n, (i) => {
    const aI = (i + floor(n / 2)) % n,
      bI = (i + ceil(n / 2)) % n
    const a = corners[aI],
      b = corners[bI]
    const center = corners[i]
    const f1 = center.to(a),
      curve = new EllipseCurve(center, f1, V3.Z.cross(f1))
    return createEdge(
      curve,
      a,
      b,
      0,
      curve.pointT(b),
      undefined,
      V3.Z.cross(f1),
      V3.Z.cross(center.to(b)),
    )
  })
}

export function round(edges: Edge[], radius: number) {
  if (eq0(radius)) {
    return edges
  }
  const corners = edges.map((edge, i) => {
    const j = (i + 1) % edges.length,
      nextEdge = edges[j]
    if (!edge.b.like(nextEdge.a)) return undefined
    const angleToNext = edge.bDir.angleTo(nextEdge.aDir)
    const c1 = edge.curve,
      c2 = nextEdge.curve
    if (c1 instanceof L3 && c2 instanceof L3) {
      const normal = c1.dir1.cross(c2.dir1)
      if (eq0(angleToNext)) return undefined

      const l1inside = normal.cross(c1.dir1),
        l2inside = normal.cross(c2.dir1)
      const l1offset = c1.transform(M4.translate(l1inside.toLength(radius)))
      const l2offset = c2.transform(M4.translate(l2inside.toLength(radius)))
      const center = l1offset.isInfoWithLine(l2offset)
      if (!center) throw new Error("tangential curves")
      const cornerA = center.plus(l1inside.toLength(-radius))
      const cornerB = center.plus(l2inside.toLength(-radius))
      const f1 = l1inside.toLength(-radius)
      const curve = new EllipseCurve(
        center,
        f1,
        normal.cross(f1).toLength(radius),
      )
      const cornerEdge = createEdge(
        curve,
        cornerA,
        cornerB,
        0,
        curve.pointT(cornerB),
        undefined,
        c1.dir1,
        c2.dir1,
      )
      return cornerEdge
    } else {
      return arbitraryCorner(edge, nextEdge, radius)
    }
  })
  const result = edges.flatMap((edge, i) => {
    const h = (i + edges.length - 1) % edges.length
    const prevCorner = corners[h],
      nextCorner = corners[i]
    if (!prevCorner && !nextCorner) {
      return edge
    }
    const [aT, a, aDir] = !prevCorner
      ? [edge.aT, edge.a, edge.aDir]
      : [edge.curve.pointT(prevCorner.b), prevCorner.b, prevCorner.bDir]
    const [bT, b, bDir] = !nextCorner
      ? [edge.bT, edge.b, edge.bDir]
      : [edge.curve.pointT(nextCorner.a), nextCorner.a, nextCorner.aDir]
    const newEdge = createEdge(edge.curve, a, b, aT, bT, undefined, aDir, bDir)
    return !nextCorner ? newEdge : [newEdge, nextCorner]
  })
  return result
}

export function arbitraryCorner(e1: Edge, e2: Edge, radius: number) {
  const c1 = e1.curve,
    c2 = e2.curve

  function f([t1, t2]: number[]) {
    const p1 = c1.at(t1),
      p2 = c2.at(t2)
    const dp1 = c1.tangentAt(t1),
      dp2 = c2.tangentAt(t2)
    const virtualPlaneNormal = dp1.cross(dp2)
    const normal1 = virtualPlaneNormal.cross(dp1).unit(),
      normal2 = virtualPlaneNormal.cross(dp2).unit()
    const dirCross = normal1.cross(normal2)
    if (virtualPlaneNormal.likeO()) {
      assert(false)
    } // lines parallel
    const p1p2 = p1.to(p2)
    // check if distance is zero (see also L3.distanceToLine)
    if (!eq0(p1p2.dot(virtualPlaneNormal))) {
      assert(false)
    }
    const dist1 = p1p2.cross(normal2).dot(dirCross) / dirCross.squared()
    const dist2 = p1p2.cross(normal1).dot(dirCross) / dirCross.squared()
    const g1 = p1.plus(normal1.times(dist1))
    const g2 = p2.plus(normal2.times(dist2))
    assert(g1.like(g2))
    return [abs(dist1) - radius, abs(dist2) - radius]
  }

  const startT1 = e1.bT - (radius * sign(e1.deltaT())) / e1.bDir.length()
  const startT2 = e2.aT + (radius * sign(e2.deltaT())) / e2.aDir.length()
  const [t1, t2] = newtonIterate(f, [startT1, startT2])
  const cornerA = e1.curve.at(t1)
  const cornerB = e2.curve.at(t2)
  const dp1 = c1.tangentAt(t1),
    dp2 = c2.tangentAt(t2)
  const virtualPlaneNormal = dp1.cross(dp2)
  const normal1 = virtualPlaneNormal.cross(dp1).unit()
  const f1 = normal1.toLength(-radius)
  const center = cornerA.minus(f1)
  const curve = new EllipseCurve(
    center,
    f1,
    virtualPlaneNormal.cross(f1).toLength(radius),
  )
  const cornerEdge = createEdge(
    curve,
    cornerA,
    cornerB,
    0,
    curve.pointT(cornerB),
    undefined,
    c1.tangentAt(t1),
    c2.tangentAt(t2),
  )
  return cornerEdge
}
