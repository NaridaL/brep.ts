import {
  AABB,
  arrayFromFunction,
  arrayRange,
  assert,
  assertf,
  assertInst,
  assertNumbers,
  assertVectors,
  callsce,
  clamp,
  eq,
  fuzzyBetween,
  int,
  le,
  lt,
  Transformable,
  Tuple3,
  V,
  V3,
} from "ts3dutils"

import { Curve, L3, P3, Surface } from "."
import { sign } from "./math"

export abstract class Edge extends Transformable {
  readonly reversed: boolean

  abstract get aDir(): V3
  abstract get bDir(): V3

  constructor(
    readonly curve: Curve,
    readonly a: V3,
    readonly b: V3,
    readonly aT: number,
    readonly bT: number,
    public flippedOf?: Edge | undefined,
    readonly name?: string,
  ) {
    super()
    assertNumbers(aT, bT)
    assert(!eq(aT, bT))
    assertVectors(a, b)
    assertf(() => curve instanceof Curve, curve)
    assertf(
      () => !curve.isValidT || (curve.isValidT(aT) && curve.isValidT(bT)),
      aT,
      bT,
      curve,
    )
    //if (curve instanceof PICurve) {
    //    assertf(() => curve.at(aT).to(a).length() < 0.1, ''+curve.at(aT)+a)
    //    assertf(() => curve.at(bT).to(b).length() < 0.1, '' + curve.at(bT) + b)
    //} else {
    assertf(
      () => curve.at(aT).like(a),
      () => "" + curve.at(aT) + a + " aT should have been " + curve.pointT(a),
    )
    assertf(
      () => curve.at(bT).like(b),
      () => "" + curve.at(bT) + b + " bT should have been " + curve.pointT(b),
    )
    //}
    assertf(
      () => fuzzyBetween(aT, curve.tMin, curve.tMax),
      aT,
      curve.tMin,
      curve.tMax,
    )
    assertf(
      () => fuzzyBetween(bT, curve.tMin, curve.tMax),
      bT,
      curve.tMin,
      curve.tMax,
    )
    assert(!a.like(b), "!a.like(b)" + a + b)
    this.aT = clamp(aT, curve.tMin, curve.tMax)
    this.bT = clamp(bT, curve.tMin, curve.tMax)
    this.reversed = this.aT > this.bT
  }

  get minT() {
    return Math.min(this.aT, this.bT)
  }

  get maxT() {
    return Math.max(this.aT, this.bT)
  }

  static isLoop(loop: Edge[]): boolean {
    return loop.every((edge, i) => edge.b.like(loop[(i + 1) % loop.length].a))
  }

  static edgesIntersect(e1: Edge, e2: Edge) {
    // TODO: still getting some NaNs here..
    assertNumbers(e1.curve.hlol, e2.curve.hlol)
    assertInst(Edge, e1, e2)
    if (e1.curve.hlol < e2.curve.hlol) {
      ;[e2, e1] = [e1, e2]
    }
    const sts = e1.curve.isInfosWithCurve(e2.curve)
    if (sts.some((info) => isNaN(info.tThis) || isNaN(info.tOther))) {
      console.log(e1.sce)
      console.log(e2.sce)
      assert(false)
    }
    return sts.some(
      /// (  e1.aT < tThis < e1.bT  )  &&  (  e2.aT < tOther < e2.bT  )
      ({ tThis, tOther }) => {
        return e1.tValueInside(tThis) && e2.tValueInside(tOther)
      },
    )
  }

  static assertLoop(edges: Edge[]): void {
    edges.forEach((edge, i) => {
      const j = (i + 1) % edges.length
      assert(
        edge.b.like(edges[j].a),
        `edges[${i}].b != edges[${j}].a (${edges[i].b.sce} != ${edges[j].a.sce})`,
      )
    })
  }

  static reversePath(path: Edge[], doReverse: boolean = true): Edge[] {
    return doReverse
      ? arrayFromFunction(path.length, (i) =>
          path[path.length - 1 - i].flipped(),
        )
      : path
  }

  abstract tangentAt(t: number): V3

  toString(): string {
    return callsce(
      "new " + this.constructor.name,
      this.curve,
      this.a,
      this.b,
      this.aT,
      this.bT,
      undefined,
      this.aDir,
      this.bDir,
    )
  }

  abstract edgeISTsWithSurface(surface: Surface): number[]

  /**
   * Returns the intersections of the edge with the plane.
   * Values are snapped to aT and bT, ie aT === t || !eq(aT, t)
   */
  abstract edgeISTsWithPlane(plane: P3): number[]

  colinearToLine(line: L3): boolean {
    return this.curve instanceof L3 && this.curve.isColinearTo(line)
  }

  tValueInside(t: number) {
    return this.aT < this.bT
      ? lt(this.aT, t) && lt(t, this.bT)
      : lt(this.bT, t) && lt(t, this.aT)
  }

  isValidT(t: number): boolean {
    return this.aT < this.bT
      ? le(this.aT, t) && le(t, this.bT)
      : le(this.bT, t) && le(t, this.aT)
  }

  clampedT(t: number): number {
    return this.aT < this.bT
      ? clamp(t, this.aT, this.bT)
      : clamp(t, this.bT, this.aT)
  }

  abstract flipped(): Edge

  /**
   * this is equals-equals. "isColinearTo" might make more sense but can't be used, because you can't get a
   * consistent hashCode for colinear curves
   * @param obj
   * @returns
   */
  equals(obj: any): boolean {
    return (
      this === obj ||
      (this.constructor == obj.constructor &&
        this.a.equals(obj.a) &&
        this.b.equals(obj.b) &&
        this.curve.equals(obj.curve))
    )
  }

  hashCode(): int {
    let hashCode = 0
    hashCode = hashCode * 31 + this.a.hashCode()
    hashCode = hashCode * 31 + this.b.hashCode()
    hashCode = hashCode * 31 + this.curve.hashCode()
    return hashCode | 0
  }

  like(edge: Edge) {
    // TODO this breaks on colinear edges,
    // TODO: what, where?
    return (
      this === edge ||
      (edge instanceof Edge &&
        this.curve.isColinearTo(edge.curve) &&
        this.a.like(edge.a) &&
        this.b.like(edge.b))
    )
  }

  /**
   * Get edge points, excluding the first one.
   */
  abstract getVerticesNo0(): V3[]

  abstract pointsCount(): int

  isCanon() {
    return !this.reversed
  }

  getCanon() {
    return this.reversed ? this.flipped() : this
  }

  overlaps(edge: Edge, noback?: boolean): boolean {
    assert(this.curve.isColinearTo(edge.curve))
    const edgeAT = this.curve.containsPoint(edge.a) && this.curve.pointT(edge.a)
    const edgeBT = this.curve.containsPoint(edge.b) && this.curve.pointT(edge.b)
    if (false === edgeAT && false === edgeBT) {
      return noback ? false : edge.overlaps(this, true)
    }
    return !(le(edge.maxT, this.minT) || le(this.maxT, edge.minT))
  }

  getAABB(): AABB {
    const min: Tuple3<number> = [Infinity, Infinity, Infinity],
      max: Tuple3<number> = [-Infinity, -Infinity, -Infinity]
    this.curve.roots().forEach((ts, dim) => {
      ts.forEach((t) => {
        if (lt(this.minT, t) && lt(t, this.maxT)) {
          min[dim] = Math.min(min[dim], this.curve.at(t).e(dim))
          max[dim] = Math.max(max[dim], this.curve.at(t).e(dim))
        }
      })
    })
    const aabb = new AABB(V(min), V(max))
    aabb.addPoint(this.a)
    aabb.addPoint(this.b)
    return aabb
  }

  length(steps: int = 1): number {
    return this.curve.arcLength(this.minT, this.maxT, steps)
  }

  abstract isCoEdge(other: Edge): boolean

  abstract points(): V3[]

  abstract split(t: number): [Edge, Edge]

  deltaT() {
    return this.bT - this.aT
  }

  deltaTSign() {
    return sign(this.bT - this.aT) as -1 | 1
  }

  atAvgT() {
    return this.curve.at((this.minT + this.maxT) / 2)
  }

  /**
   * Whether two edge loops are equal. Takes into account that two loops need not start with the same edge.
   * @param loop1
   * @param loop2
   */
  static loopsEqual(loop1: Edge[], loop2: Edge[]): boolean {
    return (
      loop1.length == loop2.length &&
      arrayRange(0, loop1.length, 1).some((offset) =>
        loop1.every((edge, i) =>
          edge.equals(loop2[(offset + i) % loop1.length]),
        ),
      )
    )
  }
}
