import {
  assert,
  assertf,
  assertVectors,
  callSource,
  int,
  M4,
  snap2,
  V3,
  snap,
} from "ts3dutils"

import { Curve, Edge, P3, PICurve, Surface } from "."

export class PCurveEdge extends Edge {
  constructor(
    curve: Curve,
    a: V3,
    b: V3,
    aT: number,
    bT: number,
    public flippedOf: PCurveEdge | undefined,
    readonly aDir: V3,
    readonly bDir: V3,
    name?: string,
  ) {
    super(curve, a, b, aT, bT, flippedOf, name)
    assertVectors(aDir, bDir)
    assertf(() => !aDir.likeO(), curve)
    assertf(() => !bDir.likeO(), curve)
    if (!(curve instanceof PICurve)) {
      // TODO
      assertf(
        () => curve.tangentAt(aT).likeOrReversed(aDir),
        "" + aT + curve.tangentAt(aT).sce + " " + aDir.sce,
      )
      assertf(
        () => curve.tangentAt(bT).likeOrReversed(bDir),
        "" + bT + curve.tangentAt(bT).sce + " " + bDir.sce,
      )
    }
    assert(
      this.reversed === this.aDir.dot(curve.tangentAt(aT)) < 0,
      aT +
        " " +
        bT +
        " " +
        curve.constructor.name +
        " " +
        this.aDir.sce +
        " " +
        this.bDir.sce +
        " " +
        curve.tangentAt(aT),
    )
    assert(
      this.reversed === this.bDir.dot(curve.tangentAt(bT)) < 0,
      aT +
        " " +
        bT +
        " " +
        curve.constructor.name +
        " " +
        this.aDir.sce +
        " " +
        this.bDir.sce +
        " " +
        curve.tangentAt(aT),
    )
  }

  static forCurveAndTs(curve: Curve, aT: number, bT: number, name?: string) {
    return new PCurveEdge(
      curve,
      curve.at(aT),
      curve.at(bT),
      aT,
      bT,
      undefined,
      aT < bT ? curve.tangentAt(aT) : curve.tangentAt(aT).negated(),
      aT < bT ? curve.tangentAt(bT) : curve.tangentAt(bT).negated(),
      name,
    )
  }

  toSource(): string {
    return callSource(
      "new PCurveEdge",
      this.curve,
      this.a,
      this.b,
      this.aT,
      this.bT,
      undefined,
      this.aDir,
      this.bDir,
      this.name,
    )
  }

  getVerticesNo0(): V3[] {
    return this.curve.calcSegmentPoints(
      this.aT,
      this.bT,
      this.a,
      this.b,
      this.reversed,
      false,
    )
  }

  pointsCount(): int {
    return this.points().length
  }

  points(): V3[] {
    return this.curve.calcSegmentPoints(
      this.aT,
      this.bT,
      this.a,
      this.b,
      this.reversed,
      true,
    )
  }

  edgeISTsWithSurface(surface: Surface): number[] {
    return this.curve
      .isTsWithSurface(surface)
      .map((edgeT) => snap2(edgeT, this.aT, this.bT))
      .filter((edgeT) => this.minT <= edgeT && edgeT <= this.maxT)
  }

  edgeISTsWithPlane(surface: P3): number[] {
    return this.curve
      .isTsWithPlane(surface)
      .map((edgeT) => snap2(edgeT, this.aT, this.bT))
      .filter((edgeT) => this.minT <= edgeT && edgeT <= this.maxT)
  }

  tangentAt(t: number): V3 {
    return !this.reversed
      ? this.curve.tangentAt(t)
      : this.curve.tangentAt(t).negated()
  }

  flipped(): PCurveEdge {
    return (
      this.flippedOf ||
      (this.flippedOf = new PCurveEdge(
        this.curve,
        this.b,
        this.a,
        this.bT,
        this.aT,
        this,
        this.bDir.negated(),
        this.aDir.negated(),
        this.name,
      ))
    )
  }

  transform(m4: M4, desc?: string): this {
    return new PCurveEdge(
      this.curve.transform(m4),
      m4.transformPoint(this.a),
      m4.transformPoint(this.b),
      this.aT,
      this.bT,
      undefined,
      m4.transformVector(this.aDir),
      m4.transformVector(this.bDir),
      "" + this.name + desc,
    ) as this
  }

  transform4(m4: M4, desc?: string): this {
    const a_ = m4.transformPoint(this.a)
    const b_ = m4.transformPoint(this.b)
    const curve_ = this.curve.transform4(m4)
    return new PCurveEdge(
      curve_,
      a_,
      b_,
      snap(curve_.pointT(a_), this.aT),
      snap(curve_.pointT(b_), this.bT),
      undefined,
      m4.transformVector(this.aDir),
      m4.transformVector(this.bDir),
      "" + this.name + desc,
    ) as this
  }

  isCoEdge(edge: Edge): boolean {
    return (
      this === edge ||
      this === edge.flippedOf ||
      (this.curve.isColinearTo(edge.curve) &&
        ((this.a.like(edge.a) && this.b.like(edge.b)) ||
          (this.a.like(edge.b) && this.b.like(edge.a))))
    )
  }
  split(t: number): [Edge, Edge] {
    const p = this.curve.at(t)
    const dir = this.tangentAt(t)
    return [
      new PCurveEdge(
        this.curve,
        this.a,
        p,
        this.aT,
        t,
        undefined,
        this.aDir,
        dir,
        this.name + "left",
      ),
      new PCurveEdge(
        this.curve,
        p,
        this.b,
        t,
        this.bT,
        undefined,
        this.aDir,
        dir,
        this.name + "right",
      ),
    ]
  }
}
