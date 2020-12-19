import {
  arrayFromFunction,
  assertf,
  assertInst,
  assertVectors,
  callsce,
  eq0,
  M4,
  snap2,
  V3,
} from "ts3dutils"

import { Edge, L3, P3, PlaneSurface, Surface } from "."

export class StraightEdge extends Edge {
  readonly tangent: V3
  readonly curve!: L3

  constructor(
    line: L3,
    a: V3,
    b: V3,
    aT: number,
    bT: number,
    public flippedOf?: StraightEdge,
    name?: string,
  ) {
    super(line, a, b, aT, bT, flippedOf, name)
    assertInst(L3, line)
    !flippedOf || assertInst(StraightEdge, flippedOf)
    !name || assertf(() => "string" === typeof name, name)
    this.tangent =
      this.aT < this.bT ? this.curve.dir1 : this.curve.dir1.negated()
  }

  get aDir() {
    return this.tangent
  }

  get bDir() {
    return this.tangent
  }

  static throughPoints(a: V3, b: V3, name?: string) {
    return new StraightEdge(
      L3.throughPoints(a, b, 0, a.to(b).length()),
      a,
      b,
      0,
      a.to(b).length(),
      undefined,
      name,
    )
  }

  /**
   * Create a list of StraightEdges from a list of vertices.
   * @param vertices
   * @param closed Whether to connect the first and last vertices. Defaults to true.
   * @returns
   */
  static chain(
    vertices: ReadonlyArray<V3>,
    closed: boolean = true,
  ): StraightEdge[] {
    const vc = vertices.length
    return arrayFromFunction(closed ? vc : vc - 1, (i) =>
      StraightEdge.throughPoints(vertices[i], vertices[(i + 1) % vc]),
    )
  }

  toSource(): string {
    return callsce(
      "new StraightEdge",
      this.curve,
      this.a,
      this.b,
      this.aT,
      this.bT,
    )
  }

  getVerticesNo0() {
    return [this.b]
  }

  pointsCount() {
    return 2
  }

  points() {
    return [this.a, this.b]
  }

  edgeISTsWithPlane(plane: P3): number[] {
    const edgeT = snap2(this.curve.isTWithPlane(plane), this.aT, this.bT)
    return this.minT <= edgeT && edgeT <= this.maxT ? [edgeT] : []
  }

  edgeISTsWithSurface(surface: Surface): number[] {
    if (surface instanceof PlaneSurface) {
      return this.edgeISTsWithPlane(surface.plane)
    } else {
      return surface
        .isTsForLine(this.curve)
        .map((edgeT) => snap2(edgeT, this.aT, this.bT))
        .filter((edgeT) => this.minT <= edgeT && edgeT <= this.maxT)
    }
  }

  tangentAt() {
    return this.tangent
  }

  flipped(): StraightEdge {
    return (
      this.flippedOf ||
      (this.flippedOf = new StraightEdge(
        this.curve,
        this.b,
        this.a,
        this.bT,
        this.aT,
        this,
        this.name,
      ))
    )
  }

  transform(m4: M4, desc?: string): this {
    const lineDir1TransLength = m4
      .transformVector2(this.curve.dir1, this.curve.anchor)
      .length()
    const curve = this.curve.transform(m4)
    const a = m4.transformPoint(this.a)
    const b = m4.transformPoint(this.b)
    return new StraightEdge(
      curve,
      a,
      b,
      m4.isNoProj() ? this.aT * lineDir1TransLength : curve.pointT(a),
      m4.isNoProj() ? this.bT * lineDir1TransLength : curve.pointT(b),
      undefined,
      "" + this.name + desc,
    ) as this
  }

  transform4(m4: M4, desc?: string): this {
    const lineDir1TransLength = m4
      .transformVector2(this.curve.dir1, this.curve.anchor)
      .length()
    const curve = this.curve.transform4(m4)
    const a = m4.transformPoint(this.a)
    const b = m4.transformPoint(this.b)
    return new StraightEdge(
      curve,
      a,
      b,
      m4.isNoProj() ? this.aT * lineDir1TransLength : curve.pointT(a),
      m4.isNoProj() ? this.bT * lineDir1TransLength : curve.pointT(b),
      undefined,
      "" + this.name + desc,
    ) as this
  }

  isCoEdge(edge: Edge): boolean {
    return (
      this === edge ||
      this === edge.flippedOf ||
      (edge.constructor instanceof StraightEdge &&
        ((this.a.like(edge.a) && this.b.like(edge.b)) ||
          (this.a.like(edge.b) && this.b.like(edge.a))))
    )
  }

  getEdgeT(p: V3): number | undefined {
    assertVectors(p)
    let edgeT = p.minus(this.curve.anchor).dot(this.curve.dir1)
    if (!eq0(this.curve.at(edgeT).distanceTo(p))) {
      return
    }
    edgeT = snap2(edgeT, this.aT, this.bT)
    return this.minT <= edgeT && edgeT <= this.maxT ? edgeT : undefined
  }

  split(t: number): [StraightEdge, StraightEdge] {
    const p = this.curve.at(t)
    return [
      new StraightEdge(
        this.curve,
        this.a,
        p,
        this.aT,
        t,
        undefined,
        this.name + "left",
      ),
      new StraightEdge(
        this.curve,
        p,
        this.b,
        t,
        this.bT,
        undefined,
        this.name + "right",
      ),
    ]
  }
}
