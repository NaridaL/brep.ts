import { SVGPathData } from 'svg-pathdata'
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
	DEG,
	eq,
	eq0,
	fuzzyBetween,
	getIntervals,
	int,
	le,
	lt,
	M4,
	MINUS,
	mod,
	newtonIterate,
	NLA_PRECISION,
	snap2,
	TAU,
	Transformable,
	Tuple3,
	V,
	V3,
} from 'ts3dutils'

import { BezierCurve, Curve, L3, P3, ParabolaCurve, PICurve, PlaneSurface, SemiEllipseCurve, Surface } from './index'

import { abs, ceil, floor, PI, sign } from './math'

export interface Edge {
	readonly aDir: V3
	readonly bDir: V3
}
export abstract class Edge extends Transformable {
	readonly reversed: boolean

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
		assertf(() => !curve.isValidT || (curve.isValidT(aT) && curve.isValidT(bT)), aT, bT, curve)
		//if (curve instanceof PICurve) {
		//    assertf(() => curve.at(aT).to(a).length() < 0.1, ''+curve.at(aT)+a)
		//    assertf(() => curve.at(bT).to(b).length() < 0.1, '' + curve.at(bT) + b)
		//} else {
		assertf(() => curve.at(aT).like(a), () => '' + curve.at(aT) + a + ' aT should have been ' + curve.pointT(a))
		assertf(() => curve.at(bT).like(b), () => '' + curve.at(bT) + b + ' bT should have been ' + curve.pointT(b))
		//}
		assertf(() => fuzzyBetween(aT, curve.tMin, curve.tMax), aT, curve.tMin, curve.tMax)
		assertf(() => fuzzyBetween(bT, curve.tMin, curve.tMax), bT, curve.tMin, curve.tMax)
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

	static forCurveAndTs(curve: Curve, aT: number = curve.tMin, bT: number = curve.tMax): Edge {
		return Edge.create(
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

	static create(
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
			return new StraightEdge(curve, a, b, aT, bT, flippedOf as StraightEdge, name)
		} else {
			return new PCurveEdge(curve, a, b, aT, bT, flippedOf as PCurveEdge, aDir, bDir, name)
		}
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
		if (sts.some(info => isNaN(info.tThis) || isNaN(info.tOther))) {
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
			assert(edge.b.like(edges[j].a), `edges[${i}].b != edges[${j}].a (${edges[i].b.sce} != ${edges[j].a.sce})`)
		})
	}

	static ngon(n: int = 3, radius: number = 1): Edge[] {
		return StraightEdge.chain(arrayFromFunction(n, i => V3.polar(radius, TAU * i / n)))
	}

	static star(pointCount: int = 5, r0: number = 1, r1: number = 0.5): Edge[] {
		const vertices = arrayFromFunction(pointCount * 2, i =>
			V3.polar(0 == i % 2 ? r0 : r1, TAU * i / pointCount / 2),
		)
		return StraightEdge.chain(vertices)
	}

	static reversePath(path: Edge[], doReverse: boolean = true): Edge[] {
		return doReverse ? arrayFromFunction(path.length, i => path[path.length - 1 - i].flipped()) : path
	}

	/**
	 * Create an axis-aligned rectangle of edges on the XY-plane with the bottom-left corner on the origin.
	 * @param width
	 * @param height
	 */
	static rect(width: number = 1, height: number = width): Edge[] {
		const vertices = [new V3(0, 0, 0), new V3(width, 0, 0), new V3(width, height, 0), new V3(0, height, 0)]
		return StraightEdge.chain(vertices)
	}

	static reuleaux(n: int = 3, radius: number = 1): Edge[] {
		assert(3 <= n)
		assert(1 == n % 2)
		const corners = arrayFromFunction(n, i => V3.polar(radius, TAU * i / n))
		return arrayFromFunction(n, i => {
			const aI = (i + floor(n / 2)) % n,
				bI = (i + ceil(n / 2)) % n
			const a = corners[aI],
				b = corners[bI]
			const center = corners[i]
			const f1 = center.to(a),
				curve = new SemiEllipseCurve(center, f1, V3.Z.cross(f1))
			return Edge.create(curve, a, b, 0, curve.pointT(b), undefined, V3.Z.cross(f1), V3.Z.cross(center.to(b)))
		})
	}

	static round(edges: Edge[], radius: number) {
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
				if (!center) throw new Error('tangential curves')
				const cornerA = center.plus(l1inside.toLength(-radius))
				const cornerB = center.plus(l2inside.toLength(-radius))
				const f1 = l1inside.toLength(-radius)
				const curve = new SemiEllipseCurve(center, f1, normal.cross(f1).toLength(radius))
				const cornerEdge = Edge.create(
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
				return Edge.arbitraryCorner(edge, nextEdge, radius)
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
			const newEdge = Edge.create(edge.curve, a, b, aT, bT, undefined, aDir, bDir)
			return !nextCorner ? newEdge : [newEdge, nextCorner]
		})
		return result
	}

	static arbitraryCorner(e1: Edge, e2: Edge, radius: number) {
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

		const startT1 = e1.bT - radius * sign(e1.deltaT()) / e1.bDir.length()
		const startT2 = e2.aT + radius * sign(e2.deltaT()) / e2.aDir.length()
		const [t1, t2] = newtonIterate(f, [startT1, startT2])
		const cornerA = e1.curve.at(t1)
		const cornerB = e2.curve.at(t2)
		const dp1 = c1.tangentAt(t1),
			dp2 = c2.tangentAt(t2)
		const virtualPlaneNormal = dp1.cross(dp2)
		const normal1 = virtualPlaneNormal.cross(dp1).unit()
		const f1 = normal1.toLength(-radius)
		const center = cornerA.minus(f1)
		const curve = new SemiEllipseCurve(center, f1, virtualPlaneNormal.cross(f1).toLength(radius))
		const cornerEdge = Edge.create(
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

	static pathFromSVG(pathString: string): Edge[] {
		let currentPos: V3 = undefined!
		const parsed: any[] = new SVGPathData(pathString)
			.toAbs()
			.normalizeHVZ()
			.sanitize(NLA_PRECISION)
			.annotateArcs().commands
		const path: Edge[] = []
		for (const c of parsed) {
			assert('x' in c && 'y' in c)
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
					const curve = ParabolaCurve.quadratic(currentPos, c1, endPos).rightAngled()
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
						.map(n => n * PI)
						.filter(stop => phiMin <= stop && stop <= phiMax)
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
						const curve = new SemiEllipseCurve(center, gtPI ? f1.negated() : f1, gtPI ? f2.negated() : f2)
						const a = phi1 == t1 ? currentPos : phi2 == t1 ? endPos : curve.at(aT)
						const b = phi1 == t2 ? currentPos : phi2 == t2 ? endPos : curve.at(bT)
						return new PCurveEdge(curve, a, b, aT, bT, undefined, curve.tangentAt(aT), curve.tangentAt(bT))
					})
					path.push(...(c.phiDelta > 0 ? edges : Edge.reversePath(edges)))
					break
				}
			}
			currentPos = endPos
		}
		return path
	}

	abstract tangentAt(t: number): V3

	toString(): string {
		return callsce(
			'new ' + this.constructor.name,
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

	split(t: number): [Edge, Edge] {
		const p = this.curve.at(t)
		const pDir = this.tangentAt(t)
		return [
			Edge.create(this.curve, this.a, p, this.aT, t, undefined, this.aDir, pDir, this.name + 'left'),
			Edge.create(this.curve, p, this.b, t, this.bT, undefined, pDir, this.bDir, this.name + 'left'),
		]
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
		return this.aT < this.bT ? lt(this.aT, t) && lt(t, this.bT) : lt(this.bT, t) && lt(t, this.aT)
	}

	isValidT(t: number): boolean {
		return this.aT < this.bT ? le(this.aT, t) && le(t, this.bT) : le(this.bT, t) && le(t, this.aT)
	}

	clampedT(t: number): number {
		return this.aT < this.bT ? clamp(t, this.aT, this.bT) : clamp(t, this.bT, this.aT)
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
			(edge instanceof Edge && this.curve.isColinearTo(edge.curve) && this.a.like(edge.a) && this.b.like(edge.b))
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
		const flipped =
			false !== edgeAT ? this.tangentAt(edgeAT).dot(edge.aDir) : this.tangentAt(edge.bT).dot(edge.bDir)
		return !(le(edge.maxT, this.minT) || le(this.maxT, edge.minT))
	}

	getAABB(): AABB {
		const min: Tuple3<number> = [Infinity, Infinity, Infinity],
			max: Tuple3<number> = [-Infinity, -Infinity, -Infinity]
		this.curve.roots().forEach((ts, dim) => {
			ts.forEach(t => {
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
			arrayRange(0, loop1.length, 1).some(offset =>
				loop1.every((edge, i) => edge.equals(loop2[(offset + i) % loop1.length])),
			)
		)
	}
}

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
			assertf(() => curve.tangentAt(aT).likeOrReversed(aDir), '' + aT + curve.tangentAt(aT).sce + ' ' + aDir.sce)
			assertf(() => curve.tangentAt(bT).likeOrReversed(bDir))
		}
		assert(
			this.reversed === this.aDir.dot(curve.tangentAt(aT)) < 0,
			aT +
				' ' +
				bT +
				' ' +
				curve.constructor.name +
				' ' +
				this.aDir.sce +
				' ' +
				this.bDir.sce +
				' ' +
				curve.tangentAt(aT),
		)
		assert(
			this.reversed === this.bDir.dot(curve.tangentAt(bT)) < 0,
			aT +
				' ' +
				bT +
				' ' +
				curve.constructor.name +
				' ' +
				this.aDir.sce +
				' ' +
				this.bDir.sce +
				' ' +
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
		return callsce(
			'new PCurveEdge',
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
		return this.curve.calcSegmentPoints(this.aT, this.bT, this.a, this.b, this.reversed, false)
	}

	pointsCount(): int {
		return this.points().length
	}

	points(): V3[] {
		return this.curve.calcSegmentPoints(this.aT, this.bT, this.a, this.b, this.reversed, true)
	}

	edgeISTsWithSurface(surface: Surface): number[] {
		return this.curve
			.isTsWithSurface(surface)
			.map(edgeT => snap2(edgeT, this.aT, this.bT))
			.filter(edgeT => this.minT <= edgeT && edgeT <= this.maxT)
	}

	edgeISTsWithPlane(surface: P3): number[] {
		return this.curve
			.isTsWithPlane(surface)
			.map(edgeT => snap2(edgeT, this.aT, this.bT))
			.filter(edgeT => this.minT <= edgeT && edgeT <= this.maxT)
	}

	tangentAt(t: number): V3 {
		return !this.reversed ? this.curve.tangentAt(t) : this.curve.tangentAt(t).negated()
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
			'' + this.name + desc,
		) as this
	}

	isCoEdge(edge: Edge): boolean {
		return (
			this === edge ||
			this === edge.flippedOf ||
			(this.curve.isColinearTo(edge.curve) &&
				((this.a.like(edge.a) && this.b.like(edge.b)) || (this.a.like(edge.b) && this.b.like(edge.a))))
		)
	}
}

export class StraightEdge extends Edge {
	readonly tangent: V3
	readonly curve!: L3

	constructor(line: L3, a: V3, b: V3, aT: number, bT: number, public flippedOf?: StraightEdge, name?: string) {
		super(line, a, b, aT, bT, flippedOf, name)
		assertInst(L3, line)
		!flippedOf || assertInst(StraightEdge, flippedOf)
		!name || assertf(() => 'string' === typeof name, name)
		assert(!a.like(b), '!a.like(b)' + a + b) // don't put in super as it will break full ellipse
		this.tangent = this.aT < this.bT ? this.curve.dir1 : this.curve.dir1.negated()
	}

	get aDir() {
		return this.tangent
	}

	get bDir() {
		return this.tangent
	}

	static throughPoints(a: V3, b: V3, name?: string) {
		return new StraightEdge(L3.throughPoints(a, b, 0, a.to(b).length()), a, b, 0, a.to(b).length(), undefined, name)
	}

	/**
	 * Create a list of StraightEdges from a list of vertices.
	 * @param vertices
	 * @param closed Whether to connect the first and last vertices. Defaults to true.
	 * @returns
	 */
	static chain(vertices: V3[], closed: boolean = true): StraightEdge[] {
		const vc = vertices.length
		return arrayFromFunction(closed ? vc : vc - 1, i =>
			StraightEdge.throughPoints(vertices[i], vertices[(i + 1) % vc]),
		)
	}

	toSource(): string {
		return callsce('new StraightEdge', this.curve, this.a, this.b, this.aT, this.bT)
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
				.map(edgeT => snap2(edgeT, this.aT, this.bT))
				.filter(edgeT => this.minT <= edgeT && edgeT <= this.maxT)
		}
	}

	tangentAt() {
		return this.tangent
	}

	flipped(): StraightEdge {
		return (
			this.flippedOf ||
			(this.flippedOf = new StraightEdge(this.curve, this.b, this.a, this.bT, this.aT, this, this.name))
		)
	}

	transform(m4: M4, desc?: string): this {
		const lineDir1TransLength = m4.transformVector(this.curve.dir1).length()
		return new StraightEdge(
			this.curve.transform(m4),
			m4.transformPoint(this.a),
			m4.transformPoint(this.b),
			this.aT * lineDir1TransLength,
			this.bT * lineDir1TransLength,
			undefined,
			'' + this.name + desc,
		) as this
	}

	isCoEdge(edge: Edge): boolean {
		return (
			this === edge ||
			this === edge.flippedOf ||
			(edge.constructor === StraightEdge &&
				((this.a.like(edge.a) && this.b.like(edge.b)) || (this.a.like(edge.b) && this.b.like(edge.a))))
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
}
