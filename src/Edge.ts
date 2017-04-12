abstract class Edge extends Transformable {
	curve: Curve
	aT: number
	bT: number
	a: V3
	b: V3
	flippedOf: Edge
	name: string
	aDir: V3
	bDir: V3
	reversed: boolean
	id: number
	aDDT: V3
	bDDT: V3

	get minT() { return Math.min(this.aT, this.bT) }

	get maxT() { return Math.max(this.aT, this.bT) }

	abstract tangentAt(t: number): V3

	constructor(curve, a, b, aT, bT, flippedOf, name) {
		super()
		this.curve = curve
		this.a = a
		this.b = b
		this.aT = aT
		this.bT = bT
        assert(!eq(aT, bT))
		this.flippedOf = flippedOf
		this.name = name
		this.reversed = this.aT > this.bT
		this.id = globalId++
	}

	toString(f?): string {
		return `new ${this.constructor.name}(${this.curve}, ${this.a}, ${this.b}, ${this.aT}, ${this.bT}, null, ${this.aDir}, ${this.bDir})`
	}

	abstract edgeISTsWithSurface(surface: Surface): number[]

	/**
	 * Returns the intersections of the edge with the plane.
	 * Values are snapped to aT and bT, ie aT === t || !NLA.eq(aT, t)
	 * @param plane
	 */
	abstract edgeISTsWithPlane(plane: P3): number[]

	colinearToLine(line): boolean {
		return this.curve instanceof L3 && this.curve.isColinearTo(line)
	}

	tValueInside(t) {
		return this.aT < this.bT
			? NLA.lt(this.aT, t) && NLA.lt(t, this.bT)
			: NLA.lt(this.bT, t) && NLA.lt(t, this.aT)
	}

	isValidT(t: number): boolean {
		return this.aT < this.bT
			? NLA.le(this.aT, t) && NLA.le(t, this.bT)
			: NLA.le(this.bT, t) && NLA.le(t, this.aT)
	}

	clampedT(t) {
		return this.aT < this.bT
			? NLA.clamp(t, this.aT, this.bT)
			: NLA.clamp(t, this.bT, this.aT)
	}


	static forCurveAndTs(curve: Curve, aT: number = curve.tMin, bT: number = curve.tMax): Edge {
		return Edge.create(curve, curve.at(aT), curve.at(bT), aT, bT, undefined,
			aT < bT ? curve.tangentAt(aT) : curve.tangentAt(aT).negated(),
			aT < bT ? curve.tangentAt(bT) : curve.tangentAt(bT).negated())
	}

	static create(curve: Curve, a: V3, b: V3, aT: number, bT: number, flippedOf: Edge, aDir: V3, bDir: V3, name?: string): Edge {
		if (curve instanceof L3) {
			return new StraightEdge(curve, a, b, aT, bT, flippedOf, name)
		} else {
			return new PCurveEdge(curve, a, b, aT, bT, flippedOf, aDir, bDir, name)
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
			[e2, e1] = [e1, e2]
		}
		let sts = e1.curve.isInfosWithCurve(e2.curve)
		if (sts.some(info => isNaN(info.tThis) || isNaN(info.tOther))) {
			console.log(e1.sce)
			console.log(e2.sce)
			assert(false)
		}
		return sts.some(
			/// (  e1.aT < tThis < e1.bT  )  &&  (  e2.aT < tOther < e2.bT  )
			({tThis, tOther}) => {
				return e1.tValueInside(tThis) && e2.tValueInside(tOther)
			})
	}

	abstract flipped(): Edge

    /**
     * this is equals-equals. "isColinearTo" might make more sense but can't be used, because you can't get a
     * consistent hashCode for colinear curves
     * @param obj
     * @returns {boolean}
     */
    equals(obj): boolean {
        return this === obj ||
            this.constructor == obj.constructor
            && this.a.equals(obj.a)
            && this.b.equals(obj.b)
            && this.curve.equals(obj.curve)
    }

    hashCode(): int {
        let hashCode = 0
        hashCode = hashCode * 31 + this.a.hashCode()
        hashCode = hashCode * 31 + this.b.hashCode()
        hashCode = hashCode * 31 + this.curve.hashCode()
        return hashCode | 0
    }

    like(edge) {
        // TODO this breaks on colinear edges,
        // TODO: what, where?
        return this === edge ||
            edge instanceof Edge &&
            this.curve.isColinearTo(edge.curve)
            && this.a.like(edge.a)
            && this.b.like(edge.b)
    }

	abstract getVerticesNo0(): V3[]

	abstract pointsCount(): int

	static assertLoop(edges: Edge[]): void {
		edges.forEach((edge, i) => {
			const j = (i + 1) % edges.length
			assert(edge.b.like(edges[j].a), `edges[${i}].b != edges[${j}].a (${edges[i].b.sce} != ${edges[j].a.sce})`)
		})
	}

	getCanon() {
		return this.reversed
			? this.flipped()
			: this
	}

	overlaps(edge: Edge, noback?: boolean): boolean {
		assert(this.curve.isColinearTo(edge.curve))
		const edgeAT = this.curve.containsPoint(edge.a) && this.curve.pointT(edge.a)
        const edgeBT = this.curve.containsPoint(edge.b) && this.curve.pointT(edge.b)
        if (false === edgeAT && false === edgeBT) {
		    return noback ? false : edge.overlaps(this, true)
        }
        const flipped = false !== edgeAT ? this.tangentAt(edgeAT).dot(edge.aDir) : this.tangentAt(edgeBT).dot(edge.bDir)
		return !(NLA.le(edgeMaxT, this.minT) || NLA.le(this.maxT, edgeMinT))
	}

	getAABB(): AABB {
        const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity]
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

	abstract isCoEdge(other: Edge): boolean

	static ngon(count: int, radius: number = 1) {
		return StraightEdge.chain(NLA.arrayFromFunction(count, i => V3.polar(radius, TAU * i / count)))
	}

	static star(pointCount: int = 5, r0: number = 1, r1: number = 0.5) {
		return StraightEdge.chain(NLA.arrayFromFunction(pointCount * 2, i => V3.polar(0 == i % 2 ? r0 : r1, TAU * i / pointCount / 2)))
	}

	static reverseLoop(loop: Edge[]) {
		return NLA.arrayFromFunction(loop.length, i => loop[loop.length - 1 - i].flipped())
	}
    static pathFromSVG(pathString: String): Edge[] {
        let currentPos
        const parsed: {type: int, x: number, y: number, relative?: boolean}[] =
	        new SVGPathData(pathString).toAbs().normalizeHVZ().sanitize(NLA_PRECISION).commands
        const path: Edge[] = []
        for (const c of parsed) {
            const endPos = ('x' in c && 'y' in c) && new V3(c.x, c.y, 0)
            switch (c.type) {
	            case SVGPathData.LINE_TO:
                    path.push(StraightEdge.throughPoints(currentPos, endPos))
                    break
                case SVGPathData.CURVE_TO: {
                    const c1 = new V3(c.x1, c.y1, 0)
                    const c2 = new V3(c.x2, c.y2, 0)
                    const curve = new BezierCurve(currentPos, c1, c2, endPos, 0, 1)
                    const edge = new PCurveEdge(curve, currentPos, endPos, 0, 1, undefined, curve.tangentAt(0), curve.tangentAt(1))
                    path.push(edge)
                    break
                }
	            case SVGPathData.QUAD_TO: {
		            const c1 = new V3(c.x1, c.y1, 0)
		            const curve = ParabolaCurve.quadratic(currentPos, c1, endPos).rightAngled()
		            const edge = new PCurveEdge(curve, currentPos, endPos, curve.tMin, curve.tMax, undefined, curve.tangentAt(curve.tMin), curve.tangentAt(curve.tMax))
		            path.push(edge)
		            break
	            }
                case SVGPathData.ARC_TO: {
	                let {rx, ry, xAxisRotation, largeArc, sweep} = c
	                rx = abs(rx)
	                ry = abs(ry)
	                if ('a' == c.code) {
		                endPos = endPos.plus(currentPos)
	                }
	                const rads = xAxisRotation * DEG
	                const midPoint = currentPos.minus(endPos).times(0.5)
	                const midPointTransformed = M4.rotationZ(-rads).transformPoint(midPoint)
	                const testValue = midPointTransformed.x ** 2 / rx ** 2 + midPointTransformed.y ** 2 / ry ** 2
	                console.log(testValue)
	                if (testValue > 1) {
		                rx *= Math.sqrt(testValue)
		                ry *= Math.sqrt(testValue)
	                }
	                const temp = (rx ** 2 * midPointTransformed.y ** 2 + ry ** 2 * midPointTransformed.x ** 2)
	                const centerTransformedScale = (largeArc != sweep ? 1 : -1) * Math.sqrt(max(0, (rx ** 2 * ry ** 2 - temp) / temp))
	                const centerTransformed = V(rx * midPointTransformed.y / ry, -ry * midPointTransformed.x / rx).times(centerTransformedScale)
	                const center = M4.rotationZ(rads).transformPoint(centerTransformed).plus(currentPos.plus(endPos).times(0.5))
	                let f1 = new V3(rx * cos(rads), rx * sin(rads), 0)
	                const f2 = new V3(ry * -sin(rads), ry * cos(rads), 0)
	                let curve = new SemiEllipseCurve(center, f1, f2)
	                let aT = curve.pointT(currentPos, PI)
	                let bT = curve.pointT(endPos, PI)
	                if (aT < bT != sweep) {
		                assert((aT != PI) || (bT != PI))
		                if (aT == PI) {
			                aT = -PI
		                } else if (bT == PI) {
			                bT = -PI
		                } else {
			                f1 = f1.negated()
			                curve = new SemiEllipseCurve(center, f1, f2)
			                aT = curve.pointT(currentPos, PI)
			                bT = curve.pointT(endPos, PI)
		                }
	                }
	                const edge = new PCurveEdge(curve, currentPos, endPos, aT, bT, undefined, curve.tangentAt(aT).times(sign(bT - aT)), curve.tangentAt(bT).times(sign(bT - aT)))
	                path.push(edge)
	                break
                }
            }
            currentPos = endPos
        }
        return path
    }

    deltaT() {
        return this.bT - this.aT
    }
}

class PCurveEdge extends Edge {

	constructor(curve, a, b, aT, bT, flippedOf, aDir, bDir, name?) {
		assertNumbers(aT, bT)
		assertVectors(a, b, aDir, bDir)
		assertf(() => curve instanceof L3 || curve instanceof Curve, curve)
		assertf(() => !curve.isValidT || curve.isValidT(aT) && curve.isValidT(bT), aT + ' ' + bT)
		assertf(() => curve.at(aT).like(a), curve.at(aT) + a)
		assertf(() => curve.at(bT).like(b), curve.at(bT) + b)
		assertf(() => curve.tangentAt(aT).likeOrReversed(aDir), '' + aT + curve.tangentAt(aT).sce + ' ' + aDir.sce)
		assertf(() => curve.tangentAt(bT).likeOrReversed(bDir))
        assertf(() => fuzzyBetween(aT, curve.tMin, curve.tMax))
        assertf(() => fuzzyBetween(bT, curve.tMin, curve.tMax))
		super(curve, a, b, clamp(aT, curve.tMin, curve.tMax), clamp(bT, curve.tMin, curve.tMax), flippedOf, name)
		this.aDir = aDir
		this.bDir = bDir
		assert(this.reversed === this.aDir.dot(curve.tangentAt(aT)) < 0, aT + ' ' + bT + ' ' + curve.constructor.name + ' ' + this.aDir.sce + ' ' + this.bDir.sce + ' ' + curve.tangentAt(aT))
		assert(this.reversed === this.bDir.dot(curve.tangentAt(bT)) < 0, aT + ' ' + bT + ' ' + curve.constructor.name + ' ' + this.aDir.sce + ' ' + this.bDir.sce + ' ' + curve.tangentAt(aT))
	}

	getVerticesNo0() {
		return this.curve.calcSegmentPoints(this.aT, this.bT, this.a, this.b, this.reversed, false)
	}

	pointsCount() {
		return this.points().length
	}

	points() {
		return this.curve.calcSegmentPoints(this.aT, this.bT, this.a, this.b, this.reversed, true)
	}

	rotViaPlane(normal, reversed) {
		let rot = this.aDir.angleRelativeNormal(this.bDir, normal)
		const counterClockWise = (normal.dot(this.curve.normal) > 0) === !this.reversed
		if (counterClockWise) {
			// counterclockwise rotation, i.e. rot > 0
			if (rot < 0) rot += 2 * Math.PI
		} else {
			if (rot > 0) rot -= 2 * Math.PI
		}
		return rot
	}

	edgeISTsWithSurface(surface) {
        return this.curve.isTsWithSurface(surface)
			.map(edgeT => NLA.snap(NLA.snap(edgeT, this.aT), this.bT))
			.filter(edgeT => this.minT <= edgeT && edgeT <= this.maxT)
	}

	edgeISTsWithPlane(surface) {
        return this.curve.isTsWithPlane(surface)
			.map(edgeT => NLA.snap(NLA.snap(edgeT, this.aT), this.bT))
			.filter(edgeT => this.minT <= edgeT && edgeT <= this.maxT)
	}

	tangentAt(t) {
		return !this.reversed ? this.curve.tangentAt(t) : this.curve.tangentAt(t).negated()
	}

	flipped() {
		return this.flippedOf || (this.flippedOf = new PCurveEdge(this.curve, this.b, this.a, this.bT, this.aT, this,
				this.bDir.negated(), this.aDir.negated(), this.name))
	}

	transform(m4, desc): this {
		return new PCurveEdge(this.curve.transform(m4), m4.transformPoint(this.a), m4.transformPoint(this.b),
			this.aT, this.bT,
			null,
			m4.transformVector(this.aDir), m4.transformVector(this.bDir), this.name + desc) as this
	}


	isCoEdge(edge) {
		return this === edge || this === edge.flippedOf ||
			this.curve.isColinearTo(edge.curve) && (
				this.a.like(edge.a) && this.b.like(edge.b)
				|| this.a.like(edge.b) && this.b.like(edge.a)
			)
	}

    static forCurveAndTs(curve: Curve, aT: number, bT: number, name?: string) {
        return new PCurveEdge(curve, curve.at(aT), curve.at(bT), aT, bT, undefined,
            aT < bT ? curve.tangentAt(aT) : curve.tangentAt(aT).negated(),
            aT < bT ? curve.tangentAt(bT) : curve.tangentAt(bT).negated(), name)
    }

	get aDDT() {
		let ddt = this.curve.ddt(this.aT)
		return this.reversed ? ddt.negated() : ddt
	}

	get bDDT() {
		let ddt = this.curve.ddt(this.bT)
		return this.reversed ? ddt.negated() : ddt
	}
}


class StraightEdge extends Edge {
	tangent: V3
	curve: L3
	// flippedOf: StraightEdge

	constructor(line: L3, a: V3, b: V3, aT: number, bT: number, flippedOf?: StraightEdge, name?: string) {
		assertInst(L3, line)
		assertNumbers(aT, bT)
		assertVectors(a, b)
		!flippedOf || assertInst(StraightEdge, flippedOf)
		!name || assertf(() => 'string' === typeof name, name)
		assert(line.at(aT).like(a), 'line.at(aT).like(a)' + aT + line + a)
		assert(line.at(bT).like(b), 'line.at(bT).like(b)' + bT + line + b)
        assert(!a.like(b), '!a.like(b)' + a + b)
		super(line, a, b, aT, bT, flippedOf, name)
		this.tangent = this.aT < this.bT ? this.curve.dir1 : this.curve.dir1.negated()
	}

    toSource() {
        //return `StraightEdge.throughPoints(${this.a}, ${this.b})`
        return `new StraightEdge(${this.curve}, ${this.a}, ${this.b}, ${this.aT}, ${this.bT})`
    }

    toString() {
        //return `StraightEdge.throughPoints(${this.a}, ${this.b})`
        return `new StraightEdge(${this.curve}, ${this.a}, ${this.b}, ${this.aT}, ${this.bT})`
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

	edgeISTsWithPlane(plane): number[] {
		let edgeT = this.curve.intersectWithPlaneLambda(plane)
		edgeT = NLA.snap(edgeT, this.aT)
		edgeT = NLA.snap(edgeT, this.bT)
		return (this.minT <= edgeT && edgeT <= this.maxT) ? [edgeT] : []
	}

	edgeISTsWithSurface(surface): number[] {
		if (surface instanceof PlaneSurface) {
			return this.edgeISTsWithPlane(surface.plane)
		} else {
			return surface.isTsForLine(this.curve)
				.map(edgeT => NLA.snap(NLA.snap(edgeT, this.aT), this.bT))
				.filter(edgeT => this.minT <= edgeT && edgeT <= this.maxT)
		}
	}

	tangentAt(p) {
		return this.tangent
	}

	flipped(): StraightEdge {
		return this.flippedOf || (this.flippedOf = new StraightEdge(this.curve, this.b, this.a, this.bT, this.aT, this, this.name))
	}

	get aDir() {
		return this.tangent
	}

	get bDir() {
		return this.tangent
	}

	transform(m4, desc): this {
	    const lineDir1TransLength = m4.transformVector(this.curve.dir1).length()
		return new StraightEdge(
			this.curve.transform(m4),
			m4.transformPoint(this.a),
			m4.transformPoint(this.b), this.aT * lineDir1TransLength, this.bT * lineDir1TransLength, null, this.name + desc) as any
	}

	isCoEdge(edge): boolean {
		return this === edge || this === edge.flippedOf || edge.constructor === StraightEdge && (
				this.a.like(edge.a) && this.b.like(edge.b)
				|| this.a.like(edge.b) && this.b.like(edge.a)
			)
	}

	getEdgeT(p: V3): number|undefined {
		assertVectors(p)
		let edgeT = p.minus(this.curve.anchor).dot(this.curve.dir1)
		if (!NLA.eq0(this.curve.at(edgeT).distanceTo(p))) {
			return
		}
		edgeT = NLA.snap(edgeT, this.aT)
		edgeT = NLA.snap(edgeT, this.bT)
		return (this.minT <= edgeT && edgeT <= this.maxT) ? edgeT : undefined
	}


	static throughPoints(a: V3, b: V3, name?: string) {
		return new StraightEdge(L3.throughPoints(a, b), a, b, 0, b.minus(a).length(), null, name)
	}

	/**
	 * Create a list of StraightEdges from a list of vertices.
	 * @param vertices
	 * @param closed Whether to connect the first and last vertices. Defaults to true.
	 * @returns
	 */
	static chain(vertices: V3[], closed: boolean = true): StraightEdge[] {
		const vc = vertices.length
		return NLA.arrayFromFunction(closed ? vc : vc - 1,
			i => StraightEdge.throughPoints(vertices[i], vertices[(i + 1) % vc]))
	}

}
StraightEdge.prototype.aDDT = V3.O
StraightEdge.prototype.bDDT = V3.O

NLA.registerClass(StraightEdge)