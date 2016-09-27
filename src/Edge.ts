abstract class Edge extends Transformable {
	curve: Curve
	aT: number
	bT: number
	a: V3
	b: V3
	flippedOf: Edge
	name: string
	aDir:V3
	bDir:V3
	reversed:boolean
	id:number


	constructor(curve, a, b, aT, bT, flippedOf, name) {
		super()
		this.curve = curve
		this.a = a
		this.b = b
		this.aT = aT
		this.bT = bT
		this.flippedOf = flippedOf
		this.name = name
		this.reversed = this.aT > this.bT
		this.id = globalId++
	}

	toString(f?):string {
		return `new ${this.constructor.name}(${this.curve}, ${this.a}, ${this.b}, ${this.aT}, ${this.bT}, null, ${this.aDir}, ${this.bDir})`
	}

	/**
	 *
	 * @param {Surface} surface
	 * @returns {number[]}
	 */
	edgeISTsWithSurface(surface) {
		assert(false, this.constructor.name)
	}

	edgeISTsWithPlane(plane) {
		assert(false, this.constructor.name)
	}

	colinearToLine(line) {
		return this.curve instanceof L3 && this.curve.isColinearTo(line)
	}

	tValueInside(t) {
		return this.aT < this.bT
			? NLA.lt(this.aT, t) && NLA.lt(t, this.bT)
			: NLA.lt(this.bT, t) && NLA.lt(t, this.aT)
	}

	clampedT(t) {
		return this.aT < this.bT
			? NLA.clamp(t, this.aT, this.bT)
			: NLA.clamp(t, this.bT, this.aT)
	}


	static forCurveAndTs(curve:Curve, aT:number, bT:number):Edge {
		return Edge.create(curve, curve.at(aT), curve.at(bT), aT, bT, undefined,
			aT < bT ? curve.tangentAt(aT) : curve.tangentAt(aT).negated(),
			aT < bT ? curve.tangentAt(bT) : curve.tangentAt(bT).negated())
	}

	static create(curve:Curve, a:V3, b:V3, aT:number, bT:number, flippedOf:Edge, aDir:V3, bDir:V3, name?:string):Edge {
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

	abstract flipped():Edge

	abstract equals(edge2: any)
}

class PCurveEdge extends Edge {
	aDir:V3
	bDir:V3
	reversed:boolean

	constructor(curve, a, b, aT, bT, flippedOf, aDir, bDir, name) {
		assertNumbers(aT, bT)
		assertVectors(a, b, aDir, bDir)
		assertf(() => curve instanceof L3 || curve instanceof Curve, curve)
		assertf(() => !curve.isValidT || curve.isValidT(aT) && curve.isValidT(bT), aT + ' ' + bT)
		assertf(() => curve.at(aT).like(a), curve.at(aT)+a)
		assertf(() => curve.at(bT).like(b), curve.at(bT)+b)
		assertf(() => curve.tangentAt(aT).likeOrReversed(aDir), curve.tangentAt(aT).sce +' '+ aDir.sce)
		assertf(() => curve.tangentAt(bT).likeOrReversed(bDir))
		super(curve, a, b, aT, bT, flippedOf, name)
		this.aDir = aDir
		this.bDir = bDir
		assert(this.reversed == this.aDir.dot(curve.tangentAt(aT)) < 0, aT+' '+bT+' '+curve.constructor.name+' '+this.aDir.sce+' '+this.bDir.sce + ' '+curve.tangentAt(aT))
	}

	getVerticesNo0() {
		return this.curve.calcSegmentPoints(this.aT, this.bT, this.a, this.b, this.reversed, false)
	}

	points() {
		return this.curve.calcSegmentPoints(this.aT, this.bT, this.a, this.b, this.reversed, true)
	}

	rotViaPlane(normal, reversed) {
		var rot = this.aDir.angleRelativeNormal(this.bDir, normal)
		var counterClockWise = (normal.dot(this.curve.normal) > 0) == !this.reversed
		if (counterClockWise) {
			// counterclockwise rotation, i.e. rot > 0
			if (rot < 0) rot += 2 * Math.PI
		} else {
			if (rot > 0) rot -= 2 * Math.PI
		}
		return rot
	}

	edgeISTsWithSurface(surface) {
		return this.curve.isTsWithSurface(surface).filter(edgeT => {
			var aT = this.aT, bT = this.bT
			edgeT = NLA.snapTo(edgeT, aT)
			edgeT = NLA.snapTo(edgeT, bT)
			if (!this.reversed) {
				if (aT < bT) {
					return aT <= edgeT && edgeT <= bT
				} else {
					return !(bT < edgeT && edgeT < aT)
				}
			} else {
				if (aT > bT) {
					return aT >= edgeT && edgeT >= bT
				} else {
					return !(bT > edgeT && edgeT > aT)
				}
			}
		})
	}

	edgeISTsWithPlane(surface) {
		return this.curve.isTsWithPlane(surface).filter(edgeT => {
			var aT = this.aT, bT = this.bT
			edgeT = NLA.snapTo(edgeT, aT)
			edgeT = NLA.snapTo(edgeT, bT)
			if (!this.reversed) {
				if (aT < bT) {
					return aT <= edgeT && edgeT <= bT
				} else {
					return !(bT < edgeT && edgeT < aT)
				}
			} else {
				if (aT > bT) {
					return aT >= edgeT && edgeT >= bT
				} else {
					return !(bT > edgeT && edgeT > aT)
				}
			}
		})
	}

	tangentAt(t) {
		return !this.reversed ? this.curve.tangentAt(t) : this.curve.tangentAt(t).negated()
	}

	flipped() {
		return this.flippedOf || (this.flippedOf = new PCurveEdge(this.curve, this.b, this.a, this.bT, this.aT, this,
				this.bDir.negated(), this.aDir.negated(), this.name))
	}

	transform(m4, desc) {
		return new PCurveEdge(this.curve.transform(m4), m4.transformPoint(this.a), m4.transformPoint(this.b),
			this.aT, this.bT,
			null,
			m4.transformVector(this.aDir), m4.transformVector(this.bDir), this.name + desc)
	}



	isCoEdge(edge) {
		return this == edge || this == edge.flippedOf ||
			this.curve.isColinearTo(edge.curve) && (
				this.a.like(edge.a) && this.b.like(edge.b)
				|| this.a.like(edge.b) && this.b.like(edge.a)
			)
	}

	like(edge) {
		// TODO this breaks on colinear edges
		return this == edge ||
			edge instanceof Edge &&
			this.curve.isColinearTo(edge.curve) &&
			this.a.like(edge.a) && this.b.like(edge.b)
	}

	static forCurveAndTs(curve:Curve, aT:number, bT:number, name?:string) {
		return new PCurveEdge(curve, curve.at(aT), curve.at(bT), aT, bT, undefined,
			aT < bT ? curve.tangentAt(aT) : curve.tangentAt(aT).negated(),
			aT < bT ? curve.tangentAt(bT) : curve.tangentAt(bT).negated(), name)
	}
}


var StraightEdge = class StraightEdge extends Edge {
	tangent:V3
	curve:L3

	constructor(line, a, b, aT, bT, flippedOf, name) {
		assertInst(L3, line)
		assertNumbers(aT, bT)
		assertVectors(a, b)
		!flippedOf || assertInst(StraightEdge, flippedOf)
		!name || assertf(() => 'string' === typeof name, name)
		assert(line.containsPoint(a), 'line.containsPoint(a)'+line+a)
		assert(line.containsPoint(b), 'line.containsPoint(b)'+line+b)
		super(line, a, b, aT, bT, flippedOf, name)
		this.tangent = this.aT < this.bT ? this.curve.dir1 : this.curve.dir1.negated()
	}

	toSource() {
		return `StraightEdge.throughPoints(${this.a}, ${this.b})`
	}

	getVerticesNo0() {
		return [this.b]
	}

	points() {
		return [this.a, this.b]
	}

	edgeISTsWithPlane(plane) {
		var minT = Math.min(this.aT, this.bT), maxT = Math.max(this.aT, this.bT)
		var edgeT = this.curve.intersectWithPlaneLambda(plane)
		edgeT = NLA.snapTo(edgeT, this.aT)
		edgeT = NLA.snapTo(edgeT, this.bT)
		return (minT <= edgeT && edgeT <= maxT) ? [edgeT] : []
	}

	edgeISTsWithSurface(surface) {
		if (surface instanceof PlaneSurface) {
			return this.edgeISTsWithPlane(surface.plane)
		} else {
			var minT = Math.min(this.aT, this.bT), maxT = Math.max(this.aT, this.bT)
			return surface.isTsForLine(/** @type {L3} */ this.curve)
				.mapFilter(edgeT => {
					edgeT = NLA.snapTo(edgeT, this.aT)
					edgeT = NLA.snapTo(edgeT, this.bT)
					if (minT <= edgeT && edgeT <= maxT) {
						return edgeT
					}
				})
		}
	}

	tangentAt(p) {
		return this.tangent
	}

	flipped() {
		return this.flippedOf || (this.flippedOf = new StraightEdge(this.curve, this.b, this.a, this.bT, this.aT, this, this.name))
	}

	get aDir() { return this.tangent }

	get bDir() { return this.tangent }

	transform(m4, desc) {
		return new StraightEdge(this.curve.transform(m4), m4.transformPoint(this.a), m4.transformPoint(this.b), this.aT, this.bT, null, this.name + desc)
	}

	isCoEdge(edge) {
		return this == edge || this == edge.flippedOf || edge.constructor == StraightEdge && (
				this.a.like(edge.a) && this.b.like(edge.b)
				|| this.a.like(edge.b) && this.b.like(edge.a)
			)
	}

	like(edge) {
		return edge.constructor == StraightEdge && this.a.like(edge.a) && this.b.like(edge.b)
	}

	equals(edge) {
		return edge.constructor == StraightEdge && this.a.equals(edge.a) && this.b.equals(edge.b)
	}

	getEdgeT(p) {
		assertVectors(p)
		var edgeT = p.minus(this.curve.anchor).dot(this.curve.dir1)
		if (!NLA.isZero(this.curve.at(edgeT).distanceTo(p))) { return }
		var minT = Math.min(this.aT, this.bT), maxT = Math.max(this.aT, this.bT)
		edgeT = NLA.snapTo(edgeT, this.aT)
		edgeT = NLA.snapTo(edgeT, this.bT)
		return (minT <= edgeT && edgeT <= maxT) ? edgeT : undefined
	}


	static throughPoints(a:V3, b:V3, name?:string) {
		return new StraightEdge(L3.throughPoints(a, b), a, b, 0, b.minus(a).length(), null, name)
	}
}
