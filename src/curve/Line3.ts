/**
 * 3-dimensional line
 */
class L3 extends Curve {

	constructor(readonly anchor: V3, // line anchor
	            readonly dir1: V3, // normalized line dir
	            tMin: number = -4096,
	            tMax: number = 4096) {
		super(tMin, tMax)
		assertVectors(anchor, dir1)
		assert(dir1.hasLength(1), 'dir must be unit' + dir1)
		assertf(() => !Number.isNaN(anchor.x))
	}

    roots(): number[][] {
        return [[], [], []]
    }

	containsPoint(p: V3): boolean {
		assertVectors(p)
		const dist = this.distanceToPoint(p)
		assertNumbers(dist)
		return eq0(dist)
	}

	likeCurve(curve: Curve): boolean {
		return this == curve ||
			Object.getPrototypeOf(curve) == L3.prototype
			&& this.anchor.like(curve.anchor)
			&& this.dir1.like(curve.dir1)
	}

	equals(obj: any): boolean {
		return this == obj ||
			Object.getPrototypeOf(obj) == L3.prototype
			&& this.anchor.equals(obj.anchor)
			&& this.dir1.equals(obj.dir1)
	}

	isColinearTo(obj: Curve): boolean {
		return obj instanceof L3
			&& this.containsPoint(obj.anchor)
			&& eq(1, Math.abs(this.dir1.dot(obj.dir1)))
	}

	distanceToLine(line: L3): number {
		assertInst(L3, line)
		if (this.isParallelToLine(line)) {
			return this.distanceToPoint(line.anchor)
		}
		const dirCross1 = this.dir1.cross(line.dir1).unit()
		const anchorDiff = this.anchor.minus(line.anchor)
		return Math.abs(anchorDiff.dot(dirCross1))
	}

	distanceToPoint(x: V3): number {
		assertVectors(x)
		// See http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
		const t = x.minus(this.anchor).dot(this.dir1)
		return this.at(t).distanceTo(x)

		//return x.minus(this.anchor).cross(x.minus(this.anchor.plus(this.dir1))).length()
	}

	asSegmentDistanceToPoint(x, sStart, sEnd) {
		let t = x.minus(this.anchor).dot(this.dir1)
		t = clamp(t, sStart, sEnd)
		return this.at(t).minus(x).length()
	}

	asSegmentDistanceToLine(line, sStart, sEnd) {
		assertInst(L3, line)
		const dirCross = this.dir1.cross(line.dir1)
		const div = dirCross.squared()
		if (eq0(div)) {
			return null
		} // lines parallel
		const anchorDiff = line.anchor.minus(this.anchor)
		// check if distance is zero (see also L3.distanceToLine)
		if (!eq0(anchorDiff.dot(dirCross.unit()))) {
			return null
		}
		let t = this.infoClosestToLine(line).t
		t = clamp(t, sStart, sEnd)
		return this.at(clamp(t, sStart, sEnd))
	}

	at(t: number): V3 {
		assertNumbers(t)
		return this.anchor.plus(this.dir1.times(t))
	}

	/**
	 * This function returns lambda for a given point x
	 *
	 * Every point x on this line is described by the equation
	 *      x = this.anchor + lambda * this.dir1 | - this.anchor
	 *      x - this.anchor = lambda * this.dir1 | DOT this.dir1
	 *      (x - this.anchor) DOT this.dir1 = lambda (dir1² is 1 as |dir1| == 1)
	 *
	 *  @param x
	 *  @returns {number}
	 */
	pointT(x: V3) {
		assertVectors(x)
		const t = x.minus(this.anchor).dot(this.dir1)
		return t
	}

	/**
	 * Returns true if the line is parallel (this.dir = line.dir || this.dir = -line.dir) to the argument.
	 */
	isParallelToLine(line: L3): boolean {
		assertInst(L3, line)
		// we know that 1 == this.dir1.length() == line.dir1.length(), we can check for parallelity simpler than
		// isParallelTo()
		return eq(1, Math.abs(this.dir1.dot(line.dir1)))
	}

	angleToLine(line: L3): number {
		assertInst(L3, line)
		return this.dir1.angleTo(line.dir1)
	}

	/**
	 *
	 * @param line
	 * @returns {boolean} If the distance between the lines is zero
	 */
	intersectsLine(line: L3): boolean {
		return eq0(this.distanceToLine(line))
	}

	isInfosWithCurve(curve: Curve) {
		assertInst(L3, curve)

		const dirCross = this.dir1.cross(curve.dir1)
		const div = dirCross.squared()
		if (eq0(div)) {
			// lines are parallel
			return []
		}
		const anchorDiff = curve.anchor.minus(this.anchor)
		if (eq0(anchorDiff.dot(dirCross))) {
			const tThis = anchorDiff.cross(curve.dir1).dot(dirCross) / div
			const tOther = anchorDiff.cross(this.dir1).dot(dirCross) / div
			const p = this.at(tThis)
			return [{tThis: tThis, tOther: tOther, p: p}]
		}
		return []
	}

	isInfoWithLine(line: L3) {
		// todo infos?
		assertInst(L3, line)
		const dirCross = this.dir1.cross(line.dir1)
		const div = dirCross.squared()
		if (eq0(div)) {
			return null
		} // lines parallel
		const anchorDiff = line.anchor.minus(this.anchor)
		// check if distance is zero (see also L3.distanceToLine)
		if (!eq0(anchorDiff.dot(dirCross.unit()))) {
			return null
		}
		const t = anchorDiff.cross(line.dir1).dot(dirCross) / div
		return this.at(t)
	}

	/**
	 * returns s and t with this.at(s) == line.at(t)
	 */
	intersectionLineST(line: L3): { s: number, t: number } {
		// the two points on two lines the closest two each other are the ones whose
		// connecting
		// TODO Where does this come from?
		// TODO: return value when no IS?
		assertInst(L3, line)
		const dirCross = this.dir1.cross(line.dir1)
		const div = dirCross.squared()
		const anchorDiff = line.anchor.minus(this.anchor)
		const s = anchorDiff.cross(this.dir1).dot(dirCross) / div
		const t = anchorDiff.cross(line.dir1).dot(dirCross) / div
		return {s: s, t: t}
		//console.log(segmentIntersectsRay, a, b, "ab", ab, "p", p, "dir", dir, s > 0 && t / div >= 0 && t / div <= 1,
		// "s", s, "t", t, "div", div)
	}

	ddt(t: number): V3 {
		return V3.O
	}

	getConstructorParameters(): any[] {
		return [this.anchor, this.dir1]
	}

	closestTToPoint(p: V3): number {
		// similar logic as pointT; we project the vector (anchor -> p) onto dir1, then add anchor back to it
		const nearestT = p.minus(this.anchor).dot(this.dir1)
		return nearestT
	}

	infoClosestToLine(line: L3): { t: number, s?: number, closest?: V3, closest2?: V3, distance: number } {
		/*
		 line = a + s*b
		 this = c + t*d

		 (this - line) * b = 0
		 (this - line) * d = 0

		 (a + s*b - c - t*d) * b = 0
		 (a + s*b - c - t*d) * d = 0

		 (a - c + s*b - t*d) * b = 0
		 (a - c + s*b - t*d) * d = 0

		 (a - c)*b + (s*b - t*d)*b = 0
		 (a - c)*d + (s*b - t*d)*d = 0

		 (a - c)*b + s*(b*b) - t*(d*b) = 0
		 (a - c)*d + s*(b*d) - t*(d*d) = 0

		 s = (t*(d*b) - (a - c)*b) / (b*b)
		 =>
		 (a - c)*d + (t*(d*b) - (a - c)*b) / (b*b)*(b*d) - t*(d*d) = 0 | * (b*b)
		 (a - c)*d * (b*b) + (t*(d*b) - (a - c)*b)*(b*d) - t*(d*d) * (b*b) = 0
		 (a - c)*d * (b*b) + t*(d*b)*(b*d) - (a - c)*b*(b*d) - t*(d*d) * (b*b) = 0
		 t = ((a - c)*b*(b*d) - (a - c)*d * (b*b)) / ((d*b)*(b*d) - (d*d) * (b*b))
		 */
		if (this.isParallelToLine(line)) {
			return {t: NaN, s: NaN, distance: this.distanceToLine(line)}
		}
		const a = line.anchor, b = line.dir1, c = this.anchor, d = this.dir1
		const bd = b.dot(d), bb = b.squared(), dd = d.squared(), amc = a.minus(c), divisor = bd * bd - dd * bb
		const t = (amc.dot(b) * bd - amc.dot(d) * bb) / divisor
		const s = (amc.dot(b) * dd - amc.dot(d) * bd) / divisor
		return {
			t: t,
			s: s,
			closest: this.at(t),
			closest2: line.at(s),
			distance: this.at(t).distanceTo(line.at(s))
		}
	}

	intersectionWithPlane(plane: P3): V3 {
		// plane: plane.normal1 * p = plane.w
		// line: p=line.point + lambda * line.dir1
		const lambda = (plane.w - plane.normal1.dot(this.anchor)) / plane.normal1.dot(this.dir1)
		const point = this.anchor.plus(this.dir1.times(lambda))
		return point
	}

	tangentAt(t: number): V3 {
		return this.dir1
	}

	intersectWithPlaneLambda(plane): number {
		// plane: plane.normal1 * p = plane.w
		// line: p=line.point + lambda * line.dir1
		const div = plane.normal1.dot(this.dir1)
		if (eq0(div)) return NaN
		const lambda = (plane.w - plane.normal1.dot(this.anchor)) / div
		return lambda
	}

	reversed() {
		return new L3(this.anchor, this.dir1.negated(), -this.tMax, -this.tMin)
	}

	isTsWithPlane(plane: P3) {
		return [this.intersectWithPlaneLambda(plane)]
	}

	flipped() {
		return new L3(this.anchor, this.dir1.negated())
	}

	transform(m4: M4): this {
		const newAnchor = m4.transformPoint(this.anchor)
		const newDir = m4.transformVector(this.dir1)
		return new L3(newAnchor, newDir.unit(), this.tMin * newDir.length(), this.tMax * newDir.length()) as this
	}

	static throughPoints(anchor: V3, b: V3, tMin?: number, tMax?: number): L3 {
		return new L3(anchor, b.minus(anchor).unit(), tMin, tMax)
	}

	static anchorDirection = (anchor: V3, dir: V3): L3 => new L3(anchor, dir.unit())


	static pointTNotNormalized(anchor, dir, x) {
		assertVectors(anchor, dir, x)
		return x.minus(anchor).dot(dir) / dir.squared()
	}

	hashCode(): int {
		return this.anchor.hashCode() * 31 + this.dir1.hashCode()
	}

	static fromPlanes(p1: P3, p2: P3): L3 {
		assertInst(P3, p1, p2)
		const dir = p1.normal1.cross(p2.normal1)
		const length = dir.length()
		if (length < 1e-10) {
			throw new Error('Parallel planes')
		}

		return p1.intersectionWithPlane(p2)
	}

	static readonly X: L3 = new L3(V3.O, V3.X)
	static readonly Y: L3 = new L3(V3.O, V3.Y)
	static readonly Z: L3 = new L3(V3.O, V3.Z)

}
L3.prototype.hlol = Curve.hlol++