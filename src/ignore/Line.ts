class Line {

	static throughPoints = (anchor, b) => new Line(anchor, b.minus(anchor).normalized())
	static anchorDirection = (anchor, direction) => new Line(anchor, direction.unit())
	dir1: Vector
	anchor: Vector

	/**
	 * n-dimensional line
	 *
	 * @param anchor
	 * @param dir1
	 * @constructor
	 */
	constructor(anchor, dir1) {
		assertVectors(anchor, dir1)
		assert(dir1.hasLength(1), 'dir must be unit')
		this.anchor = anchor
		this.dir1 = dir1
	}

	containsPoint() {
	}

	distanceToPoint(x) {
		assertVectors(x)
		// See http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
		var t = x.minus(this.anchor).dot(this.dir1)
		return this.at(t).minus(x).length()

		//return x.minus(this.anchor).cross(x.minus(this.anchor.plus(this.dir1))).length()
	}

	/**
	 * Every point x on this line is described by the equation x = this.anchor + lambda * this.dir1
	 * This function returns lambda for a given point x
	 * @param x
	 */
	pointT(x) {
		assertVectors(x)
		var t = x.minus(this.anchor).dot(this.dir1)
		return t
	}

	isParallelToLine(line) {
		assertInst(Line, line)
		// we know that 1 == this.dir1.length() == line.dir1.length(), we can check for parallelity simpler than
		// isParallelTo()
		return eq(1, Math.abs(this.dir1.dot(line.dir1)))
	}

	angleToLine(line) {
		assertInst(Line, line)
		return this.dir1.angleTo(line.dir1)
	}

	at(t: number): V3 {
		assertNumbers(t)
		return this.anchor.plus(this.dir1.times(t))
	}
}
