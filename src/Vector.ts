class Vector {
	constructor(readonly v: Float64Array) {
		assertInst(Float64Array, v)
	}

	static fromFunction(dims: int, f: (dim: int) => number): Vector {
		assertNumbers(dims)
		const e = new Float64Array(dims)
		let i = dims
		while (i--) {
			e[i] = f(i)
		}
		return new Vector(e)
	}

	static random(dims: int): Vector {
		return Vector.fromFunction(dims, (i) => Math.random())
	}

	static from(...args): Vector {
		assert(args[0] instanceof Float64Array || args.every(a => 'number' == typeof a),
			'args[0] instanceof Float64Array || args.every(a => "number" == typeof a)')
		return new Vector(args[0] instanceof Float64Array ? args[0] : Float64Array.from(args))
	}

	dim(): int {
		return this.v.length
	}

	e(index: int): number {
		if (0 > index || index >= this.v.length) {
			throw new Error('array index out of bounds')
		}
		return this.v[index]
	}

	plus(vector: Vector): Vector {
		const u = this.v, v = vector.v
		const n = new Float64Array(u.length)
		let i = u.length
		while (i--) {
			n[i] = u[i] + v[i]
		}
		return new Vector(n)
	}

	minus(vector: Vector): Vector {
		const u = this.v, v = vector.v
		const n = new Float64Array(u.length)
		let i = u.length
		while (i--) {
			n[i] = u[i] - v[i]
		}
		return new Vector(n)
	}

	times(factor: number): Vector {
		const u = this.v
		const n = new Float64Array(u.length)
		let i = u.length
		while (i--) {
			n[i] = u[i] * factor
		}
		return new Vector(n)
	}

	div(val: number): Vector {
		const u = this.v
		const n = new Float64Array(u.length)
		let i = u.length
		while (i--) {
			n[i] = u[i] / val
		}
		return new Vector(n)
	}

	dot(vector: Vector): number {
		assert(this.dim == vector.dim, 'passed vector must have the same dim')
		let result = 0
		const u = this.v, v = vector.v
		let i = u.length
		while (i--) {
			result += u[i] * v[i]
		}
		return result
	}

	cross(vector: Vector): Vector {
		assertInst(Vector, vector)
		const n = new Float64Array(3)
		n[0] = this.v[1] * vector.v[2] - this.v[2] * vector.v[1]
		n[1] = this.v[2] * vector.v[0] - this.v[0] * vector.v[2]
		n[2] = this.v[0] * vector.v[1] - this.v[1] * vector.v[0]

		return new Vector(n)
	}

	equals(obj: any): boolean {
		if (obj === this) return true
		if (obj.constructor !== Vector) return false
		if (this.v.length != obj.v.length) return false
		let i = this.v.length
		while (i--) {
			if (!eq(this.v[i], obj.v[i])) return false
		}
		return true
	}

	map(f: (el: number, dim: number) => number): Vector {
		return new Vector(this.v.map(f))
	}

	toString(roundFunction?: (x: number) => any): string {
		roundFunction = roundFunction || ((v) => +v.toFixed(6))
		return 'Vector(' + this.v.map(roundFunction).join(', ') + ')'
	}

	/*
	 get x() {
	 return this.v[0]
	 },
	 get y() {
	 return this.v[1]
	 },
	 get z() {
	 return this.v[2]
	 },
	 get w() {
	 return this.v[3]
	 },
	 */
	angleTo(vector: Vector): number {
		assertInst(Vector, vector)
		assert(!this.isZero(), '!this.isZero()')
		assert(!vector.isZero(), '!vector.isZero()')
		return Math.acos(this.dot(vector) / this.length() / vector.length())
	}

	/**
	 Returns true iff this is parallel to vector, using equals
	 Throw a DebugError
	 if vector is not a Vector or
	 if this has a length of 0 or
	 if vector has a length of 0
	 */
	isParallelTo(vector: Vector): boolean {
		assertInst(Vector, vector)
		assert(!this.isZero(), '!this.isZero()')
		assert(!vector.isZero(), '!vector.isZero()')
		// a . b takes on values of +|a|*|b| (vectors same direction) to -|a|*|b| (opposite direction)
		// in both cases the vectors are paralle, so check if abs(a . b) == |a|*|b|
		return eq(Math.sqrt(this.lengthSquared() * vector.lengthSquared()), Math.abs(this.dot(vector)))
	}

	isPerpendicularTo(vector: Vector): boolean {
		assertInst(Vector, vector)
		assert(!this.isZero(), '!this.isZero()')
		assert(!vector.isZero(), '!vector.isZero()')
		return eq0(this.dot(vector))
	}

	/**
	 Returns true iff the length of this vector is 0, as returned by NLA.isZero.
	 Definition: Vector.prototype.isZero = () => NLA.isZero(this.length())
	 */
	isZero(): boolean {
		return eq0(this.length())
	}

	/*/ Returns the length of this Vector, i.e. the euclidian norm.*/
	length(): number {
		return Math.hypot.apply(undefined, this.v)
		//return Math.sqrt(this.lengthSquared())
	}

	lengthSquared(): number {
		let result = 0
		let u = this.v
		let i = u.length
		while (i--) {
			result += u[i] * u[i]
		}
		return result
	}

	// Returns a new unit Vector (.length() === 1) with the same direction as this vector. Throws a
	// NLA_DEBUGError if this has a length of 0.
	normalized(): Vector {
		const length = this.length()
		if (eq0(length)) {
			throw new Error('cannot normalize zero vector')
		}
		return this.div(this.length())
	}

	asRowMatrix(): Matrix {
		return new Matrix(this.v.length, 1, this.v)
	}

	asColMatrix(): Matrix {
		return new Matrix(1, this.v.length, this.v)
	}

	/**
	 Returns a new Vector which is the projection of this vector onto the passed vector.
	 Examples
	 NLA.V(3, 4).projectedOn(NLA.V(1, 0)) // returns NLA.V(3, 0)
	 NLA.V(3, 4).projectedOn(NLA.V(2, 0)) // returns NLA.V(3, 0)
	 NLA.V(3, 4).projectedOn(NLA.V(-1, 0)) // returns NLA.V(-3, 0)
	 NLA.V(3, 4).projectedOn(NLA.V(0, 1)) // returns NLA.V(0, 4)
	 NLA.V(3, 4).projectedOn(NLA.V(1, 1)) // returns
	 */
	projectedOn(b: Vector): Vector {
		assertInst(Vector, b)
		// https://en.wikipedia.org/wiki/Vector_projection#Vector_projection_2
		return b.times(this.dot(b) / b.dot(b))
	}

	rejectedOn(b: Vector): Vector {
		assertInst(Vector, b)
		// https://en.wikipedia.org/wiki/Vector_projection#Vector_projection_2
		return this.minus(b.times(this.dot(b) / b.dot(b)))
	}

	/**
	 Returns true iff the length() of this vector is equal to 'length', using equals
	 E.g. NLA.V(3, 4).hasLength(5) === true
	 NLA.V(1, 1).hasLength(1) === false
	 */
	hasLength(length: number): boolean {
		assertNumbers(length)
		return eq(length, this.length())
	}

	V3(): V3 {
		//assert(this.dim() == 3)
		return new V3(this.v[0], this.v[1], this.v[2])
	}

	static Zero(dims: int): Vector {
		assertNumbers(dims)
		let i = 0
		const n = new Float64Array(dims)
		while (i--) {
			n[i] = 0
		}
		return new Vector(n)
	}

	static Unit(dims: int, dir: int): Vector {
		assertNumbers(dims, dir)
		let i = 0
		const n = new Float64Array(dims)
		while (i--) {
			n[i] = +(i == dir) // +true === 1, +false === 0
		}
		return new Vector(n)
	}
}