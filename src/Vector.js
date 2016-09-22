
/**
 *
 * @param v
 * @constructor
 * @property {Float64Array} v
 */
NLA.Vector = function (v) {
	if (NLA_DEBUG && !(v instanceof Float64Array)) {
		throw new Error("!!");
	}
	this.v = v
}
NLA.Vector.fromFunction = function(dim, f) {
	NLA.assertNumbers(dim);
	var e = new Float64Array(dim)
	var i = dim;
	while (i--) {
		e[i] = f(i)
	}
	return new NLA.Vector(e);
}
NLA.Vector.random = function (dim) {
	return NLA.Vector.fromFunction(dim, (i) => Math.random())
}
NLA.Vector.fromArguments = function (...args) {
	assert (args[0] instanceof Float64Array || args.every(a => 'number' == typeof a),
		"args[0] instanceof Float64Array || args.every(a => 'number' == typeof a)")
	return new NLA.Vector(args[0] instanceof Float64Array ? args[0] : Float64Array.from(args))
}
NLA.V = NLA.Vector.fromArguments
NLA.Vector.prototype = {
	dim: function () { return this.v.length },
	e: function (index) {
		if (0 > index || index >= this.v.length) {
			throw new Error("array index out of bounds")
		}
		return this.v[index]
	},
	plus: function (vector) {
		var u = this.v, v = vector.v
		var n = new Float64Array(u.length)
		var i = u.length
		while (i--) {
			n[i] = u[i] + v[i]
		}
		return new NLA.Vector(n)
	},
	minus: function (/** V3 */ vector) {
		var u = this.v, v = vector.v
		var n = new Float64Array(u.length)
		var i = u.length
		while (i--) {
			n[i] = u[i] - v[i]
		}
		return new NLA.Vector(n)
	},
	div: function (val) {
		var u = this.v
		var n = new Float64Array(u.length)
		var i = u.length
		while (i--) {
			n[i] = u[i] / val
		}
		return new NLA.Vector(n)
	},
	dot: function (vector) {
		assert(this.dim == vector.dim, "passed vector must have the same dim")
		var result = 0;
		var u = this.v, v = vector.v
		var i = u.length
		while (i--) {
			result += u[i] * v[i]
		}
		return result
	},
	times: function (val) {
		var u = this.v
		var n = new Float64Array(u.length)
		var i = u.length
		while (i--) {
			n[i] = u[i] * val
		}
		return new NLA.Vector(n)
	},
	cross: function(vector) {
		assertInst(NLA.Vector, vector)
		var n = new Float64Array(3)
		n[0] = this.v[1] * vector.v[2] - this.v[2] * vector.v[1]
		n[1] = this.v[2] * vector.v[0] - this.v[0] * vector.v[2]
		n[2] = this.v[0] * vector.v[1] - this.v[1] * vector.v[0]

		return new NLA.Vector(n)
	},
	equals: function (obj) {
		if (obj === this) return true
		if (!(obj.constructor != NLA.Vector)) return false
		if (this.v.length != obj.v.length) return false
		var i = this.v.length
		while (i--) {
			if (!NLA.equals(this.v[i], obj.v[i])) return false
		}
		return true
	},
	map: function(f) {
		var e = new Float64Array(this.v.length)
		var i = e.length
		while (i--) {
			e[i] = f(i)
		}
		return new NLA.Vector(e);
	},
	toString: function (roundFunction) {
		roundFunction = roundFunction || ((v) => +v.toFixed(6))
		return "NLA.Vector(" + this.v.map(roundFunction).join(", ") + ")"
	},
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
	angleTo: function (vector) {
		assertInst(NLA.Vector, vector)
		assert(!this.isZero(), "!this.isZero()")
		assert(!vector.isZero(), "!vector.isZero()")
		return Math.acos(this.dot(vector) / this.length() / vector.length())
	},
	/**
	 Returns true iff this is parallel to vector, using NLA.equals
	 Throw a DebugError
	 if vector is not a NLA.NLA.Vector or
	 if this has a length of 0 or
	 if vector has a length of 0
	 */
	isParallelTo: function (vector) {
		assertInst(NLA.Vector, vector)
		assert(!this.isZero(), "!this.isZero()")
		assert(!vector.isZero(), "!vector.isZero()")
		// a . b takes on values of +|a|*|b| (vectors same direction) to -|a|*|b| (opposite direction)
		// in both cases the vectors are paralle, so check if abs(a . b) == |a|*|b|
		return NLA.equals(Math.sqrt(this.lengthSquared() * vector.lengthSquared()), Math.abs(this.dot(vector)))
	},
	isPerpendicularTo: function (vector) {
		assertInst(NLA.Vector, vector)
		assert(!this.isZero(), "!this.isZero()")
		assert(!vector.isZero(), "!vector.isZero()")
		return NLA.isZero(this.dot(vector))
	},
	/**
	 Returns true iff the length of this vector is 0, as returned by NLA.isZero.
	 Definition: NLA.NLA.Vector.prototype.isZero = () => NLA.isZero(this.length())
	 */
	isZero: function () { return NLA.isZero(this.length()) },
	/**
	 Returns the length of this NLA.NLA.Vector, i.e. the euclidian norm.
	 */
	length: function () {
		return Math.sqrt(this.lengthSquared())
	},
	lengthSquared: function () {
		var result = 0
		var u = this.v
		var i = u.length
		while (i--) {
			result += u[i] * u[i]
		}
		return result
	},
	/**
	 Returns a new unit NLA.NLA.Vector (.length() === 1) with the same direction as this vector.
	 Throws a NLA_DEBUGError if this has a length of 0.
	 */
	normalized: function () {
		var length = this.length()
		if (NLA.isZero(length)) {
			throw new Error("cannot normalize zero vector")
		}
		return this.div(this.length())
	},
	asRowMatrix: function () {
		return new NLA.Matrix(this.v.length, 1, this.v)
	},
	asColMatrix: function () {
		return new NLA.Matrix(1, this.v.length, this.v)
	},
	/**
	 Returns a new NLA.NLA.Vector which is the projection of this vector onto the passed vector.
	 Examples
	 NLA.V(3, 4).projectedOn(NLA.V(1, 0)) // returns NLA.V(3, 0)
	 NLA.V(3, 4).projectedOn(NLA.V(2, 0)) // returns NLA.V(3, 0)
	 NLA.V(3, 4).projectedOn(NLA.V(-1, 0)) // returns NLA.V(-3, 0)
	 NLA.V(3, 4).projectedOn(NLA.V(0, 1)) // returns NLA.V(0, 4)
	 NLA.V(3, 4).projectedOn(NLA.V(1, 1)) // returns
	 */
	projectedOn: function (b) {
		assertInst(NLA.Vector, b)
		// https://en.wikipedia.org/wiki/NLA.Vector_projection#NLA.Vector_projection_2
		return b.times(this.dot(b) / b.dot(b))
	},
	rejectedOn: function (b) {
		assertInst(NLA.Vector, b)
		// https://en.wikipedia.org/wiki/NLA.Vector_projection#NLA.Vector_projection_2
		return this.minus(b.times(this.dot(b) / b.dot(b)))
	},
	/**
	 Returns true iff the length() of this vector is equal to "length", using NLA.equals
	 E.g. NLA.V(3, 4).hasLength(5) === true
	 NLA.V(1, 1).hasLength(1) === false
	 */
	hasLength: function (length) {
		NLA.assertNumbers(length)
		return NLA.equals(length, this.length())
	},
	V3: function () {
		//assert(this.dim() == 3)
		return V3.create(this.v[0], this.v[1], this.v[2])
	}
}
NLA.Vector.Zero = function (dim) {
	NLA.assertNumbers(dim)
	var i = 0
	var n = new Float64Array(dim)
	while (i--) {
		n[i] = 0
	}
	return new NLA.Vector(n)
}
NLA.Vector.Unit = function (dim, dir) {
	NLA.assertNumbers(dim, dir)
	var i = 0
	var n = new Float64Array(dim)
	while (i--) {
		n[i] = +(i == dir) // +true === 1, +false === 0
	}
	return new NLA.Vector(n)
}
