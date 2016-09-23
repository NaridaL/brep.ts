/**
 * Created by aval on 20/11/2015.
 */

/**
 * Represents a 3D vector.
 *
 * Example usage:
 *
 *     V3(1, 2, 3)
 *     V3([1, 2, 3])
 *     V3({ x: 1, y: 2, z: 3 })
 *     V3(1, 2) * assumes z=0
 *     V3([1, 2]) // assumes z=0
 *
 * @param a
 * @param b
 * @param {*=} c
 * @constructor
 * @property {number} x
 * @property {number} y
 * @property {number} z
 */
function V3(a, b, c) {
    if (arguments.length == 3) {
        return V3.create(parseFloat(a), parseFloat(b), parseFloat(c))
    } else if (arguments.length == 2) {
        return V3.create(parseFloat(a), parseFloat(b), 0)
    } else if (arguments.length == 1) {
        if (typeof(a) == "object") {
            if (a instanceof V3) {
                // immutable, so
                return a
            //} else if (a instanceof NLA.Vector2D) {
            //    this.x = a.x
            //    this.y = a.y
            //    this.z = 0
            } else if (a instanceof Array || a instanceof Float32Array || a instanceof Float64Array) {
                if (2 == a.length) {
                    return V3.create(parseFloat(a[0]), parseFloat(a[1]), 0)
                } else if (3 == a.length) {
                    return V3.create(parseFloat(a[0]), parseFloat(a[1]), parseFloat(a[2]))
                }
            } else if (('x' in a) && ('y' in a)) {
                if ('z' in a) {
                    return V3.create(parseFloat(a.x), parseFloat(a.y), parseFloat(a.z))
                } else {
                    return V3.create(parseFloat(a.x), parseFloat(a.y), 0)
                }
            }
        }
    }
    throw new Error("invalid arguments" + arguments)
}

V3.id = 0
/**
 * 
 * @param {number} x
 * @param {number} y
 * @param {number} z
 * @returns {V3}
 */
V3.create = function(x, y, z) {
    assert(arguments.length == 3)
    assertNumbers(x, y, z)
    var result = Object.create(V3.prototype, {
        x: { value: x },
        y: { value: y },
        z: { value: z },
        //id: {value: V3.id++}
    } )
    /*
    if ([525, 659].contains(result.id)) {
	    console.log("!!!!!!!!!CREATING V3 "+ result.id, result.hashCodes(), x, y, z, result.$)
    }
    */
    return result
}

V3.prototype = /** @lends V3.prototype */ {
	perturbed: function (v, delta) {
		delta = delta || NLA_PRECISION * 0.8
		return this.map(x => x + (Math.random() - 0.5) * delta)
	},
    e: function (index) {
        assert (index >= 0 && index < 3)
        return index == 0 ? this.x : (index == 1 ? this.y: this.z)
    },

	/**
	 * 
	 * @returns {V3}
	 */
    negated: function() {
        return V3.create(-this.x, -this.y, -this.z)
    },

    abs: function() {
        return V3.create(Math.abs(this.x), Math.abs(this.y), Math.abs(this.z))
    },

    /**
     *
     * @param {V3} a
     * @returns {V3}
     */
    plus: function(a) {
        assertVectors(a)
        return V3.create(this.x + a.x, this.y + a.y, this.z + a.z)
    },

	scale: function (a) {
    	return V3.create(this.x * a.x, this.y * a.y, this.z * a.z)
	},

    /**
     *
     * @param {V3} a
     * @returns {V3}
     */
    minus: function(a) {
        assertVectors(a)
        return V3.create(this.x - a.x, this.y - a.y, this.z - a.z)
    },

	/**
	 *
	 * @param {number} a
	 * @returns {V3}
	 */
    times: function(a) {
        assertNumbers(a)
        return V3.create(this.x * a, this.y * a, this.z * a)
    },

	/**
	 *
	 * @param a
	 * @returns {V3}
	 */
    div: function(a) {
        assertNumbers(a)
        return V3.create(this.x / a, this.y / a, this.z / a)
    },

    /**
     *
     * @param {V3} a
     * @returns {number}
     */
    dot: function(a) {
        assertInst(V3, a)
        return this.x * a.x + this.y * a.y + this.z * a.z
    },

    /**
     * Linearly interpolate
     * @param {V3} a
     * @param {number} t (0 <= t && t <= 1)
     * @returns {V3}
     */
    lerp: function(a, t) {
        assertVectors(a)
        assertNumbers(t)
        assert(0 <= t && t <= 1, "0 <= t && t <= 1")
        return this.plus(a.minus(this).times(t))
    },

    lengthSquared: function() {
        return this.dot(this)
    },

    /**
     *
     * @param {V3} a
     * @returns {number}
     */
    distanceTo: function(a) {
        assertVectors(a)
        //return this.minus(a).length()
	    return Math.hypot(this.x - a.x, this.y - a.y, this.z - a.z)
    },
    distanceToSquared: function(a) {
        assertVectors(a)
        return this.minus(a).lengthSquared()
    },
    // find a vector that is somewhat perpendicular to this one
    nonParallelVector: function() {
        var abs = this.abs()
        if ((abs.x <= abs.y) && (abs.x <= abs.z)) {
            return V3.X
        } else if ((abs.y <= abs.x) && (abs.y <= abs.z)) {
            return V3.Y
        } else {
            return V3.Z
        }
    },
    slerp: function (b, t) {
	    assertVectors(b)
	    NLA.assertNumbers(t)
	    assert(0 <= t && t <= 1, "0 <= t && t <= 1")
	    var sin = Math.sin
	    var omega = this.angleTo(b)
	    return this.times(sin((1 - t) * omega) / sin(omega)).plus(b.times(sin(t * omega) / sin(omega) ))
    },

    min: function(p) {
        return V3.create(
            Math.min(this.x, p.x), Math.min(this.y, p.y), Math.min(this.z, p.z))
    },

    max: function(p) {
        return V3.create(
            Math.max(this.x, p.x), Math.max(this.y, p.y), Math.max(this.z, p.z))
    },
    equals: function(v) {
        return this == v || this.x == v.x && this.y == v.y && this.z == v.z
    },

	/**
	 *
	 * The cross product is defined as:
	 * a x b = |a| * |b| * sin(phi) * n
	 * where |.| is the euclidean norm, phi is the angle between the vectors
	 * and n is a unit vector perpendicular to both a and b.
	 *
	 * The cross product is zero for parallel vectors.
	 *
	 * @param {V3} v
	 * @returns {V3}
	 */
    cross: function(v) {
        return V3.create(
            this.y * v.z - this.z * v.y,
            this.z * v.x - this.x * v.z,
            this.x * v.y - this.y * v.x
        )
    },
    unit: function() {
        return this.div(this.length())
    },
    minElement: function() {
        return Math.min(this.x, this.y, this.z)
    },
    maxElement: function() {
        return Math.max(this.x, this.y, this.z)
    },
    toAngles: function() {
        return {
            theta: Math.atan2(this.z, this.x),
            phi: Math.asin(this.y / this.length())
        }
    },
    toArray: function(n) {
        return [this.x, this.y, this.z].slice(0, n || 3)
    },
	/**
	 * Get a perpendicular vector.
	 * For vectors in the XY-Plane, returns vector rotated CCW.
	 * @returns {V3}
	 */
    getPerpendicular: function() {
        if (NLA.isZero(this.x) && NLA.isZero(this.y)) {
            if (NLA.isZero(this.z)) {
                throw new Error('zero vector');
            }
            // v is Vector(0, 0, v.z)
            return V3.Y
        }
        return V3(-this.y, this.x, 0)
    },
    dim: () => 3,
    els: function () { return [this.x, this.y, this.z] },
    angleXY: function () {
        return Math.atan2(this.y, this.x)
    },
    lengthXY: function () {
	    return Math.hypot(this.x, this.y)
	    //return Math.sqrt(this.x * this.x + this.y * this.y)
    },
	lengthSquaredXY: function () {
        return this.x * this.x + this.y * this.y
	},

    /**
     * Transform this vector element-wise by way of function f. Returns V3(f(x), f(y), f(z))
     * @param f function to apply to elements (number -> number)
     */
    map: function (f) {
        return V3.create(f(this.x), f(this.y), f(this.z))
    },
    toString: function (roundFunction) {
        roundFunction = roundFunction || NLA.defaultRoundFunction
        return "V3(" + [this.x, this.y, this.z].map(roundFunction).join(", ") + ")" //+ this.id
    },

	/**
	 *
	 * @param {V3} v
	 * @returns {number}
	 */
    angleTo: function (v) {
        assert(1 == arguments.length)
        assertVectors(v)
        assert(!this.isZero())
        assert(!v.isZero())
        return Math.acos(this.dot(v) / this.length() / v.length())
    },
	/**
	 *
	 * phi = angle between A and B
	 * alpha = angle between n and normal1
	 *
	 * A . B = ||A|| * ||B|| * cos(phi)
	 * A x B = ||A|| * ||B|| * sin(phi) * n (n = unit vector perpendicular)
	 * (A x B) . normal1 = ||A|| * ||B|| * sin(phi) * cos(alpha)
	 *
	 * @param {V3} vector
	 * @param {V3} normal1
	 * @returns {number}
	 */
    angleRelativeNormal: function (vector, normal1) {
	    assert(2 == arguments.length)
	    NLA.assertVectors (vector, normal1)
	    assert (normal1.hasLength(1))
	    assert (vector.isPerpendicularTo(normal1), "vector.isPerpendicularTo(normal1)" +vector.sce+normal1.sce)
	    assert (this.isPerpendicularTo(normal1), "this.isPerpendicularTo(normal1)" + this.dot(vector))
	    return Math.atan2(this.cross(vector).dot(normal1), this.dot(vector))
    },
    /**
     Returns true iff this is parallel to vector, i.e. this * s == vector, where s is a pos or neg number, using NLA.equals
     Throw a DebugError
     if vector is not a NLA.NLA.Vector or
     if this has a length of 0 or
     if vector has a length of 0
     */
    isParallelTo: function (vector) {
        assertVectors(vector)
        assert(!this.isZero())
        assert(!vector.isZero())
        // a . b takes on values of +|a|*|b| (vectors same direction) to -|a|*|b| (opposite direction)
        // in both cases the vectors are parallel, so check if abs(a . b) == |a|*|b|
        var dot = this.dot(vector)
        return NLA.equals(this.lengthSquared() * vector.lengthSquared(), dot * dot)
    },
    isPerpendicularTo: function (vector) {
        assertVectors(vector)
        assert(!this.isZero(), "!this.isZero()")
        assert(!vector.isZero(), "!vector.isZero()")
        return NLA.isZero(this.dot(vector))
    },
    isReverseDirTo: function (vector) {
	    assertVectors(vector)
	    assert(!this.isZero())
	    assert(!vector.isZero())
	    // a . b takes on values of +|a|*|b| (vectors same direction) to -|a|*|b| (opposite direction)
	    // in both cases the vectors are parallel, so check if abs(a . b) == |a|*|b|
	    var dot = this.dot(vector)
	    return NLA.equals(Math.sqrt(this.lengthSquared() * vector.lengthSquared()), dot)
    },

    /**
     * Returns the length of this NLA.Vector, i.e. the euclidean norm.
     *
     * Note that the partial derivatives of the euclidean norm at point x are equal to the
     * components of the normalized vector x.
     *
     * @returns {number}
     */
    length: function () {
	    return Math.hypot(this.x, this.y, this.z)
	    //return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z)
    },
    /**
     * Definition: V3.isZero == V3.like(V3.ZERO)
     * @returns {boolean}
     */
    isZero: function () {
        return this.like(V3.ZERO)
    },

	/**
	 *
	 * @param {*} obj
	 * @returns {boolean}
	 */
    like: function (obj) {
	    if (obj === this) return true
	    if (!(obj instanceof V3)) return false
	    return NLA.equals(this.x, obj.x) && NLA.equals(this.y, obj.y) && NLA.equals(this.z, obj.z)
    },

    /**
     * equivalent to this.like(v) || this.negated().like(v)
     *
     * @returns {boolean}
     */
    likeOrReversed: function (v) {
	    return NLA.equals(Math.abs(this.dot(v)), Math.sqrt(this.lengthSquared() * v.lengthSquared()))
    },

	/**
	 * Returns a new unit NLA.NLA.Vector (.length() === 1) with the same direction as this vector. Throws a
	 * NLA.DebugError if this has a length of 0.
	 *
	 * @returns {V3}
	 */
    normalized: function () {
	    assert(!this.isZero(), "cannot normalize zero vector")
	    return this.div(this.length())
    },

	/**
	 *
	 * @param {number} newLength
	 * @returns {V3}
	 */
    toLength: function (newLength) {
	    assertNumbers(newLength)
	    return this.times(newLength / this.length())
    },

	/**
	 * See also {@see #setTo} for the individual
	 *
	 * @param {V3} v
	 */
	assign: function (v) {
		assertVectors(v)
		this.x = v.x
		this.y = v.y
		this.z = v.z
	},

	/**
	 * See also {@see #assign} for the V3 version
	 *
	 * @param {number} x
	 * @param {number} y
	 * @param {number=} z
	 */
	setTo: function (x, y, z = 0) {
		this.x = x
		this.y = y
		this.z = z
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
        assertVectors(b)
        // https://en.wikipedia.org/wiki/NLA.Vector_projection#NLA.Vector_projection_2
        return b.times(this.dot(b) / b.dot(b))
    },

	/**
	 *
	 * @param b
	 * @returns {V3}
	 */
    rejectedFrom: function (b) {
        assertVectors(b)
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
    /**
     Returns the sum of the absolute values of the components of this vector.
     E.g. NLA.V(1, -2, 3) === abs(1) + abs(-2) + abs(3) === 1 + 2 + 3 === 6
     */
    absSum: function () {
        return Math.abs(this.x) + Math.abs(this.y) + Math.abs(this.z)
    },
    /**
     * returns max(|x|, |y|, |z|)
     */
    maxAbsElement: function () {
        return Math.max(Math.abs(this.x), Math.abs(this.y), Math.abs(this.z))
    },
    /**
     * returns min(|x|, |y|, |z|)
     */
    minAbsElement: function () {
        return Math.min(Math.abs(this.x), Math.abs(this.y), Math.min(this.z))
    },
    absMaxDim: function () {
        var xAbs = Math.abs(this.x), yAbs = Math.abs(this.y), zAbs = Math.abs(this.z)
        return xAbs >= yAbs ? (xAbs >= zAbs ? 0 : 2) : (yAbs >= zAbs ? 1 : 2)
    },
    absMinDim: function () {
        var xAbs = Math.abs(this.x), yAbs = Math.abs(this.y), zAbs = Math.abs(this.z)
        return xAbs < yAbs ? (xAbs < zAbs ? 0 : 2) : (yAbs < zAbs ? 1 : 2)
    },
    transform: function (m4) {
	    assertInst(M4, m4)
	    return m4.transformPoint(this)
    },
    withElement: function (dim, el) {
	    assert(["x", "y", "z"].includes(dim), ""+dim)
	    assertNumbers(el)
	    if ("x" == dim) {
		    return V3.create(el, this.y, this.z)
	    }
	    if ("y" == dim) {
		    return V3.create(this.x, el, this.z)
	    }
	    return V3.create(this.x, this.y, el)
    },
    hashCode: function() {
	    function floatHashCode(f) {
		    return ~~(f * (1 << 28))
	    }
	    return ~~((floatHashCode(this.x) * 31 + floatHashCode(this.y)) * 31 + floatHashCode(this.z))
    },
    // at most 8 hashcodes
    hashCodes: function() {
	    function floatHashCode(f) {
		    return ~~(f * (1 << 28))
	    }
	    // compare hashCode.floatHashCode
	    // the following ops are equivalent to
	    // floatHashCode((el - NLA_PRECISION) % (2 * NLA_PRECISION))
	    // this results in the hashCode for the (out of 8 possible) cube with the lowest hashCode
	    // the other 7 can be calculated by adding constants
	    var xHC = ~~(this.x * (1 << 28) - 0.5),
		    yHC = ~~(this.y * (1 << 28) - 0.5),
		    zHC = ~~(this.z * (1 << 28) - 0.5),
	    hc = ~~((xHC * 31 + yHC) * 31 + zHC)
	    return [
		    ~~(hc),
		    ~~(hc + 961),
		    ~~(hc + 31),
		    ~~(hc + 31 + 961),

		    ~~(hc + 1),
		    ~~(hc + 1 + 961),
		    ~~(hc + 1 + 31),
		    ~~(hc + 1 + 31 + 961)
	    ]
    },
    compareTo2: function (v3, precision) {
	    precision = precision || 0
	    if (!NLA.equals2(this.x, v3.x, precision)) {
		    return this.x - v3.x
	    } else if (!NLA.equals2(this.y, v3.y, precision)) {
		    return this.y - v3.y
	    } else if (!NLA.equals2(this.z, v3.z, precision)) {
		    return this.z - v3.z
	    } else {
		    return 0
	    }
    },
    compareTo: function (v3) {
	    if (this.x != v3.x) {
		    return this.x - v3.x
	    } else if (this.y != v3.y) {
		    return this.y - v3.y
	    } else {
		    return this.z - v3.z
	    }
    }
}
V3.prototype.toSource = V3.prototype.toString
V3.prototype.constructor = V3

/**
 * @returns {V3}
 */
V3.random = function () {
    return V3.create(Math.random(), Math.random(), Math.random())
}

V3.parallelCondition = function(a, b) {
	return a.dot(b) - a.length() * b.length()
}

/**
 * See http://math.stackexchange.com/questions/44689/how-to-find-a-random-axis-or-unit-vector-in-3d
 * @returns {V3}
 */
V3.randomUnit = function () {
    var zRotation = Math.random() * 2 * Math.PI
    var z = Math.random() * 2 - 1
    var zRadius = Math.sqrt(1 - z * z)
    return V3(zRadius * Math.cos(zRotation), zRadius * Math.sin(zRotation), z)
}
V3.fromAngles = function(theta, phi) {
    return V3.create(Math.cos(theta) * Math.cos(phi), Math.sin(phi), Math.sin(theta) * Math.cos(phi))
}
V3.fromFunction = function (f) {
    return V3.create(f(0), f(1), f(2))
}
V3.min = function(a, b) {
    return V3.create(Math.min(a.x, b.x), Math.min(a.y, b.y), Math.min(a.z, b.z))
}
V3.max = function(a, b) {
    return V3.create(Math.max(a.x, b.x), Math.max(a.y, b.y), Math.max(a.z, b.z))
}
V3.lerp = function(a, b, fraction) {
    return b.minus(a).times(fraction).plus(a)
}
V3.fromArray = function(a) {
    return V3.create(a[0], a[1], a[2])
}
V3.angleBetween = function(a, b) {
    return a.angleTo(b)
}
V3.zip = function (f, ...args) {
	assert(f instanceof Function)
	return V3.create(
		f.apply(undefined, args.map(x => x.x)),
		f.apply(undefined, args.map(x => x.y)),
		f.apply(undefined, args.map(x => x.z))
	)
}
V3.normalOnPoints = function (v0, v1, v2) {
	assertVectors(v0, v1, v2)
	return v1.minus(v0).cross(v2.minus(v0))
}
V3.add = function (...vs) {
	assertVectors.apply(undefined, vs)
	var x = 0, y = 0, z = 0
	var i = vs.length
	while (i--) {
		x += vs[i].x
		y += vs[i].y
		z += vs[i].z
	}
	return V3.create(x, y, z)
}
V3.sub = function (...vs) {
	assertVectors.apply(undefined, vs)
	var x = vs[0].x, y = vs[0].y, z = vs[0].z
	for (var i = 1; i < vs.length; i++) {
		x -= vs[i].x
		y -= vs[i].y
		z -= vs[i].z
	}
	return V3.create(x, y, z)
}
V3.flattenV3Array = function (v3arr, dest, srcStart, destStart, v3count) {
	//assert (v3arr.every(v3 => v3 instanceof V3), "v3arr.every(v3 => v3 instanceof V3)")
	srcStart = srcStart || 0
	destStart = destStart || 0
	v3count = v3count || (v3arr.length - srcStart)
	dest = dest || new Float32Array(3 * v3count)
	assert(dest.length - destStart >= v3count, "dest.length - destStart >= v3count")
	for (var i = 0; i < v3count; i++) {
		var v = v3arr[srcStart + i]
		dest[destStart + i * 3]     = v.x
		dest[destStart + i * 3 + 1] = v.y
		dest[destStart + i * 3 + 2] = v.z
	}
	return dest
}
V3.perturbed = function (v, delta) {
	return v.perturbed(delta)
}
V3.polar = function (radius, phi) {
	return V3.create(radius * Math.cos(phi), radius * Math.sin(phi), 0)
}
V3.areDisjoint = function (it) {
	var vSet = new NLA.CustomSet, v
	while (v = it.next().value) {
		if (!v.equals(vSet.canonicalizeLike(v))) {
			// like value already in set
			return false
		}
	}
	return true
}
/** @const */ V3.ZERO = V3.create(0, 0, 0)
/** @const */ V3.ONES = V3.create(1, 1, 1)
/** @const */ V3.X = V3.create(1, 0, 0)
/** @const */ V3.Y = V3.create(0, 1, 0)
/** @const */ V3.Z = V3.create(0, 0, 1)
/** @const */ V3.INF = V3.create(Infinity, Infinity, Infinity)
/** @const */ V3.XYZ = [V3.X, V3.Y, V3.Z]
