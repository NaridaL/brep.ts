
"use strict";
var NLA = {}

NLA.assertVectors = function () {
	if (NLA.DEBUG) {
		for (var i = 0; i < arguments.length; i++) {
			if (!(arguments[i] instanceof NLA.Vector3 || arguments[i] instanceof NLA.Vector)) {
				// Arrays.prototype.slice.call is inefficient, but it doesn't matter here
				throw new Error("NLA.assertVectors arguments[" + (i) + "] is not a vector. " + typeof arguments[i] + " == typeof " + arguments[i]);
			}
		}
	}
	return true
}
NLA.PRECISION = 1 / (1 << 27)
console.log("NLA.PRECISION", NLA.PRECISION)
/**
 *
 * @type {boolean}
 * @const
 * @define {boolean}
 */
NLA.DEBUG = true
NLA.isZero = (x) => Math.abs(x) < NLA.PRECISION
NLA.isZero2 = (x, precision) => Math.abs(x) < precision
NLA.equals = (x, y) => Math.abs(x - y) <= NLA.PRECISION
NLA.lt = (x, y) => x + NLA.PRECISION < y
NLA.gt = (x, y) => x > y + NLA.PRECISION
NLA.le = (x, y) => x <= y + NLA.PRECISION
NLA.ge = (x, y) => x + NLA.PRECISION >= y
NLA.equals2 = (x, y, precision) => Math.abs(x - y) < precision
NLA.eqAngle = (x, y) => NLA.zeroAngle(x - y)
NLA.zeroAngle = (x) => ((x % (2 * Math.PI)) + 2 * Math.PI + NLA.PRECISION) % (2 * Math.PI) < 2 * NLA.PRECISION
NLA.snapTo = (x, to) => Math.abs(x - to) < NLA.PRECISION ? to : x
NLA.canonAngle = (x) => ((x % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI)
var COLORS = {
	RD_FILL:0x9EDBF9,
	RD_STROKE:0x77B0E0,
	TS_FILL: 0xD19FE3,
	TS_STROKE: 0xA76BC2,
	PP_FILL:0xF3B6CF,
	PP_STROKE:0xEB81B4,
}
function deg2rad(angle) {
	//  discuss at: http://phpjs.org/functions/deg2rad/
	// original by: Enrique Gonzalez
	// improved by: Thomas Grainger (http://graingert.co.uk)
	//   example 1: deg2rad(45);
	//   returns 1: 0.7853981633974483

	return angle * .017453292519943295; // (angle / 180) * Math.PI;
}
function rad2deg(rad) {
	//  discuss at: http://phpjs.org/functions/deg2rad/
	// original by: Enrique Gonzalez
	// improved by: Thomas Grainger (http://graingert.co.uk)
	//   example 1: deg2rad(45);
	//   returns 1: 0.7853981633974483

	return rad / .017453292519943295; // (angle / 180) * Math.PI;
}
NLA.repeatString = function(count, str) {
	if (count == 0) {
		return ""
	}
	count *= str.length
	var halfCharLength = count / 2
	var result = str

	// double the input until it is long enough.
	while (result.length <= halfCharLength) {
		result += result
	}
	// use substring to hit the precise length target without
	// using extra memory
	return result + result.substring(0, count - result.length)
}
NLA.assertNumbers = function () {
    if (NLA.DEBUG) {
        for (var i = 0; i < arguments.length; i++) {
            if (typeof arguments[i] !== 'number') {
                // Arrays.prototype.slice.call is inefficient, but it doesn't matter here
                throw new Error("NLA.assertNumbers arguments[" + (i) + "] is not a number. " + typeof arguments[i] + " == typeof " + arguments[i]);
            }
        }
    }
	return true
}
NLA.assert = function (value, message) {
	if (NLA.DEBUG && !value) {
		throw new Error("NLA.assert failed: " + message);
	}
	return true
}
NLA.mod = function (a, b) {
	return ((a % b) + b) % b
}
NLA.arraySwap = function (arr, i, j) {
    var temp = arr[i]
    arr[i] = arr[j]
    arr[j] = temp
}
NLA.arrayCopy = function(src,sstart,dst,dstart,length) {
    dstart += length;
    length += sstart;
    while(length-- > sstart) {
        dst[--dstart] = src[length];
	}       
}
NLA.clamp = function (val, min, max) {
	NLA.assertNumbers(val, min, max)
	return Math.max(min, Math.min(max, val))
}
NLA.CustomSet = function (iterable) {
	this._map = new Map()
	this._size = 0
	if (iterable) {
		this.addAll(iterable)
	}
}
NLA.CustomSet.prototype = {
	add: function (val) {
		// you can't use this.canonicalize here, as there is no way to differentiate if val
		// is new or if val was === the exisitng value (not only .equals)
		var hashCode = val.hashCode(), bucket = this._map.get(hashCode)
		if (bucket) {
			if (bucket.some(x => x.equals(val))) { return false }
			bucket.push(val)
		} else {
			this._map.set(hashCode, [val])
		}
		this._size++
		return true
	},
	addAll: function (iterable) {
		for (var val of iterable) {
			this.add(val)
		}
	},
	canonicalize: function (val) {
		var hashCode = val.hashCode(), bucket = this._map.get(hashCode)
		if (bucket) {
			var existing = bucket.find(x => x.equals(val))
			if (existing) { return existing }
			bucket.push(val)
		} else {
			this._map.set(hashCode, [val])
		}
		this._size++
		return val
	},
	has: function (val) {
		var hashCode = val.hashCode(), bucket = this._map.get(hashCode)
		return bucket && bucket.some(x => x.equals(val))
	},
	getLike: function (val) {
		for (var hashCode of val.hashCodes()) {
			var bucket = this._map.get(hashCode)
			var canonVal = bucket && bucket.find(x => x.like(val))
			if (canonVal) return canonVal
		}
	},
	canonicalizeLike: function (val) {
		// if this.getLike(val) is defined, return it, otherwise add val and return val
		return this.getLike(val) || this.canonicalize(val)
	},
	addLike: function (val) {
		return !this.getLike(val) && this.add(val)
	},
	delete: function (val) {
		var hashCode = val.hashCode(), bucket = this._map.get(hashCode)
		if (bucket) {
			var index = bucket.findIndex(x => x.equals(val))
			if (-1 != index) {
				if (1 == bucket.size) {
					this._map.delete(bucket)
				} else {
					bucket.splice(index, 1)
				}
				this._size--
				return true
			}
		}
		return false
	},
	deleteLike: function (val) {
		for (var hashCode of val.hashCodes()) {
			var bucket = this._map.get(hashCode)
			if (bucket) {
				var index = bucket.findIndex(x => x.like(val))
				if (-1 != index) {
					var deleted = bucket[index]
					if (1 == bucket.size) {
						this._map.delete(bucket)
					} else {
						bucket.splice(index, 1)
					}
					this._size--
					return deleted
				}
			}
		}
	},
	entries: function* () {
		for (var bucket of this._map) {
			yield* bucket
		}
	},
	clear: function () {
		this._map.clear()
		this._size = 0
	},
	get size() {
		return this._size
	}
}
NLA.CustomMap = function () {

}
NLA.CustomMap.prototype = {
	set: function (key, val) {
		var hashCode = key.hashCode(), bucket = this._map.get(hashCode)
		if (bucket) {
			var pairIndex = bucket.findIndex(pair => pair.key.equals(key))
			if (-1 == pairIndex) {
				bucket.push({key: key, value: val})
			} else {
				bucket[pairIndex].value = val
				return false
			}
		} else {
			this._map.set(hashCode, [{key: key, value: val}])
		}
		this._size++
		return true
	},
	has: function (key) {
		var hashCode = key.hashCode(), bucket = this._map.get(hashCode)
		return bucket && bucket.some(x => x.key.equals(key))
	},
	getLike: function (val) {
		for (var hashCode of val.hashCodes()) {
			var bucket = this._map.get(hashCode)
			var canonVal = bucket && bucket.find(x => x.key.like(val))
			if (canonVal) return canonVal
		}
	},
	setLike: function (key, val) {
		return !this.getLike(val) && this.set(key, val)
	},
	delete: function (key) {
		var hashCode = key.hashCode(), bucket = this._map.get(hashCode)
		if (bucket) {
			var index = bucket.findIndex(x => x.key.equals(key))
			if (-1 != index) {
				if (1 == bucket.size) {
					this._map.delete(bucket)
				} else {
					bucket.splice(index, 1)
				}
				this._size--
				return true
			}
		}
		return false
	},
	deleteLike: function (key) {
		for (var hashCode of key.hashCodes()) {
			var bucket = this._map.get(hashCode)
			if (bucket) {
				var index = bucket.findIndex(x => x.key.like(key))
				if (-1 != index) {
					var deleted = bucket[index]
					if (1 == bucket.size) {
						this._map.delete(bucket)
					} else {
						bucket.splice(index, 1)
					}
					this._size--
					return deleted
				}
			}
		}
	},
	entries: function* () {
		for (var bucket of this._map) {
			yield* bucket
		}
	},
	clear: function () {
		this._map.clear()
		this._size = 0
	},
	get size() {
		return this._size
	}
}
NLA.CustomSet.prototype.values = NLA.CustomSet.prototype.entries
NLA.CustomSet.prototype[Symbol.iterator] = NLA.CustomSet.prototype.entries
NLA.CustomSet.prototype.keys = NLA.CustomSet.prototype.entries
NLA.randomColor = function () {
	return Math.floor(Math.random() * 0x1000000)
}
NLA.mapAdd = function (map, key, val) {
	var list = map.get(key)
	if (list) {
		list.push(val)
	} else {
		map.set(key, [val])
	}
}
// Add several convenience methods to the classes that support a transform() method:
NLA.addTransformationMethods = function (obj) {
	NLA.addOwnProperties(obj, NLA.tranformablePrototype)
}
String.prototype.capitalizeFirstLetter = function() {
	return this.charAt(0).toUpperCase() + this.slice(1);
}
var ARRAY_UTILITIES = {
	pushAll: function (arr) {
		Array.prototype.push.apply(this, arr)
	},
	copyStep: function(src,sstart,sstep, dst,dstart,dstep,count) {
		var srcIndex = sstart + count * sstep;
		var dIndex = dstart + count * dstep;
		while(srcIndex > sstart) {
			dst[dIndex -= dstep] = src[srcIndex -= sstep];
		}
	},
	sliceStep: function (start, step, chunkSize) {
		NLA.assertNumbers(start, step)
		chunkSize = chunkSize || 1;
		var result = new Array(Math.ceil((this.length - start) / step)); // "- start" so that chunk in the last row will also be selected, even if the row is not complete
		var index = 0;
		for (var i = 0; i < this.length; i += step) {
			for (var j = 0; j < chunkSize; j++) {
				result[index++] = this[start + i + j];
			}
		}
		return result;
	},
	concatenated: function () {
		return Array.prototype.concat.apply([], this);
	},
	min: function() {
		return Math.min.apply(null, this);
	},
	max: function() {
		return Math.max.apply(null, this);
	},
	indexWithMax: function (f) {
		var i = this.length, result = -1, maxVal = -Infinity
		while (i--) {
			var val = f(this[i])
			if (val > maxVal) {
				maxVal = val
				result = i
			}
		}
		return result
	},
	withMax: function (f) {
		var i = this.length, result = undefined, maxVal = -Infinity
		while (i--) {
			var el = this[i], val = f(el)
			if (val > maxVal) {
				maxVal = val
				result = el
			}
		}
		return result
	},
	/**
	 Returns the sum of the absolute values of the components of this vector.
	 E.g. NLA.V(1, -2, 3) === abs(1) + abs(-2) + abs(3) === 1 + 2 + 3 === 6
	 */
	absSum: function() {
		var i = this.length
		var result = 0
		while (i--) {
			result += Math.abs(this[i])
		}
		return result
	},
	sum: function () {
		var i = this.length
		var result = 0
		while (i--) {
			result += this[i]
		}
		return result
	},
	/*
	contains: function (o) {
		return -1 != this.indexOf(o)
	},
	*/
	isEmpty: function () {
		return 0 == this.length
	},
	unique: function () {
		var uniqueSet = new Set(this)
		return Array.from(uniqueSet)
	},
	remove: function (o) {
		var index = this.indexOf(o);
		if (index != -1) {
			this.splice(index, 1);
		}
	},
	removeAll: function (o) {
		var i = o.length
		while (i--) {
			this.remove(o[i]);
		}
	},
	toggle: function (o) {
		var index = this.indexOf(o);
		if (index != -1) {
			this.splice(index, 1);
		} else {
			this.push(o)
		}
	},
	binaryIndexOf: function(searchElement, cmp) {

		var minIndex = 0;
		var maxIndex = this.length - 1;
		var currentIndex;
		var currentElement;

		while (minIndex <= maxIndex) {
			currentIndex = (minIndex + maxIndex) / 2 | 0;
			currentElement = this[currentIndex];

			if (cmp(currentElement, searchElement) < 0) {
				minIndex = currentIndex + 1;
			}
			else if (cmp(currentElement, searchElement) > 0) {
				maxIndex = currentIndex - 1;
			}
			else {
				return currentIndex;
			}
		}

		return -1;
	},
	binaryInsert: function (el, cmp) {
		cmp = cmp || NLA.minus
		var minIndex = 0
		var maxIndex = this.length
		var currentIndex
		var currentElement

		while (minIndex < maxIndex) {
			currentIndex = ~~((minIndex + maxIndex) / 2)
			currentElement = this[currentIndex]

			if (cmp(currentElement, el) < 0) {
				minIndex = currentIndex + 1
			} else {
				maxIndex = currentIndex
			}
		}

		this.splice(minIndex, 0, el)
	}
}
NLA.minus = (a, b) => a - b
for (var key in ARRAY_UTILITIES) {
	NLA["array" + key.capitalizeFirstLetter()] = function (arr, ...rest) {
		ARRAY_UTILITIES[key].apply(arr, rest)
	}
}
NLA.arrayCopyStep = function(src,sstart,sstep, dst,dstart,dstep,count) {
	var srcIndex = sstart + count * sstep;
	var dIndex = dstart + count * dstep;
	while(srcIndex > sstart) {
		dst[dIndex -= dstep] = src[srcIndex -= sstep];
	}
}
NLA.arrayCopyBlocks = function (src,sstart,sstep, dst,dstart,dstep,blockSize, blockCount) {
	for (var i = 0; i < blockCount; i++) {
		NLA.arrayCopy(src, sstart + sstep * i, dst, dstart + dstep * i, blockSize)
	}
}
NLA.arrayRange = function (start, end, step) {
	NLA.assertNumbers(start, step)
	console.log(Math.ceil((end - start) / step))
	var result = new Array(Math.ceil((end - start) / step)); // "- start" so that chunk in the last row will also be selected, even if the row is not complete
	var index = 0
	for (var i = 0; i < end; i += step) {
		result[index++] = i
	}
	return result;
}
NLA.arrayFromFunction = function(length, f) {
	NLA.assertNumbers(length)
	NLA.assert("function" == typeof f)
	var a = new Array(length)
	var elIndex = length;
	while (elIndex--) {
		a[elIndex] = f(elIndex);
	}
	return a
}
function arrayEquals(arr1, arr2) {
	if (arr1.length != arr2.length) return false;
	
	var i = arr1.length
	while (i--) {
		//if (arr1[i] != (arr2.result = Math.max(result, arr[i])))
	}
	return result
}


///////// VECTOR
NLA.Vector = function (v) {
	if (NLA.DEBUG && !(v instanceof Float64Array)) {
		throw new Error("!!");
	}
	this.v = v;
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
NLA.Vector.fromArguments = function () {
	assert (arguments[0] instanceof Float64Array || arguments.every(a => 'number' == typeof a),
	"arguments[0] instanceof Float64Array || arguments.every(a => 'number' == typeof a)")
	return new NLA.Vector(arguments[0] instanceof Float64Array ? arguments[0] : Float64Array.from(arguments))
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
	minus: function (vector) {
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
		NLA.assert(this.dim == vector.dim, "passed vector must have the same dim")
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
		assert (vector instanceof NLA.Vector, "vector instanceof NLA.Vector")
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
	angleTo: function (vector) {
		assert (vector instanceof NLA.Vector, "vector instanceof NLA.Vector")
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
		assert (vector instanceof NLA.Vector, "vector instanceof NLA.Vector")
		assert(!this.isZero(), "!this.isZero()")
		assert(!vector.isZero(), "!vector.isZero()")
		// a . b takes on values of +|a|*|b| (vectors same direction) to -|a|*|b| (opposite direction)
		// in both cases the vectors are paralle, so check if abs(a . b) == |a|*|b|
		return NLA.equals(Math.sqrt(this.lengthSquared() * vector.lengthSquared()), Math.abs(this.dot(vector)))
	},
	isPerpendicularTo: function (vector) {
		assert (vector instanceof NLA.Vector, "vector instanceof NLA.Vector")
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
		Throws a NLA.DebugError if this has a length of 0.
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
		assert (b instanceof NLA.Vector, "b instanceof NLA.Vector")
		// https://en.wikipedia.org/wiki/NLA.Vector_projection#NLA.Vector_projection_2
		return b.times(this.dot(b) / b.dot(b))
	},
	rejectedOn: function (b) {
		assert (b instanceof NLA.Vector, "b instanceof NLA.Vector")
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
};
NLA.addOwnProperties = function (target, props) {
	for (var key in props) {
		if (props.hasOwnProperty(key)) {
			if (target.hasOwnProperty(key)) {
				console.warn("target ", target, " already has property ", key)
			}
			Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(props, key))
		}
	}
}
NLA.defineClass = function (name, parent, constructor, props, statics) {
	NLA.assert('function' == typeof constructor, "'function' == typeof constructor")
	constructor.prototype = NLA.defineObject(parent && parent.prototype, props)
	constructor.prototype.constructor = constructor
	constructor.prototype.name = name
	NLA.addOwnProperties(constructor, statics)
	return constructor
}
NLA.defaultRoundFunction = x => x.toFixed(2)

NLA.defineObject = function (prot, props) {
	var o = Object.create(prot || NLA.baseObject)
	NLA.addOwnProperties(o, props)
	return o
}
NLA.baseObject = NLA.defineObject(Object.prototype, {
	toSource: function () {
		return this.toString == Object.prototype.toString ? Object.prototype.toSource.apply(this) : this.toString()
	},
	get ss() {
		return this.toSource()
	}
})
NLA.forceFinite = function(val) {
	val = parseFloat(val.replace(',', '.').replace(/^[^0-9,\.\-]/, ''))
	return Number.isFinite(val) ? val : 0
}
// Closure
;(function() {
	/**
	 * Decimal adjustment of a number.
	 *
	 * @param {String}  type  The type of adjustment.
	 * @param {Number}  value The number.
	 * @param {Integer} exp   The exponent (the 10 logarithm of the adjustment base).
	 * @returns {Number} The adjusted value.
	 */
	function decimalAdjust(type, value, exp) {
		// If the exp is undefined or zero...
		if (typeof exp === 'undefined' || +exp === 0) {
			return Math[type](value);
		}
		value = +value;
		exp = +exp;
		// If the value is not a number or the exp is not an integer...
		if (isNaN(value) || !(typeof exp === 'number' && exp % 1 === 0)) {
			return NaN;
		}
		// Shift
		value = value.toString().split('e');
		value = Math[type](+(value[0] + 'e' + (value[1] ? (+value[1] - exp) : -exp)));
		// Shift back
		value = value.toString().split('e');
		return +(value[0] + 'e' + (value[1] ? (+value[1] + exp) : exp));
	}

	// Decimal round
	if (!Math.round10) {
		Math.round10 = function(value, exp) {
			return decimalAdjust('round', value, exp);
		};
	}
	// Decimal floor
	if (!Math.floor10) {
		Math.floor10 = function(value, exp) {
			return decimalAdjust('floor', value, exp);
		};
	}
	// Decimal ceil
	if (!Math.ceil10) {
		Math.ceil10 = function(value, exp) {
			return decimalAdjust('ceil', value, exp);
		};
	}
})();






/// MATRIX
NLA.Matrix = function (width, height, m) {
	NLA.assert(width * height == m.length, "uhm...")
	this.m = m;
	this.width = width;
	this.height = height;
}
NLA.Matrix.random = function (width, height) {
	NLA.assertNumbers(width, height);
	return NLA.Matrix.fromFunction(width, height, (i, j) => Math.random())
}
NLA.Matrix.fromFunction = function(width, height, f) {
	NLA.assertNumbers(width, height);
	var m = new Float64Array(height * width)
	var elIndex = height * width;
	while (elIndex--) {
		m[elIndex] = f(Math.floor(elIndex / width), elIndex % width, elIndex);
	}
	return new NLA.Matrix(width, height, m);
}
NLA.Matrix.identity = function(dim) {
	NLA.assertNumbers(dim);
	var m = new Float64Array(dim * dim)
	// Float64Array are init to 0
	var elIndex = dim * (dim + 1);
	while (elIndex) {
		elIndex -= (dim + 1)
		m[elIndex] = 1
	}
	return new NLA.Matrix(dim, dim, m);
}
NLA.Matrix.permutation = function(dim, i, k) {
	NLA.assertNumbers(dim, i, k);
	var m = new Float64Array(dim * dim)
	// Float64Array are init to 0
	var elIndex = dim * (dim + 1);
	while (elIndex) {
		elIndex -= (dim + 1)
		m[elIndex] = 1
	}
	m[i * dim + i] = 0
	m[k * dim + k] = 0
	m[i * dim + k] = 1
	m[k * dim + i] = 1
	return new NLA.Matrix(dim, dim, m);
}
NLA.Matrix.fromRowArrays = function() {
	return NLA.Matrix.fromRowArrays2(arguments)
}
NLA.Matrix.fromRowArrays2 = function(arrays) {
	if (0 == arrays.length) {
		throw new Error("cannot have 0 vector");
	}
	var height = arrays.length
	var width = arrays[0].length
	var m = new Float64Array(height * width)
	NLA.arrayCopy(arrays[0], 0, m, 0, width);
	for (var rowIndex = 1; rowIndex < height; rowIndex++) {
		if (arrays[rowIndex].length != width) {
			throw new Error("all row arrays must be the same length")
		}
		NLA.arrayCopy(arrays[rowIndex], 0, m, rowIndex * width, width)
	}
	return new NLA.Matrix(width, height, m);
}
NLA.Matrix.fromColVectors = function (colVectors) {
	return NLA.Matrix.fromColArrays(colVectors.map((v) => v.v))
}
NLA.Matrix.forWidthHeight = function(width, height) {
	return new NLA.Matrix(width, height, new Float64Array(width * height))
}
NLA.Matrix.fromColArrays = function(colArrays) {
	if (0 == colArrays.length) {
		throw new Error("cannot have 0 vector");
	}
	var width = colArrays.length
	var height = colArrays[0].length
	var m = new Float64Array(height * width)
	NLA.arrayCopyStep(colArrays[0], 0, 1, m, 0, width, height)
	for (var colIndex = 1; colIndex < width; colIndex++) {
		if (colArrays[colIndex].length != height) {
			throw new Error("all col arrays must be the same length")
		}
		NLA.arrayCopyStep(colArrays[colIndex], 0, 1, m, colIndex, width, height)
	}
	return new NLA.Matrix(width, height, m);
}
NLA.Matrix.prototype = {
    e: function (rowIndex, colIndex) {
        NLA.assertNumbers(rowIndex, colIndex)
        if (NLA.DEBUG && (rowIndex >= this.height || colIndex >= this.width)) {
            throw new Error("index " + rowIndex + ", " + colIndex + " is out of bounds (" + this.width + " x " + this.height + ")")
        }
        return this.m[rowIndex * this.width + colIndex]
    },
	setEl: function (rowIndex, colIndex, val) {
		NLA.assertNumbers(rowIndex, colIndex)
		NLA.assert(0 <= rowIndex && rowIndex < this.height, "rowIndex out of bounds " + rowIndex)
		NLA.assert(0 <= colIndex && colIndex < this.width, "colIndex out of bounds " + colIndex)
		this.m[rowIndex * this.width + colIndex] = val
	},
	toString: function (f) {
		f = f || ((v) => v.toFixed(3))
		NLA.assert(typeof f(0) == "string", "" + typeof f(0))
		var rounded = Array.prototype.slice.call(this.m).map(f);
		var colWidths = NLA.arrayFromFunction(this.width,
			(colIndex) => rounded.sliceStep(colIndex, this.width).map((x) => x.length).max())
		return NLA.arrayFromFunction(this.height,
			(rowIndex) => rounded.slice(rowIndex * this.width, (rowIndex + 1) * this.width) // select matrix row
				.map((x, colIndex) => NLA.repeatString(colWidths[colIndex] - x.length, ' ') + x) // pad numbers with spaces to col width
				.join(" ")
		).map(x => x + "\n").join(""); // join rows
    },

    row: function (rowIndex) {
        var v = new Float64Array(this.width)
        NLA.arrayCopy(this.m, rowIndex * this.width, v, 0, this.width)
        return new NLA.Vector(v)
    },
    col: function (colIndex) {
        var v = new Float64Array(this.height)
        NLA.arrayCopyStep(this.m, colIndex, this.width, v, 0, 1, this.height)
        return new NLA.Vector(v)
    },
    dim: function() {
        return {width: this.width, height: this.height};
    },
	dimString: function () {
		return this.width + "x" + this.height
	},
    equals: function (obj) {
        if (obj.constructor != NLA.Matrix) return false;
        if (this.width != obj.width || this.height != obj.height) return false;
        var elIndex = this.m.length;
        while (elIndex--) {
            if (this.m[elIndex] != obj.m[elIndex]) return false;
        }
    },
    equalsMatrix: function (matrix, precision) {
	    precision = precision || NLA.PRECISION
        if (!(matrix instanceof NLA.Matrix)) throw new Error("not a matrix");
        if (this.width != matrix.width || this.height != matrix.height) return false;
        var elIndex = this.m.length;
        while (elIndex--) {
            if (Math.abs(this.m[elIndex] - matrix.m[elIndex]) >= precision) return false;
        }
        return true
    },
	hashCode: function () {
		var result = 0
		var elIndex = this.m.length;
		while (elIndex--) {
			result = result * 31 + NLA.floatHashCode(this.m[elIndex])
		}
		return result
	},
    isOrthogonal: function () {
        return this.isSquare() && this.transposed().times(this).equalsMatrix(NLA.Matrix.identity(this.width))
    },
    luDecomposition: function() {
        NLA.assert(this.isSquare())
        var resultRowArrays = this.asRowArrays(Float64Array)
        var dim = this.width
        var lRowArrays = NLA.arrayFromFunction(dim, (row) => new Float64Array(dim))
        var pRowArrays = NLA.Matrix.identity(dim).asRowArrays(Float64Array)
        for (var colIndex = 0; colIndex < dim; colIndex++) {
            // find largest value in colIndex
            var maxAbsValue = 0
            var rowIndex = dim
            var pivotRowIndex = colIndex
            while (--rowIndex >= colIndex) {
                if (Math.abs(resultRowArrays[rowIndex][colIndex]) >= maxAbsValue) { // TODO: ">" ?
                    maxAbsValue = Math.abs(resultRowArrays[rowIndex][colIndex])
                    pivotRowIndex = rowIndex
                }
            }
            // swap row
			NLA.arraySwap(resultRowArrays, colIndex, pivotRowIndex)
			NLA.arraySwap(lRowArrays, colIndex, pivotRowIndex)
            NLA.arraySwap(pRowArrays, colIndex, pivotRowIndex)
            // subtract pivot row from all below it
            lRowArrays[colIndex][colIndex] = 1
            for (var rowIndex = colIndex + 1; rowIndex < dim; rowIndex++) {
                var l = resultRowArrays[rowIndex][colIndex] / resultRowArrays[colIndex][colIndex]
                lRowArrays[rowIndex][colIndex] = l
                // subtract pivot row * l from row "rowIndex"
                for (var colIndex2 = colIndex; colIndex2 < dim; colIndex2++) {
                    resultRowArrays[rowIndex][colIndex2] -= l * resultRowArrays[colIndex][colIndex2]
                }
            }
        }
        return {L: NLA.Matrix.fromRowArrays2(lRowArrays), U: NLA.Matrix.fromRowArrays2(resultRowArrays), P: NLA.Matrix.fromRowArrays2(pRowArrays)}
    },
	qrDecompositionGivensRotation: function () {
		function sigma (c, s) {
			if ( 0 == c) {
				return 1;
			}
			if ( Math.abs(s) < Math.abs(c)) {
				return 0.5 * Math.sign(c) * s;
			}
			return 2 * Math.sign (s) / c;
		}
		function matrixForCS(dim, i, k, c, s) {
			var m = NLA.Matrix.identity(dim)
			m.setEl(i, i, c)
			m.setEl(k, k, c)
			m.setEl(i, k, s)
			m.setEl(k, i, -s)
			return m
		}

		var qTransposed = NLA.Matrix.identity(this.height)
		for (var colIndex = 0; colIndex < this.width; colIndex++) {
			// find largest value in colIndex
			var maxAbsValue = 0
			var pivotRowIndex = colIndex
			for (var rowIndex = colIndex + 1; rowIndex < this.height; rowIndex++) {
				//console.log("row ", rowIndex, "col ", colIndex);
				var xi = this.e(colIndex, colIndex)
				var xk = this.e(rowIndex, colIndex)
				if (xk == 0) {
					continue;
				}
				var r = Math.sqrt(xi * xi + xk * xk);
				var c = xi / r
				var s = xk / r

				// apply transformation on every column:
				for (var col2 = colIndex; col2 < this.width; col2++) {
					var x1 = this.e(colIndex, col2) * c + this.e(rowIndex, col2) * s
					var x2 = this.e(rowIndex, col2) * c - this.e(colIndex, col2) * s
					this.setEl(colIndex, col2, x1)
					this.setEl(rowIndex, col2, x2)
				}
				//console.log("r ", r, "c ", c, "s ", s, "sigma", sigma(c, s));
				//console.log(this.toString(),"cs\n", matrixForCS(this.height, colIndex, rowIndex, c, s).toString())
				qTransposed = matrixForCS(this.height, colIndex, rowIndex, c, s).times(qTransposed)
			}
		}
		//console.log(qTransposed.transposed().toString(), this.toString(), qTransposed.transposed().times(this).toString())
		return {Q: qTransposed.transposed(), R: this}
	},
    isPermutation: function () {
        if (!this.isSquare()) return false
        if (this.m.some((value) => !NLA.isZero(value) && !NLA.equals(1, value))) return false
        var rows = this.asRowArrays(Array)
        if (rows.some((row) => row.filter((value) => NLA.equals(1, value)).length != 1)) return false
        var cols = this.asColArrays(Array)
        if (cols.some((col) => col.filter((value) => NLA.equals(1, value)).length != 1)) return false
        return true
    },
	isIdentity: function (precision) {
		return this.isLowerUnitriangular(precision) && this.isUpperTriangular(precision)
	},
    isUpperTriangular: function () {
        if (!this.isSquare()) return false
        for (var rowIndex = 1; rowIndex < this.height; rowIndex++) {
            for (var colIndex = 0; colIndex < rowIndex; colIndex++) {
                if (!NLA.isZero(this.m[rowIndex * this.width + colIndex])) {
                    return false
                }
            }
        }
        return true
    },
    /**
     * Returns x, so that this * x = b
	 * More efficient than calculating the inverse for few (~ <= this.height) values
     * @param y
     */
    solveLinearSystem: function (b) {
        var lup = this.luDecomposition()
        var y = lup.L.solveForwards(lup.P.timesVector(b))
        var x = lup.U.solveBackwards(y)
        return x
    },
	isLowerUnitriangular: function (precision) {
		precision = "number" == typeof precision ? precision : NLA.PRECISION
		if (!this.isSquare()) return false
		for (var rowIndex = 0; rowIndex < this.height - 1; rowIndex++) {
			for (var colIndex = rowIndex; colIndex < this.width; colIndex++) {
				var el = this.m[rowIndex * this.width + colIndex];
                if (rowIndex == colIndex ? !NLA.equals2(1, el, precision) : !NLA.isZero2(el, precision)) {
					return false
				}
			}
		}
		return true
	},
    isLowerTriangular: function () {
        if (!this.isSquare()) return false
        for (var rowIndex = 0; rowIndex < this.height - 1; rowIndex++) {
            for (var colIndex = rowIndex + 1; colIndex < this.width; colIndex++) {
                if (!NLA.isZero(this.m[rowIndex * this.width + colIndex])) {
                    return false
                }
            }
        }
        return true
    },
    solveBackwards: function (x) {
        NLA.assertVectors(x)
        NLA.assert(this.height == x.dim(), "this.height == x.dim()")
        NLA.assert(this.isUpperTriangular(), "this.isUpperTriangular()")
        var v = new Float64Array(this.width)
        var rowIndex = this.height
        while (rowIndex--) {
            var temp = x.v[rowIndex]
            for (var colIndex = rowIndex + 1; colIndex < this.width; colIndex++) {
                temp -= v[colIndex] * this.e(rowIndex, colIndex)
            }
            v[rowIndex] = temp / this.e(rowIndex, rowIndex)
        }
        return new NLA.Vector(v)
    },
    solveBackwardsMatrix: function (matrix) {
        var colVectors = new Array(matrix.width)
        var i = matrix.width
        while (i--) {
            colVectors[i] = this.solveBackwards(matrix.col(i))
        }
        return NLA.Matrix.fromColVectors(colVectors)
    },
    solveForwardsMatrix: function (matrix) {
        var colVectors = new Array(matrix.width)
        var i = matrix.width
        while (i--) {
            colVectors[i] = this.solveForwards(matrix.col(i))
        }
        return NLA.Matrix.fromColVectors(colVectors)
    },
    solveForwards: function (x) {
        NLA.assertVectors(x)
        NLA.assert(this.height == x.dim(), "°°")
        NLA.assert(this.isLowerTriangular(), "°°")
        var v = new Float64Array(this.width)
        for (var rowIndex = 0; rowIndex < this.height; rowIndex++) {
            var temp = x.v[rowIndex]
            for (var colIndex = 0; colIndex < rowIndex; colIndex++) {
                temp -= v[colIndex] * this.e(rowIndex, colIndex)
            }
            v[rowIndex] = temp / this.e(rowIndex, rowIndex)
        }
        return new NLA.Vector(v)
    },
	rank: function () {
		var {Q, R} = this.qrDecompositionGivensRotation()
		//console.log(R.toString())
		var rowIndex = this.height
		while (rowIndex-- && R.row(rowIndex).isZero()) {}
		return rowIndex + 1
	},
	rowsIndependent: function () {
		return this.height == this.rank()
	},
	colsIndependent: function () {
		var cols = NLA.arrayFromFunction(this.width, this.col.bind(this))
		return Array.from(NLA.combinations(this.width)).every(({i, j}) => i == j || !cols[i].isParallelTo(rows[j]))
	},
    asRowArrays: function (arrayConstructor) {
        arrayConstructor = arrayConstructor || Float64Array
        var rowIndex = this.height
        var result = new Array(this.height)
        while (rowIndex--) {
            result[rowIndex] = this.rowArray(rowIndex, arrayConstructor)
        }
        return result
    },
    asColArrays: function (arrayConstructor) {
        arrayConstructor = arrayConstructor || Float64Array
        var colIndex = this.width
        var result = new Array(this.width)
        while (colIndex--) {
            result[colIndex] = this.colArray(colIndex, arrayConstructor)
        }
        return result
    },
    rowArray: function (rowIndex, arrayConstructor) {
        arrayConstructor = arrayConstructor || Float64Array
        var result = new arrayConstructor(this.width)
        NLA.arrayCopy(this.m, rowIndex * this.width, result, 0, this.width)
        return result
    },
    colArray: function (colIndex, arrayConstructor) {
        arrayConstructor = arrayConstructor || Float64Array
        var result = new arrayConstructor(this.width)
        NLA.arrayCopyStep(this.m, colIndex, this.height, result, 0, 1, this.height)
        return result
    },
    subMatrix: function (firstColIndex, subWidth, firstRowIndex, subHeight) {
        if (firstColIndex + subWidth > this.width || firstRowIndex + subHeight > this.height) {
            throw new Error("inavlid params")
        }
        var m = new Float64Array(this.height)
        NLA.arrayCopyBlocks(this.m, firstColIndex, this.width, m, 0, subWidth, subHeight, subWidth)
        return new NLA.Matrix(subWidth, subHeight, m)
    },
    map: function (fn) {
        return new NLA.Matrix(this.width, this.height, this.m.map(fn));
    },

    dimEquals: function (matrix) {
        NLA.assert(matrix instanceof NLA.Matrix, "matrix is not a NLA.Matrix")
        return this.width == matrix.width && this.height == height;
    },
    inversed: function () {
        var lup = this.luDecomposition()
        var y = lup.L.solveForwardsMatrix(lup.P)
        var inverse = lup.U.solveBackwardsMatrix(y)
        return inverse
    },
    canMultiply: function (matrix) {
        NLA.assert(matrix instanceof NLA.Matrix, "matrix is not a NLA.Matrix" + matrix)
        return this.width == matrix.height
    },
    times: function (matrix) {
	    NLA.assert(matrix instanceof NLA.Matrix)
        NLA.assert(this.canMultiply(matrix), `Cannot multiply this {this.dimString()} by matrix {matrix.dimString()}`)
        var nWidth = matrix.width, nHeight = this.height, n = this.width
        var nM = new Float64Array(nWidth * nHeight)
        var nRowIndex = nHeight
        while (nRowIndex--) {
            var nColIndex = nWidth
            while (nColIndex--) {
                var result = 0
                var i = n
                while (i--) {
                    result += this.m[nRowIndex * n + i] * matrix.m[i * nWidth + nColIndex]
                }
                nM[nRowIndex * nWidth + nColIndex] = result
            }
        }
        return new NLA.Matrix(nWidth, nHeight, nM)
    },
    timesVector: function (vector) {
        NLA.assertVectors(vector)
        NLA.assert(this.width == vector.dim())
        var nHeight = this.height, n = this.width
        var nM = new Float64Array(nHeight)
        var nRowIndex = nHeight
        while (nRowIndex--) {
            var result = 0
            var i = n
            while (i--) {
                result += this.m[nRowIndex * n + i] * vector.v[i]
            }
            nM[nRowIndex] = result
        }
        return new NLA.Vector(nM)
    },
	transposed: function () {
		var tWidth = this.height, tHeight = this.width
		var tM = new Float64Array(tWidth * tHeight)
		var tRowIndex = tHeight
		while (tRowIndex--) {
			var tColIndex = tWidth
			while (tColIndex--) {
				tM[tRowIndex * tWidth + tColIndex] = this.m[tColIndex * tHeight + tRowIndex]
			}
		}
		return new NLA.Matrix(tWidth, tHeight, tM)
	},
	transpose: function () {
		var h = this.height, w = this.width, tM = this.m
		var tRowIndex = h
		while (tRowIndex--) {
			var tColIndex = Math.min(tRowIndex, w)
			while (tColIndex--) {
				console.log("col", tColIndex, "row", tRowIndex)
				var temp = tM[tRowIndex * w + tColIndex]
				tM[tRowIndex * w + tColIndex] = tM[tColIndex * h + tRowIndex]
				tM[tColIndex * h + tRowIndex] = temp
			}
		}
		this.width = h
		this.height = w
	},
    isSquare: function () {
        return this.height == this.width;
    },
    diagonal: function () {
        if (!this.isSquare()) {
            throw new Error("!!")
        }
        var v = new Float64Array(this.width)
        var elIndex = this.width * (this.width + 1)
        var vIndex = this.width
        while (vIndex--) {
            elIndex -= this.width + 1
            v[vIndex] = this.m[elIndex]
        }
        return new NLA.Vector(v)
    },
    max: function () {
        return NLA.arrayMax(this.m);
    },
    min: function () {
        return NLA.arrayMin(this.m);
    },
    maxAbsColSum: function () {
        var result = 0
        var colIndex = this.width
        while (colIndex--) {
            var absSum = 0
            var rowIndex = this.height
            while (rowIndex--) {
                absSum  += Math.abs(this.m(rowIndex * this.width + colIndex))
            }
            result = Math.max(result, absSum)
        }
        return result
    },
    maxAbsRowSum: function () {
        var result = 0
        var rowIndex = this.height
        while (rowIndex--) {
            var absSum = 0
            var colIndex = this.width
            while (colIndex--) {
                absSum += Math.abs(this.m(rowIndex * this.width + colIndex))
            }
            result = Math.max(result, absSum)
        }
        return result
    },
	getTriangularDeterminant: function () {
		assert(this.isUpperTriangular() || this.isLowerTriangular(), "not a triangular matrix")

		var product = 1
		var elIndex = this.width * (this.width + 1)
		while (elIndex) {
			elIndex -= this.width + 1
			product *= this.m[elIndex]
		}
		return product
	},
	/**
	 * Calculates the determinant by first calculating the LU decomposition. If you already have that, use U.getTriangularDeterminant()
	 * @returns {*}
	 */
	getDeterminant: function () {
		// PA = LU
		// det(A) * det(B) = det(AB)
		// det(P) == 1 (permutation matrix)
		// det(L) == 1 (main diagonal is 1s
		// =>  det(A) == det(U)
		return this.luDecomposition().U.getTriangularDeterminant()
	},
	hasFullRank: function () {
		return Math.min(this.width, this.height) == this.rank()
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

NLA.combinations = function* (n) {
	for (var i = 0; i < n; i++) {
		for (var j = i; j < n; j++) {
			yield {i: i, j: j}
		}
	}
}

NLA.Line = function (anchor, dir1) {
	NLA.assertVectors(anchor, dir1)
    console.log("sadjlkasjd", dir1)
    try {
        if (!dir1.hasLength(1)) {
            throw new Error("dir must be normalized");
        }
    } catch (e) {
        console.log(e)
    }
    this.anchor = anchor
	this.dir1 = dir1
}
NLA.Line.throughPoints = (anchor, b) => new NLA.Line(anchor, b.minus(anchor).normalized())
NLA.Line.anchorDirection = (anchor, direction) => new NLA.Line(anchor, direction.normalized())
NLA.Line.prototype = {
	containsPoint: function () {
	},
	equals: function (line) {
		NLA.assert(line instanceof NLA.Line);
        // we know that 1 == this.dir1.length() == line.dir1.length(), we can check for parallelity simpler than isParallelTo()
		return this.contains(line.anchor) && NLA.equals(1, Math.abs(this.dir1.dot(line.dir1)))
	},
	distanceToPoint: function (x) {
		NLA.assertVectors(x)
		"use strict";
		// See http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
		var t = x.minus(this.anchor).dot(this.dir1)
		return this.at(t).minus(x).length()

		//return x.minus(this.anchor).cross(x.minus(this.anchor.plus(this.dir1))).length()
	},
    /**
     * Every point x on this line is described by the equation x = this.anchor + lambda * this.dir1
     * This function returns lambda for a given point x
     * @param x
     */
    pointLambda: function (x) {
        "use strict";
        NLA.assertVectors(x)
        var t = x.minus(this.anchor).dot(this.dir1)
    },
    isParallelToLine: function (line) {
        "use strict";
        NLA.assert(line instanceof NLA.Line)
        // we know that 1 == this.dir1.length() == line.dir1.length(), we can check for parallelity simpler than isParallelTo()
        return NLA.equals(1, Math.abs(this.dir1.dot(line.dir1)))
    },
    angleToLine:  function (line) {
        "use strict";
        NLA.assert(line instanceof NLA.Line)
        return this.dir1.angleTo(line.dir1)
    },
	/**
	 *
	 * @param t
	 * @returns `this.anchor.plus(this.dir1.times(t))`
     */
	at: function (t) {
        NLA.assertNumbers(t)
		"use strict";
		return this.anchor.plus(this.dir1.times(t))
	},
	pointClosestTo2: function(obj) {
		if (obj.direction) {
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
			// obj is a line
			if (this.intersects(obj)) { return this.intersectionWith(obj); }
			if (this.isParallelTo(obj)) { return {t: NaN, s:NaN, closest: null, distance: this.distanceTo(obj)}; }
			var a = obj.anchor, b = obj.direction, c = this.anchor, d = this.direction;
			var bd = b.dot(d), bb = b.dot(b), dd = d.dot(d), amc = a.subtract(c), divisor = bd*bd - dd * bb;
			var t = (amc.dot(b)*bd - amc.dot(d)*bb) / divisor;
			var s = (amc.dot(b)*dd - amc.dot(d)*bd) / divisor;
			return {t: t, s:s, closest: this.at(t), closest2: obj.at(s), distance: this.at(t).subtract(obj.at(s)).modulus()};
		} else {
			// obj is a point
			var P = obj.elements || obj;
			if (this.contains(P)) { return Vector.create(P); }
			var A = this.anchor.elements, D = this.direction.elements;
			var D1 = D[0], D2 = D[1], D3 = D[2], A1 = A[0], A2 = A[1], A3 = A[2];
			var x = D1 * (P[1]-A2) - D2 * (P[0]-A1), y = D2 * ((P[2] || 0) - A3) - D3 * (P[1]-A2),
					z = D3 * (P[0]-A1) - D1 * ((P[2] || 0) - A3);
			var V = Vector.create([D2 * x - D3 * z, D3 * y - D1 * x, D1 * z - D2 * y]);
			var k = this.distanceTo(P) / V.modulus();
			return Vector.create([
				P[0] + V.elements[0] * k,
				P[1] + V.elements[1] * k,
				(P[2] || 0) + V.elements[2] * k
			]);
		}
	},
}
NLA.tranformablePrototype = (function() {
	var prot = {}
	prot.mirrored = function(plane) {
		return this.transform(NLA.Matrix4x4.mirroring(plane));
	}
	prot.mirroredX = function() {
		return this.mirrored(P3.YZ)
	}
	prot.mirroredY = function() {
		return this.mirrored(P3.ZX)
	}
	prot.mirroredZ = function() {
		return this.mirrored(P3.XY)
	}
	prot.translate = function(x, y, z) {
		return this.transform(NLA.Matrix4x4.translation(x, y, z))
	}
	prot.scale = function(f) {
		return this.transform(NLA.Matrix4x4.scaling(f))
	}
	prot.rotateX = function(radians) {
		return this.transform(NLA.Matrix4x4.rotationX(radians))
	}
	prot.rotateY = function(radians) {
		return this.transform(NLA.Matrix4x4.rotationY(radians))
	}
	prot.rotateZ = function(radians) {
		return this.transform(NLA.Matrix4x4.rotationZ(radians))
	}
	prot.rotate = function(rotationCenter, rotationAxis, radians) {
		return this.transform(NLA.Matrix4x4.rotation(rotationCenter, rotationAxis, radians))
	}
	prot.eulerZXZ = function(alpha, beta, gamma) {
		return this.transform(NLA.Matrix4x4.eulerZXZ(alpha, beta, gamma))
	}
	prot.project = function(plane) {
		return this.transform(NLA.Matrix4x4.projection(plane))
	}
	return prot
})();
NLA.addOwnProperties(Array.prototype, ARRAY_UTILITIES)