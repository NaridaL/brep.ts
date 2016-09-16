"use strict"
/**
 * @namespace
 */
var NLA = {}
window['NLA'] = NLA

function SCE(o) {
	return o.sce
}
function STR(o) {
	return o.str
}
Object.defineProperty(Object.prototype, 'sce', {get: function () { return this.toSource() }})
Object.defineProperty(Object.prototype, 'str', {get: function () { return this.toString() }})


var assertVectors = NLA.assertVectors = function () {
	if (NLA_DEBUG) {
		for (var i = 0; i < arguments.length; i++) {
			if (!(arguments[i] instanceof V3 || arguments[i] instanceof NLA.Vector)) {
				throw new Error("NLA.assertVectors arguments[" + (i) + "] is not a vector. " + typeof arguments[i] + " == typeof " + arguments[i]);
			}
		}
	}
	return true
}
var assertInst = NLA.assertInst = function (what, ...objs) {
	if (NLA_DEBUG) {
		for (var i = 0; i < objs.length; i++) {
			if (!(objs[i] instanceof what)) {
				throw new Error("assertInst objs[" + (i) + "] is not a "+what.prototype.constructor.name+". " + objs[i].constructor.name + objs[i])
			}
		}
	}
	return true
}
var assertNumbers = NLA.assertNumbers = function (...numbers) {
	if (NLA_DEBUG) {
		for (var i = 0; i < numbers.length; i++) {
			if (typeof numbers[i] !== 'number') {
				throw new Error("NLA.assertNumbers arguments[" + (i) + "] is not a number. " + typeof numbers[i] + " == typeof " + numbers[i]);
			}
		}
	}
	return true
}
/**
 *
 * @param value
 * @param {string} [message]
 * @returns {boolean}
 */
var assert = NLA.assert = function (value, message) {
	if (NLA_DEBUG && !value) {
		throw new Error("NLA.assert failed: " + ('function' == typeof message ? message() : message))
	}
	return true
}

var assertf = NLA.assertf = function (f, message) {
	if (!f()) {
		throw new Error("NLA.assertf failed: " + f.toString() + ('function' == typeof message ? message() : (message ? message : '')))
	}
}
NLA.CLASSES = []
NLA.registerClass = function (clazz) {
	NLA.CLASSES.push(clazz)
}

/** @const */
NLA.PRECISION = 1 / (1 << 26)

console.log("NLA.PRECISION", NLA.PRECISION)

/** @define {boolean} */
var NLA_DEBUG = true
/**
 *
 * @param {number} x
 * @returns {boolean}
 */
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
NLA.snapTo = (x, to) => Math.abs(x - to) <= NLA.PRECISION ? to : x
NLA.canonAngle = (x) => ((x % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI)

var COLORS = {
	RD_FILL:0x9EDBF9,
	RD_STROKE:0x77B0E0,
	TS_FILL: 0xD19FE3,
	TS_STROKE: 0xA76BC2,
	PP_FILL:0xF3B6CF,
	PP_STROKE:0xEB81B4,
}
var DEG = .017453292519943295
function deg2rad(angle) {
	//  discuss at: http://phpjs.org/functions/deg2rad/
	// original by: Enrique Gonzalez
	// improved by: Thomas Grainger (http://graingert.co.uk)
	//   example 1: deg2rad(45);
	//   returns 1: 0.7853981633974483

	return angle * .017453292519943295 // (angle / 180) * Math.PI;
}
function rad2deg(rad) {
	//  discuss at: http://phpjs.org/functions/deg2rad/
	// original by: Enrique Gonzalez
	// improved by: Thomas Grainger (http://graingert.co.uk)
	//   example 1: deg2rad(45);
	//   returns 1: 0.7853981633974483

	return rad / .017453292519943295 // (angle / 180) * Math.PI;
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
/**
 *
 * @param iterable
 * @constructor
 */
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
/**
 *
 * @constructor
 */
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

String.prototype.capitalizeFirstLetter = function() {
	return this.charAt(0).toUpperCase() + this.slice(1);
}
var ARRAY_UTILITIES = /** @template T
 @lends Array.<T>.prototype */ {
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

	/**
	 * Semantically identical to .map(f).filter(v => v)
	 * @template S
	 * @param {function (T, number, Array.<T>):S} f
	 * @returns {Array.<S>}
	 */
	mapFilter(f) {
		let length = this.length, result = []
		for (let i = 0; i < length; i++) {
			if (i in this) {
				let val = f(this[i], i, this)
				if (val) {
					result.push(val)
				}
			}
		}
		return result
	},

	/**
	 *
	 * @returns {Array} Array.prototype.concat.apply([], this)
	 */
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
		if (this.length == 0) { return 0 }
		var i = this.length, result = -1, maxVal = -Infinity
		while (i--) {
			var val = f(this[i], i)
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

	/**
	 * @template T
	 * @template S
	 * @this {T[]}
	 * @param {S} searchElement
	 * @param {function (T, S): number} cmp
	 * @returns {number}
	 */
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

		return -minIndex-1;
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
	},
	last: function () {
		return this[this.length - 1]
	},
	// forEachCall: function (f, ...args) {
	// 	for (let i = 0; i < this.length; i++) {
	// 		f.apply(this[i], args)
	// 	}
	// },
	// mapCall: function (f, ...args) {
	// 	return this.map(e => f.apply(e, args))
	// }
}
NLA.minus = (a, b) => a - b
for (let key in ARRAY_UTILITIES) {
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
	//console.log(Math.ceil((end - start) / step))
	var result = new Array(Math.ceil((end - start) / step)); // "- start" so that chunk in the last row will also be selected, even if the row is not complete
	var index = 0
	for (var i = 0; i < end; i += step) {
		result[index++] = i
	}
	return result;
}

/**
 * @template T
 * @param {number} length
 * @param {function (number): T} f
 * @returns {Array.<T>}
 */
NLA.arrayFromFunction = function(length, f) {
	NLA.assertNumbers(length)
	assert("function" == typeof f)
	var a = new Array(length)
	var elIndex = length;
	while (elIndex--) {
		a[elIndex] = f(elIndex);
	}
	return a
}

/**
 * @param {number[]} vals
 * @returns {number[]}
 */
NLA.fuzzyUniques = function (vals) {
	let round = val => Math.floor(val * (1 << 26)) / (1 << 26)
	let map = new Map()
	map.set(round(vals[0]), vals[0])
	for (let i = 1; i < vals.length; i++) {
		let val = vals[i], roundVal = round(val)
		let key
		if (!map.has(roundVal)
			&& !((key = map.get(roundVal - 1 / (1 << 26))) && NLA.equals(key, val))
			&& !((key = map.get(roundVal + 1 / (1 << 26))) && NLA.equals(key, val))) {
			map.set(roundVal, val)
		}
	}
	return Array.from(map.values())
}

/**
 * @template T
 * @param {T[]} vals
 * @param {function(T):number} f
 * @returns {T[]}
 */
NLA.fuzzyUniquesF = function (vals, f) {
	let round = val => Math.floor(val * (1 << 26)) / (1 << 26)
	let map = new Map()
	for (let i = 0; i < vals.length; i++) {
		let val = vals[i], roundVal = round(f(val))
		let key
		if (!map.has(roundVal)
			&& !((key = map.get(roundVal - 1 / (1 << 26))) && NLA.equals(key, f(val)))
			&& !((key = map.get(roundVal + 1 / (1 << 26))) && NLA.equals(key, f(val)))) {
			map.set(roundVal, val)
		}
	}
	return Array.from(map.values())
}


/**
 *
 * @param v
 * @constructor
 */
NLA.Vector = function (v) {
	if (NLA_DEBUG && !(v instanceof Float64Array)) {
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
};
NLA.addOwnProperties = function (target, props) {
	Object.getOwnPropertyNames(props).forEach(key => {
		//console.log(props, key)
		if (target.hasOwnProperty(key)) {
			console.warn("target ", target, " already has property ", key)
		}
		Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(props, key))
	})
}

NLA.defineClass = function (name, parent, constructor, props, statics) {
	assert('function' == typeof constructor, "'function' == typeof constructor")
	constructor.prototype = NLA.defineObject(parent && parent.prototype, props)
	constructor.prototype.constructor = constructor
	Object.defineProperty(constructor.prototype, 'name', {value: name})
	statics && NLA.addOwnProperties(constructor, statics)
	return constructor
}
NLA.defaultRoundFunction = x => x // Math.round10(x, -4)

NLA.defineObject = function (prot, props) {
	var o = Object.create(prot || NLA.BaseObject)
	NLA.addOwnProperties(o, props)
	return o
}
NLA.BaseObject = class extends Object {
	toSource() {
		return this.toString == Object.prototype.toString ? Object.prototype.toSource.apply(this) : this.toString()
	}
	get $() {
		return this.toSource()
	}
}
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
	 * @param {number}  value The number.
	 * @param {number} exp   The exponent (the 10 logarithm of the adjustment base).
	 * @returns {number} The adjusted value.
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


/**
 *
 * @param width
 * @param height
 * @param m
 * @constructor
 * @property {Float64Array} m
 */
NLA.Matrix = function (width, height, m) {
	assert(width * height == m.length, "width * height == m.length", width, height, m.length)
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
	// Float64Arrays are init to 0
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
        if (NLA_DEBUG && (rowIndex >= this.height || colIndex >= this.width)) {
            throw new Error("index " + rowIndex + ", " + colIndex + " is out of bounds (" + this.width + " x " + this.height + ")")
        }
        return this.m[rowIndex * this.width + colIndex]
    },
	setEl: function (rowIndex, colIndex, val) {
		NLA.assertNumbers(rowIndex, colIndex, val)
		assert(0 <= rowIndex && rowIndex < this.height, "rowIndex out of bounds " + rowIndex)
		assert(0 <= colIndex && colIndex < this.width, "colIndex out of bounds " + colIndex)
		this.m[rowIndex * this.width + colIndex] = val
	},
	toString: function (f) {
		f = f || ((v) => v.toFixed(3))
		assert(typeof f(0) == "string", "" + typeof f(0))
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
        assertf(() => this.isSquare(), this.dim().toSource())
        var uRowArrays = this.asRowArrays(Float64Array)
        var dim = this.width
        var lRowArrays = NLA.arrayFromFunction(dim, (row) => new Float64Array(dim))
        var pRowArrays = NLA.Matrix.identity(dim).asRowArrays(Float64Array)
	    let currentRowIndex = 0
        for (let colIndex = 0; colIndex < dim; colIndex++) {
            // find largest value in colIndex
            let maxAbsValue = 0, pivotRowIndex = undefined, numberOfNonZeroRows = 0
            for (let rowIndex = currentRowIndex; rowIndex < dim; rowIndex++) {
	            var el = uRowArrays[rowIndex][colIndex]
	            numberOfNonZeroRows += (0 != el) | 0
	            if (Math.abs(el) > maxAbsValue) {
                    maxAbsValue = Math.abs(el)
                    pivotRowIndex = rowIndex
                }
            }
            // TODO: check with NLA.isZero
            if (0 == maxAbsValue) {
            	// column contains only zeros
	            continue
            }
	        assert(undefined !== pivotRowIndex)
            // swap rows
			NLA.arraySwap(uRowArrays, currentRowIndex, pivotRowIndex)
			NLA.arraySwap(lRowArrays, currentRowIndex, pivotRowIndex)
            NLA.arraySwap(pRowArrays, currentRowIndex, pivotRowIndex)
	        lRowArrays[currentRowIndex][colIndex] = 1

	        if (1 < numberOfNonZeroRows) {
		        // subtract pivot (now current) row from all below it
		        for (let rowIndex = currentRowIndex + 1; rowIndex < dim; rowIndex++) {
			        let l = uRowArrays[rowIndex][colIndex] / uRowArrays[currentRowIndex][colIndex]
			        lRowArrays[rowIndex][colIndex] = l
			        // subtract pivot row * l from row "rowIndex"
			        for (let colIndex2 = colIndex; colIndex2 < dim; colIndex2++) {
				        uRowArrays[rowIndex][colIndex2] -= l * uRowArrays[currentRowIndex][colIndex2]
			        }
		        }
	        }
            currentRowIndex++ // this doesn't increase if pivot was zero
        }
        return {L: NLA.Matrix.fromRowArrays2(lRowArrays), U: NLA.Matrix.fromRowArrays2(uRowArrays), P: NLA.Matrix.fromRowArrays2(pRowArrays)}
    },
    gauss: function() {
        var uRowArrays = this.asRowArrays(Float64Array)
        var w = this.width, h = this.height
        var lRowArrays = NLA.arrayFromFunction(h, (row) => new Float64Array(w))
        var pRowArrays = NLA.Matrix.identity(h).asRowArrays(Float64Array)
	    let currentRowIndex = 0
        for (let colIndex = 0; colIndex < w; colIndex++) {
        	// console.log('currentRowIndex', currentRowIndex)
            // find largest value in colIndex
            let maxAbsValue = 0, pivotRowIndex = undefined, numberOfNonZeroRows = 0
            for (let rowIndex = currentRowIndex; rowIndex < h; rowIndex++) {
	            var el = uRowArrays[rowIndex][colIndex]
	            numberOfNonZeroRows += (0 != el) | 0
	            if (Math.abs(el) > maxAbsValue) {
                    maxAbsValue = Math.abs(el)
                    pivotRowIndex = rowIndex
                }
            }
            // TODO: check with NLA.isZero
            if (0 == maxAbsValue) {
            	// column contains only zeros
	            continue
            }
	        assert(undefined !== pivotRowIndex)
            // swap rows
			NLA.arraySwap(uRowArrays, currentRowIndex, pivotRowIndex)
			NLA.arraySwap(lRowArrays, currentRowIndex, pivotRowIndex)
            NLA.arraySwap(pRowArrays, currentRowIndex, pivotRowIndex)
	        lRowArrays[currentRowIndex][colIndex] = 1

	        if (1 < numberOfNonZeroRows) {
		        // subtract pivot (now current) row from all below it
		        for (let rowIndex = currentRowIndex + 1; rowIndex < h; rowIndex++) {
			        let l = uRowArrays[rowIndex][colIndex] / uRowArrays[currentRowIndex][colIndex]
			        lRowArrays[rowIndex][colIndex] = l
			        // subtract pivot row * l from row "rowIndex"
			        for (let colIndex2 = colIndex; colIndex2 < w; colIndex2++) {
				        uRowArrays[rowIndex][colIndex2] -= l * uRowArrays[currentRowIndex][colIndex2]
			        }
		        }
	        }
	        currentRowIndex++ // this doesn't increase if pivot was zero
        }
        return {L: NLA.Matrix.fromRowArrays2(lRowArrays), U: NLA.Matrix.fromRowArrays2(uRowArrays), P: NLA.Matrix.fromRowArrays2(pRowArrays)}
    },
	qrDecompositionGivensRotation: function () {
		function sigma (c, s) {
			if ( 0 == c) {
				return 1
			}
			if ( Math.abs(s) < Math.abs(c)) {
				return 0.5 * Math.sign(c) * s
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
     * @param b
     */
    solveLinearSystem: function (b) {
        var lup = this.luDecomposition()
        var y = lup.L.solveForwards(lup.P.timesVector(b))
        var x = lup.U.solveBackwards(y)
        return x
    },

	/**
	 *
	 * @param {number=} precision
	 * @returns {boolean}
	 */
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

	/**
	 *
	 * @returns {boolean}
	 */
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


	/**
	 *
	 * @param x
	 * @returns {NLA.Vector}
	 */
    solveBackwards: function (x) {
        NLA.assertVectors(x)
        assert(this.height == x.dim(), "this.height == x.dim()")
        assert(this.isUpperTriangular(), "this.isUpperTriangular()")
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
        assert(this.height == x.dim(), "°°")
        assert(this.isLowerTriangular(), "°°")
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



	/**
	 * Calculates rank of matrix.
	 * Number of linearly independant row/column vectors.
	 * Is equal to the unmber of dimensions the image of the affine transformation represented this matrix has.
	 *
	 * @returns {number} integer
	 */
	rank: function () {
		let U = this.gauss().U
		//console.log(R.toString())
		var rowIndex = this.height
		while (rowIndex-- && U.row(rowIndex).isZero()) {
			console.log("RANK" + U.row(rowIndex).toString() + U.row(rowIndex).isZero())}
		return rowIndex + 1
	},
	rowsIndependent: function () {
		return this.height == this.rank()
	},
	colsIndependent: function () {
		return this.width == this.rank()
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
        assertInst(NLA.Matrix, matrix)
        return this.width == matrix.width && this.height == matrix.height
    },

    inversed: function () {
        var lup = this.luDecomposition()
        var y = lup.L.solveForwardsMatrix(lup.P)
        var inverse = lup.U.solveBackwardsMatrix(y)
        return inverse
    },

	inversed3: function () {
		assertf(() => 3 == this.width && 3 == this.height)
		let result = new NLA.Matrix.forWidthHeight(3, 3), m = this.m, r = result.m

		r[0] = m[4]*m[8] - m[5]*m[7]
		r[1] = -m[1]*m[8] + m[2]*m[7]
		r[2] = m[1]*m[5] - m[2]*m[4]

		r[3] = -m[3]*m[8] + m[5]*m[6]
		r[4] = m[0]*m[8] - m[2]*m[6]
		r[5] = -m[0]*m[5] + m[2]*m[3]

		r[6] = m[3]*m[7] - m[4]*m[6]
		r[7] = -m[0]*m[7] + m[1]*m[6]
		r[8] = m[0]*m[4] - m[1]*m[3]

		let det = m[0]*r[0] + m[1]*r[3] + m[2]*r[6]
		var i = 9
		while (i--) { r[i] /= det }

		return result
	},

	inversed2: function () {
		assertf(() => 2 == this.width && 2 == this.height)
		let result = new NLA.Matrix.forWidthHeight(2, 2), m = this.m, r = result.m

		let det = m[0]*m[3] - m[1]*r[2]

		r[0] = m[3] / det
		r[1] = -m[2] / det

		r[2] = -m[1] / det
		r[3] = m[0] / det

		return result
	},

    canMultiply: function (matrix) {
        assertInst(NLA.Matrix, matrix)
        return this.width == matrix.height
    },
    times: function (matrix) {
	    assertInst(NLA.Matrix, matrix)
        assert(this.canMultiply(matrix), `Cannot multiply this {this.dimString()} by matrix {matrix.dimString()}`)
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
        assert(this.width == vector.dim())
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
	 * Calculates the determinant by first calculating the LU decomposition. If you already have that, use
	 * U.getTriangularDeterminant()
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

/**
 * combinations(2) will generate
 * [0,0] [0,1] [1,1] [0,2] [1,2] [2,2]
 * @param n
 */
NLA.combinations = function* (n) {
	for (var i = 0; i < n; i++) {
		for (var j = i; j < n; j++) {
			yield {i: i, j: j}
		}
	}
}
function isCCW(vertices, normal) {
	var dsa = doubleSignedArea(vertices, normal)
	assert(0 != dsa)
	return dsa < 0
}
function triangulateVertices(normal, vertices, holeStarts) {
	var absMaxDim = normal.absMaxDim(), factor = normal.e(absMaxDim) < 0 ? -1 : 1
	var contour = new Array(vertices.length * 2)
	var i = vertices.length
	/*
	 var [coord0, coord1] = [['y', 'z'], ['z', 'x'], ['x', 'y']][absMaxDim]
	 while (i--) {
	 contour[i * 2    ] = vertices[i][coord0] * factor
	 contour[i * 2 + 1] = vertices[i][coord1]
	 }
	 */

	while (i--) {
		// unroll disambiguation instead of accessing elements by string name ([coord0] etc)
		// as it confuses google closure
		switch (absMaxDim) {
			case 0:
				contour[i * 2    ] = vertices[i].y * factor
				contour[i * 2 + 1] = vertices[i].z
				break
			case 1:
				contour[i * 2    ] = vertices[i].z * factor
				contour[i * 2 + 1] = vertices[i].x
				break
			case 2:
				contour[i * 2    ] = vertices[i].x * factor
				contour[i * 2 + 1] = vertices[i].y
				break
		}
	}
	return earcut(contour, holeStarts)
}

function arrayFilterMap(arr, f) {
	var result = []
	for (var i = 0; i < arr.length; i++) {
		var v = f(arr[i])
		if (undefined !== v) {
			result.push(v)
		}
	}
	return result
}


function doubleSignedArea(vertices, normal) {
	assert(!normal.isZero(),'!normal.isZero()')
	var absMaxDim = normal.absMaxDim()
	// order is important, coord0 and coord1 must be set so that coord0, coord1 and maxDim span a right-hand coordinate system
	//var [coord0, coord1] = [['y', 'z'], ['z', 'x'], ['x', 'y']][absMaxDim]
	var doubleSignedArea = vertices.map((v0, i, vertices) => {
		var v1 = vertices[(i + 1) % vertices.length]
		//return (v1[coord0] - v0[coord0]) * (v1[coord1] + v0[coord1])
		switch (absMaxDim) {
			case 0:
				return (v1.y - v0.y) * (v1.z + v0.z)
			case 1:
				return (v1.z - v0.z) * (v1.x + v0.x)
			case 2:
				return (v1.x - v0.x) * (v1.y + v0.y)
		}
	}).reduce((a, b) => a + b)
	return NLA.snapTo(doubleSignedArea * Math.sign(normal.e(absMaxDim)), 0)
}

/**
 * solves x² + px + q = 0
 *
 * @param {number} p
 * @param {number} q
 * @returns {Array.<number>}
 */
function pqFormula(p, q) {
	// 4 times the discriminant:in
	var discriminantX4 = p * p / 4 - q
	if (discriminantX4 < -NLA.PRECISION) {
		return []
	} else if (discriminantX4 <= NLA.PRECISION) {
		return [-p/2]
	} else {
		var root = Math.sqrt(discriminantX4)
		return [-p/2 - root, -p/2 + root]
	}
}

/**
 * solves ax³ + bx² + cx + d = 0
 *
 * @param {number} a
 * @param {number} b
 * @param {number} c
 * @param {number} d
 * @returns {number[]} with 0-3 entries
 */
function solveCubicReal2(a, b, c, d) {
	if (NLA.isZero(a)) {
		if (NLA.isZero(b)) {
			return [-d/c]
		} else {
			return pqFormula(c / b, d / b)
		}
	}
	var div = a
	a = b/div
	b = c/div
	c = d/div
	var p = (3*b - a*a)/3,
		pDiv3 = p/3,
		pDiv3Pow3 = pDiv3 * pDiv3 * pDiv3,
		q = (2*a*a*a - 9*a*b + 27*c)/27,
		qDiv2 = q/2,
		discriminant = qDiv2*qDiv2 + pDiv3Pow3,
		u1,v1,x1,x2,x3
	// 18abcd - 4b³d + b²c² - 4ac³ - 27a²d²
	if (discriminant < -NLA.PRECISION / 8) {
		var r = Math.sqrt(-pDiv3Pow3),
			t = -q/(2*r),
			cosphi = t<-1 ? -1 : t>1 ? 1 : t, // clamp t to [-1;1]
			phi = Math.acos(cosphi),
			t1 = 2 * Math.cbrt(r)
		x1 = t1 * Math.cos(phi/3) - a/3
		x2 = t1 * Math.cos((phi+2*Math.PI)/3) - a/3
		x3 = t1 * Math.cos((phi+4*Math.PI)/3) - a/3
		return [x1, x2, x3]
	} else if(discriminant <= NLA.PRECISION / 8) {
		if (0 == qDiv2) {
			// TODO: compare with NLA.isZero?
			return [-a/3]
		}
		u1 = Math.cbrt(Math.abs(qDiv2))
		x1 = 2*u1 - a/3
		x2 = -u1 - a/3
		return [x1,x2]
	} else {
		var sd = Math.sqrt(discriminant)
		u1 = Math.cbrt(-qDiv2+sd)
		v1 = Math.cbrt(qDiv2+sd)
		return [u1-v1-a/3]
	}
}
// function solveCubicReal2(a, b, c, d) {
// 	if (NLA.isZero(a)) {
// 		if (NLA.isZero(b)) {
// 			return [-d/c]
// 		} else {
// 			return pqFormula(c / b, d / b)
// 		}
// 	}
// 	var div = a
// 	a = b/div
// 	b = c/div
// 	c = d/div
// 	var p = (3*b - a*a)/3,
// 		p3 = p/3,
// 		q = (2*a*a*a - 9*a*b + 27*c)/27,
// 		q2 = q/2,
// 		discriminant = q2*q2 + p3*p3*p3,
// 		u1,v1,x1,x2,x3
// 	// 18abcd - 4b³d + b²c² - 4ac³ - 27a²d²
// 	if (discriminant < -NLA.PRECISION / 8) {
// 		var mp3 = -p/3,
// 			mp33 = mp3*mp3*mp3,
// 			r = Math.sqrt( mp33 ),
// 			t = -q/(2*r),
// 			cosphi = t<-1 ? -1 : t>1 ? 1 : t,
// 			phi = Math.acos(cosphi),
// 			crtr = Math.cbrt(r),
// 			t1 = 2*crtr
// 		x1 = t1 * Math.cos(phi/3) - a/3
// 		x2 = t1 * Math.cos((phi+2*Math.PI)/3) - a/3
// 		x3 = t1 * Math.cos((phi+4*Math.PI)/3) - a/3
// 		return [x1, x2, x3]
// 	} else if(discriminant <= NLA.PRECISION / 8) {
// 		if (0 == q2) {
// 			// TODO: compare with NLA.isZero?
// 			return [-a/3]
// 		}
// 		u1 = q2 < 0 ? Math.cbrt(-q2) : -Math.cbrt(q2)
// 		x1 = 2*u1-a/3
// 		x2 = -u1 - a/3
// 		return [x1,x2]
// 	} else {
// 		var sd = Math.sqrt(discriminant)
// 		u1 = Math.cbrt(-q2+sd)
// 		v1 = Math.cbrt(q2+sd)
// 		return [u1-v1-a/3]
// 	}
// }

/**
 * @template T
 */
class Transformable extends NLA.BaseObject {
	mirrored(plane) {
		return this.transform(M4.mirroring(plane));
	}
	mirroredX() {
		return this.mirrored(P3.YZ)
	}
	mirroredY() {
		return this.mirrored(P3.ZX)
	}
	mirroredZ() {
		return this.mirrored(P3.XY)
	}
	translate(x, y, z) {
		return this.transform(M4.translation(x, y, z), `.translate(${x}, ${y}, ${z})`)
	}
	scale(x, y, z) {
		return this.transform(M4.scaling(x, y, z))
	}
	rotateX(radians) {
		return this.transform(M4.rotationX(radians), `.rotateX(${radians})`)
	}
	rotateY(radians) {
		return this.transform(M4.rotationY(radians), `.rotateY(${radians})`)
	}
	rotateZ(radians) {
		return this.transform(M4.rotationZ(radians), `.rotateZ(${radians})`)
	}
	rotate(rotationCenter, rotationAxis, radians) {
		return this.transform(M4.rotation(rotationCenter, rotationAxis, radians))
	}
	eulerZXZ(alpha, beta, gamma) {
		return this.transform(M4.eulerZXZ(alpha, beta, gamma))
	}
	project(plane) {
		return this.transform(M4.projection(plane))
	}
	projXY() {
		return this.transform(M4.projection(P3.XY))
	}
	projYZ() {
		return this.transform(M4.projection(P3.YZ))
	}
	projZX() {
		return this.transform(M4.projection(P3.ZX))
	}
	shearedX (y, z) {
		return this.transform(M4([
			1, y, z, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1]))
	}
	foo () {
		return this.transform(M4.FOO)
	}
	bar () {
		return this.transform(M4.BAR)
	}

	/**
	 * @abstract
	 * @param {M4} m4
	 * @returns {T}
	 * @method
	 */
}
NLA.addOwnProperties(Array.prototype, ARRAY_UTILITIES)
