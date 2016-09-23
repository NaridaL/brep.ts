"use strict"
/**
 * @namespace
 */
var NLA = {}
window['NLA'] = NLA


/** @define {boolean} */
const NLA_DEBUG = true

function SCE(o) {
	return o.sce
}
function STR(o) {
	return o.str
}
Object.defineProperty(Object.prototype, 'sce', {get: function () { return this.toSource() }})
Object.defineProperty(Object.prototype, 'str', {get: function () { return this.toString() }})
if (!Object.prototype.toSource) {
	Object.defineProperty(Object.prototype, 'toSource', {value: function () { return JSON.stringify(this) }})
}


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
var assert = NLA.assert = function (value, ...messages) {
	if (NLA_DEBUG && !value) {
		throw new Error("NLA.assert failed: " + messages.map(message => ('function' == typeof message ? message() : message)).join('\n'))
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
	firstMatch: function (f) {
		for (let i = 0; i < this.length; i++) {
			if (i in this) {
				let val = f(this[i])
				if (val) {
					return val
				}
			}
		}
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
 * combinations(2) will generate
 * [0,0] [0,1] [1,1] [0,2] [1,2] [2,2]
 * @param {number} n
 * @returns {{i: number, j: number}}
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

/**
 * @template T
 * @class Transformable.<T>
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
	 * @method transform
	 */
}
NLA.addOwnProperties(Array.prototype, ARRAY_UTILITIES)












/**
 *
 * @param {function(Array.<number>):Array.<number>} f
 * @param {Array.<number>} xStart
 * @param {number=} steps
 * @param {number} EPSILON
 * @returns {*}
 */
function newtonIterate(f, xStart, steps, EPSILON) {
	steps = steps || 4
	EPSILON = EPSILON || 1e-8

	let x = xStart
	for (let i = 0; i < steps; i++) {
		let fx = f(x)
		let dfdx = NLA.Matrix.jacobi(f, x, fx, EPSILON)
		assert(!dfdx.isZero())
		let dx = dfdx.solveLinearSystem(new NLA.Vector(new Float64Array(fx))).v
		assert (!isNaN(dx[0]))
		//console.log("fx / dfdx", fx / dfdx)
		for (let j = 0; j < x.length; j++) x[j] -= dx[j]
	}
	return x
}

/**
 *
 * @param {function(number):number} f
 * @param {number} xStart
 * @param {number=} steps int
 * @param {number=} EPSILON
 * @returns {*}
 */
function newtonIterate1d(f, xStart, steps, EPSILON) {
	steps = steps || 4
	EPSILON = EPSILON || 1e-8

	let x = xStart

	for (let i = 0; i < steps; i++) {
		let fx = f(x)
		let dfdx = (f(x + EPSILON) - fx) / EPSILON
		//console.log("fx / dfdx", fx / dfdx)
		x = x - fx / dfdx
	}
	return x
}
function newtonIterateWithDerivative(f, xStart, steps, df) {
	steps = steps || 4
	let x = xStart
	for (let i = 0; i < steps; i++) {
		let fx = f(x)
		let dfdx = df(x)
		//console.log("fx / dfdx", fx / dfdx)
		x = x - fx / dfdx
	}
	return x
}
/**
 *
 * @param {function (number, number):number} f1
 * @param {function (number, number):number} f2
 * @param {number} sStart
 * @param {number} tStart
 * @param {number=} steps
 * @returns {V3}
 */
function newtonIterate2d(f1, f2, sStart, tStart, steps) {
	const EPSILON = 1e-6
	steps = steps || 4
	let s = sStart, t = tStart
	do {
		/*
		 | a b |-1                   |  d -b |
		 | c d |   = 1 / (ad - bc) * | -c  a |
		 */
		var f1ts = f1(s, t), f2ts = f2(s, t)
		/*
		 let df1s = (f1(s + EPSILON, t) - f1ts) / EPSILON, df1t = (f1(s, t + EPSILON) - f1ts) / EPSILON,
		 df2s = (f2(s + EPSILON, t) - f2ts) / EPSILON, df2t = (f2(s, t + EPSILON) - f2ts) / EPSILON
		 let det = df1s * df2t - df1t * df2s
		 s = s - ( df2t * f1ts - df1t * f2ts) / det
		 t = t - (-df2s * f1ts + df1s * f2ts) / det
		 */
		// TODO: is this even more accurate?
		let df1s = (f1(s + EPSILON, t) - f1ts), df1t = (f1(s, t + EPSILON) - f1ts),
			df2s = (f2(s + EPSILON, t) - f2ts), df2t = (f2(s, t + EPSILON) - f2ts)
		let det = (df1s * df2t - df1t * df2s) / EPSILON
		let ds = ( df2t * f1ts - df1t * f2ts) / det
		let dt = (-df2s * f1ts + df1s * f2ts) / det
		s -= ds
		t -= dt
	} while (--steps && f1ts * f1ts + f2ts * f2ts > NLA.PRECISION)
	if (!steps) {
		//console.log(f1ts * f1ts + f2ts * f2ts)
		return null
	}
	return V3(s, t, 0)
}
function newtonIterate2dWithDerivatives(f, g, sStart, tStart, steps, dfds, dfdt, dgds, dgdt) {
	steps = steps || 4
	let s = sStart, t = tStart
	let eps = 1e-6
	do {
		/*
		 | a b |-1                   |  d -b |
		 | c d |   = 1 / (ad - bc) * | -c  a |
		 */
		var f1ts = f(s, t), f2ts = g(s, t)
		let df1s = dfds(s, t), df1t = dfdt(s, t),
			df2s = dgds(s, t), df2t = dgdt(s, t)
		// TODO: is this even more accurate?
		let det = df1s * df2t - df1t * df2s
		let ds = ( df2t * f1ts - df1t * f2ts) / det
		let dt = (-df2s * f1ts + df1s * f2ts) / det
		s -= ds
		t -= dt
	} while (--steps && f1ts * f1ts + f2ts * f2ts > NLA.PRECISION / 32)
	if (!steps) {
		//console.log(f1ts * f1ts + f2ts * f2ts)
		return null
	}
	return V3(s, t, 0)
}