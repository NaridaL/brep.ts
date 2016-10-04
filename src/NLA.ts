"use strict"


// I'd use a loop but it kills the type checker.
const ceil = Math.ceil
const floor = Math.floor
const abs = Math.abs
const sign = Math.sign
const atan2 = Math.atan2
const atan = Math.atan
const cos = Math.cos
const sin = Math.sin
const min = Math.min
const max = Math.max
const PI = Math.PI
const sqrt = Math.sqrt
const pow = Math.pow
const round = Math.round
const log = Math.log

var globalId = 0;
type int = number
type FloatArray = Float32Array | Float64Array | number[]

/** @define {boolean} */
const NLA_DEBUG = true
const NLA_PRECISION = 1 / (1 << 26)
console.log("NLA_PRECISION", NLA_PRECISION)

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

interface Object {
	sce:string
	str:string
	toSource():string
}
interface Array<T> {
    withMax(f: (el: T)=>number): T
    mapFilter<U>(f: (el: T, elIndex: int, array: Array<T>)=>U): U[]
    concatenated(): T
    includes(el: T): boolean
    isEmpty(): boolean
    remove(el: T): boolean
    unique(): T[]
    last(): T
    pushAll(els: T[])
    flatMap<U>(f: (el: T) => (U | U[])): U[]
    sum(): number
	max(): number
	min(): number
}

var oldConsole = undefined
function disableConsole() {
	oldConsole = console.log
	console.log = function() {}
}
function enableConsole() {
	if (oldConsole) {
		console.log = oldConsole;
	}
}


function assertVectors(...vectors:(Vector|V3)[]) {
	if (NLA_DEBUG) {
		for (var i = 0; i < arguments.length; i++) {
			if (!(arguments[i] instanceof V3 || arguments[i] instanceof NLA.Vector)) {
				throw new Error("assertVectors arguments[" + (i) + "] is not a vector. " + typeof arguments[i] + " == typeof " + arguments[i]);
			}
		}
	}
	return true
}
function assertInst(what, ...objs) {
	if (NLA_DEBUG) {
		for (var i = 0; i < objs.length; i++) {
			if (!(objs[i] instanceof what)) {
				throw new Error("assertInst objs[" + (i) + "] is not a "+what.prototype.constructor.name+". " + objs[i].constructor.name + objs[i])
			}
		}
	}
	return true
}
function assertNumbers(...numbers) {
	if (NLA_DEBUG) {
		for (var i = 0; i < numbers.length; i++) {
			if (typeof numbers[i] !== 'number') {
				throw new Error("assertNumbers arguments[" + (i) + "] is not a number. " + typeof numbers[i] + " == typeof " + numbers[i]);
			}
		}
	}
	return true
}
function assert(value:any, ...messages:(any|(() => string))[]):boolean {
	if (NLA_DEBUG && !value) {
		throw new Error("NLA.assert failed: "
			+ messages.map(message => ('function' === typeof message ? message() : message || '')).join('\n'))
	}
	return true
}
function assertNever(value?: never): never {
    throw new Error()
}

function assertf(f:() => any, ...messages:(string|(() => any))[]) {
	if (!f()) {
		throw new Error("NLA.assertf failed: " + f.toString()
			+ messages.map(message => ('function' === typeof message ? message() : message || '')).join('\n'))
	}
}

namespace NLA {
    export const eq0 = (x: number): boolean => Math.abs(x) < NLA_PRECISION
    export const eq02 = (x: number, precision: number) => Math.abs(x) < precision
    export const eq = (x: number, y: number) => Math.abs(x - y) <= NLA_PRECISION
    export const lt = (x: number, y: number): boolean => x + NLA_PRECISION < y
    export const gt = (x: number, y: number): boolean => x > y + NLA_PRECISION
    export const le = (x: number, y: number): boolean => x <= y + NLA_PRECISION
    export const ge = (x: number, y: number): boolean => x + NLA_PRECISION >= y
    export const eq2 = (x, y, precision): boolean => Math.abs(x - y) < precision
    export const eqAngle = (x: number, y: number): boolean => zeroAngle(x - y)
    export const zeroAngle = (x:number):number => ((x % (2 * Math.PI)) + 2 * Math.PI + NLA_PRECISION) % (2 * Math.PI) < 2 * NLA_PRECISION
    export const snapTo = (x:number, to:number):number => Math.abs(x - to) <= NLA_PRECISION ? to : x
    export const canonAngle = (x:number):number => ((x % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI)


    /**
	 * Decimal adjustment of a number.
	 *
	 * @param  f  The type of adjustment.
	 * @param  value The number.
	 * @param exp   The exponent (the 10 logarithm of the adjustment base).
	 * @returns {number} The adjusted value.
	 */
	function decimalAdjust(f: (x: number)=>number, value: number, exp: number): number {
		// If the exp is undefined or zero...
		if (typeof exp === 'undefined' || +exp === 0) {
			return f(value)
		}
		value = +value
		exp = +exp
		// If the value is not a number or the exp is not an integer...
		if (isNaN(value) || !(typeof exp === 'number' && exp % 1 === 0)) {
			return NaN
		}
		// Shift
		let vs = value.toString().split('e')
		value = f(+(vs[0] + 'e' + (vs[1] ? (+vs[1] - exp) : -exp)))
		// Shift back
		vs = value.toString().split('e')
		return +(vs[0] + 'e' + (vs[1] ? (+vs[1] + exp) : exp))
	}

	export const round10: (value: number, exp: number) => number = decimalAdjust.bind(undefined, Math.round)
	export const floor10: (value: number, exp: number) => number = decimalAdjust.bind(undefined, Math.floor)
	export const ceil10: (value: number, exp: number) => number = decimalAdjust.bind(undefined, Math.ceil)





	export function repeatString(count, str) {
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
	export const mod = function (a, b) {
		return ((a % b) + b) % b
	}
	export const arraySwap = function (arr, i, j) {
		var temp = arr[i]
		arr[i] = arr[j]
		arr[j] = temp
	}
	export const arrayCopy = function(src,sstart,dst,dstart,length) {
		dstart += length;
		length += sstart;
		while(length-- > sstart) {
			dst[--dstart] = src[length];
		}
	}
	export const clamp = function (val, min, max) {
		assertNumbers(val, min, max)
		return Math.max(min, Math.min(max, val))
	}

	export const randomColor = function () {
		return Math.floor(Math.random() * 0x1000000)
	}
	export function mapAdd<T, U>(map:Map<T, U[]>, key:T, val:U) {
		var list = map.get(key)
		if (list) {
			list.push(val)
		} else {
			map.set(key, [val])
		}
	}


	export const arrayCopyStep = function(src,sstart,sstep, dst,dstart,dstep,count) {
		var srcIndex = sstart + count * sstep;
		var dIndex = dstart + count * dstep;
		while(srcIndex > sstart) {
			dst[dIndex -= dstep] = src[srcIndex -= sstep];
		}
	}
	export const arrayCopyBlocks = function (src,sstart,sstep, dst,dstart,dstep,blockSize, blockCount) {
		for (var i = 0; i < blockCount; i++) {
			NLA.arrayCopy(src, sstart + sstep * i, dst, dstart + dstep * i, blockSize)
		}
	}
	export const arrayRange = function (start, end, step) {
		assertNumbers(start, step)
		//console.log(Math.ceil((end - start) / step))
		var result = new Array(Math.ceil((end - start) / step)); // "- start" so that chunk in the last row will also be selected, even if the row is not complete
		var index = 0
		for (var i = 0; i < end; i += step) {
			result[index++] = i
		}
		return result;
	}

	export const arrayFromFunction = function<T>(length:number, f:(i:number) => T):T[] {
		assertNumbers(length)
		assert("function" == typeof f)
		var a = new Array(length)
		var elIndex = length;
		while (elIndex--) {
			a[elIndex] = f(elIndex);
		}
		return a
	}

	/**
	 * @param vals
	 * @returns {number[]}
	 */
	export const fuzzyUniques = function (vals:number[]):number[] {
		let round = val => Math.floor(val * (1 << 26)) / (1 << 26)
		let map = new Map()
		map.set(round(vals[0]), vals[0])
		for (let i = 1; i < vals.length; i++) {
			let val = vals[i], roundVal = round(val)
			let key
			if (!map.has(roundVal)
				&& !((key = map.get(roundVal - 1 / (1 << 26))) && NLA.eq(key, val))
				&& !((key = map.get(roundVal + 1 / (1 << 26))) && NLA.eq(key, val))) {
				map.set(roundVal, val)
			}
		}
		return Array.from(map.values())
	}

	/**
	 * @param vals
	 * @param f
	 * @returns {T[]}
	 */
	export const fuzzyUniquesF = function<T>(vals:T[], f: (o:T) => number) {
		let round = val => Math.floor(val * (1 << 26)) / (1 << 26)
		let map = new Map()
		for (let i = 0; i < vals.length; i++) {
			let val = vals[i], roundVal = round(f(val))
			let key
			if (!map.has(roundVal)
				&& !((key = map.get(roundVal - 1 / (1 << 26))) && NLA.eq(key, f(val)))
				&& !((key = map.get(roundVal + 1 / (1 << 26))) && NLA.eq(key, f(val)))) {
				map.set(roundVal, val)
			}
		}
		return Array.from(map.values())
	}


	export const addOwnProperties = function (target, props) {
		Object.getOwnPropertyNames(props).forEach(key => {
			//console.log(props, key)
			if (target.hasOwnProperty(key)) {
				console.warn("target ", target, " already has property ", key)
			}
			Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(props, key))
		})
	}

	export const defineClass = function (name, parent, constructor, props, statics) {
		assert('function' == typeof constructor, "'function' == typeof constructor")
		constructor.prototype = NLA.defineObject(parent && parent.prototype, props)
		constructor.prototype.constructor = constructor
		Object.defineProperty(constructor.prototype, 'name', {value: name})
		statics && NLA.addOwnProperties(constructor, statics)
		return constructor
	}
	export const defaultRoundFunction = x => x // Math.round10(x, -4)

	export const defineObject = function (prot, props) {
		var o = Object.create(prot || NLA.BaseObject)
		addOwnProperties(o, props)
		return o
	}
	export const BaseObject = class extends Object {
		toSource() {
			return this.toString == Object.prototype.toString ? Object.prototype.toSource.apply(this) : this.toString()
		}
	}
	export const forceFinite = function(val) {
		val = parseFloat(val.replace(',', '.').replace(/^[^0-9,\.\-]/, ''))
		return Number.isFinite(val) ? val : 0
	}

	export const minus = (a, b) => a - b





	/**
	 * combinations(2) will generate
	 * [0,0] [0,1] [1,1] [0,2] [1,2] [2,2]
	 * @param n
	 * @returns {{i: number, j: number}}
	 */
	export const combinations = function* (n:number) {
		for (var i = 0; i < n; i++) {
			for (var j = i; j < n; j++) {
				yield {i: i, j: j}
			}
		}
	}


	export const CLASSES:any = {}
	export const registerClass = function (clazz: any) {
		NLA.CLASSES[clazz.name] = clazz
	}
}
/**
 *
 * @param x
 * @returns {boolean}
 */

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

interface String {
	capitalizeFirstLetter(): string
}
String.prototype.capitalizeFirstLetter = function() {
	return this.charAt(0).toUpperCase() + this.slice(1);
}
var ARRAY_UTILITIES =
/** @template T
@lends Array.prototype */ {
	pushAll: function (arr) {
		Array.prototype.push.apply(this, arr)
	},
	sliceStep: function (start, step, chunkSize) {
		assertNumbers(start, step)
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
	 * @param f
	 * @returns {Array.<S>}
	 */
	mapFilter<U>(f:(T) => U) {
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
	flatMap(f) {
		return Array.prototype.concat.apply([], this.map(f));
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
	 * @param searchElement
	 * @param cmp
	 * @returns {number}
	 */
	binaryIndexOf: function<S>(searchElement:S, cmp:(a:T, b:S) => number) {

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
for (let key in ARRAY_UTILITIES) {
	NLA["array" + key.capitalizeFirstLetter()] = function (arr, ...rest) {
		assert(!ARRAY_UTILITIES[key])
		ARRAY_UTILITIES[key].apply(arr, rest)
	}
}
// Closure

	;(function() {
})();



function isCCW(vertices, normal) {
	var dsa = doubleSignedArea(vertices, normal)
	assert(0 != dsa)
	return dsa < 0
}
declare function earcut(data:number[], holeIndices:number[], dim?:number): number[]
function triangulateVertices(normal, vertices, holeStarts) {
	var absMaxDim = normal.maxAbsDim(), factor = normal.e(absMaxDim) < 0 ? -1 : 1
	var contour = new Array(vertices.length * 2)
	var i = vertices.length
	/*
	 var [coord0, coord1] = [['y', 'z'], ['z', 'x'], ['x', 'y']][maxAbsDim]
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
	var absMaxDim = normal.maxAbsDim()
	// order is important, coord0 and coord1 must be set so that coord0, coord1 and maxDim span a right-hand coordinate system
	//var [coord0, coord1] = [['y', 'z'], ['z', 'x'], ['x', 'y']][maxAbsDim]
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
 * @param p
 * @param q
 * @returns {Array.<number>}
 */
function pqFormula(p:number, q:number) {
	// 4 times the discriminant:in
	var discriminantX4 = p * p / 4 - q
	if (discriminantX4 < -NLA_PRECISION) {
		return []
	} else if (discriminantX4 <= NLA_PRECISION) {
		return [-p/2]
	} else {
		var root = Math.sqrt(discriminantX4)
		return [-p/2 - root, -p/2 + root]
	}
}

/**
 * solves ax³ + bx² + cx + d = 0
 *
 * @param a
 * @param b
 * @param c
 * @param d
 * @returns {number[]} with 0-3 entries
 */
function solveCubicReal2(a:number, b:number, c:number, d:number) {
	if (NLA.eq0(a)) {
		if (NLA.eq0(b)) {
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
	if (discriminant < -NLA_PRECISION / 8) {
		var r = Math.sqrt(-pDiv3Pow3),
			t = -q/(2*r),
			cosphi = t<-1 ? -1 : t>1 ? 1 : t, // clamp t to [-1;1]
			phi = Math.acos(cosphi),
			t1 = 2 * Math.cbrt(r)
		x1 = t1 * Math.cos(phi/3) - a/3
		x2 = t1 * Math.cos((phi+2*Math.PI)/3) - a/3
		x3 = t1 * Math.cos((phi+4*Math.PI)/3) - a/3
		return [x1, x2, x3]
	} else if(discriminant <= NLA_PRECISION / 8) {
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

NLA.addOwnProperties(Array.prototype, ARRAY_UTILITIES)

let a = [1, 2, 3].mapFilter(x => x * 2)












/**
 *
 * @param f
 * @param xStart
 * @param steps
 * @param EPSILON
 * @returns {*}
 */
function newtonIterate(f: (x: number[]) => number[], xStart: number[], steps?: number, EPSILON?: number) {
	steps = steps || 4
	EPSILON = EPSILON || 1e-8

	let x = xStart
	for (let i = 0; i < steps; i++) {
		let fx = f(x)
		let dfdx = Matrix.jacobi(f, x, fx, EPSILON)
		assert(!dfdx.isZero())
		let dx = dfdx.solveLinearSystem(new NLA.Vector(new Float64Array(fx))).v
		assert(!isNaN(dx[0]))
		//console.log("fx / dfdx", fx / dfdx)
		for (let j = 0; j < x.length; j++) x[j] -= dx[j]
	}
	return x
}

function newtonIterate1d(f:(x:number) => number, xStart: number, steps?:number, EPSILON?:number):number {
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
function newtonIterateWithDerivative(f:(x:number) => number, xStart:number, steps:number, df:(x:number)=>number) {
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
 * @param f1
 * @param f2
 * @param sStart
 * @param tStart
 * @param {number=} steps
 * @returns {V3}
 */
function newtonIterate2d(f1:(s:number,t:number)=>number, f2:(s:number,t:number)=>number, sStart:number, tStart:number, steps?:number) {
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
	} while (--steps && f1ts * f1ts + f2ts * f2ts > NLA_PRECISION)
	if (!steps) {
		//console.log(f1ts * f1ts + f2ts * f2ts)
		return null
	}
	return V(s, t, 0)
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
	} while (--steps && f1ts * f1ts + f2ts * f2ts > NLA_PRECISION / 32)
	if (!steps) {
		//console.log(f1ts * f1ts + f2ts * f2ts)
		return null
	}
	return V(s, t, 0)
}