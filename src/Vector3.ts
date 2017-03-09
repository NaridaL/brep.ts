import Equalable = NLA.Equalable
class V3 implements Equalable {
	readonly x: number
	readonly y: number
	readonly z: number

	constructor(x: number, y: number, z: number) {
		//assert(arguments.length == 3)
		//assertNumbers(x, y, z)
		//Object.defineProperties(this, {
		//	x: {value: x},
		//	y: {value: y},
		//	z: {value: z},
		//	//id: {value: V3.id++}
		//})
        this.x = x
        this.y = y
        this.z = z
	}


	perturbed(delta?: number): V3 {
		delta = delta || NLA_PRECISION * 0.8
		return this.map(x => x + (Math.random() - 0.5) * delta)
	}


	e(index: int): number {
		assert(index >= 0 && index < 3)
		return 0 == index ? this.x : (1 == index ? this.y : this.z)
	}

	negated(): V3 {
		return new V3(-this.x, -this.y, -this.z)
	}

	abs(): V3 {
		return new V3(Math.abs(this.x), Math.abs(this.y), Math.abs(this.z))
	}

	plus(a: V3): V3 {
		assertVectors(a)
		return new V3(this.x + a.x, this.y + a.y, this.z + a.z)
	}

	scale(a: V3): V3 {
		return new V3(this.x * a.x, this.y * a.y, this.z * a.z)
	}

	/**
	 * See also {@link to} which is a.minus(this)
	 */
	minus(a: V3): V3 {
		assertVectors(a)
		return new V3(this.x - a.x, this.y - a.y, this.z - a.z)
	}

	to(a: V3): V3 {
		assertVectors(a)
		return a.minus(this)
	}

	times(factor: number): V3 {
		assertNumbers(factor)
		return new V3(this.x * factor, this.y * factor, this.z * factor)
	}

	div(a: number): V3 {
		assertNumbers(a)
		return new V3(this.x / a, this.y / a, this.z / a)
	}

	dot(a: V3): number {
		assertInst(V3, a)
		return this.x * a.x + this.y * a.y + this.z * a.z
	}

	/**
	 * Linearly interpolate
	 */
	lerp(b: V3, t: number): V3 {
		assertVectors(b)
		assertNumbers(t)
		return this.plus(b.minus(this).times(t))
	}

	squared(): number {
		return this.dot(this)
	}

	distanceTo(a: V3): number {
		assertVectors(a)
		//return this.minus(a).length()
		return Math.hypot(this.x - a.x, this.y - a.y, this.z - a.z)
	}

	distanceToSquared(a: V3): number {
		assertVectors(a)
		return this.minus(a).squared()
	}

	toSource(): string {
		return V3.NAMEMAP.get(this) || this.toString()
	}

	nonParallelVector(): V3 {
		const abs = this.abs()
		if ((abs.x <= abs.y) && (abs.x <= abs.z)) {
			return V3.X
		}
		else if ((abs.y <= abs.x) && (abs.y <= abs.z)) {
			return V3.Y
		}
		else {
			return V3.Z
		}
	}

	slerp(b: V3, t: number): V3 {
		assertVectors(b)
		assertNumbers(t)
		const sin = Math.sin
		const omega = this.angleTo(b)
		return this.times(sin((1 - t) * omega) / sin(omega)).plus(b.times(sin(t * omega) / sin(omega)))
	}

	min(b: V3): V3 {
		return new V3(Math.min(this.x, b.x), Math.min(this.y, b.y), Math.min(this.z, b.z))
	}

	max(b: V3): V3 {
		return new V3(Math.max(this.x, b.x), Math.max(this.y, b.y), Math.max(this.z, b.z))
	}

	equals(v: any): boolean {
		return this == v || this.x == v.x && this.y == v.y && this.z == v.z
	}

	/**
	 *
	 * The cross product is defined as:
	 * a x b = |a| * |b| * sin(phi) * n
	 * where |.| is the euclidean norm, phi is the angle between the vectors
	 * and n is a unit vector perpendicular to both a and b.
	 *
	 * The cross product is zero for parallel vectors.
	 */
	cross(v: V3): V3 {
		return new V3(this.y * v.z - this.z * v.y, this.z * v.x - this.x * v.z, this.x * v.y - this.y * v.x)
	}

	/**
	 * Documentation stub. You want {@link normalized}
	 */
	unit(): V3 { throw new Error() }

	minElement(): number {
		return Math.min(this.x, this.y, this.z)
	}

	maxElement(): number {
		return Math.max(this.x, this.y, this.z)
	}

	toArray(n: int = 3): number[] {
		return [this.x, this.y, this.z].slice(0, n)
	}

	/**
	 * Get a perpendicular vector.
	 * For vectors in the XY-Plane, returns vector rotated 90° CCW.
	 */
	getPerpendicular(): V3 {
		if (NLA.eq0(this.x) && NLA.eq0(this.y)) {
			if (NLA.eq0(this.z)) {
				throw new Error('zero vector')
			}
			// v is Vector(0, 0, v.z)
			return V3.Y
		}
		return new V3(-this.y, this.x, 0)
	}

	dim(): int {
		return 3
	}

	els(): number[] {
		return [this.x, this.y, this.z]
	}

	angleXY(): number {
		return Math.atan2(this.y, this.x)
	}

	lengthXY(): number {
		return Math.hypot(this.x, this.y)
		//return Math.sqrt(this.x * this.x + this.y * this.y)
	}

	squaredXY(): number {
		return this.x * this.x + this.y * this.y
	}

	/**
	 * Transform this vector element-wise by way of function f. Returns V3(f(x), f(y), f(z))
	 * @param f function to apply to elements (number -> number)
	 */
	map(f: (el: number) => number): V3 {
		return new V3(f(this.x), f(this.y), f(this.z))
	}

	toString(roundFunction?): string {
		roundFunction = roundFunction || NLA.defaultRoundFunction
		return 'V(' + [this.x, this.y, this.z].map(roundFunction).join(', ') + ')' //+ this.id
	}

	angleTo(b: V3): number {
		assert(1 == arguments.length)
		assertVectors(b)
		assert(!this.isZero())
		assert(!b.isZero())
		return Math.acos(this.dot(b) / this.length() / b.length())
	}

	/**
	 *
	 * phi = angle between A and B
	 * alpha = angle between n and normal1
	 *
	 * A . B = ||A|| * ||B|| * cos(phi)
	 * A x B = ||A|| * ||B|| * sin(phi) * n (n = unit vector perpendicular)
	 * (A x B) . normal1 = ||A|| * ||B|| * sin(phi) * cos(alpha)
	 */
	angleRelativeNormal(vector: V3, normal1: V3): number {
		assertf(() => 2 == arguments.length)
		assertVectors(vector, normal1)
		assertf(() => normal1.hasLength(1))
		assert(vector.isPerpendicularTo(normal1), 'vector.isPerpendicularTo(normal1)' + vector.sce + normal1.sce)
		assert(this.isPerpendicularTo(normal1), 'this.isPerpendicularTo(normal1)' + this.dot(vector))
		return Math.atan2(this.cross(vector).dot(normal1), this.dot(vector))
	}

	/**
	 Returns true iff this is parallel to vector, i.e. this * s == vector, where s is a pos or neg number, using NLA.equals
	 Throw a DebugError
	 if vector is not a NLA.NLA.Vector or
	 if this has a length of 0 or
	 if vector has a length of 0
	 */
	isParallelTo(vector: V3): boolean {
		assertVectors(vector)
		assert(!this.isZero())
		assert(!vector.isZero())
		// a . b takes on values of +|a|*|b| (vectors same direction) to -|a|*|b| (opposite direction)
		// in both cases the vectors are parallel, so check if abs(a . b) == |a|*|b|
		const dot = this.dot(vector)
		return NLA.eq(this.squared() * vector.squared(), dot * dot)
	}

	isPerpendicularTo(vector: V3): boolean {
		assertVectors(vector)
		assert(!this.isZero(), '!this.isZero()')
		assert(!vector.isZero(), '!vector.isZero()')
		return NLA.eq0(this.dot(vector))
	}

	isReverseDirTo(other: V3): boolean {
		assertVectors(other)
		assert(!this.isZero())
		assert(!other.isZero())
		// a . b takes on values of +|a|*|b| (vectors same direction) to -|a|*|b| (opposite direction)
		// in both cases the vectors are parallel, so check if abs(a . b) == |a|*|b|
		const dot = this.dot(other)
		return NLA.eq(Math.sqrt(this.squared() * other.squared()), dot)
	}

	/**
	 * Returns the length of this NLA.Vector, i.e. the euclidean norm.
	 *
	 * Note that the partial derivatives of the euclidean norm at point x are equal to the
	 * components of the normalized vector x.
	 */
	length(): number {
		return Math.hypot(this.x, this.y, this.z)
		//return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z)
	}

	/**
	 * Definition: V3.isZero == V3.like(V3.ZERO)
	 */
	isZero(): boolean {
		return this.like(V3.ZERO)
	}

	like(obj): boolean {
		if (obj === this)
			return true
		if (!(obj instanceof V3))
			return false
		return NLA.eq(this.x, obj.x) && NLA.eq(this.y, obj.y) && NLA.eq(this.z, obj.z)
	}

	/**
	 * equivalent to this.like(v) || this.negated().like(v)
	 */
	likeOrReversed(v: V3): boolean {
		return NLA.eq(Math.abs(this.dot(v)), Math.sqrt(this.squared() * v.squared()))
	}

	/**
	 * Returns a new unit NLA.NLA.Vector (.length() === 1) with the same direction as this vector. Throws a
	 * NLA.DebugError if this has a length of 0.
	 */
	normalized(): V3 {
		assert(!this.isZero(), 'cannot normalize zero vector')
		return this.div(this.length())
	}

	/**
	 * Returns a new V3 equal to this scaled so that its length is equal to newLength.
	 *
	 * Passing a negative newLength will flip the vector.
	 */
	toLength(newLength: number): V3 {
		assertNumbers(newLength)
		return this.times(newLength / this.length())
	}

	///**
	// * See also {@see #setTo} for the individual
	// *
	// * @param v
	// */
	//assign(v) {
	//	assertVectors(v)
	//	this.x = v.x
	//	this.y = v.y
	//	this.z = v.z
	//}
	//
	///**
	// * See also {@see #assign} for the V3 version
	// *
	// * @param x
	// * @param y
	// * @param z
	// */
	//setTo(x, y, z = 0) {
	//	this.x = x
	//	this.y = y
	//	this.z = z
	//}

	/**
	 Returns a new NLA.NLA.Vector which is the projection of this vector onto the passed vector.
	 Examples
	 NLA.V(3, 4).projectedOn(NLA.V(1, 0)) // returns NLA.V(3, 0)
	 NLA.V(3, 4).projectedOn(NLA.V(2, 0)) // returns NLA.V(3, 0)
	 NLA.V(3, 4).projectedOn(NLA.V(-1, 0)) // returns NLA.V(-3, 0)
	 NLA.V(3, 4).projectedOn(NLA.V(0, 1)) // returns NLA.V(0, 4)
	 NLA.V(3, 4).projectedOn(NLA.V(1, 1)) // returns
	 */
	projectedOn(b: V3): V3 {
		assertVectors(b)
		// https://en.wikipedia.org/wiki/NLA.Vector_projection#NLA.Vector_projection_2
		return b.times(this.dot(b) / b.dot(b))
	}

	rejectedFrom(b: V3): V3 {
		assertVectors(b)
		// https://en.wikipedia.org/wiki/NLA.Vector_projection#NLA.Vector_projection_2
		return this.minus(b.times(this.dot(b) / b.dot(b)))
	}

	rejectedFrom1(b1: V3): V3 {
		assertVectors(b1)
		assert(b1.hasLength(1))
		// https://en.wikipedia.org/wiki/NLA.Vector_projection#NLA.Vector_projection_2
		return this.minus(b1.times(this.dot(b1)))
	}

	/**
	 * Returns the length of this vector rejected from the unit vector b.
	 *
	 *       /|
	 * this / |    ^
	 *     /__|    | b
	 *      r
	 *  Returns length of r (r === this.rejectedFrom(b))
	 */
	rejectedLength(b: V3): number {
		assertVectors(b)
		return Math.sqrt(this.dot(this) - this.dot(b) ** 2 / b.dot(b))
	}

	/**
	 * Returns the length of this vector rejected from the unit vector b1.
	 *
	 *       /|
	 * this / |    ^
	 *     /__|    | b1
	 *      r
	 *  Returns length of r (r === this.rejectedFrom(b1))
	 */
	rejected1Length(b1: V3): number {
		assertVectors(b1)
		assert(b1.hasLength(1))
		return Math.sqrt(this.dot(this) - this.dot(b1) ** 2)
	}

	/**
	 Returns true iff the length() of this vector is equal to 'length', using NLA.eq
	 E.g. NLA.V(3, 4).hasLength(5) === true
	 NLA.V(1, 1).hasLength(1) === false
	 */
	hasLength(length: number): boolean {
		assertNumbers(length)
		return NLA.eq(length, this.length())
	}

	/**
	 Returns the sum of the absolute values of the components of this vector.
	 E.g. NLA.V(1, -2, 3) === abs(1) + abs(-2) + abs(3) === 1 + 2 + 3 === 6
	 */
	absSum(): number {
		return Math.abs(this.x) + Math.abs(this.y) + Math.abs(this.z)
	}

	/**
	 * returns max(|x|, |y|, |z|)
	 */
	maxAbsElement(): number {
		return Math.max(Math.abs(this.x), Math.abs(this.y), Math.abs(this.z))
	}

	/**
	 * returns min(|x|, |y|, |z|)
	 */
	minAbsElement(): number {
		return Math.min(Math.abs(this.x), Math.abs(this.y), Math.min(this.z))
	}

	maxAbsDim(): number {
		const xAbs = Math.abs(this.x), yAbs = Math.abs(this.y), zAbs = Math.abs(this.z)
		return xAbs >= yAbs ? (xAbs >= zAbs ? 0 : 2) : (yAbs >= zAbs ? 1 : 2)
	}

	minAbsDim(): number {
		const xAbs = Math.abs(this.x), yAbs = Math.abs(this.y), zAbs = Math.abs(this.z)
		return xAbs < yAbs ? (xAbs < zAbs ? 0 : 2) : (yAbs < zAbs ? 1 : 2)
	}

	withElement(dim: 'x' | 'y' | 'z', el: number) {
		assert(['x', 'y', 'z'].includes(dim), '' + dim)
		assertNumbers(el)
		if ('x' == dim) {
			return new V3(el, this.y, this.z)
		}
		if ('y' == dim) {
			return new V3(this.x, el, this.z)
		}
		return new V3(this.x, this.y, el)
	}

	hashCode(): int {
		function floatHashCode(f) {
			return ~~(f * (1 << 28))
		}

		return ~~((floatHashCode(this.x) * 31 + floatHashCode(this.y)) * 31 + floatHashCode(this.z))
	}

	hashCodes() {
		//function floatHashCode(f) {
		//	return ~~(f * (1 << 28))
		//}

		// compare hashCode.floatHashCode
		// the following ops are equivalent to
		// floatHashCode((el - NLA_PRECISION) % (2 * NLA_PRECISION))
		// this results in the hashCode for the (out of 8 possible) cube with the lowest hashCode
		// the other 7 can be calculated by adding constants
		const xHC = ~~(this.x * (1 << 28) - 0.5), 
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
	}

	compareTo(other: V3): number {
		if (this.x != other.x) {
			return this.x - other.x
		}
		else if (this.y != other.y) {
			return this.y - other.y
		}
		else {
			return this.z - other.z
		}
	}

	compareTo2(other: V3, precision: number = 0): number {
		if (!NLA.eq2(this.x, other.x, precision)) {
			return this.x - other.x
		}
		else if (!NLA.eq2(this.y, other.y, precision)) {
			return this.y - other.y
		}
		else if (!NLA.eq2(this.z, other.z, precision)) {
			return this.z - other.z
		}
		else {
			return 0
		}
	}


	static readonly ZERO: V3 = new V3(0, 0, 0)
	static readonly ONES: V3 = new V3(1, 1, 1)
    static readonly X: V3 = new V3(1, 0, 0)
    static readonly Y: V3 = new V3(0, 1, 0)
	static readonly Z: V3 = new V3(0, 0, 1)
    static readonly XY: V3 = new V3(1, 1, 0)
    static readonly INF: V3 = new V3(Infinity, Infinity, Infinity)
	static readonly XYZ: V3[] = [V3.X, V3.Y, V3.Z]
	
	static readonly NAMEMAP = new NLA.CustomMap<V3, string>()
		.set(V3.ZERO, 'V3.ZERO')
		.set(V3.ONES, 'V3.ONES')
		.set(V3.X, 'V3.X')
		.set(V3.Y, 'V3.Y')
		.set(V3.Z, 'V3.Z')
		.set(V3.INF, 'V3.INF')

	static random(): V3 {
		return new V3(Math.random(), Math.random(), Math.random())
	}

	static parallelCondition(a, b) {
		return a.dot(b) - a.length() * b.length()
	}

	/**
	 * See http://math.stackexchange.com/questions/44689/how-to-find-a-random-axis-or-unit-vector-in-3d
	 * @returns A random point on the unit sphere with uniform distribution across the surface.
	 */
	static randomUnit(): V3 {
		const zRotation = Math.random() * 2 * Math.PI
		const z = Math.random() * 2 - 1
		const zRadius = Math.sqrt(1 - z * z)
		return new V3(zRadius * Math.cos(zRotation), zRadius * Math.sin(zRotation), z)
	}

	static fromAngles(theta: number, phi: number): V3 {
		return new V3(Math.cos(theta) * Math.cos(phi), Math.sin(phi), Math.sin(theta) * Math.cos(phi))
	}

	static fromFunction(f: (dim: number) => number) {
		return new V3(f(0), f(1), f(2))
	}

	static min(a: V3, b: V3): V3 {
		return new V3(Math.min(a.x, b.x), Math.min(a.y, b.y), Math.min(a.z, b.z))
	}

	static max(a: V3, b: V3): V3 {
		return new V3(Math.max(a.x, b.x), Math.max(a.y, b.y), Math.max(a.z, b.z))
	}

	static lerp(a: V3, b: V3, fraction: number): V3 {
		return b.minus(a).times(fraction).plus(a)
	}

	static fromArray(a: number[]): V3 {
		return new V3(a[0], a[1], a[2])
	}

	static angleBetween(a: V3, b: V3): number {
		return a.angleTo(b)
	}

	static zip(f: (...args: number[]) => number, ...args):V3 {
		assert(f instanceof Function)
		return new V3(
			f.apply(undefined, args.map(x => x.x)),
			f.apply(undefined, args.map(x => x.y)),
			f.apply(undefined, args.map(x => x.z)))
	}

	static normalOnPoints(v0: V3, v1: V3, v2: V3): V3 {
		assertVectors(v0, v1, v2)
		return v1.minus(v0).cross(v2.minus(v0))
	}

	static add(...vs: V3[]): V3 {
		assertVectors.apply(undefined, vs)
		let x = 0, y = 0, z = 0
		let i = vs.length
		while (i--) {
			x += vs[i].x
			y += vs[i].y
			z += vs[i].z
		}
		return new V3(x, y, z)
	}

	static sub(...vs: V3[]): V3 {
		assertVectors.apply(undefined, vs)
		let x = vs[0].x, y = vs[0].y, z = vs[0].z
		let i = vs.length
		while (i--) {
			x -= vs[i].x
			y -= vs[i].y
			z -= vs[i].z
		}
		return new V3(x, y, z)
	}

	static flattenV3Array(v3arr: V3[], dest?: number[]|Float64Array|Float32Array, srcStart?: number, destStart?: number, v3count?: number) {
		//assert (v3arr.every(v3 => v3 instanceof V3), 'v3arr.every(v3 => v3 instanceof V3)')
		srcStart = srcStart || 0
		destStart = destStart || 0
		v3count = v3count || (v3arr.length - srcStart)
		dest = dest || new Float32Array(3 * v3count)
		assert(dest.length - destStart >= v3count, 'dest.length - destStart >= v3count')
		let i = v3count
		while (i--) {
			const v = v3arr[srcStart + i]
			dest[destStart + i * 3] = v.x
			dest[destStart + i * 3 + 1] = v.y
			dest[destStart + i * 3 + 2] = v.z
		}
		return dest
	}

	static perturbed(v: V3, delta?: number): V3 {
		return v.perturbed(delta)
	}

	static polar(radius: number, phi: number): V3 {
		return new V3(radius * Math.cos(phi), radius * Math.sin(phi), 0)
	}

	static areDisjoint(it: Iterable<V3>): boolean {
		const vSet = new NLA.CustomSet
		for (const v of it) {
			if (!v.equals(vSet.canonicalizeLike(v))) {
				// like value already in set
				return false
			}
		}
		return true
	}


	toAngles(): {theta: number, phi: number} {
		return {
			theta: Math.atan2(this.y, this.x),
			phi: Math.asin(this.z / this.length())
		}
	}

	/**
	 *
	 * @param longitude angle in XY plane
	 * @param latitude "height"/z dir angle
	 */
    static sphere(longitude: number, latitude: number): V3 {
        return new V3(
            Math.cos(latitude) * Math.cos(longitude),
            Math.cos(latitude) * Math.sin(longitude),
            Math.sin(latitude))
    }
}

/**
 * Utility method for creating V3s
 *
 * Example usage:
 *
 *     V(1, 2, 3)
 *     V([1, 2, 3])
 *     V({ x: 1, y: 2, z: 3 })
 *     V(1, 2) * assumes z=0
 *     V([1, 2]) // assumes z=0
 *
 */
function V(a: any, b?: any, c?: any): V3 {
	if (arguments.length == 3) {
		return new V3(parseFloat(a), parseFloat(b), parseFloat(c))
	} else if (arguments.length == 2) {
		return new V3(parseFloat(a), parseFloat(b), 0)
	} else if (arguments.length == 1) {
		if (typeof(a) == 'object') {
			if (a instanceof V3) {
				// immutable, so
				return a
			} else if (a instanceof Array || a instanceof Float32Array || a instanceof Float64Array) {
				if (2 == a.length) {
					return new V3(parseFloat(a[0]), parseFloat(a[1]), 0)
				} else if (3 == a.length) {
					return new V3(parseFloat(a[0]), parseFloat(a[1]), parseFloat(a[2]))
				}
			} else if (('x' in a) && ('y' in a)) {
				if ('z' in a) {
					return new V3(parseFloat(a.x), parseFloat(a.y), parseFloat(a.z))
				} else {
					return new V3(parseFloat(a.x), parseFloat(a.y), 0)
				}
			}
		}
	}
	throw new Error('invalid arguments' + arguments)
}

NLA.registerClass(V3)