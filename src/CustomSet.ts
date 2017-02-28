namespace NLA {
	export interface Equalable {
		hashCode(): int
		equals(x: any): boolean
		hashCodes?(): int[]
		like?(x: any): boolean
	}

	export class CustomSet<T extends Equalable> implements Set<T> {
		[Symbol.toStringTag]:'Set' = 'Set'

		forEach(callbackfn: (value: T, index: T, set: Set<T>)=>void, thisArg?: any): void {
			for (const value of this.entries()) {
				callbackfn.call(thisArg, value, value, this)
			}
		}
		_map: Map<int, T[]>
		_size: int


		constructor(iterable?: Iterable<T>) {
			this._map = new Map()
			this._size = 0
			if (iterable) {
				this.addAll(iterable)
			}
		}

		add(val: T): this {
			this.add2(val)
			return this
		}

		add2(val: T): boolean {
			// you can't use this.canonicalize here, as there is no way to differentiate if val
			// is new or if val was === the exisitng value (not only .equals)
			const hashCode = val.hashCode(), bucket = this._map.get(hashCode)
			if (bucket) {
				if (bucket.some(x => x.equals(val))) {
					return false
				}
				bucket.push(val)
			} else {
				this._map.set(hashCode, [val])
			}
			this._size++
			return true
		}

		addAll(iterable: Iterable<T>): this {
			for (const val of iterable) {
				this.add(val)
			}
			return this
		}

		canonicalize(val: T): T {
			const hashCode = val.hashCode(), bucket = this._map.get(hashCode)
			if (bucket) {
				const existing = bucket.find(x => x.equals(val))
				if (existing) {
					return existing
				}
				bucket.push(val)
			} else {
				this._map.set(hashCode, [val])
			}
			this._size++
			return val
		}

		has(val: T): boolean {
			const hashCode = val.hashCode(), bucket = this._map.get(hashCode)
			return bucket && bucket.some(x => x.equals(val))
		}

		getLike(val: T) {
			for (const hashCode of val.hashCodes()) {
				const bucket = this._map.get(hashCode)
				const canonVal = bucket && bucket.find(x => x.like(val))
				if (canonVal) return canonVal
			}
		}

		canonicalizeLike(val: T) {
			// if this.getLike(val) is defined, return it, otherwise add val and return val
			return this.getLike(val) || this.canonicalize(val)
		}

		addLike(val) {
			return !this.getLike(val) && this.add(val)
		}

		'delete'(val) {
			const hashCode = val.hashCode(), bucket = this._map.get(hashCode)
			if (bucket) {
				const index = bucket.findIndex(x => x.equals(val))
				if (-1 != index) {
					if (1 == bucket.length) {
						this._map.delete(hashCode)
					} else {
						bucket.splice(index, 1)
					}
					this._size--
					return true
				}
			}
			return false
		}

		deleteLike(val) {
			for (let hashCode of val.hashCodes()) {
				const bucket = this._map.get(hashCode)
				if (bucket) {
					const index = bucket.findIndex(x => x.like(val))
					if (-1 != index) {
						const deleted = bucket[index]
						if (1 == bucket.length) {
							this._map.delete(hashCode)
						} else {
							bucket.splice(index, 1)
						}
						this._size--
						return deleted
					}
				}
			}
		}

		*values(): IterableIterator<T> {
			for (const bucket of this._map.values()) {
				yield* bucket
			}
		}
		*entries(): IterableIterator<[T, T]> {
			for (const bucket of this._map.values()) {
				for (const value of bucket) {
					yield [value, value]
				}
			}
		}

		clear(): void {
			this._map.clear()
			this._size = 0
		}

		get size(): int {
			return this._size
		}

		toString() {
			return '{' + Array.from(this.values()).join(', ') + '}'
		}

		[Symbol.iterator] = CustomSet.prototype.values
		keys = CustomSet.prototype.values
	}

	/**
	 * Java style map.
	 */
	export class CustomMap<K extends {hashCode(): int, equals(x: any): boolean, hashCodes?():int[], like(x: any): boolean}, V> implements Map<K, V> {
		[Symbol.toStringTag]:"Map" = "Map"

		toString() {
			return '{' + Array.from(this.entries2()).map(({key, value}) => key + ':' + value).join(', ') + '}'
		}

		forEach(callbackfn: (value: V, index: K, map: Map<K, V>)=>void, thisArg?: any): void {
			for (const bucket of this._map.values()) {
				for (const {key, value} of bucket) {
					callbackfn.call(thisArg, value, key, this)
				}
			}
		}

		*keys(): IterableIterator<K> {
			for (const bucket of this._map.values()) {
				for (const {key} of bucket) {
					yield key
				}
			}
		}

		*values(): IterableIterator<V> {
			for (const bucket of this._map.values()) {
				for (const {value} of bucket) {
					yield value
				}
			}
		}
		_map: Map<int, {key: K, value: V}[]>
		_size: int

		constructor() {
			this._map = new Map()
			this._size = 0
		}

		[Symbol.iterator]() {
			return this.entries()
		}

		set(key: K, value?: V): this {
			this.set2(key, value)
			return this
		}
		set2(key: K, val: V): boolean {
			const hashCode = key.hashCode(), bucket = this._map.get(hashCode)
			if (bucket) {
				const pairIndex = bucket.findIndex(pair => pair.key.equals(key))
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
		}

		has(key: K): boolean {
			const hashCode = key.hashCode(), bucket = this._map.get(hashCode)
			return bucket && bucket.some(pair => pair.key.equals(key))
		}

		/**
		 * or undefined
		 */
		get(key: K): V {
			const
				hashCode = key.hashCode(),
				bucket = this._map.get(hashCode),
				pair = bucket && bucket.find(pair => pair.key.equals(key))
			return pair && pair.value
		}

		getLike(val) {
			for (let hashCode of val.hashCodes()) {
				const bucket = this._map.get(hashCode)
				const canonVal = bucket && bucket.find(x => x.key.like(val))
				if (canonVal) return canonVal
			}
		}

		setLike(key, val) {
			return !this.getLike(val) && this.set(key, val)
		}

		'delete'(key) {
			const hashCode = key.hashCode(), bucket = this._map.get(hashCode)
			if (bucket) {
				const index = bucket.findIndex(x => x.key.equals(key))
				if (-1 != index) {
					if (1 == bucket.length) {
						this._map.delete(hashCode)
					} else {
						bucket.splice(index, 1)
					}
					this._size--
					return true
				}
			}
			return false
		}

		deleteLike(key) {
			for (const hashCode of key.hashCodes()) {
				const bucket = this._map.get(hashCode)
				if (bucket) {
					const index = bucket.findIndex(x => x.key.like(key))
					if (-1 != index) {
						const deleted = bucket[index]
						if (1 == bucket.length) {
							this._map.delete(hashCode)
						} else {
							bucket.splice(index, 1)
						}
						this._size--
						return deleted
					}
				}
			}
		}

		*entries2(): IterableIterator<{key: K, value: V}> {
			for (const bucket of this._map.values()) {
				yield* bucket
			}
		}

		*entries(): IterableIterator<[K, V]> {
			for (const bucket of this._map.values()) {
				for (const {key, value} of bucket) {
					yield [key, value]
				}
			}
		}

		clear() {
			this._map.clear()
			this._size = 0
		}

		get size() {
			return this._size
		}
	}

	export class Pair<L extends Equalable, R extends Equalable> implements Equalable {

		constructor(public left: L, public right: R) {}

		hashCode() {
			return this.left.hashCode() * 31 + this.right.hashCode()
		}
		equals(other: any) {
			return this == other || Object.getPrototypeOf(other) == Pair.prototype && this.left.equals(other.left) && this.right.equals(other.right)
		}

		toString() {
			return '('+this.left.toString() +', '+ this.right.toString()+')'
		}
		toSource() {
			return 'new NLA.Pair('+this.left.toSource() +', '+ this.right.toSource()+')'
		}
	}
}