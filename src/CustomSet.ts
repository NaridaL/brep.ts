namespace NLA {
	export class CustomSet<T extends {hashCode(): int, equals(x: any): boolean, hashCodes?():int[], like(x: any): boolean}> {
		_map: Map<number, T[]>
		_size: int


		constructor(iterable?: Iterable<T>) {
			this._map = new Map()
			this._size = 0
			if (iterable) {
				this.addAll(iterable)
			}
		}

		add(val: T): boolean {
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

		addAll(iterable: Iterable<T>): void {
			for (var val of iterable) {
				this.add(val)
			}
		}

		canonicalize(val: T): T {
			var hashCode = val.hashCode(), bucket = this._map.get(hashCode)
			if (bucket) {
				var existing = bucket.find(x => x.equals(val))
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
			var hashCode = val.hashCode(), bucket = this._map.get(hashCode)
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

		delete(val) {
			var hashCode = val.hashCode(), bucket = this._map.get(hashCode)
			if (bucket) {
				var index = bucket.findIndex(x => x.equals(val))
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
			for (var hashCode of val.hashCodes()) {
				var bucket = this._map.get(hashCode)
				if (bucket) {
					var index = bucket.findIndex(x => x.like(val))
					if (-1 != index) {
						var deleted = bucket[index]
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

		*entries(): IterableIterator<T> {
			for (const bucket of this._map) {
				yield* bucket
			}
		}

		clear(): void {
			this._map.clear()
			this._size = 0
		}

		size(): int {
			return this._size
		}

		[Symbol.iterator] = CustomSet.prototype.entries
		values = CustomSet.prototype.entries;
		keys = CustomSet.prototype.entries;
	}

	class CustomMap {
		constructor() {

		}

		set(key, val) {
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
		}

		has(key) {
			var hashCode = key.hashCode(), bucket = this._map.get(hashCode)
			return bucket && bucket.some(x => x.key.equals(key))
		}

		getLike(val) {
			for (var hashCode of val.hashCodes()) {
				var bucket = this._map.get(hashCode)
				var canonVal = bucket && bucket.find(x => x.key.like(val))
				if (canonVal) return canonVal
			}
		}

		setLike(key, val) {
			return !this.getLike(val) && this.set(key, val)
		}

		delete(key) {
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
		}

		deleteLike(key) {
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
		}

		*entries() {
			for (var bucket of this._map) {
				yield* bucket
			}
		}

		clear() {
			this._map.clear()
			this._size = 0
		}

		size() {
			return this._size
		}
	}
}