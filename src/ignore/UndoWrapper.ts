class UndoArray<T> implements Array<T> {

	constructor(baseArray: T[]) {
	}

	length: number

	[Symbol.iterator](): IterableIterator<T> {
		return undefined
	}

	[Symbol.unscopables](): {copyWithin: boolean; entries: boolean; fill: boolean; find: boolean; findIndex: boolean; keys: boolean; values: boolean} {
		return undefined
	}

	entries(): IterableIterator<any> {
		return undefined
	}

	keys(): IterableIterator<number> {
		return undefined
	}

	values(): IterableIterator<T> {
		return undefined
	}

	find(predicate: (value: T, index: number, obj) => boolean, thisArg?: any): T {
		return undefined
	}

	findIndex(predicate: (value: T)=>boolean, thisArg?: any): number {
		return undefined
	}

	fill(value: T, start?: number, end?: number): T[] {
		return undefined
	}

	copyWithin(target: number, start: number, end?: number): T[] {
		return undefined
	}

	withMax(f: (el: T)=>number): T {
		return undefined
	}

	mapFilter<U>(f: (el: T)=>U): U[] {
		return undefined
	}

	isEmpty(): boolean {
		return undefined
	}

	remove(T): boolean {
		return undefined
	}

	unique(): T[] {
		return undefined
	}

	last(): T {
		return undefined
	}

	push(...items: T[]): number {
		Array.prototype.push.apply(this.baseArray, items.map(UndoObject.for))
	}

	pop(): T {
		return undefined
	}

	concat<U extends T[]>(items: U): T[] {
		return undefined
	}

	join(separator?: string): string {
		return undefined
	}

	reverse(): T[] {
		return undefined
	}

	shift(): T {
		return undefined
	}

	slice(start?: number, end?: number): T[] {
		return undefined
	}

	sort(compareFn?: (a: T, b: T)=>number): T[] {
		return undefined
	}

	unshift(items: T): number {
		return undefined
	}

	indexOf(searchElement: T, fromIndex?: number): number {
		return undefined
	}

	lastIndexOf(searchElement: T, fromIndex?: number): number {
		return undefined
	}

	every(callbackfn: (value: T, index: number, array: T[])=>boolean, thisArg?: any): boolean {
		return undefined
	}

	some(callbackfn: (value: T, index: number, array: T[])=>boolean, thisArg?: any): boolean {
		return undefined
	}

	forEach(callbackfn: (value: T, index: number, array: T[])=>void, thisArg?: any): void {
	}

	map<U>(callbackfn: (value: T, index: number, array: T[])=>U, thisArg?: any): U[] {
		return undefined
	}

	filter(callbackfn: (value: T, index: number, array: T[])=>boolean, thisArg?: any): T[] {
		return undefined
	}

	reduce(callbackfn: (previousValue: T, currentValue: T, currentIndex: number, array: T[])=>T, initialValue?: T): T {
		return undefined
	}

	reduce<U>(callbackfn: (previousValue: U, currentValue: T, currentIndex: number, array: T[])=>U, initialValue: U): U {
		return undefined
	}

	reduceRight(callbackfn: (previousValue: T, currentValue: T, currentIndex: number, array: T[])=>T, initialValue?: T): T {
		return undefined
	}

	reduceRight<U>(callbackfn: (previousValue: U, currentValue: T, currentIndex: number, array: T[])=>U, initialValue: U): U {
		return undefined
	}

}

const UNDO_STACK = []
//class UndoObject {
//	static from(o: any) {
//		let keys = Object.getOwnPropertyNames(o)
//	const result = Object.create(null, )
//		for (const i = 0; i < keys.length; i++) {
//			let key = keys[i]
//			Object.definePropertyconstsult, key, {
//				get: function () {
//					return o[key]
//				},
//				set: function (newValue) {
//					if (!(newValue instanceof UndoObject)) newValue = UndoObject.from(newValue)
//					UNDO_STACK.push({o: o, key: key, oldValue: o[key], newValue: newValue})
//					o[key] = UndoObject.from(newValue)
//				}
//			})}}
//		}
//		return result
//	}
//}