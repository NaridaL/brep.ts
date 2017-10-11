import {assert} from 'ts3dutils'
import {getGlobalId} from './index'

export class ClassSerializer {
	CLASS_NAMES = new Map<any, string>()
	NAME_CLASSES = new Map<string, any>()

	constructor() {
		this.addClass('Object', Object)
	}

	addClass(name: string, clazz: any) {
		this.NAME_CLASSES.set(name, clazz)
		this.CLASS_NAMES.set(clazz, name)
	}
	addNamespace(namespaceName: string, namespace: { [symbol: string]: any }) {
		Object.keys(namespace).forEach(symbol => {
			const o = namespace[symbol]
			if ('function' == typeof o && o.name) {
				this.addClass(namespaceName +'.'+ symbol, o)
			}
		})
	}

	serialize(v: any) {
		const gatherList = (v: any) => {
			if (Array.isArray(v)) {
				if (visited.has(v)) {
					if (!listMap.has(v)) {
						listMap.set(v, resultList.length)
						resultList.push(v)
					}
				} else {
					visited.add(v)
					for (let i = 0; i < v.length; i++) {
						gatherList(v[i])
					}
				}
			} else if (undefined !== v && 'object' == typeof v) {
				if (visited.has(v)) {
					if (!listMap.has(v)) {
						listMap.set(v, resultList.length)
						resultList.push(v)
					}
				} else {
					visited.add(v)
					const keys = Object.keys(v)
					for (let i = 0; i < keys.length; i++) {
						gatherList(v[keys[i]])
					}
				}
			}
		}

		const transform = (v: any, allowLinks: boolean, first?: true): any => {
			if ('string' == typeof v || 'number' == typeof v || 'boolean' == typeof v || 'undefined' == typeof v || undefined == v) {
				return v
			}
			let index
			if (allowLinks && !first && undefined !== (index = listMap.get(v))) {
				return {'#REF': index}
			}

			if (Array.isArray(v)) {
				return v.map(x => transform(x, allowLinks))
			}

			//if (mobx && mobx.isObservableArray(v)) {
			//	const result = {'#PROTO': 'ObservableArray'} as any
			//	v.forEach((val, i) => result[i] = transform(val))
			//	return result
			//}

			if ('object' == typeof v) {
				if (v.getConstructorParameters) {
					return {'#CONSTRUCTOR': this.CLASS_NAMES.get(v.constructor), '#ARGS': transform(v.getConstructorParameters(), false)}
				}
				const result = Object.prototype == v.prototype ? {} as any
					: assert(this.CLASS_NAMES.get(v.constructor),
					() => (console.log(v), v.toSource() + v.constructor.name)) && {'#PROTO': this.CLASS_NAMES.get(v.constructor)}
				const keys = Object.keys(v)
				for (let i = 0; i < keys.length; i++) {
					result[keys[i]] = transform(v[keys[i]], allowLinks)
				}
				return result
			}

			throw new Error('?' + typeof v + v.toString())
		}

		const visited = new Set()
		const listMap = new Map()
		let resultList: {}[] = []
		listMap.set(v, 0)
		resultList.push(v)
		gatherList(v)

		resultList = resultList.map(v => transform(v, true, true))
		// console.log(JSON.stringify(resultList))
		const resString = JSON.stringify(resultList)
		return resString
	}

	unserialize(string: string) {
		const fixObjects = (v: any) => {
			if (v && v.constructor === Array) {
				for (let i = 0; i < v.length; i++) {
					v[i] = fixObjects(v[i])
				}
				return v
			} else if ('object' == typeof v && undefined != v) {
				if ('#CONSTRUCTOR' in v) {
					const protoName = v['#CONSTRUCTOR'] as string
					const proto = this.NAME_CLASSES.get(protoName as string)
					assert(proto, protoName + ' Missing ')
					const args: any[] = fixObjects(v['#ARGS'])
					return new proto(...args)
				} else if ('#PROTO' in v) {
					const protoName = v['#PROTO'] as string
					const proto = this.NAME_CLASSES.get(protoName as string)
					assert(proto, protoName + ' Missing ')
					const result = Object.create(proto.prototype)

					const keys = Object.keys(v)
					for (let i = 0; i < keys.length; i++) {
						//if ('name' == keys[i]) console.log(result)
						if ('#PROTO' != keys[i]) {
							result[keys[i]] = fixObjects(v[keys[i]])
							//Object.defineProperty(result, keys[i], {
							//	value: fixObjects(v[keys[i]]),
							//	enumerable: true,
							//	writable: true,
							//	configurable: true
							//})
						}
					}
					Object.defineProperty(result, 'loadID', {value: getGlobalId(), enumerable: false, writable: false})
					return result
				} else {
					return v
				}
			} else {
				return v
			}
		}

		const linkReferences = (v: any) => {
			if (v && v.constructor === Array) {
				for (let i = 0; i < v.length; i++) {
					v[i] = linkReferences(v[i])
				}
				return v
			} else if ('object' == typeof v && undefined != v) {
				if ('#REF' in v) {
					return tree[v['#REF']]
				} else {
					const keys = Object.keys(v)
					for (let i = 0; i < keys.length; i++) {
						v[keys[i]] = linkReferences(v[keys[i]])
					}
					return v
				}
			} else {
				return v
			}
		}

		const tree = JSON.parse(string)
		// console.log(tree)
		fixObjects(tree)
		// console.log(tree)
		linkReferences(tree)
		// console.log(tree)
		return tree[0]
	}

}