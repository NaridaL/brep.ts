import { assert } from "ts3dutils"
import { getGlobalId } from "."

export function doNotSerialize(target: any, key: PropertyKey) {
  const map =
    target.__SERIALIZATION_BLACKLIST || (target.__SERIALIZATION_BLACKLIST = {})
  map[key] = "no"
}

export class ClassSerializer {
  CLASS_NAMES = new Map<any, string>()
  NAME_CLASSES = new Map<string, any>()
  private updater: ((v: any) => void) | undefined

  constructor() {
    this.addClass("Object", Object)
  }

  addClass(name: string, clazz: any) {
    if (this.NAME_CLASSES.has(name)) {
      throw new Error(name)
    }
    this.NAME_CLASSES.set(name, clazz)
    this.CLASS_NAMES.set(clazz, name)
    return this
  }

  addNamespace(namespace: { [symbol: string]: any }, namespaceName?: string) {
    Object.keys(namespace).forEach((symbol) => {
      const o = namespace[symbol]
      if ("function" == typeof o && o.name) {
        this.addClass((namespaceName ? namespaceName + "." : "") + symbol, o)
      }
    })
    return this
  }

  setUpdater(f: (v: any) => void) {
    this.updater = f
    return this
  }

  serialize(v: any) {
    return JSON.stringify(this.serializeObj(v))
  }

  serializeObj(v: any) {
    const path: string[] = []
    const gatherList = (v: any) => {
      //console.log(path.toString())
      if (
        undefined !== v &&
        v.hasOwnProperty("constructor") &&
        this.CLASS_NAMES.has(v.constructor)
      ) {
        // do nothing, this is a class/function prototype
      } else if (Array.isArray(v)) {
        if (visited.has(v)) {
          if (!listMap.has(v)) {
            listMap.set(v, resultList.length)
            resultList.push(v)
          }
        } else {
          visited.add(v)
          for (let i = 0; i < v.length; i++) {
            path.push("" + i)
            gatherList(v[i])
            path.pop()
          }
        }
      } else if (undefined !== v && "object" == typeof v) {
        if (visited.has(v)) {
          if (!listMap.has(v)) {
            listMap.set(v, resultList.length)
            resultList.push(v)
          }
        } else {
          assert(!v.__noxTarget || !visited.has(v.__noxTarget))
          assert(!v.__noxProxy || !visited.has(v.__noxProxy))
          visited.add(v)
          if (!v.getConstructorParameters) {
            for (const key of Object.keys(v).sort()) {
              if (key == "__noxProxy" || key == "__noxTarget") continue
              if (
                !v.__SERIALIZATION_BLACKLIST ||
                !v.__SERIALIZATION_BLACKLIST[key]
              ) {
                path.push(key)
                gatherList(v[key])
                path.pop()
              }
            }
          }
          path.push("proto")
          gatherList(Object.getPrototypeOf(v))
          path.pop()
        }
      }
    }

    const transform = (v: any, allowLinks: boolean, first?: true): any => {
      if (
        "string" == typeof v ||
        "number" == typeof v ||
        "boolean" == typeof v ||
        null === v
      ) {
        return v
      }
      if ("undefined" == typeof v) {
        return { "#REF": -1 }
      }
      if (
        v.hasOwnProperty("constructor") &&
        this.CLASS_NAMES.has(v.constructor)
      ) {
        return { "#REF": this.CLASS_NAMES.get(v.constructor) }
      }
      let index
      if (allowLinks && !first && undefined !== (index = listMap.get(v))) {
        return { "#REF": index }
      }

      if (Array.isArray(v)) {
        return v.map((x) => transform(x, allowLinks))
      }

      //if (mobx && mobx.isObservableArray(v)) {
      //	const result = {'#PROTO': 'ObservableArray'} as any
      //	v.forEach((val, i) => result[i] = transform(val))
      //	return result
      //}

      if ("object" == typeof v) {
        if (v.getConstructorParameters) {
          return {
            "#CONSTRUCTOR": this.CLASS_NAMES.get(v.constructor),
            "#ARGS": transform(v.getConstructorParameters(), false),
          }
        }
        const result: any = {}
        if (Object.prototype !== Object.getPrototypeOf(v)) {
          result["#PROTO"] = transform(Object.getPrototypeOf(v), allowLinks)
        }
        for (const key of Object.keys(v)) {
          if (key == "__noxProxy" || key == "__noxTarget") continue
          if (
            !v.__SERIALIZATION_BLACKLIST ||
            !v.__SERIALIZATION_BLACKLIST[key]
          ) {
            result[key] = transform(v[key], allowLinks)
          }
        }
        return result
      }

      throw new Error("?" + typeof v + v.toString())
    }

    const visited = new Set()
    const listMap = new Map()
    let resultList: {}[] = []
    listMap.set(v, 0)
    resultList.push(v)
    gatherList(v)

    resultList = resultList.map((v) => transform(v, true, true))
    return resultList
  }

  unserialize(string: string) {
    let depth = 0
    const fixObject = (v: any, onReady: (x: any) => void): void => {
      depth++
      if (depth > 100) throw new Error()
      if (v && v.constructor === Array) {
        onReady(v)
        for (let i = 0; i < v.length; i++) {
          fixObject(v[i], (x) => (v[i] = x))
        }
      } else if ("object" == typeof v && undefined != v) {
        if ("#CONSTRUCTOR" in v) {
          const protoName = v["#CONSTRUCTOR"] as string
          const proto = this.NAME_CLASSES.get(protoName as string)
          assert(proto, protoName + " Missing ")
          let args: any[] = undefined!
          fixObject(v["#ARGS"], (x) => (args = x))
          onReady(new proto(...args))
        } else if ("#REF" in v) {
          const ref = v["#REF"]
          if ("string" == typeof ref) {
            onReady(this.NAME_CLASSES.get(ref).prototype)
          } else if ("number" == typeof ref) {
            if (-1 == ref) {
              onReady(undefined)
            } else if (fixedObjects[ref]) {
              onReady(fixedObjects[ref])
            } else {
              fixObject(tree[ref], (x) => onReady((fixedObjects[ref] = x)))
            }
          }
        } else {
          let result: any
          if ("#PROTO" in v) {
            fixObject(v["#PROTO"], (x) => {
              result = Object.create(x)
              onReady(result)
            })
          } else {
            onReady((result = v))
          }

          const keys = Object.keys(v)
          for (let i = 0; i < keys.length; i++) {
            //if ('name' == keys[i]) console.log(result)
            if ("#PROTO" != keys[i]) {
              fixObject(v[keys[i]], (x) => (result[keys[i]] = x))
              //Object.defineProperty(result, keys[i], {
              //	value: fixObjects(v[keys[i]]),
              //	enumerable: true,
              //	writable: true,
              //	configurable: true
              //})
            }
          }
          Object.defineProperty(result, "loadID", {
            value: getGlobalId(),
            enumerable: false,
            writable: false,
          })
          this.updater && this.updater(result)
        }
      } else {
        onReady(v)
      }
      depth--
    }

    // const linkReferences = (v: any) => {
    // 	if (v && v.constructor === Array) {
    // 		for (let i = 0; i < v.length; i++) {
    // 			v[i] = linkReferences(v[i])
    // 		}
    // 		return v
    // 	} else if ('object' == typeof v && undefined != v) {
    // 		if ('#REF' in v) {
    // 			return tree[v['#REF']]
    // 		} else {
    // 			const keys = Object.keys(v)
    // 			for (let i = 0; i < keys.length; i++) {
    // 				v[keys[i]] = linkReferences(v[keys[i]])
    // 			}
    // 			return v
    // 		}
    // 	} else {
    // 		return v
    // 	}
    // }

    const tree = JSON.parse(string)
    // console.log(tree)
    const fixedObjects = new Array(tree.length)
    fixObject({ "#REF": 0 }, () => {})
    // console.log(tree)
    // linkReferences(tree)
    // console.log(tree)
    return fixedObjects[0]
  }
}
