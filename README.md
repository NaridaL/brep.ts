<!--- header generated automatically, don't edit --->
[![Travis](https://img.shields.io/travis/NaridaL/brep.ts.svg?style=flat-square)](https://travis-ci.org/NaridaL/brep.ts)
[![npm](https://img.shields.io/npm/v/brep.ts.svg?style=flat-square)](https://www.npmjs.com/package/brep.ts)
[![David](https://img.shields.io/david/expressjs/express.svg?style=flat-square)](https://david-dm.org/NaridaL/brep.ts)

# brep.ts
Boundary representation volume modeling in TypeScript.

## Installation
NPM:  `npm install brep.ts --save`

<!--- CONTENT-START --->
Work in progress. Allows modelling volumes a boundary representations and intersecting them. See the [interactive demo](https://naridal.github.io/brep.ts/demo.html)

## Example Usage
```ts
import {B2T} from 'brep.ts'
const sphere = B2T.sphere()
const box = B2T.box()
const result = sphere.minus(box)
const resultMesh = result.toMesh()

```
<!--- CONTENT-END --->
<!--- footer generated automatically, don't edit --->
## License
[MIT](./LICENSE)
