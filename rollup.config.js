import * as fs from "fs"
import typescriptPlugin from "@rollup/plugin-typescript"
import typescript from "typescript"
import { terser as terserPlugin } from "rollup-plugin-terser"

const pkg = JSON.parse(fs.readFileSync("package.json", "utf-8"))
export default {
  input: "src/index.ts",
  output: [
    ["es", false],
    ["es", true],
    ["umd", false],
    ["umd", true],
  ].map(([format, compress]) => ({
    format: format,
    entryFileNames: "[name].[format]" + (compress ? ".min" : "") + ".js",
    sourcemap: true,
    dir: "lib",
    name: pkg.umdGlobal,
    exports: "named",
    globals: pkg.umdGlobals,
    plugins: compress
      ? [
          terserPlugin({
            compress: {
              passes: 3,
              unsafe: true,
              ecma: 7,
            },
            toplevel: true,
            mangle: {
              properties: { regex: /^_/ },
            },
          }),
        ]
      : [],
  })),
  external: Object.keys(pkg.dependencies),
  plugins: [typescriptPlugin({ typescript })],
  onwarn: function (warning, warn) {
    // Suppress this error message... there are hundreds of them. Angular team says to ignore it.
    // https://github.com/rollup/rollup/wiki/Troubleshooting#this-is-undefined
    if ("THIS_IS_UNDEFINED" === warning.code) return
    if ("CIRCULAR_DEPENDENCY" === warning.code) {
      const m = warning.message.match(/^Circular dependency: (.*) -> .* -> .*$/)
      if (m) {
        const start = m[1]
        if (start.match(/out[/\\]index.js|src[/\\]index.ts/)) {
          // this is a loop of length three starting at the index file: don't warn
          return
        }
      }
    }

    warn(warning)
  },
}
