import sourcemapsPlugin from "rollup-plugin-sourcemaps"
import typescriptPlugin from "@rollup/plugin-typescript"
import commonjsPlugin from "@rollup/plugin-commonjs"
import nodeResolvePlugin from "@rollup/plugin-node-resolve"
import servePlugin from "rollup-plugin-serve"
import livereloadPlugin from "rollup-plugin-livereload"
import replacePlugin from "@rollup/plugin-replace"
import glsl from "rollup-plugin-glsl"
import * as typescript from "typescript"
import * as fs from "fs"

const pkg = JSON.parse(fs.readFileSync("package.json", "utf-8"))
export default {
  input: "src/viewer.tsx",
  output: {
    format: "iife",
    file: "dist/viewer.js",
    sourcemap: true,
    name: "viewer",
    globals: (moduleName) => {
      const x =
        require(moduleName + "/package.json").umdGlobal ||
        (pkg.umdGlobals && pkg.umdGlobals[moduleName])
      console.log(moduleName, x)
      return x
    },
  },
  external: [],
  plugins: [
    typescriptPlugin({
      typescript,
    }),
    sourcemapsPlugin(),
    nodeResolvePlugin({
      mainFields: ["module", "jsnext:main", "main"],

      // some package.json files have a `browser` field which
      // specifies alternative files to load for people bundling
      // for the browser. If that's you, use this option, otherwise
      // pkg.browser will be ignored
      browser: true, // Default: false

      // not all files you want to resolve are .js files
      extensions: [".js", ".json"], // Default: ['.js']

      // whether to prefer built-in modules (e.g. `fs`, `path`) or
      // local ones with the same names
      preferBuiltins: false, // Default: true

      // Lock the module search in this path (like a chroot). Module defined
      // outside this path will be mark has external
      //   jail: '/my/jail/path', // Default: '/'

      // If true, inspect resolved files to check that they are
      // ES2015 modules
      //   modulesOnly: true, // Default: false

      // Any additional options that should be passed through
      // to node-resolve
      //   customResolveOptions: {
      // 	moduleDirectory: 'js_modules'
      //   }
    }),
    commonjsPlugin({
      // non-CommonJS modules will be ignored, but you can also
      // specifically include/exclude files
      // include: 'node_modules/**',  // Default: undefined
      // exclude: [ 'node_modules/foo/**', 'node_modules/bar/**' ],  // Default: undefined
      // these values can also be regular expressions
      // include: /node_modules/

      // search for files other than .js files (must already
      // be transpiled by a previous plugin!)
      extensions: [".js", ".coffee"], // Default: [ '.js' ]

      // if true then uses of `global` won't be dealt with by this plugin
      ignoreGlobal: false, // Default: false

      // if false then skip sourceMap generation for CommonJS modules
      // sourceMap: false,  // Default: true

      // sometimes you have to leave require statements
      // unconverted. Pass an array containing the IDs
      // or a `id => boolean` function. Only use this
      // option if you know what you're doing!
      // ignore: [ 'conditional-runtime-dependency' ]
    }),
    replacePlugin({
      "process.env.NODE_ENV": JSON.stringify("development"),
    }),
    glsl({
      // By default, everything gets included
      include: "src/**/*.glslx",

      // Undefined by default
      // exclude: ['**/index.html'],

      // Source maps are on by default
      // sourceMap: false
    }),
    process.env.ROLLUP_WATCH && servePlugin({}),
    // {
    // 	contentBase: '.',
    // 	open: true,
    // 	host: 'localhost',
    // 	port: 10002
    // }
    process.env.ROLLUP_WATCH && livereloadPlugin(),
  ].filter((x) => x),
  onwarn: function (warning, warn) {
    // Suppress this error message... there are hundreds of them. Angular team says to ignore it.
    // https://github.com/rollup/rollup/wiki/Troubleshooting#this-is-undefined
    if ("THIS_IS_UNDEFINED" === warning.code) return
    if ("CIRCULAR_DEPENDENCY" === warning.code) return
    warn(warning)
  },
}
