import viewerConfig from "./rollup.viewer.config"
import uglify from "rollup-plugin-uglify"
import Visualizer from "rollup-plugin-visualizer"
import { minify } from "uglify-es"
import merge from "deepmerge"

const config = merge(viewerConfig, {
  input: "src/demo.tsx",
  output: {
    format: "iife",
    file: "dist/demo.js",
    name: "demo",
    globals: {
      "highlight.js": "hljs",
      "svg-pathdata": "svgpathdata",
    },
  },
  external: ["highlight.js"],
  // plugins: [
  // 	uglify({}, minify),
  // ],
})
export default config
