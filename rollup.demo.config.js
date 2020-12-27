import viewerConfig from "./rollup.viewer.config"
import merge from "deepmerge"
import servePlugin from "rollup-plugin-serve"

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
  plugins: [
    process.env.ROLLUP_WATCH &&
      servePlugin({
        port: 9995,
        open: true,
        openPage: "/demo.html",
      }),
  ],
})
export default config
