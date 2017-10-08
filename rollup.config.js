import sourcemaps from 'rollup-plugin-sourcemaps'
import * as fs from 'fs'

const pkg = JSON.parse(fs.readFileSync('package.json'))
export default {
	input: 'out/index.js',
	output: {format: 'umd', file: 'dist/bundle.js'},
	name: pkg.umdGlobal,
	sourcemap: true,
	external: Object.keys(pkg.dependencies),
	// globals: {'javasetmap.ts': '' },
	plugins: [
		sourcemaps()
	],
	onwarn: function (warning) {
		// Suppress this error message... there are hundreds of them. Angular team says to ignore it.
		// https://github.com/rollup/rollup/wiki/Troubleshooting#this-is-undefined
		if (warning.code === 'THIS_IS_UNDEFINED') {
			return
		}
		console.error(warning.message)
	},
}
