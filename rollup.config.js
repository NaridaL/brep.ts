import sourcemaps from 'rollup-plugin-sourcemaps'
import * as fs from 'fs'
import livereload from 'rollup-plugin-livereload'
import serve from 'rollup-plugin-serve'

const pkg = JSON.parse(fs.readFileSync('package.json'))
export default {
	input: 'out/index.js',
	output: [
		// {
		// 	format: 'cjs',
		// 	file: 'dist/bundle.js',
		// 	name: pkg.umdGlobal,
		// 	sourcemap: true,
		// 	globals: moduleName => require(moduleName + '/package.json').umdGlobal || pkg.umdGlobals && pkg.umdGlobals[moduleName],
		// },
		{
			format: 'es',
			sourcemap: true,
			file: 'dist/bundle.module.js'
		},
	],
	external: Object.keys(pkg.dependencies),
	// plugins: [
	// 	sourcemaps(),
    // process.env.ROLLUP_WATCH && serve(
    //   // {
    //   // 	contentBase: '.',
    //   // 	open: true,
    //   // 	host: 'localhost',
    //   // 	port: 10002
    //   // }
    // ),
    // process.env.ROLLUP_WATCH && livereload(),
	// ].filter(x => x), // remove undefined
	onwarn: function (warning, warn) {
		// Suppress this error message... there are hundreds of them. Angular team says to ignore it.
		// https://github.com/rollup/rollup/wiki/Troubleshooting#this-is-undefined
		if ('THIS_IS_UNDEFINED' === warning.code) return
		if ('CIRCULAR_DEPENDENCY' === warning.code) return

		warn(warning)
	},
}
