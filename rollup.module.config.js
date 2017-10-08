import sourcemaps from 'rollup-plugin-sourcemaps';
import * as fs from 'fs';

export default {
	input: 'out/index.js',
	output: {format: 'es', file: 'dist/bundle.module.js'},
	sourcemap: true,
	plugins: [
		sourcemaps()
	],
	external: Object.keys(JSON.parse(fs.readFileSync('package.json')).dependencies),
	onwarn: function (warning) {
		// Suppress this error message... there are hundreds of them. Angular team says to ignore it.
		// https://github.com/rollup/rollup/wiki/Troubleshooting#this-is-undefined
		if (warning.code === 'THIS_IS_UNDEFINED') {
			return
		}
		console.error(warning.message)
	},
};