import sourcemaps from 'rollup-plugin-sourcemaps';
export default {
	input: 'out/index.js',
	output: {format: 'es', file: 'dist/bundle.module.js'},
	sourcemap: true,
	plugins: [
		sourcemaps()
	]
};