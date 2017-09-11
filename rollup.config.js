import sourcemaps from 'rollup-plugin-sourcemaps';
export default {
	input: 'out/index.js',
	output: {format: 'umd', file: 'dist/bundle.js'},
	name: 'nla',
	sourcemap: true,
	external: ['chai'],
	plugins: [
		sourcemaps()
	]
};