import viewerConfig from './rollup.viewer.config'
import uglify from 'rollup-plugin-uglify';
import Visualizer from 'rollup-plugin-visualizer'
import { minify } from 'uglify-es'
const config = Object.assign(viewerConfig, {
	input: 'src/demo.tsx',
	output: {
		format: 'iife',
		file: 'dist/demo.js',
	},
})
config.plugins.push(
	uglify({}, minify),
	Visualizer({
		sourcemap: true
	}),
)
export default viewerConfig