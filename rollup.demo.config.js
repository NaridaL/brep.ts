import viewerConfig from './rollup.viewer.config'
Object.assign(viewerConfig, {
	input: 'src/demo.tsx',
	output: { format: 'iife', file: 'dist/demo.js' },
})
export default viewerConfig