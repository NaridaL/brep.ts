import chroma from 'chroma-js'
import classnames from 'classnames'
import { Component, h, render } from 'preact'
import { DEG, round10, V, V3 } from 'ts3dutils'
import { LightGLContext, Mesh } from 'tsgl'

import { B2, B2T, CustomPlane } from './index'
import { BREPGLContext, initMeshes, initNavigationEvents, initShaders, setupCamera } from './viewer'

declare const hljs: any
export type Demo = {
	id: string, f(...args: any[]): B2 | B2[],
	args: DemoArg[],
	container: HTMLDivElement,
	canvas: HTMLCanvasElement,
	gl: BREPGLContext,
	eye: { pos: V3, focus: V3, up: V3, zoomFactor: 8 }
	b2s: B2[]
	meshes: Mesh[]
}
const demos: Demo[] = []

export function makeDemo(id: string, f: (...args: any[]) => B2 | B2[], args: any[]) {
	render(h(DemoView, { id, f, args }), document.currentScript)
}

class DemoView extends Component<Demo, Demo> {
	constructor(demo: Demo) {
		super(demo)
		this.state = demo
	}
	onChange(e: Event) {
		const input = e.target as HTMLInputElement, demo = this.state
		const allInputs = Array.from(input.parentElement.parentElement.getElementsByTagName('input'))
		const params = allInputs.map(i => i.value)
		update(demo, params)
		this.forceUpdate()
		console.log(e)
	}

	onComponentHasMounted() {
		setupDemo(this.state)
	}
	render() {
		const demo = this.state
		const info = demo.b2s && 'faces: ' + demo.b2s.map(b2 => b2.faces.length).join('/')
			+ ' edges: ' + demo.b2s.map(b2 => b2.edgeFaces && b2.edgeFaces.size || '?').join('/')
			+ ' triangles: ' + demo.meshes.map(m => m.TRIANGLES.length / 3).join('/')
		return <div ref={r => demo.container = r as any}>
			<canvas ref={r => demo.canvas = r as any}></canvas>
			{demo.args.map(arg =>
				<span class='incont'>
					<input type='text' data-name={arg.name} value={arg.def} onChange={this.onChange} />
					{arg.name}
				</span>
			)}
			<span class='info'>{info}</span>
			<span class='navinfo'>pan: drag-mmb | rotate: drag-rmb | zoom: scroll</span>
			<a class='sourcelink' href='#'>{true ? 'show source' : 'hide source'}</a>
			<code class={classnames('src', true && 'hide')}>{demo.f.toSource()}</code>
		</div>
	}
}

type DemoArg = { step?: number, name: string, def, type: 'int' | 'angle' | 'string' | 'number' }

class InputComponent extends Component<{ ref: any, arg: DemoArg }, {}> {
	onWheel = (e: WheelEvent) => {
		const target = e.target as HTMLInputElement
		const delta = (e.shiftKey ? 0.1 : 1) * Math.sign(e.deltaY) * this.props.arg.step
		target.value = '' + round10(parseFloat(target.value) + delta, -6)
		// target.fire('change')
		e.preventDefault()
	}
	render() {
		const { arg, ref, ...atts } = this.props
		return <input {...atts} class={classnames(this.props.class, arg.step && 'scrollable')} onWheel={arg.step && this.onWheel} />
	}
}

export async function demoMain() {
	await B2T.loadFonts()
}
function setupDemo(demo: Demo) {
	//if (demoi !== 0) return
	const { id, f, args, canvas, container } = demo
	// const canvas = new MooEl('canvas', {
	// 	width: container.clientWidth,
	// 	height: container.clientHeight,
	// }) as HTMLCanvasElement
	window.addEventListener('resize', e => {
		canvas.width = container.clientWidth
		canvas.height = container.clientHeight
		gl.viewport(0, 0, canvas.width, canvas.height)
		setupCamera(demo.eye, gl)
	})
	const gl = demo.gl = BREPGLContext.create(LightGLContext.create({ canvas }))
	gl.clearColor(1.0, 1.0, 1.0, 0.0)
	gl.enable(gl.BLEND)
	gl.enable(gl.DEPTH_TEST)
	gl.enable(gl.CULL_FACE)
	canvas.oncontextmenu = () => false
	gl.depthFunc(gl.LEQUAL)
	gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA) // TODO ?!

	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
	gl.loadIdentity()
	gl.scale(10, 10, 10)

	gl.loadIdentity()

	gl.shaders = initShaders(gl)
	initMeshes(gl.meshes = {}, gl)
	demo.eye = { pos: V(10, 10, 100), focus: V(5, 5, 0), up: V3.Y, zoomFactor: 8 }
	initNavigationEvents(gl, demo.eye, () => paintDemo(demo))
	//initInfoEvents()
	//initPointInfoEvents()
	setupCamera(demo.eye, gl)

	// container.adopt(
		// demo.srcLink = new MooEl('a.sourcelink', {text: 'show source', href: '#'})
		// 	.addEvent('click', e => {
		// 		const showing = demo.srcLink.get('text') == 'hide source'
		// 		demo.srcContainer.setStyle('display', showing ? 'none' : 'block')
		// 		demo.srcLink.set('text', showing ? 'show source' : 'hide source')
		// 		return false
		// 	}),
	// )
	// hljs.highlightBlock(demo.srcContainer)
	update(demo, args.map(a => a.def))

}

const meshColorss = [
	chroma.scale(['#ffa621', '#ffd026']).mode('lab').colors(20, 'gl'),
	chroma.scale(['#ff297f', '#6636FF']).mode('lab').colors(20, 'gl'),
	chroma.scale(['#19ff66', '#1ffff2']).mode('lab').colors(20, 'gl'),
]
const demoPlanes = [
	new CustomPlane(V3.O, V3.Y, V3.Z, 'planeYZ', 0xff0000, -5, 5, -5, 5),
	new CustomPlane(V3.O, V3.X, V3.Z, 'planeZX', 0x00ff00, -5, 5, -5, 5),
	new CustomPlane(V3.O, V3.X, V3.Y, 'planeXY', 0x0000ff, -5, 5, -5, 5),
	//	sketchPlane
]
const hovering = undefined
function paintDemo(demo) {
	const gl = demo.gl
	gl.makeCurrent()
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
	gl.loadIdentity()
	demoPlanes.forEach(plane => gl.drawPlane(plane, plane.color, gl))
	gl.drawVectors(gl)
	if (!demo.meshes) return
	//viewerGL.scale(100, 100, 100)

	for (let i = 0; i < demo.meshes.length; i++) {
		const mesh = demo.meshes[i], b2 = demo.b2s[i]
		gl.pushMatrix()
		//viewerGL.translate(30, 0, 0)
		gl.projectionMatrix.m[11] -= 1 / (1 << 22) // prevent Z-fighting
		mesh.lines && gl.shaders.singleColor.uniforms({ color: chroma('#bfbfbf').gl() }).draw(mesh, gl.LINES)
		gl.projectionMatrix.m[11] += 1 / (1 << 22)

		let faceIndex = b2.faces.length
		while (faceIndex--) {

			const face = b2.faces[faceIndex]
			const faceTriangleIndexes = mesh.faceIndexes.get(face)
			gl.shaders.lighting.uniforms({
				color: hovering == face ? chroma('purple').gl() : meshColorss.emod(i).emod(faceIndex),
			}).draw(mesh, gl.TRIANGLES, faceTriangleIndexes.start, faceTriangleIndexes.count)
			//shaders.singleColor.uniforms({
			//color: hexIntToGLColor(0x0000ff)
			//}).draw(brepMesh, 'LINES')
		}

		gl.popMatrix()
	}
}

function update(demo: Demo, params) {
	const fixedParams = params.map((p, i) => {
		switch (demo.args[i].type) {
		case 'number':
		case 'int':
			return parseFloat(p)
		case 'angle':
			return parseFloat(p) * DEG
		}
		return p
	})
	//try {
	const b2s = demo.f.apply(undefined, fixedParams)
	demo.b2s = b2s instanceof Array ? b2s : [b2s]
	demo.b2s.forEach(b2 => b2.buildAdjacencies())
	//} catch (e) {
	//	console.log(e.message)
	//}
	demo.gl.makeCurrent()
	demo.meshes = demo.b2s && demo.b2s.map(b2 => b2.toMesh())
	paintDemo(demo)
}
