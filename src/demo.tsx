import chroma from 'chroma-js'
import classnames from 'classnames'
import { Component, h, render } from 'preact'
import { clamp, DEG, int, round10, TAU, V, V3 } from 'ts3dutils'
import { LightGLContext, Mesh } from 'tsgl'

import { B2, B2T, CustomPlane, Edge, P3 } from './index'
import { BREPGLContext, initMeshes, initNavigationEvents, initShaders, setupCamera } from './viewer'

const fakeB2Mesh = (false as true) && ({} as B2).toMesh()
type B2Mesh = typeof fakeB2Mesh

declare const hljs: any
export type DemoDesc = {
	id: string, f(...args: any[]): B2 | B2[],
	args: DemoArg[],
	container: HTMLDivElement,
	canvas: HTMLCanvasElement,
	gl: BREPGLContext,
	eye: { pos: V3, focus: V3, up: V3, zoomFactor: 8 }
	b2s: B2[]
	meshes: (B2Mesh & {lines?: int[]})[]
}
const demos: DemoDesc[] = []

class Demo extends Component<JSX.HTMLAttributes & {id: string, f: (...args: any[]) => B2 | B2[], args: any[]}, DemoDesc> {
	constructor(demo: DemoDesc) {
		super(demo)
		this.state = demo
		demo.args.forEach(arg => arg.value = arg.def)
	}
	onChange = (e: Event) => {
		const input = e.target as HTMLInputElement, demo = this.state
		const allInputs = Array.from(input.parentElement.parentElement.getElementsByTagName('input'))
		const params = allInputs.map(i => i.value)
		update(demo, demo.args.map(arg => arg.value))
		this.forceUpdate()
		console.log(e)
	}

	componentDidMount() {
		setupDemo(this.state)
	}
	render() {
		const demo = this.state
		const info = demo.b2s && 'faces: ' + demo.b2s.map(b2 => b2.faces.length).join('/')
			+ ' edges: ' + demo.b2s.map(b2 => b2.edgeFaces && b2.edgeFaces.size || '?').join('/')
			+ ' triangles: ' + demo.meshes.map(m => m.TRIANGLES.length / 3).join('/')
		return <div {...this.props} ref={r => demo.container = r as any} class={classnames(this.props.class, 'democontainer')} >
			<canvas ref={r => demo.canvas = r as any}></canvas>
			{demo.args.map(arg =>
				<span class='incont'>
					<InputComponent type='text' data-name={arg.name} value={arg.value} step={arg.step}
						onChange={e => {arg.value = e.target.value; this.onChange(e)}} />
					{arg.name}
				</span>
			)}
			<span class='info'>{info}</span>
			<span class='navinfo'>pan: drag-mmb | rotate: drag-rmb | zoom: scroll</span>
			<a class='sourcelink' href='#'>{true ? 'show source' : 'hide source'}</a>
			{/* <code class={classnames('src', true && 'hide')}>{demo.f.toSource()}</code> */}
		</div>
	}
}

type DemoArg = { step?: number, name: string, def: any, value: any, type: 'int' | 'angle' | 'string' | 'number' }

class InputComponent extends Component<JSX.HTMLAttributes & { step: number, onChange: (e: Event & {target: HTMLInputElement}) => void }, {}> {
	onWheel = (e: WheelEvent) => {
		const target = e.target as HTMLInputElement
		const delta = (e.shiftKey ? 0.1 : 1) * Math.sign(-e.deltaY) * this.props.step
		target.value = '' + round10(parseFloat(target.value) + delta, -6)
		target.dispatchEvent(new Event('change'))
		e.preventDefault()
	}
	render() {
		const { step, ref, ...atts } = this.props
		return <input {...atts}
			class={classnames(this.props.class, step && 'scrollable')}
			onWheel={step && this.onWheel} />
	}
}

export async function demoMain() {
	await B2T.loadFonts()
}
function setupDemo(demo: DemoDesc) {
	//if (demoi !== 0) return
	const { id, f, args, canvas, container } = demo
	// const canvas = new MooEl('canvas', {
	canvas.width = container.clientWidth,
	canvas.height = container.clientHeight,
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
function paintDemo(demo: DemoDesc) {
	const gl = demo.gl
	gl.makeCurrent()
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
	gl.loadIdentity()
	demoPlanes.forEach(plane => gl.drawPlane(plane, chroma(plane.color).gl(), false))
	// gl.drawVectors(g.Vs)
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

function update(demo: DemoDesc, params) {
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
	demo.meshes = demo.b2s && demo.b2s.map(b2 => b2.toMesh().compile())
	paintDemo(demo)
}

window.onload = async e => {
	await B2T.loadFonts()
	render(<Body />, document.body)
}

const Body = () => <div>
	<h1>BREP.TS</h1>
	<p>This library describes volumes using a boundary representation model. Basically like triangle meshes, except faces can be any shape and are not necessarily planar. This allows for fast operations while maintaining a high degree of accuracy.</p>

	<p>Once you have two volumes, you can combine them using boolean operations. For instance:</p>
	<h4>box - sphere</h4>
	<Demo id='boxminussphere' style='width: 100%; height: 500px' f={function (sphereRadius, sphereZ) {
		const box = B2T.box(10, 10, 3)
		const sphere = B2T.sphere(sphereRadius).translate(5, 5, sphereZ)
		const result = box.minus(sphere)
		return [box, sphere, result.translate(12)]
	}} args={[
		{ name: 'sphere radius', type: 'number', fix: val => clamp(val, 0.1, 1000), def: 2, step: 0.5 },
		{ name: 'sphere height', type: 'number', fix: val => clamp(val, 0.1, 1000), def: 2, step: 0.5 }]} />
	<h3>Functionality this library implements</h3>
	<ul>
		<li>Parametric curves: lines, ellipses, parabolas, hyperbolas, quadratic/cubic beziers, intersection curves between parametric and implicit surfaces.</li>
		<li>Surfaces: planes, cylinders, spheres, projected beziers. Functionality includes:
			<ul>
				<li>surface/surface intersections</li>
				<li>testing if two surfaces are coplanar</li>
				<li>testing if a surface contains a point</li>
				<li>testing if a surface contains a curve</li>
			</ul></li>
		<li>Edge: segment of a curve</li>
		<li>Faces: test line intersection/point containment</li>
		<li>BREP volumes: intersection/union/subtraction, conversion to triangle mesh</li>
	</ul>
	<h3>Generator function examples</h3>
	<h4>cylinder(radius: number = 1, height: number = 1, rads: number = TAU)</h4>
	<Demo id='cyl' class='democontainer' style='width: 50%; height: 400px' f={B2T.cylinder} args={[
		{ name: 'radius', type: 'number', fix: val => clamp(val, 0.1, 1000), def: 10, step: 1 },
		{ name: 'height', type: 'number', fix: val => clamp(val, 0.1, 1000), def: 10, step: 1 },
		{ name: 'degrees', type: 'angle', fix: val => clamp(val, 0.1, TAU), def: 360, step: 10 }
	]} />
	<h4>text(text: string, size: number, depth: number = 1, font: opentypejs.Font = defaultFont)</h4>
	<Demo id='text' style='width: 50%; height: 400px' f={B2T.text} args={[
		{ name: 'text', type: 'text', fix: val => val, def: 'foo' },
		{ name: 'size', type: 'number', fix: val => clamp(val, 0.1, 100), def: 10, step: 5 },
		{ name: 'depth', type: 'number', fix: val => clamp(val, 0.1, 100), def: 10, step: 1 }]} />
	<h3>Advanced example: cylinder - extruded rounded edges</h3>
	<Demo id='test1' class='democontainer' style='width: 50%; height: 400px'
		f={(outsideRadius, n, insideRadius, cornerRadius, depth) => {
			const cylinder = B2T.cylinder(outsideRadius, 10)

			// create an n-gon centered on the XY-plane with corners insideRadius away from the origin
			// the ngon is counter-clockwise (CCW) when viewed from +Z
			const ngon = Edge.ngon(n, insideRadius)

			// round the corners of the ngon with radius cornerRadius
			const ngonRounded = Edge.round(ngon, cornerRadius)

			// create a hole punch by extruding ngonRounded in direction +Z
			// Like triangles in opengl, faces face the direction they are CCW when viewed from.
			// Because we extrude in the direction the outline faces, the resulting is "inside-out",
			// with faces facing inwards. As we want to remove the 'punch' volume, this is not a
			// problem. Otherwise we could call volume.flipped() or extrude in -Z and volume.translate(...)
			// into the correct position.
			const punch = B2T.extrudeEdges(ngonRounded, P3.XY, V(0, 0, 10)).translate(0, 0, depth)

			// punch is already inside-out, so use volume.and(...) instead of volume.minus(...)
			const result = cylinder.and(punch)

			return [cylinder, punch, result.translate(2 + outsideRadius * 2)]
		}} args={[
			{ name: 'outside radius', type: 'number', fix: val => clamp(val, 0.1, 10), def: 5, step: 1 },
			{ name: 'corner count', type: 'int', fix: val => clamp(val, 3, 21), def: 3, step: 1 },
			{ name: 'ngon radius', type: 'number', fix: val => clamp(val, 0.1, 12), def: 4, step: 1 },
			{ name: 'corner radius', type: 'number', fix: val => clamp(val, 0.1, 5), def: 0.5, step: 0.1 },
			{ name: 'depth', type: 'number', fix: val => clamp(val, -9, 9), def: 0, step: 1 }]} />
</div>