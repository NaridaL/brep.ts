import LightGLContext = GL.LightGLContext
const demos = []
function makeDemo(id: string, f: (...args: any[]) => B2 | B2[], args) {
	demos.push({id, f, args})
}
window.onload = async function () {
	await B2T.loadFonts()
	demos.forEach((demo, demoi) => {
		//if (demoi !== 0) return
		const {id, f, args} = demo
		const container = demo.container = $(id) as HTMLDivElement
		const canvas = new MooEl('canvas', {width: container.clientWidth, height: container.clientHeight}) as HTMLCanvasElement
		$(window).addEvent('resize', () => {
			canvas.width = container.clientWidth
			canvas.height = container.clientHeight
		})
		const gl =demo.gl= GL.create({canvas})
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

		initShaders(gl.shaders = {})
		initMeshes(gl.meshes = {})
		demo.eye = {eyePos: V(10, 10, 100), eyeFocus: V(5, 5, 0), eyeUp: V3.Y, zoomFactor: 8}
		initNavigationEvents(gl, demo.eye, () => paintDemo(demo))
		//initInfoEvents()
		//initPointInfoEvents()
		setupCamera(gl, demo.eye)

		container.adopt(
			canvas,
			args.map(arg => {
				const input = new MooEl('input.canvasinput', {type: 'text', dataName: arg.name, value: arg.def, foo:eye}) as HTMLInputElement
				input.demo = demo
				input.arg = arg
				arg.step && input.addEvent('mousewheel', function (e) {
					const delta = (e.shift ? 10 : 1) * Math.sign(e.wheel) * arg.step
					this.set('value', NLA.round10(parseFloat(this.value) + delta, -6))
					this.fireEvent('change')
					e.preventDefault()
				})
				input.addEvent('change', changeHandler)
				return new MooEl('span.incont').adopt(input).appendText(arg.name)
			}),
			new MooEl('span.info'),
			new MooEl('span.navinfo', {text: 'pan: drag-mmb | rotate: drag-right | zoom: wheel'})
		)

		update(demo, args.map(a => a.def))

	})
}
const meshColorss = [
	chroma.scale(['#ffa621', '#ffd026']).mode('lab').colors(20).map(s => chroma(s).gl()),
	chroma.scale(['#ff297f', '#6636FF']).mode('lab').colors(20).map(s => chroma(s).gl()),
	chroma.scale(['#19ff66', '#1ffff2']).mode('lab').colors(20).map(s => chroma(s).gl())
]
const demoPlanes = [
	new CustomPlane(V3.O, V3.Y, V3.Z, 'planeYZ', 0xff0000, -100, 100, -100, 100),
	new CustomPlane(V3.O, V3.X, V3.Z, 'planeZX', 0x00ff00, -100, 100, -100, 100),
	new CustomPlane(V3.O, V3.X, V3.Y, 'planeXY', 0x0000ff, -100, 100, -100, 100),
	//	sketchPlane
]
function paintDemo(demo) {
	const gl = demo.gl
	gl.makeCurrent()
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
	gl.loadIdentity()
	demoPlanes.forEach(plane => drawPlane(plane, plane.color, gl))
	drawVectors(gl)
	if (!demo.meshes) return
	//gl.scale(100, 100, 100)

	for (let i = 0; i < demo.meshes.length; i++) {
		const mesh = demo.meshes[i], b2 = demo.b2s[i]
		gl.pushMatrix()
		//gl.translate(30, 0, 0)
		gl.projectionMatrix.m[11] -= 1 / (1 << 22) // prevent Z-fighting
		mesh.lines && gl.shaders.singleColor.uniforms({ color: chroma('#bfbfbf').gl() }).draw(mesh, 'LINES')
		gl.projectionMatrix.m[11] += 1 / (1 << 22)

		let faceIndex = b2.faces.length
		while (faceIndex--) {

			const face = b2.faces[faceIndex]
			const faceTriangleIndexes = mesh.faceIndexes.get(face)
			gl.shaders.lighting.uniforms({
				color: hovering == face ? chroma('purple').gl() : meshColorss[i % meshColorss.length][faceIndex % randomColors.length]
			}).draw(mesh, 'TRIANGLES', faceTriangleIndexes.start, faceTriangleIndexes.count)
			//shaders.singleColor.uniforms({
			//color: rgbToVec4(0x0000ff)
			//}).draw(brepMesh, 'LINES')
		}

		gl.popMatrix()
	}
}
function update(demo, params) {
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
	try {
		const b2s = demo.f.apply(undefined, fixedParams)
		demo.b2s = b2s instanceof Array ? b2s : [b2s]
		demo.b2s.forEach(b2 => b2.buildAdjacencies())
	} catch (e) {
		console.log(e.message)
	}
	demo.gl.makeCurrent()
	demo.meshes = demo.b2s && demo.b2s.map(b2 => b2.toMesh())
	const info = demo.b2s && 'faces: ' + demo.b2s.map(b2 => b2.faces.length).join('/')
		+ ' edges: ' + demo.b2s.map(b2 => b2.edgeFaces && b2.edgeFaces.size || '?').join('/')
		+ ' triangles: ' + demo.meshes.map(m => m.triangles.length / 3).join('/')
	demo.container.getElement('.info').set('text', info)
	paintDemo(demo)
}
function changeHandler(e) {
	const input = this
	const allInputs = input.getParent().getParent().getElements('input')
	const params = allInputs.map(i => i.value)
	update(input.demo, params)
	console.log(e)
}