function parseGetParams(str) {
    const result = {}
    str
        .split('&')
        .forEach(function (item) {
        	const splitIndex = item.indexOf('=')
            result[item.substr(0, splitIndex)] = decodeURI(item.substr(splitIndex + 1))
        })
    return result
}
let a, b, c, d, edges = [], hovering
function initB2() {
    dMesh = new GL.Mesh()
    eye.eyePos = V(1, 2, 101)
	eye.eyeFocus = V(0, 1, 0)
	eye.eyeUp = V(0, 1, 0)
    eye.zoomFactor = 8

    const gets = parseGetParams(window.location.search.substr(1) || window.location.hash.substr(1))
    'abcd'.split('').forEach(k => gets[k] && eval(k + '=' + gets[k] + ';' + k + 'Mesh = ' + k + '.toMesh()'))

    //cMesh && cMesh.computeWireframeFromFlatTriangles() && cMesh.compile()
    if (gets['points']) {
        console.log('drPs from GET')
        drPs = eval(gets['points'])
    }
	if (gets['vectors']) {
		console.log('vectors from GET')
		drVs.pushAll(eval(gets['vectors']))
	}

    if (gets['edges']) {
        console.log('edges from GET')
        dMesh = new GL.Mesh({triangles: false})
        edges = eval(gets['edges'])
        edges && dMesh.addVertexBuffer('curve1', 'curve1')
        edges.forEach(edge => {
            const points = edge.points()
            for (let i = 0; i < points.length - 1; i++) {
                dMesh.curve1.push(points[i], points[i + 1])
            }
        })
    }
    if (gets['mesh']) {
        console.log('mesh from GET')
        sMesh = eval(gets['mesh'])
        sMesh.computeWireframeFromFlatTriangles()
        //sMesh.computeNormalLines()
        sMesh.compile()
    }
    if (gets['meshes']) {
        console.log('meshes from GET')
        b2meshes = eval(gets['meshes']) || []
        b2meshes.forEach(m => m.computeWireframeFromFlatTriangles())
	    b2meshes.forEach(m => m.computeNormalLines(0.5))
	    b2meshes.forEach(m => m.compile())
        console.log('meshes from GET', b2meshes)
    }

    dMesh.compile()
}





const randomColors = chroma.scale(['#ff297f', '#6636FF']).mode('lab').colors(20).map(s => chroma(s).gl())
let aMesh: GL.Mesh & {faceIndexes: Map<Face, {start: int, count: int}>},
	bMesh: GL.Mesh & {faceIndexes: Map<Face, {start: int, count: int}>},
	cMesh: GL.Mesh & {faceIndexes: Map<Face, {start: int, count: int}>},
	dMesh: GL.Mesh & {faceIndexes: Map<Face, {start: int, count: int}>},
	sMesh: GL.Mesh,
	b2meshes: GL.Mesh[] = []
function viewerPaint() {
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
    gl.loadIdentity()

    drawVectors()

    //gl.scale(100, 100, 100)

    if (aMesh) {
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        aMesh.lines && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.PP_STROKE) }).draw(aMesh, 'LINES')
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
        shaders.lighting.uniforms({ color: rgbToVec4(COLORS.PP_FILL),
            camPos: eye.eyePos }).draw(aMesh)
    }

    if (bMesh) {
        gl.pushMatrix()
        //gl.translate(15, 0, 0)
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        bMesh.lines && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.PP_STROKE) }).draw(bMesh, 'LINES')
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
        shaders.lighting.uniforms({ color: rgbToVec4(COLORS.RD_FILL),
            camPos: eye.eyePos }).draw(bMesh)
        bMesh.edgeTangents && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.TS_STROKE) })
            .drawBuffers({LGL_Vertex: bMesh.vertexBuffers.edgeTangents}, null, gl.LINES)
        bMesh.edgeTangents2 && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.RD_STROKE) })
            .drawBuffers({LGL_Vertex: bMesh.vertexBuffers.edgeTangents2}, null, gl.LINES)
        gl.popMatrix()
    }
    if (cMesh) {
        gl.pushMatrix()
        //gl.translate(30, 0, 0)
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        cMesh.lines && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.TS_STROKE) }).draw(cMesh, 'LINES')
        gl.projectionMatrix.m[11] += 1 / (1 << 20)

        let faceIndex = c.faces.length
        while (faceIndex--) {

            const face = c.faces[faceIndex]
            const faceTriangleIndexes = cMesh.faceIndexes.get(face)
            shaders.lighting.uniforms({
                color: hovering == face ? chroma('purple').gl() : randomColors[faceIndex % randomColors.length]
            }).draw(cMesh, 'TRIANGLES', faceTriangleIndexes.start, faceTriangleIndexes.count)
            //shaders.singleColor.uniforms({
            //color: rgbToVec4(0x0000ff)
            //}).draw(brepMesh, 'LINES')
        }

        gl.popMatrix()
    }
    if (dMesh && dMesh.hasBeenCompiled) {
        gl.pushMatrix()
        //gl.scale(10, 10, 10)
        gl.projectionMatrix.m[11] -= 1 / (1 << 21) // prevent Z-fighting
        dMesh.lines && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.RD_STROKE) }).draw(dMesh, 'LINES')
        gl.projectionMatrix.m[11] += 1 / (1 << 21)
        dMesh.triangles && shaders.lighting.uniforms({ color: rgbToVec4(0xffFF00),
            camPos: eye.eyePos }).draw(dMesh)

        dMesh.curve1 && shaders.singleColor.uniforms({ color: rgbToVec4(0xff00000) })
            .drawBuffers({LGL_Vertex: dMesh.vertexBuffers.curve1}, null, gl.LINES)
        dMesh.curve2 && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.RD_STROKE) })
            .drawBuffers({LGL_Vertex: dMesh.vertexBuffers.curve2}, null, gl.LINES)

        dMesh.curve3 && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.PP_STROKE) })
            .drawBuffers({LGL_Vertex: dMesh.vertexBuffers.curve3}, null, gl.LINES)
        dMesh.curve4 && shaders.singleColor.uniforms({ color: rgbToVec4(0x00ff00) })
            .drawBuffers({LGL_Vertex: dMesh.vertexBuffers.curve4}, null, gl.LINES)
        gl.popMatrix()
    }
    if (sMesh && sMesh.hasBeenCompiled) {
        gl.pushMatrix()
        //gl.scale(10, 10, 10)
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        sMesh.lines && shaders.singleColor.uniforms({ color: rgbToVec4(0xFF6600) }).draw(sMesh, 'LINES')
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
        sMesh.triangles && shaders.lighting.uniforms({ color: rgbToVec4(0xffFF00),
            camPos: eye.eyePos }).draw(sMesh)

        sMesh.curve1 && shaders.singleColor.uniforms({ color: rgbToVec4(0xff00000) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve1}, null, gl.LINES)
        sMesh.curve2 && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.RD_STROKE) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve2}, null, gl.LINES)

        sMesh.curve3 && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.PP_STROKE) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve3}, null, gl.LINES)
        sMesh.curve4 && shaders.singleColor.uniforms({ color: rgbToVec4(0x00ff00) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve4}, null, gl.LINES)
        gl.popMatrix()
    }
    for (const sMesh of b2meshes) {
        gl.pushMatrix()
        //gl.scale(10, 10, 10)
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        sMesh.lines && shaders.singleColor.uniforms({ color: rgbToVec4(0xFF6600) }).draw(sMesh, 'LINES')
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
        sMesh.triangles && shaders.lighting.uniforms({ color: rgbToVec4(0xffFF00),
            camPos: eye.eyePos }).draw(sMesh)

        sMesh.curve1 && shaders.singleColor.uniforms({ color: rgbToVec4(0xff00000) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve1}, null, gl.LINES)
        sMesh.curve2 && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.RD_STROKE) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve2}, null, gl.LINES)

        sMesh.curve3 && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.PP_STROKE) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve3}, null, gl.LINES)
        sMesh.curve4 && shaders.singleColor.uniforms({ color: rgbToVec4(0x00ff00) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve4}, null, gl.LINES)
        gl.popMatrix()
    }

    if (hovering instanceof Edge) {
	    gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
	    drawEdge(hovering, 0x000000, 0.1)
	    gl.projectionMatrix.m[11] += 1 / (1 << 20)
    }

	drawEdge(Edge.forCurveAndTs(SemiEllipseCurve.UNIT.scale(100)), 0x000000, 0.1)

    //drPs.forEach(v => drawPoint(v, undefined, 0.3))
    drawPoints(0.05)
	b2planes.forEach(plane => drawPlane(plane, plane.color))
}

















function initInfoEvents() {
	gl.onmousemove.push(function (e) {
		const mouseLine = getMouseLine({x: e.clientX, y: e.clientY})
		const faces = [a,b,c,d].mapFilter(b2 => b2 && b2.faces).concatenated()
		const testEdges: Edge[] = [a,b,c,d].mapFilter(b2 => b2 && (b2.buildAdjacencies(), Array.from(b2.edgeFaces.keys())))
			.concatenated()
			.concat(edges)
		hovering = getHovering(mouseLine, faces, undefined, [], testEdges, 0.1, 'faces', 'edges')
		let html = '', pp
		if (hovering instanceof Edge) {
			pp = V({x: e.clientX, y: e.clientY})
			NLA.defaultRoundFunction = x => NLA.round10(x, -3)
			html = hovering.toString(x => NLA.round10(x, -3)) + ' length=' + hovering.length().toFixed(3)
		} else if (hovering instanceof Face) {
			pp = V({x: e.clientX, y: e.clientY})
			NLA.defaultRoundFunction = x => NLA.round10(x, -3)
			let area, f = hovering
			try { area = hovering.calcArea()} catch (e) {}
			html = `face surface=${f.surface.constructor.name} edges=${f.contour.length} area=${area}`
		}
		if (pp) {
			//const pSC = gl.projectionMatrix.times(gl.modelViewMatrix).transformPoint(pp)
			//const x = (pSC.x * 0.5 + 0.5) * window.innerWidth, y = (-pSC.y * 0.5 + 0.5) * window.innerHeight
			tooltipShow(html, pp.x, pp.y)
		} else {
			tooltipHide()
		}
		paintScreen()
	})
}


//var sketchPlane = new CustomPlane(V3.X, V3(1, 0, -1).unit(), V3.Y, -500, 500, -500, 500, 0xff00ff);
const b2planes = [
    new CustomPlane(V3.O, V3.Y, V3.Z, 'planeYZ', 0xff0000),
    new CustomPlane(V3.O, V3.X, V3.Z, 'planeZX', 0x00ff00),
    new CustomPlane(V3.O, V3.X, V3.Y, 'planeXY', 0x0000ff),
    //	sketchPlane
]

async function viewerMain() {
	paintScreen = viewerPaint
	await B2T.loadFonts()
	modeStack.push({})
	window.onerror = function (errorMsg, url, lineNumber, column, errorObj) {
		console.log(errorMsg, url, lineNumber, column, errorObj)
	}
	gl = GL.create({canvas: document.getElementById('testcanvas') as HTMLCanvasElement})
	gl.fullscreen()
	gl.canvas.oncontextmenu = () => false

	setupCamera(eye, gl)
	//gl.cullFace(gl.FRONT_AND_BACK);
	gl.clearColor(1.0, 1.0, 1.0, 0.0)
	gl.enable(gl.BLEND)
	gl.enable(gl.DEPTH_TEST)
	gl.enable(gl.CULL_FACE)
	gl.depthFunc(gl.LEQUAL)
	gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA) // TODO ?!

	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
	gl.loadIdentity()
	gl.scale(10, 10, 10)

	gl.loadIdentity()


	initMeshes(gl.meshes = meshes)
	initShaders(gl.shaders = shaders)
	initNavigationEvents(gl, eye, paintScreen)
	initInfoEvents()
	//initPointInfoEvents()
	initB2()
	setupCamera(eye, gl)
	paintScreen()
}