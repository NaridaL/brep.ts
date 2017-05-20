///<reference path="surface/ConicSurface.ts"/>
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
let a: B2, b: B2, c: B2, d: B2, edges: Edge[] = [], hovering
const edgeViewerColors = arrayFromFunction(20, i => chroma.random().gl())
function initB2() {
    dMesh = new Mesh()
    eye.pos = V(1, 2, 101)
	eye.focus = V(0, 1, 0)
	eye.up = V(0, 1, 0)
    eye.zoomFactor = 8

	const hash = window.location.search.substr(1) || window.location.hash.substr(1)
	//noinspection TsLint
	let i, points, vectors, mesh, hjk
	//Object.assign(window, HJK())
	eval(decodeURI(hash))
	let gets: any = {a, b, c, d, mesh, edges, points, vectors}
	hjk && (gets = HJK())
	i && Object.assign(eye, i)
	'abcd'.split('').forEach(k => gets[k] && eval(`${k}Mesh = (${k} = gets.${k}).toMesh()`))

    //cMesh && cMesh.computeWireframeFromFlatTriangles() && cMesh.compile()
    if (gets.points) {
        console.log('drPs from GET')
        drPs = gets.points
    }
	if (gets.vectors) {
		console.log('vectors from GET')
		drVs.pushAll(gets.vectors)
	}

	//cMesh && cMesh.computeNormalLines(0.1) && cMesh.compile()
	//aMesh && aMesh.computeNormalLines(0.1) && aMesh.compile()

    if (gets.edges) {
        console.log('edges from GET')
        dMesh = new Mesh({triangles: true, normals: true})
        edges = gets.edges
        edges && dMesh.addVertexBuffer('curve1', 'curve1')
        edges && dMesh.addVertexBuffer('curve1colors', 'curve1colors')
        edges.forEach((edge, edgeIndex) => {
            const points = edge.points()
            for (let i = 0; i < points.length - 1; i++) {
	            const color = edgeViewerColors[(edgeIndex + (i % 2)) % edgeViewerColors.length]
	            // const tangent = edge.tangentAt(i)
                // dMesh.curve1.push(points[i], points[i].plus(tangent.toLength(1)))
                dMesh.curve1.push(points[i], points[i + 1])
                dMesh.curve1colors.push(color, color)
            }
	        ;(edge.curve instanceof PICurve) && edge.curve.addToMesh(dMesh, 8, 0.02, 2, edge.minT, edge.maxT)
        })
	    //dMesh.computeWireframeFromFlatTriangles()
    }
    if (gets.mesh) {
	    console.log('mesh/es from GET', b2meshes)
	    b2meshes = gets.mesh instanceof Array ? gets.mesh : [gets.mesh]
	    b2meshes.forEach(m => m.computeWireframeFromFlatTriangles())
	    b2meshes.forEach(m => m.computeNormalLines(0.5))
	    b2meshes.forEach(m => m.compile())
    }

    dMesh.compile()
}





const randomColors = chroma.scale(['#ff297f', '#6636FF']).mode('lab').colors(20).map(s => chroma(s).gl())
let aMesh: Mesh & {faceIndexes?: Map<Face, {start: int, count: int}>},
	bMesh: Mesh & {faceIndexes?: Map<Face, {start: int, count: int}>},
	cMesh: Mesh & {faceIndexes?: Map<Face, {start: int, count: int}>},
	dMesh: Mesh & {faceIndexes?: Map<Face, {start: int, count: int}>},
	sMesh: Mesh,
	b2meshes: Mesh[] = []
function viewerPaint() {
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
    gl.loadIdentity()

    drawVectors()

    //gl.scale(100, 100, 100)

    if (aMesh) {
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        aMesh.lines && shaders.singleColor.uniforms({ color: hexIntToGLColor(COLORS.PP_STROKE) }).draw(aMesh, DRAW_MODES.LINES)
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
        shaders.lighting.uniforms({ color: hexIntToGLColor(COLORS.PP_FILL),
            camPos: eye.pos }).draw(aMesh)
    }

    if (bMesh) {
        gl.pushMatrix()
        //gl.translate(15, 0, 0)
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        bMesh.lines && shaders.singleColor.uniforms({ color: hexIntToGLColor(COLORS.PP_STROKE) }).draw(bMesh, DRAW_MODES.LINES)
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
        shaders.lighting.uniforms({ color: hexIntToGLColor(COLORS.RD_FILL),
            camPos: eye.pos }).draw(bMesh)
        bMesh.edgeTangents && shaders.singleColor.uniforms({ color: hexIntToGLColor(COLORS.TS_STROKE) })
            .drawBuffers({LGL_Vertex: bMesh.vertexBuffers.edgeTangents}, null, gl.LINES)
        bMesh.edgeTangents2 && shaders.singleColor.uniforms({ color: hexIntToGLColor(COLORS.RD_STROKE) })
            .drawBuffers({LGL_Vertex: bMesh.vertexBuffers.edgeTangents2}, null, gl.LINES)
        gl.popMatrix()
    }
    if (cMesh) {
        gl.pushMatrix()
        //gl.translate(30, 0, 0)
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        cMesh.lines && shaders.singleColor.uniforms({ color: hexIntToGLColor(COLORS.TS_STROKE) }).draw(cMesh, DRAW_MODES.LINES)
        gl.projectionMatrix.m[11] += 1 / (1 << 20)

        let faceIndex = c.faces.length
        while (faceIndex--) {

            const face = c.faces[faceIndex]
            const faceTriangleIndexes = cMesh.faceIndexes.get(face)
            shaders.lighting.uniforms({
                color: hovering == face ? chroma('purple').gl() : randomColors[faceIndex % randomColors.length]
            }).draw(cMesh, DRAW_MODES.TRIANGLES, faceTriangleIndexes.start, faceTriangleIndexes.count)
        }

        gl.popMatrix()
    }
    if (dMesh && dMesh.hasBeenCompiled) {
        gl.pushMatrix()
        //gl.scale(10, 10, 10)
        gl.projectionMatrix.m[11] -= 1 / (1 << 21) // prevent Z-fighting
        dMesh.lines && shaders.singleColor.uniforms({ color: hexIntToGLColor(COLORS.RD_STROKE) }).draw(dMesh, DRAW_MODES.LINES)
        gl.projectionMatrix.m[11] += 1 / (1 << 21)
        dMesh.triangles && shaders.lighting.uniforms({ color: chroma('#ff8b23').gl(),
            camPos: eye.pos }).draw(dMesh)

        dMesh.curve1 && shaders.multiColor
            .drawBuffers({LGL_Vertex: dMesh.vertexBuffers.curve1, color: dMesh.vertexBuffers.curve1colors}, null, gl.LINES)
        dMesh.curve2 && shaders.singleColor.uniforms({ color: hexIntToGLColor(COLORS.RD_STROKE) })
            .drawBuffers({LGL_Vertex: dMesh.vertexBuffers.curve2}, null, gl.LINES)

        dMesh.curve3 && shaders.singleColor.uniforms({ color: hexIntToGLColor(COLORS.PP_STROKE) })
            .drawBuffers({LGL_Vertex: dMesh.vertexBuffers.curve3}, null, gl.LINES)
        dMesh.curve4 && shaders.singleColor.uniforms({ color: hexIntToGLColor(0x00ff00) })
            .drawBuffers({LGL_Vertex: dMesh.vertexBuffers.curve4}, null, gl.LINES)
        gl.popMatrix()
    }
    if (sMesh && sMesh.hasBeenCompiled) {
        gl.pushMatrix()
        //gl.scale(10, 10, 10)
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        sMesh.lines && shaders.singleColor.uniforms({ color: hexIntToGLColor(0xFF6600) }).draw(sMesh, DRAW_MODES.LINES)
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
        sMesh.triangles && shaders.lighting.uniforms({ color: hexIntToGLColor(0xffFF00),
            camPos: eye.pos }).draw(sMesh)

        sMesh.curve1 && shaders.singleColor.uniforms({ color: hexIntToGLColor(0xff00000) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve1}, null, gl.LINES)
        sMesh.curve2 && shaders.singleColor.uniforms({ color: hexIntToGLColor(COLORS.RD_STROKE) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve2}, null, gl.LINES)

        sMesh.curve3 && shaders.singleColor.uniforms({ color: hexIntToGLColor(COLORS.PP_STROKE) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve3}, null, gl.LINES)
        sMesh.curve4 && shaders.singleColor.uniforms({ color: hexIntToGLColor(0x00ff00) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve4}, null, gl.LINES)
        gl.popMatrix()
    }
    for (const sMesh of b2meshes) {
        gl.pushMatrix()
        //gl.scale(10, 10, 10)
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        sMesh.lines && shaders.singleColor.uniforms({ color: hexIntToGLColor(0xFF6600) }).draw(sMesh, DRAW_MODES.LINES)
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
        sMesh.triangles && shaders.lighting.uniforms({ color: hexIntToGLColor(0xffFF00),
            camPos: eye.pos }).draw(sMesh)

        sMesh.curve1 && shaders.singleColor.uniforms({ color: hexIntToGLColor(0xff00000) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve1}, null, gl.LINES)
        sMesh.curve2 && shaders.singleColor.uniforms({ color: hexIntToGLColor(COLORS.RD_STROKE) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve2}, null, gl.LINES)

        sMesh.curve3 && shaders.singleColor.uniforms({ color: hexIntToGLColor(COLORS.PP_STROKE) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve3}, null, gl.LINES)
        sMesh.curve4 && shaders.singleColor.uniforms({ color: hexIntToGLColor(0x00ff00) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve4}, null, gl.LINES)
        gl.popMatrix()
    }

    if (hovering instanceof Edge) {
	    gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
	    drawEdge(hovering, 0x000000, 0.01)
	    gl.projectionMatrix.m[11] += 1 / (1 << 20)
    }
    //edges.forEach((e, i) => drawEdge(e, 0x0ff000, 0.01))

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
			pp = V(e.clientX, e.clientY)
			defaultRoundFunction = x => round10(x, -3)
			html = hovering.toString(x => round10(x, -3)) + ' length=' + hovering.length().toFixed(3)
		} else if (hovering instanceof Face) {
			pp = V(e.clientX, e.clientY)
			defaultRoundFunction = x => round10(x, -3)
			let area, f = hovering
			try { area = hovering.calcArea() } catch (e) {}
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
	gl = create({canvas: document.getElementById('testcanvas') as HTMLCanvasElement})
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
	setupCameraListener = function (eye) {
		const iSource = 'i=' + eye.toSource().replace(/[\n\r\s]+|^\(|\)$/g, '')
		const hash = window.location.hash.substr(1) || iSource
		const result = hash.match(/i=\{[^}]*\}/)
			? hash.replace(/i=\{[^}]*\}/, iSource)
			: hash + ';' + iSource
		window.history.replaceState(undefined, undefined, '#' + result)
	}
	initInfoEvents()
	//initPointInfoEvents()
	initB2()
	setupCamera(eye, gl)
	paintScreen()
}


function test2() {
	const ic = (x, y) => sin(x+y)-cos(x*y)+1
	const dids = (x, y) => y * sin(x * y) + cos(x + y)
	const didt = (x, y) => x * sin(x * y) + cos(x + y)
	const ic2 = (x, y) => (3 * x ** 2 - y ** 2) ** 2 * y ** 2 - (x ** 2 + y ** 2) ** 4
	const di2ds = (x, y) => 4* x* (9* x**2* y**2 - 3* y**4 - 2* (x**2 + y**2)**3)
	const di2dt = (x, y) => 2 * y * (-4 * (x ** 2 + y ** 2) ** 3 + (3 * x ** 2 - y ** 2) ** 2 + 2 * y ** 2 * (y ** 2 - 3 * x ** 2))
	const start = V(-3.6339970071165784, 3.5625834844534974, 0) // curvePoint(ic, V(-4, 4))
	assert(eq02(ic(start.x, start.y), 0.1))
	const bounds = (s, t) => -5 <= s && s <= 5 && -5 <= t && t <= 5
	//const curves =  Curve.breakDownIC(ic, -5, 5, -5, 5, 0.1, 0.1, 0.05, dids, didt)
	const curves =  Curve.breakDownIC(ic2, {sMin: -5, sMax: 5, tMin: -5, tMax: 5}, 0.1, 0.1, 0.02, di2ds, di2dt)
	//const curves =  Curve.breakDownIC(cassini(1, 1.02), -5, 5, -5, 5, 0.1, 0.1, 0.02)
	//const curves = mkcurves(ic, start.x, start.y, 0.05, dids, didt, bounds)
		.map(({points, tangents}, i) => {
			const curve = new ImplicitCurve(ic, points, tangents)
			return Edge.forCurveAndTs(curve.translate(5, 0, 0.1 * i))
		})
	//checkDerivate(s => ic(s, 0), s => dids(s, 0), -5, 5, 0)
	//checkDerivate(t => ic(0, t), t => dids(0, t), -5, 5, 0)
	console.log(curves.length)
	return curves

}
function cassini(a: number, c: number): (x: number, y: number) => number {
	return (x, y) => (x * x + y * y) * (x * x + y * y) - 2 * c * c * (x * x - y * y) - (a ** 4 - c ** 4)
}

/**
 * A function R² -> R with first and second derivatives.
 */
interface MathFunctionR2_R {
	(s: number, t: number): number
	readonly x: R2_R
	readonly y: R2_R
	readonly xx?: R2_R
	readonly xy?: R2_R
	readonly yy?: R2_R
}
namespace MathFunctionR2_R {
	export function forNerdamer(expression: nerdamer.ExpressionParam, args: [string, string] = ['x', 'y']): MathFunctionR2_R {
		const ndf = nerdamer(expression)
		const ndfs = nerdamer.diff(ndf, args[0])
		const ndft = nerdamer.diff(ndf, args[1])
		const f = ndf.buildFunction(args) as any
		f.x = ndfs.buildFunction(args)
		f.y = ndft.buildFunction(args)
		f.xx = nerdamer.diff(ndfs, args[0]).buildFunction(args)
		f.xy = nerdamer.diff(ndfs, args[1]).buildFunction(args)
		f.yy = nerdamer.diff(ndft, args[1]).buildFunction(args)
		return f
	}

	export function nerdamerToR2_R(expression: nerdamer.Expression, args: [string, string] = ['x', 'y']) {
		return expression.buildFunction(args)
	}

	export function forFFxFy(f: R2_R, fx: R2_R, fy: R2_R): MathFunctionR2_R {
		0;(f as any).x = fx
		0;(f as any).y = fy
		return f as any
	}
}
const cas2 = cassini(0.9, 1.02)
function HJK() {
	const x = math.compile('x^2')
	//math.derivative('x^2', 'x')
	nerdamer.setFunction('cassini', 'acxy'.split(''), '(x^2 + y^2)^2 + 2 c^2 (x^2 - y^2) - (a^4 - c^4)')
	const cassini = nerdamer('(x^2 + y^2)^2 + 2 c^2 (x^2 - y^2) - (a^4 - c^4)')
	const ndf = nerdamer('cassini(1, 1, x, y)')
	const mf = MathFunctionR2_R.forNerdamer(ndf)
	//const pf = cs.parametricFunction(), icc = ses.implicitFunction()
	//const start = V(-3.6339970071165784, 3.5625834844534974, 0) // curvePoint(ic, V(-4, 4))
	//assert(eq02(ic(start.x, start.y), 0.1))

	const bounds2={sMin: -2, sMax : 2, tMin : -2, tMax : 2}
	const {sMin, tMin, sMax, tMax} = bounds2
	const bounds = (s, t) => sMin <= s && s <= sMax && tMin <= t && t <= tMax
	//const curves =  Curve.breakDownIC(ic, -5, 5, -5, 5, 0.1, 0.1, 0.05, dids, didt)
	const {points, tangents}
		= followAlgorithm2d(mf, curvePointMF(mf, V(1, 0.5)), 0.02, bounds)
	//const edges = [Edge.forCurveAndTs(new ImplicitCurve(points, tangents).scale(10))]
	const edges =  Curve.breakDownIC(mf, bounds2, 0.1, 0.1, 0.05)
		.map(({points, tangents}, i) => Edge.forCurveAndTs(new ImplicitCurve(points, tangents, 1)).translate(0,0,i*0.1))
	//const curves =  Curve.breakDownIC(cassini(1, 1.02), -5, 5, -5, 5, 0.1, 0.1, 0.02)
	//const curves = mkcurves(ic, start.x, start.y, 0.05, dids, didt, bounds)
	console.log(edges.length)
	return {edges: edges, points: edges[0].curve.points}

}
function HJK_() {
	//return {}
	//const a = B2T.cone(12, 12).scale(0.05,0.2)
	//	.rotateZ(90*DEG)
	//	.rotateY(-90*DEG)
	//	.translate(2,0.2,0.7)
	//const b = B2T.sphere(1, undefined, PI)
	//const cone = new ConicSurface(
	//	V(2, 0.2, 1.1),
	//	V(0, 0.6, 0),
	//	V(0, 0, -2.4),
	//	V(-12, 0, 0))
	//const sphere = new SemiEllipsoidSurface(V3.O,V3.X,V3.Y,V(0, 0, -1))
	////const curves = cone.isCurvesWithSurface(sphere)
	////assert(sphere.containsCurve(curves[0]))
	//const c = a.minus(b).translate(0,3)
	//const pcs = new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -200, 200)
	const pcs = SemiCylinderSurface.UNIT.scale(0.5, 0.05, 4).translate(0.5,0,-2).flipped()
	const ses = new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z)
	const curves = ses.isCurvesWithSurface(pcs)
	return {
		//a,
		//b,
		//c,
		mesh: [pcs.toMesh(), ses.toMesh()],
		//edges: [Edge.forCurveAndTs(HyperbolaCurve.XY.shearedX(2, 3))]
	}

}