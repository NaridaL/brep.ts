import chroma, { Color } from 'chroma-js'
import nerdamer from 'nerdamer'
import { arrayFromFunction, assert, DEG, int, round10, V, V3} from 'ts3dutils'
import { GL_COLOR, GL_COLOR_BLACK, Mesh, TSGLContext } from 'tsgl'
import deepmerge from 'deepmerge'

import {
	BRep, B2T, BREPGLContext, cameraChangeListeners, COLORS, Curve, curvePointMF, CustomPlane, Edge, Face, followAlgorithm2d, getMouseLine, ImplicitCurve, initNavigationEvents, L3,
	MathFunctionR2R,
	P3,
	PICurve,
	SemiCylinderSurface,
	SemiEllipsoidSurface,
	setupCamera,
} from './index'

const eye = { pos: V(1000, 1000, 1000), focus: V3.O, up: V3.Z, zoomFactor: 1 }
const drVs: any[] = []
const b2s: BRep[] = []
const edgeViewerColors = arrayFromFunction(20, i => chroma.random().gl())
const aMeshes: (Mesh & { faceIndexes?: Map<Face, { start: int, count: int }>, TRIANGLES: int[], normals: V3[] })[] = []
//bMesh: Mesh & {faceIndexes?: Map<Face, {start: int, count: int}>},
//cMesh: Mesh & {faceIndexes?: Map<Face, {start: int, count: int}>},
let dMesh: Mesh & { faceIndexes?: Map<Face, { start: int, count: int }>, TRIANGLES: int[], normals: V3[], curve1: V3[], curve1colors: GL_COLOR[] }
let faceMesh: Mesh & {tangents: V3[], TRIANGLES: int[], normals: V3[]}
let b2meshes: (Mesh & { TRIANGLES: int[], normals: V3[] })[] = []
let hovering: {}

import * as ts3dutils from 'ts3dutils'
import * as tsgl from 'tsgl'
import * as brepts from './index'
const addMissing = (to: any, from: any) => Object.keys(from).forEach(key => 'Buffer' != key && !to[key] && (to[key] = from[key]))
// tslint:disable-next-line:class-name
export class RenderObjects {
	a: BRep = undefined
	b: BRep = undefined
	c: BRep = undefined
	d: BRep = undefined
    face: Face[]
	edges: Edge[] = []
	wireframe: boolean = false
	normallines: boolean = false
	i: any = undefined
	hjk: any = undefined
	drPs:  (V3 | { info: string, p: V3 })[] = []
	drVs: any = []
	mesh: (Mesh & { TRIANGLES: int[], normals: V3[] }) = undefined
}
const renderObjectKeys = Object.keys(new RenderObjects()) as (keyof RenderObjects)[]

declare function INIT_HTML(): void
addMissing(window, ts3dutils)
addMissing(window, tsgl)
addMissing(window, brepts)
addMissing(window, new RenderObjects())
declare global {
    interface Window extends RenderObjects {}
}
const arrayLiteralType = <T extends string>(x: T[]): T[] => x
const g = window
function objectAssignConcatArray<T, U>(a: T, b: U): T & U {
	for (const key of Object.keys(b)) {
		if (Array.isArray(g[key]) && Array.isArray(b[key])) {
			a[key].push(...b[key])
		} else if (undefined !== b[key]) {
			a[key] = b[key]
		}
	}
	return a as any
}
function initBRep() {
	eye.pos = V(1, 2, 101)
	eye.focus = V(0, 1, 0)
	eye.up = V(0, 1, 0)
	eye.zoomFactor = 8

    const htmlContext = INIT_HTML()

    const hash = window.location.search.substr(1) || window.location.hash.substr(1) || ''
    const command = decodeURIComponent(hash)
    console.log(command)
    const hashContext = new Function(`let ${renderObjectKeys.join(',')};${command};return{${renderObjectKeys.join(',')}}`)() as RenderObjects

	// hashContext last, so i value in hash wins
	objectAssignConcatArray(g, htmlContext)
	objectAssignConcatArray(g, hashContext)
	console.log(htmlContext)


	Object.assign(eye, g.i)
	// let gets: any = {a, b, c, d, mesh, edges, points, vectors}
	g.hjk && Object.assign(g, HJK())
    arrayLiteralType(['a', 'b', 'c', 'd']).forEach(k => {
        if (g[k]) {
            aMeshes.push(g[k].toMesh())
            b2s.push(g[k])
        }
    })

	for (let i = 0; i < aMeshes.length; i++) {
		aMeshes[i].computeWireframeFromFlatTriangles('wireframe')
		aMeshes[i].computeNormalLines(0.1, 'normallines')
		aMeshes[i].compile()
	}

	if (g.edges) {
		console.log('edges from GET')
		dMesh = new Mesh()
			.addIndexBuffer('TRIANGLES')
			.addVertexBuffer('normals', 'ts_Normal')
			.addVertexBuffer('curve1', 'curve1')
			.addVertexBuffer('curve1colors', 'curve1colors')
		g.edges.forEach((edge, edgeIndex) => {
			const points = edge.points()
			for (let i = 0; i < points.length - 1; i++) {
				const color = edgeViewerColors[(edgeIndex + (i % 2)) % edgeViewerColors.length]
				// const tangent = edge.tangentAt(i)
				// dMesh.curve1.push(points[i], points[i].plus(tangent.toLength(1)))
				dMesh.curve1.push(points[i], points[i + 1])
				dMesh.curve1colors.push(color, color)
			}
			edge.curve instanceof PICurve && (edge.curve as PICurve).addToMesh(dMesh, 8, 0.02, 2)
		})
		//dMesh.computeWireframeFromFlatTriangles()
	}
	if (g.mesh) {
		console.log('mesh/es from GET', b2meshes)
		b2meshes = g.mesh instanceof Array ? g.mesh : [g.mesh]
		b2meshes.forEach(m => m.computeWireframeFromFlatTriangles())
		b2meshes.forEach(m => m.computeNormalLines(0.5))
		b2meshes.forEach(m => m.compile())
	}
	if (g.face) {
        if (!g.face.length) {
            g.face = [g.face] as any
        }
        faceMesh = new Mesh()
            .addIndexBuffer('TRIANGLES')
            .addIndexBuffer('LINES')
            .addVertexBuffer('tangents', 'tangents')
            .addVertexBuffer('normals', 'ts_Normal')
        for (const face of g.face) {
            face.addToMesh(faceMesh)
            for (const edge of face.allEdges) {
                const ts = edge.curve.calcSegmentTs(edge.aT, edge.bT, edge.reversed, true)
                for (const t of ts) {
                    const p = edge.curve.at(t)
                    faceMesh.tangents.push(p, p.plus(edge.tangentAt(t)))
                }
            }
        }
        faceMesh.compile()
    }

	dMesh && dMesh.compile()
}

const meshColors: Color[][] = [
	chroma.scale(['#ff297f', '#6636FF']).mode('lab').colors(20, null as 'alpha'),
	chroma.scale(['#ffe93a', '#ff6e35']).mode('lab').colors(20, null as 'alpha'),
	chroma.scale(['#1eff33', '#4960ff']).mode('lab').colors(20, null as 'alpha'),
	chroma.scale(['#31fff8', '#2dff2a']).mode('lab').colors(20, null as 'alpha'),
]
const meshColorssGL = meshColors.map(cs => cs.map(c => c.gl()))

function viewerPaint(time: int, gl: BREPGLContext) {
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
	gl.loadIdentity()

	gl.drawVectors(drVs)
	gl.shaders.lighting.uniforms({ camPos: eye.pos })
	//gl.scale(100, 100, 100)
	for (let i = 0; i < aMeshes.length; i++) {
		const aMesh = aMeshes[i]
		gl.pushMatrix()
		//gl.translate(30, 0, 0)
		gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
		aMesh.indexBuffers.wireframe && gl.shaders.singleColor.uniforms({ color: COLORS.TS_STROKE.gl() })
			.drawBuffers(aMesh.vertexBuffers, aMesh.indexBuffers.wireframe, gl.LINES)
		aMesh.indexBuffers.normallines && gl.shaders.singleColor.uniforms({ color: COLORS.TS_STROKE.gl() })
			.drawBuffers(aMesh.vertexBuffers, aMesh.indexBuffers.normallines, gl.LINES)
		gl.shaders.singleColor.uniforms({ color: COLORS.TS_STROKE.gl() })
			.drawBuffers(aMesh.vertexBuffers, aMesh.indexBuffers.LINES, gl.LINES)
		gl.projectionMatrix.m[11] += 1 / (1 << 20)

		let faceIndex = b2s[i].faces.length
		while (faceIndex--) {
			const face = b2s[i].faces[faceIndex]
			const faceTriangleIndexes = aMesh.faceIndexes.get(face)
			gl.shaders.lighting.uniforms({
				color: hovering == face
					? meshColors.emod(i).emod(faceIndex).darken(2).gl()
					: meshColorssGL.emod(i).emod(faceIndex),
			}).draw(aMesh, gl.TRIANGLES, faceTriangleIndexes.start, faceTriangleIndexes.count)
		}

		gl.popMatrix()
	}

	if (faceMesh) {
	    gl.shaders.singleColor.uniforms({color: chroma('red').gl() })
            .drawBuffers({ts_Vertex: faceMesh.vertexBuffers['tangents']}, undefined, gl.LINES)
    }

    //if (dMesh) {
	 //   gl.shaders.multiColor.uniforms({color: COLORS.PP_STROKE.gl() }).drawBuffers({
    //        ts_Vertex: dMesh.vertexBuffers.curve1,
    //        color: dMesh.vertexBuffers.curve1colors,
    //    }, undefined, gl.LINES)
    //}

	//gl.disable(gl.CULL_FACE)
	for (const sMesh of b2meshes) {
		gl.pushMatrix()
		//gl.scale(10, 10, 10)
		gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
		sMesh.LINES && gl.shaders.singleColor.uniforms({ color: chroma('#FF6600').gl() }).draw(sMesh, gl.LINES)
		gl.projectionMatrix.m[11] += 1 / (1 << 20)
		sMesh.TRIANGLES && gl.shaders.lighting.uniforms({
			color: chroma('#ffFF00').gl(),
			camPos: eye.pos,
		}).draw(sMesh)
		gl.popMatrix()
	}
    gl.enable(gl.CULL_FACE)

	if (hovering instanceof Edge) {
		gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
		gl.drawEdge(hovering, GL_COLOR_BLACK, 2 / eye.zoomFactor)
		gl.projectionMatrix.m[11] += 1 / (1 << 20)
	}
	g.edges && g.edges.forEach((e, i) => gl.drawEdge(e, edgeViewerColors.emod(i), 0.01))

	//drPs.forEach(v => drawPoint(v, undefined, 0.3))
	g.drPs.forEach(info => gl.drawPoint(info instanceof V3 ? info : info.p, chroma('#cc0000').gl(), 5 / eye.zoomFactor))
	b2planes.forEach(plane => gl.drawPlane(plane, chroma(plane.color).gl(), hovering == plane))

    gl.begin(gl.LINE_STRIP)
    gl.color('red')
    gl.end()
}

// let meshes: any = {}

function getHovering(mouseLine: L3, faces: Face[], planes: CustomPlane[], points: V3[], edges: Edge[],
	mindist: number,
	...consider: ('faces' | 'planes' | 'sketchElements' | 'points' | 'edges' | 'features')[]): any {
	let hoverHighlight = null, nearest = Infinity

	const checkFeatures = consider.includes('features')
	assert(!checkFeatures || !consider.includes('faces'))

	function checkEl(el: any, distance: number) {
		if (distance < nearest) {
			nearest = distance
			hoverHighlight = el
		}
	}

	if (faces && (consider.includes('faces') || consider.includes('features'))) {
		for (const face of faces) {
			checkEl(checkFeatures ? face.info.feature : face, face.intersectsLine(mouseLine))
		}
	}
	if (planes && consider.includes('planes')) {
		for (const plane of planes) {
			checkEl(plane, plane.distanceTo(mouseLine, mindist))
		}
	}
	if (consider.includes('points')) {
		for (const p of points) {
			const t = mouseLine.pointT(p)
			if (mouseLine.at(t).distanceTo(p) < mindist * 1.2) {
				checkEl(p, t - 0.1)
			}
		}
	}
	if (consider.includes('edges')) {
		const projPlane = new P3(mouseLine.dir1, 0)
		const projPoint = projPlane.projectedPoint(mouseLine.anchor)
		for (const edge of edges) {
			const curve = edge.curve
			const prio = 0.05
			if (curve instanceof L3 && curve.dir1.isParallelTo(mouseLine.dir1)) {
				const d = mouseLine.distanceToPoint(edge.a)
				const t = mouseLine.pointT(edge.a)

				if (d < mindist) {
					checkEl(edge, t - prio)
				}
			} else {
				if (!(curve instanceof ImplicitCurve)) {
					const projCurve = curve.project(projPlane)
					const curveT = edge.clampedT(projCurve.closestTToPoint(projPoint))
					const p = curve.at(curveT)
					const t = mouseLine.pointT(p)
					if (projCurve.at(curveT).distanceTo(projPoint) < mindist) {
						checkEl(edge, t - prio)
					}
				}
			}
		}
	}

	return hoverHighlight
}

function initInfoEvents(paintScreen: () => {}, gl: BREPGLContext) {
	gl.canvas.addEventListener('mousemove', function (e) {
		const mouseLine = getMouseLine({ x: e.clientX, y: e.clientY }, gl)
		const faces = b2s.flatMap(b2 => b2 && b2.faces)
		const testEdges: Edge[] = [
			...b2s.flatMap(b2 => Array.from<Edge>(b2.buildAdjacencies().edgeFaces.keys())),
			...g.edges]
		hovering = getHovering(mouseLine, faces, undefined, [], testEdges, 0.1, 'faces', 'edges')
		//let html = '', pp
		//if (hovering instanceof Edge) {
		//	pp = V(e.clientX, e.clientY)
		//	defaultRoundFunction = x => round10(x, -3)
		//	html = hovering.toString(x => round10(x, -3)) + ' length=' + hovering.length().toFixed(3)
		//} else if (hovering instanceof Face) {
		//	pp = V(e.clientX, e.clientY)
		//	defaultRoundFunction = x => round10(x, -3)
		//   let area
		//   try { area = hovering.calcArea() } catch (e) {}
		//	html = `face surface=${hovering.surface.constructor.name} edges=${hovering.contour.length} area=${area}`
		//}
		//if (pp) {
		//	//const pSC = gl.projectionMatrix.times(gl.modelViewMatrix).transformPoint(pp)
		//	//const x = (pSC.x * 0.5 + 0.5) * window.innerWidth, y = (-pSC.y * 0.5 + 0.5) * window.innerHeight
		//	tooltipShow(html, pp.x, pp.y)
		//} else {
		//	tooltipHide()
		//}
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
declare var BREPTS_ROOT: string
export async function viewerMain() {
	const paintScreen = () => requestAnimationFrame(t => viewerPaint(t, gl))
	B2T.defaultFont = await B2T.loadFont(BREPTS_ROOT + '/fonts/FiraSansMedium.woff')
	window.onerror = function (errorMsg, url, lineNumber, column, errorObj) {
		console.log(errorMsg, url, lineNumber, column, errorObj)
	}
	const gl = BREPGLContext.create(TSGLContext.create({ canvas: document.getElementById('testcanvas') as HTMLCanvasElement }))
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

	initNavigationEvents(gl, eye, paintScreen)
	cameraChangeListeners.push(function (eye) {
		const round = (x: number) => round10(x, -3)
		const roundedEye = { pos: eye.pos.map(round), focus: eye.focus.map(round), up: eye.up.map(round), zoomFactor: round(eye.zoomFactor) }
		const iSource = 'i=' + roundedEye.toSource().replace(/[\n\r\s]+|^\(|\)$/g, '')
		const hash = window.location.hash.substr(1) || iSource
		const result = hash.match(/i=\{[^}]*\}/)
			? hash.replace(/i=\{[^}]*\}/, iSource)
			: hash + ';' + iSource
		window.history.replaceState(undefined, undefined, '#' + result)
	})
	// initInfoEvents(paintScreen, g l)
	//initToolTips() // hide tooltip on mouseover
	//initPointInfoEvents()
	initBRep()
	setupCamera(eye, gl)
	paintScreen()
}

function HJKl() {
	//math.derivative('x^2', 'x')
	nerdamer.setFunction('cassini', 'acxy'.split(''), '(x^2 + y^2)^2 + 2 c^2 (x^2 - y^2) - (a^4 - c^4)')
	const cassini = nerdamer('(x^2 + y^2)^2 + 2 c^2 (x^2 - y^2) - (a^4 - c^4)')
	const ndf = nerdamer('cassini(1, 1, x, y)')
	const mf = MathFunctionR2R.forNerdamer(ndf)
	//const pf = cs.parametricFunction(), icc = ses.implicitFunction()
	//const start = V(-3.6339970071165784, 3.5625834844534974, 0) // curvePoint(ic, V(-4, 4))
	//assert(eq02(ic(start.x, start.y), 0.1))

	const bounds2 = { sMin: -2, sMax: 2, tMin: -2, tMax: 2 }
	const { sMin, tMin, sMax, tMax } = bounds2
	const bounds = (s: number, t: number) => sMin <= s && s <= sMax && tMin <= t && t <= tMax
	//const curves =  Curve.breakDownIC(ic, -5, 5, -5, 5, 0.1, 0.1, 0.05, dids, didt)
	const { points, tangents }
		= followAlgorithm2d(mf, curvePointMF(mf, V(1, 0.5)), 0.05, bounds)
	// const edges = [Edge.forCurveAndTs(new ImplicitCurve(points, tangents).scale(10))]
	const edges = Curve.breakDownIC(mf, bounds2, 0.1, 0.1, 0.1)
		.map(({ points, tangents }, i) => Edge.forCurveAndTs(new ImplicitCurve(points, tangents, 1)).translate(0, 0, i * 0.1).scale(10))
	//const curves =  Curve.breakDownIC(cassini(1, 1.02), -5, 5, -5, 5, 0.1, 0.1, 0.02)
	//const curves = mkcurves(ic, start.x, start.y, 0.05, dids, didt, bounds)
	console.log(edges.length)
	return { edges: edges, points: (edges[0].curve as ImplicitCurve).points }

}

function HJK() {
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
	//const pcs = new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625,
	// -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009,
	// 1.2), 0, 1), V(0, 0, -1), 0, 1, -200, 200)
	const pcs = SemiCylinderSurface.UNIT
		.rotateZ(-40 * DEG)
		.scale(0.5, 0.05, 4)
		.translate(0.5, 0, -2)
		.flipped()
	const ses = new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z)
	const curves = ses.isCurvesWithSurface(pcs)
	return {
		//a,
		//b,
		//c,
		mesh: [pcs.toMesh(), ses.toMesh()],
		edges: curves.map(c => Edge.forCurveAndTs(c)),
		points: curves.flatMap(c => (c as ImplicitCurve).points),
	}

}
