import chroma, { Color } from 'chroma-js'
import deepmerge from 'deepmerge'
import nerdamer from 'nerdamer'
import { AABB, arrayFromFunction, assert, DEG, int, M4, round10, V, V3 } from 'ts3dutils'
import { GL_COLOR, GL_COLOR_BLACK, Mesh, TSGLContext } from 'tsgl'

import {
	B2T,
	BRep,
	BREPGLContext,
	cameraChangeListeners,
	COLORS,
	Curve,
	curvePointMF,
	CustomPlane,
	CylinderSurface,
	Edge,
	EllipsoidSurface,
	Face,
	FaceMesh,
	followAlgorithm2d,
	getMouseLine,
	ImplicitCurve,
	initNavigationEvents,
	L3,
	MathFunctionR2R,
	P3,
	PICurve,
	setupCamera,
} from './index'

const eye = { pos: V(1000, 1000, 1000), focus: V3.O, up: V3.Z, zoomFactor: 1 }
const bReps: BRep[] = []
const edgeViewerColors = ['darkorange', 'darkgreen', 'cyan'].map(c => chroma(c).gl())
let bRepMeshes: (Mesh & { faceIndexes?: Map<Face, { start: int; count: int }>; TRIANGLES: int[]; normals: V3[] })[] = []
//bMesh: Mesh & {faceIndexes?: Map<Face, {start: int, count: int}>},
//cMesh: Mesh & {faceIndexes?: Map<Face, {start: int, count: int}>},
let edgesMesh: Mesh & {
	faceIndexes?: Map<Face, { start: int; count: int }>
	TRIANGLES: int[]
	normals: V3[]
	curve1: V3[]
	curve1colors: GL_COLOR[]
}
let faceMesh: FaceMesh & { tangents: V3[] }
let meshes: (Mesh & { TRIANGLES: int[]; normals: V3[] })[] = []
let hovering: {}
const edgeDebugPoints = [],
	edgeDebugLines = []

import * as ts3dutils from 'ts3dutils'
import * as tsgl from 'tsgl'
import * as brepts from './index'
const addMissing = (to: any, from: any) =>
	Object.keys(from).forEach(key => 'Buffer' != key && !to[key] && (to[key] = from[key]))
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
	drPs: (V3 | { p: V3; info?: string; color?: string })[] = []
	drVs: { v: V3; anchor: V3; color?: GL_COLOR }[] = []
	drLines: V3[] = []
	mesh: Mesh & { TRIANGLES: int[]; normals: V3[] } = undefined
	aabbs: AABB[] = []
	paintMeshNormals = false
	paintWireframe = false
	paintCurveDebug = false
    planes = []
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
	const hashContext = new Function(
		`let ${renderObjectKeys.join(',')};${command};return{${renderObjectKeys.join(',')}}`,
	)() as RenderObjects

	// hashContext last, so i value in hash wins
	objectAssignConcatArray(g, htmlContext)
	objectAssignConcatArray(g, hashContext)
	console.log(htmlContext)

	Object.assign(eye, g.i)
	// let gets: any = {a, b, c, d, mesh, edges, points, vectors}
	// g.hjk && Object.assign(g, HJK())
	arrayLiteralType(['a', 'b', 'c', 'd']).forEach(k => {
		if (g[k]) {
			bReps.push(g[k])
		}
	})

	bRepMeshes = bReps.map(bRep => bRep.toMesh())
	bRepMeshes.forEach(mesh => {
		mesh.computeWireframeFromFlatTriangles('wireframe')
		mesh.computeNormalLines(0.1, 'normallines')
		mesh.compile()
	})

	if (g.mesh) {
		console.log('mesh/es from GET', bRepMeshes)
		meshes = g.mesh instanceof Array ? g.mesh : [g.mesh]
		meshes.forEach(mesh => {
			mesh.computeWireframeFromFlatTriangles('wireframe')
			mesh.computeNormalLines(0.1, 'normallines')
			mesh.compile()
		})
	}

	if (g.edges) {
		console.log('edges from GET')
		edgesMesh = new Mesh()
			.addIndexBuffer('TRIANGLES')
			.addVertexBuffer('normals', 'ts_Normal')
			.addVertexBuffer('curve1', 'curve1')
			.addVertexBuffer('curve1colors', 'curve1colors')
		g.edges.forEach((edge, edgeIndex) => {
			const points = edge.points()
			for (let i = 0; i < points.length - 1; i++) {
				const color = edgeViewerColors[(edgeIndex + i % 2) % edgeViewerColors.length]
				// const tangent = edge.tangentAt(i)
				// dMesh.curve1.push(points[i], points[i].plus(tangent.toLength(1)))
				edgesMesh.curve1.push(points[i], points[i + 1])
				edgesMesh.curve1colors.push(color, color)
			}
			edge.curve instanceof PICurve && (edge.curve as PICurve).addToMesh(edgesMesh, 8, 0.02, 2)

			if (edge.curve.debugInfo) {
				const { points, lines } = edge.curve.debugInfo()
				points && edgeDebugPoints.push(...points)
				lines && edgeDebugLines.push(...lines)
			}
		})
		//dMesh.computeWireframeFromFlatTriangles()
		edgesMesh.compile()
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

	g.drPs.push()
}

const meshColors: Color[][] = [
	chroma.scale(['#ff297f', '#6636FF']),
	chroma.scale(['#ffe93a', '#ff6e35']),
	chroma.scale(['#1eff33', '#4960ff']),
	chroma.scale(['#31fff8', '#2dff2a']),
].map(scale => scale.mode('lab').colors(20, null as 'alpha'))

const meshColorssGL = meshColors.map(cs => cs.map(c => c.gl()))

function viewerPaint(time: int, gl: BREPGLContext) {
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
	gl.loadIdentity()

	setupCamera(eye, gl)

	gl.drawVectors(g.drVs, 4 / eye.zoomFactor)

	g.drPs.forEach(info =>
		gl.drawPoint(
			info instanceof V3 ? info : info.p,
			info instanceof V3 || !info.color ? chroma('#cc0000').gl() : chroma(info.color).gl(),
			6 / eye.zoomFactor,
		),
	)
	drawPlanes.forEach(plane => gl.drawPlane(plane, plane.color, hovering == plane))
    g.planes.forEach(plane => gl.drawPlane(plane, plane.color, hovering == plane))

	g.aabbs.forEach(aabb => gl.drawAABB(aabb, chroma('black').gl()))

	gl.shaders.lighting.uniforms({ camPos: eye.pos })
	for (let i = 0; i < bRepMeshes.length; i++) {
		const mesh = bRepMeshes[i]
		gl.pushMatrix()
		//gl.translate(30, 0, 0)
		gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
		g.paintWireframe &&
			mesh.indexBuffers.wireframe &&
			gl.shaders.singleColor
				.uniforms({ color: COLORS.TS_STROKE.gl() })
				.drawBuffers(mesh.vertexBuffers, mesh.indexBuffers.wireframe, gl.LINES)
		g.paintMeshNormals &&
			mesh.indexBuffers.normallines &&
			gl.shaders.singleColor
				.uniforms({ color: COLORS.TS_STROKE.gl() })
				.drawBuffers(mesh.vertexBuffers, mesh.indexBuffers.normallines, gl.LINES)
		gl.shaders.singleColor
			.uniforms({ color: COLORS.TS_STROKE.gl() })
			.drawBuffers(mesh.vertexBuffers, mesh.indexBuffers.LINES, gl.LINES)
		gl.projectionMatrix.m[11] += 1 / (1 << 20)

		let faceIndex = bReps[i].faces.length
		while (faceIndex--) {
			const face = bReps[i].faces[faceIndex]
			const faceTriangleIndexes = mesh.faceIndexes.get(face)
			gl.shaders.lighting
				.uniforms({
					color:
						hovering == face
							? meshColors
									.emod(i)
									.emod(faceIndex)
									.darken(2)
									.gl()
							: meshColorssGL.emod(i).emod(faceIndex),
				})
				.draw(mesh, gl.TRIANGLES, faceTriangleIndexes.start, faceTriangleIndexes.count)
		}

		gl.popMatrix()
	}

	if (faceMesh) {
		gl.shaders.singleColor
			.uniforms({ color: chroma('red').gl() })
			.drawBuffers({ ts_Vertex: faceMesh.vertexBuffers.tangents }, undefined, gl.LINES)
	}

	for (const mesh of meshes) {
		gl.pushMatrix()
		gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
		g.paintWireframe &&
			mesh.indexBuffers.wireframe &&
			gl.shaders.singleColor
				.uniforms({ color: COLORS.TS_STROKE.gl() })
				.drawBuffers(mesh.vertexBuffers, mesh.indexBuffers.wireframe, gl.LINES)
		g.paintMeshNormals &&
			mesh.indexBuffers.normallines &&
			gl.shaders.singleColor
				.uniforms({ color: COLORS.TS_STROKE.gl() })
				.drawBuffers(mesh.vertexBuffers, mesh.indexBuffers.normallines, gl.LINES)
		gl.projectionMatrix.m[11] += 1 / (1 << 20)
		mesh.TRIANGLES &&
			gl.shaders.lighting
				.uniforms({
					color: chroma('#ffFF00').gl(),
					camPos: eye.pos,
				})
				.draw(mesh)
		gl.popMatrix()
	}

	if (hovering instanceof Edge) {
		gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
		gl.drawEdge(hovering, GL_COLOR_BLACK, 2 / eye.zoomFactor)
		gl.projectionMatrix.m[11] += 1 / (1 << 20)
	}
	g.edges.forEach((e, i) => gl.drawEdge(e, edgeViewerColors.emod(i), 3 / eye.zoomFactor))

	if (g.paintCurveDebug) {
		gl.begin(gl.LINES)
		gl.color('red')
		edgeDebugLines.forEach(x => gl.vertex(x))
		gl.end()
		edgeDebugPoints.forEach(p => gl.drawPoint(p, chroma('red').gl(), 6 / eye.zoomFactor))
	}
	if (0 !== g.drLines.length) {
		gl.begin(gl.LINES)
		g.drLines.forEach(x => {
			gl.color((x as any).color || 'red')
			gl.vertex(x)
		})
		gl.end()
	}
}

function getHovering(
	mouseLine: L3,
	faces: Face[],
	planes: CustomPlane[],
	points: V3[],
	edges: Edge[],
	mindist: number,
	...consider: ('faces' | 'planes' | 'sketchElements' | 'points' | 'edges' | 'features')[]
): any {
	let hoverHighlight = null,
		nearest = Infinity

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
	gl.canvas.addEventListener('mousemove', function(e) {
		const mouseLine = getMouseLine({ x: e.clientX, y: e.clientY }, gl)
		const faces = bReps.flatMap(b2 => b2 && b2.faces)
		const testEdges: Edge[] = [
			...bReps.flatMap(b2 => Array.from<Edge>(b2.buildAdjacencies().edgeFaces.keys())),
			...g.edges,
		]
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
const drawPlanes = [
	new CustomPlane(V3.O, V3.Y, V3.Z, 'planeYZ', chroma(0xff0000).gl()),
	new CustomPlane(V3.O, V3.X, V3.Z, 'planeZX', chroma(0x00ff00).gl()),
	new CustomPlane(V3.O, V3.X, V3.Y, 'planeXY', chroma(0x0000ff).gl()),
	//	sketchPlane
]
let paintScreen: () => void
declare var BREPTS_ROOT: string
export async function viewerMain() {
	const meshNormalsCheckbox = document.getElementById('paint-mesh-normals')
	meshNormalsCheckbox.onclick = e => {
		g.paintMeshNormals = !g.paintMeshNormals
		paintScreen()
	}

	const wireframeCheckbox = document.getElementById('paint-wireframe')
	wireframeCheckbox.onclick = e => {
		g.paintWireframe = !g.paintWireframe
		paintScreen()
	}

	const paintDebugCheckbox = document.getElementById('paint-curvedebug')
	paintDebugCheckbox.onclick = e => {
		g.paintCurveDebug = !g.paintCurveDebug
		paintScreen()
	}

	paintScreen = () => requestAnimationFrame(t => viewerPaint(t, gl))
	B2T.defaultFont = await B2T.loadFont(BREPTS_ROOT + '/fonts/FiraSansMedium.woff')
	window.onerror = function(errorMsg, url, lineNumber, column, errorObj) {
		console.log(errorMsg, url, lineNumber, column, errorObj)
	}
	const gl = BREPGLContext.create(
		TSGLContext.create({ canvas: document.getElementById('testcanvas') as HTMLCanvasElement }),
	)
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
	cameraChangeListeners.push(function(eye) {
		const round = (x: number) => round10(x, -3)
		const roundedEye = {
			pos: eye.pos.map(round),
			focus: eye.focus.map(round),
			up: eye.up.map(round),
			zoomFactor: round(eye.zoomFactor),
		}
		const iSource = 'i=' + roundedEye.toSource().replace(/[\n\r\s]+|^\(|\)$/g, '')
		const hash = window.location.hash.substr(1) || iSource
		const result = hash.match(/i=\{[^}]*\}/) ? hash.replace(/i=\{[^}]*\}/, iSource) : hash + ';' + iSource
		window.history.replaceState(undefined, undefined, '#' + result)
	})
	// initInfoEvents(paintScreen, g l)
	//initToolTips() // hide tooltip on mouseover
	//initPointInfoEvents()
	initBRep()
	setupCamera(eye, gl)
	paintScreen()
}

export function alignX(dir: number) {
	eye.focus = V3.O
	eye.pos = V(100 * dir, 0, 0)
	eye.up = V3.Z
	paintScreen()
}
export function alignY(dir: number) {
	eye.focus = V3.O
	eye.pos = V(0, 100 * dir, 0)
	eye.up = V3.Z
	paintScreen()
}
export function alignZ(dir: number) {
	eye.focus = V3.O
	eye.pos = V(0, 0, 100 * dir)
	eye.up = eye.pos.cross(V3.X).unit()
	paintScreen()
}
export function rot(angleInDeg: number) {
	eye.up = M4.rotateLine(eye.pos, eye.pos.to(eye.focus), angleInDeg * DEG).transformVector(eye.up)
	paintScreen()
}
