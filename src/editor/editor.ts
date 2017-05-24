///<reference path="Feature.ts"/>
const MooEl: ElementConstructor = Element
function formatStack(stack) {
	const re = /\s*([^@]*)@(.+):(\d+):(\d+)\s*$/mg
	let match, matches = []

	let maxWidths = [0, 0, 0, 0]
	while (match = re.exec(stack)) {
		// console.log(match)
		matches.push(match)
		maxWidths = maxWidths.map((maxWidth, i) => max(maxWidth, match[i + 1].length))
	}
	const output = matches.map(match => {
		return [
			'\t',

			match[1],
			repeatString(maxWidths[0] + 2 - match[1].length, ' '),

			match[2].replace('file:///C:/Users/aval/Desktop/cs', ''),
			repeatString(maxWidths[1] + 2 - match[2].length, ' '),

			repeatString(maxWidths[2] + 2 - match[3].length, ' '),
			match[3],

			repeatString(maxWidths[2] + 2 - match[3].length, ' '),
			match[3]

		].join('')
	}).join('\n')
	return output
}



let savedTextures = new Map()
function getTextureForString(str: string) {
	let texture
	if (texture = savedTextures.get(str)) {
		return texture
	}

	const canvas = $('textTextureCanvas') as any as HTMLCanvasElement
	const ctx = canvas.getContext('2d')

	const font = TEXT_TEXTURE_HEIGHT + 'px Anonymous Pro'

	ctx.font = font
	const textWidthPx = ctx.measureText(str).width
	const canvasWidth = 1 << ceil(log(textWidthPx) / log(2))
	canvas.width = canvasWidth
	canvas.height = TEXT_TEXTURE_HEIGHT


	ctx.font = font
	ctx.fillStyle = '#ffffff'
	ctx.textAlign = 'left'	// This determines the alignment of text, e.g. left, center, right
	ctx.textBaseline = 'top'	// This determines the baseline of the text, e.g. top, middle, bottom
	;(ctx as any).imageSmoothingEnabled = false


	ctx.fillText(str, 0, 0)


	texture = Texture.fromImage(canvas, {minFilter: gl.LINEAR_MIPMAP_LINEAR, magFilter: gl.LINEAR})
	texture.textWidthPx = textWidthPx

	savedTextures.set(str, texture)

	return texture
}
function paintLooseSeg(seg: SketchSegment, color: int) {
	const sketch = seg.sketch
	if (!sketch.plane || !sketch.sketchToWorldMatrix || !sketch.sketchToWorldMatrix.m) return

	gl.pushMatrix()
	gl.multMatrix(sketch.sketchToWorldMatrix)
	paintSeg(seg, color)
	gl.popMatrix()
}
function paintSeg(seg: SketchSegment, color?: int) {
	function drawPoint(p) {
		paintArc(p.V3(), 1.5, 3, colorFor(highlighted.includes(p) || hoverHighlight == p, selected.includes(p)))
	}

	undefined == color && colorFor(highlighted.includes(seg) || hoverHighlight == seg, selected.includes(seg))
	if (seg instanceof SketchLineSeg) {
		//console.log("seg", seg);
		//console.log("hoverHighlight.length", hoverHighlight.length);
		//ctx.beginPath();
		paintLineXY(seg.a.V3(), seg.b.V3(), color)


		//console.log(seg.a);
		drawPoint(seg.a)
		drawPoint(seg.b)
	} else if (seg instanceof SketchArc) {
		const radius = seg.radiusA()
		const angleA = seg.angleA()
        let angleB = seg.angleB()
		if (angleB <= angleA) { angleB += Math.PI * 2 }
		paintArc(seg.c.V3(), radius, 2, color, angleA, angleB)
		drawPoint(seg.a)
		drawPoint(seg.b)
		drawPoint(seg.c)
	} else if (seg instanceof SketchBezier) {
		// paintBezier(seg.points.map(p => p.V3()), 2, 0xdddddd, -2, 3)
		paintBezier(seg.points.map(p => p.V3()), 2, color, 0, 1)
		seg.points.forEach(drawPoint)
	} else if (seg instanceof SegmentEndPoint) {
		drawPoint(seg)
	} else {
		throw new Error('unknown sketch element' + seg)
	}

}
function paintSketch(sketch: Sketch) {

	if (!sketch.plane || !sketch.sketchToWorldMatrix || !sketch.sketchToWorldMatrix.m) return
	//console.log("painting segments", sketch.elements.length);
	/*ctx.clearRect (0, 0, ctx.canvas.width, ctx.canvas.height);
	 ctx.fillStyle="rgb(100, 100, 255)";
	 ctx.lineWidth=2;*/
	//console.log(sketch.sketchToWorldMatrix);
	gl.pushMatrix()
	gl.multMatrix(sketch.sketchToWorldMatrix)
	sketch.elements.forEach(seg => paintSeg(seg))

	paintConstraints(sketch)
	gl.popMatrix()
}
function colorFor(highlighted, selected) {
	return !selected
		? (!highlighted ? 0x33CCFF : 0x145266)
		: (!highlighted ? 0xFF3399 : 0x330A1E)
}
function getAllPoints(sketch) {
	return sketch.elements.flatMap(segment => segment.points)
}
function getCachedCurve(other: NameRef|any, plane: P3) {
	if (other instanceof NameRef) other = other.lastHit
	if (other instanceof Edge) return other.curve
	if (other instanceof SketchLineSeg || other instanceof SketchArc || other instanceof SketchBezier) return other.getCurve()
	if (other instanceof CustomPlane) return other.intersectionWithPlane(plane)
	if (other instanceof Face) other = other.surface
	if (other instanceof Surface) return other.isCurvesWithPlane(plane)[0]
}
function paintConstraints(sketch: Sketch) {
	let crossCount = 2, parallelCrossCount = 1
	sketch.constraints.forEach(function (cst) {
		paintDefaultColor = hoverHighlight == cst ? (!cst.error ? 0x00ff00 : 0xffff00) : (!cst.error ? 0x000000 : 0xff0000)
		switch (cst.type) {
			case 'coincident':
				const point = cst.cs[0]
				paintArc(point.V3(), 4, 1, colorFor(false, false), 0, 2 * Math.PI)
				paintArc(point.V3(), 1.5, 3, colorFor(highlighted.includes(point) || hoverHighlight == point, selected.includes(point)), 0, 2 * Math.PI)
				break
			case 'parallel': {
				const dir1 = cst.segments[0].getVectorAB().unit()
				const dir90 = dir1.getPerpendicular().unit()
				const ARR_SPACING = 5
				for (let c = 0; c < cst.segments.length; c++) {
					const line = cst.segments[c]
					const ab = line.getVectorAB()
					const abLength = ab.length()
					const ab1 = ab.unit()
					for (let i = 0; i < parallelCrossCount; i++) {
						const s = abLength / 2 - ARR_SPACING * parallelCrossCount / 2 + i * ARR_SPACING - 10
						const arrPoint = line.a.V3().plus(ab1.times(s))
						//console.log('crosspos', crossStart, crossEnd);
						paintLineXY(arrPoint.plus(dir90.times(-1)), arrPoint.plus(dir1.times(6).plus(dir90.times(5))))
						paintLineXY(arrPoint.plus(dir90.times(1)), arrPoint.plus(dir1.times(6).plus(dir90.times(-5))))
					}
				}
				parallelCrossCount++
				break
			}
			case 'perpendicular':
				if (cst.other instanceof SketchLineSeg) {
					const ab = cst.segment.getVectorAB()
					const cd = cst.other.getVectorAB()
					const intersection = cst.segment.intersection(cst.other)
					const abPos = cst.segment.pointT(intersection)
					const cdPos = cst.other.pointT(intersection)
					const abLine = ab.unit().times(0.5 < abPos ? -16 : 16)
					const cdLine = cd.unit().times(0.5 < cdPos ? -16 : 16)
					paintLineXY(intersection.plus(abLine).plus(cdLine), intersection.plus(abLine))
					paintLineXY(intersection.plus(abLine).plus(cdLine), intersection.plus(cdLine))
				} else {
					const ab = cst.segment.getVectorAB()
					const sketchLineWC = sketch.plane.intersectionWithPlane(cst.other.plane)
					const dir = sketch.worldToSketchMatrix.transformVector(sketchLineWC.dir1)
					const p = sketch.worldToSketchMatrix.transformPoint(sketchLineWC.anchor)
					const abXcd = ab.cross(dir)
					const div = abXcd.squared()
					const anchorDiff = p.minus(cst.segment.a.V3())
					const abPos = anchorDiff.cross(dir).dot(abXcd) / div
					const linePos = anchorDiff.cross(ab).dot(abXcd) / div
					const intersection = p.plus(dir.times(linePos))
					const abLine = ab.unit().times(0.5 < abPos ? -16 : 16)
					const cdLine = dir.times(16)
					paintLineXY(intersection.plus(abLine).plus(cdLine), intersection.plus(abLine))
					paintLineXY(intersection.plus(abLine).plus(cdLine), intersection.plus(cdLine))
				}
				break
			case 'colinear':
				// assume segments dont overlap
				const segments = cst.segments
				let ab = segments[0].getVectorAB().unit()
				const coord = ab.x > ab.y ? 'x' : 'y'
				if (ab[coord] < 0) {
					ab = ab.times(-1)
				}
				const scale = ab[coord]
				const offsetPoint = segments[0].a.V3().minus(ab.times(segments[0].a.V3()[coord] / scale))
				segments.sort((a, b) => a.a[coord] - b.a[coord])
				for (let i = 0; i < segments.length - (cst.fixed ? 0 : 1); i++) {
					let startS = max(segments[i].a.V3()[coord] / scale, segments[i].b.V3()[coord] / scale) + 8
					let endS = startS + 20
					paintLineXY(offsetPoint.plus(ab.times(startS)), offsetPoint.plus(ab.times(endS)), undefined, 2)

					startS = min(segments[(i + 1) % segments.length].a[coord] / scale, segments[(i + 1) % segments.length].b[coord] / scale) - 8
					endS = startS - 20
					paintLineXY(offsetPoint.plus(ab.times(startS)), offsetPoint.plus(ab.times(endS)), undefined, 2)
				}
				break
			case 'pointDistance': {
				const a = cst.cs[0].V3(), b = cst.cs[1].V3()
				const ab = b.minus(a)
				const ab1 = ab.unit()
				const abLength = ab.length()
				const ab90 = ab1.getPerpendicular()
				const texture = getTextureForString(cst.distance)
				const textLength = texture.textWidthPx / TEXT_TEXTURE_HEIGHT * 20
				paintLineXY(a.plus(ab90.times(6)), a.plus(ab90.times(22)))
				paintLineXY(b.plus(ab90.times(6)), b.plus(ab90.times(22)))

				paintLineXY(a.plus(ab90.times(14)), a.plus(ab90.times(14)).plus(ab1.times(abLength / 2 - textLength / 2 - 5)))
				paintLineXY(a.plus(ab90.times(14)).plus(ab1.times(abLength / 2 + textLength / 2 + 5)), a.plus(ab90.times(14)).plus(ab1.times(abLength)))
				const textCenter = a.plus(ab90.times(14)).plus(ab1.times(abLength / 2))
				gl.pushMatrix()
				gl.translate(textCenter)
				gl.scale(20, 20, 10)
				gl.multMatrix(M4.forSys(ab1, ab90.times(1), V3.Z))
				gl.translate(-textLength / 2 / 20, -0.5, 0)
				renderText(cst.distance, 0xff0000)
				gl.popMatrix()
				break
			}
			case 'pointLineDistance': {
				const texture = getTextureForString(cst.distance)
				const p = cst.point.V3(), ab = cst.other.getVectorAB(), a = cst.other.a.V3()
				const ab1 = ab.unit()
				const abLength = ab.length()
				const ab90 = ab1.getPerpendicular()
				const ap = p.minus(a)
				const apProj = ab90.times(-ap.dot(ab90))
				paintLineXY(a, a.plus(apProj))
				paintLineXY(a.plus(apProj), p)
				const textLength = texture.textWidthPx / TEXT_TEXTURE_HEIGHT * 20
				const textCenter = a.plus(apProj.times(1 / 2))
				gl.pushMatrix()
				gl.translate(textCenter)
				gl.scale(20, 20, 10)
				gl.multMatrix(M4.forSys(apProj.unit(), ab1, V3.Z))
				gl.translate(-textLength / 2 / 20, -0.5, 0)
				renderText(cst.distance, 0xff0000)
				gl.popMatrix()
				break
			}
			case 'pointPlaneDistance': {
				const texture = getTextureForString(cst.distance)
				const p = cst.point.V3()
				const sketchLineWC = sketch.plane.intersectionWithPlane(cst.other.plane)
				const sketchLineSC = sketchLineWC.transform(sketch.worldToSketchMatrix)
				const a = sketchLineSC.closestPointToPoint(p)
				const ap = p.minus(a)
				paintLineXY(a, p)
				const textLength = texture.textWidthPx / TEXT_TEXTURE_HEIGHT * 20
				const textCenter = a.plus(ap.times(1 / 2))
				gl.pushMatrix()
				gl.translate(textCenter)
				gl.scale(20, 20, 10)
				gl.multMatrix(M4.forSys(sketchLineSC.dir1.cross(V3.Z), sketchLineSC.dir1, V3.Z))
				gl.translate(-textLength / 2 / 20, -0.5, 0)
				renderText(cst.distance, 0xff0000)
				gl.popMatrix()
				break
			}
			case 'pointOnLine': {
				const curve = getCachedCurve(cst.other, sketch.plane)
				const t = curve.closestTToPoint(cst.point.V3())
				const tangentLength = curve.tangentAt(t).length()
				CURVE_PAINTERS[curve.constructor.name](curve, 0x000000, t - 12 / tangentLength, t - 4 / tangentLength, 3)
				CURVE_PAINTERS[curve.constructor.name](curve, 0x000000, t + 4 / tangentLength, t + 12 / tangentLength, 3)
				break
			}
			case 'angle':
				const first = cst.cs[0], second = cst.cs[1]
				const intersection = first.intersection(second)
				const startAngle = (first.angleAB() + (cst.f[0] == -1 ? PI : 0)) % (2 * PI), endAngle = startAngle + cst.value
				paintArc(intersection, 20, 1, 0x000000, startAngle, endAngle)
				break
			case 'equalLength': {
				for (let c = 0; c < cst.segments.length; c++) {
					const line = cst.segments[c]
					const ab = line.getVectorAB()
					const abLength = ab.length()
					const ab1 = ab.unit()
					const ab90 = ab.getPerpendicular().unit()
					const crossLength = 10
					const crossSpacing = 3
					for (let i = 0; i < crossCount; i++) {
						const s = abLength / 2 - crossSpacing * crossCount / 2 + i * crossSpacing + 10
						const crossStart = line.a.V3().plus(ab1.times(s)).plus(ab90.times(crossLength / 2))
						const crossEnd = line.a.V3().plus(ab1.times(s)).plus(ab90.times(-crossLength / 2))
						//console.log('crosspos', crossStart, crossEnd)
						paintLineXY(crossStart, crossEnd)
					}
				}
				crossCount++
				break
			}
			default:
				throw new Error('unknown cst ' + cst.type)
			// console.log('errror', cst.type);
		}

	})
}


function serializeFeatures(features: Feature[]) {
	return serialize(features)
	//return '[' + featureStack.map(f => f.toSource()).join(',') + ']'
}

function load(key: string) {
	$<HTMLInputElement>('saveNameInput').value = key
	window.location.hash = key
	console.log(localStorage.getItem(key))
	// console.log(unserialize(localStorage.getItem(key)))
	featureStack = unserialize(localStorage.getItem(key))
	selected = []
	hoverHighlight = null
	editingSketch = null
	updateSelected()
	updateFeatureDisplay()
	// rebuildModel()
	// let lastSketch = featureStack.slice().reverse().find(f => f instanceof Sketch)
	// lastSketch && modePush(MODES.SKETCH, lastSketch)
}
function initLoadSave() {
	const loadSelect = $('loadSelect') as HTMLSelectElement
	const keys = arrayFromFunction(localStorage.length, i => localStorage.key(i))
	keys.sort()
	keys.forEach(key => loadSelect.adopt(new MooEl('option', {html: key})))
	loadSelect.onchange = function () {
		const key = loadSelect.value
		if (key) {
			load(key)
		}
	}


	const saveButton = $<HTMLButtonElement>('saveButton')
	saveButton.onclick = function () {
		const key = $<HTMLInputElement>('saveNameInput').value
		if (null == key) {
			loadSelect.adopt(new MooEl('option', {html: key}))
		}
		localStorage.setItem(key, serializeFeatures(featureStack))
		console.log('saved ' + key)
	}
}

function isSketchEl(el: any): el is SketchSegment | SegmentEndPoint {
	return el instanceof SketchBezier || el instanceof SegmentEndPoint
		|| el instanceof SketchLineSeg || el instanceof SketchArc
}

//var sketchPlane = new CustomPlane(V3(0, 0,1), V3.X, V3.Y, -500, 500, -500, 500, 0xff00ff);
let editingSketch: Sketch, featureStack: Feature[] = []
function initModel() {

}
function directionObjectToVector(dirObj) {
	if (dirObj instanceof CustomPlane) {
		return dirObj.normal1
	} else {
		console.log(dirObj)
		throw new Error('uuuh' + dirObj)
	}
}
const modelColorCount = 8
const modelColors = chroma.scale(['#59f7ff', '#ff2860']).mode('lab').colors(modelColorCount).map(s => chroma(s).gl() as GL_COLOR)
function rebuildModel() {
	const DISABLE_CONSOLE = false
	DISABLE_CONSOLE && disableConsole()
	console.log('rebuilding model')
	//DEBUG = false


	namesPublishedBy = new Map()
	publishedObjects = new Map()
	function publish(name: string, object: any, feature: Feature = null) {
		if (!namesPublishedBy.has(name)) {
			namesPublishedBy.set(name, feature)
			publishedObjects.set(name, object)
		}

	}


	modelBREP = undefined
	brepMesh = undefined
	brepPoints = []
	brepEdges = []

	featureError = undefined

	planes = [
		new CustomPlane(V3.O, V3.Y, V3.Z, 'planeYZ', 0xffaaaa, -500, 500, -500, 500),
		new CustomPlane(V3.O, V3.Z, V3.X, 'planeZX', 0xaaffaa, -500, 500, -500, 500),
		new CustomPlane(V3.O, V3.X, V3.Y, 'planeXY', 0xaaaaff, -500, 500, -500, 500),
	]
	planes.forEach(customPlane => publish(customPlane.name, customPlane))

	const loopEnd = -1 == rebuildLimit ? featureStack.length : min(rebuildLimit, featureStack.length)
	for (let featureIndex = 0; featureIndex < loopEnd; featureIndex++) {
		const color = modelColors[featureIndex % modelColorCount]
		const feature = featureStack[featureIndex]
		if (!catchErrors) {
			feature.apply(publish, color)
		} else {
			try {
				feature.apply(publish, color)
			} catch (error) {
				featureError = {feature, error}
				updateFeatureDisplay()
				console.error(error)
				// throw error
				break
			}
		}
		publish(feature.name, feature, feature)
	}
	if (modelBREP) {
		brepMesh = modelBREP.toMesh()
		modelBREP.faces.forEach(face =>
			brepEdges.pushAll(face.getAllEdges().filter(edge => !edge.flippedOf || edge.id < edge.flippedOf.id)))
		// brepEdges.pushAll(isEdges)
		modelBREP.vertexNames && (brepPoints = Array.from(modelBREP.vertexNames.keys()))
		console.log(modelBREP.vertexNames)
		//brepMesh.computeNormalLines(10)
		//brepMesh.computeWireframeFromFlatTriangles()
		brepMesh.compile()
	}
	updateSelected()
	paintScreen()
	DISABLE_CONSOLE && enableConsole()
}
// DECLARATIONS
// GLOBALS
let setupCameraListener: (e: typeof eye) => void
let missingEls: any[] = [], namesPublishedBy: Map<string, Feature>, publishedObjects: Map<string, any>
let faces = []
let highlighted: any[] = [], selected: any[] = [], paintDefaultColor = 0x000000
let gl: LightGLContext
let modelBREP: B2, brepMesh: Mesh, brepPoints: V3[], planes: P3[], brepEdges: Edge[]
let rebuildLimit = -1, rollBackIndex = -1, featureError: {feature: Feature, error: any}
let drPs: V3[] = [], drVs: {anchor: V3, dir1: V3, color: int}[] = []

let eye = {pos: V(1000, 1000, 1000), focus: V3.O, up: V3.Z, zoomFactor: 1}
let hoverHighlight: any = undefined
// console.log = oldConsole;
let modeStack: any[] = []
let shaders: { [name: string]: Shader } = {}

function renderColor(mesh: Mesh, color: int, mode?: DRAW_MODES) {
	shaders.singleColor.uniforms({
		color: hexIntToGLColor(color)
	}).draw(mesh, mode)
}
function renderColorLines(mesh: Mesh, color: int) {
	shaders.singleColor.uniforms({
		color: hexIntToGLColor(color)
	}).draw(mesh, DRAW_MODES.LINES)
}
let TEXT_TEXTURE_HEIGHT = 128
function renderText(string: string, color: int) {
	const texture = getTextureForString(string)
	texture.bind(0)
	gl.pushMatrix()
	gl.scale(texture.width / TEXT_TEXTURE_HEIGHT, 1, 1)
	shaders.textureColor.uniforms({
		texture: 0,
		color: hexIntToGLColor(color)
	}).draw(meshes.text)
	gl.popMatrix()
}
function drawVector(vector: V3, anchor: V3, color = 0x0000ff, size = 1, _gl = gl) {
	_gl.pushMatrix()

	const vT = vector.getPerpendicular().unit()
	_gl.multMatrix(M4.forSys(vector, vT, vector.cross(vT).unit(), anchor))
	1 != size && _gl.scale(size, size, size)
	_gl.shaders.singleColor.uniforms({
		color: hexIntToGLColor(color)
	}).draw(_gl.meshes.vector)

	_gl.popMatrix()
}
function drawVectors(_gl = gl) {
	drawVector(V3.X, V3.O, 0xff0000, undefined, _gl)
	drawVector(V3.Y, V3.O, 0x00ff00, undefined, _gl)
	drawVector(V3.Z, V3.O, 0x0000ff, undefined, _gl)

	drVs.forEach(vi => drawVector(vi.dir1, vi.anchor, vi.color, undefined, _gl))
}

function drawPoint(p, color = 0x000000, size = 5) {
	gl.pushMatrix()
	gl.translate(p)
	gl.scale(size, size, size)
	shaders.singleColor.uniforms({color: hexIntToGLColor(color)}).draw(meshes.sphere1)
	gl.popMatrix()
}
function drawPoints(size: number = 2) {
	drPs.forEach(info => drawPoint(info.p || info, 0xcc0000, size))
	brepPoints && brepPoints.forEach(p =>
		drawPoint(p, hoverHighlight == p ? 0x0adfdf : (selected.includes(p) ? 0xff0000 : 0xcccc00), size))
}
let catchErrors = false
function handleCatchErrors(e) {
	catchErrors = e.target.value == 'on'
	rebuildModel()
}
function conicPainter(mode: 0 | 1 | 2, ellipse: SemiEllipseCurve, color: int, startT: number, endT: number, width = 2) {
	shaders.ellipse3d.uniforms({
		f1: ellipse.f1,
		f2: ellipse.f2,
		center: ellipse.center,
		color: hexIntToGLColor(color),
		startT: startT,
		endT: endT,
		scale: width,
		mode: mode
	}).draw(meshes.pipe)
}
const CURVE_PAINTERS: {[curveConstructorName: string]: (curve: Curve, color: int, startT: number, endT: number, width: number) => void} = {
	[SemiEllipseCurve.name]: conicPainter.bind(undefined, 0),
	[EllipseCurve.name]: conicPainter.bind(undefined, 0),
	[ParabolaCurve.name]: conicPainter.bind(undefined, 1),
	[HyperbolaCurve.name]: conicPainter.bind(undefined, 2),
	[ImplicitCurve.name](curve: ImplicitCurve, color, startT, endT, width = 2, normal = V3.Z) {
		let mesh = cachedMeshes.get(curve)
		if (!mesh) {
			mesh = new Mesh({triangles: true, normals: true, lines: false, colors: false})
			curve.addToMesh(mesh)
			mesh.compile()
			//mesh=Mesh.sphere(2)
			cachedMeshes.set(curve, mesh)
		}
		// TODO: draw only part
		//startT: startT,
		//	endT: endT,
		shaders.generic3d.uniforms({
			color: hexIntToGLColor(color),
			scale: width,
		}).draw(mesh)
	},
	[BezierCurve.name](curve: BezierCurve, color, startT, endT, width = 2, normal = V3.Z) {
		shaders.bezier3d.uniforms({
			p0: curve.p0,
			p1: curve.p1,
			p2: curve.p2,
			p3: curve.p3,
			color: hexIntToGLColor(color),
			startT: startT,
			endT: endT,
			scale: width,
			normal: normal
		}).draw(meshes.pipe)
	},
	[L3.name](curve: L3, color, startT, endT, width = 2, normal = V3.Z) {
		gl.pushMatrix()
		const a = curve.at(startT), b = curve.at(endT)
		const ab = b.minus(a), abT = ab.getPerpendicular().unit()
		const m = M4.forSys(ab, abT, ab.cross(abT).unit(), a)
		gl.multMatrix(m)
		gl.scale(1, width, width)
		shaders.singleColor.uniforms({
			color: hexIntToGLColor(color), // TODO: error checking
		}).draw(meshes.pipe)

		gl.popMatrix()
	},
}
CURVE_PAINTERS[PICurve.name] = CURVE_PAINTERS[ImplicitCurve.name]


function paintLineXY(a: V3, b: V3, color = paintDefaultColor, width = 2) {
	const ab = b.minus(a)
	if (ab.likeO()) { return }
	const abT = ab.getPerpendicular().toLength(width)
	//console.log(ab)
	gl.pushMatrix()
	gl.multMatrix(M4.forSys(ab, abT, V3.Z, a))
	renderColor(meshes.segment, color)
	gl.popMatrix()
}
function paintArc(center: V3, radius: number, width: number, color: int, startAngle = 0, endAngle = 2 * Math.PI) {
	// startAngle = isFinite(startAngle) ? startAngle : 0
	// endAngle = isFinite(endAngle) ? endAngle : 2 * Math.PI
	gl.pushMatrix()
	gl.translate(center)
	shaders.arc2.uniforms({
		color: hexIntToGLColor(color),
		offset: startAngle,
		step: endAngle - startAngle,
		radius: radius,
		width: width
	}).draw(meshes.segment)
	gl.popMatrix()
}

function paintBezier(ps: V3[], width: number, color: int, startT: number, endT: number) {
	// TODO PS AREN'T IN THE ORDER YOU EXPECT THEM TO BE!!!
	assertVectors.apply(undefined, ps)
	shaders.bezier.uniforms({
		color: hexIntToGLColor(color),
		width: width,
		p0: ps[0],
		p1: ps[2],
		p2: ps[3],
		p3: ps[1],
		startT: startT || 0,
		endT: endT || 1
	}).draw(meshes.segment)
}

function randomVec4Color(opacity = 1.0) {
	return [Math.random(), Math.random(), Math.random(), opacity]
}
function drawEdge(edge: Edge, color = 0x000000, width = 2) {
	CURVE_PAINTERS[edge.curve.constructor.name](edge.curve, color, edge.minT, edge.maxT, width)
}
const cachedMeshes = new Map()
function drawFace(face: Face, color: int) {
	let mesh = cachedMeshes.get(face)
	if (!mesh) {
		mesh = new Mesh({triangles: true, lines: true, normals: true})
		face.addToMesh(mesh)
		mesh.compile()
		cachedMeshes.set(face, mesh)
	}
	shaders.singleColor.uniforms({color: hexIntToGLColor(color)}).draw(mesh)
}
let meshes: any = {}
let ZERO_EL = {}
function drawEl(el: any, color: int) {
	if (el instanceof V3) {
		drawPoint(el, color)
	} else if (el instanceof Edge) {
		drawEdge(el, color)
	} else if (el instanceof Face) {
		drawFace(el, color)
	} else if (el instanceof CustomPlane) {
		drawPlane(el, color)
	} else if (el instanceof Curve) {
		drawEdge(Edge.forCurveAndTs(el), color)
	} else if (isSketchEl(el)) {
		paintLooseSeg(el, color)
	} else if (el == ZERO_EL) {} else {
		assert(false)
	}
}
let paintScreen = function () {
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)


	gl.loadIdentity()
	planes.forEach(plane => drawPlane(plane, plane.color))
	/*
	 gl.translate(0, 0, -5);
	 gl.rotate(30, 1, 0, 0);
	 gl.rotate(angle, 0, 1, 0);
	 */
	drawVectors()

	drawPoints()

	missingEls.forEach(el => drawEl(el, 0xff0000))

	gl.pushMatrix()
	if (brepMesh) {
		// paint faces
		let faceIndex = modelBREP.faces.length
		while (faceIndex--) {

			const face = modelBREP.faces[faceIndex]
			const faceTriangleIndexes = brepMesh.faceIndexes.get(face)
			const faceColor = hoverHighlight == face ? hexIntToGLColor(0xff00ff)
				: selected.includes(face) ? hexIntToGLColor(0x00ff45)
					: hoverHighlight == face.info.feature ? hexIntToGLColor(0x4500ff)
					: face.info.color
			shaders.lighting.uniforms({
				color: faceColor
			}).draw(brepMesh, DRAW_MODES.TRIANGLES, faceTriangleIndexes.start, faceTriangleIndexes.count)
			/*
			 shaders.singleColor.uniforms({
			 color: hexIntToGLColor(0x0000ff)
			 }).draw(brepMesh, DRAW_MODES.LINES);
			 */
		}

		// paint edges
		brepEdges.forEach(edge => drawEdge(edge, hoverHighlight == edge ? 0xff00ff : (selected.includes(edge) ? 0x00ff45 : COLORS.RD_STROKE), 2))

		gl.projectionMatrix.m[11] -= 1 / (1 << 22) // prevent Z-fighting
		shaders.singleColor.uniforms({
			color: hexIntToGLColor(COLORS.PP_STROKE)
		}).draw(brepMesh, DRAW_MODES.LINES)
		gl.projectionMatrix.m[11] += 1 / (1 << 22) // prevent Z-fighting
	}
	gl.popMatrix()

	const DZ = 0.1
	gl.projectionMatrix.m[11] -= DZ // prevent Z-fighting
	featureStack.filter(f => f instanceof Sketch && !f.hide).forEach(s => {
		if (-1 != rebuildLimit && rebuildLimit <= featureStack.indexOf(s)) return
		paintSketch(s)
	})
	gl.projectionMatrix.m[11] += DZ
}

function setupCamera(_eye: typeof eye, _gl: LightGLContext = gl) {
	const {pos, focus, up, zoomFactor} = _eye
	//console.log("pos", pos.$, "focus", focus.$, "up", up.$)
	_gl.matrixMode(_gl.PROJECTION)
	_gl.loadIdentity()
	//_gl.perspective(70, _gl.canvas.width / _gl.canvas.height, 0.1, 1000);
	const lr = _gl.canvas.width / 2 / zoomFactor
	const bt = _gl.canvas.height / 2 /zoomFactor
	_gl.ortho(-lr, lr, -bt, bt, -1e4, 1e4)
	_gl.lookAt(pos, focus, up)
	_gl.matrixMode(_gl.MODELVIEW)
	setupCameraListener && setupCameraListener(_eye)
}
function drawPlane(customPlane: CustomPlane, color: int, _gl: LightGLContext = gl) {
	_gl.pushMatrix()
	_gl.multMatrix(M4.forSys(customPlane.right, customPlane.up, customPlane.normal1))
	_gl.translate(customPlane.sMin, customPlane.rMin, customPlane.w)
	_gl.scale(customPlane.sMax - customPlane.sMin, customPlane.tMax - customPlane.rMin, 1)

	const shader = hoverHighlight == customPlane ? _gl.shaders.singleColorHighlight : _gl.shaders.singleColor

	shader.uniforms({color: hexIntToGLColor(color)}).draw(_gl.meshes.xyLinePlane, DRAW_MODES.LINES)

	_gl.popMatrix()
}

function segmentIntersectsRay(a, b, p, dir) {
	const ab = b.minus(a)
	const abXdir = ab.cross(dir)
	const div = abXdir.squared()
	const anchorDiff = p.minus(a)
	const s = anchorDiff.cross(dir).dot(abXdir) / div
	const t = anchorDiff.cross(ab).dot(abXdir) / div
	//console.log(segmentIntersectsRay, a, b, "ab", ab, "p", p, "dir", dir, s > 0 && t / div >= 0 && t / div <= 1, "s",
	// s, "t", t, "div", div)
	return t > 0 && s >= 0 && s <= 1
}
function template(templateName: string, map: {[key: string]: string}): HTMLElement {
	let html = $<HTMLScriptElement>(templateName).text || $<HTMLScriptElement>(templateName).innerHTML
	for (const key in map) {
		html = html.replace(new RegExp('\\$' + key, 'g'), map[key])
	}
	return $(new MooEl('div', {html: html.trim()}).firstChild)

}
function featureRollBack(feature: any, featureIndex: number) {
	rollBackIndex = featureIndex
}
let featureComp: Element
function updateFeatureDisplay() {
	const props = {features: featureStack}
	featureComp = preact.render(h(FeatureStackDisplay, props), $('featureDisplay'), featureComp)
	//const div = $('featureDisplay')
	//div.erase('text')
	//featureStack.forEach(function (feature, featureIndex) {
	//	if (rebuildLimit == featureIndex) {
	//		div.adopt(new MooEl('div', {text: 'REBUILD LIMIT', class: 'rebuildLimit'}))
	//		// div.adopt(<div class='rebuildLimit'>REBUILD LIMIT</div>)
	//	}
	//	const snAndNames = {
	//		[Extrude.name]: ['EXTR', MODES.EXTRUDE],
	//		[Rotate.name]: ['ROTA', MODES.ROTATE],
	//		[Sketch.name]: ['SKTC', MODES.SKETCH],
	//		[Pattern.name]: ['PTRN', MODES.PATTERN],
	//		[PlaneDefinition.name]: ['PLNE', MODES.PLANE_DEFINITION],
	//	}
	//	let newChild
	//	const [sn, mode] = snAndNames[feature.constructor.name]
	//	newChild = template('templateFeatureExtrude', {what: sn, name: feature.name})
	//	newChild.getElement('[name=edit]').onclick = function () {
	//		modePush(mode, feature)
	//	}
	//	newChild.inject(div)
	//	newChild.getElement('[name=delete]').onclick = function () {
	//		featureDelete(feature)
	//	}
	//	newChild.getElement('[name=rollBack]').onclick = function () {
	//		featureRollBack(feature, featureIndex)
	//	}
	//	newChild.featureLink = feature
	//	newChild.onmouseover = function (e) {
	//		const dependencies = featureDependencies(feature)
	//		const dependents = featureDependents(feature)
	//		div.getChildren().filter((subDiv: any) => dependencies.includes(subDiv.featureLink)).addClass('isDependedOn')
	//		div.getChildren().filter((subDiv: any) => dependents.includes(subDiv.featureLink)).addClass('hasDependents')
	//		hoverHighlight = this.featureLink
	//		paintScreen()
	//	}
	//	newChild.onmouseout = function (e) {
	//		div.getChildren().removeClass('isDependedOn').removeClass('hasDependents')
	//	}
	//	feature.hide && newChild.getElement('[name=toggleHide]').addClass('hidden')
	//	newChild.getElement('[name=toggleHide]').onclick = function () {
	//		feature.hide = !feature.hide
	//		this.toggleClass('hidden', feature.hide)
	//		paintScreen()
	//	}
	//})
}
function updateSelected() {
	let div = $('selectedElements')
	div.erase('text')
	selected.forEach(sel => {
		// TODO, if only necessary part of model is rebuilt, this probably wont be necessary:
		let target, name
		if (sel instanceof NameRef) {
			target = sel.get()
			name = sel.ref
		} else {
			target = sel
			name = sel.name || modelBREP && modelBREP.vertexNames && modelBREP.vertexNames.get(sel)
		}
		const newChild = template('template', {what: sel.constructor.name, name: name})
		if (!target) {
			newChild.addClass('notfound')
		}
		newChild.onmouseover = function (e) {
			if (target) {
				hoverHighlight = target
			} else {
				// lookup missing
				missingEls.push(sel.lastHit)
			}
			paintScreen()
		}
		newChild.onmouseout = function (e) {
			if (target) {
				hoverHighlight = null
			} else {
				missingEls.remove(sel.lastHit)
			}
			paintScreen()
		}
		newChild.getElement('.remove').onclick = function (e) {
			editingSketch.removeElement(target)
			updateSelected()
			paintScreen()
		}
		newChild.getElement('.unsel').onclick = function (e) {
			selected.remove(sel)
			updateSelected()
			paintScreen()
		}
		sel.toBrepEdge && newChild.grab(new MooEl('span', {
			text: sel.toBrepEdge().curve.toSource(x => round10(x, -3)),
			style: 'font-size: small;'
		}))
		sel.surface && newChild.grab(new MooEl('textarea', {
			text: sel.sce,
			style: 'font-size: xx-small;display:block;width:100%;'
		}))
		newChild.inject(div)
	})
	div = $('selectedConstraints')
	div.erase('text')
	if (MODES.SKETCH == modeGetCurrent()) {
		selected.flatMap((el) => editingSketch.getConstraintsFor(el)).unique().forEach(cst => {
			let newChild
			if ('pointDistance' == cst.type
				|| 'pointLineDistance' == cst.type
				|| 'pointPlaneDistance' == cst.type) {
				newChild = template('templateDistance', {name: cst.type, id: cst.id})
				newChild.getElement('.distanceInput').value = cst.distance
				newChild.getElement('.distanceInput').onchange = function (e) {
					cst.distance = e.target.value
					rebuildModel()
					paintScreen()
				}
			} else if ('angle' == cst.type) {
				newChild = template('templateAngle', {name: cst.type, id: cst.id})
				const input = newChild.getElement('.distanceInput')
				newChild.getElement('.distanceInput').value = round10(rad2deg(cst.value), -5)
				newChild.getElement('.fa').onclick = () => {
					cst.f[0] *= -1
					input.value = round(rad2deg(abs(cst.value = (-PI + cst.value) % (2 * PI))))
					paintScreen()
				}
				newChild.getElement('.fb').onclick = () => {
					cst.f[1] *= -1
					input.value = round(rad2deg(abs(cst.value = (-PI + cst.value) % (2 * PI))))
					paintScreen()
				}
				newChild.getElement('.fv').onclick = () => {
					input.value = round(rad2deg(abs(cst.value -= sign(cst.value) * 2 * PI)))
					paintScreen()
				}
				input.onchange = function (e) {
					cst.value = e.target.value * DEG
					rebuildModel()
					paintScreen()
				}
			} else {
				newChild = template('templateConstraint', {name: cst.type})

				cst.cs.forEach(function (el) {
					const subChild = template('templateConstraintSub', {what: el.constructor.name, name: el.name})
					subChild.inject(newChild)
					subChild.getElement('.removeFromConstraint').onclick = function (e) {
						removeFromConstraint(el, editingSketch, cst)
						updateSelected()
						paintScreen()
					}
					subChild.onmouseover = function (e) {
						hoverHighlight = el
						e.stopPropagation()
						paintScreen()
					}
					subChild.onmouseout = function (e) {
						hoverHighlight = null
						paintScreen()
					}
				})
			}
			newChild.getElement('.remove').onclick = function (e) {
                editingSketch.deleteConstraint(cst)
				updateSelected()
			}
			newChild.onmouseover = function (e) {
				hoverHighlight = cst
				paintScreen()
			}
			newChild.onmouseout = function (el) {
				hoverHighlight = null
				paintScreen()
			}
			newChild.inject(div)
		})
	}
}
class NameRef {
	static pool: Map<string, NameRef> = new Map()

	private constructor(public ref: string, public lastHit: any) {}

	get() {
		const hit = publishedObjects.get(this.ref)
		// let hit = planes.find(plane => plane.name == this.ref)
		// 	|| modelBREP && modelBREP.faces.find(face => face.name == this.ref)
		// 	|| modelBREP && modelBREP.faces.firstMatch(face => face.getAllEdges().find(edge => edge.name ==
		// this.ref)) || modelBREP && modelBREP.vertexNames && mapReverse(modelBREP.vertexNames, this.ref)
		if (hit) this.lastHit = hit
		return hit
	}

	getOrThrow() {
		const hit = this.get()
		if (!hit) {
			throw new Error(`could not find ${this.lastHit.constructor.name} ${this.ref}`)
		}
		return hit
	}

	getPlane() {
		return planes.find(plane => plane.name == this.ref) || modelBREP && modelBREP.faces.find(face => face.name == this.ref).plane
	}

	get plane() {
		return this.getPlane()
	}

	get name() {
		const hit = this.get()
		return hit && hit.name
	}

	static forRef(ref: string, lastHit?: any) {
        const x = NameRef.pool.get(ref)
        if (x) return x
        return new NameRef(ref, lastHit)
    }
	static forObject(o: any) {
		if (o instanceof V3) {
			return NameRef.forRef(modelBREP.vertexNames.get(o), o)
		} else {
			assert(o.name)
			return NameRef.forRef(o.name, o)
		}
	}

	static UNASSIGNED = new NameRef('UNASSIGNED', ZERO_EL)
}
NameRef.UNASSIGNED.get = function () { return undefined }
NameRef.UNASSIGNED.getOrThrow = function () { throw new Error('this nameref has never been assigned a value') }



function mapReverse(map, value) {
	const it = map.keys()
	let key
	while (key = it.next().value) {
		if (map.get(key) == value) {
			return key
		}
	}
}
function initMeshes(_meshes) {
	_meshes.sphere1 = Mesh.sphere(2)
	_meshes.segment = Mesh.plane({startY: -0.5, height: 1, detailX: 128})
	_meshes.text = Mesh.plane()
	_meshes.vector = Mesh.rotation([V3.O, V(0, 0.05, 0), V(0.8, 0.05), V(0.8, 0.1), V(1, 0)], L3.X, TAU, 16, true)
	_meshes.pipe = Mesh.rotation(arrayFromFunction(128, i => new V3(i / 127, -0.5, 0)), L3.X, TAU, 8, true)
	_meshes.xyLinePlane = Mesh.plane()
}
function initShaders(_shaders) {
	_shaders.singleColor = new Shader(vertexShaderBasic, fragmentShaderColor)
	_shaders.multiColor = new Shader(vertexShaderColor, fragmentShaderVaryingColor)
	_shaders.singleColorHighlight = new Shader(vertexShaderBasic, fragmentShaderColorHighlight)
	_shaders.textureColor = new Shader(vertexShaderTexture, fragmentShaderTextureColor)
	_shaders.arc = new Shader(vertexShaderRing, fragmentShaderColor)
	_shaders.arc2 = new Shader(vertexShaderArc, fragmentShaderColor)
	_shaders.ellipse3d = new Shader(vertexShaderConic3d, fragmentShaderColor)
	_shaders.generic3d = new Shader(vertexShaderGeneric, fragmentShaderColor)
	_shaders.bezier3d = new Shader(vertexShaderBezier3d, fragmentShaderColor)
	_shaders.bezier = new Shader(vertexShaderBezier, fragmentShaderColor)
	_shaders.lighting = new Shader(vertexShaderLighting, fragmentShaderLighting)
	_shaders.waves = new Shader(vertexShaderWaves, fragmentShaderLighting)
}
function initOtherEvents() {
	window.onkeypress = function (e) {
		if ('Delete' == e.key) {
			selected.forEach(x => editingSketch.removeElement(x))
			paintScreen()
		}
	}
	gl.canvas.addEventListener('mouseup', function (e) {
		// don't put modeGetCurrent().mouseup in local var as 'this' won't be bound correctly
		modeGetCurrent().mouseup && modeGetCurrent().mouseup(e, getMouseLine(e))
		paintScreen()
	})
	gl.canvas.addEventListener('mousedown', function (e) {
		if (1 == e.button) {
			modePop()
		}
		if (0 == e.button) {
			const mouseLine = getMouseLine(e)
			modeGetCurrent().mousedown(e, mouseLine)
		}

		return false
	})
	$('clearmode').addEvent('click', modePop)

}
function initPointInfoEvents(_gl = gl) {
	_gl.canvas.addEventListener('mousemove', function (e) {
		const mouseLine = getMouseLine({x: e.clientX, y: e.clientY})
		modeGetCurrent().mousemove && modeGetCurrent().mousemove(e, mouseLine)
		{
			let pp, html = '', closestP = Infinity
			drPs.forEach(info => {
				const p = info.p || info
				let text = p.toString(x => round10(x, -4)) + (info.p ? ' ' + info.text : '')
				text = `<li>${text}</li>`
				const dist = mouseLine.distanceToPoint(p)
				if (dist < 16) {
					if (pp && p.distanceTo(pp) < 10) {
						html += text
					} else if (dist < closestP) {
						pp = p
						closestP = dist
						html = text
					}
				}
			})
			if (pp) {
				const pSC = _gl.projectionMatrix.times(_gl.modelViewMatrix).transformPoint(pp)
				const x = (pSC.x * 0.5 + 0.5) * window.innerWidth, y = (-pSC.y * 0.5 + 0.5) * window.innerHeight
				tooltipShow(html, x, y)
			} else {
				tooltipHide()
			}
		}
		paintScreen()
	})
}
function initNavigationEvents(_gl = gl, eye: {pos: V3, focus: V3, up: V3, zoomFactor: number}, paintScreen: () => void) {
	const canvas: HTMLCanvasElement = $(_gl.canvas)
	let lastPos: V3 = V3.O
	//_gl.onmousedown.push((e) => {
	//	e.preventDefault()
	//	e.stopPropagation()
	//})
	//_gl.onmouseup.push((e) => {
	//	e.preventDefault()
	//	e.stopPropagation()
	//})
	canvas.addEventListener('mousemove', e => {
		const pagePos = V(e.pageX, e.pageY)
		const delta = lastPos.to(pagePos)
		//noinspection JSBitwiseOperatorUsage
		if (e.buttons & 4) {
			// pan
			const moveCamera = V(-delta.x * 2 / _gl.canvas.width, delta.y * 2 / _gl.canvas.height)
			const inverseProjectionMatrix = _gl.projectionMatrix.inversed()
			const worldMoveCamera = inverseProjectionMatrix.transformVector(moveCamera)
			eye.pos = eye.pos.plus(worldMoveCamera)
			eye.focus = eye.focus.plus(worldMoveCamera)
			setupCamera(eye, _gl)
			paintScreen()
		}
		// scene rotation
		//noinspection JSBitwiseOperatorUsage
		if (e.buttons & 2) {
			const rotateLR = -delta.x / 6.0 * DEG
			const rotateUD = -delta.y / 6.0 * DEG
			// rotate
			let matrix = M4.rotateLine(eye.focus, eye.up, rotateLR)
			//let horizontalRotationAxis = focus.minus(pos).cross(up)
			const horizontalRotationAxis = eye.up.cross(eye.pos.minus(eye.focus))
			matrix = matrix.times(M4.rotateLine(eye.focus, horizontalRotationAxis, rotateUD))
			eye.pos = matrix.transformPoint(eye.pos)
			eye.up = matrix.transformVector(eye.up)

			setupCamera(eye, _gl)
			paintScreen()
		}
		lastPos = pagePos
	})
	canvas.addEvent('mousewheel', function (e) {
		//console.log(e)
		eye.zoomFactor *= pow(0.9, -e.wheel)
		const targetPos = e.target.getPosition()
		const mouseCoords = {x: e.page.x - targetPos.x, y: e.page.y - targetPos.y}
		const moveCamera = V(mouseCoords.x * 2 / _gl.canvas.width - 1, -mouseCoords.y * 2 / _gl.canvas.height + 1, 0).times(1 - 1 / pow(0.9, -e.wheel))
		const inverseProjectionMatrix = _gl.projectionMatrix.inversed()
		const worldMoveCamera = inverseProjectionMatrix.transformVector(moveCamera)
		//console.log("moveCamera", moveCamera)
		//console.log("worldMoveCamera", worldMoveCamera)
		eye.pos = eye.pos.plus(worldMoveCamera)
		eye.focus = eye.focus.plus(worldMoveCamera)
		setupCamera(eye, _gl)
		paintScreen()
		e.preventDefault()
	})
}
//const sFace = B2T.rotateEdges([Edge.forCurveAndTs(EllipseCurve.XY, 0, 90 * DEG).rotateX(90 * DEG),StraightEdge.throughPoints(V3.Z, V3.X)], 45 * DEG, 'blah').faces.find(face => face.surface instanceof EllipsoidSurface)
//const face2 = B2T.extrudeEdges([Edge.forCurveAndTs(EllipseCurve.forAB(1, -1), -PI, 0), StraightEdge.throughPoints(V3.X, V3.X.negated())], P3.XY.flipped(), V3.Z, 'cyl')
//const cylface = cyl.faces.find(face => face instanceof RotationFace)//.rotateX(50 * DEG)
//assert(cylface.surface.facesOutwards())
//const cyl = B2T.extrudeEdges([Edge.forCurveAndTs(EllipseCurve.forAB(1, -1), -PI, 0), StraightEdge.throughPoints(V3.X, V3.X.negated())], P3.XY.flipped(), V3.Z, 'cyl')
function tooltipShowEl(htmlContent: string, target: HTMLElement, side: string) {
	const tipWrap = $('tip-wrap') as HTMLDivElement
	const arrow = tipWrap.getElement('#tip-arrow') as HTMLDivElement
	assert(tipWrap)
	const tip = $('tip') as HTMLDivElement
	tipWrap.setStyle('visibility', 'visible')
	tip.set('html', htmlContent)
	const targetSize = target.getSize()
	const tipSize = tipWrap.getSize()
	const targetPos = target.getPosition()
	const arrowSize = 8
	side = side || 'bottom'
	const color = tip.getStyle('background-color')
	const borderColor = ['top', 'right', 'bottom', 'left']
		.map(s => s == side ? color : 'transparent').join(' ')
	if ('left' == side || 'right' == side) {
		const tipTop = targetPos.y + tipSize.y <=window.innerHeight
			? targetPos.y
			: targetPos.y + targetSize.y - tipSize.y
		//const borderColor = 'red blue green purple'
		const arrowRight = 'left' == side ? '' + -arrowSize * 2 + 'px' : 'auto'
		const arrowLeft = 'right' == side ? '' + -arrowSize * 2 + 'px' : 'auto'
		arrow.setStyles({
			top: '' + (targetPos.y - tipTop + targetSize.y / 2 - arrowSize) + 'px',
			bottom: 'auto',
			left: arrowLeft,
			right: arrowRight,
			borderColor,
			borderWidth: arrowSize + 'px'
		})
		const tipLeft = 'left' == side
			? targetPos.x - arrowSize - tipSize.x
			: targetPos.x + targetSize.x + arrowSize
		tipWrap.setStyle('left', tipLeft + 'px')
		tipWrap.setStyle('top', tipTop + 'px')
	} else {
		const tipLeft = targetPos.x + tipSize.x <= window.innerWidth
			? targetPos.x
			: targetPos.x + targetSize.x - tipSize.x
		const tipTop = 'top' == side
			? targetPos.y - arrowSize - tipSize.y
			: targetPos.y + targetSize.y + arrowSize
		const arrowTop = 'bottom' == side ? '' + -arrowSize * 2 + 'px' : 'auto'
		const arrowBottom = 'top' == side ? '' + -arrowSize * 2 + 'px' : 'auto'
		// place below
		$('tip-arrow').setStyles({
			top: arrowTop,
			bottom: arrowBottom,
			left: '' + (targetPos.x - tipLeft + targetSize.x / 2 - arrowSize) + 'px',
			right: 'auto',
			borderColor,
			borderWidth: arrowSize + 'px'
		})
		tipWrap.setStyle('left', tipLeft + 'px')
		tipWrap.setStyle('top', tipTop + 'px')
	}
}
function initToolTips() {
	const DELAY_DEFAULT = 1000
	let timeout: number, target: HTMLElement

    window.addEventListener('mouseover', function (e) {
		let el = e.target as HTMLElement, tt
		while (!(tt = el.dataset.tooltip || el.dataset.tooltipsrc) && el != document.body) {
			el = el.getParent()  as HTMLElement
		}
		if (tt) {
			target = el
			timeout = setTimeout(() => {
				const html = el.dataset.tooltip || $(el.dataset.tooltipsrc).innerHTML
				tooltipShowEl(html, el, el.dataset.tooltipside)
			}, DELAY_DEFAULT)
		}
	}, true)
	window.addEventListener('mouseout', function (e) {
		if (e.target == target) {
			clearTimeout(timeout)
			tooltipHide()
			target = undefined
		}
	})
	$('tip-wrap').addEvent('mouseover', tooltipHide)
}
window.onload = function () {
    //const rt = RangeTree.fromArray([-2, -1, 0,2,3,4,7-1,7])
    //rt.addIntervals([
    //    {left: -1, right: 6},
    //    {left: -2, right: 3},
    //    {left: 0, right: 4},
    //    {left: 2, right: 7},
    //    {left: 3, right: 4},
    //    {left: -2, right: -1}])
    //console.log(rt.str)
    //
    //return
	modePush(MODES.DEFAULT)

	$$('.sketchControl').set('disabled', true)

	gl = create({canvas: document.getElementById('mainCanvas') as HTMLCanvasElement})
	gl.fullscreen()
	gl.canvas.oncontextmenu = () => false

	setupCamera(eye, gl)
	//gl.cullFace(gl.FRONT_AND_BACK);
	gl.clearColor(1.0, 1.0, 1.0, 0.0)
	gl.enable(gl.BLEND)
	gl.enable(gl.DEPTH_TEST)
	// gl.enable(gl.CULL_FACE)
	gl.depthFunc(gl.LEQUAL)
	gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA) // TODO ?!


	initMeshes(gl.meshes = meshes)
	initShaders(gl.shaders = shaders)
	const b = editingSketch && editingSketch.elements.find(el => el instanceof SketchBezier).toBrepEdge().curve
	//	console.log(mesh.vertices)


	initNavigationEvents(gl, eye, paintScreen)
	initOtherEvents()
	initPointInfoEvents()
	initToolTips()


	initLoadSave()

	if (window.location.hash && window.location.hash.length > 1) {
		const key = window.location.hash.substr(1)
		console.log(key)
		load(key)
	} else {
		//initModel()
	}

	const feature = new Rotate()
	rebuildModel() // necessary to init planes
	updateFeatureDisplay()
	rebuildModel() // so warning will show
	paintScreen()

	const lastInterst = featureStack.last()
	lastInterst && modePush(snAndNames[lastInterst.constructor.name][1], lastInterst)

}

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
			checkEl(plane, plane.distanceTo(mouseLine,mindist))
		}
	}
	if (consider.includes('sketchElements')) {
        0;(featureStack.filter(f => f instanceof Sketch) as Sketch[]).forEach(sketch => {
			if (!sketch.hide && sketch.plane && sketch.plane.normal1.dot(mouseLine.dir1) < -0.1) {
				// sketch plane is facing user; ensures there is an intersection

				const mouseLineIS = mouseLine.intersectionWithPlane(sketch.plane)
				const sketchCoords = sketch.worldToSketchMatrix.transformPoint(mouseLineIS)

				const closestElement = sketch.elements.concat(getAllPoints(sketch).map(p => p.canon()).unique())
					.withMax(el => {
						let d = el.distanceToCoords(sketchCoords)
						if (el instanceof SegmentEndPoint) d -= 8
						return -d
					})

				if (closestElement && closestElement.distanceToCoords(sketchCoords) < mindist) {
					// subtract 0.001 so that sketch elements have priority over things in same plane
					checkEl(closestElement, mouseLine.pointT(mouseLineIS) - 0.1)
				}
			}
		})
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
/**
 * Transforms mouse positions on the screen into a line in world coordinates.
 */
function getMouseLine(pos: {x: number, y: number}): L3 {
	const ndc1 = V(pos.x * 2 / gl.canvas.width - 1, -pos.y * 2 / gl.canvas.height + 1, 0)
	const ndc2 = V(pos.x * 2 / gl.canvas.width - 1, -pos.y * 2 / gl.canvas.height + 1, 1)
	//console.log(ndc)
	const inverseProjectionMatrix = gl.projectionMatrix.inversed()
	const s = inverseProjectionMatrix.transformPoint(ndc1)
	const dir = inverseProjectionMatrix.transformPoint(ndc2).minus(s)
	return L3.anchorDirection(s, dir)
}


function featureDelete(feature: Feature) {
	const correspondingMode = modeStack.find(mode => mode.feature && mode.feature == feature)
	correspondingMode && modeEnd(correspondingMode)

	const featureIndex = featureStack.indexOf(feature)
	featureStack.remove(feature)
	if (-1 != featureIndex && featureIndex < rebuildLimit) {
		rebuildLimit--
	}
	updateFeatureDisplay()
	rebuildModel()
}
function featureDependents(feature: Feature) {
	const featureIndex = featureStack.indexOf(feature)
	return featureStack
		.slice(featureIndex + 1)
		.filter(followingFeature => featureDependencies(followingFeature).includes(feature))
}
function featureDependencies(feature: Feature) {
	const result = new Set()
	const featureNameDependecies = feature.dependentOnNames()
	featureNameDependecies.forEach(nameRef => {
		const dependentOnFeature = namesPublishedBy.get(nameRef.ref)
		dependentOnFeature && result.add(dependentOnFeature)
	})
	return Array.from(result.values())
}
function chainComparison(diff1: number, diff2: number) {
	return diff1 != 0 ? diff1 : diff2
}
function pointsFirst(a, b) {
	//console.log((a instanceof SegmentEndPoint ? -1 : 1) - (b instanceof SegmentEndPoint ? -1 : 1));
	return (a instanceof SegmentEndPoint ? -1 : 1) - (b instanceof SegmentEndPoint ? -1 : 1)
}

declare var saveAs: (blob: Blob, name: string, no_auto_bom?: boolean) => void
function exportModel() {
	saveAs(brepMesh.toBinarySTL(), 'export.stl')
}
function clicky() {
	for (let j = 0; j < 1; j++) {
		editingSketch.gaussNewtonStep()
	}
	editingSketch.reverse()
	paintScreen()
}
// colinear, equal length or parallel
/*
 function makeGroup2(type) {
 var els = selected.filter((el) => el instanceof SketchLineSeg || ('equalLength' != type && el instanceof CustomPlane) );
 console.log("makeGroup, els", els, type);
 if (els.length >= 2) {
 var newConstraint = {"fixed": null, "type": type, "which": []};
 var oldConstraints = [];
 for (var i = 0; i < els.length; i++) {
 var el = els[i];
 var segConstraint = getGroupConstraint(el, editingSketch, type);
 if (segConstraint) {
 oldConstraints.push(segConstraint);
 if (newConstraint.fixed && segConstraint.fixed) {
 throw new Error("cannot have two fixed");
 }
 if (segConstraint.fixed) {
 newConstraint.fixed = segConstraint.fixed;
 }
 if (segConstraint != groupConstraint) {
 segConstraint.cs.forEach(function (segConstraintOther) {
 groupConstraint.cs.push(segConstraintOther);
 });
 editingSketch.constraints.remove(segConstraint);
 }
 } else {
 if (el instanceof CustomPlane) {
 if (null != groupConstraint.fixed) {
 throw new Error("cannot have two fixed");
 }
 newConstraint.fixed = el;
 } else {
 newConstraint.cs.push(el);
 }
 }
 }
 constraints.cs.push(el);
 constraints.cs.removeAll(oldConstraints);
 rebuildModel();
 paintScreen();
 }
 console.log("sketch.constraints", editingSketch.constraints);
 }
 */
function makeGroup(type) {
	const isFixed = x => x instanceof Face || x instanceof Edge || x instanceof CustomPlane
	const els = selected.filter((el) => el instanceof SketchLineSeg || ('equalLength' != type && isFixed(el)))
	if (els.length < 2) return

	const newGroup = els.flatMap(el => {
		if (isFixed(el)) el = NameRef.forObject(el)
		const c = editingSketch.getGroupConstraint(el, type)
		return c ? c.cs : el
	}).unique()
	const fixeds = newGroup.filter(el => el instanceof NameRef)
	if (1 < fixeds.length) {
		throw new Error('cannot have two fixed')
	}

	const oldConstraints = newGroup.map(el => editingSketch.getGroupConstraint(el, type))
	editingSketch.constraints.removeAll(oldConstraints)

	// add new constraint
	// fixeds[0] may be null
	const segments = newGroup.filter(x => !(x instanceof NameRef))
	const f = fixeds[0]
	editingSketch.constraints.push(new Constraint(type, f ? segments.concat(f) : segments, {
		fixed: f,
		segments: segments
	}))

	rebuildModel()
	paintScreen()
	updateSelected()
}


function makeAngle() {
	const selSegments = selected.filter((el) => el instanceof SketchLineSeg || el.plane)
	console.log(selSegments)
	if (selSegments.length == 2) {
		const segment = selSegments[0] instanceof SketchLineSeg ? selSegments[0] : selSegments[1]
		const other = selSegments[0] instanceof SketchLineSeg ? selSegments[1] : selSegments[0]
		if (!(segment instanceof SketchLineSeg)) {
			throw new Error('at least one must be a segment')
		}
		editingSketch.constraints.push(new Constraint('angle', [segment, other], {
			segment: segment,
			other: other,
			f: [1, 1],
			value: selSegments[0].angleTo(selSegments[1])
		})) //
		rebuildModel()
		paintScreen()
		updateSelected()
	}
}
function makePerpendicular() {
	const selSegments = selected.filter((el) => el instanceof SketchLineSeg || el.plane)
	if (selSegments.length == 2) {
		const segment = selSegments[0] instanceof SketchLineSeg ? selSegments[0] : selSegments[1]
		const other = selSegments[0] instanceof SketchLineSeg ? selSegments[1] : selSegments[0]
		if (!(segment instanceof SketchLineSeg)) {
			throw new Error('at least one must be a segment')
		}
		editingSketch.constraints.push(new Constraint('perpendicular', [segment, other], {
			segment: segment,
			other: other
		}))
		rebuildModel()
		paintScreen()
		updateSelected()
	}
}
function makeDistance() {
	if (2 != selected.length) return
	const point = selected[0] instanceof SegmentEndPoint ? selected[0] : selected[1]
	let other = selected[0] instanceof SegmentEndPoint ? selected[1] : selected[0]
	if (point instanceof SegmentEndPoint && (other instanceof SketchLineSeg || other instanceof SegmentEndPoint || other.plane)) {
		let newConstraint
		if (other instanceof SegmentEndPoint) {
			newConstraint = new Constraint('pointDistance', selected.slice(), {distance: round(other.distanceToCoords(point))})
		} else if (other instanceof SketchLineSeg) {
			const distance = round(other.distanceToCoords(point))
			newConstraint = new Constraint('pointLineDistance', [point, other],
				{point: point, other: other, distance: distance})
		} else {
			const distance = round(other.plane.intersectionWithPlane(editingSketch.plane)
				.transform(editingSketch.worldToSketchMatrix).distanceToPoint(point.V3()))
			other = NameRef.forObject(other)
			newConstraint = new Constraint('pointPlaneDistance', [point, other],
				{point: point, other: other, distance: distance})
			console.log(newConstraint)
		}
		editingSketch.constraints.push(newConstraint)
		rebuildModel()
		paintScreen()
		updateSelected()
		;
		($('distanceInput' + newConstraint.id) as HTMLInputElement).select()
	}
}
function selPointOnLine() {
	if (2 != selected.length) return
	const point = selected[0] instanceof SegmentEndPoint ? selected[0] : selected[1]
	let other = selected[0] instanceof SegmentEndPoint ? selected[1] : selected[0]
	if (point instanceof SegmentEndPoint && (other instanceof SketchLineSeg || other instanceof SketchBezier || other instanceof SketchArc
		|| other instanceof Face || other instanceof Edge || other instanceof CustomPlane)) {
		other = isSketchEl(other) ? other : NameRef.forObject(other)
		const newConstraint = new Constraint('pointOnLine', [point, other], {point: point, other: other})
		console.log(newConstraint)
		editingSketch.constraints.push(newConstraint)
		rebuildModel()
		paintScreen()
		updateSelected()
	}
}
function makeSelCoincident() {
	const selPoints = selected.filter(function (s) { return s instanceof SegmentEndPoint })
	if (selPoints.length >= 2) {
		for (let i = 1; i < selPoints.length; i++) {
            editingSketch.makeCoincident(selPoints[0], selPoints[i])
        }
		rebuildModel()
		paintScreen()
		updateSelected()
	}
}

function faceSketchPlane() {
	const plane = editingSketch.plane
	if (!plane) return
	const viewLine = L3.throughPoints(eye.pos, eye.focus)
	eye.focus = plane.intersectionWithLine(viewLine) || viewLine.closestPointToPoint(V3.O)
	eye.pos = eye.focus.plus(plane.normal1.times(100))
	eye.up = plane.up
	setupCamera(eye, gl)
	paintScreen()
}
function editable(str: string) {
	return function (target: any, name: string) {
		console.log(target)
	}
}
class Extrude extends Feature {

	readonly type: string = 'extrude'

	constructor(public name: string = 'extrude' + (globalId++),
	            private _segmentName: NameRef = NameRef.UNASSIGNED,
	            public start: number = 0,
	            public end: number = 100,
	            public operation: 'minus' | 'plus' = 'minus') {
		super()
		assertInst(NameRef, _segmentName)
	}

	set segmentName(sn: NameRef) {
		assertInst(NameRef, sn)
		this._segmentName = sn
	}

	get segmentName() {
		assertInst(NameRef, this._segmentName)
		return this._segmentName
	}

	dependentOnNames(): NameRef[] {
		return []
	}
	getB2(m4: M4, color?: colorstr, genFeature?: Feature) {
		const feature = this
		const loopSegment = feature.segmentName.getOrThrow()
		const loopSketch = loopSegment.sketch
		// opposite dir to plane normal1:
		const startOffset = loopSketch.plane.normal1.times(-min(feature.start, feature.end))
		const lengthOffset = m4.transformVector(loopSketch.plane.normal1.times(-Math.abs(feature.start - feature.end)))
		const startMatrix = m4.times(M4.translate(startOffset).times(loopSketch.sketchToWorldMatrix))
		let edgeLoop = loopSketch.getLoopForSegment(loopSegment)
		edgeLoop = edgeLoop.map(edge => edge.transform(startMatrix, ''))
		// TODO: test for self-intersection of edgeloop
		if (!new PlaneSurface(loopSketch.plane, loopSketch.right, loopSketch.up).edgeLoopCCW(edgeLoop)) {
			edgeLoop = edgeLoop.map(edge => edge.flipped()).reverse()
		}
		//console.log(polygonPoints.map(v =>v.$))
		const length = feature.end - feature.start
		//console.log("polypoints", polygonPoints, polygonPoints.toSource(),
		// loopSketch.plane.translated().toSource(), offsetDir.times(feature.start))
		type FI = { feature: Feature, color: number[] }
		const featureFaceInfo = { feature: genFeature, color }
		const infoFactory = FaceInfoFactory.makeStatic({ feature: genFeature, color })
		const brep = B2T.extrudeEdges(edgeLoop, loopSketch.plane.transform(startMatrix),
			lengthOffset, feature.name, undefined, infoFactory)
		return brep
	}

	apply(publish, color, m4: M4 = M4.IDENTITY, genFeature: Feature = this) {
		const feature = this
		const brep = this.getB2(m4, color, genFeature)
		brep.assertSanity()

		if (modelBREP) {
			//isEdges = modelBREP.getIntersectionEdges(brep)
			//drVs = isEdges.map(e => ({anchor: e.a, dir: e.curve.tangentAt(e.aT).unit()}))
			type FI = { feature: Feature, color: number[] }
			modelBREP = modelBREP[feature.operation](brep, new class extends FaceInfoFactory<FI> {
				constructor() {
					super()
				}

				newSubFace(original: Face, surface: Surface, contour: Edge[], holes: Edge[][] = []): FI {
					return original.info
				}
			})
			//modelBREP = B2.join([modelBREP, brep])
		} else {
			modelBREP = brep
		}

		modelBREP.faces.forEach(face => publish(face.name, face))
		modelBREP.faces.map(face => face.getAllEdges().filter(edge => !edge.flippedOf || edge.id < edge.flippedOf.id).forEach(edge => publish(edge.name, edge)))
		console.log('modelBREP.vertexNames', modelBREP.vertexNames)
		modelBREP.vertexNames && modelBREP.vertexNames.forEach((name, p) => publish(name, p, feature))
	}
}
class PlaneDefinition extends Feature {
	/*
	 if ("face" == feature.planeType && feature.faceName) {
	 var face = modelBREP.faces.find(face => face.name == feature.faceName)
	 var plane = face.surface.plane, right = plane.normal1.getPerpendicular().unit(), up = plane.normal1.cross(right)

	 var cp
	 planes.push(cp = new CustomPlane(plane.anchor.plus(plane.normal1.times(feature.offset)),
	 right, up, -500, 500, -500, 500, 0xFF4500, feature.planeId))
	 }
	 if ("immediate" == feature.planeType) {
	 var plane = eval(feature.source), right = plane.normal1.getPerpendicular().unit(), up = plane.normal1.cross(right)

	 planes.push(new CustomPlane(plane.anchor,
	 right, up, -500, 500, -500, 500, 0xFF4500, feature.planeId))
	 }
	 */
	planeId = 'plane' + globalId++
	name = 'plane' + globalId++
	planeType = 'dynamic'
	source = ''
	offset = 0
	angle = 0
	whats = []
	flipped = false

	dependentOnNames() {
		return this.whats
	}

	apply(publish, color, m4: M4 = M4.IDENTITY, genFeature: Feature = this) {
		const feature = this
		const sel = feature.whats.map(w => w.getOrThrow())
		let plane = MODES.PLANE_DEFINITION.magic(sel, feature.angle * DEG)
		let cp
		if (plane) {
			if (feature.flipped) plane = plane.flipped()

			const right = plane.normal1.getPerpendicular().unit(), up = plane.normal1.cross(right)
			planes.push(cp = new CustomPlane(plane.anchor.plus(plane.normal1.times(feature.offset)),
				right, up, feature.planeId, 0xFF4500, -500, 500, -500, 500))
		}
		if ('immediate' == feature.planeType) {
			const plane = eval(feature.source), right = plane.normal1.getPerpendicular().unit(),
				up = plane.normal1.cross(right)

			planes.push(new CustomPlane(plane.anchor,
				right, up, feature.planeId, 0xFF4500, -500, 500, -500, 500))
		}
		publish(feature.name, cp)
	}
}

class Rotate extends Feature {

	readonly type: string = 'extrude'

	constructor(public name: string = 'rotate' + (globalId++),
	            private _segmentName: NameRef = NameRef.UNASSIGNED,
	            public axis: any[] = [],
	            public start: raddd = 0,
	            public end: raddd = TAU,
	            public operation: 'minus' | 'plus' = 'minus') {
		super()
		assertInst(NameRef, _segmentName)
	}

	set segmentName(sn: NameRef) {
		assertInst(NameRef, sn)
		this._segmentName = sn
	}

	get segmentName() {
		assertInst(NameRef, this._segmentName)
		return this._segmentName
	}

	dependentOnNames(): NameRef[] {
		return []
	}

	getB2(m4: M4, color?: colorstr, genFeature?: Feature) {
		const feature = this
		const loopSegment = feature.segmentName.getOrThrow()
		const loopSketch = loopSegment.sketch
		const axis = MODES.SELECT_LINE.magic(feature.axis.map(w => w.get()))[0]
		assert(loopSketch.plane.containsLine(axis))
		const m1 = M4.forSys(axis.dir1.cross(loopSketch.plane.normal1).negated(), loopSketch.plane.normal1, axis.dir1, axis.anchor)
		const m1i = m1.inversed()
		// opposite dir to plane normal1:
		const startMatrix = M4.multiplyMultiple(m4, m1, M4.rotateZ(feature.start))
		const edgeLoopXY = loopSketch.getLoopForSegment(loopSegment)
		let edgeLoopXZ = edgeLoopXY.map(edge => edge.transform(m1i, ''))
		const rads = min(TAU, abs(feature.end - feature.start))
		// TODO: test for self-intersection of edgeloop
		if (new PlaneSurface(P3.ZX, loopSketch.right, loopSketch.up).edgeLoopCCW(edgeLoopXZ)) {
			edgeLoopXZ = Edge.reversePath(edgeLoopXZ)
		}
		//console.log(polygonPoints.map(v =>v.$))
		const length = feature.end - feature.start
		//console.log("polypoints", polygonPoints, polygonPoints.toSource(),
		// loopSketch.plane.translated().toSource(), offsetDir.times(feature.start))
		const infoFactory = FaceInfoFactory.makeStatic({ feature: genFeature, color })
		const brepXZ = B2T.rotateEdges(edgeLoopXZ, rads, feature.name, undefined, infoFactory)
		brepXZ.assertSanity()
		const brep = brepXZ.transform(startMatrix)
		return brep
	}

	apply(publish, color, m4: M4 = M4.IDENTITY, genFeature: Feature = this) {
		const feature = this
		const brep = this.getB2(m4, color, genFeature)
		brep.assertSanity()

		if (modelBREP) {
			//isEdges = modelBREP.getIntersectionEdges(brep)
			//drVs = isEdges.map(e => ({anchor: e.a, dir: e.curve.tangentAt(e.aT).unit()}))
			type FI = { feature: Feature, color: number[] }
			modelBREP = modelBREP[feature.operation](brep, new class extends FaceInfoFactory<FI> {
				constructor() {
					super()
				}

				newSubFace(original: Face, surface: Surface, contour: Edge[], holes: Edge[][] = []): FI {
					return original.info
				}
			})
			//modelBREP = B2.join([modelBREP, brep])
		} else {
			modelBREP = brep
		}

		modelBREP.faces.forEach(face => publish(face.name, face))
		modelBREP.faces.map(face => face.getAllEdges().filter(edge => !edge.flippedOf || edge.id < edge.flippedOf.id).forEach(edge => publish(edge.name, edge)))
		console.log('modelBREP.vertexNames', modelBREP.vertexNames)
		modelBREP.vertexNames && modelBREP.vertexNames.forEach((name, p) => publish(name, p, feature))
	}
}


type PatternDimensionType = {
	direction: V3,
	directionFlipped: boolean,
	count: int,
	totalLength: number,
	intervalLength: number,
	passive: 'totalLength' | 'intervalLength' | 'count'}

function onChange(callback) {
	return function (proto, name) {
		return {
			set: function (x) {
				this['_' + name] = x
				callback(this)
			},
			get: function () { return this['_' +name] }}
	}
}
class Pattern extends Feature {
	fail = this
	static readonly SN: string = 'PTRN'
	name: string = 'pattern' + globalId++
	features: Feature[] = []
	direction: V3 = []
	directionFlipped: boolean = false
	count: int = 2
	totalLength: number = 40
	intervalLength: number = 20
	passive: 'totalLength' | 'intervalLength' | 'count' = 'totalLength'


	constructor() {
		super()
	}

	apply(publish, color, m4: M4 = M4.IDENTITY, genFeature: Feature = this) {
		const feature = this
		const sel = feature.direction.map(w => w.getOrThrow())
		const dir = MODES.SELECT_DIRECTION.magic(sel)[0].times(feature.directionFlipped ? -1 : 1)
		for (const featureRef of feature.features) {
			const patFeature = featureRef.getOrThrow()
			for (let i = 0; i < feature.count; i++) {
				const offset = dir.times((i + 1) * feature.intervalLength)
				const matrix = M4.translate( offset)
				patFeature.apply(publish, color, m4.times(matrix), genFeature)
			}
		}

	}

	dependentOnNames() {
		return []
	}

	addDimension() {
		this.dimensions.push({
			direction: V3.X,
			directionFlipped: false,

			count: 2,
			totalLength: 20,
			intervalLength: 10,
			passive: 'totalLength',
		})
	}
}


function modeGetName(mode) {
	return Object.getOwnPropertyNames(MODES).find(name => MODES[name] == mode)
}
function modeUpdateDisplay() {
	$('modeBox').set('html', modeStack.map(modeGetName).join('<br>'))
}
function modePush(mode: MODE, ...args: any[]) {
	assert(mode.init && mode.end && mode.mousemove && mode.mousedown)
	mode.before && mode.before()
	const feature = args[0]
	if (mode.modeEditsFeatureAndRequiresRollback && feature) {
		rebuildLimit = featureStack.indexOf(feature) + 1
		updateFeatureDisplay()
		rebuildModel()
	}
	modeStack.push(mode)
	mode.init.apply(mode, args)
	modeUpdateDisplay()
}
function modeEnd(mode: MODE) {
	if (modeStack.includes(mode)) {
		let popped
		do {
			popped = modeStack.pop()
			popped.end()
			if (popped.modeEditsFeatureAndRequiresRollback) {
				rebuildLimit = -1
				updateFeatureDisplay()
				rebuildModel()
			}
		} while (popped != mode)
		modeUpdateDisplay()
	}
}
function modePop() {
	if (1 == modeStack.length) return
	modeEnd(modeStack.last())
}
function modeGetCurrent() {
	return modeStack.last()
}

function tooltipShow(htmlContent: string, x: number, y: number) {
	const tipWrap = $('tip-wrap') as HTMLDivElement
	assert(tipWrap)
	const tip = $('tip') as HTMLDivElement
	tipWrap.setStyle('visibility', 'visible')
	tip.set('html', htmlContent)
	const size = tipWrap.getSize()
	tipWrap.setStyle('left', (x - 16 - 16) + 'px')
	if (window.innerHeight - y < size.y) {
		// place on top
		tipWrap.setStyle('top', (y - size.y - 16) + 'px')
		$('tip-arrow').setStyles({
			top: 'auto',
			bottom: '-16px',
			borderWidth: '16px 16px 0px'
		})
	} else {
		// place below
		$('tip-arrow').setStyles({
			top: '-16px',
			bottom: 'auto',
			borderWidth: '0px 16px 16px'
		})
		tipWrap.setStyle('top', (y + 16 + 8) + 'px')
	}
}
function tooltipHide() {
	$('tip-wrap').setStyle('visibility', 'hidden')
}


function serialize(v: {}) {
	function gatherList(v: {}) {
		if (Array.isArray(v)) {
			if (visited.has(v)) {
				if (!listMap.has(v)) {
					listMap.set(v, resultList.length)
					resultList.push(v)
				}
			} else {
				visited.add(v)
				for (let i = 0; i < v.length; i++) {
					gatherList(v[i])
				}
			}
		} else if (null !== v && 'object' == typeof v) {
			if (visited.has(v)) {
				if (!listMap.has(v)) {
					listMap.set(v, resultList.length)
					resultList.push(v)
				}
			} else {
				visited.add(v)
				const keys = Object.keys(v)
				for (let i = 0; i < keys.length; i++) {
					gatherList(v[keys[i]])
				}
			}
		}
	}

	function transform(v, first?) {
		if ('string' == typeof v || 'number' == typeof v || 'boolean' == typeof v || 'undefined' == typeof v || null == v) {
			return v
		} else if (v.constructor === Array) {
			let index
			if (true !== first && undefined !== (index = listMap.get(v))) {
				return {'#REF': index}
			} else {
				return v.map(transform)
			}
		} else if ('object' == typeof v) {
			let index
			if (true !== first && undefined !== (index = listMap.get(v))) {
				return {'#REF': index}
			} else {
				const result = Object.prototype == v.prototype ? {}
					: assert(v.constructor && v.constructor.name,
					() => (console.log(v), v.toSource() + v.constructor.name)) && {'#PROTO': v.constructor.name}
				const keys = Object.keys(v)
				for (let i = 0; i < keys.length; i++) {
					result[keys[i]] = transform(v[keys[i]])
				}
				return result
			}
		} else {
			throw new Error('?' + typeof v + v.toString())
		}
	}

	const visited = new Set()
	const listMap = new Map()
    let resultList: {}[] = []
	listMap.set(v, 0)
	resultList.push(v)
	gatherList(v)
	console.log(resultList)

	resultList = resultList.map(v => transform(v, true))
	console.log(JSON.stringify(resultList))
	return JSON.stringify(resultList)
}

function unserialize(string: string) {
	function fixObjects(v) {
		if (v && v.constructor === Array) {
			for (let i = 0; i < v.length; i++) {
				v[i] = fixObjects(v[i])
			}
			return v
		} else if ('object' == typeof v && null != v) {
			if ('#PROTO' in v) {
				const protoName = v['#PROTO'] as string
				const proto = (window[protoName] as string) || CLASSES[protoName] || /^[_$a-zA-Z\xA0-\uFFFF][_$a-zA-Z0-9\xA0-\uFFFF]*$/.test(protoName) && eval(protoName)
				assert(proto, protoName + ' Missing ' + window[protoName])
				const result = Object.create((proto).prototype)
				const keys = Object.keys(v)
				for (let i = 0; i < keys.length; i++) {
					//if ('name' == keys[i]) console.log(result)
					if ('#PROTO' != keys[i]) {
						Object.defineProperty(result, keys[i], {
							value: fixObjects(v[keys[i]]),
							enumerable: true,
							writable: true,
							configurable: true
						})
					}
				}
				Object.defineProperty(result, 'loadID', {value: globalId++, enumerable: false, writable: false})
				return result
			} else {
				return v
			}
		} else {
			return v
		}
	}

	function linkReferences(v) {
		if (v && v.constructor === Array) {
			for (let i = 0; i < v.length; i++) {
				v[i] = linkReferences(v[i])
			}
			return v
		} else if ('object' == typeof v && null != v) {
			if ('#REF' in v) {
				return tree[v['#REF']]
			} else {
				const keys = Object.keys(v)
				for (let i = 0; i < keys.length; i++) {
					v[keys[i]] = linkReferences(v[keys[i]])
				}
				return v
			}
		} else {
			return v
		}
	}

	const tree = JSON.parse(string)
	// console.log(tree)
	fixObjects(tree)
	// console.log(tree)
	linkReferences(tree)
	// console.log(tree)
	return tree[0]
}
