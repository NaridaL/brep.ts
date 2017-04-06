// window.onerror = function (errorMsg, url, lineNumber, column, errorObj) {
// 	// console.error(errorMsg, url, lineNumber, column, errorObj);
// 	console.error("%c"+errorMsg, 'color: black')
// 	if (errorObj) {
// 		console.log(formatStack(errorObj.stack.replace(/^(.*NLA\.assert.*\r?\n?)+/, '')))
// 		return true
// 	} else {
// 		return false
// 	}
// }

const MooEl: ElementConstructor = Element
function formatStack(stack) {
	let re = /\s*([^@]*)@(.+):(\d+):(\d+)\s*$/mg
	let match, matches = []

	let maxWidths = [0, 0, 0, 0]
	while (match = re.exec(stack)) {
		// console.log(match)
		matches.push(match)
		maxWidths = maxWidths.map((maxWidth, i) => max(maxWidth, match[i + 1].length))
	}
	let output = matches.map(match => {
		return [
			"\t",

			match[1],
			NLA.repeatString(maxWidths[0] + 2 - match[1].length, ' '),

			match[2].replace('file:///C:/Users/aval/Desktop/cs', ''),
			NLA.repeatString(maxWidths[1] + 2 - match[2].length, ' '),

			NLA.repeatString(maxWidths[2] + 2 - match[3].length, ' '),
			match[3],

			NLA.repeatString(maxWidths[2] + 2 - match[3].length, ' '),
			match[3]

		].join('')
	}).join('\n')
	return output
}


function Point(x, y) {
	this.x = x
	this.y = y
}


function makeCoincident(p1, p2, sketch) {
	if (p1.coincidence && p2.coincidence) {
		for (var i = 0; i < p2.coincidence.length; i++) {
			p1.coincidence.cs.push(p2.coincidence.cs[i])
			p2.coincidence.cs[i].coincidence = p1.coincidence
		}
		p1.coincidence.cs.push(p2)
		p2.coincidence = p1.coincidence
	} else if (p1.coincidence) {
		p1.coincidence.cs.push(p2)
		p2.coincidence = p1.coincidence
	} else if (p2.coincidence) {
		p2.coincidence.cs.push(p1)
		p1.coincidence = p2.coincidence
	} else {
		sketch.constraints.push(p1.coincidence = p2.coincidence = new Constraint("coincident", [p1, p2]))
	}
	p2.x = p1.x
	p2.y = p1.y
}
function removeFromCoincidence(p, sketch) {
	if (!p.coincidence) return
	p.coincidence.cs.remove(p)
	if (p.coincidence.cs.length == 1) {
		sketch.constraints.remove(p.coincidence)
		p.coincidence.cs[0].coincidence = null
	}
	p.coincidence = null
}
function deleteCoincidence(coincidence, sketch) {
	sketch.constraints.remove(coincidence)
	coincidence.cs.forEach(p => p.coincidence = null)
}


var savedTextures = new Map()
function getTextureForString(str) {
	var texture
	if (texture = savedTextures.get(str)) {
		return texture
	}

	var canvas = $("textTextureCanvas") as any as HTMLCanvasElement
	var ctx = canvas.getContext('2d')

	var font = TEXT_TEXTURE_HEIGHT + "px Anonymous Pro"

	ctx.font = font
	var textWidthPx = ctx.measureText(str).width
	var canvasWidth = 1 << ceil(log(textWidthPx) / log(2))
	canvas.width = canvasWidth
	canvas.height = TEXT_TEXTURE_HEIGHT


	ctx.font = font
	ctx.fillStyle = "#ffffff"
	ctx.textAlign = "left"	// This determines the alignment of text, e.g. left, center, right
	ctx.textBaseline = "top"	// This determines the baseline of the text, e.g. top, middle, bottom
	;
	(ctx as any).imageSmoothingEnabled = false


	ctx.fillText(str, 0, 0)


	texture = GL.Texture.fromImage(canvas, {minFilter: gl.LINEAR_MIPMAP_LINEAR, magFilter: gl.LINEAR})
	texture.textWidthPx = textWidthPx

	savedTextures.set(str, texture)

	return texture
}
var zoomFactor = 1
function paintSketch(sketch) {

	if (!sketch.plane || !sketch.sketchToWorldMatrix || !sketch.sketchToWorldMatrix.m) return
	if (-1 != rebuildLimit && rebuildLimit <= featureStack.indexOf(sketch)) return
	//console.log("painting segments", sketch.elements.length);
	/*ctx.clearRect (0, 0, ctx.canvas.width, ctx.canvas.height);
	 ctx.fillStyle="rgb(100, 100, 255)";
	 ctx.lineWidth=2;*/
	//console.log(sketch.sketchToWorldMatrix);
	gl.multMatrix(sketch.sketchToWorldMatrix)
	sketch.elements.forEach(function (seg) {
		function drawPoint(p) {
			paintArc(p.V3(), 1.5, 3, colorFor(highlighted.includes(p) || hoverHighlight == p, selected.includes(p)))
		}

		let color = colorFor(highlighted.includes(seg) || hoverHighlight == seg, selected.includes(seg))
		if (seg instanceof SketchLineSeg) {
			//console.log("seg", seg);
			//console.log("hoverHighlight.length", hoverHighlight.length);
			//ctx.beginPath();
			paintLineXY(seg.a.V3(), seg.b.V3(), color)


			//console.log(seg.a);
			drawPoint(seg.a)
			drawPoint(seg.b)
		} else if (seg instanceof SketchArc) {
			var radius = seg.radiusA()
			var angleA = seg.angleA(), angleB = seg.angleB()
			if (angleB <= angleA) { angleB += Math.PI * 2 }
			paintArc(seg.c.V3(), radius, 2, color, angleA, angleB)
			drawPoint(seg.c)
		} else if (seg instanceof SketchBezier) {
			// paintBezier(seg.points.map(p => p.V3()), 2, 0xdddddd, -2, 3)
			paintBezier(seg.points.map(p => p.V3()), 2, color, 0, 1)
			seg.points.forEach(drawPoint)
		} else {
			throw new Error("unknown sketch element" + seg)
		}
		drawPoint(seg.a)
		drawPoint(seg.b)

	})

	paintConstraints(sketch)
}
function colorFor(highlighted, selected) {
	return !selected
		? (!highlighted ? 0x33CCFF : 0x145266)
		: (!highlighted ? 0xFF3399 : 0x330A1E)
}
function getAllPoints(sketch) {
	return sketch.elements.map(segment => segment.points).concatenated()
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
	var crossCount = 2, parallelCrossCount = 1
	sketch.constraints.forEach(function (cst) {
		paintDefaultColor = hoverHighlight == cst ? (!cst.error ? 0x00ff00 : 0xffff00) : (!cst.error ? 0x000000 : 0xff0000)
		switch (cst.type) {
			case 'coincident':
				let point = cst.cs[0]
				paintArc(point.V3(), 4, 1, colorFor(false, false), 0, 2 * Math.PI)
				paintArc(point.V3(), 1.5, 3, colorFor(highlighted.includes(point) || hoverHighlight == point, selected.includes(point)), 0, 2 * Math.PI)
				break
			case 'parallel': {
				let dir1 = cst.segments[0].getVectorAB().unit()
				let dir90 = dir1.getPerpendicular().unit()
				let ARR_SPACING = 5
				for (let c = 0; c < cst.segments.length; c++) {
					let line = cst.segments[c]
					let ab = line.getVectorAB()
					let abLength = ab.length()
					let ab1 = ab.unit()
					for (let i = 0; i < parallelCrossCount; i++) {
						let s = abLength / 2 - ARR_SPACING * parallelCrossCount / 2 + i * ARR_SPACING - 10
						let arrPoint = line.a.V3().plus(ab1.times(s))
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
					let ab = cst.segment.getVectorAB()
					let cd = cst.other.getVectorAB()
					let intersection = cst.segment.intersection(cst.other)
					let abPos = cst.segment.pointT(intersection)
					let cdPos = cst.other.pointT(intersection)
					let abLine = ab.unit().times(0.5 < abPos ? -16 : 16)
					let cdLine = cd.unit().times(0.5 < cdPos ? -16 : 16)
					paintLineXY(intersection.plus(abLine).plus(cdLine), intersection.plus(abLine))
					paintLineXY(intersection.plus(abLine).plus(cdLine), intersection.plus(cdLine))
				} else {
					let ab = cst.segment.getVectorAB()
					let sketchLineWC = sketch.plane.intersectionWithPlane(cst.other.plane)
					let dir = sketch.worldToSketchMatrix.transformVector(sketchLineWC.dir1)
					let p = sketch.worldToSketchMatrix.transformPoint(sketchLineWC.anchor)
					let abXcd = ab.cross(dir)
					let div = abXcd.squared()
					let anchorDiff = p.minus(cst.segment.a.V3())
					let abPos = anchorDiff.cross(dir).dot(abXcd) / div
					let linePos = anchorDiff.cross(ab).dot(abXcd) / div
					let intersection = p.plus(dir.times(linePos))
					let abLine = ab.unit().times(0.5 < abPos ? -16 : 16)
					let cdLine = dir.times(16)
					paintLineXY(intersection.plus(abLine).plus(cdLine), intersection.plus(abLine))
					paintLineXY(intersection.plus(abLine).plus(cdLine), intersection.plus(cdLine))
				}
				break
			case 'colinear':
				// assume segments dont overlap
				let segments = cst.segments
				let ab = segments[0].getVectorAB().unit()
				let coord = ab.x > ab.y ? 'x' : 'y'
				if (ab[coord] < 0) {
					ab = ab.times(-1)
				}
				let scale = ab[coord]
				let offsetPoint = segments[0].a.V3().minus(ab.times(segments[0].a.V3()[coord] / scale))
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
				let a = cst.cs[0].V3(), b = cst.cs[1].V3()
				let ab = b.minus(a)
				let ab1 = ab.unit()
				let abLength = ab.length()
				let ab90 = ab1.getPerpendicular()
				let texture = getTextureForString(cst.distance)
				let textLength = texture.textWidthPx / TEXT_TEXTURE_HEIGHT * 20
				paintLineXY(a.plus(ab90.times(6)), a.plus(ab90.times(22)))
				paintLineXY(b.plus(ab90.times(6)), b.plus(ab90.times(22)))

				paintLineXY(a.plus(ab90.times(14)), a.plus(ab90.times(14)).plus(ab1.times(abLength / 2 - textLength / 2 - 5)))
				paintLineXY(a.plus(ab90.times(14)).plus(ab1.times(abLength / 2 + textLength / 2 + 5)), a.plus(ab90.times(14)).plus(ab1.times(abLength)))
				let textCenter = a.plus(ab90.times(14)).plus(ab1.times(abLength / 2))
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
				let texture = getTextureForString(cst.distance)
				let p = cst.point.V3(), ab = cst.other.getVectorAB(), a = cst.other.a.V3()
				let ab1 = ab.unit()
				let abLength = ab.length()
				let ab90 = ab1.getPerpendicular()
				let ap = p.minus(a)
				let apProj = ab90.times(-ap.dot(ab90))
				paintLineXY(a, a.plus(apProj))
				paintLineXY(a.plus(apProj), p)
				let textLength = texture.textWidthPx / TEXT_TEXTURE_HEIGHT * 20
				let textCenter = a.plus(apProj.times(1 / 2))
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
				let texture = getTextureForString(cst.distance)
				let p = cst.point.V3()
				let sketchLineWC = sketch.plane.intersectionWithPlane(cst.other.plane)
				let sketchLineSC = sketchLineWC.transform(sketch.worldToSketchMatrix)
				let a = sketchLineSC.closestPointToPoint(p)
				let ap = p.minus(a)
				paintLineXY(a, p)
				let textLength = texture.textWidthPx / TEXT_TEXTURE_HEIGHT * 20
				let textCenter = a.plus(ap.times(1 / 2))
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
				let curve = getCachedCurve(cst.other, sketch.plane)
				let t = curve.closestTToPoint(cst.point.V3())
				let tangentLength = curve.tangentAt(t).length()
				CURVE_PAINTERS[curve.constructor.name](curve, 0x000000, t - 12 / tangentLength, t - 4 / tangentLength, 3)
				CURVE_PAINTERS[curve.constructor.name](curve, 0x000000, t + 4 / tangentLength, t + 12 / tangentLength, 3)
				break
			}
			case 'angle':
				let first = cst.cs[0], second = cst.cs[1]
				let intersection = first.intersection(second)
				let startAngle = (first.angleAB() + (cst.f[0] == -1 ? PI : 0)) % (2 * PI), endAngle = startAngle + cst.value
				paintArc(intersection, 20, 1, 0x000000, startAngle, endAngle)
				break
			case 'equalLength': {
				for (let c = 0; c < cst.segments.length; c++) {
					let line = cst.segments[c]
					let ab = line.getVectorAB()
					let abLength = ab.length()
					let ab1 = ab.unit()
					let ab90 = ab.getPerpendicular().unit()
					let crossLength = 10
					let crossSpacing = 3
					for (let i = 0; i < crossCount; i++) {
						let s = abLength / 2 - crossSpacing * crossCount / 2 + i * crossSpacing + 10
						let crossStart = line.a.V3().plus(ab1.times(s)).plus(ab90.times(crossLength / 2))
						let crossEnd = line.a.V3().plus(ab1.times(s)).plus(ab90.times(-crossLength / 2))
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

class PlaneDefinition {
	/*
	 if ("face" == feature.planeType && feature.faceName) {
	 var face = modelBREP.faces.find(face => face.name == feature.faceName)
	 var plane = face.surface.plane, right = plane.normal.getPerpendicular().unit(), up = plane.normal.cross(right)

	 var cp
	 planes.push(cp = new CustomPlane(plane.anchor.plus(plane.normal.times(feature.offset)),
	 right, up, -500, 500, -500, 500, 0xFF4500, feature.planeId))
	 }
	 if ("immediate" == feature.planeType) {
	 var plane = eval(feature.source), right = plane.normal.getPerpendicular().unit(), up = plane.normal.cross(right)

	 planes.push(new CustomPlane(plane.anchor,
	 right, up, -500, 500, -500, 500, 0xFF4500, feature.planeId))
	 }
	 */
	type = "planeDefinition"
	planeId
	planeType = 'dynamic'
	source = ''
	offset = 0
	angle = 0
	whats = []
	flipped = false

	constructor(planeId) {
		assert(planeId)
		this.planeId = planeId
	}

	toSource(line) {
		return `new PlaneDefinition()`
	}

	dependentOnNames() {
		return this.whats
	}
}
NLA.registerClass(PlaneDefinition)

function serializeFeatures(features) {
	return serialize(features)
	//return '[' + featureStack.map(f => f.toSource()).join(',') + ']'
}

function load(key) {
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
	let loadSelect = $('loadSelect') as HTMLSelectElement
	let keys = NLA.arrayFromFunction(localStorage.length, i => localStorage.key(i))
	keys.sort()
	keys.forEach(key => loadSelect.adopt(new MooEl('option', {html: key})))
	loadSelect.onchange = function () {
		let key = loadSelect.value
		if (key) {
			load(key)
		}
	}


	let saveButton = $<HTMLButtonElement>('saveButton')
	saveButton.onclick = function () {
		let key = $<HTMLInputElement>('saveNameInput').value
		if (null == key) {
			loadSelect.adopt(new MooEl('option', {html: key}))
		}
		localStorage.setItem(key, serializeFeatures(featureStack))
		console.log("saved " + key)
	}
}

function isSketchEl(el) {
	return el instanceof SketchBezier || el instanceof SegmentEndPoint
		|| el instanceof SketchLineSeg || el instanceof SketchArc
}

//var sketchPlane = new CustomPlane(V3(0, 0,1), V3.X, V3.Y, -500, 500, -500, 500, 0xff00ff);
var editingSketch: Sketch, featureStack = []
function initModel() {

}
function directionObjectToVector(dirObj) {
	if (dirObj instanceof CustomPlane) {
		return dirObj.normal
	} else {
		console.log(dirObj)
		throw new Error("uuuh" + dirObj)
	}
}
function rebuildModel() {
	let DISABLE_CONSOLE = false
	DISABLE_CONSOLE && disableConsole()
	console.log("rebuilding model")
	//NLA.DEBUG = false


	namesPublishedBy = new Map()
	publishedObjects = new Map()
	function publish(name, object) {
		if (!namesPublishedBy.has(name)) {
			namesPublishedBy.set(name, feature)
			publishedObjects.set(name, object)
		}

	}


	modelBREP = undefined
	brepMesh = undefined
	brepPoints = []
	brepEdges = []

	planes = [
		new CustomPlane(V3.O, V3.Y, V3.Z, -500, 500, -500, 500, 0xffaaaa, "planeYZ"),
		new CustomPlane(V3.O, V3.Z, V3.X, -500, 500, -500, 500, 0xaaffaa, "planeZX"),
		new CustomPlane(V3.O, V3.X, V3.Y, -500, 500, -500, 500, 0xaaaaff, "planeXY"),
	]
	planes.forEach(customPlane => publish(customPlane.name, customPlane))

	let loopEnd = -1 == rebuildLimit ? featureStack.length : min(rebuildLimit, featureStack.length)
	for (var featureIndex = 0; featureIndex < loopEnd; featureIndex++) {
		var feature = featureStack[featureIndex]
		//try {
		if (feature instanceof Sketch) {
			feature.plane = feature.planeRef.getOrThrow()
			//console.log("LENGTHS", feature.plane.right.length(), feature.plane.up.length(),
			// feature.plane.normal.length())
			feature.sketchToWorldMatrix =
				M4.forSys(feature.plane.right, feature.plane.up, feature.plane.normal, feature.plane.anchor)
			feature.worldToSketchMatrix = feature.sketchToWorldMatrix.inversed()
			feature.recalculate()
			feature.elements.forEach(el => publish(el.name, el))
		} else if (feature.type && feature.type == "extrude") {

			let loopSegment = feature.segmentName.getOrThrow()
			let loopSketch = loopSegment.sketch
			// opposite dir to plane normal:
			let startOffset = loopSketch.plane.normal.times(-min(feature.start, feature.end))
			let lengthOffset = loopSketch.plane.normal.times(-Math.abs(feature.start - feature.end))
			let startMatrix = M4.translation(startOffset).times(loopSketch.sketchToWorldMatrix)
			let edgeLoop = loopSketch.getLoopForSegment(loopSegment)
			edgeLoop = edgeLoop.map(edge => edge.transform(startMatrix, ''))
			// TODO: test for self-intersection of edgeloop
			if (!new PlaneSurface(loopSketch.plane, loopSketch.right, loopSketch.up).edgeLoopCCW(edgeLoop)) {
				edgeLoop = edgeLoop.map(edge => edge.flipped()).reverse()
			}
			//console.log(polygonPoints.map(v =>v.$))
			let length = feature.end - feature.start
			//console.log("polypoints", polygonPoints, polygonPoints.toSource(),
			// loopSketch.plane.translated().toSource(), offsetDir.times(feature.start))
			let brep = B2T.extrudeEdges(edgeLoop, loopSketch.plane.translated(startOffset), lengthOffset, feature.name)

			if (modelBREP) {
				// isEdges = modelBREP.getIntersectionEdges(brep)
				// drVs = isEdges.map(e => ({anchor: e.a, dir: e.curve.tangentAt(e.aT).unit()}))
				modelBREP = modelBREP[feature.operation](brep)
			} else {
				modelBREP = brep
			}

			modelBREP.faces.forEach(face => publish(face.name, face))
			modelBREP.faces.map(face => face.getAllEdges().filter(edge => !edge.flippedOf || edge.id < edge.flippedOf.id).forEach(edge => publish(edge.name, edge)))
			console.log('modelBREP.vertexNames', modelBREP.vertexNames)
			modelBREP.vertexNames && modelBREP.vertexNames.forEach((name, p) => publish(name, p))
		} else if (feature.type && "planeDefinition" == feature.type) {
			let sel = feature.whats.map(w => w.getOrThrow())
			let plane = MODES.PLANE_DEFINITION.magic(sel, feature.angle * DEG)
			if (plane) {
				if (feature.flipped) plane = plane.flipped()
				let right = plane.normal.getPerpendicular().unit(), up = plane.normal.cross(right)

				var cp
				planes.push(cp = new CustomPlane(plane.anchor.plus(plane.normal.times(feature.offset)),
					right, up, -500, 500, -500, 500, 0xFF4500, feature.planeId))
			}
			if ("immediate" == feature.planeType) {
				let plane = eval(feature.source), right = plane.normal.getPerpendicular().unit(), up = plane.normal.cross(right)

				planes.push(new CustomPlane(plane.anchor,
					right, up, -500, 500, -500, 500, 0xFF4500, feature.planeId))
			}
			publish(feature.planeId, cp)
		} else {
			//noinspection ExceptionCaughtLocallyJS
			throw new Error("unknown feature")
		}
		// brepMesh.computeWireframeFromFlatTriangles()
		// brepMesh.compile()

		//} catch (error) {
		//	let featureDiv = $('featureDisplay').getChildren().filter(child => child.featureLink == feature)[0]
		//
		//	if (featureDiv) {
		//		let ediv = featureDiv.getElement('[name=error]')
		//		ediv.setStyle('display', 'inline')
		//		ediv.title= error.toString() + '\n' + error.stack
		//	}
		//	console.error(error)
		//	// throw error
		//	break
		//}
	}
	if (modelBREP) {
		brepMesh = modelBREP.toMesh()
		modelBREP.faces.forEach(face =>
			brepEdges.pushAll(face.getAllEdges().filter(edge => !edge.flippedOf || edge.id < edge.flippedOf.id)))
		// brepEdges.pushAll(isEdges)
		modelBREP.vertexNames && (brepPoints = Array.from(modelBREP.vertexNames.keys()))
		console.log(modelBREP.vertexNames)
		brepMesh.computeNormalLines(10)
		//brepMesh.computeWireframeFromFlatTriangles()
		brepMesh.compile()
	}
	updateSelected()
	paintScreen()
	DISABLE_CONSOLE && enableConsole()
}
// DECLARATIONS
// GLOBALS
var missingEls = [], namesPublishedBy, publishedObjects
var faces = []
var highlighted = [], selected = [], paintDefaultColor = 0x000000
var gl: GL.LightGLContext
var modelBREP, brepMesh, brepPoints, planes, brepEdges, isEdges = []
var rebuildLimit = -1, rollBackIndex = -1
var drPs = [], drVs = []

var eyePos = V(1000, 1000, 1000), eyeFocus = V3.O, eyeUp = V3.Z
var hoverHighlight = undefined
// console.log = oldConsole;
var modeStack = []
var shaders: any = {}

function rgbToVec4(color) {
	return [(color >> 16) / 255.0, ((color >> 8) & 0xff) / 255.0, (color & 0xff) / 255.0, 1.0]
}
function renderColor(mesh, color, mode?) {
	shaders.singleColor.uniforms({
		color: rgbToVec4(color)
	}).draw(mesh, mode)
}
function renderColorLines(mesh, color) {
	shaders.singleColor.uniforms({
		color: rgbToVec4(color)
	}).draw(mesh, 'LINES')
}
var TEXT_TEXTURE_HEIGHT = 128
function renderText(string, color) {
	var texture = getTextureForString(string)
	texture.bind(0)
	gl.pushMatrix()
	gl.scale(texture.width / TEXT_TEXTURE_HEIGHT, 1, 1)
	shaders.textureColor.uniforms({
		texture: 0,
		color: rgbToVec4(color)
	}).draw(meshes.text)
	gl.popMatrix()
}
function drawVector(vector, anchor, color = 0x0000ff, size = 50) {
	gl.pushMatrix()
	let v2 = vector.getPerpendicular()
	gl.multMatrix(M4.forSys(vector, v2, vector.cross(v2), anchor))
	gl.scale(size, size, size)
	shaders.singleColor.uniforms({
		color: rgbToVec4(color)
	}).draw(meshes.vector)
	gl.popMatrix()
}
function drawVectors() {
	drawVector(V3.X, V3.O, 0xff0000)
	drawVector(V3.Y, V3.O, 0x00ff00)
	drawVector(V3.Z, V3.O, 0x0000ff)

	drVs.forEach(vi => drawVector(vi.dir1, vi.anchor, vi.color))
}

function drawPoint(p, color = 0x000000, size = 5) {
	gl.pushMatrix()
	gl.translate(p)
	gl.scale(size, size, size)
	shaders.singleColor.uniforms({color: rgbToVec4(color)}).draw(meshes.sphere1)
	gl.popMatrix()
}
function drawPoints(size: number) {
	drPs.forEach(info => drawPoint(info.p || info, 0xcc0000, size))
	brepPoints && brepPoints.forEach(p =>
		drawPoint(p, hoverHighlight == p ? 0x0adfdf : (selected.includes(p) ? 0xff0000 : 0xcccc00), 2))
}


const CURVE_PAINTERS: {[curveConstructorName: string]: (curve: Curve, color: int, startT: number, endT: number, width: number) => void} = {}
CURVE_PAINTERS[SemiEllipseCurve.name] = function paintEllipseCurve(ellipse: SemiEllipseCurve, color, startT, endT, width = 2) {
	shaders.ellipse3d.uniforms({
		f1: ellipse.f1,
		f2: ellipse.f2,
		center: ellipse.center,
		color: rgbToVec4(color),
		startT: startT,
		endT: endT,
		scale: width,
		mode: 0 // ellipse
	}).draw(meshes.pipe)
}
CURVE_PAINTERS[BezierCurve.name] = function paintBezierCurve(curve: BezierCurve, color, startT, endT, width = 2, normal = V3.Z) {
	shaders.bezier3d.uniforms({
		p0: curve.p0,
		p1: curve.p1,
		p2: curve.p2,
		p3: curve.p3,
		color: rgbToVec4(color),
		startT: startT,
		endT: endT,
		scale: width,
		normal: normal
	}).draw(meshes.pipe)
}
CURVE_PAINTERS[L3.name] = function (curve: L3, color, startT, endT, width = 2, normal = V3.Z) {
	gl.pushMatrix()
	let a = curve.at(startT), b = curve.at(endT)
	let ab = b.minus(a), abT = ab.getPerpendicular().unit()
	let m = M4.forSys(ab, abT, ab.cross(abT).unit(), a)
	gl.multMatrix(m)
	gl.scale(1, width, width)
	shaders.singleColor.uniforms({
		color: rgbToVec4(color), // TODO: error checking
	}).draw(meshes.pipe)

	gl.popMatrix()
}


function paintLineXY(a, b, color?, width?) {
	color = color || paintDefaultColor
	width = width || 2
	var ab = b.minus(a)
	if (ab.isZero()) { return }
	var abT = ab.getPerpendicular().toLength(width)
	//console.log(ab)
	gl.pushMatrix()
	gl.multMatrix(M4.forSys(ab, abT, V3.Z, a))
	renderColor(meshes.segment, color)
	gl.popMatrix()
}
function paintArc(center, radius, width, color, startAngle = 0, endAngle = 2 * Math.PI) {
	// startAngle = isFinite(startAngle) ? startAngle : 0
	// endAngle = isFinite(endAngle) ? endAngle : 2 * Math.PI
	gl.pushMatrix()
	gl.translate(center)
	shaders.arc2.uniforms({
		color: rgbToVec4(color),
		offset: startAngle,
		step: endAngle - startAngle,
		radius: radius,
		width: width
	}).draw(meshes.segment)
	gl.popMatrix()
}

function paintBezier(ps, width, color, startT, endT) {
	// TODO PS AREN'T IN THE ORDER YOU EXPECT THEM TO BE!!!
	assertVectors.apply(undefined, ps)
	shaders.bezier.uniforms({
		color: rgbToVec4(color),
		width: width,
		p0: ps[0],
		p1: ps[2],
		p2: ps[3],
		p3: ps[1],
		startT: startT || 0,
		endT: endT || 1
	}).draw(meshes.segment)
}

function randomVec4Color(opacity) {
	opacity = opacity || 1.0
	return [Math.random(), Math.random(), Math.random(), opacity]
}
function drawEdge(edge, color = 0x000000, width = 2) {
	const startT = min(edge.aT, edge.bT)
	const endT = max(edge.aT, edge.bT)
	CURVE_PAINTERS[edge.curve.constructor.name](edge.curve, color, startT, endT, width)
}
var cachedMeshes = new Map()
function drawFace(face, color) {
	let mesh = cachedMeshes.get(face)
	if (!mesh) {
		mesh = new GL.Mesh({triangles: true, lines: true, normals: true})
		face.addToMesh(mesh)
		mesh.compile()
		cachedMeshes.set(face, mesh)
	}
	shaders.singleColor.uniforms({color: rgbToVec4(color)}).draw(mesh)
}
var /**@type GL.Mesh */ mesh1, mesh2
var meshes: any = {}
var ZERO_EL = {}
function drawEl(el, color) {
	if (el instanceof V3) {
		drawPoint(el, color)
	} else if (el instanceof Edge) {
		drawEdge(el, color)
	} else if (el instanceof Face) {
		drawFace(el, color)
	} else if (el instanceof CustomPlane) {
		drawPlane(el, color)
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
	// paintEllipseCurve(EllipseCurve.forAB(200, 300), COLORS.RD_FILL, 0, 2 * Math.PI)
	// paintBezierCurve(BezierCurve.EX3D.scale(100, 100, 100), COLORS.PP_STROKE, -2, 3)


// 	if (!mesh1) {
// //		let pF = modelBREP.faces.filter(face => face.constructor == RotationFace)[1]
// 		let pf = new RotationFace(
// 			new ProjectedCurveSurface(new BezierCurve(V3(142.87578921496748, -191.46078243076332, 0),
// V(161.78547089700214, -252.13248349581008, 0), V(284.63214994898954, -163.59789158697575, 0), V(372.40411211189405,
// -210.3992206435476, 0)), V(0, 0, 1), 0, 1), [ new PCurveEdge(new BezierCurve(V3(142.87578921496748,
// -191.46078243076332, 0), V(161.78547089700214, -252.13248349581008, 0), V(284.63214994898954, -163.59789158697575,
// 0), V(372.40411211189405, -210.3992206435476, 0)), V(372.40411211189405, -210.3992206435476, 0),
// V(142.87578921496748, -191.46078243076332, 0), 1, 0, new PCurveEdge(new BezierCurve(V3(142.87578921496748,
// -191.46078243076332, 0), V(161.78547089700214, -252.13248349581008, 0), V(284.63214994898954, -163.59789158697575,
// 0), V(372.40411211189405, -210.3992206435476, 0)), V(142.87578921496748, -191.46078243076332, 0),
// V3(372.40411211189405, -210.3992206435476, 0), 0, 1, null, V3(142.87578921496748, -191.46078243076332, 0),
// V3(372.40411211189405, -210.3992206435476, 0)), V3(-372.40411211189405, 210.3992206435476, 0),
// V3(-142.87578921496748, 191.46078243076332, 0)), StraightEdge.throughPoints(V3(142.87578921496748, -191.46078243076332, 0), V3(142.87578921496748, -191.46078243076332, -100)), new PCurveEdge(new BezierCurve(V3(142.87578921496748, -191.46078243076332, -100), V3(161.78547089700214, -252.13248349581008, -100), V3(284.63214994898954, -163.59789158697575, -100), V3(372.40411211189405, -210.3992206435476, -100)), V3(142.87578921496748, -191.46078243076332, -100), V3(372.40411211189405, -210.3992206435476, -100), 0, 1, null, V3(142.87578921496748, -191.46078243076332, 0), V3(372.40411211189405, -210.3992206435476, 0)), StraightEdge.throughPoints(V3(372.40411211189405, -210.3992206435476, -100), V3(372.40411211189405, -210.3992206435476, 0))], []) sphere = new GL.Mesh({lines: true, normals: true}) sphere.faceIndexes = new Map() pf.addToMesh(sphere) sphere.computeWireframeFromFlatTriangles() console.log("sphere.lines", sphere.lines, sphere.vertices, sphere.triangles) sphere.compile() }


	gl.pushMatrix()
	!mesh1.hasBeenCompiled && mesh1.compile()
	mesh1 && (mesh1.normals ? shaders.lighting : shaders.singleColor).uniforms({
		color: rgbToVec4(0xff0000)
	}).draw(mesh1)
	mesh1 && mesh1.curve1 && shaders.singleColor.uniforms({
		color: rgbToVec4(0xff0000)
	}).drawBuffers({LGL_Vertex: mesh1.vertexBuffers.curve1}, undefined, gl.LINES)
	mesh1 && mesh1.lines && shaders.singleColor.uniforms({
		color: rgbToVec4(0x0000ff)
	}).draw(mesh1, 'LINES')

	mesh2 && (mesh2.normals ? shaders.lighting : shaders.singleColor).uniforms({
		color: rgbToVec4(0x00ff00)
	}).draw(mesh2)
	mesh2 && mesh2.lines && shaders.singleColor.uniforms({
		color: rgbToVec4(0x0ff000)
	}).draw(mesh2, 'LINES')
	gl.popMatrix()


	gl.pushMatrix()
	if (brepMesh) {
		// paint faces
		let faceIndex = modelBREP.faces.length
		while (faceIndex--) {

			let face = modelBREP.faces[faceIndex]
			let faceTriangleIndexes = brepMesh.faceIndexes.get(face)
			shaders.lighting.uniforms({
				color: rgbToVec4(hoverHighlight == face ? 0xff00ff : (selected.includes(face) ? 0x00ff45 : COLORS.RD_FILL))
			}).draw(brepMesh, 'TRIANGLES', faceTriangleIndexes.start, faceTriangleIndexes.count)
			/*
			 shaders.singleColor.uniforms({
			 color: rgbToVec4(0x0000ff)
			 }).draw(brepMesh, 'LINES');
			 */
		}

		// paint edges
		isEdges.forEach(edge => drawEdge(edge, hoverHighlight == edge ? 0xff00ff : (selected.includes(edge) ? 0x00ff45 : 0x000000), 2))
		brepEdges.forEach(edge => drawEdge(edge, hoverHighlight == edge ? 0xff00ff : (selected.includes(edge) ? 0x00ff45 : COLORS.RD_STROKE), 2))

		gl.projectionMatrix.m[11] -= 1 / (1 << 22) // prevent Z-fighting
		shaders.singleColor.uniforms({
			color: rgbToVec4(COLORS.PP_STROKE)
		}).draw(brepMesh, 'LINES')
		gl.projectionMatrix.m[11] += 1 / (1 << 22) // prevent Z-fighting
	}
	gl.popMatrix()

	const DZ = 0.1
	gl.projectionMatrix.m[11] -= DZ // prevent Z-fighting
	featureStack.filter(f => f instanceof Sketch && !f.hide).forEach(s => {
		gl.pushMatrix()
		paintSketch(s)
		gl.popMatrix()
	})
	gl.projectionMatrix.m[11] += DZ
}

function setupCamera() {
	//console.log("eyePos", eyePos.$, "eyeFocus", eyeFocus.$, "eyeUp", eyeUp.$)
	gl.matrixMode(gl.PROJECTION)
	gl.loadIdentity()
	//gl.perspective(70, gl.canvas.width / gl.canvas.height, 0.1, 1000);
	var lr = gl.canvas.width / 2 / zoomFactor
	var bt = gl.canvas.height / 2 / zoomFactor
	gl.ortho(-lr, lr, -bt, bt, -1e5, 1e5)
	gl.lookAt(eyePos, eyeFocus, eyeUp)
	gl.matrixMode(gl.MODELVIEW)
}
function drawPlane(customPlane, color) {
	gl.pushMatrix()
	gl.multMatrix(M4.forSys(customPlane.right, customPlane.up, customPlane.normal))
	gl.translate(customPlane.rightStart, customPlane.upStart, customPlane.w)
	gl.scale(customPlane.rightEnd - customPlane.rightStart, customPlane.upEnd - customPlane.upStart, 1)

	let shader = hoverHighlight == customPlane ? shaders.singleColorHighlight : shaders.singleColor

	shader.uniforms({color: rgbToVec4(color)}).draw(meshes.xyLinePlane, 'LINES')

	gl.popMatrix()
}

function segmentIntersectsRay(a, b, p, dir) {
	var ab = b.minus(a)
	var abXdir = ab.cross(dir)
	var div = abXdir.squared()
	var anchorDiff = p.minus(a)
	var s = anchorDiff.cross(dir).dot(abXdir) / div
	var t = anchorDiff.cross(ab).dot(abXdir) / div
	//console.log(segmentIntersectsRay, a, b, "ab", ab, "p", p, "dir", dir, s > 0 && t / div >= 0 && t / div <= 1, "s",
	// s, "t", t, "div", div)
	return t > 0 && s >= 0 && s <= 1
}
function template(templateName, map): HTMLElement {
	var html = $<HTMLScriptElement>(templateName).text
	for (var key in map) {
		html = html.replace(new RegExp('\\$' + key, 'g'), map[key])
	}
	return $(new MooEl('div', {html: html.trim()}).firstChild)

}
function featureRollBack(feature: any, featureIndex: number) {
	rollBackIndex = featureIndex
}
function updateFeatureDisplay() {
	var div = $('featureDisplay')
	div.erase('text')
	featureStack.forEach(function (feature, featureIndex) {
		if (rebuildLimit == featureIndex) {
			div.adopt(new MooEl('div', {text: 'REBUILD LIMIT', class: 'rebuildLimit'}))
			// div.adopt(<div class='rebuildLimit'>REBUILD LIMIT</div>)
		}
		var newChild
		if (feature.type == "extrude") {
			newChild = template('templateFeatureExtrude', {what: "EXTRUDE", name: feature.name})
			newChild.getElement('[name=edit]').onclick = function () {
				modePush(MODES.EXTRUDE, feature)
			}
		} else if (feature.type == "planeDefinition") {
			newChild = template('templateFeatureExtrude', {what: "PLANE", name: feature.planeId})
			newChild.getElement('[name=edit]').onclick = function () {
				modePush(MODES.PLANE_DEFINITION, feature)
			}
		} else if (feature instanceof Sketch) {
			newChild = template('templateFeatureExtrude', {what: "SKETCH", name: feature.name})
			newChild.getElement('[name=edit]').onclick = function () {
				modePush(MODES.SKETCH, feature)
			}
		} else {
			console.log(feature)
			throw new Error("Unknown feature" + feature.toSource() + feature.constructor.name)
		}
		newChild.inject(div)
		newChild.getElement('[name=delete]').onclick = function () {
			featureDelete(feature)
		}
		newChild.getElement('[name=rollBack]').onclick = function () {
			featureRollBack(feature, featureIndex)
		}
		newChild.featureLink = feature
		newChild.onmouseover = function (e) {
			const dependencies = featureDependencies(feature)
			const dependents = featureDependents(feature)
			div.getChildren().filter((subDiv: any) => dependencies.includes(subDiv.featureLink)).addClass('isDependedOn')
			div.getChildren().filter((subDiv: any) => dependents.includes(subDiv.featureLink)).addClass('hasDependents')
		}
		newChild.onmouseout = function (e) {
			div.getChildren().removeClass('isDependedOn').removeClass('hasDependents')
		}
		feature.hide && newChild.getElement('[name=toggleHide]').addClass('hidden')
		newChild.getElement('[name=toggleHide]').onclick = function () {
			feature.hide = !feature.hide
			this.toggleClass('hidden', feature.hide)
			paintScreen()
		}
	})
}
function updateSelected() {
	var div = $('selectedElements')
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
		var newChild = template("template", {what: sel.constructor.name, name: name})
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
		newChild.getElement(".remove").onclick = function (e) {
			editingSketch.removeElement(target)
			updateSelected()
			paintScreen()
		}
		newChild.getElement(".unsel").onclick = function (e) {
			selected.remove(sel)
			updateSelected()
			paintScreen()
		}
		sel.toBrepEdge && newChild.grab(new MooEl('span', {
			text: sel.toBrepEdge().curve.toSource(x => NLA.round10(x, -3)),
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
			var newChild
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
				var input = newChild.getElement('.distanceInput')
				newChild.getElement('.distanceInput').value = NLA.round10(rad2deg(cst.value), -5)
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
				newChild = template("templateConstraint", {name: cst.type})

				cst.cs.forEach(function (el) {
					var subChild = template("templateConstraintSub", {what: el.constructor.name, name: el.name})
					subChild.inject(newChild)
					subChild.getElement(".removeFromConstraint").onclick = function (e) {
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
			newChild.getElement(".remove").onclick = function (e) {
				deleteConstraint(editingSketch, cst)
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

	ref: string
	lastHit: any

	constructor(name, lastHit) {
		let x = NameRef.pool.get(name)
		if (x) return x
		this.ref = name
		this.lastHit = lastHit
		NameRef.pool.set(name, this)
	}

	get() {
		let hit = publishedObjects.get(this.ref)
		// let hit = planes.find(plane => plane.name == this.ref)
		// 	|| modelBREP && modelBREP.faces.find(face => face.name == this.ref)
		// 	|| modelBREP && modelBREP.faces.firstMatch(face => face.getAllEdges().find(edge => edge.name ==
		// this.ref)) || modelBREP && modelBREP.vertexNames && mapReverse(modelBREP.vertexNames, this.ref)
		if (hit) this.lastHit = hit
		return hit
	}

	getOrThrow() {
		let hit = this.get()
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

	get what() {
		return this.get() instanceof Face ? "face" : "plane"
	}

	get name() {
		let hit = this.get()
		return hit && hit.name
	}

	static forObject(o) {
		if (o instanceof V3) {
			return new NameRef(modelBREP.vertexNames.get(o), o)
		} else {
			assert(o.name)
			return new NameRef(o.name, o)
		}
	}

	static UNASSIGNED = new NameRef('UNASSIGNED', ZERO_EL)
}
NameRef.UNASSIGNED.get = function () { return undefined }
NameRef.UNASSIGNED.getOrThrow = function () { throw new Error('this nameref has never been assigned a value') }
NLA.registerClass(NameRef)


function mapReverse(map, value) {
	let it = map.keys(), key
	while (key = it.next().value) {
		if (map.get(key) == value) {
			return key
		}
	}
}
function initMeshes() {
	meshes.sphere1 = GL.Mesh.sphere(2)
	meshes.segment = GL.Mesh.plane({startY: -0.5, height: 1, detailX: 128})
	meshes.text = GL.Mesh.plane()
	meshes.vector = GL.Mesh.rotation([V3.O, V(0, 0.05, 0), V(0.8, 0.05), V(0.8, 0.1), V(1, 0)], L3.X, Math.PI * 2, 16, true)
	meshes.pipe = GL.Mesh.rotation(NLA.arrayFromFunction(128, i => new V3(i / 127, -0.5, 0)), L3.X, Math.PI * 2, 8, true)
	meshes.xyLinePlane = GL.Mesh.plane()
}
function initShaders() {
	shaders.singleColor = new GL.Shader(vertexShaderBasic, fragmentShaderColor)
	shaders.singleColorHighlight = new GL.Shader(vertexShaderBasic, fragmentShaderColorHighlight)
	shaders.textureColor = new GL.Shader(vertexShaderTexture, fragmentShaderTextureColor)
	shaders.arc = new GL.Shader(vertexShaderRing, fragmentShaderColor)
	shaders.arc2 = new GL.Shader(vertexShaderArc, fragmentShaderColor)
	shaders.ellipse3d = new GL.Shader(vertexShaderConic3d, fragmentShaderColor)
	shaders.bezier3d = new GL.Shader(vertexShaderBezier3d, fragmentShaderColor)
	shaders.bezier = new GL.Shader(vertexShaderBezier, fragmentShaderColor)
	shaders.lighting = new GL.Shader(vertexShaderLighting, fragmentShaderLighting)
	shaders.waves = new GL.Shader(vertexShaderWaves, fragmentShaderLighting)
}
function initOtherEvents() {
	window.onkeypress = function (e) {
		if ("Delete" == e.key) {
			selected.forEach(x => editingSketch.removeElement(x))
			paintScreen()
		}
	}
	gl.onmouseup.push(function (e) {
		// don't put modeGetCurrent().mouseup in local var as 'this' won't be bound correctly
		modeGetCurrent().mouseup && modeGetCurrent().mouseup(e, getMouseLine(e))
		paintScreen()
	})
	gl.onmousedown.push(function (e) {
		if (1 == e.button) {
			modePop()
		}
		if (0 == e.button) {
			const mouseLine = getMouseLine(e)
			console.log("mouseLine", getMouseLine(e).toString(x => x), "mode", modeGetCurrent())
			modeGetCurrent().mousedown(e, mouseLine)
		}

		return false
	})
	$("clearmode").addEvent('click', modePop)

}
function initPointInfoEvents() {
	gl.onmousemove.push(function (e) {
		const mouseLine = getMouseLine({x: e.clientX, y: e.clientY})
		modeGetCurrent().mousemove && modeGetCurrent().mousemove(e, mouseLine)
		{
			let pp, html = '', closestP = Infinity
			drPs.forEach(info => {
				let p = info.p || info
				let text = p.toString(x => NLA.round10(x, -4)) + (info.p ? ' ' + info.text : '')
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
				let pSC = gl.projectionMatrix.times(gl.modelViewMatrix).transformPoint(pp)
				let x = (pSC.x * 0.5 + 0.5) * window.innerWidth, y = (-pSC.y * 0.5 + 0.5) * window.innerHeight
				tooltipShow(html, x, y)
				console.log("show tt")
			} else {
				tooltipHide()
			}
		}
		paintScreen()
	})
}
function initNavigationEvents() {
	const canvas: HTMLCanvasElement = $(gl.canvas)

	gl.onmousemove.push(function (e) {

		if (e.dragging) {

			//noinspection JSBitwiseOperatorUsage
			if (e.buttons & 4) {
				// pan
				const moveCamera = V(-e.deltaX, e.deltaY, 0).times(2 / gl.canvas.width)
				const inverseProjectionMatrix = gl.projectionMatrix.inversed()
				const worldMoveCamera = inverseProjectionMatrix.transformVector(moveCamera)
				eyePos = eyePos.plus(worldMoveCamera)
				eyeFocus = eyeFocus.plus(worldMoveCamera)
				setupCamera()
			}
			// scene rotation
			//noinspection JSBitwiseOperatorUsage
			if (e.buttons & 2) {
				let rotateLR = -e.deltaX / 6.0 * DEG
				let rotateUD = -e.deltaY / 6.0 * DEG
				// rotate
				let matrix = M4.rotationLine(eyeFocus, eyeUp, rotateLR)
				//let horizontalRotationAxis = eyeFocus.minus(eyePos).cross(eyeUp)
				let horizontalRotationAxis = eyeUp.cross(eyePos.minus(eyeFocus))
				matrix = matrix.times(M4.rotationLine(eyeFocus, horizontalRotationAxis, rotateUD))
				eyePos = matrix.transformPoint(eyePos)
				eyeUp = matrix.transformVector(eyeUp)

				setupCamera()
			}
		}
		paintScreen()
	})
	canvas.addEvent('mousewheel', function (e) {
		//console.log(e)
		zoomFactor *= pow(0.9, -e.wheel)
		const mouseCoords = e.client
		const moveCamera = V(mouseCoords.x * 2 / gl.canvas.width - 1, -mouseCoords.y * 2 / gl.canvas.height + 1, 0).times(1 - 1 / pow(0.9, -e.wheel))
		const inverseProjectionMatrix = gl.projectionMatrix.inversed()
		const worldMoveCamera = inverseProjectionMatrix.transformVector(moveCamera)
		//console.log("moveCamera", moveCamera)
		//console.log("worldMoveCamera", worldMoveCamera)
		eyePos = eyePos.plus(worldMoveCamera)
		eyeFocus = eyeFocus.plus(worldMoveCamera)
		setupCamera()
		paintScreen()
	})
}
//const sFace = B2T.rotateEdges([Edge.forCurveAndTs(EllipseCurve.XY, 0, 90 * DEG).rotateX(90 * DEG),StraightEdge.throughPoints(V3.Z, V3.X)], 45 * DEG, 'blah').faces.find(face => face.surface instanceof EllipsoidSurface)
//const face2 = B2T.extrudeEdges([Edge.forCurveAndTs(EllipseCurve.forAB(1, -1), -PI, 0), StraightEdge.throughPoints(V3.X, V3.X.negated())], P3.XY.flipped(), V3.Z, 'cyl')
//const cylface = cyl.faces.find(face => face instanceof RotationFace)//.rotateX(50 * DEG)
//assert(cylface.surface.facesOutwards())
//const cyl = B2T.extrudeEdges([Edge.forCurveAndTs(EllipseCurve.forAB(1, -1), -PI, 0), StraightEdge.throughPoints(V3.X, V3.X.negated())], P3.XY.flipped(), V3.Z, 'cyl')
window.onload = function () {
    const rt = RangeTree.fromArray([-2, -1, 0,2,3,4,7-1,7])
    rt.addIntervals([
        {left: -1, right: 6},
        {left: -2, right: 3},
        {left: 0, right: 4},
        {left: 2, right: 7},
        {left: 3, right: 4},
        {left: -2, right: -1}])
    console.log(rt.str)

    return
	gl = GL.create({canvas: document.getElementById('mainCanvas') as HTMLCanvasElement})

	const sphereSlice = B2T.rotateEdges([Edge.forCurveAndTs(EllipseCurve.XY, 0, 90 * DEG).rotateX(90 * DEG),StraightEdge.throughPoints(V3.Z, V3.X)], 45 * DEG, 'blah')

	//const cyl2 = B2T.extrudeEdges([Edge.forCurveAndTs(EllipseCurve.XY, PI, -PI)], P3.XY.flipped(), V3.Z, 'cyl')

	const plane = new P3(V3.sphere(-170 * DEG, 0), cos(10 * DEG))
	const isc = SemiEllipsoidSurface.unitISCurvesWithPlane(plane)[0]
	console.log(isc, SemiEllipsoidSurface.UNIT.containsCurve(isc))

	const face2 = new RotationFace(SemiEllipsoidSurface.UNIT, [
		Edge.forCurveAndTs(SemiEllipsoidSurface.unitISCurvesWithPlane(new P3(V3.sphere(-170 * DEG, 0), cos(10 * DEG)))[0], -PI, PI)])
	const face = sphereSlice.faces[0]
	//console.log("asd", cylface.toMesh().calcVolume().area)
	//console.log("asd", cylface.calcArea())
	return
	//console.log("ASHD", ConicSurface.XY.calculateArea([
	//	Edge.forCurveAndTs(EllipseCurve.XY.translate(0, 0, 1).shearedX(0, 1), -PI, +PI)
	//]))
	//console.log(ConicSurface.XY.shearedX(0, 1).toMesh(0, 1).calcVolume().area)
	////console.log(BezierCurve.EX2D.getAreaInDirSurface(V3.Y, new PlaneSurface(P3.XY), 0, 1))
	//console.log(M4.rotationLine(V3.O, V3.X, 3).str)
	//console.log(M4.rotationLine(V3.O, V3.X.times(2), 3).str)
	//console.log("4/3 PI = " + 4 / 3 * PI)
	//const cxz = new EllipseCurve(V3.O, V3.X, V3.Z), cxy = new EllipseCurve(V3.O, V3.X, V3.Y)
	//const loop = [Edge.forCurveAndTs(cxz, -PI/2, PI/2), Edge.forCurveAndTs(cxz, PI/2, -PI/2)]
	//console.log('EllipsoidSurface.XY.zDirVolumeForLoop(loop)', EllipsoidSurface.XY.zDirVolumeForLoop(loop))
	////console.log('EllipsoidSurface.XY.zDirVolumeForLoop(loop)', EllipsoidSurface.XY.zDirVolumeForLoop([
	////	Edge.forCurveAndTs(cxz, 0, PI/2),
	////	Edge.forCurveAndTs(cxz, PI/2, 0),
	////	Edge.forCurveAndTs(cxy, 2 * PI, 0),
	////]))
	//return
	////const ec = EllipseCurve.forAB(2.23, 3.05)
	//const ec = EllipseCurve.forAB(3.05, 2.23)
	//const f = t => ec.tangentAt(t).length(), startT = 0, endT = ec.angleToT(50 * DEG), steps = 0
	//console.log(PI)
	//console.log(99.54182501646524863791320)
	//console.log(gaussLegendreQuadrature24(f, Math.floor(endT / (PI / 2)) * PI / 2, endT) + Math.floor(endT / (PI / 2)) * gaussLegendreQuadrature24(f, 0, PI / 2))
	//console.log(glqInSteps(f, startT, endT, 8))
	//console.log(integrateCurve(ec, startT, endT, steps || 1024))
	//console.log(midpointRuleQuadrature(f, startT, endT, steps || 10024))
//	// goal: group linearly dependent vectors together
//	// these represent constraints which conflict
//	// are these groups clearly defined?
//	new Float32Array(3) instanceof MimeTypeArray
//	let m = Matrix.fromRowArrays(
//		[2, 3, 0, 0, 0],
//		[1, 0, 0, 0, 0],
//		[0, 0, 1, 2, 3],
//		[0, 1, 0, 0, 0],
//		[0, 0, 1, 2, 3],
//		[0, 0, 0, 0, 3],
//		[0, 0, 0, 0, 0]
//	)
//	//m = new M4(
//	//	2, -3, -3, 4,
//	//	0, 0, 1, -1,
//	//	0, 0, 0, 1,
//	//	0, 0, 0, 0
//	//)
//	let {L, U, P} = m.gauss()
//	const height = L.height
//	let dependents = new Array(height)
//	let uRowIndex = height
//	while (uRowIndex--) {
//		let uRow = U.row(uRowIndex)
//		if (uRow.length() < NLA_PRECISION) {
//			dependents[uRowIndex] = true
//		} else {
//			break
//		}
//	}
//	let lRowIndex = height
//	while (lRowIndex--) {
//		if (dependents[lRowIndex]) {
//			let lColIndex = Math.min(lRowIndex, L.width)
//			while (lColIndex--) {
//				if (0 !== L.e(lRowIndex, lColIndex)) {
//					dependents[lColIndex] = true
//				}
//			}
//		}
//	}
//	let indexMap = P.permutationAsIndexMap()
//	console.log(dependents, dependents.map((b, index) => b && index).filter(x => x != void 0))
//	console.log(dependents, dependents.map((b, index) => b && indexMap[index]).filter(x => x != void 0))
//	console.log(P.permutationAsIndexMap())
//	console.log("m\n", m.str)
//	console.log("L\n", L.str)
//	console.log("U\n", U.str)
//	console.log("P\n", P.str)
//return
	// initTips()
	// new Request({url: 'src/testshader.glsl', async: false, onSuccess: text => console.log(text)}).send()

	modePush(MODES.DEFAULT)

	$$('.sketchControl').set('disabled', true)

	gl = GL.create({canvas: document.getElementById('mainCanvas') as HTMLCanvasElement})
	gl.fullscreen()
	gl.canvas.oncontextmenu = () => false

	setupCamera()
	//gl.cullFace(gl.FRONT_AND_BACK);
	gl.clearColor(1.0, 1.0, 1.0, 0.0)
	gl.enable(gl.BLEND)
	gl.enable(gl.DEPTH_TEST)
	// gl.enable(gl.CULL_FACE)
	gl.depthFunc(gl.LEQUAL)
	gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA) // TODO ?!


	// let s1 = new ProjectedCurveSurface(BezierCurve.EX3D.scale(100, 100, 100).project(P3.XY), V3.Z, -1, 2)
	// let s2 = new ProjectedCurveSurface(BezierCurve.EX2D.scale(100, 100, 100), V3.Z, -1,
	// 2).rotateY(Math.PI/2).translate(0, 260, 0)
	mesh1 = GL.Mesh.sphere2(8, 8).scale(100, 100, 100)
	mesh1.computeWireframe()
	// mesh2 = s2.toMesh()
	// let isc = s2.isCurvesWithSurface(s1)[0]
	// mesh1 = face.inB2().toMesh()
	// mesh1.computeWireframeFromFlatTriangles()
	// mesh1.compile()
	// let mint = -2, maxt = 2
	// drPs.push(curve.at(mint), curve.at(maxt))
	// let line = L3(V3(-1560.8950828838565, 716.07295580975, 249.61382611323648), V(0.9130103135570956,
	// -0.36545647611595106, -0.18125598308272678)) isc.debugToMesh(mesh1, 'curve1') console.log(mesh1)
	mesh1.compile()
	// drPs.pushAll(isc.rootPoints().concatenated())
	// mesh2 = curve.getAABB(mint, maxt).toMesh()
	// console.log(curve.getAABB().toSource())
	// console.log(mesh2)
	//mesh1 = mesh1.transform(M4.FOO.as3x3())

	let linksHaendig = M4.forSys(V3.X, V3.Y, V3.Z.negated())
	let v0 = V(2, 3, 4), v1 = V(-2, 3, 5), v2 = v0.cross(v1)
	let v2t = linksHaendig.transformVector(v0).cross(linksHaendig.transformVector(v1))
	let v2s = linksHaendig.inversed().transposed().transformVector(v2)
	console.log("ASKDKJALDS", v2t.dot(v2s), v2t.isParallelTo(v2s))
	initMeshes()
	initShaders()
	let b = editingSketch && editingSketch.elements.find(el => el instanceof SketchBezier).toBrepEdge().curve
	//	console.log(mesh.vertices)
	console.log("BBB", b)
	initNavigationEvents()
	initOtherEvents()


	initLoadSave()

	mesh1.compile()
	if (window.location.hash && window.location.hash.length > 1) {
		const key = window.location.hash.substr(1)
		console.log(key)
		load(key)
	} else {
		//initModel()
	}

	rebuildModel() // necessary to init planes
	updateFeatureDisplay()
	rebuildModel() // so warning will show
	paintScreen()
	let lastInterst = featureStack.slice().reverse().find(f => f instanceof Sketch)
	lastInterst && modePush(MODES.SKETCH, lastInterst)

}

function getHovering(mouseLine: L3, ...consider: ('faces' | 'planes' | 'sketchElements' | 'brepPoints' | 'brepEdges')[]) {
	let hoverHighlight = null, nearest = Infinity

	function checkEl(el, distance) {
		if (distance < nearest) {
			nearest = distance
			hoverHighlight = el
		}
	}

	if (consider.includes('faces') && modelBREP) {
		modelBREP.faces.forEach((face) => {
			checkEl(face, face.intersectsLine(mouseLine))
		})
	}
	if (consider.includes('planes')) {
		planes.forEach(plane => checkEl(plane, plane.distanceTo(mouseLine)))
	}
	if (consider.includes('sketchElements')) {
		featureStack.filter(f => f instanceof Sketch).forEach(sketch => {
			if (!sketch.hide && sketch.plane && sketch.plane.normal.dot(mouseLine.dir1) < -0.1) {
				// sketch plane is facing user; ensures there is an intersection

				let mouseLineIS = mouseLine.intersectionWithPlane(sketch.plane)
				let sketchCoords = sketch.worldToSketchMatrix.transformPoint(mouseLineIS)

				let closestElement = sketch.elements.concat(getAllPoints(sketch).map(p => p.canon()).unique())
					.withMax(el => {
						let d = el.distanceToCoords(sketchCoords)
						if (el instanceof SegmentEndPoint) d -= 8
						return -d
					})

				if (closestElement && closestElement.distanceToCoords(sketchCoords) < 16) {
					// subtract 0.001 so that sketch elements have priority over things in same plane
					checkEl(closestElement, mouseLine.pointT(mouseLineIS) - 0.1)
				}
			}
		})
	}
	if (consider.includes('brepPoints')) {
		brepPoints.forEach(p => {
			let t = mouseLine.pointT(p)
			if (mouseLine.at(t).distanceTo(p) < 20) {
				checkEl(p, t - 0.1)
			}
		})
	}
	if (consider.includes('brepEdges')) {
		let projPlane = new P3(mouseLine.dir1, 0)
		let projPoint = projPlane.projectedPoint(mouseLine.anchor)
		brepEdges.forEach(edge => {
			let curve = edge.curve
			const prio = 0.05
			if (curve instanceof L3) {
				if (curve.dir1.isParallelTo(mouseLine.dir1)) {
					let d = mouseLine.distanceToPoint(edge.a)
					let t = mouseLine.pointT(edge.a)

					if (d < 16) {
						checkEl(edge, t - prio)
					}
				} else {
					let projAnchor = projPlane.projectedPoint(curve.anchor)
					let projDir = projPlane.projectedVector(curve.dir1)
					let tCurve = projPoint.minus(projAnchor).dot(projDir) / projDir.squared()
					tCurve = edge.clampedT(tCurve)
					let p = curve.at(tCurve)
					let t = mouseLine.pointT(p)
					if (mouseLine.at(t).distanceTo(p) < 16) {
						checkEl(edge, t - prio)
					}
				}
			} else {
				let projCurve = curve.project(projPlane)
				let tCurve = projCurve.closestTToPoint(projPoint)
				tCurve = edge.clampedT(tCurve)
				let p = curve.at(tCurve)
				let t = mouseLine.pointT(p)
				if (projCurve.at(tCurve).distanceTo(projPoint) < 16) {
					checkEl(edge, t - prio)
				}
			}
		})
	}

	return hoverHighlight
}
/**
 * Transforms mouse positions on the screen into a line in world coordinates.
 * @param {{x: number, y: number}} pos
 * @returns {L3}
 */
function getMouseLine(pos) {
	let ndc1 = V(pos.x * 2 / gl.canvas.width - 1, -pos.y * 2 / gl.canvas.height + 1, 0)
	let ndc2 = V(pos.x * 2 / gl.canvas.width - 1, -pos.y * 2 / gl.canvas.height + 1, 1)
	//console.log(ndc)
	let inverseProjectionMatrix = gl.projectionMatrix.inversed()
	let s = inverseProjectionMatrix.transformPoint(ndc1)
	let dir = inverseProjectionMatrix.transformPoint(ndc2).minus(s)
	return L3.anchorDirection(s, dir)
}


function setupSelectors(el, feature, mode) {

	el.getElements('.face-select, .plane-select, .segment-select')
		.removeEvents()
		.removeClass('selecting')
		.addEvents(tooltipEvents)
		.each(el => {
			el.set('text', feature[el.dataset.featureProperty].ref)
			el.linkRef = feature[el.dataset.featureProperty]
		})
		.addEvent('mouseover', function (e) {
			let target = this.linkRef.get()
			if (target) {
				hoverHighlight = target
			} else {
				missingEls.push(this.linkRef.lastHit)
			}
			paintScreen()
		})
		.addEvent('mouseout', function (e) {
			hoverHighlight = null
			missingEls.remove(this.linkRef.lastHit)
			paintScreen()
		})

	el.getElements('.face-select')
		.addEvent('click', function (e) {
			const selector = this
			this.addClass('selecting')
			selector.set('text', 'Click on a face')
			modePush(MODES.SELECT_FACE, function (face) {
				selector.removeClass('selecting')
				selector.set('text', face.name)

				feature[selector.dataset.featureProperty] = face.name
				rebuildModel()
			})
		})

	el.getElements('.plane-select')
		.addEvent('click', function (e) {
			const selector = this
			this.addClass('selecting')
			selector.set('text', 'Click on a plane')
			modePush(MODES.PLANE_SELECT, function (planeRef) {
				console.log("plane-select callback", planeRef)
				selector.removeClass('selecting')
				selector.set('text', planeRef.ref)
				selector.linkRef = planeRef
				modeEnd(MODES.PLANE_SELECT)

				feature[selector.dataset.featureProperty] = planeRef
				rebuildModel()
			})
		})

	el.getElements('.segment-select')
		.addEvent('click', function (e) {
			const selector = this
			this.addClass('selecting')
			selector.set('text', 'Click on a sketch segment')
			modePush(MODES.SELECT_SEGMENT, function (segmentRef) {
				selector.removeClass('selecting')
				selector.set('text', segmentRef.ref)
				selector.linkRef = segmentRef

				feature[selector.dataset.featureProperty] = segmentRef
				rebuildModel()
			})
		})

	el.getElements('.string-id-input, .select-text')
		.removeEvents()
		.each(el => el.value = feature[el.dataset.featureProperty])
		.addEvent('change', function (e) {
			// TODO: check if unique
			let propName = el.dataset.featureProperty
			feature[propName] = this.value
			rebuildModel()
		})

	el.getElements('.select-select')
		.removeEvents()
		.each(el => el.value = feature[el.dataset.featureProperty])
		.addEvent('change', function (e) {
			var propName = this.dataset.featureProperty
			feature[propName] = this.value
			rebuildModel()
		})

	el.getElements('.dimension-input')
		.removeEvents()
		.each(el => el.value = feature[el.dataset.featureProperty])
		.addEvent('mousewheel', function (e) {
			console.log(e)
			let delta = (e.shift ? 1 : 10) * Math.sign(e.wheel)
			this.set('value', NLA.round10(parseFloat(this.value) + delta, -6))
			this.fireEvent('change')
		})
		.addEvent('click', function (e) {
			if (e.event.button == 1) {
				this.set('value', this.dataset.defaultValue)
				this.fireEvent('change')
			}
		})
		.addEvent('change', function (e) {
			var propName = this.dataset.featureProperty
			feature[propName] = this.value = NLA.forceFinite(this.value)
			rebuildModel()
		})

	el.getElements('.boolean-input')
		.removeEvents()
		.each(el => el.checked = feature[el.dataset.featureProperty])
		.addEvent('change', function (e) {
			var propName = this.dataset.featureProperty
			feature[propName] = this.checked
			rebuildModel()
		})

	el.getElement('[name=done]')
		.removeEvents()
		.addEvent('click', function () {
			modeEnd(mode)
		})
	el.getElement('[name=delete]')
		.removeEvents()
		.addEvent('click', function () {
			featureDelete(feature)
		})

}
function featureDelete(feature) {
	let correspondingMode = modeStack.find(mode => mode.feature && mode.feature == feature)
	correspondingMode && modeEnd(correspondingMode)

	let featureIndex = featureStack.indexOf(feature)
	featureStack.remove(feature)
	if (-1 != featureIndex && featureIndex < rebuildLimit) {
		rebuildLimit--
	}
	updateFeatureDisplay()
	rebuildModel()
}
function featureDependents(feature) {
	let featureIndex = featureStack.indexOf(feature)
	return featureStack
		.slice(featureIndex + 1)
		.filter(followingFeature => featureDependencies(followingFeature).includes(feature))
}
function featureDependencies(feature) {
	let result = new Set()
	let featureNameDependecies = feature.dependentOnNames()
	featureNameDependecies.forEach(nameRef => {
		let dependentOnFeature = namesPublishedBy.get(nameRef.ref)
		dependentOnFeature && result.add(dependentOnFeature)
	})
	return Array.from(result.values())
}
function chainComparison(diff1, diff2) {
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
	for (var j = 0; j < 1; j++) {
		editingSketch.gaussNewtonStep()
	}
	editingSketch.reverse()
	paintScreen()
}
function deleteConstraint(sketch, constraint) {
	if ("coincident" == constraint.type) {
		deleteCoincidence(constraint, sketch)
	} else {
		sketch.constraints.remove(constraint)
	}
	sketch.recalculate()
	paintScreen()
}
function getGroupConstraint(el, sketch, type) {
	return sketch.constraints.find(function (c) {
		return c.type == type && (c.fixed == el || c.segments.includes(el))
	})
}
// removes segment or plane from group constraint
function removeFromGroupConstraint(el, sketch, type) {
	var groupConstraint = getGroupConstraint(el, sketch, type)
	if (!groupConstraint) return
	groupConstraint.cs.remove(el)
	if (1 == groupConstraint.cs.length) {
		sketch.constraints.remove(groupConstraint)
	}
	groupConstraint.segments.remove(el)
	if (groupConstraint.fixed == el) {
		groupConstraint.fixed = undefined
	}
}
function removeFromConstraint(el: SketchSegment | SegmentEndPoint, sketch: Sketch, constraint: Constraint) {
	switch (constraint.type) {
		case 'coincident':
			removeFromCoincidence(el, sketch)
			break
		case 'parallel':
		case 'colinear':
		case 'equalLength':
			removeFromGroupConstraint(el, sketch, constraint.type)
			break
		case 'perpendicular':
		case 'pointDistance':
		case 'pointLineDistance':
		case 'pointOnLine':
		case 'angle':
			sketch.constraints.remove(constraint)
			break
		default:
			throw new Error('unknown constraint ' + constraint.type)
		// console.log('errror', constraint.type);
	}
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
	var els = selected.filter((el) => el instanceof SketchLineSeg || ('equalLength' != type && isFixed(el)))
	if (els.length < 2) return

	var newGroup = els.map(el => {
		if (isFixed(el)) el = NameRef.forObject(el)
		var c = getGroupConstraint(el, editingSketch, type)
		return c ? c.cs : el
	}).concatenated().unique()
	var fixeds = newGroup.filter(el => el instanceof NameRef)
	if (1 < fixeds.length) {
		throw new Error("cannot have two fixed")
	}

	var oldConstraints = newGroup.map(el => getGroupConstraint(el, editingSketch, type))
	editingSketch.constraints.removeAll(oldConstraints)

	// add new constraint
	// fixeds[0] may be null
	var segments = newGroup.filter(x => !(x instanceof NameRef))
	var f = fixeds[0]
	editingSketch.constraints.push(new Constraint(type, f ? segments.concat(f) : segments, {
		fixed: f,
		segments: segments
	}))

	rebuildModel()
	paintScreen()
	updateSelected()
}


function makeAngle() {
	var selSegments = selected.filter((el) => el instanceof SketchLineSeg || el.plane)
	console.log(selSegments)
	if (selSegments.length == 2) {
		var segment = selSegments[0] instanceof SketchLineSeg ? selSegments[0] : selSegments[1]
		var other = selSegments[0] instanceof SketchLineSeg ? selSegments[1] : selSegments[0]
		if (!(segment instanceof SketchLineSeg)) {
			throw new Error("at least one must be a segment")
		}
		editingSketch.constraints.push(new Constraint("angle", [segment, other], {
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
	var selSegments = selected.filter((el) => el instanceof SketchLineSeg || el.plane)
	if (selSegments.length == 2) {
		var segment = selSegments[0] instanceof SketchLineSeg ? selSegments[0] : selSegments[1]
		var other = selSegments[0] instanceof SketchLineSeg ? selSegments[1] : selSegments[0]
		if (!(segment instanceof SketchLineSeg)) {
			throw new Error("at least one must be a segment")
		}
		editingSketch.constraints.push(new Constraint("perpendicular", [segment, other], {
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
	var point = selected[0] instanceof SegmentEndPoint ? selected[0] : selected[1]
	var other = selected[0] instanceof SegmentEndPoint ? selected[1] : selected[0]
	if (point instanceof SegmentEndPoint && (other instanceof SketchLineSeg || other instanceof SegmentEndPoint || other.plane)) {
		var newConstraint
		if (other instanceof SegmentEndPoint) {
			newConstraint = new Constraint("pointDistance", selected.slice(), {distance: round(other.distanceToCoords(point))})
		} else if (other instanceof SketchLineSeg) {
			let distance = round(other.distanceToCoords(point))
			newConstraint = new Constraint("pointLineDistance", [point, other],
				{point: point, other: other, distance: distance})
		} else {
			let distance = round(other.plane.intersectionWithPlane(editingSketch.plane)
				.transform(editingSketch.worldToSketchMatrix).distanceToPoint(point.V3()))
			other = NameRef.forObject(other)
			newConstraint = new Constraint("pointPlaneDistance", [point, other],
				{point: point, other: other, distance: distance})
			console.log(newConstraint)
		}
		editingSketch.constraints.push(newConstraint)
		rebuildModel()
		paintScreen()
		updateSelected()
		;
		($("distanceInput" + newConstraint.id) as HTMLInputElement).select()
	}
}
interface getCurvable {
	getCurve()
}
interface Edge extends getCurvable {}
Edge.prototype.getCurve = function () {

}
function selPointOnLine() {
	if (2 != selected.length) return
	var point = selected[0] instanceof SegmentEndPoint ? selected[0] : selected[1]
	var other = selected[0] instanceof SegmentEndPoint ? selected[1] : selected[0]
	if (point instanceof SegmentEndPoint && (other instanceof SketchLineSeg || other instanceof SketchBezier || other instanceof SketchArc
		|| other instanceof Face || other instanceof Edge || other instanceof CustomPlane)) {
		other = isSketchEl(other) ? other : NameRef.forObject(other)
		var newConstraint = new Constraint("pointOnLine", [point, other], {point: point, other: other})
		console.log(newConstraint)
		editingSketch.constraints.push(newConstraint)
		rebuildModel()
		paintScreen()
		updateSelected()
	}
}
function makeSelCoincident() {
	var selPoints = selected.filter(function (s) { return s instanceof SegmentEndPoint })
	if (selPoints.length >= 2) {
		for (var i = 1; i < selPoints.length; i++) {
			makeCoincident(selPoints[0], selPoints[i], editingSketch)
		}
		rebuildModel()
		paintScreen()
		updateSelected()
	}
}

function faceSketchPlane() {
	var plane = editingSketch.plane
	var viewLine = L3.throughPoints(eyePos, eyeFocus)
	eyeFocus = plane.intersectionWithLine(viewLine) || viewLine.closestPointToPoint(V3.O)
	eyePos = eyeFocus.plus(plane.normal.times(100))
	eyeUp = plane.up
	setupCamera()
	paintScreen()
}
class Extrude {
	type: string
	name: string
	_segmentName: NameRef
	start: number
	end: number
	operation: 'minus' | 'plus'

	constructor(name, segmentName, start, end, operation) {
		this.type = "extrude"
		this.name = name || "extrude" + (globalId++)
		segmentName = segmentName || NameRef.UNASSIGNED
		assertInst(NameRef, segmentName)
		this.segmentName = segmentName
		this.start = start || 0
		this.end = end || 100
		this.operation = operation || "minus"
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
		return [this.segmentName]
	}
}
NLA.registerClass(Extrude)

class Pattern {
	name = 'pattern' + globalId++
	features
	direction = V3.X
	count = 2
	totalLength = 20
	intervalLength = 10
	mode = 0

	build() {

	}
}


function modeGetName(mode) {
	return Object.getOwnPropertyNames(MODES).find(name => MODES[name] == mode)
}
function modeUpdateDisplay() {
	$('modeBox').set('html', modeStack.map(modeGetName).join('<br>'))
}
function modePush(mode, ...args) {
	assert(mode.init && mode.end && mode.mousemove && mode.mousedown)
	mode.before && mode.before()
	let feature = args[0]
	if (mode.modeEditsFeatureAndRequiresRollback && feature) {
		rebuildLimit = featureStack.indexOf(feature) + 1
		updateFeatureDisplay()
		rebuildModel()
	}
	modeStack.push(mode)
	mode.init.apply(mode, args)
	modeUpdateDisplay()
}
function modeEnd(mode) {
	if (modeStack.includes(mode)) {
		do {
			var popped = modeStack.pop()
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

interface Mode {
	end()
	init()
	mousemove()
	mousedown()

}


function tooltipShow(htmlContent, x, y) {
	let tipWrap = $('tip-wrap') as HTMLDivElement
	assert(tipWrap)
	let tip = $('tip') as HTMLDivElement
	tipWrap.setStyle('visibility', 'visible')
	tip.set('html', htmlContent)
	let size = tipWrap.getSize()
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

let tooltipEvents = {}


function serialize(v) {
	function gatherList(v) {
		if (v && v.constructor === Array) {
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
				let keys = Object.keys(v)
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
				let result = Object.prototype == v.prototype ? {}
					: assert(v.constructor && v.constructor.name,
					() => (console.log(v), v.toSource() + v.constructor.name)) && {'#PROTO': v.constructor.name}
				let keys = Object.keys(v)
				for (let i = 0; i < keys.length; i++) {
					result[keys[i]] = transform(v[keys[i]])
				}
				return result
			}
		} else {
			throw new Error('?' + typeof v + v.toString())
		}
	}

	let visited = new Set()
	let listMap = new Map(), resultList = []
	listMap.set(v, 0)
	resultList.push(v)
	gatherList(v)
	console.log(resultList)

	resultList = resultList.map(v => transform(v, true))
	console.log(JSON.stringify(resultList))
	return JSON.stringify(resultList)
}

function unserialize(string) {
	function fixObjects(v) {
		if (v && v.constructor === Array) {
			for (let i = 0; i < v.length; i++) {
				v[i] = fixObjects(v[i])
			}
			return v
		} else if ('object' == typeof v && null != v) {
			if ('#PROTO' in v) {
				assert(window[v['#PROTO']] || NLA.CLASSES[v['#PROTO']], v['#PROTO'] + ' Missing ' + window[v['#PROTO']])
				let result = Object.create((window[v['#PROTO']] || NLA.CLASSES[v['#PROTO']]).prototype)
				let keys = Object.keys(v)
				for (let i = 0; i < keys.length; i++) {
					//if ('name' == keys[i]) console.log(result)
					if ('#PROTO' != keys[i]) {
						Object.defineProperty(result, keys[i], {
							value: fixObjects(v[keys[i]]),
							enumerable: true,
							writable: true
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
				let keys = Object.keys(v)
				for (let i = 0; i < keys.length; i++) {
					v[keys[i]] = linkReferences(v[keys[i]])
				}
				return v
			}
		} else {
			return v
		}
	}

	let tree = JSON.parse(string)
	// console.log(tree)
	fixObjects(tree)
	// console.log(tree)
	linkReferences(tree)
	// console.log(tree)
	return tree[0]
}
