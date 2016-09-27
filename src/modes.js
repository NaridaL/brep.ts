"use strict"
/**
 * Created by aval on 26.08.2016.
 */
var MODES = {}
MODES.DEFAULT = {
	init: function () {},
	end: function () {
		throw new Error("Can't end default mode")
	},
	mousemove: function (e) {},
	mousedown: function (e) {}
}
MODES.SKETCH = {
	modeEditsFeatureAndRequiresRollback: true,
	relPoss: null, // positions of sketch elements relative to mouse coords, used for DnD
	draggingElements: false,
	before: function () {
		modeEnd(MODES.SKETCH)
		modeEnd(MODES.PLANE_DEFINITION)
		modeEnd(MODES.EXTRUDE)
	},
	init: function (feature) {
		if (feature) {
			assertInst(Sketch, feature)
		} else {
			feature = new Sketch()
			featureStack.push(feature)
			updateFeatureDisplay()
		}

		this.feature = feature
		editingSketch = feature
		let div = $('sketchEditor')
		selected = []
		updateSelected()
		div.setStyle('display', 'block')
		setupSelectors(div, feature, this)



		$$('.sketchControl').set('disabled', false)


		if (NameRef.UNASSIGNED == feature.planeRef) {
			div.getElement('[data-feature-property=planeRef]').fireEvent('click')
		}
	},
	end: function () {
		editingSketch = null
		$$('.sketchControl').set('disabled', true)
		let div = $('sketchEditor')
		div.setStyle('display', 'none')
	},
	mouseup: function (/** MouseEvent */ e, mouseLine) {
		console.log('MOUSEUP', 'buttons', e.buttons, 'button', e.button, 'which', e.which)
		if (this.draggingElements) {
			console.log('draggingElements')
			rebuildModel()
		} else {
			if (e.shiftKey) {
				selected.toggle(hoverHighlight)
			} else if (hoverHighlight) {
				selected = [hoverHighlight]
			}
			updateSelected()
		}

	},
	mousemove: function (e, mouseLine) {
		// TODO: highlight sketch elements
		hoverHighlight = getHovering(mouseLine, 'sketchElements', 'planes', 'faces', 'brepPoints', 'brepEdges')

		// drag elements
		//noinspection JSBitwiseOperatorUsage
		if (e.buttons & 1) {
			let intersection = mouseLine.intersectionWithPlane(editingSketch.plane)

			if (intersection) {
				// parallel to mouseline
				let mouseSC = editingSketch.worldToSketchMatrix.transformPoint(intersection)
				for (var [point, relPos] of this.relPoss.entries()) {
					var newPos = mouseSC.plus(relPos)
					point.moveCoincidence(newPos)
				}
				this.draggingElements = true
				paintScreen();
			}
		}
	},
	mousedown: function (e, mouseLine) {
		console.log('mousedown', 'buttons', e.buttons, 'button', e.button, 'which', e.which)
		// TODO: select sketch elements
		this.draggingElements = false
		this.relPoss = new Map()
		if (hoverHighlight != null) {
			let intersection = mouseLine.intersectionWithPlane(editingSketch.plane)
			if (intersection && !GL.keys.SHIFT) {
				let mouseSC = editingSketch.worldToSketchMatrix.transformPoint(intersection)
				if (!selected.contains(hoverHighlight)) {
					selected = [hoverHighlight]
				}
				selected.filter(el => isSketchEl(el) && el.sketch == editingSketch)
					.map(el => el instanceof SegmentEndPoint ? el : el.points)
					.concatenated()
					.forEach(p => this.relPoss.set(p, p.V3().minus(mouseSC)))
			}
			updateSelected()
			console.log("selected is now ", selected);
			console.log("highlighted is now ", highlighted);
		}
	}

}
const BUTTONS = {LEFT: 0, MIDDLE: 1, RIGHT: 2, BACK: 3, FORWARD: 4}
MODES.EXTRUDE = {
	modeEditsFeatureAndRequiresRollback: true,
	feature: undefined,
	before: function () {
		modeEnd(MODES.SKETCH)
		modeEnd(MODES.PLANE_DEFINITION)
		modeEnd(MODES.EXTRUDE)
	},
	init: function (feature) {
		if (feature) {
			assertInst(Extrude, feature)
		} else {
			feature = new Extrude()
			featureStack.push(feature)
			updateFeatureDisplay()
		}
		this.feature = feature
		let div = $('extrudeEditor')
		div.setStyle('display', 'block')
		setupSelectors(div, feature, this)
		if (NameRef.UNASSIGNED == feature.segmentName) {
			div.getElement('[data-feature-property=segmentName]').fireEvent('click')
		}
	},
	end: function () {
		let div = $('extrudeEditor')
		div.setStyle('display', 'none')
	},
	mousemove: function (e) {},
	mousedown: function (e) {}

}
MODES.PLANE_DEFINITION = {
	modeEditsFeatureAndRequiresRollback: true,
	magic: function (sel, rads) {
		const CLASS_ORDER = [V3, Curve, Edge, Surface, Face]
		const sortMap = o => CLASS_ORDER.findIndex(clazz => o instanceof clazz)
		sel.sort((a, b) => sortMap(a) - sortMap(b))
		let [a, b, c] = sel
		if (sel.length == 3) {
			if (a instanceof V3 && b instanceof V3 && c instanceof V3) {
				return P3.throughPoints(a, b, c)
			}
		}
		if (2 == sel.length) {
			if (a instanceof V3 && b instanceof Face) {
				let p = /** @type V3 */ a
				let face = /** @type Face */ b
				if (face.containsPoint(p)) {
					return P3.normalOnAnchor(face.surface.normalAt(p), p)
				}
			}
			if (a instanceof V3 && (b instanceof Curve || b instanceof Edge)) {
				let p = /** @type V3 */ a
				let curve = /** @type Curve */ (b instanceof Curve ? b : b.curve)
				if (curve instanceof L3) {
					// a line and a point can define a plane one of two ways
					//  1. line as normal, point as anchor
					//  2. any two points on the line and the point define three points on the plane
					//     requires point not to be on line
					if (curve.containsPoint(p) /**|| 'pointAndNormal' == mode*/) {
						// mode = 'pointAndNormal'
						return P3.normalOnAnchor(curve.dir1, p)
					} else {
						let planeVector2 = M4.rotation(rads, curve.dir1, M4.temp0).transformVector(p.minus(curve.anchor))
						return P3.forAnchorAndPlaneVectors(curve.anchor, curve.dir1, planeVector2)
					}
				} else if (curve.containsPoint(p)) {
					// if the curve is not a line, the point needs to be on the curve so we can get a specific vector
					// TODO: closestTTOPoint instead, so we get a specific one in case of self intersection?
					return P3.normalOnAnchor(curve.tangentAt(curve.pointLambda(p)), p)
				}
			}
		}
		if (1 == sel.length) {
			if (a instanceof Face && a.surface instanceof PlaneSurface) {
				return a.surface.plane
			}
			if (a instanceof P3) {
				return a
			}
		}
	},
	before: function () {
		modeEnd(MODES.SKETCH)
		modeEnd(MODES.PLANE_DEFINITION)
		modeEnd(MODES.EXTRUDE)
	},
	init: function (feature) {
		if (feature) {
			assertInst(PlaneDefinition, feature)
			selected = feature.whats
		} else {
			feature = new PlaneDefinition('plane' + globalId++)
			featureStack.push(feature)
			updateFeatureDisplay()
		}
		updateSelected()
		this.feature = feature
		var div = $('planeDefiner')
		div.setStyle('display', 'block')
		setupSelectors(div, feature, this)
	},
	end: function () {
		var div = $('planeDefiner')
		div.setStyle('display', 'none')
	},
	mousemove: function (e, mouseLine) {
		hoverHighlight = getHovering(mouseLine, 'planes', 'faces', 'brepEdges', 'brepPoints')
	},
	mouseup: function (e) {
		console.log('e', e)
		if (BUTTONS.LEFT == e.button) {
			const nameRef = hoverHighlight && NameRef.forObject(hoverHighlight)
			if (e.shiftKey) {
				nameRef && selected.toggle(nameRef)
			} else {
				if (nameRef) {
					selected = [nameRef]
				} else {
					selected = []
				}
			}
			this.feature.whats = selected.slice()
			rebuildModel()
			updateSelected()
		}
	},
	mousedown: function (e) {}
}
MODES.PLANE_SELECT = {
	init: function (callback) {
		this.callback = callback
	},
	end: function () {},
	mousemove: function (e, mouseLine) {
		hoverHighlight = getHovering(mouseLine, 'planes', 'faces')
	},
	mousedown: function (e) {
		if (null == hoverHighlight) return
		let customPlane
		if (hoverHighlight instanceof CustomPlane) {
			customPlane = hoverHighlight
		} else {
			// TODO: create plane on face
			let planeDefinition = new PlaneDefinition()
			planeDefinition.type = 'face'
			planeDefinition.faceName = hoverHighlight.name
			featureStack.push(planeDefinition)
			this.callback(planeDefinition)
		}
		this.callback(NameRef.forObject(customPlane))
	}
}
MODES.SELECT_SEGMENT = {
	init: function (callback) {
		this.callback = callback
	},
	end: function () {},
	mousemove: function (e, mouseLine) {
		hoverHighlight = getHovering(mouseLine, 'sketchElements') // todo only sketch segments
		paintScreen()
	},
	mousedown: function (e) {
		this.callback(NameRef.forObject(hoverHighlight))
	}
}
MODES.ADD_SEGMENT = {
	before: function () {
		modeEnd(MODES.ADD_SEGMENT)
	},
	init: function (constructor) {
		this.currentAddingSegment = null
		this.constructor = constructor
	},
	end: function () {
		this.currentAddingSegment.points.forEach(p => removeFromCoincidence(p, editingSketch))
		editingSketch.elements.remove(this.currentAddingSegment);
		paintScreen();
	},
	mousemove: function (e, mouseLine) {
		let intersection = mouseLine.intersectionWithPlane(editingSketch.plane)
		if (intersection == null) {
			return;
		}
		if (!this.currentAddingSegment) {
			this.currentAddingSegment = new this.constructor(editingSketch)
			editingSketch.elements.push(this.currentAddingSegment)
			this.arcmode = 0
		}
		let mouseSC = editingSketch.worldToSketchMatrix.transformPoint(intersection)
		removeFromCoincidence(this.currentAddingSegment.points[this.arcmode], editingSketch)
		//console.log(mousePos);
		var points = getAllPoints(editingSketch);
		points.removeAll(this.currentAddingSegment.points);
		var pointDistances = points
			.map(function (point) { return {point: point, distance: point.distanceToCoords(mouseSC)} })
			.filter(function (pair) { return pair.distance < 16; })
			.sort(function (a, b) { return a.distance - b.distance; })
		if (pointDistances.length != 0 && this.currentAddingSegment) {
			makeCoincident(pointDistances[0].point, this.currentAddingSegment.points[this.arcmode], editingSketch)
			mouseSC = pointDistances[0].point;
		}
		for (let i = this.arcmode; i < this.currentAddingSegment.points.length; i++) {
			this.currentAddingSegment.points[i].x = mouseSC.x
			this.currentAddingSegment.points[i].y = mouseSC.y
		}
		/*
		 for (var i in lines) {
		 var line = lines[i];
		 //console.log("line", line);
		 if (line != this.currentAddingSegment && line.distanceTo(mousePos.x, mousePos.y) < 16) {
		 mousePos = line.getClosestPoint(mousePos.x, mousePos.y);
		 break;
		 }
		 };
			 */
	},
	mousedown: function (e, mouseLine) {
		let intersection = mouseLine.intersectionWithPlane(editingSketch.plane)
		let sketchCoords = editingSketch.worldToSketchMatrix.transformPoint(intersection)
		if (intersection == null) {
			return;
		}
		this.arcmode++
		if (this.arcmode == this.currentAddingSegment.points.length) {
			// finished adding the current segment
			this.currentAddingSegment = null
			if (SketchLineSeg == this.constructor) {
				this.mousemove(e, mouseLine)
				this.mousedown(e, mouseLine)
			}
		}
	}
}