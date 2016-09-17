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
		setupSelectors(div, feature)

		$$('.sketchControl').set('disabled', false)


		if (!feature.planeName) {
			div.getElement('[name=planeName]').fireEvent('click')
		}
	},
	end: function () {
		editingSketch = null
		$$('.sketchControl').set('disabled', true)
	},
	mouseup: function (/** MouseEvent */ e, mouseLine) {
		console.log('MOUSEUP', this.draggingElements, this)
		console.log(this.draggingElements, e.shiftKey)
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
		hoverHighlight = getHovering(mouseLine, 'sketchElements', 'planes', 'faces')

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
	deleteFeature: function () {
		featureStack.remove(this.feature)
		updateFeatureDisplay()
		rebuildModel()
		modeEnd(this)
	},
	mousedown: function (e, mouseLine) {
		// TODO: select sketch elements
		this.draggingElements = false
		this.relPoss = new Map()
		if(hoverHighlight != null) {
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
MODES.EXTRUDE = {
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
		if (!feature.segmentName) {
			div.getElement('[name=segmentName]').fireEvent('click')
		}
	},
	deleteFeature: function () {
		featureStack.remove(this.feature)
		updateFeatureDisplay()
		rebuildModel()
		modeEnd(this)
	},
	end: function () {
		let div = $('extrudeEditor')
		div.setStyle('display', 'none')
	},
	mousemove: function (e) {},
	mousedown: function (e) {}

}
MODES.PLANE_DEFINITION = {
	before: function () {
		modeEnd(MODES.SKETCH)
		modeEnd(MODES.PLANE_DEFINITION)
		modeEnd(MODES.EXTRUDE)
	},
	init: function (feature) {
		if (feature) {
			// assertInst(Extrude, feature)
		} else {
			feature = new PlaneDefinition("face")
			featureStack.push(feature)
			updateFeatureDisplay()
		}
		this.feature = feature
		var div = $("planeDefiner")
		div.setStyle("display", "block")
		setupSelectors(div, feature)
		if (feature.planeType == "face" && !feature.faceName) {
			div.getElement("[name=faceName]").fireEvent("click")
		}
	},
	end: function () {},
	mousemove: function (e) {},
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
		if (hoverHighlight instanceof CustomPlane) {
			this.callback(hoverHighlight)
		} else {
			// TODO: create plane on face
			let planeDefinition = new PlaneDefinition()
			planeDefinition.type = 'face'
			planeDefinition.faceName = hoverHighlight.name
			featureStack.push(planeDefinition)
			this.callback(planeDefinition)
		}
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
		this.callback(hoverHighlight)
	}
}
MODES.ADD_SEGMENT = {
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
		}
	}
}