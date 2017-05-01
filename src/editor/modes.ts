const MODES:any = {}
MODES.DEFAULT = {
	init: function () {},
	end: function () {
		throw new Error('Can\'t end default mode')
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
	mouseup: function (e, mouseLine) {
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
		hoverHighlight = getHovering(mouseLine, modelBREP && modelBREP.faces, planes, brepPoints, brepEdges, 16, 'sketchElements', 'planes', 'faces', 'points', 'edges')

		// drag elements
		//noinspection JSBitwiseOperatorUsage
		if (e.buttons & 1) {
			let intersection = mouseLine.intersectionWithPlane(editingSketch.plane)

			if (intersection) {
				// parallel to mouseline
				let mouseSC = editingSketch.worldToSketchMatrix.transformPoint(intersection)
				for (const [point, relPos] of this.relPoss.entries()) {
					const newPos = mouseSC.plus(relPos)
					point.moveCoincidence(newPos)
				}
				this.draggingElements = true
				paintScreen()
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
			if (intersection && !GL.KEYS.SHIFT) {
				let mouseSC = editingSketch.worldToSketchMatrix.transformPoint(intersection)
				if (!selected.includes(hoverHighlight)) {
					selected = [hoverHighlight]
				}
				selected.filter(el => isSketchEl(el) && editingSketch == el.sketch)
					.map(el => el instanceof SegmentEndPoint ? el : el.points)
					.concatenated()
					.forEach(p => this.relPoss.set(p, p.V3().minus(mouseSC)))
				console.log(this.relPoss)
			}
			updateSelected()
			console.log('selected is now ', selected)
			console.log('highlighted is now ', highlighted)
		}
	}

}
const BUTTONS = {LEFT: 0, MIDDLE: 1, RIGHT: 2, BACK: 3, FORWARD: 4}
function initFeatureMode(this: any, feature: Feature, featType: any, component) {
	if (feature) {
		assertInst(featType, feature)
	} else {
		feature = new featType()
		featureStack.push(feature)
		updateFeatureDisplay()
	}
	this.feature = feature

	const props = {
		feature,
		del: e => featureDelete(feature),
		done: e => modePop(),
		notifier: (feature, propName, value) => {
			feature[propName] = value
			preact.render(h(component, props), $('preact'), this.comp)
			rebuildModel()
		}}
	this.comp = preact.render(h(component, props), $('preact'))
	//if (NameRef.UNASSIGNED == feature.segmentName) {
	//	div.getElement('[data-feature-property=segmentName]').fireEvent('click')
	//}
}
MODES.ROTATE = {
	modeEditsFeatureAndRequiresRollback: true,
	feature: undefined,
	before: function () {
		modeEnd(MODES.SKETCH)
		modeEnd(MODES.PLANE_DEFINITION)
		modeEnd(MODES.EXTRUDE)
		modeEnd(MODES.ROTATE)
	},
	init: function (feature) {
		initFeatureMode.call(this, feature, Rotate, RotateEditor)
	},
	end: function () {
		preact.render('', $('preact'), this.comp)
	},
	mousemove: function (e) {},
	mousedown: function (e) {}
}
MODES.EXTRUDE = {
	modeEditsFeatureAndRequiresRollback: true,
	feature: undefined,
	before: function () {
		modeEnd(MODES.SKETCH)
		modeEnd(MODES.PLANE_DEFINITION)
		modeEnd(MODES.EXTRUDE)
	},
	init: function (feature) {
		initFeatureMode.call(this, feature, Extrude, ExtrudeEditor)
	},
	end: function () {
		preact.render('', $('preact'), this.comp)
	},
	mousemove: function (e) {},
	mousedown: function (e) {}

}
MODES.PATTERN = {
	modeEditsFeatureAndRequiresRollback: true,
	feature: undefined,
	before: function () {
		modeEnd(MODES.SKETCH)
		modeEnd(MODES.PLANE_DEFINITION)
		modeEnd(MODES.EXTRUDE)
		modeEnd(MODES.PATTERN)
	},
	init: function (feature) {
		if (feature) {
			assertInst(Pattern, feature)
		} else {
			feature = new Pattern()
			featureStack.push(feature)
			updateFeatureDisplay()
		}
		const div = $('patternEditor')
		div.setStyle('display', 'block')
		setupSelectors(div, feature, this)
	},
	end: function () {
		const div = $('patternEditor')
		//div.erase('text')
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
		const [a, b, c] = sel
		if (3 == sel.length) {
			if (a instanceof V3 && b instanceof V3 && c instanceof V3) {
				return P3.throughPoints(a, b, c)
			}
		}
		if (2 == sel.length) {
			if (a instanceof V3 && b instanceof Face) {
				if (b.containsPoint(a)) {
					return P3.normalOnAnchor(b.surface.normalAt(a), a)
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
					return P3.normalOnAnchor(curve.tangentAt(curve.pointT(p)), p)
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
		const div = $('planeDefiner')
		div.setStyle('display', 'block')
		setupSelectors(div, feature, this)
	},
	end: function () {
		const div = $('planeDefiner')
		div.setStyle('display', 'none')
	},
	mousemove: function (e, mouseLine) {
		hoverHighlight = getHovering(mouseLine, modelBREP.faces, planes, brepPoints, brepEdges, 16, 'planes', 'faces', 'edges', 'points')
	},
	mouseup: function (e) {
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
	init: function (callback: (NameRef) => void) {
		this.callback = callback
	},
	end: function () {},
	mousemove: function (e, mouseLine) {
		hoverHighlight = getHovering(mouseLine, modelBREP && modelBREP.faces, planes, brepPoints, brepEdges, 16, 'planes', 'faces')
	},
	mouseup: function (e) {
		if (null == hoverHighlight) return
		let customPlane
		if (hoverHighlight instanceof CustomPlane) {
			customPlane = hoverHighlight
		} else {
			// TODO: create plane on face
			const planeDefinition = new PlaneDefinition()
			planeDefinition.planeType = 'face'
			planeDefinition.whats = [NameRef.forObject(hoverHighlight)]
			featureStack.splice(featureStack.length - 1, 0, planeDefinition)
			customPlane = planeDefinition
		}
		this.callback(NameRef.forObject(customPlane))
	},
	mousedown: function (e) {},
}
MODES.SELECT_DIRECTION = {
	init: function (callback: (NameRef) => void) {
		this.callback = callback
	},
	magic: function (sel): [V3, string] {
		const r = (vector: V3, type: string): [V3, string] => [vector, type]
		const CLASS_ORDER = [V3, Edge, Surface, Face, P3]
		const sortMap = o => CLASS_ORDER.findIndex(clazz => o instanceof clazz)
		sel.sort((a, b) => sortMap(a) - sortMap(b))
		const [a, b, c] = sel
		if (3 == sel.length) {
			if (a instanceof V3 && b instanceof V3 && c instanceof V3) {
				return r(V3.normalOnPoints(a, b, c).unit(),
					'normal of implicit plane through three points')
			}
		}
		if (2 == sel.length) {
			if (a instanceof V3 && b instanceof V3) {
				return r(a.to(b).unit(), 'through two points')
			}
			if (a instanceof V3 &&
				b instanceof Edge &&
				!(b.curve instanceof L3) &&
				b.curve.containsPoint(a)) {
				const aT = b.curve.pointT(a)
				return r(b.curve.tangentAt(aT).unit(), 'tangent of edge at point')
			}
			if (a instanceof V3 &&
				b instanceof Face &&
				!(b.surface instanceof PlaneSurface) &&
				b.surface.containsPoint(a)) {
				return r(b.surface.normalAt(a).unit(), 'normal of surface at point')
			}
		}
		if (1 == sel.length) {
			if (a instanceof P3) {
				return r(a.normal1, 'plane normal')
			}
			if (a instanceof Edge && a.curve instanceof L3) {
				return r(a.curve.dir1, 'tangent of straight edge')
			}
			if (a instanceof Face && a.surface instanceof PlaneSurface) {
				return r(a.surface.plane.normal1, 'normal of planar surface')
			}
		}
	},
	end: function () {
		this.callback(selected.map(s => NameRef.forObject(s)))
	},
	mousemove: function (e, mouseLine) {
		hoverHighlight = getHovering(mouseLine, modelBREP && modelBREP.faces, planes, brepPoints, brepEdges, 16, 'planes', 'faces', 'points', 'edges')
	},
	mouseup: function (e) {
		if (BUTTONS.LEFT == e.button) {
			if (e.shiftKey) {
				selected.toggle(hoverHighlight)
			} else {
				if (hoverHighlight) {
					selected = [hoverHighlight]
				} else {
					selected = []
				}
			}
			this.callback(selected.map(s => NameRef.forObject(s)))
			updateSelected()
		}
	},
	mousedown: function (e, mouseLine) {},
}
MODES.SELECT_LINE = {
	init: function (callback: (NameRef) => void) {
		this.callback = callback
	},
	magic: function (sel): [L3, string] {
		const r = (vector: L3, type: string): [L3, string] => [vector, type]
		const CLASS_ORDER = [V3, Edge, Surface, Face, P3]
		const sortMap = o => CLASS_ORDER.findIndex(clazz => o instanceof clazz)
		sel.sort((a, b) => sortMap(a) - sortMap(b))
		const [a, b, c] = sel
		if (2 == sel.length) {
			if (a instanceof V3 && b instanceof V3) {
				return r(L3.throughPoints(a, b), 'through two points')
			}
			if (a instanceof V3 &&
				b instanceof Edge &&
				!(b.curve instanceof L3)) {
				const curve = b.curve
				if (curve.containsPoint(a)) {
					const aT = b.curve.pointT(a)
					return r(new L3(a, b.curve.tangentAt(aT).unit()), 'tangent of edge at point')
				} else {
					const t = b.curve.closestTToPoint(a)
					return r(new L3(a, b.curve.tangentAt(t).unit()), 'point and tangent at closest point on edge')
				}
			}
			if (a instanceof V3 &&
				b instanceof Face) {
				if (b.surface instanceof PlaneSurface) {
					const line = new L3(a, b.surface.plane.normal1)
					return r(line, 'point and plane normal')
				} else if (b.surface.containsPoint(a)) {
					const line = new L3(a, b.surface.normalAt(a))
					return r(line, 'point and surface normal')
				} else {
					const foot = b.surface.pointFoot(a)
					const line = L3.throughPoints(foot, a)
					return r(line, 'through point and surface point foot')
				}
			}
			if (a instanceof V3 && b instanceof P3) {
				return r(new L3(a, b.normal1), 'plane normal')
			}
		}
		if (1 == sel.length) {
			if (a instanceof Edge && a.curve instanceof L3) {
				return r(a.curve, 'straight edge')
			}
			if (a instanceof SketchLineSeg) {
				return r(a.getCurve(), 'sketch segment')
			}
			if (a instanceof L3) {
				return r(a, 'line')
			}
		}
	},
	end: function () {
		this.callback(selected.map(s => NameRef.forObject(s)))
	},
	mousemove: function (e, mouseLine) {
		hoverHighlight = getHovering(mouseLine, modelBREP && modelBREP.faces, planes, brepPoints, brepEdges, 16,
			'planes', 'faces', 'points', 'edges', 'sketchElements')
	},
	mouseup: function (e) {
		if (BUTTONS.LEFT == e.button) {
			if (e.shiftKey) {
				selected.toggle(hoverHighlight)
			} else {
				if (hoverHighlight) {
					selected = [hoverHighlight]
				} else {
					selected = []
				}
			}
			this.callback(selected.map(s => NameRef.forObject(s)))
			updateSelected()
		}
	},
	mousedown: function (e, mouseLine) {},
}
MODES.SELECT_SEGMENT = {
	init: function (callback) {
		this.callback = callback
	},
	end: function () {},
	mousemove: function (e, mouseLine) {
		hoverHighlight = getHovering(mouseLine, modelBREP && modelBREP.faces, planes, brepPoints, brepEdges, 16, 'sketchElements') // todo only sketch segments
		paintScreen()
	},
	mousedown: function (e) {
		this.callback(NameRef.forObject(hoverHighlight))
	}
}
MODES.SELECT_FEATURE = {
	init: function (callback) {
		this.callback = callback
		const __this = this
		const events = {
			click: function (e) {
				const fl = e.target.featureLink || e.target.getParent().featureLink
				if (fl) {
					e.stopPropagation()
					modePop()
					$('featureDisplay')
						.removeEvents(events)
						.removeClass('selectable')
					__this.mousedown(e)
				}
			},
		}
		$('featureDisplay')
			.addEvents(events)
			.addClass('selectable')
	},
	end: function () {},
	mousemove: function (e, mouseLine) {
		hoverHighlight = getHovering(mouseLine, modelBREP && modelBREP.faces, planes, brepPoints, brepEdges, 16,
			'features') // todo only sketch segments
		paintScreen()
	},
	mousedown: function (e) {
		if (BUTTONS.LEFT == e.button) {
			if (e.shiftKey) {
				selected.toggle(hoverHighlight)
			} else {
				if (hoverHighlight) {
					selected = [hoverHighlight]
				} else {
					selected = []
				}
			}
			this.callback(selected.map(s => NameRef.forObject(s)))
			updateSelected()
		}
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
		editingSketch.elements.remove(this.currentAddingSegment)
		paintScreen()
	},
	mousemove: function (e, mouseLine) {
		let intersection = mouseLine.intersectionWithPlane(editingSketch.plane)
		if (intersection == null) {
			return
		}
		if (!this.currentAddingSegment) {
			this.currentAddingSegment = new this.constructor(editingSketch)
			editingSketch.elements.push(this.currentAddingSegment)
			this.arcmode = 0
		}
		let mouseSC = editingSketch.worldToSketchMatrix.transformPoint(intersection)
		removeFromCoincidence(this.currentAddingSegment.points[this.arcmode], editingSketch)
		//console.log(mousePos);
		const points = getAllPoints(editingSketch)
		points.removeAll(this.currentAddingSegment.points)
		const pointDistances = points
			.map(function (point) { return {point: point, distance: point.distanceToCoords(mouseSC)} })
			.filter(function (pair) { return pair.distance < 16 })
			.sort(function (a, b) { return a.distance - b.distance })
		if (pointDistances.length != 0 && this.currentAddingSegment) {
			makeCoincident(pointDistances[0].point, this.currentAddingSegment.points[this.arcmode], editingSketch)
			mouseSC = pointDistances[0].point
		}
		for (let i = this.arcmode; i < this.currentAddingSegment.points.length; i++) {
			this.currentAddingSegment.points[i].x = mouseSC.x
			this.currentAddingSegment.points[i].y = mouseSC.y
		}
		/*
		 for (var i in lines) {
		 var line = lines[i];
		 //console.log('line', line);
		 if (line != this.currentAddingSegment && line.distanceTo(mousePos.x, mousePos.y) < 16) {
		 mousePos = line.getClosestPoint(mousePos.x, mousePos.y);
		 break;
		 }
		 };
			 */
	},
	mousedown: function (e, mouseLine) {
		const intersection = mouseLine.intersectionWithPlane(editingSketch.plane)
		const sketchCoords = editingSketch.worldToSketchMatrix.transformPoint(intersection)
		if (intersection == null) {
			return
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