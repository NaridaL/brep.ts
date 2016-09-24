"use strict";
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

// I'd use a loop but it kills the type checker.
const ceil = Math.ceil
const floor = Math.floor
const abs = Math.abs
const sign = Math.sign
const atan2 = Math.atan2
const atan = Math.atan
const cos = Math.cos
const sin = Math.sin
const min = Math.min
const max = Math.max
const PI = Math.PI
const sqrt = Math.sqrt
const pow = Math.pow
const round = Math.round
const log = Math.log


function distanceDiffForPointSquared(c) {
	return p => {
		let dist = c.distanceTo(p)
		return V3.create(
			(p.x - c.x) / dist,
			(p.y - c.y) / dist,
			(p.z - c.z) / dist
		)
	}
}

function Point(x, y) {
	this.x = x;
	this.y = y;
}
function Sketch(planeName) {
	// elements in 2D coordinates on x-y plane
	this.planeName = planeName
	this.elements = []
	this.constraints = []
	this.name = "sketch"+(globalId++)
	this.plane = null
}
Sketch.prototype = {
	getConstraintsFor: function (el) {
		return this.constraints.filter((constraint) => constraint.constrains(el))
	},
	constrainDistancePointFixedLineWC: function (point, lineWC, distance) {
		this.b.push(distance);
		var px = this.varMap.get(point), py = px + 1;
		var lineA_SC = this.worldToSketchMatrix.transformPoint(lineWC.anchor);
		var lineB_SC = lineA_SC.plus(this.worldToSketchMatrix.transformVector(lineWC.dir1));
		// console.log(lineA_SC, lineB_SC);
		if (NLA.isZero(distance)) {
			this.F.push(x => distanceLinePointSigned(lineA_SC.x, lineA_SC.y, lineB_SC.x, lineB_SC.y, x[px], x[py]))
		} else {
			this.F.push(x => distanceLinePoint(lineA_SC.x, lineA_SC.y, lineB_SC.x, lineB_SC.y, x[px], x[py]))
		}

		// console.log("calling F", lineA_SC.x, lineA_SC.y, lineB_SC.x, lineB_SC.y, x[px], x[py]);

	},
	constrainDistancePointFixedPlane: function (point, plane, distance) {
		var sketchLineWC = this.plane.intersectionWithPlane(plane);
		if (null == sketchLineWC) throw new Error("no intersection!!");
		this.constrainDistancePointFixedLineWC(point, sketchLineWC, distance);
	},
	constrainAngleSegmentToPlane: function (segment, plane, cosAngle) {
		var sketchLineWC = this.plane.intersectionWithPlane(plane);
		if (null == sketchLineWC) throw new Error("no intersection!!");
		this.constrainAngleSegmentFixedLineWC(segment, sketchLineWC, cosAngle);
	},
	constrainAngleSegmentFixedLineWC: function (segment, lineWC, cosAngle) {
		this.b.push(cosAngle);
		// indexes:
		var ia = this.varMap.get(segment.a), ib = this.varMap.get(segment.b);
		var lineDirectionSC = this.worldToSketchMatrix.transformVector(lineWC.dir1);
		this.F.push(function (x) { return angleVectors(x[ib] - x[ia], x[ib + 1] - x[ia + 1], lineDirectionSC.x, lineDirectionSC.y); })
	},
	constrainEqualDistance: function (ia, ib, ic, id) {
		this.b.push(0);
		this.F.push((x) => distance(x[ia], x[ia+1], x[ib], x[ib+1]) - distance(x[ic], x[ic+1], x[id], x[id+1]));
	},
	constrainAngleSegmentSegment: function(line1, line2, cosAngle) {
		this.b.push(cosAngle * cosAngle);
		var ia = this.varMap.get(line1.a),
			ib = this.varMap.get(line1.b),
			ic = this.varMap.get(line2.a),
			id = this.varMap.get(line2.b)
		this.F.push((x) => {
			let angle = angleABCD(
				x[ia], x[ia + 1],
				x[ib], x[ib + 1],
				x[ic], x[ic + 1],
				x[id], x[id + 1]
			)
			return angle * angle
		})
	},
	constrainDistancePointPoint: function(pA, pB, pDistance) {
		this.b.push(pDistance);
		var ia = this.varMap.get(pA), ib = this.varMap.get(pB);
		this.F.push((x) => distance(x[ia], x[ia + 1], x[ib], x[ib + 1]));
	},
	constrainDistancePointSegment: function(point, segment, distance) {
		this.b.push(distance);
		var ip = this.varMap.get(point), ia = this.varMap.get(segment.a), ib = this.varMap.get(segment.b);
		if (NLA.isZero(distance)) {
			this.F.push(x => distanceLinePointSigned(x[ia], x[ia + 1], x[ib], x[ib + 1], x[ip], x[ip + 1]))
		} else {
			this.F.push(x => distanceLinePoint(x[ia], x[ia + 1], x[ib], x[ib + 1], x[ip], x[ip + 1]))
		}
	},
	constrainDistancePointBezier: function(point, bezier, pDistance) {
		let startT = bezier.getBezierCurve().closestTToPoint(point.V3())
		let startTIndex = this.x.length
		this.x.push(startT)
		this.b.push(pDistance)
		this.b.push(0)
		var ip = this.varMap.get(point),
			ib = bezier.points.map(p => this.varMap.get(p))
		this.F.push(
		x => {
			// cos
			let {x:tPx, y:tPy} = bezierCurveAt(
				x[ib[0]], x[ib[0] + 1],
				x[ib[2]], x[ib[2] + 1],
				x[ib[3]], x[ib[3] + 1],
				x[ib[1]], x[ib[1] + 1],
				x[startTIndex]
			)
			let {x:tPdx, y:tPdy} = bezierCurveTangentAt(
				x[ib[0]], x[ib[0] + 1],
				x[ib[2]], x[ib[2] + 1],
				x[ib[3]], x[ib[3] + 1],
				x[ib[1]], x[ib[1] + 1],
				x[startTIndex]
			)
			return [
				distance(tPx, tPy, x[ip], x[ip + 1]),
				tPdx * (tPx - x[ip]) + tPdy * (tPy - x[ip + 1])]
		})
	},
	constrainAngleSegmentSegment2: function(ab, cd, f, value) {
		// value is in [-2 pi ; 2 pi]; divide by 2 to map to [-pi;pi]
		this.b.push(sin(value)/*, sin(value)*/)
		var ia = this.varMap.get(ab.a),
			ib = this.varMap.get(ab.b),
			ic = this.varMap.get(cd.a),
			id = this.varMap.get(cd.b);
		// calculate the angle for each segment relative to the x axis using atan
		// subtract angle of ab from angle of cd to get signed difference
		// two functions, one calculate sin, one cos of signed difference / 2 to map the signed difference to a point on the unit circle
		// using two functions is necessary so that the resulting functions are "continuous"
		this.F.push(
			(x) => {
				var angle = atan2(
						(x[id + 1] - x[ic + 1]) * f[1],
						(x[id]     - x[ic]    ) * f[1])
					- atan2(
						(x[ib + 1] - x[ia + 1]) * f[0],
						(x[ib]     - x[ia]    ) * f[0]);
				console.log("vaalue", rad2deg(value).toFixed(6), "angle", rad2deg(angle).toFixed(6), cos(angle), sin(angle));
				return sin(
					( atan2(
						(x[id + 1] - x[ic + 1]) * f[1],
						(x[id]     - x[ic]    ) * f[1])
					- atan2(
						(x[ib + 1] - x[ia + 1]) * f[0],
						(x[ib]     - x[ia]    ) * f[0]))) }/*,
			 (x) => sin(
			 ( atan2(
			 (x[id + 1] - x[ic + 1]) * f[1],
			 (x[id]     - x[ic]    ) * f[1])
			 - atan2(
			 (x[ib + 1] - x[ia + 1]) * f[0],
			 (x[ib]     - x[ia]    ) * f[0])))*/ )
	},
	gaussNewtonStep: function () {
		let DISABLE_CONSOLE = true
		DISABLE_CONSOLE && disableConsole();
		var {F, x, b} = this
		var Fx = NLA.Vector.fromFunction(F.length, i => F[i](x))
		console.log('F', F)
		console.log("x", x)
		console.log("Fx", Fx.toString())
		console.log("b", b)
		var jacobi = NLA.Matrix.jacobi(x => F.map(f => f(x)).concatenated(), x, Fx.v)
		console.log("jacobi\n", jacobi.toString())
		var jacobiTranspose = jacobi.transposed()
		var matrix = jacobiTranspose.times((jacobi.times(jacobiTranspose)).inversed())
		var bVector = new NLA.Vector(new Float64Array(b))
		var xDiff = matrix.timesVector(Fx.minus(bVector))
		console.log("matrix\n",matrix.toString(),"\nFx.minus(bVector)", Fx.minus(bVector).toString(), "\nxDiff", xDiff.toString())
		this.x = new NLA.Vector(new Float64Array(x)).minus(xDiff).v;
		Fx = NLA.Vector.fromFunction(F.length, i => F[i](this.x))
		DISABLE_CONSOLE && enableConsole();
		return Fx.minus(bVector)
	},

	/**
	 *
	 * @param segment
	 * @returns {Array.<Edge>}
	 */
	// TODO: check for self-intersections
	getLoopForSegment: function (segment) {
		var startPoint = segment.b;
		var currentPoint = startPoint;
		var loop = [];
		do {
			let currentSegment = currentPoint.line
			let edge = currentSegment.toBrepEdge()
			if (currentSegment.b != currentPoint) {
				edge = edge.flipped()
			}
			loop.push(edge);
			//console.log(currentPoint.coincidence.cs.filter(function (point) { point instanceof SegmentEndPoint && point != currentPoint; }));
			var otherPointsInCoincidence = currentPoint.coincidence && currentPoint.coincidence.cs.filter(function (point) {
					//console.log("point", point, point instanceof SegmentEndPoint, point != currentPoint);
					return point instanceof SegmentEndPoint && point != currentPoint; });
			if (!otherPointsInCoincidence || otherPointsInCoincidence.length != 1) {
				throw new Error("The selected segment is not part of an unambiguous loop.");
			}
			var nextSegmentCoincidencePoint = otherPointsInCoincidence[0]
			var nextSegment = nextSegmentCoincidencePoint.line
			currentPoint = nextSegment.getOtherPoint(nextSegmentCoincidencePoint);
		} while (startPoint.coincidence != currentPoint.coincidence);
		return loop
	},
	delete: function () {
		$('sketchEditor').set('display', 'none')
		featureStack.remove(this)

		rebuildModel()
		updateFeatureDisplay()
	},
	toSource: function () {
		return `(function () {
			let sketch = new Sketch('${this.planeName}')
			let els = sketch.elements = [${this.elements.map(el => el.toSource()).join(',')}]
			sketch.constraints = [${this.constraints.map(el => el.serialize(this.elements)).join(',')}]
			return sketch
		})()`
	}
}
Sketch.prototype.constructor = Sketch

function SegmentEndPoint(x, y, line) {
	this.x = x;
	this.y = y;
	this.line = line;
	this.id = globalId++;
	this.name = "segEndPoint" + (this.id)
	this.coincidence = null
}
SegmentEndPoint.prototype = {
	distanceToCoords: function (p) {
		return distance(this.x, this.y, p.x, p.y)
	},
	toString: function () {
		return "SegmentEndPoint #" + this.id;
	},
	/**
	 *
	 * @returns {V3}
	 */
	V3: function () {
		return V3(this.x, this.y, 0);
	},

	isConstrained: function (sketch) {
		return sketch.constraints.some(constraint => constraint.constrains(this) || constraint.constrains(this.line))
	},
	freeFromConstraints: function (sketch) {
		sketch.constraints.forEach(constraint => constraint.constrains(this) && removeFromConstraint(this, sketch, constraint))
	},
	/**
	 * Can be called even if already removed
	 * @param sketch
	 */
	removeFromSketch: function (sketch) {
		this.line.removeFromSketch(sketch)
	},
	moveCoincidence: function (v3) {
		if (this.coincidence) {
			this.coincidence.cs.forEach(p => {
				p.x = v3.x
				p.y = v3.y
			})
		} else {
			this.x = v3.x
			this.y = v3.y
		}
	},
	get sketch() {
		return this.line.sketch
	},
	canon: function () {
		return this.coincidence && this.coincidence.cs[0] || this
	}
}
SegmentEndPoint.fromV3 = function(segment, p) {
	return new SegmentEndPoint(p.x, p.y, segment)
}
SegmentEndPoint.prototype.constructor = SegmentEndPoint


function makeCoincident(p1, p2, sketch) {
	if (p1.coincidence && p2.coincidence) {
		for (var i = 0; i < p2.coincidence.length;  i++) {
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
		sketch.constraints.push(p1.coincidence = p2.coincidence = new Constraint("coincident", [p1,p2]))
	}
	p2.x = p1.x
	p2.y = p1.y
}
function removeFromCoincidence(p, sketch) {
	if (!p.coincidence) return;
	p.coincidence.cs.remove(p);
	if (p.coincidence.cs.length == 1) {
		sketch.constraints.remove(p.coincidence);
		p.coincidence.cs[0].coincidence = null
	}
	p.coincidence = null;
}
function deleteCoincidence(coincidence, sketch) {
	sketch.constraints.remove(coincidence)
	coincidence.cs.forEach(p => p.coincidence = null)
}
/**
 * Arc goes CCW from a to b
 * @param sketch
 * @param ax
 * @param ay
 * @param bx
 * @param by
 * @param cx
 * @param cy
 * @constructor
 */
function SketchArc(sketch, ax, ay, bx, by, cx, cy) {
	assertInst(Sketch, sketch)
	this.sketch = sketch
	this.a = new SegmentEndPoint(ax, ay, this)
	this.b = new SegmentEndPoint(bx, by, this)
	this.c = new SegmentEndPoint(cx, cy, this)
	this.points = [this.c, this.a, this.b]
	this.id = globalId++
	this.name = "SketchArc" + this.id
}
SketchArc.prototype = {
	angleA: function () {
		return this.a.V3().minus(this.c.V3()).angleXY()
	},
	angleB: function () {
		return this.b.V3().minus(this.c.V3()).angleXY()
	},
	radiusA: function () {
		return this.a.distanceToCoords(this.c)
	},
	distanceToCoords: function (coords) {
		coords = V3(coords)
		var angleA = this.angleA(), angleB = this.angleB()
		if (angleB <= angleA) { angleB += Math.PI * 2 }
		var relCoords = coords.minus(this.c.V3()), angle = relCoords.angleXY(), radius = this.radiusA()
		if (angle < 0) { angle += Math.PI * 2 }
		var angle2 = NLA.clamp(angle, angleA, angleB)
		return V3(radius * cos(angle2), radius * sin(angle2), 0).minus(relCoords).length()
	},
	flip: function () {
		[this.a, this.b] = [this.b, this.a]
	},
	getIntermediatePoints: function () {
		var result = []
		var angleA = this.angleA(), angleB = this.angleB()
		if (angleB <= angleA) { angleB += Math.PI * 2 }
		var radius = this.radiusA()
		var center = this.c.V3()
		var segmentLength = radius * (angleB - angleA), pointCount = floor(segmentLength / 10)
		var intervalAngle = (angleB - angleA) / (pointCount + 1)
		for (var i = 1; i < pointCount + 1; i++) {
			var angle = angleA + i * intervalAngle
			result.push(V3(radius * cos(angle), radius * sin(angle), 0).plus(center))
		}
		return result
	},
	getVectorCA: function () {
		return this.a.V3().minus(this.c.V3())
	},
	getVectorCB: function () {
		return this.b.V3().minus(this.c.V3())
	},
	getOtherPoint: function (p) {
		if (p == this.a) return this.b
		if (p == this.b) return this.a
		assert(false)
	},
	toBrepEdge: function () {
		let ca = this.getVectorCA()
		let curve = new EllipseCurve(this.c.V3(), ca.negated(), ca.negated().getPerpendicular())
		return new PCurveEdge(curve,
			this.a.V3(), this.b.V3(),
			-PI, curve.pointLambda(this.b.V3()),
			null,
			curve.tangentAt(-PI), curve.tangentAt(curve.pointLambda(this.b.V3())),
			this.name)
	}
}
SketchArc.prototype.constructor = SketchArc


function SketchLineSeg(sketch, x1, y1, x2, y2) {
	assertInst(Sketch, sketch)
	this.sketch = sketch
	this.points = [new SegmentEndPoint(x1, y1, this), new SegmentEndPoint(x2, y2, this)]
	this.a = this.points[0]
	this.b = this.points[1]
	this.id = globalId++
	this.name = "segment" + this.id
}
SketchLineSeg.prototype = {
	remove: function () {
		assert (editingSketch.elements.contains(this))
		this.removeFromSketch(editingSketch)
	},
	distanceToCoords: function (coords) {
		return this.distanceTo(coords.x, coords.y);
	},
	angleTo: function (segment) {
		assertInst(SketchLineSeg, segment)
		return segment.angleAB() - this.angleAB();
	},
	toString: function ()  {
		return "SketchLineSeg #" + this.id;
	},
	angleAB: function () {
		return atan2(this.b.y - this.a.y, this.b.x - this.a.x);
	},
	/**
	 * Can be called even if already removed
	 * @param sketch
	 */
	removeFromSketch: function (sketch) {
		this.a.freeFromConstraints(sketch)
		this.b.freeFromConstraints(sketch)
		sketch.elements.remove(this)
		sketch.constraints.forEach(constraint => constraint.constrains(this) && removeFromConstraint(this, sketch, constraint))
		selected.remove(this)
	},
	getOtherPoint: function (p) {
		if (p == this.a) return this.b
		if (p == this.b) return this.a
		assert(false)
	},
	pointLambda: function (v) {
		if (this.b.x - this.a.x > this.b.y - this.a.y) {
			return (v.x - this.a.x) / (this.b.x - this.a.x);
		} else {
			return (v.y - this.a.y) / (this.b.y - this.a.y);
		}
	},
	distanceTo: function (x, y) {
		var x1 = this.a.x, y1 = this.a.y, x2 = this.b.x,  y2 = this.b.y;
		var a = y1 - y2;
		var b = x2 - x1;
		var c = x2 * y1 - x1 * y2;
		var dist = Math.abs(a * x + b * y - c) / length(a, b);
		var xClosest = (b * (b * x - a * y) + a * c) / lengthSquared(a, b);
		var yClosest = (a * ( -b * x + a * y) + b * c) / lengthSquared(a, b);
		if (x1 != x2 ? isBetween(xClosest, x1, x2) : isBetween(yClosest, y1, y2)) {
			return dist;
		} else {
			if (x1 < x2 && x < x1 || x1 > x2 && x > x1
				|| x1 == x2 && (y1 < y2 && y < y2 || y1 > y2 && y > y1)) {
				//noinspection JSSuspiciousNameCombination
				return distance(x, x1, y, y1);
			} else {
				//noinspection JSSuspiciousNameCombination
				return distance(x, x2, y, y2);
			}
		}
	},
	getVectorAB: function () {
		return V3(this.b.x - this.a.x, this.b.y - this.a.y, 0);
	},
	getClosestPoint: function (x, y) {
		var x1 = this.a.x, y1 = this.a.y, x2 = this.b.x, y2 = this.b.y;
		var a = y1 - y2;
		var b = x2 - x1;
		var c = x2 * y1 - x1 * y2;
		var dist = abs(a * x + b * y - c) / length(a, b);
		var xClosest = (b * (b * x - a * y) + a * c) / lengthSquared(a, b);
		var yClosest = (a * ( -b * x + a * y) + b * c) / lengthSquared(a, b);
		if (isBetween(xClosest, x1, x2)) {
			return {"x": xClosest, "y": yClosest};
		} else {
			if (x1 < x2 && x < x1 || x1 > x2 && x > x1
				|| x1 == x2 && (y1 < y2 && y < y2 || y1 > y2 && y > y1)) {
				return {x: x1, y: y1}
			} else {
				//noinspection JSSuspiciousNameCombination
				return {x: y1, y: y2}
			}
		}
	},
	length: function () {
		return this.points[0].distanceTo(this.points[1]);
	},
	intersection: function (segment) {
		return intersection(this.a.x, this.a.y, this.b.x, this.b.y,
			segment.a.x, segment.a.y, segment.b.x, segment.b.y);
	},
	getL3: function() {
		return L3.anchorDirection(this.a.V3(), this.getVectorAB())
	},
	toBrepEdge: function() {
		return StraightEdge.throughPoints(this.a.V3(), this.b.V3(), this.name)
	}
}
SketchLineSeg.prototype.constructor = SketchLineSeg

function intersection(x1, y1, x2, y2, x3, y3, x4, y4) {
	var denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
	var xNominator = (x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4);
	var yNominator = (x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4);
	return V3(xNominator / denominator, yNominator / denominator, 0);
}

function isBetween(val, a, b) {
	if (a < b) {
		return a < val && val < b;
	} else {
		return b < val && val < a;
	}
}

function distance(x1, y1, x2, y2) {
	return length(x2 - x1, y2 - y1);
}
function distanceSquared(x1, y1, x2, y2) {
	return lengthSquared(x2 - x1, y2 - y1);
}
function bezierCurveAt(x0, y0, x1, y1, x2, y2, x3, y3, t) {
	//console.log(x0, y0, x1, y1, x2, y2, x3, y3, t)
	let s = 1 - t, c0 = s * s * s, c1 = 3 * s * s * t, c2 = 3 * s * t * t, c3 = t * t * t
	return {
		x: x0 * c0 + x1 * c1 + x2 * c2 + x3 * c3,
		y: y0 * c0 + y1 * c1 + y2 * c2 + y3 * c3
	}
}
function bezierCurveTangentAt(x0, y0, x1, y1, x2, y2, x3, y3, t) {
	//console.log(x0, y0, x1, y1, x2, y2, x3, y3, t)
	let s = 1 - t, c01 = 3 * s * s, c12 = 6 * s * t, c23 = 3 * t * t
	return {
		x: (x1 - x0) * c01 + (x2 - x1) * c12 + (x3 - x2) * c23,
		y: (y1 - y0) * c01 + (y2 - y1) * c12 + (y3 - y2) * c23,
	}
}
function angle(x1, y1, x2, y2) {
	return length(x2 - x1, y2 - y1);
}

function length(x, y) {
	return Math.sqrt(x * x + y * y);
}
function lengthSquared(x, y) {
	return x * x + y * y
}

function distanceLinePoint(x1, y1, x2, y2, x, y) {
	let a = y1 - y2;
	let b = x2 - x1;
	let c = x2 * y1 - x1 * y2;
	let dist = Math.abs(a * x + b * y - c) / length(a, b);
	return dist;
}


function distanceLinePointSigned(x1, y1, x2, y2, x, y) {
	// function needs to be differentiable around target value
	// therefore this funciton is necessary
	let a = y1 - y2;
	let b = x2 - x1;
	let c = x2 * y1 - x1 * y2;
	//TODOlet dist = Math.abs(a * x + b * y - c) / length(a, b);
	let dist = (a * x + b * y - c) / length(a, b);
	return dist;
}

function angleVectors(ax, ay, bx, by) {
//	console.log(ax, ay, bx, by);
	return Math.abs(ax * bx + ay * by) / Math.sqrt((ax * ax + ay * ay) * (bx * bx + by * by));
}

// returns angle between segments AB, CD
function angleABCD(ax, ay, bx, by, cx, cy, dx, dy) {
	return angleVectors(bx - ax, by - ay, dx - cx, dy - cy);
	// return ((bx - ax) * (dx - cx) + (by - ay) * (dy - cy)) / Math.sqrt(((bx - ax)*(bx - ax)+(by - ay)*(by - ay))*((dx - cx)*(dx - cx)+(dy - cy)*(dy - cy)))
}
// 4 vars per line: s0, s1, d0, d1
// (d0, d1)^T is normalized
// x.length = 4 * lines.length
function recalculate(sketch) {
	function pushPoint(p) {
		var varMap = sketch.varMap, x = sketch.x
		if (varMap.has(p) /*	|| !p.isConstrained()*/) {
			return;
		}
		if (p.coincidence) {
			p.coincidence.cs.forEach(function (p2) {
				varMap.set(p2, xIndex);
			});
		} else {
			varMap.set(p, xIndex);
		}
		x.push(p.x, p.y);
		xIndex += 2;
	}
	console.log("recalculating");
	Object.defineProperty(sketch, 'x', {enumerable: false, value: [], writable: true})
	Object.defineProperty(sketch, 'b', {enumerable: false, value: [], writable: true})
	Object.defineProperty(sketch, 'F', {enumerable: false, value: [], writable: true})
	Object.defineProperty(sketch, 'varMap', {enumerable: false, value: new Map(), writable: true})
	// init x to current values
	var xIndex = 0;
	sketch.elements.forEach(function (seg) {
		seg.points.forEach(pushPoint)
		if (seg instanceof SketchArc && sketch.varMap.get(seg.a) != sketch.varMap.get(seg.b)) {
			sketch.constrainEqualDistance.apply(sketch, [seg.a, seg.c, seg.b, seg.c].map(point => sketch.varMap.get(point)))
		}
	});
	console.log("varMap", sketch.varMap)
	var constraintFunctionSources = []
	sketch.constraints.forEach(function (cst) {
		//console.log(cst);
		if (cst.type == "parallel" || cst.type == "colinear" || cst.type == "equalLength") {
			if (cst.fixed) {
				for (let j = 0; j < cst.segments.length; j++) {
					if ("colinear" == cst.type) {
						sketch.constrainDistancePointFixedPlane(cst.segments[j].a, cst.fixed.plane, 0);
						sketch.constrainDistancePointFixedPlane(cst.segments[j].b, cst.fixed.plane, 0);
					} else {
						sketch.constrainAngleSegmentToPlane(cst.segments[j], cst.fixed.plane, 1);
					}
				}
			} else {
				for (let j = 1; j < cst.segments.length; j++) {
					var first = cst.segments[0], second = cst.segments[j];
					if ("parallel" == cst.type) {
						// assume that max. one element can be a line or a plane
						sketch.constrainAngleSegmentSegment(first, second, 1);
					}
					if ("colinear" == cst.type) {
						sketch.constrainDistancePointSegment(second.a, first, 0)
						sketch.constrainDistancePointSegment(second.b, first, 0)
					}
					if ("equalLength" == cst.type) {
						sketch.constrainEqualDistance.apply(sketch,
							[0, j]
								.map((segmentsIndex) => cst.segments[segmentsIndex].points)
								.concatenated()
								.map((point) => sketch.varMap.get(point)))
					}
				}
			}
		}
		if (cst.type == "pointDistance") {
			sketch.constrainDistancePointPoint(cst.cs[0], cst.cs[1], cst.distance);
		}
		if (cst.type == "pointOnLine" || cst.type == "pointLineDistance" || cst.type == "pointPlaneDistance") {
			var distance = cst.type != "pointOnLine" ? cst.distance : 0;
			if (cst.other.plane) {
				sketch.constrainDistancePointFixedPlane(cst.point, cst.other.plane, distance)
			} else if(cst.other instanceof SketchLineSeg) {
				sketch.constrainDistancePointSegment(cst.point, cst.other, distance)
			} else if (cst.other instanceof SketchBezier) {
				sketch.constrainDistancePointBezier(cst.point, cst.other, distance)
			} else if (cst.other instanceof SketchArc) {
				sketch.constrainEqualDistance.apply(sketch, [cst.other.a, cst.other.c, cst.point, cst.other.c].map(point => sketch.varMap.get(point)))
			}
		}
		if (cst.type == "angle") {
			sketch.constrainAngleSegmentSegment2(cst.cs[0], cst.cs[1], cst.f, cst.value);
		}
		if (cst.type == "perpendicular") {
			var cosValue = cst.type == "angle" ? cos(cst.value) : 0
			// assume that max. one element can be a line or a plane
			if (cst.other.plane) {
				sketch.constrainAngleSegmentToPlane(cst.segment, cst.other.plane, cosValue);
			} else {
				sketch.constrainAngleSegmentSegment(cst.segment, cst.other, cosValue);
			}
			/*
			 constrainAngle2D.apply(undefined,
			 [0, 1]
			 .map((whichIndex) => cst.constrains[whichIndex].points).concatenated()
			 .map((point) => varMap.get(point))
			 .map((pointXCoord) => [pointXCoord, pointXCoord + 1]).concatenated()
			 .concat(cst.value));
			 */
		}
//			console.log("added constraint, b:", b);
		for (var i = constraintFunctionSources.length; i < sketch.b.length; i++) {
			constraintFunctionSources.push(cst)
		}
	});
	if (sketch.b.isEmpty()) {
		return;
	}
	for (var count = 0; count < 100; count++) {
		var lastDiffs = sketch.gaussNewtonStep(), lastSize = lastDiffs.length()
		if (lastSize < NLA_PRECISION / 1000) {
			break;
		}
	}
	//enableConsole()
	console.log(`broke at ${count}, lastDiffs ${lastDiffs}, lastSize ${lastSize}`)
	// first set all of them to false
	sketch.constraints.forEach(cst => cst.error = false)
	if (lastSize < 1) {
		console.log("REVERSING")
		reverse(sketch);
	} else {
		// then mark every constraint which has a problematic function
		lastDiffs.v.forEach((el, index) => constraintFunctionSources[index].error |= !NLA.isZero(el))
	}
}


function paintLineXY(a, b, color, width) {
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
function randomVec4Color(opacity) {
	opacity = opacity || 1.0
	return [Math.random(), Math.random(), Math.random(), opacity]
}
function paintArc(center, radius, width, color, startAngle, endAngle) {
	startAngle = isFinite(startAngle) ? startAngle : 0
	endAngle = isFinite(endAngle) ? endAngle : 2 * Math.PI
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
	}).draw(meshes.segment);
}
var savedTextures = new Map();
function getTextureForString(str) {
	var texture
	if (texture = savedTextures.get(str)) {
		return texture
	}

	var canvas = $("textTextureCanvas");
	var ctx = canvas.getContext('2d');

	var font = TEXT_TEXTURE_HEIGHT + "px Anonymous Pro";

	ctx.font = font ;
	var textWidthPx = ctx.measureText(str).width;
	var canvasWidth = 1 << ceil(log(textWidthPx)/log(2));
	canvas.width = canvasWidth;
	canvas.height = TEXT_TEXTURE_HEIGHT;


	ctx.font = font;
	ctx.fillStyle = "#ffffff";
	ctx.textAlign = "left";	// This determines the alignment of text, e.g. left, center, right
	ctx.textBaseline = "top";	// This determines the baseline of the text, e.g. top, middle, bottom
	ctx.imageSmoothingEnabled = false;


	ctx.fillText(str, 0, 0)


	texture = GL.Texture.fromImage(canvas, {minFilter: gl.LINEAR_MIPMAP_LINEAR, magFilter: gl.LINEAR})
	texture.textWidthPx = textWidthPx

	savedTextures.set(str, texture)

	return texture
}
var zoomFactor = 1;
function paintSegments(sketch) {

	if (!sketch.plane) return
	//console.log("painting segments", sketch.elements.length);
	/*ctx.clearRect (0, 0, ctx.canvas.width, ctx.canvas.height);
	 ctx.fillStyle="rgb(100, 100, 255)";
	 ctx.lineWidth=2;*/
	//console.log(sketch.sketchToWorldMatrix);
	gl.multMatrix(sketch.sketchToWorldMatrix)
	sketch.elements.forEach(function (seg) {
		function drawPoint(p) {
			paintArc(p.V3(), 1.5, 3, colorFor(highlighted.contains(p) || hoverHighlight == p, selected.contains(p)))
		}
		let color = colorFor(highlighted.contains(seg) || hoverHighlight == seg, selected.contains(seg))
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

	});

	paintConstraints(sketch);
}
function colorFor(highlighted, selected) {
	return !selected
		? (!highlighted ? 0x33CCFF : 0x145266)
		: (!highlighted ? 0xFF3399 : 0x330A1E)
}
function getAllPoints(sketch) {
	return sketch.elements.map(segment => segment.points).concatenated()
}
function paintConstraints(sketch) {
	var crossCount = 2, parallelCrossCount = 1;
	sketch.constraints.forEach(function (cst) {
		paintDefaultColor = hoverHighlight == cst ? (!cst.error ? 0x00ff00 : 0xffff00) : (!cst.error ? 0x000000 : 0xff0000)
		switch (cst.type) {
			case "coincident":
				var point = cst.cs[0]
				paintArc(point.V3(), 4, 1, colorFor(false, false), 0, 2 * Math.PI)
				paintArc(point.V3(), 1.5, 3, colorFor(highlighted.contains(point) || hoverHighlight == point, selected.contains(point)), 0, 2 * Math.PI)
				break;
			case "parallel": {
				var dir1 = cst.segments[0].getVectorAB().normalized()
				var dir90 = dir1.getPerpendicular().normalized()
				var ARR_SPACING = 5
				for (var c = 0; c < cst.segments.length; c++) {
					let line = cst.segments[c]
					let ab = line.getVectorAB()
					let abLength = ab.length()
					let ab1 = ab.normalized()
					for (let i = 0; i < parallelCrossCount; i++) {
						var s = abLength / 2 - ARR_SPACING * parallelCrossCount / 2 + i * ARR_SPACING - 10
						var arrPoint = line.a.V3().plus(ab1.times(s))
						//console.log("crosspos", crossStart, crossEnd);
						paintLineXY(arrPoint.plus(dir90.times(-1)), arrPoint.plus(dir1.times(6).plus(dir90.times(5))))
						paintLineXY(arrPoint.plus(dir90.times(1)), arrPoint.plus(dir1.times(6).plus(dir90.times(-5))))
					}
				}
				parallelCrossCount++;
				break;
			}
			case "perpendicular":
				if (cst.other instanceof SketchLineSeg) {
					let ab = cst.segment.getVectorAB()
					let cd = cst.other.getVectorAB()
					let intersection = cst.segment.intersection(cst.other)
					let abPos = cst.segment.pointLambda(intersection)
					let cdPos = cst.other.pointLambda(intersection)
					let abLine = ab.normalized().times(0.5 < abPos ? -16 : 16)
					let cdLine = cd.normalized().times(0.5 < cdPos ? -16 : 16)
					paintLineXY(intersection.plus(abLine).plus(cdLine), intersection.plus(abLine))
					paintLineXY(intersection.plus(abLine).plus(cdLine), intersection.plus(cdLine))
				} else {
					let ab = cst.segment.getVectorAB()
					let sketchLineWC = sketch.plane.intersectionWithPlane(cst.other.plane)
					let dir = sketch.worldToSketchMatrix.transformVector(sketchLineWC.dir1)
					let p = sketch.worldToSketchMatrix.transformPoint(sketchLineWC.anchor)
					let abXcd = ab.cross(dir)
					let div = abXcd.lengthSquared()
					let anchorDiff = p.minus(cst.segment.a.V3())
					let abPos = anchorDiff.cross(dir).dot(abXcd) / div
					let linePos = anchorDiff.cross(ab).dot(abXcd) / div
					let intersection = p.plus(dir.times(linePos))
					let abLine = ab.normalized().times(0.5 < abPos ? -16 : 16)
					let cdLine = dir.times(16)
					paintLineXY(intersection.plus(abLine).plus(cdLine), intersection.plus(abLine))
					paintLineXY(intersection.plus(abLine).plus(cdLine), intersection.plus(cdLine))
				}
				break
			case "colinear":
				// assume segments dont overlap
				var segments = cst.segments;
				var ab = segments[0].getVectorAB().normalized();
				var coord = ab.x > ab.y ? "x" : "y";
				if (ab[coord] < 0) {
					ab = ab.times(-1);
				}
				var scale = ab[coord];
				var offsetPoint = segments[0].a.V3().minus(ab.times(segments[0].a.V3()[coord] / scale));
				segments.sort((a, b) => a.a[coord] - b.a[coord]);
				for (var i = 0; i < segments.length - (cst.fixed ? 0 : 1); i++) {
					var startS = max(segments[i].a.V3()[coord]/scale, segments[i].b.V3()[coord]/scale) + 8;
					var endS = startS + 20;
					paintLineXY(offsetPoint.plus(ab.times(startS)), offsetPoint.plus(ab.times(endS)), undefined, 2);

					startS = min(segments[(i + 1) % segments.length].a[coord]/scale, segments[(i + 1) % segments.length].b[coord]/scale) - 8;
					endS = startS - 20;
					paintLineXY(offsetPoint.plus(ab.times(startS)), offsetPoint.plus(ab.times(endS)), undefined, 2);
				}
				break;
			case "pointDistance": {
				let a = cst.cs[0].V3(), b = cst.cs[1].V3();
				let ab = b.minus(a);
				let ab1 = ab.normalized();
				let abLength = ab.length();
				let ab90 = ab1.getPerpendicular();
				let texture = getTextureForString(cst.distance);
				let textLength = texture.textWidthPx / TEXT_TEXTURE_HEIGHT * 20;
				paintLineXY(a.plus(ab90.times(6)), a.plus(ab90.times(22)));
				paintLineXY(b.plus(ab90.times(6)), b.plus(ab90.times(22)));

				paintLineXY(a.plus(ab90.times(14)), a.plus(ab90.times(14)).plus(ab1.times(abLength / 2 - textLength / 2 - 5)));
				paintLineXY(a.plus(ab90.times(14)).plus(ab1.times(abLength / 2 + textLength / 2 + 5)), a.plus(ab90.times(14)).plus(ab1.times(abLength)));
				let textCenter = a.plus(ab90.times(14)).plus(ab1.times(abLength / 2));
				gl.pushMatrix();
				gl.translate(textCenter);
				gl.scale(20, 20, 10);
				gl.multMatrix(M4.forSys(ab1, ab90.times(1), V3.Z));
				gl.translate(-textLength / 2 / 20, -0.5, 0);
				renderText(cst.distance, 0xff0000);
				gl.popMatrix();
				break;
			}
			case "pointLineDistance": {
				let texture = getTextureForString(cst.distance)
				let p = cst.point.V3(), ab = cst.other.getVectorAB(), a = cst.other.a.V3()
				let ab1 = ab.normalized()
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
				gl.multMatrix(M4.forSys(apProj.normalized(), ab1, V3.Z))
				gl.translate(-textLength / 2 / 20, -0.5, 0)
				renderText(cst.distance, 0xff0000)
				gl.popMatrix()
				break
			}
			case "pointPlaneDistance": {
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
			case "pointOnLine": {
				if (cst.other instanceof SketchBezier) {
					let bez = cst.other.getBezierCurve()
					let t = bez.closestTToPoint(cst.point.V3())
					let tangentLength = bez.tangentAt(t).length()
					const ps = cst.other.points.map(p => p.V3())
					paintBezier(ps, 3, 0x000000, t + 4 / tangentLength, t + 12 / tangentLength)
					paintBezier(ps, 3, 0x000000, t - 12 / tangentLength, t - 4 / tangentLength)
					break
				}
				if (cst.other instanceof SketchArc) {
					let angle = cst.point.V3().minus(cst.other.c.V3()).angleXY()
					let radius = cst.other.radiusA()
					paintArc(cst.other.c.V3(), radius, 3, 0x000000, angle + 4 / radius, angle + 12 / radius)
					paintArc(cst.other.c.V3(), radius, 3, 0x000000, angle - 12 / radius, angle - 4 / radius)
					break
				}
				if (cst.other.plane) break
				let ab = cst.other instanceof SketchLineSeg ? cst.other.getVectorAB() : sketch.plane.normal.cross(cst.other.plane.normal)
				let ab1 = ab.normalized()
				let p = cst.point.V3()
				paintLineXY(p.plus(ab1.times(12)), p.plus(ab1.times(4)), 0x000000, 3)
				paintLineXY(p.plus(ab1.times(-12)), p.plus(ab1.times(-4)), 0x000000, 3)
				break
			}
			case "angle":
				var first = cst.cs[0], second = cst.cs[1]
				var intersection = first.intersection(second)
				var startAngle = (first.angleAB() + (cst.f[0] == -1 ? PI : 0)) % (2 * PI), endAngle = startAngle + cst.value
				paintArc(intersection, 20, 1, 0x000000, startAngle, endAngle)
				break
			case "equalLength":
				for (let c = 0; c < cst.segments.length; c++) {
					let line = cst.segments[c];
					let ab = line.getVectorAB();
					let abLength = ab.length();
					let ab1 = ab.normalized();
					let ab90 = ab.getPerpendicular().normalized();
					let crossLength = 10;
					let crossSpacing = 3;
					for (let i = 0; i < crossCount; i++) {
						let s = abLength / 2 - crossSpacing * crossCount / 2 + i * crossSpacing + 10
						let crossStart = line.a.V3().plus(ab1.times(s)).plus(ab90.times(crossLength / 2))
						let crossEnd = line.a.V3().plus(ab1.times(s)).plus(ab90.times(-crossLength / 2))
						//console.log("crosspos", crossStart, crossEnd)
						paintLineXY(crossStart, crossEnd)
					}
				}
				crossCount++;
				break
			default:
				throw new Error("unknown cst " + cst.type)
			// console.log("errror", cst.type);
		}

	});
}

var PlaneDefinition = class PlaneDefinition {
	/*
	if ("face" == feature.planeType && feature.faceName) {
		var face = modelBREP.faces.find(face => face.name == feature.faceName)
		var plane = face.surface.plane, right = plane.normal.getPerpendicular().normalized(), up = plane.normal.cross(right)

		var cp
		planes.push(cp = new CustomPlane(plane.anchor.plus(plane.normal.times(feature.offset)),
			right, up, -500, 500, -500, 500, 0xFF4500, feature.planeId))
	}
	if ("immediate" == feature.planeType) {
		var plane = eval(feature.source), right = plane.normal.getPerpendicular().normalized(), up = plane.normal.cross(right)

		planes.push(new CustomPlane(plane.anchor,
			right, up, -500, 500, -500, 500, 0xFF4500, feature.planeId))
	}
	*/
	constructor(planeId) {
		assert(planeId)
		this.type = "planeDefinition"
		this.planeId = planeId
		this.planeType = 'dynamic'
		this.source = ''
		this.offset = 0
		this.angle = 0
		this.whats = []
	}

	toSource(line) {
		return `new PlaneDefinition()`
	}

	delete(line) {
		$('planeDefiner').set('display', 'none')
		featureStack.remove(this)

		rebuildModel()
		updateFeatureDisplay()
	}
}

function serializeFeatures(features) {
	return serialize(features)
	//return '[' + featureStack.map(f => f.toSource()).join(',') + ']'
}

function load(key) {
	$('saveNameInput').value = key
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
	let loadSelect = $('loadSelect')
	let keys = NLA.arrayFromFunction(localStorage.length, i => localStorage.key(i))
	keys.sort()
	keys.forEach(key => loadSelect.adopt(new Element('option', {html: key})))
	loadSelect.onchange = function () {
		let key = loadSelect.value
		if (key) {
			load(key)
		}
	}


	let saveButton = $('saveButton')
	saveButton.onclick = function () {
		let key = $('saveNameInput').value
		if (null == key) {
			loadSelect.adopt(new Element('option', {html: key}))
		}
		localStorage.setItem(key, serializeFeatures(featureStack))
		console.log("saved " + key)
	}
}

function isSketchEl(el) {
	return el instanceof SketchBezier || el instanceof  SegmentEndPoint
			|| el instanceof SketchLineSeg || el instanceof SketchArc
}

//var sketchPlane = new CustomPlane(V3(0, 0,1), V3.X, V3.Y, -500, 500, -500, 500, 0xff00ff);
var /** @type Sketch */ editingSketch, /** @type Array */ featureStack = []
function initModel() {


	return
	featureStack.push(new PlaneDefinition("immediate", "new P3(V3.X, 450.00086883286616)", "sketchPlane"))
	//featureStack.push({type: "planeDefinition", planeType: "immediate", source: "P3.normalOnAnchor(V3(1, 0, -1).normalized(), V3.X)", planeId: "sketchPlane", delete: PlaneDefinition.prototype.delete})
	var sketch = new Sketch("sketchPlane");
	featureStack.push(sketch);
	editingSketch = sketch
	console.log("init sketch");
	sketch.elements.push(new SketchLineSeg(sketch, 0,100,100,0));
	let sb = SketchBezier.forBezierCurve(BezierCurve.graphXY(2,-3,-3,2).scale(100, 100, 100), sketch)
	sketch.elements.push(sb)
	console.log("sb", sb)
	//sketch.elements.push(new SketchArc(0,100,100,0, 0,0))
	sketch.elements.push(new SketchLineSeg(sketch, 0,0,100 * sqrt(3) /2,50));
	//sketch.elements.push(new SketchLineSeg(100 * sqrt(3) /2,50, 100, 100));
	//sketch.elements.push(new SketchLineSeg(100,100,0,0));
	makeCoincident(sketch.elements[0].b, sketch.elements[1].a, sketch);
	//makeCoincident(sketch.elements[1].b, sketch.elements[2].a, sketch);
	//makeCoincident(sketch.elements[2].b, sketch.elements[0].a, sketch);
	selected = [sketch.elements[0].b, sketch.elements[1]];


	selected = sketch.elements.slice(0, 2)
	makePerpendicular()
	sketch.constraints[0].error = true
	//makeGroup("equalLength")
	selected = [sketch.elements[0], sketch.elements[1]]
	//makeAngle()
	//makePerpendicular()
	//makeGroup("parallel")
	//console.log("angle" ,rad2deg(sketch.elements[0].angleTo(sketch.elements[1])));
	//makeGroup("colinear");
	//sketch.constraints[0].value = deg2rad(30);
	//selected = sketch.elements.slice(1, 3);
	//makeDistance();
	//selected = [sketch.elements[0], sketch.elements[1]];
	//makeAngle();
	console.log("constraints", sketch.constraints);
	featureStack.push(new Extrude("initExtrude", sketch.elements[0].name))
	//featureStack.push({type: "planeDefinition", planeType: "face", faceName: "initExtrudewall0", offset: 0, planeId: "planeCustom1", delete: PlaneDefinition.prototype.delete})
	//featureStack.push(editingSketch = new Sketch("planeCustom1"))
	//editingSketch.elements.push(new SketchLineSeg(0,10,100,0));
	//rebuildModel();
	//selected = [editingSketch.elements[0], new NameRef("initExtruderoof").get()]
	//makeGroup("parallel")

	rebuildModel();

	/*
	 sketch.elements.push(new SketchLineSeg(10,10,10,100));
	 sketch.elements.push(new SketchLineSeg(10,10,100,10));
	 sketch.elements.push(new SketchLineSeg(10,100,100,10));
	 makeCoincident(sketch.elements[0].a, sketch.elements[1].a, sketch);
	 makeCoincident(sketch.elements[0].b, sketch.elements[2].a, sketch);
	 makeCoincident(sketch.elements[1].b, sketch.elements[2].b, sketch);
	 */
}
function directionObjectToVector(dirObj) {
	if (dirObj instanceof CustomPlane) {
		return dirObj.normal;
	} else {
		console.log(dirObj)
		throw new Error("uuuh"+ dirObj);
	}
}
function rebuildModel() {
	let DISABLE_CONSOLE = false
	console.log("rebuilding model")
	//NLA.DEBUG = false
	DISABLE_CONSOLE && disableConsole()
	planes = [
		CustomPlane(V3.ZERO, V3.Y, V3.Z, -500, 500, -500, 500, 0xffaaaa, "planeYZ"),
		CustomPlane(V3.ZERO, V3.Z, V3.X, -500, 500, -500, 500, 0xaaffaa, "planeZX"),
		CustomPlane(V3.ZERO, V3.X, V3.Y, -500, 500, -500, 500, 0xaaaaff, "planeXY"),
	];
	modelBREP = undefined
	brepMesh = undefined
	brepPoints = []
	brepEdges = []
	featureStack.forEach((feature) => {
		// try {
			if (feature instanceof Sketch) {
				feature.plane = planes.find(p => p.name == feature.planeName);
				//console.log("LENGTHS", feature.plane.right.length(), feature.plane.up.length(), feature.plane.normal.length())
				if (feature.plane) {
					feature.sketchToWorldMatrix =
						M4.forSys(feature.plane.right, feature.plane.up, feature.plane.normal, feature.plane.anchor)
					feature.worldToSketchMatrix = feature.sketchToWorldMatrix.inversed();
					recalculate(feature)
				}
			} else if (feature.type && feature.type == "extrude") {
				let loopSketch = featureStack.filter(f => f instanceof Sketch)
					.find(sketch => sketch.elements.some(el => el.name == feature.segmentName))
				if (loopSketch.plane) {
					let loopSegment = loopSketch.elements.find(el => el.name == feature.segmentName)
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
					//console.log("polypoints", polygonPoints, polygonPoints.toSource(), loopSketch.plane.translated().toSource(), offsetDir.times(feature.start))
					let brep = B2.extrudeEdges(edgeLoop, loopSketch.plane.translated(startOffset), lengthOffset, feature.name)
					if (modelBREP) {
						isEdges = modelBREP.getIntersectionEdges(brep)
						// drVs = isEdges.map(e => ({anchor: e.a, dir: e.curve.tangentAt(e.aT).normalized()}))
						// modelBREP = modelBREP[feature.operation](brep)
					} else {
						modelBREP = brep;
					}
					brepMesh = modelBREP.toMesh()
					modelBREP.faces.forEach(face =>
						brepEdges.pushAll(face.getAllEdges().filter(edge => !edge.flippedOf || edge.id < edge.flippedOf.id)))
					// brepEdges.pushAll(isEdges)
					brepPoints = Array.from(modelBREP.vertexNames.keys())
					console.log(modelBREP.vertexNames)
					// brepMesh.computeWireframeFromFlatTriangles()
					// brepMesh.compile()
				}
			} else if (feature.type && "planeDefinition" == feature.type) {
					let face = modelBREP.faces.find(face => face.name == feature.faceName)
					let sel = feature.whats.map(w => w.get())
					let plane = MODES.PLANE_DEFINITION.magic(sel, feature.angle * DEG)
					console.log('PLANE', plane, feature.whats.sce, sel.sce, feature.whats.length, sel.length)
					if (plane) {
						let right = plane.normal.getPerpendicular().normalized(), up = plane.normal.cross(right)

						var cp
						planes.push(cp = new CustomPlane(plane.anchor.plus(plane.normal.times(feature.offset)),
							right, up, -500, 500, -500, 500, 0xFF4500, feature.planeId))
					}
				if ("immediate" == feature.planeType) {
					let plane = eval(feature.source), right = plane.normal.getPerpendicular().normalized(), up = plane.normal.cross(right)

					planes.push(new CustomPlane(plane.anchor,
						right, up, -500, 500, -500, 500, 0xFF4500, feature.planeId))
				}
			} else {
				//noinspection ExceptionCaughtLocallyJS
				throw new Error("unknown feature");
			}

		// } catch (error) {
		// 	let featureDiv = $('featureDisplay').getChildren().filter(child => child.featureLink == feature)[0]
		//
		// 	if (featureDiv) {
		// 		let ediv = featureDiv.getElement('[name=error]')
		// 		ediv.setStyle('display', 'inline')
		// 		ediv.title= error.toString() + '\n' + error.stack
		// 	}
		// 	throw error
		// 	console.log(error)
		// 	// TODO stop rebuild
		// }
	})
	if (brepMesh) {
		brepMesh.computeNormalLines(10)
		//brepMesh.computeWireframeFromFlatTriangles()
		brepMesh.compile()
	}
	paintScreen()
	DISABLE_CONSOLE && enableConsole()
}
// DECLARATIONS
// GLOBALS
var faces = []
var highlighted = [], selected = [], paintDefaultColor = 0x000000
var /** @type LightGLContext */ gl;
var modelBREP, brepMesh, brepPoints, planes, brepEdges, isEdges
var drPs = [], drVs = []
var oldConsole = undefined;
var eyePos = V3(1000, 1000, 1000), eyeFocus = V3.ZERO, eyeUp = V3.Z
var hoverHighlight = undefined
// console.log = oldConsole;
var modeStack = []
var actionStack = []
var shaders = {}

function rgbToVec4(color) {
	return [(color >> 16) / 255.0, ((color >> 8) & 0xff) / 255.0, (color & 0xff) / 255.0, 1.0];
}
function renderColor(mesh, color, mode) {
	shaders.singleColor.uniforms({
		color: rgbToVec4(color)
	}).draw(mesh)
}
function renderColorLines(mesh, color) {
	shaders.singleColor.uniforms({
		color: rgbToVec4(color)
	}).draw(mesh, 'LINES');
}
var TEXT_TEXTURE_HEIGHT = 128;
function renderText(string, color) {
	var texture = getTextureForString(string);
	texture.bind(0);
	gl.pushMatrix();
	gl.scale(texture.width / TEXT_TEXTURE_HEIGHT, 1, 1);
	shaders.textureColor.uniforms({
		texture: 0,
		color: rgbToVec4(color)
	}).draw(meshes.text);
	gl.popMatrix();
}
function drawVectors() {
	gl.pushMatrix();
	gl.scale(50, 50, 50)
	shaders.singleColor.uniforms({
		color: rgbToVec4(0xff0000)
	}).draw(meshes.vector);
	gl.popMatrix();

	gl.pushMatrix();
	gl.multMatrix(M4.forSys(V3.Y, V3.Z, V3.X))
	gl.scale(50, 50, 50)
	shaders.singleColor.uniforms({
		color: rgbToVec4(0x00ff00)
	}).draw(meshes.vector);
	gl.popMatrix();

	gl.pushMatrix();
	gl.multMatrix(M4.forSys(V3.Z, V3.X, V3.Y))
	gl.scale(50, 50, 50)
	shaders.singleColor.uniforms({
		color: rgbToVec4(0x0000ff)
	}).draw(meshes.vector);
	gl.popMatrix();

	drVs.forEach(vi => {
		gl.pushMatrix()
		let v2 = vi.dir.getPerpendicular()
		gl.multMatrix(M4.forSys(vi.dir, v2, vi.dir.cross(v2), vi.anchor))
		gl.scale(50, 50, 50)
		shaders.singleColor.uniforms({
			color: rgbToVec4(vi.color || 0x0000ff)
		}).draw(meshes.vector);
		gl.popMatrix()
	})
}

function drawPoints() {
	drPs.forEach(info => {
		let p = info.p || info
		gl.pushMatrix()
		gl.translate(p)
		gl.scale(5,5,5)
		shaders.singleColor.uniforms({ color: rgbToVec4(0xcc0000)}).draw(meshes.sphere1)
		gl.popMatrix()
	})
	brepPoints && brepPoints.forEach(p => {
		const color = hoverHighlight == p ? 0x0adfdf : (selected.contains(p) ? 0xff0000 : 0xcccc00)
		gl.pushMatrix()
		gl.translate(p)
		gl.scale(2,2,2)
		shaders.singleColor.uniforms({ color: rgbToVec4(color)}).draw(meshes.sphere1)
		gl.popMatrix()
	})
}
const CURVE_PAINTERS = {}
CURVE_PAINTERS[EllipseCurve.name] = function paintEllipseCurve(ellipse, color, startT, endT, width = 2) {
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
CURVE_PAINTERS[BezierCurve.name] = function paintBezierCurve(curve, color, startT, endT, width = 2, normal = V3.Z) {
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
CURVE_PAINTERS[L3.name] = function (curve, color, startT, endT, width = 2, normal = V3.Z) {
	gl.pushMatrix()
	let a = curve.at(startT), b = curve.at(endT)
	let ab = b.minus(a), abT = ab.getPerpendicular().normalized()
	let m = M4.forSys(ab, abT, ab.cross(abT).normalized(), a)
	gl.multMatrix(m)
	gl.scale(1, width, width)
	shaders.singleColor.uniforms({
		color: rgbToVec4(color), // TODO: error checking
	}).draw(meshes.pipe)

	gl.popMatrix()
}
var /**@type GL.Mesh */ mesh1, mesh2
var meshes = {}
function paintScreen () {
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);


	gl.loadIdentity();
	drawPlanes();
	/*
	 gl.translate(0, 0, -5);
	 gl.rotate(30, 1, 0, 0);
	 gl.rotate(angle, 0, 1, 0);
	 */
	drawVectors()

	drawPoints()


	// paintEllipseCurve(EllipseCurve.forAB(200, 300), COLORS.RD_FILL, 0, 2 * Math.PI)
	// paintBezierCurve(BezierCurve.EX3D.scale(100, 100, 100), COLORS.PP_STROKE, -2, 3)


// 	if (!mesh1) {
// //		let pF = modelBREP.faces.filter(face => face.constructor == RotationFace)[1]
// 		let pf = new RotationFace(
// 			new ProjectedCurveSurface(new BezierCurve(V3(142.87578921496748, -191.46078243076332, 0), V3(161.78547089700214, -252.13248349581008, 0), V3(284.63214994898954, -163.59789158697575, 0), V3(372.40411211189405, -210.3992206435476, 0)), V3(0, 0, 1), 0, 1), [
// 			new PCurveEdge(new BezierCurve(V3(142.87578921496748, -191.46078243076332, 0), V3(161.78547089700214, -252.13248349581008, 0), V3(284.63214994898954, -163.59789158697575, 0), V3(372.40411211189405, -210.3992206435476, 0)), V3(372.40411211189405, -210.3992206435476, 0), V3(142.87578921496748, -191.46078243076332, 0), 1, 0, new PCurveEdge(new BezierCurve(V3(142.87578921496748, -191.46078243076332, 0), V3(161.78547089700214, -252.13248349581008, 0), V3(284.63214994898954, -163.59789158697575, 0), V3(372.40411211189405, -210.3992206435476, 0)), V3(142.87578921496748, -191.46078243076332, 0), V3(372.40411211189405, -210.3992206435476, 0), 0, 1, null, V3(142.87578921496748, -191.46078243076332, 0), V3(372.40411211189405, -210.3992206435476, 0)), V3(-372.40411211189405, 210.3992206435476, 0), V3(-142.87578921496748, 191.46078243076332, 0)),
// 			StraightEdge.throughPoints(V3(142.87578921496748, -191.46078243076332, 0), V3(142.87578921496748, -191.46078243076332, -100)),
// 			new PCurveEdge(new BezierCurve(V3(142.87578921496748, -191.46078243076332, -100), V3(161.78547089700214, -252.13248349581008, -100), V3(284.63214994898954, -163.59789158697575, -100), V3(372.40411211189405, -210.3992206435476, -100)), V3(142.87578921496748, -191.46078243076332, -100), V3(372.40411211189405, -210.3992206435476, -100), 0, 1, null, V3(142.87578921496748, -191.46078243076332, 0), V3(372.40411211189405, -210.3992206435476, 0)),
// 			StraightEdge.throughPoints(V3(372.40411211189405, -210.3992206435476, -100), V3(372.40411211189405, -210.3992206435476, 0))], [])
// 		sphere = new GL.Mesh({lines: true, normals: true})
// 		sphere.faceIndexes = new Map()
// 		pf.addToMesh(sphere)
// 		sphere.computeWireframeFromFlatTriangles()
// 		console.log("sphere.lines", sphere.lines, sphere.vertices, sphere.triangles)
// 		sphere.compile()
// 	}



	gl.pushMatrix()
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


	gl.pushMatrix();
	if (brepMesh) {
		// paint faces
		var faceIndex = modelBREP.faces.length;
		while (faceIndex--) {

			let face = modelBREP.faces[faceIndex]
			let faceTriangleIndexes = brepMesh.faceIndexes.get(face)
			shaders.lighting.uniforms({
				color: rgbToVec4(hoverHighlight == face ? 0xff00ff : (selected.contains(face) ? 0x00ff45 : COLORS.RD_FILL))
			}).draw(brepMesh, 'TRIANGLES', faceTriangleIndexes.start, faceTriangleIndexes.count);
			/*
			 shaders.singleColor.uniforms({
			 color: rgbToVec4(0x0000ff)
			 }).draw(brepMesh, 'LINES');
			 */
		}

		// paint edges
		brepEdges.forEach(edge => {
			const startT = min(edge.aT, edge.bT)
			const endT = max(edge.aT, edge.bT)
			const color = hoverHighlight == edge ? 0xff00ff : (selected.contains(edge) ? 0x00ff45 : 0x000000)
			CURVE_PAINTERS[edge.curve.constructor.name](edge.curve, color, startT, endT, 2)
		})
		isEdges.forEach(edge => {
			const startT = min(edge.aT, edge.bT)
			const endT = max(edge.aT, edge.bT)
			const color = hoverHighlight == edge ? 0xff00ff : (selected.contains(edge) ? 0x00ff45 : 0x000000)
			CURVE_PAINTERS[edge.curve.constructor.name](edge.curve, color, startT, endT, 2)
		})

		gl.projectionMatrix.m[11] -= 1 / (1 << 22) // prevent Z-fighting
		shaders.singleColor.uniforms({
			color: rgbToVec4(COLORS.PP_STROKE)
		}).draw(brepMesh, 'LINES');
		gl.projectionMatrix.m[11] += 1 / (1 << 22) // prevent Z-fighting
	}
	gl.popMatrix();

	const DZ = 0.1
	gl.projectionMatrix.m[11] -= DZ // prevent Z-fighting
	featureStack.filter(f => f instanceof Sketch && !f.hide).forEach(s => {
		gl.pushMatrix()
		paintSegments(s)
		gl.popMatrix()
	})
	gl.projectionMatrix.m[11] += DZ
}

function disableConsole() {
	oldConsole = console.log;
	console.log = function() {};
}
function enableConsole() {
	if (oldConsole) {
		console.log = oldConsole;
	}
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
function drawPlanes() {
	planes.forEach(function (plane) {
		gl.pushMatrix();
		gl.multMatrix(M4.forSys(plane.right, plane.up, plane.normal));
		gl.translate(plane.rightStart, plane.upStart, plane.w);
		gl.scale(plane.rightEnd - plane.rightStart, plane.upEnd - plane.upStart, 1);

		let shader = hoverHighlight == plane ? shaders.singleColorHighlight : shaders.singleColor;

		shader.uniforms({
			color: rgbToVec4(plane.color)
		}).draw(meshes.xyLinePlane, 'LINES')

		gl.popMatrix()
	})
}

function segmentIntersectsRay(a, b, p, dir) {
	var ab = b.minus(a)
	var abXdir = ab.cross(dir)
	var div = abXdir.lengthSquared()
	var anchorDiff = p.minus(a)
	var s = anchorDiff.cross(dir).dot(abXdir) / div
	var t = anchorDiff.cross(ab).dot(abXdir) / div
	//console.log(segmentIntersectsRay, a, b, "ab", ab, "p", p, "dir", dir, s > 0 && t / div >= 0 && t / div <= 1, "s", s, "t", t, "div", div)
	return t > 0 && s >= 0 && s <= 1
}
function template(templateName, map) {
	var html = $(templateName).text;
	for (var key in map) {
		html = html.replace(new RegExp('\\$'+key, 'g'), map[key]);
	}
	return $(new Element('div', {html: html.trim()}).firstChild)

}
function updateFeatureDisplay() {
	var div = $('featureDisplay')
	div.erase('text')
	featureStack.forEach(function (feature) {
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
		newChild.inject(div);
		newChild.getElement('[name=delete]').onclick = function () {
			featureStack.remove(feature)
			let correspondingMode = modeStack.find(mode => mode.feature && mode.feature == feature)
			correspondingMode && modeEnd(correspondingMode)
			updateFeatureDisplay()
		}
		newChild.featureLink = feature
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
	selected.forEach(function (sel) {
		// TODO, if only necessary part of model is rebuilt, this probably wont be necessary:
		assert(!sel.get || sel.get(), sel)
		if (sel instanceof NameRef) sel = sel.get()
		let name = sel.name || modelBREP && modelBREP.vertexNames && modelBREP.vertexNames.get(sel)
		var newChild = template("template", {what: sel.constructor.name, name: name});
		newChild.onmouseover = function (e) {
			hoverHighlight = sel;
			paintScreen();
		};
		newChild.onmouseout = function (e) {
			hoverHighlight = null
			paintScreen();
		};
		newChild.getElement(".remove").onclick = function (e) {
			(sel instanceof SegmentEndPoint ? sel.line : sel).removeFromSketch(editingSketch)
			updateSelected()
			paintScreen()
		};
		newChild.getElement(".unsel").onclick = function (e) {
			selected.remove(sel)
			updateSelected()
			paintScreen()
		};
		sel.toBrepEdge && newChild.grab(new Element('span', {text: sel.toBrepEdge().curve.toSource(x => Math.round10(x, -3)), style: 'font-size: small;'}))
		sel.surface && newChild.grab(new Element('textarea', {text: sel.sce, style: 'font-size: xx-small;display:block;width:100%;'}))
		newChild.inject(div)
	});
	div = $('selectedConstraints')
	div.erase('text');
	(MODES.SKETCH == modeGetCurrent()) && selected.map((el) => editingSketch.getConstraintsFor(el)).concatenated().unique().forEach(function (cst) {
		var newChild
		if ('pointDistance' == cst.type
			|| 'pointLineDistance' == cst.type
			|| 'pointPlaneDistance' == cst.type) {
			newChild = template('templateDistance', {name: cst.type, id: cst.id});
			newChild.getElement('.distanceInput').value = cst.distance;
			newChild.getElement('.distanceInput').onchange = function (e) {
				cst.distance = e.target.value;
				rebuildModel();
				paintScreen();
			}
		} else if ('angle' == cst.type) {
			newChild = template('templateAngle', {name: cst.type, id: cst.id});
			var input = newChild.getElement('.distanceInput');
			newChild.getElement('.distanceInput').value = round(rad2deg(cst.value, 5));
			newChild.getElement('.fa').onclick = () => { cst.f[0] *= -1; input.value = round(rad2deg(abs(cst.value = (-PI + cst.value) % (2 * PI)))); paintScreen() };
			newChild.getElement('.fb').onclick = () => { cst.f[1] *= -1; input.value = round(rad2deg(abs(cst.value = (-PI + cst.value) % (2 * PI)))); paintScreen() };
			newChild.getElement('.fv').onclick = () => { input.value = round(rad2deg(abs(cst.value -= sign(cst.value) * 2*PI))); paintScreen() };
			input.onchange = function (e) {
				cst.value = deg2rad(e.target.value);
				rebuildModel();
				paintScreen();
			}
		} else {
			newChild = template("templateConstraint", {name: cst.type});

			cst.cs.forEach(function (el) {
				var subChild = template("templateConstraintSub", {what: el.constructor.name, name: el.name});
				subChild.inject(newChild);
				subChild.getElement(".removeFromConstraint").onclick = function (e) {
					removeFromConstraint(el, editingSketch, cst)
					updateSelected();
					paintScreen();
				};
				subChild.onmouseover = function (e) {
					hoverHighlight = el
					e.stopPropagation()
					paintScreen()
				};
				subChild.onmouseout = function (e) {
					hoverHighlight = null
					paintScreen()
				};
			});
		}
		newChild.getElement(".remove").onclick = function (e) {
			deleteConstraint(editingSketch, cst)
			updateSelected();
		}
		newChild.onmouseover = function (e) {
			hoverHighlight = cst
			paintScreen()
		}
		newChild.onmouseout = function (el) {
			hoverHighlight = null
			paintScreen()
		}
		newChild.inject(div);
	});
}
function NameRef(name) {
	var x
	if (x = NameRef.pool.get(name)) return x
	this.ref = name
}
Object.defineProperty(NameRef.prototype, "what", { get: function () { return this.get() instanceof Face ? "face" : "plane" } });
Object.defineProperty(NameRef.prototype, "name", { get: function () { return this.get().name } });
Object.defineProperty(NameRef.prototype, "plane", { get: function () { return this.getPlane() } });
NameRef.pool = new Map()
NameRef.prototype.isRef = true
NameRef.prototype.get = function () {
	return planes.find(plane => plane.name == this.ref)
		|| modelBREP && modelBREP.faces.find(face => face.name == this.ref)
		|| modelBREP && modelBREP.faces.firstMatch(face => face.getAllEdges().find(edge => edge.name == this.ref))
		|| modelBREP && mapReverse(modelBREP.vertexNames, this.ref)
}
NameRef.prototype.getPlane = function () {
	return planes.find(plane => plane.name == this.ref) || modelBREP && modelBREP.faces.find(face => face.name == this.ref).plane
}
function mapReverse(map, value) {
	let it = map.keys(), key
	while (key = it.next().value) {
		if (map.get(key) == value) {
			return key
		}
	}
}
NameRef.forObject = function (o) {
	if (o instanceof V3) {
		return new NameRef(modelBREP.vertexNames.get(o))
	} else {
		assert(o.name)
		return new NameRef(o.name)
	}
}
function initTips() {
	console.log($$('[data-tooltip]'))
	let tips = new Tips($$('[data-tooltip]'), {
		title: 'data-tooltip',
		showDelay: 0,
		fixed: true
			})
}
function initMeshes() {
	meshes.sphere1 = GL.Mesh.sphere(2)
	meshes.segment = GL.Mesh.plane({startY: -0.5, height: 1, detailX: 128})
	meshes.text = GL.Mesh.plane()
	meshes.vector = GL.Mesh.rotation([V3.ZERO, V3(0, 0.05, 0), V3(0.8, 0.05), V3(0.8, 0.1), V3(1, 0)], L3.X, Math.PI * 2, 8, true)
	meshes.pipe = GL.Mesh.rotation(NLA.arrayFromFunction(128, i => V3.create(i / 127, -0.5, 0)), L3.X, Math.PI * 2, 8, true)
	meshes.xyLinePlane = GL.Mesh.plane()
}
function initShaders() {
	shaders.singleColor = new GL.Shader(vertexShaderBasic, fragmentShaderColor);
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
function main() {
	 let o = {

	 }
	 o[EllipseCurve.name] = "a"
	console.log(Object.keys(o))

	// initTips()
	// new Request({url: 'src/testshader.glsl', async: false, onSuccess: text => console.log(text)}).send()

	modePush(MODES.DEFAULT)

	$$('.sketchControl').set('disabled', true)

	gl = GL.create({});
	gl.fullscreen();
	gl.canvas.oncontextmenu = () => false;

	setupCamera();
	//gl.cullFace(gl.FRONT_AND_BACK);
	gl.clearColor(1.0, 1.0, 1.0, 0.0)
	gl.enable(gl.BLEND)
	gl.enable(gl.DEPTH_TEST)
	// gl.enable(gl.CULL_FACE)
	gl.depthFunc(gl.LEQUAL)
	gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA); // TODO ?!


	// let s1 = new ProjectedCurveSurface(BezierCurve.EX3D.scale(100, 100, 100).project(P3.XY), V3.Z, -1, 2)
	// let s2 = new ProjectedCurveSurface(BezierCurve.EX2D.scale(100, 100, 100), V3.Z, -1, 2).rotateY(Math.PI/2).translate(0, 260, 0)
	mesh1 = new GL.Mesh()
	// mesh2 = s2.toMesh()
	// let isc = s2.isCurvesWithSurface(s1)[0]
	// mesh1 = face.inB2().toMesh()
	// mesh1.computeWireframeFromFlatTriangles()
	// mesh1.compile()
	// let mint = -2, maxt = 2
	// drPs.push(curve.at(mint), curve.at(maxt))
	// let line = L3(V3(-1560.8950828838565, 716.07295580975, 249.61382611323648), V3(0.9130103135570956, -0.36545647611595106, -0.18125598308272678))
	// isc.debugToMesh(mesh1, 'curve1')
	// console.log(mesh1)
	mesh1.compile()
	// drPs.pushAll(isc.rootPoints().concatenated())
	// mesh2 = curve.getAABB(mint, maxt).toMesh()
	// console.log(curve.getAABB().toSource())
	// console.log(mesh2)
	//mesh1 = mesh1.transform(M4.FOO.as3x3())

	let linksHaendig = M4.forSys(V3.X, V3.Y, V3.Z.negated())
	let v0 = V3(2, 3, 4), v1 = V3(-2, 3, 5), v2 = v0.cross(v1)
	let v2t = linksHaendig.transformVector(v0).cross(linksHaendig.transformVector(v1))
	let v2s = linksHaendig.inversed().transposed().transformVector(v2)
	console.log("ASKDKJALDS", v2t.dot(v2s), v2t.isParallelTo(v2s))
	initMeshes()
	initShaders()
	//	console.log(mesh.vertices)
	if (window.location.hash && window.location.hash.length > 1) {
		var key = window.location.hash.substr(1)
		console.log(key)
		load(key)
	} else {
		initModel()
	}
	let b = editingSketch && editingSketch.elements.find(el => el instanceof SketchBezier).toBrepEdge().curve
	console.log("BBB", b)

	window.onkeypress = function (e) {
		if ("Delete" == e.key) {
			selected.map((x) => x instanceof SegmentEndPoint ? x.line : x).unique().forEach((x) => x.removeFromSketch(editingSketch))
			paintScreen()
		}
	}
	gl.onclick = function (e) {
		//noinspection JSBitwiseOperatorUsage
		if (1 == e.button) {
			modePop()
		}
	}
	gl.onmouseup = function (e) {
		// don't put modeGetCurrent().mouseup in local var as 'this' won't be bound correctly
		modeGetCurrent().mouseup && modeGetCurrent().mouseup(e, getMouseLine(e))
		paintScreen()
	}
	gl.onmousedown = function (e) {
		//noinspection JSBitwiseOperatorUsage
		if (!(e.buttons & 1)) {
			return
		}

		let mouseLine = getMouseLine(e)
		console.log("mouseLine", getMouseLine(e).toString(x => x), "mode", modeGetCurrent())
		modeGetCurrent().mousedown(e, mouseLine)
	}

	gl.onmousemove = function (e) {
		var mouseLine = getMouseLine(e)
		modeGetCurrent().mousemove(e, mouseLine)
		{
			let pp, html = '', closestP = Infinity
			drPs.forEach(info => {
				let p = info.p || info
				let text = p.toString(x => Math.round10(x, -4)) + (info.p ? ' ' + info.text : '')
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
			} else {
				tooltipHide()
			}
		}
		if (e.dragging) {

			//noinspection JSBitwiseOperatorUsage
			if (e.buttons & 4) {
				// pan
				let moveCamera = V3(-e.deltaX * 2 / gl.canvas.width, e.deltaY * 2 / gl.canvas.height, 0);
				let inverseProjectionMatrix = gl.projectionMatrix.inversed();
				let worldMoveCamera = inverseProjectionMatrix.transformVector(moveCamera);
				eyePos = eyePos.plus(worldMoveCamera);
				eyeFocus = eyeFocus.plus(worldMoveCamera);
				setupCamera();
			}
			// scene rotation
			//noinspection JSBitwiseOperatorUsage
			if (e.buttons & 2) {
				let rotateLR = deg2rad(-e.deltaX / 6.0);
				let rotateUD = deg2rad(-e.deltaY / 6.0);

				// rotate
				let matrix = M4.rotationLine(eyeFocus, eyeUp, rotateLR)
				//let horizontalRotationAxis = eyeFocus.minus(eyePos).cross(eyeUp)
				let horizontalRotationAxis = eyeUp.cross(eyePos.minus(eyeFocus))
				matrix = matrix.times(M4.rotationLine(eyeFocus, horizontalRotationAxis, rotateUD))
				eyePos = matrix.transformPoint(eyePos)
				eyeUp = matrix.transformVector(eyeUp)

				setupCamera();
			}
		}
		paintScreen()
	}
	$(gl.canvas).addEvent('mousewheel', function (e) {
		//console.log(e)
		zoomFactor *= pow(0.9, -e.wheel)
		var mouseCoords = e.client
		var moveCamera = V3(mouseCoords.x * 2 / gl.canvas.width - 1, -mouseCoords.y * 2 / gl.canvas.height + 1, 0).times(1 - 1 / pow(0.9, -e.wheel))
		var inverseProjectionMatrix = gl.projectionMatrix.inversed()
		var worldMoveCamera = inverseProjectionMatrix.transformVector(moveCamera)
		//console.log("moveCamera", moveCamera)
		//console.log("worldMoveCamera", worldMoveCamera)
		eyePos = eyePos.plus(worldMoveCamera)
		eyeFocus = eyeFocus.plus(worldMoveCamera)
		setupCamera()
		paintScreen()
	})
	$("clearmode").onclick = modePop



		rebuildModel() // necessary to init planes
		updateFeatureDisplay()
		rebuildModel() // so warning will show
		initLoadSave()

	mesh1.compile()
	paintScreen();
	let lastPlane = featureStack.slice().reverse().find(f => f instanceof PlaneDefinition)
	lastPlane && modePush(MODES.PLANE_DEFINITION, lastPlane)

}

/**
 *
 * @param {L3} mouseLine
 * @param {String} consider
 */
function getHovering(mouseLine, ...consider) {
	var hoverHighlight = null, nearest = Infinity
	function checkEl(el, distance) {
		if (distance < nearest) {
			nearest = distance;
			hoverHighlight = el
		}
	}
	if (consider.contains('faces') && modelBREP) {
		modelBREP.faces.forEach((face) => {
			checkEl(face, face.intersectsLine(mouseLine))
		})
	}
	if (consider.contains('planes')) {
		planes.forEach(plane => checkEl(plane, plane.distanceTo(mouseLine)))
	}
	if (consider.contains('sketchElements')) {
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
					checkEl(closestElement, mouseLine.pointLambda(mouseLineIS) - 0.1)
				}
			}
		})
	}
	if (consider.contains('brepPoints')) {
		brepPoints.forEach(p => {
			let t = mouseLine.pointLambda(p)
			if (mouseLine.at(t).distanceTo(p) < 20) {
				checkEl(p, t - 0.1)
			}
		})
	}
	if (consider.contains('brepEdges')) {
		let projPlane = P3(mouseLine.dir1, 0)
		let projPoint = projPlane.projectedPoint(mouseLine.anchor)
		brepEdges.forEach(edge => {
			let curve = edge.curve
			const prio = 0.05
			if (curve instanceof L3) {
				if (curve.dir1.isParallelTo(mouseLine.dir1)) {
					let d = mouseLine.distanceToPoint(edge.a)
					let t = mouseLine.pointLambda(edge.a)

					if (d < 16) {
						checkEl(edge, t - prio)
					}
				} else {
					let projAnchor = projPlane.projectedPoint(curve.anchor)
					let projDir = projPlane.projectedVector(curve.dir1)
					let tCurve = projPoint.minus(projAnchor).dot(projDir) / projDir.lengthSquared()
					tCurve = edge.clampedT(tCurve)
					let p = curve.at(tCurve)
					let t = mouseLine.pointLambda(p)
					if (mouseLine.at(t).distanceTo(p) < 16) {
						checkEl(edge, t - prio)
					}
				}
			} else {
				let projCurve = curve.project(projPlane)
				let tCurve = projCurve.closestTToPoint(projPoint)
				tCurve = edge.clampedT(tCurve)
				let p = curve.at(tCurve)
				let t = mouseLine.pointLambda(p)
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
	let ndc1 = V3(pos.x * 2 / gl.canvas.width - 1, -pos.y * 2 / gl.canvas.height + 1, 0)
	let ndc2 = V3(pos.x * 2 / gl.canvas.width - 1, -pos.y * 2 / gl.canvas.height + 1, 1)
	//console.log(ndc)
	let inverseProjectionMatrix = gl.projectionMatrix.inversed()
	let s = inverseProjectionMatrix.transformPoint(ndc1)
	let dir = inverseProjectionMatrix.transformPoint(ndc2).minus(s)
	return L3.anchorDirection(s, dir)
}


function setupSelectors(el, feature, mode) {

	el.getElements('.face-select')
		.removeEvents()
		.addEvent('click', function (e) {
			var el = this
			this.addClass('selecting')
			el.set('text', 'Click on a face')
			modePush(MODES.SELECT_FACE, function (face) {
				el.removeClass('selecting')
				el.set('text', face.name)

				feature[el.dataset.featureProperty] = face.name
				rebuildModel()
			})
		})
		.each(el => el.innerHTML = feature[el.dataset.featureProperty])
		.removeClass('selecting')

	el.getElements('.plane-select')
		.removeEvents()
		.addEvent('click', function (e) {
			var el = this
			this.addClass('selecting')
			el.set('text', 'Click on a plane')
			modePush(MODES.PLANE_SELECT, function (plane) {
				console.log("plane-select callback", plane)
				el.removeClass('selecting')
				el.set('text', plane.name)

				feature[el.dataset.featureProperty] = plane.name
				rebuildModel()
			})
		})
		.each(el => el.innerHTML = feature[el.dataset.featureProperty])
		.removeClass('selecting')

	el.getElements('.segment-select')
		.removeEvents()
		.addEvent('click', function (e) {
			var el = this
			this.addClass('selecting')
			el.set('text', 'Click on a sketch segment')
			modePush(MODES.SELECT_SEGMENT, function (segment) {
				el.removeClass('selecting')
				el.set('text', segment.name)
				el.segment = segment

				feature[el.dataset.featureProperty] = segment.name
				rebuildModel()
			})
		})
		.addEvent('mouseover', function (e) {
			if (el.segment) {
				hoverHighlight = segment
				paintScreen()
			}
		})
		.each(function (el) {
			el.set('text', feature[el.dataset.featureProperty])
		})
		.removeClass('selecting')

	el.getElements('.string-id-input, .select-text')
		.removeEvents()
		.addEvent('change', function (e) {
			// TODO: check if unique
			let propName = el.dataset.featureProperty
			feature[propName] = this.value
			rebuildModel()
		})
		.each(el => el.value = feature[el.dataset.featureProperty])

	el.getElements('.select-select')
		.removeEvents()
		.addEvent('change', function (e) {
			var propName = this.dataset.featureProperty
			feature[propName] = this.value
			rebuildModel()
		})
		.each(el => el.value = feature[el.dataset.featureProperty])

	el.getElements('.dimension-input')
		.removeEvents()
		.addEvent('mousewheel', function (e) {
			console.log(e)
			let delta = (e.shift ? 1 : 10) * Math.sign(e.wheel)
			this.set('value', Math.round10(parseFloat(this.value) + delta, -6))
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
		.each(el => el.value = feature[el.dataset.featureProperty])

	el.getElement('[name=done]')
		.removeEvents()
		.addEvent('click', function () {
			modeEnd(mode)
		})
	el.getElement('[name=delete]')
		.removeEvents()
		.addEvent('click', function () {
			mode.deleteFeature()
		})

}
function chainComparison(diff1, diff2) {
	return diff1 != 0 ? diff1 : diff2;
}
function reverse(sketch) {
	sketch.varMap.forEach(function (pointIndex, point) {
		point.x = sketch.x[pointIndex]
		point.y = sketch.x[pointIndex + 1]
	});
}
function pointsFirst(a, b) {
	//console.log((a instanceof SegmentEndPoint ? -1 : 1) - (b instanceof SegmentEndPoint ? -1 : 1));
	return (a instanceof SegmentEndPoint ? -1 : 1) - (b instanceof SegmentEndPoint ? -1 : 1);
}
function clicky() {
	for (var j = 0; j < 1; j++) {
		editingSketch.gaussNewtonStep()
	}
	reverse(editingSketch)
	paintScreen()
}
function deleteConstraint(sketch, constraint) {
	if ("coincident" == constraint.type) {
		deleteCoincidence(constraint, sketch)
	} else {
		sketch.constraints.remove(constraint)
	}
	recalculate(sketch)
	paintScreen()
}
function getGroupConstraint(el, sketch, type) {
	return sketch.constraints.find(function (c) { return c.type == type && (c.fixed == el || c.segments.contains(el)); });
}
// removes segment or plane from group constraint
function removeFromGroupConstraint(el, sketch, type) {
	var groupConstraint = getGroupConstraint(el, sketch, type)
	if (!groupConstraint) return;
	groupConstraint.segments.remove(el)
	if (1 == groupConstraint.segments.length) {
		sketch.constraints.remove(groupConstraint);
	}
	if (groupConstraint.fixed == el) {
		groupConstraint.fixed = undefined;
	}
	groupConstraint.cs.remove(el)
}
function removeFromConstraint(el, sketch, constraint) {
	switch (constraint.type) {
		case "coincident":
			removeFromCoincidence(el, sketch)
			break;
		case "parallel":
		case "colinear":
		case "equalLength":
			removeFromGroupConstraint(el, sketch, constraint.type)
			break;
		case "perpendicular":
		case "pointDistance":
		case "pointLineDistance":
		case "pointOnLine":
		case "angle":
			sketch.constraints.remove(constraint)
			break
		default:
			throw new Error("implement!" + constraint.type)
			//throw new Error("unknown constraint " + constraint.type);
		// console.log("errror", constraint.type);
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
	var els = selected.filter((el) => el instanceof SketchLineSeg || ('equalLength' != type && el.plane) )
	if (els.length < 2) return;

	var newGroup = els.map(el => {
		el = el.plane ? new NameRef(el.name) : el
		var c = getGroupConstraint(el, editingSketch, type);
		return c
			? (c.fixed ? c.segments.concat(c.fixed) : c.segments)
			: el }).concatenated().unique();
	var fixeds = newGroup.filter(el => el.isRef);
	if (1 < fixeds.length) {
		throw new Error("cannot have two fixed");
	}

	var oldConstraints = els.map(el => getGroupConstraint(el, editingSketch, type)).filter(x => x);
	editingSketch.constraints.removeAll(oldConstraints);

	// add new constraint
	// fixeds[0] may be null
	var segments = newGroup.filter(x => x instanceof SketchLineSeg);
	var f = fixeds[0];
	editingSketch.constraints.push(new Constraint(type, f ? segments.concat(f) : segments, {fixed: fixeds[0], segments: segments}));

	rebuildModel();
	paintScreen();
	updateSelected();
}
function Constraint(type, constrains, props) {
	if (constrains.constructor != Array) {
		throw new Error("not an array: " + constrains);
	}
	this.type = type;
	this.id = globalId++;
	this.cs = constrains;
	for (var key in props) {
		this[key] = props[key];
	}
}
Constraint.prototype = {
	constrains: function (o) {
		if (!this.cs) {
			console.log(this);
		}
		return this.cs.contains(o);
	},
	segmentOtherTypeFree(sketch, el) {

	},
	/*
	serialize: function (els) {
		function lookup (val) {
			let i = els.indexOf(val)
			return i != -1 ? 'els[' + i + ']' : serialize(val)
		}
		return serialize(this, lookup)
		return `new Constraint('${this.type}',[${els.map(lookup).join(',')}])`
	}
	*/
	//constrains: o => this.cs.contains(o)
}
Constraint.prototype.constructor = Constraint


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
	function transform(v, first) {
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
				: assert(v.constructor && v.constructor.name && v.constructor.name != "Object",
					() => (console.log(v), v.toSource() +v.constructor.name)) && {'#PROTO': v.constructor.name}
				let keys = Object.keys(v)
				for (let i = 0; i < keys.length; i++) {
					result[keys[i]] = transform(v[keys[i]])
				}
				return result
			}
		} else {
			throw new Error("?" + typeof v + v.toString())
		}
	}

	let visited = new Set()
	let listMap = new Map(), resultList = []
	listMap.set(v, 0)
	resultList.push(v)
	gatherList(v)
	console.log(resultList)

	resultList = resultList.map(v => transform(v, true))
	console.log(JSON.encode(resultList))
	return JSON.encode(resultList)
}

function unserialize(string) {
	function fixObjects(v) {
		if (v && v.constructor === Array) {
			for (let i = 0; i < v.length; i++) {
				v[i] = fixObjects(v[i])
			}
			return v
		} else if ('object' == typeof v && null != v) {
			if ("#PROTO" in v) {
				assert(window[v["#PROTO"]], v["#PROTO"] + ' Missing ' + window[v["#PROTO"]])
				let result = Object.create(window[v["#PROTO"]].prototype)
				let keys = Object.keys(v)
				for (let i = 0; i < keys.length; i++) {
					//if ("name" == keys[i]) console.log(result)
					if ("#PROTO" != keys[i]) {
						Object.defineProperty(result, keys[i], {value: fixObjects(v[keys[i]]), enumerable: true, writable: true})
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
			if ("#REF" in v) {
				return tree[v["#REF"]]
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
	let tree = JSON.decode(string, true)
	// console.log(tree)
	fixObjects(tree)
	// console.log(tree)
	linkReferences(tree)
	// console.log(tree)
	return tree[0]
}

function makeAngle() {
	var selSegments = selected.filter((el) => el instanceof SketchLineSeg || el.plane )
	console.log(selSegments);
	if (selSegments.length == 2) {
		var segment = selSegments[0] instanceof SketchLineSeg ? selSegments[0] : selSegments[1];
		var other = selSegments[0] instanceof SketchLineSeg ? selSegments[1] : selSegments[0];
		if (! (segment instanceof SketchLineSeg)) {
			throw new Error ("at least one must be a segment");
		}
		editingSketch.constraints.push(new Constraint("angle", [segment, other], {segment:segment, other:other, f: [1, 1], value: selSegments[0].angleTo(selSegments[1])})); //
		rebuildModel();
		paintScreen();
		updateSelected();
	}
}
function makePerpendicular() {
	var selSegments = selected.filter((el) => el instanceof SketchLineSeg || el.plane )
	if (selSegments.length == 2) {
		var segment = selSegments[0] instanceof SketchLineSeg ? selSegments[0] : selSegments[1];
		var other = selSegments[0] instanceof SketchLineSeg ? selSegments[1] : selSegments[0];
		if (! (segment instanceof SketchLineSeg)) {
			throw new Error ("at least one must be a segment");
		}
		editingSketch.constraints.push(new Constraint("perpendicular", [segment, other], {segment:segment, other:other}));
		rebuildModel();
		paintScreen();
		updateSelected();
	}
}
function makeDistance() {
	if (2 != selected.length) return
	var point = selected[0] instanceof SegmentEndPoint ? selected[0] : selected[1]
	var other = selected[0] instanceof SegmentEndPoint ? selected[1] : selected[0]
	if (point instanceof SegmentEndPoint && (other instanceof SketchLineSeg || other instanceof SegmentEndPoint || other.plane)) {
		var newConstraint
		if (other instanceof SegmentEndPoint) {
			newConstraint = new Constraint("pointDistance", selected.slice(), {distance:round(other.distanceToCoords(point))})
		} else if (other instanceof SketchLineSeg) {
			let distance = round(other.distanceToCoords(point))
			newConstraint = new Constraint("pointLineDistance", [point, other],
				{point:point, other:other, distance: distance})
		} else {
			let distance = round(other.plane.intersectionWithPlane(editingSketch.plane)
				.transform(editingSketch.worldToSketchMatrix).distanceToPoint(V3(point)))
			other = new NameRef(other.name)
			newConstraint = new Constraint("pointPlaneDistance", [point, other],
				{point:point, other:other, distance: distance})
			console.log(newConstraint)
		}
		editingSketch.constraints.push(newConstraint)
		rebuildModel()
		paintScreen()
		updateSelected()
		$("distanceInput"+newConstraint.id).select()
	}
}
function selPointOnLine() {
	if (2 != selected.length) return;
	var point = selected[0] instanceof SegmentEndPoint ? selected[0] : selected[1];
	var other = selected[0] instanceof SegmentEndPoint ? selected[1] : selected[0];
	if (point instanceof SegmentEndPoint && (other instanceof SketchLineSeg || other instanceof SketchBezier || other instanceof SketchArc || other.plane)) {
		other = !(other.plane) ? other : new NameRef(other.name)
		var newConstraint = new Constraint("pointOnLine", [point, other], {point:point, other:other});
		console.log(newConstraint);
		editingSketch.constraints.push(newConstraint);
		rebuildModel();
		paintScreen();
		updateSelected();
	}
}
function makeSelCoincident() {
	var selPoints = selected.filter(function (s) { return s instanceof SegmentEndPoint; });
	if (selPoints.length >= 2) {
		for (var i = 1; i < selPoints.length; i++) {
			makeCoincident(selPoints[0], selPoints[i], editingSketch);
		}
		rebuildModel();
		paintScreen();
		updateSelected();
	}
}

function faceSketchPlane() {
	var plane = editingSketch.plane
	var viewLine = L3.throughPoints(eyePos, eyeFocus)
	eyeFocus = plane.intersectionWithLine(viewLine) || viewLine.closestPointToPoint(V3.ZERO)
	eyePos = eyeFocus.plus(plane.normal.times(100));
	eyeUp = plane.up;
	setupCamera();
	paintScreen();
}
function Extrude(name, segmentName, start, end, operation) {
	this.type = "extrude"
	this.name = name || "extrude" + (globalId++)
	this.segmentName = segmentName
	this.start = start || 0
	this.end = end || 100
	this.operation = operation || "minus"
}
Extrude.prototype = {}
Extrude.prototype.constructor = Extrude
class Pattern {
	constructor(name, features) {
		this.name = name || 'pattern' + globalId
		this.features = features
		this.direction = V3.X
		this.count = 2
		this.totalLength = 20
		this.intervalLength = 10
		this.mode = 0
	}
	build() {

	}
}
function addSketch() {
	let feature = new Sketch()
	featureStack.push(feature)
	modePush(MODES.SKETCH, feature)
	updateFeatureDisplay()
}
function addPlane() {
	let feature = new PlaneDefinition("face")  // TODO?
	featureStack.push(feature)
	modePush(MODES.PLANE_DEFINITION, feature)
	updateFeatureDisplay()
}
function addExtrude() {
	let feature = new Extrude();
	featureStack.push(feature);
	modePush(MODES.EXTRUDE, feature)
	updateFeatureDisplay()
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
	modeStack.push(mode)
	mode.init.apply(mode, args)
	modeUpdateDisplay()
}
function modeEnd(mode) {
	if (modeStack.contains(mode)) {
		do {
			var popped = modeStack.pop()
			popped.end()
		} while (popped != mode)
		modeUpdateDisplay()
	}
}
function modePop() {
	modeStack.pop().end()
	modeUpdateDisplay()
}
function modeGetCurrent() {
	return modeStack.last()
}















function createArcMesh(insideRadius, outsideRadius, cornerCount) {
	var mesh = new GL.Mesh({});
	for (var i = 0; i < cornerCount; i++) {
		mesh.vertices.push(
			V3(outsideRadius, i / (cornerCount - 1), 0),
			V3(insideRadius, i / (cornerCount - 1), 0))
		if (i == cornerCount - 1) break;
		mesh.triangles.pushAll([
			0, 2, 1,
			1, 2, 3].map(offset=>(2 * i + offset) % (2 * cornerCount)))
	}
	mesh.compile();
	return mesh;
}


function tooltipShow(htmlContent, x, y) {
	let tipWrap = $('tip-wrap')
	assert(tipWrap)
	let tip = $('tip')
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
		tipWrap.setStyle('top', (y + 16+8) + 'px')
	}


}
function tooltipHide() {
	$('tip-wrap').setStyle('visibility', 'hidden')
}

