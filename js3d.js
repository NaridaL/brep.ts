"use strict";
window.onerror = function (errorMsg, url, lineNumber, column, errorObj) {
	console.log(errorMsg, url, lineNumber, column, errorObj);
}

Object.getOwnPropertyNames(Math).forEach(function (propertyName) {
	/*if (window[propertyName]) {
	 throw new Error("already exists"+propertyName)
	 }*/
	window[propertyName] = Math[propertyName];
});
var M4 = NLA.Matrix4x4, V3 = NLA.Vector3, P3 = NLA.Plane3, L3 = NLA.Line3
var DEBUG = true

/*
var S0 = 0, S1 = 1, D0 = 2, D1 = 3;
Array.prototype.concatenated = function () {
	return Array.prototype.concat.apply([], this);
}
Array.prototype.pushAll = function (arr) {
	Array.prototype.push.apply(this, arr);
};
Array.prototype.max = function() {
	return Math.max.apply(null, this);
};

Array.prototype.min = function() {
	return Math.min.apply(null, this);
};
Array.prototype.unique = function (o) {
	return this.filter(function(item, pos, array) {
		return array.indexOf(item) == pos;
	});
}
Array.prototype.isEmpty = function () {
	return 0 == this.length;
};
Array.prototype.removeAll = function (os) {
	for (var i = this.length; --i > 0;) {
		if (os.contains(this[i])) {
			this.splice(i, 1);
		}
	}
};
Array.prototype.sliceStep = function (start, step, chunkSize) {
	chunkSize = chunkSize || 1;
	var result = new Array(ceil((this.length - start) / step)); // "- start" so that chunk in the last row will also be selected, even if the row is not complete
	var index = 0;
	for (var i = 0; i < this.length; i += step) {
		for (var j = 0; j < chunkSize; j++) {
			result[index++] = this[start + i + j];
		}
	}
	return result;
}
*/
var globalId = 0;
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
}
Sketch.prototype = {
	getConstraintsFor: function (el) {
		console.log("el", el);
		return this.constraints.filter((constraint) => constraint.type != 'coincident' && constraint.constrains(el));
	},

	constrainDistancePointFixedLineWC: function (point, lineWC, distance) {
		this.b.push(distance);
		var px = this.varMap.get(point), py = px + 1;
		var lineA_SC = this.worldToSketchMatrix.transformPoint(lineWC.anchor);
		var lineB_SC = lineA_SC.plus(this.worldToSketchMatrix.transformVector(lineWC.dir1));
		console.log(lineA_SC, lineB_SC);
		this.F.push(function (x) { return distanceLinePoint(lineA_SC.x, lineA_SC.y, lineB_SC.x, lineB_SC.y, x[px], x[py]); })
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
		var ax = this.varMap.get(segment.a), ay = ax + 1, bx = this.varMap.get(segment.b), by = bx + 1;
		var lineDirectionSC = this.worldToSketchMatrix.transformVector(lineWC.dir1);
		this.F.push(function (x) { return angleVectors(x[bx] - x[ax], x[by] - x[ay], lineDirectionSC.x, lineDirectionSC.y); })
	},
	constrainEqualDistance: function (ax, ay, bx, by, cx, cy, dx, dy) {
		this.b.push(0);
		this.F.push((x) => distance(x[ax], x[ay], x[bx], x[by]) - distance(x[cx], x[cy], x[dx], x[dy]));
	},
	constrainAngleSegmentSegment: function(line1, line2, cosAngle) {
		this.b.push(cosAngle);
		var ia = this.varMap.get(line1.a),
			ib = this.varMap.get(line1.b),
			ic = this.varMap.get(line2.a),
			id = this.varMap.get(line2.b)
		this.F.push((x) => angleABCD(
			x[ia], x[ia + 1],
			x[ib], x[ib + 1],
			x[ic], x[ic + 1],
			x[id], x[id + 1]
		))
	},
	constrainDistancePointPoint: function(pA, pB, pDistance) {
		this.b.push(pDistance);
		var ia = this.varMap.get(pA), ib = this.varMap.get(pB);
		this.F.push((x) => distance(x[ia], x[ia + 1], x[ib], x[ib + 1]));
	},
	constrainDistancePointSegment: function(point, segment, distance) {
		this.b.push(distance);
		var ip = this.varMap.get(point), ia = this.varMap.get(segment.a), ib = this.varMap.get(segment.b);
		this.F.push((x) => distanceLinePoint(x[ia], x[ia + 1], x[ib], x[ib + 1], x[ip], x[ip + 1]));
	},
	constrainAngleSegmentSegment2: function(ab, cd, f, value) {
	// value is in [-2 pi ; 2 pi]; divide by 2 to map to [-pi;pi]
		this.b.push(1 / cos(value / 2), 1 / sin(value / 2));
		var ia = this.varMap.get(ab.a),
			ib = this.varMap.get(ab.b),
			ic = this.varMap.get(cd.a),
			id = this.varMap.get(cd.b);
		// calculate the angle for each segment relative to the x axis using atan
		// subtract angle of ab from angle of cd to get signed difference
		// two functions, one calculate sin, one cos of signed difference / 2 to map the signed difference to a point on the unit circle
		// using two functions is necessary so that the resulting functions are "stetig"
		this.F.push(
		(x) => {
			var angle = atan2(
					(x[id + 1] - x[ic + 1]) * f[1],
					(x[id]     - x[ic]    ) * f[1])
				- atan2(
					(x[ib + 1] - x[ia + 1]) * f[0],
					(x[ib]     - x[ia]    ) * f[0]);
			console.log("angle", angle, 1 / cos(angle / 2), sin(angle / 2));
			return 1 / cos(
					( atan2(
						(x[id + 1] - x[ic + 1]) * f[1],
						(x[id]     - x[ic]    ) * f[1])
					- atan2(
						(x[ib + 1] - x[ia + 1]) * f[0],
						(x[ib]     - x[ia]    ) * f[0])) / 2); },
		(x) => 1 / sin(
			( atan2(
				(x[id + 1] - x[ic + 1]) * f[1],
				(x[id]     - x[ic]    ) * f[1])
			- atan2(
				(x[ib + 1] - x[ia + 1]) * f[0],
				(x[ib]     - x[ia]    ) * f[0])) / 2) );
	},
	/*
	gaussNewtonStep2: function () {
		//disableConsole();
		var {F, x, b} = this
		var Fx = F.map(function (f) { return f(x)});
		console.log("x", x);
		console.log("Fx", Fx);
		console.log("b", b);
		var jacobi = $M(calcDFx(F, x, Fx));
		console.log("jacobi\n", calcDFx(F, x, Fx).toSource());
		console.log("jacobi\n", jacobi.inspect());
		var jacobiTranspose = jacobi.transpose();
		var matrix = jacobiTranspose.multiply((jacobi.multiply(jacobiTranspose)).inverse());
		console.log("inverseJacobi",matrix.inspect());
		var xDif = matrix.multiply($V(Fx).subtract($V(b)));
		console.log("xDif", xDif.inspect());
		this.x = $V(x).subtract(xDif).elements;
		Fx = F.map(function (f) { return f(x)});
		//enableConsole();
		return $V(Fx).subtract($V(b)).modulus();
	},
	*/
	gaussNewtonStep: function () {
		//disableConsole();
		var {F, x, b} = this
		var Fx = NLA.Vector.fromFunction (F.length, function (i) { return  F[i](x) })
		//console.log("x", x);
		//console.log("Fx", Fx.toString());
		//console.log("b", b);
		var jacobi = calcDFxNLA(F, x, Fx.v);
		var jacobiTranspose = jacobi.transposed();
		var matrix = jacobiTranspose.times((jacobi.times(jacobiTranspose)).inversed());
		var bVector = new NLA.Vector(new Float64Array(b))
		var xDiff = matrix.timesVector(Fx.minus(bVector));
		this.x = new NLA.Vector(new Float64Array(x)).minus(xDiff).v;
		Fx = NLA.Vector.fromFunction (F.length, i => F[i](this.x))
		//enableConsole();
		return Fx.minus(bVector).length();
	},
	getLoopForSegment: function (segment) {
		var startPoint = segment.a;
		var currentPoint = startPoint;
		var points = [];
		do {
			points.push(currentPoint);
			//console.log(currentPoint.coincidence.which.filter(function (point) { point instanceof SegmentEndPoint && point != currentPoint; }));
			var otherPointsInCoincidence = currentPoint.coincidence.which.filter(function (point) {
				//console.log("point", point, point instanceof SegmentEndPoint, point != currentPoint);
				return point instanceof SegmentEndPoint && point != currentPoint; });
			if (otherPointsInCoincidence.length != 1) {
				throw new Error("not a loop");
			}
			var nextSegmentCoincidencePoint = otherPointsInCoincidence[0]
			currentPoint = nextSegmentCoincidencePoint.line.getOtherPoint(nextSegmentCoincidencePoint);
		} while (startPoint.coincidence != currentPoint.coincidence);
		return points
	},
	delete: function () {
		$('sketchEditor').set('display', 'none')
		features.remove(this)

		rebuildModel()
		updateFeatureDisplay()
	}
};
function SegmentEndPoint(x, y, line, sketch) {
	this.x = x;
	this.y = y;
	this.line = line;
	this.id = globalId++;
	this.sketch = sketch
	//this.sceneElement = addScenePoint();
	//this.updateSceneElement();
}
SegmentEndPoint.prototype = {
	distanceToCoords: function (p) {
		return distance(this.x, this.y, p.x, p.y);
	},
	toString: function () {
		return "SegmentEndPoint #" + this.id;
	},
	V3: function () {
		return V3(this.x, this.y, 0);
	},

	isConstrained: function () {
		var line = this.line;
		return constraints.some(
			function (constraint) {
				return constraint.type == "parallel" && constraint.which.contains(line)
					|| constraint.type == "perpendicular" && (constraint.first == line || constraint.second == line)
					|| constraint.type == "pointOnLine" && (constraint.line == line || constraint.point == this);
			});
	},
	remove: function () {
		removeCoincidence(this, this.sketch);
	},
}
var globHighlight = "";
function makeCoincident(p1, p2, sketch) {
	if (p1.coincidence && p2.coincidence) {
		for (var i = 0; i < p2.coincidence.length;  i++) {
			p1.coincidence.which.push(p2.coincidence.which[i]);
			p2.coincidence.which[i].coincidence = p1.coincidence;
		}
		p1.coincidence.which.push(p2);
		p2.coincidence = p1.coincidence;
	} else if (p1.coincidence) {
		p1.coincidence.which.push(p2);
		p2.coincidence = p1.coincidence;
	} else if (p2.coincidence) {
		p2.coincidence.which.push(p1);
		p1.coincidence = p2.coincidence;
	} else {
		sketch.constraints.push(p1.coincidence = p2.coincidence = {"type":"coincident", "which":[p1,p2]});
	}
	p2.x = p1.x;
	p2.y = p1.y;
}
function removeCoincidence(p, sketch) {
	if (!p.coincidence) return;
	p.coincidence.which.remove(p);
	if (p.coincidence.which.length == 1) {
		sketch.constraints.remove(p.coincidence);
		p.coincidence.which[0].coincidence = undefined
	}
	p.coincidence = undefined;
}
function Segment(x1, y1, x2, y2) {
	this.points = [new SegmentEndPoint(x1, y1, this), new SegmentEndPoint(x2, y2, this)];
	this.a = this.points[0];
	this.b = this.points[1];
	this.id = globalId++;
	this.name = "segment" + this.id
	//this.sceneElement = addSceneSegment();
	//this.updateSceneElement();
}
Segment.prototype = {
	distanceToCoords: function (coords) {
		return this.distanceTo(coords.x, coords.y);
	},
	angleTo: function (segment) {
		return segment.angleAB() - this.angleAB();
	},
	toString: function ()  {
		return "Segment #" + this.id;
	},
	angleAB: function () {
		return atan2(this.b.y - this.a.y, this.b.x - this.a.x);
	},
	removeFromSketch: function (sketch) {
		this.a.remove();
		this.b.remove();
		sketch.elements.remove(this);
		removeFromGroupConstraint(this, sketch, "colinear");
		removeFromGroupConstraint(this, sketch, "parallel");
		removeFromGroupConstraint(this, sketch, "equalLength");
		selected.remove(this);
	},
	getOtherPoint: function (point) {
		return this.a != point ? this.a : this.b;
	},
	getS: function (v) {
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
				return distance(x, x1, y, y1);
			} else {
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
		var dist = Math.abs(a * x + b * y - c) / length(a, b);
		var xClosest = (b * (b * x - a * y) + a * c) / lengthSquared(a, b);
		var yClosest = (a * ( -b * x + a * y) + b * c) / lengthSquared(a, b);
		if (isBetween(xClosest, x1, x2)) {
			return {"x": xClosest, "y": yClosest};
		} else {
			if (x1 < x2 && x < x1 || x1 > x2 && x > x1
				|| x1 == x2 && (y1 < y2 && y < y2 || y1 > y2 && y > y1)) {
				return {x: x1, y: y1}
			} else {
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
	}
}
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
function angle(x1, y1, x2, y2) {
	return length(x2 - x1, y2 - y1);
}

function length(x, y) {
	return Math.sqrt(x * x + y * y);
}
function lengthSquared(x, y) {
	return x * x + y * y
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
			p.coincidence.which.forEach(function (p2) {
				varMap.set(p2, xIndex);
			});
		} else {
			varMap.set(p, xIndex);
		}
		x.push(p.x, p.y);
		xIndex += 2;
	};
	try {
		console.log("recalculating");
		sketch.b = [];
		sketch.x = [];
		sketch.F = [];
		sketch.varMap = new Map();
		// init x to current values
		var xIndex = 0;
		sketch.elements.forEach(function (line) {
			pushPoint(line.a);
			pushPoint(line.b);
		});
		sketch.constraints.forEach(function (c) {
			//console.log(c);
			if (c.type == "parallel" || c.type == "colinear" || c.type == "equalLength") {
				if (c.fixed) {
					for (var j = 0; j < c.segments.length; j++) {
						sketch.constrainAngleSegmentToPlane(c.segments[j], c.fixed, 1);
						if ("colinear" == c.type) {
							sketch.constrainDistancePointFixedPlane(c.segments[j].a, c.fixed, 0);
						}
					}
				} else {
					for (var j = 1; j < c.segments.length;  j++) {
						var first = c.segments[0], second = c.segments[j];
						if ("parallel" == c.type || "colinear" == c.type) {
							// assume that max. one element can be a line or a plane
							sketch.constrainAngleSegmentSegment(first, second, 1);
						}
						if ("colinear" == c.type) {
							sketch.constrainDistancePointSegment(second.a, first, 0)
						}
						if ("equalLength" == c.type) {
							sketch.constrainEqualDistance.apply(sketch,
								[0, j]
									.map((segmentsIndex) => c.segments[segmentsIndex].points).concatenated()
									.map((point) => sketch.varMap.get(point))
									.map((pointXCoord) => [pointXCoord, pointXCoord + 1]).concatenated());
						}
					}
				}
			}
			if (c.type == "pointDistance") {
				sketch.constrainDistancePointPoint(c.cs[0], c.cs[1], c.distance);
			}
			if (c.type == "pointOnLine" || c.type == "pointLineDistance") {
				var distance = c.type == "pointLineDistance" ? c.distance : 0;
				if (c.other instanceof CustomPlane) {
					sketch.constrainDistancePointFixedPlane(c.point, c.other, distance)
				} else {
					sketch.constrainDistancePointSegment(c.point, c.other, distance)
				}
			}
			if (c.type == "angle") {
				sketch.constrainAngleSegmentSegment2(c.cs[0], c.cs[1], c.f, c.value);
			}
			if (c.type == "perpendicular") {
				var cosValue = c.type == "angle" ? cos(c.value) : 0
				// assume that max. one element can be a line or a plane
				if (c.other instanceof CustomPlane) {
					sketch.constrainAngleSegmentToPlane(c.segment, c.other, cosValue);
				} else {
					sketch.constrainAngleSegmentSegment(c.segment, c.other, cosValue);
				}
				/*
				 constrainAngle2D.apply(undefined,
				 [0, 1]
				 .map((whichIndex) => c.constrains[whichIndex].points).concatenated()
				 .map((point) => varMap.get(point))
				 .map((pointXCoord) => [pointXCoord, pointXCoord + 1]).concatenated()
				 .concat(c.value));
				 */
			}
//			console.log("added constraint, b:", b);
		});
		if (sketch.b.isEmpty()) {
			return;
		}
		var lastSize = 0
		for (var count = 0; count < 2; count++) {
			if ((lastSize = sketch.gaussNewtonStep()) < 1e-6) {
				break;
			}
		}
		enableConsole()
		console.log(`broke at ${count}, lastSize ${lastSize}`)
		reverse(sketch);
	} catch (e) {
		console.log(e);
	}
}
/*
function wtfStep() {
	disableConsole();
	console.log("x", x);
	Fx = F.map(function (f) { return f(x)});
	console.log("Fx", Fx);
	console.log("b", b);
	var diffs = calcDFx2(F, x, Fx);
	console.log("diffs", diffs);
	var maxDiff = 0;
	var useXIndex;
	for (var i = 0; i < x.length; i++) {
		if (diffs[i] > maxDiff){
			maxDiff = diffs[i];
			useXIndex = i;
		}
	}
	var currLength = $V(Fx).minus($V(b)).modulus();
	var oldLength = currLength;
	var xDif = x.map(function () { return 0; });
	xDif[useXIndex] = currLength / maxDiff;
	var subValue = currLength / maxDiff;
	var newX, newLength;
	var halvedCount = -1;
	do {
		newX = x.slice();
		newX[useXIndex] -= subValue;
		subValue /= 2;
		halvedCount++;
		var newFx = F.map(function (f) { return f(newX)});
		newLength = $V(newFx).minus($V(b)).modulus();
	} while (newLength > 1e-6 && newLength >= currLength);
	x = newX;
	console.log("xdif", xDif);
	console.log("currLength", $V(Fx).minus($V(b)).modulus(), "halved", halvedCount);
	enableConsole();

	return $V(Fx).minus($V(b)).modulus();
	// DFx = DF.map(function (dfLine) { return dfLine.map(function (dff) { return dff(x) })});
	/*DFx = calcDFx(F, x, Fx);

	 var DFxM = $M(DFx);
	 console.log("DFx", DFx);

	 var xDif = [];
	 // iterate over x values
	 var minAngle = 10, changeVarIndex;
	 for (var i = 0; i < x.length; i++) {
	 var angle = $V(b).minus($V(Fx)).angleFrom($V(DFxM.col(i)));
	 console.log(DFxM.col(i), angle, i);
	 if (angle < minAngle) {
	 minAngle = angle;
	 changeVarIndex = i;
	 }
	 xDif.push(0);
	 /*
	 var sum = 0;
	 // iterate over b values / functions
	 for (var j = 0; j < b.length; j++) {
	 sum += DFx[j][i] * (Fx[j] - b[j]); // / $V(DFx[j]).modulus();
	 //console.log("sumArray", DFx[j],  lengthArray(DFx[j]));
	 //sum += (Fx[j] - b[j]) / DFx[j][i];
	 }
	 xDif.push(sum / 100);

	 }
	 xDif[changeVarIndex] = $V(b).minus($V(Fx)).modulus() - DFxM.col(i).modulus();
	 */
//}

function distanceLinePoint(x1, y1, x2, y2, x, y) {
	var a = y1 - y2;
	var b = x2 - x1;
	var c = x2 * y1 - x1 * y2;
	var dist = Math.abs(a * x + b * y - c) / length(a, b);
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
function calcDFx(F, x, Fx) {
	var DIFF = 1e-3;
	var DFx = F.map(function (a) { return new Array(x.length); });
	for (var i = 0; i < x.length; i++) {
		x[i] += DIFF;
		for (var j = 0; j < F.length; j++) {
			DFx[j][i] = (F[j](x) - Fx[j]) / DIFF;
		}
		x[i] -= DIFF;
	}
	return DFx;
}
function calcDFxNLA(F, x, Fx) {
	var DIFF = 1e-3;
	var DFx = NLA.Matrix.forWidthHeight(x.length, F.length)
	for (var colIndex = 0; colIndex < x.length; colIndex++) {
		x[colIndex] += DIFF;
		for (var rowIndex = 0; rowIndex < F.length; rowIndex++) {
			DFx.setEl(rowIndex, colIndex, (F[rowIndex](x) - Fx[rowIndex]) / DIFF)
		}
		x[colIndex] -= DIFF;
	}
	return DFx;
}
function calcDFx2(F, x, Fx) {
	var DIFF = 1e-9;
	var DFx = new Array(x.length);
	var currLength = $V(Fx).minus($V(b)).modulus();
	for (var i = 0; i < x.length; i++) {
		x[i] += DIFF;
		var newLength = $V(F.map(function (f) { return f(x); })).minus($V(b)).modulus();
		DFx[i] = (newLength - currLength) / DIFF;
		x[i] -= DIFF;
	}
	return DFx;
}

var oldConsole = undefined;
function disableConsole() {
	oldConsole = console.log;
	console.log = function() {};
}
function enableConsole() {
	if (oldConsole) {
		console.log = oldConsole;
	}
}
function paintLineXY(a, b, color, width) {
	width = width || 1;
	var ab = b.minus(a);
	if (ab.isZero()) {return; }
	var abT = ab.getPerpendicular().unit();
	//console.log(ab);
	gl.pushMatrix();
	//gl.translate(a);
	gl.multMatrix(M4.forSys(ab, abT, V3.Z, a));
	gl.scale(1, width, 1);
	renderColor(segmentMesh, color);
	gl.popMatrix();
}
function randomVec4Color(opacity) {
	opacity = opacity || 1.0
	return [Math.random(), Math.random(), Math.random(), opacity]
}
function drawArc(center, startAngle, endAngle, innerRadius, outerRadius, color) {
	gl.pushMatrix();
	gl.translate(center);
	arcShader2.uniforms({
		color: rgbToVec4(color),
		offset: startAngle,
		step: endAngle - startAngle,
		innerRadius: innerRadius,
		outerRadius: outerRadius
	}).draw(arcMesh);
	gl.popMatrix();
}
var savedTextures = {};
function getTextureForString(str) {
	if (savedTextures[str]) {
		return savedTextures[str];
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


	ctx.fillText(str, 0, 0);


	var texture = GL.Texture.fromImage(canvas, {minFilter: gl.LINEAR_MIPMAP_LINEAR, magFilter: gl.LINEAR});
	texture.textWidthPx = textWidthPx;

	savedTextures[str] = texture;

	return texture;
}
var zoomFactor = 1;
function paintSegments(sketch) {
	if (!sketch.plane) return;
	//console.log("painting segments", sketch.elements.length);
	/*ctx.clearRect (0, 0, ctx.canvas.width, ctx.canvas.height);
	 ctx.fillStyle="rgb(100, 100, 255)";
	 ctx.lineWidth=2;*/
	//console.log(sketch.sketchToWorldMatrix);
	gl.multMatrix(sketch.sketchToWorldMatrix);
	sketch.elements.forEach(function (line) {
		//console.log("line", line);
		//console.log("hoverHighlight.length", hoverHighlight.length);
		//ctx.beginPath();
		paintLineXY(line.a.V3(), line.b.V3(), colorFor(highlighted.contains(line) || hoverHighlight.contains(line), selected.contains(line)));


		//console.log(line.a);
		gl.pushMatrix();
		gl.translate(line.a.V3());
		renderColor(circleMesh, colorFor(highlighted.contains(line.a) || hoverHighlight.contains(line.a), selected.contains(line.a)));
		gl.popMatrix();

		gl.pushMatrix();
		gl.translate(line.b.V3());
		renderColor(circleMesh, colorFor(highlighted.contains(line.b) || hoverHighlight.contains(line.b), selected.contains(line.b)));
		gl.popMatrix();

	});

	paintConstraints(sketch);
	// console.log("linecount", lines.length);
}
function colorFor(highlighted, selected) {
	return !selected
		? (!highlighted ? 0x33CCFF : 0x145266)
		: (!highlighted ? 0xFF3399 : 0x330A1E)
}
function getAllPoints(sketch) {
	return [].concat.apply([], sketch.elements.map(function (segment) { return segment.points; }));
}
var highlighted = [], selected = [];
function paintConstraints(sketch) {
	var crossCount = 2;
	sketch.constraints.forEach(function (constraint) {
		switch (constraint.type) {
			case "coincident":
				gl.pushMatrix();
				gl.translate(constraint.which[0].x, constraint.which[0].y, 0);
				renderColor(ringMesh, colorFor(highlighted.contains(constraint.which[0]) || hoverHighlight.contains(constraint.which[0]), selected.contains(constraint.which[0])));
				gl.popMatrix();
				break;
			case "parallel":

				for (var c = 0; c < constraint.segments.length; c++) {
					var line = constraint.segments[c];
					var ab = line.getVectorAB();
					var abLength = ab.length();
					var abNormalized = ab.unit();
					var ab90 = ab.getPerpendicular().unit();
					var crossLength = 10;
					var crossSpacing = 3;
					for (var i = 0; i < crossCount; i++) {
						var s = abLength / 2 - crossSpacing * crossCount / 2 + i * crossSpacing;
						var crossStart = line.a.V3().plus(abNormalized.times(s)).plus(ab90.times(crossLength / 2));
						var crossEnd = line.a.V3().plus(abNormalized.times(s)).plus(ab90.times(-crossLength / 2));
						//console.log("crosspos", crossStart, crossEnd);
						paintLineXY(crossStart, crossEnd);
					}
				}
				crossCount++;
				break;
			case "perpendicular":
				/*
				var intersection
				var ab = constraint.segment.getVectorAB()
				if (constraint.other instanceof Segment) {
					intersection = constraint.segment.intersection(constraint.other);
				} else if (constraint.other instanceof CustomPlane) {
					var sketchLineWC = sketch.plane.intersectionWithPlane(constraint.other)
					var dir = sketch.worldToSketchMatrix.transformVector(sketchLineWC.dir1)
					var p = sketch.worldToSketchMatrix.transformPoint(sketchLineWC.anchor)
					var abXcd = ab.cross(dir)
					var div = abXcd.lengthSquared()
					var anchorDiff = p.minus(constraint.segment.a.V3())
					var s = anchorDiff.cross(anchorDiff).dot(abXcd) / div
					var t = anchorDiff.cross(ab).dot(abXcd) / div

					//console.log(segmentIntersectsRay, a, b, "ab", ab, "p", p, "dir", dir, s > 0 && t / div >= 0 && t / div <= 1, "s", s, "t", t, "div", div)
					return {segmentAt: s, lineAt: t}
				}
				var abPos = constraint.segment.getS(intersection);
				var cdPos = constraint.other instanceof Segment
				? constraint.other.getS(intersection)
					 : t;
				var abLine = ab.unit().times(0.5 < abPos ? -16 : 16);
				var cdLine = (constraint.other instanceof Segment ?
					constraint.other.getVectorAB().unit()
						:dir).times(0.5 < cdPos ? -16 : 16);
				paintLineXY(intersection.plus(abLine).plus(cdLine), intersection.plus(abLine));
				paintLineXY(intersection.plus(abLine).plus(cdLine), intersection.plus(cdLine));
				*/
				break;
			case "colinear":
				// assume segments dont overlap
				var segments = constraint.segments;
				var ab = segments[0].getVectorAB().unit();
				var coord = ab.x > ab.y ? "x" : "y";
				if (ab[coord] < 0) {
					ab = ab.times(-1);
				}
				var scale = ab[coord];
				var offsetPoint = segments[0].a.V3().minus(ab.times(segments[0].a.V3()[coord] / scale));
				segments.sort((a, b) => a.a[coord] - b.a[coord]);
				for (var i = 0; i < segments.length - (constraint.fixed ? 0 : 1); i++) {
					var startS = max(segments[i].a.V3()[coord]/scale, segments[i].b.V3()[coord]/scale) + 8;
					var endS = startS + 20;
					paintLineXY(offsetPoint.plus(ab.times(startS)), offsetPoint.plus(ab.times(endS)), 0x000000, 2);

					startS = min(segments[(i + 1) % segments.length].a[coord]/scale, segments[(i + 1) % segments.length].b[coord]/scale) - 8;
					endS = startS - 20;
					paintLineXY(offsetPoint.plus(ab.times(startS)), offsetPoint.plus(ab.times(endS)), 0x000000, 2);
				}
				break;
			case "pointDistance":
			{
				var a = constraint.cs[0].V3(), b = constraint.cs[1].V3();
				var ab = b.minus(a);
				var ab1 = ab.unit();
				var abLength = ab.length();
				var ab90 = ab1.getPerpendicular();
				var texture = getTextureForString(constraint.distance);
				var textLength = texture.textWidthPx / TEXT_TEXTURE_HEIGHT * 20;
				paintLineXY(a.plus(ab90.times(6)), a.plus(ab90.times(22)));
				paintLineXY(b.plus(ab90.times(6)), b.plus(ab90.times(22)));

				paintLineXY(a.plus(ab90.times(14)), a.plus(ab90.times(14)).plus(ab1.times(abLength / 2 - textLength / 2 - 5)));
				paintLineXY(a.plus(ab90.times(14)).plus(ab1.times(abLength / 2 + textLength / 2 + 5)), a.plus(ab90.times(14)).plus(ab1.times(abLength)));
				var textCenter = a.plus(ab90.times(14)).plus(ab1.times(abLength / 2));
				gl.pushMatrix();
				gl.translate(textCenter);
				gl.scale(20, 20, 10);
				gl.multMatrix(M4.forSys(ab1, ab90.times(1), V3.Z));
				gl.translate(-textLength / 2 / 20, -0.5, 0);
				renderText(constraint.distance, 0xff0000);
				gl.popMatrix();
				break;
			}
			case "pointLineDistance":
				var texture = getTextureForString(constraint.distance);
				var p = constraint.point.V3(), ab = constraint.other.getVectorAB(), a = constraint.other.a.V3();
				var ab1 = ab.unit();
				var abLength = ab.length();
				var ab90 = ab1.getPerpendicular();
				var ap = p.minus(a);
				var apProj = ab90.times(ap.dot(ab90));
				paintLineXY(a, a.plus(apProj));
				paintLineXY(a.plus(apProj), p);
				var textLength = texture.textWidthPx / TEXT_TEXTURE_HEIGHT * 20;
				var textCenter =  a.plus(apProj.times(1/2))
				gl.pushMatrix();
				gl.translate(textCenter);
				gl.scale(20, 20, 10);
				gl.multMatrix(M4.forSys(apProj.unit(), ab1, V3.Z));
				gl.translate(-textLength / 2 / 20, -0.5, 0);
				renderText(constraint.distance, 0xff0000);
				gl.popMatrix();
				break;
			case "pointOnLine":
				var o = constraint.other;
				var ab = o instanceof Segment ? constraint.other.getVectorAB() : sketch.plane.normal.cross(o.normal);
				var ab1 = ab.unit();
				var p = constraint.point.V3();
				paintLineXY(p.plus(ab1.times(16)), p.plus(ab1.times(4)), 0x000000, 2);
				paintLineXY(p.plus(ab1.times(-16)), p.plus(ab1.times(-4)), 0x000000, 2);
				break;
			case "angle":
				var first = constraint.cs[0], second = constraint.cs[1];
				var intersection = first.intersection(second);
				var startAngle = (first.angleAB() + (constraint.f[0] == -1 ? PI : 0)) % (2 * PI), endAngle = startAngle + constraint.value;
				drawArc(intersection, startAngle, endAngle, 20, 21, 0x000000);
				break;
			case "equalLength":
				break;
			default:
				throw new Error("unknown constraint " + constraint.type);
			// console.log("errror", constraint.type);
		}

	});
}
function CustomPlane(anchor2, right, up, upStart, upEnd, rightStart, rightEnd, color, name) {
	var p = P3.forAnchorAndPlaneVectors(anchor2, right, up, CustomPlane.prototype);
	p.up = up;
	p.right = right;
	p.upStart = upStart;
	p.upEnd = upEnd;
	p.rightStart = rightStart;
	p.rightEnd = rightEnd;
	p.color = color;
	p.id = globalId++;
	p.name = name
	return p;
}
CustomPlane.prototype = Object.create(P3.prototype);
CustomPlane.prototype.constructor = CustomPlane;
Object.defineProperty(CustomPlane.prototype, "plane", { get: function () { return this } });
CustomPlane.prototype.toString = function() {
	return "Plane #" + this.id;
}
CustomPlane.prototype.distanceTo = function (line) {
	try {
		return [
			L3(this.anchor.plus(this.right.times(this.rightStart)), this.up),
			L3(this.anchor.plus(this.right.times(this.rightEnd)), this.up),
			L3(this.anchor.plus(this.up.times(this.upStart)), this.right),
			L3(this.anchor.plus(this.up.times(this.upEnd)), this.right)].map(function (line2, line2Index) {
			var info = line2.pointClosestTo2(line);
			if (isNaN(info.t) // parallel lines
				|| line2Index < 2 && this.upStart <= info.t && info.t <= this.upEnd
				|| line2Index >= 2 && this.rightStart <= info.t && info.t <= this.rightEnd ) {
				return info.distance;
			} else {
				return 1e9;
			}
		}, this).min();
	}catch (e) {
		console.log(e);
	}
}

var PlaneDefinition = function() {

}
PlaneDefinition.prototype.delete = function (line) {
	$('planeDefiner').set('display', 'none')
	features.remove(this)

	rebuildModel()
	updateFeatureDisplay()
}

//var sketchPlane = new CustomPlane(V3(0, 0,1), V3.X, V3.Y, -500, 500, -500, 500, 0xff00ff);
var currentSketch, features = [], currentPlaneDefinition;
function initModel() {
	features.push({type: "planeDefinition", planeType: "immediate", source: "new P3(V3.X, 450.00086883286616)", planeId: "sketchPlane", delete: PlaneDefinition.prototype.delete})
	//features.push({type: "planeDefinition", planeType: "immediate", source: "P3.normalOnAnchor(V3(1, 0, -1).unit(), V3.X)", planeId: "sketchPlane", delete: PlaneDefinition.prototype.delete})
	var sketch = new Sketch("sketchPlane");
	features.push(sketch);
	currentSketch = sketch;
	console.log("init sketch");
	sketch.elements.push(new Segment(0,0,100,0));
	sketch.elements.push(new Segment(100,0,100 * sqrt(3) /2,50));
	sketch.elements.push(new Segment(100 * sqrt(3) /2,50, 100, 100));
	sketch.elements.push(new Segment(100,100,0,0));
	makeCoincident(sketch.elements[0].b, sketch.elements[1].a, sketch);
	makeCoincident(sketch.elements[1].b, sketch.elements[2].a, sketch);
	makeCoincident(sketch.elements[2].b, sketch.elements[3].a, sketch);
	makeCoincident(sketch.elements[3].b, sketch.elements[0].a, sketch);
	selected = [sketch.elements[0].b, sketch.elements[1]];

	selected = sketch.elements.slice(0, 2)
	makeGroup("equalLength")
	selected = [sketch.elements[0], planes[1]]
	makePerpendicular()
	//makeGroup("parallel")
	//console.log("angle" ,rad2deg(sketch.elements[0].angleTo(sketch.elements[1])));
	//makeGroup("colinear");
	//sketch.constraints[0].value = deg2rad(30);
	//selected = sketch.elements.slice(1, 3);
	//makeDistance();
	//selected = [sketch.elements[0], sketch.elements[1]];
	//makeAngle();
	console.log("constraints", sketch.constraints);
	features.push(new Extrude("initExtrude", sketch.elements[0].name))
	features.push({type: "planeDefinition", planeType: "face", faceName: "initExtrudewall0", offset: 0, planeId: "planeCustom1", delete: PlaneDefinition.prototype.delete})
	//features.push(currentSketch = new Sketch("planeCustom1"))
	rebuildModel();

	/*
	 sketch.elements.push(new Segment(10,10,10,100));
	 sketch.elements.push(new Segment(10,10,100,10));
	 sketch.elements.push(new Segment(10,100,100,10));
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
	console.log("rebuilding model")
	//NLA.DEBUG = false
	//disableConsole()
	planes = [
		CustomPlane(V3.ZERO, V3.Y, V3.Z, -500, 500, -500, 500, 0xff0000, "planeYZ"),
		CustomPlane(V3.ZERO, V3.Z, V3.X, -500, 500, -500, 500, 0x00ff00, "planeZX"),
		CustomPlane(V3.ZERO, V3.X, V3.Y, -500, 500, -500, 500, 0x0000ff, "planeXY"),
	];
	modelCSG = undefined
	csgMesh = undefined
	features.forEach((feature) => {
		if (feature instanceof Sketch) {
			feature.plane = planes.find(p => p.name == feature.planeName);
			//console.log("LENGTHS", feature.plane.right.length(), feature.plane.up.length(), feature.plane.normal.length())
			feature.sketchToWorldMatrix =
				M4.forSys(feature.plane.right, feature.plane.up, feature.plane.normal, feature.plane.anchor)
			feature.worldToSketchMatrix = feature.sketchToWorldMatrix.inversed();
			recalculate(feature)
		} else if (feature.type && feature.type == "extrude") {
			var loopSketch = features.filter(f => f instanceof Sketch)
				.find(sketch => sketch.elements.some(el => el instanceof Segment && el.name == feature.segmentName))
			var loopSegment = loopSketch.elements.find(el => el instanceof Segment && el.name == feature.segmentName)
			var polygonPoints = loopSketch.getLoopForSegment(loopSegment).map(p => V3.create(p.x, p.y, 0))
			//console.log(polygonPoints.map(v =>v.ss))
			polygonPoints = loopSketch.sketchToWorldMatrix.transformedPoints(polygonPoints)
			var length = feature.end - feature.start
			var brep = BREP.extrude(polygonPoints, loopSketch.plane,
				loopSketch.plane.normal.times(-length), feature.name)

			if (modelCSG) {
				modelCSG = modelCSG[feature.operation](brep)
			} else {
				modelCSG = brep;
			}
			csgMesh = modelCSG.toNormalMesh()
		} else if (feature.type && "planeDefinition" == feature.type) {
			if ("face" == feature.planeType && feature.faceName) {
				var face = modelCSG.faces.find(face => face.name == feature.faceName)
				var plane = face.plane, right = plane.normal.getPerpendicular().normalized(), up = plane.normal.cross(right)

				var cp
				planes.push(cp = new CustomPlane(plane.anchor.plus(plane.normal.times(feature.offset)),
					right, up, -500, 500, -500, 500, 0xFF4500, feature.planeId))
			}
			if ("immediate" == feature.planeType) {
				var plane = eval(feature.source), right = plane.normal.getPerpendicular().normalized(), up = plane.normal.cross(right)

				planes.push(new CustomPlane(plane.anchor,
					right, up, -500, 500, -500, 500, 0xFF4500, feature.planeId))
			}
		} else {
			throw new Error("unknown feature");
		}
	})
	paintScreen()
	//enableConsole()
}
var circleShaderMaterial, circlePlaneGeometry, segmentPlaneGeometry, segmentMaterial;
var currentAddingSegment, connectingConstraint = [], currentPointConstraint, hoverHighlight = [];
// console.log = oldConsole;
var cvs, ctx;
var modeStack = [], mode = ""
var isShift;
var singleColorShader, textureColorShader, singleColorShaderHighlight, arcShader, arcShader2, lightingShader;
var extrudeFeature;
function rgbToVec4(color) {
	return [(color >> 16) / 255.0, ((color >> 8) & 0xff) / 255.0, (color & 0xff) / 255.0, 1.0];
}
function renderColor(mesh, color, mode) {
	singleColorShader.uniforms({
		color: rgbToVec4(color)
	}).draw(mesh);
}
function renderColorLines(mesh, color) {
	singleColorShader.uniforms({
		color: rgbToVec4(color)
	}).draw(mesh, gl.LINES);
}
var TEXT_TEXTURE_HEIGHT = 128;
function renderText(string, color) {
	var texture = getTextureForString(string);
	texture.bind(0);
	gl.pushMatrix();
	gl.scale(texture.width / TEXT_TEXTURE_HEIGHT, 1, 1);
	textureColorShader.uniforms({
		texture: 0,
		color: rgbToVec4(color)
	}).draw(textMesh);
	gl.popMatrix();
}
function drawVectors() {
	gl.pushMatrix();
	gl.scale(50, 50, 50)
	singleColorShader.uniforms({
		color: rgbToVec4(0xff0000)
	}).draw(vectorMesh);
	gl.popMatrix();

	gl.pushMatrix();
	gl.multMatrix(M4.forSys(V3.Y, V3.Z, V3.X))
	gl.scale(50, 50, 50)
	singleColorShader.uniforms({
		color: rgbToVec4(0x00ff00)
	}).draw(vectorMesh);
	gl.popMatrix();

	gl.pushMatrix();
	gl.multMatrix(M4.forSys(V3.Z, V3.X, V3.Y))
	gl.scale(50, 50, 50)
	singleColorShader.uniforms({
		color: rgbToVec4(0x0000ff)
	}).draw(vectorMesh);
	gl.popMatrix();

}
function paintScreen () {
	try {
		gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
		gl.loadIdentity();
		gl.scale(10, 10, 10);

		gl.loadIdentity();
		drawPlanes();
		/*
		 gl.translate(0, 0, -5);
		 gl.rotate(30, 1, 0, 0);
		 gl.rotate(angle, 0, 1, 0);
		 */
		//renderColor(ringMesh, 0xffff00);
		drawVectors()


		gl.pushMatrix();
		if (csgMesh) {
			var faceIndex = modelCSG.faces.length;
			while (faceIndex--) {
				var face = modelCSG.faces[faceIndex]
				var faceTriangleIndexes = csgMesh.faceIndexes.get(face)
				lightingShader.uniforms({
					color: rgbToVec4(hFace == face ? 0xff00ff : 0xff0000)
				}).draw(csgMesh, gl.TRIANGLES, faceTriangleIndexes.start, faceTriangleIndexes.count);
				/*
				singleColorShader.uniforms({
					color: rgbToVec4(0x0000ff)
				}).draw(csgMesh, gl.LINES);
				*/
			}
		}
		gl.popMatrix();

		features.filter(f => f instanceof Sketch).forEach(s => {
			gl.pushMatrix()
			paintSegments(s)
			gl.popMatrix()
		})
	} catch(e) {
		console.log(e);
	}
}
var gl;
var circleMesh, ringMesh, segmentMesh, textMesh, xyLinePlaneMesh, cubeMesh, indexBuffer, arcMesh, csgMesh, vectorMesh;
var modelCSG;
var lastSegment;
function setupCamera() {
	gl.matrixMode(gl.PROJECTION);
	gl.loadIdentity();
	//gl.perspective(70, gl.canvas.width / gl.canvas.height, 0.1, 1000);
	var lr = gl.canvas.width / 2 / zoomFactor;
	var bt = gl.canvas.height / 2 / zoomFactor;
	gl.ortho(-lr, lr, -bt, bt, -1e6, 1e6);
	gl.lookAt(
		eyePos.x, eyePos.y, eyePos.z,
		eyeFocus.x, eyeFocus.y, eyeFocus.z,
		eyeUp.x, eyeUp.y, eyeUp.z);
	gl.matrixMode(gl.MODELVIEW);
}
function drawPlanes() {
	planes.forEach(function (plane) {
		gl.pushMatrix();
		gl.multMatrix(M4.forSys(plane.right, plane.up, plane.normal));
		gl.translate(plane.rightStart, plane.upStart, plane.w);
		gl.scale(plane.rightEnd - plane.rightStart, plane.upEnd - plane.upStart, 1);

		var shader = hoverHighlight.contains(plane) ? singleColorShaderHighlight : singleColorShader;

		shader.uniforms({
			color: rgbToVec4(plane.color)
		}).draw(xyLinePlaneMesh, gl.LINES)

		gl.popMatrix()
	})
}

var faces = [], hFace = [];
CSG.Polygon.prototype.intersectsLine = function (line) {
	//console.log(line);
	var plane = this.plane
	var lambda = line.intersectWithPlaneLambda(plane)
	var p = line.at(lambda)
	//console.log("p", p)
	var direction = this.vertices[this.vertices.length - 1].pos.minus(this.vertices[0].pos)
	//var testRay = new Line3d(p, direction)
	var i = this.vertices.length - 1;
	var inside = false;
	while (i--) {
		var a = this.vertices[i].pos, b = this.vertices[i + 1].pos
		// check if segment ab intersects testRay
		if (segmentIntersectsRay(a, b, p, direction)) {
			inside = !inside
		}
	}
	return inside ? lambda : NaN;
}
function segmentIntersectsRay2(a, b, p, dir) {
	var ab = b.minus(a)
	var coord0 = "x", coord1 = "y"
	//console.log(SYL.isZero(ab.x), SYL.isZero(ab.y), SYL.isZero(ab.x) && SYL.isZero(ab.y))
	if (SYL.isZero(ab.x) && SYL.isZero(ab.y) || SYL.isZero(dir.x) && SYL.isZero(dir.y)) {
		coord0 = "z"
	}
	console.log(segmentIntersectsRay, coord0, coord1, dir[coord0] , ab[coord1] , dir[coord1] , ab[coord0])
	var div = dir[coord0] * ab[coord1] - dir[coord1] * ab[coord0];
	var t = (ab[coord0] * (p[coord1]- a[coord1]) + ab[coord1] * (a[coord0] - p[coord0])) / div
	var s = (dir[coord0] * (p[coord1]- a[coord1]) + dir[coord1] * (a[coord0] - p[coord0])) / div
	console.log(segmentIntersectsRay, a, b, "ab", ab, "p", p, "dir", dir, s > 0 && t / div >= 0 && t / div <= 1, "s", s, "t", t, "div", div, "coord0", coord0)
	return t > 0 && s >= 0 && s <= 1
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
CSG.prototype.toMesh = function() {
	function canonEdge(i0, i1) {
		var iMin = min(i0, i1), iMax = max(i0, i1)
		return (iMin << 16) | iMax
	}
	function uncanonEdge(key) {
		return [key >> 16, key & 0xffff]
	}
	faces = [];
	console.log("^calling tomesh");
	var mesh = new GL.Mesh({ normals: true, colors: true, lines:true });
	var indexer = new GL.Indexer();
	var sharedMap = new Map();
	var polygons = this.toPolygons();
	var i = polygons.length;
	while (i--) {
		var p = polygons[i];
		var sharedPolygons = sharedMap.get(p.shared);
		if (!sharedPolygons) {
			sharedMap.set(p.shared, sharedPolygons = []);
		}
		sharedPolygons.push(p);
	}
	console.log("Map", sharedMap);
	var allValidEdges = [];
	sharedMap.forEach((polygons, shared, map) => {
		var allSharedEdges = []
		faces.push({start: mesh.triangles.length * 3, count: polygons.length * 3, shared: shared})
		polygons.forEach((polygon) => {
			var indices = polygon.vertices.map(function(vertex) {
				vertex.color = polygon.shared || [1, 1, 1];
				return indexer.add(vertex);
			});
			for (var i = 2; i < indices.length; i++) {
				mesh.triangles.push([indices[0], indices[i - 1], indices[i]]);
			}
			for (var i = 0; i < indices.length; i++) {
				var i0 = indices[i], i1 = indices[(i + 1) % indices.length];
				allSharedEdges.push(canonEdge(i0, i1))
			}
		})
		allSharedEdges.sort()
		console.log("allSharedEdges", allSharedEdges)
		var validEdges = []
		var i
		for (i = allSharedEdges.length - 2; i >= 0; i--) {
			if (allSharedEdges[i] == allSharedEdges[i + 1]) {
				// remove both
				i--;
				continue;
			}
			validEdges.push(allSharedEdges[i + 1])
		}
		if (i == -1) {
			validEdges.push(allSharedEdges[0])
		}
		console.log(validEdges)
		allValidEdges = allValidEdges.concat(validEdges	)
	})
	allValidEdges.sort()
	console.log("allValidEdges", allValidEdges)
	mesh.lines = allValidEdges.unique().map(uncanonEdge)
	console.log("indexer.unique.length", indexer.map);
	console.log(mesh.lines)
	mesh.vertices = indexer.unique.map(function(v) { return [v.pos.x, v.pos.y, v.pos.z]; });
	//mesh.normals = indexer.unique.map(function(v) { return [v.normal.x, v.normal.y, v.normal.z]; });
	mesh.colors = indexer.unique.map(function(v) { return v.color; });
	mesh.compile()
	return mesh;
}
function template(templateName, map) {
	var html = $(templateName).text;
	for (var key in map) {
		html = html.replace(new RegExp('\\$'+key, 'g'), map[key]);
	}
	return $(new Element('div', {html: html.trim()}).firstChild);
}
function updateFeatureDisplay() {
	var div = $('featureDisplay')
	div.erase('text')
	features.forEach(function (feature) {
		if (feature.type == "extrude") {
			var newChild = template('templateFeatureExtrude', {what: "EXTRUDE", name: feature.name})
			newChild.getElement('[name=edit]').onclick = function () {
				editExtrude(feature)
			}
			newChild.getElement('[name=delete]').onclick = function () {
				feature.delete()
			}
			newChild.inject(div);
		} else if (feature.type == "planeDefinition") {
			var newChild = template('templateFeatureExtrude', {what: "PLANE", name: feature.planeId})
			newChild.getElement('[name=edit]').onclick = function () {
				editPlane(feature)
			}
			newChild.getElement('[name=delete]').onclick = function () {
				feature.delete()
			}
			newChild.inject(div);
		} else if (feature instanceof Sketch) {
			var newChild = template('templateFeatureExtrude', {what: "SKETCH", name: feature.name})
			newChild.getElement('[name=edit]').onclick = function () {
				editSketch(feature)
			}
			newChild.getElement('[name=delete]').onclick = function () {
				feature.delete()
			}
			newChild.inject(div);
		}
	})
}
function updateSelected() {
	var div = $('selectedElements')
	div.erase('text')
	selected.forEach(function (sel) {
		var newChild = template("template", {what: sel.constructor.name, name: sel.name});
		newChild.onmouseover = function (e) {
			highlighted = [sel];
			paintScreen();
		};
		newChild.onmouseout = function (e) {
			highlighted = [];
			paintScreen();
		};
		newChild.getElement(".unsel").onclick = function (e) {
			(sel instanceof SegmentEndPoint ? sel.line : sel).remove();
			updateSelected();
			paintScreen();
		};
		newChild.inject(div);
	});
	div = $('selectedConstraints')
	div.erase('text')
	selected.map((el) => currentSketch.getConstraintsFor(el)).concatenated().unique().forEach(function (constraint) {
		if ("pointDistance" == constraint.type || "pointLineDistance" == constraint.type) {
			console.log("constraint.type", constraint.type);
			var newChild = template("templateDistance", {name: constraint.type, id: constraint.id});
			newChild.inject(div);
			newChild.getElement(".distanceInput").value = constraint.distance;
			newChild.getElement(".distanceInput").onchange = function (e) {
				constraint.distance = e.target.value;
				rebuildModel();
				paintScreen();
			}
		} else if ("angle" == constraint.type) {
			var newChild = template("templateAngle", {name: constraint.type, id: constraint.id});
			newChild.inject(div);
			console.log("AAANGLE");
			var input = newChild.getElement(".distanceInput");
			newChild.getElement(".distanceInput").value = round(rad2deg(constraint.value, 5));
			newChild.getElement(".fa").onclick = () => { constraint.f[0] *= -1; input.value = round(rad2deg(abs(constraint.value = (-PI + constraint.value) % (2 * PI))), 3); paintScreen() };
			newChild.getElement(".fb").onclick = () => { constraint.f[1] *= -1; input.value = round(rad2deg(abs(constraint.value = (-PI + constraint.value) % (2 * PI))), 3); paintScreen() };
			newChild.getElement(".fv").onclick = () => { input.value = round(rad2deg(abs(constraint.value -= sign(constraint.value) * 2*PI)), 3); paintScreen() };
			input.onchange = function (e) {
				constraint.value = deg2rad(e.target.value);
				rebuildModel();
				paintScreen();
			}
		} else {
			console.log("wat");
			var newChild = template("templateConstraint", {name: constraint.type});
			newChild.inject(div);

			constraint.cs.forEach(function (segment) {
				var subChild = template("templateConstraintSub", {name: segment.toString()});
				subChild.inject(newChild);
				subChild.getElement(".removeFromConstraint").onclick = function (e) {
					removeFromGroupConstraint(segment, currentSketch, constraint.type);
					updateSelected();
					paintScreen();
				};
				subChild.onmouseover = function (e) {
					highlighted = [segment];
					paintScreen();
				};
				subChild.onmouseout = function (e) {
					highlighted = [];
					paintScreen();
				};
			});
		}
		newChild.getElement(".remove").onclick = function (e) {
			deleteConstraint(currentSketch, constraint)
			updateSelected();
		}
	});
}
window.onerror = function (errorMsg, url, lineNumber, column, errorObj) {
	console.log(errorObj);
}
function ensureNumbers() {
	for (var i = 0; i < arguments.length; i++) {
		if (typeof arguments[i] !== 'number') {
			// Arrays.prototype.slice.call is inefficient, but it doesn't matter here
			throw new Error("ensureNumbers arguments[" + (i) + "] is not a number. " + typeof arguments[i] + " == typeof " + arguments[i]);
		}
	}
}
var main = function () {
	gl = GL.create({});
	gl.fullscreen();
	gl.canvas.oncontextmenu = () => false;

	gl.scaleVector = function(x, y, z) {
		gl.multMatrix(M4.scaleVector(x, y, z));
	};

	setupCamera();
	//gl.cullFace(gl.FRONT_AND_BACK);
	gl.clearColor(1.0, 1.0, 1.0, 0.0);
	gl.enable(gl.BLEND);
	gl.enable(gl.DEPTH_TEST);
	gl.enable(gl.CULL_FACE);
	gl.depthFunc(gl.LEQUAL)
	gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA); // TODO ?!

	cubeMesh = GL.Mesh.cube();
	ringMesh = createRingMesh(8, 10, 16);
	indexBuffer = {"index": createIndexBuffer(64)};
	circleMesh = createCircleMesh(6, 16);
	arcMesh = createArcMesh(0, 1.0, 19);
	segmentMesh = createRectangleMesh(0, -1, 1, 1);
	textMesh = createRectangleMesh(0, 0, 1, 1);
	vectorMesh = rotationMesh([V3.ZERO, V3(0, 0.05, 0), V3(0.8, 0.05), V3(0.8, 0.1), V3(1, 0)], L3.X, Math.PI * 2, 8, false)
	xyLinePlaneMesh = new GL.Mesh({lines: true, triangles: false});
	xyLinePlaneMesh.vertices = [[0, 0], [0, 1], [1, 1], [1, 0]];
	xyLinePlaneMesh.lines = [[0, 1], [1, 2], [2, 3], [3, 0]];
	xyLinePlaneMesh.compile();

	singleColorShader = new GL.Shader(vertexShaderBasic, fragmentShaderColor);
	singleColorShaderHighlight = new GL.Shader(vertexShaderBasic, fragmentShaderColorHighlight);
	textureColorShader = new GL.Shader(vertexShaderTextureColor, fragmentShaderTextureColor);
	arcShader = new GL.Shader(vertexShaderRing, fragmentShaderColor);
	arcShader2 = new GL.Shader(vertexShaderArc, fragmentShaderColor);
	lightingShader = new GL.Shader(vertexShaderLighting, fragmentShaderLighting);
	//	console.log(mesh.vertices);
	initModel()
	updateFeatureDisplay()

	var a = CSG.cube();
	var b = CSG.sphere({ radius: 1.35 });
	a.setColor(1, 1, 0);
	b.setColor(0, 0.5, 1);
	//testCSG = createCSG();
	paintScreen();
	window.onkeypress = function (e) {
		if ("Delete" == e.key) {
			selected.map((x) => x instanceof SegmentEndPoint ? x.line : x).unique().forEach((x) => x.remove());
			paintScreen();
		}
	}
	gl.onmousedown = function (e) {

		if (!(e.buttons & 1)) {
			return;
		}
		var mouseLine = getMouseLine(e);
		var sketch = currentSketch
		console.log("mouseLine", getMouseLine(e).toString(), "mode", mode);
		if (mode == "addsegments") {
			var intersection = mouseLine.intersectWithPlane(sketch.plane)
			console.log("INTERSECTION", intersection.ss)
			if (intersection == null) {
				return;
			}
			var sketchCoords = sketch.worldToSketchMatrix.transformPoint(intersection)
			currentPointConstraint = undefined;
			lastSegment = currentAddingSegment;
			currentAddingSegment = new Segment(sketchCoords.x, sketchCoords.y, sketchCoords.x, sketchCoords.y);
			sketch.elements.push(currentAddingSegment);
			if (lastSegment) {
				makeCoincident(lastSegment.b, currentAddingSegment.a, sketch);
			}
			//console.log("e", e);
			console.log("lines", sketch.elements);
			//console.log("constraints", constraints);
			console.log(currentAddingSegment.a);
			rebuildModel(sketch);
			console.log(currentAddingSegment.a);
		} else if (mode == "face-select") {
			var face = hFace
			if (face instanceof BREP.Face) {
				modeClickCallback(face)
			}
		} else if (mode == "plane-select") {
			var plane = hoverHighlight[0]
			if (plane instanceof CustomPlane) {
				modeClickCallback(plane)
			}
		} else if (mode == "segment-select") {
			var segment = hoverHighlight[0]
			if (segment instanceof Segment) {
				modeClickCallback(segment)
			}
		} else {
			if (!GL.keys.SHIFT) {
				selected = hoverHighlight.slice();
			} else {
				if (selected.contains(hoverHighlight[0])) {
					selected.remove(hoverHighlight[0]);
				} else {
					selected.push(hoverHighlight[0]);
				}
			}
			updateSelected();
			console.log ("selected is now ", selected);
			console.log ("highlighted is now ", highlighted);
			paintScreen();
		}
	}
	gl.onmousemove = function (e) {
		// try {
		var mouseLine = getMouseLine(e);
		var sketch = currentSketch
		var intersection = sketch.plane && mouseLine.intersectWithPlane(sketch.plane)
		if (mode == "addsegments") {
			if (!intersection) {
				return;
			}
			// parallel to mouseline
			var sketchCoords = sketch.worldToSketchMatrix.transformPoint(intersection)
			currentAddingSegment && removeCoincidence(currentAddingSegment.b, sketch);
			//console.log(mousePos);
			var points = getAllPoints(sketch);
			currentAddingSegment && points.removeAll(currentAddingSegment.points);
			var pointDistances = points
				.map(function (point) { return {point: point, distance: point.distanceToCoords(sketchCoords)}; })
				.filter(function (pair) { return pair.distance < 16; })
				.sort(function (a, b) { return a.distance - b.distance; });
			if (pointDistances.length != 0) {
				makeCoincident(pointDistances[0].point, currentAddingSegment.b, sketch);
				sketchCoords = pointDistances[0].point;
			}
			if (currentAddingSegment) {
				currentAddingSegment.b.x = sketchCoords.x;
				currentAddingSegment.b.y = sketchCoords.y;
			}
			paintScreen();
			/*
			 for (var i in lines) {
			 var line = lines[i];
			 //console.log("line", line);
			 if (line != currentAddingSegment && line.distanceTo(mousePos.x, mousePos.y) < 16) {
			 mousePos = line.getClosestPoint(mousePos.x, mousePos.y);
			 break;
			 }
			 };
			 */
		} else {
			var nearest = Infinity
			hFace = null
			if (modelCSG) {
				modelCSG.faces.forEach((face) => {
					var distance;
					if (!isNaN(distance = face.intersectsLine(mouseLine))) {
						if (distance < nearest) {
							//console.log(hFace, distance)
							nearest = distance;
							hFace = face
							//console.log(polygon.shared, polygon, "[" + polygon.vertices.map((v)=>v2as(v.pos)).join(", ") +"]",v2as(mouseLine.anchor.csg()), v2as(mouseLine.direction.csg()), distance);
						}
					}
				});
			}
			var planesWithDistance = planes.map((plane) => ({el: plane, distance: plane.distanceTo(mouseLine)}) );
			hoverHighlight = [];
			if (intersection != null) {
				var sketchCoords = sketch.worldToSketchMatrix.transformPoint(intersection)
				planesWithDistance = sketch.elements.concat(getAllPoints(sketch))
					.map(function (el) { return {"el": el, "distance": el.distanceToCoords(sketchCoords)}; })
					.concat(planesWithDistance)
			}
			var elsByDistance = planesWithDistance
				.filter(function (pair) { return pair.distance < 16; })
				.sort(function (a, b) { return chainComparison(pointsFirst(a.el, b.el), a.distance - b.distance); });
			//console.log("elsByDistance", elsByDistance);
			if (0 != elsByDistance.length) {
				hoverHighlight.push(elsByDistance[0].el);
			}
			paintScreen();
		}
		if (e.dragging) {
			if (e.buttons & 4) {
				// pan
				var moveCamera = V3(-e.deltaX * 2 / gl.canvas.width, e.deltaY * 2 / gl.canvas.height, 0);
				var inverseProjectionMatrix = gl.projectionMatrix.inversed();
				var worldMoveCamera = inverseProjectionMatrix.transformVector(moveCamera);
				eyePos = eyePos.plus(worldMoveCamera);
				eyeFocus = eyeFocus.plus(worldMoveCamera);
				setupCamera();
				paintScreen();
			}
			if (e.buttons & 2) {
				var rotateLR = deg2rad(-e.deltaX / 6.0);
				var rotateUD = deg2rad(-e.deltaY / 6.0);
				// rotate
				var matrix = M4.rotation(rotateLR, eyeUp);
				//var horizontalRotationAxis = eyeFocus.minus(eyePos).cross(eyeUp);
				var horizontalRotationAxis = eyeUp.cross(eyePos.minus(eyeFocus));
				matrix = matrix.times(M4.rotation(rotateUD, horizontalRotationAxis));
				//console.log(matrix);
				eyePos = matrix.transformPoint(eyePos);
				eyeUp = matrix.transformVector(eyeUp)
				//console.log("eyepos", eyePos);
				setupCamera();
				paintScreen();
			}
		}
	}
	$(gl.canvas).addEvent('mousewheel', function (e) {
		//console.log(e);
		zoomFactor *= pow(0.9, -e.wheel);
		var mouseCoords = e.client;
		var moveCamera = V3(mouseCoords.x * 2 / gl.canvas.width - 1, -mouseCoords.y * 2 / gl.canvas.height + 1, 0).times(1 - 1 / pow(0.9, -e.wheel));
		var inverseProjectionMatrix = gl.projectionMatrix.inversed();
		var worldMoveCamera = inverseProjectionMatrix.transformVector(moveCamera);
		//console.log("moveCamera", moveCamera);
		//console.log("worldMoveCamera", worldMoveCamera);
		eyePos = eyePos.plus(worldMoveCamera);
		eyeFocus = eyeFocus.plus(worldMoveCamera);
		setupCamera();
		paintScreen();
	});
	$("clearmode").onclick = function () {
		if (connectingConstraint) {
			removeCoincidence(currentAddingSegment.a, currentSketch);
		}
		if (currentAddingSegment.b.coincidence) {
			removeCoincidence(currentAddingSegment.b.coincidence, currentSketch);
		}
		currentSketch.elements.remove(currentAddingSegment);
		currentAddingSegment = undefined;
		console.log("mode cleared");
		paintScreen();
		mode = "";
	}
}
var eyePos = V3(0, 0, 1000), eyeFocus = V3.ZERO, eyeUp = V3.Y;
function getMouseLine(pos) {
	var ndc1 = V3(pos.x * 2 / gl.canvas.width - 1, -pos.y * 2 / gl.canvas.height + 1, 0);
	var ndc2 = V3(pos.x * 2 / gl.canvas.width - 1, -pos.y * 2 / gl.canvas.height + 1, 1);
	//console.log(ndc);
	var inverseProjectionMatrix = gl.projectionMatrix.inversed();
	var s = inverseProjectionMatrix.transformPoint(ndc1);
	var dir = inverseProjectionMatrix.transformPoint(ndc2).minus(s);
	return L3.anchorDirection(s, dir);
}
function rotationMesh(vertices, axis, total, count, close) {
	var mesh = new GL.Mesh({normals: false})
	var vc = vertices.length, vTotal = vc * count

	for (var i = 0; i < count; i++) {
		var angle = total / count * i
		var m = M4()
		M4.rotationLine(axis.anchor, axis.dir1, angle, m)
		Array.prototype.push.apply(mesh.vertices, m.transformedPoints(vertices).map(v => v.els()))

		// add triangles
		for (var j = 0; j < vc - 1; j++) {
			mesh.triangles.push([i * vc + j + 1, i * vc + j, (i + 1) * vc + j].map(x => x % vTotal))
			mesh.triangles.push([i * vc + j + 1, (i + 1) * vc + j, (i + 1) * vc + j + 1].map(x => x % vTotal))
		}
	}

	mesh.compile()
	return mesh
}
var modeClickCallback
function setupSelectors(el, feature){
	el.getElements('.face-select')
		.removeEvents()
		.addEvent('click', function (uh) {
			var prevMode = mode, prevModeClickCallback = modeClickCallback
			var el = this
			this.addClass('selecting')
			el.set('text', 'Click on a face')
			mode = 'face-select'
			modeClickCallback = function (face) {
				mode = prevMode
				modeClickCallback = prevModeClickCallback
				el.removeClass('selecting')
				el.set('text', face.name)

				feature[el.getProperty('name')] = face.name
				rebuildModel()
			}
		})
		.each(el => el.innerHTML = feature[$(el).getProperty('name')])
		.removeClass('selecting')

	el.getElements('.plane-select')
		.removeEvents()
		.addEvent('click', function (e) {
			var prevMode = mode, prevModeClickCallback = modeClickCallback
			var el = this
			this.addClass('selecting')
			el.set('text', 'Click on a plane')
			mode = 'plane-select'
			modeClickCallback = function (plane) {
				mode = prevMode
				modeClickCallback = prevModeClickCallback
				el.removeClass('selecting')
				el.set('text', plane.name)

				feature[el.getProperty('name')] = plane.name
				rebuildModel()
			}
		})
		.each(el => el.innerHTML = feature[$(el).getProperty('name')])
		.removeClass('selecting')

	el.getElements('.segment-select')
		.removeEvents()
		.addEvent('click', function (e) {
			var prevMode = mode, prevModeClickCallback = modeClickCallback
			var el = this
			this.addClass('selecting')
			el.set('text', 'Click on a sketch segment')
			mode = 'segment-select'
			modeClickCallback = function (segment) {
				mode = prevMode
				modeClickCallback = prevModeClickCallback
				el.removeClass('selecting')
				el.set('text', segment.name)

				feature[el.getProperty('name')] = segment.name
				rebuildModel()
			}
		})
		.each(function (el) {
			el.set('text', feature[$(el).getProperty('name')])
		})
		.removeClass('selecting')

	el.getElements('.string-id-input, .select-text')
		.removeEvents()
		.addEvent('change', function (e) {
			// TODO: check if unique
			var propName = this.getProperty('name')
			feature[propName] = this.value
			rebuildModel()
		})
		.each(el => el.value = feature[$(el).getProperty('name')])

	el.getElements('.select-select')
		.removeEvents()
		.addEvent('change', function (e) {
			var propName = this.getProperty('name')
			feature[propName] = this.value
			rebuildModel()
		})
		.each(el => el.value = feature[$(el).getProperty('name')])

	el.getElements('.dimension-input')
		.removeEvents()
		.addEvent('mousewheel', function (e) {
			this.set('value', Math.round10( parseFloat(this.value) + Math.sign(e.wheel), -6))
			this.fireEvent('change')
		})
		.addEvent('click', function (e) {
			if (e.event.button == 1) {
				this.set('value', this.getProperty('defaultvalue'))
				this.fireEvent('change')
			}
		})
		.addEvent('change', function (e) {
			var propName = this.getProperty('name')
			feature[propName] = this.value = NLA.forceFinite(this.value)
			rebuildModel()
		})
		.each(el => el.value = feature[$(el).getProperty('name')])

	el.getElement('[name=done]')
		.removeEvents()
		.addEvent('click', function () {
			el.setStyle('display', 'none')
		})
	el.getElement('[name=delete]')
		.removeEvents()
		.addEvent('click', function () {
			feature.delete()
		})

}
function chainComparison(diff1, diff2) {
	return diff1 != 0 ? diff1 : diff2;
}
function reverse(sketch) {
	sketch.varMap.forEach(function (pointIndex, point) {
		point.x = sketch.x[pointIndex];
		point.y = sketch.x[pointIndex + 1 ];
	});
}
function pointsFirst(a, b) {
	//console.log((a instanceof SegmentEndPoint ? -1 : 1) - (b instanceof SegmentEndPoint ? -1 : 1));
	return (a instanceof SegmentEndPoint ? -1 : 1) - (b instanceof SegmentEndPoint ? -1 : 1);
}
function clicky() {
	for (var j = 0; j < 1; j++) {
		currentSketch.gaussNewtonStep();
	}
	reverse();
	paintScreen();
}
function deleteConstraint(sketch, constraint) {
	sketch.constraints.remove(constraint)

}
function getGroupConstraint(el, sketch, type) {
	return sketch.constraints.find(function (c) { return c.type == type && (c.fixed == el || c.segments.contains(el)); });
}
// removes segment or plane from group constraint
function removeFromGroupConstraint(el, sketch, type) {
	var groupConstraint = getGroupConstraint(el, sketch, type);
	if (!groupConstraint) return;
	groupConstraint.segments.remove(el);
	if (1 == groupConstraint.segments.length) {
		sketch.constraints.remove(groupConstraint);
	}
	if (groupConstraint.fixed == el) {
		groupConstraint.fixed = undefined;
	}
}
// colinear, equal length or parallel
/*
function makeGroup2(type) {
	var els = selected.filter((el) => el instanceof Segment || ('equalLength' != type && el instanceof CustomPlane) );
	console.log("makeGroup, els", els, type);
	if (els.length >= 2) {
		var newConstraint = {"fixed": null, "type": type, "which": []};
		var oldConstraints = [];
		for (var i = 0; i < els.length; i++) {
			var el = els[i];
			var segConstraint = getGroupConstraint(el, currentSketch, type);
			if (segConstraint) {
				oldConstraints.push(segConstraint);
				if (newConstraint.fixed && segConstraint.fixed) {
					throw new Error("cannot have two fixed");
				}
				if (segConstraint.fixed) {
					newConstraint.fixed = segConstraint.fixed;
				}
				if (segConstraint != groupConstraint) {
					segConstraint.which.forEach(function (segConstraintOther) {
						groupConstraint.which.push(segConstraintOther);
					});
					currentSketch.constraints.remove(segConstraint);
				}
			} else {
				if (el instanceof CustomPlane) {
					if (null != groupConstraint.fixed) {
						throw new Error("cannot have two fixed");
					}
					newConstraint.fixed = el;
				} else {
					newConstraint.which.push(el);
				}
			}
		}
		constraints.which.push(el);
		constraints.which.removeAll(oldConstraints);
		rebuildModel();
		paintScreen();
	}
	console.log("sketch.constraints", currentSketch.constraints);
}
*/
function makeGroup(type) {
	var els = selected.filter((el) => el instanceof Segment || ('equalLength' != type && el instanceof CustomPlane) );
	if (els.length < 2) return;

	var newGroup = els.map(el => {
		var c = getGroupConstraint(el, currentSketch, type);
		return c
			? (c.fixed ? c.segments.concat(c.fixed) : c.segments)
			: el; }).concatenated().unique();
	var fixeds = newGroup.filter(el => el instanceof CustomPlane);
	if (1 < fixeds.length) {
		throw new Error("cannot have two fixed");
	}

	var oldConstraints = els.map(el => getGroupConstraint(el, currentSketch, type)).filter(x => x);
	currentSketch.constraints.removeAll(oldConstraints);

	// add new constraint
	// fixeds[0] may be null
	var segments = newGroup.filter(x => x instanceof Segment);
	var f = fixeds[0];
	currentSketch.constraints.push(new Constraint(type, f ? segments.concat(f) : segments, {fixed: fixeds[0], segments: segments}));

	rebuildModel(currentSketch);
	paintScreen();
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
	}
	//constrains: o => this.cs.contains(o)
}

function makeAngle() {
	var selSegments = selected.filter(function (s) { return s instanceof Segment || s instanceof CustomPlane; });
	console.log(selSegments);
	if (selSegments.length == 2) {
		var segment = selSegments[0] instanceof Segment ? selSegments[0] : selSegments[1];
		var other = selSegments[0] instanceof Segment ? selSegments[1] : selSegments[0];
		if (! (segment instanceof Segment)) {
			throw new Error ("at least one must be a segment");
		}
		currentSketch.constraints.push(new Constraint("angle", [segment, other], {segment:segment, other:other, f: [1, 1], value: selSegments[0].angleTo(selSegments[1])})); //
		rebuildModel();
		paintScreen();
		updateSelected();
	}
}
function makePerpendicular() {
	var selSegments = selected.filter(function (s) { return s instanceof Segment || s instanceof CustomPlane; });
	if (selSegments.length == 2) {
		var segment = selSegments[0] instanceof Segment ? selSegments[0] : selSegments[1];
		var other = selSegments[0] instanceof Segment ? selSegments[1] : selSegments[0];
		if (! (segment instanceof Segment)) {
			throw new Error ("at least one must be a segment");
		}
		currentSketch.constraints.push(new Constraint("perpendicular", [segment, other], {segment:segment, other:other}));
		rebuildModel();
		paintScreen();
	}
}
function makeDistance() {
	if (2 != selected.length) return;
	var point = selected[0] instanceof SegmentEndPoint ? selected[0] : selected[1];
	var other = selected[0] instanceof SegmentEndPoint ? selected[1] : selected[0];
	if (point instanceof SegmentEndPoint && (other instanceof Segment || other instanceof SegmentEndPoint || other instanceof CustomPlane)) {
		var newConstraint;
		if (other instanceof SegmentEndPoint) {
			newConstraint = new Constraint("pointDistance", selected.slice(), {distance:round(other.distanceToCoords(point))});
		} else {
			newConstraint = new Constraint("pointLineDistance", [point, other], {point:point, other:other, distance:round(other.distanceToCoords(point))});
			console.log(newConstraint);
		}
		currentSketch.constraints.push(newConstraint);
		rebuildModel();
		paintScreen();
		updateSelected();
		$("distanceInput"+newConstraint.id).select();
	}
}
function selPointOnLine() {
	if (2 != selected.length) return;
	var point = selected[0] instanceof SegmentEndPoint ? selected[0] : selected[1];
	var other = selected[0] instanceof SegmentEndPoint ? selected[1] : selected[0];
	if (point instanceof SegmentEndPoint && (other instanceof Segment || other instanceof CustomPlane)) {
		var newConstraint = new Constraint("pointOnLine", [point, other], {point:point, other:other});
		console.log(newConstraint);
		currentSketch.constraints.push(newConstraint);
		rebuildModel();
		paintScreen();
		updateSelected();
	}
	$("distanceInput"+newConstraint.id).select();
}
function makeSelCoincident() {
	var selPoints = selected.filter(function (s) { return s instanceof SegmentEndPoint; });
	if (selPoints.length >= 2) {
		for (var i = 1; i < selPoints.length; i++) {
			makeCoincident(selPoints[0], selPoints[i], currentSketch);
		}
		rebuildModel();
		paintScreen();
	}
}
function createRingMesh(insideRadius, outsideRadius, cornerCount) {
	var mesh = new GL.Mesh({});
	for (var i = 0; i < cornerCount; i++) {
		mesh.vertices.push(
			[outsideRadius * cos(2 * PI * i / cornerCount), outsideRadius * sin(2 * PI * i / cornerCount), 0],
			[insideRadius * cos(2 * PI * (i + 0.5) / cornerCount), insideRadius * sin(2 * PI * (i + 0.5) / cornerCount), 0]
		);
		mesh.triangles.push(
			[0, 2, 1].map(function (offset) { return (2 * i + offset) % (2 * cornerCount); }),
			[1, 2, 3].map(function (offset) { return (2 * i + offset) % (2 * cornerCount); })
		);
	}
	mesh.compile();
	return mesh;
}
function createArcMesh(insideRadius, outsideRadius, cornerCount) {
	var mesh = new GL.Mesh({});
	for (var i = 0; i < cornerCount; i++) {
		mesh.vertices.push(
			[outsideRadius, i / (cornerCount - 1), 0],
			[insideRadius, i / (cornerCount - 1), 0]
		);
		if (i == cornerCount - 1) break;
		mesh.triangles.push(
			[0, 2, 1].map(function (offset) { return (2 * i + offset) % (2 * cornerCount); }),
			[1, 2, 3].map(function (offset) { return (2 * i + offset) % (2 * cornerCount); })
		);
	}
	mesh.compile();
	return mesh;
}
function createCircleMesh(radius, cornerCount) {
	var mesh = new GL.Mesh({});
	mesh.vertices.push([0, 0, 0]);
	for (var i = 0; i < cornerCount; i++) {
		mesh.vertices.push(
			[radius * cos(2 * PI * i / cornerCount), radius * sin(2 * PI * i / cornerCount), 0]
		);
		mesh.triangles.push(
			[0, 1 + i, 1 + (i + 1) % cornerCount]
		);
	}
	mesh.compile();
	return mesh;
}
function createRectangleMesh(x0, y0, x1, y1) {
	var mesh = new GL.Mesh({});
	mesh.vertices.push(
		[x0, y0],
		[x1, y0],
		[x0, y1],
		[x1, y1]
	);
	mesh.triangles.push(
		[0, 1, 2],
		[3, 2, 1]
	);
	mesh.compile();
	return mesh;
}
function createIndexBuffer(count) {
	var buffer = new GL.Buffer(gl.ARRAY_BUFFER, Float32Array);
	for (var i = 0; i < count; i++) {
		buffer.data.push([i]);
	}
	buffer.compile();
	//console.log("buffer.buffer.length",buffer.buffer.length, buffer.buffer.spacing);
	return buffer;
}

function faceSketchPlane() {
	eyePos = eyeFocus.plus(currentSketch.plane.normal.times(100));
	eyeUp = currentSketch.plane.up;
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
Extrude.prototype.delete = function () {
	$('extrudeEditor').set('display', 'none')
	features.remove(this)

	rebuildModel()
	updateFeatureDisplay()
}
function addSketch() {
	var sketch = new Sketch();
	features.push(sketch);
	editSketch(sketch)
	updateFeatureDisplay()
}
function addPlane() {
	currentPlaneDefinition = {type: "planeDefinition", planeType: "face",
		faceName: null, offset: 0, planeId: "plane"+(globalId++), delete: PlaneDefinition.prototype.delete}
	features.push(currentPlaneDefinition)
	editPlane(currentPlaneDefinition)
	updateFeatureDisplay()
}
function addExtrude() {
	var feature = new Extrude();
	features.push(feature);
	editExtrude(feature)
	updateFeatureDisplay()
}
function editSketch(feature) {
	currentSketch = feature;
	var div = $("sketchEditor")
	div.setStyle("display", "block")
	setupSelectors(div, feature)

	if (!feature.planeName) {
		div.getElement("[name=planeName]").fireEvent("click")
	}
}
function editPlane(feature) {
	var div = $("planeDefiner")
	div.setStyle("display", "block")
	setupSelectors(div, feature)
	if (feature.planeType == "face" && !feature.faceName) {
		div.getElement("[name=faceName]").fireEvent("click")
	}
}
function editExtrude(feature) {
	var div = $("extrudeEditor")
	div.setStyle("display", "block")
	setupSelectors(div, feature)
	/*
	div.getElement("[name=operation]")
		.set("value", feature.operation)
		.addEvent("change", function () {
			feature.operation = this.value
			rebuildModel()
		})
		*/
	if (!feature.segmentName) {
		div.getElement("[name=segmentName]").fireEvent("click")
	}
}



















