var S0 = 0, S1 = 1, D0 = 2, D1 = 3;
Array.prototype.remove = function (o) {
    var index = this.indexOf(o);
    if (index != -1) {
        this.splice(index, 1);
    }
};
Array.prototype.contains = function (o) {
    return -1 != this.indexOf(o);
};
Array.prototype.isEmpty = function (o) {
    return 0 == this.length;
};
Array.prototype.removeAll = function (os) {
    for (var i = this.length; --i > 0;) {
        if (os.contains(this[i])) {
            this.splice(i, 1);
        }
    }
};
var lineId = 0;
function Point(x, y) {
	this.x = x;
	this.y = y;
}
function SegmentEndPoint(x, y, line) {
	this.x = x;
	this.y = y;
	this.line = line;
	this.distanceTo = function (p) {
		return distance(this.x, this.y, p.x, p.y);
	}
	this.$V = function () {
		return $V([this.x, this.y]);
	};
}
function Segment(x1, y1, x2, y2) {
	this.points = [new SegmentEndPoint(x1, y1, this), new SegmentEndPoint(x2, y2, this)];
	this.a = this.points[0];
	this.b = this.points[1];
	this.id = lineId++;
	this.distanceToCoords = function (coords) {
		return this.distanceTo(coords.x, coords.y);
	}
	this.distanceTo = function (x, y) {
		var x1 = this.a.x, y1 = this.a.y, x2 = this.b.x,  y2 = this.b.y;
		var a = y1 - y2;
		var b = x2 - x1;
		var c = x2 * y1 - x1 * y2;
		var dist = Math.abs(a * x + b * y - c) / length(a, b);
		var xClosest = (b * (b * x - a * y) + a * c) / lengthSquared(a, b);
		var yClosest = (a * ( -b * x + a * y) + b * c) / lengthSquared(a, b);
		if (isBetween(xClosest, x1, x2)) {
			return dist;
		} else {
			if (x1 < x2 && x < x1 || x1 > x2 && x > x1 
                || x1 == x2 && (y1 < y2 && y < y2 || y1 > y2 && y > y1)) {
				return distance(x, x1, y, y1);
			} else {
				return distance(x, x2, y, y2);
			}
		}
	}
	this.getVectorAB = function () {
		return $V([this.b.x - this.a.x, this.b.y - this.a.y]);
	}
	this.getClosestPoint = function (x, y) {
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
				return {"x": x1, "y": y1};
			} else {
				return {"x": y1, "y": y2};
			}
		}
	}
	this.length = function () {
		return this.points[0].distanceTo(this.points[1]);
	}
}

function isBetween(x, a, b) {
	if (a < b) {
		return a < x && x < b;
	} else {
		return b < x && x < a;
	}
}

function distance(x1, y1, x2, y2) {
	return length(x2 - x1, y2 - y1);
}

function length(x, y) {
	return Math.sqrt(x * x + y * y);
}
function vectorLength(v) {
	return length(v.x, v.y);
}
var lines = [];
lines.push(new Segment(100,100,200,300));
lines.push(new Segment(400,480,700,400));

function lengthSquared(x, y) {
	return x * x + y * y
}

var x, F, DF, b;
// 4 vars per line: s0, s1, d0, d1
// (d0, d1)^T is normalized
// x.length = 4 * lines.length
function recalculate() {
	/*
	lines.forEach(function (l, i) {
		b.push(1);
		F.push(function (x) { return lengthSquared(x[i * 4 + D0], x[i * 4 + D1]); });
		var DFline = [];
		for (var j = 0; j < 4 * lines.length; j++) {
			if (j == i * 4 + D0 || j == i * 4 + D1) {
				DFline.push(
					(function (constVarIndex) {
						return function (x) { return 2 * x[constVarIndex]; }
					})(j))
			} else {
				DFline.push(function (x) { return 0; });
			}
		}
		DF.push(DFline);
	});
	*/
	b = [];
	DF = [];
	x = [];
	F = [];
	constraints.forEach(function (c) {
		console.log(c);
		if (c.type == "parallel") {
			for (var j = 1; j < c.which.length;  j++) {
				constrainAngleAB(lines.indexOf(c.which[0]), lines.indexOf(c.which[j]), 1);
				/*
				b.push(0);
				var lineIndex = lines.indexOf(c.which[j]);
				F.push(function (x) {
					// console.log("lines.indexOf(c.which[j]) = " + lines.indexOf(c.which[j]));
					// console.log("c.which[j] = " + c.which[j]);
					// console.log("j = " + lineIndex);
					return x[4 * lines.indexOf(c.which[0]) + D0] - x[4 * lineIndex + D0];
				});
				b.push(0);
				F.push(function (x) {
					return x[4 * lines.indexOf(c.which[0]) + D1] - x[4 * lineIndex + D1];
				});
				*/
			}
		}
		if (c.type == "perpendicular") {
			//constrainAngle(lines.indexOf(c.first), lines.indexOf(c.second), 0);
		}
		if (c.type == "coincident") {
            var p1Index = pointIndex(c.which[0]);
			for (var j = 1; j < c.which.length;  j++) {
                var p2Index = pointIndex(c.which[j]);
				constrainCoincident(p1Index, p2Index);
			}
		}
	});
	x = [];
	// init x to current values
	lines.forEach(function (line) {
		x.push(line.a.x);
		x.push(line.a.y);
		x.push(line.b.x);
		x.push(line.b.y);
	});
    var lastSize = 0
	for (var count = 0; count < 100000; count++) {
        if ((lastSize = wtfStep()) < 1e-3) {
            enableConsole();
            console.log("breaking at ", count);
            break;
        }
	}
    enableConsole();
    console.log("lastSize ", lastSize);
	for (var j = 0; j < lines.length; j++) {
		lines[j].a.x = x[4 * j + 0];
		lines[j].a.y = x[4 * j + 1];
		lines[j].b.x = x[4 * j + 2];
		lines[j].b.y = x[4 * j + 3];
	}
}
function pointIndex(point) {
    var result = 4 * lines.indexOf(point.line);
    if (point.line.b == point) {
           result += 2;
    }
    return result;
}
function constrainEqual(aIndex, bIndex) {
    console.log("constrainEqual", aIndex,bIndex);
    b.push(0);
	F.push(function (x) { return x[aIndex] - x[bIndex]; });
	var varIndexMap = {};
	varIndexMap[aIndex] = 1;
	varIndexMap[bIndex] = -1;
	var DFline = [];	
	for (var varIndex = 0; varIndex < 4 * lines.length; varIndex++) {
		DFline.push((function (constant) {
			return constant
				? function (x) { return constant; }
				: function (x) { return 0; }
		})(varIndexMap[varIndex]));
	}
	DF.push(DFline);
}
function constrainCoincident(aIndex, bIndex) {
    constrainEqual(aIndex, bIndex);
    constrainEqual(aIndex + 1, bIndex + 1);
}
function calcDFx(F, x, Fx) {
    var DIFF = 1e-9;
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
function constrainAngleSD(lineIndex1, lineIndex2, cosAngle) {
	b.push(cosAngle);
	F.push(function (x) { return x[lineIndex1 * 4 + D0] * x[lineIndex2 * 4 + D0] + x[lineIndex1 * 4 + D1] * x[lineIndex2 * 4 + D1]});
	var varIndexMap = {};
	varIndexMap[lineIndex1 * 4 + D0] = lineIndex2 * 4 + D0;
	varIndexMap[lineIndex2 * 4 + D0] = lineIndex1 * 4 + D0;
	varIndexMap[lineIndex1 * 4 + D1] = lineIndex2 * 4 + D1;
	varIndexMap[lineIndex2 * 4 + D1] = lineIndex1 * 4 + D1;
	var DFline = [];	
	for (var varIndex = 0; varIndex < 4 * lines.length; varIndex++) {
		DFline.push((function (constVarIndex) {
			return constVarIndex
				? function (x) { return x[constVarIndex]; }
				: function (x) { return 0; }
		})(varIndexMap[varIndex]));
	}
	DF.push(DFline);
}
function sumAbsArray (arr) {
    var sum = 0;
    for (var i = 0; i < arr.length; i++) {
        sum += Math.abs(arr[i]);
    }
    return sum;
}
function lengthArray (arr) {
    var sum = 0;
    for (var i = 0; i < arr.length; i++) {
        sum += arr[i] * arr[i];
    }
    return Math.sqrt(sum);
}
function wtfStep() {
    //disableConsole();
    //console.log("x", x);
    Fx = F.map(function (f) { return f(x)});
    //console.log("Fx", Fx);
    //console.log("b", b);
    // DFx = DF.map(function (dfLine) { return dfLine.map(function (dff) { return dff(x) })});
    DFx = calcDFx(F, x, Fx);
    //console.log("DFx", DFx);

    var xDif = [];
    for (var i = 0; i < x.length; i++) {
        
        var sum = 0;
        for (var j = 0; j < b.length; j++) {
            sum += DFx[j][i] * (Fx[j] - b[j]) / lengthArray(DFx[j]);
            //console.log("sumArray", DFx[j],  lengthArray(DFx[j]));
            //	sum += (Fx[j] - b[j]) / DFx[j][i];
        }
        xDif.push(sum / 10);
    }
    //console.log("xdif", xDif);
    
    for (var i = 0; i < x.length; i++) {
        x[i] -= xDif[i];
    }
    return lengthArray(arraySub(Fx, b));
}
function arraySub(a, b) {
    return a.map(function (aEl, aI) { return aEl - b[aI]; });
}

function angleVectors(ax, ay, bx, by) {
    return (ax * bx + ay * by) / Math.sqrt((ax * ax + ay * ay) * (bx * bx + by * by));
}
	
// returns angle between segments AB, CD
function angleABCD(ax, ay, bx, by, cx, cy, dx, dy) {
    return angleVectors(bx - ax, by - ay, dx - cx, dy - cy);
	// return ((bx - ax) * (dx - cx) + (by - ay) * (dy - cy)) / Math.sqrt(((bx - ax)*(bx - ax)+(by - ay)*(by - ay))*((dx - cx)*(dx - cx)+(dy - cy)*(dy - cy)))
}
function constrainAngleAB(lineIndex1, lineIndex2, cosAngle) {
	b.push(cosAngle);
	F.push(function (x) { return angleABCD(
		x[lineIndex1 * 4 + 0],
		x[lineIndex1 * 4 + 1],
		x[lineIndex1 * 4 + 2],
		x[lineIndex1 * 4 + 3],
		x[lineIndex2 * 4 + 0],
		x[lineIndex2 * 4 + 1],
		x[lineIndex2 * 4 + 2],
		x[lineIndex2 * 4 + 3]
	)});
	
	var varIndexMap = {};
	varIndexMap[lineIndex1 * 4 + 0] = [lineIndex2 * 4 + 0, lineIndex2 * 4 + 2];
	varIndexMap[lineIndex1 * 4 + 1] = [lineIndex2 * 4 + 1, lineIndex2 * 4 + 3];
	varIndexMap[lineIndex1 * 4 + 2] = [lineIndex2 * 4 + 2, lineIndex2 * 4 + 0];
	varIndexMap[lineIndex1 * 4 + 3] = [lineIndex2 * 4 + 3, lineIndex2 * 4 + 1];
	varIndexMap[lineIndex2 * 4 + 0] = [lineIndex1 * 4 + 0, lineIndex1 * 4 + 2];
	varIndexMap[lineIndex2 * 4 + 1] = [lineIndex1 * 4 + 1, lineIndex1 * 4 + 3];
	varIndexMap[lineIndex2 * 4 + 2] = [lineIndex1 * 4 + 2, lineIndex1 * 4 + 0];
	varIndexMap[lineIndex2 * 4 + 3] = [lineIndex1 * 4 + 3, lineIndex1 * 4 + 1];
	var DFline = [];	
	for (var varIndex = 0; varIndex < 4 * lines.length; varIndex++) {
		DFline.push((function (constVarIndexes) {
			return constVarIndexes
				? function (x) { return x[constVarIndexes[0]], - x[constVarIndexes[1]]; }
				: function (x) { return 0; }
		})(varIndexMap[varIndex]));
	}
	DF.push(DFline);
}
function constrainAngleAB(lineIndex1, lineIndex2, cosAngle) {
	b.push(cosAngle);
	F.push(function (x) { return angleABCD(
		x[lineIndex1 * 4 + 0],
		x[lineIndex1 * 4 + 1],
		x[lineIndex1 * 4 + 2],
		x[lineIndex1 * 4 + 3],
		x[lineIndex2 * 4 + 0],
		x[lineIndex2 * 4 + 1],
		x[lineIndex2 * 4 + 2],
		x[lineIndex2 * 4 + 3]
	)});
	
	var varIndexMap = {};
	varIndexMap[lineIndex1 * 4 + 0] = [lineIndex2 * 4 + 0, lineIndex2 * 4 + 2];
	varIndexMap[lineIndex1 * 4 + 1] = [lineIndex2 * 4 + 1, lineIndex2 * 4 + 3];
	varIndexMap[lineIndex1 * 4 + 2] = [lineIndex2 * 4 + 2, lineIndex2 * 4 + 0];
	varIndexMap[lineIndex1 * 4 + 3] = [lineIndex2 * 4 + 3, lineIndex2 * 4 + 1];
	varIndexMap[lineIndex2 * 4 + 0] = [lineIndex1 * 4 + 0, lineIndex1 * 4 + 2];
	varIndexMap[lineIndex2 * 4 + 1] = [lineIndex1 * 4 + 1, lineIndex1 * 4 + 3];
	varIndexMap[lineIndex2 * 4 + 2] = [lineIndex1 * 4 + 2, lineIndex1 * 4 + 0];
	varIndexMap[lineIndex2 * 4 + 3] = [lineIndex1 * 4 + 3, lineIndex1 * 4 + 1];
	var DFline = [];	
	for (var varIndex = 0; varIndex < 4 * lines.length; varIndex++) {
		DFline.push((function (constVarIndexes) {
			return constVarIndexes
				? function (x) { return x[constVarIndexes[0]], - x[constVarIndexes[1]]; }
				: function (x) { return 0; }
		})(varIndexMap[varIndex]));
	}
	DF.push(DFline);
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

function sigma (c, s) {
    if ( 0 == c) {
        return 1;
    }
    if ( abs(s) < abs(c)) {
        return 0.5 * sgn(c) * s;
    }
    return 2 * sgn (s) / c;
}

function clicky2() {
	Fx = F.map(function (f) { return f(x)});
	console.log("Fx", Fx);
	DFx = DF.map(function (dfLine) { return dfLine.map(function (dff) { return dff(x) })});
	console.log("DFx", DFx);
    
    for (var row = 1; row < b.length; row++) {
        for (var col = 0; col < row; col++) {
            console.log("row ", row, "col ", col);
            var xi = DFx[col][col];
            var xk = DFx[row][col];
            if (xk == 0) {
                continue;
            }
            var r = Math.sqrt(xi * xi + xk * xk);
            var c = xi / r;
            var s = xk / r;
            
            // apply transformation on every column:
            for (var col2 = col; col2 < x.length; col2++) {
                var x1 = DFx[col][col2] * c + DFx[row][col2] * s;
                var x2 = DFx[row][col2] * c - DFx[col][col2] * s;
                DFx[col][col2] = x1;
                DFx[row][col2] = x2;
                console.log("r ", r, "c ", c, "s ", s);
                console.log("DFx ", DFx);
            }
        }
    }

}
var constraints = [		{"type": "parallel", "which": [lines[0], lines[1]]},
						{"type": "perpendicular", "first": "b", "second": "c"},
						{"type": "coincident", "which": [lines[0].b, lines[1].a]}];
function relMouseCoords(e){
    var totalOffsetX = 0;
    var totalOffsetY = 0;
    var canvasX = 0;
    var canvasY = 0;
    var currentElement = this;

    do{
        totalOffsetX += currentElement.offsetLeft - currentElement.scrollLeft;
        totalOffsetY += currentElement.offsetTop - currentElement.scrollTop;
    }
    while(currentElement = currentElement.offsetParent)

    canvasX = e.pageX - totalOffsetX;
    canvasY = e.pageY - totalOffsetY;

    return {x:canvasX, y:canvasY}
}
HTMLCanvasElement.prototype.relMouseCoords = relMouseCoords;
function paintSegments() {
    ctx.clearRect (0, 0, ctx.canvas.width, ctx.canvas.height);
	ctx.moveToPoint = function (p) {
		this.moveTo(p.x, p.y);
	}
	ctx.moveToVector = function (p) {
		this.moveTo(p.e(1), p.e(2));
	}
	ctx.lineToVector = function (p) {
		this.lineTo(p.e(1), p.e(2));
	}
	ctx.lineToPoint = function (p) {
		this.lineTo(p.x, p.y);
	}
	ctx.arcPoint = function (p, radius, start, end, clockwise) {
		this.arc(p.x, p.y, radius, start, end, clockwise);
	}
ctx.fillStyle="rgb(100, 100, 255)";
ctx.lineWidth=2;
	lines.forEach(function (line) {
		//console.log("line", line);
		console.log("hoverHighlight.length", hoverHighlight.length);
		ctx.strokeStyle=colorFor(highlighted.contains(line) || hoverHighlight.contains(line), selected.contains(line));
		ctx.beginPath();
		ctx.moveToPoint(line.a);
		ctx.lineToPoint(line.b);
		ctx.stroke();
		ctx.beginPath();
		ctx.arcPoint(line.a, 2, 0, 2 * Math.PI, false);
		ctx.arcPoint(line.b, 2, 0, 2 * Math.PI, false);
		ctx.fill();
	});
	paintConstraints();
	// console.log("linecount", lines.length);
}
function colorFor(highlighted, selected) {
	return !selected
		? (!highlighted ? "#33CCFF" : "#145266")
		: (!highlighted ? "#FF3399" : "#330A1E")
}	
function getAllPoints() {
	return [].concat.apply([], lines.map(function (segment) { return segment.points; }));
}
var highlighted = [], selected = [];
function paintConstraints() {
	var crossCount = 2;
    ctx.lineWidth=1;
	constraints.forEach(function (constraint) {
		switch (constraint.type) {
			case "coincident":
                ctx.beginPath();
				ctx.arcPoint(constraint.which[0], 4, 0, 2 * Math.PI, false);
				ctx.stroke();
				break;
			case "parallel":
				
				ctx.beginPath();
				for (var c = 0; c < constraint.which.length; c++) {
					var line = constraint.which[c];
					var ab = line.getVectorAB();
					var abLength = ab.modulus();
					var abNormalized = ab.toUnitVector();
					var abPerpendicular = $V([ -ab.e(2), ab.e(1) ]).toUnitVector();
					var crossLength = 10;
					var crossSpacing = 3;
					for (var i = 0; i < crossCount; i++) {
						var s = abLength / 2 - crossSpacing * crossCount / 2 + i * crossSpacing;
						var crossStart = line.a.$V().add(abNormalized.multiply(s)).add(abPerpendicular.multiply(crossLength / 2));
						var crossEnd = line.a.$V().add(abNormalized.multiply(s)).add(abPerpendicular.multiply(-crossLength / 2));
						//console.log("crosspos", crossStart, crossEnd);
						ctx.moveToVector(crossStart);
						ctx.lineToVector(crossEnd);
					}
					ctx.stroke();
				}
				crossCount++;
			default:
				// console.log("errror", constraint.type);
		}
		
	});
}
function normalizeVector(v) {
	var vLength = vectorLength(v);
	v.x /= vLength;
	v.y /= vLength;
}
var currentAddingSegment, connectingConstraint = [], currentPointConstraint, hoverHighlight = [];
// console.log = oldConsole;
var cvs, ctx;
var mode = "";
var isShift;
window.onload = function () {
	cvs = document.getElementById("myCanvas");
	cvs.onclick = function (e) {
		var mousePos = cvs.relMouseCoords(e);
        if (mode == "addsegments") {
			currentPointConstraint = undefined;
			lastSegment = currentAddingSegment;
            currentAddingSegment = new Segment(mousePos.x, mousePos.y, mousePos.x, mousePos.y);
            lines.push(currentAddingSegment);
            if (lastSegment) {
                constraints.push(connectingConstraint = {"type":"coincident", "which":[lastSegment.b, currentAddingSegment.a]});
            }
            //console.log("e", e);
			console.log("lines", lines);
			console.log("constraints", constraints);
			recalculate();
        } else {
			if (!isShift) {
				selected = hoverHighlight.slice();
			} else {
				if (selected.contains(hoverHighlight[0])) {
					selected.remove(hoverHighlight[0]);
				} else {
					selected.push(hoverHighlight[0]);
				}
			}
			console.log ("selected is now ", selected);
			console.log ("highlighted is now ", highlighted);
			paintSegments();
		}
	}
	window.onkeydown = function (ev) {
  var key;
  if (window.event) {
    key = window.event.keyCode;
    isShift = !!window.event.shiftKey; // typecast to boolean
  } else {
    key = ev.which;
    isShift = !!ev.shiftKey;
  }
}
	window.onkeyup = function (ev) {
  var key;
  if (window.event) {
    key = window.event.keyCode;
    isShift = !!window.event.shiftKey; // typecast to boolean
  } else {
    key = ev.which;
    isShift = !!ev.shiftKey;
  }
}
	cvs.onmousemove = function (e) {
		var mousePos = cvs.relMouseCoords(e);
        if (mode == "addsegments") {
			constraints.remove(currentPointConstraint);
            //console.log(mousePos);
			var points = getAllPoints();
			points.removeAll(currentAddingSegment.points);
			var pointDistances = points
				.map(function (point) { return {"point": point, "distance": point.distanceTo(mousePos)}; })
				.filter(function (pair) { return pair.distance < 16; })
				.sort(function (a, b) { return a.distance - b.distance; });
			if (pointDistances.length != 0) {
				currentPointConstraint = {"type": "coincident", "which": [currentAddingSegment.b, pointDistances[0].point]}
				constraints.push(currentPointConstraint);
				mousePos = pointDistances[0].point;
			}
			if (currentAddingSegment) {
				currentAddingSegment.b.x = mousePos.x;
				currentAddingSegment.b.y = mousePos.y;
			}
			paintSegments();
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
			hoverHighlight = [];
			var linesByDistance = lines
				.map(function (line) { return {"line": line, "distance": line.distanceToCoords(mousePos)}; })
				.filter(function (pair) { return pair.distance < 16; })
				.sort(function (a, b) { return a.distance - b.distance; });
			console.log("linesByDistance", linesByDistance);
			if (0 != linesByDistance.length) {
				hoverHighlight.push(linesByDistance[0].line);
			}
			paintSegments();
		}
	}
    document.getElementById("clearmode").onclick = function () {
        lines.remove(currentAddingSegment);
        constraints.remove(connectingConstraint);
        currentAddingSegment = undefined;
        console.log("mode cleared");
        paintSegments();
        mode = "";
    }
	ctx = cvs.getContext("2d");
	recalculate();
	paintSegments();
	console.log("distanceto", lines[0].distanceTo(340, 354));
	console.log("getClosestPoint", lines[0].getClosestPoint(340, 354));
	console.log("line", lines[0].a);
}
function clicky() {
	for (var j = 0; j < 1; j++) {
		wtfStep();
	}
	for (var j = 0; j < lines.length; j++) {
		lines[j].a.x = x[4 * j + 0];
		lines[j].a.y = x[4 * j + 1];
		lines[j].b.x = x[4 * j + 2];
		lines[j].b.y = x[4 * j + 3];
	}
	paintSegments();
}
function makeParallel() {
	if (selected.length == 2) {
		constraints.push({"type": "parallel", "which": selected.slice()});
		recalculate();
		paintSegments();
	}
}



