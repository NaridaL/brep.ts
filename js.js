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
}
SegmentEndPoint.prototype = {
	"distanceToCoords": function (p) {
		return distance(this.x, this.y, p.x, p.y);
	},
	"$V": function () {
		return $V([this.x, this.y]);
	},
    "isConstrained": function () {
        var line = this.line;
        return constraints.some(
            function (constraint) {
                return constraint.type == "parallel" && constraint.which.contains(line)
				|| constraint.type == "perpendicular" && (constraint.first == line || constraint.second == line)
                || constraint.type == "pointOnLine" && (constraint.line == line || constraint.point == this); 
            });
    }
}
function makeCoincident(p1, p2) {
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
        constraints.push(p1.coincidence = p2.coincidence = {"type":"coincident", "which":[p1,p2]});
    }
    p2.x = p1.x;
    p2.y = p1.y;
}
function removeCoincidence(p) {
    p.coincidence.which.remove(p);
    if (p.coincidence.which.length == 1) {
        constraints.remove(p.coincidence);
        p.coincidence.which[0].coincidence = undefined
    }
    p.coincidence = undefined;
}
function Segment(x1, y1, x2, y2) {
	this.points = [new SegmentEndPoint(x1, y1, this), new SegmentEndPoint(x2, y2, this)];
	this.a = this.points[0];
	this.b = this.points[1];
	this.id = lineId++;
}
Segment.prototype = {
	"distanceToCoords": function (coords) {
		return this.distanceTo(coords.x, coords.y);
	},
    "getS": function (v) {
        if (this.b.x - this.a.x > this.b.y - this.a.y) {
            return (v.e(1) - this.a.x) / (this.b.x - this.a.x);
        } else {
            return (v.e(2) - this.a.y) / (this.b.y - this.a.y);
        }
    },
	"distanceTo": function (x, y) {
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
	},
	"getVectorAB": function () {
		return $V([this.b.x - this.a.x, this.b.y - this.a.y]);
	},
	"getClosestPoint": function (x, y) {
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
	},
	"length": function () {
		return this.points[0].distanceTo(this.points[1]);
	},
    "intersection": function (segment) {
        return intersection(this.a.x, this.a.y, this.b.x, this.b.y,
            segment.a.x, segment.a.y, segment.b.x, segment.b.y);
    }
}
function intersection(x1, y1, x2, y2, x3, y3, x4, y4) {
    var denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    var xNominator = (x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4);
    var yNominator = (x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4);
    return $V([xNominator / denominator, yNominator / denominator]);
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
lines.push(new Segment(100,100,200,700));
lines.push(new Segment(400,480,700,400));

function lengthSquared(x, y) {
	return x * x + y * y
}

var x, F, DF, b, varMap;
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
    varMap = new Map();
	// init x to current values
    var xIndex = 0;
    function pushPoint(p) {
        if (varMap.has(p) /*    || !p.isConstrained()*/) {
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
	lines.forEach(function (line) {
        pushPoint(line.a);
        pushPoint(line.b);
        console.log("varMap.get(line.a]", varMap.get(line.a));
        console.log("varMap.get(line.b]", varMap.get(line.b));  
	});
	constraints.forEach(function (c) {
		console.log(c);
		if (c.type == "parallel") {
			for (var j = 1; j < c.which.length;  j++) {
				constrainAngleAB(c.which[0], c.which[j], 1);
			}
		}
		if (c.type == "perpendicular") {
				constrainAngleAB(c.first, c.second, 0);
		}
        if (c.type == "colinear") {
            constrainAngleAB(c.lines[0], c.lines[1], 1);
            constrainPointOnLine(
                varMap.get(c.lines[0].a), varMap.get(c.lines[0].a) + 1,
                varMap.get(c.lines[0].b), varMap.get(c.lines[0].b) + 1,
                varMap.get(c.lines[1].a), varMap.get(c.lines[1].a) + 1);
        }
		if (c.type == "pointDistance") {
			constrainDistance(c.constrains[0], c.constrains[1], c.distance);
		}
		if (c.type == "coincident") {
            /*var p1Index = pointIndex(c.which[0]);
			for (var j = 1; j < c.which.length;  j++) {
                var p2Index = pointIndex(c.which[j]);
				constrainCoincident(p1Index, p2Index);
			}*/
		}
		if (c.type == "pointOnLine") {
            constrainPointOnLine(
                varMap.get(c.line.a), varMap.get(c.line.a) + 1,
                varMap.get(c.line.b), varMap.get(c.line.b) + 1,
                varMap.get(c.point), varMap.get(c.point) + 1);
                
		}
        console.log(b);
	});
	
//lines.push(new Segment(100,100,200,700));
//lines.push(new Segment(400,480,700,400));
/*
	var xsToConstrain = [0, 1, 3, 6, 7];
	var originalFLength = F.length;
	for (var i = originalFLength; i < x.length; i++) {
		var varIndexToConstrainEqual = xsToConstrain[i - originalFLength];
		b.push(x[varIndexToConstrainEqual]);
		F.push(
			(function (xIndex) { return function (x) { return x[xIndex]; }; })(varIndexToConstrainEqual));		
	}
    */
    var lastSize = 0
	for (var count = 0; count < 0; count++) {
        if ((lastSize = newtonStep()) < 1e-3) {
            enableConsole();
            console.log("breaking at ", count);
            break;
        }
	}
    enableConsole();
    console.log("lastSize ", lastSize);
	reverse();
}
function gaussNewtonStep() {
    Fx = F.map(function (f) { return f(x)});
    console.log("x", x);
    console.log("Fx", Fx);
    console.log("b", b);
	var jacobi = $M(calcDFx(F, x, Fx));
	console.log("jacobi", jacobi);
	var jacobiTranspose = jacobi.transpose();
    var matrix = jacobiTranspose.multiply((jacobi.multiply(jacobiTranspose)).inverse());
	console.log("inverseJacobi",matrix);
    console.log($V(x).subtract($V(b)));
	var xDif = matrix.multiply($V(Fx).subtract($V(b)));
    console.log("xDif", xDif);
	x = $V(x).subtract(xDif).elements;
}
function newtonStep() {
    gaussNewtonStep(); return;
    Fx = F.map(function (f) { return f(x)});
    console.log("x", x);
    console.log("Fx", Fx);
    console.log("b", b);
	var jacobi = $M(calcDFx(F, x, Fx));
	console.log("jacobi", jacobi);
	var inverseJacobi = jacobi.inverse();
	console.log("inverseJacobi",inverseJacobi);
    console.log($V(x).subtract($V(b)));
	var xDif = inverseJacobi.multiply($V(Fx).subtract($V(b)));
    console.log("xDif", xDif);
	x = $V(x).subtract(xDif).elements;
}
function wtfStep() {
    disableConsole();
    console.log("x", x);
    Fx = F.map(function (f) { return f(x)});
    console.log("Fx", Fx);
    console.log("b", b);
    var diffs = calcDFx2(F, x, Fx);
    console.log("diffs", diffs);
	/*
    var maxDiff = 0;
    var useXIndex;
    for (var i = 0; i < x.length; i++) {
        if (diffs[i] > maxDiff){ 
            maxDiff = diffs[i];
            useXIndex = i;
        }
    }
	*/
    var currLength = $V(Fx).subtract($V(b)).modulus();
    var oldLength = currLength;
    var xDif = $V(diffs).multiply(10);
    var newX, newLength;
    var halvedCount = -1;
    do {
        newX = $V(x);
		newX = newX.subtract(xDif);
        xDif = xDif.multiply(1/2);
        halvedCount++;
        var newFx = F.map(function (f) { return f(newX.elements)});
        newLength = $V(newFx).subtract($V(b)).modulus();
    } while (newLength > 1e-6 && newLength >= currLength);
    x = newX.elements;
    console.log("xdif", xDif);
    console.log("currLength", $V(Fx).subtract($V(b)).modulus(), "halved", halvedCount);
    enableConsole();
    
    return $V(Fx).subtract($V(b)).modulus();
    // DFx = DF.map(function (dfLine) { return dfLine.map(function (dff) { return dff(x) })});
    /*DFx = calcDFx(F, x, Fx);
    
	var DFxM = $M(DFx);
    console.log("DFx", DFx);

    var xDif = [];
	// iterate over x values
	var minAngle = 10, changeVarIndex;
    for (var i = 0; i < x.length; i++) {
        var angle = $V(b).subtract($V(Fx)).angleFrom($V(DFxM.col(i)));
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
	xDif[changeVarIndex] = $V(b).subtract($V(Fx)).modulus() - DFxM.col(i).modulus();
    */
}
function constrainPointOnLine(axIndex, ayIndex, bxIndex, byIndex, cxIndex, cyIndex) {
    b.push(0);
    var indexArray = Array.slice(arguments);
    F.push(function (x) { 
        return blurgh.apply(null, indexArray.map(function (index) { return x[index]; })); });
}

function blurgh2(ax, ay, bx, by, cx, cy) {
    // x1 / y1 == x2 / y2 => x1 * y2 - x2 * y1 = 0
    var x1 = cx - ax, y1 = cy - ay, x2 = cx - bx, y2 = cy - by;
    return x1 * y2 - x2 * y1;
}
function blurgh(x1, y1, x2, y2, x, y) {
    var a = y1 - y2;
    var b = x2 - x1;
    var c = x2 * y1 - x1 * y2;
    var dist = Math.abs(a * x + b * y - c) / length(a, b);
    return dist;
}
function constrainDistance(p0, p1, pDistance) {
    b.push(1);
    F.push(function (x) { return distance(
        x[varMap.get(p0)], x[varMap.get(p0) + 1], 
        x[varMap.get(p1)], x[varMap.get(p1) + 1]) / pDistance; });
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
function calcDFx2(F, x, Fx) {
    var DIFF = 1e-9;
    var DFx = new Array(x.length);
    var currLength = $V(Fx).subtract($V(b)).modulus();
    for (var i = 0; i < x.length; i++) {
        x[i] += DIFF;
        var newLength = $V(F.map(function (f) { return f(x); })).subtract($V(b)).modulus();
        DFx[i] = (newLength - currLength) / DIFF;
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

function angleVectors(ax, ay, bx, by) {
//    console.log(ax, ay, bx, by);
    return Math.abs(ax * bx + ay * by) / Math.sqrt((ax * ax + ay * ay) * (bx * bx + by * by));
}
	
// returns angle between segments AB, CD
function angleABCD(ax, ay, bx, by, cx, cy, dx, dy) {
    return angleVectors(bx - ax, by - ay, dx - cx, dy - cy);
	// return ((bx - ax) * (dx - cx) + (by - ay) * (dy - cy)) / Math.sqrt(((bx - ax)*(bx - ax)+(by - ay)*(by - ay))*((dx - cx)*(dx - cx)+(dy - cy)*(dy - cy)))
}
/*
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
*/
function constrainAngleAB(line1, line2, cosAngle) {
	b.push(cosAngle);
	F.push(function (x) { 
        //console.log("blaaah", varMap.get(line1.b), x[varMap.get(line1.b)]);
        return angleABCD(
		x[varMap.get(line1.a)],
		x[varMap.get(line1.a) + 1],
		x[varMap.get(line1.b)],
		x[varMap.get(line1.b) + 1],
		x[varMap.get(line2.a)],
		x[varMap.get(line2.a) + 1],
		x[varMap.get(line2.b)],
		x[varMap.get(line2.b) + 1]
	); });
	/*
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
    */
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
var constraints = [		//{"type": "parallel", "which": [lines[0], lines[1]]},
						//{"type": "perpendicular", "first": lines[0], "second": lines[1]},
                        //{"type": "pointOnLine", "line": lines[0], "point": lines[1].a},
                            {"type": "colinear", "lines": [lines[0], lines[1]]},
                            {"type":"pointDistance", "constrains": [lines[1].a, lines[1].b], "distance": 400}];
// makeCoincident(lines[0].b, lines[1].a);
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
    ctx.fillStyle="rgb(100, 100, 255)";
    ctx.lineWidth=2;
	lines.forEach(function (line) {
		//console.log("line", line);
		//console.log("hoverHighlight.length", hoverHighlight.length);
		ctx.beginPath();
		ctx.moveToPoint(line.a);
		ctx.lineToPoint(line.b);
		ctx.strokeStyle=colorFor(highlighted.contains(line) || hoverHighlight.contains(line), selected.contains(line));
		ctx.stroke();
		ctx.beginPath();
		ctx.arcPoint(line.a, 2, 0, 2 * Math.PI, false);
		ctx.fillStyle=colorFor(highlighted.contains(line.a) || hoverHighlight.contains(line.a), selected.contains(line.a));
		ctx.fill();
		ctx.beginPath();
		ctx.arcPoint(line.b, 2, 0, 2 * Math.PI, false);
		ctx.fillStyle=colorFor(highlighted.contains(line.b) || hoverHighlight.contains(line.b), selected.contains(line.b));
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
function getPerpendicular(v) {
    return $V([ -v.e(2), v.e(1) ]);
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
					var ab90 = $V([ -ab.e(2), ab.e(1) ]).toUnitVector();
					var crossLength = 10;
					var crossSpacing = 3;
					for (var i = 0; i < crossCount; i++) {
						var s = abLength / 2 - crossSpacing * crossCount / 2 + i * crossSpacing;
						var crossStart = line.a.$V().add(abNormalized.multiply(s)).add(ab90.multiply(crossLength / 2));
						var crossEnd = line.a.$V().add(abNormalized.multiply(s)).add(ab90.multiply(-crossLength / 2));
						//console.log("crosspos", crossStart, crossEnd);
						ctx.moveToVector(crossStart);
						ctx.lineToVector(crossEnd);
					}
					ctx.stroke();
				}
				crossCount++;
                break;
            case "perpendicular":
                var intersection = constraint.first.intersection(constraint.second);
                var abPos = constraint.first.getS(intersection);
                var cdPos = constraint.second.getS(intersection);
                var ab = constraint.first.getVectorAB().toUnitVector().multiply(0.5 < abPos ? -16 : 16);
                var cd = constraint.second.getVectorAB().toUnitVector().multiply(0.5 < cdPos ? -16 : 16);
                ctx.beginPath();
                ctx.moveToVector(intersection.add(ab));
                ctx.lineToVector(intersection.add(ab).add(cd));
                ctx.lineToVector(intersection.add(cd));
                ctx.stroke();
                break;
            case "pointDistance":
                var a = constraint.constrains[0].$V(), b = constraint.constrains[1].$V();
                var ab = b.subtract(a);
                var ab1 = ab.toUnitVector();
                var abLength = ab.modulus();
                var ab90 = getPerpendicular(ab1);
                var textLength = 60;
                ctx.beginPath();
                ctx.moveToVector(a.add(ab90.x(6)));
                ctx.lineToVector(a.add(ab90.x(22)));
                ctx.moveToVector(b.add(ab90.x(6)));
                ctx.lineToVector(b.add(ab90.x(22)));
                
                ctx.moveToVector(a.add(ab90.x(14)));
                ctx.lineToVector(a.add(ab90.x(14)).add(ab1.x(abLength / 2 - textLength / 2)));
                ctx.moveToVector(a.add(ab90.x(14)).add(ab1.x(abLength / 2 + textLength / 2)));
                ctx.lineToVector(a.add(ab90.x(14)).add(ab1.x(abLength)));
                ctx.lineWidth = 1;
                ctx.stroke();
                ctx.font = "14px sans";
                ctx.textAlign = "center";
                ctx.textBaseline = "middle";
                var angle = Math.atan(ab.e(2) / ab.e(1));
                var textCenter = a.add(ab90.x(14)).add(ab1.x(abLength / 2))
                ctx.save();
                ctx.translate(textCenter.e(1), textCenter.e(2));
                ctx.rotate(angle);
                ctx.textAlign = "center";
                ctx.fillText(constraint.distance, 0,0);
                ctx.restore();
                break;
            case "pointOnLine": 
                var ab = constraint.line.getVectorAB();
                var ab1 = ab.toUnitVector();
                var p = constraint.point.$V();
                ctx.beginPath();
                ctx.moveToVector(p.add(ab1.x(16)));
                ctx.lineToVector(p.add(ab1.x(4)));
                ctx.moveToVector(p.add(ab1.x(-16)));
                ctx.lineToVector(p.add(ab1.x(-4)));
                ctx.lineWidth = 4;
                ctx.stroke();
                break;
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
                makeCoincident(lastSegment.b, currentAddingSegment.a);
            }
            //console.log("e", e);
			console.log("lines", lines);
			console.log("constraints", constraints);
            console.log(currentAddingSegment.a);
			recalculate();
            console.log(currentAddingSegment.a);
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
			currentAddingSegment && points.removeAll(currentAddingSegment.points);
			var pointDistances = points
				.map(function (point) { return {"point": point, "distance": point.distanceToCoords(mousePos)}; })
				.filter(function (pair) { return pair.distance < 16; })
				.sort(function (a, b) { return a.distance - b.distance; });
			if (pointDistances.length != 0) {
                //makeCoincident(currentAddingSegment.b, pointDistances[0].point);
				currentPointConstraint = pointDistances[0].point;
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
			var elsByDistance = lines.concat(getAllPoints())
				.map(function (el) { return {"el": el, "distance": el.distanceToCoords(mousePos)}; })
				.filter(function (pair) { return pair.distance < 16; })
				.sort(function (a, b) { return chainComparison(pointsBeforeLines(a.el, b.el), a.distance - b.distance); });
			//console.log("elsByDistance", elsByDistance);
			if (0 != elsByDistance.length) {
				hoverHighlight.push(elsByDistance[0].el);
			}
			paintSegments();
		}
	}
    document.getElementById("clearmode").onclick = function () {
        if (connectingConstraint) {
            removeCoincidence(currentAddingSegment.a);
        }
        if (currentAddingSegment.b.coincidence) {
            removeCoincidence(currentAddingSegment.b.coincidence);
        }
        lines.remove(currentAddingSegment);
        currentAddingSegment = undefined;
        console.log("mode cleared");
        paintSegments();
        mode = "";
    }
	ctx = cvs.getContext("2d");
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
	recalculate();
	paintSegments();
	console.log("distanceto", lines[0].distanceTo(340, 354));
	console.log("getClosestPoint", lines[0].getClosestPoint(340, 354));
	console.log("line", lines[0].a);
}
function chainComparison(diff1, diff2) {
    return diff1 != 0 ? diff1 : diff2;
}
function reverse() {
	varMap.forEach(function (pointIndex, point) {
		point.x = x[pointIndex];
		point.y = x[pointIndex + 1 ];
	});
}
function pointsBeforeLines(a, b) {
    console.log((a instanceof SegmentEndPoint ? -1 : 1) - (b instanceof SegmentEndPoint ? -1 : 1));
    return (a instanceof SegmentEndPoint ? -1 : 1) - (b instanceof SegmentEndPoint ? -1 : 1);
}
function clicky() {
	for (var j = 0; j < 1; j++) {
		newtonStep();
	}
	reverse();
	paintSegments();
}
function getParallelConstraint(segment) {
    return constraints.find(function (c) { c.type == "parallel" && c.which.contains(segment); });
}
function makeParallel() {
	var selSegments = selected.filter(function (s) { return s instanceof Segment; });
    console.log(selSegments);
    if (selSegments.length >= 2) {
        var parallelConstraint = getParallelConstraint(selSegments[0]);
        if (!parallelConstraint) {
            parallelConstraint = {"type": "parallel", "which": [selSegments[0]]};
        }
        for (var i = 1; i < selSegments.length; i++) {
            var seg = selSegments[i];
            var segConstraint = getParallelConstraint(seg);
            if (segConstraint) {
                if (segConstraint != parallelConstraint) {
                    segConstraint.which.forEach(function (segConstraintOther) {
                        parallelConstraint.which.push(segConstraintOther);
                    });
                    constraints.remove(segConstraint);
                }
            } else {
                parallelConstraint.which.push(seg);
            }   
        }
        constraints.push(parallelConstraint);
        recalculate();
        paintSegments();
    }
    console.log("parallelConstraint", parallelConstraint);
}

function makePerpendicular() {
	var selSegments = selected.filter(function (s) { return s instanceof Segment; });
    console.log(selSegments);
    if (selSegments.length == 2) {
        constraints.push({"type":"perpendicular", "first":selSegments[0], "second":selSegments[1]});
        recalculate();
        paintSegments();
    }
}
function makeSelCoincident() {
	var selPoints = selected.filter(function (s) { return s instanceof SegmentEndPoint; });
    if (selPoints.length >= 2) {
        for (var i = 1; i < selPoints.length; i++) {
            makeCoincident(selPoints[0], selPoints[i]);
        }
        recalculate();
        paintSegments();
    }
}
function selPointOnLine() {
    if (selected.length >= 2 &&
        ((selected[0] instanceof SegmentEndPoint && selected[1] instanceof Segment) 
        || (selected[1] instanceof SegmentEndPoint && selected[0] instanceof Segment))) {
        var point = selected[0] instanceof SegmentEndPoint ? selected[0] : selected[1];
        var segment = selected[0] instanceof Segment ? selected[0] : selected[1];
        constraints.push({"type":"pointOnLine", "line": segment, "point": point});
        recalculate();
        paintSegments();
    }
}
function makeSelColinear() {
	var selSegments = selected.filter(function (s) { return s instanceof Segment; });
    console.log(selSegments);
    if (selSegments.length >= 2) {
        var colinearConstraint = getColinearConstraint(selSegments[0]);
        if (!parallelConstraint) {
            colinearConstraint = {"type": "colinear", "which": [selSegments[0]]};
        }
        for (var i = 0; i < selSegments.length; i++) {
            var seg = selSegments[i];
            var segConstraint = getParallelConstraint(seg);
            if (segConstraint) {
                if (segConstraint != colinearConstraint) {
                    segConstraint.which.forEach(function (segConstraintOther) {
                        colinearConstraint.which.push(segConstraintOther);
                    });
                    constraints.remove(segConstraint);
                }
            } else {
                colinearConstraint.which.push(seg);
            }   
        }
        constraints.push(colinearConstraint);
        recalculate();
        paintSegments();
    }    
}

