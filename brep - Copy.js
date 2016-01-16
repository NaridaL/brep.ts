
["min", "max", "PI", "sqrt","pow","round"].forEach(function (propertyName) {
	/*if (window[propertyName]) {
	 throw new Error("already exists"+propertyName)
	 }*/
	window[propertyName] = Math[propertyName];
});
if (!Array.prototype.find) {
	Array.prototype.find = function(predicate) {
		if (this === null) {
			throw new TypeError('Array.prototype.find called on null or undefined');
		}
		if (typeof predicate !== 'function') {
			throw new TypeError('predicate must be a function');
		}
		var list = Object(this);
		var length = list.length >>> 0;
		var thisArg = arguments[1];
		var value;

		for (var i = 0; i < length; i++) {
			value = list[i];
			if (predicate.call(thisArg, value, i, list)) {
				return value;
			}
		}
		return undefined;
	};
}
function BREP(vertices, edges, faces) {
	this.vertices = vertices;
	this.edges = edges;
	this.faces = faces;
}
var M4 = NLA.Matrix4x4, V3 = NLA.Vector3, P3 = NLA.Plane3, L3 = NLA.Line3, cl = console.log, assert = NLA.assert
BREP.box = function (w, h, d) {
	/*
	 6 - 7
	 2 - 3 |
	 | 4 | 5
	 0 - 1
	 */
	var baseVertices = [
		V3(0, 0, 0),
		V3(0, h, 0),
		V3(w, h, 0),
		V3(w, 0, 0),
	];
	return BREP.extrude(baseVertices, P3(V3.Z.negated(), 0), V3(0, 0, d))
}
var faceId = 0
BREP.Face = function(vertices, plane) {
	assert(!plane || plane instanceof P3, "plane is not a P3")
	this.vertices = vertices
	// this.aabb
	this.plane = plane ? plane : P3.throughPoints(vertices[0], vertices[1], vertices[2])
	assert(!plane || vertices.every((v) => plane.containsPoint(v)), "all vertices must be contained within the plane")
	var ccw = isCCW(vertices, this.plane.normal)
	assert(!plane || ccw, "vertices are not ccw to passed plane")
	if (!ccw) {
		this.plane = this.plane.flipped()
	}
	this.id = faceId++
}
BREP.Face.prototype = {
	toString: function () {
		return "[" + this.vertices.map((v) => v.toString()).join(", ") + "] on " + this.plane.toString()
	},
	equals: function (face) {

		return this.vertices.length == face.vertices.length &&
			this.vertices.some((vertex, i, vs) =>
				face.vertices.every((v2, j) => v2.like(vs[(i + j) % vs.length]))
		)
	},
	containsPoint: function (point) {
		// TODO: problem if the ray hits a corner exactly
		assert(point instanceof V3, "point was"+ point.toString())
		assert(this.plane.containsPoint(point))
		var vs = this.vertices
		var direction = vs[0].minus(vs[vs.length - 1])
		//var testRay = new Line3d(p, direction)
		var i = vs.length - 1;
		var inside = false;
		while (i--) {
			var a = vs[i], b = vs[i + 1]
			// check if segment ab intersects testRay
			if (segmentIntersectsRay(a, b, point, direction)) {
				console.log("ab", a.toString(), b.toString())
				inside = !inside
			}
		}
		return inside
	},
	intersectsLine: function (line) {
		var lambda = line.intersectWithPlaneLambda(this.plane)
		var inside = this.containsPoint(line.at(lambda));
		return inside ? lambda : NaN;
	},
	transform: function (m4) {
		return new BREP.Face(m4.transformedPoints(this.vertices), this.plane.transform(m4))
	}
}
CSG.addTransformationMethodsToPrototype(V3.prototype)
CSG.addTransformationMethodsToPrototype(BREP.Face.prototype)
BREP.prototype.ss = function() {
	return "new BREP(null, null, [" + this.faces.map((face) =>
		"new BREP.Face(["+face.vertices.map((v) => v.toString(x => x)).join(",")+"], "+face.plane.toString((x) => x)+")"
	).join(",") + "])"
}
BREP.prototype.equals = function (brep) {
	return this.faces.length == brep.faces.length &&
			this.faces.every((face) => brep.faces.some((face2) => face.equals(face2)))
}
BREP.prototype.toMesh = function() {
	var vertexIndexMap = new Map()
	var canonEdges = new Set()
	function canonEdge(i0, i1) {
		var iMin = min(i0, i1), iMax = max(i0, i1)
		return (iMin << 16) | iMax
	}
	function uncanonEdge(key) {
		return [key >> 16, key & 0xffff]
	}
	function vertexIndex (vertex) {
		var index = vertexIndexMap.get(vertex)
		if (index === undefined) {
			index = mesh.vertices.length
			vertexIndexMap.set(vertex, index)
			mesh.vertices.push(vertex.toArray())
		}
		return index
	}
	faces = [];
	var mesh = new GL.Mesh({ normals: true, colors: true, lines:true });
	var allValidEdges = [];
	this.faces.forEach((face) => {
		var triangleSubIndexes = triangulateFace(face)
		for (var i = 0; i < triangleSubIndexes.length; i += 3) {
			var
				vi0 = vertexIndex(face.vertices[triangleSubIndexes[i]]),
				vi1 = vertexIndex(face.vertices[triangleSubIndexes[i + 1]]),
				vi2 = vertexIndex(face.vertices[triangleSubIndexes[i + 2]]);
			var dot = V3.normalOnPoints(
				face.vertices[triangleSubIndexes[i]],
				face.vertices[triangleSubIndexes[i + 1]],
				face.vertices[triangleSubIndexes[i + 2]]).dot(face.plane.normal)
			mesh.triangles.push(dot > 0 ? [vi0, vi1, vi2] : [vi0, vi2, vi1])
		}
		for (var i = 0; i < face.vertices.length; i++) {
			var i0 = vertexIndex(face.vertices[i]), i1 = vertexIndex(face.vertices[(i + 1) % face.vertices.length]);
			canonEdges.add(canonEdge(i0, i1))
		}
	})
	var it =canonEdges.values(), ce
	while (ce = it.next().value) mesh.lines.push(uncanonEdge(ce))
	//console.log("mesh.lines", mesh.lines)
	mesh.compile()
	return mesh;
}
function segmentIntersectsRay(a, b, p, dir) {
	var st = segmentIntersectsLineST(a, b, p, dir), s = st && st.segmentAt, t = st && st.lineAt
	console.log("st", st)
	return (st && s > 0 && s < 1 && t > 0)
}
function segmentIntersectsLineST(a, b, p, dir) {
	var ab = b.minus(a)
	var abXdir = ab.cross(dir)
	var div = abXdir.lengthSquared()
	if (NLA.isZero(div)) { return null }
	var anchorDiff = p.minus(a)
	var s = anchorDiff.cross(dir).dot(abXdir) / div
	var t = anchorDiff.cross(ab).dot(abXdir) / div

	//console.log(segmentIntersectsRay, a, b, "ab", ab, "p", p, "dir", dir, s > 0 && t / div >= 0 && t / div <= 1, "s", s, "t", t, "div", div)
	return {segmentAt: s, lineAt: t}
}
function segmentSegmentIntersection(a, b, c, d) {
	var ab = b.minus(a)
	var abXdir = ab.cross(dir)
	var div = abXdir.lengthSquared()
	var anchorDiff = p.minus(a)
	var s = anchorDiff.cross(dir).dot(abXdir) / div
	var t = anchorDiff.cross(ab).dot(abXdir) / div
}
BREP.extrude = function(baseVertices, baseFacePlane, offset) {
	assert (offset instanceof V3, "offset must be V3")
	var topVertices = baseVertices.map((v) => v.plus(offset)).reverse()
	//var topPlane = basePlane.translated(offset)
	var faces = [
		new BREP.Face(baseVertices, baseFacePlane),
		new BREP.Face(topVertices, baseFacePlane.flipped().translated(offset))]
	var m = baseVertices.length
	for (var i = 0; i < m; i++) {
		faces.push(new BREP.Face([
			baseVertices[(i + 1) % m],
			baseVertices[i],
			topVertices[m - 1 - i],
			topVertices[m - 1 - (i + 1) % m]]))
	}
	return new BREP(null, null, faces)
}
BREP.prototype.union = function (brep2) {
}
function arrayFilterMap(arr, f) {
	var result = []
	for (var i = 0; i < arr.length; i++) {
		var v = f(arr[i])
		if (undefined !== v) {
			result.push(v)
		}
	}
	return result
}
function outsideVector(edgeArray, normal) {
	return edgeArray[1].minus(edgeArray[0]).cross(normal)
}
BREP.prototype.clipped = function (brep2) {
	function facesWithEdge(edge, faces) {
		return arrayFilterMap(faces, (face) => {
			for (var i = 0; i < face.vertices.length; i++) {
				var v0 = face.vertices[i], v1 = face.vertices[(i + 1) % face.vertices.length]
				if (v0.like(edge[0]) && v1.like(edge[1]) || v0.like(edge[1]) && v1.like(edge[0])) {
					console.log("OV", outsideVector(edge, face.plane.normal))
					return {face: face, edgeVector: segmentVector(edge), ov: outsideVector(edge, face.plane.normal), reversed: !v0.like(edge[0])}
				}
			}
		})
	}
	function splitsAreaEnclosingVolume(segment, plane) {
		var ab1 = segmentVector(segment).normalized()
		var ov = ab1.cross(plane.normal)
		NLA.assert(ov.hasLength(1))
		var m = M4.forSys(ov, plane.normal.negated(), ab1).transposed()
		var relFaces = facesWithEdge(segment, brep2.faces)
		console.log("relFaces",relFaces)
		var ff = relFaces.map((face) => ({face: face, angle: m.transformVector(face.ov).angleXY()})).sort((a, b) => a.angle - b.angle)
		assert(ff.length != 0)

		//console.log("fffff", ff[0].face.reversed)
		return !ff[0].face.reversed

	}
	function getFacePlaneIntersectionSs(face2, line, facePlane) {
		var xss = []
		face2.vertices.forEach((v0, i, vertices) => {
			var v1 = vertices[(i + 1) % vertices.length], v2 = vertices[(i + 2) % vertices.length]
			console.log("segment", v0.toString(), v1.toString())
			var st = segmentIntersectsLineST(v0, v1, line.anchor, line.dir1), s = st && st.segmentAt, t = st && st.lineAt
			// st is undefined if this segment is parallel to the plane
			if (st) {
				if (NLA.equals(s, 1)) {
					// the end of this segment lies on the line
					// the start doesn't as st is undefined if the segment is parallel to the line
					// we add this point if the next segment ends on the other side of the line
					var pointDistanceSigned = facePlane.distancePointSigned(v2)
					if (NLA.isZero(pointDistanceSigned)) {
						// next segment is colinear
						var outVector = outsideVector([v1, v2], face2.plane.normal)
						if (outVector.dot(segmentVector([v0, v1])) < 0) {
							xss.push({t: line.pointLambda(v1), on: [v0, v1], p: v1, endpoint: true})
						}
					} else if (pointDistanceSigned * facePlane.distancePointSigned(v0) < 0) {
						console.log("HHHEELEO", v1.toString(), line.pointLambda(v1))
						xss.push({t: line.pointLambda(v1), on: [v0, v1], p: v1, endpoint: true})
					}
				}
				if (s - NLA.PRECISION < 0 || s + NLA.PRECISION > 1) { return }
				console.log("HHHEELE2O", v1.toString(), line.pointLambda(v1))
				xss.push({t: t, on: [v0, v1], p: line.at(t)})
			} else {
				// check if the segment is IN the plane
				if (line.containsPoint(v0)) {
					console.log("getFacePLaneINter", line.toString(), v0.toString())
					//
					if (splitsAreaEnclosingVolume([v0, v1], facePlane)) {
						// TODO: ??
						xss.push({t: line.pointLambda(v0), on: [v0, v1], p: v0, endpoint: true})
						xss.push({t: line.pointLambda(v1), on: [v0, v1], p: v1, endpoint: true})
					}
					// figure out if we need to re-add the end point of the colinear segment
					var outVector = outsideVector([v0, v1], face2.plane.normal)
					if (outVector.dot(segmentVector([v1, v2])) > 0) {
						xss.push({t: line.pointLambda(v1), on: [v0, v1], p: v1, endpoint: true})
					}

				}
			}
			console.log("partxss", xss)
		})
		return xss
	}
	function segmentVector(edgeArray) {
		return edgeArray[1].minus(edgeArray[0])
	}
	function getNextPoint(segment, currentPoint, intersections) {
		var sv = segmentVector(segment)
		var sorted = intersections.filter((intersection) => {
				return equalsEdge(segment, intersection.faceSegment) || equalsEdge(segment, intersection.projectedSegment)
			})
			.map((intersection) => ({is: intersection, t: intersection.p.minus(currentPoint).dot(sv)}))
			.filter((pair) => pair.t > 0)
			.sort((pair, pair2) => (pair.t - pair2.t))
		if (sorted.length == 0 ) {
			return segment[1]
		} else {
			return sorted[0].is
		}

	}
	function segmentColinearToEdge(segment, face) {
		var segmentLine = L3.throughPoints(segment[0], segment[1])
		return face.vertices.some((v0, i, vertices) => {
			var v1 = vertices[(i + 1) % vertices.length]
			return null == segmentIntersectsLineST(v0, v1, segmentLine.anchor, segmentLine.dir1) && segmentLine.containsPoint(v0)
		})
	}
	function equalsEdge(arr1, arr2) {
		return arr1[0].like(arr2[0]) && arr1[1].like(arr2[1])
	}
	var newFaces = []
	var edgePoints = []
	var newSegments = []
	this.faces.forEach((face) => {
		// for each face in this BREP
		var plane = face.plane
		console.log("face", face.toString())
		var projectedSegments = []
		var intersections = []
		// iterate though the segments of the projection of face2 on face's plane
		brep2.faces.forEach((face2) => {
			console.log("face2", face2.toString())
			var p2 = face2.plane
			var line = L3.fromPlanes(plane, p2)
			console.log("line", line.toString())
			var xss = getFacePlaneIntersectionSs(face2, line, plane)
			xss.sort((a, b) => a.t - b.t)
			console.log("xss", xss)
			var segmentOutsideVector = line.dir1.cross(plane.normal)
			var dot = segmentOutsideVector.dot(p2.normal)
			for (var i = 0; i < xss.length; i += 2) {
				var a = xss[i], b = xss[i + 1]
				console.log("DOT", dot, "a", a, "b", b)
				var newSegment = dot < 0 ? [a.p, b.p] : [b.p, a.p]
				projectedSegments.push(newSegment)
				// check if segment intersects an edge of face
				face.vertices.forEach((v0, i, vertices) => {
					var v1 = vertices[(i + 1) % vertices.length]
					//console.log("segment", v0.toString(), v1.toString())
					var st = segmentIntersectsLineST(v0, v1, line.anchor, line.dir1), s = st && st.segmentAt, t = st && st.lineAt
					if (st && s > 0 && s < 1) {
						if (NLA.equals(t, a.t)) {
							intersections.push({faceSegment: [v0, v1], projectedSegment: newSegment, p: a.p})
						} else if (NLA.equals(t, b.t)) {
							intersections.push({faceSegment: [v0, v1], projectedSegment: newSegment, p: b.p})
						} else if (t > a.t && t < b.t) {
							intersections.push({faceSegment: [v0, v1], projectedSegment: newSegment, p: line.at(t)})
						}
						if (t > a.t && t <= b.t) {
							/*
							if (NLA.equals(0, s)) {
								intersections.push({faceSegment: [v0, v1], projectedSegment: newSegment, p: v0})
							} else if (NLA.equals(1, s)) {
								intersections.push({faceSegment: [v0, v1], projectedSegment: newSegment, p: v1})
							} else if (s > 0 && s < 1) {
								intersections.push({faceSegment: [v0, v1], projectedSegment: newSegment, p: line.at(t)})
							}
							*/
						}
					} else {
						if (!line.containsPoint(v0)) {
						}
					}
				})
			}
		})
		console.log("face.vertices", face.vertices)
		console.log("projectedSegments\n", projectedSegments.map((segment) => segment[0].toString() + "-"+segment[1].toString()).join("\n"))
		console.log("intersections", intersections)

		if (projectedSegments.length == 0 || intersections.length == 0) {
			newFaceVertices = face.vertices.slice()
			console.log("newFaceVertices", newFaceVertices)
			newFaces.push(new BREP.Face(newFaceVertices, plane))
		}
		while (projectedSegments.some((segment) => !segment.visited) || intersections.some((intersection) => !intersection.visited)) {
			console.log("new face")
			var newFaceVertices = [], i = 0
			var currentpoint = intersections.find((is) => !is.visited), startpoint = currentpoint
			if (!currentpoint) {
				var unvisitedSegment = projectedSegments.find((is) => !is.visited)
				if (!unvisitedSegment) {
					break;
				} else {
					var p = unvisitedSegment[0]
					unvisitedSegment.visited = true
					console.log("unvisitedSegment", unvisitedSegment)
					if (face.containsPoint(p) && !segmentColinearToEdge(unvisitedSegment, face)) {
						throw new Error("implement holes"+ face.toString() + p.toString())
					} else {
						continue
					}
				}
			}
			// calculate a - b
			//
			do {
				console.log("currentpoint", currentpoint)
				if (currentpoint.faceSegment) {
					intersections.map((is) => {if (is.p.like(currentpoint.p)) is.visited = true })
					// at an intersection
					var faceEdgeOutsideVector = outsideVector(currentpoint.faceSegment, plane.normal)
					var projectedSegmentVector = segmentVector(currentpoint.projectedSegment)
					newFaceVertices.push(currentpoint.p)
					if (faceEdgeOutsideVector.dot(projectedSegmentVector) < 0) {
						// projectedSegmentVector points into the edge
						// as we want to subtract the projectedSegment s, follow that
						projectedSegments.map((ps) => {
							if (equalsEdge(ps, currentpoint.projectedSegment))  { console.log("ps",ps); ps.visited = true} })
						currentpoint = getNextPoint(currentpoint.projectedSegment, currentpoint.p, intersections)
						console.log("nxt currentpoint", currentpoint)
					} else {
						// follow segmentEdge
						currentpoint = getNextPoint(currentpoint.faceSegment, currentpoint.p, intersections)
						console.log("nxt currentpoint", currentpoint)

					}
				} else {
					// simple edge point, find next edge
					newFaceVertices.push(currentpoint)
					var currentPointFaceEdgeIndex = face.vertices.findIndex((val) => val.like(currentpoint)), nextEdgeEnd
					if (currentPointFaceEdgeIndex != -1) {
						nextEdgeEnd = face.vertices[(currentPointFaceEdgeIndex + 1) % face.vertices.length]
					} else {
						projectedSegments.map((ps) => {
							if (ps[0].like(currentpoint) || ps[1].like(currentpoint)) { console.log("ps", ps); ps.visited = true} })
						nextEdgeEnd = projectedSegments.find((val) => {
							return val[0].like(currentpoint)
						})[1]
					}
					currentpoint = getNextPoint([currentpoint, nextEdgeEnd], currentpoint, intersections)
				}
			} while (currentpoint != startpoint && i++ < 20)
			newFaces.push(new BREP.Face(newFaceVertices, plane))
		}
	})
	return new BREP(null, null, newFaces)
}
/*
 console.log("fff", xss2.map((t) => t.t).join("->"))
 // generate overlapping segments
 var i = 0, j = 0
 while (i < xss.length && j < xss2.length) {
 var a, b, c, d
 if (xss2[j].t > xss[i].t && xss2[j].t < xss[i + 1].t) {
 a = xss[i], b = xss[i + 1], c = xss2[j], d = xss2[j + 1]
 console.log("asdjlj")
 } else if (xss[i].t > xss2[j].t && xss[i].t < xss2[j + 1].t) {
 a = xss2[j], b = xss2[j + 1], c = xss[i], d = xss[i + 1]
 } else {
 // move the lowest one
 if (xss[i].t < xss2[j].t) {
 i += 2
 } else {
 j += 2
 }
 continue
 }
 if ((d.t < b.t) ^ (a == xss2[j])) {
 j += 2
 } else {
 i += 2
 }
 console.log("i", i, "j", j)
 console.log("a", a, "b", c, "c", d, "d")
 // a---------b
 //     c-----?
 // first segment in xss2 starts in the middle of the first segment of xss
 var segmentInsideVector = line.dir1.cross(p2.normal)
 var dot = segmentInsideVector.dot(p.normal)
 edgePoints.push(c)
 if (d.t < b.t) {
 // a---------b
 //     c--d
 newSegments.push({face: face, segment: dot > 0 ? [c.p, d.p] : [d.p, c.p]})
 newSegments.push({face: face2, segment: dot < 0 ? [c.p, d.p] : [d.p, c.p]})
 edgePoints.push(d)
 j += 2
 } else {
 // a---------b
 //     c--------d
 // mutual overlap
 edgePoints.push(b)
 newSegments.push({face: face, segment: dot > 0 ? [c.p, b.p] : [d.p, b.p]})
 newSegments.push({face: face2, segment: dot < 0 ? [c.p, b.p] : [d.p, b.p]})
 i += 2
 }
 var xss2 = getFacePlaneIntersectionSs(face, line)
 console.log("fff", xss.map((t) => t.t).join("->"))
 console.log("newSegments", newSegments)
 console.log("edgePoints", edgePoints)
 */


var newFaces = function () {

}
function triangulateFace(face) {
	var coords = [['y', 'z'], ['z', 'x'], ['x', 'y']][face.plane.normal.absMaxDim()]
	var contour = [].concat.apply([], face.vertices.map((vertex) => [vertex[coords[0]], vertex[coords[1]]]))
	return earcut(contour, [])
}
function isCCW(vertices, normal) {
	assert(!normal.isZero())
	var maxDim = normal.absMaxDim()
	// order is important, coord0 and coord1 must be set so that coord0, coord1 and maxDim span a right-hand coordinate system
	var [coord0, coord1] = [['y', 'z'], ['z', 'x'], ['x', 'y']][maxDim]
	var doubleSignedArea = vertices.map((v0, i, vertices) => {
		var v1 = vertices[(i + 1) % vertices.length]
		return (v1[coord0] - v0[coord0]) * (v1[coord1] + v0[coord1])
	}).reduce((a, b) => a + b)
	console.log("doubleSignedArea", doubleSignedArea)
	assert(!NLA.isZero(doubleSignedArea))
	return doubleSignedArea * Math.sign(normal.e(maxDim)) < 0
}

var zoomFactor = 10;
var eyePos = V3(0, 0, 1000), eyeFocus = V3.ZERO, eyeUp = V3.Y;
function setupCamera() {
	gl.matrixMode(gl.PROJECTION);
	gl.loadIdentity();
	//gl.perspective(70, gl.canvas.width / gl.canvas.height, 0.1, 1000);
	var lr = gl.canvas.width / 2 / zoomFactor;
	var bt = gl.canvas.height / 2 / zoomFactor;
	gl.ortho(-lr, lr, -bt, bt, -1e3, 1e5);
	gl.lookAt(
		eyePos.x, eyePos.y, eyePos.z,
		eyeFocus.x, eyeFocus.y, eyeFocus.z,
		eyeUp.x, eyeUp.y, eyeUp.z);
	gl.matrixMode(gl.MODELVIEW);
}
var v2as = (vertex) => "["+vertex.x+", "+vertex.y+", "+vertex.z+"]";
function CustomPlane(anchor2, right, up, upStart, upEnd, rightStart, rightEnd, color) {
	var p = P3.throughPoints(anchor2, right, up, CustomPlane.prototype);
	p.up = up;
	p.right = right;
	p.upStart = upStart;
	p.upEnd = upEnd;
	p.rightStart = rightStart;
	p.rightEnd = rightEnd;
	p.color = color;
	p.id = lineId++;
	return p;
}

function rgbToArr4(color) {
	return [(color >> 16) / 255.0, ((color >> 8) & 0xff) / 255.0, (color & 0xff) / 255.0, 1.0];
}
var lineId = 0
CustomPlane.prototype = Object.create(P3.prototype);
CustomPlane.prototype.constructor = CustomPlane;
CustomPlane.prototype.toString = function () {
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
				|| line2Index >= 2 && this.rightStart <= info.t && info.t <= this.rightEnd) {
				return info.distance;
			} else {
				return 1e9;
			}
		}, this).min();
	} catch (e) {
		console.log(e);
	}
}

//var sketchPlane = new CustomPlane(V3.X, V3(1, 0, -1).unit(), V3.Y, -500, 500, -500, 500, 0xff00ff);
var planes = [
	CustomPlane(V3.ZERO, V3.Y, V3.Z, -500, 500, -500, 500, 0xff0000),
	CustomPlane(V3.ZERO, V3.X, V3.Z, -500, 500, -500, 500, 0x00ff00),
	CustomPlane(V3.ZERO, V3.X, V3.Y, -500, 500, -500, 500, 0x0000ff),
	//	sketchPlane
];

function drawPlanes() {
	planes.forEach(function (plane) {
		gl.pushMatrix();
		gl.multMatrix(M4.forSys(plane.right, plane.up, plane.normal));
		gl.translate(plane.rightStart, plane.upStart, 0);
		gl.scale(plane.rightEnd - plane.rightStart, plane.upEnd - plane.upStart, 1);

		var shader = singleColorShader;

		shader.uniforms({
			color: rgbToArr4(plane.color)
		}).draw(xyLinePlaneMesh, gl.LINES);

		gl.popMatrix();
	});
}
// abcd can be in any order
BREP.tetrahedron = function (a, b, c, d) {
	var dDistance = P3.throughPoints(a, b, c).distancePointSigned(d)
	if (NLA.isZero(dDistance)) {
		throw new Error("four points are coplanar")
	}
	var faceVertices = [
		[a, b, c], [a, d, b], [b, d, c], [c, d, a]
	]
	if (dDistance > 0) {
		faceVertices.forEach((vertices) => vertices.reverse())
	}
	return new BREP(null, null, faceVertices.map((vertices) => new BREP.Face(vertices)))
}
var singleColorShader, textureColorShader, singleColorShaderHighlight, arcShader, arcShader2,xyLinePlaneMesh;

window.loadup = function () {
	window.onerror = function (errorMsg, url, lineNumber, column, errorObj) {
		console.log(errorMsg, url, lineNumber, column, errorObj);
	}
	gl = GL.create({canvas: document.getElementById("testcanvas")});
	gl.fullscreen();
	gl.canvas.oncontextmenu = () => false;

	gl.scaleVector = function (x, y, z) {
		gl.multMatrix(M4.scaleVector(x, y, z));
	};
	gl.translateVector = gl.translateV3
	gl.scaleVector = gl.scaleV3

	var a= V3.Y
	for (var i = 0; i < PI * 2; i+=PI * 2 / 16) {
		var x = V3.Y.rotateX(i)
		console.log("Math", i, x, a.cross(x).length(), a.dot(x), round(rad2deg(Math.atan2(a.cross(x).length(), a.dot(x)))))
	}

	setupCamera();
	//gl.cullFace(gl.FRONT_AND_BACK);
	gl.clearColor(1.0, 1.0, 1.0, 0.0);
	gl.enable(gl.BLEND);
	gl.enable(gl.DEPTH_TEST);
	gl.enable(gl.CULL_FACE);
	gl.depthFunc(gl.LEQUAL)
	gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA); // TODO ?!

	cubeMesh = GL.Mesh.cube();
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
	gl.loadIdentity();
	gl.scale(10, 10, 10);

	gl.loadIdentity();

	gl.onmousemove = function (e) {
		if (e.dragging) {
			if (e.buttons & 2) {
				// pan
				var moveCamera = V3(-e.deltaX * 2 / gl.canvas.width, e.deltaY * 2 / gl.canvas.height, 0);
				var inverseProjectionMatrix = gl.projectionMatrix.inversed();
				var worldMoveCamera = inverseProjectionMatrix.transformVector(moveCamera);
				eyePos = eyePos.plus(worldMoveCamera);
				eyeFocus = eyeFocus.plus(worldMoveCamera);
				setupCamera();
				paintScreen();
			}
			if (e.buttons & 4) {
				var rotateLR = -e.deltaX / 6.0;
				var rotateUD = -e.deltaY / 6.0;
				// rotate
				var matrix = M4.rotation(rotateLR, eyeUp);
				var horizontalRotationAxis = eyeFocus.minus(eyePos).cross(eyeUp);
				matrix = matrix.times(M4.rotation(rotateUD, horizontalRotationAxis));
				//console.log(matrix);
				eyePos = matrix.transformVector(eyePos);
				eyeUp = matrix.transformVector(eyeUp)
				//console.log("eyepos", eyePos);
				setupCamera();
				paintScreen();
			}
		}
		/*
		 } catch (e) {
		 console.log(e);
		 }*/
	}
	xyLinePlaneMesh = new GL.Mesh({lines: true, triangles: false});
	xyLinePlaneMesh.vertices = [[0, 0], [0, 1], [1, 1], [1, 0]];
	xyLinePlaneMesh.lines = [[0, 1], [1, 2], [2, 3], [3, 0]];
	xyLinePlaneMesh.compile();
	singleColorShader = new GL.Shader($('vertexShaderBasic').text, $('fragmentShaderColor').text);
	singleColorShaderHighlight = new GL.Shader($('vertexShaderBasic').text, $('fragmentShaderColorHighlight').text);
	textureColorShader = new GL.Shader($('vertexShaderTextureColor').text, $('fragmentShaderTextureColor').text);
	arcShader = new GL.Shader($('vertexShaderRing').text, $('fragmentShaderColor').text);
	arcShader2 = new GL.Shader($('vertexShaderArc').text, $('fragmentShaderColor').text);

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

	f0 = new BREP.Face([V3(0, 0, 0), V3(10, 0, 0), V3(10, 10, 0), V3(0, 10, 0)], P3.XY)
	f4 = new BREP.Face([V3(5, 0, 5), V3(5, 5, 5), V3(0, 5, 5), V3(0, 0, 5)], P3(V3.Z, 5))
	var boxBottom = new BREP.Face([V3(5, 0, 0), V3(5, 5, 0), V3(0, 5, 0), V3(0, 0, 0)].reverse(), P3(V3.Z.negated(), 0))
	f00 = new BREP.Face([V3(0, 0, -6), V3(10, 0, -6), V3(11, 10, 0), V3(0, 10, 0)])
	f01 = new BREP.Face([V3(10, 0, 0), V3(10, 0, -6), V3(11, 10, 0)])
	f1 = new BREP.Face([V3(5, 8, -5), V3(9, 12, 1), V3(5, 8, 5)])
	f2 = new BREP.Face([V3(5, 8, -5), V3(5, 8, 5), V3(1, 12, 2)])
	f3 = new BREP.Face([V3(5, 8, -5), V3(1, 12, 2), V3(9, 12, 1)])
	var a =new BREP(null, null, [new BREP.Face([V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0),V3(0, 10, 0)], P3(V3(0, 0, 1), 0))])
	var b = BREP.tetrahedron(V3(11, 11, 0), V3(7, 7, 2), V3(7, 9, 2), V3(7, 7, -5))
	//var a = new BREP(null, null, [ boxBottom])
	//var a = new BREP(null, null, [new BREP.Face([V3(5, 0, 5), V3(0, 0, 5), V3(0, 5, 5), V3(5, 5, 5)], P3(V3(0, 0, 1), 5))])
	//var a = BREP.box(5, 5, 5)
	//var tet = BREP.tetrahedron(V3(4, 4, 4), V3(8, 6, 0), V3(6, 7, 5), V3(1, 7, 5))
	//var b = BREP.tetrahedron(V3(2, -2, -2), V3(2, 4, -2), V3(1, -2, 4), V3(3, 4, 4))
	var concave = BREP.extrude(
		[V3(5,5,0), V3(5,8,0),V3(8,10,-2),V3(2,11,2)].reverse(),
		P3.normalOnAnchor(V3(-2, 0, -3), V3(5, 5, 0)),
		V3(4,0,6))
	//var b = new BREP(null, null, [f1, f2, f3])
	//var b = tet
	aMesh = a.toMesh()
	bMesh = b.toMesh()
	console.log("css",b.ss())


	try {
		var c = a.clipped(b)
		cMesh = c.toMesh()
		console.log("css",c.ss())
	} catch (e) {
		console.log(e)
	}
	paintScreen()
}
var aMesh, bMesh, cMesh
function paintScreen() {
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
	gl.loadIdentity();
	gl.scale(10, 10, 10);

	gl.loadIdentity();


	if (aMesh) {
		singleColorShader.uniforms({
			color: rgbToArr4(0x9edbf9)
		}).draw(aMesh);
		singleColorShader.uniforms({
			color: rgbToArr4(0x77b0e0)
		}).draw(aMesh, gl.LINES);
	}
	if (bMesh) {
		singleColorShader.uniforms({
			color: rgbToArr4(0xEE4144)
		}).draw(bMesh);
		singleColorShader.uniforms({
			color: rgbToArr4(0x1E98D3)
		}).draw(bMesh, gl.LINES);
	}
	if (cMesh) {
		singleColorShader.uniforms({
			color: rgbToArr4(0xF37033)
		}).draw(cMesh);
		singleColorShader.uniforms({
			color: rgbToArr4(0x1E98D3)
		}).draw(cMesh, gl.LINES);
	}
	drawPlanes();
}