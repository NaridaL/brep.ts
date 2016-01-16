"use strict";
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
function BREP(vertices, edges, faces, infiniteVolume) {
	this.vertices = vertices;
	this.edges = edges;
	this.faces = faces;
	this.infiniteVolume = infiniteVolume || false
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
BREP.Face = function(vertices, plane, name) {
	assert(!plane || plane instanceof P3, "plane is not a P3")
	this.vertices = vertices
	this.name = name
	// this.aabb
	this.plane = plane ? plane : P3.throughPoints(vertices[0], vertices[1], vertices[2])
	//console.log(vertices, plane)
	//console.log(plane && vertices.map(v=> plane.distanceToPointSigned(v)))
	assert(!plane || vertices.every((v) => plane.containsPoint(v)), "all vertices must be contained within the plane")
	assert(Array.from(BREP.Face.segmentsIterator(vertices)).every((seg) => !seg[0].like(seg[1])), "no degenerate segments"+vertices.map(v=>v.ss).join(","))
	var ccw = isCCW(vertices, this.plane.normal)
//	assert(!plane || ccw, "vertices are not ccw to passed plane")
	if (!ccw) {
		this.plane = this.plane.flipped()
	}
	this.id = faceId++
}
BREP.Face.segmentsIterator = function*(vertices) {
	for (var i = 0; i < vertices.length; i++) {
		yield [vertices[i], vertices[(i + 1) % vertices.length]]
	}
}
BREP.prototype.flipped = function () {
	return new BREP(null, null, this.faces.map(f => f.flipped()), !this.infiniteVolume)
}
BREP.prototype.transform = function (m4) {
	return new BREP(null, null, this.faces.map(f => f.transform(m4)), this.infiniteVolume)
}
BREP.prototype.minus = function (brep) {
	return this.clipped(brep, true, true).newAll(brep.flipped().clipped(this.flipped(), true, false))
}
BREP.prototype.plus = function (brep) {
	return this.clipped(brep, true, true).newAll(brep.clipped(this, false, true))
}
BREP.segmentIntersectsVertexLoop = function (a, b, vertices) {
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
		assert(this.plane.containsPoint(point), "this.plane.containsPoint(point)" + this.plane.distanceToPointSigned(point))
		var vs = this.vertices
		var direction = vs[0].minus(vs[vs.length - 1])
		//var testRay = new Line3d(p, direction)
		var i = vs.length - 1;
		var inside = false;
		while (i--) {
			var a = vs[i], b = vs[i + 1]
			// check if segment ab intersects testRay
			if (segmentIntersectsRay(a, b, point, direction)) {
				inside = !inside
			}
		}
		return inside
	},
	containsPoint2: function (point, includeBorder) {
		assert(point instanceof V3, "point was"+ point.toString())
		assert(this.plane.containsPoint(point))
		var vs = this.vertices.map(v => v.minus(point))
		if (vs[0].isZero()) return includeBorder // point lies on first vertex
		//var testRay = new Line3d(p, direction)
		var angle = 0
		for (var i = 0; i < vs.length; i++) {
			var av = vs[i], bv = vs[(i + 1) % vs.length]
			if (bv.isZero()) return includeBorder // point lies on vertex
			var delta = av.angleRelativeNormal(bv, this.plane.normal);
			if (NLA.eqAngle(Math.PI, delta)) return includeBorder // point lies on edge
			angle += delta
		}
		return !NLA.isZero(angle)
	},
	pointsToInside: function (p, dir, includeEdges) {
		var ff = this.vertices.map((v2, i, vertices) => {
			var vNext = vertices[(i + 1) % vertices.length]
			var vPrev = vertices[(i - 1 + vertices.length) % vertices.length]
			if (v2 == p) {
				console.log("sadas", dir.angleRelativeNormal(vNext.minus(v2), this.plane.normal).toSource(),
					vNext.minus(v2).toSource())
				return [
					{reversed: false, angle: (dir.angleRelativeNormal(vNext.minus(v2), this.plane.normal) + 2 * Math.PI + NLA.PRECISION / 2) % (2 * Math.PI)},
					{reversed: true, angle: (dir.angleRelativeNormal(vPrev.minus(v2), this.plane.normal) + 2 * Math.PI + NLA.PRECISION / 2) % (2 * Math.PI)}]
			} else {
				return []
			}
		}).concatenated().sort((a, b) => a.angle - b.angle)
		assert(ff.length != 0, "ff.length != 0")
		return NLA.isZero(ff[0].angle) ? includeEdges : ff[0].reversed
	},
	flipped: function () {
		return new BREP.Face(this.vertices.slice().reverse(), this.plane.flipped(), this.name + "flipped")
	},
	intersectsLine: function (line) {
		assert(line instanceof L3)
		var lambda = line.intersectWithPlaneLambda(this.plane)
		if (!Number.isFinite(lambda)) {
			return NaN
		}
		var inside = this.containsPoint(line.at(lambda))
		return inside ? lambda : NaN
	},
	withHole: function (holeVertices) {
		var v = this.vertices
		var fail, i, j
		for (i = 0; i < v.length; i++) {
			var a = v[i]
			for (j = 0; j < holeVertices.length; j++) {
				var b = holeVertices[j]
				fail = false
				// check if segment ab intersects any other segment
				for (var k = i + 1; k < v.length + i - 1; k++) {
					var v0 = v[k % v.length], v1 = v[(k + 1) % v.length]
					if (segmentsTouchOrIntersect(a, b, v0, v1)) {
						fail = true;
						break;
					}
				}
				if (fail) continue
				for (var k = j + 1; k < holeVertices.length + j - 1; k++) {
					var v0 = holeVertices[k % holeVertices.length], v1 = holeVertices[(k + 1) % holeVertices.length]
					if (segmentsTouchOrIntersect(a, b, v0, v1)) {
						fail = true;
						break;
					}
				}
				if (!fail) break
			}
			if (!fail) break
		}
		if (fail) {
			throw new Error("")
		}
		var cycledHole = holeVertices.slice(j).concat(holeVertices.slice(0, j + 1))
		var newVertices = v.slice(0, i + 1).concat(cycledHole, v.slice(i))
		console.log(newVertices.toSource())
		return new BREP.Face(newVertices, this.plane, this.name + "holed")
	},
	transform: function (m4) {
		return new BREP.Face(m4.transformedPoints(this.vertices), this.plane.transform(m4), this.name)
	}
}
CSG.addTransformationMethodsToPrototype(V3.prototype)
CSG.addTransformationMethodsToPrototype(BREP.Face.prototype)
CSG.addTransformationMethodsToPrototype(BREP.prototype)
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
	var mesh = new GL.Mesh({ normals: true, colors: true, lines:true });
	var allValidEdges = [];
	this.faces.forEach((face) => {
		//var faceEdges = new Map()
		var triangleSubIndexes = triangulateFace(face)
		var vs = face.vertices;
		for (var i = 0; i < triangleSubIndexes.length; i += 3) {
			var
				vi0 = vertexIndex(vs[triangleSubIndexes[i]]),
				vi1 = vertexIndex(vs[triangleSubIndexes[i + 1]]),
				vi2 = vertexIndex(vs[triangleSubIndexes[i + 2]]);
			var dot = V3.normalOnPoints(
				vs[triangleSubIndexes[i]],
				vs[triangleSubIndexes[i + 1]],
				vs[triangleSubIndexes[i + 2]]).dot(face.plane.normal)
			mesh.triangles.push(dot > 0 ? [vi0, vi1, vi2] : [vi0, vi2, vi1])
		}
		for (var i = 0; i < vs.length; i++) {
			var i0 = vertexIndex(vs[i]), i1 = vertexIndex(vs[(i + 1) % vs.length]);

			canonEdges.add(canonEdge(i0, i1))
		}
	})
	var it =canonEdges.values(), ce
	while (ce = it.next().value) mesh.lines.push(uncanonEdge(ce))
	//console.log("mesh.lines", mesh.lines)
	mesh.compile()
	return mesh;
}
BREP.prototype.toNormalMesh = function() {
	var vertexIndexMap = new Map()
	var canonEdges = new Set()
	function canonEdge(i0, i1) {
		var iMin = min(i0, i1), iMax = max(i0, i1)
		return (iMin << 16) | iMax
	}
	function uncanonEdge(key) {
		return [key >> 16, key & 0xffff]
	}
	var mesh = new GL.Mesh({ normals: true, colors: true, lines:false });
	var allValidEdges = [];
	var faceIndexes = new Map()
	this.faces.forEach((face) => {
		//var faceEdges = new Map()
		var triangleSubIndexes = triangulateFace(face)
		var vs = face.vertices;
		var startIndex = mesh.vertices.length
		faceIndexes.set(face, {start: mesh.triangles.length * 3, count: triangleSubIndexes.length})
		vs.forEach(v => {
			mesh.vertices.push(v.toArray())
			mesh.normals.push(face.plane.normal.els())
		})
		for (var i = 0; i < triangleSubIndexes.length; i += 3) {
			var
				vi0 = triangleSubIndexes[i],
				vi1 = triangleSubIndexes[i + 1],
				vi2 = triangleSubIndexes[i + 2];
			var dot = V3.normalOnPoints(
				vs[vi0],
				vs[vi1],
				vs[vi2]).dot(face.plane.normal)
			mesh.triangles.push((dot > 0 ? [vi0, vi1, vi2] : [vi0, vi2, vi1]).map(x => x + startIndex))
		}
		/*
		for (var i = 0; i < vs.length; i++) {
			var i0 = vertexIndex(vs[i]), i1 = vertexIndex(vs[(i + 1) % vs.length]);

			canonEdges.add(canonEdge(i0, i1))
		}
		*/
	})
	//console.log("mesh.lines", mesh.lines)
	mesh.faceIndexes = faceIndexes
	mesh.compile()
	return mesh;
}
function segmentIntersectsRay(a, b, p, dir) {
	var st = segmentIntersectsLineST(a, b, p, dir), s = st && st.segmentAt, t = st && st.lineAt
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
function segmentSegmentIntersectionST(a, b, c, d) {
	var ab = b.minus(a)
	var cd = d.minus(c)
	var abXcd = ab.cross(cd)
	var div = abXcd.lengthSquared()
	var ac = c.minus(a)
	var s = ac.cross(cd).dot(abXcd) / div
	var t = ac.cross(ab).dot(abXcd) / div
	return {s: s, t: t}
}
function segmentsTouchOrIntersect(a, b, c, d) {
	var {s, t} = segmentSegmentIntersectionST(a, b, c, d)
	return (NLA.equals(s, 0) || NLA.equals(s, 1) || (s > 0 && s < 1))
		 && (NLA.equals(t, 0) || NLA.equals(t, 1) || (t > 0 && t < 1))
}
BREP.extrude = function(baseVertices, baseFacePlane, offset, name) {
	assert (baseVertices.every(v => v instanceof V3))
	assert (baseFacePlane instanceof P3)
	assert (offset instanceof V3, "offset must be V3")
	if (baseFacePlane.normal.dot(offset) > 0) baseFacePlane = baseFacePlane.flipped()
	if (!isCCW(baseVertices, baseFacePlane.normal)) {
		baseVertices = baseVertices.reverse()
	}
	var topVertices = baseVertices.map((v) => v.plus(offset)).reverse()
	//var topPlane = basePlane.translated(offset)
	var faces = [
		new BREP.Face(baseVertices, baseFacePlane, name + "base"),
		new BREP.Face(topVertices, baseFacePlane.flipped().translated(offset), name + "roof")]
	var m = baseVertices.length
	for (var i = 0; i < m; i++) {
		faces.push(new BREP.Face([
			baseVertices[(i + 1) % m],
			baseVertices[i],
			topVertices[m - 1 - i],
			topVertices[m - 1 - (i + 1) % m]], undefined, name + "wall" + i))
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
BREP.prototype.newAll = function(brep2) {
	return new BREP(null, null, this.faces.concat(brep2.faces), this.infiniteVolume)
}
function facesWithEdge(edge, faces) {
	return arrayFilterMap(faces, (face) => {
		for (var i = 0; i < face.vertices.length; i++) {
			var v0 = face.vertices[i], v1 = face.vertices[(i + 1) % face.vertices.length]
			if (v0.like(edge[0]) && v1.like(edge[1]) || v0.like(edge[1]) && v1.like(edge[0])) {
				console.log("OV", outsideVector([v0, v1], face.plane.normal))
				return {face: face, edgeVector: segmentVector(edge), ov: outsideVector([v0, v1], face.plane.normal), reversed: !v0.like(edge[0])}
			}
		}
	})
}
function segmentVector(edgeArray) {
	return edgeArray[1].minus(edgeArray[0])
}
function splitsVolumeEnclosingFaces(brep, segment, dir, faceNormal, removeCoplanarSame, removeCoplanarOpposite) {
	console.log("splitsVolumeEnclosingFaces", brep, segment, dir, faceNormal)
	var ab1 = segmentVector(segment).normalized()
	var relFaces = facesWithEdge(segment, brep.faces)
	console.log("relFaces", relFaces)
	var ff = relFaces.map((face) => {
		console.log("dir", dir.ss, face.ov.negated().ss, dir.angleRelativeNormal(face.ov.negated(), ab1))
		return ({
			face: face,
			angle: (dir.angleRelativeNormal(face.ov.negated(), ab1) + 2 * Math.PI + NLA.PRECISION / 2) % (2 * Math.PI)
		})
	}).sort((a, b) => a.angle - b.angle)
	assert(ff.length != 0)

	console.log("fffff", segment, ff.toSource())
	if (NLA.isZero(ff[0].angle)) {
		var coplanarSame = ff[0].face.face.plane.normal.dot(faceNormal) > 0
		console.log("returning", coplanarSame
			? assert(removeCoplanarSame !== undefined) && removeCoplanarSame
			: assert(removeCoplanarOpposite !== undefined) && removeCoplanarOpposite
		)
		return coplanarSame
			? assert(removeCoplanarSame !== undefined) && removeCoplanarSame
			: assert(removeCoplanarOpposite !== undefined) && removeCoplanarOpposite
	} else {
		return !ff[0].face.reversed
	}
}
BREP.prototype.clipped = function (brep2, removeCoplanarSame, removeCoplanarOpposite) {
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
					var pointDistanceSigned = facePlane.distanceToPointSigned(v2)
					if (NLA.isZero(pointDistanceSigned)) {
						// next segment is colinear
						var outVector = outsideVector([v1, v2], face2.plane.normal)
						if (outVector.dot(segmentVector([v0, v1])) < 0) {
							xss.push({t: line.pointLambda(v1), on: [v0, v1], p: v1, endpoint: true})
						}
					} else if (pointDistanceSigned * facePlane.distanceToPointSigned(v0) < 0) {
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
					var normal = facePlane.normal
					console.log("brep2.infiniteVolume",brep2.infiniteVolume)
					if (splitsVolumeEnclosingFaces(brep2, [v0, v1], facePlane.projectedVector(face2.plane.normal), normal, removeCoplanarSame, removeCoplanarOpposite)
					!= splitsVolumeEnclosingFaces(brep2, [v0, v1], facePlane.projectedVector(face2.plane.normal).negated(), normal, removeCoplanarSame, removeCoplanarOpposite)) {
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
	function getNextPoint(segment, currentPoint, intersections) {
		var sv = segmentVector(segment)
		var sorted = intersections.filter((intersection) => {
				return equalsEdge(segment, intersection.faceSegment)
					|| equalsEdge(segment, intersection.projectedSegment)
			})
			.map((intersection) => ({is: intersection, t: intersection.p.minus(currentPoint).dot(sv)}))
			.filter((pair) => pair.t > 0)
			.sort((pair, pair2) => (pair.t - pair2.t))
		if (sorted.length == 0 ) {
			return segment[1]
		} else {
			return sorted[0].is.p
		}

	}
	function possDirs(currentPoint, intersections) {
		return intersections.filter((is) => {
			return is.p.like(currentPoint)
		})
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
		var addedFaces = 0
		// for each face in this BREP
		var plane = face.plane
		console.log("face", face.toString())
		var projectedSegments = []
		var intersections = []
		brep2.faces.forEach((face2) => {
			console.log("face2", face2.toString())
			var p2 = face2.plane
			if (plane.isParallelToPlane(p2)) {
				return;
			}
			var line = L3.fromPlanes(plane, p2)
			console.log("line", line.toString())
			var xss = getFacePlaneIntersectionSs(face2, line, plane)
			xss.sort((a, b) => a.t - b.t)
			console.log("xss", xss)
			var segmentOutsideVector = line.dir1.cross(plane.normal)
			var dot = segmentOutsideVector.dot(p2.normal)
			// iterate though the segments of the intersection of face2 on face's plane
			for (var i = 0; i < xss.length; i += 2) {
				var a = xss[i], b = xss[i + 1]
				var newSegment = dot < 0 ? [a.p, b.p] : [b.p, a.p], t0 = dot < 0 ? a.t : b.t, t1 = dot < 0 ? b.t : a.t
				projectedSegments.push(newSegment)
				// check if segment intersects an edge of face
				face.vertices.forEach((v0, i, vertices) => {
					var v1 = vertices[(i + 1) % vertices.length]
					//console.log("segment", v0.toString(), v1.toString())
					var st = segmentIntersectsLineST(v0, v1, line.anchor, line.dir1), s = st && st.segmentAt, t = st && st.lineAt
					if (st) {
						s = NLA.snapTo(s, 0)
						s = NLA.snapTo(s, 1)
						t = NLA.snapTo(t, t0)
						t = NLA.snapTo(t, t1)
						/*
							Invalid combinations:
								s < 0 || s > 1       intersection not on face edge
								t < a.t || t > b.t   intersection not on projected segment
								s == 0 && t == t0    both directed edges start on the same point
								s == 1 && t == t1    both directed edges end on the same point
						 */
						var point
						if (s == 0 && (t > a.t && t < b.t || t == t1)) {
							point = v0
						} else if (s > 0 && s < 1 && t >= a.t && t <= b.t) {
							point = t == a.t ? a.p : (t == b.t ? b.p : line.at(t))
						} else if (s == 1 && (t > a.t && t < b.t || t == t0)) {
							point = v1
						}
						if (point) {
							var newIs = {faceSegment: [v0, v1], projectedSegment: newSegment, p: point,
								projectedEndpoint: t == t1, edgeEndpoint: s == 1, edgeStartpoint: s == 0, projectedStartpoint: t == t0}
							intersections.push(newIs)
							var popped = false
							if (newIs.edgeEndpoint) {
								// check that the projected segment vector splits an area enclosing pair at p == faceSegment[1]
								var dirr = segmentVector(newIs.projectedSegment)
								console.log("face.pointsToInside(point, dirr)", face.pointsToInside(point, dirr, false))
								if(!face.pointsToInside(point, dirr, false)) {
									console.log("popping is")
									intersections.pop()
									popped = true
								}
							} else if (newIs.edgeStartpoint) {
								if (segmentVector(newIs.faceSegment).dot(outsideVector(newIs.projectedSegment, face.plane.normal)) > 0) {
									intersections.pop()
									popped = true
								}
							}
							/*
							if (!popped && newIs.projectedEndpoint) {
								var dirr = segmentVector(newIs.projectedSegment)
								// check that the edge segment vector splits an area enclosing pair at p == projectedSegment[1]
								if(!pointsToInsideSegmentsPoly(point, dirr)) {
									console.log("popping is")
									intersections.pop()
									popped = true
								}
							}*/
						}
					} else {
						if (!line.containsPoint(v0)) {
						}
					}
				})
			}
		})
		function pointsToInside(looseSegments, p, dir, includeEdges) {
			var ff = looseSegments.map((segment) => {
				if (segment[0].like(p) || segment[1].like(p)) {
					return [{
						reversed: segment[0].like(p),
						angle: (dir.angleRelativeNormal(segmentVector(segment), face.plane.normal) + 2 * Math.PI + NLA.PRECISION / 2) % (2 * Math.PI)}]
				} else {
					return []
				}
			}).concatenated().sort((a, b) => a.angle - b.angle)
			assert(ff.length != 0, "ff.length != 0")
			console.log(dir.toSource(), "dirr", dir.toSource(), ff)
		return NLA.isZero(ff[0].angle) ? includeEdges : ff[0].reversed
		}
		// check intersections again, now that we have all the cross-section segments
		intersections = intersections.filter(is => {
			if (is.projectedEndpoint) {
				return pointsToInside(projectedSegments, is.p, segmentVector(is.faceSegment), false)
			}
			return true
		})
		console.log("face.vertices", face.vertices)
		console.log("projectedSegments\n", projectedSegments.map((segment) => segment[0].toString() + "-"+segment[1].toString()).join("\n"))
		if (projectedSegments.length == 0) {
			if (!brep2.infiniteVolume) {
				newFaceVertices = face.vertices.slice()
				newFaces.push(new BREP.Face(newFaceVertices, plane, face.name))
			}
			return
		}

		var loopCount = 0
		while (projectedSegments.some((segment) => !segment.visited) || intersections.some((intersection) => !intersection.visited)) {
			loopCount++
			console.log("new face")
			var newFaceVertices = [], i = 0
			console.log("intersections\n", intersections.map(is => is.toSource()).join("\n"))
			var onProjected = true
			var currentPoint = intersections.find((is) => !is.visited)
			var startPoint, startNextPoint = null, prevPoint, hole = false
			if (!currentPoint) {
				var unvisitedSegment = projectedSegments.find((is) => !is.visited)
				if (!unvisitedSegment) {
					break;
				} else {
					var p = unvisitedSegment[0]
					unvisitedSegment.visited = true
					console.log("unvisitedSegment", unvisitedSegment)
					hole = true
					currentPoint = unvisitedSegment[1]
					prevPoint = unvisitedSegment[0]
					startPoint = currentPoint
				}
			} else {
				var is = currentPoint
				onProjected = is.projectedStartpoint
					|| (!is.edgeEndpoint && outsideVector(currentPoint.projectedSegment, face.plane.normal)
						.dot(segmentVector(is.faceSegment)) > 0)
				if (onProjected) {
					prevPoint = currentPoint.faceSegment[0]
				} else {
					prevPoint = currentPoint.projectedSegment[0]
				}
				currentPoint = currentPoint.p
				startPoint = currentPoint
			}
			// calculate a - b
			//
			var isNewLoop = true
			var onProjectedNext
			var flipCount = 0
			do {
				var possibleDirections = possDirs(currentPoint, intersections)
				var currentDir = currentPoint.minus(prevPoint)
				var maxAngle = -Infinity, nextSegment
				possibleDirections.forEach((is) => {
					is.visited = loopCount
					if (!is.projectedEndpoint) {
						var dir = is.projectedSegment[1].minus(currentPoint)
						var angle = currentDir.angleRelativeNormal(dir, plane.normal)
						if (angle > maxAngle) {
							maxAngle = angle
							nextSegment = is.projectedSegment
							onProjectedNext = true
						}
					}
					if (!is.edgeEndpoint) {
						var dir = is.faceSegment[1].minus(currentPoint)
						var angle = currentDir.angleRelativeNormal(dir, plane.normal)
						if (angle > maxAngle) {
							maxAngle = angle
							nextSegment = is.faceSegment
							onProjectedNext = false
						}
					}
				})

				if (!onProjected) {
					face.vertices.forEach((v0, i, vertices) => {
						var v1 = vertices[(i + 1) % vertices.length]
						if (v0.like(currentPoint)) {
							var dir = v1.minus(v0)
							var angle = currentDir.angleRelativeNormal(dir, plane.normal)
							if (angle > maxAngle) {
								maxAngle = angle
								nextSegment = [v0, v1]
								onProjectedNext = false
							}
						}
					})
				}
				if (hole && projectedSegments.some((seg) => seg[0].like(prevPoint) && seg[1].like(currentPoint) && seg.visited && seg.visited != loopCount)) {
					isNewLoop = false
					break
				}
				projectedSegments.forEach((seg) => {
					if (seg[0].like(prevPoint) && seg[1].like(currentPoint)) {
						//console.log("ps", ps);
						seg.visited = loopCount} })
				if (onProjected) {
					projectedSegments.forEach((seg) => {
						var v0 = seg[0], v1 = seg[1]
						if (v0.like(currentPoint)) {
							var dir = v1.minus(v0)
							var angle = currentDir.angleRelativeNormal(dir, plane.normal)
							if (angle > maxAngle) {
								maxAngle = angle
								nextSegment = seg
								onProjectedNext = true
								console.log("here")
							}
						}
					})
				}
				assert(nextSegment)
				console.log("nextSegment", nextSegment.toSource())
				// simple edge point, find next edge
				flipCount += +(onProjected != onProjectedNext)
				onProjected = onProjectedNext
				console.log("flipCount", flipCount)
				newFaceVertices.push(currentPoint)
				prevPoint = currentPoint
				currentPoint = getNextPoint(nextSegment, currentPoint, intersections)
				if (!startNextPoint) {
					startNextPoint = currentPoint
				}
				console.log
					("asdhasjkdhKJ","currentpoint", currentPoint.toString(),
					"startNextPoint", startNextPoint.toString(),
					"prevPoint", prevPoint.toString(),
					"startPoint",startPoint.toString())
			} while (!(prevPoint == startPoint && startNextPoint == currentPoint) && i++ < 20 || i == 0)
			if (i >= 20) {
				assert(false, "too many")
			}
			if (!isNewLoop) {
				continue
			}
			newFaceVertices.pop()
			if (!hole && flipCount > 1) {
				newFaces.push(new BREP.Face(newFaceVertices, plane, face.name))
				console.log("add face 3")
				addedFaces++
			} else {
				console.log("is hole",(face.vertices[0]).ss, new BREP.Face(newFaceVertices).containsPoint(face.vertices[0]))
				console.log("newFaceVertices",  new BREP.Face(newFaceVertices).toString())
				console.log("reverseblash",new BREP.Face(newFaceVertices).containsPoint(face.vertices[0]))
				var loop =new BREP.Face(newFaceVertices, plane, face.name)
				var loopContainsFace = face.vertices.every(v => loop.containsPoint(v))
				var faceContainsLoop = loop.vertices.every(v => face.containsPoint(v))
				console.log("loopContainsFace", loopContainsFace, "faceContainsLoop", faceContainsLoop)
				if (loopContainsFace && faceContainsLoop) {
					// perfect overlap
				} else if (loopContainsFace) {
					newFaceVertices = face.vertices.slice()
					newFaces.push(new BREP.Face(newFaceVertices, plane, face.name))
					console.log("add face 4")
				} else if (faceContainsLoop) {
					newFaces.push(face.withHole(newFaceVertices))
					console.log("add face 2")
					addedFaces++
				} else {
					// outside each other
					newFaceVertices = face.vertices.slice()
					newFaces.push(new BREP.Face(newFaceVertices, plane, face.name))
					console.log("add face 1")
				}
				console.log("HOLE", newFaceVertices.toSource())
			}
		}
		if (!addedFaces) {
		}
	})
	return new BREP(null, null, newFaces, this.infiniteVolume)
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


function testBREP() {
	try {
		var f0 = new BREP.Face([V3(0, 0, 0), V3(10, 0, 0), V3(10, 10, 0), V3(0, 10, 0)], P3.XY)
		var f4 = new BREP.Face([V3(5, 0, 5), V3(5, 5, 5), V3(0, 5, 5), V3(0, 0, 5)], P3(V3.Z, 5))
		var boxBottom = new BREP.Face([V3(5, 0, 0), V3(5, 5, 0), V3(0, 5, 0), V3(0, 0, 0)].reverse(), P3(V3.Z.negated(), 0))
		var f00 = new BREP.Face([V3(0, 0, -6), V3(10, 0, -6), V3(11, 10, 0), V3(0, 10, 0)])
		var f01 = new BREP.Face([V3(10, 0, 0), V3(10, 0, -6), V3(11, 10, 0)])
		var f1 = new BREP.Face([V3(5, 8, -5), V3(9, 12, 1), V3(5, 8, 5)])
		var f2 = new BREP.Face([V3(5, 8, -5), V3(5, 8, 5), V3(1, 12, 2)])
		var f3 = new BREP.Face([V3(5, 8, -5), V3(1, 12, 2), V3(9, 12, 1)])
		var holeVertices = [V3(2, 3, 0), V3(8, 7, 0), V3(7, 2, 0)]
		//var holeFace = a.withHole(holeVertices)
		var a = new BREP(null, null, [f0])
		var b = BREP.tetrahedron(V3(5, 10, 0), V3(15, 10, 0), V3(5, 5, 0), V3(5, 10, -2))
		var a = new BREP(null, null, [new BREP.Face([V3(0, 0, 0), V3(0, 1, 0), V3(0, 1, 5), V3(0, 0, 5)])])
		var a = new BREP(null, null, [new BREP.Face([V3(0, 0, 0), V3(0, 1, 0), V3(0, 1, 5), V3(0, 0, 5)])])
		var a = new BREP(null, null, [ new BREP.Face([V3(5, 0, 0), V3(5, 5, 0), V3(0, 5, 0), V3(0, 0, 0)].reverse(), P3(V3.Z.negated(), 0))])
		//var a = new BREP(null, null, [new BREP.Face([V3(5, 0, 5), V3(0, 0, 5), V3(0, 5, 5), V3(5, 5, 5)], P3(V3(0, 0, 1), 5))])
		var a = BREP.box(5, 5, 5)
		var b = BREP.box(1, 1, 5).translate(V3(1, 1, 1))
		//var b = BREP.tetrahedron(V3(1, 0, -1), V3(1, 0, 8), V3(1, 4, -1), V3(4, 0, -1))
		//var b = BREP.tetrahedron(V3(2, -2, -2), V3(2, 4, -2), V3(1, -2, 4), V3(3, 4, 4))
		//var b = new BREP(null, null, [new BREP.Face([V3(1, 0, 0), V3(1, 1, 0), V3(0, 1, 0), V3(0, 0, 0)])])
		var concave = BREP.extrude(
			[V3(5, 5, 0), V3(5, 8, 0), V3(8, 10, -2), V3(2, 11, 2)].reverse(),
			P3.normalOnAnchor(V3(-2, 0, -3), V3(5, 5, 0)),
			V3(4, 0, 6))
		//var b = new BREP(null, null, [f1, f2, f3])
		//var b = tet
		//aMesh = a.toMesh()
		//bMesh = b.toMesh()
		console.log("css", b.ss())


		var c = a.minus(b)//clipped(b, true, true)
		cMesh = c.toNormalMesh()
		console.log("css", c.ss())
		//bMesh = b.clipped(a).toMesh()
	} catch (e) {
		console.log(e)
	}
	paintScreen()
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
	var dDistance = P3.throughPoints(a, b, c).distanceToPointSigned(d)
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
var singleColorShader, textureColorShader, singleColorShaderHighlight, arcShader, arcShader2,xyLinePlaneMesh,gl,cubeMesh,lightingShader

window.loadup = function () {
	/*
	var start = new Date().getTime();
	var m = M4.fromFunction(Math.random)
	for (var i = 0; i < 500000; ++i) {
		var  d= m.isMirroring()
	}

	console.log(m.determinant())
	var end = new Date().getTime();
	var time = end - start;
	console.log('Execution time: ' + time);
	*/

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

	assert(NLA.equals(-PI/2, V3.Y.times(3).angleRelativeNormal(V3.Z.times(-2), V3.X)))

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
	lightingShader = new GL.Shader($('vertexShaderLighting').text, $('fragmentShaderLighting').text);

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


	testBREP();
}
var aMesh, bMesh, cMesh
function paintScreen() {
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
	gl.loadIdentity();
	gl.scale(10, 10, 10);

	gl.loadIdentity();

	console.log("painting")


	if (aMesh) {
		singleColorShader.uniforms({
			color: rgbToArr4(COLORS.PP_STROKE)
		}).draw(aMesh, gl.LINES);
		singleColorShader.uniforms({
			color: rgbToArr4(COLORS.PP_FILL)
		}).draw(aMesh);
	}
	if (bMesh) {
		singleColorShader.uniforms({
			color: rgbToArr4(COLORS.TS_STROKE)
		}).draw(bMesh, gl.LINES);
		singleColorShader.uniforms({
			color: rgbToArr4(COLORS.TS_FILL)
		}).draw(bMesh);
	}
	/*
	if (cMesh) {
		singleColorShader.uniforms({
			color: rgbToArr4(COLORS.RD_STROKE)
		}).draw(cMesh, gl.LINES);
		singleColorShader.uniforms({
			color: rgbToArr4(COLORS.RD_FILL)
		}).draw(cMesh);
	}
	*/
	if (cMesh) {
		/*singleColorShader.uniforms({
			color: rgbToArr4(COLORS.RD_STROKE)
		}).draw(cMesh, gl.LINES);*/
		lightingShader.uniforms({
			color: rgbToArr4(COLORS.RD_FILL)
		}).draw(cMesh);
	}
	drawPlanes();
}