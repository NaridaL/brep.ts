/**
 * Created by aval on 21/12/2015.
 */
var P3 = NLA.Plane3

var goldenRatio = (1 + Math.sqrt(5)) / 2, u = V3(1, goldenRatio, 0).normalized(), s = u.x, t = u.y
var icosahedronVertices = [
	V3(-s,  t,  0),
	V3( s,  t,  0),
	V3(-s, -t,  0),
	V3( s, -t,  0),

	V3( 0, -s,  t),
	V3( 0,  s,  t),
	V3( 0, -s, -t),
	V3( 0,  s, -t),

	V3( t,  0, -s),
	V3( t,  0,  s),
	V3(-t,  0, -s),
	V3(-t,  0,  s),]
var icosahedronTriangles = [
	// 5 faces around point 0
	0, 11, 5,
	0, 5, 1,
	0, 1, 7,
	0, 7, 10,
	0, 10, 11,

	// 5 adjacent faces
	1, 5, 9,
	5, 11, 4,
	11, 10, 2,
	10, 7, 6,
	7, 1, 8,

	// 5 faces around point 3
	3, 9, 4,
	3, 4, 2,
	3, 2, 6,
	3, 6, 8,
	3, 8, 9,

	// 5 adjacent faces
	4, 9, 5,
	2, 4, 11,
	6, 2, 10,
	8, 6, 7,
	9, 8, 1,
]
function genFace(a, b, c, res, vertices, triangles) {
	var totalAngle = a.angleTo(b)
	var startIndex = vertices.length
	var ab = b.minus(a), ac = c.minus(a)
	var center = V3.add(a, b, c).div(3), center1 = center.normalized()

	// area of sphere triangle = sphere area / 20 = 4/ 20 * PI * R² (R == 1)
	// area of dedocahedron face sqrt(3) / 4 * a²
	var factor = 4 / 20 * Math.PI / (Math.sqrt(3) / 4 * ab.length() * ab.length())
	console.log("factor", factor)

	for (let i = 0; i<= res; i++) {
		var tAngle = (i / res - 0.5) * totalAngle
		var t = Math.tan(tAngle) / Math.tan(totalAngle / 2) / 2 + 0.5
		t = i /res
		for (let j = 0; j<= res - i; j++) {
			var sAngle = (j / res - 0.5) * totalAngle
			var s = Math.tan(sAngle) / Math.tan(totalAngle / 2) / 2 + 0.5
			s = j / res
			var s0 = s - 1/3, t0 = t - 1/3
			// cos satz
			var dist = Math.sqrt(s0 * s0 + t0 * t0 - 2 * s0 * t0 * Math.cos(Math.PI * 2 / 3)) * factor
			var distFactor = Math.sqrt(1 - dist * dist / 4)
			var weight = Math.min(s, t, (1 - t - s)) * Math.sin(Math.PI * 2 / 3) * 3
			distFactor = distFactor * weight + (1 - weight)
			console.log(i, j, "dist",dist, distFactor, s, t, Math.min(s * Math.sin(Math.PI * 2 / 3), t * Math.sin(Math.PI * 2 / 3)) * 3)
			var heightOffset = (-dist * dist / 2 + (1 - center.length())) * weight
			var v = V3.add(center, ab.times(s0 * distFactor), ac.times(t0 * distFactor), center1.times(heightOffset)).normalized()
			//var v = V3.add(center, ab.times(s0), ac.times(t0)).normalized()
			vertices.push([v.x, v.y, v.z])
		}
	}
	// add triangles
	for (let i = 0; i < res; i++) {
		for (let j = 0; j < res - i; j++) {
			if (0 != j) {
				triangles.push([startIndex + res - i + j + 1, startIndex + j, startIndex + res - i + j ])
			}
			triangles.push([startIndex + j + 1, startIndex + j, startIndex + res - i + j + 1])
		}
		startIndex += res - i + 1
	}
}
function st(res) {
	var mesh = new GL.Mesh({ normals: false, colors: false, lines:false});
	mesh.vertices.push(icosahedronVertices)
	for (var i = 0; i < 20; i++) {
		var [a, b, c] = icosahedronTriangles.slice(i * 3, i * 3 + 3).map(index => icosahedronVertices[index])
		//console.log(icosahedronTriangles.slice(i * 3, i * 3 + 3).map(index => icosahedronVertices[index]))
		genFace(a, b, c, res, mesh.vertices, mesh.triangles)
	}

	//console.log(mesh.toSource())
	mesh.computeWireframe()
	mesh.compile()
	return mesh
}

function tesselateRecursively(a, b, c, res, vertices, triangles, ia, ib, ic, lines) {
	if (0 == res) {
		triangles.push(ia, ib, ic)
		if (ia < ib) lines.push(ia, ib)
		if (ib < ic) lines.push(ib, ic)
		if (ic < ia) lines.push(ic, ia)
	} else {
		var abMid1 = a.plus(b).normalized(), bcMid1 = b.plus(c).normalized(), caMid1 = c.plus(a).normalized()
		var iabm = vertices.length, ibcm = iabm + 1, icam = iabm + 2
		vertices.push(abMid1, bcMid1, caMid1)
		tesselateRecursively(abMid1, bcMid1, caMid1, res - 1, vertices, triangles, iabm, ibcm, icam, lines)
		tesselateRecursively(a, abMid1, caMid1, res - 1, vertices, triangles, ia, iabm, icam, lines)
		tesselateRecursively(b, bcMid1, abMid1, res - 1, vertices, triangles, ib, ibcm, iabm, lines)
		tesselateRecursively(c, caMid1, bcMid1, res - 1, vertices, triangles, ic, icam, ibcm, lines)
	}
}
function sphereMesh(subdivisions) {
	var mesh = new GL.Mesh({normals: true, colors: false, lines: true});

	//mesh.vertices.push(a, b, c)
	//tesselateRecursively(a, b, c, 1, mesh.vertices, mesh.triangles, 0, 1, 2)

	mesh.vertices.pushAll(icosahedronVertices)
	for (var i = 0; i < 20; i++) {
		var [ia, ic, ib] = icosahedronTriangles.slice(i * 3, i * 3 + 3)
		tesselateRecursively(icosahedronVertices[ia], icosahedronVertices[ic], icosahedronVertices[ib], subdivisions, mesh.vertices, mesh.triangles, ia, ic, ib, mesh.lines)
	}

	mesh.normals = mesh.vertices
	mesh.compile()
	return mesh
}











