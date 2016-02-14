"use strict";
["min", "max", "PI", "sqrt","pow","round"].forEach(function (propertyName) {
	/*if (window[propertyName]) {
	 throw new Error("already exists"+propertyName)
	 }*/
	window[propertyName] = Math[propertyName];
});
/**
 * Created by aval on 21/12/2015.
 */
var M4 = NLA.Matrix4x4, V3 = NLA.Vector3, P3 = NLA.Plane3, L3 = NLA.Line3

var golden = (1 + Math.sqrt(5)) / 2, u = V3(1, golden, 0).normalized(), s = u.x, t = u.y
var vertices = [
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
var triangles = [
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

	for (var i = 0; i<= res; i++) {
		var tAngle = (i / res - 0.5) * totalAngle
		var t = Math.tan(tAngle) / Math.tan(totalAngle / 2) / 2 + 0.5
		t = i /res
		for (var j = 0; j<= res - i; j++) {
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
	for (var i = 0; i < res; i++) {
		for (var j = 0; j < res - i; j++) {
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
	var a = V3(1, 0, 0), b = V3(0, 1, 0), c = V3(0, 0, 1)
	//genFace(a, b, c, res, mesh.vertices, mesh.triangles)
	mesh.vertices.push(vertices)
	for (var i = 0; i < 20; i++) {
		var [a, b, c] = triangles.slice(i * 3, i * 3 + 3).map(index => vertices[index])
		console.log(triangles.slice(i * 3, i * 3 + 3).map(index => vertices[index]))
		genFace(a, b, c, res, mesh.vertices, mesh.triangles)
	}

	console.log(mesh.toSource())
	mesh.computeWireframe()
	mesh.compile()
	return mesh
}

function tesselateRecursively(a, b, c, res, vertices, triangles, ia, ib, ic, lines) {
	if (0 == res) {
		triangles.push([ia, ib, ic])
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
	var a = V3(1, 0, 0), b = V3(0, 1, 0), c = V3(0, 0, 1)
	var subdividecount = 3

	//mesh.vertices.push(a, b, c)
	//tesselateRecursively(a, b, c, 1, mesh.vertices, mesh.triangles, 0, 1, 2)

	mesh.vertices.pushAll(vertices)
	for (var i = 0; i < 20; i++) {
		var [ia, ic, ib] = triangles.slice(i * 3, i * 3 + 3)
		tesselateRecursively(vertices[ia], vertices[ic], vertices[ib], subdivisions, mesh.vertices, mesh.triangles, ia, ic, ib, mesh.lines)
	}

	mesh.normals = mesh.vertices
	mesh.compile()
	return mesh
}





































// abcd can be in any order
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
	gl = GL.create();
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
	//gl.enable(gl.BLEND);
	gl.enable(gl.DEPTH_TEST);
	gl.enable(gl.CULL_FACE);
	gl.depthFunc(gl.LEQUAL)
	//gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA); // TODO ?!

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
	lightingShader = new GL.Shader($('vertexShaderLighting').text, $('fragmentShaderLighting').text);
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

	aMesh = sphereMesh(3)

	//pointMesh = st()
	paintScreen();
}
var aMesh, bMesh, cMesh,pointMesh
function paintScreen() {
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
	gl.loadIdentity();
	gl.scale(10, 10, 10);

	gl.loadIdentity();
	/*
if (false) {
	gl.begin(gl.LINES)
	gl.color(0, 0, 0, 0)
	for (var v of pointMesh.vertices) {
		gl.vertex.apply(gl, v)
	}

	gl.end()
}else if (false)
{
	singleColorShader.uniforms({
		color: [0, 0, 0, 1]
	}).draw(pointMesh, gl.TRIANGLES);

}
*/

	if (aMesh) {
		lightingShader.uniforms({
			color: rgbToArr4(COLORS.TS_FILL)
		}).draw(aMesh);

		singleColorShader.uniforms({
			color: rgbToArr4(COLORS.TS_STROKE)
		})//.draw(aMesh, gl.LINES);
	}
	if (bMesh) {
		singleColorShader.uniforms({
			color: rgbToArr4(COLORS.TS_STROKE)
		}).draw(bMesh, gl.LINES);
		singleColorShader.uniforms({
			color: rgbToArr4(COLORS.TS_FILL)
		}).draw(bMesh);
	}
	if (cMesh) {
		singleColorShader.uniforms({
			color: rgbToArr4(COLORS.RD_STROKE)
		}).draw(cMesh, gl.LINES);
		singleColorShader.uniforms({
			color: rgbToArr4(COLORS.RD_FILL)
		}).draw(cMesh);
	}
	drawPlanes();
}