// TODO: V3.create instead of V3 where necessar
function parseGetParams() {
	var result = {}
	window.location.search
		.substr(1)
		.split("&")
		.forEach(function (item) {
			var tmp = item.split("=");
			result[tmp[0]] = decodeURI(tmp[1])
		});
	return result;
}
function initB2() {
	dMesh = new GL.Mesh()
	/*
	var c1 = EllipseCurve.circle(5), c2 = EllipseCurve.circle(5, V3(3, 0))
	var test = new EllipseCurve(V3(6, 1, 0), V3(3, 1, 0), V3(4, 0, 0))
	var cyl = new CylinderSurface(new EllipseCurve(V3.ZERO, V3(5, 5, 0), V3(0, 5, 0)), V3.Z, 1)
	var ell = new CylinderSurface(new EllipseCurve(V3.ZERO, V3(5, 5, 0), V3(0, 2, 0)), V3.Z, 1).rotateX(PI/3)

	aMesh = cyl.toMesh()
	bMesh = ell.toMesh()
	c1.isPointsWithEllipse(test)
	dMesh.compile()
	*/
	eyePos = V3(0, 100, 100)
	eyeFocus = V3(0, 100, 0)
	eyeUp = V3(0, 1, 0)
	zoomFactor = 0.5

	var face = PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(10, 0), V3(10, 10), V3(0, 10)])
	// splitting contour in base position:
	var brep = B2.extrudeVertices([V3(5, 0), V3(2, 3), V3(8, 3)], P3.XY.flipped(), V3(0, 0, 10)).translate(0, 0, -2).flipped()

	aMesh = face.inB2().toMesh()
	bMesh = brep.toMesh()

	var gets = parseGetParams()
	"abc".split('').forEach(c => gets[c] && (console.log(c+" from GET: ", gets[c]), eval(c+'Mesh = '+gets[c] + '.toMesh()')))

	if (gets['points']) {
		console.log("drPs from GET")
		drPs = eval(gets['points'])
	}

	if (gets['edges']) {
		console.log("edges from GET")
		dMesh = new GL.Mesh({triangles: false})
		dMesh.addVertexBuffer('curve1', 'curve1')
		var edges = eval(gets['edges'])
		edges.forEach(edge => {
			var points = edge.points
			for (var i = 0; i < points.length - 1; i++) {
				dMesh.curve1.push(points[i], points[i + 1])
			}
		})
		console.log(dMesh.curve1)
	}

	dMesh.compile()
}






var aMesh, bMesh, cMesh, dMesh
function paintScreen2() {
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.loadIdentity();
    gl.scale(10, 10, 10);

    gl.loadIdentity();

    drawVectors()

    gl.scale(10, 10, 10)

    if (aMesh) {
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        aMesh.lines && singleColorShader.uniforms({ color: rgbToVec4(COLORS.PP_STROKE) }).draw(aMesh, 'LINES');
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
        lightingShader.uniforms({ color: rgbToVec4(COLORS.PP_FILL),
            camPos: eyePos }).draw(aMesh);
    }
    if (bMesh) {
        gl.pushMatrix()
        //gl.translate(15, 0, 0)
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        bMesh.lines && singleColorShader.uniforms({ color: rgbToVec4(COLORS.PP_STROKE) }).draw(bMesh, 'LINES');
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
        lightingShader.uniforms({ color: rgbToVec4(COLORS.RD_FILL),
            camPos: eyePos }).draw(bMesh);
        bMesh.edgeTangents && singleColorShader.uniforms({ color: rgbToVec4(COLORS.TS_STROKE) })
            .drawBuffers({gl_Vertex: bMesh.vertexBuffers.edgeTangents}, null, gl.LINES)
        bMesh.edgeTangents2 && singleColorShader.uniforms({ color: rgbToVec4(COLORS.RD_STROKE) })
            .drawBuffers({gl_Vertex: bMesh.vertexBuffers.edgeTangents2}, null, gl.LINES)
        gl.popMatrix()
    }
    if (cMesh) {
        gl.pushMatrix()
        //gl.translate(30, 0, 0)
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        cMesh.lines && singleColorShader.uniforms({ color: rgbToVec4(COLORS.TS_STROKE) }).draw(cMesh, 'LINES');
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
        lightingShader.uniforms({ color: rgbToVec4(COLORS.RD_FILL),
            camPos: eyePos }).draw(cMesh)

        gl.popMatrix()
    }
    if (dMesh) {
        gl.pushMatrix()
        gl.translate(20, 0, 0)
        //gl.scale(10, 10, 10)
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        dMesh.lines && singleColorShader.uniforms({ color: rgbToVec4(COLORS.RD_STROKE) }).draw(dMesh, 'LINES');
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
        lightingShader.uniforms({ color: rgbToVec4(0xffFF00),
            camPos: eyePos }).draw(dMesh)

        dMesh.curve1 && singleColorShader.uniforms({ color: rgbToVec4(0xff00000) })
            .drawBuffers({gl_Vertex: dMesh.vertexBuffers.curve1}, null, gl.LINES)
        dMesh.curve2 && singleColorShader.uniforms({ color: rgbToVec4(COLORS.RD_STROKE) })
            .drawBuffers({gl_Vertex: dMesh.vertexBuffers.curve2}, null, gl.LINES)

        dMesh.curve3 && singleColorShader.uniforms({ color: rgbToVec4(COLORS.PP_STROKE) })
            .drawBuffers({gl_Vertex: dMesh.vertexBuffers.curve3}, null, gl.LINES)
        dMesh.curve4 && singleColorShader.uniforms({ color: rgbToVec4(0x00ff00) })
            .drawBuffers({gl_Vertex: dMesh.vertexBuffers.curve4}, null, gl.LINES)
        gl.popMatrix()
    }

    drPs.forEach(v => {
        gl.pushMatrix()
        gl.translate(v)
        gl.scale(0.3,0.3,0.3)
        lightingShader.uniforms({color: rgbToVec4(0xabcdef)}).draw(sMesh)
        //singleColorShader.uniforms({ color: rgbToVec4(COLORS.RD_STROKE) }).draw(sMesh, 'LINES')
        gl.popMatrix()
    })
    drawPlanes();
}



















//var sketchPlane = new CustomPlane(V3.X, V3(1, 0, -1).unit(), V3.Y, -500, 500, -500, 500, 0xff00ff);
var planes = [
    CustomPlane(V3.ZERO, V3.Y, V3.Z, -500, 500, -500, 500, 0xff0000),
    CustomPlane(V3.ZERO, V3.X, V3.Z, -500, 500, -500, 500, 0x00ff00),
    CustomPlane(V3.ZERO, V3.X, V3.Y, -500, 500, -500, 500, 0x0000ff),
    //	sketchPlane
];

var singleColorShader, textureColorShader, singleColorShaderHighlight, arcShader, arcShader2,xyLinePlaneMesh,gl,cubeMesh,lightingShader, vectorMesh

var sMesh

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

    setupCamera();
    //gl.cullFace(gl.FRONT_AND_BACK);
    gl.clearColor(1.0, 1.0, 1.0, 0.0);
    gl.enable(gl.BLEND);
    gl.enable(gl.DEPTH_TEST);
    gl.enable(gl.CULL_FACE);
    gl.depthFunc(gl.LEQUAL)
    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA); // TODO ?!

    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.loadIdentity();
    gl.scale(10, 10, 10);

    gl.loadIdentity();

    gl.onmousemove = function (e) {
        if (e.dragging) {
            if (0x4 & e.buttons) {
                // pan
                var moveCamera = V3(-e.deltaX * 2 / gl.canvas.width, e.deltaY * 2 / gl.canvas.height, 0);
                var inverseProjectionMatrix = gl.projectionMatrix.inversed();
                var worldMoveCamera = inverseProjectionMatrix.transformVector(moveCamera);
                eyePos = eyePos.plus(worldMoveCamera);
                eyeFocus = eyeFocus.plus(worldMoveCamera);
                setupCamera();
                paintScreen2();
            }
            if (0x2 & e.buttons) {
                var rotateLR = deg2rad(-e.deltaX / 6.0);
                var rotateUD = deg2rad(-e.deltaY / 6.0);

                // rotate
                var matrix = M4.rotationLine(eyeFocus, eyeUp, rotateLR)
                //var horizontalRotationAxis = eyeFocus.minus(eyePos).cross(eyeUp)
                var horizontalRotationAxis = eyeUp.cross(eyePos.minus(eyeFocus))
                matrix = matrix.times(M4.rotationLine(eyeFocus, horizontalRotationAxis, rotateUD))
                eyePos = matrix.transformPoint(eyePos)
                eyeUp = matrix.transformVector(eyeUp)

                setupCamera();
                paintScreen2();
            }
        }
    }
    xyLinePlaneMesh = new GL.Mesh({lines: true, triangles: false});
    xyLinePlaneMesh.vertices = [[0, 0], [0, 1], [1, 1], [1, 0]];
    xyLinePlaneMesh.lines = [[0, 1], [1, 2], [2, 3], [3, 0]];
    xyLinePlaneMesh.compile();
    vectorMesh = GL.Mesh.rotation([V3.ZERO, V3(0, 0.05, 0), V3(0.8, 0.05), V3(0.8, 0.1), V3(1, 0)], L3.X, Math.PI * 2, 8, false, normals)
    sMesh = GL.Mesh.sphere(2)

    singleColorShader = new GL.Shader(vertexShaderBasic, fragmentShaderColor);
    singleColorShaderHighlight = new GL.Shader(vertexShaderBasic, fragmentShaderColorHighlight);
    //textureColorShader = new GL.Shader(vertexShaderTextureColor, fragmentShaderTextureColor);
    arcShader = new GL.Shader(vertexShaderRing, fragmentShaderColor);
    arcShader2 = new GL.Shader(vertexShaderArc, fragmentShaderColor);
    lightingShader = new GL.Shader(vertexShaderLighting, fragmentShaderLighting);

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
        paintScreen2();
    });
    initB2()
    setupCamera()
    paintScreen2()

}