function parseGetParams() {
    const result = {}
    window.location.search
        .substr(1)
        .split("&")
        .forEach(function (item) {
            const tmp = item.split("=")
            result[tmp[0]] = decodeURI(tmp[1])
        })
    return result
}
let a, b, c, d, edges
function initB2() {
    dMesh = new GL.Mesh()
    /*
     var c1 = SemiEllipseCurve.circle(5), c2 = SemiEllipseCurve.circle(5, V3(3, 0))
     var test = new SemiEllipseCurve(V3(6, 1, 0), V3(3, 1, 0), V3(4, 0, 0))
     var cyl = new SemiCylinderSurface(new SemiEllipseCurve(V3.O, V3(5, 5, 0), V3(0, 5, 0)), V3.Z, 1)
     var ell = new SemiCylinderSurface(new SemiEllipseCurve(V3.O, V3(5, 5, 0), V3(0, 2, 0)), V3.Z, 1).rotateX(PI/3)

     aMesh = cyl.toMesh()
     bMesh = ell.toMesh()
     c1.isPointsWithEllipse(test)
     dMesh.compile()
     */
    eyePos = V(0, 100, 100)
    eyeFocus = V(0, 100, 0)
    eyeUp = V(0, 1, 0)
    zoomFactor = 1

    const gets = parseGetParams()
    "abcd".split('').forEach(k => gets[k] && (console.log(k + '=' + gets[k] + ';' + k + 'Mesh = ' + k + '.toMesh()'), eval(k + '=' + gets[k] + ';' + k + 'Mesh = ' + k + '.toMesh()')))

    //cMesh && cMesh.computeWireframeFromFlatTriangles() && cMesh.compile()
    if (gets['points']) {
        console.log("drPs from GET")
        drPs = eval(gets['points'])
    }

    if (gets['edges']) {
        console.log("edges from GET")
        dMesh = new GL.Mesh({triangles: false})
        edges = eval(gets['edges'])
        edges && dMesh.addVertexBuffer('curve1', 'curve1')
        edges.forEach(edge => {
            const points = edge.points()
            console.log(points)
            for (let i = 0; i < points.length - 1; i++) {
                dMesh.curve1.push(points[i], points[i + 1])
            }
        })
        console.log(dMesh.curve1)
    }
    if (gets['mesh']) {
        console.log("mesh from GET")
        sMesh = eval(gets['mesh'])
        sMesh.computeWireframeFromFlatTriangles()
        //sMesh.computeNormalLines()
        sMesh.compile()
    }
    if (gets['meshes']) {
        console.log("meshes from GET")
        b2meshes = eval(gets['meshes']) || []
        b2meshes.forEach(m => m.computeWireframeFromFlatTriangles())
        b2meshes.forEach(m => m.compile())
        console.log("meshes from GET", b2meshes)
        //sMesh.computeNormalLines()
    }
    if (gets['vectors']) {
        console.log("vectors from GET")
        drVs.pushAll(eval(gets['vectors']))
    }

    dMesh.compile()
}





const randomColors = NLA.arrayFromFunction(20, i => NLA.randomColor())
let aMesh: GL.Mesh, bMesh: GL.Mesh, cMesh: GL.Mesh, dMesh: GL.Mesh, sMesh: GL.Mesh, b2meshes: GL.Mesh[] = []
paintScreen = function() {
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
    gl.loadIdentity()

    drawVectors()

    gl.scale(100, 100, 100)

    if (aMesh) {
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        aMesh.lines && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.PP_STROKE) }).draw(aMesh, 'LINES')
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
        shaders.lighting.uniforms({ color: rgbToVec4(COLORS.PP_FILL),
            camPos: eyePos }).draw(aMesh)
    }

    if (bMesh) {
        gl.pushMatrix()
        //gl.translate(15, 0, 0)
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        bMesh.lines && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.PP_STROKE) }).draw(bMesh, 'LINES')
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
        shaders.lighting.uniforms({ color: rgbToVec4(COLORS.RD_FILL),
            camPos: eyePos }).draw(bMesh)
        bMesh.edgeTangents && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.TS_STROKE) })
            .drawBuffers({LGL_Vertex: bMesh.vertexBuffers.edgeTangents}, null, gl.LINES)
        bMesh.edgeTangents2 && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.RD_STROKE) })
            .drawBuffers({LGL_Vertex: bMesh.vertexBuffers.edgeTangents2}, null, gl.LINES)
        gl.popMatrix()
    }
    if (cMesh) {
        gl.pushMatrix()
        //gl.translate(30, 0, 0)
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        cMesh.lines && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.TS_STROKE) }).draw(cMesh, 'LINES')
        gl.projectionMatrix.m[11] += 1 / (1 << 20)

        let faceIndex = c.faces.length
        while (faceIndex--) {

            const face = c.faces[faceIndex]
            const faceTriangleIndexes = cMesh.faceIndexes.get(face)
            shaders.lighting.uniforms({
                color: rgbToVec4(randomColors[faceIndex % randomColors.length])
            }).draw(cMesh, 'TRIANGLES', faceTriangleIndexes.start, faceTriangleIndexes.count)
            //shaders.singleColor.uniforms({
            //color: rgbToVec4(0x0000ff)
            //}).draw(brepMesh, 'LINES')
        }

        gl.popMatrix()
    }
    if (dMesh && dMesh.hasBeenCompiled) {
        gl.pushMatrix()
        //gl.scale(10, 10, 10)
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        dMesh.lines && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.RD_STROKE) }).draw(dMesh, 'LINES')
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
        dMesh.triangles && shaders.lighting.uniforms({ color: rgbToVec4(0xffFF00),
            camPos: eyePos }).draw(dMesh)

        dMesh.curve1 && shaders.singleColor.uniforms({ color: rgbToVec4(0xff00000) })
            .drawBuffers({LGL_Vertex: dMesh.vertexBuffers.curve1}, null, gl.LINES)
        dMesh.curve2 && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.RD_STROKE) })
            .drawBuffers({LGL_Vertex: dMesh.vertexBuffers.curve2}, null, gl.LINES)

        dMesh.curve3 && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.PP_STROKE) })
            .drawBuffers({LGL_Vertex: dMesh.vertexBuffers.curve3}, null, gl.LINES)
        dMesh.curve4 && shaders.singleColor.uniforms({ color: rgbToVec4(0x00ff00) })
            .drawBuffers({LGL_Vertex: dMesh.vertexBuffers.curve4}, null, gl.LINES)
        gl.popMatrix()
    }
    if (sMesh && sMesh.hasBeenCompiled) {
        gl.pushMatrix()
        //gl.scale(10, 10, 10)
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        sMesh.lines && shaders.singleColor.uniforms({ color: rgbToVec4(0xFF6600) }).draw(sMesh, 'LINES')
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
        sMesh.triangles && shaders.lighting.uniforms({ color: rgbToVec4(0xffFF00),
            camPos: eyePos }).draw(sMesh)

        sMesh.curve1 && shaders.singleColor.uniforms({ color: rgbToVec4(0xff00000) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve1}, null, gl.LINES)
        sMesh.curve2 && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.RD_STROKE) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve2}, null, gl.LINES)

        sMesh.curve3 && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.PP_STROKE) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve3}, null, gl.LINES)
        sMesh.curve4 && shaders.singleColor.uniforms({ color: rgbToVec4(0x00ff00) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve4}, null, gl.LINES)
        gl.popMatrix()
    }
    for (const sMesh of b2meshes) {
        gl.pushMatrix()
        //gl.scale(10, 10, 10)
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        sMesh.lines && shaders.singleColor.uniforms({ color: rgbToVec4(0xFF6600) }).draw(sMesh, 'LINES')
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
        sMesh.triangles && shaders.lighting.uniforms({ color: rgbToVec4(0xffFF00),
            camPos: eyePos }).draw(sMesh)

        sMesh.curve1 && shaders.singleColor.uniforms({ color: rgbToVec4(0xff00000) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve1}, null, gl.LINES)
        sMesh.curve2 && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.RD_STROKE) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve2}, null, gl.LINES)

        sMesh.curve3 && shaders.singleColor.uniforms({ color: rgbToVec4(COLORS.PP_STROKE) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve3}, null, gl.LINES)
        sMesh.curve4 && shaders.singleColor.uniforms({ color: rgbToVec4(0x00ff00) })
            .drawBuffers({LGL_Vertex: sMesh.vertexBuffers.curve4}, null, gl.LINES)
        gl.popMatrix()
    }

    //drPs.forEach(v => drawPoint(v, undefined, 0.3))
    drawPoints(0.05)
    planes.forEach(plane => drawPlane(plane))
}



















//var sketchPlane = new CustomPlane(V3.X, V3(1, 0, -1).unit(), V3.Y, -500, 500, -500, 500, 0xff00ff);
var planes = [
    new CustomPlane(V3.O, V3.Y, V3.Z, -500, 500, -500, 500, 0xff0000),
    new CustomPlane(V3.O, V3.X, V3.Z, -500, 500, -500, 500, 0x00ff00),
    new CustomPlane(V3.O, V3.X, V3.Y, -500, 500, -500, 500, 0x0000ff),
    //	sketchPlane
]

//const font = opentype.loadSync('fonts/Verdana.ttf')
//opentype.load('fonts/calibri.ttf', function(err, font) {
//opentype.load('fonts/FiraSans-Medium.woff', function(err, font) {
window.onload = function () {
    B2T.loadFonts().then(function () {
        modeStack.push({})
        window.onerror = function (errorMsg, url, lineNumber, column, errorObj) {
            console.log(errorMsg, url, lineNumber, column, errorObj)
        }
        gl = GL.create({canvas: document.getElementById("testcanvas") as HTMLCanvasElement})
        gl.fullscreen()
        gl.canvas.oncontextmenu = () => false

        setupCamera()
        //gl.cullFace(gl.FRONT_AND_BACK);
        gl.clearColor(1.0, 1.0, 1.0, 0.0)
        gl.enable(gl.BLEND)
        gl.enable(gl.DEPTH_TEST)
        gl.enable(gl.CULL_FACE)
        gl.depthFunc(gl.LEQUAL)
        gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA) // TODO ?!

        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
        gl.loadIdentity()
        gl.scale(10, 10, 10)

        gl.loadIdentity()


        initMeshes()
        initShaders()
        initNavigationEvents()
        initPointInfoEvents()
        initB2()
        setupCamera()
        paintScreen()
    })
}
