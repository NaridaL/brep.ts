///<reference path="surface/ConicSurface.ts"/>

///<reference path="CustomPlane.ts"/>

function parseGetParams(str: string) {
    const result: {[key:string]: string} = {}
    str
            .split('&')
            .forEach(function (item) {
                const splitIndex = item.indexOf('=')
                result[item.substr(0, splitIndex)] = decodeURI(item.substr(splitIndex + 1))
            })
    return result
}
const COLORS = {
    RD_FILL: chroma('#9EDBF9'),
    RD_STROKE: chroma('#77B0E0'),
    TS_FILL: chroma('#D19FE3'),
    TS_STROKE: chroma('#A76BC2'),
    PP_FILL: chroma('#F3B6CF'),
    PP_STROKE: chroma('#EB81B4'),
}
class BREPGLContext extends (Object as any as typeof LightGLContext) {
    shaders: typeof SHADER_TYPE_VAR

    cachedMeshes: WeakMap<any, Mesh & { TRIANGLES: int[], normals: V3[] }> = new WeakMap()

    static create(gl: LightGLContext) {
        addOwnProperties(gl, BREPGLContext.prototype)
        addOwnProperties(gl, new BREPGLContext(gl as BREPGLContext))
        return gl as BREPGLContext
    }

    constructor(gl: BREPGLContext) {
        super(gl)
        this.shaders = initShaders(gl)
        initMeshes(this.meshes, this)
    }

    drawPoint(p: V3, color: GL_COLOR = GL_COLOR_BLACK, size = 5) {
        this.pushMatrix()
        this.translate(p)
        this.scale(size, size, size)
        this.shaders.singleColor.uniforms({color: color}).draw(this.meshes.sphere1)
        this.popMatrix()
    }

    drawEdge(edge: Edge, color: GL_COLOR = GL_COLOR_BLACK, width = 2) {
        CURVE_PAINTERS[edge.curve.constructor.name](this, edge.curve, color, edge.minT, edge.maxT, width)
    }

    drawCurve(curve: Curve, color: GL_COLOR = GL_COLOR_BLACK, width = 2, tStart: number, tEnd: number) {
        CURVE_PAINTERS[curve.constructor.name](this, curve, color, tStart, tEnd, width)
    }

    drawVector(vector: V3, anchor: V3, color: GL_COLOR = GL_COLOR_BLACK, size = 1) {
        this.pushMatrix()

        const vT = vector.getPerpendicular().unit()
        this.multMatrix(M4.forSys(vector, vT, vector.cross(vT).unit(), anchor))
        1 != size && this.scale(size, size, size)
        this.shaders.singleColor.uniforms({
            color: color
        }).draw(this.meshes.vector)

        this.popMatrix()
    }

    drawVectors(drVs: {dir1: V3, anchor: V3, color: GL_COLOR}[]) {
        this.drawVector(V3.X, V3.O, chroma('red').gl(), undefined)
        this.drawVector(V3.Y, V3.O, chroma('green').gl(), undefined)
        this.drawVector(V3.Z, V3.O, chroma('blue').gl(), undefined)

        drVs.forEach(vi => this.drawVector(vi.dir1, vi.anchor, vi.color, undefined))
    }

    drawPlane(customPlane: CustomPlane, color: GL_COLOR, dotted: boolean = false) {
        this.pushMatrix()
        this.multMatrix(M4.forSys(customPlane.right, customPlane.up, customPlane.normal1))
        this.translate(customPlane.sMin, customPlane.rMin, customPlane.w)
        this.scale(customPlane.sMax - customPlane.sMin, customPlane.tMax - customPlane.rMin, 1)

        const mesh = dotted ? this.meshes.xyDottedLinePlane : this.meshes.xyLinePlane
        this.shaders.singleColor.uniforms({color: color}).draw(mesh, DRAW_MODES.LINES)

        this.popMatrix()
    }
}
function conicPainter(mode: 0 | 1 | 2, gl: BREPGLContext, ellipse: SemiEllipseCurve, color: GL_COLOR, startT: number, endT: number, width = 2) {
    gl.shaders.ellipse3d.uniforms({
        f1: ellipse.f1,
        f2: ellipse.f2,
        center: ellipse.center,
        color: color,
        startT: startT,
        endT: endT,
        scale: width,
        mode: mode
    }).draw(gl.meshes.pipe)
}
const CURVE_PAINTERS: {[curveConstructorName: string]: (gl: BREPGLContext, curve: Curve, color: GL_COLOR, startT: number, endT: number, width: number) => void} = {
    [SemiEllipseCurve.name]: conicPainter.bind(undefined, 0),
    [EllipseCurve.name]: conicPainter.bind(undefined, 0),
    [ParabolaCurve.name]: conicPainter.bind(undefined, 1),
    [HyperbolaCurve.name]: conicPainter.bind(undefined, 2),
    [ImplicitCurve.name](gl, curve: ImplicitCurve, color, startT, endT, width = 2, normal = V3.Z) {
        let mesh = gl.cachedMeshes.get(curve)
        if (!mesh) {
            mesh = new Mesh()
                .addIndexBuffer('TRIANGLES')
                .addVertexBuffer('normals', 'LGL_Normal')
            curve.addToMesh(mesh)
            mesh.compile()
            //mesh=Mesh.sphere(2)
            gl.cachedMeshes.set(curve, mesh)
        }
        // TODO: draw only part
        //startT: startT,
        //	endT: endT,
        gl.shaders.generic3d.uniforms({
            color: color,
            scale: width,
        }).draw(mesh)
    },
    [BezierCurve.name](gl, curve: BezierCurve, color, startT, endT, width = 2, normal = V3.Z) {
        gl.shaders.bezier3d.uniforms({
            p0: curve.p0,
            p1: curve.p1,
            p2: curve.p2,
            p3: curve.p3,
            color: color,
            startT: startT,
            endT: endT,
            scale: width,
            normal: normal
        }).draw(gl.meshes.pipe)
    },
    [L3.name](gl, curve: L3, color, startT, endT, width = 2, normal = V3.Z) {
        gl.pushMatrix()
        const a = curve.at(startT), b = curve.at(endT)
        const ab = b.minus(a), abT = ab.getPerpendicular().unit()
        const m = M4.forSys(ab, abT, ab.cross(abT).unit(), a)
        gl.multMatrix(m)
        gl.scale(1, width, width)
        gl.shaders.singleColor.uniforms({
            color: color, // TODO: error checking
        }).draw(gl.meshes.pipe)

        gl.popMatrix()
    },
}
CURVE_PAINTERS[PICurve.name] = CURVE_PAINTERS[ImplicitCurve.name]


let setupCameraListener: (e: typeof eye) => void
const SHADER_TYPE_VAR = (false as true) && initShaders(0 as any)
// let shaders: typeof SHADER_TYPE_VAR
let a: B2, b: B2, c: B2, d: B2, edges: Edge[] = [], hovering: any,
        wireframe: boolean = false, normallines: boolean = false, b2s: B2[] = []
let eye = {pos: V(1000, 1000, 1000), focus: V3.O, up: V3.Z, zoomFactor: 1}
let hoverHighlight: any
let drPs: (V3 | {info: string, p: V3})[] = [], drVs: any[] = []
const edgeViewerColors = arrayFromFunction(20, i => chroma.random().gl())
function initB2() {
    eye.pos = V(1, 2, 101)
    eye.focus = V(0, 1, 0)
    eye.up = V(0, 1, 0)
    eye.zoomFactor = 8

    const hash = window.location.search.substr(1) || window.location.hash.substr(1)
    //noinspection TsLint
    let i, points, vectors, mesh, hjk
    //Object.assign(window, HJK())
    eval(decodeURI(hash))
    let gets: any = {a, b, c, d, mesh, edges, points, vectors}
    hjk && (gets = HJK())
    i && Object.assign(eye, i)
    'abcd'.split('').forEach(k => gets[k] && eval(`aMeshes.push((${k} = gets.${k}).toMesh()); b2s.push(${k})`))

    //cMesh && cMesh.computeWireframeFromFlatTriangles() && cMesh.compile()
    if (gets.points) {
        console.log('drPs from GET')
        drPs = gets.points
    }
    if (gets.vectors) {
        console.log('vectors from GET')
        drVs.pushAll(gets.vectors)
    }

    for (let i = 0; i < aMeshes.length; i++) {
        aMeshes[i].computeWireframeFromFlatTriangles('wireframe')
        aMeshes[i].computeNormalLines(0.1, 'normallines')
        aMeshes[i].compile()
    }

    if (gets.edges) {
        console.log('edges from GET')
        edges = gets.edges
        dMesh = new Mesh()
            .addIndexBuffer('TRIANGLES')
            .addVertexBuffer('normals', 'LGL_Normal')
            .addVertexBuffer('curve1', 'curve1')
            .addVertexBuffer('curve1colors', 'curve1colors')
        edges.forEach((edge, edgeIndex) => {
            const points = edge.points()
            for (let i = 0; i < points.length - 1; i++) {
                const color = edgeViewerColors[(edgeIndex + (i % 2)) % edgeViewerColors.length]
                // const tangent = edge.tangentAt(i)
                // dMesh.curve1.push(points[i], points[i].plus(tangent.toLength(1)))
                dMesh.curve1.push(points[i], points[i + 1])
                dMesh.curve1colors.push(color, color)
            }
            edge.curve instanceof PICurve && (edge.curve as PICurve).addToMesh(dMesh, 8, 0.02, 2)
        })
        //dMesh.computeWireframeFromFlatTriangles()
    }
    if (gets.mesh) {
        console.log('mesh/es from GET', b2meshes)
        b2meshes = gets.mesh instanceof Array ? gets.mesh : [gets.mesh]
        b2meshes.forEach(m => m.computeWireframeFromFlatTriangles())
        b2meshes.forEach(m => m.computeNormalLines(0.5))
        b2meshes.forEach(m => m.compile())
    }

    dMesh.compile()
}



const meshColors: any[][] = [
    chroma.scale(['#ff297f', '#6636FF']).mode('lab').colors(20, null),
    chroma.scale(['#ffe93a', '#ff6e35']).mode('lab').colors(20, null),
    chroma.scale(['#1eff33', '#4960ff']).mode('lab').colors(20, null),
    chroma.scale(['#31fff8', '#2dff2a']).mode('lab').colors(20, null)
]
const meshColorssGL = meshColors.map(cs => cs.map(c => c.gl()))
let aMeshes: (Mesh & {faceIndexes?: Map<Face, {start: int, count: int}>, TRIANGLES: int[], normals: V3[]})[] = [],
        //bMesh: Mesh & {faceIndexes?: Map<Face, {start: int, count: int}>},
        //cMesh: Mesh & {faceIndexes?: Map<Face, {start: int, count: int}>},
    dMesh: Mesh & { faceIndexes?: Map<Face, { start: int, count: int }>, TRIANGLES: int[], normals: V3[], curve1: V3[], curve1colors: GL_COLOR[] },
    b2meshes: (Mesh & { TRIANGLES: int[], normals: V3[] })[] = []
function viewerPaint(time: int, gl: BREPGLContext) {
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
    gl.loadIdentity()

    gl.drawVectors(drVs)
    gl.shaders.lighting.uniforms({ camPos: eye.pos })
    //gl.scale(100, 100, 100)
    for (let i = 0; i < aMeshes.length; i++) {
        const aMesh = aMeshes[i]
        gl.pushMatrix()
        //gl.translate(30, 0, 0)
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        wireframe && gl.shaders.singleColor.uniforms({ color: COLORS.TS_STROKE.gl() })
                .drawBuffers(aMesh.vertexBuffers, aMesh.indexBuffers.wireframe, DRAW_MODES.LINES)
        normallines && gl.shaders.singleColor.uniforms({ color: COLORS.TS_STROKE.gl() })
                .drawBuffers(aMesh.vertexBuffers, aMesh.indexBuffers.normallines, DRAW_MODES.LINES)
        gl.shaders.singleColor.uniforms({ color: COLORS.TS_STROKE.gl() })
                .drawBuffers(aMesh.vertexBuffers, aMesh.indexBuffers.LINES, DRAW_MODES.LINES)
        gl.projectionMatrix.m[11] += 1 / (1 << 20)

        let faceIndex = b2s[i].faces.length
        while (faceIndex--) {
            const face = b2s[i].faces[faceIndex]
            const faceTriangleIndexes = aMesh.faceIndexes.get(face)
            gl.shaders.lighting.uniforms({
                color: hovering == face ? meshColors.emod(i).emod(faceIndex).darken(2).gl() : meshColorssGL.emod(i).emod(faceIndex)
            }).draw(aMesh, DRAW_MODES.TRIANGLES, faceTriangleIndexes.start, faceTriangleIndexes.count)
        }

        gl.popMatrix()
    }

    for (const sMesh of b2meshes) {
        gl.pushMatrix()
        //gl.scale(10, 10, 10)
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        sMesh.LINES && gl.shaders.singleColor.uniforms({ color: hexIntToGLColor(0xFF6600) }).draw(sMesh, DRAW_MODES.LINES)
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
        sMesh.TRIANGLES && gl.shaders.lighting.uniforms({ color: hexIntToGLColor(0xffFF00),
            camPos: eye.pos }).draw(sMesh)
        gl.popMatrix()
    }

    if (hovering instanceof Edge) {
        gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
        gl.drawEdge(hovering, GL_COLOR_BLACK, 2 / eye.zoomFactor)
        gl.projectionMatrix.m[11] += 1 / (1 << 20)
    }
    //edges.forEach((e, i) => drawEdge(e, 0x0ff000, 0.01))

    //drPs.forEach(v => drawPoint(v, undefined, 0.3))
    drPs.forEach(info => gl.drawPoint(info instanceof V3 ? info : info.p, hexIntToGLColor(0xcc0000), 5 / eye.zoomFactor))
    b2planes.forEach(plane => gl.drawPlane(plane, hexIntToGLColor(plane.color), hoverHighlight))

    //console.log(gl.drawCallCount)
}
// let meshes: any = {}











/**
 * Transforms position on the screen into a line in world coordinates.
 */
function getMouseLine(pos: { x: number; y: number }, _gl: LightGLContext): L3 {
    const ndc1 = V(pos.x * 2 / _gl.canvas.width - 1, -pos.y * 2 / _gl.canvas.height + 1, 0)
    const ndc2 = V(pos.x * 2 / _gl.canvas.width - 1, -pos.y * 2 / _gl.canvas.height + 1, 1)
    //console.log(ndc)
    const inverseProjectionMatrix = _gl.projectionMatrix.inversed()
    const s = inverseProjectionMatrix.transformPoint(ndc1)
    const dir = inverseProjectionMatrix.transformPoint(ndc2).minus(s)
    return L3.anchorDirection(s, dir)
}


function getHovering(mouseLine: L3, faces: Face[], planes: CustomPlane[], points: V3[], edges: Edge[],
                     mindist: number,
                     ...consider: ('faces' | 'planes' | 'sketchElements' | 'points' | 'edges' | 'features')[]): any {
    let hoverHighlight = null, nearest = Infinity

    const checkFeatures = consider.includes('features')
    assert(!checkFeatures || !consider.includes('faces'))

    function checkEl(el: any, distance: number) {
        if (distance < nearest) {
            nearest = distance
            hoverHighlight = el
        }
    }

    if (faces && (consider.includes('faces') || consider.includes('features'))) {
        for (const face of faces) {
            checkEl(checkFeatures ? face.info.feature : face, face.intersectsLine(mouseLine))
        }
    }
    if (planes && consider.includes('planes')) {
        for (const plane of planes) {
            checkEl(plane, plane.distanceTo(mouseLine,mindist))
        }
    }
    if (consider.includes('points')) {
        for (const p of points) {
            const t = mouseLine.pointT(p)
            if (mouseLine.at(t).distanceTo(p) < mindist * 1.2) {
                checkEl(p, t - 0.1)
            }
        }
    }
    if (consider.includes('edges')) {
        const projPlane = new P3(mouseLine.dir1, 0)
        const projPoint = projPlane.projectedPoint(mouseLine.anchor)
        for (const edge of edges) {
            const curve = edge.curve
            const prio = 0.05
            if (curve instanceof L3 && curve.dir1.isParallelTo(mouseLine.dir1)) {
                const d = mouseLine.distanceToPoint(edge.a)
                const t = mouseLine.pointT(edge.a)

                if (d < mindist) {
                    checkEl(edge, t - prio)
                }
            } else {
                if (!(curve instanceof ImplicitCurve)) {
                    const projCurve = curve.project(projPlane)
                    const curveT = edge.clampedT(projCurve.closestTToPoint(projPoint))
                    const p = curve.at(curveT)
                    const t = mouseLine.pointT(p)
                    if (projCurve.at(curveT).distanceTo(projPoint) < mindist) {
                        checkEl(edge, t - prio)
                    }
                }
            }
        }
    }

    return hoverHighlight
}




function initInfoEvents(paintScreen: () => {}, gl: BREPGLContext) {
    gl.canvas.addEventListener('mousemove', function (e) {
        const mouseLine = getMouseLine({x: e.clientX, y: e.clientY}, gl)
        const faces = [a,b,c,d].mapFilter(b2 => b2 && b2.faces).concatenated()
        const testEdges: Edge[] = [a,b,c,d].mapFilter(b2 => b2 && (b2.buildAdjacencies(), Array.from(b2.edgeFaces.keys())))
                .concatenated()
                .concat(edges)
        hovering = getHovering(mouseLine, faces, undefined, [], testEdges, 0.1, 'faces', 'edges')
        //let html = '', pp
        //if (hovering instanceof Edge) {
        //	pp = V(e.clientX, e.clientY)
        //	defaultRoundFunction = x => round10(x, -3)
        //	html = hovering.toString(x => round10(x, -3)) + ' length=' + hovering.length().toFixed(3)
        //} else if (hovering instanceof Face) {
        //	pp = V(e.clientX, e.clientY)
        //	defaultRoundFunction = x => round10(x, -3)
        //   let area
        //   try { area = hovering.calcArea() } catch (e) {}
        //	html = `face surface=${hovering.surface.constructor.name} edges=${hovering.contour.length} area=${area}`
        //}
        //if (pp) {
        //	//const pSC = gl.projectionMatrix.times(gl.modelViewMatrix).transformPoint(pp)
        //	//const x = (pSC.x * 0.5 + 0.5) * window.innerWidth, y = (-pSC.y * 0.5 + 0.5) * window.innerHeight
        //	tooltipShow(html, pp.x, pp.y)
        //} else {
        //	tooltipHide()
        //}
        paintScreen()
    })
}


//var sketchPlane = new CustomPlane(V3.X, V3(1, 0, -1).unit(), V3.Y, -500, 500, -500, 500, 0xff00ff);
const b2planes = [
    new CustomPlane(V3.O, V3.Y, V3.Z, 'planeYZ', 0xff0000),
    new CustomPlane(V3.O, V3.X, V3.Z, 'planeZX', 0x00ff00),
    new CustomPlane(V3.O, V3.X, V3.Y, 'planeXY', 0x0000ff),
    //	sketchPlane
]

async function viewerMain() {
    const paintScreen = () => requestAnimationFrame(t => viewerPaint(t, gl))
    await B2T.loadFonts()
    window.onerror = function (errorMsg, url, lineNumber, column, errorObj) {
        console.log(errorMsg, url, lineNumber, column, errorObj)
    }
    const gl = BREPGLContext.create(LightGLContext.create({canvas: document.getElementById('testcanvas') as HTMLCanvasElement}))
    gl.fullscreen()
    gl.canvas.oncontextmenu = () => false

    setupCamera(eye, gl)
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

    initNavigationEvents(gl, eye, paintScreen)
    setupCameraListener = function (eye) {
        const iSource = 'i=' + eye.toSource().replace(/[\n\r\s]+|^\(|\)$/g, '')
        const hash = window.location.hash.substr(1) || iSource
        const result = hash.match(/i=\{[^}]*\}/)
                ? hash.replace(/i=\{[^}]*\}/, iSource)
                : hash + ';' + iSource
        window.history.replaceState(undefined, undefined, '#' + result)
    }
    initInfoEvents(paintScreen, gl)
    //initToolTips() // hide tooltip on mouseover
    //initPointInfoEvents()
    initB2()
    setupCamera(eye, gl)
    paintScreen()
}
function initNavigationEvents(_gl: BREPGLContext, eye: {pos: V3, focus: V3, up: V3, zoomFactor: number}, paintScreen: () => void) {
    const canvas: HTMLCanvasElement = $(_gl.canvas)
    let lastPos: V3 = V3.O
    //_gl.onmousedown.push((e) => {
    //	e.preventDefault()
    //	e.stopPropagation()
    //})
    //_gl.onmouseup.push((e) => {
    //	e.preventDefault()
    //	e.stopPropagation()
    //})
    canvas.addEventListener('mousemove', e => {
        const pagePos = V(e.pageX, e.pageY)
        const delta = lastPos.to(pagePos)
        //noinspection JSBitwiseOperatorUsage
        if (e.buttons & 4) {
            // pan
            const moveCamera = V(-delta.x * 2 / _gl.canvas.width, delta.y * 2 / _gl.canvas.height)
            const inverseProjectionMatrix = _gl.projectionMatrix.inversed()
            const worldMoveCamera = inverseProjectionMatrix.transformVector(moveCamera)
            eye.pos = eye.pos.plus(worldMoveCamera)
            eye.focus = eye.focus.plus(worldMoveCamera)
            setupCamera(eye, _gl)
            paintScreen()
        }
        // scene rotation
        //noinspection JSBitwiseOperatorUsage
        if (e.buttons & 2) {
            const rotateLR = -delta.x / 6.0 * DEG
            const rotateUD = -delta.y / 6.0 * DEG
            // rotate
            let matrix = M4.rotateLine(eye.focus, eye.up, rotateLR)
            //let horizontalRotationAxis = focus.minus(pos).cross(up)
            const horizontalRotationAxis = eye.up.cross(eye.pos.minus(eye.focus))
            matrix = matrix.times(M4.rotateLine(eye.focus, horizontalRotationAxis, rotateUD))
            eye.pos = matrix.transformPoint(eye.pos)
            eye.up = matrix.transformVector(eye.up)

            setupCamera(eye, _gl)
            paintScreen()
        }
        lastPos = pagePos
    })
    canvas.addEvent('mousewheel', function (e) {
        //console.log(e)
        eye.zoomFactor *= pow(0.9, -e.wheel)
        const targetPos = e.target.getPosition()
        const mouseCoords = {x: e.page.x - targetPos.x, y: e.page.y - targetPos.y}
        const moveCamera = V(mouseCoords.x * 2 / _gl.canvas.width - 1, -mouseCoords.y * 2 / _gl.canvas.height + 1, 0).times(1 - 1 / pow(0.9, -e.wheel))
        const inverseProjectionMatrix = _gl.projectionMatrix.inversed()
        const worldMoveCamera = inverseProjectionMatrix.transformVector(moveCamera)
        //console.log("moveCamera", moveCamera)
        //console.log("worldMoveCamera", worldMoveCamera)
        eye.pos = eye.pos.plus(worldMoveCamera)
        eye.focus = eye.focus.plus(worldMoveCamera)
        setupCamera(eye, _gl)
        paintScreen()
        e.preventDefault()
    })
}
function makeDottedLinePlane(count: int = 128) {
    const mesh = new Mesh().
        addIndexBuffer('LINES')
    const OXvertices = arrayFromFunction(count, i => new V3(i / count, 0, 0))
    mesh.vertices.pushAll(OXvertices)
    mesh.vertices.pushAll(M4.forSys(V3.Y, V3.O, V3.O, V3.X).transformedPoints(OXvertices))
    mesh.vertices.pushAll(M4.forSys(V3.X.negated(), V3.O, V3.O, new V3(1, 1, 0)).transformedPoints(OXvertices))
    mesh.vertices.pushAll(M4.forSys(V3.Y.negated(), V3.O, V3.O, V3.Y).transformedPoints(OXvertices))
    mesh.LINES = arrayFromFunction(count * 4, i => i - (i >= count * 2 ? 1 : 0))
    mesh.compile()
    return mesh
}
function initMeshes(_meshes: { [name: string]: Mesh }, _gl) {
    _gl.makeCurrent()
    _meshes.sphere1 = Mesh.sphere(2)
    _meshes.segment = Mesh.plane({startY: -0.5, height: 1, detailX: 128})
    _meshes.text = Mesh.plane()
    _meshes.vector = Mesh.rotation([V3.O, V(0, 0.05, 0), V(0.8, 0.05), V(0.8, 0.1), V(1, 0)], L3.X, TAU, 16, true)
    _meshes.pipe = Mesh.rotation(arrayFromFunction(128, i => new V3(i / 127, -0.5, 0)), L3.X, TAU, 8, true)
    _meshes.xyLinePlane = Mesh.plane()
    _meshes.xyDottedLinePlane = makeDottedLinePlane()
}
function initShaders(_gl) {
    _gl.makeCurrent()
    return {
        singleColor: Shader.create(vertexShaderBasic, fragmentShaderColor),
        multiColor: Shader.create(vertexShaderColor, fragmentShaderVaryingColor),
        singleColorHighlight: Shader.create(vertexShaderBasic, fragmentShaderColorHighlight),
        textureColor: Shader.create(vertexShaderTexture, fragmentShaderTextureColor),
        arc: Shader.create(vertexShaderRing, fragmentShaderColor),
        arc2: Shader.create(vertexShaderArc, fragmentShaderColor),
        ellipse3d: Shader.create(vertexShaderConic3d, fragmentShaderColor),
        generic3d: Shader.create(vertexShaderGeneric, fragmentShaderColor),
        bezier3d: Shader.create(vertexShaderBezier3d, fragmentShaderColor),
        bezier: Shader.create(vertexShaderBezier, fragmentShaderColor),
        lighting: Shader.create(vertexShaderLighting, fragmentShaderLighting),
        waves: Shader.create(vertexShaderWaves, fragmentShaderLighting),
    }
}


function setupCamera(_eye: typeof eye, _gl: LightGLContext) {
    const {pos, focus, up, zoomFactor} = _eye
    //console.log("pos", pos.$, "focus", focus.$, "up", up.$)
    _gl.matrixMode(_gl.PROJECTION)
    _gl.loadIdentity()
    //_gl.perspective(70, _gl.canvas.width / _gl.canvas.height, 0.1, 1000);
    const lr = _gl.canvas.width / 2 / zoomFactor
    const bt = _gl.canvas.height / 2 /zoomFactor
    _gl.ortho(-lr, lr, -bt, bt, -1e4, 1e4)
    _gl.lookAt(pos, focus, up)
    _gl.matrixMode(_gl.MODELVIEW)
    setupCameraListener && setupCameraListener(_eye)
}

//
//function test2() {
//    const ic: R2_R = (x, y) => sin(x+y)-cos(x*y)+1
//    const dids: R2_R = (x, y) => y * sin(x * y) + cos(x + y)
//    const didt: R2_R = (x, y) => x * sin(x * y) + cos(x + y)
//    const ic2: R2_R = (x, y) => (3 * x ** 2 - y ** 2) ** 2 * y ** 2 - (x ** 2 + y ** 2) ** 4
//    const di2ds: R2_R = (x, y) => 4* x* (9* x**2* y**2 - 3* y**4 - 2* (x**2 + y**2)**3)
//    const di2dt: R2_R = (x, y) => 2 * y * (-4 * (x ** 2 + y ** 2) ** 3 + (3 * x ** 2 - y ** 2) ** 2 + 2 * y ** 2 * (y ** 2 - 3 * x ** 2))
//    const start = V(-3.6339970071165784, 3.5625834844534974, 0) // curvePoint(ic, V(-4, 4))
//    assert(eq02(ic(start.x, start.y), 0.1))
//    const bounds = (s: number, t: number) => -5 <= s && s <= 5 && -5 <= t && t <= 5
//    //const curves =  Curve.breakDownIC(ic, -5, 5, -5, 5, 0.1, 0.1, 0.05, dids, didt)
//    const curves =  Curve.breakDownIC(ic2, {sMin: -5, sMax: 5, tMin: -5, tMax: 5}, 0.1, 0.1, 0.02, di2ds, di2dt)
//    //const curves =  Curve.breakDownIC(cassini(1, 1.02), -5, 5, -5, 5, 0.1, 0.1, 0.02)
//    //const curves = mkcurves(ic, start.x, start.y, 0.05, dids, didt, bounds)
//            .map(({points, tangents}, i) => {
//                const curve = new ImplicitCurve(ic, points, tangents)
//                return Edge.forCurveAndTs(curve.translate(5, 0, 0.1 * i))
//            })
//    //checkDerivate(s => ic(s, 0), s => dids(s, 0), -5, 5, 0)
//    //checkDerivate(t => ic(0, t), t => dids(0, t), -5, 5, 0)
//    console.log(curves.length)
//    return curves
//
//}
function cassini(a: number, c: number): (x: number, y: number) => number {
    return (x, y) => (x * x + y * y) * (x * x + y * y) - 2 * c * c * (x * x - y * y) - (a ** 4 - c ** 4)
}

/**
 * A function R² -> R with first and second derivatives.
 */
interface MathFunctionR2R {
    (s: number, t: number): number
    readonly x: R2_R
    readonly y: R2_R
    readonly xx?: R2_R
    readonly xy?: R2_R
    readonly yy?: R2_R
}
namespace MathFunctionR2R {
    export function forNerdamer(expression: nerdamer.ExpressionParam, args: [string, string] = ['x', 'y']): MathFunctionR2R {
        const ndf = nerdamer(expression)
        const ndfs = nerdamer.diff(ndf, args[0])
        const ndft = nerdamer.diff(ndf, args[1])
        const f = ndf.buildFunction(args) as any
        f.x = ndfs.buildFunction(args)
        f.y = ndft.buildFunction(args)
        f.xx = nerdamer.diff(ndfs, args[0]).buildFunction(args)
        f.xy = nerdamer.diff(ndfs, args[1]).buildFunction(args)
        f.yy = nerdamer.diff(ndft, args[1]).buildFunction(args)
        return f
    }

    export function nerdamerToR2_R(expression: nerdamer.Expression, args: [string, string] = ['x', 'y']) {
        return expression.buildFunction(args)
    }

    export function forFFxFy(f: R2_R, fx: R2_R, fy: R2_R): MathFunctionR2R {
        ;(f as any).x = fx
        ;(f as any).y = fy
        return f as any
    }
}
const cas2 = cassini(0.9, 1.02)
function HJKl() {
    //math.derivative('x^2', 'x')
    nerdamer.setFunction('cassini', 'acxy'.split(''), '(x^2 + y^2)^2 + 2 c^2 (x^2 - y^2) - (a^4 - c^4)')
    const cassini = nerdamer('(x^2 + y^2)^2 + 2 c^2 (x^2 - y^2) - (a^4 - c^4)')
    const ndf = nerdamer('cassini(1, 1, x, y)')
    const mf = MathFunctionR2R.forNerdamer(ndf)
    //const pf = cs.parametricFunction(), icc = ses.implicitFunction()
    //const start = V(-3.6339970071165784, 3.5625834844534974, 0) // curvePoint(ic, V(-4, 4))
    //assert(eq02(ic(start.x, start.y), 0.1))

    const bounds2={sMin: -2, sMax : 2, tMin : -2, tMax : 2}
    const {sMin, tMin, sMax, tMax} = bounds2
    const bounds = (s: number, t: number) => sMin <= s && s <= sMax && tMin <= t && t <= tMax
    //const curves =  Curve.breakDownIC(ic, -5, 5, -5, 5, 0.1, 0.1, 0.05, dids, didt)
    const {points, tangents}
            = followAlgorithm2d(mf, curvePointMF(mf, V(1, 0.5)), 0.05, bounds)
    // const edges = [Edge.forCurveAndTs(new ImplicitCurve(points, tangents).scale(10))]
    const edges =  Curve.breakDownIC(mf, bounds2, 0.1, 0.1, 0.1)
            .map(({points, tangents}, i) => Edge.forCurveAndTs(new ImplicitCurve(points, tangents, 1)).translate(0,0,i*0.1).scale(10))
    //const curves =  Curve.breakDownIC(cassini(1, 1.02), -5, 5, -5, 5, 0.1, 0.1, 0.02)
    //const curves = mkcurves(ic, start.x, start.y, 0.05, dids, didt, bounds)
    console.log(edges.length)
    return {edges: edges, points: (edges[0].curve as ImplicitCurve).points}

}
function HJK() {
    //return {}
    //const a = B2T.cone(12, 12).scale(0.05,0.2)
    //	.rotateZ(90*DEG)
    //	.rotateY(-90*DEG)
    //	.translate(2,0.2,0.7)
    //const b = B2T.sphere(1, undefined, PI)
    //const cone = new ConicSurface(
    //	V(2, 0.2, 1.1),
    //	V(0, 0.6, 0),
    //	V(0, 0, -2.4),
    //	V(-12, 0, 0))
    //const sphere = new SemiEllipsoidSurface(V3.O,V3.X,V3.Y,V(0, 0, -1))
    ////const curves = cone.isCurvesWithSurface(sphere)
    ////assert(sphere.containsCurve(curves[0]))
    //const c = a.minus(b).translate(0,3)
    //const pcs = new ProjectedCurveSurface(new BezierCurve(V(-0.1040625, -0.095, 1.2), V(-0.1040625, -0.030937500000000007, 1.2), V(-0.065, 0.010000000000000009, 1.2), V(0.047968750000000004, 0.010000000000000009, 1.2), 0, 1), V(0, 0, -1), 0, 1, -200, 200)
    const pcs = SemiCylinderSurface.UNIT
            .rotateZ(-40*DEG)
            .scale(0.5, 0.05, 4)
            .translate(0.5,0,-2)
            .flipped()
    const ses = new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z)
    const curves = ses.isCurvesWithSurface(pcs)
    return {
        //a,
        //b,
        //c,
        mesh: [pcs.toMesh(), ses.toMesh()],
        edges: curves.map(c => Edge.forCurveAndTs(c)),
        points: curves.flatMap(c => (c as ImplicitCurve).points)
    }

}