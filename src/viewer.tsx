import * as chroma from "chroma.ts"
import debounce from "debounce"
import * as prettier from "prettier"
import typescriptParser from "prettier/parser-typescript"
import {
  assert,
  DEG,
  emod,
  int,
  M4,
  round10,
  V,
  V3,
  Vector,
  VV,
} from "ts3dutils"
import { GL_COLOR, GL_COLOR_BLACK, Mesh } from "tsgl"
import * as ReactDOM from "react-dom"
import * as React from "react"
import { red, green, blue } from "@material-ui/core/colors"
import {
  Button,
  Checkbox,
  CssBaseline,
  FormControlLabel,
  makeStyles,
} from "@material-ui/core"

import {
  BRep,
  BRepGLContext,
  cameraChangeListeners,
  COLORS,
  CustomPlane,
  Edge,
  Eye,
  Face,
  FaceMesh,
  getMouseLine,
  ImplicitCurve,
  initNavigationEvents,
  L3,
  P3,
  PICurve,
  setupCamera,
} from "."
import * as ts3dutils from "ts3dutils"
import * as tsgl from "tsgl"
import * as brepts from "."
import { useCallback, useEffect, useRef, useState } from "react"

const bReps: BRep[] = []
const edgeViewerColors = ["darkorange", "darkgreen", "cyan"].map((c) =>
  chroma.css(c).gl(),
)
let bRepMeshes: (Mesh & {
  faceIndexes: Map<Face, { start: int; count: int }>
  TRIANGLES: int[]
  normals: V3[]
})[] = []
//bMesh: Mesh & {faceIndexes?: Map<Face, {start: int, count: int}>},
//cMesh: Mesh & {faceIndexes?: Map<Face, {start: int, count: int}>},
let edgesMesh: Mesh & {
  faceIndexes?: Map<Face, { start: int; count: int }>
  TRIANGLES: int[]
  normals: V3[]
  curve1: V3[]
  curve1colors: GL_COLOR[]
}
let faceMesh: FaceMesh & { tangents: V3[] }
let meshes: (Mesh & { TRIANGLES: int[]; normals: V3[] })[] = []
let hovering: {}
const edgeDebugPoints: V3[] = []
const edgeDebugLines: V3[] = []

const addMissing = (to: any, from: any) =>
  Object.keys(from).forEach(
    (key) => "Buffer" != key && !to[key] && (to[key] = from[key]),
  )

// tslint:disable-next-line:class-name
export class RenderObjects {
  a: BRep | undefined = undefined
  b: BRep | undefined = undefined
  c: BRep | undefined = undefined
  d: BRep | undefined = undefined
  face: Face[] = []
  edges: Edge[] = []
  wireframe: boolean = false
  normallines: boolean = false
  i: any = undefined
  hjk: any = undefined
  drPs: (V3 | { p: V3; info?: string; color?: string })[] = []
  drVs: { v: V3; anchor: V3; color?: GL_COLOR }[] = []
  drLines: V3[] = []
  mesh: (Mesh & { TRIANGLES: int[]; normals: V3[] })[] = []
  boxes: M4[] = []
  planes: CustomPlane[] = []
}

const renderObjectKeys = Object.keys(
  new RenderObjects(),
) as (keyof RenderObjects)[]

declare function INIT_HTML(): void

addMissing(window, ts3dutils)
addMissing(window, tsgl)
addMissing(window, brepts)
addMissing(window, new RenderObjects())
declare global {
  interface Window extends RenderObjects {}
}
const arrayLiteralType = <T extends string>(x: T[]): T[] => x
const g = window

function objectAssignConcatArray<T, U>(a: T, b: U): T & U {
  for (const key of Object.keys(b)) {
    if (Array.isArray((g as any)[key]) && Array.isArray((b as any)[key])) {
      ;(a as any)[key].push(...(b as any)[key])
    } else if (undefined !== (b as any)[key]) {
      ;(a as any)[key] = (b as any)[key]
    }
  }
  return a as any
}

function initBRep() {
  const htmlContext = INIT_HTML()

  const hash =
    window.location.search.substr(1) || window.location.hash.substr(1) || ""
  const command = decodeURIComponent(hash)
  console.log(command)
  const hashContext = new Function(
    `let ${renderObjectKeys.join(
      ",",
    )};${command};return{${renderObjectKeys.join(",")}}`,
  )() as RenderObjects

  // hashContext last, so i value in hash wins
  objectAssignConcatArray(g, htmlContext)
  objectAssignConcatArray(g, hashContext)
  console.log(htmlContext)

  // let gets: any = {a, b, c, d, mesh, edges, points, vectors}
  // g.hjk && Object.assign(g, HJK())
  arrayLiteralType(["a", "b", "c", "d"]).forEach((k) => {
    const bRep = g[k]
    if (bRep) {
      bReps.push(bRep)
    }
  })

  bRepMeshes = bReps.map((bRep) => bRep.toMesh())
  bRepMeshes.forEach((mesh) => {
    mesh.computeWireframeFromFlatTriangles("wireframe")
    mesh.computeNormalLines(0.1, "normallines")
    mesh.compile()
  })

  if (g.mesh) {
    console.log("mesh/es from GET", bRepMeshes)
    meshes.forEach((mesh) => {
      mesh.computeWireframeFromFlatTriangles("wireframe")
      mesh.computeNormalLines(0.1, "normallines")
      mesh.compile()
    })
  }

  if (g.edges) {
    console.log("edges from GET")
    edgesMesh = new Mesh()
      .addIndexBuffer("TRIANGLES")
      .addVertexBuffer("normals", "ts_Normal")
      .addVertexBuffer("curve1", "curve1")
      .addVertexBuffer("curve1colors", "curve1colors")
    g.edges.forEach((edge, edgeIndex) => {
      const points = edge.points()
      for (let i = 0; i < points.length - 1; i++) {
        const color =
          edgeViewerColors[(edgeIndex + (i % 2)) % edgeViewerColors.length]
        // const tangent = edge.tangentAt(i)
        // dMesh.curve1.push(points[i], points[i].plus(tangent.toLength(1)))
        edgesMesh.curve1.push(points[i], points[i + 1])
        edgesMesh.curve1colors.push(color, color)
      }
      edge.curve instanceof PICurve &&
        (edge.curve as PICurve).addToMesh(edgesMesh, 8, 0.02, 2)

      if (edge.curve.debugInfo) {
        const { points, lines } = edge.curve.debugInfo()
        points && edgeDebugPoints.push(...points)
        lines && edgeDebugLines.push(...lines)
      }
    })
    //dMesh.computeWireframeFromFlatTriangles()
    edgesMesh.compile()
  }
  if (g.face) {
    faceMesh = new Mesh()
      .addIndexBuffer("TRIANGLES")
      .addIndexBuffer("LINES")
      .addVertexBuffer("tangents", "tangents")
      .addVertexBuffer("normals", "ts_Normal")
    for (const face of g.face) {
      face.addToMesh(faceMesh)
      for (const edge of face.allEdges) {
        const ts = edge.curve.calcSegmentTs(
          edge.aT,
          edge.bT,
          edge.reversed,
          true,
        )
        for (const t of ts) {
          const p = edge.curve.at(t)
          faceMesh.tangents.push(p, p.plus(edge.tangentAt(t)))
        }
      }
    }
    faceMesh.compile()
  }

  g.drPs.push()
}

const brepMeshColors: chroma.Color[][] = [
  chroma.scale(["#ff297f", "#6636FF"]),
  chroma.scale(["#ffe93a", "#ff6e35"]),
  chroma.scale(["#1eff33", "#4960ff"]),
  chroma.scale(["#31fff8", "#2dff2a"]),
].map((scale) => scale.mode("lab").colors(20, "color"))
const brepMeshColorssGL = brepMeshColors.map((cs) => cs.map((c) => c.gl()))
const meshColorsGL: GL_COLOR[] = chroma.scale("GnBu").colors(16, "gl")

function viewerPaint(
  time: int,
  gl: BRepGLContext & WebGL2RenderingContext,
  eye: Eye,
  options: {
    paintMeshNormals: boolean
    paintWireframe: boolean
    paintCurveDebug: boolean
    renderElement: V3 | Vector | undefined
  },
) {
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
  gl.loadIdentity()

  if (options.renderElement?.constructor === V3) {
    gl.drawPoint(options.renderElement, undefined, 12 / eye.zoomFactor)
  } else if (options.renderElement?.constructor === Vector) {
    gl.drawPoint(options.renderElement.p3(), undefined, 12 / eye.zoomFactor)
  }

  //setupCamera(eye, gl)

  gl.drawVectors(g.drVs, 4 / eye.zoomFactor)

  g.drPs.forEach((info) =>
    gl.drawPoint(
      info instanceof V3 ? info : info.p,
      info instanceof V3 || !info.color
        ? chroma.css("#cc0000").gl()
        : chroma.css(info.color).gl(),
      6 / eye.zoomFactor,
    ),
  )
  drawPlanes.forEach((plane) =>
    gl.drawPlane(plane, plane.color, hovering == plane),
  )
  g.planes.forEach((plane) =>
    gl.drawPlane(plane, plane.color, hovering == plane),
  )

  g.boxes.forEach((m4) => gl.drawBox(m4))

  gl.shaders.lighting.uniforms({ camPos: eye.pos })
  for (let i = 0; i < bRepMeshes.length; i++) {
    const mesh = bRepMeshes[i]
    gl.pushMatrix()
    //gl.translate(30, 0, 0)
    gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
    options.paintWireframe &&
      mesh.indexBuffers.wireframe &&
      gl.shaders.singleColor
        .uniforms({ color: COLORS.TS_STROKE.gl() })
        .drawBuffers(mesh.vertexBuffers, mesh.indexBuffers.wireframe, gl.LINES)
    options.paintMeshNormals &&
      mesh.indexBuffers.normallines &&
      gl.shaders.singleColor
        .uniforms({ color: COLORS.TS_STROKE.gl() })
        .drawBuffers(
          mesh.vertexBuffers,
          mesh.indexBuffers.normallines,
          gl.LINES,
        )
    gl.shaders.singleColor
      .uniforms({ color: COLORS.TS_STROKE.gl() })
      .drawBuffers(mesh.vertexBuffers, mesh.indexBuffers.LINES, gl.LINES)
    gl.projectionMatrix.m[11] += 1 / (1 << 20)

    let faceIndex = bReps[i].faces.length
    while (faceIndex--) {
      const face = bReps[i].faces[faceIndex]
      const faceTriangleIndexes = mesh.faceIndexes.get(face)!
      gl.shaders.lighting
        .uniforms({
          color:
            hovering == face
              ? emod(emod(brepMeshColors, i), faceIndex).darker(2).gl()
              : emod(emod(brepMeshColorssGL, i), faceIndex),
        })
        .draw(
          mesh,
          gl.TRIANGLES,
          faceTriangleIndexes.start,
          faceTriangleIndexes.count,
        )
    }

    gl.popMatrix()
  }

  if (faceMesh) {
    gl.shaders.singleColor
      .uniforms({ color: chroma.css("red").gl() })
      .drawBuffers(
        { ts_Vertex: faceMesh.vertexBuffers.tangents },
        undefined,
        gl.LINES,
      )
    gl.shaders.lighting
      .uniforms({
        color: chroma.css("red").gl(),
        camPos: eye.pos,
      })
      .draw(faceMesh)
  }

  for (let i = 0; i < meshes.length; i++) {
    const mesh = meshes[i]
    gl.pushMatrix()
    gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
    options.paintWireframe &&
      mesh.indexBuffers.wireframe &&
      gl.shaders.singleColor
        .uniforms({ color: COLORS.TS_STROKE.gl() })
        .drawBuffers(mesh.vertexBuffers, mesh.indexBuffers.wireframe, gl.LINES)
    options.paintMeshNormals &&
      mesh.indexBuffers.normallines &&
      gl.shaders.singleColor
        .uniforms({ color: COLORS.TS_STROKE.gl() })
        .drawBuffers(
          mesh.vertexBuffers,
          mesh.indexBuffers.normallines,
          gl.LINES,
        )
    gl.projectionMatrix.m[11] += 1 / (1 << 20)
    mesh.TRIANGLES &&
      gl.shaders.lighting
        .uniforms({
          color: emod(meshColorsGL, i),
          camPos: eye.pos,
        })
        .draw(mesh)
    gl.popMatrix()
  }

  if (hovering instanceof Edge) {
    gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
    gl.drawEdge(hovering, GL_COLOR_BLACK, 2 / eye.zoomFactor)
    gl.projectionMatrix.m[11] += 1 / (1 << 20)
  }
  g.edges.forEach((e, i) =>
    gl.drawEdge(e, emod(edgeViewerColors, i), 3 / eye.zoomFactor),
  )

  if (options.paintCurveDebug) {
    gl.begin(gl.LINES)
    gl.color("red")
    edgeDebugLines.forEach((x) => gl.vertex(x))
    gl.end()
    edgeDebugPoints.forEach((p) =>
      gl.drawPoint(p, chroma.css("red").gl(), 6 / eye.zoomFactor),
    )
  }
  if (0 !== g.drLines.length) {
    gl.begin(gl.LINES)
    g.drLines.forEach((x) => {
      gl.color((x as any).color || "red")
      gl.vertex(x)
    })
    gl.end()
  }
}

function getHovering(
  mouseLine: L3,
  faces: Face[],
  planes: CustomPlane[] | undefined,
  points: V3[],
  edges: Edge[],
  mindist: number,
  ...consider: (
    | "faces"
    | "planes"
    | "sketchElements"
    | "points"
    | "edges"
    | "features"
  )[]
): any {
  let hoverHighlight = null,
    nearest = Infinity

  const checkFeatures = consider.includes("features")
  assert(!checkFeatures || !consider.includes("faces"))

  function checkEl(el: any, distance: number) {
    if (distance < nearest) {
      nearest = distance
      hoverHighlight = el
    }
  }

  if (faces && (consider.includes("faces") || consider.includes("features"))) {
    for (const face of faces) {
      checkEl(
        checkFeatures ? face.info.feature : face,
        face.intersectsLine(mouseLine),
      )
    }
  }
  if (planes && consider.includes("planes")) {
    for (const plane of planes) {
      checkEl(plane, plane.distanceTo(mouseLine, mindist))
    }
  }
  if (consider.includes("points")) {
    for (const p of points) {
      const t = mouseLine.pointT(p)
      if (mouseLine.at(t).distanceTo(p) < mindist * 1.2) {
        checkEl(p, t - 0.1)
      }
    }
  }
  if (consider.includes("edges")) {
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

// @ts-expect-error
function initInfoEvents(paintScreen: () => {}, gl: BRepGLContext) {
  gl.canvas.addEventListener("mousemove", function (e) {
    const mouseLine = getMouseLine({ x: e.clientX, y: e.clientY }, gl)
    const faces = bReps.flatMap((b2) => b2 && b2.faces)
    const testEdges: Edge[] = [
      ...bReps.flatMap((b2) =>
        Array.from<Edge>(b2.buildAdjacencies().edgeFaces.keys()),
      ),
      ...g.edges,
    ]
    hovering = getHovering(
      mouseLine,
      faces,
      undefined,
      [],
      testEdges,
      0.1,
      "faces",
      "edges",
    )
    //let html = '', pp
    //if (hovering instanceof Edge) {
    //	pp = V(e.clientX, e.clientY)
    //	defaultRoundFunction = x => round10(x, -3)
    //	html = hovering.toString(x => round10(x, -3)) + ' length=' +
    // hovering.length().toFixed(3) } else if (hovering instanceof Face) { pp =
    // V(e.clientX, e.clientY) defaultRoundFunction = x => round10(x, -3) let
    // area try { area = hovering.calcArea() } catch (e) {} html = `face
    // surface=${hovering.surface.constructor.name}
    // edges=${hovering.contour.length} area=${area}` } if (pp) { //const pSC =
    // gl.projectionMatrix.times(gl.modelViewMatrix).transformPoint(pp) //const
    // x = (pSC.x * 0.5 + 0.5) * window.innerWidth, y = (-pSC.y * 0.5 + 0.5) *
    // window.innerHeight tooltipShow(html, pp.x, pp.y) } else { tooltipHide()
    // }
    paintScreen()
  })
}

//var sketchPlane = new CustomPlane(V3.X, V3(1, 0, -1).unit(), V3.Y, -500, 500,
// -500, 500, 0xff00ff);
const drawPlanes = [
  new CustomPlane(V3.O, V3.Y, V3.Z, "planeYZ", chroma.color(0xff0000).gl()),
  new CustomPlane(V3.O, V3.X, V3.Z, "planeZX", chroma.color(0x00ff00).gl()),
  new CustomPlane(V3.O, V3.X, V3.Y, "planeXY", chroma.color(0x0000ff).gl()),
  //	sketchPlane
]
let paintScreen: () => void
declare var BREPTS_ROOT: string

const omap = <U, V, O extends Record<string, U>>(
  obj: O,
  f: (val: U, key: string) => V,
): { [k in keyof O]: V } =>
  Object.keys(obj).reduce(
    (prev, curr) => ((prev[curr] = f(obj[curr], curr)), prev),
    {} as any,
  ) as any

const useStyles = makeStyles((theme) => ({
  tipWrap: {
    "& > *": {
      marginRight: theme.spacing(1),
    },
  },
  max: {
    zIndex: -1000,
    display: "flex",
    "& > *": {
      zIndex: 1000,
    },
  },
  ...omap({ red, green, blue }, (color: Record<string | number, string>) => ({
    backgroundColor: color[500],
    color: "#fff",
    "&:hover": {
      backgroundColor: color[700],
    },
  })),
}))

const useInfoStyles = makeStyles({
  vspan: {
    textDecorationLine: "underline",
    textDecorationStyle: "dashed",
    "&:hover": {
      backgroundColor: "#eeeeee",
    },
  },
})

function Info({
  renderObjects,
  style,
  onRenderElement,
}: {
  renderObjects: Partial<RenderObjects>
  style?: React.CSSProperties
  onRenderElement: (e: V3 | Vector) => void
}) {
  const spanHover = useCallback(
    (e: React.MouseEvent<HTMLSpanElement>) => {
      onRenderElement(
        new Function("V", "VV", "return " + e.currentTarget.dataset.what!)(
          V,
          VV,
        ),
      )
    },
    [onRenderElement],
  )

  const classes = useInfoStyles()

  return (
    <div style={style}>
      {renderObjectKeys
        .map((k) => [k, renderObjects[k]])
        .filter(([k, v]) => "i" !== k && v && (!Array.isArray(v) || v.length))
        .map(([k, v]): [string, any[], string?] =>
          Array.isArray(v) ? [k, v, k] : [k, [v]],
        )
        .flatMap(([k, vs, name]) =>
          vs.map((v, vi) => {
            let source = v?.toSource()
            try {
              source = prettier.format(source, {
                parser: "typescript",
                semi: false,
                trailingComma: "all",
                plugins: [typescriptParser],
              })
            } catch (e) {}
            const FLOAT_REGEX = /[0-9+-.e]+/
            const V_OR_VV_REGEX = new RegExp(
              "VV?\\(\\s*" +
                FLOAT_REGEX.source +
                "(?:\\s*,\\s*" +
                FLOAT_REGEX.source +
                "){2,3}\\s*\\)" +
                "|V3\\.\\w+",
              "g",
            )
            const tokenized: (string | RegExpExecArray)[] = []
            let result: RegExpExecArray | null
            let lastPos = 0
            while ((result = V_OR_VV_REGEX.exec(source)) !== null) {
              tokenized.push(source.substr(lastPos, result.index - lastPos))
              tokenized.push(result)
              lastPos = result.index + result[0].length
            }
            tokenized.push(source.substr(lastPos))
            return (
              <div key={k + vi}>
                {name ?? `${k}[${vi}]`}
                <pre>
                  <code>
                    {tokenized.map((x, xi) =>
                      "string" === typeof x ? (
                        x
                      ) : (
                        <span
                          className={classes.vspan}
                          data-what={x[0]}
                          onMouseEnter={spanHover}
                          key={xi}
                        >
                          {x[0]}
                        </span>
                      ),
                    )}
                  </code>
                </pre>
              </div>
            )
          }),
        )}
    </div>
  )
}

function Nav({
  style,
  hash,
  testName,
}: {
  style?: React.CSSProperties
  hash: string
  testName: string
}) {
  const [state, setState] = useState({} as Record<string, string>)
  useEffect(() => {
    fetch("/tests/results/index.json")
      .then((r) => r.json())
      .then(setState)
  }, [])
  return (
    <div style={style}>
      {Object.entries(state).map(([name, url]) => (
        <div key={name}>
          <a
            style={{ fontWeight: name === testName ? "bold" : "normal" }}
            href={url + hash}
          >
            {name}
          </a>
        </div>
      ))}
    </div>
  )
}

const Viewer = () => {
  const [state, setState] = useState({
    paintMeshNormals: false,
    paintWireframe: false,
    paintCurveDebug: false,
    renderElement: undefined as V3 | Vector | undefined,
  })
  const stateRef = useRef(state)
  stateRef.current = state

  const [renderObjects, setRenderObjects] = useState<RenderObjects>()

  const [hash, setHash] = useState(window.location.hash)

  const canvasRef = useRef<HTMLCanvasElement>(null)
  const glRef = useRef<BRepGLContext>()
  const eye = useRef({
    pos: V(1000, 1000, 1000),
    focus: V3.O,
    up: V3.Z,
    zoomFactor: 1,
  }).current

  useEffect(() => {
    paintScreen = () =>
      requestAnimationFrame((t) =>
        viewerPaint(t, glRef.current!, eye, stateRef.current),
      )
    // B2T.defaultFont = await B2T.loadFont(BREPTS_ROOT +
    // '/fonts/FiraSansMedium.woff')
    window.onerror = function (errorMsg, url, lineNumber, column, errorObj) {
      console.log(errorMsg, url, lineNumber, column, errorObj)
    }
    const gl = BRepGLContext.create({
      canvas: canvasRef.current!,
    })
    glRef.current = gl
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

    initNavigationEvents(gl, eye, paintScreen, (_) => V3.O)
    window.onpopstate = function () {
      const hash =
        window.location.search.substr(1) || window.location.hash.substr(1) || ""
      setHash("#" + hash)
      const command = decodeURIComponent(hash)
      const hashContext = new Function(
        `let ${renderObjectKeys.join(
          ",",
        )};${command};return{${renderObjectKeys.join(",")}}`,
      )() as RenderObjects

      Object.assign(eye, hashContext.i)
      setupCamera(eye, gl, true)
      paintScreen()
    }

    cameraChangeListeners.push(
      debounce(function (eye) {
        const round = (x: number) => round10(x, -3)
        const roundedEye = {
          pos: eye.pos.map(round),
          focus: eye.focus.map(round),
          up: eye.up.map(round),
          zoomFactor: round(eye.zoomFactor),
        }
        const iSource =
          "i=" + roundedEye.toSource().replace(/[\n\r\s]+|^\(|\)$/g, "")
        const hash = window.location.hash.substr(1) || iSource
        const result = hash.match(/i={[^}]*}/)
          ? hash.replace(/i={[^}]*}/, iSource)
          : hash + ";" + iSource
        window.history.pushState(undefined, "", "#" + result)
        setHash("#" + result)
      }, 500),
    )
    // initInfoEvents(paintScreen, g l)
    //initToolTips() // hide tooltip on mouseover
    //initPointInfoEvents()
    initBRep()
    setRenderObjects(window)
    document.title = window.TEST_NAME
    Object.assign(eye, g.i)
    setupCamera(eye, gl)
    paintScreen()
  }, [])

  useEffect(() => paintScreen(), [state])

  function alignX(dir: number) {
    eye.pos = eye.focus.plus(V(100 * dir, 0, 0))
    eye.up = V3.Z
    setupCamera(eye, glRef.current!)
    paintScreen()
  }

  function alignY(dir: number) {
    eye.pos = eye.focus.plus(V(0, 100 * dir, 0))
    eye.up = V3.Z
    setupCamera(eye, glRef.current!)
    paintScreen()
  }

  function alignZ(dir: number) {
    eye.pos = eye.focus.plus(V(0, 0, 100 * dir))
    eye.up = eye.pos.cross(V3.X).unit()
    setupCamera(eye, glRef.current!)
    paintScreen()
  }

  function rot(angleInDeg: number) {
    eye.up = M4.rotateLine(
      eye.pos,
      eye.pos.to(eye.focus),
      angleInDeg * DEG,
    ).transformVector(eye.up)
    setupCamera(eye, glRef.current!)
    paintScreen()
  }

  const classes = useStyles()

  return (
    <CssBaseline>
      <div style={{}} className={classes.max}>
        <canvas ref={canvasRef} style={{}} />
        <div
          id="tip-wrap"
          className={classes.tipWrap}
          style={{ display: "flex", flexDirection: "column" }}
        >
          <Button
            variant="contained"
            className={classes.red}
            onClick={() => alignX(1)}
          >
            +X
          </Button>
          <Button
            variant="contained"
            className={classes.red}
            onClick={() => alignX(-1)}
          >
            -X
          </Button>
          <Button
            variant="contained"
            className={classes.green}
            onClick={() => alignY(1)}
          >
            +Y
          </Button>
          <Button
            variant="contained"
            className={classes.green}
            onClick={() => alignY(-1)}
          >
            -Y
          </Button>
          <Button
            variant="contained"
            className={classes.blue}
            onClick={() => alignZ(1)}
          >
            +Z
          </Button>
          <Button
            variant="contained"
            className={classes.blue}
            onClick={() => alignZ(-1)}
          >
            -Z
          </Button>
          <Button variant="contained" onClick={() => rot(22.5)}>
            ↺
          </Button>
          <Button variant="contained" onClick={() => rot(-22.5)}>
            ↻
          </Button>

          <FormControlLabel
            control={
              <Checkbox
                checked={state.paintMeshNormals}
                onChange={() =>
                  setState((s) => ({
                    ...s,
                    paintMeshNormals: !s.paintMeshNormals,
                  }))
                }
                name="checkedB"
                color="primary"
              />
            }
            label="mesh normals"
          />
          <FormControlLabel
            control={
              <Checkbox
                checked={state.paintWireframe}
                onChange={() =>
                  setState((s) => ({
                    ...s,
                    paintWireframe: !s.paintWireframe,
                  }))
                }
                name="checkedB"
                color="primary"
              />
            }
            label="wireframe"
          />
          <FormControlLabel
            control={
              <Checkbox
                checked={state.paintCurveDebug}
                onChange={() =>
                  setState((s) => ({
                    ...s,
                    paintCurveDebug: !s.paintCurveDebug,
                  }))
                }
                name="checkedB"
                color="primary"
              />
            }
            label="curve debug"
          />
        </div>
        <Nav testName={window.TEST_NAME} hash={hash} />
        {renderObjects && (
          <Info
            renderObjects={renderObjects}
            onRenderElement={(e) =>
              setState((s) => ({
                ...s,
                renderElement: e,
              }))
            }
            style={{ marginLeft: "auto", width: 512 }}
          />
        )}
      </div>
    </CssBaseline>
  )
}

window.onload = function () {
  ReactDOM.render(<Viewer />, document.getElementById("react-root"))
}
