import * as chroma from "chroma.ts"
import classnames from "classnames"
import * as hljs from "highlight.js"
import React, { useCallback, useEffect, useRef, useState } from "react"
import ReactDOM from "react-dom"
import { assertf, clamp, DEG, emod, int, round10, TAU, V, V3 } from "ts3dutils"
import { TSGLContext } from "tsgl"

import {
  BREPGLContext,
  initMeshes,
  initNavigationEvents,
  initShaders,
  setupCamera,
} from "./BREPGLContext"
import { B2T, BRep, CustomPlane, edgeNgon, P3, round } from "."

const fakeB2Mesh = (false as true) && ({} as BRep).toMesh()
type B2Mesh = typeof fakeB2Mesh

function fixFunctionIndentation(src: string): string {
  const indent = /([ \t]*).*$/g.exec(src)![1]
  return src.replace(new RegExp("^" + indent, "gm"), "")
}

type DemoProps = {
  width: string
  height: string
  id: string
  f(...args: any[]): BRep | BRep[]
  args: DemoArg[]
}
type DemoDesc = {
  id: string
  f(...args: any[]): BRep | BRep[]
  args: DemoArg[]
  canvas: HTMLCanvasElement
  gl: BREPGLContext
  eye: { pos: V3; focus: V3; up: V3; zoomFactor: 8 }
  b2s: BRep[]
  meshes: (B2Mesh & { lines?: int[] })[]
}

function Demo({ id, f, args, height, width, ...props }: DemoProps) {
  const [state, setState] = useState(() => {
    args.forEach((arg) => (arg.value = "" + arg.def))
    return {
      showingSource: false,
      demo: ({ id, f, args, height, width } as any) as DemoDesc,
    }
  })

  const container = useRef<HTMLDivElement>(null)
  const canvas = useRef<HTMLCanvasElement>(null)
  const sourceContainer = useRef<HTMLDivElement>(null)

  const onChange = useCallback(() => {
    update(state.demo)
    setState((s) => ({ ...s }))
  }, [state, setState])

  const onClickSource = useCallback(() => {
    setState((s) => ({ ...s, showingSource: !s.showingSource }))
  }, [setState])

  useEffect(() => {
    setupDemo(canvas.current!, state.demo, container.current!)
    hljs.highlightBlock(sourceContainer.current!)
  }, [])

  const demo = state.demo
  const info =
    demo.b2s &&
    "faces: " +
      demo.b2s.map((b2) => b2.faces.length).join("/") +
      " edges: " +
      demo.b2s
        .map((b2) => (b2.edgeFaces && b2.edgeFaces.size) || "?")
        .join("/") +
      " triangles: " +
      demo.meshes.map((m) => (m ? m.TRIANGLES.length / 3 : 0)).join("/")
  return (
    <div {...props} style={{ width }} className={"democontainer"}>
      <div
        className="canvascontainer"
        ref={container}
        style={{ width: "100%", height }}
      >
        <canvas style={{ height }} ref={canvas} />
        {demo.args.map((arg) => (
          <InputComponent
            key={arg.name}
            label={arg.name}
            value={arg.value!}
            step={(arg as any).step}
            change={(e) => {
              arg.value = e
              onChange()
            }}
          />
        ))}
        <span className="info">{info}</span>
        <span className="navinfo">
          pan: drag-mmb | rotate: drag-rmb | zoom: scroll
        </span>
        <a className="sourcelink" onClick={onClickSource}>
          {state.showingSource ? "hide source" : "show source"}
        </a>
      </div>
      <code
        className={classnames(
          "src",
          !state.showingSource && "hide",
          "language-typescript",
        )}
        ref={sourceContainer}
      >
        {fixFunctionIndentation(f.toSource())}
      </code>
    </div>
  )
}

export type DemoArg =
  | {
      name: string
      def: string
      type: "string"
      fix(x: string): string
      value?: string
    }
  | {
      step?: number
      name: string
      def: number
      type: "int" | "angle" | "number"
      fix(x: number): number
      value?: string
    }

function InputComponent({
  step,
  value,
  change,
  label,
  ...atts
}: {
  step: number
  label: string
  value: string
  change(x: string): void
}) {
  const foo = (delta: number) =>
    change("" + round10(parseFloat(value) + delta, -6))

  const inputRef = useRef<HTMLInputElement>(null)

  const onWheel = (e: WheelEvent) => {
    e.preventDefault()
    const delta = (e.shiftKey ? 0.1 : 1) * Math.sign(-e.deltaY) * step
    foo(delta)
  }

  useEffect(() => {
    const input = inputRef.current!
    input.addEventListener("wheel", onWheel)
    return () => {
      input.removeEventListener("wheel", onWheel)
    }
  }, [onWheel])

  return (
    <span className="incont">
      <input
        {...atts}
        value={value}
        onChange={(e) => change(e.target.value)}
        className={classnames(step && "scrollable")}
        step={step}
        type={step ? "number" : "text"}
        ref={inputRef}
      />
      <label>{label}</label>
    </span>
  )
}

export async function demoMain() {
  await B2T.loadFonts()
}

function setupDemo(
  canvas: HTMLCanvasElement,
  demo: DemoDesc,
  container: HTMLDivElement,
) {
  canvas.width = container.clientWidth
  canvas.height = container.clientHeight
  window.addEventListener("resize", () => {
    canvas.width = container.clientWidth
    canvas.height = container.clientHeight
    gl.viewport(0, 0, canvas.width, canvas.height)
    setupCamera(demo.eye, gl)
  })
  const gl = (demo.gl = BREPGLContext.create(TSGLContext.create({ canvas })))
  gl.clearColor(1.0, 1.0, 1.0, 0.0)
  gl.enable(gl.BLEND)
  gl.enable(gl.DEPTH_TEST)
  gl.enable(gl.CULL_FACE)
  canvas.oncontextmenu = () => false
  gl.depthFunc(gl.LEQUAL)
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA)

  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
  gl.loadIdentity()
  gl.scale(10, 10, 10)

  gl.loadIdentity()

  gl.shaders = initShaders(gl)
  initMeshes((gl.meshes = {}), gl)
  demo.eye = { pos: V(10, 10, 100), focus: V(5, 5, 0), up: V3.Y, zoomFactor: 8 }
  initNavigationEvents(gl, demo.eye, () => paintDemo(demo))
  //initInfoEvents()
  //initPointInfoEvents()
  setupCamera(demo.eye, gl)

  // container.adopt(
  // demo.srcLink = new MooEl('a.sourcelink', {text: 'show source', href: '#'})
  // 	.addEvent('click', e => {
  // 		const showing = demo.srcLink.get('text') == 'hide source'
  // 		demo.srcContainer.setStyle('display', showing ? 'none' : 'block')
  // 		demo.srcLink.set('text', showing ? 'show source' : 'hide source')
  // 		return false
  // 	}),
  // )
  // hljs.highlightBlock(demo.srcContainer)
  update(demo)
}

const meshColorss = [
  chroma.scale(["#ffa621", "#ffd026"]).mode("lab").colors(20, "gl"),
  chroma.scale(["#ff297f", "#6636FF"]).mode("lab").colors(20, "gl"),
  chroma.scale(["#19ff66", "#1ffff2"]).mode("lab").colors(20, "gl"),
]
const demoPlanes = [
  new CustomPlane(
    V3.O,
    V3.Y,
    V3.Z,
    "planeYZ",
    chroma.color("red").gl(),
    -5,
    5,
    -5,
    5,
  ),
  new CustomPlane(
    V3.O,
    V3.X,
    V3.Z,
    "planeZX",
    chroma.color("green").gl(),
    -5,
    5,
    -5,
    5,
  ),
  new CustomPlane(
    V3.O,
    V3.X,
    V3.Y,
    "planeXY",
    chroma.color("blue").gl(),
    -5,
    5,
    -5,
    5,
  ),
]
const hovering: any = undefined

function paintDemo(demo: DemoDesc) {
  const gl = demo.gl
  gl.makeCurrent()
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
  gl.loadIdentity()
  demoPlanes.forEach((plane) =>
    gl.drawPlane(plane, chroma.color(plane.color).gl(), false),
  )
  // gl.drawVectors(g.Vs)
  if (!demo.meshes) return
  //viewerGL.scale(100, 100, 100)

  for (let i = 0; i < demo.meshes.length; i++) {
    const mesh = demo.meshes[i],
      b2 = demo.b2s[i]
    gl.pushMatrix()
    //viewerGL.translate(30, 0, 0)
    gl.projectionMatrix.m[11] -= 1 / (1 << 22) // prevent Z-fighting
    mesh.lines &&
      gl.shaders.singleColor
        .uniforms({ color: chroma.color("#bfbfbf").gl() })
        .draw(mesh, gl.LINES)
    gl.projectionMatrix.m[11] += 1 / (1 << 22)

    let faceIndex = b2.faces.length
    while (faceIndex--) {
      const face = b2.faces[faceIndex]
      const faceTriangleIndexes = mesh.faceIndexes.get(face)!
      gl.shaders.lighting
        .uniforms({
          color:
            hovering == face
              ? chroma.color("purple").gl()
              : emod(emod(meshColorss, i), faceIndex),
        })
        .draw(
          mesh,
          gl.TRIANGLES,
          faceTriangleIndexes.start,
          faceTriangleIndexes.count,
        )
      //shaders.singleColor.uniforms({
      //color: hexIntToGLColor(0x0000ff)
      //}).draw(brepMesh, 'LINES')
    }

    gl.popMatrix()
  }
}

function update(demo: DemoDesc) {
  const fixedParams = demo.args.map((arg) => {
    switch (arg.type) {
      case "number":
      case "int":
        return parseFloat(arg.value!)
      case "angle":
        return parseFloat(arg.value!) * DEG
    }
    return arg.value
  })
  //try {
  const b2s = demo.f.apply(undefined, fixedParams)
  demo.b2s = b2s instanceof Array ? b2s : [b2s]
  demo.b2s.forEach((b2) => b2.buildAdjacencies())
  //} catch (e) {
  //	console.log(e.message)
  //}
  demo.gl.makeCurrent()
  demo.meshes = demo.b2s.flatMap((b2, i) => {
    try {
      const m = b2.toMesh()
      assertf(() => m.faceIndexes.size == b2.faces.length)
      return [b2.toMesh().compile()]
    } catch (e) {
      console.error(`Error creating mesh from brep ${i}`, e, b2.toSource())
      return []
    }
  })
  paintDemo(demo)
}

window.onload = async () => {
  await B2T.loadFonts()
  ReactDOM.render(<Body />, document.getElementById("react-root"))
}

const Body = () => (
  <div>
    <h1>BREP.TS</h1>
    <p>
      This library describes volumes using a boundary representation model.
      Basically like triangle meshes, except faces can be any shape and are not
      necessarily planar. This allows for fast operations while maintaining a
      high degree of accuracy.
    </p>

    <p>
      Once you have two volumes, you can combine them using boolean operations.
      For instance:
    </p>
    <h4>box - sphere</h4>
    <Demo
      id="boxminussphere"
      width="100%"
      height="500px"
      f={function (sphereRadius, sphereZ) {
        const box = B2T.box(10, 10, 3)
        const sphere = B2T.sphere(sphereRadius).translate(5, 5, sphereZ)
        const result = box.minus(sphere)
        return [box, sphere, result.translate(12)]
      }}
      args={[
        {
          name: "sphere radius",
          type: "number",
          fix: (val) => clamp(val, 0.1, 1000),
          def: 2.5,
          step: 0.5,
        },
        {
          name: "sphere height",
          type: "number",
          fix: (val) => clamp(val, 0.1, 1000),
          def: 2,
          step: 0.5,
        },
      ]}
    />
    <h3>Functionality this library implements</h3>
    <ul>
      <li>
        Parametric curves: lines, ellipses, parabolas, hyperbolas,
        quadratic/cubic beziers, intersection curves between parametric and
        implicit surfaces.
      </li>
      <li>
        Surfaces: planes, cylinders, spheres, projected beziers. Functionality
        includes:
        <ul>
          <li>surface/surface intersections</li>
          <li>testing if two surfaces are coplanar</li>
          <li>testing if a surface contains a point</li>
          <li>testing if a surface contains a curve</li>
        </ul>
      </li>
      <li>Edge: segment of a curve</li>
      <li>Faces: test line intersection/point containment</li>
      <li>
        BREP volumes: intersection/union/subtraction, conversion to triangle
        mesh
      </li>
    </ul>
    <h3>Generator function examples</h3>
    <h4>
      cylinder(radius: number = 1, height: number = 1, rads: number = TAU)
    </h4>
    <Demo
      id="cyl"
      width="50%"
      height="400px"
      f={B2T.cylinder}
      args={[
        {
          name: "radius",
          type: "number",
          fix: (val) => clamp(val, 0.1, 1000),
          def: 10,
          step: 1,
        },
        {
          name: "height",
          type: "number",
          fix: (val) => clamp(val, 0.1, 1000),
          def: 10,
          step: 1,
        },
        {
          name: "degrees",
          type: "angle",
          fix: (val) => clamp(val, 0.1, TAU),
          def: 360,
          step: 10,
        },
      ]}
    />
    <h4>
      text(text: string, size: number, depth: number = 1, font: opentypejs.Font
      = defaultFont)
    </h4>
    <Demo
      id="text"
      width="50%"
      height="400px"
      f={B2T.text}
      args={[
        { name: "text", type: "string", fix: (val) => val, def: "foo" },
        {
          name: "size",
          type: "number",
          fix: (val) => clamp(val, 0.1, 100),
          def: 10,
          step: 5,
        },
        {
          name: "depth",
          type: "number",
          fix: (val) => clamp(val, 0.1, 100),
          def: 10,
          step: 1,
        },
      ]}
    />
    <h3>Advanced example: cylinder - extruded rounded edges</h3>
    <Demo
      id="test1"
      width="80%"
      height="400px"
      f={function (
        /* number */ outsideRadius,
        /* int */ n,
        /* number */ insideRadius,
        /* number */ cornerRadius,
        /* number */ depth,
      ) {
        const cylinder = B2T.cylinder(outsideRadius, 10)

        // create an n-gon centered on the XY-plane with corners insideRadius away from the origin
        // the ngon is counter-clockwise (CCW) when viewed from +Z
        const ngon = edgeNgon(n, insideRadius)

        // round the corners of the ngon with radius cornerRadius
        const ngonRounded = round(ngon, cornerRadius)

        // create a hole punch by extruding ngonRounded in direction +Z
        // Like triangles in opengl, faces face the direction they are CCW when viewed from.
        // Because we extrude in the direction the outline faces, the resulting is "inside-out",
        // with faces facing inwards. As we want to remove the 'punch' volume, this is not a
        // problem. Otherwise we could call volume.flipped() or extrude in -Z and volume.translate(...)
        // into the correct position.
        const punch = B2T.extrudeEdges(
          ngonRounded,
          P3.XY,
          V(0, 0, 10),
        ).translate(0, 0, depth)

        // punch is already inside-out, so use volume.and(...) instead of volume.minus(...)
        const result = cylinder.and(punch)

        return [cylinder, punch, result.translate(2 + outsideRadius * 2)]
      }}
      args={[
        {
          name: "outside radius",
          type: "number",
          fix: (val) => clamp(val, 0.1, 10),
          def: 5,
          step: 1,
        },
        {
          name: "corner count",
          type: "int",
          fix: (val) => clamp(val, 3, 21),
          def: 3,
          step: 1,
        },
        {
          name: "ngon radius",
          type: "number",
          fix: (val) => clamp(val, 0.1, 12),
          def: 4,
          step: 1,
        },
        {
          name: "corner radius",
          type: "number",
          fix: (val) => clamp(val, 0.1, 5),
          def: 0.5,
          step: 0.1,
        },
        {
          name: "depth",
          type: "number",
          fix: (val) => clamp(val, -9, 9),
          def: 0,
          step: 1,
        },
      ]}
    />
  </div>
)
