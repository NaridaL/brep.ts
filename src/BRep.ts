import earcut from "earcut"
import { JavaMap, JavaSet, JavaSet as CustomSet, Pair } from "javasetmap.ts"
import nerdamer from "nerdamer"
import {
  AABB,
  assert,
  assertf,
  assertInst,
  assertNever,
  assertNumbers,
  assertVectors,
  concatenated,
  eq,
  eq0,
  getLast,
  gt,
  indexWithMax,
  int,
  lt,
  M4,
  mapFilter,
  mapPush,
  newtonIterate2d,
  newtonIterateWithDerivative,
  NLA_DEBUG,
  NLA_PRECISION,
  SCE,
  snap,
  snap0,
  sum,
  TAU,
  Transformable,
  V,
  V3,
  withMax,
} from "ts3dutils"
import { Mesh } from "tsgl"

import {
  AABB2,
  Curve,
  curvePointMF,
  Edge,
  Face,
  FaceInfoFactory,
  L3,
  P3,
  ParametricSurface,
  PlaneFace,
  PointVsFace,
  R2,
  R2_R,
  Surface,
  uvInAABB2,
  createEdge,
} from "."

import { abs, sign, sqrt } from "./math"

export const EPS = 1e-5

let globalId = 0

export function getGlobalId() {
  return globalId++
}

export function addLikeSurfaceFaces(
  likeSurfaceFaces: Face[][],
  face1: Face,
  face2: Face,
) {
  // There cannot be two subgroups which will later be connected, as the "graph" of like surface faces is fully
  // connected
  for (let i = 0; i < likeSurfaceFaces.length; i++) {
    const faceGroup = likeSurfaceFaces[i]
    let foundFace1 = false,
      foundFace2 = false
    for (let j = 0; j < faceGroup.length; j++) {
      const face = faceGroup[j]
      if (face == face1) {
        foundFace1 = true
      }
      if (face == face2) {
        foundFace2 = true
      }
    }
    if (foundFace1 != foundFace2) {
      faceGroup.push(foundFace1 ? face2 : face1)
      return
    } else if (foundFace1) {
      // found both
      return
    }
  }
  // nothing found, add a new group
  likeSurfaceFaces.push([face1, face2])
}

export function assembleFaceFromLooseEdges(
  edges: Edge[],
  surface: Surface,
  originalFace: Face,
): Face {
  const visited = new Set()

  function nextStart() {
    return edges.find((edge) => !visited.has(edge))
  }

  const loops = []
  let startEdge,
    currentEdge: Edge = undefined!
  while ((startEdge = nextStart())) {
    currentEdge = startEdge
    const loop = []
    let total = 0
    do {
      visited.add(currentEdge)
      loop.push(currentEdge)
      const possibleEdges = edges.filter((edge) => currentEdge.b.like(edge.a))
      const normalAtCurrentB = surface.normalP(currentEdge.b)
      const nextEdgeIndex = indexWithMax(possibleEdges, (edge) =>
        currentEdge.bDir.angleRelativeNormal(edge.aDir, normalAtCurrentB),
      )
      currentEdge = possibleEdges[nextEdgeIndex]
    } while (startEdge != currentEdge && total++ < 200)
    assert(total != 201)
    loops.push(loop)
  }

  const assembledFaces = BRep.assembleFacesFromLoops(
    loops,
    surface,
    originalFace,
  )
  assertf(() => 1 == assembledFaces.length)
  return assembledFaces[0]
}

/**
 * ## Markdown header
 * ![foo](screenshots/Capture.PNG)
 * {@link ../screenshots/Capture.PNG}
 * find the next edge with the MAXIMUM angle
 */
export function calcNextEdgeIndex(
  currentEdge: Edge,
  possibleEdges: Edge[],
  faceNormalAtCurrentB: V3,
): int {
  let maxValue = -20,
    advanced = false,
    result = Number.MAX_SAFE_INTEGER
  const normVector = currentEdge.bDir.cross(faceNormalAtCurrentB)
  const eps = 1e-4
  const dir = sign(currentEdge.deltaT())
  const ecd = currentEdge.curve.diff(currentEdge.bT, -dir * eps).dot(normVector)
  for (let i = possibleEdges.length; i--; ) {
    const edge = possibleEdges[i]
    const angle1 = currentEdge.bDir
      .negated()
      .angleRelativeNormal(edge.aDir, faceNormalAtCurrentB)
    const angle = ((angle1 + TAU + NLA_PRECISION) % TAU) - NLA_PRECISION
    if (eq0(angle)) {
      // do advanced analysis
      if (currentEdge.curve.isColinearTo(edge.curve)) {
        continue
      }
      const edgeDir = sign(edge.deltaT())
      const iscd = edge.curve.diff(edge.aT, edgeDir * eps).dot(normVector)
      const diff = iscd - ecd
      // if diff > 0, the angle is actually ~= 0
      if (diff < 0 && (!advanced || diff > maxValue)) {
        advanced = true
        maxValue = diff
        result = i
      }
    } else if (!advanced) {
      if (gt(angle, maxValue)) {
        maxValue = angle
        result = i
      }
    }
  }
  return result == Number.MAX_SAFE_INTEGER ? 0 : result
}
export type FaceInfo = {
  face: Face
  edge: Edge
  normalAtCanonA: V3
  inside: V3
  reversed: boolean
  angle: number
}
export class BRep extends Transformable {
  static EMPTY = new BRep([], false, "BRep.EMPTY", new Map()).buildAdjacencies()
  static R3 = new BRep([], true, "BRep.R3", new Map()).buildAdjacencies()
  faces: Face[]
  infiniteVolume: boolean
  generator: string | undefined
  vertexNames: Map<V3, string> | undefined
  edgeFaces: JavaMap<Edge, FaceInfo[]> | undefined

  constructor(
    faces: Face[],
    infiniteVolume: boolean,
    generator?: string,
    vertexNames?: Map<V3, string>,
  ) {
    super()
    this.faces = faces
    assertInst(Face, ...faces)
    this.infiniteVolume = infiniteVolume
    assert(!this.infiniteVolume || true === this.infiniteVolume)
    this.generator = generator
    this.vertexNames = vertexNames
    this.edgeFaces = undefined
    //this.assertSanity()
  }

  static loop1ContainsLoop2(
    loop1: Edge[],
    ccw1: boolean,
    loop2: Edge[],
    ccw2: boolean,
    surface: Surface,
  ): boolean {
    for (const edge of loop2) {
      const loop1ContainsPoint = surface.loopContainsPoint(loop1, edge.a)
      if (PointVsFace.ON_EDGE != loop1ContainsPoint)
        return PointVsFace.INSIDE == loop1ContainsPoint
    }
    for (const edge of loop2) {
      const edgePoint = edge.curve.at(edge.aT * 0.2 + edge.bT * 0.8)
      const loop1ContainsPoint = surface.loopContainsPoint(loop1, edgePoint)
      if (PointVsFace.ON_EDGE != loop1ContainsPoint)
        return PointVsFace.INSIDE == loop1ContainsPoint
    }
    if (ccw1 != ccw2) {
      return ccw2
    }
    throw new Error(loop1.sce + loop2.sce)
  }

  static assembleFacesFromLoops(
    loops: Edge[][],
    surface: Surface,
    originalFace: Face,
    infoFactory?: FaceInfoFactory<any>,
  ): Face[] {
    type LoopInfo = { loop: Edge[]; ccw: boolean; subloops: LoopInfo[] }

    function placeRecursively(newLoopInfo: LoopInfo, loopInfos: LoopInfo[]) {
      if (loopInfos.length == 0) {
        loopInfos.push(newLoopInfo)
      } else {
        const subLoopInfo = loopInfos.find((loopInfo) =>
          BRep.loop1ContainsLoop2(
            loopInfo.loop,
            loopInfo.ccw,
            newLoopInfo.loop,
            newLoopInfo.ccw,
            surface,
          ),
        )
        if (subLoopInfo) {
          placeRecursively(newLoopInfo, subLoopInfo.subloops)
        } else {
          // newLoopInfo isnt contained by any other subLoopInfo
          for (let i = loopInfos.length; --i >= 0; ) {
            const subLoopInfo = loopInfos[i]
            //console.log("cheving subLoopInfo", surface.loopContainsPoint(newLoopInfo.edges,
            // subLoopInfo.edges[0].a))
            if (
              BRep.loop1ContainsLoop2(
                newLoopInfo.loop,
                newLoopInfo.ccw,
                subLoopInfo.loop,
                subLoopInfo.ccw,
                surface,
              )
            ) {
              newLoopInfo.subloops.push(subLoopInfo)
              loopInfos.splice(i, 1) // remove it
            }
          }
          loopInfos.push(newLoopInfo)
        }
      }
    }

    function newFacesRecursive(loopInfo: LoopInfo): void {
      // CW loops can be top level, if they are holes in the original face not contained in the new face
      if (loopInfo.ccw) {
        if (loopInfo.subloops.every((sl) => !sl.ccw)) {
          const holes = loopInfo.subloops.map((sl) => sl.loop)
          const info =
            infoFactory &&
            infoFactory.newSubFace(originalFace, surface, loopInfo.loop, holes)
          const newFace = new originalFace.constructor(
            surface,
            loopInfo.loop,
            holes,
            "genface" + getGlobalId(),
            info,
          )
          newFaces.push(newFace)
          loopInfo.subloops.forEach((sl) =>
            sl.subloops.forEach((slsl) => slsl.ccw && newFacesRecursive(slsl)),
          )
        } else {
          loopInfo.subloops.forEach((sl) => sl.ccw && newFacesRecursive(sl))
        }
      }
    }

    const newFaces: Face[] = []
    const topLevelLoops: LoopInfo[] = []
    loops.forEach((loop) =>
      placeRecursively(
        {
          loop: loop,
          ccw: surface.edgeLoopCCW(loop),
          subloops: [],
        },
        topLevelLoops,
      ),
    )
    topLevelLoops.forEach((tll) => newFacesRecursive(tll))
    return newFaces
  }

  /**
   * Create a [BRep] by concatenating the faces of other BReps. Only use this if certain that the faces of the BReps do not intersect.
   * Otherwise, use [BRep.plus].
   * @param bReps
   * @param generator
   */
  static join(bReps: BRep[], generator?: string) {
    return new BRep(
      bReps.flatMap((b2) => b2.faces),
      false,
      generator,
    )
  }

  containsPoint(p: V3, forceInsideOutside: boolean = false): boolean {
    const dirs = [
      V(-0.3920414696448526, -0.12936136783391444, -0.9108068525164064),
      V(0.6520650903544943, -0.07151288645511984, -0.7547827667692488),
      V(0.9433494201061395, -0.2402757256238473, -0.22882186797013926),
      V(0.13678704228501923, -0.04480387361087783, 0.9895867410047372),
      V(0.0662057922721913, -0.5865836917435423, 0.8071780259955845),
      V(-0.7322576567870621, -0.12953393611526787, 0.6685953061989045),
      V(0.6579719127258273, -0.012300218400456116, 0.7529420075219719),
      V(-0.5576497966736425, 0.8006695748324647, 0.2189861552871446),
    ]
    dirLoop: for (const dir of dirs) {
      const testLine = new L3(p, dir)
      let inside = this.infiniteVolume,
        minT = Infinity
      for (const face of this.faces) {
        assert(!face.surface.containsCurve(testLine))
        const ists = face.surface.isTsForLine(testLine)
        for (const t of ists) {
          const p = testLine.at(t)
          const pvf = face.containsPoint2(p)
          //assert(pvf != PointVsFace.ON_EDGE)
          !forceInsideOutside && assert(!eq0(t))
          if (t > 0) {
            if (pvf == PointVsFace.ON_EDGE) {
              continue dirLoop
            }
            if (pvf == PointVsFace.INSIDE) {
              inside = !inside
              if (t < minT) {
                minT = t
              }
            }
          }
        }
      }
      return inside
    }
    return false
  }

  withMergedFaces(): BRep {
    const likeSurfaceFaces = []
    for (let i = 0; i < this.faces.length; i++) {
      let addedToGroup = false
      for (let j = 0; j < i; j++) {
        if (this.faces[i].surface.isCoplanarTo(this.faces[j].surface)) {
          const faceGroup = likeSurfaceFaces.find((faceGroup) =>
            faceGroup.includes(this.faces[j]),
          )
          if (faceGroup) {
            faceGroup.push(this.faces[i])
            addedToGroup = true
          }
        }
      }
      !addedToGroup && likeSurfaceFaces.push([this.faces[i]])
    }

    console.log("likeSurfaceFaces", likeSurfaceFaces)
    if (likeSurfaceFaces.every((group) => group.length == 1)) return this

    const newFaces = []
    let total = 0
    for (const faceGroup of likeSurfaceFaces) {
      console.log(faceGroup)
      if (faceGroup.length == 1) {
        newFaces.push(faceGroup[0])
      } else {
        const allEdges = faceGroup.flatMap((face) => face.getAllEdges())
        for (let i = allEdges.length; i-- > 0; ) {
          for (let j = 0; j < i; j++) {
            console.log("blugh", total)
            assert(i >= 0 && j >= 0 && total++ < 500, i + " " + j + " " + total)
            if (allEdges[i].isCoEdge(allEdges[j])) {
              // remove both
              allEdges.splice(i, 1)
              allEdges.splice(j, 1)
              i--
              break
            }
          }
        }
        const newFace = assembleFaceFromLooseEdges(
          allEdges,
          faceGroup[0].surface,
          faceGroup[0],
        )
        newFaces.push(newFace)
      }
    }

    return new BRep(
      newFaces,
      this.infiniteVolume,
      this.generator && this.generator + ".withMergedFaces()",
      this.vertexNames,
    )
  }

  calculateVolume(): number {
    return sum(this.faces.map((face) => face.zDirVolume().volume))
  }

  toMesh(): Mesh & {
    faceIndexes: Map<Face, { start: int; count: int }>
    TRIANGLES: int[]
    LINES: int[]
    normals: V3[]
  } {
    const mesh = new Mesh()
      .addVertexBuffer("normals", "ts_Normal")
      .addIndexBuffer("TRIANGLES")
      .addIndexBuffer("LINES") as any
    mesh.faceIndexes = new Map()
    for (const face of this.faces) {
      const triangleStart = mesh.TRIANGLES.length
      face.addToMesh(mesh)
      mesh.faceIndexes.set(face, {
        start: triangleStart,
        count: mesh.TRIANGLES.length - triangleStart,
      })
    }
    //this.buildAdjacencies()
    //for (const edge of this.edgeFaces.keys()) {
    //
    //}
    return mesh
  }

  minus(other: BRep, infoFactory?: FaceInfoFactory<any>): BRep {
    const generator =
      this.generator &&
      other.generator &&
      this.generator + ".minus(" + other.generator + ")"
    return this.intersection(
      other.flipped(),
      true,
      true,
      generator,
      infoFactory,
    )
  }

  plus(other: BRep, infoFactory?: FaceInfoFactory<any>): BRep {
    const generator =
      this.generator &&
      other.generator &&
      this.generator + ".plus(" + other.generator + ")"
    return this.flipped()
      .intersection(other.flipped(), true, true, generator, infoFactory)
      .flipped()
  }

  and(other: BRep, infoFactory?: FaceInfoFactory<any>): BRep {
    const generator =
      this.generator &&
      other.generator &&
      this.generator + ".and(" + other.generator + ")"
    return this.intersection(other, true, true, generator, infoFactory)
  }

  xor(other: BRep, infoFactory?: FaceInfoFactory<any>): BRep {
    const generator =
      this.generator &&
      other.generator &&
      this.generator + ".xor(" + other.generator + ")"
    return new BRep(
      this.minus(other, infoFactory).faces.concat(
        other.minus(this, infoFactory).faces,
      ),
      this.infiniteVolume != other.infiniteVolume,
      generator,
    )
  }

  equals(obj: any): boolean {
    return (
      this.faces.length == obj.faces.length &&
      this.faces.every((face) =>
        (obj as BRep).faces.some((face2) => face.equals(face2)),
      )
    )
  }

  like(brep: BRep): boolean {
    return (
      this.faces.length == brep.faces.length &&
      this.faces.every((face) =>
        brep.faces.some((face2) => face.likeFace(face2)),
      )
    )
  }

  //reconstituteCoplanarFaces(likeSurfacePlanes, edgeLooseSegments, faceMap, newFaces) {
  //    likeSurfacePlanes.forEach(faceGroup => {
  //        // calculate total contours
  //        let surface = faceGroup[0].surface, bag = []
  //        faceGroup.forEach(face => {
  //            Array.prototype.push.apply(bag, faceMap(face))
  //            face.getAllEdges().forEach(edge => {
  //                let edgeSubSegments
  //                if (edgeSubSegments = edgeLooseSegments.get(edge)) {
  //                    Array.prototype.push.apply(bag, edgeSubSegments)
  //                } else {
  //                    bag.push(edge)
  //                }
  //            })
  //        })
  //        let currentEdge, loops = []
  //        while (currentEdge = bag.find(edge => !edge.visited)) {
  //            let path = []
  //            do {
  //                currentEdge.visited = true
  //                path.push(currentEdge)
  //                let possibleNextEdges = bag.filter(edge => currentEdge.b.like(edge.a))
  //                // lowest angle, i.e. the right-most next edge
  //                let nextEdgeIndex = possibleNextEdges.indexWithMax((edge, index) =>
  // -currentEdge.bDir.angleRelativeNormal(edge.aDir, surface.normalP(currentEdge.b))) currentEdge =
  // possibleNextEdges[nextEdgeIndex] } while (!currentEdge.visited) let startIndex = path.find(currentEdge) if (-1
  // != startIndex) { loops.push(path.slice(startIndex)) } } }) }

  toString(): string {
    return `new BRep([\n${this.faces.join(",\n").replace(/^/gm, "\t")}], ${
      this.infiniteVolume
    })`
  }

  getConstructorParameters() {
    return [this.faces, this.infiniteVolume]
  }

  toSource(useGenerator: boolean = true): string {
    return (
      (useGenerator && this.generator) ||
      `new BRep([\n${this.faces.map(SCE).join(",\n").replace(/^/gm, "\t")}], ${
        this.infiniteVolume
      })`
    )
  }

  /**
   * Rightmost next segment doesn't work, as the correct next segment isn't obvious from the current corner
   * alone.
   * (at least, not without extensive pre-analysis on the face edges, which shouldn't be necessary, as the
   * correct new faces are defined by the new edges already.) Leftmost edge should work. Holes which touch the
   * edge of the face will be added to the face contour.
   *
   * New segments will always be part left-er than existing ones, so no special check is required.
   *
   */
  reconstituteFaces(
    oldFaces: Face[],
    edgeSubEdges: Map<Edge, Edge[]>,
    faceMap: Map<Face, Edge[]>,
    newFaces: Face[],
    infoFactory?: FaceInfoFactory<any>,
  ): void {
    const oldFaceStatuses: Map<Face, string> = new Map()
    // reconstitute faces
    const insideEdges: Edge[] = []
    for (const face of oldFaces) {
      const usableOldEdges = face
        .getAllEdges()
        .filter((edge) => !edgeSubEdges.get(edge))
      const subEdges: Edge[] = concatenated(
        mapFilter(face.getAllEdges(), (edge) => edgeSubEdges.get(edge)),
      )
      const newEdges = faceMap.get(face) || []
      if (newEdges.length || subEdges.length) {
        oldFaceStatuses.set(face, "partial")
        const loops = []
        // new edges are definitely part of a resulting loop
        // old edges (both contour and holes) can either be part of a new loop, in which case they will already
        // have been visited when starting a loop search with a new edge, OR they can be stranded, OR they can
        // remain in their old loop
        function getNextStart() {
          return (
            newEdges.find((edge) => !visitedEdges.has(edge)) ||
            subEdges.find((edge) => !visitedEdges.has(edge)) ||
            usableOldEdges.find((edge) => !visitedEdges.has(edge))
          )
        }

        const visitedEdges = new Set()

        // search for a loop:
        let currentEdge: Edge | undefined
        while ((currentEdge = getNextStart())) {
          const startEdge = currentEdge,
            edges: Edge[] = []
          let i = 0
          // wether only new edges are used (can include looseSegments)
          do {
            visitedEdges.add(currentEdge)
            edges.push(currentEdge)
            // find next edge
            const possibleOldEdges = usableOldEdges.filter((edge) =>
              currentEdge!.b.like(edge.a),
            )
            const possibleSubEdges = subEdges.filter((edge) =>
              currentEdge!.b.like(edge.a),
            )
            const possibleNewEdges = newEdges.filter((edge) =>
              currentEdge!.b.like(edge.a),
            )
            const possibleEdges = possibleOldEdges.concat(
              possibleSubEdges,
              possibleNewEdges,
            )
            if (0 == possibleEdges.length) break
            assert(0 < possibleEdges.length, () => face.sce)
            const faceNormalAtCurrentB = face.surface.normalP(currentEdge.b)
            const nextEdgeIndex = calcNextEdgeIndex(
              currentEdge,
              possibleEdges,
              faceNormalAtCurrentB,
            )
            currentEdge = possibleEdges[nextEdgeIndex]
            if (visitedEdges.has(currentEdge)) {
              break
            }
            assert(currentEdge)
            assert(currentEdge != startEdge)
          } while (++i < 400)
          if (400 == i) {
            assert(false, "too many")
          }
          // check if we found a loop
          if (edges.length > 1 && currentEdge == startEdge) {
            loops.push(edges)
          }
        }
        const faceNewFaces = BRep.assembleFacesFromLoops(
          loops,
          face.surface,
          face,
          infoFactory,
        )
        newFaces.push(...faceNewFaces)
        const faceNewFacesEdges = faceNewFaces.flatMap((face) =>
          face.getAllEdges(),
        )
        insideEdges.push(
          ...usableOldEdges.filter((edge) => faceNewFacesEdges.includes(edge)),
        )
      }
    }
    while (insideEdges.length != 0) {
      const insideEdge = insideEdges.pop()!
      const adjacentFaces = this.edgeFaces!.get(insideEdge.getCanon())!
      adjacentFaces.forEach((info) => {
        if (!oldFaceStatuses.has(info.face)) {
          oldFaceStatuses.set(info.face, "inside")
          insideEdges.push.apply(insideEdges, info.face.getAllEdges())
        }
      })
    }
    newFaces.push(
      ...oldFaces.filter((face) => oldFaceStatuses.get(face) == "inside"),
    )
  }

  static getLooseEdgeSegments(
    edgePointInfoss: JavaMap<Edge, IntersectionPointInfo[]>,
    edgeFaces: JavaMap<Edge, FaceInfo[]>,
  ): Map<Edge, Edge[]> {
    const result = new JavaMap<Edge, Edge[]>()
    // if there are no point info, the original edge will be kept, so we should return nothing
    // otherwise, something will be returned, even if it a new edge identical to the base edge
    for (const [canonEdge, pointInfos] of edgePointInfoss) {
      if (0 == pointInfos.length) continue
      const allFaces = edgeFaces.get(canonEdge)!
      pointInfos.sort((a, b) => snap0(a.edgeT - b.edgeT) || +!!undefined)
      let startP = canonEdge.a,
        startDir = canonEdge.aDir,
        startT = canonEdge.aT,
        startInfo

      function addNewEdge(
        startInfo: IntersectionPointInfo | undefined,
        endInfo: IntersectionPointInfo | undefined,
        newEdge: Edge,
      ) {
        for (let i = 0; i < allFaces.length; i++) {
          const faceInfo = allFaces[i]
          mapPush(
            result,
            !faceInfo.reversed ? canonEdge : canonEdge.flipped(),
            !faceInfo.reversed ? newEdge : newEdge.flipped(),
          )
        }
      }

      for (let i = 0; i < pointInfos.length; i++) {
        const info = pointInfos[i]
        const pDir = canonEdge.tangentAt(info.edgeT)
        if (!eq(info.edgeT, startT)) {
          const newEdge = createEdge(
            canonEdge.curve,
            startP,
            info.p,
            startT,
            info.edgeT,
            undefined,
            startDir,
            pDir,
            "looseSegment" + getGlobalId(),
          )
          addNewEdge(startInfo, info, newEdge)
        }
        startP = info.p
        startT = info.edgeT
        startInfo = info
        startDir = pDir
      }
      if (startInfo && !eq(startT, canonEdge.bT)) {
        const newEdge = createEdge(
          canonEdge.curve,
          startP,
          canonEdge.b,
          startT,
          canonEdge.bT,
          undefined,
          startDir,
          canonEdge.bDir,
          "looseSegment" + getGlobalId(),
        )
        addNewEdge(startInfo, undefined, newEdge)
      }
    }
    return result
  }

  getIntersectionEdges(brep2: BRep) {
    const faceMap = new Map(),
      thisEdgePoints = new JavaMap<Edge, IntersectionPointInfo[]>(),
      otherEdgePoints = new JavaMap<Edge, IntersectionPointInfo[]>()

    const checkedPairs = new JavaSet<Pair<any, any>>()

    this.faces.forEach((face) => {
      //console.log('face', face.toString())
      brep2.faces.forEach((face2) => {
        //console.log('face2', face2.toString())
        face.intersectFace(
          face2,
          this,
          brep2,
          faceMap,
          thisEdgePoints,
          otherEdgePoints,
          checkedPairs,
        )
      })
    })

    return concatenated(Array.from(faceMap.values()))
  }

  shellCount(): int {
    const foundFaces = new Set<Face>()
    let face: Face | undefined,
      result = 0
    while ((face = this.faces.find((face) => !foundFaces.has(face)))) {
      result++
      const stack = [face]
      while ((face = stack.pop())) {
        // @ts-ignore
        for (const edge of face.getAllEdges()) {
          // @ts-ignore
          for (const { face: face2 } of this.edgeFaces!.get(edge.getCanon())!) {
            if (face !== face2 && !foundFaces.has(face2)) {
              foundFaces.add(face2)
              stack.push(face2)
            }
          }
        }
      }
    }
    return result
  }

  getAABB(): AABB {
    return AABB.forAABBs(this.faces.map((face) => face.getAABB()))
  }

  assertSanity(): void {
    if (!NLA_DEBUG) return
    // const allFaceEdges = this.faces.flatMap(face => face.getAllEdges())
    // for (const { i, j } of combinations(allFaceEdges.length)) {
    // const a = allFaceEdges[i],
    // 	b = allFaceEdges[j]
    // assert(i == j || !a.isCoEdge(b) || a == b || a.flippedOf == b, 'coedges not linked properly', a, b)
    // assert(
    // 	i == j ||
    // 		!a.curve.isColinearTo(b.curve) ||
    // 		(a.curve.equals(b.curve) && a.isCoEdge(b)) ||
    // 		!a.overlaps(b),
    // 	'colinear edges overlap',
    // 	a,
    // 	b,
    // )
    // }

    this.buildAdjacencies()
    for (const [canonEdge, edgeFaceInfos] of this.edgeFaces!) {
      // TODO handle curved faces
      assert(edgeFaceInfos.length % 2 == 0, () => canonEdge + edgeFaceInfos.sce)
    }
  }

  //intersection3(other: BRep, buildThis: boolean, buildOther: boolean, name?: string): BRep {
  //    this.assertSanity()
  //    other.assertSanity()
  //    this.buildAdjacencies()
  //    other.buildAdjacencies()
  //
  //    // edge / edge
  //    for (const [edge1, edge1Faces] of this.edgeFaces) {
  //        for (const [edge2, edge2Faces] of other.edgeFaces) {
  //            const curve1 = edge1.curve, curve2 = edge2.curve
  //            if (curve1.isColinearTo(curve2)) {
  //                if (edge1.overlaps(edge2)) {
  //                    // faces have a common edge
  //                    const aT = curve1.pointT(edge2.a), bT = curve1.pointT(edge2.a)
  //                    const minT = min(aT, bT), maxT = max(aT, bT)
  //                    const commonEdge = createEdge(curve1, min(edge1.minT, minT), min(edge1.maxT, maxT), )
  //                }
  //            } else if (x = curve1.isInfosWithCurve(edge2.curve)) {
  //                // edges intersect in a point
  //            }
  //        }
  //    }
  //
  //    // point / edge
  //    function pointEdge(b1, b2, has, add) {
  //        for (const v1 of this.vertFaces.keys()) {
  //            for (const edge2 of other.edgeFaces.keys()) {
  //                if (edge2.curve.containsPoint(v1)) {
  //                    const edge2T = edge2.curve.pointT(v1)
  //                    if (eq(edge2.aT, edge2T) || eq(edge2.bT, edge2T)) {
  //                        add(v1, eq(edge2.aT, edge2T) ? edge2.a : edge2.b)
  //                    }
  //                }
  //            }
  //        }
  //    }
  //    const pairs: CustomSet<[Equalable, Equalable]> = new CustomSet<[Equalable, Equalable]>()
  //    pointEdge(this, other, (a, b) => pairs.has([a, b]), (a, b) => pairs.add([a, b]))
  //    pointEdge(other, this, (b, a) => pairs.has([a, b]), (b, a) => pairs.add([a, b]))
  //
  //
  //    // point / point
  //    for (const v1 of this.vertFaces.keys()) {
  //        for (const v2 of other.vertFaces.keys()) {
  //            if (v1.like(v2)) {
  //
  //            }
  //        }
  //    }
  //
  //    for (const face1 of this.faces) {
  //        for (const face2 of other.faces) {
  //            face1.intersectFace(face2)
  //        }
  //    }
  //
  //}

  buildAdjacencies(): this & {
    edgeFaces: JavaMap<
      Edge,
      {
        face: Face
        edge: Edge
        normalAtCanonA: V3
        inside: V3
        reversed: boolean
        angle: number
      }[]
    >
  } {
    if (this.edgeFaces) return this as any

    this.edgeFaces = new JavaMap()
    for (const face of this.faces) {
      for (const edge of face.getAllEdges()) {
        const canon = edge.getCanon()
        const normalAtCanonA = face.surface.normalP(canon.a)
        const inside = normalAtCanonA.cross(
          canon == edge ? edge.aDir : edge.bDir,
        )
        mapPush(this.edgeFaces!, canon, {
          face: face,
          edge: edge,
          normalAtCanonA: normalAtCanonA,
          reversed: canon != edge,
          inside: inside,
          angle: 0,
        })
      }
    }

    for (const [canonEdge, edgeFaceInfos] of this.edgeFaces!) {
      // TODO handle curved faces
      //assert(edgeFaceInfos.length % 2 == 0, () => canonEdge + edgeFaceInfos.sce)
      const faceInfo0 = edgeFaceInfos.find((faceInfo) => faceInfo.reversed)
      if (!faceInfo0) {
        console.warn("invalid brep")
        continue
      }
      edgeFaceInfos.forEach((faceInfo) => {
        if (faceInfo != faceInfo0) {
          faceInfo.angle = faceInfo0.inside.angleRelativeNormal(
            faceInfo.inside,
            canonEdge.aDir.unit(),
          )
          if (faceInfo.angle < 0) faceInfo.angle += 2 * Math.PI
        }
      })
      edgeFaceInfos.sort((a, b) => snap(a.angle - b.angle, 0)) // TODO  || assertNever()
    }

    return this as any
  }

  /**
   * Cases for volumes A and B
   *
   *          1.  Volumes do not touch.
   *          2.  face/face Face surfaces intersect each other.
   *              implies edges going through faces.
   *              e.g. box(5, 5, 5) - box(5, 5, 5).translate(1, 1, 1)
   *          3.  face/edge Edge of A lies in a face of B
   *              implies vertices of A lying in face of B
   *              e.g. box(5, 5, 5) - box(3, 3, 3).rotateZ([0, 1, 2] * PI / 2).translate(0, 1, 1)
   *          4.  edge/edge Two edges are colinear.
   *              implies vertex of A lying in edge of B
   *           5.  vertex/edge Vertex of A lies on edge of B (but no edge/edge)
   *          6.  vertex/vertex with/without edge/edge, edge/face and face/face intersections
   *          7.  vertex lies in face
   *
   *
   *
   */
  intersection(
    other: BRep,
    buildThis: boolean,
    buildOther: boolean,
    generator?: string,
    infoFactory?: FaceInfoFactory<any>,
  ): BRep {
    this.assertSanity()
    other.assertSanity()
    this.buildAdjacencies()
    other.buildAdjacencies()

    const faceMap = new Map()
    const thisEdgePoints = new JavaMap<Edge, IntersectionPointInfo[]>(),
      otherEdgePoints = new JavaMap<Edge, IntersectionPointInfo[]>()

    const checkedPairs = new CustomSet<Pair<any, any>>()

    for (const thisFace of this.faces) {
      for (const otherFace of other.faces) {
        thisFace.intersectFace(
          otherFace,
          this,
          other,
          faceMap,
          thisEdgePoints,
          otherEdgePoints,
          checkedPairs,
        )
      }
    }
    for (const edge of thisEdgePoints.keys()) {
      assert(this.edgeFaces!.get(edge))
    }
    for (const edge of otherEdgePoints.keys()) {
      assert(other.edgeFaces!.get(edge))
    }
    const newFaces: Face[] = []

    if (
      0 == faceMap.size &&
      0 == thisEdgePoints.size &&
      0 == otherEdgePoints.size
    ) {
      const thisInOther =
        other.containsPoint(this.faces[0].contour[0].a, true) !==
        other.infiniteVolume
      const otherInThis =
        !thisInOther &&
        this.containsPoint(other.faces[0].contour[0].a) !== this.infiniteVolume
      if (thisInOther || otherInThis) {
        const [inside, outside] = thisInOther ? [this, other] : [other, this]
        if (inside.infiniteVolume) {
          if (outside.infiniteVolume) {
            return outside
          } else {
            return BRep.join([inside, outside])
          }
        } else {
          if (outside.infiniteVolume) {
            return BRep.EMPTY
          } else {
            return inside
          }
        }
      } else {
        if (this.infiniteVolume) {
          if (other.infiniteVolume) {
            return BRep.join([this, other])
          } else {
            other
          }
        } else {
          if (other.infiniteVolume) {
            return this
          } else {
            return BRep.EMPTY
          }
        }
      }

      return BRep.EMPTY
    } else {
      if (buildThis) {
        const edgeLooseSegments = BRep.getLooseEdgeSegments(
          thisEdgePoints,
          this.edgeFaces!,
        )
        // @ts-ignore
        const els = this.faces.map((face) => [
          face,
          Array.from(edgeLooseSegments.entries()).flatMap(([edge, subs]) =>
            face.getAllEdges().some((e) => e.equals(edge)) ? subs : [],
          ),
        ])
        this.reconstituteFaces(
          this.faces,
          edgeLooseSegments,
          faceMap,
          newFaces,
          infoFactory,
        )
      }
      if (buildOther) {
        const edgeLooseSegments = BRep.getLooseEdgeSegments(
          otherEdgePoints,
          other.edgeFaces!,
        )
        // @ts-ignore
        const els = other.faces.map((face) => [
          face,
          Array.from(edgeLooseSegments.entries()).flatMap(([edge, subs]) =>
            face.getAllEdges().some((e) => e.equals(edge)) ? subs : [],
          ),
        ])
        other.reconstituteFaces(
          other.faces,
          edgeLooseSegments,
          faceMap,
          newFaces,
          infoFactory,
        )
      }
    }
    //buildCoplanar && this.reconstituteCoplanarFaces(likeSurfaceFaces, edgeLooseSegments, faceMap, newFaces,
    // this.infiniteVolume, other.infiniteVolume)

    const result = new BRep(
      newFaces,
      this.infiniteVolume && other.infiniteVolume,
      generator,
    )
    //result.buildAdjacencies()
    return result
  }

  transform(m4: M4, desc?: string) {
    let vertexNames: Map<V3, string> | undefined
    if (this.vertexNames) {
      vertexNames = new Map()
      this.vertexNames.forEach((name, vertex) =>
        vertexNames!.set(m4.transformPoint(vertex), name + desc),
      )
    }
    return new BRep(
      this.faces.map((f) => f.transform(m4)),
      this.infiniteVolume,
      this.generator && desc && this.generator + desc, // if desc isn't set, the generator will be invalid
      vertexNames,
    ) as this
  }
  transform4(m4: M4, desc?: string) {
    let vertexNames: Map<V3, string> | undefined
    if (this.vertexNames) {
      vertexNames = new Map()
      this.vertexNames.forEach((name, vertex) =>
        vertexNames!.set(m4.transformPoint(vertex), name + desc),
      )
    }
    return new BRep(
      this.faces.map((f) => f.transform4(m4)),
      this.infiniteVolume,
      this.generator && desc && this.generator + desc, // if desc isn't set, the generator will be invalid
      vertexNames,
    ) as this
  }

  flipped(): BRep {
    return new BRep(
      this.faces.map((f) => f.flipped()),
      !this.infiniteVolume,
      this.generator && this.generator + ".flipped()",
      this.vertexNames,
    )
  }
}

export type IntersectionPointInfo = {
  p: V3 // intersection point
  insideDir: V3
  t: number // param on intersection curve
  edge: Edge // face edge doing the intersection
  edgeT: number
  colinear: boolean // whether edge is colinear to intersection line
  used?: boolean
}

export function dotCurve(v: V3, cDir: V3, cDDT: V3): number {
  let dot = v.dot(cDir)
  if (eq0(dot)) {
    dot = v.dot(cDDT)
  }
  assert(!eq0(dot))
  return dot
}

export function dotCurve2(
  curve: Curve,
  t: number,
  normal: V3,
  sign: number,
): number {
  assert(sign == 1 || sign == -1, sign)
  const tangentDot = curve.tangentAt(t).dot(normal)
  // if tangentDot != 0 the curve simply crosses the plane
  if (!eq0(tangentDot)) {
    return sign * tangentDot
  }
  if (curve.ddt) {
    const ddtDot = curve.ddt(t).dot(normal)
    // tangentDot == 0 ==> critical point at t, if ddtDot != 0, then it is a turning point, otherwise we can't be sure
    // and must do a numeric test
    if (!eq0(ddtDot)) {
      return ddtDot
    }
  }
  const numericDot = curve
    .at(t)
    .to(curve.at(t + sign * 4 * NLA_PRECISION))
    .dot(normal)
  assert(!(curve instanceof L3))
  return numericDot
}

export const INSIDE = 0,
  OUTSIDE = 1,
  COPLANAR_SAME = 2,
  COPLANAR_OPPOSITE = 3,
  ALONG_EDGE_OR_PLANE = 4

/**
 *
 * @param brep BREP to check
 * @param edge edge to check
 * @param dirAtEdgeA the direction vector to check
 * @param faceNormal If dirAtEdgeA doesn't split a volume, but is along a face, the returned value depends on
 *     wether that faces normal1 points in the same direction as faceNormal
 * @returns INSIDE, OUTSIDE, COPLANAR_SAME or COPLANAR_OPPOSITE
 */
//function splitsVolumeEnclosingFaces(brep: BRep, edge: Edge, dirAtEdgeA: V3, faceNormal: V3): int {
//    assert(arguments.length == 4)
//    //assert(p.equals(edge.a))
//    const ab1 = edge.aDir.unit()
//    const relFaces = facesWithEdge(edge, brep.faces) as any[]
//    relFaces.forEach(faceInfo => {
//        faceInfo.normalAtEdgeA = faceInfo.face.surface.normalP(edge.a)
//        faceInfo.edgeDirAtEdgeA = !faceInfo.reversed
//            ? faceInfo.edge.aDir
//            : faceInfo.edge.bDir
//        faceInfo.outsideVector = faceInfo.edgeDirAtEdgeA.cross(faceInfo.normalAtEdgeA)
//        faceInfo.angle = (dirAtEdgeA.angleRelativeNormal(faceInfo.outsideVector.negated(), ab1) + 2 * Math.PI +
// NLA_PRECISION / 2) % (2 * Math.PI) }) assert(relFaces.length != 0, edge.toSource()) relFaces.sort((a, b) => a.angle
// - b.angle) // assert(relFaces.length % 2 == 0, edge.toSource()) // even number of touching faces  if
// (eq0(relFaces[0].angle)) { //assert(false) todo const coplanarSame = relFaces[0].normalAtEdgeA.dot(faceNormal) > 0;
// return coplanarSame ? COPLANAR_SAME : COPLANAR_OPPOSITE } else { return !relFaces[0].reversed ? INSIDE : OUTSIDE } }
export function splitsVolumeEnclosingFaces(
  brep: BRep,
  canonEdge: Edge,
  dirAtEdgeA: V3,
  faceNormal: V3,
): int {
  assert(arguments.length == 4)
  assert(canonEdge == canonEdge.getCanon())
  //assert(p.equals(canonEdge.a))
  const edgeFaceInfos = brep.edgeFaces!.get(canonEdge) as any[]
  assertf(() => edgeFaceInfos.length % 2 == 0)
  assertf(() => brep.edgeFaces)
  const faceInfo0 = edgeFaceInfos[0]
  const aDir1 = canonEdge.aDir.unit()
  const angleToCanon =
    ((faceInfo0.inside.angleRelativeNormal(dirAtEdgeA, aDir1) +
      2 * Math.PI +
      NLA_PRECISION) %
      (2 * Math.PI)) -
    NLA_PRECISION
  const nearestFaceInfoIndex = edgeFaceInfos.findIndex((faceInfo) =>
    lt(angleToCanon, faceInfo.angle),
  )
  const nearestFaceInfo =
    edgeFaceInfos[
      nearestFaceInfoIndex == -1
        ? edgeFaceInfos.length - 1
        : nearestFaceInfoIndex - 1
    ]
  if (eq(nearestFaceInfo.angle, angleToCanon)) {
    //assert(false) todo
    const coplanarSame = nearestFaceInfo.normalAtCanonA.dot(faceNormal) > 0
    return coplanarSame ? COPLANAR_SAME : COPLANAR_OPPOSITE
  } else {
    return nearestFaceInfo.reversed ? INSIDE : OUTSIDE
  }
}

export function splitsVolumeEnclosingFacesP(
  brep: BRep,
  canonEdge: Edge,
  p: V3,
  pInside: V3,
  pFaceNormal: V3,
): int {
  assert(arguments.length == 5)
  assert(canonEdge == canonEdge.getCanon())
  //assert(p.equals(canonEdge.a))
  assertf(() => brep.edgeFaces)
  const edgeFaceInfos = brep.edgeFaces!.get(canonEdge)!
  assertf(() => edgeFaceInfos.length % 2 == 0)
  const pDir1 = canonEdge.tangentAt(canonEdge.curve.pointT(p)).unit()
  const faceInfoAngleFromPInsideNeg = (faceInfo: FaceInfo) => {
    const faceInfoPDir =
      faceInfo.edge.getCanon() == faceInfo.edge ? pDir1 : pDir1.negated()
    const faceInfoInsideAtP = faceInfo.face.surface
      .normalP(p)
      .cross(faceInfoPDir)
    const faceInfoAngleAtP = pInside.angleRelativeNormal(
      faceInfoInsideAtP,
      pDir1,
    )
    return -(((faceInfoAngleAtP + TAU + NLA_PRECISION) % TAU) - NLA_PRECISION)
  }
  const nearestFaceInfo = withMax(edgeFaceInfos, faceInfoAngleFromPInsideNeg)!
  if (eq0(faceInfoAngleFromPInsideNeg(nearestFaceInfo))) {
    //assert(false) todo
    const coplanarSame =
      nearestFaceInfo.face.surface.normalP(p).dot(pFaceNormal) > 0
    return coplanarSame ? COPLANAR_SAME : COPLANAR_OPPOSITE
  } else {
    return nearestFaceInfo.reversed ? OUTSIDE : INSIDE
  }
}

export function splitsVolumeEnclosingFacesP2(
  brep: BRep,
  canonEdge: Edge,
  p: V3,
  testCurve: Curve,
  curveT: number,
  dir: -1 | 1,
  faceNormal: V3,
): int {
  assert(canonEdge == canonEdge.getCanon())
  //assert(p.equals(canonEdge.a))
  assertf(() => brep.edgeFaces)
  const edgeFaceInfos = brep.edgeFaces!.get(canonEdge) as any[]
  assertf(() => edgeFaceInfos.length % 2 == 0)
  const pDir1 = canonEdge.tangentAt(canonEdge.curve.pointT(p)).unit()
  let pInside = testCurve.tangentAt(curveT).times(dir)
  if (pInside.isParallelTo(pDir1)) {
    pInside = testCurve
      .diff(curveT, (1e-4 * dir) / testCurve.tangentAt(curveT).length())
      .rejectedFrom(pDir1)
    pInside = pInside.div(pInside.length())
  }
  let minValue = 20,
    advanced = false,
    result = OUTSIDE
  for (const faceInfo of edgeFaceInfos) {
    const faceInfoPDir =
      faceInfo.edge.getCanon() == faceInfo.edge ? pDir1 : pDir1.negated()
    const faceInfoInsideAtP = faceInfo.face.surface
      .normalP(p)
      .cross(faceInfoPDir)
    const faceInfoAngleAtP = pInside.angleRelativeNormal(
      faceInfoInsideAtP,
      pDir1,
    )
    const angle =
      ((faceInfoAngleAtP + TAU + NLA_PRECISION) % TAU) - NLA_PRECISION
    if (eq0(angle)) {
      // do advanced analysis
      const normVector = faceInfo.face.surface.normalP(p)
      if (faceInfo.face.surface.containsCurve(testCurve)) {
        const coplanarSame = normVector.dot(faceNormal) > 0
        return coplanarSame ? COPLANAR_SAME : COPLANAR_OPPOSITE
      }
      const testPlane = P3.normalOnAnchor(pDir1, p)
      const isCurve = faceInfo.face.surface.isCurvesWithPlane(testPlane)[0]
      const isCurvePT = isCurve.pointT(p)
      const dirFactor = sign(isCurve.tangentAt(isCurvePT).dot(pInside))
      const eps = 1e-4
      const iscd = isCurve
        .at(isCurvePT)
        .to(isCurve.at(isCurvePT + dir * dirFactor * eps))
        .dot(normVector)
      const ecd = testCurve
        .at(curveT)
        .to(testCurve.at(curveT + dir * eps))
        .dot(normVector)
      const diff = (iscd - ecd) * (faceInfo.reversed ? -1 : 1)
      if (diff > 0 && (!advanced || diff < minValue)) {
        advanced = true
        minValue = diff
        result = faceInfo.reversed ? OUTSIDE : INSIDE
      }
    } else if (!advanced) {
      if (angle < minValue) {
        minValue = angle
        result = faceInfo.reversed ? OUTSIDE : INSIDE
      }
    }
  }
  return result
}

export function splitsVolumeEnclosingCone(brep: BRep, p: V3, dir: V3) {
  const testPlane = P3.forAnchorAndPlaneVectors(p, dir, dir.getPerpendicular())
  const rays = []
  for (let k = 0; k < brep.faces.length; k++) {
    const planeFace = brep.faces[k] as PlaneFace
    assertf(() => planeFace instanceof PlaneFace)
    if (planeFace.getAllEdges().some((edge) => edge.a.like(p))) {
      if (testPlane.isParallelToPlane(planeFace.surface.plane)) {
        if (planeFace.pointsToInside(p, dir) != PointVsFace.OUTSIDE) {
          return ALONG_EDGE_OR_PLANE
        }
      } else {
        const isLine = L3.fromPlanes(testPlane, planeFace.surface.plane)
        const ps = planeFace.edgeISPsWithPlane(isLine, testPlane)
        let i = 0
        while (i < ps.length) {
          const a = ps[i++],
            b = ps[i++]
          const out = a.p.like(p)
          if (out || b.p.like(p)) {
            const dir2 = out ? isLine.dir1 : isLine.dir1.negated()
            const angle =
              (dir.angleRelativeNormal(dir2, testPlane.normal1) +
                2 * Math.PI +
                NLA_PRECISION / 2) %
              (2 * Math.PI)
            rays.push({ angle: angle, out: out })
          }
        }
      }
    }
  }
  rays.sort((a, b) => a.angle - b.angle)
  //console.log("testPlane"Plane.toSource(), "rays", rays.toSource())

  if (eq0(rays[0].angle)) {
    return ALONG_EDGE_OR_PLANE
  } else {
    return rays[0].out ? OUTSIDE : INSIDE
  }
}

export function splitsVolumeEnclosingCone2(
  brep: BRep,
  p: V3,
  curve: Curve,
  curveT: number,
  fb: 1 | -1,
) {
  assert(curve.containsPoint(p))
  const pFaces = brep.faces.filter((face) =>
    face.getAllEdges().some((edge) => edge.a.like(p)),
  )
  for (let k = 0; k < pFaces.length; k++) {
    const face = pFaces[k]
    if (face.surface.containsCurve(curve)) {
      //assert(false)
      if (face.pointsToInside3(p, curve, curveT, fb) != PointVsFace.OUTSIDE) {
        return ALONG_EDGE_OR_PLANE
      }
    }
  }
  const EPS = 1e-6
  return brep.containsPoint(curve.at(curveT + fb * EPS), true)
    ? INSIDE
    : OUTSIDE
}

export function fff(
  info: {
    face: Face
    edge: Edge
    normalAtCanonA: V3
    inside: V3
    reversed: boolean
    angle: number
  },
  surface: Surface,
): int {
  const canonA = info.edge.reversed ? info.edge.b : info.edge.a
  const surfaceNormalAtCanonA = surface.normalP(canonA)
  const dot = snap0(info.inside.dot(surfaceNormalAtCanonA))
  if (0 !== dot) {
    return 0 < dot ? OUTSIDE : INSIDE
  }
  if (surface.isCoplanarTo(info.face.surface)) {
    return 0 < info.normalAtCanonA.dot(surfaceNormalAtCanonA)
      ? COPLANAR_SAME
      : COPLANAR_OPPOSITE
  }
  throw new Error()
}

export function triangulateVertices(
  normal: V3,
  vertices: V3[],
  holeStarts: int[],
) {
  const absMaxDim = normal.maxAbsDim(),
    factor = sign(normal.e(absMaxDim))
  const contour = new Float64Array(vertices.length * 2)
  let i = vertices.length
  /*
	 var [coord0, coord1] = [['y', 'z'], ['z', 'x'], ['x', 'y']][maxAbsDim]
	 while (i--) {
	 contour[i * 2    ] = vertices[i][coord0] * factor
	 contour[i * 2 + 1] = vertices[i][coord1]
	 }
	 */

  while (i--) {
    // unroll disambiguation instead of accessing elements by string name ([coord0] etc)
    // as it confuses google closure
    switch (absMaxDim) {
      case 0:
        contour[i * 2] = vertices[i].y * factor
        contour[i * 2 + 1] = vertices[i].z
        break
      case 1:
        contour[i * 2] = vertices[i].z * factor
        contour[i * 2 + 1] = vertices[i].x
        break
      case 2:
        contour[i * 2] = vertices[i].x * factor
        contour[i * 2 + 1] = vertices[i].y
        break
    }
  }
  return earcut(contour, holeStarts)
}

/**
 * Solves a quadratic system of equations of the form
 *      a * x + b * y = c
 *      x² + y² = 1
 * This can be understood as the intersection of the unit circle with a line.
 *      => y = (c - a x) / b
 *      => x² + (c - a x)² / b² = 1
 *      => x² b² + c² - 2 c a x + a² x² = b²
 *      => (a² + b²) x² - 2 a c x + (c² - b²) = 0
 *
 * a * b + (b -c) * (b + c)
 */
export function intersectionUnitCircleLine(
  a: number,
  b: number,
  c: number,
): { x1: number; y1: number; x2: number; y2: number } {
  assertNumbers(a, b, c)
  // TODO: disambiguate on a < b
  const term = sqrt(a * a + b * b - c * c)
  return {
    x1: (a * c + b * term) / (a * a + b * b),
    x2: (a * c - b * term) / (a * a + b * b),
    y1: (b * c - a * term) / (a * a + b * b),
    y2: (b * c + a * term) / (a * a + b * b),
  }
}

export function intersectionUnitCircleLine2(
  a: number,
  b: number,
  c: number,
): [number, number][] {
  assertNumbers(a, b, c)
  // TODO: disambiguate on a < b
  // cf. pqFormula
  const termSqr = snap0(a * a + b * b - c * c)
  if (termSqr < 0) {
    return []
  } else if (termSqr == 0) {
    return [[(a * c) / (a * a + b * b), (b * c) / (a * a + b * b)]]
  } else {
    const term = sqrt(termSqr)
    return [
      [
        (a * c + b * term) / (a * a + b * b),
        (b * c - a * term) / (a * a + b * b),
      ],
      [
        (a * c - b * term) / (a * a + b * b),
        (b * c + a * term) / (a * a + b * b),
      ],
    ]
  }
}

export function intersectionCircleLine(
  a: number,
  b: number,
  c: number,
  r: number,
): { x1: number; x2: number; y1: number; y2: number } {
  assertNumbers(a, b, c, r)
  const term = sqrt(r * r * (a * a + b * b) - c * c)
  return {
    x1: (a * c + b * term) / (a * a + b * b),
    x2: (a * c - b * term) / (a * a + b * b),
    y1: (b * c - a * term) / (a * a + b * b),
    y2: (b * c + a * term) / (a * a + b * b),
  }
}

/**
 * Solves a quadratic system of equations of the form
 *      a * x + b * y = c
 *      x^2 - y^2 = 1
 * This can be understood as the intersection of the unit hyperbola with a line.
 *
 * @returns with x1 >= x2 and y1 <= y2
 * a * b + (b -c) * (b + c)
 */
export function intersectionUnitHyperbolaLine(
  a: number,
  b: number,
  c: number,
): { x1: number; y1: number; x2: number; y2: number } {
  assertNumbers(a, b, c)
  const aa = a * a,
    bb = b * b,
    cc = c * c
  // TODO: disambiguate on a < b
  //var xTerm = sqrt(4*cc*aa-4*(bb-aa)*(-cc-bb))
  const xTerm = 2 * sqrt(bb * cc + bb * bb - aa * bb)
  const yTerm = sqrt(4 * cc * bb - 4 * (bb - aa) * (cc - aa))
  return {
    x1: (-2 * a * c + xTerm) / 2 / (bb - aa),
    x2: (-2 * a * c - xTerm) / 2 / (bb - aa),
    y1: (2 * b * c - yTerm) / 2 / (bb - aa),
    y2: (2 * b * c + yTerm) / 2 / (bb - aa),
  }
}

export function curvePointPP(
  ps1: ParametricSurface,
  ps2: ParametricSurface,
  startPoint: V3,
  dontCheck?: boolean,
) {
  const EPS = NLA_PRECISION / 4
  //if (!dontCheck) {
  //    const p = curvePointPP(ps1, ps2, startPoint, true).p
  //    if (!ps1.containsPoint(p)) {
  //        console.log("foo, startPoint was " + startPoint.sce)
  //        ps1.containsPoint(p)
  //    }
  //}
  let Q = startPoint
  let st1 = ps1.pointFoot(Q)
  let st2 = ps2.pointFoot(Q)
  let a, b, aNormal, bNormal, abNormalsCross
  //console.log("curvePointPP, startPoint was " + startPoint.sce)
  //console.log(Q.sce+ ',')
  let i = 16
  do {
    a = ps1.pUV(st1.x, st1.y)
    b = ps2.pUV(st2.x, st2.y)
    if (eq0(a.distanceTo(b), EPS)) break
    // drPs.push({p:a,text:'a'+j+' '+i})
    // drPs.push({p:b,text:'b'+j+' '+i})
    aNormal = ps1.normalUV(st1.x, st1.y)
    bNormal = ps2.normalUV(st2.x, st2.y)
    // next Q is the intersection of the planes
    // (x - a) * aNormal,
    // (x - b) * bNormal and
    // (x - Q) * (aNormal X bNormal)
    abNormalsCross = aNormal.cross(bNormal)
    // drVs.push({anchor: Q, dir: aNormal})
    // drVs.push({anchor: Q, dir: bNormal})
    Q = V3.add(
      bNormal.cross(abNormalsCross).times(a.dot(aNormal)),
      abNormalsCross.cross(aNormal).times(b.dot(bNormal)),
      abNormalsCross.times(abNormalsCross.dot(Q)),
    ).div(abNormalsCross.squared())
    //console.log(Q.sce+ ',')
    // feet of Q on ps1 and ps2 (closest points)
    st1 = ps1.pointFoot(Q, st1.x, st1.y)
    st2 = ps2.pointFoot(Q, st2.x, st2.y)
  } while (--i)
  //assert(ps1.containsPoint(Q), Q, ps1)
  //assert(ps2.containsPoint(Q))
  if (!eq0(a.distanceTo(b), EPS)) {
    return undefined
  }
  return { p: Q, st1: st1, st2: st2 }
}

/**
 * Follow the intersection curve of two parametric surfaces starting from a given point.
 * @param {ParametricSurface} ps1
 * @param {ParametricSurface} ps2
 * @param {number} s1Step
 * @param {number} t1Step
 * @param {number} s2Step
 * @param {number} t2Step
 * @param {number} curveStepSize
 * @return {Curve[]}
 */
export function followAlgorithmPP(
  ps1: ParametricSurface,
  ps2: ParametricSurface,
  startPoint: V3,
  curveStepSize: number,
  bounds1: (u: number, v: number) => boolean = uvInAABB2.bind(undefined, ps1),
  bounds2: (u: number, v: number) => boolean = uvInAABB2.bind(undefined, ps2),
): { points: V3[]; tangents: V3[]; st1s: V3[]; st2s: V3[] } {
  const points: V3[] = []
  const tangents: V3[] = []
  const st1s: V3[] = []
  const st2s: V3[] = []
  let Q = startPoint
  let st1 = ps1.uvP(Q)
  let st2 = ps2.uvP(Q)
  assert(ps1.pUV(st1.x, st1.y).like(Q))
  assert(st1.like(ps1.pointFoot(Q, st1.x, st1.y)))
  assert(ps2.pUV(st2.x, st2.y).like(Q))
  assert(st2.like(ps2.pointFoot(Q, st2.x, st2.y)))
  for (let i = 0; i < 1000; i++) {
    ;({ p: Q, st1, st2 } = curvePointPP(ps1, ps2, Q)!)
    assert(ps1.containsPoint(Q), Q, ps1)
    assert(ps2.containsPoint(Q))
    const aNormal = ps1.normalUV(st1.x, st1.y)
    const bNormal = ps2.normalUV(st2.x, st2.y)
    const tangent = aNormal.cross(bNormal).toLength(curveStepSize)
    tangents.push(tangent)
    points.push(Q)
    st1s.push(st1)
    st2s.push(st2)
    if (i > 4) {
      if (!bounds1(st1.x, st1.y) || !bounds2(st2.x, st2.y)) {
        break
      }
    }
    Q = Q.plus(tangent)
  }
  return { points, tangents, st1s, st2s }
}

/**
 * Iteratively calculate points on an implicit 2D curve.
 * @param ic The curve in question.
 * @param startP The point at which to start.
 * @param stepLength The step the algorithm takes. Will be the approximate distance between points.
 * @param bounds Bounds function.
 * @param endP End point. If undefined, algorithm will continue until out of bounds or back at start point.
 * @param startTangent TODO Ignore this.
 * @returns Calculated points and tangents. points[0] and tangents[0] will be startP and startTangent.
 */
export function followAlgorithm2d(
  ic: MathFunctionR2R,
  startP: V3,
  stepLength: number = 0.5,
  bounds: AABB2,
  validUV: R2<boolean>,
  endP?: V3,
  startTangent?: V3,
): {
  points: V3[]
  tangents: V3[]
} {
  assertNumbers(stepLength, ic(0, 0))
  assertVectors(startP)
  if (!startTangent) {
    startTangent = new V3(
      -ic.y(startP.x, startP.y),
      ic.x(startP.x, startP.y),
      0,
    ).toLength(stepLength)
  }
  assertVectors(startTangent)
  const points: V3[] = []
  const tangents: V3[] = []
  assert(
    eq0(ic(startP.x, startP.y), 0.01),
    "isZero(implicitCurve(startPoint.x, startPoint.y))",
    ic(startP.x, startP.y),
  )
  let i = 0,
    p = startP,
    tangent = startTangent,
    fullLoop = false
  do {
    points.push(p)
    tangents.push(tangent)
    const searchStart = p.plus(tangent)
    assert(searchStart)
    const newP = curvePointMF(ic, searchStart)
    const dfpdx = ic.x(newP.x, newP.y),
      dfpdy = ic.y(newP.x, newP.y)
    const newTangent = new V3(-dfpdy, dfpdx, 0).toLength(stepLength)
    //const reversedDir = p.minus(prevp).dot(tangent) < 0
    assert(!p.equals(newP))
    // check if we passed a singularity
    if (tangent.dot(newTangent) < 0) {
      const singularity = newtonIterate2d(ic.x, ic.y, p.x, p.y)!
      if (
        eq0(ic(singularity.x, singularity.y)) &&
        singularity.distanceTo(p) < abs(stepLength)
      ) {
        // end on this point
        points.push(singularity)
        tangents.push(p.to(singularity))
        break
      } else {
        throw new Error()
      }
    }
    // check for endP
    if (endP && p.equals(endP)) {
      break
    }
    // check if loop
    if (fullLoop) {
      if (p.distanceTo(startP) > abs(stepLength)) {
        points.pop()
        tangents.pop()
        assert(getLast(points).distanceTo(startP) <= abs(stepLength))
        break
      }
    } else {
      if (i > 4 && p.distanceTo(startP) <= abs(stepLength)) {
        fullLoop = true
      }
    }
    // check if out of bounds
    if (i > 1 && !uvInAABB2(bounds, p.x, p.y)) {
      const endP = figureOutBorderPoint(bounds, p, ic)
      points.pop()
      tangents.pop()
      if (getLast(points).distanceTo(endP) < abs(stepLength) / 2) {
        points.pop()
        tangents.pop()
      }
      const endTangent = new V3(
        -ic.y(endP.x, endP.y),
        ic.x(endP.x, endP.y),
        0,
      ).toLength(stepLength)
      points.push(endP)
      tangents.push(endTangent)
      break
    }
    if (i > 4 && !validUV(p.x, p.y)) {
      break
    }
    assert(
      eq0(ic(newP.x, newP.y), NLA_PRECISION * 2),
      p,
      newP,
      searchStart,
      ic(newP.x, newP.y),
    )
    tangent = newTangent
    p = newP
  } while (++i < 1000)
  assert(i < 1000)

  //assert(points.length > 6)
  return { points, tangents }
}

/**
 * Given a point p just outside the bounds, figure out the nearby intersection of the bounds with the ic.
 * @param bounds
 * @param p
 * @param ic
 */
function figureOutBorderPoint(bounds: AABB2, p: V3, ic: MathFunctionR2R): V3 {
  if (p.x < bounds.uMin || bounds.uMax < p.x) {
    const u = bounds.uMax < p.x ? bounds.uMax : bounds.uMin
    const v = newtonIterateWithDerivative(
      (t) => ic(u, t),
      p.y,
      4,
      (t) => ic.y(u, t),
    )
    if (uvInAABB2(bounds, u, v)) {
      return new V3(u, v, 0)
    }
  }
  if (p.y < bounds.vMin || bounds.vMax < p.y) {
    const v = bounds.vMax < p.y ? bounds.vMax : bounds.vMin
    const u = newtonIterateWithDerivative(
      (s) => ic(s, v),
      p.x,
      4,
      (s) => ic.x(s, v),
    )
    assert(uvInAABB2(bounds, u, v))
    return new V3(u, v, 0)
  }
  throw new Error(p + " " + bounds)
}

export function followAlgorithm2dAdjustable(
  ic: MathFunctionR2R,
  start: V3,
  stepLength: number = 0.5,
  bounds: (u: number, v: number) => boolean,
  endp: V3 = start,
): { points: V3[]; tangents: V3[] } {
  assertNumbers(stepLength, ic(0, 0))
  assertVectors(start)
  //assert (!startDir || startDir instanceof V3)
  const points = []
  const tangents = []
  assert(
    eq0(ic(start.x, start.y), 0.01),
    "isZero(implicitCurve(startPoint.x, startPoint.y))",
  )
  let p = start,
    prevp = p
  let i = 0
  do {
    const dfpdx = ic.x(p.x, p.y),
      dfpdy = ic.y(p.x, p.y)
    const dfpdxx = ic.xx!(p.x, p.y),
      dfpdyy = ic.yy!(p.x, p.y),
      dfpdxy = ic.xy!(p.x, p.y)
    const c2factor = abs(
      (dfpdy ** 2 * dfpdxx - 2 * dfpdx * dfpdy * dfpdxy + dfpdx ** 2 * dfpdyy) /
        (dfpdx ** 2 + dfpdy ** 2) ** 2,
    )
    const c2 = new V3(dfpdx, dfpdy, 0).times(c2factor)
    const s = 1 / 16 / c2.length()
    const tangent = new V3(-dfpdy, dfpdx, 0).unit()
    const newPStart = p.plus(tangent.times(s).plus(c2.times(s ** 2 / 2)))
    points.push(p)
    tangents.push(tangent)
    prevp = p
    const newP = curvePointMF(ic, newPStart)
    if (newP.equals(p)) {
      assertNever()
    }
    console.log(p.to(newP).length())
    p = newP

    assert(eq0(ic(p.x, p.y)))
  } while (
    i++ < 1000 &&
    (i < 4 || prevp.distanceTo(endp) > stepLength) &&
    bounds(p.x, p.y)
  )
  assert(i != 1000)
  //assert(bounds(p.x, p.y))
  const end = i < 4 || prevp.distanceTo(endp) > stepLength ? p : endp
  const endTangent = new V3(
    -ic.y(end.x, end.y),
    ic.x(end.x, end.y),
    0,
  ).toLength(stepLength)
  points.push(end)
  tangents.push(endTangent)

  //assert(points.length > 6)
  // TODO gleichmäßige Verteilung der Punkte
  return { points, tangents }
}

// both curves must be in the same s-t coordinates for this to make sense
export function intersectionICurveICurve(
  iCurve1: MathFunctionR2R,
  startParams1: V3,
  endParams1: V3,
  startDir: V3,
  stepLength: number,
  iCurve2: MathFunctionR2R,
) {
  assertNumbers(stepLength, iCurve1(0, 0), iCurve2(0, 0))
  assertVectors(startParams1, endParams1)
  assert(!startDir || startDir instanceof V3)
  const vertices = []
  assert(eq0(iCurve1(startParams1.x, startParams1.y)))
  stepLength = stepLength || 0.5
  const eps = 1e-5
  let p = startParams1,
    prevp = p // startDir ? p.minus(startDir) : p
  let i = 0
  while (i++ < 1000 && (i < 4 || p.distanceTo(endParams1) > 1.1 * stepLength)) {
    const fp = iCurve1(p.x, p.y)
    const dfpdx = (iCurve1(p.x + eps, p.y) - fp) / eps,
      dfpdy = (iCurve1(p.x, p.y + eps) - fp) / eps
    let tangent = new V3(-dfpdy, dfpdx, 0).toLength(stepLength)
    if (p.minus(prevp).dot(tangent) < 0) tangent = tangent.negated()
    prevp = p
    p = curvePointMF(iCurve1, p.plus(tangent))
    vertices.push(p)
  }
  // TODO gleichmäßige Verteilung der Punkte
  return vertices
}

export function intersectionICurveICurve2(
  iCurve1: R2_R,
  loopPoints1: V3[],
  iCurve2: R2_R,
) {
  let p = loopPoints1[0],
    val = iCurve2(p.x, p.y),
    lastVal
  const iss = []
  for (let i = 0; i < loopPoints1.length; i++) {
    lastVal = val
    p = loopPoints1[i]
    val = iCurve2(p.x, p.y)
    if (val * lastVal <= 0) {
      // TODO < ?
      iss.push(newtonIterate2d(iCurve1, iCurve2, p.x, p.y))
    }
  }
  return iss
}

// export function intersectionPCurveISurface(
// 	parametricCurve: Curve,
// 	searchStart: number,
// 	searchEnd: number,
// 	searchStep: number,
// 	implicitSurface: ImplicitSurface,
// ) {
// 	assertNumbers(searchStart, searchEnd, searchStep)
// 	const iss = []
// 	let val = implicitSurface(parametricCurve(searchStart)),
// 		lastVal
// 	for (let t = searchStart + searchStep; t <= searchEnd; t += searchStep) {
// 		lastVal = val
// 		val = implicitSurface(parametricCurve(t))
// 		if (val * lastVal <= 0) {
// 			iss.push(newtonIterate1d(t => implicitSurface(parametricCurve(t)), t))
// 		}
// 	}
// 	return iss
// }

export function cassini(
  a: number,
  c: number,
): (x: number, y: number) => number {
  return (x, y) =>
    (x * x + y * y) * (x * x + y * y) -
    2 * c * c * (x * x - y * y) -
    (a ** 4 - c ** 4)
}

/**
 * A function R² -> R with first and (optionally) second derivatives.
 */
export interface MathFunctionR2R {
  readonly x: R2_R
  readonly y: R2_R
  readonly xx?: R2_R
  readonly xy?: R2_R
  readonly yy?: R2_R

  (u: number, v: number): number
}

export namespace MathFunctionR2R {
  export function forNerdamer(
    expression: nerdamer.ExpressionParam,
    args: [string, string] = ["x", "y"],
  ): MathFunctionR2R {
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

  export function nerdamerToR2_R(
    expression: nerdamer.Expression,
    args: [string, string] = ["x", "y"],
  ) {
    return expression.buildFunction(args)
  }

  export function forFFxFy(f: R2_R, fx: R2_R, fy: R2_R): MathFunctionR2R {
    ;(f as any).x = fx
    ;(f as any).y = fy
    return f as any
  }
}
export const cas2 = cassini(0.9, 1.02)

export function arrayLerp<T>(
  lerp: (a: T, b: T, t: number) => T,
  arr: T[],
  t: number,
): T {
  if (0 === t % 1) return arr[t]
  return lerp(arr[Math.floor(t)], arr[Math.ceil(t)], t % 1)
}
