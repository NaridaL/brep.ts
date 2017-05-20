///<reference path="Face.ts"/>
let eps = 1e-5

let globalId = 0
function addLikeSurfaceFaces(likeSurfaceFaces: Face[][], face1: Face, face2: Face) {
    // There cannot be two subgroups which will later be connected, as the "graph" of like surface faces is fully
    // connected
    for (let i = 0; i < likeSurfaceFaces.length; i++) {
        let faceGroup = likeSurfaceFaces[i]
        let foundFace1 = false, foundFace2 = false
        for (let j = 0; j < faceGroup.length; j++) {
            let face = faceGroup[j]
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

function assembleFaceFromLooseEdges(edges: Edge[], surface: Surface, faceConstructor: typeof Face.constructor): Face {
    const visited = new Set()
    function nextStart() { return edges.find(edge => !visited.has(edge)) }
    const loops = []
    let startEdge, currentEdge
    while (startEdge = nextStart()) {
        currentEdge = startEdge
        const loop = []
        let total = 0
        do {
            visited.add(currentEdge)
            loop.push(currentEdge)
            const possibleEdges = edges.filter(edge => currentEdge.b.like(edge.a))
            const normalAtCurrentB = surface.normalP(currentEdge.b)
            const nextEdgeIndex = possibleEdges.indexWithMax(
                (edge, index) => currentEdge.bDir.angleRelativeNormal(edge.aDir, normalAtCurrentB))
            currentEdge = possibleEdges[nextEdgeIndex]
        } while (startEdge != currentEdge && total++ < 200)
        assert(total != 201)
        loops.push(loop)
    }


    const assembledFaces = B2.assembleFacesFromLoops(loops, surface, faceConstructor)
    assertf(() => 1 == assembledFaces.length)
    return assembledFaces[0]
}
/**
 * ## Markdown header
 * ![foo](screenshots/Capture.PNG)
 * {@link ../screenshots/Capture.PNG}
 * find the next edge with the MAXIMUM angle
 */
function calcNextEdgeIndex(currentEdge: Edge , possibleEdges: Edge[], faceNormalAtCurrentB: V3): int {
	let maxValue = -20, advanced = false, result = Number.MAX_SAFE_INTEGER
	const normVector = currentEdge.bDir.cross(faceNormalAtCurrentB)
	const eps = 1e-4
	const dir = sign(currentEdge.deltaT())
	const ecd = currentEdge.curve.diff(currentEdge.bT, - dir * eps).dot(normVector)
	for (let i = possibleEdges.length ; i--;) {
		const edge = possibleEdges[i]
		const angle1 = currentEdge.bDir.negated().angleRelativeNormal(edge.aDir, faceNormalAtCurrentB)
		const angle = (angle1 + TAU + NLA_PRECISION) % TAU - NLA_PRECISION
		if (eq0(angle)) {
			// do advanced analysis
			if (currentEdge.curve.isColinearTo(edge.curve)) {
				continue
			}
			const edgeDir = sign(edge.deltaT())
			const iscd = edge.curve.diff(edge.aT, edgeDir * eps).dot(normVector)
			const diff = (iscd - ecd)
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
	return result
}
class B2 extends Transformable {
    faces: Face[]
    infiniteVolume: boolean
    generator: string
    vertexNames: Map<V3, string>
    edgeFaces: CustomMap<Edge, {face: Face, edge: Edge, normalAtCanonA: V3, inside: V3, reversed: boolean, angle: number}[]>
    vertFaces: CustomMap<V3, Edge[]>

    constructor(faces: Face[], infiniteVolume: boolean, generator?: string, vertexNames?) {
        super()
        this.faces = faces
        assertInst.apply(undefined, [Face as any].concat(faces))
	    this.infiniteVolume = infiniteVolume
	    assert(false === this.infiniteVolume || true === this.infiniteVolume)
	    this.generator = generator
        this.vertexNames = vertexNames
        this.edgeFaces = undefined
        //this.assertSanity()
    }

    containsPoint(p: V3, forceInsideOutside: boolean = false): boolean {
	    const dirs = [V(-0.3920414696448526, -0.12936136783391444, -0.9108068525164064),V(0.6520650903544943, -0.07151288645511984, -0.7547827667692488),V(0.9433494201061395, -0.2402757256238473, -0.22882186797013926),V(0.13678704228501923, -0.04480387361087783, 0.9895867410047372),V(0.0662057922721913, -0.5865836917435423, 0.8071780259955845),V(-0.7322576567870621, -0.12953393611526787, 0.6685953061989045),V(0.6579719127258273, -0.012300218400456116, 0.7529420075219719),V(-0.5576497966736425, 0.8006695748324647, 0.2189861552871446)]
	    dirLoop: for (const dir of dirs) {
		    const testLine = new L3(p, dir)
		    let inside = this.infiniteVolume, result = false, minT = Infinity
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
    }

    withMergedFaces(): B2 {
        const likeSurfaceFaces = []
        for (let i = 0; i < this.faces.length; i++) {
            let addedToGroup = false
            for (let j = 0; j < i; j++) {
                if (this.faces[i].surface.isCoplanarTo(this.faces[j].surface)) {
                    const faceGroup = likeSurfaceFaces.find(faceGroup => faceGroup.contains(this.faces[j]))
                    if (faceGroup) {
                        faceGroup.push(this.faces[i])
                        addedToGroup = true
                    }
                }
            }
            !addedToGroup && likeSurfaceFaces.push([this.faces[i]])
        }

        console.log('likeSurfaceFaces', likeSurfaceFaces)
        if (likeSurfaceFaces.every(group => group.length == 1)) return this

        const newFaces = []
        let total = 0
        for (const faceGroup of likeSurfaceFaces) {
            console.log(faceGroup)
            if (faceGroup.length == 1) {
                newFaces.push(faceGroup[0])
            } else {
                const allEdges = faceGroup.flatMap(face => face.getAllEdges())
                for (let i = allEdges.length; i-- > 0;) {
                    for (let j = 0; j < i; j++) {
                        console.log('blugh', total)
                        assert(i >= 0 && j >= 0 && total++ < 500, i + ' '+ j+' '+ total)
                        if (allEdges[i].isCoEdge(allEdges[j])) {
                            // remove both
                            allEdges.splice(i, 1)
                            allEdges.splice(j, 1)
                            i--
                            break
                        }
                    }
                }
                const newFace = assembleFaceFromLooseEdges(allEdges, faceGroup[0].surface, faceGroup[0].constructor)
                newFaces.push(newFace)
            }
        }

        return new B2(newFaces, this.infiniteVolume, this.generator && this.generator + '.withMergedFaces()', this.vertexNames)
    }

    calculateVolume(): number {
        return this.faces.map(face => face.zDirVolume().volume).sum()
    }

    toMesh(): Mesh & {faceIndexes: Map<Face, {start: int, count: int}>} {
        const mesh = new Mesh({triangles: true, normals: true, lines: true}) as any
        mesh.faceIndexes = new Map()
        for (const face of this.faces) {
            let triangleStart = mesh.triangles.length
            face.addToMesh(mesh)
            mesh.faceIndexes.set(face, {start: triangleStart, count: mesh.triangles.length - triangleStart})
        }
        //this.buildAdjacencies()
        //for (const edge of this.edgeFaces.keys()) {
        //
        //}
	    mesh.compile()
        return mesh
    }

    minus(other: B2, infoFactory: FaceInfoFactory<any>): B2 {
        return this.intersection(other.flipped(), true, true,
	        this.generator && other.generator && this.generator + '.minus(' + other.generator + ')', infoFactory)
    }

    plus(other: B2): B2 {
        const result = this.flipped().intersection(other.flipped(), true, true).flipped()
        result.generator = this.generator && other.generator && this.generator + '.plus(' + other.generator + ')'
        return result
    }

    and(other: B2): B2 {
        return this.intersection(other, true, true,
            this.generator && other.generator && this.generator + '.and(' + other.generator + ')')
    }

    xor(other: B2): B2 {
        return new B2(
            this.minus(other).faces.concat(other.minus(this).faces),
            this.infiniteVolume != other.infiniteVolume,
            this.generator && other.generator && this.generator + '.xor(' + other.generator + ')')
    }

    equals(obj: any): boolean {
        return this.faces.length == obj.faces.length &&
            this.faces.every((face) => obj.faces.some((face2) => face.equals(face2)))
    }

    like(brep): boolean {
        return this.faces.length == brep.faces.length &&
            this.faces.every((face) => brep.faces.some((face2) => face.likeFace(face2)))
    }

    toString(): string {
        return `new B2([\n${this.faces.join(',\n').replace(/^/gm, '\t')}], ${this.infiniteVolume})`
    }

    toSource(useGenerator: boolean = true): string {
        return useGenerator && this.generator ||
            `new B2([\n${this.faces.map(SCE).join(',\n').replace(/^/gm, '\t')}], ${this.infiniteVolume})`
    }

    static loop1ContainsLoop2(loop1: Edge[], ccw1: boolean, loop2: Edge[], ccw2: boolean, surface: Surface): boolean {
        for (const edge of loop2) {
            const loop1ContainsPoint = surface.loopContainsPoint(loop1, edge.a)
            if (PointVsFace.ON_EDGE != loop1ContainsPoint) return PointVsFace.INSIDE == loop1ContainsPoint
        }
        for (const edge of loop2) {
            const edgePoint = edge.curve.at(edge.aT * 0.2 + edge.bT * 0.8)
            const loop1ContainsPoint = surface.loopContainsPoint(loop1, edgePoint)
            if (PointVsFace.ON_EDGE != loop1ContainsPoint) return PointVsFace.INSIDE == loop1ContainsPoint
        }
        if (ccw1 != ccw2) {
        	return ccw2
        }
        throw new Error(loop1.sce+ loop2.sce)
    }

    static assembleFacesFromLoops(loops: Edge[][],
                                  surface: Surface,
                                  originalFace: Face,
                                  infoFactory: FaceInfoFactory<any>): Face[] {
        type LoopInfo = {loop: Edge[], ccw: boolean, subloops: LoopInfo[]}
        function placeRecursively(newLoopInfo: LoopInfo, loopInfos: LoopInfo[]) {
            if (loopInfos.length == 0) {
                loopInfos.push(newLoopInfo)
            } else {
                const subLoopInfo = loopInfos.find(loopInfo => B2.loop1ContainsLoop2(loopInfo.loop, loopInfo.ccw, newLoopInfo.loop, newLoopInfo.ccw, surface))
                if (subLoopInfo) {
                    placeRecursively(newLoopInfo, subLoopInfo.subloops)
                } else {
                    // newLoopInfo isnt contained by any other subLoopInfo
                    for (let i = loopInfos.length; --i >= 0;) {
                        const subLoopInfo = loopInfos[i]
                        //console.log("cheving subLoopInfo", surface.loopContainsPoint(newLoopInfo.edges,
                        // subLoopInfo.edges[0].a))
                        if (B2.loop1ContainsLoop2(newLoopInfo.loop, newLoopInfo.ccw, subLoopInfo.loop, subLoopInfo.ccw, surface)) {
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
                if (loopInfo.subloops.every(sl => !sl.ccw)) {
	                const holes = loopInfo.subloops.map(sl => sl.loop)
	                const info = infoFactory && infoFactory.newSubFace(originalFace, surface, loopInfo.loop, holes)
	                const newFace = new originalFace.constructor(surface, loopInfo.loop, holes, 'genface'+globalId++, info)
                    newFaces.push(newFace)
                    loopInfo.subloops.forEach(sl => sl.subloops.forEach(slsl => slsl.ccw && newFacesRecursive(slsl)))
                } else {
                    loopInfo.subloops.forEach(sl => sl.ccw && newFacesRecursive(sl))
                }
            }
        }

        const newFaces: Face[] = []
        const topLevelLoops:LoopInfo[] = []
        loops.forEach(loop => placeRecursively({loop: loop, ccw: surface.edgeLoopCCW(loop), subloops: []}, topLevelLoops))
        topLevelLoops.forEach(tll => newFacesRecursive(tll))
        return newFaces
    }


    /**
     * Rightmost next segment doesn't work, as the correct next segment isn't obvious from the current corner
     * alone.
     * (at least, not without extensive pre-analysos on the face edges, which shouldn't be necessary, as the
     * correct new faces are defined by the new edges already.) Leftmost edge should work. Holes which touch the
     * edge of the face will be added to the face contour.
     *
     * New segments will always be part left-er than exisiting ones, so no special check is required.
     *
     */
    reconstituteFaces(oldFaces: Face[],
                      edgeSubEdges: Map<Edge, Edge[]>,
                      faceMap: Map<Face, Edge[]>,
                      newFaces: Face[],
                      infoFactory: FaceInfoFactory<any>): void {

        const oldFaceStatuses: Map<Face, string> = new Map()
        // reconstitute faces
        const insideEdges: Edge[] = []
        for (const face of oldFaces) {
            const usableOldEdges = face.getAllEdges().filter(edge => !edgeSubEdges.get(edge))
            const subEdges = face.getAllEdges().mapFilter(edge => edgeSubEdges.get(edge)).concatenated()
            const newEdges = faceMap.get(face) || []
            if (newEdges.length || subEdges.length) {
                oldFaceStatuses.set(face, 'partial')
                const loops = []
                // new edges are definitely part of a resulting loop
                // old edges (both contour and holes) can either be part of a new loop, in which case they will already
                // have been visited when starting a loop search with a new edge, OR they can be stranded, OR they can
                // remain in their old loop
                function getNextStart() {
                   return newEdges.find(edge => !visitedEdges.has(edge))
                    || subEdges.find(edge => !visitedEdges.has(edge))
                    || usableOldEdges.find(edge => !visitedEdges.has(edge))
                }
                const visitedEdges = new Set()

                // search for a loop:
                let currentEdge: Edge
                while (currentEdge = getNextStart()) {
                    const startEdge = currentEdge, edges = []
	                let i = 0
                    // wether only new edges are used (can include looseSegments)
                    do {
                        visitedEdges.add(currentEdge)
                        edges.push(currentEdge)
                        // find next edge
                        const possibleOldEdges = usableOldEdges.filter(edge => currentEdge.b.like(edge.a))
                        const possibleSubEdges = subEdges.filter(edge => currentEdge.b.like(edge.a))
                        const possibleNewEdges = newEdges.filter(edge => currentEdge.b.like(edge.a))
                        const possibleEdges = possibleOldEdges.concat(possibleSubEdges, possibleNewEdges)
                        if (0 == possibleEdges.length) break
                        assert(0 < possibleEdges.length, () => face.sce)
                        const faceNormalAtCurrentB = face.surface.normalP(currentEdge.b)
	                    const correct =  possibleEdges.indexWithMax(
		                    (edge, index) => (currentEdge.bDir.angleRelativeNormal(edge.aDir, faceNormalAtCurrentB) + NLA_PRECISION + PI) % TAU)
	                    const nextEdgeIndex = calcNextEdgeIndex(currentEdge, possibleEdges, faceNormalAtCurrentB)
                        currentEdge = possibleEdges[nextEdgeIndex]
                        if (visitedEdges.has(currentEdge)) {
                            break
                        }
                        assert(currentEdge)
                        assert(currentEdge != startEdge)
                    } while (++i < 400)
                    if (400 == i) {
                        assert(false, 'too many')
                    }
                    // check if we found a loop
                    if (edges.length > 1 && currentEdge == startEdge) {
                        loops.push(edges)
                    }
                }
	            const faceNewFaces = B2.assembleFacesFromLoops(loops, face.surface, face, infoFactory)
                newFaces.pushAll(faceNewFaces)
                const faceNewFacesEdges = faceNewFaces.flatMap(face => face.getAllEdges())
                insideEdges.pushAll(usableOldEdges.filter(edge => faceNewFacesEdges.includes(edge)))
            }
        }
	    while (insideEdges.length != 0) {
            const insideEdge = insideEdges.pop()
            const adjacentFaces = this.edgeFaces.get(insideEdge.getCanon())
            adjacentFaces.forEach(info => {
                if (!oldFaceStatuses.has(info.face)) {
                    oldFaceStatuses.set(info.face, 'inside')
                    insideEdges.push.apply(insideEdges, info.face.getAllEdges())
                }
            })
        }
        newFaces.pushAll(oldFaces.filter(face => oldFaceStatuses.get(face) == 'inside'))
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
    //                let nextEdgeIndex = possibleNextEdges.indexWithMax((edge, index) => -currentEdge.bDir.angleRelativeNormal(edge.aDir, surface.normalP(currentEdge.b)))
    //                currentEdge = possibleNextEdges[nextEdgeIndex]
    //            } while (!currentEdge.visited)
    //            let startIndex = path.find(currentEdge)
    //            if (-1 != startIndex) {
    //                loops.push(path.slice(startIndex))
    //            }
    //        }
    //    })
    //}

    getLooseEdgeSegments(
        edgePointInfoss: CustomMap<Edge, IntersectionPointInfo[]>,
        edgeFaces: CustomMap<Edge, any[]>): Map<Edge, Edge[]> {

        const result = new CustomMap<Edge, Edge[]>()
        // if there are no point info, the original edge will be kept, so we should return nothing
        // otherwise, something will be returned, even if it a new edge identical to the base edge
        for (const [canonEdge, pointInfos] of edgePointInfoss) {
            if (0 == pointInfos.length) continue
            const allFaces = edgeFaces.get(canonEdge)
            pointInfos.sort((a, b) => snap0(a.edgeT - b.edgeT) || +!!a.faces)
            let startP = canonEdge.a, startDir = canonEdge.aDir, startT = canonEdge.aT, startInfo
            function addNewEdge(startInfo, endInfo, newEdge) {
                for (let i = 0; i < allFaces.length; i++) {
                    const faceInfo = allFaces[i]
                    const startYes = !startInfo || !startInfo.faces || startInfo.faces[i]
                    const endYes = !endInfo || !endInfo.faces
                    endYes && mapPush(result,
                        !faceInfo.reversed ? canonEdge : canonEdge.flipped(),
                        !faceInfo.reversed ? newEdge : newEdge.flipped())
                }
            }
            for (let i = 0; i < pointInfos.length; i++) {
                const info = pointInfos[i]
                const pDir = canonEdge.tangentAt(info.edgeT)
                if (!eq(info.edgeT, startT) ) {
                    const newEdge = Edge.create(canonEdge.curve, startP, info.p, startT, info.edgeT, null, startDir, pDir, 'looseSegment' + globalId++)
                    addNewEdge(startInfo, info, newEdge)
                }
                startP = info.p
                startT = info.edgeT
                startInfo = info
                startDir = pDir
            }
            if (startInfo && !eq(startT, canonEdge.bT)) {
                const newEdge = Edge.create(canonEdge.curve, startP, canonEdge.b, startT, canonEdge.bT, null, startDir, canonEdge.bDir, 'looseSegment' + globalId++)
                addNewEdge(startInfo, undefined, newEdge)
            }
        }
        return result
    }

    getIntersectionEdges(brep2) {
        const faceMap = new Map(), thisEdgePoints = new CustomMap(), otherEdgePoints = new CustomMap()

        let likeSurfaceFaces = []

        this.faces.forEach(face => {
            //console.log('face', face.toString())
            brep2.faces.forEach(face2 => {
                //console.log('face2', face2.toString())
                face.intersectFace(face2, this, brep2, faceMap, thisEdgePoints, otherEdgePoints, likeSurfaceFaces)
            })
        })

        return Array.from(faceMap.values()).concatenated()

    }

    static join(b2s: B2[], generator?: string) {
        return new B2(b2s.flatMap(b2 => b2.faces), false, generator)
    }

    shellCount(): int {
        const foundFaces = new Set<Face>()
        let face, result = 0
        while (face = this.faces.find(face => !foundFaces.has(face))) {
            result++
            const stack = [face]
            while (face = stack.pop()) {
                for (const edge of face.getAllEdges()) {
                    for (const {face: face2} of this.edgeFaces.get(edge.getCanon())) {
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
        return AABB.forAABBs(this.faces.map(face => face.getAABB()))
    }

	assertSanity(): void {
        if (!NLA_DEBUG) return
		const allFaceEdges = this.faces.flatMap(face => face.getAllEdges())
		for (const {i, j} of combinations(allFaceEdges.length)) {
			const a = allFaceEdges[i], b = allFaceEdges[j]
			//assert(i == j || !a.isCoEdge(b) || a == b || a.flippedOf == b, 'coedges not linked properly', a, b)

			//assert(i == j
			//	|| !a.curve.isColinearTo(b.curve)
			//	|| (a.curve.equals(b.curve) && a.isCoEdge(b))
             //   || !a.overlaps(b), 'colinear edges overlap', a, b)
		}

		this.buildAdjacencies()
		for (const [canonEdge, edgeFaceInfos] of this.edgeFaces) {
			// TODO handle curved faces
			assert(edgeFaceInfos.length % 2 == 0, () => canonEdge + edgeFaceInfos.sce)
		}

	}


    buildAdjacencies(): void {
        if (this.edgeFaces) return

        this.edgeFaces = new CustomMap() as any
        for (const face of this.faces) {
            for (const edge of face.getAllEdges()) {
                const canon = edge.getCanon()
                const normalAtCanonA = face.surface.normalP(canon.a)
                const inside = normalAtCanonA.cross(canon == edge ? edge.aDir : edge.bDir)
                mapPush(this.edgeFaces, canon,
                    {face: face, edge: edge, normalAtCanonA: normalAtCanonA, reversed: canon != edge, inside: inside, angle: 0})
            }
        }

        for (const [canonEdge, edgeFaceInfos] of this.edgeFaces) {
	        // TODO handle curved faces
	        //assert(edgeFaceInfos.length % 2 == 0, () => canonEdge + edgeFaceInfos.sce)
            const faceInfo0 = edgeFaceInfos.find(faceInfo => faceInfo.reversed)
	        if (!faceInfo0) {
		        console.warn('invalid brep')
		        continue
	        }
            edgeFaceInfos.forEach(faceInfo => {
	            if (faceInfo != faceInfo0) {
		            faceInfo.angle = faceInfo0.inside.angleRelativeNormal(faceInfo.inside, canonEdge.aDir.unit())
		            if (faceInfo.angle < 0) faceInfo.angle += 2 * Math.PI
	            }
            })
            edgeFaceInfos.sort((a, b) => snap(a.angle - b.angle, 0)) // TODO  || assertNever()
        }
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
    intersection(other: B2, buildThis: boolean, buildOther: boolean, generator?: string, infoFactory: FaceInfoFactory<any>): B2 {
    	const link = makeLink({a: this, b: other})
        this.assertSanity()
        other.assertSanity()
        this.buildAdjacencies()
        other.buildAdjacencies()

        const faceMap = new Map()
        const thisEdgePoints = new CustomMap<Edge, IntersectionPointInfo[]>(),
	         otherEdgePoints = new CustomMap<Edge, IntersectionPointInfo[]>()

        const likeSurfaceFaces = new CustomSet()

        for (const thisFace of this.faces) {
            for (const otherFace of other.faces) {
                thisFace.intersectFace(otherFace, this, other, faceMap, thisEdgePoints, otherEdgePoints, likeSurfaceFaces)
            }
        }
        for (const edge of thisEdgePoints.keys()) {
            assert(this.edgeFaces.get(edge))
        }
        for (const edge of otherEdgePoints.keys()) {
            assert(other.edgeFaces.get(edge))
        }
        const newFaces: Face[] = []

        if (0 == faceMap.size && 0 == thisEdgePoints.size && 0 == otherEdgePoints.size) {
            const thisInOther = other.containsPoint(this.faces[0].contour[0].a, true)
            const otherInThis = !thisInOther && this.containsPoint(other.faces[0].contour[0].a)
	        return this
        } else {
            if (buildThis) {
                const edgeLooseSegments = this.getLooseEdgeSegments(thisEdgePoints, this.edgeFaces)
                const els = this.faces.map(face => [face,
                    Array.from(edgeLooseSegments.entries()).filter(([edge, subs]) => face.getAllEdges().some(e => e.equals(edge))).concatenated()])
                this.reconstituteFaces(this.faces, edgeLooseSegments, faceMap, newFaces, infoFactory)
            }
            if (buildOther) {
                const edgeLooseSegments = this.getLooseEdgeSegments(otherEdgePoints, other.edgeFaces)
                const els = other.faces.map(face => [face,
                    Array.from(edgeLooseSegments.entries())
                        .filter(([edge, subs]) => face.getAllEdges().some(e => e.equals(edge)))
                        .flatMap(([edge, subs]) => subs)])
                other.reconstituteFaces(other.faces, edgeLooseSegments, faceMap, newFaces, infoFactory)
            }
        }
        //buildCoplanar && this.reconstituteCoplanarFaces(likeSurfaceFaces, edgeLooseSegments, faceMap, newFaces,
        // this.infiniteVolume, other.infiniteVolume)

        const result = new B2(newFaces, this.infiniteVolume && other.infiniteVolume, generator)
	    //result.buildAdjacencies()
	    return result

    }

    //intersection3(other: B2, buildThis: boolean, buildOther: boolean, name?: string): B2 {
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
    //                    const commonEdge = Edge.create(curve1, min(edge1.minT, minT), min(edge1.maxT, maxT), )
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

    transform(m4: M4, desc?: string) {

        let vertexNames: Map<V3,string>
        if (this.vertexNames) {
            vertexNames = new Map()
            this.vertexNames.forEach((name, vertex) => vertexNames.set(m4.transformPoint(vertex), name + desc))
        }
        return new B2(
            this.faces.map(f => f.transform(m4)),
            this.infiniteVolume,
            this.generator && desc && this.generator + desc, // if desc isn't set, the generator will be invalid
            vertexNames
        )
    }

    flipped(): B2 {
        return new B2(
            this.faces.map(f => f.flipped()),
            !this.infiniteVolume,
            this.generator && this.generator + '.flipped()',
            this.vertexNames)
    }


    static EMPTY = new B2([], false, 'B2.EMPTY', new Map())
    static R3 = new B2([], true, 'B2.R3', new Map())
}
namespace B2 {
	export const asldk = 0
}

type IntersectionPointInfo = {
    p: V3, // intersection point
    insideDir: V3,
    t: number, // param on intersection curve
    edge: Edge, // face edge doing the intersection
    edgeT: number,
    colinear: boolean, // whether edge is colinear to intersection line
    used?: boolean }

function planeFaceEdgeISPsWithPlane(face: PlaneFace, isLine: L3, plane2: P3): IntersectionPointInfo[] {
    assert(face.surface.plane.containsLine(isLine))
    assert(plane2.containsLine(isLine))
    const plane = face.surface.plane
    const ps = []
    const loops = [face.contour].concat(face.holes)
    loops.forEach(loop => {
        const colinearEdges = loop.map((edge) => edge.colinearToLine(isLine) && -sign(edge.aDir.dot(isLine.dir1)))
        const isLineOut = isLine.dir1.cross(plane.normal1)

        loop.forEach((edge, edgeIndex, edges) => {
            const nextEdgeIndex = (edgeIndex + 1) % edges.length, nextEdge = edges[nextEdgeIndex]
            //console.log(edge.toSource()) {p:V3(2, -2.102, 0),
            if (colinearEdges[edgeIndex]) {
                // edge colinear to intersection line
                const curveAT = isLine.pointT(edge.a), curveBT = isLine.pointT(edge.b)
                // add interval for colinear segment
                ps.push(
                    {p: edge.a, insideDir: edge.aDir,           t: curveAT, edge: edge, edgeT: edge.aT, colinear: true},
                    {p: edge.b, insideDir: edge.bDir.negated(), t: curveBT, edge: edge, edgeT: edge.bT, colinear: true})
                // open next interval if necessary
                const nextSide = colinearEdges[nextEdgeIndex] || dotCurve(isLineOut, nextEdge.aDir, nextEdge.aDDT)
                if (colinearEdges[edgeIndex] * nextSide < 0) {
                    // side changes
                    ps.push({p: nextEdge.a, insideDir: edge.bDir, t: curveBT, edge: nextEdge, edgeT: nextEdge.aT, colinear: false})
                }
            } else {
                // not necessarily a straight edge, so multiple intersections are possible
                const edgeTs = edge.edgeISTsWithPlane(plane2)
                assert(edgeTs.every(t => plane2.containsPoint(edge.curve.at(t))), edgeTs)
                for (const edgeT of edgeTs) {
                    if (edgeT == edge.bT) {
                        // endpoint lies on intersection line
                        const side = -dotCurve(isLineOut, edge.bDir, edge.bDDT)
                        const nextSide = colinearEdges[nextEdgeIndex] || dotCurve(isLineOut, nextEdge.aDir, nextEdge.aDDT)
                        if (side * nextSide < 0) {
                            // next segment is not colinear and ends on different side
                            ps.push({p: edge.b, insideDir: plane2.normal1.negated(), t: isLine.pointT(edge.b), edge: edge, edgeT: edge.bT, colinear: false})
                        }
                    } else if (edgeT != edge.aT) {
                        // edge crosses intersection line, neither starts nor ends on it
                        const p = edge.curve.at(edgeT)
                        assert(plane2.containsPoint(p), edge.toString(), p, edgeT, plane2.distanceToPoint(p))
                        assert(isLine.containsPoint(p), edge.toString(), p, edgeT, isLine.distanceToPoint(p))
                        const insideDir = plane2.normal1.negated()
                        ps.push({p: p, insideDir: insideDir, t: isLine.pointT(p), edge: edge, edgeT: edgeT, colinear: false})
                    }
                }
            }
        })
    })
    // duplicate 't's are ok, as sometimes a segment needs to stop and start again
    // should be sorted so that back facing ones are first
    ps.sort((a, b) => a.t - b.t || a.insideDir.dot(isLine.dir1))
    return ps
}
function faceEdgeISPsWithSurface(face: Face, isCurve: Curve, surface2: Surface): IntersectionPointInfo[] {
    const surface = face.surface
    const loops = face.holes.concat([face.contour])
	const ps = []
	for (const loop of loops) {
		const colinearEdges: boolean[] = loop.map(edge => edge.curve.isColinearTo(isCurve))
		//const colinearSides = loop.map((edge, edgeIndex) => -1 != colinearEdges[edgeIndex]
		//            && -sign(isCurves[colinearEdges[edgeIndex]].tangentAt(edge.aT).dot(edge.aDir)))
		for (let edgeIndex = 0; edgeIndex < loop.length; edgeIndex++) {
			const edge = loop[edgeIndex]
			const nextEdgeIndex = (edgeIndex + 1) % loop.length, nextEdge = loop[nextEdgeIndex]
			//console.log(edge.toSource()) {p:V3(2, -2.102, 0),
			if (colinearEdges[edgeIndex]) {
				if (isCurve.containsPoint(edge.a)) {
					const prevEdgeIndex = (edgeIndex - 1 + loop.length) % loop.length, prevEdge = loop[prevEdgeIndex]
					const curveAT = isCurve.pointT(edge.a)
					const colinearOutA = edge.aDir.cross(surface.normalP(edge.a))
					if (!colinearEdges[prevEdgeIndex] && dotCurve2(prevEdge.curve, prevEdge.bT, colinearOutA, -sign(prevEdge.deltaT())) > 0) {
						ps.push({p: prevEdge.b, insideDir: edge.aDir.negated(), t: curveAT, edge: prevEdge, edgeT: prevEdge.bT, colinear: false})
					}
					ps.push({p: edge.a, insideDir: edge.aDir, t: curveAT, edge: edge, edgeT: edge.aT, colinear: true})
				}
				if (isCurve.containsPoint(edge.b)) {
					const curveBT = isCurve.pointT(edge.b)
					const colinearOutB = edge.bDir.cross(surface.normalP(edge.b))
					if (!colinearEdges[nextEdgeIndex] && dotCurve2(nextEdge.curve, nextEdge.aT, colinearOutB, sign(nextEdge.deltaT())) > 0) {
						ps.push({p: edge.b, insideDir: edge.bDir, t: curveBT, edge: nextEdge, edgeT: nextEdge.aT, colinear: false})
					}
					ps.push({p: edge.b, insideDir: edge.bDir.negated(), t: curveBT, edge: edge, edgeT: edge.bT, colinear: true})
				}

			} else {
				const edgeTs = edge.edgeISTsWithSurface(surface2)
				for (const edgeT of edgeTs) {
					const p = edge.curve.at(edgeT)
					if (!isCurve.containsPoint(p)) continue
					const curveT = isCurve.pointT(p)
					assert(!isNaN(curveT))
					const insideDir = edge.tangentAt(edgeT).cross(surface.normalP(p)).negated()

					const isTangent = isCurve.tangentAt(curveT)
					const dirFactor = sign(isTangent.dot(edge.curve.tangentAt(edgeT)))
					const eps = 1e-4
					const normVector = surface2.normalP(p)
					//if(!eq0(insideDir.dot(isTangent))) {
						// Edge.edgeISTsWithSurface returns snapped values, so comparison with == is ok:
						if (edgeT == edge.bT) {
							// endpoint lies on intersection line
							if (!colinearEdges[nextEdgeIndex]) {
								let thisSide = -normVector.dot(edge.bDir)
								if (eq0(thisSide)) {
									// advanced test
									const dir = -sign(edge.deltaT())
									const iscd = isCurve.at(curveT).to(isCurve.at(curveT + dir * dirFactor * eps)).dot(normVector)
									const ecd = edge.curve.at(edgeT).to(edge.curve.at(edgeT + dir * eps)).dot(normVector)
									thisSide = sign(ecd - iscd)
								}
								let nextSide = normVector.dot(nextEdge.aDir)
								if (eq0(nextSide)) {
									// advanced test
									const dirFactor = sign(snap0(isTangent.dot(nextEdge.curve.tangentAt(nextEdge.aT))))
									assert(dirFactor !== 0)
									const dir = sign(nextEdge.deltaT())
									const iscd = isCurve.at(curveT).to(isCurve.at(curveT + dir * dirFactor * eps)).dot(normVector)
									const ecd = nextEdge.curve.at(nextEdge.aT).to(nextEdge.curve.at(nextEdge.aT + dir * eps)).dot(normVector)
									nextSide = sign(ecd - iscd)
								}
								if (nextSide * thisSide < 0) {
									assert(!eq0(insideDir.dot(isTangent)))
									// next segment is not colinear and ends on different side
									ps.push({ p: edge.b, insideDir: insideDir, t: curveT, edge: edge, edgeT: edge.bT, colinear: false})
								}
							}
						} else if (edgeT != edge.aT) {
							// edge crosses/touches an intersection curve, neither starts nor ends on it
							if(eq0(insideDir.dot(isTangent))) {
								const dirFactor = sign(isTangent.dot(edge.curve.tangentAt(edgeT)))
								const eps = 1e-4
								for (const dir of [-1, 1]) {
									if (-1 == dir * dirFactor && edgeT == edge.minT ||
										1 == dir * dirFactor && edgeT == edge.maxT ||
										-1 == dir && curveT == isCurve.tMin ||
										1 == dir && curveT == isCurve.tMax) continue
									const iscd = isCurve.at(curveT).to(isCurve.at(curveT + dir * eps)).dot(insideDir)
									const ecd = edge.curve.at(edgeT).to(edge.curve.at(edgeT + dir * dirFactor * eps)).dot(insideDir)
									if (iscd - ecd > 0) {
										ps.push({
											p,
											insideDir: isTangent.times(dir),
											t: curveT,
											edge: edge,
											edgeT: edgeT,
											colinear: false
										})
									}
								}
							} else {
								ps.push({
									p: p,
									insideDir: insideDir,
									t: curveT,
									edge: edge,
									edgeT: edgeT,
									colinear: false
								})
							}
						}
					//} else {
                    //
                    //	const dirFactor = sign(isTangent.dot(edge.curve.tangentAt(edgeT)))
                    //	const eps = 1e-4
                    //	const normVector = surface2.normalP(p)
                    //	for (const dir of [-1, 1]) {
                    //		if (-1 == dir * dirFactor && edgeT == edge.minT ||
                    //			1 == dir * dirFactor && edgeT == edge.maxT ||
                    //			-1 == dir && curveT == isCurve.tMin ||
                    //			1 == dir && curveT == isCurve.tMax) continue
                    //		const iscd = isCurve.at(curveT).to(isCurve.at(curveT + dir * eps)).dot(normVector)
                    //		const ecd = edge.curve.at(edgeT).to(edge.curve.at(edgeT + dir * dirFactor *
                    // eps)).dot(normVector) if (iscd > ecd) { ps.push({p, insideDir: isTangent.times(dir * dirFactor),
                    // t: curveT, edge: edge, edgeT: edgeT, colinear: false}) } }
						//curveVsSurface(isCurve, curveT, p, surface2)
					//}
				}
			}
		}
	}
    // duplicate 't's are ok, as sometimes a segment needs to stop and start again
    // should be sorted so that back facing ones are first
	ps.sort((a, b) => a.t - b.t || a.insideDir.dot(isCurve.tangentAt(a.t)))
    return ps
}
function planeFaceEdgeISPsWithPlaneOld(face: PlaneFace, isLine: L3, plane2: P3): IntersectionPointInfo[] {
    assert(face.surface.plane.containsLine(isLine))
    assert(plane2.containsLine(isLine))
    let plane = face.surface.plane
    let ps = []
    let loops = [face.contour].concat(face.holes)
    loops.forEach(loop => {
        let colinearEdges = loop.map((edge) => edge.colinearToLine(isLine))
        let isLineOut = isLine.dir1.cross(plane.normal1)

        loop.forEach((edge, edgeIndex, edges) => {
            let nextEdgeIndex = (edgeIndex + 1) % edges.length, nextEdge = edges[nextEdgeIndex]
            //console.log(edge.toSource()) {p:V3(2, -2.102, 0),
            if (colinearEdges[edgeIndex]) {
                // edge colinear to intersection line
                let prevEdgeIndex = (edgeIndex - 1 + edges.length) % edges.length, prevEdge = edges[prevEdgeIndex]
                // if colinear, edge.curve must be a line, so colinearOut is constant
                let colinearOut = edge.aDir.cross(plane.normal1)
                let curveAT = isLine.pointT(edge.a), curveBT = isLine.pointT(edge.b)
                // close previous interval if necessary
                if (!colinearEdges[prevEdgeIndex] && dotCurve(colinearOut, prevEdge.bDir, prevEdge.bDDT) < 0) {
                    ps.push({p: prevEdge.b, insideDir: edge.aDir.negated(), t: curveAT, edge: prevEdge, edgeT: prevEdge.bT,
                        colinear: false})
                }
                // add interval for colinear segment
                ps.push(
                    {p: edge.a, insideDir: edge.aDir, t: curveAT, edge: edge, edgeT: edge.aT,
                        colinear: true},
                    {p: edge.b, insideDir: edge.bDir.negated(), t: curveBT, edge: edge, edgeT: edge.bT,
                        colinear: true})
                // open next interval if necessary
                if (!colinearEdges[nextEdgeIndex] && dotCurve(colinearOut, nextEdge.aDir, nextEdge.aDDT) > 0) {
                    ps.push({p: nextEdge.a, insideDir: edge.bDir, t: curveBT, edge: nextEdge, edgeT: nextEdge.aT, colinear: false})
                }
            } else {
                // not necessarily a straight edge, so multiple intersections are possible
                let edgeTs = edge.edgeISTsWithPlane(plane2)
                assert(edgeTs.every(t => plane2.containsPoint(edge.curve.at(t))), edgeTs)
                for (let k = 0; k < edgeTs.length; k++) {
                    let edgeT = edgeTs[k]
                    if (edgeT == edge.bT) {
                        // endpoint lies on intersection line
                        if (!colinearEdges[nextEdgeIndex]
                            && edgesDifferentSidesOfDir(isLineOut, edge, nextEdge)) {
                            // next segment is not colinear and ends on different side
                            ps.push({p: edge.b, insideDir: plane2.normal1.negated(), t: isLine.pointT(edge.b), edge: edge, edgeT: edge.bT, colinear: false})
                        }
                    } else if (edgeT != edge.aT) {
                        // edge crosses intersection line, neither starts nor ends on it
                        let p = edge.curve.at(edgeT)
                        assert(plane2.containsPoint(p), edge.toString(), p, edgeT, plane2.distanceToPoint(p))
                        assert(isLine.containsPoint(p), edge.toString(), p, edgeT, isLine.distanceToPoint(p))
                        let insideDir = plane2.normal1.negated()
                        ps.push({p: p, insideDir: insideDir, t: isLine.pointT(p), edge: edge, edgeT: edgeT, colinear: false})
                    }
                }
            }
        })
    })
    // duplicate 't's are ok, as sometimes a segment needs to stop and start again
    // should be sorted so that back facing ones are first
    ps.sort((a, b) => a.t - b.t || a.insideDir.dot(isLine.dir1))
    return ps
}
function dotCurve(v: V3, cDir: V3, cDDT: V3): number {
	let dot = v.dot(cDir)
	if (eq0(dot)) { dot = v.dot(cDDT) }
	assert(!eq0(dot))
	return dot
}
function dotCurve2(curve: Curve, t: number, normal: V3, sign: number): number {
	assert(sign == 1 || sign == -1, sign)
	const tangentDot = curve.tangentAt(t).dot(normal)
	// if tangentDot != 0 the curve simply crosses the plane
	if (!eq0(tangentDot)) { return sign * tangentDot }
	const ddtDot = curve.ddt(t).dot(normal)
	// tangentDot == 0 ==> critical point at t, if ddtDot != 0, then it is a turning point, otherwise we can't be sure
    // and must do a numeric test
	if (!eq0(ddtDot)) { return ddtDot }
	const numericDot = curve.at(t).to(curve.at(t + sign * 4 * NLA_PRECISION)).dot(normal)
	assert(!(curve instanceof L3))
	return numericDot
}
function edgesDifferentSidesOfDir(dir: V3, e1: Edge, e2: Edge): boolean {
    let factor1 = dir.dot(e1.bDir)
    if (eq0(factor1)) { factor1 = dir.dot(e1.bDDT) }
    assert(!eq0(factor1))
    let factor2 = dir.dot(e2.aDir)
    if (eq0(factor2)) { factor2 = dir.dot(e2.aDDT) }
    assert(!eq0(factor2))
    return factor1 * factor2 > 0
}


const INSIDE = 0, OUTSIDE = 1, COPLANAR_SAME = 2, COPLANAR_OPPOSITE= 3, ALONG_EDGE_OR_PLANE = 4
/**
 *
 * @param brep BREP to check
 * @param edge edge to check
 * @param dirAtEdgeA the direction vector to check
 * @param faceNormal If dirAtEdgeA doesn't split a volume, but is along a face, the returned value depends on
 *     wether that faces normal1 points in the same direction as faceNormal
 * @returns INSIDE, OUTSIDE, COPLANAR_SAME or COPLANAR_OPPOSITE
 */
//function splitsVolumeEnclosingFaces(brep: B2, edge: Edge, dirAtEdgeA: V3, faceNormal: V3): int {
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
function splitsVolumeEnclosingFaces(brep: B2, canonEdge: Edge, dirAtEdgeA: V3, faceNormal: V3): int {
    assert(arguments.length == 4)
    assert(canonEdge == canonEdge.getCanon())
    //assert(p.equals(canonEdge.a))
    const edgeFaceInfos = brep.edgeFaces.get(canonEdge) as any[]
    assertf(() => edgeFaceInfos.length % 2 == 0)
    assertf(() => brep.edgeFaces)
    const faceInfo0 = edgeFaceInfos[0]
    const aDir1 = canonEdge.aDir.unit()
    const angleToCanon = (faceInfo0.inside.angleRelativeNormal(dirAtEdgeA, aDir1) + 2 * Math.PI + NLA_PRECISION) % (2 * Math.PI) - NLA_PRECISION
    const nearestFaceInfoIndex = edgeFaceInfos.findIndex(faceInfo => lt(angleToCanon, faceInfo.angle))
	const nearestFaceInfo = edgeFaceInfos[nearestFaceInfoIndex == -1 ? edgeFaceInfos.length - 1 : nearestFaceInfoIndex - 1]
    if (eq(nearestFaceInfo.angle, angleToCanon)) {
        //assert(false) todo
        const coplanarSame = nearestFaceInfo.normalAtCanonA.dot(faceNormal) > 0
        return coplanarSame ? COPLANAR_SAME : COPLANAR_OPPOSITE
    } else {
        return nearestFaceInfo.reversed ? INSIDE : OUTSIDE
    }
}
function splitsVolumeEnclosingFacesP(brep: B2, canonEdge: Edge, p: V3, pInside: V3, faceNormal: V3): int {
    assert(arguments.length == 5)
    assert(canonEdge == canonEdge.getCanon())
    //assert(p.equals(canonEdge.a))
    assertf(() => brep.edgeFaces)
    const edgeFaceInfos = brep.edgeFaces.get(canonEdge) as any[]
    assertf(() => edgeFaceInfos.length % 2 == 0)
    const pDir1 = canonEdge.tangentAt(canonEdge.curve.pointT(p)).unit()
    const faceInfoAngleFromPInsideNeg = faceInfo => {
        const faceInfoPDir = faceInfo.edge.getCanon() == faceInfo.edge ? pDir1 : pDir1.negated()
        const faceInfoInsideAtP = faceInfo.face.surface.normalP(p).cross(faceInfoPDir)
        const faceInfoAngleAtP = pInside.angleRelativeNormal(faceInfoInsideAtP, pDir1)
        return -((faceInfoAngleAtP + TAU + NLA_PRECISION) % TAU - NLA_PRECISION)
    }
    const nearestFaceInfo = edgeFaceInfos.withMax(faceInfoAngleFromPInsideNeg)
    if (eq0(faceInfoAngleFromPInsideNeg(nearestFaceInfo))) {
        //assert(false) todo
        const coplanarSame = nearestFaceInfo.normalAtCanonA.dot(faceNormal) > 0
        return coplanarSame ? COPLANAR_SAME : COPLANAR_OPPOSITE
    } else {
        return nearestFaceInfo.reversed ? OUTSIDE : INSIDE
    }
}
function splitsVolumeEnclosingFacesP2(brep: B2, canonEdge: Edge, p: V3, testCurve: Curve, curveT: number, dir: -1 | 1, faceNormal: V3): int {
    assert(canonEdge == canonEdge.getCanon())
    //assert(p.equals(canonEdge.a))
    assertf(() => brep.edgeFaces)
    const edgeFaceInfos = brep.edgeFaces.get(canonEdge) as any[]
    assertf(() => edgeFaceInfos.length % 2 == 0)
    const pDir1 = canonEdge.tangentAt(canonEdge.curve.pointT(p)).unit()
	const pInside = testCurve.tangentAt(curveT).times(dir)
	let minValue = 20, advanced = false, result = OUTSIDE
	for (const faceInfo of edgeFaceInfos) {
        const faceInfoPDir = faceInfo.edge.getCanon() == faceInfo.edge ? pDir1 : pDir1.negated()
        const faceInfoInsideAtP = faceInfo.face.surface.normalP(p).cross(faceInfoPDir)
        const faceInfoAngleAtP = pInside.angleRelativeNormal(faceInfoInsideAtP, pDir1)
	    const angle = (faceInfoAngleAtP + TAU + NLA_PRECISION) % TAU - NLA_PRECISION
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
		    const iscd = isCurve.at(isCurvePT).to(isCurve.at(isCurvePT + dir *dirFactor* eps)).dot(normVector)
		    const ecd = testCurve.at(curveT).to(testCurve.at(curveT + dir * eps)).dot(normVector)
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
function splitsVolumeEnclosingCone(brep: B2, p: V3, dir: V3) {
    const testPlane = P3.forAnchorAndPlaneVectors(p, dir, dir.getPerpendicular())
    const rays = []
    for (let k = 0; k < brep.faces.length; k++) {
        const face = brep.faces[k] as PlaneFace
        assertf(() => face instanceof PlaneFace)
        if (face.getAllEdges().some(edge => edge.a.like(p))) {
            if (testPlane.isParallelToPlane(face.surface.plane)) {
                if (face.pointsToInside(p, dir) != PointVsFace.OUTSIDE) {
                    return ALONG_EDGE_OR_PLANE
                }
            } else {
                const isLine = L3.fromPlanes(testPlane, face.surface.plane)
                const ps = planeFaceEdgeISPsWithPlane(face, isLine, testPlane)
                let i = 0
                while (i < ps.length) {
                    const a = ps[i++], b = ps[i++]
                    const out = a.p.like(p)
                    if (out || b.p.like(p)) {
                        const dir2 = out ? isLine.dir1 : isLine.dir1.negated()
                        const angle = (dir.angleRelativeNormal(dir2, testPlane.normal1) + 2 * Math.PI + NLA_PRECISION / 2) % (2 * Math.PI)
                        rays.push({angle: angle, out: out})
                    }
                }
            }
        }
    }
    rays.sort((a, b) => a.angle - b.angle)
    //console.log("testPlane", testPlane.toSource(), "rays", rays.toSource())

    if (eq0(rays[0].angle)) {
        return ALONG_EDGE_OR_PLANE
    } else {
        return rays[0].out ? OUTSIDE : INSIDE
    }
}
function splitsVolumeEnclosingCone2(brep: B2, p: V3, curve: Curve, curveT: number, fb: 1 | -1) {
	assert(curve.containsPoint(p))
	const dir = curve.tangentAt(curveT).times(fb)
	const testPlane = P3.forAnchorAndPlaneVectors(p, dir, dir.getPerpendicular())
	const rays = []
	const pFaces = brep.faces.filter(face => face.getAllEdges().some(edge => edge.a.like(p)))
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
	return brep.containsPoint(curve.at(curveT + fb * EPS), true) ? INSIDE : OUTSIDE
}
function fff(info: {face: Face, edge: Edge, normalAtCanonA: V3, inside: V3, reversed: boolean, angle: number}, surface: Surface): int {
    const canonA = info.edge.reversed ? info.edge.b : info.edge.a
    const surfaceNormalAtCanonA = surface.normalP(canonA)
    const dot = snap0(info.inside.dot(surfaceNormalAtCanonA))
    if (0 !== dot) {
        return 0 < dot ? OUTSIDE : INSIDE
    }
    if (surface.isCoplanarTo(info.face.surface)) {
        return 0 < info.normalAtCanonA.dot(surfaceNormalAtCanonA) ? COPLANAR_SAME : COPLANAR_OPPOSITE
    }
    assert(false)
}
function makeLink(values: any) {
	return 'viewer.html#' + Object.getOwnPropertyNames(values).map(name => {
			const val = values[name]
			return name + '=' + (typeof val == 'string' ? val : val.toSource())
		}).join(';')
}
declare function earcut(data: FloatArray, holeIndices: number[], dim: int): int[]
function triangulateVertices(normal: V3, vertices: V3[], holeStarts: int[]) {
	const absMaxDim = normal.maxAbsDim(), factor = sign(normal.e(absMaxDim))
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
function intersectionUnitCircleLine(a: number, b: number, c: number): {x1: number, y1: number, x2: number, y2: number} {
	assertNumbers(a, b, c)
	// TODO: disambiguate on a < b
	const term = sqrt(a * a + b * b - c * c)
	return {
		x1: (a * c + b * term) / (a * a + b * b),
		x2: (a * c - b * term) / (a * a + b * b),
		y1: (b * c - a * term) / (a * a + b * b),
		y2: (b * c + a * term) / (a * a + b * b)
	}
}
function intersectionUnitCircleLine2(a: number, b: number, c: number): [number, number][] {
	assertNumbers(a, b, c)
	// TODO: disambiguate on a < b
	// cf. pqFormula
	const termSqr = snap0(a * a + b * b - c * c)
	const term = sqrt(termSqr)
	if (termSqr < 0) {
		return []
	} else if (termSqr == 0) {
		return [    [   (a * c + b * term) / (a * a + b * b),
				        (b * c - a * term) / (a * a + b * b)]]
	} else {
		return [    [   (a * c + b * term) / (a * a + b * b),
						(b * c - a * term) / (a * a + b * b)],
					[   (a * c - b * term) / (a * a + b * b),
						(b * c + a * term) / (a * a + b * b)]]
	}
}
function intersectionCircleLine(a: number, b: number, c: number, r: number): {x1: number, x2: number, y1: number, y2: number} {
    assertNumbers(a, b, c, r)
    const term = sqrt(r * r * (a * a + b * b) - c * c)
    return {
        x1: (a * c + b * term) / (a * a + b * b),
        x2: (a * c - b * term) / (a * a + b * b),
        y1: (b * c - a * term) / (a * a + b * b),
        y2: (b * c + a * term) / (a * a + b * b)
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
function intersectionUnitHyperbolaLine(a: number, b: number, c: number): { x1: number, y1: number, x2: number, y2: number } {
    assertNumbers(a, b, c)
    const aa = a * a, bb = b * b, cc = c * c
    // TODO: disambiguate on a < b
    //var xTerm = sqrt(4*cc*aa-4*(bb-aa)*(-cc-bb))
    const xTerm = 2 * sqrt(bb * cc + bb * bb - aa * bb)
    const yTerm = sqrt(4 * cc * bb - 4 * (bb - aa) * (cc - aa))
    return {
        x1: (-2 * a * c + xTerm) / 2 / (bb - aa),
        x2: (-2 * a * c - xTerm) / 2 / (bb - aa),
        y1: (2 * b * c - yTerm) / 2 / (bb - aa),
        y2: (2 * b * c + yTerm) / 2 / (bb - aa)
    }
}


//function followAlgorithm2d(implicitCurve: (s: number, t: number) => number,
//                           start: V3,
//                           stepLength: number = 0.5,
//                           dids: (s: number, t: number) => number,
//                           didt: (s: number, t: number) => number,
//                           bounds: (s: number, t: number) => boolean,
//                           endp: V3 = start): {points: V3[], tangents: V3[]} {
//	assertNumbers(stepLength, implicitCurve(0, 0))
//	assertVectors(start)
//	//assert (!startDir || startDir instanceof V3)
//	const points = []
//	//let tangents = tangents || []
//	const tangents = []
//	assert (eq02(implicitCurve(start.x, start.y), 0.01), 'isZero(implicitCurve(startPoint.x, startPoint.y))')
//	const eps = stepLength / 32
//	let p = start, prevp = p
//	let i = 0
//	do {
//		const dfpdx = dids(p.x, p.y), dfpdy = didt(p.x, p.y)
//		let tangent = new V3(-dfpdy, dfpdx, 0)
//		const reversedDir = p.minus(prevp).dot(tangent) < 0
//		tangent = tangent.toLength(stepLength)
//		const tangentEndPoint = p.plus(tangent)
//		points.push(p)
//		tangents.push(tangent)
//		prevp = p
//		let newP = curvePoint(implicitCurve, tangentEndPoint, dids, didt)
//		if (newP.equals(p)) {
//			assertNever()
//		}
//		p = newP
//		assert(eq0(implicitCurve(p.x, p.y)))
//	} while (i++ < 1000 && (i < 4 || prevp.distanceTo(endp) > stepLength && bounds(p.x, p.y)))
//	assert(i != 1000)
//	//assert(bounds(p.x, p.y))
//	const end = (i < 4 || prevp.distanceTo(endp) > stepLength) ? p : endp
//	const endTangent = new V3(-didt(end.x, end.y), dids(end.x, end.y), 0).toLength(stepLength)
//	points.push(end)
//	tangents.push(endTangent)
//
//	//assert(points.length > 6)
//	// TODO gleichmäßige Verteilung der Punkte
//	return {points, tangents}
//}
function followAlgorithm2d(ic: MathFunctionR2_R,
                           start: V3,
                           stepLength: number = 0.5,
                           bounds: (s: number, t: number) => boolean,
                           endp: V3 = start): {points: V3[], tangents: V3[]} {
	assertNumbers(stepLength, ic(0, 0))
	assertVectors(start)
	//assert (!startDir || startDir instanceof V3)
	const points = []
	//let tangents = tangents || []
	const tangents = []
	assert (eq02(ic(start.x, start.y), 0.01), 'isZero(implicitCurve(startPoint.x, startPoint.y))')
	const eps = stepLength / 32
	let p = start, prevp = p, prevTangent = V3.O, broke = false
	let i = 0
	do {
		const dfpdx = ic.x(p.x, p.y), dfpdy = ic.y(p.x, p.y)
		let tangent = new V3(-dfpdy, dfpdx, 0)
		const reversedDir = p.minus(prevp).dot(tangent) < 0
		tangent = tangent.toLength(stepLength)
		if (prevTangent.dot(tangent) < 0) {
			// singularity
			const singularity = newtonIterate2d(ic.x, ic.y, p.x, p.y)
			if (eq0(ic(singularity.x, singularity.y)) && singularity.distanceTo(p) < abs(stepLength)) {
                // end on this point
                points.push(singularity)
                tangents.push(p.to(singularity))
                broke = true
                break
            } else {
				throw new Error()
			}
		}
		const tangentEndPoint = p.plus(tangent)
		points.push(p)
		tangents.push(tangent)
		prevp = p
        prevTangent = tangent
		const newP = curvePointMF(ic, tangentEndPoint)
		if (newP.equals(p)) {
			assertNever()
		}
		p = newP
		assert(eq0(ic(p.x, p.y)))
	} while (++i < 1000 && (i < 4 || prevp.distanceTo(endp) > stepLength && bounds(p.x, p.y)))
	assert(i != 1000)
	//assert(bounds(p.x, p.y))
    if (!broke) {
        const end = (i < 4 || prevp.distanceTo(endp) > stepLength) ? p : endp
        const endTangent = new V3(-ic.y(end.x, end.y), ic.x(end.x, end.y), 0).toLength(stepLength)
        points.push(end)
        tangents.push(endTangent)
    }

	//assert(points.length > 6)
	// TODO gleichmäßige Verteilung der Punkte
	return {points, tangents}
}
function followAlgorithm2dAdjustable(ic: MathFunctionR2_R,
                                     start: V3,
                                     stepLength: number = 0.5,
                                     bounds: (s: number, t: number) => boolean,
                                     endp: V3 = start): {points: V3[], tangents: V3[]} {
	assertNumbers(stepLength, ic(0, 0))
	assertVectors(start)
	//assert (!startDir || startDir instanceof V3)
	const points = []
	//let tangents = tangents || []
	const tangents = []
	assert (eq02(ic(start.x, start.y), 0.01), 'isZero(implicitCurve(startPoint.x, startPoint.y))')
	const eps = stepLength / 32
	let p = start, prevp = p
	let i = 0
	do {
		const dfpdx = ic.x(p.x, p.y), dfpdy = ic.y(p.x, p.y)
        const dfpdxx = ic.xx(p.x, p.y), dfpdyy = ic.yy(p.x, p.y), dfpdxy = ic.xy(p.x, p.y)
		const c2factor = abs((dfpdy ** 2 * dfpdxx - 2 * dfpdx * dfpdy * dfpdxy + dfpdx ** 2 * dfpdyy) /
								(dfpdx ** 2 + dfpdy ** 2) ** 2)
		const c2 = new V3(dfpdx, dfpdy, 0).times(c2factor)
		const s = 1 / 16 / c2.length()
		const tangent = new V3(-dfpdy, dfpdx, 0).unit()
		const reversedDir = p.minus(prevp).dot(tangent) < 0
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
	} while (i++ < 1000 && (i < 4 || prevp.distanceTo(endp) > stepLength) && bounds(p.x, p.y))
	assert(i != 1000)
	//assert(bounds(p.x, p.y))
	const end = (i < 4 || prevp.distanceTo(endp) > stepLength) ? p : endp
	const endTangent = new V3(-ic.y(end.x, end.y), ic.x(end.x, end.y), 0).toLength(stepLength)
	points.push(end)
	tangents.push(endTangent)

	//assert(points.length > 6)
	// TODO gleichmäßige Verteilung der Punkte
	return {points, tangents}
}
// both curves must be in the same s-t coordinates for this to make sense
function intersectionICurveICurve(
    iCurve1: (s: number, t: number) => number,
    startParams1: V3,
    endParams1: V3,
    startDir,
    stepLength: number,
    iCurve2: (s: number, t: number) => number) {

    assertNumbers(stepLength, iCurve1(0, 0), iCurve2(0, 0))
    assertVectors(startParams1, endParams1)
    assert (!startDir || startDir instanceof V3)
    const vertices = []
    assert (eq0(iCurve1(startParams1.x, startParams1.y)))
    stepLength = stepLength || 0.5
    const eps = 1e-5
    let p = startParams1, prevp = p // startDir ? p.minus(startDir) : p
    let i = 0
    while (i++ < 1000 && (i < 4 || p.distanceTo(endParams1) > 1.1 * stepLength)) {
        const fp = iCurve1(p.x, p.y)
        const dfpdx = (iCurve1(p.x + eps, p.y) - fp) / eps,
            dfpdy = (iCurve1(p.x, p.y + eps) - fp) / eps
        let tangent = new V3(-dfpdy, dfpdx, 0).toLength(stepLength)
        if (p.minus(prevp).dot(tangent) < 0) tangent = tangent.negated()
        prevp = p
        p = curvePoint(iCurve1, p.plus(tangent))
        vertices.push(p)
    }
    // TODO gleichmäßige Verteilung der Punkte
    return vertices

}
function intersectionICurveICurve2(iCurve1, loopPoints1, iCurve2) {
    let p = loopPoints1[0], val = iCurve2(p.x, p.y), lastVal
    const iss = []
    for (let i = 0; i < loopPoints1.length; i++) {
        lastVal = val
        p = loopPoints1[i]
        val = iCurve2(p)
        if (val * lastVal <= 0) { // TODO < ?
            iss.push(newtonIterate2d(iCurve1, iCurve2, p.x, p.y))
        }
    }
    return iss
}

function intersectionPCurveISurface(parametricCurve, searchStart, searchEnd, searchStep, implicitSurface) {
    assertNumbers(searchStart, searchEnd, searchStep)
    const iss = []
    let val = implicitSurface(parametricCurve(searchStart)), lastVal
    for (let t = searchStart + searchStep; t <= searchEnd; t += searchStep) {
        lastVal = val
        val = implicitSurface(parametricCurve(t))
        if (val * lastVal <= 0) {
            iss.push(newtonIterate1d(t => implicitSurface(parametricCurve(t)), t))
        }
    }
    return iss
}
function intersectionICurvePSurface(f0, f1, parametricSurface) {

}



