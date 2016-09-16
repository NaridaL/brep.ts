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


var eps = 1e-5
/**
 * @constructor
 * @augments {Transformable}
 * @returns {B2}
 */
var B2 = NLA.defineClass('B2', Transformable,
	/**
	 *
	 * @param {Face[]} faces
	 * @param infiniteVolume
	 * @param generator
	 */
	function (faces, infiniteVolume, generator) {
		this.faces = faces
		assert(faces.every(f => f instanceof B2.Face), () => 'faces.every(f => f instanceof B2.Face)\n' + this.toString())
		this.infiniteVolume = !!infiniteVolume
		this.generator = generator
	},
	/** @lends {B2.prototype} */ {
		calculateVolume: function () {
			return this.faces.map(face => face.zDirVolume()).sum()
		},

		/**
		 *
		 * @returns {Mesh}
		 */
		toMesh: function () {
			var mesh = new GL.Mesh({triangles: true, normals: true, lines: true})
			mesh.faceIndexes = new Map()
			this.faces.forEach((face, i) => {
				face.addToMesh(mesh)
			})
			mesh.compile()
			return mesh
		},

		/**
		 *
		 * @param {B2} brep2
		 * @returns {B2}
		 */
		minus: function (brep2) {
			return this.intersection(brep2.flipped(), true, true)
		},

		/**
		 *
		 * @param {B2} brep2
		 * @returns {B2}
		 */
		plus: function (brep2) {
			return this.flipped().intersection(brep2.flipped(), true, true).flipped()
		},

		/**
		 *
		 * @param {B2} brep2
		 * @returns {B2}
		 */
		xor: function (brep2) {
			// TODO
			assert(false, "Not implemented yet")
		},

		/**
		 *
		 * @param {B2} brep
		 * @returns {boolean}
		 */
		equals: function (brep) {
			return this.faces.length == brep.faces.length &&
				this.faces.every((face) => brep.faces.some((face2) => face.equals(face2)))
		},

		/**
		 *
		 * @returns {string}
		 */
		toString: function () {
			return `new B2([\n${this.faces.join(',\n').replace(/^/gm, '\t')}])`
		},

		/**
		 *
		 * @returns {string}
		 */
		toSource: function () {
			return this.generator || `new B2([\n${this.faces.map(face => face.toSource()).join(',\n').replace(/^/gm, '\t')}])`
		},

		/**
		 *
		 * @param newFaces
		 * @param loops
		 * @param oldFace
		 * @returns {boolean}
		 */
		assembleFacesFromLoops: function (newFaces, loops, oldFace) {
			var surface = oldFace.surface
			console.log(loops.map(loop=> loop.join('\n')).join('\n\n'))
			var topLevelLoops = []
			function placeRecursively(newLoop, arr) {
				if (arr.length == 0) {
					arr.push(newLoop)
				} else {
					var sl = arr.find(subloop => {
						var contains = surface.edgeLoopContainsPoint(subloop.edges, newLoop.edges[0].a)
						console.log(newLoop.edges[0].a.sce, newLoop.ccw ? contains : !contains)
						return contains
					})
					console.log("here", sl)
					if (sl) {
						placeRecursively(newLoop, sl.subloops)
					} else {
						// newLoop isnt contained by any other newLoop
						for (var i = arr.length; --i >= 0; ) {
							var subloop = arr[i]
							//console.log("cheving subloop", surface.edgeLoopContainsPoint(newLoop.edges, subloop.edges[0].a))
							if (surface.edgeLoopContainsPoint(newLoop.edges, subloop.edges[0].a)) {
								newLoop.subloops.push(subloop)
								arr.splice(i, 1) // remove it
							}
						}
						arr.push(newLoop)
					}
				}
			}
			function newFacesRecursive(loop) {
				var face = new oldFace.constructor(oldFace.surface, loop.edges, loop.subloops.map(sl => sl.edges))
				loop.subloops.forEach(sl => sl.subloops.forEach(tlLoop => newFacesRecursive(tlLoop)))
				console.log(face.toString(), face.holes.toString())
				newFaces.push(face)
			}
			loops.forEach(loop => {
				var ccw = surface.edgeLoopCCW(loop)
				console.log('CCW', ccw)
				placeRecursively({edges: loop, ccw: ccw, subloops: []}, topLevelLoops)
				console.log(topLevelLoops)
			})
			if (topLevelLoops[0].ccw) {
				topLevelLoops.forEach(tlLoop => newFacesRecursive(tlLoop))
				return false
			} else {
				assert(false)
				newFacesRecursive({edges: oldFace.edges, ccw: true, subloops: topLevelLoops})
				return true
			}
		},

		/**
		 * Rightmost next segment doesn't work, as the correct next segment isn't obvious from the current corner
		 * alone.
		 * (at least, not without extensive pre-analysos on the face edges, which shouldn't be necessary, as the
		 * correct new faces are defined by the new edges already.) Leftmost edge should work. Holes which touch the
		 * edge of the face will be added to the face contour.
		 *
		 * New segments will always be part left-er than exisiting ones, so no special check is required.
		 *
		 * @param oldFaces
		 * @param edgeLooseSegments
		 * @param faceMap
		 * @param newFaces
		 */
		reconstituteFaces: function (oldFaces, edgeLooseSegments, faceMap, newFaces) {
			// reconstitute faces
			var insideEdges = []
			oldFaces.forEach(face => {
				var allOldFaceEdges = face.getAllEdges()
				console.log('reconstituting face', face.toString())
				var els = face.edges.map(edge => edgeLooseSegments.get(edge) && edgeLooseSegments.get(edge).join('\n')).map(s => '\n\n'+s).join()
				console.log('edgeLooseSegments', els)
				var newSegments = faceMap.get(face)
				if (!newSegments) {
					face.insideOutside = 'undecided'
				} else {
					face.insideOutside = 'part'
					console.log('newSegments\n', newSegments.map(e=>e.toString()).join('\n'))
					var loops = []
					var currentEdge
					var edgeCond = face instanceof B2.PlaneFace
						? (edge => edge.a.like(currentEdge.b))
						: (edge => (edge.curve == currentEdge.curve // TODO: ??
										? NLA.equals(edge.aT, currentEdge.bT)
										: edge.a.like(currentEdge.b)))
					var getNextStart = function () {
						return newSegments.find(edge => !edge.visited) || allOldFaceEdges.find(edge => !edge.visited)
					}
					allOldFaceEdges.forEach(edge => edge.visited = false)
					while (currentEdge = getNextStart()) {
						var cancelLoop = false
						var startEdge = currentEdge, edges = [], i = 0
						var looseLoop = true // uses only new segments edges
						do {
							currentEdge.visited = true
							console.log('currentEdge', currentEdge.b.sce, currentEdge.toSource())
							edges.push(currentEdge)
							// find next edge
							var possibleLooseEdges = newSegments.filter(edge => edgeCond(edge)), looseSegments
							var possibleLooseEdgesCount = possibleLooseEdges.length
							allOldFaceEdges.forEach(edge => (looseSegments = edgeLooseSegments.get(edge)) && looseSegments.forEach(
									edge => edgeCond(edge) && possibleLooseEdges.push(edge)))
							var noSegments = false
							if (possibleLooseEdgesCount == possibleLooseEdges.length) {
								allOldFaceEdges.forEach(edge => edgeCond(edge) && possibleLooseEdges.push(edge))
								noSegments = true
							}
							var index = possibleLooseEdges.indexWithMax((edge, index) => currentEdge.bDir.angleRelativeNormal(edge.aDir, face.surface.normalAt(currentEdge.b)) + (index < possibleLooseEdgesCount) * NLA.PRECISION)
							//console.log('possibleLooseEdges\n', possibleLooseEdges.map(e=>e.toString()).join('\n'), allOldFaceEdges.find(edge => edgeCond(edge)), index, possibleLooseEdges[index])
							// TODO assert(possibleLooseEdges.length < 2)
							currentEdge = possibleLooseEdges[index]
							if (currentEdge.visited) {
								console.log("breaking")
								break
							}
							if (index < possibleLooseEdgesCount) {
								possibleLooseEdges.forEach(possibleLooseEdge => possibleLooseEdge.isCoEdge(currentEdge) && (possibleLooseEdge.visited = true))
							} else {
								looseLoop = false
								noSegments && insideEdges.push(currentEdge)
							}
							assert(currentEdge)
							assert(currentEdge != startEdge)
						} while (++i < 200)
						if (20 == i) {
							assert(false, "too many")
						}
						if (currentEdge == startEdge) {
							console.log('finished loop')
							loops.push(edges)
						} else {
							insideEdges.removeAll(edges)
						}
					}
					if (this.assembleFacesFromLoops(newFaces, loops, face)) {
						insideEdges.pushAll(face.edges)
					}
				}
			})
			while (insideEdges.length != 0) {
				var edge = insideEdges.pop()
				var adjoiningFaces = facesWithEdge(edge, oldFaces)
				adjoiningFaces.forEach(info => {
					if (info.face.insideOutside == 'undecided') {
						info.face.insideOutside = 'inside'
						insideEdges.push.apply(insideEdges, info.face.edges)
					}
				})
			}
			oldFaces.forEach(face => {
				if (face.insideOutside == 'inside') {
					newFaces.push(face)
				}
			})
		},

		reconstituteCoplanarFaces: function (likeSurfacePlanes, edgeLooseSegments, faceMap, newFaces) {
			likeSurfacePlanes.forEach(faceGroup => {
				// calculate total contours
				let surface = faceGroup[0].surface, bag = []
				faceGroup.forEach(face => {
					Array.prototype.push.apply(bag, faceMap(face))
					face.getAllEdges().forEach(edge => {
						let edgeSubSegments
						if (edgeSubSegments = edgeLooseSegments.get(edge)) {
							Array.prototype.push.apply(bag, edgeSubSegments)
						} else {
							bag.push(edge)
						}
					})
				})
				let currentEdge, loops = []
				while (currentEdge = bag.find(edge => !edge.visited)) {
					let path = []
					do {
						currentEdge.visited = true
						path.push(currentEdge)
						let possibleNextEdges = bag.filter(edge => currentEdge.b.like(edge.a))
						// lowest angle, i.e. the right-most next edge
						let nextEdgeIndex = possibleNextEdges.indexWithMax((edge, index) => -currentEdge.bDir.angleRelativeNormal(edge.aDir, surface.normalAt(currentEdge.b)))
						currentEdge = possibleNextEdges[nextEdgeIndex]
					} while (!currentEdge.visited)
					let startIndex = path.find(currentEdge)
					if (-1 != startIndex) {
						loops.push(path.slice(startIndex))
					}
				}
			})
		},

		/**
		 *
		 * @param edgeMap
		 * @returns {Map}
		 */
		getLooseEdgeSegmentsInsideDirs: function (edgeMap) {
			var edgeLooseSegments = new Map()
			console.log("edgeMap", edgeMap)
			edgeMap.forEach((pointInfos, baseEdge) => {
				// TODO: make sure it works with loop
				// TODO: dont generate unnecessarry segments
				var looseSegments = []
				if (!baseEdge.reversed) {
					pointInfos.sort((a, b) => NLA.snapTo(a.edgeT - b.edgeT, 0) || (b.insideDir.dot(baseEdge.curve.dir1) - a.insideDir.dot(baseEdge.curve.dir1)))
				} else {
					pointInfos.sort((b, a) => NLA.snapTo(a.edgeT - b.edgeT, 0) || (b.insideDir.dot(baseEdge.curve.dir1) - a.insideDir.dot(baseEdge.curve.dir1)))
				}
				var startP = baseEdge.a, startDir = baseEdge.aDir, startT = baseEdge.aT, startInfo
				for (var i = 0; i < pointInfos.length; i++) {
					var info = pointInfos[i]
					assert(info.insideDir, info.toSource())
					console.log("info.insideDir.dot(baseEdge.curve.dir1)", info.insideDir.dot(baseEdge.curve.dir1))
					// ignore start, end and repeating points
					if (NLA.equals(info.edgeT, baseEdge.bT) || NLA.equals(info.edgeT, baseEdge.aT) || NLA.equals(info.edgeT, startT)) {
						continue
					}
					var pDir = baseEdge.tangentAt(info.edgeT)
					// add segment only if insideDir points backwards
					if (info.insideDir.dot(baseEdge.curve.dir1) < 0) {
						looseSegments.push(new baseEdge.constructor(baseEdge.curve, startP, info.p, startT, info.edgeT, null, startDir, pDir))
					}
					startP = info.p
					startT = info.edgeT
					startInfo = info
					startDir = pDir
				}
				if (startInfo && startInfo.insideDir.dot(baseEdge.curve.dir1) > 0) {
					looseSegments.push(new baseEdge.constructor(baseEdge.curve, startP, baseEdge.b, startT, baseEdge.bT, null, startDir, baseEdge.bDir))
				}
				if (baseEdge.b.like(V3(0, 4, 1))) {
					console.log('pointInfos', baseEdge.reversed, baseEdge.sce, pointInfos.map(pi => '\n' + pi.toSource()).join())
					console.log(looseSegments)
				}
				edgeLooseSegments.set(baseEdge, looseSegments)
			})
			return edgeLooseSegments;
		},
		/**
		 *
		 * @param edgeMap
		 * @returns {Map}
		 */
		getLooseEdgeSegments: function (edgeMap) {
			var edgeLooseSegments = new Map()
			console.log("edgeMap", edgeMap)
			edgeMap.forEach((pointInfos, baseEdge) => {
				// TODO: make sure it works with loop
				// TODO: dont generate unnecessarry segments
				var looseSegments = []
				if (!baseEdge.reversed) {
					pointInfos.sort((a, b) => NLA.snapTo(a.edgeT - b.edgeT, 0) || (b.insideDir.dot(baseEdge.curve.dir1) - a.insideDir.dot(baseEdge.curve.dir1)))
				} else {
					pointInfos.sort((b, a) => NLA.snapTo(a.edgeT - b.edgeT, 0) || (b.insideDir.dot(baseEdge.curve.dir1) - a.insideDir.dot(baseEdge.curve.dir1)))
				}
				var startP = baseEdge.a, startDir = baseEdge.aDir, startT = baseEdge.aT, startInfo
				for (var i = 0; i < pointInfos.length; i++) {
					var info = pointInfos[i]
					// ignore start, end and repeating points
					if (NLA.equals(info.edgeT, baseEdge.bT) || NLA.equals(info.edgeT, baseEdge.aT) || NLA.equals(info.edgeT, startT)) {
						continue
					}
					var pDir = baseEdge.tangentAt(info.edgeT)
					looseSegments.push(new baseEdge.constructor(baseEdge.curve, startP, info.p, startT, info.edgeT, null, startDir, pDir))
					startP = info.p
					startT = info.edgeT
					startInfo = info
					startDir = pDir
				}
				looseSegments.push(new baseEdge.constructor(baseEdge.curve, startP, baseEdge.b, startT, baseEdge.bT, null, startDir, baseEdge.bDir))
				edgeLooseSegments.set(baseEdge, looseSegments)
			})
			return edgeLooseSegments;
		},

		/**
		 *
		 * @param brep2
		 * @param buildThis
		 * @param buildBREP2
		 */
		intersection: function (brep2, buildThis, buildBREP2, buildCoplanar) {
			var faceMap = new Map(), edgeMap = new Map()

			let likeSurfaceFaces = []

			this.faces.forEach(face => {
				//console.log('face', face.toString())
				brep2.faces.forEach(face2 => {
					//console.log('face2', face2.toString())
					face.doo(face2, this, brep2, faceMap, edgeMap, likeSurfaceFaces)
				})
			})
			var newFaces = []

			/*
			TODO:
			faceMap.forEach((faceLooses, face) => {
				faceLooses.forEach(edge => {
					face.edges.forEach(faceEdge => {
						var edgeT = faceEdge.getEdgeT(edge.a)
						if (undefined !== edgeT) {
							console.log("WAARGH", edge.a.$, faceEdge.toString(), edgeT)
							NLA.mapAdd(edgeMap, faceEdge, {edgeT: edgeT, p: edge.a})
						}
					})
				})
			})
			*/
			var edgeLooseSegments = this.getLooseEdgeSegments(edgeMap);

			buildThis && this.reconstituteFaces(this.faces, edgeLooseSegments, faceMap, newFaces, brep2.infiniteVolume)
			buildBREP2 && this.reconstituteFaces(brep2.faces, edgeLooseSegments, faceMap, newFaces, this.infiniteVolume)
			//buildCoplanar && this.reconstituteCoplanarFaces(likeSurfaceFaces, edgeLooseSegments, faceMap, newFaces, this.infiniteVolume, brep2.infiniteVolume)
			if (0 == newFaces.length) {
				return null
			} else {
				return new B2(newFaces, this.infiniteVolume && brep2.infiniteVolume)
			}
		},

		/**
		 *
		 * @param {M4} m4
		 * @param {string} doo
		 * @return {B2}
		 */
		transform: function (m4, doo) {
			return new B2(this.faces.map(f => f.transform(m4)), this.infiniteVolume, this.generator && this.generator + doo)
		},

		/**
		 * @return {B2}
		 */
		flipped: function () {
			return new B2(this.faces.map(f => f.flipped()), !this.infiniteVolume, this.generator && this.generator + ".flipped()")
		}
	}
)
/**
 *
 * @param {Surface} surface
 * @param {Edge[]} contour
 * @param {Edge[][]=} holes
 * @constructor
 */
B2.Face = function (surface, contour, holes) {
	this.assertChain(contour)
	//assert(surface.edgeLoopCCW(contour), surface.toString()+contour.join("\n"))
	assert(contour.every(f => f instanceof B2.Edge), 'contour.every(f => f instanceof B2.Edge)' + contour.toSource())
	// contour.forEach(e => !surface.containsCurve(e.curve) && console.log("FAIL:"+surface.distanceToPoint(e.curve.anchor)))
	contour.forEach(e => assert(surface.containsCurve(e.curve), e + surface))
	holes && holes.forEach(hole => this.assertChain(hole))
	//holes && holes.forEach(hole => assert(!surface.edgeLoopCCW(hole)))
	assert (!holes || holes.constructor == Array, holes && holes.toString())
	this.surface = surface
	this.edges = contour // TODO refactor to contour
	this.holes = holes || []
	this.id = globalId++

}
B2.Face.prototype = NLA.defineObject(Transformable.prototype, {
	transform: function (m4) {
		var newEdges = this.edges.map(e => e.transform(m4))
		var newHoles = this.holes.map(hole => hole.map(e => e.transform(m4)))
		return new this.constructor(this.surface.transform(m4), newEdges, newHoles)
	},
	assertChain: function (edges) {
		edges.forEach((edge, i) => {
			var j = (i + 1) % edges.length
			assert(edge.b.like(edges[j].a), `edges[${i}].b != edges[${j}].a (${edges[i].b.sce} != ${edges[j].a.sce})`)
		})
	},
	flipped: function () {
		var newEdges = this.edges.map(e => e.flipped()).reverse()
		var newHoles = this.holes.map(hole => hole.map(e => e.flipped()).reverse())
		return new this.constructor(this.surface.flipped(), newEdges, newHoles)
	},
	toString: function () {
		return `new ${this.name}(${this.surface}, [${this.edges.map(e => '\n\t' + e).join()}]`
			+`${this.holes.map(hole => '\n\t\thole: ' + hole.join())})`
	},
	toSource: function () {
		return `new ${this.name}(${this.surface.toSource()}, [${this.edges.map(e => '\n\t' + e.toSource()).join(',')}], [${
			this.holes.map(hole => '['+hole.map(e => '\n\t' + e.toSource()).join(',')+']').join(',')}])`
	},
	equals: function (face) {
//TODO		assert(false)
		var edgeCount = this.edges.length

		return this.surface.equalsSurface(face.surface) &&
				this.edges.length == face.edges.length &&
				NLA.arrayRange(0, edgeCount, 1)
					.some(offset => this.edges.every((edge, i) => edge.equals(face.edges[(offset + i) % edgeCount])))
	},
	likeFace: function (face2) {
		function loopsLike(a, b) {
			return a.length == b.length &&
				NLA.arrayRange(0, a.length, 1)
					.some(offset => a.every((edge, i) => edge.likeEdge(b[(offset + i) % a.length])))

		}
		assertInst(B2.Face, face2)
		return this.surface.equalsSurface(face2.surface)
			&& this.holes.length == face2.holes.length
			&& loopsLike(this.edges, face2.edges)
			&& this.holes.every(hole => face2.holes.some(hole2 => loopsLike(hole, hole2)))
	},
	getAllEdges: function () {
		return Array.prototype.concat.apply(this.edges, this.holes)
	},
	addEdgeLines: function (mesh) {
		assert(false, "buggy, fix")
		var vertices = this.edges.map(edge => edge.getVerticesNo0()).concatenated(), mvl = mesh.vertices.length
		for (var i = 0; i < vertices.length; i++) {
			mesh.vertices.push(vertices[i])
			mesh.lines.push(mvl + i, mvl + (i + 1) % vertices.length)

		}
	},
	inB2: function () {
		return new B2([this])
	},

	/**
	 *
	 * @param {V3} p
	 * @returns {boolean}
	 */
	containsPoint: function (p) {
		assertVectors (p)
		return this.surface.edgeLoopContainsPoint(this.edges, p)
			&& !this.holes.some(hole => this.surface.edgeLoopContainsPoint(hole, p))
	},

	/**
	 *
	 * @param {L3} line
	 * @returns {number}
	 */
	intersectsLine: function (line) {
		assertInst(L3, line)
		let containedIntersections = this.surface.isPointsWithLine(line).filter(p => this.containsPoint(p))
		let nearestPoint = containedIntersections.withMax(p => -line.pointLambda(p))

		return nearestPoint ? line.pointLambda(nearestPoint) : NaN
	},
	toMesh: function () {
		let mesh = new GL.Mesh({triangles: true, normals: true, lines: true})
		mesh.faceIndexes = new Map()
		this.addToMesh(mesh)
		mesh.compile()
		return mesh
	},
	constructor: B2.Face
})


B2.Edge = function () {}
B2.Edge.prototype = NLA.defineObject(Transformable.prototype, {
	toString: function (f) {
		return `new ${this.name}(${this.curve.toString(f)}, ${this.a}, ${this.b}, ${this.aT}, ${this.bT}, null, ${this.aDir}, ${this.bDir})`
	},
	getIntersectionsWithPlane: function (p) {
		assert(false, this.name + '.isTsWithPlane')
	},
	colinearToLine: function (line) {
		return this.curve instanceof L3 && this.curve.isColinearTo(line)
	},
	tValueInside: function (t) {
		return this.aT < this.bT 
			? NLA.lt(this.aT, t) && NLA.lt(t, this.bT)
			: NLA.lt(this.bT, t) && NLA.lt(t, this.aT)
	}
})
/**
 *
 * @param {Array.<B2.Edge>} loop
 * @returns {boolean}
 */
B2.Edge.isLoop = function (loop) {
	return loop.every((edge, i) => edge.b.like(loop[(i + 1) % loop.length].a))
}
B2.Edge.edgesIntersect = function (e1, e2) {
	assertInst(B2.Edge, e1, e2)
	if (e1.hlol > e2.hlol) {
		[e2, e1] = [e1, e2]
	}
	let sts = e1.curve.isInfosWithCurve(e2.curve)
	// console.log(sts.map(SCE), e1.aT, e1.bT, e2.aT, e2.bT)
	return sts.some(
		/// (  e1.aT < tThis < e1.bT  )  &&  (  e2.aT < tOther < e2.bT  )
		({tThis, tOther}) => e1.tValueInside(tThis) && e2.tValueInside(tOther))
}

B2.PlaneFace = NLA.defineClass('B2.PlaneFace', B2.Face,
	function (planeSurface, contour, holes, name) {
		B2.Face.call(this, planeSurface, contour, holes, name)
		assertInst(PlaneSurface, planeSurface)
	},
	{
		calculateArea: function () {

		},

		addToMesh: function (mesh) {
			var mvl = mesh.vertices.length
			var normal = this.surface.plane.normal
			var vertices = this.edges.map(edge => edge.getVerticesNo0()).concatenated()
			for (var i = 0; i < vertices.length; i++) { mesh.lines.push(mvl + i, mvl + (i + 1) % vertices.length) }
			var holeStarts = []
			this.holes.forEach(hole => {
				holeStarts.push(vertices.length)
				vertices.pushAll(hole.map(edge => edge.getVerticesNo0()).concatenated())
			})
			var triangles = triangulateVertices(normal, vertices, holeStarts).map(index => index + mvl)
			mesh.faceIndexes.set(this, {start: mesh.triangles.length, count: triangles.length})
			Array.prototype.push.apply(mesh.vertices, vertices)
			Array.prototype.push.apply(mesh.triangles, triangles)
			Array.prototype.push.apply(mesh.normals, NLA.arrayFromFunction(vertices.length, i => normal))
		},
		containsPoint: function (p) {
			assertVectors (p)
			return this.surface.edgeLoopContainsPoint(this.edges, p)
			&& !this.holes.some(hole => this.surface.edgeLoopContainsPoint(hole, p))
		},
		intersectsLine: function (line) {
			assertInst(L3, line)
			var lambda = line.intersectWithPlaneLambda(this.surface.plane)
			if (!Number.isFinite(lambda)) {
				return NaN
			}
			var inside = this.containsPoint(line.at(lambda))
			return inside ? lambda : NaN
		},
		withHole: function (holeEdges) {
			return new B2.PlaneFace(this.surface, this.edges, [holeEdges])
		},
	/*
		doo2: function (face2, thisBrep, face2Brep, faceMap, edgeMap, removeCoplanarSame, removeCoplanarOpposite) {
			if (face2 instanceof B2.RotationFace) {
				// get intersection
				var newCurves = []
				// get intersections of newCurve with other edges of face and face2
				var pss = new Map(), ps1count = 0, ps2count = 0
				this.edges.forEach(edge => {
					var iss = edge.getIntersectionsWithISurface(face2.surface)
					//console.log('iss',iss, edge.toString())
					for (var i = 0; i < iss.length; i++) {
						var edgeT = iss[i], p = edge.curve.at(edgeT), newCurveT
						var newCurve = newCurves.find(curve => !isNaN(newCurveT = curve.pointLambda(p)))
						if (!newCurve) {
							newCurves.push(newCurve = new CurvePI(this.surface, face2.surface, p))
							newCurveT = newCurve.pointLambda(p)
							pss.set(newCurve, {ps1: [], ps2: [],
								thisDir: face2.surface.normalAt(p).cross(this.surface.normalAt(p)).dot(newCurve.tangentAt(newCurveT)) > 0})
							/*console.log("NEWCURVE", p.$, face2.surface.normalAt(p).cross(this.surface.normalAt(p)).$, 'nct', newCurve.tangentAt(newCurveT).$,
								face2.surface.normalAt(p).cross(this.surface.normalAt(p)).dot(newCurve.tangentAt(newCurveT)) > 0,
								'newCurveT',newCurveT)
						}
						var ov = edge.tangentAt(edgeT).cross(this.surface.normalAt(p))
						var ct = newCurve.tangentAt(newCurveT)
						console.log("ov", p.$,edge.tangentAt(edgeT).$, this.surface.normalAt(p).$, ov.$, ct.$, ov.dot(ct) > 0)
						if (ov.dot(ct) > 0) ct = ct.negated()
						pss.get(newCurve).ps1.push({p: p, insideDir: ct, t: newCurveT, edge: edge, edgeT: edgeT})
						ps1count++
					}
				})
				//console.log(new CurvePIEdge(newCurve, ps[0], ps[1], ts[0], ts[1]))
				face2.edges.forEach(edge => {
					var iss = edge.getIntersectionsWithPSurface(this.surface)
				})
				if (ps1count == 0 && ps2count == 0) {
					// faces to not intersect
					return
				}
				newCurves.forEach((newCurve, key) => {
					var {ps1, ps2, thisDir} = pss.get(newCurve)
					var segments = (newCurve instanceof L3 )
						?
						: newCurve.getIntersectionSegments(ps1, ps2)
					console.log('ps', ps1.toSource(), ps2.toSource())
					// TODO: getCanon() TODO TODO TODO
					console.log('segments', segments.toSource())
					ps1.forEach(ps => ps.used && mapAdd(edgeMap, ps.edge, ps))
					ps2.forEach(ps => ps.used && mapAdd(edgeMap, ps.edge, ps))
					segments.forEach(segment => {
						console.log('segment', segment.toString())
						mapAdd(faceMap, this, thisDir ? segment : segment.flipped())
						mapAdd(faceMap, face2, thisDir ? segment.flipped() : segment)
					})
				})
				console.log('faceMap', faceMap)
			} else if (face2 instanceof B2.PlaneFace) {
				this.dooPlaneFace(face2, thisBrep, face2Brep, faceMap, edgeMap, removeCoplanarSame, removeCoplanarOpposite)
			}
			/*
			 // get intersection
			 var newCurve = this.surface.isCurvesWithSurface(face2.surface)
			 // get intersections of newCurve with other edges of face and face2
			 var ps1 = []
			 this.edges.forEach(edge => {
			 var iss = edge.getIntersectionsWithISurface(face2.surface)
			 for (var i = 0; i < iss.length; i++) {
			 var p = edge.curve.at(iss[i])
			 var ov = edge.pointTangent(p).cross(this.surface.normalAt(p))
			 var ct = newCurve.pointTangent(p)
			 //console.log("ov", p.$,edge.pointTangent(p).$, this.surface.normalAt(p).$, ov.$, ct.$)
			 if (ov.dot(ct) > 0) ct = ct.negated()
			 ps1.push({p: p, insideDir: ct, t: NaN, edge: edge, edgeT: iss[i]})
			 }
			 })
			 //console.log(new CurvePIEdge(newCurve, ps[0], ps[1], ts[0], ts[1]))
			 var ps2 = []
			 face2.edges.forEach(edge => {
			 var iss = edge.getIntersectionsWithPSurface(this.surface)
			 })
			 if (ps1.length == 0 && ps2.length == 0) {
			 // faces to not intersect
			 return
			 }
			 console.log(ps1.toSource(), ps2)
			 var segments = newCurve.getIntersectionSegments(ps1, ps2)
			 // TODO: getCanon()
			 ps1.forEach(ps => ps.used && mapAdd(edgeMap, ps.edge, ps))
			 segments.forEach(segment => {
			 mapAdd(faceMap, this, segment.flipped())
			 mapAdd(faceMap, face2, segment)
			 })
			 }
		},*/
		doo: function (face2, thisBrep, face2Brep, faceMap, edgeMap, likeSurfaceFaces) {
			if (face2 instanceof B2.RotationFace) {
				if (this.surface.isCoplanarTo(face2.surface)) { return }

				// get intersections
				var newCurves = face2.surface.isTsWithSurface(this.surface)

				console.log(newCurves)
				if (newCurves.length == 0) {
					return
				}

				// get intersections of newCurves with other edges of face and face2
				var pss1 = getFacePlaneIntersectionSs2(thisBrep, this, newCurves, face2.surface, true, false)
				var pss2 = getFacePlaneIntersectionSs2(face2Brep, face2, newCurves, this.surface, false, false)

				newCurves.forEach((newCurve, i) => {
					var ps1 = pss1[i], ps2 = pss2[i]
					if (ps1.length == 0 || ps2.length == 0) { return }

					var ps = ps1.length != 0 ? ps1[0] : ps2[0]
					var thisDir = !(face2.surface.normalAt(ps.p).cross(this.surface.normalAt(ps.p)).dot(newCurve.tangentAt(ps.t)) > 0)

					var in1 = ps1[0].insideDir.dot(newCurve.tangentAt(ps1[0].t)) < 0
					var in2 = ps2[0].insideDir.dot(newCurve.tangentAt(ps2[0].t)) < 0
					if(newCurve.debugToMesh) {
						dMesh = new GL.Mesh()
						newCurve.debugToMesh(dMesh, 'curve2')
						dMesh.compile()
						console.log(dMesh)
					}
					console.log('iscurve', newCurve.toString(), newCurve.tangentAt(ps1[0].t).sce)
					console.log('ps1\n', ps1.map(m => m.toSource()).join('\n'), '\nps2\n', ps2.map(m => m.toSource()).join('\n'))
					var segments = newCurve instanceof L3
						? getBlug(ps1, ps2, newCurve)
						: getIntersectionSegments(ps1, ps2, in1, in2, B2.PCurveEdge, newCurve)
					// TODO: getCanon() TODO TODO TODO
					console.log('segments', segments.toSource())
					ps1.forEach(ps => ps.used && mapAdd(edgeMap, ps.edge, ps))
					ps2.forEach(ps => ps.used && mapAdd(edgeMap, ps.edge, ps))
					segments.forEach(segment => {
						console.log('segment', segment.toString())
						mapAdd(faceMap, this, thisDir ? segment : segment.flipped())
						mapAdd(faceMap, face2, thisDir ? segment.flipped() : segment)
					})
				})
				console.log('faceMap', faceMap)
			} else if (face2 instanceof B2.PlaneFace) {
				this.dooPlaneFace(face2, thisBrep, face2Brep, faceMap, edgeMap, likeSurfaceFaces)
			}
		},

		/**
		 *
		 * @param face2
		 * @param thisBrep
		 * @param face2Brep
		 * @param faceMap
		 * @param edgeMap
		 * @param likeSurfaceFaces
		 */
		dooPlaneFace: function (face2, thisBrep, face2Brep, faceMap, edgeMap, likeSurfaceFaces) {
			assertInst(B2.PlaneFace, face2)
			var face = this
			// get intersection
			var thisPlane = this.surface.plane, face2Plane = face2.surface.plane
			if (thisPlane.isParallelToPlane(face2Plane)) {
				if (thisPlane.like(face2Plane)) {
					// normal same and same location in space
					addLikeSurfaceFaces(likeSurfaceFaces, this, face2)
				}
				return
			}
			var intersectionLine = L3.fromPlanes(thisPlane, face2Plane)
			var thisDir = true
			// get intersections of newCurve with other edges of face and face2
			var ps1 = getFacePlaneIntersectionSs(thisBrep, this, intersectionLine, face2Plane)
			var ps2 = getFacePlaneIntersectionSs(face2Brep, face2, intersectionLine, thisPlane)
			console.log('ps1\n', ps1.map(m => m.toSource()).join('\n'), '\nps2\n', ps2.map(m => m.toSource()).join('\n'))
			if (ps1.length == 0 || ps2.length == 0) {
				// faces to not intersect
				return
			}
			console.log(''+thisPlane+face2Plane)

			/**
			 *
			 * @param seg generated segment
			 * @param col1 if seg is colinear to an edge of this, the edge in question
			 * @param col2 same for face2
			 * @param in1
			 * @param in2
			 * @param a
			 * @param b
			 */
			function handleGeneratedSegment(seg, col1, col2, in1, in2, a, b) {
				console.log("handling seg", col1 && col1.toSource(), col2 && col2.toSource(), seg.toSource())
				if (!col1 && !col2) {
					console.log("adding")
					NLA.mapAdd(faceMap, face, seg)
					NLA.mapAdd(faceMap, face2, seg.flipped())
				}
				if (!col1 && col2) {
					console.log('!col1 && col2')
					if (col2.aDir.cross(face2Plane.normal).dot(thisPlane.normal) > 0) {
						// NB: a new edge is inserted even though it may be the same as an old one
						// however it indicates that it intersects the other volume here, i.e. the old edge cannot
						// be counted as "inside" for purposes of reconstitution
						NLA.mapAdd(faceMap, face2, seg.flipped())

						NLA.mapAdd(edgeMap, col2, {edgeT: col2.curve.pointLambda(seg.a), p: seg.a})
						NLA.mapAdd(edgeMap, col2, {edgeT: col2.curve.pointLambda(seg.b), p: seg.b})
						let testVector = face2Plane.normal.negated().rejectedFrom(thisPlane.normal)
						console.log("testVector", testVector.sce, face2Brep.sce, splitsVolumeEnclosingFaces(face2Brep, col2, testVector, thisPlane.normal) == INSIDE)
						if (splitsVolumeEnclosingFaces(face2Brep, col2, testVector, thisPlane.normal) == INSIDE) {
							NLA.mapAdd(faceMap, face, seg)
						}
					}
				}
				if (!col2 && col1) {
					if (col1.aDir.cross(thisPlane.normal).dot(face2Plane.normal) > 0) {
						NLA.mapAdd(faceMap, face, seg)

						NLA.mapAdd(edgeMap, col1, {edgeT: col1.curve.pointLambda(seg.a), p: seg.a})
						NLA.mapAdd(edgeMap, col1, {edgeT: col1.curve.pointLambda(seg.b), p: seg.b})
						let testVector = thisPlane.normal.negated().rejectedFrom(face2Plane.normal)
						if (splitsVolumeEnclosingFaces(thisBrep, col1, testVector, face2Plane.normal) == INSIDE) {
							NLA.mapAdd(faceMap, face2, seg.flipped())
						}
					}
				}
				if (col1 && col2) {
					// faces are touching edge-edge
					let testVector1 = col1.aDir.cross(thisPlane.normal).negated()
					let testVector2 = col2.aDir.cross(face2Plane.normal).negated()
					console.log("testVector", testVector1.sce, face2Brep.sce, splitsVolumeEnclosingFaces(face2Brep, col2, testVector1, thisPlane.normal) == INSIDE)
					if (splitsVolumeEnclosingFaces(face2Brep, col2, testVector1, thisPlane.normal) == INSIDE
						&& splitsVolumeEnclosingFaces(thisBrep, col1, testVector2, face2Plane.normal) == INSIDE) {
						//NLA.mapAdd(faceMap, face2, seg.flipped())
						NLA.mapAdd(edgeMap, col1, {edgeT: col1.curve.pointLambda(seg.a), p: seg.a})
						NLA.mapAdd(edgeMap, col1, {edgeT: col1.curve.pointLambda(seg.b), p: seg.b})

						NLA.mapAdd(edgeMap, col2, {edgeT: col2.curve.pointLambda(seg.a), p: seg.a})
						NLA.mapAdd(edgeMap, col2, {edgeT: col2.curve.pointLambda(seg.b), p: seg.b})
					}
				}
			}

			/**
			 * What does this do??
			 * @param a
			 * @param b
			 * @param useA
			 * @param useB
			 */
			function handlePoint(a, b, useA, useB) {
				console.log("logging point", useA, a.toSource())
				console.log("logging point", useB, b.toSource())
				if (useA && !useB) {
					if (!a.colinear) {
						NLA.mapAdd(edgeMap, a.edge, a)
					}
					// else colinear segment ends in middle of other face, do nothing
				}
				if (useB && !useA) {
					if (!b.colinear) {
						NLA.mapAdd(edgeMap, b.edge, b)
					}
					// else colinear segment ends in middle of other face, do nothing
				}
				if (useA && useB) {
					// segment start/ends on edge/edge intersection
					if (!a.colinear && !NLA.equals(a.edgeT, a.edge.aT) && !NLA.equals(a.edgeT, a.edge.bT)) {
//						assert(a.edgeT != a.edge.aT && a.edgeT != a.edge.bT, "implement point intersections")
						// ends on a, on colinear segment b bT != a.edge.bT &&
						let testVector = a.edge.aDir.rejectedFrom(b.edge.aDir)
						let sVEF1 = splitsVolumeEnclosingFaces(face2Brep, b.edge, testVector, thisPlane.normal)
						let sVEF2 = splitsVolumeEnclosingFaces(face2Brep, b.edge, testVector.negated(), thisPlane.normal)
						if (!(INSIDE == sVEF1 && INSIDE == sVEF2
							|| (OUTSIDE == sVEF1 || COPLANAR_OPPOSITE == sVEF1) && (OUTSIDE == sVEF2 || COPLANAR_OPPOSITE == sVEF2))) {
							NLA.mapAdd(edgeMap, a.edge, a)
						}
					}
					if (!b.colinear && !NLA.equals(b.edgeT, b.edge.aT) && !NLA.equals(b.edgeT, b.edge.bT)) {
//						assert(b.edgeT != b.edge.aT && b.edgeT != b.edge.bT, "implement point intersections")
						// ends on b, on colinear segment a
						let testVector = b.edge.aDir.rejectedFrom(a.edge.aDir)
						let sVEF1 = splitsVolumeEnclosingFaces(thisBrep, a.edge, testVector, face2Plane.normal)
						let sVEF2 = splitsVolumeEnclosingFaces(thisBrep, a.edge, testVector.negated(), face2Plane.normal)
						if (!(INSIDE == sVEF1 && INSIDE == sVEF2
							|| (OUTSIDE == sVEF1 || COPLANAR_OPPOSITE == sVEF1) && (OUTSIDE == sVEF2 || COPLANAR_OPPOSITE == sVEF2))) {
							NLA.mapAdd(edgeMap, b.edge, b)
						}
					}
				}
			}

			var in1 = false, in2 = false, colinear1 = false, colinear2 = false
			var i = 0, j = 0, last, segments = []
			var startP, startDir, startT
			while (i < ps1.length || j < ps2.length) {
				var a = ps1[i], b = ps2[j]
				if (j >= ps2.length || i < ps1.length && NLA.lt(a.t, b.t)) {
					last = a
					in1 = !in1
					// ": colinear1" to remember the colinear segment, as segments are only handled once it ends (in1 false)
					colinear1 = in1 ? a.colinear && a.edge : colinear1
					i++
				} else if (i >= ps1.length || NLA.gt(a.t, b.t)) {
					last = b
					in2 = !in2
					colinear2 = in2 ? b.colinear && b.edge : colinear2
					j++
				} else {
					last = a
					in1 = !in1
					in2 = !in2
					//if (in1 == in2) {
					a.used = true
					b.used = true
					colinear1 = in1 ? a.colinear && a.edge : colinear1
					colinear2 = in2 ? b.colinear && b.edge : colinear2
					//}
					i++
					j++
				}
				console.log("as", a, b, in1, in2)
				if (startP && !(in1 && in2)) {
					// segment end
					handleGeneratedSegment(new StraightEdge(intersectionLine, startP, last.p, startT, last.t, null, startDir, last.insideDir && last.insideDir.negated()), colinear1, colinear2, in1, in2, a, b)
					startP = undefined
					last.used = true
					handlePoint(a, b, colinear1 || !!a.used, colinear2 || !!b.used)
				} else if (in1 && in2) {
					// new segment just started
					startP = last.p
					startDir = last.insideDir
					startT = last.t
					last.used = true
					handlePoint(a, b, colinear1 || !!a.used, colinear2 || !!b.used)
				}
			}
		}
	}
)

function addLikeSurfaceFaces(likeSurfaceFaces, face1, face2) {
	// There cannot be two subgroups which will later be connected, as the "graph" of like surface faces is fully connected
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

B2.PlaneFace.forVertices = function (planeSurface, vs, ...holeVs) {
	if (planeSurface instanceof P3) {
		planeSurface = new PlaneSurface(planeSurface)
	}
	assert(isCCW(vs, planeSurface.plane.normal), 'isCCW(vs, planeSurface.plane.normal)')
	var edges = vs.map((a, i, vs) => {
		var b = vs[(i + 1) % vs.length]
		return StraightEdge.throughPoints(a, b)
	})
	holeVs.forEach(vs => assert(doubleSignedArea(vs, planeSurface.plane.normal) >= 0, 'doubleSignedArea(vs, planeSurface.plane.normal) >= 0'))
	var holes = holeVs.map(hvs => hvs.map((a, i, vs) => {
		var b = vs[(i + 1) % vs.length]
		return StraightEdge.throughPoints(a, b)
	}))
	return new B2.PlaneFace(planeSurface, edges, holes)
}
function facesWithEdge(edge, faces) {
	return arrayFilterMap(faces, (face) => {
		var matchingEdge = face.getAllEdges().find(e => e.isCoEdge(edge))
		if (matchingEdge) {
			return {face: face, reversed: !edge.a.like(matchingEdge.a), angle: NaN, normalAtEdgeA: null, edge: matchingEdge}
		}
	})
}
/**
 *
 * @param brep
 * @param brepFace
 * @param line
 * @param plane2
 * @returns {Array}
 * p: intersection point
 * insideDir: currently not needed
 * t: param on intersection line
 * edge: face edge doing the intersection
 * edgeT: !!
 * colinear: whether edge is colinear to intersection line
 */
function getFacePlaneIntersectionSs(brep, brepFace, line, plane2) {
	var facePlane = brepFace.surface.plane
	var ps = []
	var loops = brepFace.holes.concat([brepFace.edges])
	loops.forEach(loop => {
		var colinearSegments = loop.map((edge) => edge.colinearToLine(line))
		var intersectionLinePerpendicular = line.dir1.cross(facePlane.normal)

		loop.forEach((edge, i, edges) => {
			var j = (i + 1) % edges.length, nextEdge = edges[j]
			//console.log(edge.toSource()) {p:V3(2, -2.102, 0),
			if (colinearSegments[i]) {
				// edge colinear to intersection line
				var h = (i + edges.length - 1) % edges.length, prevEdge = edges[h]
				var colinearOutside = edge.aDir.cross(facePlane.normal)
				if (prevEdge.bDir.dot(colinearOutside) < 0) {
					ps.push({p: prevEdge.b, insideDir: edge.aDir.negated(), t: line.pointLambda(prevEdge.b), edge: prevEdge, edgeT: prevEdge.bT,
						colinear: false})
				}
				ps.push(
					{p: edge.a, insideDir: edge.aDir, t: line.pointLambda(edge.a), edge: edge, edgeT: edge.aT,
						colinear: true},
					{p: edge.b, insideDir: edge.bDir.negated(), t: line.pointLambda(edge.b), edge: edge, edgeT: edge.bT,
						colinear: true})
				if (nextEdge.aDir.dot(colinearOutside) > 0) {
					ps.push({p: edge.b, insideDir: edge.bDir, t: line.pointLambda(edge.b), edge: edge, edgeT: prevEdge.bT,
						colinear: false})
				}
			} else {
				// not necessarily a straight edge, so multiple intersections are possible
				var edgeTs = edge.edgeISTsWithPlane(plane2)
				for (var k = 0; k < edgeTs.length; k++) {
					var edgeT = edgeTs[k]
					if (edgeT == edge.bT) {
						// endpoint lies on intersection line
						console.log('endpoint lies on intersection line',
							intersectionLinePerpendicular.dot(edge.bDir) , intersectionLinePerpendicular.dot(nextEdge.aDir),
							intersectionLinePerpendicular.dot(edge.bDir) * intersectionLinePerpendicular.dot(nextEdge.aDir), intersectionLinePerpendicular.sce,
							edge.bDir.sce, nextEdge.aDir.sce)
						if (!colinearSegments[j] && intersectionLinePerpendicular.dot(edge.bDir) * intersectionLinePerpendicular.dot(nextEdge.aDir) > 0) {
							// next segment is not colinear and ends on different side
							console.log("adding")
							ps.push({p: edge.b, insideDir: plane2.normal.negated(), t: line.pointLambda(edge.b), edge: edge, edgeT: edge.bT, colinear: false})
							//console.log('end on line, next other side')
						}
					} else if (edgeT != edge.aT) {
						// edge crosses intersection line, neither starts nor ends on it
						var p = edge.curve.at(edgeT)
						assert(line.containsPoint(p), edge.toString() + p+edgeT)
						var insideDir = plane2.normal.negated()
						ps.push({p: p, insideDir: insideDir, t: line.pointLambda(p), edge: edge, edgeT: edgeT, colinear: false})
						console.log('middle')
					}
				}
			}
		})
	})
	// duplicate 't's are ok, as sometimes a segment needs to stop and start again
	// should be sorted so that back facing ones are first
	ps.sort((a, b) => a.t - b.t || a.insideDir.dot(line.dir1))
	return ps
}
function getFacePlaneIntersectionSs2(brep, brepFace, isCurves, surface2, removeCoplanarSame, removeCoplanarOpposite) {
	var faceSurface = brepFace.surface
	var colinearSegments = brepFace.edges.map((edge) => false)
	/*
	 var colinearSegments = brepFace.edges.map((edge) => edge.curve.colinearTo(isCurve))
	var colinearSegmentsInsideCaseTrue = [], colinearSegmentsInsideCaseFalse = []
	for (var i = 0; i < brepFace.edges.length; i++) {
		if (colinearSegments[i]) {
			var edge = brepFace.edges[i]
			var surface2Normal = surface2.normalAt(edge.a)
			var testVector = faceSurface.normalAt(edge.a).rejectedFrom(surface2Normal)
			var csi1 = splitsVolumeEnclosingFaces(brep, edge, testVector, surface2Normal, removeCoplanarSame, removeCoplanarOpposite)
			var csi2 = splitsVolumeEnclosingFaces(brep, edge, testVector.negated(), surface2Normal, removeCoplanarSame, removeCoplanarOpposite)
			var a = INSIDE == csi1, b = INSIDE  == csi1 || COPLANAR_SAME == csi1
			var c = INSIDE == csi2, d = INSIDE  == csi2 || COPLANAR_SAME == csi2
			colinearSegmentsInsideCaseTrue[i] = b != d
			colinearSegmentsInsideCaseFalse[i] = a != c
		}
	}
	*/
	var pss = NLA.arrayFromFunction(isCurves.length, i => [])


	//console.log(colinearSegments, colinearSegmentsInside)
	brepFace.edges.forEach((edge, i, edges) => {
		let j = (i + 1) % edges.length, nextEdge = edges[j]
		//console.log(edge.toSource()) {p:V3(2, -2.102, 0),
		if (colinearSegments[i]) {
			assert(false)
			// edge colinear to intersection
			let outVector = edge.bDir.cross(faceSurface.normal)
			let insideNext = outVector.dot(nextEdge.aDir) > 0
			let caseA = insideNext != colinearSegmentsInsideCaseTrue[i],
				caseB = insideNext != colinearSegmentsInsideCaseFalse[i]
			let colinearSegmentOutsideVector = edge.aDir.cross(faceSurface.normal)
			let displayOnFace = colinearSegmentOutsideVector.dot(surface2.normalAt(edge.b)) > 0
			if (caseA || caseB || displayOnFace != insideNext) {
				ps.push({p: edge.b, insideDir: null, t: isCurve.pointLambda(edge.b), edge: edge, edgeT: edge.bT,
					caseA: caseA, caseB: caseB, colinear: true, hideOnFace: displayOnFace == insideNext})
				//console.log('colinear')
			}
		} else {
			let edgeTs = edge.isTsWithSurface(surface2)
			for (let k = 0; k < edgeTs.length; k++) {
				let edgeT = edgeTs[k]
				if (edgeT == edge.bT) {
					assert(false)
					let intersectionLinePerpendicular = curve.tangentAt(edge.b).cross(faceSurface.normalAt(edge.b))
					// endpoint lies on intersection isCurve
					console.log('endpoint lies on intersection isCurve',
						intersectionLinePerpendicular.dot(edge.bDir) , intersectionLinePerpendicular.dot(nextEdge.aDir),
						intersectionLinePerpendicular.dot(edge.bDir) * intersectionLinePerpendicular.dot(nextEdge.aDir), intersectionLinePerpendicular.sce,
						edge.bDir.sce, nextEdge.aDir.sce)
					if (colinearSegments[j]) {
						// next segment is colinear
						// we need to calculate if the section of the plane intersection isCurve BEFORE the colinear segment is
						// inside or outside the face. It is inside when the colinear segment out vector and the current segment vector
						// point in the same direction (dot > 0)
						// TODO: UUUH?
						let colinearSegmentOutsideVector = nextEdge.aDir.cross(faceSurface.normal)
						let insideFaceBeforeColinear = colinearSegmentOutsideVector.dot(edge.bDir) < 0
						let caseA = insideFaceBeforeColinear != colinearSegmentsInsideCaseTrue[j],
							caseB = insideFaceBeforeColinear != colinearSegmentsInsideCaseFalse[j]
						let displayOnFace = colinearSegmentOutsideVector.dot(surface2.normalAt(edge.a)) > 0
						// if the "inside-ness" changes, add intersection point
						//console.log("segment end on isCurve followed by colinear", insideFaceBeforeColinear != colinearSegmentInsideFace, nextSegmentOutsideVector)
						if (caseA || caseB || displayOnFace != insideFaceBeforeColinear) {
							ps.push({p: edge.b, insideDir: null, t: isCurve.pointLambda(edge.b), edge: edge, edgeT: edge.bT
								, caseA: caseA, caseB: caseB, colinear: true, hideOnFace: displayOnFace == insideFaceBeforeColinear})
							//console.log('next colinear')
						}
					} else if (intersectionLinePerpendicular.dot(edge.bDir) * intersectionLinePerpendicular.dot(nextEdge.aDir) > 0) {
						// next segment is not colinear and ends on different side
						ps.push({p: edge.b, insideDir: null, t: isCurve.pointLambda(edge.b), edge: edge, edgeT: edge.bT, caseA: true, caseB: true})
						//console.log('end on isCurve, next other side')
					}
				} else if (edgeT != edge.aT) {
					// edge crosses is isCurve, neither starts nor ends on it
					// TODO: figure out which curve it is on

					let onCurve = isCurves.length, isCurve
					let p = edge.curve.at(edgeT)
					while (--onCurve >= 0 && !(isCurve = isCurves[onCurve]).containsPoint(p)) {}
					console.log('edgeT', edgeT, 'p', p.sce, edge, onCurve)
					if (onCurve < 0) {
						assert (false, p.sce)
					}
					assert(isCurve.containsPoint(p))
					let ov = edge.tangentAt(edgeT).cross(faceSurface.normalAt(p))
					let newCurveT = isCurve.pointLambda(p, ov)
					let ct = isCurve.tangentAt(newCurveT)
					if (ov.dot(ct) > 0) ct = ct.negated()
					pss[onCurve].push({p: p, insideDir: ct, t: newCurveT, edge: edge, edgeT: edgeT, caseA: true, caseB: true})
					console.log('middle')
				}
			}
		}
	})
	pss.forEach(ps => ps.sort((a, b) => a.t - b.t || assert(false, a.t + ' '+b.t+' '+a.p.sce)))
	return pss
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


var INSIDE = 0, OUTSIDE = 1, COPLANAR_SAME = 2, COPLANAR_OPPOSITE= 3, ALONG_EDGE_OR_PLANE = 4
/**
 *
 * @param {B2} brep BREP to check
 * @param {B2.Edge} edge edge to check
 * @param {V3} dirAtEdgeA the direction vector to check
 * @param {V3} faceNormal If dirAtEdgeA doesn't split a volume, but is along a face, the returned value depends on
 *     wether that faces normal points in the same direction as faceNormal
 * @returns {number} INSIDE, OUTSIDE, COPLANAR_SAME or COPLANAR_OPPOSITE
 */
function splitsVolumeEnclosingFaces(brep, edge, dirAtEdgeA, faceNormal) {
	assert(arguments.length == 4)
	//assert(p.equals(edge.a))
	var ab1 = edge.aDir.normalized()
	var relFaces = facesWithEdge(edge, brep.faces)
	relFaces.forEach(faceInfo => {
		faceInfo.normalAtEdgeA = faceInfo.face.surface.normalAt(edge.a)
		faceInfo.edgeDirAtEdgeA = !faceInfo.reversed
				? faceInfo.edge.aDir
				: faceInfo.edge.bDir
		faceInfo.outsideVector = faceInfo.edgeDirAtEdgeA.cross(faceInfo.normalAtEdgeA)
		faceInfo.angle = (dirAtEdgeA.angleRelativeNormal(faceInfo.outsideVector.negated(), ab1) + 2 * Math.PI + NLA.PRECISION / 2) % (2 * Math.PI)
	})
	assert(relFaces.length != 0, edge.toSource())
	relFaces.sort((a, b) => a.angle - b.angle)
	// assert(relFaces.length % 2 == 0, edge.toSource()) // even number of touching faces

	if (NLA.isZero(relFaces[0].angle)) {
		var coplanarSame = relFaces[0].normalAtEdgeA.dot(faceNormal) > 0
		return coplanarSame ? COPLANAR_SAME : COPLANAR_OPPOSITE
	} else {
		return !relFaces[0].reversed ? INSIDE : OUTSIDE
	}
}
function splitsVolumeEnclosingCone(brep, point, dir) {
	var testPlane = P3.forAnchorAndPlaneVectors(point, dir, dir.getPerpendicular())
	var rays = []
	for (var k = 0; k < brep.faces.length; k++) {
		var face = brep.faces[k]
		if (face.getAllEdges().some(edge => edge.a.like(point))) {
			if (testPlane.isParallelToPlane(face.surface.plane)) {
				return ALONG_EDGE_OR_PLANE
			}
			var intersectionLine = L3.fromPlanes(testPlane, face.surface.plane)
			var ps = getFacePlaneIntersectionSs(null, face, intersectionLine, testPlane)
			var i = 0
			while (i < ps.length) {
				var a = ps[i++], b = ps[i++]
				var out = a.p.like(point)
				if (out || b.p.like(point)) {
					var dir2 = out ? intersectionLine.dir1 : intersectionLine.dir1.negated()
					var angle = (dir.angleRelativeNormal(dir2, testPlane.normal) + 2 * Math.PI + NLA.PRECISION / 2) % (2 * Math.PI)
					rays.push({angle: angle, out: out})
				}
			}
		}
	}
	rays.sort((a, b) => a.angle - b.angle)
	//console.log("testPlane", testPlane.toSource(), "rays", rays.toSource())

	if (NLA.isZero(rays[0].angle)) {
		return ALONG_EDGE_OR_PLANE
	} else {
		return rays[0].out ? OUTSIDE : INSIDE
	}
}
B2.RotationFace = function (rot, contour, holes, name) {
	B2.Face.call(this, rot, contour, holes, name)
	//assertInst(RotationReqFofZ, rot)
}
B2.RotationFace.prototype = NLA.defineObject(B2.Face.prototype, {
	constructor: B2.RotationFace,
	name: 'B2.RotationFace',
	unrollLoop: function (edgeLoop) {
		var vs = []
		edgeLoop.forEach((edge, e) => {
			var reverseFunc = this.surface.pointToParameterFunction()
			var hint = edge.bDir
			if (edge instanceof StraightEdge && edge.curve.dir1.isParallelTo(this.surface.dir || this.surface.dir1)) {
				hint = this.surface.normalAt(edge.b).cross(edge.bDir)
			}
			edge.getVerticesNo0().forEach(p => vs.push(reverseFunc(p, hint)))
		})
		//console.log("e2\n", vs.join("\n"), vs.length)
		return vs
	},
	addToMesh: function (mesh) {
		var closed = false
		var hSplit = 32, zSplit = 1
		var ribs = []
		var minZ = Infinity, maxZ = -Infinity
		var cmp = (a, b) => a.value - b.value
		var f = this.surface.parametricFunction()
		var normalF = this.surface.parametricNormal()
		var vertexLoops = this.holes.concat([this.edges]).map(loop => this.unrollLoop(loop))
		vertexLoops.forEach(vertexLoop => {
			vertexLoop.forEach(({x: d, y: z}) => {
				var index0 = ribs.binaryIndexOf(d, (a, b) => NLA.snapTo(a.value - b, 0))
				if (index0 < 0) {
					ribs.splice(-index0-1, 0, {value: d, left: [], right: []})
				}
				minZ = min(minZ, z)
				maxZ = max(maxZ, z)
			})
		})
		var correction = 1
		vertexLoops.forEach(vertexLoop => {
			vertexLoop.forEach((v0, i, vs) => {
				var v1 = vs[(i + 1) % vs.length], dDiff = v1.x - v0.x
				//console.log(v0.sce, v1.sce)
				if (NLA.isZero(dDiff)) {
					return
				}
				if (dDiff < 0) {
					[v0, v1] = [v1, v0]
					dDiff = -dDiff
				}
				var index0 = ribs.binaryIndexOf(v0.x, (a, b) => NLA.snapTo(a.value - b, 0))
				var index1 = ribs.binaryIndexOf(v1.x, (a, b) => NLA.snapTo(a.value - b, 0))
				ribs[index0].right.binaryInsert(v0.y)
				for (var j = (index0 + correction) % ribs.length; j != index1; j = (j + correction) % ribs.length) {
					var x = ribs[j].value
					var part = (x - v0.x) / dDiff
					var interpolated = v1.y * part + v0.y * (1 - part)
					ribs[j].left.binaryInsert(interpolated)
					ribs[j].right.binaryInsert(interpolated)
				}
				ribs[index1].left.binaryInsert(v1.y)
				//console.log(ribs.map(r=>r.toSource()).join('\n'))
			})
		})
		var vertices = [], triangles = [], normals = []
		for (let i = 0; i < ribs.length; i++) {
			let ribLeft = ribs[i], ribRight = ribs[(i + 1) % ribs.length]
			assert(ribLeft.right.length == ribRight.left.length)
			for (let j = 0; j < ribLeft.right.length; j++) {
				vertices.push(f(ribLeft.value, ribLeft.right[j]), f(ribRight.value, ribRight.left[j]))
				normals.push(normalF(ribLeft.value, ribLeft.right[j]), normalF(ribRight.value, ribRight.left[j]))
			}
		}
		//console.log(ribs.map(r=>r.toSource()).join('\n'))
		var vss = vertices.length, detailVerticesStart = vss
		var zInterval = maxZ - minZ, zStep = zInterval / zSplit
		var detailZs = NLA.arrayFromFunction(zSplit - 1, i => minZ + (1 + i) * zStep)
		for (let i = 0; i < ribs.length; i++) {
			var d = ribs[i].value
			for (let j = 0; j < detailZs.length; j++) {
				vertices.push(f(d, detailZs[j]))
				normals.push(normalF(d, detailZs[j]))
			}
		}
		//console.log('detailVerticesStart', detailVerticesStart, 'vl', vertices.length, vertices.length - detailVerticesStart, ribs.length)
		// finally, fill in the ribs
		var vsStart = 0
		var flipped2 = this.surface instanceof ProjectedCurveSurface ? true : this.surface.normalDir == 1
		//for (var i = 0; i < 1; i++) {
		var end = closed ? ribs.length : ribs.length - 1
		for (let i = 0; i < end; i++) {
			var ipp = (i + 1) % ribs.length
			let inside = false, colPos = 0, ribLeft = ribs[i], ribRight = ribs[(i + 1) % ribs.length]
			for (let j = 0; j < detailZs.length + 1; j++) {
				var detailZ = detailZs[j] || 100000
				if (!inside) {
					if (ribLeft.right[colPos] < detailZ && ribRight.left[colPos] < detailZ) {
						if (ribLeft.right[colPos + 1] < detailZ || ribRight.left[colPos + 1] < detailZ) {
							pushQuad(triangles, flipped2,
								vsStart + colPos * 2,
								vsStart + (colPos + 1) * 2,
								vsStart + colPos * 2 + 1,
								vsStart + (colPos + 1) * 2 + 1)
							colPos += 2
							if (ribLeft.right[colPos] < detailZ || ribRight.left[colPos] < detailZ) {
								j--
							}
						} else {
							pushQuad(triangles, flipped2,
								vsStart + colPos * 2,
								vsStart + colPos * 2 + 1,
								detailVerticesStart + i * detailZs.length + j,
								detailVerticesStart + ipp * detailZs.length + j)
							inside = true
							colPos++
						}
					}
				} else {
					if (ribLeft.right[colPos] < detailZ || ribRight.left[colPos] < detailZ) {
						pushQuad(triangles, flipped2,
							detailVerticesStart + i * detailZs.length + j - 1,
							detailVerticesStart + ipp * detailZs.length + j - 1,
							vsStart + colPos * 2,
							vsStart + colPos * 2 + 1)
						inside = false
						colPos++
						if (ribLeft.right[colPos] < detailZ || ribRight.left[colPos] < detailZ) {
							j--
						}
					} else {
						pushQuad(triangles, flipped2,
							detailVerticesStart + i * detailZs.length + j,
							detailVerticesStart + i * detailZs.length + j - 1,
							detailVerticesStart + ipp * detailZs.length + j,
							detailVerticesStart + ipp * detailZs.length + j - 1)
					}
				}
			}
			vsStart += ribLeft.right.length * 2
		}
		//console.log('trinagle', triangles.max(), vertices.length, triangles.length, triangles.toSource(), triangles.map(i => vertices[i].$).toSource() )
		triangles = triangles.map(index => index + mesh.vertices.length)
		//assert(normals.every(n => n.hasLength(1)), normals.find(n => !n.hasLength(1)).length() +" "+normals.findIndex(n => !n.hasLength(1)))
		mesh.faceIndexes.set(this, {start: mesh.triangles.length, count: triangles.length})
		Array.prototype.push.apply(mesh.vertices, vertices)
		Array.prototype.push.apply(mesh.triangles, triangles)
		Array.prototype.push.apply(mesh.normals, normals)
		//this.addEdgeLines(mesh)

	},
	doo: function (face) {
		if (face.surface instanceof PlaneSurface) {
			face.doo.apply(face, [this].concat(Array.from(arguments).slice(1)))
		} else {
			B2.PlaneFace.prototype.doo.apply(this, arguments)
		}
	}
})
B2.RotationFace.prototype.constructor = B2.RotationFace


B2.PCurveEdge = NLA.defineClass('B2.PCurveEdge', B2.Edge,
	function (curve, a, b, aT, bT, flippedOf, aDir, bDir) {
		assertNumbers(aT, bT)
		assertVectors(a, b, aDir, bDir)
		assertf(() => curve instanceof L3 || curve instanceof Curve, curve)
		assertf(() => !curve.isValidT || curve.isValidT(aT) && curve.isValidT(bT), aT + ' ' + bT)
		assertf(() => curve.at(aT).like(a), curve.at(aT)+a)
		assertf(() => curve.at(bT).like(b), curve.at(bT)+b)
		assertf(() => curve.tangentAt(aT).likeOrReversed(aDir), curve.tangentAt(aT).sce +' '+ aDir.sce)
		assertf(() => curve.tangentAt(bT).likeOrReversed(bDir))
		this.curve = curve
		this.a = a
		this.b = b
		this.aT = aT
		this.bT = bT
		this.aDir = aDir
		this.bDir = bDir
		this.canon = flippedOf
		this.reversed = this.aDir.dot(curve.tangentAt(aT)) < 0
		assert(this.reversed != aT < bT, aT+' '+bT+' '+curve.constructor.name+' '+this.aDir.sce+' '+this.bDir.sce + ' '+curve.tangentAt(aT))
		this.id = globalId++
	},
	{
		toSource: function() {
			return `new B2.PCurveEdge(${this.curve}, ${this.a}, ${this.b}, ${this.aT}, ${this.bT}, ${this.canon}, ${this.aDir}, ${this.bDir})`
		},
		getVerticesNo0: function () {
			return this.curve.calcSegmentPoints(this.aT, this.bT, this.a, this.b, this.reversed, false)
		},
		get points() {
			return this.curve.calcSegmentPoints(this.aT, this.bT, this.a, this.b, this.reversed, true)
		},
		rotViaPlane: function (normal, reversed) {
			var rot = this.aDir.angleRelativeNormal(this.bDir, normal)
			var counterClockWise = (normal.dot(this.curve.normal) > 0) == !this.reversed
			if (counterClockWise) {
				// counterclockwise rotation, i.e. rot > 0
				if (rot < 0) rot += 2 * PI
			} else {
				if (rot > 0) rot -= 2 * PI
			}
			return rot
		},
		edgeISTsWithSurface: function (surface) {
			return this.curve.isTsWithSurface(surface).filter(edgeT => {
				var aT = this.aT, bT = this.bT
				edgeT = NLA.snapTo(edgeT, aT)
				edgeT = NLA.snapTo(edgeT, bT)
				if (!this.reversed) {
					if (aT < bT) {
						return aT <= edgeT && edgeT <= bT
					} else {
						return !(bT < edgeT && edgeT < aT)
					}
				} else {
					if (aT > bT) {
						return aT >= edgeT && edgeT >= bT
					} else {
						return !(bT > edgeT && edgeT > aT)
					}
				}
			})
		},
		edgeISTsWithPlane: function (surface) {
			return this.curve.isTsWithPlane(surface).filter(edgeT => {
				var aT = this.aT, bT = this.bT
				edgeT = NLA.snapTo(edgeT, aT)
				edgeT = NLA.snapTo(edgeT, bT)
				if (!this.reversed) {
					if (aT < bT) {
						return aT <= edgeT && edgeT <= bT
					} else {
						return !(bT < edgeT && edgeT < aT)
					}
				} else {
					if (aT > bT) {
						return aT >= edgeT && edgeT >= bT
					} else {
						return !(bT > edgeT && edgeT > aT)
					}
				}
			})
		},
		tangentAt: function (t) {
			return !this.reversed ? this.curve.tangentAt(t) : this.curve.tangentAt(t).negated()
		},
		flipped: function () {
			return new B2.PCurveEdge(this.curve, this.b, this.a, this.bT, this.aT, this, this.bDir.negated(), this.aDir.negated())
		},
		transform: function (m4) {
			return new B2.PCurveEdge(this.curve.transform(m4), m4.transformPoint(this.a), m4.transformPoint(this.b),
				this.aT, this.bT,
				null,
				m4.transformVector(this.aDir), m4.transformVector(this.bDir))
		},
		isCoEdge: function (edge) {
			// TODO: optimization with flippedOf etc
			return this == edge ||
				this.curve.isColinearTo(edge.curve) && (
					this.a.like(edge.a) && this.b.like(edge.b)
					|| this.a.like(edge.b) && this.b.like(edge.a)
				)
		},
	}
)
B2.PCurveEdge.forCurveAndTs = function (curve, aT, bT) {
	return new B2.PCurveEdge(curve, curve.at(aT), curve.at(bT), aT, bT, undefined,
		aT < bT ? curve.tangentAt(aT) : curve.tangentAt(aT).negated(),
		aT < bT ? curve.tangentAt(bT) : curve.tangentAt(bT).negated())
}
var StraightEdge = NLA.defineClass('StraightEdge', B2.Edge,
	function (line, a, b, aT, bT, flippedOf) {
		assertNumbers(aT, bT)
		assertVectors(a, b)
		assertInst(L3, line)
		assert(line.containsPoint(a), 'line.containsPoint(a)'+line+a)
		assert(line.containsPoint(b), 'line.containsPoint(b)'+line+b)
		this.curve = line
		this.a = a || line.at(aT)
		this.b = b || line.at(bT)
		this.aT = aT
		this.bT = bT
		this.reversed = this.aT > this.bT
		this.canon = flippedOf
		this.tangent = this.aT < this.bT ? this.curve.dir1 : this.curve.dir1.negated()
		this.id = globalId++
	},
	{

		toSource: function () {
			return `StraightEdge.throughPoints(${this.a}, ${this.b})`
		},
		getVerticesNo0: function () {
			return [this.b]
		},
		get points() {
			return [this.a, this.b]
		},
		edgeISTsWithPlane: function (plane) {
			var minT = min(this.aT, this.bT), maxT = max(this.aT, this.bT)
			var edgeT = this.curve.intersectWithPlaneLambda(plane)
			edgeT = NLA.snapTo(edgeT, this.aT)
			edgeT = NLA.snapTo(edgeT, this.bT)
			return (minT <= edgeT && edgeT <= maxT) ? [edgeT] : []
		},
		edgeISTsWithSurface: function (surface) {
			if (surface instanceof PlaneSurface) {
				return this.edgeISTsWithPlane(surface.plane)
			} else if (surface instanceof CylinderSurface) {
				var minT = min(this.aT, this.bT), maxT = max(this.aT, this.bT)
				return surface.isPointsWithLine(this.curve)
					.map(p => this.curve.pointLambda(p))
					.filter(edgeT => minT <= edgeT && edgeT <= maxT)
			} else {
				assert(false)
			}
		},
		tangentAt: function (p) {
			return this.tangent
		},
		flipped: function () {
			return new StraightEdge(this.curve, this.b, this.a, this.bT, this.aT, this)
		},
		get aDir() { return this.tangent },
		get bDir() { return this.tangent },
		set aDir(x) {  },
		set bDir(x) {  },
		transform: function (m4) {
			return new StraightEdge(this.curve.transform(m4), m4.transformPoint(this.a), m4.transformPoint(this.b), this.aT, this.bT)
		},
		isCoEdge: function (edge) {
			// TODO: optimization with flippedOf etc
			return edge.constructor == StraightEdge && (
					this.a.like(edge.a) && this.b.like(edge.b)
					|| this.a.like(edge.b) && this.b.like(edge.a)
				)
		},
		likeEdge: function (edge) {
			return edge.constructor == StraightEdge && this.a.like(edge.a) && this.b.like(edge.b)
		},
		equals: function (edge) {
			return edge.constructor == StraightEdge && this.a.equals(edge.a) && this.b.equals(edge.b)
		},
		getEdgeT: function (p) {
			assertVectors(p)
			var edgeT = p.minus(this.curve.anchor).dot(this.curve.dir1)
			if (!NLA.isZero(this.curve.at(edgeT).distanceTo(p))) { return }
			var minT = min(this.aT, this.bT), maxT = max(this.aT, this.bT)
			edgeT = NLA.snapTo(edgeT, this.aT)
			edgeT = NLA.snapTo(edgeT, this.bT)
			return (minT <= edgeT && edgeT <= maxT) ? edgeT : undefined
		}
	}
)
StraightEdge.throughPoints = function (a, b) {
	return new StraightEdge(L3.throughPoints(a, b), a, b, 0, b.minus(a).length())
}

/**
 * Solves a quadratic system of equations of the form
 *      a * x + b * y = c
 *      x^2 + y^2 = 1
 * This can be understood as the intersection of the unit circle with a line.
 *
 * @param {number} a double
 * @param {number} b double
 * @param {number} c double
 * @returns {{x1: number, y1: number, x2: number, y2: number}} with x1 >= x2 and y1 <= y2
 * a * b + (b -c) * (b + c)
 */
function intersectionUnitCircleLine(a, b, c) {
	assertNumbers(a, b, c)
	// TODO: disambiguate on a < b
	var term = sqrt(a * a + b * b - c * c)
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
 * @param {number} a
 * @param {number} b
 * @param {number} c
 * @returns {{x1:number, y1:number, x2:number, y2:number}} with x1 >= x2 and y1 <= y2
 * a * b + (b -c) * (b + c)
 */
function intersectionUnitHyperbolaLine(a, b, c) {
	assertNumbers(a, b, c)
	var aa = a*a, bb = b*b, cc = c*c
	// TODO: disambiguate on a < b
	//var xTerm = sqrt(4*cc*aa-4*(bb-aa)*(-cc-bb))
	var xTerm = 2 * sqrt(bb*cc + bb*bb - aa*bb)
	var yTerm = sqrt(4*cc*bb-4*(bb-aa)*(cc-aa))
	return {
		x1: (-2 * a * c + xTerm) / 2 / (bb - aa),
		x2: (-2 * a * c - xTerm) / 2 / (bb - aa),
		y1: (2 * b * c - yTerm) / 2 / (bb - aa),
		y2: (2 * b * c + yTerm) / 2 / (bb - aa)
	}
}
/**
 *
 * @param {number} a
 * @param {number} b
 * @param {number} c
 * @param {number} r
 * @returns {{x1: number, x2: number, y1: number, y2: number}}
 */
function intersectionCircleLine(a, b, c, r) {
	assertNumbers(a, b, c, r)
	var term = sqrt(r * r * (a * a + b * b) - c * c)
	return {
		x1: (a * c + b * term) / (a * a + b * b),
		x2: (a * c - b * term) / (a * a + b * b),
		y1: (b * c - a * term) / (a * a + b * b),
		y2: (b * c + a * term) / (a * a + b * b)
	}
}

/**
 *
 * @param {Curve} curve
 * @param {number} startT
 * @param {number} endT
 * @param {number} steps integer
 * @returns {number}
 */
function integrateCurve(curve, startT, endT, steps) {
	var step = (endT - startT) / steps
	var length = 0
	var p = curve.at(startT)
	for (var i = 0, t = startT + step; i < steps; i++, t += step) {
		var next = curve.at(t)
		length += p.distanceTo(next)
		p = next
	}
	return length
}

var gaussLegendreXs = [
	-0.0640568928626056260850430826247450385909,
	0.0640568928626056260850430826247450385909,
	-0.1911188674736163091586398207570696318404,
	0.1911188674736163091586398207570696318404,
	-0.3150426796961633743867932913198102407864,
	0.3150426796961633743867932913198102407864,
	-0.4337935076260451384870842319133497124524,
	0.4337935076260451384870842319133497124524,
	-0.5454214713888395356583756172183723700107,
	0.5454214713888395356583756172183723700107,
	-0.6480936519369755692524957869107476266696,
	0.6480936519369755692524957869107476266696,
	-0.7401241915785543642438281030999784255232,
	0.7401241915785543642438281030999784255232,
	-0.8200019859739029219539498726697452080761,
	0.8200019859739029219539498726697452080761,
	-0.8864155270044010342131543419821967550873,
	0.8864155270044010342131543419821967550873,
	-0.9382745520027327585236490017087214496548,
	0.9382745520027327585236490017087214496548,
	-0.9747285559713094981983919930081690617411,
	0.9747285559713094981983919930081690617411,
	-0.9951872199970213601799974097007368118745,
	0.9951872199970213601799974097007368118745
]
var gaussLegendreWeights = [
	0.1279381953467521569740561652246953718517,
	0.1279381953467521569740561652246953718517,
	0.1258374563468282961213753825111836887264,
	0.1258374563468282961213753825111836887264,
	0.1216704729278033912044631534762624256070,
	0.1216704729278033912044631534762624256070,
	0.1155056680537256013533444839067835598622,
	0.1155056680537256013533444839067835598622,
	0.1074442701159656347825773424466062227946,
	0.1074442701159656347825773424466062227946,
	0.0976186521041138882698806644642471544279,
	0.0976186521041138882698806644642471544279,
	0.0861901615319532759171852029837426671850,
	0.0861901615319532759171852029837426671850,
	0.0733464814110803057340336152531165181193,
	0.0733464814110803057340336152531165181193,
	0.0592985849154367807463677585001085845412,
	0.0592985849154367807463677585001085845412,
	0.0442774388174198061686027482113382288593,
	0.0442774388174198061686027482113382288593,
	0.0285313886289336631813078159518782864491,
	0.0285313886289336631813078159518782864491,
	0.0123412297999871995468056670700372915759,
	0.0123412297999871995468056670700372915759
]

function gaussLegendreQuadrature24(fn, startT, endT) {
	let result = 0
	for (let i = 0; i < gaussLegendreXs.length; i++) {
		// gauss-legendre goes from -1 to 1, so we need to scale
		let t = startT + (gaussLegendreXs[i] + 1) / 2 * (endT - startT)
		result += gaussLegendreWeights[i] * fn(t)
	}
	// again, [-1,1], so div by 2
	return result / 2 * (endT - startT)
}

function curveLengthByDerivative(df, startT, endT, steps) {
	var dt = (endT - startT) / steps
	var length = 0
	for (var i = 0, t = startT + dt / 2; i < steps; i++, t += dt) {
		length += df(t) * dt
	}
	return length
}

function getIntersectionSegments (ps1, ps2, in1, in2, constructor, curve) {
	var currentSegment
	assert (!(in1 && in2), '!(in1 && in2) '+in1+' '+in2)
	console.log('in', in1, in2)
	// generate overlapping segments
	var i = 0, j = 0, last, segments = []
	var startP, startDir, startT
	while (i < ps1.length || j < ps2.length) {
		var a = ps1[i], b = ps2[j]
		if (j >= ps2.length || i < ps1.length && NLA.lt(a.t, b.t)) {
			last = a
			in1 = !in1
			i++
		} else if (i >= ps1.length || NLA.gt(a.t, b.t)) {
			last = b
			in2 = !in2
			j++
		} else {
			last = a
			in1 = !in1
			in2 = !in2
			if (in1 == in2) {
				a.used = true
				b.used = true
			}
			i++
			j++
		}
//		console.log("as", a, b, in1, in2)
		if (startP && !(in1 && in2)) {
			segments.push(new constructor(curve, startP, last.p, startT, last.t, null, startDir, last.insideDir.negated()))
			startP = undefined
			last.used = true
		} else if (in1 && in2) {
			startP = last.p
			startDir = last.insideDir
			startT = last.t
			last.used = true
		}
	}
	assert (!(in1 && in2), '!(in1 && in2) '+in1+' '+in2)
	return segments
}
function getBlug(ps1, ps2, curve) {
	// generate overlapping segments
	var in1 = false, in2 = false
	var i = 0, j = 0, last, segments = []
	var startP, startDir, startT
	// TODO : skip -><-
	while (i < ps1.length || j < ps2.length) {
		var a = ps1[i], b = ps2[j]
		if (j >= ps2.length || i < ps1.length && NLA.lt(a.t, b.t)) {
			last = a
			in1 = !in1
			i++
		} else if (i >= ps1.length || NLA.gt(a.t, b.t)) {
			last = b
			in2 = !in2
			j++
		} else {
			last = a
			in1 = !in1
			in2 = !in2
			//if (in1 == in2) {
				a.used = true
				b.used = true
			//}
			i++
			j++
		}
//		console.log("as", a, b, in1, in2)
		if (startP && !(in1 && in2)) {
			segments.push(new StraightEdge(curve, startP, last.p, startT, last.t, null, startDir, last.insideDir && last.insideDir.negated()))
			startP = undefined
			last.used = true
		} else if (in1 && in2) {
			startP = last.p
			startDir = last.insideDir
			startT = last.t
			last.used = true
		}
	}
	assert (!in1 && !in2, '!in1 && !in2 '+in1+' '+in2)
	return segments
}
/**
 *
 * @param {P3} plane
 * @param {V3=} right
 * @param {V3=} up
 * @constructor
 * @extends {Surface}
 */
function PlaneSurface(plane, right, up) {
	assertInst(P3, plane)
	this.plane = plane
	this.up = up || plane.normal.getPerpendicular().normalized()
	this.right = right || this.up.cross(this.plane.normal).normalized()
	assert(this.right.cross(this.up).like(this.plane.normal))
}
PlaneSurface.throughPoints = function (a, b, c) {
	return new PlaneSurface(P3.throughPoints(a, b, c))
}
PlaneSurface.prototype = NLA.defineObject(Surface.prototype, /** @lends PlaneSurface.prototype */ {
	isCoplanarTo: function (surface) {
		return surface instanceof PlaneSurface && this.plane.isCoplanarToPlane(surface.plane)
	},
	parametricFunction: function () {
		var matrix = M4.forSys(this.right, this.up, this.normal, this.plane.anchor)
		return function (s, t) {
			return matrix.transformPoint(V3.create(s, t, 0))
		}
	},
	implicitFunction: function () {
		return p => this.plane.distanceToPointSigned(p)
	},
	isCurvesWithISurface: function (implicitSurface) {
		assert (implicitSurface.implicitFunction, 'implicitSurface.implicitFunction')
		return new CurvePI(this, implicitSurface)
	},

	/**
	 * @inheritDoc
	 */
	isCurvesWithSurface: function (surface2) {
		// prefer other surface to be the paramteric one
		if (surface2.implicitFunction) {
			return new CurvePI(this, surface2)
		} else if (surface2.parametricFunction) {
			return new CurvePI(surface2, this)
		}
	},
	/**
	 *
	 * @param {Array.<B2.Edge>} contour
	 * @returns {boolean}
	 */
	edgeLoopCCW: function (contour) {
		var totalAngle = 0
		for (var i = 0; i < contour.length; i++) {
			var ipp = (i + 1) % contour.length
			var edge = contour[i], nextEdge = contour[ipp]
			assert(edge.b.like(nextEdge.a), "edges dont form a loop")
			if (edge.curve instanceof EllipseCurve) {
				totalAngle += edge.rotViaPlane(this.plane.normal)
				console.log(edge.toString(), edge.rotViaPlane(this.plane.normal))
			}
			totalAngle += edge.bDir.angleRelativeNormal(nextEdge.aDir, this.plane.normal)
		}
		return totalAngle > 0
	},

	/**
	 *
	 * @param {Array.<B2.Edge>} contour
	 * @param {V3} p
	 * @returns {boolean}
	 */
	edgeLoopContainsPoint: function (contour, p) {
		assert(contour)
		assertVectors (p)
		var dir = this.right.plus(this.up.times(0.123)).normalized()
		var line = L3(p, dir)
		var plane = this.plane
		var intersectionLinePerpendicular = dir.cross(plane.normal)
		var plane2 = P3.normalOnAnchor(intersectionLinePerpendicular, p)
		var colinearSegments = contour.map((edge) => edge.colinearToLine(line))
		var colinearSegmentsInside = contour.map((edge, i) => edge.aDir.dot(dir) > 0)
		var inside = false
		function logIS(p) {
			if (line.pointLambda(p) > 0) {
				inside = !inside
			}
		}
		contour.forEach((edge, i, edges) => {
			var j = (i + 1) % edges.length, nextEdge = edges[j]
			//console.log(edge.toSource()) {p:V3(2, -2.102, 0),
			if (colinearSegments[i]) {
				// edge colinear to intersection
				var outVector = edge.bDir.cross(plane.normal)
				var insideNext = outVector.dot(nextEdge.aDir) > 0
				if (colinearSegmentsInside[i] != insideNext) {
					logIS(edge.b)
				}
			} else {
				var edgeTs = edge.edgeISTsWithPlane(plane2)
				for (var k = 0; k < edgeTs.length; k++) {
					var edgeT = edgeTs[k]
					if (edgeT == edge.bT) {
						// endpoint lies on intersection line
						if (colinearSegments[j]) {
							// next segment is colinear
							// we need to calculate if the section of the plane intersection line BEFORE the colinear segment is
							// inside or outside the face. It is inside when the colinear segment out vector and the current segment vector
							// point in the same direction (dot > 0)
							var colinearSegmentOutsideVector = nextEdge.aDir.cross(plane.normal)
							var insideFaceBeforeColinear = colinearSegmentOutsideVector.dot(edge.bDir) < 0
							// if the "inside-ness" changes, add intersection point
							//console.log("segment end on line followed by colinear", insideFaceBeforeColinear != colinearSegmentInsideFace, nextSegmentOutsideVector)
							if (colinearSegmentsInside[j] != insideFaceBeforeColinear) {
								logIS(edge.b)
							}
						} else if (intersectionLinePerpendicular.dot(edge.bDir) * intersectionLinePerpendicular.dot(nextEdge.aDir) > 0) {
							logIS(edge.b)
						}
					} else if (edgeT != edge.aT) {
						// edge crosses line, neither starts nor ends on it
						logIS(edge.curve.at(edgeT))
					}
				}
			}
		})
		return inside

	},
	pointToParameterFunction: function (p) {
		var matrix = M4.forSys(this.right, this.up, this.normal, this.plane.anchor)
		var matrixInverse = matrix.inversed()
		return function (pWC) {
			return matrixInverse.transformPoint(pWC)
		}
	},
	normalAt: function (p) {
		return this.plane.normal
	},
	containsPoint: function (p) { return this.plane.containsPoint(p) },
	containsCurve: function (curve) {
		if (curve instanceof L3) {
			return this.plane.containsLine(curve)
		} else if (curve instanceof EllipseCurve) {
			return this.plane.containsPoint(curve.center) && this.plane.normal.isParallelTo(curve.normal)
		} else if (curve instanceof BezierCurve) {
			return curve.points.every(p => this.plane.containsPoint(p))
		} else {
			assert(false, edge.toString())
		}
	},
	transform: function (m4) {
		return new PlaneSurface(this.plane.transform(m4))
	},
	flipped: function () {
		return new PlaneSurface(this.plane.flipped(), this.right, this.up.negated())
	},
	toString: function () {
		return this.plane.toString()
	},
	toSource: function () {
		return `new PlaneSurface(${this.plane})`
	},
	equalsSurface: function (surface) {
		return surface instanceof PlaneSurface && this.plane.like(surface.plane)
	}
})
function curvePoint(implicitCurve, startPoint) {
	var eps = 1e-5
	var p = startPoint
	for (var i = 0; i < 4; i++) {
		var fp = implicitCurve(p.x, p.y)
		var dfpdx = (implicitCurve(p.x + eps, p.y) - fp) / eps,
			dfpdy = (implicitCurve(p.x, p.y + eps) - fp) / eps
		var scale = fp / (dfpdx * dfpdx + dfpdy * dfpdy)
		//console.log(p.$)
		p = p.minus(V3(scale * dfpdx, scale * dfpdy))
	}
	return p
}
function followAlgorithm (implicitCurve, startPoint, endPoint, stepLength, startDir, tangentEndPoints, boundsFunction) {
	NLA.assertNumbers(stepLength, implicitCurve(0, 0))
	NLA.assertVectors(startPoint, endPoint)
	assert (!startDir || startDir instanceof V3)
	var points = []
	tangentEndPoints = tangentEndPoints || []
	assert (NLA.isZero(implicitCurve(startPoint.x, startPoint.y)), 'NLA.isZero(implicitCurve(startPoint.x, startPoint.y))')
	stepLength = stepLength || 0.5
	var eps = 1e-5
	var p = startPoint, prevp = startDir ? p.minus(startDir) : p
	var i = 0
	do {
		var fp = implicitCurve(p.x, p.y)
		var dfpdx = (implicitCurve(p.x + eps, p.y) - fp) / eps,
			dfpdy = (implicitCurve(p.x, p.y + eps) - fp) / eps
		var tangent = V3.create(-dfpdy, dfpdx, 0)
		var reversedDir = p.minus(prevp).dot(tangent) < 0
		tangent = tangent.toLength(reversedDir ? -stepLength : stepLength)
		var tangentEndPoint = p.plus(tangent)
		points.push(p)
		tangentEndPoints.push(tangentEndPoint)
		prevp = p
		p = curvePoint(implicitCurve, tangentEndPoint)
	} while (i++ < 100 && (i < 4 || prevp.distanceTo(endPoint) > 1.1 * stepLength) && boundsFunction(p.x, p.x))
	// TODO gleichm¨aßige Verteilung der Punkte
	return points
}
// both curves must be in the same s-t coordinates for this to make sense
function intersectionICurveICurve(pCurve1, startParams1, endParams1, startDir, stepLength, pCurve2) {
	NLA.assertNumbers(stepLength, pCurve1(0, 0), pCurve2(0, 0))
	NLA.assertVectors(startParams1, endParams1)
	assert (!startDir || startDir instanceof V3)
	var vertices = []
	assert (NLA.isZero(pCurve1(startParams1.x, startParams1.y)))
	stepLength = stepLength || 0.5
	var eps = 1e-5
	var p = startParams1, prevp = p // startDir ? p.minus(startDir) : p
	var i = 0
	while (i++ < 1000 && (i < 4 || p.distanceTo(endParams1) > 1.1 * stepLength)) {
		var fp = pCurve1(p.x, p.y)
		var dfpdx = (pCurve1(p.x + eps, p.y) - fp) / eps,
			dfpdy = (pCurve1(p.x, p.y + eps) - fp) / eps
		var tangent = V3(-dfpdy, dfpdx, 0).toLength(stepLength)
		if (p.minus(prevp).dot(tangent) < 0) tangent = tangent.negated()
		prevp = p
		p = curvePoint(pCurve1, p.plus(tangent))
		vertices.push(p)
	}
	// TODO gleichmäßige Verteilung der Punkte
	return vertices

}
function asj(iCurve1, loopPoints1, iCurve2) {
	var p = loopPoints1[0], val = iCurve2(p.x, p.y), lastVal
	var iss = []
	for (var i = 0; i < loopPoints1.length; i++) {
		lastVal = val
		p = loopPoints1[i]
		val = iCurve2(p)
		if (val * lastVal <= 0) { // TODO < ?
			iss.push(newtonIterate2d(iCurve1, iCurve2, p.x, p.y))
		}
	}
	return iss
}
/**
 *
 * @param {function (number, number):number} f1
 * @param {function (number, number):number} f2
 * @param {number} startS
 * @param {number} startT
 * @param {number=} iterations
 * @returns {V3}
 */
function newtonIterate2d(f1, f2, startS, startT, iterations) {
	iterations = iterations || 4
	var s = startS, t = startT
	var eps = 1e-6
	do {
		/*
			| a b |-1                   |  d -b |
			| c d |   = 1 / (ad - bc) * | -c  a |
		 */
		var f1ts = f1(s, t), f2ts = f2(s, t)
		/*
		var df1s = (f1(s + eps, t) - f1ts) / eps, df1t = (f1(s, t + eps) - f1ts) / eps,
			df2s = (f2(s + eps, t) - f2ts) / eps, df2t = (f2(s, t + eps) - f2ts) / eps
		var det = df1s * df2t - df1t * df2s
		s = s - ( df2t * f1ts - df1t * f2ts) / det
		t = t - (-df2s * f1ts + df1s * f2ts) / det
		*/
		// TODO: is this even more accurate?
		var df1s = (f1(s + eps, t) - f1ts), df1t = (f1(s, t + eps) - f1ts),
			df2s = (f2(s + eps, t) - f2ts), df2t = (f2(s, t + eps) - f2ts)
		var det = (df1s * df2t - df1t * df2s) / eps
		var ds = ( df2t * f1ts - df1t * f2ts) / det
		var dt = (-df2s * f1ts + df1s * f2ts) / det
		s -= ds
		t -= dt
	} while (--iterations && f1ts * f1ts + f2ts * f2ts > NLA.PRECISION)
	if (!iterations) {
		//console.log(f1ts * f1ts + f2ts * f2ts)
	 	return null
	}
	return V3(s, t, 0)
}
function newtonIterate2dWithDerivatives(f, g, startS, startT, iterations, dfds, dfdt, dgds, dgdt) {
	iterations = iterations || 4
	var s = startS, t = startT
	var eps = 1e-6
	do {
		/*
			| a b |-1                   |  d -b |
			| c d |   = 1 / (ad - bc) * | -c  a |
		 */
		var f1ts = f(s, t), f2ts = g(s, t)
		var df1s = dfds(s, t), df1t = dfdt(s, t),
			df2s = dgds(s, t), df2t = dgdt(s, t)
		// TODO: is this even more accurate?
		var det = df1s * df2t - df1t * df2s
		var ds = ( df2t * f1ts - df1t * f2ts) / det
		var dt = (-df2s * f1ts + df1s * f2ts) / det
		s -= ds
		t -= dt
	} while (--iterations && f1ts * f1ts + f2ts * f2ts > NLA.PRECISION / 32)
	if (!iterations) {
		//console.log(f1ts * f1ts + f2ts * f2ts)
	 	return null
	}
	return V3(s, t, 0)
}
function newtonIterate(f, startValue, steps) {
	var t = startValue
	var eps = 1e-8
	for (var i = 0; i < (steps || 4); i++) {
		var ft = f(t)
		var dft = (f(t + eps) - ft) / eps
		//console.log("ft / dft", ft / dft)
		t = t - ft / dft
	}
	return t
}
function newtonIterateWithDerivative(f, startValue, steps, df) {
	var t = startValue
	for (var i = 0; i < (steps || 4); i++) {
		var ft = f(t)
		var dft = df(t)
		//console.log("ft / dft", ft / dft)
		t = t - ft / dft
	}
	return t
}
function intersectionPCurveISurface(parametricCurve, searchStart, searchEnd, searchStep, implicitSurface) {
	assertNumbers(searchStart, searchEnd, searchStep)
	var iss = []
	var val = implicitSurface(parametricCurve(searchStart)), lastVal
	for (var t = searchStart + searchStep; t <= searchEnd; t += searchStep) {
		lastVal = val
		val = implicitSurface(parametricCurve(t))
		if (val * lastVal <= 0) {
			iss.push(newtonIterate(t => implicitSurface(parametricCurve(t)), t))
		}
	}
	return iss
}
function intersectionICurvePSurface(f0, f1, parametricSurface) {

}
function blugh(f, df, ddf, start, end, da) {
	var t = start, res = []
	while (t < end) {
		res.push(t)
		var cx = t, cy = f(t),
			dcx = 1, dcy = df(t),
			ddcx = 0, ddcy = ddf(t),
			div = Math.max(0.3, Math.abs(ddcy)),
			dt = da * (1 + dcy * dcy) / div
//		console.log(t, div, dt)
		t += dt
	}
	return res
}
function cassini(a, c) {
	return (x,y) => (x*x+y*y) * (x*x+y*y) - 2 * c * c * (x * x - y * y) - (a * a * a * a - c * c * c * c)
}
// TODO: V3.create instead of V3 where necessar
var drPs = [], drVs = []
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

	var face = B2.PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(10, 0), V3(10, 10), V3(0, 10)])
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



















function CustomPlane(anchor2, right, up, upStart, upEnd, rightStart, rightEnd, color, name) {
	var p = P3.forAnchorAndPlaneVectors(anchor2, right, up, CustomPlane.prototype);
	p.up = up;
	p.right = right;
	p.upStart = upStart;
	p.upEnd = upEnd;
	p.rightStart = rightStart;
	p.rightEnd = rightEnd;
	p.color = color;
	p.id = globalId++;
	p.name = name
	return p;
}
CustomPlane.prototype = Object.create(P3.prototype);
CustomPlane.prototype.constructor = CustomPlane;
Object.defineProperty(CustomPlane.prototype, "plane", { get: function () { return this } });
CustomPlane.prototype.toString = function() {
	return "Plane #" + this.id;
}
CustomPlane.prototype.what ="Plane"
CustomPlane.prototype.distanceTo = function (line) {
	return [
		L3(this.anchor.plus(this.right.times(this.rightStart)), this.up),
		L3(this.anchor.plus(this.right.times(this.rightEnd)), this.up),
		L3(this.anchor.plus(this.up.times(this.upStart)), this.right),
		L3(this.anchor.plus(this.up.times(this.upEnd)), this.right)].map(function (line2, line2Index) {
		var info = line2.infoClosestToLine(line);
		if ((isNaN(info.t) // parallel lines
			|| line2Index < 2 && this.upStart <= info.t && info.t <= this.upEnd
			|| line2Index >= 2 && this.rightStart <= info.t && info.t <= this.rightEnd)
			&& info.distance <= 16) {
			return info.s;
		} else {
			return Infinity;
		}
	}, this).min();
}
CustomPlane.forPlane = function (plane, color, name) {
	var p = P3(plane.normal, plane.w, CustomPlane.prototype)
	p.up = plane.normal.getPerpendicular().normalized()
	p.right = p.up.cross(p.normal)
	p.upStart = -500
	p.upEnd = 500
	p.rightStart = -500
	p.rightEnd = 500
	p.color = color || NLA.randomColor()
	p.id = globalId++
	p.name = name
	return p;
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