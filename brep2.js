"use strict";
["min", "max", "PI", "sqrt","pow","round"].forEach(function (propertyName) {
	/*if (window[propertyName]) {
	 throw new Error("already exists"+propertyName)
	 }*/
	window[propertyName] = Math[propertyName];
});
NLA.addOwnProperties(window, NLA)
/**
 * Created by aval on 21/12/2015.
 */
var M4 = NLA.Matrix4x4, V3 = NLA.Vector3, P3 = NLA.Plane3, L3 = NLA.Line3


var eps = 1e-5
var B2 = NLA.defineClass('B2', null,
	function (faces, infiniteVolume) {
		this.faces = faces
		assert(faces.every(f => f instanceof B2.Face), 'faces.every(f => f instanceof B2.Face)\n' + this.toString())
		this.infiniteVolume = !!infiniteVolume
	},
	{
		toMesh: function () {
			var mesh = new GL.Mesh({triangles: true, normals: true, lines: true})
			this.faces.forEach((face, i) => {
				face.addToMesh(mesh)
			})
			mesh.compile()
			return mesh
		},
		minus: function (brep2) {
			return this.intersection(brep2.flipped(), true, true)
		},
		plus: function (brep2) {
			return this.flipped().intersection(brep2.flipped(), true, true).flipped()
		},
		equals: function (brep) {
			return this.faces.length == brep.faces.length &&
				this.faces.every((face) => brep.faces.some((face2) => face.equals(face2)))
		},
		toString: function () {
			return `new B2([\n${this.faces.join(',\n').replace(/^/gm, '\t')}])`
		},
		toSource: function () {
			return `new B2([\n${this.faces.map(face => face.toSource()).join(',\n').replace(/^/gm, '\t')}])`
		},
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
						console.log(newLoop.edges[0].a.ss, newLoop.ccw ? contains : !contains)
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
				newFacesRecursive({edges: oldFace.edges, ccw: true, subloops: topLevelLoops})
				return true
			}
		},
		reconstituteFaces: function (oldFaces, edgeLooseSegments, faceMap, newFaces) {
			// reconstitute faces
			var insideEdges = []
			oldFaces.forEach(face => {
				var allOldFaceEdges = Array.prototype.concat.apply(face.edges, face.holes)
				console.log('reconstituting face', face.toString())
				var els = face.edges.map(edge => edgeLooseSegments.get(edge) && edgeLooseSegments.get(edge).join('\n')).map(s => '\n\n'+s).join()
				console.log('edgeLooseSegments', els)
				var looseEdges = faceMap.get(face)
				if (!looseEdges) {
					face.insideOutside = 'undecided'
				} else {
					face.insideOutside = 'part'
					// other brep does not intersect this face
					// assume it is outside
					console.log('looseEdges\n', looseEdges.map(e=>e.toString()).join('\n'))
					var loops = []
					var currentEdge
					var edgeCond = face instanceof B2.PlaneFace
						? (edge => edge.a.like(currentEdge.b))
						: (edge => (edge.curve == currentEdge.curve
										? NLA.equals(edge.aT, currentEdge.bT)
										: edge.a.like(currentEdge.b)))
					while (currentEdge = looseEdges.find(edge => !edge.visited)) {
						currentEdge.visited = true
						var startEdge = currentEdge, edges = [], i = 0
						var looseLoop = true
						do {
							console.log('currentEdge', currentEdge.b.ss, currentEdge.toSource())
							edges.push(currentEdge)
							// find next edge
							var possibleLooseEdges = looseEdges.filter(edge => edgeCond(edge)), looseSegments
							var possibleLooseEdgesCount = possibleLooseEdges.length
							allOldFaceEdges.forEach(edge => (looseSegments = edgeLooseSegments.get(edge)) && looseSegments.forEach(
									edge => edgeCond(edge) && possibleLooseEdges.push(edge)))
							var noSegments = false
							if (possibleLooseEdgesCount == possibleLooseEdges.length) {
								allOldFaceEdges.forEach(edge => edgeCond(edge) && possibleLooseEdges.push(edge))
								noSegments = true
							}
							var index = possibleLooseEdges.indexWithMax((edge, index) => currentEdge.bDir.angleRelativeNormal(edge.aDir, face.surface.normalAt(currentEdge.b)) + (index < possibleLooseEdgesCount ? NLA.PRECISION : 0))
							//console.log('possibleLooseEdges\n', possibleLooseEdges.map(e=>e.toString()).join('\n'), allOldFaceEdges.find(edge => edgeCond(edge)), index, possibleLooseEdges[index])
							// TODO assert(possibleLooseEdges.length < 2)
							currentEdge = possibleLooseEdges[index]
							if (index < possibleLooseEdgesCount) {
								possibleLooseEdges.forEach(possibleLooseEdge => possibleLooseEdge.isCoEdge(currentEdge) && (possibleLooseEdge.visited = true))
							} else {
								looseLoop = false
								noSegments && insideEdges.push(currentEdge)
							}
							assert(currentEdge)
						} while (++i < 10 && currentEdge != startEdge)
						if (10 == i) {
							assert(false, "too many")
						}
						console.log('finished loop')
						loops.push(edges)
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
		intersection: function (brep2, reconThis, reconB2, removeCoplanarSame, removeCoplanarOpposite) {
			var faceMap = new Map(), edgeMap = new Map()
			this.faces.forEach(face => {
				//console.log('face', face.toString())
				brep2.faces.forEach(face2 => {
					//console.log('face2', face2.toString())
					face.doo(face2, this, brep2, faceMap, edgeMap, removeCoplanarSame, removeCoplanarOpposite)
				})
			})
			var faces = []

			/*
			TODO:
			faceMap.forEach((faceLooses, face) => {
				faceLooses.forEach(edge => {
					face.edges.forEach(faceEdge => {
						var edgeT = faceEdge.getEdgeT(edge.a)
						if (undefined !== edgeT) {
							console.log("WAARGH", edge.a.ss, faceEdge.toString(), edgeT)
							NLA.mapAdd(edgeMap, faceEdge, {edgeT: edgeT, p: edge.a})
						}
					})
				})
			})
			*/
			var edgeLooseSegments = new Map()
			edgeMap.forEach((pointInfos, baseEdge) => {
				// TODO: make sure it works with loop
				// TODO: dont generate unnecessarry segments
				var looseSegments = []
				if (!baseEdge.reversed) {
					pointInfos.sort((a, b) => a.edgeT - b.edgeT)
				} else {
					pointInfos.sort((a, b) => b.edgeT - a.edgeT)
				}
				var startP = baseEdge.a, startDir = baseEdge.aDir, startT = baseEdge.aT
				for (var i = 0; i < pointInfos.length; i++) {
					var info = pointInfos[i]
					if (NLA.equals(info.edgeT, baseEdge.bT) || NLA.equals(info.edgeT, baseEdge.aT) || NLA.equals(info.edgeT, startT)) {
						continue
					}
					var pDir = baseEdge.tangentAt(info.edgeT)
					looseSegments.push(new baseEdge.constructor(baseEdge.curve, startP, info.p, startT, info.edgeT, null, startDir, pDir))
					startP = info.p
					startT = info.edgeT
					startDir = pDir
				}
				looseSegments.push(new baseEdge.constructor(baseEdge.curve, startP, baseEdge.b, startT, baseEdge.bT, null, startDir, baseEdge.bDir))
				if (baseEdge.b.like(V3(0, 4, 1))) {
					console.log('pointInfos', baseEdge.reversed, baseEdge.ss, pointInfos.map(pi => '\n'+pi.toSource()).join())
					console.log(looseSegments)
				}
				edgeLooseSegments.set(baseEdge, looseSegments)
			})

			reconThis && this.reconstituteFaces(this.faces, edgeLooseSegments, faceMap, faces, brep2.infiniteVolume)
			reconB2 && this.reconstituteFaces(brep2.faces, edgeLooseSegments, faceMap, faces, this.infiniteVolume)
			return new B2(faces, this.infiniteVolume && brep2.infiniteVolume)
		},
		transform: function (m4) {
			return new B2(this.faces.map(f => f.transform(m4)), this.infiniteVolume)
		},
		flipped: function () {
			return new B2(this.faces.map(f => f.flipped()), !this.infiniteVolume)
		}
	}
)
NLA.addTransformationMethods(B2.prototype)
B2.Face = function (surface, contour, holes) {
	this.assertChain(contour)
	//assert(surface.edgeLoopCCW(contour), surface.toString()+contour.join("\n"))
	assert(contour.every(f => f instanceof B2.Edge), 'contour.every(f => f instanceof B2.Edge)' + contour.toSource())
	contour.forEach(e => assert(surface.containsCurve(e.curve), e))
	holes && holes.forEach(hole => this.assertChain(hole))
	//holes && holes.forEach(hole => assert(!surface.edgeLoopCCW(hole)))
	assert (!holes || holes.constructor == Array, holes && holes.toString())
	this.surface = surface
	this.edges = contour
	this.holes = holes || []
	this.id = globalId++

}
B2.Face.prototype = NLA.defineObject(null, {
	transform: function (m4) {
		var newEdges = this.edges.map(e => e.transform(m4))
		var newHoles = this.holes.map(hole => hole.map(e => e.transform(m4)))
		return new this.constructor(this.surface.transform(m4), newEdges, newHoles)
	},
	assertChain: function (edges) {
		edges.forEach((edge, i) => {
			var j = (i + 1) % edges.length
			assert(edge.b.like(edges[j].a), `edges[${i}].b != edges[${j}].a (${edges[i].b.ss} != ${edges[j].a.ss})`)
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
		return `new ${this.name}(${this.surface.toSource()}, [${this.edges.map(e => '\n\t' + e.toSource()).join(',')}])`
	},
	equals: function (face) {
		var edgeCount = this.edges.length

		return this.surface.equalsSurface(face.surface) &&
				this.edges.length == face.edges.length &&
				NLA.arrayRange(0, edgeCount, 1)
					.some(offset => this.edges.every((edge, i) => edge.equals(face.edges[(offset + i) % edgeCount])))
	},
	addEdgeLines: function (mesh) {
		assert(false, "buggy, fix")
		var vertices = this.edges.map(edge => edge.getVerticesNo0()).concatenated(), mvl = mesh.vertices.length
		for (var i = 0; i < vertices.length; i++) {
			mesh.vertices.push(vertices[i])
			mesh.lines.push(mvl + i, mvl + (i + 1) % vertices.length)

		}
	}
})
NLA.addTransformationMethods(B2.Face.prototype)
B2.FaceOnPISurface = function (surface, edges) {
	assert(surface.parametricFunction && surface.implicitFunction, 'surface.parametricFunction && surface.implicitFunction')
	this.surface = surface
	this.edges = edges
}
B2.FaceOnPISurface.prototype = NLA.defineObject(B2.Face.prototype, {
	addToMesh: function (mesh) {
	}
})
/*
Object.prototype.getIntersectionsWithPlane= function (p) {
	assert(false, this.name + '.getIntersectionsWithPlane' + this.constructor.name)
}
*/
B2.Edge = function () {}
B2.Edge.prototype = NLA.defineObject(null, {
	toString: function (f) {
		return `new ${this.name}(${this.curve.toString(f)}, ${this.a}, ${this.b}, ${this.aT}, ${this.bT}, null, ${this.aDir}, ${this.bDir})`
	},
	getIntersectionsWithPlane: function (p) {
		assert(false, this.name + '.getIntersectionsWithPlane')
	}
})
B2.PlaneFace = NLA.defineClass('B2.PlaneFace', B2.Face,
	function (planeSurface, contour, holes, name) {
		B2.Face.call(this, planeSurface, contour, holes, name)
		assert (planeSurface instanceof PlaneSurface)
	},
	{
		addToMesh: function (mesh) {
			var mvl = mesh.vertices.length
			var normal = this.surface.plane.normal
			var vertices = this.edges.map(edge => edge.getVerticesNo0()).concatenated()
			console.log(vertices.map(V3.ss).join())
			for (var i = 0; i < vertices.length; i++) { mesh.lines.push(mvl + i, mvl + (i + 1) % vertices.length) }
			var holeStarts = []
			this.holes.forEach(hole => {
				holeStarts.push(vertices.length)
				vertices.pushAll(hole.map(edge => edge.getVerticesNo0()).concatenated())
			})
			var triangles = triangulateVertices(normal, vertices, holeStarts).map(index => index + mvl)
			Array.prototype.push.apply(mesh.vertices, vertices)
			Array.prototype.push.apply(mesh.triangles, triangles)
			Array.prototype.push.apply(mesh.normals, NLA.arrayFromFunction(vertices.length, i => normal))
		},
		containsPoint: function (p) {
			assertVectors (p)
			return this.surface.edgeLoopContainsPoint(this.edges, p)
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
							/*console.log("NEWCURVE", p.ss, face2.surface.normalAt(p).cross(this.surface.normalAt(p)).ss, 'nct', newCurve.tangentAt(newCurveT).ss,
								face2.surface.normalAt(p).cross(this.surface.normalAt(p)).dot(newCurve.tangentAt(newCurveT)) > 0,
								'newCurveT',newCurveT)
						}
						var ov = edge.tangentAt(edgeT).cross(this.surface.normalAt(p))
						var ct = newCurve.tangentAt(newCurveT)
						console.log("ov", p.ss,edge.tangentAt(edgeT).ss, this.surface.normalAt(p).ss, ov.ss, ct.ss, ov.dot(ct) > 0)
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
			 var newCurve = this.surface.getIntersectionCurve(face2.surface)
			 // get intersections of newCurve with other edges of face and face2
			 var ps1 = []
			 this.edges.forEach(edge => {
			 var iss = edge.getIntersectionsWithISurface(face2.surface)
			 for (var i = 0; i < iss.length; i++) {
			 var p = edge.curve.at(iss[i])
			 var ov = edge.pointTangent(p).cross(this.surface.normalAt(p))
			 var ct = newCurve.pointTangent(p)
			 //console.log("ov", p.ss,edge.pointTangent(p).ss, this.surface.normalAt(p).ss, ov.ss, ct.ss)
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
		doo: function (face2, thisBrep, face2Brep, faceMap, edgeMap, removeCoplanarSame, removeCoplanarOpposite) {
			if (face2 instanceof B2.RotationFace) {
				if (this.surface.isCoplanarTo(face2.surface)) { return }

				// get intersections
				var newCurves = face2.surface.getIntersectionsWithSurface(this.surface)

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
					console.log('iscurve', newCurve.toString(), newCurve.tangentAt(ps1[0].t).ss)
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
				this.dooPlaneFace(face2, thisBrep, face2Brep, faceMap, edgeMap, removeCoplanarSame, removeCoplanarOpposite)
			}
		},
		dooPlaneFace: function (face2, thisBrep, face2Brep, faceMap, edgeMap, removeCoplanarSame, removeCoplanarOpposite) {
			assert(face2 instanceof B2.PlaneFace)
			// get intersection
			var thisPlane = this.surface.plane, face2Plane = face2.surface.plane
			if (thisPlane.isParallelToPlane(face2Plane)) { return }
			var intersectionLine = L3.fromPlanes(thisPlane, face2Plane)
			var thisDir = (this.surface.normalAt(null).cross(face2.surface.normalAt(null)).dot(intersectionLine.dir1) > 0)
			// get intersections of newCurve with other edges of face and face2
			var ps1 = getFacePlaneIntersectionSs(thisBrep, this, intersectionLine, face2Plane, true, false)
			var ps2 = getFacePlaneIntersectionSs(face2Brep, face2, intersectionLine, thisPlane, false, false)

			if (ps1.length == 0 || ps2.length == 0) {
				// faces to not intersect
				return
			}
			console.log(''+thisPlane+face2Plane)

			var segments1 = getBlug(ps1.filter(ps => ps.caseB && !ps.hideOnFace), ps2.filter(ps => ps.caseB), intersectionLine)
			var segments2 = getBlug(ps1.filter(ps => ps.caseA), ps2.filter(ps => ps.caseB && !ps.hideOnFace), intersectionLine)
			// TODO: getCanon() TODO TODO TODO
			if (intersectionLine.dir1.isParallelTo(V3.Z)) {
				console.log(intersectionLine.toString(), ps1[0].t, ps2[1].t, NLA.equals(ps1[0].t, ps2[1].t))
				console.log('segments', segments1.toSource(), segments2.toSource())
				console.log('ps1\n', ps1.map(m => m.toSource()).join('\n'), '\nps2\n', ps2.map(m => m.toSource()).join('\n'))
			}
			ps1.forEach(ps => ps.used && mapAdd(edgeMap, ps.edge, ps))
			ps2.forEach(ps => ps.used && mapAdd(edgeMap, ps.edge, ps))
			segments1.forEach(segment => {
				console.log('segment', segment.toString(), thisDir)
				mapAdd(faceMap, this, thisDir ? segment : segment.flipped())
			})
			segments2.forEach(segment => {
				console.log('segment', segment.toString(), thisDir)
				mapAdd(faceMap, face2, thisDir ? segment.flipped() : segment)
			})
		}
	}
)

B2.PlaneFace.forVertices = function (planeSurface, vs) {
	if (planeSurface instanceof P3) {
		planeSurface = new PlaneSurface(planeSurface)
	}
	assert(isCCW(vs, planeSurface.plane.normal), 'isCCW(vs, planeSurface.plane.normal)')
	var edges = vs.map((a, i) => {
		var b = vs[(i + 1) % vs.length]
		return StraightEdge.throughPoints(a, b)
	})
	return new B2.PlaneFace(planeSurface, edges, this.infiniteVolume)
}
function facesWithEdge(edge, faces) {
	return arrayFilterMap(faces, (face) => {
		var matchingEdge = face.edges.find(e => e.isCoEdge(edge))
		if (matchingEdge) {
			return {face: face, reversed: !edge.a.like(matchingEdge.a), angle: NaN, normalAtEdgeA: null, edge: matchingEdge}
		}
	})
}
function getFacePlaneIntersectionSs(brep, brepFace, line, plane2, removeCoplanarSame, removeCoplanarOpposite) {
	var facePlane = brepFace.surface.plane
	var ps = []
	var loops = brepFace.holes.concat([brepFace.edges])
	loops.forEach(loop => {

		var colinearSegments = loop.map((edge) => edge.colinearToLine(line))
		var testVector = plane2.projectedVector(facePlane.normal)
		var intersectionLinePerpendicular = line.dir1.cross(facePlane.normal)
		var colinearSegmentsInside = loop.map((edge, i) => colinearSegments[i] &&
			(splitsVolumeEnclosingFaces(brep, edge, testVector, plane2.normal, removeCoplanarSame, removeCoplanarOpposite)
			!= splitsVolumeEnclosingFaces(brep, edge, testVector.negated(), plane2.normal, removeCoplanarSame, removeCoplanarOpposite))
		)
		var colinearSegmentsInsideCaseTrue = [], colinearSegmentsInsideCaseFalse = []
		for (var i = 0; i < loop.length; i++) {
			if (colinearSegments[i]) {
				var edge = loop[i]
				var csi1 = splitsVolumeEnclosingFaces(brep, edge, testVector, plane2.normal, removeCoplanarSame, removeCoplanarOpposite)
				var csi2 = splitsVolumeEnclosingFaces(brep, edge, testVector.negated(), plane2.normal, removeCoplanarSame, removeCoplanarOpposite)
				var a = INSIDE == csi1, b = INSIDE  == csi1 || COPLANAR_SAME == csi1
				var c = INSIDE == csi2, d = INSIDE  == csi2 || COPLANAR_SAME == csi2
				colinearSegmentsInsideCaseTrue[i] = b != d
				colinearSegmentsInsideCaseFalse[i] = a != c
			}
		}

		loop.forEach((edge, i, edges) => {
			var j = (i + 1) % edges.length, nextEdge = edges[j]
			//console.log(edge.toSource()) {p:V3(2, -2.102, 0),
			if (colinearSegments[i]) {
				// edge colinear to intersection
				var outVector = edge.bDir.cross(facePlane.normal)
				var insideNext = outVector.dot(nextEdge.aDir) > 0
				var caseA = insideNext != colinearSegmentsInsideCaseTrue[i],
					caseB = insideNext != colinearSegmentsInsideCaseFalse[i]
				var colinearSegmentOutsideVector = edge.aDir.cross(facePlane.normal)
				var displayOnFace = colinearSegmentOutsideVector.dot(plane2.normal) > 0
				if (caseA || caseB || displayOnFace != insideNext) {
					ps.push({p: edge.b, insideDir: null, t: line.pointLambda(edge.b), edge: edge, edgeT: edge.bT,
						caseA: caseA, caseB: caseB, colinear: true, hideOnFace: displayOnFace == insideNext})
					//console.log('colinear')
				}
			} else {
				var edgeTs = edge.getIntersectionsWithPlane(plane2)
				for (var k = 0; k < edgeTs.length; k++) {
					var edgeT = edgeTs[k]
					if (edgeT == edge.bT) {
						// endpoint lies on intersection line
						console.log('endpoint lies on intersection line',
							intersectionLinePerpendicular.dot(edge.bDir) , intersectionLinePerpendicular.dot(nextEdge.aDir),
							intersectionLinePerpendicular.dot(edge.bDir) * intersectionLinePerpendicular.dot(nextEdge.aDir), intersectionLinePerpendicular.ss,
							edge.bDir.ss, nextEdge.aDir.ss)
						if (colinearSegments[j]) {
							// next segment is colinear
							// we need to calculate if the section of the plane intersection line BEFORE the colinear segment is
							// inside or outside the face. It is inside when the colinear segment out vector and the current segment vector
							// point in the same direction (dot > 0)
							// TODO: UUUH?
							var colinearSegmentOutsideVector = nextEdge.aDir.cross(facePlane.normal)
							var insideFaceBeforeColinear = colinearSegmentOutsideVector.dot(edge.bDir) < 0
							var caseA = insideFaceBeforeColinear != colinearSegmentsInsideCaseTrue[j],
								caseB = insideFaceBeforeColinear != colinearSegmentsInsideCaseFalse[j]
							var displayOnFace = colinearSegmentOutsideVector.dot(plane2.normal) > 0
							// if the "inside-ness" changes, add intersection point
							//console.log("segment end on line followed by colinear", insideFaceBeforeColinear != colinearSegmentInsideFace, nextSegmentOutsideVector)
							if (caseA || caseB || displayOnFace != insideFaceBeforeColinear) {
								ps.push({p: edge.b, insideDir: null, t: line.pointLambda(edge.b), edge: edge, edgeT: edge.bT
									, caseA: caseA, caseB: caseB, colinear: true, hideOnFace: displayOnFace == insideFaceBeforeColinear})
								//console.log('next colinear')
							}
						} else if (intersectionLinePerpendicular.dot(edge.bDir) * intersectionLinePerpendicular.dot(nextEdge.aDir) > 0) {
							// next segment is not colinear and ends on different side
							ps.push({p: edge.b, insideDir: null, t: line.pointLambda(edge.b), edge: edge, edgeT: edge.bT, caseA: true, caseB: true})
							//console.log('end on line, next other side')
						}
					} else if (edgeT != edge.aT) {
						// edge crosses is line, neither starts nor ends on it
						var p = edge.curve.at(edgeT)
						assert(line.containsPoint(p), edge.toString() + p+edgeT)
						ps.push({p: p, insideDir: null, t: line.pointLambda(p), edge: edge, edgeT: edgeT, caseA: true, caseB: true})
						console.log('middle')
					}
				}
			}
		})
	})
	ps.sort((a, b) => a.t - b.t || -a.insideDir.dot(line.dir1))
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
		var j = (i + 1) % edges.length, nextEdge = edges[j]
		//console.log(edge.toSource()) {p:V3(2, -2.102, 0),
		if (colinearSegments[i]) {
			assert(false)
			// edge colinear to intersection
			var outVector = edge.bDir.cross(faceSurface.normal)
			var insideNext = outVector.dot(nextEdge.aDir) > 0
			var caseA = insideNext != colinearSegmentsInsideCaseTrue[i],
				caseB = insideNext != colinearSegmentsInsideCaseFalse[i]
			var colinearSegmentOutsideVector = edge.aDir.cross(faceSurface.normal)
			var displayOnFace = colinearSegmentOutsideVector.dot(surface2.normalAt(edge.b)) > 0
			if (caseA || caseB || displayOnFace != insideNext) {
				ps.push({p: edge.b, insideDir: null, t: isCurve.pointLambda(edge.b), edge: edge, edgeT: edge.bT,
					caseA: caseA, caseB: caseB, colinear: true, hideOnFace: displayOnFace == insideNext})
				//console.log('colinear')
			}
		} else {
			var edgeTs = edge.getIntersectionsWithSurface(surface2)
			for (var k = 0; k < edgeTs.length; k++) {
				var edgeT = edgeTs[k]
				if (edgeT == edge.bT) {
					assert(false)
					// endpoint lies on intersection isCurve
					console.log('endpoint lies on intersection isCurve',
						intersectionLinePerpendicular.dot(edge.bDir) , intersectionLinePerpendicular.dot(nextEdge.aDir),
						intersectionLinePerpendicular.dot(edge.bDir) * intersectionLinePerpendicular.dot(nextEdge.aDir), intersectionLinePerpendicular.ss,
						edge.bDir.ss, nextEdge.aDir.ss)
					var intersectionLinePerpendicular = curve.tangentAt(edge.b).cross(faceSurface.normalAt(edge.b))
					if (colinearSegments[j]) {
						// next segment is colinear
						// we need to calculate if the section of the plane intersection isCurve BEFORE the colinear segment is
						// inside or outside the face. It is inside when the colinear segment out vector and the current segment vector
						// point in the same direction (dot > 0)
						// TODO: UUUH?
						var colinearSegmentOutsideVector = nextEdge.aDir.cross(faceSurface.normal)
						var insideFaceBeforeColinear = colinearSegmentOutsideVector.dot(edge.bDir) < 0
						var caseA = insideFaceBeforeColinear != colinearSegmentsInsideCaseTrue[j],
							caseB = insideFaceBeforeColinear != colinearSegmentsInsideCaseFalse[j]
						var displayOnFace = colinearSegmentOutsideVector.dot(surface2.normalAt(edge.a)) > 0
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

					var onCurve = isCurves.length, isCurve
					var p = edge.curve.at(edgeT)
					while (--onCurve >= 0 && !(isCurve = isCurves[onCurve]).containsPoint(p)) {}
					console.log('edgeT', edgeT, 'p', p.ss, edge, onCurve)
					if (onCurve < 0) {
						assert (false, p.ss)
					}
					assert(isCurve.containsPoint(p))
					var ov = edge.tangentAt(edgeT).cross(faceSurface.normalAt(p))
					var newCurveT = isCurve.pointLambda(p, ov)
					var ct = isCurve.tangentAt(newCurveT)
					if (ov.dot(ct) > 0) ct = ct.negated()
					pss[onCurve].push({p: p, insideDir: ct, t: newCurveT, edge: edge, edgeT: edgeT, caseA: true, caseB: true})
					console.log('middle')
				}
			}
		}
	})
	pss.forEach(ps => ps.sort((a, b) => a.t - b.t || assert(false, a.t + ' '+b.t+' '+a.p.ss)))
	return pss
}

var INSIDE = 0, OUTSIDE = 1, COPLANAR_SAME = 2, COPLANAR_OPPOSITE= 3
/**
 *
 * @param brep BREP to check
 * @param edge edge to check
 * @param dirAtEdgeA the direction vector to check
 * @param faceNormal If dirAtEdgeA doesn't split a volume, but is along a face, the returned value depends on wether
 * that faces normal points in the same direction as faceNormal
 * @param coplanarSameInside
 * @param coplanarOppositeInside
 * @returns {*}
 */
function splitsVolumeEnclosingFaces(brep, edge, dirAtEdgeA, faceNormal, coplanarSameInside, coplanarOppositeInside) {
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
	relFaces.sort((a, b) => a.angle - b.angle)
	assert(relFaces.length != 0)
	//console.log(relFaces.map(f => f.toSource()).join('\n'))

	if (NLA.isZero(relFaces[0].angle)) {
		var coplanarSame = relFaces[0].normalAtEdgeA.dot(faceNormal) > 0
		return coplanarSame ? COPLANAR_SAME : INSIDE
	} else {
		return !relFaces[0].reversed
	}
}
B2.RotationFace = function (rot, contour, holes, name) {
	B2.Face.call(this, rot, contour, holes, name)
	//assert(rot instanceof RotationReqFofZ)
}
B2.RotationFace.prototype = NLA.defineObject(B2.Face.prototype, {
	constructor: B2.RotationFace,
	name: 'B2.RotationFace',
	unrollLoop: function (edgeLoop) {
		var vs = []
		edgeLoop.forEach((edge, e) => {
			var reverseFkt = this.surface.pointToParameterFunction()
			var hint = edge.bDir
			if (edge instanceof StraightEdge && edge.curve.dir1.isParallelTo(this.surface.dir)) {
				hint = this.surface.baseEllipse.f1.times(this.surface.normalDir).cross(edge.bDir)
			}
			edge.getVerticesNo0().forEach(p => vs.push(reverseFkt(p, hint)))
			console.log(edge.ss, "\n", edge.getVerticesNo0().map(p => reverseFkt(p, hint)).join())
		})
		return vs
	},
	addToMesh: function (mesh) {
		console.log("mlsadkl")
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
				console.log(v0.ss, v1.ss)
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
				console.log(ribs.map(r=>r.toSource()).join('\n'))
			})
		})
		var vertices = [], triangles = [], normals = []
		for (var i = 0; i < ribs.length; i++) {
			var ribLeft = ribs[i], ribRight = ribs[(i + 1) % ribs.length]
			assert(ribLeft.right.length == ribRight.left.length)
			for (var j = 0; j < ribLeft.right.length; j++) {
				vertices.push(f(ribLeft.value, ribLeft.right[j]), f(ribRight.value, ribRight.left[j]))
				normals.push(normalF(ribLeft.value, ribLeft.right[j]), normalF(ribRight.value, ribRight.left[j]))
			}
		}
		console.log(ribs.map(r=>r.toSource()).join('\n'))
		var vss = vertices.length, detailVerticesStart = vss
		var zInterval = maxZ - minZ, zStep = zInterval / zSplit
		var detailZs = NLA.arrayFromFunction(zSplit - 1, i => minZ + (1 + i) * zStep)
		for (var i = 0; i < ribs.length; i++) {
			var d = ribs[i].value
			for (var j = 0; j < detailZs.length; j++) {
				vertices.push(f(d, detailZs[j]))
				normals.push(normalF(d, detailZs[j]))
			}
		}
		//console.log('detailVerticesStart', detailVerticesStart, 'vl', vertices.length, vertices.length - detailVerticesStart, ribs.length)
		// finally, fill in the ribs
		var vsStart = 0
		var flipped2 = this.surface.normalDir == 1
		//for (var i = 0; i < 1; i++) {
		var end = closed ? ribs.length : ribs.length - 1
		for (var i = 0; i < end; i++) {
			var ipp = (i + 1) % ribs.length
			var inside = false, colPos = 0, ribLeft = ribs[i], ribRight = ribs[(i + 1) % ribs.length]
			for (var j = 0; j < detailZs.length + 1; j++) {
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
		//console.log('trinagle', triangles.max(), vertices.length, triangles.length, triangles.toSource(), triangles.map(i => vertices[i].ss).toSource() )
		triangles = triangles.map(index => index + mesh.vertices.length)
		//assert(normals.every(n => n.hasLength(1)), normals.find(n => !n.hasLength(1)).length() +" "+normals.findIndex(n => !n.hasLength(1)))
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
function pushQuad(triangles, flipped, a, b, c, d) {
	if (flipped) {
		triangles.push(a, c, b,
			b, c, d)
	} else {
		triangles.push(a, b, c,
			b, d, c)
	}
}
B2.PCurveEdge = NLA.defineClass('B2.PCurveEdge', B2.Edge,
	function (curve, a, b, aT, bT, flippedOf, aDir, bDir) {
		assertNumbers(aT, bT)
		assertVectors(a, b, aDir, bDir)
		assert(curve instanceof L3 || curve instanceof EllipseCurve, curve)
		assert(!curve.isValidT || curve.isValidT(aT) && curve.isValidT(bT), aT + ' ' + bT)
		assert(curve.at(aT).like(a))
		assert(curve.at(bT).like(b))
		this.curve = curve
		this.a = a
		this.b = b
		this.aT = aT
		this.bT = bT
		this.aDir = aDir
		this.bDir = bDir
		this.canon = flippedOf
		this.reversed = this.aDir.dot(curve.tangentAt(aT)) < 0
		assert(this.reversed != aT < bT, aT+' '+bT)
		this.id = globalId++
	},
	{
		getVerticesNo0: function () {
			return this.curve.asklkjas(this.aT, this.bT, this.a, this.b, this.reversed, false)
		},
		get points() {
			return this.curve.asklkjas(this.aT, this.bT, this.a, this.b, this.reversed, true)
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
		getIntersectionsWithISurface: function (is) {
			assert (is.implicitFunction)
			var start = min(this.aT, this.bT), end = max(this.aT, this.bT)
			return intersectionPCurveISurface(t => this.curve.at(t), start, end, 0.1, is.implicitFunction())
		},
		getIntersectionsWithSurface: function (surface) {
			return this.curve.getIntersectionsWithSurface(surface).filter(edgeT => {
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
		getIntersectionsWithPlane: function (surface) {
			console.log(this.curve, this.curve.constructor.name)
			return this.curve.getIntersectionsWithPlane(surface).filter(edgeT => {
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
		colinearToLine: function (line) {
			return this.curve.equals(line)
		},
		isCoEdge: function (edge) {
			// TODO: optimization with flippedOf etc
			return this == edge ||
				this.curve.isColinearTo(edge.curve) && (
					this.a.like(edge.a) && this.b.like(edge.b)
					|| this.a.like(edge.b) && this.b.like(edge.a)
				)
		},
		likeEdge: function (edge) {
			return edge.constructor == StraightEdge && this.a.like(edge.a) && this.b.like(edge.b)
		}
	}
)
var StraightEdge = NLA.defineClass('StraightEdge', B2.Edge,
	function (line, a, b, aT, bT, flippedOf) {
		assertNumbers(aT, bT)
		assertVectors(a, b)
		assert(line instanceof L3)
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
		getIntersectionsWithISurface: function (is) {
			assert (is.implicitFunction)
			var start = min(this.aT, this.bT), end = max(this.aT, this.bT)
			return intersectionPCurveISurface(t => this.curve.at(t), start, end, 0.1, is.implicitFunction())
		},
		getIntersectionsWithPlane: function (plane) {
			var minT = min(this.aT, this.bT), maxT = max(this.aT, this.bT)
			var edgeT = this.curve.intersectWithPlaneLambda(plane)
			edgeT = NLA.snapTo(edgeT, this.aT)
			edgeT = NLA.snapTo(edgeT, this.bT)
			return (minT <= edgeT && edgeT <= maxT) ? [edgeT] : []
		},
		getIntersectionsWithSurface: function (surface) {
			if (surface instanceof PlaneSurface) {
				return this.getIntersectionsWithPlane(surface.plane)
			} else if (surface instanceof CylinderSurface) {
				var minT = min(this.aT, this.bT), maxT = max(this.aT, this.bT)
				return surface.intersectionWithLine(this.curve)
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
		colinearToLine: function (line) {
			return this.curve.equals(line)
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
 *      a^2 + b^2 = 1
 * This can be understood as the intersection of the unit circle with a line.
 * @param a double
 * @param b double
 * @param c double
 * @returns {x1, y1, x2, y2} with x1 >= x2 and y1 <= y2
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

function getIntersectionSegments (ps1, ps2, in1, in2, constructor, curve) {
	var currentSegment
	assert (!(in1 && in2), '!(in1 && in2) '+in1+' '+in2)
	console.log('in', in1, in2)
	// generate overlapping segments
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
function PlaneSurface(plane, right, up) {
	assert(plane instanceof P3)
	this.plane = plane
	this.up = up || plane.normal.getPerpendicular().normalized()
	this.right = right || this.up.cross(this.plane.normal).normalized()
	assert(this.right.cross(this.up).like(this.plane.normal))
}
PlaneSurface.throughPoints = function (a, b, c) {
	return new PlaneSurface(P3.throughPoints(a, b, c))
}
PlaneSurface.prototype = NLA.defineObject(null, {
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
	intersectionCurveWithImplicitSurface: function (implicitSurface) {
		assert (implicitSurface.implicitFunction, 'implicitSurface.implicitFunction')
		return new CurvePI(this, implicitSurface)
	},
	getIntersectionCurve: function (surface2) {
		// prefer other surface to be the paramteric one
		if (surface2.implicitFunction) {
			return new CurvePI(this, surface2)
		} else if (surface2.parametricFunction) {
			return new CurvePI(surface2, this)
		}
	},
	edgeLoopCCW: function (contour) {
		var totalAngle = 0
		for (var i = 0; i < contour.length; i++) {
			var ipp = (i + 1) % contour.length
			var edge = contour[i], nextEdge = contour[ipp]
			if (edge.curve instanceof EllipseCurve) {
				totalAngle += edge.rotViaPlane(this.plane.normal)
				console.log(edge.toString(), edge.rotViaPlane(this.plane.normal))
			}
			totalAngle += edge.bDir.angleRelativeNormal(nextEdge.aDir, this.plane.normal)
		}
		return totalAngle > 0
	},
	edgeLoopContainsPoint: function (contour, p) {
		assertVectors (p)
		var dir = this.right
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
				var edgeTs = edge.getIntersectionsWithPlane(plane2)
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
function RotationReqFofZ(l3Axis, FofR, minZ, maxZ) {
	assert(l3Axis instanceof L3)
	this.l3Axis = l3Axis
	this.FofR = FofR
	this.minZ = minZ
	this.maxZ = maxZ
}
RotationReqFofZ.prototype = {
	toMesh: function (zStart, zEnd, count) {
		var zInterval = zEnd - zStart, zStep = zInterval / (count - 1)
		var vertices = NLA.arrayFromFunction(count, i => (z = zStart + i * zStep, V3.create(this.FofR(z), 0, z)))
		var normals = NLA.arrayFromFunction(count, i => {
			var z = zStart + i * zStep
			var fz = this.FofR(z)
			var dfz = (this.FofR(z + eps) - fz) / eps
			return V3.create(1, 0, -dfz).normalized()
		})
		var z = this.l3Axis.dir1, x = z.getPerpendicular().normalized(), y = z.cross(x)
		var matrix = M4.forSys(x, y, z, this.l3Axis.anchor);
		vertices = matrix.transformedPoints(vertices)
		normals = matrix.inversed().transposed().transformedVectors(normals).map(v => v.normalized())
		return rotationMesh(vertices, this.l3Axis, 2 * Math.PI, 64, true, normals)
	},
	parametricNormal: function () {
		var z = this.l3Axis.dir1, x = z.getPerpendicular().normalized(), y = z.cross(x)
		var matrix = M4.forSys(x, y, z, this.l3Axis.anchor).inversed().transposed()
		return (d, z) => {
			var fz = this.FofR(z)
			var dfz = (this.FofR(z + eps) - fz) / eps
			return matrix.transformVector(V3.create(cos(d), sin(d), -dfz)).normalized()
		}
	},
	normalAt: function (p) {
		var pmPoint = this.pointToParameterFunction()(p)
		return this.parametricNormal()(pmPoint.x, pmPoint.y)
	},
	parametricFunction: function () {
		var z = this.l3Axis.dir1, x = z.getPerpendicular().normalized(), y = z.cross(x)
		var matrix = M4.forSys(x, y, z, this.l3Axis.anchor)
		var f = this.FofR
		return function (d, z) {
			var radius = f(z)
			return matrix.transformPoint(V3.create(radius * cos(d), radius * sin(d), z))
		}
	},
	implicitFunction: function () {
		var z = this.l3Axis.dir1, x = z.getPerpendicular().normalized(), y = z.cross(x)
		var matrix = M4.forSys(x, y, z, this.l3Axis.anchor)
		var matrixInverse = matrix.inversed()
		var f = this.FofR
		return function (pWC) {
			var p = matrixInverse.transformPoint(pWC)
			var radiusLC = Math.sqrt(p.x * p.x + p.y * p.y)
			return f(p.z) - radiusLC
		}
	},
	boundsFunction: function () {
		var z = this.l3Axis.dir1, x = z.getPerpendicular().normalized(), y = z.cross(x)
		var matrix = M4.forSys(x, y, z, this.l3Axis.anchor)
		var matrixInverse = matrix.inversed()
		var f = this.FofR, minZ = this.minZ, maxZ = this.maxZ
		return function (pWC) {
			var z = matrixInverse.transformPoint(pWC).z
			return minZ <= z && z <= maxZ
		}
	},
	pointToParameterFunction: function (p) {
		var z = this.l3Axis.dir1, x = z.getPerpendicular().normalized(), y = z.cross(x)
		var matrix = M4.forSys(x, y, z, this.l3Axis.anchor)
		var matrixInverse = matrix.inversed()
		var f = this.FofR
		return function (pWC) {
			var p = matrixInverse.transformPoint(pWC)
			return V3.create(atan2(p.y, p.x), p.z, 0)
		}
	},
	getIntersectionCurve: function (surface2) {
		// prefer other surface to be the paramteric one
		if (surface2.parametricFunction) {
			return new CurvePI(surface2, this)
		} else if (surface2.implicitFunction) {
			return new CurvePI(this, surface2)
		}
	}
}

function curvePoint(implicitCurve, startPoint) {
	var eps = 1e-5
	var p = startPoint
	for (var i = 0; i < 4; i++) {
		var fp = implicitCurve(p.x, p.y)
		var dfpdx = (implicitCurve(p.x + eps, p.y) - fp) / eps,
			dfpdy = (implicitCurve(p.x, p.y + eps) - fp) / eps
		var scale = fp / (dfpdx * dfpdx + dfpdy * dfpdy)
		//console.log(p.ss)
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
	// TODO gleichm¨aßige Verteilung der Punkte
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


function cylinderPoints (l3Axis, radius) {
	assert(l3Axis instanceof L3)
	var z = l3Axis.dir1, x = z.getPerpendicular().normalized(), y = z.cross(x)
	var matrix = M4.forSys(x.times(radius), y.times(radius), z, l3Axis.anchor)
	return function (d, z) {
		return matrix.transformPoint(V3.create(cos(d), sin(d), z))
	}
}
function cylinderImplicit(l3Axis, radius) {
	assert(l3Axis instanceof L3)
	var z = l3Axis.dir1, x = z.getPerpendicular().normalized(), y = z.cross(x)
	var matrix = M4.forSys(x.times(radius), y.times(radius), z, l3Axis.anchor)
	var matrixInverse = matrix.inversed()
	assert(matrixInverse.times(matrix).isIdentity(NLA.PRECISION))
	return function (pWC) {
		var p = matrixInverse.transformPoint(pWC)
		var radiusLC = Math.sqrt(p.x * p.x + p.y * p.y)
		return 1 - radiusLC
	}
}
function newtonIterate2d(f1, f2, startS, startT, iterations) {
	iterations = iterations || 4
	var s = startS, t = startT
	var eps = 1e-5
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
function newtonIterate(f, startValue) {
	var t = startValue
	var eps = 1e-5
	for (var i = 0; i < 4; i++) {
		var ft = f(t)
		var dft = (f(t + eps) - ft) / eps
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
var drPs = []
addGenerator()
function initB2() {
	/*
	var c1 = EllipseCurve.circle(5), c2 = EllipseCurve.circle(5, V3(3, 0))
	var test = new EllipseCurve(V3(6, 1, 0), V3(3, 1, 0), V3(4, 0, 0))
	dMesh = new GL.Mesh()
	var cyl = new CylinderSurface(new EllipseCurve(V3.ZERO, V3(5, 5, 0), V3(0, 5, 0)), V3.Z, 1)
	var ell = new CylinderSurface(new EllipseCurve(V3.ZERO, V3(5, 5, 0), V3(0, 2, 0)), V3.Z, 1).rotateX(PI/3)

	aMesh = cyl.toMesh()
	bMesh = ell.toMesh()
	c1.intersectWithEllipse(test)
	dMesh.compile()
	*/
	eyePos = V3(0, 100, 100)
	eyeFocus = V3(0, 100, 0)
	eyeUp = V3(0, 1, 0)
	zoomFactor = 0.5
	setupCamera()
	paintScreen2()

	var face = B2.PlaneFace.forVertices(new PlaneSurface(P3.XY), [V3(0, 0, 0), V3(10, 0, 0),V3(10, 10, 0), V3(4, 4,0), V3(0, 10, 0)])
	var extrusion = B2.extrudeVertices([V3(0, 0, 0), V3(10, 0, 0),V3(10, 10, 0), V3(4, 4,0), V3(0, 10, 0)], P3.XY, V3(0, 0, 4), "ex0")
	var wideBox = B2.box(10, 10, 5)
	//wideBox = new B2(wideBox.faces.slice(2, 3))
	//wideBox = B2.box(5, 5, 5).flipped()
	//var wideBox = new B2(extrusion.faces.slice(0, 1))

	// top.addToMesh(new GL.Mesh({normals: true}));bottom.addToMesh(new GL.Mesh({normals: true}))
	//tallBox = B2.box(5, 5, 10).translate(0, -5, 0).rotateX(-0.2).translate(0, 3, 0).rotateZ(-0.2)
	var tallBox = B2.box(5, 5, 10).translate(0,-1,1).flipped()
	tallBox = B2.extrudeVertices([V3(0,0),V3(1, 3), V3(3, 1)], P3.XY.flipped(), V3(0, 0, 10), 'lol')
		.translate(1,0,-1)
		//.flipped()
	/*
	bMesh.computeNormalLines(1)
	bMesh.compile()
	*/
	var a = B2.box(2, 20, 2)
	.translate(0, -10, 2)
	//.rotateX(0.1)
	.flipped()
	var curve = new EllipseCurve(V3.ZERO, V3(3, 0, 0), V3(0, 6, 0))
	a = B2.extrudeEdges([
		new B2.PCurveEdge(curve, V3(-3, 0, 0), curve.at(-3), PI, -3, null, V3.Y, curve.tangentAt(-3).negated()),
		StraightEdge.throughPoints(curve.at(-3), V3.ZERO),
		StraightEdge.throughPoints(V3.ZERO, V3(-3, 0, 0))], P3.XY.flipped(), V3(0, 0, 6))
	.translate(2, 2, 2)
	.flipped()
	//a = new B2([a.faces[4]])
	//b = new B2([new B2.PlaneFace(new PlaneSurface(P3.XY), verticesChain([V3(0, 0, 0), V3(8, 0, 0), V3(0, 8, 0)]),
	//	[verticesChain([V3(1, 1, 0), V3(1, 5, 0), V3(5, 1, 0)])])])
	//var b = B2.cylinder(4, 12)
	var b = B2.rotateEdges(verticesChain([V3(0, 0, 0), V3(6, 0, 0), V3(6, 0, 5), V3(0, 0, 5)]), 2*PI, 3)
		//.rotateX(PI/4)
		//.rotateZ(0.1)
		//.translate(7.2, 4,3.9)
		//.flipped()
	//b = new B2([b.faces[1]])
	a  = wideBox
	b = tallBox.flipped()
	console.log(b.ss)
	bMesh = b.toMesh()
	//var b = B2.puckman(10, PI/2 -0.1, 12, 'puckman').translate(-4, -1, 1)

	aMesh = a.toMesh()
	var c
	paintScreen2()
	c = a.intersection(b, true, false)
	//c = new B2([c.faces[0]])
	cMesh = c && c.toMesh()
	//dMesh = a.minus(b).toMesh()
	console.log(c.toSource())
	paintScreen2()
	//dMesh = wideBox.intersection(tallBox, false, true).flipped().toMesh()
	/*
	 var cyl = new CylinderSurface(L3.Z, 5)
	 bMesh = cyl.toMesh(0, 50)
	 var curve = cyl.getIntersectionWithPlane(cpTop)
	 bMesh.addVertexBuffer('edgeTangents', 'edgeTangents')
	 for (var t = 0; t < 2 * PI; t+=0.1) {
	 var p = curve.at(t);
	 bMesh.edgeTangents.push(p, p.plus(curve.tangentAt(t).toLength(1)))
	 bMesh.edgeTangents.push(p, p.plus(curve.normalAt(t).toLength(1)))
	 }
	 console.log(bMesh.edgeTangents.map(V3.ss))
	 bMesh.compile()
	 //disableConsole()
	 */
	//c = new B2(curveface)
	//aMesh = cMesh = null
	//cMesh.computeNormalLines(0.2);cMesh.compile()
	//cMesh = c.toMesh()
}






var aMesh, bMesh, cMesh, dMesh
function paintScreen2() {
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
	gl.loadIdentity();
	gl.scale(10, 10, 10);

	gl.loadIdentity();

	//drawVectors()

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
		lightingShader.uniforms({ color: rgbToVec4(COLORS.TS_FILL),
			camPos: eyePos }).draw(bMesh);
		bMesh.edgeTangents && singleColorShader.uniforms({ color: rgbToVec4(COLORS.TS_STROKE) })
			.drawBuffers({gl_Vertex: bMesh.vertexBuffers.edgeTangents}, null, gl.LINES)
		bMesh.edgeTangents2 && singleColorShader.uniforms({ color: rgbToVec4(COLORS.RD_STROKE) })
			.drawBuffers({gl_Vertex: bMesh.vertexBuffers.edgeTangents2}, null, gl.LINES)
		gl.popMatrix()
	}
	if (cMesh) {
		gl.pushMatrix()
		gl.translate(30, 0, 0)
		gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
		cMesh.lines && singleColorShader.uniforms({ color: rgbToVec4(COLORS.TS_STROKE) }).draw(cMesh, 'LINES');
		gl.projectionMatrix.m[11] += 1 / (1 << 20)
		lightingShader.uniforms({ color: rgbToVec4(COLORS.RD_FILL),
			camPos: eyePos }).draw(cMesh)

		gl.popMatrix()
	}
	if (dMesh) {
		gl.pushMatrix()
		gl.translate(45, 0, 0)
		//gl.scale(10, 10, 10)
		gl.projectionMatrix.m[11] -= 1 / (1 << 20) // prevent Z-fighting
		dMesh.lines && singleColorShader.uniforms({ color: rgbToVec4(COLORS.RD_STROKE) }).draw(dMesh, 'LINES');
		gl.projectionMatrix.m[11] += 1 / (1 << 20)
		lightingShader.uniforms({ color: rgbToVec4(0xffFF00),
			camPos: eyePos }).draw(dMesh)
		dMesh.curve1 && singleColorShader.uniforms({ color: rgbToVec4(COLORS.TS_STROKE) })
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
		//gl.scale(0.5,0.5,0.5)
		lightingShader.uniforms({color: rgbToVec4(NLA.randomColor())}).draw(sMesh)
		lightingShader.uniforms({color: rgbToVec4(NLA.randomColor())}).draw(sMesh)
		singleColorShader.uniforms({ color: rgbToVec4(COLORS.RD_STROKE) }).draw(sMesh, 'LINES')
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
				paintScreen2();
			}
			if (e.buttons & 2) {
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
	vectorMesh = rotationMesh([V3.ZERO, V3(0, 0.05, 0), V3(0.8, 0.05), V3(0.8, 0.1), V3(1, 0)], L3.X, Math.PI * 2, 8, false)
	sMesh = GL.Mesh.sphere2(2)

	singleColorShader = new GL.Shader(vertexShaderBasic, fragmentShaderColor);
	singleColorShaderHighlight = new GL.Shader(vertexShaderBasic, fragmentShaderColorHighlight);
	textureColorShader = new GL.Shader(vertexShaderTextureColor, fragmentShaderTextureColor);
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

}