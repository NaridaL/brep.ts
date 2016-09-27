/**
 * Created by aval on 16/02/2016.
 */

/**
 * Create a list of StraightEdges from a list of vertices.
 * @param {Array.<V3>} vertices
 * @param {boolean=} closed Whether to connect the first and last vertices. Defaults to true.
 * @returns {StraightEdge[]}
 */
function verticesChain(vertices, closed) {
	closed = 'boolean' === typeof closed ? closed : true
	var vc = vertices.length
	return NLA.arrayFromFunction(closed ? vc : vc - 1,
		i => StraightEdge.throughPoints(vertices[i], vertices[(i + 1) % vc]))
}
B2.box = function (w, h, d, name) {
	assertNumbers(w, h, d)
	assertInst('string' === typeof name)
	var baseVertices = [
		V(0, 0, 0),
		V(0, h, 0),
		V(w, h, 0),
		V(w, 0, 0)
	]
	return B2.extrudeVertices(baseVertices, P3.XY.flipped(), V(0, 0, d), name)
}

/**
 *
 * @param {number} radius
 * @param {number} rads
 * @param {number} height
 * @param {string} name
 * @returns {B2}
 */
B2.puckman = function (radius, rads, height, name) {
	// TODO: argument checking
	var circleCurve = new EllipseCurve(V3.ZERO, V(radius, 0, 0), V(0, -radius, 0))
	var a = circleCurve.at(0)
	var b = circleCurve.at(-rads)
	var edges = [
		StraightEdge.throughPoints(a, V3.ZERO),
		StraightEdge.throughPoints(V3.ZERO, b),
		new PCurveEdge(circleCurve, b, a, -rads, 0, null, circleCurve.tangentAt(-rads), circleCurve.tangentAt(0))]
	return B2.extrudeEdges(edges, P3.XY.flipped(), V(0, 0, height), name)
}

B2.registerVertexName = function (map, name, p) {
	// TODO
	if (!Array.from(map.keys()).some(p2 => p2.like(p))) {
		map.set(p, name)
	}
}

/**
 *
 * @param {Array.<Edge>} baseFaceEdges
 * @param {P3} baseFacePlane
 * @param {V3} offset
 * @param {string} name
 * @param {string=} gen
 * @returns {B2}
 */
B2.extrudeEdges = function (baseFaceEdges, baseFacePlane, offset, name, gen) {
	Array.from(NLA.combinations(baseFaceEdges.length)).forEach(({i, j}) => {
		assertf(() => !Edge.edgesIntersect(baseFaceEdges[i], baseFaceEdges[j]), baseFaceEdges[i].sce+baseFaceEdges[j].sce)
	})
	assertf(() => Edge.isLoop(baseFaceEdges))
	// TODO checks..
	if (offset.dot(baseFacePlane.normal) > 0) {
		baseFacePlane = baseFacePlane.flipped()
	}
	let vertexNames = new Map()
	let basePlaneSurface = new PlaneSurface(baseFacePlane)
	assert(basePlaneSurface.edgeLoopCCW(baseFaceEdges), "edges not CCW on baseFacePlane")
	var translationMatrix = M4.translation(offset)
	var topEdges = baseFaceEdges.map(edge => edge.transform(translationMatrix, 'top'))
	var edgeCount = baseFaceEdges.length
	var bottomFace = new PlaneFace(basePlaneSurface, baseFaceEdges, [], name + 'Bottom')
	var topFaceEdges = topEdges.map(edge => edge.flipped()).reverse()
	var topFace = new PlaneFace(new PlaneSurface(baseFacePlane.flipped().translated(offset)), topFaceEdges, [], name + 'Top')


	baseFaceEdges.forEach(edge => B2.registerVertexName(vertexNames, edge.name + 'A', edge.a))
	topFaceEdges.forEach(edge => B2.registerVertexName(vertexNames, edge.name + 'A', edge.a))

	var ribs = NLA.arrayFromFunction(edgeCount,
		i => StraightEdge.throughPoints(baseFaceEdges[i].a, topEdges[i].a, name + 'Rib' + i))

	let faces = baseFaceEdges.map((edge, i) => {
		let faceName = name + 'Wall' + i
		let j = (i + 1) % edgeCount
		let faceEdges = [baseFaceEdges[i].flipped(), ribs[i], topEdges[i], ribs[j].flipped()]
		let curve = edge.curve
		let surface = projectCurve(curve, offset, edge.reversed)
		if (edge instanceof StraightEdge) {
			return new PlaneFace(surface, faceEdges, undefined, faceName)
		} else if (curve instanceof EllipseCurve) {
			return new RotationFace(surface, faceEdges, undefined, faceName)
		} else if (curve instanceof BezierCurve) {
			return new RotationFace(surface, faceEdges, undefined, faceName)
		} else {
			assert(false, edge)
		}
	})
	faces.push(bottomFace, topFace)
	gen = gen || `B2.extrudeEdges(${baseFaceEdges.sce}, ${baseFacePlane.sce}, ${offset.sce}, ${JSON.stringify(name)})`
	return new B2(faces, false, gen, vertexNames)
}

function projectCurve(curve, offset, flipped) {
	if (curve instanceof L3) {
		var surfaceNormal = offset.cross(curve.dir1).toLength(flipped ? -1 : 1)
		return new PlaneSurface(P3.normalOnAnchor(surfaceNormal, curve.anchor))
	}
	if (curve instanceof EllipseCurve) {
		let curveDir = flipped ? offset : offset.negated()
		return new CylinderSurface(curve, curveDir.normalized())
	}
	if (curve instanceof BezierCurve) {
		let curveDir = flipped ? offset : offset.negated()
		return new ProjectedCurveSurface(curve, curveDir.normalized(), 0, 1)
	}
	assert(false, curve)
}

B2.cylinder = function (radius, height, rads) {
	return B2.rotateEdges(verticesChain([V3(0, 0, 0), V3(radius, 0, 0), V3(radius, 0, height), V3(0, 0, height)]), rads || 2 * PI)
}
B2.torus = function (rSmall, rLarge, rads) {
	assertNumbers(rSmall, rLarge, rads)
	assertf(() => rLarge > rSmall)
	return B2.rotateEdges([EllipseCurve.circle(rSmall, V(rLarge, 0, 0))], rads)
}

/**
 *
 * @param {Edge[]} edges
 * @param {number} rads
 * @returns {B2}
 */
B2.rotateEdges = function (edges, rads, name) {
	var rotationMatrix = M4.rotationZ(rads)
	var open = !NLA.equals(rads, 2 * PI)
	var endEdges = open ? edges.map(edge => edge.transform(rotationMatrix)) : edges
	var edgeCount = edges.length;
	var ribs = NLA.arrayFromFunction(edgeCount, i => {
		var a = edges[i].a, radius = a.lengthXY()
		var b = endEdges[i].a
		if (!NLA.isZero(radius)) {
			var curve = new EllipseCurve(V(0, 0, a.z), V(-radius, 0, 0), V(0, -radius, 0))
			var aT = -PI, bT = -PI + rads
			return new PCurveEdge(curve, a, b, aT, bT, null, curve.tangentAt(aT), curve.tangentAt(bT), name + 'rib' + i)
		}
	})
	var faces = edges.map((edge, i) => {
		var ipp = (i + 1) % edgeCount
		var faceEdges = [
			edge.flipped(),
			!NLA.isZero(edge.a.x) && ribs[i],
			endEdges[i],
			!NLA.isZero(edge.b.x) && ribs[ipp].flipped()].filter(x => x)
		var curve = edge.curve;
		if (edge instanceof StraightEdge) {
			var line = edge.curve
			if (line.dir1.isParallelTo(V3.Z)) {
				if (NLA.isZero(edge.a.x)) { return }
				let flipped = edge.a.z > edge.b.z
				let surface = new CylinderSurface(ribs[i].curve, !flipped ? V3.Z : V3.Z.negated())
				return new RotationFace(surface, faceEdges)
			} else if (line.dir1.isPerpendicularTo(V3.Z)) {
				let flipped = edge.a.x > edge.b.x
				let surface = new PlaneSurface(new P3(V3.Z, edge.a.z))
				if (!flipped) surface = surface.flipped()
				if (!open) {
					var hole = flipped
						? !NLA.isZero(edge.b.x) && ribs[ipp].flipped()
						: !NLA.isZero(edge.a.x) && ribs[i]
					return new PlaneFace(surface, [flipped ? ribs[i] : ribs[ipp].flipped()], hole && [[hole]])
				}
				return new PlaneFace(surface, faceEdges)
			} else {
				// apex is intersection of segment with Z-axis
				let a = edge.a, b = edge.b
				let apexZ = a.z - a.x * (b.z - a.z) / (b.x - a.x)
				let apex = V3(0, 0, apexZ)
				let flipped = edge.a.z > edge.b.z
				let surface = ConicSurface.atApexThroughEllipse(apex, ribs[a.x > b.x ? i : ipp].curve, !flipped ? 1 : -1)
				return new RotationFace(surface, faceEdges)
			}
		} if (edge.curve instanceof EllipseCurve) {
			let flipped = undefined
			let /** @type EllipseCurve */ ell = edge.curve.rightAngled()
			let f1Perp = ell.f1.isPerpendicularTo(V3.Z), f2Perp = ell.f2.isPerpendicularTo(V3.Z)
			if (L3.Z.containsPoint(ell.center) && (f1Perp || f2Perp)) {
				let f3length = f1Perp ? ell.f1.length() : ell.f2.length()
				if (flipped) { f3length *= -1 }
				let surface = new EllipsoidSurface(ell.center, ell.f1, ell.f2, ell.f1.cross(ell.f2).toLength(f3length))
				return new RotationFace(surface, faceEdges)
			}
		} else {
			assert (false, edge)
		}
	}).filter(x =>x)
	if (open) {
		var endFaceEdges = endEdges.map(edge => edge.flipped()).reverse()
		var endFace = new PlaneFace(new PlaneSurface(P3.ZX.rotateZ(rads)), endFaceEdges)
		faces.push(new PlaneFace(new PlaneSurface(P3.ZX.flipped()), edges), endFace)
	}
	return new B2(faces)
}
B2.rotStep = function (edges, rads, count) {
	var radStep = rads / count
	var open = !NLA.equals(rads, 2 * PI)
	var ribCount = !open ? count : count + 1
	var ribs = NLA.arrayFromFunction(ribCount, i => {
		if (i == 0) return edges
		var matrix = M4.rotationZ(radStep * i)
		return edges.map(edge => edge.transform(matrix))
	})
	var hs = NLA.arrayFromFunction(count, i => {
		var ipp = (i + 1) % ribCount
		return NLA.arrayFromFunction(edges.length, j => {
			if (!NLA.isZero(edges[j].a.lengthXY())) {
				return StraightEdge.throughPoints(ribs[i][j].a, ribs[ipp][j].a)
			}
		})
	})
	var faces = [], surface, face
	edges.forEach((edge, i) => {
		var ipp = (i + 1) % edges.length
		if (edge instanceof StraightEdge && edge.curve.dir1.isPerpendicularTo(V3.Z)) {
			var flipped = edge.a.x > edge.b.x;
			surface = new PlaneSurface(flipped ? new P3(V3.Z, edge.a.z) : new P3(V3.Z.negated(), -edge.a.z))
			if (open) {
				var newEdges = []
				if (!NLA.isZero(edge.a.x)) {
					newEdges.pushAll(NLA.arrayFromFunction(count, j => hs[j][i]))
				}
				newEdges.push(ribs[count][i])
				if (!NLA.isZero(edge.b.x)) {
					newEdges.pushAll(NLA.arrayFromFunction(count, j => hs[count - j - 1][ipp].flipped()))
				}
				newEdges.push(edge.flipped())
				face = new PlaneFace(surface, newEdges)
			} else {
				var contour = flipped
					? NLA.arrayFromFunction(count, j => hs[j][i])
					: NLA.arrayFromFunction(count, j => hs[count - j - 1][ipp].flipped())
				var hole
				if (flipped && !NLA.isZero(edge.b.x)) {
					hole = NLA.arrayFromFunction(count, j => hs[count - j - 1][ipp].flipped())
				} else if (!flipped && !NLA.isZero(edge.a.x)) {
					hole = NLA.arrayFromFunction(count, j => hs[j][i])
				}
				face = new PlaneFace(surface, contour, hole ? [hole] : [])
			}
			faces.push(face)
			return
		} else if (edge instanceof StraightEdge) {
			if (NLA.isZero(edge.a.lengthXY()) && NLA.isZero(edge.b.lengthXY())) {
				return
			}
		}
		for (var r = 0; r < count; r++) {
			var rpp = (r + 1) % ribCount
			var faceEdges = [ribs[r][i].flipped(), hs[r][i], ribs[rpp][i], hs[r][ipp] && hs[r][ipp].flipped()].filter(x => x)
			if (edge instanceof StraightEdge) {
				var surface = new PlaneSurface(P3.throughPoints(faceEdges[0].a, faceEdges[1].a, faceEdges[2].a))
				faces.push(new PlaneFace(surface, faceEdges))
			} else {
				assert(false, edge.toString())
			}
		}
	})
	if (open) {
		var endFaceEdges = ribs[count].map(edge => edge.flipped()).reverse()
		var endFace = new PlaneFace(new PlaneSurface(P3.XZ.rotateZ(rads)), endFaceEdges)
		faces.push(new PlaneFace(new PlaneSurface(P3.XZ.flipped()), edges), endFace)
	}
	return new B2(faces)
}

B2.extrudeVertices = function (baseVertices, baseFacePlane, offset, name, source) {
	assert(baseVertices.every(v => v instanceof V3), "baseVertices.every(v => v instanceof V3)")
	assertInst(P3, baseFacePlane)
	assertVectors(offset)
	if (baseFacePlane.normal.dot(offset) > 0) baseFacePlane = baseFacePlane.flipped()
	if (!isCCW(baseVertices, baseFacePlane.normal)) {
		baseVertices = baseVertices.reverse()
	}
	var topVertices = baseVertices.map((v) => v.plus(offset)).reverse()
	//var topPlane = basePlane.translated(offset)
	var top, bottom
	var faces = [
		bottom = PlaneFace.forVertices(new PlaneSurface(baseFacePlane), baseVertices),
		top = PlaneFace.forVertices(new PlaneSurface(baseFacePlane.flipped().translated(offset)), topVertices)]
	var m = baseVertices.length
	var ribs = NLA.arrayFromFunction(m, i => StraightEdge.throughPoints(baseVertices[i], topVertices[m - 1 - i]))
	for (var i = 0; i < m; i++) {
		var j = (i + 1) % m
		faces.push(
			new PlaneFace(
				PlaneSurface.throughPoints(baseVertices[j], baseVertices[i], topVertices[m - j - 1]),
				[bottom.edges[i].flipped(), ribs[i], top.edges[m - j - 1].flipped(), ribs[j].flipped()], [], name + "wall" + i))
	}
	let edges = verticesChain(baseVertices, true)
	source = source || `B2.extrudeVertices(${baseVertices.sce}, ${baseFacePlane.sce}, ${offset.sce}, "${name}")`
	return B2.extrudeEdges(edges, baseFacePlane, offset, name, source)
}

/**
 * Returns a tetrahedron (3 sided pyramid).
 * Faces will face outwards.
 * abcd can be in any order. The only constraint is that abcd cannot be on a common plane.
 *
 * @param {V3} a
 * @param {V3} b
 * @param {V3} c
 * @param {V3} d
 * @returns {B2}
 */
B2.tetrahedron = function (a, b, c, d) {
	assertVectors(a, b, c, d)
	let dDistance = P3.throughPoints(a, b, c).distanceToPointSigned(d)
	if (NLA.isZero(dDistance)) {
		throw new Error("four points are coplanar")
	}
	if (dDistance > 0) {
		[c, d] = [d, c]
	}
	let ab = StraightEdge.throughPoints(a, b)
	let ac = StraightEdge.throughPoints(a, c)
	let ad = StraightEdge.throughPoints(a, d)
	let bc = StraightEdge.throughPoints(b, c)
	let bd = StraightEdge.throughPoints(b, d)
	let cd = StraightEdge.throughPoints(c, d)
	let faces = [
		new PlaneFace(PlaneSurface.throughPoints(a, b, c), [ab, bc, ac.flipped()]),
		new PlaneFace(PlaneSurface.throughPoints(a, d, b), [ad, bd.flipped(), ab.flipped()]),
		new PlaneFace(PlaneSurface.throughPoints(b, d, c), [bd, cd.flipped(), bc.flipped()]),
		new PlaneFace(PlaneSurface.throughPoints(c, d, a), [cd, ad.flipped(), ac])
	]
	let gen = `B2.tetrahedron(${a.sce}, ${b.sce}, ${c.sce}, ${d.sce})`
	return new B2(faces, false, gen)
}