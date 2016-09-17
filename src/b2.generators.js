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
		V3.create(0, 0, 0),
		V3.create(0, h, 0),
		V3.create(w, h, 0),
		V3.create(w, 0, 0)
	]
	return B2.extrudeVertices(baseVertices, P3.XY.flipped(), V3.create(0, 0, d), name)
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
	var circleCurve = new EllipseCurve(V3.ZERO, V3.create(radius, 0, 0), V3.create(0, -radius, 0))
	var a = circleCurve.at(0)
	var b = circleCurve.at(-rads)
	var edges = [
		StraightEdge.throughPoints(a, V3.ZERO),
		StraightEdge.throughPoints(V3.ZERO, b),
		new B2.PCurveEdge(circleCurve, b, a, -rads, 0, null, circleCurve.tangentAt(-rads), circleCurve.tangentAt(0))]
	return B2.extrudeEdges(edges, P3.XY.flipped(), V3.create(0, 0, height), name)
}

/**
 *
 * @param {Edge[]} baseFaceEdges
 * @param {P3} baseFacePlane
 * @param {V3} offset
 * @param {string} name
 * @param {string=} source
 * @returns {B2}
 */
B2.extrudeEdges = function (baseFaceEdges, baseFacePlane, offset, name, source) {
	Array.from(NLA.combinations(baseFaceEdges.length)).forEach(({i, j}) => {
		assertf(() => !B2.Edge.edgesIntersect(baseFaceEdges[i], baseFaceEdges[j]), baseFaceEdges[i].sce+baseFaceEdges[j].sce)
	})
	assertf(() => B2.Edge.isLoop(baseFaceEdges))
	// TODO checks..
	if (offset.dot(baseFacePlane.normal) > 0) {
		baseFacePlane = baseFacePlane.flipped()
	}
	let basePlaneSurface = new PlaneSurface(baseFacePlane)
	assert(basePlaneSurface.edgeLoopCCW(baseFaceEdges), "edges not CCW on baseFacePlane")
	var translationMatrix = M4.translation(offset)
	var topEdges = baseFaceEdges.map(edge => edge.transform(translationMatrix))
	var edgeCount = baseFaceEdges.length
	var bottomFace = new B2.PlaneFace(basePlaneSurface, baseFaceEdges)
	var topFaceEdges = topEdges.map(edge => edge.flipped()).reverse()
	var topFace = new B2.PlaneFace(new PlaneSurface(baseFacePlane.flipped().translated(offset)), topFaceEdges)

	var ribs = NLA.arrayFromFunction(edgeCount,
		i => StraightEdge.throughPoints(baseFaceEdges[i].a, topEdges[i].a, name + 'Rib' + i))

	var faces = baseFaceEdges.map((edge, i) => {
		var j = (i + 1) % edgeCount
		var faceEdges = [baseFaceEdges[i].flipped(), ribs[i], topEdges[i], ribs[j].flipped()]
		var curve = edge.curve
		if (edge instanceof StraightEdge) {
			var surfaceNormal = offset.cross(edge.tangent).normalized()
			let surface = new PlaneSurface(P3.normalOnAnchor(surfaceNormal, edge.a))
			return new B2.PlaneFace(surface, faceEdges)
		} else if (curve instanceof EllipseCurve) {
			let surface = new CylinderSurface(curve, offset.normalized(), -1)
			return new B2.RotationFace(surface, faceEdges)
		} else if (curve instanceof BezierCurve) {
			let curveDir = edge.reversed ? offset : offset.negated()
			let surface = new ProjectedCurveSurface(curve, curveDir.normalized(), 0, 1)
			return new B2.RotationFace(surface, faceEdges)
		} else {
			assert(false, edge)
		}
	})
	faces.push(bottomFace, topFace)
	source = source || `B2.extrudeEdges(${baseFaceEdges.sce}, ${baseFacePlane.sce}, ${offset.sce}, "${name}")` // todo proper escaping
	return new B2(faces, false, source)
}
B2.cylinder = function (radius, height, rads) {
	return B2.rotateEdges(verticesChain([V3(0, 0, 0), V3(radius, 0, 0), V3(radius, 0, height), V3(0, 0, height)]), rads || 2 * PI)
}
B2.torus = function (rSmall, rLarge, rads) {
	assertNumbers(rSmall, rLarge, rads)
	assertf(() => rLarge > rSmall)
	return B2.rotateEdges([EllipseCurve.circle(rSmall, V3.create(rLarge, 0, 0))], rads)
}

/**
 *
 * @param {Edge[]} edges
 * @param {number} rads
 * @returns {B2}
 */
B2.rotateEdges = function (edges, rads) {
	var rotationMatrix = M4.rotationZ(rads)
	var open = !NLA.equals(rads, 2 * PI)
	var endEdges = open ? edges.map(edge => edge.transform(rotationMatrix)) : edges
	var edgeCount = edges.length;
	var ribs = NLA.arrayFromFunction(edgeCount, i => {
		var a = edges[i].a, radius = a.lengthXY()
		var b = endEdges[i].a
		if (!NLA.isZero(radius)) {
			var curve = new EllipseCurve(V3.create(0, 0, a.z), V3.create(-radius, 0, 0), V3.create(0, -radius, 0))
			var aT = -PI, bT = -PI + rads
			return new B2.PCurveEdge(curve, a, b, aT, bT, null, curve.tangentAt(aT), curve.tangentAt(bT))
		}
	})
	var faces = edges.map((edge, i) => {
		var ipp = (i + 1) % edgeCount
		var faceEdges = [
			edge.flipped(),
			!NLA.isZero(edge.a.x) && ribs[i],
			endEdges[i],
			!NLA.isZero(edge.b.x) && ribs[ipp].flipped()].filter(x => x)
		var surface
		var curve = edge.curve;
		if (edge instanceof StraightEdge) {
			var line = edge.curve
			if (line.dir1.isParallelTo(V3.Z)) {
				if (NLA.isZero(edge.a.x)) { return }
				let flipped = edge.a.z > edge.b.z
				let surface = new CylinderSurface(ribs[i].curve, V3.Z, !flipped ? 1 : -1)
				return new B2.RotationFace(surface, faceEdges)
			} else if (line.dir1.isPerpendicularTo(V3.Z)) {
				let flipped = edge.a.x > edge.b.x
				let surface = new PlaneSurface(P3(V3.Z, edge.a.z))
				if (!flipped) surface = surface.flipped()
				if (!open) {
					var hole = flipped
						? !NLA.isZero(edge.b.x) && ribs[ipp].flipped()
						: !NLA.isZero(edge.a.x) && ribs[i]
					return new B2.PlaneFace(surface, [flipped ? ribs[i] : ribs[ipp].flipped()], hole && [[hole]])
				}
				return new B2.PlaneFace(surface, faceEdges)
			} else {
				// apex is intersection of segment with Z-axis
				let a = edge.a, b = edge.b
				let apexZ = a.z - a.x * (b.z - a.z) / (b.x - a.x)
				let apex = V3(0, 0, apexZ)
				let flipped = edge.a.z > edge.b.z
				let surface = ConicSurface.atApexThroughEllipse(apex, ribs[a.x > b.x ? i : ipp].curve, !flipped ? 1 : -1)
				return new B2.RotationFace(surface, faceEdges)
			}
		} if (edge.curve instanceof EllipseCurve) {
			let flipped = undefined
			let /** @type EllipseCurve */ ell = edge.curve.rightAngled()
			let f1Perp = ell.f1.isPerpendicularTo(V3.Z), f2Perp = ell.f2.isPerpendicularTo(V3.Z)
			if (L3.Z.containsPoint(ell.center) && (f1Perp || f2Perp)) {
				let f3length = f1Perp ? ell.f1.length() : ell.f2.length()
				if (flipped) { f3length *= -1 }
				let surface = new EllipsoidSurface(ell.center, ell.f1, ell.f2, ell.f1.cross(ell.f2).toLength(f3length))
				return new B2.RotationFace(surface, faceEdges)
			}
		} else {
			assert (false, edge)
		}
	}).filter(x =>x)
	if (open) {
		var endFaceEdges = endEdges.map(edge => edge.flipped()).reverse()
		var endFace = new B2.PlaneFace(new PlaneSurface(P3.ZX.rotateZ(rads)), endFaceEdges)
		faces.push(new B2.PlaneFace(new PlaneSurface(P3.ZX.flipped()), edges), endFace)
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
			surface = new PlaneSurface(flipped ? P3(V3.Z, edge.a.z) : P3(V3.Z.negated(), -edge.a.z))
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
				face = new B2.PlaneFace(surface, newEdges)
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
				face = new B2.PlaneFace(surface, contour, hole ? [hole] : [])
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
				faces.push(new B2.PlaneFace(surface, faceEdges))
			} else {
				assert(false, edge.toString())
			}
		}
	})
	if (open) {
		var endFaceEdges = ribs[count].map(edge => edge.flipped()).reverse()
		var endFace = new B2.PlaneFace(new PlaneSurface(P3.XZ.rotateZ(rads)), endFaceEdges)
		faces.push(new B2.PlaneFace(new PlaneSurface(P3.XZ.flipped()), edges), endFace)
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
		bottom = B2.PlaneFace.forVertices(new PlaneSurface(baseFacePlane), baseVertices),
		top = B2.PlaneFace.forVertices(new PlaneSurface(baseFacePlane.flipped().translated(offset)), topVertices)]
	var m = baseVertices.length
	var ribs = NLA.arrayFromFunction(m, i => StraightEdge.throughPoints(baseVertices[i], topVertices[m - 1 - i]))
	for (var i = 0; i < m; i++) {
		var j = (i + 1) % m
		faces.push(
			new B2.PlaneFace(
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
		new B2.PlaneFace(PlaneSurface.throughPoints(a, b, c), [ab, bc, ac.flipped()]),
		new B2.PlaneFace(PlaneSurface.throughPoints(a, d, b), [ad, bd.flipped(), ab.flipped()]),
		new B2.PlaneFace(PlaneSurface.throughPoints(b, d, c), [bd, cd.flipped(), bc.flipped()]),
		new B2.PlaneFace(PlaneSurface.throughPoints(c, d, a), [cd, ad.flipped(), ac])
	]
	return new B2(faces)
}