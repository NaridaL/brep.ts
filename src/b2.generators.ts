function projectCurve(curve: Curve, offset: V3, flipped: boolean): Surface {
	if (curve instanceof L3) {
		let surfaceNormal = offset.cross(curve.dir1).toLength(flipped ? -1 : 1)
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
	assertNever()
}
function rotateCurve(curve: Curve, offset: V3, flipped: boolean): Surface {
	if (curve instanceof L3) {
		let surface
		if (curve.dir1.isParallelTo(V3.Z)) {
			if (NLA.eq0(line.anchor.x)) {
				return
			}
			let flipped = line.anchor.z > edge.b.z
			surface = new CylinderSurface(ribs[i].curve, !flipped ? V3.Z : V3.Z.negated())
		} else if (curve.dir1.isPerpendicularTo(V3.Z)) {
			let flipped = line.anchor.x > edge.b.x
			let surface = new PlaneSurface(new P3(V3.Z, line.anchor.z))
			if (!flipped) surface = surface.flipped()
			if (!open) {
				const hole = flipped
					? !NLA.eq0(edge.b.x) && ribs[ipp].flipped()
					: !NLA.eq0(line.anchor.x) && ribs[i]
				return new PlaneFace(surface, [flipped ? ribs[i] : ribs[ipp].flipped()], hole && [[hole]])
			}
			return new PlaneFace(surface, faceEdges)
		} else {
			// apex is intersection of segment with Z-axis
			let a = line.anchor, b = edge.b
			let apexZ = a.z - a.x * (b.z - a.z) / (b.x - a.x)
			let apex = new V3(0, 0, apexZ)
			let flipped = line.anchor.z > edge.b.z
			surface = ConicSurface.atApexThroughEllipse(apex, ribs[a.x > b.x ? i : ipp].curve as EllipseCurve, !flipped ? 1 : -1)
		}
		return Face.create(surface, faceEdges)
	}
	if (edge.curve instanceof EllipseCurve) {
		let flipped = undefined
		let ell = edge.curve.rightAngled()
		let f1Perp = ell.f1.isPerpendicularTo(V3.Z), f2Perp = ell.f2.isPerpendicularTo(V3.Z)
		if (L3.Z.containsPoint(ell.center) && (f1Perp || f2Perp)) {
			let f3length = f1Perp ? ell.f1.length() : ell.f2.length()
			if (flipped) {
				f3length *= -1
			}
			let surface = new EllipsoidSurface(ell.center, ell.f1, ell.f2, ell.f1.cross(ell.f2).toLength(f3length))
			return new RotationFace(surface, faceEdges)
		}
	} else {
		assert(false, edge)
	}
	if (curve instanceof L3) {
		let surfaceNormal = offset.cross(curve.dir1).toLength(flipped ? -1 : 1)
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
	assertNever()
}

namespace B2T {

    export function box(w: number = 1, h: number = 1, d: number = 1, name?: string): B2 {
		assertNumbers(w, h, d)
		assertInst('string' === typeof name)
		const baseVertices = [
			new V3(0, 0, 0),
			new V3(0, h, 0),
			new V3(w, h, 0),
			new V3(w, 0, 0)
		]
		return B2T.extrudeVertices(baseVertices, P3.XY.flipped(), new V3(0, 0, d), name, `B2T.box(${w}, ${h}, ${d}, "${name || ''}")`)
	}

	export function puckman(radius: number, rads: number, height: number, name: string): B2 {
		assertf(() => NLA.lt(0, radius))
		assertf(() => NLA.lt(0, rads) && NLA.le(rads, TAU))
		assertf(() => NLA.lt(0, height))
		const edges = StraightEdge.chain([V3.ZERO, new V3(radius, 0, 0), new V3(radius, 0, height), new V3(0, 0, height)], true)
		return B2T.rotateEdges(edges, rads, name || 'puckman' + globalId++)
	}

	export function registerVertexName(map, name, p) {
		// TODO
		if (!Array.from(map.keys()).some(p2 => p2.like(p))) {
			map.set(p, name)
		}
	}

	export function extrudeEdges(baseFaceEdges: Edge[], baseFacePlane: P3, offset: V3, name: string, gen?: string): B2 {
		Array.from(NLA.combinations(baseFaceEdges.length)).forEach(({i, j}) => {
			assertf(() => !Edge.edgesIntersect(baseFaceEdges[i], baseFaceEdges[j]), baseFaceEdges[i].sce + baseFaceEdges[j].sce)
		})
		assertf(() => Edge.isLoop(baseFaceEdges))
		// TODO checks..
		if (offset.dot(baseFacePlane.normal) > 0) {
			baseFacePlane = baseFacePlane.flipped()
		}
		let vertexNames = new Map()
		let basePlaneSurface = new PlaneSurface(baseFacePlane)
		assert(basePlaneSurface.edgeLoopCCW(baseFaceEdges), "edges not CCW on baseFacePlane")
		const translationMatrix = M4.translation(offset)
		const topEdges = baseFaceEdges.map(edge => edge.transform(translationMatrix, 'top'))
		const edgeCount = baseFaceEdges.length
		const bottomFace = new PlaneFace(basePlaneSurface, baseFaceEdges, [], name + 'Bottom')
		const topFaceEdges = topEdges.map(edge => edge.flipped()).reverse()
		const topFace = new PlaneFace(new PlaneSurface(baseFacePlane.flipped().translated(offset)), topFaceEdges, [], name + 'Top')


		baseFaceEdges.forEach(edge => B2T.registerVertexName(vertexNames, edge.name + 'A', edge.a))
		topFaceEdges.forEach(edge => B2T.registerVertexName(vertexNames, edge.name + 'A', edge.a))

		const ribs = NLA.arrayFromFunction(edgeCount,
			i => StraightEdge.throughPoints(baseFaceEdges[i].a, topEdges[i].a, name + 'Rib' + i))

		const faces = baseFaceEdges.map((edge, i) => {
			const faceName = name + 'Wall' + i
			const j = (i + 1) % edgeCount
			const faceEdges = [baseFaceEdges[i].flipped(), ribs[i], topEdges[i], ribs[j].flipped()]
			const surface = projectCurve(edge.curve, offset, edge.reversed)
			return Face.create(surface, faceEdges, undefined, faceName)
		})
		faces.push(bottomFace, topFace)
		gen = gen || `B2T.extrudeEdges(${baseFaceEdges.sce}, ${baseFacePlane.sce}, ${offset.sce}, ${JSON.stringify(name)})`
		return new B2(faces, false, gen, vertexNames)
	}


	export function cylinder(radius: number, height: number, rads: number, name: string): B2 {
		const vertices = [new V3(0, 0, 0), new V3(radius, 0, 0), new V3(radius, 0, height), new V3(0, 0, height)]
		return B2T.rotateEdges(StraightEdge.chain(vertices, true), rads || 2 * PI, name)
	}

	export function sphere(radius: number, name: string = 'sphere' + globalId++): B2 {
        const ee = PCurveEdge.forCurveAndTs(new EllipseCurve(V3.ZERO, new V3(0, 0, radius), new V3(radius, 0, 0)), -PI, 0)
        return rotateEdges([StraightEdge.throughPoints(ee.b, ee.a), ee], TAU, name)
    }

	export function torus(rSmall: number, rLarge: number, rads: number, name: string): B2 {
		assertNumbers(rSmall, rLarge, rads)
		assertf(() => rLarge > rSmall)
		const curve = EllipseCurve.circle(rSmall, new V3(rLarge, 0, 0))
		const baseEdges = [PCurveEdge.forCurveAndTs(curve, -Math.PI, 0), PCurveEdge.forCurveAndTs(curve, 0, Math.PI)]
		return B2T.rotateEdges(baseEdges, rads, name || 'torus' + globalId++)
	}
	export function torusUnsplit(rSmall: number, rLarge: number, rads: number, name: string): B2 {
		assertNumbers(rSmall, rLarge, rads)
		assertf(() => rLarge > rSmall)
		let baseEdge = PCurveEdge.forCurveAndTs(EllipseCurve.circle(rSmall, new V3(rLarge, 0, 0)), -Math.PI, Math.PI)
		return B2T.rotateEdges([baseEdge], rads, name || 'torus' + globalId++)
	}

	/**
	 * baseLoop should be CCW on XZ plane for a bounded B2
	 */
	export function rotateEdges(baseLoop: Edge[], totalRads: number, name: string): B2 {
		assert(!NLA.eq(PI, totalRads) || PI == totalRads) // URHGJ
		assertf(() => NLA.lt(0, totalRads) && NLA.le(totalRads, TAU))
		assertf(() => Edge.isLoop(baseLoop))
		const rotationSteps = ceil((totalRads - NLA_PRECISION) / PI)
		const angles = rotationSteps == 1 ? [-PI, -PI + totalRads] : [-PI, 0, totalRads - PI]
		console.log("angles", angles)
		const open = !NLA.eq(totalRads, 2 * PI)
		const ribCurves = baseLoop.map(edge =>  {
			const a = edge.a, radius = a.lengthXY()
			if (!NLA.eq0(radius)) {
				return new EllipseCurve(V(0, 0, a.z), V(-radius, 0, 0), V(0, -radius, 0))
			}
		})
		const surfaces = baseLoop.map((edge, i) => {
			const ipp = (i + 1) % baseLoop.length
			if (edge instanceof StraightEdge) {
				const line = edge.curve
				if (line.dir1.isParallelTo(V3.Z)) {
					if (NLA.eq0(edge.a.x)) {
						return
					}
					let flipped = edge.a.z > edge.b.z
					return new CylinderSurface(ribCurves[i], !flipped ? V3.Z : V3.Z.negated())
				} else if (line.dir1.isPerpendicularTo(V3.Z)) {
					let flipped = edge.a.x > edge.b.x
					let surface = new PlaneSurface(new P3(V3.Z, edge.a.z))
					if (!flipped) surface = surface.flipped()
					return surface
				} else {
					// apex is intersection of segment with Z-axis
					let a = edge.a, b = edge.b
					let apexZ = a.z - a.x * (b.z - a.z) / (b.x - a.x)
					let apex = new V3(0, 0, apexZ)
					let flipped = edge.a.z > edge.b.z
					return ConicSurface.atApexThroughEllipse(apex, ribCurves[a.x > b.x ? i : ipp] as EllipseCurve, !flipped ? 1 : -1)
				}
			}
			/*
			    at(t) = f1 * cos t + f2 sin t
			    rotated projection
			    at2(t) = V(at(t).lengthXY(), 0, at(t).z)
			    at2(t).x = sqrt((f1x cos t + f2x sin t)² + (f1y cos t + f2y sin t)²)
			    at2(t).x = sqrt((f1x² + f1y²) cos² t + (f1x f2x + f1y f2y) cos t sin t + (f2x² + f2y²)sin²t)
			    at2(t).x = sqrt((a² + b²) cos² t + (a c + b d) cos t sin t + (c² + d²)sin²t)
			    (x cos t + y sin t)² = x² cos² t + x y cos t sin t + y² sin² t
			 */
			if (edge.curve instanceof EllipseCurve) {
				let flipped = edge.a.z > edge.b.z
				let ell = edge.curve.rightAngled()
                assert(ell.normal.isPerpendicularTo(V3.Z))
                assert(L3.Z.containsPoint(ell.center))
                let width = ell.f1.length(), height = ell.f2.length()
                if (!ell.isCircular()) {
				    assert(ell.f1.isParallelTo(V3.Z) && ell.f2.isParallelTo(V3.X)
                        || ell.f2.isParallelTo(V3.Z) && ell.f1.isParallelTo(V3.X))
                    if (ell.f1.isParallelTo(V3.Z)) {
				        [width, height] = [height, width]
                    }
                }
                return EllipsoidSurface.forABC(width, (!flipped ? 1 : -1) * width, height, ell.center)
			} else {
				assert(false, edge)
			}
		})
		let stepStartEdges = baseLoop, stepEndEdges
		const faces = []
		for (let rotStep = 0; rotStep < rotationSteps; rotStep++) {
			const aT = angles[rotStep], bT = angles[rotStep + 1]
			const rotationMatrix = M4.rotationZ(angles[rotStep + 1] - PI)
			stepEndEdges = NLA.eq(TAU, bT) ? baseLoop : baseLoop.map(edge => edge.transform(rotationMatrix))
			const ribs = NLA.arrayFromFunction(baseLoop.length, i => {
				const a = stepStartEdges[i].a, radius = a.lengthXY()
				const b = stepEndEdges[i].a
				if (!NLA.eq0(radius)) {
					const curve = ribCurves[i]
					return new PCurveEdge(curve, a, b, aT, bT, null, curve.tangentAt(aT), curve.tangentAt(bT), name + 'rib' + i)
				}
			})
			for (let edgeIndex = 0; edgeIndex < baseLoop.length; edgeIndex++) {
				if (surfaces[edgeIndex]) {
					const edge = stepStartEdges[edgeIndex]
					const ipp = (edgeIndex + 1) % baseLoop.length
					const faceEdges = [
						stepStartEdges[edgeIndex].flipped(),
						!NLA.eq0(edge.a.x) && ribs[edgeIndex],
						stepEndEdges[edgeIndex],
						!NLA.eq0(edge.b.x) && ribs[ipp].flipped()].filter(x => x)
					console.log("ljl", edgeIndex, ipp, ribs, surfaces[rotStep], faceEdges)
					faces.push(Face.create(surfaces[edgeIndex], faceEdges))
				}
			}
			stepStartEdges = stepEndEdges
		}
		if (open) {
			const endFaceEdges = stepEndEdges.map(edge => edge.flipped()).reverse()
			faces.push(
				new PlaneFace(new PlaneSurface(P3.ZX.flipped()), baseLoop),
				new PlaneFace(new PlaneSurface(P3.ZX.rotateZ(totalRads)), endFaceEdges))
		}
		return new B2(faces)
	}

	/**
	 * loop should be CCW on XZ plane for a bounded B2
	 */
	export function rotateEdgesUnsplit(loop: Edge[], rads: number, name: string): B2 {
		assert(Edge.isLoop(loop))
		const rotationMatrix = M4.rotationZ(rads)
		const open = !NLA.eq(rads, 2 * PI)
		const endEdges = open ? loop.map(edge => edge.transform(rotationMatrix)) : loop
		const edgeCount = loop.length
		const ribs = NLA.arrayFromFunction(edgeCount, i => {
			const a = loop[i].a, radius = a.lengthXY()
			const b = endEdges[i].a
			if (!NLA.eq0(radius)) {
				const curve = new EllipseCurve(V(0, 0, a.z), V(-radius, 0, 0), V(0, -radius, 0))
				const aT = -PI, bT = -PI + rads
				return new PCurveEdge(curve, a, b, aT, bT, null, curve.tangentAt(aT), curve.tangentAt(bT), name + 'rib' + i)
			}
		})
		const faces = loop.map((edge, i) => {
			const ipp = (i + 1) % edgeCount
			console.log("ljl", i, ipp, ribs)
			const faceEdges = [
				edge.flipped(),
				!NLA.eq0(edge.a.x) && ribs[i],
				endEdges[i],
				!NLA.eq0(edge.b.x) && ribs[ipp].flipped()].filter(x => x)
			if (edge instanceof StraightEdge) {
				const line = edge.curve
				let surface
				if (line.dir1.isParallelTo(V3.Z)) {
					if (NLA.eq0(edge.a.x)) {
						return
					}
					let flipped = edge.a.z > edge.b.z
					surface = new CylinderSurface(ribs[i].curve, !flipped ? V3.Z : V3.Z.negated())
				} else if (line.dir1.isPerpendicularTo(V3.Z)) {
					let flipped = edge.a.x > edge.b.x
					let surface = new PlaneSurface(new P3(V3.Z, edge.a.z))
					if (!flipped) surface = surface.flipped()
					if (!open) {
						const hole = flipped
							? !NLA.eq0(edge.b.x) && ribs[ipp].flipped()
							: !NLA.eq0(edge.a.x) && ribs[i]
						return new PlaneFace(surface, [flipped ? ribs[i] : ribs[ipp].flipped()], hole && [[hole]])
					}
					return new PlaneFace(surface, faceEdges)
				} else {
					// apex is intersection of segment with Z-axis
					let a = edge.a, b = edge.b
					let apexZ = a.z - a.x * (b.z - a.z) / (b.x - a.x)
					let apex = new V3(0, 0, apexZ)
					let flipped = edge.a.z > edge.b.z
					surface = ConicSurface.atApexThroughEllipse(apex, ribs[a.x > b.x ? i : ipp].curve as EllipseCurve, !flipped ? 1 : -1)
				}
				return Face.create(surface, faceEdges)
			}
			if (edge.curve instanceof EllipseCurve) {
				let flipped = undefined
				let ell = edge.curve.rightAngled()
				let f1Perp = ell.f1.isPerpendicularTo(V3.Z), f2Perp = ell.f2.isPerpendicularTo(V3.Z)
				if (L3.Z.containsPoint(ell.center) && (f1Perp || f2Perp)) {
					let f3length = f1Perp ? ell.f1.length() : ell.f2.length()
					if (flipped) {
						f3length *= -1
					}
					let surface = new EllipsoidSurface(ell.center, ell.f1, ell.f2, ell.f1.cross(ell.f2).toLength(f3length))
					return new RotationFace(surface, faceEdges)
				}
			} else {
				assert(false, edge)
			}
		}).filter(x => x)
		if (open) {
			const endFaceEdges = endEdges.map(edge => edge.flipped()).reverse()
			faces.push(
				new PlaneFace(new PlaneSurface(P3.ZX.flipped()), loop),
				new PlaneFace(new PlaneSurface(P3.ZX.rotateZ(rads)), endFaceEdges))
		}
		return new B2(faces)
	}

	export function rotStep(edges: Edge[], totalRads: number, count: int) {
		const radStep = totalRads / count
		const open = !NLA.eq(totalRads, 2 * PI)
		const ribCount = !open ? count : count + 1
		const ribs = NLA.arrayFromFunction(ribCount, i => {
			if (i == 0) return edges
			const matrix = M4.rotationZ(radStep * i)
			return edges.map(edge => edge.transform(matrix))
		})
		const hs = NLA.arrayFromFunction(count, i => {
			const ipp = (i + 1) % ribCount
			return NLA.arrayFromFunction(edges.length, j => {
				if (!NLA.eq0(edges[j].a.lengthXY())) {
					return StraightEdge.throughPoints(ribs[i][j].a, ribs[ipp][j].a)
				}
			})
		})
		const faces = []
		let surface, face
		edges.forEach((edge, i) => {
			const ipp = (i + 1) % edges.length
			if (edge instanceof StraightEdge && edge.curve.dir1.isPerpendicularTo(V3.Z)) {
				let flipped = edge.a.x > edge.b.x
				surface = new PlaneSurface(flipped ? new P3(V3.Z, edge.a.z) : new P3(V3.Z.negated(), -edge.a.z))
				if (open) {
					const newEdges = []
					if (!NLA.eq0(edge.a.x)) {
						newEdges.pushAll(NLA.arrayFromFunction(count, j => hs[j][i]))
					}
					newEdges.push(ribs[count][i])
					if (!NLA.eq0(edge.b.x)) {
						newEdges.pushAll(NLA.arrayFromFunction(count, j => hs[count - j - 1][ipp].flipped()))
					}
					newEdges.push(edge.flipped())
					face = new PlaneFace(surface, newEdges)
				} else {
					const contour = flipped
						? NLA.arrayFromFunction(count, j => hs[j][i])
						: NLA.arrayFromFunction(count, j => hs[count - j - 1][ipp].flipped())
					let hole
					if (flipped && !NLA.eq0(edge.b.x)) {
						hole = NLA.arrayFromFunction(count, j => hs[count - j - 1][ipp].flipped())
					} else if (!flipped && !NLA.eq0(edge.a.x)) {
						hole = NLA.arrayFromFunction(count, j => hs[j][i])
					}
					face = new PlaneFace(surface, contour, hole ? [hole] : [])
				}
				faces.push(face)
				return
			} else if (edge instanceof StraightEdge) {
				if (NLA.eq0(edge.a.lengthXY()) && NLA.eq0(edge.b.lengthXY())) {
					return
				}
			}
			for (let r = 0; r < count; r++) {
				const rpp = (r + 1) % ribCount
				const faceEdges = [ribs[r][i].flipped(), hs[r][i], ribs[rpp][i], hs[r][ipp] && hs[r][ipp].flipped()].filter(x => x)
				if (edge instanceof StraightEdge) {
					const surface = new PlaneSurface(P3.throughPoints(faceEdges[0].a, faceEdges[1].a, faceEdges[2].a))
					faces.push(new PlaneFace(surface, faceEdges))
				} else {
					assert(false, edge.toString())
				}
			}
		})
		if (open) {
			const endFaceEdges = ribs[count].map(edge => edge.flipped()).reverse()
			const endFace = new PlaneFace(new PlaneSurface(P3.ZX.rotateZ(totalRads)), endFaceEdges)
			faces.push(new PlaneFace(new PlaneSurface(P3.ZX.flipped()), edges), endFace)
		}
		return new B2(faces)
	}

	export function extrudeVertices(baseVertices, baseFacePlane, offset, name?, source?) {
		assert(baseVertices.every(v => v instanceof V3), "baseVertices.every(v => v instanceof V3)")
		assertInst(P3, baseFacePlane)
		assertVectors(offset)
		if (baseFacePlane.normal.dot(offset) > 0) baseFacePlane = baseFacePlane.flipped()
		if (!isCCW(baseVertices, baseFacePlane.normal)) {
			baseVertices = baseVertices.reverse()
		}
		let topVertices = baseVertices.map((v) => v.plus(offset)).reverse()
		//let topPlane = basePlane.translated(offset)
		let top, bottom
		let faces = [
			bottom = PlaneFace.forVertices(new PlaneSurface(baseFacePlane), baseVertices),
			top = PlaneFace.forVertices(new PlaneSurface(baseFacePlane.flipped().translated(offset)), topVertices)]
		let m = baseVertices.length
		let ribs = NLA.arrayFromFunction(m, i => StraightEdge.throughPoints(baseVertices[i], topVertices[m - 1 - i]))
		for (let i = 0; i < m; i++) {
			let j = (i + 1) % m
			faces.push(
				new PlaneFace(
					PlaneSurface.throughPoints(baseVertices[j], baseVertices[i], topVertices[m - j - 1]),
					[bottom.edges[i].flipped(), ribs[i], top.edges[m - j - 1].flipped(), ribs[j].flipped()], [], name + "wall" + i))
		}
		let edges = StraightEdge.chain(baseVertices, true)
		source = source || `B2T.extrudeVertices(${baseVertices.sce}, ${baseFacePlane.sce}, ${offset.sce}, "${name}")`
		return B2T.extrudeEdges(edges, baseFacePlane, offset, name, source)
	}

	// Returns a tetrahedron (3 sided pyramid).
	// Faces will face outwards.
	// abcd can be in any order. The only constraint is that abcd cannot be on a common plane.
	export function tetrahedron(a: V3, b: V3, c: V3, d: V3, name: string = 'tetra' + globalId++): B2 {
		assertVectors(a, b, c, d)
		const dDistance = P3.throughPoints(a, b, c).distanceToPointSigned(d)
		if (NLA.eq0(dDistance)) {
			throw new Error("four points are coplanar")
		}
		if (dDistance > 0) {
			[c, d] = [d, c]
		}
		const ab = StraightEdge.throughPoints(a, b)
		const ac = StraightEdge.throughPoints(a, c)
		const ad = StraightEdge.throughPoints(a, d)
		const bc = StraightEdge.throughPoints(b, c)
		const bd = StraightEdge.throughPoints(b, d)
		const cd = StraightEdge.throughPoints(c, d)
		const faces = [
			new PlaneFace(PlaneSurface.throughPoints(a, b, c), [ab, bc, ac.flipped()], [], name + 'abc'),
			new PlaneFace(PlaneSurface.throughPoints(a, d, b), [ad, bd.flipped(), ab.flipped()], [], name + 'adb'),
			new PlaneFace(PlaneSurface.throughPoints(b, d, c), [bd, cd.flipped(), bc.flipped()], [], name + 'bdc'),
			new PlaneFace(PlaneSurface.throughPoints(c, d, a), [cd, ad.flipped(), ac], [], name + 'cda')
		]
		const gen = `B2T.tetrahedron(${a.sce}, ${b.sce}, ${c.sce}, ${d.sce})`
		return new B2(faces, false, gen)
	}
    const b = 1 / NLA.GOLDEN_RATIO, c = 2 - NLA.GOLDEN_RATIO
    export const DODECAHEDRON_VERTICES = [
        new V3( c,  0,  1),
        new V3(-c,  0,  1),
        new V3(-b,  b,  b),
        new V3( 0,  1,  c),
        new V3( b,  b,  b),
        new V3( b, -b,  b),
        new V3( 0, -1,  c),
        new V3(-b, -b,  b),
        new V3( c,  0, -1),
        new V3(-c,  0, -1),
        new V3(-b, -b, -b),
        new V3( 0, -1, -c),
        new V3( b, -b, -b),
        new V3( b,  b, -b),
        new V3( 0,  1, -c),
        new V3(-b,  b, -b),
        new V3( 1,  c,  0),
        new V3(-1,  c,  0),
        new V3(-1, -c,  0),
        new V3( 1, -c,  0)
    ].map(v => v.normalized())

    export const DODECAHEDRON_FACE_VERTICES = [
        [  4,  3,  2,  1,  0 ],
        [  7,  6,  5,  0,  1 ],
        [ 12, 11, 10,  9,  8 ],
        [ 15, 14, 13,  8,  9 ],
        [ 14,  3,  4, 16, 13 ],
        [  3, 14, 15, 17,  2 ],
        [ 11,  6,  7, 18, 10 ],
        [  6, 11, 12, 19,  5 ],
        [  4,  0,  5, 19, 16 ],
        [ 12,  8, 13, 16, 19 ],
        [ 15,  9, 10, 18, 17 ],
        [  7,  1,  2, 17, 18 ]].concatenated()

    export const OCTAHEDRON_VERTICES = [
        new V3( 1,  0,  0),
        new V3(-1,  0,  0),
        new V3( 0,  1,  0),
        new V3( 0,  -1,  0),
        new V3( 0,  0,  1),
        new V3( 0,  0,  -1)]

    export const OCTAHEDRON_FACE_VERTICES = [
        0, 2, 4,
        2, 1, 4,
        1, 3, 4,
        3, 0, 4,

        2, 0, 5,
        1, 2, 5,
        3, 1, 5,
        0, 3, 5]

    export function dodecahedron() {
        return makePlatonic(DODECAHEDRON_VERTICES, DODECAHEDRON_FACE_VERTICES, 5, 'B2T.dodecahedron()')
    }
    export function octahedron() {
        return makePlatonic(OCTAHEDRON_VERTICES, OCTAHEDRON_FACE_VERTICES, 3, 'B2T.octahedron()')
    }
	function makePlatonic(VS, FVIS, FACE_EDGE_COUNT, generator) {
	    const edgeMap = new Map(), faces = []
        for (let start = 0; start < FVIS.length; start += FACE_EDGE_COUNT) {
            console.log(start + 1, FVIS[start + 1], VS[FVIS[start + 1]])
            const surface = PlaneSurface.throughPoints(VS[FVIS[start]], VS[FVIS[start + 1]], VS[FVIS[start + 2]])
            const contour = []
	        for (let i = 0; i < FACE_EDGE_COUNT; i++) {
                const ipp = (i + 1) % FACE_EDGE_COUNT
                const iA = FVIS[start + i], iB = FVIS[start + ipp]
                const iMin = min(iA, iB), iMax = max(iA, iB), edgeID = iMin * VS.length + iMax
                let edge = edgeMap.get(edgeID)
                !edge && edgeMap.set(edgeID, edge = StraightEdge.throughPoints(VS[iMin], VS[iMax]))
                contour.push(iA < iB ? edge : edge.flipped())
            }
            faces.push(new PlaneFace(surface, contour))
        }
        return new B2(faces, false, 'B2T.dodecahedron()')
    }

	export function pyramidEdges(baseEdges: Edge[], apex: V3, name: string = 'pyramid' + globalId++): B2 {
		assertInst.apply(undefined, [Edge].concat(baseEdges))
		assertVectors(apex)

		const ribs = baseEdges.map(baseEdge => StraightEdge.throughPoints(apex, baseEdge.a))
		const faces = baseEdges.map((baseEdge, index) => {
			const faceName = name + 'Wall' + index
			const nextIndex = (index + 1) % baseEdges.length
			const faceEdges = [ribs[index], baseEdge, ribs[nextIndex].flipped()]
			const surface = undefined // TODO
			return Face.create(surface, faceEdges, undefined, faceName)
		})
		faces.push(bottomFace)
		const gen = makeGen('B2T.pyramidEdges', baseEdges, apex, name)
		return new B2(faces, false, gen, name)
	}
}

function makeGen(name: string, ...params: { toSource(): string }[]) {
	return name + '(' + params.map(p => p.toSource()).join(', ') + ')'
}