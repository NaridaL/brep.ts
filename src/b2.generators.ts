function projectCurve(curve: Curve, offset: V3, flipped: boolean): Surface {
	if (curve instanceof L3) {
		const surfaceNormal = offset.cross(curve.dir1).toLength(flipped ? -1 : 1)
		return new PlaneSurface(P3.normalOnAnchor(surfaceNormal, curve.anchor))
	}
	if (curve instanceof SemiEllipseCurve) {
		const curveDir = flipped ? offset : offset.negated()
		return new SemiCylinderSurface(curve, curveDir.unit())
	}
	if (curve instanceof BezierCurve || curve instanceof ParabolaCurve) {
		const curveDir = flipped ? offset : offset.negated()
		return new ProjectedCurveSurface(curve, curveDir.unit(), 0, 1)
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
			surface = new SemiCylinderSurface(ribs[i].curve, !flipped ? V3.Z : V3.Z.negated())
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
			surface = ConicSurface.atApexThroughEllipse(apex, ribs[a.x > b.x ? i : ipp].curve as SemiEllipseCurve, !flipped ? 1 : -1)
		}
		return Face.create(surface, faceEdges)
	}
	if (edge.curve instanceof SemiEllipseCurve) {
		let flipped = undefined
		let ell = edge.curve.rightAngled()
		let f1Perp = ell.f1.isPerpendicularTo(V3.Z), f2Perp = ell.f2.isPerpendicularTo(V3.Z)
		if (L3.Z.containsPoint(ell.center) && (f1Perp || f2Perp)) {
			let f3length = f1Perp ? ell.f1.length() : ell.f2.length()
			if (flipped) {
				f3length *= -1
			}
			let surface = new SemiEllipsoidSurface(ell.center, ell.f1, ell.f2, ell.f1.cross(ell.f2).toLength(f3length))
			return new RotationFace(surface, faceEdges)
		}
	} else {
		assert(false, edge)
	}
	if (curve instanceof L3) {
		let surfaceNormal = offset.cross(curve.dir1).toLength(flipped ? -1 : 1)
		return new PlaneSurface(P3.normalOnAnchor(surfaceNormal, curve.anchor))
	}
	if (curve instanceof SemiEllipseCurve) {
		let curveDir = flipped ? offset : offset.negated()
		return new SemiCylinderSurface(curve, curveDir.unit())
	}
	if (curve instanceof BezierCurve) {
		let curveDir = flipped ? offset : offset.negated()
		return new ProjectedCurveSurface(curve, curveDir.unit(), 0, 1)
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
	    const generator = makeGen('B2T.box', w, h, d, name)
		return B2T.extrudeVertices(baseVertices, P3.XY.flipped(), new V3(0, 0, d), name, generator)
	}

	export function puckman(radius: number, rads: raddd, height: number, name: string): B2 {
		assertf(() => NLA.lt(0, radius))
		assertf(() => NLA.lt(0, rads) && NLA.le(rads, TAU))
		assertf(() => NLA.lt(0, height))
		const edges = StraightEdge.chain([V3.O, new V3(radius, 0, 0), new V3(radius, 0, height), new V3(0, 0, height)], true)
		return B2T.rotateEdges(edges, rads, name || 'puckman' + globalId++)
	}

	export function registerVertexName(map, name, p) {
		// TODO
		if (!Array.from(map.keys()).some(p2 => p2.like(p))) {
			map.set(p, name)
		}
	}

	export function extrudeEdges(baseFaceEdges: Edge[],
	                             baseFacePlane: P3 = P3.XY,
	                             offset: V3 = V3.Z,
	                             name: string = 'extrude' + globalId++,
	                             gen?: string,
	                             infoFactory?: FaceInfoFactory<any>): B2 {
    	baseFaceEdges = fixEdges(baseFaceEdges)
		//Array.from(NLA.combinations(baseFaceEdges.length)).forEach(({i, j}) => {
		//	assertf(() => !Edge.edgesIntersect(baseFaceEdges[i], baseFaceEdges[j]), baseFaceEdges[i].sce + baseFaceEdges[j].sce)
		//})
		assertf(() => Edge.isLoop(baseFaceEdges))
		// TODO checks..
		//if (offset.dot(baseFacePlane.normal1) > 0) {
		//	baseFacePlane = baseFacePlane.flipped()
		//}
		const vertexNames = new Map()
		const basePlaneSurface = new PlaneSurface(baseFacePlane)
		//assert(basePlaneSurface.edgeLoopCCW(baseFaceEdges), 'edges not CCW on baseFacePlane')
		const translationMatrix = M4.translation(offset)
		const topEdges = baseFaceEdges.map(edge => edge.transform(translationMatrix, 'top')) as Edge[]
		const edgeCount = baseFaceEdges.length

		const bottomInfo = infoFactory && infoFactory.extrudeBottom(basePlaneSurface, baseFaceEdges)
		const bottomFace = new PlaneFace(basePlaneSurface, baseFaceEdges, [], name + 'Bottom', bottomInfo)

		const topFaceEdges = topEdges.map(edge => edge.flipped()).reverse()
		const topSurface = new PlaneSurface(baseFacePlane.flipped().translated(offset))
		const topInfo = infoFactory && infoFactory.extrudeBottom(topSurface, topFaceEdges)
		const topFace = new PlaneFace(topSurface, topFaceEdges, [], name + 'Top', topInfo)


		baseFaceEdges.forEach(edge => B2T.registerVertexName(vertexNames, edge.name + 'A', edge.a))
		topFaceEdges.forEach(edge => B2T.registerVertexName(vertexNames, edge.name + 'A', edge.a))

		const ribs = NLA.arrayFromFunction(edgeCount,
			i => StraightEdge.throughPoints(baseFaceEdges[i].a, topEdges[i].a, name + 'Rib' + i))

		const faces = baseFaceEdges.map((edge, i) => {
			const faceName = name + 'Wall' + i
			const j = (i + 1) % edgeCount
			const faceEdges = [baseFaceEdges[i].flipped(), ribs[i], topEdges[i], ribs[j].flipped()]
			const surface = projectCurve(edge.curve, offset, edge.reversed)
			const info = infoFactory && infoFactory.extrudeWall(i, surface, faceEdges)
			return Face.create(surface, faceEdges, undefined, faceName,	info)
		}) as Face[]
		faces.push(bottomFace, topFace)
		gen = gen || makeGen('B2T.extrudeEdges', baseFaceEdges, baseFacePlane, offset, name)
		return new B2(faces, false, gen, vertexNames)
	}


	export function cylinder(radius: number = 1, height: number = 1, rads: raddd = TAU, name: string = 'cylinder' + globalId++): B2 {
		const vertices = [new V3(0, 0, 0), new V3(radius, 0, 0), new V3(radius, 0, height), new V3(0, 0, height)]
		return rotateEdges(StraightEdge.chain(vertices, true), rads, name)
	}

	export function cone(radius: number = 1, height: number = 1, rads: raddd = TAU, name: string = 'cone' + globalId++): B2 {
		const vertices = [new V3(0, 0, 0), new V3(radius, 0, height), new V3(0, 0, height)]
		return rotateEdges(StraightEdge.chain(vertices, true), rads, name)
	}

    export function sphere(radius: number = 1, name: string = 'sphere' + globalId++, rot: raddd = TAU): B2 {
        const ee = PCurveEdge.create(
            new SemiEllipseCurve(V3.O, new V3(0, 0, -radius), new V3(radius, 0, 0)),
            new V3(0, 0, -radius), new V3(0, 0, radius),
            0, PI,
            undefined,
            new V3(radius, 0, 0), new V3(-radius, 0, 0))
	    const generator = makeGen('B2T.sphere', radius, name, rot)
	    return rotateEdges([StraightEdge.throughPoints(ee.b, ee.a), ee], rot, name, generator)
    }

    export function menger(res: int = 2, name: string = 'menger' + globalId++): B2 {
        let result = B2T.box(1,1,1)
        if (0 == res) return result
        const punch = B2T.box(1/3, 1/3, 2).translate(1/3, 1/3, -1/2).flipped()
        function recurse(steps: int, m4: M4) {
            result = result.and(punch.transform(m4))
            a = result
            if (steps > 1) {
                const scaled = m4.times(M4.scaling(1/3, 1/3, 1))
                for (let i = 0; i < 9; i++) {
                    if (4 == i) continue
                    recurse(steps - 1, scaled.times(M4.translation(i % 3, i / 3 | 0, 0)))
                }
            }
        }
        recurse(res, M4.IDENTITY)
        recurse(res, M4.YZX)
        recurse(res, M4.ZXY)
        return result
    }
    export function menger2(res: int = 2, name: string = 'menger' + globalId++): B2 {
	    if (0 == res) return B2T.box(1,1,1)

        const punch = B2T.box(1/3, 1/3, 2).translate(1/3, 1/3, -1/2).flipped()
        const stencilFaces = []
        function recurse(steps: int, m4: M4) {
            stencilFaces.pushAll(punch.transform(m4).faces)
            if (steps > 1) {
                const scaled = m4.times(M4.scaling(1/3, 1/3, 1))
                for (let i = 0; i < 9; i++) {
                    if (4 == i) continue
                    recurse(steps - 1, scaled.times(M4.translation(i % 3, i / 3 | 0, 0)))
                }
            }
        }
        recurse(res, M4.IDENTITY)
        const stencil = new B2(stencilFaces, true)

        return B2T.box()
            .and(stencil)
            .and(stencil.transform(M4.YZX))
            .and(stencil.transform(M4.ZXY))
    }

	export function torus(rSmall: number, rLarge: number, rads: raddd, name: string): B2 {
		assertNumbers(rSmall, rLarge, rads)
		assertf(() => rLarge > rSmall)
		const curve = SemiEllipseCurve.semicircle(rSmall, new V3(rLarge, 0, 0))
		const baseEdges = [PCurveEdge.forCurveAndTs(curve, -Math.PI, 0), PCurveEdge.forCurveAndTs(curve, 0, Math.PI)]
		return B2T.rotateEdges(baseEdges, rads, name || 'torus' + globalId++)
	}
	export function torusUnsplit(rSmall: number, rLarge: number, rads: raddd, name: string): B2 {
		assertNumbers(rSmall, rLarge, rads)
		assertf(() => rLarge > rSmall)
		const baseEdge = PCurveEdge.forCurveAndTs(SemiEllipseCurve.semicircle(rSmall, new V3(rLarge, 0, 0)), -Math.PI, Math.PI)
		return B2T.rotateEdges([baseEdge], rads, name || 'torus' + globalId++)
	}

	/**
	 * baseLoop should be CCW on XZ plane for a bounded B2
	 */
	export function rotateEdges(baseLoop: Edge[], totalRads: raddd, name: string, generator?: string): B2 {
		assert(!NLA.eq(PI, totalRads) || PI == totalRads) // URHGJ
		assertf(() => NLA.lt(0, totalRads) && NLA.le(totalRads, TAU))
        totalRads = NLA.snap(totalRads, TAU)
		assertf(() => Edge.isLoop(baseLoop))
		const rotationSteps = ceil((totalRads - NLA_PRECISION) / PI)
		const angles = rotationSteps == 1 ? [-PI, -PI + totalRads] : [-PI, 0, totalRads - PI]
		const open = !NLA.eq(totalRads, 2 * PI)
		const baseRibCurves = baseLoop.map(edge =>  {
			const a = edge.a, radius = a.lengthXY()
			if (!NLA.eq0(radius)) {
				return new SemiEllipseCurve(V(0, 0, a.z), V(radius, 0, 0), V(0, radius, 0))
			}
		})
		const baseSurfaces = baseLoop.map((edge, i) => {
			const ipp = (i + 1) % baseLoop.length
			if (edge instanceof StraightEdge) {
				const line = edge.curve
				if (line.dir1.isParallelTo(V3.Z)) {
					if (NLA.eq0(edge.a.x)) {
						return
					}
					const flipped = edge.a.z > edge.b.z
					return new SemiCylinderSurface(baseRibCurves[i], !flipped ? V3.Z : V3.Z.negated())
				} else if (line.dir1.isPerpendicularTo(V3.Z)) {
					const flipped = edge.a.x > edge.b.x
					let surface = new PlaneSurface(new P3(V3.Z, edge.a.z))
					if (!flipped) surface = surface.flipped()
					return surface
				} else {
					// apex is intersection of segment with Z-axis
					const a = edge.a, b = edge.b
					const apexZ = a.z - a.x * (b.z - a.z) / (b.x - a.x)
					const apex = new V3(0, 0, apexZ)
					const flipped = edge.a.z > edge.b.z
					return ConicSurface.atApexThroughEllipse(apex, baseRibCurves[a.x > b.x ? i : ipp] as SemiEllipseCurve, !flipped ? 1 : -1)
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
			if (edge.curve instanceof SemiEllipseCurve) {
				const flipped = edge.a.z > edge.b.z
				const ell = edge.curve.rightAngled()
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
                return SemiEllipsoidSurface.forABC(width, (!flipped ? 1 : -1) * width, height, ell.center)
			} else {
				assert(false, edge)
			}
		})
		let stepStartEdges = baseLoop, stepEndEdges
		const faces = []
		for (let rot = 0; rot < totalRads; rot += PI) {
			const aT = 0, bT = min(totalRads - rot, PI)
			const rotation = M4.rotationZ(rot + bT), rotrot = M4.rotationZ(rot)
			stepEndEdges = rot + bT == TAU ? baseLoop : baseLoop.map(edge => edge.transform(rotation))
			const ribs = NLA.arrayFromFunction(baseLoop.length, i => {
				const a = stepStartEdges[i].a, radius = a.lengthXY()
				const b = stepEndEdges[i].a
				if (!NLA.eq0(radius)) {
					const curve = 0 == rot ? baseRibCurves[i] : baseRibCurves[i].rotateZ(rot)
					return new PCurveEdge(curve, a, b, aT, bT, null, curve.tangentAt(aT), curve.tangentAt(bT), name + 'rib' + i)
				}
			})
			for (let edgeIndex = 0; edgeIndex < baseLoop.length; edgeIndex++) {
				if (baseSurfaces[edgeIndex]) {
					const edge = stepStartEdges[edgeIndex]
					const ipp = (edgeIndex + 1) % baseLoop.length
					const faceEdges = [
						stepStartEdges[edgeIndex].flipped(),
						!NLA.eq0(edge.a.x) && ribs[edgeIndex],
						stepEndEdges[edgeIndex],
						!NLA.eq0(edge.b.x) && ribs[ipp].flipped()].filter(x => x)
                    const surface = 0 == rot ? baseSurfaces[edgeIndex] : baseSurfaces[edgeIndex].rotateZ(rot)
					faces.push(Face.create(surface, faceEdges))
				}
			}
			stepStartEdges = stepEndEdges
		}
		if (open) {
			const endFaceEdges = Edge.reverseLoop(stepEndEdges)
			faces.push(
				new PlaneFace(new PlaneSurface(P3.ZX.flipped()), baseLoop),
				new PlaneFace(new PlaneSurface(P3.ZX.rotateZ(totalRads)), endFaceEdges))
		}
		const infiniteVolume = new PlaneSurface(P3.ZX).edgeLoopCCW(baseLoop)
		return new B2(faces, infiniteVolume, generator)
	}

	/**
	 * loop should be CCW on XZ plane for a bounded B2
	 */
	//export function rotateEdgesUnsplit(loop: Edge[], rads: raddd, name: string): B2 {
	//	assert(Edge.isLoop(loop))
	//	const rotationMatrix = M4.rotationZ(rads)
	//	const open = !NLA.eq(rads, 2 * PI)
	//	const endEdges = open ? loop.map(edge => edge.transform(rotationMatrix)) : loop
	//	const edgeCount = loop.length
	//	const ribs = NLA.arrayFromFunction(edgeCount, i => {
	//		const a = loop[i].a, radius = a.lengthXY()
	//		const b = endEdges[i].a
	//		if (!NLA.eq0(radius)) {
	//			const curve = new SemiEllipseCurve(V(0, 0, a.z), V(-radius, 0, 0), V(0, -radius, 0))
	//			const aT = -PI, bT = -PI + rads
	//			return new PCurveEdge(curve, a, b, aT, bT, null, curve.tangentAt(aT), curve.tangentAt(bT), name + 'rib' + i)
	//		}
	//	})
	//	const faces = loop.map((edge, i) => {
	//		const ipp = (i + 1) % edgeCount
	//		console.log('ljl', i, ipp, ribs)
	//		const faceEdges = [
	//			edge.flipped(),
	//			!NLA.eq0(edge.a.x) && ribs[i],
	//			endEdges[i],
	//			!NLA.eq0(edge.b.x) && ribs[ipp].flipped()].filter(x => x)
	//		if (edge instanceof StraightEdge) {
	//			const line = edge.curve
	//			let surface
	//			if (line.dir1.isParallelTo(V3.Z)) {
	//				if (NLA.eq0(edge.a.x)) {
	//					return
	//				}
	//				let flipped = edge.a.z > edge.b.z
	//				surface = new SemiCylinderSurface(ribs[i].curve, !flipped ? V3.Z : V3.Z.negated())
	//			} else if (line.dir1.isPerpendicularTo(V3.Z)) {
	//				let flipped = edge.a.x > edge.b.x
	//				let surface = new PlaneSurface(new P3(V3.Z, edge.a.z))
	//				if (!flipped) surface = surface.flipped()
	//				if (!open) {
	//					const hole = flipped
	//						? !NLA.eq0(edge.b.x) && ribs[ipp].flipped()
	//						: !NLA.eq0(edge.a.x) && ribs[i]
	//					return new PlaneFace(surface, [flipped ? ribs[i] : ribs[ipp].flipped()], hole && [[hole]])
	//				}
	//				return new PlaneFace(surface, faceEdges)
	//			} else {
	//				// apex is intersection of segment with Z-axis
	//				let a = edge.a, b = edge.b
	//				let apexZ = a.z - a.x * (b.z - a.z) / (b.x - a.x)
	//				let apex = new V3(0, 0, apexZ)
	//				let flipped = edge.a.z > edge.b.z
	//				surface = ConicSurface.atApexThroughEllipse(apex, ribs[a.x > b.x ? i : ipp].curve as SemiEllipseCurve, !flipped ? 1 : -1)
	//			}
	//			return Face.create(surface, faceEdges)
	//		}
	//		if (edge.curve instanceof SemiEllipseCurve) {
	//			let flipped = undefined
	//			let ell = edge.curve.rightAngled()
	//			let f1Perp = ell.f1.isPerpendicularTo(V3.Z), f2Perp = ell.f2.isPerpendicularTo(V3.Z)
	//			if (L3.Z.containsPoint(ell.center) && (f1Perp || f2Perp)) {
	//				let f3length = f1Perp ? ell.f1.length() : ell.f2.length()
	//				if (flipped) {
	//					f3length *= -1
	//				}
	//				let surface = new SemiEllipsoidSurface(ell.center, ell.f1, ell.f2, ell.f1.cross(ell.f2).toLength(f3length))
	//				return new RotationFace(surface, faceEdges)
	//			}
	//		} else {
	//			assert(false, edge)
	//		}
	//	}).filter(x => x)
	//	if (open) {
	//		const endFaceEdges = endEdges.map(edge => edge.flipped()).reverse()
	//		faces.push(
	//			new PlaneFace(new PlaneSurface(P3.ZX.flipped()), loop),
	//			new PlaneFace(new PlaneSurface(P3.ZX.rotateZ(rads)), endFaceEdges))
	//	}
	//	return new B2(faces, undefined)
	//}

	export function quaffle() {
	    //const baseK = B2T.sphere(0.2).translate(0, 0.95).flipped()
	    const baseK = B2T.box().scale(0.2).translate(0, 0.95).flipped()
        const vs = B2T.DODECAHEDRON_VERTICES.concat(
            B2T.DODECAHEDRON_FACE_VERTICES.map(fis => fis
                    .map(vi => B2T.DODECAHEDRON_VERTICES[vi])
                    .reduce((a,b) => a.plus(b), V3.O)
                    .unit()))
        vs.forEach(v => B2T.sphere().and(baseK.rotateAB(V3.Y, v)))
        const ss = new B2(vs.map(v => baseK.rotateAB(V3.Y, v).faces).concatenated(), false)
        return ss
    }

    export function extrudeFace(face, dir) {
	    return new B2(
            extrudeEdges(face.contour, face.surface.plane, dir).faces.slice(0, -2).concat(
	        face, face.translate(dir.x,dir.y,dir.z).flipped(),
	        face.holes.map(
	            hole =>
                    extrudeEdges(hole, face.surface.plane.flipped(), dir).faces.slice(0, -2)).concatenated()), false)
    }

    let defaultFont: opentypejs.Font
    export function loadFonts(): Promise<opentypejs.Font> {
		return loadFont('fonts/FiraSansMedium.woff').then(font => defaultFont = font)
    }
    const loadedFonts = new Map<string, opentypejs.Font>()
    export function loadFont(fontPath): Promise<opentypejs.Font> {
        return new Promise(function (executor, reject) {
        	const font = loadedFonts.get(fontPath)
            if (font) {
                executor(font)
            } else {
                opentype.load(fontPath, function (err, f) {
                    if (err) {
                        reject(err)
                    } else {
                        loadedFonts.set(fontPath, f)
                        executor(f)
                    }
                })
            }
        })
    }
    export function loadFontsAsync(callback) {
        if (defaultFont) {
            callback()
        } else {
            opentype.load('fonts/FiraSansMedium.woff', function (err, font) {
                if (err) {
                    throw new Error('Could not load font: ' + err)
                } else {
                    defaultFont = font
                    callback()
                }
            })
        }
    }
    export function text(text: string, size: number, depth: number = 1, font: opentypejs.Font = defaultFont) {
        const path = font.getPath(text, 0, 0, size)
        const subpaths = []
        path.commands.forEach(c => {
            if (c.type == 'M') {
                subpaths.push([])
            }
            subpaths.last().push(c)
        })
        const loops = subpaths.map(sp => {
            const path = new opentype.Path()
            path.commands = sp
            const loop = Edge.reverseLoop(Edge.pathFromSVG(path.toPathData(13))).map(e => e.mirroredY())
            assert(Edge.isLoop(loop))
            return loop
        })
        const faces = Face.assembleFacesFromLoops(loops, new PlaneSurface(P3.XY), PlaneFace)
	    const generator = `B2T.text(${text.sce}, ${size}, ${depth})`
	    const hello = B2.join(faces.map(face => B2T.extrudeFace(face, V(0,0,-depth))), generator)
        return hello

    }

    export function minorityReport() {
	    const a = B2T.sphere()
	    const b = B2T.text('LEO CROW', 64, 128).scale(0.1 / 32).translate(-0.5, -0.05, 1.2).flipped()
	    const c = B2T.sphere(0.98)
	    return a.and(b).plus(c)
    }

    export function whatever() {
        const iso = isocahedron()
        const numbersB2 = B2.join(iso.faces.map((face, i) => {
            const numberB2 = text('' + (i + 1), 0.4, -2)
            const centroid = face.contour.map(edge => edge.a).reduce((a, b) => a.plus(b), V3.O).div(3)

            const sys = M4.forSys(
                face.contour[0].aDir,
                centroid.cross(face.contour[0].aDir),
                centroid.unit(),
                centroid)
            return numberB2.transform(sys.times(M4.translation(-numberB2.getAABB().size().x / 2, -0.1, -0.04)))
        }))
        const s = sphere(0.9)
        //return iso.and(numbersB2)
        return iso.and(s).and(numbersB2)
        //return numbersB2
    }

    export function d20() {
        const iso = isocahedron()
        const numbersB2 = B2.join(iso.faces.map((face, i) => {
            const numberB2 = text('' + (i + 1), 0.4, -2)
            const centroid = face.contour.map(edge => edge.a).reduce((a, b) => a.plus(b), V3.O).div(3)

            const sys = M4.forSys(
                face.contour[0].aDir,
                centroid.cross(face.contour[0].aDir),
                centroid.unit(),
                centroid)
            return numberB2.transform(sys.times(M4.translation(-numberB2.getAABB().size().x / 2, -0.1, -0.04)))
        }))
        const s = sphere(0.9)
        //return iso.and(numbersB2)
        return iso.and(s).and(numbersB2)
        //return numbersB2
    }

	export function rotStep(edges: Edge[], totalRads: raddd, count: int) {
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

	export function fixEdges(edges: Edge[]): Edge[] {
		return edges.map(edge => {
			const c = edge.curve
			if (c instanceof EllipseCurve) {
				const splitEdges = (edge.minT < 0 && edge.maxT > 0)
					? edge.split(0)
					: [edge]
				return splitEdges.map(edge => {
					if (edge.minT >= 0) {
						return Edge.create(new SemiEllipseCurve(c.center, c.f1, c.f2, max(0, c.tMin), c.tMax),
							edge.a, edge.b,
							edge.aT, edge.bT,
							undefined,
							edge.aDir, edge.bDir,
							edge.name)
					} else {
						// "rotate" the curve
						return Edge.create(new SemiEllipseCurve(c.center, c.f1.negated(), c.f2.negated(), c.tMin + PI, min(PI, c.tMax + PI)),
							edge.a, edge.b,
							edge.aT + PI, edge.bT + PI,
							undefined,
							edge.aDir, edge.bDir,
							edge.name)
					}
				})
			}
			if (c instanceof BezierCurve) {
				if (edge.a.like(edge.b)) {
					return edge.split(lerp(edge.aT, edge.bT, 0.5))
				}
			}
			return edge
		}).concatenated() as Edge[]
	}

	export function extrudeVertices(baseVertices: V3[], baseFacePlane: P3, offset: V3, name?, generator?) {
		assert(baseVertices.every(v => v instanceof V3), 'baseVertices.every(v => v instanceof V3)')
		assertInst(P3, baseFacePlane)
		assertVectors(offset)
		if (baseFacePlane.normal1.dot(offset) > 0) baseFacePlane = baseFacePlane.flipped()
		//if (!isCCW(baseVertices, baseFacePlane.normal1)) {
		//	baseVertices = baseVertices.reverse()
		//}
		//let topVertices = baseVertices.map((v) => v.plus(offset)).reverse()
		//let topPlane = basePlane.translated(offset)
		//let top, bottom
		//let faces = [
		//	bottom = PlaneFace.forVertices(new PlaneSurface(baseFacePlane), baseVertices),
		//	top = PlaneFace.forVertices(new PlaneSurface(baseFacePlane.flipped().translated(offset)), topVertices)]
		//let m = baseVertices.length
		//let ribs = NLA.arrayFromFunction(m, i => StraightEdge.throughPoints(baseVertices[i], topVertices[m - 1 - i]))
		//for (let i = 0; i < m; i++) {
		//	let j = (i + 1) % m
		//	faces.push(
		//		new PlaneFace(
		//			PlaneSurface.throughPoints(baseVertices[j], baseVertices[i], topVertices[m - j - 1]),
		//			[bottom.contour[i].flipped(), ribs[i], top.contour[m - j - 1].flipped(), ribs[j].flipped()], [], name + 'wall' + i))
		//}
		let edges = StraightEdge.chain(baseVertices, true)
		generator = generator || makeGen('B2T.extrudeVertices', baseVertices, baseFacePlane, offset, name)
		return B2T.extrudeEdges(edges, baseFacePlane, offset, name, generator)
	}

	export function bug() {
		const good = `M22.2656,33.2348
		L25.7646,29.7638
		C26.7356,31.2938 27.2636,33.0868 27.2636,34.9698
		C27.2636,37.5888 26.2356,40.0578 24.4126,41.9118
		L21.1486,45.1468
		C19.2956,46.9698 16.8556,47.9998 14.2086,47.9998
		C11.5926,47.9998 9.1236,46.9698 7.2416,45.1178
		L2.8896,40.7638
		C-0.9634,36.9408 -0.9634,30.7048 2.8896,26.8818
		L6.1226,23.6468
		C7.9766,21.8228 10.4466,20.7928 13.0626,20.7928
		C14.9736,20.7928 16.7676,21.3228 18.3256,22.3228
		L14.8556,25.7938
		C14.2686,25.5878 13.6806,25.4708 13.0626,25.4708
		C11.6806,25.4708 10.3866,25.9988 9.4166,26.9708
		L6.1816,30.2048
		C4.1826,32.2048 4.1826,35.4698 6.1816,37.4698
		L10.5636,41.8228
		C11.5336,42.7938 12.8286,43.3228 14.2086,43.3228
		C15.5916,43.3228 16.8856,42.7938 17.8546,41.8228
		L21.0896,38.5878
		C22.5306,37.1468 22.9126,35.0588 22.2656,33.2348
		L22.2656,33.2348
		L22.2656,33.2348
		Z
		M33.6736,17.6458
		L30.3516,14.3228
		L14.0626,30.6178
		L17.3846,33.9408
		L33.6736,17.6458
		L33.6736,17.6458
		Z
		M45.1116,21.1168
		L41.9066,24.3228
		C40.0536,26.1748 37.5846,27.2048 34.9686,27.2048
		C33.0866,27.2048 31.2926,26.6768 29.7646,25.7058
		L33.2326,22.2348
		C33.7916,22.4108 34.3496,22.5298 34.9686,22.5298
		C36.3506,22.5298 37.6136,21.9998 38.5846,21.0298
		L41.8176,17.7938
		C43.8186,15.7938 43.8186,12.5288 41.8176,10.5288
		L37.4676,6.1758
		C36.4966,5.2058 35.2036,4.6758 33.8216,4.6758
		C32.4396,4.6758 31.1446,5.2058 30.2056,6.1758
		L26.9696,9.4118
		C26.0006,10.3528 25.4706,11.6468 25.4706,13.0288
		C25.4706,13.6458 25.5876,14.2348 25.7946,14.7938
		L22.3246,18.2648
		C21.3246,16.7348 20.7946,14.9108 20.7946,13.0288
		C20.7946,10.4118 21.8246,7.9408 23.6476,6.0878
		L26.8826,2.8518
		C28.7346,1.0288 31.2046,-0.0002 33.8216,-0.0002
		C36.4376,-0.0002 38.9076,1.0288 40.7596,2.8518
		L45.1116,7.2338
		C48.9626,11.0578 48.9626,17.2938 45.1116,21.1168
		L45.1116,21.1168
		L45.1116,21.1168
		Z`
		const bad = `M22.266 33.235
        l3.499-3.471
        a9.684 9.684 0 0 1 1.499 5.206 9.867 9.867 0 0 1-2.851 6.942
        l-3.264 3.235
        C19.296 46.97 16.856 48 14.209 48
        c-2.616 0-5.085-1.03-6.967-2.882
        L2.89 40.764
        c-3.853-3.823-3.853-10.06 0-13.882
        l3.233-3.235
        a9.866 9.866 0 0 1 6.94-2.854
        c1.91 0 3.705.53 5.263 1.53
        l-3.47 3.47
        a5.352 5.352 0 0 0-1.793-.322 5.106 5.106 0 0 0-3.646 1.5
        l-3.235 3.234
        a5.15 5.15 0 0 0 0 7.265
        l4.382 4.353
        c.97.97 2.265 1.5 3.645 1.5
        a5.107 5.107 0 0 0 3.646-1.5
        l3.235-3.235
        c1.44-1.441 1.823-3.53 1.176-5.353
        z
        m11.408-15.59
        l-3.322-3.322-16.29 16.295 3.323 3.323 16.289-16.295
        z
        m11.438 3.472
        l-3.205 3.206
        a9.768 9.768 0 0 1-6.938 2.882 9.674 9.674 0 0 1-5.204-1.5
        l3.468-3.47
        c.559.176 1.117.295 1.736.295 1.382 0 2.645-.53 3.616-1.5
        l3.233-3.236
        a5.147 5.147 0 0 0 0-7.265
        l-4.35-4.353
        a5.116 5.116 0 0 0-3.646-1.5
        c-1.382 0-2.677.53-3.616 1.5
        L26.97 9.412
        c-.97.94-1.5 2.235-1.5 3.617 0 .617.118 1.206.325 1.765
        l-3.47 3.47
        c-1-1.53-1.53-3.353-1.53-5.235
        a9.868 9.868 0 0 1 2.853-6.941
        l3.235-3.236
        A9.86 9.86 0 0 1 33.822 0
        a9.86 9.86 0 0 1 6.938 2.852
        l4.352 4.382
        c3.85 3.824 3.85 10.06 0 13.883
        z`
		const c = `M10 315
           L 110 215
           A 30 50 0 0 1 162.55 162.45
           L 172.55 152.45
           A 30 50 -45 0 1 215.1 109.9
           L 315 10`
		const g = new SVGPathData(bad).encode()
		const b = new SVGPathData(bad).ySymmetry(30).encode()
		return Edge.pathFromSVG(bad).concat(Edge.pathFromSVG(b).map(e=>e.translate(300)))
		//return Edge.pathFromSVG(c)
	}


	// Returns a tetrahedron (3 sided pyramid).
	// Faces will face outwards.
	// abcd can be in any order. The only constraint is that abcd cannot be on a common plane.
	export function tetrahedron(a: V3, b: V3, c: V3, d: V3, name: string = 'tetra' + globalId++): B2 {
		assertVectors(a, b, c, d)
		const dDistance = P3.throughPoints(a, b, c).distanceToPointSigned(d)
		if (NLA.eq0(dDistance)) {
			throw new Error('four points are coplanar')
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
    ].map(v => v.unit())
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
        [  7,  1,  2, 17, 18 ]]

    export const OCTAHEDRON_VERTICES = [
        new V3( 1,  0,  0),
        new V3(-1,  0,  0),
        new V3( 0,  1,  0),
        new V3( 0,  -1,  0),
        new V3( 0,  0,  1),
        new V3( 0,  0,  -1)]
    export const OCTAHEDRON_FACE_VERTICES = [
        [0, 2, 4],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],

        [2, 0, 5],
        [1, 2, 5],
        [3, 1, 5],
        [0, 3, 5]]

    const {x: s, y: t} = new V3(1, NLA.GOLDEN_RATIO, 0).unit()
    export const ISOCAHEDRON_VERTICES = [
        new V3(-s, t, 0),
        new V3(s, t, 0),
        new V3(-s, -t, 0),
        new V3(s, -t, 0),

        new V3(0, -s, t),
        new V3(0, s, t),
        new V3(0, -s, -t),
        new V3(0, s, -t),

        new V3(t, 0, -s),
        new V3(t, 0, s),
        new V3(-t, 0, -s),
        new V3(-t, 0, s)]
    export const ISOCAHEDRON_FACE_VERTICES = [
        // 5 faces around point 0
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],

        // 5 adjacent faces
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],

        // 5 faces around point 3
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],

        // 5 adjacent faces
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1]]

    export function dodecahedron() {
        return makePlatonic(DODECAHEDRON_VERTICES, DODECAHEDRON_FACE_VERTICES, 'B2T.dodecahedron()')
    }
    export function octahedron() {
        return makePlatonic(OCTAHEDRON_VERTICES, OCTAHEDRON_FACE_VERTICES, 'B2T.octahedron()')
    }
    export function isocahedron() {
        return makePlatonic(ISOCAHEDRON_VERTICES, ISOCAHEDRON_FACE_VERTICES, 'B2T.octahedron()')
    }
	function makePlatonic(VS: V3[], FVIS: int[][], generator: string) {
	    const edgeMap = new Map()
        const faces = FVIS.map(faceIndexes => {
            const surface = PlaneSurface.throughPoints(VS[faceIndexes[0]], VS[faceIndexes[1]], VS[faceIndexes[2]])
            const contour = NLA.arrayFromFunction(faceIndexes.length, i => {
                const ipp = (i + 1) % faceIndexes.length
                const iA = faceIndexes[i], iB = faceIndexes[ipp]
                const iMin = min(iA, iB), iMax = max(iA, iB), edgeID = iMin * VS.length + iMax
                let edge = edgeMap.get(edgeID)
                !edge && edgeMap.set(edgeID, edge = StraightEdge.throughPoints(VS[iMin], VS[iMax]))
                return iA < iB ? edge : edge.flipped()
            })
            return new PlaneFace(surface, contour)
        })
        return new B2(faces, false, generator)
    }

	export function pyramidEdges(baseEdges: Edge[], apex: V3, name: string = 'pyramid' + globalId++): B2 {
		assertInst.apply(undefined, [Edge].concat(baseEdges))
		assertVectors(apex)

		const ribs = baseEdges.map(baseEdge => StraightEdge.throughPoints(apex, baseEdge.a))
		const faces = baseEdges.map((baseEdge, i) => {
			const faceName = name + 'Wall' + i
			const ipp = (i + 1) % baseEdges.length
			const faceEdges = [ribs[i], baseEdge, ribs[ipp].flipped()]
			const surface = undefined // TODO
			return Face.create(surface, faceEdges, undefined, faceName)
		})
		faces.push(bottomFace)
		const generator = makeGen('B2T.pyramidEdges', baseEdges, apex, name)
		return new B2(faces, false, generator, name)
	}
}

function makeGen(name: string, ...params: any[]) {
	return name + '(' + params.map(SCE).join(',') + ')'
}