
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

namespace B2 {
    export function box(w: number, h: number, d: number, name?: string): B2 {
        assertNumbers(w, h, d)
        assertInst('string' === typeof name)
        const baseVertices = [
	        V(0, 0, 0),
	        V(0, h, 0),
	        V(w, h, 0),
	        V(w, 0, 0)
        ];
        return B2.extrudeVertices(baseVertices, P3.XY.flipped(), V(0, 0, d), name, `B2.box(${w}, ${h}, ${d}, "${name || ''}")`)
    }

    export function puckman(radius: number, rads: number, height: number, name: string): B2 {
        // TODO: argument checking
        const circleCurve = new EllipseCurve(V3.ZERO, V(radius, 0, 0), V(0, -radius, 0));
        const a = circleCurve.at(0);
        const b = circleCurve.at(-rads);
        const edges = [
	        StraightEdge.throughPoints(a, V3.ZERO),
	        StraightEdge.throughPoints(V3.ZERO, b),
	        new PCurveEdge(circleCurve, b, a, -rads, 0, null, circleCurve.tangentAt(-rads), circleCurve.tangentAt(0))];
        return B2.extrudeEdges(edges, P3.XY.flipped(), V(0, 0, height), name)
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
        const translationMatrix = M4.translation(offset);
        const topEdges = baseFaceEdges.map(edge => edge.transform(translationMatrix, 'top'));
        const edgeCount = baseFaceEdges.length;
        const bottomFace = new PlaneFace(basePlaneSurface, baseFaceEdges, [], name + 'Bottom');
        const topFaceEdges = topEdges.map(edge => edge.flipped()).reverse();
        const topFace = new PlaneFace(new PlaneSurface(baseFacePlane.flipped().translated(offset)), topFaceEdges, [], name + 'Top');


        baseFaceEdges.forEach(edge => B2.registerVertexName(vertexNames, edge.name + 'A', edge.a))
        topFaceEdges.forEach(edge => B2.registerVertexName(vertexNames, edge.name + 'A', edge.a))

        const ribs = NLA.arrayFromFunction(edgeCount,
	        i => StraightEdge.throughPoints(baseFaceEdges[i].a, topEdges[i].a, name + 'Rib' + i));

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


    export function cylinder(radius: number, height: number, rads: number, name: string): B2 {
	    const vertices = [new V3(0, 0, 0), new V3(radius, 0, 0), new V3(radius, 0, height), new V3(0, 0, height)]
	    return B2.rotateEdges(StraightEdge.chain(vertices, true), rads || 2 * PI, name)
    }

    export function torus(rSmall: number, rLarge: number, rads: number, name: string): B2 {
        assertNumbers(rSmall, rLarge, rads)
        assertf(() => rLarge > rSmall)
	    let baseEdge = PCurveEdge.forCurveAndTs(EllipseCurve.circle(rSmall, new V3(rLarge, 0, 0)), -Math.PI, Math.PI)
        return B2.rotateEdges([baseEdge], rads, name || 'torus' + globalId++)
    }

    export function rotateEdges(edges: Edge[], rads: number, name: string): B2 {
        const rotationMatrix = M4.rotationZ(rads);
        let open = !NLA.eq(rads, 2 * PI);
        const endEdges = open ? edges.map(edge => edge.transform(rotationMatrix)) : edges;
        const edgeCount = edges.length;
        const ribs = NLA.arrayFromFunction(edgeCount, i => {
	        const a = edges[i].a, radius = a.lengthXY();
	        const b = endEdges[i].a;
	        if (!NLA.eq0(radius)) {
		        const curve = new EllipseCurve(V(0, 0, a.z), V(-radius, 0, 0), V(0, -radius, 0));
		        const aT = -PI, bT = -PI + rads;
		        return new PCurveEdge(curve, a, b, aT, bT, null, curve.tangentAt(aT), curve.tangentAt(bT), name + 'rib' + i)
	        }
        });
        const faces = edges.map((edge, i) => {
	        const ipp = (i + 1) % edgeCount;
	        const faceEdges = [
		        edge.flipped(),
		        !NLA.eq0(edge.a.x) && ribs[i],
		        endEdges[i],
		        !NLA.eq0(edge.b.x) && ribs[ipp].flipped()].filter(x => x);
	        let curve = edge.curve;
	        if (edge instanceof StraightEdge) {
		        const line = edge.curve;
		        if (line.dir1.isParallelTo(V3.Z)) {
			        if (NLA.eq0(edge.a.x)) {
				        return
			        }
			        let flipped = edge.a.z > edge.b.z
			        let surface = new CylinderSurface(ribs[i].curve, !flipped ? V3.Z : V3.Z.negated())
			        return new RotationFace(surface, faceEdges)
		        } else if (line.dir1.isPerpendicularTo(V3.Z)) {
			        let flipped = edge.a.x > edge.b.x
			        let surface = new PlaneSurface(new P3(V3.Z, edge.a.z))
			        if (!flipped) surface = surface.flipped()
			        if (!open) {
				        const hole = flipped
					        ? !NLA.eq0(edge.b.x) && ribs[ipp].flipped()
					        : !NLA.eq0(edge.a.x) && ribs[i];
				        return new PlaneFace(surface, [flipped ? ribs[i] : ribs[ipp].flipped()], hole && [[hole]])
			        }
			        return new PlaneFace(surface, faceEdges)
		        } else {
			        // apex is intersection of segment with Z-axis
			        let a = edge.a, b = edge.b
			        let apexZ = a.z - a.x * (b.z - a.z) / (b.x - a.x)
			        let apex = new V3(0, 0, apexZ)
			        let flipped = edge.a.z > edge.b.z
			        let surface = ConicSurface.atApexThroughEllipse(apex, ribs[a.x > b.x ? i : ipp].curve as EllipseCurve, !flipped ? 1 : -1)
			        return new RotationFace(surface, faceEdges)
		        }
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
        }).filter(x =>x);
        if (open) {
            const endFaceEdges = endEdges.map(edge => edge.flipped()).reverse();
            const endFace = new PlaneFace(new PlaneSurface(P3.ZX.rotateZ(rads)), endFaceEdges);
            faces.push(new PlaneFace(new PlaneSurface(P3.ZX.flipped()), edges), endFace)
        }
        return new B2(faces)
    }

    export function rotStep(edges: Edge[], totalRads: number, count: int) {
        const radStep = totalRads / count;
        let open = !NLA.eq(totalRads, 2 * PI);
        const ribCount = !open ? count : count + 1;
        const ribs = NLA.arrayFromFunction(ribCount, i => {
	        if (i == 0) return edges
	        const matrix = M4.rotationZ(radStep * i);
	        return edges.map(edge => edge.transform(matrix))
        });
        const hs = NLA.arrayFromFunction(count, i => {
	        const ipp = (i + 1) % ribCount;
	        return NLA.arrayFromFunction(edges.length, j => {
		        if (!NLA.eq0(edges[j].a.lengthXY())) {
			        return StraightEdge.throughPoints(ribs[i][j].a, ribs[ipp][j].a)
		        }
	        })
        });
        const faces = []
	    let surface, face
        edges.forEach((edge, i) => {
            const ipp = (i + 1) % edges.length;
            if (edge instanceof StraightEdge && edge.curve.dir1.isPerpendicularTo(V3.Z)) {
                let flipped = edge.a.x > edge.b.x;
                surface = new PlaneSurface(flipped ? new P3(V3.Z, edge.a.z) : new P3(V3.Z.negated(), -edge.a.z))
                if (open) {
                    const newEdges = [];
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
	                    : NLA.arrayFromFunction(count, j => hs[count - j - 1][ipp].flipped());
                    let hole;
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
                const rpp = (r + 1) % ribCount;
                const faceEdges = [ribs[r][i].flipped(), hs[r][i], ribs[rpp][i], hs[r][ipp] && hs[r][ipp].flipped()].filter(x => x);
                if (edge instanceof StraightEdge) {
                    var surface = new PlaneSurface(P3.throughPoints(faceEdges[0].a, faceEdges[1].a, faceEdges[2].a))
                    faces.push(new PlaneFace(surface, faceEdges))
                } else {
                    assert(false, edge.toString())
                }
            }
        })
        if (open) {
            const endFaceEdges = ribs[count].map(edge => edge.flipped()).reverse();
            const endFace = new PlaneFace(new PlaneSurface(P3.XZ.rotateZ(totalRads)), endFaceEdges);
            faces.push(new PlaneFace(new PlaneSurface(P3.XZ.flipped()), edges), endFace)
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
	    source = source || `B2.extrudeVertices(${baseVertices.sce}, ${baseFacePlane.sce}, ${offset.sce}, "${name}")`
        return B2.extrudeEdges(edges, baseFacePlane, offset, name, source)
    }

    // Returns a tetrahedron (3 sided pyramid).
    // Faces will face outwards.
    // abcd can be in any order. The only constraint is that abcd cannot be on a common plane.
    export function tetrahedron(a: V3, b: V3, c: V3, d: V3, name:string = 'tetra' + globalId++): B2 {
        assertVectors(a, b, c, d)
        let dDistance = P3.throughPoints(a, b, c).distanceToPointSigned(d)
        if (NLA.eq0(dDistance)) {
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
            new PlaneFace(PlaneSurface.throughPoints(a, b, c), [ab, bc, ac.flipped()], [], name + 'abc'),
            new PlaneFace(PlaneSurface.throughPoints(a, d, b), [ad, bd.flipped(), ab.flipped()], [], name + 'adb'),
            new PlaneFace(PlaneSurface.throughPoints(b, d, c), [bd, cd.flipped(), bc.flipped()], [], name + 'bdc'),
            new PlaneFace(PlaneSurface.throughPoints(c, d, a), [cd, ad.flipped(), ac], [], name + 'cda')
        ]
        let gen = `B2.tetrahedron(${a.sce}, ${b.sce}, ${c.sce}, ${d.sce})`
        return new B2(faces, false, gen)
    }

	export function pyramidEdges(edges: Edge[], apex: V3, name: string = 'pyramid' + globalId++): B2 {

	}
}