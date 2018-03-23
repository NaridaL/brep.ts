import * as opentype from 'opentype.js'
import {
	arrayFromFunction,
	assert,
	assertf,
	assertInst,
	assertNumbers,
	assertVectors,
	callsce,
	eq,
	eq0,
	GOLDEN_RATIO,
	int,
	le,
	lerp,
	lt,
	M4,
	MINUS,
	raddd,
	snap,
	TAU,
	V,
	V3,
} from 'ts3dutils'

import {
	BezierCurve,
	BRep,
	ConicSurface,
	Curve,
	Edge,
	Face,
	FaceInfoFactory,
	getGlobalId,
	L3,
	P3,
	PCurveEdge,
	PlaneFace,
	PlaneSurface,
	ProjectedCurveSurface,
	RotatedCurveSurface,
	SemiCylinderSurface,
	SemiEllipseCurve,
	SemiEllipsoidSurface,
	StraightEdge,
	Surface,
	XiEtaCurve,
} from './index'

import { max, min, PI } from './math'

/**
 * Create a surface by projecting a curve in a direction.
 *
 * @param curve The curve to project.
 * @param offset The direction and distance to project curve.
 * @param flipped Whether the surface's default orientation (normal = curve tangent cross offset) should be flipped.
 */
export function projectCurve(curve: Curve, offset: V3, flipped: boolean): Surface {
	if (curve instanceof L3) {
		const surfaceNormal = offset.cross(curve.dir1).toLength(flipped ? -1 : 1)
		return new PlaneSurface(P3.normalOnAnchor(surfaceNormal, curve.anchor))
	}
	if (curve instanceof SemiEllipseCurve) {
		const curveDir = flipped ? offset : offset.negated()
		return new SemiCylinderSurface(curve, curveDir.unit(), undefined, undefined)
	}
	if (curve instanceof BezierCurve || curve instanceof XiEtaCurve) {
		const curveDir = offset.times(flipped ? 1 : -1)
		return new ProjectedCurveSurface(curve, curveDir, undefined, undefined, flipped ? 0 : -1, flipped ? 1 : 0)
	}
	throw new Error()
}

/**
 * Create a surface by projecting a curve onto a point.
 */
export function projectPointCurve(
	curve: Curve,
	tMin = curve.tMin,
	tMax = curve.tMax,
	p: V3,
	flipped: boolean,
): Surface {
	if (curve instanceof L3) {
		const up = curve.anchor.to(p).rejectedFrom(curve.dir1)
		return PlaneSurface.forAnchorAndPlaneVectors(curve.anchor, curve.dir1, up.unit(), tMin, tMax, 0, up.length())
	} else if (curve instanceof SemiEllipseCurve) {
		// flip f2 by default
		const factor = -1 * (flipped ? -1 : 1)
		return new ConicSurface(p, curve.f1.times(factor), curve.f2, p.to(curve.center), tMin, tMax, 0, 1)
	} else {
		throw new Error('projectPointCurve not implemented for ' + curve.constructor.name)
	}
}

/**
 * Create a surface by rotating a curve in the XZ-plane, with X > 0, around the Z-axis according to the right-hand rule.
 * @param curve The curve to rotate.
 * @param rotationAxis The line around which to rotate the curve.
 * @param flipped Whether the surface's default orientation (normal = curve tangent cross rotation tangent) should be
 * flipped.
 */
export function rotateCurve(
	curve: Curve,
	tMin = curve.tMin,
	tMax = curve.tMax,
	degrees: raddd,
	flipped: boolean,
): Surface {
	assertf(() => new PlaneSurface(P3.ZX).containsCurve(curve))
	if (curve instanceof L3) {
		if (curve.dir1.isParallelTo(V3.Z)) {
			if (eq0(curve.anchor.x)) {
				return undefined
			}
			const baseEllipse = new SemiEllipseCurve(
				V3.O,
				curve.anchor.xy(),
				curve.anchor.xy().getPerpendicular(),
				0,
				degrees,
			)
			// if curve.dir1 is going up (+Z), it the cylinder surface should face inwards
			const factor = (curve.dir1.z > 0 ? -1 : 1) * (flipped ? -1 : 1)
			const [zMin, zMax] = [curve.at(tMin).z * factor, curve.at(tMax).z * factor].sort(MINUS)
			return new SemiCylinderSurface(baseEllipse, V3.Z.times(factor), 0, degrees, zMin, zMax)
		}
		if (
			curve
				.at(tMin)
				.xy()
				.dot(curve.dir1) *
				curve
					.at(tMax)
					.xy()
					.dot(curve.dir1) <
			0
		) {
			throw new Error(
				'line cannot cross the Z axis in the [tMin, tMax] interval, as conic surfaces cannot have an hourglass shape.',
			)
		}
		if (curve.dir1.isPerpendicularTo(V3.Z)) {
			// if line.dir1 is pointing aways from V3.Z, then the surface should face up
			const factor = (curve.at(lerp(tMin, tMax, 0.5)).dot(curve.dir1) > 0 ? 1 : -1) * (flipped ? -1 : 1)
			return new PlaneSurface(new P3(V3.Z.times(factor), curve.anchor.z * factor))
		} else {
			// apex is intersection of segment with Z-axis
			const a = curve.at(tMin),
				b = curve.at(tMax)
			const apexZ = a.z - a.x * (b.z - a.z) / (b.x - a.x)
			const apex = new V3(0, 0, apexZ)
			const factor = -(a.x > b.x ? -1 : 1) * (flipped ? -1 : 1)
			const s = new ConicSurface(
				apex,
				new V3(curve.dir1.lengthXY(), 0, 0),
				new V3(0, curve.dir1.lengthXY(), 0),
				new V3(0, 0, (a.x > b.x ? -1 : 1) * curve.dir1.z),
				0,
				degrees,
				0,
				1,
			)
			return factor > 0 ? s : s.flipped()
		}
	}
	if (curve instanceof SemiEllipseCurve) {
		const a = curve.at(tMin),
			b = curve.at(tMax)
		const ell = curve.rightAngled()
		const f1Perp = ell.f1.isPerpendicularTo(V3.Z),
			f2Perp = ell.f2.isPerpendicularTo(V3.Z)
		if (L3.Z.containsPoint(ell.center) && (f1Perp || f2Perp)) {
			flipped = flipped == a.z > b.z
			let width = ell.f1.length(),
				height = ell.f2.length()
			if (ell.f1.isParallelTo(V3.Z)) {
				;[width, height] = [height, width]
			}
			return SemiEllipsoidSurface.forABC(width, (!flipped ? 1 : -1) * width, height, ell.center)
		} else {
			const s = new RotatedCurveSurface(curve, M4.IDENTITY, tMin, tMax)
			return s
		}
	}
	throw new Error()
}

export namespace B2T {
	/**
	 * Create a [BRep] of an axis-aligned box width starting at the origin and extending into +XYZ space.
	 * @param width x-direction size.
	 * @param height y-direction size.
	 * @param depth z-direction size.
	 * @param name
	 */
	export function box(
		width: number = 1,
		height: number = 1,
		depth: number = 1,
		name: string = 'box' + getGlobalId(),
	): BRep {
		assertNumbers(width, height, depth)
		assert('string' === typeof name)
		const baseVertices = [new V3(0, 0, 0), new V3(0, height, 0), new V3(width, height, 0), new V3(width, 0, 0)]
		const generator = callsce('B2T.box', width, height, depth, name)
		return B2T.extrudeVertices(baseVertices, P3.XY.flipped(), new V3(0, 0, depth), name, generator)
	}

	export function puckman(
		radius: number,
		rads: raddd,
		height: number,
		name: string = 'puckman' + getGlobalId(),
	): BRep {
		assertf(() => lt(0, radius))
		assertf(() => lt(0, rads) && le(rads, TAU))
		assertf(() => lt(0, height))
		const edges = StraightEdge.chain(
			[V3.O, new V3(radius, 0, 0), new V3(radius, 0, height), new V3(0, 0, height)],
			true,
		)
		return B2T.rotateEdges(edges, rads, name)
	}

	export function registerVertexName(map: Map<V3, string>, name: string, p: V3) {
		// TODO
		if (!Array.from(map.keys()).some(p2 => p2.like(p))) {
			map.set(p, name)
		}
	}

	/**
	 * Create a [BRep] by projecting a number of edges in a direction.
	 * @param baseFaceEdges
	 * @param baseFacePlane
	 * @param offset
	 * @param name
	 * @param gen
	 * @param infoFactory
	 */
	export function extrudeEdges(
		baseFaceEdges: Edge[],
		baseFacePlane: P3 = P3.XY,
		offset: V3 = V3.Z,
		name: string = 'extrude' + getGlobalId(),
		gen?: string,
		infoFactory?: FaceInfoFactory<any>,
	): BRep {
		baseFaceEdges = fixEdges(baseFaceEdges)
		//Array.from(combinations(baseFaceEdges.length)).forEach(({i, j}) => {
		//	assertf(() => !Edge.edgesIntersect(baseFaceEdges[i], baseFaceEdges[j]), baseFaceEdges[i].sce +
		// baseFaceEdges[j].sce) })
		assertf(() => Edge.isLoop(baseFaceEdges))
		// TODO checks..
		//if (offset.dot(baseFacePlane.normal1) > 0) {
		//	baseFacePlane = baseFacePlane.flipped()
		//}
		const vertexNames = new Map()
		const basePlaneSurface = new PlaneSurface(baseFacePlane)
		//assert(basePlaneSurface.edgeLoopCCW(baseFaceEdges), 'edges not CCW on baseFacePlane')
		const translationMatrix = M4.translate(offset)
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

		const ribs = arrayFromFunction(edgeCount, i =>
			StraightEdge.throughPoints(baseFaceEdges[i].a, topEdges[i].a, name + 'Rib' + i),
		)

		const faces = baseFaceEdges.map((edge, i) => {
			const faceName = name + 'Wall' + i
			const j = (i + 1) % edgeCount
			const faceEdges = [baseFaceEdges[i].flipped(), ribs[i], topEdges[i], ribs[j].flipped()]
			const surface = projectCurve(edge.curve, offset, edge.reversed)
			const info = infoFactory && infoFactory.extrudeWall(i, surface, faceEdges)
			return Face.create(surface, faceEdges, undefined, faceName, info)
		}) as Face[]
		faces.push(bottomFace, topFace)
		gen = gen || callsce('B2T.extrudeEdges', baseFaceEdges, baseFacePlane, offset, name)
		return new BRep(faces, baseFacePlane.normal1.dot(offset) > 0, gen, vertexNames)
	}

	export function cylinder(
		radius: number = 1,
		height: number = 1,
		rads: raddd = TAU,
		name: string = 'cylinder' + getGlobalId(),
	): BRep {
		const vertices = [new V3(0, 0, 0), new V3(radius, 0, 0), new V3(radius, 0, height), new V3(0, 0, height)]
		return rotateEdges(StraightEdge.chain(vertices, true), rads, name)
	}

	export function cone(
		radius: number = 1,
		height: number = 1,
		rads: raddd = TAU,
		name: string = 'cone' + getGlobalId(),
	): BRep {
		const vertices = [new V3(0, 0, 0), new V3(radius, 0, height), new V3(0, 0, height)]
		return rotateEdges(StraightEdge.chain(vertices, true), rads, name)
	}

	export function sphere(radius: number = 1, name: string = 'sphere' + getGlobalId(), rot: raddd = TAU): BRep {
		const ee = PCurveEdge.create(
			new SemiEllipseCurve(V3.O, new V3(0, 0, -radius), new V3(radius, 0, 0)),
			new V3(0, 0, -radius),
			new V3(0, 0, radius),
			0,
			PI,
			undefined,
			new V3(radius, 0, 0),
			new V3(-radius, 0, 0),
		)
		const generator = callsce('B2T.sphere', radius, name, rot)
		return rotateEdges([StraightEdge.throughPoints(ee.b, ee.a), ee], rot, name, generator)
	}

	/**
	 * Create a [[BRep]] of a menger sponge.
	 * @param res 0: just a cube, 1: every cube face has one hole, 2: 9 holes, etc
	 * @param name
	 */
	export function menger(res: int = 2, name: string = 'menger' + getGlobalId()): BRep {
		let result = B2T.box(1, 1, 1)
		if (0 == res) return result
		const punch = B2T.box(1 / 3, 1 / 3, 2)
			.translate(1 / 3, 1 / 3, -1 / 2)
			.flipped()

		function recurse(steps: int, m4: M4) {
			result = result.and(punch.transform(m4))
			if (steps > 1) {
				const scaled = m4.times(M4.scale(1 / 3, 1 / 3, 1))
				for (let i = 0; i < 9; i++) {
					if (4 == i) continue
					recurse(steps - 1, scaled.times(M4.translate(i % 3, (i / 3) | 0, 0)))
				}
			}
		}

		recurse(res, M4.IDENTITY)
		recurse(res, M4.YZX)
		recurse(res, M4.ZXY)
		return result
	}

	export function menger2(res: int = 2, name: string = 'menger' + getGlobalId()): BRep {
		if (0 == res) return B2T.box(1, 1, 1)

		const punch = B2T.box(1 / 3, 1 / 3, 2)
			.translate(1 / 3, 1 / 3, -1 / 2)
			.flipped()
		const stencilFaces: Face[] = []

		function recurse(steps: int, m4: M4) {
			stencilFaces.push(...punch.transform(m4).faces)
			if (steps > 1) {
				const scaled = m4.times(M4.scale(1 / 3, 1 / 3, 1))
				for (let i = 0; i < 9; i++) {
					if (4 == i) continue
					recurse(steps - 1, scaled.times(M4.translate(i % 3, (i / 3) | 0, 0)))
				}
			}
		}

		recurse(res, M4.IDENTITY)
		const stencil = new BRep(stencilFaces, true)

		return B2T.box()
			.and(stencil)
			.and(stencil.transform(M4.YZX))
			.and(stencil.transform(M4.ZXY))
	}

	/**
	 * Create a [BRep] of a torus.
	 * @param rSmall The radius to the surface of the torus.
	 * @param rLarge The radius from the origin to the inside of the torus.
	 * @param rads
	 * @param name
	 */
	export function torus(
		rSmall: number,
		rLarge: number,
		rads: raddd = TAU,
		name: string = 'torus' + getGlobalId(),
	): BRep {
		assertNumbers(rSmall, rLarge, rads)
		assertf(() => rLarge > rSmall)
		const curves = [
			SemiEllipseCurve.semicircle(rSmall, new V3(rLarge, 0, 0)),
			SemiEllipseCurve.semicircle(-rSmall, new V3(rLarge, 0, 0)),
		]
		const baseEdges = curves.map(c => PCurveEdge.forCurveAndTs(c, 0, Math.PI).rotateX(PI / 2))
		return B2T.rotateEdges(baseEdges, rads, name)
	}

	/**
	 * Create a [BRep] by smoothly rotating edges around Z.
	 * baseLoop should be CCW on XZ plane for a bounded BRep
	 */
	export function rotateEdges(
		baseLoop: Edge[],
		totalRads: raddd,
		name: string = 'rotateEdges' + getGlobalId(),
		generator?: string,
		infoFactory?: FaceInfoFactory<any>,
	): BRep {
		assert(baseLoop.every(e => new PlaneSurface(P3.ZX).containsCurve(e.curve)))
		assert(!eq(PI, totalRads) || PI == totalRads) // URHGJ
		assertf(() => lt(0, totalRads) && le(totalRads, TAU))
		totalRads = snap(totalRads, TAU)
		assertf(() => Edge.isLoop(baseLoop))
		const basePlane = new PlaneSurface(P3.ZX.flipped()).edgeLoopCCW(baseLoop)
			? new PlaneSurface(P3.ZX.flipped())
			: new PlaneSurface(P3.ZX)
		// const rotationSteps = ceil((totalRads - NLA_PRECISION) / PI)
		// const angles = rotationSteps == 1 ? [-PI, -PI + totalRads] : [-PI, 0, totalRads - PI]
		const open = !eq(totalRads, 2 * PI)
		const baseRibCurves = baseLoop.map(edge => {
			const a = edge.a,
				radius = a.lengthXY()
			if (!eq0(radius)) {
				return new SemiEllipseCurve(V(0, 0, a.z), V(radius, 0, 0), V(0, radius, 0))
			}
			return undefined
		})
		const baseSurfaces = baseLoop.map((edge, i) => {
			const s = rotateCurve(edge.curve, edge.minT, edge.maxT, PI, edge.deltaT() > 0)
			const t = lerp(edge.aT, edge.bT, 0.5)
			s &&
				assert(
					edge
						.tangentAt(t)
						.cross(V3.Y)
						.dot(s.normalP(edge.curve.at(t))) < 0,
				)
			return s
		})
		let stepStartEdges = baseLoop,
			stepEndEdges: Edge[]
		const faces = []
		for (let rot = 0; rot < totalRads; rot += PI) {
			const aT = 0,
				bT = min(totalRads - rot, PI)
			const rotation = M4.rotateZ(rot + bT)
			stepEndEdges = rot + bT == TAU ? baseLoop : baseLoop.map(edge => edge.transform(rotation))
			const ribs = arrayFromFunction(baseLoop.length, i => {
				const a = stepStartEdges[i].a,
					radius = a.lengthXY()
				const b = stepEndEdges[i].a
				if (!eq0(radius)) {
					const curve = 0 === rot ? baseRibCurves[i]! : baseRibCurves[i]!.rotateZ(rot)
					return new PCurveEdge(
						curve,
						a,
						b,
						aT,
						bT,
						undefined,
						curve.tangentAt(aT),
						curve.tangentAt(bT),
						name + 'rib' + i,
					)
				}
				return undefined
			})
			for (let edgeIndex = 0; edgeIndex < baseLoop.length; edgeIndex++) {
				if (baseSurfaces[edgeIndex]) {
					const edge = stepStartEdges[edgeIndex]
					const ipp = (edgeIndex + 1) % baseLoop.length
					const faceEdges = [
						stepStartEdges[edgeIndex].flipped(),
						!eq0(edge.a.x) && ribs[edgeIndex],
						stepEndEdges[edgeIndex],
						!eq0(edge.b.x) && ribs[ipp]!.flipped(),
					].filter((x: any): x is Edge => x)
					const surface = 0 === rot ? baseSurfaces[edgeIndex] : baseSurfaces[edgeIndex].rotateZ(rot)
					const info = infoFactory && infoFactory.extrudeWall(edgeIndex, surface, faceEdges, undefined)
					faces.push(Face.create(surface, faceEdges, undefined, name + 'Wall' + edgeIndex, info))
				}
			}
			stepStartEdges = stepEndEdges
		}
		if (open) {
			const endFaceEdges = Edge.reversePath(stepEndEdges!)
			const infoStart = infoFactory && infoFactory.rotationStart(basePlane, baseLoop, undefined)
			const infoEnd =
				infoFactory && infoFactory.rotationEnd(basePlane.flipped().rotateZ(totalRads), endFaceEdges, undefined)
			faces.push(
				new PlaneFace(basePlane, baseLoop, undefined, name + 'start', infoStart),
				new PlaneFace(basePlane.flipped().rotateZ(totalRads), endFaceEdges, undefined, name + 'end', infoEnd),
			)
		}
		const infiniteVolume = new PlaneSurface(P3.ZX).edgeLoopCCW(baseLoop)
		return new BRep(faces, infiniteVolume, generator)
	}

	/**
	 * loop should be CCW on XZ plane for a bounded BRep
	 */
	//export function rotateEdgesUnsplit(loop: Edge[], rads: raddd, name: string): BRep {
	//	assert(Edge.isLoop(loop))
	//	const rotationMatrix = M4.rotateZ(rads)
	//	const open = !eq(rads, 2 * PI)
	//	const endEdges = open ? loop.map(edge => edge.transform(rotationMatrix)) : loop
	//	const edgeCount = loop.length
	//	const ribs = arrayFromFunction(edgeCount, i => {
	//		const a = loop[i].a, radius = a.lengthXY()
	//		const b = endEdges[i].a
	//		if (!eq0(radius)) {
	//			const curve = new SemiEllipseCurve(V(0, 0, a.z), V(-radius, 0, 0), V(0, -radius, 0))
	//			const aT = -PI, bT = -PI + rads
	//			return new PCurveEdge(curve, a, b, aT, bT, undefined, curve.tangentAt(aT), curve.tangentAt(bT), name
	// + 'rib' + i) } }) const faces = loop.map((edge, i) => { const ipp = (i + 1) % edgeCount console.log('ljl', i,
	// ipp, ribs) const faceEdges = [ edge.flipped(), !eq0(edge.a.x) && ribs[i], endEdges[i], !eq0(edge.b.x) &&
	// ribs[ipp].flipped()].filter(x => x) if (edge instanceof StraightEdge) { const line = edge.curve let surface if
	// (line.dir1.isParallelTo(V3.Z)) { if (eq0(edge.a.x)) { return } let flipped = edge.a.z > edge.b.z surface = new
	// SemiCylinderSurface(ribs[i].curve, !flipped ? V3.Z : V3.Z.negated()) } else if
	// (line.dir1.isPerpendicularTo(V3.Z)) { let flipped = edge.a.x > edge.b.x let surface = new PlaneSurface(new
	// P3(V3.Z, edge.a.z)) if (!flipped) surface = surface.flipped() if (!open) { const hole = flipped ? !eq0(edge.b.x)
	// && ribs[ipp].flipped() : !eq0(edge.a.x) && ribs[i] return new PlaneFace(surface, [flipped ? ribs[i] :
	// ribs[ipp].flipped()], hole && [[hole]]) } return new PlaneFace(surface, faceEdges) } else { // apex is
	// intersection of segment with Z-axis let a = edge.a, b = edge.b let apexZ = a.z - a.x * (b.z - a.z) / (b.x - a.x)
	// let apex = new V3(0, 0, apexZ) let flipped = edge.a.z > edge.b.z surface =
	// ConicSurface.atApexThroughEllipse(apex, ribs[a.x > b.x ? i : ipp].curve as SemiEllipseCurve, !flipped ? 1 : -1)
	// } return Face.create(surface, faceEdges) } if (edge.curve instanceof SemiEllipseCurve) { let flipped = undefined
	// let ell = edge.curve.rightAngled() let f1Perp = ell.f1.isPerpendicularTo(V3.Z), f2Perp =
	// ell.f2.isPerpendicularTo(V3.Z) if (L3.Z.containsPoint(ell.center) && (f1Perp || f2Perp)) { let f3length = f1Perp
	// ? ell.f1.length() : ell.f2.length() if (flipped) { f3length *= -1 } let surface = new
	// SemiEllipsoidSurface(ell.center, ell.f1, ell.f2, ell.f1.cross(ell.f2).toLength(f3length)) return new
	// RotationFace(surface, faceEdges) } } else { assert(false, edge) } }).filter(x => x) if (open) { const
	// endFaceEdges = endEdges.map(edge => edge.flipped()).reverse() faces.push( new PlaneFace(new
	// PlaneSurface(P3.ZX.flipped()), loop), new PlaneFace(new PlaneSurface(P3.ZX.rotateZ(rads)), endFaceEdges)) }
	// return new BRep(faces, undefined) }

	export function quaffle() {
		const baseK = B2T.sphere(1)
			.translate(0, 1.7)
			.flipped()
		//const baseK = B2T.box().scale(0.2).translate(0, 0.95).flipped()
		// const vs = B2T.DODECAHEDRON_VERTICES.concat(
		// B2T.DODECAHEDRON_FACE_VERTICES.map(fis => fis
		// .map(vi => B2T.DODECAHEDRON_VERTICES[vi])
		// .reduce((a,b) => a.plus(b), V3.O)
		// .unit()))
		const ss = new BRep(TETRAHEDRON_VERTICES.flatMap(v => baseK.rotateAB(V3.Y, v).faces), false)
		//return ss
		return B2T.sphere().and(ss)
	}

	export function extrudeFace(face: PlaneFace, dir: V3) {
		return new BRep(
			extrudeEdges(face.contour, face.surface.plane, dir)
				.faces.slice(0, -2)
				.concat(
					face,
					face.translate(dir.x, dir.y, dir.z).flipped(),
					face.holes.flatMap(hole =>
						extrudeEdges(hole, face.surface.plane.flipped(), dir).faces.slice(0, -2),
					),
				),
			false,
		)
	}

	export let defaultFont: opentype.Font

	export function loadFonts(): Promise<opentype.Font> {
		return loadFont('fonts/FiraSansMedium.woff').then(font => (defaultFont = font))
	}

	const loadedFonts = new Map<string, opentype.Font>()

	export function loadFont(fontPath: string): Promise<opentype.Font> {
		return new Promise<opentype.Font>(function(resolve, reject) {
			const font = loadedFonts.get(fontPath)
			if (font) {
				resolve(font)
			} else {
				opentype.load(fontPath, function(err, f) {
					if (err) {
						reject(err)
					} else {
						loadedFonts.set(fontPath, f!)
						resolve(f)
					}
				})
			}
		})
	}

	export function loadFontsAsync(callback: () => void) {
		if (defaultFont) {
			callback()
		} else {
			opentype.load('fonts/FiraSansMedium.woff', function(err, font) {
				if (err) {
					throw new Error('Could not load font: ' + err)
				} else {
					defaultFont = font!
					callback()
				}
			})
		}
	}

	/**
	 * Create the [BRep] of a string rendered in a font.
	 * @param text
	 * @param size
	 * @param depth
	 * @param font An opentype.js font.
	 */
	export function text(text: string, size: number, depth: number = 1, font: opentype.Font = defaultFont) {
		const path = font.getPath(text, 0, 0, size)
		const subpaths: opentype.PathCommand[][] = []
		path.commands.forEach(c => {
			if (c.type == 'M') {
				subpaths.push([])
			}
			subpaths.last.push(c)
		})
		const loops = subpaths.map(sp => {
			const path = new opentype.Path()
			path.commands = sp
			const loop = Edge.reversePath(Edge.pathFromSVG(path.toPathData(13))).map(e => e.mirrorY())
			assert(Edge.isLoop(loop))
			return loop
		})
		const faces = Face.assembleFacesFromLoops(loops, new PlaneSurface(P3.XY), PlaneFace as any)
		const generator = callsce('B2T.text', text, size, depth)
		return BRep.join(faces.map(face => B2T.extrudeFace(face as PlaneFace, V(0, 0, -depth))), generator)
	}

	export function minorityReport() {
		const a = B2T.sphere()
		const b = B2T.text('LEO CROW', 64, 128)
			.scale(0.1 / 32)
			.translate(-0.5, -0.05, 1.2)
			.flipped()
		const c = B2T.sphere(0.98)
		return a.and(b).plus(c)
	}

	export function whatever() {
		const iso = isocahedron()
		const numbersBRep = BRep.join(
			iso.faces.map((face, i) => {
				const numberBRep = text('' + (i + 1), 0.4, -2)
				const centroid = face.contour
					.map(edge => edge.a)
					.reduce((a, b) => a.plus(b), V3.O)
					.div(3)

				const sys = M4.forSys(
					face.contour[0].aDir,
					centroid.cross(face.contour[0].aDir),
					centroid.unit(),
					centroid,
				)
				return numberBRep.transform(sys.times(M4.translate(-numberBRep.getAABB().size().x / 2, -0.1, -0.04)))
			}),
		)
		const s = sphere(0.9)
		//return iso.and(numbersBRep)
		return iso.and(s).and(numbersBRep)
		//return numbersBRep
	}

	export function whatever3() {
		const t = B2T.torus(1, 2)
		return B2T.box(5, 5, 2)
			.translate(-2.5, -2.5)
			.minus(t)
	}

	export function d20() {
		const iso = isocahedron()
		const numbersBRep = BRep.join(
			iso.faces.map((face, i) => {
				const numberBRep = text('' + (i + 1), 0.4, -2)
				const centroid = face.contour
					.map(edge => edge.a)
					.reduce((a, b) => a.plus(b), V3.O)
					.div(3)

				const sys = M4.forSys(
					face.contour[0].aDir,
					centroid.cross(face.contour[0].aDir),
					centroid.unit(),
					centroid,
				)
				return numberBRep.transform(sys.times(M4.translate(-numberBRep.getAABB().size().x / 2, -0.1, -0.04)))
			}),
		)
		const s = sphere(0.9)
		//return iso.and(numbersBRep)
		return iso.and(s).and(numbersBRep)
		//return numbersBRep
	}

	/**
	 * Create a [BRep] by rotating a number of edges in steps in the ZX plane, with X > 0
	 * around the Z axis according to the right-hand rule. The edges from each rotation step
	 * ("ribs") are connected by projecting a rib onto the next. For example, line edges
	 * will always be connected by [PlaneSurface]s
	 *
	 * @example Roundabout way of creating a cube:
	 * ```ts
	 * B2T.rotStep(StraightEdge.rect(Math.sqrt(2) / 2, 1), 360 * DEG, 4)
	 * ```
	 * @param the edges to rotate.
	 * @param totalRads The angle between the original and last rib. [0; TAU]
	 * @param count The number of steps to take. If totalRads == TAU, it must be >= 3, otherwise >= 2.
	 */
	export function rotStep(edges: Edge[], totalRads: raddd, count: int): BRep
	export function rotStep(edges: Edge[], angles: raddd[]): BRep
	export function rotStep(edges: Edge[], totalRadsOrAngles: raddd | raddd[], countO?: int): BRep {
		const angles: number[] =
			'number' === typeof totalRadsOrAngles
				? arrayFromFunction(countO!, i => (i + 1) / countO! * totalRadsOrAngles)
				: totalRadsOrAngles
		const count = angles.length
		const open = !eq(TAU, angles.last)
		const ribs = [
			edges,
			...angles.map(phi => {
				if (eq(TAU, phi)) {
					return edges
				}
				const matrix = M4.rotateZ(phi)
				return edges.map(edge => edge.transform(matrix))
			}),
		]
		const horizontalEdges = arrayFromFunction(count, i => {
			const ipp = (i + 1) % (count + 1)
			return arrayFromFunction(edges.length, j => {
				if (!eq0(edges[j].a.lengthXY())) {
					return StraightEdge.throughPoints(ribs[i][j].a, ribs[ipp][j].a)
				}
				return undefined!
			})
		})

		const faces: Face[] = []
		let face: Face
		edges.forEach((edge, i) => {
			const ipp = (i + 1) % edges.length
			// for straight edges perpendicular to the Z-axis, we only create one face.
			if (edge instanceof StraightEdge && edge.curve.dir1.isPerpendicularTo(V3.Z)) {
				const flipped = edge.a.x > edge.b.x
				const surface = new PlaneSurface(flipped ? new P3(V3.Z, edge.a.z) : new P3(V3.Z.negated(), -edge.a.z))
				if (open) {
					const faceEdges: Edge[] = []
					if (!eq0(edge.a.x)) {
						faceEdges.push(...arrayFromFunction(count, j => horizontalEdges[j][i]!))
					}
					faceEdges.push(ribs[count][i])
					if (!eq0(edge.b.x)) {
						faceEdges.push(...arrayFromFunction(count, j => horizontalEdges[count - j - 1][ipp]!.flipped()))
					}
					faceEdges.push(edge.flipped())
					face = new PlaneFace(surface, faceEdges)
				} else {
					const contour = flipped
						? arrayFromFunction(count, j => horizontalEdges[j][i])
						: arrayFromFunction(count, j => horizontalEdges[count - j - 1][ipp]!.flipped())
					let hole
					if (flipped && !eq0(edge.b.x)) {
						hole = arrayFromFunction(count, j => horizontalEdges[count - j - 1][ipp]!.flipped())
					} else if (!flipped && !eq0(edge.a.x)) {
						hole = arrayFromFunction(count, j => horizontalEdges[j][i])
					}
					face = new PlaneFace(surface, contour, hole ? [hole] : [])
				}
				faces.push(face)
				return
			} else if (edge instanceof StraightEdge) {
				if (eq0(edge.a.lengthXY()) && eq0(edge.b.lengthXY())) {
					return
				}
			}
			for (let r = 0; r < count; r++) {
				const rpp = (r + 1) % (count + 1)
				const faceEdges = [
					ribs[r][i].flipped(),
					horizontalEdges[r][i],
					ribs[rpp][i],
					horizontalEdges[r][ipp] && horizontalEdges[r][ipp]!.flipped(),
				].filter(x => x)
				let surface
				if (edge instanceof StraightEdge) {
					surface = new PlaneSurface(P3.throughPoints(faceEdges[0].a, faceEdges[1].a, faceEdges[2].a))
				} else {
					const maxX = edges[i].getAABB().max.x
					const phi = angles[r],
						prevPhi = 0 == r ? 0 : angles[r - 1]
					const offset = V3.polar(maxX, prevPhi).to(V3.polar(maxX, phi))
					surface = projectCurve(ribs[r][i].curve, offset, false)
				}
				faces.push(Face.create(surface, faceEdges))
			}
		})
		if (open) {
			const endFaceEdges = ribs[count].map(edge => edge.flipped()).reverse()
			const endFace = new PlaneFace(new PlaneSurface(P3.ZX.rotateZ(angles.last)), endFaceEdges)
			faces.push(new PlaneFace(new PlaneSurface(P3.ZX.flipped()), edges), endFace)
		}
		return new BRep(faces, new PlaneSurface(P3.ZX).edgeLoopCCW(edges))
	}

	export function fixEdges(edges: Edge[]): Edge[] {
		return edges.flatMap(edge => {
			const c = edge.curve
			if (c instanceof SemiEllipseCurve && c.tMin === -PI && c.tmax === PI) {
				const splitEdges = edge.minT < 0 && edge.maxT > 0 ? edge.split(0) : [edge]
				return splitEdges.map(edge => {
					if (edge.minT >= 0) {
						return Edge.create(
							new SemiEllipseCurve(c.center, c.f1, c.f2, max(0, c.tMin), c.tMax),
							edge.a,
							edge.b,
							edge.aT,
							edge.bT,
							undefined,
							edge.aDir,
							edge.bDir,
							edge.name,
						)
					} else {
						// "rotate" the curve
						return Edge.create(
							new SemiEllipseCurve(
								c.center,
								c.f1.negated(),
								c.f2.negated(),
								c.tMin + PI,
								min(PI, c.tMax + PI),
							),
							edge.a,
							edge.b,
							edge.aT + PI,
							edge.bT + PI,
							undefined,
							edge.aDir,
							edge.bDir,
							edge.name,
						)
					}
				})
			}
			if (c instanceof BezierCurve) {
				if (edge.a.like(edge.b)) {
					return edge.split(lerp(edge.aT, edge.bT, 0.5))
				}
			}
			return edge
		})
	}

	/**
	 * Create a [BRep] by projecting edges created by joining vertices with straight edges.
	 * @param baseVertices
	 * @param baseFacePlane
	 * @param offset
	 * @param name
	 * @param generator
	 */
	export function extrudeVertices(
		baseVertices: V3[],
		baseFacePlane: P3,
		offset: V3,
		name?: string,
		generator?: string,
	) {
		assert(baseVertices.every(v => v instanceof V3), 'baseVertices.every(v => v instanceof V3)')
		assertInst(P3, baseFacePlane)
		assertVectors(offset)
		if (baseFacePlane.normal1.dot(offset) > 0) baseFacePlane = baseFacePlane.flipped()
		const edges = StraightEdge.chain(baseVertices, true)
		generator = generator || callsce('B2T.extrudeVertices', baseVertices, baseFacePlane, offset, name)
		return B2T.extrudeEdges(edges, baseFacePlane, offset, name, generator)
	}

	/**
	 * Create a tetrahedron (3 sided pyramid) [BRep].
	 * `a`, `b`, `c` and `d` can be in any order. The only constraint is that they cannot be on a common plane.
	 * The resulting tetrahedron will always have outwards facing faces.
	 * @param a
	 * @param b
	 * @param c
	 * @param d
	 * @param name
	 */
	export function tetrahedron(a: V3, b: V3, c: V3, d: V3, name: string = 'tetra' + getGlobalId()): BRep {
		assertVectors(a, b, c, d)
		const dDistance = P3.throughPoints(a, b, c).distanceToPointSigned(d)
		if (eq0(dDistance)) {
			throw new Error('four points are coplanar')
		}
		if (dDistance > 0) {
			;[c, d] = [d, c]
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
			new PlaneFace(PlaneSurface.throughPoints(c, d, a), [cd, ad.flipped(), ac], [], name + 'cda'),
		]
		const gen = callsce('B2T.tetrahedron', a, b, c, d)
		return new BRep(faces, false, gen)
	}

	const b = 1 / GOLDEN_RATIO,
		c = 2 - GOLDEN_RATIO
	export const TETRAHEDRON_VERTICES = [
		new V3(1, 0, -1 / Math.sqrt(2)),
		new V3(-1, 0, -1 / Math.sqrt(2)),
		new V3(0, -1, 1 / Math.sqrt(2)),
		new V3(0, 1, 1 / Math.sqrt(2)),
	].map(v => v.unit())
	export const DODECAHEDRON_VERTICES = [
		new V3(c, 0, 1),
		new V3(-c, 0, 1),
		new V3(-b, b, b),
		new V3(0, 1, c),
		new V3(b, b, b),
		new V3(b, -b, b),
		new V3(0, -1, c),
		new V3(-b, -b, b),
		new V3(c, 0, -1),
		new V3(-c, 0, -1),
		new V3(-b, -b, -b),
		new V3(0, -1, -c),
		new V3(b, -b, -b),
		new V3(b, b, -b),
		new V3(0, 1, -c),
		new V3(-b, b, -b),
		new V3(1, c, 0),
		new V3(-1, c, 0),
		new V3(-1, -c, 0),
		new V3(1, -c, 0),
	].map(v => v.unit())
	export const DODECAHEDRON_FACE_VERTICES = [
		[4, 3, 2, 1, 0],
		[7, 6, 5, 0, 1],
		[12, 11, 10, 9, 8],
		[15, 14, 13, 8, 9],
		[14, 3, 4, 16, 13],
		[3, 14, 15, 17, 2],
		[11, 6, 7, 18, 10],
		[6, 11, 12, 19, 5],
		[4, 0, 5, 19, 16],
		[12, 8, 13, 16, 19],
		[15, 9, 10, 18, 17],
		[7, 1, 2, 17, 18],
	]

	export const OCTAHEDRON_VERTICES = [
		new V3(1, 0, 0),
		new V3(-1, 0, 0),
		new V3(0, 1, 0),
		new V3(0, -1, 0),
		new V3(0, 0, 1),
		new V3(0, 0, -1),
	]
	export const OCTAHEDRON_FACE_VERTICES = [
		[0, 2, 4],
		[2, 1, 4],
		[1, 3, 4],
		[3, 0, 4],

		[2, 0, 5],
		[1, 2, 5],
		[3, 1, 5],
		[0, 3, 5],
	]

	const { x: s, y: t } = new V3(1, GOLDEN_RATIO, 0).unit()
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
		new V3(-t, 0, s),
	]
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
		[9, 8, 1],
	]

	/**
	 * Create a dodecahedron [BRep]. The vertices are on the unit sphere.
	 */
	export function dodecahedron() {
		return makePlatonic(DODECAHEDRON_VERTICES, DODECAHEDRON_FACE_VERTICES, 'B2T.dodecahedron()')
	}

	/**
	 * Create an octahedron [BRep]. The vertices are on the unit sphere.
	 */
	export function octahedron() {
		return makePlatonic(OCTAHEDRON_VERTICES, OCTAHEDRON_FACE_VERTICES, 'B2T.octahedron()')
	}

	/**
	 * Create an isocahedron [BRep]. The vertices are on the unit sphere.
	 */
	export function isocahedron() {
		return makePlatonic(ISOCAHEDRON_VERTICES, ISOCAHEDRON_FACE_VERTICES, 'B2T.octahedron()')
	}

	function makePlatonic(VS: V3[], FVIS: int[][], generator: string) {
		const edgeMap = new Map()
		const faces = FVIS.map(faceIndexes => {
			const surface = PlaneSurface.throughPoints(VS[faceIndexes[0]], VS[faceIndexes[1]], VS[faceIndexes[2]])
			const contour = arrayFromFunction(faceIndexes.length, i => {
				const ipp = (i + 1) % faceIndexes.length
				const iA = faceIndexes[i],
					iB = faceIndexes[ipp]
				const iMin = min(iA, iB),
					iMax = max(iA, iB),
					edgeID = iMin * VS.length + iMax
				let edge = edgeMap.get(edgeID)
				!edge && edgeMap.set(edgeID, (edge = StraightEdge.throughPoints(VS[iMin], VS[iMax])))
				return iA < iB ? edge : edge.flipped()
			})
			return new PlaneFace(surface, contour)
		})
		return new BRep(faces, false, generator)
	}

	/**
	 * Create a [BRep] by projecting a number of edges onto a point.
	 * @param baseEdges The edges forming the base of the pyramid.
	 * @param apex The tip of the pyramid.
	 * @param name
	 */
	export function pyramidEdges(baseEdges: Edge[], apex: V3, name: string = 'pyramid' + getGlobalId()): BRep {
		assertInst(Edge, ...baseEdges)
		assertVectors(apex)

		const ribs = baseEdges.map(baseEdge => StraightEdge.throughPoints(apex, baseEdge.a))
		const faces = baseEdges.map((baseEdge, i) => {
			const faceName = name + 'Wall' + i
			const ipp = (i + 1) % baseEdges.length
			const faceEdges = [ribs[i], baseEdge, ribs[ipp].flipped()]
			const surface = projectPointCurve(baseEdge.curve, baseEdge.minT, baseEdge.maxT, apex, baseEdge.deltaT() < 0)
			return Face.create(surface, faceEdges, undefined, faceName)
		})
		const baseSurface = new PlaneSurface(P3.XY).flipped()
		const bottomFace = Face.create(baseSurface, Edge.reversePath(baseEdges))
		faces.push(bottomFace)
		const generator = callsce('B2T.pyramidEdges', baseEdges, apex, name)
		return new BRep(faces, false, generator)
	}
}
