import { Assert, outputLink, skip, suite, test, testBRepAnd, testBRepOp } from './manager'

import { JavaMap as CustomMap, JavaSet as CustomSet } from 'javasetmap.ts'
import { DEG, M4, V, V3 } from 'ts3dutils'
import {
	B2T,
	BezierCurve,
	BRep,
	Edge,
	Face,
	IntersectionPointInfo,
	L3,
	P3,
	PCurveEdge,
	PlaneFace,
	PlaneSurface,
	ProjectedCurveSurface,
	RotationFace,
	SemiCylinderSurface,
	SemiEllipseCurve,
	StraightEdge,
} from '..'

function doTest(assert: Assert, face: PlaneFace, brep2: BRep, resultEdges: Edge[], resultPoints: V3[], desc?: string) {
	if (brep2 instanceof Face) {
		brep2 = new BRep([brep2], false)
	}
	brep2.buildAdjacencies()
	outputLink(
		{
			a: new BRep([face]),
			b: brep2,
			edges: resultEdges,
			drPs: resultPoints,
		},
		desc,
	)
	const faceMap = new Map(),
		edgeMap = new Map(),
		colinearEdgePairs = new CustomSet()
	brep2.faces.forEach(face2 => {
		face.intersectPlaneFace(face2, new BRep([face]), brep2, faceMap, edgeMap, colinearEdgePairs)
	})
	const edges = faceMap.get(face) || []
	console.log(faceMap)
	assert.equal(
		edges.length,
		resultEdges.length,
		resultEdges.length + ' == resultEdges.length == edges.length' + edges.toSource(),
	)
	resultEdges.forEach(edge => {
		assert.ok(
			edges.some(edge2 => edge.like(edge2)),
			`edges.some(edge2 => edge.like(edge2)) [${edges.toSource()}] ${edge.toSource()}`,
		)
	})
	const uniquePoints = []
	face.contour.forEach(edge => {
		const em = edgeMap.get(edge)
		em &&
			em.forEach(
				info =>
					info && !uniquePoints.some(up => up.like(info.p)) && assert.ok(info.p) && uniquePoints.push(info.p),
			)
	})
	assert.equal(
		uniquePoints.length,
		resultPoints.length,
		resultPoints.length + ' == points.length == resultPoints.length' + uniquePoints.toSource(),
	)
	resultPoints.forEach(p => {
		assert.ok(
			uniquePoints.some(up => up.like(p)),
			`edges.some(edge2 => edge.like(edge2)) [${uniquePoints.toSource()}]`,
		)
	})
}

function doTestWithBrep(
	assert: Assert,
	face: Face,
	faceBrep: BRep,
	brep2: BRep,
	expectedEdges: Edge[],
	expectedPoints: V3[],
	desc: string,
	backwards?: boolean,
) {
	return
	faceBrep.buildAdjacencies()
	brep2.buildAdjacencies()
	outputLink(
		assert,
		{
			a: faceBrep,
			b: brep2,
			edges: expectedEdges,
			drPs: expectedPoints,
		},
		`expected ${desc}`,
	)
	const faceMap = new Map(),
		faceEdgePoints: Map<Edge, IntersectionPointInfo[]> = new CustomMap(),
		checkedPairs = new CustomSet()
	let x
	if (!backwards) {
		brep2.faces.forEach(face2 => {
			face.intersectFace(face2, faceBrep, brep2, faceMap, faceEdgePoints, (x = new CustomMap()), checkedPairs)
		})
	} else {
		brep2.faces.forEach(face2 => {
			face2.intersectFace(face, brep2, faceBrep, faceMap, new CustomMap(), faceEdgePoints, checkedPairs)
		})
	}
	const edges = faceMap.get(face) || []
	const edgePointInfos = face
		.getAllEdges()
		.mapFilter(e => {
			const canonEdgePoints = faceEdgePoints.get(e.getCanon())
			return canonEdgePoints
		})
		.concatenated()
	const uniquePoints = []

	edgePointInfos &&
		edgePointInfos.forEach(
			info => !uniquePoints.some(up => up.like(info.p)) && assert.ok(info.p) && uniquePoints.push(info.p),
		)
	console.log('edgePointInfos', edgePointInfos)
	outputLink(
		assert,
		{
			a: faceBrep,
			b: brep2,
			edges: edges,
			drPs: uniquePoints,
		},
		`actual ${desc}`,
	)
	console.log(faceMap)
	assert.equal(
		edges.length,
		expectedEdges.length,
		expectedEdges.length + ' == resultEdges.length == edges.length' + edges.toSource(),
	)
	expectedEdges.forEach(edge => {
		assert.ok(
			edges.some(edge2 => edge.like(edge2)),
			`edges.some(edge2 => edge.like(edge2)) [${edges.toSource()}] ${edge.toSource()}`,
		)
	})
	assert.equal(
		uniquePoints.length,
		expectedPoints.length,
		expectedPoints.length + ' == resultPoints.length == uniquePoints.length' + uniquePoints.toSource(),
	)
	expectedPoints.forEach(p => {
		assert.ok(
			uniquePoints.some(up => up.like(p)),
			`uniquePoints.some(up => up.like(p)) [${uniquePoints.toSource()}]`,
		)
	})
}

function doTest2(assert: Assert, face: Face, brep: BRep, resultFaces: Face[], desc: string) {
	if (brep instanceof Face) {
		brep = new BRep([brep])
	}
	const faceMap = new Map(),
		edgeMap = new CustomMap<Edge, IntersectionPointInfo[]>()
	const faceBrep = new BRep([face])

	outputLink(
		assert,
		{
			a: faceBrep,
			b: brep,
			c: new BRep(resultFaces).translate(20, 0, 0),
		},
		desc,
	)
	brep.faces.forEach(face2 => {
		face.intersectFace(face2, faceBrep, brep, faceMap, edgeMap, new CustomMap(), new CustomSet())
	})
	console.log('faceMap', faceMap)
	const edgeLooseSegments = BRep.getLooseEdgeSegments(edgeMap, faceBrep.edgeFaces)
	const newFaces: Face[] = []
	faceBrep.reconstituteFaces([face], edgeLooseSegments, faceMap, newFaces)
	assert.equal(newFaces.length, resultFaces.length, 'number of new faces')
	outputLink(
		assert,
		{
			a: faceBrep,
			b: brep,
			c: new BRep(newFaces).translate(20, 0, 0),
		},
		desc,
	)
	resultFaces.forEach(face => {
		assert.ok(
			newFaces.some(newFace => newFace.likeFace(face)),
			`newFaces.some(newFace => newFace.likeFace(face) ${newFaces.toSource()}`,
		)
	})
}

function doTest3(assert, face: Face, newEdges: Edge[], points: Map<Edge, V3[]>, resultFaces: Face[], desc: string) {
	const faceBrep = new BRep([face])
	assert.ok(
		true,
		`<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
href='brep2.html?a=${faceBrep.toSource()}&c=${new BRep(resultFaces).toSource()}.translate(20, 0, 0)'>${desc}</a>`,
	)
	const isps = Array.from(points.entries())
		.map(([edge, ps]) =>
			ps.map(p => ({ edge: edge, p: p, edgeT: snap(snap(edge.curve.pointT(p), edge.aT), edge.bT) })),
		)
		.concatenated()
	const edgeLooseSegments = BRep.prototype.getLooseEdgeSegments(new Map().set(face, isps))
	const newFaces = []
	BRep.reconstituteFaces([face], edgeLooseSegments, new Map().set(face, newEdges), newFaces)
	assert.ok(
		true,
		`<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
        href='brep2.html?a=${faceBrep.toSource()}&c=${new BRep(newFaces).toSource()}.translate(20, 0, 0)'>result</a>`,
	)
	assert.equal(newFaces.length, resultFaces.length, 'number of new faces')
	resultFaces.forEach(face => {
		assert.ok(
			newFaces.some(newFace => newFace.likeFace(face)),
			`newFaces.some(newFace => newFace.likeFace(face) ${newFaces.toSource()}`,
		)
	})
}

// prettier-ignore
suite('BRep 2', () => {
	test('tetrahedron intersections doTestWithBrep', assert => {
		const box = B2T.box(8, 9, 10, 'box')
		const face = box.faces.find(face => (face as PlaneFace).surface.plane.normal1.like(V3.Z))
		let testBrep
		// splitting contour in base position:

		testBrep = B2T.tetrahedron(V(-1, 1, 9), V(5, 1, 9), V(-1, -4, 14), V(2, -4, 10))
		doTestWithBrep(assert, face, box, testBrep, [], [],
			'face of tetra touches edge of test face, tetra intersects test volume (test face not part of result, points dont really matter)')

		testBrep = B2T.tetrahedron(V(-1, 1, 9), V(4, 1, 9), V(-1, -4, 14), V(-1, -4, 10)).flipped()
		doTestWithBrep(assert, face, box, testBrep, [StraightEdge.throughPoints(V(0, 0, 10), V(3, 0, 10))], [V(3, 0, 10), V(0, 0, 10)],
			'face of _flipped_ tetra touches edge of testface and also intersects main volume (expect point on edge)')

		testBrep = B2T.tetrahedron(V(-1, 0, 10), V(5, 0, 10), V(2, 0, 14), V(2, -4, 10))
		doTestWithBrep(assert, face, box, testBrep, [], [],
			'volumes touch edge-edge but no overlap (empty result volume; generated points dont matter)')

		testBrep = B2T.tetrahedron(V(-1, 0, 10), V(5, 0, 10), V(-1, -2, 8), V(-1, 2, 8)).flipped()
		doTestWithBrep(assert, face, box, testBrep, [StraightEdge.throughPoints(V(0, 0, 10), V(5, 0, 10))], [V(5, 0, 10), V(0, 0, 10)],
			'Tetrahedron is flipped, testface only touched at edge, needs point on edge as tetrahedron intersects side of volume')


		testBrep = B2T.tetrahedron(V(-1, 0, 10), V(5, 0, 10), V(-1, 0, 14), V(-1, -4, 9)).flipped()
		doTestWithBrep(assert, face, box, testBrep, [StraightEdge.throughPoints(V(0, 0, 10), V(5, 0, 10))], [V(0, 0, 10), V(5, 0, 10)],
			'volumes do not intersect, tetra is flipped, touches box edge-edge (result would be entire box)')
	})

	test('tetrahedorn is point on line', assert => {
		const box = B2T.box(8, 9, 10, 'box')
		let face = box.faces.find(face => (face as PlaneFace).surface.plane.normal1.like(V3.Z)), testBrep
		// splitting contour in base position:

		testBrep = B2T.tetrahedron(V(4, 0, 10), V(4, 3, 9), V(7, 3, 12), V(1, 3, 12))
		doTestWithBrep(assert, face, box, testBrep, [
				StraightEdge.throughPoints(V(3, 3, 10), V(4, 0, 10)),
				StraightEdge.throughPoints(V(4, 0, 10), V(5, 3, 10)),
				StraightEdge.throughPoints(V(5, 3, 10), V(3, 3, 10))], [],
			'face of tetra touches edge of test face, tetra intersects test volume (test face not part of result, points dont really matter)')

		doTestWithBrep(assert, face, box, testBrep.flipped(), [
				StraightEdge.throughPoints(V(4, 0, 10), V(3, 3, 10)),
				StraightEdge.throughPoints(V(5, 3, 10), V(4, 0, 10)),
				StraightEdge.throughPoints(V(3, 3, 10), V(5, 3, 10))], [V(4, 0, 10)],
			'face of tetra touches edge of test face, tetra intersects test volume (test face not part of result, points dont really matter)')

		doTestWithBrep(assert, face, box, testBrep.rotate(V(4, 0, 10), V3.Z, 90 * DEG), [
				StraightEdge.throughPoints(V(4, 0, 10), V(1, 1, 10)),
				StraightEdge.throughPoints(V(1, 1, 10), V(1, 0, 10))], [V(1, 0, 10), V(4, 0, 10)],
			'face of tetra touches edge of test face, tetra intersects test volume (test face not part of result, points dont really matter)')

		doTestWithBrep(assert, face, box, testBrep.rotate(V(4, 0, 10), V3.Z, 90 * DEG).flipped(), [
				StraightEdge.throughPoints(V(1, 1, 10), V(4, 0, 10)),
				StraightEdge.throughPoints(V(1, 0, 10), V(1, 1, 10))], [V(1, 0, 10), V(4, 0, 10)],
			'face of tetra touches edge of test face, tetra intersects test volume (test face not part of result, points dont really matter)')

		doTestWithBrep(assert, face, box, testBrep.rotate(V(4, 0, 10), V3.Z, 180 * DEG), [], [],
			'face of tetra touches edge of test face, tetra intersects test volume (test face not part of result, points dont really matter)')

		doTestWithBrep(assert, face, box, testBrep.rotate(V(4, 0, 10), V3.Z, 180 * DEG).flipped(), [], [],
			'face of tetra touches edge of test face, tetra intersects test volume (test face not part of result, points dont really matter)')
	})

	test('face box', assert => {
		const box = B2T.box(8, 9, 10, 'box')
		const face = box.faces.find(face => (face as PlaneFace).surface.plane.normal1.like(V3.Z))
		// splitting contour in base position:

		const testBrep = B2T.box(8, 1, 1).translate(0, 0, 9).flipped()
		doTestWithBrep(assert, face, box, testBrep, [
				StraightEdge.throughPoints(V(0, 1, 10), V(8, 1, 10))], [V(0, 1, 10), V(8, 1, 10)],
			'???')
	})

	test('face box 2', assert => {
		const box = B2T.box(1, 1, 6, 'box').flipped()
		const face = box.faces.find(face => (face as PlaneFace).surface.plane.normal1.like(V3.X))
		// splitting contour in base position:

		const testBrep = B2T.box(5, 5, 5)
		doTestWithBrep(assert, face, box, testBrep, [], [],
			'???')

	})

	test('face box 3', assert => {
		const box = B2T.box(10, 10, 5, 'box')
		const face = box.faces.find(face => (face.surface as PlaneSurface).plane.normal1.like(V3.Z))
		// splitting contour in base position:

		const testBrep = B2T.box(10, 10, 5, 'box').translate(3, 3, 0)
		// this fails, because the point is added as the bounds of the edge in face, where the edge belongs to a
		// neighborig face which isn't checked
		doTestWithBrep(assert, face, box, testBrep, [
				StraightEdge.throughPoints(V(3, 10, 5), V(3, 3, 5)),
				StraightEdge.throughPoints(V(3, 3, 5), V(10, 3, 5))], [V(3, 10, 5), V(10, 3, 5)],
			'???')

		doTestWithBrep(assert, face, box, testBrep, [], [],
			'???', true)

	})

	test('face box 4', assert => {
		const box = B2T.box(10, 10, 5, 'box').flipped()
		const face = box.faces.find(face => (face.surface as PlaneSurface).plane.normal1.z == -1)
		// splitting contour in base position:

		const testBrep = box.translate(3, 3, 0)
		doTestWithBrep(assert, face, box, testBrep, [
				StraightEdge.throughPoints(V(3, 10, 5), V(3, 3, 5)),
				StraightEdge.throughPoints(V(3, 3, 5), V(3, 10, 5)),
				StraightEdge.throughPoints(V(3, 3, 5), V(10, 3, 5)),
				StraightEdge.throughPoints(V(10, 3, 5), V(3, 3, 5))], [V(3, 10, 5), V(10, 3, 5)],
			'???')
		doTestWithBrep(assert, face, box, testBrep, [
				StraightEdge.throughPoints(V(3, 10, 5), V(3, 3, 5)),
				StraightEdge.throughPoints(V(3, 3, 5), V(10, 3, 5))], [V(3, 10, 5), V(10, 3, 5)],
			'???', true)

	})
	test('face box 5', assert => {
		const box = B2T.box(10, 10, 5, 'box').flipped()
		const face = box.faces.find(face => (face.surface as PlaneSurface).plane.normal1.z == -1)
		// splitting contour in base position:

		const testBrep = box.translate(0, 0, 5)
		doTestWithBrep(assert, face, box, testBrep, [], [],
			'???')
		doTestWithBrep(assert, face, box, testBrep, [], [],
			'???', true)

	})
	test('face box 6', assert => {
		const box = B2T.box(10, 10, 10, 'box').flipped()
		const box2 = B2T.box(10, 12, 12, 'box').flipped()

		const face = box.faces.find(face => (face.surface as PlaneSurface).plane.normal1.z == -1)
		// splitting contour in base position:

		const testBrep = box2.translate(0, 0, -2)
		doTestWithBrep(assert, face, box, testBrep, [
				new StraightEdge(new L3(V(0, 0, 10), V(0, -1, 0)), V(0, 0, 10), V(0, 10, 10), 0, -10),
				new StraightEdge(new L3(V(10, 0, 10), V(0, 1, 0)), V(10, 10, 10), V(10, 0, 10), 10, 0),
				new StraightEdge(new L3(V(0, 0, 10), V(1, 0, 0)), V(10, 0, 10), V(0, 0, 10), 10, 0)], [
				V(0, 10, 10), V(0, 0, 10), V(10, 0, 10), V(10, 10, 10),
			],
			'???')
		doTestWithBrep(assert, face, box, testBrep, [], [],
			'???', true)

	})
	test('face box 7', assert => {
		const box = B2T.box(1, 1, 1, 'box')
		const face = box.faces.find(face => (face as PlaneFace).surface.plane.normal1.like(V3.Z))
		// splitting contour in base position:

		const testBrep = B2T.box(1 / 3, 1 / 3, 1).flipped()
		doTestWithBrep(assert, face, box, testBrep, [
				new StraightEdge(new L3(V(0, 1 / 3, 1), V(1, 0, 0)), V(0, 1 / 3, 1), V(1 / 3, 1 / 3, 1), 0, 1 / 3),
				new StraightEdge(new L3(V(1 / 3, 0, 1), V(0, -1, 0)), V(1 / 3, 1 / 3, 1), V(1 / 3, 0, 1), -1 / 3, 0)],
			[V(1 / 3, 0, 1), V(0, 1 / 3, 1)],
			'???')

	})
	test('sphere face box', assert => {
		const sphere = B2T.sphere(4, 'ball')
		const face = sphere.faces.find(face => face.containsPoint(new V3(0, 4, 0)))
		// splitting contour in base position:

		const testBrep = B2T.box(10, 10, 10, 'box').translate(1, 2, 3).flipped()
		doTestWithBrep(assert, face, sphere, testBrep, [
			new PCurveEdge(
				new SemiEllipseCurve(V(1, 0, 0), V(0, 0, 3.872983346207417), V(0, 3.872983346207417, 0), 0, 3.141592653589793),
				V(1, 2, 3.3166247903554), V(1, 2.449489742783178, 3),
				0.5426391022496527, 0.684719203002283, undefined,
				V(0, 3.3166247903554003, -2), V(0, 3, -2.4494897427831783)),
			new PCurveEdge(
				new SemiEllipseCurve(V(0, 2, 0), V(0, 0, -3.4641016151377544), V(3.4641016151377544, 0, 0), 0, 3.141592653589793),
				V(1.7320508075688767, 2, 3), V(1, 2, 3.3166247903554),
				2.6179938779914944, 2.8487498818612176, undefined,
				V(-3, 0, 1.732050807568877), V(-3.3166247903554, 0, 1)),
			new PCurveEdge(
				new SemiEllipseCurve(V(0, 0, 3), V(-2.6457513110645907, 0, 0), V(0, 2.6457513110645907, 0), 0, 3.141592653589793),
				V(1, 2.449489742783178, 3), V(1.7320508075688767, 2, 3),
				1.9583930134500773, 2.284520705739662, undefined,
				V(2.449489742783178, -1, 0), V(2, -1.7320508075688765, 0)),
		], [], '???')

	})
	test('sphere face box 2', assert => {
		const sphere = B2T.sphere(5, 'ball')
		const face = sphere.faces.find(face => face.containsPoint(new V3(0, 5, 0)))
		// splitting contour in base position:

		const testBrep = B2T.box(4, 4, 4, 'box').translate(0, 1, 0).flipped()
		doTestWithBrep(assert, face, sphere, testBrep, [], [], '???')

	})

	test('coplanar things', assert => {
		const box = B2T.box(8, 9, 10, 'box')
		const face = box.faces.find(face => (face as PlaneFace).surface.plane.normal1.like(V3.Z))
		// splitting contour in base position:

		doTestWithBrep(assert, face, box,
			B2T.tetrahedron(V(-1, -1, 10), V(-1, 5, 10), V(5, -1, 10), V(-1, -1, 5)).flipped(),
			[StraightEdge.throughPoints(V(0, 4, 10), V(4, 0, 10))], [V(0, 4, 10), V(4, 0, 10)],
			'cut off corner of box with flipped tetra (anti)coplanar to test face')

		doTestWithBrep(assert, face, box,
			B2T.tetrahedron(V(-1, -1, 10), V(-1, 5, 10), V(4, 0, 10), V(-1, -1, 5)).flipped(),
			[StraightEdge.throughPoints(V(0, 4, 10), V(4, 0, 10))], [V(0, 4, 10), V(4, 0, 10)],
			'cut off corner of box with flipped tetra (anti)coplanar to test face')

		doTestWithBrep(assert, face, box,
			B2T.tetrahedron(V(1, 0, 10), V(1, 5, 10), V(4, 0, 10), V(-1, -1, 5)).flipped(),
			[StraightEdge.throughPoints(V(1, 0, 10), V(1, 5, 10)), StraightEdge.throughPoints(V(1, 5, 10), V(4, 0, 10))],
			[V(1, 0, 10), V(4, 0, 10)],
			'cut off corner of box with flipped tetra (anti)coplanar to test face')

		doTestWithBrep(assert, face, box,
			B2T.tetrahedron(V(0, 0, 10), V(0, 5, 10), V(5, 0, 10), V(0, 0, 5)).flipped(),
			[StraightEdge.throughPoints(V(0, 5, 10), V(5, 0, 10))],
			[V(0, 5, 10), V(5, 0, 10), V(0, 0, 10)],
			'cut off corner of box with flipped tetra (anti)coplanar to test face')

		doTestWithBrep(assert, face, box,
			B2T.tetrahedron(V(0, 0, 10), V(-1, 5, 10), V(5, -1, 10), V(-1, -1, 5)).flipped(),
			[StraightEdge.throughPoints(V(0, 4, 10), V(4, 0, 10))], [V(0, 4, 10), V(4, 0, 10)],
			'cut off corner of box with flipped tetra (anti)coplanar to test face')

		doTestWithBrep(assert, face, box,
			B2T.tetrahedron(V(0, 0, 10), V(1, 5, 10), V(5, 1, 10), V(0, 0, 5)).flipped(),
			[StraightEdge.throughPoints(V(0, 0, 10), V(1, 5, 10)), StraightEdge.throughPoints(V(1, 5, 10), V(5, 1, 10)), StraightEdge.throughPoints(V(5, 1, 10), V(0, 0, 10))],
			[],
			'cut hole at corner of test face')

	})
	test('coplanar things 2', assert => {
		const box = B2T.box(8, 9, 10, 'box')
		const face = box.faces.find(face => (face as PlaneFace).surface.plane.normal1.like(V3.Z))
		// splitting contour in base position:

		const testBrep = B2T.tetrahedron(V(-1, -1, 10), V(-1, 5, 10), V(5, -1, 10), V(-1, -1, 5))
		doTestWithBrep(assert, face, box, testBrep, [StraightEdge.throughPoints(V(0, 4, 10), V(4, 0, 10))], [V(0, 4, 10), V(4, 0, 10)],
			'cut off corner of box with flipped tetra (anti)coplanar to test face')

		const testBrep2 = B2T.tetrahedron(V(-1, -1, 10), V(-1, 5, 10), V(5, -1, 10), V(-1, -1, 5)).flipped()
		doTestWithBrep(assert, face, box.flipped(), testBrep2, [StraightEdge.throughPoints(V(0, 4, 10), V(4, 0, 10))], [V(0, 4, 10), V(4, 0, 10)],
			'cut off corner of box with flipped tetra (anti)coplanar to test face')

	})


	skip('assembly', assert => {
		const baseFace = PlaneFace.forVertices(P3.XY, [V(0, 0), V(10, 0), V(10, 10), V(0, 10)])

		const extrude = B2T.extrudeVertices([V(5, -1), V(2, 2), V(8, 2)], P3.XY.flipped(), V(0, 0, 10)).translate(0, 0, -2).flipped()
		const result = PlaneFace.forVertices(P3.XY, [V(0, 0), V(4, 0), V(2, 2), V(8, 2), V(6, 0), V(10, 0), V(10, 10), V(0, 10)])
		doTest2(assert, baseFace, extrude, [result], 'volume cuts edge of test face (twice)')

		const edges = StraightEdge.chain([V(5, 0), V(2, 3), V(8, 3)])
		const result2 = PlaneFace.forVertices(P3.XY, [V(0, 0), V(10, 0), V(10, 10), V(0, 10)], [V(5, 0), V(2, 3), V(8, 3)])
		doTest3(assert, baseFace, edges, new Map(), [result2], 'volume touches inside of test face edge')


		// from test case 3:4
		// V-shape spine touching middle, splits volume enclosing in both directions
		const brep = new BRep([
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).flipped(),
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).rotateX(PI / 4),
		]).rotateX(-PI / 2).flipped().translate(0, 2, 0)
		// degenerate cycle in middle of face
		const result4 = PlaneFace.forVertices(P3.XY, [V(0, 0), V(10, 0), V(10, 10), V(0, 10)], [V(5, 2), V(0, 2)])


		const extrude5 = B2T.extrudeVertices([V(0, 0), V(3, 2), V(2, 3)], P3.XY.flipped(), V(0, 0, 10)).translate(0, 0, -2).flipped()
		const result5 = PlaneFace.forVertices(P3.XY, [V(0, 0), V(2, 3), V(3, 2), V(0, 0), V(10, 0), V(10, 10), V(0, 10)])
		doTest2(assert, baseFace, extrude5, [result5], 'volume touches inside of test face corner')

		const edges6 = StraightEdge.chain([V(1, 0), V(8, 10), V(3, 0)], false)
		const points6 = new Map().set(baseFace.contour.find(edge => edge.aDir.like(V3.X)), [V(1, 0), V(3, 0)])
		const result6 = [PlaneFace.forVertices(P3.XY, [V(0, 0), V(1, 0), V(8, 10), V(3, 0), V(10, 0), V(10, 10), V(0, 10)])]
		doTest3(assert, baseFace, edges6, points6, result6, 'volume touches inside of test face corner')

		const edges7 = [StraightEdge.throughPoints(V(1, 0), V(7, 10)), StraightEdge.throughPoints(V(8, 10), V(3, 0))]
		const points7 = new Map()
			.set(baseFace.contour.find(edge => edge.aDir.like(V3.X)), [V(1, 0), V(3, 0)])
			.set(baseFace.contour.find(edge => edge.aDir.like(V3.X.negated())), [V(7, 10), V(8, 10)])
		const result7 = [PlaneFace.forVertices(P3.XY, [V(0, 0), V(1, 0), V(7, 10), V(0, 10)]), PlaneFace.forVertices(P3.XY, [V(10, 10), V(8, 10), V(3, 0), V(10, 0)])]
		doTest3(assert, baseFace, edges7, points7, result7, 'volume touches inside of test face corner')


		const edges8 = [StraightEdge.throughPoints(V(1, 0), V(8, 10)), StraightEdge.throughPoints(V(8, 10), V(1, 0))]
		const points8 = new Map()
			.set(baseFace.contour.find(edge => edge.aDir.like(V3.X)), [V(1, 0)])
			.set(baseFace.contour.find(edge => edge.aDir.like(V3.X.negated())), [V(8, 10)])
		const result8 = [PlaneFace.forVertices(P3.XY, [V(0, 0), V(1, 0), V(8, 10), V(0, 10)]), PlaneFace.forVertices(P3.XY, [V(10, 10), V(8, 10), V(1, 0), V(10, 0)])]
		doTest3(assert, baseFace, edges8, points8, result8, 'volume touches inside of test face corner')


		const edges9 = [StraightEdge.throughPoints(V(1, 0), V(8, 5)), StraightEdge.throughPoints(V(8, 5), V(1, 0))]
		const points9 = new Map().set(baseFace.contour.find(edge => edge.aDir.like(V3.X)), [V(1, 0)])
		const result9 = [PlaneFace.forVertices(P3.XY, [V(0, 0), V(1, 0), V(8, 5), V(1, 0), V(10, 0), V(10, 10), V(0, 10)])]
		doTest3(assert, baseFace, edges9, points9, result9, 'volume touches inside of test face corner')


	})

	test('assembly holes', assert => {
		const baseFace = PlaneFace.forVertices(P3.XY, [V(0, 0), V(10, 0), V(10, 10), V(0, 10)])

		const extrude7 = B2T.extrudeVertices([V(0, 0), V(5, 8), V(8, 5)], P3.XY.flipped(), V(0, 0, 10)).translate(0, 0, -2).flipped()
		const result7 = [PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 8), V(8, 5), V(0, 0), V(10, 0), V(10, 10), V(0, 10)])]
		doTest2(assert, baseFace, extrude7, result7, 'volume touches inside of test face corner')

		const extrude8 = B2T.extrudeVertices([V(1, 1), V(1, -1), V(-1, 1)], P3.XY.flipped(), V(0, 0, 10)).translate(0, 0, -2).flipped()
		const result8 = [PlaneFace.forVertices(P3.XY, [V(1, 0), V(10, 0), V(10, 10), V(0, 10), V(0, 1), V(1, 1)])]
		doTest2(assert, baseFace, extrude8, result8, 'volume touches inside of test face corner')

		const extrude9 = B2T.extrudeVertices([V(-1, -1), V(1, 1), V(1, -1)], P3.XY.flipped(), V(0, 0, 10)).translate(0, 0, -2).flipped()
		const result9 = [PlaneFace.forVertices(P3.XY, [V(1, 0), V(10, 0), V(10, 10), V(0, 10), V(0, 0), V(1, 1)])]
		doTest2(assert, baseFace, extrude9, result9, 'volume touches inside of test face corner')

		const extrude10 = B2T.extrudeVertices([V(1, 1), V(2, 2), V(2, 1)], P3.XY.flipped(), V(0, 0, 10)).translate(0, 0, -2).flipped()
		const result10 = [PlaneFace.forVertices(P3.XY, [V(10, 0), V(10, 10), V(0, 10), V(0, 0)], [V(1, 1), V(2, 2), V(2, 1)])]
		doTest2(assert, baseFace, extrude10, result10, 'adding a hole')

		const baseFace11 = PlaneFace.forVertices(P3.XY, [V(10, 0), V(10, 10), V(0, 10), V(0, 0)], [V(1, 1), V(2, 2), V(2, 1)])
		const extrude11 = B2T.extrudeVertices([V(5, 5), V(6, 6), V(6, 5)], P3.XY.flipped(), V(0, 0, 10)).translate(0, 0, -2).flipped()
		const result11 = [PlaneFace.forVertices(P3.XY, [V(10, 0), V(10, 10), V(0, 10), V(0, 0)], [V(1, 1), V(2, 2), V(2, 1)], [V(5, 5), V(6, 6), V(6, 5)])]
		doTest2(assert, baseFace11, extrude11, result11, 'adding a hole to a face which already has one')


		const baseFace12 = PlaneFace.forVertices(P3.XY, [V(10, 0), V(10, 10), V(0, 10), V(0, 0)], [V(1, 1), V(5, 5), V(5, 1)])
		const extrude12 = B2T.extrudeVertices([V(2, 1.5), V(2, 4), V(4.5, 4)], P3.XY.flipped(), V(0, 0, 10)).translate(0, 0, -2).flipped()
		const result12 = [PlaneFace.forVertices(P3.XY, [V(10, 0), V(10, 10), V(0, 10), V(0, 0)], [V(5, 5), V(5, 1), V(1, 1), V(2, 2), V(2, 4), V(4, 4)])]
		doTest2(assert, baseFace12, extrude12, result12, 'extending an existing hole')

		const baseFace13 = PlaneFace.forVertices(P3.XY, [V(10, 0), V(10, 10), V(0, 10), V(0, 0)], [V(1, 1), V(5, 5), V(5, 1)])
		const extrude13 = B2T.extrudeVertices([V(3, -1), V(4, -1), V(4, 2), V(3, 2)], P3.XY.flipped(), V(0, 0, 10)).translate(0, 0, -2).flipped()
		const result13 = [PlaneFace.forVertices(P3.XY, [V(10, 0), V(10, 10), V(0, 10), V(0, 0), V(3, 0), V(3, 1), V(1, 1), V(5, 5), V(5, 1), V(4, 1), V(4, 0)])]
		doTest2(assert, baseFace13, extrude13, result13, 'removing a hole by cutting a channel')

		const baseFace14 = PlaneFace.forVertices(P3.XY, [V(10, 0), V(10, 10), V(0, 10), V(0, 0)], [V(1, 1), V(5, 5), V(5, 1)])
		const edges14 = StraightEdge.chain([V(1, 1), V(1, 5), V(5, 5)], false)
		const points14 = new Map()
		const result14 = [PlaneFace.forVertices(P3.XY, [V(10, 0), V(10, 10), V(0, 10), V(0, 0)], [V(5, 5), V(5, 1), V(1, 1), V(1, 5)])]
		doTest3(assert, baseFace14, edges14, points14, result14, 'extending an existing hole')
	})

	test('remove half of a half-pie', assert => {
		const pie = B2T.puckman(8, 180 * DEG, 5, 'pie/2')
		const boxKnife = B2T.box(11, 10, 7, 'knife').translate(-10, -1, -1).flipped()

		const resultTopPoint = V(1, 8 * Math.sin(Math.acos(1 / 8)), 0)
		const result = B2T.extrudeEdges([
				StraightEdge.throughPoints(V(8, 0, 0), V(1, 0, 0)),
				StraightEdge.throughPoints(V(1, 0, 0), resultTopPoint),
				PCurveEdge.forCurveAndTs(SemiEllipseCurve.semicircle(8), Math.acos(1 / 8), 0)],
			P3.XY.flipped(), V(0, 0, 5), 'pie/4')
		testBRepAnd(assert, pie, boxKnife, result)


		const testFace = pie.faces.find(face => face.surface.plane && face.surface.plane.normal1.like(V3.Y.negated()) && face.contour.some(edge => edge.a.x > 1))
		const k2 = boxKnife.translate(-1, 0, 0)
		doTestWithBrep(assert, testFace, pie, k2.flipped(), [], [V3.O, V(0, 0, 5)],
			'volumes touch edge-edge but no overlap (empty result volume; generated points dont matter)')
		testBRepOp(assert, pie, k2, () => pie.minus(k2), B2T.puckman(8, 90 * DEG, 5, 'pie/4'))
	})
	test('cut hole through side of pie', assert => {
		const pie = B2T.puckman(8, 180 * DEG, 5, 'pie/2')
		const punch = B2T.box(5, 10, 3, 'knife').translate(1, -1, 1).flipped()
		const result = new BRep([
			new RotationFace(new SemiCylinderSurface(new SemiEllipseCurve(V(0, 0, 0), V(8, 0, 0), V(0, 8, 0)), V3.Z, undefined, undefined), [
				new StraightEdge(new L3(V(8, 0, 0), V(0, 0, 1)), V(8, 0, 5), V(8, 0, 0), 5, 0),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(8, 0, 0), V(0, 8, 0)), V(8, 0, 0), V(-8, 0, 0), 0, 3.141592653589793, undefined, V(0, 8, 0), V(-0, -8, 0)),
				new StraightEdge(new L3(V(-8, 0, 0), V(0, 0, 1)), V(-8, 0, 0), V(-8, 0, 5), 0, 5),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 5), V(8, 0, 0), V(0, 8, 0)), V(-8, 0, 5), V(8, 0, 5), 3.141592653589793, 0, undefined, V(0, 8, 0), V(0, -8, 0))], [[
				new StraightEdge(new L3(V(1, 7.937253933193772, 0), V(0, 0, -1)), V(1, 7.937253933193773, 4), V(1, 7.937253933193773, 1), -4, -1),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 1), V(-8, 0, 0), V(0, 8, 0)), V(1, 7.937253933193773, 1), V(6, 5.291502622129181, 1), 1.696124157962962, 2.4188584057763776, undefined, V(7.937253933193772, -1, 0), V(5.291502622129181, -6, 0)),
				new StraightEdge(new L3(V(6, 5.2915026221291805, 0), V(0, 0, 1)), V(6, 5.291502622129181, 1), V(6, 5.291502622129181, 4), 1, 4),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 4), V(8, 0, 0), V(0, 8, 0)), V(6, 5.291502622129181, 4), V(1, 7.937253933193773, 4), 0.7227342478134156, 1.4454684956268313, undefined, V(-5.2915026221291805, 6, 0), V(-7.937253933193772, 1, 0))]]),
			new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0)), [
				new StraightEdge(new L3(V(0, 0, 0), V(1, 0, 0)), V(0, 0, 0), V(8, 0, 0), 0, 8),
				new StraightEdge(new L3(V(8, 0, 0), V(0, 0, 1)), V(8, 0, 0), V(8, 0, 5), 0, 5),
				new StraightEdge(new L3(V(8, 0, 5), V(-1, 0, 0)), V(8, 0, 5), V(0, 0, 5), 0, 8),
				new StraightEdge(new L3(V(0, 0, 5), V(0, 0, -1)), V(0, 0, 5), V(0, 0, 0), 0, 5)], [[
				new StraightEdge(new L3(V(1, 0, 0), V(0, 0, 1)), V(1, 0, 1), V(1, 0, 4), 1, 4),
				new StraightEdge(new L3(V(0, 0, 4), V(1, 0, 0)), V(1, 0, 4), V(6, 0, 4), 1, 6),
				new StraightEdge(new L3(V(6, 0, 0), V(0, 0, -1)), V(6, 0, 4), V(6, 0, 1), -4, -1),
				new StraightEdge(new L3(V(0, 0, 1), V(-1, 0, 0)), V(6, 0, 1), V(1, 0, 1), -6, -1)]]),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 0)), [
				new StraightEdge(new L3(V(0, 0, 0), V(1, 0, 0)), V(8, 0, 0), V(0, 0, 0), 8, 0),
				new StraightEdge(new L3(V(0, 0, 0), V(-1, 0, 0)), V(0, 0, 0), V(-8, 0, 0), 0, 8),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 0), V(8, 0, 0), V(0, 8, 0)), V(-8, 0, 0), V(8, 0, 0), 3.141592653589793, 0, undefined, V(0, 8, 0), V(0, -8, 0))], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 5)), [
				new StraightEdge(new L3(V(8, 0, 5), V(-1, 0, 0)), V(0, 0, 5), V(8, 0, 5), 8, 0),
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 5), V(8, 0, 0), V(0, 8, 0)), V(8, 0, 5), V(-8, 0, 5), 0, 3.141592653589793, undefined, V(0, 8, 0), V(-0, -8, 0)),
				new StraightEdge(new L3(V(-8, 0, 5), V(1, -0, 0)), V(-8, 0, 5), V(0, 0, 5), 0, 8)], []),
			new PlaneFace(new PlaneSurface(new P3(V(-0, -1, 0), 0)), [
				new StraightEdge(new L3(V(0, 0, 5), V(0, 0, -1)), V(0, 0, 0), V(0, 0, 5), 5, 0),
				new StraightEdge(new L3(V(-8, 0, 5), V(1, -0, 0)), V(0, 0, 5), V(-8, 0, 5), 8, 0),
				new StraightEdge(new L3(V(-8, 0, 0), V(0, 0, 1)), V(-8, 0, 5), V(-8, 0, 0), 5, 0),
				new StraightEdge(new L3(V(0, 0, 0), V(-1, 0, 0)), V(-8, 0, 0), V(0, 0, 0), 8, 0)], []),
			new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 1)), [
				new StraightEdge(new L3(V(1, 7.937253933193772, 0), V(0, 0, -1)), V(1, 7.937253933193773, 1), V(1, 7.937253933193773, 4), -1, -4),
				new StraightEdge(new L3(V(1, -1, 4), V(0, 1, 0)), V(1, 7.937253933193773, 4), V(1, 0, 4), 8.937253933193773, 1),
				new StraightEdge(new L3(V(1, 0, 0), V(0, 0, 1)), V(1, 0, 4), V(1, 0, 1), 4, 1),
				new StraightEdge(new L3(V(1, -1, 1), V(0, 1, 0)), V(1, 0, 1), V(1, 7.937253933193773, 1), 1, 8.937253933193773)], []),
			new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), -6)), [
				new StraightEdge(new L3(V(6, 5.2915026221291805, 0), V(0, 0, 1)), V(6, 5.291502622129181, 4), V(6, 5.291502622129181, 1), 4, 1),
				new StraightEdge(new L3(V(6, 9, 1), V(0, -1, 0)), V(6, 5.291502622129181, 1), V(6, 0, 1), 3.7084973778708186, 9),
				new StraightEdge(new L3(V(6, 0, 0), V(0, 0, -1)), V(6, 0, 1), V(6, 0, 4), -1, -4),
				new StraightEdge(new L3(V(6, 9, 4), V(0, -1, 0)), V(6, 0, 4), V(6, 5.291502622129181, 4), 9, 3.7084973778708186)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 1)), [
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 1), V(-8, 0, 0), V(0, 8, 0)), V(6, 5.291502622129181, 1), V(1, 7.937253933193773, 1), 2.4188584057763776, 1.696124157962962, undefined, V(-5.291502622129181, 6, 0), V(-7.937253933193772, 1, 0)),
				new StraightEdge(new L3(V(1, -1, 1), V(0, 1, 0)), V(1, 7.937253933193773, 1), V(1, 0, 1), 8.937253933193773, 1),
				new StraightEdge(new L3(V(0, 0, 1), V(-1, 0, 0)), V(1, 0, 1), V(6, 0, 1), -1, -6),
				new StraightEdge(new L3(V(6, 9, 1), V(0, -1, 0)), V(6, 0, 1), V(6, 5.291502622129181, 1), 9, 3.7084973778708186)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), -4)), [
				new PCurveEdge(new SemiEllipseCurve(V(0, 0, 4), V(8, 0, 0), V(0, 8, 0)), V(1, 7.937253933193773, 4), V(6, 5.291502622129181, 4), 1.4454684956268313, 0.7227342478134156, undefined, V(7.937253933193772, -1, 0), V(5.2915026221291805, -6, 0)),
				new StraightEdge(new L3(V(6, 9, 4), V(0, -1, 0)), V(6, 5.291502622129181, 4), V(6, 0, 4), 3.7084973778708186, 9),
				new StraightEdge(new L3(V(0, 0, 4), V(1, 0, 0)), V(6, 0, 4), V(1, 0, 4), 6, 1),
				new StraightEdge(new L3(V(1, -1, 4), V(0, 1, 0)), V(1, 0, 4), V(1, 7.937253933193773, 4), 1, 8.937253933193773)], [])], false)
		testBRepAnd(assert, pie, punch, result)
	})

	test('minus including ProjectedCurveSurface', assert =>  {
		const s1 = new ProjectedCurveSurface(new BezierCurve(V(0, 0, 0), V(-5, 5, 0), V(15, 5, 0), V(10, 0, 0), -0.1, 1.1), V(0, 0, -1), 0, 1)
		const s2 = new ProjectedCurveSurface(new BezierCurve(V(0, 0, 0), V(-5, 5, 0), V(15, 5, 0), V(10, 0, 0), -0.1, 1.1), V(0, 0, -1), 0, 1)
		const a = B2T.extrudeEdges([
				StraightEdge.throughPoints(V(10, 0, 0), V(0, 0, 0)),
				PCurveEdge.forCurveAndTs(new BezierCurve(V3.O, V(-5, 5, 0), V(15, 5, 0), V(10, 0, 0)), 0, 1)],
			P3.XY.flipped(), V(0, 0, 5), 'a/4')
		const punch = B2T.box(5, 10, 3, 'knife').translate(1, -1, 1).flipped()

		const result = new BRep([
			new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0), V(0, 0, -1), V(1, 0, 0)), [
				new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(0, 0, 0), V(10, 0, 0), 10, 0),
				new StraightEdge(new L3(V(10, 0, 0), V(0, 0, 1)), V(10, 0, 0), V(10, 0, 5), 0, 5),
				new StraightEdge(new L3(V(10, 0, 5), V(-1, 0, 0)), V(10, 0, 5), V(0, 0, 5), 0, 10),
				new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 5), V(0, 0, 0), 5, 0)], [[
				new StraightEdge(new L3(V(1, 0, 0), V(0, 0, 1)), V(1, 0, 1), V(1, 0, 4), 1, 4),
				new StraightEdge(new L3(V(0, 0, 4), V(1, 0, 0)), V(1, 0, 4), V(6, 0, 4), 1, 6),
				new StraightEdge(new L3(V(6, 0, 0), V(0, 0, -1)), V(6, 0, 4), V(6, 0, 1), -4, -1),
				new StraightEdge(new L3(V(0, 0, 1), V(-1, 0, 0)), V(6, 0, 1), V(1, 0, 1), -6, -1)]]),
			new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0, 0, 0), V(-5, 5, 0), V(15, 5, 0), V(10, 0, 0), -0.1, 1.1), V(0, 0, -1), 0, 1, -Infinity, Infinity), [
				new PCurveEdge(new BezierCurve(V(0, 0, 0), V(-5, 5, 0), V(15, 5, 0), V(10, 0, 0), -0.1, 1.1), V(10, 0, 0), V(0, 0, 0), 1, 0, undefined, V(15, 15, 0), V(15, -15, 0)),
				new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 0), V(0, 0, 5), 0, 5),
				new PCurveEdge(new BezierCurve(V(0, 0, 5), V(-5, 5, 5), V(15, 5, 5), V(10, 0, 5), -0.1, 1.1), V(0, 0, 5), V(10, 0, 5), 0, 1, undefined, V(-15, 15, 0), V(-15, -15, 0)),
				new StraightEdge(new L3(V(10, 0, 0), V(0, 0, 1)), V(10, 0, 5), V(10, 0, 0), 5, 0)], [[
				new StraightEdge(new L3(V(0.9999999999999971, 3.1854361043800057, 0), V(0, 0, -1)), V(1, 3.185436104380006, 4), V(1, 3.185436104380006, 1), -4, -1),
				new PCurveEdge(new BezierCurve(V(0, 0, 1), V(-5, 5, 1), V(15, 5, 1), V(10, 0, 1), -0.1, 1.1), V(1, 3.185436104380006, 1), V(6, 3.720106174228432, 1), 0.3059958942668148, 0.5446421518086274, undefined, V(16.854361043800058, 5.820123171995554, 0), V(22.201061742284324, -1.339264554258822, 0)),
				new StraightEdge(new L3(V(5.999999999999993, 3.7201061742284325, 0), V(0, 0, 1)), V(6, 3.720106174228432, 1), V(6, 3.720106174228432, 4), 1, 4),
				new PCurveEdge(new BezierCurve(V(10, 0, 4), V(15, 5, 4), V(-5, 5, 4), V(0, 0, 4), -0.10000000000000009, 1.1), V(6, 3.720106174228432, 4), V(1, 3.185436104380006, 4), 0.4553578481913726, 0.6940041057331853, undefined, V(-22.20106174228432, 1.339264554258822, 0), V(-16.85436104380006, -5.820123171995558, 0))]]),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 0), V(-1, 0, 0), V(0, 1, 0)), [
				new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(10, 0, 0), V(0, 0, 0), 0, 10),
				new PCurveEdge(new BezierCurve(V(0, 0, 0), V(-5, 5, 0), V(15, 5, 0), V(10, 0, 0), -0.1, 1.1), V(0, 0, 0), V(10, 0, 0), 0, 1, undefined, V(-15, 15, 0), V(-15, -15, 0))], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 5), V(1, 0, 0), V(0, 1, 0)), [
				new PCurveEdge(new BezierCurve(V(0, 0, 5), V(-5, 5, 5), V(15, 5, 5), V(10, 0, 5), -0.1, 1.1), V(10, 0, 5), V(0, 0, 5), 1, 0, undefined, V(15, 15, 0), V(15, -15, 0)),
				new StraightEdge(new L3(V(10, 0, 5), V(-1, 0, 0)), V(0, 0, 5), V(10, 0, 5), 10, 0)], []),
			new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 1), V(0, 0, -1), V(0, 1, 0)), [
				new StraightEdge(new L3(V(1, 0, 0), V(0, 0, 1)), V(1, 0, 4), V(1, 0, 1), 4, 1),
				new StraightEdge(new L3(V(1, -1, 1), V(0, 1, 0)), V(1, 0, 1), V(1, 3.185436104380006, 1), 1, 4.185436104380006),
				new StraightEdge(new L3(V(0.9999999999999971, 3.1854361043800057, 0), V(0, 0, -1)), V(1, 3.185436104380006, 1), V(1, 3.185436104380006, 4), -1, -4),
				new StraightEdge(new L3(V(1, -1, 4), V(0, 1, 0)), V(1, 3.185436104380006, 4), V(1, 0, 4), 4.185436104380006, 1)], []),
			new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), -6), V(0, 0, -1), V(0, -1, 0)), [
				new StraightEdge(new L3(V(6, 0, 0), V(0, 0, -1)), V(6, 0, 1), V(6, 0, 4), -1, -4),
				new StraightEdge(new L3(V(6, 9, 4), V(0, -1, 0)), V(6, 0, 4), V(6, 3.720106174228432, 4), 9, 5.279893825771568),
				new StraightEdge(new L3(V(5.999999999999993, 3.7201061742284325, 0), V(0, 0, 1)), V(6, 3.720106174228432, 4), V(6, 3.720106174228432, 1), 4, 1),
				new StraightEdge(new L3(V(6, 9, 1), V(0, -1, 0)), V(6, 3.720106174228432, 1), V(6, 0, 1), 5.279893825771568, 9)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 1), V(-1, 0, 0), V(0, -1, 0)), [
				new StraightEdge(new L3(V(0, 0, 1), V(-1, 0, 0)), V(1, 0, 1), V(6, 0, 1), -1, -6),
				new StraightEdge(new L3(V(6, 9, 1), V(0, -1, 0)), V(6, 0, 1), V(6, 3.720106174228432, 1), 9, 5.279893825771568),
				new PCurveEdge(new BezierCurve(V(0, 0, 1), V(-5, 5, 1), V(15, 5, 1), V(10, 0, 1), -0.1, 1.1), V(6, 3.720106174228432, 1), V(1, 3.185436104380006, 1), 0.5446421518086274, 0.3059958942668148, undefined, V(-22.201061742284324, 1.339264554258822, 0), V(-16.854361043800058, -5.820123171995554, 0)),
				new StraightEdge(new L3(V(1, -1, 1), V(0, 1, 0)), V(1, 3.185436104380006, 1), V(1, 0, 1), 4.185436104380006, 1)], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), -4), V(1, 0, 0), V(0, -1, 0)), [
				new StraightEdge(new L3(V(0, 0, 4), V(1, 0, 0)), V(6, 0, 4), V(1, 0, 4), 6, 1),
				new StraightEdge(new L3(V(1, -1, 4), V(0, 1, 0)), V(1, 0, 4), V(1, 3.185436104380006, 4), 1, 4.185436104380006),
				new PCurveEdge(new BezierCurve(V(10, 0, 4), V(15, 5, 4), V(-5, 5, 4), V(0, 0, 4), -0.10000000000000009, 1.1), V(1, 3.185436104380006, 4), V(6, 3.720106174228432, 4), 0.6940041057331853, 0.4553578481913726, undefined, V(16.85436104380006, 5.820123171995558, 0), V(22.20106174228432, -1.339264554258822, 0)),
				new StraightEdge(new L3(V(6, 9, 4), V(0, -1, 0)), V(6, 3.720106174228432, 4), V(6, 0, 4), 5.279893825771568, 9)], [])], false)
		testBRepAnd(assert, a, punch, result)
	})

	suite('intersection of volumes where no faces intersect', () => {
		/**
		 * A in B:
		 * A - B = 0
		 * B - A = B + A flipped
		 * A + B = B
		 * intersection = A
		 *
		 * A and B outside each other
		 * A - B = A (B - A is equivalent)
		 * A + B = A + B
		 * intersection = 0
		 */

		const b = B2T.box(10, 10, 10), a = B2T.box(2, 2, 2).translate(1, 1, 1)
		test('a - b', assert => testBRepAnd(assert, a, b.flipped(), BRep.EMPTY) )
		test('b - a', assert => testBRepAnd(assert, a.flipped(), b, new BRep(b.faces.concat(a.flipped().faces), false)) )
		test('a + b', assert => testBRepAnd(assert, a.flipped(), b.flipped(), b.flipped()) )
		test('a and b', assert => testBRepAnd(assert, a, b, a) )

		const c = B2T.box(2, 2, 2).translate(20)
		test('a - c', assert => testBRepAnd(assert, a, c.flipped(), a) )
		test('a + c', assert => testBRepAnd(assert, a.flipped(), c.flipped(), new BRep(c.faces.concat(a.faces), false).flipped()) )
		test('a and c', assert => testBRepAnd(assert, a, c, BRep.EMPTY) )
	})
})
