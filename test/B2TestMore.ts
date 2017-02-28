QUnit.module('B2TestMore')

interface Assert {
	dotTest2(face: Face, brep: B2, resultFaces: Face[], desc: string)
	dotTest(face: Face, brep: B2, resultFaces: Face[], desc: string)
}
function doTest(test, face: PlaneFace, brep2: B2, resultEdges: Edge[], resultPoints: V3[], desc?: string) {
	if (brep2 instanceof Face) {
		brep2 = new B2([brep2])
	}
	brep2.buildAdjacencies()
	test.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
						href='../brep2.html?a=${new B2([face]).toSource()}&b=${brep2.toSource()}
						&edges=[${resultEdges.map(e => e.toSource()).join(',')}]
						&points=[${resultPoints.map(e => e.toSource()).join(',')}]'>${desc}</a>`)
	const faceMap = new Map(), edgeMap = new Map(), colinearEdgePairs = new NLA.CustomSet()
	brep2.faces.forEach(face2 => {
		face.intersectPlaneFace(face2, new B2([face]), brep2, faceMap, edgeMap, colinearEdgePairs)
	})
	const edges = faceMap.get(face) || []
	console.log(faceMap)
	test.equal(edges.length, resultEdges.length, resultEdges.length + ' == resultEdges.length == edges.length'+ edges.toSource())
	resultEdges.forEach(edge => {
		test.ok(edges.some(edge2 => edge.like(edge2)), `edges.some(edge2 => edge.like(edge2)) [${edges.toSource()}] ${edge.toSource()}`)
	})
	const uniquePoints = []
	face.edges.forEach(edge => {
		const em = edgeMap.get(edge)
		em && em.forEach(info => info && !uniquePoints.some(up => up.like(info.p)) && assert(info.p) && uniquePoints.push(info.p))
	})
	test.equal(uniquePoints.length, resultPoints.length, resultPoints.length + ' == points.length == resultPoints.length'+ uniquePoints.toSource())
	resultPoints.forEach(p => {
		test.ok(uniquePoints.some(up => up.like(p)), `edges.some(edge2 => edge.like(edge2)) [${uniquePoints.toSource()}]`)
	})
}
function doTestWithBrep(test, face, faceBrep, brep2, resultEdges, resultPoints, desc, backwards?) {
	if (brep2 instanceof Face) {
		brep2 = new B2([brep2])
	}
	faceBrep.buildAdjacencies()
	brep2.buildAdjacencies()
	test.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
						href='../brep2.html?a=${faceBrep.toSource()}&b=${brep2.toSource()}
						&edges=[${resultEdges.map(e => e.toSource()).join(',')}]
						&points=[${resultPoints.map(e => e.toSource()).join(',')}]'>${desc}</a>`)
	const faceMap = new Map(), faceEdgePoints = new Map(), colinearEdgePairs = new NLA.CustomSet()
	if (!backwards) {
		brep2.faces.forEach(face2 => {
			face.intersectPlaneFace(face2, faceBrep, brep2, faceMap, faceEdgePoints, colinearEdgePairs)
		})
	} else {
		brep2.faces.forEach(face2 => {
			face2.intersectPlaneFace(face, brep2, faceBrep, faceMap, faceEdgePoints, colinearEdgePairs)
		})
	}
	const edges = faceMap.get(face) || []
	console.log(faceMap)
	test.equal(edges.length, resultEdges.length, resultEdges.length + ' == resultEdges.length == edges.length'+ edges.toSource())
	resultEdges.forEach(edge => {
		test.ok(edges.some(edge2 => edge.like(edge2)), `edges.some(edge2 => edge.like(edge2)) [${edges.toSource()}] ${edge.toSource()}`)
	})
	const edgePointInfos = faceEdgePoints.get(face), uniquePoints = []
	edgePointInfos && edgePointInfos.forEach(info => !uniquePoints.some(up => up.like(info.p)) && assert(info.p) && uniquePoints.push(info.p))
	console.log('edgePointInfos', edgePointInfos)
	test.equal(uniquePoints.length, resultPoints.length, resultPoints.length + ' == resultPoints.length == uniquePoints.length'+ uniquePoints.toSource())
	resultPoints.forEach(p => {
		test.ok(uniquePoints.some(up => up.like(p)), `uniquePoints.some(up => up.like(p)) [${uniquePoints.toSource()}]`)
	})
}
QUnit.assert.doTest2 = function (face, brep, resultFaces, desc) {
	if (brep instanceof Face) {
		console.log('blah')
		brep = new B2([brep])
	}
	const faceMap = new Map(), edgeMap = new Map()
	const faceBrep = new B2([face])
	this.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
href='../brep2.html?a=${faceBrep.toSource()}&b=${brep.toSource()}&c=${new B2(resultFaces).toSource()}.translate(20, 0, 0)'>${desc}</a>`)
	brep.faces.forEach(face2 => {
		face.intersectPlaneFace(face2, faceBrep, brep, faceMap, edgeMap)
	})
	console.log('faceMap', faceMap)
	const edgeLooseSegments = B2.prototype.getLooseEdgeSegments(edgeMap)
	const newFaces = []
	B2.reconstituteFaces([face], edgeLooseSegments, faceMap, newFaces)
	this.equal(newFaces.length, resultFaces.length, 'number of new faces')
	this.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
href='../brep2.html?a=${faceBrep.toSource()}&b=${brep.toSource()}&c=${new B2(newFaces).toSource()}.translate(20, 0, 0)'>result</a>`)
	resultFaces.forEach(face => {
		this.ok(newFaces.some(newFace => newFace.likeFace(face)), `newFaces.some(newFace => newFace.likeFace(face) ${newFaces.toSource()}`)
	})
}
function doTest3(assert, face: Face, newEdges: Edge[], points: Map<Edge, V3[]>, resultFaces: Face[], desc: string) {
    const faceBrep = new B2([face])
	assert.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
href='../brep2.html?a=${faceBrep.toSource()}&c=${new B2(resultFaces).toSource()}.translate(20, 0, 0)'>${desc}</a>`)
	const isps = Array.from(points.entries()).map(([edge, ps]) =>
		ps.map(p => ({edge: edge, p: p, edgeT: NLA.snap(NLA.snap(edge.curve.pointLambda(p), edge.aT), edge.bT)}))
	).concatenated()
	const edgeLooseSegments = B2.prototype.getLooseEdgeSegments(new Map().set(face, isps))
	const newFaces = []
	B2.reconstituteFaces([face], edgeLooseSegments, new Map().set(face, newEdges), newFaces)
	assert.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank' 
        href='../brep2.html?a=${faceBrep.toSource()}&c=${new B2(newFaces).toSource()}.translate(20, 0, 0)'>result</a>`)
	assert.equal(newFaces.length, resultFaces.length, 'number of new faces')
	resultFaces.forEach(face => {
		assert.ok(newFaces.some(newFace => newFace.likeFace(face)), `newFaces.some(newFace => newFace.likeFace(face) ${newFaces.toSource()}`)
	})
}
registerTests({
	'two faces cutting each other in the middle'(assert) {
		const face = PlaneFace.forVertices(P3.XY, [V(0, 0), V(10, 0), V(10, 10), V(0, 10)])
		console.log(face, Face.prototype.rotateX)
		const face2 = face.rotateX(PI / 2).translate(1, 2, -3)
		doTest(assert, face, face2, [StraightEdge.throughPoints(V(1, 2), V(10, 2))], [V(10, 2)])
	
		const face3 = PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).rotateX(PI / 2).translate(1, 2, -3)
		doTest(assert, face, face3, [StraightEdge.throughPoints(V(1, 2), V(6, 2))], [])
	},
	
	'face touching edge of test face with its middle'(assert) {
		const face = PlaneFace.forVertices(P3.XY, [V(0, 0), V(10, 0), V(10, 10), V(0, 10)])
		const face2 = PlaneFace.forVertices(P3.XY, [V(0, 0), V(12, 0), V(12, 10), V(0, 10)]).rotateX(PI / 2).translate(-1, 0, -5)
		doTest(assert, face, face2, [], [V(0, 0), V(10, 0)]) // StraightEdge.throughPoints(V(0, 0), V(10, 0))
		doTest(assert, face, face2.flipped(), [], [])
	},
	
	'V-shape with spine touching middle of test face'(assert) {
		const face = PlaneFace.forVertices(P3.XY, [V(-10, -10), V(10, -10), V(10, 10), V(-10, 10)])
		const brep = B2T.tetrahedron(V(0, 0, 0), V(5, 0, 0), V(0, 10, 0), V(0, 10, 6))
	
		doTest(assert, face, brep.rotateX(PI / 2), [], [])
		doTest(assert, face, brep.rotateX(PI / 2).flipped(), [], [])
		doTest(assert, face, brep.rotateX(-PI / 2), [], [])
		doTest(assert, face, brep.rotateX(-PI / 2).flipped(), [StraightEdge.throughPoints(V(0, 0), V(5, 0)), StraightEdge.throughPoints(V(5, 0), V(0, 0))], [])
	},
	
	'V-shape splitting test face at spine'(assert) {
		const face = PlaneFace.forVertices(P3.XY, [V(-10, -10), V(10, -10), V(10, 10), V(-10, 10)])
		const brep = new B2([
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).flipped().rotateX(-PI / 4),
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).rotateX(PI / 4)
		])
	
		doTest(assert, face, brep, [StraightEdge.throughPoints(V(0, 0), V(5, 0))], [])
		doTest(assert, face, brep.flipped(), [StraightEdge.throughPoints(V(5, 0), V(0, 0))], [])
	},
	
	'V-shape with spine parallel to test face normal; touching test face edge'(assert) {
		const face = PlaneFace.forVertices(P3.XY, [V(-20, 0), V(20, 0), V(20, 10), V(-20, 10)])
		const brep = new B2([
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).flipped().rotateX(-PI / 4),
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).rotateX(PI / 4)
		]).rotateY(PI / 2).translate(0, 0, 2)
	
		doTest(assert, face, brep, [StraightEdge.throughPoints(V(0, 0), V(5 * sqrt(2), 5 * sqrt(2))), StraightEdge.throughPoints(V(-5 * sqrt(2), 5 * sqrt(2)), V(0, 0))], [])
		doTest(assert, face, brep.flipped(), [
			StraightEdge.throughPoints(V(5 * sqrt(2), 5 * sqrt(2)), V(0, 0)), StraightEdge.throughPoints(V(0, 0), V(-5 * sqrt(2), 5 * sqrt(2)))], [])
		doTest(assert, face, brep.rotateZ(PI), [], [])
		doTest(assert, face, brep.rotateZ(PI).flipped(), [], [])
	},
	
	'V-shape with spine parallel to test face normal; touching test face edge; splitting test face edge'(assert) {
		const face = PlaneFace.forVertices(P3.XY, [V(-20, 0), V(20, 0), V(20, 10), V(-20, 10)])
		const brep = new B2([
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).flipped().rotateX(-PI / 4),
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).rotateX(PI / 4)
		]).rotateY(PI / 2).translate(0, 0, 2)
	
		doTest(assert, face, brep.rotateZ(-PI / 2), [
			StraightEdge.throughPoints(V(5 * sqrt(2), 5 * sqrt(2)), V(0, 0))], [V(0, 0)])
		doTest(assert, face, brep.rotateZ(-PI / 2).flipped(), [
			StraightEdge.throughPoints(V(0, 0), V(5 * sqrt(2), 5 * sqrt(2)))], [V(0, 0)])
	},
	
	'V-shape with spine parallel to test face normal; one wing of V overlaps edge of test face'(assert) {
		const face = PlaneFace.forVertices(P3.XY, [V(-10, -10), V(10, -10), V(10, 10), V(-10, 10)])
		const brep = new B2([
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).flipped().rotateX(-PI / 4),
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)])
		]).rotateY(PI / 2).translate(10, 5, 2)
	
		doTest(assert, face, brep, [
			StraightEdge.throughPoints(V(5, 10), V(10, 5))], [V(10, 5), V(5, 10), V(10, 10)])
		doTest(assert, face, brep.flipped(), [
			StraightEdge.throughPoints(V(10, 5), V(5, 10))], [V(10, 5), V(5, 10)])
	},
	
	'V-shape with spine touching test face, overlapping edge of test face'(assert) {
		const face = PlaneFace.forVertices(P3.XY, [V(-10, -10), V(2, -10), V(2, 10), V(-10, 10)])
		const brep = new B2([
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).flipped(),
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).rotateX(PI / 2)
		])
	
		doTest(assert, face, brep.rotateX(PI / 4), [], [])
		doTest(assert, face, brep.rotateX(PI / 4).flipped(), [], [])
		doTest(assert, face, brep.rotateX(-PI / 4), [StraightEdge.throughPoints(V(0, 0), V(2, 0))], [V(2, 0)])
		doTest(assert, face, brep.rotateX(-PI / 4).flipped(), [StraightEdge.throughPoints(V(2, 0), V(0, 0))], [V(2, 0)])
		doTest(assert, face, brep.rotateX(-PI / 4 * 3), [], [])
		doTest(assert, face, brep.rotateX(-PI / 4 * 3).flipped(), [StraightEdge.throughPoints(V(2, 0), V(0, 0)), StraightEdge.throughPoints(V(0, 0), V(2, 0))], [])
	},
	
	'V-shape with spine touching test face, touching edge of test face'(assert) {
		const face = PlaneFace.forVertices(P3.XY, [V(-10, -10), V(5, -10), V(5, 10), V(-10, 10)])
		const brep = new B2([
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).flipped(),
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).rotateX(PI / 2)
		])
		doTest(assert, face, brep.rotateX(PI / 4), [], [])
		doTest(assert, face, brep.rotateX(PI / 4).flipped(), [], [])
		doTest(assert, face, brep.rotateX(-PI / 4), [StraightEdge.throughPoints(V(0, 0), V(5, 0))], [V(5, 0)])
		doTest(assert, face, brep.rotateX(-PI / 4).flipped(), [StraightEdge.throughPoints(V(5, 0), V(0, 0))], [V(5, 0)])
		doTest(assert, face, brep.rotateX(-PI / 4 * 3), [], [])
		doTest(assert, face, brep.rotateX(-PI / 4 * 3).flipped(), [StraightEdge.throughPoints(V(5, 0), V(0, 0)), StraightEdge.throughPoints(V(0, 0), V(5, 0))], [])
	},
	
	'V-shape with spine touching test face edge middle'(assert) {
		const face = PlaneFace.forVertices(P3.XY, [V(-10, 0), V(10, 0), V(10, 10), V(-10, 10)])
		const brep = new B2([
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).flipped(),
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).rotateX(PI / 2)
		])
	
		doTest(assert, face, brep.rotateX(PI / 4), [], [])
		doTest(assert, face, brep.rotateX(PI / 4).flipped(), [], [V(0, 0), V(5, 0)]) // StraightEdge.throughPoints(V(0, 0),
		                                                                            // V(5, 0)),
		                                                                            // StraightEdge.throughPoints(V(0, 0),
		                                                                            // V(5, 0))
		doTest(assert, face, brep.rotateX(-PI / 4), [], [V(0, 0), V(5, 0)]) // StraightEdge.throughPoints(V(0, 0), V(5, 0)),
		                                                                   // StraightEdge.throughPoints(V(0, 0), V(5, 0))
		doTest(assert, face, brep.rotateX(-PI / 4).flipped(), [], [])
		doTest(assert, face, brep.rotateX(-PI / 4 * 3), [], [])
		doTest(assert, face, brep.rotateX(-PI / 4 * 3).flipped(), [], [V(0, 0), V(5, 0)]) // StraightEdge.throughPoints(V(0,
		                                                                                 // 0), V(5, 0)),
		                                                                                 // StraightEdge.throughPoints(V(0,
		                                                                                 // 0), V(5, 0))
	},
	
	'V-shape with spine overlapping test face edge'(assert) {
		const face = PlaneFace.forVertices(P3.XY, [V(-10, 0), V(3, 0), V(3, 10), V(-10, 10)])
		const brep = new B2([
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).flipped(),
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).rotateX(PI / 2)
		])
	
		doTest(assert, face, brep.rotateX(PI / 4), [], [])
		doTest(assert, face, brep.rotateX(PI / 4).flipped(), [], [V(0, 0), V(3, 0)]) // StraightEdge.throughPoints(V(0, 0),
		                                                                            // V(3, 0)),
		                                                                            // StraightEdge.throughPoints(V(0, 0),
		                                                                            // V(3, 0))
		doTest(assert, face, brep.rotateX(-PI / 4), [], [V(0, 0), V(3, 0)]) // StraightEdge.throughPoints(V(0, 0), V(3, 0)),
		                                                                   // StraightEdge.throughPoints(V(0, 0), V(3, 0))
		doTest(assert, face, brep.rotateX(-PI / 4).flipped(), [], [])
		doTest(assert, face, brep.rotateX(-PI / 4 * 3), [], [])
		doTest(assert, face, brep.rotateX(-PI / 4 * 3).flipped(), [], [V(0, 0), V(3, 0)]) // StraightEdge.throughPoints(V(0,
		                                                                                 // 0), V(3, 0)),
		                                                                                 // StraightEdge.throughPoints(V(0,
		                                                                                 // 0), V(3, 0))
	},
	
	'V-shape with spine parallel to test face normal; touching test face corner'(assert) {
		const face = PlaneFace.forVertices(P3.XY, [V(0, 0), V(10, 0), V(10, 10), V(0, 10)])
		// splitting contour in base position:
		const brep = new B2([
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).flipped().rotateX(-PI / 8),
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).rotateX(PI / 8)
		]).rotateY(PI / 2).translate(0, 0, 2)
	
		doTest(assert, face, brep, [
			StraightEdge.throughPoints(V(0, 0), V3.polar(10, PI / 2 - PI / 8))], [])
	},
	
	'tetrahedorn intersctions doTestWithBrep'(assert) {
		const box = B2T.box(8, 9, 10, 'box')
		const face = box.faces.find(face => face.surface.plane.normal.like(V3.Z))
		let testBrep
		// splitting contour in base position:
	
		testBrep = B2T.tetrahedron(V(-1, 1, 9), V(5, 1, 9), V(-1, -4, 14), V(2, -4, 10))
		doTestWithBrep(assert, face, box, testBrep, [], [],
			'face of tetra touches edge of test face, tetra intersects test volume (test face not part of result, points dont really matter)')
	
		testBrep = B2T.tetrahedron(V(-1, 1, 9), V(4, 1, 9), V(-1, -4, 14), V(-1, -4, 10)).flipped()
		doTestWithBrep(assert, face, box, testBrep, [], [V(3, 0, 10), V(0, 0, 10)],
			'face of _flipped_ tetra touches edge of testface and also intersects main volume (expect point on edge)')
	
		testBrep = B2T.tetrahedron(V(-1, 0, 10), V(5, 0, 10), V(2, 0, 14), V(2, -4, 10))
		doTestWithBrep(assert, face, box, testBrep, [], [],
			'volumes touch edge-edge but no overlap (empty result volume; generated points dont matter)')
	
		testBrep = B2T.tetrahedron(V(-1, 0, 10), V(5, 0, 10), V(-1, -2, 8), V(-1, 2, 8)).flipped()
		doTestWithBrep(assert, face, box, testBrep, [], [V(5, 0, 10), V(0, 0, 10)],
			'Tetrahedron is flipped, testface only touched at edge, needs point on edge as tetrahedron intersects side of volume')
	
	
		testBrep = B2T.tetrahedron(V(-1, 0, 10), V(5, 0, 10), V(-1, 0, 14), V(-1, -4, 9)).flipped()
		doTestWithBrep(assert, face, box, testBrep, [], [],
			'volumes do not intersect, tetra is flipped, touches box edge-edge (result would be entire box, so no points)')
	},
	
	'tetrahedorn is point on line'(assert) {
		const box = B2T.box(8, 9, 10, 'box')
		let face = box.faces.find(face => face.surface.plane.normal.like(V3.Z)), testBrep
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
				StraightEdge.throughPoints(V(3, 3, 10), V(5, 3, 10))], [],
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
	},
	
	'face box'(assert) {
		const box = B2T.box(8, 9, 10, 'box')
		const face = box.faces.find(face => face.surface.plane.normal.like(V3.Z))
		// splitting contour in base position:
	
		const testBrep = B2T.box(8, 1, 1).translate(0, 0, 9).flipped()
		doTestWithBrep(assert, face, box, testBrep, [
				StraightEdge.throughPoints(V(0, 1, 10), V(8, 1, 10))], [V(0, 1, 10), V(8, 1, 10)],
			'???')
	},
	
	'face box 2'(assert) {
		const box = B2T.box(1, 1, 6, 'box').flipped()
		const face = box.faces.find(face => face.surface.plane.normal.like(V3.X))
		// splitting contour in base position:
	
		const testBrep = B2T.box(5, 5, 5)
		doTestWithBrep(assert, face, box, testBrep, [], [],
			'???')
	
	},
	
	'face box 3'(assert) {
		const box = B2T.box(10, 10, 5, 'box')
		const face = box.faces.find(face => face.surface.plane.normal.like(V3.Z))
		// splitting contour in base position:
	
		const testBrep = B2T.box(10, 10, 5, 'box').translate(3, 3, 0)
		doTestWithBrep(assert, face, box, testBrep, [], [],
			'???', true)
		doTestWithBrep(assert, face, box, testBrep, [
				StraightEdge.throughPoints(V(3, 10, 5), V(3, 3, 5)),
				StraightEdge.throughPoints(V(3, 3, 5), V(10, 3, 5))], [V(3, 10, 5), V(10, 3, 5)],
			'???')
	
	},

    'face box 4'(assert) {
        const box = B2T.box(10, 10, 5, 'box').flipped()
        const face = box.faces.find(face => (face.surface as PlaneSurface).plane.normal.z == -1)
        // splitting contour in base position:

        const testBrep = box.translate(3, 3, 0)
        doTestWithBrep(assert, face, box, testBrep, [], [V(3, 10, 5), V(10, 3, 5)],
            '???')
        doTestWithBrep(assert, face, box, testBrep, [
                StraightEdge.throughPoints(V(3, 10, 5), V(3, 3, 5)),
                StraightEdge.throughPoints(V(3, 3, 5), V(10, 3, 5))], [V(3, 10, 5), V(10, 3, 5)],
            '???', true)

    },
    'face box 5'(assert) {
        const box = B2T.box(10, 10, 5, 'box').flipped()
        const face = box.faces.find(face => (face.surface as PlaneSurface).plane.normal.z == -1)
        // splitting contour in base position:

        const testBrep = box.translate(0, 0, 5)
        doTestWithBrep(assert, face, box, testBrep, [], [],
            '???')
        doTestWithBrep(assert, face, box, testBrep, [], [],
            '???', true)

    },
    'face box 6'(assert) {
        const box = B2T.box(10, 10, 10, 'box').flipped()
        const box2 = B2T.box(10, 12, 12, 'box').flipped()

        const face = box.faces.find(face => (face.surface as PlaneSurface).plane.normal.z == -1)
        // splitting contour in base position:

        const testBrep = box2.translate(0, 0, -2)
        doTestWithBrep(assert, face, box, testBrep, [], [],
            '???')
        doTestWithBrep(assert, face, box, testBrep, [], [],
            '???', true)

    },

	'coplanar things'(assert) {
		const box = B2T.box(8, 9, 10, 'box')
		let face = box.faces.find(face => face.surface.plane.normal.like(V3.Z)), testBrep
		// splitting contour in base position:
	
		testBrep = B2T.tetrahedron(V(-1, -1, 10), V(-1, 5, 10), V(5, -1, 10), V(-1, -1, 5)).flipped()
		doTestWithBrep(assert, face, box, testBrep, [StraightEdge.throughPoints(V(0, 4, 10), V(4, 0, 10))], [V(0, 4, 10), V(4, 0, 10)],
			'cut off corner of box with flipped tetra (anti)coplanar to test face')
	
		testBrep = B2T.tetrahedron(V(-1, -1, 10), V(-1, 5, 10), V(4, 0, 10), V(-1, -1, 5)).flipped()
		doTestWithBrep(assert, face, box, testBrep, [StraightEdge.throughPoints(V(0, 4, 10), V(4, 0, 10))], [V(0, 4, 10), V(4, 0, 10)],
			'cut off corner of box with flipped tetra (anti)coplanar to test face')
	
		testBrep = B2T.tetrahedron(V(1, 0, 10), V(1, 5, 10), V(4, 0, 10), V(-1, -1, 5)).flipped()
		doTestWithBrep(assert, face, box, testBrep,
			[StraightEdge.throughPoints(V(1, 0, 10), V(1, 5, 10)), StraightEdge.throughPoints(V(1, 5, 10), V(4, 0, 10))],
			[V(1, 0, 10), V(4, 0, 10)],
			'cut off corner of box with flipped tetra (anti)coplanar to test face')
	
		testBrep = B2T.tetrahedron(V(0, 0, 10), V(0, 5, 10), V(5, 0, 10), V(0, 0, 5)).flipped()
		doTestWithBrep(assert, face, box, testBrep,
			[StraightEdge.throughPoints(V(0, 5, 10), V(5, 0, 10))],
			[V(0, 5, 10), V(5, 0, 10), V(0, 0, 10)],
			'cut off corner of box with flipped tetra (anti)coplanar to test face')
	
		testBrep = B2T.tetrahedron(V(0, 0, 10), V(-1, 5, 10), V(5, -1, 10), V(-1, -1, 5)).flipped()
		doTestWithBrep(assert, face, box, testBrep, [StraightEdge.throughPoints(V(0, 4, 10), V(4, 0, 10))], [V(0, 4, 10), V(4, 0, 10)],
			'cut off corner of box with flipped tetra (anti)coplanar to test face')
	
		testBrep = B2T.tetrahedron(V(0, 0, 10), V(1, 5, 10), V(5, 1, 10), V(0, 0, 5)).flipped()
		doTestWithBrep(assert, face, box, testBrep,
			[StraightEdge.throughPoints(V(0, 0, 10), V(1, 5, 10)), StraightEdge.throughPoints(V(1, 5, 10), V(5, 1, 10)), StraightEdge.throughPoints(V(5, 1, 10), V(0, 0, 10))],
			[],
			'cut hole at corner of test face')
	
	},
	'coplanar things 2'(assert) {
		const box = B2T.box(8, 9, 10, 'box')
		let face = box.faces.find(face => face.surface.plane.normal.like(V3.Z)), testBrep
		// splitting contour in base position:
	
		testBrep = B2T.tetrahedron(V(-1, -1, 10), V(-1, 5, 10), V(5, -1, 10), V(-1, -1, 5))
		doTestWithBrep(assert, face, box, testBrep, [StraightEdge.throughPoints(V(0, 4, 10), V(4, 0, 10))], [V(0, 4, 10), V(4, 0, 10)],
			'cut off corner of box with flipped tetra (anti)coplanar to test face')
	
		testBrep = B2T.tetrahedron(V(-1, -1, 10), V(-1, 5, 10), V(5, -1, 10), V(-1, -1, 5)).flipped()
		doTestWithBrep(assert, face, box.flipped(), testBrep, [StraightEdge.throughPoints(V(0, 4, 10), V(4, 0, 10))], [V(0, 4, 10), V(4, 0, 10)],
			'cut off corner of box with flipped tetra (anti)coplanar to test face')
	
	},
	
	
	'assembly'(assert) {
		const baseFace = PlaneFace.forVertices(P3.XY, [V(0, 0), V(10, 0), V(10, 10), V(0, 10)])
	
		const extrude = B2T.extrudeVertices([V(5, -1), V(2, 2), V(8, 2)], P3.XY.flipped(), V(0, 0, 10)).translate(0, 0, -2).flipped()
		const result = PlaneFace.forVertices(P3.XY, [V(0, 0), V(4, 0), V(2, 2), V(8, 2), V(6, 0), V(10, 0), V(10, 10), V(0, 10)])
		assert.doTest2(baseFace, extrude, [result], 'volume cuts edge of test face (twice)')
	
		const edges = StraightEdge.chain([V(5, 0), V(2, 3), V(8, 3)])
		const result2 = PlaneFace.forVertices(P3.XY, [V(0, 0), V(10, 0), V(10, 10), V(0, 10)], [V(5, 0), V(2, 3), V(8, 3)])
		doTest3(assert, baseFace, edges, new Map(), [result2], 'volume touches inside of test face edge')
	
	
		// from test case 3:4
		// V-shape spine touching middle, splits volume enclosing in both directions
		let brep = new B2([
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).flipped(),
			PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 0), V(5, 10), V(0, 10)]).rotateX(PI / 4)
		]).rotateX(-PI / 2).flipped().translate(0, 2, 0)
		// degenerate cycle in middle of face
		let result4 = PlaneFace.forVertices(P3.XY, [V(0, 0), V(10, 0), V(10, 10), V(0, 10)], [V(5, 2), V(0, 2)])
	
	
		const extrude5 = B2T.extrudeVertices([V(0, 0), V(3, 2), V(2, 3)], P3.XY.flipped(), V(0, 0, 10)).translate(0, 0, -2).flipped()
		const result5 = PlaneFace.forVertices(P3.XY, [V(0, 0), V(2, 3), V(3, 2), V(0, 0), V(10, 0), V(10, 10), V(0, 10)])
		assert.doTest2(baseFace, extrude5, [result5], 'volume touches inside of test face corner')
	
		const edges6 = StraightEdge.chain([V(1, 0), V(8, 10), V(3, 0)], false)
		const points6 = new Map().set(baseFace.edges.find(edge => edge.aDir.like(V3.X)), [V(1, 0), V(3, 0)])
		const result6 = [PlaneFace.forVertices(P3.XY, [V(0, 0), V(1, 0), V(8, 10), V(3, 0), V(10, 0), V(10, 10), V(0, 10)])]
		doTest3(assert, baseFace, edges6, points6, result6, 'volume touches inside of test face corner')

		const edges7 = [StraightEdge.throughPoints(V(1, 0), V(7, 10)), StraightEdge.throughPoints(V(8, 10), V(3, 0))]
		const points7 = new Map()
			.set(baseFace.edges.find(edge => edge.aDir.like(V3.X)), [V(1, 0), V(3, 0)])
			.set(baseFace.edges.find(edge => edge.aDir.like(V3.X.negated())), [V(7, 10), V(8, 10)])
		const result7 = [PlaneFace.forVertices(P3.XY, [V(0, 0), V(1, 0), V(7, 10), V(0, 10)]), PlaneFace.forVertices(P3.XY, [V(10, 10), V(8, 10), V(3, 0), V(10, 0)])]
		doTest3(assert, baseFace, edges7, points7, result7, 'volume touches inside of test face corner')


		const edges8 = [StraightEdge.throughPoints(V(1, 0), V(8, 10)), StraightEdge.throughPoints(V(8, 10), V(1, 0))]
		const points8 = new Map()
			.set(baseFace.edges.find(edge => edge.aDir.like(V3.X)), [V(1, 0)])
			.set(baseFace.edges.find(edge => edge.aDir.like(V3.X.negated())), [V(8, 10)])
		const result8 = [PlaneFace.forVertices(P3.XY, [V(0, 0), V(1, 0), V(8, 10), V(0, 10)]), PlaneFace.forVertices(P3.XY, [V(10, 10), V(8, 10), V(1, 0), V(10, 0)])]
		doTest3(assert, baseFace, edges8, points8, result8, 'volume touches inside of test face corner')


		const edges9 = [StraightEdge.throughPoints(V(1, 0), V(8, 5)), StraightEdge.throughPoints(V(8, 5), V(1, 0))]
		const points9 = new Map().set(baseFace.edges.find(edge => edge.aDir.like(V3.X)), [V(1, 0)])
		const result9 = [PlaneFace.forVertices(P3.XY, [V(0, 0), V(1, 0), V(8, 5), V(1, 0), V(10, 0), V(10, 10), V(0, 10)])]
		doTest3(assert, baseFace, edges9, points9, result9, 'volume touches inside of test face corner')


	},
	
	'assembly holes'(assert) {
		const baseFace = PlaneFace.forVertices(P3.XY, [V(0, 0), V(10, 0), V(10, 10), V(0, 10)])
	
		const extrude7 = B2T.extrudeVertices([V(0, 0), V(5, 8), V(8, 5)], P3.XY.flipped(), V(0, 0, 10)).translate(0, 0, -2).flipped()
		const result7 = [PlaneFace.forVertices(P3.XY, [V(0, 0), V(5, 8), V(8, 5), V(0, 0), V(10, 0), V(10, 10), V(0, 10)])]
		assert.doTest2(baseFace, extrude7, result7, 'volume touches inside of test face corner')
	
		const extrude8 = B2T.extrudeVertices([V(1, 1), V(1, -1), V(-1, 1)], P3.XY.flipped(), V(0, 0, 10)).translate(0, 0, -2).flipped()
		const result8 = [PlaneFace.forVertices(P3.XY, [V(1, 0), V(10, 0), V(10, 10), V(0, 10), V(0, 1), V(1, 1)])]
		assert.doTest2(baseFace, extrude8, result8, 'volume touches inside of test face corner')
	
		const extrude9 = B2T.extrudeVertices([V(-1, -1), V(1, 1), V(1, -1)], P3.XY.flipped(), V(0, 0, 10)).translate(0, 0, -2).flipped()
		const result9 = [PlaneFace.forVertices(P3.XY, [V(1, 0), V(10, 0), V(10, 10), V(0, 10), V(0, 0), V(1, 1)])]
		assert.doTest2(baseFace, extrude9, result9, 'volume touches inside of test face corner')
	
		const extrude10 = B2T.extrudeVertices([V(1, 1), V(2, 2), V(2, 1)], P3.XY.flipped(), V(0, 0, 10)).translate(0, 0, -2).flipped()
		const result10 = [PlaneFace.forVertices(P3.XY, [V(10, 0), V(10, 10), V(0, 10), V(0, 0)], [V(1, 1), V(2, 2), V(2, 1)])]
		assert.doTest2(baseFace, extrude10, result10, 'adding a hole')
	
		const baseFace11 = PlaneFace.forVertices(P3.XY, [V(10, 0), V(10, 10), V(0, 10), V(0, 0)], [V(1, 1), V(2, 2), V(2, 1)])
		const extrude11 = B2T.extrudeVertices([V(5, 5), V(6, 6), V(6, 5)], P3.XY.flipped(), V(0, 0, 10)).translate(0, 0, -2).flipped()
		const result11 = [PlaneFace.forVertices(P3.XY, [V(10, 0), V(10, 10), V(0, 10), V(0, 0)], [V(1, 1), V(2, 2), V(2, 1)], [V(5, 5), V(6, 6), V(6, 5)])]
		assert.doTest2(baseFace11, extrude11, result11, 'adding a hole to a face which already has one')
	
	
		let baseFace12 = PlaneFace.forVertices(P3.XY, [V(10, 0), V(10, 10), V(0, 10), V(0, 0)], [V(1, 1), V(5, 5), V(5, 1)])
		const extrude12 = B2T.extrudeVertices([V(2, 1.5), V(2, 4), V(4.5, 4)], P3.XY.flipped(), V(0, 0, 10)).translate(0, 0, -2).flipped()
		const result12 = [PlaneFace.forVertices(P3.XY, [V(10, 0), V(10, 10), V(0, 10), V(0, 0)], [V(5, 5), V(5, 1), V(1, 1), V(2, 2), V(2, 4), V(4, 4)])]
		assert.doTest2(baseFace12, extrude12, result12, 'extending an existing hole')
	
		const baseFace13 = PlaneFace.forVertices(P3.XY, [V(10, 0), V(10, 10), V(0, 10), V(0, 0)], [V(1, 1), V(5, 5), V(5, 1)])
		const extrude13 = B2T.extrudeVertices([V(3, -1), V(4, -1), V(4, 2), V(3, 2)], P3.XY.flipped(), V(0, 0, 10)).translate(0, 0, -2).flipped()
		const result13 = [PlaneFace.forVertices(P3.XY, [V(10, 0), V(10, 10), V(0, 10), V(0, 0), V(3, 0), V(3, 1), V(1, 1), V(5, 5), V(5, 1), V(4, 1), V(4, 0)])]
		assert.doTest2(baseFace13, extrude13, result13, 'removing a hole by cutting a channel')
	
		const baseFace14 = PlaneFace.forVertices(P3.XY, [V(10, 0), V(10, 10), V(0, 10), V(0, 0)], [V(1, 1), V(5, 5), V(5, 1)])
		const edges14 = StraightEdge.chain([V(1, 1), V(1, 5), V(5, 5)], false)
		const points14 = new Map()
		const result14 = [PlaneFace.forVertices(P3.XY, [V(10, 0), V(10, 10), V(0, 10), V(0, 0)], [V(5, 5), V(5, 1), V(1, 1), V(1, 5)])]
		doTest3(assert, baseFace14, edges14, points14, result14, 'extending an existing hole')
	},
	
	'B2.prototype.minus remove half of a half-pie'(assert) {
		let pie = B2T.puckman(8, 180 * DEG, 5, 'pie/2')
		let boxKnife = B2T.box(11, 10, 7, 'knife').translate(-10, -1, -1)
	
		let resultTopPoint = V(1, 8 * Math.sin(Math.acos(1 / 8)), 0)
		let result = B2T.extrudeEdges([
				StraightEdge.throughPoints(V(8, 0, 0), V(1, 0, 0)),
				StraightEdge.throughPoints(V(1, 0, 0), resultTopPoint),
				PCurveEdge.forCurveAndTs(EllipseCurve.circle(8), Math.acos(1 / 8), 0)],
			P3.XY.flipped(), V(0, 0, 5), 'pie/4')
		assert.b2Equal(pie, boxKnife, pie.minus(boxKnife), result)
	
	
		{
			const testFace = pie.faces.find(face => face.surface.plane && face.surface.plane.normal.like(V3.Y.negated()) && face.edges.some(edge => edge.a.x > 1))
			const k2 = boxKnife.translate(-1, 0, 0)
			doTestWithBrep(assert, testFace, pie, k2.flipped(), [], [V3.ZERO, V(0, 0, 5)],
				'volumes touch edge-edge but no overlap (empty result volume; generated points dont matter)')
			assert.b2Equal(pie, k2, pie.minus(k2), B2T.puckman(8, 90 * DEG, 5, 'pie/4'))
		}
	},
	'B2.prototype.minus cut hole through side of pie'(assert) {
		let pie = B2T.puckman(8, 180 * DEG, 5, 'pie/2')
		let punch = B2T.box(5, 10, 3, 'knife').translate(1, -1, 1)
		let result = new B2([
			new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0)), [
				StraightEdge.throughPoints(V(0, 0, 0), V(8, 0, 0)),
				StraightEdge.throughPoints(V(8, 0, 0), V(8, 0, 5)),
				StraightEdge.throughPoints(V(8, 0, 5), V(0, 0, 5)),
				StraightEdge.throughPoints(V(0, 0, 5), V(0, 0, 0))], [[
				StraightEdge.throughPoints(V(1, 0, 1), V(1, 0, 4)),
				StraightEdge.throughPoints(V(1, 0, 4), V(6, 0, 4)),
				StraightEdge.throughPoints(V(6, 0, 4), V(6, 0, 1)),
				StraightEdge.throughPoints(V(6, 0, 1), V(1, 0, 1))]]),
			new RotationFace(new CylinderSurface(new EllipseCurve(V(0, 0, 0), V(8, 0, 0), V(0, -8, 0)), V(0, 0, -1)), [
				new PCurveEdge(new EllipseCurve(V(0, 0, 0), V(8, 0, 0), V(0, -8, 0)), V(8, 0, 0), V(-8, 9.797174393178826e-16, 0), 0, -3.141592653589793, null, V(0, 8, 0), V(-9.797174393178826e-16, -8, 0)),
				StraightEdge.throughPoints(V(-8, 9.797174393178826e-16, 0), V(-8, 9.797174393178826e-16, 5)),
				new PCurveEdge(new EllipseCurve(V(0, 0, 5), V(8, 0, 0), V(0, -8, 0)), V(-8, 9.797174393178826e-16, 5), V(8, 0, 5), -3.141592653589793, 0, null, V(9.797174393178826e-16, 8, 0), V(0, -8, 0)),
				StraightEdge.throughPoints(V(8, 0, 5), V(8, 0, 0))], [[
				StraightEdge.throughPoints(V(1, 7.937253933193773, 4), V(1, 7.937253933193773, 1)),
				new PCurveEdge(new EllipseCurve(V(0, 0, 1), V(8, 0, 0), V(0, -8, 0)), V(1, 7.937253933193773, 1), V(6, 5.291502622129181, 1), -1.4454684956268313, -0.7227342478134156, null, V(7.937253933193772, -0.9999999999999991, 0), V(5.2915026221291805, -6, 0)),
				StraightEdge.throughPoints(V(6, 5.291502622129181, 1), V(6, 5.291502622129181, 4)),
				new PCurveEdge(new EllipseCurve(V(0, 0, 4), V(8, 0, 0), V(0, 8, 0)), V(6, 5.291502622129181, 4), V(1, 7.937253933193773, 4), 0.7227342478134156, 1.4454684956268313, null, V(-5.2915026221291805, 6, 0), V(-7.937253933193772, 0.9999999999999991, 0))]]),
			new PlaneFace(new PlaneSurface(new P3(V(-1.2246467991473532e-16, -1, 0), 0)), [
				StraightEdge.throughPoints(V(-8, 9.797174393178826e-16, 0), V(0, 0, 0)),
				StraightEdge.throughPoints(V(0, 0, 0), V(0, 0, 5)),
				StraightEdge.throughPoints(V(0, 0, 5), V(-8, 9.797174393178826e-16, 5)),
				StraightEdge.throughPoints(V(-8, 9.797174393178826e-16, 5), V(-8, 9.797174393178826e-16, 0))], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 0)), [
				StraightEdge.throughPoints(V(8, 0, 0), V(0, 0, 0)),
				StraightEdge.throughPoints(V(0, 0, 0), V(-8, 9.797174393178826e-16, 0)),
				new PCurveEdge(new EllipseCurve(V(0, 0, 0), V(8, 0, 0), V(0, -8, 0)), V(-8, 9.797174393178826e-16, 0), V(8, 0, 0), -3.141592653589793, 0, null, V(9.797174393178826e-16, 8, 0), V(0, -8, 0))], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 5)), [
				new PCurveEdge(new EllipseCurve(V(0, 0, 5), V(8, 0, 0), V(0, -8, 0)), V(8, 0, 5), V(-8, 9.797174393178826e-16, 5), 0, -3.141592653589793, null, V(0, 8, 0), V(-9.797174393178826e-16, -8, 0)),
				StraightEdge.throughPoints(V(-8, 9.797174393178826e-16, 5), V(0, 0, 5)),
				StraightEdge.throughPoints(V(0, 0, 5), V(8, 0, 5))], []),
			new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 1)), [
				StraightEdge.throughPoints(V(1, 0, 4), V(1, 0, 1)),
				StraightEdge.throughPoints(V(1, 0, 1), V(1, 7.937253933193773, 1)),
				StraightEdge.throughPoints(V(1, 7.937253933193773, 1), V(1, 7.937253933193773, 4)),
				StraightEdge.throughPoints(V(1, 7.937253933193773, 4), V(1, 0, 4))], []),
			new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), -6)), [
				StraightEdge.throughPoints(V(6, 0, 1), V(6, 0, 4)),
				StraightEdge.throughPoints(V(6, 0, 4), V(6, 5.291502622129181, 4)),
				StraightEdge.throughPoints(V(6, 5.291502622129181, 4), V(6, 5.291502622129181, 1)),
				StraightEdge.throughPoints(V(6, 5.291502622129181, 1), V(6, 0, 1))], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 1)), [
				StraightEdge.throughPoints(V(1, 0, 1), V(6, 0, 1)),
				StraightEdge.throughPoints(V(6, 0, 1), V(6, 5.291502622129181, 1)),
				new PCurveEdge(new EllipseCurve(V(0, 0, 1), V(8, 0, 0), V(0, -8, 0)), V(6, 5.291502622129181, 1), V(1, 7.937253933193773, 1), -0.7227342478134156, -1.4454684956268313, null, V(-5.2915026221291805, 6, 0), V(-7.937253933193772, 0.9999999999999991, 0)),
				StraightEdge.throughPoints(V(1, 7.937253933193773, 1), V(1, 0, 1))], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), -4)), [
				StraightEdge.throughPoints(V(6, 0, 4), V(1, 0, 4)),
				StraightEdge.throughPoints(V(1, 0, 4), V(1, 7.937253933193773, 4)),
				new PCurveEdge(new EllipseCurve(V(0, 0, 4), V(8, 0, 0), V(0, 8, 0)), V(1, 7.937253933193773, 4), V(6, 5.291502622129181, 4), 1.4454684956268313, 0.7227342478134156, null, V(7.937253933193772, -0.9999999999999991, 0), V(5.2915026221291805, -6, 0)),
				StraightEdge.throughPoints(V(6, 5.291502622129181, 4), V(6, 0, 4))], [])], false)
	
		assert.b2Equal(pie, punch, pie.minus(punch), result)
		console.log(pie.minus(punch).sce)
	},
	
	'B2.prototype.minus including ProjectedCurveSurface'(assert) {
		let s1 = new ProjectedCurveSurface(new BezierCurve(V(0, 0, 0), V(-5, 5, 0), V(15, 5, 0), V(10, 0, 0), -0.1, 1.1), V(0, 0, -1), 0, 1)
		let s2 = new ProjectedCurveSurface(new BezierCurve(V(0, 0, 0), V(-5, 5, 0), V(15, 5, 0), V(10, 0, 0), -0.1, 1.1), V(0, 0, -1), 0, 1)
		let a = B2T.extrudeEdges([
				StraightEdge.throughPoints(V(10, 0, 0), V(0, 0, 0)),
				PCurveEdge.forCurveAndTs(new BezierCurve(V3.ZERO, V(-5, 5, 0), V(15, 5, 0), V(10, 0, 0)), 0, 1)],
			P3.XY.flipped(), V(0, 0, 5), 'a/4')
		let punch = B2T.box(5, 10, 3, 'knife').translate(1, -1, 1)
	
		let result = new B2([
			new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0)), [
				StraightEdge.throughPoints(V(0, 0, 0), V(10, 0, 0)),
				StraightEdge.throughPoints(V(10, 0, 0), V(10, 0, 5)),
				StraightEdge.throughPoints(V(10, 0, 5), V(0, 0, 5)),
				StraightEdge.throughPoints(V(0, 0, 5), V(0, 0, 0))], [[
				StraightEdge.throughPoints(V(1, 0, 1), V(1, 0, 4)),
				StraightEdge.throughPoints(V(1, 0, 4), V(6, 0, 4)),
				StraightEdge.throughPoints(V(6, 0, 4), V(6, 0, 1)),
				StraightEdge.throughPoints(V(6, 0, 1), V(1, 0, 1))]]),
			new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(0, 0, 0), V(-5, 5, 0), V(15, 5, 0), V(10, 0, 0), -0.1, 1.1), V(0, 0, -1), 0, 1), [
				new PCurveEdge(new BezierCurve(V(0, 0, 0), V(-5, 5, 0), V(15, 5, 0), V(10, 0, 0), -0.1, 1.1), V(10, 0, 0), V(0, 0, 0), 1, 0, null, V(15, 15, 0), V(15, -15, 0)),
				StraightEdge.throughPoints(V(0, 0, 0), V(0, 0, 5)),
				new PCurveEdge(new BezierCurve(V(0, 0, 5), V(-5, 5, 5), V(15, 5, 5), V(10, 0, 5), -0.1, 1.1), V(0, 0, 5), V(10, 0, 5), 0, 1, null, V(-15, 15, 0), V(-15, -15, 0)),
				StraightEdge.throughPoints(V(10, 0, 5), V(10, 0, 0))], [[
				StraightEdge.throughPoints(V(1, 3.185436104380006, 4), V(1, 3.185436104380006, 1)),
				new PCurveEdge(new BezierCurve(V(0, 0, 1), V(-5, 5, 1), V(15, 5, 1), V(10, 0, 1), -0.1, 1.1), V(1, 3.185436104380006, 1), V(6, 3.720106174228432, 1), 0.3059958942668147, 0.5446421518086273, null, V(16.85436104380006, 5.820123171995558, 0), V(22.201061742284324, -1.3392645542588184, 0)),
				StraightEdge.throughPoints(V(6, 3.720106174228432, 1), V(6, 3.720106174228432, 4)),
				new PCurveEdge(new BezierCurve(V(10, 0, 4), V(15, 5, 4), V(-5, 5, 4), V(0, 0, 4), -1.1, 0.1), V(6, 3.720106174228432, 4), V(1, 3.185436104380006, 4), 0.45535784819137265, 0.6940041057331853, null, V(-22.201061742284324, 1.3392645542588197, 0), V(-16.85436104380006, -5.820123171995558, 0))]]),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 0)), [
				StraightEdge.throughPoints(V(10, 0, 0), V(0, 0, 0)),
				new PCurveEdge(new BezierCurve(V(0, 0, 0), V(-5, 5, 0), V(15, 5, 0), V(10, 0, 0), -0.1, 1.1), V(0, 0, 0), V(10, 0, 0), 0, 1, null, V(-15, 15, 0), V(-15, -15, 0))], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 5)), [
				new PCurveEdge(new BezierCurve(V(0, 0, 5), V(-5, 5, 5), V(15, 5, 5), V(10, 0, 5), -0.1, 1.1), V(10, 0, 5), V(0, 0, 5), 1, 0, null, V(15, 15, 0), V(15, -15, 0)),
				StraightEdge.throughPoints(V(0, 0, 5), V(10, 0, 5))], []),
			new PlaneFace(new PlaneSurface(new P3(V(1, 0, 0), 1)), [
				StraightEdge.throughPoints(V(1, 0, 4), V(1, 0, 1)),
				StraightEdge.throughPoints(V(1, 0, 1), V(1, 3.185436104380006, 1)),
				StraightEdge.throughPoints(V(1, 3.185436104380006, 1), V(1, 3.185436104380006, 4)),
				StraightEdge.throughPoints(V(1, 3.185436104380006, 4), V(1, 0, 4))], []),
			new PlaneFace(new PlaneSurface(new P3(V(-1, 0, 0), -6)), [
				StraightEdge.throughPoints(V(6, 0, 1), V(6, 0, 4)),
				StraightEdge.throughPoints(V(6, 0, 4), V(6, 3.720106174228432, 4)),
				StraightEdge.throughPoints(V(6, 3.720106174228432, 4), V(6, 3.720106174228432, 1)),
				StraightEdge.throughPoints(V(6, 3.720106174228432, 1), V(6, 0, 1))], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, 1), 1)), [
				StraightEdge.throughPoints(V(1, 0, 1), V(6, 0, 1)),
				StraightEdge.throughPoints(V(6, 0, 1), V(6, 3.720106174228432, 1)),
				new PCurveEdge(new BezierCurve(V(0, 0, 1), V(-5, 5, 1), V(15, 5, 1), V(10, 0, 1), -0.1, 1.1), V(6, 3.720106174228432, 1), V(1, 3.185436104380006, 1), 0.5446421518086273, 0.3059958942668147, null, V(-22.201061742284324, 1.3392645542588184, 0), V(-16.85436104380006, -5.820123171995558, 0)),
				StraightEdge.throughPoints(V(1, 3.185436104380006, 1), V(1, 0, 1))], []),
			new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), -4)), [
				StraightEdge.throughPoints(V(6, 0, 4), V(1, 0, 4)),
				StraightEdge.throughPoints(V(1, 0, 4), V(1, 3.185436104380006, 4)),
				new PCurveEdge(new BezierCurve(V(10, 0, 4), V(15, 5, 4), V(-5, 5, 4), V(0, 0, 4), -1.1, 0.1), V(1, 3.185436104380006, 4), V(6, 3.720106174228432, 4), 0.6940041057331853, 0.45535784819137265, null, V(16.85436104380006, 5.820123171995558, 0), V(22.201061742284324, -1.3392645542588197, 0)),
				StraightEdge.throughPoints(V(6, 3.720106174228432, 4), V(6, 0, 4))], [])], false)
	
		assert.b2Equal(a, punch, a.minus(punch), result)
	}
})