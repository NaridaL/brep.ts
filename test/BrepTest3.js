QUnit.module('brep3')

/**
 *
 * @param {Face} face
 * @param {B2} brep2
 * @param {Edge[]} resultEdges
 * @param {V3[]} resultPoints
 * @param {string=} desc
 */
QUnit.assert.doTest = function (face, brep2, resultEdges, resultPoints, desc) {
	if (brep2 instanceof Face) {
		brep2 = new B2([brep2])
	}
	this.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
						href='file:///C:/Users/aval/Desktop/cs/brep2.html?a=${new B2([face]).toSource()}&b=${brep2.toSource()}
						&edges=[${resultEdges.map(e => e.toSource()).join(',')}]
						&points=[${resultPoints.map(e => e.toSource()).join(',')}]'>${desc}</a>`)
	var faceMap = new Map(), edgeMap = new Map()
	brep2.faces.forEach(face2 => {
		face.intersectPlaneFace(face2, new B2([face]), brep2, faceMap, edgeMap)
	})
	var edges = faceMap.get(face) || []
	console.log(faceMap)
	this.equal(edges.length, resultEdges.length, 'resultEdges.length == edges.length'+edges.toSource())
	resultEdges.forEach(edge => {
		this.ok(edges.some(edge2 => edge.like(edge2)), `edges.some(edge2 => edge.like(edge2)) [${edges.toSource()}] ${edge.toSource()}`)
	})
	var uniquePoints = []
	face.edges.forEach(edge => {
		var em = edgeMap.get(edge)
		em && em.forEach(info => info && !uniquePoints.some(up => up.like(info.p)) && assert(info.p) && uniquePoints.push(info.p))
	})
	this.equal(uniquePoints.length, resultPoints.length, 'points.length == resultPoints.length'+uniquePoints.toSource())
	resultPoints.forEach(p => {
		this.ok(uniquePoints.some(up => up.like(p)), `edges.some(edge2 => edge.like(edge2)) [${uniquePoints.toSource()}]`)
	})
}
QUnit.assert.doTestWithBrep = function (face, faceBrep, brep2, resultEdges, resultPoints, desc) {
	if (brep2 instanceof Face) {
		brep2 = new B2([brep2])
	}
	this.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
						href='file:///C:/Users/aval/Desktop/cs/brep2.html?a=${faceBrep.toSource()}&b=${brep2.toSource()}
						&edges=[${resultEdges.map(e => e.toSource()).join(',')}]
						&points=[${resultPoints.map(e => e.toSource()).join(',')}]'>${desc}</a>`)
	var faceMap = new Map(), edgeMap = new Map()
	brep2.faces.forEach(face2 => {
		face.intersectPlaneFace(face2, faceBrep, brep2, faceMap, edgeMap)
	})
	var edges = faceMap.get(face) || []
	console.log(faceMap)
	this.equal(edges.length, resultEdges.length, 'resultEdges.length == edges.length'+edges.toSource())
	resultEdges.forEach(edge => {
		this.ok(edges.some(edge2 => edge.like(edge2)), `edges.some(edge2 => edge.like(edge2)) [${edges.toSource()}] ${edge.toSource()}`)
	})
	var uniquePoints = []
	face.edges.forEach(edge => {
		var em = edgeMap.get(edge)
		em && em.forEach(info => info && !uniquePoints.some(up => up.like(info.p)) && assert(info.p) && uniquePoints.push(info.p))
	})
	this.equal(uniquePoints.length, resultPoints.length, 'points.length == resultPoints.length'+uniquePoints.toSource())
	resultPoints.forEach(p => {
		this.ok(uniquePoints.some(up => up.like(p)), `uniquePoints.some(up => up.like(p)) [${uniquePoints.toSource()}]`)
	})
}

QUnit.assert.doTest2 = function (face, brep, resultFaces, desc) {
	if (brep instanceof Face) {
		console.log("blah")
		brep = new B2([brep])
	}
	var faceMap = new Map(), edgeMap = new Map()
	var faceBrep = new B2([face]);
	this.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
href='file:///C:/Users/aval/Desktop/cs/brep2.html?a=${faceBrep.toSource()}&b=${brep.toSource()}&c=${new B2(resultFaces).toSource()}.translate(20, 0, 0)'>${desc}</a>`)
	brep.faces.forEach(face2 => {
		face.intersectPlaneFace(face2, faceBrep, brep, faceMap, edgeMap)
	})
	console.log('faceMap', faceMap)
	var edgeLooseSegments = B2.prototype.getLooseEdgeSegments(edgeMap)
	var newFaces = []
	B2.prototype.reconstituteFaces([face], edgeLooseSegments, faceMap, newFaces)
	this.equal(newFaces.length, resultFaces.length, "number of new faces")
	this.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
href='file:///C:/Users/aval/Desktop/cs/brep2.html?a=${faceBrep.toSource()}&b=${brep.toSource()}&c=${new B2(newFaces).toSource()}.translate(20, 0, 0)'>result</a>`)
	resultFaces.forEach(face => {
		this.ok(newFaces.some(newFace => newFace.likeFace(face)), `newFaces.some(newFace => newFace.likeFace(face) ${newFaces.toSource()}`)
	})
}
QUnit.test( "two faces cutting each other in the middle", function( assert ) {
	var face = PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(10, 0), V3(10, 10), V3(0, 10)])
	console.log(face, Face.prototype.rotateX)
	var face2 = face.rotateX(PI / 2).translate(1, 2, -3)
	assert.doTest(face, face2, [StraightEdge.throughPoints(V3(1, 2), V3(10, 2))], [V3(10, 2)])

	var face3 = PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).rotateX(PI / 2).translate(1, 2, -3)
	assert.doTest(face, face3, [StraightEdge.throughPoints(V3(1, 2), V3(6, 2))], [])
});
QUnit.test( "face touching edge of test face with its middle", function( assert ) {
	var face = PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(10, 0), V3(10, 10), V3(0, 10)])
	var face2 = PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(12, 0), V3(12, 10), V3(0, 10)]).rotateX(PI/2).translate(-1, 0, -5)
	assert.doTest(face, face2, [], [V3(0, 0), V3(10, 0)]) // StraightEdge.throughPoints(V3(0, 0), V3(10, 0))
	assert.doTest(face, face2.flipped(), [], [])
});
QUnit.test( "testing V-shape with spine touching middle of test face", function( assert ) {
	var face = PlaneFace.forVertices(P3.XY, [V3(-10, -10), V3(10, -10), V3(10, 10), V3(-10, 10)])
	var brep = new B2([
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).flipped(),
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).rotateX(PI/4)
	])

	assert.doTest(face, brep.rotateX(PI/2), [], [])
	assert.doTest(face, brep.rotateX(PI/2).flipped(), [], [])
	assert.doTest(face, brep.rotateX(-PI/2), [], [])
	assert.doTest(face, brep.rotateX(-PI/2).flipped(), [StraightEdge.throughPoints(V3(0, 0), V3(5, 0)), StraightEdge.throughPoints(V3(5, 0), V3(0, 0))], [])
});
QUnit.test( "V-shape splitting test face at spine", function( assert ) {
	var face = PlaneFace.forVertices(P3.XY, [V3(-10, -10), V3(10, -10), V3(10, 10), V3(-10, 10)])
	var brep = new B2([
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).flipped().rotateX(-PI/4),
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).rotateX(PI/4)
	])

	assert.doTest(face, brep, [StraightEdge.throughPoints(V3(0, 0), V3(5, 0))], [])
	assert.doTest(face, brep.flipped(), [StraightEdge.throughPoints(V3(5, 0), V3(0, 0))], [])
});
QUnit.test( "V-shape with spine parallel to test face normal; touching test face edge", function( assert ) {
	var face = PlaneFace.forVertices(P3.XY, [V3(-20, 0), V3(20, 0), V3(20, 10), V3(-20, 10)])
	var brep = new B2([
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).flipped().rotateX(-PI/4),
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).rotateX(PI/4)
	]).rotateY(PI/2).translate(0, 0, 2)

	assert.doTest(face, brep, [StraightEdge.throughPoints(V3(0, 0), V3(5 * sqrt(2), 5 * sqrt(2))), StraightEdge.throughPoints(V3(-5 * sqrt(2), 5 * sqrt(2)), V3(0, 0))], [])
	assert.doTest(face, brep.flipped(), [
		StraightEdge.throughPoints(V3(5 * sqrt(2), 5 * sqrt(2)), V3(0, 0)), StraightEdge.throughPoints(V3(0, 0), V3(-5 * sqrt(2), 5 * sqrt(2)))], [])
	assert.doTest(face, brep.rotateZ(PI), [], [])
	assert.doTest(face, brep.rotateZ(PI).flipped(), [], [])
});
QUnit.test( "V-shape with spine parallel to test face normal; touching test face edge; splitting test face edge", function( assert ) {
	var face = PlaneFace.forVertices(P3.XY, [V3(-20, 0), V3(20, 0), V3(20, 10), V3(-20, 10)])
	var brep = new B2([
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).flipped().rotateX(-PI/4),
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).rotateX(PI/4)
	]).rotateY(PI/2).translate(0, 0, 2)

	assert.doTest(face, brep.rotateZ(-PI/2), [
		StraightEdge.throughPoints(V3(5 * sqrt(2), 5 * sqrt(2)), V3(0, 0))], [V3(0, 0)])
	assert.doTest(face, brep.rotateZ(-PI/2).flipped(), [
		StraightEdge.throughPoints(V3(0, 0), V3(5 * sqrt(2), 5 * sqrt(2)))], [V3(0, 0)])
});
QUnit.test( "V-shape with spine parallel to test face normal; one wing of V overlaps edge of test face", function( assert ) {
	var face = PlaneFace.forVertices(P3.XY, [V3(-10, -10), V3(10, -10), V3(10, 10), V3(-10, 10)])
	var brep = new B2([
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).flipped().rotateX(-PI/4),
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)])
	]).rotateY(PI/2).translate(10, 5, 2)

	assert.doTest(face, brep, [
		StraightEdge.throughPoints(V3(5, 10), V3(10, 5))], [V3(10, 5), V3(5, 10), V3(10, 10)])
	assert.doTest(face, brep.flipped(), [
		StraightEdge.throughPoints(V3(10, 5), V3(5, 10))], [V3(10, 5), V3(5, 10)])
});
QUnit.test( "testing V-shape with spine touching test face, overlapping edge of test face", function( assert ) {
	var face = PlaneFace.forVertices(P3.XY, [V3(-10, -10), V3(2, -10), V3(2, 10), V3(-10, 10)])
	var brep = new B2([
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).flipped(),
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).rotateX(PI/2)
	])

	assert.doTest(face, brep.rotateX(PI/4), [], [])
	assert.doTest(face, brep.rotateX(PI/4).flipped(), [], [])
	assert.doTest(face, brep.rotateX(-PI/4), [StraightEdge.throughPoints(V3(0, 0), V3(2, 0))], [V3(2, 0)])
	assert.doTest(face, brep.rotateX(-PI/4).flipped(), [StraightEdge.throughPoints(V3(2, 0), V3(0, 0))], [V3(2, 0)])
	assert.doTest(face, brep.rotateX(-PI/4 * 3), [], [])
	assert.doTest(face, brep.rotateX(-PI/4 * 3).flipped(), [StraightEdge.throughPoints(V3(2, 0), V3(0, 0)), StraightEdge.throughPoints(V3(0, 0), V3(2, 0))], [])
});
QUnit.test( "testing V-shape with spine touching test face, touching edge of test face", function( assert ) {
	var face = PlaneFace.forVertices(P3.XY, [V3(-10, -10), V3(5, -10), V3(5, 10), V3(-10, 10)])
	var brep = new B2([
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).flipped(),
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).rotateX(PI/2)
	])

	assert.doTest(face, brep.rotateX(PI/4), [], [])
	assert.doTest(face, brep.rotateX(PI/4).flipped(), [], [])
	assert.doTest(face, brep.rotateX(-PI/4), [StraightEdge.throughPoints(V3(0, 0), V3(5, 0))], [V3(5, 0)])
	assert.doTest(face, brep.rotateX(-PI/4).flipped(), [StraightEdge.throughPoints(V3(5, 0), V3(0, 0))], [V3(5, 0)])
	assert.doTest(face, brep.rotateX(-PI/4 * 3), [], [])
	assert.doTest(face, brep.rotateX(-PI/4 * 3).flipped(), [StraightEdge.throughPoints(V3(5, 0), V3(0, 0)), StraightEdge.throughPoints(V3(0, 0), V3(5, 0))], [])
});
QUnit.test( "testing V-shape with spine touching test face edge middle", function( assert ) {
	var face = PlaneFace.forVertices(P3.XY, [V3(-10, 0), V3(10, 0), V3(10, 10), V3(-10, 10)])
	var brep = new B2([
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).flipped(),
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).rotateX(PI/2)
	])

	assert.doTest(face, brep.rotateX(PI/4), [], [])
	assert.doTest(face, brep.rotateX(PI/4).flipped(), [], [V3(0, 0), V3(5, 0)]) // StraightEdge.throughPoints(V3(0, 0), V3(5, 0)), StraightEdge.throughPoints(V3(0, 0), V3(5, 0))
	assert.doTest(face, brep.rotateX(-PI/4), [], [V3(0, 0), V3(5, 0)]) // StraightEdge.throughPoints(V3(0, 0), V3(5, 0)), StraightEdge.throughPoints(V3(0, 0), V3(5, 0))
	assert.doTest(face, brep.rotateX(-PI/4).flipped(), [], [])
	assert.doTest(face, brep.rotateX(-PI/4 * 3), [], [])
	assert.doTest(face, brep.rotateX(-PI/4 * 3).flipped(), [], [V3(0, 0), V3(5, 0)]) // StraightEdge.throughPoints(V3(0, 0), V3(5, 0)), StraightEdge.throughPoints(V3(0, 0), V3(5, 0))
});
QUnit.test( "testing V-shape with spine overlapping test face edge", function( assert ) {
	var face = PlaneFace.forVertices(P3.XY, [V3(-10, 0), V3(3, 0), V3(3, 10), V3(-10, 10)])
	var brep = new B2([
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).flipped(),
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).rotateX(PI/2)
	])

	assert.doTest(face, brep.rotateX(PI/4), [], [])
	assert.doTest(face, brep.rotateX(PI/4).flipped(), [], [V3(0, 0), V3(3, 0)]) // StraightEdge.throughPoints(V3(0, 0), V3(3, 0)), StraightEdge.throughPoints(V3(0, 0), V3(3, 0))
	assert.doTest(face, brep.rotateX(-PI/4), [], [V3(0, 0), V3(3, 0)]) // StraightEdge.throughPoints(V3(0, 0), V3(3, 0)), StraightEdge.throughPoints(V3(0, 0), V3(3, 0))
	assert.doTest(face, brep.rotateX(-PI/4).flipped(), [], [])
	assert.doTest(face, brep.rotateX(-PI/4 * 3), [], [])
	assert.doTest(face, brep.rotateX(-PI/4 * 3).flipped(), [], [V3(0, 0), V3(3, 0)]) // StraightEdge.throughPoints(V3(0, 0), V3(3, 0)), StraightEdge.throughPoints(V3(0, 0), V3(3, 0))
});
QUnit.test( "V-shape with spine parallel to test face normal; touching test face corner", function( assert ) {
	var face = PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(10, 0), V3(10, 10), V3(0, 10)])
	// splitting contour in base position:
	var brep = new B2([
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).flipped().rotateX(-PI/8),
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).rotateX(PI/8)
	]).rotateY(PI/2).translate(0, 0, 2)

	assert.doTest(face, brep, [
		StraightEdge.throughPoints(V3(0, 0), V3.polar(10, PI/2-PI/8))], [V3(0, 0)])
});
QUnit.test( "V-shape with spine parallel to test face normal; touching test face corner", function( assert ) {
	var face = PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(10, 0), V3(10, 10), V3(0, 10)])
	// splitting contour in base position:
	var brep = new B2([
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).flipped().rotateX(-PI/8),
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).rotateX(PI/8)
	]).rotateY(PI/2).translate(0, 0, 2)

	assert.doTest(face, brep, [
		StraightEdge.throughPoints(V3(0, 0), V3.polar(10, PI/2-PI/8))], [V3(0, 0)])
});
QUnit.test( " ABC V-shape with spine parallel to test face normal; touching test face corner", function( assert ) {
	var box = B2.box(8, 9, 10, "box")
	var face = box.faces[1], testBrep
	// splitting contour in base position:

	testBrep = B2.tetrahedron(V3(-1, 1, 9), V3(5, 1, 9), V3(-1, -4, 14), V3(2, -4, 10))
	assert.doTestWithBrep(face, box, testBrep, [], [],
		"face of tetra touches edge of test face, tetra intersects test volume (test face not part of result, points dont really matter)")

	testBrep = B2.tetrahedron(V3(-1, 1, 9), V3(4, 1, 9), V3(-1, -4, 14), V3(-1, -4, 10)).flipped()
	assert.doTestWithBrep(face, box, testBrep, [StraightEdge.throughPoints(V3(0, 0, 10), V3(3, 0 , 10))], [V3(3, 0, 10), V3(0, 0, 10)],
		"face of _flipped_ tetra touches edge of testface and also intersects main volume (expect point on edge)")

	testBrep = B2.tetrahedron(V3(-1, 0, 10), V3(5, 0, 10), V3(2, 0, 14), V3(2, -4, 10))
	assert.doTestWithBrep(face, box, testBrep, [], [],
		"volumes touch edge-edge but no overlap (empty result volume; generated points dont matter)")

	testBrep = B2.tetrahedron(V3(-1, 0, 10), V3(5, 0, 10), V3(-1, -2, 8), V3(-1, 2, 8)).flipped()
	assert.doTestWithBrep(face, box, testBrep, [], [V3(5, 0, 10), V3(0, 0, 10)],
		"Tetrahedron is flipped, testface only touched at edge, needs point on edge as tetrahedron intersects side of volume")


	testBrep = B2.tetrahedron(V3(-1, 0, 10), V3(5, 0, 10), V3(-1, 0, 14), V3(-1, -4, 9)).flipped()
	assert.doTestWithBrep(face, box, testBrep, [], [],
		"volumes do not intersect, tetra is flipped, touches box edge-edge (result would be entire box, so no points)")
});
QUnit.test( " coplanar things", function( assert ) {
	var box = B2.box(8, 9, 10, "box")
	var face = box.faces[1], testBrep
	// splitting contour in base position:

	testBrep = B2.tetrahedron(V3(-1, -1, 10), V3(-1, 5, 10), V3(5, -1, 10), V3(-1, -1, 5)).flipped()
	assert.doTestWithBrep(face, box, testBrep, [StraightEdge.throughPoints(V3(0, 4, 10), V3(4, 0, 10))], [V3(0, 4, 10), V3(4, 0, 10)],
		"cut off corner of box with flipped tetra (anti)coplanar to test face")

	testBrep = B2.tetrahedron(V3(-1, -1, 10), V3(-1, 5, 10), V3(4, 0, 10), V3(-1, -1, 5)).flipped()
	assert.doTestWithBrep(face, box, testBrep, [StraightEdge.throughPoints(V3(0, 4, 10), V3(4, 0, 10))], [V3(0, 4, 10), V3(4, 0, 10)],
		"cut off corner of box with flipped tetra (anti)coplanar to test face")

	testBrep = B2.tetrahedron(V3(1, 0, 10), V3(1, 5, 10), V3(4, 0, 10), V3(-1, -1, 5)).flipped()
	assert.doTestWithBrep(face, box, testBrep,
		[StraightEdge.throughPoints(V3(1, 0, 10), V3(1, 5, 10)), StraightEdge.throughPoints(V3(1, 5, 10), V3(4, 0, 10))],
		[V3(1, 0, 10), V3(4, 0, 10)],
		"cut off corner of box with flipped tetra (anti)coplanar to test face")

	testBrep = B2.tetrahedron(V3(0, 0, 10), V3(0, 5, 10), V3(5, 0, 10), V3(0, 0, 5)).flipped()
	assert.doTestWithBrep(face, box, testBrep,
		[StraightEdge.throughPoints(V3(0, 5, 10), V3(5, 0, 10))],
		[V3(0, 5, 10), V3(5, 0, 10)],
		"cut off corner of box with flipped tetra (anti)coplanar to test face")

	testBrep = B2.tetrahedron(V3(0, 0, 10), V3(-1, 5, 10), V3(5, -1, 10), V3(-1, -1, 5)).flipped()
	assert.doTestWithBrep(face, box, testBrep, [StraightEdge.throughPoints(V3(0, 4, 10), V3(4, 0, 10))], [V3(0, 4, 10), V3(4, 0, 10)],
		"cut off corner of box with flipped tetra (anti)coplanar to test face")

	testBrep = B2.tetrahedron(V3(0, 0, 10), V3(1, 5, 10), V3(5, 1, 10), V3(0, 0, 5)).flipped()
	assert.doTestWithBrep(face, box, testBrep,
		[StraightEdge.throughPoints(V3(0, 0, 10), V3(1, 5, 10)), StraightEdge.throughPoints(V3(1, 5, 10), V3(5, 1, 10)), StraightEdge.throughPoints(V3(5, 1, 10), V3(0, 0, 10))],
		[],
		"cut hole at corner of test face")

});
QUnit.test( " coplanar things 2", function( assert ) {
	var box = B2.box(8, 9, 10, "box")
	var face = box.faces[1], testBrep
	// splitting contour in base position:

	testBrep = B2.tetrahedron(V3(-1, -1, 10), V3(-1, 5, 10), V3(5, -1, 10), V3(-1, -1, 5))
	assert.doTestWithBrep(face, box, testBrep, [StraightEdge.throughPoints(V3(0, 4, 10), V3(4, 0, 10))], [V3(0, 4, 10), V3(4, 0, 10)],
		"cut off corner of box with flipped tetra (anti)coplanar to test face")

	testBrep = B2.tetrahedron(V3(-1, -1, 10), V3(-1, 5, 10), V3(5, -1, 10), V3(-1, -1, 5)).flipped()
	assert.doTestWithBrep(face, box.flipped(), testBrep, [StraightEdge.throughPoints(V3(0, 4, 10), V3(4, 0, 10))], [V3(0, 4, 10), V3(4, 0, 10)],
		"cut off corner of box with flipped tetra (anti)coplanar to test face")

});




QUnit.test( "test assembly", function( assert ) {
	var baseFace = PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(10, 0), V3(10, 10), V3(0, 10)])

	var extrude = B2.extrudeVertices([V3(5, -1), V3(2, 2), V3(8, 2)], P3.XY.flipped(), V3(0, 0, 10)).translate(0, 0, -2).flipped()
	var result = PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(4, 0), V3(2, 2), V3(8, 2), V3(6, 0), V3(10, 0), V3(10, 10), V3(0, 10)])
	assert.doTest2(baseFace, extrude, [result], "volume cuts edge of test face (twice)")

	var extrude2 = B2.extrudeVertices([V3(5, 0), V3(2, 3), V3(8, 3)], P3.XY.flipped(), V3(0, 0, 10)).translate(0, 0, -2).flipped()
	var result2 = PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(10, 0), V3(10, 10), V3(0, 10)], [V3(5, 0), V3(2, 3), V3(8, 3)])
	assert.doTest2(baseFace, extrude2, [result2], "volume touches inside of test face edge")


	// from test case 3:4
	// V-shape spine touching middle, splits volume enclosing in both directions
	var brep = new B2([
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).flipped(),
		PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 0), V3(5, 10), V3(0, 10)]).rotateX(PI/4)
	]).rotateX(-PI/2).flipped().translate(0, 2, 0)
	// degenerate cycle in middle of face
	var result4 = PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(10, 0), V3(10, 10), V3(0, 10)], [V3(5, 2), V3(0, 2)])


	var extrude5 = B2.extrudeVertices([V3(0, 0), V3(3, 2), V3(2, 3)], P3.XY.flipped(), V3(0, 0, 10)).translate(0, 0, -2).flipped()
	var result5 = PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(2, 3), V3(3, 2), V3(0, 0), V3(10, 0), V3(10, 10), V3(0, 10)])
	assert.doTest2(baseFace, extrude5, [result5], "volume touches inside of test face corner")

	var extrude6 = B2.extrudeVertices([V3(1, 0), V3(3, 0), V3(8, 10)], P3.XY.flipped(), V3(0, 0, 10)).translate(0, 0, -2).flipped()
	var result6 = [PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(1, 0), V3(8, 10), V3(3, 0), V3(10, 0), V3(10, 10), V3(0, 10)])]
	assert.doTest2(baseFace, extrude6, result6, "volume touches inside of test face corner")

	var extrude7 = B2.extrudeVertices([V3(1, 0), V3(3, 0), V3(8, 10), V3(7, 10)], P3.XY.flipped(), V3(0, 0, 10)).translate(0, 0, -2).flipped()
	var result7 = [PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(1, 0), V3(7, 10), V3(0, 10)]), PlaneFace.forVertices(P3.XY, [V3(10, 10), V3(8, 10), V3(3, 0), V3(10, 0)])]
	assert.doTest2(baseFace, extrude7, result7, "volume touches inside of test face corner")


});

QUnit.test( "test assembly holes", function( assert ) {
	var baseFace = PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(10, 0), V3(10, 10), V3(0, 10)])

	var extrude7 = B2.extrudeVertices([V3(0, 0), V3(5, 8), V3(8, 5)], P3.XY.flipped(), V3(0, 0, 10)).translate(0, 0, -2).flipped()
	var result7 = [PlaneFace.forVertices(P3.XY, [V3(0, 0), V3(5, 8), V3(8, 5), V3(0, 0), V3(10, 0), V3(10, 10), V3(0, 10)])]
	assert.doTest2(baseFace, extrude7, result7, "volume touches inside of test face corner")

	var extrude8 = B2.extrudeVertices([V3(1, 1), V3(1, -1), V3(-1, 1)], P3.XY.flipped(), V3(0, 0, 10)).translate(0, 0, -2).flipped()
	var result8 = [PlaneFace.forVertices(P3.XY, [V3(1, 0), V3(10, 0), V3(10, 10), V3(0, 10), V3(0, 1), V3(1, 1)])]
	assert.doTest2(baseFace, extrude8, result8, "volume touches inside of test face corner")

	var extrude9 = B2.extrudeVertices([V3(-1, -1), V3(1, 1), V3(1, -1)], P3.XY.flipped(), V3(0, 0, 10)).translate(0, 0, -2).flipped()
	var result9 = [PlaneFace.forVertices(P3.XY, [V3(1, 0), V3(10, 0), V3(10, 10), V3(0, 10), V3(0, 0), V3(1, 1)])]
	assert.doTest2(baseFace, extrude9, result9, "volume touches inside of test face corner")

	var extrude10 = B2.extrudeVertices([V3(1, 1), V3(2, 2), V3(2, 1)], P3.XY.flipped(), V3(0, 0, 10)).translate(0, 0, -2).flipped()
	var result10 = [PlaneFace.forVertices(P3.XY, [V3(10, 0), V3(10, 10), V3(0, 10), V3(0, 0)], [V3(1, 1), V3(2, 2), V3(2, 1)])]
	assert.doTest2(baseFace, extrude10, result10, "adding a hole")

	var baseFace11 = PlaneFace.forVertices(P3.XY, [V3(10, 0), V3(10, 10), V3(0, 10), V3(0, 0)], [V3(1, 1), V3(2, 2), V3(2, 1)])
	var extrude11 = B2.extrudeVertices([V3(5, 5), V3(6, 6), V3(6, 5)], P3.XY.flipped(), V3(0, 0, 10)).translate(0, 0, -2).flipped()
	var result11 = [PlaneFace.forVertices(P3.XY, [V3(10, 0), V3(10, 10), V3(0, 10), V3(0, 0)], [V3(1, 1), V3(2, 2), V3(2, 1)], [V3(5, 5), V3(6, 6), V3(6, 5)])]
	assert.doTest2(baseFace11, extrude11, result11, "adding a hole to a face which already has one")


	let baseFace12 = PlaneFace.forVertices(P3.XY, [V3(10, 0), V3(10, 10), V3(0, 10), V3(0, 0)], [V3(1, 1), V3(5, 5), V3(5, 1)])
	var extrude12 = B2.extrudeVertices([V3(2, 1.5), V3(2, 4), V3(4.5, 4)], P3.XY.flipped(), V3(0, 0, 10)).translate(0, 0, -2).flipped()
	var result12 = [PlaneFace.forVertices(P3.XY, [V3(10, 0), V3(10, 10), V3(0, 10), V3(0, 0)], [V3(5, 5), V3(5, 1), V3(1, 1), V3(2, 2), V3(2, 4), V3(4, 4)])]
	assert.doTest2(baseFace12, extrude12, result12, "extending an existing hole")

	var baseFace13 = PlaneFace.forVertices(P3.XY, [V3(10, 0), V3(10, 10), V3(0, 10), V3(0, 0)], [V3(1, 1), V3(5, 5), V3(5, 1)])
	var extrude13 = B2.extrudeVertices([V3(3, -1), V3(4, -1), V3(4, 2), V3(3, 2)], P3.XY.flipped(), V3(0, 0, 10)).translate(0, 0, -2).flipped()
	var result13 = [PlaneFace.forVertices(P3.XY, [V3(10, 0), V3(10, 10), V3(0, 10), V3(0, 0), V3(3, 0), V3(3, 1), V3(1, 1), V3(5, 5), V3(5, 1), V3(4, 1), V3(4, 0)])]
	assert.doTest2(baseFace13, extrude13, result13, "removing a hole by cutting a channel")

	var baseFace14 = PlaneFace.forVertices(P3.XY, [V3(10, 0), V3(10, 10), V3(0, 10), V3(0, 0)], [V3(1, 1), V3(5, 5), V3(5, 1)])
	var extrude14 = B2.extrudeVertices([V3(1, 1), V3(1, 5), V3(5, 5)], P3.XY.flipped(), V3(0, 0, 10)).translate(0, 0, -2).flipped()
	var result14 = [PlaneFace.forVertices(P3.XY, [V3(10, 0), V3(10, 10), V3(0, 10), V3(0, 0)], [V3(5, 5), V3(5, 1), V3(1, 1), V3(1, 5)])]
	assert.doTest2(baseFace14, extrude14, result14, "extending an existing hole")
});
QUnit.test( "B2.prototype.minus remove half of a half-pie", function( assert ) {
	let pie = B2.puckman(8, 180*DEG, 5, 'pie/2')
	let boxKnife = B2.box(11, 10, 7, 'knife').translate(-10, -1, -1)

	let resultTopPoint = V3(1,8*Math.sin(Math.acos(1/8)), 0)
	let result = B2.extrudeEdges([
			StraightEdge.throughPoints(V3(8, 0, 0), V3(1, 0, 0)),
			StraightEdge.throughPoints(V3(1, 0, 0), resultTopPoint),
			PCurveEdge.forCurveAndTs(EllipseCurve.circle(8), Math.acos(1/8), 0)],
		P3.XY.flipped(), V3(0,0,5), 'pie/4')
	assert.b2Equal(pie, boxKnife, pie.minus(boxKnife), result)

	assert.b2Equal(pie, boxKnife, pie.minus(boxKnife.translate(-1, 0, 0)), B2.puckman(8, 90*DEG, 5, 'pie/4'))
});
QUnit.test( "B2.prototype.minus cut hole through side of pie", function( assert ) {
	let pie = B2.puckman(8, 180*DEG, 5, 'pie/2')
	let punch = B2.box(5, 10, 3, 'knife').translate(1, -1, 1)
	let result = new B2([
		new PlaneFace(new PlaneSurface(P3(V3(0, -1, 0), 0)), [
			StraightEdge.throughPoints(V3(0, 0, 0), V3(8, 0, 0)),
			StraightEdge.throughPoints(V3(8, 0, 0), V3(8, 0, 5)),
			StraightEdge.throughPoints(V3(8, 0, 5), V3(0, 0, 5)),
			StraightEdge.throughPoints(V3(0, 0, 5), V3(0, 0, 0))], [[
			StraightEdge.throughPoints(V3(1, 0, 1), V3(1, 0, 4)),
			StraightEdge.throughPoints(V3(1, 0, 4), V3(6, 0, 4)),
			StraightEdge.throughPoints(V3(6, 0, 4), V3(6, 0, 1)),
			StraightEdge.throughPoints(V3(6, 0, 1), V3(1, 0, 1))]]),
		new RotationFace(new CylinderSurface(new EllipseCurve(V3(0, 0, 0), V3(8, 0, 0), V3(0, -8, 0)), V3(0, 0, -1)), [
			new PCurveEdge(new EllipseCurve(V3(0, 0, 0), V3(8, 0, 0), V3(0, -8, 0)), V3(8, 0, 0), V3(-8, 9.797174393178826e-16, 0), 0, -3.141592653589793, null, V3(0, 8, 0), V3(-9.797174393178826e-16, -8, 0)),
			StraightEdge.throughPoints(V3(-8, 9.797174393178826e-16, 0), V3(-8, 9.797174393178826e-16, 5)),
			new PCurveEdge(new EllipseCurve(V3(0, 0, 5), V3(8, 0, 0), V3(0, -8, 0)), V3(-8, 9.797174393178826e-16, 5), V3(8, 0, 5), -3.141592653589793, 0, null, V3(9.797174393178826e-16, 8, 0), V3(0, -8, 0)),
			StraightEdge.throughPoints(V3(8, 0, 5), V3(8, 0, 0))], [[
			StraightEdge.throughPoints(V3(1, 7.937253933193773, 4), V3(1, 7.937253933193773, 1)),
			new PCurveEdge(new EllipseCurve(V3(0, 0, 1), V3(8, 0, 0), V3(0, -8, 0)), V3(1, 7.937253933193773, 1), V3(6, 5.291502622129181, 1), -1.4454684956268313, -0.7227342478134156, null, V3(7.937253933193772, -0.9999999999999991, 0), V3(5.2915026221291805, -6, 0)),
			StraightEdge.throughPoints(V3(6, 5.291502622129181, 1), V3(6, 5.291502622129181, 4)),
			new PCurveEdge(new EllipseCurve(V3(0, 0, 4), V3(8, 0, 0), V3(0, 8, 0)), V3(6, 5.291502622129181, 4), V3(1, 7.937253933193773, 4), 0.7227342478134156, 1.4454684956268313, null, V3(-5.2915026221291805, 6, 0), V3(-7.937253933193772, 0.9999999999999991, 0))]]),
		new PlaneFace(new PlaneSurface(P3(V3(-1.2246467991473532e-16, -1, 0), 0)), [
			StraightEdge.throughPoints(V3(-8, 9.797174393178826e-16, 0), V3(0, 0, 0)),
			StraightEdge.throughPoints(V3(0, 0, 0), V3(0, 0, 5)),
			StraightEdge.throughPoints(V3(0, 0, 5), V3(-8, 9.797174393178826e-16, 5)),
			StraightEdge.throughPoints(V3(-8, 9.797174393178826e-16, 5), V3(-8, 9.797174393178826e-16, 0))], []),
		new PlaneFace(new PlaneSurface(P3(V3(0, 0, -1), 0)), [
			StraightEdge.throughPoints(V3(8, 0, 0), V3(0, 0, 0)),
			StraightEdge.throughPoints(V3(0, 0, 0), V3(-8, 9.797174393178826e-16, 0)),
			new PCurveEdge(new EllipseCurve(V3(0, 0, 0), V3(8, 0, 0), V3(0, -8, 0)), V3(-8, 9.797174393178826e-16, 0), V3(8, 0, 0), -3.141592653589793, 0, null, V3(9.797174393178826e-16, 8, 0), V3(0, -8, 0))], []),
		new PlaneFace(new PlaneSurface(P3(V3(0, 0, 1), 5)), [
			new PCurveEdge(new EllipseCurve(V3(0, 0, 5), V3(8, 0, 0), V3(0, -8, 0)), V3(8, 0, 5), V3(-8, 9.797174393178826e-16, 5), 0, -3.141592653589793, null, V3(0, 8, 0), V3(-9.797174393178826e-16, -8, 0)),
			StraightEdge.throughPoints(V3(-8, 9.797174393178826e-16, 5), V3(0, 0, 5)),
			StraightEdge.throughPoints(V3(0, 0, 5), V3(8, 0, 5))], []),
		new PlaneFace(new PlaneSurface(P3(V3(1, 0, 0), 1)), [
			StraightEdge.throughPoints(V3(1, 0, 4), V3(1, 0, 1)),
			StraightEdge.throughPoints(V3(1, 0, 1), V3(1, 7.937253933193773, 1)),
			StraightEdge.throughPoints(V3(1, 7.937253933193773, 1), V3(1, 7.937253933193773, 4)),
			StraightEdge.throughPoints(V3(1, 7.937253933193773, 4), V3(1, 0, 4))], []),
		new PlaneFace(new PlaneSurface(P3(V3(-1, 0, 0), -6)), [
			StraightEdge.throughPoints(V3(6, 0, 1), V3(6, 0, 4)),
			StraightEdge.throughPoints(V3(6, 0, 4), V3(6, 5.291502622129181, 4)),
			StraightEdge.throughPoints(V3(6, 5.291502622129181, 4), V3(6, 5.291502622129181, 1)),
			StraightEdge.throughPoints(V3(6, 5.291502622129181, 1), V3(6, 0, 1))], []),
		new PlaneFace(new PlaneSurface(P3(V3(0, 0, 1), 1)), [
			StraightEdge.throughPoints(V3(1, 0, 1), V3(6, 0, 1)),
			StraightEdge.throughPoints(V3(6, 0, 1), V3(6, 5.291502622129181, 1)),
			new PCurveEdge(new EllipseCurve(V3(0, 0, 1), V3(8, 0, 0), V3(0, -8, 0)), V3(6, 5.291502622129181, 1), V3(1, 7.937253933193773, 1), -0.7227342478134156, -1.4454684956268313, null, V3(-5.2915026221291805, 6, 0), V3(-7.937253933193772, 0.9999999999999991, 0)),
			StraightEdge.throughPoints(V3(1, 7.937253933193773, 1), V3(1, 0, 1))], []),
		new PlaneFace(new PlaneSurface(P3(V3(0, 0, -1), -4)), [
			StraightEdge.throughPoints(V3(6, 0, 4), V3(1, 0, 4)),
			StraightEdge.throughPoints(V3(1, 0, 4), V3(1, 7.937253933193773, 4)),
			new PCurveEdge(new EllipseCurve(V3(0, 0, 4), V3(8, 0, 0), V3(0, 8, 0)), V3(1, 7.937253933193773, 4), V3(6, 5.291502622129181, 4), 1.4454684956268313, 0.7227342478134156, null, V3(7.937253933193772, -0.9999999999999991, 0), V3(5.2915026221291805, -6, 0)),
			StraightEdge.throughPoints(V3(6, 5.291502622129181, 4), V3(6, 0, 4))], [])])

	assert.b2Equal(pie, punch, pie.minus(punch), result)
	console.log(pie.minus(punch).sce)
});

QUnit.test( "B2.prototype.minus including ProjectedCurveSurface", function( assert ) {
	let s1 = new ProjectedCurveSurface(new BezierCurve(V3(0, 0, 0), V3(-5, 5, 0), V3(15, 5, 0), V3(10, 0, 0), -0.1, 1.1), V3(0, 0, -1), 0, 1)
	let s2 = new ProjectedCurveSurface(new BezierCurve(V3(0, 0, 0), V3(-5, 5, 0), V3(15, 5, 0), V3(10, 0, 0), -0.1, 1.1), V3(0, 0, -1), 0, 1)
	let a = B2.extrudeEdges([
			StraightEdge.throughPoints(V3(10, 0, 0), V3(0, 0, 0)),
			PCurveEdge.forCurveAndTs(new BezierCurve(V3.ZERO, V3(-5, 5, 0), V3(15, 5, 0), V3(10, 0, 0)), 0, 1)],
		P3.XY.flipped(), V3(0,0,5), 'a/4')
	let punch = B2.box(5, 10, 3, 'knife').translate(1, -1, 1)

	let result = new B2([
		new PlaneFace(new PlaneSurface(P3(V3(0, -1, 0), 0)), [
			StraightEdge.throughPoints(V3(0, 0, 0), V3(10, 0, 0)),
			StraightEdge.throughPoints(V3(10, 0, 0), V3(10, 0, 5)),
			StraightEdge.throughPoints(V3(10, 0, 5), V3(0, 0, 5)),
			StraightEdge.throughPoints(V3(0, 0, 5), V3(0, 0, 0))], [[
			StraightEdge.throughPoints(V3(1, 0, 1), V3(1, 0, 4)),
			StraightEdge.throughPoints(V3(1, 0, 4), V3(6, 0, 4)),
			StraightEdge.throughPoints(V3(6, 0, 4), V3(6, 0, 1)),
			StraightEdge.throughPoints(V3(6, 0, 1), V3(1, 0, 1))]]),
		new RotationFace(new ProjectedCurveSurface(new BezierCurve(V3(0, 0, 0), V3(-5, 5, 0), V3(15, 5, 0), V3(10, 0, 0), -0.1, 1.1), V3(0, 0, -1), 0, 1), [
			new PCurveEdge(new BezierCurve(V3(0, 0, 0), V3(-5, 5, 0), V3(15, 5, 0), V3(10, 0, 0), -0.1, 1.1), V3(10, 0, 0), V3(0, 0, 0), 1, 0, null, V3(15, 15, 0), V3(15, -15, 0)),
			StraightEdge.throughPoints(V3(0, 0, 0), V3(0, 0, 5)),
			new PCurveEdge(new BezierCurve(V3(0, 0, 5), V3(-5, 5, 5), V3(15, 5, 5), V3(10, 0, 5), -0.1, 1.1), V3(0, 0, 5), V3(10, 0, 5), 0, 1, null, V3(-15, 15, 0), V3(-15, -15, 0)),
			StraightEdge.throughPoints(V3(10, 0, 5), V3(10, 0, 0))], [[
			StraightEdge.throughPoints(V3(1, 3.185436104380006, 4), V3(1, 3.185436104380006, 1)),
			new PCurveEdge(new BezierCurve(V3(0, 0, 1), V3(-5, 5, 1), V3(15, 5, 1), V3(10, 0, 1), -0.1, 1.1), V3(1, 3.185436104380006, 1), V3(6, 3.720106174228432, 1), 0.3059958942668147, 0.5446421518086273, null, V3(16.85436104380006, 5.820123171995558, 0), V3(22.201061742284324, -1.3392645542588184, 0)),
			StraightEdge.throughPoints(V3(6, 3.720106174228432, 1), V3(6, 3.720106174228432, 4)),
			new PCurveEdge(new BezierCurve(V3(10, 0, 4), V3(15, 5, 4), V3(-5, 5, 4), V3(0, 0, 4), -1.1, 0.1), V3(6, 3.720106174228432, 4), V3(1, 3.185436104380006, 4), 0.45535784819137265, 0.6940041057331853, null, V3(-22.201061742284324, 1.3392645542588197, 0), V3(-16.85436104380006, -5.820123171995558, 0))]]),
		new PlaneFace(new PlaneSurface(P3(V3(0, 0, -1), 0)), [
			StraightEdge.throughPoints(V3(10, 0, 0), V3(0, 0, 0)),
			new PCurveEdge(new BezierCurve(V3(0, 0, 0), V3(-5, 5, 0), V3(15, 5, 0), V3(10, 0, 0), -0.1, 1.1), V3(0, 0, 0), V3(10, 0, 0), 0, 1, null, V3(-15, 15, 0), V3(-15, -15, 0))], []),
		new PlaneFace(new PlaneSurface(P3(V3(0, 0, 1), 5)), [
			new PCurveEdge(new BezierCurve(V3(0, 0, 5), V3(-5, 5, 5), V3(15, 5, 5), V3(10, 0, 5), -0.1, 1.1), V3(10, 0, 5), V3(0, 0, 5), 1, 0, null, V3(15, 15, 0), V3(15, -15, 0)),
			StraightEdge.throughPoints(V3(0, 0, 5), V3(10, 0, 5))], []),
		new PlaneFace(new PlaneSurface(P3(V3(1, 0, 0), 1)), [
			StraightEdge.throughPoints(V3(1, 0, 4), V3(1, 0, 1)),
			StraightEdge.throughPoints(V3(1, 0, 1), V3(1, 3.185436104380006, 1)),
			StraightEdge.throughPoints(V3(1, 3.185436104380006, 1), V3(1, 3.185436104380006, 4)),
			StraightEdge.throughPoints(V3(1, 3.185436104380006, 4), V3(1, 0, 4))], []),
		new PlaneFace(new PlaneSurface(P3(V3(-1, 0, 0), -6)), [
			StraightEdge.throughPoints(V3(6, 0, 1), V3(6, 0, 4)),
			StraightEdge.throughPoints(V3(6, 0, 4), V3(6, 3.720106174228432, 4)),
			StraightEdge.throughPoints(V3(6, 3.720106174228432, 4), V3(6, 3.720106174228432, 1)),
			StraightEdge.throughPoints(V3(6, 3.720106174228432, 1), V3(6, 0, 1))], []),
		new PlaneFace(new PlaneSurface(P3(V3(0, 0, 1), 1)), [
			StraightEdge.throughPoints(V3(1, 0, 1), V3(6, 0, 1)),
			StraightEdge.throughPoints(V3(6, 0, 1), V3(6, 3.720106174228432, 1)),
			new PCurveEdge(new BezierCurve(V3(0, 0, 1), V3(-5, 5, 1), V3(15, 5, 1), V3(10, 0, 1), -0.1, 1.1), V3(6, 3.720106174228432, 1), V3(1, 3.185436104380006, 1), 0.5446421518086273, 0.3059958942668147, null, V3(-22.201061742284324, 1.3392645542588184, 0), V3(-16.85436104380006, -5.820123171995558, 0)),
			StraightEdge.throughPoints(V3(1, 3.185436104380006, 1), V3(1, 0, 1))], []),
		new PlaneFace(new PlaneSurface(P3(V3(0, 0, -1), -4)), [
			StraightEdge.throughPoints(V3(6, 0, 4), V3(1, 0, 4)),
			StraightEdge.throughPoints(V3(1, 0, 4), V3(1, 3.185436104380006, 4)),
			new PCurveEdge(new BezierCurve(V3(10, 0, 4), V3(15, 5, 4), V3(-5, 5, 4), V3(0, 0, 4), -1.1, 0.1), V3(1, 3.185436104380006, 4), V3(6, 3.720106174228432, 4), 0.6940041057331853, 0.45535784819137265, null, V3(16.85436104380006, 5.820123171995558, 0), V3(22.201061742284324, -1.3392645542588197, 0)),
			StraightEdge.throughPoints(V3(6, 3.720106174228432, 4), V3(6, 0, 4))], [])])

	assert.b2Equal(a, punch, a.minus(punch), result)
});