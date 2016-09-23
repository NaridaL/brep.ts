QUnit.module('brep2')


QUnit.test( "Face.equals", function( assert ) {
	var a = PlaneFace.forVertices(P3(V3(0, 0, 1), 0), [V3(0, 0, 0), V3(10, 0, 0), V3(10, 10, 0), V3(0, 10, 0)])
	var b = PlaneFace.forVertices(P3(V3(0, 0, 1), 0), [V3(0, 10, 0), V3(0, 0, 0), V3(10, 0, 0), V3(10, 10, 0)])
	var c = PlaneFace.forVertices(P3(V3(0, 0, -1), 0), [V3(0, 0, 0), V3(10, 0, 0), V3(10, 10, 0), V3(0, 10, 0)].slice().reverse())
	assert.ok(a.equals(a))

	assert.ok(a.equals(b))
	assert.ok(b.equals(a))

	assert.notOk(a.equals(c))
	assert.notOk(c.equals(a))

	assert.notOk(b.equals(c))
	assert.notOk(c.equals(b))
});
QUnit.test( "B2.equals", function( assert ) {
	var a = B2.tetrahedron(V3(5, 5, 5), V3(5, 5, -5), V3(10, 12, 1), V3(0, 12, 1))
	var b = B2.tetrahedron(V3(5, 5, 5), V3(10, 12, 1), V3(5, 5, -5), V3(0, 12, 1))
	var c = B2.tetrahedron(V3(5, 5, 5), V3(12, 12, 1), V3(5, 5, -5), V3(0, 12, 1))

	assert.ok(a.equals(a))

	assert.ok(a.equals(b))
	assert.ok(b.equals(a))

	assert.notOk(a.equals(c))
	assert.notOk(c.equals(a))

	assert.notOk(b.equals(c))
	assert.notOk(c.equals(b))
});
QUnit.test( "B2.tetrahedron", function( assert ) {
	var result = new B2([
		PlaneFace.forVertices(P3(V3(0.8320502943378437, -0.5547001962252291, 0), 1.386750490563073), [V3(5, 5, 5), V3(5, 5, -5), V3(9, 11, 1)]),
		PlaneFace.forVertices(P3(V3(-0.7071067811865475, -0.7071067811865475, 0), -7.071067811865475), [V3(5, 5, 5), V3(1, 9, 0), V3(5, 5, -5)]),
		PlaneFace.forVertices(P3(V3(-0.1003911722115382, 0.7362019295512802, -0.6692744814102547), 6.525426193749983), [V3(5, 5, -5), V3(1, 9, 0), V3(9, 11, 1)]),
		PlaneFace.forVertices(P3(V3(-0.25177250044296223, 0.6474150011390457, 0.7193500012656064), 5.57496250980845), [V3(9, 11, 1), V3(1, 9, 0), V3(5, 5, 5)])])
	var a = B2.tetrahedron(V3(5, 5, 5), V3(5, 5, -5), V3(9, 11, 1), V3(1, 9, 0))
	assert.B2equals(a, result)
});


QUnit.test( "Face.prototype.containsPoint", function( assert ) {
	var a = PlaneFace.forVertices(P3(V3(0, 0, 1), 0), [V3(0, 0, 0), V3(10, 0, 0), V3(10, 10, 0), V3(0, 10, 0)])
	assert.ok(a.containsPoint(V3(5, 5, 0)))
	assert.notOk(a.containsPoint(V3(11, 5, 0)))


	var b = PlaneFace.forVertices(P3(V3(0, 0, 1), 0), [V3(0, 10, 0), V3(0, 0, 0), V3(5, 0, 0), V3(6, 5, 0)])
	assert.ok(b.containsPoint(V3(2, 5, 0)))

	var c = PlaneFace.forVertices(P3(V3(0, 0, 1), 0), [V3(0, 10, 0), V3(0, 0, 0), V3(5, 0, 0), V3(6, 5, 0), V3(10, 5, 0)])
	assert.ok(c.containsPoint(V3(2, 5, 0)))



	a = a.rotateZ(deg2rad(30))
	var m = M4.rotationZ(deg2rad(30))
	console.log(a.toString())
	assert.ok(a.containsPoint(m.transformPoint(V3(5, 5, 0))))
	assert.notOk(a.containsPoint(m.transformPoint(V3(-5, 5, 0))))

	b = b.rotateZ(deg2rad(30))
	assert.ok(b.containsPoint(m.transformPoint(V3(2, 5, 0))))

	c = c.rotateZ(deg2rad(30))
	assert.ok(c.containsPoint(m.transformPoint(V3(2, 5, 0))))
});

QUnit.test( "Face.prototype.containsPoint", function( assert ) {
	var a = PlaneFace.forVertices(P3(V3(0, 0, 1), 0), [V3(0, 0, 0), V3(10, 0, 0), V3(10, 10, 0), V3(0, 10, 0)])
	assert.notOk(a.containsPoint(V3(-0.00000001, 11, 0)))

	assert.ok(PlaneFace.forVertices(P3(V3(0, 0, -1), 0), [V3(-1, -10, 0), V3(0, 25, 0), V3(25, 0, 0)]).containsPoint(V3(0, 0, 0)))
});
QUnit.test( "Face.withHole", function( assert ) {
	var a = PlaneFace.forVertices(P3(V3(0, 0, 1), 0), [V3(0, 0, 0), V3(10, 0, 0), V3(10, 10, 0), V3(0, 10, 0)])
	var holeVertices = [V3(2, 3, 0),V3(8, 7, 0),V3(7, 2, 0)]


	assert.notOk(a.containsPoint(V3(-0.00000001, 11, 0)))

});
QUnit.test( "segmentTouchOrIntersect", function( assert ) {
	assert.ok(segmentsTouchOrIntersect(V3.ZERO, V3(3, 2.2), V3(2, 2.8), V3(2.8, 2.0)))

});
QUnit.test( "NLA.eqAngle", function( assert ) {
	assert.ok(NLA.zeroAngle(0))
	assert.ok(NLA.zeroAngle(- NLA_PRECISION / 2))
	assert.ok(NLA.zeroAngle(2 * Math.PI - NLA_PRECISION / 2))
	assert.ok(NLA.zeroAngle(2 * Math.PI + NLA_PRECISION / 2))
	assert.ok(NLA.eqAngle(-Math.PI, Math.PI))
	assert.ok(NLA.eqAngle(0, 2 * Math.PI))
	assert.ok(NLA.eqAngle(0, 2 * Math.PI - NLA_PRECISION / 2))
	assert.ok(NLA.eqAngle(0, 2 * Math.PI + NLA_PRECISION / 2))
	assert.notOk(NLA.eqAngle(-Math.PI, 2 * Math.PI))
	assert.notOk(NLA.eqAngle(0, Math.PI))

});

QUnit.test( "angleRelativeNormal", function( assert ) {
	assert.fuzzyEquals(V3.X.angleRelativeNormal(V3.Y, V3.Z), Math.PI / 2 )
	assert.fuzzyEquals(V3.Y.angleRelativeNormal(V3(32, Math.sqrt(2), -Math.sqrt(2)), V3.X), -Math.PI / 4 )
	assert.fuzzyEquals(V3(-0.1, 1, 0).angleRelativeNormal(V3(0.0, 0, -1), V3.X), -Math.PI / 2)
	assert.fuzzyEquals(V3(-0, 1, 0).angleRelativeNormal(V3(2, 0, -1), V3.X), -Math.PI / 2)
});

QUnit.test( "splitsVolumeEnclosingFaces", function( assert ) {
	var brep = B2.tetrahedron(V3(0, 0,0),V3(10,0,0),V3(0,10,0),V3(0,0,10))
	// pointing into tetrahedon
	var edge = StraightEdge.throughPoints
	assert.ok(INSIDE == splitsVolumeEnclosingFaces(brep, edge(V3(0, 0,0),V3(10,0,0)), V3(0,1,1), V3(0,-1,1)))
	assert.ok(INSIDE == splitsVolumeEnclosingFaces(brep, edge(V3(0, 0,0),V3(10,0,0)), V3(0,1,1), V3(0,1,-1)))
	// pointing out of tetrahedon
	assert.notOk(INSIDE == splitsVolumeEnclosingFaces(brep, edge(V3(0, 0,0),V3(10,0,0)), V3(0,-1,0), V3(0,1,1)))

	assert.notOk(INSIDE == splitsVolumeEnclosingFaces(brep, edge(V3(0, 0,0),V3(10,0,0)), V3(0,-1,-1), V3(0,-1,1)))
	assert.notOk(INSIDE == splitsVolumeEnclosingFaces(brep, edge(V3(0, 0,0),V3(10,0,0)), V3(0,-1,-1), V3(0,1,-1)))
});


QUnit.test( "splitsVolumeEnclosingFaces 2", function( assert ) {
	var brep = B2.tetrahedron(V3(0, 0,0),V3(10,0,0),V3(0,10,0),V3(0,0,10))
	// pointing out of tetrahedon
	assert.equal(splitsVolumeEnclosingFaces(brep, StraightEdge.throughPoints(V3(0, 0,0),V3(10,0,0)), V3(0,1,0), V3(0,0,1)), COPLANAR_OPPOSITE)
	assert.equal(splitsVolumeEnclosingFaces(brep, StraightEdge.throughPoints(V3(0, 0,0),V3(10,0,0)), V3(0,1,0), V3(0,0,-1)), COPLANAR_SAME)
	assert.equal(splitsVolumeEnclosingFaces(brep, StraightEdge.throughPoints(V3(0, 0,0),V3(10,0,0)), V3(0,0,1), V3(0,1,0)), COPLANAR_OPPOSITE)
});
QUnit.test( "splitsVolumeEnclosingFaces 3", function( assert ) {
	var brep = B2.box(5, 5, 5).flipped()
	console.log(brep.sce)

	assert.equal(splitsVolumeEnclosingFaces(brep, StraightEdge.throughPoints(V3(0, 5,0),V3(0,0,0)), V3(0,0,-1), V3(1,0,0)), INSIDE)
	assert.equal(splitsVolumeEnclosingFaces(brep, StraightEdge.throughPoints(V3(0,0,0),V3(0, 5,0)), V3(0,0,-1), V3(1,0,0)), INSIDE)

	assert.equal(splitsVolumeEnclosingFaces(brep, StraightEdge.throughPoints(V3(0,0,0),V3(0, 5,0)), V3(0,0,1), V3(-1,0,0)), COPLANAR_OPPOSITE)
	assert.equal(splitsVolumeEnclosingFaces(brep, StraightEdge.throughPoints(V3(0,0,0),V3(0, 5,0)), V3(0,0,1), V3(1,0,0)), COPLANAR_SAME)
});
QUnit.test( "splitsVolumeEnclosingCone", function( assert ) {
	var brep = B2.box(5, 5, 5)

	assert.equal(splitsVolumeEnclosingCone(brep, V3.ZERO, V3(1, 1, 1)), INSIDE)
	assert.equal(splitsVolumeEnclosingCone(brep, V3.ZERO, V3(-1, 1, 1)), OUTSIDE)
	assert.equal(splitsVolumeEnclosingCone(brep, V3.ZERO, V3.X), ALONG_EDGE_OR_PLANE)
	assert.equal(splitsVolumeEnclosingCone(brep, V3.ZERO, V3(1, 1, 0)), ALONG_EDGE_OR_PLANE)
});

QUnit.test( "pointsToInside", function( assert ) {
	var face = PlaneFace.forVertices(P3(V3(0, 0, 1), 0), [V3(0, 0, 0), V3(10, 0, 0), V3(10, 10, 0), V3(0, 10, 0)])
	var v = V3(10, 10, 0)
	assert.ok(face.pointsToInside(v, V3(-1, -1, 0)))
	assert.notOk(face.pointsToInside(v, V3(1, 1, 0)))
	assert.notOk(face.pointsToInside(v, V3(-1, 0, 0)))
	assert.notOk(face.pointsToInside(v, V3(0, -1, 0)))
});
QUnit.assert.compareV3arraysLike = function (actual, expected, message) {
	this.push(expected.every((v, i) => v.like(actual[i])), actual.toSource(), expected.toSource(), message)
}
QUnit.test( "planeFaceEdgeISPsWithPlane", function(assert ) {
	var brep = B2.extrudeVertices([
			V3(3.318793683290006e-10, 9.12934465653568e-10, 3.8212177478055603e-10), // 0 0 0
			V3(9.999999999857609, 4.561232401125957e-13, 5.033029171552123e-10), // 10 0 0
			V3(2.984222284459465e-10, 9.999999999852161, -1.0461618780772852e-9)], // 0 10 0
		P3(V3(0, 0, 1), 0), V3(-1.8192047548394726e-10, -3.0150747681012967e-10, -5.0000000009123795), "ex0")
	/*var result = planeFaceEdgeISPsWithPlane(brep,
	 new BREP.Face([
	 V3(9.999999999675689, -3.010513535700171e-10, -5.000000000409076), // 10 0 -5
	 V3(1.4995889284505332e-10, 6.114269888434384e-10, -5.000000000530258), // 0 0 -5
	 V3(3.318793683290006e-10, 9.12934465653568e-10, 3.8212177478055603e-10), // 0 0 0
	 V3(9.999999999857609, 4.561232401125957e-13, 5.033029171552123e-10)], // 10 0 0
	 P3(V3(9.12478342464039e-11, 1, -6.03014953543423e-11), 9.129344656608087e-10)), // 0 1 0
	 L3(V3(-1.3833878355530264e-10, 6.114269894465992e-10, -4.999999990964091), V3(-1, 9.12478342480723e-11, 2.7667756772219476e-11)),
	 P3(V3(2.766775686256173e-11, 9.90075577448337e-10, 1), -4.999999990964091),
	 true, true)
	 assert.deepEqual(result, [])*/
	console.log(brep.faces[2].toSource())
	let line = L3.X.translate(0, 0, -1)
	var result = planeFaceEdgeISPsWithPlane(brep, brep.faces[2], L3.Y.translate(0, 0, -1), P3.XY.translate(0, 0, -1), true, true, new NLA.CustomSet()).map(is => is.p)
	assert.compareV3arraysLike(result, [V3(0, 0, -1), V3(0, 10, -1)])
	var result = planeFaceEdgeISPsWithPlane(brep, brep.faces[2], L3.Y, P3.XY, true, true, new NLA.CustomSet()).map(is => is.p)
	assert.compareV3arraysLike(result, [V3(0, 0, 0), V3(0, 10, 0)])
	var result = planeFaceEdgeISPsWithPlane(brep.translate(0, 0, 10), brep.translate(0, 0, 10).faces[2], L3.Y.translate(0, 0, 6), P3.XY.translate(0, 0, 6), true, true, new NLA.CustomSet()).map(is => is.p)
	assert.compareV3arraysLike(result, [V3(0, 0, 6), V3(0, 10, 6)])
});

QUnit.test( "B2.prototype.minus B2.box(5, 5, 5).minus(B2.box(1, 1, 6))", function( assert ) {
	var a = B2.box(5, 5, 5)
	var b = B2.box(1, 1, 6)
	var result = new B2([
		PlaneFace.forVertices(P3(V3(0, 0, -1), 0), [V3(0, 1, 0),V3(0, 5, 0),V3(5, 5, 0),V3(5, 0, 0),V3(1, 0, 0),V3(1, 1, 0)]),
		PlaneFace.forVertices(P3(V3(0, 0, 1), 5), [V3(0, 1, 5),V3(1, 1, 5),V3(1, 0, 5),V3(5, 0, 5),V3(5, 5, 5),V3(0, 5, 5)]),
		PlaneFace.forVertices(P3(V3(-1, 0, 0), 0), [V3(0, 1, 0),V3(0, 1, 5),V3(0, 5, 5),V3(0, 5, 0)]),
		PlaneFace.forVertices(P3(V3(0, 1, 0), 5), [V3(5, 5, 0),V3(0, 5, 0),V3(0, 5, 5),V3(5, 5, 5)]),
		PlaneFace.forVertices(P3(V3(1, 0, 0), 5), [V3(5, 0, 0),V3(5, 5, 0),V3(5, 5, 5),V3(5, 0, 5)]),
		PlaneFace.forVertices(P3(V3(0, -1, 0), 0), [V3(1, 0, 0),V3(5, 0, 0),V3(5, 0, 5),V3(1, 0, 5)]),
		PlaneFace.forVertices(P3(V3(0, -1, 0), -1), [V3(0, 1, 0),V3(1, 1, 0),V3(1, 1, 5),V3(0, 1, 5)]),
		PlaneFace.forVertices(P3(V3(-1, 0, 0), -1), [V3(1, 1, 0),V3(1, 0, 0),V3(1, 0, 5),V3(1, 1, 5)])])
	assert.b2Equal(a, b, a.minus(b), result)
});
QUnit.test( "B2.prototype.minus B2.box(5, 5, 5).minus(B2.box(1, 1, 5))", function( assert ) {
	var a = B2.box(5, 5, 5)
	var b = B2.box(1, 1, 5)
	var result = new B2([
		PlaneFace.forVertices(P3(V3(0, 0, -1), 0), [V3(0, 1, 0),V3(0, 5, 0),V3(5, 5, 0),V3(5, 0, 0),V3(1, 0, 0),V3(1, 1, 0)]),
		PlaneFace.forVertices(P3(V3(0, 0, 1), 5), [V3(0, 1, 5),V3(1, 1, 5),V3(1, 0, 5),V3(5, 0, 5),V3(5, 5, 5),V3(0, 5, 5)]),
		PlaneFace.forVertices(P3(V3(-1, 0, 0), 0), [V3(0, 1, 0),V3(0, 1, 5),V3(0, 5, 5),V3(0, 5, 0)]),
		PlaneFace.forVertices(P3(V3(0, 1, 0), 5), [V3(5, 5, 0),V3(0, 5, 0),V3(0, 5, 5),V3(5, 5, 5)]),
		PlaneFace.forVertices(P3(V3(1, 0, 0), 5), [V3(5, 0, 0),V3(5, 5, 0),V3(5, 5, 5),V3(5, 0, 5)]),
		PlaneFace.forVertices(P3(V3(0, -1, 0), 0), [V3(1, 0, 0),V3(5, 0, 0),V3(5, 0, 5),V3(1, 0, 5)]),
		PlaneFace.forVertices(P3(V3(0, -1, 0), -1), [V3(0, 1, 0),V3(1, 1, 0),V3(1, 1, 5),V3(0, 1, 5)]),
		PlaneFace.forVertices(P3(V3(-1, 0, 0), -1), [V3(1, 1, 0),V3(1, 0, 0),V3(1, 0, 5),V3(1, 1, 5)])])
	assert.b2Equal(a, b, a.minus(b), result)
});
QUnit.test( "B2 edge/face intersection", function( assert ) {
	var wideBox = B2.box(10, 10, 5)
	var extrusion = B2.extrudeVertices([V3(1,0),V3(0, 3), V3(-2, 5), V3(5, 5)], P3.XY.flipped(), V3(0, 0, 10), 'lol')
		.translate(0,1,1)
	var result = new B2([
		new PlaneFace(new PlaneSurface(P3(V3(0, 0, 1), 5)), [
			StraightEdge.throughPoints(V3(0, 3.999999999999999, 5), V3(0, 0, 5)),
			StraightEdge.throughPoints(V3(0, 0, 5), V3(10, 0, 5)),
			StraightEdge.throughPoints(V3(10, 0, 5), V3(10, 10, 5)),
			StraightEdge.throughPoints(V3(10, 10, 5), V3(0, 10, 5)),
			StraightEdge.throughPoints(V3(0, 10, 5), V3(0, 6, 5)),
			StraightEdge.throughPoints(V3(0, 6, 5), V3(5, 6, 5)),
			StraightEdge.throughPoints(V3(5, 6, 5), V3(1, 1, 5)),
			StraightEdge.throughPoints(V3(1, 1, 5), V3(0, 3.999999999999999, 5))]),
		new PlaneFace(new PlaneSurface(P3(V3(-1, 0, 0), 0)), [
			StraightEdge.throughPoints(V3(0, 6, 1), V3(0, 6, 5)),
			StraightEdge.throughPoints(V3(0, 6, 5), V3(0, 10, 5)),
			StraightEdge.throughPoints(V3(0, 10, 5), V3(0, 10, 0)),
			StraightEdge.throughPoints(V3(0, 10, 0), V3(0, 0, 0)),
			StraightEdge.throughPoints(V3(0, 0, 0), V3(0, 0, 5)),
			StraightEdge.throughPoints(V3(0, 0, 5), V3(0, 4.000000000000002, 5)),
			StraightEdge.throughPoints(V3(0, 3.999999999999999, 5), V3(0, 4, 1)),
			StraightEdge.throughPoints(V3(0, 4, 1), V3(0, 6, 1))]),
		new PlaneFace(new PlaneSurface(P3(V3(0, 0, -1), 0)), [
			StraightEdge.throughPoints(V3(0, 0, 0), V3(0, 10, 0)),
			StraightEdge.throughPoints(V3(0, 10, 0), V3(10, 10, 0)),
			StraightEdge.throughPoints(V3(10, 10, 0), V3(10, 0, 0)),
			StraightEdge.throughPoints(V3(10, 0, 0), V3(0, 0, 0))]),
		new PlaneFace(new PlaneSurface(P3(V3(0, 1, 0), 10)), [
			StraightEdge.throughPoints(V3(10, 10, 0), V3(0, 10, 0)),
			StraightEdge.throughPoints(V3(0, 10, 0), V3(0, 10, 5)),
			StraightEdge.throughPoints(V3(0, 10, 5), V3(10, 10, 5)),
			StraightEdge.throughPoints(V3(10, 10, 5), V3(10, 10, 0))]),
		new PlaneFace(new PlaneSurface(P3(V3(1, 0, 0), 10)), [
			StraightEdge.throughPoints(V3(10, 0, 0), V3(10, 10, 0)),
			StraightEdge.throughPoints(V3(10, 10, 0), V3(10, 10, 5)),
			StraightEdge.throughPoints(V3(10, 10, 5), V3(10, 0, 5)),
			StraightEdge.throughPoints(V3(10, 0, 5), V3(10, 0, 0))]),
		new PlaneFace(new PlaneSurface(P3(V3(0, -1, 0), 0)), [
			StraightEdge.throughPoints(V3(0, 0, 0), V3(10, 0, 0)),
			StraightEdge.throughPoints(V3(10, 0, 0), V3(10, 0, 5)),
			StraightEdge.throughPoints(V3(10, 0, 5), V3(0, 0, 5)),
			StraightEdge.throughPoints(V3(0, 0, 5), V3(0, 0, 0))]),
		new PlaneFace(new PlaneSurface(P3(V3(0, 0, -1), -1)), [
			StraightEdge.throughPoints(V3(0, 4, 1), V3(-2, 6, 1)),
			StraightEdge.throughPoints(V3(-2, 6, 1), V3(0, 6, 1)),
			StraightEdge.throughPoints(V3(0, 6, 1), V3(0, 4, 1))]),
		new PlaneFace(new PlaneSurface(P3(V3(-0.9486832980505139, -0.316227766016838, 0), -1.2649110640673518)), [
			StraightEdge.throughPoints(V3(1, 1, 5), V3(1, 1, 11)),
			StraightEdge.throughPoints(V3(1, 1, 11), V3(0, 4, 11)),
			StraightEdge.throughPoints(V3(0, 4, 11), V3(0, 4, 5)),
			StraightEdge.throughPoints(V3(0, 3.999999999999999, 5), V3(1, 1, 5))]),
		new PlaneFace(new PlaneSurface(P3(V3(-0.7071067811865477, -0.7071067811865472, 0), -2.8284271247461894)), [
			StraightEdge.throughPoints(V3(0, 4, 5), V3(0, 4, 11)),
			StraightEdge.throughPoints(V3(0, 4, 11), V3(-2, 6, 11)),
			StraightEdge.throughPoints(V3(-2, 6, 11), V3(-2, 6, 1)),
			StraightEdge.throughPoints(V3(-2, 6, 1), V3(0, 4, 1)),
			StraightEdge.throughPoints(V3(0, 4, 1), V3(0, 4.000000000000002, 5))]),
		new PlaneFace(new PlaneSurface(P3(V3(0, 1, 0), 6)), [
			StraightEdge.throughPoints(V3(0, 6, 5), V3(0, 6, 1)),
			StraightEdge.throughPoints(V3(0, 6, 1), V3(-2, 6, 1)),
			StraightEdge.throughPoints(V3(-2, 6, 1), V3(-2, 6, 11)),
			StraightEdge.throughPoints(V3(-2, 6, 11), V3(5, 6, 11)),
			StraightEdge.throughPoints(V3(5, 6, 11), V3(5, 6, 5)),
			StraightEdge.throughPoints(V3(5, 6, 5), V3(0, 6, 5))]),
		new PlaneFace(new PlaneSurface(P3(V3(0.7808688094430303, -0.6246950475544244, 0), 0.15617376188860604)), [
			StraightEdge.throughPoints(V3(5, 6, 5), V3(5, 6, 11)),
			StraightEdge.throughPoints(V3(5, 6, 11), V3(1, 1, 11)),
			StraightEdge.throughPoints(V3(1, 1, 11), V3(1, 1, 5)),
			StraightEdge.throughPoints(V3(1, 1, 5), V3(5, 6, 5))]),
		new PlaneFace(new PlaneSurface(P3(V3(0, 0, 1), 11)), [
			StraightEdge.throughPoints(V3(5, 6, 11), V3(-2, 6, 11)),
			StraightEdge.throughPoints(V3(-2, 6, 11), V3(0, 4, 11)),
			StraightEdge.throughPoints(V3(0, 4, 11), V3(1, 1, 11)),
			StraightEdge.throughPoints(V3(1, 1, 11), V3(5, 6, 11))])])
	assert.b2Equal(wideBox, extrusion, wideBox.plus(extrusion), result)
});
QUnit.test( "EllipseCurve", function( assert ) {
	var a = B2.box(5, 5, 5)
	var b = B2.box(1, 1, 5)
	var result = new B2([
		PlaneFace.forVertices(P3(V3(0, 0, -1), 0), [V3(0, 1, 0),V3(0, 5, 0),V3(5, 5, 0),V3(5, 0, 0),V3(1, 0, 0),V3(1, 1, 0)]),
		PlaneFace.forVertices(P3(V3(0, 0, 1), 5), [V3(0, 1, 5),V3(1, 1, 5),V3(1, 0, 5),V3(5, 0, 5),V3(5, 5, 5),V3(0, 5, 5)]),
		PlaneFace.forVertices(P3(V3(-1, 0, 0), 0), [V3(0, 1, 0),V3(0, 1, 5),V3(0, 5, 5),V3(0, 5, 0)]),
		PlaneFace.forVertices(P3(V3(0, 1, 0), 5), [V3(5, 5, 0),V3(0, 5, 0),V3(0, 5, 5),V3(5, 5, 5)]),
		PlaneFace.forVertices(P3(V3(1, 0, 0), 5), [V3(5, 0, 0),V3(5, 5, 0),V3(5, 5, 5),V3(5, 0, 5)]),
		PlaneFace.forVertices(P3(V3(0, -1, 0), 0), [V3(1, 0, 0),V3(5, 0, 0),V3(5, 0, 5),V3(1, 0, 5)]),
		PlaneFace.forVertices(P3(V3(0, -1, 0), -1), [V3(0, 1, 0),V3(1, 1, 0),V3(1, 1, 5),V3(0, 1, 5)]),
		PlaneFace.forVertices(P3(V3(-1, 0, 0), -1), [V3(1, 1, 0),V3(1, 0, 0),V3(1, 0, 5),V3(1, 1, 5)])])
	assert.B2equals(a.minus(b), result)
});
QUnit.test( "intersectionUnitCircleLine", function( assert ) {
	// y = -x + 1 => x + y = 1
	assert.deepEqual(intersectionUnitCircleLine(1, 1, 1), {x1: 1, x2: 0, y1:0, y2: 1})
});
QUnit.test( "intersectionCircleLine", function( assert ) {
	// y = -x + 2 => x + y = 2
	assert.deepEqual(intersectionCircleLine(1, 1, 2, 2), {x1: 2, x2: 0, y1:0, y2: 2})
});
QUnit.test( "CylinderSurface.intersectionLine", function( assert ) {
	var cylSurface = CylinderSurface.cylinder(5)
	var line = L3.throughPoints(V3(10, 0, 0), V3(-10, 2, 10))
	let isPoints = cylSurface.isTsForLine(line).map(line.at, line)

	assert.equal(isPoints.length, 2, 'no of points')
	assert.notOk(isPoints[0].like(isPoints[1]))

	assert.ok(cylSurface.containsPoint(isPoints[0]))
	assert.ok(cylSurface.containsPoint(isPoints[1]))

	assert.ok(line.containsPoint(isPoints[0]), line.distanceToPoint(isPoints[0]))
	assert.ok(line.containsPoint(isPoints[1]), line.distanceToPoint(isPoints[1]))
});

QUnit.test( "CylinderSurface.intersectionLine 2", function( assert ) {
	var cylSurface = new CylinderSurface(new EllipseCurve(V3.ZERO, V3(8, 0, 0), V3(0, 5, 0)), V3.Z)
	var line = L3.throughPoints(V3(10, 0, 0), V3(-10, 2, 10))
	var isPoints = cylSurface.isTsForLine(line).map(line.at, line)
	console.log(isPoints.toSource())
	assert.equal(isPoints.length, 2, 'no of points')
	assert.notOk(isPoints[0].like(isPoints[1]))

	assert.ok(cylSurface.containsPoint(isPoints[0]))
	assert.ok(cylSurface.containsPoint(isPoints[1]))

	assert.ok(line.containsPoint(isPoints[0]), line.distanceToPoint(isPoints[0]))
	assert.ok(line.containsPoint(isPoints[1]), line.distanceToPoint(isPoints[1]))
});
function testEllipseIntersections(assert, e1, e2, count) {
	var intersections = e1.isInfosWithEllipse(e2).map(info => info.p)
	assert.ok(intersections.length == count, `intersections.length == count: ${intersections.length} == ${count}`)
	intersections.forEach((is, i) => {
		assert.ok(intersections.every((is2, j) => j == i || !is.like(is2)), is.sce+' is not unique '+intersections)
		assert.ok(e1.containsPoint(is), `e1.containsPoint(is): ${e1.toSource()}.containsPoint(${is.sce})`)
		assert.ok(e2.containsPoint(is), `e2.containsPoint(is): ${e1.toSource()}.containsPoint(${is.sce})`)
	})
}
QUnit.test( "Ellipse.isPointsWithEllipse", function(assert ) {
	var c1 = EllipseCurve.circle(5), c2 = EllipseCurve.circle(5, V3(3, 0))
	testEllipseIntersections(assert, c1, c2, 2)
	var verticalEllipse = new EllipseCurve(V3(2, 0), V3(1, 1), V3(1, 10))
	testEllipseIntersections(assert, c1, verticalEllipse, 4)
	var verticalEllipse2 = new EllipseCurve(V3(10, 2), V3(1, 1), V3(1, 10))
	testEllipseIntersections(assert, c1, verticalEllipse2, 0)
	var smallEllipse = EllipseCurve.forAB(2, 3)
	testEllipseIntersections(assert, c1, smallEllipse, 0)
	var test = new EllipseCurve(V3(6, 1, 0), V3(3, 1, 0), V3(4, 0, 0))
	testEllipseIntersections(assert, c1, test, 2)
});
QUnit.test( "cylinder surface inss", function( assert ) {
	var cyl = CylinderSurface.cylinder(5)
	var ell = new CylinderSurface(new EllipseCurve(V3(6, 1, 4), V3(3, 1, 4), V3(4, 0, 0)), V3.Z)
	var iss = cyl.isCurvesWithSurface(ell)
	assert.ok(iss.length == 2)
	assert.ok(cyl.containsPoint(iss[0].anchor))
	assert.ok(ell.containsPoint(iss[0].anchor), ell.implicitFunction()(iss[0].anchor))
});
QUnit.test( "B2.box(10, 10, 5) + B2.box(10, 10, 5).translate(0, 0, 5) touching coplanar faces; result contains both", function( assert ) {
	var box = B2.box(10, 10, 5)
	assert.b2Equal(box, box.translate(0, 0, 5), box.plus(box.translate(0, 0, 5)), B2.box(10, 10, 10))
});
QUnit.test( "B2.box(10, 10, 5) + B2.box(10, 5, 5).translate(0, 0, 5) + B2.box(5, 5, 5).translate(5, 5, 5)", function( assert ) {
	var box = B2.box(10, 10, 5), box2 = B2.box(10, 5, 5), box3 = B2.box(5, 5, 5)
	assert.b2Equal(box.plus(box2.translate(0, 0, 5)), box3.translate(5, 5, 5), box.plus(box2.translate(0, 0, 5)).plus(box3.translate(5, 5, 5)),
		B2.box(10, 10, 10).minus(box3.translate(0, 5, 5)))
});

QUnit.test( "B2.box(10, 10, 5) && B2.box(10, 10, 5).translate(3, 3, 0) overlapping faces; result contains intersection of faces", function( assert ) {
	let box = B2.box(10, 10, 5)
	assert.b2Equal(box, box.translate(3, 3, 0), box.intersection(box.translate(3, 3, 0), true, true), B2.box(7, 7, 5).translate(3, 3, 0))
});
QUnit.test( "B2.box(10, 10, 5) + B2.box(10, 10, 5).translate(3, 3, 0) overlapping faces; result contains union of faces", function( assert ) {
	let box = B2.box(10, 10, 5)
	assert.b2Equal(box, box.translate(3, 3, 0), box.plus(box.translate(3, 3, 0), true, true), B2.box(10, 10, 10))
});


QUnit.test( "B2.box(10, 10, 5) + B2.box(10, 10, 5).translate(3, 3, 0) overlapping faces; result contains union of faces", function( assert ) {
	let box = B2.box(10, 10, 5), box2 = B2.box(4, 10, 2).translate(2, 0, 3)
	assert.b2Equal(box, box2, box.minus(box2), box)
});


QUnit.test( "serialization", function( assert ) {
	let a = {a: 2, b: 3}
	assert.equal(unserialize(serialize(a)).toString(), a.toString())

	a.c = a

	let a2 = unserialize(serialize(a))
	assert.equal(a2.a, 2)
	assert.equal(a2.b, 3)
	assert.equal(a2.c, a2)


	a = [1, 2, 3]
	assert.equal(unserialize(serialize(a)).toString(), a.toString())
});