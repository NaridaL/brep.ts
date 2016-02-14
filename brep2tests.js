QUnit.assert.B2equals = function(actual, expected, message) {
	if (!(actual instanceof B2)) {
		this.push(false, actual, null, "actual is not a B2")
		return
	}
	console.log(actual)
	this.push(actual.equals(expected), actual.toString(), expected.toString(), message)
}
QUnit.assert.fuzzyEquals = function(actual, expected, message) {
	this.push(NLA.equals(actual, expected), actual, expected, message)
}


QUnit.test( "B2.Face.equals", function( assert ) {
	var a = B2.PlaneFace.forVertices(P3(V3(0, 0, 1), 0), [V3(0, 0, 0), V3(10, 0, 0), V3(10, 10, 0), V3(0, 10, 0)])
	var b = B2.PlaneFace.forVertices(P3(V3(0, 0, 1), 0), [V3(0, 10, 0), V3(0, 0, 0), V3(10, 0, 0), V3(10, 10, 0)])
	var c = B2.PlaneFace.forVertices(P3(V3(0, 0, -1), 0), [V3(0, 0, 0), V3(10, 0, 0), V3(10, 10, 0), V3(0, 10, 0)].slice().reverse())
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
QUnit.test( "B2.tetrahedon", function( assert ) {
	var result = new B2([
		B2.PlaneFace.forVertices(P3(V3(0.8320502943378437, -0.5547001962252291, 0), 1.386750490563073), [V3(5, 5, 5), V3(5, 5, -5), V3(9, 11, 1)]),
		B2.PlaneFace.forVertices(P3(V3(-0.7071067811865475, -0.7071067811865475, 0), -7.071067811865475), [V3(5, 5, 5), V3(1, 9, 0), V3(5, 5, -5)]),
		B2.PlaneFace.forVertices(P3(V3(-0.1003911722115382, 0.7362019295512802, -0.6692744814102547), 6.525426193749983), [V3(5, 5, -5), V3(1, 9, 0), V3(9, 11, 1)]),
		B2.PlaneFace.forVertices(P3(V3(-0.25177250044296223, 0.6474150011390457, 0.7193500012656064), 5.57496250980845), [V3(9, 11, 1), V3(1, 9, 0), V3(5, 5, 5)])])
	var a = B2.tetrahedron(V3(5, 5, 5), V3(5, 5, -5), V3(9, 11, 1), V3(1, 9, 0))
	assert.B2equals(a, result)
});


QUnit.test( "B2.Face.prototype.containsPoint", function( assert ) {
	var a = B2.PlaneFace.forVertices(P3(V3(0, 0, 1), 0), [V3(0, 0, 0), V3(10, 0, 0), V3(10, 10, 0), V3(0, 10, 0)])
	assert.ok(a.containsPoint(V3(5, 5, 0)))
	assert.notOk(a.containsPoint(V3(11, 5, 0)))


	var b = B2.PlaneFace.forVertices(P3(V3(0, 0, 1), 0), [V3(0, 10, 0), V3(0, 0, 0), V3(5, 0, 0), V3(6, 5, 0)])
	assert.ok(b.containsPoint(V3(2, 5, 0)))

	var c = B2.PlaneFace.forVertices(P3(V3(0, 0, 1), 0), [V3(0, 10, 0), V3(0, 0, 0), V3(5, 0, 0), V3(6, 5, 0), V3(10, 5, 0)])
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

QUnit.test( "B2.Face.prototype.containsPoint", function( assert ) {
	var a = B2.PlaneFace.forVertices(P3(V3(0, 0, 1), 0), [V3(0, 0, 0), V3(10, 0, 0), V3(10, 10, 0), V3(0, 10, 0)])
	assert.notOk(a.containsPoint(V3(-0.00000001, 11, 0)))

	assert.ok(B2.PlaneFace.forVertices(P3(V3(0, 0, -1), 0), [V3(-1, -10, 0), V3(0, 25, 0), V3(25, 0, 0)]).containsPoint(V3(0, 0, 0)))
});
QUnit.test( "B2.Face.withHole", function( assert ) {
	var a = B2.PlaneFace.forVertices(P3(V3(0, 0, 1), 0), [V3(0, 0, 0), V3(10, 0, 0), V3(10, 10, 0), V3(0, 10, 0)])
	var holeVertices = [V3(2, 3, 0),V3(8, 7, 0),V3(7, 2, 0)]


	assert.notOk(a.containsPoint(V3(-0.00000001, 11, 0)))

});
QUnit.test( "segmentTouchOrIntersect", function( assert ) {
	assert.ok(segmentsTouchOrIntersect(V3.ZERO, V3(3, 2.2), V3(2, 2.8), V3(2.8, 2.0)))

});
QUnit.test( "NLA.eqAngle", function( assert ) {
	assert.ok(NLA.zeroAngle(0))
	assert.ok(NLA.zeroAngle(- NLA.PRECISION / 2))
	assert.ok(NLA.zeroAngle(2 * Math.PI - NLA.PRECISION / 2))
	assert.ok(NLA.zeroAngle(2 * Math.PI + NLA.PRECISION / 2))
	assert.ok(NLA.eqAngle(-Math.PI, Math.PI))
	assert.ok(NLA.eqAngle(0, 2 * Math.PI))
	assert.ok(NLA.eqAngle(0, 2 * Math.PI - NLA.PRECISION / 2))
	assert.ok(NLA.eqAngle(0, 2 * Math.PI + NLA.PRECISION / 2))
	assert.notOk(NLA.eqAngle(-Math.PI, 2 * Math.PI))
	assert.notOk(NLA.eqAngle(0, Math.PI))

});

QUnit.test( "angleRelativeNormal", function( assert ) {
	assert.fuzzyEquals(V3.Y.angleRelativeNormal(V3(32, Math.sqrt(2), -Math.sqrt(2)), V3.X), -Math.PI / 4 )
	assert.fuzzyEquals(V3(-0.1, 1, 0).angleRelativeNormal(V3(0.0, 0, -1), V3.X), -Math.PI / 2)
});

QUnit.test( "splitsVolumeEnclosingFaces", function( assert ) {
	var brep = B2.tetrahedron(V3(0, 0,0),V3(10,0,0),V3(0,10,0),V3(0,0,10))
	// pointing into tetrahedon
	var edge = StraightEdge.throughPoints
	assert.ok(splitsVolumeEnclosingFaces(brep, edge(V3(0, 0,0),V3(10,0,0)), V3(0,1,1), V3(0,-1,1)))
	assert.ok(splitsVolumeEnclosingFaces(brep, edge(V3(0, 0,0),V3(10,0,0)), V3(0,1,1), V3(0,1,-1)))
	// pointing out of tetrahedon
	assert.notOk(splitsVolumeEnclosingFaces(brep, edge(V3(0, 0,0),V3(10,0,0)), V3(0,-1,0), V3(0,1,1)))

	assert.notOk(splitsVolumeEnclosingFaces(brep, edge(V3(0, 0,0),V3(10,0,0)), V3(0,-1,-1), V3(0,-1,1)))
	assert.notOk(splitsVolumeEnclosingFaces(brep, edge(V3(0, 0,0),V3(10,0,0)), V3(0,-1,-1), V3(0,1,-1)))
});


QUnit.test( "splitsVolumeEnclosingFaces 2", function( assert ) {
	var brep = B2.tetrahedron(V3(0, 0,0),V3(10,0,0),V3(0,10,0),V3(0,0,10))
	// pointing out of tetrahedon
	assert.ok(splitsVolumeEnclosingFaces(brep, StraightEdge.throughPoints(V3(0, 0,0),V3(10,0,0)), V3(0,1,0), V3(0,0,1), false, true))
	assert.notOk(splitsVolumeEnclosingFaces(brep, StraightEdge.throughPoints(V3(0, 0,0),V3(10,0,0)), V3(0,1,0), V3(0,0,-1), false, true))
	assert.ok(splitsVolumeEnclosingFaces(brep, StraightEdge.throughPoints(V3(0, 0,0),V3(10,0,0)), V3(0,0,1), V3(0,1,0), false, true))
});
QUnit.test( "splitsVolumeEnclosingFaces 3", function( assert ) {
	var brep = B2.box(5, 5, 5).flipped()
	console.log(brep.ss)
	// pointing out of tetrahedon
	assert.ok(splitsVolumeEnclosingFaces(brep, StraightEdge.throughPoints(V3(0, 5,0),V3(0,0,0)), V3(0,0,-1), V3(1,0,0)))
	assert.ok(splitsVolumeEnclosingFaces(brep, StraightEdge.throughPoints(V3(0,0,0),V3(0, 5,0)), V3(0,0,-1), V3(1,0,0)))

	assert.ok(splitsVolumeEnclosingFaces(brep, StraightEdge.throughPoints(V3(0,0,0),V3(0, 5,0)), V3(0,0,1), V3(-1,0,0), false, true))
	assert.notOk(splitsVolumeEnclosingFaces(brep, StraightEdge.throughPoints(V3(0,0,0),V3(0, 5,0)), V3(0,0,1), V3(1,0,0), false, true))
});

QUnit.test( "pointsToInside", function( assert ) {
	var face = B2.PlaneFace.forVertices(P3(V3(0, 0, 1), 0), [V3(0, 0, 0), V3(10, 0, 0), V3(10, 10, 0), V3(0, 10, 0)])
	var v = face.vertices[2]
	assert.ok(face.pointsToInside(v, V3(-1, -1, 0)))
	assert.notOk(face.pointsToInside(v, V3(1, 1, 0)))
	assert.notOk(face.pointsToInside(v, V3(-1, 0, 0)))
	assert.notOk(face.pointsToInside(v, V3(0, -1, 0)))
});
QUnit.assert.compareV3arraysLike = function (actual, expected, message) {
	this.push(expected.every((v, i) => v.like(actual[i])), actual.toSource(), expected.toSource(), message)
}
QUnit.test( "getFacePlaneIntersectionSs", function( assert ) {
	var brep = B2.extrudeVertices([
			V3(3.318793683290006e-10, 9.12934465653568e-10, 3.8212177478055603e-10), // 0 0 0
			V3(9.999999999857609, 4.561232401125957e-13, 5.033029171552123e-10), // 10 0 0
			V3(2.984222284459465e-10, 9.999999999852161, -1.0461618780772852e-9)], // 0 10 0
		P3(V3(0, 0, 1), 0), V3(-1.8192047548394726e-10, -3.0150747681012967e-10, -5.0000000009123795), "ex0")
	/*var result = getFacePlaneIntersectionSs(brep,
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
	var result = getFacePlaneIntersectionSs(brep, brep.faces[2], L3.X.translate(0, 0, -1), P3.XY.translate(0, 0, -1), true, true, new NLA.CustomSet()).map(is => is.p)
	assert.compareV3arraysLike(result, [V3(0, 0, -1), V3(10, 0, -1)])
	var result = getFacePlaneIntersectionSs(brep, brep.faces[2], L3.X, P3.XY, true, true, new NLA.CustomSet()).map(is => is.p)
	assert.compareV3arraysLike(result, [V3(0, 0, 0), V3(10, 0, 0)])
	var result = getFacePlaneIntersectionSs(brep.translate(0, 0, 10), brep.translate(0, 0, 10).faces[2], L3.X.translate(0, 0, 6), P3.XY.translate(0, 0, 6), true, true, new NLA.CustomSet()).map(is => is.p)
	assert.compareV3arraysLike(result, [V3(0, 0, 6), V3(10, 0, 6)])
});

QUnit.test( "B2.prototype.minus B2.box(5, 5, 5).minus(B2.box(1, 1, 6))", function( assert ) {
	var a = B2.box(5, 5, 5)
	var b = B2.box(1, 1, 6)
	var result = new B2([
		B2.PlaneFace.forVertices(P3(V3(0, 0, -1), 0), [V3(0, 1, 0),V3(0, 5, 0),V3(5, 5, 0),V3(5, 0, 0),V3(1, 0, 0),V3(1, 1, 0)]),
		B2.PlaneFace.forVertices(P3(V3(0, 0, 1), 5), [V3(0, 1, 5),V3(1, 1, 5),V3(1, 0, 5),V3(5, 0, 5),V3(5, 5, 5),V3(0, 5, 5)]),
		B2.PlaneFace.forVertices(P3(V3(-1, 0, 0), 0), [V3(0, 1, 0),V3(0, 1, 5),V3(0, 5, 5),V3(0, 5, 0)]),
		B2.PlaneFace.forVertices(P3(V3(0, 1, 0), 5), [V3(5, 5, 0),V3(0, 5, 0),V3(0, 5, 5),V3(5, 5, 5)]),
		B2.PlaneFace.forVertices(P3(V3(1, 0, 0), 5), [V3(5, 0, 0),V3(5, 5, 0),V3(5, 5, 5),V3(5, 0, 5)]),
		B2.PlaneFace.forVertices(P3(V3(0, -1, 0), 0), [V3(1, 0, 0),V3(5, 0, 0),V3(5, 0, 5),V3(1, 0, 5)]),
		B2.PlaneFace.forVertices(P3(V3(0, -1, 0), -1), [V3(0, 1, 0),V3(1, 1, 0),V3(1, 1, 5),V3(0, 1, 5)]),
		B2.PlaneFace.forVertices(P3(V3(-1, 0, 0), -1), [V3(1, 1, 0),V3(1, 0, 0),V3(1, 0, 5),V3(1, 1, 5)])])
	assert.B2equals(a.minus(b), result)
});
QUnit.test( "B2.prototype.minus B2.box(5, 5, 5).minus(B2.box(1, 1, 5))", function( assert ) {
	var a = B2.box(5, 5, 5)
	var b = B2.box(1, 1, 5)
	var result = new B2([
		B2.PlaneFace.forVertices(P3(V3(0, 0, -1), 0), [V3(0, 1, 0),V3(0, 5, 0),V3(5, 5, 0),V3(5, 0, 0),V3(1, 0, 0),V3(1, 1, 0)]),
		B2.PlaneFace.forVertices(P3(V3(0, 0, 1), 5), [V3(0, 1, 5),V3(1, 1, 5),V3(1, 0, 5),V3(5, 0, 5),V3(5, 5, 5),V3(0, 5, 5)]),
		B2.PlaneFace.forVertices(P3(V3(-1, 0, 0), 0), [V3(0, 1, 0),V3(0, 1, 5),V3(0, 5, 5),V3(0, 5, 0)]),
		B2.PlaneFace.forVertices(P3(V3(0, 1, 0), 5), [V3(5, 5, 0),V3(0, 5, 0),V3(0, 5, 5),V3(5, 5, 5)]),
		B2.PlaneFace.forVertices(P3(V3(1, 0, 0), 5), [V3(5, 0, 0),V3(5, 5, 0),V3(5, 5, 5),V3(5, 0, 5)]),
		B2.PlaneFace.forVertices(P3(V3(0, -1, 0), 0), [V3(1, 0, 0),V3(5, 0, 0),V3(5, 0, 5),V3(1, 0, 5)]),
		B2.PlaneFace.forVertices(P3(V3(0, -1, 0), -1), [V3(0, 1, 0),V3(1, 1, 0),V3(1, 1, 5),V3(0, 1, 5)]),
		B2.PlaneFace.forVertices(P3(V3(-1, 0, 0), -1), [V3(1, 1, 0),V3(1, 0, 0),V3(1, 0, 5),V3(1, 1, 5)])])
	assert.B2equals(a.minus(b), result)
});
QUnit.test( "EllipseCurve", function( assert ) {
	var a = B2.box(5, 5, 5)
	var b = B2.box(1, 1, 5)
	var result = new B2([
		B2.PlaneFace.forVertices(P3(V3(0, 0, -1), 0), [V3(0, 1, 0),V3(0, 5, 0),V3(5, 5, 0),V3(5, 0, 0),V3(1, 0, 0),V3(1, 1, 0)]),
		B2.PlaneFace.forVertices(P3(V3(0, 0, 1), 5), [V3(0, 1, 5),V3(1, 1, 5),V3(1, 0, 5),V3(5, 0, 5),V3(5, 5, 5),V3(0, 5, 5)]),
		B2.PlaneFace.forVertices(P3(V3(-1, 0, 0), 0), [V3(0, 1, 0),V3(0, 1, 5),V3(0, 5, 5),V3(0, 5, 0)]),
		B2.PlaneFace.forVertices(P3(V3(0, 1, 0), 5), [V3(5, 5, 0),V3(0, 5, 0),V3(0, 5, 5),V3(5, 5, 5)]),
		B2.PlaneFace.forVertices(P3(V3(1, 0, 0), 5), [V3(5, 0, 0),V3(5, 5, 0),V3(5, 5, 5),V3(5, 0, 5)]),
		B2.PlaneFace.forVertices(P3(V3(0, -1, 0), 0), [V3(1, 0, 0),V3(5, 0, 0),V3(5, 0, 5),V3(1, 0, 5)]),
		B2.PlaneFace.forVertices(P3(V3(0, -1, 0), -1), [V3(0, 1, 0),V3(1, 1, 0),V3(1, 1, 5),V3(0, 1, 5)]),
		B2.PlaneFace.forVertices(P3(V3(-1, 0, 0), -1), [V3(1, 1, 0),V3(1, 0, 0),V3(1, 0, 5),V3(1, 1, 5)])])
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