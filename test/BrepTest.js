
QUnit.module('brep')

QUnit.test( "BREP.Face.equals", function( assert ) {
	var a =new BREP.Face([V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0),V3(0, 10, 0)], P3(V3(0, 0, 1), 0))
	var b =new BREP.Face([V3(0, 10, 0), V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0)], P3(V3(0, 0, 1), 0))
	var c =new BREP.Face([V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0),V3(0, 10, 0)].slice().reverse(), P3(V3(0, 0, -1), 0))

	assert.ok(a.equals(a))

	assert.ok(a.equals(b))
	assert.ok(b.equals(a))

	assert.notOk(a.equals(c))
	assert.notOk(c.equals(a))

	assert.notOk(b.equals(c))
	assert.notOk(c.equals(b))
});

QUnit.test( "BREP.equals", function( assert ) {
	var a = BREP.tetrahedron(V3(5, 5, 5), V3(5, 5, -5), V3(10, 12, 1), V3(0, 12, 1))
	var b = BREP.tetrahedron(V3(5, 5, 5), V3(10, 12, 1), V3(5, 5, -5), V3(0, 12, 1))
	var c = BREP.tetrahedron(V3(5, 5, 5), V3(12, 12, 1), V3(5, 5, -5), V3(0, 12, 1))

	assert.ok(a.equals(a))

	assert.ok(a.equals(b))
	assert.ok(b.equals(a))

	assert.notOk(a.equals(c))
	assert.notOk(c.equals(a))

	assert.notOk(b.equals(c))
	assert.notOk(c.equals(b))
});
QUnit.test( "BREP.isCCW", function( assert ) {
	var vertices =[V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0),V3(0, 10, 0)]
	assert.ok(isCCW(vertices, V3(0, 0, 1)))
	assert.notOk(isCCW(vertices, V3(0, 0, -1)))
});
QUnit.assert.brepEquals = function(actual, expected, message) {
	this.push(actual.equals(expected), actual.ss(), expected.ss(), message)
}
QUnit.assert.fuzzyEquals = function(actual, expected, message) {
	this.push(NLA.eq(actual, expected), actual, expected, message)
}
QUnit.test( "BREP.tetrahedon", function( assert ) {
	var result = new BREP(null, null, [
		new BREP.Face([V3(5, 5, 5),V3(5, 5, -5),V3(9, 11, 1)], P3(V3(0.8320502943378437, -0.5547001962252291, 0), 1.386750490563073)),
		new BREP.Face([V3(5, 5, 5),V3(1, 9, 0),V3(5, 5, -5)], P3(V3(-0.7071067811865475, -0.7071067811865475, 0), -7.071067811865475)),
		new BREP.Face([V3(5, 5, -5),V3(1, 9, 0),V3(9, 11, 1)], P3(V3(-0.1003911722115382, 0.7362019295512802, -0.6692744814102547), 6.525426193749983)),
		new BREP.Face([V3(9, 11, 1),V3(1, 9, 0),V3(5, 5, 5)], P3(V3(-0.25177250044296223, 0.6474150011390457, 0.7193500012656064), 5.57496250980845))])
	var a = BREP.tetrahedron(V3(5, 5, 5), V3(5, 5, -5), V3(9, 11, 1), V3(1, 9, 0))
	assert.brepEquals(a, result)
});

QUnit.test( "BREP.prototype.clipped() 1: no vertex/? intersections, no parallel faces", function( assert ) {
	var a =new BREP(null, null, [new BREP.Face([V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0),V3(0, 10, 0)], P3(V3(0, 0, 1), 0))])
	var b = BREP.tetrahedron(V3(5, 5, 5), V3(5, 5, -5), V3(10, 12, 1), V3(0, 12, 1))
	var result = new BREP(null, null, [new BREP.Face([
		V3(8.57142857142857, 9.999999999999998, 0),
		V3(4.999999999999999, 4.999999999999999, 0),
		V3(1.4285714285714297, 9.999999999999998, 0),
		V3(0, 10, 0),V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0)], P3(V3(0, 0, 1), 0))])
	assert.brepEquals(a.clipped(b), result)
});

QUnit.test( "BREP.prototype.clipped() 2: like 1, but clipping BREP vertices are on a plane (not inside polygon)", function( assert ) {
	var a =new BREP(null, null, [new BREP.Face([V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0),V3(0, 10, 0)], P3(V3(0, 0, 1), 0))])
	var b = BREP.tetrahedron(V3(5, 5, 5), V3(5, 5, -5), V3(10, 12, 0), V3(0, 12, 0))
	var result = new BREP(null, null, [new BREP.Face([V3(8.57142857142857, 9.999999999999998, 0),V3(4.999999999999999, 4.999999999999999, 0),V3(1.4285714285714297, 9.999999999999998, 0),V3(0, 10, 0),V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0)], P3(V3(0, 0, 1), 0))])
	assert.brepEquals(a.clipped(b), result)
});


QUnit.test( "BREP.prototype.clipped() 3: clipping BREP vertex is on face edge", function( assert ) {
	var a =new BREP(null, null, [new BREP.Face([V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0),V3(0, 10, 0)], P3(V3(0, 0, 1), 0))])
	var b = BREP.tetrahedron(V3(5, 5, 5), V3(5, 5, -5), V3(9, 10, 0), V3(0, 12, 1))
	var result = new BREP(null, null, [new BREP.Face([V3(9, 10, 0),V3(4.999999999999999, 5, 0),V3(1.4285714285714297, 9.999999999999998, 0),V3(0, 10, 0),V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0)], P3(V3(0, 0, 1), 0))])
	assert.brepEquals(a.clipped(b), result)
});

QUnit.test( "BREP.prototype.clipped() 4: clipping BREP edge is on face edge (+ non-zero overlap)", function( assert ) {
	var a =new BREP(null, null, [new BREP.Face([V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0),V3(0, 10, 0)], P3(V3(0, 0, 1), 0))])
	var b = BREP.tetrahedron(V3(5, 5, 5), V3(5, 5, -5), V3(9, 10, 0), V3(1, 10, 0))
	var result = new BREP(null, null, [new BREP.Face([V3(9, 10, 0),V3(4.999999999999999, 5, 0),V3(1, 10, 0),V3(0, 10, 0),V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0)], P3(V3(0, 0, 1), 0))])
	assert.brepEquals(a.clipped(b), result)
});

QUnit.test( "BREP.prototype.clipped() 5: projected edge is on face edge (+ non-zero overlap)", function( assert ) {
	var a =new BREP(null, null, [new BREP.Face([V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0),V3(0, 10, 0)], P3(V3(0, 0, 1), 0))])
	var b = BREP.tetrahedron(V3(5, 5, 5), V3(5, 5, -5), V3(9, 11, 1), V3(1, 10, 0))
	var result = new BREP(null, null, [new BREP.Face([V3(8.333333333333332, 9.999999999999998, 0),V3(5, 4.999999999999999, 0),V3(1, 10, 0),V3(0, 10, 0),V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0)], P3(V3(0, 0, 1), 0))])
	assert.brepEquals(a.clipped(b), result)
});

QUnit.test( "BREP.prototype.clipped() 6: projected poly completely inside face, projected vertex touches inside of face", function( assert ) {
	var a =new BREP(null, null, [new BREP.Face([V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0),V3(0, 10, 0)], P3(V3(0, 0, 1), 0))])
	var b = BREP.tetrahedron(V3(5, 5, 5), V3(5, 5, -5), V3(9, 11, 1), V3(1, 9, 0))
	var result = new BREP(null, null, [new BREP.Face([V3(8.333333333333332, 9.999999999999998, 0),V3(5, 4.999999999999999, 0),V3(1, 9, 0),V3(8.333333333333332, 9.999999999999998, 0),V3(0, 10, 0),V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0)], P3(V3(0, 0, 1), 0))])
	assert.brepEquals(a.clipped(b), result)
});

QUnit.test( "BREP.prototype.clipped() 7", function( assert ) {
	var a =new BREP(null, null, [new BREP.Face([V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0),V3(0, 10, 0)], P3(V3(0, 0, 1), 0))])
	var b = BREP.tetrahedron(V3(5, 5, 0), V3(8, 8, 0), V3(10, 10, 5), V3(1, 9, 5))
	var result = a
	assert.brepEquals(a.clipped(b), result)
});
QUnit.test( "BREP.prototype.clipped() 8", function( assert ) {
	var a =new BREP(null, null, [new BREP.Face([V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0),V3(0, 10, 0)], P3(V3(0, 0, 1), 0))])
	var b = BREP.tetrahedron(V3(-1, 4, 0), V3(-1, 6, 0), V3(5, 12, 5), V3(5, 12, -5))
	var result = new BREP(null, null, [
			new BREP.Face([V3(3.0000000000000004, 10, 0),V3(0, 10, 0),V3(0, 6.999999999999999, 0)], P3(V3(0, 0, 1), 0)),
			new BREP.Face([V3(3.5, 10, 0),V3(0, 5.333333333333334, 0),V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0)], P3(V3(0, 0, 1), 0))])

	assert.brepEquals(a.clipped(b), result)
});
QUnit.test( "BREP.prototype.clipped() 9 overlapping vertices", function( assert ) {
	var a =new BREP(null, null, [new BREP.Face([V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0),V3(0, 10, 0)], P3(V3(0, 0, 1), 0))])
	var b = BREP.tetrahedron(V3(10, 10, 0), V3(6, 7, 1), V3(7, 9, 1), V3(6, 7, -5)).newAll(
		BREP.tetrahedron(V3(10, 10, 0), V3(7, 6, 1), V3(9,7, 1), V3(7, 6, -5))
	)
	var result = new BREP(null, null, [new BREP.Face([V3(10, 10, 0),V3(8.666666666666668, 6.833333333333332, 0),V3(7, 6, 0),V3(10, 10, 0),V3(6, 7.000000000000001, 0),V3(6.833333333333334, 8.666666666666666, 0),V3(10, 10, 0),V3(0, 10, 0),V3(0, 0, 0),V3(10, 0, 0)], P3(V3(0, 0, 1), 0))])
	assert.brepEquals(a.clipped(b), result)
});
QUnit.test( "BREP.prototype.clipped() 10 no overlap", function( assert ) {
	var b = BREP.tetrahedron(V3(2, 2, 8), V3(3, 3, 8), V3(3, 2, -2), V3(2, 3, -2)).translate(5, 5, 0)
	var a = new BREP(null, null, [ new BREP.Face([V3(5, 0, 0), V3(5, 5, 0), V3(0, 5, 0), V3(0, 0, 0)].reverse(), P3(V3.Z.negated(), 0))])
	assert.brepEquals(a.clipped(b), a)
});
QUnit.test( "BREP.prototype.clipped() 11 a face corner touches inside of projected polygon edge", function( assert ) {
	var a =new BREP(null, null, [ new BREP.Face([V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0),V3(0, 10, 0)], P3(V3(0, 0, 1), 0))])
	var b = BREP.tetrahedron(V3(11, 1, 0), V3(9, -1, 0), V3(5, 5, 0), V3(5,5, 10))
	var result = new BREP(null, null, [
		new BREP.Face([V3(8.333333333333334, 0, 0),V3(5, 5, 0),V3(10, 1.6666666666666687, 0),V3(10, 10, 0),V3(0, 10, 0),V3(0, 0, 0)], P3(V3(0, 0, 1), 0))])
	assert.brepEquals(a.clipped(b, false, true), result)
});
QUnit.test( "BREP.prototype.clipped() 12, projected poly on inside of face, edge lies on face edge in the middle", function( assert ) {
	var b = BREP.tetrahedron(V3(1, 0, -1), V3(1, 0, 8), V3(1, 4, -1), V3(4, 0, -1))
	var a = new BREP(null, null, [ new BREP.Face([V3(5, 0, 0), V3(5, 5, 0), V3(0, 5, 0), V3(0, 0, 0)].reverse(), P3(V3.Z.negated(), 0))])
	var result = new BREP(null, null, [new BREP.Face([V3(1, 0, 0),V3(0, 0, 0),V3(0, 5, 0),V3(5, 5, 0),V3(5, 0, 0),V3(3.6666666666666665, -3.537972128266059e-16, 0),V3(1, 3.555555555555555, 0)], P3(V3(0, 0, -1), 0))])
	assert.brepEquals(a.clipped(b), result)
});
QUnit.test( "BREP.prototype.clipped() 13", function( assert ) {

	var a = BREP.extrude([V3(0, 0, 0), V3(10, 0, 0), V3(0, 10, 0)], P3.XY, V3(0, 0, -5), "ex0").flipped()
	//a = new BREP(null, null, [a.faces[0]])
	var b = BREP.extrude([V3(0, 0, -6), V3(0, 4, -6), V3(0, 0, -2)], P3.YZ.flipped().translate(0, 0, 0), V3(13, 0, 0), "ex1").flipped()
	b = new BREP(null, null, [b.faces[2]])

	var c = b.clipped(a)
	var result = new BREP(null, null, [new BREP.Face([V3(0, 2.999999999999999, -5),V3(7, 3, -5),V3(10, 0, -2.0000000000000004),V3(0, 0, -2)], P3(V3(0, -0.7071067811865476, -0.7071067811865476), 1.4142135623730954))])
	assert.brepEquals(b.clipped(a), result)
});
QUnit.test( "BREP.prototype.minus BREP.box(5, 5, 5).minus(BREP.box(1, 1, 6))", function( assert ) {
	var a = BREP.box(5, 5, 5)
	var b = BREP.box(1, 1, 6)
	var result = new BREP(null, null, [
		new BREP.Face([V3(0, 1, 0),V3(0, 5, 0),V3(5, 5, 0),V3(5, 0, 0),V3(1, 0, 0),V3(1, 1, 0)], P3(V3(0, 0, -1), 0)),
		new BREP.Face([V3(0, 1, 5),V3(1, 1, 5),V3(1, 0, 5),V3(5, 0, 5),V3(5, 5, 5),V3(0, 5, 5)], P3(V3(0, 0, 1), 5)),
		new BREP.Face([V3(0, 1, 0),V3(0, 1, 5),V3(0, 5, 5),V3(0, 5, 0)], P3(V3(-1, 0, 0), 0)),
		new BREP.Face([V3(5, 5, 0),V3(0, 5, 0),V3(0, 5, 5),V3(5, 5, 5)], P3(V3(0, 1, 0), 5)),
		new BREP.Face([V3(5, 0, 0),V3(5, 5, 0),V3(5, 5, 5),V3(5, 0, 5)], P3(V3(1, 0, 0), 5)),
		new BREP.Face([V3(1, 0, 0),V3(5, 0, 0),V3(5, 0, 5),V3(1, 0, 5)], P3(V3(0, -1, 0), 0)),
		new BREP.Face([V3(0, 1, 0),V3(1, 1, 0),V3(1, 1, 5),V3(0, 1, 5)], P3(V3(0, -1, 0), -1)),
		new BREP.Face([V3(1, 1, 0),V3(1, 0, 0),V3(1, 0, 5),V3(1, 1, 5)], P3(V3(-1, 0, 0), -1))])
	assert.brepEquals(a.minus(b), result)
});
QUnit.test( "BREP.prototype.minus BREP.box(5, 5, 5).minus(BREP.box(1, 1, 5))", function( assert ) {
	var a = BREP.box(5, 5, 5)
	var b = BREP.box(1, 1, 5)
	var result = new BREP(null, null, [
		new BREP.Face([V3(0, 1, 0),V3(0, 5, 0),V3(5, 5, 0),V3(5, 0, 0),V3(1, 0, 0),V3(1, 1, 0)], P3(V3(0, 0, -1), 0)),
		new BREP.Face([V3(0, 1, 5),V3(1, 1, 5),V3(1, 0, 5),V3(5, 0, 5),V3(5, 5, 5),V3(0, 5, 5)], P3(V3(0, 0, 1), 5)),
		new BREP.Face([V3(0, 1, 0),V3(0, 1, 5),V3(0, 5, 5),V3(0, 5, 0)], P3(V3(-1, 0, 0), 0)),
		new BREP.Face([V3(5, 5, 0),V3(0, 5, 0),V3(0, 5, 5),V3(5, 5, 5)], P3(V3(0, 1, 0), 5)),
		new BREP.Face([V3(5, 0, 0),V3(5, 5, 0),V3(5, 5, 5),V3(5, 0, 5)], P3(V3(1, 0, 0), 5)),
		new BREP.Face([V3(1, 0, 0),V3(5, 0, 0),V3(5, 0, 5),V3(1, 0, 5)], P3(V3(0, -1, 0), 0)),
		new BREP.Face([V3(0, 1, 0),V3(1, 1, 0),V3(1, 1, 5),V3(0, 1, 5)], P3(V3(0, -1, 0), -1)),
		new BREP.Face([V3(1, 1, 0),V3(1, 0, 0),V3(1, 0, 5),V3(1, 1, 5)], P3(V3(-1, 0, 0), -1))])
	assert.brepEquals(a.minus(b), result)
});


QUnit.test( "BREP.Face.prototype.containsPoint", function( assert ) {
	var a =new BREP.Face([V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0),V3(0, 10, 0)], P3(V3(0, 0, 1), 0))

	assert.ok(a.containsPoint(V3(5, 5, 0)))
	assert.notOk(a.containsPoint(V3(11, 5, 0)))


	var b = new BREP.Face([V3(0, 10, 0), V3(0, 0, 0),V3(5, 0, 0), V3(6, 5, 0)], P3(V3(0, 0, 1), 0))
	assert.ok(b.containsPoint(V3(2, 5, 0)))

	var c = new BREP.Face([V3(0, 10, 0), V3(0, 0, 0),V3(5, 0, 0), V3(6, 5, 0), V3(10, 5, 0)], P3(V3(0, 0, 1), 0))
	assert.ok(c.containsPoint(V3(2, 5, 0)))



	a = a.rotateZ(deg2rad(30))
	console.log(a.toString())
	assert.ok(a.containsPoint(V3(5, 5, 0).rotateZ(deg2rad(30))))
	assert.notOk(a.containsPoint(V3(-5, 5, 0).rotateZ(deg2rad(30))))

	b = b.rotateZ(deg2rad(30))
	assert.ok(b.containsPoint(V3(2, 5, 0).rotateZ(deg2rad(30))))

	c = c.rotateZ(deg2rad(30))
	assert.ok(c.containsPoint(V3(2, 5, 0).rotateZ(deg2rad(30))))
});

QUnit.test( "BREP.Face.prototype.containsPoint", function( assert ) {
	var a =new BREP.Face([V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0),V3(0, 10, 0)], P3(V3(0, 0, 1), 0))

	assert.notOk(a.containsPoint(V3(-0.00000001, 11, 0)))

	assert.ok(new BREP.Face([V3(-1, -10, 0), V3(0, 25, 0), V3(25, 0, 0)]).containsPoint(V3(0,0,0)))

});
QUnit.test( "BREP.Face.withHole", function( assert ) {
	var a =new BREP.Face([V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0),V3(0, 10, 0)], P3(V3(0, 0, 1), 0))
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
});

QUnit.test( "splitsVolumeEnclosingFaces", function( assert ) {
	var brep = BREP.tetrahedron(V3(0, 0,0),V3(10,0,0),V3(0,10,0),V3(0,0,10))
	// pointing into tetrahedon
	assert.ok(INSIDE == splitsVolumeEnclosingFaces(brep, [V3(0, 0,0),V3(10,0,0)], V3(0,1,1), V3(0,-1,1)))
	assert.ok(INSIDE == splitsVolumeEnclosingFaces(brep, [V3(0, 0,0),V3(10,0,0)], V3(0,1,1), V3(0,1,-1)))
	// pointing out of tetrahedon
	assert.notOk(INSIDE == splitsVolumeEnclosingFaces(brep, [V3(0, 0,0),V3(10,0,0)], V3(0,-1,0), V3(0,1,1)))

	assert.notOk(INSIDE == splitsVolumeEnclosingFaces(brep, [V3(0, 0,0),V3(10,0,0)], V3(0,-1,-1), V3(0,-1,1)))
	assert.notOk(INSIDE == splitsVolumeEnclosingFaces(brep, [V3(0, 0,0),V3(10,0,0)], V3(0,-1,-1), V3(0,1,-1)))
});


QUnit.test( "splitsVolumeEnclosingFaces 2", function( assert ) {
	var brep = BREP.tetrahedron(V3(0, 0,0),V3(10,0,0),V3(0,10,0),V3(0,0,10))
	// pointing out of tetrahedon
	assert.ok(splitsVolumeEnclosingFaces(brep, [V3(0, 0,0),V3(10,0,0)], V3(0,1,0), V3(0,0,1), false, true))
	assert.notOk(splitsVolumeEnclosingFaces(brep, [V3(0, 0,0),V3(10,0,0)], V3(0,1,0), V3(0,0,-1), false, true))
	assert.ok(splitsVolumeEnclosingFaces(brep, [V3(0, 0,0),V3(10,0,0)], V3(0,0,1), V3(0,1,0), false, true))
});
QUnit.test( "splitsVolumeEnclosingFaces 3", function( assert ) {
	var brep = BREP.box(5, 5, 5).flipped()
	// pointing out of tetrahedon
	assert.ok(splitsVolumeEnclosingFaces(brep, [V3(0, 5,0),V3(0,0,0)], V3(0,0,-1), V3(1,0,0)))
	assert.ok(splitsVolumeEnclosingFaces(brep, [V3(0,0,0),V3(0, 5,0)], V3(0,0,-1), V3(1,0,0)))

	assert.ok(splitsVolumeEnclosingFaces(brep, [V3(0,0,0),V3(0, 5,0)], V3(0,0,1), V3(-1,0,0), false, true))
	assert.notOk(splitsVolumeEnclosingFaces(brep, [V3(0,0,0),V3(0, 5,0)], V3(0,0,1), V3(1,0,0), false, true))
});

QUnit.test( "pointsToInside", function( assert ) {
	var face = new BREP.Face([V3(0, 0, 0),V3(10, 0, 0),V3(10, 10, 0),V3(0, 10, 0)], P3(V3(0, 0, 1), 0))
	var v = face.vertices[2]
	assert.ok(face.pointsToInside(v, V3(-1, -1, 0)))
	assert.notOk(face.pointsToInside(v, V3(1, 1, 0)))
	assert.notOk(face.pointsToInside(v, V3(-1, 0, 0)))
	assert.notOk(face.pointsToInside(v, V3(0, -1, 0)))
});
QUnit.test( "planeFaceEdgeISPsWithPlane 2", function( assert ) {
	var brep = BREP.extrude([
			V3(3.318793683290006e-10, 9.12934465653568e-10, 3.8212177478055603e-10), // 0 0 0
			V3(9.999999999857609, 4.561232401125957e-13, 5.033029171552123e-10), // 10 0 0
			V3(2.984222284459465e-10, 9.999999999852161, -1.0461618780772852e-9)], // 0 10 0
		P3(V3(0, 0, 1), 0), V3(-1.8192047548394726e-10, -3.0150747681012967e-10, -5.0000000009123795), "ex0").flipped()
	var result = planeFaceEdgeISPsWithPlane(brep,
	 new BREP.Face([
	 V3(9.999999999675689, -3.010513535700171e-10, -5.000000000409076), // 10 0 -5
	 V3(1.4995889284505332e-10, 6.114269888434384e-10, -5.000000000530258), // 0 0 -5
	 V3(3.318793683290006e-10, 9.12934465653568e-10, 3.8212177478055603e-10), // 0 0 0
	 V3(9.999999999857609, 4.561232401125957e-13, 5.033029171552123e-10)], // 10 0 0
	 P3(V3(9.12478342464039e-11, 1, -6.03014953543423e-11), 9.129344656608087e-10)), // 0 1 0
	 L3(V3(-1.3833878355530264e-10, 6.114269894465992e-10, -4.999999990964091), V3(-1, 9.12478342480723e-11, 2.7667756772219476e-11)),
	 P3(V3(2.766775686256173e-11, 9.90075577448337e-10, 1), -4.999999990964091),
	 true, true, new NLA.CustomSet()).map(is => is.p)
	assert.compareV3arraysLike(result, [])
	console.log(brep.faces[2])
	var result2 = planeFaceEdgeISPsWithPlane(brep, brep.faces[2], L3.X, P3.XY, true, true, new NLA.CustomSet())
	assert.compareV3arraysLike(result2, [])
	console.log(result2)
});
QUnit.assert.compareV3arraysLike = function (actual, expected, message) {
	this.push(expected.every((v, i) => v.like(actual[i])), actual.toSource(), expected.toSource(), message)
}
QUnit.test( "planeFaceEdgeISPsWithPlane", function(assert ) {
	var brep = BREP.extrude([
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
	console.log(brep.faces[2].vertices.map(v => v.toString(v => v.toFixed(3))).toSource())
	var result = planeFaceEdgeISPsWithPlane(brep, brep.faces[2], L3.X.translate(0, 0, -1), P3.XY.translate(0, 0, -1), true, true, new NLA.CustomSet()).map(is => is.p)
	assert.compareV3arraysLike(result, [V3(0, 0, -1), V3(10, 0, -1)])
	var result = planeFaceEdgeISPsWithPlane(brep, brep.faces[2], L3.X, P3.XY, true, true, new NLA.CustomSet()).map(is => is.p)
	assert.compareV3arraysLike(result, [V3(0, 0, 0), V3(10, 0, 0)])
	var result = planeFaceEdgeISPsWithPlane(brep.translate(0, 0, 10), brep.translate(0, 0, 10).faces[2], L3.X.translate(0, 0, 6), P3.XY.translate(0, 0, 6), true, true, new NLA.CustomSet()).map(is => is.p)
	assert.compareV3arraysLike(result, [V3(0, 0, 6), V3(10, 0, 6)])
});

