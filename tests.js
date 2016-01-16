var a = new GL.Vector(1, 1, 0);
var b = new GL.Vector(-2, -2, 0);
QUnit.test( "Vector.isParallelTo", function( assert ) {
  assert.ok( a.isParallelTo(b) );
});
QUnit.test( "undefined multiplication", function( assert ) {
	var boo = {};
  assert.ok( 1 == 1 + boo.far );
});
QUnit.test( "segmentIntersectsRay", function( assert ) {
	//var boo = VectorM.imm(1, 2, 3)
	assert.ok( "segmentIntersectsRay", segmentIntersectsRay(new CSG.Vector3D(0, 0, 0), new CSG.Vector3D(10, 0, 0), new CSG.Vector3D(5, 10, 0), new CSG.Vector3D(0, -5, 0)));
	assert.ok( "segmentIntersectsRay", !segmentIntersectsRay(new CSG.Vector3D(0, 0, 0), new CSG.Vector3D(10, 0, 0), new CSG.Vector3D(5, 10, 0), new CSG.Vector3D(0, 5, 0)));
	assert.ok( "segmentIntersectsRay", !segmentIntersectsRay(new CSG.Vector3D(0, 0, 0), new CSG.Vector3D(0, 10, 0), new CSG.Vector3D(0, 5, 10), new CSG.Vector3D(0, 0, 5)));
});
QUnit.test( "polygon intersectsLine", function( assert ) {
	//var boo = VectorM.imm(1, 2, 3)
	assert.ok(
		new CSG.Polygon([new CSG.Vector3D(10, 0, 0), new CSG.Vector3D(0, 20, 0), new CSG.Vector3D(0, 0, 30)].map((v) => new CSG.Vertex(v)))
			.intersectsLine(new CSG.Line3D(new CSG.Vector3D(0, 0, 0), new CSG.Vector3D(1, 1, 1))));
	assert.notOk(
		new CSG.Polygon([new CSG.Vector3D(10, 0, 0), new CSG.Vector3D(0, 20, 0), new CSG.Vector3D(0, 0, 30)].map((v) => new CSG.Vertex(v)))
			.intersectsLine(new CSG.Line3D(new CSG.Vector3D(11, 0, 0), new CSG.Vector3D(1, 1, 1))));
	assert.notOk(
		new CSG.Polygon([new CSG.Vector3D(0, 0, 0), new CSG.Vector3D(0, 20, 0), new CSG.Vector3D(0, 0, 30)].map((v) => new CSG.Vertex(v)))
			.intersectsLine(new CSG.Line3D(new CSG.Vector3D(11, 0, 0), new CSG.Vector3D(0, 1, 0))));
	assert.notOk(
		new CSG.Polygon([new CSG.Vector3D(100, 0, 0), new CSG.Vector3D(0, 0, 0), new CSG.Vector3D(87, 50, 0)].map((v) => new CSG.Vertex(v)))
			.intersectsLine(new CSG.Line3D(new CSG.Vector3D(107, -113, 0), new CSG.Vector3D(0, 0, -1))));
	assert.notOk(
		new CSG.Polygon([new CSG.Vector3D(87, 0, 50), new CSG.Vector3D(100, 0, 50), new CSG.Vector3D(100, 50, 0)].map((v) => new CSG.Vertex(v)))
			.intersectsLine(new CSG.Line3D(new CSG.Vector3D(267, 178, 990), new CSG.Vector3D(-0.402, -0.380, -0.833))));
	assert.notOk(
		new CSG.Polygon([new CSG.Vector3D(87, 50, 50), new CSG.Vector3D(100, 0, 50), new CSG.Vector3D(100, 0, 0)].map((v) => new CSG.Vertex(v)))
			.intersectsLine(new CSG.Line3D(new CSG.Vector3D(-355, 313, 955), new CSG.Vector3D(0.567, -0.41, -0.714))));
	assert.ok(
		new CSG.Polygon([new CSG.Vector3D(100, 0, 0), new CSG.Vector3D(0, 0, 0), new CSG.Vector3D(87, 50, 0)].map((v) => new CSG.Vertex(v)))
			.intersectsLine(new CSG.Line3D(new CSG.Vector3D(174, 182, -970), new CSG.Vector3D(-0.109, -0.169, 0.979))));
	assert.ok(
		new CSG.Polygon([new CSG.Vector3D(87, 50, 50), new CSG.Vector3D(100, 0, 0), new CSG.Vector3D(87, 50, 0)].map((v) => new CSG.Vertex(v)))
			.intersectsLine(new CSG.Line3D(new CSG.Vector3D(778, 508, 387), new CSG.Vector3D(-0.746, -0.525, -0.409))));
	var polyWall1_0 = new CSG.Polygon([[100, 0, 50], [0, 0, 50], [0, 0, 0]].map(a => new CSG.Vertex(new CSG.Vector3D(a))));
	var polyWall1_1 = new CSG.Polygon([[100, 0, 50], [0, 0, 0], [100, 0, 0]].map(a => new CSG.Vertex(new CSG.Vector3D(a))));
	assert.ok(
		polyWall1_0.intersectsLine(new CSG.Line3D(new CSG.Vector3D([196.06951107727852, -749.6852119549773, 633.861156159236]), new CSG.Vector3D([-0.18628252773580262, 0.7545516732654738, -0.6292460506293492] ))));
	assert.notOk(
		polyWall1_1.intersectsLine(new CSG.Line3D(new CSG.Vector3D([196.06951107727852, -749.6852119549773, 633.861156159236]), new CSG.Vector3D([-0.18628252773580262, 0.7545516732654738, -0.6292460506293492] ))));
});
QUnit.test( "polygon intersectsLine", function( assert ) {
	assert.notOk(
		new CSG.Polygon([new CSG.Vector3D(100, 0, 0), new CSG.Vector3D(0, 0, 0), new CSG.Vector3D(87, 50, 0)].map((v) => new CSG.Vertex(v)))
			.intersectsLine(new CSG.Line3D(new CSG.Vector3D(350, 134, -931), new CSG.Vector3D(-0.285, -0.163, 0.945))));
});
