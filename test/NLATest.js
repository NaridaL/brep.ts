window.onerror = function (errorMsg, url, lineNumber, column, errorObj) {
	console.log(errorMsg, url, lineNumber, column, errorObj);
}
QUnit.assert.matrixEquals = function(actual, expected, message, precision) {
	this.push(actual.equalsMatrix(expected, precision), actual.toString(), expected.toString(), message)
}
QUnit.assert.matrixEquivalent = function(actual, expected, message, precision) {
	this.push(actual.normalized2().equalsMatrix(expected.normalized2(), precision), actual.toString(), expected.toString(), message)
}
QUnit.assert.V3equals = function(actual, expected, message, precision) {
	this.push(false, actual.toString(), expected.toString(), message)
}
QUnit.assert.V3like = function(actual, expected, message, precision) {
	this.push(expected.like(actual), actual.toString(), expected.toString(), (message ?message:'V3like')+ "; |dv| = " + expected.distanceTo(actual))
}

QUnit.module('NLA')
/**
 *
 * @param name
 * @param {function (Object, M4)} test
 * @param {...M4} what
 *
 */
QUnit.testDifferentSystems = function (name, test, ...what) {
	if (!what.length) {
		what = [M4.IDENTITY, M4.FOO]
	}
	QUnit.test(name, assert => {
		what.forEach(m => {
			console.log(`TESTING '${name}' WITH '${m.name || m.toString()}`)
			assert.push(true,undefined,undefined, `TESTING '${name}' WITH '${m.name || m.toString()}'`)
			test(assert, m)
		})
	})
}
QUnit.test( "Vector.isParallelTo", function( assert ) {
	assert.equal(1, new Line(V3.ZERO, V3.X).distanceToPoint(V(1, 1, 0)))
});
QUnit.test( "NLA.eq etc", function( assert ) {
	assert.notOk(NLA.lt(2, 2 + NLA_PRECISION / 2))
	assert.notOk(NLA.lt(2, 2 - NLA_PRECISION / 2))
	assert.ok(NLA.le(2, 2 + NLA_PRECISION / 2))
	assert.ok(NLA.le(2, 2 - NLA_PRECISION / 2))

	assert.notOk(NLA.gt(2, 2 + NLA_PRECISION / 2))
	assert.notOk(NLA.gt(2, 2 - NLA_PRECISION / 2))
	assert.ok(NLA.ge(2, 2 + NLA_PRECISION / 2))
	assert.ok(NLA.ge(2, 2 - NLA_PRECISION / 2))

	assert.ok(NLA.lt(2, 3))
	assert.ok(NLA.gt(3, 2))
	assert.ok(NLA.le(2, 3))
	assert.ok(NLA.ge(3, 2))

	assert.ok(NLA.eq(2, 2 + NLA_PRECISION) == !NLA.gt(2, 2 + NLA_PRECISION))
	assert.ok(NLA.eq(2, 2 + NLA_PRECISION) == NLA.ge(2, 2 + NLA_PRECISION))
});
QUnit.test( "arrayCopyStep", function( assert ) {
	var a = [1, 2, 3, 4, 5, 6, 7, 8]
	var b = [-1, -2, -3, -4]
	NLA.arrayCopyStep(b, 0, 1, a, 1, 2, 3)
	assert.deepEqual(a, [1, -1, 3, -2, 5, -3, 7, 8])
});
QUnit.test( "arrayCopy", function( assert ) {
	var a = [1, 2, 3, 4, 5, 6, 7, 8]
	var b = [-1, -2, -3, -4]
	NLA.arrayCopy(b, 1, a, 2, 2)
	assert.deepEqual(a, [1, 2, -2, -3, 5, 6, 7, 8])
});
QUnit.test( "matrix rowArrays", function( assert ) {
	var a = Matrix.fromRowArrays([6, 3],[4, 3])
	assert.deepEqual(a.rowArray(0, Array), [6, 3])
	assert.deepEqual(a.rowArray(1, Array), [4, 3])
	assert.deepEqual(a.asRowArrays(Array), [[6, 3],[4, 3]])
});
QUnit.test( "matrix transposed", function( assert ) {
	var a = Matrix.fromRowArrays([6, 3],[4, 3])
	assert.deepEqual(a.transposed().asRowArrays(Array), [[6, 4],[3, 3]])
});
QUnit.test( "matrix transpose", function( assert ) {
	var a = Matrix.fromRowArrays([6, 3],[4, 3])
	a.transpose()
	assert.deepEqual(a.asRowArrays(Array), [[6, 4],[3, 3]])

	a = Matrix.fromRowArrays([6, 3, 4, 3])
	a.transpose()
	assert.deepEqual(a.asRowArrays(Array), [[6],[3],[4],[3]])
	a.transpose()
	assert.deepEqual(a.asRowArrays(Array), [[6, 3, 4, 3]])
});
QUnit.test( "Matrix.prototype.times", function( assert ) {
	var a = Matrix.fromRowArrays([6, 3],[4, 3])
	assert.deepEqual(Matrix.identity(2).times(a).asRowArrays(Array), [[6, 3],[4, 3]])
	assert.deepEqual(a.times(Matrix.identity(2)).asRowArrays(Array), [[6, 3],[4, 3]])
});
QUnit.test( "Matrix.identity", function( assert ) {
	var a = Matrix.identity(2)
	assert.deepEqual(a.asRowArrays(Array), [[1, 0],[0, 1]])
});
QUnit.test( "Matrix.prototype.rowsIndependent", function( assert ) {
	var a = Matrix.fromRowArrays([6, 3],[4, 3])
	var b = Matrix.fromRowArrays([6, 3],[12, 6])

	// all rows on plane through origin with normal = V(1, -1, 0)
	var c = Matrix.fromRowArrays([1, -1, 1],[1, 1, 1],[-1, 0, -1])
	assert.ok(a.rowsIndependent())
	assert.notOk(b.rowsIndependent())
	assert.notOk(c.rowsIndependent())
});
QUnit.test( "Matrix.prototype.rank", function( assert ) {
	var a = Matrix.fromRowArrays([6, 3],[4, 3])
	var b = Matrix.fromRowArrays([6, 3],[12, 6])

	// all rows on plane through origin with normal = V(1, -1, 0)
	var c = Matrix.fromRowArrays([1, -1, 1],[1, 1, 1],[-1, 0, -1])
	assert.equal(a.rank(), 2)
	assert.equal(b.rank(), 1)
	assert.equal(c.rank(), 2)

	var d = Matrix.fromRowArrays([1, 1, 0, 2],[-1, -1, 0, -2])
	assert.equal(d.rank(), 1)
	assert.equal(d.transposed().rank(), 1)


	let e = new M4(
		-60.16756109919886, 1, 1, 0,
		3, -56.16756109919886, 8, 0,
		21, 34, -5.167561099198863, 0,
		0, 0, 0, 1)
	console.log(e.rank())
	console.log(e.determinant())

});
QUnit.test( "LU Decomposition", function( assert ) {
	var a = Matrix.fromRowArrays([6, 3],[4, 3])
	var aLU = a.luDecomposition()
	assert.deepEqual(aLU.L.asRowArrays(Array), [[1, 0],[4/6, 1]])
	assert.deepEqual(aLU.U.asRowArrays(Array), [[6, 3],[0, 1]])
	assert.deepEqual(aLU.P.asRowArrays(Array), [[1, 0],[0, 1]])
	assert.matrixEquals(aLU.P.times(a), aLU.L.times(aLU.U))
});
QUnit.test( "LU Decomposition 2", function( assert ) {
	var a = Matrix.fromFunction(8, 8, (i, j) => Math.round((Math.random() - 0.5) * 4096))
	//var a = Matrix.fromRowArrays2([[-1636, 1740, -516], [-708, 403, 1986], [-1256, -1493, 996]])
	assert.ok(a.rowsIndependent())
	var aLU = a.luDecomposition()
	assert.ok(aLU.P.isPermutation())
	assert.ok(aLU.L.isLowerTriangular())
	assert.ok(aLU.U.isUpperTriangular())
	assert.matrixEquals(aLU.P.times(a), aLU.L.times(aLU.U))
});
QUnit.test( "LU Decomposition 3", function( assert ) {
	var a = Matrix.fromRowArrays([0,1,1],[1,1,1],[1,2,3])
	var aLU = a.luDecomposition()
	assert.deepEqual(aLU.U.asRowArrays(Array), [[1,1,1],[0,1,1],[0,0,1]])
	assert.deepEqual(aLU.L.asRowArrays(Array), [[1,0,0],[0,1,0],[1,1,1]])
	assert.deepEqual(aLU.P.asRowArrays(Array), [[0,1,0],[1,0,0],[0,0,1]])
	assert.matrixEquals(aLU.P.times(a), aLU.L.times(aLU.U))
});
QUnit.test( "Matrix.isOrthogonal", function( assert ) {
	var a = Matrix.identity(4)
	assert.ok(a.isOrthogonal())
	var b = Matrix.fromRowArrays([Math.sqrt(2) / 2, Math.sqrt(2) / 2],[-Math.sqrt(2) / 2, Math.sqrt(2) / 2])
	assert.ok(a.isOrthogonal())
});
QUnit.test( "Matrix.prototype.solveLinearSystem", function( assert ) {
	var a = Matrix.fromRowArrays([0,1,1],[1,1,1],[1,2,3])
	var b = NLA.Vector.from(1, 2, 3)
	var x = a.solveLinearSystem(b)
	assert.push(x.equals(NLA.Vector.from(1, 1, 0)), x, NLA.Vector.from(1, 1, 0))
	assert.push(a.timesVector(x).equals(b), a.timesVector(x), b)
});
QUnit.test( "Matrix.prototype.inverse", function( assert ) {
	var a = Matrix.fromRowArrays([0,1,1],[1,1,1],[1,2,3])
	var aInverse = a.inversed()
	console.log(aInverse.toString())
	assert.matrixEquals(a.times(aInverse), Matrix.identity(3))
});
QUnit.test( "Matrix.prototype.inverse 2", function( assert ) {
	var a = Matrix.random(8, 8)
	var aInverse = a.inversed()
	assert.matrixEquals(a.times(aInverse), Matrix.identity(8))
});
QUnit.test( "Matrix.prototype.inverse 3", function( assert ) {
	var a = new Matrix(1, 1, new Float64Array([5]))
	var aInverse = a.inversed()
	console.log(a.luDecomposition())
	assert.equal(aInverse.m[0], 1/5)
});
QUnit.test( "Matrix.prototype.qrDecompositionGivensRotation", function( assert ) {
	var sqrt = Math.sqrt
	var m = Matrix.fromRowArrays([3,5], [0,2], [0,0], [4,5])
	var {Q, R} = m.qrDecompositionGivensRotation()
	assert.matrixEquals(Q, Matrix.fromRowArrays(
		[3/5, 4/5/sqrt(5), 0, -8/5/sqrt(5)],
		[0, 2/sqrt(5), 0, 1/sqrt(5)],
		[0,0,1,0],
		[4/5, -3/5/sqrt(5), 0, 6/5/sqrt(5)]
	))
	assert.ok(Q.isOrthogonal())
	assert.matrixEquals(R, Matrix.fromRowArrays(
		[5,7],
		[0, sqrt(5)],
		[0,0],
		[0,0]
	))
	assert.matrixEquals(Q.times(R), Matrix.fromRowArrays([3,5], [0,2], [0,0], [4,5]))
});

QUnit.test( "Plane3.prototype.projectedVector", function( assert ) {
	var p = new P3(V(0,0,1), 2)
	assert.ok(V(1, 1, 0).like(p.projectedVector(V(1,1,1))))
});
QUnit.test( "Line3.prototype.distanceToLine", function( assert ) {
	assert.equal(L3.X.distanceToLine(new L3(V3.Z, V3.Y)), 1)
	assert.equal(L3.X.distanceToLine(new L3(V3.Z, V3.X)), 1)
});
QUnit.test( "Plane3.prototype.transformed", function( assert ) {
	var p = new P3(V(0,0,1), 2)
	assert.ok(P3.XY.like(P3.XY.transform(M4.identity())))
});
QUnit.test( "Matrix4x4.prototype.isAxisAligned", function( assert ) {
	assert.ok(M4.rotationX(Math.PI / 2).isAxisAligned())
	console.log(M4.rotationX(Math.PI / 4).toString())
	console.log(false + true + true)
	assert.notOk(M4.rotationX(Math.PI / 4).isAxisAligned())
});
QUnit.test( "Matrix4x4.prototype.rotationLine", function( assert ) {
	assert.matrixEquals(M4.rotationLine(V3.ZERO, V3.X, 1), M4.rotationX(1))
	assert.matrixEquals(M4.rotationLine(V3.ZERO, V3.Y, 1), M4.rotationY(1))
	assert.matrixEquals(M4.rotationLine(V3.ZERO, V3.Z, 1), M4.rotationZ(1))

	var a = V(1, 2, 3), PI = Math.PI;
	assert.matrixEquals(
		M4.rotationLine(a, V(1, 1, 0).normalized(), 1),
		M4.multiplyMultiple(M4.translation(a), M4.rotationZ(PI / 4), M4.rotationX(1), M4.rotationZ(-PI / 4), M4.translation(a.negated())),
		"",
		1e-6)
});
QUnit.test( "Matrix4x4.multiplyMultiple", function( assert ) {
	assert.matrixEquals(M4.multiplyMultiple(M4.rotationX(1), M4.rotationZ(1)), M4.rotationX(1).times(M4.rotationZ(1)))
});
QUnit.test( "Plane3.prototype.intersectionWithPlane", function( assert ) {
	assert.ok(P3.XY.intersectionWithPlane(P3.ZX).isColinearTo(L3.X))
	assert.ok(P3.ZX.intersectionWithPlane(P3.XY).isColinearTo(L3.X))
	assert.notOk(P3.ZX.intersectionWithPlane(P3.XY).isColinearTo(L3.Y))
});
QUnit.test( "Line3.prototype.isTsForLine", function(assert ) {
	console.log(L3.X.isInfoWithLine(new L3(V(1, 1, 0), V3.Y)).sce)
	assert.ok(L3.X.isInfoWithLine(new L3(V(1, 1, 0), V3.Y)).equals(V3.X))
});
QUnit.test( "V3.prototype.zip", function( assert ) {
	var a = V(1, 2, 3), b = V(4, 5, 6)
	assert.ok(V3.zip((a, b) => a + 3 * b, a, b).equals(a.plus(b.times(3))))
});
QUnit.test( "NLA.magic", function( assert ) {
	assert.expect(0)
	var a = V(1, 2, 3), b = V(4, 5, 6)
	//assert.ok(NLA.magic("a b c s => abs(a) x b .. c + 3 + ")(a, b, c, 3).equals(a.abs().cross(b).dot(c) + 3))
});
QUnit.test( "AABB", function( assert ) {
	var a = new AABB(V3.ZERO, V(20, 10, 30))
	var b = new AABB(V3.ZERO, V(1, 1, 1))
	assert.ok(a.likeAABB(a))
	assert.notOk(a.likeAABB(a.translate(10, 0, 0)))
	assert.ok(a.withoutAABB(b).likeAABB(new AABB(V(0, 0, 1), V(20, 10, 30))))
});
QUnit.test( "V3.areDisjoint", function(assert ) {

	assert.ok(V3.areDisjoint([V3.ZERO, V3.X, V3.Y].entries()))
	assert.ok(V3.areDisjoint([V3.ZERO, V3.X, V3.X, V3.Y].entries())) // same value twice is ok, as same reference
	assert.notOk(V3.areDisjoint([V3.ZERO, V3.X, V(0, 0, 0), V3.Y].entries())) // not ok as V3.ZERO != V(0, 0, 0)
	assert.notOk(V3.areDisjoint([V3.ZERO, V3.X, V(NLA_PRECISION / 2, 0, 0), V3.Y].entries())) // not ok as !V3.ZERO.like(V(...))
	assert.ok(V3.areDisjoint([V(NLA_PRECISION * -0.7, 0, 0), V(NLA_PRECISION * 0.7, 0, 0)].entries())) // not ok as V3.ZERO != V(0, 0, 0)
});
QUnit.test( "V3.areDisjoint2", function(assert ) {
	console.log(~~2147483657)
	var s = new NLA.CustomSet()
	var a = V(0, 2.7499999999999996, -5), b = V(0, 2.749999999999999, -5)
	s.canonicalizeLike(a)
	console.log(s._map, a.like(b), a.hashCodes(), b.hashCodes(), a.hashCode(), b.hashCode())
	assert.ok(s.canonicalizeLike(b) == a)
});
QUnit.test( "NLA.arrayBinaryInsert", function(assert ) {
	var arr = [1, 2, 3, 4]
	NLA.arrayBinaryInsert(arr, 2.5, (a, b) => a - b)
	assert.deepEqual(arr, [1, 2, 2.5, 3, 4])

	var arr2 = []
	NLA.arrayBinaryInsert(arr2, -2, NLA.minus)
	NLA.arrayBinaryInsert(arr2, 5, NLA.minus)
	assert.deepEqual(arr2, [-2, 5])
});
QUnit.test( "NLA.arrayBinaryIndexOf", function(assert ) {
	assert.equal([1, 2, 2.5, 3, 4].binaryIndexOf(3, (a, b) => a - b), 3)
	assert.equal([1, 2, 2.5, 3, 4].binaryIndexOf(2.6, (a, b) => a - b), -3-1)
});
QUnit.test( "newtonIterate2d", function (assert) {
	var res = newtonIterate2d((s, t) => s - 2, (s, t) => t - 4, 5, 5)
	assert.push(res.like(V(2, 4)), res.sce, V(2, 4).sce)
});

QUnit.test( "NLA.M4.projection", function (assert) {
	var plane = new P3(V(1, 2, 3).normalized(), 5);
	var proj = M4.projection(plane)
	console.log(proj.transformPoint(V(2, 4, 6)))
	assert.V3like(proj.transformPoint(V(2, 4, 6)), plane.anchor)
	assert.V3like(proj.transformVector(V(2, 4, 6)), V3.ZERO)
	var p2 = V(3, 5, 22)
	assert.ok(proj.transformPoint(p2).minus(p2).isParallelTo(plane.normal))
	assert.ok(plane.containsPoint(proj.transformPoint(p2)))
	assert.ok(proj.transformVector(p2).minus(p2).isParallelTo(plane.normal))
	assert.ok(proj.transformVector(p2).isPerpendicularTo(plane.normal))
});
QUnit.test( "NLA.M4.projection 2", function (assert) {
	[V(1, 1, 1), V(0, 0, -1)].forEach(dir => {
		var plane = new P3(V(1, 2, 3).normalized(), 5);
		var proj = M4.projection(plane, dir)
		var p2 = V(3, 5, 22)
		assert.ok(proj.transformPoint(p2).minus(p2).isParallelTo(dir))
		assert.ok(plane.containsPoint(proj.transformPoint(p2)))
		assert.ok(proj.transformVector(p2).minus(p2).isParallelTo(dir))
		assert.ok(proj.transformVector(p2).isPerpendicularTo(plane.normal))
		console.log(proj.transformPoint(p2).sce)
		console.log(proj.str)
	})
});
QUnit.test( "NLA.M4.isIdentity", function (assert) {
	assert.ok(M4.identity().isIdentity())
	assert.notOk(M4.scaling(1, 2, 3).isIdentity())
});
QUnit.test( "ParabolaCurve", function (assert) {
	var curve = new ParabolaCurve(V(1, 1), V(4, 1, -2), V(1, 10, 2))
	assert.ok(curve.containsPoint(curve.at(0)))
	assert.ok(curve.containsPoint(curve.at(1)))
	assert.ok(curve.containsPoint(curve.at(-1)))
	var plane = new P3(V(2,7,1).normalized(), 10)
	var iss = curve.isTsWithPlane(plane)
	assert.equal(iss.length, 2)
	assert.ok(plane.containsPoint(curve.at(iss[0])), plane.distanceToPointSigned(curve.at(iss[0])))
	assert.ok(plane.containsPoint(curve.at(iss[1])), plane.distanceToPointSigned(curve.at(iss[1])))


	var curveRA = curve.rightAngled()
	NLA.arrayRange(-10, 10, 1).forEach(t => assert.ok(curveRA.containsPoint(curve.at(t))))

	//var curve = ParabolaCurve.forAB(10, 20)
	var startT = -2, endT = 3, steps = 1000
	console.log(integrateCurve(curve, startT, endT, 1000))
	console.log(curve.arcLength(startT, endT))
});
QUnit.test( "BezierCurve", function (assert) {
	var curve = BezierCurve.graphXY(2,-3,-3,2)//.rotateZ(PI/3)
	//NLA.arrayRange(-1, 1, 0.1).forEach(t => assert.ok(NLA.eq(curve.at(t).x, t)))
	//curve.pointLambda(V3.ZERO)//(V(1,-2))
	NLA.arrayRange(-1, 1, 0.1).forEach(t => assert.push(NLA.eq(t, curve.pointLambda(curve.at(t))), t, curve.pointLambda(curve.at(t))))
});
QUnit.test( "BezierCurve.distanceToPoint", function (assert) {
	var curve = BezierCurve.graphXY(0,0,0,1)//.rotateZ(PI/3)
//        assert.ok(NLA.eq2(curve.distanceToPoint(V(0.5, 0)), 1, NLA_PRECISION))

	let curve2 = BezierCurve.graphXY(2,-3,-3,2)
	let p = V(0.5, 0.2)
	let closestT = curve2.closestTToPoint(p)
	let pDist = curve2.at(closestT).distanceTo(p)
	const EPS = NLA_PRECISION
	assert.push(pDist < curve2.at(closestT - EPS).distanceTo(p), curve2.at(closestT - EPS).distanceTo(p), "> " + pDist, "" + (pDist - curve2.at(closestT - EPS).distanceTo(p)) + "larger")
	assert.push(pDist < curve2.at(closestT + EPS).distanceTo(p), curve2.at(closestT + EPS).distanceTo(p), "> " + pDist)

	let curve3 = BezierCurve.graphXY(2,-3,-3,2).scale(100, 100, 100)
	let p3 = V(71, -65, 0)
	let closestT3 = curve3.closestTToPoint(p3)
	let pDist3 = curve3.at(closestT3).distanceTo(p3)
	assert.push(pDist3 < curve3.at(closestT3 - EPS).distanceTo(p3), curve3.at(closestT3 - EPS).distanceTo(p3), "> " + pDist3, "" + (pDist3 - curve3.at(closestT3 - EPS).distanceTo(p3)) + "larger")
	assert.push(pDist3 < curve3.at(closestT3 + EPS).distanceTo(p3), curve3.at(closestT3 + EPS).distanceTo(p3), "> " + pDist3)

});
QUnit.test( "EllipseCurve.distanceToPoint", function (assert) {
	let curve = EllipseCurve.forAB(10, 15)
	let p = V(12, 12)
	let closestT = curve.closestTToPoint(p)
	let pDist = curve.at(closestT).distanceTo(p)
	const EPS = 0.001
	assert.push(pDist < curve.at(closestT - EPS).distanceTo(p), curve.at(closestT - EPS).distanceTo(p), "> " + pDist, "" + (pDist - curve.at(closestT - EPS).distanceTo(p)) + "larger")
	assert.push(pDist < curve.at(closestT + EPS).distanceTo(p), curve.at(closestT + EPS).distanceTo(p), "> " + pDist)
});
QUnit.test( "BezierCurve.isInfosWithLine", function (assert) {
	console.log(solveCubicReal2(1, 0, 0, 0))
	let curve = BezierCurve.graphXY(2,-3,-3,2,   -3, 4)
	let line = new L3(V3.Y, V3.X)
	let isInfos = curve.isInfosWithLine(line.anchor, line.dir1)
	assert.equal(isInfos.length, 3)
	isInfos.forEach(info => {
		let p = info.p
		assert.ok(line.containsPoint(p))
		assert.ok(curve.containsPoint(p))
	})


	let line2 = new L3(V(0, 2, 1), V3.Z)
	let isInfos2 = curve.isInfosWithLine(line2.anchor, line2.dir1)
	assert.equal(isInfos2.length, 1)
	assert.ok(V(0, 2, 0).like(isInfos2[0].p))


	let line3 = new L3(V3.Z, V3.X)
	assert.equal(curve.isInfosWithLine(line3.anchor, line3.dir1).length, 0)

});
QUnit.test( "BezierCurve.isTsWithPlane", function (assert) {
	let curve = BezierCurve.graphXY(2,-3,-3,2)
	let plane = new P3(V(0, 1, 1).normalized(), 1)

	let iss = curve.isTsWithPlane(plane)
	assert.equal(iss.length, 3)
	iss.forEach(t => {
		assert.ok(plane.containsPoint(curve.at(t)))
	})
});

QUnit.test( "solveCubicReal2()", function (assert) {
	assert.deepEqual(solveCubicReal2(0, 1, 0, 0), [0])
	assert.deepEqual(solveCubicReal2(1, 0, 0, 0), [0])
});

QUnit.test( "M4.projectPlanePoint()", function (assert) {
	var m4 = M4.projectPlanePoint(V3.Z.negated(), P3.XY)
	assert.V3like(m4.transformPoint(V(4, 0, 1)), V(2, 0, 0))
	assert.V3like(m4.transformPoint(V(4, 8, 1)), V(2, 4, 0))
	assert.V3like(m4.transformPoint(V(4, 8, 2)), V(4/3, 8/3, 0))
	assert.matrixEquivalent(
		M4.projectPlanePoint(M4.FOO.transformPoint(V3.Z.negated()), P3.XY.transform(M4.FOO)),
		M4.multiplyMultiple(M4.FOO, m4, M4.BAR))

});
QUnit.test( "ConicSurface.coplanar", function (assert) {
	var unitCone = ConicSurface.UNIT
	assert.ok(unitCone.matrix.isIdentity(), 'unitCone.matrix.isIdentity()')
	assert.V3like(unitCone.parametricFunction()(0,3), V(3, 0, 3))
	var ellipseAtZ3 = EllipseCurve.UNIT.scale(3, 3, 3).translate(0,0,3)
	var planeAtZ3 = P3.XY.translate(0, 0, 3)
	var issAtZ3 = unitCone.isCurvesWithPlane(planeAtZ3)
	assert.equal(issAtZ3.length, 1)
	assert.push(ellipseAtZ3.isColinearTo(issAtZ3[0]), issAtZ3.toString(), ellipseAtZ3.toString())
	assert.ok(unitCone.containsEllipse(ellipseAtZ3))


	var scaledUnit = ConicSurface.UNIT.scale(2, 2, 1)
	assert.notOk(scaledUnit.isCoplanarTo(unitCone))
	assert.notOk(unitCone.isCoplanarTo(scaledUnit))
	var ell1 = unitCone.isCurvesWithPlane(new P3(V(2, 3, 10).normalized(), 10))[0]
	assert.ok(unitCone.containsEllipse(ell1), 'unitCone.containsEllipse(ell1)')
	var ell2 = unitCone.isCurvesWithPlane(new P3(V(1, 1, 2).normalized(), 4))[0]
	var ell1Cone = ConicSurface.atApexThroughEllipse(V3.ZERO, ell1, 1)
	var ell2Cone = ConicSurface.atApexThroughEllipse(V3.ZERO, ell2, 1)
	console.log(ell1Cone   )
	assert.ok(unitCone.isCoplanarTo(ell1Cone))
	assert.ok(unitCone.isCoplanarTo(ell2Cone))
	assert.ok(ell1Cone.isCoplanarTo(ell2Cone))
	assert.ok(ell2Cone.isCoplanarTo(ell1Cone))
	assert.ok(ell1Cone.foo().isCoplanarTo(ell2Cone.foo()))
});
QUnit.test( "ConicSurface.containsParabola", function (assert) {
	var unitCone = ConicSurface.UNIT
	var pb = unitCone.isCurvesWithPlane(new P3(V(1,0,1).normalized(), 4))[0]
	assert.ok(unitCone.containsParabola(pb))

	var c2 = unitCone.shearedX(2, 3)
	var pb2 = c2.isCurvesWithPlane(new P3(V(1,0,1).normalized(), 4).shearedX(2, 3))[0]
	assert.ok(c2.containsParabola(pb2))
});
QUnit.test( "Hyperbola", function (assert) {
	var hb = HyperbolaCurve.UNIT
	NLA.arrayRange(-10, 10, 1).forEach(t => assert.ok(hb.containsPoint(hb.at(t))))

	var hbSheared = hb.shearedX(2, 3)
	assert.notOk(hbSheared.isOrthogonal())
	var hbScaledRA = hbSheared.rightAngled()
	NLA.arrayRange(-10, 10, 1).forEach(t => assert.ok(hbScaledRA.containsPoint(hbSheared.at(t))))

	assert.deepEqual(intersectionUnitHyperbolaLine(1,0,2), {x1:2,y1:sqrt(3),x2:2,y2:-sqrt(3)})
});
QUnit.testDifferentSystems( "ProjectedCurveSurface", function (assert, m4) {
	let pp = V(0.5, 1)
	let curve = BezierCurve.graphXY(2,-3,-3,2)
	let pcs = new ProjectedCurveSurface(curve, V3.Z).transform(m4)
	let p = pcs.parametricFunction()(pp.x, pp.y)
	console.log(p.sce, pcs.pointToParameterFunction()(p))
	assert.V3like(pcs.pointToParameterFunction()(p), pp, "ptpf(pcs.pf(pp)) == pp")
});
QUnit.testDifferentSystems( "ProjectedCurveSurface Face line intersection test", function (assert, m4) {
	let curve = BezierCurve.graphXY(2,-3,-3,2)
	let edge = PCurveEdge.forCurveAndTs(curve, 0, 1)
	let edges = [
		edge,
		StraightEdge.throughPoints(curve.at(1), curve.at(1).plus(V(0, 0, 10))),
		edge.flipped().transform(M4.translation(0, 0, 10)),
		StraightEdge.throughPoints(curve.at(0).plus(V(0, 0, 10)), curve.at(0))]
	let surface = new ProjectedCurveSurface(curve, V3.Z)
	let face = new Face(surface, edges).transform(m4)
	let line = new L3(V3.Z, V3.X).transform(m4)
	let d = face.intersectsLine(line)
	assert.ok(d)
});
QUnit.test( "ProjectedCurveSurface Face containsPoint", function (assert, m4) {
	let face = new RotationFace(new ProjectedCurveSurface(new BezierCurve(V(142.87578921496748, -191.46078243076332, 0), V(161.78547089700214, -252.13248349581008, 0), V(284.63214994898954, -163.59789158697575, 0), V(372.40411211189405, -210.3992206435476, 0), -3, 4), V(0, 0, 1), -3, 4), [
		PCurveEdge.forCurveAndTs(
			new BezierCurve(V(142.87578921496748, -191.46078243076332, 0), V(161.78547089700214, -252.13248349581008, 0), V(284.63214994898954, -163.59789158697575, 0), V(372.40411211189405, -210.3992206435476, 0)), 1, 0),
		StraightEdge.throughPoints(V(142.87578921496748, -191.46078243076332, 0), V(142.87578921496748, -191.46078243076332, -100)),
		PCurveEdge.forCurveAndTs(new BezierCurve(V(142.87578921496748, -191.46078243076332, -100), V(161.78547089700214, -252.13248349581008, -100), V(284.63214994898954, -163.59789158697575, -100), V(372.40411211189405, -210.3992206435476, -100)), 0, 1),
		StraightEdge.throughPoints(V(372.40411211189405, -210.3992206435476, -100), V(372.40411211189405, -210.3992206435476, 0))], [])
	let line = new L3(V(1241.5987, -1214.1894, 38.9886), V(-0.6705, 0.7386, -0.0696).normalized())
	let ists = face.surface.isTsForLine(line)
	assert.equal(ists.length, 3)
	ists.forEach(t => {
		let p = line.at(t)
		assert.ok(face.surface.containsPoint(p))
	})
//		let p = V(1192.4056247755673, -1243.899135769775, 220.80458903468156)
//		assert.ok(face.surface.containsPoint(p))
});
QUnit.test( "CylinderSurface Face containsPoint", function (assert, m4) {
	let face = new RotationFace(new CylinderSurface(new EllipseCurve(V(73.03583314037537, -69.86032483338774, 0), V(-24.176861672352132, -146.16681457389276, 0), V(146.16681457389276, -24.176861672352132, 0)), V(0, 0, 1)), [
		new PCurveEdge(new EllipseCurve(V(73.03583314037537, -69.86032483338774, 0), V(-24.176861672352132, -146.16681457389276, 0), V(146.16681457389276, -24.176861672352132, 0)), V(219.75148278474705, -90.44615667066816, 0), V(97.2126948127275, 76.30648974050503, 0), 1.5953170840348225, -3.141592653589793, null, V(-20.58583183728038, -146.71564964437164, 0), V(146.16681457389276, -24.176861672352114, 0)),
		StraightEdge.throughPoints(V(97.2126948127275, 76.30648974050503, 0), V(97.2126948127275, 76.30648974050503, -100)),
		new PCurveEdge(new EllipseCurve(V(73.03583314037537, -69.86032483338774, -100), V(-24.176861672352132, -146.16681457389276, 0), V(146.16681457389276, -24.176861672352132, 0)), V(97.2126948127275, 76.30648974050503, -100), V(219.75148278474705, -90.44615667066816, -100), -3.141592653589793, 1.5953170840348225, null, V(-146.16681457389276, 24.176861672352114, 0), V(20.58583183728038, 146.71564964437164, 0)),
		StraightEdge.throughPoints(V(219.75148278474705, -90.44615667066816, -100), V(219.75148278474705, -90.44615667066816, 0))], [])
	// let line = new L3(V(-1344.04574670165, 826.5930889273866, 720.915318266099), V(0.776732950940391, -0.43614824442447003, -0.45437939192802856))
	let line = new L3(V(-1560.8950828838565, 716.07295580975, 249.61382611323648), V(0.9130103135570956, -0.36545647611595106, -0.18125598308272678))
	let face2 = new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 100)), [
		new PCurveEdge(new EllipseCurve(V(73.03583314037537, -69.86032483338774, -100), V(-24.176861672352132, -146.16681457389276, 0), V(146.16681457389276, -24.176861672352132, 0)), V(219.75148278474705, -90.44615667066816, -100), V(97.2126948127275, 76.30648974050503, -100), 1.5953170840348225, -3.141592653589793, null, V(-20.58583183728038, -146.71564964437164, 0), V(146.16681457389276, -24.176861672352114, 0)),
		StraightEdge.throughPoints(V(97.2126948127275, 76.30648974050503, -100), V(275.99999999999966, 255.99999999999972, -100)),
		new PCurveEdge(new BezierCurve(V(219.75148278474705, -90.44615667066816, -100), V(-82.00000000000018, -138.00000000000023, -100), V(539.9999999999997, 225.9999999999997, -100), V(275.99999999999966, 255.99999999999972, -100), -0.1, 1.1), V(275.99999999999966, 255.99999999999972, -100), V(219.75148278474705, -90.44615667066816, -100), 1, 0, null, V(792, -90.00000000000009, 0), V(905.2544483542417, 142.6615299879962, 0))], [])
	console.log(face.intersectsLine(line), face.surface.isTsForLine(line), face.surface.isTsForLine(line))
	const t = line.intersectWithPlaneLambda(face2.surface.plane);
	console.log(face2.intersectsLine(line), t, line.at(t).sce)
	assert.ok(face.intersectsLine(line))
});
QUnit.testDifferentSystems( "Matrix4x4 eigenValues and eigenVectors", function (assert, /** M4*/ m4) {
	let eigenValues = m4.realEigenValues3()
	console.log(eigenValues)
//		assert.equal(eigenValues.length, 3)
	eigenValues.forEach((eigenValue, i)=> {
		assert.ok(NLA.eq0(M4.IDENTITY.timesScalar(-eigenValue).plus(m4.as3x3()).determinant()))
		//assert.ok(ei)
	})
	let eigenVectors = m4.realEigenVectors3()
	console.log(eigenVectors)
	eigenVectors.forEach((eigenVector, i) => {
		m4.isNormal() && assert.ok(eigenVector.isPerpendicularTo(eigenVectors[(i+1)%eigenVectors.length]), `eigenVector${i}.isPerpendicularTo(eigenVector${(i+1)%eigenVectors.length})`)
		assert.ok(!eigenVector.isZero(), `'!eigenVector${i}.isZero()` + !eigenVector.isZero())
		assert.ok(eigenVector.isParallelTo(m4.transformVector(eigenVector)), `eigenVector${i}.isParallelTo(m4.transformVector(eigenVector))`)
	})

}, new M4(
	3,2,4,0,
	2,0,2,0,
	4,2,3,0,
	0,0,0,1));

QUnit.testDifferentSystems( "EllipseCurve.isPointsWithBezier()", function (assert, /** M4*/ m4) {
	let ell = new EllipseCurve(V(-223.34900663163222, -176.63214006755936, 0), V(-169.5891804980124, -35.54247345835796, 0), V(35.54247345835796, -169.5891804980124, 0))
	let bez = new BezierCurve(V(-267.6481190901419, -368.37017217006473, 0), V(563.959389388763, 94.96018577817034, 0), V(-1110.7787051488917, -95.8394860073627, 0), V(-59.14331799274822, -299.7830665459221, 0))
	let isPoints = ell.isPointsWithBezier(bez)
	assert.equal(isPoints.length, 6)
	console.log(isPoints.join())
	isPoints.forEach(p => {
		assert.ok(ell.containsPoint(p), `ell.distanceToPoint(${p}) = ${ell.distanceToPoint(p)}`)
		assert.ok(bez.containsPoint(p), `bez.distanceToPoint(${p}) = ${bez.distanceToPoint(p)}`)
	})
})

QUnit.test( "BezierCurve.isPointsWithBezier()", function (assert, /** M4*/ m4) {
	let curve1 = BezierCurve.graphXY(2,-3,-3,2, -2, 3)
	let curve2 = curve1.transform(M4.rotationLine(V(0.5, 0), V3.Z, PI/2))
	let isInfos = curve1.isInfosWithBezier(curve2)
	assert.equal(isInfos.length, 9)
	console.log(isInfos.map(SCE))
	isInfos.forEach(info => {
		let p = info.p
		assert.ok(curve1.containsPoint(p), `curve1.distanceToPoint(${p}) = ${curve1.distanceToPoint(p, -2, 3)}`)
		assert.ok(curve2.containsPoint(p), `curve2.distanceToPoint(${p}) = ${curve2.distanceToPoint(p, -2, 3)}`)
	})

	// test self-intersections
	;[new BezierCurve(V(133, 205, 0), V(33, 240, 0), V(63, 168, 0), V(151, 231, 0)),
		new BezierCurve(V3.ZERO, V(10, 10), V(-9, 10), V3.X)].forEach(curve => {
		assert.push(true, undefined, undefined, 'Testing ' + curve.sce)
		let isInfos = curve.isInfosWithBezier(curve, 0,1,0,1)
		assert.equal(isInfos.length, 1)
		console.log(isInfos.map(SCE))
		isInfos.forEach(info => {
			let p = info.p
			assert.ok(NLA.eq0(curve.distanceToPoint(p)), `curve.distanceToPoint(${p}) = ${curve.distanceToPoint(p, -2, 3)}`)
		})
	})
})
QUnit.test( "EllipseCurve.isInfosWithBezier2D()", function (assert, /** M4*/ m4) {
	let ell = EllipseCurve.forAB(3, 1)
	let bez = BezierCurve.graphXY(2,-3,-3,2, -2, 3)
	let isInfos = ell.isInfosWithBezier2D(bez)
	assert.equal(isInfos.length, 6)
	console.log(isInfos.map(SCE))
	isInfos.forEach(info => {
		let p = info.p
		assert.ok(ell.containsPoint(p), `curve1.distanceToPoint(${p}) = ${ell.distanceToPoint(p)}`)
		assert.ok(bez.containsPoint(p), `curve2.distanceToPoint(${p}) = ${bez.distanceToPoint(p, -2, 3)}`)
	})

})

QUnit.test( 'Edge.edgesIntersects' , function (assert) {
	let curve1 = BezierCurve.graphXY(2,-3,-3,2, -2, 3)
	let curve2 = curve1.transform(M4.rotationLine(V(0.5, 0), V3.Z, PI/2))
	let edge1 = PCurveEdge.forCurveAndTs(curve1, 0, 1)
	let edge2 = PCurveEdge.forCurveAndTs(curve2, 0, 1)
	assert.ok(Edge.edgesIntersect(edge1, edge2))
	assert.notOk(Edge.edgesIntersect(edge1, edge1.translate(10, 0, 0)))
	assert.notOk(Edge.edgesIntersect(edge1, edge2.translate(10, 0, 0)))


})
QUnit.testDifferentSystems( 'BezierCurve.isTsWithSurface(CylinderSurface)' , function (assert, m4) {
	let bez = BezierCurve.graphXY(2,-3,-3,2, -2, 3).rotateX(15 * DEG).translate(0, 0, 100).transform(m4)
	let cyl = new CylinderSurface(EllipseCurve.forAB(4, 1).rotateY(10 * DEG), V3.Z).transform(m4)
	let ts = bez.isTsWithSurface(cyl)
	assert.equal(6, ts.length)
	ts.forEach(t => {
		const p = bez.at(t)
		assert.ok(cyl.containsPoint(p), cyl.implicitFunction()(p))
	})
})

QUnit.test( 'BezierCurve.pointLambda' , function (assert) {
	let curve = new BezierCurve(
		V( 92.48132002394416, 253.35277539335377, 0),
		V( 99.18055157018783, 225.4322156490681, 0),
		V( 63.52151476563836, 168.59279980361327, 0),
		V(151.89049176954802, 231.21343792479922, 0))
	let p = V(90.8280915025532, 214.7764313721318, 0)
	assert.ok(NLA.eq0(curve.distanceToPoint(p)))
	assert.ok(isFinite(curve.pointLambda(p)))


})

QUnit.testDifferentSystems('EllipseCurve.getAreaInDir', function (assert, m4) {
		let k = 1;
		[
			{right: V3.X, up: V3.Y, s: 0, t: PI, result: PI / 2, c: V(0, 4 / 3 / PI)},
			{right: V3.X, up: V3.Y, s: PI, t: 0, result: -PI / 2, c: V(0, 4 / 3 / PI)},
			{right: V3.X, up: V3.Y, s: -PI / 2, t: PI / 2, result: PI / 2, c: V(4 / 3 / PI, 0)},
			{right: V3.X, up: V3.Y, s: -PI, t: 0, result: PI / 2, c: V(0, -4 / 3 / PI)},
			// let "down" be X
			{right: V3.Y, up: V3.X.negated(), s: 0, t: PI, result: PI / 2, c: V(0, 4 / 3 / PI)},
			{right: V3.Y, up: V3.X.negated(), s: -PI, t: 0, result: PI / 2, c: V(0, -4 / 3 / PI)}
		].forEach(test => {
			[0, 4].forEach(yDiff => {
				let r = m4.transformVector(test.right)
				let areaFactor = m4.transformVector(V3.X).cross(m4.transformVector(V3.Y)).length()
				console.log(areaFactor)
				const ell = EllipseCurve.UNIT.translate(0, yDiff, 0).transform(m4)
				let up = m4.transformVector(test.up).normalized()
				let offsetArea = yDiff * ((1 - cos(test.t)) - (1 - cos(test.s))) * test.up.dot(V3.Y)
				const totalArea = test.result + offsetArea
				const expectedArea = totalArea * areaFactor
				let result = ell.getAreaInDir(r, up, test.s, test.t)
				let offsetCentroid = V((cos(test.t) + cos(test.s)) / 2, yDiff / 2)
				let movedCentroid = test.c.plus(V(0, yDiff))
				const expectedCentroid = m4.transformPoint(movedCentroid.times(test.result).plus(offsetCentroid.times(offsetArea)).div(totalArea))
				console.log(test.t, test.s, 1 - cos(test.t), 1 - cos(test.s))
				console.log(test.c.times(test.result).str, offsetCentroid.str, offsetArea, offsetCentroid.times(offsetArea).str, test.c.times(test.result).plus(offsetCentroid.times(offsetArea)).str, totalArea, expectedCentroid.str   )
				assert.fuzzyEquals(
					result.area,
					expectedArea,
					`yDiff: ${yDiff}, ${test.sce}, offsetArea: ${offsetArea}, expected: ${expectedArea}, ${areaFactor * offsetArea}`)

				assert.fuzzyEquals(result.centroid.x, expectedCentroid.x, 'cx ' + result.centroid.x)
				assert.fuzzyEquals(result.centroid.y, expectedCentroid.y, 'cy ' + result.centroid.y)
				// if (!k--) throw new Error()
			})
		})
	},
	M4.IDENTITY,
	M4.rotationZ(45 * DEG),
	M4.forRows(V(1, 2, 3), V3.Y, V3.Z),
	M4.FOO.as3x3())

QUnit.test( 'M4.gauss' , function (assert) {
	let m = new M4(
		0, 0, 1, -1,
		2, -3, -3, 4,
		0, 0, 0, 0,
		0, 0, 0, 1
	)
	let mGauss = new M4(
		2, -3, -3, 4,
		0, 0, 1, -1,
		0, 0, 0, 1,
		0, 0, 0, 0
	)
	assert.matrixEquals(m.gauss().U, mGauss)

	console.log(M4.FOO.gauss().U.str)
})

QUnit.test('asd', function (assert) {
	const surface = new CylinderSurface(new EllipseCurve(V(0, 0, 0), V(8, 0, 0), V(0, -8, 0)), V(0, 0, -1))
	const loop = [
		StraightEdge.throughPoints(V(1, 7.937253933193773, 4), V(1, 7.937253933193773, 1)),
		new PCurveEdge(new EllipseCurve(V(0, 0, 1), V(8, 0, 0), V(0, -8, 0)), V(1, 7.937253933193773, 1), V(6, 5.291502622129181, 1), -1.4454684956268313, -0.7227342478134156, null, V(7.937253933193772, -0.9999999999999991, 0), V(5.2915026221291805, -6, 0)),
		StraightEdge.throughPoints(V(6, 5.291502622129181, 1), V(6, 5.291502622129181, 4)),
		new PCurveEdge(new EllipseCurve(V(0, 0, 4), V(8, 0, 0), V(0, 8, 0)), V(6, 5.291502622129181, 4), V(1, 7.937253933193773, 4), 0.7227342478134156, 1.4454684956268313, null, V(-5.2915026221291805, 6, 0), V(-7.937253933193772, 0.9999999999999991, 0))
	]

	assert.equal(surface.edgeLoopContainsPoint(loop, V(8, 0, 0)), PointVsFace.OUTSIDE)
	assert.equal(surface.edgeLoopContainsPoint(loop, V(1, 7.937253933193773, 3)), PointVsFace.ON_EDGE)
})