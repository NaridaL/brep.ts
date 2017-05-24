//window.onerror = function (errorMsg, url, lineNumber, column, errorObj) {
//	console.log(errorMsg, url, lineNumber, column, errorObj)
//}
QUnit.assert.matrixEquals = function (actual, expected, message, precision) {
	this.push(actual.equalsMatrix(expected, precision), actual.toString(), expected.toString(), message)
}
QUnit.assert.matrixEquivalent = function (actual, expected, message, precision) {
	this.push(actual.normalized2().equalsMatrix(expected.normalized2(), precision), actual.toString(), expected.toString(), message)
}
QUnit.assert.V3equals = function (actual, expected, message, precision) {
	this.push(false, actual.toString(), expected.toString(), message)
}
QUnit.assert.V3like = function (actual, expected, message, precision) {
	this.push(expected.like(actual), actual.toString(), expected.toString(), (message ? message : 'V3like') + '; |dv| = ' + expected.distanceTo(actual))
}

QUnit.module('NLA')
function testDifferentSystems(name: string, test: (assert: Assert, m4: M4) => void, ...what: M4[]) {
	if (!what.length) {
		what = [M4.IDENTITY, M4.FOO]
	}
    what.forEach((m, i) => {
        QUnit.test(name + '_' + i, assert => {
            console.log(`TESTING '${name}' WITH '${m.toString()}`)
            assert.push(true, undefined, undefined, `<html>'${name}' WITH <pre style="font-weight: bold">${m.toString()}</pre>`)
            test(assert, m)
        })
    })
}

testDifferentSystems('Matrix4x4 eigenValues and eigenVectors', function (assert, m4: M4) {
	const eigenValues = m4.realEigenValues3()
	console.log(eigenValues)
//		assert.equal(eigenValues.length, 3)
	eigenValues.forEach((eigenValue, i)=> {
		assert.ok(eq0(M4.IDENTITY.scale(-eigenValue).plus(m4).determinant()))
		//assert.ok(ei)
	})
	const eigenVectors = m4.realEigenVectors3()
	console.log(eigenVectors)
	eigenVectors.forEach((eigenVector, i) => {
		m4.isSymmetric() && assert.ok(eigenVector.isPerpendicularTo(eigenVectors[(i + 1) % eigenVectors.length]), `eigenVector${i}.isPerpendicularTo(eigenVector${(i + 1) % eigenVectors.length})`)
		assert.ok(!eigenVector.likeO(), `'!eigenVector${i}.isZero()` + !eigenVector.likeO())
		assert.ok(eigenVector.isParallelTo(m4.transformVector(eigenVector)), `eigenVector${i}.isParallelTo(m4.transformVector(eigenVector))`)
	})

}, new M4(
	3, 2, 4, 0,
	2, 0, 2, 0,
	4, 2, 3, 0,
	0, 0, 0, 1), M4.FOO.as3x3())

testDifferentSystems('SemiEllipseCurve.isPointsWithBezier()', function (assert, m4: M4) {
	const ell = new SemiEllipseCurve(V(-223.34900663163222, -176.63214006755936, 0), V(-169.5891804980124, -35.54247345835796, 0), V(35.54247345835796, -169.5891804980124, 0))
    const bez = new BezierCurve(V(-267.6481190901419, -368.37017217006473, 0), V(563.959389388763, 94.96018577817034, 0), V(-1110.7787051488917, -95.8394860073627, 0), V(-59.14331799274822, -299.7830665459221, 0))
    const isPoints = ell.isInfosWithBezier(bez).map(info => info.p)
	assert.equal(isPoints.length, 3)
	isPoints.forEach(p => {
		assert.ok(ell.containsPoint(p), `ell.distanceToPoint(${p}) = ${ell.distanceToPoint(p)}`)
		assert.ok(bez.containsPoint(p), `bez.distanceToPoint(${p}) = ${bez.distanceToPoint(p)}`)
	})
})

testDifferentSystems('M4.svd3()', function (assert, m4: M4) {

	const {U, SIGMA, VSTAR} = m4.svd3()
	const U_SIGMA_VSTAR = M4.multiplyMultiple(U, SIGMA, VSTAR)
	assert.ok(SIGMA.isDiagonal(), 'SIGMA.isDiagonal()')
	assert.ok(U.isOrthogonal(), 'U.isOrthogonal()\n' + U.str + '\n' + U.times(U.transposed()).str)
	assert.ok(VSTAR.isOrthogonal(), 'VSTAR.isOrthogonal()\n' + VSTAR.str + '\n' + VSTAR.times(VSTAR.transposed()).str)
	assert.push(U_SIGMA_VSTAR.likeM4(m4.as3x3()), U_SIGMA_VSTAR.str, m4.as3x3().str, 'U_SIGMA_VSTAR == A')
})

//testDifferentSystems('SemiCylinderSurface.calculateArea', function (assert, m4) {
//	const surface = SemiCylinderSurface.UNIT.transform(m4)
//	// loop which is 1 high and goes around a quarter of the cylinder
//	const loop = [
//		StraightEdge.throughPoints(V(1, 0, 1), V(1, 0, 0)),
//		Edge.forCurveAndTs(SemiEllipseCurve.UNIT, 0, PI / 2),
//		StraightEdge.throughPoints(V(0, 1, 0), V(0, 1, 1)),
//		Edge.forCurveAndTs(SemiEllipseCurve.UNIT.translate(0, 0, 1), PI / 2, 0)].map(edge => edge.transform(m4))
//	const face = Face.create(surface, loop)
//	linkB2(assert, `mesh=${face.sce}.scale(100, 100, 100).toMesh()`)
//	const area = face.calcArea()
//	if (m4.isOrthogonal()) {
//		assert.push(eq(area, PI/2), area, PI / 2)
//	} else {
//		const expectedArea = face.toMesh().calcVolume().area
//		assert.push(eq2(area, expectedArea, 0.1), area, expectedArea)
//	}
//
//
//	const loopReverse = Edge.reversePath(loop)
//	const holeArea = surface.calculateArea(loopReverse)
//	if (m4.isOrthogonal()) {
//		assert.push(eq(holeArea, -PI/2), area, -PI / 2)
//	} else {
//		const expectedArea = face.toMesh().calcVolume().area
//		assert.push(eq2(holeArea, -expectedArea, 0.1), holeArea, -expectedArea)
//	}
//
//	const flippedSurfaceArea = surface.flipped().calculateArea(loop)
//	if (m4.isOrthogonal()) {
//		assert.push(eq(holeArea, -PI/2), area, -PI / 2)
//	} else {
//		const expectedArea = face.toMesh().calcVolume().area
//		assert.push(eq2(holeArea, -expectedArea, 0.1), holeArea, -expectedArea)
//	}
//
//
//	{
//		const loop = [
//			StraightEdge.throughPoints(V(1, 0, 1), V(1, 0, 0)),
//			Edge.forCurveAndTs(SemiEllipseCurve.forAB(1, 1), 0, PI/2),
//			Edge.forCurveAndTs(new SemiEllipseCurve(V3.O, V(1, 0, 1), V(0, 1, 0)), PI / 2, 0)].map(edge => edge.transform(m4))
//		const face = Face.create(surface, loop)
//		linkB2(assert, `mesh=${face.sce}.scale(100, 100, 100).toMesh()`)
//		const area = face.calcArea()
//		if (m4.isOrthogonal()) {
//			assert.push(eq(area, 1), area, 1)
//		} else {
//			const expectedArea = face.toMesh().calcVolume().area
//			assert.push(eq2(area, expectedArea, 0.1), area, expectedArea)
//		}
//
//
//		const loopReverse = Edge.reversePath(loop)
//		const holeArea = surface.calculateArea(loopReverse)
//		if (m4.isOrthogonal()) {
//			assert.push(eq(holeArea, -1), area, -1)
//		} else {
//			const expectedArea = face.toMesh().calcVolume().area
//			console.log('expectedArea', expectedArea)
//			console.log('expectedArea', eval(face.sce).toMesh().calcVolume().area)
//			assert.push(eq2(holeArea, -expectedArea, 0.1), holeArea, -expectedArea)
//		}
//	}
//}, M4.IDENTITY, M4.translate(0, 0, 4).times(M4.rotateY(10 * DEG)), M4.scale(1, 2, 1))
//
//testDifferentSystems('EllipsoidSurface.calculateArea', function (assert, m4) {
//	const surface = SemiEllipsoidSurface.UNIT.transform(m4)
//	// loop which is 1 high and goes around a quarter of the cylinder
//	const loop = [
//		Edge.forCurveAndTs(SemiEllipseCurve.UNIT, 10 * DEG, 40 * DEG),
//		Edge.forCurveAndTs(new SemiEllipseCurve(V3.O, V3.sphere(40 * DEG, 0), V3.Z), 0, PI / 2),
//		Edge.forCurveAndTs(new SemiEllipseCurve(V3.O, V3.sphere(10 * DEG, 0), V3.Z), PI / 2, 0)].map(edge => edge.transform(m4))
//	const face = Face.create(surface, loop)
//	assert.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
//						href='brep2.html?mesh=${face.sce}.scale(100, 100, 100).toMesh()'>view</a>`)
//	const area = face.calcArea()
//	const expectedArea = face.toMesh().calcVolume().area
//	assert.push(eq2(area, expectedArea, 0.1), area, expectedArea)
//
//
//	//const loopReverse = Edge.reversePath(loop)
//	//const holeArea = surface.calculateArea(loopReverse)
//	//if (m4.isOrthogonal()) {
//	//	assert.push(eq(holeArea, -PI/2), area, -PI / 2)
//	//} else {
//	//	const expectedArea = face.toMesh().calcVolume().area
//	//	assert.push(eq2(holeArea, -expectedArea, 0.1), holeArea, -expectedArea)
//	//}
//	//
//	//const flippedSurfaceArea = surface.flipped().calculateArea(loop)
//	//if (m4.isOrthogonal()) {
//	//	assert.push(eq(holeArea, -PI/2), area, -PI / 2)
//	//} else {
//	//	const expectedArea = face.toMesh().calcVolume().area
//	//	assert.push(eq2(holeArea, -expectedArea, 0.1), holeArea, -expectedArea)
//	//}
//	//
//	//
//	//{
//	//	const loop = [
//	//		StraightEdge.throughPoints(V(1, 0, 1), V(1, 0, 0)),
//	//		Edge.forCurveAndTs(SemiEllipseCurve.forAB(1, -1), 0, -PI / 2),
//	//		Edge.forCurveAndTs(new SemiEllipseCurve(V3.O, V(0, 1, 0), V(1, 0, 1)), 0, PI / 2)].map(edge => edge.transform(m4))
//	//	const face = Face.create(surface, loop)
//	//	assert.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
//	//					href='viewer.html?mesh=${face.sce}.scale(100, 100, 100).toMesh()'>view</a>`)
//	//	const area = face.calcArea()
//	//	if (m4.isOrthogonal()) {
//	//		assert.push(eq(area, 1), area, 1)
//	//	} else {
//	//		const expectedArea = face.toMesh().calcVolume().area
//	//		assert.push(eq2(area, expectedArea, 0.1), area, expectedArea)
//	//	}
//	//
//	//
//	//	const loopReverse = Edge.reversePath(loop)
//	//	const holeArea = surface.calculateArea(loopReverse)
//	//	if (m4.isOrthogonal()) {
//	//		assert.push(eq(holeArea, -1), area, -1)
//	//	} else {
//	//		const expectedArea = face.toMesh().calcVolume().area
//	//		console.log("expectedArea", expectedArea)
//	//		console.log("expectedArea", eval(face.sce).toMesh().calcVolume().area)
//	//		assert.push(eq2(holeArea, -expectedArea, 0.1), holeArea, -expectedArea)
//	//	}
//	//}
//}, M4.IDENTITY, M4.translate(0, 0, 4).times(M4.rotateY(10 * DEG)), M4.scale(1, 1, 2))

testDifferentSystems('SemiEllipseCurve.getAreaInDir', function (assert, m4) {
		const k = 1;
		[
			{right: V3.X, up: V3.Y, s: 0, t: PI, result: PI / 2, c: V(0, 4 / 3 / PI)},
			{right: V3.X, up: V3.Y, s: PI, t: 0, result: -PI / 2, c: V(0, 4 / 3 / PI)},
			{right: V3.X, up: V3.Y, s: -PI / 2, t: PI / 2, result: PI / 2, c: V(4 / 3 / PI, 0)},
			{right: V3.X, up: V3.Y, s: -PI, t: 0, result: PI / 2, c: V(0, -4 / 3 / PI)},
			// let 'down' be X
			{right: V3.Y, up: V3.X.negated(), s: 0, t: PI, result: PI / 2, c: V(0, 4 / 3 / PI)},
			{right: V3.Y, up: V3.X.negated(), s: -PI, t: 0, result: PI / 2, c: V(0, -4 / 3 / PI)}
		].forEach(test => {
			[0, 4].forEach(yDiff => {
				const r = m4.transformVector(test.right)
				const areaFactor = m4.transformVector(V3.X).cross(m4.transformVector(V3.Y)).length()
				console.log(areaFactor)
				const ell = SemiEllipseCurve.UNIT.translate(0, yDiff, 0).transform(m4)
				const up = m4.transformVector(test.up).unit()
				const offsetArea = yDiff * ((1 - cos(test.t)) - (1 - cos(test.s))) * test.up.dot(V3.Y)
				const totalArea = test.result + offsetArea
				const expectedArea = totalArea * areaFactor
				const result = ell.getAreaInDir(r, up, test.s, test.t)
				const offsetCentroid = V((cos(test.t) + cos(test.s)) / 2, yDiff / 2)
				const movedCentroid = test.c.plus(V(0, yDiff))
				const expectedCentroid = m4.transformPoint(movedCentroid.times(test.result).plus(offsetCentroid.times(offsetArea)).div(totalArea))
				console.log(test.t, test.s, 1 - cos(test.t), 1 - cos(test.s))
				console.log(test.c.times(test.result).str, offsetCentroid.str, offsetArea, offsetCentroid.times(offsetArea).str, test.c.times(test.result).plus(offsetCentroid.times(offsetArea)).str, totalArea, expectedCentroid.str)
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
	M4.rotateZ(45 * DEG),
	M4.forRows(V(1, 2, 3), V3.Y, V3.Z),
	M4.FOO.as3x3())

registerTests({

	'Edge.edgesIntersects'(assert) {
		const curve1 = BezierCurve.graphXY(2, -3, -3, 2, -2, 3)
		const curve2 = curve1.transform(M4.rotateLine(V(0.5, 0), V3.Z, PI / 2))
		const edge1 = PCurveEdge.forCurveAndTs(curve1, 0, 1)
		const edge2 = PCurveEdge.forCurveAndTs(curve2, 0, 1)
		assert.ok(Edge.edgesIntersect(edge1, edge2))
		assert.notOk(Edge.edgesIntersect(edge1, edge1.translate(10, 0, 0)))
		assert.notOk(Edge.edgesIntersect(edge1, edge2.translate(10, 0, 0)))
	},
	'eqAngle'(assert) {
		assert.ok(zeroAngle(0))
		assert.ok(zeroAngle(-NLA_PRECISION / 2))
		assert.ok(zeroAngle(2 * Math.PI - NLA_PRECISION / 2))
		assert.ok(zeroAngle(2 * Math.PI + NLA_PRECISION / 2))
		assert.ok(eqAngle(-Math.PI, Math.PI))
		assert.ok(eqAngle(0, 2 * Math.PI))
		assert.ok(eqAngle(0, 2 * Math.PI - NLA_PRECISION / 2))
		assert.ok(eqAngle(0, 2 * Math.PI + NLA_PRECISION / 2))
		assert.notOk(eqAngle(-Math.PI, 2 * Math.PI))
		assert.notOk(eqAngle(0, Math.PI))
	},
	'eq etc'(assert) {
		assert.notOk(lt(2, 2 + NLA_PRECISION / 2))
		assert.notOk(lt(2, 2 - NLA_PRECISION / 2))
		assert.ok(le(2, 2 + NLA_PRECISION / 2))
		assert.ok(le(2, 2 - NLA_PRECISION / 2))

		assert.notOk(gt(2, 2 + NLA_PRECISION / 2))
		assert.notOk(gt(2, 2 - NLA_PRECISION / 2))
		assert.ok(ge(2, 2 + NLA_PRECISION / 2))
		assert.ok(ge(2, 2 - NLA_PRECISION / 2))

		assert.ok(lt(2, 3))
		assert.ok(gt(3, 2))
		assert.ok(le(2, 3))
		assert.ok(ge(3, 2))

		assert.ok(eq(2, 2 + NLA_PRECISION) == !gt(2, 2 + NLA_PRECISION))
		assert.ok(eq(2, 2 + NLA_PRECISION) == ge(2, 2 + NLA_PRECISION))
	},
	'arrayCopyStep'(assert) {
		const a = [1, 2, 3, 4, 5, 6, 7, 8]
		const b = [-1, -2, -3, -4]
		arrayCopyStep(b, 0, 1, a, 1, 2, 3)
		assert.deepEqual(a, [1, -1, 3, -2, 5, -3, 7, 8])
	},
	'arrayCopy'(assert) {
		const a = [1, 2, 3, 4, 5, 6, 7, 8]
		const b = [-1, -2, -3, -4]
		arrayCopy(b, 1, a, 2, 2)
		assert.deepEqual(a, [1, 2, -2, -3, 5, 6, 7, 8])
	},
	'LU Decomposition'(assert) {
		const a = Matrix.fromRowArrays([6, 3], [4, 3])
		const aLU = a.luDecomposition()
		assert.deepEqual(aLU.L.asRowArrays(Array), [[1, 0], [4 / 6, 1]])
		assert.deepEqual(aLU.U.asRowArrays(Array), [[6, 3], [0, 1]])
		assert.deepEqual(aLU.P.asRowArrays(Array), [[1, 0], [0, 1]])
		assert.matrixEquals(aLU.P.times(a), aLU.L.times(aLU.U))
	},
	'LU Decomposition 2'(assert) {
		const a = Matrix.fromFunction(8, 8, (i, j) => Math.round((Math.random() - 0.5) * 4096))
		//var a = Matrix.fromRowArrays2([[-1636, 1740, -516], [-708, 403, 1986], [-1256, -1493, 996]])
		assert.ok(a.rowsIndependent())
		const aLU = a.luDecomposition()
		assert.ok(aLU.P.isPermutation())
		assert.ok(aLU.L.isLowerTriangular())
		assert.ok(aLU.U.isUpperTriangular())
		assert.matrixEquals(aLU.P.times(a), aLU.L.times(aLU.U))
	},
	'LU Decomposition 3'(assert) {
		const a = Matrix.fromRowArrays([0, 1, 1], [1, 1, 1], [1, 2, 3])
		const aLU = a.luDecomposition()
		assert.deepEqual(aLU.U.asRowArrays(Array), [[1, 1, 1], [0, 1, 1], [0, 0, 1]])
		assert.deepEqual(aLU.L.asRowArrays(Array), [[1, 0, 0], [0, 1, 0], [1, 1, 1]])
		assert.deepEqual(aLU.P.asRowArrays(Array), [[0, 1, 0], [1, 0, 0], [0, 0, 1]])
		assert.matrixEquals(aLU.P.times(a), aLU.L.times(aLU.U))
	},

	'Plane3.prototype.projectedVector'(assert) {
		const p = new P3(V(0, 0, 1), 2)
		assert.ok(V(1, 1, 0).like(p.projectedVector(V(1, 1, 1))))
	},
	'Line3.prototype.distanceToLine'(assert) {
		assert.equal(L3.X.distanceToLine(new L3(V3.Z, V3.Y)), 1)
		assert.equal(L3.X.distanceToLine(new L3(V3.Z, V3.X)), 1)
	},
	'Plane3.prototype.transformed'(assert) {
		const p = new P3(V(0, 0, 1), 2)
		assert.ok(P3.XY.like(P3.XY.transform(M4.identity())))
	},
	'Plane3.prototype.intersectionWithPlane'(assert) {
		assert.ok(P3.XY.intersectionWithPlane(P3.ZX).isColinearTo(L3.X))
		assert.ok(P3.ZX.intersectionWithPlane(P3.XY).isColinearTo(L3.X))
		assert.notOk(P3.ZX.intersectionWithPlane(P3.XY).isColinearTo(L3.Y))
	},
	'Line3.prototype.isTsForLine'(assert) {
		console.log(L3.X.isInfoWithLine(new L3(V(1, 1, 0), V3.Y)).sce)
		assert.ok(L3.X.isInfoWithLine(new L3(V(1, 1, 0), V3.Y)).equals(V3.X))
	},
	'V3.prototype.zip'(assert) {
		const a = V(1, 2, 3), b = V(4, 5, 6)
		assert.ok(V3.zip((a, b) => a + 3 * b, a, b).equals(a.plus(b.times(3))))
	},
	'magic'(assert) {
		assert.expect(0)
		const a = V(1, 2, 3), b = V(4, 5, 6)
		//assert.ok(magic('a b c s => abs(a) x b .. c + 3 + ')(a, b, c, 3).equals(a.abs().cross(b).dot(c) + 3))
	},
	'AABB'(assert) {
		const a = new AABB(V3.O, V(20, 10, 30))
		const b = new AABB(V3.O, V(1, 1, 1))
		assert.ok(a.likeAABB(a))
		assert.notOk(a.likeAABB(a.translate(10, 0, 0)))
		assert.ok(a.withoutAABB(b).likeAABB(new AABB(V(0, 0, 1), V(20, 10, 30))))
	},
// QUnit.test( 'V3.areDisjoint', function(assert ) {
// 	assert.ok(V3.areDisjoint([V3.O, V3.X, V3.Y]))
// 	assert.ok(V3.areDisjoint([V3.O, V3.X, V3.X, V3.Y])) // same value twice is ok, as same reference
// 	assert.notOk(V3.areDisjoint([V3.O, V3.X, V(0, 0, 0), V3.Y])) // not ok as V3.O != V(0, 0, 0)
// 	assert.notOk(V3.areDisjoint([V3.O, V3.X, V(NLA_PRECISION / 2, 0, 0), V3.Y])) // not ok as !V3.O.like(V(...))
// 	assert.ok(V3.areDisjoint([V(NLA_PRECISION * -0.7, 0, 0), V(NLA_PRECISION * 0.7, 0, 0)]))
// })
	'V3.areDisjoint2'(assert) {
		console.log(~~2147483657)
		const s = new CustomSet()
		const a = V(0, 2.7499999999999996, -5), b = V(0, 2.749999999999999, -5)
		s.canonicalizeLike(a)
		console.log(s._map, a.like(b), a.hashCodes(), b.hashCodes(), a.hashCode(), b.hashCode())
		assert.ok(s.canonicalizeLike(b) == a)
	},
	'arrayBinaryInsert'(assert) {
		const arr = [1, 2, 3, 4]
		arr.binaryInsert(2.5, (a, b) => a - b)
		assert.deepEqual(arr, [1, 2, 2.5, 3, 4])

		const arr2 = []
		arr2.binaryInsert(-2)
		arr2.binaryInsert(5)
		assert.deepEqual(arr2, [-2, 5])
	},
	'arrayBinaryIndexOf'(assert) {
		assert.equal([1, 2, 2.5, 3, 4].binaryIndexOf(3, (a, b) => a - b), 3)
		assert.equal([1, 2, 2.5, 3, 4].binaryIndexOf(2.6, (a, b) => a - b), -3 - 1)
	},
	'newtonIterate2d'(assert) {
		const res = newtonIterate2d((s, t) => s - 2, (s, t) => t - 4, 5, 5)
		assert.push(res.like(V(2, 4)), res.sce, V(2, 4).sce)
	},

	//'SemiEllipseCurve.getVolZAnd'(assert) {
	//
	//	assert.equal(SemiEllipseCurve.UNIT.getVolZAnd(V3.Z, -PI, PI).volume, 0)
	//	assert.equal(SemiEllipseCurve.UNIT.rotateY(90 * DEG).translate(1, 0, 0).getVolZAnd(V3.Z, -PI, PI).volume, PI)
	//},
    'solveCubicReal2()'(assert) {
        assert.deepEqual(solveCubicReal2(0, 1, 0, 0), [0])
        assert.deepEqual(solveCubicReal2(1, 0, 0, 0), [0])
    },
    'solveCubicReal2() 2'(assert) {
        const [a, b, c, d] = [-2.38, -11.88, +30.9, -16.64]
        const results = solveCubicReal2(a,b,c,d)
        assert.equal(results.length, 2)
        results.forEach(x => {
            const y = a * x ** 3 + b * x ** 2 + c * x + d
            assert.push(eq0(y), y, 0)
        })
    },
	'intersectionUnitCircleLine'(assert) {
		// y = -x + 1 => x + y = 1
		assert.deepEqual(intersectionUnitCircleLine(1, 1, 1), {x1: 1, x2: 0, y1: 0, y2: 1})
	},
	'intersectionCircleLine'(assert) {
		// y = -x + 2 => x + y = 2
		assert.deepEqual(intersectionCircleLine(1, 1, 2, 2), {x1: 2, x2: 0, y1: 0, y2: 2})
	},
	//'SemiCylinderSurface Face containsPoint'(assert) {
	//	let face = new RotationFace(new SemiCylinderSurface(new SemiEllipseCurve(V(73.03583314037537, -69.86032483338774, 0), V(-24.176861672352132, -146.16681457389276, 0), V(146.16681457389276, -24.176861672352132, 0)), V(0, 0, 1)), [
	//		new PCurveEdge(new SemiEllipseCurve(V(73.03583314037537, -69.86032483338774, 0), V(-24.176861672352132, -146.16681457389276, 0), V(146.16681457389276, -24.176861672352132, 0)), V(219.75148278474705, -90.44615667066816, 0), V(97.2126948127275, 76.30648974050503, 0), 1.5953170840348225, -3.141592653589793, null, V(-20.58583183728038, -146.71564964437164, 0), V(146.16681457389276, -24.176861672352114, 0)),
	//		StraightEdge.throughPoints(V(97.2126948127275, 76.30648974050503, 0), V(97.2126948127275, 76.30648974050503, -100)),
	//		new PCurveEdge(new SemiEllipseCurve(V(73.03583314037537, -69.86032483338774, -100), V(-24.176861672352132, -146.16681457389276, 0), V(146.16681457389276, -24.176861672352132, 0)), V(97.2126948127275, 76.30648974050503, -100), V(219.75148278474705, -90.44615667066816, -100), -3.141592653589793, 1.5953170840348225, null, V(-146.16681457389276, 24.176861672352114, 0), V(20.58583183728038, 146.71564964437164, 0)),
	//		StraightEdge.throughPoints(V(219.75148278474705, -90.44615667066816, -100), V(219.75148278474705, -90.44615667066816, 0))], [])
	//	// let line = new L3(V(-1344.04574670165, 826.5930889273866, 720.915318266099), V(0.776732950940391, -0.43614824442447003, -0.45437939192802856))
	//	let line = new L3(V(-1560.8950828838565, 716.07295580975, 249.61382611323648), V(0.9130103135570956, -0.36545647611595106, -0.18125598308272678))
	//	let face2 = new PlaneFace(new PlaneSurface(new P3(V(0, 0, -1), 100)), [
	//		new PCurveEdge(new SemiEllipseCurve(V(73.03583314037537, -69.86032483338774, -100), V(-24.176861672352132, -146.16681457389276, 0), V(146.16681457389276, -24.176861672352132, 0)), V(219.75148278474705, -90.44615667066816, -100), V(97.2126948127275, 76.30648974050503, -100), 1.5953170840348225, -3.141592653589793, null, V(-20.58583183728038, -146.71564964437164, 0), V(146.16681457389276, -24.176861672352114, 0)),
	//		StraightEdge.throughPoints(V(97.2126948127275, 76.30648974050503, -100), V(275.99999999999966, 255.99999999999972, -100)),
	//		new PCurveEdge(new BezierCurve(V(219.75148278474705, -90.44615667066816, -100), V(-82.00000000000018, -138.00000000000023, -100), V(539.9999999999997, 225.9999999999997, -100), V(275.99999999999966, 255.99999999999972, -100), -0.1, 1.1), V(275.99999999999966, 255.99999999999972, -100), V(219.75148278474705, -90.44615667066816, -100), 1, 0, null, V(792, -90.00000000000009, 0), V(905.2544483542417, 142.6615299879962, 0))], [])
	//	console.log(face.intersectsLine(line), face.surface.isTsForLine(line), face.surface.isTsForLine(line))
	//	const t = line.isTWithPlane(face2.surface.plane);
	//	console.log(face2.intersectsLine(line), t, line.at(t).sce)
	//	assert.ok(face.intersectsLine(line))
	//},
	'serialization'(assert) {
		let a = {a: 2, b: 3}
		assert.equal(unserialize(serialize(a)).toString(), a.toString())

		a.c = a

		const a2 = unserialize(serialize(a))
		assert.equal(a2.a, 2)
		assert.equal(a2.b, 3)
		assert.equal(a2.c, a2)


		a = [1, 2, 3]
		assert.equal(unserialize(serialize(a)).toString(), a.toString())
	},
	'PlaneSurface.loopContainsPoint'(assert) {
		const loop = StraightEdge.chain([V(0, 0), V(10, 0), V(10, 10), V(0, 10)], true)
		assert.equal(new PlaneSurface(P3.XY).loopContainsPoint(loop, V(8, 10)), PointVsFace.ON_EDGE)
	},
	'angleRelativeNormal'(assert) {
		assert.fuzzyEquals(V3.X.angleRelativeNormal(V3.Y, V3.Z), Math.PI / 2)
		assert.fuzzyEquals(V3.X.angleRelativeNormal(V3.Y.negated(), V3.Z), -Math.PI / 2)
		// assert.fuzzyEquals(V3.Y.angleRelativeNormal(V(32, Math.sqrt(2), -Math.sqrt(2)), V3.X), -Math.PI / 4 )
		// assert.fuzzyEquals(V(-0.1, 1, 0).angleRelativeNormal(V(0.0, 0, -1), V3.X), -Math.PI / 2)
		// assert.fuzzyEquals(V(-0, 1, 0).angleRelativeNormal(V(2, 0, -1), V3.X), -Math.PI / 2)
	},
	'sumInPlaceTree'(assert) {
		assert.equal([].sumInPlaceTree(), 0)
		assert.equal([1].sumInPlaceTree(), 1)
		assert.equal([1, 2].sumInPlaceTree(), 3)
		assert.equal([1, 2, 3].sumInPlaceTree(), 6)
		assert.equal([1, 2, 3, 4].sumInPlaceTree(), 10)
	},
	'PlaneSurface.loopContainsPoint 2'(assert) {
		const loop = [
			new StraightEdge(new L3(V(2, 10, 0), V3.Z), V(2, 10, 3), V(2, 10, 5), 3, 5),
			new StraightEdge(new L3(V(0, 10, 5), V3.X), V(2, 10, 5), V(0, 10, 5), 2, 0),
			new StraightEdge(new L3(V(0, 10, 0), V3.Z), V(0, 10, 5), V(0, 10, 0), 5, 0),
			new StraightEdge(new L3(V(0, 10, 0), V3.X), V(0, 10, 0), V(10, 10, 0), 0, 10),
			new StraightEdge(new L3(V(10, 10, 0), V3.Z), V(10, 10, 0), V(10, 10, 5), 0, 5),
			new StraightEdge(new L3(V(0, 10, 5), V3.X), V(10, 10, 5), V(6, 10, 5), 10, 6),
			new StraightEdge(new L3(V(6, 10, 0), V(0, 0, -1)), V(6, 10, 5), V(6, 10, 3), -5, -3),
			new StraightEdge(new L3(V(0, 10, 3), V(-1, 0, 0)), V(6, 10, 3), V(2, 10, 3), -6, -2)
		]
		const p = V(6, 10, 3)
		testLoopContainsPoint(assert, new PlaneSurface(new P3(V(0,-1,0),-10)), loop, p, PointVsFace.ON_EDGE)
	},
    //'EllipsoidSurface.splitOnPlaneLoop'(assert) {
    //    //const es = SemiEllipsoidSurface.UNIT
    //    const a = V3.sphere(30 * DEG, 70 * DEG), z = a.z, xy = a.lengthXY(), center = V(0, 0, z), f1 = V(a.x, a.y, 0), f2 = V(-a.y, a.x)
    //    const curve = new SemiEllipseCurve(center, f1, f2)
    //    const seamCurve = SemiEllipseCurve.UNIT.rotateX(-PI / 2)
    //    const edge = Edge.forCurveAndTs(curve, -PI, PI)
    //    assert.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
		//				href='viewer.html?mesh=${es.sce}.toMesh()&points=[V(-5, 1, -1)]&edges=[${edge.str}]'>view</a>`)
    //    const [front, back] = SemiEllipsoidSurface.splitOnPlaneLoop([edge], true)
    //
    //    assert.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
		//				href='viewer.html?mesh=${es.sce}.toMesh()&points=[V(-5, 1, -1)]&edges=${back.sce}'>view</a>`)
    //    console.log(front, back)
    //    const expectedFront = []
    //    const expectedBack = [Edge.forCurveAndTs(curve, -120 * DEG, 60 * DEG), Edge.forCurveAndTs(seamCurve)]
    //},
    'V3.inverseLerp'(assert) {
	    const a = V(0.1, 0.2, 0.3), b = V(3, 2, 1)
        assert.equal(0, V3.inverseLerp(a, b, a))
        assert.equal(1, V3.inverseLerp(a, b, b))
        assert.equal(0.5, V3.inverseLerp(a, b, a.plus(b).div(2)))
    },
})