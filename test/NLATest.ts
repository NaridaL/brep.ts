import round10 = NLA.round10
import eq2 = NLA.eq2
window.onerror = function (errorMsg, url, lineNumber, column, errorObj) {
	console.log(errorMsg, url, lineNumber, column, errorObj)
}
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
QUnit.testDifferentSystems = function (name, test, ...what) {
	if (!what.length) {
		what = [M4.IDENTITY, M4.FOO]
	}
	QUnit.test(name, assert => {
		what.forEach(m => {
			console.log(`TESTING '${name}' WITH '${m.name || m.toString()}`)
			assert.push(true, undefined, undefined, `TESTING '${name}' WITH '${m.name || m.toString()}'`)
			test(assert, m)
		})
	})
}

QUnit.testDifferentSystems('ProjectedCurveSurface', function (assert, m4) {
	let pp = V(0.5, 1)
	let curve = BezierCurve.graphXY(2, -3, -3, 2)
	let pcs = new ProjectedCurveSurface(curve, V3.Z).transform(m4)
	let p = pcs.parametricFunction()(pp.x, pp.y)
	console.log(p.sce, pcs.pointToParameterFunction()(p))
	assert.V3like(pcs.pointToParameterFunction()(p), pp, 'ptpf(pcs.pf(pp)) == pp')
})
QUnit.testDifferentSystems('ProjectedCurveSurface Face line intersection test', function (assert, m4) {
	let curve = BezierCurve.graphXY(2, -3, -3, 2)
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
})
QUnit.testDifferentSystems('Matrix4x4 eigenValues and eigenVectors', function (assert, /** M4*/ m4) {
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
		/*m4.isNormal() && */assert.ok(eigenVector.isPerpendicularTo(eigenVectors[(i + 1) % eigenVectors.length]), `eigenVector${i}.isPerpendicularTo(eigenVector${(i + 1) % eigenVectors.length})`)
		assert.ok(!eigenVector.isZero(), `'!eigenVector${i}.isZero()` + !eigenVector.isZero())
		assert.ok(eigenVector.isParallelTo(m4.transformVector(eigenVector)), `eigenVector${i}.isParallelTo(m4.transformVector(eigenVector))`)
	})

}, new M4(
	3, 2, 4, 0,
	2, 0, 2, 0,
	4, 2, 3, 0,
	0, 0, 0, 1), M4.FOO)

QUnit.testDifferentSystems('SemiEllipseCurve.isPointsWithBezier()', function (assert, /** M4*/ m4) {
	const ell = new SemiEllipseCurve(V(-223.34900663163222, -176.63214006755936, 0), V(-169.5891804980124, -35.54247345835796, 0), V(35.54247345835796, -169.5891804980124, 0))
    const bez = new BezierCurve(V(-267.6481190901419, -368.37017217006473, 0), V(563.959389388763, 94.96018577817034, 0), V(-1110.7787051488917, -95.8394860073627, 0), V(-59.14331799274822, -299.7830665459221, 0))
    const isPoints = ell.isPointsWithBezier(bez)
	assert.equal(isPoints.length, 3)
	isPoints.forEach(p => {
		assert.ok(ell.containsPoint(p), `ell.distanceToPoint(${p}) = ${ell.distanceToPoint(p)}`)
		assert.ok(bez.containsPoint(p), `bez.distanceToPoint(${p}) = ${bez.distanceToPoint(p)}`)
	})
})

QUnit.testDifferentSystems('M4.svd3()', function (assert, /** M4*/ m4) {

	const {U, SIGMA, VSTAR} = m4.svd3()
	const U_SIGMA_VSTAR = M4.multiplyMultiple(U, SIGMA, VSTAR)
	assert.ok(SIGMA.isDiagonal(), "SIGMA.isDiagonal()")
	assert.ok(U.isOrthogonal(), "U.isOrthogonal()\n" + U.str + "\n" + U.times(U.transposed()).str)
	assert.ok(VSTAR.isOrthogonal(), "VSTAR.isOrthogonal()\n" + VSTAR.str + "\n" + VSTAR.times(VSTAR.transposed()).str)
	assert.push(U_SIGMA_VSTAR.likeM4(m4.as3x3()), U_SIGMA_VSTAR.str, m4.as3x3().str, "U_SIGMA_VSTAR == A")
})

QUnit.testDifferentSystems('BezierCurve.isTsWithSurface(SemiCylinderSurface)', function (assert, m4) {
	let bez = BezierCurve.graphXY(2, -3, -3, 2, -2, 3).rotateX(15 * DEG).translate(0, 0, 100).transform(m4)
	let cyl = new SemiCylinderSurface(SemiEllipseCurve.forAB(4, 1).rotateY(10 * DEG), V3.Z).transform(m4)
	let ts = bez.isTsWithSurface(cyl)
	assert.equal(ts.length, 3)
	ts.forEach(t => {
		const p = bez.at(t)
		assert.ok(cyl.containsPoint(p), cyl.implicitFunction()(p))
	})
})
QUnit.testDifferentSystems('SemiCylinderSurface.calculateArea', function (assert, m4) {
	const surface = SemiCylinderSurface.UNIT.transform(m4)
	// loop which is 1 high and goes around a quarter of the cylinder
	const loop = [
		StraightEdge.throughPoints(V(1, 0, 1), V(1, 0, 0)),
		Edge.forCurveAndTs(SemiEllipseCurve.UNIT, 0, PI / 2),
		StraightEdge.throughPoints(V(0, 1, 0), V(0, 1, 1)),
		Edge.forCurveAndTs(SemiEllipseCurve.UNIT.translate(0, 0, 1), PI / 2, 0)].map(edge => edge.transform(m4))
	const face = Face.create(surface, loop)
	linkB2(assert, `mesh=${face.sce}.scale(100, 100, 100).toMesh()`)
	const area = face.calcArea()
	if (m4.isOrthogonal()) {
		assert.push(NLA.eq(area, PI/2), area, PI / 2)
	} else {
		const expectedArea = face.toMesh().calcVolume().area
		assert.push(NLA.eq2(area, expectedArea, 0.1), area, expectedArea)
	}


	const loopReverse = Edge.reverseLoop(loop)
	const holeArea = surface.calculateArea(loopReverse)
	if (m4.isOrthogonal()) {
		assert.push(NLA.eq(holeArea, -PI/2), area, -PI / 2)
	} else {
		const expectedArea = face.toMesh().calcVolume().area
		assert.push(NLA.eq2(holeArea, -expectedArea, 0.1), holeArea, -expectedArea)
	}

	const flippedSurfaceArea = surface.flipped().calculateArea(loop)
	if (m4.isOrthogonal()) {
		assert.push(NLA.eq(holeArea, -PI/2), area, -PI / 2)
	} else {
		const expectedArea = face.toMesh().calcVolume().area
		assert.push(NLA.eq2(holeArea, -expectedArea, 0.1), holeArea, -expectedArea)
	}


	{
		const loop = [
			StraightEdge.throughPoints(V(1, 0, 1), V(1, 0, 0)),
			Edge.forCurveAndTs(SemiEllipseCurve.forAB(1, 1), 0, PI/2),
			Edge.forCurveAndTs(new SemiEllipseCurve(V3.O, V(1, 0, 1), V(0, 1, 0)), PI / 2, 0)].map(edge => edge.transform(m4))
		const face = Face.create(surface, loop)
		linkB2(assert, `mesh=${face.sce}.scale(100, 100, 100).toMesh()`)
		const area = face.calcArea()
		if (m4.isOrthogonal()) {
			assert.push(NLA.eq(area, 1), area, 1)
		} else {
			const expectedArea = face.toMesh().calcVolume().area
			assert.push(NLA.eq2(area, expectedArea, 0.1), area, expectedArea)
		}


		const loopReverse = Edge.reverseLoop(loop)
		const holeArea = surface.calculateArea(loopReverse)
		if (m4.isOrthogonal()) {
			assert.push(NLA.eq(holeArea, -1), area, -1)
		} else {
			const expectedArea = face.toMesh().calcVolume().area
			console.log("expectedArea", expectedArea)
			console.log("expectedArea", eval(face.sce).toMesh().calcVolume().area)
			assert.push(NLA.eq2(holeArea, -expectedArea, 0.1), holeArea, -expectedArea)
		}
	}
}, M4.IDENTITY, M4.translation(0, 0, 4).times(M4.rotationY(10 * DEG)), M4.scaling(1, 2, 1))

QUnit.testDifferentSystems('EllipsoidSurface.calculateArea', function (assert, m4) {
	const surface = SemiEllipsoidSurface.UNIT.transform(m4)
	// loop which is 1 high and goes around a quarter of the cylinder
	const loop = [
		Edge.forCurveAndTs(SemiEllipseCurve.UNIT, 10 * DEG, 40 * DEG),
		Edge.forCurveAndTs(new SemiEllipseCurve(V3.O, V3.sphere(40 * DEG, 0), V3.Z), 0, PI / 2),
		Edge.forCurveAndTs(new SemiEllipseCurve(V3.O, V3.sphere(10 * DEG, 0), V3.Z), PI / 2, 0)].map(edge => edge.transform(m4))
	const face = Face.create(surface, loop)
	assert.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
						href='brep2.html?mesh=${face.sce}.scale(100, 100, 100).toMesh()'>view</a>`)
	const area = face.calcArea()
	const expectedArea = face.toMesh().calcVolume().area
	assert.push(NLA.eq2(area, expectedArea, 0.1), area, expectedArea)


	//const loopReverse = Edge.reverseLoop(loop)
	//const holeArea = surface.calculateArea(loopReverse)
	//if (m4.isOrthogonal()) {
	//	assert.push(NLA.eq(holeArea, -PI/2), area, -PI / 2)
	//} else {
	//	const expectedArea = face.toMesh().calcVolume().area
	//	assert.push(NLA.eq2(holeArea, -expectedArea, 0.1), holeArea, -expectedArea)
	//}
	//
	//const flippedSurfaceArea = surface.flipped().calculateArea(loop)
	//if (m4.isOrthogonal()) {
	//	assert.push(NLA.eq(holeArea, -PI/2), area, -PI / 2)
	//} else {
	//	const expectedArea = face.toMesh().calcVolume().area
	//	assert.push(NLA.eq2(holeArea, -expectedArea, 0.1), holeArea, -expectedArea)
	//}
	//
	//
	//{
	//	const loop = [
	//		StraightEdge.throughPoints(V(1, 0, 1), V(1, 0, 0)),
	//		Edge.forCurveAndTs(SemiEllipseCurve.forAB(1, -1), 0, -PI / 2),
	//		Edge.forCurveAndTs(new SemiEllipseCurve(V3.O, V(0, 1, 0), V(1, 0, 1)), 0, PI / 2)].map(edge => edge.transform(m4))
	//	const face = Face.create(surface, loop)
	//	assert.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
	//					href='brep2.html?mesh=${face.sce}.scale(100, 100, 100).toMesh()'>view</a>`)
	//	const area = face.calcArea()
	//	if (m4.isOrthogonal()) {
	//		assert.push(NLA.eq(area, 1), area, 1)
	//	} else {
	//		const expectedArea = face.toMesh().calcVolume().area
	//		assert.push(NLA.eq2(area, expectedArea, 0.1), area, expectedArea)
	//	}
	//
	//
	//	const loopReverse = Edge.reverseLoop(loop)
	//	const holeArea = surface.calculateArea(loopReverse)
	//	if (m4.isOrthogonal()) {
	//		assert.push(NLA.eq(holeArea, -1), area, -1)
	//	} else {
	//		const expectedArea = face.toMesh().calcVolume().area
	//		console.log("expectedArea", expectedArea)
	//		console.log("expectedArea", eval(face.sce).toMesh().calcVolume().area)
	//		assert.push(NLA.eq2(holeArea, -expectedArea, 0.1), holeArea, -expectedArea)
	//	}
	//}
}, M4.IDENTITY, M4.translation(0, 0, 4).times(M4.rotationY(10 * DEG)), M4.scaling(1, 1, 2))

QUnit.testDifferentSystems('SemiEllipseCurve.getAreaInDir', function (assert, m4) {
		let k = 1;
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
				let r = m4.transformVector(test.right)
				let areaFactor = m4.transformVector(V3.X).cross(m4.transformVector(V3.Y)).length()
				console.log(areaFactor)
				const ell = SemiEllipseCurve.UNIT.translate(0, yDiff, 0).transform(m4)
				let up = m4.transformVector(test.up).unit()
				let offsetArea = yDiff * ((1 - cos(test.t)) - (1 - cos(test.s))) * test.up.dot(V3.Y)
				const totalArea = test.result + offsetArea
				const expectedArea = totalArea * areaFactor
				let result = ell.getAreaInDir(r, up, test.s, test.t)
				let offsetCentroid = V((cos(test.t) + cos(test.s)) / 2, yDiff / 2)
				let movedCentroid = test.c.plus(V(0, yDiff))
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
	M4.rotationZ(45 * DEG),
	M4.forRows(V(1, 2, 3), V3.Y, V3.Z),
	M4.FOO.as3x3())

registerTests({

	'Edge.edgesIntersects'(assert) {
		let curve1 = BezierCurve.graphXY(2, -3, -3, 2, -2, 3)
		let curve2 = curve1.transform(M4.rotationLine(V(0.5, 0), V3.Z, PI / 2))
		let edge1 = PCurveEdge.forCurveAndTs(curve1, 0, 1)
		let edge2 = PCurveEdge.forCurveAndTs(curve2, 0, 1)
		assert.ok(Edge.edgesIntersect(edge1, edge2))
		assert.notOk(Edge.edgesIntersect(edge1, edge1.translate(10, 0, 0)))
		assert.notOk(Edge.edgesIntersect(edge1, edge2.translate(10, 0, 0)))
	},

	'M4.gauss'(assert) {
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
	},

	'NLA.eqAngle'(assert) {
		assert.ok(NLA.zeroAngle(0))
		assert.ok(NLA.zeroAngle(-NLA_PRECISION / 2))
		assert.ok(NLA.zeroAngle(2 * Math.PI - NLA_PRECISION / 2))
		assert.ok(NLA.zeroAngle(2 * Math.PI + NLA_PRECISION / 2))
		assert.ok(NLA.eqAngle(-Math.PI, Math.PI))
		assert.ok(NLA.eqAngle(0, 2 * Math.PI))
		assert.ok(NLA.eqAngle(0, 2 * Math.PI - NLA_PRECISION / 2))
		assert.ok(NLA.eqAngle(0, 2 * Math.PI + NLA_PRECISION / 2))
		assert.notOk(NLA.eqAngle(-Math.PI, 2 * Math.PI))
		assert.notOk(NLA.eqAngle(0, Math.PI))
	},

	'Vector.isParallelTo'(assert) {
		assert.equal(new Line(V3.O, V3.X).distanceToPoint(V(1, 1, 0)), 1)
	},
	'NLA.eq etc'(assert) {
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
	},
	'arrayCopyStep'(assert) {
		const a = [1, 2, 3, 4, 5, 6, 7, 8]
		const b = [-1, -2, -3, -4]
		NLA.arrayCopyStep(b, 0, 1, a, 1, 2, 3)
		assert.deepEqual(a, [1, -1, 3, -2, 5, -3, 7, 8])
	},
	'arrayCopy'(assert) {
		const a = [1, 2, 3, 4, 5, 6, 7, 8]
		const b = [-1, -2, -3, -4]
		NLA.arrayCopy(b, 1, a, 2, 2)
		assert.deepEqual(a, [1, 2, -2, -3, 5, 6, 7, 8])
	},
	'matrix rowArrays'(assert) {
		const a = Matrix.fromRowArrays([6, 3], [4, 3])
		assert.deepEqual(a.rowArray(0, Array), [6, 3])
		assert.deepEqual(a.rowArray(1, Array), [4, 3])
		assert.deepEqual(a.asRowArrays(Array), [[6, 3], [4, 3]])
	},
	'matrix transposed'(assert) {
		const a = Matrix.fromRowArrays([6, 3], [4, 3])
		assert.deepEqual(a.transposed().asRowArrays(Array), [[6, 4], [3, 3]])
	},
	'matrix transpose'(assert) {
		let a = Matrix.fromRowArrays([6, 3], [4, 3])
		a.transpose()
		assert.deepEqual(a.asRowArrays(Array), [[6, 4], [3, 3]])

		a = Matrix.fromRowArrays([6, 3, 4, 3])
		a.transpose()
		assert.deepEqual(a.asRowArrays(Array), [[6], [3], [4], [3]])
		a.transpose()
		assert.deepEqual(a.asRowArrays(Array), [[6, 3, 4, 3]])
	},
	'Matrix.prototype.times'(assert) {
		const a = Matrix.fromRowArrays([6, 3], [4, 3])
		assert.deepEqual(Matrix.identity(2).times(a).asRowArrays(Array), [[6, 3], [4, 3]])
		assert.deepEqual(a.times(Matrix.identity(2)).asRowArrays(Array), [[6, 3], [4, 3]])
	},
	'Matrix.identity'(assert) {
		const a = Matrix.identity(2)
		assert.deepEqual(a.asRowArrays(Array), [[1, 0], [0, 1]])
	},
	'Matrix.prototype.rowsIndependent'(assert) {
		const a = Matrix.fromRowArrays([6, 3], [4, 3])
		const b = Matrix.fromRowArrays([6, 3], [12, 6])

		// all rows on plane through origin with normal = V(1, -1, 0)
		const c = Matrix.fromRowArrays([1, -1, 1], [1, 1, 1], [-1, 0, -1])
		assert.ok(a.rowsIndependent())
		assert.notOk(b.rowsIndependent())
		assert.notOk(c.rowsIndependent())
	},
	'Matrix.prototype.rank'(assert) {
		const a = Matrix.fromRowArrays([6, 3], [4, 3])
		const b = Matrix.fromRowArrays([6, 3], [12, 6])

		// all rows on plane through origin with normal = V(1, -1, 0)
		const c = Matrix.fromRowArrays([1, -1, 1], [1, 1, 1], [-1, 0, -1])
		assert.equal(a.rank(), 2)
		assert.equal(b.rank(), 1)
		assert.equal(c.rank(), 2)

		const d = Matrix.fromRowArrays([1, 1, 0, 2], [-1, -1, 0, -2])
		assert.equal(d.rank(), 1)
		assert.equal(d.transposed().rank(), 1)


		let e = new M4(
			-60.16756109919886, 1, 1, 0,
			3, -56.16756109919886, 8, 0,
			21, 34, -5.167561099198863, 0,
			0, 0, 0, 1)
		console.log(e.rank())
		console.log(e.determinant())

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
	'Matrix.isOrthogonal'(assert) {
		const a = Matrix.identity(4)
		assert.ok(a.isOrthogonal())
		let b = Matrix.fromRowArrays([Math.sqrt(2) / 2, Math.sqrt(2) / 2], [-Math.sqrt(2) / 2, Math.sqrt(2) / 2])
		assert.ok(a.isOrthogonal())
	},
	'Matrix.prototype.solveLinearSystem'(assert) {
		const a = Matrix.fromRowArrays([0, 1, 1], [1, 1, 1], [1, 2, 3])
		const b = NLA.Vector.from(1, 2, 3)
		const x = a.solveLinearSystem(b)
		assert.push(x.equals(NLA.Vector.from(1, 1, 0)), x, NLA.Vector.from(1, 1, 0))
		assert.push(a.timesVector(x).equals(b), a.timesVector(x), b)
	},
	'Matrix.prototype.inverse'(assert) {
		const a = Matrix.fromRowArrays([0, 1, 1], [1, 1, 1], [1, 2, 3])
		const aInverse = a.inversed()
		console.log(aInverse.toString())
		assert.matrixEquals(a.times(aInverse), Matrix.identity(3))
	},
	'Matrix.prototype.inverse 2'(assert) {
		const a = Matrix.random(8, 8)
		const aInverse = a.inversed()
		assert.matrixEquals(a.times(aInverse), Matrix.identity(8))
	},
	'Matrix.prototype.inverse 3'(assert) {
		const a = new Matrix(1, 1, new Float64Array([5]))
		const aInverse = a.inversed()
		console.log(a.luDecomposition())
		assert.equal(aInverse.m[0], 1 / 5)
	},
	'Matrix.prototype.qrDecompositionGivensRotation'(assert) {
		const sqrt = Math.sqrt
		let m = Matrix.fromRowArrays([3, 5], [0, 2], [0, 0], [4, 5])
		const {Q, R} = m.qrDecompositionGivensRotation()
		assert.matrixEquals(Q, Matrix.fromRowArrays(
			[3 / 5, 4 / 5 / sqrt(5), 0, -8 / 5 / sqrt(5)],
			[0, 2 / sqrt(5), 0, 1 / sqrt(5)],
			[0, 0, 1, 0],
			[4 / 5, -3 / 5 / sqrt(5), 0, 6 / 5 / sqrt(5)]
		))
		assert.ok(Q.isOrthogonal())
		assert.matrixEquals(R, Matrix.fromRowArrays(
			[5, 7],
			[0, sqrt(5)],
			[0, 0],
			[0, 0]
		))
		assert.matrixEquals(Q.times(R), Matrix.fromRowArrays([3, 5], [0, 2], [0, 0], [4, 5]))
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
		let p = new P3(V(0, 0, 1), 2)
		assert.ok(P3.XY.like(P3.XY.transform(M4.identity())))
	},
	'Matrix4x4.prototype.isAxisAligned'(assert) {
		assert.ok(M4.rotationX(Math.PI / 2).isAxisAligned())
		console.log(M4.rotationX(Math.PI / 4).toString())
		console.log(false + true + true)
		assert.notOk(M4.rotationX(Math.PI / 4).isAxisAligned())
	},
	'Matrix4x4.prototype.rotationLine'(assert) {
		assert.matrixEquals(M4.rotationLine(V3.O, V3.X, 1), M4.rotationX(1))
		assert.matrixEquals(M4.rotationLine(V3.O, V3.Y, 1), M4.rotationY(1))
		assert.matrixEquals(M4.rotationLine(V3.O, V3.Z, 1), M4.rotationZ(1))

		const a = V(1, 2, 3), PI = Math.PI
		assert.matrixEquals(
			M4.rotationLine(a, V(1, 1, 0).unit(), 1),
			M4.multiplyMultiple(M4.translation(a), M4.rotationZ(PI / 4), M4.rotationX(1), M4.rotationZ(-PI / 4), M4.translation(a.negated())),
			'',
			1e-6)
	},
	'Matrix4x4.multiplyMultiple'(assert) {
		assert.matrixEquals(M4.multiplyMultiple(M4.rotationX(1), M4.rotationZ(1)), M4.rotationX(1).times(M4.rotationZ(1)))
	},
	'NLA.M4.projection'(assert) {
		const plane = new P3(V(1, 2, 3).unit(), 5)
		const proj = M4.projection(plane)
		console.log(proj.transformPoint(V(2, 4, 6)))
		assert.V3like(proj.transformPoint(V(2, 4, 6)), plane.anchor)
		assert.V3like(proj.transformVector(V(2, 4, 6)), V3.O)
		const p2 = V(3, 5, 22)
		assert.ok(proj.transformPoint(p2).minus(p2).isParallelTo(plane.normal))
		assert.ok(plane.containsPoint(proj.transformPoint(p2)))
		assert.ok(proj.transformVector(p2).minus(p2).isParallelTo(plane.normal))
		assert.ok(proj.transformVector(p2).isPerpendicularTo(plane.normal))
	},
	'NLA.M4.projection 2'(assert) {
		[V(1, 1, 1), V(0, 0, -1)].forEach(dir => {
			const plane = new P3(V(1, 2, 3).unit(), 5)
			const proj = M4.projection(plane, dir)
			const p2 = V(3, 5, 22)
			assert.ok(proj.transformPoint(p2).minus(p2).isParallelTo(dir))
			assert.ok(plane.containsPoint(proj.transformPoint(p2)))
			assert.ok(proj.transformVector(p2).minus(p2).isParallelTo(dir))
			assert.ok(proj.transformVector(p2).isPerpendicularTo(plane.normal))
			console.log(proj.transformPoint(p2).sce)
			console.log(proj.str)
		})
	},
	'NLA.M4.isIdentity'(assert) {
		assert.ok(M4.identity().isIdentity())
		assert.notOk(M4.scaling(1, 2, 3).isIdentity())
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
	'NLA.magic'(assert) {
		assert.expect(0)
		let a = V(1, 2, 3), b = V(4, 5, 6)
		//assert.ok(NLA.magic('a b c s => abs(a) x b .. c + 3 + ')(a, b, c, 3).equals(a.abs().cross(b).dot(c) + 3))
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
		const s = new NLA.CustomSet()
		const a = V(0, 2.7499999999999996, -5), b = V(0, 2.749999999999999, -5)
		s.canonicalizeLike(a)
		console.log(s._map, a.like(b), a.hashCodes(), b.hashCodes(), a.hashCode(), b.hashCode())
		assert.ok(s.canonicalizeLike(b) == a)
	},
	'NLA.arrayBinaryInsert'(assert) {
		const arr = [1, 2, 3, 4]
		NLA.arrayBinaryInsert(arr, 2.5, (a, b) => a - b)
		assert.deepEqual(arr, [1, 2, 2.5, 3, 4])

		const arr2 = []
		NLA.arrayBinaryInsert(arr2, -2, NLA.minus)
		NLA.arrayBinaryInsert(arr2, 5, NLA.minus)
		assert.deepEqual(arr2, [-2, 5])
	},
	'NLA.arrayBinaryIndexOf'(assert) {
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
        const results2 = solveCubic(a,b,c,d)
        assert.equal(results.length, 2)
        results.forEach(x => {
            const y = a * x ** 3 + b * x ** 2 + c * x + d
            assert.push(eq0(y), y, 0)
        })
    },
	'M4.projectPlanePoint()'(assert) {
		const m4 = M4.projectPlanePoint(V3.Z.negated(), P3.XY)
		assert.V3like(m4.transformPoint(V(4, 0, 1)), V(2, 0, 0))
		assert.V3like(m4.transformPoint(V(4, 8, 1)), V(2, 4, 0))
		assert.V3like(m4.transformPoint(V(4, 8, 2)), V(4 / 3, 8 / 3, 0))
		assert.matrixEquivalent(
			M4.projectPlanePoint(M4.FOO.transformPoint(V3.Z.negated()), P3.XY.transform(M4.FOO)),
			M4.multiplyMultiple(M4.FOO, m4, M4.BAR))

	},
	'ConicSurface.coplanar'(assert) {
		const unitCone = ConicSurface.UNIT
		assert.ok(unitCone.matrix.isIdentity(), 'unitCone.matrix.isIdentity()')
		assert.V3like(unitCone.parametricFunction()(0, 3), V(3, 0, 3))
		const ellipseAtZ3 = SemiEllipseCurve.UNIT.scale(3, 3, 3).translate(0, 0, 3)
		const planeAtZ3 = P3.XY.translate(0, 0, 3)
		const issAtZ3 = unitCone.isCurvesWithPlane(planeAtZ3)
		assert.equal(issAtZ3.length, 1)
		assert.push(ellipseAtZ3.isColinearTo(issAtZ3[0]), issAtZ3.toString(), ellipseAtZ3.toString())
		assert.ok(unitCone.containsEllipse(ellipseAtZ3))


		const scaledUnit = ConicSurface.UNIT.scale(2, 2, 1)
		assert.notOk(scaledUnit.isCoplanarTo(unitCone))
		assert.notOk(unitCone.isCoplanarTo(scaledUnit))
		const ell1 = unitCone.isCurvesWithPlane(new P3(V(2, 3, 10).unit(), 10))[0]
		assert.ok(unitCone.containsEllipse(ell1), 'unitCone.containsEllipse(ell1)')
		const ell2 = unitCone.isCurvesWithPlane(new P3(V(1, 1, 2).unit(), 4))[0]
		const ell1Cone = ConicSurface.atApexThroughEllipse(V3.O, ell1, 1)
		const ell2Cone = ConicSurface.atApexThroughEllipse(V3.O, ell2, 1)
		console.log(ell1Cone)
		assert.ok(unitCone.isCoplanarTo(ell1Cone))
		assert.ok(unitCone.isCoplanarTo(ell2Cone))
		assert.ok(ell1Cone.isCoplanarTo(ell2Cone))
		assert.ok(ell2Cone.isCoplanarTo(ell1Cone))
		assert.ok(ell1Cone.foo().isCoplanarTo(ell2Cone.foo()))
	},
	'ConicSurface.containsParabola'(assert) {
		const unitCone = ConicSurface.UNIT
		const pb = unitCone.isCurvesWithPlane(new P3(V(1, 0, 1).unit(), 4))[0]
		assert.ok(unitCone.containsParabola(pb))

		const c2 = unitCone.shearedX(2, 3)
		const pb2 = c2.isCurvesWithPlane(new P3(V(1, 0, 1).unit(), 4).shearedX(2, 3))[0]
		assert.ok(c2.containsParabola(pb2))
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
	//	const t = line.intersectWithPlaneLambda(face2.surface.plane);
	//	console.log(face2.intersectsLine(line), t, line.at(t).sce)
	//	assert.ok(face.intersectsLine(line))
	//},
	'SemiCylinderSurface.intersectionLine'(assert) {
		const cylSurface = SemiCylinderSurface.semicylinder(5)
		const line = L3.throughPoints(V(10, 0, 0), V(-10, 2, 10))
		const isPoints = cylSurface.isTsForLine(line).map(line.at, line)

		assert.equal(isPoints.length, 2, 'no of points')
		assert.notOk(isPoints[0].like(isPoints[1]))

		assert.ok(cylSurface.containsPoint(isPoints[0]))
		assert.ok(cylSurface.containsPoint(isPoints[1]))

		assert.ok(line.containsPoint(isPoints[0]), '' + line.distanceToPoint(isPoints[0]))
		assert.ok(line.containsPoint(isPoints[1]), '' + line.distanceToPoint(isPoints[1]))
	},
	'SemiCylinderSurface.zDirVolume'(assert) {
		const face = B2T.extrudeEdges([Edge.forCurveAndTs(SemiEllipseCurve.forAB(-1, 1), 0, PI), StraightEdge.throughPoints(V3.X, V3.X.negated())], P3.XY.flipped(), V3.Z, 'cyl').faces.find(face => face.surface instanceof SemiCylinderSurface)
		const face2 = B2T.extrudeEdges([Edge.forCurveAndTs(SemiEllipseCurve.UNIT, PI, 0), StraightEdge.throughPoints(V3.X, V3.X.negated())], P3.XY.flipped(), V3.Z, 'cyl').faces.find(face => face.surface instanceof SemiCylinderSurface)
		const face3 = B2T.extrudeEdges([Edge.forCurveAndTs(SemiEllipseCurve.UNIT, PI, 0).rotateY(-80 * DEG), StraightEdge.throughPoints(V3.X, V3.X.negated()).rotateY(-80 * DEG)], P3.XY.flipped().rotateY(-80 * DEG), new V3(-10, -1, 0).unit(), 'cyl').faces.find(face => face.surface instanceof SemiCylinderSurface)
		const modface = face.rotateY(-45 * DEG).translate(1, 0, 2)
		const e0 = modface.contour[0].project(new P3(modface.surface.dir1, 0))
		const face4 = Face.create(modface.surface, [e0, StraightEdge.throughPoints(e0.b, modface.contour[2].a), modface.contour[2], StraightEdge.throughPoints(modface.contour[2].b, e0.a)])

		testZDirVolume(assert, face)
		testZDirVolume(assert, face.rotateY(-45 * DEG).translate(1, 0, 2))
		testZDirVolume(assert, face.rotateY(90 * DEG).translate(1, 0, 2))

		testZDirVolume(assert, face2)
		testZDirVolume(assert, face2.rotateY(-45 * DEG).translate(1, 0, 2))
		testZDirVolume(assert, face2.rotateY(90 * DEG).translate(1, 0, 2))

		testZDirVolume(assert, face3)
		testZDirVolume(assert, face3.translate(1, 0, 2))

		testZDirVolume(assert, face4)
		testZDirVolume(assert, face4.translate(1, 0, 2))
	},
	'SemiCylinderSurface.loopContainsPoint'(assert) {
		const surface = new SemiCylinderSurface(new SemiEllipseCurve(V(0, 0, 0), V(8, 0, 0), V(0, 8, 0)), V(0, 0, -1))
		const loop = [
			StraightEdge.throughPoints(V(1, 7.937253933193773, 4), V(1, 7.937253933193773, 1)),
			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 1), V(8, 0, 0), V(0, 8, 0)), V(1, 7.937253933193773, 1), V(6, 5.291502622129181, 1), 1.4454684956268313, 0.7227342478134156, null, V(7.937253933193772, -0.9999999999999991, 0), V(5.2915026221291805, -6, 0)),
			StraightEdge.throughPoints(V(6, 5.291502622129181, 1), V(6, 5.291502622129181, 4)),
			new PCurveEdge(new SemiEllipseCurve(V(0, 0, 4), V(8, 0, 0), V(0, 8, 0)), V(6, 5.291502622129181, 4), V(1, 7.937253933193773, 4), 0.7227342478134156, 1.4454684956268313, null, V(-5.2915026221291805, 6, 0), V(-7.937253933193772, 0.9999999999999991, 0))
		]
		testLoopContainsPoint(assert, surface, loop, V(8, 0, 0), PointVsFace.OUTSIDE)
		testLoopContainsPoint(assert, surface, loop, V(1, 7.937253933193773, 3), PointVsFace.ON_EDGE)
	},
    'EllipsoidSurface.loopContainsPoint'(assert) {
        const testFace = B2T.rotateEdges([
            Edge.forCurveAndTs(SemiEllipseCurve.UNIT, 0, 90 * DEG).rotateX(90 * DEG),
            StraightEdge.throughPoints(V3.Z, V3.X)], 45 * DEG, 'blah')
            .faces.find(face => face.surface instanceof SemiEllipsoidSurface)

        const p1 = V3.sphere(10 * DEG, 10 * DEG)
	    testLoopContainsPoint(assert, testFace.surface, testFace.contour, p1, PointVsFace.INSIDE)
        const p2 = V3.sphere(10 * DEG, -10 * DEG)
	    testLoopContainsPoint(assert, testFace.surface, testFace.contour, p2, PointVsFace.OUTSIDE)
	    testLoopContainsPoint(assert,
		    testFace.surface.foo(),
		    testFace.contour.map(edge => edge.foo()),
		    M4.FOO.transformPoint(p1),
		    PointVsFace.INSIDE)
    },
    'EllipsoidSurface.loopContainsPoint 2'(assert) {
        const testFace = B2T.sphere(1).faces[0]

        const p1 = new V3(0, 1, 0)
	    testLoopContainsPoint(assert, testFace.surface, testFace.contour, p1, PointVsFace.INSIDE)
    },
	'serialization'(assert) {
		let a = {a: 2, b: 3}
		assert.equal(unserialize(serialize(a)).toString(), a.toString())

		a.c = a

		let a2 = unserialize(serialize(a))
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
	'EllipsoidSurface.mainAxes'(assert) {
		const es = new SemiEllipsoidSurface(V3.O, V(5, 0, -1), V(5, 1, 1), V(5, -1, 1))
		assert.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
						href='brep2.html?mesh=${es.sce}.toMesh()&points=[V(-5, 1, -1)]&edges=[
						StraightEdge.throughPoints(V(-5, 1, -1), ${V(-5, 1, -1).plus(V(-9.660064978873681e-7, 0.999999999962889, 0.000008560880502697739).times(100)).sce}),
						StraightEdge.throughPoints(V(-5, 1, -1), ${V(-5, 1, -1).plus(V(0.14427746420619014, -0.6925318281897137, -1.413919149220665).times(100)).sce})]'>view</a>`)
		const esn = es.mainAxes()
		assert.ok(esn.f1.isPerpendicularTo(esn.f2))
		assert.ok(esn.f2.isPerpendicularTo(esn.f3))
		assert.ok(esn.f3.isPerpendicularTo(esn.f1))
		assert.ok(es.inverseMatrix.times(esn.matrix).isOrthogonal())
	},
    //'EllipsoidSurface.splitOnPlaneLoop'(assert) {
    //    //const es = SemiEllipsoidSurface.UNIT
    //    const a = V3.sphere(30 * DEG, 70 * DEG), z = a.z, xy = a.lengthXY(), center = V(0, 0, z), f1 = V(a.x, a.y, 0), f2 = V(-a.y, a.x)
    //    const curve = new SemiEllipseCurve(center, f1, f2)
    //    const seamCurve = SemiEllipseCurve.UNIT.rotateX(-PI / 2)
    //    const edge = Edge.forCurveAndTs(curve, -PI, PI)
    //    assert.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
		//				href='brep2.html?mesh=${es.sce}.toMesh()&points=[V(-5, 1, -1)]&edges=[${edge.str}]'>view</a>`)
    //    const [front, back] = SemiEllipsoidSurface.splitOnPlaneLoop([edge], true)
    //
    //    assert.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
		//				href='brep2.html?mesh=${es.sce}.toMesh()&points=[V(-5, 1, -1)]&edges=${back.sce}'>view</a>`)
    //    console.log(front, back)
    //    const expectedFront = []
    //    const expectedBack = [Edge.forCurveAndTs(curve, -120 * DEG, 60 * DEG), Edge.forCurveAndTs(seamCurve)]
    //},
    'EllipsoidSurface.isCurvesWithPlane'(assert) {
        const es = SemiEllipsoidSurface.sphere(5)
        testISCurves(assert, es, new PlaneSurface(new P3(V(0, -1, 0.1).unit(), 4)), 0)
        testISCurves(assert, es, new PlaneSurface(new P3(V(0, -1, 0.1).unit(), 4)).flipped(), 0)
        testISCurves(assert, es, new PlaneSurface(new P3(V3.sphere(-PI/2, sin(3/5)), 4)), 0)
        testISCurves(assert, es, new PlaneSurface(new P3(V(0, -0.1, 1).unit(), 4)), 1)
        testISCurves(assert, es, new PlaneSurface(new P3(V(0, 0, 1), 4)), 1)
        testISCurves(assert, es, new PlaneSurface(new P3(V(0, 0.1, 1).unit(), 4)), 2)
        testISCurves(assert, es, new PlaneSurface(new P3(V(0, 1, 0.1).unit(), 4)), 2)
        testISCurves(assert, es, new PlaneSurface(P3.XY), 1)
        testISCurves(assert, es, new PlaneSurface(P3.XY.flipped()), 1)

        // slices perpendicular to V3.Y
        testISCurves(assert, es, new PlaneSurface(P3.ZX), 2)
        testISCurves(assert, es, new PlaneSurface(P3.ZX.flipped()), 2)
        testISCurves(assert, es, new PlaneSurface(P3.ZX).translate(0, 2, 0), 2)
        testISCurves(assert, es, new PlaneSurface(P3.ZX).translate(0, -2, 0), 0)

        // slices perpendicular to V3.X
        testISCurves(assert, es, new PlaneSurface(P3.YZ), 1)
    },
    'V3.inverseLerp'(assert) {
	    const a = V(0.1, 0.2, 0.3), b = V(3, 2, 1)
        assert.equal(0, V3.inverseLerp(a, b, a))
        assert.equal(1, V3.inverseLerp(a, b, b))
        assert.equal(0.5, V3.inverseLerp(a, b, a.plus(b).div(2)))
    },
})