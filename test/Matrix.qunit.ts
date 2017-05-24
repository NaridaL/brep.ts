{
	QUnit.module('Matrix')
	registerTests({
		'Matrix.isOrthogonal'(assert) {
			const a = Matrix.identity(4)
			assert.ok(a.isOrthogonal())
			let b = Matrix.fromRowArrays([Math.sqrt(2) / 2, Math.sqrt(2) / 2], [-Math.sqrt(2) / 2, Math.sqrt(2) / 2])
			assert.ok(a.isOrthogonal())
		},
		'Matrix.prototype.solveLinearSystem'(assert) {
			const a = Matrix.fromRowArrays([0, 1, 1], [1, 1, 1], [1, 2, 3])
			const b = Vector.from(1, 2, 3)
			const x = a.solveLinearSystem(b)
			assert.push(x.equals(Vector.from(1, 1, 0)), x, Vector.from(1, 1, 0))
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

			// all rows on plane through origin with normal1 = V(1, -1, 0)
			const c = Matrix.fromRowArrays([1, -1, 1], [1, 1, 1], [-1, 0, -1])
			assert.ok(a.rowsIndependent())
			assert.notOk(b.rowsIndependent())
			assert.notOk(c.rowsIndependent())
		},
		'Matrix.prototype.rank'(assert) {
			const a = Matrix.fromRowArrays([6, 3], [4, 3])
			const b = Matrix.fromRowArrays([6, 3], [12, 6])

			// all rows on plane through origin with normal1 = V(1, -1, 0)
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
		'M4.prototype.isAxisAligned'(assert) {
			assert.ok(M4.rotateX(Math.PI / 2).isAxisAligned())
			console.log(M4.rotateX(Math.PI / 4).toString())
			console.log(false + true + true)
			assert.notOk(M4.rotateX(Math.PI / 4).isAxisAligned())
		},
		'M4.prototype.rotationLine'(assert) {
			assert.matrixEquals(M4.rotateLine(V3.O, V3.X, 1), M4.rotateX(1))
			assert.matrixEquals(M4.rotateLine(V3.O, V3.Y, 1), M4.rotateY(1))
			assert.matrixEquals(M4.rotateLine(V3.O, V3.Z, 1), M4.rotateZ(1))

			const a = V(1, 2, 3), PI = Math.PI
			assert.matrixEquals(
				M4.rotateLine(a, V(1, 1, 0).unit(), 1),
				M4.multiplyMultiple(M4.translate(a), M4.rotateZ(PI / 4), M4.rotateX(1), M4.rotateZ(-PI / 4), M4.translate(a.negated())),
				'',
				1e-6)
		},
		'M4.multiplyMultiple'(assert) {
			assert.matrixEquals(M4.multiplyMultiple(M4.rotateX(1), M4.rotateZ(1)), M4.rotateX(1).times(M4.rotateZ(1)))
		},
		'M4.projection'(assert) {
			const plane = new P3(V(1, 2, 3).unit(), 5)
			const proj = M4.project(plane)
			console.log(proj.transformPoint(V(2, 4, 6)))
			assert.V3like(proj.transformPoint(V(2, 4, 6)), plane.anchor)
			assert.V3like(proj.transformVector(V(2, 4, 6)), V3.O)
			const p2 = V(3, 5, 22)
			assert.ok(proj.transformPoint(p2).minus(p2).isParallelTo(plane.normal1))
			assert.ok(plane.containsPoint(proj.transformPoint(p2)))
			assert.ok(proj.transformVector(p2).minus(p2).isParallelTo(plane.normal1))
			assert.ok(proj.transformVector(p2).isPerpendicularTo(plane.normal1))
		},
		'M4.projection 2'(assert) {
			[V(1, 1, 1), V(0, 0, -1)].forEach(dir => {
				const plane = new P3(V(1, 2, 3).unit(), 5)
				const proj = M4.project(plane, dir)
				const p2 = V(3, 5, 22)
				assert.ok(proj.transformPoint(p2).minus(p2).isParallelTo(dir))
				assert.ok(plane.containsPoint(proj.transformPoint(p2)))
				assert.ok(proj.transformVector(p2).minus(p2).isParallelTo(dir))
				assert.ok(proj.transformVector(p2).isPerpendicularTo(plane.normal1))
				console.log(proj.transformPoint(p2).sce)
				console.log(proj.str)
			})
		},
		'M4.isIdentity'(assert) {
			assert.ok(M4.identity().isIdentity())
			assert.notOk(M4.scale(1, 2, 3).isIdentity())
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
		'M4.projectPlanePoint()'(assert) {
			const m4 = M4.projectPlanePoint(V3.Z.negated(), P3.XY)
			assert.V3like(m4.transformPoint(V(4, 0, 1)), V(2, 0, 0))
			assert.V3like(m4.transformPoint(V(4, 8, 1)), V(2, 4, 0))
			assert.V3like(m4.transformPoint(V(4, 8, 2)), V(4 / 3, 8 / 3, 0))
			assert.matrixEquivalent(
				M4.projectPlanePoint(M4.FOO.transformPoint(V3.Z.negated()), P3.XY.transform(M4.FOO)),
				M4.multiplyMultiple(M4.FOO, m4, M4.BAR))

		},


	})
}
