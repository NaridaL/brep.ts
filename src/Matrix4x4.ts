
import eq = NLA.eq
class M4 extends Matrix implements Transformable {
	m: Float64Array

	/**
	 * This constructor takes 16 arguments in row-major order, which can be passed individually, as a list, or even as
	 * four lists, one for each row. If the arguments are omitted then the identity matrix is constructed instead.
	 */
	constructor(...var_args: (number|number[])[]) {
		let m
		if (0 == arguments.length) {
			m = new Float64Array(16)
		} else {
			const flattened = Array.prototype.concat.apply([], arguments)
			assert(flattened.length == 16, "flattened.length == 16" + flattened.length)
			m = new Float64Array(flattened)
		}
		super(4, 4, m)
		let o = Object.create(M4.prototype)
		Object.defineProperty(o, 'm', {value: m})
		return o
	}

	/**
	 * Returns a new M4 which is equal to the inverse of this.
	 */
	inversed(): M4 {
		return M4.inverse(this)
	}

	/**
	 * Matrix trace is defined as the sum of the elements of the main diagonal.
	 */
	trace(): number {
		return this.m[0] + this.m[5] + this.m[10] + this.m[15]
	}

	as3x3(): M4 {
		let result = M4.copy(this), m = result.m

		m[3] = m[7] = m[11] = m[12] = m[13] = m[14] = 0
		m[15] = 1
		return result
	}

	transform(m4: M4): M4 {
		return m4.times(this)
	}

	realEigenValues3(): number[] {
		let m = this.m
		assert(0 == m[12] && 0 == m[13] && 0 == m[14])
		// determinant of (this - λI):
		// | a-λ  b   c  |
		// |  d  e-λ  f  | = -λ^3 + λ^2 (a+e+i) + λ (-a e-a i+b d+c g-e i+f h) + a(ei - fh) - b(di - fg) + c(dh - eg)
		// |  g   h  i-λ |

		let [a, b, c, ,
			d, e, f, ,
			g, h, i] = m
		// det(this - λI) = -λ^3 +λ^2 (a+e+i) + λ (-a e-a i-b d+c g-e i+f h)+ (a e i-a f h-b d i+b f g+c d h-c e g)
		let s = -1
		let t = a + e + i // equivalent to trace of matrix
		let u = -a * e - a * i + b * d + c * g - e * i + f * h // equivalent to 1/2 (trace(this²) - trace²(A))
		let w = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g) // equivalent to matrix determinant

		console.log(s, t, u, w)
		return solveCubicReal2(s, t, u, w)

	}

	realEigenVectors3():V3[] {
		let eigenValues = this.realEigenValues3()
		let this3x3 = this.times(M4.IDENTITY3)
		console.log(this.toString())
		console.log(this3x3.toString())
		let mats = eigenValues.map(ev => M4.IDENTITY3.timesScalar(-ev).plus(this3x3))
		console.log(mats.map(m=>m.determinant3()))
		console.log(mats.map(m=>'' + m.toString(v=>'' + v)).join("\n\n"))
		console.log(mats.map(m=>'' + m.gauss().U.toString(v=>'' + v)).join("\n\n"))
		console.log("mats.map(m=>m.rank())", mats.map(m=>m.rank()))
		if (1 == eigenValues.length) {
			console.log(mats[0].toString())
			assertf(() => 0 == mats[0].rank())
			// col vectors
			return NLA.arrayFromFunction(3, col => new V3(this.m[col], this.m[4 + col], this.m[8 + col]))
		}
		if (2 == eigenValues.length) {
			// one matrix should have rank 1, the other rank 2
			if (1 == mats[0].rank()) {
				mats = [mats[1], mats[0]]
			}
			assertf(() => 2 == mats[0].rank())
			assertf(() => 1 == mats[1].rank())

			// mat[0] has rank 2, mat[1] has rank 1
			let gauss0 = mats[0].gauss().U
			let eigenVector0 = gauss0.row(0).cross(gauss0.row(1)).V3().normalized()

			let planeNormal = mats[1].gauss().U.row(0).V3()
			let eigenVector1 = planeNormal.getPerpendicular().normalized()
			let eigenVector2 = eigenVector0.cross(eigenVector1).rejectedFrom(planeNormal)

			return [eigenVector0, eigenVector1, eigenVector2]
		}
		if (3 == eigenValues.length) {
			mats.forEach((mat, i) => assert(2 == mat.rank(), i + ': ' + mat.rank()))
			// the (A - lambda I) matrices map to a plane. This means, that there is an entire line in R³ which maps to the point V3.ZERO
			return mats.map(mat => {
				let gauss = mat.gauss().U
				return gauss.row(0).cross(gauss.row(1)).V3().normalized()
			})
		}
	}

	/**
	 * U * SIGMA * VSTAR = this
	 * U and VSTAR are orthogonal matrices
	 * SIGMA is a diagonal matrix
	 */
	svd3(): {U: M4, SIGMA: M4, VSTAR: M4} {
		function matrixForCS(i, k, c, s) {
			const m = M4.identity()
			m.setEl(i, i, c)
			m.setEl(k, k, c)
			m.setEl(i, k, s)
			m.setEl(k, i, -s)
			return m
		}

		const A = this.as3x3()
		let S = A.transposed().times(A), V = M4.identity()
		console.log(S.str)
		for (let it = 0; it < 16; it++) {
			console.log("blahg\n", V.times(S).times(V.transposed()).str)
			assert(V.times(S).times(V.transposed()).likeM4(A.transposed().times(A)),
				V.times(S).times(V.transposed()).str,
				A.transposed().times(A).str)
			let maxOffDiagonal = 0, maxOffDiagonalIndex = 1, j = 10
			while (j--) {
				const val = Math.abs(S.m[j])
				if (j % 4 != Math.floor(j / 4) && val > maxOffDiagonal) {
					maxOffDiagonal = val
					maxOffDiagonalIndex = j
				}
			}

			const i = Math.floor(maxOffDiagonalIndex / 4), k = maxOffDiagonalIndex % 4
			const a_ii = S.m[5 * i], a_kk = S.m[5 * k], a_ik = S.m[maxOffDiagonalIndex]
			const phi = a_ii === a_kk ? PI / 4 : Math.atan(2 * a_ik / (a_ii - a_kk)) / 2
			console.log(maxOffDiagonalIndex, i, k, "phi", phi)
			const cos = Math.cos(phi), sin = Math.sin(phi)
			const givensRotation = matrixForCS(i, k, cos, -sin)
			assert(givensRotation.transposed().times(givensRotation).isIdentity())
			console.log(givensRotation.str)
			V = V.times(givensRotation)
			S = M4.multiplyMultiple(givensRotation.transposed(), S, givensRotation)
			console.log(S.str)
		}

		const sigma = S.map((el, elIndex) => elIndex % 5 == 0 ? Math.sqrt(el) : 0) as M4
		return {U: M4.multiplyMultiple(A, V, sigma.map((el, elIndex) => elIndex % 5 == 0 ? 1 / el : 0)), SIGMA: sigma, VSTAR: V.transposed()}
	}

	map(fn: (el: number, elIndex: number, array: Float64Array) => number): M4 {
		return M4.fromFunction((x, y, i) => fn(this.m[i], i, this.m))
	}

	likeM4(m4): boolean {
		assertInst(M4, m4)
		return this.m.every((el, index) => NLA.eq(el, m4.m[index]))
	}

	/**
	 * Returns a new M4 equal to the transpose of this.
	 */
	transposed(): M4 {
		return M4.transpose(this)
	}

	/**
	 * Returns a new M4 which equal to (this * matrix) (in that order)
	 */
	times(matrix: M4): M4 {
		return M4.multiply(this, matrix)
	}

	/**
	 * Transforms the vector as a point with a w coordinate of 1. This means translations will have an effect, for
	 * example.
	 */
	transformPoint(v: V3): V3 {
		assertVectors(v)
		const m = this.m
		const vx = v.x, vy = v.y, vz = v.z, vw = 1
		const x = vx * m[0] + vy * m[1] + vz * m[2] + vw * m[3]
		const y = vx * m[4] + vy * m[5] + vz * m[6] + vw * m[7]
		const z = vx * m[8] + vy * m[9] + vz * m[10] + vw * m[11]
		const w = vx * m[12] + vy * m[13] + vz * m[14] + vw * m[15]
		// scale such that fourth element becomes 1:
		return new V3(x / w, y / w, z / w)
	}

	/**
	 * Transforms the vector as a vector with a w coordinate of 0. This means translations will have no effect, for
	 * example. Will throw an exception if the calculated w component != 0. This occurs for example when attempting
	 * to transform a vector with a perspective matrix.
	 */
	transformVector(v: V3):V3 {
		assertVectors(v)
		const m = this.m
		const w = v.x * m[12] + v.y * m[13] + v.z * m[14]
		assert(w === 0, () => 'w != 0 needs to be true for this to make sense (w =' + w + this.str)
		return new V3(m[0] * v.x + m[1] * v.y + m[2] * v.z, m[4] * v.x + m[5] * v.y + m[6] * v.z, m[8] * v.x + m[9] * v.y + m[10] * v.z)
	}

	transformedPoints(vs: V3[]): V3[] {
		return vs.map(v => this.transformPoint(v))
	}

	transformedVectors(vs: V3[]): V3[] {
		return vs.map(v => this.transformVector(v))
	}

	plus(m: M4): M4 {
		let r = new M4(), i = 16
		while (i--) r.m[i] = this.m[i] + m.m[i]
		return r
	}

	minus(m: M4): M4 {
		let r = new M4(), i = 16
		while (i--) r.m[i] = this.m[i] - m.m[i]
		return r
	}

	timesScalar(factor: number): M4 {
		let r = new M4(), i = 16
		while (i--) r.m[i] = this.m[i] * factor
		return r
	}

	divScalar(scalar: number): M4 {
		let r = new M4(), i = 16
		while (i--) r.m[i] = this.m[i] / scalar
		return r
	}

	copy(): M4 {
		return M4.copy(this)
	}

	isRegular(): boolean {
		return !NLA.eq0(this.determinant())
	}


	isAxisAligned():boolean {
		const m = this.m
		return (1 >= +!NLA.eq0(m[0]) + +!NLA.eq0(m[1]) + +!NLA.eq0(m[2]))
			&& (1 >= +!NLA.eq0(m[4]) + +!NLA.eq0(m[5]) + +!NLA.eq0(m[6]))
			&& (1 >= +!NLA.eq0(m[8]) + +!NLA.eq0(m[9]) + +!NLA.eq0(m[10]))
	}

	/**
	 * A matrix M is orthogonal iff M * M^T = I
	 * I being the identity matrix.
	 *
	 * @returns {boolean} If this matrix is orthogonal or very close to it. Comparison of the identity matrix and
	 * this * this^T is done with {@link #likeM4}
	 */
	isOrthogonal(): boolean {
		// return this.transposed().times(this).likeM4(M4.IDENTITY)
		M4.transpose(this, M4.temp0)
		M4.multiply(this, M4.temp0, M4.temp1)
		return M4.IDENTITY.likeM4(M4.temp1)
	}

	/**
	 * A matrix M is symmetric iff M == M^T
	 * I being the identity matrix.
	 *
	 * @returns {boolean} If this matrix is symmetric or very close to it. Comparison of the identity matrix and
	 * this * this^T is done with {@link #likeM4}
	 */
	isSymmetric(): boolean {
		M4.transpose(this, M4.temp0)
		return this.likeM4(M4.temp0)
	}

	/**
	 * A matrix M is normal iff M * M^-T == M^T * M TODO: ^-T?
	 * I being the identity matrix.
	 *
	 * @returns {boolean} If this matrix is symmetric or very close to it. Comparison of the identity matrix and
	 * this * this^T is done with {@link #likeM4}
	 */
	isNormal(): boolean {
		M4.transpose(this, M4.temp0) // temp0 = this^-T
		M4.multiply(this, M4.temp0, M4.temp1) // temp1 = this * this^-T
		M4.multiply(M4.temp0, this, M4.temp2) // temp2 = this^-T * this
		return M4.temp1.likeM4(M4.temp2)
	}

	/**
	 * Determinant of matrix.
	 *
	 * Notes:
	 *      For matrices A and B
	 *      det(A * B) = det(A) * det(B)
	 *      det(A^-1) = 1 / det(A)
	 */
	determinant():number {
		/*
		 | a b c d |
		 | e f g h |
		 | i j k l |
		 | m n o p |
		 */
		const $ = this.m,
			a = $[0], b = $[1], c = $[2], d = $[3],
			e = $[4], f = $[5], g = $[6], h = $[7],
			i = $[8], j = $[9], k = $[10], l = $[11],
			m = $[12], n = $[13], o = $[14], p = $[15],
			klop = k * p - l * o, jlnp = j * p - l * n, jkno = j * o - k * n,
			ilmp = i * p - l * m, ikmo = i * o - k * m, ijmn = i * n - j * m
		return (
			a * (f * klop - g * jlnp + h * jkno)
			- b * (e * klop - g * ilmp + h * ikmo)
			+ c * (e * jlnp - f * ilmp + h * ijmn)
			- d * (e * jkno - f * ikmo + g * ijmn))
	}

	determinant3(): number {
		const [a, b, c, ,
			   d, e, f, ,
			   g, h, i] = this.m
		const det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
		return det
	}

	/**
	 * determine whether this matrix is a mirroring transformation
	 */
	isMirroring():boolean {
		/*
		 var u = V(this.m[0], this.m[4], this.m[8])
		 var v = V(this.m[1], this.m[5], this.m[9])
		 var w = V(this.m[2], this.m[6], this.m[10])

		 // for a true orthogonal, non-mirrored base, u.cross(v) == w
		 // If they have an opposite direction then we are mirroring
		 var mirrorvalue = u.cross(v).dot(w)
		 var ismirror = (mirrorvalue < 0)
		 return ismirror
		 */

		return this.determinant() < 0 // TODO: also valid for 4x4?

	}

	/**
	 * Get the translation part of this matrix, i.e. the result of this.transformVector(V3.ZERO)
	 */
	getTranslation():V3 {
		const m = this.m, w = m[15]
		return new V3(m[3] / w, m[7] / w, m[11] / w)
	}

	/**
	 * Returns this matrix scaled so that the determinant is 1.
	 * det(c * A) = (c ** n) * det(A) for n x n matrices,
	 * so we need to divide by the 4th root of the determinant
	 * @returns {M4}
	 */
	normalized():M4 {
		const det = this.determinant()
		return 1 == det ? this : this.divScalar(Math.sqrt(Math.sqrt(det)))
	}

	/**
	 * Returns this matrix scaled so that the determinant is 1.
	 * det(c * A) = (c ** n) * det(A) for n x n matrices,
	 * so we need to divide by the 4th root of the determinant
	 * @returns {M4}
	 */
	normalized2() {
		let div = this.m[15]
		return 1 == div ? this : this.divScalar(Math.sqrt(Math.sqrt(div)))
	}

	/**
	 * Returns if the matrix has the following form (within NLA_PRECISION):
	 * a b c 0
	 * c d e 0
	 * f g h 0
	 * 0 0 0 1
	 */
	is3x3(): boolean {
		const m = this.m
		return NLA.eq(1, m[15])
			&& NLA.eq0(m[12]) && NLA.eq0(m[13]) && NLA.eq0(m[14])
			&& NLA.eq0(m[3]) && NLA.eq0(m[7]) && NLA.eq0(m[11])
	}

	isIdentity(): boolean {
		return this.m.every((val, i) => (i / 4 | 0) == (i % 4) ? NLA.eq(1, val) : NLA.eq0(val))
	}

	toString(f?: (number) => string): string {
		f = f || ((v) => v.toFixed(6).replace(/(0|\.)(?=0*$)/g, " ").toString())
		assert(typeof f(0) == "string", "" + typeof f(0))
		// slice this.m to convert it to an Array (from TypeArray)
		const rounded = Array.prototype.slice.call(this.m).map(f)
		const colWidths = [0, 1, 2, 3].map((colIndex) => rounded.sliceStep(colIndex, 4).map((x) => x.length).max())
		return [0, 1, 2, 3].map(
			(rowIndex) => rounded
				.slice(rowIndex * 4, rowIndex * 4 + 4) // select matrix row
				.map((x, colIndex) => NLA.repeatString(colWidths[colIndex] - x.length, ' ') + x) // pad numbers with spaces to col width
				.join(" ")
		).join("\n") // join rows
	}

	/**
	 * Returns the matrix that when multiplied with `matrix` results in the
	 * identity matrix. You can optionally pass an existing matrix in `result`
	 * to avoid allocating a new matrix. This implementation is from the Mesa
	 * OpenGL function `__gluInvertMatrixd()` found in `project.c`.
	 */
	static inverse(matrix: M4, result?: M4): M4 {
		assertInst(M4, matrix)
		!result || assertInst(M4, result)
		assert(matrix != result, "matrix != result")
		result = result || new M4()
		const m = matrix.m, r = result.m

		// first compute transposed cofactor matrix:
		// cofactor of an element is the determinant of the 3x3 matrix gained by removing the column and row belonging to the element
		r[0] = m[5] * m[10] * m[15] - m[5] * m[14] * m[11] - m[6] * m[9] * m[15] + m[6] * m[13] * m[11] + m[7] * m[9] * m[14] - m[7] * m[13] * m[10]
		r[1] = -m[1] * m[10] * m[15] + m[1] * m[14] * m[11] + m[2] * m[9] * m[15] - m[2] * m[13] * m[11] - m[3] * m[9] * m[14] + m[3] * m[13] * m[10]
		r[2] = m[1] * m[6] * m[15] - m[1] * m[14] * m[7] - m[2] * m[5] * m[15] + m[2] * m[13] * m[7] + m[3] * m[5] * m[14] - m[3] * m[13] * m[6]
		r[3] = -m[1] * m[6] * m[11] + m[1] * m[10] * m[7] + m[2] * m[5] * m[11] - m[2] * m[9] * m[7] - m[3] * m[5] * m[10] + m[3] * m[9] * m[6]

		r[4] = -m[4] * m[10] * m[15] + m[4] * m[14] * m[11] + m[6] * m[8] * m[15] - m[6] * m[12] * m[11] - m[7] * m[8] * m[14] + m[7] * m[12] * m[10]
		r[5] = m[0] * m[10] * m[15] - m[0] * m[14] * m[11] - m[2] * m[8] * m[15] + m[2] * m[12] * m[11] + m[3] * m[8] * m[14] - m[3] * m[12] * m[10]
		r[6] = -m[0] * m[6] * m[15] + m[0] * m[14] * m[7] + m[2] * m[4] * m[15] - m[2] * m[12] * m[7] - m[3] * m[4] * m[14] + m[3] * m[12] * m[6]
		r[7] = m[0] * m[6] * m[11] - m[0] * m[10] * m[7] - m[2] * m[4] * m[11] + m[2] * m[8] * m[7] + m[3] * m[4] * m[10] - m[3] * m[8] * m[6]

		r[8] = m[4] * m[9] * m[15] - m[4] * m[13] * m[11] - m[5] * m[8] * m[15] + m[5] * m[12] * m[11] + m[7] * m[8] * m[13] - m[7] * m[12] * m[9]
		r[9] = -m[0] * m[9] * m[15] + m[0] * m[13] * m[11] + m[1] * m[8] * m[15] - m[1] * m[12] * m[11] - m[3] * m[8] * m[13] + m[3] * m[12] * m[9]
		r[10] = m[0] * m[5] * m[15] - m[0] * m[13] * m[7] - m[1] * m[4] * m[15] + m[1] * m[12] * m[7] + m[3] * m[4] * m[13] - m[3] * m[12] * m[5]
		r[11] = -m[0] * m[5] * m[11] + m[0] * m[9] * m[7] + m[1] * m[4] * m[11] - m[1] * m[8] * m[7] - m[3] * m[4] * m[9] + m[3] * m[8] * m[5]

		r[12] = -m[4] * m[9] * m[14] + m[4] * m[13] * m[10] + m[5] * m[8] * m[14] - m[5] * m[12] * m[10] - m[6] * m[8] * m[13] + m[6] * m[12] * m[9]
		r[13] = m[0] * m[9] * m[14] - m[0] * m[13] * m[10] - m[1] * m[8] * m[14] + m[1] * m[12] * m[10] + m[2] * m[8] * m[13] - m[2] * m[12] * m[9]
		r[14] = -m[0] * m[5] * m[14] + m[0] * m[13] * m[6] + m[1] * m[4] * m[14] - m[1] * m[12] * m[6] - m[2] * m[4] * m[13] + m[2] * m[12] * m[5]
		r[15] = m[0] * m[5] * m[10] - m[0] * m[9] * m[6] - m[1] * m[4] * m[10] + m[1] * m[8] * m[6] + m[2] * m[4] * m[9] - m[2] * m[8] * m[5]

		// calculate determinant using laplace expansion (cf https://en.wikipedia.org/wiki/Laplace_expansion),
		// as we already have the cofactors. We multiply a column by a row as the cofactor matrix is transposed.
		const det = m[0] * r[0] + m[1] * r[4] + m[2] * r[8] + m[3] * r[12]
		// assert(!NLA.isZero(det), "det may not be zero, i.e. the matrix is not invertible")
		let i = 16
		while (i--) {
			r[i] /= det
		}
		return result
	}

	/**
	 * Returns `matrix`, exchanging columns for rows. You can optionally pass an
	 * existing matrix in `result` to avoid allocating a new matrix.
	 */
	static transpose(matrix: M4, result?: M4): M4 {
		assertInst(M4, matrix)
		!result || assertInst(M4, result)
		assert(matrix != result, "matrix != result" + matrix + result)
		result = result || new M4()
		const m = matrix.m, r = result.m
		r[0] = m[0]
		r[1] = m[4]
		r[2] = m[8]
		r[3] = m[12]
		r[4] = m[1]
		r[5] = m[5]
		r[6] = m[9]
		r[7] = m[13]
		r[8] = m[2]
		r[9] = m[6]
		r[10] = m[10]
		r[11] = m[14]
		r[12] = m[3]
		r[13] = m[7]
		r[14] = m[11]
		r[15] = m[15]
		return result
	}

	/**
	 * Returns the concatenation of the transforms for `left` and `right`.
	 */
	static multiply(left: M4, right: M4, result?: M4):M4 {
		assertInst(M4, left, right)
		!result || assertInst(M4, result)
		assert(left != result, "left != result")
		assert(right != result, "right != result")
		result = result || new M4()
		const a = left.m, b = right.m, r = result.m

		r[0] = a[0] * b[0] + a[1] * b[4] + (a[2] * b[8] + a[3] * b[12])
		r[1] = a[0] * b[1] + a[1] * b[5] + (a[2] * b[9] + a[3] * b[13])
		r[2] = a[0] * b[2] + a[1] * b[6] + (a[2] * b[10] + a[3] * b[14])
		r[3] = a[0] * b[3] + a[1] * b[7] + (a[2] * b[11] + a[3] * b[15])

		r[4] = a[4] * b[0] + a[5] * b[4] + (a[6] * b[8] + a[7] * b[12])
		r[5] = a[4] * b[1] + a[5] * b[5] + (a[6] * b[9] + a[7] * b[13])
		r[6] = a[4] * b[2] + a[5] * b[6] + (a[6] * b[10] + a[7] * b[14])
		r[7] = a[4] * b[3] + a[5] * b[7] + (a[6] * b[11] + a[7] * b[15])

		r[8] = a[8] * b[0] + a[9] * b[4] + (a[10] * b[8] + a[11] * b[12])
		r[9] = a[8] * b[1] + a[9] * b[5] + (a[10] * b[9] + a[11] * b[13])
		r[10] = a[8] * b[2] + a[9] * b[6] + (a[10] * b[10] + a[11] * b[14])
		r[11] = a[8] * b[3] + a[9] * b[7] + (a[10] * b[11] + a[11] * b[15])

		r[12] = a[12] * b[0] + a[13] * b[4] + (a[14] * b[8] + a[15] * b[12])
		r[13] = a[12] * b[1] + a[13] * b[5] + (a[14] * b[9] + a[15] * b[13])
		r[14] = a[12] * b[2] + a[13] * b[6] + (a[14] * b[10] + a[15] * b[14])
		r[15] = a[12] * b[3] + a[13] * b[7] + (a[14] * b[11] + a[15] * b[15])

		return result
	}

	static copy(src: M4, result?: M4) {
		assertInst(M4, src)
		!result || assertInst(M4, result)
		assert(result != src, "result != src")
		result = result || new M4()
		const s = src.m, d = result.m
		let i = 16
		while (i--) {
			d[i] = s[i]
		}
		return result
	}

	/**
	 *
	 * @param e0
	 * @param e1
	 * @param e2
	 * @param origin
	 * @returns {M4}
	 */
	static forSys(e0: V3, e1: V3, e2?: V3, origin?: V3):M4 {
		assertVectors(e0, e1)
		!e2 || assertVectors(e2)
		!origin || assertVectors(origin)

		e2 = e2 || e0.cross(e1)
		origin = origin || V3.ZERO
		return new M4(
			e0.x, e1.x, e2.x, origin.x,
			e0.y, e1.y, e2.y, origin.y,
			e0.z, e1.z, e2.z, origin.z,
			0, 0, 0, 1)
	}

	static forRows(n0: V3, n1: V3, n2: V3, n3?: V3):M4 {
		assertVectors(n0, n1, n2)
		!n3 || assertVectors(n2)
		n3 = n3 || V3.ZERO
		return new M4(
			n0.x, n0.y, n0.z, 0,
			n1.x, n1.y, n1.z, 0,
			n2.x, n2.y, n2.z, 0,
			n3.x, n3.y, n3.z, 1)
	}

	/**
	 * Returns an identity matrix. You can optionally pass an existing matrix in `result` to avoid allocating a new
	 * matrix. This emulates the OpenGL function `glLoadIdentity()`
	 *
	 * Unless initializing a matrix to be modified, use M4.IDENTITY
	 */
	static identity(result?: M4):M4 {
		!result || assertInst(M4, result)
		result = result || new M4()
		const m = result.m
		m[0] = m[5] = m[10] = m[15] = 1
		m[1] = m[2] = m[3] = m[4] = m[6] = m[7] = m[8] = m[9] = m[11] = m[12] = m[13] = m[14] = 0

		return result
	}

	/**
	 * Creates a new M4 initialized by a user defined callback function
	 *
	 * @param f signature: (elRow, elCol, elIndex) =>
	 *     el, where elIndex is the row-major index, i.e. eLindex == elRow * 4 + elCol
	 * @param result
	 */
	static fromFunction(f: (elRow: number, elCol: number, elIndex: number) => number, result?: M4): M4 {
		assert(typeof f == "function", 'typeof f == "function"' + typeof f)
		!result || assertInst(M4, result)
		result = result || new M4()
		const m = result.m
		let i = 16
		while (i--) {
			m[i] = f(Math.floor(i / 4), i % 4, i)
		}
		return result
	}

// ### GL.Matrix.perspective(fov, aspect, near, far[, result])
//
// Returns a perspective transform matrix, which makes far away objects appear
// smaller than nearby objects. The `aspect` argument should be the width
// divided by the height of your viewport and `fov` is the top-to-bottom angle
// of the field of view in degrees. You can optionally pass an existing matrix
// in `result` to avoid allocating a new matrix. This emulates the OpenGL
// function `gluPerspective()`.
	/**
	 *
	 * @param fov in degrees
	 * @param aspect Aspect ration, width/height of viewport
	 * @param near
	 * @param far
	 * @param {M4} [result]
	 * @returns {M4}
	 */
	static perspective(fov: number, aspect: number, near: number, far: number, result): M4 {
		!result || assertInst(M4, result)
		assertNumbers(fov, aspect, near, far)
		let y = Math.tan(fov * Math.PI / 360) * near
		let x = y * aspect
		return M4.frustum(-x, x, -y, y, near, far, result)
	}

// ### GL.Matrix.frustum(left, right, bottom, top, near, far[, result])
//
// Sets up a viewing frustum, which is shaped like a truncated pyramid with the
// camera where the point of the pyramid would be. You can optionally pass an
// existing matrix in `result` to avoid allocating a new matrix. This emulates
// the OpenGL function `glFrustum()`.
	static frustum(left: number, right: number, bottom: number, top: number, near: number, far: number, result?: M4): M4 {
		assertNumbers(left, right, bottom, top, near, far)
		assert(0 < near, "0 < near")
		assert(near < far, "near < far")
		!result || assertInst(M4, result)
		result = result || new M4()
		const m = result.m

		m[0] = 2 * near / (right - left)
		m[1] = 0
		m[2] = (right + left) / (right - left)
		m[3] = 0

		m[4] = 0
		m[5] = 2 * near / (top - bottom)
		m[6] = (top + bottom) / (top - bottom)
		m[7] = 0

		m[8] = 0
		m[9] = 0
		m[10] = -(far + near) / (far - near)
		m[11] = -2 * far * near / (far - near)

		m[12] = 0
		m[13] = 0
		m[14] = -1
		m[15] = 0

		return result
	}

	/**
	 * Returns a new M4 representing the a projection through/towards a point onto a plane.
	 */
	static projectPlanePoint(p: V3, plane: P3, result?: M4): M4 {
		assertVectors(p)
		assertInst(P3, plane)
		!result || assertInst(M4, result)
		result = result || new M4()
		const m = result.m
		const n = plane.normal, w = plane.w
		const np = n.dot(p)

		m[0] = p.x * n.x + w - np
		m[1] = p.x * n.y
		m[2] = p.x * n.z
		m[3] = -w * p.x

		m[4] = p.y * n.x
		m[5] = p.y * n.y + w - np
		m[6] = p.y * n.z
		m[7] = -w * p.y

		m[8] = p.z * n.x
		m[9] = p.z * n.y
		m[10] = p.z * n.z + w - np
		m[11] = -w * p.z

		m[12] = n.x
		m[13] = n.y
		m[14] = n.z
		m[15] = -np

		return result
	}


	/**
	 * Orthographic/orthogonal projection. Transforms the cuboid with the dimensions X: [left right] Y: [bottom, top]
	 * Z: [near far] to the cuboid X: [-1 1] Y [-1 1] Z [-1, 1]
	 */
	static ortho(left: number, right: number, bottom: number, top: number, near: number, far: number, result?: M4): M4 {
		assertNumbers(left, right, bottom, top, near, far)
		!result || assertInst(M4, result)
		result = result || new M4()
		const m = result.m

		m[0] = 2 / (right - left)
		m[1] = 0
		m[2] = 0
		m[3] = -(right + left) / (right - left)

		m[4] = 0
		m[5] = 2 / (top - bottom)
		m[6] = 0
		m[7] = -(top + bottom) / (top - bottom)

		m[8] = 0
		m[9] = 0
		m[10] = -2 / (far - near)
		m[11] = -(far + near) / (far - near)

		m[12] = 0
		m[13] = 0
		m[14] = 0
		m[15] = 1

		return result
	}

	/**
	 * This emulates the OpenGL function `glScale()`. You can optionally pass an existing matrix in `result` to avoid
	 * allocating a new matrix.
	 */
	static scaling(x: number, y: number, z: number, result?: M4): M4
	static scaling(v: V3, result?: M4): M4
	static scaling(scale: number, result?: M4): M4
	static scaling(x, y?, z?, result?): M4 {
		if (1 == arguments.length || 2 == arguments.length) {
			result = y
            if ('number' === typeof x) {
                y = x
                z = x
            }
			if (x instanceof V3) {
				y = x.y
				z = x.z
				x = x.x
			}
		} else {
			assert(3 == arguments.length || 4 == arguments.length, "3 == arguments.length || 4 == arguments.length")
			assertNumbers(x, y, z)
		}
		!result || assertInst(M4, result)

		result = result || new M4()
		const m = result.m

		m[0] = x
		m[1] = 0
		m[2] = 0
		m[3] = 0

		m[4] = 0
		m[5] = y
		m[6] = 0
		m[7] = 0

		m[8] = 0
		m[9] = 0
		m[10] = z
		m[11] = 0

		m[12] = 0
		m[13] = 0
		m[14] = 0
		m[15] = 1

		return result
	}

	/**
	 * This emulates the OpenGL function `glTranslate()`. You can optionally pass
	 * an existing matrix in `result` to avoid allocating a new matrix.
	 */
	static translation(x: number, y: number, z: number, result?: M4): M4
	static translation(v: V3, result?: M4): M4
	static translation(scale: number, result?: M4): M4
	static translation(x, y?, z?, result?): M4 {
		if (1 == arguments.length || 2 == arguments.length) {
			assertVectors(x)
			result = y
			if (x instanceof V3) {
				y = x.y
				z = x.z
				x = x.x
			} else if ('number' === typeof x) {
				y = x
				z = x
			}
		} else {
			assert(3 == arguments.length || 4 == arguments.length, "3 == arguments.length || 4 == arguments.length")
			assertNumbers(x, y, z)
		}
		!result || assertInst(M4, result)

		result = result || new M4()
		const m = result.m

		m[0] = 1
		m[1] = 0
		m[2] = 0
		m[3] = x

		m[4] = 0
		m[5] = 1
		m[6] = 0
		m[7] = y

		m[8] = 0
		m[9] = 0
		m[10] = 1
		m[11] = z

		m[12] = 0
		m[13] = 0
		m[14] = 0
		m[15] = 1

		return result
	}

	/**
	 * Returns a matrix that rotates by `a` degrees around the vector (x, y, z). You can optionally pass an existing
	 * matrix in `result` to avoid allocating a new matrix. This emulates the OpenGL function `glRotate()`.
	 */
	static rotation(radians: number, x: number, y: number, z: number, result?: M4): M4
	static rotation(radians: number, v: V3, result?: M4): M4
	static rotation(radians, x, y?, z?, result?): M4 {
		if (2 == arguments.length || 3 == arguments.length) {
			assertVectors(x)
			assertNumbers(radians)
			result = y
			y = x.y
			z = x.z
			x = x.x
		} else {
			assert(4 == arguments.length || 5 == arguments.length, "4 == arguments.length || 5 == arguments.length")
			assertNumbers(radians, x, y, z)
		}
		!result || assertInst(M4, result)
		// TODO
		assert(!V(x, y, z).isZero(), "!V(x, y, z).isZero()")

		result = result || new M4()
		const m = result.m

		const d = Math.sqrt(x * x + y * y + z * z)
		x /= d
		y /= d
		z /= d
		const cos = Math.cos(radians), sin = Math.sin(radians), t = 1 - cos

		m[0] = x * x * t + cos
		m[1] = x * y * t - z * sin
		m[2] = x * z * t + y * sin
		m[3] = 0

		m[4] = y * x * t + z * sin
		m[5] = y * y * t + cos
		m[6] = y * z * t - x * sin
		m[7] = 0

		m[8] = z * x * t - y * sin
		m[9] = z * y * t + x * sin
		m[10] = z * z * t + cos
		m[11] = 0

		m[12] = 0
		m[13] = 0
		m[14] = 0
		m[15] = 1

		return result
	}

	/**
	 * Returns a matrix that puts the camera at the eye point `ex, ey, ez` looking
	 * toward the center point `cx, cy, cz` with an up direction of `ux, uy, uz`.
	 * You can optionally pass an existing matrix in `result` to avoid allocating
	 * a new matrix. This emulates the OpenGL function `gluLookAt()`.
	 */
	static lookAt(eye: V3, focus: V3, up: V3, result?: M4): M4 {
		assert(3 == arguments.length || 4 == arguments.length, "3 == arguments.length || 4 == arguments.length")
		assertVectors(eye, focus, up)
		!result || assertInst(M4, result)

		result = result || new M4()
		const m = result.m

		const f = eye.minus(focus).normalized()
		const s = up.cross(f).normalized()
		const t = f.cross(s).normalized()

		m[0] = s.x
		m[1] = s.y
		m[2] = s.z
		m[3] = -s.dot(eye)

		m[4] = t.x
		m[5] = t.y
		m[6] = t.z
		m[7] = -t.dot(eye)

		m[8] = f.x
		m[9] = f.y
		m[10] = f.z
		m[11] = -f.dot(eye)

		m[12] = 0
		m[13] = 0
		m[14] = 0
		m[15] = 1

		return result
	}

	/**
	 * Create a rotation matrix for rotating around the X axis
	 */
	static rotationX(radians: number): M4 {
		assertNumbers(radians)
		const sin = Math.sin(radians), cos = Math.cos(radians)
		const els = [
			1, 0, 0, 0, 0, cos, -sin, 0, 0, sin, cos, 0, 0, 0, 0, 1
		]
		return new M4(els)
	}

	/**
	 * Create a rotation matrix for rotating around the Y axis
	 */
	static rotationY(radians: number): M4 {
		const sin = Math.sin(radians), cos = Math.cos(radians)
		const els = [
			cos, 0, sin, 0, 0, 1, 0, 0, -sin, 0, cos, 0, 0, 0, 0, 1
		]
		return new M4(els)
	}

	/**
	 * Create a rotation matrix for rotating around the Z axis
	 */
	static rotationZ(radians: number) {
		const sin = Math.sin(radians), cos = Math.cos(radians)
		const els = [
			cos, -sin, 0, 0, sin, cos, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1
		]
		return new M4(els)
	}

	/**
	 * New rotation matrix such that result.transformVector(a).isParallelTo(b) through smallest rotation.
	 * Performs no scaling.
	 */
	static rotationAB(a: V3, b: V3, result?: M4): M4 {
		// see http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
		assertVectors(a, b)
		!result || assertInst(M4, result)
		const rotationAxis = a.cross(b), rotationAxisLength = rotationAxis.length()
		if (eq0(rotationAxisLength)) {
			return M4.identity(result)
		}
		const radians = Math.atan2(rotationAxisLength, a.dot(b))
		return M4.rotationLine(V3.ZERO, rotationAxis, radians, result)
	}

	/**
	 * Matrix for rotation about arbitrary line defined by an anchor point and direction.
	 * rotationAxis does not need to be normalized
	 */
	static rotationLine(rotationAnchor: V3, rotationAxis: V3, radians: number, result?: M4): M4 {
		// see http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
		assertVectors(rotationAnchor, rotationAxis)
		assertNumbers(radians)
		!result || assertInst(M4, result)
		result = result || new M4()
		rotationAxis = rotationAxis.normalized()

		const ax = rotationAnchor.x, ay = rotationAnchor.y, az = rotationAnchor.z,
			  dx = rotationAxis.x, dy = rotationAxis.y, dz = rotationAxis.z
		const m = result.m, cos = Math.cos(radians), sin = Math.sin(radians)

		m[0] = dx * dx + (dy * dy + dz * dz) * cos
		m[1] = dx * dy * (1 - cos) - dz * sin
		m[2] = dx * dz * (1 - cos) + dy * sin
		m[3] = (ax * (dy * dy + dz * dz) - dx * (ay * dy + az * dz)) * (1 - cos) + (ay * dz - az * dy) * sin

		m[4] = dx * dy * (1 - cos) + dz * sin
		m[5] = dy * dy + (dx * dx + dz * dz) * cos
		m[6] = dy * dz * (1 - cos) - dx * sin
		m[7] = (ay * (dx * dx + dz * dz) - dy * (ax * dx + az * dz)) * (1 - cos) + (az * dx - ax * dz) * sin

		m[8] = dx * dz * (1 - cos) - dy * sin
		m[9] = dy * dz * (1 - cos) + dx * sin
		m[10] = dz * dz + (dx * dx + dy * dy) * cos
		m[11] = (az * (dx * dx + dy * dy) - dz * (ax * dx + ay * dy)) * (1 - cos) + (ax * dy - ay * dx) * sin

		m[12] = 0
		m[13] = 0
		m[14] = 0
		m[15] = 1

		return result
	}

	/**
	 * Create an affine matrix for mirroring into an arbitrary plane:
	 */
	static mirroring(plane: P3, result?: M4): M4 {
		assertInst(P3, plane)
		!result || assertInst(M4, result)
		const nx = plane.normal.x
		const ny = plane.normal.y
		const nz = plane.normal.z
		const w = plane.w
		result = result || new M4()
		const m = result.m

		m[0] = 1.0 - 2.0 * nx * nx
		m[1] = -2.0 * ny * nx
		m[2] = -2.0 * nz * nx
		m[3] = 0

		m[4] = -2.0 * nx * ny
		m[5] = 1.0 - 2.0 * ny * ny
		m[6] = -2.0 * nz * ny
		m[7] = 0

		m[8] = -2.0 * nx * nz
		m[9] = -2.0 * ny * nz
		m[10] = 1.0 - 2.0 * nz * nz
		m[11] = 0

		m[12] = 2.0 * nx * w
		m[13] = 2.0 * ny * w
		m[14] = 2.0 * nz * w
		m[15] = 1
		return result
	}

	/**
	 *
	 * @param plane
	 * @param dir Projection direction. Optional, if not specified plane normal will be used.
	 * @param result
	 */
	static projection(plane: P3, dir?: V3, result?: M4): M4 {
		// TODO: doc
		/**
		 * plane.normal DOT (p + lambda * dir) = w (1)
		 * extract lambda:
		 * plane.normal DOT p + lambda * plane.normal DOT dir = w
		 * lambda = (w - plane.normal DOT p) / plane.normal DOT dir
		 * result = p + lambda * dir
		 * result = p + dir * (w - plane.normal DOT p) / plane.normal DOT dir
		 * result =  w * dir / (plane.normal DOT dir) + p - plane.normal DOT p * dir / (plane.normal DOT dir) *
		 *

		 a + d * (w - n . a) / (nd)
		 a + dw - d * na
		 */
		assertInst(P3, plane)
		!dir || assertVectors(dir)
		!result || assertInst(M4, result)
		dir = dir || plane.normal
		const w = plane.w
		result = result || new M4()
		const m = result.m
		const nd = plane.normal.dot(dir)
		const {x: nx, y: ny, z: nz} = plane.normal
		const {x: dx, y: dy, z: dz} = dir.div(nd)
		/*
		 rejectedFrom: return this.minus(b.times(this.dot(b) / b.dot(b)))
		 return M4.forSys(
		 V3.X.rejectedFrom(plane.normal),
		 V3.Y.rejectedFrom(plane.normal),
		 V3.Z.rejectedFrom(plane.normal),
		 plane.anchor,
		 result
		 )
		 */

		m[0] = 1.0 - nx * dx
		m[1] = -ny * dx
		m[2] = -nz * dx
		m[3] = dx * w

		m[4] = -nx * dy
		m[5] = 1.0 - ny * dy
		m[6] = -nz * dy
		m[7] = dy * w

		m[8] = -nx * dz
		m[9] = -ny * dz
		m[10] = 1.0 - nz * dz
		m[11] = dz * w

		m[12] = 0
		m[13] = 0
		m[14] = 0
		m[15] = 1

		return result
	}

	static lineProjection(line: L3, result?: M4): M4 {
		assertInst(L3, line)
		!result || assertInst(M4, result)
		const ax = line.anchor.x, ay = line.anchor.y, az = line.anchor.z
		const dx = line.dir1.x, dy = line.dir1.y, dz = line.dir1.z
		result = result || new M4()
		const m = result.m

		/*
		 projectedOn: return b.times(this.dot(b) / b.dot(b))
		 */

		m[0] = dx * dx
		m[1] = dx * dy
		m[2] = dx * dz
		m[3] = ax

		m[4] = dy * dx
		m[5] = dy * dy
		m[6] = dy * dz
		m[7] = ay

		m[8] = dz * dx
		m[9] = dz * dy
		m[10] = dz * dz
		m[11] = az

		m[12] = 0
		m[13] = 0
		m[14] = 0
		m[15] = 1
		return result
	}

	static multiplyMultiple(...m4s: M4[]): M4 {
		if (0 == m4s.length) return M4.identity()
		let temp = M4.identity(), result = m4s[0].copy()
		for (let i = 1; i < m4s.length; i++) {
			M4.multiply(result, m4s[i], temp)

			;[temp, result] = [result, temp]
		}
		return result
	}

	static pointInversion(p: V3, result?: M4): M4 {
		assertVectors(p)
		!result || assertInst(M4, result)
		result = result || new M4()
		const m = result.m

		m[0] = -1
		m[1] = 0
		m[2] = 0
		m[3] = 2 * p.x

		m[4] = 0
		m[5] = -1
		m[6] = 0
		m[7] = 2 * p.y

		m[8] = 0
		m[9] = 0
		m[10] = -1
		m[11] = 2 * p.z

		m[12] = 0
		m[13] = 0
		m[14] = 0
		m[15] = 1
		return result
	}

	private static readonly temp0 = new M4()
	private static readonly temp1 = new M4()
	private static readonly temp2 = new M4()

	/**
	 * A simple (consists of integers), regular, non-orthogonal matrix, useful mainly for testing.
	 * M4.BAR = M4.FOO.inverse()
	 */
	static readonly FOO = new M4(
		0, 1, 1, 2,
		0.3, 0.4, 0.8, 13,
		2.1, 3.4, 5.5, 8.9,
		0, 0, 0, 1)

	static readonly BAR = M4.FOO.inversed()

	static readonly IDENTITY = M4.identity()
	static readonly YZX = M4.forSys(V3.Y, V3.Z, V3.X)
	static readonly ZXY = M4.forSys(V3.Z, V3.X, V3.Y)

	static IDENTITY3 = new M4(
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 0
	)

	xyAreaFactor(): number {
		return this.transformVector(V3.X).cross(this.transformVector(V3.Y)).length()
	}

	getX(): V3 {
		return this.transformVector(V3.X)
	}
	getY(): V3 {
		return this.transformVector(V3.Y)
	}
	getZ(): V3 {
		return this.transformVector(V3.Z)
	}
}
M4.prototype.height = 4
M4.prototype.width = 4
NLA.registerClass(M4)