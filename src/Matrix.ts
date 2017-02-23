import Vector = NLA.Vector
class Matrix implements Equalable {
	m: Float64Array
	width: number
	height: number

	constructor(width: int, height: int, m: Float64Array) {
		assert(width * height == m.length, "width * height == m.length", width, height, m.length)
		this.m = m
		this.width = width
		this.height = height
	}

	static random(width, height): Matrix {
		assertNumbers(width, height)
		return Matrix.fromFunction(width, height, (i, j) => Math.random())
	}

	static fromFunction(width, height, f): Matrix {
		assertNumbers(width, height)
		const m = new Float64Array(height * width)
		let elIndex = height * width
		while (elIndex--) {
			m[elIndex] = f(Math.floor(elIndex / width), elIndex % width, elIndex)
		}
		return new Matrix(width, height, m)
	}

	static identity(dim: number): Matrix {
		assertNumbers(dim)
		const m = new Float64Array(dim * dim)
		// Float64Arrays are init to 0
		let elIndex = dim * (dim + 1)
		while (elIndex) {
			elIndex -= (dim + 1)
			m[elIndex] = 1
		}
		return new Matrix(dim, dim, m)
	}

	static permutation(dim, i, k): Matrix {
		assertNumbers(dim, i, k)
		const m = new Float64Array(dim * dim)
		// Float64Array are init to 0
		let elIndex = dim * (dim + 1)
		while (elIndex) {
			elIndex -= (dim + 1)
			m[elIndex] = 1
		}
		m[i * dim + i] = 0
		m[k * dim + k] = 0
		m[i * dim + k] = 1
		m[k * dim + i] = 1
		return new Matrix(dim, dim, m)
	}

	static fromRowArrays(...args:number[][]): Matrix {
		return Matrix.fromRowArrays2(arguments)
	}

	static fromRowArrays2(arrays) {
		if (0 == arrays.length) {
			throw new Error("cannot have 0 vector")
		}
		const height = arrays.length
		const width = arrays[0].length
		const m = new Float64Array(height * width)
		NLA.arrayCopy(arrays[0], 0, m, 0, width)
		for (let rowIndex = 1; rowIndex < height; rowIndex++) {
			if (arrays[rowIndex].length != width) {
				throw new Error("all row arrays must be the same length")
			}
			NLA.arrayCopy(arrays[rowIndex], 0, m, rowIndex * width, width)
		}
		return new Matrix(width, height, m)
	}

	static fromColVectors(colVectors): Matrix {
		return Matrix.fromColArrays(colVectors.map((v) => v.v))
	}

	static forWidthHeight(width: int, height: int): Matrix {
		return new Matrix(width, height, new Float64Array(width * height))
	}

	static fromColArrays(colArrays): Matrix {
		if (0 == colArrays.length) {
			throw new Error("cannot have 0 vector")
		}
		const width = colArrays.length
		const height = colArrays[0].length
		const m = new Float64Array(height * width)
		NLA.arrayCopyStep(colArrays[0], 0, 1, m, 0, width, height)
		for (let colIndex = 1; colIndex < width; colIndex++) {
			if (colArrays[colIndex].length != height) {
				throw new Error("all col arrays must be the same length")
			}
			NLA.arrayCopyStep(colArrays[colIndex], 0, 1, m, colIndex, width, height)
		}
		return new Matrix(width, height, m)
	}


	e(rowIndex: number, colIndex: number): number {
		assertNumbers(rowIndex, colIndex)
		if (NLA_DEBUG && (rowIndex >= this.height || colIndex >= this.width)) {
			throw new Error("index " + rowIndex + ", " + colIndex + " is out of bounds (" + this.width + " x " + this.height + ")")
		}
		return this.m[rowIndex * this.width + colIndex]
	}

	setEl(rowIndex: number, colIndex: number, val: number): void {
		assertNumbers(rowIndex, colIndex, val)
		assert(0 <= rowIndex && rowIndex < this.height, "rowIndex out of bounds " + rowIndex)
		assert(0 <= colIndex && colIndex < this.width, "colIndex out of bounds " + colIndex)
		this.m[rowIndex * this.width + colIndex] = val
	}

	toString(f?: (el: number) => string): string {
		f = f || ((v) => v.toFixed(6))
		assert(typeof f(0) == "string", "" + typeof f(0))
		const rounded = Array.prototype.slice.call(this.m).map(f)
		const colWidths = NLA.arrayFromFunction(this.width,
			(colIndex) => rounded.sliceStep(colIndex, this.width).map((x) => x.length).max())
		return NLA.arrayFromFunction(this.height,
			(rowIndex) => rounded.slice(rowIndex * this.width, (rowIndex + 1) * this.width) // select matrix row
				.map((x, colIndex) => NLA.repeatString(colWidths[colIndex] - x.length, ' ') + x) // pad numbers with spaces to col width
				.join("  ")
		).map(x => x + "\n").join(""); // join rows
	}

	row(rowIndex): Vector {
		const v = new Float64Array(this.width)
		NLA.arrayCopy(this.m, rowIndex * this.width, v, 0, this.width)
		return new Vector(v)
	}

	col(colIndex): Vector {
		const v = new Float64Array(this.height)
		NLA.arrayCopyStep(this.m, colIndex, this.width, v, 0, 1, this.height)
		return new Vector(v)
	}

	dim(): {width: int, height: int} {
		return {width: this.width, height: this.height}
	}

	dimString(): string {
		return this.width + "x" + this.height
	}

	equals(obj): boolean {
		if (obj.constructor != Matrix) return false
		if (this.width != obj.width || this.height != obj.height) return false
		let elIndex = this.m.length
		while (elIndex--) {
			if (this.m[elIndex] != obj.m[elIndex]) return false
		}
	}

	equalsMatrix(matrix: Matrix, precision?: number): boolean {
		precision = precision || NLA_PRECISION
		if (!(matrix instanceof Matrix)) throw new Error("not a matrix")
		if (this.width != matrix.width || this.height != matrix.height) return false
		let elIndex = this.m.length
		while (elIndex--) {
			if (Math.abs(this.m[elIndex] - matrix.m[elIndex]) >= precision) return false
		}
		return true
	}

	hashCode(): int {
		let result = 0
		let elIndex = this.m.length
		while (elIndex--) {
			result = result * 31 + NLA.floatHashCode(this.m[elIndex])
		}
		return result
	}

	isZero(): boolean {
		let elIndex = this.m.length
		while (elIndex--) {
			if (!NLA.eq0(this.m[elIndex])) {
				return false
			}
		}
		return true
	}

	isOrthogonal(): boolean {
		return this.isSquare() && this.transposed().times(this).equalsMatrix(Matrix.identity(this.width))
	}

	/**
	 * Returns L, U, P such that L * U == P * this
	 */
	luDecomposition(): {L: Matrix, U: Matrix, P: Matrix} {
		assertf(() => this.isSquare(), this.dim().toSource())
		let dim = this.width
		let uRowArrays = this.asRowArrays(Float64Array)
		let lRowArrays = NLA.arrayFromFunction(dim, (row) => new Float64Array(dim))
		let pRowArrays = Matrix.identity(dim).asRowArrays(Float64Array)
		let currentRowIndex = 0
		for (let colIndex = 0; colIndex < dim; colIndex++) {
			// find largest value in colIndex
			let maxAbsValue = 0, pivotRowIndex = undefined, numberOfNonZeroRows: number = 0
			for (let rowIndex = currentRowIndex; rowIndex < dim; rowIndex++) {
				const el: number = uRowArrays[rowIndex][colIndex]
				numberOfNonZeroRows += +(0 != el)
				if (Math.abs(el) > maxAbsValue) {
					maxAbsValue = Math.abs(el)
					pivotRowIndex = rowIndex
				}
			}
			// TODO: check with NLA.isZero
			if (0 == maxAbsValue) {
				// column contains only zeros
				continue
			}
			assert(undefined !== pivotRowIndex)
			// swap rows
			NLA.arraySwap(uRowArrays, currentRowIndex, pivotRowIndex)
			NLA.arraySwap(lRowArrays, currentRowIndex, pivotRowIndex)
			NLA.arraySwap(pRowArrays, currentRowIndex, pivotRowIndex)
			lRowArrays[colIndex][colIndex] = 1

			if (1 < numberOfNonZeroRows) {
				// subtract pivot (now current) row from all below it
				for (let rowIndex = currentRowIndex + 1; rowIndex < dim; rowIndex++) {
					let l = uRowArrays[rowIndex][colIndex] / uRowArrays[currentRowIndex][colIndex]
					lRowArrays[rowIndex][colIndex] = l
					// subtract pivot row * l from row "rowIndex"
					for (let colIndex2 = colIndex; colIndex2 < dim; colIndex2++) {
						uRowArrays[rowIndex][colIndex2] -= l * uRowArrays[currentRowIndex][colIndex2]
					}
				}
			}
			currentRowIndex++ // this doesn't increase if pivot was zero
		}
		return {
			L: Matrix.fromRowArrays2(lRowArrays),
			U: Matrix.fromRowArrays2(uRowArrays),
			P: Matrix.fromRowArrays2(pRowArrays)
		}
	}

	gauss(): {L: Matrix, U: Matrix, P: Matrix} {
		let width = this.width, height = this.height
		let uRowArrays = this.asRowArrays(Float64Array)
		let lRowArrays = NLA.arrayFromFunction(height, (row) => new Float64Array(width))
		let pRowArrays = Matrix.identity(height).asRowArrays(Float64Array)
		let currentRowIndex = 0
		for (let colIndex = 0; colIndex < width; colIndex++) {
			// console.log('currentRowIndex', currentRowIndex)	// find largest value in colIndex
			let maxAbsValue = 0, pivotRowIndex = undefined, numberOfNonZeroRows = 0
			for (let rowIndex = currentRowIndex; rowIndex < height; rowIndex++) {
				const el = uRowArrays[rowIndex][colIndex]
				numberOfNonZeroRows += +(0 != el)
				if (Math.abs(el) > maxAbsValue) {
					maxAbsValue = Math.abs(el)
					pivotRowIndex = rowIndex
				}
			}
			// TODO: check with NLA.isZero
			if (0 == maxAbsValue) {
				// column contains only zeros
				continue
			}
			assert(undefined !== pivotRowIndex)
			// swap rows
			NLA.arraySwap(uRowArrays, currentRowIndex, pivotRowIndex)
			NLA.arraySwap(lRowArrays, currentRowIndex, pivotRowIndex)
			NLA.arraySwap(pRowArrays, currentRowIndex, pivotRowIndex)
			lRowArrays[currentRowIndex][colIndex] = 1

			if (1 < numberOfNonZeroRows) {
				// subtract pivot (now current) row from all below it
				for (let rowIndex = currentRowIndex + 1; rowIndex < height; rowIndex++) {
					let l = uRowArrays[rowIndex][colIndex] / uRowArrays[currentRowIndex][colIndex]
					lRowArrays[rowIndex][colIndex] = l
					// subtract pivot row * l from row "rowIndex"
					for (let colIndex2 = colIndex; colIndex2 < width; colIndex2++) {
						uRowArrays[rowIndex][colIndex2] -= l * uRowArrays[currentRowIndex][colIndex2]
					}
				}
			}
			currentRowIndex++ // this doesn't increase if pivot was zero
		}
		return {
			L: Matrix.fromRowArrays2(lRowArrays),
			U: Matrix.fromRowArrays2(uRowArrays),
			P: Matrix.fromRowArrays2(pRowArrays)
		}
	}

	qrDecompositionGivensRotation(): {Q: Matrix, R: Matrix} {
		function sigma(c, s) {
			if (0 == c) {
				return 1
			}
			if (Math.abs(s) < Math.abs(c)) {
				return 0.5 * Math.sign(c) * s
			}
			return 2 * Math.sign(s) / c
		}

		function matrixForCS(dim, i, k, c, s) {
			const m = Matrix.identity(dim)
			m.setEl(i, i, c)
			m.setEl(k, k, c)
			m.setEl(i, k, s)
			m.setEl(k, i, -s)
			return m
		}

		let qTransposed = Matrix.identity(this.height)
		for (let colIndex = 0; colIndex < this.width; colIndex++) {
			// find largest value in colIndex
			for (let rowIndex = colIndex + 1; rowIndex < this.height; rowIndex++) {
				//console.log("row ", rowIndex, "col ", colIndex)
				const xi = this.e(colIndex, colIndex)
				const xk = this.e(rowIndex, colIndex)
				if (xk == 0) {
					continue
				}
				const r = Math.sqrt(xi * xi + xk * xk)
				const c = xi / r
				const s = xk / r

				// apply transformation on every column:
				for (let col2 = colIndex; col2 < this.width; col2++) {
					const x1 = this.e(colIndex, col2) * c + this.e(rowIndex, col2) * s
					const x2 = this.e(rowIndex, col2) * c - this.e(colIndex, col2) * s
					this.setEl(colIndex, col2, x1)
					this.setEl(rowIndex, col2, x2)
				}
				//console.log("r ", r, "c ", c, "s ", s, "sigma", sigma(c, s))
				//console.log(this.toString(),"cs\n", matrixForCS(this.height, colIndex, rowIndex, c, s).toString())
				qTransposed = matrixForCS(this.height, colIndex, rowIndex, c, s).times(qTransposed)
			}
		}
		//console.log(qTransposed.transposed().toString(), this.toString(), qTransposed.transposed().times(this).toString())
		return {Q: qTransposed.transposed(), R: this}
	}

	isPermutation(): boolean {
		if (!this.isSquare()) return false
		if (this.m.some((value) => !NLA.eq0(value) && !NLA.eq(1, value))) return false

		const rows = this.asRowArrays(Array)
		if (rows.some((row) => row.filter((value) => NLA.eq(1, value)).length != 1)) return false

		const cols = this.asColArrays(Array)
		if (cols.some((col) => col.filter((value) => NLA.eq(1, value)).length != 1)) return false

		return true
	}

	isDiagonal(precision?: number): boolean {
		let i = this.m.length
		while (i--) {
			if (0 !== i % (this.width + 1) && !NLA.eq0(this.m[i])) return false
		}
		return true
	}

	isIdentity(precision?: number): boolean {
		return this.isLowerUnitriangular(precision) && this.isUpperTriangular(precision)
	}

	isUpperTriangular(precision?: number) {
		precision = "number" == typeof precision ? precision : NLA_PRECISION
		if (!this.isSquare()) return false
		for (let rowIndex = 1; rowIndex < this.height; rowIndex++) {
			for (let colIndex = 0; colIndex < rowIndex; colIndex++) {
				if (!NLA.eq02(this.m[rowIndex * this.width + colIndex], precision)) {
					return false
				}
			}
		}
		return true
	}

	/**
	 * Returns x, so that this * x = b
	 * More efficient than calculating the inverse for few (~ <= this.height) values
	 */
	solveLinearSystem(b: Vector): Vector {
		const lup = this.luDecomposition()
		// console.log(lup.L.toString())
		// console.log(lup.U.toString())
		// console.log(lup.P.toString())
		const y = lup.L.solveForwards(lup.P.timesVector(b))
		const x = lup.U.solveBackwards(y)
		return x
	}

	isLowerUnitriangular(precision?: number): boolean {
		precision = "number" == typeof precision ? precision : NLA_PRECISION
		if (!this.isSquare()) return false
		for (let rowIndex = 0; rowIndex < this.height - 1; rowIndex++) {
			for (let colIndex = rowIndex; colIndex < this.width; colIndex++) {
				const el = this.m[rowIndex * this.width + colIndex]
				if (rowIndex == colIndex ? !NLA.eq2(1, el, precision) : !NLA.eq02(el, precision)) {
					return false
				}
			}
		}
		return true
	}

	isLowerTriangular(): boolean {
		if (!this.isSquare()) return false
		for (let rowIndex = 0; rowIndex < this.height - 1; rowIndex++) {
			for (let colIndex = rowIndex + 1; colIndex < this.width; colIndex++) {
				if (!NLA.eq0(this.m[rowIndex * this.width + colIndex])) {
					return false
				}
			}
		}
		return true
	}

	solveBackwards(x: Vector): Vector {
		assertVectors(x)
		assert(this.height == x.dim(), "this.height == x.dim()")
		assert(this.isUpperTriangular(), "this.isUpperTriangular()\n" + this.str)
		const v = new Float64Array(this.width)
		let rowIndex = this.height
		while (rowIndex--) {
			let temp = x.v[rowIndex]
			for (let colIndex = rowIndex + 1; colIndex < this.width; colIndex++) {
				temp -= v[colIndex] * this.e(rowIndex, colIndex)
			}
			v[rowIndex] = temp / this.e(rowIndex, rowIndex)
		}
		return new Vector(v)
	}

	solveBackwardsMatrix(matrix: Matrix): Matrix {
		const colVectors = new Array(matrix.width)
		let i = matrix.width
		while (i--) {
			colVectors[i] = this.solveBackwards(matrix.col(i))
		}
		return Matrix.fromColVectors(colVectors)
	}

	solveForwardsMatrix(matrix: Matrix): Matrix {
		const colVectors = new Array(matrix.width)
		let i = matrix.width
		while (i--) {
			colVectors[i] = this.solveForwards(matrix.col(i))
		}
		return Matrix.fromColVectors(colVectors)
	}

	solveForwards(x: Vector): Vector {
		assertVectors(x)
		assert(this.height == x.dim(), "this.height == x.dim()")
		assertf(() => this.isLowerTriangular(), this.toString())
		const v = new Float64Array(this.width)
		for (let rowIndex = 0; rowIndex < this.height; rowIndex++) {
			let temp = x.v[rowIndex]
			for (let colIndex = 0; colIndex < rowIndex; colIndex++) {
				temp -= v[colIndex] * this.e(rowIndex, colIndex)
			}
			v[rowIndex] = temp / this.e(rowIndex, rowIndex)
		}
		return new Vector(v)
	}


	/**
	 * Calculates rank of matrix.
	 * Number of linearly independant row/column vectors.
	 * Is equal to the unmber of dimensions the image of the affine transformation represented this matrix has.
	 */
	rank(): int {
		let U = this.gauss().U
		//console.log(R.toString())
		let rowIndex = this.height
		while (rowIndex-- && U.row(rowIndex).isZero()) {
			console.log("RANK" + U.row(rowIndex).toString() + U.row(rowIndex).isZero())
		}
		return rowIndex + 1
	}

	rowsIndependent(): boolean {
		return this.height == this.rank()
	}

	colsIndependent(): boolean {
		return this.width == this.rank()
	}

	asRowArrays<T extends FloatArray>(arrayConstructor: new (length: int) => T): T[] {
		arrayConstructor = arrayConstructor || Float64Array
		let rowIndex = this.height
		const result = new Array(this.height)
		while (rowIndex--) {
			result[rowIndex] = this.rowArray(rowIndex, arrayConstructor)
		}
		return result
	}

	asColArrays<T extends FloatArray>(arrayConstructor: new (length: int) => T): T[] {
		arrayConstructor = arrayConstructor || Float64Array
		let colIndex = this.width
		const result = new Array(this.width)
		while (colIndex--) {
			result[colIndex] = this.colArray(colIndex, arrayConstructor)
		}
		return result
	}

	rowArray<T extends FloatArray>(rowIndex, arrayConstructor: new (length: int) => T): T {
		arrayConstructor = arrayConstructor || Float64Array
		const result = new arrayConstructor(this.width)
		NLA.arrayCopy(this.m, rowIndex * this.width, result, 0, this.width)
		return result
	}

	colArray<T extends FloatArray>(colIndex, arrayConstructor: new (length: int) => T): T {
		arrayConstructor = arrayConstructor || Float64Array
		const result = new arrayConstructor(this.width)
		NLA.arrayCopyStep(this.m, colIndex, this.height, result, 0, 1, this.height)
		return result
	}

	subMatrix(firstColIndex: int, subWidth: int, firstRowIndex: int, subHeight: int): Matrix {
		assert(firstColIndex + subWidth > this.width || firstRowIndex + subHeight > this.height)
		const m = new Float64Array(this.height)
		NLA.arrayCopyBlocks(this.m, firstColIndex, this.width, m, 0, subWidth, subHeight, subWidth)
		return new Matrix(subWidth, subHeight, m)
	}

	map(fn: (el: number, elIndex: number, array: Float64Array) => number): Matrix {
		return new Matrix(this.width, this.height, this.m.map(fn))
	}

	dimEquals(matrix: Matrix): boolean {
		assertInst(Matrix, matrix)
		return this.width == matrix.width && this.height == matrix.height
	}

	inversed(): Matrix {
		let lup = this.luDecomposition()
		let y = lup.L.solveForwardsMatrix(lup.P)
		console.log(y)
		let inverse = lup.U.solveBackwardsMatrix(y)
		return inverse
	}

	inversed3(): Matrix {
		assertf(() => 3 == this.width && 3 == this.height)
		let result = Matrix.forWidthHeight(3, 3), m = this.m, r = result.m

		r[0] = m[4] * m[8] - m[5] * m[7]
		r[1] = -m[1] * m[8] + m[2] * m[7]
		r[2] = m[1] * m[5] - m[2] * m[4]

		r[3] = -m[3] * m[8] + m[5] * m[6]
		r[4] = m[0] * m[8] - m[2] * m[6]
		r[5] = -m[0] * m[5] + m[2] * m[3]

		r[6] = m[3] * m[7] - m[4] * m[6]
		r[7] = -m[0] * m[7] + m[1] * m[6]
		r[8] = m[0] * m[4] - m[1] * m[3]

		let det = m[0] * r[0] + m[1] * r[3] + m[2] * r[6]
		let i = 9
		while (i--) {
			r[i] /= det
		}

		return result
	}

	inversed2(): Matrix {
		assertf(() => 2 == this.width && 2 == this.height)
		let result = Matrix.forWidthHeight(2, 2), m = this.m, r = result.m

		let det = m[0] * m[3] - m[1] * r[2]

		r[0] = m[3] / det
		r[1] = -m[2] / det

		r[2] = -m[1] / det
		r[3] = m[0] / det

		return result
	}

	canMultiply(matrix: Matrix): boolean {
		assertInst(Matrix, matrix)
		return this.width == matrix.height
	}

	times(matrix: Matrix): Matrix {
		assertInst(Matrix, matrix)
		assert(this.canMultiply(matrix), `Cannot multiply this {this.dimString()} by matrix {matrix.dimString()}`)
		const nWidth = matrix.width, nHeight = this.height, n = this.width
		const nM = new Float64Array(nWidth * nHeight)
		let nRowIndex = nHeight
		while (nRowIndex--) {
			let nColIndex = nWidth
			while (nColIndex--) {
				let result = 0
				let i = n
				while (i--) {
					result += this.m[nRowIndex * n + i] * matrix.m[i * nWidth + nColIndex]
				}
				nM[nRowIndex * nWidth + nColIndex] = result
			}
		}
		return new Matrix(nWidth, nHeight, nM)
	}

	timesVector(v: Vector): Vector {
		assertVectors(v)
		assert(this.width == v.dim())
		const nHeight = this.height, n = this.width
		const nM = new Float64Array(nHeight)
		let nRowIndex = nHeight
		while (nRowIndex--) {
			let result = 0
			let i = n
			while (i--) {
				result += this.m[nRowIndex * n + i] * v.v[i]
			}
			nM[nRowIndex] = result
		}
		return new Vector(nM)
	}

	transposed(): Matrix {
		const tWidth = this.height, tHeight = this.width
		const tM = new Float64Array(tWidth * tHeight)
		let tRowIndex = tHeight
		while (tRowIndex--) {
			let tColIndex = tWidth
			while (tColIndex--) {
				tM[tRowIndex * tWidth + tColIndex] = this.m[tColIndex * tHeight + tRowIndex]
			}
		}
		return new Matrix(tWidth, tHeight, tM)
	}

	/**
	 * In-place transpose.
	 */
	transpose(): void {
		const h = this.height, w = this.width, tM = this.m
		let tRowIndex = h
		while (tRowIndex--) {
			let tColIndex = Math.min(tRowIndex, w)
			while (tColIndex--) {
				console.log("col", tColIndex, "row", tRowIndex)
				const temp = tM[tRowIndex * w + tColIndex]
				tM[tRowIndex * w + tColIndex] = tM[tColIndex * h + tRowIndex]
				tM[tColIndex * h + tRowIndex] = temp
			}
		}
		this.width = h
		this.height = w
	}

	isSquare(): boolean {
		return this.height == this.width
	}

	diagonal(): Vector {
		if (!this.isSquare()) {
			throw new Error("!!")
		}
		const v = new Float64Array(this.width)
		let elIndex = this.width * (this.width + 1)
		let vIndex = this.width
		while (vIndex--) {
			elIndex -= this.width + 1
			v[vIndex] = this.m[elIndex]
		}
		return new Vector(v)
	}

	maxEl(): number {
		return Math.max.apply(undefined, this.m)
	}

	minEl(): number {
		return Math.min.apply(undefined, this.m)
	}

	maxAbsColSum(): number {
		let result = 0
		let colIndex = this.width
		while (colIndex--) {
			let absSum = 0
			let rowIndex = this.height
			while (rowIndex--) {
				absSum += Math.abs(this.m[rowIndex * this.width + colIndex])
			}
			result = Math.max(result, absSum)
		}
		return result
	}

	maxAbsRowSum(): number {
		let result = 0
		let rowIndex = this.height
		while (rowIndex--) {
			let absSum = 0
			let colIndex = this.width
			while (colIndex--) {
				absSum += Math.abs(this.m[rowIndex * this.width + colIndex])
			}
			result = Math.max(result, absSum)
		}
		return result
	}

	getTriangularDeterminant(): number {
		assert(this.isUpperTriangular() || this.isLowerTriangular(), "not a triangular matrix")

		let product = 1
		let elIndex = this.width * (this.width + 1)
		while (elIndex) {
			elIndex -= this.width + 1
			product *= this.m[elIndex]
		}
		return product
	}

	/**
	 * Calculates the determinant by first calculating the LU decomposition. If you already have that, use
	 * U.getTriangularDeterminant()
	 */
	getDeterminant(): number {
		// PA = LU
		// det(A) * det(B) = det(A * B)
		// det(P) == 1 (permutation matrix)
		// det(L) == 1 (main diagonal is 1s
		// =>  det(A) == det(U)
		return this.luDecomposition().U.getTriangularDeterminant()
	}

	hasFullRank(): boolean {
		return Math.min(this.width, this.height) == this.rank()
	}

	permutationAsIndexMap():int[] {
		assertf(() => this.isPermutation())
		let result = new Array(this.height)
		let i = this.height
		while (i--) {
			let searchIndexStart = i * this.width
			let searchIndex = searchIndexStart
			while (this.m[searchIndex] < 0.5) searchIndex++
			result[i] = searchIndex - searchIndexStart
		}
		return result
	}

	getDependentRowIndexes(gauss: {L: Matrix, U: Matrix, P: Matrix} = this.gauss()): int[] {
		let {L, U, P} = gauss
		let dependents = new Array(this.height)
		let uRowIndex = this.height
		while (uRowIndex--) {
			let uRow = U.row(uRowIndex)
			if (uRow.length() < NLA_PRECISION) {
				dependents[uRowIndex] = true
			} else {
				break
			}
		}
		let lRowIndex = this.height
		while (lRowIndex--) {
			if (dependents[lRowIndex]) {
				let lColIndex = Math.min(lRowIndex, this.width)
				while (lColIndex--) {
					if (0 !== L.e(lRowIndex, lColIndex)) {
						dependents[lColIndex] = true
					}
				}
			}
		}
		console.log("m\n", this.toString(x => '' + x))
		console.log("L\n", L.toString(x => '' + x))
		console.log("U\n", U.toString(x => '' + x))
		console.log("P\n", P.toString(x => '' + x))
		let indexMap = P.permutationAsIndexMap()
		let dependentRowIndexes = dependents.map((b, index) => b && indexMap[index]).filter(x => x != void 0)
		return dependentRowIndexes
	}


	/**
	 * Numerically calculate all the partial derivates of f at x0.
	 *
	 *
	 * @param f
	 * @param x0
	 * @param fx0 f(x0), pass it if you have it already
	 * @param EPSILON
	 */
	static jacobi(f: (x: number[]) => number[], x0: number[], fx0: number[], EPSILON?: number): Matrix {
		EPSILON = EPSILON || 1e-6
		fx0 = fx0 || f(x0)
		let jacobi = Matrix.forWidthHeight(x0.length, fx0.length)
		for (let colIndex = 0; colIndex < x0.length; colIndex++) {
			x0[colIndex] += EPSILON
			let fx = f(x0)
			for (let rowIndex = 0; rowIndex < fx0.length; rowIndex++) {
				jacobi.setEl(rowIndex, colIndex, (fx[rowIndex] - fx0[rowIndex]) / EPSILON)
			}
			x0[colIndex] -= EPSILON
		}
		return jacobi
	}
}