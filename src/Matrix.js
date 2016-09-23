
/**
 *
 * @param width
 * @param height
 * @param m
 * @constructor
 * @property {Float64Array} m
 */
NLA.Matrix = function (width, height, m) {
	assert(width * height == m.length, "width * height == m.length", width, height, m.length)
	this.m = m;
	this.width = width;
	this.height = height;
}
NLA.Matrix.random = function (width, height) {
	NLA.assertNumbers(width, height);
	return NLA.Matrix.fromFunction(width, height, (i, j) => Math.random())
}
NLA.Matrix.fromFunction = function(width, height, f) {
	NLA.assertNumbers(width, height);
	var m = new Float64Array(height * width)
	var elIndex = height * width;
	while (elIndex--) {
		m[elIndex] = f(Math.floor(elIndex / width), elIndex % width, elIndex);
	}
	return new NLA.Matrix(width, height, m);
}
NLA.Matrix.identity = function(dim) {
	NLA.assertNumbers(dim);
	var m = new Float64Array(dim * dim)
	// Float64Arrays are init to 0
	var elIndex = dim * (dim + 1);
	while (elIndex) {
		elIndex -= (dim + 1)
		m[elIndex] = 1
	}
	return new NLA.Matrix(dim, dim, m);
}
NLA.Matrix.permutation = function(dim, i, k) {
	NLA.assertNumbers(dim, i, k);
	var m = new Float64Array(dim * dim)
	// Float64Array are init to 0
	var elIndex = dim * (dim + 1);
	while (elIndex) {
		elIndex -= (dim + 1)
		m[elIndex] = 1
	}
	m[i * dim + i] = 0
	m[k * dim + k] = 0
	m[i * dim + k] = 1
	m[k * dim + i] = 1
	return new NLA.Matrix(dim, dim, m);
}
NLA.Matrix.fromRowArrays = function() {
	return NLA.Matrix.fromRowArrays2(arguments)
}
NLA.Matrix.fromRowArrays2 = function(arrays) {
	if (0 == arrays.length) {
		throw new Error("cannot have 0 vector");
	}
	var height = arrays.length
	var width = arrays[0].length
	var m = new Float64Array(height * width)
	NLA.arrayCopy(arrays[0], 0, m, 0, width);
	for (var rowIndex = 1; rowIndex < height; rowIndex++) {
		if (arrays[rowIndex].length != width) {
			throw new Error("all row arrays must be the same length")
		}
		NLA.arrayCopy(arrays[rowIndex], 0, m, rowIndex * width, width)
	}
	return new NLA.Matrix(width, height, m);
}
NLA.Matrix.fromColVectors = function (colVectors) {
	return NLA.Matrix.fromColArrays(colVectors.map((v) => v.v))
}
NLA.Matrix.forWidthHeight = function(width, height) {
	return new NLA.Matrix(width, height, new Float64Array(width * height))
}
NLA.Matrix.fromColArrays = function(colArrays) {
	if (0 == colArrays.length) {
		throw new Error("cannot have 0 vector");
	}
	var width = colArrays.length
	var height = colArrays[0].length
	var m = new Float64Array(height * width)
	NLA.arrayCopyStep(colArrays[0], 0, 1, m, 0, width, height)
	for (var colIndex = 1; colIndex < width; colIndex++) {
		if (colArrays[colIndex].length != height) {
			throw new Error("all col arrays must be the same length")
		}
		NLA.arrayCopyStep(colArrays[colIndex], 0, 1, m, colIndex, width, height)
	}
	return new NLA.Matrix(width, height, m);
}
NLA.Matrix.prototype = {
	e: function (rowIndex, colIndex) {
		NLA.assertNumbers(rowIndex, colIndex)
		if (NLA_DEBUG && (rowIndex >= this.height || colIndex >= this.width)) {
			throw new Error("index " + rowIndex + ", " + colIndex + " is out of bounds (" + this.width + " x " + this.height + ")")
		}
		return this.m[rowIndex * this.width + colIndex]
	},
	setEl: function (rowIndex, colIndex, val) {
		NLA.assertNumbers(rowIndex, colIndex, val)
		assert(0 <= rowIndex && rowIndex < this.height, "rowIndex out of bounds " + rowIndex)
		assert(0 <= colIndex && colIndex < this.width, "colIndex out of bounds " + colIndex)
		this.m[rowIndex * this.width + colIndex] = val
	},
	toString: function (f) {
		f = f || ((v) => v.toFixed(6))
		assert(typeof f(0) == "string", "" + typeof f(0))
		var rounded = Array.prototype.slice.call(this.m).map(f);
		var colWidths = NLA.arrayFromFunction(this.width,
			(colIndex) => rounded.sliceStep(colIndex, this.width).map((x) => x.length).max())
		return NLA.arrayFromFunction(this.height,
			(rowIndex) => rounded.slice(rowIndex * this.width, (rowIndex + 1) * this.width) // select matrix row
				.map((x, colIndex) => NLA.repeatString(colWidths[colIndex] - x.length, ' ') + x) // pad numbers with spaces to col width
				.join(" ")
		).map(x => x + "\n").join(""); // join rows
	},

	row: function (rowIndex) {
		var v = new Float64Array(this.width)
		NLA.arrayCopy(this.m, rowIndex * this.width, v, 0, this.width)
		return new NLA.Vector(v)
	},
	col: function (colIndex) {
		var v = new Float64Array(this.height)
		NLA.arrayCopyStep(this.m, colIndex, this.width, v, 0, 1, this.height)
		return new NLA.Vector(v)
	},
	dim: function() {
		return {width: this.width, height: this.height};
	},
	dimString: function () {
		return this.width + "x" + this.height
	},
	equals: function (obj) {
		if (obj.constructor != NLA.Matrix) return false;
		if (this.width != obj.width || this.height != obj.height) return false;
		var elIndex = this.m.length;
		while (elIndex--) {
			if (this.m[elIndex] != obj.m[elIndex]) return false;
		}
	},
	equalsMatrix: function (matrix, precision) {
		precision = precision || NLA_PRECISION
		if (!(matrix instanceof NLA.Matrix)) throw new Error("not a matrix");
		if (this.width != matrix.width || this.height != matrix.height) return false;
		var elIndex = this.m.length;
		while (elIndex--) {
			if (Math.abs(this.m[elIndex] - matrix.m[elIndex]) >= precision) return false;
		}
		return true
	},
	hashCode: function () {
		var result = 0
		var elIndex = this.m.length;
		while (elIndex--) {
			result = result * 31 + NLA.floatHashCode(this.m[elIndex])
		}
		return result
	},

	/**
	 *
	 * @returns {boolean}
	 */
	isZero() {
		var elIndex = this.m.length;
		while (elIndex--) {
			if (!NLA.isZero(this.m[elIndex])) {
				return false
			}
		}
		return true
	},
	isOrthogonal: function () {
		return this.isSquare() && this.transposed().times(this).equalsMatrix(NLA.Matrix.identity(this.width))
	},
	luDecomposition: function() {
		assertf(() => this.isSquare(), this.dim().toSource())
		var uRowArrays = this.asRowArrays(Float64Array)
		var dim = this.width
		var lRowArrays = NLA.arrayFromFunction(dim, (row) => new Float64Array(dim))
		var pRowArrays = NLA.Matrix.identity(dim).asRowArrays(Float64Array)
		let currentRowIndex = 0
		for (let colIndex = 0; colIndex < dim; colIndex++) {
			// find largest value in colIndex
			let maxAbsValue = 0, pivotRowIndex = undefined, numberOfNonZeroRows = 0
			for (let rowIndex = currentRowIndex; rowIndex < dim; rowIndex++) {
				var el = uRowArrays[rowIndex][colIndex]
				numberOfNonZeroRows += (0 != el) | 0
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
		return {L: NLA.Matrix.fromRowArrays2(lRowArrays), U: NLA.Matrix.fromRowArrays2(uRowArrays), P: NLA.Matrix.fromRowArrays2(pRowArrays)}
	},
	gauss: function() {
		var uRowArrays = this.asRowArrays(Float64Array)
		var w = this.width, h = this.height
		var lRowArrays = NLA.arrayFromFunction(h, (row) => new Float64Array(w))
		var pRowArrays = NLA.Matrix.identity(h).asRowArrays(Float64Array)
		let currentRowIndex = 0
		for (let colIndex = 0; colIndex < w; colIndex++) {
			// console.log('currentRowIndex', currentRowIndex)
			// find largest value in colIndex
			let maxAbsValue = 0, pivotRowIndex = undefined, numberOfNonZeroRows = 0
			for (let rowIndex = currentRowIndex; rowIndex < h; rowIndex++) {
				var el = uRowArrays[rowIndex][colIndex]
				numberOfNonZeroRows += (0 != el) | 0
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
				for (let rowIndex = currentRowIndex + 1; rowIndex < h; rowIndex++) {
					let l = uRowArrays[rowIndex][colIndex] / uRowArrays[currentRowIndex][colIndex]
					lRowArrays[rowIndex][colIndex] = l
					// subtract pivot row * l from row "rowIndex"
					for (let colIndex2 = colIndex; colIndex2 < w; colIndex2++) {
						uRowArrays[rowIndex][colIndex2] -= l * uRowArrays[currentRowIndex][colIndex2]
					}
				}
			}
			currentRowIndex++ // this doesn't increase if pivot was zero
		}
		return {L: NLA.Matrix.fromRowArrays2(lRowArrays), U: NLA.Matrix.fromRowArrays2(uRowArrays), P: NLA.Matrix.fromRowArrays2(pRowArrays)}
	},
	qrDecompositionGivensRotation: function () {
		function sigma (c, s) {
			if ( 0 == c) {
				return 1
			}
			if (Math.abs(s) < Math.abs(c)) {
				return 0.5 * Math.sign(c) * s
			}
			return 2 * Math.sign (s) / c;
		}
		function matrixForCS(dim, i, k, c, s) {
			var m = NLA.Matrix.identity(dim)
			m.setEl(i, i, c)
			m.setEl(k, k, c)
			m.setEl(i, k, s)
			m.setEl(k, i, -s)
			return m
		}

		var qTransposed = NLA.Matrix.identity(this.height)
		for (var colIndex = 0; colIndex < this.width; colIndex++) {
			// find largest value in colIndex
			for (var rowIndex = colIndex + 1; rowIndex < this.height; rowIndex++) {
				//console.log("row ", rowIndex, "col ", colIndex);
				var xi = this.e(colIndex, colIndex)
				var xk = this.e(rowIndex, colIndex)
				if (xk == 0) {
					continue;
				}
				var r = Math.sqrt(xi * xi + xk * xk);
				var c = xi / r
				var s = xk / r

				// apply transformation on every column:
				for (var col2 = colIndex; col2 < this.width; col2++) {
					var x1 = this.e(colIndex, col2) * c + this.e(rowIndex, col2) * s
					var x2 = this.e(rowIndex, col2) * c - this.e(colIndex, col2) * s
					this.setEl(colIndex, col2, x1)
					this.setEl(rowIndex, col2, x2)
				}
				//console.log("r ", r, "c ", c, "s ", s, "sigma", sigma(c, s));
				//console.log(this.toString(),"cs\n", matrixForCS(this.height, colIndex, rowIndex, c, s).toString())
				qTransposed = matrixForCS(this.height, colIndex, rowIndex, c, s).times(qTransposed)
			}
		}
		//console.log(qTransposed.transposed().toString(), this.toString(), qTransposed.transposed().times(this).toString())
		return {Q: qTransposed.transposed(), R: this}
	},
	isPermutation: function () {
		if (!this.isSquare()) return false
		if (this.m.some((value) => !NLA.isZero(value) && !NLA.equals(1, value))) return false

		var rows = this.asRowArrays(Array)
		if (rows.some((row) => row.filter((value) => NLA.equals(1, value)).length != 1)) return false

		var cols = this.asColArrays(Array)
		if (cols.some((col) => col.filter((value) => NLA.equals(1, value)).length != 1)) return false

		return true
	},
	isIdentity: function (precision) {
		return this.isLowerUnitriangular(precision) && this.isUpperTriangular(precision)
	},
	isUpperTriangular: function () {
		if (!this.isSquare()) return false
		for (var rowIndex = 1; rowIndex < this.height; rowIndex++) {
			for (var colIndex = 0; colIndex < rowIndex; colIndex++) {
				if (!NLA.isZero(this.m[rowIndex * this.width + colIndex])) {
					return false
				}
			}
		}
		return true
	},

	/**
	 * Returns x, so that this * x = b
	 * More efficient than calculating the inverse for few (~ <= this.height) values
	 * @param {NLA.Vector} b
	 * @returns {NLA.Vector}
	 */
	solveLinearSystem: function (b) {
		var lup = this.luDecomposition()
		// console.log(lup.L.toString())
		// console.log(lup.U.toString())
		// console.log(lup.P.toString())
		var y = lup.L.solveForwards(lup.P.timesVector(b))
		var x = lup.U.solveBackwards(y)
		return x
	},

	/**
	 *
	 * @param {number=} precision
	 * @returns {boolean}
	 */
	isLowerUnitriangular: function (precision) {
		precision = "number" == typeof precision ? precision : NLA_PRECISION
		if (!this.isSquare()) return false
		for (var rowIndex = 0; rowIndex < this.height - 1; rowIndex++) {
			for (var colIndex = rowIndex; colIndex < this.width; colIndex++) {
				var el = this.m[rowIndex * this.width + colIndex];
				if (rowIndex == colIndex ? !NLA.equals2(1, el, precision) : !NLA.isZero2(el, precision)) {
					return false
				}
			}
		}
		return true
	},

	/**
	 *
	 * @returns {boolean}
	 */
	isLowerTriangular: function () {
		if (!this.isSquare()) return false
		for (var rowIndex = 0; rowIndex < this.height - 1; rowIndex++) {
			for (var colIndex = rowIndex + 1; colIndex < this.width; colIndex++) {
				if (!NLA.isZero(this.m[rowIndex * this.width + colIndex])) {
					return false
				}
			}
		}
		return true
	},


	/**
	 *
	 * @param x
	 * @returns {NLA.Vector}
	 */
	solveBackwards: function (x) {
		NLA.assertVectors(x)
		assert(this.height == x.dim(), "this.height == x.dim()")
		assert(this.isUpperTriangular(), "this.isUpperTriangular()\n" + this.str)
		var v = new Float64Array(this.width)
		var rowIndex = this.height
		while (rowIndex--) {
			var temp = x.v[rowIndex]
			for (var colIndex = rowIndex + 1; colIndex < this.width; colIndex++) {
				temp -= v[colIndex] * this.e(rowIndex, colIndex)
			}
			v[rowIndex] = temp / this.e(rowIndex, rowIndex)
		}
		return new NLA.Vector(v)
	},
	solveBackwardsMatrix: function (matrix) {
		var colVectors = new Array(matrix.width)
		var i = matrix.width
		while (i--) {
			colVectors[i] = this.solveBackwards(matrix.col(i))
		}
		return NLA.Matrix.fromColVectors(colVectors)
	},
	solveForwardsMatrix: function (matrix) {
		var colVectors = new Array(matrix.width)
		var i = matrix.width
		while (i--) {
			colVectors[i] = this.solveForwards(matrix.col(i))
		}
		return NLA.Matrix.fromColVectors(colVectors)
	},
	solveForwards: function (x) {
		NLA.assertVectors(x)
		assert(this.height == x.dim(), "°°")
		assertf(() => this.isLowerTriangular(), this.toString())
		var v = new Float64Array(this.width)
		for (var rowIndex = 0; rowIndex < this.height; rowIndex++) {
			var temp = x.v[rowIndex]
			for (var colIndex = 0; colIndex < rowIndex; colIndex++) {
				temp -= v[colIndex] * this.e(rowIndex, colIndex)
			}
			v[rowIndex] = temp / this.e(rowIndex, rowIndex)
		}
		return new NLA.Vector(v)
	},



	/**
	 * Calculates rank of matrix.
	 * Number of linearly independant row/column vectors.
	 * Is equal to the unmber of dimensions the image of the affine transformation represented this matrix has.
	 *
	 * @returns {number} integer
	 */
	rank: function () {
		let U = this.gauss().U
		//console.log(R.toString())
		var rowIndex = this.height
		while (rowIndex-- && U.row(rowIndex).isZero()) {
			console.log("RANK" + U.row(rowIndex).toString() + U.row(rowIndex).isZero())}
		return rowIndex + 1
	},
	rowsIndependent: function () {
		return this.height == this.rank()
	},
	colsIndependent: function () {
		return this.width == this.rank()
	},
	asRowArrays: function (arrayConstructor) {
		arrayConstructor = arrayConstructor || Float64Array
		var rowIndex = this.height
		var result = new Array(this.height)
		while (rowIndex--) {
			result[rowIndex] = this.rowArray(rowIndex, arrayConstructor)
		}
		return result
	},
	asColArrays: function (arrayConstructor) {
		arrayConstructor = arrayConstructor || Float64Array
		var colIndex = this.width
		var result = new Array(this.width)
		while (colIndex--) {
			result[colIndex] = this.colArray(colIndex, arrayConstructor)
		}
		return result
	},
	rowArray: function (rowIndex, arrayConstructor) {
		arrayConstructor = arrayConstructor || Float64Array
		var result = new arrayConstructor(this.width)
		NLA.arrayCopy(this.m, rowIndex * this.width, result, 0, this.width)
		return result
	},
	colArray: function (colIndex, arrayConstructor) {
		arrayConstructor = arrayConstructor || Float64Array
		var result = new arrayConstructor(this.width)
		NLA.arrayCopyStep(this.m, colIndex, this.height, result, 0, 1, this.height)
		return result
	},
	subMatrix: function (firstColIndex, subWidth, firstRowIndex, subHeight) {
		if (firstColIndex + subWidth > this.width || firstRowIndex + subHeight > this.height) {
			throw new Error("inavlid params")
		}
		var m = new Float64Array(this.height)
		NLA.arrayCopyBlocks(this.m, firstColIndex, this.width, m, 0, subWidth, subHeight, subWidth)
		return new NLA.Matrix(subWidth, subHeight, m)
	},
	map: function (fn) {
		return new NLA.Matrix(this.width, this.height, this.m.map(fn));
	},

	dimEquals: function (matrix) {
		assertInst(NLA.Matrix, matrix)
		return this.width == matrix.width && this.height == matrix.height
	},

	inversed: function () {
		var lup = this.luDecomposition()
		var y = lup.L.solveForwardsMatrix(lup.P)
		var inverse = lup.U.solveBackwardsMatrix(y)
		return inverse
	},

	inversed3: function () {
		assertf(() => 3 == this.width && 3 == this.height)
		let result = new NLA.Matrix.forWidthHeight(3, 3), m = this.m, r = result.m

		r[0] = m[4]*m[8] - m[5]*m[7]
		r[1] = -m[1]*m[8] + m[2]*m[7]
		r[2] = m[1]*m[5] - m[2]*m[4]

		r[3] = -m[3]*m[8] + m[5]*m[6]
		r[4] = m[0]*m[8] - m[2]*m[6]
		r[5] = -m[0]*m[5] + m[2]*m[3]

		r[6] = m[3]*m[7] - m[4]*m[6]
		r[7] = -m[0]*m[7] + m[1]*m[6]
		r[8] = m[0]*m[4] - m[1]*m[3]

		let det = m[0]*r[0] + m[1]*r[3] + m[2]*r[6]
		var i = 9
		while (i--) { r[i] /= det }

		return result
	},

	inversed2: function () {
		assertf(() => 2 == this.width && 2 == this.height)
		let result = new NLA.Matrix.forWidthHeight(2, 2), m = this.m, r = result.m

		let det = m[0]*m[3] - m[1]*r[2]

		r[0] = m[3] / det
		r[1] = -m[2] / det

		r[2] = -m[1] / det
		r[3] = m[0] / det

		return result
	},

	canMultiply: function (matrix) {
		assertInst(NLA.Matrix, matrix)
		return this.width == matrix.height
	},
	times: function (matrix) {
		assertInst(NLA.Matrix, matrix)
		assert(this.canMultiply(matrix), `Cannot multiply this {this.dimString()} by matrix {matrix.dimString()}`)
		var nWidth = matrix.width, nHeight = this.height, n = this.width
		var nM = new Float64Array(nWidth * nHeight)
		var nRowIndex = nHeight
		while (nRowIndex--) {
			var nColIndex = nWidth
			while (nColIndex--) {
				var result = 0
				var i = n
				while (i--) {
					result += this.m[nRowIndex * n + i] * matrix.m[i * nWidth + nColIndex]
				}
				nM[nRowIndex * nWidth + nColIndex] = result
			}
		}
		return new NLA.Matrix(nWidth, nHeight, nM)
	},
	timesVector: function (vector) {
		NLA.assertVectors(vector)
		assert(this.width == vector.dim())
		var nHeight = this.height, n = this.width
		var nM = new Float64Array(nHeight)
		var nRowIndex = nHeight
		while (nRowIndex--) {
			var result = 0
			var i = n
			while (i--) {
				result += this.m[nRowIndex * n + i] * vector.v[i]
			}
			nM[nRowIndex] = result
		}
		return new NLA.Vector(nM)
	},
	transposed: function () {
		var tWidth = this.height, tHeight = this.width
		var tM = new Float64Array(tWidth * tHeight)
		var tRowIndex = tHeight
		while (tRowIndex--) {
			var tColIndex = tWidth
			while (tColIndex--) {
				tM[tRowIndex * tWidth + tColIndex] = this.m[tColIndex * tHeight + tRowIndex]
			}
		}
		return new NLA.Matrix(tWidth, tHeight, tM)
	},
	transpose: function () {
		var h = this.height, w = this.width, tM = this.m
		var tRowIndex = h
		while (tRowIndex--) {
			var tColIndex = Math.min(tRowIndex, w)
			while (tColIndex--) {
				console.log("col", tColIndex, "row", tRowIndex)
				var temp = tM[tRowIndex * w + tColIndex]
				tM[tRowIndex * w + tColIndex] = tM[tColIndex * h + tRowIndex]
				tM[tColIndex * h + tRowIndex] = temp
			}
		}
		this.width = h
		this.height = w
	},
	isSquare: function () {
		return this.height == this.width;
	},
	diagonal: function () {
		if (!this.isSquare()) {
			throw new Error("!!")
		}
		var v = new Float64Array(this.width)
		var elIndex = this.width * (this.width + 1)
		var vIndex = this.width
		while (vIndex--) {
			elIndex -= this.width + 1
			v[vIndex] = this.m[elIndex]
		}
		return new NLA.Vector(v)
	},
	max: function () {
		return NLA.arrayMax(this.m);
	},
	min: function () {
		return NLA.arrayMin(this.m);
	},
	maxAbsColSum: function () {
		var result = 0
		var colIndex = this.width
		while (colIndex--) {
			var absSum = 0
			var rowIndex = this.height
			while (rowIndex--) {
				absSum  += Math.abs(this.m(rowIndex * this.width + colIndex))
			}
			result = Math.max(result, absSum)
		}
		return result
	},
	maxAbsRowSum: function () {
		var result = 0
		var rowIndex = this.height
		while (rowIndex--) {
			var absSum = 0
			var colIndex = this.width
			while (colIndex--) {
				absSum += Math.abs(this.m(rowIndex * this.width + colIndex))
			}
			result = Math.max(result, absSum)
		}
		return result
	},
	getTriangularDeterminant: function () {
		assert(this.isUpperTriangular() || this.isLowerTriangular(), "not a triangular matrix")

		var product = 1
		var elIndex = this.width * (this.width + 1)
		while (elIndex) {
			elIndex -= this.width + 1
			product *= this.m[elIndex]
		}
		return product
	},
	/**
	 * Calculates the determinant by first calculating the LU decomposition. If you already have that, use
	 * U.getTriangularDeterminant()
	 * @returns {*}
	 */
	getDeterminant: function () {
		// PA = LU
		// det(A) * det(B) = det(AB)
		// det(P) == 1 (permutation matrix)
		// det(L) == 1 (main diagonal is 1s
		// =>  det(A) == det(U)
		return this.luDecomposition().U.getTriangularDeterminant()
	},
	hasFullRank: function () {
		return Math.min(this.width, this.height) == this.rank()
	}
}

/**
 * Numerically calculate all the partial derivates of f at x0.
 *
 *
 * @param {function(Array.<number>):Array.<number>} f
 * @param {Array.<number>} x0
 * @param {Array.<number>=} fx0 f(x0), pass it if you have it already
 * @param {number=} EPSILON
 */
NLA.Matrix.jacobi = function (f, x0, fx0, EPSILON) {
	EPSILON = EPSILON || 1e-8
	fx0 = fx0 || f(x0)
	let jacobi = NLA.Matrix.forWidthHeight(x0.length, fx0.length)
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