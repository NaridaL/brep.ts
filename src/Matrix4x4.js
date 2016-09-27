/**
* Created by aval on 20/11/2015.
*/
"use strict";
if (!V3) {0
	throw new Error("Define NLA.V3 first")
}

/**
 * This constructor takes 16 arguments in row-major order, which can be passed individually, as a list, or even as
 * four lists, one for each row. If the arguments are omitted then the identity matrix is constructed instead.
 *
 * @returns {M4}
 * @constructor
 * @property {Float64Array} m
 * @param {...number|number[]} var_args
 */
function M4(...var_args) {
    let m
    if (0 == arguments.length) {
        m = new Float64Array(16)
    } else {
        var flattened = Array.prototype.concat.apply([], arguments)
        assert(flattened.length == 16, "flattened.length == 16")
        //noinspection JSCheckFunctionSignatures
        m = new Float64Array(flattened)
    }
    let o = Object.create(M4.prototype)
    Object.defineProperty(o, 'm', {value: m})
    return o
}

M4.prototype = NLA.defineObject(NLA.Matrix.prototype, /** @lends {M4.prototype} */ {
    constructor: M4,
    height: 4,
    width: 4,

	/**
	 * Returns a new M4 which is equal to the inverse of this.
	 *
	 * @returns {M4}
	 */
    inversed: function() {
        return M4.inverse(this);
    },

	/**
	 * Matrix trace is defined as the sum of the elements of the main diagonal.
	 *
	 * @returns {number}
	 */
    trace: function () {
	    return this.m[0] + this.m[5] + this.m[10] + this.m[15]
    },

    /**
     * @returns {M4}
     */
    as3x3: function () {
		let result = M4.copy(this), m = result.m

	    m[3] = m[7] = m[11] = m[12] = m[13] = m[14] = 0
	    m[15] = 1
	    return result
    },

    realEigenValues3: function () {
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
	    let u = -a*e - a*i + b*d + c*g - e*i + f*h // equivalent to 1/2 (trace(this²) - trace²(A))
	    let w = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g) // equivalent to matrix determinant

	    console.log(s, t, u, w)
	    return solveCubicReal2(s, t, u, w)

    },

    realEigenVectors3: function () {
        let eigenValues = this.realEigenValues3()
	    let this3x3 = this.times(M4.IDENTITY3)
	    console.log(this.toString())
	    console.log(this3x3.toString())
	    let mats = eigenValues.map(ev => M4.IDENTITY3.timesScalar(-ev).plus(this3x3))
	    console.log(mats.map(m=>m.determinant3()))
	    console.log(mats.map(m=>''+m.toString(v=>''+v)).join("\n\n"))
	    console.log(mats.map(m=>''+m.gauss().U.toString(v=>''+v)).join("\n\n"))
	    console.log("mats.map(m=>m.rank())", mats.map(m=>m.rank()))
	    if (1 == eigenValues.length) {
	        console.log(mats[0].toString())
            assertf(() => 0 == mats[0].rank())
		    // col vectors
		    return NLA.arrayFromFunction(3, col => V3.create(this.m[col], this.m[4 + col], this.m[8 + col]))
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
		    mats.forEach((mat, i) => assert(2 == mat.rank(), i+': '+mat.rank()))
		    // the (A - lambda I) matrices map to a plane. This means, that there is an entire line in R³ which maps to the point V3.ZERO
			return mats.map(mat => {
				let gauss = mat.gauss().U
				return gauss.row(0).cross(gauss.row(1)).V3().normalized()
			})
	    }
    },

	/**
	 *
	 * @param m4
	 * @returns {boolean}
	 */
    likeM4: function (m4) {
	    assertInst(M4, m4)
	    return this.m.every((el, index) => NLA.equals(el, m4.m[index]))
    },

	/**
	 * Returns a new M4 equal to the transpose of this.
	 * @returns {M4}
	 */
    transposed: function() {
        return M4.transpose(this);
    },

	/**
	 * Returns a new M4 which equal to (this * matrix) (in that order)
	 * @param {M4} matrix
	 * @returns {M4}
	 */
    times: function(matrix) {
        return M4.multiply(this, matrix);
    },

	/**
	 * Transforms the vector as a point with a w coordinate of 1. This means translations will have an effect, for
	 * example.
	 *
	 * @param {V3} v
	 * @returns {V3}
	 */
    transformPoint: function(v) {
	    assertVectors(v)
	    var m = this.m;
	    var vx = v.x, vy = v.y, vz = v.z, vw = 1
	    var x = vx * m[ 0] + vy * m[ 1] + vz * m[ 2] + vw * m[ 3]
	    var y = vx * m[ 4] + vy * m[ 5] + vz * m[ 6] + vw * m[ 7]
	    var z = vx * m[ 8] + vy * m[ 9] + vz * m[10] + vw * m[11]
	    var w = vx * m[12] + vy * m[13] + vz * m[14] + vw * m[15]
	    // scale such that fourth element becomes 1:
	    return V3.create(x / w, y / w, z / w);
    },

	/**
	 * Transforms the vector as a vector with a w coordinate of 0. This means translations will have no effect, for
	 * example. Will throw an exception if the calculated w component != 0. This occurs for example when attempting
	 * to transform a vector with a perspective matrix.
	 *
	 * @param {V3} v
	 * @returns {V3}
	 */
    transformVector: function(v) {
        assertVectors(v)
        var m = this.m;
        var w = v.x * m[12] + v.y * m[13] + v.z * m[14]
        assert(w == 0, 'w != 0 needs to be true for this to make sense (w ='+w)
        return V3.create(
            m[0] * v.x + m[1] * v.y + m[2] * v.z,
            m[4] * v.x + m[5] * v.y + m[6] * v.z,
            m[8] * v.x + m[9] * v.y + m[10] * v.z
        );
    },

	/**
	 *
	 * @param {Array<V3>} vs
	 * @returns {Array<V3>}
	 */
    transformedPoints: function(vs) {
	    return vs.map(v => this.transformPoint(v))
    },

	/**
	 *
	 * @param {Array<V3>} vs
	 * @returns {Array<V3>}
	 */
    transformedVectors: function(vs) {
	    return vs.map(v => this.transformVector(v))
    },

	/**
	 *
	 * @param {M4} m
	 * @returns {M4}
	 */
    plus: function(m) {
        var r = M4()
        for (var i = 0; i < 16; i++) {
            r.m[i] = this.m[i] + m.m[i];
        }
        return r;
    },

	/**
	 *
	 * @param {M4} m
	 * @returns {M4}
	 */
    minus: function(m) {
        var r = M4()
        for (var i = 0; i < 16; i++) {
            r.m[i] = this.m[i] - m.m[i];
        }
        return r;
    },

	/**
	 *
	 * @param {number} scalar
	 * @returns {M4}
	 */
    timesScalar: function (scalar) {
	    var r = M4()
	    for (var i = 0; i < 16; i++) {
		    r.m[i] = this.m[i] * scalar
	    }
	    return r;
    },

	/**
	 *
	 * @param {number} scalar
	 * @returns {M4}
	 */
    divScalar: function (scalar) {
	    var r = M4()
	    for (var i = 0; i < 16; i++) {
		    r.m[i] = this.m[i] / scalar
	    }
	    return r;
    },

    /**
     *
     * @returns {M4}
     */
    clone: function() {
	    return M4.copy(this)
    },

    /**
     *
     * @returns {M4}
     */
    copy: function() {
	    return M4.copy(this)
    },

	/**
	 *
	 * @returns {boolean}
	 */
    isRegular: function () {
	    return !NLA.isZero(this.determinant())
    },


    /**
     *
     * @returns {boolean}
     */
	isAxisAligned: function () {
		var $ = this.m
		return (1 >= !NLA.isZero($[0]) + !NLA.isZero($[1]) + !NLA.isZero($[2]))
			&& (1 >= !NLA.isZero($[4]) + !NLA.isZero($[5]) + !NLA.isZero($[6]))
			&& (1 >= !NLA.isZero($[8]) + !NLA.isZero($[9]) + !NLA.isZero($[10]))
	},

	/**
	 * A matrix M is orthogonal iff M * M^T = I
	 * I being the identity matrix.
	 *
	 * @returns {boolean} If this matrix is orthogonal or very close to it. Comparison of the identity matrix and
	 * this * this^T is done with {@link #likeM4}
	 */
    isOrthogonal: function () {
		// return this.transposed().times(this).likeM4(M4.IDENTITY)
		M4.transpose(this, M4.temp0)
	    M4.multiply(this, M4.temp0, M4.temp1)
	    return M4.IDENTITY.likeM4(M4.temp1)
    },

	/**
	 * A matrix M is symmetric iff M == M^T
	 * I being the identity matrix.
	 *
	 * @returns {boolean} If this matrix is symmetric or very close to it. Comparison of the identity matrix and
	 * this * this^T is done with {@link #likeM4}
	 */
	isSymmetric: function () {
		M4.transpose(this, M4.temp0)
		return this.likeM4(M4.temp0)
	},

	/**
	 * A matrix M is normal iff M * M^-T == M^T * M
	 * I being the identity matrix.
	 *
	 * @returns {boolean} If this matrix is symmetric or very close to it. Comparison of the identity matrix and
	 * this * this^T is done with {@link #likeM4}
	 */
	isNormal: function () {
		M4.transpose(this, M4.temp0) // temp0 = this^-T
		M4.multiply(this, M4.temp0, M4.temp1) // temp1 = this * this^-T
		M4.multiply(M4.temp0, this, M4.temp2) // temp2 = this^-T * this
		return M4.temp1.likeM4(M4.temp2)
	},

	/**
     * Determinant of matrix.
     *
     * Notes:
     *      For matrices A and B
     *      det(AB) = det(A) * det(B)
     *      det(A^-1) = 1 / det(A)
     *
     * @returns {number}
     */
    determinant: function() {
        /*
         | a b c d |
         | e f g h |
         | i j k l |
         | m n o p |
        */
        var $ = this.m,
            a = $[ 0], b = $[ 1], c = $[ 2], d = $[ 3],
            e = $[ 4], f = $[ 5], g = $[ 6], h = $[ 7],
            i = $[ 8], j = $[ 9], k = $[10], l = $[11],
            m = $[12], n = $[13], o = $[14], p = $[15],
            klop = k * p - l * o, jlnp = j * p - l * n, jkno = j * o - k * n,
            ilmp = i * p - l * m, ikmo = i * o - k * m, ijmn = i * n - j * m;
        return (
              a * (f * klop - g * jlnp + h * jkno)
            - b * (e * klop - g * ilmp + h * ikmo)
            + c * (e * jlnp - f * ilmp + h * ijmn)
            - d * (e * jkno - f * ikmo + g * ijmn))
    },

	determinant3: function() {
		let [a, b, c, ,
			d, e, f, ,
			g, h, i] = this.m
		let det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)
		return det
	},

	/**
	 * determine whether this matrix is a mirroring transformation
	 * @returns {boolean}
	 */
    isMirroring: function() {
        /*
        var u = V3.create(this.m[0], this.m[4], this.m[8]);
        var v = V3.create(this.m[1], this.m[5], this.m[9]);
        var w = V3.create(this.m[2], this.m[6], this.m[10]);

        // for a true orthogonal, non-mirrored base, u.cross(v) == w
        // If they have an opposite direction then we are mirroring
        var mirrorvalue = u.cross(v).dot(w);
        var ismirror = (mirrorvalue < 0);
         return ismirror;
        */

        return this.determinant() < 0 // TODO: also valid for 4x4?

    },

    /**
     * Get the translation part of this matrix, i.e. the result of this.transformVector(V3.ZERO)
     * @returns {V3}
     */
    getTranslation: function () {
	    var m = this.m, w = m[15]
	    return V3.create(m[3] / w, m[7] / w, m[11] / w)
    },

    /**
     * Returns this matrix scaled so that the determinant is 1.
     * det(c * A) = (c ** n) * det(A) for n x n matrices,
     * so we need to divide by the 4th root of the determinant
     * @returns {M4}
     */
    normalized: function () {
	    var det = this.determinant()
		return 1 == det ? this : this.divScalar(Math.sqrt(Math.sqrt(det)))
    },
    /**
     * Returns this matrix scaled so that the determinant is 1.
     * det(c * A) = (c ** n) * det(A) for n x n matrices,
     * so we need to divide by the 4th root of the determinant
     * @returns {M4}
     */
    normalized2: function () {
	    let div = this.m[15]
		return 1 == div ? this : this.divScalar(Math.sqrt(Math.sqrt(div)))
    },

    /**
     * Returns if the matrix has the following form (within NLA_PRECISION):
     * a b c 0
     * c d e 0
     * f g h 0
     * 0 0 0 1
     * @returns {boolean}
     */
    is3x3: function () {
	    var m = this.m
	    return NLA.equals(1, m[15])
		    && NLA.isZero(m[12]) && NLA.isZero(m[13]) && NLA.isZero(m[14])
		    && NLA.isZero(m[3]) && NLA.isZero(m[7]) && NLA.isZero(m[11])
    },

    /**
     *
     * @returns {boolean}
     */
    isIdentity: function () {
	    return this.m.every((val, i) => (i / 4 | 0) == (i % 4) ? NLA.equals(1, val) : NLA.isZero(val))
    },

    /**
     *
     * @param {function(number):string=} f
     * @returns {string|*}
     */
    toString: function (f) {
        f = f || ((v) => v.toFixed(6).replace(/(0|\.)(?=0*$)/g, " ").toString())
	    assert(typeof f(0) == "string", "" + typeof f(0))
	    // slice this.m to convert it to an Array (from TypeArray)
	    var rounded = Array.prototype.slice.call(this.m).map(f);
	    var colWidths = [0, 1, 2, 3].map((colIndex) => rounded.sliceStep(colIndex, 4).map((x) => x.length).max());
	    return [0, 1, 2, 3].map(
		    (rowIndex) => rounded
			    .slice(rowIndex * 4, rowIndex * 4 + 4) // select matrix row
			    .map((x, colIndex) => NLA.repeatString(colWidths[colIndex] - x.length, ' ') + x) // pad numbers with spaces to col width
			    .join(" ")
	    ).join("\n"); // join rows
    }
})

/**
 * Returns the matrix that when multiplied with `matrix` results in the
 * identity matrix. You can optionally pass an existing matrix in `result`
 * to avoid allocating a new matrix. This implementation is from the Mesa
 * OpenGL function `__gluInvertMatrixd()` found in `project.c`.
 *
 * @param {M4} matrix
 * @param {M4=} result
 * @returns {M4}
 */
M4.inverse = function(matrix, result) {
    assertInst(M4, matrix)
    assert(!result || result instanceof M4, "!result || result instanceof M4")
    assert(matrix != result, "matrix != result")
    result = result || M4()
    var m = matrix.m, r = result.m

	// first compute transposed cofactor matrix:
	// cofactor of an element is the determinant of the 3x3 matrix gained by removing the column and row belonging to the element
    r[0] = m[5]*m[10]*m[15] - m[5]*m[14]*m[11] - m[6]*m[9]*m[15] + m[6]*m[13]*m[11] + m[7]*m[9]*m[14] - m[7]*m[13]*m[10]
    r[1] = -m[1]*m[10]*m[15] + m[1]*m[14]*m[11] + m[2]*m[9]*m[15] - m[2]*m[13]*m[11] - m[3]*m[9]*m[14] + m[3]*m[13]*m[10]
    r[2] = m[1]*m[6]*m[15] - m[1]*m[14]*m[7] - m[2]*m[5]*m[15] + m[2]*m[13]*m[7] + m[3]*m[5]*m[14] - m[3]*m[13]*m[6]
    r[3] = -m[1]*m[6]*m[11] + m[1]*m[10]*m[7] + m[2]*m[5]*m[11] - m[2]*m[9]*m[7] - m[3]*m[5]*m[10] + m[3]*m[9]*m[6]

    r[4] = -m[4]*m[10]*m[15] + m[4]*m[14]*m[11] + m[6]*m[8]*m[15] - m[6]*m[12]*m[11] - m[7]*m[8]*m[14] + m[7]*m[12]*m[10]
    r[5] = m[0]*m[10]*m[15] - m[0]*m[14]*m[11] - m[2]*m[8]*m[15] + m[2]*m[12]*m[11] + m[3]*m[8]*m[14] - m[3]*m[12]*m[10]
    r[6] = -m[0]*m[6]*m[15] + m[0]*m[14]*m[7] + m[2]*m[4]*m[15] - m[2]*m[12]*m[7] - m[3]*m[4]*m[14] + m[3]*m[12]*m[6]
    r[7] = m[0]*m[6]*m[11] - m[0]*m[10]*m[7] - m[2]*m[4]*m[11] + m[2]*m[8]*m[7] + m[3]*m[4]*m[10] - m[3]*m[8]*m[6]

    r[8] = m[4]*m[9]*m[15] - m[4]*m[13]*m[11] - m[5]*m[8]*m[15] + m[5]*m[12]*m[11] + m[7]*m[8]*m[13] - m[7]*m[12]*m[9]
    r[9] = -m[0]*m[9]*m[15] + m[0]*m[13]*m[11] + m[1]*m[8]*m[15] - m[1]*m[12]*m[11] - m[3]*m[8]*m[13] + m[3]*m[12]*m[9]
    r[10] = m[0]*m[5]*m[15] - m[0]*m[13]*m[7] - m[1]*m[4]*m[15] + m[1]*m[12]*m[7] + m[3]*m[4]*m[13] - m[3]*m[12]*m[5]
    r[11] = -m[0]*m[5]*m[11] + m[0]*m[9]*m[7] + m[1]*m[4]*m[11] - m[1]*m[8]*m[7] - m[3]*m[4]*m[9] + m[3]*m[8]*m[5]

    r[12] = -m[4]*m[9]*m[14] + m[4]*m[13]*m[10] + m[5]*m[8]*m[14] - m[5]*m[12]*m[10] - m[6]*m[8]*m[13] + m[6]*m[12]*m[9]
    r[13] = m[0]*m[9]*m[14] - m[0]*m[13]*m[10] - m[1]*m[8]*m[14] + m[1]*m[12]*m[10] + m[2]*m[8]*m[13] - m[2]*m[12]*m[9]
    r[14] = -m[0]*m[5]*m[14] + m[0]*m[13]*m[6] + m[1]*m[4]*m[14] - m[1]*m[12]*m[6] - m[2]*m[4]*m[13] + m[2]*m[12]*m[5]
    r[15] = m[0]*m[5]*m[10] - m[0]*m[9]*m[6] - m[1]*m[4]*m[10] + m[1]*m[8]*m[6] + m[2]*m[4]*m[9] - m[2]*m[8]*m[5]

	// calculate determinant using laplace expansion (cf https://en.wikipedia.org/wiki/Laplace_expansion),
	// as we already have the cofactors. We multiply a column by a row as the cofactor matrix is transposed.
    var det = m[0]*r[0] + m[1]*r[4] + m[2]*r[8] + m[3]*r[12]
   // assert(!NLA.isZero(det), "det may not be zero, i.e. the matrix is not invertible")
    var i = 16
    while (i--) { r[i] /= det }
    return result
}

/**
 * Returns `matrix`, exchanging columns for rows. You can optionally pass an
 * existing matrix in `result` to avoid allocating a new matrix.
 * @param {M4} matrix
 * @param {M4=} result
 * @returns {M4}
 */
M4.transpose = function(matrix, result) {
    assertInst(M4, matrix)
    assert(!result || result instanceof M4, "!result || result instanceof M4")
    assert(matrix != result, "matrix != result")
    result = result || M4();
    var m = matrix.m, r = result.m;
    r[ 0] = m[0]; r[ 1] = m[4]; r[ 2] = m[ 8]; r[ 3] = m[12];
    r[ 4] = m[1]; r[ 5] = m[5]; r[ 6] = m[ 9]; r[ 7] = m[13];
    r[ 8] = m[2]; r[ 9] = m[6]; r[10] = m[10]; r[11] = m[14];
    r[12] = m[3]; r[13] = m[7]; r[14] = m[11]; r[15] = m[15];
    return result;
};

/**
 * Returns the concatenation of the transforms for `left` and `right`.
 *
 * @param {M4} left
 * @param {M4} right
 * @param {M4=} result
 * @returns {M4}
 */
M4.multiply = function(left, right, result) {
    assertInst(M4, left, right)
    assert(!result || result instanceof M4, "!result || result instanceof M4")
    assert(left != result, "left != result")
    assert(right != result, "right != result")
    result = result || M4();
    var a = left.m, b = right.m, r = result.m;

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

    return result;
};

/**
 *
 * @param {M4} src
 * @param {M4=} result
 * @returns {M4}
 */
M4.copy = function(src, result) {
	assertInst(M4, src)
	assert(!result || result instanceof M4, "!result || result instanceof M4")
	assert(result != src, "result != src")
	result = result || M4();
	var s = src.m, d = result.m
	var i = 16
	while (i--) {
		d[i] = s[i]
	}
	return result
}

/**
 *
 * @param {V3} e0
 * @param {V3} e1
 * @param {V3} e2
 * @param {V3=} origin
 * @returns {M4}
 */
M4.forSys = function (e0, e1, e2, origin) {
    assertVectors(e0)
    assertVectors(e1)
    assert(!e2 || e2 instanceof V3, "!e2 || e2 instanceof V3")
    assert(!origin || origin instanceof V3, "!origin || origin instanceof V3")

    e2 = e2 || e0.cross(e1);
    origin = origin || V3.ZERO
    return M4(
        e0.x, e1.x, e2.x, origin.x,
        e0.y, e1.y, e2.y, origin.y,
        e0.z, e1.z, e2.z, origin.z,
        0,	     0,    0, 1);
}

/**
 *
 * @param {V3} n0
 * @param {V3} n1
 * @param {V3} n2
 * @param {V3=} n3
 * @returns {M4}
 */
M4.forRows = function (n0, n1, n2, n3) {
    assertVectors(n0, n1, n2)
	!n3 || assertVectors(n2)
	n3 = n3 || V3.ZERO
    return M4(
	    n0.x, n0.y, n0.z, 0,
	    n1.x, n1.y, n1.z, 0,
	    n2.x, n2.y, n2.z, 0,
	    n3.x, n3.y, n3.z, 1);
}

/**
 * Returns an identity matrix. You can optionally pass an existing matrix in `result` to avoid allocating a new
 * matrix. This emulates the OpenGL function `glLoadIdentity()`
 *
 * Unless initializing a matrix to be modified, use M4.IDENTITY
 * @param {M4=} result
 * @returns {M4}
 */
M4.identity = function(result) {
    assert(!result || result instanceof M4, "!result || result instanceof M4")
    result = result || M4();
    var m = result.m;
    m[0] = m[5] = m[10] = m[15] = 1;
    m[1] = m[2] = m[3] = m[4] = m[6] = m[7] = m[8] = m[9] = m[11] = m[12] = m[13] = m[14] = 0;

    return result;
};

/**
 * Creates a new M4 initialized by a user defined callback function
 *
 * @param {function(elRow: number, elCol: number, elIndex: number): number} f signature: (elRow, elCol, elIndex) =>
 *     el, where elIndex is the row-major index, i.e. elRow * 4 + elCol
 * @param {M4=} result
 * @returns {M4}
 */
M4.fromFunction = function (f, result) {
	assert(typeof f == "function", 'typeof f == "function"' + typeof f)
	assert(!result || result instanceof M4, "!result || result instanceof M4")
	result = result || M4()
	var m = result.m;
	var i = 16
	while (i--) {
		m[i] = f(Math.floor(i / 4), i % 4, i);
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
 * @param {number} fov in degrees
 * @param {number} aspect Aspect ration, width/height of viewport
 * @param {number} near
 * @param {number} far
 * @param {M4} [result]
 * @returns {M4}
 */
M4.perspective = function(fov, aspect, near, far, result) {
    assert(!result || result instanceof M4, "!result || result instanceof M4")
    assertNumbers(fov, aspect, near, far)
    var y = Math.tan(fov * Math.PI / 360) * near;
    var x = y * aspect;
    return M4.frustum(-x, x, -y, y, near, far, result);
};

// ### GL.Matrix.frustum(left, right, bottom, top, near, far[, result])
//
// Sets up a viewing frustum, which is shaped like a truncated pyramid with the
// camera where the point of the pyramid would be. You can optionally pass an
// existing matrix in `result` to avoid allocating a new matrix. This emulates
// the OpenGL function `glFrustum()`.
/**
 *
 * @param {number} l
 * @param {number} r
 * @param {number} b
 * @param {number} t
 * @param {number} n
 * @param {number} f
 * @param {M4} [result]
 * @returns {M4}
 */
M4.frustum = function(l, r, b, t, n, f, result) {
    assertNumbers(l, r, b, t, n, f)
    assert(0 < n, "0 < n")
    assert(n < f, "n < f")
    assert(!result || result instanceof M4, "!result || result instanceof M4")
    result = result || M4()
    var m = result.m

    m[0] = 2 * n / (r - l)
    m[1] = 0
    m[2] = (r + l) / (r - l)
    m[3] = 0

    m[4] = 0
    m[5] = 2 * n / (t - b)
    m[6] = (t + b) / (t - b)
    m[7] = 0

    m[8] = 0
    m[9] = 0
    m[10] = -(f + n) / (f - n)
    m[11] = -2 * f * n / (f - n)

    m[12] = 0
    m[13] = 0
    m[14] = -1
    m[15] = 0

    return result
};

/**
 * Returns a new M4 representing the a projection through/towards a point onto a plane.
 *
 * @param {V3} p
 * @param {P3} plane
 * @param {M4=} result
 * @returns {M4}
 */
M4.projectPlanePoint = function(p, plane, result) {
    assertVectors(p)
    assertInst(P3, plane)
    console.log('p', p.sce, plane)
    assert(!result || result instanceof M4, "!result || result instanceof M4")
    result = result || M4()
    var m = result.m
    var n = plane.normal, w = plane.w
    var np = n.dot(p)

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

    console.log(m)
    return result
};


/**
 * Orthographic/orthogonal projection. Transforms the cuboid with the dimensions X: [left; right] Y: [bottom, top]
 * Z: [near; far] to the cuboid X: [-1; 1] Y[-1; 1] Z [-1, 1]
 * @param {number} left
 * @param {number} right
 * @param {number} bottom
 * @param {number} top
 * @param {number} near
 * @param {number} far
 * @param {M4=} result
 * @returns {M4}
 */
M4.ortho = function(left, right, bottom, top, near, far, result) {
    assertNumbers(left, right, bottom, top, near, far)
    assert(!result || result instanceof M4, "!result || result instanceof M4")
    result = result || M4()
    var m = result.m

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

    return result;
};

/**
 * This emulates the OpenGL function `glScale()`. You can optionally pass an existing matrix in `result` to avoid
 * allocating a new matrix.
 *
 * @param {number|V3} x
 * @param {number|M4} [y]
 * @param {number|undefined} [z]
 * @param {M4|undefined} [result]
 * @returns {M4}
 */
M4.scaling = function(x, y, z, result) {
    if (1 == arguments.length || 2 == arguments.length) {
        assertVectors(x)
        result = y
        y = x.y
        z = x.z
        x = x.x
    } else {
        assert(3 == arguments.length || 4 == arguments.length, "3 == arguments.length || 4 == arguments.length")
        assertNumbers(x, y, z)
    }
    assert(!result || result instanceof M4, "!result || result instanceof M4")

    result = result || M4();
    var m = result.m;

    m[0] = x;
    m[1] = 0;
    m[2] = 0;
    m[3] = 0;

    m[4] = 0;
    m[5] = y;
    m[6] = 0;
    m[7] = 0;

    m[8] = 0;
    m[9] = 0;
    m[10] = z;
    m[11] = 0;

    m[12] = 0;
    m[13] = 0;
    m[14] = 0;
    m[15] = 1;

    return result;
};

/**
 * This emulates the OpenGL function `glTranslate()`. You can optionally pass
 * an existing matrix in `result` to avoid allocating a new matrix.
 *
 * @param {number|V3} x
 * @param {number|M4=} y
 * @param {number=} z
 * @param {M4=} result optional matrix to store result into
 * @returns {M4}
 */
M4.translation = function(x, y, z, result) {
    if (1 == arguments.length || 2 == arguments.length) {
        assertVectors(x)
        result = y
        y = x.y
        z = x.z
        x = x.x
    } else {
        assert(3 == arguments.length || 4 == arguments.length, "3 == arguments.length || 4 == arguments.length")
        assertNumbers(x, y, z)
    }
    assert(!result || result instanceof M4, "!result || result instanceof M4")

    result = result || M4();
    var m = result.m;

    m[0] = 1;
    m[1] = 0;
    m[2] = 0;
    m[3] = x;

    m[4] = 0;
    m[5] = 1;
    m[6] = 0;
    m[7] = y;

    m[8] = 0;
    m[9] = 0;
    m[10] = 1;
    m[11] = z;

    m[12] = 0;
    m[13] = 0;
    m[14] = 0;
    m[15] = 1;

    return result;
};

// ### GL.Matrix.rotate(a, x, y, z[, result])
//

/**
 * Returns a matrix that rotates by `a` degrees around the vector (x, y, z). You can optionally pass an existing
 * matrix in `result` to avoid allocating a new matrix. This emulates the OpenGL function `glRotate()`.
 * @param {number} radians
 * @param {V3|number} x
 * @param {number|M4=} y
 * @param {number=} z
 * @param {M4=} result
 * @returns {M4}
 */
M4.rotation = function(radians, x, y, z, result) {
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
    assert(!result || result instanceof M4, "!result || result instanceof M4")
    // TODO
    assert(!V3(x, y, z).isZero(), "!V3(x, y, z).isZero()")

    result = result || M4();
    var m = result.m;

    var d = Math.sqrt(x*x + y*y + z*z);
    x /= d; y /= d; z /= d;
    var cos = Math.cos(radians), sin = Math.sin(radians), t = 1 - cos;

    m[0] = x * x * t + cos;
    m[1] = x * y * t - z * sin;
    m[2] = x * z * t + y * sin;
    m[3] = 0;

    m[4] = y * x * t + z * sin;
    m[5] = y * y * t + cos;
    m[6] = y * z * t - x * sin;
    m[7] = 0;

    m[8] = z * x * t - y * sin;
    m[9] = z * y * t + x * sin;
    m[10] = z * z * t + cos;
    m[11] = 0;

    m[12] = 0;
    m[13] = 0;
    m[14] = 0;
    m[15] = 1;

    return result;
};

/**
 * Returns a matrix that puts the camera at the eye point `ex, ey, ez` looking
 * toward the center point `cx, cy, cz` with an up direction of `ux, uy, uz`.
 * You can optionally pass an existing matrix in `result` to avoid allocating
 * a new matrix. This emulates the OpenGL function `gluLookAt()`.
 *
 * @param {V3} e eye coordinates
 * @param {V3} c center coordinates
 * @param {V3} u up coordinates
 * @param {M4=} result optional matrix to store result into
 * @returns {M4}
 */
M4.lookAt = function(e, c, u, result) {
    assert(3 == arguments.length || 4 == arguments.length, "3 == arguments.length || 4 == arguments.length")
    NLA.assertVectors(e, c, u)
    assert(!result || result instanceof M4, "!result || result instanceof M4")

    result = result || M4();
    var m = result.m;

    var f = e.minus(c).unit()
    var s = u.cross(f).unit()
    var t = f.cross(s).unit()

    m[0] = s.x;
    m[1] = s.y;
    m[2] = s.z;
    m[3] = -s.dot(e);

    m[4] = t.x;
    m[5] = t.y;
    m[6] = t.z;
    m[7] = -t.dot(e);

    m[8] = f.x;
    m[9] = f.y;
    m[10] = f.z;
    m[11] = -f.dot(e);

    m[12] = 0;
    m[13] = 0;
    m[14] = 0;
    m[15] = 1;

    return result;
};

/**
 * Create a rotation matrix for rotating around the X axis
 * @param {number} radians
 * @returns {M4}
 */
M4.rotationX = function(radians) {
    assertNumbers(radians)
    var cos = Math.cos(radians);
    var sin = Math.sin(radians);
    var els = [
        1, 0, 0, 0, 0, cos, -sin, 0, 0, sin, cos, 0, 0, 0, 0, 1
    ];
    return M4(els);
};

/**
 * Create a rotation matrix for rotating around the Y axis
 * @param {number} radians
 * @returns {M4}
 */
M4.rotationY = function(radians) {
    var cos = Math.cos(radians);
    var sin = Math.sin(radians);
    var els = [
        cos, 0, sin, 0, 0, 1, 0, 0, -sin, 0, cos, 0, 0, 0, 0, 1
    ];
    return M4(els);
};

/**
 * Create a rotation matrix for rotating around the Z axis
 * @param {number} radians
 * @returns {M4}
 */
M4.rotationZ = function(radians) {
    var cos = Math.cos(radians);
    var sin = Math.sin(radians);
    var els = [
        cos, -sin, 0, 0, sin, cos, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1
    ];
    return M4(els);
};

/**
 * Matrix for rotation about arbitrary line defined by a point and axis
 * @param {V3} rotationCenter
 * @param {V3} rotationAxis
 * @param {number} radians
 * @param {M4=} result
 * @returns {*|M4}
 */
M4.rotationLine = function(rotationCenter, rotationAxis, radians, result) {
    // see http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
    assertVectors(rotationCenter)
    assertVectors(rotationAxis)
    assert(!result || result instanceof M4, "!result || result instanceof M4")
    result = result || M4()
    rotationAxis = rotationAxis.normalized()

    var a = rotationCenter.x, b = rotationCenter.y, c = rotationCenter.z,
        u = rotationAxis.x, v = rotationAxis.y, w = rotationAxis.z
    var m = result.m, cos = Math.cos(radians), sin = Math.sin(radians)

    m[0] = u * u + (v * v + w * w) * cos
    m[1] = u * v * (1 - cos) - w * sin
    m[2] = u * w * (1 - cos) + v * sin
    m[3] = (a * (v * v + w * w) - u * (b * v + c * w)) * (1 - cos) + (b * w - c * v) * sin

    m[4] = u * v * (1 - cos) + w * sin
    m[5] = v * v + (u * u + w * w) * cos
    m[6] = v * w * (1 - cos) - u * sin
    m[7] = (b * (u * u + w * w) - v * (a * u + c * w)) * (1 - cos) + (c * u - a * w) * sin

    m[8] = u * w * (1 - cos) - v * sin
    m[9] = v * w * (1 - cos) + u * sin
    m[10] = w * w + (u * u + v * v) * cos
    m[11] = (c * (u * u + v * v) - w * (a * u + b * v)) * (1 - cos) + (a * v - b * u) * sin

    m[12] = 0
    m[13] = 0
    m[14] = 0
    m[15] = 1
    return result
    /*
    var rotationPlane = NLA.Plane.fromNormalAndPoint(rotationAxis, rotationCenter);
    var orthobasis = new NLA.OrthoNormalBasis(rotationPlane);
    var transformation = M4.translation(rotationCenter.negated());
    transformation = transformation.multiply(orthobasis.getProjectionMatrix());
    transformation = transformation.multiply(M4.rotationZ(radians));
    transformation = transformation.multiply(orthobasis.getInverseProjectionMatrix());
    transformation = transformation.multiply(M4.translation(rotationCenter));
    return transformation;
    */
};

/**
 * Create an affine matrix for mirroring into an arbitrary plane:
 *
 * @param {P3} plane
 * @param {M4=} result
 * @returns {M4}
 */
M4.mirroring = function(plane, result) {
    assertInst(P3, plane)
    assert(!result || result instanceof M4, "!result || result instanceof M4")
    var nx = plane.normal.x;
    var ny = plane.normal.y;
    var nz = plane.normal.z;
    var w = plane.w;
    result = result || M4()
    var m = result.m

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
};

/**
 *
 * @param {P3} plane
 * @param {V3} [dir] Projection direction. Optional, if not specified plane normal will be used.
 * @param {M4=} result
 * @returns {M4}
 */
M4.projection = function (plane, dir, result) {
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
	assert(!dir || dir instanceof V3, "!dir || dir instanceof V3")
	assert(!result || result instanceof M4, "!result || result instanceof M4")
	dir = dir || plane.normal
	var {x: nx, y: ny, z: nz} = plane.normal
	var {x: dx, y: dy, z: dz} = dir
	var w = plane.w;
	result = result || M4()
	var m = result.m
	var nd = plane.normal.dot(dir)
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

	dx /= nd
	dy /= nd
	dz /= nd

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

/**
 *
 * @param line
 * @param result
 */
M4.lineProjection = function (line, result) {
	assertInst(L3, line)
	assert(!result || result instanceof M4, "!result || result instanceof M4")
	var dx = line.dir1.x, dy = line.dir1.y, dz = line.dir1.z
	var ax = line.anchor.x, ay = line.anchor.y, az = line.anchor.z
	result = result || M4()
	var m = result.m

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
}
M4.multiplyMultiple = function () {
	var m4s = arguments
	if (0 == m4s.length) return M4.identity()
	var temp = M4(), result = m4s[0].copy()
	for (var i = 1; i < m4s.length; i++) {
        M4.multiply(result, m4s[i], temp); // <-- ; required
		[temp, result] = [result, temp]
	}
	return result
}

/**
 *
 * @param {V3} p
 * @param {M4=} result
 * @returns {M4}
 */
M4.pointInversion = function(p, result) {
	assertVectors(p)
	assert(!result || result instanceof M4, "!result || result instanceof M4")
	result = result || M4()
	var m = result.m

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

/**
 * Temp variables used
 * @type {M4}
 */
M4.temp0 = M4()
M4.temp1 = M4()
M4.temp2 = M4()

/**
 * A simple (consists of integers), regular, non-orthogonal matrix, useful mainly for testing.
 * M4.BAR = M4.FOO.inverse()
 *
 * @const
 * @type {M4}
 */
M4.FOO = M4(
	0, 1, 1, 2,
	0.3, 0.4, 0.8, 13,
	2.1, 3.4, 5.5, 8.9,
	0, 0, 0, 1)
M4.FOO.name = 'M4.FOO'
/**
 * Inverse of M4.FOO
 * @const
 * @type {M4}
 */
M4.BAR = M4.FOO.inversed()
M4.BAR.name = 'M4.BAR'

/**
 * @const
 * @type {M4}
 */
M4.IDENTITY = M4.identity()
M4.IDENTITY.name = 'M4.IDENTITY'

M4.IDENTITY3 = M4(
	1, 0, 0, 0,
	0, 1, 0, 0,
	0, 0, 1, 0,
	0, 0, 0, 0
)
M4.IDENTITY3.name = 'M4.IDENTITY3'