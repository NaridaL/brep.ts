/**
 * Created by aval on 20/11/2015.
 */
"use strict";
if (!NLA.Vector3) {
    throw new Error("Define NLA.V3 first")
}
(function (module) {
    var V3 = NLA.Vector3

    var assert = NLA.assert, assertNumbers = NLA.assertNumbers
    // This constructor takes 16 arguments in row-major order, which can be passed
    // individually, as a list, or even as four lists, one for each row. If the
    // arguments are omitted then the identity matrix is constructed instead.
    var M4 = function() {
	    var m
        if (0 == arguments.length) {
	        m = new Float64Array(16)
        } else {
	        var flattened = Array.prototype.concat.apply([], arguments)
	        assert(flattened.length == 16, "flattened.length == 16")
	        m = new Float64Array(flattened);
        }
	    var o = Object.create(M4.prototype)
	    Object.defineProperty(o, "m", {value: m})
	    return o
    }
	M4.prototype = Object.create(NLA.Matrix.prototype)

    var M4prot = {
        // ### .inverse()
        //
        // Returns the matrix that when multiplied with this matrix results in the
        // identity matrix.
        inversed: function() {
            return M4.inverse(this);
        },
	    trace: function () {
		    return this.m[0] + this.m[5] + this.m[10] + this.m[15]
	    },
	    likeMatrix4: function (m4) {
		    assert(m4 instanceof M4, "m4 instanceof M4")
		    return this.m.every((el, index) => NLA.equals(el, m4.m[index]))
	    },

        // ### .transpose()
        //
        // Returns this matrix, exchanging columns for rows.
        transposed: function() {
            return M4.transpose(this);
        },

        // ### .multiply(matrix)
        //
        // Returns the concatenation of the transforms for this matrix and `matrix`.
        // This emulates the OpenGL function `glMultMatrix()`.
        times: function(matrix) {
            return M4.multiply(this, matrix);
        },

        // ### .transformPoint(point)
        //
        // Transforms the vector as a point with a w coordinate of 1. This
        // means translations will have an effect, for example.
	    transformPoint: function(v) {
		    assert(v instanceof V3, "v instanceof V3")
		    var m = this.m;
		    var v0 = v.x;
		    var v1 = v.y;
		    var v2 = v.z;
		    var v3 = 1;
		    var x = v0 * m[ 0] + v1 * m[ 1] + v2 * m[ 2] + v3 * m[ 3];
		    var y = v0 * m[ 4] + v1 * m[ 5] + v2 * m[ 6] + v3 * m[ 7];
		    var z = v0 * m[ 8] + v1 * m[ 9] + v2 * m[10] + v3 * m[11];
		    var w = v0 * m[12] + v1 * m[13] + v2 * m[14] + v3 * m[15];
		    // scale such that fourth element becomes 1:
		    return V3.create(x / w, y / w, z / w);
	    },

        // ### .transformPoint(vector)
        //
        // Transforms the vector as a vector with a w coordinate of 0. This
        // means translations will have no effect, for example.
        transformVector: function(v) {
	        assert(v instanceof V3, "v instanceof V3")
            var m = this.m;
            return V3.create(
                m[0] * v.x + m[1] * v.y + m[2] * v.z,
                m[4] * v.x + m[5] * v.y + m[6] * v.z,
                m[8] * v.x + m[9] * v.y + m[10] * v.z
            );
        },
	    transformedPoints: function(vs) {
		    return vs.map(v => this.transformPoint(v))
	    },
	    transformedVectors: function(vs) {
		    return vs.map(v => this.transformVector(v))
	    },

        plus: function(m) {
            var r = M4()
            for (var i = 0; i < 16; i++) {
                r.m[i] = this.m[i] + m.m[i];
            }
            return r;
        },

        minus: function(m) {
            var r = M4()
            for (var i = 0; i < 16; i++) {
                r.m[i] = this.m[i] - m.m[i];
            }
            return r;
        },

	    clone: function() {
		    return M4.copy(this)
	    },
	    copy: function() {
		    return M4.copy(this)
	    },
	    isRegular: function () {
		    return !NLA.isZero(this.determinant())
	    },
		isAxisAligned: function () {
			var $ = this.m
			return (1 >= !NLA.isZero($[0]) + !NLA.isZero($[1]) + !NLA.isZero($[2]))
				&& (1 >= !NLA.isZero($[4]) + !NLA.isZero($[5]) + !NLA.isZero($[6]))
				&& (1 >= !NLA.isZero($[8]) + !NLA.isZero($[9]) + !NLA.isZero($[10]))
		},

        // Multiply a NLA.V3D (interpreted as 3 column, 1 row) by this matrix
        // (result = v*M)
        // Fourth element is taken as 1
        leftMultiplyVector3: function(v) {
            var v0 = v.x;
            var v1 = v.y;
            var v2 = v.z;
            var v3 = 1;
            var x = v0 * this.m[0] + v1 * this.m[4] + v2 * this.m[8] + v3 * this.m[12];
            var y = v0 * this.m[1] + v1 * this.m[5] + v2 * this.m[9] + v3 * this.m[13];
            var z = v0 * this.m[2] + v1 * this.m[6] + v2 * this.m[10] + v3 * this.m[14];
            var w = v0 * this.m[3] + v1 * this.m[7] + v2 * this.m[11] + v3 * this.m[15];
            // scale such that fourth element becomes 1:
            return V3.create(x / w, y / w, z / w);
        },

        // Right multiply the matrix by a NLA.Vector2D (interpreted as 2 row, 1 column)
        // (result = M*v)
        // Fourth element is taken as 1
        rightMultiply1x2Vector: function(v) {
            var v0 = v.x;
            var v1 = v.y;
            var v2 = 0;
            var v3 = 1;
            var x = v0 * this.m[ 0] + v1 * this.m[ 1] + v2 * this.m[ 2] + v3 * this.m[ 3];
            var y = v0 * this.m[ 4] + v1 * this.m[ 5] + v2 * this.m[ 6] + v3 * this.m[ 7];
            var w = v0 * this.m[12] + v1 * this.m[13] + v2 * this.m[14] + v3 * this.m[15];
            // scale such that w becomes 1:
            return new NLA.Vector2D(x / w, y / w);
        },

        // Multiply a NLA.Vector2D (interpreted as 2 column, 1 row) by this matrix
        // (result = v*M)
        // Fourth element is taken as 1
        leftMultiply1x2Vector: function(v) {
            var v0 = v.x;
            var v1 = v.y;
            var v2 = 0;
            var v3 = 1;
            var x = v0 * this.m[0] + v1 * this.m[4] + v2 * this.m[8] + v3 * this.m[12];
            var y = v0 * this.m[1] + v1 * this.m[5] + v2 * this.m[9] + v3 * this.m[13];
            var z = v0 * this.m[2] + v1 * this.m[6] + v2 * this.m[10] + v3 * this.m[14];
            var w = v0 * this.m[3] + v1 * this.m[7] + v2 * this.m[11] + v3 * this.m[15];
            // scale such that fourth element becomes 1:
            return new NLA.Vector2D(x / w, y / w);
        },
	    isOrthogonal: function () {
		    var t = M4.temp0, s = M4.temp1
		    M4.transpose(this, t)
		    M4.multiply(this, t, s)
		    M4.identity(t)
		    return t.likeMatrix4(s)
	    },
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
        // determine whether this matrix is a mirroring transformation
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

	        return this.determinant() < 0
        },

	    /**
	     * Returns if the matrix has the following form:
	     * a b c 0
	     * c d e 0
	     * f g h 0
	     * 0 0 0 1
	     * @returns {*}
	     */
	    is3x3: function () {
		    var m = this.m
			return NLA.equals(1, m[15])
		        && NLA.isZero(m[12]) && NLA.isZero(m[13]) && NLA.isZero(m[14])
		        && NLA.isZero(m[3]) && NLA.isZero(m[7]) && NLA.isZero(m[11])
	    },
	    toString: function (f) {
	        f = f || ((v) => v.toFixed(6).replace(/(0|\.)(?=0*$)/g, " ").toString())
		    assert(typeof f(0) == "string", "" + typeof f(0))
		    // slice this.m to convert it to an Array (from TypeArray)
		    var rounded = Array.prototype.slice.call(this.m).map(f);
		    var colWidths = [0, 1, 2, 3].map((colIndex) => rounded.sliceStep(colIndex, 4).map((x) => x.length).max());
		    return [0, 1, 2, 3].map(
			    (rowIndex) => rounded
				    .slice(rowIndex * 4, rowIndex * 4 + 4) // select matrix row
				    .map((x, colIndex) => NLA.repeatChar(colWidths[colIndex] - x.length, ' ') + x) // pad numbers with spaces to col width
				    .join(" ")
		    ).join("\n"); // join rows
	    }
    };
	for (var prop in M4prot) {
		M4.prototype[prop] = M4prot[prop]
	}
    // ### GL.Matrix.inverse(matrix[, result])
    //
    // Returns the matrix that when multiplied with `matrix` results in the
    // identity matrix. You can optionally pass an existing matrix in `result`
    // to avoid allocating a new matrix. This implementation is from the Mesa
    // OpenGL function `__gluInvertMatrixd()` found in `project.c`.
    M4.inverse = function(matrix, result) {
	    assert(matrix instanceof M4, "matrix instanceof M4")
	    assert(!result || result instanceof M4, "!result || result instanceof M4")
	    assert(matrix != result, "matrix != result")
        result = result || M4();
        var m = matrix.m, r = result.m;

        r[0] = m[5]*m[10]*m[15] - m[5]*m[14]*m[11] - m[6]*m[9]*m[15] + m[6]*m[13]*m[11] + m[7]*m[9]*m[14] - m[7]*m[13]*m[10];
        r[1] = -m[1]*m[10]*m[15] + m[1]*m[14]*m[11] + m[2]*m[9]*m[15] - m[2]*m[13]*m[11] - m[3]*m[9]*m[14] + m[3]*m[13]*m[10];
        r[2] = m[1]*m[6]*m[15] - m[1]*m[14]*m[7] - m[2]*m[5]*m[15] + m[2]*m[13]*m[7] + m[3]*m[5]*m[14] - m[3]*m[13]*m[6];
        r[3] = -m[1]*m[6]*m[11] + m[1]*m[10]*m[7] + m[2]*m[5]*m[11] - m[2]*m[9]*m[7] - m[3]*m[5]*m[10] + m[3]*m[9]*m[6];

        r[4] = -m[4]*m[10]*m[15] + m[4]*m[14]*m[11] + m[6]*m[8]*m[15] - m[6]*m[12]*m[11] - m[7]*m[8]*m[14] + m[7]*m[12]*m[10];
        r[5] = m[0]*m[10]*m[15] - m[0]*m[14]*m[11] - m[2]*m[8]*m[15] + m[2]*m[12]*m[11] + m[3]*m[8]*m[14] - m[3]*m[12]*m[10];
        r[6] = -m[0]*m[6]*m[15] + m[0]*m[14]*m[7] + m[2]*m[4]*m[15] - m[2]*m[12]*m[7] - m[3]*m[4]*m[14] + m[3]*m[12]*m[6];
        r[7] = m[0]*m[6]*m[11] - m[0]*m[10]*m[7] - m[2]*m[4]*m[11] + m[2]*m[8]*m[7] + m[3]*m[4]*m[10] - m[3]*m[8]*m[6];

        r[8] = m[4]*m[9]*m[15] - m[4]*m[13]*m[11] - m[5]*m[8]*m[15] + m[5]*m[12]*m[11] + m[7]*m[8]*m[13] - m[7]*m[12]*m[9];
        r[9] = -m[0]*m[9]*m[15] + m[0]*m[13]*m[11] + m[1]*m[8]*m[15] - m[1]*m[12]*m[11] - m[3]*m[8]*m[13] + m[3]*m[12]*m[9];
        r[10] = m[0]*m[5]*m[15] - m[0]*m[13]*m[7] - m[1]*m[4]*m[15] + m[1]*m[12]*m[7] + m[3]*m[4]*m[13] - m[3]*m[12]*m[5];
        r[11] = -m[0]*m[5]*m[11] + m[0]*m[9]*m[7] + m[1]*m[4]*m[11] - m[1]*m[8]*m[7] - m[3]*m[4]*m[9] + m[3]*m[8]*m[5];

        r[12] = -m[4]*m[9]*m[14] + m[4]*m[13]*m[10] + m[5]*m[8]*m[14] - m[5]*m[12]*m[10] - m[6]*m[8]*m[13] + m[6]*m[12]*m[9];
        r[13] = m[0]*m[9]*m[14] - m[0]*m[13]*m[10] - m[1]*m[8]*m[14] + m[1]*m[12]*m[10] + m[2]*m[8]*m[13] - m[2]*m[12]*m[9];
        r[14] = -m[0]*m[5]*m[14] + m[0]*m[13]*m[6] + m[1]*m[4]*m[14] - m[1]*m[12]*m[6] - m[2]*m[4]*m[13] + m[2]*m[12]*m[5];
        r[15] = m[0]*m[5]*m[10] - m[0]*m[9]*m[6] - m[1]*m[4]*m[10] + m[1]*m[8]*m[6] + m[2]*m[4]*m[9] - m[2]*m[8]*m[5];

        var det = m[0]*r[0] + m[1]*r[4] + m[2]*r[8] + m[3]*r[12];
	   // assert(!NLA.isZero(det), "det may not be zero, i.e. the matrix is not invertible")
        var i = 16
        while (i--) { r[i] /= det }
        return result;
    };

    // ### GL.Matrix.transpose(matrix[, result])
    //
    // Returns `matrix`, exchanging columns for rows. You can optionally pass an
    // existing matrix in `result` to avoid allocating a new matrix.
    M4.transpose = function(matrix, result) {
	    assert(matrix instanceof M4, "matrix instanceof M4")
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

    // ### GL.Matrix.multiply(left, right[, result])
    //
    // Returns the concatenation of the transforms for `left` and `right`. You can
    // optionally pass an existing matrix in `result` to avoid allocating a new
    // matrix. This emulates the OpenGL function `glMultMatrix()`.
    M4.multiply = function(left, right, result) {
	    assert(left instanceof M4, "left instanceof M4")
	    assert(right instanceof M4, "right instanceof M4")
	    assert(!result || result instanceof M4, "!result || result instanceof M4")
	    assert(left != result, "left != result")
	    assert(right != result, "right != result")
        result = result || M4();
        var a = left.m, b = right.m, r = result.m;

        r[0] = a[0] * b[0] + a[1] * b[4] + a[2] * b[8] + a[3] * b[12];
        r[1] = a[0] * b[1] + a[1] * b[5] + a[2] * b[9] + a[3] * b[13];
        r[2] = a[0] * b[2] + a[1] * b[6] + a[2] * b[10] + a[3] * b[14];
        r[3] = a[0] * b[3] + a[1] * b[7] + a[2] * b[11] + a[3] * b[15];

        r[4] = a[4] * b[0] + a[5] * b[4] + a[6] * b[8] + a[7] * b[12];
        r[5] = a[4] * b[1] + a[5] * b[5] + a[6] * b[9] + a[7] * b[13];
        r[6] = a[4] * b[2] + a[5] * b[6] + a[6] * b[10] + a[7] * b[14];
        r[7] = a[4] * b[3] + a[5] * b[7] + a[6] * b[11] + a[7] * b[15];

        r[8] = a[8] * b[0] + a[9] * b[4] + a[10] * b[8] + a[11] * b[12];
        r[9] = a[8] * b[1] + a[9] * b[5] + a[10] * b[9] + a[11] * b[13];
        r[10] = a[8] * b[2] + a[9] * b[6] + a[10] * b[10] + a[11] * b[14];
        r[11] = a[8] * b[3] + a[9] * b[7] + a[10] * b[11] + a[11] * b[15];

        r[12] = a[12] * b[0] + a[13] * b[4] + a[14] * b[8] + a[15] * b[12];
        r[13] = a[12] * b[1] + a[13] * b[5] + a[14] * b[9] + a[15] * b[13];
        r[14] = a[12] * b[2] + a[13] * b[6] + a[14] * b[10] + a[15] * b[14];
        r[15] = a[12] * b[3] + a[13] * b[7] + a[14] * b[11] + a[15] * b[15];

        return result;
    };
	M4.copy = function(src, result) {
		assert(src instanceof M4, "src instanceof M4")
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
    M4.forSys = function (e0, e1, e2, origin) {
	    assert(e0 instanceof V3, "e0 instanceof V3")
	    assert(e1 instanceof V3, "e1 instanceof V3")
	    assert(!e2 || e2 instanceof V3, "!e2 || e2 instanceof V3")

        e2 = e2 || e0.cross(e1);
	    origin = origin || V3.ZERO
        return M4(
            e0.x, e1.x, e2.x, origin.x,
            e0.y, e1.y, e2.y, origin.y,
            e0.z, e1.z, e2.z, origin.z,
            0,	     0,    0, 1);
    }

    // ### GL.Matrix.identity([result])
    //
    // Returns an identity matrix. You can optionally pass an existing matrix in
    // `result` to avoid allocating a new matrix. This emulates the OpenGL function
    // `glLoadIdentity()`.
    M4.identity = function(result) {
	    assert(!result || result instanceof M4, "!result || result instanceof M4")
        result = result || M4();
        var m = result.m;
        m[0] = m[5] = m[10] = m[15] = 1;
        m[1] = m[2] = m[3] = m[4] = m[6] = m[7] = m[8] = m[9] = m[11] = m[12] = m[13] = m[14] = 0;

        return result;
    };

	M4.fromFunction = function (f, result) {
		assert(typeof f == "function", 'typeof f == "function"' + typeof f)
		assert(!result || result instanceof M4, "!result || result instanceof M4")
		result = result || M4()
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
    M4.perspective = function(fov, aspect, near, far, result) {
	    assert(!result || result instanceof M4, "!result || result instanceof M4")
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
    M4.frustum = function(l, r, b, t, n, f, result) {
	    assertNumbers(l, r, b, t, n, f)
	    assert(!result || result instanceof M4, "!result || result instanceof M4")
        result = result || M4();
        var m = result.m;

        m[0] = 2 * n / (r - l);
        m[1] = 0;
        m[2] = (r + l) / (r - l);
        m[3] = 0;

        m[4] = 0;
        m[5] = 2 * n / (t - b);
        m[6] = (t + b) / (t - b);
        m[7] = 0;

        m[8] = 0;
        m[9] = 0;
        m[10] = -(f + n) / (f - n);
        m[11] = -2 * f * n / (f - n);

        m[12] = 0;
        m[13] = 0;
        m[14] = -1;
        m[15] = 0;

        return result;
    };

	/**
	 * Orthographic/orthogonal projection. Transforms the cuboid with the dimensions X: [left; right] Y: [bottom, top] Z: [near; far]
	 * to the cuboid X: [-1; 1] Y[-1; 1] Z [-1, 1]
	 * @param left
	 * @param right
	 * @param bottom
	 * @param top
	 * @param near
	 * @param far
	 * @param result
	 * @returns {*}
	 */
    // ### GL.Matrix.ortho(left, right, bottom, top, near, far[, result])
    //
    // Returns an orthographic projection, in which objects are the same size no
    // matter how far away or nearby they are. You can optionally pass an existing
    // matrix in `result` to avoid allocating a new matrix. This emulates the OpenGL
    // function `glOrtho()`.
    M4.ortho = function(left, right, bottom, top, near, far, result) {
	    assertNumbers(left, right, bottom, top, near, far)
	    assert(!result || result instanceof M4, "!result || result instanceof M4")
        result = result || M4();
        var m = result.m;

        m[0] = 2 / (right - left);
        m[1] = 0;
        m[2] = 0;
        m[3] = -(right + left) / (right - left);

        m[4] = 0;
        m[5] = 2 / (top - bottom);
        m[6] = 0;
        m[7] = -(top + bottom) / (top - bottom);

        m[8] = 0;
        m[9] = 0;
        m[10] = -2 / (far - near);
        m[11] = -(far + near) / (far - near);

        m[12] = 0;
        m[13] = 0;
        m[14] = 0;
        m[15] = 1;

        return result;
    };

    // ### GL.Matrix.scale(x, y, z[, result])
    //
    // This emulates the OpenGL function `glScale()`. You can optionally pass an
    // existing matrix in `result` to avoid allocating a new matrix.
    M4.scaling = function(x, y, z, result) {
        if (1 == arguments.length || 2 == arguments.length) {
            assert(x instanceof V3, "x instanceof V3")
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

    // ### GL.Matrix.translation(x, y, z[, result])
    //
    // This emulates the OpenGL function `glTranslate()`. You can optionally pass
    // an existing matrix in `result` to avoid allocating a new matrix.
    M4.translation = function(x, y, z, result) {
        if (1 == arguments.length || 2 == arguments.length) {
            assert(x instanceof V3, "" + x)
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
    // Returns a matrix that rotates by `a` degrees around the vector `x, y, z`.
    // You can optionally pass an existing matrix in `result` to avoid allocating
    // a new matrix. This emulates the OpenGL function `glRotate()`.
    M4.rotation = function(radians, x, y, z, result) {
        if (2 == arguments.length || 3 == arguments.length) {
            assert(x instanceof V3, "x instanceof V3")
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
        var c = Math.cos(radians), s = Math.sin(radians), t = 1 - c;

        m[0] = x * x * t + c;
        m[1] = x * y * t - z * s;
        m[2] = x * z * t + y * s;
        m[3] = 0;

        m[4] = y * x * t + z * s;
        m[5] = y * y * t + c;
        m[6] = y * z * t - x * s;
        m[7] = 0;

        m[8] = z * x * t - y * s;
        m[9] = z * y * t + x * s;
        m[10] = z * z * t + c;
        m[11] = 0;

        m[12] = 0;
        m[13] = 0;
        m[14] = 0;
        m[15] = 1;

        return result;
    };

    // ### GL.Matrix.lookAt(ex, ey, ez, cx, cy, cz, ux, uy, uz[, result])
    //
    // Returns a matrix that puts the camera at the eye point `ex, ey, ez` looking
    // toward the center point `cx, cy, cz` with an up direction of `ux, uy, uz`.
    // You can optionally pass an existing matrix in `result` to avoid allocating
    // a new matrix. This emulates the OpenGL function `gluLookAt()`.
    M4.lookAt = function(ex, ey, ez, cx, cy, cz, ux, uy, uz, result) {
        if (3 == arguments.length || 4 == arguments.length) {
            var eye = ex, center = ey, up = ez, result = cx
            assert(eye instanceof V3, "eye instanceof V3")
            assert(center instanceof V3, "center instanceof V3")
            assert(up instanceof V3, "up instanceof V3")
            return M4.lookAt(eye.x, eye.y, eye.z, center.x, center.y, center.z, up.x, up.y, up.z, result)
        } else {
            assert(9 == arguments.length || 10 == arguments.length, "9 == arguments.length || 10 == arguments.length")
            assertNumbers(ex, ey, ez, cx, cy, cz, ux, uy, uz)
        }
        assert(!result || result instanceof M4, "!result || result instanceof M4")

        result = result || M4();
        var m = result.m;

        var e = V3.create(ex, ey, ez);
        var c = V3.create(cx, cy, cz);
        var u = V3.create(ux, uy, uz);
        var f = e.minus(c).unit();
        var s = u.cross(f).unit();
        var t = f.cross(s).unit();

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

    // Create a rotation matrix for rotating around the x axis
    M4.rotationX = function(radians) {
        assertNumbers(radians)
        var cos = Math.cos(radians);
        var sin = Math.sin(radians);
        var els = [
            1, 0, 0, 0, 0, cos, -sin, 0, 0, sin, cos, 0, 0, 0, 0, 1
        ];
        return M4(els);
    };

    // Create a rotation matrix for rotating around the y axis
    M4.rotationY = function(radians) {
        var cos = Math.cos(radians);
        var sin = Math.sin(radians);
        var els = [
            cos, 0, sin, 0, 0, 1, 0, 0, -sin, 0, cos, 0, 0, 0, 0, 1
        ];
        return M4(els);
    };

    // Create a rotation matrix for rotating around the z axis
    M4.rotationZ = function(radians) {
        var cos = Math.cos(radians);
        var sin = Math.sin(radians);
        var els = [
            cos, -sin, 0, 0, sin, cos, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1
        ];
        return M4(els);
    };

    // Matrix for rotation about arbitrary point and axis
    M4.rotationLine = function(rotationCenter, rotationAxis, radians, result) {
	    // see http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
	    assert(rotationCenter instanceof V3, "rotationCenter instanceof V3")
	    assert(rotationAxis instanceof V3, "rotationAxis instanceof V3")
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


    // Create an affine matrix for mirroring into an arbitrary plane:
    M4.mirroring = function(plane, result) {
	    assert(plane instanceof P3, "plane instanceof P3")
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
	M4.pointInversion = function(p, result) {
		assert(p instanceof V3, "p instanceof V3")
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
	M4.temp0 = M4()
	M4.temp1 = M4()
    module.Matrix4x4 = M4

})(NLA)