/**
 * Created by aval on 20/11/2015.
 */

/**
 * Oriented plane, i.e. splits R^3 in half, with one half being "in front" of the plane.
 * Leads to multiple comparisons: isCoplanarToPlane returns if the plane occupies the same space,
 * like returns if the plane occupies the same space and has the same orientation
 *
 * Points x on the plane fulfill the equation: normal DOT x = w
 *
 * @param {V3} normal1 normalized plane normal
 * @param {number} w signed (rel to normal1) distance from the origin
 * @param {Object=} prototype object to set as prototype of the new Plane3 object. Defaults to Plane3.prototype
 * @alias NLA.Plane3
 * @constructor NLA.Plane3
 * @augments {NLA.Transformable}
 * @property {V3} normal
 * @property {number} w
 */
var P3 = function (normal1, w, prototype) {
    NLA.assertVectors(normal1)
    NLA.assertNumbers(w)
    assert(normal1.hasLength(1), "normal1.hasLength(1)")
    var p = Object.create(prototype || P3.prototype)
    p.w = w
	p.normal = normal1
	return p
}
/**
 * @alias NLA.Plane3.throughPoints
 * @param {V3} a
 * @param {V3} b
 * @param {V3} c
 * @param prototype
 * @returns {NLA.Plane3}
 */
P3.throughPoints = function (a, b, c, prototype) {
    NLA.assertVectors(a, b, c)
    var n1 = b.minus(a).cross(c.minus(a)).normalized();
    return P3(n1, n1.dot(a), prototype)
}

/**
 * @alias NLA.Plane3.normalOnAnchor
 * @param normal
 * @param anchor
 * @param prototype
 * @returns {NLA.Plane3}
 */
P3.normalOnAnchor = function (normal, anchor, prototype) {
    NLA.assertVectors(normal, anchor)
    var n1 = normal.normalized()
    return P3(n1, n1.dot(anchor), prototype)
}

P3.forAxisIntercepts = function (x0, y0, z0) {
	NLA.assertNumbers(x0, y0, z0)
	return P3(V3.create(1/x0, 1/y0, 1/z0))
}

/**
 * @alias NLA.Plane3.forAnchorAndPlaneVectors
 * @param anchor
 * @param v0
 * @param v1
 * @param prototype
 * @returns {NLA.Plane3}
 */
P3.forAnchorAndPlaneVectors = function (anchor, v0, v1, prototype) {
    NLA.assertVectors(anchor, v0, v1)
    return P3.normalOnAnchor(v0.cross(v1), anchor, prototype)
}


P3.prototype = /** @lends {NLA.Plane3.prototype} */ {
	/**
	 *
	 * @returns {V3}
	 */
	axisIntercepts: function () {
		let w = this.w, n = this.normal
		return V3.create(w / n.x, w / n.y, w / n.z)
	},
    get anchor () {
        return this.normal.times(this.w)
    },
    // Returns true iff the plane occupies the same space as the argument
    isCoplanarToPlane: function (plane) {
        assertInst(P3, plane)
        return this.like(plane) || this.likeFlipped(plane)
    },
    like: function (plane) {
	    assertInst(P3, plane)
	    return NLA.equals(this.w, plane.w) && this.normal.like(plane.normal)
    },
    likeFlipped: function (plane) {
	    assertInst(P3, plane)
	    return NLA.equals(this.w, -plane.w) && this.normal.like(plane.normal.negated())
    },

	/**
	 * True iff plane.normal is equal to this.normal or it's negation.
	 *
	 * @param {NLA.Plane3} plane
	 */
    isParallelToPlane: function (plane) {
        assertInst(P3, plane)
        return NLA.equals(1, Math.abs(this.normal.dot(plane.normal)))
    },

    isParallelToLine: function (line) {
        assertInst(L3, line)
        return NLA.isZero(this.normal.dot(line.dir1))
    },

    isPerpendicularToLine: function (line) {
        assertInst(L3, line)
        // this.normal || line.dir1
        return NLA.equals(1, Math.abs(this.normal.dot(line.normal)))
    },

    isPerpendicularToPlane: function (plane) {
        assertInst(P3, plane)
        return NLA.isZero(this.normal.dot(plane.normal))
    },
    toString: function (roundFunction) {
	    roundFunction = roundFunction || (v => v) //((v) => +v.toFixed(3))
	    return "P3("+this.normal.toString(roundFunction) + ", " + roundFunction(this.w) +")"
    },
    translated: function (offset) {
	    return P3(this.normal, this.w + offset.dot(this.normal))
    },

    transform: function(m4) {
	    var mirror = m4.isMirroring()
	    // get two vectors in the plane:
	    var u = this.normal.getPerpendicular()
	    var v = u.cross(this.normal)
	    // get 3 points in the plane:
	    var p1 = m4.transformPoint(this.anchor),
	        p2 = m4.transformPoint(this.anchor.plus(v)),
		    p3 = m4.transformPoint(this.anchor.plus(u))
	    // and create a new plane from the transformed points:
	    return P3.throughPoints(p1, !mirror ? p2 : p3, !mirror ? p3 : p2)
    },
    distanceToLine: function (line) {
	    assertInst(L3, line)
	    if (!this.isParallelToLine(line)) {
		    return this.distanceToPoint(line.anchor)
	    } else {
		    return 0
	    }
    },
    containsPoint: function (x) {
        NLA.assertVectors(x)
        return NLA.equals(this.w, this.normal.dot(x))
    },
    containsLine: function (line) {
        assertInst(L3, line)
        return this.containsPoint(line.anchor) && this.isParallelToLine(line)
    },
    distanceToPointSigned: function (point) {
	    NLA.assertInst(V3, point)
	    return this.normal.dot(point) - this.w
    },
    distanceToPoint: function (point) {
	    NLA.assertInst(V3, point)
	    return Math.abs(this.normal.dot(point) - this.w)
    },

    intersectionWithLine: function (line) {
		line.intersectionWithPlane(this)
    },

	/**
	 *
	 * @param plane
	 * @returns {L3}
	 */
    intersectionWithPlane: function (plane) {
        /*

            this: n0 * x = w0
            plane: n1 * x = w1
            plane perpendicular to both which goes through origin:
                n2 := n0 X x1
                n2 * x = 0
        */
        assertInst(P3, plane)
        assert(!this.isParallelToPlane(plane), "!this.isParallelToPlane(plane)")
        /*
        var n0 = this.normal, n1 = plane.normal, n2 = n0.cross(n1).normalized(), m = M4.forSys(n0, n1, n2)
        var x0 = this.anchor, x1 = plane.anchor, x2 = V3.ZERO
        var p = n2.times(x2.dot(n2))
	        .plus(n1.cross(n2).times(x0.dot(n0)))
	        .plus(n2.cross(n0).times(x1.dot(n1)))
	            .div(m.determinant())
        */
        var n0 = this.normal, n1 = plane.normal, n2 = n0.cross(n1).normalized()
        var p = M4.forRows(n0, n1, n2).inversed().transformVector(V3.create(this.w, plane.w, 0))
		return L3(p, n2)
    },

	/**
	 * Returns the point in the plane closest to the given point
	 *
	 * @param {V3} x
	 * @returns {V3}
	 */
    projectedPoint: function (x) {
	    // See http://math.stackexchange.com/questions/444968/project-a-point-in-3d-on-a-given-plane
	    // p = x - ((x - planeAnchor) * normal) * normal
	    return x.minus(this.normal.times(x.minus(this.anchor).dot(this.normal)))
    },

	/**
	 *
	 * @param {V3} x
	 * @returns {V3}
	 */
    projectedVector: function (x) {
	    // See V3.rejectedFrom. Simplified, as this.normal.length() == 1
	    return x.minus(this.normal.times(x.dot(this.normal)))
    },

	/**
	 *
	 * @returns {NLA.Plane3}
	 */
    flipped: function () {
	    return P3(this.normal.negated(), -this.w)
    }

}
NLA.addTransformationMethods(P3.prototype)
 // X-Y-Z planes
/** @const */ P3.YZ = P3(V3.X, 0)
/** @const */ P3.ZX = P3(V3.Y, 0)
/** @const */ P3.XY = P3(V3.Z, 0)

NLA.Plane3 = P3
