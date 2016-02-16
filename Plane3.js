/**
 * Created by aval on 20/11/2015.
 */
"use strict";
if (!NLA.Vector3) {
    throw new Error("Define NLA.V3 first")
}
(function (module) {
    var V3 = NLA.Vector3

	/**
	 * Oriented plane, i.e. splits R^3 in half, with one half being "in front" of the plane.
	 * Leads to multiple comparisons: isCoplanarToPlane returns if the plane occupies the same space,
	 * like returns if the plane occupies the same space and has the same orientation
	 *
	 * Points x on the plane fulfill the equation: normal DOT x = w
	 *
	 * @param normal1 normalized plane normal
	 * @param w signed (rel to normal1) distance from the origin
	 * @param prototype object to set as prototype of the new Plane3 object. Defaults to Plane3.prototype
	 * @returns new Plane3 object
	 * @constructor
	 */
    var P3 = function (normal1, w, prototype) {
        NLA.assertVectors(normal1)
        NLA.assertNumbers(w)
        NLA.assert(normal1.hasLength(1), "normal1.hasLength(1)")
        var p = Object.create(prototype || P3.prototype)
        p.normal = normal1
        p.w = w
        return p
    }
    P3.throughPoints = function (a, b, c, prototype) {
        NLA.assertVectors(a, b, c)
        var n1 = b.minus(a).cross(c.minus(a)).normalized();
        return P3(n1, n1.dot(a), prototype)
    }
    P3.normalOnAnchor = function (normal, anchor, prototype) {
        NLA.assertVectors(normal, anchor)
        var n1 = normal.normalized()
        return P3(n1, n1.dot(anchor), prototype)
    }
    P3.forAnchorAndPlaneVectors = function (anchor, v0, v1, prototype) {
        NLA.assertVectors(anchor, v0, v1)
        return P3.normalOnAnchor(v0.cross(v1), anchor, prototype)
    }
    P3.prototype = {

        get anchor () {
            return this.normal.times(this.w)
        },
        // Returns true iff the plane occupies the same space as the argument
        isCoplanarToPlane: function (plane) {
	        NLA.assert(plane instanceof P3, "plane instanceof P3")
            return this.like(plane) || this.likeFlipped(plane)
        },
	    like: function (plane) {
		    NLA.assert(plane instanceof P3, "plane instanceof P3")
		    return NLA.equals(this.w, plane.w) && this.normal.like(plane.normal)
	    },
	    likeFlipped: function (plane) {
		    NLA.assert(plane instanceof P3, "plane instanceof P3")
		    return NLA.equals(this.w, -plane.w) && this.normal.like(plane.normal.negated())
	    },
        isParallelToPlane: function (plane) {
            NLA.assert(plane instanceof P3, "plane instanceof P3")
            return NLA.equals(1, Math.abs(this.normal.dot(plane.normal)))
        },

        isParallelToLine: function (line) {
            NLA.assert(line instanceof NLA.Line3)
            return NLA.isZero(this.normal.dot(line.dir1))
        },

        isPerpendicularToLine: function (line) {
            NLA.assert(line instanceof NLA.Line3)
            // this.normal || line.dir1
            return NLA.equals(1, Math.abs(this.normal.dot(line.normal)))
        },

        isPerpendicularToPlane: function (plane) {
            NLA.assert(plane instanceof P3)
            return NLA.isZero(this.normal.dot(plane.normal))
        },
	    toString: function (roundFunction) {
		    roundFunction = roundFunction || (v => v) //((v) => +v.toFixed(3))
		    return "P3("+this.normal.toString(roundFunction) + ", " + roundFunction(this.w) +")"
	    },
	    translated: function (offset) {
		    return P3(this.normal, offset.dot(this.normal))
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
		    NLA.assert(line instanceof NLA.Line3)
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
            NLA.assert(line instanceof NLA.Line3)
            return this.containsPoint(line.anchor) && this.isParallelToLine(line)
        },
	    distanceToPointSigned: function (point) {
		    NLA.assert (point instanceof V3)
		    return this.normal.dot(point) - this.w
	    },
	    distanceToPoint: function (point) {
		    NLA.assert (point instanceof V3)
		    return Math.abs(this.normal.dot(point) - this.w)
	    },

        intersectionWithLine: function (line) {
			line.intersectionWithPlane(this)
        },

        intersectionWithPlane: function (plane) {
            /*

                this: n0 * x = w0
                plane: n1 * x = w1
                plane perpendicular to both which goes through origin:
                    n2 := n0 X x1
                    n2 * x = 0
            */
	        assert(plane instanceof P3, "plane instanceof P3")
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
			return NLA.Line3(p, n2)
        },

        // Returns the point in the plane closest to the given point
	    projectedPoint: function (x) {
		    // See http://math.stackexchange.com/questions/444968/project-a-point-in-3d-on-a-given-plane
		    // p = x - ((x - planeAnchor) * normal) * normal
		    return x.minus(this.normal.times(x.minus(this.anchor()).dot(this.normal)))
	    },
	    projectedVector: function (x) {
		    // See Vector3.rejectedFrom. Simplified, as this.normal.length() == 1
		    return x.minus(this.normal.times(x.dot(this.normal)))
	    },

	    flipped: function () {
		    return P3(this.normal.negated(), -this.w)
	    }

    }
	NLA.addTransformationMethods(P3.prototype)
     // X-Y-Z planes
     P3.YZ = P3.ZY = P3(V3.X, 0)
     P3.ZX = P3.XZ = P3(V3.Y, 0)
     P3.XY = P3.YX = P3(V3.Z, 0)
    module.Plane3 = P3

})(NLA)