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
		    return P3(this.normal, this.anchor.plus(offset).dot(this.normal))
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

        // Returns the plane's distance from the given object (point, line or plane)
        distanceTo: function (obj) {
	        throw new Error("uhhh")
            if (this.intersects(obj) || this.contains(obj)) {
                return 0;
            }
            if (obj.anchor) {
                // obj is a plane or line
                var A = this.anchor.elements, B = obj.anchor.elements, N = this.normal.elements;
                return Math.abs((A[0] - B[0]) * N[0] + (A[1] - B[1]) * N[1] + (A[2] - B[2]) * N[2]);
            } else {
                // obj is a point
                var P = obj.elements || obj;
                var A = this.anchor.elements, N = this.normal.elements;
                return Math.abs((A[0] - P[0]) * N[0] + (A[1] - P[1]) * N[1] + (A[2] - (P[2] || 0)) * N[2]);
            }
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
			line.intersectWithPlane(this)
        },

        intersectionWithPlane: function (plane) {
            /*
             (x - a) * n = 0 // this
             (x - b) * m = 0 // plane
             x * n - a * n = 0 | * m0
             x * m - b * m = 0 | * n0

             (x1 * n1 + x2 * n2 - a * n) * m0 - (x1 * m1 + x2 * m2 - b * m) * n0= 0
             x1 * (n1 * m0 - m1 * n0) + x2 * (n2 * m0 - m2 * n0) -a * n * m0 + b * m * n0 | * (dir1 * m0 - m1 * dir0)

             x * dir - a * dir = 0
             ...
             x1 * (dir1 * m0 - m1 * dir0) + x2 * (dir2 * m0 - m2 * dir0) -a * dir * m0 + b * m * dir0 | * (n1 * m0 - m1 * n0)
             (x2 * (n2 * m0 - m2 * n0) -a * n * m0 + b * m * n0)  * (dir1 * m0 - m1 * dir0)-(x2 * (dir2 * m0 - m2 * dir0) -a * dir * m0 + b * m * dir0) * (n1 * m0 - m1 * n0)
             x2 * ((n2 * m0 - m2 * n0) * (dir1 * m0 - m1 * dir0) - (dir2 * m0 - m2 * dir0) * (n1 * m0 - m1 * n0))
             x2 * ((n2 * m0 * dir1 * m0 - m2 * n0 * dir1 * m0 - n2 * m0 * m1 * dir0) - (dir2 * m0 * n1 * m0 - m2 * dir0 * n1 * m0 - m1 * n0 * dir2 * m0))
             x2 * m0 (n2 * dir1 * m0 - m2 * n0 * dir1 - n2 * m1 * dir0 - dir2 * m0 * n1 + m2 * dir0 * n1 -m1 * n0 * dir2)
             x2 * m0 (dir0 * (- m1 * n2 + n1 * m2) + dir1 *(m0 * n2 -n0 * m2) + dir2 * (- m0 * n1 -m1 * n0)
             x2 * m0 * dir * (m X n)

            var direction = this.normal.cross(plane.normal), anchor;
            var f1 = this.normal.y * plane.normal.x - plane.normal.y * this.normal.x, f2 = this.normal.z * plane.normal.x - plane.normal.z * this.normal.x;
            var goal = plane.anchor.dot(plane.normal) * this.normal.x - this.anchor.dot(this.normal) * plane.normal.x;
            */
	        assert(plane instanceof P3, "plane instanceof P3")
	        assert(!this.isParallelToPlane(plane), "!this.isParallelToPlane(plane)")
	        var n0 = this.normal, n1 = plane.normal, n2 = n0.cross(n1).normalized(), m = M4.forSys(n0, n1, n2)
	        var x0 = this.anchor, x1 = plane.anchor, x2 = V3.ZERO
	        var p = n2.times(x2.dot(n2))
		        .plus(n1.cross(n2).times(x0.dot(n0)))
		        .plus(n2.cross(n0).times(x1.dot(n1)))
		            .div(m.determinant())
			return NLA.Line3(p, n2)
        },

        // Returns the point in the plane closest to the given point
	    projectedPoint: function (x) {
		    // See http://math.stackexchange.com/questions/444968/project-a-point-in-3d-on-a-given-plane
		    // p = x - ((x - planeAnchor) * normal) * normal
		    return x.minus(this.normal.times(x.minus(this.anchor()).dot(this.normal)))
	    },
	    projectedVector: function (x) {
		    return x.minus(this.normal.times(x.dot(this.normal)))
	    },

        // Returns the reflection of the plane in the given point, line or Plane3.
        reflectionIn: function (obj) {
            if (obj.normal) {
                // obj is a plane
                var A = this.anchor.elements, N = this.normal.elements;
                var A1 = A[0], A2 = A[1], A3 = A[2], N1 = N[0], N2 = N[1], N3 = N[2];
                var newA = this.anchor.reflectionIn(obj).elements;
                // Add the plane's normal to its anchor, then mirror that in the other plane
                var AN1 = A1 + N1, AN2 = A2 + N2, AN3 = A3 + N3;
                var Q = obj.pointClosestTo([AN1, AN2, AN3]).elements;
                var newN = [Q[0] + (Q[0] - AN1) - newA[0], Q[1] + (Q[1] - AN2) - newA[1], Q[2] + (Q[2] - AN3) - newA[2]];
                return Plane3.create(newA, newN);
            } else if (obj.direction) {
                // obj is a line
                return this.rotate(Math.PI, obj);
            } else {
                // obj is a point
                var P = obj.elements || obj;
                return Plane3.create(this.anchor.reflectionIn([P[0], P[1], (P[2] || 0)]), this.normal);
            }
        },
	    flipped: function () {
		    return P3(this.normal.negated(), -this.w)
	    }

    }
     // X-Y-Z planes
     P3.YZ = P3.ZY = P3(V3.X, 0)
     P3.ZX = P3.XZ = P3(V3.Y, 0)
     P3.XY = P3.YX = P3(V3.Z, 0)
    module.Plane3 = P3

})(NLA)