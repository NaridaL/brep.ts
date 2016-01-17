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
            return this.like(plane) || this.likeFlipped(plane)
        },
	    like: function (plane) {
		    NLA.assert(plane instanceof P3)
		    return NLA.equals(this.w, plane.w) && this.normal.like(plane.normal)
	    },
	    likeFlipped: function (plane) {
		    NLA.assert(plane instanceof P3)
		    return NLA.equals(this.w, -plane.w) && this.normal.like(plane.normal.negated())
	    },
        isParallelToPlane: function (plane) {
            NLA.assert(plane instanceof P3)
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
		    roundFunction = roundFunction || ((v) => +v.toFixed(3))
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
        // Returns the unique intersection with the argument, if one exists. The result
        // will be a vector if a line is supplied, and a line if a plane is supplied.
        intersectionWith: function (obj) {
            if (!this.intersects(obj)) {
                return null;
            }
            if (obj.direction) {
                // obj is a line
                var A = obj.anchor.elements, D = obj.direction.elements,
                    P = this.anchor.elements, N = this.normal.elements;
                var multiplier = (N[0] * (P[0] - A[0]) + N[1] * (P[1] - A[1]) + N[2] * (P[2] - A[2])) / (N[0] * D[0] + N[1] * D[1] + N[2] * D[2]);
                return Vector.create([A[0] + D[0] * multiplier, A[1] + D[1] * multiplier, A[2] + D[2] * multiplier]);
            } else if (obj.normal) {
                // obj is a plane
                var direction = this.normal.cross(obj.normal).toUnitVector();
                // To find an anchor point, we find one co-ordinate that has a value
                // of zero somewhere on the intersection, and remember which one we picked
                var N = this.normal.elements, A = this.anchor.elements,
                    O = obj.normal.elements, B = obj.anchor.elements;
                var solver = Matrix.Zero(2, 2), i = 0;
                while (solver.isSingular()) {
                    i++;
                    solver = Matrix.create([
                        [N[i % 3], N[(i + 1) % 3]],
                        [O[i % 3], O[(i + 1) % 3]]
                    ]);
                }
                // Then we solve the simultaneous equations in the remaining dimensions
                var inverse = solver.inverse().elements;
                var x = N[0] * A[0] + N[1] * A[1] + N[2] * A[2];
                var y = O[0] * B[0] + O[1] * B[1] + O[2] * B[2];
                var intersection = [
                    inverse[0][0] * x + inverse[0][1] * y,
                    inverse[1][0] * x + inverse[1][1] * y
                ];
                var anchor = [];
                for (var j = 1; j <= 3; j++) {
                    // This formula picks the right element from intersection by
                    // cycling depending on which element we set to zero above
                    anchor.push((i == j) ? 0 : intersection[(j + (5 - i) % 3) % 3]);
                }
                console.log("WAAAH", anchor, direction);
                return Line.create(new GL.Vector(anchor[0], anchor[1], anchor[2]), direction);
            }
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
            if (SYL.isZero(f1)) {
                if (SYL.isZero(f2)) {
                    anchor = new GL.Vector(plane.anchor.dot(plane.normal) / plane.normal.x, 0, 0);
                } else {

                }
            }
             */
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

        // Returns a copy of the plane, rotated by t radians about the given line
        // See notes on Line#rotate.
        rotate: function (t, line) {
            var R = Matrix.Rotation(t, line.direction).elements;
            var C = line.pointClosestTo(this.anchor).elements;
            var A = this.anchor.elements, N = this.normal.elements;
            var C1 = C[0], C2 = C[1], C3 = C[2], A1 = A[0], A2 = A[1], A3 = A[2];
            var x = A1 - C1, y = A2 - C2, z = A3 - C3;
            return Plane3.create([
                C1 + R[0][0] * x + R[0][1] * y + R[0][2] * z,
                C2 + R[1][0] * x + R[1][1] * y + R[1][2] * z,
                C3 + R[2][0] * x + R[2][1] * y + R[2][2] * z
            ], [
                R[0][0] * N[0] + R[0][1] * N[1] + R[0][2] * N[2],
                R[1][0] * N[0] + R[1][1] * N[1] + R[1][2] * N[2],
                R[2][0] * N[0] + R[2][1] * N[1] + R[2][2] * N[2]
            ]);
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

    };
     // X-Y-Z planes
     P3.YZ = P3.ZY = P3(V3.X, 0)
     P3.ZX = P3.XZ = P3(V3.Y, 0)
     P3.XY = P3.YX = P3(V3.Z, 0)
    module.Plane3 = P3

})(NLA)