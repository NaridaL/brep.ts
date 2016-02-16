/**
 * Created by aval on 20/11/2015.
 */
"use strict";
if (!NLA.Vector3) {
    throw new Error("Define NLA.V3 first")
}
(function (module) {
    var V3 = NLA.Vector3, P3 = NLA.Plane3
    var assert = NLA.assert, assertNumbers = NLA.assertNumbers
    

    var L3 = function (anchor, dir1) {
        NLA.assertVectors(anchor, dir1)
        assert(dir1.hasLength(1), "dir must be normalized" + dir1)
        var l = Object.create(L3.prototype)
        l.anchor = anchor
        l.dir1 = dir1
        return l
    }
    L3.throughPoints = (anchor, b) => L3(anchor, b.minus(anchor).normalized())
    L3.anchorDirection = (anchor, direction) => L3(anchor, direction.normalized())

    L3.prototype = {


        // Returns true if the line is parallel to the argument. Here, 'parallel to'
        // means that the argument's direction is either parallel or antiparallel to
        // the line's own direction. A line is parallel to a plane if the two do not
        // have a unique intersection.
        isParallelTo: function (obj) {
            if (obj.normal) {
                return obj.isParallelTo(this);
            }
            return this.dir1.isParallelTo(obj.dir1);
        },
        containsPoint: function (point) {
	        NLA.assertVectors(point)
            var dist = this.distanceToPoint(point);
            assertNumbers(dist)
            return NLA.isZero(dist)
        },
        equals: function (line) {
            NLA.assert(line instanceof L3);
            // we know that 1 == this.dir1.length() == line.dir1.length(), we can check for parallelity simpler than isParallelTo()
            return this.containsPoint(line.anchor) && NLA.equals(1, Math.abs(this.dir1.dot(line.dir1)))
        },

	    // Returns the line's perpendicular distance from the argument,
	    // which can be a point, a line or a plane
	    distanceTo: function (obj) {
		    if (obj.normal) {
			    return obj.distanceTo(this);
		    }
		    if (obj.dir1) {
			    // obj is a line
			    return this.distanceToLine(obj)
			    /*
			    if (this.isParallelTo(obj)) {
				    return this.distanceTo(obj.anchor);
			    }
			    var N = this.dir1.cross(obj.dir1).toUnitVector().elements;
			    var A = this.anchor.elements, B = obj.anchor.elements;
			    return Math.abs((A[0] - B[0]) * N[0] + (A[1] - B[1]) * N[1] + (A[2] - B[2]) * N[2]);
			    */
		    } else {
			    // obj is a point
			    return this.distanceToPoint(obj)
			    /*
			     // TODO what is it doing?
			     var P = obj.elements || obj;
			     var A = this.anchor.elements, D = this.dir1.elements;
			     var PA1 = P[0] - A[0], PA2 = P[1] - A[1], PA3 = (P[2] || 0) - A[2];
			     var modPA = Math.sqrt(PA1 * PA1 + PA2 * PA2 + PA3 * PA3);
			     if (modPA === 0) return 0;
			     // Assumes direction vector is normalized
			     var cosTheta = (PA1 * D[0] + PA2 * D[1] + PA3 * D[2]) / modPA;
			     var sin2 = 1 - cosTheta * cosTheta;
			     return Math.abs(modPA * Math.sqrt(sin2 < 0 ? 0 : sin2));
			     */
		    }
	    },

	    distanceToLine: function (line) {
			assert(line instanceof L3)
		    if (this.isParallelToLine(line)) {
			    return this.distanceToPoint(line.anchor)
		    }
		    var cross1 = this.dir1.cross(line.dir1).unit()
		    return Math.abs(this.anchor.minus(line.anchor).dot(cross1))
	    },

        distanceToPoint: function (x) {
            NLA.assertVectors(x)
            // See http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
            var t = x.minus(this.anchor).dot(this.dir1)
            return this.at(t).minus(x).length()

            //return x.minus(this.anchor).cross(x.minus(this.anchor.plus(this.dir1))).length()
        },
	    asSegmentDistanceToPoint: function (x, sStart, sEnd) {
		    var t = x.minus(this.anchor).dot(this.dir1)
		    t = NLA.clamp(t, sStart, sEnd)
		    return this.at(t).minus(x).length()
	    },
	    asSegmentDistanceToLine: function (line, sStart, sEnd) {
		    assert (line instanceof  L3)
		    var dirCross = this.dir1.cross(line.dir1)
		    var div = dirCross.lengthSquared()
		    if (NLA.isZero(div)) { return null } // lines parallel
		    var anchorDiff = line.anchor.minus(this.anchor)
		    // check if distance is zero (see also L3.distanceToLine)
		    if (!NLA.isZero(anchorDiff.dot(dirCross.normalized()))) { return null }
		    var t = this.pointClosestTo2(line).t
		    t = NLA.clamp(t, sStart, sEnd)
		    return this.at(NLA.clamp(t, sStart, sEnd))
	    },

	    at: function (lambda) {
		    assertNumbers(lambda)
		    return this.anchor.plus(this.dir1.times(lambda))
	    },
        /**
         * Every point x on this line is described by the equation x = this.anchor + lambda * this.dir1
         * This function returns lambda for a given point x
         * @param x
         */
        pointLambda: function (x) {
            NLA.assertVectors(x)
            var t = x.minus(this.anchor).dot(this.dir1)
	        return t
        },

        isParallelToLine: function (line) {
            NLA.assert(line instanceof L3)
            // we know that 1 == this.dir1.length() == line.dir1.length(), we can check for parallelity simpler than isParallelTo()
            return NLA.equals(1, Math.abs(this.dir1.dot(line.dir1)))
        },
        angleToLine: function (line) {
            NLA.assert(line instanceof L3)
            return this.dir1.angleTo(line.dir1)
        },

        // Returns true iff the argument is a point on the line
        contains: function (point) {
            var dist = this.distanceTo(point);
            return (dist !== null && dist <= NLA.PRECISION);
        },

        // Returns true iff the line lies in the given plane
        liesIn: function (plane) {
            return Plane3.contains(this);
        },

        // Returns true iff the line has a unique point of intersection with the argument
        intersects: function (obj) {
            if (obj.normal) {
                return obj.intersects(this);
            }
            return (!this.isParallelTo(obj) && this.distanceTo(obj) <= NLA.PRECISION);
        },

        // Returns the unique intersection point with the argument, if one exists
        intersectionWith: function (obj) {
            if (obj.normal) {
                return obj.intersectionWith(this);
            }
            if (!this.intersects(obj)) {
                return null;
            }
            var P = this.anchor.elements, X = this.dir1.elements,
                Q = obj.anchor.elements, Y = obj.dir1.elements;
            var X1 = X[0], X2 = X[1], X3 = X[2], Y1 = Y[0], Y2 = Y[1], Y3 = Y[2];
            var PsubQ1 = P[0] - Q[0], PsubQ2 = P[1] - Q[1], PsubQ3 = P[2] - Q[2];
            var XdotQsubP = -X1 * PsubQ1 - X2 * PsubQ2 - X3 * PsubQ3;
            var YdotPsubQ = Y1 * PsubQ1 + Y2 * PsubQ2 + Y3 * PsubQ3;
            var XdotX = X1 * X1 + X2 * X2 + X3 * X3;
            var YdotY = Y1 * Y1 + Y2 * Y2 + Y3 * Y3;
            var XdotY = X1 * Y1 + X2 * Y2 + X3 * Y3;
            var k = (XdotQsubP * YdotY / XdotX + XdotY * YdotPsubQ) / (YdotY - XdotY * XdotY);
            return Vector.create([P[0] + k * X1, P[1] + k * X2, P[2] + k * X3]);
        },
	    intersectionWithLine: function (line) {
		    assert (line instanceof  L3)
		    var dirCross = this.dir1.cross(line.dir1)
		    var div = dirCross.lengthSquared()
		    if (NLA.isZero(div)) { return null } // lines parallel
		    var anchorDiff = line.anchor.minus(this.anchor)
		    // check if distance is zero (see also L3.distanceToLine)
		    if (!NLA.isZero(anchorDiff.dot(dirCross.normalized()))) { return null }
		    var t = anchorDiff.cross(line.dir1).dot(dirCross) / div
		    return this.at(t)
	    },
	    intersectionLineST: function (line) {
		    // TODO Where does this come from?
		    assert (line instanceof  L3)
		    var abXdir = this.dir1.cross(line.dir1)
		    var div = abXdir.lengthSquared()
		    var anchorDiff = line.anchor.minus(this.anchor)
		    var s = anchorDiff.cross(this.dir1).dot(abXdir) / div
		    var t = anchorDiff.cross(line.dir1).dot(abXdir) / div
		    return {s: s, t: t}
		    //console.log(segmentIntersectsRay, a, b, "ab", ab, "p", p, "dir", dir, s > 0 && t / div >= 0 && t / div <= 1, "s", s, "t", t, "div", div)
	    },
	    toString: function (roundFunction) {
		    roundFunction = roundFunction || ((v) => +v.toFixed(4))
		    return "L3("+this.anchor.toString(roundFunction) + ", "+this.dir1.toString(roundFunction) +")"
	    },

        // Returns the point on the line that is closest to the given point or line
        pointClosestTo: function (obj) {
            if (obj.dir1) {
                // obj is a line
                if (this.isParallelTo(obj)) {
                    return null;
                }
                var D = this.dir1.elements, E = obj.dir1.elements;
                var D1 = D[0], D2 = D[1], D3 = D[2], E1 = E[0], E2 = E[1], E3 = E[2];
                // Create plane containing obj and the shared normal and intersect this with it
                // Thank you: http://www.cgafaq.info/wiki/Line-line_distance
                var x = (D3 * E1 - D1 * E3), y = (D1 * E2 - D2 * E1), z = (D2 * E3 - D3 * E2);
                var N = Vector.create([x * E3 - y * E2, y * E1 - z * E3, z * E2 - x * E1]);
                var P = Plane3.create(obj.anchor, N);
                return P.intersectionWith(this);
            } else {
                // obj is a point
                var P = obj.elements || obj;
                if (this.contains(P)) {
                    return Vector.create(P);
                }
                var A = this.anchor.elements, D = this.dir1.elements;
                var D1 = D[0], D2 = D[1], D3 = D[2], A1 = A[0], A2 = A[1], A3 = A[2];
                var x = D1 * (P[1] - A2) - D2 * (P[0] - A1), y = D2 * ((P[2] || 0) - A3) - D3 * (P[1] - A2),
                    z = D3 * (P[0] - A1) - D1 * ((P[2] || 0) - A3);
                var V = Vector.create([D2 * x - D3 * z, D3 * y - D1 * x, D1 * z - D2 * y]);
                var k = this.distanceTo(P) / V.modulus();
                return Vector.create([
                    P[0] + V.elements[0] * k,
                    P[1] + V.elements[1] * k,
                    (P[2] || 0) + V.elements[2] * k
                ]);
            }
        },
        pointClosestTo2: function (obj) {
            if (obj.dir1) {
                /*
                 line = a + s*b
                 this = c + t*d

                 (this - line) * b = 0
                 (this - line) * d = 0

                 (a + s*b - c - t*d) * b = 0
                 (a + s*b - c - t*d) * d = 0

                 (a - c + s*b - t*d) * b = 0
                 (a - c + s*b - t*d) * d = 0

                 (a - c)*b + (s*b - t*d)*b = 0
                 (a - c)*d + (s*b - t*d)*d = 0

                 (a - c)*b + s*(b*b) - t*(d*b) = 0
                 (a - c)*d + s*(b*d) - t*(d*d) = 0

                 s = (t*(d*b) - (a - c)*b) / (b*b)
                 =>
                 (a - c)*d + (t*(d*b) - (a - c)*b) / (b*b)*(b*d) - t*(d*d) = 0 | * (b*b)
                 (a - c)*d * (b*b) + (t*(d*b) - (a - c)*b)*(b*d) - t*(d*d) * (b*b) = 0
                 (a - c)*d * (b*b) + t*(d*b)*(b*d) - (a - c)*b*(b*d) - t*(d*d) * (b*b) = 0
                 t = ((a - c)*b*(b*d) - (a - c)*d * (b*b)) / ((d*b)*(b*d) - (d*d) * (b*b))
                 */
                // obj is a line
                if (this.intersects(obj)) {
                    return this.intersectionWithLine(obj);
                }
                if (this.isParallelTo(obj)) {
                    return {t: NaN, s: NaN, closest: null, distance: this.distanceTo(obj)};
                }
                var a = obj.anchor, b = obj.dir1, c = this.anchor, d = this.dir1;
                var bd = b.dot(d), bb = b.lengthSquared(), dd = d.lengthSquared(), amc = a.minus(c), divisor = bd * bd - dd * bb;
                var t = (amc.dot(b) * bd - amc.dot(d) * bb) / divisor;
                var s = (amc.dot(b) * dd - amc.dot(d) * bd) / divisor;
                return {
                    t: t,
                    s: s,
                    closest: this.at(t),
                    closest2: obj.at(s),
                    distance: this.at(t).minus(obj.at(s)).length()
                };
            } else {
                // obj is a point
                var P = obj.elements || obj;
                if (this.contains(P)) {
                    return Vector.create(P);
                }
                var A = this.anchor.elements, D = this.dir1.elements;
                var D1 = D[0], D2 = D[1], D3 = D[2], A1 = A[0], A2 = A[1], A3 = A[2];
                var x = D1 * (P[1] - A2) - D2 * (P[0] - A1), y = D2 * ((P[2] || 0) - A3) - D3 * (P[1] - A2),
                    z = D3 * (P[0] - A1) - D1 * ((P[2] || 0) - A3);
                var V = Vector.create([D2 * x - D3 * z, D3 * y - D1 * x, D1 * z - D2 * y]);
                var k = this.distanceTo(P) / V.modulus();
                return Vector.create([
                    P[0] + V.elements[0] * k,
                    P[1] + V.elements[1] * k,
                    (P[2] || 0) + V.elements[2] * k
                ]);
            }
        },

        intersectionWithPlane: function (plane) {
            // plane: plane.normal * p = plane.w
            // line: p=line.point + lambda * line.dir1
            var lambda = (plane.w - plane.normal.dot(this.anchor)) / plane.normal.dot(this.dir1);
            var point = this.anchor.plus(this.dir1.times(lambda));
            return point;
        },
	    tangentAt: function () {
		    return this.dir1
	    },
        intersectWithPlaneLambda: function (plane) {
            // plane: plane.normal * p = plane.w
            // line: p=line.point + lambda * line.dir1
	        var div = plane.normal.dot(this.dir1)
	        if (NLA.isZero(div)) return NaN
            var lambda = (plane.w - plane.normal.dot(this.anchor)) / div
            return lambda
        },
	    getIntersectionsWithPlane: function (plane) {
		    return [this.intersectWithPlaneLambda(plane)]
	    },
	    flipped: function () {
		    return L3(this.anchor, this.dir1.negated())
	    },

        transform: function (m4) {
            var newAnchor = m4.transformPoint(this.anchor)
            var newDir = m4.transformVector(this.dir1)
            return L3(newAnchor, newDir.normalized())
        },

	    closestPointOnLine: function (point) {
		    point = V3(point);
		    var t = point.minus(this.anchor).dot(this.dir1) / this.dir1.dot(this.dir1);
		    var closestpoint = this.anchor.plus(this.dir1.times(t));
		    return closestpoint;
	    },
		projectedOnPlane: function (plane) {
			assert(plane instanceof P3)
			return L3(plane.projectedPoint(this.anchor), plane.projectedVector(this.dir1).normalized())
		}

    }
	NLA.addTransformationMethods(L3.prototype)

    L3.fromPlanes = function (p1, p2) {
        assert(p1 instanceof P3)
        assert(p2 instanceof P3)
        var direction = p1.normal.cross(p2.normal);
        var l = direction.length();
        if (l < 1e-10) {
            throw new Error("Parallel planes");
        }
        direction = direction.times(1.0 / l);

	    return p1.intersectionWithPlane(p2)
    }

    Object.defineProperties(L3, {
        X: { value: L3(V3.ZERO, V3.X) },
        Y: { value: L3(V3.ZERO, V3.Y) },
        Z: { value: L3(V3.ZERO, V3.Z) },
    })

    module.Line3 = L3
})(NLA)