/**
 * Created by aval on 20/11/2015.
 */
/**
 * 3-dimensional line
 * 
 * @param anchor
 * @param dir1
 * @class L3
 * @constructor
 * @property {V3} dir1 Normalized direction of the line.
 * @property {V3} anchor Anchor of the line.
 * @returns L3
 */
function L3(anchor, dir1) {
    NLA.assertVectors(anchor, dir1)
    assert(dir1.hasLength(1), "dir must be normalized" + dir1)
	assert(!Number.isNaN(anchor.x))
    var l = Object.create(L3.prototype)
    l.dir1 = dir1
	l.anchor = anchor
	return l
}

/**
 *
 * @param {V3} anchor
 * @param {V3} b
 * @returns {L3}
 */
L3.throughPoints = (anchor, b) => L3(anchor, b.minus(anchor).normalized())
/**
 *
 * @param {V3} anchor
 * @param {V3} direction
 * @returns {L3}
 */
L3.anchorDirection = (anchor, direction) => L3(anchor, direction.normalized())

L3.prototype = {
	constructor: L3,

	/**
	 * Returns true if the line is parallel to the argument. Here, 'parallel to' means that the argument's direction is
	 * either parallel or antiparallel to the line's own direction. A line is parallel to a plane if the two do not
	 * have a unique intersection.
	 *
	 * @param {L3|NLA.Plane3} obj
	 * @returns {boolean}
	 */
    isParallelTo: function (obj) {
        if (obj.normal) {
            return obj.isParallelTo(this)
        }
        return this.dir1.isParallelTo(obj.dir1)
    },


	/**
	 *
	 * @param {V3} point
	 * @returns {boolean}
	 */
    containsPoint: function (point) {
        NLA.assertVectors(point)
        var dist = this.distanceToPoint(point);
        assertNumbers(dist)
        return NLA.isZero(dist)
    },
    equals: function (line) {
        NLA.assertInst(L3, line);
        // we know that 1 == this.dir1.length() == line.dir1.length(), we can check for parallelity simpler than isParallelTo()
        return this.containsPoint(line.anchor) && NLA.equals(1, Math.abs(this.dir1.dot(line.dir1)))
    },

    // Returns the line's perpendicular distance from the argument,
    // which can be a point, a line or a plane
	/**
	 *
	 * @param {NLA.Plane3|L3|V3} obj
	 * @returns {number}
	 */
    distanceTo: function (obj) {
	    if (obj.normal) {
		    // obj is a plane
		    return obj.distanceTo(this);
	    }
	    if (obj.dir1) {
		    // obj is a line
		    return this.distanceToLine(obj)
	    } else {
		    // obj is a point
		    return this.distanceToPoint(obj)
	    }
    },

	/**
	 *
	 * @param {L3} line
	 * @returns {number}
	 */
    distanceToLine: function (line) {
		assertInst(L3, line)
	    if (this.isParallelToLine(line)) {
		    return this.distanceToPoint(line.anchor)
	    }
	    var cross1 = this.dir1.cross(line.dir1).unit()
	    return Math.abs(this.anchor.minus(line.anchor).dot(cross1))
    },

    /**
     *
     * @param {V3} x
     * @returns {number}
     */
    distanceToPoint: function (x) {
        NLA.assertVectors(x)
        // See http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
        var t = x.minus(this.anchor).dot(this.dir1)
        return this.at(t).distanceTo(x)

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

    /**
     *
     * @param {number} lambda
     * @returns {V3}
     */
    at: function (lambda) {
	    assertNumbers(lambda)
	    return this.anchor.plus(this.dir1.times(lambda))
    },
    /**
     * This function returns lambda for a given point x
     *
     * Every point x on this line is described by the equation
     *      x = this.anchor + lambda * this.dir1 | - this.anchor
     *      x - this.anchor = lambda * this.dir1 | DOT this.dir1
     *      (x - this.anchor) DOT this.dir1 = lambda (dir1² is 1 as |dir1| == 1)
     *
     *  @param {V3} x
     *  @returns {number}
     */
    pointLambda: function (x) {
        NLA.assertVectors(x)
        var t = x.minus(this.anchor).dot(this.dir1)
        return t
    },

    isParallelToLine: function (line) {
        NLA.assertInst(L3, line)
        // we know that 1 == this.dir1.length() == line.dir1.length(), we can check for parallelity simpler than isParallelTo()
        return NLA.equals(1, Math.abs(this.dir1.dot(line.dir1)))
    },
    angleToLine: function (line) {
        NLA.assertInst(L3, line)
        return this.dir1.angleTo(line.dir1)
    },

    // Returns true iff the line lies in the given plane
    liesIn: function (plane) {
        return P3.contains(this);
    },

    // Returns true iff the line has a unique point of intersection with the argument
    intersects: function (obj) {
        if (obj.normal) {
            return obj.intersects(this);
        }
        return (!this.isParallelTo(obj) && this.distanceTo(obj) <= NLA.PRECISION);
    },

	/**
	 *
	 * @param {L3} line
	 * @returns {V3}
	 */
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

	/**
	 * returns s and t with this.at(s) == line.at(t)
	 *
	 * @param {L3} line
	 * @returns {{s: number, t: number}}
	 */
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

	/**
	 * Returns the point on the line that is closest to the given point or line
	 * @param {L3|P3} obj
	 * @returns {V3}
	 */
    pointClosestTo: function (obj) {
        if (obj.dir1) {
            // obj is a line
            if (this.isParallelTo(obj)) {
                return null;
            }
            return this.pointClosestTo2(obj).closest
        } else {
            // obj is a point
            // similar logic as pointLambda; we project the vector (anchor -> p) onto dir1, then add anchor back to it
            let nearestT = p.minus(this.anchor).dot(this.dir1)
			return this.at(nearestT)
        }
    },

	/**
	 *
	 * @param {L3|P3|V3} obj
	 * @returns {{t: number, s?: number, closest?: V3, closest2?: V3, distance: number}}
	 */
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
                return {t: NaN, s: NaN, distance: this.distanceTo(obj)};
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
                distance: this.at(t).distanceTo(obj.at(s))
            };
        } else {
            // obj is a point
            /** @type V3 */ let p = obj
            let nearestT = p.minus(this.anchor).dot(this.dir1)
            let closest = this.at(nearestT)
            return {t: nearestT, closest: closest, distance: closest.distanceTo(p)}
        }
    },

    /**
     *
     * @param {NLA.Plane3} plane
     * @returns {V3}
     */
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
    isTsWithPlane: function (plane) {
	    return [this.intersectWithPlaneLambda(plane)]
    },
    flipped: function () {
	    return L3(this.anchor, this.dir1.negated())
    },

	/**
	 *
	 * @param {M4} m4
	 * @returns {L3}
	 */
    transform: function (m4) {
        var newAnchor = m4.transformPoint(this.anchor)
        var newDir = m4.transformVector(this.dir1)
        return L3(newAnchor, newDir.normalized())
    },

	projectedOnPlane: function (plane) {
		assertInst(P3, plane)
		return L3(plane.projectedPoint(this.anchor), plane.projectedVector(this.dir1).normalized())
	},
    debugToMesh: function(mesh, bufferName) {
	    mesh[bufferName] || mesh.addVertexBuffer(bufferName, bufferName)
	    mesh[bufferName].push(this.at(-100), this.at(100))
    }

}
NLA.addTransformationMethods(L3.prototype)

L3.pointLambdaNotNormalized = function (anchor, dir, x) {
	NLA.assertVectors(anchor, dir, x)
	return x.minus(anchor).dot(dir) / dir.lengthSquared()
}

/**
 *
 * @param {NLA.Plane3} p1
 * @param {NLA.Plane3} p2
 * @returns {V3}
 */
L3.fromPlanes = function (p1, p2) {
    assertInst(P3, p1, p2)
    var direction = p1.normal.cross(p2.normal);
    var l = direction.length();
    if (l < 1e-10) {
        throw new Error("Parallel planes");
    }

    return p1.intersectionWithPlane(p2)
}
/** @type {L3} */
L3.X = L3(V3.ZERO, V3.X)
/** @type {L3} */
L3.Y = L3(V3.ZERO, V3.Y)
/** @type {L3} */
L3.Z = L3(V3.ZERO, V3.Z)