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
 * @augments Curve
 */
function L3(anchor, dir1) {
    NLA.assertVectors(anchor, dir1)
    assert(dir1.hasLength(1), "dir must be normalized" + dir1)
	assertf(() => !Number.isNaN(anchor.x))
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
 * @param {V3} dir
 * @returns {L3}
 */
L3.anchorDirection = (anchor, dir) => L3(anchor, dir.normalized())

L3.prototype = Object.create(Curve.prototype)

NLA.addOwnProperties(L3.prototype, /** @lends L3.prototype */ {
	hlol: Curve.hlol++,
	constructor: L3,

    containsPoint: function (point) {
        NLA.assertVectors(point)
        var dist = this.distanceToPoint(point);
        assertNumbers(dist)
        return NLA.isZero(dist)
    },

    likeCurve: function (curve) {
    	return this.anchor.like(curve.anchor) && this.dir1.like(curve.dir1)
    },

	isColinearTo: function (obj) {
    	return obj instanceof L3
		    && this.containsPoint(obj.anchor)
		    && NLA.equals(1, Math.abs(this.dir1.dot(obj.dir1)))
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
	    var dirCross1 = this.dir1.cross(line.dir1).normalized()
		var anchorDiff = this.anchor.minus(line.anchor)
		return Math.abs(anchorDiff.dot(dirCross1))
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
	    assertInst(L3, line)
	    var dirCross = this.dir1.cross(line.dir1)
	    var div = dirCross.lengthSquared()
	    if (NLA.isZero(div)) { return null } // lines parallel
	    var anchorDiff = line.anchor.minus(this.anchor)
	    // check if distance is zero (see also L3.distanceToLine)
	    if (!NLA.isZero(anchorDiff.dot(dirCross.normalized()))) { return null }
	    var t = this.infoClosestToLine(line).t
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

	/**
	 * Returns true if the line is parallel (this.dir = line.dir || this.dir = -line.dir) to the argument.
	 *
	 * @param {L3} line
	 * @returns {boolean}
	 */
    isParallelToLine: function (line) {
        NLA.assertInst(L3, line)
        // we know that 1 == this.dir1.length() == line.dir1.length(), we can check for parallelity simpler than isParallelTo()
        return NLA.equals(1, Math.abs(this.dir1.dot(line.dir1)))
    },

	/**
	 *
	 * @param {V3} line
	 * @returns {number}
	 */
    angleToLine: function (line) {
        NLA.assertInst(L3, line)
        return this.dir1.angleTo(line.dir1)
    },

	/**
	 *
	 * @param {L3} line
	 * @returns {boolean} If the distance between the lines is zero
	 */
    intersectsLine: function (line) {
        return this.distanceToLine(line) <= NLA.PRECISION;
    },

	isInfosWithCurve: function (line) {
    	assertInst(L3, line)

		let dirCross = this.dir1.cross(line.dir1)
		let div = dirCross.lengthSquared()
		let anchorDiff = line.anchor.minus(this.anchor)
		if (NLA.isZero(anchorDiff.dot(dirCross))) {
			let tThis = anchorDiff.cross(this.dir1).dot(dirCross) / div
			let tOther = anchorDiff.cross(line.dir1).dot(dirCross) / div
			let p = this.at(tThis)
			return [{tThis: tThis, tOther: tOther, p: p}]
		}
		return []
	},

	/**
	 *
	 * @param {L3} line
	 * @returns {V3}
	 */
    isInfoWithLine: function (line) {
	    assertInst (L3, line)
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
    	// the two points on two lines the closest two each other are the ones whose
		// connecting
	    // TODO Where does this come from?
		// TODO: return value when no IS?
	    assertInst(L3, line)
	    var dirCross = this.dir1.cross(line.dir1)
	    var div = dirCross.lengthSquared()
	    var anchorDiff = line.anchor.minus(this.anchor)
	    var s = anchorDiff.cross(this.dir1).dot(dirCross) / div
	    var t = anchorDiff.cross(line.dir1).dot(dirCross) / div
	    return {s: s, t: t}
	    //console.log(segmentIntersectsRay, a, b, "ab", ab, "p", p, "dir", dir, s > 0 && t / div >= 0 && t / div <= 1, "s", s, "t", t, "div", div)
    },

	ddt: function (t) {
    	return V3.ZERO
	},

    toString: function (roundFunction) {
	    roundFunction = roundFunction || ((v) => +v.toFixed(4))
	    return "L3("+this.anchor.toString(roundFunction) + ", "+this.dir1.toString(roundFunction) +")"
    },

	/**
	 * Returns the point on the line that is closest to the given point.
	 *
	 * @param {V3} p
	 * @returns {V3}
	 */
    closestPointToPoint: function (p) {
		// similar logic as pointLambda; we project the vector (anchor -> p) onto dir1, then add anchor back to it
		let nearestT = p.minus(this.anchor).dot(this.dir1)
		return this.at(nearestT)
    },

	/**
	 *
	 * @param {L3} obj
	 * @returns {{t: number, s?: number, closest?: V3, closest2?: V3, distance: number}}
	 */
	infoClosestToLine: function (obj) {
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
		if (this.isParallelToLine(obj)) {
			return {t: NaN, s: NaN, distance: this.distanceToLine(obj)};
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
	},

    /**
     *
     * @param {P3} plane
     * @returns {V3}
     */
    intersectionWithPlane: function (plane) {
        // plane: plane.normal * p = plane.w
        // line: p=line.point + lambda * line.dir1
        var lambda = (plane.w - plane.normal.dot(this.anchor)) / plane.normal.dot(this.dir1);
        var point = this.anchor.plus(this.dir1.times(lambda));
        return point;
    },

    tangentAt: function (t) {
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

})

L3.prototype.hlol = Curve.hlol++

L3.pointLambdaNotNormalized = function (anchor, dir, x) {
	NLA.assertVectors(anchor, dir, x)
	return x.minus(anchor).dot(dir) / dir.lengthSquared()
}

/**
 *
 * @param {P3} p1
 * @param {P3} p2
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