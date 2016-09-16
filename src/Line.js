/**
 * Created by aval on 16/03/2016.
 */
/**
 * n-dimensional line
 *
 * @param anchor
 * @param dir1
 * @constructor
 */
NLA.Line = function (anchor, dir1) {
	NLA.assertVectors(anchor, dir1)
	console.log("sadjlkasjd", dir1)
	assert(dir1.hasLength(1), "dir must be normalized")
	this.anchor = anchor
	this.dir1 = dir1
}
NLA.Line.throughPoints = (anchor, b) => new NLA.Line(anchor, b.minus(anchor).normalized())
NLA.Line.anchorDirection = (anchor, direction) => new NLA.Line(anchor, direction.normalized())
NLA.Line.prototype = {
	containsPoint: function () {
	},
	equals: function (line) {
	    assertInst(NLA.Line, line);
		// we know that 1 == this.dir1.length() == line.dir1.length(), we can check for parallelity simpler than isParallelTo()
		return this.contains(line.anchor) && NLA.equals(1, Math.abs(this.dir1.dot(line.dir1)))
	},
	distanceToPoint: function (x) {
		NLA.assertVectors(x)
		// See http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
		var t = x.minus(this.anchor).dot(this.dir1)
		return this.at(t).minus(x).length()

		//return x.minus(this.anchor).cross(x.minus(this.anchor.plus(this.dir1))).length()
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
		assertInst(NLA.Line, line)
		// we know that 1 == this.dir1.length() == line.dir1.length(), we can check for parallelity simpler than isParallelTo()
		return NLA.equals(1, Math.abs(this.dir1.dot(line.dir1)))
	},
	angleToLine:  function (line) {
		assertInst(NLA.Line, line)
		return this.dir1.angleTo(line.dir1)
	},
	/**
	 *
	 * @param t
	 * @returns `this.anchor.plus(this.dir1.times(t))`
	 */
	at: function (t) {
		NLA.assertNumbers(t)
		return this.anchor.plus(this.dir1.times(t))
	},
	/**
	 *
	 * @param {L3|V3} obj
	 * @returns {*}
	 */
	pointClosestTo2: function(obj) {
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
			if (this.intersects(obj)) { return this.intersectionWith(obj); }
			if (this.isParallelTo(obj)) { return {t: NaN, s:NaN, closest: null, distance: this.distanceTo(obj)}; }
			var a = obj.anchor, b = obj.direction, c = this.anchor, d = this.direction;
			var bd = b.dot(d), bb = b.dot(b), dd = d.dot(d), amc = a.subtract(c), divisor = bd*bd - dd * bb;
			var t = (amc.dot(b)*bd - amc.dot(d)*bb) / divisor;
			var s = (amc.dot(b)*dd - amc.dot(d)*bd) / divisor;
			return {t: t, s:s, closest: this.at(t), closest2: obj.at(s), distance: this.at(t).subtract(obj.at(s)).modulus()};
		} else {
			/** @type V3 */ let p = obj
			if (this.contains(P)) { return Vector.create(P); }
			var A = this.anchor.elements, D = this.direction.elements;
			var D1 = D[0], D2 = D[1], D3 = D[2], A1 = A[0], A2 = A[1], A3 = A[2];
			var x = D1 * (P[1]-A2) - D2 * (P[0]-A1), y = D2 * ((P[2] || 0) - A3) - D3 * (P[1]-A2),
				z = D3 * (P[0]-A1) - D1 * ((P[2] || 0) - A3);
			var V = Vector.create([D2 * x - D3 * z, D3 * y - D1 * x, D1 * z - D2 * y]);
			var k = this.distanceTo(P) / V.modulus();
			return Vector.create([
				P[0] + V.elements[0] * k,
				P[1] + V.elements[1] * k,
				(P[2] || 0) + V.elements[2] * k
			]);
		}
	},
}