
var assert = NLA.assert, V3 = NLA.Vector3, M4 = NLA.Matrix4x4
function AABB(min = V3.INF, max = V3.INF.negated()) {
	assert(min instanceof V3, "min instanceof V3" + min.toSource())
	assert(max instanceof V3)
	this.min = min
	this.max = max
}
AABB.prototype = {
	addPoint: function (p) {
		assert(p instanceof V3)
		this.min = this.min.min(p)
		this.max = this.max.max(p)
		return this
	},
	addPoints: function (ps) {
		ps.forEach(p => this.addPoint(p))
		return this
	},
	addAABB: function (aabb) {
		assert(aabb instanceof AABB)
		this.addPoint(aabb.min)
		this.addPoint(aabb.max)
	},
	/**
	 * Returns the largest AABB contained in this which doesn't overlap with aabb
	 * @param aabb
	 */
	withoutAABB: function (aabb) {
		assert(aabb instanceof AABB)
		var min, max
		var volume = this.volume(), size = this.size()
		var remainingVolume = -Infinity
		for (var i = 0; i < 3; i++) {
			var dim = ["x", "y", "z"][i]
			var cond = aabb.min[dim] - this.min[dim] > this.max[dim] - aabb.max[dim]
			var dimMin = cond ? this.min[dim] : Math.max(this.min[dim], aabb.max[dim])
			var dimMax = !cond ? this.max[dim] : Math.min(this.max[dim], aabb.min[dim])
			var newRemainingVolume = (dimMax - dimMin) * volume / size[dim]
			if (newRemainingVolume > remainingVolume) {
				remainingVolume = newRemainingVolume
				min = this.min.withElement(dim, dimMin)
				max = this.max.withElement(dim, dimMax)
			}
		}
		return new AABB(min, max)
	},
	getIntersectionAABB: function (aabb) {
		assert(aabb instanceof AABB)
		return new AABB(this.min.max(aabb.min), this.max.min(aabb.max))
	},
	intersectsAABB: function (aabb) {
		assert(aabb instanceof AABB)
		return !(
			this.min.x > aabb.max.x || this.max.x < aabb.min.x
			|| this.min.y > aabb.max.y || this.max.y < aabb.min.y
			|| this.min.z > aabb.max.z || this.max.z < aabb.min.z)
	},
	containsPoint: function (p) {
		assert(p instanceof V3)
		return this.min.x <= p.x && this.min.y <= p.y && this.min.z <= p.z
			 && this.max.x >= p.x && this.max.y >= p.y && this.max.z >= p.z
	},
	containsSphere: function (center, radius) {
		assert(center instanceof V3)
		NLA.assertNumbers(radius)
		return this.distanceToPoint(center) > radius
	},
	intersectsSphere: function (center, radius) {
		assert(center instanceof V3)
		NLA.assertNumbers(radius)
		return this.distanceToPoint(center) <= radius
	},
	distanceToPoint: function (p) {
		assert(p instanceof V3)
		var x = p.x, y = p.y, z = p.z
		var min = this.min, max = this.max
		if (this.containsPoint(p)) {
			return Math.max(
				min.x - x, x - max.x,
				min.y - y, y - max.y,
				min.z - z, z - max.z)
		}
		return p.distanceToPoint(V3(
			NLA.clamp(x, min.x, max.x),
			NLA.clamp(y, min.y, max.y),
			NLA.clamp(z, min.z, max.z)))
	},
	containsAABB: function (aabb) {
		assert(aabb instanceof AABB)
		return this.containsPoint(aabb.min) && this.containsPoint(aabb.max)
	},
	likeAABB: function (aabb) {
		assert(aabb instanceof AABB)
		return this.min.like(aabb.min) && this.max.like(aabb.max)
	},
	intersectsLine: function (l3) {
		assert (l3 instanceof L3)
		var maxDim = l3.dir1.absMaxDim()
		var [coord0, coord1] = [['y', 'z'], ['z', 'x'], ['x', 'y']][maxDim]
		var s0 = (this.min[maxDim] - l3.anchor[maxDim]) / l3.dir1[maxDim]
		var s1 = (this.max[maxDim] - l3.anchor[maxDim]) / l3.dir1[maxDim]
		var sMin = Math.min(s0, s1)
		var sMax = Math.max(s0, s1)
		var c = l3.dir1[coord0] * l3.anchor[coord1] - l3.anchor[coord0] * l3.dir1[coord1]
		function lineSide(pCoord0, pCoord1) {
			return l3.dir1[coord1] * pCoord0 + l3.dir1[coord0] * pCoord1 + c
		}
		var sideBL = lineSide()
	},
	hasVolume: function () {
		return this.min.x <= this.max.x && this.min.y <= this.max.y && this.min.z <= this.max.z
	},
	volume: function () {
		if (!this.hasVolume()) {
			return -1
		}
		var v = this.max.minus(this.min)
		return v.x * v.y * v.z
	},
	size: function () {
		return this.max.minus(this.min)
	},
	getCenter: function () {
		return this.min.plus(this.max).div(2)
	},
	transform: function (m4) {
		assert(m4 instanceof M4, "m4 instanceof M4")
		assert(m4.isAxisAligned())
		var aabb = new AABB()
		aabb.addPoint(m4.transformPoint(this.min))
		aabb.addPoint(m4.transformPoint(this.max))
		return aabb
	},
	ofTransformed: function (m4) {
		assert(m4 instanceof M4, "m4 instanceof M4")
		var aabb = new AABB()
		aabb.addPoints(m4.transformedPoints(this.corners()))
		return aabb
	},
	corners: function () {
		var min = this.min, max = this.max
		return [
			min,
			V3(min.x, min.y, max.z),
			V3(min.x, max.y, min.z),
			V3(min.x, max.y, max.z),

			V3(max.x, min.y, min.z),
			V3(max.x, min.y, max.z),
			V3(max.x, max.y, min.z),
			max
		]
	},

	toString: function () {
		return "new AABB("+this.min.toString()+", "+this.max.toString()+")"
	}

	
}
NLA.addTransformationMethodsToPrototype(AABB.prototype)
AABB.intersections = function (aabbs) {

}