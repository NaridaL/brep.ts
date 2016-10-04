class AABB extends Transformable {
	min: V3
	max: V3

	constructor(min: V3 = V3.INF, max: V3 = V3.INF.negated()) {
		assertVectors(min, max)
		super()
		this.min = min
		this.max = max
	}

	addPoint(p: V3): this {
		assertVectors(p)
		this.min = this.min.min(p)
		this.max = this.max.max(p)
		return this
	}

	addPoints(ps: V3[]): this {
		ps.forEach(p => this.addPoint(p))
		return this
	}

	addAABB(aabb: AABB): this {
		assertInst(AABB, aabb)
		this.addPoint(aabb.min)
		this.addPoint(aabb.max)
		return this
	}

	/**
	 * Returns the largest AABB contained in this which doesn't overlap with aabb
	 * @param aabb
	 */
	withoutAABB(aabb) {
		assertInst(AABB, aabb)
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
	}

	getIntersectionAABB(aabb) {
		assertInst(AABB, aabb)
		return new AABB(this.min.max(aabb.min), this.max.min(aabb.max))
	}

	touchesAABB(aabb) {
		assertInst(AABB, aabb)
		return !(
		this.min.x > aabb.max.x || this.max.x < aabb.min.x
		|| this.min.y > aabb.max.y || this.max.y < aabb.min.y
		|| this.min.z > aabb.max.z || this.max.z < aabb.min.z)
	}

	intersectsAABB(aabb) {
		assertInst(AABB, aabb)
		return !(
		this.min.x >= aabb.max.x || this.max.x <= aabb.min.x
		|| this.min.y >= aabb.max.y || this.max.y <= aabb.min.y
		|| this.min.z >= aabb.max.z || this.max.z <= aabb.min.z)
	}

	intersectsAABB2d(aabb) {
		assertInst(AABB, aabb)
		return !(
		this.min.x >= aabb.max.x || this.max.x <= aabb.min.x
		|| this.min.y >= aabb.max.y || this.max.y <= aabb.min.y)
	}

	containsPoint(p) {
		assertVectors(p)
		return this.min.x <= p.x && this.min.y <= p.y && this.min.z <= p.z
			&& this.max.x >= p.x && this.max.y >= p.y && this.max.z >= p.z
	}

	containsSphere(center, radius) {
		assertVectors(center)
		assertNumbers(radius)
		return this.distanceToPoint(center) > radius
	}

	intersectsSphere(center, radius) {
		assertVectors(center)
		assertNumbers(radius)
		return this.distanceToPoint(center) <= radius
	}

	distanceToPoint(p) {
		assertVectors(p)
		var x = p.x, y = p.y, z = p.z
		var min = this.min, max = this.max
		if (this.containsPoint(p)) {
			return Math.max(
				min.x - x, x - max.x,
				min.y - y, y - max.y,
				min.z - z, z - max.z)
		}
		return p.distanceToPoint(new V3(
			NLA.clamp(x, min.x, max.x),
			NLA.clamp(y, min.y, max.y),
			NLA.clamp(z, min.z, max.z)))
	}

	containsAABB(aabb) {
		assertInst(AABB, aabb)
		return this.containsPoint(aabb.min) && this.containsPoint(aabb.max)
	}

	likeAABB(aabb) {
		assertInst(AABB, aabb)
		return this.min.like(aabb.min) && this.max.like(aabb.max)
	}

	intersectsLine(l3) {
		assertInst(L3, l3)
		var maxDim = l3.dir1.maxAbsDim()
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
	}

	hasVolume() {
		return this.min.x <= this.max.x && this.min.y <= this.max.y && this.min.z <= this.max.z
	}

	volume() {
		if (!this.hasVolume()) {
			return -1
		}
		var v = this.max.minus(this.min)
		return v.x * v.y * v.z
	}

	size() {
		return this.max.minus(this.min)
	}

	getCenter() {
		return this.min.plus(this.max).div(2)
	}

	transform(m4) {
		assertInst(M4, m4)
		assert(m4.isAxisAligned())
		var aabb = new AABB()
		aabb.addPoint(m4.transformPoint(this.min))
		aabb.addPoint(m4.transformPoint(this.max))
		return aabb
	}

	ofTransformed(m4) {
		assertInst(M4, m4)
		var aabb = new AABB()
		aabb.addPoints(m4.transformedPoints(this.corners()))
		return aabb
	}

	corners() {
		var min = this.min, max = this.max
		return [
			min,
			new V3(min.x, min.y, max.z),
			new V3(min.x, max.y, min.z),
			new V3(min.x, max.y, max.z),

			new V3(max.x, min.y, min.z),
			new V3(max.x, min.y, max.z),
			new V3(max.x, max.y, min.z),
			max
		]
	}

	toString() {
		return "new AABB(" + this.min.toString() + ", " + this.max.toString() + ")"
	}

	toSource() {
		return this.toString()
	}

	toMesh() {
		let matrix = M4.multiplyMultiple(
			M4.translation(this.min),
			M4.scaling(this.size().max(new V3(NLA_PRECISION, NLA_PRECISION, NLA_PRECISION))))
		console.log(matrix.str)
		console.log(matrix.inversed().transposed().str)
		let mesh = GL.Mesh.cube()
		console.log(mesh)
		// mesh.vertices = this.corners()
		mesh = mesh.transform(matrix)
		console.log(matrix.transformedPoints(mesh.vertices))
		mesh.compile()

		return mesh
	}
}
