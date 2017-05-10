
QUnit.module('brep')
QUnit.assert.brepEquals = function(actual, expected, message) {
	this.push(actual.equals(expected), actual.ss(), expected.ss(), message)
}
QUnit.assert.fuzzyEquals = function(actual, expected, message) {
	this.push(eq(actual, expected), actual, expected, message)
}

registerTests({

	'B2.prototype.clipped() 1: no vertex/? intersections, no parallel faces'(assert) {
		const a = new B2(null, null, [new B2.Face([V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)], P3.XY)])
		const b = B2T.tetrahedron(V(5, 5, 5), V(5, 5, -5), V(10, 12, 1), V(0, 12, 1))
		const result = new B2(null, null, [new B2.Face([
			V(8.57142857142857, 9.999999999999998, 0),
			V(4.999999999999999, 4.999999999999999, 0),
			V(1.4285714285714297, 9.999999999999998, 0),
			V(0, 10, 0), V(0, 0, 0), V(10, 0, 0), V(10, 10, 0)], P3.XY)])
		assert.brepEquals(a.clipped(b), result)
	},

	'B2.prototype.clipped() 2: like 1, but clipping B2 vertices are on a plane (not inside polygon)'(assert) {
		const a = new B2(null, null, [new B2.Face([V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)], P3.XY)])
		const b = B2T.tetrahedron(V(5, 5, 5), V(5, 5, -5), V(10, 12, 0), V(0, 12, 0))
		const result = new B2(null, null, [new B2.Face([V(8.57142857142857, 9.999999999999998, 0), V(4.999999999999999, 4.999999999999999, 0), V(1.4285714285714297, 9.999999999999998, 0), V(0, 10, 0), V(0, 0, 0), V(10, 0, 0), V(10, 10, 0)], P3.XY)])
		assert.brepEquals(a.clipped(b), result)
	},


	'B2.prototype.clipped() 3: clipping B2 vertex is on face edge'(assert) {
		const a = new B2(null, null, [new B2.Face([V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)], P3.XY)])
		const b = B2T.tetrahedron(V(5, 5, 5), V(5, 5, -5), V(9, 10, 0), V(0, 12, 1))
		const result = new B2(null, null, [new B2.Face([V(9, 10, 0), V(4.999999999999999, 5, 0), V(1.4285714285714297, 9.999999999999998, 0), V(0, 10, 0), V(0, 0, 0), V(10, 0, 0), V(10, 10, 0)], P3.XY)])
		assert.brepEquals(a.clipped(b), result)
	},

	'B2.prototype.clipped() 4: clipping B2 edge is on face edge (+ non-zero overlap)'(assert) {
		const a = new B2(null, null, [new B2.Face([V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)], P3.XY)])
		const b = B2T.tetrahedron(V(5, 5, 5), V(5, 5, -5), V(9, 10, 0), V(1, 10, 0))
		const result = new B2(null, null, [new B2.Face([V(9, 10, 0), V(4.999999999999999, 5, 0), V(1, 10, 0), V(0, 10, 0), V(0, 0, 0), V(10, 0, 0), V(10, 10, 0)], P3.XY)])
		assert.brepEquals(a.clipped(b), result)
	},

	'B2.prototype.clipped() 5: projected edge is on face edge (+ non-zero overlap)'(assert) {
		const a = new B2(null, null, [new B2.Face([V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)], P3.XY)])
		const b = B2T.tetrahedron(V(5, 5, 5), V(5, 5, -5), V(9, 11, 1), V(1, 10, 0))
		const result = new B2(null, null, [new B2.Face([V(8.333333333333332, 9.999999999999998, 0), V(5, 4.999999999999999, 0), V(1, 10, 0), V(0, 10, 0), V(0, 0, 0), V(10, 0, 0), V(10, 10, 0)], P3.XY)])
		assert.brepEquals(a.clipped(b), result)
	},

	'B2.prototype.clipped() 6: projected poly completely inside face, projected vertex touches inside of face'(assert) {
		const a = new B2(null, null, [new B2.Face([V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)], P3.XY)])
		const b = B2T.tetrahedron(V(5, 5, 5), V(5, 5, -5), V(9, 11, 1), V(1, 9, 0))
		const result = new B2(null, null, [new B2.Face([V(8.333333333333332, 9.999999999999998, 0), V(5, 4.999999999999999, 0), V(1, 9, 0), V(8.333333333333332, 9.999999999999998, 0), V(0, 10, 0), V(0, 0, 0), V(10, 0, 0), V(10, 10, 0)], P3.XY)])
		assert.brepEquals(a.clipped(b), result)
	},

	'B2.prototype.clipped() 7'(assert) {
		const a = new B2(null, null, [new B2.Face([V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)], P3.XY)])
		const b = B2T.tetrahedron(V(5, 5, 0), V(8, 8, 0), V(10, 10, 5), V(1, 9, 5))
		const result = a
		assert.brepEquals(a.clipped(b), result)
	},
	'B2.prototype.clipped() 8'(assert) {
		const a = new B2(null, null, [new B2.Face([V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)], P3.XY)])
		const b = B2T.tetrahedron(V(-1, 4, 0), V(-1, 6, 0), V(5, 12, 5), V(5, 12, -5))
		const result = new B2(null, null, [
			new B2.Face([V(3.0000000000000004, 10, 0), V(0, 10, 0), V(0, 6.999999999999999, 0)], P3.XY),
			new B2.Face([V(3.5, 10, 0), V(0, 5.333333333333334, 0), V(0, 0, 0), V(10, 0, 0), V(10, 10, 0)], P3.XY)])

		assert.brepEquals(a.clipped(b), result)
	},
	'B2.prototype.clipped() 9 overlapping vertices'(assert) {
		const a = new B2(null, null, [new B2.Face([V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)], P3.XY)])
		const b = B2T.tetrahedron(V(10, 10, 0), V(6, 7, 1), V(7, 9, 1), V(6, 7, -5)).newAll(
			B2T.tetrahedron(V(10, 10, 0), V(7, 6, 1), V(9, 7, 1), V(7, 6, -5))
		)
		const result = new B2(null, null, [new B2.Face([V(10, 10, 0), V(8.666666666666668, 6.833333333333332, 0), V(7, 6, 0), V(10, 10, 0), V(6, 7.000000000000001, 0), V(6.833333333333334, 8.666666666666666, 0), V(10, 10, 0), V(0, 10, 0), V(0, 0, 0), V(10, 0, 0)], P3.XY)])
		assert.brepEquals(a.clipped(b), result)
	},
	'B2.prototype.clipped() 10 no overlap'(assert) {
		const b = B2T.tetrahedron(V(2, 2, 8), V(3, 3, 8), V(3, 2, -2), V(2, 3, -2)).translate(5, 5, 0)
		const a = new B2(null, null, [new B2.Face([V(5, 0, 0), V(5, 5, 0), V(0, 5, 0), V(0, 0, 0)].reverse(), new P3(V3.Z.negated(), 0))])
		assert.brepEquals(a.clipped(b), a)
	},
	'B2.prototype.clipped() 11 a face corner touches inside of projected polygon edge'(assert) {
		const a = new B2(null, null, [new B2.Face([V(0, 0, 0), V(10, 0, 0), V(10, 10, 0), V(0, 10, 0)], P3.XY)])
		const b = B2T.tetrahedron(V(11, 1, 0), V(9, -1, 0), V(5, 5, 0), V(5, 5, 10))
		const result = new B2(null, null, [
			new B2.Face([V(8.333333333333334, 0, 0), V(5, 5, 0), V(10, 1.6666666666666687, 0), V(10, 10, 0), V(0, 10, 0), V(0, 0, 0)], P3.XY)])
		assert.brepEquals(a.clipped(b, false, true), result)
	},
	'B2.prototype.clipped() 12, projected poly on inside of face, edge lies on face edge in the middle'(assert) {
		const b = B2T.tetrahedron(V(1, 0, -1), V(1, 0, 8), V(1, 4, -1), V(4, 0, -1))
		const a = new B2(null, null, [new B2.Face([V(5, 0, 0), V(5, 5, 0), V(0, 5, 0), V(0, 0, 0)].reverse(), new P3(V3.Z.negated(), 0))])
		const result = new B2(null, null, [new B2.Face([V(1, 0, 0), V(0, 0, 0), V(0, 5, 0), V(5, 5, 0), V(5, 0, 0), V(3.6666666666666665, -3.537972128266059e-16, 0), V(1, 3.555555555555555, 0)], new P3(V(0, 0, -1), 0))])
		assert.brepEquals(a.clipped(b), result)
	},
	'B2.prototype.clipped() 13'(assert) {

		const a = B2.extrude([V(0, 0, 0), V(10, 0, 0), V(0, 10, 0)], P3.XY, V(0, 0, -5), 'ex0').flipped()
		//a = new B2(null, null, [a.faces[0]])
		let b = B2.extrude([V(0, 0, -6), V(0, 4, -6), V(0, 0, -2)], P3.YZ.flipped().translate(0, 0, 0), V(13, 0, 0), 'ex1').flipped()
		b = new B2(null, null, [b.faces[2]])

		let c = b.clipped(a)
		const result = new B2(null, null, [new B2.Face([V(0, 2.999999999999999, -5), V(7, 3, -5), V(10, 0, -2.0000000000000004), V(0, 0, -2)], new P3(V(0, -0.7071067811865476, -0.7071067811865476), 1.4142135623730954))])
		assert.brepEquals(b.clipped(a), result)
	},


})