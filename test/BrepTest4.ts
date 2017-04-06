QUnit.module('BrepTest4')




QUnit.assert.testIntersectFace = function (face, brep2, resultEdges, resultPoints, desc) {
	if (brep2 instanceof Face) {
		brep2 = new B2([brep2])
	}
	this.ok(true, `<html><a style='color: #0000ff text-decoration: underline' target='blank'
						href='file:///C:/Users/aval/Desktop/cs/brep2.html?a=${new B2([face]).toSource()}&b=${brep2.toSource()}
						&edges=[${resultEdges.map(e => e.toSource()).join(',')}]
						&points=[${resultPoints.map(e => e.toSource()).join(',')}]'>${desc}</a>`)
	const faceMap = new Map(), edgeMap = new Map()
	brep2.faces.forEach(face2 => {
		face.intersectPlaneFace(face2, new B2([face]), brep2, faceMap, edgeMap)
	})
	const edges = faceMap.get(face) || []
	console.log(faceMap)
	this.equal(edges.length, resultEdges.length, 'resultEdges.length == edges.length'+edges.toSource())
	resultEdges.forEach(edge => {
		this.ok(edges.some(edge2 => edge.like(edge2)), `edges.some(edge2 => edge.like(edge2)) [${edges.toSource()}] ${edge.toSource()}`)
	})
	const uniquePoints = []
	face.contour.forEach(edge => {
		const em = edgeMap.get(edge)
		em && em.forEach(info => info && !uniquePoints.some(up => up.like(info.p)) && assert(info.p) && uniquePoints.push(info.p))
	})
	this.equal(uniquePoints.length, resultPoints.length, 'points.length == resultPoints.length'+uniquePoints.toSource())
	resultPoints.forEach(p => {
		this.ok(uniquePoints.some(up => up.like(p)), `edges.some(edge2 => edge.like(edge2)) [${uniquePoints.toSource()}]`)
	})
}
{
    new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(10, 0, 0), V(8, 0, 0), 0, 1.9999999999999984)[{
        face: new PlaneFace(new PlaneSurface(new P3(V(0, -1, 0), 0)), [
            new StraightEdge(new L3(V(0, 0, 0), V(0.7071067811865476, 0, 0.7071067811865475)), V(0, 0, 0), V(3.381197846482994, 0, 3.3811978464829933), 0, 4.781735851562952),
            new StraightEdge(new L3(V(2.791322082992153, 0, 3.813016875511774), V(0.8068982213550735, 0, -0.5906904945688721)), V(3.381197846482994, 0, 3.3811978464829933), V(8, 0, 0), 0.7310411002024879, 6.45518577084059),
            new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(8, 0, 0), V(10, 0, 0), 1.9999999999999984, 0),
            new StraightEdge(new L3(V(10, 0, 0), V(0, 0, 1)), V(10, 0, 0), V(10, 0, 10), 0, 10),
            new StraightEdge(new L3(V(10, 0, 10), V(-1, 0, 0)), V(10, 0, 10), V(0, 0, 10), 0, 10),
            new StraightEdge(new L3(V(0, 0, 0), V(0, 0, 1)), V(0, 0, 10), V(0, 0, 0), 10, 0)], []),
        edge: new StraightEdge(new L3(V(10, 0, 0), V(-1, 0, 0)), V(8, 0, 0), V(10, 0, 0), 1.9999999999999984, 0),
        normalAtCanonA: V(0, -1, 0),
        reversed: true,
        inside: V3.Z,
        angle: 0
    }]
}
{
	QUnit.module('B2.assembleFacesFromLoops')
	function ccwSquare(x, y, width, height) { return StraightEdge.chain([V(x, y), V(x + width, y), V(x + width, y + height), V(x, y + height)], true) }
	function reversedLoop(loop) { return loop.map(edge => edge.flipped()).reverse() }
	const surface = new PlaneSurface(P3.XY)
	// 1(3(2), 4)
	const CCW1 = ccwSquare(0, 0, 10, 10) // outer loop
	const CCW2 = ccwSquare(5, 5, 1, 1), CW2 = reversedLoop(CCW2) // small in middle
	const CCW3 = ccwSquare(4, 4, 3, 3), CW3 = reversedLoop(CCW3) // around small in middle
	const CCW4 = ccwSquare(1, 1, 1, 1), CW4 = reversedLoop(CCW4)
	const CCW5 = ccwSquare(10, 10, 1, 1)

	function test(name, loops, expected) {
		QUnit.test(name, function (assert) {
			const actual = B2.assembleFacesFromLoops(loops, surface, PlaneFace)
			expected = expected.length ? expected : [expected]
			assert.equal(actual.length, expected.length)
			for (let i = 0; i < expected.length; i++) {
				assert.ok(actual[i].likeFace(expected[i]))

			}
		})
	}

	test('CCW: CCW', [CCW1], new PlaneFace(surface, CCW1, []))

	test('CCW1(CW2): CCW1(CW2)', [CCW1, CW2], new PlaneFace(surface, CCW1, [CW2]))

	test('CCW1(CCW4, CW2): CCW4', [CCW1, CCW4, CW2], new PlaneFace(surface, CCW4, []))

	test('CCW1(CCW3(CW2)): CCW3(CW2)', [CCW1, CCW3, CW2], new PlaneFace(surface, CCW3, [CW2]))

	test('CCW1, CCW5: CCW1, CCW5', [CCW1, CCW5], [new PlaneFace(surface, CCW1, []), new PlaneFace(surface, CCW5, [])])

	test('CCW1(CW4, CW2): CCW1(CW4, CW2)', [CCW1, CW4, CW2], new PlaneFace(surface, CCW1, [CW2, CW4]))

}