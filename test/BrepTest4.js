QUnit.module('BrepTest4')


/**
 * Cases
 *          1. Volumes do not touch.
 *          2. face/face Face surfaces intersect each other. Implies edges going through faces.
 *             e.g. box(5, 5, 5) - box(5, 5, 5).translate(1, 1, 1)
 *          3. face/edge Edge of one volume lies in a face of the other
 *             e.g. box(5, 5, 5) - box(3, 3, 3).rotateZ([0, 1, 2] * PI / 2).translate(0, 1, 1)
 *          4. edge/edge Two edges are colinear.
 *
 *
 *
 */

/**
 *
 * @param face
 * @param brep2
 * @param resultEdges
 * @param resultPoints
 * @param desc
 */


QUnit.assert.testIntersectFace = function (face, brep2, resultEdges, resultPoints, desc) {
	if (brep2 instanceof Face) {
		brep2 = new B2([brep2])
	}
	this.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
						href='file:///C:/Users/aval/Desktop/cs/brep2.html?a=${new B2([face]).toSource()}&b=${brep2.toSource()}
						&edges=[${resultEdges.map(e => e.toSource()).join(',')}]
						&points=[${resultPoints.map(e => e.toSource()).join(',')}]'>${desc}</a>`)
	var faceMap = new Map(), edgeMap = new Map()
	brep2.faces.forEach(face2 => {
		face.intersectPlaneFace(face2, new B2([face]), brep2, faceMap, edgeMap)
	})
	var edges = faceMap.get(face) || []
	console.log(faceMap)
	this.equal(edges.length, resultEdges.length, 'resultEdges.length == edges.length'+edges.toSource())
	resultEdges.forEach(edge => {
		this.ok(edges.some(edge2 => edge.like(edge2)), `edges.some(edge2 => edge.like(edge2)) [${edges.toSource()}] ${edge.toSource()}`)
	})
	var uniquePoints = []
	face.edges.forEach(edge => {
		var em = edgeMap.get(edge)
		em && em.forEach(info => info && !uniquePoints.some(up => up.like(info.p)) && assert(info.p) && uniquePoints.push(info.p))
	})
	this.equal(uniquePoints.length, resultPoints.length, 'points.length == resultPoints.length'+uniquePoints.toSource())
	resultPoints.forEach(p => {
		this.ok(uniquePoints.some(up => up.like(p)), `edges.some(edge2 => edge.like(edge2)) [${uniquePoints.toSource()}]`)
	})
}