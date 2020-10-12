import { Assert, outputLink, suite, test, testCurveCentralProjection } from './manager'

import { arrayFromFunction, arraySamples, lerp, M4, NLA_PRECISION, V, V3, Vector } from 'ts3dutils'
import { Mesh } from 'tsgl'
import { Curve, Edge, L3, P3, CustomPlane } from '..'

suite('L3', () => {
	only.test('isInfosWithLine', assert => {
	test('isInfosWithLine', assert => {
		const l1 = new L3(V3.Y, V3.X)
		const res = l1.isInfosWithLine(V3.X, V(0, 2))[0]
		assert.equal(res.tThis, 1)
		assert.equal(res.tOther, 0.5)
		assert.v3like(res.p, V(1, 1))
	})
	test('distanceToLine', assert => {
		assert.equal(L3.X.distanceToLine(new L3(V3.Z, V3.Y)), 1)
		assert.equal(L3.X.distanceToLine(new L3(V3.Z, V3.X)), 1)
	})
	test('isTsForLine', assert => {
		console.log(L3.X.isInfoWithLine(new L3(V(1, 1, 0), V3.Y)).sce)
		assert.ok(L3.X.isInfoWithLine(new L3(V(1, 1, 0), V3.Y)).equals(V3.X))
	})

	function testCurvePerspectiveTransform(assert: Assert, curve: Curve, m4: M4) {
		Mesh.prototype.compile = function() {
			return this
		}
		const cubeMesh = Mesh.cube().transform(M4.scale(2).translate(-1, -1, 2))
		const cubeLines = cubeMesh.LINES.map(i => cubeMesh.vertices[i])
		const minv = m4.inversed()
		const curveTransformed = (curve as any).transform4(m4) as Curve
		outputLink(assert, {
			edges: [Edge.forCurveAndTs(curve), Edge.forCurveAndTs(curveTransformed)],
			drLines: [...cubeLines, ...m4.transformedPoints(cubeLines)],
		})
	}

	test('transform4', assert => {
		const m4 = M4.perspective(45, 1, 1, 5)
		const line = new L3(V3.X, V3.Z, -10, -1)
		const lineT = line.transform4(m4)
		assert.ok(lineT.containsPoint(m4.transformPoint(line.at(-1))))
		assert.ok(lineT.containsPoint(m4.transformPoint(line.at(-10))))
	})

	test('transform4 central projection', assert => {
		const l = new L3(V(1, 1, 0), V(1, 1, 0.3).unit(), 0.1, 10)
		testCurveCentralProjection(assert, l)
	})
})
