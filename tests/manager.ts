declare const require: any
try {
	;(global as any).WebGLRenderingContext = {}
	//const mock = require('mock-require')
	//mock('tsgl', {})
} catch (e) {}

export * from 'ts3dutils/tests/manager'
import {
	arrayFromFunction,
	assert,
	DEG,
	eq,
	eq0,
	glqInSteps,
	int,
	lerp,
	M4,
	NLA_PRECISION,
	raddd,
	toSource,
	V,
	V3,
	arraySamples,
} from 'ts3dutils'
import { Assert, test } from 'ts3dutils/tests/manager'
import { RenderObjects } from '../src/viewer'

import slug from 'slug'

function sanitizeFilename(s: string) {
	return slug(s.replace(/-/g, 'minus').replace(/\+/g, 'plus'), '_')
}

import {
	BRep,
	ConicSurface,
	Curve,
	Edge,
	EllipseCurve,
	Face,
	ImplicitSurface,
	L3,
	P3,
	ParametricSurface,
	PlaneSurface,
	PointVsFace,
	rotateCurve,
	StraightEdge,
	Surface,
	CustomPlane,
} from '..'

import * as fs from 'fs'

export function testCurveCentralProjection(assert: Assert, curve: Curve) {
	const m4 = M4.projectPlanePoint(V3.O, new P3(V3.Z, 1))
	const pls = arraySamples(curve.tMin, curve.tMax, 16).flatMap(t => {
		const p = curve.at(t)
		return [p, m4.transformPoint(p)]
	})
	let curveTransformed
	try {
		curveTransformed = (curve as any).transform4(m4) as Curve
	} catch (e) {}
	outputLink(assert, {
		edges: [Edge.forCurveAndTs(curve), ...(curveTransformed ? [Edge.forCurveAndTs(curveTransformed)] : [])],
		drPs: pls,
		drLines: pls,
		planes: [new CustomPlane(V3.Z, V3.X, V3.Y, 'Z=1', [0, 0, 0, 1], -100, 100, -100, 100)],
	})
	if (curveTransformed) {
		arraySamples(curve.tMin, curve.tMax, 4).forEach(t => {
			assert.ok(curveTransformed.containsPoint(m4.transformPoint(curve.at(t))))
		})
	}
}

export function b2equals(assert: Assert, actual: BRep, expected: BRep, message = '') {
	if (!(actual instanceof BRep)) {
		assert.push(false, typeof actual, BRep, 'actual is not a BRep')
		return
	}

	assert.push(
		actual.faces.length == expected.faces.length,
		actual.toSource(false),
		expected.toSource(false),
		`no of faces ${actual.faces.length} != ${expected.faces.length}`,
	)

	actual.faces.forEach(face => {
		if (!expected.faces.some(expectedFace => expectedFace.likeFace(face))) {
			assert.push(
				false,
				actual.toSource(false),
				expected.toSource(false),
				'Unexpected face in result:' + face.toSource(),
			)
		}
	})
}

export function bRepEqual(assert: Assert, actual: BRep, expected: BRep, message = '') {
	let actualTranslated
	//try {
	const x = expected.getAABB().max.x - actual.getAABB().min.x
	actualTranslated = actual.translate(x === -Infinity ? 0 : x + 1)
	//} catch (e) { }
	outputLink(assert, {
		a: expected,
		b: actualTranslated,
	})
	b2equals(assert, actual, expected)
}

export function testBRepAnd(assert: Assert, a: BRep, b: BRep, expected: BRep, message?: string) {
	return testBRepOp(assert, a, b, () => a.and(b), expected, message)
}

export function testBRepOp(
	assert: Assert,
	a: BRep,
	b: BRep,
	calculateActual: () => BRep,
	expected: BRep,
	message?: string,
) {
	let actual
	try {
		actual = calculateActual()
	} finally {
		if (actual) {
			const abWidth = a
				.getAABB()
				.addAABB(b.getAABB())
				.size().x
			outputLink(
				assert,
				{
					a,
					b,
					c: actual.translate(abWidth + 1).toSource(false),
					d: expected.translate(2 * (abWidth + 1)).toSource(false),
				},
				message,
			)
			b2equals(assert, actual, expected)
		} else {
			outputLink(assert, { a, b })
		}
	}
}

export function registerTests(o: { [key: string]: (assert: Assert) => void })
export function registerTests(moduleName: string, o: { [key: string]: (assert: Assert) => void })
export function registerTests(moduleName: any, o?: any) {
	if ('string' == typeof moduleName) {
		QUnit.module(moduleName)
	} else {
		o = moduleName
	}
	for (const key in o) {
		QUnit.test(key, o[key])
	}
}

export function makeLink(values: any) {
	return Object.getOwnPropertyNames(values)
		.map(name => {
			const val = values[name]
			return name + '=' + (typeof val == 'string' ? val : val.toSource())
		})
		.join(';')
}

export function outputLink(
	assert: Assert,
	values: { [K in keyof RenderObjects]?: RenderObjects[K] | string },
	msg = 'view',
) {
	const script =
		'TEST_NAME = ' +
		assert.getTestName().toSource() +
		'\n' +
		Object.getOwnPropertyNames(values)
			.map(name => {
				const val = values[name]
				return 'const ' + name + ' = ' + (typeof val == 'string' ? val : toSource(val))
			})
			.join('\n') +
		'\n' +
		'return {' +
		Object.keys(values).join(',') +
		'}'
	const o = sanitizeFilename(assert.getTestName() + '_' + msg) + '.html'
	fs.writeFileSync(__dirname + '/results/' + o, demoFile.replace('/*INSERT*/', script), 'utf8')
	assert.link('http://localhost:10001/tests/results/' + o, msg)
}

const demoFile = fs.readFileSync(__dirname + '/../viewer.html', 'utf8')

export function testISCurves(assert: Assert, surface1: Surface | P3, surface2: Surface | P3, curveCount: int) {
	surface1 instanceof P3 && (surface1 = new PlaneSurface(surface1))
	surface2 instanceof P3 && (surface2 = new PlaneSurface(surface2))
	let isCurves
	try {
		isCurves = surface1.isCurvesWithSurface(surface2)
	} finally {
		if (isCurves) {
			outputLink(assert, {
				mesh: `[${surface1}.toMesh(), ${surface2}.toMesh()]`,
				edges: isCurves.map(c => Edge.forCurveAndTs(c)),
			})

			assert.equal(isCurves.length, curveCount, 'number of curves = ' + curveCount)
			for (const curve of isCurves) {
				assert.ok(
					surface1.containsCurve(curve),
					'surface1.containsCurve(curve) ' + surface1.str + ' ' + curve.str,
				)
				assert.ok(
					surface2.containsCurve(curve),
					'surface2.containsCurve(curve) ' + surface2.str + ' ' + curve.str,
				)
				const t = lerp(curve.tMin, curve.tMax, 0.5),
					p = curve.at(t),
					dp = curve.tangentAt(t)
				assert.ok(surface1.containsPoint(p), 'surface1.containsPoint(curve.at(curve.tMin))')
				assert.ok(surface2.containsPoint(p), 'surface2.containsPoint(curve.at(curve.tMin))')

				const pN1 = surface1.normalP(p)
				const pN2 = surface2.normalP(p)
				const expectedTangent = pN1.cross(pN2)
				// expectedTangent can be zero if the surfaces just touch and dont cross each other
				//!expectedTangent.likeZero() && assert.ok(expectedTangent.isParallelTo(dp),
				// 'pN1.cross(pN2).isParallelTo(dp)') !expectedTangent.likeZero() && assert.ok(expectedTangent.dot(dp)
				// > 0, 'pN1.cross(pN2).dot(dp) > 0')
			}
		} else {
			outputLink(assert, {
				mesh: `[${surface1}.toMesh(), ${surface2}.toMesh()]`,
			})
			assert.ok(false, 'no isCurves returned: ' + isCurves)
		}
	}
	return isCurves
}

export function testPointT(assert: Assert, curve: Curve, p: V3, expectedT?: number, precision?: number) {
	outputLink(assert, {
		edges: [Edge.forCurveAndTs(curve)],
		drPs: [p],
	})
	const actualT = curve.pointT(p)
	if (undefined !== expectedT) {
		if (isNaN(expectedT)) {
			assert.push(isNaN(actualT), curve.pointT(p), expectedT, 'testing pointT')
		} else {
			assert.fuzzyEqual(actualT, expectedT, 'testing pointT', precision)
		}
	}
	if (!isNaN(actualT)) {
		assert.v3like(curve.at(actualT), p)
	}
}

/**
 * Tests that the passed loop is CCW on the passed surface, and that the reversed loop is CW.
 * @param assert
 * @param surface
 * @param loop
 */
export function testLoopCCW(assert: Assert, surface: Surface, loop: Edge[]) {
	const points = [loop[0].a, loop[0].atAvgT()]
	outputLink(
		assert,
		{
			mesh: surface.sce + '.toMesh()',
			edges: loop,
			drPs: points,
		},
		'testLoopCCW',
	)
	assert.ok(surface.edgeLoopCCW(loop))
	assert.ok(!surface.edgeLoopCCW(Edge.reversePath(loop)))
}

export function surfaceVolumeAndAreaTests(face: Face, msg = 'face', expectedVolume?: number) {
	const flippedFace = face.flipped()
	const faceMesh = face.toMesh()
	const faceMeshVol = faceMesh.calcVolume()

	test(msg + ' area', assert => {
		outputLink(assert, { mesh: face.toSource() + '.toMesh()', face: face })
		const actualArea = face.calcArea()
		const expectedArea = faceMeshVol.area
		assert.fuzzyEqual(actualArea, expectedArea, undefined, 0.05)
		console.log('OK! actual = ' + actualArea + ', expected = ' + expectedArea)
	})
	test(msg + ' flipped() area', assert => {
		outputLink(assert, { mesh: face.flipped().toSource() + '.toMesh()', face: face.flipped() })
		const actualArea = flippedFace.calcArea()
		const expectedArea = faceMeshVol.area
		assert.fuzzyEqual(actualArea, expectedArea, undefined, 0.05)
		console.log('OK! actual = ' + actualArea + ', expected = ' + expectedArea)
	})
	test(msg + ' volume', assert => {
		outputLink(assert, { mesh: face.toSource() + '.toMesh()', edges: face.allEdges })
		faceMesh.calcVolume()
		const actual = face.zDirVolume().volume,
			expected = undefined === expectedVolume ? faceMeshVol.volume : expectedVolume
		assert.fuzzyEqual(actual, expected, undefined, 0.05)
		console.log('OK! actual = ' + actual + ', expected = ' + expected + ', |dv| = ' + (actual - expected))
	})
	test(msg + ' flipped() volume', assert => {
		outputLink(assert, { mesh: flippedFace.toSource() + '.toMesh()', edges: flippedFace.allEdges })
		const actual = flippedFace.zDirVolume().volume
		const expected = undefined === expectedVolume ? -faceMeshVol.volume : -expectedVolume
		assert.fuzzyEqual(actual, expected, undefined, 0.05)
		console.log('OK! actual = ' + actual + ', expected = ' + expected + ', |dv| = ' + (actual - expected))
	})
	test(msg + ' centroid', assert => {
		const actual = flippedFace.zDirVolume()
		const expected = faceMeshVol.centroid
		outputLink(assert, {
			mesh: face.toSource() + '.toMesh()',
			drPs: [expected, actual.centroid],
			edges: face.getAllEdges(),
		})
		if (!eq0(actual.volume)) {
			// centroid doesn't make sense when volume is 0
			assert.v3like(actual.centroid, expected, undefined, expected.length() / 100)
			console.log(
				'OK! actual = ' +
					actual.centroid +
					', expected = ' +
					expected +
					', |dv| = ' +
					actual.centroid.distanceTo(expected),
			)
		}
	})
	test(msg + ' flipped() centroid', assert => {
		const actual = flippedFace.zDirVolume()
		const expected = faceMeshVol.centroid
		outputLink(assert, {
			mesh: face.toSource() + '.toMesh()',
			drPs: [expected, actual.centroid],
		})
		if (!eq0(actual.volume)) {
			// centroid doesn't make sense when volume is 0
			assert.v3like(actual.centroid, expected, undefined, expected.length() / 100)
		}
	})
	test(msg + ' aabb', assert => {
		const expected = faceMesh.getAABB()
		const actual = face.getAABB()
		assert.v3like(actual.min, expected.min, 'aabb.min', 0.1)
		assert.v3like(actual.max, expected.max, 'aabb.max', 0.1)
	})
}

export function testCurve(assert: Assert, curve: Curve, checkTangents = true, msg?: string) {
	const edge = Edge.forCurveAndTs(curve)
	const debugInfo = curve.debugInfo && curve.debugInfo()
	outputLink(
		assert,
		{
			edges: [
				edge,
				...(debugInfo && debugInfo.lines
					? arrayFromFunction(debugInfo.lines.length / 2, i => {
							const a = debugInfo.lines[i * 2],
								b = debugInfo.lines[i * 2 + 1]
							return !a.like(b) && StraightEdge.throughPoints(a, b)
					  }).filter(x => x)
					: []),
			],
			drPs: [edge.a, edge.b, ...((debugInfo && debugInfo.points) || [])],
		},
		msg,
	)
	const STEPS = 12
	for (let i = 0; i < STEPS; i++) {
		const t = lerp(curve.tMin, curve.tMax, i / (STEPS - 1))
		const p = curve.at(t)
		// check that pointT and containsPoint behave as expected
		assert.push(eq(t, curve.pointT(p)), curve.pointT(p), t, 't eq pointT(at(t)) for ' + t)
		assert.ok(curve.containsPoint(p), `containsPoint(at(t == ${t}) == ${p})`)

		// check that tangentAt() behaves correctly
		if (checkTangents) {
			const eps = t != curve.tMax ? 1e-8 : -1e-8
			const expectedTangent = curve
				.at(t + eps)
				.minus(p)
				.div(eps)
			const actualTangent = curve.tangentAt(t)
			assert.v3like(actualTangent, expectedTangent, undefined, 1e-3)
		}
	}

	// test curve length
	if (curve.arcLength !== Curve.prototype.arcLength) {
		const expected = glqInSteps(t => curve.tangentAt(t).length(), curve.tMin, curve.tMax, 4)
		const actual = curve.arcLength(curve.tMin, curve.tMax)
		assert.push(
			eq(expected, actual, 1e-6),
			expected,
			actual,
			'curve should have same length as the numerically calculated value',
		)
	}
}

export function suiteSurface(surface: Surface) {
	if (ParametricSurface.is(surface)) {
		test('ParametricSurface', assert => testParametricSurface(assert, surface))
	}
	if (ImplicitSurface.is(surface)) {
		test('ImplicitSurface', assert => testImplicitSurface(assert, surface))
	}
}
export function testParametricSurface(assert: Assert, surf: ParametricSurface) {
	const debugInfo = surf.debugInfo && surf.debugInfo()
	outputLink(
		assert,
		{
			mesh: `[${surf}.toMesh()]`,
			edges: [
				...(debugInfo && debugInfo.lines
					? arrayFromFunction(debugInfo.lines.length / 2, i => {
							const a = debugInfo.lines[i * 2],
								b = debugInfo.lines[i * 2 + 1]
							return !a.like(b) && StraightEdge.throughPoints(a, b)
					  }).filter(x => x)
					: []),
			],
			drPs: [...((debugInfo && debugInfo.points) || [])],
		},
		'view',
	)

	assert.ok(ParametricSurface.is(surf))

	// test equals
	const clone = new surf.constructor(...surf.getConstructorParameters())
	assert.ok(clone.equals(surf), 'clone.equals(surf)')
	assert.ok(clone.isCoplanarTo(surf), 'clone.isParallelTo(surf)')
	assert.ok(clone.like(surf), 'clone.like(surf)')

	const params = [V(0, 0), V(0, 1), V(1, 0), V(1, 1), V(0.25, 0.25), V(0.6, 0.25), V(0.25, 0.6), V(0.6, 0.7)].map(
		pm => new V3(lerp(surf.uMin, surf.uMax, pm.x), lerp(surf.vMin, surf.vMax, pm.y), 0),
	)
	const points = params.map(({ x, y }) => surf.pUV(x, y))
	const psFlipped = surf.flipped()
	for (let i = 0; i < points.length; i++) {
		const p = points[i]
		const pm = params[i]
		assert.ok(surf.containsPoint(p))

		// test dpdu and dpdv
		const eps = 0.0001
		const dpdu = surf.dpdu()(pm.x, pm.y)
		const dpdv = surf.dpdv()(pm.x, pm.y)
		const dpduNumeric = p.to(surf.pUV(pm.x + eps, pm.y)).div(eps)
		const dpdvNumeric = p.to(surf.pUV(pm.x, pm.y + eps)).div(eps)
		assert.v3like(dpduNumeric, dpdu, 'dpdu', 0.01)
		assert.v3like(dpdvNumeric, dpdv, 'dpdv', 0.01)
		const pmNormal = surf.normalUV(pm.x, pm.y)
		assert.ok(pmNormal.hasLength(1), 'pmNormal.hasLength(1)')
		const dpduXdpdv = dpdu.cross(dpdv)
		if (ImplicitSurface.is(surf)) {
			assert.ok(eq0(surf.implicitFunction()(p)))
		}
		const pm2 = surf.uvP(p)
		const pNormal = surf.normalP(p)
		const psFlippedUV = psFlipped.uvP(p)
		const psFlippedNormal = psFlipped.normalP(p)
		if (!dpdu.likeO() && !dpdv.likeO()) {
			assert.v3like(pm2, pm, 'pm == uvP(pUV(pm))')
			assert.v3like(
				psFlipped.pUV(psFlippedUV.x, psFlippedUV.y),
				p,
				'psFlipped.pUV(psFlippedUV.x, psFlippedUV.y) == p',
			)

			assert.v3like(pmNormal, pNormal)

			const computedNormal = dpduXdpdv.unit()
			assert.ok(computedNormal.angleTo(pNormal) < 5 * DEG)

			assert.v3like(psFlippedNormal, pNormal.negated(), 'pNormal == -psFlippedNormal')
		}

		// test pointFoot:
		const offsetPoint = p.plus(pNormal.toLength(0.5))
		const actualFoot = surf.pointFoot(offsetPoint)
		// the original params may not actually be the closest.
		//if (!actualFoot.like(pm) && surf.pUV(actualFoot.x, actualFoot.y).distanceTo(p) > 0.5) {
		//    assert.v3like(actualFoot, pm, 'distance to foot ' + surf.pUV(actualFoot.x, actualFoot.y).distanceTo(p))
		//}
		const offsetPoint2 = p.plus(pNormal.toLength(-0.5))
		const actualFoot2 = surf.pointFoot(offsetPoint2)
		// the original params may not actually be the closest.
		//if (!actualFoot2.like(pm) && surf.pUV(actualFoot2.x, actualFoot2.y).distanceTo(p) > 0.5) {
		//    assert.v3like(actualFoot2, pm, 'distance to foot ' + surf.pUV(actualFoot2.x, actualFoot2.y).distanceTo(p))
		//}
	}
	const matrices = [M4.mirror(P3.XY), M4.mirror(P3.YZ), M4.mirror(P3.ZX)]
	for (let mI = 0; mI < matrices.length; mI++) {
		const m = matrices[mI]
		for (let i = 0; i < points.length; i++) {
			const dpdu = surf.dpdu()(params[i].x, params[i].y)
			const dpdv = surf.dpdv()(params[i].x, params[i].y)
			const p = points[i],
				pNormal = surf.normalP(p)
			const normalMatrix = m
				.as3x3()
				.inversed()
				.transposed()
			const mNormal = normalMatrix.transformVector(pNormal)
			const mP = m.transformPoint(p)
			const mSurface = surf.transform(m)

			if (!dpdu.likeO() && !dpdv.likeO()) {
				assert.v3like(mSurface.normalP(mP), mNormal)
			}

			assert.ok(mSurface.containsPoint(mP))

			//const mPSFlipped = mSurface.flipped()
			//assert.ok(mPSFlipped.normalP(mP).negated().like(mNormal))
			//assert(mPSFlipped.normalP(mP).negated().like(mNormal))
		}
	}
}

export function testContainsCurve(assert: Assert, surface: Surface, curve: Curve, expected = true, msg?: string) {
	outputLink(
		assert,
		{
			mesh: surface.sce + '.toMesh()',
			edges: [Edge.forCurveAndTs(curve)],
		},
		msg,
	)
	assert.equal(surface.containsCurve(curve), expected, 'surface contains curve')
}

export function rotateEdge(edge: Edge, angle: raddd) {
	const surface = rotateCurve(edge.curve, undefined, undefined, angle, edge.deltaT() > 0)
	const edges = [
		edge,
		Edge.forCurveAndTs(EllipseCurve.semicircle(edge.b.lengthXY(), V(0, 0, edge.b.z), 0, angle)),
		edge.rotateZ(angle).flipped(),
		Edge.forCurveAndTs(EllipseCurve.semicircle(edge.a.lengthXY(), V(0, 0, edge.a.z), 0, angle)).flipped(),
	]
	return Face.create(surface, edges)
}

function testImplicitSurface(assert: Assert, surface: ImplicitSurface) {
	const EPS = 1e-8
	const testPoints = [
		V3.O.plus(V(0.2, 0, 0)), // V3.O fails on ellipsoidSurface
		V3.Y,
		V3.X,
		V3.Z.plus(V(0.2, 0, 0)),
		V3.XY,
		V3.XYZ,
		new V3(10, 10, 10),
		new V3(5, 6, 7),
	]
	for (const testPoint of testPoints) {
		const i = surface.implicitFunction()(testPoint)
		const didpGuess = testPoint.map((el, dim) => {
			const i2 = surface.implicitFunction()(testPoint.plus(V3.O.withElement(dim, EPS)))
			return (i2 - i) / EPS
		})
		const didp = surface.didp(testPoint)
		assert.v3like(didp, didpGuess, undefined, 1e-5)
	}
}

export function testCurveISInfos(
	assert: Assert,
	c1: Curve,
	c2: Curve,
	expectedCount: int,
	msg: string = 'view',
	f = 'isInfosWithCurve',
) {
	let intersections
	try {
		intersections = c1[f](c2).map(info => info.p)
		outputLink(assert, { edges: [c1, c2].map(c => Edge.forCurveAndTs(c)), drPs: intersections }, msg)
		assert.equal(
			intersections.length,
			expectedCount,
			`intersections.length == count: ${intersections.length} == ${expectedCount}`,
		)
		intersections.forEach((is, i) => {
			assert.ok(
				intersections.every((is2, j) => j == i || !is.like(is2)),
				is.sce + ' is not unique ' + intersections,
			)
			assert.ok(c1.containsPoint(is), `e1.containsPoint(is): ${c1.toSource()}.containsPoint(${is.sce},`)
			assert.ok(c2.containsPoint(is), `e2.containsPoint(is): ${c1.toSource()}.containsPoint(${is.sce},`)
		})
	} finally {
		!intersections && outputLink(assert, { edges: [c1, c2].map(c => Edge.forCurveAndTs(c)) }, msg)
	}
}

/**
 * Test intersections of a curve with a surface.
 * @param assert
 * @param curve
 * @param surface
 * @param tCount
 */
export function testISTs(assert: Assert, curve: Curve, surface: Surface | P3, tCount: int) {
	surface instanceof P3 && (surface = new PlaneSurface(surface))
	let ists
	try {
		ists = curve instanceof L3 ? surface.isTsForLine(curve) : curve.isTsWithSurface(surface)
		assert.equal(ists.length, tCount, 'number of isps = ' + tCount)
		for (const t of ists) {
			const p = curve.at(t)
			assert.ok(
				surface.containsPoint(p),
				'surface.containsPoint(p) ' +
					surface.str +
					' ' +
					p.str +
					' t: ' +
					t +
					(ImplicitSurface.is(surface) ? ' dist: ' + surface.implicitFunction()(p) : ''),
			)
		}
	} finally {
		outputLink(
			assert,
			{
				mesh: `[${surface}.toMesh()]`,
				edges: [Edge.forCurveAndTs(curve)],
				drPs: ists ? ists.map(t => curve.at(t)) : [],
			},
			(ists && ists.join(', ')) || 'view',
		)
	}
}

export function linkBRep(assert: Assert, hash: string, message = 'view') {
	const escapedHash = encodeURIComponent(hash.replace(/, /g, ',').replace(/(\n|\t)+/g, ''))
		.replace(/\(/g, '%28')
		.replace(/\)/g, '%29')
	assert.link('http://localhost:10001/viewer.html#' + escapedHash, message)
}

export function testLoopContainsPoint(assert: Assert, surface: Surface, loop: Edge[], p: V3, result: PointVsFace) {
	const ccwLoop = surface.edgeLoopCCW(loop) ? loop : Edge.reversePath(loop)
	outputLink(assert, {
		mesh: Face.create(surface, loop).sce + '.toMesh()',
		drPs: [p],
	})
	assert.equal(surface.loopContainsPoint(loop, p), result)
}

export const FONT_PATH = __dirname + '/../fonts'

export function testCurvesColinear(test: Assert, curve1: Curve, curve2: Curve): void {
	test.ok(curve1.isColinearTo(curve2))
	const t = (curve1.tMin + curve1.tMax) / 2
	test.notOk(
		curve1
			.translate(
				curve1
					.tangentAt(t)
					.getPerpendicular()
					.unit(),
			)
			.isColinearTo(curve2),
	)
	outputLink(test, { edges: [curve1, curve2].map(c => Edge.forCurveAndTs(c)) })
	for (let i = 0; i < 10; i++) {
		const t = lerp(curve1.tMin, curve1.tMax, i / 9)
		if (!curve2.containsPoint(curve1.at(t))) {
			test.ok(false)
		}
	}
}
