declare const require: any
try {
	(global as any).WebGLRenderingContext = {}
	//const mock = require('mock-require')
	//mock('tsgl', {})
} catch (e) { }

export * from 'ts3dutils/tests/manager'
import { arrayFromFunction, assert, DEG, eq, eq0, eq, glqInSteps, int, lerp, M4, V, V3, NLA_PRECISION } from 'ts3dutils'
import { Assert, test } from 'ts3dutils/tests/manager'

import slug from 'slug'
function sanitizeFilename(s: string) {
	return slug(s.replace(/-/g, 'minus').replace(/\+/g, 'plus'), '_')
}

import {
	BRep, ConicSurface, Curve, Edge, Face, ImplicitSurface, L3, P3, ParametricSurface, PlaneSurface,
	PointVsFace, Surface,
} from '..'

import * as fs from 'fs'

export function b2equals(assert: Assert, actual, expected, message = '') {
	if (!(actual instanceof BRep)) {
		assert.push(false, typeof actual, BRep, 'actual is not a BRep')
		return
	}

	assert.equal(actual.faces.length, expected.faces.length, 'no of faces')

	actual.faces.forEach(face => {
		if (!expected.faces.some(expectedFace => expectedFace.likeFace(face))) {
			assert.push(false, actual.toSource(false), expected.toSource(false), 'Unexpected face in result:' + face.toSource())
		}
	})
}

export function b2EqualAnd(test, a: BRep, b: BRep, expected: BRep) {
	return b2Equal(test, a, b, () => a.and(b), expected)
}
export function b2Equal(test, a: BRep, b: BRep, calculateActual: () => BRep, expected: BRep) {
	let actual
	try {
		actual = calculateActual()
	} finally {
		if (actual) {
			const abWidth = a.getAABB().addAABB(b.getAABB()).size().x
			linkB3(test, {
				a,
				b,
				c: actual.translate(abWidth + 1).toSource(false),
				d: expected.translate(2 * (abWidth + 1)).toSource(false),
			})
			b2equals(test, actual, expected)
		} else {
			linkB3(test, { a, b })
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
	return Object.getOwnPropertyNames(values).map(name => {
		const val = values[name]
		return name + '=' + (typeof val == 'string' ? val : val.toSource())
	}).join(';')
}

export function linkB3(assert: Assert, values, msg = 'view') {
	const script = 'TEST_NAME = ' + assert.getTestName().toSource() + '\n' +
		'FROM_SCRIPT = true\n' +
		Object.getOwnPropertyNames(values).map(name => {
			const val = values[name]
			return name + ' = ' + (typeof val == 'string' ? val : val.toSource())
		}).join('\n')
	const o = sanitizeFilename(assert.getTestName() + '_' + msg) + '.html'
	fs.writeFileSync('results/' + o, demoFile.replace('/*INSERT*/', script), 'utf8')
	// linkBRep(assert, makeLink(values), msg)
	assert.link('http://localhost:10001/tests/results/' + o)
}
const demoFile = fs.readFileSync('../viewer.html', 'utf8')
export function testISCurves(assert: Assert, surface1: Surface | P3, surface2: Surface | P3, curveCount: int) {
	surface1 instanceof P3 && (surface1 = new PlaneSurface(surface1))
	surface2 instanceof P3 && (surface2 = new PlaneSurface(surface2))
	let isCurves
	try {
		isCurves = surface1.isCurvesWithSurface(surface2)
	} finally {
		if (isCurves) {
		    linkB3(assert, {mesh: `[${surface1}.toMesh(), ${surface2}.toMesh()]`, edges: isCurves.map(c => Edge.forCurveAndTs(c))})
//			const script =
//				`var FROM_SCRIPT = true
//mesh=[${surface1}.toMesh(), ${surface2}.toMesh()]
//edges=${isCurves.map(c => Edge.forCurveAndTs(c)).sce}`
//			fs.writeFileSync('results/' + sanitizeFilename(assert.getTestName()) + '.html', demoFile.replace('/*INSERT*/', script), 'utf8')
//			linkBRep(assert, `mesh=[${surface1}.toMesh(), ${surface2}.toMesh()];edges=${isCurves.map(c => Edge.forCurveAndTs(c)).sce}`)

			assert.equal(isCurves.length, curveCount, 'number of curves = ' + curveCount)
			for (const curve of isCurves) {
				assert.ok(surface1.containsCurve(curve), 'surface1.containsCurve(curve) ' + surface1.str + ' ' + curve.str)
				assert.ok(surface2.containsCurve(curve), 'surface2.containsCurve(curve) ' + surface2.str + ' ' + curve.str)
				const t = curve.tMin || 0.2, p = curve.at(t), dp = curve.tangentAt(t)
				assert.ok(surface1.containsPoint(p), 'surface1.containsPoint(curve.at(curve.sMin))')
				assert.ok(surface2.containsPoint(p), 'surface2.containsPoint(curve.at(curve.tMax))')

				const pN1 = surface1.normalP(p)
				const pN2 = surface2.normalP(p)
				const expectedTangent = pN1.cross(pN2)
				// expectedTangent can be zero if the surfaces just touch and dont cross each other
				//!expectedTangent.likeZero() && assert.ok(expectedTangent.isParallelTo(dp),
				// 'pN1.cross(pN2).isParallelTo(dp)') !expectedTangent.likeZero() && assert.ok(expectedTangent.dot(dp)
				// > 0, 'pN1.cross(pN2).dot(dp) > 0')
			}
		} else {
			linkBRep(assert, `mesh=[${surface1}.toMesh(), ${surface2}.toMesh()]`)
		}
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
	linkB3(assert, {
	    mesh: surface.sce + '.toMesh()',
	    edges: loop,
        points: points,
    }, 'testLoopCCW')
	assert.ok(surface.edgeLoopCCW(loop))
	assert.ok(!surface.edgeLoopCCW(Edge.reversePath(loop)))
}

export function testZDirVolumeAndArea(assert: Assert, face: Face) {
	linkBRep(assert, `mesh=${face.sce}.toMesh()`)
	const faceMeshVol = face.toMesh().calcVolume()
	const actual = face.zDirVolume().volume, expected = faceMeshVol.volume
	assert.push(eq(actual, expected, 0.1), actual, expected, 'diff = ' + (actual - expected))

	const actualArea = face.calcArea()
	const expectedArea = faceMeshVol.area
	assert.push(eq(actualArea, expectedArea, 0.1), actualArea, expectedArea, 'diff = ' + (actualArea - expectedArea))
}

export function surfaceVolumeAndArea(face: Face) {
    return () => {
        const flippedFace = face.flipped()
        const faceMeshVol = face.toMesh().calcVolume()

        test('face volume', assert => {
            linkB3(assert, {mesh: face.toSource()+'.toMesh()'})
            const actual = face.zDirVolume().volume, expected = faceMeshVol.volume
            assert.push(eq(actual, expected, 0.1), actual, expected, 'diff = ' + (actual - expected))
        })
        test('face.flipped() volume', assert => {
            const actual = flippedFace.zDirVolume().volume, expected = -faceMeshVol.volume
            assert.push(eq(actual, expected, 0.1), actual, expected, 'diff = ' + (actual - expected))
        })
        test('face area', assert => {
            const actualArea = face.calcArea()
            const expectedArea = faceMeshVol.area
            assert.push(eq(actualArea, expectedArea, 0.1), actualArea, expectedArea, 'diff = ' + (actualArea - expectedArea))
        })
        test('face.flipped() area', assert => {
            const actualArea = flippedFace.calcArea()
            const expectedArea = faceMeshVol.area
            assert.push(eq(actualArea, expectedArea, 0.1), actualArea, expectedArea, 'diff = ' + (actualArea - expectedArea))
        })
    }
}

export function testCurve(ass: Assert, curve: Curve) {
	const STEPS = 12
	arrayFromFunction(STEPS, i => {
		const t = lerp(curve.tMin, curve.tMax, i / (STEPS - 1))
		const p = curve.at(t)
		// check that pointT and containsPoint behave as expected
		ass.push(
			eq(t, curve.pointT(p)),
			curve.pointT(p),
			t,
			't eq pointT(at(t)) for ' + t,
		)
		ass.ok(curve.containsPoint(p), `containsPoint(at(t == ${t}) == ${p})`)
	})

	// test curve length
	if (curve.arcLength !== Curve.prototype.arcLength) {
		const expected = glqInSteps(t => curve.tangentAt(t).length(), curve.tMin, curve.tMax, 4)
		const actual = curve.arcLength(curve.tMin, curve.tMax)
		ass.push(
			eq(expected, actual, 1e-6),
			expected,
			actual,
			'curve should have same length as the numericaly calculated value',
		)
	}
}

export function testParametricSurface(ass: Assert, surf: ParametricSurface) {
	linkBRep(ass, `mesh=[${surf}.toMesh()]`, 'view')
	const params = [V(0.25, 0.25), V(0.6, 0.25), V(0.25, 0.6), V(0.6, 0.7)]
		.map(pm => new V3(lerp(surf.sMin, surf.sMax, pm.x), lerp(surf.tMin, surf.tMax, pm.y), 0))
	const points = params.map(({ x, y }) => surf.pST(x, y))
	const psFlipped = surf.flipped()
	for (let i = 0; i < points.length; i++) {
		const p = points[i], pNormal = surf.normalP(p)
		const pm = params[i]
		ass.ok(surf.containsPoint(p))
		assert(surf.containsPoint(p))

		const psFlippedNormal = psFlipped.normalP(p)
		ass.ok(psFlippedNormal.negated().like(pNormal))
		assert(psFlippedNormal.negated().like(pNormal))

		const pm2 = surf.stP(p)
		ass.ok(pm.like(pm2))
		assert(pm.like(pm2))

		if (ParametricSurface.is(surf)) {
			const eps = 0.0001
			const dpds = surf.dpds()(pm.x, pm.y)
			const dpdt = surf.dpdt()(pm.x, pm.y)
			const dpdsNumeric = p.to(surf.pST(pm.x + eps, pm.y)).div(eps)
			const dpdtNumeric = p.to(surf.pST(pm.x, pm.y + eps)).div(eps)
			assert(dpdsNumeric.angleTo(dpds) < 5 * DEG)
			assert(dpdtNumeric.angleTo(dpdt) < 5 * DEG)
			const pmNormal = surf.normalST(pm.x, pm.y)
			assert(pmNormal.hasLength(1))
			ass.ok(pNormal.like(pmNormal))
			assert(pNormal.like(pmNormal))
			const computedNormal = dpds.cross(dpdt).unit()
			assert(computedNormal.angleTo(pNormal) < 5 * DEG)
		}
		if (ImplicitSurface.is(surf)) {
			ass.ok(eq0(surf.implicitFunction()(p)))
			assert(eq0(surf.implicitFunction()(p)))
		}
	}
	const matrices = [M4.mirror(P3.XY), M4.mirror(P3.YZ), M4.mirror(P3.ZX)]
	for (let mI = 0; mI < matrices.length; mI++) {
		const m = matrices[mI]
		for (let i = 0; i < points.length; i++) {
			const p = points[i], pNormal = surf.normalP(p)
			const normalMatrix = m.as3x3().inversed().transposed()
			const mNormal = normalMatrix.transformVector(pNormal)
			const mP = m.transformPoint(p)
			const mSurface = surf.transform(m)
			ass.ok(mSurface.normalP(mP).like(mNormal))
			assert(mSurface.normalP(mP).like(mNormal))

			ass.ok(mSurface.containsPoint(mP))
			assert(mSurface.containsPoint(mP))


			//const mPSFlipped = mSurface.flipped()
			//ass.ok(mPSFlipped.normalP(mP).negated().like(mNormal))
			//assert(mPSFlipped.normalP(mP).negated().like(mNormal))
		}
	}
}

export function testImplicitSurface(t: Assert, surface: ImplicitSurface) {
    const EPS = 1e-8
    const testPoints = [
        V3.O.plus(V(0.2, 0, 0)), // V3.O fails on ellipsoidSurface
        V3.Y,
        V3.X,
        V3.Z.plus(V(0.2, 0, 0)),
        V3.XY,
        V3.XYZ,
        new V3(10, 10, 10),
        new V3(5, 6, 7)
    ]
    for (const testPoint of testPoints) {
        const i = surface.implicitFunction()(testPoint)
        const didpGuess = testPoint.map((el, dim) => {
            const i2 = surface.implicitFunction()(testPoint.plus(V3.O.withElement(dim, EPS)))
            return (i2 - i) / EPS
        })
        const didp = surface.didp(testPoint)
        t.push(didpGuess.to(didp).length() < 1e-6, didp, didpGuess, `actual: ${didp} guess: ${didpGuess} p: ${testPoint}`)
    }
}

export function testCurveISInfos(assert: Assert, c1: Curve, c2: Curve, count, f = 'isInfosWithCurve') {
	const intersections = c1[f](c2).map(info => info.p)
	linkB3(assert, { edges: [c1, c2].map(c => Edge.forCurveAndTs(c)), points: intersections }, `view`)
	assert.equal(intersections.length, count, `intersections.length == count: ${intersections.length} == ${count}`)
	intersections.forEach((is, i) => {
		assert.ok(intersections.every((is2, j) => j == i || !is.like(is2)), is.sce + ' is not unique ' + intersections)
		assert.ok(c1.containsPoint(is), `e1.containsPoint(is): ${c1.toSource()}.containsPoint(${is.sce},`)
		assert.ok(c2.containsPoint(is), `e2.containsPoint(is): ${c1.toSource()}.containsPoint(${is.sce},`)
	})
}

export function testISTs(assert: Assert, curve: Curve, surface: Surface | P3, tCount: int) {
	surface instanceof P3 && (surface = new PlaneSurface(surface))
	const ists = curve instanceof L3 ? surface.isTsForLine(curve) : curve.isTsWithSurface(surface)
	const points = ists.map(t => curve.at(t))
	linkBRep(assert, `mesh=[${surface}.toMesh()];edges=[${Edge.forCurveAndTs(curve, curve.tMin, curve.tMax)}];points=${points.sce}`,
		ists.join(', ') || 'view')
	assert.equal(ists.length, tCount, 'number of isps = ' + tCount)
	for (const t of ists) {
		const p = curve.at(t)
		assert.ok(surface.containsPoint(p), 'surface.containsPoint(p) ' + surface.str + ' ' + p.str
			+ ' t: ' + t
			+ (ImplicitSurface.is(surface) ? ' dist: ' + surface.implicitFunction()(p) : ''))
	}
}

export function linkBRep(assert: Assert, hash: string, message = 'view') {
	const escapedHash = encodeURIComponent(hash.replace(/, /g, ',').replace(/(\n|\t)+/g, '')).replace(/\(/g, '%28').replace(/\)/g, '%29')
	assert.link('http://localhost:10001/viewer.html#' + escapedHash, message)
}

export function testLoopContainsPoint(assert: Assert, surface: Surface, loop: Edge[], p: V3, result: PointVsFace) {
	const ccwLoop = surface.edgeLoopCCW(loop) ? loop : Edge.reversePath(loop)
	linkBRep(assert, `mesh=[${Face.create(surface, loop).sce}.toMesh()];points=[${p.sce}]`)
	assert.equal(surface.loopContainsPoint(loop, p), result)
}