function b2equals(assert: Assert, actual, expected, message = '') {
	if (!(actual instanceof B2)) {
		assert.push(false, actual, null, 'actual is not a B2')
		return
	}

	assert.equal(actual.faces.length, expected.faces.length, 'no of faces')

	actual.faces.forEach(face => {
		if (!expected.faces.some(expectedFace => expectedFace.likeFace(face))) {
			assert.ok(false, 'Unexpected face in result:' + face.toSource())
		}
	})
}

QUnit.assert.fuzzyEquals = function(actual, expected, message) {
	this.push(eq(actual, expected), actual, expected, message)
}

function b2Equal(test, a, b, actual, expected) {

	linkB2(test, `a=${a.toSource()};b=${b.toSource()};c=${expected.translate(20, 0, 0).toSource(false)}'`, 'expected')
	linkB2(test, `a=${a.toSource()};b=${b.toSource()};c=${actual.translate(20, 0, 0).toSource(false)}`, 'actual')
	b2equals(test, actual, expected)
}
function b2EqualAnd(test, a: B2, b: B2, expected: B2) {
	let actual
	try {
		actual = a.and(b)
	} finally {
		if (actual) {
			const abWidth = a.getAABB().addAABB(b.getAABB()).size().x
			linkB3(test, {
				a, b, c: actual.translate(abWidth + 1).toSource(false), d: expected.translate(2 * (abWidth + 1)).toSource(false)})
			b2equals(test, actual, expected)
		} else {
			linkB3(test, {a, b})
		}
	}
}


QUnit.assert.V3ArraysLike = function (actual, expected, message) {
	this.push(expected.every((v, i) => v.like(actual[i])), actual.toSource(), expected.toSource(), message)
}


function registerTests(o: { [key: string]: (assert: Assert) => void }) {
	for (const key in o) {
		QUnit.test(key, o[key])
	}
}
function linkB2(assert: Assert, link, msg = 'view') {
	//link = link.replace(/, /g, ',').replace(/[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?/g, (numberStr) => {
	//	const f = parseFloat(numberStr), rd = round10(f, -7)
	//	return eq(f, rd) ? rd : f
	//})
	assert.ok(true, `<html><a href='${link}'>${msg}</a>`)
}
function linkB3(assert: Assert, values, msg = 'view') {
	linkB2(assert, makeLink(values), msg)
}
function testISCurves(assert: Assert, surface1: Surface | P3, surface2: Surface | P3, curveCount: int) {
	surface1 instanceof P3 && (surface1 = new PlaneSurface(surface1))
	surface2 instanceof P3 && (surface2 = new PlaneSurface(surface2))
	let isCurves
	try {
		isCurves = surface1.isCurvesWithSurface(surface2)
	} finally {
		if (isCurves) {
			linkB2(assert, `./viewer.html#mesh=[${surface1}.toMesh(), ${surface2}.toMesh()];edges=${isCurves.map(c => Edge.forCurveAndTs(c)).sce}`)

			assert.equal(isCurves.length, curveCount, 'number of curves = ' +  curveCount)
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
				//!expectedTangent.likeZero() && assert.ok(expectedTangent.isParallelTo(dp), 'pN1.cross(pN2).isParallelTo(dp)')
				//!expectedTangent.likeZero() && assert.ok(expectedTangent.dot(dp) > 0, 'pN1.cross(pN2).dot(dp) > 0')
			}
		} else {
			linkB2(assert, `./viewer.html#mesh=[${surface1}.toMesh(), ${surface2}.toMesh()]`)
		}
	}
}
function testLoopCCW(assert: Assert, surface: ConicSurface, loop: Edge[]) {
	const points = [loop[0].a, loop[0].atAvgT()]
	linkB2(assert, `viewer.html#mesh=${surface.sce}.toMesh();edges=${loop.toSource()};points=${points.sce}`)
	assert.push(surface.edgeLoopCCW(loop))
	assert.push(!surface.edgeLoopCCW(Edge.reversePath(loop)))
}
function testZDirVolume(assert: Assert, face: Face) {
	linkB2(assert, `mesh=${face.sce}.toMesh()`)
	const actual = face.zDirVolume().volume, expected = face.toMesh().calcVolume().volume
	assert.push(eq2(actual, expected, 0.1), actual, expected, 'diff = ' + (actual - expected))
}
function testCurve(ass: Assert, curve: Curve) {
	const STEPS = 12
	arrayFromFunction(STEPS, i => {
		const t = lerp(curve.tMin, curve.tMax, i / (STEPS - 1))
		const p = curve.at(t)
		// check that pointT and containsPoint behave as expected
		ass.pushResult({
			result: eq(t, curve.pointT(p)),
			actual: curve.pointT(p),
			expected: t,
			message: 't eq pointT(at(t) for ' + t})
		ass.ok(curve.containsPoint(p), `containsPoint(at(t = ${t}) = ${p})`)
	})

	// test curve length
	if (curve.arcLength !== Curve.prototype.arcLength) {
		const expected = glqInSteps(t => curve.tangentAt(t).length(), curve.tMin, curve.tMax, 4)
		const actual = curve.arcLength(curve.tMin, curve.tMax)
		ass.pushResult({result: eq2(expected, actual, 1e-6), expected, actual, message: 'curve should have same length as the numericaly calculated value'})
	}
}

function testParametricSurface(ass: Assert, surf: ParametricSurface) {
	linkB2(ass, `viewer.html#mesh=[${surf}.toMesh()]`, 'view')
	const params = [V(0.25, 0.25), V(0.6, 0.25), V(0.25, 0.6), V(0.6, 0.7)]
		.map(pm => new V3(lerp(surf.sMin, surf.sMax, pm.x), lerp(surf.tMin, surf.tMax, pm.y), 0))
	const points = params.map(({x, y}) => surf.pST(x, y))
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
	const matrices = [M4.mirroring(P3.XY), M4.mirroring(P3.YZ), M4.mirroring(P3.ZX)]
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
function testCurveISInfos(assert: Assert, c1: Curve, c2: Curve, count, f = 'isInfosWithCurve') {
	const intersections = c1[f](c2).map(info => info.p)
	linkB3(assert, {edges: [c1, c2].map(c => Edge.forCurveAndTs(c)), points: intersections}, `view ${c1} ${c2}`)
	assert.equal(intersections.length, count, `intersections.length == count: ${intersections.length} == ${count}`)
	intersections.forEach((is, i) => {
		assert.ok(intersections.every((is2, j) => j == i || !is.like(is2)), is.sce + ' is not unique ' + intersections)
		assert.ok(c1.containsPoint(is), `e1.containsPoint(is): ${c1.toSource()}.containsPoint(${is.sce},`)
		assert.ok(c2.containsPoint(is), `e2.containsPoint(is): ${c1.toSource()}.containsPoint(${is.sce},`)
	})
}
function testISTs(assert: Assert, curve: Curve, surface: Surface | P3, tCount: int) {
	surface instanceof P3 && (surface = new PlaneSurface(surface))
	const ists = curve instanceof L3 ? surface.isTsForLine(curve) : curve.isTsWithSurface(surface)
	const points = ists.map(t => curve.at(t))
	linkB2(assert, `viewer.html#mesh=[${surface}.toMesh()];edges=[${Edge.forCurveAndTs(curve, curve.tMin, curve.tMax)}];points=${points.sce}`,
		ists.join(', ') || 'view')
	assert.equal(ists.length, tCount, 'number of isps = ' +  tCount)
	for (const t of ists) {
		const p = curve.at(t)
		assert.ok(surface.containsPoint(p), 'surface.containsPoint(p) ' + surface.str + ' ' + p.str
			+ ' t: ' + t
			+ (surface.implicitFunction ? ' dist: ' + surface.implicitFunction()(p) : ''))
	}
}

function testLoopContainsPoint(assert: Assert, surface: Surface, loop: Edge[], p: V3, result: PointVsFace) {
	const ccwLoop = surface.edgeLoopCCW(loop) ? loop : Edge.reversePath(loop)
	linkB2(assert, `mesh=[${Face.create(surface, loop).sce}.toMesh()];points=[${p.sce}]`)
	assert.equal(surface.loopContainsPoint(loop, p), result)
}