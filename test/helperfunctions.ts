QUnit.assert.B2equals = function(actual, expected, message) {
	if (!(actual instanceof B2)) {
		this.push(false, actual, null, "actual is not a B2")
		return
	}
	console.log(actual)

	this.equal(actual.faces.length, expected.faces.length, "no of faces")

	actual.faces.forEach(face => {
		if (!expected.faces.some(expectedFace => expectedFace.likeFace(face))) {
			this.ok(false, "Unexpected face in result:" + face.toSource())
		}
	})
}

QUnit.assert.fuzzyEquals = function(actual, expected, message) {
	this.push(NLA.eq(actual, expected), actual, expected, message)
}

function b2Equal(test, a, b, actual, expected) {

    test.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
						href='brep2.html?a=${a.toSource()}&b=${b.toSource()}&c=${expected.translate(20, 0, 0).toSource()}'>expected</a>`)
    test.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
						href='brep2.html?a=${a.toSource()}&b=${b.toSource()}&c=${actual.translate(20, 0, 0).toSource()}'>actual</a>`)
    test.B2equals(actual, expected)
}


QUnit.assert.V3ArraysLike = function (actual, expected, message) {
	this.push(expected.every((v, i) => v.like(actual[i])), actual.toSource(), expected.toSource(), message)
}


function registerTests(o: { [key: string]: (assert: Assert) => void }) {
	for (const key in o) {
		QUnit.test(key, o[key])
	}
}
