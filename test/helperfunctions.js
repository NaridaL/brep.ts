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
	this.push(NLA.equals(actual, expected), actual, expected, message)
}

QUnit.assert.b2Equal = function (a, b, actual, expected) {

	this.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
						href='file:///C:/Users/aval/Desktop/cs/brep2.html?a=${a.toSource()}&b=${b.toSource()}&c=${expected.translate(20, 0, 0).toSource()}'>link</a>`)
	this.ok(true, `<html><a style='color: #0000ff; text-decoration: underline;' target='blank'
						href='file:///C:/Users/aval/Desktop/cs/brep2.html?a=${a.toSource()}&b=${b.toSource()}&c=${actual.translate(20, 0, 0).toSource()}'>link</a>`)
	this.B2equals(actual, expected)
}