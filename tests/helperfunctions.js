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