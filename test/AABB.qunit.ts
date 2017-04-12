QUnit.module('AABB')
QUnit.test('intersectLines', function (assert) {
	const aabb = new AABB(V3.O, V3.XYZ)
	assert.notOk(aabb.intersectsLine(new L3(V(2, 0, 0), V3.Z)))
	assert.ok(aabb.intersectsLine(new L3(V(1, 0, 0), V3.Z)))
	assert.ok(aabb.intersectsLine(new L3(V(0.5, 0, 0), V3.Z)))
	assert.ok(aabb.intersectsLine(new L3(V(0.5, 0.5, 0), V3.Z)))
	assert.ok(aabb.intersectsLine(new L3(V3.XYZ.negated(), V3.XYZ)))
})