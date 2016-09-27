"use strict"

/**
 *
 * @param anchor2
 * @param right
 * @param up
 * @param upStart
 * @param upEnd
 * @param rightStart
 * @param rightEnd
 * @param color
 * @param name
 * @returns {P3}
 * @constructor
 */
function CustomPlane(anchor2, right, up, upStart, upEnd, rightStart, rightEnd, color, name) {
	var p = P3.forAnchorAndPlaneVectors(anchor2, right, up, CustomPlane.prototype)
	p.up = up
	p.right = right
	p.upStart = upStart
	p.upEnd = upEnd
	p.rightStart = rightStart
	p.rightEnd = rightEnd
	p.color = color
	p.id = globalId++
	p.name = name
	return p
}
CustomPlane.prototype = Object.create(P3.prototype);
CustomPlane.prototype.constructor = CustomPlane;
Object.defineProperty(CustomPlane.prototype, "plane", { get: function () { return this } });
CustomPlane.prototype.toString = function() {
	return "Plane #" + this.id;
}
CustomPlane.prototype.what ="Plane"
CustomPlane.prototype.distanceTo = function (line) {
	return [
		L3(this.anchor.plus(this.right.times(this.rightStart)), this.up),
		L3(this.anchor.plus(this.right.times(this.rightEnd)), this.up),
		L3(this.anchor.plus(this.up.times(this.upStart)), this.right),
		L3(this.anchor.plus(this.up.times(this.upEnd)), this.right)].map(function (line2, line2Index) {
		var info = line2.infoClosestToLine(line);
		if ((isNaN(info.t) // parallel lines
			|| line2Index < 2 && this.upStart <= info.t && info.t <= this.upEnd
			|| line2Index >= 2 && this.rightStart <= info.t && info.t <= this.rightEnd)
			&& info.distance <= 16) {
			return info.s;
		} else {
			return Infinity;
		}
	}, this).min();
}
CustomPlane.forPlane = function (plane, color, name) {
	var p = P3(plane.normal, plane.w, CustomPlane.prototype)
	p.up = plane.normal.getPerpendicular().normalized()
	p.right = p.up.cross(p.normal)
	p.upStart = -500
	p.upEnd = 500
	p.rightStart = -500
	p.rightEnd = 500
	p.color = color || NLA.randomColor()
	p.id = globalId++
	p.name = name
	return p;
}
