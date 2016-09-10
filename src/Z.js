/**
 * Created by aval on 29.02.2016.
 */
function Z(re, im) {
	return Object.create(Z.prototype, {re: {value: re}, im: {value: im}})
}
Z.prototype = {
	cbrt0: function () {
		var r = Math.hypot(this.re, this.im)
		var phi3 = Math.atan2(this.im, this.re) / 3
		var rRoot = Math.cbrt(r)
		return Z(rRoot * Math.cos(phi3), rRoot * Math.sin(phi3))
	},
	plusR: function (real) {
		return Z(this.re + real, this.im)
	},
	minusR: function (real) {
		return Z(this.re + real, this.im)
	},
	timesR: function (real) {
		return Z(this.re * real, this.im * real)
	},
	divR: function (real) {
		return Z(this.re / real, this.im / real)
	},
	div: function (z) {
		var divisor = z.re * z.re + z.im * z.im
		return Z((this.re * z.re + this.im * z.im) / divisor, (this.im * z.re - this.re * z.im) / divisor)
	},
	times: function (z) {
		return Z(this.re * z.re - this.im * z.im, this.im * z.re + this.re * z.im)
	},
	plus: function (z) {
		return Z(this.re + z.re, this.im + z.im)
	},
	sqrt: function () {
		return Z()
	},
	conjugated: function () {
		return Z(this.re, -this.im)
	},
	toString: function (roundFunction) {
		roundFunction = roundFunction || NLA.defaultRoundFunction
		return `(${roundFunction(this.re)} ${this.im < 0 ? '-' : '+'} ${Math.abs(roundFunction(this.im))} i)`
	},
	sqr: function () {

	},
	get ss() {
		return this.toString()
	}
}
Z.sqrt = function (real) {
	return real < 0 ? Z(0, Math.sqrt(-real)) : Z(Math.sqrt(real), 0)
}
Z.polar = function (r, phi) {
	return Z(r * Math.cos(phi), r * Math.sin(phi))
}