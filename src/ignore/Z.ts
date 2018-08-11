class Complex {
	constructor(readonly re: number, readonly im: number) {}

	plus(zOrR: Complex | number): Complex {
		const z = 'number' === typeof zOrR ? Z(zOrR, 0) : zOrR
		return Z(this.re + z.re, this.im + z.im)
	}

	minus(zOrR: Complex | number): Complex {
		const z = 'number' === typeof zOrR ? Z(zOrR, 0) : zOrR
		return Z(this.re + z.re, this.im - z.re)
	}

	times(zOrR: Complex | number): Complex {
		return 'number' === typeof zOrR
			? Z(this.re * zOrR, this.im * zOrR)
			: Z(this.re * zOrR.re - this.im * zOrR.im, this.im * zOrR.re + this.re * zOrR.im)
	}

	div(zOrR: Complex | number): Complex {
		const z = 'number' === typeof zOrR ? Z(zOrR, 0) : zOrR
		const divisor = z.re * z.re + z.im * z.im
		return Z((this.re * z.re + this.im * z.im) / divisor, (this.im * z.re - this.re * z.im) / divisor)
	}

	conjugated() {
		return Z(this.re, -this.im)
	}

	negated(): Complex {
		return Z(-this.re, this.im)
	}

	toString(roundFunction) {
		roundFunction = roundFunction || defaultRoundFunction
		return `(${roundFunction(this.re)} ${this.im < 0 ? '-' : '+'} ${Math.abs(roundFunction(this.im))} i)`
	}

	sqr(): Complex {
		return Z(this.re * this.re - this.im * this.im, 2 * this.im * this.re)
	}

	sqrt() {
		const f = Math.SQRT1_2
		const r = Math.hypot(this.re, this.im)
		return Z(f * Math.sqrt(r + this.re), f * Math.sqrt(r - this.re) * Math.sign(this.im))
	}

	cbrt(): Complex {
		const r = Math.hypot(this.re, this.im)
		const phi3 = Math.atan2(this.im, this.re) / 3
		const rCbrt = Math.cbrt(r)
		return Z(rCbrt * Math.cos(phi3), rCbrt * Math.sin(phi3))
	}
}

function Z(re: number, im: number) {
	return new Complex(re, im)
}

namespace Z {
	export function sqrt(real: number): Complex {
		return real < 0 ? Z(0, Math.sqrt(-real)) : Z(Math.sqrt(real), 0)
	}

	export function polar(r: number, phi: number): Complex {
		return Z(r * Math.cos(phi), r * Math.sin(phi))
	}
}
