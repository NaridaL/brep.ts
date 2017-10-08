class Complex {
	constructor(readonly re: number, readonly im: number) {}

	cbrt(): Complex {
		const r = Math.hypot(this.re, this.im)
		const phi3 = Math.atan2(this.im, this.re) / 3
		const rRoot = Math.cbrt(r)
		return Z(rRoot * Math.cos(phi3), rRoot * Math.sin(phi3))
	}

	plusR(real: number): Complex {
		return Z(this.re + real, this.im)
	}

	minusR(real: number): Complex {
		return Z(this.re + real, this.im)
	}

	timesR(real: number): Complex {
		return Z(this.re * real, this.im * real)
	}

	divR(real: number): Complex {
		return Z(this.re / real, this.im / real)
	}

	div(z: Complex): Complex {
		const divisor = z.re * z.re + z.im * z.im
		return Z((this.re * z.re + this.im * z.im) / divisor, (this.im * z.re - this.re * z.im) / divisor)
	}

	times(z: Complex | number): Complex {
		return 'number' == typeof z
			? Z(this.re * z, this.im)
			: Z(this.re * z.re - this.im * z.im, this.im * z.re + this.re * z.im)
	}

	plus(z: Complex | number): Complex {
		return Z(this.re + z.re, this.im + z.im)
	}

	sqrt() {
		return Z()
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
}


function Z(re: number, im: number) { return new Complex(re, im) }

namespace Z {
	export function sqrt(real: number): Complex {
		return real < 0 ? Z(0, Math.sqrt(-real)) : Z(Math.sqrt(real), 0)
	}

	export function polar(r: number, phi: number): Complex {
		return Z(r * Math.cos(phi), r * Math.sin(phi))
	}
}