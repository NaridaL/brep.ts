abstract class Transformable extends NLA.BaseObject {
	mirrored(plane): this {
		return this.transform(M4.mirroring(plane));
	}

	mirroredX(): this {
		return this.mirrored(P3.YZ)
	}

	mirroredY(): this {
		return this.mirrored(P3.ZX)
	}

	mirroredZ(): this {
		return this.mirrored(P3.XY)
	}

	translate(x, y, z): this {
		return this.transform(M4.translation(x, y, z), `.translate(${x}, ${y}, ${z})`)
	}

	scale(x, y, z): this {
		return this.transform(M4.scaling(x, y, z))
	}

	rotateX(radians): this {
		return this.transform(M4.rotationX(radians), `.rotateX(${radians})`)
	}

	rotateY(radians): this {
		return this.transform(M4.rotationY(radians), `.rotateY(${radians})`)
	}

	rotateZ(radians): this {
		return this.transform(M4.rotationZ(radians), `.rotateZ(${radians})`)
	}

	rotate(rotationCenter, rotationAxis, radians): this {
		console.log(M4.rotationLine(rotationCenter, rotationAxis, radians).str)
		return this.transform(M4.rotationLine(rotationCenter, rotationAxis, radians),
			`.rotate(${rotationCenter.sce}, ${rotationAxis.sce}, ${radians})`)
	}

	eulerZXZ(alpha, beta, gamma): this {
		return this.transform(M4.eulerZXZ(alpha, beta, gamma))
	}

	project(plane): this {
		return this.transform(M4.projection(plane))
	}

	projXY(): this {
		return this.transform(M4.projection(P3.XY))
	}

	projYZ(): this {
		return this.transform(M4.projection(P3.YZ))
	}

	projZX(): this {
		return this.transform(M4.projection(P3.ZX))
	}

	shearedX(y, z): this {
		return this.transform(new M4([
			1, y, z, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1]))
	}

	foo(): this {
		return this.transform(M4.FOO)
	}

	bar(): this {
		return this.transform(M4.BAR)
	}

	abstract transform(m4: M4, desc?: string): this

}
