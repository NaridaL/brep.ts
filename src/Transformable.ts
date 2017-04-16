abstract class Transformable extends NLA.BaseObject {
	mirrored(plane: P3): this {
		return this.transform(M4.mirroring(plane)) as this
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

	translate(x: number, y: number = 0, z: number = 0): this {
		return this.transform(M4.translation(x, y, z), `.translate(${x}, ${y}, ${z})`) as this
	}

	scale(f): this
	scale(x: number, y: number, z: number): this
	scale(x, y?, z?): this {
        return 1 == arguments.length
            ? this.transform(M4.scaling(x), `.scale(${x})`) as this
            : this.transform(M4.scaling(x, y, z), `.scale(${x}, ${y}, ${z})`) as this
	}

	rotateX(radians: raddd): this {
		return this.transform(M4.rotationX(radians), `.rotateX(${radians})`) as this
	}

	rotateY(radians: raddd): this {
		return this.transform(M4.rotationY(radians), `.rotateY(${radians})`) as this
	}

	rotateZ(radians: raddd): this {
		return this.transform(M4.rotationZ(radians), `.rotateZ(${radians})`) as this
	}

    rotate(rotationCenter: V3, rotationAxis: V3, radians: raddd): this {
        console.log(M4.rotationLine(rotationCenter, rotationAxis, radians).str)
        return this.transform(M4.rotationLine(rotationCenter, rotationAxis, radians),
            `.rotate(${rotationCenter.sce}, ${rotationAxis.sce}, ${radians})`) as this
    }

    rotateAB(from: V3, to: V3): this {
        return this.transform(M4.rotationAB(from, to), `.rotateAB(${from.sce}, ${to.sce})`) as this
    }

	eulerZXZ(alpha: raddd, beta: raddd, gamma: raddd): this {
		return this.transform(M4.eulerZXZ(alpha, beta, gamma)) as this
	}

	project(plane: P3): this {
		return this.transform(M4.projection(plane)) as this
	}

	projXY(): this {
		return this.transform(M4.projection(P3.XY)) as this
	}

	projYZ(): this {
		return this.transform(M4.projection(P3.YZ)) as this
	}

	projZX(): this {
		return this.transform(M4.projection(P3.ZX)) as this
	}

	shearedX(y: number, z: number): this {
		return this.transform(new M4([
			1, y, z, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1])) as this
	}

	foo(): this {
		return this.transform(M4.FOO) as this
	}

	bar(): this {
		return this.transform(M4.BAR) as this
	}

	abstract transform(m4: M4, desc?: string): Transformable

	visit(visitor: {[name: string]: (o: Transformable, ...args: any[]) => any}, ...args: any[]) {
		let proto = Object.getPrototypeOf(this)
		// walk up the prototype chain until we find a defined function in o
		while (!visitor.hasOwnProperty(proto.constructor.name) && proto !== Transformable.prototype) {
			proto = Object.getPrototypeOf(proto)
		}
		if (visitor.hasOwnProperty(proto.constructor.name)) {
			return visitor[proto.constructor.name](this, ...args)
		} else {
			throw new Error('No implementation for ' + this.constructor.name)
		}
	}
}
