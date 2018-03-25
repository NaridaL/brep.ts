import {
	assert,
	assertInst,
	assertVectors,
	eq,
	eq0,
	getIntervals,
	M4,
	newtonIterate,
	pqFormula,
	TAU,
	V3,
} from 'ts3dutils'

import {
	Curve,
	Edge,
	HyperbolaCurve,
	ImplicitSurface,
	L3,
	P3,
	ParabolaCurve,
	ParametricSurface,
	PlaneSurface,
	PointVsFace,
	SemiEllipseCurve,
	Surface,
} from '../index'

import { abs, cos, PI, sign, sin, sqrt } from '../math'

export class ConicSurface extends ParametricSurface implements ImplicitSurface {
	pointFoot(pWC: V3, ss?: number, st?: number): V3 {
		if (undefined === ss || undefined === st) {
			// similar to stP
			const pLC = this.matrixInverse.transformPoint(pWC)
			const angle = pLC.angleXY()
			if (undefined === ss) {
				ss = angle < -PI / 2 ? angle + TAU : angle
			}
			if (undefined === st) {
				st = pLC.z + (pLC.lengthXY() - pLC.z) * sqrt(2) / 2
			}
		}
		const f = ([s, t]: number[]) => {
			const pSTToPWC = this.pST(s, t).to(pWC)
			return [this.dpds()(s, t).dot(pSTToPWC), this.dpdt()(s).dot(pSTToPWC)]
		}
		const { 0: x, 1: y } = newtonIterate(f, [ss, st])
		return new V3(x, y, 0)
	}
	/**
	 * Unit cone. x² + y² = z², 0 <= z
	 */
	static readonly UNIT = new ConicSurface(V3.O, V3.X, V3.Y, V3.Z)
	readonly matrix: M4
	readonly matrixInverse: M4
	readonly pLCNormalWCMatrix: M4
	readonly normalDir: number // -1 | 1

	/**
	 * returns new cone C = {apex + f1 * z * cos(d) + f2 * z * sin(d) + f3 * z | -PI <= d <= PI, 0 <= z}
	 * @param f1
	 * @param f2
	 * @param dir Direction in which the cone opens. The ellipse spanned by f1, f2 is contained at (apex + f1).
	 */
	constructor(
		readonly center: V3,
		readonly f1: V3,
		readonly f2: V3,
		readonly dir: V3,
		sMin: number = 0,
		sMax: number = PI,
		tMin: number = 0,
		tMax: number = 16,
	) {
		super(sMin, sMax, tMin, tMax)
		assertVectors(center, f1, f2, dir)
		assert(0 <= tMin)
		this.matrix = M4.forSys(f1, f2, dir, center)
		this.matrixInverse = this.matrix.inversed()
		this.normalDir = sign(this.f1.cross(this.f2).dot(this.dir))
		this.pLCNormalWCMatrix = this.matrix
			.as3x3()
			.inversed()
			.transposed()
			.scale(this.normalDir)
	}

	get apex() {
		return this.center
	}

	static atApexThroughEllipse(
		apex: V3,
		ellipse: SemiEllipseCurve,
		sMin?: number,
		sMax?: number,
		tMin?: number,
		tMax?: number,
	): ConicSurface {
		assertVectors(apex)
		assertInst(SemiEllipseCurve, ellipse)
		return new ConicSurface(apex, ellipse.f1, ellipse.f2, apex.to(ellipse.center), sMin, sMax, tMin, tMax)
	}

	static unitISLineTs(anchor: V3, dir: V3): number[] {
		const { x: ax, y: ay, z: az } = anchor
		const { x: dx, y: dy, z: dz } = dir

		// this cone: x² + y² = z²
		// line: p = anchor + t * dir1
		// split line equation into 3 component equations, insert into cone equation
		// transform to form (a t² + b t + c = 0) and solve with pqFormula
		const a = dx * dx + dy * dy - dz * dz
		const b = 2 * (ax * dx + ay * dy - az * dz)
		const c = ax * ax + ay * ay - az * az
		// cone only defined for 0 <= z, so filter invalid values
		return pqFormula(b / a, c / a).filter(t => 0 < az + t * dz)
	}

	// calculate intersection of plane ax + cz = d and cone x² + y² = z²
	static unitISPlane(a: number, c: number, d: number): Curve[] {
		if (eq0(c)) {
			// plane is "vertical", i.e. parallel to Y and Z axes
			assert(!eq0(a)) // normal would be zero, which is invalid
			// z² - y² = d²/a²
			if (eq0(d)) {
				// d = 0 => z² - y² = 0 => z² = y² => z = y
				// plane goes through origin/V3.O
				return [
					new L3(V3.O, new V3(0, -sqrt(2) / 2, -sqrt(2) / 2), undefined, 0),
					new L3(V3.O, new V3(0, -sqrt(2) / 2, sqrt(2) / 2), 0),
				]
			} else {
				// hyperbola
				const center = new V3(d / a, 0, 0)
				const f1 = new V3(0, 0, abs(d / a)) // abs, because we always want the hyperbola to be pointing up
				const f2 = new V3(0, d / a, 0)
				return [new HyperbolaCurve(center, f1, f2)]
			}
		} else {
			// c != 0
			const aa = a * a,
				cc = c * c
			if (eq0(d)) {
				// ax + cz = d => x = d - cz / a => x² = d² - 2cdz/a + c²z²/a²
				// x² + y² = z²
				// => d² - 2cdz/a + c²z²/a² + y² = z²

				if (eq(aa, cc)) {
					return [new L3(V3.O, new V3(c, 0, -a).unit())]
				} else if (aa < cc) {
					throw new Error('intersection is single point V3.O')
				} else if (aa > cc) {
					return [
						new L3(V3.O, new V3(c, sqrt(aa - cc), -a).unit()),
						new L3(V3.O, new V3(c, -sqrt(aa - cc), -a).unit()),
					]
				}
			} else {
				if (eq(aa, cc)) {
					// parabola
					const parabolaVertex = new V3(d / 2 / a, 0, d / 2 / c)
					const parabolaVertexTangentPoint = new V3(d / 2 / a, d / c, d / 2 / c)
					const p2 = new V3(0, 0, d / c)
					const f2 = p2.minus(parabolaVertex)
					return [
						new ParabolaCurve(
							parabolaVertex,
							parabolaVertexTangentPoint.minus(parabolaVertex),
							f2.z < 0 ? f2.negated() : f2,
						),
					]
				} else if (aa < cc) {
					// ellipse
					const center = new V3(-a * d / (cc - aa), 0, d * c / (cc - aa))
					if (center.z < 0) {
						return []
					}
					const p1 = new V3(d / (a - c), 0, -d / (a - c))
					const p2 = new V3(-a * d / (cc - aa), d / sqrt(cc - aa), d * c / (cc - aa))
					return [new SemiEllipseCurve(center, center.to(p1), center.to(p2), -PI, PI)]
				} else if (aa > cc) {
					// hyperbola
					const center = new V3(-a * d / (cc - aa), 0, d * c / (cc - aa))
					const p1 = new V3(d / (a - c), 0, -d / (a - c))
					const p2 = new V3(-a * d / (cc - aa), d / sqrt(aa - cc), d * c / (cc - aa))
					const f1 = center.to(p1)
					return [new HyperbolaCurve(center, f1.z > 0 ? f1 : f1.negated(), center.to(p2))]
				}
			}
		}
		throw new Error('???')
	}

	equals(obj: any): boolean {
		return (
			this == obj ||
			(Object.getPrototypeOf(this) == Object.getPrototypeOf(obj) &&
				this.center.equals(obj.center) &&
				this.f1.equals(obj.f1) &&
				this.f2.equals(obj.f2) &&
				this.dir.equals(obj.dir))
		)
	}

	like(object: any): boolean {
		if (!this.isCoplanarTo(object)) return false
		// normals need to point in the same direction (outwards or inwards) for both
		return this.normalDir == object.normalDir
	}

	getVectors() {
		return [
			{ anchor: this.center, dir1: this.dir },
			{ anchor: this.center.plus(this.dir), dir1: this.f1 },
			{ anchor: this.center.plus(this.dir), dir1: this.f2 },
		]
	}

	getSeamPlane(): P3 {
		return P3.forAnchorAndPlaneVectors(this.center, this.f1, this.dir)
	}

	loopContainsPoint(contour: Edge[], p: V3): PointVsFace {
		assertVectors(p)
		const line = this.center.like(p)
			? new L3(p, this.matrix.transformVector(new V3(0, 1, 1)).unit())
			: L3.throughPoints(p, this.apex)
		const lineOut = line.dir1.cross(this.dir)

		return Surface.loopContainsPointGeneral(contour, p, line, lineOut)
	}

	getConstructorParameters(): any[] {
		return [this.center, this.f1, this.f2, this.dir, this.sMin, this.sMax, this.tMin, this.tMax]
	}

	isTsForLine(line: L3): number[] {
		// transforming line manually has advantage that dir1 will not be renormalized,
		// meaning that calculated values t for lineLC are directly transferable to line
		const anchorLC = this.matrixInverse.transformPoint(line.anchor)
		const dirLC = this.matrixInverse.transformVector(line.dir1)
		return ConicSurface.unitISLineTs(anchorLC, dirLC)
	}

	/**
	 * Interestingly, two cones don't need to have parallel dirs to be coplanar.
	 */
	isCoplanarTo(surface: Surface): boolean {
		if (this === surface) return true
		if (!(surface instanceof ConicSurface) || !this.apex.like(surface.apex)) return false
		// at this point apexes are equal
		return this.containsEllipse(new SemiEllipseCurve(surface.center.plus(surface.dir), surface.f1, surface.f2))
	}

	containsEllipse(ellipse: SemiEllipseCurve): boolean {
		const ellipseLC = ellipse.transform(this.matrixInverse)
		if (ellipseLC.center.z < 0) {
			return false
		}
		const { f1, f2 } = ellipseLC.rightAngled()
		const p1 = ellipseLC.center.plus(f1),
			p2 = ellipseLC.center.plus(f2)
		// check if both endpoints are on the cone's surface
		// and that one main axis is perpendicular to the Z-axis
		return eq(p1.x ** 2 + p1.y ** 2, p1.z ** 2) && eq(p2.x ** 2 + p2.y ** 2, p2.z ** 2) && (eq0(f1.z) || eq0(f2.z))
	}

	containsLine(line: L3): boolean {
		const lineLC = line.transform(this.matrixInverse)
		const d = lineLC.dir1
		return lineLC.containsPoint(V3.O) && eq(d.x * d.x + d.y * d.y, d.z * d.z)
	}

	containsParabola(curve: ParabolaCurve): boolean {
		assertInst(ParabolaCurve, curve)
		const curveLC = curve.transform(this.matrixInverse)
		if (curveLC.center.z < 0 || curveLC.f2.z < 0) {
			return false
		}
		const { center, f1, f2 } = curveLC.rightAngled()
		// check if center is on the surface,
		// that tangent is perpendicular to the Z-axis
		// and that "y" axis is parallel to surface
		return (
			eq(center.x * center.x + center.y * center.y, center.z * center.z) &&
			eq0(f1.z) &&
			eq(f2.x * f2.x + f2.y * f2.y, f2.z * f2.z)
		)
	}

	containsHyperbola(curve: HyperbolaCurve): boolean {
		assertInst(HyperbolaCurve, curve)
		return true
		const curveLC = curve.transform(this.matrixInverse)
		if (curveLC.center.z < 0 || curveLC.f2.z < 0) {
			return false
		}
		const { center, f1, f2 } = curveLC.rightAngled()
		// check if center is on the surface,
		// that tangent is perpendicular to the Z-axis
		return true
		return eq(center.x * center.x + center.y * center.y, center.z * center.z) && eq0(f1.z)
	}

	containsCurve(curve: Curve): boolean {
		if (curve instanceof SemiEllipseCurve) {
			return this.containsEllipse(curve)
		} else if (curve instanceof L3) {
			return this.containsLine(curve)
		} else if (curve instanceof HyperbolaCurve) {
			return this.containsHyperbola(curve)
		} else if (curve instanceof ParabolaCurve) {
			return this.containsParabola(curve)
		} else {
			return super.containsCurve(curve)
		}
	}

	transform(m4: M4): this {
		return new ConicSurface(
			m4.transformPoint(this.center),
			m4.transformVector(this.f1).times(m4.isMirroring() ? -1 : 1),
			m4.transformVector(this.f2),
			m4.transformVector(this.dir),
		) as this
	}

	flipped(): this {
		return new ConicSurface(this.center, this.f1.negated(), this.f2, this.dir) as this
	}

	normalSTFunc(): (s: number, t: number) => V3 {
		const { f1, f2 } = this,
			f3 = this.dir
		return (d, z) => {
			return f2
				.cross(f1)
				.plus(f2.cross(f3.times(Math.cos(d))))
				.plus(f3.cross(f1.times(Math.sin(d))))
				.unit()
		}
	}

	normalP(p: V3): V3 {
		//TODO assert(!p.like(this.center))
		const pLC = this.matrixInverse.transformPoint(p)
		return this.normalSTFunc()(pLC.angleXY(), pLC.z)
	}

	pSTFunc(): (s: number, t: number) => V3 {
		return (s, t) => {
			// center + f1 t cos s + f2 t sin s + t dir
			const resultLC = new V3(t * cos(s), t * sin(s), t)
			return this.matrix.transformPoint(resultLC)
		}
	}

	dpds(): (s: number, t: number) => V3 {
		return (s, t) => {
			const resultLC = new V3(t * -sin(s), t * cos(s), 0)
			return this.matrix.transformVector(resultLC)
		}
	}

	dpdt(): (s: number) => V3 {
		return s => {
			const resultLC = new V3(cos(s), sin(s), 1)
			return this.matrix.transformVector(resultLC)
		}
	}

	implicitFunction(): (pWC: V3) => number {
		return pWC => {
			const pLC = this.matrixInverse.transformPoint(pWC)
			const radiusLC = pLC.lengthXY()
			return this.normalDir * (radiusLC - pLC.z)
		}
	}

	didp(pWC: V3): V3 {
		const pLC = this.matrixInverse.transformPoint(pWC)
		return this.pLCNormalWCMatrix.transformVector(
			pLC
				.xy()
				.unit()
				.withElement('z', -1)
				.times(this.normalDir),
		)
	}

	containsPoint(p: V3) {
		return eq0(this.implicitFunction()(p))
	}

	stP(pWC: V3) {
		const pLC = this.matrixInverse.transformPoint(pWC)
		const angle = pLC.angleXY()
		return new V3(angle < -PI / 2 ? angle + TAU : angle, pLC.z, 0)
	}

	isCurvesWithSurface(surface: Surface): Curve[] {
		if (surface instanceof PlaneSurface) {
			return this.isCurvesWithPlane(surface.plane)
		} else if (ImplicitSurface.is(surface)) {
			return ParametricSurface.isCurvesParametricImplicitSurface(
				this,
				surface,
				0.1,
				0.1 / this.dir.length(),
				0.02,
			)
		}
		return super.isCurvesWithSurface(surface)
	}

	getCenterLine() {
		return new L3(this.center, this.dir)
	}

	isCurvesWithPlane(plane: P3): Curve[] {
		assertInst(P3, plane)
		const planeLC = plane.transform(this.matrixInverse)
		const planeNormal = planeLC.normal1
		const c = planeNormal.z
		/** "rotate" plane normal1 when passing to {@link ConicSurface.unitISPlane} so that
		 *  y-component of normal1 is 0 */
		const a = planeNormal.lengthXY()
		const d = planeLC.w
		// generated curves need to be rotated back before transforming to world coordinates
		const rotationMatrix = M4.rotateZ(planeNormal.angleXY())
		const wcMatrix = eq0(planeNormal.lengthXY()) ? this.matrix : this.matrix.times(rotationMatrix)
		return ConicSurface.unitISPlane(a, c, d).flatMap<Curve>(curve => {
			const curveWC = curve.transform(wcMatrix)
			if (curve instanceof SemiEllipseCurve) {
				const curveLC = curve.transform(rotationMatrix)
				const ts = curveLC.isTsWithPlane(P3.ZX)
				const intervals = getIntervals(ts, -PI, PI).filter(([a, b]) => curveLC.at((a + b) / 2).y > 0)
				return intervals.flatMap(([a, b]) => (curveWC as SemiEllipseCurve).split(a, b))
			}
			const p = curveWC.at(0.2)
			return this.normalP(p)
				.cross(plane.normal1)
				.dot(curveWC.tangentAt(0.2)) > 0
				? curveWC
				: curveWC.reversed()
		})
	}
}

ConicSurface.prototype.uStep = PI / 16
ConicSurface.prototype.vStep = 256
