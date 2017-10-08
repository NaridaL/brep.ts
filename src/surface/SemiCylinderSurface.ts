import {int, Tuple3, Tuple2, Tuple4,AABB,DEG,EllipticE,EllipticF,GOLDEN_RATIO,M4,MINUS,Matrix,NLA_DEBUG,NLA_PRECISION,P3XY,P3YZ,P3ZX,SCE,STR,TAU,Transformable,V,V3,Vector,addOwnProperties,arithmeticGeometricMean,arrayCopy,arrayCopyBlocks,arrayCopyStep,arrayFromFunction,arrayRange,arraySwap,assert,assertInst,assertNever,assertNumbers,assertVectors,assertf,between,bisect,callsce,canonAngle,ceil10,checkDerivate,clamp,combinations,defaultRoundFunction,disableConsole,doubleSignedArea,enableConsole,eq,eq0,eq02,eq2,eqAngle,floatHashCode,floor10,forceFinite,fuzzyBetween,fuzzyUniques,fuzzyUniquesF,gaussLegendre24Weights,gaussLegendre24Xs,gaussLegendreQuadrature24,ge,getIntervals,getRoots,glq24_11,glqInSteps,gt,hasConstructor,isCCW,le,lerp,lt,mapPush,midpointRuleQuadrature,mod,newtonIterate,newtonIterate1d,newtonIterate2d,newtonIterate2dWithDerivatives,newtonIterateSmart,newtonIterateWithDerivative,numberToStr,pqFormula,rad2deg,randomColor,repeatString,round10,snap,snap0,snap2,snapEPS,solveCubicReal2,time,zeroAngle} from 'ts3dutils'
import {Mesh, pushQuad} from 'tsgl'

import {Curve, P3, Surface, ProjectedCurveSurface, L3, ISInfo, ParametricSurface, ImplicitSurface,
    Edge,
    HyperbolaCurve,
    SemiEllipseCurve,
    ParabolaCurve,
    EllipseCurve, PointVsFace, PICurve} from '../index'

const {PI, cos, sin, min, max, tan, sign, ceil, floor, abs, sqrt, pow, atan2, round} = Math


export class SemiCylinderSurface extends ProjectedCurveSurface {
	readonly matrix: M4
	readonly inverseMatrix: M4
	readonly normalMatrix: M4
	readonly normalDir: number
	readonly baseCurve: SemiEllipseCurve

	constructor(baseCurve: SemiEllipseCurve,
	            dir1: V3,
	            sMin: number,
	            sMax: number,
	            zMin = -Infinity,
	            zMax = Infinity) {
		super(baseCurve, dir1, sMin, sMax, zMin, zMax)
		assertInst(SemiEllipseCurve, baseCurve)
		//assert(!baseCurve.normal1.isPerpendicularTo(dir1), !baseCurve.normal1.isPerpendicularTo(dir1))
		this.matrix = M4.forSys(baseCurve.f1, baseCurve.f2, dir1, baseCurve.center)
		this.inverseMatrix = this.matrix.inversed()
		this.normalDir = sign(this.baseCurve.normal.dot(this.dir))
		this.normalMatrix = this.matrix.as3x3().inversed().transposed().scale(this.normalDir)
		if (!sema) {
			sema = true
			assert(this)
			const src = this.toSource()
			const n = eval(src)
			sema = false
		}
	}

	getConstructorParameters(): any[] {
		return [this.baseCurve, this.dir, this.sMin, this.sMax, this.tMin, this.tMax]
	}

	normalP(p: V3): V3 {
		return this.normalMatrix.transformVector(this.inverseMatrix.transformPoint(p).xy()).unit()
	}

	loopContainsPoint(loop: Edge[], p: V3): PointVsFace {
		assertVectors(p)
        if (!this.containsPoint(p)) return OUTSIDE
		// create plane that goes through cylinder seam
		const line = new L3(p, this.dir.unit())
		const seamBase = this.baseCurve.at(PI)
		const lineOut = this.dir.cross(this.normalP(p))
		return Surface.loopContainsPointGeneral(loop, p, line, lineOut)
	}

	isTsForLine(line: L3) {
		assertInst(L3, line)
		// transforming line manually has advantage that dir1 will not be renormalized,
		// meaning that calculated values t for localLine are directly transferable to line
		const dirLC = this.inverseMatrix.transformVector(line.dir1)
		if (dirLC.isParallelTo(V3.Z)) {
			// line is parallel to this.dir
			return []
		}
		const anchorLC = this.inverseMatrix.transformPoint(line.anchor)
		assert(!SemiCylinderSurface.unitISLineTs(anchorLC, dirLC).length || !isNaN(SemiCylinderSurface.unitISLineTs(anchorLC, dirLC)[0]), 'sad ' + dirLC)
		return SemiCylinderSurface.unitISLineTs(anchorLC, dirLC)
	}

	isCoplanarTo(surface: Surface): surface is SemiCylinderSurface {
		return this == surface ||
			hasConstructor(surface, SemiCylinderSurface)
			&& this.dir.isParallelTo(surface.dir)
			&& this.containsSemiEllipse(surface.baseCurve, false)
	}

	like(surface: Surface): boolean {
		if (!this.isCoplanarTo(surface)) return false
		// normals need to point in the same direction (outwards or inwards) for both
		const thisFacesOut = 0 < this.baseCurve.normal.dot(this.dir)
		const objectFacesOut = 0 < surface.baseCurve.normal.dot(surface.dir)
		return thisFacesOut == objectFacesOut
	}

	containsSemiEllipse(ellipse: SemiEllipseCurve, checkAABB: boolean = true) {
		const projEllipse = ellipse.transform(M4.project(this.baseCurve.getPlane(), this.dir))
		return this.baseCurve == ellipse || this.baseCurve.isColinearTo(projEllipse) &&
			(!checkAABB || le(0, ellipse.transform(this.inverseMatrix).getAABB().min.y))
	}

	containsCurve(curve: Curve) {
		if (curve instanceof L3) {
			return this.containsLine(curve)
		} else if (curve instanceof SemiEllipseCurve) {
			return this.containsSemiEllipse(curve)
		} else if (curve instanceof BezierCurve) {
			return false
		} else {
			return super.containsCurve(curve)
		}
	}

	implicitFunction() {
		return (pWC) => {
			const pLC = this.inverseMatrix.transformPoint(pWC)
			const radiusLC = pLC.lengthXY()
			const normalDir = Math.sign(this.baseCurve.normal.dot(this.dir))
			return normalDir * (1 - radiusLC)
		}
	}

	containsPoint(pWC: V3): boolean {
	    const pLC = this.inverseMatrix.transformPoint(pWC)
		return SemiEllipseCurve.XYLCValid(pLC)
	}

	stP(pWC: V3): V3 {
        assert(arguments.length == 1)
        const pLC = this.inverseMatrix.transformPoint(pWC)
		const u = SemiEllipseCurve.XYLCPointT(pLC)
		return new V3(u, pLC.z, 0)
	}

	isCurvesWithSurface(surface2: Surface): Curve[] {
		if (surface2 instanceof PlaneSurface) {
			return this.isCurvesWithPlane(surface2.plane)
		} else if (surface2 instanceof SemiCylinderSurface) {
			if (surface2.dir.isParallelTo(this.dir)) {
				const projEllipse = surface2.baseCurve.transform(M4.project(this.baseCurve.getPlane(), this.dir))
				return this.baseCurve.isInfosWithEllipse(projEllipse).map(info => {
					const lineDir = sign(this.normalP(info.p).cross(surface2.normalP(info.p)).dot(this.dir)) || 1
					return new L3(info.p, this.dir.times(lineDir))
				})
			} else if (eq0(this.getCenterLine().distanceToLine(surface2.getCenterLine()))) {
				assert(false)
			} else {
				assert(false)
			}
		}
	}

	getCenterLine(): L3 {
		return new L3(this.baseCurve.center, this.dir)
	}

	static semicylinder(radius: number): SemiCylinderSurface {
		return new SemiCylinderSurface(new SemiEllipseCurve(V3.O, new V3(radius, 0, 0), new V3(0, radius, 0)), V3.Z, undefined, undefined)
	}

	/**
	 *
	 * @param anchorLC
	 * @param dirLC not necessarily unit
	 */
	static unitISLineTs(anchorLC: V3, dirLC: V3): number[] {
		const {x: ax, y: ay} = anchorLC
		const {x: dx, y: dy} = dirLC

		// this cylinder: x² + y² = 1
		// line: p = anchorLC + t * dirLC
		// split line equation into 3 component equations, insert into cylinder equation
		// x = ax + t * dx
		// y = ay + t * dy
		// (ax² + 2 ax t dx + t²dx²) + (ay² + 2 ay t dy + t²dy²) = 1
		// transform to form (a t² + b t + c = 0) and solve with pqFormula
		const a = dx ** 2 + dy ** 2
		const b = 2 * (ax * dx + ay * dy)
		const c = ax ** 2 + ay ** 2 - 1
		return pqFormula(b / a, c / a).filter(t => SemiEllipseCurve.XYLCValid(new V3(ax + dx * t, ay + dy * t, 0)))
	}

	facesOutwards(): boolean {
		return this.baseCurve.normal.dot(this.dir) > 0
	}

	getSeamPlane(): P3 {
		let normal = this.baseCurve.f1.cross(this.dir)
		normal = normal.times(-sign(normal.dot(this.baseCurve.f2)))
		return P3.normalOnAnchor(normal, this.baseCurve.center)
	}

	clipCurves(curves: Curve[]): Curve[] {
		return curves.flatMap(curve => curve.clipPlane(this.getSeamPlane()))
	}

	static readonly UNIT = new SemiCylinderSurface(SemiEllipseCurve.UNIT, V3.Z, undefined, undefined, 0, 1)
}
SemiCylinderSurface.prototype.uStep = TAU  / 32
SemiCylinderSurface.prototype.vStep = 256










