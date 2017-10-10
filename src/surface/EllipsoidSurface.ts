import {
	arrayFromFunction, assert, assertf, assertInst, assertNever, assertNumbers, assertVectors, between, eq, eq0,
	gaussLegendreQuadrature24, ge, glqInSteps, int, le, lt, M4, MINUS, NLA_PRECISION, pqFormula, V3,
} from 'ts3dutils'
import {Mesh} from 'tsgl'

import {
	Curve, CylinderSurface, Edge, EllipseCurve, ImplicitSurface, L3, P3, ParametricSurface, PlaneSurface, PointVsFace,
	ProjectedCurveSurface, Surface, dotCurve,
} from '../index'

const {PI, cos, sin, abs} = Math

export class EllipsoidSurface extends ParametricSurface implements ImplicitSurface {
	static readonly UNIT = new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z)
	readonly matrix: M4
	readonly inverseMatrix: M4
	readonly normalMatrix: M4
	readonly pLCNormalWCMatrix: M4
	readonly pWCNormalWCMatrix: M4
	readonly normalDir: number // -1 | 1

	constructor(readonly center: V3,
				readonly f1: V3,
				readonly f2: V3,
				readonly f3: V3) {
		super()
		assertVectors(center, f1, f2, f3)
		this.matrix = M4.forSys(f1, f2, f3, center)
		this.inverseMatrix = this.matrix.inversed()
		this.normalDir = sign(this.f1.cross(this.f2).dot(this.f3))
		this.pLCNormalWCMatrix = this.matrix.as3x3().inversed().transposed().scale(this.normalDir)
		this.pWCNormalWCMatrix = this.pLCNormalWCMatrix.times(this.inverseMatrix)
	}

	/**
	 * unit sphere: x² + y² + z² = 1
	 * line: p = anchor + t * dir |^2
	 * p² = (anchor + t * dir)^2
	 * 1 == (anchor + t * dir)^2
	 * 1 == anchor DOT anchor + 2 * anchor * t * dir + t² * dir DOT dir
	 */
	static unitISTsWithLine(anchor: V3, dir: V3): number[] {
		// for 0 = a t² + b t + c
		const a = dir.dot(dir)
		const b = 2 * anchor.dot(dir)
		const c = anchor.dot(anchor) - 1
		return pqFormula(b / a, c / a)
	}

	/**
	 * unit sphere: x² + y² + z² = 1
	 * plane: normal1 DOT p = w
	 */
	static unitISCurvesWithPlane(plane: P3): EllipseCurve[] {
		assertInst(P3, plane)
		let distPlaneCenter = Math.abs(plane.w)
		if (lt(distPlaneCenter, 1)) {
			// result is a circle
			// radius of circle: imagine right angled triangle (origin -> center of intersection circle -> point on
			// intersection circle) pythagoras: 1² == distPlaneCenter² + isCircleRadius² => isCircleRadius == sqrt(1 -
			// distPlaneCenter²)
			const isCircleRadius = Math.sqrt(1 - distPlaneCenter * distPlaneCenter)
			const center = plane.anchor
			const f1 = plane.normal1.getPerpendicular().toLength(isCircleRadius)
			const f2 = plane.normal1.cross(f1)
			return [new EllipseCurve(plane.anchor, f1, f2)]
		} else {
			return []
		}
	}

	static sphere(radius: number, center?: V3): EllipsoidSurface {
		assertNumbers(radius)
		center && assertVectors(center)
		return new EllipsoidSurface(center || V3.O, new V3(radius, 0, 0), new V3(0, radius, 0), new V3(0, 0, radius))
	}

	/**
	 * x²/a² + y²/b² + z²/c² = 1
	 */
	static forABC(a: number, b: number, c: number, center?: V3): EllipsoidSurface {
		return new EllipsoidSurface(center || V3.O, new V3(a, 0, 0), new V3(0, b, 0), new V3(0, 0, c))
	}

	static calculateAreaSpheroid(a: V3, b: V3, c: V3, edges: Edge[]): number {
		assertf(() => a.isPerpendicularTo(b))
		assertf(() => b.isPerpendicularTo(c))
		assertf(() => c.isPerpendicularTo(a))

		// handling discontinuities:
		// option 1: check for intersections with baseline, if there are any integrate parts separetely
		// "rotate" the edge so that there are no overlaps
		const matrix = M4.forSys(a, b, c), inverseMatrix = matrix.inversed()
		const circleRadius = a.length()
		const c1 = c.unit()
		const totalArea = edges.map(edge => {
			if (edge.curve instanceof EllipseCurve) {
				const f = (t) => {
					const at = edge.curve.at(t), tangent = edge.tangentAt(t)
					const localAt = inverseMatrix.transformPoint(at)
					const angleXY = localAt.angleXY()
					const arcLength = angleXY * circleRadius * Math.sqrt(1 + localAt.z ** 2)
					const scaling = Math.sqrt(1 + c1.dot(tangent) ** 2)
					return arcLength * scaling
				}
				const val = glqInSteps(f, edge.aT, edge.bT, 1)
				console.log('edge', edge, val)
				return val
			} else {
				assertNever()
			}
		}).sum()


		return totalArea
	}

	like(obj: any): boolean {
		return this.isCoplanarTo(obj) && this.isInsideOut() == obj.isInsideOut()
	}

	edgeLoopCCW(loop: Edge[]): boolean {
		throw new Error()
	}

	rootPoints() {

	}

	getConstructorParameters(): any[] {
		return [this.center, this.f1, this.f2, this.f3]
	}

	equals(obj: any): boolean {
		return this == obj ||
			Object.getPrototypeOf(obj) == this.constructor.prototype
			&& this.matrix.equals(obj.matrix)
	}

	isCurvesWithPlane(plane: P3): Curve[] {
		const planeLC = plane.transform(this.inverseMatrix)
		return EllipsoidSurface.unitISCurvesWithPlane(planeLC).map(c => c.transform(this.matrix))
	}

	isCurvesWithSurface(surface: Surface) {
		if (surface instanceof PlaneSurface) {
			return this.isCurvesWithPlane(surface.plane)
		} else if (surface instanceof CylinderSurface) {
			if (surface.dir.isParallelTo(this.dir1)) {
				const ellipseProjected = surface.baseCurve.transform(M4.project(this.baseEllipse.getPlane(), this.dir1))
				return this.baseEllipse.isInfosWithEllipse(ellipseProjected).map(info => new L3(info.p, this.dir1))
			} else if (eq0(this.getCenterLine().distanceToLine(surface.getCenterLine()))) {
				assert(false)
			} else {
				assert(false)
			}
		} else if (surface instanceof ProjectedCurveSurface) {
			const surfaceLC = surface.transform(this.inverseMatrix)
			const baseCurveLC = surfaceLC.baseCurve.project(new P3(surfaceLC.dir, 0))
			const ists = baseCurveLC.isTsWithSurface(EllipsoidSurface.UNIT)
			const insideIntervals = iii(ists, EllipsoidSurface.UNIT, baseCurveLC)
			const curves = insideIntervals.flatMap(ii => {
				const aLine = new L3(baseCurveLC.at(ii[0]), surfaceLC.dir)
				const a = EllipsoidSurface.UNIT.isTsForLine(aLine).map(t => aLine.at(t))
				const bLine = new L3(baseCurveLC.at(ii[1]), surfaceLC.dir)
				const b = EllipsoidSurface.UNIT.isTsForLine(bLine).map(t => bLine.at(t))
				return [0, 1].map(i => {
					let aP = a[i] || a[0], bP = b[i] || b[0]
					0 !== i && ([aP, bP] = [bP, aP])
					assert(EllipsoidSurface.UNIT.containsPoint(aP))
					assert(EllipsoidSurface.UNIT.containsPoint(bP))
					return PICurve.forStartEnd(surface, this.asEllipsoidSurface(), aP, bP)
				})
			})
			const f = (t) => baseCurveLC.at(t).length() - 1
			const fRoots = undefined
			return curves
		}
	}

	isTsForLine(line) {
		assertInst(L3, line)
		// transforming line manually has advantage that dir1 will not be renormalized,
		// meaning that calculated values t for localLine are directly transferable to line
		const anchorLC = this.inverseMatrix.transformPoint(line.anchor)
		const dirLC = this.inverseMatrix.transformVector(line.dir1)
		return EllipsoidSurface.unitISTsWithLine(anchorLC, dirLC)
	}

	isCoplanarTo(surface) {
		if (this === surface) return true
		if (surface.constructor !== EllipsoidSurface) return false
		if (!this.center.like(surface.center)) return false
		if (this.isSphere()) return surface.isSphere() && eq(this.f1.length(), this.f2.length())

		const localOtherMatrix = this.inverseMatrix.times(surface.matrix)
		// Ellipsoid with matrix localOtherMatrix is unit sphere iff localOtherMatrix is orthogonal
		return localOtherMatrix.is3x3() && localOtherMatrix.isOrthogonal()
	}

	containsEllipse(ellipse: EllipseCurve): boolean {
		const localEllipse = ellipse.transform(this.inverseMatrix)
		const distLocalEllipseCenter = localEllipse.center.length()
		const correctRadius = Math.sqrt(1 - distLocalEllipseCenter * distLocalEllipseCenter)
		return lt(distLocalEllipseCenter, 1) && localEllipse.isCircular() && localEllipse.f1.hasLength(correctRadius)
	}

	containsCurve(curve: Curve): boolean {
		if (curve instanceof EllipseCurve) {
			return this.containsEllipse(curve)
		} else {
			return super.containsCurve(curve)
		}
	}

	transform(m4: M4) {
		return new EllipsoidSurface(
			m4.transformPoint(this.center),
			m4.transformVector(this.f1),
			m4.transformVector(this.f2),
			m4.transformVector(this.f3))
	}

	isInsideOut(): boolean {
		return this.f1.cross(this.f2).dot(this.f3) < 0
	}

	flipped(): EllipsoidSurface {
		return new EllipsoidSurface(
			this.center,
			this.f1,
			this.f2,
			this.f3.negated())
	}

	toMesh(subdivisions: int = 3): Mesh {
		return Mesh.sphere(subdivisions).transform(this.matrix)
		// let mesh = new Mesh({triangles: true, lines: false, normals: true})
		// let pf = this.pSTFunc()
		// let pn = this.normalSTFunc()
		// let aCount = 32, bCount = 16, vTotal = aCount * bCount
		// for (let i = 0, a = -PI; i < aCount; i++, a += 2 * PI / aCount) {
		// 	for (let j = 0, b = -Math.PI / 2; j < bCount; j++, b += Math.PI / (bCount - 1)) {
		// 		mesh.vertices.push(pf(a, b))
		// 		mesh.normals.push(pn(a, b))
		// 		j != (bCount - 1) && pushQuad(mesh.triangles, true,
		// 			i * bCount + j, i * bCount + j + 1,
		// 			((i + 1) * bCount + j) % vTotal, ((i + 1) * bCount + j + 1) % vTotal)
		// 	}
		// }
		// mesh.compile()
		// return mesh
	}

	normalSTFunc() {
		// ugh
		// paramtric ellipsoid point q(a, b)
		// normal1 == (dq(a, b) / da) X (dq(a, b) / db) (Cross product of partial derivatives
		// normal1 == cos b * (f2 X f3 * cos b * cos a + f3 X f1 * cos b * sin a + f1 X f2 * sin b)
		return (a, b) => {
			let {f1, f2, f3} = this
			let normal = f2.cross(f3).times(Math.cos(b) * Math.cos(a))
				.plus(f3.cross(f1).times(Math.cos(b) * Math.sin(a)))
				.plus(f1.cross(f2).times(Math.sin(b)))
				//.times(Math.cos(b))
				.unit()
			return normal
		}
	}

	normalP(p) {
		return this.normalMatrix.transformVector(this.inverseMatrix.transformPoint(p)).unit()
	}

	normalST(s, t) {
		return this.normalMatrix.transformVector(V3.sphere(s, t))
	}

	pST(s: number, t: number): V3 {
		return this.matrix.transformPoint(V3.sphere(s, t))
	}

	//   d/dp (this.implicitFunction(p)) =
	// = d/dp (this.inverseMatrix.transformPoint(p).length() - 1)
	// = d/dp (this.inverseMatrix.transformPoint(p) * this.inverseMatrix.transformPoint(pWC).unit()

	dpds(): (s: number, t: number) => V3 {
		return (s: number, t: number) => this.matrix.transformVector(new V3(
			sin(s) * -cos(t),
			cos(s) * cos(t),
			0))
	}

	dpdt(): (s: number, t: number) => V3 {
		return (s: number, t: number) => this.matrix.transformVector(new V3(
			sin(t) * -cos(s),
			-sin(s) * sin(t),
			cos(t)))
	}

	stPFunc() {
		return (pWC: V3, hint?) => {
			const pLC = this.inverseMatrix.transformPoint(pWC)
			let alpha = pLC.angleXY()
			if (abs(alpha) > Math.PI - NLA_PRECISION) {
				assert(hint == -PI || hint == PI)
				alpha = hint
			}
			let beta = Math.asin(pLC.z)
			return new V3(alpha, beta, 0)
		}
	}

	isSphere(): boolean {
		return eq(this.f1.length(), this.f2.length())
			&& eq(this.f2.length(), this.f3.length())
			&& eq(this.f3.length(), this.f1.length())
			&& this.f1.isPerpendicularTo(this.f2)
			&& this.f2.isPerpendicularTo(this.f3)
			&& this.f3.isPerpendicularTo(this.f1)
	}

	isVerticalSpheroid(): boolean {
		return eq(this.f1.length(), this.f2.length())
			&& this.f1.isPerpendicularTo(this.f2)
			&& this.f2.isPerpendicularTo(this.f3)
			&& this.f3.isPerpendicularTo(this.f1)
	}

	implicitFunction() {
		return (pWC: V3) => {
			const pLC = this.inverseMatrix.transformPoint(pWC)
			return pLC.length() - 1
		}
	}

	// = this.inverseMatrix.transformPoint(this.inverseMatrix.transformPoint(pWC).unit())
	didp(pWC: V3) {
		const pLC = this.inverseMatrix.transformPoint(pWC)
		return this.inverseMatrix.transformVector(pLC.unit())
	}

	mainAxes(): EllipsoidSurface {
		// q(a, b) = f1 cos a cos b + f2 sin a cos b + f3 sin b
		// q(s, t, u) = s * f1 + t * f2 + u * f3 with s² + t² + u² = 1
		// (del q(a, b) / del a) = f1 (-sin a) cos b  + f2 cos a cos b
		// (del q(a, b) / del b) = f1 cos a (-sin b) + f2 sin a (-sin b) + f2 cos b
		// del q(s, t, u) / del a = -t f1 + s f2
		// (del q(a, b) / del a) DOT q(a, b) == 0
		// (f1 (-sin a) cos b  + f2 cos a cos b) DOT (f1 cos a cos b + f2 sin a cos b + f2 sin b) == 0
		// (del q(a, b) / del b) DOT q(a, b) == 0
		// (f1 cos a (-sin b) + f2 sin a (-sin b) + f2 cos b) DOT (f1 cos a cos b + f2 sin a cos b + f2 sin b) == 0

		// Solve[
		// (f1 (-sin a) cos b  + f2 cos a cos b) * (f1 cos a cos b + f2 sin a cos b + f2 sin b) = 0,
		// (f1 cos a (-sin b) + f2 sin a (-sin b) + f2 cos b) * (f1 cos a cos b + f2 sin a cos b + f2 sin b) = 0}, a, b]
		const {f1, f2, f3} = this

		if (eq0(f1.dot(f2)) && eq0(f2.dot(f3)) && eq0(f3.dot(f1))) {
			return this
		}

		//const f = ([a, b], x?) => {
		//    const sinA = Math.sin(a), cosA = Math.cos(a), sinB = Math.sin(b), cosB = Math.cos(b)
		//    const centerToP = V3.add(f1.times(cosA * cosB), f2.times(sinA * cosB), f3.times(sinB))
		//    const centerToPdelA = f1.times(-sinA * cosB).plus(f2.times(cosA * cosB))
		//    const centerToPdelB = V3.add(f1.times(cosA * -sinB), f2.times(sinA * -sinB), f3.times(cosB))
		//    x && console.log(centerToP.sce, centerToPdelA.sce, centerToPdelB.sce)
		//    return [centerToP.dot(centerToPdelA), centerToP.dot(centerToPdelB)]
		//}
		//const mainF1Params = newtonIterate(f, [0, 0], 8), mainF1 = this.pSTFunc()(mainF1Params[0], mainF1Params[1])
		//console.log(f(mainF1Params, 1).sce)
		//const mainF2Params = newtonIterate(f, this.stPFunc()(f2.rejectedFrom(mainF1)).toArray(2), 8),
		//   mainF2 = this.pSTFunc()(mainF2Params[0], mainF2Params[1])
		//console.log(this.normalSTFunc()(mainF2Params[0], mainF2Params[1]).sce)
		//assert(mainF1.isPerpendicularTo(mainF2), mainF1, mainF2, mainF1.dot(mainF2), mainF1Params)
		//const mainF3Params = this.stPFunc()(mainF1.cross(mainF2)), mainF3 = this.pSTFunc()(mainF3Params[0],
		// mainF3Params[1]) return new EllipsoidSurface(this.center, mainF1, mainF2, mainF3)

		const {U, SIGMA} = this.matrix.svd3()
		assert(SIGMA.isDiagonal())
		assert(U.isOrthogonal())
		const U_SIGMA = U.times(SIGMA)
		// column vectors of U_SIGMA
		const [mainF1, mainF2, mainF3] = arrayFromFunction(3, i => new V3(U_SIGMA.m[i], U_SIGMA.m[i + 4], U_SIGMA.m[i + 8]))
		return new EllipsoidSurface(this.center, mainF1, mainF2, mainF3)
	}

	containsPoint(p: V3): boolean {
		return eq0(this.implicitFunction()(p))
	}

	boundsFunction() {
		return (a, b) => between(b, -PI, PI)
	}

	volume(): number {
		return 4 / 3 * Math.PI * this.f1.dot(this.f2.cross(this.f3))
	}

	// TODO: also a commented out test
	//static splitOnPlaneLoop(loop: Edge[], ccw: boolean): [Edge[], Edge[]] {
	//const seamPlane = P3.ZX, seamSurface = new PlaneSurface(seamPlane)
	//const frontParts = [], backParts = [], iss = []
	//const colinearEdges = loop.map((edge) => seamSurface.containsCurve(edge.curve))
	//// a colinear edge is in front when
	//// ccw is true
	//// the edge curve is CCW on the seamPlane
	//// the edge is the same dir as the curve (bT > aT)
	//const colinearEdgesSide = loop.map((edge, i) => colinearEdges[i] &&
	//		(ccw ? 1 : -1) * seamPlane.normal1.dot(edge.curve.normal1) * (edge.bT - edge.aT))
	//
	//for (let edgeIndex = 0; edgeIndex < loop.length; edgeIndex++) {
	//	const edge = loop[edgeIndex]
	//	const nextEdgeIndex = (edgeIndex + 1) % loop.length, nextEdge = loop[nextEdgeIndex]
	//	//console.log(edge.toSource()) {p:V(2, -2.102, 0),
	//	if (colinearEdges[edgeIndex]) {
	//		const nextSide = colinearEdges[nextEdgeIndex] ? colinearEdgesSide[nextEdgeIndex]
	//			: dotCurve2(nextEdge.curve, nextEdge.aT, seamPlane.normal1, nextEdge.bT - nextEdge.aT)
	//		if (nextSide * colinearEdgesSide[edgeIndex] < 0) {
	//			iss.push({p: edge.b, t: 0, out: nextSide > 0})
	//		}
	//		(colinearEdgesSide[edgeIndex] > 0 ? frontParts : backParts).push(edge)
	//	} else {
	//		const f = sign(edge.bT - edge.aT)
	//		const ists = edge.edgeISTsWithPlane(seamPlane).sort((a, b) => f * (a - b))
	//		let prevT = edge.aT,
	//			prevP = edge.a,
	//			prevDir = edge.aDir,
	//			prevSide = snap0(seamPlane.distanceToPointSigned(edge.a)) || dotCurve2(edge.curve, edge.aT, V3.Y,
	// f) for (let i = 0; i < ists.length; i++) { const t = ists[i] if (edge.aT == t || edge.bT == t) { edge.bT ==
	// t && iss.push({p: edge.b, t: 0, out: true}) continue } const nextSide = dotCurve2(edge.curve, t, V3.Y, 1) if
	// (prevSide * nextSide < 0) { // switches sides, so: const newP = edge.curve.at(t) const newDir =
	// edge.tangentAt(t) const newEdge = Edge.create(edge.curve, prevP, newP, prevT, t, undefined, prevDir, newDir)
	// ;(prevSide > 0 ? frontParts : backParts).push(newEdge) iss.push({p: newP, t: 0, out: nextSide > 0}) prevP =
	// newP prevDir = newDir prevT = t prevSide = nextSide } } const lastEdge = Edge.create(edge.curve, prevP,
	// edge.b, prevT, edge.bT, undefined, prevDir, edge.bDir) ;(prevSide > 0 ? frontParts :
	// backParts).push(lastEdge) } } iss.forEach(is => is.t = V3.X.negated().angleRelativeNormal(is.p, V3.Y))
	// iss.sort((a, b) => a.t - b.t) let i = ccw == iss[0].out ? 1 : 0 const curve = new EllipseCurve(V3.O,
	// V3.X.negated(), V3.Z) //if (1 == i) {
	//    	//frontParts.push(
	//    	//	Edge.create(curve, V3.Y.negated(), iss[0].p, -PI, iss[0].t, undefined, V3.Z.negated(),
	// curve.tangentAt(iss[0].t)),
	////        Edge.create(curve, iss.last.p, V3.Y.negated(), iss.last.t, PI, undefined,
	// curve.tangentAt(iss.last.t), V3.Z.negated())) //} for (let i = ccw == iss[0].out ? 1 : 0; i < iss.length; i
	// += 2) {
	//    	let is0 = iss[i], is1 = iss[(i + 1) % iss.length]
	//	if (lt(is0.t, -PI) && lt(-PI, is1.t)) {
	//    		iss.splice(i + 1, 0, is1 = {p: V3.Y.negated(), t: -PI, out: true}, {p: V3.Y.negated(), t: -PI, out:
	// true})
	//	} else if (lt(is0.t, PI) && lt(PI, is1.t)) {
	//		iss.splice(i + 1, 0, is1 = {p: V3.Y, t: -PI, out: true}, {p: V3.Y, t: PI, out: true})
	//	}
	//	const edge = Edge.create(curve, is0.p, is1.p, is0.t, is1.t, undefined,
	//		curve.tangentAt(is0.t).times(sign(is1.t - is0.t)),
	//		curve.tangentAt(is1.t).times(sign(is1.t - is0.t)))
	//	frontParts.push(edge)
	//	backParts.push(edge.flipped())
	//}
	//return [frontParts, backParts]
	//}

	loopContainsPoint(loop: Edge[], p: V3): PointVsFace {
		assertVectors(p)
		const testLine = new EllipseCurve(
			this.center,
			this.matrix.transformVector(this.inverseMatrix.transformPoint(p).withElement('z', 0).unit()),
			this.f3)
		const pT = testLine.pointT(p)


		const lineOut = testLine.normal
		const testPlane = P3.normalOnAnchor(testLine.normal, p)
		const colinearEdges = loop.map((edge) => edge.curve.isColinearTo(testLine))
		let inside = false

		function logIS(isP) {
			const isT = testLine.pointT(isP)
			if (eq(pT, isT)) {
				return true
			} else if (pT < isT && le(isT, PI)) {
				inside = !inside
			}
		}

		for (let edgeIndex = 0; edgeIndex < loop.length; edgeIndex++) {
			const edge = loop[edgeIndex]
			const nextEdgeIndex = (edgeIndex + 1) % loop.length, nextEdge = loop[nextEdgeIndex]
			//console.log(edge.toSource()) {p:V(2, -2.102, 0),
			if (colinearEdges[edgeIndex]) {
				const lineAT = testLine.pointT(edge.a), lineBT = testLine.pointT(edge.b)
				if (le(Math.min(lineAT, lineBT), pT) && ge(pT, Math.max(lineAT, lineBT))) {
					return PointVsFace.ON_EDGE
				}
				// edge colinear to intersection
				const nextInside = colinearEdges[nextEdgeIndex] || dotCurve(lineOut, nextEdge.aDir, nextEdge.aDDT) < 0
				if (nextInside) {
					if (logIS(edge.b)) return PointVsFace.ON_EDGE
				}
			} else {
				for (const edgeT of edge.edgeISTsWithPlane(testPlane)) {
					if (edgeT == edge.bT) {
						if (!testLine.containsPoint(edge.b)) continue
						// endpoint lies on intersection testLine
						const edgeInside = dotCurve(lineOut, edge.bDir, edge.bDDT) < 0
						const nextInside = colinearEdges[nextEdgeIndex] || dotCurve(lineOut, nextEdge.aDir, nextEdge.aDDT) < 0
						if (edgeInside != nextInside) {
							if (logIS(edge.b)) return PointVsFace.ON_EDGE
						}
					} else if (edgeT != edge.aT) {
						const p = edge.curve.at(edgeT)
						if (!testLine.containsPoint(p)) continue
						// edge crosses testLine, neither starts nor ends on it
						if (logIS(p)) return PointVsFace.ON_EDGE
						// TODO: tangents?
					}
				}
			}
		}
		return inside ? PointVsFace.INSIDE : PointVsFace.OUTSIDE

	}

	surfaceAreaApprox(): number {
		// See https://en.wikipedia.org/wiki/Ellipsoid#Surface_area
		const mainAxes = this.mainAxes(),
			a = mainAxes.f1.length(),
			b = mainAxes.f2.length(),
			c = mainAxes.f3.length()
		const p = 1.6075
		return 4 * PI * Math.pow((Math.pow(a * b, p) + Math.pow(b * c, p) + Math.pow(c * a, p)) / 3, 1 / p)
	}

	surfaceArea(): number {
		// See https://en.wikipedia.org/wiki/Ellipsoid#Surface_area
		const mainAxes = this.mainAxes(),
			f1l = mainAxes.f1.length(),
			f2l = mainAxes.f2.length(),
			f3l = mainAxes.f3.length(),
			[c, b, a] = [f1l, f2l, f3l].sort(MINUS)

		// https://en.wikipedia.org/w/index.php?title=Spheroid&oldid=761246800#Area
		function spheroidArea(a, c) {
			if (c < a) {
				const eccentricity2 = 1 - c ** 2 / a ** 2
				const eccentricity = Math.sqrt(eccentricity2)
				return 2 * PI * a ** 2 * (1 + (1 - eccentricity2) / Math.sqrt(eccentricity) * Math.atanh(eccentricity))
			} else {
				const eccentricity = Math.sqrt(1 - a ** 2 / c ** 2)
				return 2 * PI * a ** 2 * (1 + c / a / eccentricity * Math.asin(eccentricity))
			}
		}

		if (eq(a, b)) {
			return spheroidArea(a, c)
		} else if (eq(b, c)) {
			return spheroidArea(b, a)
		} else if (eq(c, a)) {
			return spheroidArea(c, b)
		}

		const phi = Math.acos(c / a)
		const k2 = a ** 2 * (b ** 2 - c ** 2) / (b ** 2 * (a ** 2 - c ** 2)), k = Math.sqrt(k2)
		const incompleteEllipticInt1 = gaussLegendreQuadrature24(phi => Math.pow(1 - k2 * Math.sin(phi) ** 2, -0.5), 0, phi)
		const incompleteEllipticInt2 = gaussLegendreQuadrature24(phi => Math.pow(1 - k2 * Math.sin(phi) ** 2, 0.5), 0, phi)
		return 2 * PI * c ** 2 + 2 * PI * a * b / Math.sin(phi) * (incompleteEllipticInt2 * Math.sin(phi) ** 2 + incompleteEllipticInt1 * Math.cos(phi) ** 2)
	}

	getSeamPlane(): P3 {
		return P3.forAnchorAndPlaneVectors(this.center, this.f1, this.f3)
	}
}

EllipsoidSurface.prototype.uStep = PI / 32
EllipsoidSurface.prototype.vStep = PI / 32
