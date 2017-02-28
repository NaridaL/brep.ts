type SketchSegment = SketchBezier | SketchArc | SketchLineSeg

class Sketch {
	planeRef: NameRef
	elements: (SketchLineSeg | SketchBezier | SketchArc)[]
	constraints: (Constraint|any)[]
	name: string
	plane: CustomPlane
	hide: boolean

	F: ((x: number[]) => number[]|number)[]
	x: number[]
	b: number[]
	varMap: Map<SegmentEndPoint, int>
	worldToSketchMatrix: M4
	sketchToWorldMatrix: M4

	constructor() {
		// elements in 2D coordinates on x-y plane
		this.planeRef = NameRef.UNASSIGNED
		this.elements = []
		this.constraints = []
		this.name = "sketch" + (globalId++)
		this.plane = null
	}


	/**
	 * Can be called even if already removed
	 */
	removeSegment(seg: SketchSegment) {
		seg.points.forEach(p => p.freeFromConstraints(this))
		this.elements.remove(seg)
		this.constraints.forEach(constraint => constraint.constrains(seg) && removeFromConstraint(seg, this, constraint))
		selected.remove(seg)
	}

	removeElement(el: SketchSegment | SegmentEndPoint) {
		this.removeSegment((el as any).line || el)
	}

	dependentOnNames():NameRef[] {
		return Array.prototype.concat.apply(
			[this.planeRef],
			this.constraints.map(constraint => constraint.cs.filter(c => c instanceof NameRef)))
	}

	getConstraintsFor(el) {
		return this.constraints.filter((constraint) => constraint.constrains(el))
	}


	constrainDistancePointFixedCurveWC(point: SegmentEndPoint, curve: Curve, distance: number): void {
		if (curve instanceof L3) {
			this.constrainDistancePointFixedLineWC(point, curve, distance)
		} else if (curve instanceof EllipseCurve) {
			this.constrainDistancePointFixedEllipseWC(point, curve, distance)
		} else if (curve instanceof BezierCurve) {
			this.constrainDistancePointFixedGeneralCurveWC(point, curve, distance)
		} else {
			throw new Error('Cannot constrain distance point to ' + curve.constructor.name)
		}
	}

	constrainDistancePointFixedEllipseWC(point: SegmentEndPoint, curveWC: EllipseCurve, dist: number) {
		let ellSC = curveWC.transform(this.worldToSketchMatrix).project(P3.XY)
		let ip = this.varMap.get(point)
		console.log("p on ell")
		if (ellSC.isCircular()) {
			console.log("p circ")
			if (NLA.eq0(dist)) {
				this.constrainDistancePointFixedPoint(point, ellSC.center.x, ellSC.center.y, ellSC.f1.length())
			} else {
				this.b.push(dist)
				this.F.push(x => Math.abs(distance(x[ip], x[ip + 1], ellSC.center.x, ellSC.center.y) - dist))
			}
		} else {
			this.constrainDistancePointFixedGeneralCurveWC(point, curveWC, dist)
		}
	}
	constrainDistancePointFixedGeneralCurveWC(point: SegmentEndPoint, curveWC: Curve, dist: number) {
		let curveSC = curveWC.transform(this.worldToSketchMatrix).project(P3.XY)
		let ip = this.varMap.get(point)
		console.log("p not curc")
		let tStart = curveSC.closestTToPoint(point.V3())
		let tIndex = this.x.length
		this.x.push(tStart)
		this.b.push(dist) // target for distance
		this.b.push(0) // target for dot product
		this.F.push(x => {
			let p = new V3(x[ip], x[ip + 1], 0)
			let t = x[tIndex]
			let tPoint = curveSC.at(t)
			let tTangent = curveSC.tangentAt(t)
			let tPointToP = p.minus(tPoint)
			return [tPointToP.length(), tPointToP.dot(tTangent)]
		})

	}

	constrainDistancePointFixedLineWC(point, lineWC, distance) {
		this.b.push(distance)
		var px = this.varMap.get(point), py = px + 1
		var lineA_SC = this.worldToSketchMatrix.transformPoint(lineWC.anchor)
		var lineB_SC = lineA_SC.plus(this.worldToSketchMatrix.transformVector(lineWC.dir1))
		// console.log(lineA_SC, lineB_SC);
		if (NLA.eq0(distance)) {
			this.F.push(x => distanceLinePointSigned(lineA_SC.x, lineA_SC.y, lineB_SC.x, lineB_SC.y, x[px], x[py]))
		} else {
			this.F.push(x => distanceLinePoint(lineA_SC.x, lineA_SC.y, lineB_SC.x, lineB_SC.y, x[px], x[py]))
		}

		// console.log("calling F", lineA_SC.x, lineA_SC.y, lineB_SC.x, lineB_SC.y, x[px], x[py]);

	}

	constrainAngleSegmentToFixedPlane(segment, plane, cosAngle) {
		var sketchLineWC = this.plane.intersectionWithPlane(plane)
		if (null == sketchLineWC) throw new Error("no intersection!!")
		this.constrainAngleSegmentFixedLineWC(segment, sketchLineWC, cosAngle)
	}

	constrainAngleSegmentFixedLineWC(segment, lineWC, cosAngle) {
		console.log("WARGH")
		// indexes:
		var ia = this.varMap.get(segment.a), ib = this.varMap.get(segment.b)
		var lineDirSC = this.worldToSketchMatrix.transformVector(lineWC.dir1)
		if (1 == cosAngle || -1 == cosAngle) {
			this.b.push(0)
			this.F.push(
				x => Math.sin(Math.atan2(x[ib + 1] - x[ia + 1], x[ib] - x[ia]) - Math.atan2(lineDirSC.y, lineDirSC.x)))
		} else {
			this.b.push(cosAngle)
			this.F.push(x => angleVectors(x[ib] - x[ia], x[ib + 1] - x[ia + 1], lineDirSC.x, lineDirSC.y))
		}
	}

	constrainEqualDistance(ia, ib, ic, id) {
		this.b.push(0)
		this.F.push((x) => distance(x[ia], x[ia + 1], x[ib], x[ib + 1]) - distance(x[ic], x[ic + 1], x[id], x[id + 1]))
	}

	constrainAngleSegmentSegment(line1, line2, cosAngle) {
		this.b.push(cosAngle * cosAngle)
		var ia = this.varMap.get(line1.a),
			ib = this.varMap.get(line1.b),
			ic = this.varMap.get(line2.a),
			id = this.varMap.get(line2.b)
		this.F.push((x) => {
			let angle = angleABCD(
				x[ia], x[ia + 1],
				x[ib], x[ib + 1],
				x[ic], x[ic + 1],
				x[id], x[id + 1]
			)
			return angle * angle
		})
	}

	constrainDistancePointPoint(pA: SegmentEndPoint, pB: SegmentEndPoint, pDistance: number) {
		this.b.push(pDistance)
		var ia = this.varMap.get(pA), ib = this.varMap.get(pB)
		this.F.push((x) => distance(x[ia], x[ia + 1], x[ib], x[ib + 1]))
	}

	constrainDistancePointFixedPoint(pA: SegmentEndPoint, bx: number, by: number, pDistance: number) {
		this.b.push(pDistance)
		var ia = this.varMap.get(pA)
		this.F.push((x) => distance(x[ia], x[ia + 1], bx, by))
	}

	constrainDistancePointSegment(point, segment, distance) {
		this.b.push(distance)
		var ip = this.varMap.get(point), ia = this.varMap.get(segment.a), ib = this.varMap.get(segment.b)
		if (NLA.eq0(distance)) {
			this.F.push(x => distanceLinePointSigned(x[ia], x[ia + 1], x[ib], x[ib + 1], x[ip], x[ip + 1]))
		} else {
			this.F.push(x => distanceLinePoint(x[ia], x[ia + 1], x[ib], x[ib + 1], x[ip], x[ip + 1]))
		}
	}

	constrainDistancePointBezier(point, bezier, pDistance) {
		let startT = bezier.getCurve().closestTToPoint(point.V3())
		let startTIndex = this.x.length
		this.x.push(startT)
		this.b.push(pDistance)
		this.b.push(0)
		var ip = this.varMap.get(point),
			ib = bezier.points.map(p => this.varMap.get(p))
		this.F.push(
			x => {
				let {x:tPx, y:tPy} = bezierCurveAt(
					x[ib[0]], x[ib[0] + 1],
					x[ib[2]], x[ib[2] + 1],
					x[ib[3]], x[ib[3] + 1],
					x[ib[1]], x[ib[1] + 1],
					x[startTIndex]
				)
				let {x:tPdx, y:tPdy} = bezierCurveTangentAt(
					x[ib[0]], x[ib[0] + 1],
					x[ib[2]], x[ib[2] + 1],
					x[ib[3]], x[ib[3] + 1],
					x[ib[1]], x[ib[1] + 1],
					x[startTIndex]
				)
				return [
					distance(tPx, tPy, x[ip], x[ip + 1]),
					tPdx * (tPx - x[ip]) + tPdy * (tPy - x[ip + 1])]
			})
	}

	constrainAngleSegmentSegment2(ab, cd, f, value) {
		// value is in [-2 pi ; 2 pi]; divide by 2 to map to [-pi;pi]
		this.b.push(sin(value)/*, sin(value)*/)
		var ia = this.varMap.get(ab.a),
			ib = this.varMap.get(ab.b),
			ic = this.varMap.get(cd.a),
			id = this.varMap.get(cd.b)
		// calculate the angle for each segment relative to the x axis using atan
		// subtract angle of ab from angle of cd to get signed difference
		// two functions, one calculate sin, one cos of signed difference / 2 to map the signed difference to a point on the unit circle
		// using two functions is necessary so that the resulting functions are "continuous"
		this.F.push(
			(x) => {
				const angle = atan2(
						(x[id + 1] - x[ic + 1]) * f[1],
						(x[id] - x[ic]    ) * f[1])
					- atan2(
						(x[ib + 1] - x[ia + 1]) * f[0],
						(x[ib] - x[ia]    ) * f[0])
				console.log("vaalue", rad2deg(value).toFixed(6), "angle", rad2deg(angle).toFixed(6), cos(angle), sin(angle))
				return sin(angle)
			}/*,
			 (x) => sin(
			 ( atan2(
			 (x[id + 1] - x[ic + 1]) * f[1],
			 (x[id]     - x[ic]    ) * f[1])
			 - atan2(
			 (x[ib + 1] - x[ia + 1]) * f[0],
			 (x[ib]     - x[ia]    ) * f[0])))*/)
	}

	gaussNewtonStep() {
		let DISABLE_CONSOLE = true
		DISABLE_CONSOLE && disableConsole()
		var {F, x, b} = this
		const flatF = x => F.flatMap(f => f(x))
		var Fx = new NLA.Vector(new Float64Array(flatF(x)))
		console.log('F', F)
		console.log("x", x)
		console.log("Fx", Fx.toString())
		console.log("b", b)
		var jacobi = Matrix.jacobi(x => F.flatMap(f => f(x)), x, Fx.v, 1e-6)
		console.log("jacobi\n", jacobi.toString(x=>''+x))
		let jacobiDependentRowIndexes = jacobi.getDependentRowIndexes()
		if (0 != jacobiDependentRowIndexes.length) {
			let error:any = new Error()
			error.jacobiDependentRowIndexes = jacobiDependentRowIndexes
			throw error
		}
		var jacobiTranspose = jacobi.transposed()
		console.log((jacobi.times(jacobiTranspose)).str)
		console.log((jacobi.times(jacobiTranspose)).inversed().str)
		var matrix = jacobiTranspose.times((jacobi.times(jacobiTranspose)).inversed())
		var bVector = new NLA.Vector(new Float64Array(b))
		var xDiff = matrix.timesVector(Fx.minus(bVector))
		console.log("matrix\n", matrix.toString(), "\nFx.minus(bVector)", Fx.minus(bVector).toString(), "\nxDiff", xDiff.toString())
		this.x = new NLA.Vector(new Float64Array(x)).minus(xDiff).v
		Fx = new NLA.Vector(new Float64Array(flatF(x)))
		DISABLE_CONSOLE && enableConsole()
		return Fx.minus(bVector)
	}

	getLoopForSegment(segment):Edge[] {
		var startPoint = segment.b
		var currentPoint = startPoint
		var loop = []
		do {
			let currentSegment = currentPoint.line
			let edge = currentSegment.toBrepEdge()
			if (currentSegment.b != currentPoint) {
				edge = edge.flipped()
			}
			loop.push(edge)
			//console.log(currentPoint.coincidence.cs.filter(function (point) { point instanceof SegmentEndPoint && point != currentPoint; }));
			var otherPointsInCoincidence = currentPoint.coincidence && currentPoint.coincidence.cs.filter(function (point) {
					//console.log("point", point, point instanceof SegmentEndPoint, point != currentPoint);
					return point instanceof SegmentEndPoint && point != currentPoint
				})
			if (!otherPointsInCoincidence || otherPointsInCoincidence.length != 1) {
				throw new Error("The selected segment is not part of an unambiguous loop.")
			}
			var nextSegmentCoincidencePoint = otherPointsInCoincidence[0]
			var nextSegment = nextSegmentCoincidencePoint.line
			currentPoint = nextSegment.getOtherPoint(nextSegmentCoincidencePoint)
		} while (startPoint.coincidence != currentPoint.coincidence)
		return loop
	}

	toSource() {
		return `(function () {
				let sketch = new Sketch('${this.planeRef.ref}')
				let els = sketch.elements = [${this.elements.map(el => el.toSource()).join(',')}]
				sketch.constraints = [${this.constraints.map(el => el.serialize(this.elements)).join(',')}]
				return sketch
			})()`
	}


	recalculate() {
		let sketch = this
		function pushPoint(p) {
			var varMap = sketch.varMap, x = sketch.x
			if (varMap.has(p) /*	|| !p.isConstrained()*/) {
				return
			}
			if (p.coincidence) {
				p.coincidence.cs.forEach(function (p2) {
					varMap.set(p2, xIndex)
				})
			} else {
				varMap.set(p, xIndex)
			}
			x.push(p.x, p.y)
			xIndex += 2
		}
		console.log("recalculating")
		Object.defineProperty(sketch, 'x', {enumerable: false, value: [], writable: true})
		Object.defineProperty(sketch, 'b', {enumerable: false, value: [], writable: true})
		Object.defineProperty(sketch, 'F', {enumerable: false, value: [], writable: true})
		Object.defineProperty(sketch, 'varMap', {enumerable: false, value: new Map(), writable: true})
		// init x to current values
		var xIndex = 0
		sketch.elements.forEach(function (seg) {
			seg.points.forEach(pushPoint)
			if (seg instanceof SketchArc && sketch.varMap.get(seg.a) != sketch.varMap.get(seg.b)) {
				sketch.constrainEqualDistance.apply(sketch, [seg.a, seg.c, seg.b, seg.c].map(point => sketch.varMap.get(point)))
			}
		})
		console.log("varMap", sketch.varMap)
		var fIndexToConstraint = []
		sketch.constraints.forEach(cst => {
			//console.log(cst);
			if (cst.type == "parallel" || cst.type == "colinear" || cst.type == "equalLength") {
				if (cst.fixed) {
					let fixed = cst.fixed.getOrThrow(), lineWC
					if (fixed instanceof Face || fixed instanceof CustomPlane) {
						let plane = fixed instanceof Face ? fixed.surface.plane : fixed.plane
						lineWC = this.plane.intersectionWithPlane(plane)
						if (null == lineWC) throw new Error("no intersection!!")
					} else {
						lineWC = fixed.curve
						if (lineWC.dir1.isParallelTo(this.plane.normal)) {
							if (null == lineWC) throw new Error("line perp to sketch plane")
						}
					}
					for (let j = 0; j < cst.segments.length; j++) {
						if ("colinear" == cst.type) {
							sketch.constrainDistancePointFixedLineWC(cst.segments[j].a, lineWC, 0)
							sketch.constrainDistancePointFixedLineWC(cst.segments[j].b, lineWC, 0)
						} else {
							sketch.constrainAngleSegmentFixedLineWC(cst.segments[j], lineWC, 1)
						}
					}
				} else {
					for (let j = 1; j < cst.segments.length; j++) {
						let first = cst.segments[0], second = cst.segments[j]
						if ("parallel" == cst.type) {
							// assume that max. one element can be a line or a plane
							sketch.constrainAngleSegmentSegment(first, second, 1)
						}
						if ("colinear" == cst.type) {
							sketch.constrainDistancePointSegment(second.a, first, 0)
							sketch.constrainDistancePointSegment(second.b, first, 0)
						}
						if ("equalLength" == cst.type) {
							sketch.constrainEqualDistance.apply(sketch,
								[0, j]
									.map((segmentsIndex) => cst.segments[segmentsIndex].points)
									.concatenated()
									.map((point) => sketch.varMap.get(point)))
						}
					}
				}
			}
			if (cst.type == "pointDistance") {
				sketch.constrainDistancePointPoint(cst.cs[0], cst.cs[1], cst.distance)
			}
			if (cst.type == "pointOnLine" || cst.type == "pointLineDistance" || cst.type == "pointPlaneDistance") {
				let distance = cst.type != "pointOnLine" ? cst.distance : 0
				if (cst.other instanceof NameRef) {
					let fixed = cst.other.getOrThrow()
					if (fixed instanceof Edge) {
						let curve = fixed.curve
						sketch.constrainDistancePointFixedCurveWC(cst.point, curve, distance)
					} else if (fixed instanceof CustomPlane || fixed instanceof Face && fixed.surface instanceof PlaneSurface) {
						let plane = fixed instanceof Face ? fixed.surface.plane : fixed.plane
						let sketchLineWC = sketch.plane.intersectionWithPlane(plane)
						if (null == sketchLineWC) throw new Error("no intersection!!")
						sketch.constrainDistancePointFixedLineWC(cst.point, sketchLineWC, distance)
					} else if (fixed instanceof Face) {
						let surface = fixed.surface
						let isCurves = surface.isCurvesWithPlane(sketch.plane)
						if (isCurves.length == 0) throw new Error('No intersection!')
						if (isCurves.length != 1) throw new Error('Cannot constrain distance point to face/surface with multiple intersection curves')
						let isCurve = isCurves[0]
						sketch.constrainDistancePointFixedCurveWC(cst.point, isCurve, distance)
					} else {
						assert(false)
					}
				} else if(cst.other instanceof SketchLineSeg) {
					sketch.constrainDistancePointSegment(cst.point, cst.other, distance)
				} else if (cst.other instanceof SketchBezier) {
					sketch.constrainDistancePointBezier(cst.point, cst.other, distance)
				} else if (cst.other instanceof SketchArc) {
					sketch.constrainEqualDistance.apply(sketch, [cst.other.a, cst.other.c, cst.point, cst.other.c].map(point => sketch.varMap.get(point)))
				} else {
					assert(false)
				}
			}
			if (cst.type == "angle") {
				sketch.constrainAngleSegmentSegment2(cst.cs[0], cst.cs[1], cst.f, cst.value)
			}
			if (cst.type == "perpendicular") {
				let cosValue = cst.type == "angle" ? cos(cst.value) : 0
				// assume that max. one element can be a line or a plane
				if (cst.other.plane) {
					sketch.constrainAngleSegmentToFixedPlane(cst.segment, cst.other.plane, cosValue)
				} else {
					sketch.constrainAngleSegmentSegment(cst.segment, cst.other, cosValue)
				}
				/*
				 constrainAngle2D.apply(undefined,
				 [0, 1]
				 .map((whichIndex) => cst.constrains[whichIndex].points).concatenated()
				 .map((point) => varMap.get(point))
				 .map((pointXCoord) => [pointXCoord, pointXCoord + 1]).concatenated()
				 .concat(cst.value));
				 */
			}
//			console.log("added constraint, b:", b);
			for (let i = fIndexToConstraint.length; i < sketch.b.length; i++) {
				fIndexToConstraint.push(cst)
			}
		})

		let fxIndexToConstraint = []
		let subFxs = this.F.map(f => f(this.x))
		subFxs.forEach((subFx, fIndex) => {
			let subFxSize = subFx.length || 1
			while (subFxSize--) fxIndexToConstraint.push(fIndexToConstraint[fIndex])
		})



		if (sketch.b.isEmpty()) {
			return
		}
		let count, lastDiffs, lastSize
		try {
			for (count = 0; count < 5; count++) {
				lastDiffs = sketch.gaussNewtonStep()
				lastSize = lastDiffs.length()
				if (lastSize < NLA_PRECISION / 1000) {
					break
				}
			}
			//enableConsole()
			console.log(`broke at ${count}, lastDiffs ${lastDiffs}, lastSize ${lastSize}`)
			// first set all of them to false
			sketch.constraints.forEach(cst => cst.error = false)
			if (lastSize < 0.00001) {
				sketch.reverse()
			} else {
				// then mark every constraint which has a problematic function
				lastDiffs.v.forEach((el, index) => fxIndexToConstraint[index].error |= +!NLA.eq0(el))
			}
		} catch (e) {
			if (!e.jacobiDependentRowIndexes) throw e // wrong type of error
			console.log("PROBLEMS", e.jacobiDependentRowIndexes)
			e.jacobiDependentRowIndexes.forEach(fxIndex => fxIndexToConstraint[fxIndex].error = true)
		}
	}
	reverse() {
		console.log("REVERSING")
		this.varMap.forEach((pointIndex, point) => {
			point.x = this.x[pointIndex]
			point.y = this.x[pointIndex + 1]
		})
	}
}
NLA.registerClass(Sketch)

class SegmentEndPoint {
	x: number
	y: number
	line: SketchSegment
	id: number
	name: string
	coincidence: SegmentEndPoint[] | null

	constructor(x, y, line) {
		this.x = x
		this.y = y
		this.line = line
		this.id = globalId++
		this.name = "segEndPoint" + (this.id)
		this.coincidence = null
	}

	distanceToCoords(p) {
		return distance(this.x, this.y, p.x, p.y)
	}

	toString() {
		return "SegmentEndPoint #" + this.id
	}

	V3() {
		return V(this.x, this.y, 0)
	}

	isConstrained(sketch) {
		return sketch.constraints.some(constraint => constraint.constrains(this) || constraint.constrains(this.line))
	}

	freeFromConstraints(sketch) {
		sketch.constraints.forEach(constraint => constraint.constrains(this) && removeFromConstraint(this, sketch, constraint))
	}

	/**
	 * Can be called even if already removed
	 * @param sketch
	 */
	removeFromSketch(sketch) {
		sketch.removeSegment(this.line)
	}

	moveCoincidence(v3) {
		if (this.coincidence) {
			this.coincidence.cs.forEach(p => {
				p.x = v3.x
				p.y = v3.y
			})
		} else {
			this.x = v3.x
			this.y = v3.y
		}
	}

	get sketch() {
		return this.line.sketch
	}

	canon() {
		return this.coincidence && this.coincidence.cs[0] || this
	}

	static fromV3(segment, p) {
		return new SegmentEndPoint(p.x, p.y, segment)
	}

}
NLA.registerClass(SegmentEndPoint)

class SketchArc {
	a: SegmentEndPoint
	b: SegmentEndPoint
	c: SegmentEndPoint
	sketch: Sketch
	points: SegmentEndPoint[]
	id: number
	name: string

	/**
	 * Arc goes CCW from a to b
	 * @param sketch
	 * @param ax
	 * @param ay
	 * @param bx
	 * @param by
	 * @param cx
	 * @param cy
	 * @constructor
	 */
	constructor(sketch, ax?, ay?, bx?, by?, cx?, cy?) {
		assertInst(Sketch, sketch)
		this.sketch = sketch
		this.a = new SegmentEndPoint(ax, ay, this)
		this.b = new SegmentEndPoint(bx, by, this)
		this.c = new SegmentEndPoint(cx, cy, this)
		this.points = [this.c, this.a, this.b]
		this.id = globalId++
		this.name = "SketchArc" + this.id
	}

	angleA() {
		return this.a.V3().minus(this.c.V3()).angleXY()
	}

	angleB() {
		return this.b.V3().minus(this.c.V3()).angleXY()
	}

	radiusA() {
		return this.a.distanceToCoords(this.c)
	}

	distanceToCoords(coords) {
		coords = V(coords)
		var angleA = this.angleA(), angleB = this.angleB()
		if (angleB <= angleA) {
			angleB += Math.PI * 2
		}
		var relCoords = coords.minus(this.c.V3()), angle = relCoords.angleXY(), radius = this.radiusA()
		if (angle < 0) {
			angle += Math.PI * 2
		}
		var angle2 = NLA.clamp(angle, angleA, angleB)
		return V(radius * cos(angle2), radius * sin(angle2), 0).minus(relCoords).length()
	}

	flip() {
		[this.a, this.b] = [this.b, this.a]
	}

	getIntermediatePoints() {
		var result = []
		var angleA = this.angleA(), angleB = this.angleB()
		if (angleB <= angleA) {
			angleB += Math.PI * 2
		}
		var radius = this.radiusA()
		var center = this.c.V3()
		var segmentLength = radius * (angleB - angleA), pointCount = floor(segmentLength / 10)
		var intervalAngle = (angleB - angleA) / (pointCount + 1)
		for (var i = 1; i < pointCount + 1; i++) {
			var angle = angleA + i * intervalAngle
			result.push(center.plus(V3.polar(radius, angle)))
		}
		return result
	}

	getVectorCA() {
		return this.a.V3().minus(this.c.V3())
	}

	getVectorCB() {
		return this.b.V3().minus(this.c.V3())
	}

	getOtherPoint(p) {
		if (p == this.a) return this.b
		if (p == this.b) return this.a
		assert(false)
	}

	toBrepEdge() {
		let curve = this.getCurve()
		return new PCurveEdge(curve,
			this.a.V3(), this.b.V3(),
			-PI, curve.pointLambda(this.b.V3()),
			null,
			curve.tangentAt(-PI), curve.tangentAt(curve.pointLambda(this.b.V3())),
			this.name + 'Edge')
	}

	getCurve() {
		let ca = this.getVectorCA()
		return new EllipseCurve(this.c.V3(), ca.negated(), ca.negated().getPerpendicular())
	}
}
NLA.registerClass(SketchArc)

class SketchLineSeg {
	sketch: Sketch
	a: SegmentEndPoint
	b: SegmentEndPoint
	points: SegmentEndPoint[]
	id: number
	name: string

	constructor(sketch, x1?, y1?, x2?, y2?) {
		assertInst(Sketch, sketch)
		this.sketch = sketch
		this.points = [new SegmentEndPoint(x1, y1, this), new SegmentEndPoint(x2, y2, this)]
		this.a = this.points[0]
		this.b = this.points[1]
		this.id = globalId++
		this.name = "segment" + this.id
	}

	remove() {
		assert (editingSketch.elements.includes(this))
		this.removeFromSketch(editingSketch)
	}

	distanceToCoords(coords) {
		return this.distanceTo(coords.x, coords.y)
	}

	angleTo(segment) {
		assertInst(SketchLineSeg, segment)
		return segment.angleAB() - this.angleAB()
	}

	toString() {
		return "SketchLineSeg #" + this.id
	}

	angleAB() {
		return atan2(this.b.y - this.a.y, this.b.x - this.a.x)
	}


	getOtherPoint(p) {
		if (p == this.a) return this.b
		if (p == this.b) return this.a
		assert(false)
	}

	pointLambda(v) {
		if (this.b.x - this.a.x > this.b.y - this.a.y) {
			return (v.x - this.a.x) / (this.b.x - this.a.x)
		} else {
			return (v.y - this.a.y) / (this.b.y - this.a.y)
		}
	}

	distanceTo(x, y) {
		var x1 = this.a.x, y1 = this.a.y, x2 = this.b.x,  y2 = this.b.y
		var a = y1 - y2
		var b = x2 - x1
		var c = x2 * y1 - x1 * y2
		var dist = Math.abs(a * x + b * y - c) / length(a, b)
		var xClosest = (b * (b * x - a * y) + a * c) / lengthSquared(a, b)
		var yClosest = (a * ( -b * x + a * y) + b * c) / lengthSquared(a, b)
		if (x1 != x2 ? isBetween(xClosest, x1, x2) : isBetween(yClosest, y1, y2)) {
			return dist
		} else {
			if (x1 < x2 && x < x1 || x1 > x2 && x > x1
				|| x1 == x2 && (y1 < y2 && y < y2 || y1 > y2 && y > y1)) {
				//noinspection JSSuspiciousNameCombination
				return distance(x, x1, y, y1)
			} else {
				//noinspection JSSuspiciousNameCombination
				return distance(x, x2, y, y2)
			}
		}
	}

	getVectorAB() {
		return V(this.b.x - this.a.x, this.b.y - this.a.y, 0)
	}

	getClosestPoint(x, y) {
		var x1 = this.a.x, y1 = this.a.y, x2 = this.b.x, y2 = this.b.y
		var a = y1 - y2
		var b = x2 - x1
		var c = x2 * y1 - x1 * y2
		var dist = abs(a * x + b * y - c) / length(a, b)
		var xClosest = (b * (b * x - a * y) + a * c) / lengthSquared(a, b)
		var yClosest = (a * ( -b * x + a * y) + b * c) / lengthSquared(a, b)
		if (isBetween(xClosest, x1, x2)) {
			return {"x": xClosest, "y": yClosest}
		} else {
			if (x1 < x2 && x < x1 || x1 > x2 && x > x1
				|| x1 == x2 && (y1 < y2 && y < y2 || y1 > y2 && y > y1)) {
				return {x: x1, y: y1}
			} else {
				//noinspection JSSuspiciousNameCombination
				return {x: y1, y: y2}
			}
		}
	}

	length() {
		return this.points[0].distanceTo(this.points[1])
	}

	intersection(segment) {
		return intersection(this.a.x, this.a.y, this.b.x, this.b.y,
			segment.a.x, segment.a.y, segment.b.x, segment.b.y)
	}

	getCurve() {
		return L3.anchorDirection(this.a.V3(), this.getVectorAB())
	}

	toBrepEdge() {
		return StraightEdge.throughPoints(this.a.V3(), this.b.V3(), this.name + 'Edge')
	}
}
NLA.registerClass(SketchLineSeg)

class Constraint {
	type: string
	id: int
	cs: any[]

	constructor(type, constrains, props?) {
		if (constrains.constructor != Array) {
			throw new Error("not an array: " + constrains)
		}
		this.type = type
		this.id = globalId++
		this.cs = constrains
		for (var key in props) {
			this[key] = props[key]
		}

	}

	constrains(o) {
		if (!this.cs) {
			console.log(this)
		}
		return this.cs.includes(o)
	}

	segmentOtherTypeFree(sketch, el) {

	}
}
NLA.registerClass(Constraint)


//////////////////////////////////
///////// functions! /////////////
//////////////////////////////////
function intersection(x1, y1, x2, y2, x3, y3, x4, y4) {
	var denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
	var xNominator = (x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)
	var yNominator = (x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)
	return V(xNominator / denominator, yNominator / denominator, 0)
}

function isBetween(val, a, b) {
	if (a < b) {
		return a < val && val < b
	} else {
		return b < val && val < a
	}
}

function distance(x1, y1, x2, y2) {
	return length(x2 - x1, y2 - y1)
}
function distanceSquared(x1, y1, x2, y2) {
	return lengthSquared(x2 - x1, y2 - y1)
}
function bezierCurveAt(x0, y0, x1, y1, x2, y2, x3, y3, t) {
	//console.log(x0, y0, x1, y1, x2, y2, x3, y3, t)
	let s = 1 - t, c0 = s * s * s, c1 = 3 * s * s * t, c2 = 3 * s * t * t, c3 = t * t * t
	return {
		x: x0 * c0 + x1 * c1 + x2 * c2 + x3 * c3,
		y: y0 * c0 + y1 * c1 + y2 * c2 + y3 * c3
	}
}
function bezierCurveTangentAt(x0, y0, x1, y1, x2, y2, x3, y3, t) {
	//console.log(x0, y0, x1, y1, x2, y2, x3, y3, t)
	let s = 1 - t, c01 = 3 * s * s, c12 = 6 * s * t, c23 = 3 * t * t
	return {
		x: (x1 - x0) * c01 + (x2 - x1) * c12 + (x3 - x2) * c23,
		y: (y1 - y0) * c01 + (y2 - y1) * c12 + (y3 - y2) * c23,
	}
}
function length(x, y) {
	return Math.sqrt(x * x + y * y)
}
function alpha(x1, y1, x2, y2) {
	return length(x2 - x1, y2 - y1)
}

function lengthSquared(x, y) {
	return x * x + y * y
}

function distanceLinePoint(x1, y1, x2, y2, x, y) {
	let a = y1 - y2
	let b = x2 - x1
	let c = x2 * y1 - x1 * y2
	let dist = Math.abs(a * x + b * y - c) / length(a, b)
	return dist
}


function distanceLinePointSigned(x1, y1, x2, y2, x, y) {
	// function needs to be differentiable around target value
	// therefore this funciton is necessary
	let a = y1 - y2
	let b = x2 - x1
	let c = x2 * y1 - x1 * y2
	//TODOlet dist = Math.abs(a * x + b * y - c) / length(a, b);
	let dist = (a * x + b * y - c) / length(a, b)
	return dist
}

function angleVectors(ax, ay, bx, by) {
//	console.log(ax, ay, bx, by);
	return Math.abs(ax * bx + ay * by) / Math.sqrt((ax * ax + ay * ay) * (bx * bx + by * by))
}

// returns angle between segments AB, CD
function angleABCD(ax, ay, bx, by, cx, cy, dx, dy) {
	return angleVectors(bx - ax, by - ay, dx - cx, dy - cy)
	// return ((bx - ax) * (dx - cx) + (by - ay) * (dy - cy)) / Math.sqrt(((bx - ax)*(bx - ax)+(by - ay)*(by - ay))*((dx - cx)*(dx - cx)+(dy - cy)*(dy - cy)))
}