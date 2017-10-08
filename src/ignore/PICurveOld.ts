const STEP_SIZE = 1

class PICurveOld extends Curve {
	parametricSurface: Surface
	implicitSurface: Surface
	startPoint: V3
	isLoop: boolean

	constructor(parametricSurface: Surface, implicitSurface: Surface, startPoint: V3, endPoint?: V3) {
		super()
		assert(parametricSurface.parametricFunction, 'parametricSurface.parametricFunction')
		assert(implicitSurface.implicitFunction, 'implicitSurface.implicitFunction')
		this.parametricSurface = parametricSurface
		this.implicitSurface = implicitSurface
		if (!startPoint) {
			const pmPoint = curvePoint(this.implicitCurve(), V(1, 1, 0))
			this.startPoint = this.parametricSurface.parametricFunction()(pmPoint.x, pmPoint.y)
		} else {
			this.startPoint = startPoint
		}
		this.isLoop = false
		this.calcPoints(this.startPoint)
	}

	implicitCurve() {
		const pF = this.parametricSurface.parametricFunction()
		const iF = this.implicitSurface.implicitFunction()
		return function (s, t) {
			return iF(pF(s, t))
		}
	}

	containsPoint(p) {
		assertVectors(p)
		return this.parametricSurface.containsPoint(p) && eq0(this.implicitSurface.implicitFunction()(p))
	}

	getVerticesNo0() {
		function sliceCyclic(arr, start, end) {
			if (start <= end) {
				return arr.slice(start, end)
			} else {
				return arr.slice(start).concat(arr.slice(0, end))
			}
		}

		// TODOOO
		let start, end, arr
		if (!this.canon) {
			start = Math.floor(this.aT + 1)
			end = ceil(this.bT)
			arr = sliceCyclic(this.curve.points, start, end)
		} else {
			start = Math.floor(this.bT + 1)
			end = ceil(this.aT)
			arr = sliceCyclic(this.curve.points, start, end)
			console.log('this.canon', !!this.canon, arr.length, start, end, this.aT)
			arr.reverse()
		}
		arr.push(this.b)
		return arr
	}

	calcPoints() {
		if (!this.points) {
			const pF = this.parametricSurface.parametricFunction()
			const iF = this.implicitSurface.implicitFunction()
			const iBounds = this.implicitSurface.boundsFunction()
			const curveFunction = (s, t) => iF(pF(s, t))
			const pTPF = this.parametricSurface.pointToParameterFunction()
			const startParams = pTPF(this.startPoint)
			this.pmTangentEndPoints = []
			this.pmPoints = followAlgorithm(curveFunction, startParams, startParams, STEP_SIZE, undefined,
				this.pmTangentEndPoints, (s, t) => iBounds(pF(s, t)))
			this.isLoop = this.pmPoints[0].distanceTo(this.pmPoints[this.pmPoints.length - 1]) < STEP_SIZE * 1.1
			this.startT = 0
			if (!this.isLoop) {
				// the curve starting at curveStartPoint is not closed, so we need to find curve points in the other
				// direction until out of bounds
				const pmTangent0 = this.pmTangentEndPoints[0].minus(this.pmPoints[0])
				let pmTangentEndPoints2 = []
				const pmPoints2 = followAlgorithm(curveFunction, startParams, startParams, STEP_SIZE, pmTangent0.negated(),
					pmTangentEndPoints2, (s, t) => iBounds(pF(s, t)))
				pmTangentEndPoints2 = pmTangentEndPoints2.map((ep, i) => pmPoints2[i].times(2).minus(ep))
				this.startT = pmPoints2.length
				pmPoints2.reverse()
				pmPoints2.pop()
				this.pmPoints = pmPoints2.concat(this.pmPoints)
				pmTangentEndPoints2.reverse()
				pmTangentEndPoints2.pop()
				this.pmTangentEndPoints = pmTangentEndPoints2.concat(this.pmTangentEndPoints)
			}
			this.points = this.pmPoints.map(({x: d, y: z}) => pF(d, z))
			this.tangents = this.pmTangentEndPoints.map(({x: d, y: z}, i, ps) => pF(d, z).minus(this.points[i]),
			)
			//console.log('this.points', this.points.map(v => v.$).join(", "))
			this.startTangent = this.tangentAt(this.startT)
		}
	}

	pointTangent(point) {
		assertVectors(point)
		assert(this.containsPoint(point), 'this.containsPoint(point)' + this.containsPoint(point))
		this.calcPoints(point)
		const pIndex = this.pointT(point)
		return this.tangents[pIndex]
	}

	tangentAt(t: number): V3 {
		return this.tangents[Math.round(t)]
	}

	pointT(point) {
		assertVectors(point)
		assert(this.containsPoint(point), 'this.containsPoint(p)')
		const pmPoint = this.parametricSurface.pointToParameterFunction()(point)
		let ps = this.points, pmps = this.pmPoints, t = 0, prevDistance, pmDistance = pmPoint.distanceTo(pmps[0])
		while (pmDistance > STEP_SIZE && t < ps.length - 1) { // TODO -1?
			//console.log(t, pmps[t].$, pmDistance)
			t += Math.min(1, Math.round(pmDistance / STEP_SIZE / 2))
			pmDistance = pmPoint.distanceTo(pmps[t])
		}
		if (t >= ps.length - 1) {
			// point is not on this curve
			return NaN
		}
		if (ps[t].like(point)) return t
		const nextT = (t + 1) % ps.length, prevT = (t + ps.length - 1) % ps.length
		if (ps[nextT].distanceTo(point) < ps[prevT].distanceTo(point)) {
			return t + 0.4
		} else {
			return t - 0.4
		}
	}
}