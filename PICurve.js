/**
 * Created by aval on 19/02/2016.
 */
function CurvePI(parametricSurface, implicitSurface, startPoint) {
	assert (parametricSurface.parametricFunction, 'parametricSurface.parametricFunction')
	assert(implicitSurface.implicitFunction, 'implicitSurface.implicitFunction')
	this.parametricSurface = parametricSurface
	this.implicitSurface = implicitSurface
	if (!startPoint) {
		var pmPoint = curvePoint(this.implicitCurve(), V3(1, 1, 0))
		this.startPoint = this.parametricSurface.parametricFunction()(pmPoint.x, pmPoint.y)
	} else {
		this.startPoint = startPoint
	}
	this.isLoop = false
	this.calcPoints(this.startPoint)
}
var STEP_SIZE = 1
CurvePI.prototype = NLA.defineObject(null, {
	implicitCurve: function () {
		var pF = this.parametricSurface.parametricFunction()
		var iF = this.implicitSurface.implicitFunction()
		return function (s, t) {
			return iF(pF(s, t))
		}
	},
	containsPoint: function (p) {
		assertVectors(p)
		return this.parametricSurface.containsPoint(p) && isZero(this.implicitSurface.implicitFunction()(p))
	},
	getVerticesNo0: function () {
		function sliceCyclic(arr, start, end) {
			if (start <= end) {
				return arr.slice(start, end)
			} else {
				return arr.slice(start).concat(arr.slice(0, end))
			}
		}

		// TODOOO
		if (!this.canon) {
			var start = floor(this.aT + 1), end = ceil(this.bT)
			var arr = sliceCyclic(this.curve.points, start, end)
		} else {
			var start = floor(this.bT + 1), end = ceil(this.aT)
			var arr = sliceCyclic(this.curve.points, start, end)
			console.log("this.canon", !!this.canon, arr.length, start, end, this.aT)
			arr.reverse()
		}
		arr.push(this.b)
		return arr
	},
	calcPoints: function (curveStartPoint) {
		if (!this.points) {
			var pF = this.parametricSurface.parametricFunction()
			var iF = this.implicitSurface.implicitFunction()
			var iBounds = this.implicitSurface.boundsFunction()
			var curveFunction = (s, t) => iF(pF(s, t))
			var pTPF = this.parametricSurface.pointToParameterFunction()
			var startParams = pTPF(this.startPoint)
			this.pmTangentEndPoints = []
			this.pmPoints = followAlgorithm(curveFunction, startParams, startParams, STEP_SIZE, null,
				this.pmTangentEndPoints, (s, t) => iBounds(pF(s, t)))
			this.isLoop = this.pmPoints[0].distanceTo(this.pmPoints[this.pmPoints.length - 1]) < STEP_SIZE * 1.1
			this.startT = 0
			if (!this.isLoop) {
				// the curve starting at curveStartPoint is not closed, so we need to find curve points in the other
				// direction until out of bounds
				var pmTangent0 = this.pmTangentEndPoints[0].minus(this.pmPoints[0])
				var pmTangentEndPoints2 = []
				var pmPoints2 = followAlgorithm(curveFunction, startParams, startParams, STEP_SIZE, pmTangent0.negated(),
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
			this.tangents = this.pmTangentEndPoints.map(
				({x: d, y: z}, i, ps) => pF(d, z).minus(this.points[i]))
			//console.log('this.points', this.points.map(v => v.ss).join(", "))
			this.startTangent = this.tangentAt(this.startT)
		}
	},
	pointTangent: function (point) {
		assertVectors(point)
		assert(this.containsPoint(point), 'this.containsPoint(point)'+this.containsPoint(point))
		this.calcPoints(point)
		var pIndex = this.pointLambda(point)
		return this.tangents[pIndex]
	},
	tangentAt: function (t) {
		return this.tangents[Math.round(t)]
	},
	pointLambda: function (point) {
		assertVectors(point)
		assert(this.containsPoint(point), 'this.containsPoint(p)')
		var pmPoint = this.parametricSurface.pointToParameterFunction()(point)
		var ps = this.points, pmps = this.pmPoints, t = 0, prevDistance, pmDistance = pmPoint.distanceTo(pmps[0])
		while (pmDistance > STEP_SIZE && t < ps.length - 1) { // TODO -1?
			//console.log(t, pmps[t].ss, pmDistance)
			t += Math.min(1, Math.round(pmDistance / STEP_SIZE / 2))
			pmDistance = pmPoint.distanceTo(pmps[t])
		}
		if (t >= ps.length - 1) {
			// point is not on this curve
			return NaN
		}
		if (ps[t].like(point)) return t
		var nextT = (t + 1) % ps.length, prevT = (t + ps.length - 1) % ps.length
		if (ps[nextT].distanceTo(point) < ps[prevT].distanceTo(point)) {
			return t + 0.4
		} else {
			return t - 0.4
		}
	}
})
function CurvePILoop(curve, startPoint) {
	assert(curve instanceof CurvePI)
	this.curve = curve
	assert(this.curve.isLoop)
	this.a = this.b = this.startPoint = startPoint
}
CurvePILoop.prototype = NLA.defineObject(B2.Edge.prototype, {
	getVerticesNo0: function () {
		this.curve.calcPoints()
		this.points = this.curve.points
		return this.curve.points
	},
	getIntersectionsWithPSurface: function (pSurface) {
		assert (pSurface.parametricFunction)
	},
	tangentAt: function (p) {
		return this.curve.tangentAt(p)
	},
	isCCW: function (normal) {
		var step = floor(this.points.length / 4), verts = NLA.arrayFromFunction(4, i => this.points[step * i])
		return isCCW(verts, normal)
	}
})