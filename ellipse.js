/**
 * Created by aval on 14/02/2016.
 */

"use strict"
var EllipseCurve = NLA.defineClass('EllipseCurve', null,
	/** @constructor */
	function (center, f1, f2) {
		assertVectors(center, f1, f2)
		this.center = center
		this.f1 = f1
		this.f2 = f2
		this.normal = f1.cross(f2).normalized()
		this.matrix = M4.forSys(f1, f2, this.normal, center)
		this.inverseMatrix = this.matrix.inversed()
	},
	{
		toString: function (f) {
			return `new EllipseCurve(${this.center}, ${this.f1}, ${this.f2})`
		},
		isValidT: (t) => -PI <= t && t <= PI,
		asklkjas: function (aT, bT, a, b, reversed, includeFirst) {
			var split = 4 * 17, inc = 2 * PI / split
			var verts = []
			if (includeFirst) verts.push(a)
			console.log("revrs", reversed)
			if (!reversed) {
				assert(aT < bT)
				var start = ceil((aT + NLA.PRECISION) / inc)
				var end = floor((bT - NLA.PRECISION) / inc)
				console.log(aT, bT, start, end, inc)
				for (var i = start; i <= end; i++) {
					verts.push(this.at(i * inc))
				}
			} else {
				assert(bT < aT)
				var start = floor((aT - NLA.PRECISION) / inc)
				var end = ceil((bT + NLA.PRECISION) / inc)
				for (var i = start; i >= end; i--) {
					verts.push(this.at(i * inc))
				}
			}
			verts.push(b)
			return verts
		},
		at: function (t) {
			return this.center.plus(this.f1.times(cos(t))).plus(this.f2.times(sin(t)))
		},
		at2: function (xi, eta) {
			return this.center.plus(this.f1.times(xi)).plus(this.f2.times(eta))
		},
		tangentAt: function (t) {
			assertNumbers(t)
			return this.f2.times(cos(t)).minus(this.f1.times(sin(t))).normalized()
		},
		tangentAt2: function (xi, eta) {
			return this.f2.times(xi).minus(this.f1.times(eta)).normalized()
		},
		isCircular: function () {
			return NLA.equals(this.f1.length(), this.f2.length())
		},
		equals: function (curve) {
			return curve.constructor == EllipseCurve
				&& this.center.like(curve.center)
				&& this.f1.like(curve.f1)
				&& this.f2.like(curve.f2)
		},
		isColinearTo: function (curve) {
			if (curve.constructor != EllipseCurve) { return false }
			if (!this.center.like(curve.center)) { return false }
			if (this == curve) { return true }
			if (this.isCircular()) {
				return curve.isCircular() && NLA.equals(this.f1.length(), curve.f1.length())
			} else {
				var {f1: f1, f2: f2} = this.mainAxes(), {f1: c1, f2: c2} = curve.mainAxes()
				return NLA.equals(f1.lengthSquared(), abs(f1.dot(c1)))
					&& NLA.equals(f2.lengthSquared(), abs(f2.dot(c2)))
			}
		},
		normalAt: function (t) {
			return this.tangentAt(t).cross(this.normal)
		},
		pointLambda: function (p, hint) {
			assertVectors(p)
			var p2 = this.inverseMatrix.transformPoint(p)
			var angle = p2.angleXY()
			if (angle < -Math.PI + NLA.PRECISION || angle > Math.PI - NLA.PRECISION) {
				return hint.dot(this.f2) < 0
					? Math.PI
					: -Math.PI
			}
			return angle
		},
		isOrthogonal: function (p) {
			return this.f1.isPerpendicularTo(this.f2)
		},
		mainAxes: function () {
			var f1 = this.f1, f2 = this.f2, a = f1.dot(f2), b = f2.lengthSquared() - f1.lengthSquared()
			if (NLA.isZero(a)) {
				return this
			}
			var g1 = 2 * a, g2 = b +     sqrt(b * b + 4 * a * a)
			// TODO
			// g1 * xi + g2 * eta = 0 (1)
			// xi² + eta² = 1         (2)
			// (1) => eta = -g1 * xi / g2
			// => xi² + g1² * xi² / g2² = 1 = xi² * (1 + g1² / g2²)
			// => xi = +-sqrt(1 / (1 + g1² / g2²))
			// => eta = +-sqrt(1 / (1 + g2² / g1²))
			/*
			 var xi = -sqrt(1 / (1 + g1 * g1 / g2 / g2))
			 var eta = sqrt(1 / (1 + g2 * g2 / g1 / g1))
			 return {f1: f1.times(xi).plus(f2.times(eta)), f2: f1.times(-eta).plus(f2.times(xi))}
			 */
			var {x1: xi, y1: eta} = intersectionUnitCircleLine(g1, g2, 0)
			return {f1: f1.times(xi).plus(f2.times(eta)), f2: f1.times(-eta).plus(f2.times(xi))}

		},
		circumference: function (p) {
			assert(false, "implementation?")
		},
		circumferenceApproximate: function (p) {
			// approximate circumference by Ramanujan
			// https://en.wikipedia.org/wiki/Ellipse#Circumference
			var {f1, f2} = this.mainAxes(), a = f1.length(), b = f2.length()
			var h = (a - b) * (a - b) / (a + b) / (a + b)
			return PI * (a + b) * (1 + 3 * h / (10 + sqrt(4 - 3 * h)))
		},
		transform: function (m4) {
			return new EllipseCurve(m4.transformPoint(this.center), m4.transformVector(this.f1), m4.transformVector(this.f2))
		},
		rightAngled: function () {
			var {f1, f2} = this.mainAxes()
			return new EllipseCurve(this.center, f1, f2)
		},
		projectedOn: function (plane) {
			assert (plane instanceof P3)
			return new EllipseCurve()
		},
		debugToMesh: function (mesh, bufferName) {
			mesh.addVertexBuffer(bufferName, bufferName)
			for (var t = 0; t < 2 * PI; t+=0.1) {
				var p = this.at(t);
				mesh[bufferName].push(p, p.plus(this.tangentAt(t).toLength(1)))
				mesh[bufferName].push(p, p.plus(this.normalAt(t).toLength(1)))
			}
			mesh[bufferName].push(this.center, this.center.plus(this.f1.times(1.2)))
			mesh[bufferName].push(this.center, this.center.plus(this.f2))
			mesh[bufferName].push(this.center, this.center.plus(this.normal))
		},
		getIntersectionsWithSurface: function (surface) {
			if (surface instanceof PlaneSurface) {
				return this.getIntersectionsWithPlane(surface.plane)
			} else {
				assert (false)
			}
		},
		getIntersectionsWithPlane: function (plane) {
			assert (plane instanceof P3)
			/*
			 this: x = center + f1 * cos t + f2 * sin t  (1)
			 plane:
			 n := plane.normal
			 n DOT x == plane.w           (2)
			 plane defined by f1/f2
			 x = center + f1 * xi + f2 * eta         (3)
			 intersection plane and planef1/f2:
			 insert (3) into (2):
			 n DOT center + n DOT f1 * xi + n DOT f2 * eta = plane.w | -n DOT center
			 n DOT f1 * xi + n DOT f2 * eta = plane.w - n DOT center (4)
			 points on ellipse have additional condition
			 eta * eta + xi * xi = 1 (5)
			 g1 := n DOT f1
			 g2 := n DOT f2
			 g3 := w - n DOT center
			 solve system (5)/(6)
			 g1 * eta + g2 * eta = g3 (6)
			 */
			if (plane.normal.isParallelTo(this.normal)) {
				return []
			}
			var
				n = plane.normal, w = plane.w,
				center = this.center, f1 = this.f1, f2 = this.f2,
				g1 = n.dot(f1), g2 = n.dot(f2), g3 = w - n.dot(center)

			var {x1: xi1, y1: eta1, x2: xi2, y2: eta2} = intersectionUnitCircleLine(g1, g2, g3)
			if (isNaN(xi1)) {
				return []
			}
			return [atan2(eta1, xi1), atan2(eta2, xi2)]

		},
		getPlane: function () {
			return P3.normalOnAnchor(this.normal, this.center)
		},
		containsPoint: function (p) {
			var localP = this.inverseMatrix.transformPoint(p)
			return NLA.equals(1, localP.lengthXY()) && NLA.isZero(localP.z)
		},
		intersectWithEllipse: function (ellipse) {
			if (this.normal.isParallelTo(ellipse.normal) && NLA.isZero(this.center.minus(ellipse.center).dot(ellipse.normal))) {
				// ellipses are coplanar
				var localEllipse = ellipse.transform(this.inverseMatrix)

			} else {
				assert(false)
			}
		}
	},
	{
		forAB: function (a, b) {
			return new EllipseCurve(V3.ZERO, V3(a, 0, 0), V3(0, b, 0))
		}
	}
)
var CylinderSurface = NLA.defineClass('CylinderSurface', null,
	/** @constructor */
	function (baseEllipse, dir, normalDir) {
		assertVectors(dir)
		assert(1 == normalDir || -1 == normalDir, "normalDir == 1 || normalDir == -1" + normalDir)
		assert(baseEllipse instanceof EllipseCurve)
		//assert(!baseEllipse.normal.isPerpendicularTo(dir), !baseEllipse.normal.isPerpendicularTo(dir))
		assert(dir.hasLength(1))
		this.baseEllipse = baseEllipse
		this.dir = dir
		this.normalDir = normalDir
		this.matrix = M4.forSys(baseEllipse.f1, baseEllipse.f2, dir, baseEllipse.center)
		this.inverseMatrix = this.matrix.inversed()
	},
	{
		edgeLoopContainsPoint: function (contour, p) {
			assertVectors(p)
			var line = L3(p, this.dir)
			// create plane that goes through cylinder seam
			var seamBase = this.baseEllipse.at(PI)
			var intersectionLinePerpendicular = this.dir.cross(p.minus(seamBase))
			var plane2 = P3.normalOnAnchor(intersectionLinePerpendicular, p)
			var colinearSegments = contour.map((edge) => edge.colinearToLine(line))
			var colinearSegmentsInside = contour.map((edge, i) => edge.aDir.dot(this.dir) > 0)
			var inside = false

			function logIS(p) {
				if (line.pointLambda(p) > 0) {
					inside = !inside
				}
			}

			contour.forEach((edge, i, edges) => {
				var j = (i + 1) % edges.length, nextEdge = edges[j]
				//console.log(edge.toSource()) {p:V3(2, -2.102, 0),
				if (colinearSegments[i]) {
					// edge colinear to intersection
					var outVector = edge.bDir.cross(plane.normal)
					var insideNext = outVector.dot(nextEdge.aDir) > 0
					if (colinearSegmentsInside[i] != insideNext) {
						logIS(edge.b)
					}
				} else {
					var edgeTs = edge.getIntersectionsWithPlane(plane2)
					for (var k = 0; k < edgeTs.length; k++) {
						var edgeT = edgeTs[k]
						if (edgeT == edge.bT) {
							// endpoint lies on intersection line
							if (colinearSegments[j]) {
								// next segment is colinear
								// we need to calculate if the section of the plane intersection line BEFORE the colinear segment is
								// inside or outside the face. It is inside when the colinear segment out vector and the current segment vector
								// point in the same direction (dot > 0)
								var colinearSegmentOutsideVector = nextEdge.aDir.cross(plane.normal)
								var insideFaceBeforeColinear = colinearSegmentOutsideVector.dot(edge.bDir) < 0
								// if the "inside-ness" changes, add intersection point
								//console.log("segment end on line followed by colinear", insideFaceBeforeColinear != colinearSegmentInsideFace, nextSegmentOutsideVector)
								if (colinearSegmentsInside[j] != insideFaceBeforeColinear) {
									logIS(edge.b)
								}
							} else if (intersectionLinePerpendicular.dot(edge.bDir) * intersectionLinePerpendicular.dot(nextEdge.aDir) > 0) {
								logIS(edge.b)
							}
						} else if (edgeT != edge.aT) {
							// edge crosses line, neither starts nor ends on it
							logIS(edge.curve.at(edgeT))
						}
					}
				}
			})
			return inside

		},
		toString: function () {
			return "CylinderSurface"
		},
		intersectionWithLine: function (line) {
			// fun fun fun
			var lineLocal = line.transform(this.inverseMatrix)
			var {x: dx, y : dy, z: dz} = lineLocal.dir1
			var {x: ax, y : ay, z: az} = lineLocal.anchor
			var a = - dy, b = dx, c = a * ax + b * ay
			var {x1, y1, x2, y2} = intersectionUnitCircleLine(a, b, c)
			var localLambda1 = ((x1 - ax) * dx + (y1 - ay) * dy) / (dx * dx + dy * dy)
			var localLambda2 = ((x2 - ax) * dx + (y2 - ay) * dy) / (dx * dx + dy * dy)
			return this.matrix.transformedPoints([
				V3.create(x1, y1, az + localLambda1 * dz),
				V3.create(x2, y2, az + localLambda2 * dz)])
		},
		coplanarTo: function (surface) {
			return this == surface ||
				surface instanceof CylinderSurface
				&& this.containsEllipse(surface.baseEllipse)
				&& this.dir.isParallelTo(surface.dir)
		},
		containsEllipse: function (ellipse) {
			var ellipseProjected = ellipse.transform(M4.projection(this.baseEllipse.getPlane(), this.dir))
			return this == ellipse || this.baseEllipse.isColinearTo(ellipseProjected)
		},
		containsLine: function (line) {
			return this.dir.isParallelTo(line.dir1) && this.containsPoint(line.anchor)
		},
		containsCurve: function (curve) {
			if (curve instanceof EllipseCurve) {
				return this.containsEllipse(curve)
			} else if (curve instanceof L3) {
				return this.containsLine(curve)
			} else {
				assert(false)
			}
		},
		transform: function (m4) {
			return new CylinderSurface(
				this.baseEllipse.transform(m4),
				m4.transformVector(this.dir),
				this.normalDir)
		},
		flipped: function () {
			return new CylinderSurface(
				this.baseEllipse,
				this.dir,
				-this.normalDir)
		},
		toMesh: function (zStart, zEnd) {
			zStart = zStart || -100
			zEnd = zEnd || 100
			var mesh = new GL.Mesh({triangles: true, lines: false, normals: true})
			var pF = this.parametricFunction(), pN = this.parametricNormal()
			var split = 4 * 3, inc = 2 * PI / split
			var c = split * 2
			for (var i = 0; i < split; i++) {
				var v = pF(i * inc, zStart)
				mesh.vertices.push(pF(i * inc, zStart), pF(i * inc, zEnd))
				pushQuad(mesh.triangles, 2 * i, (2 * i + 2) % c, (2 * i + 1), (2 * i + 3) % c)
				var normal = pN(i * inc, 0)
				mesh.normals.push(normal, normal)
			}
			console.log(mesh)
			//mesh.computeNormalLi00nes()
			mesh.compile()
			return mesh
		},
		parametricNormal: function () {
			return (d, z) => {
				return this.baseEllipse.normalAt(d).rejectedFrom(this.dir).toLength(this.normalDir)
			}
		},
		normalAt: function (p) {
			var localP = this.inverseMatrix.transformPoint(p)
			return this.parametricNormal()(localP.angleXY(), localP.z)
		},
		parametricFunction: function () {
			return (d, z) => {
				return this.baseEllipse.at(d).plus(this.dir.times(z))
			}
		},
		implicitFunction: function () {
			return (pWC) => {
				var p = this.inverseMatrix.transformPoint(pWC)
				var radiusLC = p.lengthXY()
				return this.normalDir * (1 - radiusLC)
			}
		},
		containsPoint: function (p) {
			return NLA.isZero(this.implicitFunction()(p))
		},
		boundsFunction: function () {
			assert(false)
		},
		pointToParameterFunction: function () {
			return (pWC, hint) => {
				var p2 = this.inverseMatrix.transformPoint(pWC)
				var angle = p2.angleXY()
				if (angle < -Math.PI + NLA.PRECISION || angle > Math.PI - NLA.PRECISION) {
					angle = hint.dot(this.baseEllipse.f2) < 0
						? Math.PI
						: -Math.PI
				}
				return V3.create(angle, p2.z, 0)
			}
		},
		getIntersectionsWithSurface: function (surface2) {
			if (surface2 instanceof PlaneSurface) {
				return this.getIntersectionsWithPlane(surface2.plane)
			}
		},
		getIntersectionsWithPlane: function (plane) {
			assert(plane instanceof P3)
			if (this.dir.isPerpendicularTo(plane.normal)) {
				var eTs = this.baseEllipse.getIntersectionsWithPlane(plane)
				return eTs.map(t => L3(this.baseEllipse.at(t), this.dir))
			} else {
				return [this.baseEllipse.transform(M4.projection(plane, this.dir))]
			}
		},
		edgeLoopCCW: function (contour) {
			if (contour.length < 3) {
				var totalAngle = 0
				for (var i = 0; i < contour.length; i++) {
					var ipp = (i + 1) % contour.length
					var edge = contour[i], nextEdge = contour[ipp]
					if (edge.curve instanceof EllipseCurve) {
						totalAngle += edge.rotViaPlane(this.plane.normal)
						console.log(edge.toString(), edge.rotViaPlane(this.plane.normal))
					}
					totalAngle += edge.bDir.angleRelativeNormal(nextEdge.aDir, this.plane.normal)
				}
				return totalAngle > 0
			} else {
				var ptpF = this.pointToParameterFunction()
				return isCCW(contour.map(e => ptpF(e.a)), V3.create(0, 0, this.normalDir))
			}
		},
	},
	{
		cyl: function (radius) {
			return new CylinderSurface(new EllipseCurve(V3.ZERO, V3(radius, 0, 0), V3(0, radius, 0)), V3.Z, 1)
		}
	}
)
