'use strict';

Object.defineProperty(exports, '__esModule', { value: true });

function _interopDefault (ex) { return (ex && (typeof ex === 'object') && 'default' in ex) ? ex['default'] : ex; }

var ts3dutils = require('ts3dutils');
var tsgl = require('tsgl');
var opentype = require('opentype.js');
var svgPathdata = require('svg-pathdata');
var javasetmap_ts = require('javasetmap.ts');
var earcut = _interopDefault(require('earcut'));
var nerdamer = _interopDefault(require('nerdamer'));
var chroma = _interopDefault(require('chroma-js'));

const { ceil, floor, abs: abs$1 } = Math;
class Curve extends ts3dutils.Transformable {
    constructor(tMin, tMax) {
        super();
        this.tMin = tMin;
        this.tMax = tMax;
        ts3dutils.assertNumbers(tMin, tMax);
        ts3dutils.assert('number' == typeof tMin && !isNaN(tMin));
        ts3dutils.assert('number' == typeof tMax && !isNaN(tMax));
        ts3dutils.assert(tMin < tMax);
    }
    static integrate(curve, startT, endT, steps) {
        const step = (endT - startT) / steps;
        let length = 0;
        let p = curve.at(startT);
        let i = 0, t = startT + step;
        for (; i < steps; i++, t += step) {
            const next = curve.at(t);
            length += p.distanceTo(next);
            p = next;
        }
        return length;
    }
    static ispsRecursive(curve1, tMin, tMax, curve2, sMin, sMax) {
        // the recursive function finds good approximates for the intersection points
        // curve1 function uses newton iteration to improve the result as much as possible
        function handleStartTS(startT, startS) {
            if (!result.some(info => ts3dutils.eq(info.tThis, startT) && ts3dutils.eq(info.tOther, startS))) {
                const f1 = (t, s) => curve1.tangentAt(t).dot(curve1.at(t).minus(curve2.at(s)));
                const f2 = (t, s) => curve2.tangentAt(s).dot(curve1.at(t).minus(curve2.at(s)));
                // f = (b1, b2, t1, t2) = b1.tangentAt(t1).dot(b1.at(t1).minus(b2.at(t2)))
                const dfdt1 = (b1, b2, t1, t2) => b1.ddt(t1).dot(b1.at(t1).minus(b2.at(t2))) + (b1.tangentAt(t1).squared());
                const dfdt2 = (b1, b2, t1, t2) => -b1.tangentAt(t1).dot(b2.tangentAt(t2));
                const ni = ts3dutils.newtonIterate2dWithDerivatives(f1, f2, startT, startS, 16, dfdt1.bind(undefined, curve1, curve2), dfdt2.bind(undefined, curve1, curve2), (t, s) => -dfdt2(curve2, curve1, s, t), (t, s) => -dfdt1(curve2, curve1, s, t));
                ts3dutils.assert(isFinite(ni.x));
                ts3dutils.assert(isFinite(ni.y));
                if (ni == undefined)
                    console.log(startT, startS, curve1.sce, curve2.sce);
                result.push({ tThis: ni.x, tOther: ni.y, p: curve1.at(ni.x) });
            }
        }
        // returns whether an intersection was immediately found (i.e. without further recursion)
        function findRecursive(tMin, tMax, sMin, sMax, curve1AABB, curve2AABB, depth = 0) {
            const EPS = ts3dutils.NLA_PRECISION;
            if (curve1AABB.fuzzyTouchesAABB(curve2AABB)) {
                const tMid = (tMin + tMax) / 2;
                const sMid = (sMin + sMax) / 2;
                if (Math.abs(tMax - tMin) < EPS || Math.abs(sMax - sMin) < EPS) {
                    handleStartTS(tMid, sMid);
                    return true;
                }
                else {
                    const curve1AABBleft = curve1.getAABB(tMin, tMid);
                    const curve2AABBleft = curve2.getAABB(sMin, sMid);
                    let curve1AABBright, curve2AABBright;
                    // if one of the following calls immediately finds an intersection, we don't want to call the others
                    // as that will lead to the same intersection being output multiple times
                    findRecursive(tMin, tMid, sMin, sMid, curve1AABBleft, curve2AABBleft, depth + 1)
                        || findRecursive(tMin, tMid, sMid, sMax, curve1AABBleft, curve2AABBright = curve2.getAABB(sMid, sMax), depth + 1)
                        || findRecursive(tMid, tMax, sMin, sMid, curve1AABBright = curve1.getAABB(tMid, tMax), curve2AABBleft, depth + 1)
                        || findRecursive(tMid, tMax, sMid, sMax, curve1AABBright, curve2AABBright, depth + 1);
                }
            }
            return false;
        }
        const result = [];
        findRecursive(tMin, tMax, sMin, sMax, curve1.getAABB(tMin, tMax), curve2.getAABB(sMin, sMax));
        return ts3dutils.fuzzyUniquesF(result, info => info.tThis);
    }
    static breakDownIC(implicitCurve, { sMin, sMax, tMin, tMax }, sStep, tStep, stepSize, dids, didt) {
        const bounds = (s, t) => sMin <= s && s <= sMax && tMin <= t && t <= tMax;
        const deltaS = sMax - sMin, deltaT = tMax - tMin;
        const sRes = ceil(deltaS / sStep), tRes = ceil(deltaT / tStep);
        const grid = new Array(sRes * tRes).fill(0);
        ts3dutils.arrayFromFunction(tRes, i => grid.slice(sRes * i, sRes * (i + 1)).map(v => v ? 'X' : '_').join('')).join('\n');
        const at = (i, j) => grid[j * sRes + i];
        const set = (i, j) => 0 <= i && i < sRes && 0 <= j && j < tRes && (grid[j * sRes + i] = 1);
        const result = [];
        const logTable = [];
        for (let i = 0; i < sRes; i++) {
            search: for (let j = 0; j < tRes; j++) {
                if (at(i, j))
                    continue;
                set(i, j);
                let s = sMin + (i + 0.5) * sStep, t = tMin + (j + 0.5) * tStep;
                const startS = s, startT = t;
                // basically curvePoint
                for (let k = 0; k < 8; k++) {
                    const fp = implicitCurve(s, t);
                    const dfpdx = implicitCurve.x(s, t), dfpdy = implicitCurve.y(s, t);
                    if (0 == dfpdx * dfpdx + dfpdy * dfpdy) {
                        // top of a hill, keep looking
                        continue search;
                    }
                    const scale = fp / (dfpdx * dfpdx + dfpdy * dfpdy);
                    s -= scale * dfpdx;
                    t -= scale * dfpdy;
                }
                const li = floor((s - sMin) / sStep), lj = floor((t - tMin) / tStep);
                logTable.push({
                    i,
                    j,
                    li,
                    lj,
                    startS,
                    startT,
                    s,
                    t,
                    'bounds(s, t)': bounds(s, t),
                    'ic(s,t)': implicitCurve(s, t),
                });
                if (!(i == li && j == lj) && at(li, lj)) {
                    continue search;
                }
                set(li, lj);
                // s, t are now good starting coordinates to use follow algo
                if (bounds(s, t) && ts3dutils.eq0(implicitCurve(s, t))) {
                    console.log(ts3dutils.V(s, t).sce);
                    const subresult = mkcurves(implicitCurve, s, t, stepSize, implicitCurve.x, implicitCurve.y, bounds);
                    for (const curvedata of subresult) {
                        ts3dutils.assert(curvedata.points.length > 2);
                        for (const { x, y } of curvedata.points) {
                            const lif = (x - sMin) / sStep, ljf = (y - tMin) / tStep;
                            set((lif - 0.5) | 0, (ljf - 0.5) | 0);
                            set((lif - 0.5) | 0, (ljf + 0.5) | 0);
                            set((lif + 0.5) | 0, (ljf - 0.5) | 0);
                            set((lif + 0.5) | 0, (ljf + 0.5) | 0);
                        }
                    }
                    result.push(...subresult);
                }
            }
        }
        //console.table(logTable)
        for (const { points } of result) {
            for (let i = 0; i < points.length - 1; i++) {
                ts3dutils.assert(!points[i].equals(points[i + 1]));
            }
        }
        return result;
    }
    toString() {
        return this.toSource();
    }
    toSource(rounder = x => x) {
        return ts3dutils.callsce.call(undefined, 'new ' + this.constructor.name, ...this.getConstructorParameters());
    }
    withBounds(tMin = this.tMin, tMax = this.tMax) {
        ts3dutils.assert(this.tMin <= tMin && tMin <= this.tMax);
        ts3dutils.assert(this.tMin <= tMax && tMax <= this.tMax);
        ts3dutils.assert(this.tMin <= tMax && tMax <= this.tMax);
        return new this.constructor(...this.getConstructorParameters().slice(0, -2), tMin, tMax);
    }
    /**
     * The point on the line that is closest to the given point.
     */
    closestPointToPoint(p) {
        return this.at(this.closestTToPoint(p));
    }
    isValidT(t) {
        return ts3dutils.le(this.tMin, t) && ts3dutils.le(t, this.tMax);
    }
    diff(t, eps) {
        return this.at(t).to(this.at(t + eps));
    }
    closestTToPoint(p, tStart) {
        // this.at(t) has minimal distance to p when this.tangentAt(t) is perpendicular to
        // the vector between this.at(t) and p. This is the case iff the dot product of the two is 0.
        // f = (this.at(t) - p) . (this.tangentAt(t)
        // df = this.tangentAt(t) . this.tangentAt(t) + (this.at(t) - p) . this.ddt(t)
        //    = this.tangentAt(t)² + (this.at(t) - p) . this.ddt(t)
        const f = (t) => this.at(t).minus(p).dot(this.tangentAt(t)); // 5th degree polynomial
        const df = (t) => this.tangentAt(t).squared() + (this.at(t).minus(p).dot(this.ddt(t)));
        const STEPS = 32;
        const startT = undefined !== tStart
            ? tStart
            : ts3dutils.arrayFromFunction(STEPS, i => this.tMin + (this.tMax - this.tMin) * i / STEPS)
                .withMax(t => -this.at(t).distanceTo(p));
        return ts3dutils.newtonIterateWithDerivative(f, startT, 16, df);
    }
    /**
     * So different edges on the same curve do not have different vertices, they are always generated
     * on fixed points this.at(k * this.tIncrement), with k taking integer values
     *
     */
    calcSegmentPoints(aT, bT, a, b, reversed, includeFirst) {
        ts3dutils.assert(this.tIncrement, 'tIncrement not defined on ' + this);
        const inc = this.tIncrement;
        const points = [];
        if (includeFirst)
            points.push(a);
        ts3dutils.assert(reversed != aT < bT);
        if (aT < bT) {
            const start = Math.ceil((aT + ts3dutils.NLA_PRECISION) / inc);
            const end = Math.floor((bT - ts3dutils.NLA_PRECISION) / inc);
            for (let i = start; i <= end; i++) {
                points.push(this.at(i * inc));
            }
        }
        else {
            const start = Math.floor((aT - ts3dutils.NLA_PRECISION) / inc);
            const end = Math.ceil((bT + ts3dutils.NLA_PRECISION) / inc);
            for (let i = start; i >= end; i--) {
                points.push(this.at(i * inc));
            }
        }
        points.push(b);
        return points;
    }
    /**
     *
     * @param p
     * @param tStart Defines interval with tEnd in which a start value for t will be searched.
     * Result is not necessarily in this interval.
     * @param tEnd
     */
    distanceToPoint(p, tStart, tEnd) {
        const closestT = this.closestTToPoint(p, tStart, tEnd);
        return this.at(closestT).distanceTo(p);
    }
    asSegmentDistanceToPoint(p, tStart, tEnd) {
        let t = this.closestTToPoint(p, tStart, tEnd);
        t = ts3dutils.clamp(t, tStart, tEnd);
        return this.at(t).distanceTo(p);
    }
    /**
     * Behavior when curves are colinear: self intersections
     */
    isInfosWithCurve(curve) {
        return Curve.ispsRecursive(this, this.tMin, this.tMax, curve, curve.tMin, curve.tMax);
    }
    arcLength(startT, endT, steps = 1) {
        ts3dutils.assert(startT < endT, 'startT < endT');
        return ts3dutils.glqInSteps(t => this.tangentAt(t).length(), startT, endT, steps);
    }
    getAABB(tMin = this.tMin, tMax = this.tMax) {
        tMin = isFinite(tMin) ? tMin : this.tMin;
        tMax = isFinite(tMax) ? tMax : this.tMax;
        const tMinAt = this.at(tMin), tMaxAt = this.at(tMax);
        const roots = this.roots();
        const mins = new Array(3), maxs = new Array(3);
        for (let dim = 0; dim < 3; dim++) {
            const tRoots = roots[dim];
            mins[dim] = Math.min(tMinAt.e(dim), tMaxAt.e(dim));
            maxs[dim] = Math.max(tMinAt.e(dim), tMaxAt.e(dim));
            for (const tRoot of tRoots) {
                if (tMin < tRoot && tRoot < tMax) {
                    mins[dim] = Math.min(mins[dim], this.at(tRoot).e(dim));
                    maxs[dim] = Math.max(maxs[dim], this.at(tRoot).e(dim));
                }
            }
        }
        return new ts3dutils.AABB(ts3dutils.V3.fromArray(mins), ts3dutils.V3.fromArray(maxs));
    }
    reversed() {
        throw new Error();
    }
    clipPlane(plane) {
        const ists = this.isTsWithPlane(plane).filter(ist => this.tMin <= ist && ist <= this.tMax);
        return ts3dutils.getIntervals(ists, this.tMin, this.tMax).mapFilter(([a, b]) => {
            const midT = (a + b) / 2;
            return !ts3dutils.eq(a, b) && plane.distanceToPointSigned(this.at(midT)) < 0 && this.withBounds(a, b);
        });
    }
}
Curve.hlol = 0;
function mkcurves(implicitCurve, sStart, tStart, stepSize, dids, didt, bounds) {
    const start = ts3dutils.V(sStart, tStart);
    // checkDerivate(s => implicitCurve(s, 0), s => dids(s, 0), -1, 1, 0)
    // checkDerivate(t => implicitCurve(0, t), t => didt(0, t), -1, 1, 0)
    const { points, tangents } = followAlgorithm2d(implicitCurve, start, stepSize, bounds);
    if (points[0].distanceTo(points.last) < stepSize && points.length > 2) {
        // this is a loop: split it
        for (let i = 0; i < points.length - 1; i++) {
            ts3dutils.assert(!points[i].equals(points[i + 1]));
        }
        const half = floor(points.length / 2);
        const points1 = points.slice(0, half), points2 = points.slice(half - 1, points.length);
        const tangents1 = tangents.slice(0, half), tangents2 = tangents.slice(half - 1, tangents.length);
        tangents2[tangents2.length - 1] = tangents1[0];
        points2[tangents2.length - 1] = points1[0];
        for (let i = 0; i < points1.length - 1; i++) {
            ts3dutils.assert(!points1[i].equals(points1[i + 1]));
        }
        for (let i = 0; i < points2.length - 1; i++) {
            ts3dutils.assert(!points2[i].equals(points2[i + 1]));
        }
        return [{ points: points1, tangents: tangents1 }, { points: points2, tangents: tangents2 }];
    }
    else {
        // not a loop: check in the other direction
        const { points: reversePoints, tangents: reverseTangents } = followAlgorithm2d(implicitCurve, start, -stepSize, bounds);
        const result = followAlgorithm2d(implicitCurve, reversePoints.last, stepSize, bounds, undefined, reverseTangents.last.negated());
        ts3dutils.assert(result.points.length > 2);
        return [result];
    }
}
function curvePoint(implicitCurve, startPoint, dids, didt) {
    let p = startPoint;
    for (let i = 0; i < 8; i++) {
        const fp = implicitCurve(p.x, p.y);
        const dfpdx = dids(p.x, p.y), dfpdy = didt(p.x, p.y);
        const scale = fp / (dfpdx * dfpdx + dfpdy * dfpdy);
        //console.log(p.$)
        p = p.minus(new ts3dutils.V3(scale * dfpdx, scale * dfpdy, 0));
    }
    return p;
}
function curvePointMF(mf, startPoint, steps = 8, eps = 1 / (1 << 30)) {
    let p = startPoint;
    for (let i = 0; i < steps; i++) {
        const fp = mf(p.x, p.y);
        const dfpdx = mf.x(p.x, p.y), dfpdy = mf.y(p.x, p.y);
        const scale = fp / (dfpdx * dfpdx + dfpdy * dfpdy);
        //console.log(p.$)
        p = p.minus(new ts3dutils.V3(scale * dfpdx, scale * dfpdy, 0));
        if (abs$1(fp) <= eps)
            break;
    }
    return p;
}

const { PI } = Math;
class XiEtaCurve extends Curve {
    constructor(center, f1, f2, tMin = -PI, tMax = PI) {
        super(tMin, tMax);
        this.center = center;
        this.f1 = f1;
        this.f2 = f2;
        this.tMin = tMin;
        this.tMax = tMax;
        ts3dutils.assertVectors(center, f1, f2);
        this.normal = f1.cross(f2);
        if (!this.normal.likeO()) {
            this.normal = this.normal.unit();
            this.matrix = ts3dutils.M4.forSys(f1, f2, this.normal, center);
            this.inverseMatrix = this.matrix.inversed();
        }
        else {
            this.matrix = ts3dutils.M4.forSys(f1, f2, f1.unit(), center);
            const f1p = f1.getPerpendicular();
            this.inverseMatrix = new ts3dutils.M4(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1).times(ts3dutils.M4.forSys(f1, f1p, f1.cross(f1p), center).inversed());
        }
    }
    static magic(a, b, c) {
        throw new Error('abstract');
    }
    /**
     * Returns a new EllipseCurve representing an ellipse parallel to the XY-plane
     * with semi-major/minor axes parallel t the X and Y axes and of length a and b.
     *
     * @param a length of the axis parallel to X axis
     * @param b length of the axis parallel to Y axis
     * @param center Defaults to V3.O
     */
    static forAB(a, b, center = ts3dutils.V3.O) {
        return new this(center, ts3dutils.V(a, 0, 0), ts3dutils.V(0, b, 0));
    }
    static XYLCValid(pLC) {
        throw new Error('abstract');
    }
    static XYLCPointT(pLC) {
        throw new Error('abstract');
    }
    static unitIsInfosWithLine(anchorLC, dirLC, anchorWC, dirWC) {
        throw new Error('abstract');
    }
    addToMesh(mesh, res = 4, radius = 0, pointStep = 1) {
        const baseNormals = ts3dutils.arrayFromFunction(res, i => ts3dutils.V3.polar(1, ts3dutils.TAU * i / res));
        const baseVertices = ts3dutils.arrayFromFunction(res, i => ts3dutils.V3.polar(radius, ts3dutils.TAU * i / res));
        const inc = this.tIncrement;
        const start = Math.ceil((this.tMin + ts3dutils.NLA_PRECISION) / inc);
        const end = Math.floor((this.tMax - ts3dutils.NLA_PRECISION) / inc);
        for (let i = start; i <= end; i += pointStep) {
            const t = i * inc;
            const start = mesh.vertices.length;
            if (0 !== i) {
                for (let j = 0; j < res; j++) {
                    tsgl.pushQuad(mesh.TRIANGLES, true, start - res + j, start + j, start - res + (j + 1) % res, start + (j + 1) % res);
                }
            }
            const point = this.at(t), tangent = this.tangentAt(t);
            const matrix = ts3dutils.M4.forSys(this.normal, tangent.cross(this.normal), tangent, point);
            mesh.normals.push(...matrix.transformedVectors(baseNormals));
            mesh.vertices.push(...matrix.transformedPoints(baseVertices));
        }
    }
    getConstructorParameters() {
        return [this.center, this.f1, this.f2, this.tMin, this.tMax];
    }
    isInfosWithCurve(curve) {
        if (curve instanceof L3$1) {
            return this.isInfosWithLine(curve.anchor, curve.dir1, this.tMin, this.tMax, curve.tMin, curve.tMax);
        }
        if (curve instanceof BezierCurve) {
            return this.isInfosWithBezier(curve);
        }
        if (curve instanceof XiEtaCurve) {
            if (!this.normal.isParallelTo(curve.normal)) {
                return this.isTsWithPlane(curve.getPlane()).mapFilter(tThis => {
                    const p = this.at(tThis);
                    if (curve.containsPoint(p)) {
                        return { tThis, tOther: curve.pointT(p), p };
                    }
                });
            }
        }
        return super.isInfosWithCurve(curve);
    }
    transform(m4) {
        return new this.constructor(m4.transformPoint(this.center), m4.transformVector(this.f1), m4.transformVector(this.f2), this.tMin, this.tMax);
    }
    equals(obj) {
        return this == obj ||
            obj.constructor == this.constructor
                && this.center.equals(obj.center)
                && this.f1.equals(obj.f1)
                && this.f2.equals(obj.f2);
    }
    hashCode() {
        let hashCode = 0;
        hashCode = hashCode * 31 + this.center.hashCode();
        hashCode = hashCode * 31 + this.f1.hashCode();
        hashCode = hashCode * 31 + this.f2.hashCode();
        return hashCode | 0;
    }
    likeCurve(curve) {
        return ts3dutils.hasConstructor(curve, this.constructor)
            && this.center.like(curve.center)
            && this.f1.like(curve.f1)
            && this.f2.like(curve.f2);
    }
    normalP(t) {
        return this.tangentAt(t).cross(this.normal);
    }
    getPlane() {
        return P3.normalOnAnchor(this.normal, this.center);
    }
    isTsWithPlane(plane) {
        ts3dutils.assertInst(P3, plane);
        /*
         this: x = center + f1 * cos t + f2 * sin t  (1)
         plane:
         n := plane.normal1
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
         g1 * xi + g2 * eta = g3 (6)
         */
        if (plane.normal1.isParallelTo(this.normal)) {
            return [];
        }
        const n = plane.normal1, w = plane.w, center = this.center, f1 = this.f1, f2 = this.f2, g1 = n.dot(f1), g2 = n.dot(f2), g3 = w - n.dot(center);
        return this.constructor.magic(g1, g2, g3);
    }
    pointT(p) {
        ts3dutils.assertVectors(p);
        const pLC = this.inverseMatrix.transformPoint(p);
        return this.constructor.XYLCPointT(pLC);
    }
    containsPoint(p) {
        const pLC = this.inverseMatrix.transformPoint(p);
        return ts3dutils.eq0(pLC.z) && this.constructor.XYLCValid(pLC);
    }
    isInfosWithLine(anchorWC, dirWC, tMin, tMax, lineMin = -100000, lineMax = 100000) {
        const anchorLC = this.inverseMatrix.transformPoint(anchorWC);
        const dirLC = this.inverseMatrix.transformVector(dirWC);
        if (ts3dutils.eq0(dirLC.z)) {
            // local line parallel to XY-plane
            if (ts3dutils.eq0(anchorLC.z)) {
                // local line lies in XY-plane
                return this.constructor.unitIsInfosWithLine(anchorLC, dirLC, anchorWC, dirWC);
            }
        }
        else {
            // if the line intersects the XY-plane in a single point, there can be an intersection there
            // find point, then check if distance from circle = 1
            const otherTAtZ0 = anchorLC.z / dirLC.z;
            const isp = dirLC.times(otherTAtZ0).plus(anchorLC);
            if (this.constructor.XYLCValid(isp)) {
                // point lies on unit circle
                return [{
                        tThis: this.constructor.XYLCPointT(isp),
                        tOther: otherTAtZ0,
                        p: anchorWC.plus(dirWC.times(otherTAtZ0)),
                    }];
            }
        }
        return [];
    }
    isTsWithSurface(surface) {
        if (surface instanceof PlaneSurface$1) {
            return this.isTsWithPlane(surface.plane);
        }
        else if (surface instanceof SemiEllipsoidSurface) {
            const isEllipse = surface.asEllipsoidSurface().isCurvesWithSurface(new PlaneSurface$1(this.getPlane()));
            if (isEllipse.length < 1)
                return [];
            const possibleInfos = this.isInfosWithCurve(isEllipse[0]);
            return possibleInfos.filter(info => surface.containsPoint(info.p)).map(info => info.tThis);
        }
        else if (surface instanceof ProjectedCurveSurface ||
            surface instanceof EllipsoidSurface ||
            surface instanceof ConicSurface) {
            return surface.isCurvesWithPlane(this.getPlane())
                .flatMap(curve => this.isInfosWithCurve(curve))
                .map(info => info.tThis);
        }
        else {
            throw new Error();
        }
    }
    isInfosWithBezier(bezierWC) {
        const bezierLC = bezierWC.transform(this.inverseMatrix);
        if (new PlaneSurface$1(P3.XY).containsCurve(bezierLC)) {
            return this.isInfosWithBezier2D(bezierWC);
        }
        else {
            const infos = bezierLC.isTsWithPlane(P3.XY).mapFilter(tOther => {
                const pLC = bezierLC.at(tOther);
                if (this.constructor.XYLCValid(pLC)) {
                    return { tOther: tOther, p: bezierWC.at(tOther), tThis: this.constructor.XYLCPointT(pLC) };
                }
            });
            return infos;
        }
    }
    isInfosWithBezier2D(bezierWC, sMin, sMax) {
        sMin = isFinite(sMin) ? sMin : bezierWC.tMin;
        sMax = isFinite(sMax) ? sMax : bezierWC.tMax;
        ts3dutils.assertf(() => 0 < Math.PI);
        ts3dutils.assertf(() => sMin < sMax);
        return Curve.ispsRecursive(this, this.tMin, this.tMax, bezierWC, sMin, sMax);
    }
    isOrthogonal() {
        return this.f1.isPerpendicularTo(this.f2);
    }
    at2(xi, eta) {
        ts3dutils.assertNumbers(xi, eta);
        // center + f1 xi + f2 eta
        return this.center.plus(this.f1.times(xi)).plus(this.f2.times(eta));
    }
    debugToMesh(mesh, bufferName) {
        mesh[bufferName] || mesh.addVertexBuffer(bufferName, bufferName);
        for (let t = 0; t < Math.PI; t += 0.1) {
            const p = this.at(t);
            mesh[bufferName].push(p, p.plus(this.tangentAt(t).toLength(1)));
            mesh[bufferName].push(p, p.plus(this.normalP(t).toLength(1)));
        }
        mesh[bufferName].push(this.center, this.center.plus(this.f1.times(1.2)));
        mesh[bufferName].push(this.center, this.center.plus(this.f2));
        mesh[bufferName].push(this.center, this.center.plus(this.normal));
    }
}

const { ceil: ceil$1, floor: floor$1 } = Math;
class ImplicitCurve extends Curve {
    constructor(points, tangents, dir = 1, generator, tMin = (1 == dir ? 0 : -(points.length - 1)), tMax = (1 == dir ? points.length - 1 : 0)) {
        super(tMin, tMax);
        this.points = points;
        this.tangents = tangents;
        this.dir = dir;
        this.generator = generator;
        ts3dutils.assert(points.length > 2);
        ts3dutils.assert(0 <= tMin && tMin <= points.length - 1);
        ts3dutils.assert(0 <= tMax && tMax <= points.length - 1);
    }
    likeCurve(curve) {
        throw new Error('Method not implemented.');
    }
    toSource(rounder = x => x) {
        return this.generator || super.toSource(rounder);
    }
    containsPoint(p) {
        ts3dutils.assertVectors(p);
        return !isNaN(this.pointT(p));
    }
    equals(obj) {
        return this == obj ||
            Object.getPrototypeOf(obj) == PICurve.prototype
                && this.points[0].equals(obj.points[0])
                && this.tangents[0].equals(obj.tangents[0]);
    }
    hashCode() {
        return [this.points[0], this.tangents[0]].hashCode();
    }
    tangentP(pWC) {
        ts3dutils.assertVectors(pWC);
        ts3dutils.assert(this.containsPoint(pWC), 'this.containsPoint(pWC)' + this.containsPoint(pWC));
        const t = this.pointT(pWC);
        return this.tangentAt(t);
    }
    tangentAt(t) {
        t = ts3dutils.clamp(t, this.tMin, this.tMax);
        return ts3dutils.V3.lerp(this.tangents[floor$1(t)], this.tangents[ceil$1(t)], t % 1);
    }
    at(t) {
        ts3dutils.assert(!isNaN(t));
        return ts3dutils.V3.lerp(this.points[floor$1(t)], this.points[ceil$1(t)], t % 1);
    }
    getConstructorParameters() {
        return [];
    }
    transform(m4) {
        return new ImplicitCurve(m4.transformedPoints(this.points), m4.transformedVectors(this.tangents));
    }
    roots() {
        const allTs = ts3dutils.arrayRange(0, this.points.length);
        return [allTs, allTs, allTs];
    }
    addToMesh(mesh, res = 4, radius = 0, pointStep = 1) {
        const baseNormals = ts3dutils.arrayFromFunction(res, i => ts3dutils.V3.polar(1, ts3dutils.TAU * i / res));
        const baseVertices = ts3dutils.arrayFromFunction(res, i => ts3dutils.V3.polar(radius, ts3dutils.TAU * i / res));
        let prevTangent = ts3dutils.V3.Z, prevMatrix = ts3dutils.M4.IDENTITY;
        for (let i = ceil$1(this.tMin); i < floor$1(this.tMax); i += pointStep) {
            const start = mesh.vertices.length;
            if (ceil$1(this.tMin) !== i) {
                for (let j = 0; j < res; j++) {
                    tsgl.pushQuad(mesh.TRIANGLES, true, start - res + j, start + j, start - res + (j + 1) % res, start + (j + 1) % res);
                }
            }
            const point = this.points[i], tangent = this.tangents[i];
            const tangentMatrix = ts3dutils.M4.rotateAB(prevTangent, tangent).times(prevMatrix);
            mesh.normals.push(...tangentMatrix.transformedVectors(baseNormals));
            const baseMatrix = ts3dutils.M4.translate(point).times(tangentMatrix);
            mesh.vertices.push(...baseMatrix.transformedPoints(baseVertices));
            prevTangent = tangent;
            prevMatrix = tangentMatrix;
        }
    }
}
ImplicitCurve.prototype.tIncrement = 1;

const { PI: PI$1, abs: abs$2, sin, cos } = Math;
class BezierCurve extends Curve {
    constructor(p0, p1, p2, p3, tMin = -0.1, tMax = 1.1) {
        super(tMin, tMax);
        ts3dutils.assertVectors(p0, p1, p2, p3);
        ts3dutils.assert(isFinite(tMin) && isFinite(tMax));
        //assert(!L3.throughPoints(p0, p3).containsPoint(p1) || !L3.throughPoints(p0, p3).containsPoint(p2))
        this.p0 = p0;
        this.p1 = p1;
        this.p2 = p2;
        this.p3 = p3;
    }
    get points() {
        return [this.p0, this.p1, this.p2, this.p3];
    }
    /**
     * Returns a curve with curve.at(x) == V(x, ax³ + bx² + cx + d, 0)
     */
    static graphXY(a, b, c, d, tMin, tMax) {
        // d = p0y
        // c = -3 p0y + 3 p1y => p1y = c/3 + p0y
        // b = 3 p0y - 6 p1y + 3 p2y => p2y = b/3 - p0y + 2 p1y
        // a = -p0y + 3 p1y -3 p2y + p3y => p3y = a + p0y - 3 p1y + 3 p2y
        const p0y = d;
        const p1y = c / 3 + p0y;
        const p2y = b / 3 - p0y + 2 * p1y;
        const p3y = a + p0y - 3 * p1y + 3 * p2y;
        return new BezierCurve(ts3dutils.V(0, p0y), ts3dutils.V(1 / 3, p1y), ts3dutils.V(2 / 3, p2y), ts3dutils.V(1, p3y), tMin, tMax);
    }
    static quadratic(a, b, c, tMin = 0, tMax = 1) {
        const line = L3$1.throughPoints(a, c);
        if (line.containsPoint(b)) {
            return line;
        }
        else {
            // p1 = 1/3 a + 2/3 b
            // p2 = 1/3 c + 2/3 b
            return new BezierCurve(a, b.times(2).plus(a).div(3), b.times(2).plus(c).div(3), c, tMin, tMax);
        }
    }
    /**
     * Returns a bezier curve which approximates a CCW unit circle arc starting at V3.X of angle phi
     * phi <= PI / 2 is recommended
     *
     * Formula from here: https://pomax.github.io/bezierinfo/#circles_cubic
     */
    static approximateUnitArc(phi) {
        const f = 4 / 3 * Math.tan(phi / 4);
        return new BezierCurve(ts3dutils.V3.X, new ts3dutils.V3(1, f, 0), new ts3dutils.V3(cos(phi) + f * sin(phi), sin(phi) - f * cos(phi), 0), ts3dutils.V3.sphere(phi, 0), 0, 1);
    }
    static testEdges() {
        const curve2 = BezierCurve.graphXY(2, -3, -3, 2, 0.6, 2);
        const items = curve2.magic().map(c => Edge.forCurveAndTs(c).translate(3));
        console.log(items.length);
        return [Edge.forCurveAndTs(curve2)].concat(items);
    }
    getConstructorParameters() {
        return [this.p0, this.p1, this.p2, this.p3, this.tMin, this.tMax];
    }
    at(t) {
        // = s^3 p0 + 3 s^2 t p1 + 3 s t^2 p2 + t^3 p3
        ts3dutils.assertNumbers(t);
        const p0 = this.p0, p1 = this.p1, p2 = this.p2, p3 = this.p3;
        const s = 1 - t, c0 = s * s * s, c1 = 3 * s * s * t, c2 = 3 * s * t * t, c3 = t * t * t;
        return new ts3dutils.V3(p0.x * c0 + p1.x * c1 + p2.x * c2 + p3.x * c3, p0.y * c0 + p1.y * c1 + p2.y * c2 + p3.y * c3, p0.z * c0 + p1.z * c1 + p2.z * c2 + p3.z * c3);
    }
    /**
     * s := (1 - t)
     * at(t) := s³ p0 + 3 s² t p1 + 3 s t² p2 + t³ p3
     * tangent(t) := 3 s² (p1 - p0) + 6 s t (p2 - p1) + 3 t² (p3 - p2)
     *            := 3 (1 - t)² (p1 - p0) + 6 (1 - t) t (p2 - p1) + 3 t² (p3 - p2)
     *            := 3 (1 - 2 t + t²) (p1 - p0) + 6 (t - t²) (p2 - p1) + 3 t² (p3 - p2)
     *            := (3 (p3 - p2) - 6 (p2 - p1) + 3 (p1 - p0)) t²*
     *                + (-6 (p1 - p0) + (p2 - p1)) t
     *                + 3 (p1 - p0)
     */
    tangentAt(t) {
        ts3dutils.assertNumbers(t);
        const p0 = this.p0, p1 = this.p1, p2 = this.p2, p3 = this.p3;
        const s = 1 - t, c01 = 3 * s * s, c12 = 6 * s * t, c23 = 3 * t * t;
        return new ts3dutils.V3((p1.x - p0.x) * c01 + (p2.x - p1.x) * c12 + (p3.x - p2.x) * c23, (p1.y - p0.y) * c01 + (p2.y - p1.y) * c12 + (p3.y - p2.y) * c23, (p1.z - p0.z) * c01 + (p2.z - p1.z) * c12 + (p3.z - p2.z) * c23);
    }
    ddt(t) {
        ts3dutils.assertNumbers(t);
        const p0 = this.p0, p1 = this.p1, p2 = this.p2, p3 = this.p3;
        const c012 = 6 * (1 - t), c123 = 6 * t;
        return new ts3dutils.V3((p2.x - 2 * p1.x + p0.x) * c012 + (p3.x - 2 * p2.x + p1.x) * c123, (p2.y - 2 * p1.y + p0.y) * c012 + (p3.y - 2 * p2.y + p1.y) * c123, (p2.z - 2 * p1.z + p0.z) * c012 + (p3.z - 2 * p2.z + p1.z) * c123);
    }
    normalP(t) {
        const tangent = this.tangentAt(t);
        const rot = tangent.cross(this.ddt(t));
        return rot.cross(tangent);
    }
    isTsWithPlane(plane) {
        ts3dutils.assertInst(P3, plane);
        /*
         We are solving for t:
         n := plane.normal1
         this.at(t) DOT n == plane.w // according to plane definition
         (a t³ + b t² + c t + d) DOT n == plane.w // bezier curve as cubic equation
         (a DOT n) t³ + (b DOT n) t³ + (c DOT n) t + d DOT n - plane.w == 0 // multiply out DOT n, minus plane.w
         */
        const { p0, p1, p2, p3 } = this;
        const n = plane.normal1;
        const a = p1.minus(p2).times(3).minus(p0).plus(p3);
        const b = p0.plus(p2).times(3).minus(p1.times(6));
        const c = p1.minus(p0).times(3);
        const d = p0;
        return ts3dutils.solveCubicReal2(a.dot(n), b.dot(n), c.dot(n), d.dot(n) - plane.w)
            .filter(t => ts3dutils.between(t, this.tMin, this.tMax));
    }
    isTsWithSurface(surface) {
        if (surface instanceof PlaneSurface$1) {
            return this.isTsWithPlane(surface.plane);
        }
        if (surface instanceof SemiCylinderSurface) {
            const projPlane = new P3(surface.dir.unit(), 0);
            const projThis = this.project(projPlane);
            const projEllipse = surface.baseCurve.project(projPlane);
            return projEllipse.isInfosWithBezier2D(projThis).map(info => info.tOther);
        }
        if (surface instanceof ProjectedCurveSurface) {
            const projPlane = new P3(surface.dir.unit(), 0);
            const projThis = this.project(projPlane);
            const projEllipse = surface.baseCurve.project(projPlane);
            return projEllipse.isInfosWithCurve(projThis).map(info => info.tOther);
        }
        if (surface instanceof EllipsoidSurface) {
            const thisOC = this.transform(surface.inverseMatrix);
            const f = (t) => thisOC.at(t).length() - 1;
            const df = (t) => thisOC.at(t).unit().dot(thisOC.tangentAt(t));
            const stepSize = 1 / (1 << 11);
            const STEPS = (this.tMax - this.tMin) / stepSize;
            const result = [];
            for (let startT = this.tMin; startT <= this.tMax; startT += stepSize) {
                const dt = stepSize * thisOC.tangentAt(startT).length();
                if (abs$2(f(startT)) <= dt) {
                    //const t = newtonIterate1d(f, startT, 16)
                    let t = ts3dutils.newtonIterateWithDerivative(f, startT, 16, df);
                    if (!ts3dutils.eq0(f(t)) || ts3dutils.eq0(df(t))) {
                        t = ts3dutils.newtonIterate1d(df, startT, 16);
                        //if (f(a) * f(b) < 0) {
                        //    t = bisect(f, a, b, 16)
                        //} else if (df(a) * df(b) < 0) {
                        //    t = bisect(df, a, b, 16)
                        //}
                    }
                    if (ts3dutils.eq0(f(t)) && !result.some(r => ts3dutils.eq(r, t))) {
                        result.push(t);
                    }
                }
            }
            return result;
        }
        if (surface instanceof SemiEllipsoidSurface) {
            return this.isTsWithSurface(surface.asEllipsoidSurface()).filter(t => surface.containsPoint(this.at(t)));
        }
        throw new Error();
    }
    likeCurve(curve) {
        return this == curve ||
            ts3dutils.hasConstructor(curve, BezierCurve)
                && this.p0.like(curve.p0)
                && this.p1.like(curve.p1)
                && this.p2.like(curve.p2)
                && this.p3.like(curve.p3);
    }
    equals(obj) {
        return this == obj ||
            ts3dutils.hasConstructor(obj, BezierCurve)
                && this.p0.equals(obj.p0)
                && this.p1.equals(obj.p1)
                && this.p2.equals(obj.p2)
                && this.p3.equals(obj.p3);
    }
    hashCode() {
        let hashCode = 0;
        hashCode = hashCode * 31 + this.p0.hashCode();
        hashCode = hashCode * 31 + this.p1.hashCode();
        hashCode = hashCode * 31 + this.p2.hashCode();
        hashCode = hashCode * 31 + this.p3.hashCode();
        return hashCode | 0;
    }
    /**
     * Checks if this curve is colinear to the passed curve, i.e.
     * for every t:number there exists a s:number with this.at(t) = curve.at(s)
     */
    isColinearTo(curve) {
        if (this === curve || this.likeCurve(curve))
            return true;
        if (!(curve instanceof BezierCurve))
            return false;
        // first, find out where/if curve.p0 and curve.p3 are on this
        // then split this at curve.p0 --> curve.p3 to compare points p1 and p2
        let curveP0T, curveP3T;
        // assign in if condition to exploit short-circuit
        if (isNaN(curveP0T = this.pointT(curve.p0)) || isNaN(curveP3T = this.pointT(curve.p3))) {
            return false;
        }
        let thisSplit;
        if (ts3dutils.eq(1, curveP0T)) {
            // this.split(curveP0T).right is degenerate in this case, so we need to handle it separately
            // this.split(curveP3T): 0 --> curveP3T --> 1
            // .right: curveP3T --> 1
            // .reversed(): 1 --> curveP3T
            thisSplit = this.split(curveP3T)[1].reversed();
        }
        else {
            // curveP3T describes the point on this
            // adjust it so it describes the same point on this.split(curveP0T).right
            // this:                       0           p0t        p3t      1
            //                             |            |          |       |
            // this.split(curveP0T).right:              0        p3tad     1
            const curveP3Tadjusted = (curveP3T - curveP0T) / (1 - curveP0T);
            thisSplit = this.split(curveP0T)[1].split(curveP3Tadjusted)[0];
        }
        return curve.likeCurve(thisSplit);
    }
    reversed() {
        return new BezierCurve(this.p3, this.p2, this.p1, this.p0, 1 - this.tMax, 1 - this.tMin);
    }
    getCoefficients() {
        const { p0, p1, p2, p3 } = this;
        // calculate cubic equation coefficients
        // a t³ + b t² + c t + d = 0
        // multiplying out the cubic Bézier curve equation gives:
        // a = -p0 + 3 p1 - 3 p2 + p3
        // b = 3 p0 - 6 p1 + 3 p2
        // c = -3 p0 + 3 p1
        // d = p0 - p
        const a = p1.minus(p2).times(3).minus(p0).plus(p3);
        const b = p0.plus(p2).times(3).minus(p1.times(6));
        const c = p1.minus(p0).times(3);
        const d = p0;
        return [a, b, c, d];
    }
    tangentCoefficients() {
        const { p0, p1, p2, p3 } = this;
        const p01 = p1.minus(p0), p12 = p2.minus(p1), p23 = p3.minus(p2);
        const a = p01.plus(p23).times(3).minus(p12.times(6));
        const b = p12.minus(p01).times(6);
        const c = p01.times(3);
        return [ts3dutils.V3.O, a, b, c];
    }
    pointT(p) {
        return this.closestTToPoint(p);
    }
    pointT3(p) {
        const { p0, p1, p2, p3 } = this;
        // calculate cubic equation coefficients
        // a t³ + b t² + c t + d = 0
        // multiplying out the cubic Bézier curve equation gives:
        // a = -p0 + 3 p1 - 3 p2 + p3
        // b = 3 p0 - 6 p1 + 3 p2
        // c = -3 p0 + 3 p1
        // d = p0 - p
        const a = p1.minus(p2).times(3).minus(p0).plus(p3);
        const b = p0.plus(p2).times(3).minus(p1.times(6));
        const c = p1.minus(p0).times(3);
        const d = p0.minus(p);
        // a t³ + b t² + c t + d = 0 is 3 cubic equations, some of which can be degenerate
        const maxDim = ts3dutils.NLA_PRECISION < a.maxAbsElement() ? a.maxAbsDim()
            : ts3dutils.NLA_PRECISION < b.maxAbsElement() ? b.maxAbsDim()
                : ts3dutils.NLA_PRECISION < c.maxAbsElement() ? c.maxAbsDim()
                    : ts3dutils.assertNever();
        const results = ts3dutils.solveCubicReal2(a.e(maxDim), b.e(maxDim), c.e(maxDim), d.e(maxDim)).filter(t => this.at(t).like(p));
        if (0 == results.length)
            return NaN;
        if (1 == results.length)
            return results[0];
        ts3dutils.assert(false, 'multiple intersection ' + this.toString() + p.sce);
    }
    pointT2(p) {
        const { p0, p1, p2, p3 } = this;
        // calculate cubic equation coefficients
        // a t³ + b t² + c t + d = 0
        // multiplying out the cubic Bézier curve equation gives:
        // a = -p0 + 3 p1 - 3 p2 + p3
        // b = 3 p0 - 6 p1 + 3 p2
        // c = -3 p0 + 3 p1
        // d = p0 - p
        const a = p1.minus(p2).times(3).minus(p0).plus(p3).els();
        const b = p0.plus(p2).times(3).minus(p1.times(6)).els();
        const c = p1.minus(p0).times(3).els();
        const d = p0.minus(p).els();
        let results = undefined;
        // assume passed point is on curve and that curve does not self-intersect,
        // i.e. there is exactly one correct result for t
        // try to find a single result in the x-dimension, if multiple are found,
        // filter them by checking the other dimensions
        for (let dim = 0; dim < 3; dim++) {
            if (ts3dutils.eq0(a[dim]) && ts3dutils.eq0(b[dim]) && ts3dutils.eq0(c[dim])) {
                // for case x:
                // ax == bx == cx == 0 => x(t) = dx
                // x value is constant
                // if x == 0 for all t, this does not limit the result, otherwise, there is no result, i.e
                // the passed point is not on the curve
                if (!ts3dutils.eq0(d[dim]))
                    return NaN;
            }
            else {
                const newResults = ts3dutils.solveCubicReal2(a[dim], b[dim], c[dim], d[dim]);
                if (0 == newResults.length)
                    return NaN;
                if (1 == newResults.length)
                    return newResults[0];
                if (results) {
                    results = results.filter(t => newResults.some(t2 => ts3dutils.eq(t, t2)));
                    if (0 == results.length)
                        return NaN;
                    if (1 == results.length)
                        return results[0];
                }
                else {
                    results = newResults;
                }
            }
        }
        ts3dutils.assert(false, 'multiple intersection ' + results + this.toString() + p.sce);
    }
    transform(m4) {
        return new BezierCurve(m4.transformPoint(this.p0), m4.transformPoint(this.p1), m4.transformPoint(this.p2), m4.transformPoint(this.p3), this.tMin, this.tMax);
    }
    isClosed() {
        return this.p0.like(this.p3);
    }
    isQuadratic() {
        return this.p1.like(this.p2);
    }
    debugToMesh(mesh, bufferName) {
        const result = mesh.addVertexBuffer(bufferName, bufferName);
        for (let t = -2; t <= 2; t += 0.01) {
            const p = this.at(t);
            result[bufferName].push(p, p.plus(this.tangentAt(t).toLength(1)));
            result[bufferName].push(p, p.plus(this.normalP(t).toLength(1)));
        }
        result[bufferName].push(this.p0, this.p1);
        result[bufferName].push(this.p1, this.p2);
        result[bufferName].push(this.p2, this.p3);
    }
    split(t) {
        // do de Casteljau's algorithm at t, the resulting points are the points needed to create 2 new curves
        const s = (1 - t);
        const { p0, p1, p2, p3 } = this;
        /*
        p3 // n3
        b01 = s p0 + t p1
        b11 = s p1 + t p2
        b21 = s p2 + t p3 // n2
        b02 = s b01 + t b11
        b12 = s b11 + t b21 // n1
        b03 = s b02 + t b12 // n0

        c01 =
        */
        const b01 = p0.times(s).plus(p1.times(t)), b11 = p1.times(s).plus(p2.times(t)), b21 = p2.times(s).plus(p3.times(t));
        const b02 = b01.times(s).plus(b11.times(t)), b12 = b11.times(s).plus(b21.times(t));
        const b03 = b02.times(s).plus(b12.times(t));
        return [new BezierCurve(p0, b01, b02, b03), new BezierCurve(b03, b12, b21, p3)];
    }
    containsPoint(p) {
        return isFinite(this.pointT(p));
    }
    roots() {
        /**
         *            := (3 (p3 - p2) - 6 (p2 - p1) + 3 (p1 - p0)) t²*
         *                + (-6 (p1 - p0) + 6 (p2 - p1)) t
         *                + 3 (p1 - p0)
         *                */
        const { p0, p1, p2, p3 } = this;
        const p01 = p1.minus(p0), p12 = p2.minus(p1), p23 = p3.minus(p2);
        const a = p01.plus(p23).times(3).minus(p12.times(6));
        const b = p12.minus(p01).times(6);
        const c = p01.times(3);
        return ts3dutils.arrayFromFunction(3, dim => ts3dutils.solveCubicReal2(0, a.e(dim), b.e(dim), c.e(dim)));
    }
    isInfosWithLine(anchorWC, dirWC, tMin, tMax, lineMin = -100000, lineMax = 100000) {
        const dirLength = dirWC.length();
        // TODO: no:
        let result = Curve.ispsRecursive(this, this.tMin, this.tMax, new L3$1(anchorWC, dirWC.unit()), lineMin, lineMax);
        result = ts3dutils.fuzzyUniquesF(result, info => info.tOther);
        result.forEach(info => (info.tOther /= dirLength));
        return result;
        // looking for this.at(t) == line.at(s)
        // this.at(t).x == anchorWC.x + dirWC.x * s
        // (this.at(t).x - anchorWC.x) / dirWC.x == s (analogue for y and z) (1x, 1y, 1z)
        // (1x) - (1y):
        // (this.at(t).x - anchorWC.x) / dirWC.x - (this.at(t).y - anchorWC.y) / dirWC.y == 0
        // (this.at(t).x - anchorWC.x) * dirWC.y - (this.at(t).y - anchorWC.y) * dirWC.x == 0 (2)
        // cubic equation params (see #pointT):
        const { p0, p1, p2, p3 } = this;
        const a = p1.minus(p2).times(3).minus(p0).plus(p3);
        const b = p0.plus(p2).times(3).minus(p1.times(6));
        const c = p1.minus(p0).times(3);
        const d = p0;
        // modifier cubic equation stP to get (1)
        // const w = a.x * dirWC.y - a.y * dirWC.x
        // const x = b.x * dirWC.y - b.y * dirWC.x
        // const y = c.x * dirWC.y - c.y * dirWC.x
        // const z = (d.x - anchorWC.x) * dirWC.y - (d.y - anchorWC.y) * dirWC.x
        // the above version doesn't work for dirWC.x == dirWC.y == 0, so:
        const absMinDim = dirWC.minAbsDim();
        const [coord0, coord1] = [[1, 2], [2, 0], [0, 1]][absMinDim];
        const w = a.e(coord0) * dirWC.e(coord1) - a.e(coord1) * dirWC.e(coord0);
        const x = b.e(coord0) * dirWC.e(coord1) - b.e(coord1) * dirWC.e(coord0);
        const y = c.e(coord0) * dirWC.e(coord1) - c.e(coord1) * dirWC.e(coord0);
        const z = (d.e(coord0) - anchorWC.e(coord0)) * dirWC.e(coord1) - (d.e(coord1) - anchorWC.e(coord1)) * dirWC.e(coord0);
        tMin = isFinite(tMin) ? tMin : this.tMin;
        tMax = isFinite(tMax) ? tMax : this.tMax;
        // we ignored a dimension in the previous step, so we need to check it too
        return ts3dutils.solveCubicReal2(w, x, y, z).mapFilter(tThis => {
            if (tMin <= tThis && tThis <= tMax) {
                const p = this.at(tThis);
                // console.log(t*t*t*w+t*t*x+t*y+z, dirWC.length())
                const s = p.minus(anchorWC).dot(dirWC) / dirWC.dot(dirWC);
                const lineAtS = dirWC.times(s).plus(anchorWC);
                if (lineAtS.like(p))
                    return { tThis: tThis, tOther: s, p: p };
            }
        });
    }
    closestPointToLine(line, tMin, tMax) {
        // (this(t)-line(s)) * line.dir == 0 (1)
        // (this(t)-line(s)) * this.tangentAt(t) == 0 (2)
        // this(t) * line.dir - line(s) * line.dir == 0
        // this(t) * line.dir - line.anchor * line.dir - s line.dir * line.dir == 0
        // this(t) * line.dir - line.anchor * line.dir == s (3)
        // insert (3) in (2)
        // (this(t)-line(this(t) * line.dir - line.anchor * line.dir)) * this.tangentAt(t) == 0 (4)
        // (4) is a 5th degree polynomial, solve numerically
        tMin = isFinite(tMin) ? tMin : this.tMin;
        tMax = isFinite(tMax) ? tMax : this.tMax;
        const anchorDotDir1 = line.anchor.dot(line.dir1);
        const f = (t) => {
            const atT = this.at(t);
            return (atT.minus(line.at(atT.dot(line.dir1) - anchorDotDir1))).dot(this.tangentAt(t));
        };
        const STEPS = 32;
        const startT = ts3dutils.arrayFromFunction(STEPS, i => tMin + (tMax - tMin) * i / STEPS).withMax(t => -f(t));
        return ts3dutils.newtonIterate1d(f, startT, 8);
    }
    /**
     *
     * @param bezier
     * @param tMin
     * @param tMax
     * @param sMin
     * @param {number=} sMax
     * @returns
     */
    isInfosWithBezie3(bezier, tMin, tMax, sMin, sMax) {
        const handleStartTS = (startT, startS) => {
            if (!result.some(info => ts3dutils.eq(info.tThis, startT) && ts3dutils.eq(info.tOther, startS))) {
                const f1 = (t, s) => this.tangentAt(t).dot(this.at(t).minus(bezier.at(s)));
                const f2 = (t, s) => bezier.tangentAt(s).dot(this.at(t).minus(bezier.at(s)));
                // f = (b1, b2, t1, t2) = b1.tangentAt(t1).dot(b1.at(t1).minus(b2.at(t2)))
                const fdt1 = (b1, b2, t1, t2) => b1.ddt(t1).dot(b1.at(t1).minus(b2.at(t2))) + (b1.tangentAt(t1).squared());
                const fdt2 = (b1, b2, t1, t2) => -b1.tangentAt(t1).dot(b2.tangentAt(t2));
                const ni = ts3dutils.newtonIterate2dWithDerivatives(f1, f2, startT, startS, 16, fdt1.bind(undefined, this, bezier), fdt2.bind(undefined, this, bezier), (t, s) => -fdt2(bezier, this, s, t), (t, s) => -fdt1(bezier, this, s, t));
                result.push({ tThis: ni.x, tOther: ni.y, p: this.at(ni.x) });
            }
        };
        tMin = 'number' == typeof tMin && isFinite(tMin) ? tMin : this.tMin;
        tMax = 'number' == typeof tMax && isFinite(tMax) ? tMax : this.tMax;
        sMin = 'number' == typeof sMin && isFinite(sMin) ? sMin : bezier.tMin;
        sMax = 'number' == typeof sMax && isFinite(sMax) ? sMax : bezier.tMax;
        // stack of indices:
        const indices = [tMin, tMax, sMin, sMax];
        const tMid = (tMin + tMax) / 2;
        const sMid = (sMin + sMax) / 2;
        const aabbs = [this.getAABB(tMin, tMid), this.getAABB(tMid, tMax), bezier.getAABB(sMin, sMin), bezier.getAABB(sMid, sMax)];
        const result = [];
        while (indices.length) {
            const i = indices.length - 4;
            const tMin = indices[i], tMax = indices[i + 1], sMin = indices[i + 2], sMax = indices[i + 3];
            indices.length -= 4;
            const thisAABB = this.getAABB(tMin, tMax);
            const otherAABB = bezier.getAABB(sMin, sMax);
            // console.log(tMin, tMax, sMin, sMax, thisAABB.sce, otherAABB.sce)
            if (thisAABB && otherAABB && thisAABB.intersectsAABB2d(otherAABB)) {
                const tMid = (tMin + tMax) / 2;
                const sMid = (sMin + sMax) / 2;
                const EPS = 0.00001;
                if (tMax - tMin < EPS || sMax - sMin < EPS) {
                    console.log(tMin, tMax, sMin, sMax);
                    console.log(thisAABB.sce);
                    console.log(otherAABB.sce);
                    console.log(tMid, sMid);
                    handleStartTS(tMid, sMid);
                }
                else {
                    Array.prototype.push.call(indices, tMin, tMid, sMin, sMid, tMin, tMid, sMid, sMax, tMid, tMax, sMin, sMid, tMid, tMax, sMid, sMax);
                }
            }
        }
        return result;
    }
    isInfosWithBezier(bezier, tMin, tMax, sMin, sMax) {
        tMin = 'number' == typeof tMin && isFinite(tMin) ? tMin : this.tMin;
        tMax = 'number' == typeof tMax && isFinite(tMax) ? tMax : this.tMax;
        sMin = 'number' == typeof sMin && isFinite(sMin) ? sMin : bezier.tMin;
        sMax = 'number' == typeof sMax && isFinite(sMax) ? sMax : bezier.tMax;
        ts3dutils.assertf(() => tMin < tMax);
        ts3dutils.assertf(() => sMin < sMax);
        const result = [];
        const likeCurves = this.likeCurve(bezier), colinearCurves = this.isColinearTo(bezier);
        if (likeCurves || colinearCurves) {
            if (!likeCurves) {
                // only colinear
                // recalculate sMin and sMax so they are valid on this, from then on we can ignore bezier
                sMin = this.pointT(bezier.at(sMin));
                sMax = this.pointT(bezier.at(sMax));
            }
            tMin = Math.min(tMin, sMin);
            tMax = Math.max(tMax, sMax);
            const splits = ts3dutils.fuzzyUniques(this.roots().concatenated().filter(isFinite).concat([tMin, tMax])).sort(ts3dutils.MINUS);
            //const aabbs = arrayFromFunction(splits.length - 1, i => this.getAABB(splits[i], splits[i + 1]))
            Array.from(ts3dutils.combinations(splits.length - 1)).forEach(({ i, j }) => {
                // adjacent curves can't intersect
                if (Math.abs(i - j) > 2) {
                    // console.log(splits[i], splits[i + 1], splits[j], splits[j + 1], aabbs[i], aabbs[j])
                    //findRecursive(splits[i], splits[i + 1], splits[j], splits[j + 1], aabbs[i], aabbs[j])
                    result.push(...Curve.ispsRecursive(this, splits[i], splits[i + 1], bezier, splits[j], splits[j + 1]));
                }
            });
        }
        else {
            return Curve.ispsRecursive(this, tMin, tMax, bezier, sMin, sMax);
        }
        return result;
    }
    selfIntersectionsInfo() {
        return this.isInfosWithBezier(this);
    }
    isInfosWithCurve(curve) {
        if (curve instanceof L3$1) {
            return this.isInfosWithLine(curve.anchor, curve.dir1, curve.tMin, curve.tMax);
        }
        if (curve instanceof BezierCurve) {
            return this.isInfosWithBezier(curve);
        }
        return curve.isInfosWithCurve(this).map(({ tThis, tOther, p }) => ({ tThis: tOther, tOther: tThis, p }));
    }
    getAreaInDirSurface(dir1, surface, aT, bT) {
        ts3dutils.assertf(() => dir1.hasLength(1));
        // INT[aT; bT] at(t) * dir1 * tangentAt(t).rejectedFrom(dir1) dt
        const f = (t) => {
            const tangent = this.tangentAt(t);
            const at = this.at(t);
            const outsideVector = tangent.cross(surface.normalP(at));
            const sign = Math.sign(outsideVector.dot(dir1));
            return at.dot(dir1) * tangent.rejected1Length(dir1) * sign;
            //return this.at(t).dot(dir1) * tangent.minus(dir1.times(tangent.dot(dir1))).length()
        };
        const cx = (t) => {
            const height = this.at(t).dot(dir1);
            //console.log(t, this.at(t).minus(dir1.times(height / 2)).sce, f(t))
            return this.at(t).minus(dir1.times(height / 2));
        };
        const area = ts3dutils.gaussLegendreQuadrature24(f, aT, bT);
        const x = ts3dutils.V3.add.apply(undefined, ts3dutils.arrayFromFunction(24, i => {
            const t = aT + (ts3dutils.gaussLegendre24Xs[i] + 1) / 2 * (bT - aT);
            return cx(t).times(ts3dutils.gaussLegendre24Weights[i] * f(t));
        })).div(2 * (bT - aT) * area);
        return { area: area, centroid: x };
    }
    magic(t0 = this.tMin, t1 = this.tMax, result = []) {
        const splits = 20;
        const ts = ts3dutils.arrayFromFunction(splits, i => ts3dutils.lerp(t0, t1, i / (splits - 1)));
        const ps = ts.map(t => this.at(t));
        const ns = ts.map(t => this.normalP(t).unit());
        const f = (ns) => {
            const ls = ts.map((t, i) => new L3$1(ps[i], ns[i].unit()));
            const isInfos = ts3dutils.arrayFromFunction(splits - 1, i => {
                const j = i + 1;
                const li = ls[i], lj = ls[j];
                return li.infoClosestToLine(lj);
            });
            const a = isInfos.map(isInfo => isInfo.s - isInfo.t);
            const centers = isInfos.map(isInfo => ts3dutils.V3.lerp(isInfo.closest, isInfo.closest2, 0.5));
            const b = ts3dutils.arrayFromFunction(splits - 1, i => {
                const tMid = ts3dutils.lerp(ts[i], ts[i + 1], 0.5);
                const pMid = this.at(tMid);
                return Math.pow(pMid.distanceTo(centers[i]), 0.5);
            });
            return a.concat(b);
        };
        const startX = ts3dutils.V3.packXY(ns);
        const ff = (xs) => {
            return f(ts3dutils.V3.unpackXY(xs));
        };
        const x = new ts3dutils.Vector(new Float64Array(startX));
        for (let i = 0; i < 2; i++) {
            const Fx = new ts3dutils.Vector(new Float64Array(ff(x.v)));
            console.log(Fx.v);
            const jacobi = ts3dutils.Matrix.jacobi(ff, x.v);
            console.log('jacobi\n', jacobi.toString(x => '' + x));
            const jacobiDependentRowIndexes = jacobi.getDependentRowIndexes();
            //if (0 != jacobiDependentRowIndexes.length) {
            //	const error:any = new Error()
            //	error.jacobiDependentRowIndexes = jacobiDependentRowIndexes
            //	throw error
            //}
            const jacobiTranspose = jacobi.transposed();
            console.log((jacobi.times(jacobiTranspose)).str);
            console.log((jacobi.times(jacobiTranspose)).inversed().str);
            const matrix = jacobiTranspose.times((jacobi.times(jacobiTranspose)).inversed());
            const xDiff = matrix.timesVector(Fx);
            x = x.minus(xDiff);
        }
        const ns2 = ts3dutils.V3.unpackXY(x.v);
        const ls2 = ts3dutils.arrayFromFunction(splits, i => new L3$1(ps[i], ns2[i].unit()));
        const curves = ts3dutils.arrayFromFunction(splits - 1, i => {
            const j = i + 1;
            const li = ls2[i], lj = ls2[j];
            const isInfo = li.infoClosestToLine(lj);
            return EllipseCurve.circleForCenter2P(isInfo.closest, ps[i], ps[j], isInfo.s);
        });
        return curves;
    }
    magic2(t0 = this.tMin, t1 = this.tMax, result = []) {
        const max3d = 0.01, eps = 0.01;
        const a = this.at(t0), b = this.at(t1);
        const aN = this.normalP(t0).unit(), bN = this.normalP(t1).unit();
        const aL = new L3$1(a, aN), bL = new L3$1(b, bN);
        const isInfo = aL.infoClosestToLine(bL);
        if (isInfo.s < 0 || isInfo.t < 0
            || isInfo.distance > max3d
            || !ts3dutils.eq(isInfo.s, isInfo.t, eps)) {
        }
        else {
            const centerPoint = ts3dutils.V3.lerp(isInfo.closest, isInfo.closest2, 0.5);
            const testT1 = ts3dutils.lerp(t0, t1, 1 / 2), testP1 = this.at(testT1);
            const testT2 = ts3dutils.lerp(t0, t1, 2 / 3), testP2 = this.at(testT2);
            const radius = (isInfo.s + isInfo.t) / 2;
            if (ts3dutils.eq(centerPoint.distanceTo(testP1), radius, eps)) {
                const newCurve = EllipseCurve.circleForCenter2P(centerPoint, a, b, radius);
                result.push(newCurve);
                return result;
            }
        }
        const tMid = (t0 + t1) / 2;
        this.magic(t0, tMid, result);
        this.magic(tMid, t1, result);
        return result;
    }
}
/**
 * https://en.wikipedia.org/wiki/Cubic_function#/media/File:Graph_of_cubic_polynomial.svg
 */
BezierCurve.EX2D = BezierCurve.graphXY(2, -3, -3, 2);
BezierCurve.EX3D = new BezierCurve(ts3dutils.V3.O, ts3dutils.V(-0.1, -1, 1), ts3dutils.V(1.1, 1, 1), ts3dutils.V3.X);
BezierCurve.QUARTER_CIRCLE = BezierCurve.approximateUnitArc(PI$1 / 2);
BezierCurve.prototype.hlol = Curve.hlol++;
BezierCurve.prototype.tIncrement = 1 / 80;

const { PI: PI$2, cos: cos$1, sin: sin$1, min, max, tan, sign, ceil: ceil$2, floor: floor$2, abs: abs$3, sqrt, pow, atan2, round } = Math;
class EllipseCurve extends XiEtaCurve {
    constructor(center, f1, f2, tMin = -PI$2, tMax = PI$2) {
        super(center, f1, f2, tMin, tMax);
        ts3dutils.assert(EllipseCurve.isValidT(tMin));
        ts3dutils.assert(EllipseCurve.isValidT(tMax));
    }
    static isValidT(t) {
        return -Math.PI <= t && t <= Math.PI;
    }
    static XYLCValid(pLC) {
        return ts3dutils.eq(1, pLC.lengthXY());
    }
    /**
     * @param hint +-PI, whichever is correct
     */
    static XYLCPointT(pLC, hint) {
        const angle = pLC.angleXY();
        if (angle < -Math.PI + ts3dutils.NLA_PRECISION || angle > Math.PI - ts3dutils.NLA_PRECISION) {
            ts3dutils.assert(isFinite(hint));
            return Math.sign(hint) * Math.PI;
        }
        return angle;
    }
    static magic(a, b, c) {
        const isLC = intersectionUnitCircleLine2(a, b, c);
        return isLC.map(([xi, eta]) => Math.atan2(eta, xi));
    }
    static unitIsInfosWithLine(anchorLC, dirLC, anchorWC, dirWC) {
        // ell: x² + y² = 1 = p²
        // line(t) = anchor + t dir
        // anchor² - 1 + 2 t dir anchor + t² dir² = 0
        const pqDiv = dirLC.dot(dirLC);
        const lineTs = ts3dutils.pqFormula(2 * dirLC.dot(anchorLC) / pqDiv, (anchorLC.dot(anchorLC) - 1) / pqDiv);
        return lineTs.map(tOther => ({
            tThis: Math.atan2(anchorLC.y + tOther * dirLC.y, anchorLC.x + tOther * dirLC.x),
            tOther: tOther,
            p: L3$1.at(anchorWC, dirWC, tOther),
        }));
    }
    /**
     * Returns a new EllipseCurve representing a circle parallel to the XY-plane.`
     */
    static circle(radius, center = ts3dutils.V3.O) {
        return new EllipseCurve(center, new ts3dutils.V3(radius, 0, 0), new ts3dutils.V3(0, radius, 0));
    }
    static circleForCenter2P(center, a, b, radius) {
        const f1 = center.to(a);
        const normal = f1.cross(center.to(b));
        const f2 = normal.cross(f1).toLength(f1.length());
        const tMax = f1.angleTo(center.to(b));
        return new EllipseCurve(center, f1, f2, 0, tMax);
    }
    // TODO: there'S alsoa commented out test
    getVolZAnd(dir1, tStart, tEnd) {
        // let p = at(t)
        // integrate area [p -> plane.projectPoint(p)] to x axis...
        // INTEGRATE[tStart, tEnd] fp(this.at(t)) dt
        function fp(p) {
            const p0ToP = dir1.times(dir1.dot(p));
            const area = p0ToP.lengthXY() * (p.z - p0ToP.z / 2);
            return area;
        }
        const f = (t) => fp(this.at(t)) * this.tangentAt(t).cross(this.normal).unit().z;
        return { volume: ts3dutils.glqInSteps(f, tStart, tEnd, 4), centroid: undefined };
    }
    getAreaInDir(right, up, tStart, tEnd) {
        //assertf(() => tStart < tEnd)
        ts3dutils.assertf(() => right.isPerpendicularTo(this.normal));
        ts3dutils.assertf(() => up.isPerpendicularTo(this.normal));
        //assertf(() => EllipseCurve.isValidT(tStart), tStart)
        //assertf(() => EllipseCurve.isValidT(tEnd), tEnd)
        const upLC = this.inverseMatrix.transformVector(up);
        const rightLC = upLC.cross(ts3dutils.V3.Z);
        const normTStart = tStart - rightLC.angleXY();
        const normTEnd = tEnd - rightLC.angleXY();
        const transformedOriginY = this.inverseMatrix.getTranslation().dot(upLC.unit());
        //console.log(upLC.str, rightLC.str, normTStart, normTEnd, 'upLC.length()', upLC.length())
        //console.log('transformedOriginY', transformedOriginY)
        //assertf(() => upLC.hasLength(1), upLC.length())
        const fPi = Math.PI / 4;
        // integral of sqrt(1 - x²) from 0 to cos(t)
        // Basically, we want
        // INTEGRAL[cos(t); PI/2] sqrt(1 - x²) dx
        // INTEGRAL[PI/2: cos(t)] -sqrt(1 - x²) dx
        // = INTEGRAL[cos(0); cos(t)] -sqrt(1 - x²) dx
        // = INTEGRAL[0; t] -sqrt(1 - cos²(t)) * -sin(t) dt
        // = INTEGRAL[0; t] -sin(t) * -sin(t) dt
        // = INTEGRAL[0; t] sin²(t) dt (partial integration / wolfram alpha)
        // = (1/2 * (t - sin(t) * cos(t)))[0; t] (this form has the distinct advantage of being defined everywhere)
        function fArea(t) { return (t - Math.sin(t) * Math.cos(t)) / 2; }
        // for the centroid, we want
        // cx = 1 / area * INTEGRAL[cos(t); PI/2] x * f(x) dx
        // cx = 1 / area * INTEGRAL[cos(t); PI/2] x * sqrt(1 - x²) dx
        // cx = 1 / area * INTEGRAL[cos(0); cos(t)] x * -sqrt(1 - x²) dx
        // ...
        // cx = 1 / area * INTEGRAL[0; t] cos(t) * sin²(t) dt // WA
        // cx = 1 / area * (sin^3(t) / 3)[0; t]
        function cxTimesArea(t) { return Math.pow(Math.sin(t), 3) / 3; }
        // cy = 1 / area * INTEGRAL[cos(t); PI/2] f²(x) / 2 dx
        // cy = 1 / area * INTEGRAL[cos(0); cos(t)] -(1 - x²) / 2 dx
        // cy = 1 / area * INTEGRAL[0; t] (cos²(t) - 1) * -sin(t) / 2 dt
        // cy = 1 / area * (cos (3 * t) - 9 * cos(t)) / 24 )[0; t]
        function cyTimesArea(t) { return (Math.cos(3 * t) - 9 * Math.cos(t)) / 24; }
        const restArea = -transformedOriginY * (-Math.cos(normTEnd) + Math.cos(normTStart));
        const area = fArea(normTEnd) - fArea(normTStart) + restArea;
        const cxt = (cxTimesArea(normTEnd) - cxTimesArea(normTStart) + -transformedOriginY * (-Math.cos(normTEnd) - Math.cos(normTStart)) / 2 * restArea) / area;
        const cyt = (cyTimesArea(normTEnd) - cyTimesArea(normTStart) - -transformedOriginY / 2 * restArea) / area;
        const factor = this.matrix.xyAreaFactor(); // * upLC.length()
        //console.log('fctor', factor, 'area', area, 'resultarea', area* factor)
        ts3dutils.assert(!ts3dutils.eq0(factor));
        return {
            area: area * factor,
            centroid: this.matrix.transformPoint(ts3dutils.M4.rotateZ(rightLC.angleXY()).transformPoint(new ts3dutils.V3(cxt, cyt, 0))),
        };
    }
    at(t) {
        // = center + f1 cos t + f2 sin t
        return this.center.plus(this.f1.times(Math.cos(t))).plus(this.f2.times(Math.sin(t)));
    }
    tangentAt(t) {
        ts3dutils.assertNumbers(t);
        // f2 cos(t) - f1 sin(t)
        return this.f2.times(Math.cos(t)).minus(this.f1.times(Math.sin(t)));
    }
    ddt(t) {
        ts3dutils.assertNumbers(t);
        return this.f2.times(-Math.sin(t)).minus(this.f1.times(Math.cos(t)));
    }
    tangentAt2(xi, eta) {
        return this.f2.times(xi).minus(this.f1.times(eta));
    }
    isCircular() {
        return ts3dutils.eq(this.f1.length(), this.f2.length()) && this.f1.isPerpendicularTo(this.f2);
    }
    reversed() {
        return new this.constructor(this.center, this.f1, this.f2.negated(), -this.tMax, -this.tMin);
    }
    isColinearTo(curve) {
        if (!ts3dutils.hasConstructor(curve, EllipseCurve))
            return false;
        if (!this.center.like(curve.center)) {
            return false;
        }
        if (this == curve) {
            return true;
        }
        if (this.isCircular()) {
            return curve.isCircular() && ts3dutils.eq(this.f1.length(), curve.f1.length()) && this.normal.isParallelTo(curve.normal);
        }
        else {
            let { f1: f1, f2: f2 } = this.rightAngled(), { f1: c1, f2: c2 } = curve.rightAngled();
            if (f1.length() > f2.length()) {
                [f1, f2] = [f2, f1];
            }
            if (c1.length() > c2.length()) {
                [c1, c2] = [c2, c1];
            }
            return ts3dutils.eq(f1.squared(), Math.abs(f1.dot(c1)))
                && ts3dutils.eq(f2.squared(), Math.abs(f2.dot(c2)));
        }
    }
    eccentricity() {
        const mainAxes = this.rightAngled();
        const f1length = mainAxes.f1.length(), f2length = mainAxes.f1.length();
        const [a, b] = f1length > f2length ? [f1length, f2length] : [f2length, f1length];
        return Math.sqrt(1 - b * b / a / a);
    }
    circumference() {
        return this.arcLength(-Math.PI, Math.PI);
    }
    arcLength(startT, endT, steps) {
        ts3dutils.assert(startT < endT, 'startT < endT');
        if (this.isCircular()) {
            return this.f1.length() * (endT - startT);
        }
        return super.arcLength(startT, endT, steps);
    }
    circumferenceApproximate() {
        // approximate circumference by Ramanujan
        // https://en.wikipedia.org/wiki/Ellipse#Circumference
        const { f1, f2 } = this.rightAngled(), a = f1.length(), b = f2.length();
        const h = Math.pow((a - b), 2) / Math.pow((a + b), 2);
        return Math.PI * (a + b) * (1 + 3 * h / (10 + Math.sqrt(4 - 3 * h)));
    }
    rightAngled() {
        const f1 = this.f1, f2 = this.f2, a = f1.dot(f2), b = f2.squared() - f1.squared();
        if (ts3dutils.eq0(a)) {
            return this;
        }
        const g1 = 2 * a, g2 = b + Math.sqrt(b * b + 4 * a * a);
        const { x1: xi, y1: eta } = intersectionUnitCircleLine(g1, g2, 0);
        return new EllipseCurve(this.center, f1.times(xi).plus(f2.times(eta)), f1.times(-eta).plus(f2.times(xi)));
    }
    isInfosWithEllipse(ellipse) {
        if (this.normal.isParallelTo(ellipse.normal) && ts3dutils.eq0(this.center.minus(ellipse.center).dot(ellipse.normal))) {
            // ellipses are coplanar
            const ellipseLCRA = ellipse.transform(this.inverseMatrix).rightAngled();
            const r1 = ellipseLCRA.f1.lengthXY(), r2 = ellipseLCRA.f2.lengthXY(), centerDist = ellipseLCRA.center.lengthXY();
            const rMin = min(r1, r2), rMax = max(r1, r2);
            if (ts3dutils.lt(centerDist + rMax, 1) || // entirely inside unit circle
                ts3dutils.lt(1, centerDist - rMax) || // entirely outside unit circle
                ts3dutils.lt(1, rMin - centerDist) || // contains unit circle
                ts3dutils.eq(1, r1) && ts3dutils.eq(1, r2) && ts3dutils.eq0(centerDist) // also unit circle, return no IS
            ) {
                return [];
            }
            const f = (t) => ellipseLCRA.at(t).lengthXY() - 1;
            const df = (t) => ellipseLCRA.at(t).xy().dot(ellipseLCRA.tangentAt(t)) / ellipseLCRA.at(t).lengthXY();
            ts3dutils.checkDerivate(f, df, -PI$2, PI$2, 1);
            const ts = [];
            const tsvs = ts3dutils.arrayRange(-4 / 5 * PI$2, PI$2, PI$2 / 4).map(startT => [startT, df(startT), ts3dutils.newtonIterateSmart(f, startT, 16, df, 1e-4), f(ts3dutils.newtonIterateSmart(f, startT, 16, df, 1e-4))]);
            for (let startT = -4 / 5 * PI$2; startT < PI$2; startT += PI$2 / 4) {
                let t = ts3dutils.newtonIterateSmart(f, startT, 16, df, 1e-4);
                ts3dutils.le(t, -PI$2) && (t += ts3dutils.TAU);
                ts3dutils.assert(!isNaN(t));
                if (ellipseLCRA.isValidT(t) && ts3dutils.eq0(f(t)) && !ts.some(r => ts3dutils.eq(t, r))) {
                    ts.push(t);
                }
            }
            return ts.map(raT => {
                const p = this.matrix.transformPoint(ellipseLCRA.at(raT));
                return { tThis: this.pointT(p), tOther: ellipse.pointT(p, PI$2), p };
            });
            //const angle = ellipseLCRA.f1.angleXY()
            //const aSqr = ellipseLCRA.f1.squared(), bSqr = ellipseLCRA.f2.squared()
            //const a = Math.sqrt(aSqr), b = Math.sqrt(bSqr)
            //const {x: centerX, y: centerY} = ellipseLCRA.center
            //const rotCenterX = centerX * Math.cos(-angle) + centerY * -Math.sin(-angle)
            //const rotCenterY = centerX * Math.sin(-angle) + centerY * Math.cos(-angle)
            //const rotCenter = V(rotCenterX, rotCenterY)
            //const f = t => {
            //	const lex = Math.cos(t) - rotCenterX, ley = Math.sin(t) - rotCenterY
            //	return lex * lex / aSqr + ley * ley / bSqr - 1
            //}
            //const f2 = (x, y) => (x * x + y * y - 1)
            //const f3 = (x, y) => ((x - rotCenterX) * (x - rotCenterX) / aSqr + (y - rotCenterY) * (y - rotCenterY) /
            // bSqr - 1) const results = [] const resetMatrix = this.matrix.times(M4.rotateZ(angle)) for (let startT =
            // Math.PI / 4; startT < 2 * Math.PI; startT += Math.PI / 2) { const startP = EllipseCurve.XY.at(startT)
            // const p = newtonIterate2d(f3, f2, startP.x, startP.y, 10) if (p && !results.some(r => r.like(p))) {
            // results.push(p) } } const rotEl = new EllipseCurve(rotCenter, V(a, 0, 0), V(0, b, 0)) return
            // results.map(pLC => { const p = resetMatrix.transformPoint(pLC) return {tThis: this.pointT(p, PI),
            // tOther: ellipse.pointT(p, PI), p} })
        }
        else {
            return this.isTsWithPlane(ellipse.getPlane()).mapFilter(t => {
                const p = this.at(t);
                if (ellipse.containsPoint(p)) {
                    return { tThis: t, tOther: ellipse.pointT(p), p };
                }
            });
        }
    }
    isInfosWithCurve(curve) {
        if (curve instanceof EllipseCurve) {
            return this.isInfosWithEllipse(curve);
        }
        return super.isInfosWithCurve(curve);
    }
    roots() {
        // tangent(t) = f2 cos t - f1 sin t
        // solve for each dimension separately
        // tangent(eta, xi) = f2 eta - f1 xi
        return ts3dutils.arrayFromFunction(3, dim => {
            const a = this.f2.e(dim), b = -this.f1.e(dim);
            const { x1, y1, x2, y2 } = intersectionUnitCircleLine(a, b, 0);
            return [Math.atan2(y1, x1), Math.atan2(y2, x2)];
        });
    }
    closestTToPoint(p, tStart) {
        // (at(t) - p) * tangentAt(t) = 0
        // (xi f1 + eta f2 + q) * (xi f2 - eta f1) = 0
        // xi eta (f2^2-f1^2) + xi f2 q - eta² f1 f2 + xi² f1 f2 - eta f1 q = 0
        //  (xi² - eta²) f1 f2 + xi eta (f2^2-f1^2) + xi f2 q - eta f1 q = 0
        // atan2 of p is a good first approximation for the searched t
        const startT = this.inverseMatrix.transformPoint(p).angleXY();
        const pRelCenter = p.minus(this.center);
        const f = (t) => this.tangentAt(t).dot(this.f1.times(Math.cos(t)).plus(this.f2.times(Math.sin(t))).minus(pRelCenter));
        return ts3dutils.newtonIterate1d(f, startT);
    }
    area() {
        // see
        // https://upload.wikimedia.org/wikipedia/commons/thumb/4/4e/Cross_product_parallelogram.svg/220px-Cross_product_parallelogram.svg.png
        return Math.PI * this.f1.cross(this.f2).length();
    }
    angleToT(phi) {
        // atan2(y, x) = phi
        const phiDir = this.f1.unit().times(Math.cos(phi)).plus(this.f2.rejectedFrom(this.f1).unit().times(Math.sin(phi)));
        const dirLC = this.inverseMatrix.transformVector(phiDir);
        return dirLC.angleXY();
    }
}
EllipseCurve.XY = new EllipseCurve(ts3dutils.V3.O, ts3dutils.V3.X, ts3dutils.V3.Y);
EllipseCurve.prototype.hlol = Curve.hlol++;
EllipseCurve.prototype.tIncrement = 2 * Math.PI / (4 * 800);

const { PI: PI$3, cos: cos$2, sin: sin$2, min: min$1, max: max$1, tan: tan$1, sign: sign$1, ceil: ceil$3, floor: floor$3, abs: abs$4, sqrt: sqrt$1, pow: pow$1, atan2: atan2$1, round: round$1 } = Math;
/**
 * x² - y² = 1
 *
 */
class HyperbolaCurve extends XiEtaCurve {
    constructor(center, f1, f2, tMin = -7, tMax = 7) {
        super(center, f1, f2, tMin, tMax);
    }
    static XYLCValid(pLC) {
        return pLC.x > 0 && ts3dutils.eq(1, pLC.x * pLC.x - pLC.y * pLC.y);
    }
    static XYLCPointT(pLC) {
        return Math.asinh(pLC.y);
    }
    /**
     * http://www.wolframalpha.com/input/?i=x%C2%B2-y%C2%B2%3D1,ax%2Bby%3Dc
     * Minor empiric test shows asinh(eta) consistently gets more accurate results than atanh(eta/xi)
     */
    static magic(a, b, c) {
        if (ts3dutils.eq0(b)) {
            const sqrtVal = ts3dutils.snap0(Math.pow(c, 2) / Math.pow(a, 2) - 1);
            if (sqrtVal < 0 || c * a < 0) {
                return [];
            }
            else if (sqrtVal == 0) {
                return [0];
            }
            const eta1 = Math.sqrt(sqrtVal);
            return [-Math.asinh(eta1), Math.asinh(eta1)];
        }
        else if (ts3dutils.eq(abs$4(a), abs$4(b))) {
            if (ts3dutils.le(c * a, 0)) {
                return [];
            }
            const eta = sign$1(a * b) * (Math.pow(c, 2) - Math.pow(a, 2)) / 2 / a / c;
            return [Math.asinh(eta)];
        }
        else {
            const sqrtVal = ts3dutils.snap0(Math.pow(b, 2) * (-(Math.pow(a, 2)) + Math.pow(b, 2) + Math.pow(c, 2)));
            if (sqrtVal < 0) {
                return [];
            }
            const xi1 = (a * c - Math.sqrt(sqrtVal)) / (Math.pow(a, 2) - Math.pow(b, 2));
            const xi2 = (a * c + Math.sqrt(sqrtVal)) / (Math.pow(a, 2) - Math.pow(b, 2));
            const eta1 = (Math.pow(b, 2) * c - a * Math.sqrt(sqrtVal)) / (b * (Math.pow(b, 2) - Math.pow(a, 2)));
            const eta2 = (Math.pow(b, 2) * c + a * Math.sqrt(sqrtVal)) / (b * (Math.pow(b, 2) - Math.pow(a, 2)));
            return [xi1 > 0 && Math.asinh(eta1), xi2 > 0 && Math.asinh(eta2)].filter((x) => x !== false);
        }
    }
    at(t) {
        ts3dutils.assertNumbers(t);
        // = center + f1 cosh t + f2 sinh t
        return this.center.plus(this.f1.times(Math.cosh(t))).plus(this.f2.times(Math.sinh(t)));
    }
    tangentAt(t) {
        ts3dutils.assertNumbers(t);
        // = f1 sinh t + f2 cosh t
        return this.f1.times(Math.sinh(t)).plus(this.f2.times(Math.cosh(t)));
    }
    tangentAt2(xi, eta) {
        ts3dutils.assertNumbers(xi, eta);
        // = f1 eta + f2 xi
        return this.f1.times(eta).plus(this.f2.times(xi));
    }
    ddt(t) {
        ts3dutils.assertNumbers(t);
        return this.f1.times(Math.cosh(t)).plus(this.f2.times(Math.sinh(t)));
    }
    isColinearTo(curve) {
        if (!ts3dutils.hasConstructor(curve, HyperbolaCurve))
            return false;
        if (!curve.center || !this.center.like(curve.center)) {
            return false;
        }
        if (this === curve) {
            return true;
        }
        const { f1: f1, f2: f2 } = this.rightAngled(), { f1: c1, f2: c2 } = curve.rightAngled();
        return ts3dutils.eq(f1.squared(), Math.abs(f1.dot(c1)))
            && ts3dutils.eq(f2.squared(), Math.abs(f2.dot(c2)));
    }
    reversed() {
        return new HyperbolaCurve(this.center, this.f1, this.f2.negated(), -this.tMax, -this.tMin);
    }
    rightAngled() {
        const f1 = this.f1, f2 = this.f2, a = f1.dot(f2), b = f2.squared() + f1.squared();
        if (ts3dutils.eq0(a)) {
            return this;
        }
        const g1 = 2 * a, g2 = b + Math.sqrt(b * b - 4 * a * a);
        const { x1: xi, y1: eta } = intersectionUnitHyperbolaLine(g1, g2, 0);
        return new HyperbolaCurve(this.center, f1.times(xi).plus(f2.times(eta)), f1.times(eta).plus(f2.times(xi)));
    }
    eccentricity() {
        const mainAxes = this.rightAngled();
        const f1length = mainAxes.f1.length(), f2length = mainAxes.f1.length();
        const [a, b] = f1length > f2length ? [f1length, f2length] : [f2length, f1length];
        return Math.sqrt(1 + b * b / a / a);
    }
    roots() {
        // tangent(t) = f1 sinh t + f2 cosh t = 0
        // tangentAt2(xi, eta) = f1 eta + f2 xi = V3.O
        // xi² - eta² = 1 (by def for hyperbola)
        return ts3dutils.arrayFromFunction(3, dim => {
            const a = this.f2.e(dim), b = this.f1.e(dim);
            return HyperbolaCurve.magic(a, b, 0);
        });
    }
}
HyperbolaCurve.XY = new HyperbolaCurve(ts3dutils.V3.O, ts3dutils.V3.X, ts3dutils.V3.Y);
HyperbolaCurve.prototype.tIncrement = PI$3 / 16;

/**
 * 3-dimensional line
 */
class L3$1 extends Curve {
    constructor(anchor, // line anchor
    dir1, // normalized line dir
    tMin = -4096, tMax = 4096) {
        super(tMin, tMax);
        this.anchor = anchor;
        this.dir1 = dir1;
        ts3dutils.assertVectors(anchor, dir1);
        ts3dutils.assert(dir1.hasLength(1), 'dir must be unit' + dir1);
        ts3dutils.assertf(() => !Number.isNaN(anchor.x));
    }
    static throughPoints(anchor, b, tMin, tMax) {
        return new L3$1(anchor, b.minus(anchor).unit(), tMin, tMax);
    }
    static pointT(anchor, dir, x) {
        ts3dutils.assertVectors(anchor, dir, x);
        return x.minus(anchor).dot(dir) / dir.squared();
    }
    static at(anchor, dir, t) {
        return anchor.plus(dir.times(t));
    }
    static fromPlanes(p1, p2) {
        ts3dutils.assertInst(P3, p1, p2);
        const dir = p1.normal1.cross(p2.normal1);
        const length = dir.length();
        if (length < 1e-10) {
            throw new Error('Parallel planes');
        }
        return p1.intersectionWithPlane(p2);
    }
    static containsPoint(anchor, dir, p) {
        const closestT = L3$1.pointT(anchor, dir, p);
        const distance = L3$1.at(anchor, dir, closestT).distanceTo(p);
        return ts3dutils.eq0(distance);
    }
    addToMesh(mesh, res = 4, radius = 0, pointStep = 1, tMin = this.tMin, tMax = this.tMax) {
        const baseNormals = ts3dutils.arrayFromFunction(res, i => ts3dutils.V3.polar(1, ts3dutils.TAU * i / res));
        const baseVertices = ts3dutils.arrayFromFunction(res, i => ts3dutils.V3.polar(radius, ts3dutils.TAU * i / res));
        const inc = this.tIncrement;
        const start = Math.ceil((this.tMin + ts3dutils.NLA_PRECISION) / inc);
        const end = Math.floor((this.tMax - ts3dutils.NLA_PRECISION) / inc);
        for (let i = 0; i <= 1; i += 1) {
            const start = mesh.vertices.length;
            if (0 !== i) {
                for (let j = 0; j < res; j++) {
                    tsgl.pushQuad(mesh.TRIANGLES, true, start - res + j, start + j, start - res + (j + 1) % res, start + (j + 1) % res);
                }
            }
            const t = 0 == i ? tMin : tMax;
            const point = this.at(t), tangent = this.dir1, x = this.dir1.getPerpendicular();
            const matrix = ts3dutils.M4.forSys(x, this.dir1.cross(x), this.dir1, point);
            mesh.normals.push(...matrix.transformedVectors(baseNormals));
            mesh.vertices.push(...matrix.transformedPoints(baseVertices));
        }
    }
    roots() {
        return [[], [], []];
    }
    containsPoint(p) {
        ts3dutils.assertVectors(p);
        const dist = this.distanceToPoint(p);
        ts3dutils.assertNumbers(dist);
        return ts3dutils.eq0(dist);
    }
    likeCurve(curve) {
        return this == curve ||
            ts3dutils.hasConstructor(curve, L3$1)
                && this.anchor.like(curve.anchor)
                && this.dir1.like(curve.dir1);
    }
    equals(obj) {
        return this == obj ||
            Object.getPrototypeOf(obj) == L3$1.prototype
                && this.anchor.equals(obj.anchor)
                && this.dir1.equals(obj.dir1);
    }
    isColinearTo(obj) {
        return obj instanceof L3$1
            && this.containsPoint(obj.anchor)
            && ts3dutils.eq(1, Math.abs(this.dir1.dot(obj.dir1)));
    }
    distanceToLine(line) {
        ts3dutils.assertInst(L3$1, line);
        if (this.isParallelToLine(line)) {
            return this.distanceToPoint(line.anchor);
        }
        const dirCross1 = this.dir1.cross(line.dir1).unit();
        const anchorDiff = this.anchor.minus(line.anchor);
        return Math.abs(anchorDiff.dot(dirCross1));
    }
    distanceToPoint(x) {
        ts3dutils.assertVectors(x);
        // See http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
        const t = x.minus(this.anchor).dot(this.dir1);
        return this.at(t).distanceTo(x);
        //return x.minus(this.anchor).cross(x.minus(this.anchor.plus(this.dir1))).length()
    }
    asSegmentDistanceToPoint(x, sStart, sEnd) {
        let t = x.minus(this.anchor).dot(this.dir1);
        t = ts3dutils.clamp(t, sStart, sEnd);
        return this.at(t).minus(x).length();
    }
    asSegmentDistanceToLine(line, sStart, sEnd) {
        ts3dutils.assertInst(L3$1, line);
        const dirCross = this.dir1.cross(line.dir1);
        const div = dirCross.squared();
        if (ts3dutils.eq0(div)) {
            return undefined;
        } // lines parallel
        const anchorDiff = line.anchor.minus(this.anchor);
        // check if distance is zero (see also L3.distanceToLine)
        if (!ts3dutils.eq0(anchorDiff.dot(dirCross.unit()))) {
            return undefined;
        }
        let t = this.infoClosestToLine(line).t;
        t = ts3dutils.clamp(t, sStart, sEnd);
        return this.at(ts3dutils.clamp(t, sStart, sEnd));
    }
    at(t) {
        ts3dutils.assertNumbers(t);
        return this.anchor.plus(this.dir1.times(t));
    }
    /**
     * This function returns lambda for a given point x
     *
     * Every point x on this line is described by the equation
     *      x = this.anchor + lambda * this.dir1 | - this.anchor
     *      x - this.anchor = lambda * this.dir1 | DOT this.dir1
     *      (x - this.anchor) DOT this.dir1 = lambda (dir1² is 1 as |dir1| == 1)
     *
     *  @param x
     *  @returns
     */
    pointT(x) {
        ts3dutils.assertVectors(x);
        const t = x.minus(this.anchor).dot(this.dir1);
        return t;
    }
    /**
     * Returns true if the line is parallel (this.dir = line.dir || this.dir = -line.dir) to the argument.
     */
    isParallelToLine(line) {
        ts3dutils.assertInst(L3$1, line);
        // we know that 1 == this.dir1.length() == line.dir1.length(), we can check for parallelity simpler than
        // isParallelTo()
        return ts3dutils.eq(1, Math.abs(this.dir1.dot(line.dir1)));
    }
    angleToLine(line) {
        ts3dutils.assertInst(L3$1, line);
        return this.dir1.angleTo(line.dir1);
    }
    /**
     *
     * @param line
     * @returns {boolean} If the distance between the lines is zero
     */
    intersectsLine(line) {
        return ts3dutils.eq0(this.distanceToLine(line));
    }
    isInfosWithCurve(curve) {
        if (curve instanceof L3$1) {
            const dirCross = this.dir1.cross(curve.dir1);
            const div = dirCross.squared();
            if (ts3dutils.eq0(div)) {
                // lines are parallel
                return [];
            }
            const anchorDiff = curve.anchor.minus(this.anchor);
            if (ts3dutils.eq0(anchorDiff.dot(dirCross))) {
                const tThis = anchorDiff.cross(curve.dir1).dot(dirCross) / div;
                const tOther = anchorDiff.cross(this.dir1).dot(dirCross) / div;
                const p = this.at(tThis);
                return [{ tThis: tThis, tOther: tOther, p: p }];
            }
            return [];
        }
        throw new Error();
    }
    isInfoWithLine(line) {
        // todo infos?
        ts3dutils.assertInst(L3$1, line);
        const dirCross = this.dir1.cross(line.dir1);
        const div = dirCross.squared();
        if (ts3dutils.eq0(div)) {
            return undefined;
        } // lines parallel
        const anchorDiff = line.anchor.minus(this.anchor);
        // check if distance is zero (see also L3.distanceToLine)
        if (!ts3dutils.eq0(anchorDiff.dot(dirCross.unit()))) {
            return undefined;
        }
        const t = anchorDiff.cross(line.dir1).dot(dirCross) / div;
        return this.at(t);
    }
    /**
     * returns s and t with this.at(s) == line.at(t)
     */
    intersectionLineST(line) {
        // the two points on two lines the closest two each other are the ones whose
        // connecting
        // TODO Where does this come from?
        // TODO: return value when no IS?
        ts3dutils.assertInst(L3$1, line);
        const dirCross = this.dir1.cross(line.dir1);
        const div = dirCross.squared();
        const anchorDiff = line.anchor.minus(this.anchor);
        const s = anchorDiff.cross(this.dir1).dot(dirCross) / div;
        const t = anchorDiff.cross(line.dir1).dot(dirCross) / div;
        return { s: s, t: t };
        //console.log(segmentIntersectsRay, a, b, "ab", ab, "p", p, "dir", dir, s > 0 && t / div >= 0 && t / div <= 1,
        // "s", s, "t", t, "div", div)
    }
    ddt(t) {
        return ts3dutils.V3.O;
    }
    getConstructorParameters() {
        return [this.anchor, this.dir1];
    }
    closestTToPoint(p) {
        // similar logic as pointT; we project the vector (anchor -> p) onto dir1, then add anchor back to it
        const nearestT = p.minus(this.anchor).dot(this.dir1);
        return nearestT;
    }
    infoClosestToLine(line) {
        /*
         line = a + s*b
         this = c + t*d

         (this - line) * b = 0
         (this - line) * d = 0

         (a + s*b - c - t*d) * b = 0
         (a + s*b - c - t*d) * d = 0

         (a - c + s*b - t*d) * b = 0
         (a - c + s*b - t*d) * d = 0

         (a - c)*b + (s*b - t*d)*b = 0
         (a - c)*d + (s*b - t*d)*d = 0

         (a - c)*b + s*(b*b) - t*(d*b) = 0
         (a - c)*d + s*(b*d) - t*(d*d) = 0

         s = (t*(d*b) - (a - c)*b) / (b*b)
         =>
         (a - c)*d + (t*(d*b) - (a - c)*b) / (b*b)*(b*d) - t*(d*d) = 0 | * (b*b)
         (a - c)*d * (b*b) + (t*(d*b) - (a - c)*b)*(b*d) - t*(d*d) * (b*b) = 0
         (a - c)*d * (b*b) + t*(d*b)*(b*d) - (a - c)*b*(b*d) - t*(d*d) * (b*b) = 0
         t = ((a - c)*b*(b*d) - (a - c)*d * (b*b)) / ((d*b)*(b*d) - (d*d) * (b*b))
         */
        if (this.isParallelToLine(line)) {
            return { t: NaN, s: NaN, distance: this.distanceToLine(line) };
        }
        const a = line.anchor, b = line.dir1, c = this.anchor, d = this.dir1;
        const bd = b.dot(d), bb = b.squared(), dd = d.squared(), amc = a.minus(c), divisor = bd * bd - dd * bb;
        const t = (amc.dot(b) * bd - amc.dot(d) * bb) / divisor;
        const s = (amc.dot(b) * dd - amc.dot(d) * bd) / divisor;
        return {
            t: t,
            s: s,
            closest: this.at(t),
            closest2: line.at(s),
            distance: this.at(t).distanceTo(line.at(s)),
        };
    }
    intersectionWithPlane(plane) {
        // plane: plane.normal1 * p = plane.w
        // line: p=line.point + lambda * line.dir1
        const lambda = (plane.w - plane.normal1.dot(this.anchor)) / plane.normal1.dot(this.dir1);
        const point = this.anchor.plus(this.dir1.times(lambda));
        return point;
    }
    tangentAt(t) {
        return this.dir1;
    }
    isTWithPlane(plane) {
        // plane: plane.normal1 * p = plane.w
        // line: p=line.point + lambda * line.dir1
        const div = plane.normal1.dot(this.dir1);
        if (ts3dutils.eq0(div))
            return NaN;
        const lambda = (plane.w - plane.normal1.dot(this.anchor)) / div;
        return lambda;
    }
    reversed() {
        return new L3$1(this.anchor, this.dir1.negated(), -this.tMax, -this.tMin);
    }
    isTsWithPlane(plane) {
        return [this.isTWithPlane(plane)];
    }
    flipped() {
        return new L3$1(this.anchor, this.dir1.negated());
    }
    transform(m4) {
        const newAnchor = m4.transformPoint(this.anchor);
        const newDir = m4.transformVector(this.dir1);
        return new L3$1(newAnchor, newDir.unit(), this.tMin * newDir.length(), this.tMax * newDir.length());
    }
    hashCode() {
        return this.anchor.hashCode() * 31 + this.dir1.hashCode();
    }
}
L3$1.anchorDirection = (anchor, dir) => new L3$1(anchor, dir.unit());
L3$1.X = new L3$1(ts3dutils.V3.O, ts3dutils.V3.X);
L3$1.Y = new L3$1(ts3dutils.V3.O, ts3dutils.V3.Y);
L3$1.Z = new L3$1(ts3dutils.V3.O, ts3dutils.V3.Z);
L3$1.prototype.hlol = Curve.hlol++;

const { floor: floor$4, abs: abs$5, ceil: ceil$4, min: min$2, max: max$2 } = Math;
class PICurve$1 extends ImplicitCurve {
    constructor(points, tangents, parametricSurface, implicitSurface, pmPoints, pmTangents, stepSize, dir = 1, generator, tMin, tMax) {
        super(points, tangents, dir, generator, tMin, tMax);
        this.parametricSurface = parametricSurface;
        this.implicitSurface = implicitSurface;
        this.pmPoints = pmPoints;
        this.pmTangents = pmTangents;
        this.stepSize = stepSize;
        ts3dutils.assert(Array.isArray(pmPoints));
        ts3dutils.assert(dir == 1);
        ts3dutils.assert(stepSize <= 1);
        const pf = parametricSurface.pSTFunc();
        const dpds = parametricSurface.dpds();
        const dpdt = parametricSurface.dpdt();
        const didp = implicitSurface.didp.bind(implicitSurface);
        this.dids = (s, t) => didp(pf(s, t)).dot(dpds(s, t));
        this.didt = (s, t) => didp(pf(s, t)).dot(dpdt(s, t));
        for (let i = 0; i < points.length - 1; i++) {
            ts3dutils.assert(!points[i].equals(points[i + 1]));
            //assert(parametricSurface.pST(pmPoints[i].x, pmPoints[i].y).equals(points[i]))
        }
    }
    static forParametricStartEnd(ps, is, pmStart, pmEnd, stepSize = 0.02, startPMTangent, tMin, tMax) {
        const pFunc = ps.pSTFunc(), iFunc = is.implicitFunction();
        const dpds = ps.dpds();
        const dpdt = ps.dpdt();
        const didp = is.didp.bind(is);
        const mf = exports.MathFunctionR2R.forFFxFy((x, y) => iFunc(pFunc(x, y)), (s, t) => didp(pFunc(s, t)).dot(dpds(s, t)), (s, t) => didp(pFunc(s, t)).dot(dpdt(s, t)));
        const { points, tangents } = followAlgorithm2d(mf, pmStart, stepSize, ps.bounds.bind(ps), pmEnd, startPMTangent);
        return PICurve$1.forParametricPointsTangents(ps, is, points, tangents, stepSize, 1, tMin, tMax);
    }
    //	assert(!startPoint.like(endPoint))
    //	assert(ParametricSurface.is(parametricSurface))
    //	assert(ImplicitSurface.is(implicitSurface))
    //	this.parametricSurface = parametricSurface
    //	this.implicitSurface = implicitSurface
    //	if (!startPoint) {
    //	const pmPoint = curvePoint(this.implicitCurve(), V(1, 1, 0))
    //	this.startPoint = this.parametricSurface.pSTFunc()(pmPoint.x, pmPoint.y)
    //} else {
    //	this.startPoint = startPoint
    //}
    //this.endPoint = endPoint
    //this.dir = dir
    //this.isLoop = false
    //try {
    //	this.calcPoints(startPoint, endPoint)
    //	this.startPoint = startPoint
    //	this.endPoint = endPoint
    //} catch (e) {
    //	this.calcPoints(this.endPoint, this.startPoint)
    //	this.startPoint = endPoint
    //	this.endPoint = startPoint
    //}
    //this.tMin = 0
    //this.tMax = this.points.length - 1
    static forStartEnd(ps, is, start, end, stepSize = 0.02, startTangent, min, max) {
        const startPM = ps.stP(start);
        const dpds = ps.dpds()(startPM.x, startPM.y), dpdt = ps.dpdt()(startPM.x, startPM.y);
        const startPMTangent = startTangent && ts3dutils.M4.forSys(dpds, dpdt).inversed().transformVector(startTangent);
        // assert(dpds.times(startPMTangent.x).plus(dpdt.times(startPMTangent.y)).like(startTangent))
        const curve = PICurve$1.forParametricStartEnd(ps, is, startPM, ps.stP(end), stepSize, startPMTangent);
        return curve.withBounds(min && curve.pointT(min), max && curve.pointT(max));
    }
    static forParametricPointsTangents(ps, is, pmPoints, pmTangents, stepSize, dir = 1, tMin, tMax) {
        const pFunc = ps.pSTFunc(), iFunc = is.implicitFunction();
        const dpds = ps.dpds();
        const dpdt = ps.dpdt();
        const points = pmPoints.map(({ x, y }) => pFunc(x, y));
        const tangents = pmPoints.map(({ x: s, y: t }, i) => {
            const ds = dpds(s, t);
            const dt = dpdt(s, t);
            return ds.times(pmTangents[i].x).plus(dt.times(pmTangents[i].y));
            //const p = points[i]
            //return cs.normalP(p).cross(ses.normalP(p))
            //	.toLength(ds.times(pmTangents[i].x).plus(dt.times(pmTangents[i].y)).length())
        });
        return new PICurve$1(points, tangents, ps, is, pmPoints, pmTangents, stepSize, dir, undefined, tMin, tMax);
    }
    getConstructorParameters() {
        return [this.points, this.tangents,
            this.parametricSurface, this.implicitSurface,
            this.pmPoints, this.pmTangents,
            this.stepSize, this.dir,
            this.generator, this.tMin, this.tMax];
    }
    reversed() {
        ts3dutils.assertNever();
        return new PICurve$1(this.parametricSurface, this.implicitSurface, this.endPoint, this.startPoint, -this.dir);
    }
    implicitCurve() {
        const pF = this.parametricSurface.pSTFunc();
        const iF = this.implicitSurface.implicitFunction();
        return function (s, t) {
            return iF(pF(s, t));
        };
    }
    isColinearTo(curve) {
        if (curve instanceof PICurve$1) {
            if (this.equals(curve)) {
                return true;
            }
            if (this.parametricSurface.isCoplanarTo(curve.parametricSurface) && this.implicitSurface.isCoplanarTo(curve.implicitSurface)) {
            }
            return false;
            ts3dutils.assertNever();
        }
        else {
            return false;
        }
    }
    //getVerticesNo0() {
    //
    //	// TODO
    //	let start, end, arr
    //	if (!this.canon) {
    //		start = Math.floor(this.aT + 1)
    //		end = ceil(this.bT)
    //		arr = sliceCyclic(this.curve.points, start, end)
    //	} else {
    //		start = Math.floor(this.bT + 1)
    //		end = ceil(this.aT)
    //		arr = sliceCyclic(this.curve.points, start, end)
    //		console.log("this.canon", !!this.canon, arr.length, start, end, this.aT)
    //		arr.reverse()
    //	}
    //	arr.push(this.b)
    //	return arr
    //}
    containsPoint(p) {
        ts3dutils.assertVectors(p);
        const t = this.pointT(p);
        return !isNaN(t) && this.isValidT(t);
    }
    equals(obj) {
        return Object.getPrototypeOf(obj) == PICurve$1.prototype
            && this.parametricSurface.equals(obj.parametricSurface)
            && this.implicitSurface.equals(obj.implicitSurface)
            && this.points[0].equals(obj.points[0])
            && this.tangents[0].equals(obj.tangents[0])
            && this.dir === obj.dir;
    }
    hashCode() {
        let hashCode = 0;
        hashCode = hashCode * 31 + this.parametricSurface.hashCode();
        hashCode = hashCode * 31 + this.implicitSurface.hashCode();
        hashCode = hashCode * 31 + this.points[0].hashCode();
        hashCode = hashCode * 31 + this.tangents[0].hashCode();
        return hashCode | 0;
    }
    tangentP(point) {
        ts3dutils.assertVectors(point);
        ts3dutils.assert(this.containsPoint(point), 'this.containsPoint(point)' + this.containsPoint(point));
        const t = this.pointT(point);
        return this.tangentAt(t);
    }
    tangentAt(t) {
        return ts3dutils.V3.lerp(this.tangents[floor$4(t)], this.tangents[ceil$4(t)], t % 1);
    }
    at(t) {
        // assert(!isNaN(t))
        // const pointParams = this.stT(t)
        // const result = this.parametricSurface.pSTFunc()(pointParams.x, pointParams.y)
        // // assert(eq(t, this.pointT(result)))
        // return result
        ts3dutils.assert(!isNaN(t));
        if (t % 1 == 0)
            return this.points[t];
        const startParams = ts3dutils.V3.lerp(this.pmPoints[floor$4(t)], this.pmPoints[ceil$4(t)], t % 1);
        return this.closestPointToParams(startParams);
    }
    stT(t) {
        ts3dutils.assert(!isNaN(t));
        if (t % 1 == 0)
            return this.points[t];
        const startParams = ts3dutils.V3.lerp(this.pmPoints[floor$4(t)], this.pmPoints[ceil$4(t)], t % 1);
        return curvePoint(this.implicitCurve(), startParams, this.dids, this.didt);
    }
    closestTToPoint(p, tStart) {
        return 0;
    }
    closestPointToParams(startParams) {
        const pointParams = curvePoint(this.implicitCurve(), startParams, this.dids, this.didt);
        return this.parametricSurface.pSTFunc()(pointParams.x, pointParams.y);
    }
    isTsWithSurface(surface) {
        if (surface instanceof PlaneSurface$1) {
            return this.isTsWithPlane(surface.plane);
        }
        else if (surface instanceof EllipsoidSurface || surface instanceof SemiEllipsoidSurface) {
            const ps = this.parametricSurface, is = this.implicitSurface;
            if (ps instanceof ProjectedCurveSurface && is instanceof SemiEllipsoidSurface) {
                const iscs = is.isCurvesWithSurface(surface);
                const points = iscs.flatMap(isc => isc.isTsWithSurface(ps).map(t => isc.at(t)));
                const ts = ts3dutils.fuzzyUniques(points.map(p => this.pointT(p)));
                return ts.filter(t => !isNaN(t) && this.isValidT(t));
            }
        }
        throw new Error();
    }
    isTsWithPlane(plane) {
        ts3dutils.assertInst(P3, plane);
        const ps = this.parametricSurface, is = this.implicitSurface;
        const pscs = ps.isCurvesWithPlane(plane);
        const iscs = is.isCurvesWithPlane(plane);
        const infos = iscs.flatMap(isc => pscs.flatMap(psc => isc.isInfosWithCurve(psc)));
        const ts = ts3dutils.fuzzyUniques(infos.map(info => this.pointT(info.p)));
        return ts.filter(t => !isNaN(t) && this.isValidT(t));
    }
    pointT(p) {
        ts3dutils.assertVectors(p);
        if (!this.parametricSurface.containsPoint(p) || !this.implicitSurface.containsPoint(p)) {
            return NaN;
        }
        const pmPoint = this.parametricSurface.stPFunc()(p);
        const ps = this.points, pmps = this.pmPoints;
        let t = 0, prevDistance, pmDistance = pmPoint.distanceTo(pmps[0]);
        while (pmDistance > abs$5(this.stepSize) && t < ps.length - 1) {
            //console.log(t, pmps[t].$, pmDistance)
            t = min$2(pmps.length - 1, t + max$2(1, Math.round(pmDistance / abs$5(this.stepSize) / 2 / 2)));
            pmDistance = pmPoint.distanceTo(pmps[t]);
        }
        // if (t < this.pmPoints.length - 1 && pmDistance > pmPoint.distanceTo(pmps[t + 1])) {
        //     t++
        // }
        if (pmDistance > abs$5(this.stepSize) * 1.1) {
            // p is not on this curve
            return NaN;
        }
        if (t == ps.length - 1) {
            t--;
        }
        if (ps[t].like(p))
            return t;
        if (ps[t + 1].like(p))
            return t + 1;
        const startT = t + ts3dutils.V3.inverseLerp(ps[t], ps[t + 1], p);
        if (startT)
            return ts3dutils.newtonIterate1d(t => this.at(t).distanceTo(p), startT, 2);
    }
    transform(m4) {
        const dirFactor = m4.isMirroring() ? -1 : 1;
        return PICurve$1.forStartEnd(this.parametricSurface.transform(m4), this.implicitSurface.transform(m4), m4.transformPoint(this.points[0]), m4.transformPoint(this.points.last), this.stepSize * dirFactor, m4.transformVector(this.tangents[0]), m4.transformPoint(this.at(this.tMin)), m4.transformPoint(this.at(this.tMax)));
        //return PICurve.forParametricStartEnd(
        //	this.parametricSurface.transform(m4),
        //	this.implicitSurface.transform(m4),
        //	this.pmPoints[0],
        //	this.pmPoints.last,
        //	this.stepSize,
        //	this.dir,
        //	this.tMin,
        //	this.tMax)
        // TODO: pass transformed points?
        //return new PICurve(
        //	m4.transformedPoints(this.points),
        //	m4.transformedVectors(this.tangents),
        //    this.parametricSurface.transform(m4),
        //   this.implicitSurface.transform(m4),
        //   this.pmPoints,
        //   this.pmTangents,
        //this.stepSize,
        //   this.dir,
        //this.generator,
        //this.tMin, this.tMax)
    }
    roots() {
        const allTs = ts3dutils.arrayRange(0, this.points.length);
        return [allTs, allTs, allTs];
    }
    toSource(rounder = x => x) {
        const result = ts3dutils.callsce('PICurve.forParametricStartEnd', this.parametricSurface, this.implicitSurface, this.pmPoints[0], this.pmPoints.last, this.stepSize, this.pmTangents[0], this.tMin, this.tMax);
        return result;
    }
}
PICurve$1.prototype.tIncrement = 1;

/**
 * eta = xi²
 */
class ParabolaCurve extends XiEtaCurve {
    constructor(center, f1, f2, tMin = -10, tMax = 10) {
        super(center, f1, f2, tMin, tMax);
    }
    static eccentricity() {
        return 1;
    }
    static unitIsInfosWithLine(anchorLC, dirLC, anchorWC, dirWC) {
        // para: x² = y
        // line(t) = anchor + t dir
        // (ax + t dx)² = ay + t dy
        // ax² + t ax dx + t² dx² = ay + t dy
        // t² dx² + t (ax dx + dy) + ay² + ay = 0
        const pqDiv = Math.pow(dirLC.x, 2);
        const lineTs = ts3dutils.pqFormula((anchorLC.x * dirLC.x + dirLC.y) / pqDiv, (Math.pow(anchorLC.x, 2) + anchorLC.y) / pqDiv);
        return lineTs.filter(tOther => ts3dutils.le(0, anchorLC.y + tOther * dirLC.y))
            .map(tOther => ({
            tThis: dirLC.x * tOther + anchorLC.x,
            tOther: tOther,
            p: L3$1.at(anchorWC, dirWC, tOther),
        }));
    }
    static magic(a, b, c) {
        /*
         solve system (5)/(6)
         g1 * xi + g2 * eta = g3 (6)
         g1 * xi + g2 * xi * xi = g3
         xi² + xi * g1/g2 - g3/g2 = 0
         */
        return ts3dutils.pqFormula(a / b, -c / b);
    }
    static XYLCValid(pLC) {
        return ts3dutils.eq(Math.pow(pLC.x, 2), pLC.y);
    }
    static XYLCPointT(pLC) {
        return pLC.x;
    }
    static quadratic(a, b, c) {
        // (1 - t)² a + 2 * t * (1 - t) b + t² c
        // (1 -2t +t²)a + (2t -2t²) b + t² c
        // = t²(a - 2b + c) + t (-2a + 2b) + a
        // (2t - 2) a + (1 - 2t) b + 2t c = t(2a + 2b - 2c) - 2a + b
        // 2 a + -2 b + 2 c
        const f2 = a.plus(c).minus(b.times(2));
        const f1 = b.minus(a).times(2);
        const center = a;
        return new ParabolaCurve(center, f1, f2, 0, 1);
    }
    at(t) {
        // center + f1 t + f2 t²
        return this.center.plus(this.f1.times(t)).plus(this.f2.times(t * t));
    }
    tangentAt(t) {
        ts3dutils.assertNumbers(t);
        // f1 + f2 2 t
        return this.f1.plus(this.f2.times(2 * t));
    }
    ddt(t) {
        ts3dutils.assertNumbers(t);
        return this.f2.times(2);
    }
    tangentAt2(xi, eta) {
        ts3dutils.assertNumbers(xi, eta);
        return this.f1.plus(this.f2.times(2 * eta));
    }
    reversed() {
        return new this.constructor(this.center, this.f1.negated(), this.f2, -this.tMax, -this.tMin);
    }
    /**
     * tangent: f1 + 2 * t * f2 = 0
     * t = -f1 / 2 / f2 (for individual dimensions)
     */
    roots() {
        const dimRoots = (dim) => ts3dutils.eq0(this.f2.e(dim)) ? [] : [-this.f1.e(dim) / 2 / this.f2.e(dim)];
        return ts3dutils.arrayFromFunction(3, dimRoots);
    }
    isColinearTo(curve) {
        if (!ts3dutils.hasConstructor(curve, ParabolaCurve))
            return false;
        const thisRA = this.rightAngled(), curveRA = curve.rightAngled();
        return thisRA.center.like(curveRA.center)
            && thisRA.f2.like(curveRA.f2)
            && thisRA.f1.likeOrReversed(curveRA.f1);
    }
    rightAngled() {
        // looking for vertex of parabola
        // this is the point where the tangent is perpendicular to the main axis (f2)
        // tangent = f1 + f2 * 2 * t0
        // f2 DOT (f1 + f2 * 2 * t0) == 0
        // f1 DOT f2 + f2 DOT f2 * 2 * t0 == 0
        // t0 == -(f1 DOT f2) / (f2 DOT f2 * 2)
        const f1 = this.f1, f2 = this.f2;
        const f1DOTf2 = f1.dot(f2);
        if (ts3dutils.eq0(f1DOTf2) && f1.hasLength(1)) {
            return this;
        }
        const t0 = -f1DOTf2 / f2.squared() / 2;
        // we need to rearange tMin/tMax
        // tMin' = pointT(at(tMin)) =
        const raCenter = this.at(t0);
        const raF1 = this.tangentAt(t0), raF1Length = raF1.length(), raF11 = raF1.unit();
        const repos = (t) => this.at(t).minus(raCenter).dot(raF11);
        return new ParabolaCurve(raCenter, raF11, f2.div(Math.pow(raF1Length, 2)), repos(this.tMin), repos(this.tMax));
    }
    arcLength(startT, endT) {
        let f1 = this.f1;
        const f2 = this.f2;
        const f1DOTf2 = f1.dot(f2);
        let t0 = 0;
        if (!ts3dutils.eq0(f1DOTf2)) {
            t0 = -f1DOTf2 / f2.squared() / 2;
            f1 = f1.plus(f2.times(2 * t0));
        }
        const f1Length = f1.length();
        const a = f2.length() / f1Length;
        function F(x) {
            return Math.asinh(a * 2 * x) / 4 / a + x * Math.sqrt(1 + a * a * 4 * x * x) / 2;
        }
        return f1Length * (F(endT - t0) - F(startT - t0));
    }
    asBezier() {
        return BezierCurve.quadratic(this.at(-1), new L3$1(this.at(-1), this.tangentAt(-1).unit()).isInfoWithLine(new L3$1(this.at(1), this.tangentAt(1).unit())), this.at(1));
    }
}
ParabolaCurve.XY = new ParabolaCurve(ts3dutils.V3.O, ts3dutils.V3.X, ts3dutils.V3.Y);
ParabolaCurve.YZ = new ParabolaCurve(ts3dutils.V3.O, ts3dutils.V3.Y, ts3dutils.V3.Z);
ParabolaCurve.ZX = new ParabolaCurve(ts3dutils.V3.O, ts3dutils.V3.Z, ts3dutils.V3.X);
ParabolaCurve.prototype.tIncrement = 1 / 32;

const { PI: PI$4, min: min$3, max: max$3 } = Math;
class SemiEllipseCurve extends XiEtaCurve {
    constructor(center, f1, f2, tMin = 0, tMax = PI$4) {
        super(center, f1, f2, tMin, tMax);
        ts3dutils.assert(0 <= this.tMin && this.tMin < PI$4);
        ts3dutils.assert(0 < this.tMax && this.tMax <= PI$4);
    }
    static XYLCValid(pLC) {
        const { x, y } = pLC;
        return ts3dutils.le(0, y) && ts3dutils.eq0(Math.pow(x, 2) + Math.pow(y, 2) - 1);
    }
    static XYLCPointT(pLC) {
        // assert(le(0, pLC.y))
        const angle = Math.atan2(pLC.y, pLC.x);
        return angle < -PI$4 / 2 ? angle + ts3dutils.TAU : angle; // 0 ? (assert(eq0(angle) || eq(PI, abs(angle))), abs(angle)) :
        // angle
    }
    static magic(a, b, c) {
        const isLC = intersectionUnitCircleLine2(a, b, c);
        const result = [];
        for (const [xi, eta] of isLC) {
            ts3dutils.le(0, eta) && result.push(SemiEllipseCurve.XYLCPointT(new ts3dutils.V3(xi, eta, 0)));
        }
        return result;
    }
    static unitIsInfosWithLine(anchorLC, dirLC, anchorWC, dirWC) {
        // ell: x² + y² = 1 = p²
        // line(t) = anchor + t dir
        // anchor² - 1 + 2 t dir anchor + t² dir² = 0
        const pqDiv = dirLC.squared();
        const lineTs = ts3dutils.pqFormula(2 * dirLC.dot(anchorLC) / pqDiv, (anchorLC.squared() - 1) / pqDiv);
        return lineTs.filter(tOther => ts3dutils.le(0, anchorLC.y + tOther * dirLC.y))
            .map(tOther => ({
            tThis: SemiEllipseCurve.XYLCPointT(dirLC.times(tOther).plus(anchorLC)),
            tOther: tOther,
            p: L3$1.at(anchorWC, dirWC, tOther),
        }));
    }
    /**
     * Returns a new SemiEllipseCurve representing a circle parallel to the XY-plane.`
     */
    static semicircle(radius, center = ts3dutils.V3.O) {
        return new SemiEllipseCurve(center, new ts3dutils.V3(radius, 0, 0), new ts3dutils.V3(0, radius, 0));
    }
    static fromEllipse(curve, tMin, tMax) {
        return [
            tMin < 0 && new SemiEllipseCurve(curve.center, curve.f1.negated(), curve.f2.negated(), tMin + PI$4, min$3(0, tMax) + PI$4),
            tMax > 0 && new SemiEllipseCurve(curve.center, curve.f1, curve.f2, max$3(0, tMin), tMax),
        ].filter(x => x);
    }
    getAreaInDir(right, up, tStart, tEnd) {
        return EllipseCurve.prototype.getAreaInDir.call(this, right, up, tStart, tEnd);
    }
    at(t) {
        ts3dutils.assertNumbers(t);
        //assert(this.isValidT(t))
        // center + f1 cos t + f2 sin t
        return this.center.plus(this.f1.times(Math.cos(t))).plus(this.f2.times(Math.sin(t)));
    }
    tangentAt(t) {
        ts3dutils.assertNumbers(t);
        //assert(this.isValidT(t))
        // f2 cos(t) - f1 sin(t)
        return this.f2.times(Math.cos(t)).minus(this.f1.times(Math.sin(t)));
    }
    ddt(t) {
        ts3dutils.assertNumbers(t);
        ts3dutils.assert(this.isValidT(t));
        return this.f2.times(-Math.sin(t)).minus(this.f1.times(Math.cos(t)));
    }
    isCircular() {
        return ts3dutils.eq(this.f1.length(), this.f2.length()) && this.f1.isPerpendicularTo(this.f2);
    }
    isColinearTo(curve) {
        if (!((x) => x.constructor == this.constructor)(curve)) {
            return false;
        }
        if (!ts3dutils.hasConstructor(curve, SemiEllipseCurve))
            return false;
        if (!this.center.like(curve.center)) {
            return false;
        }
        if (this == curve) {
            return true;
        }
        if (this.isCircular()) {
            return curve.isCircular() && ts3dutils.eq(this.f1.length(), curve.f1.length()) && this.normal.isParallelTo(curve.normal);
        }
        else {
            let { f1: f1, f2: f2 } = this.rightAngled(), { f1: c1, f2: c2 } = curve.rightAngled();
            if (f1.length() > f2.length()) {
                [f1, f2] = [f2, f1];
            }
            if (c1.length() > c2.length()) {
                [c1, c2] = [c2, c1];
            }
            return ts3dutils.eq(f1.squared(), Math.abs(f1.dot(c1)))
                && ts3dutils.eq(f2.squared(), Math.abs(f2.dot(c2)));
        }
    }
    isValidT(t) {
        return ts3dutils.le(0, t) && ts3dutils.le(t, PI$4);
    }
    pointT(p) {
        ts3dutils.assertVectors(p);
        ts3dutils.assert(this.containsPoint(p));
        const pLC = this.inverseMatrix.transformPoint(p);
        const t = SemiEllipseCurve.XYLCPointT(pLC);
        ts3dutils.assert(this.isValidT(t));
        return t;
    }
    reversed() {
        return new SemiEllipseCurve(this.center, this.f1.negated(), this.f2, PI$4 - this.tMax, PI$4 - this.tMin);
    }
    eccentricity() {
        const mainAxes = this.rightAngled();
        const f1length = mainAxes.f1.length(), f2length = mainAxes.f1.length();
        const [a, b] = f1length > f2length ? [f1length, f2length] : [f2length, f1length];
        return Math.sqrt(1 - b * b / a / a);
    }
    circumference() {
        return this.arcLength(-Math.PI, Math.PI);
    }
    arcLength(startT, endT, steps = 2) {
        ts3dutils.assert(startT < endT, 'startT < endT');
        const f1Length = this.f1.length();
        if (ts3dutils.eq(f1Length, this.f2.length())) {
            return f1Length * (endT - startT);
        }
        return super.arcLength(startT, endT, steps);
    }
    circumferenceApproximate() {
        // approximate circumference by Ramanujan
        // https://en.wikipedia.org/wiki/Ellipse#Circumference
        const { f1, f2 } = this.rightAngled(), a = f1.length(), b = f2.length();
        const h = (a - b) * (a - b) / (a + b) / (a + b); // (a - b)² / (a + b)²
        return Math.PI * (a + b) * (1 + 3 * h / (10 + Math.sqrt(4 - 3 * h)));
    }
    /**
     * Radii of the ellipse are described by
     * q(phi) = f1 * cos(phi) + f2 * sin(phi)
     * or q(xi, eta) = f1 * xi + f2 * eta (1) with the added condition
     * xi² + eta² = 1 (2)
     * we want to find the radius where the corresponding tangent is perpendicular
     * tangent: q'(phi) = f1 * -sin(phi) + f2 * cos(phi)
     * tangent: q'(xi, eta) = f1 * -eta + f2 * xi
     * perpendicular when: q'(xi, eta) DOT q(xi, eta) = 0
     * (f1 * -eta + f2 * xi) DOT (f1 * xi + f2 * eta) = 0
     * DOT is distributive:
     * f1² * (-eta * xi) + f1 * f2 * (-eta² + xi²) + f2² * (xi * eta) = 0
     * (f2² - f1²) * (eta * xi) + f1 * f2 * (-eta² + xi²) = 0
     * a * (xi² - eta²) + b * xi * eta = 0 (2)
     * with a = f1 * f2, b = f2² - f1²
     * => (xi/eta)² + xi/eta * b/a + 1 = 0 (divide by a * eta²)
     * xi/eta = b/a/2 +- sqrt(b²/a²/4 - 1) | * 2*a*eta
     * 2 * a * xi = eta * (b +- sqrt(b² - 4 * a²))
     * g1 * xi - g2 * eta = 0 (3)
     * with g1 = 2 * a, g2 = b +- sqrt(b² - 4 * a²)
     * Solve (3), (2) with intersectionUnitCircleLine
     */
    rightAngled() {
        const f1 = this.f1, f2 = this.f2, a = f1.dot(f2), b = f2.squared() - f1.squared();
        if (ts3dutils.eq0(a)) {
            return this;
        }
        const g1 = 2 * a, g2 = b + Math.sqrt(b * b + 4 * a * a);
        const { x1: xi, y1: eta } = intersectionUnitCircleLine(g1, g2, 0);
        const f1RA = f1.times(xi).plus(f2.times(eta));
        const f2RA = f1.times(-eta).plus(f2.times(xi));
        return new SemiEllipseCurve(this.center, f1RA, f2RA);
    }
    asEllipse() {
        return new EllipseCurve(this.center, this.f1, this.f2, this.tMin, this.tMax);
    }
    isInfosWithEllipse(ellipse) {
        if (this.normal.isParallelTo(ellipse.normal) && ts3dutils.eq0(this.center.minus(ellipse.center).dot(ellipse.normal))) {
            ellipse instanceof SemiEllipseCurve && (ellipse = ellipse.asEllipse());
            return this.asEllipse().isInfosWithCurve(ellipse).filter(info => this.isValidT(info.tThis) && ellipse.isValidT(info.tOther));
        }
        else {
            return this.isTsWithPlane(P3.normalOnAnchor(ellipse.normal.unit(), ellipse.center)).mapFilter(t => {
                const p = this.at(t);
                if (ellipse.containsPoint(p)) {
                    return { tThis: t, tOther: ellipse.pointT(p), p };
                }
            });
        }
    }
    isInfosWithCurve(curve) {
        if (curve instanceof SemiEllipseCurve || curve instanceof EllipseCurve) {
            return this.isInfosWithEllipse(curve);
        }
        return super.isInfosWithCurve(curve);
    }
    roots() {
        // tangent(t) = f2 cos t - f1 sin t
        // solve for each dimension separately
        // tangent(eta, xi) = f2 eta - f1 xi
        return ts3dutils.arrayFromFunction(3, dim => {
            const a = this.f2.e(dim), b = -this.f1.e(dim);
            const { x1, y1, x2, y2 } = intersectionUnitCircleLine(a, b, 0);
            return [Math.atan2(y1, x1), Math.atan2(y2, x2)];
        });
    }
    closestTToPoint(p) {
        // (at(t) - p) * tangentAt(t) = 0
        // (xi f1 + eta f2 + q) * (xi f2 - eta f1) = 0
        // xi eta (f2^2-f1^2) + xi f2 q - eta² f1 f2 + xi² f1 f2 - eta f1 q = 0
        //  (xi² - eta²) f1 f2 + xi eta (f2^2-f1^2) + xi f2 q - eta f1 q = 0
        // atan2 of p is a good first approximation for the searched t
        const startT = this.inverseMatrix.transformPoint(p).angleXY();
        const pRelCenter = p.minus(this.center);
        const f = (t) => this.tangentAt(t).dot(this.f1.times(Math.cos(t)).plus(this.f2.times(Math.sin(t))).minus(pRelCenter));
        return ts3dutils.newtonIterate1d(f, startT);
    }
    area() {
        // see
        // https://upload.wikimedia.org/wikipedia/commons/thumb/4/4e/Cross_product_parallelogram.svg/220px-Cross_product_parallelogram.svg.png
        return Math.PI * this.f1.cross(this.f2).length();
    }
    angleToT(phi) {
        // atan2(y, x) = phi
        const phiDir = this.f1.unit().times(Math.cos(phi)).plus(this.f2.rejectedFrom(this.f1).unit().times(Math.sin(phi)));
        const localDir = this.inverseMatrix.transformVector(phiDir);
        return localDir.angleXY();
    }
}
SemiEllipseCurve.UNIT = new SemiEllipseCurve(ts3dutils.V3.O, ts3dutils.V3.X, ts3dutils.V3.Y);
SemiEllipseCurve.prototype.hlol = Curve.hlol++;
SemiEllipseCurve.prototype.tIncrement = 2 * Math.PI / (4 * 32);

class P3 extends ts3dutils.Transformable {
    /**
     * Oriented plane, i.e. splits R^3 in half, with one half being "in front" of the plane.
     * Leads to multiple comparisons: isCoplanarToPlane returns if the plane occupies the same space,
     * like returns if the plane occupies the same space and has the same orientation
     *
     * Points x on the plane fulfill the equation: normal1 DOT x = w
     *
     * @param normal1 unit plane normal1
     * @param w signed (rel to normal1) distance from the origin
     */
    constructor(normal1, w = 0) {
        super();
        this.normal1 = normal1;
        this.w = w;
        ts3dutils.assertVectors(normal1);
        ts3dutils.assertNumbers(w);
        ts3dutils.assert(normal1.hasLength(1), 'normal1.hasLength(1)' + normal1);
    }
    get anchor() {
        return this.normal1.times(this.w);
    }
    static throughPoints(a, b, c) {
        ts3dutils.assertVectors(a, b, c);
        const n1 = b.minus(a).cross(c.minus(a)).unit();
        return new P3(n1, n1.dot(a));
    }
    static normalOnAnchor(normal, anchor) {
        ts3dutils.assertVectors(normal, anchor);
        const n1 = normal.unit();
        return new P3(n1, n1.dot(anchor));
    }
    /**
     * x/x0 + y/y0 + y/y0 = 1
     *
     */
    static forAxisIntercepts(x0, y0, z0) {
        ts3dutils.assertNumbers(x0, y0, z0);
        const normal = new ts3dutils.V3(1 / x0, 1 / y0, 1 / z0);
        return new P3(normal.unit(), normal.length());
    }
    static forAnchorAndPlaneVectors(anchor, v0, v1) {
        ts3dutils.assertVectors(anchor, v0, v1);
        return P3.normalOnAnchor(v0.cross(v1), anchor);
    }
    axisIntercepts() {
        const w = this.w, n = this.normal1;
        return new ts3dutils.V3(w / n.x, w / n.y, w / n.z);
    }
    isCoplanarToPlane(plane) {
        ts3dutils.assertInst(P3, plane);
        return this.like(plane) || this.likeFlipped(plane);
    }
    like(plane) {
        ts3dutils.assertInst(P3, plane);
        return ts3dutils.eq(this.w, plane.w) && this.normal1.like(plane.normal1);
    }
    likeFlipped(plane) {
        ts3dutils.assertInst(P3, plane);
        return ts3dutils.eq(this.w, -plane.w) && this.normal1.like(plane.normal1.negated());
    }
    /**
     * True iff plane.normal1 is equal to this.normal1 or it's negation.
     *
     */
    isParallelToPlane(plane) {
        ts3dutils.assertInst(P3, plane);
        return ts3dutils.eq(1, Math.abs(this.normal1.dot(plane.normal1)));
    }
    isParallelToLine(line) {
        ts3dutils.assertInst(L3$1, line);
        return ts3dutils.eq0(this.normal1.dot(line.dir1));
    }
    isPerpendicularToLine(line) {
        ts3dutils.assertInst(L3$1, line);
        // this.normal1 || line.dir1
        return ts3dutils.eq(1, Math.abs(this.normal1.dot(line.dir1)));
    }
    isPerpendicularToPlane(plane) {
        ts3dutils.assertInst(P3, plane);
        return ts3dutils.eq0(this.normal1.dot(plane.normal1));
    }
    toSource(rounder) {
        return ts3dutils.callsce('new P3', this.normal1, this.w);
    }
    translated(offset) {
        return new P3(this.normal1, this.w + offset.dot(this.normal1));
    }
    transform(m4) {
        const mirror = m4.isMirroring();
        // get two vectors in the plane:
        const u = this.normal1.getPerpendicular();
        const v = u.cross(this.normal1);
        // get 3 points in the plane:
        const p1 = m4.transformPoint(this.anchor), p2 = m4.transformPoint(this.anchor.plus(v)), p3 = m4.transformPoint(this.anchor.plus(u));
        // and create a new plane from the transformed points:
        return P3.throughPoints(p1, !mirror ? p2 : p3, !mirror ? p3 : p2);
    }
    distanceToLine(line) {
        ts3dutils.assertInst(L3$1, line);
        if (!this.isParallelToLine(line)) {
            return this.distanceToPoint(line.anchor);
        }
        else {
            return 0;
        }
    }
    containsPoint(x) {
        ts3dutils.assertVectors(x);
        return ts3dutils.eq(this.w, this.normal1.dot(x));
    }
    containsLine(line) {
        ts3dutils.assertInst(L3$1, line);
        return this.containsPoint(line.anchor) && this.isParallelToLine(line);
    }
    distanceToPointSigned(point) {
        ts3dutils.assertInst(ts3dutils.V3, point);
        return this.normal1.dot(point) - this.w;
    }
    distanceToPoint(point) {
        ts3dutils.assertInst(ts3dutils.V3, point);
        return Math.abs(this.normal1.dot(point) - this.w);
    }
    intersectionWithLine(line) {
        return line.intersectionWithPlane(this);
    }
    intersectionWithPlane(plane) {
        /*

         this: n0 * x = w0
         plane: n1 * x = w1
         plane perpendicular to both which goes through origin:
         n2 := n0 X x1
         n2 * x = 0
         */
        ts3dutils.assertInst(P3, plane);
        ts3dutils.assert(!this.isParallelToPlane(plane), '!this.isParallelToPlane(plane)');
        /*
         var n0 = this.normal1, n1 = plane.normal1, n2 = n0.cross(n1).unit(), m = M4.forSys(n0, n1, n2)
         var x0 = this.anchor, x1 = plane.anchor, x2 = V3.O
         var p = n2.times(x2.dot(n2))
         .plus(n1.cross(n2).times(x0.dot(n0)))
         .plus(n2.cross(n0).times(x1.dot(n1)))
         .div(m.determinant())
         */
        const n0 = this.normal1, n1 = plane.normal1, n2 = n0.cross(n1).unit();
        const p = ts3dutils.M4.forRows(n0, n1, n2).inversed().transformVector(new ts3dutils.V3(this.w, plane.w, 0));
        return new L3$1(p, n2);
    }
    /**
     * Returns the point in the plane closest to the given point
     *
     */
    projectedPoint(x) {
        // See http://math.stackexchange.com/questions/444968/project-a-point-in-3d-on-a-given-plane
        // p = x - ((x - planeAnchor) * normal1) * normal1
        return x.minus(this.normal1.times(x.minus(this.anchor).dot(this.normal1)));
    }
    projectedVector(x) {
        // See V3.rejectedFrom. Simplified, as this.normal1.length() == 1
        return x.minus(this.normal1.times(x.dot(this.normal1)));
    }
    flipped() {
        return new P3(this.normal1.negated(), -this.w);
    }
    containsCurve(curve) {
        if (curve instanceof L3$1) {
            return this.containsLine(curve);
        }
        else if (curve instanceof SemiEllipseCurve ||
            curve instanceof EllipseCurve ||
            curve instanceof HyperbolaCurve ||
            curve instanceof ParabolaCurve) {
            return this.containsPoint(curve.center) && this.normal1.isParallelTo(curve.normal);
        }
        else if (curve instanceof BezierCurve) {
            return curve.points.every(p => this.containsPoint(p));
        }
        else {
            throw new Error('' + curve);
        }
    }
    hashCode() {
        return this.normal1.hashCode() * 31 | 0 + ts3dutils.floatHashCode(this.w);
    }
}
P3.YZ = new P3(ts3dutils.V3.X, 0);
P3.ZX = new P3(ts3dutils.V3.Y, 0);
P3.XY = new P3(ts3dutils.V3.Z, 0);

const { ceil: ceil$5, floor: floor$5 } = Math;
class Surface extends ts3dutils.Transformable {
    static loopContainsPointGeneral(loop, p, testLine, lineOut) {
        const testPlane = P3.normalOnAnchor(lineOut, p);
        // edges colinear to the testing line; these will always be counted as "inside" relative to the testing line
        const colinearEdges = loop.map((edge) => edge.colinearToLine(testLine));
        let inside = false;
        function logIS(isP) {
            const isT = testLine.pointT(isP);
            if (ts3dutils.eq0(isT)) {
                return true;
            }
            else if (isT > 0) {
                inside = !inside;
            }
        }
        for (let edgeIndex = 0; edgeIndex < loop.length; edgeIndex++) {
            const edge = loop[edgeIndex];
            const nextEdgeIndex = (edgeIndex + 1) % loop.length, nextEdge = loop[nextEdgeIndex];
            //console.log(edge.toSource()) {p:V(2, -2.102, 0),
            if (colinearEdges[edgeIndex]) {
                const lineAT = testLine.pointT(edge.a), lineBT = testLine.pointT(edge.b);
                if (Math.min(lineAT, lineBT) <= ts3dutils.NLA_PRECISION && -ts3dutils.NLA_PRECISION <= Math.max(lineAT, lineBT)) {
                    return exports.PointVsFace.ON_EDGE;
                }
                // edge colinear to intersection
                const nextInside = colinearEdges[nextEdgeIndex] || dotCurve(lineOut, nextEdge.aDir, nextEdge.aDDT) < 0;
                if (!nextInside) {
                    if (logIS(edge.b))
                        return exports.PointVsFace.ON_EDGE;
                }
            }
            else {
                for (const edgeT of edge.edgeISTsWithPlane(testPlane)) {
                    if (edgeT == edge.bT) {
                        if (!testLine.containsPoint(edge.b))
                            continue;
                        // endpoint lies on intersection line
                        if (edge.b.like(p)) {
                            // TODO: refactor, dont check for different sides, just logIs everything
                            return exports.PointVsFace.ON_EDGE;
                        }
                        const edgeInside = dotCurve(lineOut, edge.bDir, edge.bDDT) > 0;
                        const nextInside = colinearEdges[nextEdgeIndex] || dotCurve(lineOut, nextEdge.aDir, nextEdge.aDDT) < 0;
                        if (edgeInside != nextInside) {
                            if (logIS(edge.b))
                                return exports.PointVsFace.ON_EDGE;
                        }
                    }
                    else if (edgeT != edge.aT) {
                        const p = edge.curve.at(edgeT);
                        if (!testLine.containsPoint(p))
                            continue;
                        // edge crosses line, neither starts nor ends on it
                        if (logIS(p))
                            return exports.PointVsFace.ON_EDGE;
                        // TODO: tangents?
                    }
                }
            }
        }
        return inside ? exports.PointVsFace.INSIDE : exports.PointVsFace.OUTSIDE;
    }
    toString() {
        return this.toSource();
    }
    toSource(rounder = x => x) {
        return ts3dutils.callsce.call(undefined, 'new ' + this.constructor.name, ...this.getConstructorParameters());
    }
    isCurvesWithSurface(surface) {
        return surface.isCurvesWithSurface(this).map(curve => curve.reversed());
    }
    containsCurve(curve) {
        if (curve instanceof ImplicitCurve) {
            for (let i = ceil$5(curve.tMin); i <= floor$5(curve.tMax); i++) {
                if (!this.containsPoint(curve.points[i])) {
                    return false;
                }
            }
            return true;
        }
        else {
            return false;
        }
    }
    flipped2(doFlip) {
        return doFlip ? this.flipped() : this;
    }
    clipCurves(curves) {
        return curves;
    }
    hashCode() {
        return this.getConstructorParameters().hashCode();
    }
    zDirVolume(allEdges) {
        return this.visit(ZDirVolumeVisitor, allEdges);
    }
    calculateArea(allEdges) {
        return this.visit(CalculateAreaVisitor, allEdges);
    }
}

(function (PointVsFace) {
    PointVsFace[PointVsFace["INSIDE"] = 0] = "INSIDE";
    PointVsFace[PointVsFace["OUTSIDE"] = 1] = "OUTSIDE";
    PointVsFace[PointVsFace["ON_EDGE"] = 2] = "ON_EDGE";
})(exports.PointVsFace || (exports.PointVsFace = {}));

const { ceil: ceil$6, min: min$4 } = Math;
class ParametricSurface extends Surface {
    static isCurvesParametricImplicitSurface(ps, is, sStep, tStep = sStep, curveStepSize) {
        const pf = ps.pSTFunc(), icc = is.implicitFunction();
        const dpds = ps.dpds();
        const dpdt = ps.dpdt();
        const didp = is.didp.bind(is);
        const ist = (x, y) => icc(pf(x, y));
        const dids = (s, t) => didp(pf(s, t)).dot(dpds(s, t));
        const didt = (s, t) => didp(pf(s, t)).dot(dpdt(s, t));
        const mf = exports.MathFunctionR2R.forFFxFy(ist, dids, didt);
        const curves = Curve.breakDownIC(mf, ps, sStep, tStep, curveStepSize, dids, didt)
            .map(({ points, tangents }, i) => PICurve$1.forParametricPointsTangents(ps, is, points, tangents, curveStepSize));
        return curves;
    }
    static is(obj) {
        return obj.pSTFunc;
    }
    pST(s, t) {
        return this.pSTFunc()(s, t);
    }
    pSTFunc() {
        return this.pST.bind(this);
    }
    stP(pWC) {
        return this.stPFunc()(pWC);
    }
    stPFunc() {
        return this.stP.bind(this);
    }
    bounds(s, t) {
        return this.sMin <= s && s <= this.sMax && this.tMin <= t && t <= this.tMax;
    }
    /**
     * Positive values are inside bounds.
     */
    boundsSigned(s, t) {
        return min$4(s - this.sMin, this.sMax - s, t - this.tMin, this.tMax - t);
    }
    normalP(p) {
        const pmPoint = this.stPFunc()(p);
        return this.normalST(pmPoint.x, pmPoint.y);
    }
    normalSTFunc() {
        return this.normalST.bind(this);
    }
    normalST(s, t) {
        return this.normalSTFunc()(s, t);
    }
    parametersValid(s, t) {
        return ts3dutils.between(s, this.sMin, this.sMax) && ts3dutils.between(t, this.tMin, this.tMax);
    }
    pointFoot(pWC, ss, st) {
        throw new Error();
    }
    toMesh() {
        ts3dutils.assert(isFinite(this.tMin) && isFinite(this.tMax) && isFinite(this.sMin) && isFinite(this.sMax));
        return tsgl.Mesh.parametric(this.pSTFunc(), this.normalSTFunc(), this.sMin, this.sMax, this.tMin, this.tMax, ceil$6((this.sMax - this.sMin) / this.uStep), ceil$6((this.tMax - this.tMin) / this.vStep));
    }
    isCurvesWithImplicitSurface(is, sStep, tStep, stepSize) {
        return ParametricSurface.isCurvesParametricImplicitSurface(this, is, sStep, tStep, stepSize);
    }
}
class ImplicitSurface extends Surface {
    static is(obj) {
        return obj.implicitFunction;
    }
}

const { PI: PI$5, cos: cos$3, sin: sin$3, min: min$5, max: max$4, tan: tan$2, ceil: ceil$7, floor: floor$6, abs: abs$6, sqrt: sqrt$2, pow: pow$2, atan2: atan2$2, round: round$2, sign: sign$2 } = Math;
class ConicSurface extends ParametricSurface {
    /**
     * returns new cone C = {apex + f1 * z * cos(d) + f2 * z * sin(d) + f3 * z | -PI <= d <= PI, 0 <= z}
     * @param f1
     * @param f2
     * @param dir Direction in which the cone opens. The ellipse spanned by f1, f2 is contained at (apex + f1).
     */
    constructor(center, f1, f2, dir) {
        super();
        this.center = center;
        this.f1 = f1;
        this.f2 = f2;
        this.dir = dir;
        ts3dutils.assertVectors(center, f1, f2, dir);
        this.matrix = ts3dutils.M4.forSys(f1, f2, dir, center);
        this.inverseMatrix = this.matrix.inversed();
        this.normalDir = sign$2(this.f1.cross(this.f2).dot(this.dir));
        this.normalMatrix = this.matrix.as3x3().inversed().transposed().scale(this.normalDir);
    }
    get apex() {
        return this.center;
    }
    static atApexThroughEllipse(apex, ellipse) {
        ts3dutils.assertVectors(apex);
        ts3dutils.assertInst(SemiEllipseCurve, ellipse);
        return new ConicSurface(apex, ellipse.f1, ellipse.f2, apex.to(ellipse.center));
    }
    static unitISLineTs(anchor, dir) {
        const { x: ax, y: ay, z: az } = anchor;
        const { x: dx, y: dy, z: dz } = dir;
        // this cone: x² + y² = z²
        // line: p = anchor + t * dir1
        // split line equation into 3 component equations, insert into cone equation
        // transform to form (a t² + b t + c = 0) and solve with pqFormula
        const a = dx * dx + dy * dy - dz * dz;
        const b = 2 * (ax * dx + ay * dy - az * dz);
        const c = ax * ax + ay * ay - az * az;
        // cone only defined for 0 <= z, so filter invalid values
        return ts3dutils.pqFormula(b / a, c / a).filter(t => 0 < az + t * dz);
    }
    // calculate intersection of plane ax + cz = d and cone x² + y² = z²
    static unitISPlane(a, c, d) {
        if (ts3dutils.eq0(c)) {
            // plane is "vertical", i.e. parallel to Y and Z axes
            ts3dutils.assert(!ts3dutils.eq0(a)); // normal would be zero, which is invalid
            // z² - y² = d²/a²
            if (ts3dutils.eq0(d)) {
                // d = 0 => z² - y² = 0 => z² = y² => z = y
                // plane goes through origin/V3.O
                return [new L3$1(ts3dutils.V3.O, new ts3dutils.V3(0, -sqrt$2(2) / 2, -sqrt$2(2) / 2), undefined, 0),
                    new L3$1(ts3dutils.V3.O, new ts3dutils.V3(0, -sqrt$2(2) / 2, sqrt$2(2) / 2), 0)];
            }
            else {
                // hyperbola
                const center = new ts3dutils.V3(d / a, 0, 0);
                const f1 = new ts3dutils.V3(0, 0, abs$6(d / a)); // abs, because we always want the hyperbola to be pointing up
                const f2 = new ts3dutils.V3(0, d / a, 0);
                return [new HyperbolaCurve(center, f1, f2)];
            }
        }
        else {
            // c != 0
            const aa = a * a, cc = c * c;
            if (ts3dutils.eq0(d)) {
                // ax + cz = d => x = d - cz / a => x² = d² - 2cdz/a + c²z²/a²
                // x² + y² = z²
                // => d² - 2cdz/a + c²z²/a² + y² = z²
                if (ts3dutils.eq(aa, cc)) {
                    return [new L3$1(ts3dutils.V3.O, new ts3dutils.V3(c, 0, -a).unit())];
                }
                else if (aa < cc) {
                    ts3dutils.assert(false, 'intersection is single point V3.O');
                }
                else if (aa > cc) {
                    return [new L3$1(ts3dutils.V3.O, new ts3dutils.V3(c, sqrt$2(aa - cc), -a).unit()),
                        new L3$1(ts3dutils.V3.O, new ts3dutils.V3(c, -sqrt$2(aa - cc), -a).unit())];
                }
            }
            else {
                if (ts3dutils.eq(aa, cc)) {
                    // parabola
                    const parabolaVertex = new ts3dutils.V3(d / 2 / a, 0, d / 2 / c);
                    const parabolaVertexTangentPoint = new ts3dutils.V3(d / 2 / a, d / c, d / 2 / c);
                    const p2 = new ts3dutils.V3(0, 0, d / c);
                    const f2 = p2.minus(parabolaVertex);
                    return [new ParabolaCurve(parabolaVertex, parabolaVertexTangentPoint.minus(parabolaVertex), f2.z < 0
                            ? f2.negated()
                            : f2)];
                }
                else if (aa < cc) {
                    // ellipse
                    const center = new ts3dutils.V3(-a * d / (cc - aa), 0, d * c / (cc - aa));
                    if (center.z < 0) {
                        return [];
                    }
                    const p1 = new ts3dutils.V3(d / (a - c), 0, -d / (a - c));
                    const p2 = new ts3dutils.V3(-a * d / (cc - aa), d / sqrt$2(cc - aa), d * c / (cc - aa));
                    return [new EllipseCurve(center, center.to(p1), center.to(p2))];
                }
                else if (aa > cc) {
                    // hyperbola
                    const center = new ts3dutils.V3(-a * d / (cc - aa), 0, d * c / (cc - aa));
                    const p1 = new ts3dutils.V3(d / (a - c), 0, -d / (a - c));
                    const p2 = new ts3dutils.V3(-a * d / (cc - aa), d / sqrt$2(aa - cc), d * c / (cc - aa));
                    const f1 = center.to(p1);
                    return [new HyperbolaCurve(center, f1.z > 0 ? f1 : f1.negated(), center.to(p2))];
                }
            }
        }
    }
    equals(obj) {
        return this == obj ||
            Object.getPrototypeOf(this) == Object.getPrototypeOf(obj)
                && this.center.equals(obj.center)
                && this.f1.equals(obj.f1)
                && this.f2.equals(obj.f2)
                && this.dir.equals(obj.dir);
    }
    like(object) {
        if (!this.isCoplanarTo(object))
            return false;
        // normals need to point in the same direction (outwards or inwards) for both
        return this.normalDir == object.normalDir;
    }
    getVectors() {
        return [{ anchor: this.center, dir1: this.dir },
            { anchor: this.center.plus(this.dir), dir1: this.f1 },
            { anchor: this.center.plus(this.dir), dir1: this.f2 }];
    }
    getSeamPlane() {
        return P3.forAnchorAndPlaneVectors(this.center, this.f1, this.dir);
    }
    loopContainsPoint(contour, p) {
        ts3dutils.assertVectors(p);
        const line = this.center.like(p)
            ? new L3$1(p, this.matrix.transformVector(new ts3dutils.V3(0, 1, 1)).unit())
            : L3$1.throughPoints(p, this.apex);
        const lineOut = line.dir1.cross(this.dir);
        return Surface.loopContainsPointGeneral(contour, p, line, lineOut);
    }
    getConstructorParameters() {
        return [this.center, this.f1, this.f2, this.dir];
    }
    isTsForLine(line) {
        // transforming line manually has advantage that dir1 will not be renormalized,
        // meaning that calculated values t for lineLC are directly transferable to line
        const anchorLC = this.inverseMatrix.transformPoint(line.anchor);
        const dirLC = this.inverseMatrix.transformVector(line.dir1);
        return ConicSurface.unitISLineTs(anchorLC, dirLC);
    }
    /**
     * Interestingly, two cones don't need to have parallel dirs to be coplanar.
     */
    isCoplanarTo(surface) {
        if (this === surface)
            return true;
        if (!(surface instanceof ConicSurface) || !this.apex.like(surface.apex))
            return false;
        // at this point apexes are equal
        return this.containsEllipse(new SemiEllipseCurve(surface.center.plus(surface.dir), surface.f1, surface.f2));
    }
    containsEllipse(ellipse) {
        const ellipseLC = ellipse.transform(this.inverseMatrix);
        if (ellipseLC.center.z < 0) {
            return false;
        }
        const { f1, f2 } = ellipseLC.rightAngled();
        const p1 = ellipseLC.center.plus(f1), p2 = ellipseLC.center.plus(f2);
        // check if both endpoints are on the cone's surface
        // and that one main axis is perpendicular to the Z-axis
        return ts3dutils.eq(Math.pow(p1.x, 2) + Math.pow(p1.y, 2), Math.pow(p1.z, 2))
            && ts3dutils.eq(Math.pow(p2.x, 2) + Math.pow(p2.y, 2), Math.pow(p2.z, 2))
            && (ts3dutils.eq0(f1.z) || ts3dutils.eq0(f2.z));
    }
    containsLine(line) {
        const lineLC = line.transform(this.inverseMatrix);
        const d = lineLC.dir1;
        return lineLC.containsPoint(ts3dutils.V3.O) && ts3dutils.eq(d.x * d.x + d.y * d.y, d.z * d.z);
    }
    containsParabola(curve) {
        ts3dutils.assertInst(ParabolaCurve, curve);
        const curveLC = curve.transform(this.inverseMatrix);
        if (curveLC.center.z < 0 || curveLC.f2.z < 0) {
            return false;
        }
        const { center, f1, f2 } = curveLC.rightAngled();
        // check if center is on the surface,
        // that tangent is perpendicular to the Z-axis
        // and that "y" axis is parallel to surface
        return ts3dutils.eq(center.x * center.x + center.y * center.y, center.z * center.z)
            && ts3dutils.eq0(f1.z)
            && ts3dutils.eq(f2.x * f2.x + f2.y * f2.y, f2.z * f2.z);
    }
    containsHyperbola(curve) {
        ts3dutils.assertInst(HyperbolaCurve, curve);
        return true;
        const curveLC = curve.transform(this.inverseMatrix);
        if (curveLC.center.z < 0 || curveLC.f2.z < 0) {
            return false;
        }
        const { center, f1, f2 } = curveLC.rightAngled();
        // check if center is on the surface,
        // that tangent is perpendicular to the Z-axis
        return true;
        return ts3dutils.eq(center.x * center.x + center.y * center.y, center.z * center.z)
            && ts3dutils.eq0(f1.z);
    }
    containsCurve(curve) {
        if (curve instanceof SemiEllipseCurve) {
            return this.containsEllipse(curve);
        }
        else if (curve instanceof L3$1) {
            return this.containsLine(curve);
        }
        else if (curve instanceof HyperbolaCurve) {
            return this.containsHyperbola(curve);
        }
        else if (curve instanceof ParabolaCurve) {
            return this.containsParabola(curve);
        }
        else {
            return super.containsCurve(curve);
        }
    }
    transform(m4) {
        return new ConicSurface(m4.transformPoint(this.center), m4.transformVector(this.f1).times(m4.isMirroring() ? -1 : 1), m4.transformVector(this.f2), m4.transformVector(this.dir));
    }
    rightAngled() {
        // TODO
    }
    flipped() {
        return new ConicSurface(this.center, this.f1.negated(), this.f2, this.dir);
    }
    normalSTFunc() {
        const { f1, f2 } = this, f3 = this.dir;
        return (d, z) => {
            return f2.cross(f1).plus(f2.cross(f3.times(Math.cos(d)))).plus(f3.cross(f1.times(Math.sin(d)))).unit();
        };
    }
    normalP(p) {
        //TODO assert(!p.like(this.center))
        const pLC = this.inverseMatrix.transformPoint(p);
        return this.normalSTFunc()(pLC.angleXY(), pLC.z);
    }
    pSTFunc() {
        return (s, t) => {
            // center + f1 t cos s + f2 t sin s + t dir
            return this.matrix.transformPoint(new ts3dutils.V3(t * cos$3(s), t * sin$3(s), t));
        };
    }
    dpds() {
        return (s, t) => {
            const resultLC = new ts3dutils.V3(t * -sin$3(s), t * cos$3(s), 0);
            return this.matrix.transformVector(resultLC);
        };
    }
    dpdt() {
        return (s, t) => {
            const resultLC = new ts3dutils.V3(cos$3(s), sin$3(s), 1);
            return this.matrix.transformVector(resultLC);
        };
    }
    implicitFunction() {
        return pWC => {
            const pLC = this.inverseMatrix.transformPoint(pWC);
            const radiusLC = pLC.lengthXY();
            return this.normalDir * (radiusLC - pLC.z);
        };
    }
    containsPoint(p) {
        return ts3dutils.eq0(this.implicitFunction()(p));
    }
    boundsFunction() {
        ts3dutils.assert(false);
    }
    stP(pWC) {
        const pLC = this.inverseMatrix.transformPoint(pWC);
        const angle = pLC.angleXY();
        return new ts3dutils.V3(angle < -PI$5 / 2 ? angle + ts3dutils.TAU : angle, pLC.z, 0);
    }
    isCurvesWithSurface(surface) {
        if (surface instanceof PlaneSurface$1) {
            return this.isCurvesWithPlane(surface.plane);
        }
        else if (ImplicitSurface.is(surface)) {
            return ParametricSurface.isCurvesParametricImplicitSurface(this, surface, 0.1, 0.1 / this.dir.length(), 0.02);
        }
        return super.isCurvesWithSurface(surface);
    }
    getCenterLine() {
        return new L3$1(this.center, this.dir);
    }
    isCurvesWithPlane(plane) {
        ts3dutils.assertInst(P3, plane);
        const planeLC = plane.transform(this.inverseMatrix);
        const planeNormal = planeLC.normal1;
        const c = planeNormal.z;
        /** "rotate" plane normal1 when passing to {@link ConicSurface.unitISPlane} so that
         *  y-component of normal1 is 0 */
        const a = planeNormal.lengthXY();
        const d = planeLC.w;
        // generated curves need to be rotated back before transforming to world coordinates
        const rotationMatrix = ts3dutils.M4.rotateZ(planeNormal.angleXY());
        const wcMatrix = ts3dutils.eq0(planeNormal.lengthXY())
            ? this.matrix
            : this.matrix.times(rotationMatrix);
        return ConicSurface.unitISPlane(a, c, d).flatMap(curve => {
            const curveWC = curve.transform(wcMatrix);
            if (curve instanceof EllipseCurve) {
                const curveLC = curve.transform(rotationMatrix);
                const ts = curveLC.isTsWithPlane(P3.ZX);
                const intervals = ts3dutils.getIntervals(ts, -PI$5, PI$5).filter(([a, b]) => curveLC.at((a + b) / 2).y > 0);
                return intervals.flatMap(([a, b]) => SemiEllipseCurve.fromEllipse(curveWC, a, b));
            }
            const p = curveWC.at(0.2);
            return this.normalP(p).cross(plane.normal1).dot(curveWC.tangentAt(0.2)) > 0
                ? curveWC : curveWC.reversed();
        });
    }
    edgeLoopCCW(contour) {
        const ptpF = this.stPFunc();
        return ts3dutils.isCCW(contour.flatMap(e => e.getVerticesNo0()).map(v => ptpF(v)), ts3dutils.V3.Z);
    }
}
/**
 * Unit cone. x² + y² = z², 0 <= z
 */
ConicSurface.UNIT = new ConicSurface(ts3dutils.V3.O, ts3dutils.V3.X, ts3dutils.V3.Y, ts3dutils.V3.Z);
ConicSurface.prototype.uStep = PI$5 / 16;
ConicSurface.prototype.vStep = 256;
ConicSurface.prototype.sMin = 0;
ConicSurface.prototype.sMax = PI$5;
ConicSurface.prototype.tMin = 0;
ConicSurface.prototype.tMax = 16;

const { PI: PI$6, cos: cos$4, sin: sin$4, abs: abs$7, sign: sign$3 } = Math;
class EllipsoidSurface extends ParametricSurface {
    constructor(center, f1, f2, f3) {
        super();
        this.center = center;
        this.f1 = f1;
        this.f2 = f2;
        this.f3 = f3;
        ts3dutils.assertVectors(center, f1, f2, f3);
        this.matrix = ts3dutils.M4.forSys(f1, f2, f3, center);
        this.inverseMatrix = this.matrix.inversed();
        this.normalDir = sign$3(this.f1.cross(this.f2).dot(this.f3));
        this.pLCNormalWCMatrix = this.matrix.as3x3().inversed().transposed().scale(this.normalDir);
        this.pWCNormalWCMatrix = this.pLCNormalWCMatrix.times(this.inverseMatrix);
    }
    /**
     * unit sphere: x² + y² + z² = 1
     * line: p = anchor + t * dir |^2
     * p² = (anchor + t * dir)^2
     * 1 == (anchor + t * dir)^2
     * 1 == anchor DOT anchor + 2 * anchor * t * dir + t² * dir DOT dir
     */
    static unitISTsWithLine(anchor, dir) {
        // for 0 = a t² + b t + c
        const a = dir.dot(dir);
        const b = 2 * anchor.dot(dir);
        const c = anchor.dot(anchor) - 1;
        return ts3dutils.pqFormula(b / a, c / a);
    }
    /**
     * unit sphere: x² + y² + z² = 1
     * plane: normal1 DOT p = w
     */
    static unitISCurvesWithPlane(plane) {
        ts3dutils.assertInst(P3, plane);
        let distPlaneCenter = Math.abs(plane.w);
        if (ts3dutils.lt(distPlaneCenter, 1)) {
            // result is a circle
            // radius of circle: imagine right angled triangle (origin -> center of intersection circle -> point on
            // intersection circle) pythagoras: 1² == distPlaneCenter² + isCircleRadius² => isCircleRadius == sqrt(1 -
            // distPlaneCenter²)
            const isCircleRadius = Math.sqrt(1 - distPlaneCenter * distPlaneCenter);
            const center = plane.anchor;
            const f1 = plane.normal1.getPerpendicular().toLength(isCircleRadius);
            const f2 = plane.normal1.cross(f1);
            return [new EllipseCurve(plane.anchor, f1, f2)];
        }
        else {
            return [];
        }
    }
    static sphere(radius, center) {
        ts3dutils.assertNumbers(radius);
        center && ts3dutils.assertVectors(center);
        return new EllipsoidSurface(center || ts3dutils.V3.O, new ts3dutils.V3(radius, 0, 0), new ts3dutils.V3(0, radius, 0), new ts3dutils.V3(0, 0, radius));
    }
    /**
     * x²/a² + y²/b² + z²/c² = 1
     */
    static forABC(a, b, c, center) {
        return new EllipsoidSurface(center || ts3dutils.V3.O, new ts3dutils.V3(a, 0, 0), new ts3dutils.V3(0, b, 0), new ts3dutils.V3(0, 0, c));
    }
    static calculateAreaSpheroid(a, b, c, edges) {
        ts3dutils.assertf(() => a.isPerpendicularTo(b));
        ts3dutils.assertf(() => b.isPerpendicularTo(c));
        ts3dutils.assertf(() => c.isPerpendicularTo(a));
        // handling discontinuities:
        // option 1: check for intersections with baseline, if there are any integrate parts separetely
        // "rotate" the edge so that there are no overlaps
        const matrix = ts3dutils.M4.forSys(a, b, c), inverseMatrix = matrix.inversed();
        const circleRadius = a.length();
        const c1 = c.unit();
        const totalArea = edges.map(edge => {
            if (edge.curve instanceof EllipseCurve) {
                const f = (t) => {
                    const at = edge.curve.at(t), tangent = edge.tangentAt(t);
                    const localAt = inverseMatrix.transformPoint(at);
                    const angleXY = localAt.angleXY();
                    const arcLength = angleXY * circleRadius * Math.sqrt(1 + Math.pow(localAt.z, 2));
                    const scaling = Math.sqrt(1 + Math.pow(c1.dot(tangent), 2));
                    return arcLength * scaling;
                };
                const val = ts3dutils.glqInSteps(f, edge.aT, edge.bT, 1);
                console.log('edge', edge, val);
                return val;
            }
            else {
                ts3dutils.assertNever();
            }
        }).sum();
        return totalArea;
    }
    like(obj) {
        return this.isCoplanarTo(obj) && this.isInsideOut() == obj.isInsideOut();
    }
    edgeLoopCCW(loop) {
        throw new Error();
    }
    rootPoints() {
    }
    getConstructorParameters() {
        return [this.center, this.f1, this.f2, this.f3];
    }
    equals(obj) {
        return this == obj ||
            Object.getPrototypeOf(obj) == this.constructor.prototype
                && this.matrix.equals(obj.matrix);
    }
    isCurvesWithPlane(plane) {
        const planeLC = plane.transform(this.inverseMatrix);
        return EllipsoidSurface.unitISCurvesWithPlane(planeLC).map(c => c.transform(this.matrix));
    }
    isCurvesWithSurface(surface) {
        if (surface instanceof PlaneSurface$1) {
            return this.isCurvesWithPlane(surface.plane);
        }
        else if (surface instanceof CylinderSurface) {
            if (surface.dir.isParallelTo(this.dir1)) {
                const ellipseProjected = surface.baseCurve.transform(ts3dutils.M4.project(this.baseEllipse.getPlane(), this.dir1));
                return this.baseEllipse.isInfosWithEllipse(ellipseProjected).map(info => new L3$1(info.p, this.dir1));
            }
            else if (ts3dutils.eq0(this.getCenterLine().distanceToLine(surface.getCenterLine()))) {
                ts3dutils.assert(false);
            }
            else {
                ts3dutils.assert(false);
            }
        }
        else if (surface instanceof ProjectedCurveSurface) {
            const surfaceLC = surface.transform(this.inverseMatrix);
            const baseCurveLC = surfaceLC.baseCurve.project(new P3(surfaceLC.dir, 0));
            const ists = baseCurveLC.isTsWithSurface(EllipsoidSurface.UNIT);
            const insideIntervals = iii(ists, EllipsoidSurface.UNIT, baseCurveLC);
            const curves = insideIntervals.flatMap(ii => {
                const aLine = new L3$1(baseCurveLC.at(ii[0]), surfaceLC.dir);
                const a = EllipsoidSurface.UNIT.isTsForLine(aLine).map(t => aLine.at(t));
                const bLine = new L3$1(baseCurveLC.at(ii[1]), surfaceLC.dir);
                const b = EllipsoidSurface.UNIT.isTsForLine(bLine).map(t => bLine.at(t));
                return [0, 1].map(i => {
                    let aP = a[i] || a[0], bP = b[i] || b[0];
                    0 !== i && ([aP, bP] = [bP, aP]);
                    ts3dutils.assert(EllipsoidSurface.UNIT.containsPoint(aP));
                    ts3dutils.assert(EllipsoidSurface.UNIT.containsPoint(bP));
                    return PICurve.forStartEnd(surface, this.asEllipsoidSurface(), aP, bP);
                });
            });
            return curves;
        }
    }
    isTsForLine(line) {
        ts3dutils.assertInst(L3$1, line);
        // transforming line manually has advantage that dir1 will not be renormalized,
        // meaning that calculated values t for localLine are directly transferable to line
        const anchorLC = this.inverseMatrix.transformPoint(line.anchor);
        const dirLC = this.inverseMatrix.transformVector(line.dir1);
        return EllipsoidSurface.unitISTsWithLine(anchorLC, dirLC);
    }
    isCoplanarTo(surface) {
        if (this === surface)
            return true;
        if (surface.constructor !== EllipsoidSurface)
            return false;
        if (!this.center.like(surface.center))
            return false;
        if (this.isSphere())
            return surface.isSphere() && ts3dutils.eq(this.f1.length(), this.f2.length());
        const localOtherMatrix = this.inverseMatrix.times(surface.matrix);
        // Ellipsoid with matrix localOtherMatrix is unit sphere iff localOtherMatrix is orthogonal
        return localOtherMatrix.is3x3() && localOtherMatrix.isOrthogonal();
    }
    containsEllipse(ellipse) {
        const localEllipse = ellipse.transform(this.inverseMatrix);
        const distLocalEllipseCenter = localEllipse.center.length();
        const correctRadius = Math.sqrt(1 - distLocalEllipseCenter * distLocalEllipseCenter);
        return ts3dutils.lt(distLocalEllipseCenter, 1) && localEllipse.isCircular() && localEllipse.f1.hasLength(correctRadius);
    }
    containsCurve(curve) {
        if (curve instanceof EllipseCurve) {
            return this.containsEllipse(curve);
        }
        else {
            return super.containsCurve(curve);
        }
    }
    transform(m4) {
        return new EllipsoidSurface(m4.transformPoint(this.center), m4.transformVector(this.f1), m4.transformVector(this.f2), m4.transformVector(this.f3));
    }
    isInsideOut() {
        return this.f1.cross(this.f2).dot(this.f3) < 0;
    }
    flipped() {
        return new EllipsoidSurface(this.center, this.f1, this.f2, this.f3.negated());
    }
    toMesh(subdivisions = 3) {
        return tsgl.Mesh.sphere(subdivisions).transform(this.matrix);
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
            let { f1, f2, f3 } = this;
            let normal = f2.cross(f3).times(Math.cos(b) * Math.cos(a))
                .plus(f3.cross(f1).times(Math.cos(b) * Math.sin(a)))
                .plus(f1.cross(f2).times(Math.sin(b)))
                .unit();
            return normal;
        };
    }
    normalP(p) {
        return this.normalMatrix.transformVector(this.inverseMatrix.transformPoint(p)).unit();
    }
    normalST(s, t) {
        return this.normalMatrix.transformVector(ts3dutils.V3.sphere(s, t));
    }
    pST(s, t) {
        return this.matrix.transformPoint(ts3dutils.V3.sphere(s, t));
    }
    //   d/dp (this.implicitFunction(p)) =
    // = d/dp (this.inverseMatrix.transformPoint(p).length() - 1)
    // = d/dp (this.inverseMatrix.transformPoint(p) * this.inverseMatrix.transformPoint(pWC).unit()
    dpds() {
        return (s, t) => this.matrix.transformVector(new ts3dutils.V3(sin$4(s) * -cos$4(t), cos$4(s) * cos$4(t), 0));
    }
    dpdt() {
        return (s, t) => this.matrix.transformVector(new ts3dutils.V3(sin$4(t) * -cos$4(s), -sin$4(s) * sin$4(t), cos$4(t)));
    }
    stPFunc() {
        return (pWC, hint) => {
            const pLC = this.inverseMatrix.transformPoint(pWC);
            let alpha = pLC.angleXY();
            if (abs$7(alpha) > Math.PI - ts3dutils.NLA_PRECISION) {
                ts3dutils.assert(hint == -PI$6 || hint == PI$6);
                alpha = hint;
            }
            let beta = Math.asin(pLC.z);
            return new ts3dutils.V3(alpha, beta, 0);
        };
    }
    isSphere() {
        return ts3dutils.eq(this.f1.length(), this.f2.length())
            && ts3dutils.eq(this.f2.length(), this.f3.length())
            && ts3dutils.eq(this.f3.length(), this.f1.length())
            && this.f1.isPerpendicularTo(this.f2)
            && this.f2.isPerpendicularTo(this.f3)
            && this.f3.isPerpendicularTo(this.f1);
    }
    isVerticalSpheroid() {
        return ts3dutils.eq(this.f1.length(), this.f2.length())
            && this.f1.isPerpendicularTo(this.f2)
            && this.f2.isPerpendicularTo(this.f3)
            && this.f3.isPerpendicularTo(this.f1);
    }
    implicitFunction() {
        return (pWC) => {
            const pLC = this.inverseMatrix.transformPoint(pWC);
            return pLC.length() - 1;
        };
    }
    // = this.inverseMatrix.transformPoint(this.inverseMatrix.transformPoint(pWC).unit())
    didp(pWC) {
        const pLC = this.inverseMatrix.transformPoint(pWC);
        return this.inverseMatrix.transformVector(pLC.unit());
    }
    mainAxes() {
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
        const { f1, f2, f3 } = this;
        if (ts3dutils.eq0(f1.dot(f2)) && ts3dutils.eq0(f2.dot(f3)) && ts3dutils.eq0(f3.dot(f1))) {
            return this;
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
        const { U, SIGMA } = this.matrix.svd3();
        ts3dutils.assert(SIGMA.isDiagonal());
        ts3dutils.assert(U.isOrthogonal());
        const U_SIGMA = U.times(SIGMA);
        // column vectors of U_SIGMA
        const [mainF1, mainF2, mainF3] = ts3dutils.arrayFromFunction(3, i => new ts3dutils.V3(U_SIGMA.m[i], U_SIGMA.m[i + 4], U_SIGMA.m[i + 8]));
        return new EllipsoidSurface(this.center, mainF1, mainF2, mainF3);
    }
    containsPoint(p) {
        return ts3dutils.eq0(this.implicitFunction()(p));
    }
    boundsFunction() {
        return (a, b) => ts3dutils.between(b, -PI$6, PI$6);
    }
    volume() {
        return 4 / 3 * Math.PI * this.f1.dot(this.f2.cross(this.f3));
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
    loopContainsPoint(loop, p) {
        ts3dutils.assertVectors(p);
        const testLine = new EllipseCurve(this.center, this.matrix.transformVector(this.inverseMatrix.transformPoint(p).withElement('z', 0).unit()), this.f3);
        const pT = testLine.pointT(p);
        const lineOut = testLine.normal;
        const testPlane = P3.normalOnAnchor(testLine.normal, p);
        const colinearEdges = loop.map((edge) => edge.curve.isColinearTo(testLine));
        let inside = false;
        function logIS(isP) {
            const isT = testLine.pointT(isP);
            if (ts3dutils.eq(pT, isT)) {
                return true;
            }
            else if (pT < isT && ts3dutils.le(isT, PI$6)) {
                inside = !inside;
            }
        }
        for (let edgeIndex = 0; edgeIndex < loop.length; edgeIndex++) {
            const edge = loop[edgeIndex];
            const nextEdgeIndex = (edgeIndex + 1) % loop.length, nextEdge = loop[nextEdgeIndex];
            //console.log(edge.toSource()) {p:V(2, -2.102, 0),
            if (colinearEdges[edgeIndex]) {
                const lineAT = testLine.pointT(edge.a), lineBT = testLine.pointT(edge.b);
                if (ts3dutils.le(Math.min(lineAT, lineBT), pT) && ts3dutils.ge(pT, Math.max(lineAT, lineBT))) {
                    return exports.PointVsFace.ON_EDGE;
                }
                // edge colinear to intersection
                const nextInside = colinearEdges[nextEdgeIndex] || dotCurve(lineOut, nextEdge.aDir, nextEdge.aDDT) < 0;
                if (nextInside) {
                    if (logIS(edge.b))
                        return exports.PointVsFace.ON_EDGE;
                }
            }
            else {
                for (const edgeT of edge.edgeISTsWithPlane(testPlane)) {
                    if (edgeT == edge.bT) {
                        if (!testLine.containsPoint(edge.b))
                            continue;
                        // endpoint lies on intersection testLine
                        const edgeInside = dotCurve(lineOut, edge.bDir, edge.bDDT) < 0;
                        const nextInside = colinearEdges[nextEdgeIndex] || dotCurve(lineOut, nextEdge.aDir, nextEdge.aDDT) < 0;
                        if (edgeInside != nextInside) {
                            if (logIS(edge.b))
                                return exports.PointVsFace.ON_EDGE;
                        }
                    }
                    else if (edgeT != edge.aT) {
                        const p = edge.curve.at(edgeT);
                        if (!testLine.containsPoint(p))
                            continue;
                        // edge crosses testLine, neither starts nor ends on it
                        if (logIS(p))
                            return exports.PointVsFace.ON_EDGE;
                        // TODO: tangents?
                    }
                }
            }
        }
        return inside ? exports.PointVsFace.INSIDE : exports.PointVsFace.OUTSIDE;
    }
    surfaceAreaApprox() {
        // See https://en.wikipedia.org/wiki/Ellipsoid#Surface_area
        const mainAxes = this.mainAxes(), a = mainAxes.f1.length(), b = mainAxes.f2.length(), c = mainAxes.f3.length();
        const p = 1.6075;
        return 4 * PI$6 * Math.pow((Math.pow(a * b, p) + Math.pow(b * c, p) + Math.pow(c * a, p)) / 3, 1 / p);
    }
    surfaceArea() {
        // See https://en.wikipedia.org/wiki/Ellipsoid#Surface_area
        const mainAxes = this.mainAxes(), f1l = mainAxes.f1.length(), f2l = mainAxes.f2.length(), f3l = mainAxes.f3.length(), [c, b, a] = [f1l, f2l, f3l].sort(ts3dutils.MINUS);
        // https://en.wikipedia.org/w/index.php?title=Spheroid&oldid=761246800#Area
        function spheroidArea(a, c) {
            if (c < a) {
                const eccentricity2 = 1 - Math.pow(c, 2) / Math.pow(a, 2);
                const eccentricity = Math.sqrt(eccentricity2);
                return 2 * PI$6 * Math.pow(a, 2) * (1 + (1 - eccentricity2) / Math.sqrt(eccentricity) * Math.atanh(eccentricity));
            }
            else {
                const eccentricity = Math.sqrt(1 - Math.pow(a, 2) / Math.pow(c, 2));
                return 2 * PI$6 * Math.pow(a, 2) * (1 + c / a / eccentricity * Math.asin(eccentricity));
            }
        }
        if (ts3dutils.eq(a, b)) {
            return spheroidArea(a, c);
        }
        else if (ts3dutils.eq(b, c)) {
            return spheroidArea(b, a);
        }
        else if (ts3dutils.eq(c, a)) {
            return spheroidArea(c, b);
        }
        const phi = Math.acos(c / a);
        const k2 = Math.pow(a, 2) * (Math.pow(b, 2) - Math.pow(c, 2)) / (Math.pow(b, 2) * (Math.pow(a, 2) - Math.pow(c, 2)));
        const incompleteEllipticInt1 = ts3dutils.gaussLegendreQuadrature24(phi => Math.pow(1 - k2 * Math.pow(Math.sin(phi), 2), -0.5), 0, phi);
        const incompleteEllipticInt2 = ts3dutils.gaussLegendreQuadrature24(phi => Math.pow(1 - k2 * Math.pow(Math.sin(phi), 2), 0.5), 0, phi);
        return 2 * PI$6 * Math.pow(c, 2) + 2 * PI$6 * a * b / Math.sin(phi) * (incompleteEllipticInt2 * Math.pow(Math.sin(phi), 2) + incompleteEllipticInt1 * Math.pow(Math.cos(phi), 2));
    }
    getSeamPlane() {
        return P3.forAnchorAndPlaneVectors(this.center, this.f1, this.f3);
    }
}
EllipsoidSurface.UNIT = new EllipsoidSurface(ts3dutils.V3.O, ts3dutils.V3.X, ts3dutils.V3.Y, ts3dutils.V3.Z);
EllipsoidSurface.prototype.uStep = PI$6 / 32;
EllipsoidSurface.prototype.vStep = PI$6 / 32;

const { sign: sign$4 } = Math;
/**
 * Surface normal1 is (t, z) => this.baseCurve.tangentAt(t) X this.dir
 * Choose dir appropriately to select surface orientation.
 */
class ProjectedCurveSurface extends ParametricSurface {
    constructor(baseCurve, dir, sMin = baseCurve.tMin, sMax = baseCurve.tMax, tMin = -100, tMax = 100) {
        super();
        this.baseCurve = baseCurve;
        this.dir = dir;
        this.sMin = sMin;
        this.sMax = sMax;
        this.tMin = tMin;
        this.tMax = tMax;
        ts3dutils.assertInst(Curve, baseCurve);
        ts3dutils.assertInst(ts3dutils.V3, dir);
        ts3dutils.assertNumbers(sMin, sMax, tMin, tMax);
        ts3dutils.assert(sMin < sMax);
        ts3dutils.assert(tMin < tMax);
    }
    getConstructorParameters() {
        return [this.baseCurve, this.dir, this.sMin, this.sMax, this.tMin, this.tMax];
    }
    equals(obj) {
        return this == obj ||
            Object.getPrototypeOf(this) == Object.getPrototypeOf(obj)
                && this.dir.equals(obj.dir)
                && this.baseCurve.equals(obj.baseCurve);
    }
    hashCode() {
        return [this.dir, this.baseCurve].hashCode();
    }
    containsLine(line) {
        return this.dir.isParallelTo(line.dir1) && this.containsPoint(line.anchor);
    }
    dpds() {
        return (s, t) => this.baseCurve.tangentAt(s);
    }
    dpdt() {
        return (s, t) => this.dir;
    }
    normalST(s, t) {
        return this.baseCurve.tangentAt(s).cross(this.dir).unit();
    }
    pST(s, t) {
        return this.baseCurve.at(s).plus(this.dir.times(t));
    }
    pointFoot(pWC, ss, st) {
        const basePlane = new P3(this.dir, 0);
        const projCurve = this.baseCurve.project(basePlane);
        const projPoint = basePlane.projectedPoint(pWC);
        const t = projCurve.closestTToPoint(projPoint, ss);
        const z = pWC.minus(this.baseCurve.at(t)).dot(this.dir);
        return new ts3dutils.V3(t, z, 0);
    }
    stPFunc() {
        const projPlane = new P3(this.dir.unit(), 0);
        const projBaseCurve = this.baseCurve.project(projPlane);
        return (pWC) => {
            const projPoint = projPlane.projectedPoint(pWC);
            const t = projBaseCurve.pointT(projPoint);
            const z = L3$1.pointT(this.baseCurve.at(t), this.dir, pWC);
            return new ts3dutils.V3(t, z, 0);
        };
    }
    isCurvesWithPlane(plane) {
        ts3dutils.assertInst(P3, plane);
        if (this.dir.isPerpendicularTo(plane.normal1)) {
            const ts = this.baseCurve.isTsWithPlane(plane);
            return ts.map(t => {
                const l3dir = 0 < this.baseCurve.tangentAt(t).dot(plane.normal1)
                    ? this.dir
                    : this.dir.negated();
                return new L3$1(this.baseCurve.at(t), l3dir.unit());
            });
        }
        else {
            let projCurve = this.baseCurve.transform(ts3dutils.M4.project(plane, this.dir));
            if (this.dir.dot(plane.normal1) > 0) {
                // we need to flip the ellipse so the tangent is correct
                projCurve = projCurve.reversed();
            }
            return [projCurve];
        }
    }
    isCurvesWithSurface(surface) {
        if (surface instanceof PlaneSurface$1) {
            return this.isCurvesWithPlane(surface.plane);
        }
        if (surface instanceof ProjectedCurveSurface) {
            const dir1 = surface.dir;
            if (this.dir.isParallelTo(dir1)) {
                const otherCurve = surface.baseCurve;
                const infos = this.baseCurve.isInfosWithCurve(otherCurve);
                return infos.map(info => {
                    const correctDir = this.normalP(info.p).cross(surface.normalP(info.p));
                    return new L3$1(info.p, dir1.times(sign$4(correctDir.dot(dir1))));
                });
            }
            if (surface instanceof ProjectedCurveSurface) {
                const line = new L3$1(this.baseCurve.at(0.5), this.dir);
                const startPoint = line.at(surface.isTsForLine(line)[0]);
                console.log(startPoint);
                return [new PPCurve(this, surface, startPoint)];
                // const testVector = this.dir.cross(surface.dir).unit()
                // // look for points on surface.baseCurve where tangent DOT testVector == 0
                // const abcd1 = surface.baseCurve.tangentCoefficients().map(c => c.dot(testVector))
                // const ts1 = solveCubicReal2.apply(undefined, abcd1).concat(surface.sMin, surface.sMax)
                // const abcd2 = this.baseCurve.tangentCoefficients().map(c => c.dot(testVector))
                // const ts2 = solveCubicReal2.apply(undefined, abcd2)
                // const tt1 = ts1.map(t => surface.baseCurve.at(t).dot(testVector))
                // const tt2 = ts1.map(t => surface.baseCurve.at(t).dot(testVector))
                // console.log(ts1, ts2, tt1, tt2)
                // ts1.forEach(t => drPs.push(surface.baseCurve.at(t)))
                // ts2.forEach(t => drPs.push(this.baseCurve.at(t)))
                // return
            }
        }
        if (surface instanceof SemiEllipsoidSurface) {
            return surface.isCurvesWithSurface(this);
        }
        ts3dutils.assertNever();
    }
    containsPoint(pWC) {
        const uv = this.stPFunc()(pWC);
        return this.pSTFunc()(uv.x, uv.y).like(pWC);
    }
    containsCurve(curve) {
        if (curve instanceof L3$1) {
            return this.dir.isParallelTo(curve.dir1) && this.containsPoint(curve.anchor);
        }
        if (curve instanceof PICurve$1) {
            return super.containsCurve(curve);
        }
        // project baseCurve and test curve onto a common plane and check if the curves are alike
        const projPlane = new P3(this.dir.unit(), 0);
        const projBaseCurve = this.baseCurve.project(projPlane);
        const projCurve = curve.project(projPlane);
        return projBaseCurve.isColinearTo(projCurve);
    }
    isCoplanarTo(surface) {
        return this == surface ||
            ts3dutils.hasConstructor(surface, ProjectedCurveSurface)
                && this.dir.isParallelTo(surface.dir)
                && this.containsCurve(surface.baseCurve);
    }
    like(object) {
        if (!this.isCoplanarTo(object))
            return false;
        // normals need to point in the same direction (outwards or inwards) for both
        const p00 = this.pSTFunc()(0, 0);
        const thisNormal = this.normalSTFunc()(0, 0);
        const otherNormal = object.normalP(p00);
        return 0 < thisNormal.dot(otherNormal);
    }
    loopContainsPoint(loop, p) {
        ts3dutils.assertVectors(p);
        ts3dutils.assert(isFinite(p.x), p.y, p.z);
        const line = new L3$1(p, this.dir.unit());
        const ptpf = this.stPFunc();
        const pp = ptpf(p);
        if (isNaN(pp.x)) {
            console.log(this.sce, p.sce);
            ts3dutils.assert(false);
        }
        const lineOut = this.baseCurve.tangentAt(pp.x).rejectedFrom(this.dir);
        return Surface.loopContainsPointGeneral(loop, p, line, lineOut);
    }
    edgeLoopCCW(loop) {
        if (loop.length < 56) {
            let totalAngle = 0;
            for (let i = 0; i < loop.length; i++) {
                const ipp = (i + 1) % loop.length;
                const edge = loop[i], nextEdge = loop[ipp];
                totalAngle += edge.bDir.angleRelativeNormal(nextEdge.aDir, this.normalP(edge.b));
            }
            return totalAngle > 0;
        }
        else {
            const ptpF = this.stPFunc();
            return ts3dutils.isCCW(loop.map(e => ptpF(e.a)), ts3dutils.V3.Z);
        }
    }
    transform(m4) {
        const f = m4.isMirroring() ? -1 : 1;
        return new this.constructor(this.baseCurve.transform(m4), m4.transformVector(this.dir).times(f), this.sMin, this.sMax, 1 == f ? this.tMin : -this.tMax, 1 == f ? this.tMax : -this.tMin);
    }
    isTsForLine(line) {
        ts3dutils.assertInst(L3$1, line);
        const projPlane = new P3(this.dir.unit(), 0);
        const projDir = projPlane.projectedVector(line.dir1);
        if (projDir.likeO()) {
            // line is parallel to this.dir
            return [];
        }
        const projAnchor = projPlane.projectedPoint(line.anchor);
        const projBaseCurve = this.baseCurve.project(projPlane);
        return projBaseCurve
            .isInfosWithLine(projAnchor, projDir, this.sMin, this.sMax, line.tMin, line.tMax)
            .map(info => info.tOther);
    }
    flipped() {
        return new this.constructor(this.baseCurve, this.dir.negated(), this.sMin, this.sMax, -this.tMax, -this.tMin);
    }
}
ProjectedCurveSurface.prototype.uStep = 1 / 40;
ProjectedCurveSurface.prototype.vStep = 256;

const { PI: PI$7 } = Math;
class CylinderSurface extends ProjectedCurveSurface {
    constructor(baseEllipse, dir1, zMin = -Infinity, zMax = Infinity) {
        super(baseEllipse, dir1, undefined, undefined, zMin, zMax);
        ts3dutils.assert(2 == arguments.length);
        ts3dutils.assertVectors(dir1);
        ts3dutils.assertInst(EllipseCurve, baseEllipse);
        //assert(!baseCurve.normal1.isPerpendicularTo(dir1), !baseCurve.normal1.isPerpendicularTo(dir1))
        ts3dutils.assert(dir1.hasLength(1));
        this.matrix = ts3dutils.M4.forSys(baseEllipse.f1, baseEllipse.f2, dir1, baseEllipse.center);
        this.inverseMatrix = this.matrix.inversed();
    }
    static cylinder(radius) {
        return new CylinderSurface(new EllipseCurve(ts3dutils.V3.O, new ts3dutils.V3(radius, 0, 0), new ts3dutils.V3(0, radius, 0)), ts3dutils.V3.Z);
    }
    /**
     *
     * @param anchor
     * @param dir not necessarily unit
     */
    static unitISLineTs(anchor, dir) {
        const { x: ax, y: ay } = anchor;
        const { x: dx, y: dy } = dir;
        // this cylinder: x² + y² = 1
        // line: p = anchor + t * dir
        // split line equation into 3 component equations, insert into cylinder equation
        // x = ax + t * dx
        // y = ay + t * dy
        // (ax² + 2 ax t dx + t²dx²) + (ay² + 2 ay t dy + t²dy²) = 1
        // transform to form (a t² + b t + c = 0) and solve with pqFormula
        const a = Math.pow(dx, 2) + Math.pow(dy, 2);
        const b = 2 * (ax * dx + ay * dy);
        const c = Math.pow(ax, 2) + Math.pow(ay, 2) - 1;
        return ts3dutils.pqFormula(b / a, c / a);
    }
    getConstructorParameters() {
        return [this.baseCurve, this.dir];
    }
    loopContainsPoint(loop, p) {
        ts3dutils.assertVectors(p);
        // create plane that goes through cylinder seam
        const line = new L3$1(p, this.dir);
        const seamBase = this.baseCurve.at(PI$7);
        const lineOut = this.dir.cross(this.normalP(p));
        return Surface.loopContainsPointGeneral(loop, p, line, lineOut);
    }
    isTsForLine(line) {
        ts3dutils.assertInst(L3$1, line);
        // transforming line manually has advantage that dir1 will not be renormalized,
        // meaning that calculated values t for localLine are directly transferable to line
        const localDir = this.inverseMatrix.transformVector(line.dir1);
        if (localDir.isParallelTo(ts3dutils.V3.Z)) {
            // line is parallel to this.dir
            return [];
        }
        const localAnchor = this.inverseMatrix.transformPoint(line.anchor);
        ts3dutils.assert(!CylinderSurface.unitISLineTs(localAnchor, localDir).length || !isNaN(CylinderSurface.unitISLineTs(localAnchor, localDir)[0]), 'sad ' + localDir);
        return CylinderSurface.unitISLineTs(localAnchor, localDir);
    }
    isCoplanarTo(surface) {
        return this == surface ||
            surface instanceof CylinderSurface
                && this.dir.isParallelTo(surface.dir)
                && this.containsEllipse(surface.baseCurve);
    }
    like(object) {
        if (!this.isCoplanarTo(object))
            return false;
        // normals need to point in the same direction (outwards or inwards) for both
        const thisFacesOut = 0 < this.baseCurve.normal.dot(this.dir);
        const objectFacesOut = 0 < object.baseCurve.normal.dot(object.dir);
        return thisFacesOut == objectFacesOut;
    }
    containsEllipse(ellipse) {
        const ellipseProjected = ellipse.transform(ts3dutils.M4.project(this.baseCurve.getPlane(), this.dir));
        return this.baseCurve == ellipse || this.baseCurve.isColinearTo(ellipseProjected);
    }
    containsCurve(curve) {
        if (curve instanceof EllipseCurve) {
            return this.containsEllipse(curve);
        }
        else if (curve instanceof L3$1) {
            return this.containsLine(curve);
        }
        else if (curve instanceof SemiEllipseCurve) {
            return this.containsEllipse(curve);
        }
        else {
            ts3dutils.assert(false);
        }
    }
    normalP(p) {
        const pLC = this.inverseMatrix.transformPoint(p);
        return this.normalSTFunc()(pLC.angleXY(), pLC.z);
    }
    implicitFunction() {
        return (pWC) => {
            const p = this.inverseMatrix.transformPoint(pWC);
            const radiusLC = p.lengthXY();
            const normalDir = Math.sign(this.baseCurve.normal.dot(this.dir));
            return normalDir * (1 - radiusLC);
        };
    }
    containsPoint(p) {
        return ts3dutils.eq0(this.implicitFunction()(p));
    }
    pointToParameterFunction() {
        return (pWC, hint) => {
            const pLC = this.inverseMatrix.transformPoint(pWC);
            let angle = pLC.angleXY();
            if (abs(angle) > Math.PI - ts3dutils.NLA_PRECISION) {
                ts3dutils.assert(hint == -PI$7 || hint == PI$7);
                angle = hint;
            }
            return new ts3dutils.V3(angle, pLC.z, 0);
        };
    }
    isCurvesWithSurface(surface) {
        if (surface instanceof PlaneSurface) {
            return this.isCurvesWithPlane(surface.plane);
        }
        else if (surface instanceof CylinderSurface) {
            if (surface.dir.isParallelTo(this.dir)) {
                const projEllipse = surface.baseCurve.transform(ts3dutils.M4.project(this.baseCurve.getPlane(), this.dir));
                return this.baseCurve.isInfosWithEllipse(projEllipse).map(info => new L3$1(info.p, this.dir));
            }
            else if (ts3dutils.eq0(this.getCenterLine().distanceToLine(surface.getCenterLine()))) {
                ts3dutils.assert(false);
            }
            else {
                ts3dutils.assert(false);
            }
        }
    }
    getCenterLine() {
        return new L3$1(this.baseCurve.center, this.dir);
    }
    edgeLoopCCW(loop) {
        if (loop.length < 56) {
            let totalAngle = 0;
            for (let i = 0; i < loop.length; i++) {
                const ipp = (i + 1) % loop.length;
                const edge = loop[i], nextEdge = loop[ipp];
                totalAngle += edge.bDir.angleRelativeNormal(nextEdge.aDir, this.normalP(edge.b));
            }
            return totalAngle > 0;
        }
        else {
            const ptpF = this.stPFunc();
            return ts3dutils.isCCW(loop.map(e => ptpF(e.a)), ts3dutils.V3.Z);
        }
    }
    facesOutwards() {
        return this.baseCurve.normal.dot(this.dir) > 0;
    }
    getSeamPlane() {
        return P3.forAnchorAndPlaneVectors(this.baseCurve.center, this.baseCurve.f1, this.dir);
    }
}
CylinderSurface.UNIT = new CylinderSurface(EllipseCurve.XY, ts3dutils.V3.Z);
CylinderSurface.prototype.uStep = ts3dutils.TAU / 128;
CylinderSurface.prototype.vStep = 256;

const { PI: PI$8, cos: cos$5, sin: sin$5, min: min$6, max: max$5, sign: sign$5, tan: tan$3, ceil: ceil$8, floor: floor$7, abs: abs$8, sqrt: sqrt$3, pow: pow$3, atan2: atan2$3, round: round$3 } = Math;
/**
 * Rotation surface with r = f(z)
 */
class RotationREqFOfZ extends ParametricSurface {
    constructor(matrix, rt, // r(z)
    tMin, tMax, normalDir, drdz = z => (rt(z + EPS) - rt(z)) / EPS) {
        super();
        this.matrix = matrix;
        this.rt = rt;
        this.tMin = tMin;
        this.tMax = tMax;
        this.normalDir = normalDir;
        this.drdz = drdz;
        ts3dutils.assertInst(ts3dutils.M4, matrix);
        ts3dutils.assert(matrix.isNoProj());
        ts3dutils.assert(1 == normalDir || -1 == normalDir);
        this.matrixInverse = matrix.inversed();
    }
    getConstructorParameters() {
        return [this.matrix, this.rt, this.tMin, this.tMax, this.normalDir, this.drdz];
    }
    flipped() {
        return new RotationREqFOfZ(this.matrix, this.rt, this.tMin, this.tMax, -this.normalDir, this.drdz);
    }
    transform(m4) {
        return new RotationREqFOfZ(m4.times(this.matrix), this.rt, this.tMin, this.tMax, this.normalDir, this.drdz);
    }
    containsPoint(p) {
        return ts3dutils.eq0(this.implicitFunction()(p));
    }
    pSTFunc() {
        return (d, z) => {
            const radius = this.rt(z);
            return this.matrix.transformPoint(ts3dutils.V3.polar(radius, d, z));
        };
    }
    dpds() {
        return (s, t) => {
            const radius = this.rt(t);
            return this.matrix.transformVector(new ts3dutils.V3(radius * -sin$5(s), radius * cos$5(s), 0));
        };
    }
    /**
     * new V3(f(z) * cos d, f(z) * sin d, z)
     */
    dpdt() {
        return (s, t) => {
            const drdt = this.drdz(t);
            return this.matrix.transformVector(new ts3dutils.V3(drdt * cos$5(s), drdt * sin$5(s), 1));
        };
    }
    normalSTFunc() {
        /**
         * (radius * -sin(s), radius * cos(s), 0) X (drds * cos(s), drds * sin(s), 1)
         * =(radius * cos(s)*1,
         * -radius * -sin(s)*1,
         * radius * -sin(s)* drds * sin(s)- radius * cos(s)*drds * cos(s))
         * div by radius
         * => (cos s, sin s, -drds * (sin² + cos²))
         */
        const matrix = this.matrix.inversed().transposed();
        return (d, z) => {
            const drdz = this.drdz(z);
            return matrix.transformVector(ts3dutils.V3.polar(1, d, -drdz)).toLength(this.normalDir);
        };
    }
    implicitFunction() {
        return (pWC) => {
            const pLC = this.matrixInverse.transformPoint(pWC);
            const radiusLC = pLC.lengthXY();
            return this.rt(pLC.z) - radiusLC;
        };
    }
    stPFunc() {
        return (pWC) => {
            const pLC = this.matrixInverse.transformPoint(pWC);
            return new ts3dutils.V3(atan2$3(pLC.y, pLC.x), pLC.z, 0);
        };
    }
}
Object.assign(RotationREqFOfZ.prototype, ImplicitSurface.prototype);
RotationREqFOfZ.prototype.sMin = 0;
RotationREqFOfZ.prototype.sMax = PI$8;

const { PI: PI$9, cos: cos$6, sin: sin$6, min: min$7, max: max$6, tan: tan$4, sign: sign$6, ceil: ceil$9, floor: floor$8, abs: abs$9, sqrt: sqrt$4, pow: pow$4, atan2: atan2$4, round: round$4 } = Math;
class SemiCylinderSurface extends ProjectedCurveSurface {
    constructor(baseCurve, dir1, sMin, sMax, zMin = -Infinity, zMax = Infinity) {
        super(baseCurve, dir1, sMin, sMax, zMin, zMax);
        ts3dutils.assertInst(SemiEllipseCurve, baseCurve);
        //assert(!baseCurve.normal1.isPerpendicularTo(dir1), !baseCurve.normal1.isPerpendicularTo(dir1))
        this.matrix = ts3dutils.M4.forSys(baseCurve.f1, baseCurve.f2, dir1, baseCurve.center);
        this.inverseMatrix = this.matrix.inversed();
        this.normalDir = sign$6(this.baseCurve.normal.dot(this.dir));
        this.normalMatrix = this.matrix.as3x3().inversed().transposed().scale(this.normalDir);
    }
    static semicylinder(radius) {
        return new SemiCylinderSurface(new SemiEllipseCurve(ts3dutils.V3.O, new ts3dutils.V3(radius, 0, 0), new ts3dutils.V3(0, radius, 0)), ts3dutils.V3.Z, undefined, undefined);
    }
    /**
     *
     * @param anchorLC
     * @param dirLC not necessarily unit
     */
    static unitISLineTs(anchorLC, dirLC) {
        const { x: ax, y: ay } = anchorLC;
        const { x: dx, y: dy } = dirLC;
        // this cylinder: x² + y² = 1
        // line: p = anchorLC + t * dirLC
        // split line equation into 3 component equations, insert into cylinder equation
        // x = ax + t * dx
        // y = ay + t * dy
        // (ax² + 2 ax t dx + t²dx²) + (ay² + 2 ay t dy + t²dy²) = 1
        // transform to form (a t² + b t + c = 0) and solve with pqFormula
        const a = Math.pow(dx, 2) + Math.pow(dy, 2);
        const b = 2 * (ax * dx + ay * dy);
        const c = Math.pow(ax, 2) + Math.pow(ay, 2) - 1;
        return ts3dutils.pqFormula(b / a, c / a).filter(t => SemiEllipseCurve.XYLCValid(new ts3dutils.V3(ax + dx * t, ay + dy * t, 0)));
    }
    getConstructorParameters() {
        return [this.baseCurve, this.dir, this.sMin, this.sMax, this.tMin, this.tMax];
    }
    normalP(p) {
        return this.normalMatrix.transformVector(this.inverseMatrix.transformPoint(p).xy()).unit();
    }
    loopContainsPoint(loop, p) {
        ts3dutils.assertVectors(p);
        if (!this.containsPoint(p))
            return OUTSIDE;
        // create plane that goes through cylinder seam
        const line = new L3$1(p, this.dir.unit());
        const seamBase = this.baseCurve.at(PI$9);
        const lineOut = this.dir.cross(this.normalP(p));
        return Surface.loopContainsPointGeneral(loop, p, line, lineOut);
    }
    isTsForLine(line) {
        ts3dutils.assertInst(L3$1, line);
        // transforming line manually has advantage that dir1 will not be renormalized,
        // meaning that calculated values t for localLine are directly transferable to line
        const dirLC = this.inverseMatrix.transformVector(line.dir1);
        if (dirLC.isParallelTo(ts3dutils.V3.Z)) {
            // line is parallel to this.dir
            return [];
        }
        const anchorLC = this.inverseMatrix.transformPoint(line.anchor);
        ts3dutils.assert(!SemiCylinderSurface.unitISLineTs(anchorLC, dirLC).length || !isNaN(SemiCylinderSurface.unitISLineTs(anchorLC, dirLC)[0]), 'sad ' + dirLC);
        return SemiCylinderSurface.unitISLineTs(anchorLC, dirLC);
    }
    isCoplanarTo(surface) {
        return this == surface ||
            ts3dutils.hasConstructor(surface, SemiCylinderSurface)
                && this.dir.isParallelTo(surface.dir)
                && this.containsSemiEllipse(surface.baseCurve, false);
    }
    like(surface) {
        if (!this.isCoplanarTo(surface))
            return false;
        // normals need to point in the same direction (outwards or inwards) for both
        const thisFacesOut = 0 < this.baseCurve.normal.dot(this.dir);
        const objectFacesOut = 0 < surface.baseCurve.normal.dot(surface.dir);
        return thisFacesOut == objectFacesOut;
    }
    containsSemiEllipse(ellipse, checkAABB = true) {
        const projEllipse = ellipse.transform(ts3dutils.M4.project(this.baseCurve.getPlane(), this.dir));
        return this.baseCurve == ellipse || this.baseCurve.isColinearTo(projEllipse) &&
            (!checkAABB || ts3dutils.le(0, ellipse.transform(this.inverseMatrix).getAABB().min.y));
    }
    containsCurve(curve) {
        if (curve instanceof L3$1) {
            return this.containsLine(curve);
        }
        else if (curve instanceof SemiEllipseCurve) {
            return this.containsSemiEllipse(curve);
        }
        else if (curve instanceof BezierCurve) {
            return false;
        }
        else {
            return super.containsCurve(curve);
        }
    }
    implicitFunction() {
        return (pWC) => {
            const pLC = this.inverseMatrix.transformPoint(pWC);
            const radiusLC = pLC.lengthXY();
            const normalDir = Math.sign(this.baseCurve.normal.dot(this.dir));
            return normalDir * (1 - radiusLC);
        };
    }
    containsPoint(pWC) {
        const pLC = this.inverseMatrix.transformPoint(pWC);
        return SemiEllipseCurve.XYLCValid(pLC);
    }
    stP(pWC) {
        ts3dutils.assert(arguments.length == 1);
        const pLC = this.inverseMatrix.transformPoint(pWC);
        const u = SemiEllipseCurve.XYLCPointT(pLC);
        return new ts3dutils.V3(u, pLC.z, 0);
    }
    isCurvesWithSurface(surface2) {
        if (surface2 instanceof PlaneSurface$1) {
            return this.isCurvesWithPlane(surface2.plane);
        }
        else if (surface2 instanceof SemiCylinderSurface) {
            if (surface2.dir.isParallelTo(this.dir)) {
                const projEllipse = surface2.baseCurve.transform(ts3dutils.M4.project(this.baseCurve.getPlane(), this.dir));
                return this.baseCurve.isInfosWithEllipse(projEllipse).map(info => {
                    const lineDir = sign$6(this.normalP(info.p).cross(surface2.normalP(info.p)).dot(this.dir)) || 1;
                    return new L3$1(info.p, this.dir.times(lineDir));
                });
            }
            else if (ts3dutils.eq0(this.getCenterLine().distanceToLine(surface2.getCenterLine()))) {
                ts3dutils.assert(false);
            }
            else {
                ts3dutils.assert(false);
            }
        }
    }
    getCenterLine() {
        return new L3$1(this.baseCurve.center, this.dir);
    }
    facesOutwards() {
        return this.baseCurve.normal.dot(this.dir) > 0;
    }
    getSeamPlane() {
        let normal = this.baseCurve.f1.cross(this.dir);
        normal = normal.times(-sign$6(normal.dot(this.baseCurve.f2)));
        return P3.normalOnAnchor(normal, this.baseCurve.center);
    }
    clipCurves(curves) {
        return curves.flatMap(curve => curve.clipPlane(this.getSeamPlane()));
    }
}
SemiCylinderSurface.UNIT = new SemiCylinderSurface(SemiEllipseCurve.UNIT, ts3dutils.V3.Z, undefined, undefined, 0, 1);
SemiCylinderSurface.prototype.uStep = ts3dutils.TAU / 32;
SemiCylinderSurface.prototype.vStep = 256;

const { PI: PI$10, min: min$8, max: max$7, sign: sign$7, abs: abs$10, sqrt: sqrt$5 } = Math;
class SemiEllipsoidSurface extends EllipsoidSurface {
    constructor(center, f1, f2, f3) {
        super(center, f1, f2, f3);
        this.center = center;
        this.f1 = f1;
        this.f2 = f2;
        this.f3 = f3;
        ts3dutils.assertVectors(center, f1, f2, f3);
        this.matrix = ts3dutils.M4.forSys(f1, f2, f3, center);
        this.inverseMatrix = this.matrix.inversed();
        this.normalDir = sign$7(this.f1.cross(this.f2).dot(this.f3));
        this.pLCNormalWCMatrix = this.matrix.as3x3().inversed().transposed().scale(this.normalDir);
        this.pWCNormalWCMatrix = this.pLCNormalWCMatrix.times(this.inverseMatrix);
    }
    static unitArea(contour) {
        const totalArea = contour.map(edge => {
            if (edge.curve instanceof PICurve$1) {
                const points = edge.curve.calcSegmentPoints(edge.aT, edge.bT, edge.a, edge.b, edge.aT > edge.bT, true);
                let sum = 0;
                for (let i = 0; i < points.length - 1; i++) {
                    const p = points[i], ppp = points[i + 1];
                    sum += (abs$10(p.angleXY()) + abs$10(ppp.angleXY())) / 2 * (ppp.z - p.z);
                }
                return sum;
            }
            else if (edge.curve instanceof SemiEllipseCurve) {
                const f = (t) => {
                    const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t);
                    const angleXY = abs$10(at.angleXY());
                    //const arcLength = angleXY * Math.sqrt(1 - at.z ** 2) ( == at.lengthXY())
                    //const scaling = tangent.z / at.lengthXY()
                    return angleXY * tangent.z;
                };
                const val = ts3dutils.glqInSteps(f, edge.aT, edge.bT, 1);
                return val;
            }
            else {
                ts3dutils.assertNever();
            }
        }).sum();
        return totalArea;
    }
    /**
     * unit sphere: x² + y² + z² = 1
     * line: p = anchor + t * dir |^2
     * p² = (anchor + t * dir)^2
     * 1 == (anchor + t * dir)^2
     * 1 == anchor DOT anchor + 2 * anchor * t * dir + t² * dir DOT dir
     */
    static unitISTsWithLine(anchor, dir) {
        // for 0 = a t² + b t + c
        const a = dir.dot(dir);
        const b = 2 * anchor.dot(dir);
        const c = anchor.dot(anchor) - 1;
        return ts3dutils.pqFormula(b / a, c / a).filter(t => ts3dutils.le(0, anchor.y + t * dir.y));
    }
    /**
     * unit sphere: x² + y² + z² = 1
     * plane: normal1 DOT p = w
     */
    static unitISCurvesWithPlane(plane) {
        const distPlaneCenter = Math.abs(plane.w);
        if (ts3dutils.lt(distPlaneCenter, 1)) {
            // result is a circle
            // radius of circle: imagine right angled triangle (origin -> center of intersection circle -> point on
            // intersection circle) pythagoras: 1² == distPlaneCenter² + isCircleRadius² => isCircleRadius == sqrt(1 -
            // distPlaneCenter²)
            const isCircleRadius = Math.sqrt(1 - Math.pow(distPlaneCenter, 2));
            const anchorY = plane.normal1.y * plane.w;
            const d = abs$10(distPlaneCenter * isCircleRadius);
            if (ts3dutils.le(anchorY, -d) && !ts3dutils.eq0(distPlaneCenter)) {
                return [];
            }
            else if (ts3dutils.le(anchorY, 0) && !plane.normal1.isParallelTo(ts3dutils.V3.Y)) {
                let f1 = plane.normal1.isParallelTo(ts3dutils.V3.Y) ? ts3dutils.V3.Z : plane.normal1.cross(ts3dutils.V3.Y).toLength(isCircleRadius);
                const f2 = f1.cross(plane.normal1);
                const minEta = -anchorY / f2.y, minT = max$7(0, Math.asin(minEta));
                return [new SemiEllipseCurve(plane.anchor, f1, f2, minT, PI$10 - minT)];
            }
            else {
                const f2 = (plane.normal1.isParallelTo(ts3dutils.V3.Y)
                    ? ts3dutils.V3.X
                    : plane.normal1.cross(ts3dutils.V3.Y)).toLength(isCircleRadius);
                const f1 = f2.cross(plane.normal1);
                const minXi = ts3dutils.eq0(f1.y) ? -1 : -anchorY / f1.y, maxT = Math.acos(max$7(-1, minXi - ts3dutils.NLA_PRECISION));
                return [new SemiEllipseCurve(plane.anchor, f1.negated(), f2, PI$10 - maxT, PI$10),
                    new SemiEllipseCurve(plane.anchor, f1, f2.negated(), 0, maxT)];
            }
        }
        else {
            return [];
        }
    }
    static unitISCurvesWithEllipsoidSurface(surface) {
        if (surface.isSphere()) {
            const surfaceRadius = surface.f1.length();
            const surfaceCenterDist = surface.center.length();
            if (ts3dutils.le(1, surfaceCenterDist - surfaceRadius) || ts3dutils.le(surfaceCenterDist + surfaceRadius, 1) || ts3dutils.le(surfaceCenterDist - surfaceRadius, -1)) {
                return [];
            }
            else {
                // origin, surface.center and points on the intersection curves form a triangle.
                // the height on the segment origin - surface.center is the radius of the is curves
                // the distance from the origin to the lot point is the distance to the intersection plane
                function heron(a, b, c) {
                    const p = (a + b + c) / 2;
                    return sqrt$5(p * (p - a) * (p - b) * (p - c));
                }
                const triangleArea = heron(1, surfaceRadius, surfaceCenterDist);
                const radius = triangleArea * 2 / surfaceCenterDist;
                const isCurvesCenterDist = sign$7(1 + Math.pow(surfaceCenterDist, 2) - Math.pow(surfaceRadius, 2)) * sqrt$5(1 - Math.pow(radius, 2));
                const plane = new P3(surface.center.unit(), isCurvesCenterDist);
                return SemiEllipsoidSurface.unitISCurvesWithPlane(plane.flipped());
            }
        }
        ts3dutils.assertNever();
    }
    static unitISCurvesWithSemiCylinderSurface(surface) {
        if (new L3$1(surface.baseCurve.center, surface.dir).containsPoint(ts3dutils.V3.O)) {
            const projEllipse = surface.baseCurve.transform(ts3dutils.M4.project(new P3(surface.dir, 0)));
            const f1Length = projEllipse.f1.length(), f2Length = projEllipse.f2.length();
            if (ts3dutils.lt(1, min$8(f1Length, f2Length)))
                return [];
            if (projEllipse.isCircular()) {
                const distISCurveCenter = Math.sqrt(1 - Math.pow(min$8(1, f1Length), 2));
                const isCurveCenter = (surface.dir.y < 0 ? surface.dir.negated() : surface.dir).times(distISCurveCenter);
                // isCurve.at(t).y = isCurveCenter.y + projEllipse.f1.y * cos(t) + projEllipse.f2.y * sin(t) = 0
                return [new SemiEllipseCurve(isCurveCenter, projEllipse.f1, projEllipse.f2)];
            }
        }
        ts3dutils.assert(false);
    }
    static sphere(radius, center = ts3dutils.V3.O) {
        ts3dutils.assertNumbers(radius);
        return new SemiEllipsoidSurface(center, new ts3dutils.V3(radius, 0, 0), new ts3dutils.V3(0, radius, 0), new ts3dutils.V3(0, 0, radius));
    }
    /**
     * x²/a² + y²/b² + z²/c² = 1
     */
    static forABC(a, b, c, center = ts3dutils.V3.O) {
        return new SemiEllipsoidSurface(center, new ts3dutils.V3(a, 0, 0), new ts3dutils.V3(0, b, 0), new ts3dutils.V3(0, 0, c));
    }
    static calculateAreaSpheroid(a, b, c, edges) {
        ts3dutils.assertf(() => a.isPerpendicularTo(b));
        ts3dutils.assertf(() => b.isPerpendicularTo(c));
        ts3dutils.assertf(() => c.isPerpendicularTo(a));
        // handling discontinuities:
        // option 1: check for intersections with baseline, if there are any integrate parts separetely
        // "rotate" the edge so that there are no overlaps
        const matrix = ts3dutils.M4.forSys(a, b, c), inverseMatrix = matrix.inversed();
        const circleRadius = a.length();
        const c1 = c.unit();
        const totalArea = edges.map(edge => {
            if (edge.curve instanceof SemiEllipseCurve) {
                const f = (t) => {
                    const at = edge.curve.at(t), tangent = edge.tangentAt(t);
                    const localAt = inverseMatrix.transformPoint(at);
                    const angleXY = localAt.angleXY();
                    const arcLength = angleXY * circleRadius * Math.sqrt(1 + Math.pow(localAt.z, 2));
                    const scaling = Math.sqrt(1 + Math.pow(c1.dot(tangent), 2));
                    return arcLength * scaling;
                };
                const val = ts3dutils.glqInSteps(f, edge.aT, edge.bT, 1);
                return val;
            }
            else {
                ts3dutils.assertNever();
            }
        }).sum();
        return totalArea;
    }
    equals(obj) {
        return this == obj ||
            Object.getPrototypeOf(obj) == this.constructor.prototype
                && this.matrix.equals(obj.matrix);
    }
    edgeLoopCCW(loop) {
        return SemiEllipsoidSurface.unitArea(loop.map(edge => edge.transform(this.inverseMatrix))) > 0;
        //let totalAngle = 0
        //for (let i = 0; i < contour.length; i++) {
        //    const ipp = (i + 1) % contour.length
        //    const edge = contour[i], nextEdge = contour[ipp]
        //    totalAngle += edge.bDir.angleRelativeNormal(nextEdge.aDir, this.normalP(edge.b))
        //}
        //return le(0, totalAngle)
    }
    like(object) {
        if (!this.isCoplanarTo(object))
            return false;
        // normals need to point in the same direction (outwards or inwards) for both
        return this.matrix.determinant3() * object.matrix.determinant3() > 0;
    }
    rootPoints() {
    }
    toMesh() {
        return ParametricSurface.prototype.toMesh.call(this);
    }
    getConstructorParameters() {
        return [this.center, this.f1, this.f2, this.f3];
    }
    clipCurves(curves) {
        return curves.flatMap(curve => curve.clipPlane(this.getSeamPlane()));
    }
    isCurvesWithPCS(surface) {
        let curves2 = ParametricSurface.isCurvesParametricImplicitSurface(surface, this, 0.1, 0.1 / surface.dir.length(), 0.05);
        curves2 = this.clipCurves(curves2);
        curves2 = surface.clipCurves(curves2);
        return curves2;
        const surfaceLC = surface.transform(this.inverseMatrix);
        //const lcMinZ0RelO =
        const baseCurveLC = surfaceLC.baseCurve.project(new P3(surfaceLC.dir, 0));
        const ists = baseCurveLC.isTsWithSurface(EllipsoidSurface.UNIT);
        const insideIntervals = ts3dutils.getIntervals(ists, baseCurveLC.tMin, baseCurveLC.tMax)
            .filter(([a, b]) => baseCurveLC.at((a + b) / 2).length() < 1);
        const projectedCurves = [0, 1].map(id => {
            return (t) => {
                const atSqr = ts3dutils.snap(baseCurveLC.at(t).squared(), 1);
                const lineISTs = sqrt$5(1 - atSqr);
                //assert(!isNaN(lineISTs))
                return ts3dutils.eq0(lineISTs)
                    ? baseCurveLC.at(t)
                    : baseCurveLC.at(t).plus(surfaceLC.dir.times(sign$7(id - 0.5) * lineISTs));
            };
        });
        const dProjectedCurves = [0, 1].map(id => {
            return (t) => {
                // d/dt sqrt(1 - baseCurveLC.at(t).squared())
                // = -1/2 * 1/sqrt(1 - baseCurveLC.at(t).squared()) * -2*baseCurveLC.at(t) * baseCurveLC.tangentAt(t)
                const atSqr = ts3dutils.snap(baseCurveLC.at(t).squared(), 1);
                const lineISTs = baseCurveLC.at(t).times(-1 / sqrt$5(1 - atSqr)).dot(baseCurveLC.tangentAt(t));
                //assert(!isNaN(lineISTs))
                return baseCurveLC.tangentAt(t).plus(surfaceLC.dir.times(sign$7(id - 0.5) * lineISTs));
            };
        });
        //const f2 = t => sqrt(1 - baseCurveLC.at(t).squared())
        //const df2 = t => baseCurveLC.at(t).times(-1 / sqrt(1 -
        // baseCurveLC.at(t).squared())).dot(baseCurveLC.tangentAt(t)) checkDerivate(f2, df2, 0.31, 0.60)
        const curves = [];
        for (const [aT, bT] of insideIntervals) {
            //const aLine = new L3(baseCurveLC.at(aT), surfaceLC.dir1)
            //const a = EllipsoidSurface.UNIT.isTsForLine(aLine).map(t => aLine.at(t))
            //const bLine = new L3(baseCurveLC.at(bT), surfaceLC.dir1)
            //const b = EllipsoidSurface.UNIT.isTsForLine(bLine).map(t => bLine.at(t))
            for (const i of [0, 1]) {
                const f = (t) => projectedCurves[i](t).y;
                const df = (t) => dProjectedCurves[i](t).y;
                ts3dutils.checkDerivate(f, df, aT + 0.1, bT - 0.1);
                const tsAtY0 = ts3dutils.getRoots(f, aT + ts3dutils.NLA_PRECISION, bT - ts3dutils.NLA_PRECISION, 1 / (1 << 11), df);
                const ii2 = ts3dutils.getIntervals(tsAtY0, aT, bT).filter(([a, b]) => f((a + b) / 2) > 0);
                for (const [aT2, bT2] of ii2) {
                    let aP = projectedCurves[i](aT2), bP = projectedCurves[i](bT2);
                    0 === i && ([aP, bP] = [bP, aP]);
                    ts3dutils.assert(EllipsoidSurface.UNIT.containsPoint(aP));
                    ts3dutils.assert(EllipsoidSurface.UNIT.containsPoint(bP));
                    curves.push(PICurve$1.forStartEnd(surface, this, this.matrix.transformPoint(bP), this.matrix.transformPoint(aP), undefined, 1));
                }
            }
        }
        return surface.clipCurves(curves);
    }
    isCurvesWithSurface(surface) {
        if (surface instanceof PlaneSurface$1) {
            return this.isCurvesWithPlane(surface.plane);
        }
        else if (surface instanceof SemiCylinderSurface) {
            return this.isCurvesWithSemiCylinderSurface(surface);
        }
        else if (surface instanceof SemiEllipsoidSurface) {
            const surfaceLC = surface.transform(this.inverseMatrix);
            const curves = SemiEllipsoidSurface.unitISCurvesWithEllipsoidSurface(surfaceLC)
                .map(c => c.transform(this.matrix));
            return surface.clipCurves(curves);
        }
        else if (surface instanceof ProjectedCurveSurface) {
            return this.isCurvesWithPCS(surface);
        }
        else if (surface instanceof ParametricSurface) {
            let curves2 = ParametricSurface.isCurvesParametricImplicitSurface(surface, this, 0.1, 0.1, 0.05);
            curves2 = this.clipCurves(curves2);
            curves2 = surface.clipCurves(curves2);
            return curves2;
        }
        else {
            ts3dutils.assert(false);
        }
    }
    isCurvesWithPlane(plane) {
        const planeLC = plane.transform(this.inverseMatrix);
        return SemiEllipsoidSurface.unitISCurvesWithPlane(planeLC).map(c => c.transform(this.matrix));
    }
    isCurvesWithSemiCylinderSurface(surface) {
        if (L3$1.containsPoint(surface.baseCurve.center, surface.dir, this.center)) {
            ts3dutils.assert(this.isSphere());
            const ellipseProjected = surface.baseCurve.transform(ts3dutils.M4.project(surface.baseCurve.getPlane(), surface.dir));
            if (ellipseProjected.isCircular()) {
                const thisRadius = this.f1.length();
                const surfaceRadius = ellipseProjected.f1.length();
                // sphereRadius² = distanceISFromCenter² + isRadius²
                if (ts3dutils.eq(thisRadius, surfaceRadius)) {
                    // return
                }
                else if (surfaceRadius < thisRadius) {
                }
                ts3dutils.assert(false);
            }
        }
        return this.isCurvesWithPCS(surface);
    }
    isTsForLine(line) {
        ts3dutils.assertInst(L3$1, line);
        // transforming line manually has advantage that dir1 will not be renormalized,
        // meaning that calculated values t for localLine are directly transferable to line
        const anchorLC = this.inverseMatrix.transformPoint(line.anchor);
        const dirLC = this.inverseMatrix.transformVector(line.dir1);
        return SemiEllipsoidSurface.unitISTsWithLine(anchorLC, dirLC);
    }
    isCoplanarTo(surface) {
        if (this === surface)
            return true;
        if (!ts3dutils.hasConstructor(surface, SemiEllipsoidSurface))
            return false;
        if (!this.center.like(surface.center))
            return false;
        if (this.isSphere())
            return surface.isSphere() && ts3dutils.eq(this.f1.length(), this.f2.length());
        const otherMatrixLC = this.inverseMatrix.times(surface.matrix);
        // Ellipsoid with matrix otherMatrixLC is unit sphere iff otherMatrixLC is orthogonal
        return otherMatrixLC.is3x3() && otherMatrixLC.isOrthogonal();
    }
    containsEllipse(ellipse) {
        const ellipseLC = ellipse.transform(this.inverseMatrix);
        const distEllipseLCCenter = ellipseLC.center.length();
        const correctRadius = Math.sqrt(1 - Math.pow(distEllipseLCCenter, 2));
        return ts3dutils.lt(distEllipseLCCenter, 1)
            && ellipseLC.isCircular()
            && ellipseLC.f1.hasLength(correctRadius);
        //&& le(0, ellipseLC.getAABB().min.y)
    }
    containsCurve(curve) {
        if (curve instanceof SemiEllipseCurve) {
            return this.containsEllipse(curve);
        }
        else {
            return super.containsCurve(curve);
        }
    }
    transform(m4) {
        return new SemiEllipsoidSurface(m4.transformPoint(this.center), m4.transformVector(this.f1), m4.transformVector(this.f2), m4.transformVector(this.f3).times(m4.isMirroring() ? -1 : 1));
    }
    isInsideOut() {
        return this.f1.cross(this.f2).dot(this.f3) < 0;
    }
    //implicitFunction() {
    //    return (pWC) => {
    //        const pLC = this.inverseMatrix.transformPoint(pWC)
    //        return (pLC.y > 0
    //            ? pLC.length() - 1
    //            : (-pLC.y + Math.hypot(pLC.x, pLC.z) - 1)) * this.normalDir
    //    }
    //}
    //didp(pWC) {
    //    const pLC = this.inverseMatrix.transformPoint(pWC)
    //    const didpLC = (pLC.y > 0
    //                ? pLC.unit()
    //                : V(pLC.x / Math.hypot(pLC.x, pLC.z), -1, pLC.z / Math.hypot(pLC.x, pLC.z))).times(this.normalDir)
    //    return this.inverseMatrix.transformVector(didpLC)
    //}
    flipped() {
        return new SemiEllipsoidSurface(this.center, this.f1, this.f2, this.f3.negated());
    }
    normalSTFunc() {
        // ugh
        // paramtric ellipsoid point q(a, b)
        // normal1 == (dq(a, b) / da) X (dq(a, b) / db) (Cross product of partial derivatives
        // normal1 == cos b * (f2 X f3 * cos b * cos a + f3 X f1 * cos b * sin a + f1 X f2 * sin b)
        return (a, b) => {
            const { f1, f2, f3 } = this;
            const normal = f2.cross(f3).times(Math.cos(b) * Math.cos(a))
                .plus(f3.cross(f1).times(Math.cos(b) * Math.sin(a)))
                .plus(f1.cross(f2).times(Math.sin(b)))
                .unit();
            return normal;
        };
    }
    normalP(p) {
        return this.pLCNormalWCMatrix.transformVector(this.inverseMatrix.transformPoint(p)).unit();
    }
    normalST(s, t) {
        return this.pLCNormalWCMatrix.transformVector(ts3dutils.V3.sphere(s, t)).unit();
    }
    stPFunc() {
        return (pWC) => {
            const pLC = this.inverseMatrix.transformPoint(pWC);
            const alpha = abs$10(pLC.angleXY());
            const beta = Math.asin(ts3dutils.clamp(pLC.z, -1, 1));
            ts3dutils.assert(isFinite(alpha));
            ts3dutils.assert(isFinite(beta));
            return new ts3dutils.V3(alpha, beta, 0);
        };
    }
    pSTFunc() {
        // this(a, b) = f1 cos a cos b + f2 sin a cos b + f2 sin b
        return (alpha, beta) => {
            return this.matrix.transformPoint(ts3dutils.V3.sphere(alpha, beta));
        };
    }
    isSphere() {
        return ts3dutils.eq(this.f1.length(), this.f2.length())
            && ts3dutils.eq(this.f2.length(), this.f3.length())
            && ts3dutils.eq(this.f3.length(), this.f1.length())
            && this.f1.isPerpendicularTo(this.f2)
            && this.f2.isPerpendicularTo(this.f3)
            && this.f3.isPerpendicularTo(this.f1);
    }
    isVerticalSpheroid() {
        return ts3dutils.eq(this.f1.length(), this.f2.length())
            && this.f1.isPerpendicularTo(this.f2)
            && this.f2.isPerpendicularTo(this.f3)
            && this.f3.isPerpendicularTo(this.f1);
    }
    mainAxes() {
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
        const { f1, f2, f3 } = this;
        if (ts3dutils.eq0(f1.dot(f2)) && ts3dutils.eq0(f2.dot(f3)) && ts3dutils.eq0(f3.dot(f1))) {
            return this;
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
        const { U, SIGMA } = this.matrix.svd3();
        ts3dutils.assert(SIGMA.isDiagonal());
        ts3dutils.assert(U.isOrthogonal());
        const U_SIGMA = U.times(SIGMA);
        // column vectors of U_SIGMA
        const [mainF1, mainF2, mainF3] = ts3dutils.arrayFromFunction(3, i => new ts3dutils.V3(U_SIGMA.m[i], U_SIGMA.m[i + 4], U_SIGMA.m[i + 8]));
        return new SemiEllipsoidSurface(this.center, mainF1, mainF2, mainF3);
    }
    containsPoint(p) {
        return ts3dutils.eq0(this.implicitFunction()(p));
    }
    boundsFunction() {
        return (a, b) => ts3dutils.between(a, 0, PI$10) && ts3dutils.between(b, -PI$10, PI$10);
    }
    volume() {
        return 4 / 3 * Math.PI * this.f1.dot(this.f2.cross(this.f3));
    }
    loopContainsPoint(loop, p) {
        if (!this.containsPoint(p))
            return exports.PointVsFace.OUTSIDE;
        ts3dutils.assertVectors(p);
        const pLCXY = this.inverseMatrix.transformPoint(p).withElement('z', 0);
        const testLine = new SemiEllipseCurve(this.center, this.f3, pLCXY.likeO() ? this.f2 : this.matrix.transformVector(pLCXY.unit()));
        const pT = testLine.pointT(p);
        if (P3.normalOnAnchor(this.f2.unit(), this.center).containsPoint(p)) {
            let edgeT;
            return loop.some(edge => edge.curve.containsPoint(p) && ts3dutils.le(edge.minT, edgeT = edge.curve.pointT(p)) && ts3dutils.le(edgeT, edge.maxT))
                ? exports.PointVsFace.ON_EDGE
                : exports.PointVsFace.OUTSIDE;
        }
        const lineOut = testLine.normal;
        const testPlane = P3.normalOnAnchor(testLine.normal, p);
        const colinearEdges = loop.map((edge) => testLine.isColinearTo(edge.curve));
        let inside = false;
        function logIS(isP) {
            const isT = testLine.pointT(isP);
            if (ts3dutils.eq(pT, isT)) {
                return true;
            }
            else if (pT < isT && ts3dutils.le(isT, PI$10)) {
                inside = !inside;
            }
        }
        for (let edgeIndex = 0; edgeIndex < loop.length; edgeIndex++) {
            const edge = loop[edgeIndex];
            const nextEdgeIndex = (edgeIndex + 1) % loop.length, nextEdge = loop[nextEdgeIndex];
            //console.log(edge.toSource()) {p:V(2, -2.102, 0),
            if (colinearEdges[edgeIndex]) {
                let edgeT;
                if (edge.curve.containsPoint(p) && ts3dutils.le(edge.minT, edgeT = edge.curve.pointT(p)) && ts3dutils.le(edgeT, edge.maxT)) {
                    return exports.PointVsFace.ON_EDGE;
                }
                // edge colinear to intersection
                const nextInside = colinearEdges[nextEdgeIndex] || dotCurve(lineOut, nextEdge.aDir, nextEdge.aDDT) < 0;
                if (!nextInside && testLine.containsPoint(edge.b)) {
                    if (logIS(edge.b))
                        return exports.PointVsFace.ON_EDGE;
                }
            }
            else {
                for (const edgeT of edge.edgeISTsWithPlane(testPlane)) {
                    if (edgeT == edge.bT) {
                        if (!testLine.containsPoint(edge.b))
                            continue;
                        // endpoint lies on intersection testLine
                        const edgeInside = dotCurve2(edge.curve, edge.bT, lineOut, -sign$7(edge.deltaT())) < 0; // TODO:
                        // bDDT
                        // negated?
                        const nextInside = colinearEdges[nextEdgeIndex] || dotCurve(lineOut, nextEdge.aDir, nextEdge.aDDT) < 0;
                        if (edgeInside != nextInside) {
                            if (logIS(edge.b))
                                return exports.PointVsFace.ON_EDGE;
                        }
                    }
                    else if (edgeT != edge.aT) {
                        const p = edge.curve.at(edgeT);
                        if (!testLine.containsPoint(p))
                            continue;
                        // edge crosses testLine, neither starts nor ends on it
                        if (logIS(p))
                            return exports.PointVsFace.ON_EDGE;
                        // TODO: tangents?
                    }
                }
            }
        }
        return inside ? exports.PointVsFace.INSIDE : exports.PointVsFace.OUTSIDE;
    }
    zDirVolumeForLoop2(loop) {
        const angles = this.inverseMatrix.getZ().toAngles();
        const T = ts3dutils.M4.rotateY(-angles.theta).times(ts3dutils.M4.rotateZ(-angles.phi)).times(this.inverseMatrix);
        const rot90x = ts3dutils.M4.rotateX(PI$10 / 2);
        let totalVolume = 0;
        ts3dutils.assert(ts3dutils.V3.X.isParallelTo(T.transformVector(ts3dutils.V3.Z)));
        //const zDistanceFactor = toT.transformVector(V3.Z).length()
        loop.map(edge => edge.transform(T)).forEach((edge, edgeIndex, edges) => {
            const nextEdgeIndex = (edgeIndex + 1) % edges.length, nextEdge = edges[nextEdgeIndex];
            function f(t) {
                const at2d = edge.curve.at(t).withElement('x', 0);
                const result = 1 / 3 * (1 - (Math.pow(at2d.y, 2) + Math.pow(at2d.z, 2))) * edge.tangentAt(t).dot(rot90x.transformVector(at2d.unit()));
                return result;
            }
            //if (edge.)
            if (edge.b.like(ts3dutils.V3.X)) {
                const angleDiff = (edge.bDir.angleRelativeNormal(nextEdge.aDir, ts3dutils.V3.X) + 2 * PI$10) % (2 * PI$10);
                totalVolume += 2 / 3 * angleDiff;
            }
            if (edge.b.like(ts3dutils.V3.X.negated())) {
                const angleDiff = (edge.bDir.angleRelativeNormal(nextEdge.aDir, ts3dutils.V3.X) + 2 * PI$10) % (2 * PI$10);
                totalVolume += 2 / 3 * angleDiff;
            }
            const volume = ts3dutils.gaussLegendreQuadrature24(f, edge.aT, edge.bT);
            totalVolume += volume;
        });
        return totalVolume * this.f1.dot(this.f2.cross(this.f3));
    }
    surfaceAreaApprox() {
        // See https://en.wikipedia.org/wiki/Ellipsoid#Surface_area
        const mainAxes = this.mainAxes(), a = mainAxes.f1.length(), b = mainAxes.f2.length(), c = mainAxes.f3.length();
        const p = 1.6075;
        return 4 * PI$10 * Math.pow((Math.pow(a * b, p) + Math.pow(b * c, p) + Math.pow(c * a, p)) / 3, 1 / p);
    }
    surfaceArea() {
        // See https://en.wikipedia.org/wiki/Ellipsoid#Surface_area
        const mainAxes = this.mainAxes(), f1l = mainAxes.f1.length(), f2l = mainAxes.f2.length(), f3l = mainAxes.f3.length(), [c, b, a] = [f1l, f2l, f3l].sort(ts3dutils.MINUS);
        // https://en.wikipedia.org/w/index.php?title=Spheroid&oldid=761246800#Area
        function spheroidArea(a, c) {
            if (c < a) {
                const eccentricity2 = 1 - Math.pow(c, 2) / Math.pow(a, 2);
                const eccentricity = Math.sqrt(eccentricity2);
                return 2 * PI$10 * Math.pow(a, 2) * (1 + (1 - eccentricity2) / Math.sqrt(eccentricity) * Math.atanh(eccentricity));
            }
            else {
                const eccentricity = Math.sqrt(1 - Math.pow(a, 2) / Math.pow(c, 2));
                return 2 * PI$10 * Math.pow(a, 2) * (1 + c / a / eccentricity * Math.asin(eccentricity));
            }
        }
        if (ts3dutils.eq(a, b)) {
            return spheroidArea(a, c);
        }
        else if (ts3dutils.eq(b, c)) {
            return spheroidArea(b, a);
        }
        else if (ts3dutils.eq(c, a)) {
            return spheroidArea(c, b);
        }
        const phi = Math.acos(c / a);
        const k2 = Math.pow(a, 2) * (Math.pow(b, 2) - Math.pow(c, 2)) / (Math.pow(b, 2) * (Math.pow(a, 2) - Math.pow(c, 2)));
        const incompleteEllipticInt1 = ts3dutils.gaussLegendreQuadrature24(phi => Math.pow(1 - k2 * Math.pow(Math.sin(phi), 2), -0.5), 0, phi);
        const incompleteEllipticInt2 = ts3dutils.gaussLegendreQuadrature24(phi => Math.pow(1 - k2 * Math.pow(Math.sin(phi), 2), 0.5), 0, phi);
        return 2 * PI$10 * Math.pow(c, 2) + 2 * PI$10 * a * b / Math.sin(phi) * (incompleteEllipticInt2 * Math.pow(Math.sin(phi), 2) + incompleteEllipticInt1 * Math.pow(Math.cos(phi), 2));
    }
    getSeamPlane() {
        const plane = P3.forAnchorAndPlaneVectors(this.center, this.f1, this.f3);
        return plane.normal1.dot(this.f2) < 0 ? plane : plane.flipped();
    }
    asEllipsoidSurface() {
        return new EllipsoidSurface(this.center, this.f1, this.f2, this.f3);
    }
    getExtremePoints() {
        ts3dutils.assert(this.isSphere());
        const thisRadius = this.f1.length();
        // points on the edge of the hemisphere don't need to be included, because if they can at most be on the edge
        // of a face hemisphere can be orientated anyway, so dot with this.f2 to make sure they are "inside"
        return [ts3dutils.V3.X, ts3dutils.V3.X.negated(), ts3dutils.V3.Y, ts3dutils.V3.Y.negated(), ts3dutils.V3.Z, ts3dutils.V3.Z.negated()]
            .filter(p => ts3dutils.lt(0, p.dot(this.f2)))
            .map(p => p.times(thisRadius).plus(this.center));
    }
}
SemiEllipsoidSurface.UNIT = new SemiEllipsoidSurface(ts3dutils.V3.O, ts3dutils.V3.X, ts3dutils.V3.Y, ts3dutils.V3.Z);
SemiEllipsoidSurface.prototype.uStep = PI$10 / 16;
SemiEllipsoidSurface.prototype.vStep = PI$10 / 16;
SemiEllipsoidSurface.prototype.sMin = 0;
SemiEllipsoidSurface.prototype.sMax = PI$10;
SemiEllipsoidSurface.prototype.tMin = -PI$10 / 2;
SemiEllipsoidSurface.prototype.tMax = PI$10 / 2;

class PlaneSurface$1 extends ParametricSurface {
    constructor(plane, right = plane.normal1.getPerpendicular().unit(), up = plane.normal1.cross(right).unit(), sMin = -100, sMax = 100, tMin = -100, tMax = 100) {
        super();
        this.plane = plane;
        this.right = right;
        this.up = up;
        this.sMin = sMin;
        this.sMax = sMax;
        this.tMin = tMin;
        this.tMax = tMax;
        ts3dutils.assertInst(P3, plane);
        ts3dutils.assert(this.right.cross(this.up).like(this.plane.normal1));
        this.matrix = ts3dutils.M4.forSys(right, up, plane.normal1, plane.anchor);
    }
    toSource(rounder = x => x) {
        return ts3dutils.callsce.call(undefined, 'new PlaneSurface', ...this.getConstructorParameters());
    }
    static throughPoints(a, b, c) {
        return new PlaneSurface$1(P3.throughPoints(a, b, c));
    }
    isCoplanarTo(surface) {
        return surface instanceof PlaneSurface$1 && this.plane.isCoplanarToPlane(surface.plane);
    }
    isTsForLine(line) {
        return line.isTsWithPlane(this.plane);
    }
    like(surface) {
        return surface instanceof PlaneSurface$1 && this.plane.like(surface.plane);
    }
    pST(s, t) {
        return this.matrix.transformPoint(new ts3dutils.V3(s, t, 0));
    }
    implicitFunction() {
        return p => this.plane.distanceToPointSigned(p);
    }
    isCurvesWithSurface(surface2) {
        if (surface2 instanceof PlaneSurface$1) {
            return this.isCurvesWithPlane(surface2.plane);
        }
        return super.isCurvesWithSurface(surface2);
    }
    isCurvesWithPlane(plane) {
        if (this.plane.isParallelToPlane(plane)) {
            return [];
        }
        return [this.plane.intersectionWithPlane(plane)];
    }
    edgeLoopCCW(contour) {
        return ts3dutils.isCCW(contour.flatMap(edge => edge.points()), this.plane.normal1);
        // let totalAngle = 0
        // for (let i = 0; i < contour.length; i++) {
        // 	const ipp = (i + 1) % contour.length
        // 	const edge = contour[i], nextEdge = contour[ipp]
        // 	assert(edge.b.like(nextEdge.a), 'edges dont form a loop')
        // 	if (edge.curve instanceof SemiEllipseCurve) {
        // 		totalAngle += edge.rotViaPlane(this.plane.normal1)
        // 		// console.log(edge.toString(), edge.rotViaPlane(this.plane.normal1))
        // 	}
        // 	totalAngle += edge.bDir.angleRelativeNormal(nextEdge.aDir, this.plane.normal1)
        // }
        // const result = totalAngle > 0
        // const result2 = PlaneFace.prototype.calculateArea.apply({surface: this, contour: contour}).area > 0
        // //assert (result == result2)
        // return result2
    }
    loopContainsPoint(loop, p) {
        const dir = this.right.plus(this.up.times(0.123)).unit();
        const line = new L3$1(p, dir);
        const lineOut = dir.cross(this.plane.normal1);
        return Surface.loopContainsPointGeneral(loop, p, line, lineOut);
    }
    stPFunc() {
        const matrixInverse = this.matrix.inversed();
        return function (pWC) {
            return matrixInverse.transformPoint(pWC);
        };
    }
    pointFoot(pWC) {
        return this.stP(pWC);
    }
    normalP(pWC) {
        return this.plane.normal1;
    }
    containsPoint(p) {
        return this.plane.containsPoint(p);
    }
    containsCurve(curve) {
        return this.plane.containsCurve(curve);
    }
    transform(m4) {
        return new PlaneSurface$1(this.plane.transform(m4));
    }
    flipped() {
        return new PlaneSurface$1(this.plane.flipped(), this.right, this.up.negated());
    }
    getConstructorParameters() {
        return [this.plane, this.right, this.up];
    }
    toMesh(xMin = -10, xMax = 10, yMin = -10, yMax = 10) {
        const mesh = new tsgl.Mesh()
            .addIndexBuffer('TRIANGLES')
            .addVertexBuffer('normals', 'LGL_Normal');
        const matrix = ts3dutils.M4.forSys(this.right, this.up, this.plane.normal1, this.plane.anchor);
        mesh.vertices = [ts3dutils.V(xMin, yMin), ts3dutils.V(xMax, yMin), ts3dutils.V(xMin, yMax), ts3dutils.V(xMax, yMax)].map(p => matrix.transformPoint(p));
        mesh.normals = ts3dutils.arrayFromFunction(4, i => this.plane.normal1);
        tsgl.pushQuad(mesh.TRIANGLES, false, 0, 1, 2, 3);
        mesh.compile();
        return mesh;
    }
    dpds() {
        return () => this.right;
    }
    dpdt() {
        return () => this.up;
    }
    equals(obj) {
        return undefined;
    }
    didp(pWC) {
        return this.plane.normal1;
    }
}

const { PI: PI$11 } = Math;
const ZDirVolumeVisitor = {
    /**
     * at(t)
     * |\                                    ^
     * | \ at(t).projectedOn(dir)             \  dir
     * |  \                                    \
     * |   \ at(t).rejectedFrom(dir)
     * |   |
     * |___|
     *        z = 0
     *
     *
     * A = ((at(t) + at(t).rejectedFrom(dir)) / 2).z * at(t).projectedOn(dir).lengthXY()
     * scaling = tangentAt(t) DOT dir.cross(V3.Z).unit()
     */
    [ConicSurface.name](allEdges) {
        // INT[edge.at; edge.bT] (at(t) DOT dir) * (at(t) - at(t).projectedOn(dir) / 2).z
        const totalVolume = allEdges.map(edge => {
            if (edge.curve instanceof SemiEllipseCurve || edge.curve instanceof HyperbolaCurve || edge.curve instanceof ParabolaCurve) {
                const f = (t) => {
                    const at = edge.curve.at(t), tangent = edge.tangentAt(t);
                    return (at.z + at.rejectedFrom(this.dir).z) / 2 * at.projectedOn(this.dir).lengthXY() *
                        tangent.dot(ts3dutils.V3.Z.cross(this.dir).unit());
                };
                // ellipse with normal1 parallel to dir need to be counted negatively so CCW faces result in a positive
                // area
                const sign = edge.curve instanceof SemiEllipseCurve
                    ? -Math.sign(edge.curve.normal.dot(this.dir))
                    : -Math.sign(this.center.to(edge.curve.center).cross(edge.curve.f1).dot(this.dir));
                const val = ts3dutils.glqInSteps(f, edge.aT, edge.bT, 1);
                return val * sign;
            }
            else if (edge.curve instanceof L3) {
                return 0;
            }
            else {
                ts3dutils.assertNever();
            }
        }).sum();
        return { volume: totalVolume * Math.sign(this.normal.dot(this.dir)) };
    },
    [PlaneSurface$1.name]() {
        const { centroid, area } = this.calculateArea();
        return {
            volume: this.surface.plane.normal1.z * centroid.z * area,
            centroid: new ts3dutils.V3(centroid.x, centroid.y, centroid.z / 2),
        };
    },
    /**
     * at(t)
     * |\                                    ^
     * | \ at(t).projectedOn(dir1)            \  dir1
     * |  \                                    \
     * |   \ at(t).rejectedFrom(dir1)
     * |   |
     * |___|
     *        z = 0
     *
     *
     * A = ((at(t) + at(t).rejectedFrom(dir1)) / 2).z * at(t).projectedOn(dir1).lengthXY()
     * scaling = tangentAt(t) DOT dir1.cross(V3.Z).unit()
     */
    [SemiCylinderSurface.name](allEdges) {
        if (ts3dutils.V3.Z.cross(this.dir).likeO())
            return { volume: 0 };
        // the tangent needs to be projected onto a vector which is perpendicular to the volume-slices
        const scalingVector = this.dir.cross(ts3dutils.V3.Z).unit();
        // the length of the base of the trapezoid is calculated by dotting with the baseVector
        const baseVector = this.dir.rejectedFrom(ts3dutils.V3.Z).unit();
        // INT[edge.at; edge.bT] (at(t) DOT dir1) * (at(t) - at(t).projectedOn(dir) / 2).z
        const totalArea = edges.map(edge => {
            if (edge.curve instanceof SemiEllipseCurve) {
                const f = (t) => {
                    // use curve.tangent not edge.tangent, reverse edges are handled by the integration boundaries
                    const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t);
                    const area = (at.z + at.rejectedFrom(this.dir).z) / 2 * at.projectedOn(this.dir).dot(baseVector);
                    const scale = tangent.dot(scalingVector);
                    //assert(Math.sign(scale) == Math.sign(this.normalP(at).dot(V3.Z)), this.normalP(at).dot(V3.Z))
                    //console.log(
                    //	"", t,
                    //	",", area,
                    //	",", scale,
                    //	"atz", at.z)
                    return area * scale;
                };
                // ellipse with normal parallel to dir1 need to be counted negatively so CCW faces result in a positive
                // area
                const sign = -Math.sign(edge.curve.normal.dot(this.dir));
                const val = ts3dutils.glqInSteps(f, edge.aT, edge.bT, 1);
                return val * sign;
            }
            else if (edge.curve instanceof L3) {
                return 0;
            }
            else {
                ts3dutils.assertNever();
            }
        }).sum();
        return { volume: totalArea * Math.sign(this.baseCurve.normal.dot(this.dir)) };
    },
    /**
     * at(t)
     * |\                                    ^
     * | \ at(t).projectedOn(dir1)            \  dir1
     * |  \                                    \
     * |   \ at(t).rejectedFrom(dir1)
     * |   |
     * |___|
     *        z = 0
     *
     *
     * A = ((at(t) + at(t).rejectedFrom(dir1)) / 2).z * at(t).projectedOn(dir1).lengthXY()
     * scaling = tangentAt(t) DOT dir1.cross(V3.Z).unit()
     */
    [CylinderSurface.name](edges) {
        if (ts3dutils.V3.Z.cross(this.dir).likeO())
            return { volume: 0 };
        // the tangent needs to be projected onto a vector which is perpendicular to the volume-slices
        const scalingVector = this.dir.cross(ts3dutils.V3.Z).unit();
        // the length of the base of the trapezoid is calculated by dotting with the baseVector
        const baseVector = this.dir.rejectedFrom(ts3dutils.V3.Z).unit();
        // INT[edge.at; edge.bT] (at(t) DOT dir1) * (at(t) - at(t).projectedOn(dir) / 2).z
        console.log('scalingVector', scalingVector.sce);
        const totalArea = edges.map(edge => {
            if (edge.curve instanceof EllipseCurve) {
                const f = (t) => {
                    // use curve.tangent not edge.tangent, reverse edges are handled by the integration boundaries
                    const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t);
                    const area = (at.z + at.rejectedFrom(this.dir).z) / 2 * at.projectedOn(this.dir).dot(baseVector);
                    const scale = tangent.dot(scalingVector);
                    //assert(Math.sign(scale) == Math.sign(this.normalP(at).dot(V3.Z)), this.normalP(at).dot(V3.Z))
                    //console.log(
                    //	"", t,
                    //	",", area,
                    //	",", scale,
                    //	"atz", at.z)
                    return area * scale;
                };
                // ellipse with normal parallel to dir1 need to be counted negatively so CCW faces result in a positive
                // area
                const sign = -Math.sign(edge.curve.normal.dot(this.dir));
                const val = ts3dutils.glqInSteps(f, edge.aT, edge.bT, 1);
                console.log('edge', edge, val, sign);
                return val * sign;
            }
            else if (edge.curve instanceof L3) {
                return 0;
            }
            else {
                ts3dutils.assertNever();
            }
        }).sum();
        return { volume: totalArea * Math.sign(this.baseCurve.normal.dot(this.dir)) };
    },
    // volume does scale linearly, so this can be done in the local coordinate system
    // first transform edges with inverse matrix
    // then rotate everything edges so the original world Z dir again points in Z dir
    // now we have a problem because edges which originally  did not cross the seam plane can now be anywhere
    // we need to split the transformed loop along the local seam plane
    // and then sum the zDir volumes of the resulting loops
    [EllipsoidSurface.name](loop) {
        const angles = this.inverseMatrix.transformVector(ts3dutils.V3.Z).toAngles();
        const T = ts3dutils.M4.rotateAB(this.inverseMatrix.transformVector(ts3dutils.V3.Z), ts3dutils.V3.Z).times(ts3dutils.M4.rotateZ(-angles.phi)).times(this.inverseMatrix);
        function calc(loop) {
            let totalVolume = 0;
            assert(ts3dutils.V3.Z.isParallelTo(T.transformVector(ts3dutils.V3.Z)));
            //const zDistanceFactor = toT.transformVector(V3.Z).length()
            loop.map(edge => edge.transform(T)).forEach((edge, edgeIndex, edges) => {
                const nextEdgeIndex = (edgeIndex + 1) % edges.length, nextEdge = edges[nextEdgeIndex];
                function f(t) {
                    const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t);
                    const r = at.lengthXY();
                    const at2d = at.withElement('z', 0);
                    const angleAdjusted = (at.angleXY() + ts3dutils.TAU - ts3dutils.NLA_PRECISION) % ts3dutils.TAU + ts3dutils.NLA_PRECISION;
                    const result = angleAdjusted * Math.sqrt(1 - r * r) * r * Math.abs(tangent.dot(at2d.unit())) * Math.sign(tangent.z);
                    //console.log("at2d", at2d.sce, "result", result, 'angle', angleAdjusted, '
                    // edge.tangentAt(t).dot(at2d.unit())', edge.tangentAt(t).dot(at2d.unit()))
                    return result;
                }
                const volume = ts3dutils.gaussLegendreQuadrature24(f, edge.aT, edge.bT);
                console.log('edge', edge, 'volume', volume);
                totalVolume += volume;
            });
            return totalVolume;
        }
        const [front, back] = EllipsoidSurface.splitOnPlaneLoop(loop.map(edge => edge.transform(T)), ccw);
        const localVolume = calc(front, PI$11) + calc(back, -PI$11);
        return { area: localVolume * this.f1.dot(this.f2.cross(this.f3)), centroid: undefined };
    },
    zDirVolumeForLoop2(loop) {
        const angles = this.inverseMatrix.getZ().toAngles();
        const T = ts3dutils.M4.rotateY(-angles.theta).times(ts3dutils.M4.rotateZ(-angles.phi)).times(this.inverseMatrix);
        const rot90x = ts3dutils.M4.rotateX(PI$11 / 2);
        let totalVolume = 0;
        assert(ts3dutils.V3.X.isParallelTo(T.transformVector(ts3dutils.V3.Z)));
        //const zDistanceFactor = toT.transformVector(V3.Z).length()
        loop.map(edge => edge.transform(T)).forEach((edge, edgeIndex, edges) => {
            const nextEdgeIndex = (edgeIndex + 1) % edges.length, nextEdge = edges[nextEdgeIndex];
            function f(t) {
                const at2d = edge.curve.at(t).withElement('x', 0);
                const result = 1 / 3 * (1 - (Math.pow(at2d.y, 2) + Math.pow(at2d.z, 2))) * edge.tangentAt(t).dot(rot90x.transformVector(at2d.unit()));
                console.log('at2d', at2d.sce, 'result', result);
                return result;
            }
            //if (edge.)
            if (edge.b.like(ts3dutils.V3.X)) {
                const angleDiff = (edge.bDir.angleRelativeNormal(nextEdge.aDir, ts3dutils.V3.X) + 2 * PI$11) % (2 * PI$11);
                totalVolume += 2 / 3 * angleDiff;
                console.log('xaa');
            }
            if (edge.b.like(ts3dutils.V3.X.negated())) {
                const angleDiff = (edge.bDir.angleRelativeNormal(nextEdge.aDir, ts3dutils.V3.X) + 2 * PI$11) % (2 * PI$11);
                totalVolume += 2 / 3 * angleDiff;
                console.log('xbb');
            }
            const volume = ts3dutils.gaussLegendreQuadrature24(f, edge.aT, edge.bT);
            console.log('edge', edge, 'volume', volume);
            totalVolume += volume;
        });
        return totalVolume * this.f1.dot(this.f2.cross(this.f3));
    },
    // volume does scale linearly, so this can be done in the local coordinate system
    // first transform edges with inverse matrix
    // then rotate everything edges so the original world Z dir again points in Z dir
    // now we have a problem because edges which originally  did not cross the seam plane can now be anywhere
    // we need to split the transformed loop along the local seam plane
    // and then sum the zDir volumes of the resulting loops
    [SemiEllipsoidSurface.name](loop) {
        const angles = this.inverseMatrix.transformVector(ts3dutils.V3.Z).toAngles();
        const T = ts3dutils.M4.rotateAB(this.inverseMatrix.transformVector(ts3dutils.V3.Z), ts3dutils.V3.Z).times(ts3dutils.M4.rotateZ(-angles.phi)).times(this.inverseMatrix);
        function calc(loop) {
            let totalVolume = 0;
            assert(ts3dutils.V3.Z.isParallelTo(T.transformVector(ts3dutils.V3.Z)));
            //const zDistanceFactor = toT.transformVector(V3.Z).length()
            loop.map(edge => edge.transform(T)).forEach((edge, edgeIndex, edges) => {
                const nextEdgeIndex = (edgeIndex + 1) % edges.length, nextEdge = edges[nextEdgeIndex];
                function f(t) {
                    const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t);
                    const r = at.lengthXY();
                    const at2d = at.withElement('z', 0);
                    const angleAdjusted = (at.angleXY() + ts3dutils.TAU - ts3dutils.NLA_PRECISION) % ts3dutils.TAU + ts3dutils.NLA_PRECISION;
                    const result = angleAdjusted * Math.sqrt(1 - r * r) * r * Math.abs(tangent.dot(at2d.unit())) * Math.sign(tangent.z);
                    //console.log("at2d", at2d.sce, "result", result, 'angle', angleAdjusted, '
                    // edge.tangentAt(t).dot(at2d.unit())', edge.tangentAt(t).dot(at2d.unit()))
                    return result;
                }
                const volume = ts3dutils.gaussLegendreQuadrature24(f, edge.aT, edge.bT);
                totalVolume += volume;
            });
            return totalVolume;
        }
        const [front, back] = SemiEllipsoidSurface.splitOnPlaneLoop(loop.map(edge => edge.transform(T)), ccw);
        const localVolume = calc(front, PI$11) + calc(back, -PI$11);
        return { volume: localVolume * this.f1.dot(this.f2.cross(this.f3)), centroid: undefined };
    },
};

const { PI: PI$12 } = Math;
const CalculateAreaVisitor = {
    [ConicSurface.name](edges) {
        // calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
        const totalArea = edges.map(edge => {
            if (edge.curve instanceof SemiEllipseCurve || edge.curve instanceof HyperbolaCurve || edge.curve instanceof ParabolaCurve) {
                const f = (t) => {
                    const at = edge.curve.at(t), tangent = edge.tangentAt(t);
                    return at.minus(this.center).cross(tangent.rejectedFrom(this.dir)).length() / 2;
                };
                // ellipse with normal1 parallel to dir1 need to be counted negatively so CCW faces result in a
                // positive area hyperbola normal1 can be perpendicular to
                const sign = edge.curve instanceof SemiEllipseCurve
                    ? -Math.sign(edge.curve.normal.dot(this.dir))
                    : -Math.sign(this.center.to(edge.curve.center).cross(edge.curve.f1).dot(this.dir));
                return ts3dutils.glqInSteps(f, edge.aT, edge.bT, 4) * sign;
            }
            else if (edge.curve instanceof L3$1) {
                return 0;
            }
            else {
                ts3dutils.assertNever();
            }
        }).sum();
        // if the cylinder faces inwards, CCW faces will have been CW, so we need to reverse that here
        // Math.abs is not an option as "holes" may also be passed
        return totalArea * Math.sign(this.normal.dot(this.dir));
    },
    [PlaneSurface$1.name](edges) {
        let centroid = ts3dutils.V3.O, tcs = 0, tct = 0, totalArea = 0;
        let r1 = this.surface.right, u1 = this.surface.up;
        for (const edge of edges) {
            let edgeCentroid, edgeArea, centroidS, centroidT;
            if (edge instanceof StraightEdge) {
                const midPoint = edge.a.lerp(edge.b, 0.5);
                edgeCentroid = new ts3dutils.V3(midPoint.x, centroid.y, centroid.z / 2);
                centroidS = midPoint.dot(r1) / 2;
                centroidT = midPoint.dot(u1);
                const edgeLength = edge.a.distanceTo(edge.b);
                edgeArea = edgeLength * edge.curve.dir1.dot(r1);
                edgeArea = (edge.a.dot(u1) + edge.b.dot(u1)) / 2 * edge.b.to(edge.a).dot(r1);
            }
            else {
                let curve = edge.curve;
                if (curve instanceof SemiEllipseCurve) {
                    let info = curve.getAreaInDir(r1, u1, edge.aT, edge.bT);
                    edgeArea = info.area;
                    let parametricCentroid = this.surface.stPFunc()(info.centroid);
                    centroidS = parametricCentroid.x;
                    centroidT = parametricCentroid.y;
                }
                else if (curve instanceof BezierCurve) {
                    edgeArea = curve.getAreaInDirSurface(u1, this.surface, edge.aT, edge.bT).area;
                }
                else {
                    ts3dutils.assertNever();
                }
            }
            tcs += edgeArea * centroidS;
            tct += edgeArea * centroidT;
            totalArea += edgeArea;
        }
        centroid = r1.times(tcs).plus(u1.times(tct));
        ts3dutils.assert(isFinite(totalArea));
        return { area: totalArea, centroid: centroid };
    },
    /**
     * Calculating the surface area of a projected ellipse is analogous to the circumference of the ellipse
     * ==> Elliptic integrals/numeric calculation is necessary
     */
    [CylinderSurface.name](edges) {
        // calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
        const totalArea = edges.map(edge => {
            if (edge.curve instanceof EllipseCurve) {
                const f = (t) => {
                    const at = edge.curve.at(t), tangent = edge.tangentAt(t);
                    return at.dot(this.dir) * tangent.rejected1Length(this.dir);
                };
                // ellipse with normal parallel to dir1 need to be counted negatively so CCW faces result in a positive
                // area
                const sign = -Math.sign(edge.curve.normal.dot(this.dir));
                const val = ts3dutils.glqInSteps(f, edge.aT, edge.bT, 4);
                console.log('edge', edge, val);
                return val * sign;
            }
            else if (edge.curve instanceof L3$1) {
                return 0;
            }
            else {
                ts3dutils.assertNever();
            }
        }).sum();
        // if the cylinder faces inwards, CCW faces will have been CW, so we need to reverse that here
        // Math.abs is not an option as "holes" may also be passed
        return totalArea * Math.sign(this.baseCurve.normal.dot(this.dir));
    },
    [EllipsoidSurface.name](edges, canApproximate = true) {
        ts3dutils.assert(this.isVerticalSpheroid());
        const { f1, f2, f3 } = this;
        // calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
        const circleRadius = f1.length();
        const f31 = f3.unit();
        const totalArea = edges.map(edge => {
            if (edge.curve instanceof EllipseCurve) {
                const f = (t) => {
                    const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t);
                    const localAt = this.inverseMatrix.transformPoint(at);
                    let angleXY = localAt.angleXY();
                    if (ts3dutils.eq(Math.abs(angleXY), PI$12)) {
                        if (edge.curve.normal.isParallelTo(this.f2)) {
                            angleXY = PI$12 * -Math.sign((edge.bT - edge.aT) * edge.curve.normal.dot(this.f2));
                        }
                        else {
                            angleXY = PI$12 * dotCurve(this.f2, tangent, edge.curve.ddt(t));
                        }
                        console.log(angleXY);
                    }
                    const arcLength = angleXY * circleRadius * Math.sqrt(1 - Math.pow(localAt.z, 2));
                    const dotter = this.matrix.transformVector(new ts3dutils.V3(-localAt.z * localAt.x / localAt.lengthXY(), -localAt.z * localAt.y / localAt.lengthXY(), localAt.lengthXY())).unit();
                    const df3 = tangent.dot(f31);
                    //const scaling = df3 / localAt.lengthXY()
                    const scaling = dotter.dot(tangent);
                    //console.log(t, at.str, arcLength, scaling)
                    return arcLength * scaling;
                };
                const val = ts3dutils.glqInSteps(f, edge.aT, edge.bT, 1);
                console.log('edge', edge, val);
                return val;
            }
            else {
                ts3dutils.assertNever();
            }
        }).sum();
        return totalArea * Math.sign(this.f1.cross(this.f2).dot(this.f3));
    },
    [SemiCylinderSurface.name](edges, canApproximate = true) {
        ts3dutils.assert(this.isVerticalSpheroid());
        const { f1, f2, f3 } = this;
        // calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
        const circleRadius = f1.length();
        const f31 = f3.unit();
        const totalArea = edges.map(edge => {
            if (edge.curve instanceof SemiEllipseCurve) {
                const f = (t) => {
                    const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t);
                    const localAt = this.inverseMatrix.transformPoint(at);
                    let angleXY = localAt.angleXY();
                    if (ts3dutils.eq(Math.abs(angleXY), PI$12)) {
                        if (edge.curve.normal.isParallelTo(this.f2)) {
                            angleXY = PI$12 * -Math.sign((edge.bT - edge.aT) * edge.curve.normal.dot(this.f2));
                        }
                        else {
                            angleXY = PI$12 * dotCurve(this.f2, tangent, edge.curve.ddt(t));
                        }
                    }
                    const arcLength = angleXY * circleRadius * Math.sqrt(1 - Math.pow(localAt.z, 2));
                    const dotter = this.matrix.transformVector(new ts3dutils.V3(-localAt.z * localAt.x / localAt.lengthXY(), -localAt.z * localAt.y / localAt.lengthXY(), localAt.lengthXY())).unit();
                    const df3 = tangent.dot(f31);
                    //const scaling = df3 / localAt.lengthXY()
                    const scaling = dotter.dot(tangent);
                    //console.log(t, at.str, arcLength, scaling)
                    return arcLength * scaling;
                };
                const val = ts3dutils.glqInSteps(f, edge.aT, edge.bT, 1);
                return val;
            }
            else {
                ts3dutils.assertNever();
            }
        }).sum();
        return totalArea * Math.sign(this.f1.cross(this.f2).dot(this.f3));
    },
    [ProjectedCurveSurface.name](edges) {
        // calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
        const totalArea = edges.map(edge => {
            if (edge.curve instanceof SemiEllipseCurve) {
                const f = (t) => {
                    const at = edge.curve.at(t), tangent = edge.tangentAt(t);
                    return at.dot(this.dir) * tangent.rejected1Length(this.dir);
                };
                // ellipse with normal1 parallel to dir1 need to be counted negatively so CCW faces result in a
                // positive area
                const sign = -Math.sign(edge.curve.normal.dot(this.dir));
                const val = ts3dutils.glqInSteps(f, edge.aT, edge.bT, 4);
                return val * sign;
            }
            else if (edge.curve instanceof L3$1) {
                return 0;
            }
            else {
                ts3dutils.assertNever();
            }
        }).sum();
        // if the cylinder faces inwards, CCW faces will have been CW, so we need to reverse that here
        // Math.abs is not an option as "holes" may also be passed
        return totalArea * Math.sign(this.baseCurve.normal.dot(this.dir));
    },
    /**
     * Calculating the surface area of a projected ellipse is analogous to the circumference of the ellipse
     * ==> Elliptic integrals/numeric calculation is necessary
     */
    [SemiCylinderSurface.name](edges) {
        // calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
        const totalArea = edges.map(edge => {
            if (edge.curve instanceof SemiEllipseCurve) {
                const f = (t) => {
                    const at = edge.curve.at(t), tangent = edge.tangentAt(t);
                    return at.dot(this.dir) * tangent.rejected1Length(this.dir);
                };
                // ellipse with normal parallel to dir1 need to be counted negatively so CCW faces result in a positive
                // area
                const sign = -Math.sign(edge.curve.normal.dot(this.dir));
                const val = ts3dutils.glqInSteps(f, edge.aT, edge.bT, 4);
                return val * sign;
            }
            else if (edge.curve instanceof L3$1) {
                return 0;
            }
            else {
                ts3dutils.assertNever();
            }
        }).sum();
        // if the cylinder faces inwards, CCW faces will have been CW, so we need to reverse that here
        // Math.abs is not an option as "holes" may also be passed
        return totalArea * Math.sign(this.baseCurve.normal.dot(this.dir));
    },
};

const { PI: PI$13, min: min$9, max: max$8, ceil: ceil$10 } = Math;
function projectCurve(curve, offset, flipped) {
    if (curve instanceof L3$1) {
        const surfaceNormal = offset.cross(curve.dir1).toLength(flipped ? -1 : 1);
        return new PlaneSurface$1(P3.normalOnAnchor(surfaceNormal, curve.anchor));
    }
    if (curve instanceof SemiEllipseCurve) {
        const curveDir = flipped ? offset : offset.negated();
        return new SemiCylinderSurface(curve, curveDir.unit(), undefined, undefined);
    }
    if (curve instanceof BezierCurve || curve instanceof XiEtaCurve) {
        const curveDir = offset.times(flipped ? 1 : -1);
        return new ProjectedCurveSurface(curve, curveDir, 0, 1, flipped ? 0 : -1, flipped ? 1 : 0);
    }
    throw new Error();
}

(function (B2T) {
    function box(w = 1, h = 1, d = 1, name) {
        ts3dutils.assertNumbers(w, h, d);
        ts3dutils.assertInst('string' === typeof name);
        const baseVertices = [
            new ts3dutils.V3(0, 0, 0),
            new ts3dutils.V3(0, h, 0),
            new ts3dutils.V3(w, h, 0),
            new ts3dutils.V3(w, 0, 0),
        ];
        const generator = ts3dutils.callsce('B2T.box', w, h, d, name);
        return B2T.extrudeVertices(baseVertices, P3.XY.flipped(), new ts3dutils.V3(0, 0, d), name, generator);
    }
    B2T.box = box;
    function puckman(radius, rads, height, name) {
        ts3dutils.assertf(() => ts3dutils.lt(0, radius));
        ts3dutils.assertf(() => ts3dutils.lt(0, rads) && ts3dutils.le(rads, ts3dutils.TAU));
        ts3dutils.assertf(() => ts3dutils.lt(0, height));
        const edges = StraightEdge.chain([ts3dutils.V3.O, new ts3dutils.V3(radius, 0, 0), new ts3dutils.V3(radius, 0, height), new ts3dutils.V3(0, 0, height)], true);
        return B2T.rotateEdges(edges, rads, name || 'puckman' + getGlobalId());
    }
    B2T.puckman = puckman;
    function registerVertexName(map, name, p) {
        // TODO
        if (!Array.from(map.keys()).some(p2 => p2.like(p))) {
            map.set(p, name);
        }
    }
    B2T.registerVertexName = registerVertexName;
    function extrudeEdges(baseFaceEdges, baseFacePlane = P3.XY, offset = ts3dutils.V3.Z, name = 'extrude' + getGlobalId(), gen, infoFactory) {
        baseFaceEdges = fixEdges(baseFaceEdges);
        //Array.from(combinations(baseFaceEdges.length)).forEach(({i, j}) => {
        //	assertf(() => !Edge.edgesIntersect(baseFaceEdges[i], baseFaceEdges[j]), baseFaceEdges[i].sce +
        // baseFaceEdges[j].sce) })
        ts3dutils.assertf(() => Edge.isLoop(baseFaceEdges));
        // TODO checks..
        //if (offset.dot(baseFacePlane.normal1) > 0) {
        //	baseFacePlane = baseFacePlane.flipped()
        //}
        const vertexNames = new Map();
        const basePlaneSurface = new PlaneSurface$1(baseFacePlane);
        //assert(basePlaneSurface.edgeLoopCCW(baseFaceEdges), 'edges not CCW on baseFacePlane')
        const translationMatrix = ts3dutils.M4.translate(offset);
        const topEdges = baseFaceEdges.map(edge => edge.transform(translationMatrix, 'top'));
        const edgeCount = baseFaceEdges.length;
        const bottomInfo = infoFactory && infoFactory.extrudeBottom(basePlaneSurface, baseFaceEdges);
        const bottomFace = new PlaneFace(basePlaneSurface, baseFaceEdges, [], name + 'Bottom', bottomInfo);
        const topFaceEdges = topEdges.map(edge => edge.flipped()).reverse();
        const topSurface = new PlaneSurface$1(baseFacePlane.flipped().translated(offset));
        const topInfo = infoFactory && infoFactory.extrudeBottom(topSurface, topFaceEdges);
        const topFace = new PlaneFace(topSurface, topFaceEdges, [], name + 'Top', topInfo);
        baseFaceEdges.forEach(edge => B2T.registerVertexName(vertexNames, edge.name + 'A', edge.a));
        topFaceEdges.forEach(edge => B2T.registerVertexName(vertexNames, edge.name + 'A', edge.a));
        const ribs = ts3dutils.arrayFromFunction(edgeCount, i => StraightEdge.throughPoints(baseFaceEdges[i].a, topEdges[i].a, name + 'Rib' + i));
        const faces = baseFaceEdges.map((edge, i) => {
            const faceName = name + 'Wall' + i;
            const j = (i + 1) % edgeCount;
            const faceEdges = [baseFaceEdges[i].flipped(), ribs[i], topEdges[i], ribs[j].flipped()];
            const surface = projectCurve(edge.curve, offset, edge.reversed);
            const info = infoFactory && infoFactory.extrudeWall(i, surface, faceEdges);
            return Face.create(surface, faceEdges, undefined, faceName, info);
        });
        faces.push(bottomFace, topFace);
        gen = gen || ts3dutils.callsce('B2T.extrudeEdges', baseFaceEdges, baseFacePlane, offset, name);
        return new B2(faces, baseFacePlane.normal1.dot(offset) > 0, gen, vertexNames);
    }
    B2T.extrudeEdges = extrudeEdges;
    function cylinder(radius = 1, height = 1, rads = ts3dutils.TAU, name = 'cylinder' + getGlobalId()) {
        const vertices = [new ts3dutils.V3(0, 0, 0), new ts3dutils.V3(radius, 0, 0), new ts3dutils.V3(radius, 0, height), new ts3dutils.V3(0, 0, height)];
        return rotateEdges(StraightEdge.chain(vertices, true), rads, name);
    }
    B2T.cylinder = cylinder;
    function cone(radius = 1, height = 1, rads = ts3dutils.TAU, name = 'cone' + getGlobalId()) {
        const vertices = [new ts3dutils.V3(0, 0, 0), new ts3dutils.V3(radius, 0, height), new ts3dutils.V3(0, 0, height)];
        return rotateEdges(StraightEdge.chain(vertices, true), rads, name);
    }
    B2T.cone = cone;
    function sphere(radius = 1, name = 'sphere' + getGlobalId(), rot = ts3dutils.TAU) {
        const ee = PCurveEdge.create(new SemiEllipseCurve(ts3dutils.V3.O, new ts3dutils.V3(0, 0, -radius), new ts3dutils.V3(radius, 0, 0)), new ts3dutils.V3(0, 0, -radius), new ts3dutils.V3(0, 0, radius), 0, PI$13, undefined, new ts3dutils.V3(radius, 0, 0), new ts3dutils.V3(-radius, 0, 0));
        const generator = ts3dutils.callsce('B2T.sphere', radius, name, rot);
        return rotateEdges([StraightEdge.throughPoints(ee.b, ee.a), ee], rot, name, generator);
    }
    B2T.sphere = sphere;
    function menger(res = 2, name = 'menger' + getGlobalId()) {
        let result = B2T.box(1, 1, 1);
        if (0 == res)
            return result;
        const punch = B2T.box(1 / 3, 1 / 3, 2).translate(1 / 3, 1 / 3, -1 / 2).flipped();
        function recurse(steps, m4) {
            result = result.and(punch.transform(m4));
            a = result;
            if (steps > 1) {
                const scaled = m4.times(ts3dutils.M4.scale(1 / 3, 1 / 3, 1));
                for (let i = 0; i < 9; i++) {
                    if (4 == i)
                        continue;
                    recurse(steps - 1, scaled.times(ts3dutils.M4.translate(i % 3, i / 3 | 0, 0)));
                }
            }
        }
        recurse(res, ts3dutils.M4.IDENTITY);
        recurse(res, ts3dutils.M4.YZX);
        recurse(res, ts3dutils.M4.ZXY);
        return result;
    }
    B2T.menger = menger;
    function menger2(res = 2, name = 'menger' + getGlobalId()) {
        if (0 == res)
            return B2T.box(1, 1, 1);
        const punch = B2T.box(1 / 3, 1 / 3, 2).translate(1 / 3, 1 / 3, -1 / 2).flipped();
        const stencilFaces = [];
        function recurse(steps, m4) {
            stencilFaces.push(...punch.transform(m4).faces);
            if (steps > 1) {
                const scaled = m4.times(ts3dutils.M4.scale(1 / 3, 1 / 3, 1));
                for (let i = 0; i < 9; i++) {
                    if (4 == i)
                        continue;
                    recurse(steps - 1, scaled.times(ts3dutils.M4.translate(i % 3, i / 3 | 0, 0)));
                }
            }
        }
        recurse(res, ts3dutils.M4.IDENTITY);
        const stencil = new B2(stencilFaces, true);
        return B2T.box()
            .and(stencil)
            .and(stencil.transform(ts3dutils.M4.YZX))
            .and(stencil.transform(ts3dutils.M4.ZXY));
    }
    B2T.menger2 = menger2;
    function torus(rSmall, rLarge, rads, name) {
        ts3dutils.assertNumbers(rSmall, rLarge, rads);
        ts3dutils.assertf(() => rLarge > rSmall);
        const curve = SemiEllipseCurve.semicircle(rSmall, new ts3dutils.V3(rLarge, 0, 0));
        const baseEdges = [PCurveEdge.forCurveAndTs(curve, -Math.PI, 0), PCurveEdge.forCurveAndTs(curve, 0, Math.PI)];
        return B2T.rotateEdges(baseEdges, rads, name || 'torus' + getGlobalId());
    }
    B2T.torus = torus;
    function torusUnsplit(rSmall, rLarge, rads, name) {
        ts3dutils.assertNumbers(rSmall, rLarge, rads);
        ts3dutils.assertf(() => rLarge > rSmall);
        const baseEdge = PCurveEdge.forCurveAndTs(SemiEllipseCurve.semicircle(rSmall, new ts3dutils.V3(rLarge, 0, 0)), -Math.PI, Math.PI);
        return B2T.rotateEdges([baseEdge], rads, name || 'torus' + getGlobalId());
    }
    B2T.torusUnsplit = torusUnsplit;
    /**
     * baseLoop should be CCW on XZ plane for a bounded B2
     */
    function rotateEdges(baseLoop, totalRads, name, generator, infoFactory) {
        ts3dutils.assert(!ts3dutils.eq(PI$13, totalRads) || PI$13 == totalRads); // URHGJ
        ts3dutils.assertf(() => ts3dutils.lt(0, totalRads) && ts3dutils.le(totalRads, ts3dutils.TAU));
        totalRads = ts3dutils.snap(totalRads, ts3dutils.TAU);
        const basePlane = new PlaneSurface$1(P3.ZX.flipped()).edgeLoopCCW(baseLoop)
            ? new PlaneSurface$1(P3.ZX.flipped())
            : new PlaneSurface$1(P3.ZX);
        ts3dutils.assertf(() => Edge.isLoop(baseLoop));
        const rotationSteps = ceil$10((totalRads - ts3dutils.NLA_PRECISION) / PI$13);
        const open = !ts3dutils.eq(totalRads, 2 * PI$13);
        const baseRibCurves = baseLoop.map(edge => {
            const a = edge.a, radius = a.lengthXY();
            if (!ts3dutils.eq0(radius)) {
                return new SemiEllipseCurve(ts3dutils.V(0, 0, a.z), ts3dutils.V(radius, 0, 0), ts3dutils.V(0, radius, 0));
            }
        });
        const baseSurfaces = baseLoop.map((edge, i) => {
            const ipp = (i + 1) % baseLoop.length;
            if (edge instanceof StraightEdge) {
                const line = edge.curve;
                if (line.dir1.isParallelTo(ts3dutils.V3.Z)) {
                    if (ts3dutils.eq0(edge.a.x)) {
                        return;
                    }
                    const flipped = edge.a.z > edge.b.z;
                    const [tMin, tMax] = [0, edge.b.z - edge.a.z].sort(ts3dutils.MINUS);
                    return new SemiCylinderSurface(baseRibCurves[i], !flipped ? ts3dutils.V3.Z : ts3dutils.V3.Z.negated(), undefined, undefined, tMin, tMax);
                }
                else if (line.dir1.isPerpendicularTo(ts3dutils.V3.Z)) {
                    const flipped = edge.a.x > edge.b.x;
                    let surface = new PlaneSurface$1(new P3(ts3dutils.V3.Z, edge.a.z));
                    if (!flipped)
                        surface = surface.flipped();
                    return surface;
                }
                else {
                    // apex is intersection of segment with Z-axis
                    const a = edge.a, b = edge.b;
                    const apexZ = a.z - a.x * (b.z - a.z) / (b.x - a.x);
                    const apex = new ts3dutils.V3(0, 0, apexZ);
                    const flipped = edge.a.z > edge.b.z;
                    const base = baseRibCurves[a.x > b.x ? i : ipp];
                    const surface = ConicSurface.atApexThroughEllipse(apex, base);
                    return flipped != (-1 == surface.normalDir) ? surface.flipped() : surface;
                }
            }
            /*
                at(t) = f1 * cos t + f2 sin t
                rotated projection
                at2(t) = V(at(t).lengthXY(), 0, at(t).z)
                at2(t).x = sqrt((f1x cos t + f2x sin t)² + (f1y cos t + f2y sin t)²)
                at2(t).x = sqrt((f1x² + f1y²) cos² t + (f1x f2x + f1y f2y) cos t sin t + (f2x² + f2y²)sin²t)
                at2(t).x = sqrt((a² + b²) cos² t + (a c + b d) cos t sin t + (c² + d²)sin²t)
                (x cos t + y sin t)² = x² cos² t + x y cos t sin t + y² sin² t
             */
            if (edge.curve instanceof SemiEllipseCurve) {
                const flipped = edge.a.z > edge.b.z;
                const ell = edge.curve.rightAngled();
                ts3dutils.assert(ell.normal.isPerpendicularTo(ts3dutils.V3.Z));
                ts3dutils.assert(L3$1.Z.containsPoint(ell.center));
                let width = ell.f1.length(), height = ell.f2.length();
                if (!ell.isCircular()) {
                    ts3dutils.assert(ell.f1.isParallelTo(ts3dutils.V3.Z) && ell.f2.isParallelTo(ts3dutils.V3.X)
                        || ell.f2.isParallelTo(ts3dutils.V3.Z) && ell.f1.isParallelTo(ts3dutils.V3.X));
                    if (ell.f1.isParallelTo(ts3dutils.V3.Z)) {
                        [width, height] = [height, width];
                    }
                }
                return SemiEllipsoidSurface.forABC(width, (!flipped ? 1 : -1) * width, height, ell.center);
            }
            else {
                ts3dutils.assert(false, edge);
            }
        });
        let stepStartEdges = baseLoop, stepEndEdges;
        const faces = [];
        for (let rot = 0; rot < totalRads; rot += PI$13) {
            const aT = 0, bT = min$9(totalRads - rot, PI$13);
            const rotation = ts3dutils.M4.rotateZ(rot + bT), rotrot = ts3dutils.M4.rotateZ(rot);
            stepEndEdges = rot + bT == ts3dutils.TAU ? baseLoop : baseLoop.map(edge => edge.transform(rotation));
            const ribs = ts3dutils.arrayFromFunction(baseLoop.length, i => {
                const a = stepStartEdges[i].a, radius = a.lengthXY();
                const b = stepEndEdges[i].a;
                if (!ts3dutils.eq0(radius)) {
                    const curve = 0 == rot ? baseRibCurves[i] : baseRibCurves[i].rotateZ(rot);
                    return new PCurveEdge(curve, a, b, aT, bT, undefined, curve.tangentAt(aT), curve.tangentAt(bT), name + 'rib' + i);
                }
            });
            for (let edgeIndex = 0; edgeIndex < baseLoop.length; edgeIndex++) {
                if (baseSurfaces[edgeIndex]) {
                    const edge = stepStartEdges[edgeIndex];
                    const ipp = (edgeIndex + 1) % baseLoop.length;
                    const faceEdges = [
                        stepStartEdges[edgeIndex].flipped(),
                        !ts3dutils.eq0(edge.a.x) && ribs[edgeIndex],
                        stepEndEdges[edgeIndex],
                        !ts3dutils.eq0(edge.b.x) && ribs[ipp].flipped()
                    ].filter(x => x);
                    const surface = 0 == rot ? baseSurfaces[edgeIndex] : baseSurfaces[edgeIndex].rotateZ(rot);
                    const info = infoFactory && infoFactory.extrudeWall(edgeIndex, surface, faceEdges, undefined);
                    faces.push(Face.create(surface, faceEdges, undefined, name + 'Wall' + edgeIndex, info));
                }
            }
            stepStartEdges = stepEndEdges;
        }
        if (open) {
            const endFaceEdges = Edge.reversePath(stepEndEdges);
            const infoStart = infoFactory && infoFactory.rotationStart(basePlane, baseLoop, undefined);
            const infoEnd = infoFactory && infoFactory.rotationEnd(basePlane.flipped().rotateZ(totalRads), endFaceEdges, undefined);
            faces.push(new PlaneFace(basePlane, baseLoop, undefined, name + 'start', infoStart), new PlaneFace(basePlane.flipped().rotateZ(totalRads), endFaceEdges, undefined, name + 'end', infoEnd));
        }
        const infiniteVolume = new PlaneSurface$1(P3.ZX).edgeLoopCCW(baseLoop);
        return new B2(faces, infiniteVolume, generator);
    }
    B2T.rotateEdges = rotateEdges;
    /**
     * loop should be CCW on XZ plane for a bounded B2
     */
    //export function rotateEdgesUnsplit(loop: Edge[], rads: raddd, name: string): B2 {
    //	assert(Edge.isLoop(loop))
    //	const rotationMatrix = M4.rotateZ(rads)
    //	const open = !eq(rads, 2 * PI)
    //	const endEdges = open ? loop.map(edge => edge.transform(rotationMatrix)) : loop
    //	const edgeCount = loop.length
    //	const ribs = arrayFromFunction(edgeCount, i => {
    //		const a = loop[i].a, radius = a.lengthXY()
    //		const b = endEdges[i].a
    //		if (!eq0(radius)) {
    //			const curve = new SemiEllipseCurve(V(0, 0, a.z), V(-radius, 0, 0), V(0, -radius, 0))
    //			const aT = -PI, bT = -PI + rads
    //			return new PCurveEdge(curve, a, b, aT, bT, undefined, curve.tangentAt(aT), curve.tangentAt(bT), name
    // + 'rib' + i) } }) const faces = loop.map((edge, i) => { const ipp = (i + 1) % edgeCount console.log('ljl', i,
    // ipp, ribs) const faceEdges = [ edge.flipped(), !eq0(edge.a.x) && ribs[i], endEdges[i], !eq0(edge.b.x) &&
    // ribs[ipp].flipped()].filter(x => x) if (edge instanceof StraightEdge) { const line = edge.curve let surface if
    // (line.dir1.isParallelTo(V3.Z)) { if (eq0(edge.a.x)) { return } let flipped = edge.a.z > edge.b.z surface = new
    // SemiCylinderSurface(ribs[i].curve, !flipped ? V3.Z : V3.Z.negated()) } else if
    // (line.dir1.isPerpendicularTo(V3.Z)) { let flipped = edge.a.x > edge.b.x let surface = new PlaneSurface(new
    // P3(V3.Z, edge.a.z)) if (!flipped) surface = surface.flipped() if (!open) { const hole = flipped ? !eq0(edge.b.x)
    // && ribs[ipp].flipped() : !eq0(edge.a.x) && ribs[i] return new PlaneFace(surface, [flipped ? ribs[i] :
    // ribs[ipp].flipped()], hole && [[hole]]) } return new PlaneFace(surface, faceEdges) } else { // apex is
    // intersection of segment with Z-axis let a = edge.a, b = edge.b let apexZ = a.z - a.x * (b.z - a.z) / (b.x - a.x)
    // let apex = new V3(0, 0, apexZ) let flipped = edge.a.z > edge.b.z surface =
    // ConicSurface.atApexThroughEllipse(apex, ribs[a.x > b.x ? i : ipp].curve as SemiEllipseCurve, !flipped ? 1 : -1)
    // } return Face.create(surface, faceEdges) } if (edge.curve instanceof SemiEllipseCurve) { let flipped = undefined
    // let ell = edge.curve.rightAngled() let f1Perp = ell.f1.isPerpendicularTo(V3.Z), f2Perp =
    // ell.f2.isPerpendicularTo(V3.Z) if (L3.Z.containsPoint(ell.center) && (f1Perp || f2Perp)) { let f3length = f1Perp
    // ? ell.f1.length() : ell.f2.length() if (flipped) { f3length *= -1 } let surface = new
    // SemiEllipsoidSurface(ell.center, ell.f1, ell.f2, ell.f1.cross(ell.f2).toLength(f3length)) return new
    // RotationFace(surface, faceEdges) } } else { assert(false, edge) } }).filter(x => x) if (open) { const
    // endFaceEdges = endEdges.map(edge => edge.flipped()).reverse() faces.push( new PlaneFace(new
    // PlaneSurface(P3.ZX.flipped()), loop), new PlaneFace(new PlaneSurface(P3.ZX.rotateZ(rads)), endFaceEdges)) }
    // return new B2(faces, undefined) }
    function quaffle() {
        const baseK = B2T.sphere(1).translate(0, 1.7).flipped();
        //const baseK = B2T.box().scale(0.2).translate(0, 0.95).flipped()
        // const vs = B2T.DODECAHEDRON_VERTICES.concat(
        // B2T.DODECAHEDRON_FACE_VERTICES.map(fis => fis
        // .map(vi => B2T.DODECAHEDRON_VERTICES[vi])
        // .reduce((a,b) => a.plus(b), V3.O)
        // .unit()))
        const ss = new B2(B2T.TETRAHEDRON_VERTICES.flatMap(v => baseK.rotateAB(ts3dutils.V3.Y, v).faces), false);
        //return ss
        return B2T.sphere().and(ss);
    }
    B2T.quaffle = quaffle;
    function extrudeFace(face, dir) {
        return new B2(extrudeEdges(face.contour, face.surface.plane, dir).faces.slice(0, -2).concat(face, face.translate(dir.x, dir.y, dir.z).flipped(), face.holes.flatMap(hole => extrudeEdges(hole, face.surface.plane.flipped(), dir).faces.slice(0, -2))), false);
    }
    B2T.extrudeFace = extrudeFace;
    function loadFonts() {
        return loadFont('fonts/FiraSansMedium.woff').then(font => B2T.defaultFont = font);
    }
    B2T.loadFonts = loadFonts;
    const loadedFonts = new Map();
    function loadFont(fontPath) {
        return new Promise(function (executor, reject) {
            const font = loadedFonts.get(fontPath);
            if (font) {
                executor(font);
            }
            else {
                opentype.load(fontPath, function (err, f) {
                    if (err) {
                        reject(err);
                    }
                    else {
                        loadedFonts.set(fontPath, f);
                        executor(f);
                    }
                });
            }
        });
    }
    B2T.loadFont = loadFont;
    function loadFontsAsync(callback) {
        if (B2T.defaultFont) {
            callback();
        }
        else {
            opentype.load('fonts/FiraSansMedium.woff', function (err, font) {
                if (err) {
                    throw new Error('Could not load font: ' + err);
                }
                else {
                    B2T.defaultFont = font;
                    callback();
                }
            });
        }
    }
    B2T.loadFontsAsync = loadFontsAsync;
    function text(text, size, depth = 1, font = B2T.defaultFont) {
        const path = font.getPath(text, 0, 0, size);
        const subpaths = [];
        path.commands.forEach(c => {
            if (c.type == 'M') {
                subpaths.push([]);
            }
            subpaths.last.push(c);
        });
        const loops = subpaths.map(sp => {
            const path = new opentype.Path();
            path.commands = sp;
            const loop = Edge.reversePath(Edge.pathFromSVG(path.toPathData(13))).map(e => e.mirrorY());
            ts3dutils.assert(Edge.isLoop(loop));
            return loop;
        });
        const faces = Face.assembleFacesFromLoops(loops, new PlaneSurface$1(P3.XY), PlaneFace);
        const generator = `B2T.text(${text.sce}, ${size}, ${depth})`;
        const hello = B2.join(faces.map(face => B2T.extrudeFace(face, ts3dutils.V(0, 0, -depth))), generator);
        return hello;
    }
    B2T.text = text;
    function minorityReport() {
        const a = B2T.sphere();
        const b = B2T.text('LEO CROW', 64, 128).scale(0.1 / 32).translate(-0.5, -0.05, 1.2).flipped();
        const c = B2T.sphere(0.98);
        return a.and(b).plus(c);
    }
    B2T.minorityReport = minorityReport;
    function whatever() {
        const iso = isocahedron();
        const numbersB2 = B2.join(iso.faces.map((face, i) => {
            const numberB2 = text('' + (i + 1), 0.4, -2);
            const centroid = face.contour.map(edge => edge.a).reduce((a, b) => a.plus(b), ts3dutils.V3.O).div(3);
            const sys = ts3dutils.M4.forSys(face.contour[0].aDir, centroid.cross(face.contour[0].aDir), centroid.unit(), centroid);
            return numberB2.transform(sys.times(ts3dutils.M4.translate(-numberB2.getAABB().size().x / 2, -0.1, -0.04)));
        }));
        const s = sphere(0.9);
        //return iso.and(numbersB2)
        return iso.and(s).and(numbersB2);
        //return numbersB2
    }
    B2T.whatever = whatever;
    function d20() {
        const iso = isocahedron();
        const numbersB2 = B2.join(iso.faces.map((face, i) => {
            const numberB2 = text('' + (i + 1), 0.4, -2);
            const centroid = face.contour.map(edge => edge.a).reduce((a, b) => a.plus(b), ts3dutils.V3.O).div(3);
            const sys = ts3dutils.M4.forSys(face.contour[0].aDir, centroid.cross(face.contour[0].aDir), centroid.unit(), centroid);
            return numberB2.transform(sys.times(ts3dutils.M4.translate(-numberB2.getAABB().size().x / 2, -0.1, -0.04)));
        }));
        const s = sphere(0.9);
        //return iso.and(numbersB2)
        return iso.and(s).and(numbersB2);
        //return numbersB2
    }
    B2T.d20 = d20;
    function rotStep(edges, totalRads, count) {
        const radStep = totalRads / count;
        const open = !ts3dutils.eq(totalRads, 2 * PI$13);
        const ribCount = !open ? count : count + 1;
        const ribs = ts3dutils.arrayFromFunction(ribCount, i => {
            if (i == 0)
                return edges;
            const matrix = ts3dutils.M4.rotateZ(radStep * i);
            return edges.map(edge => edge.transform(matrix));
        });
        const horizontalEdges = ts3dutils.arrayFromFunction(count, i => {
            const ipp = (i + 1) % ribCount;
            return ts3dutils.arrayFromFunction(edges.length, j => {
                if (!ts3dutils.eq0(edges[j].a.lengthXY())) {
                    return StraightEdge.throughPoints(ribs[i][j].a, ribs[ipp][j].a);
                }
            });
        });
        const faces = [];
        let surface, face;
        edges.forEach((edge, i) => {
            const ipp = (i + 1) % edges.length;
            const projDir = ts3dutils.V3.O;
            const surface = projectCurve(ribs[r][i], projDir, ribs[r][i].deltaT() < 0);
            if (edge instanceof StraightEdge && edge.curve.dir1.isPerpendicularTo(ts3dutils.V3.Z)) {
                const flipped = edge.a.x > edge.b.x;
                surface = new PlaneSurface$1(flipped ? new P3(ts3dutils.V3.Z, edge.a.z) : new P3(ts3dutils.V3.Z.negated(), -edge.a.z));
                if (open) {
                    const newEdges = [];
                    if (!ts3dutils.eq0(edge.a.x)) {
                        newEdges.push(...ts3dutils.arrayFromFunction(count, j => horizontalEdges[j][i]));
                    }
                    newEdges.push(ribs[count][i]);
                    if (!ts3dutils.eq0(edge.b.x)) {
                        newEdges.push(...ts3dutils.arrayFromFunction(count, j => horizontalEdges[count - j - 1][ipp].flipped()));
                    }
                    newEdges.push(edge.flipped());
                    face = new PlaneFace(surface, newEdges);
                }
                else {
                    const contour = flipped
                        ? ts3dutils.arrayFromFunction(count, j => horizontalEdges[j][i])
                        : ts3dutils.arrayFromFunction(count, j => horizontalEdges[count - j - 1][ipp].flipped());
                    let hole;
                    if (flipped && !ts3dutils.eq0(edge.b.x)) {
                        hole = ts3dutils.arrayFromFunction(count, j => horizontalEdges[count - j - 1][ipp].flipped());
                    }
                    else if (!flipped && !ts3dutils.eq0(edge.a.x)) {
                        hole = ts3dutils.arrayFromFunction(count, j => horizontalEdges[j][i]);
                    }
                    face = new PlaneFace(surface, contour, hole ? [hole] : []);
                }
                faces.push(face);
                return;
            }
            else if (edge instanceof StraightEdge) {
                if (ts3dutils.eq0(edge.a.lengthXY()) && ts3dutils.eq0(edge.b.lengthXY())) {
                    return;
                }
            }
            for (let r = 0; r < count; r++) {
                const rpp = (r + 1) % ribCount;
                const faceEdges = [ribs[r][i].flipped(), horizontalEdges[r][i], ribs[rpp][i], horizontalEdges[r][ipp] && horizontalEdges[r][ipp].flipped()].filter(x => x);
                if (edge instanceof StraightEdge) {
                    const surface = new PlaneSurface$1(P3.throughPoints(faceEdges[0].a, faceEdges[1].a, faceEdges[2].a));
                    faces.push(new PlaneFace(surface, faceEdges));
                }
                else {
                    ts3dutils.assert(false, edge.toString());
                }
            }
        });
        if (open) {
            const endFaceEdges = ribs[count].map(edge => edge.flipped()).reverse();
            const endFace = new PlaneFace(new PlaneSurface$1(P3.ZX.rotateZ(totalRads)), endFaceEdges);
            faces.push(new PlaneFace(new PlaneSurface$1(P3.ZX.flipped()), edges), endFace);
        }
        return new B2(faces);
    }
    B2T.rotStep = rotStep;
    function fixEdges(edges) {
        return edges.flatMap(edge => {
            const c = edge.curve;
            if (c instanceof EllipseCurve) {
                const splitEdges = (edge.minT < 0 && edge.maxT > 0)
                    ? edge.split(0)
                    : [edge];
                return splitEdges.map(edge => {
                    if (edge.minT >= 0) {
                        return Edge.create(new SemiEllipseCurve(c.center, c.f1, c.f2, max$8(0, c.tMin), c.tMax), edge.a, edge.b, edge.aT, edge.bT, undefined, edge.aDir, edge.bDir, edge.name);
                    }
                    else {
                        // "rotate" the curve
                        return Edge.create(new SemiEllipseCurve(c.center, c.f1.negated(), c.f2.negated(), c.tMin + PI$13, min$9(PI$13, c.tMax + PI$13)), edge.a, edge.b, edge.aT + PI$13, edge.bT + PI$13, undefined, edge.aDir, edge.bDir, edge.name);
                    }
                });
            }
            if (c instanceof BezierCurve) {
                if (edge.a.like(edge.b)) {
                    return edge.split(ts3dutils.lerp(edge.aT, edge.bT, 0.5));
                }
            }
            return edge;
        });
    }
    B2T.fixEdges = fixEdges;
    function extrudeVertices(baseVertices, baseFacePlane, offset, name, generator) {
        ts3dutils.assert(baseVertices.every(v => v instanceof ts3dutils.V3), 'baseVertices.every(v => v instanceof V3)');
        ts3dutils.assertInst(P3, baseFacePlane);
        ts3dutils.assertVectors(offset);
        if (baseFacePlane.normal1.dot(offset) > 0)
            baseFacePlane = baseFacePlane.flipped();
        //if (!isCCW(baseVertices, baseFacePlane.normal1)) {
        //	baseVertices = baseVertices.reverse()
        //}
        //let topVertices = baseVertices.map((v) => v.plus(offset)).reverse()
        //let topPlane = basePlane.translated(offset)
        //let top, bottom
        //let faces = [
        //	bottom = PlaneFace.forVertices(new PlaneSurface(baseFacePlane), baseVertices),
        //	top = PlaneFace.forVertices(new PlaneSurface(baseFacePlane.flipped().translated(offset)), topVertices)]
        //let m = baseVertices.length
        //let ribs = arrayFromFunction(m, i => StraightEdge.throughPoints(baseVertices[i], topVertices[m - 1 - i]))
        //for (let i = 0; i < m; i++) {
        //	let j = (i + 1) % m
        //	faces.push(
        //		new PlaneFace(
        //			PlaneSurface.throughPoints(baseVertices[j], baseVertices[i], topVertices[m - j - 1]),
        //			[bottom.contour[i].flipped(), ribs[i], top.contour[m - j - 1].flipped(), ribs[j].flipped()], [],
        // name + 'wall' + i)) }
        const edges = StraightEdge.chain(baseVertices, true);
        generator = generator || ts3dutils.callsce('B2T.extrudeVertices', baseVertices, baseFacePlane, offset, name);
        return B2T.extrudeEdges(edges, baseFacePlane, offset, name, generator);
    }
    B2T.extrudeVertices = extrudeVertices;
    // Returns a tetrahedron (3 sided pyramid).
    // Faces will face outwards.
    // abcd can be in any order. The only constraint is that abcd cannot be on a common plane.
    function tetrahedron(a, b, c, d, name = 'tetra' + getGlobalId()) {
        ts3dutils.assertVectors(a, b, c, d);
        const dDistance = P3.throughPoints(a, b, c).distanceToPointSigned(d);
        if (ts3dutils.eq0(dDistance)) {
            throw new Error('four points are coplanar');
        }
        if (dDistance > 0) {
            [c, d] = [d, c];
        }
        const ab = StraightEdge.throughPoints(a, b);
        const ac = StraightEdge.throughPoints(a, c);
        const ad = StraightEdge.throughPoints(a, d);
        const bc = StraightEdge.throughPoints(b, c);
        const bd = StraightEdge.throughPoints(b, d);
        const cd = StraightEdge.throughPoints(c, d);
        const faces = [
            new PlaneFace(PlaneSurface$1.throughPoints(a, b, c), [ab, bc, ac.flipped()], [], name + 'abc'),
            new PlaneFace(PlaneSurface$1.throughPoints(a, d, b), [ad, bd.flipped(), ab.flipped()], [], name + 'adb'),
            new PlaneFace(PlaneSurface$1.throughPoints(b, d, c), [bd, cd.flipped(), bc.flipped()], [], name + 'bdc'),
            new PlaneFace(PlaneSurface$1.throughPoints(c, d, a), [cd, ad.flipped(), ac], [], name + 'cda'),
        ];
        const gen = `B2T.tetrahedron(${a.sce}, ${b.sce}, ${c.sce}, ${d.sce})`;
        return new B2(faces, false, gen);
    }
    B2T.tetrahedron = tetrahedron;
    const b = 1 / ts3dutils.GOLDEN_RATIO, c = 2 - ts3dutils.GOLDEN_RATIO;
    B2T.TETRAHEDRON_VERTICES = [
        new ts3dutils.V3(1, 0, -1 / Math.sqrt(2)),
        new ts3dutils.V3(-1, 0, -1 / Math.sqrt(2)),
        new ts3dutils.V3(0, -1, 1 / Math.sqrt(2)),
        new ts3dutils.V3(0, 1, 1 / Math.sqrt(2)),
    ].map(v => v.unit());
    B2T.DODECAHEDRON_VERTICES = [
        new ts3dutils.V3(c, 0, 1),
        new ts3dutils.V3(-c, 0, 1),
        new ts3dutils.V3(-b, b, b),
        new ts3dutils.V3(0, 1, c),
        new ts3dutils.V3(b, b, b),
        new ts3dutils.V3(b, -b, b),
        new ts3dutils.V3(0, -1, c),
        new ts3dutils.V3(-b, -b, b),
        new ts3dutils.V3(c, 0, -1),
        new ts3dutils.V3(-c, 0, -1),
        new ts3dutils.V3(-b, -b, -b),
        new ts3dutils.V3(0, -1, -c),
        new ts3dutils.V3(b, -b, -b),
        new ts3dutils.V3(b, b, -b),
        new ts3dutils.V3(0, 1, -c),
        new ts3dutils.V3(-b, b, -b),
        new ts3dutils.V3(1, c, 0),
        new ts3dutils.V3(-1, c, 0),
        new ts3dutils.V3(-1, -c, 0),
        new ts3dutils.V3(1, -c, 0),
    ].map(v => v.unit());
    B2T.DODECAHEDRON_FACE_VERTICES = [
        [4, 3, 2, 1, 0],
        [7, 6, 5, 0, 1],
        [12, 11, 10, 9, 8],
        [15, 14, 13, 8, 9],
        [14, 3, 4, 16, 13],
        [3, 14, 15, 17, 2],
        [11, 6, 7, 18, 10],
        [6, 11, 12, 19, 5],
        [4, 0, 5, 19, 16],
        [12, 8, 13, 16, 19],
        [15, 9, 10, 18, 17],
        [7, 1, 2, 17, 18]
    ];
    B2T.OCTAHEDRON_VERTICES = [
        new ts3dutils.V3(1, 0, 0),
        new ts3dutils.V3(-1, 0, 0),
        new ts3dutils.V3(0, 1, 0),
        new ts3dutils.V3(0, -1, 0),
        new ts3dutils.V3(0, 0, 1),
        new ts3dutils.V3(0, 0, -1)
    ];
    B2T.OCTAHEDRON_FACE_VERTICES = [
        [0, 2, 4],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],
        [2, 0, 5],
        [1, 2, 5],
        [3, 1, 5],
        [0, 3, 5]
    ];
    const { x: s, y: t } = new ts3dutils.V3(1, ts3dutils.GOLDEN_RATIO, 0).unit();
    B2T.ISOCAHEDRON_VERTICES = [
        new ts3dutils.V3(-s, t, 0),
        new ts3dutils.V3(s, t, 0),
        new ts3dutils.V3(-s, -t, 0),
        new ts3dutils.V3(s, -t, 0),
        new ts3dutils.V3(0, -s, t),
        new ts3dutils.V3(0, s, t),
        new ts3dutils.V3(0, -s, -t),
        new ts3dutils.V3(0, s, -t),
        new ts3dutils.V3(t, 0, -s),
        new ts3dutils.V3(t, 0, s),
        new ts3dutils.V3(-t, 0, -s),
        new ts3dutils.V3(-t, 0, s)
    ];
    B2T.ISOCAHEDRON_FACE_VERTICES = [
        // 5 faces around point 0
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        // 5 adjacent faces
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        // 5 faces around point 3
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        // 5 adjacent faces
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1]
    ];
    function dodecahedron() {
        return makePlatonic(B2T.DODECAHEDRON_VERTICES, B2T.DODECAHEDRON_FACE_VERTICES, 'B2T.dodecahedron()');
    }
    B2T.dodecahedron = dodecahedron;
    function octahedron() {
        return makePlatonic(B2T.OCTAHEDRON_VERTICES, B2T.OCTAHEDRON_FACE_VERTICES, 'B2T.octahedron()');
    }
    B2T.octahedron = octahedron;
    function isocahedron() {
        return makePlatonic(B2T.ISOCAHEDRON_VERTICES, B2T.ISOCAHEDRON_FACE_VERTICES, 'B2T.octahedron()');
    }
    B2T.isocahedron = isocahedron;
    function makePlatonic(VS, FVIS, generator) {
        const edgeMap = new Map();
        const faces = FVIS.map(faceIndexes => {
            const surface = PlaneSurface$1.throughPoints(VS[faceIndexes[0]], VS[faceIndexes[1]], VS[faceIndexes[2]]);
            const contour = ts3dutils.arrayFromFunction(faceIndexes.length, i => {
                const ipp = (i + 1) % faceIndexes.length;
                const iA = faceIndexes[i], iB = faceIndexes[ipp];
                const iMin = min$9(iA, iB), iMax = max$8(iA, iB), edgeID = iMin * VS.length + iMax;
                let edge = edgeMap.get(edgeID);
                !edge && edgeMap.set(edgeID, edge = StraightEdge.throughPoints(VS[iMin], VS[iMax]));
                return iA < iB ? edge : edge.flipped();
            });
            return new PlaneFace(surface, contour);
        });
        return new B2(faces, false, generator);
    }
    function pyramidEdges(baseEdges, apex, name = 'pyramid' + getGlobalId()) {
        ts3dutils.assertInst(Edge, ...baseEdges);
        ts3dutils.assertVectors(apex);
        const ribs = baseEdges.map(baseEdge => StraightEdge.throughPoints(apex, baseEdge.a));
        const faces = baseEdges.map((baseEdge, i) => {
            const faceName = name + 'Wall' + i;
            const ipp = (i + 1) % baseEdges.length;
            const faceEdges = [ribs[i], baseEdge, ribs[ipp].flipped()];
            const surface = undefined; // TODO
            return Face.create(surface, faceEdges, undefined, faceName);
        });
        const bottomFace = Face.create(baseSurface, baseEdges);
        faces.push(bottomFace);
        const generator = ts3dutils.callsce('B2T.pyramidEdges', baseEdges, apex, name);
        return new B2(faces, false, generator, name);
    }
    B2T.pyramidEdges = pyramidEdges;
})(exports.B2T || (exports.B2T = {}));

class CustomPlane extends P3 {
    constructor(anchor, right, up, name, color = ts3dutils.randomColor(), rightStart = -500, rightEnd = 500, upStart = -500, upEnd = 500) {
        const { normal1, w } = P3.forAnchorAndPlaneVectors(anchor, right, up);
        super(normal1, w);
        this.up = up;
        this.right = right;
        this.sMin = rightStart;
        this.sMax = rightEnd;
        this.tMin = upStart;
        this.tMax = upEnd;
        this.color = color;
        this.name = name;
    }
    get plane() { return this; }
    toPlaneSurface() {
        return new PlaneSurface$1(this, this.right, this.up);
    }
    static forPlane(plane, color, name) {
        //assert(!name)
        const up = plane.normal1.getPerpendicular().unit(), right = up.cross(plane.normal1);
        return new CustomPlane(plane.anchor, right, up, name, color);
    }
    static fromPlaneSurface(surface) {
        return new CustomPlane(surface.plane.anchor, surface.right, surface.up, 'genCustomPlane' + getGlobalId());
    }
    distanceTo(line, mindist) {
        return [
            new L3$1(this.anchor.plus(this.right.times(this.sMin)), this.up),
            new L3$1(this.anchor.plus(this.right.times(this.sMax)), this.up),
            new L3$1(this.anchor.plus(this.up.times(this.tMin)), this.right),
            new L3$1(this.anchor.plus(this.up.times(this.tMax)), this.right)
        ].map((line2, line2Index) => {
            const info = line2.infoClosestToLine(line);
            if ((isNaN(info.t) // parallel LINES
                || line2Index < 2 && this.tMin <= info.t && info.t <= this.tMax
                || line2Index >= 2 && this.sMin <= info.t && info.t <= this.sMax)
                && info.distance <= mindist) {
                return info.s;
            }
            else {
                return Infinity;
            }
        }).min();
    }
    distanceTo2(line, mindist) {
        return [
            new L3$1(this.anchor.plus(this.right.times(this.sMin)), this.up),
            new L3$1(this.anchor.plus(this.right.times(this.sMax)), this.up),
            new L3$1(this.anchor.plus(this.up.times(this.tMin)), this.right),
            new L3$1(this.anchor.plus(this.up.times(this.tMax)), this.right)
        ].map((line2, line2Index) => {
            const info = line2.infoClosestToLine(line);
            if ((isNaN(info.t) // parallel LINES
                || line2Index < 2 && this.tMin <= info.t && info.t <= this.tMax
                || line2Index >= 2 && this.sMin <= info.t && info.t <= this.sMax)
                && info.distance <= mindist) {
                return info.distance;
            }
            else {
                return Infinity;
            }
        }).min();
    }
}

const { PI: PI$14, sign: sign$8, ceil: ceil$11, floor: floor$9, abs: abs$11 } = Math;
class Edge extends ts3dutils.Transformable {
    constructor(curve, a, b, aT, bT, flippedOf, name) {
        super();
        this.curve = curve;
        this.a = a;
        this.b = b;
        this.aT = aT;
        this.bT = bT;
        this.flippedOf = flippedOf;
        this.name = name;
        ts3dutils.assertNumbers(aT, bT);
        ts3dutils.assert(!ts3dutils.eq(aT, bT));
        ts3dutils.assertVectors(a, b);
        ts3dutils.assertf(() => curve instanceof Curve, curve);
        ts3dutils.assertf(() => !curve.isValidT || curve.isValidT(aT) && curve.isValidT(bT), aT + ' ' + bT);
        ts3dutils.assertf(() => curve.at(aT).like(a), +a);
        ts3dutils.assertf(() => curve.at(bT).like(b), '' + curve.at(bT) + b);
        ts3dutils.assertf(() => ts3dutils.fuzzyBetween(aT, curve.tMin, curve.tMax));
        ts3dutils.assertf(() => ts3dutils.fuzzyBetween(bT, curve.tMin, curve.tMax));
        this.aT = ts3dutils.clamp(aT, curve.tMin, curve.tMax);
        this.bT = ts3dutils.clamp(bT, curve.tMin, curve.tMax);
        this.reversed = this.aT > this.bT;
    }
    get minT() { return Math.min(this.aT, this.bT); }
    get maxT() { return Math.max(this.aT, this.bT); }
    static forCurveAndTs(curve, aT = curve.tMin, bT = curve.tMax) {
        return Edge.create(curve, curve.at(aT), curve.at(bT), aT, bT, undefined, aT < bT ? curve.tangentAt(aT) : curve.tangentAt(aT).negated(), aT < bT ? curve.tangentAt(bT) : curve.tangentAt(bT).negated());
    }
    static create(curve, a, b, aT, bT, flippedOf, aDir, bDir, name) {
        if (curve instanceof L3$1) {
            return new StraightEdge(curve, a, b, aT, bT, flippedOf, name);
        }
        else {
            return new PCurveEdge(curve, a, b, aT, bT, flippedOf, aDir, bDir, name);
        }
    }
    static isLoop(loop) {
        return loop.every((edge, i) => edge.b.like(loop[(i + 1) % loop.length].a));
    }
    static edgesIntersect(e1, e2) {
        // TODO: still getting some NaNs here..
        ts3dutils.assertNumbers(e1.curve.hlol, e2.curve.hlol);
        ts3dutils.assertInst(Edge, e1, e2);
        if (e1.curve.hlol < e2.curve.hlol) {
            [e2, e1] = [e1, e2];
        }
        const sts = e1.curve.isInfosWithCurve(e2.curve);
        if (sts.some(info => isNaN(info.tThis) || isNaN(info.tOther))) {
            console.log(e1.sce);
            console.log(e2.sce);
            ts3dutils.assert(false);
        }
        return sts.some(
        /// (  e1.aT < tThis < e1.bT  )  &&  (  e2.aT < tOther < e2.bT  )
        ({ tThis, tOther }) => {
            return e1.tValueInside(tThis) && e2.tValueInside(tOther);
        });
    }
    static assertLoop(edges) {
        edges.forEach((edge, i) => {
            const j = (i + 1) % edges.length;
            ts3dutils.assert(edge.b.like(edges[j].a), `edges[${i}].b != edges[${j}].a (${edges[i].b.sce} != ${edges[j].a.sce})`);
        });
    }
    static ngon(n = 3, radius = 1) {
        return StraightEdge.chain(ts3dutils.arrayFromFunction(n, i => ts3dutils.V3.polar(radius, ts3dutils.TAU * i / n)));
    }
    static star(pointCount = 5, r0 = 1, r1 = 0.5) {
        const vertices = ts3dutils.arrayFromFunction(pointCount * 2, i => ts3dutils.V3.polar(0 == i % 2
            ? r0
            : r1, ts3dutils.TAU * i / pointCount / 2));
        return StraightEdge.chain(vertices);
    }
    static reversePath(path, doReverse = true) {
        return doReverse ? ts3dutils.arrayFromFunction(path.length, i => path[path.length - 1 - i].flipped()) : path;
    }
    static rect(width = 1, height = width) {
        const vertices = [new ts3dutils.V3(0, 0, 0), new ts3dutils.V3(width, 0, 0), new ts3dutils.V3(width, height, 0), new ts3dutils.V3(0, height, 0)];
        return StraightEdge.chain(vertices);
    }
    static reuleaux(n = 3, radius = 1) {
        ts3dutils.assert(3 <= n);
        ts3dutils.assert(1 == n % 2);
        const corners = ts3dutils.arrayFromFunction(n, i => ts3dutils.V3.polar(radius, ts3dutils.TAU * i / n));
        return ts3dutils.arrayFromFunction(n, i => {
            const aI = (i + floor$9(n / 2)) % n, bI = (i + ceil$11(n / 2)) % n;
            const a = corners[aI], b = corners[bI];
            const center = corners[i];
            const f1 = center.to(a), curve = new SemiEllipseCurve(center, f1, ts3dutils.V3.Z.cross(f1));
            return Edge.create(curve, a, b, 0, curve.pointT(b), undefined, ts3dutils.V3.Z.cross(f1), ts3dutils.V3.Z.cross(center.to(b)));
        });
    }
    static round(edges, radius) {
        if (ts3dutils.eq0(radius)) {
            return edges;
        }
        const corners = edges.map((edge, i) => {
            const j = (i + 1) % edges.length, nextEdge = edges[j];
            if (!edge.b.like(nextEdge.a))
                return;
            const angleToNext = edge.bDir.angleTo(nextEdge.aDir);
            const c1 = edge.curve, c2 = nextEdge.curve;
            if (c1 instanceof L3$1 && c2 instanceof L3$1) {
                const normal = c1.dir1.cross(c2.dir1);
                if (ts3dutils.eq0(angleToNext))
                    return;
                const l1inside = normal.cross(c1.dir1), l2inside = normal.cross(c2.dir1);
                const l1offset = c1.transform(ts3dutils.M4.translate(l1inside.toLength(radius)));
                const l2offset = c2.transform(ts3dutils.M4.translate(l2inside.toLength(radius)));
                const center = l1offset.isInfoWithLine(l2offset);
                if (!center)
                    throw new Error('tangential curves');
                const cornerA = center.plus(l1inside.toLength(-radius));
                const cornerB = center.plus(l2inside.toLength(-radius));
                const f1 = l1inside.toLength(-radius);
                const curve = new SemiEllipseCurve(center, f1, normal.cross(f1).toLength(radius));
                const cornerEdge = Edge.create(curve, cornerA, cornerB, 0, curve.pointT(cornerB), undefined, c1.dir1, c2.dir1);
                return cornerEdge;
            }
            else {
                return Edge.arbitraryCorner(edge, nextEdge, radius);
            }
        });
        const result = edges.flatMap((edge, i) => {
            const h = (i + edges.length - 1) % edges.length, j = (i + 1) % edges.length;
            const prevCorner = corners[h], nextCorner = corners[i];
            if (!prevCorner && !nextCorner) {
                return edge;
            }
            const [aT, a, aDir] = !prevCorner
                ? [edge.aT, edge.a, edge.aDir]
                : [edge.curve.pointT(prevCorner.b), prevCorner.b, prevCorner.bDir];
            const [bT, b, bDir] = !nextCorner
                ? [edge.bT, edge.b, edge.bDir]
                : [edge.curve.pointT(nextCorner.a), nextCorner.a, nextCorner.aDir];
            const newEdge = Edge.create(edge.curve, a, b, aT, bT, undefined, aDir, bDir);
            return !nextCorner ? newEdge : [newEdge, nextCorner];
        });
        return result;
    }
    static arbitraryCorner(e1, e2, radius) {
        const c1 = e1.curve, c2 = e2.curve;
        function f([t1, t2]) {
            const p1 = c1.at(t1), p2 = c2.at(t2);
            const dp1 = c1.tangentAt(t1), dp2 = c2.tangentAt(t2);
            const virtualPlaneNormal = dp1.cross(dp2);
            const normal1 = virtualPlaneNormal.cross(dp1).unit(), normal2 = virtualPlaneNormal.cross(dp2).unit();
            const dirCross = normal1.cross(normal2);
            if (virtualPlaneNormal.likeO()) {
                ts3dutils.assert(false);
            } // lines parallel
            const p1p2 = p1.to(p2);
            // check if distance is zero (see also L3.distanceToLine)
            if (!ts3dutils.eq0(p1p2.dot(virtualPlaneNormal))) {
                ts3dutils.assert(false);
            }
            const l1 = new L3$1(p1, normal1), l2 = new L3$1(p2, normal2);
            const uh = l1.infoClosestToLine(l2), uh2 = l1.isInfoWithLine(l2);
            const dist1 = p1p2.cross(normal2).dot(dirCross) / dirCross.squared();
            const dist2 = p1p2.cross(normal1).dot(dirCross) / dirCross.squared();
            const g1 = p1.plus(normal1.times(dist1));
            const g2 = p2.plus(normal2.times(dist2));
            ts3dutils.assert(g1.like(g2));
            return [abs$11(dist1) - radius, abs$11(dist2) - radius];
        }
        const startT1 = e1.bT - radius * sign$8(e1.deltaT()) / e1.bDir.length();
        const startT2 = e2.aT + radius * sign$8(e2.deltaT()) / e2.aDir.length();
        const [t1, t2] = ts3dutils.newtonIterate(f, [startT1, startT2]);
        const cornerA = e1.curve.at(t1);
        const cornerB = e2.curve.at(t2);
        const p1 = c1.at(t1), p2 = c2.at(t2);
        const dp1 = c1.tangentAt(t1), dp2 = c2.tangentAt(t2);
        const virtualPlaneNormal = dp1.cross(dp2);
        const normal1 = virtualPlaneNormal.cross(dp1).unit(), normal2 = virtualPlaneNormal.cross(dp2).unit();
        const f1 = normal1.toLength(-radius);
        const center = cornerA.minus(f1);
        const curve = new SemiEllipseCurve(center, f1, virtualPlaneNormal.cross(f1).toLength(radius));
        const cornerEdge = Edge.create(curve, cornerA, cornerB, 0, curve.pointT(cornerB), undefined, c1.tangentAt(t1), c2.tangentAt(t2));
        return cornerEdge;
    }
    static pathFromSVG(pathString) {
        let currentPos = undefined;
        const parsed = new svgPathdata.SVGPathData(pathString).toAbs().normalizeHVZ().sanitize(ts3dutils.NLA_PRECISION).annotateArcs().commands;
        const path = [];
        for (const c of parsed) {
            ts3dutils.assert('x' in c && 'y' in c);
            const endPos = new ts3dutils.V3(c.x, c.y, 0);
            switch (c.type) {
                case svgPathdata.SVGPathData.LINE_TO:
                    path.push(StraightEdge.throughPoints(currentPos, endPos));
                    break;
                case svgPathdata.SVGPathData.CURVE_TO: {
                    const c1 = new ts3dutils.V3(c.x1, c.y1, 0);
                    const c2 = new ts3dutils.V3(c.x2, c.y2, 0);
                    const curve = new BezierCurve(currentPos, c1, c2, endPos, 0, 1);
                    const edge = new PCurveEdge(curve, currentPos, endPos, 0, 1, undefined, curve.tangentAt(0), curve.tangentAt(1));
                    path.push(edge);
                    break;
                }
                case svgPathdata.SVGPathData.QUAD_TO: {
                    const c1 = new ts3dutils.V3(c.x1, c.y1, 0);
                    const curve = ParabolaCurve.quadratic(currentPos, c1, endPos).rightAngled();
                    const edge = new PCurveEdge(curve, currentPos, endPos, curve.tMin, curve.tMax, undefined, curve.tangentAt(curve.tMin), curve.tangentAt(curve.tMax));
                    path.push(edge);
                    break;
                }
                case svgPathdata.SVGPathData.ARC: {
                    const phi1 = c.phi1 * ts3dutils.DEG, phi2 = c.phi2 * ts3dutils.DEG, [phiMin, phiMax] = [phi1, phi2].sort(ts3dutils.MINUS);
                    const stops = ts3dutils.arrayRange(-3, 4, 1).map(n => n * PI$14).filter(stop => phiMin <= stop && stop <= phiMax);
                    const center = ts3dutils.V(c.cX, c.cY);
                    const f1 = ts3dutils.V3.polar(c.rX, c.xRot * ts3dutils.DEG);
                    const f2 = ts3dutils.V3.polar(c.rY, c.xRot * ts3dutils.DEG + Math.PI / 2);
                    const edges = ts3dutils.getIntervals(stops, phiMin, phiMax).map(([t1, t2]) => {
                        const deltaT = t2 - t1;
                        const t1_ = ts3dutils.mod(t1, ts3dutils.TAU);
                        const t2_ = t1_ + deltaT;
                        ts3dutils.assert(t1_ >= 0 == t2_ >= 0);
                        const gtPI = t1_ > PI$14 || t2_ > PI$14;
                        const aT = gtPI ? t1_ - PI$14 : t1_;
                        const bT = gtPI ? t2_ - PI$14 : t2_;
                        const curve = new SemiEllipseCurve(center, gtPI ? f1.negated() : f1, gtPI ? f2.negated() : f2);
                        const a = phi1 == t1 ? currentPos : phi2 == t1 ? endPos : curve.at(aT);
                        const b = phi1 == t2 ? currentPos : phi2 == t2 ? endPos : curve.at(bT);
                        return new PCurveEdge(curve, a, b, aT, bT, undefined, curve.tangentAt(aT), curve.tangentAt(bT));
                    });
                    path.push(...c.phiDelta > 0 ? edges : Edge.reversePath(edges));
                    break;
                }
            }
            currentPos = endPos;
        }
        return path;
    }
    toString() {
        return ts3dutils.callsce('new ' + this.constructor.name, this.curve, this.a, this.b, this.aT, this.bT, undefined, this.aDir, this.bDir);
    }
    split(t) {
        const p = this.curve.at(t);
        const pDir = this.tangentAt(t);
        return [
            Edge.create(this.curve, this.a, p, this.aT, t, undefined, this.aDir, pDir, this.name + 'left'),
            Edge.create(this.curve, p, this.b, t, this.bT, undefined, pDir, this.bDir, this.name + 'left')
        ];
    }
    colinearToLine(line) {
        return this.curve instanceof L3$1 && this.curve.isColinearTo(line);
    }
    tValueInside(t) {
        return this.aT < this.bT
            ? ts3dutils.lt(this.aT, t) && ts3dutils.lt(t, this.bT)
            : ts3dutils.lt(this.bT, t) && ts3dutils.lt(t, this.aT);
    }
    isValidT(t) {
        return this.aT < this.bT
            ? ts3dutils.le(this.aT, t) && ts3dutils.le(t, this.bT)
            : ts3dutils.le(this.bT, t) && ts3dutils.le(t, this.aT);
    }
    clampedT(t) {
        return this.aT < this.bT
            ? ts3dutils.clamp(t, this.aT, this.bT)
            : ts3dutils.clamp(t, this.bT, this.aT);
    }
    /**
     * this is equals-equals. "isColinearTo" might make more sense but can't be used, because you can't get a
     * consistent hashCode for colinear curves
     * @param obj
     * @returns
     */
    equals(obj) {
        return this === obj ||
            this.constructor == obj.constructor
                && this.a.equals(obj.a)
                && this.b.equals(obj.b)
                && this.curve.equals(obj.curve);
    }
    hashCode() {
        let hashCode = 0;
        hashCode = hashCode * 31 + this.a.hashCode();
        hashCode = hashCode * 31 + this.b.hashCode();
        hashCode = hashCode * 31 + this.curve.hashCode();
        return hashCode | 0;
    }
    like(edge) {
        // TODO this breaks on colinear edges,
        // TODO: what, where?
        return this === edge ||
            edge instanceof Edge &&
                this.curve.isColinearTo(edge.curve)
                && this.a.like(edge.a)
                && this.b.like(edge.b);
    }
    isCanon() {
        return !this.reversed;
    }
    getCanon() {
        return this.reversed
            ? this.flipped()
            : this;
    }
    overlaps(edge, noback) {
        ts3dutils.assert(this.curve.isColinearTo(edge.curve));
        const edgeAT = this.curve.containsPoint(edge.a) && this.curve.pointT(edge.a);
        const edgeBT = this.curve.containsPoint(edge.b) && this.curve.pointT(edge.b);
        if (false === edgeAT && false === edgeBT) {
            return noback ? false : edge.overlaps(this, true);
        }
        const flipped = false !== edgeAT ? this.tangentAt(edgeAT).dot(edge.aDir) : this.tangentAt(edgeBT).dot(edge.bDir);
        return !(ts3dutils.le(edgeMaxT, this.minT) || ts3dutils.le(this.maxT, edgeMinT));
    }
    getAABB() {
        const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
        this.curve.roots().forEach((ts, dim) => {
            ts.forEach(t => {
                if (ts3dutils.lt(this.minT, t) && ts3dutils.lt(t, this.maxT)) {
                    min[dim] = Math.min(min[dim], this.curve.at(t).e(dim));
                    max[dim] = Math.max(max[dim], this.curve.at(t).e(dim));
                }
            });
        });
        const aabb = new ts3dutils.AABB(ts3dutils.V(min), ts3dutils.V(max));
        aabb.addPoint(this.a);
        aabb.addPoint(this.b);
        return aabb;
    }
    length(steps = 1) {
        return this.curve.arcLength(this.minT, this.maxT, steps);
    }
    deltaT() {
        return this.bT - this.aT;
    }
    atAvgT() {
        return this.curve.at((this.minT + this.maxT) / 2);
    }
}
class PCurveEdge extends Edge {
    constructor(curve, a, b, aT, bT, flippedOf, aDir, bDir, name) {
        super(curve, a, b, aT, bT, flippedOf, name);
        this.flippedOf = flippedOf;
        this.aDir = aDir;
        this.bDir = bDir;
        ts3dutils.assertVectors(aDir, bDir);
        ts3dutils.assertf(() => !aDir.likeO(), curve);
        ts3dutils.assertf(() => !bDir.likeO(), curve);
        if (!(curve instanceof PICurve$1)) {
            // TODO
            ts3dutils.assertf(() => curve.tangentAt(aT).likeOrReversed(aDir), '' + aT + curve.tangentAt(aT).sce + ' ' + aDir.sce);
            ts3dutils.assertf(() => curve.tangentAt(bT).likeOrReversed(bDir));
        }
        ts3dutils.assert(this.reversed === this.aDir.dot(curve.tangentAt(aT)) < 0, aT + ' ' + bT + ' ' + curve.constructor.name + ' ' + this.aDir.sce + ' ' + this.bDir.sce + ' ' + curve.tangentAt(aT));
        ts3dutils.assert(this.reversed === this.bDir.dot(curve.tangentAt(bT)) < 0, aT + ' ' + bT + ' ' + curve.constructor.name + ' ' + this.aDir.sce + ' ' + this.bDir.sce + ' ' + curve.tangentAt(aT));
    }
    static forCurveAndTs(curve, aT, bT, name) {
        return new PCurveEdge(curve, curve.at(aT), curve.at(bT), aT, bT, undefined, aT < bT ? curve.tangentAt(aT) : curve.tangentAt(aT).negated(), aT < bT ? curve.tangentAt(bT) : curve.tangentAt(bT).negated(), name);
    }
    toSource() {
        return ts3dutils.callsce('new PCurveEdge', this.curve, this.a, this.b, this.aT, this.bT, undefined, this.aDir, this.bDir, this.name);
    }
    getVerticesNo0() {
        return this.curve.calcSegmentPoints(this.aT, this.bT, this.a, this.b, this.reversed, false);
    }
    pointsCount() {
        return this.points().length;
    }
    points() {
        return this.curve.calcSegmentPoints(this.aT, this.bT, this.a, this.b, this.reversed, true);
    }
    rotViaPlane(normal, reversed) {
        let rot = this.aDir.angleRelativeNormal(this.bDir, normal);
        const counterClockWise = (normal.dot(this.curve.normal) > 0) === !this.reversed;
        if (counterClockWise) {
            // counterclockwise rotation, i.e. rot > 0
            if (rot < 0)
                rot += 2 * Math.PI;
        }
        else {
            if (rot > 0)
                rot -= 2 * Math.PI;
        }
        return rot;
    }
    edgeISTsWithSurface(surface) {
        return this.curve.isTsWithSurface(surface)
            .map(edgeT => ts3dutils.snap2(edgeT, this.aT, this.bT))
            .filter(edgeT => this.minT <= edgeT && edgeT <= this.maxT);
    }
    edgeISTsWithPlane(surface) {
        return this.curve.isTsWithPlane(surface)
            .map(edgeT => ts3dutils.snap2(edgeT, this.aT, this.bT))
            .filter(edgeT => this.minT <= edgeT && edgeT <= this.maxT);
    }
    tangentAt(t) {
        return !this.reversed ? this.curve.tangentAt(t) : this.curve.tangentAt(t).negated();
    }
    flipped() {
        return this.flippedOf || (this.flippedOf = new PCurveEdge(this.curve, this.b, this.a, this.bT, this.aT, this, this.bDir.negated(), this.aDir.negated(), this.name));
    }
    transform(m4, desc) {
        return new PCurveEdge(this.curve.transform(m4), m4.transformPoint(this.a), m4.transformPoint(this.b), this.aT, this.bT, undefined, m4.transformVector(this.aDir), m4.transformVector(this.bDir), '' + this.name + desc);
    }
    isCoEdge(edge) {
        return this === edge || this === edge.flippedOf ||
            this.curve.isColinearTo(edge.curve) && (this.a.like(edge.a) && this.b.like(edge.b)
                || this.a.like(edge.b) && this.b.like(edge.a));
    }
}
class StraightEdge extends Edge {
    // flippedOf: StraightEdge
    constructor(line, a, b, aT, bT, flippedOf, name) {
        super(line, a, b, aT, bT, flippedOf, name);
        this.flippedOf = flippedOf;
        ts3dutils.assertInst(L3$1, line);
        !flippedOf || ts3dutils.assertInst(StraightEdge, flippedOf);
        !name || ts3dutils.assertf(() => 'string' === typeof name, name);
        ts3dutils.assert(!a.like(b), '!a.like(b)' + a + b); // don't put in super as it will break full ellipse
        this.tangent = this.aT < this.bT ? this.curve.dir1 : this.curve.dir1.negated();
    }
    get aDir() {
        return this.tangent;
    }
    get bDir() {
        return this.tangent;
    }
    static throughPoints(a, b, name) {
        return new StraightEdge(L3$1.throughPoints(a, b, 0, a.to(b).length()), a, b, 0, a.to(b).length(), undefined, name);
    }
    /**
     * Create a list of StraightEdges from a list of vertices.
     * @param vertices
     * @param closed Whether to connect the first and last vertices. Defaults to true.
     * @returns
     */
    static chain(vertices, closed = true) {
        const vc = vertices.length;
        return ts3dutils.arrayFromFunction(closed ? vc : vc - 1, i => StraightEdge.throughPoints(vertices[i], vertices[(i + 1) % vc]));
    }
    toSource() {
        return ts3dutils.callsce('new StraightEdge', this.curve, this.a, this.b, this.aT, this.bT);
    }
    getVerticesNo0() {
        return [this.b];
    }
    pointsCount() {
        return 2;
    }
    points() {
        return [this.a, this.b];
    }
    edgeISTsWithPlane(plane) {
        const edgeT = ts3dutils.snap2(this.curve.isTWithPlane(plane), this.aT, this.bT);
        return (this.minT <= edgeT && edgeT <= this.maxT) ? [edgeT] : [];
    }
    edgeISTsWithSurface(surface) {
        if (surface instanceof PlaneSurface$1) {
            return this.edgeISTsWithPlane(surface.plane);
        }
        else {
            return surface.isTsForLine(this.curve)
                .map(edgeT => ts3dutils.snap2(edgeT, this.aT, this.bT))
                .filter(edgeT => this.minT <= edgeT && edgeT <= this.maxT);
        }
    }
    tangentAt() {
        return this.tangent;
    }
    flipped() {
        return this.flippedOf || (this.flippedOf = new StraightEdge(this.curve, this.b, this.a, this.bT, this.aT, this, this.name));
    }
    transform(m4, desc) {
        const lineDir1TransLength = m4.transformVector(this.curve.dir1).length();
        return new StraightEdge(this.curve.transform(m4), m4.transformPoint(this.a), m4.transformPoint(this.b), this.aT * lineDir1TransLength, this.bT * lineDir1TransLength, undefined, '' + this.name + desc);
    }
    isCoEdge(edge) {
        return this === edge || this === edge.flippedOf || edge.constructor === StraightEdge && (this.a.like(edge.a) && this.b.like(edge.b)
            || this.a.like(edge.b) && this.b.like(edge.a));
    }
    getEdgeT(p) {
        ts3dutils.assertVectors(p);
        let edgeT = p.minus(this.curve.anchor).dot(this.curve.dir1);
        if (!ts3dutils.eq0(this.curve.at(edgeT).distanceTo(p))) {
            return;
        }
        edgeT = ts3dutils.snap2(edgeT, this.aT, this.bT);
        return (this.minT <= edgeT && edgeT <= this.maxT) ? edgeT : undefined;
    }
}

/**
 * Created by aval on 19.04.2017.
 */
class FaceInfoFactory {
    static makeStatic(staticInfo) {
        return new class extends FaceInfoFactory {
            constructor() {
                super();
            }
            info(surface, contour, holes) {
                return staticInfo;
            }
        };
    }
    info(surface, contour, holes) {
        throw new Error('no default implementation');
    }
    extrudeBottom(surface, contour, holes = []) {
        return this.info(surface, contour, holes);
    }
    extrudeTop(surface, contour, holes = []) {
        return this.info(surface, contour, holes);
    }
    extrudeWall(index, surface, contour, holes = []) {
        return this.info(surface, contour, holes);
    }
    rotationWall(index, surface, contour, holes = []) {
        return this.info(surface, contour, holes);
    }
    rotationStart(surface, contour, holes = []) {
        return this.info(surface, contour, holes);
    }
    rotationEnd(surface, contour, holes = []) {
        return this.info(surface, contour, holes);
    }
    newSubFace(original, surface, contour, holes = []) {
        return original.info;
    }
    transform(original, m4, desc, surface, contour, holes = []) {
        return original.info;
    }
}

const { PI: PI$15, min: min$10, max: max$9, sign: sign$9, ceil: ceil$12, floor: floor$10, abs: abs$12 } = Math;
class Face extends ts3dutils.Transformable {
    constructor(surface, contour, holes = [], name, info) {
        super();
        this.surface = surface;
        this.contour = contour;
        this.holes = holes;
        this.name = name;
        this.info = info;
        //assert(name)
        Edge.assertLoop(contour);
        ts3dutils.assert(contour.every(f => f instanceof Edge), () => 'contour.every(f => f instanceof Edge)' + contour);
        // contour.forEach(e => !surface.containsCurve(e.curve) &&
        // console.log('FAIL:'+surface.distanceToPoint(e.curve.anchor)))
        contour.forEach(e => {
            ts3dutils.assert(surface.containsCurve(e.curve), 'edge not in surface ' + e + surface);
        });
        ts3dutils.assert(surface.edgeLoopCCW(contour), surface.toString() + contour.join('\n'));
        holes && holes.forEach(hole => Edge.assertLoop(hole));
        holes && holes.forEach(hole => ts3dutils.assert(!surface.edgeLoopCCW(hole)));
        ts3dutils.assert(!holes || holes.constructor == Array, holes && holes.toString());
        this.allEdges = Array.prototype.concat.apply(this.contour, this.holes);
    }
    static assembleFacesFromLoops(loops, surface, faceConstructor) {
        function placeRecursively(newLoopInfo, loopInfos) {
            if (loopInfos.length == 0) {
                loopInfos.push(newLoopInfo);
            }
            else {
                const subLoopInfo = loopInfos.find(loopInfo => B2.loop1ContainsLoop2(loopInfo.loop, loopInfo.ccw, newLoopInfo.loop, newLoopInfo.ccw, surface));
                if (subLoopInfo) {
                    placeRecursively(newLoopInfo, subLoopInfo.subloops);
                }
                else {
                    // newLoopInfo isnt contained by any other subLoopInfo
                    for (let i = loopInfos.length; --i >= 0;) {
                        const subLoopInfo = loopInfos[i];
                        //console.log('cheving subLoopInfo', surface.loopContainsPoint(newLoopInfo.edges,
                        // subLoopInfo.edges[0].a))
                        if (B2.loop1ContainsLoop2(newLoopInfo.loop, newLoopInfo.ccw, subLoopInfo.loop, subLoopInfo.ccw, surface)) {
                            newLoopInfo.subloops.push(subLoopInfo);
                            loopInfos.splice(i, 1); // remove it
                        }
                    }
                    loopInfos.push(newLoopInfo);
                }
            }
        }
        function newFacesRecursive(loopInfo) {
            newFaces.push(new faceConstructor(surface, loopInfo.ccw ? loopInfo.loop : Edge.reversePath(loopInfo.loop), loopInfo.subloops.map(sl => sl.ccw ? Edge.reversePath(sl.loop) : sl.loop)));
            loopInfo.subloops.forEach(sl => sl.subloops.forEach(sl2 => newFacesRecursive(sl2)));
        }
        const newFaces = [];
        const topLevelLoops = [];
        loops.forEach(loop => placeRecursively({
            loop: loop,
            ccw: surface.edgeLoopCCW(loop),
            subloops: [],
        }, topLevelLoops));
        topLevelLoops.forEach(tll => newFacesRecursive(tll));
        return newFaces;
    }
    //fromLoops(loops: Edge[][], surface: Surface) {
    //	type LoopInfo = {loop: Edge[], ccw: boolean, subloops: LoopInfo[]}
    //	function placeRecursively(newLoopInfo: LoopInfo, loopInfos: LoopInfo[]) {
    //		if (loopInfos.length == 0) {
    //			loopInfos.push(newLoopInfo)
    //		} else {
    //			const subLoopInfo = loopInfos.find(loopInfo => B2.loop1ContainsLoop2(loopInfo.loop, loopInfo.ccw,
    // newLoopInfo.loop, newLoopInfo.ccw, surface)) if (subLoopInfo) { placeRecursively(newLoopInfo,
    // subLoopInfo.subloops) } else { // newLoopInfo isnt contained by any other subLoopInfo for (let i =
    // loopInfos.length; --i >= 0;) { const subLoopInfo = loopInfos[i] //console.log('cheving subLoopInfo',
    // surface.loopContainsPoint(newLoopInfo.edges, subLoopInfo.edges[0].a)) if
    // (B2.loop1ContainsLoop2(newLoopInfo.loop, subLoopInfo.loop, surface)) { newLoopInfo.subloops.push(subLoopInfo)
    // loopInfos.splice(i, 1) // remove it } } loopInfos.push(newLoopInfo) } } }  function newFacesRecursive(loopInfo:
    // LoopInfo): void { // CW loops can be top level, if they are holes in the original face not contained in the new
    // face if (loopInfo.ccw) { if (loopInfo.subloops.every(sl => !sl.ccw)) { const newFace = new
    // faceConstructor(surface, loopInfo.loop, loopInfo.subloops.map(sl => sl.loop)) newFaces.push(newFace)
    // loopInfo.subloops.forEach(sl => sl.subloops.forEach(slsl => slsl.ccw && newFacesRecursive(slsl))) } else {
    // loopInfo.subloops.forEach(sl => sl.ccw && newFacesRecursive(sl)) } } }  const newFaces: Face[] = [] const
    // topLevelLoops:LoopInfo[] = [] loops.forEach(loop => placeRecursively({loop: loop, ccw:
    // surface.edgeLoopCCW(loop), subloops: []}, topLevelLoops)) topLevelLoops.forEach(tll => newFacesRecursive(tll))
    // return newFaces }
    static create(surface, faceEdges, holes, faceName, info) {
        return surface instanceof PlaneSurface$1
            ? new PlaneFace(surface, faceEdges, holes, faceName, info)
            : new RotationFace(surface, faceEdges, holes, faceName, info);
    }
    intersectFace(face2, thisBrep, face2Brep, faceMap, thisEdgePoints, otherEdgePoints, checkedPairs) {
        //thisEdgePoints = {
        //   get(key) {
        //       return _thisEdgePoints.get(key)
        //    },
        //    set(key, value) {
        //       assert(thisBrep.edgeFaces.get(key))
        //        _thisEdgePoints.set(key, value)
        //    }
        //}
        function hasPair(a, b) {
            return checkedPairs.has(new javasetmap_ts.Pair(a, b));
        }
        function addPair(a, b) {
            return checkedPairs.add(new javasetmap_ts.Pair(a, b));
        }
        /**
         * @param newEdge generated segment
         * @param col1 if newEdge is colinear to an edge of this, the edge in question
         * @param col2 same for face2
         * @return whether new edge was added.
         */
        function handleNewEdge(newEdge, col1, col2) {
            if (!col1 && !col2) {
                let correctDir = face.surface.normalP(newEdge.a).cross(face2.surface.normalP(newEdge.a));
                if (correctDir.likeO()) {
                    const t = ts3dutils.lerp(newEdge.aT, newEdge.bT, 1 / ts3dutils.GOLDEN_RATIO), p = newEdge.curve.at(t);
                    correctDir = face.surface.normalP(p).cross(face2.surface.normalP(p));
                }
                if (!correctDir.likeO()) {
                    if (correctDir.dot(newEdge.aDir) < 0) {
                        newEdge = newEdge.flipped();
                    }
                    ts3dutils.mapPush(faceMap, face, newEdge);
                    ts3dutils.mapPush(faceMap, face2, newEdge.flipped());
                }
                else {
                    const p = newEdge.a;
                    const plane = P3.normalOnAnchor(newEdge.aDir, p);
                    const up = face.surface.normalP(p);
                    const sameDir = up.dot(face2.surface.normalP(p)) > 0;
                    const canonDir = plane.normal1.cross(up);
                    const curve = face.surface.isCurvesWithPlane(plane)[0], curveT = curve.pointT(p), curveDir = sign$9(canonDir.dot(curve.tangentAt(curveT)));
                    const curve2 = face2.surface.isCurvesWithPlane(plane)[0], curve2T = curve2.pointT(p), curve2Dir = sign$9(canonDir.dot(curve.tangentAt(curve2T)));
                    const foo = curve.diff(curveT, EPS * curveDir).dot(up);
                    const foo2 = curve2.diff(curve2T, EPS * curve2Dir).dot(up);
                    if (foo2 < foo) {
                        ts3dutils.mapPush(faceMap, face2, sameDir ? newEdge.flipped() : newEdge);
                    }
                    if (up.dot(face2.surface.normalP(p)) < 0 == foo2 < foo) {
                        ts3dutils.mapPush(faceMap, face, newEdge.flipped());
                    }
                    const bar = curve.diff(curveT, EPS * curveDir).dot(up);
                    const bar2 = curve2.diff(curve2T, EPS * curve2Dir).dot(up);
                    if (bar2 < bar) {
                        ts3dutils.mapPush(faceMap, face2, sameDir ? newEdge : newEdge.flipped());
                    }
                    if (sameDir != bar2 < bar) {
                        ts3dutils.mapPush(faceMap, face, newEdge);
                    }
                }
                return true;
            }
            function handleEdgeInFace(col1, col2, face, face2, thisBrep, face2Brep, coplanarSameIsInside, has, add) {
                if (col1 && !col2) {
                    if (hasPair(col1.getCanon(), face2))
                        return false;
                    //add(col1.getCanon(), face2)
                    const surface2 = face2.surface;
                    // NB: a new edge is inserted even though it may be the same as an old one
                    // however it indicates that it intersects the other volume here, i.e. the old edge cannot
                    // be counted as 'inside' for purposes of reconstitution
                    thisBrep.edgeFaces.get(col1.getCanon()).forEach(faceInfo => {
                        //const dot = snap0(surface2.normal1.dot(faceInfo.inside))
                        //if (dot == 0 ? !coplanarSameIsInside : dot < 0) {
                        const pointsInsideFace = fff(faceInfo, face2.surface);
                        const edgeInside = pointsInsideFace == INSIDE || !coplanarSameIsInside && pointsInsideFace == COPLANAR_SAME;
                        const pushEdge = faceInfo.edge.tangentAt(faceInfo.edge.curve.pointT(newEdge.a)).like(newEdge.aDir)
                            ? newEdge
                            : newEdge.flipped();
                        ts3dutils.assert(faceInfo.edge.tangentAt(faceInfo.edge.curve.pointT(pushEdge.a)).like(pushEdge.aDir));
                        edgeInside && ts3dutils.mapPush(faceMap, faceInfo.face, pushEdge);
                    });
                    const surface2NormalAtNewEdgeA = surface2.normalP(newEdge.a);
                    const newEdgeInside = surface2NormalAtNewEdgeA.cross(newEdge.aDir);
                    const sVEF1 = splitsVolumeEnclosingFacesP(thisBrep, col1.getCanon(), newEdge.a, newEdgeInside, surface2NormalAtNewEdgeA);
                    let addNewEdge, addNewEdgeFlipped;
                    if (addNewEdge = sVEF1 == INSIDE || coplanarSameIsInside && sVEF1 == COPLANAR_SAME) {
                        ts3dutils.mapPush(faceMap, face2, newEdge);
                    }
                    const sVEF2 = splitsVolumeEnclosingFacesP(thisBrep, col1.getCanon(), newEdge.a, newEdgeInside.negated(), surface2NormalAtNewEdgeA);
                    if (addNewEdgeFlipped = sVEF2 == INSIDE || coplanarSameIsInside && sVEF2 == COPLANAR_SAME) {
                        ts3dutils.mapPush(faceMap, face2, newEdge.flipped());
                    }
                    if (addNewEdge || addNewEdgeFlipped || sVEF1 == COPLANAR_SAME && sVEF2 == INSIDE || sVEF2 == COPLANAR_SAME && sVEF1 == INSIDE) {
                        return true;
                    }
                }
            }
            const c1 = handleEdgeInFace(col1, col2, face, face2, thisBrep, face2Brep, false, hasPair, addPair);
            const c2 = handleEdgeInFace(col2, col1, face2, face, face2Brep, thisBrep, true, (a, b) => hasPair(b, a), (a, b) => addPair(b, a));
            if (c1 || c2)
                return true;
            if (col1 && col2) {
                if (hasPair(col1.getCanon(), col2.getCanon()))
                    return false;
                addPair(col1.getCanon(), col2.getCanon());
                function handleColinearEdgeFaces(col1, col2, thisBrep, face2Brep, coplanarSameIsInside, thisEdgePoints, has, add) {
                    // not entirely sure for what i had the dirInsides in?
                    //const aDirNegatedInside = (newEdge.a.like(col2.a) || newEdge.a.like(col2.b)) &&
                    // splitsVolumeEnclosingCone(face2Brep, newEdge.a, newEdge.aDir.negated()) == INSIDE const
                    // bDirInside = (newEdge.b.like(col2.a) || newEdge.b.like(col2.b)) &&
                    // splitsVolumeEnclosingCone(face2Brep, newEdge.b, newEdge.bDir) == INSIDE
                    for (const faceInfo of thisBrep.edgeFaces.get(col1.getCanon())) {
                        const sVEF = splitsVolumeEnclosingFaces(face2Brep, col2.getCanon(), faceInfo.inside, faceInfo.normalAtCanonA);
                        const edgeInside = sVEF == INSIDE || coplanarSameIsInside && sVEF == COPLANAR_SAME;
                        const pushEdge = (faceInfo.edge.aDir.like(newEdge.aDir)) ? newEdge : newEdge.flipped();
                        if (edgeInside) {
                            ts3dutils.mapPush(faceMap, faceInfo.face, pushEdge);
                            const aT = col1.getCanon().curve.pointT(newEdge.a);
                            if (!ts3dutils.eq(aT, col1.aT) && !ts3dutils.eq(aT, col1.bT)) {
                                // newEdge.a is in center of col1
                                if (splitsVolumeEnclosingCone2(face2Brep, newEdge.a, newEdge.curve, newEdge.aT, -Math.sign(newEdge.deltaT())) == INSIDE) {
                                    ts3dutils.mapPush(thisEdgePoints, col1.getCanon(), {
                                        p: newEdge.a,
                                        edgeT: aT
                                    });
                                }
                            }
                            const bT = col1.getCanon().curve.pointT(newEdge.b);
                            if (!ts3dutils.eq(bT, col1.aT) && !ts3dutils.eq(bT, col1.bT)) {
                                if (splitsVolumeEnclosingCone2(face2Brep, newEdge.b, newEdge.curve, newEdge.bT, Math.sign(newEdge.deltaT())) == INSIDE) {
                                    ts3dutils.mapPush(thisEdgePoints, col1.getCanon(), {
                                        p: newEdge.b,
                                        edgeT: bT
                                    });
                                }
                            }
                            
                        }
                    }
                }
                handleColinearEdgeFaces(col1, col2, thisBrep, face2Brep, true, thisEdgePoints, hasPair, addPair);
                handleColinearEdgeFaces(col2, col1, face2Brep, thisBrep, false, otherEdgePoints, (a, b) => hasPair(b, a), (a, b) => addPair(b, a));
                return false;
            }
        }
        // what needs to be generated: new edges on face
        // points on edges where they are cut by faces so that sub edges will be generated for loops
        // points on ends of edges where the edge will be an edge in the new volume where it goes from A to B
        //         you don't want those to be marked as 'inside', otherwise invalid faces will be added
        // if a face cuts a corner, nothing needs to be done, as that alone does not limit what adjacent faces will be
        function handleEndPoint(a, b, newEdge) {
            // ends in the middle of b's face
            if (a && !b) {
                if (!a.colinear && a.edgeT != a.edge.aT && a.edgeT != a.edge.bT) {
                    ts3dutils.mapPush(thisEdgePoints, a.edge.getCanon(), a);
                    ts3dutils.assert(a.edge.isValidT(a.edgeT));
                }
                // else colinear segment ends in middle of other face, do nothing
            }
            // ends in the middle of a's face
            if (b && !a) {
                if (!b.colinear && b.edgeT != b.edge.aT && b.edgeT != b.edge.bT) {
                    ts3dutils.mapPush(otherEdgePoints, b.edge.getCanon(), b);
                    ts3dutils.assert(b.edge.isValidT(b.edgeT));
                }
                // else colinear segment ends in middle of other face, do nothing
            }
            if (a && b) {
                ts3dutils.assert(a.colinear || b.colinear || ts3dutils.eq(a.t, b.t));
                // if a or b is colinear the correct points will already have been added to the edge by handleNewEdge
                // segment starts/ends on edge/edge intersection
                function foo(a, b, face, face2, thisPlane, face2Plane, thisBrep, face2Brep, first, thisEdgePoints) {
                    if (!a.colinear && a.edgeT != a.edge.aT && a.edgeT != a.edge.bT) {
                        //if (!hasPair(a.edge.getCanon(), b.edge.getCanon())) {
                        addPair(a.edge.getCanon(), b.edge.getCanon());
                        // ends on a, on colinear segment b bT != a.edge.bT &&
                        // b can be colinear, so edgeT == aT is possible
                        if (a.p.like(b.edge.a) || a.p.like(b.edge.b)) {
                            const corner = a.p.like(b.edge.a) ? b.edge.a : b.edge.b;
                            // face2brep corner on edge
                            const sVEC1 = splitsVolumeEnclosingCone2(face2Brep, corner, a.edge.curve, a.edgeT, 1);
                            const sVEC2 = splitsVolumeEnclosingCone2(face2Brep, corner, a.edge.curve, a.edgeT, -1);
                            // if either of these return ALONG_EDGE_OR_PLANE, then the breps share a colinear edge
                            if (INSIDE == sVEC1 || INSIDE == sVEC2) {
                                ts3dutils.mapPush(thisEdgePoints, a.edge.getCanon(), a);
                                ts3dutils.assert(a.edge.isValidT(a.edgeT));
                            }
                        }
                        else {
                            // edge / edge center intersection
                            // todo: is this even necessary considering we add edges anyway? i think so...
                            // const testVector =
                            // a.edge.tangentAt(a.edgeT).rejectedFrom(b.edge.tangentAt(b.edge.curve.pointT(a.p)))
                            // assert(!testVector.likeO())
                            const sVEF1 = splitsVolumeEnclosingFacesP2(face2Brep, b.edge.getCanon(), a.p, a.edge.curve, a.edgeT, 1, thisPlane.normalP(a.p));
                            const sVEF2 = splitsVolumeEnclosingFacesP2(face2Brep, b.edge.getCanon(), a.p, a.edge.curve, a.edgeT, -1, thisPlane.normalP(a.p));
                            if (INSIDE == sVEF1 || INSIDE == sVEF2) {
                                ts3dutils.mapPush(thisEdgePoints, a.edge.getCanon(), a);
                                ts3dutils.assert(a.edge.isValidT(a.edgeT));
                            }
                        }
                        //}
                    }
                }
                foo(a, b, face, face2, surface, surface2, thisBrep, face2Brep, true, thisEdgePoints);
                foo(b, a, face2, face, surface2, surface, face2Brep, thisBrep, false, otherEdgePoints);
            }
        }
        ts3dutils.assertInst(Face, face2);
        const face = this;
        const surface = face.surface, surface2 = face2.surface;
        if (!this.getAABB().fuzzyTouchesAABB(face2.getAABB())) {
            return;
        }
        if (surface.isCoplanarTo(surface2)) {
            return;
        }
        const isCurves = surface.isCurvesWithSurface(surface2);
        if (0 == isCurves.length) {
            return;
        }
        for (const isCurve of isCurves) {
            const t = (isCurve.tMin + isCurve.tMax) / 2, p = isCurve.at(t), dp = isCurve.tangentAt(t);
            const normal1 = surface.normalP(p), normal2 = surface2.normalP(p), dp2 = normal1.cross(normal2);
            ts3dutils.assert(surface.containsCurve(isCurve));
            ts3dutils.assert(surface2.containsCurve(isCurve));
            if (!dp2.likeO()) {
                //assert(dp2.dot(dp) > 0)
                // TODO assert(dp2.isParallelTo(dp))
            }
        }
        for (let isCurveIndex = 0; isCurveIndex < isCurves.length; isCurveIndex++) {
            // get intersections of newCurve with other edges of face and face2
            const isCurve = isCurves[isCurveIndex];
            const ps1 = face.edgeISPsWithSurface(isCurve, face2.surface);
            const ps2 = face2.edgeISPsWithSurface(isCurve, face.surface);
            // for non-endless curves, e.g. ellipses, the intersections of the faces can be non-zero, even if one of
            // the faces doesn't register any points on the curve. For example, if a cylinder is cut entirely by a
            // plane face (all its edges around the cylinder), then the face will contain the entire curve and
            // 'ps' for the plane face will be empty
            // TODO: behavior when curves touch face?
            // !! start in does depend on insideDir... TODO
            ts3dutils.assertf(() => (0 == ps1.length) || !ts3dutils.eq0(ps1[0].insideDir.dot(isCurve.tangentAt(ps1[0].t))), () => ps1[0].insideDir.dot(isCurve.tangentAt(ps1[0].t)));
            ts3dutils.assertf(() => (0 == ps2.length) || !ts3dutils.eq0(ps2[0].insideDir.dot(isCurve.tangentAt(ps2[0].t))), () => ps2[0].insideDir.dot(isCurve.tangentAt(ps2[0].t)));
            function startsInside(ps, face) {
                if (0 == ps.length) {
                    return isFinite(isCurve.tMin) && face.containsPoint2(isCurve.at(isCurve.tMin)) == exports.PointVsFace.INSIDE;
                }
                else {
                    return ps[0].insideDir.dot(isCurve.tangentAt(ps[0].t)) < 0;
                }
            }
            // they can't both be empty currently
            // they can't both start 'inside'
            let in1 = startsInside(ps1, face);
            let in2 = startsInside(ps2, face2);
            if (0 == ps1.length && !in1 || 0 == ps2.length && !in2) {
                continue;
            }
            //assert(!in1 || !in2)
            let col1, col2;
            let i = 0, j = 0, last;
            let startP = in1 && in2 && isCurve.at(isCurve.tMin), startDir, startT = isCurve.tMin, startA, startB;
            while (i < ps1.length || j < ps2.length) {
                ts3dutils.assert(i <= ps1.length);
                ts3dutils.assert(j <= ps2.length);
                const a = ps1[i], b = ps2[j];
                ts3dutils.assert(a || b);
                if (j == ps2.length || i < ps1.length && ts3dutils.lt(a.t, b.t)) {
                    last = a;
                    in1 = !in1;
                    a.used = true;
                    in1 && (col1 = a.colinear && a);
                    i++;
                }
                else if (i == ps1.length || ts3dutils.gt(a.t, b.t)) {
                    last = b;
                    b.used = true;
                    in2 = !in2;
                    in2 && (col2 = b.colinear && b);
                    j++;
                }
                else {
                    last = a;
                    a.used = true;
                    b.used = true;
                    in1 = !in1;
                    in2 = !in2;
                    //if (in1 == in2) {
                    in1 && (col1 = a.colinear && a);
                    in2 && (col2 = b.colinear && b);
                    //}
                    i++;
                    j++;
                }
                if (startP && !(in1 && in2)) {
                    // segment end
                    startDir = isCurve.tangentAt(startT);
                    if (ts3dutils.eq(startT, last.t)) {
                        startP = undefined;
                        continue;
                    }
                    ts3dutils.assert(ts3dutils.lt(startT, last.t));
                    startT > last.t && (startDir = startDir.negated());
                    let endDir = isCurve.tangentAt(last.t);
                    startT > last.t && (endDir = endDir.negated());
                    const newEdge = Edge.create(isCurve, startP, last.p, startT, last.t, undefined, startDir, endDir, 'genseg' + getGlobalId());
                    startP = undefined;
                    if (handleNewEdge(newEdge, col1 && col1.edge, col2 && col2.edge)) {
                        handleEndPoint(startA || col1, startB || col2, newEdge);
                        handleEndPoint(a && a.used && a || col1, b && b.used && b || col2, newEdge);
                    }
                }
                else if (in1 && in2) {
                    // new segment just started
                    startP = last.p;
                    startDir = last.insideDir;
                    startT = last.t;
                    startA = a && a.used && a;
                    startB = b && b.used && b;
                }
            }
            if (in1 && in2 && startT !== isCurve.tMax) {
                const endT = isCurve.tMax;
                startDir = isCurve.tangentAt(startT);
                startT > endT && (startDir = startDir.negated());
                let endDir = isCurve.tangentAt(endT);
                startT > endT && (endDir = endDir.negated());
                const newEdge = Edge.create(isCurve, startP, isCurve.at(endT), startT, endT, undefined, startDir, endDir, 'genseg' + getGlobalId());
                if (handleNewEdge(newEdge, col1 && col1.edge, col2 && col2.edge)) {
                    handleEndPoint(startA || col1, startB || col2, newEdge);
                }
            }
        }
        face.getAllEdges().forEach(edge => {
            checkedPairs.add(new javasetmap_ts.Pair(edge.getCanon(), face2));
        });
        face2.getAllEdges().forEach(edge => {
            checkedPairs.add(new javasetmap_ts.Pair(edge.getCanon(), face));
        });
    }
    edgeISPsWithSurface(isCurve, surface2) {
        const face = this;
        const surface = face.surface;
        const loops = face.holes.concat([face.contour]);
        const ps = [];
        for (const loop of loops) {
            const colinearEdges = loop.map(edge => edge.curve.isColinearTo(isCurve));
            //const colinearSides = loop.map((edge, edgeIndex) => -1 != colinearEdges[edgeIndex]
            //            && -sign(isCurves[colinearEdges[edgeIndex]].tangentAt(edge.aT).dot(edge.aDir)))
            for (let edgeIndex = 0; edgeIndex < loop.length; edgeIndex++) {
                const edge = loop[edgeIndex];
                const nextEdgeIndex = (edgeIndex + 1) % loop.length, nextEdge = loop[nextEdgeIndex];
                //console.log(edge.toSource()) {p:V3(2, -2.102, 0),
                if (colinearEdges[edgeIndex]) {
                    if (isCurve.containsPoint(edge.a)) {
                        const prevEdgeIndex = (edgeIndex - 1 + loop.length) % loop.length, prevEdge = loop[prevEdgeIndex];
                        const curveAT = isCurve.pointT(edge.a);
                        const colinearOutA = edge.aDir.cross(surface.normalP(edge.a));
                        if (!colinearEdges[prevEdgeIndex] && dotCurve2(prevEdge.curve, prevEdge.bT, colinearOutA, -sign$9(prevEdge.deltaT())) > 0) {
                            ps.push({
                                p: prevEdge.b,
                                insideDir: edge.aDir.negated(),
                                t: curveAT,
                                edge: prevEdge,
                                edgeT: prevEdge.bT,
                                colinear: false,
                            });
                        }
                        ps.push({
                            p: edge.a,
                            insideDir: edge.aDir,
                            t: curveAT,
                            edge: edge,
                            edgeT: edge.aT,
                            colinear: true,
                        });
                    }
                    if (isCurve.containsPoint(edge.b)) {
                        const curveBT = isCurve.pointT(edge.b);
                        const colinearOutB = edge.bDir.cross(surface.normalP(edge.b));
                        if (!colinearEdges[nextEdgeIndex] && dotCurve2(nextEdge.curve, nextEdge.aT, colinearOutB, sign$9(nextEdge.deltaT())) > 0) {
                            ps.push({
                                p: edge.b,
                                insideDir: edge.bDir,
                                t: curveBT,
                                edge: nextEdge,
                                edgeT: nextEdge.aT,
                                colinear: false,
                            });
                        }
                        ps.push({
                            p: edge.b,
                            insideDir: edge.bDir.negated(),
                            t: curveBT,
                            edge: edge,
                            edgeT: edge.bT,
                            colinear: true,
                        });
                    }
                }
                else {
                    const edgeTs = edge.edgeISTsWithSurface(surface2);
                    for (const edgeT of edgeTs) {
                        const p = edge.curve.at(edgeT);
                        if (!isCurve.containsPoint(p))
                            continue;
                        const curveT = isCurve.pointT(p);
                        ts3dutils.assert(!isNaN(curveT));
                        const insideDir = edge.tangentAt(edgeT).cross(surface.normalP(p)).negated();
                        const isTangent = isCurve.tangentAt(curveT);
                        const dirFactor = sign$9(isTangent.dot(edge.curve.tangentAt(edgeT)));
                        const normVector = surface2.normalP(p);
                        //if(!eq0(insideDir.dot(isTangent))) {
                        // Edge.edgeISTsWithSurface returns snapped values, so comparison with == is ok:
                        if (edgeT == edge.bT) {
                            // endpoint lies on intersection line
                            if (!colinearEdges[nextEdgeIndex]) {
                                if (!ts3dutils.eq(curveT, isCurve.tMax)) {
                                    const pointsToInside = this.pointsToInside3(edge.b, isCurve, curveT, 1);
                                    ts3dutils.assert(pointsToInside != exports.PointVsFace.ON_EDGE);
                                    if (exports.PointVsFace.INSIDE == pointsToInside) {
                                        ps.push({
                                            p: edge.b,
                                            insideDir: isTangent,
                                            t: curveT,
                                            edge: edge,
                                            edgeT: edge.bT,
                                            colinear: false,
                                        });
                                    }
                                }
                                if (!ts3dutils.eq(curveT, isCurve.tMin)) {
                                    const pointsToInside = this.pointsToInside3(edge.b, isCurve, curveT, -1);
                                    ts3dutils.assert(pointsToInside != exports.PointVsFace.ON_EDGE);
                                    if (exports.PointVsFace.INSIDE == pointsToInside) {
                                        ps.push({
                                            p: edge.b,
                                            insideDir: isTangent.negated(),
                                            t: curveT,
                                            edge: edge,
                                            edgeT: edge.bT,
                                            colinear: false,
                                        });
                                    }
                                }
                                //let thisSide = -normVector.dot(edge.bDir)
                                //if (eq0(thisSide)) {
                                //    // advanced test
                                //    const dir = -sign(edge.deltaT())
                                //    const iscd = isCurve.at(curveT).to(isCurve.at(curveT + dir * dirFactor *
                                // eps)).dot(normVector) const ecd = edge.curve.at(edgeT).to(edge.curve.at(edgeT + dir
                                // * eps)).dot(normVector) thisSide = sign(ecd - iscd) } let nextSide =
                                // normVector.dot(nextEdge.aDir) if (eq0(nextSide)) { // advanced test const dirFactor
                                // = sign(snap0(isTangent.dot(nextEdge.curve.tangentAt(nextEdge.aT)))) assert(dirFactor
                                // !== 0) const dir = sign(nextEdge.deltaT()) const iscd =
                                // isCurve.at(curveT).to(isCurve.at(curveT + dir * dirFactor * eps)).dot(normVector)
                                // const ecd = nextEdge.curve.at(nextEdge.aT).to(nextEdge.curve.at(nextEdge.aT + dir *
                                // eps)).dot(normVector) nextSide = sign(ecd - iscd) } if (nextSide < 0 || thisSide <
                                // 0) { assert(!eq0(insideDir.dot(isTangent))) // next segment is not colinear and ends
                                // on different side ps.push({ p: edge.b, insideDir: insideDir, t: curveT, edge: edge,
                                // edgeT: edge.bT, colinear: false}) }
                            }
                        }
                        else if (edgeT != edge.aT) {
                            // edge crosses/touches an intersection curve, neither starts nor ends on it
                            if (ts3dutils.eq0(insideDir.dot(isTangent))) {
                                const dirFactor = sign$9(isTangent.dot(edge.curve.tangentAt(edgeT)));
                                const eps = 1e-4;
                                for (const dir of [-1, 1]) {
                                    if (-1 == dir * dirFactor && edgeT == edge.minT ||
                                        1 == dir * dirFactor && edgeT == edge.maxT ||
                                        -1 == dir && curveT == isCurve.tMin ||
                                        1 == dir && curveT == isCurve.tMax)
                                        continue;
                                    const iscd = isCurve.at(curveT).to(isCurve.at(curveT + dir * eps)).dot(insideDir);
                                    const ecd = edge.curve.at(edgeT).to(edge.curve.at(edgeT + dir * dirFactor * eps)).dot(insideDir);
                                    if (iscd - ecd > 0) {
                                        ps.push({
                                            p,
                                            insideDir: isTangent.times(dir),
                                            t: curveT,
                                            edge: edge,
                                            edgeT: edgeT,
                                            colinear: false,
                                        });
                                    }
                                }
                            }
                            else {
                                ps.push({
                                    p: p,
                                    insideDir: insideDir,
                                    t: curveT,
                                    edge: edge,
                                    edgeT: edgeT,
                                    colinear: false,
                                });
                            }
                        }
                        //} else {
                        //
                        //	const dirFactor = sign(isTangent.dot(edge.curve.tangentAt(edgeT)))
                        //	const eps = 1e-4
                        //	const normVector = surface2.normalP(p)
                        //	for (const dir of [-1, 1]) {
                        //		if (-1 == dir * dirFactor && edgeT == edge.minT ||
                        //			1 == dir * dirFactor && edgeT == edge.maxT ||
                        //			-1 == dir && curveT == isCurve.tMin ||
                        //			1 == dir && curveT == isCurve.tMax) continue
                        //		const iscd = isCurve.at(curveT).to(isCurve.at(curveT + dir * eps)).dot(normVector)
                        //		const ecd = edge.curve.at(edgeT).to(edge.curve.at(edgeT + dir * dirFactor *
                        // eps)).dot(normVector) if (iscd > ecd) { ps.push({p, insideDir: isTangent.times(dir *
                        // dirFactor), t: curveT, edge: edge, edgeT: edgeT, colinear: false}) } }
                        // curveVsSurface(isCurve, curveT, p, surface2) }
                    }
                }
            }
        }
        // duplicate 't's are ok, as sometimes a segment needs to stop and start again
        // should be sorted so that back facing ones are first
        ps.sort((a, b) => a.t - b.t || a.insideDir.dot(isCurve.tangentAt(a.t)));
        return ps;
    }
    transform(m4) {
        const mirroring = m4.isMirroring();
        const newEdges = Edge.reversePath(this.contour.map(e => e.transform(m4)), mirroring);
        const newHoles = this.holes.map(hole => Edge.reversePath(hole.map(e => e.transform(m4)), mirroring));
        return new this.constructor(this.surface.transform(m4), newEdges, newHoles, this.name, this.info);
    }
    flipped() {
        const newEdges = this.contour.map(e => e.flipped()).reverse();
        const newHoles = this.holes.map(hole => hole.map(e => e.flipped()).reverse());
        return new this.constructor(this.surface.flipped(), newEdges, newHoles, this.name, this.info);
    }
    toString() {
        return 'new ' + this.constructor.name + '(' + this.surface + ', [' + this.contour.map(e => '\n\t' + e).join() + ']'
            + this.holes.map(hole => '\n\t\thole: ' + hole.join()) + ')';
    }
    toSource() {
        return `new ${this.constructor.name}(${this.surface.toSource()}, [${this.contour.map(e => '\n\t' + e.toSource()).join(',')}], [${this.holes.map(hole => '[' + hole.map(e => '\n\t' + e.toSource()).join(',') + ']').join(',')}])`;
    }
    equals(obj) {
        function loopsEqual(a, b) {
            return a.length == b.length &&
                ts3dutils.arrayRange(0, a.length, 1)
                    .some(offset => a.every((edge, i) => edge.equals(b[(offset + i) % a.length])));
        }
        return this == obj ||
            Object.getPrototypeOf(this) == Object.getPrototypeOf(obj)
                && this.holes.length == obj.holes.length
                && loopsEqual(this.contour, obj.contour)
                && this.holes.every(hole => obj.holes.some(hole2 => loopsEqual(hole, hole2)));
    }
    hashCode() {
        function arrayHashCode(array) {
            let hashCode = 0;
            for (const val of array) {
                hashCode = hashCode * 31 + val | 0;
            }
            return hashCode;
        }
        function loopHashCode(loop) { return arrayHashCode(loop.map(edge => edge.hashCode()).sort(ts3dutils.MINUS)); }
        let hashCode = 0;
        hashCode = hashCode * 31 + arrayHashCode(this.holes.map(loop => loopHashCode(loop)).sort(ts3dutils.MINUS)) | 0;
        hashCode = hashCode * 31 + loopHashCode(this.contour) | 0;
        hashCode = hashCode * 31 + this.surface.hashCode() | 0;
        return hashCode;
    }
    likeFace(face2) {
        function loopsLike(a, b) {
            return a.length == b.length &&
                ts3dutils.arrayRange(0, a.length, 1)
                    .some(offset => a.every((edge, i) => edge.like(b[(offset + i) % a.length])));
        }
        ts3dutils.assertInst(Face, face2);
        return this.surface.like(face2.surface)
            && this.holes.length == face2.holes.length
            && loopsLike(this.contour, face2.contour)
            && this.holes.every(hole => face2.holes.some(hole2 => loopsLike(hole, hole2)));
    }
    getAllEdges() {
        return this.allEdges;
    }
    addEdgeLines(mesh) {
        ts3dutils.assert(false, 'buggy, fix');
        const vertices = this.contour.flatMap(edge => edge.getVerticesNo0()), mvl = mesh.vertices.length;
        for (let i = 0; i < vertices.length; i++) {
            mesh.vertices.push(vertices[i]);
            mesh.LINES.push(mvl + i, mvl + (i + 1) % vertices.length);
        }
    }
    containsPoint(p) {
        ts3dutils.assertVectors(p);
        return this.surface.loopContainsPoint(this.contour, p) != exports.PointVsFace.OUTSIDE
            && !this.holes.some(hole => this.surface.loopContainsPoint(hole, p) != exports.PointVsFace.OUTSIDE);
    }
    containsPoint2(p) {
        ts3dutils.assertVectors(p);
        const contourContainsPoint = this.surface.loopContainsPoint(this.contour, p);
        if (contourContainsPoint != exports.PointVsFace.INSIDE)
            return contourContainsPoint;
        for (const hole of this.holes) {
            const loopContainsPoint = this.surface.loopContainsPoint(hole, p);
            if (loopContainsPoint != exports.PointVsFace.OUTSIDE) {
                return loopContainsPoint == exports.PointVsFace.ON_EDGE ? exports.PointVsFace.ON_EDGE : exports.PointVsFace.OUTSIDE;
            }
        }
        return exports.PointVsFace.INSIDE;
    }
    /**
     *
     * @param line
     * @returns t param of the line if there is an intersection, NaN otherwise
     */
    intersectsLine(line) {
        ts3dutils.assertInst(L3$1, line);
        if (!this.getAABB().intersectsLine(line))
            return NaN;
        const containedIntersectionsTs = this.surface.isTsForLine(line).filter(t => this.containsPoint(line.at(t)));
        const nearestPointT = containedIntersectionsTs.withMax(t => -t);
        return undefined != nearestPointT ? nearestPointT : NaN;
    }
    toMesh() {
        const mesh = new tsgl.Mesh()
            .addIndexBuffer('TRIANGLES')
            .addIndexBuffer('LINES')
            .addVertexBuffer('normals', 'LGL_Normal');
        this.addToMesh(mesh);
        //mesh.compile()
        return mesh;
    }
    zDirVolume() {
        return this.surface.zDirVolume(this.getAllEdges());
    }
    calcArea() {
        return this.surface.calculateArea(this.getAllEdges());
    }
    getLoops() {
        return this.holes.concat(this.contour);
    }
    getAABB() {
        return this.aabb || (this.aabb = ts3dutils.AABB.forAABBs(this.contour.map(e => e.getAABB())));
    }
    pointsToInside3(p, curve, curveT, dir) {
        const eps = 1e-6;
        const normal = this.surface.normalP(p);
        const curveTangent = curve.tangentAt(curveT).times(dir);
        const up = normal.cross(curveTangent);
        const ecd = curve.at(curveT).to(curve.at(curveT + dir * eps)).dot(up);
        let minValue = Infinity, result, advanced = false;
        for (const edge of this.getAllEdges()) {
            const aEqP = edge.a.like(p), bEqP = edge.b.like(p);
            ts3dutils.assert(aEqP == edge.a.like(p));
            ts3dutils.assert(bEqP == edge.b.like(p));
            if (!aEqP && !bEqP)
                continue;
            const edgeTangent = aEqP ? edge.aDir : edge.bDir.negated();
            const angle = curveTangent.angleRelativeNormal(edgeTangent, normal);
            if (ts3dutils.eq0(angle)) {
                if (curve.isColinearTo(edge.curve)) {
                    return exports.PointVsFace.ON_EDGE;
                }
                const edgeT = aEqP ? edge.aT : edge.bT;
                const edgeDir = (aEqP ? 1 : -1) * sign$9(edge.deltaT());
                const iscd = edge.curve.diff(edgeT, edgeDir * eps).dot(up);
                //const iscd = edge.curve.at(edgeT).to(curve.at(edgeT + edgeDir * eps)).dot(up)
                const diff = iscd - ecd;
                if (diff > 0 && (!advanced || diff < minValue)) {
                    advanced = true;
                    minValue = diff;
                    result = aEqP ? exports.PointVsFace.OUTSIDE : exports.PointVsFace.INSIDE;
                }
            }
            else if (!advanced) {
                const angle2 = (angle + ts3dutils.TAU) % ts3dutils.TAU;
                if (angle2 < minValue) {
                    minValue = angle2;
                    result = aEqP ? exports.PointVsFace.OUTSIDE : exports.PointVsFace.INSIDE;
                }
            }
        }
        if (result == undefined)
            throw new Error();
        return result;
    }
    pointsToInside2(p, dir) {
        return this.pointsToInside3(p, L3$1.anchorDirection(p, dir), 0, 1);
        //const normal = this.surface.normalP(p)
        //let minAngle = Infinity, inOut = false
        //function test(v, b) {
        //	const angle = (dir.angleRelativeNormal(v, normal) + TAU + NLA_PRECISION / 2) % TAU
        //	if (angle <= 2 * NLA_PRECISION) {
        //		return true
        //	}
        //	if (angle < minAngle) {
        //		minAngle = angle
        //		inOut = b
        //	}
        //}
        //for (const edge of this.getAllEdges()) {
        //	assert(edge.a.equals(p) || !edge.a.like(p))
        //	assert(edge.b.equals(p) || !edge.b.like(p))
        //	if (edge.a.equals(p) && test(edge.aDir, false)) return PointVsFace.ON_EDGE
        //	if (edge.b.equals(p) && test(edge.bDir.negated(), true)) return PointVsFace.ON_EDGE
        //}
        //return inOut ? PointVsFace.INSIDE : PointVsFace.OUTSIDE
    }
}
class PlaneFace extends Face {
    constructor(p, contour, holes, name, info) {
        ts3dutils.assert(p instanceof P3 || p instanceof PlaneSurface$1);
        super(p instanceof P3 ? new PlaneSurface$1(p) : p, contour, holes, name, info);
    }
    static forVertices(planeSurface, vs, ...holeVss) {
        const _planeSurface = planeSurface instanceof P3 ? new PlaneSurface$1(planeSurface) : planeSurface;
        ts3dutils.assert(ts3dutils.isCCW(vs, _planeSurface.plane.normal1), 'isCCW(vs, planeSurface.plane.normal1)');
        const edges = StraightEdge.chain(vs);
        holeVss.forEach(vs => ts3dutils.assert(ts3dutils.doubleSignedArea(vs, _planeSurface.plane.normal1) >= 0, 'doubleSignedArea(vs, planeSurface.plane.normal1) >= 0'));
        const holes = holeVss.map(hvs => StraightEdge.chain(hvs));
        return new PlaneFace(planeSurface, edges, holes);
    }
    addToMesh(mesh) {
        const mvl = mesh.vertices.length;
        const normal = this.surface.plane.normal1;
        const vertices = this.contour.flatMap(edge => edge.getVerticesNo0());
        for (let i = 0; i < vertices.length; i++) {
            mesh.LINES.push(mvl + i, mvl + (i + 1) % vertices.length);
        }
        const holeStarts = [];
        this.holes.forEach(hole => {
            holeStarts.push(vertices.length);
            vertices.push(...hole.flatMap(edge => edge.getVerticesNo0()));
        });
        const triangles = triangulateVertices(normal, vertices, holeStarts).map(index => index + mvl);
        Array.prototype.push.apply(mesh.vertices, vertices);
        Array.prototype.push.apply(mesh.TRIANGLES, triangles);
        Array.prototype.push.apply(mesh.normals, ts3dutils.arrayFromFunction(vertices.length, () => normal));
    }
    intersectsLine(line) {
        ts3dutils.assertInst(L3$1, line);
        const lambda = line.isTWithPlane(this.surface.plane);
        if (!Number.isFinite(lambda)) {
            return NaN;
        }
        const inside = this.containsPoint(line.at(lambda));
        return inside ? lambda : NaN;
    }
    //intersectPlaneFace(face2: PlaneFace,
    //                   thisBrep: B2,
    //                   face2Brep: B2,
    //                   faceMap: Map<Face, Edge[]>,
    //                   thisEdgePoints: CustomMap<Edge, { edge: Edge, edgeT: number, p: V3, passEdge?: Edge }[]>,
    //                   otherEdgePoints: CustomMap<Edge, { edge: Edge, edgeT: number, p: V3, passEdge?: Edge }[]>,
    //                   checkedPairs: CustomSet<Pair<Equalable, Equalable>>) {
    //	assertInst(CustomMap, thisEdgePoints, otherEdgePoints)
    //
    //	function hasPair(a: Equalable, b: Equalable) {
    //		return checkedPairs.has(new Pair(a, b))
    //	}
    //	function addPair(a: Equalable, b: Equalable) {
    //		return checkedPairs.add(new Pair(a, b))
    //	}
    //
    //	/**
    //	 * @param newEdge generated segment
    //	 * @param col1 if newEdge is colinear to an edge of this, the edge in question
    //	 * @param col2 same for face2
    //	 */
    //	function handleNewEdge(newEdge: StraightEdge, col1: Edge, col2: Edge) {
    //		if (!col1 && !col2) {
    //			mapPush(faceMap, face, newEdge)
    //			mapPush(faceMap, face2, newEdge.flipped())
    //			return true
    //		}
    //		function handleEdgeInFace(col1, col2, face, face2, thisBrep, face2Brep, coplanarSameIsInside: boolean,
    // has, add) { if (col1 && !col2) { if (hasPair(col1.getCanon(), face2)) return  //add(col1.getCanon(), face2)
    // const face2Plane = face2.surface.plane  // NB: a new edge is inserted even though it may be the same as an old
    // one // however it indicates that it intersects the other volume here, i.e. the old edge cannot // be counted as
    // 'inside' for purposes of reconstitution thisBrep.edgeFaces.get(col1.getCanon()).forEach(faceInfo => { //const
    // dot = snap0(face2Plane.normal1.dot(faceInfo.inside)) //if (dot == 0 ? !coplanarSameIsInside : dot < 0) { const
    // pointsInsideFace = fff(faceInfo, face2.surface) const edgeInside = pointsInsideFace == INSIDE ||
    // !coplanarSameIsInside && pointsInsideFace == COPLANAR_SAME const pushEdge =
    // (faceInfo.edge.aDir.like(newEdge.aDir)) ? newEdge : newEdge.flipped()
    // assert(faceInfo.edge.aDir.like(pushEdge.aDir)) edgeInside && mapPush(faceMap, faceInfo.face, pushEdge) })  const
    // newEdgeInside = face2Plane.normal1.cross(newEdge.aDir) const sVEF1 = splitsVolumeEnclosingFaces(thisBrep,
    // col1.getCanon(), newEdgeInside, face2Plane.normal1) let addNewEdge, addNewEdgeFlipped if (addNewEdge = sVEF1 ==
    // INSIDE || coplanarSameIsInside && sVEF1 == COPLANAR_SAME) { mapPush(faceMap, face2, newEdge) } const sVEF2 =
    // splitsVolumeEnclosingFaces(thisBrep, col1.getCanon(), newEdgeInside.negated(), face2Plane.normal1) if
    // (addNewEdgeFlipped = sVEF2 == INSIDE || coplanarSameIsInside && sVEF2 == COPLANAR_SAME) { mapPush(faceMap,
    // face2, newEdge.flipped()) } if (addNewEdge || addNewEdgeFlipped || sVEF1 == COPLANAR_SAME && sVEF2 == INSIDE ||
    // sVEF2 == COPLANAR_SAME && sVEF1 == INSIDE) { return true } } } const c1 = handleEdgeInFace(col1, col2, face,
    // face2, thisBrep, face2Brep, false, hasPair, addPair) const c2 = handleEdgeInFace(col2, col1, face2, face,
    // face2Brep, thisBrep, true, (a, b) => hasPair(b, a), (a, b) => addPair(b, a)) if (c1 || c2) return true  if (col1
    // && col2) { if (hasPair(col1.getCanon(), col2.getCanon())) return  addPair(col1.getCanon(), col2.getCanon())
    // function handleColinearEdgeFaces(col1, col2, thisBrep, face2Brep, coplanarSameIsInside: boolean, thisEdgePoints,
    // has, add) { // not entirely sure for what i had the dirInsides in? //const aDirNegatedInside =
    // (newEdge.a.like(col2.a) || newEdge.a.like(col2.b)) && splitsVolumeEnclosingCone(face2Brep, newEdge.a,
    // newEdge.aDir.negated()) == INSIDE //const bDirInside = (newEdge.b.like(col2.a) || newEdge.b.like(col2.b)) &&
    // splitsVolumeEnclosingCone(face2Brep, newEdge.b, newEdge.bDir) == INSIDE
    // thisBrep.edgeFaces.get(col1.getCanon()).forEach(faceInfo => { const sVEF = splitsVolumeEnclosingFaces(face2Brep,
    // col2.getCanon(), faceInfo.inside, faceInfo.normalAtCanonA) const edgeInside = sVEF == INSIDE ||
    // coplanarSameIsInside && sVEF == COPLANAR_SAME const pushEdge = (faceInfo.edge.aDir.like(newEdge.aDir)) ? newEdge
    // : newEdge.flipped() edgeInside && mapPush(faceMap, faceInfo.face, pushEdge) }) } handleColinearEdgeFaces(col1,
    // col2, thisBrep, face2Brep, true, thisEdgePoints, hasPair, addPair) handleColinearEdgeFaces(col2, col1,
    // face2Brep, thisBrep, false, otherEdgePoints, (a, b) => hasPair(b, a), (a, b) => addPair(b, a)) } }   // what
    // needs to be generated: new edges on face // points on edges where they are cut by faces so that sub edges will
    // be generated for loops // points on ends of edges where the edge will be an edge in the new volume where it goes
    // from A to B //         you don't want thos to be marked as 'inside', otherwise invalid faces will be added // if
    // a face cuts a corner, nothings needs to be done, as that alone does not limit what adjacent faces will be
    // function handleEndPoint(a: IntersectionPointInfo, b: IntersectionPointInfo, newEdge: Edge) { // ends in the
    // middle of b's face if (a && !b) { if (!a.colinear && a.edgeT != a.edge.aT && a.edgeT != a.edge.bT) {
    // mapPush(thisEdgePoints, a.edge.getCanon(), a) assert(a.edge.isValidT(a.edgeT)) } // else colinear segment ends
    // in middle of other face, do nothing } // ends in the middle of a's face if (b && !a) { if (!b.colinear &&
    // b.edgeT != b.edge.aT && b.edgeT != b.edge.bT) { mapPush(otherEdgePoints, b.edge.getCanon(), b)
    // assert(b.edge.isValidT(b.edgeT)) } // else colinear segment ends in middle of other face, do nothing } if (a &&
    // b) { // if a or b is colinear the correct points will already have been added to the edge by handleNewEdge //
    // segment starts/ends on edge/edge intersection function foo(a, b, face, face2, thisPlane, face2Plane, thisBrep,
    // face2Brep, first, thisEdgePoints) { if (!a.colinear && a.edgeT != a.edge.aT && a.edgeT != a.edge.bT) { if
    // (!hasPair(a.edge.getCanon(), b.edge.getCanon())) { addPair(a.edge.getCanon(), b.edge.getCanon()) // ends on a,
    // on colinear segment b bT != a.edge.bT && // b can be colinear, so edgeT == aT is possible if (a.p.like(b.edge.a)
    // || a.p.like(b.edge.b)) { const corner = a.p.like(b.edge.a) ? b.edge.a : b.edge.b // face2brep corner on edge
    // const sVEC1 = splitsVolumeEnclosingCone(face2Brep, corner, a.edge.aDir) const sVEC2 =
    // splitsVolumeEnclosingCone(face2Brep, corner, a.edge.aDir.negated()) // if either of these return
    // ALONG_EDGE_OR_PLANE, then the breps share a colinear edge  if (INSIDE == sVEC1 || INSIDE == sVEC2) {
    // mapPush(thisEdgePoints, a.edge.getCanon(), a) assert(a.edge.isValidT(a.edgeT)) } } else { // edge / edge center
    // intersection const aEdgeDir = a.edge.tangentAt(a.edgeT) const bEdgeDir = b.edge.tangentAt(b.edgeT) const
    // testVector = aEdgeDir.rejectedFrom(bEdgeDir) assert(!testVector.likeO()) const sVEF1 =
    // splitsVolumeEnclosingFaces(face2Brep, b.edge.getCanon(), testVector, thisPlane.normal1) const sVEF2 =
    // splitsVolumeEnclosingFaces(face2Brep, b.edge.getCanon(), testVector.negated(), thisPlane.normal1) if (INSIDE ==
    // sVEF1 || INSIDE == sVEF2) { mapPush(thisEdgePoints, a.edge.getCanon(), a) assert(a.edge.isValidT(a.edgeT)) } } }
    // } }  foo(a, b, face, face2, thisPlane, face2Plane, thisBrep, face2Brep, true, thisEdgePoints) foo(b, a, face2,
    // face, face2Plane, thisPlane, face2Brep, thisBrep, false, otherEdgePoints)  } }   assertInst(PlaneFace, face2)
    // const face: PlaneFace = this // get intersection const thisPlane = this.surface.plane, face2Plane =
    // face2.surface.plane if (thisPlane.isParallelToPlane(face2Plane)) { if (thisPlane.like(face2Plane)) { // normal1
    // same and same location in space // addLikeSurfaceFaces(likeSurfaceFaces, this, face2) } return } const isLine =
    // L3.fromPlanes(thisPlane, face2Plane) // get intersections of newCurve with other edges of face and face2 const
    // ps1 = planeFaceEdgeISPsWithPlane(face, isLine, face2Plane) const ps2 = planeFaceEdgeISPsWithPlane(face2, isLine,
    // thisPlane) if (ps1.length == 0 || ps2.length == 0) { // faces to not intersect return }  let col1:
    // IntersectionPointInfo, col2: IntersectionPointInfo let in1 = false, in2 = false let i = 0, j = 0, last let
    // startP, startDir, startT, startA, startB while (i < ps1.length || j < ps2.length) { assert(i <= ps1.length)
    // assert(j <= ps2.length) const a = ps1[i], b = ps2[j] assert(a || b) if (j == ps2.length || i < ps1.length &&
    // lt(a.t, b.t)) { last = a in1 = !in1 a.used = true in1 && (col1 = a.colinear && a) i++ } else if (i == ps1.length
    // || gt(a.t, b.t)) { last = b in2 = !in2 b.used = true in2 && (col2 = b.colinear && b) j++ } else { // TODO: this
    // will break if 3 points on the same t last = a in1 = !in1 in2 = !in2 //if (in1 == in2) { a.used = true b.used =
    // true in1 && (col1 = a.colinear && a) in2 && (col2 = b.colinear && b) //} i++ j++ } if (startP && !(in1 && in2))
    // { // segment end const newEdge = new StraightEdge(isLine, startP, last.p, startT, last.t, undefined, 'genseg' +
    // getGlobalId()) startP = undefined last.used = true if (handleNewEdge(newEdge, col1 && col1.edge, col2 &&
    // col2.edge)) { handleEndPoint(startA || col1, startB || col2, newEdge) handleEndPoint(a && a.used && a || col1, b
    // && b.used && b || col2, newEdge) } } else if (in1 && in2) { // new segment just started startP = last.p startDir
    // = last.insideDir startT = last.t startA = a && a.used && a startB = b && b.used && b } if (!in1 && a && last ==
    // a && a.colinear) { checkedPairs.add(new Pair(a.edge.getCanon(), face2)) } if (!in2 && b && (last == b || b.used)
    // && b.colinear) { checkedPairs.add(new Pair(b.edge.getCanon(), face)) } } }
    withHole(holeEdges) {
        return new PlaneFace(this.surface, this.contour, [holeEdges]);
    }
    pointsToInside(p, dir) {
        return this.containsPoint2(p.plus(dir.times(ts3dutils.NLA_PRECISION * 8)));
    }
    edgeISPsWithPlane(isLine, plane2) {
        const face = this;
        ts3dutils.assert(face.surface.plane.containsLine(isLine));
        ts3dutils.assert(plane2.containsLine(isLine));
        const plane = face.surface.plane;
        const ps = [];
        const loops = [face.contour].concat(face.holes);
        loops.forEach(loop => {
            const colinearEdges = loop.map((edge) => edge.colinearToLine(isLine) && -sign$9(edge.aDir.dot(isLine.dir1)));
            const isLineOut = isLine.dir1.cross(plane.normal1);
            loop.forEach((edge, edgeIndex, edges) => {
                const nextEdgeIndex = (edgeIndex + 1) % edges.length, nextEdge = edges[nextEdgeIndex];
                //console.log(edge.toSource()) {p:V3(2, -2.102, 0),
                if (colinearEdges[edgeIndex]) {
                    // edge colinear to intersection line
                    const curveAT = isLine.pointT(edge.a), curveBT = isLine.pointT(edge.b);
                    // add interval for colinear segment
                    ps.push({ p: edge.a, insideDir: edge.aDir, t: curveAT, edge: edge, edgeT: edge.aT, colinear: true }, {
                        p: edge.b,
                        insideDir: edge.bDir.negated(),
                        t: curveBT,
                        edge: edge,
                        edgeT: edge.bT,
                        colinear: true,
                    });
                    // open next interval if necessary
                    const nextSide = colinearEdges[nextEdgeIndex] || dotCurve(isLineOut, nextEdge.aDir, nextEdge.aDDT);
                    if (colinearEdges[edgeIndex] * nextSide < 0) {
                        // side changes
                        ps.push({
                            p: nextEdge.a,
                            insideDir: edge.bDir,
                            t: curveBT,
                            edge: nextEdge,
                            edgeT: nextEdge.aT,
                            colinear: false,
                        });
                    }
                }
                else {
                    // not necessarily a straight edge, so multiple intersections are possible
                    const edgeTs = edge.edgeISTsWithPlane(plane2);
                    ts3dutils.assert(edgeTs.every(t => plane2.containsPoint(edge.curve.at(t))), edgeTs);
                    for (const edgeT of edgeTs) {
                        if (edgeT == edge.bT) {
                            // endpoint lies on intersection line
                            const side = -dotCurve(isLineOut, edge.bDir, edge.bDDT);
                            const nextSide = colinearEdges[nextEdgeIndex] || dotCurve(isLineOut, nextEdge.aDir, nextEdge.aDDT);
                            if (side * nextSide < 0) {
                                // next segment is not colinear and ends on different side
                                ps.push({
                                    p: edge.b,
                                    insideDir: plane2.normal1.negated(),
                                    t: isLine.pointT(edge.b),
                                    edge: edge,
                                    edgeT: edge.bT,
                                    colinear: false,
                                });
                            }
                        }
                        else if (edgeT != edge.aT) {
                            // edge crosses intersection line, neither starts nor ends on it
                            const p = edge.curve.at(edgeT);
                            ts3dutils.assert(plane2.containsPoint(p), edge.toString(), p, edgeT, plane2.distanceToPoint(p));
                            ts3dutils.assert(isLine.containsPoint(p), edge.toString(), p, edgeT, isLine.distanceToPoint(p));
                            const insideDir = plane2.normal1.negated();
                            ps.push({
                                p: p,
                                insideDir: insideDir,
                                t: isLine.pointT(p),
                                edge: edge,
                                edgeT: edgeT,
                                colinear: false,
                            });
                        }
                    }
                }
            });
        });
        // duplicate 't's are ok, as sometimes a segment needs to stop and start again
        // should be sorted so that back facing ones are first
        ps.sort((a, b) => a.t - b.t || a.insideDir.dot(isLine.dir1));
        return ps;
    }
}
class RotationFace extends Face {
    constructor(rot, contour, holes, name, info) {
        super(rot, contour, holes, name, info);
    }
    static loopDoesNotCrossPlane(loop, seamPlane) {
        let side = 0;
        // returns true if d is on the other side as previous calls
        function checkSide(d) {
            if (side == 0) {
                side = d;
            }
            else {
                return !side || side * d < 0;
            }
        }
        for (const edge of loop) {
            const ts = edge.edgeISTsWithPlane(seamPlane);
            if (ts.length == 0) {
                if (!(edge.curve instanceof L3$1) && checkSide(seamPlane.distanceToPointSigned(edge.a)))
                    return false;
            }
            else {
                for (const t of ts) {
                    // TODO: this part probably should be in a separate function
                    // check 'backwards' only if if aT != t
                    if (edge.aT != t) {
                        if (checkSide(dotCurve2(edge.curve, t, seamPlane.normal1, -sign$9(edge.bT - edge.aT))))
                            return false;
                    }
                    if (edge.bT != t) {
                        if (checkSide(dotCurve2(edge.curve, t, seamPlane.normal1, sign$9(edge.bT - edge.aT))))
                            return false;
                    }
                }
            }
        }
        return true;
    }
    getAABB() {
        if (this.aabb)
            return this.aabb;
        if (this.surface instanceof SemiEllipsoidSurface || this.surface instanceof EllipsoidSurface) {
            this.aabb = ts3dutils.AABB.forAABBs(this.contour.map(e => e.getAABB()));
            this.aabb.addPoints(this.surface.getExtremePoints().filter(p => this.containsPoint(p)));
            return this.aabb;
        }
        else {
            return super.getAABB();
        }
    }
    getCanonSeamU() {
        const stPFunc = this.surface.stPFunc();
        for (const edge of this.contour) {
            // check edge.a
            let u = stPFunc(edge.a, PI$15).x;
            // if u is not PI, or ~0, return its sign
            if (u != PI$15 && !ts3dutils.eq0(u)) {
                return sign$9(u) * PI$15;
            }
            // check midpoint between edge.a and edge.b
            u = stPFunc(edge.curve.at((edge.aT + edge.bT) / 2), PI$15).x;
            if (u != PI$15 && !ts3dutils.eq0(u)) {
                return sign$9(u) * PI$15;
            }
        }
        const localEdge = this.contour[0].transform(this.surface.inverseMatrix);
        if (P3.ZX.containsCurve(localEdge.curve)) {
            const insideVector = localEdge.a.cross(localEdge.aDir);
            return sign$9(insideVector.dot(ts3dutils.V3.Y)) * PI$15;
        }
        ts3dutils.assert(false, 'Couldn\'t find canon seam u');
    }
    unrollLoop(edgeLoop) {
        const vs = [];
        const reverseFunc = this.surface.stPFunc();
        const verticesNo0s = edgeLoop.map(edge => edge.getVerticesNo0());
        const startEdgeIndex = verticesNo0s.findIndex(edgeVertices => !ts3dutils.eq(reverseFunc(edgeVertices[0], Math.PI).x, Math.PI));
        ts3dutils.assert(-1 != startEdgeIndex);
        // console.log(startEdgeIndex)
        let hint = Math.PI;
        for (let i = 0; i < edgeLoop.length; i++) {
            const edgeIndex = (i + startEdgeIndex) % edgeLoop.length;
            for (let j = 0; j < verticesNo0s[edgeIndex].length; j++) {
                const p = verticesNo0s[edgeIndex][j];
                const localP = reverseFunc(p, hint);
                if (Math.abs(localP.x) < Math.PI - ts3dutils.NLA_PRECISION) {
                    // update hint
                    hint = localP.x;
                }
                // console.log(hint, p.sce, localP.sce)
                vs.push(localP);
            }
        }
        edgeLoop.forEach((edge, e) => {
            let hint = edge.bDir;
            if (edge instanceof StraightEdge && edge.curve.dir1.isParallelTo(this.surface.dir || this.surface.dir1)) {
                hint = this.surface.normalP(edge.b).cross(edge.bDir);
            }
            edge.getVerticesNo0().forEach(p => {
                vs.push(reverseFunc(p, hint));
            });
        });
        console.log('vs\n', vs.join('\n'), vs.length);
        return vs;
    }
    /**
     * f1 cos t + f2 sin t
     * tan(phi) = sin / cos
     *          = (f1x cos t + f2x sin t) / (f1y cos t + f2y sin t)
     *
     *          = (-f1x sin t + f2x cos t) / (-f1y sin t + f2y cos t)
     */
    unrollEllipsoidLoops(edgeLoops, uStep, vStep) {
        const verticesST = [], vertices = [], loopStarts = [];
        const ellipsoid = this.surface;
        const ptpf = ellipsoid.stPFunc();
        const testDegeneratePoint = ellipsoid instanceof SemiEllipsoidSurface
            ? (nextStart) => nextStart.like(ellipsoid.center.plus(ellipsoid.f3)) || nextStart.like(ellipsoid.center.minus(ellipsoid.f3))
            : (nextStart) => nextStart.like(this.surface.center);
        for (const edgeLoop of edgeLoops) {
            loopStarts.push(verticesST.length);
            // console.log(startEdgeIndex)
            const hint = this.getCanonSeamU();
            for (let i = 0; i < edgeLoop.length; i++) {
                const ipp = (i + 1) % edgeLoop.length;
                const verticesNo0 = edgeLoop[i].getVerticesNo0();
                vertices.push(...verticesNo0);
                verticesST.push(...verticesNo0.map(v => ptpf(v)));
                const nextStart = edgeLoop[ipp].a;
                //console.log('BLAH', nextStart.str, ellipsoid.center.plus(ellipsoid.f3).str)
                if (testDegeneratePoint(nextStart)) {
                    const bDirLC = ellipsoid.inverseMatrix.transformVector(edgeLoop[i].bDir), aDirLC = ellipsoid.inverseMatrix.transformVector(edgeLoop[ipp].aDir);
                    let inAngle = Math.atan2(-bDirLC.y, -bDirLC.x);
                    if (abs$12(inAngle) > Math.PI - ts3dutils.NLA_PRECISION) {
                        ts3dutils.assert(hint == -PI$15 || hint == PI$15);
                        inAngle = hint;
                    }
                    let outAngle = Math.atan2(aDirLC.y, aDirLC.x);
                    if (abs$12(outAngle) > Math.PI - ts3dutils.NLA_PRECISION) {
                        ts3dutils.assert(hint == -PI$15 || hint == PI$15);
                        outAngle = hint;
                    }
                    const stLast = verticesST.pop();
                    verticesST.push(new ts3dutils.V3(inAngle, stLast.y, 0), new ts3dutils.V3(outAngle, stLast.y, 0));
                    vertices.push(vertices.last);
                }
                verticesST.forEach(({ x: u, y: v }) => {
                    ts3dutils.assert(isFinite(u));
                    ts3dutils.assert(isFinite(v));
                });
            }
        }
        let normals;
        if (this.surface instanceof EllipsoidSurface) {
            normals = vertices.map(v => ellipsoid.normalP(v));
        }
        else {
            const pN = ellipsoid.normalSTFunc();
            normals = verticesST.map(({ x, y }) => pN(x, y));
        }
        ts3dutils.assert(vertices.length == vertices.length);
        //console.log(verticesST.map(v => v.str).join('\n'))
        return {
            verticesUV: verticesST.map(vST => new ts3dutils.V3(vST.x / uStep, vST.y / vStep, 0)),
            vertices: vertices,
            normals: normals,
            loopStarts: loopStarts,
        };
    }
    unrollCylinderLoops(loops, uStep, vStep) {
        const vertexLoops = loops.map(loop => loop.flatMap(edge => edge.getVerticesNo0()));
        const surface = this.surface;
        const vertices = vertexLoops.concatenated();
        // this.unrollLoop(loop).map(v => new V3(v.x / uStep, v.y / vStep, 0)))
        const loopStarts = vertexLoops.reduce((arr, loop) => (arr.push(arr.last + loop.length), arr), [0]);
        const stPFunc = surface.stPFunc();
        const verticesST = vertices.map(v => stPFunc(v));
        const verticesUV = verticesST.map(st => new ts3dutils.V3(st.x / uStep, st.y / vStep, 0));
        const normalST = surface.normalSTFunc();
        const normals = verticesST.map(({ x, y }) => normalST(x, y));
        return { verticesUV: verticesUV, vertices: vertices, normals: normals, loopStarts: loopStarts };
    }
    /**
     * at(s, t) = new V3(s cos t, s sin t, t + )
     *
     * x = 0
     *
     * s cos t = 0
     * ==> s = 0 || cos t = 0
     * ==> L3.Z || V3(0, +-s, k * 2 pi)
     *
     * x = c
     * s cos t = c
     * ==> V3(c, c sin t / cos t = c tan t, t)
     * ==> V3(c, c t, arctan t)
     *
     *
     * x . n = w
     *      s cos t nx + s sin t ny + t nz = w
     *      s = (w - t nz) / (cos t nx + sub t ny)
     * ==> V3(
     *          cos t (w - t nz) / (cos t nx + sin t ny)
     *          sin t (w - t nz) / (cos t nx + sin t ny)
     *          t)
     *
     *  ==> V3(
     *          (w - z arctan t) / (x + t y)
     *          (w - z arctan t) / (y + x / t)
     *          arctan t)
     *
     *
     *
     */
    addToMesh(mesh, uStep = this.surface.uStep, vStep = this.surface.vStep) {
        ts3dutils.assertf(() => uStep > 0 && vStep > 0, uStep, vStep, 'Surface: ' + this.surface);
        const triangles = [];
        const pIJFunc = (i, j) => this.surface.pSTFunc()(i * uStep, j * vStep);
        const normalIJFunc = (i, j) => this.surface.normalSTFunc()(i * uStep, j * vStep);
        const loops = [this.contour].concat(this.holes);
        const { vertices, verticesUV, normals, loopStarts } = this.surface instanceof SemiEllipsoidSurface || this.surface instanceof ConicSurface
            ? this.unrollEllipsoidLoops(loops, uStep, vStep)
            : this.unrollCylinderLoops(loops, uStep, vStep);
        loopStarts.push(vertices.length);
        for (let vertexLoopIndex = 0; vertexLoopIndex < loops.length; vertexLoopIndex++) {
            const vertexLoopStart = loopStarts[vertexLoopIndex];
            const vertexLoopLength = loopStarts[vertexLoopIndex + 1] - vertexLoopStart;
            const base = mesh.vertices.length + loopStarts[vertexLoopIndex];
            for (let i = 0; i < vertexLoopLength; i++) {
                mesh.LINES.push(base + i, base + (i + 1) % vertexLoopLength);
            }
        }
        ts3dutils.disableConsole();
        let minU = Infinity, maxU = -Infinity, minV = Infinity, maxV = -Infinity;
        //console.log('surface', this.surface.str)
        //console.log(verticesUV)
        //drPs.push(...verticesUV.map((v, i) => ({p: vertices[i], text: `${i} uv: ${v.toString(x => round10(x,
        // -4))}`})))
        verticesUV.forEach(({ x: u, y: v }) => {
            ts3dutils.assert(isFinite(u));
            ts3dutils.assert(isFinite(v));
            minU = min$10(minU, u);
            maxU = max$9(maxU, u);
            minV = min$10(minV, v);
            maxV = max$9(maxV, v);
        });
        if (ParametricSurface.is(this.surface)) {
            ts3dutils.assert(this.surface.boundsSigned(minU * uStep, minV * vStep) > -ts3dutils.NLA_PRECISION);
            ts3dutils.assert(this.surface.boundsSigned(maxU * uStep, maxV * vStep) > -ts3dutils.NLA_PRECISION);
        }
        const uOffset = floor$10(minU + ts3dutils.NLA_PRECISION), vOffset = floor$10(minV + ts3dutils.NLA_PRECISION);
        const uRes = ceil$12(maxU - ts3dutils.NLA_PRECISION) - uOffset, vRes = ceil$12(maxV - ts3dutils.NLA_PRECISION) - vOffset;
        console.log(uStep, vStep, uRes, vRes);
        if (uRes == 1 && vRes == 1) {
            // triangulate this face as if it were a plane
            const polyTriangles = triangulateVertices(ts3dutils.V3.Z, verticesUV, loopStarts.slice(1, 1 + this.holes.length));
            triangles.push(...polyTriangles);
        }
        else {
            const partss = new Array(uRes * vRes);
            function fixUpPart(part, baseU, baseV) {
                ts3dutils.assert(baseU < uRes && baseV < vRes, `${baseU}, ${baseV}, ${uRes}, ${vRes}`);
                console.log('complete part', part, baseU, baseV);
                //console.trace()
                ts3dutils.assert(part.length);
                const cellU = baseU + uOffset, cellV = baseV + vOffset;
                for (const index of part) {
                    ts3dutils.assert(ts3dutils.le(cellU, verticesUV[index].x) && ts3dutils.le(verticesUV[index].x, cellU + 1), `${index} ${verticesUV[index].str} ${cellU} ${cellU}`);
                    ts3dutils.assert(ts3dutils.le(cellV, verticesUV[index].y) && ts3dutils.le(verticesUV[index].y, cellV + 1));
                }
                const pos = baseV * uRes + baseU;
                (partss[pos] || (partss[pos] = [])).push(part);
                //const outline = partss[pos] || (partss[pos] = [minU + baseU * uStep, minV + baseV * vStep, minU +
                // (baseU + 1) * uStep, minV + (baseV + 1) * vStep])
            }
            // 'some' instead of forEach so we can return out of the entire function if this.edges crosses no borders
            // and
            for (let vertexLoopIndex = 0; vertexLoopIndex < loops.length; vertexLoopIndex++) {
                let part = undefined, firstPart, firstPartBaseU = -1, firstPartBaseV = -1;
                let lastBaseV = -1, lastBaseU = -1;
                let partCount = 0;
                const vertexLoopStart = loopStarts[vertexLoopIndex];
                const vertexLoopLength = loopStarts[vertexLoopIndex + 1] - vertexLoopStart;
                for (let vlvi = 0; vlvi < vertexLoopLength; vlvi++) {
                    const vx0index = vertexLoopStart + vlvi, vx0 = verticesUV[vx0index];
                    const vx1index = vertexLoopStart + (vlvi + 1) % vertexLoopLength, vx1 = verticesUV[vx1index];
                    //console.log('dask', vx0index, vx1index)
                    const vx01 = vx0.to(vx1);
                    ts3dutils.assert(vx0);
                    const di = vx01.x, dj = vx01.y;
                    let vxIndex = vx0index, vx = vx0, currentT = 0;
                    let whileLimit = 400;
                    while (--whileLimit) {
                        const vxu = vx.x, vxv = vx.y;
                        // points which are on a grid line are assigned to the cell into which they are going (+
                        // NLA_PRECISION * sign(di)) if they are parallel to the gridline (eq0(di)), they belong the
                        // the cell for which they are a CCW boundary
                        const baseU = floor$10(vxu + (!ts3dutils.eq0(di) ? sign$9(di) : -sign$9(dj)) * ts3dutils.NLA_PRECISION) - uOffset;
                        const baseV = floor$10(vxv + (!ts3dutils.eq0(dj) ? sign$9(dj) : sign$9(di)) * ts3dutils.NLA_PRECISION) - vOffset;
                        ts3dutils.assert(baseU < uRes && baseV < vRes, `${baseU}, ${baseV}, ${uRes}, ${vRes}`);
                        // figure out the next intersection with a gridline:
                        // iNext is the positive horizontal distance to the next vertical gridline
                        const iNext = ceil$12(sign$9(di) * vxu + ts3dutils.NLA_PRECISION) - sign$9(di) * vxu;
                        const jNext = ceil$12(sign$9(dj) * vxv + ts3dutils.NLA_PRECISION) - sign$9(dj) * vxv;
                        const iNextT = currentT + iNext / abs$12(di);
                        const jNextT = currentT + jNext / abs$12(dj);
                        //console.log(vxIndex, vx.str, 'vij', vxu, vxv, 'd', di, dj, 'ijNext', iNext, jNext, 'nextT',
                        // iNextT, jNextT)
                        if (lastBaseU != baseU || lastBaseV != baseV) {
                            if (part) {
                                if (!firstPart) {
                                    firstPart = part;
                                    firstPartBaseU = lastBaseU;
                                    firstPartBaseV = lastBaseV;
                                }
                                else {
                                    partCount++;
                                    fixUpPart(part, lastBaseU, lastBaseV);
                                }
                            }
                            part = [vxIndex];
                        }
                        lastBaseU = baseU;
                        lastBaseV = baseV;
                        currentT = min$10(iNextT, jNextT);
                        if (ts3dutils.ge(currentT, 1)) {
                            //console.log('breaking ', vx1index)
                            part.push(vx1index);
                            break;
                        }
                        else {
                            const nextPoint = vx0.lerp(vx1, currentT);
                            const nextPointIndex = addVertex(nextPoint.x, nextPoint.y);
                            //console.log('pushing ', nextPointIndex)
                            part.push(nextPointIndex);
                            vx = nextPoint;
                            vxIndex = nextPointIndex;
                        }
                    }
                    ts3dutils.assert(whileLimit, 'whileLimit');
                }
                if (0 == partCount) {
                    // complete loop
                    ts3dutils.assert(false, 'found a hole, try increasing resolution');
                }
                // at this point, the firstPart hasn't been added, and the last part also hasn't been added
                // either they belong to the same cell, or not
                if (firstPartBaseU == lastBaseU && firstPartBaseV == lastBaseV) {
                    part.pop();
                    fixUpPart(part.concat(firstPart), lastBaseU, lastBaseV);
                }
                else {
                    fixUpPart(firstPart, firstPartBaseU, firstPartBaseV);
                    fixUpPart(part, lastBaseU, lastBaseV);
                }
                console.log('firstPart', firstPart);
            }
            console.log('calculated parts', partss);
            const fieldVertexIndices = new Array((uRes + 1) * (vRes + 1));
            function addVertex(u, v) {
                verticesUV.push(new ts3dutils.V3(u, v, 0));
                normals.push(normalIJFunc(u, v));
                return vertices.push(pIJFunc(u, v)) - 1;
            }
            function getGridVertexIndex(i, j) {
                const index = j * (uRes + 1) + i;
                return fieldVertexIndices[index] || (fieldVertexIndices[index] = addVertex(i + uOffset, j + vOffset));
            }
            for (let col = 0; col < uRes; col++) {
                let inside = false;
                for (let row = 0; row < vRes; row++) {
                    const pos = row * uRes + col;
                    const fieldU = uOffset + col, fieldV = vOffset + row;
                    const parts = partss[pos];
                    if (!parts) {
                        if (inside) {
                            tsgl.pushQuad(triangles, false, getGridVertexIndex(col, row), getGridVertexIndex(col + 1, row), getGridVertexIndex(col, row + 1), getGridVertexIndex(col + 1, row + 1));
                        }
                    }
                    else {
                        // assemble the field with segments in in
                        function opos(index) {
                            const p = verticesUV[index], u1 = p.x - fieldU, v1 = p.y - fieldV;
                            ts3dutils.assert(-ts3dutils.NLA_PRECISION < u1 && u1 < 1 + ts3dutils.NLA_PRECISION && -ts3dutils.NLA_PRECISION < v1 && v1 < 1 + ts3dutils.NLA_PRECISION, 'oob u1 v1 ' + u1 + ' ' + v1 + ' ' + index + ' ' + p.str + 'IF THIS FAILS check canonSeamU is correct');
                            return v1 < u1 ? u1 + v1 : 4 - u1 - v1;
                        }
                        while (parts.length) {
                            const outline = [];
                            const startPart = parts[0];
                            ts3dutils.assert(startPart.length > 0);
                            let currentPart = startPart;
                            do {
                                outline.push(...currentPart);
                                const currentPartEndOpos = opos(currentPart.last);
                                const nextPartIndex = parts.indexWithMax(part => -ts3dutils.mod(opos(part[0]) - currentPartEndOpos, 4));
                                const nextPart = parts.removeIndex(nextPartIndex);
                                let currentOpos = currentPartEndOpos;
                                const nextPartStartOpos = opos(nextPart[0]) > currentOpos
                                    ? opos(nextPart[0])
                                    : opos(nextPart[0]) + 4;
                                let nextOpos = ceil$12(currentOpos + ts3dutils.NLA_PRECISION);
                                let flipping = ts3dutils.eq0((currentOpos + ts3dutils.NLA_PRECISION) % 1 - ts3dutils.NLA_PRECISION);
                                //inside = inside != (!eq0(currentOpos % 1) && currentOpos % 2 < 1)
                                while (ts3dutils.lt(nextOpos, nextPartStartOpos)) {
                                    switch (nextOpos % 4) {
                                        case 0:
                                            outline.push(getGridVertexIndex(col, row));
                                            break;
                                        case 1:
                                            inside = inside != flipping;
                                            outline.push(getGridVertexIndex(col + 1, row));
                                            break;
                                        case 2:
                                            outline.push(getGridVertexIndex(col + 1, row + 1));
                                            break;
                                        case 3:
                                            inside = inside != flipping;
                                            outline.push(getGridVertexIndex(col, row + 1));
                                            break;
                                    }
                                    flipping = true;
                                    nextOpos++;
                                }
                                // if the next loop would have completed a top or bottom segment
                                inside = inside != (flipping && nextOpos % 2 == 1 && ts3dutils.eq(nextOpos, nextPartStartOpos));
                                currentOpos = nextOpos;
                                currentPart = nextPart;
                            } while (currentPart != startPart);
                            // triangulate outline
                            if (outline.length == 3) {
                                // its just a triangle
                                triangles.push(...outline);
                            }
                            else {
                                const polyTriangles = triangulateVertices(ts3dutils.V3.Z, outline.map(i => verticesUV[i]), []).map(i => outline[i]);
                                triangles.push(...polyTriangles);
                            }
                            //console.log('outline', col, row, outline)
                        }
                    }
                }
            }
        }
        //console.log('trinagle', triangles.max(), vertices.length, triangles.length, triangles.toSource(),
        // triangles.map(col => vertices[col].$).toSource() ) assert(normals.every(n => n.hasLength(1)), normals.find(n
        // => !n.hasLength(1)).length() +' '+normals.findIndex(n => !n.hasLength(1)))
        Array.prototype.push.apply(mesh.TRIANGLES, triangles.map(index => index + mesh.vertices.length));
        Array.prototype.push.apply(mesh.vertices, vertices);
        Array.prototype.push.apply(mesh.normals, normals);
        //this.addEdgeLines(mesh)
        ts3dutils.enableConsole();
    }
    addToMesh2(mesh) {
        const closed = false;
        const hSplit = 12800, zSplit = 8;
        const ribs = [];
        let minZ = Infinity, maxZ = -Infinity;
        //let cmp = (a, b) => a.value - b.value
        const f = this.surface.pSTFunc();
        const normalF = this.surface.normalSTFunc();
        const vertexLoops = this.holes.concat([this.contour]).map(loop => this.unrollLoop(loop));
        vertexLoops.forEach(vertexLoop => {
            vertexLoop.forEach(({ x: d, y: z }) => {
                const index0 = ribs.binaryIndexOf(d, (a, b) => ts3dutils.snap(a.value - b, 0));
                if (index0 < 0) {
                    ribs.splice(-index0 - 1, 0, { value: d, left: [], right: [] });
                }
                minZ = min$10(minZ, z);
                maxZ = max$9(maxZ, z);
            });
        });
        console.log('zzzs', minZ, maxZ, vertexLoops[0].toSource().replace(/\), /g, ',\n'));
        const correction = 1;
        vertexLoops.forEach(vertexLoop => {
            vertexLoop.forEach((v0, i, vs) => {
                let v1 = vs[(i + 1) % vs.length], dDiff = v1.x - v0.x;
                //console.log(v0.sce, v1.sce)
                if (ts3dutils.eq0(dDiff)) {
                    return;
                }
                if (dDiff < 0) {
                    [v0, v1] = [v1, v0];
                    dDiff = -dDiff;
                }
                const index0 = ribs.binaryIndexOf(v0.x, (a, b) => ts3dutils.snap(a.value - b, 0));
                const index1 = ribs.binaryIndexOf(v1.x, (a, b) => ts3dutils.snap(a.value - b, 0));
                ribs[index0].right.binaryInsert(v0.y);
                for (let j = (index0 + correction) % ribs.length; j != index1; j = (j + correction) % ribs.length) {
                    const x = ribs[j].value;
                    const part = (x - v0.x) / dDiff;
                    const interpolated = v1.y * part + v0.y * (1 - part);
                    ribs[j].left.binaryInsert(interpolated);
                    ribs[j].right.binaryInsert(interpolated);
                }
                ribs[index1].left.binaryInsert(v1.y);
                // console.log(ribs.map(r=>r.toSource()).join('\n'))
            });
        });
        const vertices = [], triangles0 = [], normals = [];
        for (let i = 0; i < ribs.length; i++) {
            const ribLeft = ribs[i], ribRight = ribs[(i + 1) % ribs.length];
            ts3dutils.assert(ribLeft.right.length == ribRight.left.length);
            for (let j = 0; j < ribLeft.right.length; j++) {
                vertices.push(f(ribLeft.value, ribLeft.right[j]), f(ribRight.value, ribRight.left[j]));
                normals.push(normalF(ribLeft.value, ribLeft.right[j]), normalF(ribRight.value, ribRight.left[j]));
            }
        }
        //console.log(ribs.map(r=>r.toSource()).join('\n'))
        const vss = vertices.length, detailVerticesStart = vss;
        const zInterval = maxZ - minZ, zStep = zInterval / zSplit;
        const detailZs = ts3dutils.arrayFromFunction(zSplit - 1, i => minZ + (1 + i) * zStep);
        console.log('detailsZs', detailZs);
        for (let i = 0; i < ribs.length; i++) {
            const d = ribs[i].value;
            for (let j = 0; j < detailZs.length; j++) {
                vertices.push(f(d, detailZs[j]));
                normals.push(normalF(d, detailZs[j]));
            }
        }
        // console.log('detailVerticesStart', detailVerticesStart, 'vl', vertices.length, vertices.length -
        // detailVerticesStart, ribs.length) finally, fill in the ribs
        let vsStart = 0;
        const flipped2 = true;
        //for (var i = 0; i < 1; i++) {
        const end = closed ? ribs.length : ribs.length - 1;
        for (let i = 0; i < end; i++) {
            const ipp = (i + 1) % ribs.length;
            let inside = false, colPos = 0;
            const ribLeft = ribs[i], ribRight = ribs[(i + 1) % ribs.length];
            for (let j = 0; j < detailZs.length + 1; j++) {
                const detailZ = detailZs[j] || 100000;
                if (!inside) {
                    if (ribLeft.right[colPos] < detailZ && ribRight.left[colPos] < detailZ) {
                        if (ribLeft.right[colPos + 1] < detailZ || ribRight.left[colPos + 1] < detailZ) {
                            tsgl.pushQuad(triangles0, flipped2, vsStart + colPos * 2, vsStart + (colPos + 1) * 2, vsStart + colPos * 2 + 1, vsStart + (colPos + 1) * 2 + 1);
                            colPos += 2;
                            if (ribLeft.right[colPos] < detailZ || ribRight.left[colPos] < detailZ) {
                                j--;
                            }
                        }
                        else {
                            tsgl.pushQuad(triangles0, flipped2, vsStart + colPos * 2, vsStart + colPos * 2 + 1, detailVerticesStart + i * detailZs.length + j, detailVerticesStart + ipp * detailZs.length + j);
                            inside = true;
                            colPos++;
                        }
                    }
                }
                else {
                    if (ribLeft.right[colPos] < detailZ || ribRight.left[colPos] < detailZ) {
                        tsgl.pushQuad(triangles0, flipped2, detailVerticesStart + i * detailZs.length + j - 1, detailVerticesStart + ipp * detailZs.length + j - 1, vsStart + colPos * 2, vsStart + colPos * 2 + 1);
                        inside = false;
                        colPos++;
                        if (ribLeft.right[colPos] < detailZ || ribRight.left[colPos] < detailZ) {
                            j--;
                        }
                    }
                    else {
                        tsgl.pushQuad(triangles0, flipped2, detailVerticesStart + i * detailZs.length + j, detailVerticesStart + i * detailZs.length + j - 1, detailVerticesStart + ipp * detailZs.length + j, detailVerticesStart + ipp * detailZs.length + j - 1);
                    }
                }
            }
            vsStart += ribLeft.right.length * 2;
        }
        //console.log('trinagle', triangles0.max(), vertices.length, triangles0.length, triangles0.toSource(),
        // triangles0.map(i => vertices[i].$).toSource() )
        const triangles = triangles0.map(index => index + mesh.vertices.length);
        //assert(normals.every(n => n.hasLength(1)), normals.find(n => !n.hasLength(1)).length() +'
        // '+normals.findIndex(n => !n.hasLength(1)))
        Array.prototype.push.apply(mesh.vertices, vertices);
        Array.prototype.push.apply(mesh.TRIANGLES, triangles);
        Array.prototype.push.apply(mesh.normals, normals);
        //this.addEdgeLines(mesh)
    }
}

/// <reference path="earcut.d.ts" />
const { PI: PI$16, sign: sign$10, abs: abs$13, sqrt: sqrt$6 } = Math;
const EPS = 1e-5;
let globalId = 0;
function getGlobalId() {
    return globalId++;
}
function addLikeSurfaceFaces(likeSurfaceFaces, face1, face2) {
    // There cannot be two subgroups which will later be connected, as the "graph" of like surface faces is fully
    // connected
    for (let i = 0; i < likeSurfaceFaces.length; i++) {
        const faceGroup = likeSurfaceFaces[i];
        let foundFace1 = false, foundFace2 = false;
        for (let j = 0; j < faceGroup.length; j++) {
            const face = faceGroup[j];
            if (face == face1) {
                foundFace1 = true;
            }
            if (face == face2) {
                foundFace2 = true;
            }
        }
        if (foundFace1 != foundFace2) {
            faceGroup.push(foundFace1 ? face2 : face1);
            return;
        }
        else if (foundFace1) {
            // found both
            return;
        }
    }
    // nothing found, add a new group
    likeSurfaceFaces.push([face1, face2]);
}
function assembleFaceFromLooseEdges(edges, surface, faceConstructor) {
    const visited = new Set();
    function nextStart() { return edges.find(edge => !visited.has(edge)); }
    const loops = [];
    let startEdge, currentEdge;
    while (startEdge = nextStart()) {
        currentEdge = startEdge;
        const loop = [];
        let total = 0;
        do {
            visited.add(currentEdge);
            loop.push(currentEdge);
            const possibleEdges = edges.filter(edge => currentEdge.b.like(edge.a));
            const normalAtCurrentB = surface.normalP(currentEdge.b);
            const nextEdgeIndex = possibleEdges.indexWithMax((edge, index) => currentEdge.bDir.angleRelativeNormal(edge.aDir, normalAtCurrentB));
            currentEdge = possibleEdges[nextEdgeIndex];
        } while (startEdge != currentEdge && total++ < 200);
        ts3dutils.assert(total != 201);
        loops.push(loop);
    }
    const assembledFaces = B2.assembleFacesFromLoops(loops, surface, faceConstructor);
    ts3dutils.assertf(() => 1 == assembledFaces.length);
    return assembledFaces[0];
}
/**
 * ## Markdown header
 * ![foo](screenshots/Capture.PNG)
 * {@link ../screenshots/Capture.PNG}
 * find the next edge with the MAXIMUM angle
 */
function calcNextEdgeIndex(currentEdge, possibleEdges, faceNormalAtCurrentB) {
    let maxValue = -20, advanced = false, result = Number.MAX_SAFE_INTEGER;
    const normVector = currentEdge.bDir.cross(faceNormalAtCurrentB);
    const eps = 1e-4;
    const dir = sign$10(currentEdge.deltaT());
    const ecd = currentEdge.curve.diff(currentEdge.bT, -dir * eps).dot(normVector);
    for (let i = possibleEdges.length; i--;) {
        const edge = possibleEdges[i];
        const angle1 = currentEdge.bDir.negated().angleRelativeNormal(edge.aDir, faceNormalAtCurrentB);
        const angle = (angle1 + ts3dutils.TAU + ts3dutils.NLA_PRECISION) % ts3dutils.TAU - ts3dutils.NLA_PRECISION;
        if (ts3dutils.eq0(angle)) {
            // do advanced analysis
            if (currentEdge.curve.isColinearTo(edge.curve)) {
                continue;
            }
            const edgeDir = sign$10(edge.deltaT());
            const iscd = edge.curve.diff(edge.aT, edgeDir * eps).dot(normVector);
            const diff = (iscd - ecd);
            // if diff > 0, the angle is actually ~= 0
            if (diff < 0 && (!advanced || diff > maxValue)) {
                advanced = true;
                maxValue = diff;
                result = i;
            }
        }
        else if (!advanced) {
            if (ts3dutils.gt(angle, maxValue)) {
                maxValue = angle;
                result = i;
            }
        }
    }
    return result == Number.MAX_SAFE_INTEGER ? 0 : result;
}
class B2 extends ts3dutils.Transformable {
    constructor(faces, infiniteVolume, generator, vertexNames) {
        super();
        this.faces = faces;
        ts3dutils.assertInst(Face, ...faces);
        this.infiniteVolume = infiniteVolume;
        ts3dutils.assert(false === this.infiniteVolume || true === this.infiniteVolume);
        this.generator = generator;
        this.vertexNames = vertexNames;
        this.edgeFaces = undefined;
        //this.assertSanity()
    }
    static loop1ContainsLoop2(loop1, ccw1, loop2, ccw2, surface) {
        for (const edge of loop2) {
            const loop1ContainsPoint = surface.loopContainsPoint(loop1, edge.a);
            if (exports.PointVsFace.ON_EDGE != loop1ContainsPoint)
                return exports.PointVsFace.INSIDE == loop1ContainsPoint;
        }
        for (const edge of loop2) {
            const edgePoint = edge.curve.at(edge.aT * 0.2 + edge.bT * 0.8);
            const loop1ContainsPoint = surface.loopContainsPoint(loop1, edgePoint);
            if (exports.PointVsFace.ON_EDGE != loop1ContainsPoint)
                return exports.PointVsFace.INSIDE == loop1ContainsPoint;
        }
        if (ccw1 != ccw2) {
            return ccw2;
        }
        throw new Error(loop1.sce + loop2.sce);
    }
    static assembleFacesFromLoops(loops, surface, originalFace, infoFactory) {
        function placeRecursively(newLoopInfo, loopInfos) {
            if (loopInfos.length == 0) {
                loopInfos.push(newLoopInfo);
            }
            else {
                const subLoopInfo = loopInfos.find(loopInfo => B2.loop1ContainsLoop2(loopInfo.loop, loopInfo.ccw, newLoopInfo.loop, newLoopInfo.ccw, surface));
                if (subLoopInfo) {
                    placeRecursively(newLoopInfo, subLoopInfo.subloops);
                }
                else {
                    // newLoopInfo isnt contained by any other subLoopInfo
                    for (let i = loopInfos.length; --i >= 0;) {
                        const subLoopInfo = loopInfos[i];
                        //console.log("cheving subLoopInfo", surface.loopContainsPoint(newLoopInfo.edges,
                        // subLoopInfo.edges[0].a))
                        if (B2.loop1ContainsLoop2(newLoopInfo.loop, newLoopInfo.ccw, subLoopInfo.loop, subLoopInfo.ccw, surface)) {
                            newLoopInfo.subloops.push(subLoopInfo);
                            loopInfos.splice(i, 1); // remove it
                        }
                    }
                    loopInfos.push(newLoopInfo);
                }
            }
        }
        function newFacesRecursive(loopInfo) {
            // CW loops can be top level, if they are holes in the original face not contained in the new face
            if (loopInfo.ccw) {
                if (loopInfo.subloops.every(sl => !sl.ccw)) {
                    const holes = loopInfo.subloops.map(sl => sl.loop);
                    const info = infoFactory && infoFactory.newSubFace(originalFace, surface, loopInfo.loop, holes);
                    const newFace = new originalFace.constructor(surface, loopInfo.loop, holes, 'genface' + getGlobalId(), info);
                    newFaces.push(newFace);
                    loopInfo.subloops.forEach(sl => sl.subloops.forEach(slsl => slsl.ccw && newFacesRecursive(slsl)));
                }
                else {
                    loopInfo.subloops.forEach(sl => sl.ccw && newFacesRecursive(sl));
                }
            }
        }
        const newFaces = [];
        const topLevelLoops = [];
        loops.forEach(loop => placeRecursively({
            loop: loop,
            ccw: surface.edgeLoopCCW(loop),
            subloops: [],
        }, topLevelLoops));
        topLevelLoops.forEach(tll => newFacesRecursive(tll));
        return newFaces;
    }
    static join(b2s, generator) {
        return new B2(b2s.flatMap(b2 => b2.faces), false, generator);
    }
    containsPoint(p, forceInsideOutside = false) {
        const dirs = [ts3dutils.V(-0.3920414696448526, -0.12936136783391444, -0.9108068525164064), ts3dutils.V(0.6520650903544943, -0.07151288645511984, -0.7547827667692488), ts3dutils.V(0.9433494201061395, -0.2402757256238473, -0.22882186797013926), ts3dutils.V(0.13678704228501923, -0.04480387361087783, 0.9895867410047372), ts3dutils.V(0.0662057922721913, -0.5865836917435423, 0.8071780259955845), ts3dutils.V(-0.7322576567870621, -0.12953393611526787, 0.6685953061989045), ts3dutils.V(0.6579719127258273, -0.012300218400456116, 0.7529420075219719), ts3dutils.V(-0.5576497966736425, 0.8006695748324647, 0.2189861552871446)];
        dirLoop: for (const dir of dirs) {
            const testLine = new L3$1(p, dir);
            let inside = this.infiniteVolume;
            for (const face of this.faces) {
                ts3dutils.assert(!face.surface.containsCurve(testLine));
                const ists = face.surface.isTsForLine(testLine);
                for (const t of ists) {
                    const p = testLine.at(t);
                    const pvf = face.containsPoint2(p);
                    //assert(pvf != PointVsFace.ON_EDGE)
                    !forceInsideOutside && ts3dutils.assert(!ts3dutils.eq0(t));
                    if (t > 0) {
                        if (pvf == exports.PointVsFace.ON_EDGE) {
                            continue dirLoop;
                        }
                        if (pvf == exports.PointVsFace.INSIDE) {
                            inside = !inside;
                            
                        }
                    }
                }
            }
            return inside;
        }
        return false;
    }
    withMergedFaces() {
        const likeSurfaceFaces = [];
        for (let i = 0; i < this.faces.length; i++) {
            let addedToGroup = false;
            for (let j = 0; j < i; j++) {
                if (this.faces[i].surface.isCoplanarTo(this.faces[j].surface)) {
                    const faceGroup = likeSurfaceFaces.find(faceGroup => faceGroup.includes(this.faces[j]));
                    if (faceGroup) {
                        faceGroup.push(this.faces[i]);
                        addedToGroup = true;
                    }
                }
            }
            !addedToGroup && likeSurfaceFaces.push([this.faces[i]]);
        }
        console.log('likeSurfaceFaces', likeSurfaceFaces);
        if (likeSurfaceFaces.every(group => group.length == 1))
            return this;
        const newFaces = [];
        let total = 0;
        for (const faceGroup of likeSurfaceFaces) {
            console.log(faceGroup);
            if (faceGroup.length == 1) {
                newFaces.push(faceGroup[0]);
            }
            else {
                const allEdges = faceGroup.flatMap(face => face.getAllEdges());
                for (let i = allEdges.length; i-- > 0;) {
                    for (let j = 0; j < i; j++) {
                        console.log('blugh', total);
                        ts3dutils.assert(i >= 0 && j >= 0 && total++ < 500, i + ' ' + j + ' ' + total);
                        if (allEdges[i].isCoEdge(allEdges[j])) {
                            // remove both
                            allEdges.splice(i, 1);
                            allEdges.splice(j, 1);
                            i--;
                            break;
                        }
                    }
                }
                const newFace = assembleFaceFromLooseEdges(allEdges, faceGroup[0].surface, faceGroup[0].constructor);
                newFaces.push(newFace);
            }
        }
        return new B2(newFaces, this.infiniteVolume, this.generator && this.generator + '.withMergedFaces()', this.vertexNames);
    }
    calculateVolume() {
        return this.faces.map(face => face.zDirVolume().volume).sum();
    }
    toMesh() {
        const mesh = new tsgl.Mesh()
            .addVertexBuffer('normals', 'LGL_Normal')
            .addIndexBuffer('TRIANGLES')
            .addIndexBuffer('LINES');
        mesh.faceIndexes = new Map();
        for (const face of this.faces) {
            const triangleStart = mesh.TRIANGLES.length;
            face.addToMesh(mesh);
            mesh.faceIndexes.set(face, { start: triangleStart, count: mesh.TRIANGLES.length - triangleStart });
        }
        //this.buildAdjacencies()
        //for (const edge of this.edgeFaces.keys()) {
        //
        //}
        return mesh;
    }
    minus(other, infoFactory) {
        const generator = this.generator && other.generator && this.generator + '.minus(' + other.generator + ')';
        return this.intersection(other.flipped(), true, true, generator, infoFactory);
    }
    plus(other, infoFactory) {
        const generator = this.generator && other.generator && ts3dutils.callsce(this.generator + '.plus', other.generator);
        return this.flipped().intersection(other.flipped(), true, true, generator, infoFactory).flipped();
    }
    and(other, infoFactory) {
        const generator = this.generator && other.generator && ts3dutils.callsce(this.generator + '.and', other.generator);
        return this.intersection(other, true, true, generator, infoFactory);
    }
    xor(other, infoFactory) {
        const s = this.generator && other.generator && ts3dutils.callsce(this.generator + '.xor', other.generator);
        return new B2(this.minus(other).faces.concat(other.minus(this).faces), this.infiniteVolume != other.infiniteVolume, s);
    }
    equals(obj) {
        return this.faces.length == obj.faces.length &&
            this.faces.every((face) => obj.faces.some((face2) => face.equals(face2)));
    }
    like(brep) {
        return this.faces.length == brep.faces.length &&
            this.faces.every((face) => brep.faces.some((face2) => face.likeFace(face2)));
    }
    //reconstituteCoplanarFaces(likeSurfacePlanes, edgeLooseSegments, faceMap, newFaces) {
    //    likeSurfacePlanes.forEach(faceGroup => {
    //        // calculate total contours
    //        let surface = faceGroup[0].surface, bag = []
    //        faceGroup.forEach(face => {
    //            Array.prototype.push.apply(bag, faceMap(face))
    //            face.getAllEdges().forEach(edge => {
    //                let edgeSubSegments
    //                if (edgeSubSegments = edgeLooseSegments.get(edge)) {
    //                    Array.prototype.push.apply(bag, edgeSubSegments)
    //                } else {
    //                    bag.push(edge)
    //                }
    //            })
    //        })
    //        let currentEdge, loops = []
    //        while (currentEdge = bag.find(edge => !edge.visited)) {
    //            let path = []
    //            do {
    //                currentEdge.visited = true
    //                path.push(currentEdge)
    //                let possibleNextEdges = bag.filter(edge => currentEdge.b.like(edge.a))
    //                // lowest angle, i.e. the right-most next edge
    //                let nextEdgeIndex = possibleNextEdges.indexWithMax((edge, index) =>
    // -currentEdge.bDir.angleRelativeNormal(edge.aDir, surface.normalP(currentEdge.b))) currentEdge =
    // possibleNextEdges[nextEdgeIndex] } while (!currentEdge.visited) let startIndex = path.find(currentEdge) if (-1
    // != startIndex) { loops.push(path.slice(startIndex)) } } }) }
    toString() {
        return `new B2([\n${this.faces.join(',\n').replace(/^/gm, '\t')}], ${this.infiniteVolume})`;
    }
    getConstructorParameters() {
        return [this.faces, this.infiniteVolume];
    }
    toSource(useGenerator = true) {
        return useGenerator && this.generator ||
            `new B2([\n${this.faces.map(ts3dutils.SCE).join(',\n').replace(/^/gm, '\t')}], ${this.infiniteVolume})`;
    }
    /**
     * Rightmost next segment doesn't work, as the correct next segment isn't obvious from the current corner
     * alone.
     * (at least, not without extensive pre-analysis on the face edges, which shouldn't be necessary, as the
     * correct new faces are defined by the new edges already.) Leftmost edge should work. Holes which touch the
     * edge of the face will be added to the face contour.
     *
     * New segments will always be part left-er than existing ones, so no special check is required.
     *
     */
    reconstituteFaces(oldFaces, edgeSubEdges, faceMap, newFaces, infoFactory) {
        const oldFaceStatuses = new Map();
        // reconstitute faces
        const insideEdges = [];
        for (const face of oldFaces) {
            const usableOldEdges = face.getAllEdges().filter(edge => !edgeSubEdges.get(edge));
            const subEdges = face.getAllEdges().mapFilter(edge => edgeSubEdges.get(edge)).concatenated();
            const newEdges = faceMap.get(face) || [];
            if (newEdges.length || subEdges.length) {
                oldFaceStatuses.set(face, 'partial');
                const loops = [];
                // new edges are definitely part of a resulting loop
                // old edges (both contour and holes) can either be part of a new loop, in which case they will already
                // have been visited when starting a loop search with a new edge, OR they can be stranded, OR they can
                // remain in their old loop
                function getNextStart() {
                    return newEdges.find(edge => !visitedEdges.has(edge))
                        || subEdges.find(edge => !visitedEdges.has(edge))
                        || usableOldEdges.find(edge => !visitedEdges.has(edge));
                }
                const visitedEdges = new Set();
                // search for a loop:
                let currentEdge;
                while (currentEdge = getNextStart()) {
                    const startEdge = currentEdge, edges = [];
                    let i = 0;
                    // wether only new edges are used (can include looseSegments)
                    do {
                        visitedEdges.add(currentEdge);
                        edges.push(currentEdge);
                        // find next edge
                        const possibleOldEdges = usableOldEdges.filter(edge => currentEdge.b.like(edge.a));
                        const possibleSubEdges = subEdges.filter(edge => currentEdge.b.like(edge.a));
                        const possibleNewEdges = newEdges.filter(edge => currentEdge.b.like(edge.a));
                        const possibleEdges = possibleOldEdges.concat(possibleSubEdges, possibleNewEdges);
                        if (0 == possibleEdges.length)
                            break;
                        ts3dutils.assert(0 < possibleEdges.length, () => face.sce);
                        const faceNormalAtCurrentB = face.surface.normalP(currentEdge.b);
                        const correct = possibleEdges.indexWithMax((edge, index) => (currentEdge.bDir.angleRelativeNormal(edge.aDir, faceNormalAtCurrentB) + ts3dutils.NLA_PRECISION + PI$16) % ts3dutils.TAU);
                        const nextEdgeIndex = calcNextEdgeIndex(currentEdge, possibleEdges, faceNormalAtCurrentB);
                        currentEdge = possibleEdges[nextEdgeIndex];
                        if (visitedEdges.has(currentEdge)) {
                            break;
                        }
                        ts3dutils.assert(currentEdge);
                        ts3dutils.assert(currentEdge != startEdge);
                    } while (++i < 400);
                    if (400 == i) {
                        ts3dutils.assert(false, 'too many');
                    }
                    // check if we found a loop
                    if (edges.length > 1 && currentEdge == startEdge) {
                        loops.push(edges);
                    }
                }
                const faceNewFaces = B2.assembleFacesFromLoops(loops, face.surface, face, infoFactory);
                newFaces.push(...faceNewFaces);
                const faceNewFacesEdges = faceNewFaces.flatMap(face => face.getAllEdges());
                insideEdges.push(...usableOldEdges.filter(edge => faceNewFacesEdges.includes(edge)));
            }
        }
        while (insideEdges.length != 0) {
            const insideEdge = insideEdges.pop();
            const adjacentFaces = this.edgeFaces.get(insideEdge.getCanon());
            adjacentFaces.forEach(info => {
                if (!oldFaceStatuses.has(info.face)) {
                    oldFaceStatuses.set(info.face, 'inside');
                    insideEdges.push.apply(insideEdges, info.face.getAllEdges());
                }
            });
        }
        newFaces.push(...oldFaces.filter(face => oldFaceStatuses.get(face) == 'inside'));
    }
    getLooseEdgeSegments(edgePointInfoss, edgeFaces) {
        const result = new javasetmap_ts.JavaMap();
        // if there are no point info, the original edge will be kept, so we should return nothing
        // otherwise, something will be returned, even if it a new edge identical to the base edge
        for (const [canonEdge, pointInfos] of edgePointInfoss) {
            if (0 == pointInfos.length)
                continue;
            const allFaces = edgeFaces.get(canonEdge);
            pointInfos.sort((a, b) => ts3dutils.snap0(a.edgeT - b.edgeT) || +!!a.faces);
            let startP = canonEdge.a, startDir = canonEdge.aDir, startT = canonEdge.aT, startInfo;
            function addNewEdge(startInfo, endInfo, newEdge) {
                for (let i = 0; i < allFaces.length; i++) {
                    const faceInfo = allFaces[i];
                    const startYes = !startInfo || !startInfo.faces || startInfo.faces[i];
                    const endYes = !endInfo || !endInfo.faces;
                    endYes && ts3dutils.mapPush(result, !faceInfo.reversed ? canonEdge : canonEdge.flipped(), !faceInfo.reversed ? newEdge : newEdge.flipped());
                }
            }
            for (let i = 0; i < pointInfos.length; i++) {
                const info = pointInfos[i];
                const pDir = canonEdge.tangentAt(info.edgeT);
                if (!ts3dutils.eq(info.edgeT, startT)) {
                    const newEdge = Edge.create(canonEdge.curve, startP, info.p, startT, info.edgeT, undefined, startDir, pDir, 'looseSegment' + getGlobalId());
                    addNewEdge(startInfo, info, newEdge);
                }
                startP = info.p;
                startT = info.edgeT;
                startInfo = info;
                startDir = pDir;
            }
            if (startInfo && !ts3dutils.eq(startT, canonEdge.bT)) {
                const newEdge = Edge.create(canonEdge.curve, startP, canonEdge.b, startT, canonEdge.bT, undefined, startDir, canonEdge.bDir, 'looseSegment' + getGlobalId());
                addNewEdge(startInfo, undefined, newEdge);
            }
        }
        return result;
    }
    getIntersectionEdges(brep2) {
        const faceMap = new Map(), thisEdgePoints = new javasetmap_ts.JavaMap(), otherEdgePoints = new javasetmap_ts.JavaMap();
        const likeSurfaceFaces = [];
        this.faces.forEach(face => {
            //console.log('face', face.toString())
            brep2.faces.forEach(face2 => {
                //console.log('face2', face2.toString())
                face.intersectFace(face2, this, brep2, faceMap, thisEdgePoints, otherEdgePoints, likeSurfaceFaces);
            });
        });
        return Array.from(faceMap.values()).concatenated();
    }
    shellCount() {
        const foundFaces = new Set();
        let face, result = 0;
        while (face = this.faces.find(face => !foundFaces.has(face))) {
            result++;
            const stack = [face];
            while (face = stack.pop()) {
                for (const edge of face.getAllEdges()) {
                    for (const { face: face2 } of this.edgeFaces.get(edge.getCanon())) {
                        if (face !== face2 && !foundFaces.has(face2)) {
                            foundFaces.add(face2);
                            stack.push(face2);
                        }
                    }
                }
            }
        }
        return result;
    }
    getAABB() {
        return ts3dutils.AABB.forAABBs(this.faces.map(face => face.getAABB()));
    }
    assertSanity() {
        if (!ts3dutils.NLA_DEBUG)
            return;
        const allFaceEdges = this.faces.flatMap(face => face.getAllEdges());
        for (const { i, j } of ts3dutils.combinations(allFaceEdges.length)) {
            const a = allFaceEdges[i], b = allFaceEdges[j];
            //assert(i == j || !a.isCoEdge(b) || a == b || a.flippedOf == b, 'coedges not linked properly', a, b)
            //assert(i == j
            //	|| !a.curve.isColinearTo(b.curve)
            //	|| (a.curve.equals(b.curve) && a.isCoEdge(b))
            //   || !a.overlaps(b), 'colinear edges overlap', a, b)
        }
        this.buildAdjacencies();
        for (const [canonEdge, edgeFaceInfos] of this.edgeFaces) {
            // TODO handle curved faces
            ts3dutils.assert(edgeFaceInfos.length % 2 == 0, () => canonEdge + edgeFaceInfos.sce);
        }
    }
    //intersection3(other: B2, buildThis: boolean, buildOther: boolean, name?: string): B2 {
    //    this.assertSanity()
    //    other.assertSanity()
    //    this.buildAdjacencies()
    //    other.buildAdjacencies()
    //
    //    // edge / edge
    //    for (const [edge1, edge1Faces] of this.edgeFaces) {
    //        for (const [edge2, edge2Faces] of other.edgeFaces) {
    //            const curve1 = edge1.curve, curve2 = edge2.curve
    //            if (curve1.isColinearTo(curve2)) {
    //                if (edge1.overlaps(edge2)) {
    //                    // faces have a common edge
    //                    const aT = curve1.pointT(edge2.a), bT = curve1.pointT(edge2.a)
    //                    const minT = min(aT, bT), maxT = max(aT, bT)
    //                    const commonEdge = Edge.create(curve1, min(edge1.minT, minT), min(edge1.maxT, maxT), )
    //                }
    //            } else if (x = curve1.isInfosWithCurve(edge2.curve)) {
    //                // edges intersect in a point
    //            }
    //        }
    //    }
    //
    //    // point / edge
    //    function pointEdge(b1, b2, has, add) {
    //        for (const v1 of this.vertFaces.keys()) {
    //            for (const edge2 of other.edgeFaces.keys()) {
    //                if (edge2.curve.containsPoint(v1)) {
    //                    const edge2T = edge2.curve.pointT(v1)
    //                    if (eq(edge2.aT, edge2T) || eq(edge2.bT, edge2T)) {
    //                        add(v1, eq(edge2.aT, edge2T) ? edge2.a : edge2.b)
    //                    }
    //                }
    //            }
    //        }
    //    }
    //    const pairs: CustomSet<[Equalable, Equalable]> = new CustomSet<[Equalable, Equalable]>()
    //    pointEdge(this, other, (a, b) => pairs.has([a, b]), (a, b) => pairs.add([a, b]))
    //    pointEdge(other, this, (b, a) => pairs.has([a, b]), (b, a) => pairs.add([a, b]))
    //
    //
    //    // point / point
    //    for (const v1 of this.vertFaces.keys()) {
    //        for (const v2 of other.vertFaces.keys()) {
    //            if (v1.like(v2)) {
    //
    //            }
    //        }
    //    }
    //
    //    for (const face1 of this.faces) {
    //        for (const face2 of other.faces) {
    //            face1.intersectFace(face2)
    //        }
    //    }
    //
    //}
    buildAdjacencies() {
        if (this.edgeFaces)
            return this;
        this.edgeFaces = new javasetmap_ts.JavaMap();
        for (const face of this.faces) {
            for (const edge of face.getAllEdges()) {
                const canon = edge.getCanon();
                const normalAtCanonA = face.surface.normalP(canon.a);
                const inside = normalAtCanonA.cross(canon == edge ? edge.aDir : edge.bDir);
                ts3dutils.mapPush(this.edgeFaces, canon, {
                    face: face,
                    edge: edge,
                    normalAtCanonA: normalAtCanonA,
                    reversed: canon != edge,
                    inside: inside,
                    angle: 0,
                });
            }
        }
        for (const [canonEdge, edgeFaceInfos] of this.edgeFaces) {
            // TODO handle curved faces
            //assert(edgeFaceInfos.length % 2 == 0, () => canonEdge + edgeFaceInfos.sce)
            const faceInfo0 = edgeFaceInfos.find(faceInfo => faceInfo.reversed);
            if (!faceInfo0) {
                console.warn('invalid brep');
                continue;
            }
            edgeFaceInfos.forEach(faceInfo => {
                if (faceInfo != faceInfo0) {
                    faceInfo.angle = faceInfo0.inside.angleRelativeNormal(faceInfo.inside, canonEdge.aDir.unit());
                    if (faceInfo.angle < 0)
                        faceInfo.angle += 2 * Math.PI;
                }
            });
            edgeFaceInfos.sort((a, b) => ts3dutils.snap(a.angle - b.angle, 0)); // TODO  || assertNever()
        }
        return this;
    }
    /**
     * Cases for volumes A and B
     *
     *          1.  Volumes do not touch.
     *          2.  face/face Face surfaces intersect each other.
     *              implies edges going through faces.
     *              e.g. box(5, 5, 5) - box(5, 5, 5).translate(1, 1, 1)
     *          3.  face/edge Edge of A lies in a face of B
     *              implies vertices of A lying in face of B
     *              e.g. box(5, 5, 5) - box(3, 3, 3).rotateZ([0, 1, 2] * PI / 2).translate(0, 1, 1)
     *          4.  edge/edge Two edges are colinear.
     *              implies vertex of A lying in edge of B
     *           5.  vertex/edge Vertex of A lies on edge of B (but no edge/edge)
     *          6.  vertex/vertex with/without edge/edge, edge/face and face/face intersections
     *          7.  vertex lies in face
     *
     *
     *
     */
    intersection(other, buildThis, buildOther, generator, infoFactory) {
        this.assertSanity();
        other.assertSanity();
        this.buildAdjacencies();
        other.buildAdjacencies();
        const faceMap = new Map();
        const thisEdgePoints = new javasetmap_ts.JavaMap(), otherEdgePoints = new javasetmap_ts.JavaMap();
        const checkedPairs = new javasetmap_ts.JavaSet();
        for (const thisFace of this.faces) {
            for (const otherFace of other.faces) {
                thisFace.intersectFace(otherFace, this, other, faceMap, thisEdgePoints, otherEdgePoints, checkedPairs);
            }
        }
        for (const edge of thisEdgePoints.keys()) {
            ts3dutils.assert(this.edgeFaces.get(edge));
        }
        for (const edge of otherEdgePoints.keys()) {
            ts3dutils.assert(other.edgeFaces.get(edge));
        }
        const newFaces = [];
        if (0 == faceMap.size && 0 == thisEdgePoints.size && 0 == otherEdgePoints.size) {
            const thisInOther = other.containsPoint(this.faces[0].contour[0].a, true);
            const otherInThis = !thisInOther && this.containsPoint(other.faces[0].contour[0].a);
            return this;
        }
        else {
            if (buildThis) {
                const edgeLooseSegments = this.getLooseEdgeSegments(thisEdgePoints, this.edgeFaces);
                //noinspection JSUnusedLocalSymbols
                const els = this.faces.map(face => [face,
                    Array.from(edgeLooseSegments.entries()).filter(([edge, subs]) => face.getAllEdges().some(e => e.equals(edge))).concatenated()]);
                this.reconstituteFaces(this.faces, edgeLooseSegments, faceMap, newFaces, infoFactory);
            }
            if (buildOther) {
                const edgeLooseSegments = this.getLooseEdgeSegments(otherEdgePoints, other.edgeFaces);
                //noinspection JSUnusedLocalSymbols
                const els = other.faces.map(face => [face,
                    Array.from(edgeLooseSegments.entries())
                        .filter(([edge, subs]) => face.getAllEdges().some(e => e.equals(edge)))
                        .flatMap(([edge, subs]) => subs)]);
                other.reconstituteFaces(other.faces, edgeLooseSegments, faceMap, newFaces, infoFactory);
            }
        }
        //buildCoplanar && this.reconstituteCoplanarFaces(likeSurfaceFaces, edgeLooseSegments, faceMap, newFaces,
        // this.infiniteVolume, other.infiniteVolume)
        const result = new B2(newFaces, this.infiniteVolume && other.infiniteVolume, generator);
        //result.buildAdjacencies()
        return result;
    }
    transform(m4, desc) {
        let vertexNames;
        if (this.vertexNames) {
            vertexNames = new Map();
            this.vertexNames.forEach((name, vertex) => vertexNames.set(m4.transformPoint(vertex), name + desc));
        }
        return new B2(this.faces.map(f => f.transform(m4)), this.infiniteVolume, this.generator && desc && this.generator + desc, // if desc isn't set, the generator will be invalid
        vertexNames);
    }
    flipped() {
        return new B2(this.faces.map(f => f.flipped()), !this.infiniteVolume, this.generator && this.generator + '.flipped()', this.vertexNames);
    }
}
B2.EMPTY = new B2([], false, 'B2.EMPTY', new Map()).buildAdjacencies();
B2.R3 = new B2([], true, 'B2.R3', new Map()).buildAdjacencies();
function dotCurve(v, cDir, cDDT) {
    let dot = v.dot(cDir);
    if (ts3dutils.eq0(dot)) {
        dot = v.dot(cDDT);
    }
    ts3dutils.assert(!ts3dutils.eq0(dot));
    return dot;
}
function dotCurve2(curve, t, normal, sign) {
    ts3dutils.assert(sign == 1 || sign == -1, sign);
    const tangentDot = curve.tangentAt(t).dot(normal);
    // if tangentDot != 0 the curve simply crosses the plane
    if (!ts3dutils.eq0(tangentDot)) {
        return sign * tangentDot;
    }
    const ddtDot = curve.ddt(t).dot(normal);
    // tangentDot == 0 ==> critical point at t, if ddtDot != 0, then it is a turning point, otherwise we can't be sure
    // and must do a numeric test
    if (!ts3dutils.eq0(ddtDot)) {
        return ddtDot;
    }
    const numericDot = curve.at(t).to(curve.at(t + sign * 4 * ts3dutils.NLA_PRECISION)).dot(normal);
    ts3dutils.assert(!(curve instanceof L3$1));
    return numericDot;
}
const INSIDE = 0;
const OUTSIDE$1 = 1;
const COPLANAR_SAME = 2;
const COPLANAR_OPPOSITE = 3;
const ALONG_EDGE_OR_PLANE = 4;
/**
 *
 * @param brep BREP to check
 * @param edge edge to check
 * @param dirAtEdgeA the direction vector to check
 * @param faceNormal If dirAtEdgeA doesn't split a volume, but is along a face, the returned value depends on
 *     wether that faces normal1 points in the same direction as faceNormal
 * @returns INSIDE, OUTSIDE, COPLANAR_SAME or COPLANAR_OPPOSITE
 */
//function splitsVolumeEnclosingFaces(brep: B2, edge: Edge, dirAtEdgeA: V3, faceNormal: V3): int {
//    assert(arguments.length == 4)
//    //assert(p.equals(edge.a))
//    const ab1 = edge.aDir.unit()
//    const relFaces = facesWithEdge(edge, brep.faces) as any[]
//    relFaces.forEach(faceInfo => {
//        faceInfo.normalAtEdgeA = faceInfo.face.surface.normalP(edge.a)
//        faceInfo.edgeDirAtEdgeA = !faceInfo.reversed
//            ? faceInfo.edge.aDir
//            : faceInfo.edge.bDir
//        faceInfo.outsideVector = faceInfo.edgeDirAtEdgeA.cross(faceInfo.normalAtEdgeA)
//        faceInfo.angle = (dirAtEdgeA.angleRelativeNormal(faceInfo.outsideVector.negated(), ab1) + 2 * Math.PI +
// NLA_PRECISION / 2) % (2 * Math.PI) }) assert(relFaces.length != 0, edge.toSource()) relFaces.sort((a, b) => a.angle
// - b.angle) // assert(relFaces.length % 2 == 0, edge.toSource()) // even number of touching faces  if
// (eq0(relFaces[0].angle)) { //assert(false) todo const coplanarSame = relFaces[0].normalAtEdgeA.dot(faceNormal) > 0;
// return coplanarSame ? COPLANAR_SAME : COPLANAR_OPPOSITE } else { return !relFaces[0].reversed ? INSIDE : OUTSIDE } }
function splitsVolumeEnclosingFaces(brep, canonEdge, dirAtEdgeA, faceNormal) {
    ts3dutils.assert(arguments.length == 4);
    ts3dutils.assert(canonEdge == canonEdge.getCanon());
    //assert(p.equals(canonEdge.a))
    const edgeFaceInfos = brep.edgeFaces.get(canonEdge);
    ts3dutils.assertf(() => edgeFaceInfos.length % 2 == 0);
    ts3dutils.assertf(() => brep.edgeFaces);
    const faceInfo0 = edgeFaceInfos[0];
    const aDir1 = canonEdge.aDir.unit();
    const angleToCanon = (faceInfo0.inside.angleRelativeNormal(dirAtEdgeA, aDir1) + 2 * Math.PI + ts3dutils.NLA_PRECISION) % (2 * Math.PI) - ts3dutils.NLA_PRECISION;
    const nearestFaceInfoIndex = edgeFaceInfos.findIndex(faceInfo => ts3dutils.lt(angleToCanon, faceInfo.angle));
    const nearestFaceInfo = edgeFaceInfos[nearestFaceInfoIndex == -1
        ? edgeFaceInfos.length - 1
        : nearestFaceInfoIndex - 1];
    if (ts3dutils.eq(nearestFaceInfo.angle, angleToCanon)) {
        //assert(false) todo
        const coplanarSame = nearestFaceInfo.normalAtCanonA.dot(faceNormal) > 0;
        return coplanarSame ? COPLANAR_SAME : COPLANAR_OPPOSITE;
    }
    else {
        return nearestFaceInfo.reversed ? INSIDE : OUTSIDE$1;
    }
}
function splitsVolumeEnclosingFacesP(brep, canonEdge, p, pInside, faceNormal) {
    ts3dutils.assert(arguments.length == 5);
    ts3dutils.assert(canonEdge == canonEdge.getCanon());
    //assert(p.equals(canonEdge.a))
    ts3dutils.assertf(() => brep.edgeFaces);
    const edgeFaceInfos = brep.edgeFaces.get(canonEdge);
    ts3dutils.assertf(() => edgeFaceInfos.length % 2 == 0);
    const pDir1 = canonEdge.tangentAt(canonEdge.curve.pointT(p)).unit();
    const faceInfoAngleFromPInsideNeg = faceInfo => {
        const faceInfoPDir = faceInfo.edge.getCanon() == faceInfo.edge ? pDir1 : pDir1.negated();
        const faceInfoInsideAtP = faceInfo.face.surface.normalP(p).cross(faceInfoPDir);
        const faceInfoAngleAtP = pInside.angleRelativeNormal(faceInfoInsideAtP, pDir1);
        return -((faceInfoAngleAtP + ts3dutils.TAU + ts3dutils.NLA_PRECISION) % ts3dutils.TAU - ts3dutils.NLA_PRECISION);
    };
    const nearestFaceInfo = edgeFaceInfos.withMax(faceInfoAngleFromPInsideNeg);
    if (ts3dutils.eq0(faceInfoAngleFromPInsideNeg(nearestFaceInfo))) {
        //assert(false) todo
        const coplanarSame = nearestFaceInfo.normalAtCanonA.dot(faceNormal) > 0;
        return coplanarSame ? COPLANAR_SAME : COPLANAR_OPPOSITE;
    }
    else {
        return nearestFaceInfo.reversed ? OUTSIDE$1 : INSIDE;
    }
}
function splitsVolumeEnclosingFacesP2(brep, canonEdge, p, testCurve, curveT, dir, faceNormal) {
    ts3dutils.assert(canonEdge == canonEdge.getCanon());
    //assert(p.equals(canonEdge.a))
    ts3dutils.assertf(() => brep.edgeFaces);
    const edgeFaceInfos = brep.edgeFaces.get(canonEdge);
    ts3dutils.assertf(() => edgeFaceInfos.length % 2 == 0);
    const pDir1 = canonEdge.tangentAt(canonEdge.curve.pointT(p)).unit();
    let pInside = testCurve.tangentAt(curveT).times(dir);
    if (pInside.isParallelTo(pDir1)) {
        pInside = testCurve.diff(curveT, 1e-4 * dir / testCurve.tangentAt(curveT).length()).rejectedFrom(pDir1);
        pInside = pInside.div(pInside.length());
    }
    let minValue = 20, advanced = false, result = OUTSIDE$1;
    for (const faceInfo of edgeFaceInfos) {
        const faceInfoPDir = faceInfo.edge.getCanon() == faceInfo.edge ? pDir1 : pDir1.negated();
        const faceInfoInsideAtP = faceInfo.face.surface.normalP(p).cross(faceInfoPDir);
        const faceInfoAngleAtP = pInside.angleRelativeNormal(faceInfoInsideAtP, pDir1);
        const angle = (faceInfoAngleAtP + ts3dutils.TAU + ts3dutils.NLA_PRECISION) % ts3dutils.TAU - ts3dutils.NLA_PRECISION;
        if (ts3dutils.eq0(angle)) {
            // do advanced analysis
            const normVector = faceInfo.face.surface.normalP(p);
            if (faceInfo.face.surface.containsCurve(testCurve)) {
                const coplanarSame = normVector.dot(faceNormal) > 0;
                return coplanarSame ? COPLANAR_SAME : COPLANAR_OPPOSITE;
            }
            const testPlane = P3.normalOnAnchor(pDir1, p);
            const isCurve = faceInfo.face.surface.isCurvesWithPlane(testPlane)[0];
            const isCurvePT = isCurve.pointT(p);
            const dirFactor = sign$10(isCurve.tangentAt(isCurvePT).dot(pInside));
            const eps = 1e-4;
            const iscd = isCurve.at(isCurvePT).to(isCurve.at(isCurvePT + dir * dirFactor * eps)).dot(normVector);
            const ecd = testCurve.at(curveT).to(testCurve.at(curveT + dir * eps)).dot(normVector);
            const diff = (iscd - ecd) * (faceInfo.reversed ? -1 : 1);
            if (diff > 0 && (!advanced || diff < minValue)) {
                advanced = true;
                minValue = diff;
                result = faceInfo.reversed ? OUTSIDE$1 : INSIDE;
            }
        }
        else if (!advanced) {
            if (angle < minValue) {
                minValue = angle;
                result = faceInfo.reversed ? OUTSIDE$1 : INSIDE;
            }
        }
    }
    return result;
}
function splitsVolumeEnclosingCone(brep, p, dir) {
    const testPlane = P3.forAnchorAndPlaneVectors(p, dir, dir.getPerpendicular());
    const rays = [];
    for (let k = 0; k < brep.faces.length; k++) {
        const planeFace = brep.faces[k];
        ts3dutils.assertf(() => planeFace instanceof PlaneFace);
        if (planeFace.getAllEdges().some(edge => edge.a.like(p))) {
            if (testPlane.isParallelToPlane(planeFace.surface.plane)) {
                if (planeFace.pointsToInside(p, dir) != exports.PointVsFace.OUTSIDE) {
                    return ALONG_EDGE_OR_PLANE;
                }
            }
            else {
                const isLine = L3$1.fromPlanes(testPlane, planeFace.surface.plane);
                const ps = planeFace.edgeISPsWithPlane(isLine, testPlane);
                let i = 0;
                while (i < ps.length) {
                    const a = ps[i++], b = ps[i++];
                    const out = a.p.like(p);
                    if (out || b.p.like(p)) {
                        const dir2 = out ? isLine.dir1 : isLine.dir1.negated();
                        const angle = (dir.angleRelativeNormal(dir2, testPlane.normal1) + 2 * Math.PI + ts3dutils.NLA_PRECISION / 2) % (2 * Math.PI);
                        rays.push({ angle: angle, out: out });
                    }
                }
            }
        }
    }
    rays.sort((a, b) => a.angle - b.angle);
    //console.log("testPlane", testPlane.toSource(), "rays", rays.toSource())
    if (ts3dutils.eq0(rays[0].angle)) {
        return ALONG_EDGE_OR_PLANE;
    }
    else {
        return rays[0].out ? OUTSIDE$1 : INSIDE;
    }
}
function splitsVolumeEnclosingCone2(brep, p, curve, curveT, fb) {
    ts3dutils.assert(curve.containsPoint(p));
    const dir = curve.tangentAt(curveT).times(fb);
    const testPlane = P3.forAnchorAndPlaneVectors(p, dir, dir.getPerpendicular());
    const pFaces = brep.faces.filter(face => face.getAllEdges().some(edge => edge.a.like(p)));
    for (let k = 0; k < pFaces.length; k++) {
        const face = pFaces[k];
        if (face.surface.containsCurve(curve)) {
            //assert(false)
            if (face.pointsToInside3(p, curve, curveT, fb) != exports.PointVsFace.OUTSIDE) {
                return ALONG_EDGE_OR_PLANE;
            }
        }
    }
    const EPS = 1e-6;
    return brep.containsPoint(curve.at(curveT + fb * EPS), true) ? INSIDE : OUTSIDE$1;
}
function fff(info, surface) {
    const canonA = info.edge.reversed ? info.edge.b : info.edge.a;
    const surfaceNormalAtCanonA = surface.normalP(canonA);
    const dot = ts3dutils.snap0(info.inside.dot(surfaceNormalAtCanonA));
    if (0 !== dot) {
        return 0 < dot ? OUTSIDE$1 : INSIDE;
    }
    if (surface.isCoplanarTo(info.face.surface)) {
        return 0 < info.normalAtCanonA.dot(surfaceNormalAtCanonA) ? COPLANAR_SAME : COPLANAR_OPPOSITE;
    }
    ts3dutils.assert(false);
}
function triangulateVertices(normal, vertices, holeStarts) {
    const absMaxDim = normal.maxAbsDim(), factor = sign$10(normal.e(absMaxDim));
    const contour = new Float64Array(vertices.length * 2);
    let i = vertices.length;
    /*
     var [coord0, coord1] = [['y', 'z'], ['z', 'x'], ['x', 'y']][maxAbsDim]
     while (i--) {
     contour[i * 2    ] = vertices[i][coord0] * factor
     contour[i * 2 + 1] = vertices[i][coord1]
     }
     */
    while (i--) {
        // unroll disambiguation instead of accessing elements by string name ([coord0] etc)
        // as it confuses google closure
        switch (absMaxDim) {
            case 0:
                contour[i * 2] = vertices[i].y * factor;
                contour[i * 2 + 1] = vertices[i].z;
                break;
            case 1:
                contour[i * 2] = vertices[i].z * factor;
                contour[i * 2 + 1] = vertices[i].x;
                break;
            case 2:
                contour[i * 2] = vertices[i].x * factor;
                contour[i * 2 + 1] = vertices[i].y;
                break;
        }
    }
    return earcut(contour, holeStarts);
}
/**
 * Solves a quadratic system of equations of the form
 *      a * x + b * y = c
 *      x² + y² = 1
 * This can be understood as the intersection of the unit circle with a line.
 *      => y = (c - a x) / b
 *      => x² + (c - a x)² / b² = 1
 *      => x² b² + c² - 2 c a x + a² x² = b²
 *      => (a² + b²) x² - 2 a c x + (c² - b²) = 0
 *
 * a * b + (b -c) * (b + c)
 */
function intersectionUnitCircleLine(a, b, c) {
    ts3dutils.assertNumbers(a, b, c);
    // TODO: disambiguate on a < b
    const term = sqrt$6(a * a + b * b - c * c);
    return {
        x1: (a * c + b * term) / (a * a + b * b),
        x2: (a * c - b * term) / (a * a + b * b),
        y1: (b * c - a * term) / (a * a + b * b),
        y2: (b * c + a * term) / (a * a + b * b),
    };
}
function intersectionUnitCircleLine2(a, b, c) {
    ts3dutils.assertNumbers(a, b, c);
    // TODO: disambiguate on a < b
    // cf. pqFormula
    const termSqr = ts3dutils.snap0(a * a + b * b - c * c);
    if (termSqr < 0) {
        return [];
    }
    else if (termSqr == 0) {
        return [[(a * c) / (a * a + b * b),
                (b * c) / (a * a + b * b)]];
    }
    else {
        const term = sqrt$6(termSqr);
        return [[(a * c + b * term) / (a * a + b * b),
                (b * c - a * term) / (a * a + b * b)],
            [(a * c - b * term) / (a * a + b * b),
                (b * c + a * term) / (a * a + b * b)]];
    }
}
function intersectionCircleLine(a, b, c, r) {
    ts3dutils.assertNumbers(a, b, c, r);
    const term = sqrt$6(r * r * (a * a + b * b) - c * c);
    return {
        x1: (a * c + b * term) / (a * a + b * b),
        x2: (a * c - b * term) / (a * a + b * b),
        y1: (b * c - a * term) / (a * a + b * b),
        y2: (b * c + a * term) / (a * a + b * b),
    };
}
/**
 * Solves a quadratic system of equations of the form
 *      a * x + b * y = c
 *      x^2 - y^2 = 1
 * This can be understood as the intersection of the unit hyperbola with a line.
 *
 * @returns with x1 >= x2 and y1 <= y2
 * a * b + (b -c) * (b + c)
 */
function intersectionUnitHyperbolaLine(a, b, c) {
    ts3dutils.assertNumbers(a, b, c);
    const aa = a * a, bb = b * b, cc = c * c;
    // TODO: disambiguate on a < b
    //var xTerm = sqrt(4*cc*aa-4*(bb-aa)*(-cc-bb))
    const xTerm = 2 * sqrt$6(bb * cc + bb * bb - aa * bb);
    const yTerm = sqrt$6(4 * cc * bb - 4 * (bb - aa) * (cc - aa));
    return {
        x1: (-2 * a * c + xTerm) / 2 / (bb - aa),
        x2: (-2 * a * c - xTerm) / 2 / (bb - aa),
        y1: (2 * b * c - yTerm) / 2 / (bb - aa),
        y2: (2 * b * c + yTerm) / 2 / (bb - aa),
    };
}
function followAlgorithm2d(ic, startP, stepLength = 0.5, bounds, endP = startP, startTangent) {
    ts3dutils.assertNumbers(stepLength, ic(0, 0));
    ts3dutils.assertVectors(startP);
    if (!startTangent) {
        startTangent = new ts3dutils.V3(-ic.y(startP.x, startP.y), ic.x(startP.x, startP.y), 0).toLength(stepLength);
    }
    ts3dutils.assertVectors(startTangent);
    const points = [];
    const tangents = [];
    ts3dutils.assert(ts3dutils.eq0(ic(startP.x, startP.y), 0.01), 'isZero(implicitCurve(startPoint.x, startPoint.y))');
    let i = 0, p = startP, tangent = startTangent;
    do {
        points.push(p);
        tangents.push(tangent);
        const searchStart = p.plus(tangent);
        ts3dutils.assert(searchStart);
        const newP = curvePointMF(ic, searchStart);
        const dfpdx = ic.x(newP.x, newP.y), dfpdy = ic.y(newP.x, newP.y);
        const newTangent = new ts3dutils.V3(-dfpdy, dfpdx, 0).toLength(stepLength);
        //const reversedDir = p.minus(prevp).dot(tangent) < 0
        if (p.equals(newP)) {
            ts3dutils.assertNever();
        }
        // check if we passed a singularity
        if (tangent.dot(newTangent) < 0) {
            const singularity = ts3dutils.newtonIterate2d(ic.x, ic.y, p.x, p.y);
            if (ts3dutils.eq0(ic(singularity.x, singularity.y)) && singularity.distanceTo(p) < abs$13(stepLength)) {
                // end on this point
                points.push(singularity);
                tangents.push(p.to(singularity));
                break;
            }
            else {
                throw new Error();
            }
        }
        if (i > 4) {
            if (!bounds(p.x, p.y)) {
                break;
            }
            // full loop or arrived at end
            if (p.distanceTo(endP) < stepLength) {
                points.push(endP);
                const endTangent = new ts3dutils.V3(-ic.y(endP.x, endP.y), ic.x(endP.x, endP.y), 0).toLength(stepLength);
                tangents.push(endTangent);
                break;
            }
        }
        ts3dutils.assert(ts3dutils.eq0(ic(newP.x, newP.y), ts3dutils.NLA_PRECISION * 2), p, newP, searchStart);
        tangent = newTangent;
        p = newP;
    } while (++i < 1000);
    ts3dutils.assert(i < 1000);
    //assert(points.length > 6)
    return { points, tangents };
}
function followAlgorithm2dAdjustable(ic, start, stepLength = 0.5, bounds, endp = start) {
    ts3dutils.assertNumbers(stepLength, ic(0, 0));
    ts3dutils.assertVectors(start);
    //assert (!startDir || startDir instanceof V3)
    const points = [];
    const tangents = [];
    ts3dutils.assert(ts3dutils.eq0(ic(start.x, start.y), 0.01), 'isZero(implicitCurve(startPoint.x, startPoint.y))');
    let p = start, prevp = p;
    let i = 0;
    do {
        const dfpdx = ic.x(p.x, p.y), dfpdy = ic.y(p.x, p.y);
        const dfpdxx = ic.xx(p.x, p.y), dfpdyy = ic.yy(p.x, p.y), dfpdxy = ic.xy(p.x, p.y);
        const c2factor = abs$13((Math.pow(dfpdy, 2) * dfpdxx - 2 * dfpdx * dfpdy * dfpdxy + Math.pow(dfpdx, 2) * dfpdyy) /
            Math.pow((Math.pow(dfpdx, 2) + Math.pow(dfpdy, 2)), 2));
        const c2 = new ts3dutils.V3(dfpdx, dfpdy, 0).times(c2factor);
        const s = 1 / 16 / c2.length();
        const tangent = new ts3dutils.V3(-dfpdy, dfpdx, 0).unit();
        const reversedDir = p.minus(prevp).dot(tangent) < 0;
        const newPStart = p.plus(tangent.times(s).plus(c2.times(Math.pow(s, 2) / 2)));
        points.push(p);
        tangents.push(tangent);
        prevp = p;
        const newP = curvePointMF(ic, newPStart);
        if (newP.equals(p)) {
            ts3dutils.assertNever();
        }
        console.log(p.to(newP).length());
        p = newP;
        ts3dutils.assert(ts3dutils.eq0(ic(p.x, p.y)));
    } while (i++ < 1000 && (i < 4 || prevp.distanceTo(endp) > stepLength) && bounds(p.x, p.y));
    ts3dutils.assert(i != 1000);
    //assert(bounds(p.x, p.y))
    const end = (i < 4 || prevp.distanceTo(endp) > stepLength) ? p : endp;
    const endTangent = new ts3dutils.V3(-ic.y(end.x, end.y), ic.x(end.x, end.y), 0).toLength(stepLength);
    points.push(end);
    tangents.push(endTangent);
    //assert(points.length > 6)
    // TODO gleichmäßige Verteilung der Punkte
    return { points, tangents };
}
// both curves must be in the same s-t coordinates for this to make sense
function intersectionICurveICurve(iCurve1, startParams1, endParams1, startDir, stepLength, iCurve2) {
    ts3dutils.assertNumbers(stepLength, iCurve1(0, 0), iCurve2(0, 0));
    ts3dutils.assertVectors(startParams1, endParams1);
    ts3dutils.assert(!startDir || startDir instanceof ts3dutils.V3);
    const vertices = [];
    ts3dutils.assert(ts3dutils.eq0(iCurve1(startParams1.x, startParams1.y)));
    stepLength = stepLength || 0.5;
    const eps = 1e-5;
    let p = startParams1, prevp = p; // startDir ? p.minus(startDir) : p
    let i = 0;
    while (i++ < 1000 && (i < 4 || p.distanceTo(endParams1) > 1.1 * stepLength)) {
        const fp = iCurve1(p.x, p.y);
        const dfpdx = (iCurve1(p.x + eps, p.y) - fp) / eps, dfpdy = (iCurve1(p.x, p.y + eps) - fp) / eps;
        let tangent = new ts3dutils.V3(-dfpdy, dfpdx, 0).toLength(stepLength);
        if (p.minus(prevp).dot(tangent) < 0)
            tangent = tangent.negated();
        prevp = p;
        p = curvePoint(iCurve1, p.plus(tangent));
        vertices.push(p);
    }
    // TODO gleichmäßige Verteilung der Punkte
    return vertices;
}
function intersectionICurveICurve2(iCurve1, loopPoints1, iCurve2) {
    let p = loopPoints1[0], val = iCurve2(p.x, p.y), lastVal;
    const iss = [];
    for (let i = 0; i < loopPoints1.length; i++) {
        lastVal = val;
        p = loopPoints1[i];
        val = iCurve2(p);
        if (val * lastVal <= 0) {
            iss.push(ts3dutils.newtonIterate2d(iCurve1, iCurve2, p.x, p.y));
        }
    }
    return iss;
}
//export function intersectionPCurveISurface(parametricCurve: ParametricCurve, searchStart: number, searchEnd: number, searchStep: number, implicitSurface) {
//	assertNumbers(searchStart, searchEnd, searchStep)
//	const iss = []
//	let val = implicitSurface(parametricCurve(searchStart)), lastVal
//	for (let t = searchStart + searchStep; t <= searchEnd; t += searchStep) {
//		lastVal = val
//		val = implicitSurface(parametricCurve(t))
//		if (val * lastVal <= 0) {
//			iss.push(newtonIterate1d(t => implicitSurface(parametricCurve(t)), t))
//		}
//	}
//	return iss
//}
function intersectionICurvePSurface(f0, f1, parametricSurface) {
}
//
//function test2() {
//    const ic: R2_R = (x, y) => sin(x+y)-cos(x*y)+1
//    const dids: R2_R = (x, y) => y * sin(x * y) + cos(x + y)
//    const didt: R2_R = (x, y) => x * sin(x * y) + cos(x + y)
//    const ic2: R2_R = (x, y) => (3 * x ** 2 - y ** 2) ** 2 * y ** 2 - (x ** 2 + y ** 2) ** 4
//    const di2ds: R2_R = (x, y) => 4* x* (9* x**2* y**2 - 3* y**4 - 2* (x**2 + y**2)**3)
//    const di2dt: R2_R = (x, y) => 2 * y * (-4 * (x ** 2 + y ** 2) ** 3 + (3 * x ** 2 - y ** 2) ** 2 + 2 * y ** 2 * (y
// ** 2 - 3 * x ** 2)) const start = V(-3.6339970071165784, 3.5625834844534974, 0) // curvePoint(ic, V(-4, 4))
// assert(eq02(ic(start.x, start.y), 0.1)) const bounds = (s: number, t: number) => -5 <= s && s <= 5 && -5 <= t && t
// <= 5 //const curves =  Curve.breakDownIC(ic, -5, 5, -5, 5, 0.1, 0.1, 0.05, dids, didt) const curves =
// Curve.breakDownIC(ic2, {sMin: -5, sMax: 5, tMin: -5, tMax: 5}, 0.1, 0.1, 0.02, di2ds, di2dt) //const curves =
// Curve.breakDownIC(cassini(1, 1.02), -5, 5, -5, 5, 0.1, 0.1, 0.02) //const curves = mkcurves(ic, start.x, start.y,
// 0.05, dids, didt, bounds) .map(({points, tangents}, i) => { const curve = new ImplicitCurve(ic, points, tangents)
// return Edge.forCurveAndTs(curve.translate(5, 0, 0.1 * i)) }) //checkDerivate(s => ic(s, 0), s => dids(s, 0), -5, 5,
// 0) //checkDerivate(t => ic(0, t), t => dids(0, t), -5, 5, 0) console.log(curves.length) return curves  }
function cassini(a, c) {
    return (x, y) => (x * x + y * y) * (x * x + y * y) - 2 * c * c * (x * x - y * y) - (Math.pow(a, 4) - Math.pow(c, 4));
}

(function (MathFunctionR2R) {
    function forNerdamer(expression, args = ['x', 'y']) {
        const ndf = nerdamer(expression);
        const ndfs = nerdamer.diff(ndf, args[0]);
        const ndft = nerdamer.diff(ndf, args[1]);
        const f = ndf.buildFunction(args);
        f.x = ndfs.buildFunction(args);
        f.y = ndft.buildFunction(args);
        f.xx = nerdamer.diff(ndfs, args[0]).buildFunction(args);
        f.xy = nerdamer.diff(ndfs, args[1]).buildFunction(args);
        f.yy = nerdamer.diff(ndft, args[1]).buildFunction(args);
        return f;
    }
    MathFunctionR2R.forNerdamer = forNerdamer;
    function nerdamerToR2_R(expression, args = ['x', 'y']) {
        return expression.buildFunction(args);
    }
    MathFunctionR2R.nerdamerToR2_R = nerdamerToR2_R;
    function forFFxFy(f, fx, fy) {
        f.x = fx;
        f.y = fy;
        return f;
    }
    MathFunctionR2R.forFFxFy = forFFxFy;
})(exports.MathFunctionR2R || (exports.MathFunctionR2R = {}));
const cas2 = cassini(0.9, 1.02);

function doNotSerialize(target, key) {
    const map = target.__SERIALIZATION_BLACKLIST || (target.__SERIALIZATION_BLACKLIST = {});
    map[key] = 'no';
}
class ClassSerializer {
    constructor() {
        this.CLASS_NAMES = new Map();
        this.NAME_CLASSES = new Map();
        this.addClass('Object', Object);
    }
    addClass(name, clazz) {
        if (this.NAME_CLASSES.has(name)) {
            throw new Error(name);
        }
        this.NAME_CLASSES.set(name, clazz);
        this.CLASS_NAMES.set(clazz, name);
        return this;
    }
    addNamespace(namespace, namespaceName) {
        Object.keys(namespace).forEach(symbol => {
            const o = namespace[symbol];
            if ('function' == typeof o && o.name) {
                this.addClass((namespaceName ? namespaceName + '.' : '') + symbol, o);
            }
        });
        return this;
    }
    setUpdater(f) {
        this.updater = f;
        return this;
    }
    serialize(v) {
        return JSON.stringify(this.serializeObj(v));
    }
    serializeObj(v) {
        const path = [];
        const gatherList = (v) => {
            //console.log(path.toString())
            if (undefined !== v && v.hasOwnProperty('constructor') && this.CLASS_NAMES.has(v.constructor)) {
                // do nothing, this is a class/function prototype
            }
            else if (Array.isArray(v)) {
                if (visited.has(v)) {
                    if (!listMap.has(v)) {
                        listMap.set(v, resultList.length);
                        resultList.push(v);
                    }
                }
                else {
                    visited.add(v);
                    for (let i = 0; i < v.length; i++) {
                        path.push('' + i);
                        gatherList(v[i]);
                        path.pop();
                    }
                }
            }
            else if (undefined !== v && 'object' == typeof v) {
                if (visited.has(v)) {
                    if (!listMap.has(v)) {
                        listMap.set(v, resultList.length);
                        resultList.push(v);
                    }
                }
                else {
                    ts3dutils.assert(!v.__noxTarget || !visited.has(v.__noxTarget));
                    ts3dutils.assert(!v.__noxProxy || !visited.has(v.__noxProxy));
                    visited.add(v);
                    if (!v.getConstructorParameters) {
                        for (const key of Object.keys(v).sort()) {
                            if (key == '__noxProxy' || key == '__noxTarget')
                                continue;
                            if (!v.__SERIALIZATION_BLACKLIST || !v.__SERIALIZATION_BLACKLIST[key]) {
                                path.push(key);
                                gatherList(v[key]);
                                path.pop();
                            }
                        }
                    }
                    path.push('proto');
                    gatherList(Object.getPrototypeOf(v));
                    path.pop();
                }
            }
        };
        const transform = (v, allowLinks, first) => {
            if ('string' == typeof v || 'number' == typeof v || 'boolean' == typeof v || null === v) {
                return v;
            }
            if ('undefined' == typeof v) {
                return { '#REF': -1 };
            }
            if (v.hasOwnProperty('constructor') && this.CLASS_NAMES.has(v.constructor)) {
                return { '#REF': this.CLASS_NAMES.get(v.constructor) };
            }
            let index;
            if (allowLinks && !first && undefined !== (index = listMap.get(v))) {
                return { '#REF': index };
            }
            if (Array.isArray(v)) {
                return v.map(x => transform(x, allowLinks));
            }
            //if (mobx && mobx.isObservableArray(v)) {
            //	const result = {'#PROTO': 'ObservableArray'} as any
            //	v.forEach((val, i) => result[i] = transform(val))
            //	return result
            //}
            if ('object' == typeof v) {
                if (v.getConstructorParameters) {
                    return {
                        '#CONSTRUCTOR': this.CLASS_NAMES.get(v.constructor),
                        '#ARGS': transform(v.getConstructorParameters(), false)
                    };
                }
                const result = {};
                if (Object.prototype !== Object.getPrototypeOf(v)) {
                    result['#PROTO'] = transform(Object.getPrototypeOf(v), allowLinks);
                }
                for (const key of Object.keys(v)) {
                    if (key == '__noxProxy' || key == '__noxTarget')
                        continue;
                    if (!v.__SERIALIZATION_BLACKLIST || !v.__SERIALIZATION_BLACKLIST[key]) {
                        result[key] = transform(v[key], allowLinks);
                    }
                }
                return result;
            }
            throw new Error('?' + typeof v + v.toString());
        };
        const visited = new Set();
        const listMap = new Map();
        let resultList = [];
        listMap.set(v, 0);
        resultList.push(v);
        gatherList(v);
        resultList = resultList.map(v => transform(v, true, true));
        return resultList;
    }
    unserialize(string) {
        let depth = 0;
        const fixObject = (v, onReady) => {
            depth++;
            if (depth > 100)
                throw new Error();
            if (v && v.constructor === Array) {
                onReady(v);
                for (let i = 0; i < v.length; i++) {
                    fixObject(v[i], x => v[i] = x);
                }
            }
            else if ('object' == typeof v && undefined != v) {
                if ('#CONSTRUCTOR' in v) {
                    const protoName = v['#CONSTRUCTOR'];
                    const proto = this.NAME_CLASSES.get(protoName);
                    ts3dutils.assert(proto, protoName + ' Missing ');
                    let args;
                    fixObject(v['#ARGS'], x => args = x);
                    onReady(new proto(...args));
                }
                else if ('#REF' in v) {
                    const ref = v['#REF'];
                    if ('string' == typeof ref) {
                        onReady(this.NAME_CLASSES.get(ref).prototype);
                    }
                    else if ('number' == typeof ref) {
                        if (-1 == ref) {
                            onReady(undefined);
                        }
                        else if (fixedObjects[ref]) {
                            onReady(fixedObjects[ref]);
                        }
                        else {
                            fixObject(tree[ref], x => onReady(fixedObjects[ref] = x));
                        }
                    }
                }
                else {
                    let result;
                    if ('#PROTO' in v) {
                        fixObject(v['#PROTO'], x => {
                            result = Object.create(x);
                            onReady(result);
                        });
                    }
                    else {
                        onReady(result = v);
                    }
                    const keys = Object.keys(v);
                    for (let i = 0; i < keys.length; i++) {
                        //if ('name' == keys[i]) console.log(result)
                        if ('#PROTO' != keys[i]) {
                            fixObject(v[keys[i]], x => result[keys[i]] = x);
                            //Object.defineProperty(result, keys[i], {
                            //	value: fixObjects(v[keys[i]]),
                            //	enumerable: true,
                            //	writable: true,
                            //	configurable: true
                            //})
                        }
                    }
                    Object.defineProperty(result, 'loadID', { value: getGlobalId(), enumerable: false, writable: false });
                    this.updater && this.updater(result);
                }
            }
            else {
                onReady(v);
            }
            depth--;
        };
        // const linkReferences = (v: any) => {
        // 	if (v && v.constructor === Array) {
        // 		for (let i = 0; i < v.length; i++) {
        // 			v[i] = linkReferences(v[i])
        // 		}
        // 		return v
        // 	} else if ('object' == typeof v && undefined != v) {
        // 		if ('#REF' in v) {
        // 			return tree[v['#REF']]
        // 		} else {
        // 			const keys = Object.keys(v)
        // 			for (let i = 0; i < keys.length; i++) {
        // 				v[keys[i]] = linkReferences(v[keys[i]])
        // 			}
        // 			return v
        // 		}
        // 	} else {
        // 		return v
        // 	}
        // }
        const tree = JSON.parse(string);
        // console.log(tree)
        const fixedObjects = new Array(tree.length);
        fixObject({ '#REF': 0 }, () => { });
        // console.log(tree)
        // linkReferences(tree)
        // console.log(tree)
        return fixedObjects[0];
    }
}

const fragmentShaderLighting = `
	precision highp float;
	uniform vec4 color;
	uniform vec3 camPos;
	varying vec3 normal;
	varying vec4 vPosition;
	void main() {
		vec3 normal1 = normalize(normal);
		vec3 lightPos = vec3(1000, 2000, 4000);
		vec3 lightDir = normalize(vPosition.xyz - lightPos);
        vec3 reflectionDirection = reflect(lightDir, normal1);
        vec3 eyeDirection = normalize(camPos.xyz-vPosition.xyz);
        float uMaterialShininess = 256.0;
		float specularLightWeighting = pow(max(dot(reflectionDirection, eyeDirection), 0.0), uMaterialShininess);
		float lightIntensity = 0.6 + 0.2 * max(0.0, -dot(lightDir, normal1)) + 0.2*specularLightWeighting;
		gl_FragColor = vec4(vec3(color) * lightIntensity, 1);
	}
`;
const vertexShaderLighting = `
	uniform mat4 LGL_ModelViewProjectionMatrix;
	uniform mat4 LGL_ModelViewMatrix;
	attribute vec4 LGL_Vertex;
	uniform mat3 LGL_NormalMatrix;
	attribute vec3 LGL_Normal;
	uniform vec4 color;
	varying vec3 normal;
	varying vec4 vPosition;
	void main() {
		gl_Position = LGL_ModelViewProjectionMatrix * LGL_Vertex;
        vPosition = LGL_ModelViewMatrix * LGL_Vertex;
		normal = normalize(LGL_NormalMatrix * LGL_Normal);
	}
`;
const vertexShaderWaves = `
	uniform mat4 LGL_ModelViewProjectionMatrix;
	uniform mat4 LGL_ModelViewMatrix;
	attribute vec4 LGL_Vertex;
	uniform mat3 LGL_NormalMatrix;
	attribute vec3 LGL_Normal;
	uniform vec4 color;
	varying vec3 normal;
	varying vec4 vPosition;
	void main() {
		normal = normalize(LGL_NormalMatrix * LGL_Normal);
		float offset = mod  (((LGL_Vertex.x + LGL_Vertex.y + LGL_Vertex.z) * 31.0), 20.0) - 10.0;
		vec4 modPos = LGL_Vertex + vec4(normal * offset, 0);
		gl_Position = LGL_ModelViewProjectionMatrix * modPos;
        vPosition = LGL_ModelViewMatrix * modPos;
	}
`;


const vertexShaderBasic = `
	uniform mat4 LGL_ModelViewProjectionMatrix;
	attribute vec4 LGL_Vertex;
	void main() {
		gl_Position = LGL_ModelViewProjectionMatrix * LGL_Vertex;
	}
`;
const vertexShaderColor = `
	uniform mat4 LGL_ModelViewProjectionMatrix;
	attribute vec4 LGL_Vertex;
	attribute vec4 color;
	varying vec4 fragColor;
	void main() {
		gl_Position = LGL_ModelViewProjectionMatrix * LGL_Vertex;
		fragColor = color;
	}
`;
const vertexShaderArc = `
	uniform mat4 LGL_ModelViewProjectionMatrix;
	attribute vec4 LGL_Vertex;
	uniform float step, offset;
	uniform float radius, width;
	void main() {
		float r = radius;
		float t = offset + LGL_Vertex.x * step;
		float pRadius = r - LGL_Vertex.y * width;
		vec4 p = vec4(pRadius * cos(t), pRadius * sin(t), 0, 1);
		gl_Position = LGL_ModelViewProjectionMatrix * p;
}
`;
const vertexShaderConic3d = `
	uniform mat4 LGL_ModelViewProjectionMatrix;
	attribute vec4 LGL_Vertex;
	uniform float startT, endT, scale;
	uniform vec3 center, f1, f2;
	uniform int mode;
	float sinh(float x) { return (exp(x) - exp(-x)) / 2.0; }
	float cosh(float x) { return (exp(x) + exp(-x)) / 2.0; }
	void main() {
		float t = startT + LGL_Vertex.x * (endT - startT);

		vec3 normal = normalize(cross(f1, f2));

		vec3 p, tangent;
		if (0 == mode) { // ellipse
			p = center + f1 * cos(t) + f2 * sin(t);
			tangent = f1 * -sin(t) + f2 * cos(t);
		}
		if (1 == mode) { // parabola
			p = center + f1 * t + f2 * t * t;
			tangent = f1 + f2 * t;
		}
		if (2 == mode) { // hyperbola
			p = center + f1 * cosh(t) + f2 * sinh(t);
			tangent = f1 * sinh(t) + f2 * cosh(t);
		}
		vec3 outDir = normalize(cross(normal, tangent));
		vec3 p2 = p + scale * (outDir * LGL_Vertex.y + normal * LGL_Vertex.z);
		gl_Position = LGL_ModelViewProjectionMatrix * vec4(p2, 1);
	}
`;
const vertexShaderBezier = `
    // calculates a bezier curve using LGL_Vertex.x as the (t) parameter of the curve
	uniform mat4 LGL_ModelViewProjectionMatrix;
	attribute vec4 LGL_Vertex;
	uniform float width, startT, endT;
	uniform vec3 p0, p1, p2, p3;
	void main() {
		// LGL_Vertex.y is in [0, 1]
		float t = startT + LGL_Vertex.x * (endT - startT), s = 1.0 - t;
		float c0 = s * s * s, c1 = 3.0 * s * s * t, c2 = 3.0 * s * t * t, c3 = t * t * t;
		vec3 pPos = p0 * c0 + p1 * c1 + p2 * c2 + p3 * c3;
		float c01 = 3.0 * s * s, c12 = 6.0 * s * t, c23 = 3.0 * t * t;
		vec3 pTangent = (p1 - p0) * c01 + (p2 - p1) * c12 + (p3 - p2) * c23;
		vec3 pNormal = normalize(vec3(pTangent.y, -pTangent.x, 0));
		vec4 p = vec4(pPos - LGL_Vertex.y * width * pNormal, 1);
		gl_Position = LGL_ModelViewProjectionMatrix * p;
	}
`;
const vertexShaderBezier3d = `
    // calculates a bezier curve using LGL_Vertex.x as the (t) parameter of the curve
	uniform float scale, startT, endT;
	uniform vec3 ps[4];
	uniform vec3 p0, p1, p2, p3, normal;
	uniform mat4 LGL_ModelViewProjectionMatrix;
	attribute vec4 LGL_Vertex;
	void main() {
		// LGL_Vertex.y is in [0, 1]
		vec3 p5 = ps[0];
		float t = startT + LGL_Vertex.x * (endT - startT), s = 1.0 - t;
		float c0 = s * s * s, c1 = 3.0 * s * s * t, c2 = 3.0 * s * t * t, c3 = t * t * t;
		vec3 p = p0 * c0 + p1 * c1 + p2 * c2 + p3 * c3;
		float c01 = 3.0 * s * s, c12 = 6.0 * s * t, c23 = 3.0 * t * t;
		vec3 pTangent = (p1 - p0) * c01 + (p2 - p1) * c12 + (p3 - p2) * c23;
		vec3 outDir = normalize(cross(normal, pTangent));
		vec3 correctNormal = normalize(cross(pTangent, outDir));
		vec3 p2 = p + scale * (outDir * LGL_Vertex.y + correctNormal * LGL_Vertex.z);
		gl_Position = LGL_ModelViewProjectionMatrix * vec4(p2, 1);
	}
`;
const vertexShaderGeneric = `
	uniform float scale;
	uniform mat4 LGL_ModelViewProjectionMatrix;
	attribute vec4 LGL_Vertex;
	uniform mat3 LGL_NormalMatrix;
	attribute vec3 LGL_Normal;
	void main() {
		vec3 normal = normalize(LGL_NormalMatrix * LGL_Normal);
		vec4 vertexPos = LGL_Vertex + vec4(normal * scale, 0);
		gl_Position = LGL_ModelViewProjectionMatrix * vertexPos;
	}
`;
const vertexShaderRing = `
	#define M_PI 3.1415926535897932384626433832795
	uniform float step;
	uniform float innerRadius, outerRadius;
	attribute float index;
	uniform mat4 LGL_ModelViewProjectionMatrix;
	attribute vec4 LGL_Vertex;
	void main() {
		gl_Position = LGL_ModelViewProjectionMatrix * vec4(index, index, index, 1);
		float id = atan(LGL_Vertex.x, LGL_Vertex.y) / M_PI  * 32.0;
		float radius = mod(id, 2.0) < 1.0 ? outerRadius : innerRadius;
		gl_Position = LGL_ModelViewProjectionMatrix * vec4(radius * cos(index * step), radius * sin(index * step), 0, 1);
	}
`;
const fragmentShaderColor = `
	precision highp float;
	uniform vec4 color;
	void main() {
		gl_FragColor = color;
	}
`;
const fragmentShaderVaryingColor = `
	precision highp float;
	varying vec4 fragColor;
	void main() {
		gl_FragColor = fragColor;
	}
`;
const fragmentShaderColorHighlight = `
	precision highp float;
	uniform vec4 color;
	void main() {
		float diagonal = (gl_FragCoord.x + 2.0 * gl_FragCoord.y);
		if (mod(diagonal, 50.0) > 40.0) { // mod(diagonal, 2.0) > 1.0
			discard;
			//gl_FragColor = color + vec4(0.2,0.2,0.2,0);
		} else {
			gl_FragColor = color - vec4(0.2,0.2,0.2,0);
		}
	}
`;
const vertexShaderTexture = `
	varying vec2 texturePos;
	attribute vec4 LGL_Vertex;
	uniform mat4 LGL_ModelViewProjectionMatrix;
	void main() {
		texturePos = LGL_Vertex.xy;
		gl_Position = LGL_ModelViewProjectionMatrix * LGL_Vertex;
	}
`;
const fragmentShaderTextureColor = `
	precision highp float;
	varying vec2 texturePos;
	uniform vec4 color;
	uniform sampler2D texture;
	void main() {
		gl_FragColor = texture2D(texture, texturePos) * color;
	}
`;

const { pow: pow$5, sign: sign$11 } = Math;
function parseGetParams(str) {
    const result = {};
    str
        .split('&')
        .forEach(function (item) {
        const splitIndex = item.indexOf('=');
        if (-1 == splitIndex) {
            result[item] = item;
        }
        else {
            result[item.substr(0, splitIndex)] = decodeURI(item.substr(splitIndex + 1));
        }
    });
    return result;
}
const COLORS = {
    RD_FILL: chroma('#9EDBF9'),
    RD_STROKE: chroma('#77B0E0'),
    TS_FILL: chroma('#D19FE3'),
    TS_STROKE: chroma('#A76BC2'),
    PP_FILL: chroma('#F3B6CF'),
    PP_STROKE: chroma('#EB81B4'),
};
class BREPGLContext {
    constructor(gl) {
        this.cachedMeshes = new WeakMap();
        this.shaders = initShaders(gl);
        initMeshes(this.meshes = {}, gl);
    }
    static create(gl) {
        ts3dutils.addOwnProperties(gl, BREPGLContext.prototype);
        ts3dutils.addOwnProperties(gl, new BREPGLContext(gl));
        return gl;
    }
    drawPoint(p, color = tsgl.GL_COLOR_BLACK, size = 5) {
        this.pushMatrix();
        this.translate(p);
        this.scale(size, size, size);
        this.shaders.singleColor.uniforms({ color: color }).draw(this.meshes.sphere1);
        this.popMatrix();
    }
    drawEdge(edge, color = tsgl.GL_COLOR_BLACK, width = 2) {
        CURVE_PAINTERS[edge.curve.constructor.name](this, edge.curve, color, edge.minT, edge.maxT, width);
    }
    drawCurve(curve, color = tsgl.GL_COLOR_BLACK, width = 2, tStart, tEnd) {
        CURVE_PAINTERS[curve.constructor.name](this, curve, color, tStart, tEnd, width);
    }
    drawVector(vector, anchor, color = tsgl.GL_COLOR_BLACK, size = 1) {
        this.pushMatrix();
        const vT = vector.getPerpendicular().unit();
        this.multMatrix(ts3dutils.M4.forSys(vector, vT, vector.cross(vT).unit(), anchor));
        1 != size && this.scale(size, size, size);
        this.shaders.singleColor.uniforms({
            color: color,
        }).draw(this.meshes.vector);
        this.popMatrix();
    }
    drawVectors(drVs) {
        this.drawVector(ts3dutils.V3.X, ts3dutils.V3.O, chroma('red').gl(), undefined);
        this.drawVector(ts3dutils.V3.Y, ts3dutils.V3.O, chroma('green').gl(), undefined);
        this.drawVector(ts3dutils.V3.Z, ts3dutils.V3.O, chroma('blue').gl(), undefined);
        drVs.forEach(vi => this.drawVector(vi.dir1, vi.anchor, vi.color, undefined));
    }
    drawPlane(customPlane, color, dotted = false) {
        this.pushMatrix();
        this.multMatrix(ts3dutils.M4.forSys(customPlane.right, customPlane.up, customPlane.normal1));
        this.translate(customPlane.sMin, customPlane.tMin, customPlane.w);
        this.scale(customPlane.sMax - customPlane.sMin, customPlane.tMax - customPlane.tMin, 1);
        const mesh = dotted ? this.meshes.xyDottedLinePlane : this.meshes.xyLinePlane;
        this.shaders.singleColor.uniforms({ color: color }).draw(mesh, tsgl.DRAW_MODES.LINES);
        this.popMatrix();
    }
}
function conicPainter(mode, gl, ellipse, color, startT, endT, width = 2) {
    gl.shaders.ellipse3d.uniforms({
        f1: ellipse.f1,
        f2: ellipse.f2,
        center: ellipse.center,
        color: color,
        startT: startT,
        endT: endT,
        scale: width,
        mode: mode,
    }).draw(gl.meshes.pipe);
}
const CURVE_PAINTERS = {
    [SemiEllipseCurve.name]: conicPainter.bind(undefined, 0),
    [EllipseCurve.name]: conicPainter.bind(undefined, 0),
    [ParabolaCurve.name]: conicPainter.bind(undefined, 1),
    [HyperbolaCurve.name]: conicPainter.bind(undefined, 2),
    [ImplicitCurve.name](gl, curve, color, startT, endT, width = 2, normal = ts3dutils.V3.Z) {
        let mesh = gl.cachedMeshes.get(curve);
        if (!mesh) {
            mesh = new tsgl.Mesh()
                .addIndexBuffer('TRIANGLES')
                .addVertexBuffer('normals', 'LGL_Normal');
            curve.addToMesh(mesh);
            mesh.compile();
            //mesh=Mesh.sphere(2)
            gl.cachedMeshes.set(curve, mesh);
        }
        // TODO: draw only part
        //startT: startT,
        //	endT: endT,
        gl.shaders.generic3d.uniforms({
            color: color,
            scale: width,
        }).draw(mesh);
    },
    [BezierCurve.name](gl, curve, color, startT, endT, width = 2, normal = ts3dutils.V3.Z) {
        gl.shaders.bezier3d.uniforms({
            p0: curve.p0,
            p1: curve.p1,
            p2: curve.p2,
            p3: curve.p3,
            color: color,
            startT: startT,
            endT: endT,
            scale: width,
            normal: normal,
        }).draw(gl.meshes.pipe);
    },
    [L3$1.name](gl, curve, color, startT, endT, width = 2, normal = ts3dutils.V3.Z) {
        gl.pushMatrix();
        const a = curve.at(startT), b = curve.at(endT);
        const ab = b.minus(a), abT = ab.getPerpendicular().unit();
        const m = ts3dutils.M4.forSys(ab, abT, ab.cross(abT).unit(), a);
        gl.multMatrix(m);
        gl.scale(1, width, width);
        gl.shaders.singleColor.uniforms({
            color: color,
        }).draw(gl.meshes.pipe);
        gl.popMatrix();
    },
};
CURVE_PAINTERS[PICurve$1.name] = CURVE_PAINTERS[ImplicitCurve.name];
function initMeshes(_meshes, _gl) {
    _gl.makeCurrent();
    _meshes.sphere1 = tsgl.Mesh.sphere(2);
    _meshes.segment = tsgl.Mesh.plane({ startY: -0.5, height: 1, detailX: 128 });
    _meshes.text = tsgl.Mesh.plane();
    _meshes.vector = tsgl.Mesh.rotation([ts3dutils.V3.O, ts3dutils.V(0, 0.05, 0), ts3dutils.V(0.8, 0.05), ts3dutils.V(0.8, 0.1), ts3dutils.V(1, 0)], L3$1.X, ts3dutils.TAU, 16, true);
    _meshes.pipe = tsgl.Mesh.rotation(ts3dutils.arrayFromFunction(128, i => new ts3dutils.V3(i / 127, -0.5, 0)), L3$1.X, ts3dutils.TAU, 8, true);
    _meshes.xyLinePlane = tsgl.Mesh.plane();
    _meshes.xyDottedLinePlane = makeDottedLinePlane();
}
function initShaders(_gl) {
    _gl.makeCurrent();
    return {
        singleColor: tsgl.Shader.create(vertexShaderBasic, fragmentShaderColor),
        multiColor: tsgl.Shader.create(vertexShaderColor, fragmentShaderVaryingColor),
        singleColorHighlight: tsgl.Shader.create(vertexShaderBasic, fragmentShaderColorHighlight),
        textureColor: tsgl.Shader.create(vertexShaderTexture, fragmentShaderTextureColor),
        arc: tsgl.Shader.create(vertexShaderRing, fragmentShaderColor),
        arc2: tsgl.Shader.create(vertexShaderArc, fragmentShaderColor),
        ellipse3d: tsgl.Shader.create(vertexShaderConic3d, fragmentShaderColor),
        generic3d: tsgl.Shader.create(vertexShaderGeneric, fragmentShaderColor),
        bezier3d: tsgl.Shader.create(vertexShaderBezier3d, fragmentShaderColor),
        bezier: tsgl.Shader.create(vertexShaderBezier, fragmentShaderColor),
        lighting: tsgl.Shader.create(vertexShaderLighting, fragmentShaderLighting),
        waves: tsgl.Shader.create(vertexShaderWaves, fragmentShaderLighting),
    };
}
function makeDottedLinePlane(count = 128) {
    const mesh = new tsgl.Mesh().addIndexBuffer('LINES');
    const OXvertices = ts3dutils.arrayFromFunction(count, i => new ts3dutils.V3(i / count, 0, 0));
    mesh.vertices.push(...OXvertices);
    mesh.vertices.push(...ts3dutils.M4.forSys(ts3dutils.V3.Y, ts3dutils.V3.O, ts3dutils.V3.O, ts3dutils.V3.X).transformedPoints(OXvertices));
    mesh.vertices.push(...ts3dutils.M4.forSys(ts3dutils.V3.X.negated(), ts3dutils.V3.O, ts3dutils.V3.O, new ts3dutils.V3(1, 1, 0)).transformedPoints(OXvertices));
    mesh.vertices.push(...ts3dutils.M4.forSys(ts3dutils.V3.Y.negated(), ts3dutils.V3.O, ts3dutils.V3.O, ts3dutils.V3.Y).transformedPoints(OXvertices));
    mesh.LINES = ts3dutils.arrayFromFunction(count * 4, i => i - (i >= count * 2 ? 1 : 0));
    mesh.compile();
    return mesh;
}
function initNavigationEvents(_gl, eye, paintScreen) {
    const canvas = _gl.canvas;
    let lastPos = ts3dutils.V3.O;
    //_gl.onmousedown.push((e) => {
    //	e.preventDefault()
    //	e.stopPropagation()
    //})
    //_gl.onmouseup.push((e) => {
    //	e.preventDefault()
    //	e.stopPropagation()
    //})
    canvas.addEventListener('mousemove', e => {
        const pagePos = ts3dutils.V(e.pageX, e.pageY);
        const delta = lastPos.to(pagePos);
        //noinspection JSBitwiseOperatorUsage
        if (e.buttons & 4) {
            // pan
            const moveCamera = ts3dutils.V(-delta.x * 2 / _gl.canvas.width, delta.y * 2 / _gl.canvas.height);
            const inverseProjectionMatrix = _gl.projectionMatrix.inversed();
            const worldMoveCamera = inverseProjectionMatrix.transformVector(moveCamera);
            eye.pos = eye.pos.plus(worldMoveCamera);
            eye.focus = eye.focus.plus(worldMoveCamera);
            setupCamera(eye, _gl);
            paintScreen();
        }
        // scene rotation
        //noinspection JSBitwiseOperatorUsage
        if (e.buttons & 2) {
            const rotateLR = -delta.x / 6.0 * ts3dutils.DEG;
            const rotateUD = -delta.y / 6.0 * ts3dutils.DEG;
            // rotate
            let matrix = ts3dutils.M4.rotateLine(eye.focus, eye.up, rotateLR);
            //let horizontalRotationAxis = focus.minus(pos).cross(up)
            const horizontalRotationAxis = eye.up.cross(eye.pos.minus(eye.focus));
            matrix = matrix.times(ts3dutils.M4.rotateLine(eye.focus, horizontalRotationAxis, rotateUD));
            eye.pos = matrix.transformPoint(eye.pos);
            eye.up = matrix.transformVector(eye.up);
            setupCamera(eye, _gl);
            paintScreen();
        }
        lastPos = pagePos;
    });
    canvas.addEventListener('wheel', function (e) {
        // zoom
        const wheelY = -sign$11(e.deltaY) * 2;
        // console.log(e.deltaY, e.deltaX)
        eye.zoomFactor *= pow$5(0.9, -wheelY);
        const mouseCoordsOnCanvas = getPosOnTarget(e);
        const mousePosFrustrum = ts3dutils.V(mouseCoordsOnCanvas.x * 2 / _gl.canvas.offsetWidth - 1, -mouseCoordsOnCanvas.y * 2 / _gl.canvas.offsetHeight + 1, 0);
        const moveCamera = mousePosFrustrum.times(1 - 1 / pow$5(0.9, -wheelY));
        const inverseProjectionMatrix = _gl.projectionMatrix.inversed();
        const worldMoveCamera = inverseProjectionMatrix.transformVector(moveCamera);
        //console.log("moveCamera", moveCamera)
        //console.log("worldMoveCamera", worldMoveCamera)
        eye.pos = eye.pos.plus(worldMoveCamera);
        eye.focus = eye.focus.plus(worldMoveCamera);
        // tilt
        const mousePosWC = inverseProjectionMatrix.transformPoint(mousePosFrustrum);
        const tiltMatrix = ts3dutils.M4.rotateLine(mousePosWC, eye.pos.to(eye.focus), -sign$11(e.deltaX) * 10 * ts3dutils.DEG);
        eye.up = tiltMatrix.transformVector(eye.up);
        eye.pos = tiltMatrix.transformPoint(eye.pos);
        eye.focus = tiltMatrix.transformPoint(eye.focus);
        setupCamera(eye, _gl);
        paintScreen();
        e.preventDefault();
    });
}
/**
 * Transforms position on the screen into a line in world coordinates.
 */
function getMouseLine(pos, _gl) {
    const ndc1 = ts3dutils.V(pos.x * 2 / _gl.canvas.width - 1, -pos.y * 2 / _gl.canvas.height + 1, 0);
    const ndc2 = ts3dutils.V(pos.x * 2 / _gl.canvas.width - 1, -pos.y * 2 / _gl.canvas.height + 1, 1);
    //console.log(ndc)
    const inverseProjectionMatrix = _gl.projectionMatrix.inversed();
    const s = inverseProjectionMatrix.transformPoint(ndc1);
    const dir = inverseProjectionMatrix.transformPoint(ndc2).minus(s);
    return L3$1.anchorDirection(s, dir);
}
function getPosOnTarget(e) {
    const target = e.target;
    const targetRect = target.getBoundingClientRect();
    const mouseCoordsOnElement = {
        x: e.clientX - targetRect.left,
        y: e.clientY - targetRect.top
    };
    return mouseCoordsOnElement;
}
function setupCamera(_eye, _gl) {
    const { pos, focus, up, zoomFactor } = _eye;
    //console.log("pos", pos.$, "focus", focus.$, "up", up.$)
    _gl.matrixMode(_gl.PROJECTION);
    _gl.loadIdentity();
    //_gl.perspective(70, _gl.canvas.width / _gl.canvas.height, 0.1, 1000);
    const lr = _gl.canvas.width / 2 / zoomFactor;
    const bt = _gl.canvas.height / 2 / zoomFactor;
    _gl.ortho(-lr, lr, -bt, bt, -1e4, 1e4);
    _gl.lookAt(pos, focus, up);
    _gl.matrixMode(_gl.MODELVIEW);
    cameraChangeListeners.forEach(l => l(_eye));
}
const cameraChangeListeners = [];
const SHADERS_TYPE_VAR = false && initShaders(0);
// let shaders: typeof SHADERS_TYPE_VAR
// declare let a: B2, b: B2, c: B2, d: B2, edges: Edge[] = [], hovering: any,
// 	, normallines: boolean = false, b2s: B2[] = []
// const

exports.Curve = Curve;
exports.curvePoint = curvePoint;
exports.curvePointMF = curvePointMF;
exports.XiEtaCurve = XiEtaCurve;
exports.ImplicitCurve = ImplicitCurve;
exports.BezierCurve = BezierCurve;
exports.EllipseCurve = EllipseCurve;
exports.HyperbolaCurve = HyperbolaCurve;
exports.L3 = L3$1;
exports.PICurve = PICurve$1;
exports.ParabolaCurve = ParabolaCurve;
exports.SemiEllipseCurve = SemiEllipseCurve;
exports.P3 = P3;
exports.Surface = Surface;
exports.ParametricSurface = ParametricSurface;
exports.ImplicitSurface = ImplicitSurface;
exports.ConicSurface = ConicSurface;
exports.EllipsoidSurface = EllipsoidSurface;
exports.ProjectedCurveSurface = ProjectedCurveSurface;
exports.CylinderSurface = CylinderSurface;
exports.RotationREqFOfZ = RotationREqFOfZ;
exports.SemiCylinderSurface = SemiCylinderSurface;
exports.SemiEllipsoidSurface = SemiEllipsoidSurface;
exports.PlaneSurface = PlaneSurface$1;
exports.ZDirVolumeVisitor = ZDirVolumeVisitor;
exports.CalculateAreaVisitor = CalculateAreaVisitor;
exports.CustomPlane = CustomPlane;
exports.Edge = Edge;
exports.PCurveEdge = PCurveEdge;
exports.StraightEdge = StraightEdge;
exports.FaceInfoFactory = FaceInfoFactory;
exports.Face = Face;
exports.PlaneFace = PlaneFace;
exports.RotationFace = RotationFace;
exports.EPS = EPS;
exports.getGlobalId = getGlobalId;
exports.addLikeSurfaceFaces = addLikeSurfaceFaces;
exports.assembleFaceFromLooseEdges = assembleFaceFromLooseEdges;
exports.calcNextEdgeIndex = calcNextEdgeIndex;
exports.B2 = B2;
exports.dotCurve = dotCurve;
exports.dotCurve2 = dotCurve2;
exports.INSIDE = INSIDE;
exports.OUTSIDE = OUTSIDE$1;
exports.COPLANAR_SAME = COPLANAR_SAME;
exports.COPLANAR_OPPOSITE = COPLANAR_OPPOSITE;
exports.ALONG_EDGE_OR_PLANE = ALONG_EDGE_OR_PLANE;
exports.splitsVolumeEnclosingFaces = splitsVolumeEnclosingFaces;
exports.splitsVolumeEnclosingFacesP = splitsVolumeEnclosingFacesP;
exports.splitsVolumeEnclosingFacesP2 = splitsVolumeEnclosingFacesP2;
exports.splitsVolumeEnclosingCone = splitsVolumeEnclosingCone;
exports.splitsVolumeEnclosingCone2 = splitsVolumeEnclosingCone2;
exports.fff = fff;
exports.triangulateVertices = triangulateVertices;
exports.intersectionUnitCircleLine = intersectionUnitCircleLine;
exports.intersectionUnitCircleLine2 = intersectionUnitCircleLine2;
exports.intersectionCircleLine = intersectionCircleLine;
exports.intersectionUnitHyperbolaLine = intersectionUnitHyperbolaLine;
exports.followAlgorithm2d = followAlgorithm2d;
exports.followAlgorithm2dAdjustable = followAlgorithm2dAdjustable;
exports.intersectionICurveICurve = intersectionICurveICurve;
exports.intersectionICurveICurve2 = intersectionICurveICurve2;
exports.intersectionICurvePSurface = intersectionICurvePSurface;
exports.cassini = cassini;
exports.cas2 = cas2;
exports.doNotSerialize = doNotSerialize;
exports.ClassSerializer = ClassSerializer;
exports.parseGetParams = parseGetParams;
exports.COLORS = COLORS;
exports.BREPGLContext = BREPGLContext;
exports.CURVE_PAINTERS = CURVE_PAINTERS;
exports.initMeshes = initMeshes;
exports.initShaders = initShaders;
exports.initNavigationEvents = initNavigationEvents;
exports.getMouseLine = getMouseLine;
exports.getPosOnTarget = getPosOnTarget;
exports.setupCamera = setupCamera;
exports.cameraChangeListeners = cameraChangeListeners;
exports.SHADERS_TYPE_VAR = SHADERS_TYPE_VAR;
//# sourceMappingURL=bundle.js.map
