import { Transformable, assertNumbers, assert, fuzzyUniquesF, eq0, callsce, le, arrayFromFunction, newtonIterateWithDerivative, clamp, AABB, V3, newtonIterate1d, eq, glqInSteps, hasConstructor, getIntervals, V, NLA_PRECISION, newtonIterate2dWithDerivatives, assertVectors, M4, TAU, mapFilter, assertInst, DEG, arrayRange, bisect, solveCubicReal2, between, assertNever, assertf, fuzzyUniques, MINUS, combinations, lerp, snap0, lt, VV, newtonIterate, pqFormula, fuzzyBetween, checkDerivate, newtonIterateSmart, NLA_DEBUG, firstUnsorted, Vector, arraySamples, newtonIterateWithDerivative2, vArrGet, snap, floatHashCode, arrayEquals, arrayHashCode, isCCW, getRoots, gaussLegendreQuadrature24, toSource, sliceStep, gaussLegendre24Xs, gaussLegendre24Weights, GOLDEN_RATIO, getLast, min as min$1, snap2, PI as PI$1, mod, gt, doubleSignedArea, concatenated, disableConsole, ge, indexWithMax, bagRemoveIndex, enableConsole, binaryIndexOf, binaryInsert, mapPush, sum, SCE, withMax, newtonIterate2d, addOwnProperties } from 'ts3dutils';
import { pushQuad, Mesh, GL_COLOR_BLACK, Shader } from 'tsgl';
import { P3 as P3$1, dotCurve2 as dotCurve2$1, PPCurve as PPCurve$1, ImplicitCurve as ImplicitCurve$1, ZDirVolumeVisitor as ZDirVolumeVisitor$1, CalculateAreaVisitor as CalculateAreaVisitor$1, Surface as Surface$1, MathFunctionR2R as MathFunctionR2R$1, Curve as Curve$1, PICurve as PICurve$1, breakDownPPCurves as breakDownPPCurves$1, ParametricSurface as ParametricSurface$1, EllipseCurve as EllipseCurve$1, L3 as L3$1, HyperbolaCurve as HyperbolaCurve$1, ParabolaCurve as ParabolaCurve$1, CylinderSurface as CylinderSurface$1, PlaneSurface as PlaneSurface$1, ImplicitSurface as ImplicitSurface$1, NURBS as NURBS$1, NURBSSurface as NURBSSurface$1, intersectionUnitCircleLine2 as intersectionUnitCircleLine2$1, ProjectedCurveSurface as ProjectedCurveSurface$1, OUTSIDE as OUTSIDE$1, BezierCurve as BezierCurve$1, PointVsFace as PointVsFace$1, Edge as Edge$1, getExtremePointsHelper as getExtremePointsHelper$1 } from '..';
import { load, Path } from 'opentype.js';
import chroma from 'chroma-js';
import { SVGPathData } from 'svg-pathdata';
import { Pair, JavaMap, JavaSet } from 'javasetmap.ts';
import earcut from 'earcut';
import nerdamer from 'nerdamer';

const { abs, acos, acosh, asin, asinh, atan, atanh, atan2, ceil, cbrt, expm1, clz32, cos, cosh, exp, floor, fround, hypot, imul, log, log1p, log2, log10, max, min, pow, random, round, sign, sin, sinh, sqrt, tan, tanh, trunc, E, LN10, LN2, LOG10E, LOG2E, PI, SQRT1_2, SQRT2, } = Math;

let insideIsInfosWithCurve = false;
class Curve extends Transformable {
    constructor(tMin, tMax) {
        super();
        this.tMin = tMin;
        this.tMax = tMax;
        assertNumbers(tMin, tMax);
        assert("number" == typeof tMin && !isNaN(tMin));
        assert("number" == typeof tMax && !isNaN(tMax));
        assert(tMin < tMax, "tMin < tMax " + tMin + " < " + tMax);
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
            if (!result.some((info) => eq(info.tThis, startT) && eq(info.tOther, startS))) {
                const f1 = (t, s) => curve1.tangentAt(t).dot(curve1.at(t).minus(curve2.at(s)));
                const f2 = (t, s) => curve2.tangentAt(s).dot(curve1.at(t).minus(curve2.at(s)));
                // f = (b1, b2, t1, t2) = b1.tangentAt(t1).dot(b1.at(t1).minus(b2.at(t2)))
                const dfdt1 = (b1, b2, t1, t2) => b1.ddt(t1).dot(b1.at(t1).minus(b2.at(t2))) +
                    b1.tangentAt(t1).squared();
                const dfdt2 = (b1, b2, t1, t2) => -b1.tangentAt(t1).dot(b2.tangentAt(t2));
                const ni = newtonIterate2dWithDerivatives(f1, f2, startT, startS, 16, dfdt1.bind(undefined, curve1, curve2), dfdt2.bind(undefined, curve1, curve2), (t, s) => -dfdt2(curve2, curve1, s, t), (t, s) => -dfdt1(curve2, curve1, s, t));
                assert(isFinite(ni.x));
                assert(isFinite(ni.y));
                if (ni == undefined)
                    console.log(startT, startS, curve1.sce, curve2.sce);
                result.push({ tThis: ni.x, tOther: ni.y, p: curve1.at(ni.x) });
            }
        }
        // returns whether an intersection was immediately found (i.e. without further recursion)
        function findRecursive(tMin, tMax, sMin, sMax, curve1AABB, curve2AABB, depth = 0) {
            const EPS = NLA_PRECISION;
            if (curve1AABB.touchesAABBfuzzy(curve2AABB)) {
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
                    findRecursive(tMin, tMid, sMin, sMid, curve1AABBleft, curve2AABBleft, depth + 1) ||
                        findRecursive(tMin, tMid, sMid, sMax, curve1AABBleft, (curve2AABBright = curve2.getAABB(sMid, sMax)), depth + 1) ||
                        findRecursive(tMid, tMax, sMin, sMid, (curve1AABBright = curve1.getAABB(tMid, tMax)), curve2AABBleft, depth + 1) ||
                        findRecursive(tMid, tMax, sMid, sMax, curve1AABBright, curve2AABBright, depth + 1);
                }
            }
            return false;
        }
        const result = [];
        findRecursive(tMin, tMax, sMin, sMax, curve1.getAABB(tMin, tMax), curve2.getAABB(sMin, sMax));
        return fuzzyUniquesF(result, (info) => info.tThis);
    }
    /**
     * Searches a 2d area for (an) implicit curve(s).
     * @param implicitCurve
     * @param bounds Defines area to search.
     * @param uStep Granularity of search in s-direction.
     * @param vStep Granularity of search in t-direction.
     * @param stepSize step size to take along the curve
     * @return
     */
    static breakDownIC(implicitCurve, bounds, uStep, vStep, stepSize, validUV) {
        //undefined == didu && (didu = (u, v) => (implicitCurve(u + EPS, v) - implicitCurve(u, v)) / EPS)
        //undefined == didv && (didv = (u, v) => (implicitCurve(u, v + EPS) - implicitCurve(u, v)) / EPS)
        const { uMin, uMax, vMin, vMax } = bounds;
        const deltaS = uMax - uMin, deltaT = vMax - vMin;
        const sRes = ceil(deltaS / uStep), tRes = ceil(deltaT / vStep);
        const grid = new Array(sRes * tRes).fill(0);
        // const printGrid = () =>
        // 	console.log(
        // 		arrayFromFunction(tRes, i =>
        // 			grid
        // 				.slice(sRes * i, sRes * (i + 1))
        // 				.map(v => (v ? 'X' : '_'))
        // 				.join(''),
        // 		).join('\n'),
        // 	)
        const get = (i, j) => grid[j * sRes + i];
        const set = (i, j) => 0 <= i && i < sRes && 0 <= j && j < tRes && (grid[j * sRes + i] = 1);
        const result = [];
        const logTable = [];
        for (let i = 0; i < sRes; i++) {
            search: for (let j = 0; j < tRes; j++) {
                if (get(i, j))
                    continue;
                set(i, j);
                let u = uMin + (i + 0.5) * uStep, v = vMin + (j + 0.5) * vStep;
                const startS = u, startT = v;
                // basically curvePoint
                for (let k = 0; k < 8; k++) {
                    const fp = implicitCurve(u, v);
                    const dfpdx = implicitCurve.x(u, v), dfpdy = implicitCurve.y(u, v);
                    if (0 === Math.pow(dfpdx, 2) + Math.pow(dfpdy, 2)) {
                        // top of a hill, keep looking
                        continue search;
                    }
                    const scale = fp / (Math.pow(dfpdx, 2) + Math.pow(dfpdy, 2));
                    u -= scale * dfpdx;
                    v -= scale * dfpdy;
                }
                const li = floor((u - uMin) / uStep), lj = floor((v - vMin) / vStep);
                logTable.push({
                    i,
                    j,
                    li,
                    lj,
                    startS,
                    startT,
                    u,
                    v,
                    "bounds(u, v)": uvInAABB2(bounds, u, v),
                    "ic(s,t)": implicitCurve(u, v),
                });
                if (!(i == li && j == lj) && get(li, lj)) {
                    continue search;
                }
                set(li, lj);
                // u, v are now good starting coordinates to use follow algorithm
                if (uvInAABB2(bounds, u, v) &&
                    validUV(u, v) &&
                    eq0(implicitCurve(u, v))) {
                    const subResult = mkcurves(implicitCurve, u, v, stepSize, bounds, validUV);
                    for (const curveData of subResult) {
                        assert(curveData.points.length > 2);
                        for (const { x, y } of curveData.points) {
                            const lif = (x - uMin) / uStep, ljf = (y - vMin) / vStep;
                            set((lif - 0.5) | 0, (ljf - 0.5) | 0);
                            set((lif - 0.5) | 0, (ljf + 0.5) | 0);
                            set((lif + 0.5) | 0, (ljf - 0.5) | 0);
                            set((lif + 0.5) | 0, (ljf + 0.5) | 0);
                        }
                    }
                    //printGrid()
                    result.push(...subResult);
                }
            }
        }
        // console.table(logTable)
        for (const { points } of result) {
            for (let i = 0; i < points.length - 1; i++) {
                assert(!points[i].equals(points[i + 1]));
            }
        }
        return result;
    }
    toString() {
        return this.toSource();
    }
    toSource(rounder = (x) => x) {
        return callsce.call(undefined, "new " + this.constructor.name, ...this.getConstructorParameters(), this.tMin, this.tMax);
    }
    withBounds(tMin = this.tMin, tMax = this.tMax) {
        //assert(this.tMin <= tMin && tMin <= this.tMax)
        //assert(this.tMin <= tMax && tMax <= this.tMax)
        return new this.constructor(...this.getConstructorParameters(), tMin, tMax);
    }
    /**
     * The point on the line that is closest to the given point.
     */
    closestPointToPoint(p) {
        return this.at(this.closestTToPoint(p));
    }
    isValidT(t) {
        return le(this.tMin, t) && le(t, this.tMax);
    }
    diff(t, eps) {
        return this.at(t).to(this.at(t + eps));
    }
    // TODO: tmin/tmax first
    closestTToPoint(p, tStart, tMin = this.tMin, tMax = this.tMax) {
        // this.at(t) has minimal distance to p when this.tangentAt(t) is perpendicular to
        // the vector between this.at(t) and p. This is the case iff the dot product of the two is 0.
        // f = (this.at(t) - p) . (this.tangentAt(t)
        // df = this.tangentAt(t) . this.tangentAt(t) + (this.at(t) - p) . this.ddt(t)
        //    = this.tangentAt(t)² + (this.at(t) - p) . this.ddt(t)
        const f = (t) => this.at(t).minus(p).dot(this.tangentAt(t)); // 5th degree polynomial
        const df = (t) => this.tangentAt(t).squared() + this.at(t).minus(p).dot(this.ddt(t));
        //checkDerivate(f, df, tMin, tMax)
        const STEPS = 32;
        if (undefined === tStart) {
            tStart = arrayFromFunction(STEPS, (i) => tMin + ((tMax - tMin) * i) / (STEPS - 1)).withMax((t) => -this.at(t).distanceTo(p));
        }
        return newtonIterateWithDerivative(f, tStart, 16, df);
    }
    /**
     * So different edges on the same curve do not have different vertices, they are always generated
     * on fixed points this.at(k * this.tIncrement), with k taking integer values
     *
     */
    calcSegmentPoints(aT, bT, a, b, reversed, includeFirst) {
        assert(this.tIncrement, "tIncrement not defined on " + this);
        const inc = this.tIncrement;
        const result = [];
        if (includeFirst)
            result.push(a);
        assert(reversed != aT < bT);
        if (aT < bT) {
            const start = Math.ceil((aT + NLA_PRECISION) / inc);
            const end = Math.floor((bT - NLA_PRECISION) / inc);
            for (let i = start; i <= end; i++) {
                result.push(this.at(i * inc));
            }
        }
        else {
            const start = Math.floor((aT - NLA_PRECISION) / inc);
            const end = Math.ceil((bT + NLA_PRECISION) / inc);
            for (let i = start; i >= end; i--) {
                result.push(this.at(i * inc));
            }
        }
        result.push(b);
        return result;
    }
    calcSegmentTs(aT, bT, reversed, includeFirst) {
        assert(this.tIncrement, "tIncrement not defined on " + this);
        const inc = this.tIncrement;
        const result = [];
        if (includeFirst)
            result.push(aT);
        assert(reversed != aT < bT);
        if (aT < bT) {
            const start = Math.ceil((aT + NLA_PRECISION) / inc);
            const end = Math.floor((bT - NLA_PRECISION) / inc);
            for (let i = start; i <= end; i++) {
                result.push(i * inc);
            }
        }
        else {
            const start = Math.floor((aT - NLA_PRECISION) / inc);
            const end = Math.ceil((bT + NLA_PRECISION) / inc);
            for (let i = start; i >= end; i--) {
                result.push(i * inc);
            }
        }
        result.push(bT);
        return result;
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
        t = clamp(t, tStart, tEnd);
        return this.at(t).distanceTo(p);
    }
    /**
     * Behavior when curves are colinear: self intersections
     */
    isInfosWithCurve(curve) {
        if (insideIsInfosWithCurve) {
            return Curve.ispsRecursive(this, this.tMin, this.tMax, curve, curve.tMin, curve.tMax);
        }
        else {
            try {
                insideIsInfosWithCurve = true;
                const infos = curve.isInfosWithCurve(this);
                return infos.map((info) => {
                    assert(info);
                    const { tThis, tOther, p } = info;
                    return { tOther: tThis, tThis: tOther, p };
                });
            }
            finally {
                insideIsInfosWithCurve = false;
            }
        }
    }
    isTsWithSurface(surface) {
        if (surface instanceof PlaneSurface) {
            return this.isTsWithPlane(surface.plane);
        }
        if (surface instanceof ProjectedCurveSurface) {
            const projPlane = new P3(surface.dir.unit(), 0);
            const projThis = this.project(projPlane);
            const projEllipse = surface.baseCurve.project(projPlane);
            return projEllipse.isInfosWithCurve(projThis).map((info) => info.tOther);
        }
        if (surface instanceof EllipsoidSurface) {
            const thisOC = this.transform(surface.matrixInverse);
            if (!thisOC.getAABB().touchesAABBfuzzy(new AABB(V3.XYZ.negated(), V3.XYZ))) {
                return [];
            }
            const f = (t) => thisOC.at(t).length() - 1;
            const df = (t) => thisOC.at(t).unit().dot(thisOC.tangentAt(t));
            const stepSize = 1 / (1 << 11);
            const result = [];
            for (let startT = this.tMin; startT <= this.tMax; startT += stepSize) {
                const dt = stepSize * thisOC.tangentAt(startT).length();
                if (abs(f(startT)) <= dt) {
                    //const t = newtonIterate1d(f, startT, 16)
                    let t = newtonIterateWithDerivative(f, startT, 16, df);
                    if (!eq0(f(t)) || eq0(df(t))) {
                        t = newtonIterate1d(df, startT, 16);
                        //if (f(a) * f(b) < 0) {
                        //    t = bisect(f, a, b, 16)
                        //} else if (df(a) * df(b) < 0) {
                        //    t = bisect(df, a, b, 16)
                        //}
                    }
                    if (eq0(f(t)) && !result.some((r) => eq(r, t))) {
                        result.push(t);
                    }
                }
            }
            return result.filter((t) => surface.containsPoint(this.at(t)));
        }
        throw new Error();
    }
    arcLength(startT, endT, steps = 1) {
        assert(startT < endT, "startT < endT");
        return glqInSteps((t) => this.tangentAt(t).length(), startT, endT, steps);
    }
    equals(obj) {
        if (this === obj)
            return true;
        return (hasConstructor(obj, this.constructor) &&
            this.getConstructorParameters().equals(obj.getConstructorParameters()));
    }
    hashCode() {
        return this.getConstructorParameters().hashCode();
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
        return new AABB(V3.fromArray(mins), V3.fromArray(maxs));
    }
    reversed() {
        throw new Error();
    }
    clipPlane(plane) {
        const ists = this.isTsWithPlane(plane).filter((ist) => this.tMin <= ist && ist <= this.tMax);
        return getIntervals(ists, this.tMin, this.tMax).mapFilter(([a, b]) => {
            const midT = (a + b) / 2;
            return (!eq(a, b) &&
                plane.distanceToPointSigned(this.at(midT)) < 0 &&
                this.withBounds(a, b));
        });
    }
}
Curve.hlol = 0;
function mkcurves(implicitCurve, sStart, tStart, stepSize, bounds, validUV) {
    const start = V(sStart, tStart);
    assert(stepSize > 0);
    // checkDerivate(s => implicitCurve(s, 0), s => didu(s, 0), -1, 1, 0)
    // checkDerivate(t => implicitCurve(0, t), t => didv(0, t), -1, 1, 0)
    const { points, tangents } = followAlgorithm2d(implicitCurve, start, stepSize, bounds, validUV);
    if (points.length > 4 && points[0].distanceTo(points.last) <= abs(stepSize)) {
        // this is a loop: split it
        for (let i = 0; i < points.length - 1; i++) {
            assert(!points[i].equals(points[i + 1]));
        }
        const half = floor(points.length / 2);
        const points1 = points.slice(0, half), points2 = points.slice(half - 1, points.length);
        const tangents1 = tangents.slice(0, half), tangents2 = tangents.slice(half - 1, tangents.length);
        //tangents2[tangents2.length - 1] = tangents1[0]
        //points2[tangents2.length - 1] = points1[0]
        for (let i = 0; i < points1.length - 1; i++) {
            assert(!points1[i].equals(points1[i + 1]));
        }
        for (let i = 0; i < points2.length - 1; i++) {
            assert(!points2[i].equals(points2[i + 1]));
        }
        return [
            { points: points1, tangents: tangents1 },
            { points: points2, tangents: tangents2 },
        ];
    }
    else {
        // not a loop: check in the other direction
        const { points: reversePoints, tangents: reverseTangents, } = followAlgorithm2d(implicitCurve, start, -stepSize, bounds, validUV);
        const result = followAlgorithm2d(implicitCurve, reversePoints.last, stepSize, bounds, validUV, undefined, reverseTangents.last.negated());
        assert(result.points.length > 2);
        return [result];
    }
}
function breakDownPPCurves(ps1, ps2, uStep, vStep, stepSize) {
    const { uMin, uMax, vMin, vMax } = ps1;
    const bounds = uvInAABB2.bind(undefined, ps1);
    const bounds2 = uvInAABB2.bind(undefined, ps2);
    const deltaU = uMax - uMin, deltaV = vMax - vMin;
    const sRes = ceil(deltaU / uStep), tRes = ceil(deltaV / vStep);
    const grid = new Array(sRes * tRes).fill(0);
    //const printGrid = () => console.log(arrayFromFunction(tRes, i => grid.slice(sRes * i, sRes * (i + 1)).map(v => v ? 'X' : '_').join('')).join('\n'))
    const at = (i, j) => grid[j * sRes + i];
    const set = (i, j) => 0 <= i && i < sRes && 0 <= j && j < tRes && (grid[j * sRes + i] = 1);
    const result = [];
    const logTable = [];
    for (let i = 0; i < sRes; i++) {
        search: for (let j = 0; j < tRes; j++) {
            if (at(i, j))
                continue;
            set(i, j);
            const startU = uMin + (i + 0.5) * uStep, startV = vMin + (j + 0.5) * vStep;
            // assume point is valid, currently (TODO)
            const curvePointPPResult = curvePointPP(ps1, ps2, ps1.pUV(startU, startV));
            if (undefined === curvePointPPResult) {
                continue search;
            }
            const { p: startP, st1: { x: u, y: v }, st2: { x: u2, y: v2 }, } = curvePointPPResult;
            const li = floor((u - uMin) / uStep), lj = floor((v - vMin) / vStep);
            logTable.push({
                i,
                j,
                li,
                lj,
                startU,
                startV,
                u,
                v,
                "bounds(u, v)": bounds(u, v),
            });
            if (!(i == li && j == lj) && at(li, lj)) {
                continue search;
            }
            set(li, lj);
            // u, v are now good starting coordinates to use follow algorithm
            if (bounds(u, v) && bounds2(u2, v2)) {
                console.log(V(u, v).sce);
                const subResult = mkPPCurves(ps1, ps2, startP, stepSize, bounds, bounds2);
                for (const curveData of subResult) {
                    assert(curveData.st1s.length > 2);
                    for (const { x, y } of curveData.st1s) {
                        const lif = (x - uMin) / uStep, ljf = (y - vMin) / vStep;
                        set((lif - 0.5) | 0, (ljf - 0.5) | 0);
                        set((lif - 0.5) | 0, (ljf + 0.5) | 0);
                        set((lif + 0.5) | 0, (ljf - 0.5) | 0);
                        set((lif + 0.5) | 0, (ljf + 0.5) | 0);
                    }
                }
                //printGrid()
                result.push(...subResult);
            }
        }
    }
    console.table(logTable);
    for (const { points } of result) {
        for (let i = 0; i < points.length - 1; i++) {
            assert(!points[i].equals(points[i + 1]));
        }
    }
    return result.map(({ points, tangents, st1s }) => {
        return new PPCurve(points, tangents, ps1, ps2, st1s, undefined, stepSize, 1);
    });
}
function mkPPCurves(ps1, ps2, startPoint, stepSize, bounds1, bounds2) {
    // checkDerivate(s => implicitCurve(s, 0), s => didu(s, 0), -1, 1, 0)
    // checkDerivate(t => implicitCurve(0, t), t => didv(0, t), -1, 1, 0)
    const { points, tangents, st1s } = followAlgorithmPP(ps1, ps2, startPoint, stepSize, bounds1, bounds2);
    if (points[0].distanceTo(points.last) < stepSize && points.length > 2) {
        // this is a loop: split it
        for (let i = 0; i < points.length - 1; i++) {
            assert(!points[i].equals(points[i + 1]));
        }
        const half = floor(points.length / 2);
        const points1 = points.slice(0, half), points2 = points.slice(half - 1, points.length);
        const tangents1 = tangents.slice(0, half), tangents2 = tangents.slice(half - 1, tangents.length);
        const st1s1 = st1s.slice(0, half), st1s2 = st1s.slice(half - 1, tangents.length);
        tangents2[tangents2.length - 1] = tangents1[0];
        points2[tangents2.length - 1] = points1[0];
        st1s2[tangents2.length - 1] = st1s1[0];
        for (let i = 0; i < points1.length - 1; i++) {
            assert(!points1[i].equals(points1[i + 1]));
        }
        for (let i = 0; i < points2.length - 1; i++) {
            assert(!points2[i].equals(points2[i + 1]));
        }
        return [
            { points: points1, tangents: tangents1, st1s: st1s1 },
            { points: points2, tangents: tangents2, st1s: st1s2 },
        ];
    }
    else {
        // not a loop: check in the other direction
        const { points: reversePoints } = followAlgorithmPP(ps1, ps2, startPoint, -stepSize, bounds1, bounds2);
        const result = followAlgorithmPP(ps1, ps2, reversePoints.last, stepSize, bounds1, bounds2);
        assert(result.points.length > 2);
        return [result];
    }
}
function AABB2(uMin, uMax, vMin, vMax) {
    return { uMin, uMax, vMin, vMax };
}
function uvInAABB2(aabb2, u, v) {
    return (aabb2.uMin <= u && u <= aabb2.uMax && aabb2.vMin <= v && v <= aabb2.vMax);
}
function curvePoint(implicitCurve, startPoint, didu, didv) {
    let p = startPoint;
    for (let i = 0; i < 8; i++) {
        const fp = implicitCurve(p.x, p.y);
        const dfpdx = didu(p.x, p.y), dfpdy = didv(p.x, p.y);
        const scale = fp / (dfpdx * dfpdx + dfpdy * dfpdy);
        p = p.minus(new V3(scale * dfpdx, scale * dfpdy, 0));
    }
    return p;
}
function curvePointMF(mf, startPoint, steps = 8, eps = 1 / (1 << 30)) {
    let p = startPoint;
    for (let i = 0; i < steps; i++) {
        const fp = mf(p.x, p.y);
        const dfpdx = mf.x(p.x, p.y), dfpdy = mf.y(p.x, p.y);
        const scale = fp / (dfpdx * dfpdx + dfpdy * dfpdy);
        p = p.minus(new V3(scale * dfpdx, scale * dfpdy, 0));
        if (abs(fp) <= eps)
            break;
    }
    return p;
}

class XiEtaCurve extends Curve {
    constructor(center, f1, f2, tMin, tMax) {
        super(tMin, tMax);
        this.center = center;
        this.f1 = f1;
        this.f2 = f2;
        this.tMin = tMin;
        this.tMax = tMax;
        assertVectors(center, f1, f2);
        this.normal = f1.cross(f2);
        if (!this.normal.likeO()) {
            this.normal = this.normal.unit();
            this.matrix = M4.forSys(f1, f2, this.normal, center);
            this.matrixInverse = this.matrix.inversed();
        }
        else {
            this.matrix = M4.forSys(f1, f2, f1.unit(), center);
            const f1p = f1.getPerpendicular();
            // prettier-ignore
            this.matrixInverse = new M4(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1).times(M4.forSys(f1, f1p, f1.cross(f1p), center).inversed());
        }
    }
    /**
     * Intersection of the unit curve with the line ax + by = c.
     */
    static intersectionUnitLine(a, b, c, tMin, tMax) {
        throw new Error("abstract");
    }
    /**
     * Returns a new EllipseCurve representing an ellipse parallel to the XY-plane
     * with semi-major/minor axes parallel t the X and Y axes.
     *
     * @param a length of the axis parallel to X axis.
     * @param b length of the axis parallel to Y axis.
     * @param center center of the ellipse.
     */
    static forAB(a, b, center = V3.O) {
        return new this(center, V(a, 0, 0), V(0, b, 0));
    }
    static XYLCValid(pLC) {
        throw new Error("abstract");
    }
    static XYLCPointT(pLC, tMin, tMax) {
        throw new Error("abstract");
    }
    static unitIsInfosWithLine(anchorLC, dirLC, anchorWC, dirWC, tMin, tMax) {
        throw new Error("abstract");
    }
    addToMesh(mesh, res = 4, radius = 0, pointStep = 1) {
        const baseNormals = arrayFromFunction(res, (i) => V3.polar(1, (TAU * i) / res));
        const baseVertices = arrayFromFunction(res, (i) => V3.polar(radius, (TAU * i) / res));
        const inc = this.tIncrement;
        const start = Math.ceil((this.tMin + NLA_PRECISION) / inc);
        const end = Math.floor((this.tMax - NLA_PRECISION) / inc);
        for (let i = start; i <= end; i += pointStep) {
            const t = i * inc;
            const start = mesh.vertices.length;
            if (0 !== i) {
                for (let j = 0; j < res; j++) {
                    pushQuad(mesh.TRIANGLES, true, start - res + j, start + j, start - res + ((j + 1) % res), start + ((j + 1) % res));
                }
            }
            const point = this.at(t), tangent = this.tangentAt(t);
            const matrix = M4.forSys(this.normal, tangent.cross(this.normal), tangent, point);
            mesh.normals.push(...matrix.transformedVectors(baseNormals));
            mesh.vertices.push(...matrix.transformedPoints(baseVertices));
        }
    }
    getConstructorParameters() {
        return [this.center, this.f1, this.f2];
    }
    isInfosWithCurve(curve) {
        if (curve instanceof L3) {
            return this.isInfosWithLine(curve.anchor, curve.dir1, this.tMin, this.tMax, curve.tMin, curve.tMax);
        }
        if (curve instanceof BezierCurve) {
            return this.isInfosWithBezier(curve);
        }
        if (curve instanceof XiEtaCurve) {
            if (!this.normal.isParallelTo(curve.normal)) {
                return mapFilter(this.isTsWithPlane(curve.getPlane()), (tThis) => {
                    const p = this.at(tThis);
                    if (curve.containsPoint(p)) {
                        return { tThis, tOther: curve.pointT(p), p };
                    }
                    return undefined;
                });
            }
        }
        return super.isInfosWithCurve(curve);
    }
    transform(m4) {
        return new this.constructor(m4.transformPoint(this.center), m4.transformVector(this.f1), m4.transformVector(this.f2), this.tMin, this.tMax);
    }
    equals(obj) {
        return (this == obj ||
            (undefined != obj &&
                this.constructor == obj.constructor &&
                this.center.equals(obj.center) &&
                this.f1.equals(obj.f1) &&
                this.f2.equals(obj.f2)));
    }
    hashCode() {
        let hashCode = 0;
        hashCode = hashCode * 31 + this.center.hashCode();
        hashCode = hashCode * 31 + this.f1.hashCode();
        hashCode = hashCode * 31 + this.f2.hashCode();
        return hashCode | 0;
    }
    likeCurve(curve) {
        return (hasConstructor(curve, this.constructor) &&
            this.center.like(curve.center) &&
            this.f1.like(curve.f1) &&
            this.f2.like(curve.f2));
    }
    normalP(t) {
        return this.tangentAt(t).cross(this.normal);
    }
    getPlane() {
        return P3.normalOnAnchor(this.normal, this.center);
    }
    isTsWithPlane(planeWC) {
        assertInst(P3, planeWC);
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
        if (planeWC.normal1.isParallelTo(this.normal)) {
            return [];
        }
        const n = planeWC.normal1, w = planeWC.w, center = this.center, f1 = this.f1, f2 = this.f2, g1 = n.dot(f1), g2 = n.dot(f2), g3 = w - n.dot(center);
        return this.constructor.intersectionUnitLine(g1, g2, g3, this.tMin, this.tMax);
    }
    pointT(p) {
        assertVectors(p);
        const pLC = this.matrixInverse.transformPoint(p);
        return this.constructor.XYLCPointT(pLC);
    }
    containsPoint(p) {
        const pLC = this.matrixInverse.transformPoint(p);
        return (eq0(pLC.z) &&
            this.isValidT(this.constructor.XYLCPointT(pLC, this.tMin, this.tMax)));
    }
    isInfosWithLine(anchorWC, dirWC, tMin = this.tMin, tMax = this.tMax, lineMin = -100000, lineMax = 100000) {
        const anchorLC = this.matrixInverse.transformPoint(anchorWC);
        const dirLC = this.matrixInverse.transformVector(dirWC);
        if (eq0(dirLC.z)) {
            // local line parallel to XY-plane
            if (eq0(anchorLC.z)) {
                // local line lies in XY-plane
                return this.constructor.unitIsInfosWithLine(anchorLC, dirLC, anchorWC, dirWC, tMin, tMax);
            }
        }
        else {
            // if the line intersects the XY-plane in a single point, there can be an intersection there
            // find point, then check if distance from circle = 1
            const otherTAtZ0 = anchorLC.z / dirLC.z;
            const isp = dirLC.times(otherTAtZ0).plus(anchorLC);
            if (this.constructor.XYLCValid(isp)) {
                // point lies on unit circle
                return [
                    {
                        tThis: this.constructor.XYLCPointT(isp),
                        tOther: otherTAtZ0,
                        p: anchorWC.plus(dirWC.times(otherTAtZ0)),
                    },
                ];
            }
        }
        return [];
    }
    isTsWithSurface(surface) {
        if (surface instanceof PlaneSurface) {
            return this.isTsWithPlane(surface.plane);
        }
        else if (surface instanceof EllipsoidSurface) {
            const isEllipses = surface.isCurvesWithPlane(this.getPlane());
            return isEllipses
                .flatMap((isEllipse) => this.isInfosWithCurve(isEllipse))
                .filter((info) => surface.containsPoint(info.p))
                .map((info) => info.tThis);
        }
        else if (surface instanceof ProjectedCurveSurface ||
            surface instanceof ConicSurface) {
            return surface
                .isCurvesWithPlane(this.getPlane())
                .flatMap((curve) => this.isInfosWithCurve(curve))
                .map((info) => info.tThis);
        }
        else {
            throw new Error();
        }
    }
    isInfosWithBezier(bezierWC) {
        const bezierLC = bezierWC.transform(this.matrixInverse);
        if (new PlaneSurface(P3.XY).containsCurve(bezierLC)) {
            return this.isInfosWithBezier2D(bezierWC);
        }
        else {
            const infos = mapFilter(bezierLC.isTsWithPlane(P3.XY), (tOther) => {
                const pLC = bezierLC.at(tOther);
                if (this.constructor.XYLCValid(pLC)) {
                    return {
                        tOther: tOther,
                        p: bezierWC.at(tOther),
                        tThis: this.constructor.XYLCPointT(pLC),
                    };
                }
                return undefined;
            });
            return infos;
        }
    }
    isInfosWithBezier2D(bezierWC, sMin = bezierWC.tMin, sMax = bezierWC.tMax) {
        return Curve.ispsRecursive(this, this.tMin, this.tMax, bezierWC, sMin, sMax);
    }
    isOrthogonal() {
        return this.f1.isPerpendicularTo(this.f2);
    }
    at2(xi, eta) {
        assertNumbers(xi, eta);
        // center + f1 xi + f2 eta
        return this.center.plus(this.f1.times(xi)).plus(this.f2.times(eta));
    }
    debugInfo() {
        return {
            points: [
                this.center,
                this.at2(0.5, 0),
                this.at2(0, 1 / 3),
                this.at2(0, 2 / 3),
            ],
            lines: [this.center, this.at2(0, 1), this.center, this.at2(1, 0)],
        };
    }
}
/**
 * Transforms the unit 4d parabola P(t) = t² (0, 1, 0, 0) + t (1, 0, 0, 0) + (0, 0, 0, 1) using m and projects the
 * result into 3d. This is used for the transform4 implementation of conics. The parabola may not cross the vanishing
 * plane of m in the interval [tMin, tMax], as that would result in discontinuities.
 */
function parabola4Projection(m, tMin, tMax) {
    return HyperbolaCurve.XY.rotateZ(45 * DEG);
}

class ImplicitCurve extends Curve {
    constructor(points, tangents, dir = 1, generator, tMin = 1 == dir ? 0 : -(points.length - 1), tMax = 1 == dir ? points.length - 1 : 0) {
        super(tMin, tMax);
        this.points = points;
        this.tangents = tangents;
        this.dir = dir;
        this.generator = generator;
        assert(points.length > 2);
        assert(0 <= tMin && tMin <= points.length - 1, tMin, points.length);
        assert(0 <= tMax && tMax <= points.length - 1, tMax, points.length);
    }
    likeCurve(curve) {
        throw new Error("Method not implemented.");
    }
    toSource(rounder = (x) => x) {
        return this.generator || super.toSource(rounder);
    }
    containsPoint(p) {
        assertVectors(p);
        return !isNaN(this.pointT(p));
    }
    equals(obj) {
        return (this == obj ||
            (Object.getPrototypeOf(obj) == PICurve.prototype &&
                this.points[0].equals(obj.points[0]) &&
                this.tangents[0].equals(obj.tangents[0])));
    }
    hashCode() {
        return [this.points[0], this.tangents[0]].hashCode();
    }
    tangentP(pWC) {
        assertVectors(pWC);
        assert(this.containsPoint(pWC), "this.containsPoint(pWC)" + this.containsPoint(pWC));
        const t = this.pointT(pWC);
        return this.tangentAt(t);
    }
    tangentAt(t) {
        t = clamp(t, this.tMin, this.tMax);
        return V3.lerp(this.tangents[floor(t)], this.tangents[ceil(t)], t % 1);
    }
    at(t) {
        assert(isFinite(t));
        return V3.lerp(this.points[floor(t)], this.points[ceil(t)], t % 1);
    }
    getConstructorParameters() {
        throw new Error();
    }
    roots() {
        const allTs = arrayRange(0, this.points.length);
        return [allTs, allTs, allTs];
    }
    /**
     * @param mesh
     * @param res
     * @param radius default to 0. Use the shader to achieve dynamic scaling.
     * @param pointStep
     */
    addToMesh(mesh, res = 4, radius = 0, pointStep = 1) {
        const baseNormals = arrayFromFunction(res, (i) => V3.polar(1, (TAU * i) / res));
        const baseVertices = arrayFromFunction(res, (i) => V3.polar(radius, (TAU * i) / res));
        let prevTangent = V3.Z, prevMatrix = M4.IDENTITY;
        for (let i = 0; i < this.points.length; i += pointStep) {
            const start = mesh.vertices.length;
            if (0 !== i) {
                for (let j = 0; j < res; j++) {
                    pushQuad(mesh.TRIANGLES, true, start - res + j, start + j, start - res + ((j + 1) % res), start + ((j + 1) % res));
                }
            }
            const point = this.points[i], tangent = this.tangents[i];
            const tangentMatrix = M4.rotateAB(prevTangent, tangent).times(prevMatrix);
            mesh.normals.push(...tangentMatrix.transformedVectors(baseNormals));
            const baseMatrix = M4.translate(point).times(tangentMatrix);
            mesh.vertices.push(...baseMatrix.transformedPoints(baseVertices));
            prevTangent = tangent;
            prevMatrix = tangentMatrix;
        }
    }
    rootsApprox() {
        const roots = [[], [], []];
        const points = this.points;
        let lastDiff = points[1].minus(points[0]);
        for (let i = 2; i < points.length; i++) {
            const diff = points[i].minus(points[i - 1]);
            for (let dim = 0; dim < 3; dim++) {
                if (Math.sign(lastDiff.e(dim)) != Math.sign(diff.e(dim))) {
                    roots[dim].push(i);
                }
            }
            lastDiff = diff;
        }
        return roots;
    }
    pointT(pWC) {
        const startT = arrayRange(floor(this.tMin), ceil(this.tMax), 1).withMax((t) => -pWC.distanceTo(this.points[t]));
        if (undefined === startT)
            throw new Error();
        if (this.points[startT].like(pWC))
            return startT;
        const a = max(0, startT - 1), b = min(this.points.length - 1, startT + 1);
        const tangent = this.tangentAt(startT);
        const f = (t) => this.at(t).to(pWC).dot(tangent);
        // const df = (t: number) => -this.tangentAt(clamp(t, 0, this.points.length - 1)).dot(tangent)
        //checkDerivate(f, df, 0, this.points.length - 2, 3)
        const t = bisect(f, a, b, 32);
        if (!isFinite(t) || !eq0(this.at(t).distanceTo(pWC))) {
            return NaN;
        }
        return t;
    }
}
ImplicitCurve.prototype.tIncrement = 1;
/**
 * isInfosWithLine for an ImplicitCurve defined as the intersection of two surfaces.
 */
function surfaceIsICurveIsInfosWithLine(surface1, surface2, anchorWC, dirWC, tMin, tMax, lineMin, lineMax) {
    const line = new L3(anchorWC, dirWC.unit());
    const psTs = surface1.isTsForLine(line);
    const isTs = surface2.isTsForLine(line);
    const commonTs = psTs.filter((psT) => isTs.some((isT) => eq(psT, isT)));
    const commonTInfos = commonTs.map((t) => ({
        tThis: 0,
        tOther: t / dirWC.length(),
        p: line.at(t),
    }));
    const result = commonTInfos.filter((info) => this.containsPoint(info.p));
    result.forEach((info) => (info.tThis = this.pointT(info.p)));
}

/**
 * Bezier curve with degree 3.
 */
class BezierCurve extends Curve {
    constructor(p0, p1, p2, p3, tMin = -0.1, tMax = 1.1) {
        super(tMin, tMax);
        assertVectors(p0, p1, p2, p3);
        assert(isFinite(tMin) && isFinite(tMax));
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
        return new BezierCurve(V(0, p0y), V(1 / 3, p1y), V(2 / 3, p2y), V(1, p3y), tMin, tMax);
    }
    static quadratic(a, b, c, tMin = 0, tMax = 1) {
        const line = L3.throughPoints(a, c);
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
        const f = (4 / 3) * Math.tan(phi / 4);
        return new BezierCurve(V3.X, new V3(1, f, 0), new V3(cos(phi) + f * sin(phi), sin(phi) - f * cos(phi), 0), V3.sphere(phi, 0), 0, 1);
    }
    getConstructorParameters() {
        return [this.p0, this.p1, this.p2, this.p3];
    }
    at(t) {
        // = s^3 p0 + 3 s^2 t p1 + 3 s t^2 p2 + t^3 p3
        assertNumbers(t);
        const p0 = this.p0, p1 = this.p1, p2 = this.p2, p3 = this.p3;
        const s = 1 - t, c0 = s * s * s, c1 = 3 * s * s * t, c2 = 3 * s * t * t, c3 = t * t * t;
        return new V3(p0.x * c0 + p1.x * c1 + p2.x * c2 + p3.x * c3, p0.y * c0 + p1.y * c1 + p2.y * c2 + p3.y * c3, p0.z * c0 + p1.z * c1 + p2.z * c2 + p3.z * c3);
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
        assertNumbers(t);
        const p0 = this.p0, p1 = this.p1, p2 = this.p2, p3 = this.p3;
        const s = 1 - t, c01 = 3 * s * s, c12 = 6 * s * t, c23 = 3 * t * t;
        return new V3((p1.x - p0.x) * c01 + (p2.x - p1.x) * c12 + (p3.x - p2.x) * c23, (p1.y - p0.y) * c01 + (p2.y - p1.y) * c12 + (p3.y - p2.y) * c23, (p1.z - p0.z) * c01 + (p2.z - p1.z) * c12 + (p3.z - p2.z) * c23);
    }
    ddt(t) {
        assertNumbers(t);
        const p0 = this.p0, p1 = this.p1, p2 = this.p2, p3 = this.p3;
        const c012 = 6 * (1 - t), c123 = 6 * t;
        return new V3((p2.x - 2 * p1.x + p0.x) * c012 + (p3.x - 2 * p2.x + p1.x) * c123, (p2.y - 2 * p1.y + p0.y) * c012 + (p3.y - 2 * p2.y + p1.y) * c123, (p2.z - 2 * p1.z + p0.z) * c012 + (p3.z - 2 * p2.z + p1.z) * c123);
    }
    normalP(t) {
        const tangent = this.tangentAt(t);
        const rot = tangent.cross(this.ddt(t));
        return rot.cross(tangent);
    }
    isTsWithPlane(planeWC) {
        assertInst(P3, planeWC);
        /*
             We are solving for t:
             n := plane.normal1
             this.at(t) DOT n == plane.w // according to plane definition
             (a t³ + b t² + c t + d) DOT n == plane.w // bezier curve as cubic equation
             (a DOT n) t³ + (b DOT n) t³ + (c DOT n) t + d DOT n - plane.w == 0 // multiply out DOT n, minus plane.w
             */
        const { p0, p1, p2, p3 } = this;
        const n = planeWC.normal1;
        const a = p1.minus(p2).times(3).minus(p0).plus(p3);
        const b = p0.plus(p2).times(3).minus(p1.times(6));
        const c = p1.minus(p0).times(3);
        const d = p0;
        return solveCubicReal2(a.dot(n), b.dot(n), c.dot(n), d.dot(n) - planeWC.w).filter((t) => between(t, this.tMin, this.tMax));
    }
    isTsWithSurface(surfaceWC) {
        if (surfaceWC instanceof CylinderSurface) {
            const projPlane = new P3(surfaceWC.dir.unit(), 0);
            const projThis = this.project(projPlane);
            const projEllipse = surfaceWC.baseCurve.project(projPlane);
            return projEllipse
                .isInfosWithBezier2D(projThis)
                .map((info) => info.tOther);
        }
        return super.isTsWithSurface(surfaceWC);
    }
    likeCurve(curve) {
        return (this == curve ||
            (hasConstructor(curve, BezierCurve) &&
                this.p0.like(curve.p0) &&
                this.p1.like(curve.p1) &&
                this.p2.like(curve.p2) &&
                this.p3.like(curve.p3)));
    }
    equals(obj) {
        return (this == obj ||
            (hasConstructor(obj, BezierCurve) &&
                this.p0.equals(obj.p0) &&
                this.p1.equals(obj.p1) &&
                this.p2.equals(obj.p2) &&
                this.p3.equals(obj.p3)));
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
        if (isNaN((curveP0T = this.pointT(curve.p0))) ||
            isNaN((curveP3T = this.pointT(curve.p3)))) {
            return false;
        }
        let thisSplit;
        if (eq(1, curveP0T)) {
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
    selectPart(t0, t1) {
        const t1Adjusted = (t1 - t0) / (1 - t0);
        return this.split(t0)[1].split(t1Adjusted)[0];
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
        return [V3.O, a, b, c];
    }
    pointT2(p, tMin = this.tMin, tMax = this.tMax) {
        const t = this.closestTToPoint(p, undefined, tMin, tMax);
        assert(this.at(t).like(p));
        return t;
    }
    pointT(p) {
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
        const maxDim = NLA_PRECISION < a.maxAbsElement()
            ? a.maxAbsDim()
            : NLA_PRECISION < b.maxAbsElement()
                ? b.maxAbsDim()
                : NLA_PRECISION < c.maxAbsElement()
                    ? c.maxAbsDim()
                    : assertNever();
        const results = solveCubicReal2(a.e(maxDim), b.e(maxDim), c.e(maxDim), d.e(maxDim)).filter((t) => this.at(t).like(p));
        if (0 == results.length)
            return NaN;
        if (1 == results.length)
            return results[0];
        throw new Error("multiple intersection " + this.toString() + p.sce);
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
            if (eq0(a[dim]) && eq0(b[dim]) && eq0(c[dim])) {
                // for case x:
                // ax == bx == cx == 0 => x(t) = dx
                // x value is constant
                // if x == 0 for all t, this does not limit the result, otherwise, there is no result, i.e
                // the passed point is not on the curve
                if (!eq0(d[dim]))
                    return NaN;
            }
            else {
                const newResults = solveCubicReal2(a[dim], b[dim], c[dim], d[dim]);
                if (0 == newResults.length)
                    return NaN;
                if (1 == newResults.length)
                    return newResults[0];
                if (results) {
                    results = results.filter((t) => newResults.some((t2) => eq(t, t2)));
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
        throw new Error("multiple intersection " + results + this.toString() + p.sce);
    }
    transform(m4) {
        // perspective projection turn bezier curve into rational spline
        assert(m4.isNoProj(), m4.str);
        return new BezierCurve(m4.transformPoint(this.p0), m4.transformPoint(this.p1), m4.transformPoint(this.p2), m4.transformPoint(this.p3), this.tMin, this.tMax);
    }
    isClosed() {
        return this.p0.like(this.p3);
    }
    isQuadratic() {
        return this.p0.lerp(this.p1, 1.5).like(this.p3.lerp(this.p2, 1.5));
    }
    debugInfo() {
        return {
            lines: [0, 1, 1, 2, 2, 3].map((i) => this.points[i]),
            points: this.points,
        };
    }
    split(t) {
        // do de Casteljau's algorithm at t, the resulting points are the points needed to create 2 new curves
        const s = 1 - t;
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
        return [
            new BezierCurve(p0, b01, b02, b03),
            new BezierCurve(b03, b12, b21, p3),
        ];
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
        return arrayFromFunction(3, (dim) => solveCubicReal2(0, a.e(dim), b.e(dim), c.e(dim)));
    }
    isInfosWithLine(anchorWC, dirWC, tMin, tMax, lineMin = -100000, lineMax = 100000) {
        // const dirLength = dirWC.length()
        // // TODO: no:
        // let result = Curve.ispsRecursive(this, this.tMin, this.tMax, new L3(anchorWC, dirWC.unit()), lineMin, lineMax)
        // result = fuzzyUniquesF(result, info => info.tOther)
        // result.forEach(info => (info.tOther /= dirLength))
        // return result
        // looking for this.at(t) == line.at(s)
        // this.at(t).x == anchorWC.x + dirWC.x * s
        // (this.at(t).x - anchorWC.x) / dirWC.x == s (analogue for y and z) (1x, 1y, 1z)
        // (1x) - (1y):
        // (this.at(t).x - anchorWC.x) / dirWC.x - (this.at(t).y - anchorWC.y) / dirWC.y == 0
        // (this.at(t).x - anchorWC.x) * dirWC.y - (this.at(t).y - anchorWC.y) * dirWC.x == 0 (2)
        // cubic equation params (see #pointT):
        const { p0, p1, p2, p3 } = this;
        const a = p1.minus(p2).times(3).minus(p0).plus(p3);
        const v1 = V3.UNITS[a.minAbsDim()];
        const testPlane = P3.forAnchorAndPlaneVectors(anchorWC, dirWC, v1.isParallelTo(dirWC) ? a : v1);
        return this.isTsWithPlane(testPlane)
            .map((tThis) => {
            const p = this.at(tThis);
            return { tThis, tOther: L3.pointT(anchorWC, dirWC, p), p };
        })
            .filter((info) => L3.containsPoint(anchorWC, dirWC, info.p));
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
            return atT
                .minus(line.at(atT.dot(line.dir1) - anchorDotDir1))
                .dot(this.tangentAt(t));
        };
        const STEPS = 32;
        const startT = arrayFromFunction(STEPS, (i) => tMin + ((tMax - tMin) * i) / STEPS).withMax((t) => -f(t));
        return newtonIterate1d(f, startT, 8);
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
    isInfosWithBezier3(bezier, tMin, tMax, sMin, sMax) {
        const handleStartTS = (startT, startS) => {
            if (!result.some((info) => eq(info.tThis, startT) && eq(info.tOther, startS))) {
                const f1 = (t, s) => this.tangentAt(t).dot(this.at(t).minus(bezier.at(s)));
                const f2 = (t, s) => bezier.tangentAt(s).dot(this.at(t).minus(bezier.at(s)));
                // f = (b1, b2, t1, t2) = b1.tangentAt(t1).dot(b1.at(t1).minus(b2.at(t2)))
                const fdt1 = (b1, b2, t1, t2) => b1.ddt(t1).dot(b1.at(t1).minus(b2.at(t2))) +
                    b1.tangentAt(t1).squared();
                const fdt2 = (b1, b2, t1, t2) => -b1.tangentAt(t1).dot(b2.tangentAt(t2));
                const ni = newtonIterate2dWithDerivatives(f1, f2, startT, startS, 16, fdt1.bind(undefined, this, bezier), fdt2.bind(undefined, this, bezier), (t, s) => -fdt2(bezier, this, s, t), (t, s) => -fdt1(bezier, this, s, t));
                result.push({ tThis: ni.x, tOther: ni.y, p: this.at(ni.x) });
            }
        };
        tMin = undefined !== tMin ? tMin : this.tMin;
        tMax = undefined !== tMax ? tMax : this.tMax;
        sMin = undefined !== sMin ? sMin : bezier.tMin;
        sMax = undefined !== sMax ? sMax : bezier.tMax;
        // stack of indices:
        const indices = [tMin, tMax, sMin, sMax];
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
                    indices.push(tMin, tMid, sMin, sMid, tMin, tMid, sMid, sMax, tMid, tMax, sMin, sMid, tMid, tMax, sMid, sMax);
                }
            }
        }
        return result;
    }
    isInfosWithBezier(bezier, tMin, tMax, sMin, sMax) {
        tMin = undefined !== tMin ? tMin : this.tMin;
        tMax = undefined !== tMax ? tMax : this.tMax;
        sMin = undefined !== sMin ? sMin : bezier.tMin;
        sMax = undefined !== sMax ? sMax : bezier.tMax;
        assertf(() => tMin < tMax);
        assertf(() => sMin < sMax);
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
            const splits = fuzzyUniques(this.roots().concatenated().filter(isFinite).concat([tMin, tMax])).sort(MINUS);
            //const aabbs = arrayFromFunction(splits.length - 1, i => this.getAABB(splits[i], splits[i + 1]))
            Array.from(combinations(splits.length - 1)).forEach(({ i, j }) => {
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
        if (curve instanceof L3) {
            return this.isInfosWithLine(curve.anchor, curve.dir1, curve.tMin, curve.tMax);
        }
        if (curve instanceof BezierCurve) {
            return this.isInfosWithBezier(curve);
        }
        return curve
            .isInfosWithCurve(this)
            .map(({ tThis, tOther, p }) => ({ tThis: tOther, tOther: tThis, p }));
    }
    /**
     * Approximate this bezier curve with a number of circular segments. This curve is recursively split in half until
     * segments are close enough (relative error < REL_ERR in two test points) to an arc which goes through the start,
     * end and mid points of the segment.
     * @returns each EllipseCurve is circular and their tMin and tMax respectively define their start and end points.
     * @param t0 Start parameter of segment which should be approximated.
     * @param t1 End parameter of segment which should be approximated.
     * @param REL_ERROR max allowable relative error.
     * @param result Resulting circle arcs are stored in this array. Mainly used by the recursion.
     */
    circleApprox(t0 = this.tMin, t1 = this.tMax, REL_ERROR = 1 / 1024, result = []) {
        const a = this.at(t0), b = this.at(t1), tMid = (t0 + t1) / 2, pMid = this.at(tMid), abLine = L3.throughPoints(a, b);
        if (!abLine.containsPoint(pMid) &&
            between(abLine.pointT(pMid), 0, abLine.pointT(b))) {
            const arc = EllipseCurve.circleThroughPoints(a, pMid, b), arcRadius = arc.f1.length(), pTest1 = this.at(lerp(t0, t1, 0.25)), pTest2 = this.at(lerp(t0, t1, 0.75));
            if (abs(arc.center.distanceTo(pTest1) / arcRadius - 1) <= REL_ERROR &&
                abs(arc.center.distanceTo(pTest2) / arcRadius - 1) <= REL_ERROR) {
                result.push(arc);
                return result;
            }
        }
        this.circleApprox(t0, tMid, REL_ERROR, result);
        this.circleApprox(tMid, t1, REL_ERROR, result);
        return result;
    }
}
/**
 * https://en.wikipedia.org/wiki/Cubic_function#/media/File:Graph_of_cubic_polynomial.svg
 */
BezierCurve.EX2D = BezierCurve.graphXY(2, -3, -3, 2);
BezierCurve.EX3D = new BezierCurve(V3.O, V(-0.1, -1, 1), V(1.1, 1, 1), V3.X);
BezierCurve.QUARTER_CIRCLE = BezierCurve.approximateUnitArc(PI / 2);
BezierCurve.prototype.hlol = Curve.hlol++;
BezierCurve.prototype.tIncrement = 1 / 80;

/**
 * x² - y² = 1
 * C(t) = center + f1 * cosh(t) + f2 * sinh(t)
 */
class HyperbolaCurve extends XiEtaCurve {
    constructor(center, f1, f2, tMin = -7, tMax = 7) {
        super(center, f1, f2, tMin, tMax);
    }
    static XYLCValid(pLC) {
        return pLC.x > 0 && eq(1, pLC.x * pLC.x - pLC.y * pLC.y);
    }
    static XYLCPointT(pLC) {
        return Math.asinh(pLC.y);
    }
    /**
     * http://www.wolframalpha.com/input/?i=x%C2%BRep-y%C2%BRep%3D1,ax%2Bby%3Dc
     * Minor empiric test shows asinh(eta) consistently gets more accurate results than atanh(eta/xi)
     */
    static intersectionUnitLine(a, b, c) {
        if (eq0(b)) {
            const sqrtVal = snap0(Math.pow(c, 2) / Math.pow(a, 2) - 1);
            if (sqrtVal < 0 || c * a < 0) {
                return [];
            }
            else if (sqrtVal == 0) {
                return [0];
            }
            const eta1 = Math.sqrt(sqrtVal);
            return [-Math.asinh(eta1), Math.asinh(eta1)];
        }
        else if (eq(abs(a), abs(b))) {
            if (le(c * a, 0)) {
                return [];
            }
            const eta = (sign(a * b) * (Math.pow(c, 2) - Math.pow(a, 2))) / 2 / a / c;
            return [Math.asinh(eta)];
        }
        else {
            const sqrtVal = snap0(Math.pow(b, 2) * (-(Math.pow(a, 2)) + Math.pow(b, 2) + Math.pow(c, 2)));
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
        assertNumbers(t);
        // = center + f1 cosh t + f2 sinh t
        return this.center
            .plus(this.f1.times(Math.cosh(t)))
            .plus(this.f2.times(Math.sinh(t)));
    }
    toString() {
        return `${this.center} + ${this.f1} * cosh(t) + ${this.f2} * sinh(t)`;
    }
    tangentAt(t) {
        assertNumbers(t);
        // = f1 sinh t + f2 cosh t
        return this.f1.times(Math.sinh(t)).plus(this.f2.times(Math.cosh(t)));
    }
    tangentAt2(xi, eta) {
        assertNumbers(xi, eta);
        // = f1 eta + f2 xi
        return this.f1.times(eta).plus(this.f2.times(xi));
    }
    ddt(t) {
        assertNumbers(t);
        return this.f1.times(Math.cosh(t)).plus(this.f2.times(Math.sinh(t)));
    }
    isColinearTo(curve) {
        if (!hasConstructor(curve, HyperbolaCurve))
            return false;
        if (!curve.center || !this.center.like(curve.center)) {
            return false;
        }
        if (this === curve) {
            return true;
        }
        const { f1: f1, f2: f2 } = this.rightAngled(), { f1: c1, f2: c2 } = curve.rightAngled();
        return (eq(f1.squared(), Math.abs(f1.dot(c1))) &&
            eq(f2.squared(), Math.abs(f2.dot(c2))));
    }
    reversed() {
        return new HyperbolaCurve(this.center, this.f1, this.f2.negated(), -this.tMax, -this.tMin);
    }
    rightAngled() {
        const f1 = this.f1, f2 = this.f2, a = f1.dot(f2), b = f2.squared() + f1.squared();
        if (eq0(a)) {
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
        return Math.sqrt(1 + (b * b) / a / a);
    }
    roots() {
        // tangent(t) = f1 sinh t + f2 cosh t = 0
        // tangentAt2(xi, eta) = f1 eta + f2 xi = V3.O
        // xi² - eta² = 1 (by def for hyperbola)
        return arrayFromFunction(3, (dim) => {
            const a = this.f2.e(dim), b = this.f1.e(dim);
            return HyperbolaCurve.intersectionUnitLine(a, b, 0);
        });
    }
    transform4(m4) {
        const tMap = (t) => sign(t) * min(10, sqrt(-(1 - cosh(t)) / (1 + cosh(t))));
        // prettier-ignore
        const parabolaToUnitHyperbola = new M4(0, 1, 0, 1, 2, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 1);
        return parabola4Projection(M4.product(m4, this.matrix, parabolaToUnitHyperbola), tMap(this.tMin), tMap(this.tMax));
    }
}
HyperbolaCurve.XY = new HyperbolaCurve(V3.O, V3.X, V3.Y);
HyperbolaCurve.prototype.tIncrement = PI / 16;

/**
 * A 3-dimensional line. Defined by an anchor and a normalized direction vector.
 */
class L3 extends Curve {
    constructor(anchor, // line anchor
    dir1, // normalized line dir
    tMin = -4096, tMax = 4096) {
        super(tMin, tMax);
        this.anchor = anchor;
        this.dir1 = dir1;
        assertVectors(anchor, dir1);
        assert(dir1.hasLength(1), "dir must be unit" + dir1);
        assertf(() => !Number.isNaN(anchor.x));
    }
    isTsWithSurface(surface) {
        return surface.isTsForLine(this);
    }
    static throughPoints(anchor, b, tMin = 0, tMax) {
        const dir = b.minus(anchor);
        return new L3(anchor, dir.unit(), tMin, undefined !== tMax ? tMax : dir.length());
    }
    static anchorDirection(anchor, dir, min = 0, max = dir.length()) {
        const dir1 = dir.unit();
        return new L3(anchor, dir1, "number" == typeof min ? min : min.minus(anchor).dot(dir1), "number" == typeof max ? max : max.minus(anchor).dot(dir1));
    }
    static pointT(anchor, dir, x) {
        assertVectors(anchor, dir, x);
        return x.minus(anchor).dot(dir) / dir.squared();
    }
    static at(anchor, dir, t) {
        return anchor.plus(dir.times(t));
    }
    /**
     * Create new line which is the intersection of two planes. Throws error if planes are parallel.
     * @param plane1
     * @param plane2
     */
    static fromPlanes(plane1, plane2) {
        assertInst(P3, plane1, plane2);
        const dir = plane1.normal1.cross(plane2.normal1);
        const length = dir.length();
        if (length < 1e-10) {
            throw new Error("Parallel planes");
        }
        return plane1.intersectionWithPlane(plane2);
    }
    static containsPoint(anchor, dir, p) {
        const closestT = L3.pointT(anchor, dir, p);
        const distance = L3.at(anchor, dir, closestT).distanceTo(p);
        return eq0(distance);
    }
    roots() {
        return [[], [], []];
    }
    containsPoint(p) {
        assertVectors(p);
        const dist = this.distanceToPoint(p);
        assertNumbers(dist);
        return eq0(dist);
    }
    likeCurve(curve) {
        return (this == curve ||
            (hasConstructor(curve, L3) &&
                this.anchor.like(curve.anchor) &&
                this.dir1.like(curve.dir1)));
    }
    equals(obj) {
        return (this == obj ||
            (Object.getPrototypeOf(obj) == L3.prototype &&
                this.anchor.equals(obj.anchor) &&
                this.dir1.equals(obj.dir1)));
    }
    isColinearTo(obj) {
        return (obj instanceof L3 &&
            this.containsPoint(obj.anchor) &&
            eq(1, Math.abs(this.dir1.dot(obj.dir1))));
    }
    distanceToLine(line) {
        assertInst(L3, line);
        if (this.isParallelToLine(line)) {
            return this.distanceToPoint(line.anchor);
        }
        const dirCross1 = this.dir1.cross(line.dir1).unit();
        const anchorDiff = this.anchor.minus(line.anchor);
        return Math.abs(anchorDiff.dot(dirCross1));
    }
    distanceToPoint(x) {
        assertVectors(x);
        // See http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
        const t = x.minus(this.anchor).dot(this.dir1);
        return this.at(t).distanceTo(x);
        //return x.minus(this.anchor).cross(x.minus(this.anchor.plus(this.dir1))).length()
    }
    asSegmentDistanceToPoint(x, sStart, sEnd) {
        let t = x.minus(this.anchor).dot(this.dir1);
        t = clamp(t, sStart, sEnd);
        return this.at(t).minus(x).length();
    }
    asSegmentDistanceToLine(line, sStart, sEnd) {
        assertInst(L3, line);
        const dirCross = this.dir1.cross(line.dir1);
        const div = dirCross.squared();
        if (eq0(div)) {
            return undefined;
        } // lines parallel
        const anchorDiff = line.anchor.minus(this.anchor);
        // check if distance is zero (see also L3.distanceToLine)
        if (!eq0(anchorDiff.dot(dirCross.unit()))) {
            return undefined;
        }
        let t = this.infoClosestToLine(line).t;
        t = clamp(t, sStart, sEnd);
        return this.at(clamp(t, sStart, sEnd));
    }
    at(t) {
        assertNumbers(t);
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
        assertVectors(x);
        const t = x.minus(this.anchor).dot(this.dir1);
        return t;
    }
    /**
     * Returns true if the line is parallel (this.dir = line.dir || this.dir = -line.dir) to the argument.
     */
    isParallelToLine(line) {
        assertInst(L3, line);
        // we know that 1 == this.dir1.length() == line.dir1.length(), we can check for parallelity simpler than
        // isParallelTo()
        return eq(1, Math.abs(this.dir1.dot(line.dir1)));
    }
    angleToLine(line) {
        assertInst(L3, line);
        return this.dir1.angleTo(line.dir1);
    }
    /**
     *
     * @param line
     * @returns {boolean} If the distance between the lines is zero
     */
    intersectsLine(line) {
        return eq0(this.distanceToLine(line));
    }
    isInfosWithCurve(curve) {
        if (curve instanceof L3) {
            return this.isInfosWithLine(curve.anchor, curve.dir1);
        }
        return super.isInfosWithCurve(curve);
    }
    isInfosWithLine(anchorWC, dirWC) {
        const dirCross = this.dir1.cross(dirWC);
        const div = dirCross.squared();
        if (eq0(div)) {
            // lines are parallel
            return [];
        }
        const anchorDiff = anchorWC.minus(this.anchor);
        if (eq0(anchorDiff.dot(dirCross))) {
            const tThis = anchorDiff.cross(dirWC).dot(dirCross) / div;
            const tOther = anchorDiff.cross(this.dir1).dot(dirCross) / div;
            const p = this.at(tThis);
            return [{ tThis: tThis, tOther: tOther, p: p }];
        }
        return [];
    }
    isInfoWithLine(line) {
        // todo infos?
        assertInst(L3, line);
        const dirCross = this.dir1.cross(line.dir1);
        const div = dirCross.squared();
        if (eq0(div)) {
            return undefined;
        } // lines parallel
        const anchorDiff = line.anchor.minus(this.anchor);
        // check if distance is zero (see also L3.distanceToLine)
        if (!eq0(anchorDiff.dot(dirCross.unit()))) {
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
        assertInst(L3, line);
        const dirCross = this.dir1.cross(line.dir1);
        const div = dirCross.squared();
        const anchorDiff = line.anchor.minus(this.anchor);
        const s = anchorDiff.cross(this.dir1).dot(dirCross) / div;
        const t = anchorDiff.cross(line.dir1).dot(dirCross) / div;
        return { s: s, t: t };
        // console.log(segmentIntersectsRay, a, b, "ab", ab, "p", p, "dir", dir, s > 0 && t / div >= 0 && t / div <= 1,
        // "s", s, "t", t, "div", div)
    }
    ddt() {
        return V3.O;
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
        const bd = b.dot(d), bb = b.squared(), dd = d.squared(), ca = a.minus(c), divisor = bd * bd - dd * bb;
        const t = (ca.dot(b) * bd - ca.dot(d) * bb) / divisor;
        const s = (ca.dot(b) * dd - ca.dot(d) * bd) / divisor;
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
    tangentAt() {
        return this.dir1;
    }
    isTWithPlane(plane) {
        // plane: plane.normal1 * p = plane.w
        // line: p=line.point + lambda * line.dir1
        const div = plane.normal1.dot(this.dir1);
        if (eq0(div))
            return NaN;
        const lambda = (plane.w - plane.normal1.dot(this.anchor)) / div;
        return lambda;
    }
    reversed() {
        return new L3(this.anchor, this.dir1.negated(), -this.tMax, -this.tMin);
    }
    isTsWithPlane(planeWC) {
        const t = this.isTWithPlane(planeWC);
        return isNaN(t) ? [] : [t];
    }
    flipped() {
        return new L3(this.anchor, this.dir1.negated());
    }
    transform(m4) {
        const newAnchor = m4.transformPoint(this.anchor);
        const newDir = m4.transformVector(this.dir1);
        return new L3(newAnchor, newDir.unit(), this.tMin * newDir.length(), this.tMax * newDir.length());
    }
    transform4(m4) {
        const vanishingPlane = P3.vanishingPlane(m4);
        if (!vanishingPlane)
            return this.transform(m4);
        const pMin = this.at(this.tMin);
        const pMax = this.at(this.tMax);
        if (le(vanishingPlane.distanceToPointSigned(pMin), 0) ||
            le(vanishingPlane.distanceToPointSigned(pMax), 0)) {
            throw new Error("line must be in front of vanishingPlane in [tMin, tMax]");
        }
        const anchor = lt(0, vanishingPlane.distanceToPointSigned(this.anchor))
            ? this.anchor
            : this.at((this.tMin + this.tMax) / 2);
        const transformedAnchor = m4.timesVector(VV(anchor.x, anchor.y, anchor.z, 1));
        const transformedVector = m4.timesVector(VV(this.dir1.x, this.dir1.y, this.dir1.z, 0));
        const newDir = transformedVector
            .times(transformedAnchor.w)
            .minus(transformedAnchor.times(transformedVector.w))
            .V3();
        const newAnchor = transformedAnchor.p3();
        return L3.anchorDirection(newAnchor, newDir, m4.transformPoint(pMin), m4.transformPoint(pMax));
    }
    hashCode() {
        return this.anchor.hashCode() * 31 + this.dir1.hashCode();
    }
}
L3.X = new L3(V3.O, V3.X);
L3.Y = new L3(V3.O, V3.Y);
L3.Z = new L3(V3.O, V3.Z);
L3.prototype.hlol = Curve.hlol++;
L3.prototype.tIncrement = 256;

class PICurve extends ImplicitCurve {
    constructor(points, tangents, parametricSurface, implicitSurface, pmPoints, pmTangents, stepSize, dir = 1, generator, tMin, tMax) {
        super(points, tangents, dir, generator, tMin, tMax);
        this.parametricSurface = parametricSurface;
        this.implicitSurface = implicitSurface;
        this.pmPoints = pmPoints;
        this.pmTangents = pmTangents;
        this.stepSize = stepSize;
        assert(Array.isArray(pmPoints));
        assert(dir == 1);
        assert(stepSize <= 1);
        const pf = parametricSurface.pUVFunc();
        const dpdu = parametricSurface.dpdu();
        const dpdv = parametricSurface.dpdv();
        const didp = implicitSurface.didp.bind(implicitSurface);
        this.didu = (u, v) => didp(pf(u, v)).dot(dpdu(u, v));
        this.didv = (u, v) => didp(pf(u, v)).dot(dpdv(u, v));
        for (let i = 0; i < points.length - 1; i++) {
            assert(!points[i].equals(points[i + 1]));
            //assert(parametricSurface.pUV(pmPoints[i].x, pmPoints[i].y).equals(points[i]))
        }
        {
            const ps = this.parametricSurface;
            const is = implicitSurface;
            const pFunc = ps.pUVFunc(), iFunc = is.implicitFunction();
            const dpdu = ps.dpdu();
            const dpdv = ps.dpdv();
            const didp = is.didp.bind(is);
            const mf = MathFunctionR2R.forFFxFy((x, y) => iFunc(pFunc(x, y)), (u, v) => didp(pFunc(u, v)).dot(dpdu(u, v)), (u, v) => didp(pFunc(u, v)).dot(dpdv(u, v)));
            const { points } = followAlgorithm2d(mf, this.pmPoints[0], stepSize, ps, (u, v) => is.containsPoint(pFunc(u, v)), this.pmPoints.last, this.pmTangents[0]);
            if (points.length !== this.points.length) {
                followAlgorithm2d(mf, this.pmPoints[0], stepSize, ps, (u, v) => is.containsPoint(pFunc(u, v)), this.pmPoints.last, this.pmTangents[0]);
            }
            assert(points.length == this.points.length, points.length, this.points.length);
        }
    }
    static forParametricStartEnd(ps, is, pmStart, pmEnd, stepSize = 0.02, startPMTangent, tMin, tMax) {
        const pFunc = ps.pUVFunc(), iFunc = is.implicitFunction();
        const dpdu = ps.dpdu();
        const dpdv = ps.dpdv();
        const didp = is.didp.bind(is);
        const mf = MathFunctionR2R.forFFxFy((x, y) => iFunc(pFunc(x, y)), (u, v) => didp(pFunc(u, v)).dot(dpdu(u, v)), (u, v) => didp(pFunc(u, v)).dot(dpdv(u, v)));
        const { points, tangents } = followAlgorithm2d(mf, pmStart, stepSize, ps, (u, v) => is.containsPoint(pFunc(u, v)), pmEnd, startPMTangent);
        return PICurve.forParametricPointsTangents(ps, is, points, tangents, stepSize, 1, tMin, tMax);
    }
    static forStartEnd(ps, is, start, end, stepSize = 0.02, startTangent, min, max) {
        const startPM = ps.uvP(start);
        const dpdu = ps.dpdu()(startPM.x, startPM.y), dpdv = ps.dpdv()(startPM.x, startPM.y);
        const startPMTangent = startTangent &&
            M4.forSys(dpdu, dpdv).inversed().transformVector(startTangent);
        // assert(dpdu.times(startPMTangent.x).plus(dpdv.times(startPMTangent.y)).like(startTangent))
        const curve = PICurve.forParametricStartEnd(ps, is, startPM, ps.uvP(end), stepSize, startPMTangent);
        return curve.withBounds(min && curve.pointT(min), max && curve.pointT(max));
    }
    static forParametricPointsTangents(ps, is, pmPoints, pmTangents, stepSize, dir = 1, tMin, tMax) {
        const pFunc = ps.pUVFunc(), dpdu = ps.dpdu();
        const dpdv = ps.dpdv();
        const points = pmPoints.map(({ x, y }) => pFunc(x, y));
        const tangents = pmPoints.map(({ x: u, y: v }, i) => {
            const ds = dpdu(u, v);
            const dt = dpdv(u, v);
            return ds.times(pmTangents[i].x).plus(dt.times(pmTangents[i].y));
            //const p = points[i]
            //return cs.normalP(p).cross(ses.normalP(p))
            //	.toLength(ds.times(pmTangents[i].x).plus(dt.times(pmTangents[i].y)).length())
        });
        return new PICurve(points, tangents, ps, is, pmPoints, pmTangents, stepSize, dir, undefined, tMin, tMax);
    }
    getConstructorParameters() {
        return [
            this.points,
            this.tangents,
            this.parametricSurface,
            this.implicitSurface,
            this.pmPoints,
            this.pmTangents,
            this.stepSize,
            this.dir,
            this.generator,
        ];
    }
    implicitCurve() {
        const pF = this.parametricSurface.pUVFunc();
        const iF = this.implicitSurface.implicitFunction();
        return (u, v) => iF(pF(u, v));
    }
    isColinearTo(curve) {
        if (curve instanceof PICurve) {
            if (this.equals(curve)) {
                return true;
            }
            if (this.parametricSurface.isCoplanarTo(curve.parametricSurface) &&
                this.implicitSurface.isCoplanarTo(curve.implicitSurface)) ;
            return false;
        }
        else {
            return false;
        }
    }
    containsPoint(p) {
        assertVectors(p);
        const t = this.pointT(p);
        return !isNaN(t) && this.isValidT(t);
    }
    equals(obj) {
        return (Object.getPrototypeOf(obj) == PICurve.prototype &&
            this.parametricSurface.equals(obj.parametricSurface) &&
            this.implicitSurface.equals(obj.implicitSurface) &&
            this.points[0].equals(obj.points[0]) &&
            this.tangents[0].equals(obj.tangents[0]) &&
            this.dir === obj.dir);
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
        assertVectors(point);
        assert(this.containsPoint(point), "this.containsPoint(point)");
        const t = this.pointT(point);
        return this.tangentAt(t);
    }
    tangentAt(t) {
        assert(!isNaN(t));
        if (0 === t % 1)
            return this.tangents[t];
        const uv = this.uvT(t);
        const uvTangent = new V3(-this.didv(uv.x, uv.y), this.didu(uv.x, uv.y), 0).toLength(this.stepSize);
        const du = this.parametricSurface.dpdu()(uv.x, uv.y);
        const dv = this.parametricSurface.dpdv()(uv.x, uv.y);
        return du.times(uvTangent.x).plus(dv.times(uvTangent.y));
    }
    at(t) {
        assert(!isNaN(t));
        if (0 === t % 1)
            return this.points[t];
        const startParams = V3.lerp(this.pmPoints[floor(t)], this.pmPoints[ceil(t)], t % 1);
        return this.closestPointToParams(startParams);
    }
    uvT(t) {
        assert(!isNaN(t));
        if (0 === t % 1)
            return this.pmPoints[t];
        const startParams = V3.lerp(this.pmPoints[floor(t)], this.pmPoints[ceil(t)], t % 1);
        return curvePoint(this.implicitCurve(), startParams, this.didu, this.didv);
    }
    closestTToPoint(p, tStart) {
        // TODO
        return 0;
    }
    closestPointToParams(startUV) {
        const pointParams = curvePoint(this.implicitCurve(), startUV, this.didu, this.didv);
        return this.parametricSurface.pUVFunc()(pointParams.x, pointParams.y);
    }
    isTsWithSurface(surface) {
        if (surface instanceof EllipsoidSurface) {
            const pS = this.parametricSurface, iS = this.implicitSurface;
            if (pS instanceof ProjectedCurveSurface &&
                iS instanceof EllipsoidSurface) {
                const iscs = iS.isCurvesWithSurface(surface);
                const points = iscs.flatMap((isc) => isc.isTsWithSurface(pS).map((t) => isc.at(t)));
                const ts = fuzzyUniques(points.map((p) => this.pointT(p)));
                return ts.filter((t) => !isNaN(t) && this.isValidT(t));
            }
        }
        else if (ImplicitSurface.is(surface)) {
            const result = [];
            const iF = surface.implicitFunction();
            let prevSignedDistance = iF(this.points[0]);
            for (let i = 1; i < this.points.length; i++) {
                const point = this.points[i];
                const signedDistance = iF(point);
                if (prevSignedDistance * signedDistance <= 0) {
                    const pF = this.parametricSurface.pUVFunc();
                    const dpdu = this.parametricSurface.dpdu();
                    const dpdv = this.parametricSurface.dpdv();
                    const startUV = this.pmPoints[abs(prevSignedDistance) < abs(signedDistance) ? i - 1 : i];
                    const isUV = newtonIterate2dWithDerivatives(this.implicitCurve(), (u, v) => iF(pF(u, v)), startUV.x, startUV.y, 4, this.didu, this.didv, (u, v) => dpdu(u, v).dot(surface.didp(pF(u, v))), (u, v) => dpdv(u, v).dot(surface.didp(pF(u, v))));
                    result.push(this.pointT(this.parametricSurface.pUV(isUV.x, isUV.y)));
                }
                prevSignedDistance = signedDistance;
            }
            return result;
        }
        throw new Error();
    }
    isTsWithPlane(planeWC) {
        return this.isTsWithSurface(new PlaneSurface(planeWC));
        // version which intersects the plane with the defining surfaces of this PICurve, but this causes
        // issues when they are PICurves too:
        // assertInst(P3, planeWC)
        // const ps = this.parametricSurface,
        // 	is = this.implicitSurface
        // const pscs = ps.isCurvesWithPlane(planeWC)
        // const iscs = is.isCurvesWithPlane(planeWC)
        // const infos = iscs.flatMap(isc => pscs.flatMap(psc => isc.isInfosWithCurve(psc)))
        // const ts = fuzzyUniques(infos.map(info => this.pointT(info.p)))
        // return ts.filter(t => !isNaN(t) && this.isValidT(t))
    }
    pointT(p) {
        assertVectors(p);
        if (!this.parametricSurface.containsPoint(p) ||
            !this.implicitSurface.containsPoint(p)) {
            return NaN;
        }
        const pmPoint = this.parametricSurface.uvPFunc()(p);
        const ps = this.points, pmps = this.pmPoints;
        let t = 0, pmDistance = pmPoint.distanceTo(pmps[0]);
        while (pmDistance > abs(this.stepSize) && t < ps.length - 1) {
            // TODO -1?
            //console.log(t, pmps[t].$, pmDistance)
            t = min(pmps.length - 1, t + max(1, Math.round(pmDistance / abs(this.stepSize) / 2 / 2)));
            pmDistance = pmPoint.distanceTo(pmps[t]);
        }
        // if (t < this.pmPoints.length - 1 && pmDistance > pmPoint.distanceTo(pmps[t + 1])) {
        //     t++
        // }
        if (pmDistance > abs(this.stepSize) * 1.1) {
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
        const startT = arrayRange(floor(this.tMin), ceil(this.tMax), 1).withMax((t) => -pmPoint.distanceTo(pmps[t]));
        if (undefined === startT)
            throw new Error();
        if (ps[startT].like(p))
            return startT;
        //const [a, b] = 0 === startT
        //    ? [0, 1]
        //    : this.points.length - 1 === startT
        //        ? [startT - 1, startT]
        //        : pmPoint.distanceTo(pmps[startT - 1]) < pmPoint.distanceTo(pmps[startT + 1])
        //            ? [startT - 1, startT]
        //            : [startT, startT + 1]
        const a = max(0, startT - 1), b = min(this.points.length - 1, startT + 1);
        const tangent = this.tangentAt(startT);
        const f = (t) => this.at(clamp(t, 0, this.points.length - 1))
            .to(p)
            .dot(tangent);
        // const df = (t: number) => -this.tangentAt(clamp(t, 0, this.points.length - 1)).dot(tangent)
        //checkDerivate(f, df, 0, this.points.length - 2, 3)
        // 8 steps necessary because df can currently be way off
        t = bisect(f, a, b, 32);
        if (!isFinite(t) || this.at(t).distanceTo(p) > abs(this.stepSize)) {
            return NaN;
        }
        return t;
    }
    transform(m4) {
        const dirFactor = m4.isMirroring() ? -1 : 1;
        return PICurve.forStartEnd(this.parametricSurface.transform(m4), this.implicitSurface.transform(m4), m4.transformPoint(this.points[0]), m4.transformPoint(this.points.last), this.stepSize * dirFactor, m4.transformVector(this.tangents[0]), m4.transformPoint(this.at(this.tMin)), m4.transformPoint(this.at(this.tMax)));
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
        const allTs = arrayRange(0, this.points.length);
        return [allTs, allTs, allTs];
    }
    isInfosWithLine(anchorWC, dirWC, tMin, tMax, lineMin, lineMax) {
        return surfaceIsICurveIsInfosWithLine.call(this, anchorWC, dirWC, tMin, tMax, lineMin, lineMax);
    }
    toSource(rounder = (x) => x) {
        const result = callsce("PICurve.forParametricStartEnd", this.parametricSurface, this.implicitSurface, this.pmPoints[0], this.pmPoints.last, this.stepSize, this.pmTangents[0], this.tMin, this.tMax);
        return result;
    }
}
PICurve.prototype.tIncrement = 1;

class PPCurve extends ImplicitCurve {
    constructor(points, tangents, parametricSurface1, parametricSurface2, st1s, pmTangents, stepSize, dir = 1, generator, tMin, tMax) {
        super(points, tangents, dir, generator, tMin, tMax);
        this.parametricSurface1 = parametricSurface1;
        this.parametricSurface2 = parametricSurface2;
        this.st1s = st1s;
        this.pmTangents = pmTangents;
        this.stepSize = stepSize;
        assert(ParametricSurface.is(parametricSurface1));
        assert(ParametricSurface.is(parametricSurface2));
        assert(Array.isArray(st1s));
        assert(dir == 1);
        assert(stepSize <= 1);
    }
    at(t) {
        assert(!isNaN(t));
        if (0 === t % 1)
            return this.points[t];
        const startPoint = V3.lerp(this.points[floor(t)], this.points[ceil(t)], t % 1);
        return curvePointPP(this.parametricSurface1, this.parametricSurface2, startPoint).p;
    }
    isColinearTo(curve) {
        if (curve instanceof PPCurve) {
            if (this.equals(curve)) {
                return true;
            }
            if (this.parametricSurface1.isCoplanarTo(curve.parametricSurface1) &&
                this.parametricSurface1.isCoplanarTo(curve.parametricSurface2)) ;
            return false;
        }
        else {
            return false;
        }
    }
    containsPoint(p) {
        assertVectors(p);
        // TODO: wrong, as there could be another curve
        return (this.parametricSurface1.containsPoint(p) &&
            this.parametricSurface2.containsPoint(p) &&
            !isNaN(this.pointT(p)));
    }
    rootPoints() {
        const pF1 = this.parametricSurface1.pUVFunc();
        const pF2 = this.parametricSurface2.pUVFunc();
        const pN1 = this.parametricSurface1.normalUVFunc();
        const pN2 = this.parametricSurface2.normalUVFunc();
        const rootsApprox = this.rootsApprox();
        const results = [[], [], []];
        for (let dim = 0; dim < 3; dim++) {
            for (let i = 0; i < rootsApprox[dim].length; i++) {
                const lambda = rootsApprox[dim][i];
                const p = this.at(lambda);
                assert(this.parametricSurface1.containsPoint(p));
                const pp1 = this.parametricSurface1.uvP(p);
                const { x: u, y: v } = this.parametricSurface2.uvP(p);
                const startValues = [pp1.x, pp1.y, u, v];
                function f(vals) {
                    const [u1, v1, u2, v2] = vals;
                    const diff = pF1(u1, v1).minus(pF2(u2, v2));
                    const n1 = pN1(u1, v1);
                    const n2 = pN2(u2, v2);
                    const tangent = n1.cross(n2);
                    return [diff.x, diff.y, diff.z, tangent.e(dim)];
                }
                const pps = newtonIterate(f, startValues, 8);
                // assert(pF1(pps[0], pps[1]).like(pF2(pps[2], pps[3])),
                // 	pF1(pps[0], pps[1]).sce + pF2(pps[2], pps[3]).sce)
                const result = pF1(pps[0], pps[1]);
                results[dim].push(result);
            }
        }
        return results;
    }
    roots() {
        return this.rootPoints().map((ps) => ps.map((p) => this.pointT(p)));
    }
    pointTangent(pWC) {
        assertVectors(pWC);
        assert(this.containsPoint(pWC), "this.containsPoint(pWC)");
        const n1 = this.parametricSurface1.normalP(pWC);
        const n2 = this.parametricSurface2.normalP(pWC);
        return n1.cross(n2);
    }
    transform(m4) {
        return new PPCurve(m4.transformedPoints(this.points), m4.transformedVectors(this.tangents), this.parametricSurface1.transform(m4), this.parametricSurface2.transform(m4), this.st1s, undefined, this.stepSize, this.dir, undefined);
    }
    toSource() {
        return callsce("PPCurve.forStartEnd", this.parametricSurface1, this.parametricSurface2, this.points[0], this.points.last, this.stepSize);
    }
    static forStartEnd(ps1, ps2, startPoint, end, stepSize = 0.02) {
        const { points, tangents, st1s } = followAlgorithmPP(ps1, ps2, startPoint, stepSize);
        return new PPCurve(points, tangents, ps1, ps2, st1s, undefined, stepSize, 1);
    }
    isInfosWithLine(anchorWC, dirWC, tMin, tMax, lineMin, lineMax) {
        return surfaceIsICurveIsInfosWithLine.call(this, anchorWC, dirWC, tMin, tMax, lineMin, lineMax);
    }
    isTsWithSurface(surface) {
        if (ImplicitSurface.is(surface)) {
            const result = [];
            const iF = surface.implicitFunction();
            const pUV1 = this.parametricSurface1.pUVFunc();
            const pUV2 = this.parametricSurface2.pUVFunc();
            let prevSignedDistance = iF(this.points[0]);
            for (let i = 1; i < this.points.length; i++) {
                const point = this.points[i];
                const signedDistance = iF(point);
                if (prevSignedDistance * signedDistance <= 0) {
                    const startIndex = abs(prevSignedDistance) < abs(signedDistance) ? i - 1 : i;
                    const startPoint = this.points[startIndex];
                    const startUV1 = this.st1s[startIndex];
                    const startUV2 = this.parametricSurface2.uvP(startPoint);
                    const isSTUV = newtonIterate(([u1, v1, u2, v2]) => {
                        const ps1p = pUV1(u1, v1);
                        const ps2p = pUV2(u2, v2);
                        return [...ps1p.to(ps2p), iF(ps1p)];
                    }, [startUV1.x, startUV1.y, startUV2.x, startUV2.y]);
                    result.push(this.pointT(this.parametricSurface1.pUV(isSTUV[0], isSTUV[1])));
                }
                prevSignedDistance = signedDistance;
            }
            return result;
        }
        throw new Error("Method not implemented.");
    }
    isTsWithPlane(planeWC) {
        return this.isTsWithSurface(new PlaneSurface(planeWC));
    }
}

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
        const lineTs = pqFormula((anchorLC.x * dirLC.x + dirLC.y) / pqDiv, (Math.pow(anchorLC.x, 2) + anchorLC.y) / pqDiv);
        return lineTs
            .filter((tOther) => le(0, anchorLC.y + tOther * dirLC.y))
            .map((tOther) => ({
            tThis: dirLC.x * tOther + anchorLC.x,
            tOther: tOther,
            p: L3.at(anchorWC, dirWC, tOther),
        }));
    }
    static intersectionUnitLine(a, b, c) {
        /*
             solve system (5)/(6)
             g1 * xi + g2 * eta = g3 (6)
             g1 * xi + g2 * xi * xi = g3
             xi² + xi * g1/g2 - g3/g2 = 0
             */
        return pqFormula(a / b, -c / b);
    }
    static XYLCValid(pLC) {
        return eq(Math.pow(pLC.x, 2), pLC.y);
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
        assertNumbers(t);
        // f1 + f2 2 t
        return this.f1.plus(this.f2.times(2 * t));
    }
    ddt(t) {
        assertNumbers(t);
        return this.f2.times(2);
    }
    tangentAt2(xi, eta) {
        assertNumbers(xi, eta);
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
        const dimRoots = (dim) => eq0(this.f2.e(dim)) ? [] : [-this.f1.e(dim) / 2 / this.f2.e(dim)];
        return arrayFromFunction(3, dimRoots);
    }
    isColinearTo(curve) {
        if (!hasConstructor(curve, ParabolaCurve))
            return false;
        const thisRA = this.rightAngled(), curveRA = curve.rightAngled();
        return (thisRA.center.like(curveRA.center) &&
            thisRA.f2.like(curveRA.f2) &&
            thisRA.f1.likeOrReversed(curveRA.f1));
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
        if (eq0(f1DOTf2) && f1.hasLength(1)) {
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
        if (!eq0(f1DOTf2)) {
            t0 = -f1DOTf2 / f2.squared() / 2;
            f1 = f1.plus(f2.times(2 * t0));
        }
        const f1Length = f1.length();
        const a = f2.length() / f1Length;
        function F(x) {
            return (Math.asinh(a * 2 * x) / 4 / a +
                (x * Math.sqrt(1 + a * a * 4 * x * x)) / 2);
        }
        return f1Length * (F(endT - t0) - F(startT - t0));
    }
    transform4(m4) {
        return parabola4Projection(this.matrix.transform(m4), this.tMin, this.tMax);
    }
    asBezier() {
        return BezierCurve.quadratic(this.at(-1), new L3(this.at(-1), this.tangentAt(-1).unit()).isInfoWithLine(new L3(this.at(1), this.tangentAt(1).unit())), this.at(1));
    }
    /**
     * Returns new ParabolaCurve that has its center point at this.at(t0)
     * @param t0
     */
    recenter(t0) {
        // this.at(t) = f2 t² + f1 t + center
        // c2.at(t) = f2 (t + t0)² + f1 (t + t0) + center
        // c2.at(t) = f2 (t² + 2 t0 t + t0²) + f1 (t + t0) + center
        // c2.at(t) = f2 t² + (f1 + 2 f2 t0) t + center + f2 t0² + f1 t0
        return new ParabolaCurve(this.at(t0), this.f1.plus(this.f2.times(2 * t0)), this.f2);
    }
}
ParabolaCurve.XY = new ParabolaCurve(V3.O, V3.X, V3.Y);
ParabolaCurve.YZ = new ParabolaCurve(V3.O, V3.Y, V3.Z);
ParabolaCurve.ZX = new ParabolaCurve(V3.O, V3.Z, V3.X);
ParabolaCurve.prototype.tIncrement = 1 / 32;

class EllipseCurve extends XiEtaCurve {
    constructor(center, f1, f2, tMin = 0, tMax = PI) {
        super(center, f1, f2, tMin, tMax);
        assert(-PI <= this.tMin && this.tMin < PI);
        assert(-PI < this.tMax && this.tMax <= PI);
    }
    static andFixTs(center, f1, f2, tMin = 0, tMax = PI) {
        if (-PI <= tMin && tMax <= PI) {
            return new EllipseCurve(center, f1, f2, tMin, tMax);
        }
        if (0 <= tMin && tMax <= TAU) {
            return new EllipseCurve(center, f1.negated(), f2.negated(), tMin - PI, tMax - PI);
        }
        if (-TAU <= tMin && tMax <= 0) {
            return new EllipseCurve(center, f1.negated(), f2.negated(), tMin + PI, tMax + PI);
        }
        throw new Error("Method not implemented.");
    }
    static XYLCValid(pLC) {
        const { x, y } = pLC;
        return eq0(Math.pow(x, 2) + Math.pow(y, 2) - 1);
    }
    static XYLCPointT(pLC, tMin, tMax) {
        assertNumbers(tMin, tMax);
        const t = atan2(pLC.y, pLC.x);
        const lowSplitter = lerp(tMin, tMax - TAU, 0.5);
        if (t < lowSplitter) {
            return t + TAU;
        }
        const highSplitter = lerp(tMax, tMin + TAU, 0.5);
        if (t > highSplitter) {
            return t - TAU;
        }
        return t;
    }
    static intersectionUnitLine(a, b, c, tMin, tMax) {
        const isLC = intersectionUnitCircleLine2(a, b, c);
        const result = [];
        for (const [xi, eta] of isLC) {
            const t = EllipseCurve.XYLCPointT(new V3(xi, eta, 0), tMin, tMax);
            fuzzyBetween(t, tMin, tMax) && result.push(t);
        }
        return result;
    }
    static unitIsInfosWithLine(anchorLC, dirLC, anchorWC, dirWC, tMin, tMax) {
        // ell: x² + y² = 1 = p²
        // line(t) = anchor + t dir
        // anchor² - 1 + 2 t dir anchor + t² dir² = 0
        const pqDiv = dirLC.squared();
        const lineTs = pqFormula((2 * dirLC.dot(anchorLC)) / pqDiv, (anchorLC.squared() - 1) / pqDiv);
        return lineTs
            .filter((tOther) => le(0, anchorLC.y + tOther * dirLC.y))
            .map((tOther) => ({
            tThis: EllipseCurve.XYLCPointT(dirLC.times(tOther).plus(anchorLC), tMin, tMax),
            tOther: tOther,
            p: L3.at(anchorWC, dirWC, tOther),
        }));
    }
    /**
     * Returns a new EllipseCurve representing a circle parallel to the XY-plane.`
     */
    static semicircle(radius, center = V3.O, tMin, tMax) {
        return new EllipseCurve(center, new V3(radius, 0, 0), new V3(0, radius, 0), tMin, tMax);
    }
    static circleForCenter2P(center, a, b, radius, tMin, tMax) {
        const f1 = center.to(a);
        const normal = f1.cross(center.to(b));
        const f2 = normal.cross(f1).toLength(f1.length());
        return new EllipseCurve(center, f1, f2, undefined !== tMin ? tMin : 0, undefined !== tMax ? tMax : f1.angleTo(center.to(b)));
    }
    split(tMin = this.tMin, tMax = this.tMax) {
        const result = [];
        tMin < 0 &&
            result.push(new EllipseCurve(this.center, this.f1.negated(), this.f2.negated(), tMin + PI, min(0, tMax) + PI));
        tMax > 0 &&
            result.push(new EllipseCurve(this.center, this.f1, this.f2, max(0, tMin), tMax));
        return result;
    }
    static forAB(a, b, center = V3.O) {
        return super.forAB(a, b, center);
    }
    /**
     * Create a circle curve which has a, b and c on it. a, b, c can't be on a straight line.
     * tMin defaults to 0, tMax defaults to the value for c
     */
    static circleThroughPoints(a, b, c, tMin = 0, tMax) {
        assertf(() => !L3.throughPoints(a, c).containsPoint(b));
        const normal = a.to(b).cross(b.to(c));
        const center = new L3(a.lerp(b, 0.5), normal.cross(a.to(b)).unit()).isInfoWithLine(new L3(b.lerp(c, 0.5), normal.cross(b.to(c)).unit()));
        const f1 = center.to(a).negated();
        return new EllipseCurve(center, f1, normal.unit().cross(f1), -PI, undefined === tMax
            ? f1.angleRelativeNormal(center.to(c), normal.unit())
            : tMax);
    }
    getAreaInDir(right, up, tStart, tEnd) {
        //assertf(() => tStart < tEnd)
        assertf(() => right.isPerpendicularTo(this.normal));
        assertf(() => up.isPerpendicularTo(this.normal));
        //assertf(() => EllipseCurve.isValidT(tStart), tStart)
        //assertf(() => EllipseCurve.isValidT(tEnd), tEnd)
        const upLC = this.matrixInverse.transformVector(up);
        const rightLC = upLC.cross(V3.Z);
        const normTStart = tStart - rightLC.angleXY();
        const normTEnd = tEnd - rightLC.angleXY();
        const transformedOriginY = this.matrixInverse
            .getTranslation()
            .dot(upLC.unit());
        // integral of sqrt(1 - x²) from 0 to cos(t)
        // Basically, we want
        // INTEGRAL[cos(t); PI/2] sqrt(1 - x²) dx
        // INTEGRAL[PI/2: cos(t)] -sqrt(1 - x²) dx
        // = INTEGRAL[cos(0); cos(t)] -sqrt(1 - x²) dx
        // = INTEGRAL[0; t] -sqrt(1 - cos²(t)) * -sin(t) dt
        // = INTEGRAL[0; t] -sin(t) * -sin(t) dt
        // = INTEGRAL[0; t] sin²(t) dt (partial integration / wolfram alpha)
        // = (1/2 * (t - sin(t) * cos(t)))[0; t] (this form has the distinct advantage of being defined everywhere)
        function fArea(t) {
            return (t - Math.sin(t) * Math.cos(t)) / 2;
        }
        // for the centroid, we want
        // cx = 1 / area * INTEGRAL[cos(t); PI/2] x * f(x) dx
        // cx = 1 / area * INTEGRAL[cos(t); PI/2] x * sqrt(1 - x²) dx
        // cx = 1 / area * INTEGRAL[cos(0); cos(t)] x * -sqrt(1 - x²) dx
        // ...
        // cx = 1 / area * INTEGRAL[0; t] cos(t) * sin²(t) dt // WA
        // cx = 1 / area * (sin^3(t) / 3)[0; t]
        function cxTimesArea(t) {
            return Math.pow(Math.sin(t), 3) / 3;
        }
        // cy = 1 / area * INTEGRAL[cos(t); PI/2] f²(x) / 2 dx
        // cy = 1 / area * INTEGRAL[cos(0); cos(t)] -(1 - x²) / 2 dx
        // cy = 1 / area * INTEGRAL[0; t] (cos²(t) - 1) * -sin(t) / 2 dt
        // cy = 1 / area * (cos (3 * t) - 9 * cos(t)) / 24 )[0; t]
        function cyTimesArea(t) {
            return (Math.cos(3 * t) - 9 * Math.cos(t)) / 24;
        }
        const restArea = -transformedOriginY * (-Math.cos(normTEnd) + Math.cos(normTStart));
        const area = fArea(normTEnd) - fArea(normTStart) + restArea;
        const cxt = (cxTimesArea(normTEnd) -
            cxTimesArea(normTStart) +
            ((-transformedOriginY * (-Math.cos(normTEnd) - Math.cos(normTStart))) /
                2) *
                restArea) /
            area;
        const cyt = (cyTimesArea(normTEnd) -
            cyTimesArea(normTStart) -
            (-transformedOriginY / 2) * restArea) /
            area;
        const factor = this.matrix.xyAreaFactor(); // * upLC.length()
        //console.log('fctor', factor, 'area', area, 'resultarea', area* factor)
        assert(!eq0(factor));
        return {
            area: area * factor,
            centroid: this.matrix.transformPoint(M4.rotateZ(rightLC.angleXY()).transformPoint(new V3(cxt, cyt, 0))),
        };
    }
    at(t) {
        assertNumbers(t);
        //assert(this.isValidT(t))
        // = center + f1 cos t + f2 sin t
        return this.center
            .plus(this.f1.times(Math.cos(t)))
            .plus(this.f2.times(Math.sin(t)));
    }
    tangentAt(t) {
        assertNumbers(t);
        //assert(this.isValidT(t))
        // ) f2 cos(t) - f1 sin(t)
        return this.f2.times(Math.cos(t)).minus(this.f1.times(Math.sin(t)));
    }
    ddt(t) {
        assertNumbers(t);
        assert(this.isValidT(t));
        return this.f2.times(-Math.sin(t)).minus(this.f1.times(Math.cos(t)));
    }
    tangentAt2(xi, eta) {
        return this.f2.times(xi).minus(this.f1.times(eta));
    }
    isCircular() {
        return (eq(this.f1.length(), this.f2.length()) &&
            this.f1.isPerpendicularTo(this.f2));
    }
    isColinearTo(curve) {
        if (!hasConstructor(curve, EllipseCurve))
            return false;
        if (!this.center.like(curve.center)) {
            return false;
        }
        if (this == curve) {
            return true;
        }
        if (this.isCircular()) {
            return (curve.isCircular() &&
                eq(this.f1.length(), curve.f1.length()) &&
                this.normal.isParallelTo(curve.normal));
        }
        else {
            let { f1: f1, f2: f2 } = this.rightAngled(), { f1: c1, f2: c2 } = curve.rightAngled();
            if (f1.length() > f2.length()) {
                [f1, f2] = [f2, f1];
            }
            if (c1.length() > c2.length()) {
                [c1, c2] = [c2, c1];
            }
            return (eq(f1.squared(), Math.abs(f1.dot(c1))) &&
                eq(f2.squared(), Math.abs(f2.dot(c2))));
        }
    }
    pointT(pWC) {
        assertVectors(pWC);
        assert(this.containsPoint(pWC));
        const pLC = this.matrixInverse.transformPoint(pWC);
        const t = EllipseCurve.XYLCPointT(pLC, this.tMin, this.tMax);
        assert(this.isValidT(t));
        return t;
    }
    reversed() {
        return new EllipseCurve(this.center, this.f1.negated(), this.f2, PI - this.tMax, PI - this.tMin);
    }
    eccentricity() {
        const mainAxes = this.rightAngled();
        const f1length = mainAxes.f1.length(), f2length = mainAxes.f1.length();
        const [a, b] = f1length > f2length ? [f1length, f2length] : [f2length, f1length];
        return Math.sqrt(1 - (b * b) / a / a);
    }
    circumference() {
        return this.arcLength(-Math.PI, Math.PI);
    }
    arcLength(tStart = this.tMin, tEnd = this.tMax, steps = 2) {
        assert(tStart < tEnd, "startT < endT");
        const f1Length = this.f1.length();
        if (eq(f1Length, this.f2.length())) {
            return f1Length * (tEnd - tStart);
        }
        return super.arcLength(tStart, tEnd, steps);
    }
    circumferenceApproximate() {
        // approximate circumference by Ramanujan
        // https://en.wikipedia.org/wiki/Ellipse#Circumference
        const { f1, f2 } = this.rightAngled(), a = f1.length(), b = f2.length();
        const h = Math.pow((a - b), 2) / Math.pow((a + b), 2);
        return Math.PI * (a + b) * (1 + (3 * h) / (10 + Math.sqrt(4 - 3 * h)));
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
        if (eq0(a)) {
            return this;
        }
        const g1 = 2 * a, g2 = b + Math.sqrt(b * b + 4 * a * a);
        const { x1: xi, y1: eta } = intersectionUnitCircleLine(g1, g2, 0);
        const f1RA = f1.times(xi).plus(f2.times(eta));
        const f2RA = f1.times(-eta).plus(f2.times(xi));
        return new EllipseCurve(this.center, f1RA, f2RA, -PI, PI);
    }
    isInfosWithEllipse(ellipse) {
        if (this.normal.isParallelTo(ellipse.normal) &&
            eq0(this.center.minus(ellipse.center).dot(ellipse.normal))) {
            // ellipses are coplanar
            const ellipseLCRA = ellipse.transform(this.matrixInverse).rightAngled();
            const r1 = ellipseLCRA.f1.lengthXY(), r2 = ellipseLCRA.f2.lengthXY(), centerDist = ellipseLCRA.center.lengthXY();
            const rMin = min(r1, r2), rMax = max(r1, r2);
            if (lt(centerDist + rMax, 1) || // entirely inside unit circle
                lt(1, centerDist - rMax) || // entirely outside unit circle
                lt(1, rMin - centerDist) || // contains unit circle
                (eq(1, r1) && eq(1, r2) && eq0(centerDist)) // also unit circle, return no IS
            ) {
                return [];
            }
            const f = (t) => ellipseLCRA.at(t).lengthXY() - 1;
            const df = (t) => ellipseLCRA.at(t).xy().dot(ellipseLCRA.tangentAt(t)) /
                ellipseLCRA.at(t).lengthXY();
            checkDerivate(f, df, -PI, PI, 1);
            const ellipseLCRATs = [];
            for (let startT = (-4 / 5) * PI; startT < PI; startT += PI / 4) {
                let t = newtonIterateSmart(f, startT, 16, df, 1e-4);
                le(t, -PI) && (t += TAU);
                assert(!isNaN(t));
                if (between(t, -PI, PI) &&
                    eq0(f(t)) &&
                    !ellipseLCRATs.some((r) => eq(t, r))) {
                    ellipseLCRATs.push(t);
                }
            }
            const result = [];
            for (const ellipseLCRAT of ellipseLCRATs) {
                const p = this.matrix.transformPoint(ellipseLCRA.at(ellipseLCRAT));
                if (this.containsPoint(p) && ellipse.containsPoint(p)) {
                    result.push({ tThis: this.pointT(p), tOther: ellipse.pointT(p), p });
                }
            }
            return result;
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
            return this.isTsWithPlane(P3.normalOnAnchor(ellipse.normal.unit(), ellipse.center)).mapFilter((t) => {
                const p = this.at(t);
                if (ellipse.containsPoint(p)) {
                    return { tThis: t, tOther: ellipse.pointT(p), p };
                }
                return undefined;
            });
        }
    }
    isInfosWithCurve(curve) {
        if (curve instanceof EllipseCurve) {
            return this.isInfosWithEllipse(curve);
        }
        return super.isInfosWithCurve(curve);
    }
    transform4(m4) {
        const tMap = (t) => sign(t) * sqrt((1 - cos(t)) / (1 + cos(t)));
        // prettier-ignore
        const parabolaToUnitEllipse = new M4(0, -1, 0, 1, 2, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1);
        return parabola4Projection(M4.product(m4, this.matrix, parabolaToUnitEllipse), tMap(this.tMin), tMap(this.tMax));
    }
    roots() {
        // tangent(t) = f2 cos t - f1 sin t
        // solve for each dimension separately
        // tangent(eta, xi) = f2 eta - f1 xi
        return arrayFromFunction(3, (dim) => {
            const a = this.f2.e(dim), b = -this.f1.e(dim);
            return intersectionUnitCircleLine2(a, b, 0)
                .map(([xi, eta]) => Math.atan2(eta, xi))
                .filter((t) => this.isValidT(t));
        });
    }
    closestTToPoint(p, tStart) {
        // (at(t) - p) * tangentAt(t) = 0
        // (xi f1 + eta f2 + q) * (xi f2 - eta f1) = 0
        // xi eta (f2^2-f1^2) + xi f2 q - eta² f1 f2 + xi² f1 f2 - eta f1 q = 0
        //  (xi² - eta²) f1 f2 + xi eta (f2^2-f1^2) + xi f2 q - eta f1 q = 0
        // atan2 of p is a good first approximation for the searched t
        tStart = tStart || this.matrixInverse.transformPoint(p).angleXY();
        const pRelCenter = p.minus(this.center);
        const f = (t) => this.tangentAt(t).dot(this.f1
            .times(Math.cos(t))
            .plus(this.f2.times(Math.sin(t)))
            .minus(pRelCenter));
        return newtonIterate1d(f, tStart, 8);
    }
    area() {
        // see
        // https://upload.wikimedia.org/wikipedia/commons/thumb/4/4e/Cross_product_parallelogram.svg/220px-Cross_product_parallelogram.svg.png
        return Math.PI * this.f1.cross(this.f2).length();
    }
    angleToT(phi) {
        // atan2(y, x) = phi
        const phiDir = this.f1
            .unit()
            .times(Math.cos(phi))
            .plus(this.f2.rejectedFrom(this.f1).unit().times(Math.sin(phi)));
        const dirLC = this.matrixInverse.transformVector(phiDir);
        return dirLC.angleXY();
    }
}
EllipseCurve.UNIT = new EllipseCurve(V3.O, V3.X, V3.Y);
EllipseCurve.prototype.hlol = Curve.hlol++;
EllipseCurve.prototype.tIncrement = (2 * Math.PI) / (4 * 32);

/**
 * Non-Uniform Rational B-Spline implementation.
 *
 * See https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/ for a good reference.
 *
 *
 */
class NURBS extends Curve {
    constructor(
    /**
     * The control points of the NURBS curve, as 4D homogeneous coordinates.
     */
    points, 
    /**
     * The degree of the NURBS curve. Must be at least 1 (linear).
     */
    degree, knots = NURBS.openUniformKnots(points.length, degree), tMin = knots[degree], tMax = knots[knots.length - degree - 1]) {
        super(tMin, tMax);
        this.points = points;
        this.degree = degree;
        this.knots = knots;
        const knotsLength = points.length + degree + 1;
        NLA_DEBUG && Object.freeze(points);
        NLA_DEBUG && Object.freeze(knots);
        assert(knots.length === knotsLength, "bad knot vector length: expected " +
            knotsLength +
            " (degree = " +
            degree +
            " pcount = " +
            points.length +
            "), but was " +
            knots.length);
        assert(knots[degree] <= tMin);
        assert(tMax <= knots[knots.length - degree - 1]);
        for (let i = 0; i < points.length; i++) {
            assert(points[i].dim() == 4);
        }
        assert(degree >= 1, "degree must be at least 1 (linear)");
        assert(degree % 1 == 0);
        assert(-1 == firstUnsorted(knots, MINUS), "knot values must be in ascending order");
    }
    getConstructorParameters() {
        return [this.points, this.degree, this.knots];
    }
    at4(t) {
        NLA_DEBUG && assert(between(t, this.tMin, this.tMax), t);
        const { points, degree, knots } = this;
        // find s (the spline segment) for the [t] value provided
        const s = this.tInterval(t);
        const v = Vector.pack(points, new Float64Array((degree + 1) * 4), s - degree, 0, degree + 1);
        for (let level = 0; level < degree; level++) {
            // build level l of the pyramid
            for (let i = degree; i > level; i--) {
                const alpha = (t - knots[i + s - degree]) /
                    (knots[i + s - level] - knots[i + s - degree]);
                // interpolate each component
                for (let dim = 0; dim < 4; dim++) {
                    v[i * 4 + dim] =
                        (1 - alpha) * v[(i - 1) * 4 + dim] + alpha * v[i * 4 + dim];
                }
            }
        }
        return new Vector(v.slice(degree * 4, (degree + 1) * 4));
    }
    at(t) {
        return this.at4(t).p3();
    }
    /*
      d(k, i, t) = a(i, k, t) * d(k - 1, i, t) + (1 - a(i, k, t)) * d(k - 1, i - 1, t)
      a(i, k, t) = (t - knots[i]) / (knots[i + 1 + n - k] - knots[i])
      a'(i, k, t) = 1 / (knots[i + 1 + n - k] - knots[i])
  
      d/dt =  a(i, k, t) * d'(k - 1, i, t) + a'(i, k, t) * d(k - 1, i, t)
      + (1 - a(i, k, t)) * d'(k - 1, i - 1, t) + a'(i, k, t) * d(k - 1, i - 1, t)
  */
    ptDtDdt4(t) {
        const { points, degree, knots } = this;
        // find s (the spline segment) for the [t] value provided
        const s = this.tInterval(t);
        const v = Vector.pack(points, new Float64Array((degree + 1) * 4), s - degree, 0, degree + 1);
        let ddt = Vector.Zero(4), derivative;
        for (let level = 0; level < degree; level++) {
            if (level == degree - 2) {
                // see https://www.globalspec.com/reference/61012/203279/10-8-derivatives
                const a = new Vector(v.slice(degree * 4, (degree + 1) * 4));
                const b = new Vector(v.slice((degree - 1) * 4, degree * 4));
                const c = new Vector(v.slice((degree - 2) * 4, (degree - 1) * 4));
                function step(k, i, dkMinus1iMinus1, dkMinus1i) {
                    return dkMinus1i
                        .minus(dkMinus1iMinus1)
                        .times(k / (knots[i + degree - k] - knots[i - 1]));
                }
                ddt = step(degree, s + 1, step(degree - 1, s + 1, a, b), step(degree - 1, s, b, c));
            }
            if (level == degree - 1) {
                const a = new Vector(v.slice(degree * 4, (degree + 1) * 4));
                const b = new Vector(v.slice((degree - 1) * 4, degree * 4));
                derivative = b.minus(a).times(degree / (knots[s] - knots[s + 1]));
            }
            for (let i = degree; i > level; i--) {
                const alpha = (t - knots[i + s - degree]) /
                    (knots[i + s - level] - knots[i + s - degree]);
                // interpolate each component
                for (let dim = 0; dim < 4; dim++) {
                    v[i * 4 + dim] =
                        (1 - alpha) * v[(i - 1) * 4 + dim] + alpha * v[i * 4 + dim];
                }
            }
        }
        const p = new Vector(v.slice(degree * 4, degree * 4 + 4));
        return [p, derivative, ddt];
    }
    tangentAt(t) {
        // x(t) = xw(t) / w(t)
        // quotient rule
        const [p, derivative] = this.ptDtDdt4(t);
        const expected = derivative
            .times(p.w)
            .minus(p.times(derivative.w))
            .div(Math.pow(p.w, 2))
            .V3();
        return expected;
    }
    ddt(t) {
        const [p, dt, ddt] = this.ptDtDdt4(t);
        // =(-w(t) x(t) w''(t) - 2 w(t) w'(t) x'(t) + 2 x(t) w'(t)^2 + w(t)^2 x''(t))/w(t)^3
        // =(x(t) ((-w(t)) w''(t) + 2 w'(t)^2) - x'(t) 2 w(t) w'(t) + x''(t) w(t)^2 )/w(t)^3
        // prettier-ignore
        return Vector.add(p.times(-p.w * ddt.w + 2 * Math.pow(dt.w, 2)), dt.times(-2 * p.w * dt.w), ddt.times(Math.pow(p.w, 2))).div(Math.pow(p.w, 3)).V3();
    }
    ptDtDdt(t) {
        const [pt, dt4, ddt4] = this.ptDtDdt4(t);
        return [
            pt.p3(),
            dt4
                .times(pt.w)
                .minus(pt.times(dt4.w))
                .div(Math.pow(pt.w, 2))
                .V3(),
            Vector.add(pt.times(-pt.w * ddt4.w + 2 * Math.pow(dt4.w, 2)), //
            dt4.times(-2 * pt.w * dt4.w), ddt4.times(Math.pow(pt.w, 2)))
                .div(Math.pow(pt.w, 3))
                .V3(),
        ];
    }
    pointT(pWC) {
        return this.closestTToPoint(pWC);
    }
    closestTToPoint(p, tStart, tMin = this.tMin, tMax = this.tMax) {
        // this.at(t) has minimal distance to p when this.tangentAt(t) is perpendicular to
        // the vector between this.at(t) and p. This is the case iff the dot product of the two is 0.
        // f = (this.at(t) - p) . (this.tangentAt(t)
        // df = this.tangentAt(t) . this.tangentAt(t) + (this.at(t) - p) . this.ddt(t)
        //    = this.tangentAt(t)² + (this.at(t) - p) . this.ddt(t)
        const f = (t) => {
            const [pt, dt, ddt] = this.ptDtDdt(t);
            return [pt.minus(p).dot(dt), dt.squared() + pt.minus(p).dot(ddt)];
        };
        //checkDerivate(f, df, tMin, tMax)
        const STEPS = 32;
        if (undefined === tStart) {
            tStart = arraySamples(tMin, tMax, STEPS).withMax((t) => -this.at(t).distanceTo(p));
        }
        const result = newtonIterateWithDerivative2(f, tStart, 8, this.tMin, this.tMax);
        //assert(undefined !== result)
        return result;
    }
    containsPoint(pWC) {
        const tGuess = this.closestTToPoint(pWC);
        return undefined === tGuess ? false : this.at(tGuess).like(pWC);
    }
    derivate() {
        const k = this.degree;
        const ps = arrayFromFunction(this.points.length - 1, (i) => this.points[i]
            .to(this.points[i + 1])
            .times(k / (this.knots[i + k + 1] - this.knots[i + 1])));
        return new NURBS(ps, this.degree - 1, this.knots.slice(1, -1), this.tMin, this.tMax);
    }
    /**
     * Create a new NURBS of equal degree with the added knot [newKnot]. New NURBS will have one additional control
     * point.
     */
    withKnot(newKnot) {
        assert(between(newKnot, this.tMin, this.tMax));
        const k = this.tInterval(newKnot);
        const { knots, points, degree } = this;
        const insertPoints = arrayFromFunction(this.degree, (j) => {
            const i = k - degree + 1 + j;
            const aiNumerator = newKnot - knots[i];
            // 0/0 defined as 0:
            const ai = aiNumerator == 0 ? 0 : aiNumerator / (knots[i + degree] - knots[i]);
            assert(between(ai, 0, 1));
            return Vector.lerp(points[i - 1], points[i], ai);
        });
        const newPoints = points.slice();
        newPoints.splice(k - degree + 1, degree - 1, ...insertPoints);
        const newKnots = knots.slice();
        newKnots.splice(k + 1, 0, newKnot);
        return new NURBS(newPoints, degree, newKnots, this.tMin, this.tMax);
    }
    removeKnot(t) {
        const { knots, points, degree } = this;
        let k = this.tInterval(t), s = 0; // s = multiplicity of the knot
        while (knots[k + 1] == t) {
            k++;
            s++;
        }
        if (s == 0)
            throw new Error("There is no knot " + t + "!");
        // the points which were relevant when inserting were (k - p - 1) to (k - 1). (- 1) because the current k has
        // been increased by one due to the insertion.
        // p - 1 points were replaced by p points, hence we need to generate the original p - 1 point, + 1 to check if
        // this transformation is valid.
        const insertPoints = [points[k - degree - 1]];
        const oldKnots = knots.slice();
        oldKnots.splice(k, 1);
        for (let i = k - degree; i <= k - s; i++) {
            const alphaInv = (oldKnots[i + degree] - oldKnots[i]) / (t - oldKnots[i]);
            const oldPoint = Vector.lerp(insertPoints.last, points[i], alphaInv);
            insertPoints.push(oldPoint);
        }
        if (insertPoints.last.like(points[k + 1 - s])) {
            const oldPoints = points.slice();
            oldPoints.splice(k - degree - 1, degree - s + 3, ...insertPoints);
            return new NURBS(oldPoints, degree, oldKnots);
        }
        return undefined;
    }
    static openUniformKnots(pointCount, degree, tMin = 0, tMax = 1) {
        const knotsLength = pointCount + degree + 1;
        return arrayFromFunction(knotsLength, (i) => {
            if (i <= degree) {
                return tMin;
            }
            else if (i >= knotsLength - degree - 1) {
                return tMax;
            }
            else {
                return lerp(tMin, tMax, (i - degree) / (knotsLength - degree * 2 - 1));
            }
        });
    }
    static bezierKnots(degree, tMin = 0, tMax = 1) {
        const result = new Array((degree + 1) * 2);
        for (let i = 0; i < degree + 1; i++) {
            result[i] = tMin;
            result[degree + 1 + i] = tMax;
        }
        return result;
    }
    static fromBezier(bezier) {
        const bezier01 = bezier.selectPart(bezier.tMin, bezier.tMax);
        return NURBS.Bezier(bezier01.points);
    }
    static Bezier(points, tMin = 0, tMax = 1) {
        return new NURBS(points.map((p) => p instanceof V3 ? new Vector(new Float64Array([p.x, p.y, p.z, 1])) : p), points.length - 1, arrayFromFunction(points.length * 2, (i) => (i < points.length ? 0 : 1)), tMin, tMax);
    }
    static fromHyperbola(hyperbola, tMin = hyperbola.tMin, tMax = hyperbola.tMax) {
        const p0 = HyperbolaCurve.XY.at(tMin);
        const p2 = HyperbolaCurve.XY.at(tMax);
        const p1 = new V3((sinh(tMin) - sinh(tMax)) / sinh(tMin - tMax), (cosh(tMin) - cosh(tMax)) / sinh(tMin - tMax), 0);
        // M: midpoint between p0 and p2
        // X: intersection of line through p1 and M and unit hyperbola
        // result.at(1/2) = X
        // result.at(1/2) = (1/4 p0 + 1/2 p1 w + 1/4 p2) / (1/4 + 1/ 2 w + 1/4)
        // result.at(1/2) = (1/2 p0 + p1 w + 1/2 p2) / (1 + w)
        // result.at(1/2) = (M + p1 w) / (1 + w) = X
        // => w * (p1 - X) = (X - M)
        // as p1, X and M are all on the same line, we can solve this equation with only the x
        const M = p0.lerp(p2, 0.5);
        const Xx = 1 / sqrt(1 - Math.pow((M.y / M.x), 2));
        const w = (Xx - M.x) / (p1.x - Xx);
        return NURBS.fromV3s([p0, p1, p2], 2, undefined, [1, w, 1]).transform(hyperbola.matrix);
    }
    static fromParabola(parabola) {
        return NURBS.fromBezier(parabola.asBezier());
    }
    static fromEllipse(ellipse) {
        const unitSemiEllipse = new NURBS([
            VV(1, 0, 0, 1),
            VV(1, 1, 0, 1).times(SQRT1_2),
            VV(0, 1, 0, 1),
            VV(-1, 1, 0, 1).times(SQRT1_2),
            VV(-1, 0, 0, 1),
            VV(-1, -1, 0, 1).times(SQRT1_2),
            VV(0, -1, 0, 1),
        ], 2, [0, 0, 0, PI / 2, PI / 2, PI, PI, (3 * PI) / 2, (3 * PI) / 2, 2 * PI]);
        return unitSemiEllipse.transform(ellipse.matrix);
    }
    /**
     * Create a new NURBS from V3s, with optional weights.
     * @param points
     * @param degree
     * @param knots
     * @param weights
     */
    static fromV3s(points, degree, knots, weights = arrayFromFunction(points.length, () => 1)) {
        assert(points.length == weights.length);
        return new NURBS(points.map((p, i) => Vector.fromV3AndWeight(p, weights[i])), degree, knots);
    }
    isUniform(precision = 0) {
        const intervals = arrayFromFunction(this.knots.length - 1, (i) => this.knots[i + 1] - this.knots[i]);
        const [min, max] = minAndMax(intervals);
        return eq(min, max, precision);
    }
    /**
     * NURBS is a B spline if control points all have the same weight.
     */
    isBSpline(precision = 0) {
        const [minWeight, maxWeight] = minAndMax(this.points.map((p) => p.w));
        return eq(minWeight, maxWeight, precision);
    }
    /**
     * Whether this is a (rational) bezier curve.
     */
    isBezier(precision = 0) {
        if (this.degree + 1 != this.points.length)
            return false;
        const [min0, max0] = minAndMax(this.knots, 0, this.degree + 1);
        if (!eq(min0, max0, precision))
            return false;
        const [min1, max1] = minAndMax(this.knots, this.degree + 1);
        if (!eq(min1, max1, precision))
            return false;
        return true;
    }
    /**
     * Splits NURBS curve into rational bezier curves.
     * See https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/subdivision.html
     */
    getSegments() {
        const { knots, points, degree } = this;
        const result = [];
        const v = Vector.pack(points, new Float64Array(points.length * 4));
        const vectorFromV = (i) => new Vector(v.slice(i * 4, (i + 1) * 4));
        let k = degree + 1; // k = knot index we are duplicating
        while (k < knots.length - degree - 1) {
            const t = knots[k];
            const prevKnot = knots[k - 1];
            let s = 1; // s = multiplicity of the knot
            while (knots[k + 1] == t) {
                k++;
                s++;
            }
            const newNURBSPoints = new Array(degree + 1);
            // the first s + 1 points are identical to the current curve
            for (let i = 0; i < s + 1; i++) {
                newNURBSPoints[i] = vectorFromV(k - degree - s + i);
            }
            // we need to have multiplicity degree, so insert (degree - s) times
            for (let level = 1; level <= degree - s; level++) {
                for (let i = k - degree; i <= k - s - level; i++) {
                    const alpha = (t - prevKnot) / (knots[i + degree + 1] - prevKnot);
                    for (let dim = 0; dim < 4; dim++) {
                        v[i * 4 + dim] =
                            (1 - alpha) * v[i * 4 + dim] + alpha * v[(i + 1) * 4 + dim];
                    }
                }
                newNURBSPoints[s + level] = vectorFromV(k - degree);
            }
            const newNURBSKnots = arrayFromFunction((degree + 1) * 2, (i) => i < degree + 1 ? knots[k - s] : t);
            result.push(new NURBS(newNURBSPoints, degree, newNURBSKnots));
            k++;
        }
        // last curve
        const newNURBSPoints = arrayFromFunction(degree + 1, (i) => vectorFromV(points.length - degree - 1 + i));
        const newNURBSKnots = arrayFromFunction((degree + 1) * 2, (i) => i < degree + 1 ? knots[k - 1] : knots[k]);
        result.push(new NURBS(newNURBSPoints, degree, newNURBSKnots));
        return result;
    }
    split(t) {
        const { knots, points, degree } = this;
        assert(le(this.tMin, t) && le(t, this.tMax));
        let k = this.tInterval(t), s = 0; // s = multiplicity of the knot
        while (knots[k + 1] == t) {
            k++;
            s++;
        }
        const vectorFromV = (i) => new Vector(v.slice(i * 4, (i + 1) * 4));
        const leftPoints = new Array(k + 1 - s);
        // the first k + s + 1 points are identical to the current curve
        for (let i = 0; i < k + s - degree + 1; i++) {
            leftPoints[i] = this.points[i];
        }
        const rightPointsLength = points.length - (k - degree);
        const v = Vector.pack(points, new Float64Array(rightPointsLength * 4), k - degree);
        // we need to have multiplicity degree, so insert (degree - s) times
        for (let level = 1; level <= degree - s; level++) {
            for (let i = k - degree; i <= k - s - level; i++) {
                const alpha = (t - knots[i + level]) / (knots[i + degree + 1] - knots[i + level]);
                const j = i - (k - degree);
                for (let dim = 0; dim < 4; dim++) {
                    v[j * 4 + dim] =
                        (1 - alpha) * v[j * 4 + dim] + alpha * v[(j + 1) * 4 + dim];
                }
            }
            leftPoints[k - degree + level] = vectorFromV(0);
        }
        const leftKnots = knots.slice(0, k + degree + 2 - s);
        for (let i = 0; i < degree - s + 1; i++) {
            leftKnots[k - s + 1 + i] = t;
        }
        const rightKnots = knots.slice(k - degree);
        for (let i = 0; i < degree + 1; i++) {
            rightKnots[i] = t;
        }
        const rightPoints = arrayFromFunction(rightPointsLength, (i) => vArrGet(v, 4, i));
        return [
            new NURBS(leftPoints, degree, leftKnots),
            new NURBS(rightPoints, degree, rightKnots),
        ];
    }
    simplify() {
        assert(this.isBezier());
        if (3 == this.degree && this.isBSpline()) {
            return new BezierCurve(this.points[0].p3(), this.points[1].p3(), this.points[2].p3(), this.points[3].p3(), this.tMin, this.tMax);
        }
        else if (2 == this.degree) {
            const [P0, P1, P2] = this.points;
            const [p0, p1, p2] = this.points.map((p) => p.p3());
            const c = NURBS.simplifyUnit2(P0.w, P1.w, P2.w).transform(M4.forSys(p1.to(p0), p1.to(p2), undefined, p1));
            const [tMin, tMax] = [c.pointT(p0), c.pointT(p2)].sort();
            return c.withBounds(snap(tMin, c.tMin), snap(tMax, c.tMax));
        }
        else if (1 == this.degree) {
            return L3.throughPoints(this.points[0].p3(), this.points[1].p3());
        }
        else {
            return this;
        }
    }
    static simplifyUnit2(w0, w1, w2) {
        // see https://math.stackexchange.com/a/2794874/230980
        const delta = w0 * w2 - Math.pow(w1, 2);
        const cxy = (w0 * w2) / 2 / delta;
        const center = new V3(cxy, cxy, 0);
        const k = (Math.pow(w1, 2) + delta - 2 * w1 * sqrt(abs(delta))) / 2 / delta;
        const p = V3.X;
        const q = new V3(k, cxy, 0);
        // const q = new V3(cxy, k, 0)
        if (eq0(delta)) {
            return new ParabolaCurve(new V3(1 / 4, 1 / 4, 0), new V3(1, -1, 0), new V3(1, 1, 0), -0.5, 0.5);
        }
        else if (delta < 0) {
            // hyperbola
            return new HyperbolaCurve(center, center.to(p), center.to(q));
        }
        else {
            // ellipse
            return new EllipseCurve(center, center.to(p), center.to(q), 0);
        }
    }
    elevateDegreeBezier() {
        assert(this.isBezier());
        const newPoints = new Array(this.points.length + 1);
        newPoints[0] = this.points[0];
        newPoints[this.points.length] = this.points[this.points.length - 1];
        for (let i = 1; i < this.points.length; i++) {
            newPoints[i] = Vector.lerp(this.points[i], this.points[i - 1], i / (this.degree + 1));
        }
        const newKnots = NURBS.bezierKnots(this.degree + 1, this.knots[0], this.knots[this.degree + 1]);
        return new NURBS(newPoints, this.degree + 1, newKnots, this.tMin, this.tMax);
    }
    elevateDegree() {
        const segmentsElevated = this.getSegments().map((b) => b.elevateDegreeBezier());
        // stitch together the segments
        const newPoints = new Array(2 + segmentsElevated.length * this.degree);
        newPoints[0] = segmentsElevated[0].points[0];
        newPoints.last = segmentsElevated.last.points.last;
        for (let i = 0; i < segmentsElevated.length; i++) {
            for (let pi = 1; pi < segmentsElevated[i].points.length - 1; pi++) {
                newPoints[i * (segmentsElevated[0].points.length - 2) + pi] =
                    segmentsElevated[i].points[pi];
            }
        }
        const newKnots = new Array(newPoints.length + this.degree + 2);
        for (let i = 0; i < this.degree + 2; i++) {
            newKnots[i] = this.knots[0];
        }
        for (let i = 0; i < segmentsElevated.length; i++) {
            for (let pi = 1; pi < segmentsElevated[i].points.length - 1; pi++) {
                newKnots[i * (segmentsElevated[0].points.length - 2) + pi + this.degree + 1] = segmentsElevated[i].knots.last;
            }
        }
        newKnots[newKnots.length - 1] = this.knots.last;
        newKnots[newKnots.length - 2] = this.knots.last;
        let result = new NURBS(newPoints, this.degree + 1, newKnots, this.tMin, this.tMax);
        for (let i = 0; i < segmentsElevated.length - 1; i++) {
            let optimization;
            while ((optimization = result.removeKnot(segmentsElevated[i].knots.last))) {
                result = optimization;
            }
        }
        return result;
    }
    transform(m4) {
        return this.transform4(m4);
    }
    transform4(m4) {
        return new NURBS(this.points.map((p) => m4.timesVector(p)), this.degree, this.knots, this.tMin, this.tMax);
    }
    /**
     * Returns the index of the interval which contains the value t.
     */
    tInterval(t) {
        const { degree, knots } = this;
        for (let s = degree; s < knots.length - 1 - degree; s++) {
            if (t >= knots[s] && t <= knots[s + 1]) {
                return s;
            }
        }
        throw new Error(t + " " + knots);
    }
    static UnitCircle(sections = 2, tMin = 0, tMax = PI) {
        const dt = tMax - tMin;
        const tStep = dt / sections;
        const w = sin(PI / 2 - tStep / 2);
        console.log(tStep / 2 / DEG);
        // cos
        const r = 1 / cos(tStep / 2);
        const points = arrayFromFunction(sections * 2 + 1, (i) => {
            const t = lerp(tMin, tMax, i / 2 / sections);
            if (i % 2 == 0) {
                // control point on circle
                return VV(cos(t), sin(t), 0, 1);
            }
            else {
                return VV(r * w * cos(t), r * w * sin(t), 0, w);
            }
        });
        const knots = [];
        knots.push(tMin, tMin, tMin);
        for (let i = 0; i < sections - 1; i++) {
            const knot = lerp(tMin, tMax, (i + 1) / sections);
            knots.push(knot, knot);
        }
        knots.push(tMax, tMax, tMax);
        return new NURBS(points, 2, knots);
    }
    debugInfo() {
        return {
            points: [
                ...this.knots.slice(this.degree, -this.degree).map((t) => this.at(t)),
                ...this.points.map((p) => p.p3()),
            ],
            lines: this.points.flatMap((p, i, ps) => ps[i + 1] ? [p.p3(), ps[i + 1].p3()] : []),
        };
    }
    isTsWithPlane(planeWC) {
        const { knots, degree, points } = this;
        const controlPointTs = [
            knots[degree],
            ...points
                .slice(1, -1)
                .map((p, i) => this.closestTToPoint(p.p3(), undefined, knots[i + 3], knots[i + degree])),
            knots[knots.length - degree - 1],
        ];
        const result = [];
        for (let i = 0; i < this.points.length - 1; i++) {
            const findClosest = (startT) => {
                console.log("startT", startT);
                // try {
                const f = (t) => {
                    const [p, dt] = this.ptDtDdt(t);
                    return [planeWC.distanceToPointSigned(p), planeWC.normal1.dot(dt)];
                };
                let t = newtonIterateWithDerivative2(f, startT, 8, this.tMin, this.tMax);
                let [distanceAtT, distanceDtAtT] = undefined === t ? [undefined, undefined] : f(t);
                if (t === undefined || !eq0(distanceAtT) || eq0(distanceDtAtT)) {
                    t = newtonIterateWithDerivative2((t) => {
                        const [, dt, ddt] = this.ptDtDdt(t);
                        return [planeWC.normal1.dot(dt), planeWC.normal1.dot(ddt)];
                    }, startT, 8, this.tMin, this.tMax);
                }
                [distanceAtT, distanceDtAtT] = undefined === t ? [] : f(t);
                if (undefined !== t &&
                    eq0(distanceAtT) &&
                    !result.some((r) => eq(r, t))) {
                    result.push(t);
                }
            };
            const a = this.points[i].p3();
            const b = this.points[i + 1].p3();
            const ad = snap0(planeWC.distanceToPointSigned(a));
            const bd = snap0(planeWC.distanceToPointSigned(b));
            if (ad * bd < 0) {
                const startT = lerp(controlPointTs[i], controlPointTs[i + 1], ad / (ad - bd));
                findClosest(startT);
            }
            else if (0 == bd) {
                findClosest(this.closestTToPoint(b, controlPointTs[i + 1]));
            }
        }
        return result;
    }
    isInfosWithCurve(curveWC) {
        if (curveWC instanceof L3) {
            return this.isInfosWithLine(curveWC.anchor, curveWC.dir1);
        }
        return super.isInfosWithCurve(curveWC);
    }
    isInfosWithLine(anchor, dir) {
        const thisPlane = P3.fromPoints(this.points.map((p) => p.p3()));
        const l = L3.anchorDirection(anchor, dir);
        const maxDistanceToPlane = this.points
            .map((p) => thisPlane.distanceToPoint(p.p3()))
            .max();
        const thisIsPlanar = eq0(maxDistanceToPlane);
        if (thisIsPlanar && !thisPlane.containsLine(l)) {
            const [t] = l.isTsWithPlane(thisPlane);
            if (undefined === t)
                return [];
            const p = l.at(t);
            return this.containsPoint(p)
                ? [{ tThis: this.pointT(p), tOther: L3.pointT(anchor, dir, p), p }]
                : [];
        }
        else {
            const thisTs = this.isTsWithPlane(P3.normalOnAnchor(thisPlane.normal1.cross(dir), anchor));
            const infos = thisTs.map((tThis) => {
                const p = this.at(tThis);
                return { tThis, tOther: L3.pointT(anchor, dir, p), p };
            });
            return thisIsPlanar
                ? infos
                : infos.filter((info) => L3.containsPoint(anchor, dir, info.p));
        }
    }
    roots() {
        console.log(this.tMin, this.tMax);
        arraySamples(this.tMin, this.tMax, 30).forEach((t) => {
            console.log(t + "," + this.tangentAt(t).z);
        });
        const result = [[], [], []];
        for (let i = 0; i < this.points.length - 1; i++) {
            const findClosest = (startT, d) => {
                console.log("d", d, "startT", startT);
                // try {
                const root = newtonIterateWithDerivative2((t) => {
                    const [, dt, ddt] = this.ptDtDdt(t);
                    return [dt.e(d), ddt.e(d)];
                }, startT, 8, this.tMin, this.tMax);
                if (undefined !== root) {
                    result[d].push(root);
                }
                console.log("d", d, "startT", startT, "root", root);
            };
            const a = this.points[i].p3();
            const b = this.points[i + 1].p3();
            const ab = a.to(b);
            for (let d = 0; d < 3; d++) {
                if (0 !== i && eq0(ab.e(d))) {
                    const startT = lerp(this.knots[i], this.knots[i + this.degree + 2], 0.5);
                    findClosest(startT, d);
                }
                else if (i < this.points.length - 2) {
                    const bc = b.to(this.points[i + 2].p3());
                    if (!eq0(bc.e(d)) && ab.e(d) * bc.e(d) < 0) {
                        findClosest(this.closestTToPoint(b, this.guessTClosestToControlPoint(i + 1)), d);
                    }
                }
            }
        }
        console.log(result);
        return result;
    }
    //getAABB() {
    //	return new AABB().addPoints(this.points.map(p => p.p3()))
    //}
    /**
     * Rough approximation of t param for points closest to control point.
     */
    guessTClosestToControlPoint(pointIndex) {
        return lerp(this.knots[pointIndex], this.knots[pointIndex + this.degree + 1], 0.5);
    }
    likeCurve(curve) {
        return (this == curve ||
            (hasConstructor(curve, NURBS) &&
                this.degree === curve.degree &&
                this.points.every((p, i) => p.like(curve.points[i])) &&
                this.knots.every((k, i) => eq(k, curve.knots[i]))));
    }
    isColinearTo(curve) {
        throw new Error("This doesn't even make sense.");
    }
}
NURBS.EX2D = NURBS.fromV3s([
    V(51, 141),
    V(11, 76),
    V(29, 32),
    V(46, 102),
    V(74, 148),
    V(189, 107),
    V(56, 10),
    V(206, 10),
    V(211, 98),
    V(195, 141),
    V(139, 148),
], 4);
NURBS.EX3D = new NURBS([
    VV(94, 0, -34, 1),
    VV(69, 57, 45, 0.5),
    VV(-20, 44, 91, 1),
    VV(-89, -13, 47, 0.5),
    VV(-56, -97, -7, 1),
    VV(34, -83, -54, 0.5),
    VV(112, -53, 16, 1),
    VV(79, 30, 70, 0.5),
    VV(-2, -9, 141, 1),
    VV(-80, -40, 72, 0.5),
    VV(-38, -150, 43, 1),
    VV(43, -110, -29, 0.5),
    VV(130, -106, 65, 1),
], 2, [-12, -12, -12, -8, -8, -4, -4, 0, 0, 4, 4, 8, 8, 12, 12, 12]);
NURBS.prototype.tIncrement = 1 / 128;
function minAndMax(arr, start = 0, end = arr.length) {
    let min = Infinity, max = -Infinity;
    for (let i = start; i < end; i++) {
        if (min > arr[i])
            min = arr[i];
        if (max < arr[i])
            max = arr[i];
    }
    return [min, max];
}

/**
 * Plane x DOT this.normal1 = this.w
 */
class P3 extends Transformable {
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
        assertVectors(normal1);
        assertNumbers(w);
        assert(normal1.hasLength(1), "normal1.hasLength(1)" + normal1);
    }
    get anchor() {
        return this.normal1.times(this.w);
    }
    static throughPoints(a, b, c) {
        assertVectors(a, b, c);
        const n1 = b.minus(a).cross(c.minus(a)).unit();
        return new P3(n1, n1.dot(a));
    }
    static normalOnAnchor(normal, anchor) {
        assertVectors(normal, anchor);
        const n1 = normal.unit();
        return new this(n1, n1.dot(anchor));
    }
    /**
     * Create a plane which intersects the X, Y and Z axes at the specified offsets.
     * x/x0 + y/y0 + y/y0 = 1
     */
    static forAxisIntercepts(x0, y0, z0) {
        assertNumbers(x0, y0, z0);
        const normal = new V3(1 / x0, 1 / y0, 1 / z0);
        return new P3(normal.unit(), normal.length());
    }
    /**
     * Create a plane containing `anchor` and extending in directions `v0` and `v1`.
     * `v0` and `v1` may not be parallel.
     * @param anchor
     * @param v0
     * @param v1
     */
    static forAnchorAndPlaneVectors(anchor, v0, v1) {
        assertVectors(anchor, v0, v1);
        assert(!v0.isParallelTo(v1));
        return this.normalOnAnchor(v0.cross(v1), anchor);
    }
    /**
     * Create a plane which contains botha point and a line. The point may not lie on the line.
     * @param p
     * @param line
     */
    static forPointAndLine(p, line) {
        return this.forAnchorAndPlaneVectors(line.anchor, line.dir1, line.anchor.to(p));
    }
    /**
     * ax + by + cz + d = 0
     */
    static forABCD(a, b, c, d) {
        const normalLength = Math.hypot(a, b, c);
        if (eq0(normalLength))
            return undefined;
        return new P3(new V3(a / normalLength, b / normalLength, c / normalLength), -d / normalLength);
    }
    static vanishingPlane(m4) {
        return P3.forABCD(m4.m[12], m4.m[13], m4.m[14], m4.m[15]);
    }
    static forAABB(aabb, distance = 0) {
        return [
            new P3(V3.X, aabb.max.x + distance),
            new P3(V3.X.negated(), -aabb.min.x - distance),
            new P3(V3.Y, aabb.max.y + distance),
            new P3(V3.Y.negated(), -aabb.min.y - distance),
            new P3(V3.Z, aabb.max.z + distance),
            new P3(V3.Z.negated(), -aabb.min.z - distance),
        ];
    }
    // Fit a plane to a collection of points.
    // Fast, and accurate to within a few degrees.
    // Returns None if the points do not span a plane.
    static fromPoints(points) {
        const n = points.length;
        if (n < 3) {
            return undefined;
        }
        const centroid = V3.add(...points).div(n);
        // Calculate full 3x3 covariance matrix, excluding symmetries:
        let xx = 0.0;
        let xy = 0.0;
        let xz = 0.0;
        let yy = 0.0;
        let yz = 0.0;
        let zz = 0.0;
        for (const p of points) {
            const r = p.minus(centroid);
            xx += r.x * r.x;
            xy += r.x * r.y;
            xz += r.x * r.z;
            yy += r.y * r.y;
            yz += r.y * r.z;
            zz += r.z * r.z;
        }
        xx /= n;
        xy /= n;
        xz /= n;
        yy /= n;
        yz /= n;
        zz /= n;
        let weighted_dir = V3.O;
        {
            const det_x = yy * zz - yz * yz;
            const axis_dir = new V3(det_x, xz * yz - xy * zz, xy * yz - xz * yy);
            let weight = det_x * det_x;
            if (weighted_dir.dot(axis_dir) < 0.0) {
                weight = -weight;
            }
            weighted_dir = weighted_dir.plus(axis_dir.times(weight));
        }
        {
            const det_y = xx * zz - xz * xz;
            const axis_dir = new V3(xz * yz - xy * zz, det_y, xy * xz - yz * xx);
            let weight = det_y * det_y;
            if (weighted_dir.dot(axis_dir) < 0.0) {
                weight = -weight;
            }
            weighted_dir = weighted_dir.plus(axis_dir.times(weight));
        }
        {
            const det_z = xx * yy - xy * xy;
            const axis_dir = new V3(xy * yz - xz * yy, xy * xz - yz * xx, det_z);
            let weight = det_z * det_z;
            if (weighted_dir.dot(axis_dir) < 0.0) {
                weight = -weight;
            }
            weighted_dir = weighted_dir.plus(axis_dir.times(weight));
        }
        const normal = weighted_dir.unit();
        return P3.normalOnAnchor(normal, centroid);
    }
    axisIntercepts() {
        const w = this.w, n = this.normal1;
        return new V3(w / n.x, w / n.y, w / n.z);
    }
    isCoplanarToPlane(plane) {
        assertInst(P3, plane);
        return this.like(plane) || this.likeFlipped(plane);
    }
    like(plane) {
        assertInst(P3, plane);
        return eq(this.w, plane.w) && this.normal1.like(plane.normal1);
    }
    likeFlipped(plane) {
        assertInst(P3, plane);
        return eq(this.w, -plane.w) && this.normal1.like(plane.normal1.negated());
    }
    /**
     * True iff plane.normal1 is equal to this.normal1 or it's negation.
     *
     */
    isParallelToPlane(plane) {
        assertInst(P3, plane);
        return eq(1, Math.abs(this.normal1.dot(plane.normal1)));
    }
    isParallelToLine(line) {
        assertInst(L3, line);
        return eq0(this.normal1.dot(line.dir1));
    }
    isPerpendicularToLine(line) {
        assertInst(L3, line);
        // this.normal1 || line.dir1
        return eq(1, Math.abs(this.normal1.dot(line.dir1)));
    }
    isPerpendicularToPlane(plane) {
        assertInst(P3, plane);
        return eq0(this.normal1.dot(plane.normal1));
    }
    toSource() {
        return callsce("new P3", this.normal1, this.w);
    }
    translated(offset) {
        return new P3(this.normal1, this.w + offset.dot(this.normal1));
    }
    transform(m4) {
        // See https://stackoverflow.com/questions/7685495/transforming-a-3d-plane-using-a-4x4-matrix
        // See http://www.songho.ca/opengl/gl_normaltransform.html
        // with homogeneous coordinates, the hessian normal form of this plane is
        // (p, 1) * (normal1, -w) = 0
        // transformation: (m4^-1 * (p, 1)) DOT (normal1, -w) = 0
        // => (p, 1) DOT ((m4^-T) * (normal1, -w)) = 0
        // (validity of the above transformation is easily seen by expanding the matrix multiplication and dot product)
        // hence, (newNormal, newW) = (m4^-T) * (normal1, -w)
        // we divide both newNormal and newW by newNormal.length() to normalize the normal vector
        const m4InversedTransposed = M4.transpose(M4.inverse(m4, M4.temp0), M4.temp1);
        const [nx, ny, nz] = this.normal1;
        const newNormal = m4InversedTransposed.timesVector(VV(nx, ny, nz, -this.w));
        return P3.forABCD(newNormal.x, newNormal.y, newNormal.z, newNormal.w);
    }
    distanceToLine(line) {
        assertInst(L3, line);
        if (!this.isParallelToLine(line)) {
            return this.distanceToPoint(line.anchor);
        }
        else {
            return 0;
        }
    }
    containsPoint(x) {
        assertVectors(x);
        return eq(this.w, this.normal1.dot(x));
    }
    containsLine(line) {
        assertInst(L3, line);
        return this.containsPoint(line.anchor) && this.isParallelToLine(line);
    }
    distanceToPointSigned(point) {
        assertInst(V3, point);
        return this.normal1.dot(point) - this.w;
    }
    distanceToPoint(point) {
        assertInst(V3, point);
        return Math.abs(this.normal1.dot(point) - this.w);
    }
    intersectionWithLine(line) {
        return line.intersectionWithPlane(this);
    }
    intersectionWithPlane(plane) {
        assertInst(P3, plane);
        /*
    
             this: n0 * x = w0
             plane: n1 * x = w1
             plane perpendicular to both which goes through origin:
             n2 := n0 X x1
             n2 * x = 0
             */
        if (this.isParallelToPlane(plane)) {
            return undefined;
        }
        /*
             var n0 = this.normal1, n1 = plane.normal1, n2 = n0.cross(n1).unit(), m = M4.forSys(n0, n1, n2)
             var x0 = this.anchor, x1 = plane.anchor, x2 = V3.O
             var p = n2.times(x2.dot(n2))
             .plus(n1.cross(n2).times(x0.dot(n0)))
             .plus(n2.cross(n0).times(x1.dot(n1)))
             .div(m.determinant())
             */
        const n0 = this.normal1, n1 = plane.normal1, n2 = n0.cross(n1).unit();
        const p = M4.forRows(n0, n1, n2)
            .inversed()
            .transformVector(new V3(this.w, plane.w, 0));
        return new L3(p, n2);
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
        if (curve instanceof L3) {
            return this.containsLine(curve);
        }
        else if (curve instanceof EllipseCurve ||
            curve instanceof HyperbolaCurve ||
            curve instanceof ParabolaCurve) {
            return (this.containsPoint(curve.center) &&
                this.normal1.isParallelTo(curve.normal));
        }
        else if (curve instanceof BezierCurve) {
            return curve.points.every((p) => this.containsPoint(p));
        }
        else {
            throw new Error("" + curve);
        }
    }
    equals(obj) {
        return (hasConstructor(obj, P3) &&
            this.normal1.equals(obj.normal1) &&
            this.w == obj.w);
    }
    hashCode() {
        return (this.normal1.hashCode() * 31) | (0 + floatHashCode(this.w));
    }
}
P3.YZ = new P3(V3.X, 0);
P3.ZX = new P3(V3.Y, 0);
P3.XY = new P3(V3.Z, 0);

class Surface extends Transformable {
    static loopContainsPointGeneral(loop, pWC, testLine, lineOut) {
        const testPlane = P3$1.normalOnAnchor(lineOut, pWC);
        // edges colinear to the testing line; these will always be counted as "inside" relative to the testing line
        const colinearEdges = loop.map((edge) => edge.colinearToLine(testLine));
        let inside = false;
        function logIS(isP) {
            const isT = testLine.pointT(isP);
            if (eq0(isT)) {
                return true;
            }
            else if (isT > 0) {
                inside = !inside;
            }
            return false;
        }
        for (let edgeIndex = 0; edgeIndex < loop.length; edgeIndex++) {
            const edge = loop[edgeIndex];
            const nextEdgeIndex = (edgeIndex + 1) % loop.length, nextEdge = loop[nextEdgeIndex];
            //console.log(edge.toSource()) {p:V(2, -2.102, 0),
            if (colinearEdges[edgeIndex]) {
                const lineAT = testLine.pointT(edge.a), lineBT = testLine.pointT(edge.b);
                if (Math.min(lineAT, lineBT) <= NLA_PRECISION &&
                    -NLA_PRECISION <= Math.max(lineAT, lineBT)) {
                    return PointVsFace.ON_EDGE;
                }
                // edge colinear to intersection
                const nextInside = colinearEdges[nextEdgeIndex] ||
                    dotCurve2$1(nextEdge.curve, nextEdge.aT, lineOut, sign(nextEdge.deltaT())) < 0;
                if (!nextInside) {
                    if (logIS(edge.b))
                        return PointVsFace.ON_EDGE;
                }
            }
            else {
                for (const edgeT of edge.edgeISTsWithPlane(testPlane)) {
                    if (edgeT == edge.bT) {
                        if (!testLine.containsPoint(edge.b))
                            continue;
                        // endpoint lies on intersection line
                        if (edge.b.like(pWC)) {
                            // TODO: refactor, dont check for different sides, just logIs everything
                            return PointVsFace.ON_EDGE;
                        }
                        const edgeInside = dotCurve2$1(edge.curve, edge.bT, lineOut, -sign(edge.deltaT())) < 0;
                        const nextInside = colinearEdges[nextEdgeIndex] ||
                            dotCurve2$1(nextEdge.curve, nextEdge.aT, lineOut, sign(nextEdge.deltaT())) < 0;
                        if (edgeInside != nextInside) {
                            if (logIS(edge.b))
                                return PointVsFace.ON_EDGE;
                        }
                    }
                    else if (edgeT != edge.aT) {
                        const p = edge.curve.at(edgeT);
                        if (!testLine.containsPoint(p))
                            continue;
                        // edge crosses line, neither starts nor ends on it
                        if (logIS(p))
                            return PointVsFace.ON_EDGE;
                        // TODO: tangents?
                    }
                }
            }
        }
        return inside ? PointVsFace.INSIDE : PointVsFace.OUTSIDE;
    }
    static loopContainsPointEllipse(loop, pWC, testLine, pWCT) {
        const lineOut = testLine.normal;
        const testPlane = P3$1.normalOnAnchor(testLine.normal, pWC);
        const colinearEdges = loop.map((edge) => testLine.isColinearTo(edge.curve));
        let inside = false;
        if (undefined === pWCT) {
            pWCT = testLine.pointT(pWC);
        }
        const pT = pWCT;
        function logIS(isP) {
            const isT = testLine.pointT(isP);
            if (eq(pT, isT)) {
                return true;
            }
            else if (pT < isT && le(isT, PI)) {
                inside = !inside;
            }
            return false;
        }
        for (let edgeIndex = 0; edgeIndex < loop.length; edgeIndex++) {
            const edge = loop[edgeIndex];
            const nextEdgeIndex = (edgeIndex + 1) % loop.length, nextEdge = loop[nextEdgeIndex];
            //console.log(edge.toSource()) {p:V(2, -2.102, 0),
            if (colinearEdges[edgeIndex]) {
                let edgeT;
                if (edge.curve.containsPoint(pWC) &&
                    le(edge.minT, (edgeT = edge.curve.pointT(pWC))) &&
                    le(edgeT, edge.maxT)) {
                    return PointVsFace.ON_EDGE;
                }
                // edge colinear to intersection
                const nextInside = colinearEdges[nextEdgeIndex] ||
                    dotCurve2$1(nextEdge.curve, nextEdge.aT, lineOut, sign(nextEdge.deltaT())) < 0;
                if (!nextInside && testLine.containsPoint(edge.b)) {
                    if (logIS(edge.b))
                        return PointVsFace.ON_EDGE;
                }
            }
            else {
                for (const edgeT of edge.edgeISTsWithPlane(testPlane)) {
                    if (edgeT == edge.bT) {
                        if (!testLine.containsPoint(edge.b))
                            continue;
                        // endpoint lies on intersection testLine
                        const edgeInside = dotCurve2$1(edge.curve, edge.bT, lineOut, -sign(edge.deltaT())) < 0;
                        const nextInside = colinearEdges[nextEdgeIndex] ||
                            dotCurve2$1(nextEdge.curve, nextEdge.aT, lineOut, sign(nextEdge.deltaT())) < 0;
                        if (edgeInside != nextInside) {
                            if (logIS(edge.b))
                                return PointVsFace.ON_EDGE;
                        }
                    }
                    else if (edgeT != edge.aT) {
                        const p = edge.curve.at(edgeT);
                        if (!testLine.containsPoint(p))
                            continue;
                        // edge crosses testLine, neither starts nor ends on it
                        if (logIS(p))
                            return PointVsFace.ON_EDGE;
                        // TODO: tangents?
                    }
                }
            }
        }
        return inside ? PointVsFace.INSIDE : PointVsFace.OUTSIDE;
    }
    toString() {
        return this.toSource();
    }
    toSource(rounder = (x) => x) {
        return callsce.call(undefined, "new " + this.constructor.name, ...this.getConstructorParameters());
    }
    /**
     * Return points which would touch AABB. Doesnt include borders due to paramtetric bounds, for example.
     */
    getExtremePoints() {
        return [];
    }
    isCurvesWithSurface(surface) {
        return surface.isCurvesWithSurface(this); //.map(curve => curve.reversed())
    }
    containsCurve(curve) {
        if (curve instanceof PPCurve$1) {
            if (this.equals(curve.parametricSurface1) ||
                this.equals(curve.parametricSurface2)) {
                return true;
            }
        }
        if (curve instanceof ImplicitCurve$1) {
            for (let i = ceil(curve.tMin) + 1; i <= floor(curve.tMax) - 1; i++) {
                if (!this.containsPoint(curve.points[i])) {
                    return false;
                }
            }
            return true;
        }
        return false;
    }
    flipped2(doFlip) {
        return doFlip ? this.flipped() : this;
    }
    clipCurves(curves) {
        return curves;
    }
    equals(obj) {
        return (this === obj ||
            (this.constructor === obj.constructor &&
                arrayEquals(this.getConstructorParameters(), obj.getConstructorParameters())));
    }
    hashCode() {
        return arrayHashCode(this.getConstructorParameters());
    }
    zDirVolume(allEdges) {
        return this.visit(ZDirVolumeVisitor$1, allEdges);
    }
    calculateArea(allEdges) {
        return this.visit(CalculateAreaVisitor$1, allEdges);
    }
}
var PointVsFace;
(function (PointVsFace) {
    PointVsFace[PointVsFace["INSIDE"] = 0] = "INSIDE";
    PointVsFace[PointVsFace["OUTSIDE"] = 1] = "OUTSIDE";
    PointVsFace[PointVsFace["ON_EDGE"] = 2] = "ON_EDGE";
})(PointVsFace || (PointVsFace = {}));
class ImplicitSurface extends Surface {
    static is(obj) {
        return obj.implicitFunction && obj.didp;
    }
}

class ParametricSurface extends Surface$1 {
    constructor(uMin, uMax, vMin, vMax) {
        super();
        this.uMin = uMin;
        this.uMax = uMax;
        this.vMin = vMin;
        this.vMax = vMax;
        assertNumbers(uMin, uMax, vMin, vMax);
        assert(uMin < uMax);
        assert(vMin < vMax);
        assert(((x) => x[x.length - 4])(this.getConstructorParameters()) == this.uMin, this.getConstructorParameters(), this.uMin);
    }
    static isCurvesParametricImplicitSurface(ps, is, uStep, vStep = uStep, curveStepSize) {
        const pf = ps.pUVFunc(), icc = is.implicitFunction();
        const dpdu = ps.dpdu();
        const dpdv = ps.dpdv();
        const didp = is.didp.bind(is);
        const ist = (x, y) => icc(pf(x, y));
        const didu = (u, v) => didp(pf(u, v)).dot(dpdu(u, v));
        const didv = (u, v) => didp(pf(u, v)).dot(dpdv(u, v));
        const mf = MathFunctionR2R$1.forFFxFy(ist, didu, didv);
        const curves = Curve$1.breakDownIC(mf, ps, uStep, vStep, curveStepSize, (u, v) => is.containsPoint(pf(u, v))).map(({ points, tangents }, i) => PICurve$1.forParametricPointsTangents(ps, is, points, tangents, curveStepSize));
        return curves;
    }
    static isCurvesParametricParametricSurface(ps1, ps2, s1Step, t1Step = s1Step, curveStepSize) {
        return breakDownPPCurves$1(ps1, ps2, s1Step, t1Step, curveStepSize);
    }
    static is(obj) {
        return obj.pUVFunc;
    }
    pUV(u, v) {
        return this.pUVFunc()(u, v);
    }
    pUVFunc() {
        return this.pUV.bind(this);
    }
    uvP(pWC) {
        return this.uvPFunc()(pWC);
    }
    uvPFunc() {
        return this.uvP.bind(this);
    }
    bounds(u, v) {
        return this.uMin <= u && u <= this.uMax && this.vMin <= v && v <= this.vMax;
    }
    /**
     * Positive values are inside bounds.
     */
    boundsSigned(u, v) {
        return min(u - this.uMin, this.uMax - u, v - this.vMin, this.vMax - v);
    }
    normalP(p) {
        const pmPoint = this.uvPFunc()(p);
        return this.normalUV(pmPoint.x, pmPoint.y);
    }
    normalUVFunc() {
        return this.normalUV.bind(this);
    }
    normalUV(u, v) {
        return this.normalUVFunc()(u, v);
    }
    parametersValid(u, v) {
        return between(u, this.uMin, this.uMax) && between(v, this.vMin, this.vMax);
    }
    toMesh(uStep = this.uStep, vStep = this.vStep) {
        assert(isFinite(this.vMin) &&
            isFinite(this.vMax) &&
            isFinite(this.uMin) &&
            isFinite(this.uMax));
        assert(isFinite(uStep) && isFinite(vStep));
        return Mesh.parametric(this.pUVFunc(), this.normalUVFunc(), this.uMin, this.uMax, this.vMin, this.vMax, ceil((this.uMax - this.uMin) / uStep), ceil((this.vMax - this.vMin) / vStep));
    }
    isCurvesWithImplicitSurface(is, uStep, vStep, stepSize) {
        return ParametricSurface.isCurvesParametricImplicitSurface(this, is, uStep, vStep, stepSize);
    }
    edgeLoopCCW(contour) {
        const ptpF = this.uvPFunc();
        return isCCW(contour.flatMap((e) => e.getVerticesNo0()).map((v) => ptpF(v)), V3.Z);
    }
    like(object) {
        if (!this.isCoplanarTo(object))
            return false;
        // normals need to point in the same direction (outwards or inwards) for both
        const pSMinTMin = this.pUVFunc()(this.uMin, this.vMin);
        const thisNormal = this.normalUVFunc()(this.uMin, this.vMin);
        const otherNormal = object.normalP(pSMinTMin);
        return 0 < thisNormal.dot(otherNormal);
    }
    getApproxAABB() {
        const result = new AABB();
        result.addPoints(this.getExtremePoints());
        const ps = [V(0, 0), V(0, 1), V(1, 0), V(1, 1), V(0.5, 0.5)].map((p) => this.pUV(lerp(this.uMin, this.uMax, p.x), lerp(this.vMin, this.vMax, p.y)));
        result.addPoints(ps);
        return result;
    }
}

class ConicSurface extends ParametricSurface$1 {
    /**
     * returns new cone C = {apex + f1 * z * cos(d) + f2 * z * sin(d) + f3 * z | -PI <= d <= PI, 0 <= z}
     * @param f1
     * @param f2
     * @param dir Direction in which the cone opens. The ellipse spanned by f1, f2 is contained at (apex + f1).
     */
    constructor(center, f1, f2, dir, uMin = 0, uMax = PI, vMin = 0, vMax = 16) {
        super(uMin, uMax, vMin, vMax);
        this.center = center;
        this.f1 = f1;
        this.f2 = f2;
        this.dir = dir;
        assertVectors(center, f1, f2, dir);
        assert(0 <= vMin);
        this.matrix = M4.forSys(f1, f2, dir, center);
        this.matrixInverse = this.matrix.inversed();
        this.normalDir = sign(this.f1.cross(this.f2).dot(this.dir));
        this.pLCNormalWCMatrix = this.matrix
            .as3x3()
            .inversed()
            .transposed()
            .scale(this.normalDir);
    }
    pointFoot(pWC, startU, startV) {
        if (undefined === startU || undefined === startV) {
            // similar to uvP
            const pLC = this.matrixInverse.transformPoint(pWC);
            const angle = pLC.angleXY();
            if (undefined === startU) {
                startU = angle < -PI / 2 ? angle + TAU : angle;
            }
            if (undefined === startV) {
                startV = pLC.z + (pLC.lengthXY() - pLC.z) * SQRT1_2;
            }
        }
        const f = ([u, v]) => {
            const pUVToPWC = this.pUV(u, v).to(pWC);
            return [this.dpdu()(u, v).dot(pUVToPWC), this.dpdv()(u).dot(pUVToPWC)];
        };
        const { 0: x, 1: y } = newtonIterate(f, [startU, startV]);
        return new V3(x, y, 0);
    }
    get apex() {
        return this.center;
    }
    static atApexThroughEllipse(apex, ellipse, uMin, uMax, vMin, vMax) {
        assertVectors(apex);
        assertInst(EllipseCurve$1, ellipse);
        return new ConicSurface(apex, ellipse.f1, ellipse.f2, apex.to(ellipse.center), uMin, uMax, vMin, vMax);
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
        return pqFormula(b / a, c / a).filter((t) => 0 < az + t * dz);
    }
    // calculate intersection of plane ax + cz = d and cone x² + y² = z²
    static unitISPlane(a, c, d) {
        if (eq0(c)) {
            // plane is "vertical", i.e. parallel to Y and Z axes
            assert(!eq0(a)); // normal would be zero, which is invalid
            // z² - y² = d²/a²
            if (eq0(d)) {
                // d = 0 => z² - y² = 0 => z² = y² => z = y
                // plane goes through origin/V3.O
                return [
                    new L3$1(V3.O, new V3(0, -SQRT1_2, -SQRT1_2), undefined, 0),
                    new L3$1(V3.O, new V3(0, -SQRT1_2, SQRT1_2), 0),
                ];
            }
            else {
                // hyperbola
                const center = new V3(d / a, 0, 0);
                const f1 = new V3(0, 0, abs(d / a)); // abs, because we always want the hyperbola to be pointing up
                const f2 = new V3(0, d / a, 0);
                return [new HyperbolaCurve$1(center, f1, f2)];
            }
        }
        else {
            // c != 0
            const aa = a * a, cc = c * c;
            if (eq0(d)) {
                // ax + cz = d => x = d - cz / a => x² = d² - 2cdz/a + c²z²/a²
                // x² + y² = z²
                // => d² - 2cdz/a + c²z²/a² + y² = z²
                if (eq(aa, cc)) {
                    return [new L3$1(V3.O, new V3(c, 0, -a).unit())];
                }
                else if (aa < cc) {
                    throw new Error("intersection is single point V3.O");
                }
                else if (aa > cc) {
                    return [
                        new L3$1(V3.O, new V3(c, sqrt(aa - cc), -a).unit()),
                        new L3$1(V3.O, new V3(c, -sqrt(aa - cc), -a).unit()),
                    ];
                }
            }
            else {
                if (eq(aa, cc)) {
                    // parabola
                    const parabolaVertex = new V3(d / 2 / a, 0, d / 2 / c);
                    const parabolaVertexTangentPoint = new V3(d / 2 / a, d / c, d / 2 / c);
                    const p2 = new V3(0, 0, d / c);
                    const f2 = p2.minus(parabolaVertex);
                    return [
                        new ParabolaCurve$1(parabolaVertex, parabolaVertexTangentPoint.minus(parabolaVertex), f2.z < 0 ? f2.negated() : f2),
                    ];
                }
                else if (aa < cc) {
                    // ellipse
                    const center = new V3((-a * d) / (cc - aa), 0, (d * c) / (cc - aa));
                    if (center.z < 0) {
                        return [];
                    }
                    const p1 = new V3(d / (a - c), 0, -d / (a - c));
                    const p2 = new V3((-a * d) / (cc - aa), d / sqrt(cc - aa), (d * c) / (cc - aa));
                    return [
                        new EllipseCurve$1(center, center.to(p1), center.to(p2), -PI, PI),
                    ];
                }
                else if (aa > cc) {
                    // hyperbola
                    const center = new V3((-a * d) / (cc - aa), 0, (d * c) / (cc - aa));
                    // const p1 = new V3(d / (a - c), 0, -d / (a - c))
                    // const p2 = new V3(-a * d / (cc - aa), d / sqrt(aa - cc), d * c / (cc - aa))
                    // const f1 = center.to(p1)
                    const f1 = new V3((d * c) / (aa - cc), 0, (-d * a) / (aa - cc));
                    const f2 = new V3(0, d / sqrt(aa - cc), 0);
                    return [new HyperbolaCurve$1(center, f1.z > 0 ? f1 : f1.negated(), f2)];
                }
            }
        }
        throw new Error("???");
    }
    equals(obj) {
        return (this == obj ||
            (Object.getPrototypeOf(this) == Object.getPrototypeOf(obj) &&
                this.center.equals(obj.center) &&
                this.f1.equals(obj.f1) &&
                this.f2.equals(obj.f2) &&
                this.dir.equals(obj.dir)));
    }
    like(object) {
        if (!this.isCoplanarTo(object))
            return false;
        // normals need to point in the same direction (outwards or inwards) for both
        return this.normalDir == object.normalDir;
    }
    getVectors() {
        return [
            { anchor: this.center, dir1: this.dir },
            { anchor: this.center.plus(this.dir), dir1: this.f1 },
            { anchor: this.center.plus(this.dir), dir1: this.f2 },
        ];
    }
    getSeamPlane() {
        return P3$1.forAnchorAndPlaneVectors(this.center, this.f1, this.dir);
    }
    loopContainsPoint(contour, p) {
        assertVectors(p);
        const line = this.center.like(p)
            ? new L3$1(p, this.matrix.transformVector(new V3(0, 1, 1)).unit())
            : L3$1.throughPoints(p, this.apex);
        const lineOut = line.dir1.cross(this.dir);
        return Surface$1.loopContainsPointGeneral(contour, p, line, lineOut);
    }
    getConstructorParameters() {
        return [
            this.center,
            this.f1,
            this.f2,
            this.dir,
            this.uMin,
            this.uMax,
            this.vMin,
            this.vMax,
        ];
    }
    isTsForLine(line) {
        // transforming line manually has advantage that dir1 will not be renormalized,
        // meaning that calculated values t for lineLC are directly transferable to line
        const anchorLC = this.matrixInverse.transformPoint(line.anchor);
        const dirLC = this.matrixInverse.transformVector(line.dir1);
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
        return this.containsEllipse(new EllipseCurve$1(surface.center.plus(surface.dir), surface.f1, surface.f2));
    }
    containsEllipse(ellipse) {
        const ellipseLC = ellipse.transform(this.matrixInverse);
        if (ellipseLC.center.z < 0) {
            return false;
        }
        const { f1, f2 } = ellipseLC.rightAngled();
        const p1 = ellipseLC.center.plus(f1), p2 = ellipseLC.center.plus(f2);
        // check if both endpoints are on the cone's surface
        // and that one main axis is perpendicular to the Z-axis
        return (eq(Math.pow(p1.x, 2) + Math.pow(p1.y, 2), Math.pow(p1.z, 2)) &&
            eq(Math.pow(p2.x, 2) + Math.pow(p2.y, 2), Math.pow(p2.z, 2)) &&
            (eq0(f1.z) || eq0(f2.z)));
    }
    containsLine(line) {
        const lineLC = line.transform(this.matrixInverse);
        const d = lineLC.dir1;
        return lineLC.containsPoint(V3.O) && eq(d.x * d.x + d.y * d.y, d.z * d.z);
    }
    containsParabola(curve) {
        assertInst(ParabolaCurve$1, curve);
        const curveLC = curve.transform(this.matrixInverse);
        if (curveLC.center.z < 0 || curveLC.f2.z < 0) {
            return false;
        }
        const { center, f1, f2 } = curveLC.rightAngled();
        // check if center is on the surface,
        // that tangent is perpendicular to the Z-axis
        // and that "y" axis is parallel to surface
        return (eq(center.x * center.x + center.y * center.y, center.z * center.z) &&
            eq0(f1.z) &&
            eq(f2.x * f2.x + f2.y * f2.y, f2.z * f2.z));
    }
    containsHyperbola(curve) {
        // calculate intersection of plane ax + cz = 1 and cone x² + y² = z²
        // const center = new V3(-a / (cc - aa), 0, 1 / (cc - aa))
        // const p1 = new V3(1 / (a - c), 0, -1 / (a - c))
        // const p2 = new V3(-a / (cc - aa), 1 / sqrt(aa - cc), 1 / (cc - aa))
        // const f1 = new V3(1 * c / (aa - cc), 0, -a / (aa - cc) )
        // const f2 = new V3(0, 1 / sqrt(aa - cc), 0)
        assertInst(HyperbolaCurve$1, curve);
        const curveLC = curve.transform(this.matrixInverse).rightAngled();
        const centerXY = curveLC.center.xy();
        if (centerXY.likeO()) {
            return false;
        }
        const rot = centerXY.angleXY();
        const { center, f1, f2 } = curveLC.rotateZ(-rot);
        // s = a / (aa - cc)
        // t = -c / (aa - cc)
        // s + t = 1 / (a + c)
        // s - t = 1 / (a - c)
        // (s + t)(s - t) = (ss - tt) = 1 / (aa - cc)
        // u = 1 / sqrt(aa - cc) = sqrt(ss - tt)
        // check if center is on the surface,
        // that tangent is perpendicular to the Z-axis
        return (f1.z > 0 &&
            eq(center.x, f1.z) &&
            eq(center.z, f1.x) &&
            eq0(center.y) &&
            eq0(f1.y) &&
            eq(sqrt(abs(Math.pow(center.x, 2) - Math.pow(center.z, 2))), abs(f2.y)) &&
            eq0(f2.x) &&
            eq0(f2.z));
    }
    containsCurve(curve) {
        if (curve instanceof EllipseCurve$1) {
            return this.containsEllipse(curve);
        }
        else if (curve instanceof L3$1) {
            return this.containsLine(curve);
        }
        else if (curve instanceof HyperbolaCurve$1) {
            return this.containsHyperbola(curve);
        }
        else if (curve instanceof ParabolaCurve$1) {
            return this.containsParabola(curve);
        }
        else {
            return super.containsCurve(curve);
        }
    }
    transform(m4) {
        return new ConicSurface(m4.transformPoint(this.center), m4.transformVector(this.f1).times(m4.isMirroring() ? -1 : 1), m4.transformVector(this.f2), m4.transformVector(this.dir), this.uMin, this.uMax, this.vMin, this.vMax);
    }
    transform4(m4) {
        const transformedApex = m4.timesVector(Vector.fromV3AndWeight(this.center, 1));
        const isometricZ = (z) => new EllipseCurve$1(new V3(0, 0, z), new V3(z, 0, 0), new V3(0, z, 0));
        if (!eq0(transformedApex.w)) {
            // sMin doesn't change, but tMin does...
            const c = m4.transformPoint(this.center), f1 = m4
                .transformVector2(this.f1, this.center)
                .times(m4.isMirroring() ? -1 : 1), f2 = m4.transformVector2(this.f2, this.center), dir = m4.transformVector2(this.dir, this.center);
            const matrixInv = M4.forSys(f1, f2, dir, c).inversed();
            const aabb = isometricZ(this.vMin)
                .transform4(matrixInv.times(m4.times(this.matrix)))
                .getAABB()
                .addAABB(isometricZ(this.vMax)
                .transform4(matrixInv.times(m4.times(this.matrix)))
                .getAABB());
            return new ConicSurface(c, f1, f2, dir, this.uMin, this.uMax, aabb.min.z, aabb.max.z);
        }
        else {
            const dir = transformedApex.V3();
            const baseCurve = isometricZ(this.vMin).transform4(m4.times(this.matrix));
            const matrixInv = M4.forSys(baseCurve.f1, baseCurve.f2, dir.unit(), baseCurve.center).inversed();
            const aabb = isometricZ(this.vMax)
                .transform4(matrixInv.times(m4.times(this.matrix)))
                .getAABB();
            return new CylinderSurface$1(baseCurve, dir.unit(), this.uMin, this.uMax, min(0, aabb.min.z, aabb.max.z), max(0, aabb.min.z, aabb.max.z));
        }
    }
    flipped() {
        return new ConicSurface(this.center, this.f1.negated(), this.f2, this.dir);
    }
    normalUVFunc() {
        const { f1, f2 } = this, f3 = this.dir;
        return (d, _z) => {
            return f2
                .cross(f1)
                .plus(f2.cross(f3.times(Math.cos(d))))
                .plus(f3.cross(f1.times(Math.sin(d))))
                .unit();
        };
    }
    normalP(p) {
        //TODO assert(!p.like(this.center))
        const pLC = this.matrixInverse.transformPoint(p);
        return this.normalUVFunc()(pLC.angleXY(), pLC.z);
    }
    pUVFunc() {
        return (u, v) => {
            // center + f1 v cos u + f2 v sin u + v dir
            const resultLC = new V3(v * cos(u), v * sin(u), v);
            return this.matrix.transformPoint(resultLC);
        };
    }
    dpdu() {
        return (u, v) => {
            const resultLC = new V3(v * -sin(u), v * cos(u), 0);
            return this.matrix.transformVector(resultLC);
        };
    }
    dpdv() {
        return (s) => {
            const resultLC = new V3(cos(s), sin(s), 1);
            return this.matrix.transformVector(resultLC);
        };
    }
    implicitFunction() {
        return (pWC) => {
            const pLC = this.matrixInverse.transformPoint(pWC);
            const radiusLC = pLC.lengthXY();
            return this.normalDir * (radiusLC - pLC.z);
        };
    }
    didp(pWC) {
        const pLC = this.matrixInverse.transformPoint(pWC);
        return this.pLCNormalWCMatrix.transformVector(pLC.xy().unit().withElement("z", -1).times(this.normalDir));
    }
    containsPoint(p) {
        return eq0(this.implicitFunction()(p));
    }
    uvP(pWC) {
        const pLC = this.matrixInverse.transformPoint(pWC);
        const angle = pLC.angleXY();
        return new V3(angle < -PI / 2 ? angle + TAU : angle, pLC.z, 0);
    }
    isCurvesWithSurface(surface) {
        if (surface instanceof PlaneSurface$1) {
            return this.isCurvesWithPlane(surface.plane);
        }
        else if (ImplicitSurface$1.is(surface)) {
            return ParametricSurface$1.isCurvesParametricImplicitSurface(this, surface, 0.1, 0.1 / this.dir.length(), 0.02);
        }
        return super.isCurvesWithSurface(surface);
    }
    getCenterLine() {
        return new L3$1(this.center, this.dir);
    }
    isCurvesWithPlane(plane) {
        assertInst(P3$1, plane);
        const planeLC = plane.transform(this.matrixInverse);
        const planeNormal = planeLC.normal1;
        const c = planeNormal.z;
        /** "rotate" plane normal1 when passing to {@link ConicSurface.unitISPlane} so that
         *  y-component of normal1 is 0 */
        const a = planeNormal.lengthXY();
        const d = planeLC.w;
        // generated curves need to be rotated back before transforming to world coordinates
        const rotationMatrix = M4.rotateZ(planeNormal.angleXY());
        const wcMatrix = eq0(planeNormal.lengthXY())
            ? this.matrix
            : this.matrix.times(rotationMatrix);
        return ConicSurface.unitISPlane(a, c, d).flatMap((curve) => {
            const curveWC = curve.transform(wcMatrix);
            if (curve instanceof EllipseCurve$1) {
                const curveLC = curve.transform(rotationMatrix);
                const ts = curveLC.isTsWithPlane(P3$1.ZX);
                const intervals = getIntervals(ts, -PI, PI).filter(([a, b]) => curveLC.at((a + b) / 2).y > 0);
                return intervals.flatMap(([a, b]) => curveWC.split(a, b));
            }
            const p = curveWC.at(0.2);
            return this.normalP(p).cross(plane.normal1).dot(curveWC.tangentAt(0.2)) >
                0
                ? curveWC
                : curveWC.reversed();
        });
    }
    debugInfo() {
        return {
            ps: [this.center],
            lines: [
                this.center,
                this.center.plus(this.f1),
                this.center.plus(this.f2),
                this.center.plus(this.dir),
            ],
        };
    }
}
/**
 * Unit cone. x² + y² = z², 0 <= z
 */
ConicSurface.UNIT = new ConicSurface(V3.O, V3.X, V3.Y, V3.Z);
ConicSurface.prototype.uStep = PI / 16;
ConicSurface.prototype.vStep = 256;

// }
// [].bar()
/**
 * Surface normal1 is (t, z) => this.baseCurve.tangentAt(t) X this.dir
 * Choose dir appropriately to select surface orientation.
 */
class ProjectedCurveSurface extends ParametricSurface {
    constructor(baseCurve, dir, uMin = baseCurve.tMin, uMax = baseCurve.tMax, vMin = -100, vMax = 100) {
        super(uMin, uMax, vMin, vMax);
        this.baseCurve = baseCurve;
        this.dir = dir;
        assertInst(Curve, baseCurve);
        assertInst(V3, dir);
        assert(uMin < uMax);
        assert(vMin < vMax);
    }
    getConstructorParameters() {
        return [
            this.baseCurve,
            this.dir,
            this.uMin,
            this.uMax,
            this.vMin,
            this.vMax,
        ];
    }
    equals(obj) {
        return (this == obj ||
            (Object.getPrototypeOf(this) == Object.getPrototypeOf(obj) &&
                this.dir.equals(obj.dir) &&
                this.baseCurve.equals(obj.baseCurve)));
    }
    hashCode() {
        return [this.dir, this.baseCurve].hashCode();
    }
    containsLine(line) {
        return this.dir.isParallelTo(line.dir1) && this.containsPoint(line.anchor);
    }
    dpdu() {
        return (u, v) => this.baseCurve.tangentAt(u);
    }
    dpdv() {
        return (u, v) => this.dir;
    }
    normalUV(u, v) {
        return this.baseCurve.tangentAt(u).cross(this.dir).unit();
    }
    pUV(u, v) {
        return this.baseCurve.at(u).plus(this.dir.times(v));
    }
    pointFoot(pWC, ss) {
        const basePlane = new P3(this.dir.unit(), 0);
        const projCurve = this.baseCurve.project(basePlane);
        const projPoint = basePlane.projectedPoint(pWC);
        const t = projCurve.closestTToPoint(projPoint, ss, this.uMin, this.uMax);
        const z = L3.pointT(this.baseCurve.at(t), this.dir, pWC);
        return new V3(t, z, 0);
    }
    uvPFunc() {
        const projPlane = new P3(this.dir.unit(), 0);
        const projBaseCurve = this.baseCurve.project(projPlane);
        return (pWC) => {
            const projPoint = projPlane.projectedPoint(pWC);
            assertNumbers(this.uMin);
            const t = projBaseCurve.pointT(projPoint, this.uMin, this.uMax);
            const z = L3.pointT(this.baseCurve.at(t), this.dir, pWC);
            return new V3(t, z, 0);
        };
    }
    isCurvesWithPlane(plane) {
        assertInst(P3, plane);
        if (this.dir.isPerpendicularTo(plane.normal1)) {
            const ts = this.baseCurve.isTsWithPlane(plane);
            return ts.map((t) => {
                const l3dir = 0 < this.baseCurve.tangentAt(t).dot(plane.normal1)
                    ? this.dir
                    : this.dir.negated();
                return new L3(this.baseCurve.at(t), l3dir.unit());
            });
        }
        else {
            let projCurve = this.baseCurve.transform(M4.project(plane, this.dir));
            if (this.dir.dot(plane.normal1) > 0) {
                // we need to flip the ellipse so the tangent is correct
                projCurve = projCurve.reversed();
            }
            return [projCurve];
        }
    }
    isCurvesWithSurface(surface) {
        if (surface instanceof PlaneSurface) {
            return this.isCurvesWithPlane(surface.plane);
        }
        if (surface instanceof ProjectedCurveSurface) {
            const dir1 = surface.dir;
            if (this.dir.isParallelTo(dir1)) {
                const ts = surface.baseCurve.isTsWithSurface(this);
                return ts.map((t) => {
                    const p = surface.baseCurve.at(t);
                    const correctDir = this.normalP(p).cross(surface.normalP(p));
                    return new L3(p, dir1.times(sign(correctDir.dot(dir1))));
                });
            }
            else if (ImplicitSurface.is(surface)) {
                let curves2 = ParametricSurface.isCurvesParametricImplicitSurface(this, surface, 0.1, 0.1 / surface.dir.length(), 0.05);
                curves2 = surface.clipCurves(curves2);
                return curves2;
            }
            else {
                let curves2 = ParametricSurface.isCurvesParametricParametricSurface(this, surface, 0.05, 0.1 / surface.dir.length(), 0.05);
                curves2 = this.clipCurves(curves2);
                curves2 = surface.clipCurves(curves2);
                return curves2;
            }
        }
        if (surface instanceof EllipsoidSurface) {
            return surface.isCurvesWithSurface(this);
        }
        return super.isCurvesWithSurface(surface);
    }
    containsPoint(pWC) {
        const uv = this.uvPFunc()(pWC);
        return this.pUVFunc()(uv.x, uv.y).like(pWC);
    }
    containsCurve(curve) {
        if (curve instanceof L3) {
            return (this.dir.isParallelTo(curve.dir1) && this.containsPoint(curve.anchor));
        }
        if (curve instanceof ImplicitCurve) {
            return super.containsCurve(curve);
        }
        // project baseCurve and test curve onto a common plane and check if the curves are alike
        const projPlane = new P3(this.dir.unit(), 0);
        const projBaseCurve = this.baseCurve.project(projPlane);
        const projCurve = curve.project(projPlane);
        return projBaseCurve.isColinearTo(projCurve);
    }
    isCoplanarTo(surface) {
        return (this == surface ||
            (hasConstructor(surface, ProjectedCurveSurface) &&
                this.dir.isParallelTo(surface.dir) &&
                this.containsCurve(surface.baseCurve)));
    }
    like(object) {
        if (!this.isCoplanarTo(object))
            return false;
        // normals need to point in the same direction (outwards or inwards) for both
        const p00 = this.pUVFunc()(0, 0);
        const thisNormal = this.normalUVFunc()(0, 0);
        const otherNormal = object.normalP(p00);
        return 0 < thisNormal.dot(otherNormal);
    }
    loopContainsPoint(loop, p) {
        assertVectors(p);
        assert(isFinite(p.x), p.y, p.z);
        const line = new L3(p, this.dir.unit());
        const ptpf = this.uvPFunc();
        const pp = ptpf(p);
        if (isNaN(pp.x)) {
            console.log(this.sce, p.sce);
            assert(false);
        }
        const lineOut = this.baseCurve.tangentAt(pp.x).rejectedFrom(this.dir);
        return Surface.loopContainsPointGeneral(loop, p, line, lineOut);
    }
    transform(m4) {
        const f = m4.isMirroring() ? -1 : 1;
        return new this.constructor(this.baseCurve.transform(m4), m4.transformVector(this.dir).times(f), this.uMin, this.uMax, 1 == f ? this.vMin : -this.vMax, 1 == f ? this.vMax : -this.vMin);
    }
    transform4(m4) {
        const vp = m4.vanishingPoint(this.dir);
        if (!vp) {
            const f = m4.isMirroring() ? -1 : 1;
            return new this.constructor(this.baseCurve.transform4(m4), m4.normalized().transformVector(this.dir).times(f), undefined, undefined, 1 == f ? this.tMin : -this.tMax, 1 == f ? this.tMax : -this.tMin);
        }
        const curveT = this.baseCurve.transform4(m4);
        if (curveT instanceof EllipseCurve) {
            console.log(vp.sce, curveT.sce);
            return ConicSurface.atApexThroughEllipse(vp, m4.isMirroring() ? curveT : curveT.reversed(), this.sMin, this.sMax, 1, 2);
        }
        return new PointProjectedSurface(curveT, vp, P3.throughPoints(curveT.at(curveT.tMin), curveT.at((curveT.tMin + curveT.tMax) / 2), curveT.at(curveT.tMax)), 1, this.sMin, this.sMax, 1, 2);
    }
    isTsForLine(line) {
        assertInst(L3, line);
        const projPlane = new P3(this.dir.unit(), 0);
        const projDir = projPlane.projectedVector(line.dir1);
        if (projDir.likeO()) {
            // line is parallel to this.dir
            return [];
        }
        const projAnchor = projPlane.projectedPoint(line.anchor);
        const projBaseCurve = this.baseCurve.project(projPlane);
        return projBaseCurve
            .isInfosWithLine(projAnchor, projDir, this.uMin, this.uMax, line.tMin, line.tMax)
            .map((info) => info.tOther);
    }
    flipped() {
        return new this.constructor(this.baseCurve, this.dir.negated(), this.uMin, this.uMax, -this.vMax, -this.vMin);
    }
}
ProjectedCurveSurface.prototype.uStep = 1 / 128;
ProjectedCurveSurface.prototype.vStep = 256;

/**
 * Rotation surface with r = f(z)
 */
class RotatedCurveSurface extends ParametricSurface$1 {
    constructor(curve, matrix = M4.IDENTITY, uMin = 0, uMax = PI, vMin = curve.tMin, vMax = curve.tMax) {
        // d/dz (r(z))
        super(uMin, uMax, vMin, vMax);
        this.curve = curve;
        this.matrix = matrix;
        assertInst(M4, matrix);
        assert(matrix.isNoProj());
        assert(eq0(curve.at(vMin).y));
        this.matrixInverse = matrix.inversed();
        this.vStep = this.curve.tIncrement;
    }
    getConstructorParameters() {
        return [this.curve, this.matrix, this.uMin, this.uMax, this.vMin, this.vMax];
    }
    flipped() {
        return new RotatedCurveSurface(this.curve, this.matrix.times(M4.mirror(P3$1.YZ)), this.uMin, this.uMax, this.vMin, this.vMax);
    }
    transform(m4) {
        return new RotatedCurveSurface(this.curve, m4.isMirroring()
            ? m4.times(this.matrix).times(M4.mirror(P3$1.YZ))
            : m4.times(this.matrix), this.uMin, this.uMax, this.vMin, this.vMax);
    }
    containsPoint(pWC) {
        const pLC = this.matrixInverse.transformPoint(pWC);
        const radius = pLC.lengthXY();
        return this.curve.containsPoint(new V3(radius, 0, pLC.z));
    }
    pUVFunc() {
        return (u, v) => {
            const { x: radius, z: z } = this.curve.at(v);
            return this.matrix.transformPoint(V3.polar(radius, u, z));
        };
    }
    dpdu() {
        return (u, v) => {
            const radius = this.curve.at(v).x;
            const resultLC = new V3(radius * -sin(u), radius * cos(u), 0);
            return this.matrix.transformVector(resultLC);
        };
    }
    dpdv() {
        return (u, v) => {
            const { x: drdt, z: dzdt } = this.curve.tangentAt(v);
            return this.matrix.transformVector(V3.polar(drdt, u, dzdt));
        };
    }
    normalUVFunc() {
        const matrix = this.matrix.inversed().transposed().as3x3();
        const normalLength = this.matrix.isMirroring() ? -1 : 1;
        return (u, v) => {
            const { x: drdt, z: dzdt } = this.curve.tangentAt(v);
            return matrix
                .transformVector(V3.polar(dzdt, u, -drdt))
                .toLength(normalLength);
        };
    }
    uvPFunc() {
        return (pWC) => {
            const pLC = this.matrixInverse.transformPoint(pWC);
            const angle = EllipseCurve$1.XYLCPointT(pLC, this.uMin, this.uMax);
            const radius = pLC.lengthXY();
            return new V3(angle, this.curve.pointT(new V3(radius, 0, pLC.z)), 0);
        };
    }
    pointFoot(pWC, startS, startT) {
        const pLC = this.matrixInverse.transformPoint(pWC);
        const angle = abs(pLC.angleXY());
        const radius = pLC.lengthXY();
        return new V3(angle, this.curve.closestTToPoint(new V3(radius, 0, pLC.z)), 0);
    }
    isTsForLine(line) {
        const anchorLC = this.matrixInverse.transformPoint(line.anchor);
        const dirLC = this.matrixInverse.transformVector(line.dir1);
        if (dirLC.isParallelTo(V3.Z)) {
            if (!fuzzyBetween(anchorLC.angleXY(), this.uMin, this.uMax))
                return [];
            return this.curve
                .isInfosWithLine(new V3(anchorLC.lengthXY(), 0, anchorLC.z), dirLC)
                .map((info) => info.tOther);
        }
        else if (L3$1.containsPoint(anchorLC.xy(), dirLC.xy(), V3.O)) {
            // line goes through Z axis
            const dotter = dirLC.xy().unit();
            return [
                ...this.curve.isInfosWithLine(new V3(dotter.dot(anchorLC), 0, anchorLC.z), new V3(dotter.dot(dirLC), 0, dirLC.z)),
                ...this.curve.isInfosWithLine(new V3(-dotter.dot(anchorLC), 0, anchorLC.z), new V3(-dotter.dot(dirLC), 0, dirLC.z)),
            ]
                .map((info) => info.tOther)
                .filter((t) => fuzzyBetween(L3$1.at(anchorLC, dirLC, t).angleXY(), this.uMin, this.uMax));
        }
        else if (dirLC.isPerpendicularTo(V3.Z)) {
            const secs = this.isCurvesWithPlaneLC(new P3$1(V3.Z, anchorLC.z));
            if (!secs)
                return [];
            return secs.flatMap((sec) => sec.isInfosWithLine(anchorLC, dirLC).map((info) => info.tOther));
        }
        else {
            // transform into hyperbola
            // f(t) = V(((ax + t dx)² + (ay + t dy)²) ** 1/2, 0, az + t dz)
            // f(t) = V((ax² + 2 ax t dx + t² dx² + ay² + 2 ay t dy + t² dy²) ** 1/2, 0, az + t dz)
            // f(t) = V((t² (dx² + dy²) + 2 t (ax dx + ay dy) + ax² + ay²) ** 1/2, 0, az + t * dz)
            // (anchorLC.xy + t * dirLC.xy) * dir.xy = 0
            // t * dirLC.xy² = -anchorLC.xy * dirLC.xy
            const closestTToZ = -anchorLC.xy().dot(dirLC.xy()) / dirLC.xy().squared();
            const closestPointToZ = L3$1.at(anchorLC, dirLC, closestTToZ);
            const scaleX = closestPointToZ.lengthXY();
            const lineGradientWC = dirLC.z / dirLC.lengthXY();
            const scaleZ = scaleX * lineGradientWC;
            const hc = HyperbolaCurve$1.XY.transform(M4.rotateX(90 * DEG)
                .scale(scaleX, 0, scaleZ)
                .translate(0, 0, closestPointToZ.z));
            const infos = hc.isInfosWithCurve(this.curve);
            return infos
                .map((info) => (info.p.z - anchorLC.z) / dirLC.z)
                .filter((t) => fuzzyBetween(L3$1.at(anchorLC, dirLC, t).angleXY(), this.uMin, this.uMax));
        }
    }
    isCurvesWithPlaneLC(planeLC) {
        if (planeLC.normal1.isParallelTo(V3.Z)) {
            return this.curve.isTsWithPlane(planeLC).map((t) => {
                const { x: radius } = this.curve.at(t);
                return new EllipseCurve$1(new V3(0, 0, planeLC.w), new V3(radius, 0, 0), new V3(0, radius, 0), this.uMin, this.uMax).transform(this.matrix);
            });
        }
        else if (planeLC.normal1.isPerpendicularTo(V3.Z) &&
            planeLC.containsPoint(V3.O)) {
            return [
                this.curve
                    .rotateZ(V3.Y.angleRelativeNormal(planeLC.normal1, V3.Z))
                    .transform(this.matrix),
            ];
        }
        return undefined;
    }
    isCurvesWithPlane(plane) {
        const planeLC = plane.transform(this.matrixInverse);
        const planeLCCurves = this.isCurvesWithPlaneLC(planeLC);
        if (planeLCCurves) {
            return planeLCCurves.map((curve) => curve.transform(this.matrix));
        }
        else {
            return ParametricSurface$1.isCurvesParametricImplicitSurface(this, new PlaneSurface$1(plane), 0.05, 0.05, 0.02);
        }
    }
    loopContainsPoint(loop, pWC) {
        const pLC = this.matrixInverse.transformPoint(pWC);
        const angle = EllipseCurve$1.XYLCPointT(pLC, this.uMin, this.uMax);
        const testCurveLC = EllipseCurve$1.semicircle(pLC.lengthXY(), new V3(0, 0, pLC.z));
        const testCurveWC = testCurveLC.transform(this.matrix);
        return Surface$1.loopContainsPointEllipse(loop, pWC, testCurveWC, angle);
    }
    isCoplanarTo(surface) {
        if (this === surface)
            return true;
        if (!hasConstructor(surface, RotatedCurveSurface))
            return false;
        const surfaceLCToThisLC = this.matrixInverse.times(surface.matrix);
        assert(!surfaceLCToThisLC.X.xy().likeO());
        const zRotation = surfaceLCToThisLC.X.angleXY();
        return surface.curve
            .transform(M4.rotateZ(-zRotation).times(surfaceLCToThisLC))
            .isColinearTo(this.curve);
    }
    isCurvesWithSurface(surface) {
        if (surface instanceof PlaneSurface$1) {
            return this.isCurvesWithPlane(surface.plane);
        }
        return super.isCurvesWithSurface(surface);
    }
    containsCurve(curve) {
        if (curve.constructor == this.curve.constructor) {
            const curveLC = curve.transform(this.matrixInverse);
            // find a point on curveLC which isn't on the Z-axis
            const t = [0, 0.5, 1]
                .map((x) => lerp(curveLC.tMin, curveLC.tMax, x))
                .withMax((t) => curveLC.at(t).lengthXY());
            const angle = curveLC.at(t).angleXY();
            const curveLCRotated = curveLC.rotateZ(-angle);
            if (this.curve.isColinearTo(curveLCRotated)) {
                return true;
            }
        }
        if (curve instanceof EllipseCurve$1) {
            const curveLC = curve.transform(this.matrixInverse);
            if (curveLC.normal.isParallelTo(V3.Z)) {
                return (curveLC.isCircular() &&
                    this.curve.containsPoint(new V3(curveLC.f1.length(), 0, curveLC.center.z)));
            }
            return false;
        }
        return super.containsCurve(curve);
    }
    getExtremePoints() {
        return getExtremePointsHelper.call(this, this.curve);
    }
    asNURBSSurface() {
        // y = 0 for baseNURBS
        const baseNURBS = NURBS$1.fromEllipse(this.curve);
        const rotationNURBS = NURBS$1.UnitCircle(2, this.tMin, this.tMax);
        return new NURBSSurface$1(rotationNURBS.points.flatMap((rv) => baseNURBS.points.map((b) => this.matrix.timesVector(VV(rv.x * b.x, rv.y * b.x, b.z * rv.w, rv.w * b.w)))), baseNURBS.knots, rotationNURBS.knots, baseNURBS.degree, rotationNURBS.degree, baseNURBS.tMin, baseNURBS.tMax, rotationNURBS.tMin, rotationNURBS.tMax);
    }
}
RotatedCurveSurface.prototype.uStep = EllipseCurve$1.prototype.tIncrement;
function getExtremePointsHelper(curve) {
    // this logic comes from EllipseCurve.roots
    const f1 = this.matrix.X;
    const f2 = this.matrix.Y;
    return [0, 1, 2].flatMap((dim) => {
        const a = f2.e(dim), b = -f1.e(dim);
        const xiEtas = eq0(a) && eq0(b) ? [[1, 0]] : intersectionUnitCircleLine2$1(a, b, 0);
        return xiEtas.flatMap(([xi, eta]) => {
            const u = Math.atan2(eta, xi);
            if (!(lt(this.uMin, u) && lt(u, this.uMax)))
                return [];
            const testCurve = curve.transform(this.matrix.times(M4.rotateZ(u)));
            return testCurve.roots()[dim].map((v) => this.pUV(u, v));
        });
    });
}

class CylinderSurface extends ProjectedCurveSurface$1 {
    // @ts-ignore
    // readonly baseCurve: EllipseCurve
    constructor(baseCurve, dir1, uMin = baseCurve.tMin, uMax = baseCurve.tMax, zMin = -Infinity, zMax = Infinity) {
        super(baseCurve, dir1, uMin, uMax, zMin, zMax);
        this.baseCurve = baseCurve;
        assertInst(EllipseCurve$1, baseCurve);
        //assert(!baseCurve.normal1.isPerpendicularTo(dir1), !baseCurve.normal1.isPerpendicularTo(dir1))
        this.matrix = M4.forSys(baseCurve.f1, baseCurve.f2, dir1, baseCurve.center);
        this.matrixInverse = this.matrix.inversed();
        this.normalDir = sign(this.baseCurve.normal.dot(this.dir));
        this.pLCNormalWCMatrix = this.matrix
            .as3x3()
            .inversed()
            .transposed()
            .scale(this.normalDir);
        this.pWCNormalWCMatrix = this.pLCNormalWCMatrix.times(this.matrixInverse);
    }
    static semicylinder(radius, sMin, sMax, tMin, tMax) {
        return new CylinderSurface(new EllipseCurve$1(V3.O, new V3(radius, 0, 0), new V3(0, radius, 0)), V3.Z, sMin, sMax, tMin, tMax);
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
        return pqFormula(b / a, c / a).filter((t) => EllipseCurve$1.XYLCValid(new V3(ax + dx * t, ay + dy * t, 0)));
    }
    normalP(p) {
        return this.pLCNormalWCMatrix
            .transformVector(this.matrixInverse.transformPoint(p).xy())
            .unit();
    }
    loopContainsPoint(loop, p) {
        assertVectors(p);
        if (!this.containsPoint(p))
            return OUTSIDE$1;
        const line = new L3$1(p, this.dir.unit());
        const lineOut = this.dir.cross(this.normalP(p));
        return Surface$1.loopContainsPointGeneral(loop, p, line, lineOut);
    }
    isTsForLine(line) {
        assertInst(L3$1, line);
        // transforming line manually has advantage that dir1 will not be renormalized,
        // meaning that calculated values t for localLine are directly transferable to line
        const dirLC = this.matrixInverse.transformVector(line.dir1);
        if (dirLC.isParallelTo(V3.Z)) {
            // line is parallel to this.dir
            return [];
        }
        const anchorLC = this.matrixInverse.transformPoint(line.anchor);
        assert(!CylinderSurface.unitISLineTs(anchorLC, dirLC).length ||
            !isNaN(CylinderSurface.unitISLineTs(anchorLC, dirLC)[0]), "sad " + dirLC);
        return CylinderSurface.unitISLineTs(anchorLC, dirLC);
    }
    isCoplanarTo(surface) {
        return (this == surface ||
            (hasConstructor(surface, CylinderSurface) &&
                this.dir.isParallelTo(surface.dir) &&
                this.containsEllipse(surface.baseCurve, false)));
    }
    like(surface) {
        if (!this.isCoplanarTo(surface))
            return false;
        // normals need to point in the same direction (outwards or inwards) for both
        const thisFacesOut = 0 < this.baseCurve.normal.dot(this.dir);
        const objectFacesOut = 0 < surface.baseCurve.normal.dot(surface.dir);
        return thisFacesOut == objectFacesOut;
    }
    containsEllipse(ellipse, checkAABB = true) {
        const projEllipse = ellipse.transform(M4.project(this.baseCurve.getPlane(), this.dir));
        return this.baseCurve == ellipse || this.baseCurve.isColinearTo(projEllipse);
    }
    containsCurve(curve) {
        if (curve instanceof L3$1) {
            return this.containsLine(curve);
        }
        else if (curve instanceof EllipseCurve$1) {
            return this.containsEllipse(curve);
        }
        else if (curve instanceof BezierCurve$1) {
            return false;
        }
        else {
            return super.containsCurve(curve);
        }
    }
    implicitFunction() {
        return (pWC) => {
            const pLC = this.matrixInverse.transformPoint(pWC);
            return (pLC.lengthXY() - 1) * this.normalDir;
        };
    }
    didp(pWC) {
        const pLC = this.matrixInverse.transformPoint(pWC);
        const pLCLengthXY = pLC.lengthXY();
        const didpLC = new V3(pLC.x / pLCLengthXY, pLC.y / pLCLengthXY, 0);
        return this.pLCNormalWCMatrix.transformVector(didpLC);
    }
    containsPoint(pWC) {
        const pLC = this.matrixInverse.transformPoint(pWC);
        return this.baseCurve.isValidT(EllipseCurve$1.XYLCPointT(pLC, this.uMin, this.uMax));
    }
    uvP(pWC) {
        assert(arguments.length == 1);
        const pLC = this.matrixInverse.transformPoint(pWC);
        const u = EllipseCurve$1.XYLCPointT(pLC, this.vMin, this.vMax);
        return new V3(u, pLC.z, 0);
    }
    isCurvesWithSurface(surface2) {
        if (surface2 instanceof ProjectedCurveSurface$1) {
            if (surface2.dir.isParallelTo(this.dir)) {
                const projectedCurve = surface2.baseCurve.transform(M4.project(this.baseCurve.getPlane(), this.dir));
                return this.baseCurve.isInfosWithCurve(projectedCurve).map((info) => {
                    const lineDir = sign(this.normalP(info.p)
                        .cross(surface2.normalP(info.p))
                        .dot(this.dir)) || 1;
                    return new L3$1(info.p, this.dir.times(lineDir));
                });
            }
        }
        if (surface2 instanceof CylinderSurface) {
            if (eq0(this.getCenterLine().distanceToLine(surface2.getCenterLine()))) {
                throw new Error();
            }
        }
        return super.isCurvesWithSurface(surface2);
    }
    getCenterLine() {
        return new L3$1(this.baseCurve.center, this.dir);
    }
    facesOutwards() {
        return this.baseCurve.normal.dot(this.dir) > 0;
    }
    getSeamPlane() {
        let normal = this.baseCurve.f1.cross(this.dir);
        normal = normal.times(-sign(normal.dot(this.baseCurve.f2)));
        return P3$1.normalOnAnchor(normal, this.baseCurve.center);
    }
    clipCurves(curves) {
        return curves.flatMap((curve) => curve.clipPlane(this.getSeamPlane()));
    }
}
CylinderSurface.UNIT = new CylinderSurface(EllipseCurve$1.UNIT, V3.Z, undefined, undefined, 0, 1);
CylinderSurface.prototype.uStep = TAU / 32;
CylinderSurface.prototype.vStep = 256;

class EllipsoidSurface extends ParametricSurface$1 {
    constructor(center, f1, f2, f3, uMin = 0, uMax = PI, vMin = -PI / 2, vMax = PI / 2) {
        super(uMin, uMax, vMin, vMax);
        this.center = center;
        this.f1 = f1;
        this.f2 = f2;
        this.f3 = f3;
        assert(0 <= uMin && uMin <= PI);
        assert(0 <= uMax && uMax <= PI);
        assert(-PI / 2 <= vMin && vMin <= PI / 2);
        assert(-PI / 2 <= vMax && vMax <= PI / 2);
        assertVectors(center, f1, f2, f3);
        this.matrix = M4.forSys(f1, f2, f3, center);
        this.matrixInverse = this.matrix.inversed();
        this.normalDir = sign(this.f1.cross(this.f2).dot(this.f3));
        this.pLCNormalWCMatrix = this.matrix
            .as3x3()
            .inversed()
            .transposed()
            .scale(this.normalDir);
        this.pWCNormalWCMatrix = this.pLCNormalWCMatrix.times(this.matrixInverse);
    }
    static unitArea(contour) {
        const totalArea = contour
            .map((edge) => {
            if (edge.curve instanceof PICurve$1) {
                const points = edge.curve.calcSegmentPoints(edge.aT, edge.bT, edge.a, edge.b, edge.aT > edge.bT, true);
                let sum = 0;
                for (let i = 0; i < points.length - 1; i++) {
                    const p = points[i], ppp = points[i + 1];
                    sum += ((abs(p.angleXY()) + abs(ppp.angleXY())) / 2) * (ppp.z - p.z);
                }
                return sum;
            }
            else if (edge.curve instanceof EllipseCurve$1) {
                const f = (t) => {
                    const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t);
                    const angleXY = abs(at.angleXY());
                    //const arcLength = angleXY * Math.sqrt(1 - at.z ** 2) ( == at.lengthXY())
                    //const scaling = tangent.z / at.lengthXY()
                    return angleXY * tangent.z;
                };
                const val = glqInSteps(f, edge.aT, edge.bT, 1);
                return val;
            }
            else {
                throw new Error();
            }
        })
            .sum();
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
        return pqFormula(b / a, c / a).filter((t) => le(0, anchor.y + t * dir.y));
    }
    /**
     * unit sphere: x² + y² + z² = 1
     * plane: normal1 DOT p = w
     */
    static unitISCurvesWithPlane(plane) {
        const distPlaneCenter = Math.abs(plane.w);
        if (lt(distPlaneCenter, 1)) {
            // result is a circle
            // radius of circle: imagine right angled triangle (origin -> center of intersection circle -> point on
            // intersection circle) pythagoras: 1² == distPlaneCenter² + isCircleRadius² => isCircleRadius == sqrt(1 -
            // distPlaneCenter²)
            const isCircleRadius = Math.sqrt(1 - Math.pow(distPlaneCenter, 2));
            const anchorY = plane.normal1.y * plane.w;
            const d = abs(distPlaneCenter * isCircleRadius);
            if (le(anchorY, -d) && !eq0(distPlaneCenter)) {
                return [];
            }
            else if (le(anchorY, 0) && !plane.normal1.isParallelTo(V3.Y)) {
                const f1 = plane.normal1.isParallelTo(V3.Y)
                    ? V3.Z
                    : plane.normal1.cross(V3.Y).toLength(isCircleRadius);
                const f2 = f1.cross(plane.normal1);
                const minEta = -anchorY / f2.y, minT = max(0, Math.asin(minEta));
                return [new EllipseCurve$1(plane.anchor, f1, f2, minT, PI - minT)];
            }
            else {
                const f2 = (plane.normal1.isParallelTo(V3.Y)
                    ? V3.X
                    : plane.normal1.cross(V3.Y)).toLength(isCircleRadius);
                const f1 = f2.cross(plane.normal1);
                const minXi = eq0(f1.y) ? -1 : -anchorY / f1.y, maxT = Math.acos(max(-1, minXi - NLA_PRECISION));
                return [
                    new EllipseCurve$1(plane.anchor, f1.negated(), f2, PI - maxT, PI),
                    new EllipseCurve$1(plane.anchor, f1, f2.negated(), 0, maxT),
                ];
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
            if (le(1, surfaceCenterDist - surfaceRadius) ||
                le(surfaceCenterDist + surfaceRadius, 1) ||
                le(surfaceCenterDist - surfaceRadius, -1)) {
                return [];
            }
            else {
                // origin, surface.center and points on the intersection curves form a triangle.
                // the height on the segment origin - surface.center is the radius of the is curves
                // the distance from the origin to the lot point is the distance to the intersection plane
                function heron(a, b, c) {
                    const p = (a + b + c) / 2;
                    return sqrt(p * (p - a) * (p - b) * (p - c));
                }
                const triangleArea = heron(1, surfaceRadius, surfaceCenterDist);
                const radius = (triangleArea * 2) / surfaceCenterDist;
                const isCurvesCenterDist = sign(1 + Math.pow(surfaceCenterDist, 2) - Math.pow(surfaceRadius, 2)) *
                    sqrt(1 - Math.pow(radius, 2));
                const plane = new P3$1(surface.center.unit(), isCurvesCenterDist);
                return EllipsoidSurface.unitISCurvesWithPlane(plane.flipped());
            }
        }
        throw new Error();
    }
    static unitISCurvesWithCylinderSurface(surface) {
        if (new L3$1(surface.baseCurve.center, surface.dir).containsPoint(V3.O)) {
            const projEllipse = surface.baseCurve.transform(M4.project(new P3$1(surface.dir, 0)));
            const f1Length = projEllipse.f1.length(), f2Length = projEllipse.f2.length();
            if (lt(1, min(f1Length, f2Length)))
                return [];
            if (projEllipse.isCircular()) {
                const distISCurveCenter = Math.sqrt(1 - Math.pow(min(1, f1Length), 2));
                const isCurveCenter = (surface.dir.y < 0
                    ? surface.dir.negated()
                    : surface.dir).times(distISCurveCenter);
                // isCurve.at(t).y = isCurveCenter.y + projEllipse.f1.y * cos(t) + projEllipse.f2.y * sin(t) = 0
                return [new EllipseCurve$1(isCurveCenter, projEllipse.f1, projEllipse.f2)];
            }
        }
        throw new Error();
    }
    static sphere(radius, center = V3.O) {
        assertNumbers(radius);
        return new EllipsoidSurface(center, new V3(radius, 0, 0), new V3(0, radius, 0), new V3(0, 0, radius));
    }
    /**
     * x²/a² + y²/b² + z²/c² = 1
     */
    static forABC(a, b, c, center = V3.O) {
        return new EllipsoidSurface(center, new V3(a, 0, 0), new V3(0, b, 0), new V3(0, 0, c));
    }
    static calculateAreaSpheroid(a, b, c, edges) {
        assertf(() => a.isPerpendicularTo(b));
        assertf(() => b.isPerpendicularTo(c));
        assertf(() => c.isPerpendicularTo(a));
        // handling discontinuities:
        // option 1: check for intersections with baseline, if there are any integrate parts separetely
        // "rotate" the edge so that there are no overlaps
        const matrix = M4.forSys(a, b, c), matrixInverse = matrix.inversed();
        const circleRadius = a.length();
        const c1 = c.unit();
        const totalArea = edges
            .map((edge) => {
            if (edge.curve instanceof EllipseCurve$1) {
                const f = (t) => {
                    const at = edge.curve.at(t), tangent = edge.tangentAt(t);
                    const localAt = matrixInverse.transformPoint(at);
                    const angleXY = localAt.angleXY();
                    const arcLength = angleXY * circleRadius * Math.sqrt(1 + Math.pow(localAt.z, 2));
                    const scaling = Math.sqrt(1 + Math.pow(c1.dot(tangent), 2));
                    return arcLength * scaling;
                };
                const val = glqInSteps(f, edge.aT, edge.bT, 1);
                return val;
            }
            else {
                throw new Error();
            }
        })
            .sum();
        return totalArea;
    }
    getConstructorParameters() {
        return [
            this.center,
            this.f1,
            this.f2,
            this.f3,
            this.uMin,
            this.uMax,
            this.vMin,
            this.vMax,
        ];
    }
    equals(obj) {
        return (this == obj ||
            (Object.getPrototypeOf(obj) == this.constructor.prototype &&
                this.matrix.equals(obj.matrix)));
    }
    edgeLoopCCW(loop) {
        return (EllipsoidSurface.unitArea(loop.map((edge) => edge.transform(this.matrixInverse))) > 0);
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
    rootPoints() { }
    toMesh() {
        return ParametricSurface$1.prototype.toMesh.call(this);
    }
    clipCurves(curves) {
        return curves.flatMap((curve) => curve.clipPlane(this.getSeamPlane()));
    }
    dpdu() {
        // dp(u, v) = new V3(cos(t) * cos(s), cos(t) * sin(s), sin(t)
        return (u, v) => this.matrix.transformVector(new V3(cos(v) * -sin(u), cos(v) * cos(u), 0));
    }
    dpdv() {
        return (u, v) => this.matrix.transformVector(new V3(-sin(v) * cos(u), -sin(v) * sin(u), cos(v)));
    }
    isCurvesWithPCS(surface) {
        let curves2 = ParametricSurface$1.isCurvesParametricImplicitSurface(surface, this, 0.1, 0.1 / surface.dir.length(), 0.05);
        curves2 = this.clipCurves(curves2);
        return curves2;
    }
    isCurvesWithPCSSmart(surface) {
        const surfaceLC = surface.transform(this.matrixInverse);
        //const lcMinZ0RelO =
        const baseCurveLC = surfaceLC.baseCurve.project(new P3$1(surfaceLC.dir, 0));
        const ists = baseCurveLC.isTsWithSurface(EllipsoidSurface.UNIT);
        const insideIntervals = getIntervals(ists, baseCurveLC.tMin, baseCurveLC.tMax).filter(([a, b]) => baseCurveLC.at((a + b) / 2).length() < 1);
        const projectedCurves = [0, 1].map((id) => {
            return (t) => {
                const atSqr = snap(baseCurveLC.at(t).squared(), 1);
                const lineISTs = /* +- */ sqrt(1 - atSqr);
                //assert(!isNaN(lineISTs))
                return eq0(lineISTs)
                    ? baseCurveLC.at(t)
                    : baseCurveLC
                        .at(t)
                        .plus(surfaceLC.dir.times(sign(id - 0.5) * lineISTs));
            };
        });
        const dProjectedCurves = [0, 1].map((id) => {
            return (t) => {
                // d/dt sqrt(1 - baseCurveLC.at(t).squared())
                // = -1/2 * 1/sqrt(1 - baseCurveLC.at(t).squared()) * -2*baseCurveLC.at(t) * baseCurveLC.tangentAt(t)
                const atSqr = snap(baseCurveLC.at(t).squared(), 1);
                const lineISTs = /* +- */ baseCurveLC
                    .at(t)
                    .times(-1 / sqrt(1 - atSqr))
                    .dot(baseCurveLC.tangentAt(t));
                //assert(!isNaN(lineISTs))
                return baseCurveLC
                    .tangentAt(t)
                    .plus(surfaceLC.dir.times(sign(id - 0.5) * lineISTs));
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
                checkDerivate(f, df, aT + 0.1, bT - 0.1);
                const tsAtY0 = getRoots(f, aT + NLA_PRECISION, bT - NLA_PRECISION, 1 / (1 << 11), df);
                const ii2 = getIntervals(tsAtY0, aT, bT).filter(([a, b]) => f((a + b) / 2) > 0);
                for (const [aT2, bT2] of ii2) {
                    let aP = projectedCurves[i](aT2), bP = projectedCurves[i](bT2);
                    0 === i && ([aP, bP] = [bP, aP]);
                    assert(EllipsoidSurface.UNIT.containsPoint(aP));
                    assert(EllipsoidSurface.UNIT.containsPoint(bP));
                    curves.push(PICurve$1.forStartEnd(surface, this, this.matrix.transformPoint(bP), this.matrix.transformPoint(aP), undefined));
                }
            }
        }
        return surface.clipCurves(curves);
    }
    isCurvesWithSurface(surface) {
        if (surface instanceof PlaneSurface$1) {
            return this.isCurvesWithPlane(surface.plane);
        }
        else if (surface instanceof CylinderSurface$1) {
            return this.isCurvesWithCylinderSurface(surface);
        }
        else if (surface instanceof EllipsoidSurface) {
            const surfaceLC = surface.transform(this.matrixInverse);
            const curves = EllipsoidSurface.unitISCurvesWithEllipsoidSurface(surfaceLC).map((c) => c.transform(this.matrix));
            return surface.clipCurves(curves);
        }
        else if (surface instanceof ProjectedCurveSurface$1) {
            return this.isCurvesWithPCS(surface);
        }
        else if (surface instanceof ParametricSurface$1) {
            let curves2 = ParametricSurface$1.isCurvesParametricImplicitSurface(surface, this, 0.1, 0.1, 0.05);
            curves2 = this.clipCurves(curves2);
            curves2 = surface.clipCurves(curves2);
            return curves2;
        }
        else {
            throw new Error();
        }
    }
    isCurvesWithPlane(plane) {
        const planeLC = plane.transform(this.matrixInverse);
        return EllipsoidSurface.unitISCurvesWithPlane(planeLC).map((c) => c.transform(this.matrix));
    }
    isCurvesWithCylinderSurface(surface) {
        if (L3$1.containsPoint(surface.baseCurve.center, surface.dir, this.center)) {
            assert(this.isSphere());
            const ellipseProjected = surface.baseCurve.transform(M4.project(surface.baseCurve.getPlane(), surface.dir));
            if (ellipseProjected.isCircular()) {
                const thisRadius = this.f1.length();
                const surfaceRadius = ellipseProjected.f1.length();
                // sphereRadius² = distanceISFromCenter² + isRadius²
                if (eq(thisRadius, surfaceRadius)) ;
                assert(false);
            }
        }
        return this.isCurvesWithPCS(surface);
    }
    isTsForLine(line) {
        assertInst(L3$1, line);
        // transforming line manually has advantage that dir1 will not be renormalized,
        // meaning that calculated values t for localLine are directly transferable to line
        const anchorLC = this.matrixInverse.transformPoint(line.anchor);
        const dirLC = this.matrixInverse.transformVector(line.dir1);
        return EllipsoidSurface.unitISTsWithLine(anchorLC, dirLC);
    }
    isCoplanarTo(surface) {
        if (this === surface)
            return true;
        if (!hasConstructor(surface, EllipsoidSurface))
            return false;
        if (!this.center.like(surface.center))
            return false;
        if (this.isSphere())
            return surface.isSphere() && eq(this.f1.length(), this.f2.length());
        const otherMatrixLC = this.matrixInverse.times(surface.matrix);
        // Ellipsoid with matrix otherMatrixLC is unit sphere iff otherMatrixLC is orthogonal
        return otherMatrixLC.like3x3() && otherMatrixLC.isOrthogonal();
    }
    containsEllipse(ellipse) {
        const ellipseLC = ellipse.transform(this.matrixInverse);
        const distEllipseLCCenter = ellipseLC.center.length();
        const correctRadius = Math.sqrt(1 - Math.pow(distEllipseLCCenter, 2));
        return (lt(distEllipseLCCenter, 1) &&
            ellipseLC.isCircular() &&
            ellipseLC.f1.hasLength(correctRadius));
        //&& le(0, ellipseLC.getAABB().min.y)
    }
    containsCurve(curve) {
        if (curve instanceof EllipseCurve$1) {
            return this.containsEllipse(curve);
        }
        else {
            return super.containsCurve(curve);
        }
    }
    transform(m4) {
        assert(m4.isNoProj(), () => m4.sce);
        return new EllipsoidSurface(m4.transformPoint(this.center), m4.transformVector(this.f1), m4.transformVector(this.f2), m4.transformVector(this.f3).times(m4.isMirroring() ? -1 : 1));
    }
    transform4(m4) {
        console.log("transform4");
        const resultMatrix = m4.times(this.matrix);
        console.log(resultMatrix.toString());
        const scaleDir = V(resultMatrix.m[12], resultMatrix.m[13], resultMatrix.m[14]);
        // need to find parameters where scaleDir is parallel to the normal
        const pLC = this.pLCNormalWCMatrix.inversed().transformPoint(scaleDir);
        const s = pLC.angleXY();
        const t = Math.asin(clamp(pLC.z, -1, 1));
        const fa = resultMatrix.transformPoint(scaleDir.unit());
        const fb = resultMatrix.transformPoint(scaleDir.unit().negated());
        const newCenter = V3.lerp(fa, fb, 0.5);
        console.log(scaleDir.sce, s, t, fa, fb, "newCenter", newCenter.sce);
        return new EllipsoidSurface(newCenter, m4.transformVector2(this.f1, this.center), m4.transformVector2(this.f2, this.center), m4
            .transformVector2(this.f3, this.center)
            .times(m4.isMirroring() ? -1 : 1));
    }
    isInsideOut() {
        return this.f1.cross(this.f2).dot(this.f3) < 0;
    }
    flipped() {
        return new EllipsoidSurface(this.center, this.f1, this.f2, this.f3.negated(), this.uMin, this.uMax, -this.vMax, -this.vMin);
    }
    normalUVFunc() {
        // ugh
        // paramtric ellipsoid point q(a, b)
        // normal1 == (dq(a, b) / da) X (dq(a, b) / db) (cross product of partial derivatives)
        // normal1 == cos b * (f2 X f3 * cos b * cos a + f3 X f1 * cos b * sin a + f1 X f2 * sin b)
        return (a, b) => {
            const { f1, f2, f3 } = this;
            const normal = f2
                .cross(f3)
                .times(Math.cos(b) * Math.cos(a))
                .plus(f3.cross(f1).times(Math.cos(b) * Math.sin(a)))
                .plus(f1.cross(f2).times(Math.sin(b)))
                //.times(Math.cos(b))
                .unit();
            return normal;
        };
    }
    normalP(p) {
        return this.pLCNormalWCMatrix
            .transformVector(this.matrixInverse.transformPoint(p))
            .unit();
    }
    normalUV(u, v) {
        return this.pLCNormalWCMatrix.transformVector(V3.sphere(u, v)).unit();
    }
    uvPFunc() {
        return (pWC) => {
            const pLC = this.matrixInverse.transformPoint(pWC);
            const alpha = abs(pLC.angleXY());
            const beta = Math.asin(clamp(pLC.z, -1, 1));
            assert(isFinite(alpha));
            assert(isFinite(beta));
            return new V3(alpha, beta, 0);
        };
    }
    pUVFunc() {
        // this(a, b) = f1 cos a cos b + f2 sin a cos b + f2 sin b
        return (alpha, beta) => {
            return this.matrix.transformPoint(V3.sphere(alpha, beta));
        };
    }
    isSphere() {
        return (eq(this.f1.length(), this.f2.length()) &&
            eq(this.f2.length(), this.f3.length()) &&
            eq(this.f3.length(), this.f1.length()) &&
            this.f1.isPerpendicularTo(this.f2) &&
            this.f2.isPerpendicularTo(this.f3) &&
            this.f3.isPerpendicularTo(this.f1));
    }
    isVerticalSpheroid() {
        return (eq(this.f1.length(), this.f2.length()) &&
            this.f1.isPerpendicularTo(this.f2) &&
            this.f2.isPerpendicularTo(this.f3) &&
            this.f3.isPerpendicularTo(this.f1));
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
        if (eq0(f1.dot(f2)) && eq0(f2.dot(f3)) && eq0(f3.dot(f1))) {
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
        //const mainF1Params = newtonIterate(f, [0, 0], 8), mainF1 = this.pUVFunc()(mainF1Params[0], mainF1Params[1])
        //console.log(f(mainF1Params, 1).sce)
        //const mainF2Params = newtonIterate(f, this.uvPFunc()(f2.rejectedFrom(mainF1)).toArray(2), 8),
        //   mainF2 = this.pUVFunc()(mainF2Params[0], mainF2Params[1])
        //console.log(this.normalUVFunc()(mainF2Params[0], mainF2Params[1]).sce)
        //assert(mainF1.isPerpendicularTo(mainF2), mainF1, mainF2, mainF1.dot(mainF2), mainF1Params)
        //const mainF3Params = this.uvPFunc()(mainF1.cross(mainF2)), mainF3 = this.pUVFunc()(mainF3Params[0],
        // mainF3Params[1]) return new EllipsoidSurface(this.center, mainF1, mainF2, mainF3)
        const { U, SIGMA } = this.matrix.svd3();
        assert(SIGMA.isDiagonal());
        assert(U.isOrthogonal());
        const U_SIGMA = U.times(SIGMA);
        // column vectors of U_SIGMA
        const [mainF1, mainF2, mainF3] = arrayFromFunction(3, (i) => new V3(U_SIGMA.m[i], U_SIGMA.m[i + 4], U_SIGMA.m[i + 8]));
        return new EllipsoidSurface(this.center, mainF1, mainF2, mainF3);
    }
    containsPoint(p) {
        return eq0(this.implicitFunction()(p));
    }
    boundsFunction() {
        return (a, b) => between(a, 0, PI) && between(b, -PI, PI);
    }
    volume() {
        return (4 / 3) * Math.PI * this.f1.dot(this.f2.cross(this.f3));
    }
    loopContainsPoint(loop, pWC) {
        if (!this.containsPoint(pWC))
            return PointVsFace$1.OUTSIDE;
        assertVectors(pWC);
        assert(Edge$1.isLoop(loop));
        const pLCXY = this.matrixInverse.transformPoint(pWC).xy();
        const testLine = new EllipseCurve$1(this.center, this.f3, pLCXY.likeO() ? this.f2 : this.matrix.transformVector(pLCXY.unit()));
        if (P3$1.normalOnAnchor(this.f2.unit(), this.center).containsPoint(pWC)) {
            return loop.some((edge) => edge.curve.containsPoint(pWC) &&
                fuzzyBetween(edge.curve.pointT(pWC), edge.minT, edge.maxT))
                ? PointVsFace$1.ON_EDGE
                : PointVsFace$1.OUTSIDE;
        }
        return Surface$1.loopContainsPointEllipse(loop, pWC, testLine);
    }
    surfaceAreaApprox() {
        // See https://en.wikipedia.org/wiki/Ellipsoid#Surface_area
        const mainAxes = this.mainAxes(), a = mainAxes.f1.length(), b = mainAxes.f2.length(), c = mainAxes.f3.length();
        const p = 1.6075;
        return (4 *
            PI *
            Math.pow((Math.pow(a * b, p) + Math.pow(b * c, p) + Math.pow(c * a, p)) / 3, 1 / p));
    }
    surfaceArea() {
        // See https://en.wikipedia.org/wiki/Ellipsoid#Surface_area
        const mainAxes = this.mainAxes(), f1l = mainAxes.f1.length(), f2l = mainAxes.f2.length(), f3l = mainAxes.f3.length(), [c, b, a] = [f1l, f2l, f3l].sort(MINUS);
        // https://en.wikipedia.org/w/index.php?title=Spheroid&oldid=761246800#Area
        function spheroidArea(a, c) {
            if (c < a) {
                const eccentricity2 = 1 - Math.pow(c, 2) / Math.pow(a, 2);
                const eccentricity = Math.sqrt(eccentricity2);
                return (2 *
                    PI *
                    Math.pow(a, 2) *
                    (1 +
                        ((1 - eccentricity2) / Math.sqrt(eccentricity)) *
                            Math.atanh(eccentricity)));
            }
            else {
                const eccentricity = Math.sqrt(1 - Math.pow(a, 2) / Math.pow(c, 2));
                return (2 *
                    PI *
                    Math.pow(a, 2) *
                    (1 + (c / a / eccentricity) * Math.asin(eccentricity)));
            }
        }
        if (eq(a, b)) {
            return spheroidArea(a, c);
        }
        else if (eq(b, c)) {
            return spheroidArea(b, a);
        }
        else if (eq(c, a)) {
            return spheroidArea(c, b);
        }
        const phi = Math.acos(c / a);
        const kk = (Math.pow(a, 2) * (Math.pow(b, 2) - Math.pow(c, 2))) / (Math.pow(b, 2) * (Math.pow(a, 2) - Math.pow(c, 2)));
        const incompleteEllipticInt1 = gaussLegendreQuadrature24((phi) => Math.pow(1 - kk * Math.pow(Math.sin(phi), 2), -0.5), 0, phi);
        const incompleteEllipticInt2 = gaussLegendreQuadrature24((phi) => Math.pow(1 - kk * Math.pow(Math.sin(phi), 2), 0.5), 0, phi);
        return ((2 * PI * Math.pow(c, 2) + (2 * PI * a * b) / Math.sin(phi)) *
            (incompleteEllipticInt2 * Math.pow(Math.sin(phi), 2) +
                incompleteEllipticInt1 * Math.pow(Math.cos(phi), 2)));
    }
    getSeamPlane() {
        const plane = P3$1.forAnchorAndPlaneVectors(this.center, this.f1, this.f3);
        return plane.normal1.dot(this.f2) < 0 ? plane : plane.flipped();
    }
    getExtremePoints() {
        return getExtremePointsHelper$1.call(this, new EllipseCurve$1(V3.O, V3.X, V3.Z, -PI / 2, PI / 2));
    }
    pointFoot(pWC, startS, startT) {
        console.log(pWC.sce);
        if (undefined === startS || undefined === startT) {
            let pLC1 = this.matrixInverse.transformPoint(pWC).unit();
            if (pLC1.y < 0)
                pLC1 = pLC1.negated();
            ({ x: startS, y: startT } = EllipsoidSurface.UNIT.uvP(pLC1));
        }
        const dpdu = this.dpdu();
        const dpdv = this.dpdv();
        const [u, v] = newtonIterate(([u, v]) => {
            const p = this.pUV(u, v);
            console.log([p, p.plus(dpdu(u, v)), p, p.plus(dpdv(u, v))].map(toSource).join() +
                ",");
            const pUVToPWC = this.pUV(u, v).to(pWC);
            return [pUVToPWC.dot(dpdu(u, v)), pUVToPWC.dot(dpdv(u, v))];
        }, [startS, startT], 8, undefined, 0.1);
        return new V3(u, v, 0);
    }
    implicitFunction() {
        return (pWC) => {
            const pLC = this.matrixInverse.transformPoint(pWC);
            return (pLC.length() - 1) * this.normalDir;
        };
    }
    // = this.inverseMatrix.transformPoint(this.inverseMatrix.transformPoint(pWC).unit())
    didp(pWC) {
        // i(pWC) = this.inverseMatrix.transformPoint(pWC).length() - 1
        // chain diff rule
        const pLC = this.matrixInverse.transformPoint(pWC);
        return this.pLCNormalWCMatrix.transformVector(pLC.unit()); //.times(this.normalDir)
    }
    /*+
     * An ellipsoid remains an ellipsoid after a perspective transform (as long as it does not intersect the vanishing
     * plane. This transforms a matrix with a perspective component into one which would return an identical ellipsoid,
     * but with no perspective component.
     */
    static unitTransform4(m) {
        m.m[15] !== 1 && (m = m.divScalar(m.m[15]));
        // X * P = m => X = m * P^-1
        // prettier-ignore
        const Pinv = new M4(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, -m.m[12], -m.m[13], -m.m[14], 1);
        const pn = new V3(m.m[12], m.m[13], m.m[14]), pw = m.m[15];
        const pwSqrMinusPnSqr = Math.pow(pw, 2) - pn.squared();
        if (lt(pwSqrMinusPnSqr, 0)) {
            throw new Error("vanishing plane intersects unit sphere");
        }
        const c = pn.div(-pwSqrMinusPnSqr);
        const scale = pn.times((pw * pn.length()) / (pn.squared() * -pwSqrMinusPnSqr));
        const scale1 = pw / -pwSqrMinusPnSqr;
        const scale2 = 1 / sqrt(pwSqrMinusPnSqr);
        const rotNX = M4.forSys(pn.unit(), pn.getPerpendicular().unit());
        return M4.product(m, Pinv, M4.translate(c), rotNX, M4.scale(scale1, scale2, scale2), rotNX.transposed());
    }
}
EllipsoidSurface.UNIT = new EllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z);
EllipsoidSurface.prototype.uStep = PI / 32;
EllipsoidSurface.prototype.vStep = PI / 32;

class PlaneSurface extends ParametricSurface$1 {
    constructor(plane, right = plane.normal1.getPerpendicular().unit(), up = plane.normal1.cross(right).unit(), uMin = -100, uMax = 100, vMin = -100, vMax = 100) {
        super(uMin, uMax, vMin, vMax);
        this.plane = plane;
        this.right = right;
        this.up = up;
        assertInst(P3$1, plane);
        assert(this.right.cross(this.up).like(this.plane.normal1));
        this.matrix = M4.forSys(right, up, plane.normal1, plane.anchor);
    }
    toSource(rounder = (x) => x) {
        return callsce.call(undefined, "new PlaneSurface", ...this.getConstructorParameters());
    }
    static throughPoints(a, b, c) {
        return new PlaneSurface(P3$1.throughPoints(a, b, c));
    }
    static forAnchorAndPlaneVectors(anchor, v0, v1, uMin, uMax, vMin, vMax) {
        return new PlaneSurface(P3$1.forAnchorAndPlaneVectors(anchor, v0, v1), v0, v1, uMin, uMax, vMin, vMax);
    }
    isCoplanarTo(surface) {
        return (hasConstructor(surface, PlaneSurface) &&
            this.plane.isCoplanarToPlane(surface.plane));
    }
    isTsForLine(line) {
        return line.isTsWithPlane(this.plane);
    }
    like(surface) {
        return (hasConstructor(surface, PlaneSurface) && this.plane.like(surface.plane));
    }
    pUV(u, v) {
        return this.matrix.transformPoint(new V3(u, v, 0));
    }
    implicitFunction() {
        return (p) => this.plane.distanceToPointSigned(p);
    }
    isCurvesWithSurface(surface2) {
        if (surface2 instanceof PlaneSurface) {
            return this.isCurvesWithPlane(surface2.plane);
        }
        return super.isCurvesWithSurface(surface2);
    }
    isCurvesWithPlane(plane) {
        const result = this.plane.intersectionWithPlane(plane);
        return result ? [result] : [];
    }
    edgeLoopCCW(contour) {
        assert(Edge$1.isLoop(contour), "isLoop");
        return isCCW(contour.flatMap((edge) => edge.points()), this.plane.normal1);
    }
    loopContainsPoint(loop, p) {
        const dir = this.right.plus(this.up.times(0.123)).unit();
        const line = new L3$1(p, dir);
        const lineOut = dir.cross(this.plane.normal1);
        return Surface$1.loopContainsPointGeneral(loop, p, line, lineOut);
    }
    uvPFunc() {
        const matrixInverse = this.matrix.inversed();
        return function (pWC) {
            return matrixInverse.transformPoint(pWC);
        };
    }
    pointFoot(pWC) {
        return this.uvP(pWC);
    }
    normalP(pWC) {
        return this.plane.normal1;
    }
    containsPoint(p) {
        return this.plane.containsPoint(p);
    }
    containsCurve(curve) {
        return curve instanceof ImplicitCurve$1
            ? super.containsCurve(curve)
            : this.plane.containsCurve(curve);
    }
    transform(m4) {
        return new PlaneSurface(this.plane.transform(m4));
    }
    transform4(m4) {
        return new PlaneSurface(this.plane.transform(m4));
    }
    flipped() {
        return new PlaneSurface(this.plane.flipped(), this.right, this.up.negated());
    }
    getConstructorParameters() {
        return [
            this.plane,
            this.right,
            this.up,
            this.uMin,
            this.uMax,
            this.vMin,
            this.vMax,
        ];
    }
    dpdu() {
        return () => this.right;
    }
    dpdv() {
        return () => this.up;
    }
    didp(pWC) {
        return this.plane.normal1;
    }
    normalUV() {
        return this.plane.normal1;
    }
}
PlaneSurface.prototype.uStep = 1e6;
PlaneSurface.prototype.vStep = 1e6;

class PointProjectedSurface extends ParametricSurface$1 {
    constructor(curve, apex, curvePlane, normalDir = 1, uMin = curve.tMin, uMax = curve.tMax, vMin = 0, vMax = 16) {
        super(uMin, uMax, vMin, vMax);
        this.curve = curve;
        this.apex = apex;
        this.curvePlane = curvePlane;
        this.normalDir = normalDir;
        assertInst(Curve$1, curve);
        assert(!(curve instanceof L3$1), "use PlaneSurface instead");
        assert(!(curve instanceof EllipseCurve$1), "use ConicSurface instead");
        assert(!(curve instanceof ImplicitCurve$1), "this just seems like a terrible idea");
        assert(new PlaneSurface$1(curvePlane).containsCurve(curve));
        assertVectors(apex);
        assert(0 <= vMin);
        this.planeProjectionMatrix = M4.projectPlanePoint(apex, curvePlane);
        this.uStep = curve.tIncrement;
    }
    pointFoot(pWC, ss, st) {
        if (undefined === ss || undefined === st) {
            // similar to stP
            if (undefined === ss) {
                ss = pWC.like(this.apex)
                    ? 0
                    : this.curve.closestTToPoint(this.planeProjectionMatrix.transformPoint(pWC)) * this.normalDir;
            }
            if (undefined === st) {
                st = V3.inverseLerp(this.apex, this.curve.at(ss), pWC);
            }
        }
        const f = ([s, t]) => {
            const pSTToPWC = this.pST(s, t).to(pWC);
            return [this.dpds()(s, t).dot(pSTToPWC), this.dpdt()(s).dot(pSTToPWC)];
        };
        const { 0: x, 1: y } = newtonIterate(f, [ss, st]);
        return new V3(x, y, 0);
    }
    getConstructorParameters() {
        return [
            this.curve,
            this.apex,
            this.curvePlane,
            this.normalDir,
            this.uMin,
            this.uMax,
            this.vMin,
            this.vMax,
        ];
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
        return pqFormula(b / a, c / a).filter((t) => 0 < az + t * dz);
    }
    equals(obj) {
        return (this == obj ||
            (hasConstructor(obj, PointProjectedSurface) &&
                this.curve.equals(obj.curve) &&
                this.apex.equals(this.apex)));
    }
    like(object) {
        if (!this.isCoplanarTo(object))
            return false;
        // normals need to point in the same direction (outwards or inwards) for both
        return this.normalDir == object.normalDir;
    }
    loopContainsPoint(contour, p) {
        assertVectors(p);
        const line = this.apex.like(p)
            ? new L3$1(p, this.apex.to(this.curve.at(this.curve.tMin)).unit())
            : L3$1.throughPoints(p, this.apex);
        const lineOut = line.dir1.cross(this.curvePlane.normal1);
        return Surface$1.loopContainsPointGeneral(contour, p, line, lineOut);
    }
    isTsForLine(line) {
        // transforming line manually has advantage that dir1 will not be renormalized,
        // meaning that calculated values t for lineLC are directly transferable to line
        const anchorPlane = this.planeProjectionMatrix.transformPoint(line.anchor);
        const anchor2Plane = this.planeProjectionMatrix.transformPoint(line.anchor.plus(line.dir1));
        if (anchorPlane.like(anchor2Plane)) {
            // line projects onto a point in plane.
            // there are either no or infinite intersection points
            return [];
        }
        return this.curve
            .isInfosWithLine(anchorPlane, anchorPlane.to(anchor2Plane), undefined, undefined, line.tMin, line.tMax)
            .map((info) => info.tOther);
    }
    /**
     * Interestingly, two cones don't need to have parallel dirs to be coplanar.
     */
    isCoplanarTo(surface) {
        if (this === surface)
            return true;
        if (!(surface instanceof PointProjectedSurface) ||
            !this.apex.like(surface.apex))
            return false;
        // at this point apexes are equal
        return this.containsCurve(surface.curve);
    }
    containsLine(line) {
        if (this.curvePlane.isParallelToLine(line)) {
            return false;
        }
        if (!line.containsPoint(this.apex)) {
            return false;
        }
        const p = this.curvePlane.intersectionWithLine(line);
        return this.curve.containsPoint(p);
    }
    containsCurve(curve) {
        if (curve instanceof L3$1) {
            return this.containsLine(curve);
        }
        else if (!(curve instanceof ImplicitCurve$1)) {
            const otherCurveOnThisPlane = curve.transform(this.planeProjectionMatrix);
            return this.curve.isColinearTo(otherCurveOnThisPlane);
        }
        else {
            return super.containsCurve(curve);
        }
    }
    transform(m4) {
        return new PointProjectedSurface(this.curve.transform(m4), m4.transformPoint(this.apex), this.curvePlane.transform(m4), (m4.isMirroring() ? -1 : 1) * this.normalDir, this.uMin, this.uMax, this.vMin, this.vMax);
    }
    flipped() {
        return new PointProjectedSurface(this.curve, this.apex, this.curvePlane, -this.normalDir, -this.uMax, -this.uMin, this.vMin, this.vMax);
    }
    normalSTFunc() {
        const dpdt = this.dpdt();
        return (s, t) => this.curve
            .tangentAt(s * this.normalDir)
            .times(this.normalDir)
            .cross(dpdt(s))
            .unit();
    }
    pSTFunc() {
        return (s, t) => {
            return this.apex.lerp(this.curve.at(s * this.normalDir), t);
        };
    }
    dpds() {
        return (s, t) => {
            return this.curve.tangentAt(s * this.normalDir).times(t * this.normalDir);
        };
    }
    dpdt() {
        return (s) => {
            return this.apex.to(this.curve.at(s * this.normalDir));
        };
    }
    containsPoint(pWC) {
        return (this.apex.like(pWC) ||
            this.curve.containsPoint(this.planeProjectionMatrix.transformPoint(pWC)));
    }
    stP(pWC) {
        const s = pWC.like(this.apex)
            ? 0
            : this.curve.pointT(this.planeProjectionMatrix.transformPoint(pWC));
        const t = V3.inverseLerp(this.apex, this.curve.at(s), pWC);
        return new V3(s * this.normalDir, t, 0);
    }
    isCurvesWithSurface(surface) {
        if (surface instanceof PlaneSurface$1) {
            return this.isCurvesWithPlane(surface.plane);
        }
        else if (ImplicitSurface$1.is(surface)) {
            return ParametricSurface$1.isCurvesParametricImplicitSurface(this, surface, 0.1, 0.1 / this.curvePlane.distanceToPoint(this.apex), 0.02);
        }
        return super.isCurvesWithSurface(surface);
    }
    isCurvesWithPlane(plane) {
        if (plane.containsPoint(this.apex)) {
            if (plane.isParallelToPlane(this.curvePlane)) {
                return [];
            }
            return this.curve
                .isTsWithPlane(plane)
                .map((t) => L3$1.throughPoints(this.apex, this.curve.at(t)));
        }
        return [this.curve.transform(M4.projectPlanePoint(this.apex, plane))];
    }
}
PointProjectedSurface.prototype.vStep = 256;

class NURBSSurface extends ParametricSurface$1 {
    constructor(
    /**
     * Control points in u-major order. I.e. the first pointCountU points are a NURBS.
     */
    points, knotsU, knotsV, degreeU, degreeV, uMin = knotsU[degreeU], uMax = knotsU[knotsU.length - degreeU - 1], vMin = knotsV[degreeV], vMax = knotsV[knotsV.length - degreeV - 1]) {
        super(uMin, uMax, vMin, vMax);
        this.points = points;
        this.knotsU = knotsU;
        this.knotsV = knotsV;
        this.degreeU = degreeU;
        this.degreeV = degreeV;
        const pointCountU = knotsU.length - 1 - degreeU;
        const pointCountV = knotsV.length - 1 - degreeV;
        assert(pointCountU * pointCountV == points.length);
        assert(degreeU <= degreeV, "degreeU <= degreeV");
        assert(-1 === firstUnsorted(knotsU, MINUS), "knot values must be in ascending order");
        assert(-1 === firstUnsorted(knotsV, MINUS), "knot values must be in ascending order");
    }
    getConstructorParameters() {
        return [
            this.points,
            this.knotsU,
            this.knotsV,
            this.degreeU,
            this.degreeV,
            this.uMin,
            this.uMax,
            this.vMin,
            this.vMax,
        ];
    }
    transform(m4) {
        return this.transform4(m4);
    }
    transform4(m4) {
        return new NURBSSurface(this.points.map((p) => m4.timesVector(p)), this.knotsU, this.knotsV, this.degreeU, this.degreeV, this.uMin, this.uMax, this.vMin, this.vMax);
    }
    pUV(u, v) {
        return this.isoparametricU(u).at(v);
    }
    dpdu() {
        return (u, v) => this.isoparametricV(v).tangentAt(u);
    }
    dpdv() {
        return (u, v) => this.isoparametricU(u).tangentAt(v);
    }
    normalUV(u, v) {
        const normal = this.dpdu()(u, v).cross(this.dpdv()(u, v));
        return normal.likeO() ? V3.X : normal.unit();
    }
    isoparametricU(u) {
        const pointCountU = this.knotsU.length - 1 - this.degreeU;
        const pointCountV = this.knotsV.length - 1 - this.degreeV;
        return new NURBS$1(arrayFromFunction(pointCountV, (i) => {
            return deBoor(this.points.slice(i * pointCountU, (i + 1) * pointCountU), this.degreeU, this.knotsU, u);
        }), this.degreeV, this.knotsV, this.vMin, this.vMax);
    }
    isoparametricV(v) {
        const pointCountU = this.knotsU.length - 1 - this.degreeU;
        return new NURBS$1(arrayFromFunction(pointCountU, (i) => {
            return deBoor(sliceStep(this.points, i, this.points.length, pointCountU, 1), this.degreeV, this.knotsV, v);
        }), this.degreeU, this.knotsU, this.uMin, this.uMax);
    }
    debugInfo() {
        const pointCountU = this.knotsU.length - 1 - this.degreeU;
        const pointCountV = this.knotsV.length - 1 - this.degreeV;
        const grid = [];
        for (let u = 0; u < pointCountU; u++) {
            for (let v = 0; v < pointCountV; v++) {
                const i = v * pointCountU + u;
                if (u < pointCountU - 1) {
                    const j = v * pointCountU + u + 1;
                    grid.push(this.points[i].p3(), this.points[j].p3());
                }
                if (v < pointCountV - 1) {
                    const j = (v + 1) * pointCountU + u;
                    grid.push(this.points[i].p3(), this.points[j].p3());
                }
            }
        }
        return { points: this.points.map((p) => p.p3()), lines: grid };
    }
    flipped() {
        const pointCountU = this.knotsU.length - 1 - this.degreeU;
        return new NURBSSurface(arrayFromFunction(this.points.length, (i) => {
            const u = i % pointCountU;
            return this.points[i - u + (pointCountU - u - 1)];
        }), this.knotsU.map((x) => -x).reverse(), this.knotsV, this.degreeU, this.degreeV, -this.uMax, -this.uMin, this.vMin, this.vMax);
    }
}
NURBSSurface.prototype.uStep = 1 / 8;
NURBSSurface.prototype.vStep = 1 / 8;
function getInterval(degree, knots, t) {
    for (let s = degree; s < knots.length - 1 - degree; s++) {
        if (t >= knots[s] && t <= knots[s + 1]) {
            return s;
        }
    }
    throw new Error(t + " " + knots);
}
function deBoor(points, degree, knots, t) {
    // find s (the spline segment) for the [t] value provided
    const s = getInterval(degree, knots, t);
    const v = Vector.pack(points, new Float64Array(points.length * 4));
    // l (level) goes from 1 to the curve degree + 1
    for (let l = 1; l <= degree; l++) {
        // build level l of the pyramid
        for (let i = s; i > s - degree - 1 + l; i--) {
            const alpha = (t - knots[i]) / (knots[i + degree + 1 - l] - knots[i]);
            // interpolate each component
            for (let d = 0; d < 4; d++) {
                v[i * 4 + d] = (1 - alpha) * v[(i - 1) * 4 + d] + alpha * v[i * 4 + d];
            }
        }
    }
    return new Vector(v.slice(s * 4, s * 4 + 4));
}

/**
 * In general: the z-dir shadow volume of a face is the integral: SURFACE_INTEGRAL[p in face] (normal(p).z * p.z) dp
 * In general: the centroid of the z-dir shadow volume of a face is the integral:
 *     SURFACE_INTEGRAL[p in face] ((p schur (1, 1, 0.5)) * normal(p).z * p.z) dp
 *     dividing the z component by 2 is usually done at the very end
 */
const ZDirVolumeVisitor = {
    [ConicSurface.name](edges) {
        console.log(this);
        const dpdu = this.dpdu();
        const dpdv = this.dpdv();
        // INT[edge.at; edge.bT] (at(t) DOT dir) * (at(t) - at(t).projectedOn(dir) / 2).z dt
        const totalVolume = edges
            .map((edgeWC) => {
            const curveWC = edgeWC.curve;
            if (curveWC instanceof EllipseCurve ||
                curveWC instanceof HyperbolaCurve ||
                curveWC instanceof ParabolaCurve) {
                const f = (curveT) => {
                    const at = curveWC.at(curveT), tangentWC = curveWC.tangentAt(curveT);
                    const uvOfPWC = this.uvP(at);
                    // INTEGRATE [0; atUV.y] (dpdu(atUV.x, t) X dpdv(atUV.x)).z * pUV(atUV.x, t).z dt
                    // dpdu(u, v) === t * dpdu(s, 1)
                    // => INTEGRATE [0; atUV.y] (t * dpdu(atUV.x, 1) X dpdv(atUV.x)).z * pUV(atUV.x, t).z dt
                    // => (dpdu(atUV.x, 1) X dpdv(atUV.x)).z * INTEGRATE [0; atUV.y] t * pUV(atUV.x, t).z dt
                    // pUV(u, v) === t * (pUV(s, 1) - center) + center
                    // => (dpdu(atUV.x, 1) X dpdv(atUV.x)).z
                    //      * INTEGRATE [0; atUV.y] t² * (pUV(atUV.x, t) - center).z + t * center.z dt
                    // => (dpdu(atUV.x, 1) X dpdv(atUV.x)).z
                    //      * INTEGRATE [0; atUV.y] t² * (pUV(atUV.x, t) - center).z + t * center.z dt
                    // => (dpdu(atUV.x, 1) X dpdv(atUV.x)).z
                    //      * (1/3 t³ pUV(atUV.x, 1).z + 1/2 t² center.z)[0; atUV.y]
                    const du = -M4.forSys(dpdu(uvOfPWC.x, uvOfPWC.y), dpdv(uvOfPWC.x))
                        .inversed()
                        .transformVector(tangentWC).x;
                    const factor = (Math.pow(uvOfPWC.y, 3) / 3) *
                        (this.pUV(uvOfPWC.x, 1).z - this.center.z) +
                        (Math.pow(uvOfPWC.y, 2) / 2) * this.center.z;
                    const actual = dpdu(uvOfPWC.x, factor).cross(dpdv(uvOfPWC.x)).z;
                    return actual * du;
                };
                const val = glqInSteps(f, edgeWC.aT, edgeWC.bT, 1);
                return val;
            }
            else if (curveWC instanceof L3) {
                return 0;
            }
            else {
                throw new Error();
            }
        })
            .sum();
        const centroidZX2Parts = edges.map((edgeWC) => {
            const curveWC = edgeWC.curve;
            if (curveWC instanceof EllipseCurve ||
                curveWC instanceof HyperbolaCurve ||
                curveWC instanceof ParabolaCurve) {
                const f = (curveT) => {
                    const at = curveWC.at(curveT), tangentWC = curveWC.tangentAt(curveT);
                    const uvOfPWC = this.uvP(at);
                    // INTEGRATE [0; atUV.y] dpdu(atUV.x, t) X dpdv(atUV.x, t) * pUV(atUV.x, t).z dt
                    // dpdv is constant with respect to t
                    // => (dpdu(atUV.x, t) X dpdv(atUV.x, t)).z
                    //      * (INTEGRATE [0; atUV.y] t * pUV(atUV.x, t) * pUV(atUV.x, t).z dt)
                    // dpdu(u, v) === t * dpdu(s, 1)
                    // pUV(u, v) === t * (pUV(s, 1) - center) + center
                    // INTEGRATE [0; atUV.y] t * pUV(atUV.x, t) * pUV(atUV.x, t).z dt
                    // = INTEGRATE [0; atUV.y] t *
                    //                         (t * (pUV(s, 1) - center) + center) *
                    //                         (t (pUV(s, 1) - center).z + center.z) dt
                    // = INTEGRATE [0; atUV.y] t³ (pUV(s, 1) - center) * (pUV(s, 1) - center).z
                    //                       + t² ((pUV(s, 1) - center) * center.z + (pUV(s, 1) - center).z * center)
                    //                       + t center center.z dt
                    // = (1/4 t^4 (pUV(s, 1) - center) * (pUV(s, 1) - center).z
                    //   (1/3 t³ ((pUV(s, 1) - center) * center.z + (pUV(s, 1) - center).z * center)
                    //   (1/2 t² center center.z dt)[0; atUV.y]
                    const pUVS1V = this.pUV(uvOfPWC.x, 1).minus(this.center);
                    const factor = V3.add(pUVS1V.times((1 / 4) * Math.pow(uvOfPWC.y, 4) * pUVS1V.z +
                        (1 / 3) * Math.pow(uvOfPWC.y, 3) * this.center.z), this.center.times((1 / 3) * Math.pow(uvOfPWC.y, 3) * pUVS1V.z +
                        (1 / 2) * Math.pow(uvOfPWC.y, 2) * this.center.z));
                    const partialCentroid = factor.times(dpdu(uvOfPWC.x, 1).cross(dpdv(uvOfPWC.x)).z);
                    const ds = -M4.forSys(dpdu(uvOfPWC.x, uvOfPWC.y), dpdv(uvOfPWC.x))
                        .inversed()
                        .transformVector(tangentWC).x;
                    return partialCentroid.times(ds);
                };
                return glqV3(f, edgeWC.aT, edgeWC.bT);
            }
            else if (curveWC instanceof L3) {
                return V3.O;
            }
            else {
                throw new Error();
            }
        });
        const centroid = V3.add(...centroidZX2Parts)
            .schur(new V3(1, 1, 0.5))
            .div(totalVolume);
        return { volume: totalVolume, centroid: centroid };
    },
    [PlaneSurface.name](edges) {
        const r1 = this.right;
        const u1 = this.up;
        const c = this.plane.anchor;
        assert(r1.hasLength(1));
        assert(u1.hasLength(1));
        assert(r1.isPerpendicularTo(u1));
        const volumeAndCentroidZX2Parts = edges.map((edgeWC) => {
            const curveWC = edgeWC.curve;
            if (curveWC instanceof L3) {
                // split shadow volume into two triangle shadow volumes and use the same logic as for mesh triangles:
                function triangleShadowVolumeAndCentroid(a, b, c) {
                    const ab = b.minus(a), ac = c.minus(a);
                    const normal = ab.cross(ac);
                    const faceCentroid = V3.add(a, b, c).div(3);
                    return [
                        (faceCentroid.z * normal.z) / 2,
                        V3.add(a.times(2 * a.z + b.z + c.z), b.times(a.z + 2 * b.z + c.z), c.times(a.z + b.z + 2 * c.z)).times(normal.z),
                    ];
                }
                const a = edgeWC.a, b = edgeWC.b;
                const as = a.dot(r1);
                const bs = b.dot(r1);
                const aBase = this.pUV(as, 0);
                const bBase = this.pUV(bs, 0);
                const [v1, c1] = triangleShadowVolumeAndCentroid(a, b, aBase);
                const [v2, c2] = triangleShadowVolumeAndCentroid(bBase, aBase, b);
                return [v1 + v2, c1.plus(c2).div(24)];
            }
            else if (curveWC instanceof ImplicitCurve) {
                throw new Error();
            }
            else {
                const sliceAreaAndCentroidZX2TimesDs = (curveT) => {
                    const p = curveWC.at(curveT);
                    const s = p.dot(r1);
                    const t = p.dot(u1);
                    const area = t * c.z + s * t * r1.z + (1 / 2) * Math.pow(t, 2) * u1.z;
                    const ds = -curveWC.tangentAt(curveT).dot(r1);
                    return [
                        area * ds,
                        ...V3.add(c.times(area), r1.times(c.z * s * t + r1.z * Math.pow(s, 2) * t + (1 / 2) * s * Math.pow(t, 2) * u1.z), u1.times((1 / 2) * c.z * Math.pow(t, 2) +
                            (1 / 2) * r1.z * s * Math.pow(t, 2) +
                            (1 / 3) * Math.pow(t, 3) * u1.z)).times(ds),
                    ];
                };
                const [vol, cx, cy, cz] = glqArray(sliceAreaAndCentroidZX2TimesDs, edgeWC.aT, edgeWC.bT, 4);
                return [
                    vol * this.plane.normal1.z,
                    new V3(cx, cy, cz).times(this.plane.normal1.z),
                ];
            }
        });
        return mergeVolumeAndCentroidZX2Parts(volumeAndCentroidZX2Parts);
    },
    /**
     * Generic implementation.
     */
    [ParametricSurface.name](edges) {
        const dpdu = this.dpdu();
        const dpdv = this.dpdv();
        const volume = edges.map((edgeWC) => {
            const curveWC = edgeWC.curve;
            if (curveWC instanceof ImplicitCurve) {
                throw new Error();
            }
            else {
                const sliceAreaAndCentroidZX2TimesDs = (curveT) => {
                    // use curve.tangent not edge.tangent, reverse edges are handled by the integration boundaries
                    const pWC = curveWC.at(curveT), tangentWC = curveWC.tangentAt(curveT);
                    const uvOfPWC = this.uvP(pWC);
                    const slice = (t) => {
                        const p = this.pUV(uvOfPWC.x, t);
                        const normal = dpdu(uvOfPWC.x, t).cross(dpdv(uvOfPWC.x, t));
                        return p.z * normal.z;
                    };
                    const sliceIntegral0ToPWCT = glqInSteps(slice, 0, uvOfPWC.y, 1);
                    // const dt = tangentWC.dot(scalingVector)
                    const dt = -M4.forSys(dpdu(uvOfPWC.x, uvOfPWC.y), dpdv(uvOfPWC.x, uvOfPWC.y))
                        .inversed()
                        .transformVector(tangentWC).x;
                    const sliceAreaTimesDs = sliceIntegral0ToPWCT * dt;
                    const slice2 = (t) => {
                        const p = this.pUV(uvOfPWC.x, t);
                        const normal = dpdu(uvOfPWC.x, t).cross(dpdv(uvOfPWC.x, t));
                        return p.times(p.z * normal.z);
                    };
                    const sliceIntegral0ToPWCT2 = glqV3(slice2, 0, uvOfPWC.y);
                    // const dt = tangentWC.dot(scalingVector)
                    const sliceCentroidZX2TimesDs = sliceIntegral0ToPWCT2.times(dt);
                    return [sliceAreaTimesDs, ...sliceCentroidZX2TimesDs.toArray()];
                };
                const [vol, cx, cy, cz] = glqArray(sliceAreaAndCentroidZX2TimesDs, edgeWC.aT, edgeWC.bT, 4);
                return [vol, new V3(cx, cy, cz)];
            }
        });
        return mergeVolumeAndCentroidZX2Parts(volume);
    },
    /**
     * at(t)
     * |\                                    ^
     * | \ at(t).projectedOn(dir1)            \  dir1
     * |  \                                    \
     * |   \ at(t).rejectedFrom(dir1) = b
     * |   |
     * |___|
     *        z = 0
     *
     *
     * A = ((at(t) + at(t).rejectedFrom(dir1)) / 2).z * at(t).projectedOn(dir1).lengthXY()
     * scaling = tangentAt(t) DOT dir1.cross(V3.Z).unit()
     */
    [ProjectedCurveSurface.name](edges) {
        if (V3.Z.cross(this.dir).likeO())
            return { volume: 0, centroid: V3.O };
        // normalize this.dir so it always points up
        const upDir1 = this.dir.toLength(Math.sign(this.dir.z) || 1);
        const scalingVector = V3.Z.cross(upDir1).unit();
        // the length of the base of the trapezoid is calculated by dotting with the baseVector
        const baseVector = upDir1.rejectedFrom(V3.Z).unit();
        // INT[edge.at; edge.bT] (at(t) DOT dir1) * (at(t) - at(t).projectedOn(dir) / 2).z
        const volume = edges.map((edgeWC) => {
            if (edgeWC.curve instanceof L3) {
                return [0, V3.O];
            }
            else if (edgeWC.curve instanceof ImplicitCurve) {
                return [0, V3.O];
                // 	const { points, tangents } = edgeWC.curve
                // 	const minT = edgeWC.minT,
                // 		maxT = edgeWC.maxT
                // 	let sum = 0
                // 	const start = Math.ceil(minT + NLA_PRECISION)
                // 	const end = Math.floor(maxT - NLA_PRECISION)
                // 	for (let i = start; i <= end; i++) {
                // 		const at = points[i],
                // 			tangent = tangents[i]
                // 		const area = (at.z + at.rejectedFrom1(upDir1).z) / 2 * at.projectedOn(upDir1).dot(baseVector)
                // 		const scale = tangent.dot(scalingVector)
                // 		sum += area * scale
                // 	}
                // 	const f = (t: number) => {
                // 		const at = edgeWC.curve.at(t),
                // 			tangent = edgeWC.curve.tangentAt(t)
                // 		const area = (at.z + at.rejectedFrom1(upDir1).z) / 2 * at.projectedOn(upDir1).dot(baseVector)
                // 		const scale = tangent.dot(scalingVector)
                // 		return area * scale
                // 	}
                // 	sum += f(minT) * (start - minT - 0.5)
                // 	sum += f(maxT) * (maxT - end - 0.5)
                // 	return sum * Math.sign(edgeWC.deltaT())
            }
            else {
                const f = (curveT) => {
                    // use curve.tangent not edge.tangent, reverse edges are handled by the integration boundaries
                    const at = edgeWC.curve.at(curveT), tangent = edgeWC.curve.tangentAt(curveT);
                    const b = at.rejectedFrom1(upDir1);
                    const area = (at.z * b.to(at).dot(baseVector)) / 2 +
                        (b.z * b.to(at).dot(baseVector)) / 2;
                    const areaCentroidA = V3.add(at.xy(), b, at).times((at.z * b.to(at).dot(baseVector)) / 2 / 3);
                    const areaCentroidB = V3.add(at.xy(), b, b.xy()).times((b.z * b.to(at).dot(baseVector)) / 2 / 3);
                    const scale = tangent.dot(scalingVector);
                    return [
                        area * scale,
                        ...areaCentroidA.plus(areaCentroidB).times(scale).schur(V(1, 1, 2)),
                    ];
                };
                const [vol, cx, cy, cz] = glqArray(f, edgeWC.aT, edgeWC.bT, 4);
                return [vol, new V3(cx, cy, cz)];
            }
        });
        return mergeVolumeAndCentroidZX2Parts(volume);
    },
    // volume does scale linearly, so this could be done in the local coordinate system
    // however, shear matrices lead to point-to-plane distances having to be calculated along a vector other than
    // the plane normal
    [RotatedCurveSurface.name](edges) {
        const dpdu = this.dpdu();
        const dpdv = this.dpdv();
        const totalVolume = edges
            .map((edgeWC) => {
            const curveWC = edgeWC.curve;
            const f = (curveT) => {
                const pWC = curveWC.at(curveT), tangentWC = curveWC.tangentAt(curveT);
                const uvOfPWC = this.uvP(pWC);
                const pLC = this.matrixInverse.transformPoint(pWC);
                const dpdvAtS0 = this instanceof RotatedCurveSurface
                    ? this.curve.tangentAt(uvOfPWC.y)
                    : V(-pLC.z, 0, pLC.lengthXY());
                // const slice = (phi: number) => {
                // 	const p = this.pUV(phi, uvOfPWC.y)
                // 	const normal = dpdu(phi, uvOfPWC.y).cross(dpdv(phi, uvOfPWC.y))
                // 	return p.z * normal.z
                // }
                // const z = this.curve.at(uvOfPWC.y).z
                // const r = this.curve.at(uvOfPWC.y).lengthXY()
                // const pz =
                // 	this.f1.z * r * cos(s) +
                // 	this.f2.z * r * sin(s) +
                // 	this.f3.z * z +
                // 	this.center.z
                // const dpdux = this.f1.x * r * -sin(s) + this.f2.x * r * cos(s)
                // const dpduy = this.f1.y * r * -sin(s) + this.f2.y * r * cos(s)
                // const dpdvx = this.f1.x * dr * cos(s) + this.f2.x * dr * sin(s) + this.f3.x * dz
                // const dpdvy = this.f1.y * dr * cos(s) + this.f2.y * dr * sin(s) + this.f3.y * dz
                // const normalz = dpdux * dpdvy - dpduy * dpdvx
                // result = pz * normalz
                const r = pLC.lengthXY(), z = pLC.z;
                const dr = dpdvAtS0.x;
                const dz = dpdvAtS0.z;
                const a = this.matrix.X.z * r, b = this.matrix.Y.z * r, c = this.matrix.Z.z * z + this.matrix.O.z;
                const t0 = (this.matrix.X.x * this.matrix.Y.y -
                    this.matrix.X.y * this.matrix.Y.x) *
                    r *
                    dr;
                const t1 = (this.matrix.Y.x * this.matrix.X.y -
                    this.matrix.Y.y * this.matrix.X.x) *
                    r *
                    dr;
                const t2 = (this.matrix.X.x * this.matrix.X.y -
                    this.matrix.X.y * this.matrix.X.x) *
                    r *
                    dr;
                const t3 = (this.matrix.Y.x * this.matrix.Y.y -
                    this.matrix.Y.y * this.matrix.Y.x) *
                    r *
                    dr;
                const t4 = (this.matrix.Y.x * this.matrix.Z.y -
                    this.matrix.Y.y * this.matrix.Z.x) *
                    r *
                    dz;
                const t5 = (this.matrix.X.x * this.matrix.Z.y -
                    this.matrix.X.y * this.matrix.Z.x) *
                    r *
                    dz;
                const sliceIntegral = (p) => {
                    return ((6 * (c * (-t0 + t1) + a * t4 - b * t5) * p +
                        3 *
                            (3 * b * t0 - b * t1 + a * (t2 - t3) + 4 * c * t5) *
                            cos(p) +
                        3 *
                            (3 * a * t1 - a * t0 - b * (t2 - t3) + 4 * c * t4) *
                            sin(p) +
                        3 * (a * t5 - b * t4 + c * (t2 - t3)) * cos(2 * p) +
                        3 * (a * t4 + b * t5 + c * (t0 + t1)) * sin(2 * p) +
                        (a * (t2 - t3) - b * (t0 + t1)) * cos(3 * p) +
                        (a * (t0 + t1) + b * (t2 - t3)) * sin(3 * p)) /
                        12);
                };
                const dt = M4.forSys(dpdu(uvOfPWC.x, uvOfPWC.y), dpdv(uvOfPWC.x, uvOfPWC.y))
                    .inversed()
                    .transformVector(tangentWC).y;
                const sliceIntegral0ToPWCS = sliceIntegral(uvOfPWC.x); //- sliceIntegral(0) //(always 0)
                const result = sliceIntegral0ToPWCS * dt;
                return result;
            };
            return gaussLegendreQuadrature24(f, edgeWC.aT, edgeWC.bT);
        })
            .sum();
        // calc centroid:
        const centroidZX2Parts = edges.map((edgeWC) => {
            const f = (curveT) => {
                const curveWC = edgeWC.curve;
                const pWC = curveWC.at(curveT), tangentWC = curveWC.tangentAt(curveT);
                const uvOfPWC = this.uvP(pWC);
                const slice = (phi) => {
                    const p = this.pUV(phi, uvOfPWC.y);
                    const normal = dpdu(phi, uvOfPWC.y).cross(dpdv(phi, uvOfPWC.y));
                    return p.times(p.z * normal.z);
                };
                const sliceIntegral0ToPWCS = glqV3(slice, 0, uvOfPWC.x);
                const dt = M4.forSys(dpdu(uvOfPWC.x, uvOfPWC.y), dpdv(uvOfPWC.x, uvOfPWC.y))
                    .inversed()
                    .transformVector(tangentWC).y;
                const result = sliceIntegral0ToPWCS.times(dt);
                return result;
            };
            return glqV3(f, edgeWC.aT, edgeWC.bT);
        });
        const centroid = V3.add(...centroidZX2Parts)
            .schur(new V3(1, 1, 0.5))
            .div(totalVolume);
        return { volume: totalVolume, centroid: centroid };
    },
};
ZDirVolumeVisitor[EllipsoidSurface.name] =
    ZDirVolumeVisitor[RotatedCurveSurface.name];
function glqV3(f, startT, endT) {
    return gaussLegendre24Xs
        .reduce((val, currVal, index) => {
        const x = startT + ((currVal + 1) / 2) * (endT - startT);
        return val.plus(f(x).times(gaussLegendre24Weights[index]));
    }, V3.O)
        .times((endT - startT) / 2);
}
function glqArray(f, startT, endT, numEls = 3) {
    const result = new Array(numEls).fill(0);
    for (let i = 0; i < 24; i++) {
        const x = startT + ((gaussLegendre24Xs[i] + 1) / 2) * (endT - startT);
        const fx = f(x);
        for (let j = 0; j < numEls; j++) {
            result[j] += fx[j] * gaussLegendre24Weights[i];
        }
    }
    for (let j = 0; j < numEls; j++) {
        result[j] *= (endT - startT) / 2;
    }
    return result;
}
function mergeVolumeAndCentroidZX2Parts(volumeAndCentroidZX2Parts) {
    const volume = volumeAndCentroidZX2Parts.reduce((result, [volume]) => result + volume, 0);
    const weightedCentroid = V3.add(...volumeAndCentroidZX2Parts.map(([, centroidZX2]) => centroidZX2)).schur(new V3(1, 1, 0.5));
    return { volume, centroid: weightedCentroid.div(volume) };
}

const CalculateAreaVisitor = {
    [ConicSurface.name](edges) {
        const dpdu = this.dpdu();
        const dpdv = this.dpdv();
        // calculation cannot be done in local coordinate system, as the area doesnt scale proportionally
        const totalArea = edges
            .map((edge) => {
            if (edge.curve instanceof EllipseCurve ||
                edge.curve instanceof HyperbolaCurve ||
                edge.curve instanceof ParabolaCurve) {
                const f = (t) => {
                    const at = edge.curve.at(t), tangentWC = edge.tangentAt(t);
                    const uvOfPWC = this.uvP(at);
                    // INTEGRATE [0; atUV.y]
                    //   dpdu(atUV.x, t) X dpdv(atUV.x, t)
                    // dt
                    // dpdv is constant with respect to t
                    // => dpdv(atUV.x, 0) X (INTEGRATE [0; atUV.y] dpdu(atUV.x, t) dt)
                    // dpdu(u, v) === v * dpdu(u, 1)
                    // => dpdv(atUV.x, 0) X (1/2 t² dpdu(atUV.x, 1))[0; atUV.y]
                    // => dpdv(atUV.x, 0) X dpdu(atUV.x, atUV.y² / 2)
                    const du = -M4.forSys(dpdu(uvOfPWC.x, uvOfPWC.y), dpdv(uvOfPWC.x))
                        .inversed()
                        .transformVector(tangentWC).x;
                    return (dpdu(uvOfPWC.x, Math.pow(uvOfPWC.y, 2) / 2)
                        .cross(dpdv(uvOfPWC.x))
                        .length() * du);
                };
                return glqInSteps(f, edge.aT, edge.bT, 1);
            }
            else if (edge.curve instanceof L3) {
                return 0;
            }
            else {
                throw new Error();
            }
        })
            .sum();
        return totalArea * this.normalDir;
    },
    [PlaneSurface.name](edges) {
        let totalArea = 0;
        const r1 = this.right, u1 = this.up;
        for (const edge of edges) {
            let edgeArea;
            const curve = edge.curve;
            if (curve instanceof L3) {
                edgeArea =
                    ((edge.a.dot(u1) + edge.b.dot(u1)) / 2) * edge.b.to(edge.a).dot(r1);
            }
            else if (curve instanceof EllipseCurve) {
                // INTEGRATE[aT; bT] (curve.at(t) * u1) * (tangent(t) * r1) dt
                // INTEGRATE[aT; bT] (u1 f1 cos t + u1 f2 sin t + u1 c) * (r1 f1 (-sin t) + r1 f2 cos t) dt
                const { f1, f2, center } = curve;
                const a = u1.dot(f1), b = u1.dot(f2), c = u1.dot(center), d = r1.dot(f1), e = r1.dot(f2);
                function fArea(t) {
                    return (0.25 *
                        (2 * (-b * d + a * e) * t +
                            4 * c * d * cos(t) +
                            4 * c * e * sin(t) +
                            (a * d - b * e) * cos(2 * t) +
                            (b * d + a * e) * sin(2 * t)));
                }
                edgeArea = -(fArea(edge.bT) - fArea(edge.aT));
            }
            else if (curve instanceof ImplicitCurve) {
                throw new Error("implement for implicitCurve");
            }
            else {
                const dir1 = u1;
                assertf(() => dir1.hasLength(1));
                // INT[aT; bT] at(t) * dir1 * tangentAt(t).rejectedFrom(dir1) dt
                const f = (curveT) => {
                    const at = curve.at(curveT);
                    const tangent = curve.tangentAt(curveT);
                    const ds = r1.dot(tangent);
                    const t = u1.dot(at);
                    return ds * t;
                };
                edgeArea = glqInSteps(f, edge.aT, edge.bT, 3);
            }
            totalArea += edgeArea;
        }
        assert(isFinite(totalArea));
        return totalArea;
    },
    [RotatedCurveSurface.name](edges, canApproximate = true) {
        const f1 = this.matrix.X, f2 = this.matrix.Y, f3 = this.matrix.Z;
        const likeVerticalSpheroid = eq(f1.length(), f2.length()) &&
            f1.isPerpendicularTo(f2) &&
            f2.isPerpendicularTo(f3) &&
            f3.isPerpendicularTo(f1);
        const areaParts = edges.map((edgeWC, ei) => {
            console.log("edge", ei, edgeWC.sce);
            const curveWC = edgeWC.curve;
            if (edgeWC.curve instanceof ImplicitCurve) {
                throw new Error();
            }
            else {
                if (likeVerticalSpheroid) {
                    const f = (curveT) => {
                        const pWC = curveWC.at(curveT), tangent = curveWC.tangentAt(curveT);
                        const pLC = this.matrixInverse.transformPoint(pWC);
                        const { x: angleXY, y: t } = this.uvP(pWC);
                        const arcRadius = this.matrix.transformVector(pLC.xy()).length();
                        const arcLength = angleXY * arcRadius;
                        const dpdv = this.dpdv()(angleXY, t).unit();
                        const scaling = dpdv.dot(tangent);
                        return arcLength * scaling;
                    };
                    return glqInSteps(f, edgeWC.aT, edgeWC.bT, 1);
                }
                else {
                    const dpdu = this.dpdu(), dpdv = this.dpdv();
                    const f2 = (curveT) => {
                        const pWC = curveWC.at(curveT), tangentWC = curveWC.tangentAt(curveT);
                        const uvPWC = this.uvP(pWC);
                        const slice = (phi) => {
                            //return this.dpdu()(phi, st.y).length() * this.dpdv()(phi, st.y).length()
                            return dpdu(phi, uvPWC.y).cross(dpdv(phi, uvPWC.y)).length();
                        };
                        // we need to do a coordinate transform from curveT to dt, as that is what we are integrating
                        const dt = M4.forSys(dpdu(uvPWC.x, uvPWC.y), dpdv(uvPWC.x, uvPWC.y))
                            .inversed()
                            .transformVector(tangentWC).y;
                        return glqInSteps(slice, 0, uvPWC.x, 1) * dt;
                    };
                    return glqInSteps(f2, edgeWC.aT, edgeWC.bT, 1);
                }
            }
        });
        return areaParts.sum();
    },
    [ProjectedCurveSurface.name](edges) {
        // calculation cannot be done in local coordinate system, as the area doesn't scale proportionally
        const thisDir1 = this.dir.unit();
        const totalArea = edges
            .map((edge) => {
            if (edge.curve instanceof L3) {
                return 0;
            }
            else if (edge.curve instanceof ImplicitCurve) {
                const { points, tangents } = edge.curve;
                const minT = edge.minT, maxT = edge.maxT;
                let sum = 0;
                const start = ceil(minT + NLA_PRECISION);
                const end = floor(maxT - NLA_PRECISION);
                for (let i = start; i <= end; i++) {
                    const at = points[i], tangent = tangents[i]; //.toLength(edge.curve.stepSize)
                    const scaling = this.normalP(at).cross(thisDir1).unit().dot(tangent);
                    sum += at.dot(thisDir1) * scaling;
                }
                const f = (t) => {
                    const at = edge.curve.at(t), tangent = edge.curve.tangentAt(t);
                    const scaling = this.normalP(at).cross(thisDir1).unit().dot(tangent);
                    return at.dot(thisDir1) * scaling;
                };
                sum += f(minT) * (start - minT - 0.5);
                sum += f(maxT) * (maxT - end - 0.5);
                return sum * sign(edge.deltaT());
            }
            else {
                const f = (t) => {
                    const at = edge.curve.at(t);
                    const tangent = edge.tangentAt(t);
                    const scaling = tangent.rejected1Length(thisDir1);
                    return at.dot(thisDir1) * scaling;
                };
                const val = glqInSteps(f, edge.aT, edge.bT, 1);
                const sign = Math.sign(this.normalP(edge.a)
                    .cross(this.dir)
                    .dot(edge.curve.tangentAt(edge.aT)));
                assert(0 !== sign);
                return val * sign;
            }
        })
            .sum();
        console.log("totalArea", totalArea);
        return totalArea;
    },
};
CalculateAreaVisitor[EllipsoidSurface.name] =
    CalculateAreaVisitor[RotatedCurveSurface.name];

/**
 * Create a surface by projecting a curve in a direction.
 *
 * @param curve The curve to project.
 * @param offset The direction and distance to project curve.
 * @param flipped Whether the surface's default orientation (normal = curve tangent cross offset) should be flipped.
 */
function projectCurve(curve, offset, flipped) {
    if (curve instanceof L3) {
        const surfaceNormal = offset.cross(curve.dir1).toLength(flipped ? -1 : 1);
        return new PlaneSurface(P3.normalOnAnchor(surfaceNormal, curve.anchor));
    }
    if (curve instanceof EllipseCurve) {
        const curveDir = flipped ? offset : offset.negated();
        return new CylinderSurface(curve, curveDir.unit(), undefined, undefined);
    }
    if (curve instanceof BezierCurve || curve instanceof XiEtaCurve) {
        const curveDir = offset.times(flipped ? 1 : -1);
        return new ProjectedCurveSurface(curve, curveDir, undefined, undefined, flipped ? 0 : -1, flipped ? 1 : 0);
    }
    throw new Error();
}
/**
 * Create a surface by projecting a curve onto a point.
 */
function projectPointCurve(curve, tMin = curve.tMin, tMax = curve.tMax, p, flipped) {
    if (curve instanceof L3) {
        const up = curve.anchor.to(p).rejectedFrom(curve.dir1);
        return PlaneSurface.forAnchorAndPlaneVectors(curve.anchor, curve.dir1, up.unit(), tMin, tMax, 0, up.length());
    }
    else if (curve instanceof EllipseCurve) {
        // flip f2 by default
        const factor = -1 * (flipped ? -1 : 1);
        return new ConicSurface(p, curve.f1.times(factor), curve.f2, p.to(curve.center), tMin, tMax, 0, 1);
    }
    else {
        throw new Error("projectPointCurve not implemented for " + curve.constructor.name);
    }
}
/**
 * Create a surface by rotating a curve in the XZ-plane, with X > 0, around the Z-axis according to the right-hand rule.
 * @param curve The curve to rotate.
 * @param tMin The minimum value for t for which the surface should be defined.
 * @param tMax The maximum value for t for which the surface should be defined.
 * @param angle How much the curve should be rotated. sMin/sMax will be be 0/angle.
 * @param flipped Whether the surface's default orientation (normal = curve tangent cross rotation tangent) should be
 * flipped.
 */
function rotateCurve(curve, tMin = curve.tMin, tMax = curve.tMax, angle, flipped) {
    assertf(() => new PlaneSurface(P3.ZX).containsCurve(curve));
    if (curve instanceof L3) {
        if (curve.dir1.isParallelTo(V3.Z)) {
            if (eq0(curve.anchor.x)) {
                return undefined;
                //throw new Error('Cannot rotate curve colinear to Z axis.')
            }
            const baseEllipse = new EllipseCurve(V3.O, curve.anchor.xy(), curve.anchor.xy().getPerpendicular(), 0, angle);
            // if curve.dir1 is going up (+Z), it the cylinder surface should face inwards
            const factor = (curve.dir1.z > 0 ? -1 : 1) * (flipped ? -1 : 1);
            const [zMin, zMax] = [
                curve.at(tMin).z * factor,
                curve.at(tMax).z * factor,
            ].sort(MINUS);
            return new CylinderSurface(baseEllipse, V3.Z.times(factor), 0, angle, zMin, zMax);
        }
        if (curve.at(tMin).xy().dot(curve.dir1) *
            curve.at(tMax).xy().dot(curve.dir1) <
            0) {
            throw new Error("line cannot cross the Z axis in the [tMin, tMax] interval, as conic surfaces cannot have an hourglass shape.");
        }
        if (curve.dir1.isPerpendicularTo(V3.Z)) {
            // if line.dir1 is pointing aways from V3.Z, then the surface should face up
            const factor = (curve.at(lerp(tMin, tMax, 0.5)).dot(curve.dir1) > 0 ? 1 : -1) *
                (flipped ? -1 : 1);
            return new PlaneSurface(new P3(V3.Z.times(factor), curve.anchor.z * factor));
        }
        else {
            // apex is intersection of segment with Z-axis
            const a = curve.at(tMin), b = curve.at(tMax);
            const apexZ = a.z - (a.x * (b.z - a.z)) / (b.x - a.x);
            const apex = new V3(0, 0, apexZ);
            const factor = -(a.x > b.x ? -1 : 1) * (flipped ? -1 : 1);
            const s = new ConicSurface(apex, new V3(curve.dir1.lengthXY(), 0, 0), new V3(0, curve.dir1.lengthXY(), 0), new V3(0, 0, (a.x > b.x ? -1 : 1) * curve.dir1.z), 0, angle, 0, 1);
            return factor > 0 ? s : s.flipped();
        }
    }
    if (curve instanceof EllipseCurve) {
        const a = curve.at(tMin), b = curve.at(tMax);
        const ell = curve.rightAngled();
        const f1Perp = ell.f1.isPerpendicularTo(V3.Z), f2Perp = ell.f2.isPerpendicularTo(V3.Z);
        if (L3.Z.containsPoint(ell.center) && (f1Perp || f2Perp)) {
            flipped = flipped == a.z > b.z;
            let width = ell.f1.length(), height = ell.f2.length();
            if (ell.f1.isParallelTo(V3.Z)) {
                [width, height] = [height, width];
            }
            return EllipsoidSurface.forABC(width, (!flipped ? 1 : -1) * width, height, ell.center);
        }
        else {
            const s = new RotatedCurveSurface(curve, M4.IDENTITY, tMin, tMax);
            return s;
        }
    }
    throw new Error();
}
var B2T;
(function (B2T) {
    /**
     * Create a [BRep] of an axis-aligned box width starting at the origin and extending into +XYZ space.
     * @param width x-direction size.
     * @param height y-direction size.
     * @param depth z-direction size.
     * @param name
     */
    function box(width = 1, height = 1, depth = 1, name = "box" + getGlobalId()) {
        assertNumbers(width, height, depth);
        assert("string" === typeof name);
        const baseVertices = [
            new V3(0, 0, 0),
            new V3(0, height, 0),
            new V3(width, height, 0),
            new V3(width, 0, 0),
        ];
        const generator = callsce("B2T.box", width, height, depth, name);
        return B2T.extrudeVertices(baseVertices, P3.XY.flipped(), new V3(0, 0, depth), name, generator);
    }
    B2T.box = box;
    function puckman(radius, rads, height, name = "puckman" + getGlobalId()) {
        assertf(() => lt(0, radius));
        assertf(() => lt(0, rads) && le(rads, TAU));
        assertf(() => lt(0, height));
        const edges = StraightEdge.chain([
            V3.O,
            new V3(radius, 0, 0),
            new V3(radius, 0, height),
            new V3(0, 0, height),
        ], true);
        return B2T.rotateEdges(edges, rads, name);
    }
    B2T.puckman = puckman;
    function registerVertexName(map, name, p) {
        // TODO
        if (!Array.from(map.keys()).some((p2) => p2.like(p))) {
            map.set(p, name);
        }
    }
    B2T.registerVertexName = registerVertexName;
    /**
     * Create a [BRep] by projecting a number of edges in a direction.
     * @param baseFaceEdges
     * @param baseFacePlane
     * @param offset
     * @param name
     * @param gen
     * @param infoFactory
     */
    function extrudeEdges(baseFaceEdges, baseFacePlane = P3.XY, offset = V3.Z, name = "extrude" + getGlobalId(), gen, infoFactory) {
        baseFaceEdges = fixEdges(baseFaceEdges);
        //Array.from(combinations(baseFaceEdges.length)).forEach(({i, j}) => {
        //	assertf(() => !Edge.edgesIntersect(baseFaceEdges[i], baseFaceEdges[j]), baseFaceEdges[i].sce +
        // baseFaceEdges[j].sce) })
        assertf(() => Edge.isLoop(baseFaceEdges));
        // TODO checks..
        //if (offset.dot(baseFacePlane.normal1) > 0) {
        //	baseFacePlane = baseFacePlane.flipped()
        //}
        const vertexNames = new Map();
        const basePlaneSurface = new PlaneSurface(baseFacePlane);
        //assert(basePlaneSurface.edgeLoopCCW(baseFaceEdges), 'edges not CCW on baseFacePlane')
        const translationMatrix = M4.translate(offset);
        const topEdges = baseFaceEdges.map((edge) => edge.transform(translationMatrix, "top"));
        const edgeCount = baseFaceEdges.length;
        const bottomInfo = infoFactory && infoFactory.extrudeBottom(basePlaneSurface, baseFaceEdges);
        const bottomFace = new PlaneFace(basePlaneSurface, baseFaceEdges, [], name + "Bottom", bottomInfo);
        const topFaceEdges = topEdges.map((edge) => edge.flipped()).reverse();
        const topSurface = new PlaneSurface(baseFacePlane.flipped().translated(offset));
        const topInfo = infoFactory && infoFactory.extrudeBottom(topSurface, topFaceEdges);
        const topFace = new PlaneFace(topSurface, topFaceEdges, [], name + "Top", topInfo);
        baseFaceEdges.forEach((edge) => B2T.registerVertexName(vertexNames, edge.name + "A", edge.a));
        topFaceEdges.forEach((edge) => B2T.registerVertexName(vertexNames, edge.name + "A", edge.a));
        const ribs = arrayFromFunction(edgeCount, (i) => StraightEdge.throughPoints(baseFaceEdges[i].a, topEdges[i].a, name + "Rib" + i));
        const faces = baseFaceEdges.map((edge, i) => {
            const faceName = name + "Wall" + i;
            const j = (i + 1) % edgeCount;
            const faceEdges = [
                baseFaceEdges[i].flipped(),
                ribs[i],
                topEdges[i],
                ribs[j].flipped(),
            ];
            const surface = projectCurve(edge.curve, offset, edge.reversed);
            const info = infoFactory && infoFactory.extrudeWall(i, surface, faceEdges);
            return Face.create(surface, faceEdges, undefined, faceName, info);
        });
        faces.push(bottomFace, topFace);
        gen =
            gen ||
                callsce("B2T.extrudeEdges", baseFaceEdges, baseFacePlane, offset, name);
        return new BRep(faces, baseFacePlane.normal1.dot(offset) > 0, gen, vertexNames);
    }
    B2T.extrudeEdges = extrudeEdges;
    function cylinder(radius = 1, height = 1, rads = TAU, name = "cylinder" + getGlobalId()) {
        const vertices = [
            new V3(0, 0, 0),
            new V3(radius, 0, 0),
            new V3(radius, 0, height),
            new V3(0, 0, height),
        ];
        return rotateEdges(StraightEdge.chain(vertices, true), rads, name);
    }
    B2T.cylinder = cylinder;
    function cone(radius = 1, height = 1, rads = TAU, name = "cone" + getGlobalId()) {
        const vertices = [
            new V3(0, 0, 0),
            new V3(radius, 0, height),
            new V3(0, 0, height),
        ];
        return rotateEdges(StraightEdge.chain(vertices, true), rads, name);
    }
    B2T.cone = cone;
    function sphere(radius = 1, name = "sphere" + getGlobalId(), rot = TAU) {
        const ee = new PCurveEdge(new EllipseCurve(V3.O, new V3(0, 0, -radius), new V3(radius, 0, 0)), new V3(0, 0, -radius), new V3(0, 0, radius), 0, PI, undefined, new V3(radius, 0, 0), new V3(-radius, 0, 0));
        const generator = callsce("B2T.sphere", radius, name, rot);
        return rotateEdges([StraightEdge.throughPoints(ee.b, ee.a), ee], rot, name, generator);
    }
    B2T.sphere = sphere;
    /**
     * Create a [[BRep]] of a menger sponge.
     * @param res 0: just a cube, 1: every cube face has one hole, 2: 9 holes, etc
     * @param name
     */
    function menger(res = 2, name = "menger" + getGlobalId()) {
        let result = B2T.box(1, 1, 1);
        if (0 == res)
            return result;
        const punch = B2T.box(1 / 3, 1 / 3, 2)
            .translate(1 / 3, 1 / 3, -1 / 2)
            .flipped();
        function recurse(steps, m4) {
            result = result.and(punch.transform(m4));
            if (steps > 1) {
                const scaled = m4.times(M4.scale(1 / 3, 1 / 3, 1));
                for (let i = 0; i < 9; i++) {
                    if (4 == i)
                        continue;
                    recurse(steps - 1, scaled.times(M4.translate(i % 3, (i / 3) | 0, 0)));
                }
            }
        }
        recurse(res, M4.IDENTITY);
        recurse(res, M4.YZX);
        recurse(res, M4.ZXY);
        return result;
    }
    B2T.menger = menger;
    function menger2(res = 2, name = "menger" + getGlobalId()) {
        if (0 == res)
            return B2T.box(1, 1, 1);
        const punch = B2T.box(1 / 3, 1 / 3, 2)
            .translate(1 / 3, 1 / 3, -1 / 2)
            .flipped();
        const stencilFaces = [];
        function recurse(steps, m4) {
            stencilFaces.push(...punch.transform(m4).faces);
            if (steps > 1) {
                const scaled = m4.times(M4.scale(1 / 3, 1 / 3, 1));
                for (let i = 0; i < 9; i++) {
                    if (4 == i)
                        continue;
                    recurse(steps - 1, scaled.times(M4.translate(i % 3, (i / 3) | 0, 0)));
                }
            }
        }
        recurse(res, M4.IDENTITY);
        const stencil = new BRep(stencilFaces, true);
        return B2T.box(1, 1, 1, name)
            .and(stencil)
            .and(stencil.transform(M4.YZX))
            .and(stencil.transform(M4.ZXY));
    }
    B2T.menger2 = menger2;
    /**
     * Create a [BRep] of a torus.
     * @param rSmall The radius to the surface of the torus.
     * @param rLarge The radius from the origin to the inside of the torus.
     * @param rads
     * @param name
     */
    function torus(rSmall, rLarge, rads = TAU, name = "torus" + getGlobalId()) {
        assertNumbers(rSmall, rLarge, rads);
        assertf(() => rLarge > rSmall);
        const curves = [
            EllipseCurve.semicircle(rSmall, new V3(rLarge, 0, 0)),
            EllipseCurve.semicircle(-rSmall, new V3(rLarge, 0, 0)),
        ];
        const baseEdges = curves.map((c) => PCurveEdge.forCurveAndTs(c, 0, Math.PI).rotateX(PI / 2));
        return B2T.rotateEdges(baseEdges, rads, name);
    }
    B2T.torus = torus;
    /**
     * Create a [BRep] by smoothly rotating edges around Z.
     * baseLoop should be CCW on XZ plane for a bounded BRep
     */
    function rotateEdges(baseLoop, totalRads, name = "rotateEdges" + getGlobalId(), generator, infoFactory) {
        assert(baseLoop.every((e) => new PlaneSurface(P3.ZX).containsCurve(e.curve)));
        assert(!eq(PI, totalRads) || PI == totalRads); // URHGJ
        assertf(() => lt(0, totalRads) && le(totalRads, TAU));
        totalRads = snap(totalRads, TAU);
        assertf(() => Edge.isLoop(baseLoop));
        const basePlane = new PlaneSurface(P3.ZX.flipped()).edgeLoopCCW(baseLoop)
            ? new PlaneSurface(P3.ZX.flipped())
            : new PlaneSurface(P3.ZX);
        // const rotationSteps = ceil((totalRads - NLA_PRECISION) / PI)
        // const angles = rotationSteps == 1 ? [-PI, -PI + totalRads] : [-PI, 0, totalRads - PI]
        const open = !eq(totalRads, 2 * PI);
        const baseRibCurves = baseLoop.map((edge) => {
            const a = edge.a, radius = a.lengthXY();
            if (!eq0(radius)) {
                return new EllipseCurve(V(0, 0, a.z), V(radius, 0, 0), V(0, radius, 0));
            }
            return undefined;
        });
        const baseSurfaces = baseLoop.map((edge) => {
            const s = rotateCurve(edge.curve, edge.minT, edge.maxT, PI, edge.deltaT() > 0);
            const t = lerp(edge.aT, edge.bT, 0.5);
            s &&
                assert(edge
                    .tangentAt(t)
                    .cross(V3.Y)
                    .dot(s.normalP(edge.curve.at(t))) < 0);
            return s;
        });
        let stepStartEdges = baseLoop, stepEndEdges;
        const faces = [];
        for (let rot = 0; rot < totalRads; rot += PI) {
            const aT = 0, bT = min(totalRads - rot, PI);
            const rotation = M4.rotateZ(rot + bT);
            stepEndEdges =
                rot + bT == TAU
                    ? baseLoop
                    : baseLoop.map((edge) => edge.transform(rotation));
            const ribs = arrayFromFunction(baseLoop.length, (i) => {
                const a = stepStartEdges[i].a, radius = a.lengthXY();
                const b = stepEndEdges[i].a;
                if (!eq0(radius)) {
                    const curve = 0 === rot ? baseRibCurves[i] : baseRibCurves[i].rotateZ(rot);
                    return new PCurveEdge(curve, a, b, aT, bT, undefined, curve.tangentAt(aT), curve.tangentAt(bT), name + "rib" + i);
                }
                return undefined;
            });
            for (let edgeIndex = 0; edgeIndex < baseLoop.length; edgeIndex++) {
                if (baseSurfaces[edgeIndex]) {
                    const edge = stepStartEdges[edgeIndex];
                    const ipp = (edgeIndex + 1) % baseLoop.length;
                    const faceEdges = [
                        stepStartEdges[edgeIndex].flipped(),
                        !eq0(edge.a.x) && ribs[edgeIndex],
                        stepEndEdges[edgeIndex],
                        !eq0(edge.b.x) && ribs[ipp].flipped(),
                    ].filter((x) => x);
                    const surface = 0 === rot
                        ? baseSurfaces[edgeIndex]
                        : baseSurfaces[edgeIndex].rotateZ(rot);
                    const info = infoFactory &&
                        infoFactory.extrudeWall(edgeIndex, surface, faceEdges, undefined);
                    faces.push(Face.create(surface, faceEdges, undefined, name + "Wall" + edgeIndex, info));
                }
            }
            stepStartEdges = stepEndEdges;
        }
        if (open) {
            const endFaceEdges = Edge.reversePath(stepEndEdges);
            const infoStart = infoFactory && infoFactory.rotationStart(basePlane, baseLoop, undefined);
            const infoEnd = infoFactory &&
                infoFactory.rotationEnd(basePlane.flipped().rotateZ(totalRads), endFaceEdges, undefined);
            faces.push(new PlaneFace(basePlane, baseLoop, undefined, name + "start", infoStart), new PlaneFace(basePlane.flipped().rotateZ(totalRads), endFaceEdges, undefined, name + "end", infoEnd));
        }
        const infiniteVolume = new PlaneSurface(P3.ZX).edgeLoopCCW(baseLoop);
        return new BRep(faces, infiniteVolume, generator);
    }
    B2T.rotateEdges = rotateEdges;
    /**
     * loop should be CCW on XZ plane for a bounded BRep
     */
    //export function rotateEdgesUnsplit(loop: Edge[], rads: raddd, name: string): BRep {
    //	assert(Edge.isLoop(loop))
    //	const rotationMatrix = M4.rotateZ(rads)
    //	const open = !eq(rads, 2 * PI)
    //	const endEdges = open ? loop.map(edge => edge.transform(rotationMatrix)) : loop
    //	const edgeCount = loop.length
    //	const ribs = arrayFromFunction(edgeCount, i => {
    //		const a = loop[i].a, radius = a.lengthXY()
    //		const b = endEdges[i].a
    //		if (!eq0(radius)) {
    //			const curve = new EllipseCurve(V(0, 0, a.z), V(-radius, 0, 0), V(0, -radius, 0))
    //			const aT = -PI, bT = -PI + rads
    //			return new PCurveEdge(curve, a, b, aT, bT, undefined, curve.tangentAt(aT), curve.tangentAt(bT), name
    // + 'rib' + i) } }) const faces = loop.map((edge, i) => { const ipp = (i + 1) % edgeCount console.log('ljl', i,
    // ipp, ribs) const faceEdges = [ edge.flipped(), !eq0(edge.a.x) && ribs[i], endEdges[i], !eq0(edge.b.x) &&
    // ribs[ipp].flipped()].filter(x => x) if (edge instanceof StraightEdge) { const line = edge.curve let surface if
    // (line.dir1.isParallelTo(V3.Z)) { if (eq0(edge.a.x)) { return } let flipped = edge.a.z > edge.b.z surface = new
    // CylinderSurface(ribs[i].curve, !flipped ? V3.Z : V3.Z.negated()) } else if
    // (line.dir1.isPerpendicularTo(V3.Z)) { let flipped = edge.a.x > edge.b.x let surface = new PlaneSurface(new
    // P3(V3.Z, edge.a.z)) if (!flipped) surface = surface.flipped() if (!open) { const hole = flipped ? !eq0(edge.b.x)
    // && ribs[ipp].flipped() : !eq0(edge.a.x) && ribs[i] return new PlaneFace(surface, [flipped ? ribs[i] :
    // ribs[ipp].flipped()], hole && [[hole]]) } return new PlaneFace(surface, faceEdges) } else { // apex is
    // intersection of segment with Z-axis let a = edge.a, b = edge.b let apexZ = a.z - a.x * (b.z - a.z) / (b.x - a.x)
    // let apex = new V3(0, 0, apexZ) let flipped = edge.a.z > edge.b.z surface =
    // ConicSurface.atApexThroughEllipse(apex, ribs[a.x > b.x ? i : ipp].curve as EllipseCurve, !flipped ? 1 : -1)
    // } return Face.create(surface, faceEdges) } if (edge.curve instanceof EllipseCurve) { let flipped = undefined
    // let ell = edge.curve.rightAngled() let f1Perp = ell.f1.isPerpendicularTo(V3.Z), f2Perp =
    // ell.f2.isPerpendicularTo(V3.Z) if (L3.Z.containsPoint(ell.center) && (f1Perp || f2Perp)) { let f3length = f1Perp
    // ? ell.f1.length() : ell.f2.length() if (flipped) { f3length *= -1 } let surface = new
    // EllipsoidSurface(ell.center, ell.f1, ell.f2, ell.f1.cross(ell.f2).toLength(f3length)) return new
    // RotationFace(surface, faceEdges) } } else { assert(false, edge) } }).filter(x => x) if (open) { const
    // endFaceEdges = endEdges.map(edge => edge.flipped()).reverse() faces.push( new PlaneFace(new
    // PlaneSurface(P3.ZX.flipped()), loop), new PlaneFace(new PlaneSurface(P3.ZX.rotateZ(rads)), endFaceEdges)) }
    // return new BRep(faces, undefined) }
    function quaffle() {
        const baseK = B2T.sphere(1).translate(0, 1.7).flipped();
        //const baseK = B2T.box().scale(0.2).translate(0, 0.95).flipped()
        // const vs = B2T.DODECAHEDRON_VERTICES.concat(
        // B2T.DODECAHEDRON_FACE_VERTICES.map(fis => fis
        // .map(vi => B2T.DODECAHEDRON_VERTICES[vi])
        // .reduce((a,b) => a.plus(b), V3.O)
        // .unit()))
        const ss = new BRep(B2T.TETRAHEDRON_VERTICES.flatMap((v) => baseK.rotateAB(V3.Y, v).faces), false);
        //return ss
        return B2T.sphere().and(ss);
    }
    B2T.quaffle = quaffle;
    function extrudeFace(face, dir) {
        return new BRep(extrudeEdges(face.contour, face.surface.plane, dir)
            .faces.slice(0, -2)
            .concat(face, face.translate(dir.x, dir.y, dir.z).flipped(), face.holes.flatMap((hole) => extrudeEdges(hole, face.surface.plane.flipped(), dir).faces.slice(0, -2))), false);
    }
    B2T.extrudeFace = extrudeFace;
    function loadFonts() {
        return loadFont("fonts/FiraSansMedium.woff").then((font) => (B2T.defaultFont = font));
    }
    B2T.loadFonts = loadFonts;
    const loadedFonts = new Map();
    function loadFont(fontPath) {
        return new Promise(function (resolve, reject) {
            const font = loadedFonts.get(fontPath);
            if (font) {
                resolve(font);
            }
            else {
                load(fontPath, function (err, f) {
                    if (err) {
                        reject(err);
                    }
                    else {
                        loadedFonts.set(fontPath, f);
                        resolve(f);
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
            load("fonts/FiraSansMedium.woff", function (err, font) {
                if (err) {
                    throw new Error("Could not load font: " + err);
                }
                else {
                    B2T.defaultFont = font;
                    callback();
                }
            });
        }
    }
    B2T.loadFontsAsync = loadFontsAsync;
    /**
     * Create the [BRep] of a string rendered in a font.
     * @param text
     * @param size
     * @param depth
     * @param font An opentype.js font.
     */
    function text(text, size, depth = 1, font = B2T.defaultFont) {
        const path = font.getPath(text, 0, 0, size);
        const subpaths = [];
        path.commands.forEach((c) => {
            if (c.type == "M") {
                subpaths.push([]);
            }
            subpaths.last.push(c);
        });
        const loops = subpaths.map((sp) => {
            const path = new Path();
            path.commands = sp;
            const loop = Edge.reversePath(edgePathFromSVG(path.toPathData(13))).map((e) => e.mirrorY());
            assert(Edge.isLoop(loop));
            return loop;
        });
        const faces = Face.assembleFacesFromLoops(loops, new PlaneSurface(P3.XY), PlaneFace);
        const generator = callsce("B2T.text", text, size, depth);
        return BRep.join(faces.map((face) => B2T.extrudeFace(face, V(0, 0, -depth))), generator);
    }
    B2T.text = text;
    function minorityReport() {
        const a = B2T.sphere();
        const b = B2T.text("LEO CROW", 64, 128)
            .scale(0.1 / 32)
            .translate(-0.5, -0.05, 1.2)
            .flipped();
        const c = B2T.sphere(0.98);
        return a.and(b).plus(c);
    }
    B2T.minorityReport = minorityReport;
    function whatever() {
        const iso = icosahedron();
        const numbersBRep = BRep.join(iso.faces.map((face, i) => {
            const numberBRep = text("" + (i + 1), 0.4, -2);
            const centroid = face.contour
                .map((edge) => edge.a)
                .reduce((a, b) => a.plus(b), V3.O)
                .div(3);
            const sys = M4.forSys(face.contour[0].aDir, centroid.cross(face.contour[0].aDir), centroid.unit(), centroid);
            return numberBRep.transform(sys.times(M4.translate(-numberBRep.getAABB().size().x / 2, -0.1, -0.04)));
        }));
        const s = sphere(0.9);
        //return iso.and(numbersBRep)
        return iso.and(s).and(numbersBRep);
        //return numbersBRep
    }
    B2T.whatever = whatever;
    function whatever3() {
        const t = B2T.torus(1, 2);
        return B2T.box(5, 5, 2).translate(-2.5, -2.5).minus(t);
    }
    B2T.whatever3 = whatever3;
    function d20() {
        const iso = icosahedron();
        const numbersBRep = BRep.join(iso.faces.map((face, i) => {
            const numberBRep = text("" + (i + 1), 0.4, -2);
            const centroid = face.contour
                .map((edge) => edge.a)
                .reduce((a, b) => a.plus(b), V3.O)
                .div(3);
            const sys = M4.forSys(face.contour[0].aDir, centroid.cross(face.contour[0].aDir), centroid.unit(), centroid);
            return numberBRep.transform(sys.times(M4.translate(-numberBRep.getAABB().size().x / 2, -0.1, -0.04)));
        }));
        const s = sphere(0.9);
        //return iso.and(numbersBRep)
        return iso.and(s).and(numbersBRep);
        //return numbersBRep
    }
    B2T.d20 = d20;
    function rotStep(edges, totalRadsOrAngles, countO) {
        const angles = "number" === typeof totalRadsOrAngles
            ? arrayFromFunction(countO, (i) => ((i + 1) / countO) * totalRadsOrAngles)
            : totalRadsOrAngles;
        const count = angles.length;
        const open = !eq(TAU, angles.last);
        const ribs = [
            edges,
            ...angles.map((phi) => {
                if (eq(TAU, phi)) {
                    return edges;
                }
                const matrix = M4.rotateZ(phi);
                return edges.map((edge) => edge.transform(matrix));
            }),
        ];
        const horizontalEdges = arrayFromFunction(count, (i) => {
            const ipp = (i + 1) % (count + 1);
            return arrayFromFunction(edges.length, (j) => {
                if (!eq0(edges[j].a.lengthXY())) {
                    return StraightEdge.throughPoints(ribs[i][j].a, ribs[ipp][j].a);
                }
                return undefined;
            });
        });
        const faces = [];
        let face;
        edges.forEach((edge, i) => {
            const ipp = (i + 1) % edges.length;
            // for straight edges perpendicular to the Z-axis, we only create one face.
            if (edge instanceof StraightEdge &&
                edge.curve.dir1.isPerpendicularTo(V3.Z)) {
                const flipped = edge.a.x > edge.b.x;
                const surface = new PlaneSurface(flipped ? new P3(V3.Z, edge.a.z) : new P3(V3.Z.negated(), -edge.a.z));
                if (open) {
                    const faceEdges = [];
                    if (!eq0(edge.a.x)) {
                        faceEdges.push(...arrayFromFunction(count, (j) => horizontalEdges[j][i]));
                    }
                    faceEdges.push(ribs[count][i]);
                    if (!eq0(edge.b.x)) {
                        faceEdges.push(...arrayFromFunction(count, (j) => horizontalEdges[count - j - 1][ipp].flipped()));
                    }
                    faceEdges.push(edge.flipped());
                    face = new PlaneFace(surface, faceEdges);
                }
                else {
                    const contour = flipped
                        ? arrayFromFunction(count, (j) => horizontalEdges[j][i])
                        : arrayFromFunction(count, (j) => horizontalEdges[count - j - 1][ipp].flipped());
                    let hole;
                    if (flipped && !eq0(edge.b.x)) {
                        hole = arrayFromFunction(count, (j) => horizontalEdges[count - j - 1][ipp].flipped());
                    }
                    else if (!flipped && !eq0(edge.a.x)) {
                        hole = arrayFromFunction(count, (j) => horizontalEdges[j][i]);
                    }
                    face = new PlaneFace(surface, contour, hole ? [hole] : []);
                }
                faces.push(face);
                return;
            }
            else if (edge instanceof StraightEdge) {
                if (eq0(edge.a.lengthXY()) && eq0(edge.b.lengthXY())) {
                    return;
                }
            }
            for (let r = 0; r < count; r++) {
                const rpp = (r + 1) % (count + 1);
                const faceEdges = [
                    ribs[r][i].flipped(),
                    horizontalEdges[r][i],
                    ribs[rpp][i],
                    horizontalEdges[r][ipp] && horizontalEdges[r][ipp].flipped(),
                ].filter((x) => x);
                let surface;
                if (edge instanceof StraightEdge) {
                    surface = new PlaneSurface(P3.throughPoints(faceEdges[0].a, faceEdges[1].a, faceEdges[2].a));
                }
                else {
                    const maxX = edges[i].getAABB().max.x;
                    const phi = angles[r], prevPhi = 0 == r ? 0 : angles[r - 1];
                    const offset = V3.polar(maxX, prevPhi).to(V3.polar(maxX, phi));
                    surface = projectCurve(ribs[r][i].curve, offset, false);
                }
                faces.push(Face.create(surface, faceEdges));
            }
        });
        if (open) {
            const endFaceEdges = ribs[count].map((edge) => edge.flipped()).reverse();
            const endFace = new PlaneFace(new PlaneSurface(P3.ZX.rotateZ(getLast(angles))), endFaceEdges);
            faces.push(new PlaneFace(new PlaneSurface(P3.ZX.flipped()), edges), endFace);
        }
        return new BRep(faces, new PlaneSurface(P3.ZX).edgeLoopCCW(edges));
    }
    B2T.rotStep = rotStep;
    function fixEdges(edges) {
        return edges.flatMap((edge) => {
            const c = edge.curve;
            if (c instanceof EllipseCurve && c.tMin === -PI && c.tMax === PI) {
                const splitEdges = edge.minT < 0 && edge.maxT > 0 ? edge.split(0) : [edge];
                return splitEdges.map((edge) => {
                    if (edge.minT >= 0) {
                        return createEdge(new EllipseCurve(c.center, c.f1, c.f2, max(0, c.tMin), c.tMax), edge.a, edge.b, edge.aT, edge.bT, undefined, edge.aDir, edge.bDir, edge.name);
                    }
                    else {
                        // "rotate" the curve
                        return createEdge(new EllipseCurve(c.center, c.f1.negated(), c.f2.negated(), c.tMin + PI, min(PI, c.tMax + PI)), edge.a, edge.b, edge.aT + PI, edge.bT + PI, undefined, edge.aDir, edge.bDir, edge.name);
                    }
                });
            }
            if (c instanceof BezierCurve) {
                if (edge.a.like(edge.b)) {
                    return edge.split(lerp(edge.aT, edge.bT, 0.5));
                }
            }
            return edge;
        });
    }
    B2T.fixEdges = fixEdges;
    /**
     * Create a [BRep] by projecting edges created by joining vertices with straight edges.
     * @param baseVertices
     * @param baseFacePlane
     * @param offset
     * @param name
     * @param generator
     */
    function extrudeVertices(baseVertices, baseFacePlane, offset, name, generator) {
        assert(baseVertices.every((v) => v instanceof V3), "baseVertices.every(v => v instanceof V3)");
        assertInst(P3, baseFacePlane);
        assertVectors(offset);
        if (baseFacePlane.normal1.dot(offset) > 0)
            baseFacePlane = baseFacePlane.flipped();
        const edges = StraightEdge.chain(baseVertices, true);
        generator =
            generator ||
                callsce("B2T.extrudeVertices", baseVertices, baseFacePlane, offset, name);
        return B2T.extrudeEdges(edges, baseFacePlane, offset, name, generator);
    }
    B2T.extrudeVertices = extrudeVertices;
    /**
     * Create a tetrahedron (3 sided pyramid) [BRep].
     * `a`, `b`, `c` and `d` can be in any order. The only constraint is that they cannot be on a common plane.
     * The resulting tetrahedron will always have outwards facing faces.
     * @param a
     * @param b
     * @param c
     * @param d
     * @param name
     */
    function tetrahedron(a, b, c, d, name = "tetra" + getGlobalId()) {
        assertVectors(a, b, c, d);
        const dDistance = P3.throughPoints(a, b, c).distanceToPointSigned(d);
        if (eq0(dDistance)) {
            throw new Error("four points are coplanar");
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
            new PlaneFace(PlaneSurface.throughPoints(a, b, c), [ab, bc, ac.flipped()], [], name + "abc"),
            new PlaneFace(PlaneSurface.throughPoints(a, d, b), [ad, bd.flipped(), ab.flipped()], [], name + "adb"),
            new PlaneFace(PlaneSurface.throughPoints(b, d, c), [bd, cd.flipped(), bc.flipped()], [], name + "bdc"),
            new PlaneFace(PlaneSurface.throughPoints(c, d, a), [cd, ad.flipped(), ac], [], name + "cda"),
        ];
        const gen = callsce("B2T.tetrahedron", a, b, c, d);
        return new BRep(faces, false, gen);
    }
    B2T.tetrahedron = tetrahedron;
    const b = 1 / GOLDEN_RATIO, c = 2 - GOLDEN_RATIO;
    B2T.TETRAHEDRON_VERTICES = [
        new V3(1, 0, -SQRT1_2),
        new V3(-1, 0, -SQRT1_2),
        new V3(0, -1, SQRT1_2),
        new V3(0, 1, SQRT1_2),
    ].map((v) => v.unit());
    B2T.DODECAHEDRON_VERTICES = [
        new V3(c, 0, 1),
        new V3(-c, 0, 1),
        new V3(-b, b, b),
        new V3(0, 1, c),
        new V3(b, b, b),
        new V3(b, -b, b),
        new V3(0, -1, c),
        new V3(-b, -b, b),
        new V3(c, 0, -1),
        new V3(-c, 0, -1),
        new V3(-b, -b, -b),
        new V3(0, -1, -c),
        new V3(b, -b, -b),
        new V3(b, b, -b),
        new V3(0, 1, -c),
        new V3(-b, b, -b),
        new V3(1, c, 0),
        new V3(-1, c, 0),
        new V3(-1, -c, 0),
        new V3(1, -c, 0),
    ].map((v) => v.unit());
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
        [7, 1, 2, 17, 18],
    ];
    B2T.OCTAHEDRON_VERTICES = [
        new V3(1, 0, 0),
        new V3(-1, 0, 0),
        new V3(0, 1, 0),
        new V3(0, -1, 0),
        new V3(0, 0, 1),
        new V3(0, 0, -1),
    ];
    B2T.OCTAHEDRON_FACE_VERTICES = [
        [0, 2, 4],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],
        [2, 0, 5],
        [1, 2, 5],
        [3, 1, 5],
        [0, 3, 5],
    ];
    const { x: s, y: t } = new V3(1, GOLDEN_RATIO, 0).unit();
    B2T.ICOSAHEDRON_VERTICES = [
        new V3(-s, t, 0),
        new V3(s, t, 0),
        new V3(-s, -t, 0),
        new V3(s, -t, 0),
        new V3(0, -s, t),
        new V3(0, s, t),
        new V3(0, -s, -t),
        new V3(0, s, -t),
        new V3(t, 0, -s),
        new V3(t, 0, s),
        new V3(-t, 0, -s),
        new V3(-t, 0, s),
    ];
    B2T.ICOSAHEDRON_FACE_VERTICES = [
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
        [9, 8, 1],
    ];
    /**
     * Create a dodecahedron [BRep]. The vertices are on the unit sphere.
     */
    function dodecahedron() {
        return makePlatonic(B2T.DODECAHEDRON_VERTICES, B2T.DODECAHEDRON_FACE_VERTICES, "B2T.dodecahedron()");
    }
    B2T.dodecahedron = dodecahedron;
    /**
     * Create an octahedron [BRep]. The vertices are on the unit sphere.
     */
    function octahedron() {
        return makePlatonic(B2T.OCTAHEDRON_VERTICES, B2T.OCTAHEDRON_FACE_VERTICES, "B2T.octahedron()");
    }
    B2T.octahedron = octahedron;
    /**
     * Create an icosahedron [BRep]. The vertices are on the unit sphere.
     */
    function icosahedron() {
        return makePlatonic(B2T.ICOSAHEDRON_VERTICES, B2T.ICOSAHEDRON_FACE_VERTICES, "B2T.icosahedron()");
    }
    B2T.icosahedron = icosahedron;
    function makePlatonic(VS, FVIS, generator) {
        const edgeMap = new Map();
        const faces = FVIS.map((faceIndexes) => {
            const surface = PlaneSurface.throughPoints(VS[faceIndexes[0]], VS[faceIndexes[1]], VS[faceIndexes[2]]);
            const contour = arrayFromFunction(faceIndexes.length, (i) => {
                const ipp = (i + 1) % faceIndexes.length;
                const iA = faceIndexes[i], iB = faceIndexes[ipp];
                const iMin = min(iA, iB), iMax = max(iA, iB), edgeID = iMin * VS.length + iMax;
                let edge = edgeMap.get(edgeID);
                !edge &&
                    edgeMap.set(edgeID, (edge = StraightEdge.throughPoints(VS[iMin], VS[iMax])));
                return iA < iB ? edge : edge.flipped();
            });
            return new PlaneFace(surface, contour);
        });
        return new BRep(faces, false, generator);
    }
    /**
     * Create a [BRep] by projecting a number of edges onto a point.
     * @param baseEdges The edges forming the base of the pyramid.
     * @param apex The tip of the pyramid.
     * @param name
     */
    function pyramidEdges(baseEdges, apex, name = "pyramid" + getGlobalId()) {
        assertInst(Edge, ...baseEdges);
        assertVectors(apex);
        const ribs = baseEdges.map((baseEdge) => StraightEdge.throughPoints(apex, baseEdge.a));
        const faces = baseEdges.map((baseEdge, i) => {
            const faceName = name + "Wall" + i;
            const ipp = (i + 1) % baseEdges.length;
            const faceEdges = [ribs[i], baseEdge, ribs[ipp].flipped()];
            const surface = projectPointCurve(baseEdge.curve, baseEdge.minT, baseEdge.maxT, apex, baseEdge.deltaT() < 0);
            return Face.create(surface, faceEdges, undefined, faceName);
        });
        const baseSurface = new PlaneSurface(P3.XY).flipped();
        const bottomFace = Face.create(baseSurface, Edge.reversePath(baseEdges));
        faces.push(bottomFace);
        const generator = callsce("B2T.pyramidEdges", baseEdges, apex, name);
        return new BRep(faces, false, generator);
    }
    B2T.pyramidEdges = pyramidEdges;
    function fromBPT(bpt) {
        const lineRegex = /.+/g;
        const readLine = () => lineRegex.exec(bpt)[0];
        const readLineNumbers = () => readLine()
            .trim()
            .split(/\s+/)
            .map((s) => parseFloat(s));
        const numOfPatches = parseInt(readLine());
        const faces = arrayFromFunction(numOfPatches, () => {
            const [pointsUCount, pointsVCount] = readLineNumbers();
            const points = Array.from({ length: (pointsUCount + 1) * (pointsVCount + 1) }, () => VV(...readLineNumbers(), 1));
            const surface = new NURBSSurface(points, NURBS.bezierKnots(pointsUCount), NURBS.bezierKnots(pointsVCount), pointsUCount, pointsVCount, 0, 1, 0, 1);
            return surface;
        });
        return faces;
    }
    B2T.fromBPT = fromBPT;
})(B2T || (B2T = {}));

class CustomPlane extends P3 {
    constructor(anchor, right, up, name = "CustomPlane" + getGlobalId(), color = chroma.random().gl(), uMin = -500, uMax = 500, vMin = -500, vMax = 500) {
        const { normal1, w } = P3.forAnchorAndPlaneVectors(anchor, right, up);
        super(normal1, w);
        this.up = up;
        this.right = right;
        this.uMin = uMin;
        this.uMax = uMax;
        this.vMin = vMin;
        this.vMax = vMax;
        this.name = name;
        this.color = color;
    }
    get plane() {
        return this;
    }
    toPlaneSurface() {
        return new PlaneSurface(this, this.right, this.up);
    }
    toSource() {
        return callsce("new CustomPlane", this.anchor, this.right, this.up, this.name, this.color, this.uMin, this.uMax, this.vMin, this.vMax);
    }
    static forPlane(plane, color = GL_COLOR_BLACK, name) {
        //assert(!name)
        const up = plane.normal1.getPerpendicular().unit(), right = up.cross(plane.normal1);
        return new CustomPlane(plane.anchor, right, up, name, color);
    }
    static fromPlaneSurface(surface) {
        return new CustomPlane(surface.plane.anchor, surface.right, surface.up, "genCustomPlane" + getGlobalId());
    }
    distanceTo(line, mindist) {
        return min$1([
            new L3(this.anchor.plus(this.right.times(this.uMin)), this.up),
            new L3(this.anchor.plus(this.right.times(this.uMax)), this.up),
            new L3(this.anchor.plus(this.up.times(this.vMin)), this.right),
            new L3(this.anchor.plus(this.up.times(this.vMax)), this.right),
        ]
            .map((line2, line2Index) => {
            const info = line2.infoClosestToLine(line);
            if ((isNaN(info.t) || // parallel LINES
                (line2Index < 2 && this.vMin <= info.t && info.t <= this.vMax) ||
                (line2Index >= 2 && this.uMin <= info.t && info.t <= this.uMax)) &&
                info.distance <= mindist) {
                return info.s;
            }
            else {
                return Infinity;
            }
        }));
    }
    distanceTo2(line, mindist) {
        return min$1([
            new L3(this.anchor.plus(this.right.times(this.uMin)), this.up),
            new L3(this.anchor.plus(this.right.times(this.uMax)), this.up),
            new L3(this.anchor.plus(this.up.times(this.vMin)), this.right),
            new L3(this.anchor.plus(this.up.times(this.vMax)), this.right),
        ]
            .map((line2, line2Index) => {
            const info = line2.infoClosestToLine(line);
            if ((isNaN(info.t) || // parallel LINES
                (line2Index < 2 && this.vMin <= info.t && info.t <= this.vMax) ||
                (line2Index >= 2 && this.uMin <= info.t && info.t <= this.uMax)) &&
                info.distance <= mindist) {
                return info.distance;
            }
            else {
                return Infinity;
            }
        }));
    }
}

class Edge extends Transformable {
    constructor(curve, a, b, aT, bT, flippedOf, name) {
        super();
        this.curve = curve;
        this.a = a;
        this.b = b;
        this.aT = aT;
        this.bT = bT;
        this.flippedOf = flippedOf;
        this.name = name;
        assertNumbers(aT, bT);
        assert(!eq(aT, bT));
        assertVectors(a, b);
        assertf(() => curve instanceof Curve, curve);
        assertf(() => !curve.isValidT || (curve.isValidT(aT) && curve.isValidT(bT)), aT, bT, curve);
        //if (curve instanceof PICurve) {
        //    assertf(() => curve.at(aT).to(a).length() < 0.1, ''+curve.at(aT)+a)
        //    assertf(() => curve.at(bT).to(b).length() < 0.1, '' + curve.at(bT) + b)
        //} else {
        assertf(() => curve.at(aT).like(a), () => "" + curve.at(aT) + a + " aT should have been " + curve.pointT(a));
        assertf(() => curve.at(bT).like(b), () => "" + curve.at(bT) + b + " bT should have been " + curve.pointT(b));
        //}
        assertf(() => fuzzyBetween(aT, curve.tMin, curve.tMax), aT, curve.tMin, curve.tMax);
        assertf(() => fuzzyBetween(bT, curve.tMin, curve.tMax), bT, curve.tMin, curve.tMax);
        assert(!a.like(b), "!a.like(b)" + a + b);
        this.aT = clamp(aT, curve.tMin, curve.tMax);
        this.bT = clamp(bT, curve.tMin, curve.tMax);
        this.reversed = this.aT > this.bT;
    }
    get minT() {
        return Math.min(this.aT, this.bT);
    }
    get maxT() {
        return Math.max(this.aT, this.bT);
    }
    static isLoop(loop) {
        return loop.every((edge, i) => edge.b.like(loop[(i + 1) % loop.length].a));
    }
    static edgesIntersect(e1, e2) {
        // TODO: still getting some NaNs here..
        assertNumbers(e1.curve.hlol, e2.curve.hlol);
        assertInst(Edge, e1, e2);
        if (e1.curve.hlol < e2.curve.hlol) {
            [e2, e1] = [e1, e2];
        }
        const sts = e1.curve.isInfosWithCurve(e2.curve);
        if (sts.some((info) => isNaN(info.tThis) || isNaN(info.tOther))) {
            console.log(e1.sce);
            console.log(e2.sce);
            assert(false);
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
            assert(edge.b.like(edges[j].a), `edges[${i}].b != edges[${j}].a (${edges[i].b.sce} != ${edges[j].a.sce})`);
        });
    }
    static reversePath(path, doReverse = true) {
        return doReverse
            ? arrayFromFunction(path.length, (i) => path[path.length - 1 - i].flipped())
            : path;
    }
    toString() {
        return callsce("new " + this.constructor.name, this.curve, this.a, this.b, this.aT, this.bT, undefined, this.aDir, this.bDir);
    }
    colinearToLine(line) {
        return this.curve instanceof L3 && this.curve.isColinearTo(line);
    }
    tValueInside(t) {
        return this.aT < this.bT
            ? lt(this.aT, t) && lt(t, this.bT)
            : lt(this.bT, t) && lt(t, this.aT);
    }
    isValidT(t) {
        return this.aT < this.bT
            ? le(this.aT, t) && le(t, this.bT)
            : le(this.bT, t) && le(t, this.aT);
    }
    clampedT(t) {
        return this.aT < this.bT
            ? clamp(t, this.aT, this.bT)
            : clamp(t, this.bT, this.aT);
    }
    /**
     * this is equals-equals. "isColinearTo" might make more sense but can't be used, because you can't get a
     * consistent hashCode for colinear curves
     * @param obj
     * @returns
     */
    equals(obj) {
        return (this === obj ||
            (this.constructor == obj.constructor &&
                this.a.equals(obj.a) &&
                this.b.equals(obj.b) &&
                this.curve.equals(obj.curve)));
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
        return (this === edge ||
            (edge instanceof Edge &&
                this.curve.isColinearTo(edge.curve) &&
                this.a.like(edge.a) &&
                this.b.like(edge.b)));
    }
    isCanon() {
        return !this.reversed;
    }
    getCanon() {
        return this.reversed ? this.flipped() : this;
    }
    overlaps(edge, noback) {
        assert(this.curve.isColinearTo(edge.curve));
        const edgeAT = this.curve.containsPoint(edge.a) && this.curve.pointT(edge.a);
        const edgeBT = this.curve.containsPoint(edge.b) && this.curve.pointT(edge.b);
        if (false === edgeAT && false === edgeBT) {
            return noback ? false : edge.overlaps(this, true);
        }
        return !(le(edge.maxT, this.minT) || le(this.maxT, edge.minT));
    }
    getAABB() {
        const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
        this.curve.roots().forEach((ts, dim) => {
            ts.forEach((t) => {
                if (lt(this.minT, t) && lt(t, this.maxT)) {
                    min[dim] = Math.min(min[dim], this.curve.at(t).e(dim));
                    max[dim] = Math.max(max[dim], this.curve.at(t).e(dim));
                }
            });
        });
        const aabb = new AABB(V(min), V(max));
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
    deltaTSign() {
        return sign(this.bT - this.aT);
    }
    atAvgT() {
        return this.curve.at((this.minT + this.maxT) / 2);
    }
    /**
     * Whether two edge loops are equal. Takes into account that two loops need not start with the same edge.
     * @param loop1
     * @param loop2
     */
    static loopsEqual(loop1, loop2) {
        return (loop1.length == loop2.length &&
            arrayRange(0, loop1.length, 1).some((offset) => loop1.every((edge, i) => edge.equals(loop2[(offset + i) % loop1.length]))));
    }
}

class StraightEdge extends Edge {
    constructor(line, a, b, aT, bT, flippedOf, name) {
        super(line, a, b, aT, bT, flippedOf, name);
        this.flippedOf = flippedOf;
        assertInst(L3, line);
        !flippedOf || assertInst(StraightEdge, flippedOf);
        !name || assertf(() => "string" === typeof name, name);
        this.tangent =
            this.aT < this.bT ? this.curve.dir1 : this.curve.dir1.negated();
    }
    get aDir() {
        return this.tangent;
    }
    get bDir() {
        return this.tangent;
    }
    static throughPoints(a, b, name) {
        return new StraightEdge(L3.throughPoints(a, b, 0, a.to(b).length()), a, b, 0, a.to(b).length(), undefined, name);
    }
    /**
     * Create a list of StraightEdges from a list of vertices.
     * @param vertices
     * @param closed Whether to connect the first and last vertices. Defaults to true.
     * @returns
     */
    static chain(vertices, closed = true) {
        const vc = vertices.length;
        return arrayFromFunction(closed ? vc : vc - 1, (i) => StraightEdge.throughPoints(vertices[i], vertices[(i + 1) % vc]));
    }
    toSource() {
        return callsce("new StraightEdge", this.curve, this.a, this.b, this.aT, this.bT);
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
        const edgeT = snap2(this.curve.isTWithPlane(plane), this.aT, this.bT);
        return this.minT <= edgeT && edgeT <= this.maxT ? [edgeT] : [];
    }
    edgeISTsWithSurface(surface) {
        if (surface instanceof PlaneSurface) {
            return this.edgeISTsWithPlane(surface.plane);
        }
        else {
            return surface
                .isTsForLine(this.curve)
                .map((edgeT) => snap2(edgeT, this.aT, this.bT))
                .filter((edgeT) => this.minT <= edgeT && edgeT <= this.maxT);
        }
    }
    tangentAt() {
        return this.tangent;
    }
    flipped() {
        return (this.flippedOf ||
            (this.flippedOf = new StraightEdge(this.curve, this.b, this.a, this.bT, this.aT, this, this.name)));
    }
    transform(m4, desc) {
        const lineDir1TransLength = m4
            .transformVector2(this.curve.dir1, this.curve.anchor)
            .length();
        const curve = this.curve.transform(m4);
        const a = m4.transformPoint(this.a);
        const b = m4.transformPoint(this.b);
        return new StraightEdge(curve, a, b, m4.isNoProj() ? this.aT * lineDir1TransLength : curve.pointT(a), m4.isNoProj() ? this.bT * lineDir1TransLength : curve.pointT(b), undefined, "" + this.name + desc);
    }
    transform4(m4, desc) {
        const lineDir1TransLength = m4
            .transformVector2(this.curve.dir1, this.curve.anchor)
            .length();
        const curve = this.curve.transform4(m4);
        const a = m4.transformPoint(this.a);
        const b = m4.transformPoint(this.b);
        return new StraightEdge(curve, a, b, m4.isNoProj() ? this.aT * lineDir1TransLength : curve.pointT(a), m4.isNoProj() ? this.bT * lineDir1TransLength : curve.pointT(b), undefined, "" + this.name + desc);
    }
    isCoEdge(edge) {
        return (this === edge ||
            this === edge.flippedOf ||
            (edge.constructor instanceof StraightEdge &&
                ((this.a.like(edge.a) && this.b.like(edge.b)) ||
                    (this.a.like(edge.b) && this.b.like(edge.a)))));
    }
    getEdgeT(p) {
        assertVectors(p);
        let edgeT = p.minus(this.curve.anchor).dot(this.curve.dir1);
        if (!eq0(this.curve.at(edgeT).distanceTo(p))) {
            return;
        }
        edgeT = snap2(edgeT, this.aT, this.bT);
        return this.minT <= edgeT && edgeT <= this.maxT ? edgeT : undefined;
    }
    split(t) {
        const p = this.curve.at(t);
        return [
            new StraightEdge(this.curve, this.a, p, this.aT, t, undefined, this.name + "left"),
            new StraightEdge(this.curve, p, this.b, t, this.bT, undefined, this.name + "right"),
        ];
    }
}

class PCurveEdge extends Edge {
    constructor(curve, a, b, aT, bT, flippedOf, aDir, bDir, name) {
        super(curve, a, b, aT, bT, flippedOf, name);
        this.flippedOf = flippedOf;
        this.aDir = aDir;
        this.bDir = bDir;
        assertVectors(aDir, bDir);
        assertf(() => !aDir.likeO(), curve);
        assertf(() => !bDir.likeO(), curve);
        if (!(curve instanceof PICurve)) {
            // TODO
            assertf(() => curve.tangentAt(aT).likeOrReversed(aDir), "" + aT + curve.tangentAt(aT).sce + " " + aDir.sce);
            assertf(() => curve.tangentAt(bT).likeOrReversed(bDir), "" + bT + curve.tangentAt(bT).sce + " " + bDir.sce);
        }
        assert(this.reversed === this.aDir.dot(curve.tangentAt(aT)) < 0, aT +
            " " +
            bT +
            " " +
            curve.constructor.name +
            " " +
            this.aDir.sce +
            " " +
            this.bDir.sce +
            " " +
            curve.tangentAt(aT));
        assert(this.reversed === this.bDir.dot(curve.tangentAt(bT)) < 0, aT +
            " " +
            bT +
            " " +
            curve.constructor.name +
            " " +
            this.aDir.sce +
            " " +
            this.bDir.sce +
            " " +
            curve.tangentAt(aT));
    }
    static forCurveAndTs(curve, aT, bT, name) {
        return new PCurveEdge(curve, curve.at(aT), curve.at(bT), aT, bT, undefined, aT < bT ? curve.tangentAt(aT) : curve.tangentAt(aT).negated(), aT < bT ? curve.tangentAt(bT) : curve.tangentAt(bT).negated(), name);
    }
    toSource() {
        return callsce("new PCurveEdge", this.curve, this.a, this.b, this.aT, this.bT, undefined, this.aDir, this.bDir, this.name);
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
    edgeISTsWithSurface(surface) {
        return this.curve
            .isTsWithSurface(surface)
            .map((edgeT) => snap2(edgeT, this.aT, this.bT))
            .filter((edgeT) => this.minT <= edgeT && edgeT <= this.maxT);
    }
    edgeISTsWithPlane(surface) {
        return this.curve
            .isTsWithPlane(surface)
            .map((edgeT) => snap2(edgeT, this.aT, this.bT))
            .filter((edgeT) => this.minT <= edgeT && edgeT <= this.maxT);
    }
    tangentAt(t) {
        return !this.reversed
            ? this.curve.tangentAt(t)
            : this.curve.tangentAt(t).negated();
    }
    flipped() {
        return (this.flippedOf ||
            (this.flippedOf = new PCurveEdge(this.curve, this.b, this.a, this.bT, this.aT, this, this.bDir.negated(), this.aDir.negated(), this.name)));
    }
    transform(m4, desc) {
        return new PCurveEdge(this.curve.transform(m4), m4.transformPoint(this.a), m4.transformPoint(this.b), this.aT, this.bT, undefined, m4.transformVector(this.aDir), m4.transformVector(this.bDir), "" + this.name + desc);
    }
    transform4(m4, desc) {
        const a_ = m4.transformPoint(this.a);
        const b_ = m4.transformPoint(this.b);
        const curve_ = this.curve.transform4(m4);
        return new PCurveEdge(curve_, a_, b_, snap(curve_.pointT(a_), this.aT), snap(curve_.pointT(b_), this.bT), undefined, m4.transformVector(this.aDir), m4.transformVector(this.bDir), "" + this.name + desc);
    }
    isCoEdge(edge) {
        return (this === edge ||
            this === edge.flippedOf ||
            (this.curve.isColinearTo(edge.curve) &&
                ((this.a.like(edge.a) && this.b.like(edge.b)) ||
                    (this.a.like(edge.b) && this.b.like(edge.a)))));
    }
    split(t) {
        const p = this.curve.at(t);
        const dir = this.tangentAt(t);
        return [
            new PCurveEdge(this.curve, this.a, p, this.aT, t, undefined, this.aDir, dir, this.name + "left"),
            new PCurveEdge(this.curve, p, this.b, t, this.bT, undefined, this.aDir, dir, this.name + "right"),
        ];
    }
}

function edgePathFromSVG(pathString) {
    let currentPos = undefined;
    const parsed = new SVGPathData(pathString)
        .toAbs()
        .normalizeHVZ()
        .sanitize(NLA_PRECISION)
        .annotateArcs().commands;
    const path = [];
    for (const c of parsed) {
        assert("x" in c && "y" in c);
        const endPos = new V3(c.x, c.y, 0);
        switch (c.type) {
            case SVGPathData.LINE_TO:
                path.push(StraightEdge.throughPoints(currentPos, endPos));
                break;
            case SVGPathData.CURVE_TO: {
                const c1 = new V3(c.x1, c.y1, 0);
                const c2 = new V3(c.x2, c.y2, 0);
                const curve = new BezierCurve(currentPos, c1, c2, endPos, 0, 1);
                const edge = new PCurveEdge(curve, currentPos, endPos, 0, 1, undefined, curve.tangentAt(0), curve.tangentAt(1));
                path.push(edge);
                break;
            }
            case SVGPathData.QUAD_TO: {
                const c1 = new V3(c.x1, c.y1, 0);
                const curve = ParabolaCurve.quadratic(currentPos, c1, endPos).rightAngled();
                const edge = new PCurveEdge(curve, currentPos, endPos, curve.tMin, curve.tMax, undefined, curve.tangentAt(curve.tMin), curve.tangentAt(curve.tMax));
                path.push(edge);
                break;
            }
            case SVGPathData.ARC: {
                const phi1 = c.phi1 * DEG, phi2 = c.phi2 * DEG, [phiMin, phiMax] = [phi1, phi2].sort(MINUS);
                const stops = arrayRange(-3, 4, 1)
                    .map((n) => n * PI$1)
                    .filter((stop) => phiMin <= stop && stop <= phiMax);
                const center = V(c.cX, c.cY);
                const f1 = V3.polar(c.rX, c.xRot * DEG);
                const f2 = V3.polar(c.rY, c.xRot * DEG + Math.PI / 2);
                const edges = getIntervals(stops, phiMin, phiMax).map(([t1, t2]) => {
                    const deltaT = t2 - t1;
                    const t1_ = mod(t1, TAU);
                    const t2_ = t1_ + deltaT;
                    assert(t1_ >= 0 == t2_ >= 0);
                    const gtPI = t1_ > PI$1 || t2_ > PI$1;
                    const aT = gtPI ? t1_ - PI$1 : t1_;
                    const bT = gtPI ? t2_ - PI$1 : t2_;
                    const curve = new EllipseCurve(center, gtPI ? f1.negated() : f1, gtPI ? f2.negated() : f2);
                    const a = phi1 == t1 ? currentPos : phi2 == t1 ? endPos : curve.at(aT);
                    const b = phi1 == t2 ? currentPos : phi2 == t2 ? endPos : curve.at(bT);
                    return new PCurveEdge(curve, a, b, aT, bT, undefined, curve.tangentAt(aT), curve.tangentAt(bT));
                });
                path.push(...(c.phiDelta > 0 ? edges : Edge.reversePath(edges)));
                break;
            }
        }
        currentPos = endPos;
    }
    return path;
}
/**
 * Create an axis-aligned rectangle of edges on the XY-plane with the bottom-left corner on the origin.
 * @param width
 * @param height
 */
function edgeRect(width = 1, height = width) {
    const vertices = [
        new V3(0, 0, 0),
        new V3(width, 0, 0),
        new V3(width, height, 0),
        new V3(0, height, 0),
    ];
    return StraightEdge.chain(vertices);
}
function ngon(n = 3, radius = 1) {
    return StraightEdge.chain(arrayFromFunction(n, (i) => V3.polar(radius, (TAU * i) / n)));
}
function star(pointCount = 5, r0 = 1, r1 = 0.5) {
    const vertices = arrayFromFunction(pointCount * 2, (i) => V3.polar(0 == i % 2 ? r0 : r1, (TAU * i) / pointCount / 2));
    return StraightEdge.chain(vertices);
}
function createEdge(curve, a, b, aT, bT, flippedOf, aDir, bDir, name) {
    if (curve instanceof L3) {
        return new StraightEdge(curve, a, b, aT, bT, flippedOf, name);
    }
    else {
        return new PCurveEdge(curve, a, b, aT, bT, flippedOf, aDir, bDir, name);
    }
}
function edgeForCurveAndTs(curve, aT = curve.tMin, bT = curve.tMax) {
    return createEdge(curve, curve.at(aT), curve.at(bT), aT, bT, undefined, aT < bT ? curve.tangentAt(aT) : curve.tangentAt(aT).negated(), aT < bT ? curve.tangentAt(bT) : curve.tangentAt(bT).negated());
}
function reuleaux(n = 3, radius = 1) {
    assert(3 <= n);
    assert(1 == n % 2);
    const corners = arrayFromFunction(n, (i) => V3.polar(radius, (TAU * i) / n));
    return arrayFromFunction(n, (i) => {
        const aI = (i + floor(n / 2)) % n, bI = (i + ceil(n / 2)) % n;
        const a = corners[aI], b = corners[bI];
        const center = corners[i];
        const f1 = center.to(a), curve = new EllipseCurve(center, f1, V3.Z.cross(f1));
        return createEdge(curve, a, b, 0, curve.pointT(b), undefined, V3.Z.cross(f1), V3.Z.cross(center.to(b)));
    });
}
function round$1(edges, radius) {
    if (eq0(radius)) {
        return edges;
    }
    const corners = edges.map((edge, i) => {
        const j = (i + 1) % edges.length, nextEdge = edges[j];
        if (!edge.b.like(nextEdge.a))
            return undefined;
        const angleToNext = edge.bDir.angleTo(nextEdge.aDir);
        const c1 = edge.curve, c2 = nextEdge.curve;
        if (c1 instanceof L3 && c2 instanceof L3) {
            const normal = c1.dir1.cross(c2.dir1);
            if (eq0(angleToNext))
                return undefined;
            const l1inside = normal.cross(c1.dir1), l2inside = normal.cross(c2.dir1);
            const l1offset = c1.transform(M4.translate(l1inside.toLength(radius)));
            const l2offset = c2.transform(M4.translate(l2inside.toLength(radius)));
            const center = l1offset.isInfoWithLine(l2offset);
            if (!center)
                throw new Error("tangential curves");
            const cornerA = center.plus(l1inside.toLength(-radius));
            const cornerB = center.plus(l2inside.toLength(-radius));
            const f1 = l1inside.toLength(-radius);
            const curve = new EllipseCurve(center, f1, normal.cross(f1).toLength(radius));
            const cornerEdge = createEdge(curve, cornerA, cornerB, 0, curve.pointT(cornerB), undefined, c1.dir1, c2.dir1);
            return cornerEdge;
        }
        else {
            return arbitraryCorner(edge, nextEdge, radius);
        }
    });
    const result = edges.flatMap((edge, i) => {
        const h = (i + edges.length - 1) % edges.length;
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
        const newEdge = createEdge(edge.curve, a, b, aT, bT, undefined, aDir, bDir);
        return !nextCorner ? newEdge : [newEdge, nextCorner];
    });
    return result;
}
function arbitraryCorner(e1, e2, radius) {
    const c1 = e1.curve, c2 = e2.curve;
    function f([t1, t2]) {
        const p1 = c1.at(t1), p2 = c2.at(t2);
        const dp1 = c1.tangentAt(t1), dp2 = c2.tangentAt(t2);
        const virtualPlaneNormal = dp1.cross(dp2);
        const normal1 = virtualPlaneNormal.cross(dp1).unit(), normal2 = virtualPlaneNormal.cross(dp2).unit();
        const dirCross = normal1.cross(normal2);
        if (virtualPlaneNormal.likeO()) {
            assert(false);
        } // lines parallel
        const p1p2 = p1.to(p2);
        // check if distance is zero (see also L3.distanceToLine)
        if (!eq0(p1p2.dot(virtualPlaneNormal))) {
            assert(false);
        }
        const dist1 = p1p2.cross(normal2).dot(dirCross) / dirCross.squared();
        const dist2 = p1p2.cross(normal1).dot(dirCross) / dirCross.squared();
        const g1 = p1.plus(normal1.times(dist1));
        const g2 = p2.plus(normal2.times(dist2));
        assert(g1.like(g2));
        return [abs(dist1) - radius, abs(dist2) - radius];
    }
    const startT1 = e1.bT - (radius * sign(e1.deltaT())) / e1.bDir.length();
    const startT2 = e2.aT + (radius * sign(e2.deltaT())) / e2.aDir.length();
    const [t1, t2] = newtonIterate(f, [startT1, startT2]);
    const cornerA = e1.curve.at(t1);
    const cornerB = e2.curve.at(t2);
    const dp1 = c1.tangentAt(t1), dp2 = c2.tangentAt(t2);
    const virtualPlaneNormal = dp1.cross(dp2);
    const normal1 = virtualPlaneNormal.cross(dp1).unit();
    const f1 = normal1.toLength(-radius);
    const center = cornerA.minus(f1);
    const curve = new EllipseCurve(center, f1, virtualPlaneNormal.cross(f1).toLength(radius));
    const cornerEdge = createEdge(curve, cornerA, cornerB, 0, curve.pointT(cornerB), undefined, c1.tangentAt(t1), c2.tangentAt(t2));
    return cornerEdge;
}

/**
 * Created by aval on 19.04.2017.
 */
class FaceInfoFactory {
    static makeStatic(staticInfo) {
        return new (class extends FaceInfoFactory {
            constructor() {
                super();
            }
            info(surface, contour, holes) {
                return staticInfo;
            }
        })();
    }
    info(surface, contour, holes) {
        throw new Error("no default implementation");
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

class Face extends Transformable {
    constructor(surface, contour, holes = [], name, info) {
        super();
        this.surface = surface;
        this.contour = contour;
        this.holes = holes;
        this.name = name;
        this.info = info;
        this.aabb = undefined;
        //assert(name)
        Edge.assertLoop(contour);
        assert(contour.every((f) => f instanceof Edge), () => "contour.every(f => f instanceof Edge)" + contour);
        // contour.forEach(e => !surface.containsCurve(e.curve) &&
        // console.log('FAIL:'+surface.distanceToPoint(e.curve.anchor)))
        //contour.forEach(e => {
        //	assert(surface.containsCurve(e.curve), 'edge not in surface ' + e + surface)
        //})
        //assert(surface.edgeLoopCCW(contour), surface.toString() + contour.join('\n'))
        holes && holes.forEach((hole) => Edge.assertLoop(hole));
        holes && holes.forEach((hole) => assert(!surface.edgeLoopCCW(hole)));
        assert(!holes || holes.constructor == Array, holes && holes.toString());
        this.allEdges = Array.prototype.concat.apply(this.contour, this.holes);
    }
    static assembleFacesFromLoops(loops, surface, faceConstructor) {
        function placeRecursively(newLoopInfo, loopInfos) {
            if (loopInfos.length == 0) {
                loopInfos.push(newLoopInfo);
            }
            else {
                const subLoopInfo = loopInfos.find((loopInfo) => BRep.loop1ContainsLoop2(loopInfo.loop, loopInfo.ccw, newLoopInfo.loop, newLoopInfo.ccw, surface));
                if (subLoopInfo) {
                    placeRecursively(newLoopInfo, subLoopInfo.subloops);
                }
                else {
                    // newLoopInfo isnt contained by any other subLoopInfo
                    for (let i = loopInfos.length; --i >= 0;) {
                        const subLoopInfo = loopInfos[i];
                        //console.log('cheving subLoopInfo', surface.loopContainsPoint(newLoopInfo.edges,
                        // subLoopInfo.edges[0].a))
                        if (BRep.loop1ContainsLoop2(newLoopInfo.loop, newLoopInfo.ccw, subLoopInfo.loop, subLoopInfo.ccw, surface)) {
                            newLoopInfo.subloops.push(subLoopInfo);
                            loopInfos.splice(i, 1); // remove it
                        }
                    }
                    loopInfos.push(newLoopInfo);
                }
            }
        }
        function newFacesRecursive(loopInfo) {
            newFaces.push(new faceConstructor(surface, loopInfo.ccw ? loopInfo.loop : Edge.reversePath(loopInfo.loop), loopInfo.subloops.map((sl) => sl.ccw ? Edge.reversePath(sl.loop) : sl.loop)));
            loopInfo.subloops.forEach((sl) => sl.subloops.forEach((sl2) => newFacesRecursive(sl2)));
        }
        const newFaces = [];
        const topLevelLoops = [];
        loops.forEach((loop) => placeRecursively({
            loop: loop,
            ccw: surface.edgeLoopCCW(loop),
            subloops: [],
        }, topLevelLoops));
        topLevelLoops.forEach((tll) => newFacesRecursive(tll));
        return newFaces;
    }
    //fromLoops(loops: Edge[][], surface: Surface) {
    //	type LoopInfo = {loop: Edge[], ccw: boolean, subloops: LoopInfo[]}
    //	function placeRecursively(newLoopInfo: LoopInfo, loopInfos: LoopInfo[]) {
    //		if (loopInfos.length == 0) {
    //			loopInfos.push(newLoopInfo)
    //		} else {
    //			const subLoopInfo = loopInfos.find(loopInfo => BRep.loop1ContainsLoop2(loopInfo.loop, loopInfo.ccw,
    // newLoopInfo.loop, newLoopInfo.ccw, surface)) if (subLoopInfo) { placeRecursively(newLoopInfo,
    // subLoopInfo.subloops) } else { // newLoopInfo isnt contained by any other subLoopInfo for (let i =
    // loopInfos.length; --i >= 0;) { const subLoopInfo = loopInfos[i] //console.log('cheving subLoopInfo',
    // surface.loopContainsPoint(newLoopInfo.edges, subLoopInfo.edges[0].a)) if
    // (BRep.loop1ContainsLoop2(newLoopInfo.loop, subLoopInfo.loop, surface)) { newLoopInfo.subloops.push(subLoopInfo)
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
        return surface instanceof PlaneSurface
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
            return checkedPairs.has(new Pair(a, b));
        }
        function addPair(a, b) {
            return checkedPairs.add(new Pair(a, b));
        }
        /**
         * @param newEdge generated segment
         * @param col1 if newEdge is colinear to an edge of this, the edge in question
         * @param col2 same for face2
         * @return whether new edge was added.
         */
        function handleNewEdge(newEdge, col1, col2) {
            if (!col1 && !col2) {
                let correctDir = face.surface
                    .normalP(newEdge.a)
                    .cross(face2.surface.normalP(newEdge.a));
                if (correctDir.likeO()) {
                    const t = lerp(newEdge.aT, newEdge.bT, 1 / GOLDEN_RATIO), p = newEdge.curve.at(t);
                    correctDir = face.surface.normalP(p).cross(face2.surface.normalP(p));
                }
                if (!correctDir.likeO()) {
                    if (correctDir.dot(newEdge.aDir) < 0) {
                        newEdge = newEdge.flipped();
                    }
                    mapPush(faceMap, face, newEdge);
                    mapPush(faceMap, face2, newEdge.flipped());
                }
                else {
                    const p = newEdge.a;
                    const plane = P3.normalOnAnchor(newEdge.aDir, p);
                    const up = face.surface.normalP(p);
                    const sameDir = up.dot(face2.surface.normalP(p)) > 0;
                    const canonDir = plane.normal1.cross(up);
                    const curve = face.surface.isCurvesWithPlane(plane)[0], curveT = curve.pointT(p), curveDir = sign(canonDir.dot(curve.tangentAt(curveT)));
                    const curve2 = face2.surface.isCurvesWithPlane(plane)[0], curve2T = curve2.pointT(p), curve2Dir = sign(canonDir.dot(curve.tangentAt(curve2T)));
                    const foo = curve.diff(curveT, EPS * curveDir).dot(up);
                    const foo2 = curve2.diff(curve2T, EPS * curve2Dir).dot(up);
                    if (foo2 < foo) {
                        mapPush(faceMap, face2, sameDir ? newEdge.flipped() : newEdge);
                    }
                    if (up.dot(face2.surface.normalP(p)) < 0 == foo2 < foo) {
                        mapPush(faceMap, face, newEdge.flipped());
                    }
                    const bar = curve.diff(curveT, EPS * curveDir).dot(up);
                    const bar2 = curve2.diff(curve2T, EPS * curve2Dir).dot(up);
                    if (bar2 < bar) {
                        mapPush(faceMap, face2, sameDir ? newEdge : newEdge.flipped());
                    }
                    if (sameDir != bar2 < bar) {
                        mapPush(faceMap, face, newEdge);
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
                    thisBrep.edgeFaces.get(col1.getCanon()).forEach((faceInfo) => {
                        //const dot = snap0(surface2.normal1.dot(faceInfo.inside))
                        //if (dot == 0 ? !coplanarSameIsInside : dot < 0) {
                        const pointsInsideFace = fff(faceInfo, face2.surface);
                        const edgeInside = pointsInsideFace == INSIDE ||
                            (!coplanarSameIsInside && pointsInsideFace == COPLANAR_SAME);
                        const pushEdge = faceInfo.edge
                            .tangentAt(faceInfo.edge.curve.pointT(newEdge.a))
                            .like(newEdge.aDir)
                            ? newEdge
                            : newEdge.flipped();
                        console.log(newEdge.sce);
                        assert(faceInfo.edge
                            .tangentAt(faceInfo.edge.curve.pointT(pushEdge.a))
                            .like(pushEdge.aDir));
                        edgeInside && mapPush(faceMap, faceInfo.face, pushEdge);
                    });
                    const surface2NormalAtNewEdgeA = surface2.normalP(newEdge.a);
                    const newEdgeInside = surface2NormalAtNewEdgeA.cross(newEdge.aDir);
                    const sVEF1 = splitsVolumeEnclosingFacesP(thisBrep, col1.getCanon(), newEdge.a, newEdgeInside, surface2NormalAtNewEdgeA);
                    let addNewEdge, addNewEdgeFlipped;
                    if ((addNewEdge =
                        sVEF1 == INSIDE ||
                            (coplanarSameIsInside && sVEF1 == COPLANAR_SAME))) {
                        mapPush(faceMap, face2, newEdge);
                    }
                    const sVEF2 = splitsVolumeEnclosingFacesP(thisBrep, col1.getCanon(), newEdge.a, newEdgeInside.negated(), surface2NormalAtNewEdgeA);
                    if ((addNewEdgeFlipped =
                        sVEF2 == INSIDE ||
                            (coplanarSameIsInside && sVEF2 == COPLANAR_SAME))) {
                        mapPush(faceMap, face2, newEdge.flipped());
                    }
                    if (addNewEdge ||
                        addNewEdgeFlipped ||
                        (sVEF1 == COPLANAR_SAME && sVEF2 == INSIDE) ||
                        (sVEF2 == COPLANAR_SAME && sVEF1 == INSIDE)) {
                        return true;
                    }
                }
                return false;
            }
            const c1 = handleEdgeInFace(col1, col2, face, face2, thisBrep, face2Brep, false);
            const c2 = handleEdgeInFace(col2, col1, face2, face, face2Brep, thisBrep, true);
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
                        const edgeInside = sVEF == INSIDE || (coplanarSameIsInside && sVEF == COPLANAR_SAME);
                        const pushEdge = faceInfo.edge.aDir.like(newEdge.aDir)
                            ? newEdge
                            : newEdge.flipped();
                        if (edgeInside) {
                            mapPush(faceMap, faceInfo.face, pushEdge);
                            const aT = col1.getCanon().curve.pointT(newEdge.a);
                            if (!eq(aT, col1.aT) && !eq(aT, col1.bT)) {
                                // newEdge.a is in center of col1
                                if (splitsVolumeEnclosingCone2(face2Brep, newEdge.a, newEdge.curve, newEdge.aT, -Math.sign(newEdge.deltaT())) == INSIDE) {
                                    mapPush(thisEdgePoints, col1.getCanon(), {
                                        p: newEdge.a,
                                        edgeT: aT,
                                    });
                                }
                            }
                            const bT = col1.getCanon().curve.pointT(newEdge.b);
                            if (!eq(bT, col1.aT) && !eq(bT, col1.bT)) {
                                if (splitsVolumeEnclosingCone2(face2Brep, newEdge.b, newEdge.curve, newEdge.bT, Math.sign(newEdge.deltaT())) == INSIDE) {
                                    mapPush(thisEdgePoints, col1.getCanon(), {
                                        p: newEdge.b,
                                        edgeT: bT,
                                    });
                                }
                            }
                        }
                    }
                }
                handleColinearEdgeFaces(col1, col2, thisBrep, face2Brep, true, thisEdgePoints);
                handleColinearEdgeFaces(col2, col1, face2Brep, thisBrep, false, otherEdgePoints);
                return false;
            }
            return false;
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
                    mapPush(thisEdgePoints, a.edge.getCanon(), a);
                    assert(a.edge.isValidT(a.edgeT));
                }
                // else colinear segment ends in middle of other face, do nothing
            }
            // ends in the middle of a's face
            if (b && !a) {
                if (!b.colinear && b.edgeT != b.edge.aT && b.edgeT != b.edge.bT) {
                    mapPush(otherEdgePoints, b.edge.getCanon(), b);
                    assert(b.edge.isValidT(b.edgeT));
                }
                // else colinear segment ends in middle of other face, do nothing
            }
            if (a && b) {
                assert(a.colinear || b.colinear || eq(a.t, b.t));
                // if a or b is colinear the correct points will already have been added to the edge by handleNewEdge
                // segment starts/ends on edge/edge intersection
                function handleAB(a, b, face, face2, thisPlane, face2Plane, thisBrep, face2Brep, first, thisEdgePoints) {
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
                                mapPush(thisEdgePoints, a.edge.getCanon(), a);
                                assert(a.edge.isValidT(a.edgeT));
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
                            if (INSIDE == sVEF1 ||
                                (first && COPLANAR_SAME == sVEF1) ||
                                INSIDE == sVEF2 ||
                                (first && COPLANAR_SAME == sVEF2)) {
                                mapPush(thisEdgePoints, a.edge.getCanon(), a);
                                assert(a.edge.isValidT(a.edgeT));
                            }
                        }
                        //}
                    }
                }
                handleAB(a, b, face, face2, surface, surface2, thisBrep, face2Brep, true, thisEdgePoints);
                handleAB(b, a, face2, face, surface2, surface, face2Brep, thisBrep, false, otherEdgePoints);
            }
        }
        assertInst(Face, face2);
        const face = this;
        const surface = face.surface, surface2 = face2.surface;
        if (!this.getAABB().touchesAABBfuzzy(face2.getAABB())) {
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
            assert(surface.containsCurve(isCurve));
            assert(surface2.containsCurve(isCurve));
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
            assertf(() => 0 == ps1.length ||
                !eq0(ps1[0].insideDir.dot(isCurve.tangentAt(ps1[0].t))), () => ps1[0].insideDir.dot(isCurve.tangentAt(ps1[0].t)));
            assertf(() => 0 == ps2.length ||
                !eq0(ps2[0].insideDir.dot(isCurve.tangentAt(ps2[0].t))), () => ps2[0].insideDir.dot(isCurve.tangentAt(ps2[0].t)));
            function startsInside(ps, face) {
                if (0 == ps.length) {
                    return (isFinite(isCurve.tMin) &&
                        face.containsPoint2(isCurve.at(isCurve.tMin)) == PointVsFace.INSIDE);
                }
                else {
                    return ps[0].insideDir.dot(isCurve.tangentAt(ps[0].t)) < 0;
                }
            }
            // they can't both be empty currently
            // they can't both start 'inside'
            let in1 = startsInside(ps1, face);
            let in2 = startsInside(ps2, face2);
            if ((0 == ps1.length && !in1) || (0 == ps2.length && !in2)) {
                continue;
            }
            //assert(!in1 || !in2)
            let col1, col2;
            let i = 0, j = 0, last;
            let startP = in1 && in2 ? isCurve.at(isCurve.tMin) : undefined, startDir, startT = isCurve.tMin, startA, startB;
            while (i < ps1.length || j < ps2.length) {
                assert(i <= ps1.length);
                assert(j <= ps2.length);
                const a = ps1[i], b = ps2[j];
                assert(a || b);
                if (j == ps2.length || (i < ps1.length && lt(a.t, b.t))) {
                    last = a;
                    in1 = !in1;
                    a.used = true;
                    col1 = a.colinear ? a : undefined;
                    i++;
                }
                else if (i == ps1.length || gt(a.t, b.t)) {
                    last = b;
                    b.used = true;
                    in2 = !in2;
                    col2 = b.colinear ? b : undefined;
                    j++;
                }
                else {
                    last = a;
                    a.used = true;
                    b.used = true;
                    in1 = !in1;
                    in2 = !in2;
                    //if (in1 == in2) {
                    col1 = a.colinear ? a : undefined;
                    col2 = b.colinear ? b : undefined;
                    //}
                    i++;
                    j++;
                }
                if (startP && !(in1 && in2)) {
                    // segment end
                    startDir = isCurve.tangentAt(startT);
                    if (eq(startT, last.t)) {
                        startP = undefined;
                        continue;
                    }
                    assert(lt(startT, last.t));
                    startT > last.t && (startDir = startDir.negated());
                    let endDir = isCurve.tangentAt(last.t);
                    startT > last.t && (endDir = endDir.negated());
                    const newEdge = createEdge(isCurve, startP, last.p, startT, last.t, undefined, startDir, endDir, "genseg" + getGlobalId());
                    startP = undefined;
                    if (handleNewEdge(newEdge, col1 && col1.edge, col2 && col2.edge)) {
                        handleEndPoint(startA || col1, startB || col2);
                        handleEndPoint((a && a.used && a) || col1, (b && b.used && b) || col2);
                    }
                }
                else if (in1 && in2) {
                    // new segment just started
                    startP = last.p;
                    startDir = last.insideDir;
                    startT = last.t;
                    startA = a && a.used ? a : undefined;
                    startB = b && b.used ? b : undefined;
                }
            }
            if (in1 && in2 && startT !== isCurve.tMax) {
                const endT = isCurve.tMax;
                startDir = isCurve.tangentAt(startT);
                startT > endT && (startDir = startDir.negated());
                let endDir = isCurve.tangentAt(endT);
                startT > endT && (endDir = endDir.negated());
                const newEdge = createEdge(isCurve, startP, isCurve.at(endT), startT, endT, undefined, startDir, endDir, "genseg" + getGlobalId());
                if (handleNewEdge(newEdge, col1 && col1.edge, col2 && col2.edge)) {
                    handleEndPoint(startA || col1, startB || col2);
                }
            }
        }
        face.getAllEdges().forEach((edge) => {
            checkedPairs.add(new Pair(edge.getCanon(), face2));
        });
        face2.getAllEdges().forEach((edge) => {
            checkedPairs.add(new Pair(edge.getCanon(), face));
        });
    }
    edgeISPsWithSurface(isCurve, surface2) {
        const face = this;
        const surface = face.surface;
        const loops = face.holes.concat([face.contour]);
        const ps = [];
        for (const loop of loops) {
            const colinearEdges = loop.map((edge) => edge.curve.isColinearTo(isCurve));
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
                        if (!colinearEdges[prevEdgeIndex] &&
                            dotCurve2(prevEdge.curve, prevEdge.bT, colinearOutA, -sign(prevEdge.deltaT())) > 0) {
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
                        if (!colinearEdges[nextEdgeIndex] &&
                            dotCurve2(nextEdge.curve, nextEdge.aT, colinearOutB, sign(nextEdge.deltaT())) > 0) {
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
                        assert(!isNaN(curveT));
                        const insideDir = edge
                            .tangentAt(edgeT)
                            .cross(surface.normalP(p))
                            .negated();
                        const isTangent = isCurve.tangentAt(curveT);
                        //if(!eq0(insideDir.dot(isTangent))) {
                        // Edge.edgeISTsWithSurface returns snapped values, so comparison with == is ok:
                        if (edgeT == edge.bT) {
                            // endpoint lies on intersection line
                            if (!colinearEdges[nextEdgeIndex]) {
                                if (!eq(curveT, isCurve.tMax)) {
                                    const pointsToInside = this.pointsToInside3(edge.b, isCurve, curveT, 1);
                                    assert(pointsToInside != PointVsFace.ON_EDGE);
                                    if (PointVsFace.INSIDE == pointsToInside) {
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
                                if (!eq(curveT, isCurve.tMin)) {
                                    const pointsToInside = this.pointsToInside3(edge.b, isCurve, curveT, -1);
                                    assert(pointsToInside != PointVsFace.ON_EDGE);
                                    if (PointVsFace.INSIDE == pointsToInside) {
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
                            if (eq0(insideDir.dot(isTangent))) {
                                const dirFactor = sign(isTangent.dot(edge.curve.tangentAt(edgeT)));
                                const eps = 1e-4;
                                for (const dir of [-1, 1]) {
                                    if ((-1 == dir * dirFactor && edgeT == edge.minT) ||
                                        (1 == dir * dirFactor && edgeT == edge.maxT) ||
                                        (-1 == dir && curveT == isCurve.tMin) ||
                                        (1 == dir && curveT == isCurve.tMax))
                                        continue;
                                    const iscd = isCurve
                                        .at(curveT)
                                        .to(isCurve.at(curveT + dir * eps))
                                        .dot(insideDir);
                                    const ecd = edge.curve
                                        .at(edgeT)
                                        .to(edge.curve.at(edgeT + dir * dirFactor * eps))
                                        .dot(insideDir);
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
        const newEdges = Edge.reversePath(this.contour.map((e) => e.transform(m4)), mirroring);
        const newHoles = this.holes.map((hole) => Edge.reversePath(hole.map((e) => e.transform(m4)), mirroring));
        return new this.constructor(this.surface.transform(m4), newEdges, newHoles, this.name, this.info);
    }
    transform4(m4) {
        const mirroring = m4.isMirroring();
        const newEdges = Edge.reversePath(this.contour.map((e) => e.transform4(m4)), mirroring);
        const newHoles = this.holes.map((hole) => Edge.reversePath(hole.map((e) => e.transform4(m4)), mirroring));
        return new this.constructor(this.surface.transform4(m4), newEdges, newHoles, this.name, this.info);
    }
    flipped() {
        const newEdges = this.contour.map((e) => e.flipped()).reverse();
        const newHoles = this.holes.map((hole) => hole.map((e) => e.flipped()).reverse());
        return new this.constructor(this.surface.flipped(), newEdges, newHoles, this.name, this.info);
    }
    toString() {
        return ("new " +
            this.constructor.name +
            "(" +
            this.surface +
            ", [" +
            this.contour.map((e) => "\n\t" + e).join() +
            "]" +
            this.holes.map((hole) => "\n\t\thole: " + hole.join()) +
            ")");
    }
    toSource() {
        return ("new " +
            this.constructor.name +
            "(" +
            this.surface.toSource() +
            ", [" +
            this.contour.map((e) => "\n\t" + e.toSource() + ",").join("") +
            "], [" +
            this.holes
                .map((hole) => "[" + hole.map((e) => "\n\t" + e.toSource() + ",").join("") + "]")
                .join(",") +
            "])");
    }
    equals(obj) {
        return (this == obj ||
            (Object.getPrototypeOf(this) == Object.getPrototypeOf(obj) &&
                this.holes.length == obj.holes.length &&
                Edge.loopsEqual(this.contour, obj.contour) &&
                this.holes.every((hole) => obj.holes.some((hole2) => Edge.loopsEqual(hole, hole2)))));
    }
    hashCode() {
        function arrayHashCode(array) {
            let hashCode = 0;
            for (const val of array) {
                hashCode = (hashCode * 31 + val) | 0;
            }
            return hashCode;
        }
        function loopHashCode(loop) {
            return arrayHashCode(loop.map((edge) => edge.hashCode()).sort(MINUS));
        }
        let hashCode = 0;
        hashCode =
            (hashCode * 31 +
                arrayHashCode(this.holes.map((loop) => loopHashCode(loop)).sort(MINUS))) |
                0;
        hashCode = (hashCode * 31 + loopHashCode(this.contour)) | 0;
        hashCode = (hashCode * 31 + this.surface.hashCode()) | 0;
        return hashCode;
    }
    likeFace(face2) {
        function loopsLike(a, b) {
            return (a.length == b.length &&
                arrayRange(0, a.length, 1).some((offset) => a.every((edge, i) => edge.like(b[(offset + i) % a.length]))));
        }
        assertInst(Face, face2);
        return (this.surface.like(face2.surface) &&
            this.holes.length == face2.holes.length &&
            loopsLike(this.contour, face2.contour) &&
            this.holes.every((hole) => face2.holes.some((hole2) => loopsLike(hole, hole2))));
    }
    getAllEdges() {
        return this.allEdges;
    }
    addEdgeLines(mesh) {
        assert(false, "buggy, fix");
        const vertices = this.contour.flatMap((edge) => edge.getVerticesNo0()), mvl = mesh.vertices.length;
        for (let i = 0; i < vertices.length; i++) {
            mesh.vertices.push(vertices[i]);
            mesh.LINES.push(mvl + i, mvl + ((i + 1) % vertices.length));
        }
    }
    containsPoint(p) {
        assertVectors(p);
        return (this.surface.loopContainsPoint(this.contour, p) != PointVsFace.OUTSIDE &&
            !this.holes.some((hole) => this.surface.loopContainsPoint(hole, p) != PointVsFace.OUTSIDE));
    }
    containsPoint2(p) {
        assertVectors(p);
        const contourContainsPoint = this.surface.loopContainsPoint(this.contour, p);
        if (contourContainsPoint != PointVsFace.INSIDE)
            return contourContainsPoint;
        for (const hole of this.holes) {
            const loopContainsPoint = this.surface.loopContainsPoint(hole, p);
            if (loopContainsPoint != PointVsFace.OUTSIDE) {
                return loopContainsPoint == PointVsFace.ON_EDGE
                    ? PointVsFace.ON_EDGE
                    : PointVsFace.OUTSIDE;
            }
        }
        return PointVsFace.INSIDE;
    }
    /**
     *
     * @param line
     * @returns t param of the line if there is an intersection, NaN otherwise
     */
    intersectsLine(line) {
        assertInst(L3, line);
        if (!this.getAABB().intersectsLine(line))
            return NaN;
        const containedIntersectionsTs = this.surface
            .isTsForLine(line)
            .filter((t) => this.containsPoint(line.at(t)));
        const nearestPointT = containedIntersectionsTs.withMax((t) => -t);
        return undefined != nearestPointT ? nearestPointT : NaN;
    }
    toMesh() {
        const mesh = new Mesh()
            .addIndexBuffer("TRIANGLES")
            .addIndexBuffer("LINES")
            .addVertexBuffer("normals", "ts_Normal");
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
        return [this.contour, ...this.holes];
    }
    getAABB() {
        return (this.aabb ||
            (this.aabb = AABB.forAABBs(this.contour.map((e) => e.getAABB()))));
    }
    pointsToInside3(p, curve, curveT, dir) {
        const eps = 1e-6;
        const normal = this.surface.normalP(p);
        const curveTangent = curve.tangentAt(curveT).times(dir);
        const up = normal.cross(curveTangent);
        const ecd = curve
            .at(curveT)
            .to(curve.at(curveT + dir * eps))
            .dot(up);
        let minValue = Infinity, result, advanced = false;
        for (const edge of this.getAllEdges()) {
            const aEqP = edge.a.like(p), bEqP = edge.b.like(p);
            assert(aEqP == edge.a.like(p));
            assert(bEqP == edge.b.like(p));
            if (!aEqP && !bEqP)
                continue;
            const edgeTangent = aEqP ? edge.aDir : edge.bDir.negated();
            const angle = curveTangent.angleRelativeNormal(edgeTangent, normal);
            if (eq0(angle)) {
                if (curve.isColinearTo(edge.curve)) {
                    return PointVsFace.ON_EDGE;
                }
                const edgeT = aEqP ? edge.aT : edge.bT;
                const edgeDir = (aEqP ? 1 : -1) * sign(edge.deltaT());
                const iscd = edge.curve.diff(edgeT, edgeDir * eps).dot(up);
                //const iscd = edge.curve.at(edgeT).to(curve.at(edgeT + edgeDir * eps)).dot(up)
                const diff = iscd - ecd;
                if (diff > 0 && (!advanced || diff < minValue)) {
                    advanced = true;
                    minValue = diff;
                    result = aEqP ? PointVsFace.OUTSIDE : PointVsFace.INSIDE;
                }
            }
            else if (!advanced) {
                const angle2 = (angle + TAU) % TAU;
                if (angle2 < minValue) {
                    minValue = angle2;
                    result = aEqP ? PointVsFace.OUTSIDE : PointVsFace.INSIDE;
                }
            }
        }
        if (result == undefined)
            throw new Error();
        return result;
    }
    pointsToInside2(p, dir) {
        return this.pointsToInside3(p, L3.anchorDirection(p, dir), 0, 1);
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
        assert(p instanceof P3 || p instanceof PlaneSurface);
        super(p instanceof P3 ? new PlaneSurface(p) : p, contour, holes, name, info);
    }
    static forVertices(planeSurface, vs, ...holeVss) {
        const _planeSurface = planeSurface instanceof P3 ? new PlaneSurface(planeSurface) : planeSurface;
        assert(isCCW(vs, _planeSurface.plane.normal1), "isCCW(vs, planeSurface.plane.normal1)");
        const edges = StraightEdge.chain(vs);
        holeVss.forEach((vs) => assert(doubleSignedArea(vs, _planeSurface.plane.normal1) >= 0, "doubleSignedArea(vs, planeSurface.plane.normal1) >= 0"));
        const holes = holeVss.map((hvs) => StraightEdge.chain(hvs));
        return new PlaneFace(planeSurface, edges, holes);
    }
    addToMesh(mesh) {
        const mvl = mesh.vertices.length;
        const normal = this.surface.plane.normal1;
        const vertices = this.contour.flatMap((edge) => edge.getVerticesNo0());
        for (let i = 0; i < vertices.length; i++) {
            mesh.LINES.push(mvl + i, mvl + ((i + 1) % vertices.length));
        }
        const holeStarts = [];
        this.holes.forEach((hole) => {
            holeStarts.push(vertices.length);
            vertices.push(...hole.flatMap((edge) => edge.getVerticesNo0()));
        });
        const triangles = triangulateVertices(normal, vertices, holeStarts).map((index) => index + mvl);
        Array.prototype.push.apply(mesh.vertices, vertices);
        Array.prototype.push.apply(mesh.TRIANGLES, triangles);
        Array.prototype.push.apply(mesh.normals, arrayFromFunction(vertices.length, () => normal));
    }
    intersectsLine(line) {
        assertInst(L3, line);
        const lambda = line.isTWithPlane(this.surface.plane);
        if (!Number.isFinite(lambda)) {
            return NaN;
        }
        const inside = this.containsPoint(line.at(lambda));
        return inside ? lambda : NaN;
    }
    //intersectPlaneFace(face2: PlaneFace,
    //                   thisBrep: BRep,
    //                   face2Brep: BRep,
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
        return this.containsPoint2(p.plus(dir.times(NLA_PRECISION * 8)));
    }
    edgeISPsWithPlane(isLine, plane2) {
        const face = this;
        assert(face.surface.plane.containsLine(isLine));
        assert(plane2.containsLine(isLine));
        const plane = face.surface.plane;
        const ps = [];
        const loops = [face.contour].concat(face.holes);
        loops.forEach((loop) => {
            const colinearEdges = loop.map((edge) => edge.colinearToLine(isLine) && -sign(edge.aDir.dot(isLine.dir1)));
            const isLineOut = isLine.dir1.cross(plane.normal1);
            loop.forEach((edge, edgeIndex, edges) => {
                const nextEdgeIndex = (edgeIndex + 1) % edges.length, nextEdge = edges[nextEdgeIndex], colinearEdge = colinearEdges[edgeIndex];
                //console.log(edge.toSource()) {p:V3(2, -2.102, 0),
                if (colinearEdge) {
                    // edge colinear to intersection line
                    const curveAT = isLine.pointT(edge.a), curveBT = isLine.pointT(edge.b);
                    // add interval for colinear segment
                    ps.push({
                        p: edge.a,
                        insideDir: edge.aDir,
                        t: curveAT,
                        edge: edge,
                        edgeT: edge.aT,
                        colinear: true,
                    }, {
                        p: edge.b,
                        insideDir: edge.bDir.negated(),
                        t: curveBT,
                        edge: edge,
                        edgeT: edge.bT,
                        colinear: true,
                    });
                    // open next interval if necessary
                    const nextSide = colinearEdges[nextEdgeIndex] ||
                        dotCurve2(nextEdge.curve, nextEdge.aT, isLineOut, nextEdge.deltaTSign());
                    if (colinearEdge * nextSide < 0) {
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
                    assert(edgeTs.every((t) => plane2.containsPoint(edge.curve.at(t))), edgeTs);
                    for (const edgeT of edgeTs) {
                        if (edgeT == edge.bT) {
                            // endpoint lies on intersection line
                            const side = dotCurve2(edge.curve, edge.bT, isLineOut, -edge.deltaTSign());
                            const nextSide = colinearEdges[nextEdgeIndex] ||
                                dotCurve2(nextEdge.curve, nextEdge.aT, isLineOut, nextEdge.deltaTSign());
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
                            assert(plane2.containsPoint(p), edge.toString(), p, edgeT, plane2.distanceToPoint(p));
                            assert(isLine.containsPoint(p), edge.toString(), p, edgeT, isLine.distanceToPoint(p));
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
            return false;
        }
        for (const edge of loop) {
            const ts = edge.edgeISTsWithPlane(seamPlane);
            if (ts.length == 0) {
                if (!(edge.curve instanceof L3) &&
                    checkSide(seamPlane.distanceToPointSigned(edge.a)))
                    return false;
            }
            else {
                for (const t of ts) {
                    // TODO: this part probably should be in a separate function
                    // check 'backwards' only if if aT != t
                    if (edge.aT != t) {
                        if (checkSide(dotCurve2(edge.curve, t, seamPlane.normal1, -edge.deltaTSign())))
                            return false;
                    }
                    if (edge.bT != t) {
                        if (checkSide(dotCurve2(edge.curve, t, seamPlane.normal1, edge.deltaTSign())))
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
        this.aabb = AABB.forAABBs(this.contour.map((e) => e.getAABB()));
        this.aabb.addPoints(this.surface.getExtremePoints().filter((p) => this.containsPoint(p)));
        return this.aabb;
    }
    unrollLoop(edgeLoop) {
        const vs = [];
        const uvP = this.surface.uvPFunc();
        const verticesNo0s = edgeLoop.map((edge) => edge.getVerticesNo0());
        const startEdgeIndex = verticesNo0s.findIndex((edgeVertices) => !eq(uvP(edgeVertices[0]).x, Math.PI));
        assert(-1 != startEdgeIndex);
        // console.log(startEdgeIndex)
        for (let i = 0; i < edgeLoop.length; i++) {
            const edgeIndex = (i + startEdgeIndex) % edgeLoop.length;
            for (let j = 0; j < verticesNo0s[edgeIndex].length; j++) {
                const p = verticesNo0s[edgeIndex][j];
                const localP = uvP(p);
                // console.log(hint, p.sce, localP.sce)
                vs.push(localP);
            }
        }
        edgeLoop.forEach((edge) => {
            edge.getVerticesNo0().forEach((p) => {
                vs.push(uvP(p));
            });
        });
        console.log("vs\n", vs.join("\n"), vs.length);
        return vs;
    }
    /**
     * f1 cos t + f2 sin t
     * tan(phi) = sin / cos
     *          = (f1x cos t + f2x sin t) / (f1y cos t + f2y sin t)
     *
     *          = (-f1x sin t + f2x cos t) / (-f1y sin t + f2y cos t)
     */
    unrollEllipsoidLoops(edgeLoops) {
        const verticesUV = [], vertices = [], loopStarts = [];
        const ellipsoid = this.surface;
        const ptpf = ellipsoid.uvPFunc();
        const testDegeneratePoint = ellipsoid instanceof EllipsoidSurface
            ? (nextStart) => nextStart.like(ellipsoid.center.plus(ellipsoid.f3)) ||
                nextStart.like(ellipsoid.center.minus(ellipsoid.f3))
            : (nextStart) => nextStart.like(this.surface.center);
        for (const edgeLoop of edgeLoops) {
            loopStarts.push(verticesUV.length);
            // console.log(startEdgeIndex)
            for (let i = 0; i < edgeLoop.length; i++) {
                const ipp = (i + 1) % edgeLoop.length;
                const verticesNo0 = edgeLoop[i].getVerticesNo0();
                vertices.push(...verticesNo0);
                verticesUV.push(...verticesNo0.map((v) => ptpf(v)));
                const nextStart = edgeLoop[ipp].a;
                //console.log('BLAH', nextStart.str, ellipsoid.center.plus(ellipsoid.f3).str)
                if (testDegeneratePoint(nextStart)) {
                    const bDirLC = ellipsoid.matrixInverse.transformVector(edgeLoop[i].bDir), aDirLC = ellipsoid.matrixInverse.transformVector(edgeLoop[ipp].aDir);
                    const inAngle = Math.atan2(-bDirLC.y, -bDirLC.x);
                    const outAngle = Math.atan2(aDirLC.y, aDirLC.x);
                    const stLast = verticesUV.pop();
                    verticesUV.push(new V3(inAngle, stLast.y, 0), new V3(outAngle, stLast.y, 0));
                    vertices.push(vertices.last);
                }
                verticesUV.forEach(({ u, v }) => {
                    assert(isFinite(u));
                    assert(isFinite(v));
                });
            }
        }
        let normals;
        if (this.surface instanceof EllipsoidSurface) {
            normals = vertices.map((v) => ellipsoid.normalP(v));
        }
        else {
            const normalUV = ellipsoid.normalUVFunc();
            normals = verticesUV.map(({ u, v }) => normalUV(u, v));
        }
        assert(vertices.length == vertices.length);
        //console.log(verticesUV.map(v => v.str).join('\n'))
        return {
            verticesUV: verticesUV,
            vertices: vertices,
            normals: normals,
            loopStarts: loopStarts,
        };
    }
    unrollCylinderLoops(loops) {
        const vertexLoops = loops.map((loop) => loop.flatMap((edge) => edge.getVerticesNo0()));
        const surface = this.surface;
        const vertices = concatenated(vertexLoops);
        // this.unrollLoop(loop).map(v => new V3(v.x / uStep, v.y / vStep, 0)))
        const loopStarts = vertexLoops.reduce((arr, loop) => (arr.push(arr.last + loop.length), arr), [0]);
        const uvPFunc = surface.uvPFunc();
        const verticesUV = vertices.map((v) => uvPFunc(v));
        const uvN = surface.normalUVFunc();
        const normals = verticesUV.map(({ u, v }) => uvN(u, v));
        return {
            verticesUV: verticesUV,
            vertices: vertices,
            normals: normals,
            loopStarts: loopStarts,
        };
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
        assertf(() => uStep > 0 && vStep > 0, uStep, vStep, "Surface: " + this.surface);
        const triangles = [];
        const pMN = (m, n) => this.surface.pUVFunc()(m * uStep, n * vStep);
        const normalMN = (m, n) => this.surface.normalUVFunc()(m * uStep, n * vStep);
        const loops = this.getLoops();
        const { vertices, verticesUV, normals, loopStarts } = this.surface instanceof EllipsoidSurface ||
            this.surface instanceof ConicSurface
            ? this.unrollEllipsoidLoops(loops)
            : this.unrollCylinderLoops(loops);
        loopStarts.push(vertices.length);
        const verticesMN = verticesUV.map(({ u, v }) => new V3(u / uStep, v / vStep, 0));
        for (let vertexLoopIndex = 0; vertexLoopIndex < loops.length; vertexLoopIndex++) {
            const vertexLoopStart = loopStarts[vertexLoopIndex];
            const vertexLoopLength = loopStarts[vertexLoopIndex + 1] - vertexLoopStart;
            const base = mesh.vertices.length + loopStarts[vertexLoopIndex];
            for (let i = 0; i < vertexLoopLength; i++) {
                mesh.LINES.push(base + i, base + ((i + 1) % vertexLoopLength));
            }
        }
        disableConsole();
        let minM = Infinity, maxM = -Infinity, minN = Infinity, maxN = -Infinity;
        //console.log('surface', this.surface.str)
        //console.log(verticesMN)
        //drPs.push(...verticesMN.map((v, i) => ({p: vertices[i], text: `${i} uv: ${v.toString(x => round10(x,
        // -4))}`})))
        verticesMN.forEach(([m, n]) => {
            assert(isFinite(m));
            assert(isFinite(n));
            minM = min(minM, m);
            maxM = max(maxM, m);
            minN = min(minN, n);
            maxN = max(maxN, n);
        });
        if (ParametricSurface.is(this.surface)) ;
        const mOffset = floor(minM + NLA_PRECISION), nOffset = floor(minN + NLA_PRECISION);
        const mRes = ceil(maxM - NLA_PRECISION) - mOffset, nRes = ceil(maxN - NLA_PRECISION) - nOffset;
        console.log(uStep, vStep, mRes, nRes);
        if (mRes == 1 && nRes == 1) {
            // triangulate this face as if it were a plane
            const polyTriangles = triangulateVertices(V3.Z, verticesMN, loopStarts.slice(1, 1 + this.holes.length));
            triangles.push(...polyTriangles);
        }
        else {
            const partss = new Array(mRes * nRes);
            function fixUpPart(part, baseM, baseN) {
                assert(baseM < mRes && baseN < nRes, `${baseM}, ${baseN}, ${mRes}, ${nRes}`);
                console.log("complete part", part, baseM, baseN);
                //console.trace()
                assert(part.length);
                const cellM = baseM + mOffset, cellN = baseN + nOffset;
                for (const index of part) {
                    assert(le(cellM, verticesMN[index].x) &&
                        le(verticesMN[index].x, cellM + 1), `${index} ${verticesMN[index].str} ${cellM} ${cellM}`);
                    assert(le(cellN, verticesMN[index].y) &&
                        le(verticesMN[index].y, cellN + 1));
                }
                const pos = baseN * mRes + baseM;
                (partss[pos] || (partss[pos] = [])).push(part);
                //const outline = partss[pos] || (partss[pos] = [minM + baseM * uStep, minN + baseN * vStep, minM +
                // (baseM + 1) * uStep, minN + (baseN + 1) * vStep])
            }
            // 'some' instead of forEach so we can return out of the entire function if this.edges crosses no borders
            // and
            for (let vertexLoopIndex = 0; vertexLoopIndex < loops.length; vertexLoopIndex++) {
                let part = undefined;
                let firstPart = undefined;
                let firstPartBaseM = -1;
                let firstPartBaseN = -1;
                let lastBaseM = -1, lastBaseN = -1;
                let partCount = 0;
                const vertexLoopStart = loopStarts[vertexLoopIndex];
                const vertexLoopLength = loopStarts[vertexLoopIndex + 1] - vertexLoopStart;
                for (let vlvi = 0; vlvi < vertexLoopLength; vlvi++) {
                    const vx0index = vertexLoopStart + vlvi, vx0 = verticesMN[vx0index];
                    const vx1index = vertexLoopStart + ((vlvi + 1) % vertexLoopLength), vx1 = verticesMN[vx1index];
                    //console.log('dask', vx0index, vx1index)
                    const vx01 = vx0.to(vx1);
                    assert(vx0);
                    const di = vx01.x, dj = vx01.y;
                    let vxIndex = vx0index, vx = vx0, currentT = 0;
                    let whileLimit = 400;
                    while (--whileLimit) {
                        // points which are on a grid line are assigned to the cell into which they are going (+
                        // NLA_PRECISION * sign(di)) if they are parallel to the gridline (eq0(di)), they belong the
                        // the cell for which they are a CCW boundary
                        const baseM = floor(vx.u + (!eq0(di) ? sign(di) : -sign(dj)) * NLA_PRECISION) -
                            mOffset;
                        const baseN = floor(vx.v + (!eq0(dj) ? sign(dj) : sign(di)) * NLA_PRECISION) -
                            nOffset;
                        assert(baseM < mRes && baseN < nRes, `${baseM}, ${baseN}, ${mRes}, ${nRes}`);
                        // figure out the next intersection with a gridline:
                        // iNext is the positive horizontal distance to the next vertical gridline
                        const iNext = ceil(sign(di) * vx.u + NLA_PRECISION) - sign(di) * vx.u;
                        const jNext = ceil(sign(dj) * vx.v + NLA_PRECISION) - sign(dj) * vx.v;
                        const iNextT = currentT + iNext / abs(di);
                        const jNextT = currentT + jNext / abs(dj);
                        //console.log(vxIndex, vx.str, 'vij', vx.u, vx.v, 'd', di, dj, 'ijNext', iNext, jNext, 'nextT',
                        // iNextT, jNextT)
                        if (lastBaseM != baseM || lastBaseN != baseN) {
                            if (part) {
                                if (!firstPart) {
                                    firstPart = part;
                                    firstPartBaseM = lastBaseM;
                                    firstPartBaseN = lastBaseN;
                                }
                                else {
                                    partCount++;
                                    fixUpPart(part, lastBaseM, lastBaseN);
                                }
                            }
                            part = [vxIndex];
                        }
                        lastBaseM = baseM;
                        lastBaseN = baseN;
                        currentT = min(iNextT, jNextT);
                        if (ge(currentT, 1)) {
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
                    assert(whileLimit, "whileLimit");
                }
                if (0 == partCount) {
                    // complete loop
                    assert(false, "found a hole, try increasing resolution");
                }
                // at this point, the firstPart hasn't been added, and the last part also hasn't been added
                // either they belong to the same cell, or not
                if (firstPartBaseM == lastBaseM && firstPartBaseN == lastBaseN) {
                    part.pop();
                    fixUpPart(part.concat(firstPart), lastBaseM, lastBaseN);
                }
                else {
                    fixUpPart(firstPart, firstPartBaseM, firstPartBaseN);
                    fixUpPart(part, lastBaseM, lastBaseN);
                }
                console.log("firstPart", firstPart);
            }
            console.log("calculated parts", partss);
            const fieldVertexIndices = new Array((mRes + 1) * (nRes + 1));
            function addVertex(m, n) {
                verticesMN.push(new V3(m, n, 0));
                normals.push(normalMN(m, n));
                return vertices.push(pMN(m, n)) - 1;
            }
            function getGridVertexIndex(i, j) {
                const index = j * (mRes + 1) + i;
                return (fieldVertexIndices[index] ||
                    (fieldVertexIndices[index] = addVertex(i + mOffset, j + nOffset)));
            }
            for (let col = 0; col < mRes; col++) {
                let inside = false;
                for (let row = 0; row < nRes; row++) {
                    const pos = row * mRes + col;
                    const fieldU = mOffset + col, fieldV = nOffset + row;
                    const parts = partss[pos];
                    if (!parts) {
                        if (inside) {
                            pushQuad(triangles, false, getGridVertexIndex(col, row), getGridVertexIndex(col + 1, row), getGridVertexIndex(col, row + 1), getGridVertexIndex(col + 1, row + 1));
                        }
                    }
                    else {
                        // assemble the field with segments in in
                        function opos(index) {
                            const p = verticesMN[index], u1 = p.x - fieldU, v1 = p.y - fieldV;
                            assert(-NLA_PRECISION < u1 &&
                                u1 < 1 + NLA_PRECISION &&
                                -NLA_PRECISION < v1 &&
                                v1 < 1 + NLA_PRECISION, "oob u1 v1 " +
                                u1 +
                                " " +
                                v1 +
                                " " +
                                index +
                                " " +
                                p.str +
                                "IF THIS FAILS check canonSeamU is correct");
                            return v1 < u1 ? u1 + v1 : 4 - u1 - v1;
                        }
                        while (parts.length) {
                            const outline = [];
                            const startPart = parts[0];
                            assert(startPart.length > 0);
                            let currentPart = startPart;
                            do {
                                outline.push(...currentPart);
                                const currentPartEndOpos = opos(currentPart.last);
                                const nextPartIndex = indexWithMax(parts, (part) => -mod(opos(part[0]) - currentPartEndOpos, 4));
                                const nextPart = bagRemoveIndex(parts, nextPartIndex);
                                let currentOpos = currentPartEndOpos;
                                const nextPartStartOpos = opos(nextPart[0]) > currentOpos
                                    ? opos(nextPart[0])
                                    : opos(nextPart[0]) + 4;
                                let nextOpos = ceil(currentOpos + NLA_PRECISION);
                                let flipping = eq0(((currentOpos + NLA_PRECISION) % 1) - NLA_PRECISION);
                                //inside = inside != (!eq0(currentOpos % 1) && currentOpos % 2 < 1)
                                while (lt(nextOpos, nextPartStartOpos)) {
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
                                inside =
                                    inside !=
                                        (flipping &&
                                            nextOpos % 2 == 1 &&
                                            eq(nextOpos, nextPartStartOpos));
                                currentOpos = nextOpos;
                                currentPart = nextPart;
                            } while (currentPart != startPart);
                            // triangulate outline
                            if (outline.length == 3) {
                                // its just a triangle
                                triangles.push(...outline);
                            }
                            else {
                                const polyTriangles = triangulateVertices(V3.Z, outline.map((i) => verticesMN[i]), []).map((i) => outline[i]);
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
        Array.prototype.push.apply(mesh.TRIANGLES, triangles.map((index) => index + mesh.vertices.length));
        Array.prototype.push.apply(mesh.vertices, vertices);
        Array.prototype.push.apply(mesh.normals, normals);
        //this.addEdgeLines(mesh)
        enableConsole();
    }
    addToMesh2(mesh) {
        const zSplit = 8;
        const ribs = [];
        let minZ = Infinity, maxZ = -Infinity;
        //let cmp = (a, b) => a.value - b.value
        const f = this.surface.pUVFunc();
        const normalF = this.surface.normalUVFunc();
        const vertexLoops = this.holes
            .concat([this.contour])
            .map((loop) => this.unrollLoop(loop));
        vertexLoops.forEach((vertexLoop) => {
            vertexLoop.forEach(({ x: d, y: z }) => {
                const index0 = binaryIndexOf(ribs, d, (a, b) => snap(a.value - b, 0));
                if (index0 < 0) {
                    ribs.splice(-index0 - 1, 0, { value: d, left: [], right: [] });
                }
                minZ = min(minZ, z);
                maxZ = max(maxZ, z);
            });
        });
        console.log("zzzs", minZ, maxZ, vertexLoops[0].toSource().replace(/\), /g, ",\n"));
        const correction = 1;
        vertexLoops.forEach((vertexLoop) => {
            vertexLoop.forEach((v0, i, vs) => {
                let v1 = vs[(i + 1) % vs.length], dDiff = v1.x - v0.x;
                //console.log(v0.sce, v1.sce)
                if (eq0(dDiff)) {
                    return;
                }
                if (dDiff < 0) {
                    [v0, v1] = [v1, v0];
                    dDiff = -dDiff;
                }
                const index0 = binaryIndexOf(ribs, v0.x, (a, b) => snap(a.value - b, 0));
                const index1 = binaryIndexOf(ribs, v1.x, (a, b) => snap(a.value - b, 0));
                binaryInsert(ribs[index0].right, v0.y);
                for (let j = (index0 + correction) % ribs.length; j != index1; j = (j + correction) % ribs.length) {
                    const x = ribs[j].value;
                    const part = (x - v0.x) / dDiff;
                    const interpolated = v1.y * part + v0.y * (1 - part);
                    binaryInsert(ribs[j].left, interpolated);
                    binaryInsert(ribs[j].right, interpolated);
                }
                binaryInsert(ribs[index1].left, v1.y);
                // console.log(ribs.map(r=>r.toSource()).join('\n'))
            });
        });
        const vertices = [], triangles0 = [], normals = [];
        for (let i = 0; i < ribs.length; i++) {
            const ribLeft = ribs[i], ribRight = ribs[(i + 1) % ribs.length];
            assert(ribLeft.right.length == ribRight.left.length);
            for (let j = 0; j < ribLeft.right.length; j++) {
                vertices.push(f(ribLeft.value, ribLeft.right[j]), f(ribRight.value, ribRight.left[j]));
                normals.push(normalF(ribLeft.value, ribLeft.right[j]), normalF(ribRight.value, ribRight.left[j]));
            }
        }
        //console.log(ribs.map(r=>r.toSource()).join('\n'))
        const vss = vertices.length, detailVerticesStart = vss;
        const zInterval = maxZ - minZ, zStep = zInterval / zSplit;
        const detailZs = arrayFromFunction(zSplit - 1, (i) => minZ + (1 + i) * zStep);
        console.log("detailsZs", detailZs);
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
        const end =  ribs.length - 1;
        for (let i = 0; i < end; i++) {
            const ipp = (i + 1) % ribs.length;
            let inside = false, colPos = 0;
            const ribLeft = ribs[i], ribRight = ribs[(i + 1) % ribs.length];
            for (let j = 0; j < detailZs.length + 1; j++) {
                const detailZ = detailZs[j] || 100000;
                if (!inside) {
                    if (ribLeft.right[colPos] < detailZ &&
                        ribRight.left[colPos] < detailZ) {
                        if (ribLeft.right[colPos + 1] < detailZ ||
                            ribRight.left[colPos + 1] < detailZ) {
                            pushQuad(triangles0, flipped2, vsStart + colPos * 2, vsStart + (colPos + 1) * 2, vsStart + colPos * 2 + 1, vsStart + (colPos + 1) * 2 + 1);
                            colPos += 2;
                            if (ribLeft.right[colPos] < detailZ ||
                                ribRight.left[colPos] < detailZ) {
                                j--;
                            }
                        }
                        else {
                            pushQuad(triangles0, flipped2, vsStart + colPos * 2, vsStart + colPos * 2 + 1, detailVerticesStart + i * detailZs.length + j, detailVerticesStart + ipp * detailZs.length + j);
                            inside = true;
                            colPos++;
                        }
                    }
                }
                else {
                    if (ribLeft.right[colPos] < detailZ ||
                        ribRight.left[colPos] < detailZ) {
                        pushQuad(triangles0, flipped2, detailVerticesStart + i * detailZs.length + j - 1, detailVerticesStart + ipp * detailZs.length + j - 1, vsStart + colPos * 2, vsStart + colPos * 2 + 1);
                        inside = false;
                        colPos++;
                        if (ribLeft.right[colPos] < detailZ ||
                            ribRight.left[colPos] < detailZ) {
                            j--;
                        }
                    }
                    else {
                        pushQuad(triangles0, flipped2, detailVerticesStart + i * detailZs.length + j, detailVerticesStart + i * detailZs.length + j - 1, detailVerticesStart + ipp * detailZs.length + j, detailVerticesStart + ipp * detailZs.length + j - 1);
                    }
                }
            }
            vsStart += ribLeft.right.length * 2;
        }
        //console.log('trinagle', triangles0.max(), vertices.length, triangles0.length, triangles0.toSource(),
        // triangles0.map(i => vertices[i].$).toSource() )
        const triangles = triangles0.map((index) => index + mesh.vertices.length);
        //assert(normals.every(n => n.hasLength(1)), normals.find(n => !n.hasLength(1)).length() +'
        // '+normals.findIndex(n => !n.hasLength(1)))
        Array.prototype.push.apply(mesh.vertices, vertices);
        Array.prototype.push.apply(mesh.TRIANGLES, triangles);
        Array.prototype.push.apply(mesh.normals, normals);
        //this.addEdgeLines(mesh)
    }
}

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
function assembleFaceFromLooseEdges(edges, surface, originalFace) {
    const visited = new Set();
    function nextStart() {
        return edges.find((edge) => !visited.has(edge));
    }
    const loops = [];
    let startEdge, currentEdge = undefined;
    while ((startEdge = nextStart())) {
        currentEdge = startEdge;
        const loop = [];
        let total = 0;
        do {
            visited.add(currentEdge);
            loop.push(currentEdge);
            const possibleEdges = edges.filter((edge) => currentEdge.b.like(edge.a));
            const normalAtCurrentB = surface.normalP(currentEdge.b);
            const nextEdgeIndex = indexWithMax(possibleEdges, (edge) => currentEdge.bDir.angleRelativeNormal(edge.aDir, normalAtCurrentB));
            currentEdge = possibleEdges[nextEdgeIndex];
        } while (startEdge != currentEdge && total++ < 200);
        assert(total != 201);
        loops.push(loop);
    }
    const assembledFaces = BRep.assembleFacesFromLoops(loops, surface, originalFace);
    assertf(() => 1 == assembledFaces.length);
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
    const dir = sign(currentEdge.deltaT());
    const ecd = currentEdge.curve.diff(currentEdge.bT, -dir * eps).dot(normVector);
    for (let i = possibleEdges.length; i--;) {
        const edge = possibleEdges[i];
        const angle1 = currentEdge.bDir
            .negated()
            .angleRelativeNormal(edge.aDir, faceNormalAtCurrentB);
        const angle = ((angle1 + TAU + NLA_PRECISION) % TAU) - NLA_PRECISION;
        if (eq0(angle)) {
            // do advanced analysis
            if (currentEdge.curve.isColinearTo(edge.curve)) {
                continue;
            }
            const edgeDir = sign(edge.deltaT());
            const iscd = edge.curve.diff(edge.aT, edgeDir * eps).dot(normVector);
            const diff = iscd - ecd;
            // if diff > 0, the angle is actually ~= 0
            if (diff < 0 && (!advanced || diff > maxValue)) {
                advanced = true;
                maxValue = diff;
                result = i;
            }
        }
        else if (!advanced) {
            if (gt(angle, maxValue)) {
                maxValue = angle;
                result = i;
            }
        }
    }
    return result == Number.MAX_SAFE_INTEGER ? 0 : result;
}
class BRep extends Transformable {
    constructor(faces, infiniteVolume, generator, vertexNames) {
        super();
        this.faces = faces;
        assertInst(Face, ...faces);
        this.infiniteVolume = infiniteVolume;
        assert(!this.infiniteVolume || true === this.infiniteVolume);
        this.generator = generator;
        this.vertexNames = vertexNames;
        this.edgeFaces = undefined;
        //this.assertSanity()
    }
    static loop1ContainsLoop2(loop1, ccw1, loop2, ccw2, surface) {
        for (const edge of loop2) {
            const loop1ContainsPoint = surface.loopContainsPoint(loop1, edge.a);
            if (PointVsFace.ON_EDGE != loop1ContainsPoint)
                return PointVsFace.INSIDE == loop1ContainsPoint;
        }
        for (const edge of loop2) {
            const edgePoint = edge.curve.at(edge.aT * 0.2 + edge.bT * 0.8);
            const loop1ContainsPoint = surface.loopContainsPoint(loop1, edgePoint);
            if (PointVsFace.ON_EDGE != loop1ContainsPoint)
                return PointVsFace.INSIDE == loop1ContainsPoint;
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
                const subLoopInfo = loopInfos.find((loopInfo) => BRep.loop1ContainsLoop2(loopInfo.loop, loopInfo.ccw, newLoopInfo.loop, newLoopInfo.ccw, surface));
                if (subLoopInfo) {
                    placeRecursively(newLoopInfo, subLoopInfo.subloops);
                }
                else {
                    // newLoopInfo isnt contained by any other subLoopInfo
                    for (let i = loopInfos.length; --i >= 0;) {
                        const subLoopInfo = loopInfos[i];
                        //console.log("cheving subLoopInfo", surface.loopContainsPoint(newLoopInfo.edges,
                        // subLoopInfo.edges[0].a))
                        if (BRep.loop1ContainsLoop2(newLoopInfo.loop, newLoopInfo.ccw, subLoopInfo.loop, subLoopInfo.ccw, surface)) {
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
                if (loopInfo.subloops.every((sl) => !sl.ccw)) {
                    const holes = loopInfo.subloops.map((sl) => sl.loop);
                    const info = infoFactory &&
                        infoFactory.newSubFace(originalFace, surface, loopInfo.loop, holes);
                    const newFace = new originalFace.constructor(surface, loopInfo.loop, holes, "genface" + getGlobalId(), info);
                    newFaces.push(newFace);
                    loopInfo.subloops.forEach((sl) => sl.subloops.forEach((slsl) => slsl.ccw && newFacesRecursive(slsl)));
                }
                else {
                    loopInfo.subloops.forEach((sl) => sl.ccw && newFacesRecursive(sl));
                }
            }
        }
        const newFaces = [];
        const topLevelLoops = [];
        loops.forEach((loop) => placeRecursively({
            loop: loop,
            ccw: surface.edgeLoopCCW(loop),
            subloops: [],
        }, topLevelLoops));
        topLevelLoops.forEach((tll) => newFacesRecursive(tll));
        return newFaces;
    }
    /**
     * Create a [BRep] by concatenating the faces of other BReps. Only use this if certain that the faces of the BReps do not intersect.
     * Otherwise, use [BRep.plus].
     * @param bReps
     * @param generator
     */
    static join(bReps, generator) {
        return new BRep(bReps.flatMap((b2) => b2.faces), false, generator);
    }
    containsPoint(p, forceInsideOutside = false) {
        const dirs = [
            V(-0.3920414696448526, -0.12936136783391444, -0.9108068525164064),
            V(0.6520650903544943, -0.07151288645511984, -0.7547827667692488),
            V(0.9433494201061395, -0.2402757256238473, -0.22882186797013926),
            V(0.13678704228501923, -0.04480387361087783, 0.9895867410047372),
            V(0.0662057922721913, -0.5865836917435423, 0.8071780259955845),
            V(-0.7322576567870621, -0.12953393611526787, 0.6685953061989045),
            V(0.6579719127258273, -0.012300218400456116, 0.7529420075219719),
            V(-0.5576497966736425, 0.8006695748324647, 0.2189861552871446),
        ];
        dirLoop: for (const dir of dirs) {
            const testLine = new L3(p, dir);
            let inside = this.infiniteVolume;
            for (const face of this.faces) {
                assert(!face.surface.containsCurve(testLine));
                const ists = face.surface.isTsForLine(testLine);
                for (const t of ists) {
                    const p = testLine.at(t);
                    const pvf = face.containsPoint2(p);
                    //assert(pvf != PointVsFace.ON_EDGE)
                    !forceInsideOutside && assert(!eq0(t));
                    if (t > 0) {
                        if (pvf == PointVsFace.ON_EDGE) {
                            continue dirLoop;
                        }
                        if (pvf == PointVsFace.INSIDE) {
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
                    const faceGroup = likeSurfaceFaces.find((faceGroup) => faceGroup.includes(this.faces[j]));
                    if (faceGroup) {
                        faceGroup.push(this.faces[i]);
                        addedToGroup = true;
                    }
                }
            }
            !addedToGroup && likeSurfaceFaces.push([this.faces[i]]);
        }
        console.log("likeSurfaceFaces", likeSurfaceFaces);
        if (likeSurfaceFaces.every((group) => group.length == 1))
            return this;
        const newFaces = [];
        let total = 0;
        for (const faceGroup of likeSurfaceFaces) {
            console.log(faceGroup);
            if (faceGroup.length == 1) {
                newFaces.push(faceGroup[0]);
            }
            else {
                const allEdges = faceGroup.flatMap((face) => face.getAllEdges());
                for (let i = allEdges.length; i-- > 0;) {
                    for (let j = 0; j < i; j++) {
                        console.log("blugh", total);
                        assert(i >= 0 && j >= 0 && total++ < 500, i + " " + j + " " + total);
                        if (allEdges[i].isCoEdge(allEdges[j])) {
                            // remove both
                            allEdges.splice(i, 1);
                            allEdges.splice(j, 1);
                            i--;
                            break;
                        }
                    }
                }
                const newFace = assembleFaceFromLooseEdges(allEdges, faceGroup[0].surface, faceGroup[0]);
                newFaces.push(newFace);
            }
        }
        return new BRep(newFaces, this.infiniteVolume, this.generator && this.generator + ".withMergedFaces()", this.vertexNames);
    }
    calculateVolume() {
        return sum(this.faces.map((face) => face.zDirVolume().volume));
    }
    toMesh() {
        const mesh = new Mesh()
            .addVertexBuffer("normals", "ts_Normal")
            .addIndexBuffer("TRIANGLES")
            .addIndexBuffer("LINES");
        mesh.faceIndexes = new Map();
        for (const face of this.faces) {
            const triangleStart = mesh.TRIANGLES.length;
            face.addToMesh(mesh);
            mesh.faceIndexes.set(face, {
                start: triangleStart,
                count: mesh.TRIANGLES.length - triangleStart,
            });
        }
        //this.buildAdjacencies()
        //for (const edge of this.edgeFaces.keys()) {
        //
        //}
        return mesh;
    }
    minus(other, infoFactory) {
        const generator = this.generator &&
            other.generator &&
            this.generator + ".minus(" + other.generator + ")";
        return this.intersection(other.flipped(), true, true, generator, infoFactory);
    }
    plus(other, infoFactory) {
        const generator = this.generator &&
            other.generator &&
            this.generator + ".plus(" + other.generator + ")";
        return this.flipped()
            .intersection(other.flipped(), true, true, generator, infoFactory)
            .flipped();
    }
    and(other, infoFactory) {
        const generator = this.generator &&
            other.generator &&
            this.generator + ".and(" + other.generator + ")";
        return this.intersection(other, true, true, generator, infoFactory);
    }
    xor(other, infoFactory) {
        const generator = this.generator &&
            other.generator &&
            this.generator + ".xor(" + other.generator + ")";
        return new BRep(this.minus(other, infoFactory).faces.concat(other.minus(this, infoFactory).faces), this.infiniteVolume != other.infiniteVolume, generator);
    }
    equals(obj) {
        return (this.faces.length == obj.faces.length &&
            this.faces.every((face) => obj.faces.some((face2) => face.equals(face2))));
    }
    like(brep) {
        return (this.faces.length == brep.faces.length &&
            this.faces.every((face) => brep.faces.some((face2) => face.likeFace(face2))));
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
        return `new BRep([\n${this.faces.join(",\n").replace(/^/gm, "\t")}], ${this.infiniteVolume})`;
    }
    getConstructorParameters() {
        return [this.faces, this.infiniteVolume];
    }
    toSource(useGenerator = true) {
        return ((useGenerator && this.generator) ||
            `new BRep([\n${this.faces.map(SCE).join(",\n").replace(/^/gm, "\t")}], ${this.infiniteVolume})`);
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
            const usableOldEdges = face
                .getAllEdges()
                .filter((edge) => !edgeSubEdges.get(edge));
            const subEdges = face
                .getAllEdges()
                .mapFilter((edge) => edgeSubEdges.get(edge))
                .concatenated();
            const newEdges = faceMap.get(face) || [];
            if (newEdges.length || subEdges.length) {
                oldFaceStatuses.set(face, "partial");
                const loops = [];
                // new edges are definitely part of a resulting loop
                // old edges (both contour and holes) can either be part of a new loop, in which case they will already
                // have been visited when starting a loop search with a new edge, OR they can be stranded, OR they can
                // remain in their old loop
                function getNextStart() {
                    return (newEdges.find((edge) => !visitedEdges.has(edge)) ||
                        subEdges.find((edge) => !visitedEdges.has(edge)) ||
                        usableOldEdges.find((edge) => !visitedEdges.has(edge)));
                }
                const visitedEdges = new Set();
                // search for a loop:
                let currentEdge;
                while ((currentEdge = getNextStart())) {
                    const startEdge = currentEdge, edges = [];
                    let i = 0;
                    // wether only new edges are used (can include looseSegments)
                    do {
                        visitedEdges.add(currentEdge);
                        edges.push(currentEdge);
                        // find next edge
                        const possibleOldEdges = usableOldEdges.filter((edge) => currentEdge.b.like(edge.a));
                        const possibleSubEdges = subEdges.filter((edge) => currentEdge.b.like(edge.a));
                        const possibleNewEdges = newEdges.filter((edge) => currentEdge.b.like(edge.a));
                        const possibleEdges = possibleOldEdges.concat(possibleSubEdges, possibleNewEdges);
                        if (0 == possibleEdges.length)
                            break;
                        assert(0 < possibleEdges.length, () => face.sce);
                        const faceNormalAtCurrentB = face.surface.normalP(currentEdge.b);
                        const nextEdgeIndex = calcNextEdgeIndex(currentEdge, possibleEdges, faceNormalAtCurrentB);
                        currentEdge = possibleEdges[nextEdgeIndex];
                        if (visitedEdges.has(currentEdge)) {
                            break;
                        }
                        assert(currentEdge);
                        assert(currentEdge != startEdge);
                    } while (++i < 400);
                    if (400 == i) {
                        assert(false, "too many");
                    }
                    // check if we found a loop
                    if (edges.length > 1 && currentEdge == startEdge) {
                        loops.push(edges);
                    }
                }
                const faceNewFaces = BRep.assembleFacesFromLoops(loops, face.surface, face, infoFactory);
                newFaces.push(...faceNewFaces);
                const faceNewFacesEdges = faceNewFaces.flatMap((face) => face.getAllEdges());
                insideEdges.push(...usableOldEdges.filter((edge) => faceNewFacesEdges.includes(edge)));
            }
        }
        while (insideEdges.length != 0) {
            const insideEdge = insideEdges.pop();
            const adjacentFaces = this.edgeFaces.get(insideEdge.getCanon());
            adjacentFaces.forEach((info) => {
                if (!oldFaceStatuses.has(info.face)) {
                    oldFaceStatuses.set(info.face, "inside");
                    insideEdges.push.apply(insideEdges, info.face.getAllEdges());
                }
            });
        }
        newFaces.push(...oldFaces.filter((face) => oldFaceStatuses.get(face) == "inside"));
    }
    static getLooseEdgeSegments(edgePointInfoss, edgeFaces) {
        const result = new JavaMap();
        // if there are no point info, the original edge will be kept, so we should return nothing
        // otherwise, something will be returned, even if it a new edge identical to the base edge
        for (const [canonEdge, pointInfos] of edgePointInfoss) {
            if (0 == pointInfos.length)
                continue;
            const allFaces = edgeFaces.get(canonEdge);
            pointInfos.sort((a, b) => snap0(a.edgeT - b.edgeT) || +!!undefined);
            let startP = canonEdge.a, startDir = canonEdge.aDir, startT = canonEdge.aT, startInfo;
            function addNewEdge(startInfo, endInfo, newEdge) {
                for (let i = 0; i < allFaces.length; i++) {
                    const faceInfo = allFaces[i];
                    mapPush(result, !faceInfo.reversed ? canonEdge : canonEdge.flipped(), !faceInfo.reversed ? newEdge : newEdge.flipped());
                }
            }
            for (let i = 0; i < pointInfos.length; i++) {
                const info = pointInfos[i];
                const pDir = canonEdge.tangentAt(info.edgeT);
                if (!eq(info.edgeT, startT)) {
                    const newEdge = createEdge(canonEdge.curve, startP, info.p, startT, info.edgeT, undefined, startDir, pDir, "looseSegment" + getGlobalId());
                    addNewEdge(startInfo, info, newEdge);
                }
                startP = info.p;
                startT = info.edgeT;
                startInfo = info;
                startDir = pDir;
            }
            if (startInfo && !eq(startT, canonEdge.bT)) {
                const newEdge = createEdge(canonEdge.curve, startP, canonEdge.b, startT, canonEdge.bT, undefined, startDir, canonEdge.bDir, "looseSegment" + getGlobalId());
                addNewEdge(startInfo, undefined, newEdge);
            }
        }
        return result;
    }
    getIntersectionEdges(brep2) {
        const faceMap = new Map(), thisEdgePoints = new JavaMap(), otherEdgePoints = new JavaMap();
        const checkedPairs = new JavaSet();
        this.faces.forEach((face) => {
            //console.log('face', face.toString())
            brep2.faces.forEach((face2) => {
                //console.log('face2', face2.toString())
                face.intersectFace(face2, this, brep2, faceMap, thisEdgePoints, otherEdgePoints, checkedPairs);
            });
        });
        return Array.from(faceMap.values()).concatenated();
    }
    shellCount() {
        const foundFaces = new Set();
        let face, result = 0;
        while ((face = this.faces.find((face) => !foundFaces.has(face)))) {
            result++;
            const stack = [face];
            while ((face = stack.pop())) {
                // @ts-ignore
                for (const edge of face.getAllEdges()) {
                    // @ts-ignore
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
        return AABB.forAABBs(this.faces.map((face) => face.getAABB()));
    }
    assertSanity() {
        if (!NLA_DEBUG)
            return;
        // const allFaceEdges = this.faces.flatMap(face => face.getAllEdges())
        // for (const { i, j } of combinations(allFaceEdges.length)) {
        // const a = allFaceEdges[i],
        // 	b = allFaceEdges[j]
        // assert(i == j || !a.isCoEdge(b) || a == b || a.flippedOf == b, 'coedges not linked properly', a, b)
        // assert(
        // 	i == j ||
        // 		!a.curve.isColinearTo(b.curve) ||
        // 		(a.curve.equals(b.curve) && a.isCoEdge(b)) ||
        // 		!a.overlaps(b),
        // 	'colinear edges overlap',
        // 	a,
        // 	b,
        // )
        // }
        this.buildAdjacencies();
        for (const [canonEdge, edgeFaceInfos] of this.edgeFaces) {
            // TODO handle curved faces
            assert(edgeFaceInfos.length % 2 == 0, () => canonEdge + edgeFaceInfos.sce);
        }
    }
    //intersection3(other: BRep, buildThis: boolean, buildOther: boolean, name?: string): BRep {
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
    //                    const commonEdge = createEdge(curve1, min(edge1.minT, minT), min(edge1.maxT, maxT), )
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
        this.edgeFaces = new JavaMap();
        for (const face of this.faces) {
            for (const edge of face.getAllEdges()) {
                const canon = edge.getCanon();
                const normalAtCanonA = face.surface.normalP(canon.a);
                const inside = normalAtCanonA.cross(canon == edge ? edge.aDir : edge.bDir);
                mapPush(this.edgeFaces, canon, {
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
            const faceInfo0 = edgeFaceInfos.find((faceInfo) => faceInfo.reversed);
            if (!faceInfo0) {
                console.warn("invalid brep");
                continue;
            }
            edgeFaceInfos.forEach((faceInfo) => {
                if (faceInfo != faceInfo0) {
                    faceInfo.angle = faceInfo0.inside.angleRelativeNormal(faceInfo.inside, canonEdge.aDir.unit());
                    if (faceInfo.angle < 0)
                        faceInfo.angle += 2 * Math.PI;
                }
            });
            edgeFaceInfos.sort((a, b) => snap(a.angle - b.angle, 0)); // TODO  || assertNever()
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
        const thisEdgePoints = new JavaMap(), otherEdgePoints = new JavaMap();
        const checkedPairs = new JavaSet();
        for (const thisFace of this.faces) {
            for (const otherFace of other.faces) {
                thisFace.intersectFace(otherFace, this, other, faceMap, thisEdgePoints, otherEdgePoints, checkedPairs);
            }
        }
        for (const edge of thisEdgePoints.keys()) {
            assert(this.edgeFaces.get(edge));
        }
        for (const edge of otherEdgePoints.keys()) {
            assert(other.edgeFaces.get(edge));
        }
        const newFaces = [];
        if (0 == faceMap.size &&
            0 == thisEdgePoints.size &&
            0 == otherEdgePoints.size) {
            const thisInOther = other.containsPoint(this.faces[0].contour[0].a, true) !==
                other.infiniteVolume;
            const otherInThis = !thisInOther &&
                this.containsPoint(other.faces[0].contour[0].a) !== this.infiniteVolume;
            if (thisInOther || otherInThis) {
                const [inside, outside] = thisInOther ? [this, other] : [other, this];
                if (inside.infiniteVolume) {
                    if (outside.infiniteVolume) {
                        return outside;
                    }
                    else {
                        return BRep.join([inside, outside]);
                    }
                }
                else {
                    if (outside.infiniteVolume) {
                        return BRep.EMPTY;
                    }
                    else {
                        return inside;
                    }
                }
            }
            else {
                if (this.infiniteVolume) {
                    if (other.infiniteVolume) {
                        return BRep.join([this, other]);
                    }
                }
                else {
                    if (other.infiniteVolume) {
                        return this;
                    }
                    else {
                        return BRep.EMPTY;
                    }
                }
            }
            return BRep.EMPTY;
        }
        else {
            if (buildThis) {
                const edgeLooseSegments = BRep.getLooseEdgeSegments(thisEdgePoints, this.edgeFaces);
                // @ts-ignore
                const els = this.faces.map((face) => [
                    face,
                    Array.from(edgeLooseSegments.entries()).flatMap(([edge, subs]) => face.getAllEdges().some((e) => e.equals(edge)) ? subs : []),
                ]);
                this.reconstituteFaces(this.faces, edgeLooseSegments, faceMap, newFaces, infoFactory);
            }
            if (buildOther) {
                const edgeLooseSegments = BRep.getLooseEdgeSegments(otherEdgePoints, other.edgeFaces);
                // @ts-ignore
                const els = other.faces.map((face) => [
                    face,
                    Array.from(edgeLooseSegments.entries()).flatMap(([edge, subs]) => face.getAllEdges().some((e) => e.equals(edge)) ? subs : []),
                ]);
                other.reconstituteFaces(other.faces, edgeLooseSegments, faceMap, newFaces, infoFactory);
            }
        }
        //buildCoplanar && this.reconstituteCoplanarFaces(likeSurfaceFaces, edgeLooseSegments, faceMap, newFaces,
        // this.infiniteVolume, other.infiniteVolume)
        const result = new BRep(newFaces, this.infiniteVolume && other.infiniteVolume, generator);
        //result.buildAdjacencies()
        return result;
    }
    transform(m4, desc) {
        let vertexNames;
        if (this.vertexNames) {
            vertexNames = new Map();
            this.vertexNames.forEach((name, vertex) => vertexNames.set(m4.transformPoint(vertex), name + desc));
        }
        return new BRep(this.faces.map((f) => f.transform(m4)), this.infiniteVolume, this.generator && desc && this.generator + desc, // if desc isn't set, the generator will be invalid
        vertexNames);
    }
    transform4(m4, desc) {
        let vertexNames;
        if (this.vertexNames) {
            vertexNames = new Map();
            this.vertexNames.forEach((name, vertex) => vertexNames.set(m4.transformPoint(vertex), name + desc));
        }
        return new BRep(this.faces.map((f) => f.transform4(m4)), this.infiniteVolume, this.generator && desc && this.generator + desc, // if desc isn't set, the generator will be invalid
        vertexNames);
    }
    flipped() {
        return new BRep(this.faces.map((f) => f.flipped()), !this.infiniteVolume, this.generator && this.generator + ".flipped()", this.vertexNames);
    }
}
BRep.EMPTY = new BRep([], false, "BRep.EMPTY", new Map()).buildAdjacencies();
BRep.R3 = new BRep([], true, "BRep.R3", new Map()).buildAdjacencies();
function dotCurve(v, cDir, cDDT) {
    let dot = v.dot(cDir);
    if (eq0(dot)) {
        dot = v.dot(cDDT);
    }
    assert(!eq0(dot));
    return dot;
}
function dotCurve2(curve, t, normal, sign) {
    assert(sign == 1 || sign == -1, sign);
    const tangentDot = curve.tangentAt(t).dot(normal);
    // if tangentDot != 0 the curve simply crosses the plane
    if (!eq0(tangentDot)) {
        return sign * tangentDot;
    }
    if (curve.ddt) {
        const ddtDot = curve.ddt(t).dot(normal);
        // tangentDot == 0 ==> critical point at t, if ddtDot != 0, then it is a turning point, otherwise we can't be sure
        // and must do a numeric test
        if (!eq0(ddtDot)) {
            return ddtDot;
        }
    }
    const numericDot = curve
        .at(t)
        .to(curve.at(t + sign * 4 * NLA_PRECISION))
        .dot(normal);
    assert(!(curve instanceof L3));
    return numericDot;
}
const INSIDE = 0, OUTSIDE = 1, COPLANAR_SAME = 2, COPLANAR_OPPOSITE = 3, ALONG_EDGE_OR_PLANE = 4;
/**
 *
 * @param brep BREP to check
 * @param edge edge to check
 * @param dirAtEdgeA the direction vector to check
 * @param faceNormal If dirAtEdgeA doesn't split a volume, but is along a face, the returned value depends on
 *     wether that faces normal1 points in the same direction as faceNormal
 * @returns INSIDE, OUTSIDE, COPLANAR_SAME or COPLANAR_OPPOSITE
 */
//function splitsVolumeEnclosingFaces(brep: BRep, edge: Edge, dirAtEdgeA: V3, faceNormal: V3): int {
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
    assert(arguments.length == 4);
    assert(canonEdge == canonEdge.getCanon());
    //assert(p.equals(canonEdge.a))
    const edgeFaceInfos = brep.edgeFaces.get(canonEdge);
    assertf(() => edgeFaceInfos.length % 2 == 0);
    assertf(() => brep.edgeFaces);
    const faceInfo0 = edgeFaceInfos[0];
    const aDir1 = canonEdge.aDir.unit();
    const angleToCanon = ((faceInfo0.inside.angleRelativeNormal(dirAtEdgeA, aDir1) +
        2 * Math.PI +
        NLA_PRECISION) %
        (2 * Math.PI)) -
        NLA_PRECISION;
    const nearestFaceInfoIndex = edgeFaceInfos.findIndex((faceInfo) => lt(angleToCanon, faceInfo.angle));
    const nearestFaceInfo = edgeFaceInfos[nearestFaceInfoIndex == -1
        ? edgeFaceInfos.length - 1
        : nearestFaceInfoIndex - 1];
    if (eq(nearestFaceInfo.angle, angleToCanon)) {
        //assert(false) todo
        const coplanarSame = nearestFaceInfo.normalAtCanonA.dot(faceNormal) > 0;
        return coplanarSame ? COPLANAR_SAME : COPLANAR_OPPOSITE;
    }
    else {
        return nearestFaceInfo.reversed ? INSIDE : OUTSIDE;
    }
}
function splitsVolumeEnclosingFacesP(brep, canonEdge, p, pInside, pFaceNormal) {
    assert(arguments.length == 5);
    assert(canonEdge == canonEdge.getCanon());
    //assert(p.equals(canonEdge.a))
    assertf(() => brep.edgeFaces);
    const edgeFaceInfos = brep.edgeFaces.get(canonEdge);
    assertf(() => edgeFaceInfos.length % 2 == 0);
    const pDir1 = canonEdge.tangentAt(canonEdge.curve.pointT(p)).unit();
    const faceInfoAngleFromPInsideNeg = (faceInfo) => {
        const faceInfoPDir = faceInfo.edge.getCanon() == faceInfo.edge ? pDir1 : pDir1.negated();
        const faceInfoInsideAtP = faceInfo.face.surface
            .normalP(p)
            .cross(faceInfoPDir);
        const faceInfoAngleAtP = pInside.angleRelativeNormal(faceInfoInsideAtP, pDir1);
        return -(((faceInfoAngleAtP + TAU + NLA_PRECISION) % TAU) - NLA_PRECISION);
    };
    const nearestFaceInfo = withMax(edgeFaceInfos, faceInfoAngleFromPInsideNeg);
    if (eq0(faceInfoAngleFromPInsideNeg(nearestFaceInfo))) {
        //assert(false) todo
        const coplanarSame = nearestFaceInfo.face.surface.normalP(p).dot(pFaceNormal) > 0;
        return coplanarSame ? COPLANAR_SAME : COPLANAR_OPPOSITE;
    }
    else {
        return nearestFaceInfo.reversed ? OUTSIDE : INSIDE;
    }
}
function splitsVolumeEnclosingFacesP2(brep, canonEdge, p, testCurve, curveT, dir, faceNormal) {
    assert(canonEdge == canonEdge.getCanon());
    //assert(p.equals(canonEdge.a))
    assertf(() => brep.edgeFaces);
    const edgeFaceInfos = brep.edgeFaces.get(canonEdge);
    assertf(() => edgeFaceInfos.length % 2 == 0);
    const pDir1 = canonEdge.tangentAt(canonEdge.curve.pointT(p)).unit();
    let pInside = testCurve.tangentAt(curveT).times(dir);
    if (pInside.isParallelTo(pDir1)) {
        pInside = testCurve
            .diff(curveT, (1e-4 * dir) / testCurve.tangentAt(curveT).length())
            .rejectedFrom(pDir1);
        pInside = pInside.div(pInside.length());
    }
    let minValue = 20, advanced = false, result = OUTSIDE;
    for (const faceInfo of edgeFaceInfos) {
        const faceInfoPDir = faceInfo.edge.getCanon() == faceInfo.edge ? pDir1 : pDir1.negated();
        const faceInfoInsideAtP = faceInfo.face.surface
            .normalP(p)
            .cross(faceInfoPDir);
        const faceInfoAngleAtP = pInside.angleRelativeNormal(faceInfoInsideAtP, pDir1);
        const angle = ((faceInfoAngleAtP + TAU + NLA_PRECISION) % TAU) - NLA_PRECISION;
        if (eq0(angle)) {
            // do advanced analysis
            const normVector = faceInfo.face.surface.normalP(p);
            if (faceInfo.face.surface.containsCurve(testCurve)) {
                const coplanarSame = normVector.dot(faceNormal) > 0;
                return coplanarSame ? COPLANAR_SAME : COPLANAR_OPPOSITE;
            }
            const testPlane = P3.normalOnAnchor(pDir1, p);
            const isCurve = faceInfo.face.surface.isCurvesWithPlane(testPlane)[0];
            const isCurvePT = isCurve.pointT(p);
            const dirFactor = sign(isCurve.tangentAt(isCurvePT).dot(pInside));
            const eps = 1e-4;
            const iscd = isCurve
                .at(isCurvePT)
                .to(isCurve.at(isCurvePT + dir * dirFactor * eps))
                .dot(normVector);
            const ecd = testCurve
                .at(curveT)
                .to(testCurve.at(curveT + dir * eps))
                .dot(normVector);
            const diff = (iscd - ecd) * (faceInfo.reversed ? -1 : 1);
            if (diff > 0 && (!advanced || diff < minValue)) {
                advanced = true;
                minValue = diff;
                result = faceInfo.reversed ? OUTSIDE : INSIDE;
            }
        }
        else if (!advanced) {
            if (angle < minValue) {
                minValue = angle;
                result = faceInfo.reversed ? OUTSIDE : INSIDE;
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
        assertf(() => planeFace instanceof PlaneFace);
        if (planeFace.getAllEdges().some((edge) => edge.a.like(p))) {
            if (testPlane.isParallelToPlane(planeFace.surface.plane)) {
                if (planeFace.pointsToInside(p, dir) != PointVsFace.OUTSIDE) {
                    return ALONG_EDGE_OR_PLANE;
                }
            }
            else {
                const isLine = L3.fromPlanes(testPlane, planeFace.surface.plane);
                const ps = planeFace.edgeISPsWithPlane(isLine, testPlane);
                let i = 0;
                while (i < ps.length) {
                    const a = ps[i++], b = ps[i++];
                    const out = a.p.like(p);
                    if (out || b.p.like(p)) {
                        const dir2 = out ? isLine.dir1 : isLine.dir1.negated();
                        const angle = (dir.angleRelativeNormal(dir2, testPlane.normal1) +
                            2 * Math.PI +
                            NLA_PRECISION / 2) %
                            (2 * Math.PI);
                        rays.push({ angle: angle, out: out });
                    }
                }
            }
        }
    }
    rays.sort((a, b) => a.angle - b.angle);
    //console.log("testPlane", testPlane.toSource(), "rays", rays.toSource())
    if (eq0(rays[0].angle)) {
        return ALONG_EDGE_OR_PLANE;
    }
    else {
        return rays[0].out ? OUTSIDE : INSIDE;
    }
}
function splitsVolumeEnclosingCone2(brep, p, curve, curveT, fb) {
    assert(curve.containsPoint(p));
    const pFaces = brep.faces.filter((face) => face.getAllEdges().some((edge) => edge.a.like(p)));
    for (let k = 0; k < pFaces.length; k++) {
        const face = pFaces[k];
        if (face.surface.containsCurve(curve)) {
            //assert(false)
            if (face.pointsToInside3(p, curve, curveT, fb) != PointVsFace.OUTSIDE) {
                return ALONG_EDGE_OR_PLANE;
            }
        }
    }
    const EPS = 1e-6;
    return brep.containsPoint(curve.at(curveT + fb * EPS), true)
        ? INSIDE
        : OUTSIDE;
}
function fff(info, surface) {
    const canonA = info.edge.reversed ? info.edge.b : info.edge.a;
    const surfaceNormalAtCanonA = surface.normalP(canonA);
    const dot = snap0(info.inside.dot(surfaceNormalAtCanonA));
    if (0 !== dot) {
        return 0 < dot ? OUTSIDE : INSIDE;
    }
    if (surface.isCoplanarTo(info.face.surface)) {
        return 0 < info.normalAtCanonA.dot(surfaceNormalAtCanonA)
            ? COPLANAR_SAME
            : COPLANAR_OPPOSITE;
    }
    throw new Error();
}
function triangulateVertices(normal, vertices, holeStarts) {
    const absMaxDim = normal.maxAbsDim(), factor = sign(normal.e(absMaxDim));
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
    assertNumbers(a, b, c);
    // TODO: disambiguate on a < b
    const term = sqrt(a * a + b * b - c * c);
    return {
        x1: (a * c + b * term) / (a * a + b * b),
        x2: (a * c - b * term) / (a * a + b * b),
        y1: (b * c - a * term) / (a * a + b * b),
        y2: (b * c + a * term) / (a * a + b * b),
    };
}
function intersectionUnitCircleLine2(a, b, c) {
    assertNumbers(a, b, c);
    // TODO: disambiguate on a < b
    // cf. pqFormula
    const termSqr = snap0(a * a + b * b - c * c);
    if (termSqr < 0) {
        return [];
    }
    else if (termSqr == 0) {
        return [[(a * c) / (a * a + b * b), (b * c) / (a * a + b * b)]];
    }
    else {
        const term = sqrt(termSqr);
        return [
            [
                (a * c + b * term) / (a * a + b * b),
                (b * c - a * term) / (a * a + b * b),
            ],
            [
                (a * c - b * term) / (a * a + b * b),
                (b * c + a * term) / (a * a + b * b),
            ],
        ];
    }
}
function intersectionCircleLine(a, b, c, r) {
    assertNumbers(a, b, c, r);
    const term = sqrt(r * r * (a * a + b * b) - c * c);
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
    assertNumbers(a, b, c);
    const aa = a * a, bb = b * b, cc = c * c;
    // TODO: disambiguate on a < b
    //var xTerm = sqrt(4*cc*aa-4*(bb-aa)*(-cc-bb))
    const xTerm = 2 * sqrt(bb * cc + bb * bb - aa * bb);
    const yTerm = sqrt(4 * cc * bb - 4 * (bb - aa) * (cc - aa));
    return {
        x1: (-2 * a * c + xTerm) / 2 / (bb - aa),
        x2: (-2 * a * c - xTerm) / 2 / (bb - aa),
        y1: (2 * b * c - yTerm) / 2 / (bb - aa),
        y2: (2 * b * c + yTerm) / 2 / (bb - aa),
    };
}
function curvePointPP(ps1, ps2, startPoint, dontCheck) {
    const EPS = NLA_PRECISION / 4;
    //if (!dontCheck) {
    //    const p = curvePointPP(ps1, ps2, startPoint, true).p
    //    if (!ps1.containsPoint(p)) {
    //        console.log("foo, startPoint was " + startPoint.sce)
    //        ps1.containsPoint(p)
    //    }
    //}
    let Q = startPoint;
    let st1 = ps1.pointFoot(Q);
    let st2 = ps2.pointFoot(Q);
    let a, b, aNormal, bNormal, abNormalsCross;
    //console.log("curvePointPP, startPoint was " + startPoint.sce)
    //console.log(Q.sce+ ',')
    let i = 16;
    do {
        a = ps1.pUV(st1.x, st1.y);
        b = ps2.pUV(st2.x, st2.y);
        if (eq0(a.distanceTo(b), EPS))
            break;
        // drPs.push({p:a,text:'a'+j+' '+i})
        // drPs.push({p:b,text:'b'+j+' '+i})
        aNormal = ps1.normalUV(st1.x, st1.y);
        bNormal = ps2.normalUV(st2.x, st2.y);
        // next Q is the intersection of the planes
        // (x - a) * aNormal,
        // (x - b) * bNormal and
        // (x - Q) * (aNormal X bNormal)
        abNormalsCross = aNormal.cross(bNormal);
        // drVs.push({anchor: Q, dir: aNormal})
        // drVs.push({anchor: Q, dir: bNormal})
        Q = V3.add(bNormal.cross(abNormalsCross).times(a.dot(aNormal)), abNormalsCross.cross(aNormal).times(b.dot(bNormal)), abNormalsCross.times(abNormalsCross.dot(Q))).div(abNormalsCross.squared());
        //console.log(Q.sce+ ',')
        // feet of Q on ps1 and ps2 (closest points)
        st1 = ps1.pointFoot(Q, st1.x, st1.y);
        st2 = ps2.pointFoot(Q, st2.x, st2.y);
    } while (--i);
    //assert(ps1.containsPoint(Q), Q, ps1)
    //assert(ps2.containsPoint(Q))
    if (!eq0(a.distanceTo(b), EPS)) {
        return undefined;
    }
    return { p: Q, st1: st1, st2: st2 };
}
/**
 * Follow the intersection curve of two parametric surfaces starting from a given point.
 * @param {ParametricSurface} ps1
 * @param {ParametricSurface} ps2
 * @param {number} s1Step
 * @param {number} t1Step
 * @param {number} s2Step
 * @param {number} t2Step
 * @param {number} curveStepSize
 * @return {Curve[]}
 */
function followAlgorithmPP(ps1, ps2, startPoint, curveStepSize, bounds1 = uvInAABB2.bind(undefined, ps1), bounds2 = uvInAABB2.bind(undefined, ps2)) {
    const points = [];
    const tangents = [];
    const st1s = [];
    const st2s = [];
    let Q = startPoint;
    let st1 = ps1.uvP(Q);
    let st2 = ps2.uvP(Q);
    assert(ps1.pUV(st1.x, st1.y).like(Q));
    assert(st1.like(ps1.pointFoot(Q, st1.x, st1.y)));
    assert(st2.like(ps2.pointFoot(Q, st2.x, st2.y)));
    assert(ps2.pUV(st2.x, st2.y).like(Q));
    for (let i = 0; i < 1000; i++) {
        ({ p: Q, st1, st2 } = curvePointPP(ps1, ps2, Q));
        assert(ps1.containsPoint(Q), Q, ps1);
        assert(ps2.containsPoint(Q));
        const aNormal = ps1.normalUV(st1.x, st1.y);
        const bNormal = ps2.normalUV(st2.x, st2.y);
        const tangent = aNormal.cross(bNormal).toLength(curveStepSize);
        tangents.push(tangent);
        points.push(Q);
        st1s.push(st1);
        st2s.push(st2);
        if (i > 4) {
            if (!bounds1(st1.x, st1.y) || !bounds2(st2.x, st2.y)) {
                break;
            }
        }
        Q = Q.plus(tangent);
    }
    return { points, tangents, st1s, st2s };
}
/**
 * Iteratively calculate points on an implicit 2D curve.
 * @param ic The curve in question.
 * @param startP The point at which to start.
 * @param stepLength The step the algorithm takes. Will be the approximate distance between points.
 * @param bounds Bounds function.
 * @param endP End point. If undefined, algorithm will continue until out of bounds or back at start point.
 * @param startTangent TODO Ignore this.
 * @returns Calculated points and tangents. points[0] and tangents[0] will be startP and startTangent.
 */
function followAlgorithm2d(ic, startP, stepLength = 0.5, bounds, validUV, endP, startTangent) {
    assertNumbers(stepLength, ic(0, 0));
    assertVectors(startP);
    if (!startTangent) {
        startTangent = new V3(-ic.y(startP.x, startP.y), ic.x(startP.x, startP.y), 0).toLength(stepLength);
    }
    assertVectors(startTangent);
    const points = [];
    const tangents = [];
    assert(eq0(ic(startP.x, startP.y), 0.01), "isZero(implicitCurve(startPoint.x, startPoint.y))", ic(startP.x, startP.y));
    let i = 0, p = startP, tangent = startTangent, fullLoop = false;
    do {
        points.push(p);
        tangents.push(tangent);
        const searchStart = p.plus(tangent);
        assert(searchStart);
        const newP = curvePointMF(ic, searchStart);
        const dfpdx = ic.x(newP.x, newP.y), dfpdy = ic.y(newP.x, newP.y);
        const newTangent = new V3(-dfpdy, dfpdx, 0).toLength(stepLength);
        //const reversedDir = p.minus(prevp).dot(tangent) < 0
        assert(!p.equals(newP));
        // check if we passed a singularity
        if (tangent.dot(newTangent) < 0) {
            const singularity = newtonIterate2d(ic.x, ic.y, p.x, p.y);
            if (eq0(ic(singularity.x, singularity.y)) &&
                singularity.distanceTo(p) < abs(stepLength)) {
                // end on this point
                points.push(singularity);
                tangents.push(p.to(singularity));
                break;
            }
            else {
                throw new Error();
            }
        }
        // check for endP
        if (endP && p.equals(endP)) {
            break;
        }
        // check if loop
        if (fullLoop) {
            if (p.distanceTo(startP) > abs(stepLength)) {
                points.pop();
                tangents.pop();
                assert(points.last.distanceTo(startP) <= abs(stepLength));
                break;
            }
        }
        else {
            if (i > 4 && p.distanceTo(startP) <= abs(stepLength)) {
                fullLoop = true;
            }
        }
        // check if out of bounds
        if (i > 1 && !uvInAABB2(bounds, p.x, p.y)) {
            const endP = figureOutBorderPoint(bounds, p, ic);
            points.pop();
            tangents.pop();
            if (points.last.distanceTo(endP) < abs(stepLength) / 2) {
                points.pop();
                tangents.pop();
            }
            const endTangent = new V3(-ic.y(endP.x, endP.y), ic.x(endP.x, endP.y), 0).toLength(stepLength);
            points.push(endP);
            tangents.push(endTangent);
            break;
        }
        if (i > 4 && !validUV(p.x, p.y)) {
            break;
        }
        assert(eq0(ic(newP.x, newP.y), NLA_PRECISION * 2), p, newP, searchStart, ic(newP.x, newP.y));
        tangent = newTangent;
        p = newP;
    } while (++i < 1000);
    assert(i < 1000);
    //assert(points.length > 6)
    return { points, tangents };
}
/**
 * Given a point p just outside the bounds, figure out the nearby intersection of the bounds with the ic.
 * @param bounds
 * @param p
 * @param ic
 */
function figureOutBorderPoint(bounds, p, ic) {
    if (p.x < bounds.uMin || bounds.uMax < p.x) {
        const u = bounds.uMax < p.x ? bounds.uMax : bounds.uMin;
        const v = newtonIterateWithDerivative((t) => ic(u, t), p.y, 4, (t) => ic.y(u, t));
        if (uvInAABB2(bounds, u, v)) {
            return new V3(u, v, 0);
        }
    }
    if (p.y < bounds.vMin || bounds.vMax < p.y) {
        const v = bounds.vMax < p.y ? bounds.vMax : bounds.vMin;
        const u = newtonIterateWithDerivative((s) => ic(s, v), p.x, 4, (s) => ic.x(s, v));
        assert(uvInAABB2(bounds, u, v));
        return new V3(u, v, 0);
    }
    throw new Error(p + " " + bounds);
}
function followAlgorithm2dAdjustable(ic, start, stepLength = 0.5, bounds, endp = start) {
    assertNumbers(stepLength, ic(0, 0));
    assertVectors(start);
    //assert (!startDir || startDir instanceof V3)
    const points = [];
    const tangents = [];
    assert(eq0(ic(start.x, start.y), 0.01), "isZero(implicitCurve(startPoint.x, startPoint.y))");
    let p = start, prevp = p;
    let i = 0;
    do {
        const dfpdx = ic.x(p.x, p.y), dfpdy = ic.y(p.x, p.y);
        const dfpdxx = ic.xx(p.x, p.y), dfpdyy = ic.yy(p.x, p.y), dfpdxy = ic.xy(p.x, p.y);
        const c2factor = abs((Math.pow(dfpdy, 2) * dfpdxx - 2 * dfpdx * dfpdy * dfpdxy + Math.pow(dfpdx, 2) * dfpdyy) /
            Math.pow((Math.pow(dfpdx, 2) + Math.pow(dfpdy, 2)), 2));
        const c2 = new V3(dfpdx, dfpdy, 0).times(c2factor);
        const s = 1 / 16 / c2.length();
        const tangent = new V3(-dfpdy, dfpdx, 0).unit();
        const newPStart = p.plus(tangent.times(s).plus(c2.times(Math.pow(s, 2) / 2)));
        points.push(p);
        tangents.push(tangent);
        prevp = p;
        const newP = curvePointMF(ic, newPStart);
        if (newP.equals(p)) {
            assertNever();
        }
        console.log(p.to(newP).length());
        p = newP;
        assert(eq0(ic(p.x, p.y)));
    } while (i++ < 1000 &&
        (i < 4 || prevp.distanceTo(endp) > stepLength) &&
        bounds(p.x, p.y));
    assert(i != 1000);
    //assert(bounds(p.x, p.y))
    const end = i < 4 || prevp.distanceTo(endp) > stepLength ? p : endp;
    const endTangent = new V3(-ic.y(end.x, end.y), ic.x(end.x, end.y), 0).toLength(stepLength);
    points.push(end);
    tangents.push(endTangent);
    //assert(points.length > 6)
    // TODO gleichmäßige Verteilung der Punkte
    return { points, tangents };
}
// both curves must be in the same s-t coordinates for this to make sense
function intersectionICurveICurve(iCurve1, startParams1, endParams1, startDir, stepLength, iCurve2) {
    assertNumbers(stepLength, iCurve1(0, 0), iCurve2(0, 0));
    assertVectors(startParams1, endParams1);
    assert(!startDir || startDir instanceof V3);
    const vertices = [];
    assert(eq0(iCurve1(startParams1.x, startParams1.y)));
    stepLength = stepLength || 0.5;
    const eps = 1e-5;
    let p = startParams1, prevp = p; // startDir ? p.minus(startDir) : p
    let i = 0;
    while (i++ < 1000 && (i < 4 || p.distanceTo(endParams1) > 1.1 * stepLength)) {
        const fp = iCurve1(p.x, p.y);
        const dfpdx = (iCurve1(p.x + eps, p.y) - fp) / eps, dfpdy = (iCurve1(p.x, p.y + eps) - fp) / eps;
        let tangent = new V3(-dfpdy, dfpdx, 0).toLength(stepLength);
        if (p.minus(prevp).dot(tangent) < 0)
            tangent = tangent.negated();
        prevp = p;
        p = curvePointMF(iCurve1, p.plus(tangent));
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
        val = iCurve2(p.x, p.y);
        if (val * lastVal <= 0) {
            // TODO < ?
            iss.push(newtonIterate2d(iCurve1, iCurve2, p.x, p.y));
        }
    }
    return iss;
}
// export function intersectionPCurveISurface(
// 	parametricCurve: Curve,
// 	searchStart: number,
// 	searchEnd: number,
// 	searchStep: number,
// 	implicitSurface: ImplicitSurface,
// ) {
// 	assertNumbers(searchStart, searchEnd, searchStep)
// 	const iss = []
// 	let val = implicitSurface(parametricCurve(searchStart)),
// 		lastVal
// 	for (let t = searchStart + searchStep; t <= searchEnd; t += searchStep) {
// 		lastVal = val
// 		val = implicitSurface(parametricCurve(t))
// 		if (val * lastVal <= 0) {
// 			iss.push(newtonIterate1d(t => implicitSurface(parametricCurve(t)), t))
// 		}
// 	}
// 	return iss
// }
function cassini(a, c) {
    return (x, y) => (x * x + y * y) * (x * x + y * y) -
        2 * c * c * (x * x - y * y) -
        (Math.pow(a, 4) - Math.pow(c, 4));
}
var MathFunctionR2R;
(function (MathFunctionR2R) {
    function forNerdamer(expression, args = ["x", "y"]) {
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
    function nerdamerToR2_R(expression, args = ["x", "y"]) {
        return expression.buildFunction(args);
    }
    MathFunctionR2R.nerdamerToR2_R = nerdamerToR2_R;
    function forFFxFy(f, fx, fy) {
        f.x = fx;
        f.y = fy;
        return f;
    }
    MathFunctionR2R.forFFxFy = forFFxFy;
})(MathFunctionR2R || (MathFunctionR2R = {}));
const cas2 = cassini(0.9, 1.02);
function arrayLerp(lerp, arr, t) {
    if (0 === t % 1)
        return arr[t];
    return lerp(arr[Math.floor(t)], arr[Math.ceil(t)], t % 1);
}

function doNotSerialize(target, key) {
    const map = target.__SERIALIZATION_BLACKLIST || (target.__SERIALIZATION_BLACKLIST = {});
    map[key] = "no";
}
class ClassSerializer {
    constructor() {
        this.CLASS_NAMES = new Map();
        this.NAME_CLASSES = new Map();
        this.addClass("Object", Object);
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
        Object.keys(namespace).forEach((symbol) => {
            const o = namespace[symbol];
            if ("function" == typeof o && o.name) {
                this.addClass((namespaceName ? namespaceName + "." : "") + symbol, o);
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
        const gatherList = (v) => {
            //console.log(path.toString())
            if (undefined !== v &&
                v.hasOwnProperty("constructor") &&
                this.CLASS_NAMES.has(v.constructor)) ;
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
                        gatherList(v[i]);
                    }
                }
            }
            else if (undefined !== v && "object" == typeof v) {
                if (visited.has(v)) {
                    if (!listMap.has(v)) {
                        listMap.set(v, resultList.length);
                        resultList.push(v);
                    }
                }
                else {
                    assert(!v.__noxTarget || !visited.has(v.__noxTarget));
                    assert(!v.__noxProxy || !visited.has(v.__noxProxy));
                    visited.add(v);
                    if (!v.getConstructorParameters) {
                        for (const key of Object.keys(v).sort()) {
                            if (key == "__noxProxy" || key == "__noxTarget")
                                continue;
                            if (!v.__SERIALIZATION_BLACKLIST ||
                                !v.__SERIALIZATION_BLACKLIST[key]) {
                                gatherList(v[key]);
                            }
                        }
                    }
                    gatherList(Object.getPrototypeOf(v));
                }
            }
        };
        const transform = (v, allowLinks, first) => {
            if ("string" == typeof v ||
                "number" == typeof v ||
                "boolean" == typeof v ||
                null === v) {
                return v;
            }
            if ("undefined" == typeof v) {
                return { "#REF": -1 };
            }
            if (v.hasOwnProperty("constructor") &&
                this.CLASS_NAMES.has(v.constructor)) {
                return { "#REF": this.CLASS_NAMES.get(v.constructor) };
            }
            let index;
            if (allowLinks && !first && undefined !== (index = listMap.get(v))) {
                return { "#REF": index };
            }
            if (Array.isArray(v)) {
                return v.map((x) => transform(x, allowLinks));
            }
            //if (mobx && mobx.isObservableArray(v)) {
            //	const result = {'#PROTO': 'ObservableArray'} as any
            //	v.forEach((val, i) => result[i] = transform(val))
            //	return result
            //}
            if ("object" == typeof v) {
                if (v.getConstructorParameters) {
                    return {
                        "#CONSTRUCTOR": this.CLASS_NAMES.get(v.constructor),
                        "#ARGS": transform(v.getConstructorParameters(), false),
                    };
                }
                const result = {};
                if (Object.prototype !== Object.getPrototypeOf(v)) {
                    result["#PROTO"] = transform(Object.getPrototypeOf(v), allowLinks);
                }
                for (const key of Object.keys(v)) {
                    if (key == "__noxProxy" || key == "__noxTarget")
                        continue;
                    if (!v.__SERIALIZATION_BLACKLIST ||
                        !v.__SERIALIZATION_BLACKLIST[key]) {
                        result[key] = transform(v[key], allowLinks);
                    }
                }
                return result;
            }
            throw new Error("?" + typeof v + v.toString());
        };
        const visited = new Set();
        const listMap = new Map();
        let resultList = [];
        listMap.set(v, 0);
        resultList.push(v);
        gatherList(v);
        resultList = resultList.map((v) => transform(v, true, true));
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
                    fixObject(v[i], (x) => (v[i] = x));
                }
            }
            else if ("object" == typeof v && undefined != v) {
                if ("#CONSTRUCTOR" in v) {
                    const protoName = v["#CONSTRUCTOR"];
                    const proto = this.NAME_CLASSES.get(protoName);
                    assert(proto, protoName + " Missing ");
                    let args = undefined;
                    fixObject(v["#ARGS"], (x) => (args = x));
                    onReady(new proto(...args));
                }
                else if ("#REF" in v) {
                    const ref = v["#REF"];
                    if ("string" == typeof ref) {
                        onReady(this.NAME_CLASSES.get(ref).prototype);
                    }
                    else if ("number" == typeof ref) {
                        if (-1 == ref) {
                            onReady(undefined);
                        }
                        else if (fixedObjects[ref]) {
                            onReady(fixedObjects[ref]);
                        }
                        else {
                            fixObject(tree[ref], (x) => onReady((fixedObjects[ref] = x)));
                        }
                    }
                }
                else {
                    let result;
                    if ("#PROTO" in v) {
                        fixObject(v["#PROTO"], (x) => {
                            result = Object.create(x);
                            onReady(result);
                        });
                    }
                    else {
                        onReady((result = v));
                    }
                    const keys = Object.keys(v);
                    for (let i = 0; i < keys.length; i++) {
                        //if ('name' == keys[i]) console.log(result)
                        if ("#PROTO" != keys[i]) {
                            fixObject(v[keys[i]], (x) => (result[keys[i]] = x));
                            //Object.defineProperty(result, keys[i], {
                            //	value: fixObjects(v[keys[i]]),
                            //	enumerable: true,
                            //	writable: true,
                            //	configurable: true
                            //})
                        }
                    }
                    Object.defineProperty(result, "loadID", {
                        value: getGlobalId(),
                        enumerable: false,
                        writable: false,
                    });
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
        fixObject({ "#REF": 0 }, () => { });
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
	uniform mat4 ts_ModelViewProjectionMatrix;
	uniform mat4 ts_ModelViewMatrix;
	attribute vec4 ts_Vertex;
	uniform mat3 ts_NormalMatrix;
	attribute vec3 ts_Normal;
	uniform vec4 color;
	varying vec3 normal;
	varying vec4 vPosition;
	void main() {
		gl_Position = ts_ModelViewProjectionMatrix * ts_Vertex;
        vPosition = ts_ModelViewMatrix * ts_Vertex;
		normal = normalize(ts_NormalMatrix * ts_Normal);
	}
`;
const vertexShaderWaves = `
	uniform mat4 ts_ModelViewProjectionMatrix;
	uniform mat4 ts_ModelViewMatrix;
	attribute vec4 ts_Vertex;
	uniform mat3 ts_NormalMatrix;
	attribute vec3 ts_Normal;
	uniform vec4 color;
	varying vec3 normal;
	varying vec4 vPosition;
	void main() {
		normal = normalize(ts_NormalMatrix * ts_Normal);
		float offset = mod  (((ts_Vertex.x + ts_Vertex.y + ts_Vertex.z) * 31.0), 20.0) - 10.0;
		vec4 modPos = ts_Vertex + vec4(normal * offset, 0);
		gl_Position = ts_ModelViewProjectionMatrix * modPos;
        vPosition = ts_ModelViewMatrix * modPos;
	}
`;
const vertexShaderBasic = `
	uniform mat4 ts_ModelViewProjectionMatrix;
	attribute vec4 ts_Vertex;
	void main() {
		gl_Position = ts_ModelViewProjectionMatrix * ts_Vertex;
	}
`;
const vertexShaderColor = `
	uniform mat4 ts_ModelViewProjectionMatrix;
	attribute vec4 ts_Vertex;
	attribute vec4 ts_Color;
	varying vec4 fragColor;
	void main() {
		gl_Position = ts_ModelViewProjectionMatrix * ts_Vertex;
		fragColor = ts_Color;
	}
`;
const vertexShaderArc = `
	uniform mat4 ts_ModelViewProjectionMatrix;
	attribute vec4 ts_Vertex;
	uniform float step, offset;
	uniform float radius, width;
	void main() {
		float r = radius;
		float t = offset + ts_Vertex.x * step;
		float pRadius = r - ts_Vertex.y * width;
		vec4 p = vec4(pRadius * cos(t), pRadius * sin(t), 0, 1);
		gl_Position = ts_ModelViewProjectionMatrix * p;
}
`;
const vertexShaderConic3d = `
	uniform mat4 ts_ModelViewProjectionMatrix;
	attribute vec4 ts_Vertex;
	uniform float startT, endT, scale;
	uniform vec3 center, f1, f2;
	uniform int mode;
	float sinh(float x) { return (exp(x) - exp(-x)) / 2.0; }
	float cosh(float x) { return (exp(x) + exp(-x)) / 2.0; }
	void main() {
		float t = startT + ts_Vertex.x * (endT - startT);

		vec3 normal = normalize(cross(f1, f2));

		vec3 p, tangent;
		if (0 == mode) { // ellipse
			p = center + f1 * cos(t) + f2 * sin(t);
			tangent = f1 * -sin(t) + f2 * cos(t);
		}
		if (1 == mode) { // parabola
			p = center + f1 * t + f2 * t * t;
			tangent = f1 + 2.0 * f2 * t;
		}
		if (2 == mode) { // hyperbola
			p = center + f1 * cosh(t) + f2 * sinh(t);
			tangent = f1 * sinh(t) + f2 * cosh(t);
		}
		vec3 outDir = normalize(cross(normal, tangent));
		vec3 p2 = p + scale * (outDir * ts_Vertex.y + normal * ts_Vertex.z);
		gl_Position = ts_ModelViewProjectionMatrix * vec4(p2, 1);
	}
`;
const vertexShaderNURBS = `#version 300 es
	uniform mat4 ts_ModelViewProjectionMatrix;
	in vec4 ts_Vertex;
	uniform float startT, endT, scale;
	uniform vec4 points[32];
	uniform int pointCount, degree;
	uniform float knots[40];
	uniform vec3 normal;
	const int MIN_DEGREE = 1;
	const int MAX_DEGREE = 6;
	
	int tInterval(float t) {
		for (int s = degree; s < 40 - 1 - degree; s++) {
			if (t >= knots[s] && t <= knots[s + 1]) {
				return s;
			}
		}
	}
	
	vec4 stepp(int k, int i, vec4 dkMinus1iMinus1, vec4 dkMinus1i) {
	    return dkMinus1i - dkMinus1iMinus1 * float(k) / (knots[i + degree - k] - knots[i - 1]);
	}
	
	void main() {
		// ts_Vertex.x is in [0, 1]
		float t = startT + ts_Vertex.x * (endT - startT);
		
		int s = tInterval(t);
		
		vec4 v[MAX_DEGREE + 1];
		for (int i = 0; i < degree + 1; i++) {
		    v[i] = points[s - degree + i];
		}
		
		vec4 pTangent4, ddt4 = vec4(0, 0, 1, 0);
		for (int level = 0; level < degree; level++) {
			if (level == degree - 2) {
				// see https://www.globalspec.com/reference/61012/203279/10-8-derivatives
				vec4 a = v[degree];
				vec4 b = v[degree - 1];
				vec4 c = v[degree - 2];
				ddt4 = stepp(degree, s + 1, stepp(degree - 1, s + 1, a, b), stepp(degree - 1, s, b, c));
			}
			if (level == degree - 1) {
				vec4 a = v[degree];
				vec4 b = v[degree - 1];
				pTangent4 = (b - a) * (float(degree) / (knots[s] - knots[s + 1]));
			}
			for (int i = degree; i > level; i--) {
				float alpha = (t - knots[i + s - degree]) / (knots[i + s - level] - knots[i + s - degree]);

				// interpolate each component
                v[i] = (1.0 - alpha) * v[i - 1] + alpha * v[i];
			}
		}
		
		vec4 p4 = v[degree];
		
		vec3 p = p4.xyz / p4.w;
		vec3 pTangent = ((pTangent4.xyz * p4.w) - (p4.xyz * pTangent4.w)) / (p4.w * p4.w);
		vec3 ddt = (
		    p4.xyz * (-p4.w * ddt4.w + 2.0 * pow(pTangent4.w, 2.0))
		    + pTangent4.xyz * (-2.0 * p4.w * pTangent4.w) 
		    + ddt4.xyz * pow(p4.w, 2.0)
        ) / pow(p4.w, 3.0);
		
		vec3 outDir = normalize(cross(ddt, pTangent));
		vec3 correctNormal = normalize(cross(pTangent, outDir));
		vec3 p2 = p + scale * (outDir * ts_Vertex.y + correctNormal * ts_Vertex.z);
		gl_Position = ts_ModelViewProjectionMatrix * vec4(p2, 1);
    }
`;
const vertexShaderBezier = `
    // calculates a bezier curve using ts_Vertex.x as the (t) parameter of the curve
	uniform mat4 ts_ModelViewProjectionMatrix;
	attribute vec4 ts_Vertex;
	uniform float width, startT, endT;
	uniform vec3 p0, p1, p2, p3;
	void main() {
		// ts_Vertex.x is in [0, 1]
		float t = startT + ts_Vertex.x * (endT - startT), s = 1.0 - t;
		float c0 = s * s * s, c1 = 3.0 * s * s * t, c2 = 3.0 * s * t * t, c3 = t * t * t;
		vec3 pPos = p0 * c0 + p1 * c1 + p2 * c2 + p3 * c3;
		float c01 = 3.0 * s * s, c12 = 6.0 * s * t, c23 = 3.0 * t * t;
		vec3 pTangent = (p1 - p0) * c01 + (p2 - p1) * c12 + (p3 - p2) * c23;
		vec3 pNormal = normalize(vec3(pTangent.y, -pTangent.x, 0));
		vec4 p = vec4(pPos - ts_Vertex.y * width * pNormal, 1);
		gl_Position = ts_ModelViewProjectionMatrix * p;
	}
`;
const vertexShaderBezier3d = `
    precision highp float;
    // calculates a bezier curve using ts_Vertex.x as the (t) parameter of the curve
	uniform float scale, startT, endT;
	uniform vec3 ps[4];
	uniform vec3 p0, p1, p2, p3, normal;
	uniform mat4 ts_ModelViewProjectionMatrix;
	attribute vec4 ts_Vertex;
	void main() {
		// ts_Vertex.y is in [0, 1]
		vec3 p5 = ps[0];
		float t = startT * (1.0 - ts_Vertex.x) + endT * ts_Vertex.x, s = 1.0 - t;
		float c0 = s * s * s, 
		      c1 = 3.0 * s * s * t, 
		      c2 = 3.0 * s * t * t, c3 = t * t * t;
		vec3 p = (p0 * c0 + p1 * c1) + (p2 * c2 + p3 * c3);
		float c01 = 3.0 * s * s, 
		      c12 = 6.0 * s * t, 
		      c23 = 3.0 * t * t;
		vec3 pTangent = (p1 - p0) * c01 + (p2 - p1) * c12 + (p3 - p2) * c23;
		vec3 outDir = normalize(cross(normal, pTangent));
		vec3 correctNormal = normalize(cross(pTangent, outDir));
		vec3 p2 = p + scale * (outDir * ts_Vertex.y + correctNormal * ts_Vertex.z);
		gl_Position = ts_ModelViewProjectionMatrix * vec4(p2, 1);
	}
`;
const vertexShaderGeneric = `
	uniform float scale;
	uniform mat4 ts_ModelViewProjectionMatrix;
	attribute vec4 ts_Vertex;
	uniform mat3 ts_NormalMatrix;
	attribute vec3 ts_Normal;
	void main() {
		vec3 normal = normalize(ts_NormalMatrix * ts_Normal);
		vec4 vertexPos = ts_Vertex + vec4(normal * scale, 0);
		gl_Position = ts_ModelViewProjectionMatrix * vertexPos;
	}
`;
const vertexShaderRing = `
	#define M_PI 3.1415926535897932384626433832795
	uniform float step;
	uniform float innerRadius, outerRadius;
	attribute float index;
	uniform mat4 ts_ModelViewProjectionMatrix;
	attribute vec4 ts_Vertex;
	void main() {
		gl_Position = ts_ModelViewProjectionMatrix * vec4(index, index, index, 1);
		float id = atan(ts_Vertex.x, ts_Vertex.y) / M_PI  * 32.0;
		float radius = mod(id, 2.0) < 1.0 ? outerRadius : innerRadius;
		gl_Position = ts_ModelViewProjectionMatrix * vec4(radius * cos(index * step), radius * sin(index * step), 0, 1);
	}
`;
const fragmentShaderColor = `
	precision highp float;
	uniform vec4 color;
	void main() {
		gl_FragColor = color;
	}
`;
const fragmentShaderColor3 = `#version 300 es
	precision highp float;
	uniform vec4 color;
	out vec4 fragColor;
	void main() {
		fragColor = color;
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
	attribute vec4 ts_Vertex;
	uniform mat4 ts_ModelViewProjectionMatrix;
	void main() {
		texturePos = ts_Vertex.xy;
		gl_Position = ts_ModelViewProjectionMatrix * ts_Vertex;
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

function parseGetParams(str) {
    const result = {};
    str.split("&").forEach(function (item) {
        const splitIndex = item.indexOf("=");
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
    RD_FILL: chroma("#9EDBF9"),
    RD_STROKE: chroma("#77B0E0"),
    TS_FILL: chroma("#D19FE3"),
    TS_STROKE: chroma("#A76BC2"),
    PP_FILL: chroma("#F3B6CF"),
    PP_STROKE: chroma("#EB81B4"),
};
class BREPGLContext {
    constructor(gl) {
        this.cachedMeshes = new WeakMap();
        this.shaders = initShaders(gl);
        initMeshes((this.meshes = {}), gl);
    }
    static create(gl) {
        addOwnProperties(gl, BREPGLContext.prototype);
        addOwnProperties(gl, new BREPGLContext(gl));
        return gl;
    }
    drawPoint(p, color = GL_COLOR_BLACK, size = 5) {
        this.pushMatrix();
        this.translate(p);
        this.scale(size / 2, size / 2, size / 2);
        this.shaders.singleColor
            .uniforms({ color: color })
            .draw(this.meshes.sphere1);
        this.popMatrix();
    }
    drawEdge(edge, color = GL_COLOR_BLACK, width = 2) {
        CURVE_PAINTERS[edge.curve.constructor.name](this, edge.curve, color, edge.minT, edge.maxT, width);
    }
    drawCurve(curve, color = GL_COLOR_BLACK, width = 2, tStart, tEnd) {
        CURVE_PAINTERS[curve.constructor.name](this, curve, color, tStart, tEnd, width);
    }
    drawVector(vector, anchor, color = GL_COLOR_BLACK, size = 1) {
        if (vector.likeO())
            return;
        this.pushMatrix();
        const headLength = size * 4;
        if (headLength > vector.length())
            return;
        const vT = vector.getPerpendicular().unit();
        this.multMatrix(M4.forSys(vector.unit(), vT, vector.cross(vT).unit(), anchor));
        this.scale(vector.length() - headLength, size / 2, size / 2);
        this.shaders.singleColor
            .uniforms({
            color: color,
        })
            .draw(this.meshes.vectorShaft);
        this.scale(1 / (vector.length() - headLength), 1, 1);
        this.translate(vector.length() - headLength, 0, 0);
        this.scale(headLength / 2, 1, 1);
        this.shaders.singleColor.draw(this.meshes.vectorHead);
        this.popMatrix();
    }
    drawVectors(drVs, size = undefined) {
        this.drawVector(V3.X, V3.O, chroma("red").gl(), size);
        this.drawVector(V3.Y, V3.O, chroma("green").gl(), size);
        this.drawVector(V3.Z, V3.O, chroma("blue").gl(), size);
        drVs.forEach((vi) => this.drawVector(vi.v, vi.anchor, vi.color, size));
    }
    drawPlane(customPlane, color, dotted = false) {
        this.pushMatrix();
        this.multMatrix(M4.forSys(customPlane.right, customPlane.up, customPlane.normal1));
        this.translate(customPlane.uMin, customPlane.vMin, customPlane.w);
        this.scale(customPlane.uMax - customPlane.uMin, customPlane.vMax - customPlane.vMin, 1);
        const mesh = dotted
            ? this.meshes.xyDottedLinePlane
            : this.meshes.xyLinePlane;
        this.shaders.singleColor.uniforms({ color: color }).draw(mesh, this.LINES);
        this.popMatrix();
    }
    drawBox(m4, color) {
        this.pushMatrix();
        this.multMatrix(m4);
        if (color) {
            this.shaders.singleColor
                .uniforms({ color: color })
                .draw(this.meshes.cube, this.LINES);
        }
        else {
            this.shaders.multiColor.draw(this.meshes.cube, this.LINES);
        }
        this.popMatrix();
    }
}
function conicPainter(mode, gl, ellipse, color, startT, endT, width = 2) {
    gl.shaders.ellipse3d
        .uniforms({
        f1: ellipse.f1,
        f2: ellipse.f2,
        center: ellipse.center,
        color: color,
        startT: startT,
        endT: endT,
        scale: width,
        mode: mode,
    })
        .draw(gl.meshes.pipe);
}
const CURVE_PAINTERS = {
    [EllipseCurve.name]: conicPainter.bind(undefined, 0),
    [ParabolaCurve.name]: conicPainter.bind(undefined, 1),
    [HyperbolaCurve.name]: conicPainter.bind(undefined, 2),
    [ImplicitCurve.name](gl, curve, color, startT, endT, width = 2) {
        let mesh = gl.cachedMeshes.get(curve);
        const RES = 4;
        if (!mesh) {
            mesh = new Mesh()
                .addIndexBuffer("TRIANGLES")
                .addVertexBuffer("normals", "ts_Normal");
            curve.addToMesh(mesh, RES);
            mesh.compile();
            gl.cachedMeshes.set(curve, mesh);
        }
        const startIndex = ceil(startT);
        const endIndex = floor(endT);
        if (startIndex <= endIndex) {
            const indexFactor = 2 * // no of triangles per face
                RES * // no of faces
                3; // no of indexes per triangle
            gl.shaders.generic3d
                .uniforms({
                color: color,
                scale: width,
            })
                .draw(mesh, gl.TRIANGLES, startIndex * indexFactor, (floor(endT) - startIndex) * indexFactor);
            if (startT % 1 !== 0) {
                const p = curve.at(startT);
                gl.pushMatrix();
                const m = M4.forSys(p.to(curve.points[startIndex]), mesh.normals[startIndex * RES].toLength(width), mesh.normals[startIndex * RES + 1].toLength(width), p);
                gl.multMatrix(m);
                gl.shaders.singleColor
                    .uniforms({ color: color })
                    .draw(gl.meshes.pipeSegmentForICurve);
                console.log(gl.meshes.pipeSegmentForICurve);
                gl.popMatrix();
            }
            if (endT % 1 !== 0) {
                const p = curve.at(endT);
                gl.pushMatrix();
                const m = M4.forSys(curve.points[endIndex].to(p), mesh.normals[endIndex * RES].toLength(width), mesh.normals[endIndex * RES + 1].toLength(width), curve.points[endIndex]);
                gl.multMatrix(m);
                gl.shaders.singleColor
                    .uniforms({ color: color })
                    .draw(gl.meshes.pipeSegmentForICurve);
                gl.popMatrix();
            }
        }
        else {
            const p1 = curve.at(startT);
            const p2 = curve.at(endT);
            gl.pushMatrix();
            const v0 = p1.to(p2), v1 = v0.getPerpendicular().toLength(width), v2 = v0.cross(v1).toLength(width);
            const m = M4.forSys(v0, v1, v2, p1);
            gl.multMatrix(m);
            gl.shaders.singleColor
                .uniforms({ color: color })
                .draw(gl.meshes.pipeSegmentForICurve);
            gl.popMatrix();
        }
    },
    [BezierCurve.name](gl, curve, color, startT, endT, width = 2, normal = V3.Z) {
        gl.shaders.bezier3d
            .uniforms({
            p0: curve.p0,
            p1: curve.p1,
            p2: curve.p2,
            p3: curve.p3,
            color: color,
            startT: startT,
            endT: endT,
            scale: width,
            normal: normal,
        })
            .draw(gl.meshes.pipe);
    },
    [NURBS.name](gl, curve, color, startT, endT, width = 2, normal = V3.Z) {
        gl.shaders.nurbs
            .uniforms({
            "points[0]": Vector.pack(curve.points),
            degree: curve.degree,
            "knots[0]": curve.knots,
            color: color,
            startT: startT,
            endT: endT,
            scale: width,
            normal: normal,
        })
            .draw(gl.meshes.pipe);
    },
    [L3.name](gl, curve, color, startT, endT, width = 2, normal = V3.Z) {
        gl.pushMatrix();
        const a = curve.at(startT), b = curve.at(endT);
        const ab = b.minus(a), abT = ab.getPerpendicular().unit();
        const m = M4.forSys(ab, abT, ab.cross(abT).unit(), a);
        gl.multMatrix(m);
        gl.scale(1, width, width);
        gl.shaders.singleColor
            .uniforms({
            color: color,
        })
            .draw(gl.meshes.pipe);
        gl.popMatrix();
    },
};
CURVE_PAINTERS[PICurve.name] = CURVE_PAINTERS[ImplicitCurve.name];
CURVE_PAINTERS[PPCurve.name] = CURVE_PAINTERS[ImplicitCurve.name];
function initMeshes(_meshes, _gl) {
    _gl.makeCurrent();
    _meshes.cube = (() => {
        const cube = B2T.box().toMesh().addVertexBuffer("colors", "ts_Color");
        cube.colors = cube.vertices.map((p) => [p.x, p.y, p.z, 1].map((x) => x * 0.9));
        cube.compile();
        return cube;
    })();
    _meshes.sphere1 = Mesh.sphere(2);
    _meshes.segment = Mesh.plane({ startY: -0.5, height: 1, detailX: 128 });
    _meshes.text = Mesh.plane();
    _meshes.vector = Mesh.rotation([V3.O, V(0, 0.05, 0), V(0.8, 0.05), V(0.8, 0.1), V(1, 0)], L3.X, TAU, 16, true);
    _meshes.vectorShaft = Mesh.rotation([V3.O, V3.Y, V3.XY], L3.X, TAU, 8, true);
    _meshes.vectorHead = Mesh.rotation([V3.Y, V(0, 2, 0), V(2, 0, 0)], L3.X, TAU, 8, true);
    _meshes.pipe = Mesh.rotation(arrayFromFunction(512, (i, l) => new V3(i / (l - 1), -0.5, 0)), L3.X, TAU, 8, true);
    _meshes.xyLinePlane = Mesh.plane();
    _meshes.xyDottedLinePlane = makeDottedLinePlane();
    _meshes.pipeSegmentForICurve = Mesh.offsetVertices(M4.rotateY(90 * DEG).transformedPoints(arrayFromFunction(4, (i) => V3.polar(1, (TAU * i) / 4))), V3.X, true);
}
function initShaders(_gl) {
    _gl.makeCurrent();
    return {
        singleColor: Shader.create(vertexShaderBasic, fragmentShaderColor),
        multiColor: Shader.create(vertexShaderColor, fragmentShaderVaryingColor),
        singleColorHighlight: Shader.create(vertexShaderBasic, fragmentShaderColorHighlight),
        textureColor: Shader.create(vertexShaderTexture, fragmentShaderTextureColor),
        arc: Shader.create(vertexShaderRing, fragmentShaderColor),
        arc2: Shader.create(vertexShaderArc, fragmentShaderColor),
        ellipse3d: Shader.create(vertexShaderConic3d, fragmentShaderColor),
        generic3d: Shader.create(vertexShaderGeneric, fragmentShaderColor),
        bezier3d: Shader.create(vertexShaderBezier3d, fragmentShaderColor),
        nurbs: Shader.create(vertexShaderNURBS, fragmentShaderColor3),
        bezier: Shader.create(vertexShaderBezier, fragmentShaderColor),
        lighting: Shader.create(vertexShaderLighting, fragmentShaderLighting),
        waves: Shader.create(vertexShaderWaves, fragmentShaderLighting),
    };
}
function makeDottedLinePlane(count = 128) {
    const mesh = new Mesh().addIndexBuffer("LINES");
    const OXvertices = arrayFromFunction(count, (i) => new V3(i / count, 0, 0));
    mesh.vertices.push(...OXvertices);
    mesh.vertices.push(...M4.forSys(V3.Y, V3.O, V3.O, V3.X).transformedPoints(OXvertices));
    mesh.vertices.push(...M4.forSys(V3.X.negated(), V3.O, V3.O, new V3(1, 1, 0)).transformedPoints(OXvertices));
    mesh.vertices.push(...M4.forSys(V3.Y.negated(), V3.O, V3.O, V3.Y).transformedPoints(OXvertices));
    mesh.LINES = arrayFromFunction(count * 4, (i) => i - (i >= count * 2 ? 1 : 0));
    mesh.compile();
    return mesh;
}
function initNavigationEvents(_gl, eye, paintScreen) {
    const canvas = _gl.canvas;
    let lastPos = V3.O;
    //_gl.onmousedown.push((e) => {
    //	e.preventDefault()
    //	e.stopPropagation()
    //})
    //_gl.onmouseup.push((e) => {
    //	e.preventDefault()
    //	e.stopPropagation()
    //})
    canvas.addEventListener("mousemove", (e) => {
        const pagePos = V(e.pageX, e.pageY);
        const delta = lastPos.to(pagePos);
        //noinspection JSBitwiseOperatorUsage
        if (e.buttons & 4) {
            // pan
            const moveCamera = V((-delta.x * 2) / _gl.canvas.width, (delta.y * 2) / _gl.canvas.height);
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
            const rotateLR = (-delta.x / 6.0) * DEG;
            const rotateUD = (-delta.y / 6.0) * DEG;
            // rotate
            let matrix = M4.rotateLine(eye.focus, eye.up, rotateLR);
            //let horizontalRotationAxis = focus.minus(pos).cross(up)
            const horizontalRotationAxis = eye.up.cross(eye.pos.minus(eye.focus));
            matrix = matrix.times(M4.rotateLine(eye.focus, horizontalRotationAxis, rotateUD));
            eye.pos = matrix.transformPoint(eye.pos);
            eye.up = matrix.transformVector(eye.up);
            setupCamera(eye, _gl);
            paintScreen();
        }
        lastPos = pagePos;
    });
    canvas.addEventListener("wheel", (e) => {
        // zoom
        const wheelY = -sign(e.deltaY) * 2;
        // console.log(e.deltaY, e.deltaX)
        eye.zoomFactor *= pow(0.9, -wheelY);
        const mouseCoordsOnCanvas = getPosOnTarget(e);
        const mousePosFrustrum = V((mouseCoordsOnCanvas.x * 2) / _gl.canvas.offsetWidth - 1, (-mouseCoordsOnCanvas.y * 2) / _gl.canvas.offsetHeight + 1, 0);
        const moveCamera = mousePosFrustrum.times(1 - 1 / pow(0.9, -wheelY));
        const inverseProjectionMatrix = _gl.projectionMatrix.inversed();
        const worldMoveCamera = inverseProjectionMatrix.transformVector(moveCamera);
        //console.log("moveCamera", moveCamera)
        //console.log("worldMoveCamera", worldMoveCamera)
        eye.pos = eye.pos.plus(worldMoveCamera);
        eye.focus = eye.focus.plus(worldMoveCamera);
        // tilt
        const mousePosWC = inverseProjectionMatrix.transformPoint(mousePosFrustrum);
        const tiltMatrix = M4.rotateLine(mousePosWC, eye.pos.to(eye.focus), -sign(e.deltaX) * 10 * DEG);
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
    const ndc1 = V((pos.x * 2) / _gl.canvas.width - 1, (-pos.y * 2) / _gl.canvas.height + 1, 0);
    const ndc2 = V((pos.x * 2) / _gl.canvas.width - 1, (-pos.y * 2) / _gl.canvas.height + 1, 1);
    //console.log(ndc)
    const inverseProjectionMatrix = _gl.projectionMatrix.inversed();
    const s = inverseProjectionMatrix.transformPoint(ndc1);
    const dir = inverseProjectionMatrix.transformPoint(ndc2).minus(s);
    return L3.anchorDirection(s, dir);
}
function getPosOnTarget(e) {
    const target = e.target;
    const targetRect = target.getBoundingClientRect();
    const mouseCoordsOnElement = {
        x: e.clientX - targetRect.left,
        y: e.clientY - targetRect.top,
    };
    return mouseCoordsOnElement;
}
function setupCamera(_eye, _gl, suppressEvents = false) {
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
    !suppressEvents && cameraChangeListeners.forEach((l) => l(_eye));
}
const cameraChangeListeners = [];
const SHADERS_TYPE_VAR = false ;
// let shaders: typeof SHADERS_TYPE_VAR
// declare let a: BRep, b: BRep, c: BRep, d: BRep, edges: Edge[] = [], hovering: any,
// 	, normallines: boolean = false, b2s: BRep[] = []
// const

class Quaternion {
    constructor(s, x, y, z) {
        this.s = s;
        this.x = x;
        this.y = y;
        this.z = z;
    }
    static axis(axis, rotation) {
        assertf(() => axis.hasLength(1));
        return new Quaternion(cos(rotation / 2), sin(rotation / 2) * axis.x, sin(rotation / 2) * axis.y, sin(rotation / 2) * axis.z);
    }
    static of(s, x, y, z) {
        return new Quaternion(s, x, y, z);
    }
    plus(q) {
        return new Quaternion(this.s + q.s, this.x + q.x, this.y + q.y, this.z + q.z);
    }
    times(q) {
        return "number" == typeof q
            ? new Quaternion(q * this.s, q * this.x, q * this.y, q * this.z)
            : new Quaternion(this.s * q.s - (this.x * q.x + this.y * q.y + this.z * q.z), this.y * q.z - this.z * q.y + this.s * q.x + q.s * this.x, this.z * q.x - this.x * q.z + this.s * q.y + q.s * this.y, this.x * q.y - this.y * q.x + this.s * q.z + q.s * this.z);
    }
    conjugated() {
        return new Quaternion(this.s, -this.x, -this.y, -this.z);
    }
    length() {
        return Math.hypot(this.s, this.x, this.y, this.z);
    }
    norm() {
        return Math.pow(this.s, 2) + Math.pow(this.x, 2) + (Math.pow(this.y, 2) + Math.pow(this.z, 2));
    }
    unit() {
        const l = this.length();
        return new Quaternion(this.s / l, this.x / l, this.y / l, this.z / l);
    }
    inverse() {
        return this.conjugated().times(1 / this.norm());
    }
    toM4() {
        assertf(() => eq(1, this.length()));
        const { s, x, y, z } = this;
        // prettier-ignore
        return new M4([
            1 - 2 * (y * y + z * z), 2 * (x * y - z * s), 2 * (x * z + y * s), 0,
            2 * (x * y + z * s), 1 - 2 * (x * x + z * z), 2 * (y * z - x * s), 0,
            2 * (x * z - y * s), 2 * (y * z + x * s), 1 - 2 * (x * x + y * y), 0,
            0, 0, 0, 1,
        ]);
    }
    static fromRotation(m4) {
        const sqrtTracePlus1 = Math.sqrt(m4.trace() + 1);
        const f = 1 / (2 * sqrtTracePlus1);
        return new Quaternion(sqrtTracePlus1 / 2, f * (m4.e(2, 1) - m4.e(1, 2)), f * (m4.e(0, 2) - m4.e(2, 0)), f * (m4.e(1, 0) - m4.e(0, 1)));
    }
    rotatePoint(p) {
        const v = this.times(Quaternion.of(1, p.x, p.y, p.z)).times(this.conjugated());
        return new V3(v.x, v.y, v.z);
    }
    like(q, precision) {
        return (eq(this.s, q.s, precision) &&
            eq(this.x, q.x, precision) &&
            eq(this.y, q.y, precision) &&
            eq(this.z, q.z, precision));
    }
    equals(q) {
        return (this == q ||
            (q instanceof Quaternion &&
                this.s == q.s &&
                this.x == q.x &&
                this.y == q.y &&
                this.z == q.z));
    }
    hashCode() {
        let hashCode = 0;
        hashCode = (hashCode * 31 + floatHashCode(this.s)) | 0;
        hashCode = (hashCode * 31 + floatHashCode(this.x)) | 0;
        hashCode = (hashCode * 31 + floatHashCode(this.y)) | 0;
        hashCode = (hashCode * 31 + floatHashCode(this.z)) | 0;
        return hashCode;
    }
    slerp(b, f) {
        assertf(() => eq(1, this.length()));
        assertf(() => eq(1, b.length()));
        const a = this;
        let dot = a.s * b.s + a.x * b.x + a.y * b.y + a.z * b.z;
        if (dot < 0) {
            dot = -dot;
            b = b.times(-1);
            console.log("dot < 0");
        }
        const DOT_THRESHOLD = 0.9995;
        if (dot > DOT_THRESHOLD) {
            // If the inputs are too close for comfort, linearly interpolate
            // and normalize the result.
            return a
                .times(1 - f)
                .plus(b.times(f))
                .unit();
        }
        // Since dot is in range [0, DOT_THRESHOLD], acos is safe
        const theta0 = acos(dot); // theta_0 = angle between input vectors
        const theta = theta0 * f; // theta = angle between v0 and result
        const s0 = cos(theta) - (dot * sin(theta)) / sin(theta0); // == sin(theta_0 - theta) / sin(theta_0)
        const s1 = sin(theta) / sin(theta0);
        console.log(s0, s1, a.times(s0), b.times(s1));
        return a.times(s0).plus(b.times(s1));
    }
}
Quaternion.O = new Quaternion(1, 0, 0, 0);

export { AABB2, ALONG_EDGE_OR_PLANE, B2T, BREPGLContext, BRep, BezierCurve, COLORS, COPLANAR_OPPOSITE, COPLANAR_SAME, CURVE_PAINTERS, CalculateAreaVisitor, ClassSerializer, ConicSurface, Curve, CustomPlane, CylinderSurface, EPS, Edge, EllipseCurve, EllipsoidSurface, Face, FaceInfoFactory, HyperbolaCurve, INSIDE, ImplicitCurve, ImplicitSurface, L3, MathFunctionR2R, NURBS, NURBSSurface, OUTSIDE, P3, PCurveEdge, PICurve, PPCurve, ParabolaCurve, ParametricSurface, PlaneFace, PlaneSurface, PointProjectedSurface, PointVsFace, ProjectedCurveSurface, Quaternion, RotatedCurveSurface, RotationFace, SHADERS_TYPE_VAR, StraightEdge, Surface, XiEtaCurve, ZDirVolumeVisitor, addLikeSurfaceFaces, arbitraryCorner, arrayLerp, assembleFaceFromLooseEdges, breakDownPPCurves, calcNextEdgeIndex, cameraChangeListeners, cas2, cassini, createEdge, curvePoint, curvePointMF, curvePointPP, doNotSerialize, dotCurve, dotCurve2, edgeForCurveAndTs, edgePathFromSVG, edgeRect, fff, followAlgorithm2d, followAlgorithm2dAdjustable, followAlgorithmPP, getExtremePointsHelper, getGlobalId, getMouseLine, getPosOnTarget, glqArray, glqV3, initMeshes, initNavigationEvents, initShaders, intersectionCircleLine, intersectionICurveICurve, intersectionICurveICurve2, intersectionUnitCircleLine, intersectionUnitCircleLine2, intersectionUnitHyperbolaLine, ngon, parabola4Projection, parseGetParams, projectCurve, projectPointCurve, reuleaux, rotateCurve, round$1 as round, setupCamera, splitsVolumeEnclosingCone, splitsVolumeEnclosingCone2, splitsVolumeEnclosingFaces, splitsVolumeEnclosingFacesP, splitsVolumeEnclosingFacesP2, star, surfaceIsICurveIsInfosWithLine, triangulateVertices, uvInAABB2 };
//# sourceMappingURL=bundle.module.js.map
