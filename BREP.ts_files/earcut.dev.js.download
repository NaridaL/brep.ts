(function(f){if(typeof exports==="object"&&typeof module!=="undefined"){module.exports=f()}else if(typeof define==="function"&&define.amd){define([],f)}else{var g;if(typeof window!=="undefined"){g=window}else if(typeof global!=="undefined"){g=global}else if(typeof self!=="undefined"){g=self}else{g=this}g.earcut = f()}})(function(){var define,module,exports;return (function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
'use strict';

module.exports = earcut;

function earcut(data, holeIndices, dim) {

    dim = dim || 2;

    var hasHoles = holeIndices && holeIndices.length,
        outerLen = hasHoles ? holeIndices[0] * dim : data.length,
        outerNode = linkedList(data, 0, outerLen, dim, true),
        triangles = [];

    if (!outerNode) return triangles;

    var minX, minY, maxX, maxY, x, y, size;

    if (hasHoles) outerNode = eliminateHoles(data, holeIndices, outerNode, dim);

    // if the shape is not too simple, we'll use z-order curve hash later; calculate polygon bbox
    if (data.length > 80 * dim) {
        minX = maxX = data[0];
        minY = maxY = data[1];

        for (var i = dim; i < outerLen; i += dim) {
            x = data[i];
            y = data[i + 1];
            if (x < minX) minX = x;
            if (y < minY) minY = y;
            if (x > maxX) maxX = x;
            if (y > maxY) maxY = y;
        }

        // minX, minY and size are later used to transform coords into integers for z-order calculation
        size = Math.max(maxX - minX, maxY - minY);
    }

    earcutLinked(outerNode, triangles, dim, minX, minY, size);

    return triangles;
}

// create a circular doubly linked list from polygon points in the specified winding order
function linkedList(data, start, end, dim, clockwise) {
    var sum = 0,
        i, j, last;

    // calculate original winding order of a polygon ring
    for (i = start, j = end - dim; i < end; i += dim) {
        sum += (data[j] - data[i]) * (data[i + 1] + data[j + 1]);
        j = i;
    }

    // link points into circular doubly-linked list in the specified winding order
    if (clockwise === (sum > 0)) {
        for (i = start; i < end; i += dim) last = insertNode(i, data[i], data[i + 1], last);
    } else {
        for (i = end - dim; i >= start; i -= dim) last = insertNode(i, data[i], data[i + 1], last);
    }

    return last;
}

// eliminate colinear or duplicate points
function filterPoints(start, end) {
    if (!start) return start;
    if (!end) end = start;

    var p = start,
        again;
    do {
        again = false;

        if (!p.steiner && (equals(p, p.next) || area(p.prev, p, p.next) === 0)) {
            removeNode(p);
            p = end = p.prev;
            if (p === p.next) return null;
            again = true;

        } else {
            p = p.next;
        }
    } while (again || p !== end);

    return end;
}

// main ear slicing loop which triangulates a polygon (given as a linked list)
function earcutLinked(ear, triangles, dim, minX, minY, size, pass) {
    if (!ear) return;

    // interlink polygon nodes in z-order
    if (!pass && size) indexCurve(ear, minX, minY, size);

    var stop = ear,
        prev, next;

    // iterate through ears, slicing them one by one
    while (ear.prev !== ear.next) {
        prev = ear.prev;
        next = ear.next;

        if (size ? isEarHashed(ear, minX, minY, size) : isEar(ear)) {
            // cut off the triangle
            triangles.push(prev.i / dim);
            triangles.push(ear.i / dim);
            triangles.push(next.i / dim);

            removeNode(ear);

            // skipping the next vertice leads to less sliver triangles
            ear = next.next;
            stop = next.next;

            continue;
        }

        ear = next;

        // if we looped through the whole remaining polygon and can't find any more ears
        if (ear === stop) {
            // try filtering points and slicing again
            if (!pass) {
                earcutLinked(filterPoints(ear), triangles, dim, minX, minY, size, 1);

            // if this didn't work, try curing all small self-intersections locally
            } else if (pass === 1) {
                ear = cureLocalIntersections(ear, triangles, dim);
                earcutLinked(ear, triangles, dim, minX, minY, size, 2);

            // as a last resort, try splitting the remaining polygon into two
            } else if (pass === 2) {
                splitEarcut(ear, triangles, dim, minX, minY, size);
            }

            break;
        }
    }
}

// check whether a polygon node forms a valid ear with adjacent nodes
function isEar(ear) {
    var a = ear.prev,
        b = ear,
        c = ear.next;

    if (area(a, b, c) >= 0) return false; // reflex, can't be an ear

    // now make sure we don't have other points inside the potential ear
    var p = ear.next.next;

    while (p !== ear.prev) {
        if (pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) &&
            area(p.prev, p, p.next) >= 0) return false;
        p = p.next;
    }

    return true;
}

function isEarHashed(ear, minX, minY, size) {
    var a = ear.prev,
        b = ear,
        c = ear.next;

    if (area(a, b, c) >= 0) return false; // reflex, can't be an ear

    // triangle bbox; min & max are calculated like this for speed
    var minTX = a.x < b.x ? (a.x < c.x ? a.x : c.x) : (b.x < c.x ? b.x : c.x),
        minTY = a.y < b.y ? (a.y < c.y ? a.y : c.y) : (b.y < c.y ? b.y : c.y),
        maxTX = a.x > b.x ? (a.x > c.x ? a.x : c.x) : (b.x > c.x ? b.x : c.x),
        maxTY = a.y > b.y ? (a.y > c.y ? a.y : c.y) : (b.y > c.y ? b.y : c.y);

    // z-order range for the current triangle bbox;
    var minZ = zOrder(minTX, minTY, minX, minY, size),
        maxZ = zOrder(maxTX, maxTY, minX, minY, size);

    // first look for points inside the triangle in increasing z-order
    var p = ear.nextZ;

    while (p && p.z <= maxZ) {
        if (p !== ear.prev && p !== ear.next &&
            pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) &&
            area(p.prev, p, p.next) >= 0) return false;
        p = p.nextZ;
    }

    // then look for points in decreasing z-order
    p = ear.prevZ;

    while (p && p.z >= minZ) {
        if (p !== ear.prev && p !== ear.next &&
            pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) &&
            area(p.prev, p, p.next) >= 0) return false;
        p = p.prevZ;
    }

    return true;
}

// go through all polygon nodes and cure small local self-intersections
function cureLocalIntersections(start, triangles, dim) {
    var p = start;
    do {
        var a = p.prev,
            b = p.next.next;

        // a self-intersection where edge (v[i-1],v[i]) intersects (v[i+1],v[i+2])
        if (intersects(a, p, p.next, b) && locallyInside(a, b) && locallyInside(b, a)) {

            triangles.push(a.i / dim);
            triangles.push(p.i / dim);
            triangles.push(b.i / dim);

            // remove two nodes involved
            removeNode(p);
            removeNode(p.next);

            p = start = b;
        }
        p = p.next;
    } while (p !== start);

    return p;
}

// try splitting polygon into two and triangulate them independently
function splitEarcut(start, triangles, dim, minX, minY, size) {
    // look for a valid diagonal that divides the polygon into two
    var a = start;
    do {
        var b = a.next.next;
        while (b !== a.prev) {
            if (a.i !== b.i && isValidDiagonal(a, b)) {
                // split the polygon in two by the diagonal
                var c = splitPolygon(a, b);

                // filter colinear points around the cuts
                a = filterPoints(a, a.next);
                c = filterPoints(c, c.next);

                // run earcut on each half
                earcutLinked(a, triangles, dim, minX, minY, size);
                earcutLinked(c, triangles, dim, minX, minY, size);
                return;
            }
            b = b.next;
        }
        a = a.next;
    } while (a !== start);
}

// link every hole into the outer loop, producing a single-ring polygon without holes
function eliminateHoles(data, holeIndices, outerNode, dim) {
    var queue = [],
        i, len, start, end, list;

    for (i = 0, len = holeIndices.length; i < len; i++) {
        start = holeIndices[i] * dim;
        end = i < len - 1 ? holeIndices[i + 1] * dim : data.length;
        list = linkedList(data, start, end, dim, false);
        if (list === list.next) list.steiner = true;
        queue.push(getLeftmost(list));
    }

    queue.sort(compareX);

    // process holes from left to right
    for (i = 0; i < queue.length; i++) {
        eliminateHole(queue[i], outerNode);
        outerNode = filterPoints(outerNode, outerNode.next);
    }

    return outerNode;
}

function compareX(a, b) {
    return a.x - b.x;
}

// find a bridge between vertices that connects hole with an outer ring and and link it
function eliminateHole(hole, outerNode) {
    outerNode = findHoleBridge(hole, outerNode);
    if (outerNode) {
        var b = splitPolygon(outerNode, hole);
        filterPoints(b, b.next);
    }
}

// David Eberly's algorithm for finding a bridge between hole and outer polygon
function findHoleBridge(hole, outerNode) {
    var p = outerNode,
        hx = hole.x,
        hy = hole.y,
        qx = -Infinity,
        m;

    // find a segment intersected by a ray from the hole's leftmost point to the left;
    // segment's endpoint with lesser x will be potential connection point
    do {
        if (hy <= p.y && hy >= p.next.y) {
            var x = p.x + (hy - p.y) * (p.next.x - p.x) / (p.next.y - p.y);
            if (x <= hx && x > qx) {
                qx = x;
                m = p.x < p.next.x ? p : p.next;
            }
        }
        p = p.next;
    } while (p !== outerNode);

    if (!m) return null;

    // look for points inside the triangle of hole point, segment intersection and endpoint;
    // if there are no points found, we have a valid connection;
    // otherwise choose the point of the minimum angle with the ray as connection point

    var stop = m,
        tanMin = Infinity,
        tan;

    p = m.next;

    while (p !== stop) {
        if (hx >= p.x && p.x >= m.x &&
                pointInTriangle(hy < m.y ? hx : qx, hy, m.x, m.y, hy < m.y ? qx : hx, hy, p.x, p.y)) {

            tan = Math.abs(hy - p.y) / (hx - p.x); // tangential

            if ((tan < tanMin || (tan === tanMin && p.x > m.x)) && locallyInside(p, hole)) {
                m = p;
                tanMin = tan;
            }
        }

        p = p.next;
    }

    return m;
}

// interlink polygon nodes in z-order
function indexCurve(start, minX, minY, size) {
    var p = start;
    do {
        if (p.z === null) p.z = zOrder(p.x, p.y, minX, minY, size);
        p.prevZ = p.prev;
        p.nextZ = p.next;
        p = p.next;
    } while (p !== start);

    p.prevZ.nextZ = null;
    p.prevZ = null;

    sortLinked(p);
}

// Simon Tatham's linked list merge sort algorithm
// http://www.chiark.greenend.org.uk/~sgtatham/algorithms/listsort.html
function sortLinked(list) {
    var i, p, q, e, tail, numMerges, pSize, qSize,
        inSize = 1;

    do {
        p = list;
        list = null;
        tail = null;
        numMerges = 0;

        while (p) {
            numMerges++;
            q = p;
            pSize = 0;
            for (i = 0; i < inSize; i++) {
                pSize++;
                q = q.nextZ;
                if (!q) break;
            }

            qSize = inSize;

            while (pSize > 0 || (qSize > 0 && q)) {

                if (pSize === 0) {
                    e = q;
                    q = q.nextZ;
                    qSize--;
                } else if (qSize === 0 || !q) {
                    e = p;
                    p = p.nextZ;
                    pSize--;
                } else if (p.z <= q.z) {
                    e = p;
                    p = p.nextZ;
                    pSize--;
                } else {
                    e = q;
                    q = q.nextZ;
                    qSize--;
                }

                if (tail) tail.nextZ = e;
                else list = e;

                e.prevZ = tail;
                tail = e;
            }

            p = q;
        }

        tail.nextZ = null;
        inSize *= 2;

    } while (numMerges > 1);

    return list;
}

// z-order of a point given coords and size of the data bounding box
function zOrder(x, y, minX, minY, size) {
    // coords are transformed into non-negative 15-bit integer range
    x = 32767 * (x - minX) / size;
    y = 32767 * (y - minY) / size;

    x = (x | (x << 8)) & 0x00FF00FF;
    x = (x | (x << 4)) & 0x0F0F0F0F;
    x = (x | (x << 2)) & 0x33333333;
    x = (x | (x << 1)) & 0x55555555;

    y = (y | (y << 8)) & 0x00FF00FF;
    y = (y | (y << 4)) & 0x0F0F0F0F;
    y = (y | (y << 2)) & 0x33333333;
    y = (y | (y << 1)) & 0x55555555;

    return x | (y << 1);
}

// find the leftmost node of a polygon ring
function getLeftmost(start) {
    var p = start,
        leftmost = start;
    do {
        if (p.x < leftmost.x) leftmost = p;
        p = p.next;
    } while (p !== start);

    return leftmost;
}

// check if a point lies within a convex triangle
function pointInTriangle(ax, ay, bx, by, cx, cy, px, py) {
    return (cx - px) * (ay - py) - (ax - px) * (cy - py) >= 0 &&
           (ax - px) * (by - py) - (bx - px) * (ay - py) >= 0 &&
           (bx - px) * (cy - py) - (cx - px) * (by - py) >= 0;
}

// check if a diagonal between two polygon nodes is valid (lies in polygon interior)
function isValidDiagonal(a, b) {
    return equals(a, b) || a.next.i !== b.i && a.prev.i !== b.i && !intersectsPolygon(a, b) &&
           locallyInside(a, b) && locallyInside(b, a) && middleInside(a, b);
}

// signed area of a triangle
function area(p, q, r) {
    return (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
}

// check if two points are equal
function equals(p1, p2) {
    return p1.x === p2.x && p1.y === p2.y;
}

// check if two segments intersect
function intersects(p1, q1, p2, q2) {
    return area(p1, q1, p2) > 0 !== area(p1, q1, q2) > 0 &&
           area(p2, q2, p1) > 0 !== area(p2, q2, q1) > 0;
}

// check if a polygon diagonal intersects any polygon segments
function intersectsPolygon(a, b) {
    var p = a;
    do {
        if (p.i !== a.i && p.next.i !== a.i && p.i !== b.i && p.next.i !== b.i &&
                intersects(p, p.next, a, b)) return true;
        p = p.next;
    } while (p !== a);

    return false;
}

// check if a polygon diagonal is locally inside the polygon
function locallyInside(a, b) {
    return area(a.prev, a, a.next) < 0 ?
        area(a, b, a.next) >= 0 && area(a, a.prev, b) >= 0 :
        area(a, b, a.prev) < 0 || area(a, a.next, b) < 0;
}

// check if the middle point of a polygon diagonal is inside the polygon
function middleInside(a, b) {
    var p = a,
        inside = false,
        px = (a.x + b.x) / 2,
        py = (a.y + b.y) / 2;
    do {
        if (((p.y > py) !== (p.next.y > py)) && (px < (p.next.x - p.x) * (py - p.y) / (p.next.y - p.y) + p.x))
            inside = !inside;
        p = p.next;
    } while (p !== a);

    return inside;
}

// link two polygon vertices with a bridge; if the vertices belong to the same ring, it splits polygon into two;
// if one belongs to the outer ring and another to a hole, it merges it into a single ring
function splitPolygon(a, b) {
    var a2 = new Node(a.i, a.x, a.y),
        b2 = new Node(b.i, b.x, b.y),
        an = a.next,
        bp = b.prev;

    a.next = b;
    b.prev = a;

    a2.next = an;
    an.prev = a2;

    b2.next = a2;
    a2.prev = b2;

    bp.next = b2;
    b2.prev = bp;

    return b2;
}

// create a node and optionally link it with previous one (in a circular doubly linked list)
function insertNode(i, x, y, last) {
    var p = new Node(i, x, y);

    if (!last) {
        p.prev = p;
        p.next = p;

    } else {
        p.next = last.next;
        p.prev = last;
        last.next.prev = p;
        last.next = p;
    }
    return p;
}

function removeNode(p) {
    p.next.prev = p.prev;
    p.prev.next = p.next;

    if (p.prevZ) p.prevZ.nextZ = p.nextZ;
    if (p.nextZ) p.nextZ.prevZ = p.prevZ;
}

function Node(i, x, y) {
    // vertice index in coordinates array
    this.i = i;

    // vertex coordinates
    this.x = x;
    this.y = y;

    // previous and next vertice nodes in a polygon ring
    this.prev = null;
    this.next = null;

    // z-order curve value
    this.z = null;

    // previous and next nodes in z-order
    this.prevZ = null;
    this.nextZ = null;

    // indicates whether this is a steiner point
    this.steiner = false;
}

},{}]},{},[1])(1)
});
//# sourceMappingURL=data:application/json;charset=utf-8;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIm5vZGVfbW9kdWxlcy9icm93c2VyaWZ5L25vZGVfbW9kdWxlcy9icm93c2VyLXBhY2svX3ByZWx1ZGUuanMiLCJzcmMvZWFyY3V0LmpzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBO0FDQUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EiLCJmaWxlIjoiZ2VuZXJhdGVkLmpzIiwic291cmNlUm9vdCI6IiIsInNvdXJjZXNDb250ZW50IjpbIihmdW5jdGlvbiBlKHQsbixyKXtmdW5jdGlvbiBzKG8sdSl7aWYoIW5bb10pe2lmKCF0W29dKXt2YXIgYT10eXBlb2YgcmVxdWlyZT09XCJmdW5jdGlvblwiJiZyZXF1aXJlO2lmKCF1JiZhKXJldHVybiBhKG8sITApO2lmKGkpcmV0dXJuIGkobywhMCk7dmFyIGY9bmV3IEVycm9yKFwiQ2Fubm90IGZpbmQgbW9kdWxlICdcIitvK1wiJ1wiKTt0aHJvdyBmLmNvZGU9XCJNT0RVTEVfTk9UX0ZPVU5EXCIsZn12YXIgbD1uW29dPXtleHBvcnRzOnt9fTt0W29dWzBdLmNhbGwobC5leHBvcnRzLGZ1bmN0aW9uKGUpe3ZhciBuPXRbb11bMV1bZV07cmV0dXJuIHMobj9uOmUpfSxsLGwuZXhwb3J0cyxlLHQsbixyKX1yZXR1cm4gbltvXS5leHBvcnRzfXZhciBpPXR5cGVvZiByZXF1aXJlPT1cImZ1bmN0aW9uXCImJnJlcXVpcmU7Zm9yKHZhciBvPTA7bzxyLmxlbmd0aDtvKyspcyhyW29dKTtyZXR1cm4gc30pIiwiJ3VzZSBzdHJpY3QnO1xuXG5tb2R1bGUuZXhwb3J0cyA9IGVhcmN1dDtcblxuZnVuY3Rpb24gZWFyY3V0KGRhdGEsIGhvbGVJbmRpY2VzLCBkaW0pIHtcblxuICAgIGRpbSA9IGRpbSB8fCAyO1xuXG4gICAgdmFyIGhhc0hvbGVzID0gaG9sZUluZGljZXMgJiYgaG9sZUluZGljZXMubGVuZ3RoLFxuICAgICAgICBvdXRlckxlbiA9IGhhc0hvbGVzID8gaG9sZUluZGljZXNbMF0gKiBkaW0gOiBkYXRhLmxlbmd0aCxcbiAgICAgICAgb3V0ZXJOb2RlID0gbGlua2VkTGlzdChkYXRhLCAwLCBvdXRlckxlbiwgZGltLCB0cnVlKSxcbiAgICAgICAgdHJpYW5nbGVzID0gW107XG5cbiAgICBpZiAoIW91dGVyTm9kZSkgcmV0dXJuIHRyaWFuZ2xlcztcblxuICAgIHZhciBtaW5YLCBtaW5ZLCBtYXhYLCBtYXhZLCB4LCB5LCBzaXplO1xuXG4gICAgaWYgKGhhc0hvbGVzKSBvdXRlck5vZGUgPSBlbGltaW5hdGVIb2xlcyhkYXRhLCBob2xlSW5kaWNlcywgb3V0ZXJOb2RlLCBkaW0pO1xuXG4gICAgLy8gaWYgdGhlIHNoYXBlIGlzIG5vdCB0b28gc2ltcGxlLCB3ZSdsbCB1c2Ugei1vcmRlciBjdXJ2ZSBoYXNoIGxhdGVyOyBjYWxjdWxhdGUgcG9seWdvbiBiYm94XG4gICAgaWYgKGRhdGEubGVuZ3RoID4gODAgKiBkaW0pIHtcbiAgICAgICAgbWluWCA9IG1heFggPSBkYXRhWzBdO1xuICAgICAgICBtaW5ZID0gbWF4WSA9IGRhdGFbMV07XG5cbiAgICAgICAgZm9yICh2YXIgaSA9IGRpbTsgaSA8IG91dGVyTGVuOyBpICs9IGRpbSkge1xuICAgICAgICAgICAgeCA9IGRhdGFbaV07XG4gICAgICAgICAgICB5ID0gZGF0YVtpICsgMV07XG4gICAgICAgICAgICBpZiAoeCA8IG1pblgpIG1pblggPSB4O1xuICAgICAgICAgICAgaWYgKHkgPCBtaW5ZKSBtaW5ZID0geTtcbiAgICAgICAgICAgIGlmICh4ID4gbWF4WCkgbWF4WCA9IHg7XG4gICAgICAgICAgICBpZiAoeSA+IG1heFkpIG1heFkgPSB5O1xuICAgICAgICB9XG5cbiAgICAgICAgLy8gbWluWCwgbWluWSBhbmQgc2l6ZSBhcmUgbGF0ZXIgdXNlZCB0byB0cmFuc2Zvcm0gY29vcmRzIGludG8gaW50ZWdlcnMgZm9yIHotb3JkZXIgY2FsY3VsYXRpb25cbiAgICAgICAgc2l6ZSA9IE1hdGgubWF4KG1heFggLSBtaW5YLCBtYXhZIC0gbWluWSk7XG4gICAgfVxuXG4gICAgZWFyY3V0TGlua2VkKG91dGVyTm9kZSwgdHJpYW5nbGVzLCBkaW0sIG1pblgsIG1pblksIHNpemUpO1xuXG4gICAgcmV0dXJuIHRyaWFuZ2xlcztcbn1cblxuLy8gY3JlYXRlIGEgY2lyY3VsYXIgZG91Ymx5IGxpbmtlZCBsaXN0IGZyb20gcG9seWdvbiBwb2ludHMgaW4gdGhlIHNwZWNpZmllZCB3aW5kaW5nIG9yZGVyXG5mdW5jdGlvbiBsaW5rZWRMaXN0KGRhdGEsIHN0YXJ0LCBlbmQsIGRpbSwgY2xvY2t3aXNlKSB7XG4gICAgdmFyIHN1bSA9IDAsXG4gICAgICAgIGksIGosIGxhc3Q7XG5cbiAgICAvLyBjYWxjdWxhdGUgb3JpZ2luYWwgd2luZGluZyBvcmRlciBvZiBhIHBvbHlnb24gcmluZ1xuICAgIGZvciAoaSA9IHN0YXJ0LCBqID0gZW5kIC0gZGltOyBpIDwgZW5kOyBpICs9IGRpbSkge1xuICAgICAgICBzdW0gKz0gKGRhdGFbal0gLSBkYXRhW2ldKSAqIChkYXRhW2kgKyAxXSArIGRhdGFbaiArIDFdKTtcbiAgICAgICAgaiA9IGk7XG4gICAgfVxuXG4gICAgLy8gbGluayBwb2ludHMgaW50byBjaXJjdWxhciBkb3VibHktbGlua2VkIGxpc3QgaW4gdGhlIHNwZWNpZmllZCB3aW5kaW5nIG9yZGVyXG4gICAgaWYgKGNsb2Nrd2lzZSA9PT0gKHN1bSA+IDApKSB7XG4gICAgICAgIGZvciAoaSA9IHN0YXJ0OyBpIDwgZW5kOyBpICs9IGRpbSkgbGFzdCA9IGluc2VydE5vZGUoaSwgZGF0YVtpXSwgZGF0YVtpICsgMV0sIGxhc3QpO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIGZvciAoaSA9IGVuZCAtIGRpbTsgaSA+PSBzdGFydDsgaSAtPSBkaW0pIGxhc3QgPSBpbnNlcnROb2RlKGksIGRhdGFbaV0sIGRhdGFbaSArIDFdLCBsYXN0KTtcbiAgICB9XG5cbiAgICByZXR1cm4gbGFzdDtcbn1cblxuLy8gZWxpbWluYXRlIGNvbGluZWFyIG9yIGR1cGxpY2F0ZSBwb2ludHNcbmZ1bmN0aW9uIGZpbHRlclBvaW50cyhzdGFydCwgZW5kKSB7XG4gICAgaWYgKCFzdGFydCkgcmV0dXJuIHN0YXJ0O1xuICAgIGlmICghZW5kKSBlbmQgPSBzdGFydDtcblxuICAgIHZhciBwID0gc3RhcnQsXG4gICAgICAgIGFnYWluO1xuICAgIGRvIHtcbiAgICAgICAgYWdhaW4gPSBmYWxzZTtcblxuICAgICAgICBpZiAoIXAuc3RlaW5lciAmJiAoZXF1YWxzKHAsIHAubmV4dCkgfHwgYXJlYShwLnByZXYsIHAsIHAubmV4dCkgPT09IDApKSB7XG4gICAgICAgICAgICByZW1vdmVOb2RlKHApO1xuICAgICAgICAgICAgcCA9IGVuZCA9IHAucHJldjtcbiAgICAgICAgICAgIGlmIChwID09PSBwLm5leHQpIHJldHVybiBudWxsO1xuICAgICAgICAgICAgYWdhaW4gPSB0cnVlO1xuXG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBwID0gcC5uZXh0O1xuICAgICAgICB9XG4gICAgfSB3aGlsZSAoYWdhaW4gfHwgcCAhPT0gZW5kKTtcblxuICAgIHJldHVybiBlbmQ7XG59XG5cbi8vIG1haW4gZWFyIHNsaWNpbmcgbG9vcCB3aGljaCB0cmlhbmd1bGF0ZXMgYSBwb2x5Z29uIChnaXZlbiBhcyBhIGxpbmtlZCBsaXN0KVxuZnVuY3Rpb24gZWFyY3V0TGlua2VkKGVhciwgdHJpYW5nbGVzLCBkaW0sIG1pblgsIG1pblksIHNpemUsIHBhc3MpIHtcbiAgICBpZiAoIWVhcikgcmV0dXJuO1xuXG4gICAgLy8gaW50ZXJsaW5rIHBvbHlnb24gbm9kZXMgaW4gei1vcmRlclxuICAgIGlmICghcGFzcyAmJiBzaXplKSBpbmRleEN1cnZlKGVhciwgbWluWCwgbWluWSwgc2l6ZSk7XG5cbiAgICB2YXIgc3RvcCA9IGVhcixcbiAgICAgICAgcHJldiwgbmV4dDtcblxuICAgIC8vIGl0ZXJhdGUgdGhyb3VnaCBlYXJzLCBzbGljaW5nIHRoZW0gb25lIGJ5IG9uZVxuICAgIHdoaWxlIChlYXIucHJldiAhPT0gZWFyLm5leHQpIHtcbiAgICAgICAgcHJldiA9IGVhci5wcmV2O1xuICAgICAgICBuZXh0ID0gZWFyLm5leHQ7XG5cbiAgICAgICAgaWYgKHNpemUgPyBpc0Vhckhhc2hlZChlYXIsIG1pblgsIG1pblksIHNpemUpIDogaXNFYXIoZWFyKSkge1xuICAgICAgICAgICAgLy8gY3V0IG9mZiB0aGUgdHJpYW5nbGVcbiAgICAgICAgICAgIHRyaWFuZ2xlcy5wdXNoKHByZXYuaSAvIGRpbSk7XG4gICAgICAgICAgICB0cmlhbmdsZXMucHVzaChlYXIuaSAvIGRpbSk7XG4gICAgICAgICAgICB0cmlhbmdsZXMucHVzaChuZXh0LmkgLyBkaW0pO1xuXG4gICAgICAgICAgICByZW1vdmVOb2RlKGVhcik7XG5cbiAgICAgICAgICAgIC8vIHNraXBwaW5nIHRoZSBuZXh0IHZlcnRpY2UgbGVhZHMgdG8gbGVzcyBzbGl2ZXIgdHJpYW5nbGVzXG4gICAgICAgICAgICBlYXIgPSBuZXh0Lm5leHQ7XG4gICAgICAgICAgICBzdG9wID0gbmV4dC5uZXh0O1xuXG4gICAgICAgICAgICBjb250aW51ZTtcbiAgICAgICAgfVxuXG4gICAgICAgIGVhciA9IG5leHQ7XG5cbiAgICAgICAgLy8gaWYgd2UgbG9vcGVkIHRocm91Z2ggdGhlIHdob2xlIHJlbWFpbmluZyBwb2x5Z29uIGFuZCBjYW4ndCBmaW5kIGFueSBtb3JlIGVhcnNcbiAgICAgICAgaWYgKGVhciA9PT0gc3RvcCkge1xuICAgICAgICAgICAgLy8gdHJ5IGZpbHRlcmluZyBwb2ludHMgYW5kIHNsaWNpbmcgYWdhaW5cbiAgICAgICAgICAgIGlmICghcGFzcykge1xuICAgICAgICAgICAgICAgIGVhcmN1dExpbmtlZChmaWx0ZXJQb2ludHMoZWFyKSwgdHJpYW5nbGVzLCBkaW0sIG1pblgsIG1pblksIHNpemUsIDEpO1xuXG4gICAgICAgICAgICAvLyBpZiB0aGlzIGRpZG4ndCB3b3JrLCB0cnkgY3VyaW5nIGFsbCBzbWFsbCBzZWxmLWludGVyc2VjdGlvbnMgbG9jYWxseVxuICAgICAgICAgICAgfSBlbHNlIGlmIChwYXNzID09PSAxKSB7XG4gICAgICAgICAgICAgICAgZWFyID0gY3VyZUxvY2FsSW50ZXJzZWN0aW9ucyhlYXIsIHRyaWFuZ2xlcywgZGltKTtcbiAgICAgICAgICAgICAgICBlYXJjdXRMaW5rZWQoZWFyLCB0cmlhbmdsZXMsIGRpbSwgbWluWCwgbWluWSwgc2l6ZSwgMik7XG5cbiAgICAgICAgICAgIC8vIGFzIGEgbGFzdCByZXNvcnQsIHRyeSBzcGxpdHRpbmcgdGhlIHJlbWFpbmluZyBwb2x5Z29uIGludG8gdHdvXG4gICAgICAgICAgICB9IGVsc2UgaWYgKHBhc3MgPT09IDIpIHtcbiAgICAgICAgICAgICAgICBzcGxpdEVhcmN1dChlYXIsIHRyaWFuZ2xlcywgZGltLCBtaW5YLCBtaW5ZLCBzaXplKTtcbiAgICAgICAgICAgIH1cblxuICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgIH1cbiAgICB9XG59XG5cbi8vIGNoZWNrIHdoZXRoZXIgYSBwb2x5Z29uIG5vZGUgZm9ybXMgYSB2YWxpZCBlYXIgd2l0aCBhZGphY2VudCBub2Rlc1xuZnVuY3Rpb24gaXNFYXIoZWFyKSB7XG4gICAgdmFyIGEgPSBlYXIucHJldixcbiAgICAgICAgYiA9IGVhcixcbiAgICAgICAgYyA9IGVhci5uZXh0O1xuXG4gICAgaWYgKGFyZWEoYSwgYiwgYykgPj0gMCkgcmV0dXJuIGZhbHNlOyAvLyByZWZsZXgsIGNhbid0IGJlIGFuIGVhclxuXG4gICAgLy8gbm93IG1ha2Ugc3VyZSB3ZSBkb24ndCBoYXZlIG90aGVyIHBvaW50cyBpbnNpZGUgdGhlIHBvdGVudGlhbCBlYXJcbiAgICB2YXIgcCA9IGVhci5uZXh0Lm5leHQ7XG5cbiAgICB3aGlsZSAocCAhPT0gZWFyLnByZXYpIHtcbiAgICAgICAgaWYgKHBvaW50SW5UcmlhbmdsZShhLngsIGEueSwgYi54LCBiLnksIGMueCwgYy55LCBwLngsIHAueSkgJiZcbiAgICAgICAgICAgIGFyZWEocC5wcmV2LCBwLCBwLm5leHQpID49IDApIHJldHVybiBmYWxzZTtcbiAgICAgICAgcCA9IHAubmV4dDtcbiAgICB9XG5cbiAgICByZXR1cm4gdHJ1ZTtcbn1cblxuZnVuY3Rpb24gaXNFYXJIYXNoZWQoZWFyLCBtaW5YLCBtaW5ZLCBzaXplKSB7XG4gICAgdmFyIGEgPSBlYXIucHJldixcbiAgICAgICAgYiA9IGVhcixcbiAgICAgICAgYyA9IGVhci5uZXh0O1xuXG4gICAgaWYgKGFyZWEoYSwgYiwgYykgPj0gMCkgcmV0dXJuIGZhbHNlOyAvLyByZWZsZXgsIGNhbid0IGJlIGFuIGVhclxuXG4gICAgLy8gdHJpYW5nbGUgYmJveDsgbWluICYgbWF4IGFyZSBjYWxjdWxhdGVkIGxpa2UgdGhpcyBmb3Igc3BlZWRcbiAgICB2YXIgbWluVFggPSBhLnggPCBiLnggPyAoYS54IDwgYy54ID8gYS54IDogYy54KSA6IChiLnggPCBjLnggPyBiLnggOiBjLngpLFxuICAgICAgICBtaW5UWSA9IGEueSA8IGIueSA/IChhLnkgPCBjLnkgPyBhLnkgOiBjLnkpIDogKGIueSA8IGMueSA/IGIueSA6IGMueSksXG4gICAgICAgIG1heFRYID0gYS54ID4gYi54ID8gKGEueCA+IGMueCA/IGEueCA6IGMueCkgOiAoYi54ID4gYy54ID8gYi54IDogYy54KSxcbiAgICAgICAgbWF4VFkgPSBhLnkgPiBiLnkgPyAoYS55ID4gYy55ID8gYS55IDogYy55KSA6IChiLnkgPiBjLnkgPyBiLnkgOiBjLnkpO1xuXG4gICAgLy8gei1vcmRlciByYW5nZSBmb3IgdGhlIGN1cnJlbnQgdHJpYW5nbGUgYmJveDtcbiAgICB2YXIgbWluWiA9IHpPcmRlcihtaW5UWCwgbWluVFksIG1pblgsIG1pblksIHNpemUpLFxuICAgICAgICBtYXhaID0gek9yZGVyKG1heFRYLCBtYXhUWSwgbWluWCwgbWluWSwgc2l6ZSk7XG5cbiAgICAvLyBmaXJzdCBsb29rIGZvciBwb2ludHMgaW5zaWRlIHRoZSB0cmlhbmdsZSBpbiBpbmNyZWFzaW5nIHotb3JkZXJcbiAgICB2YXIgcCA9IGVhci5uZXh0WjtcblxuICAgIHdoaWxlIChwICYmIHAueiA8PSBtYXhaKSB7XG4gICAgICAgIGlmIChwICE9PSBlYXIucHJldiAmJiBwICE9PSBlYXIubmV4dCAmJlxuICAgICAgICAgICAgcG9pbnRJblRyaWFuZ2xlKGEueCwgYS55LCBiLngsIGIueSwgYy54LCBjLnksIHAueCwgcC55KSAmJlxuICAgICAgICAgICAgYXJlYShwLnByZXYsIHAsIHAubmV4dCkgPj0gMCkgcmV0dXJuIGZhbHNlO1xuICAgICAgICBwID0gcC5uZXh0WjtcbiAgICB9XG5cbiAgICAvLyB0aGVuIGxvb2sgZm9yIHBvaW50cyBpbiBkZWNyZWFzaW5nIHotb3JkZXJcbiAgICBwID0gZWFyLnByZXZaO1xuXG4gICAgd2hpbGUgKHAgJiYgcC56ID49IG1pblopIHtcbiAgICAgICAgaWYgKHAgIT09IGVhci5wcmV2ICYmIHAgIT09IGVhci5uZXh0ICYmXG4gICAgICAgICAgICBwb2ludEluVHJpYW5nbGUoYS54LCBhLnksIGIueCwgYi55LCBjLngsIGMueSwgcC54LCBwLnkpICYmXG4gICAgICAgICAgICBhcmVhKHAucHJldiwgcCwgcC5uZXh0KSA+PSAwKSByZXR1cm4gZmFsc2U7XG4gICAgICAgIHAgPSBwLnByZXZaO1xuICAgIH1cblxuICAgIHJldHVybiB0cnVlO1xufVxuXG4vLyBnbyB0aHJvdWdoIGFsbCBwb2x5Z29uIG5vZGVzIGFuZCBjdXJlIHNtYWxsIGxvY2FsIHNlbGYtaW50ZXJzZWN0aW9uc1xuZnVuY3Rpb24gY3VyZUxvY2FsSW50ZXJzZWN0aW9ucyhzdGFydCwgdHJpYW5nbGVzLCBkaW0pIHtcbiAgICB2YXIgcCA9IHN0YXJ0O1xuICAgIGRvIHtcbiAgICAgICAgdmFyIGEgPSBwLnByZXYsXG4gICAgICAgICAgICBiID0gcC5uZXh0Lm5leHQ7XG5cbiAgICAgICAgLy8gYSBzZWxmLWludGVyc2VjdGlvbiB3aGVyZSBlZGdlICh2W2ktMV0sdltpXSkgaW50ZXJzZWN0cyAodltpKzFdLHZbaSsyXSlcbiAgICAgICAgaWYgKGludGVyc2VjdHMoYSwgcCwgcC5uZXh0LCBiKSAmJiBsb2NhbGx5SW5zaWRlKGEsIGIpICYmIGxvY2FsbHlJbnNpZGUoYiwgYSkpIHtcblxuICAgICAgICAgICAgdHJpYW5nbGVzLnB1c2goYS5pIC8gZGltKTtcbiAgICAgICAgICAgIHRyaWFuZ2xlcy5wdXNoKHAuaSAvIGRpbSk7XG4gICAgICAgICAgICB0cmlhbmdsZXMucHVzaChiLmkgLyBkaW0pO1xuXG4gICAgICAgICAgICAvLyByZW1vdmUgdHdvIG5vZGVzIGludm9sdmVkXG4gICAgICAgICAgICByZW1vdmVOb2RlKHApO1xuICAgICAgICAgICAgcmVtb3ZlTm9kZShwLm5leHQpO1xuXG4gICAgICAgICAgICBwID0gc3RhcnQgPSBiO1xuICAgICAgICB9XG4gICAgICAgIHAgPSBwLm5leHQ7XG4gICAgfSB3aGlsZSAocCAhPT0gc3RhcnQpO1xuXG4gICAgcmV0dXJuIHA7XG59XG5cbi8vIHRyeSBzcGxpdHRpbmcgcG9seWdvbiBpbnRvIHR3byBhbmQgdHJpYW5ndWxhdGUgdGhlbSBpbmRlcGVuZGVudGx5XG5mdW5jdGlvbiBzcGxpdEVhcmN1dChzdGFydCwgdHJpYW5nbGVzLCBkaW0sIG1pblgsIG1pblksIHNpemUpIHtcbiAgICAvLyBsb29rIGZvciBhIHZhbGlkIGRpYWdvbmFsIHRoYXQgZGl2aWRlcyB0aGUgcG9seWdvbiBpbnRvIHR3b1xuICAgIHZhciBhID0gc3RhcnQ7XG4gICAgZG8ge1xuICAgICAgICB2YXIgYiA9IGEubmV4dC5uZXh0O1xuICAgICAgICB3aGlsZSAoYiAhPT0gYS5wcmV2KSB7XG4gICAgICAgICAgICBpZiAoYS5pICE9PSBiLmkgJiYgaXNWYWxpZERpYWdvbmFsKGEsIGIpKSB7XG4gICAgICAgICAgICAgICAgLy8gc3BsaXQgdGhlIHBvbHlnb24gaW4gdHdvIGJ5IHRoZSBkaWFnb25hbFxuICAgICAgICAgICAgICAgIHZhciBjID0gc3BsaXRQb2x5Z29uKGEsIGIpO1xuXG4gICAgICAgICAgICAgICAgLy8gZmlsdGVyIGNvbGluZWFyIHBvaW50cyBhcm91bmQgdGhlIGN1dHNcbiAgICAgICAgICAgICAgICBhID0gZmlsdGVyUG9pbnRzKGEsIGEubmV4dCk7XG4gICAgICAgICAgICAgICAgYyA9IGZpbHRlclBvaW50cyhjLCBjLm5leHQpO1xuXG4gICAgICAgICAgICAgICAgLy8gcnVuIGVhcmN1dCBvbiBlYWNoIGhhbGZcbiAgICAgICAgICAgICAgICBlYXJjdXRMaW5rZWQoYSwgdHJpYW5nbGVzLCBkaW0sIG1pblgsIG1pblksIHNpemUpO1xuICAgICAgICAgICAgICAgIGVhcmN1dExpbmtlZChjLCB0cmlhbmdsZXMsIGRpbSwgbWluWCwgbWluWSwgc2l6ZSk7XG4gICAgICAgICAgICAgICAgcmV0dXJuO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgYiA9IGIubmV4dDtcbiAgICAgICAgfVxuICAgICAgICBhID0gYS5uZXh0O1xuICAgIH0gd2hpbGUgKGEgIT09IHN0YXJ0KTtcbn1cblxuLy8gbGluayBldmVyeSBob2xlIGludG8gdGhlIG91dGVyIGxvb3AsIHByb2R1Y2luZyBhIHNpbmdsZS1yaW5nIHBvbHlnb24gd2l0aG91dCBob2xlc1xuZnVuY3Rpb24gZWxpbWluYXRlSG9sZXMoZGF0YSwgaG9sZUluZGljZXMsIG91dGVyTm9kZSwgZGltKSB7XG4gICAgdmFyIHF1ZXVlID0gW10sXG4gICAgICAgIGksIGxlbiwgc3RhcnQsIGVuZCwgbGlzdDtcblxuICAgIGZvciAoaSA9IDAsIGxlbiA9IGhvbGVJbmRpY2VzLmxlbmd0aDsgaSA8IGxlbjsgaSsrKSB7XG4gICAgICAgIHN0YXJ0ID0gaG9sZUluZGljZXNbaV0gKiBkaW07XG4gICAgICAgIGVuZCA9IGkgPCBsZW4gLSAxID8gaG9sZUluZGljZXNbaSArIDFdICogZGltIDogZGF0YS5sZW5ndGg7XG4gICAgICAgIGxpc3QgPSBsaW5rZWRMaXN0KGRhdGEsIHN0YXJ0LCBlbmQsIGRpbSwgZmFsc2UpO1xuICAgICAgICBpZiAobGlzdCA9PT0gbGlzdC5uZXh0KSBsaXN0LnN0ZWluZXIgPSB0cnVlO1xuICAgICAgICBxdWV1ZS5wdXNoKGdldExlZnRtb3N0KGxpc3QpKTtcbiAgICB9XG5cbiAgICBxdWV1ZS5zb3J0KGNvbXBhcmVYKTtcblxuICAgIC8vIHByb2Nlc3MgaG9sZXMgZnJvbSBsZWZ0IHRvIHJpZ2h0XG4gICAgZm9yIChpID0gMDsgaSA8IHF1ZXVlLmxlbmd0aDsgaSsrKSB7XG4gICAgICAgIGVsaW1pbmF0ZUhvbGUocXVldWVbaV0sIG91dGVyTm9kZSk7XG4gICAgICAgIG91dGVyTm9kZSA9IGZpbHRlclBvaW50cyhvdXRlck5vZGUsIG91dGVyTm9kZS5uZXh0KTtcbiAgICB9XG5cbiAgICByZXR1cm4gb3V0ZXJOb2RlO1xufVxuXG5mdW5jdGlvbiBjb21wYXJlWChhLCBiKSB7XG4gICAgcmV0dXJuIGEueCAtIGIueDtcbn1cblxuLy8gZmluZCBhIGJyaWRnZSBiZXR3ZWVuIHZlcnRpY2VzIHRoYXQgY29ubmVjdHMgaG9sZSB3aXRoIGFuIG91dGVyIHJpbmcgYW5kIGFuZCBsaW5rIGl0XG5mdW5jdGlvbiBlbGltaW5hdGVIb2xlKGhvbGUsIG91dGVyTm9kZSkge1xuICAgIG91dGVyTm9kZSA9IGZpbmRIb2xlQnJpZGdlKGhvbGUsIG91dGVyTm9kZSk7XG4gICAgaWYgKG91dGVyTm9kZSkge1xuICAgICAgICB2YXIgYiA9IHNwbGl0UG9seWdvbihvdXRlck5vZGUsIGhvbGUpO1xuICAgICAgICBmaWx0ZXJQb2ludHMoYiwgYi5uZXh0KTtcbiAgICB9XG59XG5cbi8vIERhdmlkIEViZXJseSdzIGFsZ29yaXRobSBmb3IgZmluZGluZyBhIGJyaWRnZSBiZXR3ZWVuIGhvbGUgYW5kIG91dGVyIHBvbHlnb25cbmZ1bmN0aW9uIGZpbmRIb2xlQnJpZGdlKGhvbGUsIG91dGVyTm9kZSkge1xuICAgIHZhciBwID0gb3V0ZXJOb2RlLFxuICAgICAgICBoeCA9IGhvbGUueCxcbiAgICAgICAgaHkgPSBob2xlLnksXG4gICAgICAgIHF4ID0gLUluZmluaXR5LFxuICAgICAgICBtO1xuXG4gICAgLy8gZmluZCBhIHNlZ21lbnQgaW50ZXJzZWN0ZWQgYnkgYSByYXkgZnJvbSB0aGUgaG9sZSdzIGxlZnRtb3N0IHBvaW50IHRvIHRoZSBsZWZ0O1xuICAgIC8vIHNlZ21lbnQncyBlbmRwb2ludCB3aXRoIGxlc3NlciB4IHdpbGwgYmUgcG90ZW50aWFsIGNvbm5lY3Rpb24gcG9pbnRcbiAgICBkbyB7XG4gICAgICAgIGlmIChoeSA8PSBwLnkgJiYgaHkgPj0gcC5uZXh0LnkpIHtcbiAgICAgICAgICAgIHZhciB4ID0gcC54ICsgKGh5IC0gcC55KSAqIChwLm5leHQueCAtIHAueCkgLyAocC5uZXh0LnkgLSBwLnkpO1xuICAgICAgICAgICAgaWYgKHggPD0gaHggJiYgeCA+IHF4KSB7XG4gICAgICAgICAgICAgICAgcXggPSB4O1xuICAgICAgICAgICAgICAgIG0gPSBwLnggPCBwLm5leHQueCA/IHAgOiBwLm5leHQ7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgcCA9IHAubmV4dDtcbiAgICB9IHdoaWxlIChwICE9PSBvdXRlck5vZGUpO1xuXG4gICAgaWYgKCFtKSByZXR1cm4gbnVsbDtcblxuICAgIC8vIGxvb2sgZm9yIHBvaW50cyBpbnNpZGUgdGhlIHRyaWFuZ2xlIG9mIGhvbGUgcG9pbnQsIHNlZ21lbnQgaW50ZXJzZWN0aW9uIGFuZCBlbmRwb2ludDtcbiAgICAvLyBpZiB0aGVyZSBhcmUgbm8gcG9pbnRzIGZvdW5kLCB3ZSBoYXZlIGEgdmFsaWQgY29ubmVjdGlvbjtcbiAgICAvLyBvdGhlcndpc2UgY2hvb3NlIHRoZSBwb2ludCBvZiB0aGUgbWluaW11bSBhbmdsZSB3aXRoIHRoZSByYXkgYXMgY29ubmVjdGlvbiBwb2ludFxuXG4gICAgdmFyIHN0b3AgPSBtLFxuICAgICAgICB0YW5NaW4gPSBJbmZpbml0eSxcbiAgICAgICAgdGFuO1xuXG4gICAgcCA9IG0ubmV4dDtcblxuICAgIHdoaWxlIChwICE9PSBzdG9wKSB7XG4gICAgICAgIGlmIChoeCA+PSBwLnggJiYgcC54ID49IG0ueCAmJlxuICAgICAgICAgICAgICAgIHBvaW50SW5UcmlhbmdsZShoeSA8IG0ueSA/IGh4IDogcXgsIGh5LCBtLngsIG0ueSwgaHkgPCBtLnkgPyBxeCA6IGh4LCBoeSwgcC54LCBwLnkpKSB7XG5cbiAgICAgICAgICAgIHRhbiA9IE1hdGguYWJzKGh5IC0gcC55KSAvIChoeCAtIHAueCk7IC8vIHRhbmdlbnRpYWxcblxuICAgICAgICAgICAgaWYgKCh0YW4gPCB0YW5NaW4gfHwgKHRhbiA9PT0gdGFuTWluICYmIHAueCA+IG0ueCkpICYmIGxvY2FsbHlJbnNpZGUocCwgaG9sZSkpIHtcbiAgICAgICAgICAgICAgICBtID0gcDtcbiAgICAgICAgICAgICAgICB0YW5NaW4gPSB0YW47XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICBwID0gcC5uZXh0O1xuICAgIH1cblxuICAgIHJldHVybiBtO1xufVxuXG4vLyBpbnRlcmxpbmsgcG9seWdvbiBub2RlcyBpbiB6LW9yZGVyXG5mdW5jdGlvbiBpbmRleEN1cnZlKHN0YXJ0LCBtaW5YLCBtaW5ZLCBzaXplKSB7XG4gICAgdmFyIHAgPSBzdGFydDtcbiAgICBkbyB7XG4gICAgICAgIGlmIChwLnogPT09IG51bGwpIHAueiA9IHpPcmRlcihwLngsIHAueSwgbWluWCwgbWluWSwgc2l6ZSk7XG4gICAgICAgIHAucHJldlogPSBwLnByZXY7XG4gICAgICAgIHAubmV4dFogPSBwLm5leHQ7XG4gICAgICAgIHAgPSBwLm5leHQ7XG4gICAgfSB3aGlsZSAocCAhPT0gc3RhcnQpO1xuXG4gICAgcC5wcmV2Wi5uZXh0WiA9IG51bGw7XG4gICAgcC5wcmV2WiA9IG51bGw7XG5cbiAgICBzb3J0TGlua2VkKHApO1xufVxuXG4vLyBTaW1vbiBUYXRoYW0ncyBsaW5rZWQgbGlzdCBtZXJnZSBzb3J0IGFsZ29yaXRobVxuLy8gaHR0cDovL3d3dy5jaGlhcmsuZ3JlZW5lbmQub3JnLnVrL35zZ3RhdGhhbS9hbGdvcml0aG1zL2xpc3Rzb3J0Lmh0bWxcbmZ1bmN0aW9uIHNvcnRMaW5rZWQobGlzdCkge1xuICAgIHZhciBpLCBwLCBxLCBlLCB0YWlsLCBudW1NZXJnZXMsIHBTaXplLCBxU2l6ZSxcbiAgICAgICAgaW5TaXplID0gMTtcblxuICAgIGRvIHtcbiAgICAgICAgcCA9IGxpc3Q7XG4gICAgICAgIGxpc3QgPSBudWxsO1xuICAgICAgICB0YWlsID0gbnVsbDtcbiAgICAgICAgbnVtTWVyZ2VzID0gMDtcblxuICAgICAgICB3aGlsZSAocCkge1xuICAgICAgICAgICAgbnVtTWVyZ2VzKys7XG4gICAgICAgICAgICBxID0gcDtcbiAgICAgICAgICAgIHBTaXplID0gMDtcbiAgICAgICAgICAgIGZvciAoaSA9IDA7IGkgPCBpblNpemU7IGkrKykge1xuICAgICAgICAgICAgICAgIHBTaXplKys7XG4gICAgICAgICAgICAgICAgcSA9IHEubmV4dFo7XG4gICAgICAgICAgICAgICAgaWYgKCFxKSBicmVhaztcbiAgICAgICAgICAgIH1cblxuICAgICAgICAgICAgcVNpemUgPSBpblNpemU7XG5cbiAgICAgICAgICAgIHdoaWxlIChwU2l6ZSA+IDAgfHwgKHFTaXplID4gMCAmJiBxKSkge1xuXG4gICAgICAgICAgICAgICAgaWYgKHBTaXplID09PSAwKSB7XG4gICAgICAgICAgICAgICAgICAgIGUgPSBxO1xuICAgICAgICAgICAgICAgICAgICBxID0gcS5uZXh0WjtcbiAgICAgICAgICAgICAgICAgICAgcVNpemUtLTtcbiAgICAgICAgICAgICAgICB9IGVsc2UgaWYgKHFTaXplID09PSAwIHx8ICFxKSB7XG4gICAgICAgICAgICAgICAgICAgIGUgPSBwO1xuICAgICAgICAgICAgICAgICAgICBwID0gcC5uZXh0WjtcbiAgICAgICAgICAgICAgICAgICAgcFNpemUtLTtcbiAgICAgICAgICAgICAgICB9IGVsc2UgaWYgKHAueiA8PSBxLnopIHtcbiAgICAgICAgICAgICAgICAgICAgZSA9IHA7XG4gICAgICAgICAgICAgICAgICAgIHAgPSBwLm5leHRaO1xuICAgICAgICAgICAgICAgICAgICBwU2l6ZS0tO1xuICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgIGUgPSBxO1xuICAgICAgICAgICAgICAgICAgICBxID0gcS5uZXh0WjtcbiAgICAgICAgICAgICAgICAgICAgcVNpemUtLTtcbiAgICAgICAgICAgICAgICB9XG5cbiAgICAgICAgICAgICAgICBpZiAodGFpbCkgdGFpbC5uZXh0WiA9IGU7XG4gICAgICAgICAgICAgICAgZWxzZSBsaXN0ID0gZTtcblxuICAgICAgICAgICAgICAgIGUucHJldlogPSB0YWlsO1xuICAgICAgICAgICAgICAgIHRhaWwgPSBlO1xuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICBwID0gcTtcbiAgICAgICAgfVxuXG4gICAgICAgIHRhaWwubmV4dFogPSBudWxsO1xuICAgICAgICBpblNpemUgKj0gMjtcblxuICAgIH0gd2hpbGUgKG51bU1lcmdlcyA+IDEpO1xuXG4gICAgcmV0dXJuIGxpc3Q7XG59XG5cbi8vIHotb3JkZXIgb2YgYSBwb2ludCBnaXZlbiBjb29yZHMgYW5kIHNpemUgb2YgdGhlIGRhdGEgYm91bmRpbmcgYm94XG5mdW5jdGlvbiB6T3JkZXIoeCwgeSwgbWluWCwgbWluWSwgc2l6ZSkge1xuICAgIC8vIGNvb3JkcyBhcmUgdHJhbnNmb3JtZWQgaW50byBub24tbmVnYXRpdmUgMTUtYml0IGludGVnZXIgcmFuZ2VcbiAgICB4ID0gMzI3NjcgKiAoeCAtIG1pblgpIC8gc2l6ZTtcbiAgICB5ID0gMzI3NjcgKiAoeSAtIG1pblkpIC8gc2l6ZTtcblxuICAgIHggPSAoeCB8ICh4IDw8IDgpKSAmIDB4MDBGRjAwRkY7XG4gICAgeCA9ICh4IHwgKHggPDwgNCkpICYgMHgwRjBGMEYwRjtcbiAgICB4ID0gKHggfCAoeCA8PCAyKSkgJiAweDMzMzMzMzMzO1xuICAgIHggPSAoeCB8ICh4IDw8IDEpKSAmIDB4NTU1NTU1NTU7XG5cbiAgICB5ID0gKHkgfCAoeSA8PCA4KSkgJiAweDAwRkYwMEZGO1xuICAgIHkgPSAoeSB8ICh5IDw8IDQpKSAmIDB4MEYwRjBGMEY7XG4gICAgeSA9ICh5IHwgKHkgPDwgMikpICYgMHgzMzMzMzMzMztcbiAgICB5ID0gKHkgfCAoeSA8PCAxKSkgJiAweDU1NTU1NTU1O1xuXG4gICAgcmV0dXJuIHggfCAoeSA8PCAxKTtcbn1cblxuLy8gZmluZCB0aGUgbGVmdG1vc3Qgbm9kZSBvZiBhIHBvbHlnb24gcmluZ1xuZnVuY3Rpb24gZ2V0TGVmdG1vc3Qoc3RhcnQpIHtcbiAgICB2YXIgcCA9IHN0YXJ0LFxuICAgICAgICBsZWZ0bW9zdCA9IHN0YXJ0O1xuICAgIGRvIHtcbiAgICAgICAgaWYgKHAueCA8IGxlZnRtb3N0LngpIGxlZnRtb3N0ID0gcDtcbiAgICAgICAgcCA9IHAubmV4dDtcbiAgICB9IHdoaWxlIChwICE9PSBzdGFydCk7XG5cbiAgICByZXR1cm4gbGVmdG1vc3Q7XG59XG5cbi8vIGNoZWNrIGlmIGEgcG9pbnQgbGllcyB3aXRoaW4gYSBjb252ZXggdHJpYW5nbGVcbmZ1bmN0aW9uIHBvaW50SW5UcmlhbmdsZShheCwgYXksIGJ4LCBieSwgY3gsIGN5LCBweCwgcHkpIHtcbiAgICByZXR1cm4gKGN4IC0gcHgpICogKGF5IC0gcHkpIC0gKGF4IC0gcHgpICogKGN5IC0gcHkpID49IDAgJiZcbiAgICAgICAgICAgKGF4IC0gcHgpICogKGJ5IC0gcHkpIC0gKGJ4IC0gcHgpICogKGF5IC0gcHkpID49IDAgJiZcbiAgICAgICAgICAgKGJ4IC0gcHgpICogKGN5IC0gcHkpIC0gKGN4IC0gcHgpICogKGJ5IC0gcHkpID49IDA7XG59XG5cbi8vIGNoZWNrIGlmIGEgZGlhZ29uYWwgYmV0d2VlbiB0d28gcG9seWdvbiBub2RlcyBpcyB2YWxpZCAobGllcyBpbiBwb2x5Z29uIGludGVyaW9yKVxuZnVuY3Rpb24gaXNWYWxpZERpYWdvbmFsKGEsIGIpIHtcbiAgICByZXR1cm4gZXF1YWxzKGEsIGIpIHx8IGEubmV4dC5pICE9PSBiLmkgJiYgYS5wcmV2LmkgIT09IGIuaSAmJiAhaW50ZXJzZWN0c1BvbHlnb24oYSwgYikgJiZcbiAgICAgICAgICAgbG9jYWxseUluc2lkZShhLCBiKSAmJiBsb2NhbGx5SW5zaWRlKGIsIGEpICYmIG1pZGRsZUluc2lkZShhLCBiKTtcbn1cblxuLy8gc2lnbmVkIGFyZWEgb2YgYSB0cmlhbmdsZVxuZnVuY3Rpb24gYXJlYShwLCBxLCByKSB7XG4gICAgcmV0dXJuIChxLnkgLSBwLnkpICogKHIueCAtIHEueCkgLSAocS54IC0gcC54KSAqIChyLnkgLSBxLnkpO1xufVxuXG4vLyBjaGVjayBpZiB0d28gcG9pbnRzIGFyZSBlcXVhbFxuZnVuY3Rpb24gZXF1YWxzKHAxLCBwMikge1xuICAgIHJldHVybiBwMS54ID09PSBwMi54ICYmIHAxLnkgPT09IHAyLnk7XG59XG5cbi8vIGNoZWNrIGlmIHR3byBzZWdtZW50cyBpbnRlcnNlY3RcbmZ1bmN0aW9uIGludGVyc2VjdHMocDEsIHExLCBwMiwgcTIpIHtcbiAgICByZXR1cm4gYXJlYShwMSwgcTEsIHAyKSA+IDAgIT09IGFyZWEocDEsIHExLCBxMikgPiAwICYmXG4gICAgICAgICAgIGFyZWEocDIsIHEyLCBwMSkgPiAwICE9PSBhcmVhKHAyLCBxMiwgcTEpID4gMDtcbn1cblxuLy8gY2hlY2sgaWYgYSBwb2x5Z29uIGRpYWdvbmFsIGludGVyc2VjdHMgYW55IHBvbHlnb24gc2VnbWVudHNcbmZ1bmN0aW9uIGludGVyc2VjdHNQb2x5Z29uKGEsIGIpIHtcbiAgICB2YXIgcCA9IGE7XG4gICAgZG8ge1xuICAgICAgICBpZiAocC5pICE9PSBhLmkgJiYgcC5uZXh0LmkgIT09IGEuaSAmJiBwLmkgIT09IGIuaSAmJiBwLm5leHQuaSAhPT0gYi5pICYmXG4gICAgICAgICAgICAgICAgaW50ZXJzZWN0cyhwLCBwLm5leHQsIGEsIGIpKSByZXR1cm4gdHJ1ZTtcbiAgICAgICAgcCA9IHAubmV4dDtcbiAgICB9IHdoaWxlIChwICE9PSBhKTtcblxuICAgIHJldHVybiBmYWxzZTtcbn1cblxuLy8gY2hlY2sgaWYgYSBwb2x5Z29uIGRpYWdvbmFsIGlzIGxvY2FsbHkgaW5zaWRlIHRoZSBwb2x5Z29uXG5mdW5jdGlvbiBsb2NhbGx5SW5zaWRlKGEsIGIpIHtcbiAgICByZXR1cm4gYXJlYShhLnByZXYsIGEsIGEubmV4dCkgPCAwID9cbiAgICAgICAgYXJlYShhLCBiLCBhLm5leHQpID49IDAgJiYgYXJlYShhLCBhLnByZXYsIGIpID49IDAgOlxuICAgICAgICBhcmVhKGEsIGIsIGEucHJldikgPCAwIHx8IGFyZWEoYSwgYS5uZXh0LCBiKSA8IDA7XG59XG5cbi8vIGNoZWNrIGlmIHRoZSBtaWRkbGUgcG9pbnQgb2YgYSBwb2x5Z29uIGRpYWdvbmFsIGlzIGluc2lkZSB0aGUgcG9seWdvblxuZnVuY3Rpb24gbWlkZGxlSW5zaWRlKGEsIGIpIHtcbiAgICB2YXIgcCA9IGEsXG4gICAgICAgIGluc2lkZSA9IGZhbHNlLFxuICAgICAgICBweCA9IChhLnggKyBiLngpIC8gMixcbiAgICAgICAgcHkgPSAoYS55ICsgYi55KSAvIDI7XG4gICAgZG8ge1xuICAgICAgICBpZiAoKChwLnkgPiBweSkgIT09IChwLm5leHQueSA+IHB5KSkgJiYgKHB4IDwgKHAubmV4dC54IC0gcC54KSAqIChweSAtIHAueSkgLyAocC5uZXh0LnkgLSBwLnkpICsgcC54KSlcbiAgICAgICAgICAgIGluc2lkZSA9ICFpbnNpZGU7XG4gICAgICAgIHAgPSBwLm5leHQ7XG4gICAgfSB3aGlsZSAocCAhPT0gYSk7XG5cbiAgICByZXR1cm4gaW5zaWRlO1xufVxuXG4vLyBsaW5rIHR3byBwb2x5Z29uIHZlcnRpY2VzIHdpdGggYSBicmlkZ2U7IGlmIHRoZSB2ZXJ0aWNlcyBiZWxvbmcgdG8gdGhlIHNhbWUgcmluZywgaXQgc3BsaXRzIHBvbHlnb24gaW50byB0d287XG4vLyBpZiBvbmUgYmVsb25ncyB0byB0aGUgb3V0ZXIgcmluZyBhbmQgYW5vdGhlciB0byBhIGhvbGUsIGl0IG1lcmdlcyBpdCBpbnRvIGEgc2luZ2xlIHJpbmdcbmZ1bmN0aW9uIHNwbGl0UG9seWdvbihhLCBiKSB7XG4gICAgdmFyIGEyID0gbmV3IE5vZGUoYS5pLCBhLngsIGEueSksXG4gICAgICAgIGIyID0gbmV3IE5vZGUoYi5pLCBiLngsIGIueSksXG4gICAgICAgIGFuID0gYS5uZXh0LFxuICAgICAgICBicCA9IGIucHJldjtcblxuICAgIGEubmV4dCA9IGI7XG4gICAgYi5wcmV2ID0gYTtcblxuICAgIGEyLm5leHQgPSBhbjtcbiAgICBhbi5wcmV2ID0gYTI7XG5cbiAgICBiMi5uZXh0ID0gYTI7XG4gICAgYTIucHJldiA9IGIyO1xuXG4gICAgYnAubmV4dCA9IGIyO1xuICAgIGIyLnByZXYgPSBicDtcblxuICAgIHJldHVybiBiMjtcbn1cblxuLy8gY3JlYXRlIGEgbm9kZSBhbmQgb3B0aW9uYWxseSBsaW5rIGl0IHdpdGggcHJldmlvdXMgb25lIChpbiBhIGNpcmN1bGFyIGRvdWJseSBsaW5rZWQgbGlzdClcbmZ1bmN0aW9uIGluc2VydE5vZGUoaSwgeCwgeSwgbGFzdCkge1xuICAgIHZhciBwID0gbmV3IE5vZGUoaSwgeCwgeSk7XG5cbiAgICBpZiAoIWxhc3QpIHtcbiAgICAgICAgcC5wcmV2ID0gcDtcbiAgICAgICAgcC5uZXh0ID0gcDtcblxuICAgIH0gZWxzZSB7XG4gICAgICAgIHAubmV4dCA9IGxhc3QubmV4dDtcbiAgICAgICAgcC5wcmV2ID0gbGFzdDtcbiAgICAgICAgbGFzdC5uZXh0LnByZXYgPSBwO1xuICAgICAgICBsYXN0Lm5leHQgPSBwO1xuICAgIH1cbiAgICByZXR1cm4gcDtcbn1cblxuZnVuY3Rpb24gcmVtb3ZlTm9kZShwKSB7XG4gICAgcC5uZXh0LnByZXYgPSBwLnByZXY7XG4gICAgcC5wcmV2Lm5leHQgPSBwLm5leHQ7XG5cbiAgICBpZiAocC5wcmV2WikgcC5wcmV2Wi5uZXh0WiA9IHAubmV4dFo7XG4gICAgaWYgKHAubmV4dFopIHAubmV4dFoucHJldlogPSBwLnByZXZaO1xufVxuXG5mdW5jdGlvbiBOb2RlKGksIHgsIHkpIHtcbiAgICAvLyB2ZXJ0aWNlIGluZGV4IGluIGNvb3JkaW5hdGVzIGFycmF5XG4gICAgdGhpcy5pID0gaTtcblxuICAgIC8vIHZlcnRleCBjb29yZGluYXRlc1xuICAgIHRoaXMueCA9IHg7XG4gICAgdGhpcy55ID0geTtcblxuICAgIC8vIHByZXZpb3VzIGFuZCBuZXh0IHZlcnRpY2Ugbm9kZXMgaW4gYSBwb2x5Z29uIHJpbmdcbiAgICB0aGlzLnByZXYgPSBudWxsO1xuICAgIHRoaXMubmV4dCA9IG51bGw7XG5cbiAgICAvLyB6LW9yZGVyIGN1cnZlIHZhbHVlXG4gICAgdGhpcy56ID0gbnVsbDtcblxuICAgIC8vIHByZXZpb3VzIGFuZCBuZXh0IG5vZGVzIGluIHotb3JkZXJcbiAgICB0aGlzLnByZXZaID0gbnVsbDtcbiAgICB0aGlzLm5leHRaID0gbnVsbDtcblxuICAgIC8vIGluZGljYXRlcyB3aGV0aGVyIHRoaXMgaXMgYSBzdGVpbmVyIHBvaW50XG4gICAgdGhpcy5zdGVpbmVyID0gZmFsc2U7XG59XG4iXX0=
