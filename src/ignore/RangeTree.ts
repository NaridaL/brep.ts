class RangeTree {


    nextIntervals: RangeTree = null
    prevIntervals: RangeTree = null
    intervals: {left: number, right: number}[] = []
    recIntCount: int = 0

    constructor(public s: number,
                public left: RangeTree | Number,
                public right: RangeTree | Number,) {
    }

    static fromArray(values: number[], start = 0, end = values.length): RangeTree | Number {
        if (end - start == 1) {
            return values[start]
        } else {
            const p = start + Math.floor((end - start) / 2), s = values[p - 1]
            return new RangeTree(s, RangeTree.fromArray(values, start, p), RangeTree.fromArray(values, p, end))
        }
    }

    //addIntervals(intervals) {
    //    intervals.forEach(int => addTo(this, int))
    //    function addTo(rt, int) {
    //        rt.recIntCount++
    //        if (int.right < rt.s) {
    //            addTo(rt.left, int)
    //        } else if (int.left > rt.s) {
    //            addTo(rt.right, int)
    //        } else {
    //            rt.intervals.push(int)
    //        }
    //    }
    //}

    removeInterval(interval: {left: number, right: number}) {
        let x: any
        while (true) {
            if (interval.right < x.s) {
                x.recIntCount--
                x = x.left
            } else if (interval.left > x.s) {
                x.recIntCount--
                x = x.right
            } else {
                x.recIntCount--
                x.intervals.remove(interval)
                if (0 == x.recIntCount) {
                    if (x.prev) {
                        x.prev.next = x.next
                    }
                    if (x.next) {
                        x.next.prev = x.prev
                    } else {
                        // x is last in list
                        last = x.prev
                    }
                }
                break
            }
        }
    }

    addIntervals(intervals: {left: number, right: number}[]) {
        let last
        function recurse(rt: RangeTree, intervals: {left: number, right: number}[]) {
            rt.recIntCount = intervals.length
            const intsLeft = [], intsRight = []
            for (const int of intervals) {
                (int.right < rt.s ? intsLeft : rt.s < int.left ? intsRight : rt.intervals).push(int)
            }
            intsLeft.length && recurse(rt.left as RangeTree, intsLeft)
            if (rt.intervals) {
                if (last) {
                    rt.prevIntervals = last
                    last.nextIntervals = rt
                }
                last = rt
            }
            intsRight.length && recurse(rt.right as RangeTree, intsRight)
        }
        recurse(this, intervals)
    }

    getIntervals(interval) {
        const u = this

        function getPath(a, b, u) {
            const path = []
            while (true) {
                path.push(u)
                if (b < u.s) {
                    u = u.left
                } else if (u.s < a) {
                    u = u.right
                } else {
                    break
                }
            }
            return path
        }

        function addIntersections(nodes, result) {
            nodes.forEach(node => {
                if (node.s < interval.left) {
                    node.intervals.forEach(int => {
                        if (int.right >= interval.left) {
                            result.push(int)
                        }
                    })
                } else if (interval.right < node.s) {
                    node.intervals.forEach(int => {
                        if (int.right >= interval.left) {
                            result.push(int)
                        }
                    })
                } else {
                    Array.prototype.push.apply(result, node.intervals)
                }
            })
        }

        const P1 = getPath(this, interval.left, interval.right)
        const u1 = P1.last
        const P2 = getPath(u1.left, interval.left, interval.left)
        const P3 = getPath(u1.left, interval.right, interval.right)
        const result = []
        addIntersections(P1, result)
        addIntersections(P2, result)
        addIntersections(P3, result)

    }

    toString() {
        return `${this.s} next: ${this.nextIntervals && this.nextIntervals.s}, prev: ${this.prevIntervals && this.prevIntervals.s} `
            + ` ${this.recIntCount}/` + this.intervals.toSource()
            + ('\n' + this.left + '\n' + this.right).replace(/^/gm, '\t')
    }
}