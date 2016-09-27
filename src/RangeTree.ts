
function RangeTree(vals, start, end) {
	start = start || 0
	end = end || vals.length
	var m = (end - start)
	if (m == 1) {
		return Number(vals[start])
	} else {
		var p = start + Math.floor(m / 2), s = (vals[p-1] + vals[p]) / 2
		this.s = s
		this.left = new RangeTree(vals, start, p)
		this.right = new RangeTree(vals, p, end)
		this.ints = []
		this.containsIntervals = 0
		this.next = this.prev = null
	}
}
RangeTree.prototype.addIntervals = function (intervals) {
	intervals.forEach(int => addTo(this, int))
	function addTo(x, int) {
		if (int.right < x.s) {
			x.containsIntervals++
			addTo(x.left, int)
		} else if (int.left > x.s) {
			x.containsIntervals++
			addTo(x.right, int)
		} else {
			x.containsIntervals++
			x.ints.push(int)
		}
	}
}
var last = null
RangeTree.prototype.removeInterval = function (interval) {
	while (true) {
		if (interval.right < x.s) {
			x.containsIntervals--
			x = x.left
		} else if (interval.left > x.s) {
			x.containsIntervals--
			x = x.right
		} else {
			x.containsIntervals--
			x.ints.remove(interval)
			if (0 == x.containsIntervals) {
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
RangeTree.prototype.addInterval = function (interval) {
	var x = this
	function findRightOfRecursiveInOrder (x, s) {
		if (!x.containsIntervals) return null
		var res
		if (x.s > s && (res = findRightOfRecursiveInOrder(x.left, s))) {
			return res
		}
		if (x.s > s && 0 != x.ints.length) return x
		if (res = findRightOfRecursiveInOrder(x.right, s)) {
			return res
		}

	}
	while (true) {
		if (interval.right < x.s) {
			x.containsIntervals++
			x = x.left
		} else if (interval.left > x.s) {
			x.containsIntervals++
			x = x.right
		} else {
			x.containsIntervals++
			x.ints.push(interval)
			if (1 == x.containsIntervals) {
				var nextNode = findRightOfRecursiveInOrder(this, x.s)
				if (nextNode) {
					if (nextNode.prev) {
						x.prev = nextNode.prev
						x.prev.next = x
					}
					nextNode.prev = x
					x.next = nextNode
				} else {
					// this is the (new) last node in the list
					if (last && last != x) {
						x.prev = last
						last.next = x
					}
					last = x
				}
			}
			break
		}
	}
}
RangeTree.prototype.getIntervals = function (interval) {
	var u = this
	function getPath(a, b, u) {
		var P = []
		while (true) {
			P.push(u)
			if (b < u.s) {
				u = u.left
			} else if (u.s < a) {
				u = u.right
			} else {
				break
			}
		}
		return P
	}
	function addIntersections(nodes, result) {
		nodes.forEach(node => {
			if (node.s < interval.left) {
				node.ints.forEach(int => {
					if (int.right >= interval.left) {
						result.push(int)
					}
				})
			} else if (interval.right < node.s) {
				node.ints.forEach(int => {
					if (int.right >= interval.left) {
						result.push(int)
					}
				})
			} else {
				Array.prototype.push.apply(result, node.ints)
			}
		})
	}
	var P1 = getPath(this, interval.left, interval.right)
	var P2 = getPath(u1.left, interval.left, interval.left)
	var P3 = getPath(u1.left, interval.right, interval.right)
	var result = []
	addIntersections(P1, result)
	addIntersections(P2, result)
	addIntersections(P3, result)

}
RangeTree.prototype.toString = function (a) {
	a = a || 0
	return `${this.s} next: ${this.next && this.next.s}, prev: ${this.prev && this.prev.s} `
		+` ${this.containsIntervals}/`+ this.ints.toSource()
		+ "\n" + NLA.repeatString(a + 1, "  ") + (Number.isFinite(this.left) ? this.left : this.left.toString(a + 1))
		+ "\n" + NLA.repeatString(a + 1, "  ") + (Number.isFinite(this.right) ? this.right : this.right.toString(a + 1))
}
/*
 var rt = new RangeTree([-2, -1, 0,2,3,4,7-1,7]);
 rt.addInterval({left: -1, right: 6})
 rt.addInterval({left: -2, right: 3})
 rt.addInterval({left: 0, right: 4})
 rt.addInterval({left: 2, right: 7})
 rt.addInterval({left: 3, right: 4})
 rt.addInterval({left: -2, right: -1})
 console.log(rt.toString())
 NLA.addTransformationMethods(AABB.prototype)
 AABB.intersections = function (aabbs) {

 }
 */