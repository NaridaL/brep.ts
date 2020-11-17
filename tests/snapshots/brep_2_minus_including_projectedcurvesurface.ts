B2T.extrudeEdges(
  [
    new StraightEdge(
      new L3(V(10, 0, 0), V(-1, 0, 0), 0, 10),
      V(10, 0, 0),
      V3.O,
      0,
      10
    ),
    new PCurveEdge(
      new BezierCurve(V3.O, V(-5, 5, 0), V(15, 5, 0), V(10, 0, 0), -0.1, 1.1),
      V3.O,
      V(10, 0, 0),
      0,
      1,
      undefined,
      V(-15, 15, 0),
      V(-15, -15, 0),
      undefined
    ),
  ],
  new P3(V(0, 0, -1), 0),
  V(0, 0, 5),
  "a/4"
).and(B2T.box(5, 10, 3, "knife").translate(1, -1, 1).flipped())
