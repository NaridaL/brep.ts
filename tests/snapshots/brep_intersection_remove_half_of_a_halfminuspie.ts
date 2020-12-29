B2T.extrudeEdges(
  [
    new StraightEdge(
      new L3(V(8, 0, 0), V(-1, 0, 0), 0, 7),
      V(8, 0, 0),
      V3.X,
      0,
      7,
    ),
    new StraightEdge(
      new L3(V3.X, V3.Y, 0, 7.937253933193771),
      V3.X,
      V(1, 7.937253933193771, 0),
      0,
      7.937253933193771,
    ),
    new PCurveEdge(
      new EllipseCurve(V3.O, V(8, 0, 0), V(0, 8, 0), 0, 3.141592653589793),
      V(1.0000000000000009, 7.937253933193771, 0),
      V(8, 0, 0),
      1.445468495626831,
      0,
      undefined,
      V(7.937253933193771, -1.0000000000000009, 0),
      V(0, -8, 0),
      undefined,
    ),
  ],
  new P3(V(0, 0, -1), 0),
  V(0, 0, 5),
  "pie/4",
)
