new BRep(
  [
    new PlaneFace(
      new PlaneSurface(
        new P3(V(0, -1, 0), 0),
        V3.X,
        V3.Z,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(1, 1.2246467991473532e-16, 0), V(0, 0, -1), -4096, 4096),
          V3.X,
          V(1, 0, 2),
          0,
          -2
        ),
        new StraightEdge(
          new L3(V(4, 0, 2), V(-1, 0, 0), 0, 4),
          V(1, 0, 2),
          V(0, 0, 2),
          3,
          4
        ),
        new StraightEdge(new L3(V3.O, V3.Z, 0, 2), V(0, 0, 2), V3.O, 2, 0),
        new StraightEdge(
          new L3(V(4, 0, 0), V(-1, 0, 0), 0, 4),
          V3.O,
          V3.X,
          4,
          3
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(0, -1, 0), 0),
        V3.X,
        V3.Z,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(3, 0, 0), V3.Z, -4096, 4096),
          V(3, 0, 2),
          V(3, 0, 0),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(4, 0, 0), V(-1, 0, 0), 0, 4),
          V(3, 0, 0),
          V(4, 0, 0),
          1,
          0
        ),
        new StraightEdge(
          new L3(V(4, 0, 0), V3.Z, 0, 2),
          V(4, 0, 0),
          V(4, 0, 2),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(4, 0, 2), V(-1, 0, 0), 0, 4),
          V(4, 0, 2),
          V(3, 0, 2),
          0,
          1
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(0, 0, -1), 0),
        V3.Y,
        V3.X,
        -100,
        100,
        -100,
        100
      ),
      [
        new PCurveEdge(
          new EllipseCurve(V(2, 0, 0), V(-1, 0, 0), V3.Y, 0, 3.141592653589793),
          V(3, 0, 0),
          V3.X,
          3.141592653589793,
          0,
          undefined,
          V(-1.2246467991473532e-16, 1, 0),
          V(0, -1, 0),
          "genseg1020"
        ),
        new StraightEdge(
          new L3(V(4, 0, 0), V(-1, 0, 0), 0, 4),
          V3.X,
          V3.O,
          3,
          4
        ),
        new StraightEdge(new L3(V3.O, V3.Y, 0, 2), V3.O, V(0, 2, 0), 0, 2),
        new StraightEdge(
          new L3(V(0, 2, 0), V3.X, 0, 4),
          V(0, 2, 0),
          V(4, 2, 0),
          0,
          4
        ),
        new StraightEdge(
          new L3(V(4, 2, 0), V(0, -1, 0), 0, 2),
          V(4, 2, 0),
          V(4, 0, 0),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(4, 0, 0), V(-1, 0, 0), 0, 4),
          V(4, 0, 0),
          V(3, 0, 0),
          0,
          1
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Z, 2),
        V3.Y,
        V(-1, 0, 0),
        -100,
        100,
        -100,
        100
      ),
      [
        new PCurveEdge(
          new EllipseCurve(V(2, 0, 2), V3.X, V3.Y, 0, 3.141592653589793),
          V(1, 0, 2),
          V(3, 0, 2),
          3.141592653589793,
          0,
          undefined,
          V(1.2246467991473532e-16, 1, 0),
          V(0, -1, 0),
          "genseg1023"
        ),
        new StraightEdge(
          new L3(V(4, 0, 2), V(-1, 0, 0), 0, 4),
          V(3, 0, 2),
          V(4, 0, 2),
          1,
          0
        ),
        new StraightEdge(
          new L3(V(4, 2, 2), V(0, -1, 0), 0, 2),
          V(4, 0, 2),
          V(4, 2, 2),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(0, 2, 2), V3.X, 0, 4),
          V(4, 2, 2),
          V(0, 2, 2),
          4,
          0
        ),
        new StraightEdge(
          new L3(V(0, 0, 2), V3.Y, 0, 2),
          V(0, 2, 2),
          V(0, 0, 2),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(4, 0, 2), V(-1, 0, 0), 0, 4),
          V(0, 0, 2),
          V(1, 0, 2),
          4,
          3
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(-1, 0, 0), 0),
        V(0, -1, 0),
        V3.Z,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(new L3(V3.O, V3.Y, 0, 2), V(0, 2, 0), V3.O, 2, 0),
        new StraightEdge(new L3(V3.O, V3.Z, 0, 2), V3.O, V(0, 0, 2), 0, 2),
        new StraightEdge(
          new L3(V(0, 0, 2), V3.Y, 0, 2),
          V(0, 0, 2),
          V(0, 2, 2),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(0, 2, 0), V3.Z, 0, 2),
          V(0, 2, 2),
          V(0, 2, 0),
          2,
          0
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Y, 2),
        V(-1, 0, 0),
        V3.Z,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(0, 2, 0), V3.X, 0, 4),
          V(4, 2, 0),
          V(0, 2, 0),
          4,
          0
        ),
        new StraightEdge(
          new L3(V(0, 2, 0), V3.Z, 0, 2),
          V(0, 2, 0),
          V(0, 2, 2),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(0, 2, 2), V3.X, 0, 4),
          V(0, 2, 2),
          V(4, 2, 2),
          0,
          4
        ),
        new StraightEdge(
          new L3(V(4, 2, 0), V3.Z, 0, 2),
          V(4, 2, 2),
          V(4, 2, 0),
          2,
          0
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(new P3(V3.X, 4), V3.Y, V3.Z, -100, 100, -100, 100),
      [
        new StraightEdge(
          new L3(V(4, 2, 0), V(0, -1, 0), 0, 2),
          V(4, 0, 0),
          V(4, 2, 0),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(4, 2, 0), V3.Z, 0, 2),
          V(4, 2, 0),
          V(4, 2, 2),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(4, 2, 2), V(0, -1, 0), 0, 2),
          V(4, 2, 2),
          V(4, 0, 2),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(4, 0, 0), V3.Z, 0, 2),
          V(4, 0, 2),
          V(4, 0, 0),
          2,
          0
        ),
      ],
      []
    ),
    new RotationFace(
      new CylinderSurface(
        new EllipseCurve(V(2, 0, 0), V3.X, V3.Y, 0, 3.141592653589793),
        V(0, 0, -1),
        0,
        3.141592653589793,
        -2,
        0
      ),
      [
        new StraightEdge(
          new L3(V(1, 1.2246467991473532e-16, 0), V(0, 0, -1), -4096, 4096),
          V(1, 0, 2),
          V3.X,
          -2,
          0
        ),
        new PCurveEdge(
          new EllipseCurve(V(2, 0, 0), V(-1, 0, 0), V3.Y, 0, 3.141592653589793),
          V3.X,
          V(3, 0, 0),
          0,
          3.141592653589793,
          undefined,
          V3.Y,
          V(1.2246467991473532e-16, -1, 0),
          "genseg1020"
        ),
        new StraightEdge(
          new L3(V(3, 0, 0), V3.Z, -4096, 4096),
          V(3, 0, 0),
          V(3, 0, 2),
          0,
          2
        ),
        new PCurveEdge(
          new EllipseCurve(V(2, 0, 2), V3.X, V3.Y, 0, 3.141592653589793),
          V(3, 0, 2),
          V(1, 0, 2),
          0,
          3.141592653589793,
          undefined,
          V3.Y,
          V(-1.2246467991473532e-16, -1, 0),
          "genseg1023"
        ),
      ],
      []
    ),
  ],
  false
)
