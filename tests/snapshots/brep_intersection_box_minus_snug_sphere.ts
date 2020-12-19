new BRep(
  [
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
        new StraightEdge(new L3(V3.O, V3.Y, -4096, 4096), V3.Y, V3.O, 1, 0),
        new StraightEdge(
          new L3(V3.O, V3.Z, -4096, 4096),
          V3.O,
          V(0, 0, 3),
          0,
          3
        ),
        new StraightEdge(
          new L3(V(0, 0, 3), V3.Y, -4096, 4096),
          V(0, 0, 3),
          V(0, 2, 3),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(0, 2, 0), V3.Z, -4096, 4096),
          V(0, 2, 3),
          V(0, 2, 0),
          3,
          0
        ),
        new StraightEdge(
          new L3(V3.O, V3.Y, -4096, 4096),
          V(0, 2, 0),
          V3.Y,
          2,
          1
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
          new L3(V(0, 2, 0), V3.X, -4096, 4096),
          V(1, 2, 0),
          V(0, 2, 0),
          1,
          0
        ),
        new StraightEdge(
          new L3(V(0, 2, 0), V3.Z, -4096, 4096),
          V(0, 2, 0),
          V(0, 2, 3),
          0,
          3
        ),
        new StraightEdge(
          new L3(V(0, 2, 3), V3.X, -4096, 4096),
          V(0, 2, 3),
          V(2, 2, 3),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(2, 2, 0), V3.Z, -4096, 4096),
          V(2, 2, 3),
          V(2, 2, 0),
          3,
          0
        ),
        new StraightEdge(
          new L3(V(0, 2, 0), V3.X, -4096, 4096),
          V(2, 2, 0),
          V(1, 2, 0),
          2,
          1
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(new P3(V3.X, 2), V3.Y, V3.Z, -100, 100, -100, 100),
      [
        new StraightEdge(
          new L3(V(2, 2, 0), V(0, -1, 0), -4096, 4096),
          V(2, 1, 0),
          V(2, 2, 0),
          1,
          0
        ),
        new StraightEdge(
          new L3(V(2, 2, 0), V3.Z, -4096, 4096),
          V(2, 2, 0),
          V(2, 2, 3),
          0,
          3
        ),
        new StraightEdge(
          new L3(V(2, 2, 3), V(0, -1, 0), -4096, 4096),
          V(2, 2, 3),
          V(2, 0, 3),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(2, 0, 0), V3.Z, -4096, 4096),
          V(2, 0, 3),
          V(2, 0, 0),
          3,
          0
        ),
        new StraightEdge(
          new L3(V(2, 2, 0), V(0, -1, 0), -4096, 4096),
          V(2, 0, 0),
          V(2, 1, 0),
          2,
          1
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
          new L3(V(2, 0, 0), V(-1, 0, 0), -4096, 4096),
          V3.X,
          V(2, 0, 0),
          0.9999999999999999,
          0
        ),
        new StraightEdge(
          new L3(V(2, 0, 0), V3.Z, -4096, 4096),
          V(2, 0, 0),
          V(2, 0, 3),
          0,
          3
        ),
        new StraightEdge(
          new L3(V(2, 0, 3), V(-1, 0, 0), -4096, 4096),
          V(2, 0, 3),
          V(0, 0, 3),
          0,
          2
        ),
        new StraightEdge(
          new L3(V3.O, V3.Z, -4096, 4096),
          V(0, 0, 3),
          V3.O,
          3,
          0
        ),
        new StraightEdge(
          new L3(V(2, 0, 0), V(-1, 0, 0), -4096, 4096),
          V3.O,
          V3.X,
          2,
          0.9999999999999999
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
          new EllipseCurve(V(1, 1, 0), V(-1, 0, 0), V3.Y, 0, 3.141592653589793),
          V(1, 2, 0),
          V(0, 1.0000000000000002, -6.123233995736766e-17),
          1.5707963267948966,
          2.220446049250313e-16,
          undefined,
          V(-1, -6.123233995736766e-17, 0),
          V(-2.220446049250313e-16, -1, 0),
          undefined
        ),
        new StraightEdge(
          new L3(V3.O, V3.Y, -4096, 4096),
          V3.Y,
          V(0, 2, 0),
          1,
          2
        ),
        new StraightEdge(
          new L3(V(0, 2, 0), V3.X, -4096, 4096),
          V(0, 2, 0),
          V(1, 2, 0),
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
          new EllipseCurve(V(1, 1, 0), V(-1, 0, 0), V3.Y, 0, 3.141592653589793),
          V(2, 1, -6.123233995736766e-17),
          V(1, 2, 0),
          3.141592653589793,
          1.5707963267948966,
          undefined,
          V(-1.2246467991473532e-16, 1, 0),
          V(-1, -6.123233995736766e-17, 0),
          undefined
        ),
        new StraightEdge(
          new L3(V(0, 2, 0), V3.X, -4096, 4096),
          V(1, 2, 0),
          V(2, 2, 0),
          1,
          2
        ),
        new StraightEdge(
          new L3(V(2, 2, 0), V(0, -1, 0), -4096, 4096),
          V(2, 2, 0),
          V(2, 1, 0),
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
          new EllipseCurve(
            V(1, 1, 0),
            V(1, -1.2246467991473532e-16, 0),
            V(-1.2246467991473532e-16, -1, 0),
            0,
            3.141592653589793
          ),
          V3.X,
          V(2, 1, -6.123233995736766e-17),
          1.5707963267948966,
          0,
          undefined,
          V(1, -6.123233995736766e-17, 0),
          V(1.2246467991473532e-16, 1, 0),
          undefined
        ),
        new StraightEdge(
          new L3(V(2, 2, 0), V(0, -1, 0), -4096, 4096),
          V(2, 1, 0),
          V(2, 0, 0),
          1,
          2
        ),
        new StraightEdge(
          new L3(V(2, 0, 0), V(-1, 0, 0), -4096, 4096),
          V(2, 0, 0),
          V3.X,
          0,
          0.9999999999999999
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
          new EllipseCurve(
            V(1, 1, 0),
            V(1, -1.2246467991473532e-16, 0),
            V(-1.2246467991473532e-16, -1, 0),
            0,
            3.141592653589793
          ),
          V(0, 1.0000000000000002, -6.123233995736766e-17),
          V3.X,
          3.141592653589793,
          1.5707963267948966,
          undefined,
          V(0, -1, 0),
          V(1, -6.123233995736766e-17, 0),
          undefined
        ),
        new StraightEdge(
          new L3(V(2, 0, 0), V(-1, 0, 0), -4096, 4096),
          V3.X,
          V3.O,
          0.9999999999999999,
          2
        ),
        new StraightEdge(new L3(V3.O, V3.Y, -4096, 4096), V3.O, V3.Y, 0, 1),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Z, 3),
        V3.Y,
        V(-1, 0, 0),
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(2, 0, 3), V(-1, 0, 0), -4096, 4096),
          V(0, 0, 3),
          V(2, 0, 3),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(2, 2, 3), V(0, -1, 0), -4096, 4096),
          V(2, 0, 3),
          V(2, 2, 3),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(0, 2, 3), V3.X, -4096, 4096),
          V(2, 2, 3),
          V(0, 2, 3),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(0, 0, 3), V3.Y, -4096, 4096),
          V(0, 2, 3),
          V(0, 0, 3),
          2,
          0
        ),
      ],
      []
    ),
    new RotationFace(
      new EllipsoidSurface(
        V(1, 1, 0),
        V3.X,
        V3.Y,
        V(0, 0, -1),
        0,
        3.141592653589793,
        -1.5707963267948966,
        1.5707963267948966
      ),
      [
        new PCurveEdge(
          new EllipseCurve(V(1, 1, 0), V(-1, 0, 0), V3.Y, 0, 3.141592653589793),
          V(0, 1.0000000000000002, -6.123233995736766e-17),
          V(1, 2, 0),
          2.220446049250313e-16,
          1.5707963267948966,
          undefined,
          V(2.220446049250313e-16, 1, 0),
          V(1, 6.123233995736766e-17, 0),
          undefined
        ),
        new PCurveEdge(
          new EllipseCurve(V(1, 1, 0), V(-1, 0, 0), V3.Y, 0, 3.141592653589793),
          V(1, 2, 0),
          V(2, 1, -6.123233995736766e-17),
          1.5707963267948966,
          3.141592653589793,
          undefined,
          V(1, 6.123233995736766e-17, 0),
          V(1.2246467991473532e-16, -1, 0),
          undefined
        ),
        new PCurveEdge(
          new EllipseCurve(V(1, 1, 0), V(0, 0, -1), V3.X, 0, 3.141592653589793),
          V(2, 1, -6.123233995736766e-17),
          V3.XYZ,
          1.5707963267948966,
          3.141592653589793,
          undefined,
          V(6.123233995736766e-17, 0, 1),
          V(-1, 0, 0),
          undefined
        ),
        new PCurveEdge(
          new EllipseCurve(
            V(1, 1, 0),
            V(0, 0, -1),
            V(-1, 1.2246467991473532e-16, 0),
            0,
            3.141592653589793
          ),
          V3.XYZ,
          V(0, 1.0000000000000002, -6.123233995736766e-17),
          3.141592653589793,
          1.5707963267948966,
          undefined,
          V(-1, 1.2246467991473532e-16, 0),
          V(6.123233995736766e-17, -7.498798913309288e-33, -1),
          undefined
        ),
      ],
      []
    ),
    new RotationFace(
      new EllipsoidSurface(
        V(1, 1, 0),
        V(-1, 1.2246467991473532e-16, 0),
        V(-1.2246467991473532e-16, -1, 0),
        V(0, 0, -1),
        0,
        3.141592653589793,
        -1.5707963267948966,
        1.5707963267948966
      ),
      [
        new PCurveEdge(
          new EllipseCurve(
            V(1, 1, 0),
            V(1, -1.2246467991473532e-16, 0),
            V(-1.2246467991473532e-16, -1, 0),
            0,
            3.141592653589793
          ),
          V(2, 1, -6.123233995736766e-17),
          V3.X,
          0,
          1.5707963267948966,
          undefined,
          V(-1.2246467991473532e-16, -1, 0),
          V(-1, 6.123233995736766e-17, 0),
          undefined
        ),
        new PCurveEdge(
          new EllipseCurve(
            V(1, 1, 0),
            V(1, -1.2246467991473532e-16, 0),
            V(-1.2246467991473532e-16, -1, 0),
            0,
            3.141592653589793
          ),
          V3.X,
          V(0, 1.0000000000000002, -6.123233995736766e-17),
          1.5707963267948966,
          3.141592653589793,
          undefined,
          V(-1, 6.123233995736766e-17, 0),
          V3.Y,
          undefined
        ),
        new PCurveEdge(
          new EllipseCurve(
            V(1, 1, 0),
            V(0, 0, -1),
            V(-1, 1.2246467991473532e-16, 0),
            0,
            3.141592653589793
          ),
          V(0, 1.0000000000000002, -6.123233995736766e-17),
          V3.XYZ,
          1.5707963267948966,
          3.141592653589793,
          undefined,
          V(-6.123233995736766e-17, 7.498798913309288e-33, 1),
          V(1, -1.2246467991473532e-16, 0),
          undefined
        ),
        new PCurveEdge(
          new EllipseCurve(V(1, 1, 0), V(0, 0, -1), V3.X, 0, 3.141592653589793),
          V3.XYZ,
          V(2, 1, -6.123233995736766e-17),
          3.141592653589793,
          1.5707963267948966,
          undefined,
          V3.X,
          V(-6.123233995736766e-17, 0, -1),
          undefined
        ),
      ],
      []
    ),
  ],
  false
)
