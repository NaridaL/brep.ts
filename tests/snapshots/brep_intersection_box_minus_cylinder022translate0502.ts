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
          new L3(V(0.5, 0, 0), V(0, 0, -1), -4096, 4096),
          V(0.5, 0, 0),
          V(0.5, 0, 1),
          0,
          -1
        ),
        new StraightEdge(
          new L3(V(1, 0, 1), V(-1, 0, 0), -4096, 4096),
          V(0.5, 0, 1),
          V3.Z,
          0.5,
          1
        ),
        new StraightEdge(new L3(V3.O, V3.Z, -4096, 4096), V3.Z, V3.O, 1, 0),
        new StraightEdge(
          new L3(V3.X, V(-1, 0, 0), -4096, 4096),
          V3.O,
          V(0.5, 0, 0),
          1,
          0.5
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
          new L3(V(0.5, 0, 0), V(0, 0, -1), -4096, 4096),
          V(0.5, 0, 1),
          V(0.5, 0, 0),
          -1,
          0
        ),
        new StraightEdge(
          new L3(V3.X, V(-1, 0, 0), -4096, 4096),
          V(0.5, 0, 0),
          V3.X,
          0.5,
          0
        ),
        new StraightEdge(
          new L3(V3.X, V3.Z, -4096, 4096),
          V3.X,
          V(1, 0, 1),
          0,
          1
        ),
        new StraightEdge(
          new L3(V(1, 0, 1), V(-1, 0, 0), -4096, 4096),
          V(1, 0, 1),
          V(0.5, 0, 1),
          0,
          0.5
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
            V(0.5, 0.2, 0),
            V(0.2, 0, 0),
            V(0, 0.2, 0),
            0,
            3.141592653589793
          ),
          V(0.7, 0.2, 0),
          V(0.3, 0.2, 0),
          0,
          3.141592653589793,
          undefined,
          V(0, 0.2, 0),
          V(0, -0.2, 0),
          "genseg91"
        ),
        new PCurveEdge(
          new EllipseCurve(
            V(0.5, 0.2, 0),
            V(-0.2, 0, 0),
            V(0, -0.2, 0),
            0,
            3.141592653589793
          ),
          V(0.3, 0.2, 0),
          V(0.5, 0, 0),
          0,
          1.5707963267948966,
          undefined,
          V(0, -0.2, 0),
          V(0.2, 0, 0),
          "genseg92"
        ),
        new StraightEdge(
          new L3(V3.X, V(-1, 0, 0), -4096, 4096),
          V(0.5, 0, 0),
          V3.O,
          0.5,
          1
        ),
        new StraightEdge(new L3(V3.O, V3.Y, -4096, 4096), V3.O, V3.Y, 0, 1),
        new StraightEdge(
          new L3(V3.Y, V3.X, -4096, 4096),
          V3.Y,
          V(1, 1, 0),
          0,
          1
        ),
        new StraightEdge(
          new L3(V(1, 1, 0), V(0, -1, 0), -4096, 4096),
          V(1, 1, 0),
          V3.X,
          0,
          1
        ),
        new StraightEdge(
          new L3(V3.X, V(-1, 0, 0), -4096, 4096),
          V3.X,
          V(0.5, 0, 0),
          0,
          0.5
        ),
        new PCurveEdge(
          new EllipseCurve(
            V(0.5, 0.2, 0),
            V(-0.2, 0, 0),
            V(0, -0.2, 0),
            0,
            3.141592653589793
          ),
          V(0.5, 0, 0),
          V(0.7, 0.2, 0),
          1.5707963267948966,
          3.141592653589793,
          undefined,
          V(0.2, 0, 0),
          V(0, 0.2, 0),
          "genseg93"
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Z, 1),
        V3.Y,
        V(-1, 0, 0),
        -100,
        100,
        -100,
        100
      ),
      [
        new PCurveEdge(
          new EllipseCurve(
            V(0.5, 0.2, 1),
            V(-0.2, 0, 0),
            V(0, 0.2, 0),
            0,
            3.141592653589793
          ),
          V(0.3, 0.2, 1),
          V(0.7, 0.2, 1),
          0,
          3.141592653589793,
          undefined,
          V(0, 0.2, 0),
          V(0, -0.2, 0),
          "genseg94"
        ),
        new PCurveEdge(
          new EllipseCurve(
            V(0.5, 0.2, 1),
            V(0.2, 0, 0),
            V(0, -0.2, 0),
            0,
            3.141592653589793
          ),
          V(0.7, 0.2, 1),
          V(0.5, 0, 1),
          0,
          1.5707963267948966,
          undefined,
          V(0, -0.2, 0),
          V(-0.2, 0, 0),
          "genseg95"
        ),
        new StraightEdge(
          new L3(V(1, 0, 1), V(-1, 0, 0), -4096, 4096),
          V(0.5, 0, 1),
          V(1, 0, 1),
          0.5,
          0
        ),
        new StraightEdge(
          new L3(V3.XYZ, V(0, -1, 0), -4096, 4096),
          V(1, 0, 1),
          V3.XYZ,
          1,
          0
        ),
        new StraightEdge(
          new L3(V(0, 1, 1), V3.X, -4096, 4096),
          V3.XYZ,
          V(0, 1, 1),
          1,
          0
        ),
        new StraightEdge(
          new L3(V3.Z, V3.Y, -4096, 4096),
          V(0, 1, 1),
          V3.Z,
          1,
          0
        ),
        new StraightEdge(
          new L3(V(1, 0, 1), V(-1, 0, 0), -4096, 4096),
          V3.Z,
          V(0.5, 0, 1),
          1,
          0.5
        ),
        new PCurveEdge(
          new EllipseCurve(
            V(0.5, 0.2, 1),
            V(0.2, 0, 0),
            V(0, -0.2, 0),
            0,
            3.141592653589793
          ),
          V(0.5, 0, 1),
          V(0.3, 0.2, 1),
          1.5707963267948966,
          3.141592653589793,
          undefined,
          V(-0.2, 0, 0),
          V(0, 0.2, 0),
          "genseg96"
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
        new StraightEdge(new L3(V3.O, V3.Y, -4096, 4096), V3.Y, V3.O, 1, 0),
        new StraightEdge(new L3(V3.O, V3.Z, -4096, 4096), V3.O, V3.Z, 0, 1),
        new StraightEdge(
          new L3(V3.Z, V3.Y, -4096, 4096),
          V3.Z,
          V(0, 1, 1),
          0,
          1
        ),
        new StraightEdge(
          new L3(V3.Y, V3.Z, -4096, 4096),
          V(0, 1, 1),
          V3.Y,
          1,
          0
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Y, 1),
        V(-1, 0, 0),
        V3.Z,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V3.Y, V3.X, -4096, 4096),
          V(1, 1, 0),
          V3.Y,
          1,
          0
        ),
        new StraightEdge(
          new L3(V3.Y, V3.Z, -4096, 4096),
          V3.Y,
          V(0, 1, 1),
          0,
          1
        ),
        new StraightEdge(
          new L3(V(0, 1, 1), V3.X, -4096, 4096),
          V(0, 1, 1),
          V3.XYZ,
          0,
          1
        ),
        new StraightEdge(
          new L3(V(1, 1, 0), V3.Z, -4096, 4096),
          V3.XYZ,
          V(1, 1, 0),
          1,
          0
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(new P3(V3.X, 1), V3.Y, V3.Z, -100, 100, -100, 100),
      [
        new StraightEdge(
          new L3(V(1, 1, 0), V(0, -1, 0), -4096, 4096),
          V3.X,
          V(1, 1, 0),
          1,
          0
        ),
        new StraightEdge(
          new L3(V(1, 1, 0), V3.Z, -4096, 4096),
          V(1, 1, 0),
          V3.XYZ,
          0,
          1
        ),
        new StraightEdge(
          new L3(V3.XYZ, V(0, -1, 0), -4096, 4096),
          V3.XYZ,
          V(1, 0, 1),
          0,
          1
        ),
        new StraightEdge(
          new L3(V3.X, V3.Z, -4096, 4096),
          V(1, 0, 1),
          V3.X,
          1,
          0
        ),
      ],
      []
    ),
    new RotationFace(
      new CylinderSurface(
        new EllipseCurve(
          V(0.5, 0.2, 0),
          V(0.2, 0, 0),
          V(0, 0.2, 0),
          0,
          3.141592653589793
        ),
        V(0, 0, -1),
        0,
        3.141592653589793,
        -2,
        0
      ),
      [
        new PCurveEdge(
          new EllipseCurve(
            V(0.5, 0.2, 0),
            V(0.2, 0, 0),
            V(0, 0.2, 0),
            0,
            3.141592653589793
          ),
          V(0.3, 0.2, 0),
          V(0.7, 0.2, 0),
          3.141592653589793,
          0,
          undefined,
          V(0, 0.2, 0),
          V(0, -0.2, 0),
          "genseg91"
        ),
        new StraightEdge(
          new L3(V(0.7, 0.2, 0), V3.Z, -4096, 4096),
          V(0.7, 0.2, 0),
          V(0.7, 0.2, 1),
          0,
          1
        ),
        new PCurveEdge(
          new EllipseCurve(
            V(0.5, 0.2, 1),
            V(-0.2, 0, 0),
            V(0, 0.2, 0),
            0,
            3.141592653589793
          ),
          V(0.7, 0.2, 1),
          V(0.3, 0.2, 1),
          3.141592653589793,
          0,
          undefined,
          V(0, 0.2, 0),
          V(0, -0.2, 0),
          "genseg94"
        ),
        new StraightEdge(
          new L3(V(0.3, 0.2, 0), V3.Z, -4096, 4096),
          V(0.3, 0.2, 1),
          V(0.3, 0.2, 0),
          1,
          0
        ),
      ],
      []
    ),
    new RotationFace(
      new CylinderSurface(
        new EllipseCurve(
          V(0.5, 0.2, 0),
          V(-0.2, 0, 0),
          V(0, -0.2, 0),
          0,
          3.141592653589793
        ),
        V(0, 0, -1),
        0,
        3.141592653589793,
        -2,
        0
      ),
      [
        new StraightEdge(
          new L3(V(0.5, 0, 0), V(0, 0, -1), -4096, 4096),
          V(0.5, 0, 1),
          V(0.5, 0, 0),
          -1,
          0
        ),
        new PCurveEdge(
          new EllipseCurve(
            V(0.5, 0.2, 0),
            V(-0.2, 0, 0),
            V(0, -0.2, 0),
            0,
            3.141592653589793
          ),
          V(0.5, 0, 0),
          V(0.3, 0.2, 0),
          1.5707963267948966,
          0,
          undefined,
          V(-0.2, 0, 0),
          V(0, 0.2, 0),
          "genseg92"
        ),
        new StraightEdge(
          new L3(V(0.3, 0.2, 0), V3.Z, -4096, 4096),
          V(0.3, 0.2, 0),
          V(0.3, 0.2, 1),
          0,
          1
        ),
        new PCurveEdge(
          new EllipseCurve(
            V(0.5, 0.2, 1),
            V(0.2, 0, 0),
            V(0, -0.2, 0),
            0,
            3.141592653589793
          ),
          V(0.3, 0.2, 1),
          V(0.5, 0, 1),
          3.141592653589793,
          1.5707963267948966,
          undefined,
          V(0, -0.2, 0),
          V(0.2, 0, 0),
          "genseg96"
        ),
      ],
      []
    ),
    new RotationFace(
      new CylinderSurface(
        new EllipseCurve(
          V(0.5, 0.2, 0),
          V(-0.2, 0, 0),
          V(0, -0.2, 0),
          0,
          3.141592653589793
        ),
        V(0, 0, -1),
        0,
        3.141592653589793,
        -2,
        0
      ),
      [
        new StraightEdge(
          new L3(V(0.5, 0, 0), V(0, 0, -1), -4096, 4096),
          V(0.5, 0, 0),
          V(0.5, 0, 1),
          0,
          -1
        ),
        new PCurveEdge(
          new EllipseCurve(
            V(0.5, 0.2, 1),
            V(0.2, 0, 0),
            V(0, -0.2, 0),
            0,
            3.141592653589793
          ),
          V(0.5, 0, 1),
          V(0.7, 0.2, 1),
          1.5707963267948966,
          0,
          undefined,
          V(0.2, 0, 0),
          V(0, 0.2, 0),
          "genseg95"
        ),
        new StraightEdge(
          new L3(V(0.7, 0.2, 0), V3.Z, -4096, 4096),
          V(0.7, 0.2, 1),
          V(0.7, 0.2, 0),
          1,
          0
        ),
        new PCurveEdge(
          new EllipseCurve(
            V(0.5, 0.2, 0),
            V(-0.2, 0, 0),
            V(0, -0.2, 0),
            0,
            3.141592653589793
          ),
          V(0.7, 0.2, 0),
          V(0.5, 0, 0),
          3.141592653589793,
          1.5707963267948966,
          undefined,
          V(0, -0.2, 0),
          V(-0.2, 0, 0),
          "genseg93"
        ),
      ],
      []
    ),
  ],
  false
)
