new BRep(
  [
    new RotationFace(
      new EllipsoidSurface(
        V3.O,
        V3.X,
        V3.Y,
        V3.Z,
        0,
        3.141592653589793,
        -1.5707963267948966,
        1.5707963267948966,
      ),
      [
        new PCurveEdge(
          PICurve.forParametricStartEnd(
            new CylinderSurface(
              new EllipseCurve(
                V(0.5, 0, -2),
                V(0.5, 0, 0),
                V(0, 0.1, 0),
                0,
                3.141592653589793,
              ),
              V(0, 0, -4),
              0,
              3.141592653589793,
              -1,
              0,
            ),
            new EllipsoidSurface(
              V3.O,
              V3.X,
              V3.Y,
              V3.Z,
              0,
              3.141592653589793,
              -1.5707963267948966,
              1.5707963267948966,
            ),
            V(1.2731963988416202e-20, -0.5, 0),
            V(3.141592653589793, -0.75, 0),
            0.05,
            V(0.005120444534810256, -0.0008960753763405105, 0),
            0,
            64,
          ),
          V3.Z,
          V(1, 0, -6.123233995736766e-17),
          64,
          0,
          undefined,
          V(
            3.061616997868383e-18,
            0.005000000000000001,
            -6.123233995736767e-20,
          ),
          V(
            3.2596657710943364e-23,
            -0.0005120444534810256,
            -0.003584301505362042,
          ),
          "genseg839",
        ),
        new PCurveEdge(
          PICurve.forParametricStartEnd(
            new CylinderSurface(
              new EllipseCurve(
                V(0.5, 0, -2),
                V(0.5, 0, 0),
                V(0, 0.1, 0),
                0,
                3.141592653589793,
              ),
              V(0, 0, -4),
              0,
              3.141592653589793,
              -1,
              0,
            ),
            new EllipsoidSurface(
              V3.O,
              V3.X,
              V3.Y,
              V3.Z,
              0,
              3.141592653589793,
              -1.5707963267948966,
              1.5707963267948966,
            ),
            V(3.141592653589793, -0.25, 0),
            V(9.899697571034635e-21, -0.5, 0),
            0.05,
            V(-0.05, -1.5308084989341917e-20, 0),
            0,
            64,
          ),
          V(1, 0, -6.123233995736766e-17),
          V(0, 0, -1),
          64,
          0,
          undefined,
          V(
            -3.800184612919155e-23,
            0.0007677375163536417,
            -0.0053741300224980915,
          ),
          V(
            -3.061616997868383e-18,
            -0.005000000000000001,
            -6.123233995736767e-20,
          ),
          "genseg840",
        ),
        new PCurveEdge(
          new EllipseCurve(
            V3.O,
            V(0, 0, -1),
            V(-1, 1.2246467991473532e-16, 0),
            0,
            3.141592653589793,
          ),
          V(0, 0, -1),
          V3.Z,
          0,
          3.141592653589793,
          undefined,
          V(-1, 1.2246467991473532e-16, 0),
          V(1, -1.2246467991473532e-16, 0),
          "undefinedundefined",
        ),
      ],
      [],
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(-1.2246467991473532e-16, -1, 0), 0),
        V(1, -1.2246467991473532e-16, 0),
        V3.Z,
        -100,
        100,
        -100,
        100,
      ),
      [
        new StraightEdge(
          new L3(V(0, 1.2246467991473533e-17, -2), V(0, 0, -1), -4096, 4096),
          V(0, 0, -1),
          V3.Z,
          -1,
          -3,
        ),
        new PCurveEdge(
          new EllipseCurve(
            V3.O,
            V(0, 0, -1),
            V(-1, 1.2246467991473532e-16, 0),
            0,
            3.141592653589793,
          ),
          V3.Z,
          V(0, 0, -1),
          3.141592653589793,
          0,
          undefined,
          V(-1, 1.2246467991473532e-16, 0),
          V(1, -1.2246467991473532e-16, 0),
          "undefinedundefined",
        ),
      ],
      [],
    ),
    new RotationFace(
      new CylinderSurface(
        new EllipseCurve(
          V(0.5, 0, -2),
          V(0.5, 0, 0),
          V(0, 0.1, 0),
          0,
          3.141592653589793,
        ),
        V(0, 0, -4),
        0,
        3.141592653589793,
        -1,
        0,
      ),
      [
        new PCurveEdge(
          PICurve.forParametricStartEnd(
            new CylinderSurface(
              new EllipseCurve(
                V(0.5, 0, -2),
                V(0.5, 0, 0),
                V(0, 0.1, 0),
                0,
                3.141592653589793,
              ),
              V(0, 0, -4),
              0,
              3.141592653589793,
              -1,
              0,
            ),
            new EllipsoidSurface(
              V3.O,
              V3.X,
              V3.Y,
              V3.Z,
              0,
              3.141592653589793,
              -1.5707963267948966,
              1.5707963267948966,
            ),
            V(1.2731963988416202e-20, -0.5, 0),
            V(3.141592653589793, -0.75, 0),
            0.05,
            V(0.005120444534810256, -0.0008960753763405105, 0),
            0,
            64,
          ),
          V(1, 0, -6.123233995736766e-17),
          V3.Z,
          0,
          64,
          undefined,
          V(
            -3.2596657710943364e-23,
            0.0005120444534810256,
            0.003584301505362042,
          ),
          V(
            -3.061616997868383e-18,
            -0.005000000000000001,
            6.123233995736767e-20,
          ),
          "genseg839",
        ),
        new StraightEdge(
          new L3(V(0, 1.2246467991473533e-17, -2), V(0, 0, -1), -4096, 4096),
          V3.Z,
          V(0, 0, -1),
          -3,
          -1,
        ),
        new PCurveEdge(
          PICurve.forParametricStartEnd(
            new CylinderSurface(
              new EllipseCurve(
                V(0.5, 0, -2),
                V(0.5, 0, 0),
                V(0, 0.1, 0),
                0,
                3.141592653589793,
              ),
              V(0, 0, -4),
              0,
              3.141592653589793,
              -1,
              0,
            ),
            new EllipsoidSurface(
              V3.O,
              V3.X,
              V3.Y,
              V3.Z,
              0,
              3.141592653589793,
              -1.5707963267948966,
              1.5707963267948966,
            ),
            V(3.141592653589793, -0.25, 0),
            V(9.899697571034635e-21, -0.5, 0),
            0.05,
            V(-0.05, -1.5308084989341917e-20, 0),
            0,
            64,
          ),
          V(0, 0, -1),
          V(1, 0, -6.123233995736766e-17),
          0,
          64,
          undefined,
          V(3.061616997868383e-18, 0.005000000000000001, 6.123233995736767e-20),
          V(
            3.800184612919155e-23,
            -0.0007677375163536417,
            0.0053741300224980915,
          ),
          "genseg840",
        ),
      ],
      [],
    ),
  ],
  false,
)
