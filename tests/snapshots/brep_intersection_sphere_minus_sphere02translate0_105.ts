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
        1.5707963267948966
      ),
      [
        new PCurveEdge(
          new EllipseCurve(V3.O, V(0, 0, -1), V3.X, 0, 3.141592653589793),
          V3.Z,
          V(0, 0, -1),
          3.141592653589793,
          0,
          undefined,
          V3.X,
          V(-1, 0, 0),
          undefined
        ),
        new PCurveEdge(
          new EllipseCurve(
            V3.O,
            V(0, 0, -1),
            V(-1, 1.2246467991473532e-16, 0),
            0,
            3.141592653589793
          ),
          V(0, 0, -1),
          V3.Z,
          0,
          3.141592653589793,
          undefined,
          V(-1, 1.2246467991473532e-16, 0),
          V(1, -1.2246467991473532e-16, 0),
          undefined
        ),
      ],
      [
        [
          new PCurveEdge(
            new EllipseCurve(
              V(-0.567040442954097, 0.5670404429540966, 0.5670404429540967),
              V(
                -0.07680647559859467,
                -0.15361295119718937,
                0.07680647559859467
              ),
              V(0.1330327180870652, 0, 0.1330327180870652),
              0,
              3.141592653589793
            ),
            V(-0.4902339673555023, 0.720653394151286, 0.49023396735550206),
            V(-0.6438469185526916, 0.41342749175690724, 0.6438469185526914),
            3.141592653589793,
            0,
            undefined,
            V(0.1330327180870652, -1.8812160899121655e-17, 0.1330327180870652),
            V(-0.1330327180870652, 0, -0.1330327180870652),
            "genseg2undefinedundefined"
          ),
          new PCurveEdge(
            new EllipseCurve(
              V(-0.567040442954097, 0.5670404429540966, 0.5670404429540967),
              V(0.07680647559859467, 0.15361295119718937, -0.07680647559859467),
              V(-0.1330327180870652, 0, -0.1330327180870652),
              0,
              3.141592653589793
            ),
            V(-0.6438469185526916, 0.41342749175690724, 0.6438469185526914),
            V(-0.4902339673555023, 0.720653394151286, 0.49023396735550206),
            3.141592653589793,
            0,
            undefined,
            V(-0.1330327180870652, 1.8812160899121655e-17, -0.1330327180870652),
            V(0.1330327180870652, 0, 0.1330327180870652),
            "genseg3undefinedundefined"
          ),
        ],
      ]
    ),
    new RotationFace(
      new EllipsoidSurface(
        V3.O,
        V(-1, 1.2246467991473532e-16, 0),
        V(-1.2246467991473532e-16, -1, 0),
        V3.Z,
        0,
        3.141592653589793,
        -1.5707963267948966,
        1.5707963267948966
      ),
      [
        new PCurveEdge(
          new EllipseCurve(
            V3.O,
            V(0, 0, -1),
            V(-1, 1.2246467991473532e-16, 0),
            0,
            3.141592653589793
          ),
          V3.Z,
          V(0, 0, -1),
          3.141592653589793,
          0,
          undefined,
          V(-1, 1.2246467991473532e-16, 0),
          V(1, -1.2246467991473532e-16, 0),
          undefined
        ),
        new PCurveEdge(
          new EllipseCurve(V3.O, V(0, 0, -1), V3.X, 0, 3.141592653589793),
          V(0, 0, -1),
          V3.Z,
          0,
          3.141592653589793,
          undefined,
          V3.X,
          V(-1, 0, 0),
          undefined
        ),
      ],
      []
    ),
    new RotationFace(
      new EllipsoidSurface(
        V(-0.6062177826491073, 0.606217782649107, 0.6062177826491071),
        V(-0.15773502691896257, -0.11547005383792514, -0.0422649730810374),
        V(0.11547005383792514, -0.11547005383792515, -0.11547005383792515),
        V(-0.04226497308103742, 0.11547005383792515, -0.15773502691896255),
        0,
        3.141592653589793,
        -1.5707963267948966,
        1.5707963267948966
      ),
      [
        new PCurveEdge(
          new EllipseCurve(
            V(-0.567040442954097, 0.5670404429540966, 0.5670404429540967),
            V(-0.07680647559859467, -0.15361295119718937, 0.07680647559859467),
            V(0.1330327180870652, 0, 0.1330327180870652),
            0,
            3.141592653589793
          ),
          V(-0.6438469185526916, 0.41342749175690724, 0.6438469185526914),
          V(-0.4902339673555023, 0.720653394151286, 0.49023396735550206),
          0,
          3.141592653589793,
          undefined,
          V(0.1330327180870652, 0, 0.1330327180870652),
          V(-0.1330327180870652, 1.8812160899121655e-17, -0.1330327180870652),
          "genseg2undefinedundefined"
        ),
        new PCurveEdge(
          new EllipseCurve(
            V(-0.567040442954097, 0.5670404429540966, 0.5670404429540967),
            V(0.07680647559859467, 0.15361295119718937, -0.07680647559859467),
            V(-0.1330327180870652, 0, -0.1330327180870652),
            0,
            3.141592653589793
          ),
          V(-0.4902339673555023, 0.720653394151286, 0.49023396735550206),
          V(-0.6438469185526916, 0.41342749175690724, 0.6438469185526914),
          0,
          3.141592653589793,
          undefined,
          V(-0.1330327180870652, 0, -0.1330327180870652),
          V(0.1330327180870652, -1.8812160899121655e-17, 0.1330327180870652),
          "genseg3undefinedundefined"
        ),
      ],
      []
    ),
  ],
  false
)
