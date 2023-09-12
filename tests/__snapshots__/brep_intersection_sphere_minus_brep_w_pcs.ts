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
          new EllipseCurve(V3.O, V(0, 0, -1), V3.X, 0, 3.141592653589793),
          V3.Z,
          V(0, 0, -1),
          3.141592653589793,
          0,
          undefined,
          V3.X,
          V(-1, 0, 0),
          undefined,
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
      [
        [
          new PCurveEdge(
            PICurve.forParametricStartEnd(
              new ProjectedCurveSurface(
                new BezierCurve(
                  V(0.30000000000000004, 0.1, 0.4),
                  V(
                    0.30000000000000004,
                    0.10962695751594148,
                    0.5100366279236725,
                  ),
                  V(
                    0.2104569499661587,
                    0.11743114854953163,
                    0.5992389396183492,
                  ),
                  V(
                    0.10000000000000002,
                    0.11743114854953163,
                    0.5992389396183492,
                  ),
                  0,
                  1,
                ),
                V(0, 1.992389396183491, -0.17431148549531628),
                0,
                1,
                0,
                1,
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
              V(1, 0.3617394511707665, 0),
              V(0, 0.39895466399129986, 0),
              0.05,
              V(-0.04998895788866387, -0.0010507564919570807, 0),
              0,
              20,
            ),
            V(0.30000000000000004, 0.8948730420942135, 0.33045761987439176),
            V(0.1, 0.8381569952434027, 0.5361835985225124),
            20,
            0,
            undefined,
            V(0, -0.006327513258486222, 0.017134787331172165),
            V(
              -0.016564883461105677,
              0.002093516092546251,
              -0.00018315892500688608,
            ),
            "genseg798",
          ),
          new PCurveEdge(
            new EllipseCurve(
              V(0.3535967267557321, 0.030817985353536557, 0.352251184456656),
              V(0.03780987199589105, -0.8643439107956783, 0.03766599401783423),
              V(-0.6111789986634738, 0, 0.6135136031482741),
              1.535133980799257,
              3.141592653589793,
            ),
            V(0.1, 0.8381569952434027, 0.5361835985225124),
            V(0.30000000000000004, 0.8948730420942135, 0.33045761987439176),
            2.7763814662109167,
            3.1157389336980312,
            undefined,
            V(0.5573670738448863, 0.30869745583599784, -0.5865036533639291),
            V(0.6099973325422855, 0.02234401598712118, -0.6142822713509875),
            "genseg799",
          ),
        ],
      ],
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
        1.5707963267948966,
      ),
      [
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
        new PCurveEdge(
          new EllipseCurve(V3.O, V(0, 0, -1), V3.X, 0, 3.141592653589793),
          V(0, 0, -1),
          V3.Z,
          0,
          3.141592653589793,
          undefined,
          V3.X,
          V(-1, 0, 0),
          undefined,
        ),
      ],
      [],
    ),
    new RotationFace(
      new ProjectedCurveSurface(
        new BezierCurve(
          V(0.30000000000000004, 0.1, 0.4),
          V(0.30000000000000004, 0.10962695751594148, 0.5100366279236725),
          V(0.2104569499661587, 0.11743114854953163, 0.5992389396183492),
          V(0.10000000000000002, 0.11743114854953163, 0.5992389396183492),
          0,
          1,
        ),
        V(0, 1.992389396183491, -0.17431148549531628),
        0,
        1,
        0,
        1,
      ),
      [
        new PCurveEdge(
          PICurve.forParametricStartEnd(
            new ProjectedCurveSurface(
              new BezierCurve(
                V(0.30000000000000004, 0.1, 0.4),
                V(0.30000000000000004, 0.10962695751594148, 0.5100366279236725),
                V(0.2104569499661587, 0.11743114854953163, 0.5992389396183492),
                V(0.10000000000000002, 0.11743114854953163, 0.5992389396183492),
                0,
                1,
              ),
              V(0, 1.992389396183491, -0.17431148549531628),
              0,
              1,
              0,
              1,
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
            V(1, 0.3617394511707665, 0),
            V(0, 0.39895466399129986, 0),
            0.05,
            V(-0.04998895788866387, -0.0010507564919570807, 0),
            0,
            20,
          ),
          V(0.1, 0.8381569952434027, 0.5361835985225124),
          V(0.30000000000000004, 0.8948730420942135, 0.33045761987439176),
          0,
          20,
          undefined,
          V(
            0.016564883461105677,
            -0.002093516092546251,
            0.00018315892500688608,
          ),
          V(0, 0.006327513258486222, -0.017134787331172165),
          "genseg798",
        ),
        new StraightEdge(
          new L3(
            V(0.30000000000000004, 0.1, 0.4),
            V(0, 0.9961946980917455, -0.08715574274765814),
            0,
            2,
          ),
          V(0.30000000000000004, 0.8948730420942135, 0.33045761987439176),
          V(0.30000000000000004, 0.1, 0.4),
          0.7979093279825997,
          0,
        ),
        new PCurveEdge(
          new BezierCurve(
            V(0.30000000000000004, 0.1, 0.4),
            V(0.30000000000000004, 0.10962695751594148, 0.5100366279236725),
            V(0.2104569499661587, 0.11743114854953163, 0.5992389396183492),
            V(0.10000000000000002, 0.11743114854953163, 0.5992389396183492),
            0,
            1,
          ),
          V(0.30000000000000004, 0.1, 0.4),
          V(0.10000000000000002, 0.11743114854953163, 0.5992389396183492),
          0,
          1,
          undefined,
          V(0, 0.028880872547824416, 0.33010988377101746),
          V(-0.33137084989847604, 0, 0),
          "undefinedundefinedundefinedundefined",
        ),
        new StraightEdge(
          new L3(
            V(0.1, 0.11743114854953163, 0.5992389396183492),
            V(0, 0.9961946980917455, -0.08715574274765814),
            0,
            2,
          ),
          V(0.1, 0.11743114854953163, 0.5992389396183492),
          V(0.1, 0.8381569952434027, 0.5361835985225124),
          0,
          0.7234789023415331,
        ),
      ],
      [],
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(
          V(0.7071067811865476, 0.06162841671621934, 0.7044160264027587),
          0.5000612865886898,
        ),
        V(0.08682659386424758, -0.9962234400966147, 0),
        V(0.701755757082144, 0.06116204423593943, -0.7097873355780224),
        -100,
        100,
        -100,
        100,
      ),
      [
        new PCurveEdge(
          new EllipseCurve(
            V(0.3535967267557321, 0.030817985353536557, 0.352251184456656),
            V(0.03780987199589105, -0.8643439107956783, 0.03766599401783423),
            V(-0.6111789986634738, 0, 0.6135136031482741),
            1.535133980799257,
            3.141592653589793,
          ),
          V(0.30000000000000004, 0.8948730420942135, 0.33045761987439176),
          V(0.1, 0.8381569952434027, 0.5361835985225124),
          3.1157389336980312,
          2.7763814662109167,
          undefined,
          V(-0.6099973325422855, -0.02234401598712118, 0.6142822713509875),
          V(-0.5573670738448863, -0.30869745583599784, 0.5865036533639291),
          "genseg799",
        ),
        new StraightEdge(
          new L3(
            V(0.1, 0.11743114854953163, 0.5992389396183492),
            V(0, 0.9961946980917455, -0.08715574274765814),
            0,
            2,
          ),
          V(0.1, 0.8381569952434027, 0.5361835985225124),
          V(0.1, 0.11743114854953163, 0.5992389396183492),
          0.7234789023415331,
          0,
        ),
        new StraightEdge(
          new L3(
            V(0.1, 0.11743114854953163, 0.5992389396183492),
            V(0.7071067811865476, -0.06162841671621934, -0.7044160264027587),
            0,
            0.282842712474619,
          ),
          V(0.1, 0.11743114854953163, 0.5992389396183492),
          V(0.30000000000000004, 0.1, 0.4),
          0,
          0.282842712474619,
        ),
        new StraightEdge(
          new L3(
            V(0.30000000000000004, 0.1, 0.4),
            V(0, 0.9961946980917455, -0.08715574274765814),
            0,
            2,
          ),
          V(0.30000000000000004, 0.1, 0.4),
          V(0.30000000000000004, 0.8948730420942135, 0.33045761987439176),
          0,
          0.7979093279825997,
        ),
      ],
      [],
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(
          V(0, 0.9961946980917455, -0.08715574274765814),
          0.0647571727101113,
        ),
        V3.X,
        V(0, -0.08715574274765814, -0.9961946980917455),
        -100,
        100,
        -100,
        100,
      ),
      [
        new StraightEdge(
          new L3(
            V(0.1, 0.11743114854953163, 0.5992389396183492),
            V(0.7071067811865476, -0.06162841671621934, -0.7044160264027587),
            0,
            0.282842712474619,
          ),
          V(0.30000000000000004, 0.1, 0.4),
          V(0.1, 0.11743114854953163, 0.5992389396183492),
          0.282842712474619,
          0,
        ),
        new PCurveEdge(
          new BezierCurve(
            V(0.30000000000000004, 0.1, 0.4),
            V(0.30000000000000004, 0.10962695751594148, 0.5100366279236725),
            V(0.2104569499661587, 0.11743114854953163, 0.5992389396183492),
            V(0.10000000000000002, 0.11743114854953163, 0.5992389396183492),
            0,
            1,
          ),
          V(0.10000000000000002, 0.11743114854953163, 0.5992389396183492),
          V(0.30000000000000004, 0.1, 0.4),
          1,
          0,
          undefined,
          V(0.33137084989847604, 0, 0),
          V(0, -0.028880872547824416, -0.33010988377101746),
          "undefinedundefinedundefinedundefined",
        ),
      ],
      [],
    ),
  ],
  false,
)