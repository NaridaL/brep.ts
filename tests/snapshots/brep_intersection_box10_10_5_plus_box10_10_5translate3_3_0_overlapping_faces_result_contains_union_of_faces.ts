new BRep(
  [
    new PlaneFace(
      new PlaneSurface(
        new P3(V(0, -1, 0), -10),
        V3.X,
        V3.Z,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(3, 10, 0), V3.Z, -4096, 4096),
          V(3, 10, 0),
          V(3, 10, 5),
          0,
          5
        ),
        new StraightEdge(
          new L3(V(0, 10, 5), V3.X, -4096, 4096),
          V(3, 10, 5),
          V(0, 10, 5),
          3,
          0
        ),
        new StraightEdge(
          new L3(V(0, 10, 0), V3.Z, -4096, 4096),
          V(0, 10, 5),
          V(0, 10, 0),
          5,
          0
        ),
        new StraightEdge(
          new L3(V(0, 10, 0), V3.X, -4096, 4096),
          V(0, 10, 0),
          V(3, 10, 0),
          0,
          3
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(-1, 0, 0), -10),
        V(0, -1, 0),
        V3.Z,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(10, 3, 0), V(0, 0, -1), -4096, 4096),
          V(10, 3, 5),
          V(10, 3, 0),
          -5,
          0
        ),
        new StraightEdge(
          new L3(V(10, 10, 0), V(0, -1, 0), -4096, 4096),
          V(10, 3, 0),
          V(10, 0, 0),
          7,
          10
        ),
        new StraightEdge(
          new L3(V(10, 0, 0), V3.Z, -4096, 4096),
          V(10, 0, 0),
          V(10, 0, 5),
          0,
          5
        ),
        new StraightEdge(
          new L3(V(10, 10, 5), V(0, -1, 0), -4096, 4096),
          V(10, 0, 5),
          V(10, 3, 5),
          10,
          7
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Z, 0),
        V3.Y,
        V(-1, 0, 0),
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(3, 0, 0), V3.Y, -4096, 4096),
          V(3, 3, 0),
          V(3, 10, 0),
          3,
          10
        ),
        new StraightEdge(
          new L3(V(0, 10, 0), V3.X, -4096, 4096),
          V(3, 10, 0),
          V(0, 10, 0),
          3,
          0
        ),
        new StraightEdge(
          new L3(V3.O, V3.Y, -4096, 4096),
          V(0, 10, 0),
          V3.O,
          10,
          0
        ),
        new StraightEdge(
          new L3(V(10, 0, 0), V(-1, 0, 0), -4096, 4096),
          V3.O,
          V(10, 0, 0),
          10,
          0
        ),
        new StraightEdge(
          new L3(V(10, 10, 0), V(0, -1, 0), -4096, 4096),
          V(10, 0, 0),
          V(10, 3, 0),
          10,
          7
        ),
        new StraightEdge(
          new L3(V(0, 3, 0), V(-1, 0, 0), -4096, 4096),
          V(10, 3, 0),
          V(3, 3, 0),
          -10,
          -3
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Z, 0),
        V3.Y,
        V(-1, 0, 0),
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(3, 0, 0), V3.Y, -4096, 4096),
          V(3, 10, 0),
          V(3, 3, 0),
          10,
          3
        ),
        new StraightEdge(
          new L3(V(0, 3, 0), V(-1, 0, 0), -4096, 4096),
          V(3, 3, 0),
          V(10, 3, 0),
          -3,
          -10
        ),
        new StraightEdge(
          new L3(V(10, 0, 0), V3.Y, -4096, 4096),
          V(10, 3, 0),
          V(10, 10, 0),
          3,
          10
        ),
        new StraightEdge(
          new L3(V(0, 10, 0), V(-1, 0, 0), -4096, 4096),
          V(10, 10, 0),
          V(3, 10, 0),
          -10,
          -3
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(0, 0, -1), -5),
        V3.Y,
        V3.X,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(3, 0, 5), V(0, -1, 0), -4096, 4096),
          V(3, 10, 5),
          V(3, 3, 5),
          -10,
          -3
        ),
        new StraightEdge(
          new L3(V(0, 3, 5), V3.X, -4096, 4096),
          V(3, 3, 5),
          V(10, 3, 5),
          3,
          10
        ),
        new StraightEdge(
          new L3(V(10, 10, 5), V(0, -1, 0), -4096, 4096),
          V(10, 3, 5),
          V(10, 0, 5),
          7,
          10
        ),
        new StraightEdge(
          new L3(V(10, 0, 5), V(-1, 0, 0), -4096, 4096),
          V(10, 0, 5),
          V(0, 0, 5),
          0,
          10
        ),
        new StraightEdge(
          new L3(V(0, 0, 5), V3.Y, -4096, 4096),
          V(0, 0, 5),
          V(0, 10, 5),
          0,
          10
        ),
        new StraightEdge(
          new L3(V(0, 10, 5), V3.X, -4096, 4096),
          V(0, 10, 5),
          V(3, 10, 5),
          0,
          3
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(0, 0, -1), -5),
        V3.Y,
        V3.X,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(3, 0, 5), V(0, -1, 0), -4096, 4096),
          V(3, 3, 5),
          V(3, 10, 5),
          -3,
          -10
        ),
        new StraightEdge(
          new L3(V(0, 10, 5), V3.X, -4096, 4096),
          V(3, 10, 5),
          V(10, 10, 5),
          3,
          10
        ),
        new StraightEdge(
          new L3(V(10, 0, 5), V(0, -1, 0), -4096, 4096),
          V(10, 10, 5),
          V(10, 3, 5),
          -10,
          -3
        ),
        new StraightEdge(
          new L3(V(0, 3, 5), V3.X, -4096, 4096),
          V(10, 3, 5),
          V(3, 3, 5),
          10,
          3
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(new P3(V3.X, 0), V3.Y, V3.Z, -100, 100, -100, 100),
      [
        new StraightEdge(
          new L3(V(0, 10, 0), V3.Z, -4096, 4096),
          V(0, 10, 0),
          V(0, 10, 5),
          0,
          5
        ),
        new StraightEdge(
          new L3(V(0, 0, 5), V3.Y, -4096, 4096),
          V(0, 10, 5),
          V(0, 0, 5),
          10,
          0
        ),
        new StraightEdge(
          new L3(V3.O, V3.Z, -4096, 4096),
          V(0, 0, 5),
          V3.O,
          5,
          0
        ),
        new StraightEdge(
          new L3(V3.O, V3.Y, -4096, 4096),
          V3.O,
          V(0, 10, 0),
          0,
          10
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Y, 0),
        V(-1, 0, 0),
        V3.Z,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V3.O, V3.Z, -4096, 4096),
          V3.O,
          V(0, 0, 5),
          0,
          5
        ),
        new StraightEdge(
          new L3(V(10, 0, 5), V(-1, 0, 0), -4096, 4096),
          V(0, 0, 5),
          V(10, 0, 5),
          10,
          0
        ),
        new StraightEdge(
          new L3(V(10, 0, 0), V3.Z, -4096, 4096),
          V(10, 0, 5),
          V(10, 0, 0),
          5,
          0
        ),
        new StraightEdge(
          new L3(V(10, 0, 0), V(-1, 0, 0), -4096, 4096),
          V(10, 0, 0),
          V3.O,
          0,
          10
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(new P3(V3.X, 3), V3.Y, V3.Z, -100, 100, -100, 100),
      [
        new StraightEdge(
          new L3(V(3, 10, 0), V3.Z, -4096, 4096),
          V(3, 10, 5),
          V(3, 10, 0),
          5,
          0
        ),
        new StraightEdge(
          new L3(V(3, 3, 0), V3.Y, -4096, 4096),
          V(3, 10, 0),
          V(3, 13, 0),
          7,
          10
        ),
        new StraightEdge(
          new L3(V(3, 13, 0), V3.Z, -4096, 4096),
          V(3, 13, 0),
          V(3, 13, 5),
          0,
          5
        ),
        new StraightEdge(
          new L3(V(3, 3, 5), V3.Y, -4096, 4096),
          V(3, 13, 5),
          V(3, 10, 5),
          10,
          7
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Y, 3),
        V(-1, 0, 0),
        V3.Z,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(10, 3, 0), V(0, 0, -1), -4096, 4096),
          V(10, 3, 0),
          V(10, 3, 5),
          0,
          -5
        ),
        new StraightEdge(
          new L3(V(13, 3, 5), V(-1, 0, 0), -4096, 4096),
          V(10, 3, 5),
          V(13, 3, 5),
          3,
          0
        ),
        new StraightEdge(
          new L3(V(13, 3, 0), V3.Z, -4096, 4096),
          V(13, 3, 5),
          V(13, 3, 0),
          5,
          0
        ),
        new StraightEdge(
          new L3(V(13, 3, 0), V(-1, 0, 0), -4096, 4096),
          V(13, 3, 0),
          V(10, 3, 0),
          0,
          3
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Z, 0),
        V3.Y,
        V(-1, 0, 0),
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(0, 10, 0), V(-1, 0, 0), -4096, 4096),
          V(3, 10, 0),
          V(10, 10, 0),
          -3,
          -10
        ),
        new StraightEdge(
          new L3(V(10, 0, 0), V3.Y, -4096, 4096),
          V(10, 10, 0),
          V(10, 3, 0),
          10,
          3
        ),
        new StraightEdge(
          new L3(V(13, 3, 0), V(-1, 0, 0), -4096, 4096),
          V(10, 3, 0),
          V(13, 3, 0),
          3,
          0
        ),
        new StraightEdge(
          new L3(V(13, 13, 0), V(0, -1, 0), -4096, 4096),
          V(13, 3, 0),
          V(13, 13, 0),
          10,
          0
        ),
        new StraightEdge(
          new L3(V(3, 13, 0), V3.X, -4096, 4096),
          V(13, 13, 0),
          V(3, 13, 0),
          10,
          0
        ),
        new StraightEdge(
          new L3(V(3, 3, 0), V3.Y, -4096, 4096),
          V(3, 13, 0),
          V(3, 10, 0),
          10,
          7
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(0, 0, -1), -5),
        V3.Y,
        V3.X,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(0, 10, 5), V3.X, -4096, 4096),
          V(10, 10, 5),
          V(3, 10, 5),
          10,
          3
        ),
        new StraightEdge(
          new L3(V(3, 3, 5), V3.Y, -4096, 4096),
          V(3, 10, 5),
          V(3, 13, 5),
          7,
          10
        ),
        new StraightEdge(
          new L3(V(3, 13, 5), V3.X, -4096, 4096),
          V(3, 13, 5),
          V(13, 13, 5),
          0,
          10
        ),
        new StraightEdge(
          new L3(V(13, 13, 5), V(0, -1, 0), -4096, 4096),
          V(13, 13, 5),
          V(13, 3, 5),
          0,
          10
        ),
        new StraightEdge(
          new L3(V(13, 3, 5), V(-1, 0, 0), -4096, 4096),
          V(13, 3, 5),
          V(10, 3, 5),
          0,
          3
        ),
        new StraightEdge(
          new L3(V(10, 0, 5), V(0, -1, 0), -4096, 4096),
          V(10, 3, 5),
          V(10, 10, 5),
          -3,
          -10
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(0, -1, 0), -13),
        V3.X,
        V3.Z,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(13, 13, 0), V3.Z, -4096, 4096),
          V(13, 13, 0),
          V(13, 13, 5),
          0,
          5
        ),
        new StraightEdge(
          new L3(V(3, 13, 5), V3.X, -4096, 4096),
          V(13, 13, 5),
          V(3, 13, 5),
          10,
          0
        ),
        new StraightEdge(
          new L3(V(3, 13, 0), V3.Z, -4096, 4096),
          V(3, 13, 5),
          V(3, 13, 0),
          5,
          0
        ),
        new StraightEdge(
          new L3(V(3, 13, 0), V3.X, -4096, 4096),
          V(3, 13, 0),
          V(13, 13, 0),
          0,
          10
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(-1, 0, 0), -13),
        V(0, -1, 0),
        V3.Z,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(13, 3, 0), V3.Z, -4096, 4096),
          V(13, 3, 0),
          V(13, 3, 5),
          0,
          5
        ),
        new StraightEdge(
          new L3(V(13, 13, 5), V(0, -1, 0), -4096, 4096),
          V(13, 3, 5),
          V(13, 13, 5),
          10,
          0
        ),
        new StraightEdge(
          new L3(V(13, 13, 0), V3.Z, -4096, 4096),
          V(13, 13, 5),
          V(13, 13, 0),
          5,
          0
        ),
        new StraightEdge(
          new L3(V(13, 13, 0), V(0, -1, 0), -4096, 4096),
          V(13, 13, 0),
          V(13, 3, 0),
          0,
          10
        ),
      ],
      []
    ),
  ],
  true
)
