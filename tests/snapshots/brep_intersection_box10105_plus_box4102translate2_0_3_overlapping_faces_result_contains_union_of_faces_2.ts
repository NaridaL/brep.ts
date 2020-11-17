new BRep(
  [
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Y, 10),
        V(-1, 0, 0),
        V3.Z,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(0, 10, 3), V(-1, 0, 0), -4096, 4096),
          V(2, 10, 3),
          V(6, 10, 3),
          -2,
          -6
        ),
        new StraightEdge(
          new L3(V(6, 10, 0), V(0, 0, -1), -4096, 4096),
          V(6, 10, 3),
          V(6, 10, 5),
          -3,
          -5
        ),
        new StraightEdge(
          new L3(V(0, 10, 5), V3.X, -4096, 4096),
          V(6, 10, 5),
          V(10, 10, 5),
          6,
          10
        ),
        new StraightEdge(
          new L3(V(10, 10, 0), V3.Z, -4096, 4096),
          V(10, 10, 5),
          V(10, 10, 0),
          5,
          0
        ),
        new StraightEdge(
          new L3(V(0, 10, 0), V3.X, -4096, 4096),
          V(10, 10, 0),
          V(0, 10, 0),
          10,
          0
        ),
        new StraightEdge(
          new L3(V(0, 10, 0), V3.Z, -4096, 4096),
          V(0, 10, 0),
          V(0, 10, 5),
          0,
          5
        ),
        new StraightEdge(
          new L3(V(0, 10, 5), V3.X, -4096, 4096),
          V(0, 10, 5),
          V(2, 10, 5),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(2, 10, 0), V3.Z, -4096, 4096),
          V(2, 10, 5),
          V(2, 10, 3),
          5,
          3
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Y, 10),
        V(-1, 0, 0),
        V3.Z,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(0, 10, 5), V3.X, -4096, 4096),
          V(2, 10, 5),
          V(6, 10, 5),
          2,
          6
        ),
        new StraightEdge(
          new L3(V(6, 10, 0), V(0, 0, -1), -4096, 4096),
          V(6, 10, 5),
          V(6, 10, 3),
          -5,
          -3
        ),
        new StraightEdge(
          new L3(V(0, 10, 3), V(-1, 0, 0), -4096, 4096),
          V(6, 10, 3),
          V(2, 10, 3),
          -6,
          -2
        ),
        new StraightEdge(
          new L3(V(2, 10, 0), V3.Z, -4096, 4096),
          V(2, 10, 3),
          V(2, 10, 5),
          3,
          5
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
          new L3(V(10, 0, 5), V(-1, 0, 0), -4096, 4096),
          V(2, 0, 5),
          V(0, 0, 5),
          8,
          10
        ),
        new StraightEdge(
          new L3(V3.O, V3.Z, -4096, 4096),
          V(0, 0, 5),
          V3.O,
          5,
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
          new L3(V(10, 0, 0), V3.Z, -4096, 4096),
          V(10, 0, 0),
          V(10, 0, 5),
          0,
          5
        ),
        new StraightEdge(
          new L3(V(10, 0, 5), V(-1, 0, 0), -4096, 4096),
          V(10, 0, 5),
          V(6, 0, 5),
          0,
          4
        ),
        new StraightEdge(
          new L3(V(6, 0, 0), V3.Z, -4096, 4096),
          V(6, 0, 5),
          V(6, 0, 3),
          5,
          3
        ),
        new StraightEdge(
          new L3(V(0, 0, 3), V3.X, -4096, 4096),
          V(6, 0, 3),
          V(2, 0, 3),
          6,
          2
        ),
        new StraightEdge(
          new L3(V(2, 0, 0), V(0, 0, -1), -4096, 4096),
          V(2, 0, 3),
          V(2, 0, 5),
          -3,
          -5
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
          new L3(V(0, 0, 3), V3.X, -4096, 4096),
          V(2, 0, 3),
          V(6, 0, 3),
          2,
          6
        ),
        new StraightEdge(
          new L3(V(6, 0, 0), V3.Z, -4096, 4096),
          V(6, 0, 3),
          V(6, 0, 5),
          3,
          5
        ),
        new StraightEdge(
          new L3(V(0, 0, 5), V(-1, 0, 0), -4096, 4096),
          V(6, 0, 5),
          V(2, 0, 5),
          -6,
          -2
        ),
        new StraightEdge(
          new L3(V(2, 0, 0), V(0, 0, -1), -4096, 4096),
          V(2, 0, 5),
          V(2, 0, 3),
          -5,
          -3
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Z, 5),
        V3.Y,
        V(-1, 0, 0),
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(0, 10, 5), V3.X, -4096, 4096),
          V(2, 10, 5),
          V(0, 10, 5),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(0, 0, 5), V3.Y, -4096, 4096),
          V(0, 10, 5),
          V(0, 0, 5),
          10,
          0
        ),
        new StraightEdge(
          new L3(V(10, 0, 5), V(-1, 0, 0), -4096, 4096),
          V(0, 0, 5),
          V(2, 0, 5),
          10,
          8
        ),
        new StraightEdge(
          new L3(V(2, 0, 5), V(0, -1, 0), -4096, 4096),
          V(2, 0, 5),
          V(2, 10, 5),
          0,
          -10
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Z, 5),
        V3.Y,
        V(-1, 0, 0),
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(0, 0, 5), V(-1, 0, 0), -4096, 4096),
          V(2, 0, 5),
          V(6, 0, 5),
          -2,
          -6
        ),
        new StraightEdge(
          new L3(V(6, 0, 5), V3.Y, -4096, 4096),
          V(6, 0, 5),
          V(6, 10, 5),
          0,
          10
        ),
        new StraightEdge(
          new L3(V(0, 10, 5), V3.X, -4096, 4096),
          V(6, 10, 5),
          V(2, 10, 5),
          6,
          2
        ),
        new StraightEdge(
          new L3(V(2, 0, 5), V(0, -1, 0), -4096, 4096),
          V(2, 10, 5),
          V(2, 0, 5),
          -10,
          0
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Z, 5),
        V3.Y,
        V(-1, 0, 0),
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(10, 0, 5), V(-1, 0, 0), -4096, 4096),
          V(6, 0, 5),
          V(10, 0, 5),
          4,
          0
        ),
        new StraightEdge(
          new L3(V(10, 10, 5), V(0, -1, 0), -4096, 4096),
          V(10, 0, 5),
          V(10, 10, 5),
          10,
          0
        ),
        new StraightEdge(
          new L3(V(0, 10, 5), V3.X, -4096, 4096),
          V(10, 10, 5),
          V(6, 10, 5),
          10,
          6
        ),
        new StraightEdge(
          new L3(V(6, 0, 5), V3.Y, -4096, 4096),
          V(6, 10, 5),
          V(6, 0, 5),
          10,
          0
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
        new StraightEdge(
          new L3(V3.O, V3.Y, -4096, 4096),
          V(0, 10, 0),
          V3.O,
          10,
          0
        ),
        new StraightEdge(
          new L3(V3.O, V3.Z, -4096, 4096),
          V3.O,
          V(0, 0, 5),
          0,
          5
        ),
        new StraightEdge(
          new L3(V(0, 0, 5), V3.Y, -4096, 4096),
          V(0, 0, 5),
          V(0, 10, 5),
          0,
          10
        ),
        new StraightEdge(
          new L3(V(0, 10, 0), V3.Z, -4096, 4096),
          V(0, 10, 5),
          V(0, 10, 0),
          5,
          0
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(new P3(V3.X, 10), V3.Y, V3.Z, -100, 100, -100, 100),
      [
        new StraightEdge(
          new L3(V(10, 10, 0), V(0, -1, 0), -4096, 4096),
          V(10, 0, 0),
          V(10, 10, 0),
          10,
          0
        ),
        new StraightEdge(
          new L3(V(10, 10, 0), V3.Z, -4096, 4096),
          V(10, 10, 0),
          V(10, 10, 5),
          0,
          5
        ),
        new StraightEdge(
          new L3(V(10, 10, 5), V(0, -1, 0), -4096, 4096),
          V(10, 10, 5),
          V(10, 0, 5),
          0,
          10
        ),
        new StraightEdge(
          new L3(V(10, 0, 0), V3.Z, -4096, 4096),
          V(10, 0, 5),
          V(10, 0, 0),
          5,
          0
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
        new StraightEdge(
          new L3(V3.O, V3.Y, -4096, 4096),
          V3.O,
          V(0, 10, 0),
          0,
          10
        ),
        new StraightEdge(
          new L3(V(0, 10, 0), V3.X, -4096, 4096),
          V(0, 10, 0),
          V(10, 10, 0),
          0,
          10
        ),
        new StraightEdge(
          new L3(V(10, 10, 0), V(0, -1, 0), -4096, 4096),
          V(10, 10, 0),
          V(10, 0, 0),
          0,
          10
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
  ],
  false
)
