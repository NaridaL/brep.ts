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
          V(0, 0, 10),
          0,
          10
        ),
        new StraightEdge(
          new L3(V(0, 0, 10), V3.Y, -4096, 4096),
          V(0, 0, 10),
          V(0, 10, 10),
          0,
          10
        ),
        new StraightEdge(
          new L3(V(0, 10, 0), V3.Z, -4096, 4096),
          V(0, 10, 10),
          V(0, 10, 0),
          10,
          0
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(new P3(V3.X, 10), V3.Y, V3.Z, -100, 100, -100, 100),
      [
        new StraightEdge(
          new L3(V(10, 0, 0), V(0, -1, 0), -4096, 4096),
          V(10, 0, 0),
          V(10, 10, 0),
          0,
          -10
        ),
        new StraightEdge(
          new L3(V(10, 10, 0), V(0, 0, -1), -4096, 4096),
          V(10, 10, 0),
          V(10, 10, 10),
          0,
          -10
        ),
        new StraightEdge(
          new L3(V(10, 0, 10), V(0, -1, 0), -4096, 4096),
          V(10, 10, 10),
          V(10, 0, 10),
          -10,
          0
        ),
        new StraightEdge(
          new L3(V(10, 0, 0), V(0, 0, -1), -4096, 4096),
          V(10, 0, 10),
          V(10, 0, 0),
          -10,
          0
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
          new L3(V3.O, V(-1, 0, 0), -4096, 4096),
          V3.O,
          V(10, 0, 0),
          0,
          -10
        ),
        new StraightEdge(
          new L3(V(10, 0, 0), V(0, 0, -1), -4096, 4096),
          V(10, 0, 0),
          V(10, 0, 10),
          0,
          -10
        ),
        new StraightEdge(
          new L3(V(0, 0, 10), V(-1, 0, 0), -4096, 4096),
          V(10, 0, 10),
          V(0, 0, 10),
          -10,
          0
        ),
        new StraightEdge(
          new L3(V3.O, V3.Z, -4096, 4096),
          V(0, 0, 10),
          V3.O,
          10,
          0
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Z, 10),
        V3.Y,
        V(-1, 0, 0),
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(0, 0, 10), V(-1, 0, 0), -4096, 4096),
          V(0, 0, 10),
          V(10, 0, 10),
          0,
          -10
        ),
        new StraightEdge(
          new L3(V(10, 0, 10), V(0, -1, 0), -4096, 4096),
          V(10, 0, 10),
          V(10, 10, 10),
          0,
          -10
        ),
        new StraightEdge(
          new L3(V(0, 10, 10), V3.X, -4096, 4096),
          V(10, 10, 10),
          V(0, 10, 10),
          10,
          0
        ),
        new StraightEdge(
          new L3(V(0, 0, 10), V3.Y, -4096, 4096),
          V(0, 10, 10),
          V(0, 0, 10),
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
          new L3(V(0, 0, 10), V3.Y, -4096, 4096),
          V(0, 10, 10),
          V(0, 12, 10),
          10,
          12
        ),
        new StraightEdge(
          new L3(V(0, 12, -2), V3.Z, -4096, 4096),
          V(0, 12, 10),
          V(0, 12, -2),
          12,
          0
        ),
        new StraightEdge(
          new L3(V(0, 0, -2), V3.Y, -4096, 4096),
          V(0, 12, -2),
          V(0, 0, -2),
          12,
          0
        ),
        new StraightEdge(
          new L3(V(0, 0, -2), V3.Z, -4096, 4096),
          V(0, 0, -2),
          V3.O,
          0,
          2
        ),
        new StraightEdge(
          new L3(V3.O, V3.Y, -4096, 4096),
          V3.O,
          V(0, 10, 0),
          0,
          10
        ),
        new StraightEdge(
          new L3(V(0, 10, 0), V3.Z, -4096, 4096),
          V(0, 10, 0),
          V(0, 10, 10),
          0,
          10
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(new P3(V3.X, 10), V3.Y, V3.Z, -100, 100, -100, 100),
      [
        new StraightEdge(
          new L3(V(10, 0, 0), V(0, -1, 0), -4096, 4096),
          V(10, 10, 0),
          V(10, 0, 0),
          -10,
          0
        ),
        new StraightEdge(
          new L3(V(10, 0, -2), V3.Z, -4096, 4096),
          V(10, 0, 0),
          V(10, 0, -2),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(10, 12, -2), V(0, -1, 0), -4096, 4096),
          V(10, 0, -2),
          V(10, 12, -2),
          12,
          0
        ),
        new StraightEdge(
          new L3(V(10, 12, -2), V3.Z, -4096, 4096),
          V(10, 12, -2),
          V(10, 12, 10),
          0,
          12
        ),
        new StraightEdge(
          new L3(V(10, 12, 10), V(0, -1, 0), -4096, 4096),
          V(10, 12, 10),
          V(10, 10, 10),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(10, 10, 0), V(0, 0, -1), -4096, 4096),
          V(10, 10, 10),
          V(10, 10, 0),
          -10,
          0
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
          new L3(V(0, 0, -2), V3.Z, -4096, 4096),
          V3.O,
          V(0, 0, -2),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(10, 0, -2), V(-1, 0, 0), -4096, 4096),
          V(0, 0, -2),
          V(10, 0, -2),
          10,
          0
        ),
        new StraightEdge(
          new L3(V(10, 0, -2), V3.Z, -4096, 4096),
          V(10, 0, -2),
          V(10, 0, 0),
          0,
          2
        ),
        new StraightEdge(
          new L3(V3.O, V(-1, 0, 0), -4096, 4096),
          V(10, 0, 0),
          V3.O,
          -10,
          0
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Z, 10),
        V3.Y,
        V(-1, 0, 0),
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(10, 12, 10), V(0, -1, 0), -4096, 4096),
          V(10, 10, 10),
          V(10, 12, 10),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(0, 12, 10), V3.X, -4096, 4096),
          V(10, 12, 10),
          V(0, 12, 10),
          10,
          0
        ),
        new StraightEdge(
          new L3(V(0, 0, 10), V3.Y, -4096, 4096),
          V(0, 12, 10),
          V(0, 10, 10),
          12,
          10
        ),
        new StraightEdge(
          new L3(V(0, 10, 10), V3.X, -4096, 4096),
          V(0, 10, 10),
          V(10, 10, 10),
          0,
          10
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Y, 12),
        V(-1, 0, 0),
        V3.Z,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(0, 12, -2), V3.X, -4096, 4096),
          V(10, 12, -2),
          V(0, 12, -2),
          10,
          0
        ),
        new StraightEdge(
          new L3(V(0, 12, -2), V3.Z, -4096, 4096),
          V(0, 12, -2),
          V(0, 12, 10),
          0,
          12
        ),
        new StraightEdge(
          new L3(V(0, 12, 10), V3.X, -4096, 4096),
          V(0, 12, 10),
          V(10, 12, 10),
          0,
          10
        ),
        new StraightEdge(
          new L3(V(10, 12, -2), V3.Z, -4096, 4096),
          V(10, 12, 10),
          V(10, 12, -2),
          12,
          0
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(0, 0, -1), 2),
        V3.Y,
        V3.X,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(0, 0, -2), V3.Y, -4096, 4096),
          V(0, 0, -2),
          V(0, 12, -2),
          0,
          12
        ),
        new StraightEdge(
          new L3(V(0, 12, -2), V3.X, -4096, 4096),
          V(0, 12, -2),
          V(10, 12, -2),
          0,
          10
        ),
        new StraightEdge(
          new L3(V(10, 12, -2), V(0, -1, 0), -4096, 4096),
          V(10, 12, -2),
          V(10, 0, -2),
          0,
          12
        ),
        new StraightEdge(
          new L3(V(10, 0, -2), V(-1, 0, 0), -4096, 4096),
          V(10, 0, -2),
          V(0, 0, -2),
          0,
          10
        ),
      ],
      []
    ),
  ],
  false
)
