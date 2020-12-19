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
          new L3(V(3, 10, 0), V3.Z, -4096, 4096),
          V(3, 10, 0),
          V(3, 10, 4),
          0,
          4
        ),
        new StraightEdge(
          new L3(V(0, 10, 4), V3.X, -4096, 4096),
          V(3, 10, 4),
          V(10, 10, 4),
          3,
          10
        ),
        new StraightEdge(
          new L3(V(10, 10, 0), V3.Z, 0, 5),
          V(10, 10, 4),
          V(10, 10, 0),
          4,
          0
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
      new PlaneSurface(new P3(V3.X, 10), V3.Y, V3.Z, -100, 100, -100, 100),
      [
        new StraightEdge(
          new L3(V(10, 3, 0), V(0, 0, -1), -4096, 4096),
          V(10, 3, 4),
          V(10, 3, 0),
          -4,
          0
        ),
        new StraightEdge(
          new L3(V(10, 0, 0), V3.Y, -4096, 4096),
          V(10, 3, 0),
          V(10, 10, 0),
          3,
          10
        ),
        new StraightEdge(
          new L3(V(10, 10, 0), V3.Z, 0, 5),
          V(10, 10, 0),
          V(10, 10, 4),
          0,
          4
        ),
        new StraightEdge(
          new L3(V(10, 0, 4), V(0, -1, 0), -4096, 4096),
          V(10, 10, 4),
          V(10, 3, 4),
          -10,
          -3
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
          new L3(V(0, 3, 0), V(-1, 0, 0), -4096, 4096),
          V(10, 3, 0),
          V(3, 3, 0),
          -10,
          -3
        ),
        new StraightEdge(
          new L3(V(3, 0, 0), V3.Y, -4096, 4096),
          V(3, 3, 0),
          V(3, 10, 0),
          3,
          10
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(-1, 0, 0), -3),
        V(0, -1, 0),
        V3.Z,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(3, 10, 0), V3.Z, -4096, 4096),
          V(3, 10, 4),
          V(3, 10, 0),
          4,
          0
        ),
        new StraightEdge(
          new L3(V(3, 0, 0), V3.Y, -4096, 4096),
          V(3, 10, 0),
          V(3, 3, 0),
          10,
          3
        ),
        new StraightEdge(
          new L3(V(3, 3, 0), V3.Z, 0, 4),
          V(3, 3, 0),
          V(3, 3, 4),
          0,
          4
        ),
        new StraightEdge(
          new L3(V(3, 3, 4), V3.Y, 0, 10),
          V(3, 3, 4),
          V(3, 10, 4),
          0,
          7
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(0, -1, 0), -3),
        V3.X,
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
          V(10, 3, 4),
          0,
          -4
        ),
        new StraightEdge(
          new L3(V(13, 3, 4), V(-1, 0, 0), 0, 10),
          V(10, 3, 4),
          V(3, 3, 4),
          3,
          10
        ),
        new StraightEdge(
          new L3(V(3, 3, 0), V3.Z, 0, 4),
          V(3, 3, 4),
          V(3, 3, 0),
          4,
          0
        ),
        new StraightEdge(
          new L3(V(0, 3, 0), V(-1, 0, 0), -4096, 4096),
          V(3, 3, 0),
          V(10, 3, 0),
          -3,
          -10
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Z, 4),
        V3.Y,
        V(-1, 0, 0),
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(0, 10, 4), V3.X, -4096, 4096),
          V(10, 10, 4),
          V(3, 10, 4),
          10,
          3
        ),
        new StraightEdge(
          new L3(V(3, 3, 4), V3.Y, 0, 10),
          V(3, 10, 4),
          V(3, 3, 4),
          7,
          0
        ),
        new StraightEdge(
          new L3(V(13, 3, 4), V(-1, 0, 0), 0, 10),
          V(3, 3, 4),
          V(10, 3, 4),
          10,
          3
        ),
        new StraightEdge(
          new L3(V(10, 0, 4), V(0, -1, 0), -4096, 4096),
          V(10, 3, 4),
          V(10, 10, 4),
          -3,
          -10
        ),
      ],
      []
    ),
  ],
  false
)
