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
        100,
      ),
      [
        new StraightEdge(
          new L3(V3.Y, V3.Z, -4096, 4096),
          V3.Y,
          V(0, 1, 5),
          0,
          5,
        ),
        new StraightEdge(
          new L3(V(0, 0, 5), V3.Y, -4096, 4096),
          V(0, 1, 5),
          V(0, 5, 5),
          1,
          5,
        ),
        new StraightEdge(
          new L3(V(0, 5, 0), V3.Z, -4096, 4096),
          V(0, 5, 5),
          V(0, 5, 0),
          5,
          0,
        ),
        new StraightEdge(
          new L3(V3.O, V3.Y, -4096, 4096),
          V(0, 5, 0),
          V3.Y,
          5,
          1,
        ),
      ],
      [],
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(0, -1, 0), 0),
        V3.X,
        V3.Z,
        -100,
        100,
        -100,
        100,
      ),
      [
        new StraightEdge(
          new L3(V3.X, V(0, 0, -1), -4096, 4096),
          V(1, 0, 5),
          V3.X,
          -5,
          0,
        ),
        new StraightEdge(
          new L3(V(5, 0, 0), V(-1, 0, 0), -4096, 4096),
          V3.X,
          V(5, 0, 0),
          4,
          0,
        ),
        new StraightEdge(
          new L3(V(5, 0, 0), V3.Z, -4096, 4096),
          V(5, 0, 0),
          V(5, 0, 5),
          0,
          5,
        ),
        new StraightEdge(
          new L3(V(5, 0, 5), V(-1, 0, 0), -4096, 4096),
          V(5, 0, 5),
          V(1, 0, 5),
          0,
          4,
        ),
      ],
      [],
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(0, 0, -1), 0),
        V3.Y,
        V3.X,
        -100,
        100,
        -100,
        100,
      ),
      [
        new StraightEdge(
          new L3(V3.Y, V(-1, 0, 0), -4096, 4096),
          V(1, 1, 0),
          V3.Y,
          -1,
          0,
        ),
        new StraightEdge(
          new L3(V3.O, V3.Y, -4096, 4096),
          V3.Y,
          V(0, 5, 0),
          1,
          5,
        ),
        new StraightEdge(
          new L3(V(0, 5, 0), V3.X, -4096, 4096),
          V(0, 5, 0),
          V(5, 5, 0),
          0,
          5,
        ),
        new StraightEdge(
          new L3(V(5, 5, 0), V(0, -1, 0), -4096, 4096),
          V(5, 5, 0),
          V(5, 0, 0),
          0,
          5,
        ),
        new StraightEdge(
          new L3(V(5, 0, 0), V(-1, 0, 0), -4096, 4096),
          V(5, 0, 0),
          V3.X,
          0,
          4,
        ),
        new StraightEdge(
          new L3(V3.X, V3.Y, -4096, 4096),
          V3.X,
          V(1, 1, 0),
          0,
          1,
        ),
      ],
      [],
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Z, 5),
        V3.Y,
        V(-1, 0, 0),
        -100,
        100,
        -100,
        100,
      ),
      [
        new StraightEdge(
          new L3(V(0, 1, 5), V3.X, -4096, 4096),
          V(0, 1, 5),
          V(1, 1, 5),
          0,
          1,
        ),
        new StraightEdge(
          new L3(V(1, 0, 5), V(0, -1, 0), -4096, 4096),
          V(1, 1, 5),
          V(1, 0, 5),
          -1,
          0,
        ),
        new StraightEdge(
          new L3(V(5, 0, 5), V(-1, 0, 0), -4096, 4096),
          V(1, 0, 5),
          V(5, 0, 5),
          4,
          0,
        ),
        new StraightEdge(
          new L3(V(5, 5, 5), V(0, -1, 0), -4096, 4096),
          V(5, 0, 5),
          V(5, 5, 5),
          5,
          0,
        ),
        new StraightEdge(
          new L3(V(0, 5, 5), V3.X, -4096, 4096),
          V(5, 5, 5),
          V(0, 5, 5),
          5,
          0,
        ),
        new StraightEdge(
          new L3(V(0, 0, 5), V3.Y, -4096, 4096),
          V(0, 5, 5),
          V(0, 1, 5),
          5,
          1,
        ),
      ],
      [],
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Y, 5),
        V(-1, 0, 0),
        V3.Z,
        -100,
        100,
        -100,
        100,
      ),
      [
        new StraightEdge(
          new L3(V(0, 5, 0), V3.X, -4096, 4096),
          V(5, 5, 0),
          V(0, 5, 0),
          5,
          0,
        ),
        new StraightEdge(
          new L3(V(0, 5, 0), V3.Z, -4096, 4096),
          V(0, 5, 0),
          V(0, 5, 5),
          0,
          5,
        ),
        new StraightEdge(
          new L3(V(0, 5, 5), V3.X, -4096, 4096),
          V(0, 5, 5),
          V(5, 5, 5),
          0,
          5,
        ),
        new StraightEdge(
          new L3(V(5, 5, 0), V3.Z, -4096, 4096),
          V(5, 5, 5),
          V(5, 5, 0),
          5,
          0,
        ),
      ],
      [],
    ),
    new PlaneFace(
      new PlaneSurface(new P3(V3.X, 5), V3.Y, V3.Z, -100, 100, -100, 100),
      [
        new StraightEdge(
          new L3(V(5, 5, 0), V(0, -1, 0), -4096, 4096),
          V(5, 0, 0),
          V(5, 5, 0),
          5,
          0,
        ),
        new StraightEdge(
          new L3(V(5, 5, 0), V3.Z, -4096, 4096),
          V(5, 5, 0),
          V(5, 5, 5),
          0,
          5,
        ),
        new StraightEdge(
          new L3(V(5, 5, 5), V(0, -1, 0), -4096, 4096),
          V(5, 5, 5),
          V(5, 0, 5),
          0,
          5,
        ),
        new StraightEdge(
          new L3(V(5, 0, 0), V3.Z, -4096, 4096),
          V(5, 0, 5),
          V(5, 0, 0),
          5,
          0,
        ),
      ],
      [],
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(0, -1, 0), -1),
        V3.X,
        V3.Z,
        -100,
        100,
        -100,
        100,
      ),
      [
        new StraightEdge(
          new L3(V(0, 1, 5), V3.X, -4096, 4096),
          V(1, 1, 5),
          V(0, 1, 5),
          1,
          0,
        ),
        new StraightEdge(
          new L3(V3.Y, V3.Z, -4096, 4096),
          V(0, 1, 5),
          V3.Y,
          5,
          0,
        ),
        new StraightEdge(
          new L3(V3.Y, V(-1, 0, 0), -4096, 4096),
          V3.Y,
          V(1, 1, 0),
          0,
          -1,
        ),
        new StraightEdge(
          new L3(V(1, 1, 0), V3.Z, -4096, 4096),
          V(1, 1, 0),
          V(1, 1, 5),
          0,
          5,
        ),
      ],
      [],
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(-1, 0, 0), -1),
        V(0, -1, 0),
        V3.Z,
        -100,
        100,
        -100,
        100,
      ),
      [
        new StraightEdge(
          new L3(V(1, 0, 5), V(0, -1, 0), -4096, 4096),
          V(1, 0, 5),
          V(1, 1, 5),
          0,
          -1,
        ),
        new StraightEdge(
          new L3(V(1, 1, 0), V3.Z, -4096, 4096),
          V(1, 1, 5),
          V(1, 1, 0),
          5,
          0,
        ),
        new StraightEdge(
          new L3(V3.X, V3.Y, -4096, 4096),
          V(1, 1, 0),
          V3.X,
          1,
          0,
        ),
        new StraightEdge(
          new L3(V3.X, V(0, 0, -1), -4096, 4096),
          V3.X,
          V(1, 0, 5),
          0,
          -5,
        ),
      ],
      [],
    ),
  ],
  false,
)