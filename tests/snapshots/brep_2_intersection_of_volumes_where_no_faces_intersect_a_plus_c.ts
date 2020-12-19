new BRep(
  [
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.X, 1),
        V(0, -1, 0),
        V(0, 0, -1),
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(1, 3, 1), V3.Z, 0, 2),
          V(1, 3, 1),
          V(1, 3, 3),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(1, 1, 3), V3.Y, 0, 2),
          V(1, 3, 3),
          V(1, 1, 3),
          2,
          0
        ),
        new StraightEdge(new L3(V3.XYZ, V3.Z, 0, 2), V(1, 1, 3), V3.XYZ, 2, 0),
        new StraightEdge(new L3(V3.XYZ, V3.Y, 0, 2), V3.XYZ, V(1, 3, 1), 0, 2),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(0, -1, 0), -3),
        V(-1, 0, 0),
        V(0, 0, -1),
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(3, 3, 1), V3.Z, 0, 2),
          V(3, 3, 1),
          V(3, 3, 3),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(1, 3, 3), V3.X, 0, 2),
          V(3, 3, 3),
          V(1, 3, 3),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(1, 3, 1), V3.Z, 0, 2),
          V(1, 3, 3),
          V(1, 3, 1),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(1, 3, 1), V3.X, 0, 2),
          V(1, 3, 1),
          V(3, 3, 1),
          0,
          2
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(-1, 0, 0), -3),
        V3.Y,
        V(0, 0, -1),
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(3, 1, 1), V3.Z, 0, 2),
          V(3, 1, 1),
          V(3, 1, 3),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(3, 3, 3), V(0, -1, 0), 0, 2),
          V(3, 1, 3),
          V(3, 3, 3),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(3, 3, 1), V3.Z, 0, 2),
          V(3, 3, 3),
          V(3, 3, 1),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(3, 3, 1), V(0, -1, 0), 0, 2),
          V(3, 3, 1),
          V(3, 1, 1),
          0,
          2
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Y, 1),
        V3.X,
        V(0, 0, -1),
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(new L3(V3.XYZ, V3.Z, 0, 2), V3.XYZ, V(1, 1, 3), 0, 2),
        new StraightEdge(
          new L3(V(3, 1, 3), V(-1, 0, 0), 0, 2),
          V(1, 1, 3),
          V(3, 1, 3),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(3, 1, 1), V3.Z, 0, 2),
          V(3, 1, 3),
          V(3, 1, 1),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(3, 1, 1), V(-1, 0, 0), 0, 2),
          V(3, 1, 1),
          V3.XYZ,
          0,
          2
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
        new StraightEdge(
          new L3(V(3, 1, 1), V(-1, 0, 0), 0, 2),
          V3.XYZ,
          V(3, 1, 1),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(3, 3, 1), V(0, -1, 0), 0, 2),
          V(3, 1, 1),
          V(3, 3, 1),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(1, 3, 1), V3.X, 0, 2),
          V(3, 3, 1),
          V(1, 3, 1),
          2,
          0
        ),
        new StraightEdge(new L3(V3.XYZ, V3.Y, 0, 2), V(1, 3, 1), V3.XYZ, 2, 0),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(0, 0, -1), -3),
        V3.Y,
        V3.X,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(1, 1, 3), V3.Y, 0, 2),
          V(1, 1, 3),
          V(1, 3, 3),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(1, 3, 3), V3.X, 0, 2),
          V(1, 3, 3),
          V(3, 3, 3),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(3, 3, 3), V(0, -1, 0), 0, 2),
          V(3, 3, 3),
          V(3, 1, 3),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(3, 1, 3), V(-1, 0, 0), 0, 2),
          V(3, 1, 3),
          V(1, 1, 3),
          0,
          2
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.X, 20),
        V(0, -1, 0),
        V(0, 0, -1),
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(20, 2, 0), V3.Z, 0, 2),
          V(20, 2, 0),
          V(20, 2, 2),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(20, 0, 2), V3.Y, 0, 2),
          V(20, 2, 2),
          V(20, 0, 2),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(20, 0, 0), V3.Z, 0, 2),
          V(20, 0, 2),
          V(20, 0, 0),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(20, 0, 0), V3.Y, 0, 2),
          V(20, 0, 0),
          V(20, 2, 0),
          0,
          2
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(0, -1, 0), -2),
        V(-1, 0, 0),
        V(0, 0, -1),
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(22, 2, 0), V3.Z, 0, 2),
          V(22, 2, 0),
          V(22, 2, 2),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(20, 2, 2), V3.X, 0, 2),
          V(22, 2, 2),
          V(20, 2, 2),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(20, 2, 0), V3.Z, 0, 2),
          V(20, 2, 2),
          V(20, 2, 0),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(20, 2, 0), V3.X, 0, 2),
          V(20, 2, 0),
          V(22, 2, 0),
          0,
          2
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(-1, 0, 0), -22),
        V3.Y,
        V(0, 0, -1),
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(22, 0, 0), V3.Z, 0, 2),
          V(22, 0, 0),
          V(22, 0, 2),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(22, 2, 2), V(0, -1, 0), 0, 2),
          V(22, 0, 2),
          V(22, 2, 2),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(22, 2, 0), V3.Z, 0, 2),
          V(22, 2, 2),
          V(22, 2, 0),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(22, 2, 0), V(0, -1, 0), 0, 2),
          V(22, 2, 0),
          V(22, 0, 0),
          0,
          2
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V3.Y, 0),
        V3.X,
        V(0, 0, -1),
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(20, 0, 0), V3.Z, 0, 2),
          V(20, 0, 0),
          V(20, 0, 2),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(22, 0, 2), V(-1, 0, 0), 0, 2),
          V(20, 0, 2),
          V(22, 0, 2),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(22, 0, 0), V3.Z, 0, 2),
          V(22, 0, 2),
          V(22, 0, 0),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(22, 0, 0), V(-1, 0, 0), 0, 2),
          V(22, 0, 0),
          V(20, 0, 0),
          0,
          2
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
          new L3(V(22, 0, 0), V(-1, 0, 0), 0, 2),
          V(20, 0, 0),
          V(22, 0, 0),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(22, 2, 0), V(0, -1, 0), 0, 2),
          V(22, 0, 0),
          V(22, 2, 0),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(20, 2, 0), V3.X, 0, 2),
          V(22, 2, 0),
          V(20, 2, 0),
          2,
          0
        ),
        new StraightEdge(
          new L3(V(20, 0, 0), V3.Y, 0, 2),
          V(20, 2, 0),
          V(20, 0, 0),
          2,
          0
        ),
      ],
      []
    ),
    new PlaneFace(
      new PlaneSurface(
        new P3(V(0, 0, -1), -2),
        V3.Y,
        V3.X,
        -100,
        100,
        -100,
        100
      ),
      [
        new StraightEdge(
          new L3(V(20, 0, 2), V3.Y, 0, 2),
          V(20, 0, 2),
          V(20, 2, 2),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(20, 2, 2), V3.X, 0, 2),
          V(20, 2, 2),
          V(22, 2, 2),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(22, 2, 2), V(0, -1, 0), 0, 2),
          V(22, 2, 2),
          V(22, 0, 2),
          0,
          2
        ),
        new StraightEdge(
          new L3(V(22, 0, 2), V(-1, 0, 0), 0, 2),
          V(22, 0, 2),
          V(20, 0, 2),
          0,
          2
        ),
      ],
      []
    ),
  ],
  false
)
