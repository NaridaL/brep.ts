import {
  arrayCopyBlocks,
  assertf,
  getLast,
  int,
  M4,
  Matrix,
  V3,
  Vector,
} from "ts3dutils"

type DHParams = { theta: number; d: number; r: number; alpha: number }

function DH({ theta, d, r, alpha }: DHParams): M4 {
  return M4.product(M4.rotateZ(theta), M4.translate(r, 0, d), M4.rotateX(alpha))
}

function MatrixsetSub(this: Matrix, x: int, y: int, sub: Matrix) {
  assertf(() => x + sub.width <= this.width)
  assertf(() => y + sub.height <= this.height)
  arrayCopyBlocks(
    sub.m,
    0,
    sub.width,
    this.m,
    y * this.width + x,
    this.width,
    sub.width,
    sub.height,
  )
}

function MatrixfromMultiRows(...rows: (Vector | V3 | Matrix | number)[][]) {
  const width: (x: Vector | V3 | Matrix | number) => int = (x) =>
    x instanceof Matrix ? x.width : 1
  const height: (x: Vector | V3 | Matrix | number) => int = (x) =>
    x instanceof Matrix
      ? x.height
      : x instanceof Vector
      ? x.dim()
      : x instanceof V3
      ? 3
      : 1
  const result = Matrix.forWidthHeight(
    rows[0].reduce<int>((sum: int, val) => sum + width(val), 0),
    rows.reduce((sum: int, val) => sum + height(val[0]), 0),
  )
  let y = 0
  for (const row of rows) {
    let x = 0
    for (const val of row) {
      MatrixsetSub.call(result, x, y, val)
      x += height(val)
    }
    y += height(row[0])
  }
  return result
}

function JacobiCol(paramss: DHParams[], isRot: boolean[]) {
  const stack = [M4.IDENTITY]
  let curr = stack[0]
  for (const params of paramss) {
    curr = curr.times(DH(params))
    stack.push(curr)
  }
  const finalPos = getLast(stack).getTranslation()
  const moves: V3[] = []
  const rots: V3[] = []
  for (let i = 0; i < paramss.length; i++) {
    moves.push(isRot ? stack[i].Z : V3.O)
    rots.push(
      isRot
        ? stack[i].Z.cross(finalPos.minus(stack[i].getTranslation()))
        : stack[i].Z,
    )
  }
  MatrixfromMultiRows(moves, rots)
  if (isRot) {
    return
  }
}
