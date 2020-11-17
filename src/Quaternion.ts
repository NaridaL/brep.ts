import { assertf, eq, floatHashCode, M4, raddd, V3 } from "ts3dutils"

import { acos, cos, sin } from "./math"

export class Quaternion {
  private constructor(
    public readonly s: number,
    public readonly x: number,
    public readonly y: number,
    public readonly z: number,
  ) {}

  static O = new Quaternion(1, 0, 0, 0)

  static axis(axis: V3, rotation: raddd) {
    assertf(() => axis.hasLength(1))
    return new Quaternion(
      cos(rotation / 2),
      sin(rotation / 2) * axis.x,
      sin(rotation / 2) * axis.y,
      sin(rotation / 2) * axis.z,
    )
  }

  static of(s: number, x: number, y: number, z: number) {
    return new Quaternion(s, x, y, z)
  }

  public plus(q: Quaternion) {
    return new Quaternion(
      this.s + q.s,
      this.x + q.x,
      this.y + q.y,
      this.z + q.z,
    )
  }

  public times(q: Quaternion | number) {
    return "number" == typeof q
      ? new Quaternion(q * this.s, q * this.x, q * this.y, q * this.z)
      : new Quaternion(
          this.s * q.s - (this.x * q.x + this.y * q.y + this.z * q.z),
          this.y * q.z - this.z * q.y + this.s * q.x + q.s * this.x,
          this.z * q.x - this.x * q.z + this.s * q.y + q.s * this.y,
          this.x * q.y - this.y * q.x + this.s * q.z + q.s * this.z,
        )
  }

  public conjugated() {
    return new Quaternion(this.s, -this.x, -this.y, -this.z)
  }

  public length() {
    return Math.hypot(this.s, this.x, this.y, this.z)
  }

  public norm() {
    return this.s ** 2 + this.x ** 2 + (this.y ** 2 + this.z ** 2)
  }

  public unit() {
    const l = this.length()
    return new Quaternion(this.s / l, this.x / l, this.y / l, this.z / l)
  }

  public inverse() {
    return this.conjugated().times(1 / this.norm())
  }

  public toM4() {
    assertf(() => eq(1, this.length()))
    const { s, x, y, z } = this
    // prettier-ignore
    return new M4([
      1 - 2 * (y * y + z * z), 2 * (x * y - z * s), 2 * (x * z + y * s), 0,
      2 * (x * y + z * s), 1 - 2 * (x * x + z * z), 2 * (y * z - x * s), 0,
      2 * (x * z - y * s), 2 * (y * z + x * s), 1 - 2 * (x * x + y * y), 0,
      0, 0, 0, 1,
    ])
  }

  static fromRotation(m4: M4) {
    const sqrtTracePlus1 = Math.sqrt(m4.trace() + 1)
    const f = 1 / (2 * sqrtTracePlus1)
    return new Quaternion(
      sqrtTracePlus1 / 2,
      f * (m4.e(2, 1) - m4.e(1, 2)),
      f * (m4.e(0, 2) - m4.e(2, 0)),
      f * (m4.e(1, 0) - m4.e(0, 1)),
    )
  }

  public rotatePoint(p: V3) {
    const v = this.times(Quaternion.of(1, p.x, p.y, p.z)).times(
      this.conjugated(),
    )
    return new V3(v.x, v.y, v.z)
  }

  public like(q: Quaternion, precision?: number) {
    return (
      eq(this.s, q.s, precision) &&
      eq(this.x, q.x, precision) &&
      eq(this.y, q.y, precision) &&
      eq(this.z, q.z, precision)
    )
  }

  public equals(q: any) {
    return (
      this == q ||
      (q instanceof Quaternion &&
        this.s == q.s &&
        this.x == q.x &&
        this.y == q.y &&
        this.z == q.z)
    )
  }

  public hashCode() {
    let hashCode = 0
    hashCode = (hashCode * 31 + floatHashCode(this.s)) | 0
    hashCode = (hashCode * 31 + floatHashCode(this.x)) | 0
    hashCode = (hashCode * 31 + floatHashCode(this.y)) | 0
    hashCode = (hashCode * 31 + floatHashCode(this.z)) | 0
    return hashCode
  }

  public slerp(b: Quaternion, f: number) {
    assertf(() => eq(1, this.length()))
    assertf(() => eq(1, b.length()))
    const a = this
    let dot = a.s * b.s + a.x * b.x + a.y * b.y + a.z * b.z

    if (dot < 0) {
      dot = -dot
      b = b.times(-1)
      console.log("dot < 0")
    }

    const DOT_THRESHOLD = 0.9995
    if (dot > DOT_THRESHOLD) {
      // If the inputs are too close for comfort, linearly interpolate
      // and normalize the result.
      return a
        .times(1 - f)
        .plus(b.times(f))
        .unit()
    }

    // Since dot is in range [0, DOT_THRESHOLD], acos is safe
    const theta0 = acos(dot) // theta_0 = angle between input vectors
    const theta = theta0 * f // theta = angle between v0 and result

    const s0 = cos(theta) - (dot * sin(theta)) / sin(theta0) // == sin(theta_0 - theta) / sin(theta_0)
    const s1 = sin(theta) / sin(theta0)
    console.log(s0, s1, a.times(s0), b.times(s1))
    return a.times(s0).plus(b.times(s1))
  }

  toArray(): [number, number, number, number] {
    return [this.s, this.x, this.y, this.z]
  }
}
