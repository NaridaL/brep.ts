import { V3 } from "ts3dutils"
import { Surface } from ".."

export abstract class ImplicitSurface extends Surface {
  static is(obj: any): obj is ImplicitSurface {
    return obj.implicitFunction && obj.didp
  }

  abstract implicitFunction(): (pWC: V3) => number

  /**
   * partial derivatives of this.implicitFunction in point pWC
   * @param pWC
   * @return
   */
  abstract didp(pWC: V3): V3
}
