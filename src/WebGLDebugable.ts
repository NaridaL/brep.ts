import { V3 } from "ts3dutils"

export interface WebGLDebugable {
  debugInfo?(): {
    ps?: V3[]
    lines?: V3[]
  }
}
