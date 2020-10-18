import { V3 } from "ts3dutils"

export interface WebGLDebugable {
  debugInfo?(): {
    points?: V3[]
    lines?: V3[]
  }
}
