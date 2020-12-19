function drawPSI(
  ctx: CanvasRenderingContext2D,
  ps: ParametricSurface,
  is: ImplicitSurface,
) {
  const pf = ps.pSTFunc(),
    icc = is.implicitFunction()
  const implicitCurve = (x: number, y: number) => icc(pf(x, y))
  const sStep = 0.1,
    tStep = 0.05
  const dpds = ps.dpds()
  const dpdt = ps.dpdt()
  const ist = (x: number, y: number) => icc(pf(x, y))
  const dids = (s: number, t: number) => didp(pf(s, t)).dot(dpds(s, t))
  const didt = (s: number, t: number) => didp(pf(s, t)).dot(dpdt(s, t))
  const mf = MathFunctionR2R.forFFxFy(ist, dids, didt)
  const didp = is.didp.bind(is)
  let { sMin, tMin, sMax, tMax } = ps
  //sMin -= 1

  const deltaS = sMax - sMin,
    deltaT = tMax - tMin
  const sRes = ceil(deltaS / sStep),
    tRes = ceil(deltaT / tStep)
  const result: { points: V3[]; tangents: V3[] }[] = []
  const bounds2 = (s: number, t: number) => pf(s, t).y > 0
  const test = MathFunctionR2R.forNerdamer("x * y * (-x - y + 1.5) - 0.0")
  const flower = MathFunctionR2R.forNerdamer(
    "(3x^2 - y^2)^2 y^2 - (x^2 +y^2)^4",
  )
  const implicit = MathFunctionR2R.forNerdamer("x^3 + 2x - 3x* y - y^2")
  nerdamer.setFunction(
    "cassini",
    "acxy".split(""),
    "(x^2 + y^2)^2 + 2 c^2 (x^2 - y^2) - (a^4 - c^4)",
  )
  const cassini = MathFunctionR2R.forNerdamer("cassini(1,1,x,y)")
  const heart = MathFunctionR2R.forNerdamer("(x^2+(-y)^2-1)^3-x^2 (-y)^3")
  nerdamer.setConstant("pi", PI)
  const grid = MathFunctionR2R.forNerdamer("cos(pi* y) - cos(pi *x)")
  const what = mf
  // draw3(ctx, implicitCurve, dids, didt, ps, sStep, tStep){sMin: -2, sMax: 2, tMin: -1.7, tMax: 1.3}
  draw3(ctx, what, what.x, what.y, ps, 0.07, 0.07)
}
function draw3(ctx, implicitCurve, dids, didt, bounds, sStep, tStep) {
  ctx.save()
  draw(ctx, implicitCurve, bounds, sStep, tStep)
  ctx.restore()

  // ctx.save()
  // ctx.translate(0, 300)
  // draw(ctx, dids, bounds, sStep, tStep)
  // ctx.restore()

  // ctx.save()
  // ctx.translate(0, 600)
  // draw(ctx, didt, bounds, sStep, tStep)
  // ctx.restore()
}
const colorScale: any = chroma
  .scale(["white", "red", "green", "white"])
  .domain([0, 127, 128, 255])
const colors = arrayFromFunction(256, (i) => colorScale(i))
const colorsBright = colors.map((c) => c.darken(0.5))
function draw(
  ctx: CanvasRenderingContext2D,
  ic: MathFunctionR2R,
  bounds: { sMin; tMin; sMax; tMax },
  sStep,
  tStep,
) {
  let { sMin, tMin, sMax, tMax } = bounds
  const deltaS = sMax - sMin,
    deltaT = tMax - tMin
  const sRes = ceil(deltaS / sStep),
    tRes = ceil(deltaT / tStep)
  const scale = 16
  const grid = new Array(sRes * tRes).fill(0)
  arrayFromFunction(tRes, (i) =>
    grid
      .slice(sRes * i, sRes * (i + 1))
      .map((v) => (v ? "X" : "_"))
      .join(""),
  ).join("\n")
  ctx.scale(scale, scale)
  const at = (i: int, j: int) => grid[j * sRes + i]
  const set = (i: int, j: int) =>
    0 <= i && i < sRes && 0 <= j && j < tRes && (grid[j * sRes + i] = 1)
  const logTable = []
  for (let i = 0; i < sRes; i++) {
    search: for (let j = 0; j < tRes; j++) {
      // if (at(i, j)) continue
      // set(i, j)
      ctx.fillStyle = (i + j) % 2 == 0 ? "#ff86ee" : "#ee9269"
      ctx.fillRect(i, j, 1, 1)
      //let s = sMin + i * sStep, t = tMin + j * tStep
      //const startS = s, startT = t
      //// basically curvePoint
      //for (let k = 0; k < 8; k++) {
      //	const fp = implicitCurve(s, t)
      //	const dfpdx = dids(s, t), dfpdy = didt(s, t)
      //	if (0 == dfpdx * dfpdx + dfpdy * dfpdy) {
      //		// top of a hill, keep looking
      //		continue search
      //	}
      //	const scale = fp / (dfpdx * dfpdx + dfpdy * dfpdy)
      //	s -= scale * dfpdx
      //	t -= scale * dfpdy
      //}
      //const li = floor((s - sMin) / sStep), lj = floor((t - tMin) / tStep)
      //logTable.push({
      //	i,
      //	j,
      //	li,
      //	lj,
      //	startS,
      //	startT,
      //	s,
      //	t,
      //	'bounds(s, t)': bounds(s, t),
      //	'ic(s,t)': implicitCurve(s, t)
      //})
      //if (!(i == li && j == lj) && at(li, lj)) {
      //	continue search
      //}
      ////if (0 <= li && li < sRes && 0 <= lj && lj < tRes) {
      ////		set(li, lj)
      ////	}
      //// s, t are now good starting coordinates to use follow algo
      //if (bounds(s, t) && eq0(implicitCurve(s, t))) {
      //	console.log(V(s, t).sce)
      //	const subresult = mkcurves(implicitCurve, s, t, stepSize, dids, didt, bounds)
      //	for (const curvedata of subresult) {
      //		for (const {x, y} of curvedata.points) {
      //			const lif = (x - sMin) / sStep, ljf = (y - tMin) / tStep
      //			set((lif - 0.5) | 0, (ljf - 0.5) | 0)
      //			set((lif - 0.5) | 0, (ljf + 0.5) | 0)
      //			set((lif + 0.5) | 0, (ljf - 0.5) | 0)
      //			set((lif + 0.5) | 0, (ljf + 0.5) | 0)
      //		}
      //	}
      //	result.push(...subresult)
      //}
    }
  }

  for (let i = 0; i < sRes; i += 1 / scale) {
    for (let j = 0; j < tRes; j += 1 / scale) {
      let s = sMin + i * sStep,
        t = tMin + j * tStep
      const icv = ic(s, t)
      //console.log(icv)
      const grid = ((i | 0) + (j | 0)) % 2 == 0
      // const valueColor = colorScale(icv)
      // ctx.fillStyle = grid ? valueColor.brighten() : valueColor
      const colorIndex = clamp(round((icv + 0.5) * 255), 0, 255)
      ctx.fillStyle = grid ? colorsBright[colorIndex] : colors[colorIndex]
      //ctx.fillStyle = bounds2(s, t) ? (grid ? valueColor.brighten() : valueColor) : 'black'

      //ctx.fillStyle = 'red'
      ctx.fillRect(i, j, 1 / scale, 1 / scale)
    }
  }

  for (let i = 0; i < sRes; i += 1) {
    for (let j = 0; j < tRes; j += 1) {
      let s = sMin + (i + 0.5) * sStep,
        t = tMin + (j + 0.5) * tStep
      const icv = ic(s, t)
      const gradient = V(ic.x(s, t), ic.y(s, t)).toLength(0.4)
      ctx.lineWidth = 1 / scale
      ctx.moveTo(i + 0.5 - gradient.x, j + 0.5 - gradient.y)
      ctx.lineTo(i + 0.5 + gradient.x, j + 0.5 + gradient.y)
    }
  }
  ctx.stroke()

  ctx.fillStyle = "black"
  ctx.font = "1px Consolas"
  for (let i = 0; i < sRes; i += 1) {
    let s = sMin + i * sStep
    ctx.save()
    ctx.translate(i, tRes)
    ctx.rotate(90 * DEG)
    ctx.fillText("s " + round10(s, -4), 0, 0)
    ctx.restore()
  }
  for (let j = 0; j < tRes; j += 1) {
    let t = tMin + j * tStep
    ctx.fillText("t " + round10(t, -4), sRes, j + 1)
  }
}
function loadingMain() {
  const canvas = document.getElementById("c") as HTMLCanvasElement
  const ctx = canvas.getContext("2d")
  canvas.height = window.innerHeight
  canvas.width = window.innerWidth
  let frameCount = 0
  const is = new SemiEllipsoidSurface(V3.O, V3.X, V3.Y, V3.Z)
  const ps = new SemiCylinderSurface(
    new SemiEllipseCurve(
      V(0.5, 0, -2),
      V(0.5, 0, 0),
      V(0, 0.05, 0),
      0,
      3.141592653589793,
    ),
    V(0, 0, -1),
    0,
    3.141592653589793,
    -4,
    0,
  )
  const stepSize = 0.02
  //const is = SemiEllipsoidSurface.UNIT
  //const ps = new ProjectedCurveSurface(new BezierCurve(V(0.30000000000000004, -0.1, 1.2), V(0.30000000000000004,
  // 0.010456949966158674, 1.2), V(0.2104569499661587, 0.1, 1.2), V(0.10000000000000002, 0.1, 1.2), 0, 1), V(0, 0,
  // 2), 0, 1, -1, 0)

  const ps3 = SemiCylinderSurface.UNIT.rotateZ(-20 * DEG)
    .scale(0.5, 0.05, 4)
    .translate(0.5 - 0.0, 0, -2)
    .flipped()
  requestAnimationFrame((time) => drawPSI(ctx, ps, is))
}
