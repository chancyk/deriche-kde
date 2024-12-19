## Deriche KDE
##
## This code is a port of this implementation of the Deriche approximation KDE algorithm found here:
## https://github.com/liuzh-buaa/fast_kde
##
## Which is itself a port of Javascript implementation... which was the implementation of a paper
## that was based on a port of the Deriche computer vision kernel implemented in C by Getreuer.
##
## The original Javascript implementation is here:
## https://github.com/uwdata/fast-kde

import std/[math, strformat]

import derichekde/common


type
  DericheConfig* = object
    sigma*: float
    negative*: bool
    a*: array[0..4, float]
    b_causal*: array[0..3, float]
    b_anticausal*: array[0..4, float]
    sum_causal*: float
    sum_anticausal*: float


proc deriche_causal_coeff*(
  a_out: var array[0..4, float],
  b_out: var array[0..3, float],
  sigma: float
) =
  let K = 4

  # Coefficients from Getreuer's implementation
  let alpha = [
    0.84, 1.8675, 0.84, -1.8675,
    -0.34015, -0.1299, -0.34015, 0.1299
  ]

  let x1 = exp(-1.783 / sigma)
  let x2 = exp(-1.723 / sigma)
  let y1 = 0.6318 / sigma
  let y2 = 1.997 / sigma

  let beta = [
    -x1*cos(y1), x1*sin(y1), -x1*cos(-y1), x1*sin(-y1),
    -x2*cos(y2), x2*sin(y2), -x2*cos(-y2), x2*sin(-y2)
  ]

  let denom = sigma * 2.5066282746310007

  var b: array[0..7, float]
  var a: array[0..9, float]

  # Initialize
  b[0] = alpha[0]
  b[1] = alpha[1]
  for i in 2..7:
    b[i] = 0

  a[0] = 1
  a[1] = 0
  a[2] = beta[0]
  a[3] = beta[1]
  for i in 4..9:
    a[i] = 0

  for k in countup(2, 6, 2):
    b[k] = beta[k]*b[k-2] - beta[k+1]*b[k-1]
    b[k+1] = beta[k]*b[k-1] + beta[k+1]*b[k-2]
    for j in countdown(k-2, 2, 2):
      b[j] += beta[k]*b[j-2] - beta[k+1]*b[j-1]
      b[j+1] += beta[k]*b[j-1] + beta[k+1]*b[j-2]

    # Apply alpha coefficients
    for j in countup(0, k, 2):
      b[j] += alpha[k]*a[j] - alpha[k+1]*a[j+1]
      b[j+1] += alpha[k]*a[j+1] + alpha[k+1]*a[j]

    a[k+2] = beta[k]*a[k] - beta[k+1]*a[k+1]
    a[k+3] = beta[k]*a[k+1] + beta[k+1]*a[k]
    for j in countdown(k, 2, 2):
      a[j] += beta[k]*a[j-2] - beta[k+1]*a[j-1]
      a[j+1] += beta[k]*a[j-1] + beta[k+1]*a[j-2]

  for k in 0..<K:
    let j = k * 2
    b_out[k] = b[j] / denom
    a_out[k+1] = a[j+2]


proc deriche_config*(sigma: float, negative = false): DericheConfig =
  var a: array[0..4, float] = [0,0,0,0,0]
  var bc: array[0..3, float] = [0,0,0,0]

  deriche_causal_coeff(a, bc, sigma)

  let ba = [
    0.0,
    bc[1] - a[1]*bc[0],
    bc[2] - a[2]*bc[0],
    bc[3] - a[3]*bc[0],
    -a[4]*bc[0]
  ]

  let accum_denom = 1.0 + a[1] + a[2] + a[3] + a[4]
  let sum_causal = (bc[0] + bc[1] + bc[2] + bc[3]) / accum_denom
  let sum_anticausal = (ba[1] + ba[2] + ba[3] + ba[4]) / accum_denom

  result = DericheConfig(
    sigma: sigma,
    negative: negative,
    a: a,
    b_causal: bc,
    b_anticausal: ba,
    sum_causal: sum_causal,
    sum_anticausal: sum_anticausal
  )


proc deriche_init_zero_pad*(
  dest: var openArray[float],
  src: openArray[float],
  N: int, stride: int,
  b: openArray[float], p: int,
  a: openArray[float], q: int,
  sum_value: float,
  h: var openArray[float],
  sigma: float, tol: float = 0.5
): void =
  # Compute q taps of impulse response
  # q is expected to be 4 from the given code
  let strideAbs = abs(stride)
  let strideN = strideAbs * N
  var off = if stride < 0: strideN + stride else: 0

  # Compute first q taps of impulse response h
  for n in 0..<q:
    h[n] = if n <= p: b[n] else: 0
    for m in 1..min(q, n):
      h[n] -= a[m]*h[n - m]

  # dest_m = sum_{n=1}^m h_{m-n} src_n
  for m in 0..<q:
    dest[m] = 0
    for n in 1..m:
      let i = off + stride*n
      if i >= 0 and i < strideN:
        dest[m] += h[m - n]*src[i]

  # Now pad with zero (or repeated) values from src boundary
  # The code in Python accumulates with infinite zero padding
  # Actually it reuses 'cur' = src[off]
  var cur = src[off]
  var sum_val = sum_value
  let max_iter = int(ceil(sigma * 10.0))
  for n in 0 ..< max_iter:
    for m in 0 ..< q:
      dest[m] += h[m] * cur

    sum_val -= abs(h[0])
    if sum_val <= tol:
      break

    # Compute next impulse tap h_{n+q}
    var next_h = if n+q <= p: b[n+q] else: 0.0
    for m in 1..q:
      next_h -= a[m] * h[q - m]

    # Shift h array
    for m in 0 ..< q-1:
      h[m] = h[m+1]
    h[q-1] = next_h


proc deriche_conv_1d*(
  c: DericheConfig,
  src: openArray[float],
  N: int, stride: int,
  y_causal: var openArray[float],
  y_anticausal: var openArray[float],
  h: var openArray[float],
  d: var openArray[float]
) =
  let stride2 = stride * 2
  let stride3 = stride * 3
  let stride4 = stride * 4
  let strideN = stride * N

  # Initialize causal filter
  deriche_init_zero_pad(
    y_causal, src, N, stride,
    c.b_causal, 3,
    c.a, 4,
    c.sum_causal, h,
    c.sigma
  )

  # Filter interior using a 4th order filter
  var i = stride4
  for n in 4 ..< N:
    y_causal[n] =
      c.b_causal[0] * src[i] +
      c.b_causal[1] * src[i - stride] +
      c.b_causal[2] * src[i - stride2] +
      c.b_causal[3] * src[i - stride3] -
      c.a[1] * y_causal[n - 1] -
      c.a[2] * y_causal[n - 2] -
      c.a[3] * y_causal[n - 3] -
      c.a[4] * y_causal[n - 4]
    i += stride

  # Initialize anticausal filter
  deriche_init_zero_pad(
    y_anticausal, src, N, -stride,
    c.b_anticausal, 4,
    c.a, 4,
    c.sum_anticausal, h,
    c.sigma
  )

  i = strideN - stride*5
  for n in 4 ..< N:
    y_anticausal[n] =
      c.b_anticausal[1] * src[i + stride] +
      c.b_anticausal[2] * src[i + stride2] +
      c.b_anticausal[3] * src[i + stride3] +
      c.b_anticausal[4] * src[i + stride4] -
      c.a[1] * y_anticausal[n - 1] -
      c.a[2] * y_anticausal[n - 2] -
      c.a[3] * y_anticausal[n - 3] -
      c.a[4] * y_anticausal[n - 4]
    i -= stride

  # Combine causal and anticausal
  i = 0
  if c.negative:
    for n in 0..<N:
      d[i] = y_causal[n] + y_anticausal[N - n - 1]
      i += stride
  else:
    for n in 0..<N:
      d[i] = max(0.0, y_causal[n] + y_anticausal[N - n - 1])
      i += stride


proc density_1d*(
  data: var openArray[float],
  extent: tuple[lo: float, hi: float] = (0.0, 0.0),
  weight: seq[float] = @[],
  bandwidth: float = NaN,
  adjust: float = 1.0,
  pad: int = 3,
  bins: int = 512
): tuple[density: seq[float], lo: float, hi: float] =

  var w: seq[float]
  if weight.len == 0:
    # If no weight provided, use uniform weights
    w = newSeq[float](data.len)
    let uw = 1.0 / float(data.len)
    for i in 0 ..< data.len:
      w[i] = uw
  else:
    w = weight

  var bw = bandwidth
  if isNaN(bw):
    bw = adjust * nrd(data)

  var lo, hi: float
  if extent.lo == 0.0 and extent.hi == 0.0:
    let (elo, ehi) = density_extent(data, float(pad) * bw)
    lo = elo
    hi = ehi
  else:
    (lo, hi) = extent

  let grid = bin_1d(data, w, lo, hi, bins)
  let delta = (hi - lo) / float(bins - 1)
  let neg = has_negative(grid)

  let config = deriche_config(bw / delta, neg)

  # Prepare arrays for deriche_conv1d if needed:
  var y_causal = newSeq[float](bins)
  var y_anticausal = newSeq[float](bins)
  var h = newSeq[float](5)
  var d = newSeq[float](bins)

  deriche_conv_1d(config, grid, bins, 1, y_causal, y_anticausal, h, d)
  return (d, lo, hi)
