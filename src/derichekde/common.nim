import std/[math, algorithm]


proc density_extent*(data: openArray[float], pad: float = 0.0): tuple[lo: float, hi: float] =
  return (min(data) - pad, max(data) + pad)


proc stdev(values: openArray[float]): float =
  let n = values.len
  if n < 2:
    return NaN
  var count = 0
  var mean = 0.0
  var sum_val = 0.0
  for i in 0..<n:
    count += 1
    let value = values[i]
    let delta = value - mean
    mean += delta / float(count)
    sum_val += delta*(value - mean)
  return sqrt(sum_val / float(count - 1))


proc quantile(values: openArray[float], p: float): float =
  let n = values.len
  if n == 0:
    return NaN
  if p <= 0 or n < 2:
    return values[0]
  if p >= 1:
    return values[n-1]

  let i0 = (float(n) - 1.0) * p
  let i1 = floor(i0).int
  return values[i1] + (values[i1+1] - values[i1])*(i0 - float(i1))


proc nrd*(data: var openArray[float]): float =
  # Sort the data in-place
  sort(data)
  let sd = stdev(data)
  let q1 = quantile(data, 0.25)
  let q3 = quantile(data, 0.75)
  let n = data.len.float
  let h = (q3 - q1) / 1.34

  # Replicate the pythonic "or" logic:
  # v = min(sd,h) or sd or abs(q1) or 1
  var v = min(sd, h)
  if v == 0.0:
    v = sd
    if v == 0.0:
      v = abs(q1)
      if v == 0.0:
        v = 1.0

  return 1.06 * v * pow(n, -0.2)


proc bin_1d*(
  data: openArray[float], weight: openArray[float],
  lo: float, hi: float, n: int
): seq[float] =
  # Initialize the grid with zeros
  var grid = newSeq[float](n)
  let delta = (hi - lo) / float(n - 1)

  for i in 0 ..< data.len:
    let p = (data[i] - lo) / delta
    let u = int(floor(p))
    let v = u + 1

    if u >= 0 and u < n and v < n:
      grid[u] += (float(v) - p) * weight[i]
      grid[v] += (p - float(u)) * weight[i]
    elif u == -1:
      # v = 0 in this case
      if v < n:
        grid[v] += (p - float(u)) * weight[i]
    elif v == n:
      # u = n-1 in this case
      if u >= 0 and u < n:
        grid[u] += (float(v) - p) * weight[i]

  return grid


proc has_negative*(values: openArray[float]): bool =
  for value in values:
    if value < 0:
      return true
  return false
