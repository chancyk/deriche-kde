import derichekde


block:
  var data = @[1.2, 2.3, 23, 40, 50, 60, 70, 60, 50, 400, 700, 1000]
  let (density, lo, hi) = density_1d(data, extent=(0.0, 1000.0), bandwidth=2.0, bins=512)

  doAssert density.len == 512
  doAssert lo == 0.0
  doAssert hi == 1000.0
