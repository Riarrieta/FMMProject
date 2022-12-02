using FMMProject
using FMMProject: Point2D,compute_bounding_box
using Plots

N = 1000
L = 23
points = rand(Point2D,N) .* L
box = compute_bounding_box(points)

## plot
scatter([p[1] for p in points],[p[2] for p in points])
scatter!([box.lc[1]],[box.lc[2]])
scatter!([box.uc[1]],[box.uc[2]])