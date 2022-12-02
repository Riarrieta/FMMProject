using Test
using FMMProject
import FMMProject: Point2D,
                   Point3D,
                   Box,
                   center,
                   halfside,
                   corners

@testset "Box" begin
    ce = Point2D(0.5,0.5)
    a  = 0.5
    box = Box(ce,a)
    @test center(box) == ce
    @test halfside(box) == a
    @test corners(box) == (Point2D(0,0),Point2D(1,1))

    ce = Point3D(0.5,0.5,0.5)
    box = Box(ce,a)
    @test center(box) == ce
    @test halfside(box) == a
    @test corners(box) == (Point3D(0,0,0),Point3D(1,1,1))

end