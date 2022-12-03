using Test
using Random
using FMMProject
import FMMProject: Point2D,
                   Point3D,
                   Point,
                   Box,
                   center,
                   halfside,
                   corners,
                   dist,
                   compute_bounding_box,
                   isinbox,
                   split_points,
                   split_box,
                   minus1exp

Random.seed!(1)
@testset "Utils" begin

    for p in -5:5
        @test minus1exp(p) == (-1)^p
    end

    @testset "Box" begin
        for N in [2,3]
            ce = Point{N}(ntuple(i->0.5,N))
            a  = 0.5
            box = Box(ce,a)
            @test center(box) == ce
            @test halfside(box) == a
            @test corners(box) == (Point{N}(ntuple(i->0,N)),Point{N}(ntuple(i->1,N)))
    
            d = rand()
            dvec = ntuple(i->i*d,N)
            box2 = Box(ce .+ dvec,a)
            @test dist(box,box2) ≈ N*d
    
            M = 26
            L = 23
            points = rand(Point{N},M) .* L
            box = compute_bounding_box(points)
            @test L/2 <  box.a < L
            for p in points
                @test isinbox(p,box)
            end
    
            children_points = split_points(points,box)
            children_boxes = split_box(box)
            for (plist,b) in zip(children_points,children_boxes)
                for p in plist
                    @test isinbox(p,b)
                end
            end
        end
    end
end
