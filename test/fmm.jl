using Test
using Random
using FMMProject
import FMMProject: nlevels,
                   eachlevel,
                   children,
                   parent
Random.seed!(1)

@testset "FMM Laplace 2D" begin
    npoints = 53
    L = 3  # max level
    P = 3  # interaction rank
    points = rand(Point2D,npoints)

    fmm = FMMLaplace2D(points,P,L);
    @test L == nlevels(fmm)
    for (l,level) in eachlevel(fmm)
        @test length(level) == 4^(l-1)
        for tree in level
            for child in children(tree)
                @test parent(child) === tree
            end
        end
    end

end
