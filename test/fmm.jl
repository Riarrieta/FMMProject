using Test
using Random
using LinearAlgebra
using FMMProject
import FMMProject: nlevels,
                   eachlevel,
                   maxlevel,
                   children,
                   parent,
                   neighbor_list,
                   interaction_list,
                   Laplace2D
Random.seed!(1)

@testset "FMM Laplace 2D" begin
    npoints = 53
    L = 4  # max level
    P = 3  # interaction rank
    points = rand(Point2D,npoints)
    copy_points = deepcopy(points)

    fmm = FMMLaplace2D(points,P,L);
    @test L == nlevels(fmm) == maxlevel(fmm)
    points_tree = 0
    for (l,level) in eachlevel(fmm)
        @test length(level) == 4^(l-1)
        for tree in level
            for child in children(tree)
                @test parent(child) === tree
            end
            @test length(tree.Tofo) == length(children(tree))
            @test length(tree.Tifo) == length(interaction_list(tree))
            if l < maxlevel(fmm) # if not leaf
                @test !isempty(children(tree))
                @test isempty(neighbor_list(tree))
                @test isempty(tree.points)
                @test isempty(tree.points_indices)
                @test size(tree.Tofs) == (0,0)
                @test size(tree.Ttfi) == (0,0)
                if l ≤ 2  # if root or second level
                    @test isempty(interaction_list(tree))
                else
                    @test !isempty(interaction_list(tree))
                end
            else  # if leaf
                @test isempty(children(tree))
                @test !isempty(interaction_list(tree))
                @test !isempty(neighbor_list(tree))
                @test length(tree.points) == length(tree.points_indices)
                @test size(tree.Tofs) == (P,length(tree.points))
                @test size(tree.Ttfi) == (length(tree.points),P)
                @test tree ∉ neighbor_list(tree)
                @test tree ∉ interaction_list(tree)
                for (i,p) in zip(tree.points_indices,tree.points)
                    @test fmm.points[i] == p
                end
                points_tree += length(tree.points)
            end
        end
    end
    @test points_tree == npoints
    
    # Forward map
    K = Laplace2D
    charges = ones(npoints)
    copy_charges = deepcopy(charges)
    pot = K(points)*charges
    pot_approx = fmm*charges

    @test pot_approx == fmm*charges
    @test points == copy_points
    @test copy_charges == charges

    @test norm(pot_approx-pot,Inf)/norm(pot,Inf) < 1  # change error
end
