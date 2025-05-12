using GraphLab
using Test

@testset verbose = true "GraphLab Tests" begin

    A, coords = GraphLab.build_adjacency("network")

    @testset "Coordinate Partitioning" begin
        res = GraphLab.part_coordinate(A, coords)
        @test length(res) == size(A, 1)
        @test all(x -> !isnan(x), res)
    end

    @testset "Testing Inertial Partitioning" begin
        res = GraphLab.part_inertial(A, coords) 
        @test length(res) == size(A, 1) 
        @test all(x -> !isnan(x), res)  
    end

    @testset "Testing Spectral Partitioning" begin
        res = GraphLab.part_spectral(A)
        @test length(res) == size(A, 1)
        @test all(x -> !isnan(x), res) 
    end

    @testset "Testing Metis Partitioning" begin
        res = GraphLab.part_metis(A,2,:KWAY)
        @test length(res) == size(A, 1)
        @test all(x -> !isnan(x), res) 
    end

    @testset "Testing Random Sphere Partitioning" begin
        res = GraphLab.part_randsphere(A, coords)
        @test length(res) == size(A, 1)
        @test all(x -> !isnan(x), res)
    end

    @testset "Testing Geometric Spectral Partitioning" begin
        res = GraphLab.part_geospectral(A)
        @test length(res) == size(A, 1)
        @test all(x -> !isnan(x), res)
    end

    @testset "Testing Space-filling Curve Partitioning" begin
        res = GraphLab. part_sfc(A, coords)
        @test length(res) == size(A, 1)
        @test all(x -> !isnan(x), res)
    end

    @testset "Testing Adaptive Space-filling Curve Partitioning" begin
        res = GraphLab. part_adaptive_sfc(A, coords)
        @test length(res) == size(A, 1)
        @test all(x -> !isnan(x), res)
    end


    @testset "Testing Recursive Bisection: Coordinate" begin
        res = GraphLab.recursive_bisection(GraphLab.part_coordinate, 3, A, coords)
        @test length(res) == size(A, 1)
        @test all(x -> !isnan(x), res) 
    end

    @testset "Testing Recursive Bisection: Inertial" begin
        res = GraphLab.recursive_bisection(GraphLab.part_inertial, 3, A, coords)
        @test length(res) == size(A, 1)
        @test all(x -> !isnan(x), res) 
    end

    @testset "Testing Recursive Bisection: Spectral" begin
        res = GraphLab.recursive_bisection(GraphLab.part_spectral, 3, A)
        @test length(res) == size(A, 1)
        @test all(x -> !isnan(x), res) 
    end

    @testset "Testing Nested Dissection" begin
        res = GraphLab.nested_dissection(A, part_coordinate; coords=coords)
        @test length(res) == size(A, 1)
        @test all(x -> !isnan(x), res) 
        @test all(res .>= 1) && all(rest .<= size(A, 1))
    end
end