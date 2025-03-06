using GraphPartitioning
using Test

@testset verbose = true "GraphPartitioning Tests" begin

    A, coords = GraphPartitioning.build_adjacency("network")

    @testset "Coordinate Partitioning" begin
        res = GraphPartitioning.part_coordinate(A, coords)
        @test length(res) == size(A, 1)
        @test all(x -> !isnan(x), res)
    end

    @testset "Testing Inertial Partitioning" begin
        res = GraphPartitioning.part_inertial(A, coords) 
        @test length(res) == size(A, 1) 
        @test all(x -> !isnan(x), res)  
    end

    @testset "Testing Spectral Partitioning" begin
        res = GraphPartitioning.part_spectral(A)
        @test length(res) == size(A, 1)
        @test all(x -> !isnan(x), res) 
    end

    @testset "Testing Metis Partitioning" begin
        res = GraphPartitioning.part_metis(A,2,:KWAY)
        @test length(res) == size(A, 1)
        @test all(x -> !isnan(x), res) 
    end

    @testset "Testing Recursive Bisection: Coordinate" begin
        res = GraphPartitioning.recursive_bisection(GraphPartitioning.part_coordinate, 3, A, coords)
        @test length(res) == size(A, 1)
        @test all(x -> !isnan(x), res) 
    end

    @testset "Testing Recursive Bisection: Inertial" begin
        res = GraphPartitioning.recursive_bisection(GraphPartitioning.part_inertial, 3, A, coords)
        @test length(res) == size(A, 1)
        @test all(x -> !isnan(x), res) 
    end

    @testset "Testing Recursive Bisection: Spectral" begin
        res = GraphPartitioning.recursive_bisection(GraphPartitioning.part_spectral, 3, A)
        @test length(res) == size(A, 1)
        @test all(x -> !isnan(x), res) 
    end

end