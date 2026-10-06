using Test
using CubedSphere
using CubedSphere.SphericalGeometry
using Documenter
using StaticArrays: SVector

# Bytes allocated by `f(args...)`, after compiling it
allocations(f, args...) = (f(args...); @allocated f(args...))

B_Rancic_correct = [
    0.00000000000000,
    0.67698819751739,
    0.11847293456554,
    0.05317178134668,
    0.02965810434052,
    0.01912447304028,
    0.01342565621117,
    0.00998873323180,
    0.00774868996406,
    0.00620346979888,
    0.00509010874883,
    0.00425981184328,
    0.00362308956077,
    0.00312341468940,
    0.00272360948942,
    0.00239838086555,
    0.00213001905118,
    0.00190581316131,
    0.00171644156404,
    0.00155493768255,
    0.00141600715207,
    0.00129556597754,
    0.00119042140226,
    0.00109804711790,
    0.00101642216628,
    0.00094391366522,
    0.00087919021224,
    0.00082115710311,
    0.00076890728775,
    0.00072168382969,
    0.00067885087750
]

@testset "CubedSphere" begin
    @testset "Rančić et al. (1996) Taylor coefficients" begin
        for k in eachindex(B_Rancic_correct)
            @test CubedSphere.B_Rancic[k] ≈ B_Rancic_correct[k]
        end
    end

    @testset "Rančić et al. (1996) inverse mapping" begin
        ξ = collect(0:0.1:1)
        η = collect(0:0.1:1)

        ξ′ = 0ξ
        η′ = 0η

        # transform from cube -> sphere -> cube
        for j in eachindex(η), i in eachindex(ξ)
            ξ′[i], η′[j] = conformal_cubed_sphere_inverse_mapping(conformal_cubed_sphere_mapping(ξ[i], η[j])...)
        end

        @test ξ ≈ ξ′ && η ≈ η′
    end

    @test_throws ArgumentError conformal_cubed_sphere_mapping(2, 0.5)
    @test_throws ArgumentError conformal_cubed_sphere_mapping(0.5, -2)

    @testset "Spherical geometry" begin
        a₁, a₂, a₃, a₄ = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (1 / √2, 1 / √2, 0.0)

        @testset "Same areas for tuples, Vectors and SVectors" begin
            for vertex in (identity, collect, SVector)
                @test spherical_area_triangle(vertex.((a₁, a₂, a₃))...) ≈ π / 2
                @test spherical_area_quadrilateral(vertex.((a₁, a₂, a₃, a₄))...) ≈ π / 2
            end
        end

        @testset "Areas don't allocate" begin
            @test allocations(spherical_area_triangle, a₁, a₂, a₃) == 0
            @test allocations(spherical_area_quadrilateral, a₁, a₂, a₃, a₄) == 0
            @test allocations(spherical_area_quadrilateral, SVector.((a₁, a₂, a₃, a₄))...) == 0
        end

        @testset "Quadrilateral vertices are SVectors" begin
            X = rand(3, 3)
            Y = rand(3, 3)
            Z = rand(3, 3)
            vertices = spherical_quadrilateral_vertices(X, Y, Z, 1, 2)
            @test vertices isa NTuple{4, SVector{3, Float64}}
            @test vertices[3] == SVector(X[2, 3], Y[2, 3], Z[2, 3])
            @test allocations(spherical_quadrilateral_vertices, X, Y, Z, 1, 2) == 0
        end
    end
end

@time @testset "Doctests" begin
    doctest(CubedSphere)
end
