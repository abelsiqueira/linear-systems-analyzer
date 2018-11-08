include("analyzer.jl")

function test(A, PC, PN, PCT, PNT)
    C, N, CT, NT = findSubSpaces(A)

    @testset "Matrices matching" begin
        @test isapprox(norm(PC, 1), norm(C, 1), atol=ℯ^(-12))
        @test isapprox(norm(PN, 1), norm(N, 1), atol=ℯ^(-12))
        @test isapprox(norm(PCT, 1), norm(CT, 1), atol=ℯ^(-12))
        @test isapprox(norm(PNT, 1), norm(NT, 1), atol=ℯ^(-12))
    end
end

test([1 3 4; 2 2 4; 3 1 4], [1.0 3.0; 2.0 2.0; 3.0 1.0], [-1.0; -1.0; 1.0],
        [14.0 10.0; 10.0 14.0; 24.0 24.0], [0.0; 0.0; 0.0])
test([1 2; 2 1; 3 3], [1.0 2.0; 2.0 1.0; 3.0 3.0], [0.0, 0.0], [14.0 13.0; 13.0 14.0],
        [0.0, 0.0])
