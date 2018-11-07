include("analyzer.jl")

function test(A, PC, PN, PCT, PNT)
    #=
    U = reduceMatrix(A)
    C = findColumnSpace(A)
    N = findNullspace(A)
    @test isapprox(norm(PC, 1),  norm(C, 1), atol=10^-12)
    @test isapprox(norm(PN, 1), norm(N, 1), atol=10^-12)

    UT = reduceMatrixT(A)
    CT = findColumnSpaceT(A)
    NT = findNullspaceT(A)
    @test isapprox(norm(PCT, 1),  norm(CT, 1), atol=10^-12)
    @test isapprox(norm(PNT, 1), norm(NT, 1), atol=10^-12)
    =#

    C, N, CT, NT = findSubSpaces(A)

    @testset "Matrix match" begin
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


#=
function testAll()

    A = [1 2; 2 1; 3 3]
    test(A, [1 2; 2 1; 3 3], [0; 0], [1 2; 2 1], [-1; -1; 1])

    A = [1 2 3; 2 1 3]
    test(A, [2 1; 1 2], [-1; -1; 1], [1 2; 2 1; 3 3], [0; 0])

    A = [1 2; 2 1]
    test(A, [1 2; 2 1], [0; 0], [1 2; 2 1], [0; 0])

    A = [1 2 3; 2 1 3; -5 2 -3]
    test(A, [1 2; 2 1; -5 2], [-1; -1; 1], [1 2; 2 1; 3 3], [-3; 4; 1])

end
=#



#A = [1 2 3; 2 1 3; -5 2 -3]



#N = findNullspace(A)
#NT = findNullspaceT(A)
#writedlm("matrix_N.txt", N, " ")
#writedlm("matrix_NT.txt", NT, " ")
