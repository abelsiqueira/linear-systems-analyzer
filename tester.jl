include("analyzer.jl")

function test(A, PC, PN, PCT, PNT)
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
end

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

A = [1 2 3; 2 1 3; -5 2 -3]

N = findNullspace(A)
NT = findNullspaceT(A)
writedlm("matrix_N.txt", N, " ")
writedlm("matrix_NT.txt", NT, " ")

testAll()
