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

println("-----------------------------------------------------------------")
# Sobra baixo -
test(Float64[1 2; 2 1; 3 3], [1 2; 2 1; 3 3], [0; 0], [14 13; 13 14], [-1, -1, 1])
# Sobra direita -
test(Float64[1 2 3; 2 1 3], [1 2; 2 1], [-1; -1; 1], [5 4; 4 5; 9 9], [0; 0])
# Sobra baixo e direita -
test(Float64[1 3 4; 2 2 4; 3 1 4], [1 3; 2 2; 3 1], [-1; -1; 1], [14 10; 10 14; 24 24],
        [1; -2; 1])
# Sobra mais baixo -
test(Float64[1 1 2; 3 5 8; 5 3 8; 4 4 8], [1 1; 3 5; 5 3; 4 4], [-1; -1; 1],
        [51 47; 47 51; 98 98], [-8 -4; 1 0; 1 0; 0 1])
# Sobra mais direita * faltam testar Im(A^T) -
test(Float64[3 -1 1 2; -1 3 5 2; 1 10 21 11], [3 -1; -1 3; 1 10], [-1 -1; -2 -1; 1 0; 0 1],
        [11 4; 4 110; 19 224; 15 114], [-1.625; -3.875; 1])
# Sobra mais baixo e troca linhas
test(Float64[1 3.5 4.5 2.5; 0 0 0 0; 3.5 1 4.5 -2.5; 2 3 5 1], [1 3.5; 0 0; 3.5 1; 2 3],
        [-1 1; -1 -1; 1 0; 0 1], [17.25 13; 13 22.25; 30.25 35.25; -4.25 9.25],
        [0 -0.755556; 1 0; 0 -0.355556; 0 1])
# Completo
test(Float64[1 2; 2 1], [1 2; 2 1], [0, 0], [5 4; 4 5], [0; 0])

println("-----------------------------------------------------------------")
