using DelimitedFiles

#function printMatrix(fileName, M)
#    writedlm(fileName, M)
#end

function exchangeRows(M, i1, i2)
    aux = M[i1, :]
    M[i1, :] = M[i2, :]
    M[i2, :] = aux

    return M
end

# Return the index of the first non null element finded in or under row i in
# column j. If there isn't pivot, return -1
function pivotIndex(M, i, j)
    m, n = size(M)

    for cont = i:m
        if M[cont, j] != 0
            return cont
        end
        return -1
    end
end

function reduceMatrix(M)
    m, n = size(M)

    i = j = 1
    while (i <= m) && (j <= n)
        row = pivotIndex(M, i, j)
        if row != -1
            if row != i
                exchangeRows(M, row, i)
            end
            for cont = (i + 1):m
                M[cont, :] -= (M[cont, j] / M[row, j]) * M[row, :]
            end
            i += 1
            j += 1
        else
            j += 1
        end
    end

    return M
end

function pivotsMarker(M)
    m, n = size(M)

    pivots = zeros(n)
    i = j = 1
    while (i <= m) && (j <= n)
        if pivotIndex(M, i, j) != -1
            pivots[j] = 1
            i += 1
            j += 1
            continue
        end
        j += 1
    end

    return pivots
end

function findRank(M)
    M = reduceMatrix(M)
    pivots = pivotsMarker(M)
    dimens = size(pivots)

    rank = 0
    for i = 1:dimens[1]
        if pivots[i] == 1
            rank += 1
        end
    end

    return rank
end

function findColumnSpace(M)
    pivots = pivotsMarker(M)
    rank = findRank(M)
    m, n = size(M)

    C = zeros(Float64, m, rank)
    for i = 1:n
        if pivots[i] == 1
            C[:, i] = M[:, i]
        end
    end

    return C
end

# Return value of ((inner product between v1 an v2) - (v1[index] * v2[index]))
function product(v1, v2, index)
    dimens = size(v1)

    sum = 0
    for i = 1:dimens[1]
        if i != index
            sum += v1[i] * v2[i]
        end
    end

    return sum
end

function findNullspace(M)
    U = reduceMatrix(M)
    m, n = size(U)
    rank = findRank(U)
    pivots = pivotsMarker(U)
    N = zeros(n, n - rank)

    if rank != n
        column = 1 # Column where will be put the x vector finded in each iteration

        for i = 1:n
            if pivots[i] == 0
                x = zeros(Float64, n, 1)
                x[i] = 1

                cont = 0
                j = n
                while j >= 1
                    if pivots[j] == 1
                        x[j] = (-1 * product(U[rank - cont, :], x, j)) / U[rank - cont, j]
                        cont += 1
                    end

                    j -= 1
                end
                N[:, column] = x
                column += 1
            end
        end
    end

    return N
end

function testAll(U, C, N)
    return @testset "tests" begin
               @test [1 2 3 0 -3 -3] == U
               @test [1 2 0 -3] == C
               @test [-1 -1 1] == N
           end;
end


cd("/home/leonardo/githubProjects/linear-systems-analyzer")
A = readdlm("matrix_A.txt", ' ', Float64)
U = reduceMatrix(A)
C = findColumnSpace(A)
N = findNullspace(A)
#writedlm("matrix_U.txt", U, " ")
writedlm("matrix_U.txt", U, " ")
writedlm("matrix_C.txt", C, " ")
writedlm("matrix_N.txt", N, " ")
#io =  open("resultado.txt", "w+")
#testAll(U, C, N)
#close(io)
