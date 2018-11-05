using DelimitedFiles, Test, LinearAlgebra

function exchangeRows(M, i1, i2)
    M2 = copy
    aux = M2[i1, :]
    M2[i1, :] = M2[i2, :]
    M2[i2, :] = aux

    return M2
end
#=
function exchangeColumns(M, j1, j2)
    M2 = copy(M)
    aux = M2[:, j1]
    M2[:, j1] = M2[:, j2]
    M2[:, j2] = aux

    return M2
end
=#
function pivotIndex(M, i, j)
    m, n = size(M)

    for cont = i:m
        if M[cont, j] != 0
            return cont
        end
    end
    return -1
end
#=
function pivotIndexT(M, i, j)
    m, n = size(M)

    for cont = j:n
        if M[i, cont] != 0
            return cont
        end
        return -1
    end
end
=#
function reduceMatrix(M)
    m, n = size(M)
    pivots = zeros(Int8, n)

    U = copy(M)
    E = Matrix{Float64}(I, m, m)

    i = j = 1
    while (i <= m) && (j <= n)
        row = pivotIndex(U, i, j)
        if row != -1
            pivots[j] = 1
            if row != i
                exchangeRows(M, row, i)
                exchangeRows(E, row, i)
            end
            for cont = (i + 1):m
                mult = U[cont, j] / U[row, j]
                U[cont, :] -= mult * U[row, :]
                E[cont, :] -= mult * E[row, :]
            end
            i += 1
            j += 1
        else
            j += 1
        end
    end

    rank = 0
    for i = 1:n
        if pivots[i] == 1
            rank += 1
        end
    end

    return U, pivots, rank, E
end
#=
function reduceMatrixT(M)
    m, n = size(M)
    pivots = zeros(n)

    U = copy(M)
    i = j = 1
    while (i <= m) && (j <= n)
        column = pivotIndex(M, i, j)
        if column != -1
            pivots[column] = 1
            if column != i
                exchangeColumns(M, column, j)
            end
            for cont = (j + 1):n
                U[:, cont] -= (U[i, cont] / U[i,column]) * U[:, column]
            end
            i += 1
            j += 1
        else
            i += 1
        end
    end

    rank = 0
    for i = 1:n
        if pivots[i] == 1
            rank += 1
        end
    end

    return U, pivots, rank
end
=#
function findColumnSpace(M)
    U, pivots, rank = reduceMatrix(M)
    m, n = size(M)

    if rank == 0
        return zeros(Float64, m)
    end

    C = zeros(Float64, m, rank)
    cont = 1
    for i = 1:n
        if pivots[i] == 1
            C[:, cont] = M[:, i]
            cont += 1
        end
    end

    return C
end
#=
function findColumnSpaceT(M)
    U, pivots, rank = reduceMatrixT(M)
    m, n = size(M)

    if rank == 0
        return zeros(Float64, n)
    end

    C = zeros(Float64, n, rank)
    for i = 1:m
        if pivots[i] == 1
            C[:, i] = M[i, :]
        end
    end

    return C
end
=#

function findColumnSpaceT(M)
    C = findColumnSpace(M)
    CT = M' * C

    return CT

end

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
    U, pivots, rank = reduceMatrix(M)
    m, n = size(U)
    N = zeros(Float64, n, n - rank)

    if rank == n
        return zeros(Float64, n)
    end

    column = 1

    for i = 1:n
        if pivots[i] == 0
            x = zeros(Float64, n)
            x[i] = 1

            cont = 0
            j = n
            while j >= 1
                if pivots[j] == 1
                    x[j] = (-product(U[rank - cont, :], x, j)) / U[rank - cont, j]
                    cont += 1
                end

                j -= 1
            end
            N[:, column] = x
            column += 1
        end
    end

    return N
end
#=
function findNullspaceT(M)
    U, pivots = reduceMatrixT(M)
    m, n = size(U)
    rank = m
    N = zeros(m, m - rank)

    if rank == m
        return zeros(Float64, m)
    end

    row = 1

    for i = 1:m
        if pivots[i] == 0
            x = zeros(Float64, m, 1)
            x[i] = 1

            cont = 0
            j = m
            while j >= 1
                if pivots[j] == 1
                    x[j] = (-1 * product(U[:, rank - cont], x, j)) / U[j, rank - cont]
                    cont += 1
                end

                j -= 1
            end
            N[:, row] = x
            row += 1
        end
    end

    return N
end
=#

function findNullSpaceT(m, n, E, rank)
    NT = zeros(Float64, n, (m - rank))

    for cont = 1:(m - rank)
        NT[:, cont] = E[rank + cont, :]
    end

    return NT
end

A = [1 3 4; 2 2 4; 3 1 4]
U, pivots, r, E = reduceMatrix(A)
C = findColumnSpace(A)
N = findNullspace(A)
CT = findColumnSpaceT(A)
NT = findNullSpaceT(3, 3, E, r)
println("-------------------")
println("A = ", A)
println("A^T", A')
println("U = ", U)
println("pivots = ", pivots)
println("Im(A) = ", C)
println("Nu(A) = ", N)
println("Im(A^T) = ", CT)
println("NU(A^T) = ", NT)
println("E = ", E)
