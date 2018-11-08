using DelimitedFiles, Test, LinearAlgebra

function exchangeRows(M, i1, i2)
    aux = M[i1, :]
    M[i1, :] = M[i2, :]
    M[i2, :] = aux
end

function pivotIndex(M, i, j)
    m, n = size(M)

    for cont = i:m
        if M[cont, j] != 0
            return cont
        end
    end
    return 0
end

function reduceMatrix(A)
    m, n = size(A)
    U = copy(A)
    E = Matrix{Float64}(I, m, m)
    pivots = zeros(Int8, n)

    i = j = 1
    while (i <= m) && (j <= n)
        row = pivotIndex(U, i, j)
        if row != 0
            pivots[j] = 1
            if row != i
                exchangeRows(U, row, i)
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

    return U, E, pivots, rank
end

function findColumnSpace(A, pivots, rank)
    m, n = size(A)

    if rank == 0
        return zeros(Float64, m)
    end

    C = zeros(Float64, m, rank)
    cont = 1
    for i = 1:n
        if pivots[i] == 1
            C[:, cont] = A[:, i]
            cont += 1
        end
    end

    return C
end

function findColumnSpaceT(A, C)
    CT = A' * C

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

function findNullspace(U, pivots, rank)
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

function findNullSpaceT(A, E, rank)
    m, n = size(A)
    NT = zeros(Float64, n, (m - rank))

    for cont = 1:(m - rank)
        NT[:, cont] = E[rank + cont, :]
    end

    return NT
end

function findSubSpaces(A)
    U, E, pivots, rank = reduceMatrix(A)
    C = findColumnSpace(A, pivots, rank)
    N = findNullspace(U, pivots, rank)
    CT = findColumnSpaceT(A, C)
    NT = findNullspace(A, E, rank)

    return C, N, CT, NT
end
