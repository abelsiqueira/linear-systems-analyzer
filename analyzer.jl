using DelimitedFiles

function findDimensions(M)
    dimens = size(M)
    m = dimens[1]
    n = dimens[2]

    return m, n
end

function printMatrix(fileName, M)
    writedlm(fileName, M)
end

function exchangeRows(M, i1, i2)
    aux = M[i1, :]
    M[i1, :] = M[i2, :]
    M[i2, :] = aux

    return M
end

function reduceMatrix(M)
    m, n = findDimensions(M)

    row = 1
    column = 1
    while (row <= m) && (column <= n)
        for cont = (row + 1):m
            if M[row, column] == 0
                if M[cont, column] != 0
                    M = exchangeRows(M, row, cont)
                end
                continue
            end
            M[cont, :] -= (M[cont, column] / M[row, column]) * M[row, :]
        end
        if M[row, column] == 0
            column += 1
            continue
        end
        row += 1
        column += 1
    end

    return  M
end

function isPivot(M, row, column)
    m, n = findDimensions(M)

    if M[row, column] == 0
        return false
    end
    for i = (row + 1):m
        if M[i, column] != 0
            return false
        end
    end
    return true
end

function pivotsMarker(M)
    m, n = findDimensions(M)

    pivots = zeros(Float64, n, 1)
    row = 1
    column = 1
    first = true
    while (row <= m) && (column <= n)
        if isPivot(M, row, column)
            pivots[column] = 1
            row += 1
            column += 1
            continue
        end
        column += 1
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
    dimens = size(pivots)

    io = open("bases_C-U.txt", "w+")
    for i = 1:dimens[1]
        if pivots[i] == 1
            print(io, M[:, i], " ")
        end
    end
    close(io)
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
    M = reduceMatrix(M)
    m, n = findDimensions(M)
    rank = findRank(M)
    pivots = pivotsMarker(M)
    x = zeros(Float64, n, 1)

    io = open("bases_N-A.txt", "w+")
    if rank == n
        print(io, x, " ")
        close(io)
        return
    end

    for i = 1:n
        if pivots[i] == 0
            x[i] = 1

            cont = 0
            j = n
            while j >= 1
                if pivots[j] == 1
                    x[j] = (-1 * product(M[rank - cont, :], x, j)) / M[rank - cont, j]
                    cont += 1
                end

                j -= 1
            end
            print(io, x, " ")
        end
    end
    close(io)
end

A = readdlm("matrix_A.txt", ' ', Float64)
findNullspace(A)
U = reduceMatrix(A)
writedlm("matrix_U.txt", U, " ")
findColumnSpace(U)
