# Importação da biblioteca numpy
import numpy as np

# Definição da função "codeCoord" com dados de entrada das funções "readmesh" e "nodeindex"
def codeCoord(Dx, Dy, NEX, NEY, connect, X, Y, Nnodes):

    # Determinação das coordenadas X e Y dos nós da estrutura - Caso de 8 nós por elementos
    if Nnodes == 8:
        for ix in range(NEX):
            for jy in range(NEY):
                xlocal = []
                ylocal = []
                iele = (ix) * NEY + jy
                xlocal.append((ix) * Dx)
                xlocal.append(xlocal[0] + Dx)
                xlocal.append(xlocal[1])
                xlocal.append(xlocal[0])
                xlocal.append(1 / 2 * (xlocal[0] + xlocal[1]))
                xlocal.append(1 / 2 * (xlocal[1] + xlocal[2]))
                xlocal.append(1 / 2 * (xlocal[2] + xlocal[3]))
                xlocal.append(1 / 2 * (xlocal[0] + xlocal[3]))

                ylocal.append((jy) * Dy)
                ylocal.append(ylocal[0])
                ylocal.append(ylocal[0] + Dy)
                ylocal.append(ylocal[2])
                ylocal.append(ylocal[0])
                ylocal.append(1 / 2 * (ylocal[2] + ylocal[0]))
                ylocal.append(ylocal[2])
                ylocal.append(ylocal[5])

                for ilnode in range(Nnodes):
                    ignode = int(connect[ilnode, iele])
                    X[ignode, 0] = xlocal[ilnode]
                    Y[ignode, 0] = ylocal[ilnode]

    # Determinação das coordenadas X e Y dos nós da estrutura - Caso de 4 nós por elementos
    elif Nnodes == 4:
        for ix in range(NEX):
            for jy in range(NEY):
                xlocal = []
                ylocal = []
                iele = (ix) * NEY + jy
                xlocal.append((ix) * Dx)
                xlocal.append(xlocal[0] + Dx)
                xlocal.append(xlocal[1])
                xlocal.append(xlocal[0])

                ylocal.append((jy) * Dy)
                ylocal.append(ylocal[0])
                ylocal.append(ylocal[0] + Dy)
                ylocal.append(ylocal[2])

                for ilnode in range(Nnodes):
                    ignode = int(connect[ilnode, iele])
                    X[ignode, 0] = xlocal[ilnode]
                    Y[ignode, 0] = ylocal[ilnode]

    # Tratamento de erro para eventual entrada de dados incompatível com a análise
    else:
        print("Erro na quantidade de nós por elemento. Só é tolerado 4 ou 8, porém o valor adotado é igual a:", Nnodes)

    return X, Y

# Definição da função "Coord" com dados de entrada das funções "readmesh", "nodeindex" e "codeCoord"
def Coord(nodes, Dx, Dy, NEX, NEY, connect, Nnodes):

    X = np.zeros((nodes, 1))
    Y = np.zeros((nodes, 1))
    [X, Y] = codeCoord(Dx, Dy, NEX, NEY, connect, X, Y, Nnodes)

    # Geração da matriz de coordenadas
    XY = np.zeros((len(X), 2))
    XY[:, 0] = X[:, 0]
    XY[:, 1] = Y[:, 0]

    return X, Y, XY