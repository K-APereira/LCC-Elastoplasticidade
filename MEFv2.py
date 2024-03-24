import numpy as np

# Definição da função "constMtrx" com dados de entrada da função "readmesh"
def constMtrx(planestress, E, v):

    # Geração da matriz constitutiva nula
    C = np.matrix(np.zeros((3, 3)))

    # Definição da matriz constitutiva para o Estado Plano de Tensões
    if planestress == 1:
        const = (E / (1 - (v ** 2)))
        C[0, 0] = 1 * const
        C[1, 0] = v * const
        C[2, 0] = 0 * const
        C[0, 1] = v * const
        C[1, 1] = 1 * const
        C[2, 1] = 0 * const
        C[0, 2] = 0 * const
        C[1, 2] = 0 * const
        C[2, 2] = ((1 - v) / 2) * const

    # Definição da matriz constitutiva para o Estado Plano de Deformações
    if planestress == 2:
        const = (E * (1 - v) / ((1 + v) * (1 - 2 * v)))
        C[0, 0] = 1 * const
        C[1, 0] = (v/(1-v)) * const
        C[2, 0] = 0 * const
        C[0, 1] = (v/(1-v)) * const
        C[1, 1] = 1 * const
        C[2, 1] = 0 * const
        C[0, 2] = 0 * const
        C[1, 2] = 0 * const
        C[2, 2] = ((1 - 2 * v)/(2 * (1 - v))) * const

    # Saída de dados da função "nodeindex"
    return C


# Definição da função "PtsGauss1d" com dado de entrada NGP da função "readmesh"
def PtsGauss1d(ordem):

    # Geração do vetor nulo de pontos notáveis básicos do Método de Gauss
    csi = np.matrix(np.zeros((ordem, 1)))

    # Geração do vetor nulo de pesos básicos do Método de Gauss
    w = np.matrix(np.zeros((ordem, 1)))

    # Definição dos pontos notáveis e pesos básicos do método para NGP = 1
    if ordem == 1:
        csi[0] = 0
        w[0] = 2

    # Definição dos pontos notáveis e pesos básicos do método para NGP = 2
    if ordem == 2:
        csi[0] = -1 / np.sqrt(3)
        csi[1] = -csi[0]

        w[0] = 1
        w[1] = 1

    # Definição dos pontos notáveis e pesos básicos do método para NGP = 3
    if ordem == 3:
        csi[0] = -0.774596669241483
        csi[1] = 0
        csi[2] = -csi[0]

        w[0] = 5 / 9
        w[1] = 8 / 9
        w[2] = w[2]

    # Definição dos pontos notáveis e pesos básicos do método para NGP = 4
    if ordem == 4:
        csi[0] = -0.861136311594053
        csi[1] = -0.339981043584856
        csi[2] = 0.339981043584856
        csi[3] = 0.861136311594053

        w[0] = 0.347854845137454
        w[1] = 0.652145154862546
        w[2] = 0.652145154862546
        w[3] = 0.347854845137454

    # Saída de dados da função "PtsGauss1d"
    return csi, w


# Definição da função "DerivFuncFormaNat" com dados de entrada das funções "readmesh" e "PtsGauss1d"
def DerivFuncFormaNat(csi, eta, Nnodes):

    # Geração de matriz linha nula contendo as derivadas parciais das funções de forma em coordenadas naturais
    dNdcsi = np.matrix(np.zeros((1, Nnodes)))
    dNdeta = np.matrix(np.zeros((1, Nnodes)))

    # Definição da matriz linha de derivadas parciais das funções de forma para o caso de elemento com 8 nós
    if Nnodes == 8:
        dPet = 1 + eta
        dMet = 1 - eta
        dPxs = 1 + csi
        dMxs = 1 - csi

        dNdcsi[0, 4] = -csi * dMet
        dNdcsi[0, 5] = 0.50 * dMet * dPet
        dNdcsi[0, 6] = -csi * dPet
        dNdcsi[0, 7] = -0.50 * dMet * dPet

        dNdcsi[0, 0] = -0.25 * dMet - 0.5 * (dNdcsi[0, 4] + dNdcsi[0, 7])
        dNdcsi[0, 1] = 0.25 * dMet - 0.5 * (dNdcsi[0, 4] + dNdcsi[0, 5])
        dNdcsi[0, 2] = 0.25 * dPet - 0.5 * (dNdcsi[0, 5] + dNdcsi[0, 6])
        dNdcsi[0, 3] = -0.25 * dPet - 0.5 * (dNdcsi[0, 6] + dNdcsi[0, 7])

        dNdeta[0, 4] = -0.50 * dMxs * dPxs
        dNdeta[0, 5] = -eta * dPxs
        dNdeta[0, 6] = 0.50 * dMxs * dPxs
        dNdeta[0, 7] = -eta * dMxs

        dNdeta[0, 0] = -0.25 * dMxs - 0.5 * (dNdeta[0, 4] + dNdeta[0, 7])
        dNdeta[0, 1] = -0.25 * dPxs - 0.5 * (dNdeta[0, 4] + dNdeta[0, 5])
        dNdeta[0, 2] = 0.25 * dPxs - 0.5 * (dNdeta[0, 5] + dNdeta[0, 6])
        dNdeta[0, 3] = 0.25 * dMxs - 0.5 * (dNdeta[0, 6] + dNdeta[0, 7])

    # Definição da matriz linha de derivadas parciais das funções de forma para o caso de elemento com 4 nós
    elif Nnodes == 4:
        dPet = 1 + eta
        dMet = 1 - eta
        dPxs = 1 + csi
        dMxs = 1 - csi

        dNdcsi[0, 0] = -0.25 * dMet
        dNdcsi[0, 1] = 0.25 * dMet
        dNdcsi[0, 2] = 0.25 * dPet
        dNdcsi[0, 3] = -0.25 * dPet

        dNdeta[0, 0] = -0.25 * dMxs
        dNdeta[0, 1] = -0.25 * dPxs
        dNdeta[0, 2] = 0.25 * dPxs
        dNdeta[0, 3] = 0.25 * dMxs

    # Tratamento de erro para eventual entrada de dados incompatível com a análise
    else:
        print("Erro na variável Nnodes, pois a quantidade nós atribuída é:", Nnodes, "\n", "Só pode ser atribuído valor igual a 4 ou 8")

    # Saída de dados da função "DerivFuncFormaNat"
    return dNdcsi, dNdeta


# Definição da função "MatDefDesloc" com dados de entrada das funções "readmesh" e "DerivFuncFormaNat"
def MatDefDesloc(dNdx, dNdy, Nnodes):

    # Definição da função B que relaciona deformações e deslocamentos nodais para o caso de elemento com 8 nós
    if Nnodes == 8:
        B = np.matrix([[dNdx[0, 0], 0, dNdx[0, 1], 0, dNdx[0, 2], 0, dNdx[0, 3], 0, dNdx[0, 4], 0, dNdx[0, 5], 0, dNdx[0, 6], 0, dNdx[0, 7], 0],
                       [0, dNdy[0, 0], 0, dNdy[0, 1], 0, dNdy[0, 2], 0, dNdy[0, 3], 0, dNdy[0, 4], 0, dNdy[0, 5], 0, dNdy[0, 6], 0, dNdy[0, 7]],
         [dNdy[0, 0], dNdx[0, 0], dNdy[0, 1], dNdx[0, 1], dNdy[0, 2], dNdx[0, 2], dNdy[0, 3], dNdx[0, 3], dNdy[0, 4], dNdx[0, 4], dNdy[0, 5], dNdx[0, 5], dNdy[0, 6], dNdx[0, 6], dNdy[0, 7], dNdx[0, 7]]])

    # Definição da função B que relaciona deformações e deslocamentos nodais para o caso de elemento com 4 nós
    elif Nnodes == 4:
        B = np.matrix([[dNdx[0, 0], 0, dNdx[0, 1], 0, dNdx[0, 2], 0, dNdx[0, 3], 0],
                       [0, dNdy[0, 0], 0, dNdy[0, 1], 0, dNdy[0, 2], 0, dNdy[0, 3]],
         [dNdy[0, 0], dNdx[0, 0], dNdy[0, 1], dNdx[0, 1], dNdy[0, 2], dNdx[0, 2], dNdy[0, 3], dNdx[0, 3]]])

    # Tratamento de erro para eventual entrada de dados incompatível com a análise
    else:
        print("Erro na variável Nnodes, pois a quantidade nós atribuída é:", Nnodes, "\n", "Só pode ser atribuído valor igual a 4 ou 8")

    # Saída de dados da função "MatDefDesloc"
    return B


# Definição da função "MatJacobiana2D" com dados de entrada das funções "DerivFuncFormaNat" e "Coord"
def MatJacobiana2D(dNdcsi, dNdeta, coords):

    xelem = np.matrix(coords[:, 0])
    yelem = np.matrix(coords[:, 1])

    dxdcsi = dNdcsi * xelem
    dxdeta = dNdeta * xelem
    dydcsi = dNdcsi * yelem
    dydeta = dNdeta * yelem

    # Definição da matriz Jacobiana
    Jac = np.matrix([[dxdcsi[0, 0], dydcsi[0, 0]],
                    [dxdeta[0, 0], dydeta[0, 0]]])

    # Saída de dados da função "MatDefDesloc"
    return Jac


# Definição da função "derivFuncFormaCart" com dados de entrada das funções "DerivFuncFormaNat", "MatJacobiana2D" e "readmesh"
def derivFuncFormaCart(dNdcsi, dNdeta, Jac, Nnodes):

    # Geração de matriz linha nula que irá receber as derivadas parciais das funções de forma em coordenadas naturais
    dNnat = np.matrix(np.zeros((2, Nnodes)))

    for i in range(Nnodes):

        dNnat[0, i] = dNdcsi[0, i]
        dNnat[1, i] = dNdeta[0, i]

    # Definição da inversa da matriz Jacobiana
    invJac = np.linalg.inv(Jac)

    # Definição da matriz de derivadas parciais das funções de forma em coordenadas cartesianas
    dNdcart = invJac * dNnat

    # Definição dos vetores das derivadas parciais das funções de forma das coordenadas cartesianas em X e Y
    dNdx = dNdcart[0, :]
    dNdy = dNdcart[1, :]

    # Saída de dados da função "derivFuncFormaCart"
    return dNdx, dNdy


# Definição da função "get_globalK" com dados de entrada das funções "readmesh", "nodeindex", "Coord", "dofdrive" e "MEF_ep"
def get_globalK(Nnodes, NGP, t, NELE, dofelem, assmtrx, restrs, NDoF, X, Y, ifPlast,
                csi, eta, w, Cel, sigma_total, connect, f_int, dD):

    K = np.matrix(np.zeros((NDoF, NDoF)))  # Definição da matriz de rigidez global nula
    XYelem = np.matrix(np.zeros((Nnodes, 2)))  # Definição da matriz de coordenadas do elemento
    dDelem = np.matrix(np.zeros((dofelem, 1)))

    # Cálculo do número de pontos de integração
    numPts = NGP * NGP

    for iele in range(NELE):

        Kelem = np.matrix(np.zeros((dofelem, dofelem)))
        fe_int = np.matrix(np.zeros((dofelem, 1)))  # Força interna no elemento

        for idof in range(dofelem):  # Faz o assembly do vetor de deslocamento do elemento

            dDelem[idof, 0] = dD[int(assmtrx[idof, iele]), 0]

        i_sigma = iele * numPts

        for i in range(Nnodes):

            XYelem[i, 0] = X[int(connect[i, iele]), 0]
            XYelem[i, 1] = Y[int(connect[i, iele]), 0]

        for i in range(numPts):

        # Obtém as derivadas das funções de forma em relação às coordenadas naturais
            [dNdcsi, dNdeta] = DerivFuncFormaNat(csi[i], eta[i], Nnodes)

        # Obtém a matriz Jacobiana
            Jac = MatJacobiana2D(dNdcsi, dNdeta, XYelem)

        # Obtém as derivadas das funções de forma em relação às coordenadas cartesianas
            [dNdx, dNdy] = derivFuncFormaCart(dNdcsi, dNdeta, Jac, Nnodes)

        # Obtém a matriz B, que relaciona deslocamentos nodais com deformações em um ponto do elemento
            B = MatDefDesloc(dNdx, dNdy, Nnodes)

        # Calcula o determinante da matriz Jacobiana
            J = np.linalg.det(Jac)

        # Calcula a matriz B transposta
            Bt = B.T

        # Verifica estado e calcula a matriz constitutiva
            if ifPlast[iele, i] == 1:

                sigma_elem = sigma_total[:, i_sigma + i]
                sigma_xx = sigma_elem[0, 0]
                sigma_yy = sigma_elem[1, 0]
                sigma_xy = sigma_elem[2, 0]
                df = np.matrix([[(2 * sigma_xx - sigma_yy) / 3], [(2 * sigma_yy - sigma_xx) / 3], [2 * sigma_xy]])
                Cep = Cel - (Cel * df * df.T * Cel)/(df.T * Cel * df)
                C = Cep

            elif ifPlast[iele, i] == 0:
                C = Cel

            if iele == 0 and i == 0:

                print("TENSÃO ANTERIOR",'-'*57, "\n", sigma_total[:, i_sigma + i], "\n")
                print("DELTA_SIGMA",'-'*61, "\n", C * B * dDelem, "\n")

        # Calcula a matriz de rigidez através de integração numérica
            Kelem += Bt * C * B * float(J * t * w[i, 0])
            fe_int += Bt * sigma_total[:, i_sigma + i] * float(J * t * w[i, 0])

        for ildof in range(dofelem):
            igdof = int(assmtrx[ildof, iele])

            for jldof in range(dofelem):
                jgdof = int(assmtrx[jldof, iele])
                K[igdof, jgdof] += Kelem[ildof, jldof]

        for ildof in range(dofelem):
            igdof = int(assmtrx[ildof, iele])
            f_int[igdof, 0] += fe_int[ildof, 0]

    for i in range(len(restrs)):

        if restrs[i,1] != 0:
            igdof = int(2 * restrs[i,0])
            K[igdof, igdof] = 1
            f_int[igdof] = 0

            if igdof == 0:
                K[igdof, 1: NDoF] = 0
                K[1: NDoF, igdof] = 0

            else:
                K[igdof, 0:igdof] = 0
                K[igdof, (igdof + 1): NDoF + 1] = 0
                K[0: igdof, igdof] = 0
                K[(igdof + 1): NDoF, igdof] = 0

        if restrs[i,2] != 0:
            igdof = int(2 * restrs[i,0] + 1)
            K[igdof, igdof] = 1
            f_int[igdof] = 0

            if igdof == 0:
                K[igdof, 1: NDoF] = 0
                K[1: NDoF, igdof] = 0

            else:
                K[igdof, 0:igdof] = 0
                K[igdof, (igdof + 1): NDoF + 1] = 0
                K[0: igdof, igdof] = 0
                K[(igdof + 1): NDoF, igdof] = 0


    # Saída de dados da função "get_globalK"
    return K, f_int


# Definição da função "MEF_ep" com dados de entrada das funções "readmesh", "nodeindex" e "dofdrive"
def MEF_ep(NDoF, nstep, NELE, connect, Nnodes, NGP, X, Y, dofelem, t, v, E, planestress, assmtrx, fy, forces, restrs):

    F = np.matrix(np.zeros((NDoF, 1)))  # Vetor força total
    f_ext = np.matrix(np.zeros((NDoF, 1)))  # Vetor força externa de cada passo
    sigma_total = np.matrix(np.zeros((3, NGP * NGP * NELE)))  # Tensores das tensões
    D = np.matrix(np.zeros((NDoF, 1)))  # Vetor dos deslocamentos
    dD = np.matrix(np.zeros((NDoF, 1)))  # Vetor de incrementos de deslocamentos
    dDelem = np.matrix(np.zeros((dofelem, 1)))  # Deslocamentos do elemento
    XYelem = np.matrix(np.zeros((Nnodes, 2)))  # Coordenadas do elemento
    ifPlast = np.zeros((NELE, NGP * NGP))  # Matriz booleana que guarda o estado em cada pt de Gauss em cada elemento
    tolD = 0.001  # Tolerância de incremento de deslocamento -----------------------------------------------------------
    j2Elem = np.zeros(NELE)

    # Montagem do vetor Força


    for i in range(len(forces)):

        x1 = X[int(forces[i, 1]), 0]
        x2 = X[int(forces[i, 2]), 0]
        deltaX = abs(x1 - x2)

        y1 = Y[int(forces[i, 1]), 0]
        y2 = Y[int(forces[i, 2]), 0]
        deltaY = abs(y1 - y2)

        r = (((deltaX)**2)+((deltaY)**2))**0.5

        if forces[i,3] != 0: # Teste força na horizontal
            igdof = int(2 * forces[i, 1]) # Nó inicial
            Pnode = (forces[i,3] * r) / 2
            F[igdof, 0] += Pnode

            igdof = int(2 * forces[i, 2]) # Nó final
            F[igdof, 0] += Pnode


        if forces[i,4] != 0: # Teste força na vertical
            igdof = int(2 * forces[i, 1] + 1) # Nó inicial
            Pnode = (forces[i,4] * r) / 2
            F[igdof, 0] += Pnode

            igdof = int(2 * forces[i, 2] + 1) # Nó final
            F[igdof, 0] += Pnode


    # Obtenção dos pontos de Gauss
    [csi1, w1] = PtsGauss1d(NGP)
    k2 = 0
    csi = np.matrix(np.zeros((NGP * NGP, 1)))
    eta = np.matrix(np.zeros((NGP * NGP, 1)))
    w = np.matrix(np.zeros((NGP * NGP, 1)))
    for i in range(NGP):
        for j in range(NGP):
            csi[k2] = csi1[j]
            eta[k2] = csi1[i]
            w[k2] = w1[i] * w1[j]
            k2 += 1

    # Inicializa os parâmetros para o looping
    f_incr = F / nstep
    ky = fy/np.sqrt(3)  # Tensão de escoamento do material em estado de cisalhamento puro
    Cel = constMtrx(planestress, E, v)  # Matriz constituiva elástica


    # Inicia Newthon-Raphson para cálculo dos deslocamentos
    for istep in range(nstep):
        f_ext += f_incr
        count = 0
        maxdD = 11 #--------------------------------------------------------------------------------------------

        while maxdD > tolD:
            f_int = np.matrix(np.zeros((NDoF, 1)))
            [K, f_int] = get_globalK(Nnodes, NGP, t, NELE, dofelem, assmtrx, restrs, NDoF, X, Y, ifPlast,
                csi, eta, w, Cel, sigma_total, connect, f_int, dD)
            b = f_ext - f_int
            dD = np.linalg.inv(K) * b
            #print(np.max(abs(dD)))

            #print('MATRIZ DE RIGIDEZ GLOBAL','-'*48, "\n", K, '\n')
            #print('VETOR DE CARGAS','-'*57, "\n", b, '\n')
            #print('VETOR DE DESLOCAMENTOS NODAIS','-'*43, "\n", dD, '\n')

            for iele in range(NELE):
                i_sigma = iele * NGP * NGP  # ID do elemento analisado na matriz dos tensores
                for idof in range(dofelem):  # Faz o assembly do vetor de deslocamento do elemento
                    dDelem[idof, 0] = dD[int(assmtrx[idof, iele]), 0]

                for i in range(Nnodes):
                    XYelem[i, 0] = X[int(connect[i, iele]), 0]
                    XYelem[i, 1] = Y[int(connect[i, iele]), 0]

                for m in range(NGP * NGP):

                    # Obtém as derivadas das funções de forma em relação às coordenadas naturais
                    [dNdcsi, dNdeta] = DerivFuncFormaNat(csi[m], eta[m], Nnodes)

                    # Obtém a matriz Jacobiana
                    Jac = MatJacobiana2D(dNdcsi, dNdeta, XYelem)

                    # Obtém as derivadas das funções de forma em relação às coordenadas cartesianas
                    [dNdx, dNdy] = derivFuncFormaCart(dNdcsi, dNdeta, Jac, Nnodes)

                    # Obtém a matriz B, que relaciona deslocamentos nodais com deformações em um ponto do elemento
                    B = MatDefDesloc(dNdx, dNdy, Nnodes)

                    # Calcula as deformações e tensões
                    if ifPlast[iele, m] == 1:
                        sigma_elem = sigma_total[:, i_sigma + m]
                        sigma_xx = sigma_elem[0, 0]
                        sigma_yy = sigma_elem[1, 0]
                        sigma_xy = sigma_elem[2, 0]
                        df = np.matrix(
                            [[(2 * sigma_xx - sigma_yy) / 3], [(2 * sigma_yy - sigma_xx) / 3], [2 * sigma_xy]])
                        Cep = Cel - (Cel * df * df.T * Cel) / (df.T * Cel * df)
                        C = Cep

                    elif ifPlast[iele, m] == 0:
                        C = Cel

                    strain = B * dDelem
                    d_sigma = C * strain

                    # Atualiza as tensões
                    sigma_teste = sigma_total[:, i_sigma + m] + d_sigma
                    sigma_xx = sigma_teste[0, 0]
                    sigma_yy = sigma_teste[1, 0]
                    sigma_xy = sigma_teste[2, 0]

                    # Verifica se o mateial escoou
                    J2 = 1/6 * ((sigma_xx - sigma_yy) ** 2 + sigma_yy ** 2 + sigma_xx ** 2) + sigma_xy ** 2 # Calcula o segundo invariante das tensões
                    f_yield = J2 - (ky**2) # ORIGINAL: f_yield = np.sqrt(J2) - ky
                    j2Elem[iele] += np.sqrt(J2)

                    if f_yield > 0:
                        sigma_trial = sigma_total[:, i_sigma + m]
                        k_trial = np.sqrt(J2)  # Calcula o k "teste"
                        sigma_mean = np.matrix([[(sigma_xx + sigma_yy) / 3], [(sigma_xx + sigma_yy) / 3], [0]])
                        dev = sigma_trial - sigma_mean # ORIGINAL: dev = sigma_trial - np.matrix([[(2 * sigma_xx + sigma_yy) / 3], [(sigma_xx + 2 * sigma_yy) / 3], [0]])
                        aux = ky/k_trial * dev
                        sigma_total[:, i_sigma + m] = aux + sigma_mean # ORIGINAL: sigma_total[:, i_sigma + m] = aux + np.matrix([[(2 * sigma_xx + sigma_yy) / 3], [(sigma_xx + 2 * sigma_yy) / 3], [0]])
                        ifPlast[iele, m] = 1

                    elif f_yield <= 0:
                        sigma_total[:, i_sigma + m] = sigma_teste
                        ifPlast[iele, m] = 0

            D += dD

            if np.any(ifPlast) == False:
                print("STATUS: Passo", istep, "sem plastificação em algum ponto. Looping", count, "encerrado", '\n')
                break

            if np.any(ifPlast) == True:
                print("STATUS: Passo", istep, "com plastificação em algum ponto. Looping", count, "continuado", '\n')
                continue

            maxdD = abs(np.max(dD))

            if count > 20:
                print("STATUS: Looping encerrado por limite de interações atingido", '\n')
                break

            count += 1
            print(count)

    data = open("j2.txt", "w")
    for j2 in j2Elem:
        j2 /= 4
        data.write(f"{str(j2)}\n")
    data.close()

    # Saída de dados da função "MEF_ep"
    return D, sigma_total, j2Elem

