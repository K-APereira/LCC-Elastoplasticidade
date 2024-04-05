# Importação da biblioteca numpy
import numpy as np


# Definição da função "readmesh" com arquivo de entrada "filename"
def readmesh(FEMEP, NEUTRO):

    #------DADOS FEMEP-------

    # Atribuição de variável ao arquivo de entrada aberto
    inpfile = open(FEMEP)

    # Início de processo iterativo para leitura de cada linha do arquivo de entrada aberto
    for line in inpfile:

        # Definição dos nós do modelo
        if line.find("%NODE.COORD") > -1:
            mesh_nodes = int(inpfile.readline())  # Quantidade de nós do modelo
            coord_nodes = np.zeros((mesh_nodes, 3))  # [Número do nó, Coordenada X, Coordenada Y]
            for k in range(mesh_nodes):
                [nodeID, coord_x, coord_y, coord_z] = inpfile.readline().split()
                coord_nodes[k, 0] = int(nodeID) - 1 # ID do nó já com "Zero index" corrigido
                coord_nodes[k, 1] = float(coord_x) # Coordenada X
                coord_nodes[k, 2] = float(coord_y) #Coordenada Y

        # Definição das condições de contorno do modelo
        elif line.find("%NODE.SUPPORT") > -1:
            num_restrs = int(inpfile.readline())  # Quantidade de nós do modelo com restrições
            restrs = np.zeros((num_restrs, 3))  # Matriz que armazena as restrições dos nós do modelo
            for k in range(num_restrs):
                [nodeID, dx, dy, unknown1, unknown2, unknown3, rz] = inpfile.readline().split()
                restrs[k, 0] = int(nodeID) - 1  # Armazena a ID do nó já com "Zero index" corrigido
                restrs[k, 1] = int(dx)  # Armazena se há restrição do nó na direção X
                restrs[k, 2] = int(dy)  # Armazena se há restrição do nó na direção Y

        # Definição das propriedades do material (Módulo de elast., Coef. de Poisson e Espess.)
        elif line.find("%MATERIAL.ISOTROPIC") > -1:
            num_isomat = int(inpfile.readline())  # Quantidade de materiais isotrópicos distintos
            props = np.zeros((num_isomat, 4))  # 4 propriedades por elemento (E, v, t, fy)
            for k in range(num_isomat):
                [isomatID, E, v] = inpfile.readline().split()  # [ID do material isotr., E, v]
                props[k, 0] = float(E)
                props[k, 1] = float(v)
        elif line.find("%THICKNESS") > -1:
            num_thick = int(inpfile.readline())  # Quantidade de espessuras distintas para os elementos
            for k in range(num_thick):
                [thicknessID, t] = inpfile.readline().split()  # [ID da espessura, Espessura]
                props[k, 2] = float(t)

        elif line.find("%ELEMENT.Q4") > -1:
            elem_nodes = 4
            num_elem = int(inpfile.readline())  # Número de elementos da malha de EF
            connect = np.zeros((num_elem, elem_nodes))  # Definição da matriz conectividade
            for k in range(num_elem):
                [nodeID, unk1, unk2, unk3, left_down, right_down, right_top, left_top] = inpfile.readline().split()
                connect[k, 0] = int(left_down) - 1
                connect[k, 1] = int(right_down) - 1
                connect[k, 2] = int(right_top) - 1
                connect[k, 3] = int(left_top) - 1
            connect = connect.T

        elif line.find("%ELEMENT.Q8") > -1:
            elem_nodes = 8
            num_elem = int(inpfile.readline())  # Número de elementos da malha de EF
            connect = np.zeros((num_elem, elem_nodes))  # Definição da matriz conectividade
            for k in range(num_elem):
                [nodeID, unk1, unk2, unk3, left_top, left_center, left_bottom, center_bottom, right_bottom,
                 right_center, right_top, center_top ] = inpfile.readline().split()
                connect[k, 0] = int(left_top) - 1
                connect[k, 1] = int(left_center) - 1
                connect[k, 2] = int(left_bottom) - 1
                connect[k, 3] = int(center_bottom) - 1
                connect[k, 4] = int(right_bottom) - 1
                connect[k, 5] = int(right_center) - 1
                connect[k, 6] = int(right_top) - 1
                connect[k, 7] = int(center_top) - 1
            connect = connect.T



        elif line.find("%LOAD.CASE.LINE.FORCE.UNIFORM") > -1:
            force_elem = int(inpfile.readline())  # Número de elementos com carga uniformemente distribuída
            forces = np.zeros((force_elem, 5))  # Matriz que armazena as direções de aplicação das forças nos nós
            for k in range(force_elem):

                """1ª COLUNA: ID do elemento carregado;
		        2ª COLUNA: ID global do primeiro nó na aresta carregada;
		        3ª COLUNA: ID global do segundo nó na aresta carregada;
		        5ª COLUNA: direção X (magnitude associada a cada aresta [N/mm])
		        6ª COLUNA: direção Y (magnitude associada a cada aresta [N/mm])"""

                [elemID, node1ID, node2ID, unk4, Fx, Fy, unk5] = inpfile.readline().split()
                forces[k, 0] = int(elemID) - 1  # Armazena a ID do elemento já com "Zero index" corrigido
                forces[k, 1] = int(node1ID) - 1  # Armazena a ID do 1º nó já com "Zero index" corrigido
                forces[k, 2] = int(node2ID) - 1  # Armazena a ID do 2º nó já com "Zero index" corrigido
                forces[k, 3] = float(Fx) # Magnitude associada a força em x de cada aresta
                forces[k, 4] = float(Fy) # Magnitude associada a força em y de cada aresta

    inpfile.close()


    #------DADOS AUSENTES NO FEMEP-------

    # Atribuição de variável ao arquivo de entrada aberto
    inpfile = open(NEUTRO)

    # Início de processo iterativo para leitura de cada linha do arquivo de entrada aberto
    for line in inpfile:

        # Definição das propriedades do material (Tensão de escoamento)
        if line.find("%YIELD.STRESS") > -1:
            num_escoam = int(inpfile.readline())  # Quantidade de tensões de escoamento do modelo
            for k in range(num_escoam):
                fy = inpfile.readline() # Definição da tensão de escoamento dos materiais
                props[k, 3] = float(fy)

        # Definição das parcelas da força externa
        elif line.find("%NSTEP") > -1:
            nstep = int(inpfile.readline())

        # Definição do tipo de estado plano em análise
        elif line.find("%PLANESTRESS") > -1:
            planestress = int(inpfile.readline())

        # Definição do número de pontos de integração de Gauss
        elif line.find("%NGP") > -1:
            NGP = int(inpfile.readline())

    inpfile.close()


    #------DEFINIÇÃO DE VARIÁVEIS INDEPENDENTES-------

    # Definição do número de graus de liberdade (g.d.l.) por nó do elemento
    Dofnode = 2

    # Definição do número de g.d.l. por elemento local
    dofelem = int(Dofnode * elem_nodes)

    # Definição do número de g.d.l. globais
    NDoF = int(Dofnode * mesh_nodes)

    # Definição dos comprimentos das arestas com carga uniformemente distribuída
    """length_edges_forces = []
    for k in range(forces.shape[0]):
        deltaX = abs(coord_nodes[int(forces[k,1]),1] - coord_nodes[int(forces[k,2]),1])
        deltaY = abs(coord_nodes[int(forces[k,1]),2] - coord_nodes[int(forces[k,2]),2])
        length_edges_forces.append([deltaX, deltaY])
    length_edges_forces = np.array(length_edges_forces)"""


    #print(mesh_nodes,'\n', coord_nodes,'\n', num_restrs,'\n', restrs,'\n', num_isomat,'\n', props,'\n', num_thick,'\n',
          #elem_nodes,'\n', num_elem,'\n', connect,'\n', dofelem,'\n', NDoF,'\n', forces,'\n', nstep,'\n', planestress,'\n', NGP)

    return mesh_nodes,coord_nodes,num_restrs,restrs,num_isomat,props,num_thick,elem_nodes,num_elem,connect,Dofnode\
        ,dofelem,NDoF,forces,nstep,planestress,NGP

#-----------------------------------x=readmesh('exemplocook_planeStress.txt', "Model_Tk.nf")

# Definição da função "dofdrive" com dados de entrada das funções "readmesh" e "nodeindex"
def dofdrive(Nnodes, NELE, DomNodeID, dofelem):

    # Definição da Matriz de Montagem (Assembly) nula
    assmtrx = np.zeros((dofelem, NELE))   #Gera matriz ASSMtrx

    # Preenchimento da Matriz de Montagem
    for iele in range(NELE):
        ildof = 0
        for ilnode in range(Nnodes):
            ignode = DomNodeID[ilnode, iele]
            ildof = ildof + 2
            assmtrx[ildof - 2, iele] = 2 * ignode
            assmtrx[ildof - 1, iele] = 2 * ignode + 1

    # Saída de dados da função "dofdrive"
    return assmtrx