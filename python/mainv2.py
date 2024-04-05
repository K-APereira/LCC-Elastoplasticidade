# Importação das biblioteca e funções
import timeit
from readmeshv2 import *
from viewv1 import *
from MEFv2 import *
#from windowv1 import *

# Início da contagem do tempo de execução do código
inicio = timeit.default_timer()

# Janela para entrada de dados
#janela_v1()

# Atribuição de variáveis de acordo com a saída de dados da função "readmesh"

[NODES,coord_nodes,num_restrs,restrs,num_isomat,props,num_thick,elem_nodes,NELE,connect,Dofnode\
        ,dofelem,NDoF,forces,nstep,planestress,NGP] = readmesh('Exemplo 3.3.2 - 256 elementos.txt', "Model_Tk.nf")


# Atribuição de variáveis de acordo com a saída de dados da função "dofdrive"
assmtrx = dofdrive(elem_nodes, NELE, connect, dofelem)

# Visualização da saída de dados da função "dofdrive" para usuário
#print('-'*26, 'MATRIZ DE MONTAGEM','-'*26)
#print(assmtrx,'\n')

XY = coord_nodes[:,[1,2]]
X = coord_nodes[:,[1]]
Y = coord_nodes[:,[2]]

#print(connect[1, 0])
#print(type(connect[0, 0]))

# Atribuição de variáveis de acordo com a saída de dados da função "viewmesh"
#viewmesh(XY, connect, elem_nodes)

# Atribuição de variáveis auxiliares de acordo com a saída de dados da função "readmesh"
E = props[0, 0]
v = props[0, 1]
t = props[0, 2]
fy = props[0, 3]

# Atribuição de variáveis de acordo com a saída de dados da função "PtsGauss1d"
[csi, w] = PtsGauss1d(NGP)

# Atribuição de variáveis de acordo com a saída de dados da função "MEF_ep"
[D, sigma_total, j2Elem] = MEF_ep(NDoF, nstep, NELE, connect, elem_nodes, NGP, X, Y, dofelem, t, v, E, planestress, assmtrx, fy, forces, restrs)
print(D, sigma_total)

# Atribuição de variáveis de acordo com a saída de dados da função "viewmesh" e "MEF_ep"
#viewdeformedmesh(XY, connect, D, 10,elem_nodes)

viewcolormesh(XY, connect, elem_nodes, j2Elem)

# Fim da contagem do tempo de execução do código
fim = timeit.default_timer()

print("Tempo de execução: {} s".format(fim-inicio))





