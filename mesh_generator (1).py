import numpy as np
import cv2
from viewv1 import *

FILENAME = "imagem3x3"

#Lendo imagem
img = cv2.imread(FILENAME + '.tiff',0)
img=np.array(img,dtype=np.uint8)

print(img.shape)
reduction_factor = 1

img_reduzida=np.zeros(shape=(img.shape[0]//reduction_factor,img.shape[1]//reduction_factor))

reduction=np.zeros(shape=(reduction_factor,reduction_factor))
for i in range(img_reduzida.shape[0]):
    for j in range(img_reduzida.shape[1]):
        reduction[:,:] = img[i*reduction_factor:(i+1)*reduction_factor,j*reduction_factor:(j+1)*reduction_factor]
        vals, counts = np.unique(reduction,return_counts=True)
        mode = np.argmax(counts)
        img_reduzida[i,j] = vals[mode]

img = img_reduzida
print(img.shape)
void = 255

#Tamanho do pixel em dimensoes reais
voxel_size=1

num_elems_pre=img.shape #num de elems=dim das img
num_nodes_pre=[num_elems_pre[0]+1,num_elems_pre[1]+1] #num de nodes

#Verifica posicao dos nodes
nodes=np.zeros(shape=(num_nodes_pre),dtype=np.uint16)
for i in range(num_elems_pre[0]):
    for j in range(num_elems_pre[1]):
        if(img[j,i]!=void):
            nodes[j,i]=True
            nodes[j,i+1]=True
            nodes[j+1,i]=True
            nodes[j+1,i+1]=True



#Gera coord_nodes
nNodes=0
for i in range(num_nodes_pre[0]):
    for j in range(num_nodes_pre[1]):
        if(nodes[j,i]):
            nNodes+=1
coord_nodes = np.zeros((nNodes, 3))

num_nodes=0
for i in range(num_nodes_pre[0]):
    for j in range(num_nodes_pre[1]):
        if(nodes[j,i]):
            nodes[j,i]=num_nodes
            coord_nodes[num_nodes, 0] = num_nodes
            coord_nodes[num_nodes, 1] = i * voxel_size
            coord_nodes[num_nodes, 2] = (num_nodes_pre[1] - 1 - j) * voxel_size
            num_nodes+=1

# Gera connect
nElems=0
for i in range(num_elems_pre[0]):
    for j in range(num_elems_pre[1]):
        if(img[j,i]!=void):
            nElems+=1
connect = np.zeros((nElems, 4))

print(nElems)

#Gera elems
num_elems=0
for i in range(num_elems_pre[0]):
    for j in range(num_elems_pre[1]):
        if(img[j,i]!=void):
            connect[num_elems, 0] = nodes[j,i]
            connect[num_elems, 1] = nodes[j+1,i]
            connect[num_elems, 2] = nodes[j+1,i+1]
            connect[num_elems, 3] = nodes[j,i+1]
            num_elems+=1
            #print(num_elems,nodes[j,i],nodes[j+1,i],nodes[j+1,i+1],nodes[j,i+1])
connect = connect.T

file = open(FILENAME + '.txt', 'w')
file.write("%HEADER.ANALYSIS\n"
            "'plane_stress'\n\n"
            "%NODE\n"
            f"{nNodes}\n\n"
            "%NODE.COORD\n"
            f"{nNodes}\n"
            )
for node in coord_nodes:
    file.write(f"{int(node[0] + 1)} {node[1]} {node[2]} 0.0\n")
viewmesh(coord_nodes[:,[1,2]], connect, 4)
print("Quantidade de nós com restrição em X: ")
nRestrsX = int(input())
restrs = np.zeros((nRestrsX, 3))
print("Digite o Id do primeiro nó restrito: ")
firstRestId = int(input())
for i in range(nRestrsX):
    restrs[i, 0] = firstRestId + i
    restrs[i, 1] = 1
    restrs[i, 2] = 1
file.write("\n%NODE.SUPPORT\n"
           f"{nRestrsX}\n"
           )
for restr in restrs:
    file.write(f"{int(restr[0])} {int(restr[1])} {int(restr[2])} 0 0 0 0\n")


file.write("\n%MATERIAL\n"
           f"1\n\n"
           "%MATERIAL.LABEL\n"
           f"1\n"
           f"1 'M1'\n"
           "\n%MATERIAL.ISOTROPIC\n"
           "1\n"
           )

print(f"Digite o módulo de elasticidade do material: ")
E = float(input())
print(f"Digite o coeficiente de Poisson do material: ")
v = float(input())
print(f"Digite a espessura do material: ")
thickness = float(input())
file.write(f"1 {E} {v}\n"
           "\n%THICKNESS\n1\n"
           )
file.write(f"1 {thickness}\n")
file.write("\n%ELEMENT\n"
            f"{nElems}\n\n"
            "%ELEMENT.Q4\n"
            f"{nElems}\n")
for k in range(nElems):
    file.write(f"{k+1} 1 1 3 {int(connect[0, k] + 1)} {int(connect[1, k] + 1)} {int(connect[2, k]+1)} {int(connect[3, k]+1)}\n")  
    
file.write("\n%LOAD.CASE.LINE.FORCE.UNIFORM")
print("Digite a quantidade de elementos com cargas distribuídas: ")
nLoads = int(input())
print("Digite o eixo da carga distribuida (X ou Y): ")
axis = input()
file.write(f"\n{nLoads}\n")

if axis == 'X':
    print("Digite o ID do menor elemento com carga distribuída em X: ")
    elemId = int(input())
    print("Digite o valor da carga distribuída em X: ")
    load = float(input())
    for i in range(nLoads):
        file.write(f"{elemId + i} {int(connect[3,(elemId + i - 1)] + 1)} {int(connect[2, (elemId + i - 1)] + 1)} 0 {load} 0 0\n")

if axis == 'Y':
    print("Digite o ID do menor elemento com carga distribuída em Y: ")
    elemId = int(input())
    print("Digite o valor da carga distribuída em Y: ")
    load = float(input())
    for i in range(nLoads):
        file.write(f"{elemId + i} {int(connect[2,(elemId + i - 1)] + 1)} {int(connect[1, (elemId + i - 1)] + 1)} 0 0 {load} 0\n")

file.write("\n%END")
file.close()