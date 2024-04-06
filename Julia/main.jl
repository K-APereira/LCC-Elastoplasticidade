#importing modules
include("MEF.jl")
using .MEF
include("readmesh.jl")
using .meshReader

function main(meshInput)

    # Reading data from presented Json
    (N_Nodes,NodesCoord,N_Restrs,Restrs,Props,N_NodesInElem,N_Elems,Connect,DoFNode,
    DoFElem,N_DoF,Forces,N_Steps,PlaneStressOrStrain,NGP) = meshReader.readMesh(meshInput)
    
    assmtrx = meshReader.Create_AssembMtrx(N_NodesInElem, N_Elems, Connect, DoFNode)

    print(Forces)

    (D, sigma_total, j2Elem) = MEF.Mef_ep(N_DoF, N_Steps, N_Elems, Connect, N_NodesInElem, NGP, NodesCoord, DoFNode, Props, PlaneStressOrStrain, assmtrx, Forces, Restrs)
end


########################################
##### MAIN #####
########################################
# main("Exemplo1.json")
#run with timer
@time main("Julia/Exemplo1.json")