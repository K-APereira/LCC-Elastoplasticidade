#importing modules
include("MEF.jl")
using .MEF
include("readmesh.jl")
using .meshReader
include("view.jl")
using .view

function main(meshInput,savepath = nothing; pretty_save = false)
 
    # Reading data from presented Json
    (N_Nodes,NodesCoord,N_Restrs,Restrs,Props,N_NodesInElem,N_Elems,Connect,DoFNode,
    DoFElem,N_DoF,Forces,N_Steps,PlaneStressOrStrain,NGP) = meshReader.readMesh(meshInput)
    

    @time (D, sigma_total, j2Elem) = MEF.FEM_Ep(N_DoF, N_Steps, N_Elems, Connect, N_NodesInElem, NGP, NodesCoord, DoFElem, Props, PlaneStressOrStrain, DoFNode, Forces, Restrs)

    
    
    # if isnothing(savepath)
    #     println("\nD:")
    #     println(D)
    #     println("\nSigma_total:")
    #     println(sigma_total)
    #     println("\nJ2Elem:")
    #     println(j2Elem)
    # else
    #     meshReader.write_results(savepath,("D",D),("sigma_total",sigma_total),("j2Elem",j2Elem);pretty = pretty_save)
    #     println("\nResults saved in: ", savepath)
    # end

    view.viewMesh(NodesCoord, Connect, D, 10, j2Elem)
end


########################################
##### MAIN #####
######################################## 
#run with timer
# @time main(raw"Julia\Exemplos\Exemplo 1.1 - 100 elementos.json",raw"Julia\Exemplos\Resultado2 Exemplo 1.1 - 100 elementos.json",pretty_save = true)
@time main(raw"Julia\imagem3x3.json")