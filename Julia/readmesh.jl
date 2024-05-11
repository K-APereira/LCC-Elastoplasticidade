module meshReader
    import JSON3

    # function that read data from given Json
    function readMesh(inputJson)
        JsonData = JSON3.read(inputJson)

        ###### MATERIAL
        # Read the materials properties (Young's module, poisson's ratio, thickness and Yield stress)
        Props = JsonData["Props"]

        ###### MODEL
        # Read the nodes coordinates (X,Y)
        NodesCoord = copy(JsonData["NodesCoord"])
        N_Nodes = length(NodesCoord)

        # Read the elements connection (left_down, right_down, right_top, left_top)
        Connect = copy(JsonData["Elements"])
        N_Elems = length(Connect)
        N_NodesInElem = length(Connect[1])

        # Read the nodes restrictions (NodeId, bool X lock, bool Y lock)
        Restrs = copy(JsonData["Restrs"])
        N_Restrs = length(Restrs)

        # Read the loads (ElemId, NodeId1, NodeId2, Force in X, Force in Y)
        Forces = copy(JsonData["Forces"])

        # Read if the model is using Plane stress or Plane strain
        PlaneStressOrStrain = JsonData["PlaneStressOrStrain"]

        # Read the degrees of freedom (DoF) for node, the DoF for element and the model's total DoF
        DoFNode = JsonData["DoFNode"]
        DoFElem = DoFNode*N_NodesInElem
        N_DoF = DoFNode*N_Nodes

        ###### METHOD
        N_Steps = JsonData["N_Steps"]

        # Read the number of Gauss points to be used
        NGP = JsonData["NGP"]


        return (N_Nodes,NodesCoord,N_Restrs,Restrs,Props,N_NodesInElem,N_Elems,Connect,DoFNode,
        DoFElem,N_DoF,Forces,N_Steps,PlaneStressOrStrain,NGP)
    end


    function write_results(filepath, results...;pretty = false)
        data = Dict()
        for i in results
            data[i[1]] = i[2]
        end
        file = open(filepath,"w")
        if pretty
            pretty_json(data,file)
        else
            JSON3.write(file,data)
        end
        close(file)
    end

    function pretty_json(data,f,indent=0)
        if data isa Dict
            write(f,"\n","    "^indent,"{")
            for (key,value) in data
                write(f,"\n","    "^(indent+1),"\"",key,"\"",": ")
                pretty_json(value,f,indent+1)
                if last(collect(keys(data))) != key
                    write(f,",")
                end
                write(f,"\n")
            end
            write(f,"    "^indent,"}")
        elseif data isa Array
            write(f,"\n","    "^indent,"[")
            for i in data[1: length(data)-1]
                pretty_json(i,f,indent+1)
                write(f,", ")
            end
            if length(data)>0
                pretty_json(data[length(data)],f,indent+1)
            end
            write(f,"]")
        else
            if data isa String
                write(f,"\"",data,"\"")
            else
                write(f,string(data))
            end
        end
    end
end