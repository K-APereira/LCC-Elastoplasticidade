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
        NodesCoord = JsonData["NodesCoord"]
        N_Nodes = length(NodesCoord)

        # Read the elements connection (left_down, right_down, right_top, left_top)
        Connect = JsonData["Elements"]
        N_Elems = length(Connect)
        N_NodesInElem = length(Connect[1])

        # Read the nodes restrictions (NodeId, bool X lock, bool Y lock)
        Restrs = JsonData["Restrs"]
        N_Restrs = length(Restrs)

        # Read the loads (ElemId, NodeId1, NodeId2, Force in X, Force in Y)
        Forces = JsonData["Forces"]

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

    # function that creates the Assembly matrix
    # Assembly matrix generates global Id for DoF for each node
    function Create_AssembMtrx(N_NodesInElem, N_Elems, Connect, DoFNode)
        
        # Create the Assembly matrix with zeros
        AssembMtrx = zeros(Int,(N_Elems,N_NodesInElem*DoFNode))


        # Fills the matrix with global Id for each DoF for element
        for i_elem in 1:N_Elems
            ildof = 1
            for i_node in 1:N_NodesInElem
                # For element, for node in this element, get the node Id
                IdNode = Connect[i_elem][i_node]
                
                # For each node, alocate a global Id for its degrees of freedom
                for i_dof in 1:DoFNode
                    AssembMtrx[i_elem,ildof] = DoFNode*(IdNode-1) + i_dof
                    ildof = ildof+1
                end
            end
        end

        return AssembMtrx
    end
end
