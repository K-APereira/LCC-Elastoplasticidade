module MEF
    function Mef_ep(N_DoF, N_Steps, N_Elems, Connect, N_NodesInElem, NGP, NodesCoord, DoFNode, Props, PlaneStressOrStrain, assmtrx, Forces, Restrs)
        
        # Setting matrix for function's use

        F = zeros((N_DoF,1)) # total force for node
        f_ext = zeros((N_DoF,1)) # external force for node in each step
        total_sigma = zeros((3, NGP * NGP * N_Elems)) # stress tensor
        D = zeros((N_DoF,1)) # total displacement
        dD = zeros((N_DoF,1)) # displacement for step
        dD_Elem = zeros((DoFNode,1)) # elements' displacement for step
        XY_Elem = zeros((N_NodesInElem, 2)) # elemets' coords
        if_Plast = zeros(Bool,(N_Elems,NGP*NGP)) # state of each gauss point for elem
        tolD = 0.001 # displacement increment tolerance
        J2Elem = zeros(N_Elems) # von Mises yield criterion


        # filling the total force vector

        x1 = 0
        for i in eachindex(Forces)

            # calculates the area where the force is being applied
            x1 = NodesCoord[Int(Forces[i][2])][1]
            x2 = NodesCoord[Int(Forces[i][3])][1]
            deltaX = abs(x1-x2)

            y1 = NodesCoord[Int(Forces[i][2])][2]
            y2 = NodesCoord[Int(Forces[i][3])][2]
            deltaY = abs(y1-y2)

            r = ((deltaX^2)+(deltaY^2))^0.5

            # if there is force in X, add to the vector F
            if Forces[i][4] != 0
                F_Node = (Forces[i][4] * r)/2 # half force for each node

                Node = 2 * Forces[i][2] - 1
                F[Node,1] += F_Node
                Node = 2 * Forces[i][3] - 1
                F[Node,1] += F_Node
                
            end

            # if there is force in Y, add to the vector F
            if Forces[i][5] != 0
                F_Node = (Forces[i][5] * r)/2 # half force for each node

                Node = 2 * Forces[i][2]
                F[Node,1] += F_Node
                Node = 2 * Forces[i][3]
                F[Node,1] += F_Node
                
            end

        end


        return D, total_sigma,F
    end
end