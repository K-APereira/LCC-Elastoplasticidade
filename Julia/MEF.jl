module MEF

    # function that calculates the gauss points to integrate in FEM_Ep
    function Gauss_Pts(NGP)

        # setting the vector of gauss points' locations and weights
        location = zeros(NGP)
        weight = zeros(NGP)

        # filling the vectors for each number of gauss points
        if NGP == 1
            location[1] = 0
            weight[1] = 2

        elseif NGP == 2
            location[1] = -1/sqrt(3)
            location[2] = -location[1]

            weight[1] = 1
            weight[2] = 1

        elseif NGP == 3
            location[1] = -0.774596669241483
            location[2] = 0
            location[3] = -location[1]

            weight[1] = 5/9
            weight[2] = 8/9
            weight[3] = 8/9

        elseif NGP == 4
            location[1] = -0.861136311594053
            location[2] = -0.339981043584856
            location[3] = 0.339981043584856
            location[3] = 0.861136311594053

            weight[1] = 0.347854845137454
            weight[2] = 0.652145154862546
            weight[3] = 0.652145154862546
            weight[4] = 0.347854845137454

        end

        return location, weight
    end
    
    # function that generates the constitutive matrix for the material and state
    function ConstMtrx(PlaneStressOrStrain, Props)

        E = Props["E"]
        v = Props["v"]

        C = zeros((3,3))

        # select state between Plane Stress or Plane Strain

        if PlaneStressOrStrain == "Plane Stress"
            aux = (E/(1-(v^2)))
            C[1,1] = aux
            C[1,2] = v*aux
            C[1,3] = 0
            C[2,1] = v*aux
            C[2,2] = aux
            C[2,3] = 0
            C[3,1] = 0
            C[3,2] = aux
            C[3,3] = ((1-v)/2)*aux

        elseif PlaneStressOrStrain == "Plane Strain"
            aux = (E * (1 - v) / ((1 + v) * (1 - 2 * v)))
            C[1,1] = aux
            C[1,2] = (v/(1-v))*aux
            C[1,3] = 0
            C[2,1] = (v/(1-v))*aux
            C[2,2] = aux
            C[2,3] = 0
            C[3,1] = 0
            C[3,2] = 0
            C[3,3] = ((1 - 2 * v)/(2 * (1 - v)))*aux
            
        end

        return C
    end

    function DerivNatShapeFunc(csi, eta,N_NodesInElem)
        dNdcsi = zeros(N_NodesInElem)
        dNdeta = zeros(N_NodesInElem)

        if N_NodesInElem==8
            dPet = 1 + eta
            dMet = 1 - eta
            dPxs = 1 + csi
            dMxs = 1 - csi

            dNdcsi[5] = -csi * dMet
            dNdcsi[6] = 0.50 * dMet * dPet
            dNdcsi[7] = -csi * dPet
            dNdcsi[8] = -0.50 * dMet * dPet

            dNdcsi[1] = -0.25 * dMet - 0.5 * (dNdcsi[5] + dNdcsi[8])
            dNdcsi[2] = 0.25 * dMet - 0.5 * (dNdcsi[5] + dNdcsi[6])
            dNdcsi[3] = 0.25 * dPet - 0.5 * (dNdcsi[6] + dNdcsi[7])
            dNdcsi[4] = -0.25 * dPet - 0.5 * (dNdcsi[7] + dNdcsi[8])

            dNdeta[5] = -0.50 * dMxs * dPxs
            dNdeta[6] = -eta * dPxs
            dNdeta[7] = 0.50 * dMxs * dPxs
            dNdeta[8] = -eta * dMxs

            dNdeta[1] = -0.25 * dMxs - 0.5 * (dNdeta[5] + dNdeta[8])
            dNdeta[2] = -0.25 * dPxs - 0.5 * (dNdeta[5] + dNdeta[6])
            dNdeta[3] = 0.25 * dPxs - 0.5 * (dNdeta[6] + dNdeta[7])
            dNdeta[4] = 0.25 * dMxs - 0.5 * (dNdeta[7] + dNdeta[8])

        elseif Nnodes == 4
            dPet = 1 + eta
            dMet = 1 - eta
            dPxs = 1 + csi
            dMxs = 1 - csi
    
            dNdcsi[1] = -0.25 * dMet
            dNdcsi[2] = 0.25 * dMet
            dNdcsi[3] = 0.25 * dPet
            dNdcsi[4] = -0.25 * dPet
    
            dNdeta[1] = -0.25 * dMxs
            dNdeta[2] = -0.25 * dPxs
            dNdeta[3] = 0.25 * dPxs
            dNdeta[4] = 0.25 * dMxs
        end

        return dNdcsi,dNdeta
    
    end

    # function that returns the global stiffness matrix and the internal forces
    function Get_GlobalK(N_NodesInElem, NGP, Props, N_Elems, DoFElem, assmtrx, Restrs, N_DoF, NodesCoord, if_Plast,
        csi, eta, w, Cel, total_sigma, Connect, f_int, dD)
        
        K = zeros((N_DoF,N_DoF)) #stiffness matrix
        XY_Elem = zeros(N_NodesInElem, 2) # elemets' coords
        dD_Elem = zeros(DoFElem) # displacement of the nodes in one element
        N_points = NGP*NGP # number of integration points

        for i_elem in 1:N_Elems

            K_elem = zeros((DoFElem,DoFElem)) # element's K
            f_int_elem = zeros(DoFElem) # element's internal force
            i_sigma = (i_elem-1)*N_points

            # get the element i_elem nodes coords
            for i_node in 1:N_NodesInElem
                XY_Elem[i,:] = NodesCoord[int(Connect[i_elem,i_node]),:]
            end

            for i in 1:N_points
                # get the shape form functions derivates in natural coords
                dNdcsi, dNdeta = DerivNatShapeFunc(csi, eta,N_NodesInElem)
            end
        end

        return K, f_int
    end

    # function that calculates the Finite Elements method in elastoplastics conditions
    function FEM_Ep(N_DoF, N_Steps, N_Elems, Connect, N_NodesInElem, NGP, NodesCoord, DoFElem, Props, PlaneStressOrStrain, assmtrx, Forces, Restrs)
        
        # Setting matrix for function's use

        F = zeros((N_DoF)) # total force for node
        f_ext = zeros((N_DoF)) # external force for node in each step
        total_sigma = zeros((3, NGP * NGP * N_Elems)) # stress tensor
        D = zeros((N_DoF)) # total displacement
        dD = zeros((N_DoF)) # displacement for step
        dD_Elem = zeros((DoFElem)) # elements' displacement for step
        XY_Elem = zeros((N_NodesInElem, 2)) # elemets' coords
        if_Plast = zeros(Bool,(N_Elems,NGP*NGP)) # state of each gauss point for elem
        tolD = 0.001 # displacement increment tolerance
        J2Elem = zeros(N_Elems) # von Mises yield criterion
        YieldStress = Props["YieldStress"]


        # filling the total force vector
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
                F[Node] += F_Node
                Node = 2 * Forces[i][3] - 1
                F[Node] += F_Node
                
            end

            # if there is force in Y, add to the vector F
            if Forces[i][5] != 0
                F_Node = (Forces[i][5] * r)/2 # half force for each node

                Node = 2 * Forces[i][2]
                F[Node] += F_Node
                Node = 2 * Forces[i][3]
                F[Node] += F_Node
                
            end

        end

        # Getting Gauss points locations and weights
        (GP_locations, GP_weights) = Gauss_Pts(NGP)
        k2 = 0
        csi = zeros(NGP*NGP)
        eta = zeros(NGP*NGP)
        w = zeros(NGP*NGP)
        for i in 1:NGP
            for j in 1:NGP
                csi[k2] = GP_locations[j]
                eta[k2] = GP_locations[i]
                w[k2] = GP_weights[i]*GP_weights[j]
                k2 += 1
            end
        end

        # loop parameters start
        f_incr = F / N_Steps # incremental force for each step
        ky = YieldStress / sqrt(3) # material's yield stress in pure shear (used in von Mises yield criterion)
        Cel = ConstMtrx(PlaneStressOrStrain, Props)


        # start Newthon-Raphson method to calculate the displacement
        for i_step in 1:N_Steps
            f_ext += f_incr
            count = 0
            maxdD = 11

            while maxdD > tolD
                # internal forces f_int
                f_int = zeros(N_DoF)
                # calculate the stiffness matrix K and the internal forces f_int
                (K, f_int) = Get_GlobalK(N_NodesInElem, NGP, Props, N_Elems, DoFElem, assmtrx, Restrs, N_DoF, NodesCoord, if_Plast,
                csi, eta, w, Cel, total_sigma, Connect, f_int, dD)
                b = f_ext - f_int
                dD = inv(K)*b

                maxdD = 0
            
            end
        end


        return D, total_sigma, dD
        end
    end