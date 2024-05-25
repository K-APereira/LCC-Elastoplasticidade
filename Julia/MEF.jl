module MEF
using LinearAlgebra, FastGaussQuadrature

    # function that calculates the gauss points to integrate in FEM_Ep
    # function Gauss_Pts(NGP)
    #     # setting the vector of gauss points' locations and weights
    #     location = zeros(NGP)
    #     weight = zeros(NGP)

    #     # filling the vectors for each number of gauss points
    #     if NGP == 1
    #         location[1] = 0
    #         weight[1] = 2

    #     elseif NGP == 2
    #         location[1] = -1/sqrt(3)
    #         location[2] = -location[1]

    #         weight[1] = 1
    #         weight[2] = 1

    #     elseif NGP == 3
    #         location[1] = -0.774596669241483
    #         location[2] = 0
    #         location[3] = -location[1]

    #         weight[1] = 5/9
    #         weight[2] = 8/9
    #         weight[3] = 5/9

    #     elseif NGP == 4
    #         location[1] = -0.861136311594053
    #         location[2] = -0.339981043584856
    #         location[3] = 0.339981043584856
    #         location[3] = 0.861136311594053

    #         weight[1] = 0.347854845137454
    #         weight[2] = 0.652145154862546
    #         weight[3] = 0.652145154862546
    #         weight[4] = 0.347854845137454

    #     end
    #     return location, weight
    # end
    
    # function that generates the constitutive matrix for the material and state
    function ConstMtrx(PlaneStressOrStrain, Props)

        E = Props["E"]
        v = Props["v"]

        C = zeros((3,3))

        # select state between Plane Stress or Plane Strain

        if PlaneStressOrStrain == "PlaneStress"
            aux = (E/(1-(v^2)))
            C[1,1] = aux
            C[1,2] = v*aux
            C[2,1] = v*aux
            C[2,2] = aux
            C[3,3] = ((1-v)/2)*aux

        elseif PlaneStressOrStrain == "PlaneStrain"
            aux = (E * (1 - v) / ((1 + v) * (1 - 2 * v)))
            C[1,1] = aux
            C[1,2] = (v/(1-v))*aux
            C[2,1] = (v/(1-v))*aux
            C[2,2] = aux
            C[3,3] = ((1 - 2 * v)/(2 * (1 - v)))*aux
            
        end

        return C
    end

    # function that get the isoparamétrics elements derivates
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

        elseif N_NodesInElem == 4
            dPet = (1 + eta)*0.25
            dMet = (1 - eta)*0.25
            dPxs = (1 + csi)*0.25
            dMxs = (1 - csi)*0.25
    
            dNdcsi[1] = -dMet
            dNdcsi[2] = dMet
            dNdcsi[3] = dPet
            dNdcsi[4] = -dPet
    
            dNdeta[1] = -dMxs
            dNdeta[2] = -dPxs
            dNdeta[3] = dPxs
            dNdeta[4] = dMxs
        end

        return dNdcsi,dNdeta
    
    end

    # function that calculates the Jacobian matrix between natural and cartesian coords
    function InvJacMat(dNdcsi, dNdeta, XY_Elem)
        xelem = XY_Elem[:,1]
        yelem = XY_Elem[:,2]

        dxdcsi = LinearAlgebra.dot(dNdcsi,xelem)
        dxdeta = LinearAlgebra.dot(dNdeta, xelem)
        dydcsi = LinearAlgebra.dot(dNdcsi, yelem)
        dydeta = LinearAlgebra.dot(dNdeta, yelem)

        detJ = dxdcsi*dydeta-dydcsi*dxdeta

        InvJac = [dydeta -dydcsi;
                  -dxdeta dxdcsi]/detJ

        return InvJac,detJ
    end

    # function that transforms the natural coords derivates in cartesian derivates
    function DerivCartShapeFunc(dNdcsi, dNdeta, InvJac, N_NodesInElem)
        dNnat = zeros((2,N_NodesInElem))

        for i in 1:N_NodesInElem
            dNnat[1,i] = dNdcsi[i]
            dNnat[2,i] = dNdeta[i]
        end

        dNdcart = InvJac * dNnat

        return dNdcart
    end

    function MatDefDesloc(dNdcart, N_NodesInElem)
        B = zeros((3,N_NodesInElem*2))
        
        for i in 1:N_NodesInElem
            B[1,2*i-1] = dNdcart[1,i]
            B[2,2*i] = dNdcart[2,i]
            B[3,2*i-1] = dNdcart[2,i]
            B[3,2*i] = dNdcart[1,i]
        end
        return B
    end

    # function that returns the global stiffness matrix and the internal forces
    function Get_GlobalK(N_NodesInElem, NGP, Props, N_Elems, DoFElem, DoFNode, Restrs, N_DoF, NodesCoord, if_Plast,
        csi, eta, w, Cel, sigma_total, Connect, f_ext)
        
        K = zeros((N_DoF,N_DoF)) #stiffness matrix
        XY_Elem = zeros(N_NodesInElem, 2) # elemets' coords
        N_points = NGP*NGP # number of integration points
        f_int = zeros(N_DoF) 

        for i_elem in 1:N_Elems

            K_elem = zeros((DoFElem,DoFElem)) # element's K
            f_int_elem = zeros(DoFElem) # element's internal force
            i_sigma = (i_elem-1)*N_points

            mapvec = zeros(Int32,DoFElem)
            @simd for d in 1:DoFElem
                mapvec[d] = (Connect[i_elem][(d-1)÷DoFNode+1]-1)*DoFNode+((d-1)%DoFNode+1)
            end

            # get the element i_elem nodes coords
            @simd for i_node in 1:N_NodesInElem
                XY_Elem[i_node,:] = NodesCoord[(Connect[i_elem][i_node])]
            end

            for integ_points in 1:N_points
                # get the shape form functions derivates in natural coords
                dNdcsi, dNdeta = DerivNatShapeFunc(csi[integ_points], eta[integ_points],N_NodesInElem)
                
                #get the jacobian matrix and its determinant
                InvJac, J = InvJacMat(dNdcsi, dNdeta, XY_Elem)
                
                dNdcart = DerivCartShapeFunc(dNdcsi, dNdeta, InvJac, N_NodesInElem)

                B = MatDefDesloc(dNdcart, N_NodesInElem)
               
                # J = LinearAlgebra.det(Jac)

                Bt = transpose(B)

                if if_Plast[i_elem,integ_points] == 1

                    sigma_xx = sigma_total[1,i_sigma+integ_points]
                    sigma_yy = sigma_total[2,i_sigma+integ_points]
                    sigma_xy = sigma_total[3,i_sigma+integ_points]
                    df = [(2*sigma_xx-sigma_yy)/3, (2*sigma_yy - sigma_xx)/3, 2*sigma_xy]
                    C = Cel - (Cel * df * transpose(df) * Cel)/(transpose(df) * Cel * df)
                else
                    C = Cel
                end

                K_elem += Bt* C * B * J * Props["Thickness"] * w[integ_points]
                f_int_elem +=  Bt * sigma_total[:, i_sigma+integ_points] * J * Props["Thickness"] * w[integ_points]
                
            end

            for i_dof in 1:DoFElem
                @simd for j_dof in 1:DoFElem
                    K[mapvec[i_dof],mapvec[j_dof]] += K_elem[i_dof,j_dof]
                end
                f_int[mapvec[i_dof]] += f_int_elem[i_dof]
            end
        end
        
        for i_rest in Restrs
            # restricted Node
            rest_node = i_rest[1]
            # if there is a restriction on i_dof
            for i_dof in 1:DoFNode
                if i_rest[i_dof+1] == 1
                    @simd for i in 1:N_DoF
                        K[(rest_node-1)*DoFNode+i_dof,i] = zero(Float64)
                    end
                    K[(rest_node-1)*DoFNode+i_dof,(rest_node-1)*DoFNode+i_dof] = 1
                    f_int[(rest_node-1)*DoFNode+i_dof] = f_ext[(rest_node-1)*DoFNode+i_dof]
                    # aux = maximum(K[(rest_node-1)*DoFNode+i_dof,:])
                    # K[(rest_node-1)*DoFNode+i_dof,(rest_node-1)*DoFNode+i_dof] = aux*10^16
                    # f_int[(rest_node-1)*DoFNode+i_dof] = f_ext[(rest_node-1)*DoFNode+i_dof]
                end
            end
        end
        return K, f_int
    end

    # function that calculates the Finite Elements method in elastoplastics conditions
    function FEM_Ep(N_DoF, N_Steps, N_Elems, Connect, N_NodesInElem, NGP, NodesCoord, DoFElem, Props, PlaneStressOrStrain, DoFNode, Forces, Restrs)
        
        # Setting matrix for function's use

        F = zeros(Float64,N_DoF) # total force for node
        f_ext = zeros(Float64,N_DoF) # external force for node in each step
        sigma_total = zeros(Float64,(3, NGP * NGP * N_Elems)) # stress tensor
        D = zeros(Float64,(N_DoF)) # total displacement
        dD = zeros(Float64,(N_DoF)) # displacement for step
        dD_Elem = zeros(Float64,(DoFElem)) # elements' displacement for step
        XY_Elem = zeros(Float64,(N_NodesInElem, 2)) # elemets' coords
        if_Plast = zeros(Bool,(N_Elems,NGP*NGP)) # state of each gauss point for elem
        tolD = 0.001 # displacement increment tolerance
        J2Elem = zeros(Float64,N_Elems) # von Mises yield criterion
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

                i_dof = DoFNode * Forces[i][2] - 1
                F[i_dof] += F_Node
                i_dof = DoFNode * Forces[i][3] - 1
                F[i_dof] += F_Node
                
            end

            # if there is force in Y, add to the vector F
            if Forces[i][5] != 0
                F_Node = (Forces[i][5] * r)/2 # half force for each node

                i_dof = DoFNode * Forces[i][2]
                F[i_dof] += F_Node
                i_dof = DoFNode * Forces[i][3]
                F[i_dof] += F_Node
                
            end

        end

        # Getting Gauss points locations and weights
        (GP_locations, GP_weights) = gausslegendre(NGP)
        k2 = 1
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

        for i_step in 1:N_Steps
            f_ext += f_incr
            count = 0
            maxdD = 11

            # start Newthon-Raphson method to calculate the displacement
            while maxdD > tolD
                
                # calculate the stiffness matrix K and the internal forces f_int
                (K, f_int) = Get_GlobalK(N_NodesInElem, NGP, Props, N_Elems, DoFElem, DoFNode, Restrs, N_DoF, NodesCoord, if_Plast,
                csi, eta, w, Cel, sigma_total, Connect, f_ext)
                b = f_ext - f_int
                dD = inv(K)*b

                maxdD = 0

                N_points = NGP*NGP # number of integration points

                for i_elem in 1:N_Elems

                    i_sigma = (i_elem-1)*N_points


                    @simd for i_dof in 1:DoFElem
                        dD_Elem[i_dof] = dD[(Connect[i_elem][(i_dof-1)÷DoFNode+1]-1)*DoFNode+((i_dof-1)%DoFNode+1)]
                    end
        
                    # get the element i_elem nodes coords
                    @simd for i_node in 1:N_NodesInElem
                        XY_Elem[i_node,:] = NodesCoord[(Connect[i_elem][i_node])]
                    end
        
                    for integ_points in 1:N_points
                        # get the shape form functions derivates in natural coords
                        dNdcsi, dNdeta = DerivNatShapeFunc(csi[integ_points], eta[integ_points],N_NodesInElem)
        
                        #get the jacobian matrix and its determinant
                        InvJac, J = InvJacMat(dNdcsi, dNdeta, XY_Elem)
        
                        dNdcart = DerivCartShapeFunc(dNdcsi, dNdeta, InvJac, N_NodesInElem)
        
                        B = MatDefDesloc(dNdcart, N_NodesInElem)
        
                        if if_Plast[i_elem,integ_points] == 1
        
                            sigma_xx = sigma_total[1,i_sigma+integ_points]
                            sigma_yy = sigma_total[2,i_sigma+integ_points]
                            sigma_xy = sigma_total[3,i_sigma+integ_points]
                            df = [(2*sigma_xx-sigma_yy)/3, (2*sigma_yy - sigma_xx)/3, 2*sigma_xy]
                            C = Cel - (Cel * df * transpose(df) * Cel)/(transpose(df) * Cel * df)
                        else
                            C = Cel
                        end

                        d_sigma = C * B * dD_Elem

                        # sigma_xx = sigma_total[1,i_sigma+integ_points] + d_sigma[1]
                        # sigma_yy = sigma_total[2,i_sigma+integ_points] + d_sigma[2]
                        # sigma_xy = sigma_total[3,i_sigma+integ_points] + d_sigma[3]
                        sigma_xx,sigma_yy,sigma_xy = sigma_total[:,i_sigma+integ_points] + d_sigma

                        # von Mises criterion
                        J2 = 1/6 * ((sigma_xx-sigma_yy)^2 + sigma_yy^2 + sigma_xx^2) + sigma_xy^2
                        f_yield = J2 - ky^2
                        J2Elem[i_elem] += sqrt(J2)

                        if f_yield > 0
                            sigma_trial = sigma_total[:,i_sigma+integ_points]
                            k_trial = sqrt(J2)
                            sigma_mean = [(sigma_xx + sigma_yy) / 3, (sigma_xx + sigma_yy) / 3, 0]
                            dev = sigma_trial-sigma_mean
                            aux = ky/k_trial*dev
                            sigma_total[:, i_sigma+integ_points] = aux + sigma_mean
                            if_Plast[i_elem,integ_points] = 1
                        else
                            sigma_total[:,i_sigma + integ_points] += d_sigma
                            if_Plast[i_elem,integ_points] = 0
                        end
                    end
                end
                
                D += dD
                maxdD = max(maximum(dD),minimum(dD))
                count +=1
                if count>20
                    println("Possible fracture, ending loops")
                    exit()
                end
            end
        end


        return D, sigma_total, J2Elem
        end
    end