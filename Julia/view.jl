module view
    using GLMakie
   


    function viewdeformedmesh(f,NodesCoord, Connect, D, Scale)
        Axis(f[1,1])
        edges=Set{Tuple{Int,Int}}()

        for e in Connect 
            lengthElem=length(e)

            for iN in 1:lengthElem-1 
                push!(edges,Tuple([min(e[iN],e[iN+1]),max(e[iN],e[iN+1])]))
            end
            push!(edges,Tuple([min(e[end],e[1]),max(e[end],e[1])]))
        end

        for e in edges
            lines!([NodesCoord[e[1]][1],NodesCoord[e[2]][1]],   [NodesCoord[e[1]][2],NodesCoord[e[2]][2]],color=:blue)

            d1=[D[2*e[1]-1],D[2*e[1]]] .* Scale
            d2=[D[2*e[2]-1],D[2*e[2]]] .* Scale

            x1=NodesCoord[e[1]][1] + d1[1]
            x2=NodesCoord[e[2]][1] + d2[1]

            y1=NodesCoord[e[1]][2] + d1[2]
            y2=NodesCoord[e[2]][2] + d2[2]

            lines!([x1,x2],[y1,y2], linestyle = :dot, color=:red)
        end
    end

    function viewcolormesh(f,NodesCoord, Connect, j2Elem)
        Axis(f[1,2])
        minimo=minimum(j2Elem)
        maximo=maximum(j2Elem)
        
        for (index, value) in enumerate(Connect)
            edges=Point2f[]
            for p in value 
                push!(edges,Point2f(NodesCoord[p][1],NodesCoord[p][2]))
            end
            poly!(edges, color = j2Elem[index],colormap=:jet,colorrange=(minimo,maximo) , strokecolor = :black, strokewidth = 2)

        end

        
    end

    function viewMesh(NodesCoord, Connect, D, Scale, j2Elem)
        GLMakie.activate!(title = "Custom title", fxaa = false)
        f=Figure()

        viewdeformedmesh(f,NodesCoord, Connect, D, Scale)

        viewcolormesh(f,NodesCoord, Connect, j2Elem)

        f
    end
end