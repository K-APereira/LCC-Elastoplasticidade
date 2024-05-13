import numpy as np

def fwriteKeyVals(data, f, indent=0):
    if isinstance(data, list) or isinstance(data, np.ndarray):
        f.write( "\n" + "    " * indent + "[" )  
        for i in range(len(data) ):  
            if ( i == 0):
                f.write( "[" )
            else:
                f.write( "    " * indent + " [" )
            for j in range(len(data[0])):  
                f.write( "%f" % data[i][j] )
                f.write( "," ) if j != len(data[0])-1 else (f.write( "]," ) if i != len(data)-1 else f.write( "]" ))  
            f.write( "\n" ) if i != len(data)-1 else f.write( "]" )  
    elif isinstance(data, dict):
        f.write( "\n" + "    " * indent + "{" )
        for k, v in data.items():
            f.write( "\n" + "    " * indent + "\"" + k + "\"" + ": ")
            fwriteKeyVals(v, f, indent + 1)
            if( list(data.keys())[-1] != k):
                 f.write( "," )
        f.write( "\n" + "    " * indent + "}" )
    else:
        # f.write("\"" + str(data) + "\"")
        if type(data) == str:
            f.write("\"" + data + "\"")
        else:
            print(data,file=f, end="")

# import json
# f = open( r"Julia\Exemplo2.json" , "r" )
# f_write = open(r"Julia\teste.json","w")
# x = json.load(f)
# fwriteKeyVals(x, f_write,indent=1)

######################################
from readmeshv2 import *
[NODES,coord_nodes,num_restrs,restrs,num_isomat,props,num_thick,elem_nodes,NELE,connect,Dofnode\
        ,dofelem,NDoF,forces,nstep,planestress,NGP] = readmesh(r"python\Exemplo 1.1 - 100 elementos.txt", "python/Model_Tk.nf")

if planestress == 1:
    pls = "PlaneStress"
else:
    pls = "PlaneStrain"
from prettyjson import *
import json
data = {
    "Props":{
        "E": props[0][0],
        "v": props[0][1],
        "Thickness": props[0][2],
        "YieldStress": props[0][3]
    },

    "PlaneStressOrStrain": pls,

    "DoFNode": Dofnode,

    "NodesCoord":coord_nodes[:,1:3],

    "Elements":connect.T+1,

    "Restrs":restrs+[1,0,0],

    "Forces":forces+[1,1,1,0,0],

    "N_Steps": nstep,

    "NGP": NGP
    
}
file = open(r"python\Exemplo 1.1 - 100 elementos.json","w")
fwriteKeyVals(data,file,1)
file.close()
######################################