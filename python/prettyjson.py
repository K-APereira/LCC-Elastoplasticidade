def fwriteKeyVals(data, f, indent=0):
    if isinstance(data, list):
        f.write( "\n" + "    " * indent + "[" )  
        for i in range(len(data) ):  
            if ( i == 0):
                f.write( "[" )
            else:
                f.write( "    " * indent + " [" )
            for j in range(len(data[0])):  
                f.write( "%3d" % data[i][j] )
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

import json
f = open( r"Julia\Exemplo2.json" , "r" )
f_write = open(r"Julia\teste.json","w")
x = json.load(f)
fwriteKeyVals(x, f_write,indent=1)