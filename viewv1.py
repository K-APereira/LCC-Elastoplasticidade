# Importação da biblioteca matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Definição da função "viewmesh" com dados de entrada das funções "readmesh", "nodeindex" e "Coord"
def viewmesh(coords, connect, Nnodes):

    # Configurações de plotagem
    color = 'cornflowerblue'
    linestyle = '-'
    bbox = {'fc': '0.8', 'pad': 0.05}
    props = {'ha': 'left', 'va': 'bottom', 'bbox': bbox}
    count = 1
    fontsize = 0
    #fontsize = 6

    for e in range(len(connect[0])):
        x = []
        y = []
        xm = 0
        ym = 0

        for n in range(len(connect) - (Nnodes - 4)):
            x.append(coords[int(connect[n, e]), 0])
            y.append(coords[int(connect[n, e]), 1])
            xm += coords[int(connect[n, e]), 0]
            ym += coords[int(connect[n, e]), 1]
        count += 1
        x.append(coords[int(connect[n - 3, e]), 0])
        y.append(coords[int(connect[n - 3, e]), 1])
        xm += coords[int(connect[n - 3, e]), 0]
        ym += coords[int(connect[n - 3, e]), 1]
        plt.plot(x, y, linestyle=linestyle, color=color, linewidth=1)
    color = 'black'
    plt.plot(coords[:, 0], coords[:, 1], '.', color=color)
    xmin = min(coords[:, 0])
    xmax = max(coords[:, 0])
    ymin = min(coords[:, 1])
    ymax = max(coords[:, 1])
    count = 1
    for n in coords:
        plt.text(n[0] + 0.001, n[1] + 0.001, str(count), fontsize=fontsize)
        count += 1

    # Configurações para plotagem com diversas escalas de eixos
    plt.grid(True)
    plt.xlabel('x - axis')
    plt.ylabel('y - axis')
    plt.xlim(xmin - 2, xmax + 2)
    plt.ylim(ymin - 2, ymax + 6)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.draw()
    plt.show()

# Definição da função "viewmesh" com dados de entrada das funções "readmesh", "nodeindex", "Coord" e (...)
def viewdeformedmesh(coords, connect, disp, scale, Nnodes):

    # Configurações de plotagem
    color = 'cornflowerblue'
    linestyle = '-'
    bbox = {'fc': '0.8', 'pad': 2}
    props = {'ha': 'left', 'va': 'bottom', 'bbox': bbox}
    count = 1
    fontsize = 0
    #fontsize = 6
    for e in range(len(connect[0])):
        x = []
        y = []
        xm = 0
        ym = 0

        for n in range(len(connect) - (Nnodes - 4)):
            x.append(coords[int(connect[n, e]), 0])
            y.append(coords[int(connect[n, e]), 1])
            xm += coords[int(connect[n, e]), 0]
            ym += coords[int(connect[n, e]), 1]
            plt.plot(x, y, linestyle=linestyle, color=color, linewidth=2)
        count += 1
        x.append(coords[int(connect[n - 3, e]), 0])
        y.append(coords[int(connect[n - 3, e]), 1])
        xm += coords[int(connect[n - 3, e]), 0]
        ym += coords[int(connect[n - 3, e]), 1]
        plt.plot(x, y, linestyle=linestyle, color=color, linewidth=1)
    '''color = 'black'
    plt.plot(coords[:, 0], coords[:, 1], 'o', color=color)
    count = 1
    for n in coords:
        plt.text(n[0] + 0.001, n[1] + 0.001, str(count), fontsize=fontsize)
        count += 1'''

    color = 'red'
    linestyle = ':'
    count = 1
    for e in range(len(connect[0])):
        x = []
        y = []
        for n in range(len(connect) - (Nnodes - 4)):
            x.append(coords[int(connect[n, e]), 0] + scale * disp[2 * int(connect[n, e]), 0])
            y.append(coords[int(connect[n, e]), 1] + scale * disp[2 * int(connect[n, e]) + 1, 0])
            xm += coords[int(connect[n, e]), 0]
            ym += coords[int(connect[n, e]), 1]
            plt.plot(x, y, linestyle=linestyle, color=color, linewidth=2)
        count += 1

        x.append(coords[int(connect[n - 3, e]), 0] + scale * disp[2 * int(connect[n - 3, e]), 0])
        y.append(coords[int(connect[n - 3, e]), 1] + scale * disp[2 * int(connect[n - 3, e]) + 1, 0])
        xm += coords[int(connect[n - 3, e]), 0]
        ym += coords[int(connect[n - 3, e]), 1]
        plt.plot(x, y, linestyle=linestyle, color=color, linewidth=1)

    # Ajuste da região do gráfico para limites X e Y de acordo com cada geometria
    dispx = disp[0:len(disp):2, 0]
    dispy = disp[1:len(disp):2, ]
    xlist = []
    ylist = []
    for i in range(len(coords)):
        xlist.append(coords[i, 0] + scale * dispx[i, 0])
        ylist.append(coords[i, 1] + scale * dispy[i, 0])
    xmin = min(xlist)
    xmax = max(xlist)
    ymin = min(ylist)
    ymax = max(ylist)

    plt.grid(True)
    plt.xlabel('x - axis')
    plt.ylabel('y - axis')
    plt.xlim(xmin - 2, xmax + 2)
    plt.ylim(ymin - 2, ymax + 2)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.draw()
    plt.show()

def viewcolormesh(coords, connect, Nnodes, j2Elem):
    # Configurações de plotagem
    color = 'cornflowerblue'
    linestyle = '-'
    bbox = {'fc': '0.8', 'pad': 0.05}
    props = {'ha': 'left', 'va': 'bottom', 'bbox': bbox}
    count = 1
    fontsize = 0
    fontsize = 6    

    # Normalize
    norm = mcolors.Normalize(vmin=j2Elem.min(), vmax=j2Elem.max())

    for e in range(len(connect[0])):
        x = []
        y = []
        xm = 0
        ym = 0

        for n in range(len(connect) - (Nnodes - 4)):
            x.append(coords[int(connect[n, e]), 0])
            y.append(coords[int(connect[n, e]), 1])
            xm += coords[int(connect[n, e]), 0]
            ym += coords[int(connect[n, e]), 1]
        count += 1
        x.append(coords[int(connect[n - 3, e]), 0])
        y.append(coords[int(connect[n - 3, e]), 1])
        xm += coords[int(connect[n - 3, e]), 0]
        ym += coords[int(connect[n - 3, e]), 1]
        plt.plot(x, y, linestyle=linestyle, color=color, linewidth=1)
        plt.fill(x, y, color=plt.cm.jet(norm(j2Elem[e])), edgecolor='black')
    
    color = 'black'
    plt.plot(coords[:, 0], coords[:, 1], '.', color=color)
    xmin = min(coords[:, 0])
    xmax = max(coords[:, 0])
    ymin = min(coords[:, 1])
    ymax = max(coords[:, 1])
    count = 1
    for n in coords:
        plt.text(n[0] + 0.001, n[1] + 0.001, str(count), fontsize=fontsize)
        count += 1

    # Configurações para plotagem com diversas escalas de eixos
    plt.grid(True)
    plt.xlabel('x - axis')
    plt.ylabel('y - axis')
    plt.xlim(xmin - 2, xmax + 2)
    plt.ylim(ymin - 2, ymax + 6)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.draw()
    plt.show()