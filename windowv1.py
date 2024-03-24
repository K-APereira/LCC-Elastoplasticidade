from tkinter import *
from readmeshv1 import *
from Coordv1 import *
from viewv1 import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

def clicar():
    L = entrada_1.get()
    H = entrada_2.get()
    t = entrada_3.get()
    E = entrada_4.get()
    fy = entrada_5.get()
    v = entrada_6.get()
    NEX = entrada_7.get()
    NEY = entrada_8.get()
    NODES = var_1.get()
    RX_B1 = var_2.get()
    RX_B2 = var_3.get()
    RX_B3 = var_4.get()
    RX_B4 = var_5.get()
    RY_B1 = var_6.get()
    RY_B2 = var_7.get()
    RY_B3 = var_8.get()
    RY_B4 = var_9.get()
    FX_B1 = var_10.get()
    FX_B2 = var_11.get()
    FX_B3 = var_12.get()
    FX_B4 = var_13.get()
    FY_B1 = var_14.get()
    FY_B2 = var_15.get()
    FY_B3 = var_16.get()
    FY_B4 = var_17.get()
    MAGNITUDE = entrada_9.get()
    NSTEP = entrada_10.get()
    NGP = var_18.get()
    PLANESTRESS = var_19.get()

    fileModel = open("Model_TK.nf", 'w')

    fileModel.write("%TITLE\n")
    fileModel.write("1D MODEL\n")
    fileModel.write("\n")

    fileModel.write("%GEOMETRY\n")
    fileModel.write("%L - mm\n")
    fileModel.write(L)
    fileModel.write("\n")

    fileModel.write("%H - mm\n")
    fileModel.write(H)
    fileModel.write("\n")

    fileModel.write("\n%PROPERTIES\n")
    fileModel.write("1\n")
    fileModel.write(E)
    fileModel.write(" ")
    fileModel.write(v)
    fileModel.write(" ")
    fileModel.write(t)
    fileModel.write(" ")
    fileModel.write(fy)
    fileModel.write(" ")
    fileModel.write("\n")
    fileModel.write("% MPa  -   mm MPa\n")

    fileModel.write("\n%MESH\n")
    fileModel.write("%NODES\n")
    fileModel.write(NODES)
    fileModel.write("\n")
    fileModel.write("%NEX\n")
    fileModel.write(NEX)
    fileModel.write("\n")
    fileModel.write("%NEY\n")
    fileModel.write(NEY)
    fileModel.write("\n")

    fileModel.write("\n%RESTRAINS (0 - livre, 1 - fixo)\n")
    fileModel.write("4\n")
    fileModel.write("1")
    fileModel.write("  ")
    fileModel.write(RX_B1)
    fileModel.write("  ")
    fileModel.write(RY_B1)
    fileModel.write("\n")
    fileModel.write("2")
    fileModel.write("  ")
    fileModel.write(RX_B2)
    fileModel.write("  ")
    fileModel.write(RY_B2)
    fileModel.write("\n")
    fileModel.write("3")
    fileModel.write("  ")
    fileModel.write(RX_B3)
    fileModel.write("  ")
    fileModel.write(RY_B3)
    fileModel.write("\n")
    fileModel.write("4")
    fileModel.write("  ")
    fileModel.write(RX_B4)
    fileModel.write("  ")
    fileModel.write(RY_B4)
    fileModel.write("\n")
    fileModel.write("%  X  Y\n")

    fileModel.write("\n%FORCES (0 -  livre, 1 - axial)\n")
    fileModel.write("4\n")
    fileModel.write("1")
    fileModel.write("  ")
    fileModel.write(FX_B1)
    fileModel.write("  ")
    fileModel.write(FY_B1)
    fileModel.write("\n")
    fileModel.write("2")
    fileModel.write("  ")
    fileModel.write(FX_B2)
    fileModel.write("  ")
    fileModel.write(FY_B2)
    fileModel.write("\n")
    fileModel.write("3")
    fileModel.write("  ")
    fileModel.write(FX_B3)
    fileModel.write("  ")
    fileModel.write(FY_B3)
    fileModel.write("\n")
    fileModel.write("4")
    fileModel.write("  ")
    fileModel.write(FX_B4)
    fileModel.write("  ")
    fileModel.write(FY_B4)
    fileModel.write("\n")
    fileModel.write("%  X  Y\n")

    fileModel.write("\n%MAGNITUDE\n")
    fileModel.write(MAGNITUDE)
    fileModel.write("\n")
    fileModel.write("% N (NEGATIVE - COMPRESSION; POSITIVE - TRACTION)\n")

    fileModel.write("\n%NSTEP\n")
    fileModel.write(NSTEP)
    fileModel.write("\n")

    fileModel.write("\n%NGP\n")
    fileModel.write(NGP)
    fileModel.write("\n")

    fileModel.write("\n%PLANESTRESS\n")
    fileModel.write(PLANESTRESS)
    fileModel.write("\n")
    fileModel.write("% 1 = Estado plano de tensoes\n")
    fileModel.write("% 2 = Estado plano de deformacoes\n")

    fileModel.close()

    [L, H, props, NEY, NEX, NELE, Dofnode, restrs, forces, P, planestress, Nnodes, Dx, Dy, NGP, dofelem,
     nstep] = readmesh("Model_Tk.nf")
    [NODES, connect, BC1node, BC2node, BC3node, BC4node, NDoF] = nodeindex(NEX, NEY, NELE, Nnodes, Dofnode)
    assmtrx = dofdrive(Nnodes, NELE, connect, dofelem)
    [X, Y, XY] = Coord(NODES, Dx, Dy, NEX, NEY, connect, Nnodes)
    plot_v1(XY, connect, Nnodes)

def plot_v1(coords, connect, Nnodes):

    # Configurações de plotagem
    color = 'cornflowerblue'
    linestyle = '-'
    bbox = {'fc': '0.8', 'pad': 0.05}
    props = {'ha': 'left', 'va': 'bottom', 'bbox': bbox}
    count = 1
    fontsize = 6
    fig = Figure(figsize=(4, 4))
    a = fig.add_subplot()

    for e in range(len(connect[0])):
        x = []
        y = []
        xm = 0
        ym = 0

        for n in range(len(connect) - (Nnodes - 4)):
            x.append(coords[connect[n, e], 0])
            y.append(coords[connect[n, e], 1])
            xm += coords[connect[n, e], 0]
            ym += coords[connect[n, e], 1]
        count += 1
        x.append(coords[connect[n - 3, e], 0])
        y.append(coords[connect[n - 3, e], 1])
        xm += coords[connect[n - 3, e], 0]
        ym += coords[connect[n - 3, e], 1]
        a.plot(x, y, linestyle=linestyle, color=color, linewidth=1)

    color = 'black'

    a.plot(coords[:, 0], coords[:, 1], '.', color=color)

    xmin = min(coords[:, 0])
    xmax = max(coords[:, 0])
    ymin = min(coords[:, 1])
    ymax = max(coords[:, 1])
    count = 1

    for n in coords:
        a.text(n[0] + 0.001, n[1] + 0.001, str(count), fontsize=fontsize)
        count += 1

    # Configurações para plotagem com diversas escalas de eixos
    a.grid(True)
    a.set_xlabel("Eixo X")
    a.set_ylabel("Eixo Y")
    a.set_xlim(xmin - 2, xmax + 2)
    a.set_ylim(ymin - 2, ymax + 6)
    a.set_title("Modelo Discretizado", fontsize=15)

    canvas = FigureCanvasTkAgg(fig)
    canvas.get_tk_widget().grid(row=0,column=1)
    canvas.draw()

def janela_v1():

    #Criação da interface
    app = Tk()
    app.title("Análise de chapas metálicas com comportamento elastoplástico")
    app.geometry('2000x650')

    # Interface de Pré-processamento
    # Criando o quadro e legenda para o quadro do Pré-Processamento
    J1 = LabelFrame(app)
    J1.grid(row=0, column=0, padx=10, pady=10)

    label_1 = Label(J1, text='Pré-Processamento', font=('Verdana', 20))
    label_1.grid(row=0, column=0, columnspan=3)

    # Criando o quadro e legendas para o quadro de entrada de geometria da chapa do Pré-Processamento
    J2 = LabelFrame(J1)
    J2.grid(row=1, column=0, padx=10, pady=10)

    label_2 = Label(J2, text='Geometria da Chapa', font=('Verdana', 15))
    label_2.grid(row=0, column=0, ipady=10, columnspan=6)

    # Criando as legendas e entradas para os dados de entrada de geometria do Pré-Processamento
    label_3 = Label(J2, text='L [mm]', font=('Verdana', 12)) # Comprimento
    label_3.grid(row=1, column=0, ipady=10)

    global entrada_1
    entrada_1 = Entry(J2, width=5, font=('Verdana', 10), justify='center')
    entrada_1.grid(row=1, column=1, padx=10)

    label_4 = Label(J2, text='H [mm]', font=('Verdana', 12)) # Altura
    label_4.grid(row=1, column=2, ipady=10)

    global entrada_2
    entrada_2 = Entry(J2, width=5, font=('Verdana', 10), justify='center')
    entrada_2.grid(row=1, column=3, padx=10)

    label_5 = Label(J2, text='t [mm]', font=('Verdana', 12)) # Espessura
    label_5.grid(row=1, column=4, ipady=10)

    global entrada_3
    entrada_3 = Entry(J2, width=5, font=('Verdana', 10), justify='center')
    entrada_3.grid(row=1, column=5, padx=10)

    # Criando o quadro e legendas para o quadro de entrada de propriedades do Pré-Processamento
    J3 = LabelFrame(J1)
    J3.grid(row=1, column=1, padx=10, pady=10)

    label_6 = Label(J3, text='Propriedades do Material', font=('Verdana', 15))
    label_6.grid(row=0, column=0, ipady=10, columnspan=6)

    # Criando as legendas e entradas para os dados de entrada de propriedades do Pré-Processamento
    label_7 = Label(J3, text='E [MPa]', font=('Verdana', 12)) # Módulo de elasticidade
    label_7.grid(row=1, column=0, ipady=10)

    global entrada_4
    entrada_4 = Entry(J3, width=10, font=('Verdana', 10), justify='center')
    entrada_4.grid(row=1, column=1, padx=10)

    label_8 = Label(J3, text='fy [MPa]', font=('Verdana', 12)) # Tensão de escoamento
    label_8.grid(row=1, column=2, ipady=10)

    global entrada_5
    entrada_5 = Entry(J3, width=5, font=('Verdana', 10), justify='center')
    entrada_5.grid(row=1, column=3, padx=10)

    label_9 = Label(J3, text='v', font=('Verdana', 12)) # Coeficiente de Poisson
    label_9.grid(row=1, column=4, ipady=10)

    global entrada_6
    entrada_6 = Entry(J3, width=5, font=('Verdana', 10), justify='center')
    entrada_6.grid(row=1, column=5, padx=10)

    # Criando o quadro e legendas para o quadro de entrada da malha de EF do Pré-Processamento
    J4 = LabelFrame(J1)
    J4.grid(row=1, column=2, padx=10, pady=10)

    label_10 = Label(J4, text='Propriedades da Malha de EF', font=('Verdana', 15))
    label_10.grid(row=0, column=0, ipady=10, columnspan=6)

    # Criando as legendas e entradas para os dados de entrada da malha de EF do Pré-Processamento
    label_11 = Label(J4, text='Elementos em X', font=('Verdana', 12)) # Elementos em X
    label_11.grid(row=1, column=0, ipady=10)

    global entrada_7
    entrada_7 = Entry(J4, width=5, font=('Verdana', 10), justify='center')
    entrada_7.grid(row=1, column=1, padx=10)

    label_12 = Label(J4, text='Elementos em Y', font=('Verdana', 12)) # Elementos em Y
    label_12.grid(row=1, column=2, ipady=10)

    global entrada_8
    entrada_8 = Entry(J4, width=5, font=('Verdana', 10), justify='center')
    entrada_8.grid(row=1, column=3, padx=10)

    label_13 = Label(J4, text='Nº Nós por Elemento', font=('Verdana', 12)) # Elementos em Y
    label_13.grid(row=1, column=4, ipady=10)

    global var_1
    var_1 = StringVar()
    var_1.set('4')
    botao_1 = Radiobutton(J4, text='4', font=('Verdana', 12), variable=var_1, value='4')
    botao_1.grid(row=1, column=5)
    botao_2 = Radiobutton(J4, text='8', font=('Verdana', 12), variable=var_1, value='8')
    botao_2.grid(row=2, column=5)

    # Criando o quadro e legendas para o quadro de entrada das CC do Pré-Processamento
    J5 = LabelFrame(J1)
    J5.grid(row=2, column=0, padx=10, pady=10, columnspan=2)

    label_14 = Label(J5, text='Condições de Contorno', font=('Verdana', 15))
    label_14.grid(row=0, column=0, ipady=10, columnspan=14)

    # Criando legendas para visualização das chapa e suas bordas para entrada das CC do Pré-Processamento
    label_15 = Label(J5, text='Borda 4', font=('Verdana', 12)) # Borda 4
    label_15.grid(row=5, column=0)

    label_16 = Label(J5, text='  Borda 3', font=('Verdana', 12)) # Borda 3
    label_16.grid(row=2, column=1, columnspan=4)

    label_17 = Label(J5, text='-', font=('Verdana', 14)) # _
    label_17.grid(row=3, column=1)

    label_18 = Label(J5, text='-', font=('Verdana', 14)) # _
    label_18.grid(row=3, column=2)

    label_19 = Label(J5, text='-', font=('Verdana', 14)) # _
    label_19.grid(row=3, column=3)

    label_20 = Label(J5, text='-', font=('Verdana', 14)) # _
    label_20.grid(row=3, column=4)

    label_21 = Label(J5, text='|', font=('Verdana', 14)) # |
    label_21.grid(row=4, column=1)

    label_22 = Label(J5, text='|', font=('Verdana', 14)) # |
    label_22.grid(row=5, column=1)

    label_23 = Label(J5, text='|', font=('Verdana', 14)) # |
    label_23.grid(row=6, column=1)

    label_24 = Label(J5, text='-', font=('Verdana', 14)) # _
    label_24.grid(row=7, column=1)

    label_25 = Label(J5, text='-', font=('Verdana', 14)) # _
    label_25.grid(row=7, column=2)

    label_26 = Label(J5, text='-', font=('Verdana', 14)) # _
    label_26.grid(row=7, column=3)

    label_27 = Label(J5, text='-', font=('Verdana', 14)) # _
    label_27.grid(row=7, column=4)

    label_28 = Label(J5, text='|', font=('Verdana', 14)) # |
    label_28.grid(row=4, column=4)

    label_29 = Label(J5, text='|', font=('Verdana', 14)) # |
    label_29.grid(row=5, column=4)

    label_30 = Label(J5, text='|', font=('Verdana', 14)) # |
    label_30.grid(row=6, column=4)

    label_31 = Label(J5, text='Borda 2', font=('Verdana', 12)) # Borda 2
    label_31.grid(row=5, column=5)

    label_32 = Label(J5, text='  Borda 1', font=('Verdana', 12)) # Borda 1
    label_32.grid(row=8, column=1, columnspan=4)

    # Criando as legendas e checkboxs para os dados de entrada da malha de EF do Pré-Processamento

    label_33 = Label(J5, text='Apoios em X:', font=('Verdana', 12)) # Apoios em X
    label_33.grid(row=1, column=6, columnspan=2, pady=10)

    global var_2
    var_2 = StringVar()
    var_2.set('0')
    botao_3 = Checkbutton(J5, text='Borda 1', font=('Verdana', 10), variable=var_2) # Apoios em X - Borda 1
    botao_3.grid(row=2, column=6)

    global var_3
    var_3 = StringVar()
    var_3.set('0')
    botao_4 = Checkbutton(J5, text='Borda 2', font=('Verdana', 10), variable=var_3) # Apoios em X - Borda 2
    botao_4.grid(row=4, column=6)

    global var_4
    var_4 = StringVar()
    var_4.set('0')
    botao_5 = Checkbutton(J5, text='Borda 3', font=('Verdana', 10), variable=var_4) # Apoios em X - Borda 3
    botao_5.grid(row=6, column=6)

    global var_5
    var_5 = StringVar()
    var_5.set('0')
    botao_6 = Checkbutton(J5, text='Borda 4', font=('Verdana', 10), variable=var_5) # Apoios em X - Borda 4
    botao_6.grid(row=8, column=6)

    label_34 = Label(J5, text='Apoios em Y:', font=('Verdana', 12)) # Apoios em Y
    label_34.grid(row=1, column=8, columnspan=2, padx=10)

    global var_6
    var_6 = StringVar()
    var_6.set('0')
    botao_7 = Checkbutton(J5, text='Borda 1', font=('Verdana', 10), variable=var_6) # Apoios em Y - Borda 1
    botao_7.grid(row=2, column=8)

    global var_7
    var_7 = StringVar()
    var_7.set('0')
    botao_8 = Checkbutton(J5, text='Borda 2', font=('Verdana', 10), variable=var_7) # Apoios em Y - Borda 2
    botao_8.grid(row=4, column=8)

    global var_8
    var_8 = StringVar()
    var_8.set('0')
    botao_9 = Checkbutton(J5, text='Borda 3', font=('Verdana', 10), variable=var_8) # Apoios em Y - Borda 3
    botao_9.grid(row=6, column=8)

    global var_9
    var_9 = StringVar()
    var_9.set('0')
    botao_10 = Checkbutton(J5, text='Borda 4', font=('Verdana', 10), variable=var_9) # Apoios em Y - Borda 4
    botao_10.grid(row=8, column=8)

    label_35 = Label(J5, text='Forças em X:', font=('Verdana', 12)) # Forças em X
    label_35.grid(row=1, column=10, columnspan=2, padx=10)

    global var_10
    var_10 = StringVar()
    var_10.set('0')
    botao_11 = Checkbutton(J5, text='Borda 1', font=('Verdana', 10), variable=var_10) # Forças em X - Borda 1
    botao_11.grid(row=2, column=10)

    global var_11
    var_11 = StringVar()
    var_11.set('0')
    botao_12 = Checkbutton(J5, text='Borda 2', font=('Verdana', 10), variable=var_11) # Forças em X - Borda 2
    botao_12.grid(row=4, column=10)

    global var_12
    var_12 = StringVar()
    var_12.set('0')
    botao_13 = Checkbutton(J5, text='Borda 3', font=('Verdana', 10), variable=var_12) # Forças em X - Borda 3
    botao_13.grid(row=6, column=10)

    global var_13
    var_13 = StringVar()
    var_13.set('0')
    botao_14 = Checkbutton(J5, text='Borda 4', font=('Verdana', 10), variable=var_13) # Forças em X - Borda 4
    botao_14.grid(row=8, column=10)

    label_36 = Label(J5, text='Forças em Y:', font=('Verdana', 12)) # Forças em Y
    label_36.grid(row=1, column=12, columnspan=2, padx=5)

    global var_14
    var_14 = StringVar()
    var_14.set('0')
    botao_15 = Checkbutton(J5, text='Borda 1', font=('Verdana', 10), variable=var_14) # Forças em Y - Borda 1
    botao_15.grid(row=2, column=12)

    global var_15
    var_15 = StringVar()
    var_15.set('0')
    botao_16 = Checkbutton(J5, text='Borda 2', font=('Verdana', 10), variable=var_15) # Forças em Y - Borda 2
    botao_16.grid(row=4, column=12)

    global var_16
    var_16 = StringVar()
    var_16.set('0')
    botao_17 = Checkbutton(J5, text='Borda 3', font=('Verdana', 10), variable=var_16) # Forças em Y - Borda 3
    botao_17.grid(row=6, column=12)

    global var_17
    var_17 = StringVar()
    var_17.set('0')
    botao_18 = Checkbutton(J5, text='Borda 4', font=('Verdana', 10), variable=var_17) # Forças em Y - Borda 4
    botao_18.grid(row=8, column=12)

    # Criando as legendas e entradas para os dados de entrada das CC do Pré-Processamento
    label_37 = Label(J5, text='Magnitude [N]', font=('Verdana', 12)) # Magnitude das Forças
    label_37.grid(row=9, column=12, ipady=10, ipadx=10)

    global entrada_9
    entrada_9 = Entry(J5, width=10, font=('Verdana', 10), justify='center')
    entrada_9.grid(row=9, column=10, pady=10, columnspan=2)

    # Criando o quadro e legendas para o quadro de entrada dos Parâmetros da Análise de EF do Pré-Processamento
    J6 = LabelFrame(J1)
    J6.grid(row=2, column=2, padx=10, pady=10, columnspan=2)

    label_38 = Label(J6, text='Parâmetros da Análise de EF', font=('Verdana', 15))
    label_38.grid(row=0, column=0, ipady=10, columnspan=2)

    # Criando as legendas e entradas para os dados de entrada dos Parâmetros da Análise de EF do Pré-Processamento
    label_39 = Label(J6, text='Nº Incrementos', font=('Verdana', 12)) # Número de Passos
    label_39.grid(row=1, column=0, ipady=10)

    global entrada_10
    entrada_10 = Entry(J6, width=5, font=('Verdana', 10), justify='center')
    entrada_10.grid(row=1, column=1, padx=10)

    # Criando as legendas e checkboxs para os dados de entrada dos Parâmetros da Análise de EF do Pré-Processamento

    label_40 = Label(J6, text='Nº Pontos de Gauss', font=('Verdana', 12)) # NGP
    label_40.grid(row=2, column=0, ipady=10, ipadx=10)

    global var_18
    var_18 = StringVar()
    var_18.set('1')
    botao_19 = Radiobutton(J6, text='1', font=('Verdana', 12), variable=var_18, value='1') # NGP = 1
    botao_19.grid(row=2, column=1)
    botao_20 = Radiobutton(J6, text='2', font=('Verdana', 12), variable=var_18, value='2') # NGP = 2
    botao_20.grid(row=3, column=1)
    botao_21 = Radiobutton(J6, text='3', font=('Verdana', 12), variable=var_18, value='3') # NGP = 3
    botao_21.grid(row=4, column=1)
    botao_22 = Radiobutton(J6, text='4', font=('Verdana', 12), variable=var_18, value='4') # NGP = 4
    botao_22.grid(row=5, column=1)

    label_41 = Label(J6, text='Estado Plano', font=('Verdana', 12)) # Estado Plano
    label_41.grid(row=6, column=0, ipady=10)

    global var_19
    var_19 = StringVar()
    var_19.set('1')
    botao_23 = Radiobutton(J6, text='Tensões', font=('Verdana', 12), variable=var_19, value='1') # Tensões
    botao_23.grid(row=6, column=1)
    botao_24 = Radiobutton(J6, text='Deformações', font=('Verdana', 12), variable=var_19, value='2') # Deformações
    botao_24.grid(row=7, column=1)

    # Criando botões para Execução do Pré-Processamento
    botao_25 = Button(J1, text='Executar', font=('Verdana', 15), command=clicar)
    botao_25.grid(row=3, column=2, pady=10)

    app.mainloop()