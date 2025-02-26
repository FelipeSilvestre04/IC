import turtle 
import random
from geradores import Area, Gerador, posicionar, on_click, ponto_dentro_poligono
from matriz import Matriz
from botao import Botao
import pyautogui
import math
from itertools import product
import numpy as np
import time
import csv
import numba

def transformar_coordenadas_inversa(x, y, x_min, x_max, y_min, y_max, width, height):
    x_scaled = (x / width) * (x_max - x_min) + x_min
    y_scaled = ((height - y) / height) * (y_max - y_min) + y_min

    return x_scaled, y_scaled

def criar_matrizes(lista):
    lista_matrizes = []
    # Percorrendo a lista de trás para frente
    for peca in reversed(lista):
        NumeroLados = int(peca.split(',')[0])
        TamanhoLado = int(float(peca.split(',')[2]))
        Altura = int(int(float(peca.split(',')[3]))*2)
        # Criando a matriz para a peça
        matriz = criar_matriz_poligono(NumeroLados, TamanhoLado, Altura)
        lista_matrizes.append(matriz)
    return lista_matrizes


def criar_matriz_poligono(NumeroLados, TamanhoLado, Altura):
    # Criando uma matriz de zeros
    matriz = np.zeros((Altura, Altura))

    # Preenchendo a matriz para formar o polígono
    for i in range(Altura):
        for j in range(Altura):
            if i == 0 or i == TamanhoLado - 1 or j == 0 or j == NumeroLados - 1:
                matriz[i][j] = 1
    matriz = matriz.astype(np.uint8)
    return matriz


def verifica_ponto(x,y,lista):
    if len(lista) > 0:
        for pol in lista:
            if ponto_dentro_poligono(x,y,pol):
                return False
        return True
    else:
        return True

def exact_fit(geraçao, matriz,t):
    t = turtle.Turtle()
    turtle.tracer(0)
    t.color("black")
    t.penup()
    t.goto(-900,300)
    t.pendown()
    inicio = time.time()
    tempo_anterior = 0
    # Percorra cada peça na ordem predefinida

    i = 0
    while i < base:
        j = 0
        while j < altura:
            fim = time.time()
            tempo_execucao = fim - inicio 
            # Atualize o tempo a cada segundo
            if int(tempo_execucao) != int(tempo_anterior):
                tempo(tempo_execucao)
                tempo_anterior = int(tempo_execucao)
            x_transformado, y_transformado = transformar_coordenadas_inversa(i, j, x, 250, y, -100, base, altura)    
            if verifica_ponto(x_transformado,y_transformado, Gerador._CordenadasProibidas):     
                on_click2(x_transformado, y_transformado)
            #draw_piece(x_transformado, y_transformado)
            if Gerador.Colocou:
                i, j = 0, 0
                break
            j+=10
        i += 10
    return tempo_execucao



              



def check_fit(tecido, peca, x, y):
    altura_peca, largura_peca = peca.shape
    altura_tecido, largura_tecido = tecido.shape

    # Verificando se a peça cabe no tecido
    if x + altura_peca > altura_tecido or y + largura_peca > largura_tecido:
        return False

    # Verificando se a peça sobrepõe algum polígono existente no tecido
    for i in range(altura_peca):
        for j in range(largura_peca):
            if tecido[x + i][y + j] == 1 and peca[i][j] == 1:
                return False

    return True

text_turtle = turtle.Turtle()
fecho_turtle = turtle.Turtle()
def on_click2(x,y):
    on_click(x,y,matriz2,text_turtle, fecho_turtle, base, altura, greedys, Bot)


def raio_circunscrito(n, s):
    r = s / (2 * math.sin(math.pi / n))
    return r

def draw_piece(x,y):
    # Pega a posição atual do mouse
    
    if len(Gerador._geracao) > 0:
        shape_name = Gerador._geracao[-1]
        n = int(shape_name.split(',')[0])  
        s = int(shape_name.split(',')[2])  
        turtle.tracer(0)
        turtle.clear()
        turtle.penup()
        turtle.goto(x,y)
        turtle.pendown()
        turtle.right(90)  #
        turtle.right(180/n)  
        turtle.forward(raio_circunscrito(n, s))  
        turtle.left(180/n)  
        turtle.left(90)
        turtle.color('green')
        turtle.begin_fill()
        for _ in range(int(shape_name.split(',')[0])):
                turtle.forward(int(shape_name.split(',')[2]))
                turtle.left(360 / int(shape_name.split(',')[0]))
        turtle.end_fill()
        turtle.hideturtle()
        turtle.tracer(1)

def draw_piece2(x,y):
    # Pega a posição atual do mouse
    t = turtle.Turtle()
    t.clear()
    if len(Gerador._geracao) > 0:
        shape_name = Gerador._geracao[-1]
        n = int(shape_name.split(',')[0])  
        s = int(shape_name.split(',')[2])
        turtle.tracer(0)
        t.clear()
        t.penup()
        t.goto(x,y)
        t.pendown()
        t.right(90)  #
        t.right(180/n)  
        t.forward(raio_circunscrito(n, s))  
        t.left(180/n)  
        t.left(90)
        t.color('green')
        t.begin_fill()
        for _ in range(int(shape_name.split(',')[0])):
                t.forward(int(shape_name.split(',')[2]))
                t.left(360 / int(shape_name.split(',')[0]))
        t.end_fill()
        t.hideturtle()
        turtle.tracer(1)
#print(Gerador._geracao)

# Loop principal


def tempo(time):
    # Limpe a escrita anterior antes de escrever o novo tempo
    t.clear()
    t.write(f"{int(time/60)} min {round(time%60,0)} s",font=("Arial", 16, "normal"))
    t.hideturtle()
    turtle.update()
t = turtle.Turtle()
with open('Heuristica Normal Simples.csv', 'a', newline='') as arquivo:
        writer = csv.writer(arquivo)
        writer.writerow(['Heuristicas - Pequenas:'])

for i in range(10):
    Area.Last_Reward = []
    Area.AreaTotal = []
    Area.AreaGridN = [0.00,0.00,0.00,0.00,0.00,0.00]
    Area.AreaOcupada = 0
    Area.CordenadasProibidas = []
    Area.CordenadasGreedy = []
    Area.AreasGredys = []
    Area.GreedysOcupados = []
    Area.Q_total = -6

    Gerador.Colocou = False
    Gerador.N = False
    Gerador._count = 0
    Gerador._geracaopop = []
    Gerador._geracao = []
    Gerador._turtles = []
    Gerador._CordenadasProibidas = []
    Gerador._desenhando = False
    Gerador._poligonos_posicionados = []
    Gerador._Desenahdos = []
    Gerador._PolEGrid = []

    turtle.clearscreen()
    base = 600
    altura = 400
    x = -350
    y = 300
    greedys = 6
    turtle.Turtle()
    turtle.tracer(0)
    turtle.speed(200)
    tecido = Area(base, altura, greedys)
    matriz2 = Matriz(base, altura, x, y)
    tecido.area(x, y)
    tecido.greedy(x, y)
    isdrawing = False
    posicionar(130, "P")
    Bot = Botao(350, 300, 200, 150)
    turtle.tracer(1)


    t = turtle.Turtle()
    turtle.tracer(0)
    t.color("black")
    t.penup()
    t.goto(-925,250)
    t.pendown()
    t.hideturtle()

    loop = True
    inicio = time.time()
    tempo_anterior = 0

    fim = time.time()
    tempo_execucao = fim - inicio 
        # Atualize o tempo a cada segundo
    if int(tempo_execucao) != int(tempo_anterior):
        tempo(tempo_execucao)
        tempo_anterior = int(tempo_execucao)
    tempoF = exact_fit(Gerador._geracao, matriz2, t)
    with open('Heuristica Normal Simples.csv', 'a', newline='') as arquivo:
        writer = csv.writer(arquivo)
        writer.writerow([i+1, (round(Area.AreaOcupada/(base*altura)*100, 1)), int(tempoF) ])
    turtle.reset()

with open('Heuristica Normal Simples.csv', 'a', newline='') as arquivo:
        writer = csv.writer(arquivo)
        writer.writerow(['Heuristicas - Mistas:'])
for i in range(10):
    Area.Last_Reward = []
    Area.AreaTotal = []
    Area.AreaGridN = [0.00,0.00,0.00,0.00,0.00,0.00]
    Area.AreaOcupada = 0
    Area.CordenadasProibidas = []
    Area.CordenadasGreedy = []
    Area.AreasGredys = []
    Area.GreedysOcupados = []
    Area.Q_total = -6

    Gerador.Colocou = False
    Gerador.N = False
    Gerador._count = 0
    Gerador._geracaopop = []
    Gerador._geracao = []
    Gerador._turtles = []
    Gerador._CordenadasProibidas = []
    Gerador._desenhando = False
    Gerador._poligonos_posicionados = []
    Gerador._Desenahdos = []
    Gerador._PolEGrid = []

    turtle.clearscreen()
    base = 600
    altura = 400
    x = -350
    y = 300
    greedys = 6
    turtle.Turtle()
    turtle.tracer(0)
    turtle.speed(200)
    tecido = Area(base, altura, greedys)
    matriz2 = Matriz(base, altura, x, y)
    tecido.area(x, y)
    tecido.greedy(x, y)
    isdrawing = False
    posicionar(80, "M")
    Bot = Botao(350, 300, 200, 150)
    turtle.tracer(1)


    t = turtle.Turtle()
    turtle.tracer(0)
    t.color("black")
    t.penup()
    t.goto(-925,250)
    t.pendown()
    t.hideturtle()

    loop = True
    inicio = time.time()
    tempo_anterior = 0

    fim = time.time()
    tempo_execucao = fim - inicio 
        # Atualize o tempo a cada segundo
    if int(tempo_execucao) != int(tempo_anterior):
        tempo(tempo_execucao)
        tempo_anterior = int(tempo_execucao)
    tempoF = exact_fit(Gerador._geracao, matriz2)
    with open('Heuristica Normal Simples.csv', 'a', newline='') as arquivo:
        writer = csv.writer(arquivo)
        writer.writerow([i+1, (round(Area.AreaOcupada/(base*altura)*100, 1)), int(tempoF) ])
    turtle.reset()

with open('Heuristica Normal Simples.csv', 'a', newline='') as arquivo:
        writer = csv.writer(arquivo)
        writer.writerow(['Heuristicas - Grandes:'])
for i in range(10):
    Area.Last_Reward = []
    Area.AreaTotal = []
    Area.AreaGridN = [0.00,0.00,0.00,0.00,0.00,0.00]
    Area.AreaOcupada = 0
    Area.CordenadasProibidas = []
    Area.CordenadasGreedy = []
    Area.AreasGredys = []
    Area.GreedysOcupados = []
    Area.Q_total = -6

    Gerador.Colocou = False
    Gerador.N = False
    Gerador._count = 0
    Gerador._geracaopop = []
    Gerador._geracao = []
    Gerador._turtles = []
    Gerador._CordenadasProibidas = []
    Gerador._desenhando = False
    Gerador._poligonos_posicionados = []
    Gerador._Desenahdos = []
    Gerador._PolEGrid = []

    turtle.clearscreen()
    base = 600
    altura = 400
    x = -350
    y = 300
    greedys = 6
    turtle.Turtle()
    turtle.tracer(0)
    turtle.speed(200)
    tecido = Area(base, altura, greedys)
    matriz2 = Matriz(base, altura, x, y)
    tecido.area(x, y)
    tecido.greedy(x, y)
    isdrawing = False
    posicionar(60, "G")
    Bot = Botao(350, 300, 200, 150)
    turtle.tracer(1)


    t = turtle.Turtle()
    turtle.tracer(0)
    t.color("black")
    t.penup()
    t.goto(-925,250)
    t.pendown()
    t.hideturtle()

    loop = True
    inicio = time.time()
    tempo_anterior = 0

    fim = time.time()
    tempo_execucao = fim - inicio 
        # Atualize o tempo a cada segundo
    if int(tempo_execucao) != int(tempo_anterior):
        tempo(tempo_execucao)
        tempo_anterior = int(tempo_execucao)
    tempoF = exact_fit(Gerador._geracao, matriz2)
    with open('Heuristica Normal Simples.csv', 'a', newline='') as arquivo:
        writer = csv.writer(arquivo)
        writer.writerow([i+1, (round(Area.AreaOcupada/(base*altura)*100, 1)), int(tempoF) ])
    turtle.reset()








