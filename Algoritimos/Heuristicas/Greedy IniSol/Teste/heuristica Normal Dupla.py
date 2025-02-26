import turtle 
import random
from geradores import Area, Gerador, posicionar, on_click,ponto_dentro_poligono
from matriz import Matriz
from botao import Botao
import pyautogui
import math
from itertools import product
import numpy as np
import time
import csv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt36
from PIL import Image, ImageDraw
import math
import numba

class PoligonoInscrito():
    def __init__(self, num_lados, tamanho_lado):
        self._num_lados = num_lados
        self._tamanho_lado = tamanho_lado
        if self._num_lados == 3:
            self._altura = int(self._tamanho_lado * math.sqrt(3) / 2)
            self._base = self._tamanho_lado
        elif self._num_lados == 4:
            self._base = self._altura = self._tamanho_lado
        else:
            self._raio = tamanho_lado / (2 * math.sin(math.pi / num_lados))
            self._base = self._altura = int(2 * self._raio)
        self._img = Image.new('L', (self._base, self._altura), 0)
        self._matriz = np.array(self._img) 


    def _criar_poligono(self):
        angulo = 2 * math.pi / self._num_lados
        if self._num_lados == 4:
            ImageDraw.Draw(self._img).rectangle([(0, 0), (self._base, self._altura)], outline=1, fill=1)
        else:
            if self._num_lados == 6:
                angulo_inicial = 0
            else:
                angulo_inicial = math.pi / 2 if self._num_lados % 2 == 0 else math.pi / 2 - math.pi / self._num_lados
            if self._num_lados <= 4:
                raio = self._tamanho_lado / 2
            else:
                raio = self._raio
            poligono = [(math.cos(i * angulo + angulo_inicial) * raio, math.sin(i * angulo + angulo_inicial) * raio) for i in range(self._num_lados)]
            if self._num_lados == 3:
                poligono_transformado = [(x + self._base // 2, y + self._altura - self._tamanho_lado//3) for x, y in poligono]
            else:
                poligono_transformado = [(x + self._base // 2, y + self._altura // 2) for x, y in poligono]
            ImageDraw.Draw(self._img).polygon(poligono_transformado, outline=1, fill=1)
        self._matriz = np.array(self._img)

    def _plot(self):
        self._ax.clear()
        self._ax.imshow(self._matriz)
        plt.pause(0.001)


def transformar_coordenadas_inversa(x, y, x_min, x_max, y_min, y_max, width, height):
    x_scaled = (x / width) * (x_max - x_min) + x_min
    y_scaled = ((height - y) / height) * (y_max - y_min) + y_min

    return x_scaled, y_scaled

def criar_matrizes(lista):
    lista_matrizes = []
    j = 0
    for peca in lista:
        print(Gerador._geracao[j])
        j+=1
        NumeroLados = int(peca.split(',')[0])
        TamanhoLado = int(float(peca.split(',')[2]))
        Altura = int(int(float(peca.split(',')[3]))*2)
        mat = PoligonoInscrito(NumeroLados, TamanhoLado)
        mat._criar_poligono()
        #mat._plot()
        lista_matrizes.append(mat._matriz)
    return lista_matrizes













    
import numpy as np

import math

import numpy as np
#@numba.jit
import numpy as np

import numpy as np

import numpy as np

import numpy as np

import numpy as np
from scipy.ndimage import distance_transform_edt

import numpy as np

import numpy as np

import numpy as np
def raio_circunscrito(n, s):
    r = s / (2 * math.sin(math.pi / n))
    return r

def check_fit(matriz, numero, tamanho, area, apotema, x, y):
    if numero == 3 or numero == 4:
        submatriz = matriz[y-apotema:, int(x-tamanho//2):int(x+tamanho//2)+1]
    if numero == 5:
        distancia = tamanho * ((1+math.sqrt(5))/2)
        submatriz = matriz[y-apotema:, int(x-distancia//2):int(x+distancia//2)+1]
    if numero == 6:
        raio = raio_circunscrito(numero, tamanho)
        submatriz = matriz[y-apotema:, int(x-raio//2):int(x+raio//2)+1]

    area_disponivel = submatriz[np.where(submatriz == 0)].size
    if area_disponivel >= area:
        return True
    else:
        return False



import math

def calcula_pontos_poligono(x, y, n, s):
    if n == 3:  # Triângulo
        apotema = s / (2 * math.tan(math.pi / n))
        return (int(x-(s//2)),int(y-apotema)),(int(x+(s//2)),int(y-apotema)),(x,int(y+(2*apotema)))
    if n == 4:  # Quadrado
        apotema = s//2
        return (x+apotema,y+apotema),(x-apotema,y+apotema),(x+apotema,y-apotema),(x-apotema,y-apotema)

    if n == 5:  # Pentágono
        raio = raio_circunscrito(n,s)
        apotema = s / (2 * math.tan(math.pi / n))
        diagonal = 2 * s * math.sin(math.pi / 5)
        altura = raio - (math.sqrt((s*s)-((diagonal//2) * (diagonal//2))))

        return (x-(s//2), y-apotema), (x+(s//2), y-apotema),(x,y+raio)

    if n == 6:  # Hexágono
        return (x, y+s), (x-s*math.sin(math.radians(30)), y+s*math.cos(math.radians(30))), (x-s, y), (x-s*math.sin(math.radians(30)), y-s*math.cos(math.radians(30))), (x, y-s), (x+s*math.sin(math.radians(30)), y-s*math.cos(math.radians(30))), (x+s, y), (x+s*math.sin(math.radians(30)), y+s*math.cos(math.radians(30)))



def verifica_ponto(x2, y2, lista):
    x1,y1 = transformar_coordenadas_inversa(x2, y2, x, 250, y, -100, base, altura)
    tecido = [(x,y),(x+base,y),(x+base,y-altura),(x,y-altura)]
    if not ponto_dentro_poligono(x1,y1,tecido):
        return True
    else:
        if len(lista) > 0:
            for pol in lista:
                if ponto_dentro_poligono(x1, y1, pol):
                    return True
            return False
        else:
            return False


def verifica_ponto_poligono(x, y, n, s, lista_cordenada):
    poligono = []
    poligono = calcula_pontos_poligono(x,y,n,s)
    for ponto in poligono:
        if verifica_ponto(ponto[0],ponto[1],lista_cordenada):
            return False
    return True








def exact_fit(geraçao, matriz):
    t = turtle.Turtle()
    turtle.tracer(0)
    t.color("black")
    t.penup()
    t.goto(-900,300)
    t.pendown()
    inicio = time.time()
    tempo_anterior = 0
    lista = criar_matrizes(geraçao)
    # Percorra cada peça na ordem predefinida
    while len(Gerador._geracao) > 0:
        i = 0
        i2 = 600
        while i < base:
            j = 0
            j2 = 400
            while j < altura:
                Gerador.Colocou = False
                fim = time.time()
                tempo_execucao = fim - inicio 

                if int(tempo_execucao) != int(tempo_anterior):
                    tempo(tempo_execucao)
                    tempo_anterior = int(tempo_execucao)
                    
                x_transformado, y_transformado = transformar_coordenadas_inversa(i, j, x, 250, y, -100, base, altura)
                x2_transformado, y2_transformado = transformar_coordenadas_inversa(i2, j2, x, 250, y, -100, base, altura)
                #draw_piece(x_transformado,y_transformado)
                if (len(Gerador._geracao)>0):
                    if not verifica_ponto(i,j,Gerador._CordenadasProibidas):
                        if verifica_ponto_poligono(i,j,int(float(Gerador._geracao[-1].split(',')[0])),int(float(Gerador._geracao[-1].split(',')[2])), Gerador._CordenadasProibidas):
                            #if check_fit(matriz,int(float(Gerador._geracao[-1].split(',')[0])),int(float(Gerador._geracao[-1].split(',')[2])),int(float(Gerador._geracao[-1].split(',')[1])),int(float(Gerador._geracao[-1].split(',')[3])),i,j):
                            #draw_piece(x_transformado,y_transformado)
                            on_click2(x_transformado, y_transformado)
                else:
                    break

                if len(Gerador._geracao)>0:
                    if Gerador.Colocou:
                        i, j = int(float(Gerador._geracao[-1].split(',')[3])), int(float(Gerador._geracao[-1].split(',')[3]))+1
                        i2, j2 = 600-int(float(Gerador._geracao[-1].split(',')[3])), 400-int(float(Gerador._geracao[-1].split(',')[3]))+1
                        break

                if (len(Gerador._geracao)>0):
                    if not verifica_ponto(i2,j2,Gerador._CordenadasProibidas):
                        if verifica_ponto_poligono(i2,j2,int(float(Gerador._geracao[-1].split(',')[0])),int(float(Gerador._geracao[-1].split(',')[2])), Gerador._CordenadasProibidas):
                            #if check_fit(matriz,int(float(Gerador._geracao[-1].split(',')[0])),int(float(Gerador._geracao[-1].split(',')[2])),int(float(Gerador._geracao[-1].split(',')[1])),int(float(Gerador._geracao[-1].split(',')[3])),i,j):
                            #draw_piece(x_transformado,y_transformado)
                            on_click2(x2_transformado, y2_transformado)
                else:
                    break

                if len(Gerador._geracao)>0:
                    if Gerador.Colocou:
                        i, j = int(float(Gerador._geracao[-1].split(',')[3])), int(float(Gerador._geracao[-1].split(',')[3]))+1
                        i2, j2 = 600-int(float(Gerador._geracao[-1].split(',')[3])), 400-int(float(Gerador._geracao[-1].split(',')[3]))+1
                        break
                j+=10
                j2-=10
            i += 10
            i2 -= 10
        if len(Gerador._geracao)>0:
            Gerador._geracao.pop()
    return tempo_execucao



    

text_turtle = turtle.Turtle()
fecho_turtle = turtle.Turtle()
def on_click2(x,y):
    on_click(x,y,matriz2,text_turtle, fecho_turtle, base, altura, greedys, Bot)




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


with open('Heuristica Normal Dupla.csv', 'a', newline='') as arquivo:
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
    posicionar(100, "M", True)
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
    tempoF = exact_fit(Gerador._geracao, np.array((matriz2._matriz)))
    with open('Heuristica Normal Dupla.csv', 'a', newline='') as arquivo:
        writer = csv.writer(arquivo)
        writer.writerow([i+1, (round(Area.AreaOcupada/(base*altura)*100, 1)), int(tempoF) ])
    turtle.reset()

with open('Heuristica Normal Dupla.csv', 'a', newline='') as arquivo:
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
    posicionar(60, "G", True)
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
    with open('Heuristica Normal Dupla.csv', 'a', newline='') as arquivo:
        writer = csv.writer(arquivo)
        writer.writerow([i+1, (round(Area.AreaOcupada/(base*altura)*100, 1)), int(tempoF) ])
    turtle.reset()


with open('Heuristica Normal Dupla.csv', 'a', newline='') as arquivo:
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
    tempoF = exact_fit(Gerador._geracao, matriz2)
    with open('Heuristica Normal Dupla.csv', 'a', newline='') as arquivo:
        writer = csv.writer(arquivo)
        writer.writerow([i+1, (round(Area.AreaOcupada/(base*altura)*100, 1)), int(tempoF) ])
    turtle.reset()





