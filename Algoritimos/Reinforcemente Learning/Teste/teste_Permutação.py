from Aprendizado_Permutação import Net
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
import torch
import torch.nn as nn
import torch.optim as optim



def transformar_coordenadas_inversa(x, y, x_min, x_max, y_min, y_max, width, height):
    x_scaled = (x / width) * (x_max - x_min) + x_min
    y_scaled = ((height - y) / height) * (y_max - y_min) + y_min

    return x_scaled, y_scaled


def raio_circunscrito(n, s):
    r = s / (2 * math.sin(math.pi / n))
    return r


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
    #lista = criar_matrizes(geraçao)
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


model = Net(4,400,500,100)
optimizer = optim.SGD(model.parameters(), lr=0.01)

# Treine a rede neural
for epoch in range(100):  # Suponha que estamos fazendo 100 épocas de treinamento
    optimizer.zero_grad()
    # Suponha que Gerador._geracao seja sua lista de strings
    lista_de_strings = Gerador._geracao
    lista_de_listas = [s.split(',') for s in lista_de_strings]
    lista_de_tensores = [torch.tensor([float(num) for num in lista]) for lista in lista_de_listas]  # Zero os gradientes
    outputs = model.forward(lista_de_tensores)  # Faça uma passagem para a frente
    Gerador._geracao = outputs.copy
    exact_fit(Gerador._geracao, np.array(matriz2._matriz))
    loss = Area.Last_Reward[-1]  # Calcule a perda
    loss.backward()  # Faça uma passagem para trás
    optimizer.step()
