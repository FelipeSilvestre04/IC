import turtle 
from geradores import Area, Gerador, posicionar, on_click
from matriz import Matriz
from botao import Botao
import math
import time
import gymnasium as gym
from gymnasium import spaces
import numpy as np
from stable_baselines3 import PPO
from stable_baselines3.common.env_checker import check_env
import matplotlib.pyplot as plt
from stable_baselines3.common.callbacks import BaseCallback
import turtle
from geradores import Area
from botao import Botao
import math
import numpy as np
from scipy.spatial import ConvexHull
import time
import random
import csv
import copy
import matplotlib.pyplot as plt
from matriz_vertices import Matriz
from shapely.geometry import Polygon
import os
from stable_baselines3.common.callbacks import EvalCallback
from stable_baselines3 import PPO
from stable_baselines3.common.vec_env import SubprocVecEnv
import tensorflow as tf
import torch
import torch.nn as nn
from stable_baselines3.common.torch_layers import BaseFeaturesExtractor


def ponto_dentro_poligono(x, y, poligono, tag):
    def ponto_na_aresta(px, py, p1, p2):
        # Verifica se o ponto (px, py) está na aresta definida por p1 e p2
        if p1[1] == p2[1]:  # Verifica horizontal
            if p1[1] == py and min(p1[0], p2[0]) <= px <= max(p1[0], p2[0]):
                return True
        if p1[0] == p2[0]:  # Verifica vertical
            if p1[0] == px and min(p1[1], p2[1]) <= py <= max(p1[1], p2[1]):
                return True
        # Verifica diagonal
        if min(p1[0], p2[0]) <= px <= max(p1[0], p2[0]) and min(p1[1], p2[1]) <= py <= max(p1[1], p2[1]):
            if (p2[0] - p1[0]) * (py - p1[1]) == (px - p1[0]) * (p2[1] - p1[1]):
                return True
        return False

    n = len(poligono)
    dentro = False

    p1x, p1y = poligono[0]
    for i in range(n + 1):
        p2x, p2y = poligono[i % n]
        
        # Verifica se o ponto está na borda do polígono
        if ponto_na_aresta(x, y, (p1x, p1y), (p2x, p2y)):
            return tag
        
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xints:
                        dentro = not dentro
        p1x, p1y = p2x, p2y

    return dentro
def CordenadaProibida(x, y):
    global CordenadasProibidas
    if (not ponto_dentro_poligono(x, y, Area.CordenadasProibidas, True)):
        return True
    else:
        for pol in CordenadasProibidas:
            if ponto_dentro_poligono(x, y, pol,False):
                return True
 
        return False 
def mdc(a, b):
    while b:
        a, b = b, a % b
    return a
def pontos_entre_vertices(vertices):
    def mdc(a, b):
        # Função auxiliar para calcular o MDC (Máximo Divisor Comum)
        while b:
            a, b = b, a % b
        return a

    pontos = []
  
    # Adiciona o último vértice ao início para conectar o último ao primeiro
    vertices = vertices + [vertices[0]]

    for i in range(len(vertices) - 1):
        x1, y1 = vertices[i]
        x2, y2 = vertices[i + 1]
        # Calcula as diferenças
        dx = x2 - x1
        dy = y2 - y1
        # Calcula o MDC das diferenças
        passo_mdc = mdc(abs(dx), abs(dy))
        # Calcula os incrementos
        inc_x = dx // passo_mdc if passo_mdc != 0 else 0
        inc_y = dy // passo_mdc if passo_mdc != 0 else 0
        # Adiciona os pontos à lista
        for j in range(0, passo_mdc + 1):  # Inclui o ponto final
            ponto = (x1 + j * inc_x, y1 + j * inc_y)
            pontos.append(ponto)

    return pontos

def rotate_point(x, y, angle):
    rad = math.radians(angle)
    x_new = x * math.cos(rad) - y * math.sin(rad)
    y_new = x * math.sin(rad) + y * math.cos(rad)
    return x_new, y_new
def draw_piece(i,j, grau,l):
    global lista_poligonos
    if l < len(lista_poligonos):
        x, y = i,j
        pol = lista_poligonos[l]
        turtle.tracer(0)
        turtle.clear()
        turtle.begin_fill()
        turtle.fillcolor("green")
        turtle.penup()
        x_rot, y_rot = rotate_point(pol[0][0], pol[0][1], grau)
        turtle.goto(x+x_rot*Escala, y+y_rot*Escala)
        turtle.pendown()
        for cor in pol:

            x_rot, y_rot = rotate_point(cor[0], cor[1], grau)
            turtle.goto(x+x_rot*Escala, y+y_rot*Escala)
        turtle.goto(x+x_rot*Escala, y+y_rot*Escala)
        turtle.end_fill()
        turtle.hideturtle()
        turtle.tracer(1)
def calcular_area(vertices):
    num_vertices = len(vertices)
    area = 0

    for i in range(num_vertices):
        x1, y1 = vertices[i]
        x2, y2 = vertices[(i + 1) % num_vertices]  # Próximo vértice (considerando o último vértice conectado ao primeiro)

        area += (x1 * y2) - (x2 * y1)

    return abs(area) / 2
def variacao_total_eixo_x(poligonos):
    # Inicializa min_x e max_x com valores extremos para garantir que serão substituídos
    min_x = float('inf')
    max_x = float('-inf')
    
    for poligono in poligonos:
        for ponto in poligono:
            x, _ = ponto
            # Atualiza min_x e max_x se encontrar um valor menor ou maior
            min_x = min(min_x, x)
            max_x = max(max_x, x)
    
    # Calcula a variação total do eixo X para todos os polígonos
    variacao_total = max_x - min_x
    return variacao_total
def desenhar_fecho(vertices_fecho, p):
    turtle.tracer(0)
    p.clear()
   
    p.speed(5)

    
    # Levantar a caneta antes de mover para a primeira coordenada
    p.penup()
    p.goto(vertices_fecho[0])
    
    # Abaixar a caneta e começar a desenhar
    p.pendown()
    for cordenada in vertices_fecho[1:]:
        p.goto(cordenada)
    
    # Fechar o polígono voltando ao primeiro vértice
    p.goto(vertices_fecho[0])
    
    # Esconder a tartaruga e terminar o desenho
    p.hideturtle()
    turtle.tracer(1)
def calcular_fecho_e_area(poligonos):
    # Concatenar todos os vértices de todos os polígonos
    vertices = np.concatenate(poligonos, axis=0)
    
    # Calcular o fecho convexo
    fecho = ConvexHull(vertices)
    
    # Obter os vértices do fecho convexo
    vertices_fecho = vertices[fecho.vertices]
    
    # Calcular a área do fecho convexo
    area = fecho.volume
    
    return vertices_fecho, area
def area_desperdiçada_fecho(area_fecho):
    area_ocupada = 0
    Area.AreaOcupada = 0
    #for poligono in Gerador._poligonos_posicionados:
        #area_ocupada += float(poligono.split(',')[1])
    area_fecho = area_fecho - area_ocupada
    Area.AreaOcupada = area_ocupada
    return area_fecho
def ler_poligonos(arquivo):
    with open(arquivo, 'r') as f:
        conteudo = f.read().strip()

    # Divide o conteúdo em linhas
    linhas = conteudo.split('\n')

    # Lê o número total de polígonos
    num_poligonos = int(linhas[0].strip())
    #print(f"Número total de polígonos: {num_poligonos}")

    poligonos = []
    i = 1  # Começa a leitura a partir da segunda linha

    while i < len(linhas):
        # Verifica se a linha não está vazia
        if linhas[i].strip():
            try:
                # Lê o número de vértices
                num_vertices = int(linhas[i].strip())
                #print(f"Lendo polígono com {num_vertices} vértices")  # Depuração
                i += 1

                # Lê os vértices
                vertices = []
                for _ in range(num_vertices):
                    # Verifica se a linha não está vazia
                    while i < len(linhas) and not linhas[i].strip():
                        i += 1
                    if i < len(linhas):
                        coords = linhas[i].strip().split()
                        if len(coords) != 2:
                            raise ValueError(f"Esperado 2 valores por linha, mas obteve {len(coords)}: '{linhas[i].strip()}'")
                        x, y = map(float, coords)
                        vertices.append((x, y))
                        i += 1
                    else:
                        raise ValueError(f"Esperado {num_vertices} vértices, mas o arquivo terminou prematuramente.")

                poligonos.append(vertices)
            except ValueError as ve:
                print(f"Erro ao processar a linha {i}: {linhas[i].strip()} - {ve}")
                i += 1
        else:
            i += 1
    if num_poligonos == len(poligonos):
        pass
        #print(f'Todos os {num_poligonos} poligonos foram lidos com sucesso!')
    return poligonos
def tratar_lista(lista_poligonos):
    nova_lista = []
    nova_lista_completa = []
    for pol in lista_poligonos:
        novo_pol = []
        for cor in pol:
            novo_pol.append((int(cor[0]*Escala),int(cor[1]*Escala)))
        pol_completo = pontos_entre_vertices(novo_pol)
        nova_lista_completa.append(pol_completo)
        nova_lista.append(novo_pol)

    return nova_lista, nova_lista_completa
def transformar_coordenadas_inversa(x, y, x_min, x_max, y_min, y_max, width, height):
    x_scaled = (x / width) * (x_max - x_min) + x_min
    y_scaled = ((height - y) / height) * (y_max - y_min) + y_min

    return x_scaled, y_scaled
def calcular_largura(coordenadas):
    x_valores = [coord[0] for coord in coordenadas]
    largura = max(x_valores) - min(x_valores)
    return int(largura)
def calcular_larguraY(coordenadas):
    x_valores = [coord[1] for coord in coordenadas]
    largura = max(x_valores) - min(x_valores)
    return int(largura)
def calcular_centroide(pontos):
    # Número de pontos no polígono
    n = len(pontos)
    
    # Inicialização das somas
    soma_area = 0
    soma_cx = 0
    soma_cy = 0
    
    for i in range(n):
        x0, y0 = pontos[i]
        x1, y1 = pontos[(i + 1) % n]
        
        # Cálculo do termo comum (x_i * y_{i+1} - x_{i+1} * y_i)
        termo = x0 * y1 - x1 * y0
        
        # Somas para área e centroides
        soma_area += termo
        soma_cx += (x0 + x1) * termo
        soma_cy += (y0 + y1) * termo
    
    # Área do polígono
    area = soma_area / 2
    
    # Coordenadas do centroide
    cx = soma_cx / (6 * area)
    cy = soma_cy / (6 * area)
    
    return cx, cy
def tempo(time):
    # Limpe a escrita anterior antes de escrever o novo tempo
    t.clear()
    t.write(f"{int(time/60)} min {round(time%60,0)} s",font=("Arial", 16, "normal"))
    t.hideturtle()
    turtle.update()
def segmentos_se_cruzam(p1, p2, q1, q2):
    def orientacao(a, b, c):
        # Calcula a orientação do triplet (a, b, c)
        val = (b[1] - a[1]) * (c[0] - b[0]) - (b[0] - a[0]) * (c[1] - b[1])
        if val == 0:
            return 0  # colinear
        elif val > 0:
            return 1  # horário
        else:
            return 2  # anti-horário

    def na_segmento(a, b, c):
        # Verifica se o ponto c está no segmento de linha ab
        if min(a[0], b[0]) <= c[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= c[1] <= max(a[1], b[1]):
            return True
        return False

    # Encontra as 4 orientações necessárias para a verificação geral e especial
    o1 = orientacao(p1, p2, q1)
    o2 = orientacao(p1, p2, q2)
    o3 = orientacao(q1, q2, p1)
    o4 = orientacao(q1, q2, p2)

    # Caso geral
    if o1 != o2 and o3 != o4:
        return True

    # Casos especiais
    if o1 == 0 and na_segmento(p1, p2, q1):
        return True
    if o2 == 0 and na_segmento(p1, p2, q2):
        return True
    if o3 == 0 and na_segmento(q1, q2, p1):
        return True
    if o4 == 0 and na_segmento(q1, q2, p2):
        return True

    return False

def poligonos_se_sobrepoem(poli1, poli2):
    n1 = len(poli1)
    n2 = len(poli2)

    # Verifica se algum vértice de poli1 está dentro de poli2
    for p in poli1:
        if ponto_dentro_poligono(p[0], p[1], poli2, False):
            return True

    # Verifica se algum vértice de poli2 está dentro de poli1
    for p in poli2:
        if ponto_dentro_poligono(p[0], p[1], poli1, False):
            return True

    # Verifica se alguma aresta de poli1 se cruza com alguma aresta de poli2
    for i in range(n1):
        for j in range(n2):
            if segmentos_se_cruzam(poli1[i], poli1[(i + 1) % n1], poli2[j], poli2[(j + 1) % n2]):
                return True

    return False
def process_position(i, j, l, grau,base, altura):
        global lista_colunas, lista_area_fechos, fecho, AreaPreenchida
        cordenada_disponivel = True
        poli = nova_lista[l]
        x_transformado, y_transformado = transformar_coordenadas_inversa(i, j, X, X+base, Y, Y-altura, base, altura)

        
        pontos = [(x_transformado + int(rotate_point(cor[0], cor[1], grau)[0]), 
                   y_transformado + int(rotate_point(cor[0], cor[1], grau)[1])) 
                  for cor in poli]
        
        for cor in pontos:
            if cor[0] > Area.CordenadasProibidas[1][0]:
                return False
            if cor[1] > Area.CordenadasProibidas[0][1]:
                return False

        
        if any(not ponto_dentro_poligono(cor[0],cor[1],Area.CordenadasProibidas,True)for cor in pontos):
            return False
        for pol in CordenadasProibidas:
            if any(ponto_dentro_poligono(cor[0],cor[1],pontos, False) for cor in pol):
                return False

        
        if any(CordenadaProibida(cor[0], cor[1]) for cor in pontos):
            return False

        x3, y3 = calcular_centroide(pontos)
        if CordenadaProibida(x3, y3):
            return False


        
        pol = nova_lista_completa[l]
        pontos_pol = [(x_transformado + int(rotate_point(cor[0], cor[1], grau)[0]), 
                       y_transformado + int(rotate_point(cor[0], cor[1], grau)[1])) 
                      for cor in pol]
        
        if any(CordenadaProibida(cor[0], cor[1]) for cor in pontos_pol):
            return False
        
        return True
    
def normalized_distance(point1, point2, rect_width, rect_height):

    distance = np.linalg.norm(np.array(point1) - np.array(point2))
    max_distance = np.linalg.norm(np.array([rect_width, rect_height]))
    normalized_dist = distance / max_distance
    
    return normalized_dist

def on_click2(x, y, grau,fecho_turtle, dados_turtle, botao, n,colunas,base,altura,adicional2,verificado):
    global CordenadasProibidas
    global lista_colunas
    global fecho
    global AreaPreenchida
    global turtles
    global l
    global lista_area_fechos
    global recompensa
    if len(lista_poligonos) > n:
            posição_disponivel = True
            pol = nova_lista_completa[n]
            
            pontos_pol = []
            for cor in pol:
                x_rot, y_rot = rotate_point(cor[0], cor[1], grau)
                pontos_pol.append((int(x) + int(x_rot),int(y) + int(y_rot)))
                if CordenadaProibida(int(x) + int(x_rot),int(y) + int(y_rot)):
                    posição_disponivel = False
                    
            
            for pol in CordenadasProibidas:
                if any(ponto_dentro_poligono(cor[0],cor[1],pontos_pol, False) for cor in pol):
                    posição_disponivel = False
            if posição_disponivel:
                print("passou")

                fecho_turtle.clear()
                dados_turtle.clear()
                poli = nova_lista[n]
                p = turtle.Turtle()
                p.begin_fill()
                p.fillcolor("orange")
                pol_rotação = []
                turtle.tracer(0)
                # Rotaciona o primeiro ponto e move a tartaruga para lá
                x_rot, y_rot = rotate_point(int(poli[0][0]), int(poli[0][1]), grau)
                p.penup()
                p.goto(x+x_rot, y+y_rot)
                p.pendown()
                # Rotaciona e desenha os restantes pontos do polígono
                for cor in poli:
                    x_rot, y_rot = rotate_point(int(cor[0]), int(cor[1]), grau)
                    p.goto(x+x_rot, y+y_rot)
                    pol_rotação.append((x+x_rot, y+y_rot))

                # Fecha o polígono após a rotação
                x_rot, y_rot = rotate_point(int(poli[0][0]), int(poli[0][1]), grau)
                p.goto(x+x_rot, y+y_rot)
                p.end_fill()
                p.hideturtle()
                turtle.tracer(1)
                # Remove o polígono da lista após desenhá-lo
                Matriz_np._criar_poligono(pontos_pol)
                #Matriz_np._plot()
                CordenadasProibidas.append(pontos_pol)
                fecho.append(pontos_pol)
                area_poligono = calcular_area(pontos_pol)
                
                AreaPreenchida += area_poligono
                vertices_fecho = []
                vertices_fecho, area_fecho = calcular_fecho_e_area(fecho)
                #area_desperdiçada = area_desperdiçada_fecho(area_fecho)
                    #print(f"Area desperdicada do fecho convexo: {round(area_desperdiçada, 2)}")

                
                turtle.tracer(0)
                #SumQ_local = calcular_Qs_Local(Area.GreedysOcupados, Area.CordenadasGreedy, Area.AreaOcupada, matriz._matriz, base, altura, greedys, p)
                #punicao = calcular_punição_fecho(area_desperdiçada, len(Area.GreedysOcupados), p, 40000)


                area_total = base*altura  
                dados_turtle.penup()
                dados_turtle.goto(-925,200)
                dados_turtle.pendown()
                função = (((((AreaPreenchida/(variacao_total_eixo_x(fecho)*altura))*0.5)**(0))*((AreaPreenchida/area_fecho)**(1))))

                if len(fecho)>0:
                    função = (((((AreaPreenchida/(variacao_total_eixo_x(fecho)*altura))*0.5)**(0))*((AreaPreenchida/area_fecho)**(1))))        
                    dados_turtle.write(f"Preenchimento Fecho = {round(função, 2)*100}%", font=("Arial",16, "normal"))
                else:
                    dados_turtle.write(f"Preenchimento Fecho = {0}%", font=("Arial",16, "normal"))
                dados_turtle.penup()
                dados_turtle.goto(-925,350)
                dados_turtle.pendown()
                
            
                dados_turtle.write(f"Preenchimento total: {(round(((max(AreaPreenchida,(sum(lista_area_fechos)+AreaPreenchida)))/area_total)*100,2))}%", font=("Arial",16, "normal"))
                dados_turtle.penup()
                dados_turtle.goto(-925,300)
                dados_turtle.pendown()
                dados_turtle.write(f"Peças Disponiveis: {len(lista_poligonos)-1}", font=("Arial",16, "normal"))
                dados_turtle.penup()
                dados_turtle.goto(-925,250)
                dados_turtle.pendown()
                dados_turtle.write(f"Ultima recompensa: {recompensa}", font=("Arial",16, "normal"))
                dados_turtle.hideturtle()
                
                desenhar_fecho(vertices_fecho, fecho_turtle)
                turtle.tracer(1)
                turtles.append([p,lista_poligonos[n],nova_lista[n],nova_lista_completa[n]])
                lista_poligonos.pop(n)
                nova_lista.pop(n)
                nova_lista_completa.pop(n)
                l = 0

                return area_fecho
def transformar_coordenadas_inversa(x, y, x_min, x_max, y_min, y_max, width, height):
    x_scaled = (x / width) * (x_max - x_min) + x_min
    y_scaled = ((height - y) / height) * (y_max - y_min) + y_min

    return x_scaled, y_scaled



import numpy as np
from shapely.geometry import Polygon

from shapely.geometry import Polygon
import numpy as np

def calculate_max_values(features):
    """Calcula os valores máximos para normalização com base nas características fornecidas."""
    max_values = [
        max(f[0] for f in features),
        max(f[1] for f in features),
        max(f[2] for f in features),
        max(f[3] for f in features),
        max(f[4] for f in features),
        max(f[5] for f in features),
        max(f[6] for f in features)
    ]
    return max_values


def polygon_features(lista_poligonos):
    features = []
    
    for coords in lista_poligonos:
        # Agrupar as coordenadas em pares (x, y)
        
        # Criar o polígono usando Shapely
        poly = Polygon(coords)
        
        # Calcular o número de lados
        num_lados = len(coords)
        
        # Calcular a área
        area = poly.area
        
        # Calcular o perímetro
        perimetro = poly.length
        
        # Calcular o centroide
        centroide = poly.centroid
        centroide_x, centroide_y = centroide.x, centroide.y
        
        # Calcular a bounding box
        minx, miny, maxx, maxy = poly.bounds
        width = maxx - minx
        height = maxy - miny
        
        # Adicionar as características à lista
        features.append([
            num_lados, 
            area, 
            perimetro, 
            centroide_x, 
            centroide_y, 
            width, 
            height
        ])
    
    return features




class env2KP(gym.Env):
    def __init__(self, Base=None,Altura=None,x=-200,y=200,peças='fu',draw=False,escala=None,plot=False,passo=1,peso1_coluna = 2, peso2_coluna = 1, peso3_coluna = 2, peso1_pecas = 0, peso2_pecas = 1, verificacao_avancada = False):
        super(env2KP, self).__init__()
        self.fecho = 0
        self.pecas = peças
        self.erro = 0
        self.epoca = 0
        self.acao = 0
        self.tempo_anterior = 0
        self.recompensa = -6
        self.preenchimentos = []
        self.recompensas = []
        self.MediaPreenchimentos = []
        self.MediaRecompensas = []
        self.dist = 0
        
        global X,Y,tempoF,grau,CordenadasProibidas,lista_colunas,fecho,lista_area_fechos,AreaPreenchida,l,fecho_turtle,dados_turtle,turtles,lista_poligonos,t,inicio,tempo_anterior,fim,Bot,nova_lista,nova_lista_completa,base,altura
        global Escala,Matriz_np,Plot,Passo,dado,dado1,bases
        tempoF = 0
        inicio = time.time()
        tempo_anterior = 0

        fim = time.time()
        self.Plot = plot
        self.X = x
        self.Y = y
        X = x
        Y = y
        self.primeiro_poligono = []
        
        self.instancias(Base,Altura,escala)
        self.area = self.base * self.altura

        instancias = '/home/fsilvestre/Cutting_Stock_Problem/'+self.pecas + '.dat'
        lista_poligonos = ler_poligonos(instancias)
        #random.shuffle(lista_poligonos)
        self.numero_pecas = len(lista_poligonos)

            # Espaço de observação: um dicionário que inclui a matriz da área de corte e a lista de peças restantes
        self.observation_space = spaces.Dict({
            'cutting_area': spaces.Box(low=0, high=1, shape=(self.altura, self.base), dtype=np.int8),
            'pieces': spaces.Box(low=-1, high=255, shape=(len(lista_poligonos), 7), dtype=np.int8)
        })

            
            # Espaço de ação: escolha de uma peça e sua posição na área de corte
        self.action_space = spaces.MultiDiscrete([self.base,self.altura,len(lista_poligonos), 4])
        self.desconto = 0.2
        self.done = False
        
        features = []
        self.pecas_colocadas = []
        self.Mpecas_colocadas = []
        for coords in lista_poligonos:
            # Criar o polígono usando Shapely
            poly = Polygon(coords)
            
            # Calcular o número de lados
            num_lados = len(coords)
            
            # Calcular a área
            area = poly.area
            
            # Calcular o perímetro
            perimetro = poly.length
            
            # Calcular o centroide
            centroide = poly.centroid
            centroide_x, centroide_y = centroide.x, centroide.y
            
            # Calcular a bounding box
            minx, miny, maxx, maxy = poly.bounds
            width = maxx - minx
            height = maxy - miny
            
            # Adicionar as características à lista
            features.append([
                num_lados, 
                area, 
                perimetro, 
                centroide_x, 
                centroide_y, 
                width, 
                height
            ])
        
        # Calcular os valores máximos dinamicamente
        self.max_values = calculate_max_values(features)
     
    def reset(self, seed=None):
        global X,Y,tempoF,grau,CordenadasProibidas,lista_colunas,fecho,lista_area_fechos,AreaPreenchida,l,fecho_turtle,dados_turtle,turtles,lista_poligonos,t,inicio,tempo_anterior,fim,Bot,nova_lista,nova_lista_completa,base,altura
        global Escala,Matriz_np,Plot,Passo,dado,dado1,bases, PC1, PC2, PC3, PP1, PP2, va,acao 

        self.recompensa = -6
        self.epoca+=1
        self.acao = 0
        self.done = False
        if seed is not None:
            np.random.seed(seed)
        Area.CordenadasProibidas = [[self.X,self.Y],[self.X+self.base,self.Y],[self.X+self.base,self.Y-self.altura],[self.X,self.Y-self.altura]]

        grau = 0
        CordenadasProibidas = []
        lista_colunas = []
        fecho = []
        lista_area_fechos = []
        deperdicio = 0
        area_fecho = 0
        AreaPreenchida = 0
        l = 0
        Matriz_np = Matriz(self.base,self.altura,self.X,self.Y)

        instancias = '/home/fsilvestre/Cutting_Stock_Problem/fu.dat'
        lista_poligonos = ler_poligonos(instancias)
        #random.shuffle(lista_poligonos)
        nova_lista, nova_lista_completa = tratar_lista(lista_poligonos)
        
        cutting_area = np.array((Matriz_np._matriz), dtype=np.int8)

        # Atualiza o estado de observação como um dicionário
        pieces = (polygon_features(lista_poligonos))
        diferença = self.numero_pecas - len(pieces)
        for _ in range(diferença):
            pieces.append([-1,-1,-1,-1,-1,-1,-1])
        
        pieces = np.array((pieces), dtype=np.int8)
        self.fr = False
        
            
        self.state = {
            'cutting_area': cutting_area,
            'pieces': pieces,
        }
            
        return self.state, {}

    
    
    def step(self, action):
        self.acao += 1
        global acao, recompensa, fig, ax, Matriz_np
        recompensa = self.recompensa

        n = action[2]
        grau_indice = action[3]
        x = action[0]
        y = action[1]

            
        if n < len(lista_poligonos):
            posição_disponivel, pontos_pol, primeira_peça = self.verificar_posicao_disponivel(n, x, y, grau_indice)
            if primeira_peça is not None:
                self.primeiro_poligono.append(primeira_peça)
            if posição_disponivel:
                self.remover_peca(n)
                self.atualizar_estado(pontos_pol)
                return self.verificar_final(0,n)
            else:
                self.remover_peca(n)
                self.atualizar_estado(None)
                return self.verificar_final(1,pontos_pol)
        else:
            self.remover_peca(-1)
            self.atualizar_estado(None)
            return self.verificar_final(2,n)

    def verificar_posicao_disponivel(self, n, x, y, grau_indice):
        pol = nova_lista_completa[n]
        graus = [0, 90, 180, 270]
        grau = graus[grau_indice]
        x1, y1 = transformar_coordenadas_inversa(x, y, self.X, self.X + self.base, self.Y, self.Y - self.altura, self.base, self.altura)

        pontos_pol = [(int(x1) + int(rotate_point(cor[0], cor[1], grau)[0]), int(y1) + int(rotate_point(cor[0], cor[1], grau)[1])) for cor in pol]
        pontos_primeiro = None
        if not CordenadasProibidas:
            poli = nova_lista[n]
            pontos_primeiro = [(int(x1) + int(rotate_point(cor[0], cor[1], grau)[0]), int(y1) + int(rotate_point(cor[0], cor[1], grau)[1])) for cor in poli]
        if any(CordenadaProibida(px, py) for px, py in pontos_pol):
            return False, pontos_pol, None

        if any(any(ponto_dentro_poligono(cor[0], cor[1], pontos_pol, False) for cor in pol) for pol in CordenadasProibidas):
            return False, pontos_pol,None

        return True, pontos_pol, pontos_primeiro

    def atualizar_estado(self, pontos_pol):
        if pontos_pol is not None:
            Matriz_np._criar_poligono(pontos_pol)
            CordenadasProibidas.append(pontos_pol)
        areaPreenchida = np.sum(Matriz_np._matriz)
        
        if areaPreenchida > 0:
            vertices_fecho, area_fecho = calcular_fecho_e_area(CordenadasProibidas)
        else:
            area_fecho = 0
        
        cutting_area = np.array((Matriz_np._matriz), dtype=np.int8)
        pieces = self.preencher_pecas()
        self.state = {'cutting_area': cutting_area, 'pieces': pieces}
        self.fecho = area_fecho

    def preencher_pecas(self):
        pieces = polygon_features(lista_poligonos)
        diferença = self.numero_pecas - len(pieces)
        for _ in range(diferença):
            pieces.append([-1, -1, -1, -1, -1, -1, -1])
        return np.array((pieces), dtype=np.int8)

    def remover_peca(self,n):
        lista_poligonos.pop(n)
        nova_lista.pop(n)
        nova_lista_completa.pop(n)

    def verificar_final(self, n,peca):
        global CordenadasProibidas
        areaPreenchida = np.sum(Matriz_np._matriz)
        area_total = self.base*self.altura
        if areaPreenchida > 0:
            vertices_fecho, area_fecho = calcular_fecho_e_area(CordenadasProibidas)
        else:
            area_fecho = 0
        self.fecho = area_fecho
        desperdicio = area_fecho - areaPreenchida

        if not lista_poligonos:
            self.done = True
           
            
        if n == 0:
            recompensa = round(self.calcular_recompensa_intermediaria(areaPreenchida, len(CordenadasProibidas), self.numero_pecas),2)
            if not self.done:
                self.recompensas.append(round(recompensa,2))
                
            elif self.done:
                self.preenchimentos.append((areaPreenchida / (self.base * self.altura)) * 100)
                self.pecas_colocadas.append(len(CordenadasProibidas))
                
                if (self.epoca % 100) == 0:
                    self.calcular_medias()
                
                Matriz_np.fechar_janela()
                plt.show(block=False)

                self.recompensas.append(round(recompensa,2))
                ultimos_12_elementos = self.recompensas[-12:]

                media = sum(ultimos_12_elementos) / len(ultimos_12_elementos)
                print(f"{round((areaPreenchida / (self.base * self.altura)) * 100, 2)}%, {len(CordenadasProibidas)}/12, R:{ultimos_12_elementos}, SR:{round(sum(ultimos_12_elementos,2))}")
                self.dist = 0
        
                
                
        elif n ==1:
            area = calcular_area(peca)
            recompensa = round(-(3 * (area/area_total)),2)
            #recompensa = max(recompensa, -1)
            #recompensa = 0
            if not self.done:

                self.recompensas.append(recompensa)
            
            elif self.done:
                self.preenchimentos.append((areaPreenchida / (self.base * self.altura)) * 100)
                self.pecas_colocadas.append(len(CordenadasProibidas))
                
                if (self.epoca % 100) == 0:
                    self.calcular_medias()
                
                Matriz_np.fechar_janela()
                plt.show(block=False)

                self.recompensas.append(round(recompensa,2))
                ultimos_12_elementos = self.recompensas[-12:]

                media = sum(ultimos_12_elementos) / len(ultimos_12_elementos)
                print(f"{round((areaPreenchida / (self.base * self.altura)) * 100, 2)}%, {len(CordenadasProibidas)}/12, R:{ultimos_12_elementos}, SR:{round(sum(ultimos_12_elementos),2)}")
                self.dist = 0
        
                
                
        elif n ==2:
            recompensa = -1
            if not self.done:
                self.recompensas.append(recompensa)
                return self.state, float(recompensa), self.done, False, {}
            elif self.done:
                self.preenchimentos.append((areaPreenchida / (self.base * self.altura)) * 100)
                self.pecas_colocadas.append(len(CordenadasProibidas))
                
                if (self.epoca % 100) == 0:
                    self.calcular_medias()
                
                Matriz_np.fechar_janela()
                plt.show(block=False)


                self.recompensas.append(round(recompensa,2))
                ultimos_12_elementos = self.recompensas[-12:]

                media = sum(ultimos_12_elementos) / len(ultimos_12_elementos)
                print(f"{round((areaPreenchida / (self.base * self.altura)) * 100, 2)}%, {len(CordenadasProibidas)}/12, R:{ultimos_12_elementos}, SR:{round(sum(ultimos_12_elementos),2)}")
     
        
                
        return self.state, float(recompensa), self.done, False, {}

    def calcular_recompensa_intermediaria(self, areaPreenchida, num_pecas_colocadas, total_pecas):
        preenchimento_total = areaPreenchida / (self.base * self.altura)
        eficiencia_fecho = areaPreenchida / self.fecho
        
        bonus_pecas = (num_pecas_colocadas / total_pecas) ** 2
        
     
        preenchimento_ajustado = preenchimento_total ** 0.5
        eficiencia_ajustada = eficiencia_fecho ** 0.5
        
        ultimo_pol = CordenadasProibidas[-1]
        area = calcular_area(ultimo_pol)
        recompensa = (0.6 * preenchimento_ajustado + 
                    0.3 * eficiencia_ajustada + 
                    0.1 * bonus_pecas)
        

        #return 2 * recompensa - 1
        #return preenchimento_total - (2 * (1 - eficiencia_fecho))
        #return (area/self.area) - (1-eficiencia_fecho)
        #return preenchimento_total
        #return eficiencia_fecho
        return (area/self.area)
        #return 5 * (math.e**(preenchimento_total - 1) - 0.6)

    def calcular_medias(self):

        print(f"Media dos ultimos 100 preenchimentos : {round((sum(self.preenchimentos)) / (len(self.preenchimentos)), 2)}%, {round((sum(self.pecas_colocadas)) / (len(self.pecas_colocadas)), 2)}/12")
        print(f"Media das ultimas 100 recompensas: {round((sum(self.recompensas)) / (len(self.recompensas)), 2)}")
        self.MediaPreenchimentos.append(round((sum(self.preenchimentos)) / (len(self.preenchimentos)), 2))
        self.MediaRecompensas.append(round((sum(self.recompensas)) / (len(self.recompensas)), 2))
        self.Mpecas_colocadas.append(round((sum(self.pecas_colocadas)) / (len(self.pecas_colocadas)), 2))
 
        
        if len(self.MediaPreenchimentos) > 10:
            self.MediaPreenchimentos.pop(0)
            self.MediaRecompensas.pop(0)
            self.Mpecas_colocadas.pop(0)
        
        print(f"ultimas 10 preenchimentos{self.MediaPreenchimentos}")
        print(f"ultimas 10 recompensas{self.MediaRecompensas}")
        self.preenchimentos = []
        self.recompensas = []
        self.pecas_colocadas = []
                    
            
    
    def instancias(self,Base, Altura, escala):
        global Escala
        if self.pecas == 'fu':
            lados = [102,117]
            i = random.randint(0, 1)
            j = 1 - i
            self.base = 102
            self.altura = 117
            Escala = 3
            
        elif self.pecas == 'jackobs1':
            self.base = 400
            self.altura = 130
            Escala = 10
            
        elif self.pecas == 'jackobs2':
            self.base = 350
            self.altura = 141
            Escala = 5
            
        elif self.pecas == 'shapes0':
            self.base = 315
            self.altura = 200
            Escala = 5
            
        elif self.pecas == 'shapes1':
            self.base = 160
            self.altura = 236
            Escala = 4
            
        elif self.pecas == 'shapes2':
            self.base = 270
            self.altura = 150
            Escala = 10

        elif self.pecas == 'dighe1':
            self.base = 414
            self.altura = 300
            Escala = 3
            
        elif self.pecas == 'dighe2':
            self.base = 300
            self.altura = 405
            Escala = 3

        elif self.pecas == 'albano':
            self.base = 203
            self.altura = 98
            Escala = (1/50)

        elif self.pecas == 'dagli':
            self.base = 300
            self.altura = 325
            Escala = 5

        elif self.pecas == 'mao':
            self.base = 206
            self.altura = 255
            Escala = (1/10)

        elif self.pecas == 'marques':
            self.base = 166
            self.altura = 208
            Escala = 2
            
        elif self.pecas == 'shirts':
            self.base = 80
            self.altura = 126
            Escala = 2
            
        elif self.pecas == 'swim':
            self.base = 114
            self.altura = 130
            Escala = (1/50)
            
        elif self.pecas == 'trousers':
            self.base = 245
            self.altura = 79
            Escala = 1

        else:
            self.base = Base
            self.altura = Altura    
            Escala = escala

    def seed(self, seed=None):
        self.np_random, seed = gym.utils.seeding.np_random(seed)
        return [seed]
    
    def calcular_recompensa(self):
        global CordenadasProibidas
        n = 0
        sub_base = self.base//3
        sub_altura = self.altura//2
        sub_area = sub_base * sub_altura
        preenchimento = np.sum(Matriz_np._matriz) 
        g1,g2,g3,g4,g5,g6 = dividir_em_submatrizes(Matriz_np._matriz, self.base, self.altura)
        
        for g in [g1,g2,g3,g4,g5,g6]:
            if g > 0:
                n+=1


        pt = g1+g2+g3+g4+g5+g6
    
        if self.fecho>0:
            preenchimento_fecho = (float(preenchimento)/float(self.fecho))
        else:
            preenchimento_fecho = 0
        if n >0:  
            return float((preenchimento/(self.altura*self.base)))
        else:
            return float(0)
        
        
        
def dividir_em_submatrizes(matriz, base, altura):
    altura_sub = altura // 2
    base_sub = base // 3
    
    preenchimentos = []
    
    # Itera sobre as submatrizes
    for i in range(2):  # Para as linhas
        for j in range(3):  # Para as colunas
            # Extrai a submatriz
            submatriz = matriz[i*altura_sub:(i+1)*altura_sub, j*base_sub:(j+1)*base_sub]
            
            # Calcula o preenchimento da submatriz
            preenchimento = np.sum(submatriz) / (altura_sub * base_sub)
            preenchimentos.append(preenchimento)
    
    return preenchimentos

import matplotlib.pyplot as plt
fig, ax = plt.subplots()

def plot2(scores, mean_scores, fig, ax):
    # Limpa os eixos
    ax.clear()

    # Configura o gráfico
    ax.set_title('Training...')
    ax.set_xlabel('Épocas')
    ax.set_ylabel('Preenchimento, em %')
    ax.plot(mean_scores, label='Preenchimento')
    ax.plot(scores, label='Recompensa')
    ax.set_ylim(ymin=0)
    ax.text(len(scores)-1, scores[-1], str(scores[-1]))
    ax.text(len(mean_scores)-1, mean_scores[-1], str(mean_scores[-1]))
    ax.legend()
    fig.canvas.draw_idle()
    # Atualiza a figura
    
log_dir = "./ppo_csp_tensor/"

if not os.path.exists(log_dir):
    os.makedirs(log_dir)

class TensorboardCallback(BaseCallback):
    def __init__(self, verbose=0):
        super(TensorboardCallback, self).__init__(verbose)

    def _on_step(self) -> bool:
        # Acessa o ambiente e extrai a variável personalizada
        if self.locals.get('env') and hasattr(self.locals['env'], 'filling'):
            filling = self.locals['env'].preenchimentos
            if isinstance(filling, list) and filling:
                # Calcula a média da lista, por exemplo
                avg_filling = sum(filling) / len(filling)
                self.logger.record('filling_avg', avg_filling)
            else:
                # Caso filling não seja uma lista ou esteja vazio
                self.logger.record('filling_avg', 0)
        
        return True

def make_env(env_class, rank, seed=0):
    def _init():
        env = env_class()
        env.seed(seed + rank)
        return env
    return _init

def linear_schedule(initial_value, final):
    def func(progress):
        if initial_value * (progress) <= final:
            return final
        else:
            return initial_value * (progress)
    return func

class CustomCNNExtractor(BaseFeaturesExtractor):
    def __init__(self, observation_space, features_dim=256):
        super(CustomCNNExtractor, self).__init__(observation_space, features_dim)
        
        # CNN para a área de corte
        self.cnn_cutting_area = nn.Sequential(
            nn.Conv2d(1, 32, kernel_size=3, stride=1, padding=1),
            nn.ReLU(),
            nn.Conv2d(32, 64, kernel_size=3, stride=1, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2),
            nn.Conv2d(64, 128, kernel_size=3, stride=1, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2),
            nn.Flatten()
        )

        # Ajustar as camadas da MLP para as `pieces`
        n_pieces_features = observation_space['pieces'].shape[1] * observation_space['pieces'].shape[0]

        self.mlp_pieces = nn.Sequential(
            nn.Linear(n_pieces_features, 128),
            nn.ReLU(),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Flatten()
        )

        # Combinar as saídas das duas redes
        with torch.no_grad():
            n_flatten_cutting_area = self.cnn_cutting_area(
                torch.as_tensor(observation_space['cutting_area'].sample()[None, None]).float()
            ).shape[1]

            n_flatten_pieces = self.mlp_pieces(
                torch.as_tensor(observation_space['pieces'].sample().flatten()[None]).float()
            ).shape[1]

        self.linear = nn.Sequential(
            nn.Linear(n_flatten_cutting_area + n_flatten_pieces, features_dim),
            nn.ReLU()
        )

    def forward(self, observations):
        cutting_area = observations['cutting_area'].unsqueeze(1)  # Adicionar canal de cor para a CNN
        pieces = observations['pieces'].view(observations['pieces'].shape[0], -1)  # Achatar as peças
        
        cutting_area_features = self.cnn_cutting_area(cutting_area)
        pieces_features = self.mlp_pieces(pieces)
        
        combined_features = torch.cat((cutting_area_features, pieces_features), dim=1)
        return self.linear(combined_features)

if __name__ == '__main__':
    env = env2KP
    num_envs = 32
    envs = SubprocVecEnv([make_env(env, i) for i in range(num_envs)])


    policy_kwargs = dict(
        features_extractor_class=CustomCNNExtractor,
        net_arch=dict(pi=[64, 64], vf=[64, 64]),
        activation_fn=nn.ReLU
    )

    model = PPO('MultiInputPolicy', envs, verbose=1, tensorboard_log=log_dir,
                n_steps=60,  # Ajustado para ser múltiplo de batch_size
                batch_size=120,  # Aumentado para melhor estabilidade
                n_epochs=15,  # Reduzido para treinamento mais rápido
                gamma=0.95,  # Aumentado para considerar recompensas de longo prazo
                learning_rate=linear_schedule(5e-4, 0),
                ent_coef=0.01,
                #policy_kwargs=policy_kwargs,
                clip_range=linear_schedule(0.4, 0))




    model.learn(total_timesteps=10_000_000)
    caminho_para_salvar = "/home/fsilvestre/Cutting_Stock_Problem"
    os.makedirs(os.path.dirname(caminho_para_salvar), exist_ok=True)
    model.save(caminho_para_salvar)

    

  

