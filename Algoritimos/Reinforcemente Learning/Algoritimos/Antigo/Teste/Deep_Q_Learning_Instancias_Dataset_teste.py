import turtle 
from geradores import Area, Gerador, posicionar, on_click
from matriz import Matriz
from botao import Botao
import pyautogui
import math
import time
import gymnasium as gym
from gymnasium import spaces
import numpy as np
from stable_baselines3 import DQN
from stable_baselines3.common.env_checker import check_env
import matplotlib.pyplot as plt
from stable_baselines3.common.callbacks import BaseCallback
import turtle
from geradores import Area
from botao import Botao
import pyautogui
import math
import numpy as np
from scipy.spatial import ConvexHull
import time
import random
import csv
import copy
import matplotlib.pyplot as plt
from IPython import display
from matriz_vertices import Matriz

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
    print(f"Número total de polígonos: {num_poligonos}")

    poligonos = []
    i = 1  # Começa a leitura a partir da segunda linha

    while i < len(linhas):
        # Verifica se a linha não está vazia
        if linhas[i].strip():
            try:
                # Lê o número de vértices
                num_vertices = int(linhas[i].strip())
                print(f"Lendo polígono com {num_vertices} vértices")  # Depuração
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
        print(f'Todos os {num_poligonos} poligonos foram lidos com sucesso!')
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
def on_click2(x, y, grau,fecho_turtle, dados_turtle, botao, n,colunas,base,altura,adicional2,verificado):
    global CordenadasProibidas
    global lista_colunas
    global fecho
    global AreaPreenchida
    global turtles
    global l
    global lista_area_fechos
    global recompensa
    if len(nova_lista) > n:
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
                lista_poligonos.append([-1])
                l = 0

                return area_fecho
def transformar_coordenadas_inversa(x, y, x_min, x_max, y_min, y_max, width, height):
    x_scaled = (x / width) * (x_max - x_min) + x_min
    y_scaled = ((height - y) / height) * (y_max - y_min) + y_min

    return x_scaled, y_scaled



    
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
        
        self.instancias(Base,Altura,escala)

        instancias = self.pecas + '.dat'
        lista_poligonos = ler_poligonos(instancias)
        random.shuffle(lista_poligonos)

        max_vertices = max(len(p) for p in lista_poligonos)

            # Espaço de observação: um dicionário que inclui a matriz da área de corte e a lista de peças restantes
        self.observation_space = spaces.Dict({
            'cutting_area': spaces.Box(low=0, high=1, shape=(self.altura, self.base), dtype=np.float32),
            'pieces': spaces.Box(low=0, high=1, shape=(len(lista_poligonos),max_vertices, 2), dtype=np.float32)  # Usando apenas um polígono por vez
        })
            
            # Espaço de ação: escolha de uma peça e sua posição na área de corte
        self.action_space = spaces.Discrete(len(lista_poligonos) * 4 * self.base * self.altura)

        self.done = False
     
    def reset(self, seed=None):
        global X,Y,tempoF,grau,CordenadasProibidas,lista_colunas,fecho,lista_area_fechos,AreaPreenchida,l,fecho_turtle,dados_turtle,turtles,lista_poligonos,t,inicio,tempo_anterior,fim,Bot,nova_lista,nova_lista_completa,base,altura
        global Escala,Matriz_np,Plot,Passo,dado,dado1,bases, PC1, PC2, PC3, PP1, PP2, va,acao 

        self.recompensa = -6
        self.epoca+=1
        self.acao = 0
        print("reset")
        self.done = False
        if seed is not None:
            np.random.seed(seed)
        Area.CordenadasProibidas = []
        turtle.clearscreen()
        turtle.reset()
        nome = turtle.Turtle()
        acao = turtle.Turtle()
        turtle.tracer(0)
        nome.speed(100)
        nome.penup()
        nome.goto(-100,400)
        nome.pendown
        nome.write(f"{self.pecas}", font=("Arial",16, "normal"))
        nome.penup()
        nome.goto(-150,350)
        nome.pendown
        nome.write(f"epoca: {self.epoca}", font=("Arial",16, "normal"))
        nome.hideturtle()

        
        grau = 0
        CordenadasProibidas = []
        lista_colunas = []
        fecho = []
        lista_area_fechos = []
        AreaPreenchida = 0
        l = 0

        fecho_turtle = turtle.Turtle()  
        dados_turtle = turtle.Turtle()
        dado = turtle.Turtle()
        dado1 = turtle.Turtle()
        turtles = []

        t = turtle.Turtle()
        turtle.tracer(0)
        t.penup()
        t.goto(-900,375)
        t.pendown()
        turtle.tracer(1)

        Bot = Botao((self.X+self.base+200), 300, 200, 150)

        turtle.tracer(0)
        
        Matriz_np = Matriz(self.base,self.altura,self.X,self.Y)
        area = Area(self.base,self.altura,6)
        area.area(self.X, self.Y)
        area.greedy(self.X, self.Y)
        instancias = self.pecas + '.dat'
        lista_poligonos = ler_poligonos(instancias)
        #random.shuffle(lista_poligonos)
        nova_lista, nova_lista_completa = tratar_lista(lista_poligonos)

        max_vertices = 4
        pieces = np.zeros((len(lista_poligonos), max_vertices, 2), dtype=np.float32)
        for i, poligono in enumerate(lista_poligonos):
            num_vertices = len(poligono)
            # Normaliza as coordenadas dos vértices entre 0 e 1
            normalized_poligono = np.array(poligono) / [self.base, self.altura]
            
            # Preenche com zeros se tiver menos vértices que max_vertices
            if num_vertices < max_vertices:
                padded_poligono = np.zeros((max_vertices, 2), dtype=np.float32)
                padded_poligono[:num_vertices] = normalized_poligono
                pieces[i] = padded_poligono
            else:
                pieces[i] = normalized_poligono[:max_vertices]

        
        cutting_area = Matriz_np._matriz

            # Retorna o estado de observação como um dicionário
        self.state = {
            'cutting_area': cutting_area,
            'pieces': np.array(pieces)
        }
            
        return self.state, {}

    
    
    def step(self, action):
        global acao,recompensa,fig,ax,Matriz_np
        recompensa = self.recompensa
        acao.clear()
        self.acao += 1
        turtle.tracer(0)
        acao.penup()
        acao.goto(-50,350)
        acao.pendown
        acao.write(f"ação: {self.acao}", font=("Arial",16, "normal"))
        acao.hideturtle()
        fim = time.time()
        tempo_execucao = fim - inicio

        if int(tempo_execucao) != int(self.tempo_anterior):
            tempo(tempo_execucao)
            self.tempo_anterior = int(tempo_execucao)
            
        n = action // (4 * self.base * self.altura)
        remaining = action % (4 * self.base * self.altura)
        grau_indice = remaining // (self.base * self.altura)
        remaining = remaining % (self.base * self.altura)
        x = remaining // self.altura
        y = remaining % self.altura
        
        if n>=len(lista_poligonos):
            n = 0
        graus = [0,90,180,270]
        grau = graus[grau_indice]
        x_transformado, y_transformado = transformar_coordenadas_inversa(x, y, self.X, self.X + self.base, self.Y, self.Y - self.altura, self.base, self.altura)
        for i in range(10):
            draw_piece(x_transformado,y_transformado,grau,n)
        tamanho = len(nova_lista)
        res = on_click2(x_transformado, y_transformado, grau, fecho_turtle, dados_turtle, Bot, n, 3, self.base, self.altura, 100, False)
        if res is not None:
            self.fecho = res

        
        cutting_area = Matriz_np._matriz
        max_vertices = 4
        pieces = np.zeros((len(lista_poligonos), max_vertices, 2), dtype=np.float32)
        for i, poligono in enumerate(lista_poligonos):
            num_vertices = len(poligono)
            # Normaliza as coordenadas dos vértices entre 0 e 1
            normalized_poligono = np.array(poligono) / [self.base, self.altura]
            
            # Preenche com zeros se tiver menos vértices que max_vertices
            if num_vertices < max_vertices:
                padded_poligono = np.zeros((max_vertices, 2), dtype=np.float32)
                padded_poligono[:num_vertices] = normalized_poligono
                pieces[i] = padded_poligono
            else:
                pieces[i] = normalized_poligono[:max_vertices]

        # Atualiza o estado de observação como um dicionário
        self.state = {
            'cutting_area': cutting_area,
            'pieces': np.array(pieces)  # Mantém o tamanho fixo
        }

        reward = self.calcular_recompensa()
        print(reward)
        self.recompensa = reward
        
        if lista_poligonos:
            if len(nova_lista) == tamanho:
                #self.erro += 1
                #if self.erro == 5:
                    self.erro = 0
                    lista_poligonos.pop(n)
                    nova_lista.pop(n)
                    nova_lista_completa.pop(n)
                    lista_poligonos.append([-1])
                    return self.state, float(-10), self.done, False, {}
            else:
                return self.state, float(reward), self.done, False, {}
                 
        else:
            self.done = True
            if self.epoca >0:
                # Exemplo de uso:
                #plt.ion()  # Liga o modo interativo
                areaPreenchida = np.sum(Matriz_np._matriz) 
                self.preenchimentos.append((areaPreenchida/(self.base*self.altura))*100)
                self.recompensas.append(-1 * self.recompensa)
                plot2(self.recompensas, self.preenchimentos ,fig, ax)
                #plt.ioff()  # Desliga o modo interativo
                Matriz_np.fechar_janela()
                plt.show(block=False)
            return self.state, float(reward), self.done, False, {}
            
    
    def instancias(self,Base, Altura, escala):
        global Escala
        if self.pecas == 'fu':
            self.base = 340
            self.altura = 380
            Escala = 10
            
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

    def calcular_recompensa(self):
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
    
        if n >0:  
            return float((-6 + (pt)) - (float(self.fecho-preenchimento)/float(sub_area*n)))
        else:
            return -6
        
        
        
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
    



env = env2KP()
model = DQN('MultiInputPolicy' ,env,learning_rate=0.003,buffer_size=1000,exploration_final_eps=0.1, verbose=1)
model.learn(total_timesteps=10000)
from stable_baselines3.common.env_checker import check_env
check_env(env)

    

