import sys
import os

# Adiciona o diretório raiz do projeto ao sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))

import numpy as np
from PIL import Image, ImageDraw
import matplotlib.pyplot as plt
import turtle
import math
from scipy.spatial import ConvexHull

# Ajuste o caminho de importação para o módulo Botao
from Cutting_Stock_Problem.Ambiente.Main.botao import Botao

#import pyautogui

class CSP():
    def __init__(self,dataset='fu',Base=None,Altura=None,Escala=None,render=False,plot=True, x=-200, y=200):
        
        self.x = x
        self.y = y
        
        self.render = render
        self.plot = plot
        
        self.pecas = dataset
        if Base == None and Altura == None:
            self.instancias(Base,Altura,Escala)
        else:
            self.base = Base
            self.altura = Altura
            self.Escala = Escala
        self.area = self.base * self.altura

        
        
        self.lista_original = ler_poligonos('/home/fsilvestre/Cutting_Stock_Problem/' + self.pecas +'.dat')
        #self.lista_original = [ajustar_poligono(pol) for pol in self.lista_original]

        self.nova_lista_original, self.nova_lista_completa_original = tratar_lista(self.lista_original, self.Escala)
        #self.nova_lista_original = [ajustar_poligono(pol) for pol in self.nova_lista_original]
        #self.nova_lista_completa_original = [ajustar_poligono(pol) for pol in self.nova_lista_completa_original]

        self.max_pecas = len(self.lista_original)
        
        self.lista = ler_poligonos('/home/fsilvestre/Cutting_Stock_Problem/' + self.pecas +'.dat')
        #self.lista = [ajustar_poligono(pol) for pol in self.lista]

        self.nova_lista, self.nova_lista_completa = tratar_lista(self.lista, self.Escala)
        

        #self.nova_lista = [ajustar_poligono(pol) for pol in self.nova_lista]
        #self.nova_lista_completa = [ajustar_poligono(pol) for pol in self.nova_lista_completa]

        self.pecas_posicionadas = []
        self.indices_pecas_posicionadas = []
        self.lista_removida = []
        self.cordenadas_pecas = []
        
        self.cordenadas_area = ([self.x, self.y],[(self.x + self.base), self.y],[(self.x + self.base), (self.y - self.altura)],[self.x, (self.y - self.altura)])
        if self.plot:
            self._img = Image.new('L', (self.base, self.altura), 0)
            self._matriz = np.array(self._img)
            self._fig, self._ax = plt.subplots()
        
        if self.render:
            turtle.tracer(0)
            self.botao_remover = Botao(self.x + 350, self.y + 300, 200, 150)
            self.turtle_fecho = turtle.Turtle()
            self.turtle_dados = turtle.Turtle()
            self.lista_turtle = []
            self.renderizar()
            
        self.area_ocupada = 0
        self.fecho = 0
        self.vertices_fecho = []
            
    def acao(self,peca,x,y,grau_indice,flip,verificado=False):
        if peca < len(self.lista):
            if not verificado: 
                #print("Nao Verificado")       
                pol_verificacao = self.nova_lista_completa[peca]
                pol_posicionar = self.nova_lista[peca]
                if flip:
                    pol_posicionar = flip_polygon(pol_posicionar)
                    pol_verificacao = flip_polygon(pol_verificacao)
                graus = [0, 90, 180, 270]
                grau = graus[grau_indice]
        
                pontos_verificacao = [(int(x) + int(rotate_point(cor[0], cor[1], grau)[0]), int(y) + int(rotate_point(cor[0], cor[1], grau)[1])) for cor in pol_verificacao]
                pontos_posicionar = [(int(x) + int(rotate_point(cor[0], cor[1], grau)[0]), int(y) + int(rotate_point(cor[0], cor[1], grau)[1])) for cor in pol_posicionar]
                
                if self.posicao_disponivel(pontos_verificacao):
                    self.indices_pecas_posicionadas.append([x,y,grau_indice,self.nova_lista_original.index(pol_posicionar)])
                    self.pecas_posicionadas.append(pontos_posicionar)
                    if self.render:
                        p = turtle.Turtle()
                        self.posicionar_turtle(pontos_posicionar,p)
                        p.hideturtle()
                        self.lista_turtle.append(p)
                    self.remover_lista(peca) 
                    self.atualizar_matriz()
                    self.atualizar_dados()
            else:
                    #print("entrou")
                    pol_posicionar = self.nova_lista[peca]
                    if flip:
                        pol_posicionar = flip_polygon(pol_posicionar)

                    graus = [0, 90, 180, 270]
                    grau = graus[grau_indice]
            
                    pontos_posicionar = [(int(x) + int(rotate_point(cor[0], cor[1], grau)[0]), int(y) + int(rotate_point(cor[0], cor[1], grau)[1])) for cor in pol_posicionar]
                    
                    if True:
                        self.indices_pecas_posicionadas.append([x,y,grau_indice,self.nova_lista_original.index(pol_posicionar)])
                        self.pecas_posicionadas.append(pontos_posicionar)
                        if self.render:
                            p = turtle.Turtle()
                            self.posicionar_turtle(pontos_posicionar,p)
                            p.hideturtle()
                            self.lista_turtle.append(p)
                        self.remover_lista(peca) 
                        self.atualizar_matriz()
                        self.atualizar_dados()

        else:
            print("Error")
            return False
     

    def atualizar_dados(self):
        self.area_ocupada = sum([calcular_area(pol) for pol in self.pecas_posicionadas])
        self.vertices_fecho, self.fecho = calcular_fecho_e_area(self.pecas_posicionadas)
        self.atualizar_matriz()
        
        if self.render:
            turtle.tracer(0)
            self.turtle_fecho.clear()
            self.turtle_dados.clear()
            
            self.turtle_dados.penup()
            self.turtle_dados.goto(-925,200)
            self.turtle_dados.pendown()
            self.turtle_dados.write(f"Area Desperdiçada = {round(((self.fecho - self.area_ocupada)/self.area)*100, 2)}", font=("Arial",16, "normal"))
            self.turtle_dados.penup()
            self.turtle_dados.goto(-925,350)
            self.turtle_dados.pendown()

            self.turtle_dados.write(f"Preenchimento total: {(round((self.area_ocupada/self.area)*100,2))}%", font=("Arial",16, "normal"))
            self.turtle_dados.penup()
            self.turtle_dados.goto(-925,300)
            self.turtle_dados.pendown()
            self.turtle_dados.write(f"Peças Disponiveis: {len(self.lista)}", font=("Arial",16, "normal"))
            self.turtle_dados.hideturtle()
            desenhar_fecho(self.vertices_fecho, self.turtle_fecho)
              
    def posicionar_turtle(self, poli, p):
        p.begin_fill()
        p.fillcolor("orange")
        turtle.tracer(0)
        p.penup()
        p.goto(poli[0])
        p.pendown()                
        for cor in poli:
            p.goto(cor)

        p.goto(poli[0])
        p.end_fill()
        p.hideturtle()

        
    def CordenadaProibida(self,x,y):
        if (not ponto_dentro_poligono(x, y, self.cordenadas_area, True)):
            return True
        else:
            for pol in self.pecas_posicionadas:
                if ponto_dentro_poligono(x, y, pol,False):
                    return True
    
            return False 
    
    def remover_lista(self,peca):
        self.lista_removida.append([self.lista[peca], self.nova_lista[peca], self.nova_lista_completa[peca]])
        
        self.lista.pop(peca)
        self.nova_lista.pop(peca)
        self.nova_lista_completa.pop(peca)
  
    def remover_da_area(self):
        
        peca = -1
        if self.pecas_posicionadas:
            self.pecas_posicionadas.pop(peca)
            self.indices_pecas_posicionadas.pop(peca)
            
            self.lista.append(self.lista_removida[peca][0])           
            self.nova_lista.append(self.lista_removida[peca][1])
            self.nova_lista_completa.append(self.lista_removida[peca][2])
            
            self.restaurar_ordem_original()
            
            self.lista_removida.pop(peca)
            
            if self.render:
                self.lista_turtle[peca].clear()
                self.lista_turtle.pop(peca)
                
            self.atualizar_dados()
        
            
    def restaurar_ordem_original(self):
        # Função auxiliar para lidar com itens repetidos
        def restaurar_lista(lista_original, lista_atual):
            lista_restaurada = []
            # Copiar lista atual para contagem correta de repetidos
            lista_atual_copy = lista_atual[:]
            
            # Percorrer lista_original de trás para frente
            for item in reversed(lista_original):
                if item in lista_atual_copy:
                    # Inserir o item no início da lista restaurada para manter a ordem
                    lista_restaurada.insert(0, item)
                    lista_atual_copy.remove(item)  # Remove uma ocorrência do item
                    
            return lista_restaurada

        # Reorganizando self.lista com base em self.lista_original, lidando com repetições
        self.lista = restaurar_lista(self.lista_original, self.lista)
        
        # Reorganizando self.nova_lista com base em self.nova_lista_original
        self.nova_lista = restaurar_lista(self.nova_lista_original, self.nova_lista)
        
        # Reorganizando self.nova_lista_completa com base em self.nova_lista_completa_original
        self.nova_lista_completa = restaurar_lista(self.nova_lista_completa_original, self.nova_lista_completa)
      
    def posicao_disponivel(self,pontos_pol):

        pontos = pontos_no_poligono(pontos_pol)

        if any(self.CordenadaProibida(px, py) for px, py in pontos):
            #print(1)
            return False
        
        if any(self.CordenadaProibida(px, py) for px, py in pontos_pol):
            #print(2)
            return False

        if any(any(ponto_dentro_poligono(cor[0], cor[1], pontos_pol, False) for cor in pol) for pol in self.pecas_posicionadas):
            #print(3)
            return False
        
        #for pol in self.pecas_posicionadas:
            #p = pontos_no_poligono(pol)
            #for x,y in p:
                #if ponto_dentro_poligono(x,y,pontos_pol, False):
                    #print(1)
                    #return False        
        
        return True
  
    def renderizar(self):
        turtle.tracer(0)
        t  = turtle.Turtle()
        t.color('red')
        t.begin_fill()
        t.penup()
        t.goto(self.x, self.y)
        t.pendown()
        
        for cor in self.cordenadas_area:
            t.goto(cor)

        t.goto(self.cordenadas_area[0])
        t.end_fill()
        t.hideturtle()

    def _plot(self):
        self._ax.clear()
        self._ax.imshow(self._matriz)
        plt.pause(0.001)
        
    def atualizar_matriz(self):          
        self._img = Image.new('L', (self.base, self.altura), 0)
        for poligono in self.pecas_posicionadas:
            poligono_transformado = [transformar_coordenadas(x, y, self.x, (self.x + self.base), (self.y - self.altura), self.y, self.base, self.altura) for x, y in poligono]
            ImageDraw.Draw(self._img).polygon(poligono_transformado, outline=1, fill=1)
        self._matriz = np.array(self._img)
        
        if self.plot:
            self._plot()
        
    def instancias(self,Base, Altura, escala):
        if self.pecas == 'fu':
            self.base = 340
            self.altura = 380
            self.Escala = 10
            
        elif self.pecas == 'jackobs1':
            self.base = 400
            self.altura = 130
            self.Escala = 10
            
        elif self.pecas == 'jackobs2':
            self.base = 350
            self.altura = 141
            self.Escala = 5
            
        elif self.pecas == 'shapes0':
            self.base = 315
            self.altura = 200
            self.Escala = 5
            
        elif self.pecas == 'shapes1':
            self.base = 160
            self.altura = 236
            self.Escala = 4
            
        elif self.pecas == 'shapes2':
            self.base = 270
            self.altura = 150
            self.Escala = 10

        elif self.pecas == 'dighe1':
            self.base = 415
            self.altura = 300
            self.Escala = 3
            
        elif self.pecas == 'dighe2':
            self.base = 300
            self.altura = 405
            self.Escala = 3

        elif self.pecas == 'albano':
            self.base = 203
            self.altura = 98
            self.Escala = (1/50)

        elif self.pecas == 'dagli':
            self.base = 300
            self.altura = 325
            self.Escala = 5

        elif self.pecas == 'mao':
            self.base = 206
            self.altura = 255
            self.Escala = (1/10)

        elif self.pecas == 'marques':
            self.base = 166
            self.altura = 208
            self.Escala = 2
            
        elif self.pecas == 'shirts':
            self.base = 80
            self.altura = 126
            self.Escala = 2
            
        elif self.pecas == 'swim':
            self.base = 114
            self.altura = 130
            self.Escala = (1/50)
            
        elif self.pecas == 'trousers':
            self.base = 245
            self.altura = 79
            self.Escala = 1

        else:
            self.base = Base
            self.altura = Altura    
            self.Escala = escala

    def draw_click(self,x1,y1, grau_indice,n,flip):
        graus = [0,90,180,270]
        grau = graus[grau_indice]
        x, y = x1,y1
        x_escala = 0
        y_escala = 0
        if len(self.lista) > 0 and n < len(self.lista):          
            pol = self.lista[n]
            if flip:
                pol = flip_polygon(pol)
            turtle.tracer(0)
            turtle.clear()
            turtle.begin_fill()
            turtle.fillcolor("green")
            turtle.penup()
            x_rot, y_rot = rotate_point(pol[0][0], pol[0][1], grau)
            turtle.goto(x+x_rot*self.Escala-x_escala, -y+y_rot*self.Escala+y_escala)
            turtle.pendown()
            for cor in pol:

                x_rot, y_rot = rotate_point(cor[0], cor[1], grau)
                turtle.goto(x+x_rot*self.Escala-x_escala, -y+y_rot*self.Escala+y_escala)
            turtle.goto(x+x_rot*self.Escala-x_escala, -y+y_rot*self.Escala+y_escala)
            turtle.end_fill()
            turtle.hideturtle()
            turtle.tracer(1)

      
    def mais_g(self,x,y):
        self.g += 1
        if self.g == 4:
            self.g = 0
        
    def mais_n(self,x,y):
        self.n += 1
        if self.n == len(self.nova_lista):
            self.n=0
         
    def on_click(self,x,y):
 
        if self.render:
            if self.botao_remover.clicou(x,y):
                self.n = 0
                self.remover_da_area()
                return
        
        numero_pecas = len(self.lista)
        self.acao(self.n,x,y,self.g, False)
        if numero_pecas > len(self.lista):
            self.n=0
            
    def click(self):
        self.loop = True
        self.n = 0
        self.g = 0
        while self.loop:
    
            turtle.onscreenclick(self.mais_g,3)
            turtle.onscreenclick(self.mais_n,2)
            turtle.onscreenclick(self.on_click)
            #x1,y1 = pyautogui.position()
            x1 = 963
            y1 = 533
             
            self.draw_click(x1,y1,self.g,self.n, False)
            
            if len(self.nova_lista)==0:
                self.loop = False
                turtle.clear()
                break

def flip_polygon(polygon, axis='y'):
    """
    Flipa o polígono ao longo do eixo x ou y.
    
    Args:
        polygon: Lista de tuplas representando os vértices [(x1, y1), (x2, y2), ...].
        axis: Eixo de reflexão, 'x' para flipar ao longo do eixo x, ou 'y' para o eixo y.
        
    Returns:
        Lista de vértices do polígono flipado.
    """
    if axis == 'x':
        # Reflete ao longo do eixo x: inverte o sinal da coordenada y
        return [(x, -y) for x, y in polygon]
    elif axis == 'y':
        # Reflete ao longo do eixo y: inverte o sinal da coordenada x
        return [(-x, y) for x, y in polygon]
    else:
        raise ValueError("Eixo inválido. Use 'x' ou 'y'.")

from shapely.geometry import Point, Polygon
import random

def pontos_no_poligono(pontos_poligono, num_pontos=100):
    # Criar o objeto Polígono a partir dos pontos
    poligono = Polygon(pontos_poligono)

    # Ponto central (mesmo cálculo que já foi feito)
    xs = [x for x, y in pontos_poligono]
    ys = [y for x, y in pontos_poligono]
    centro_x = sum(xs) / len(xs)
    centro_y = sum(ys) / len(ys)
    ponto_central = (centro_x, centro_y)

    # Gerar pontos aleatórios dentro do polígono
    pontos_internos = []
    min_x, min_y, max_x, max_y = poligono.bounds
    while len(pontos_internos) < num_pontos:
        # Gerar coordenadas aleatórias dentro dos limites do polígono
        random_point = Point(random.uniform(min_x, max_x), random.uniform(min_y, max_y))
        if poligono.contains(random_point):
            pontos_internos.append((random_point.x, random_point.y))

    return [ponto_central] + pontos_internos


def ajustar_poligono(poligono):
    # Encontrar o primeiro vértice
    if poligono[0][0] == 0 and poligono[0][1] == 0:
        return poligono
    else:
        x_min = poligono[0][0]
        y_min = poligono[0][1]
        
        # Transladar todos os vértices para garantir que o primeiro ponto seja (0,0)
        poligono_ajustado = [(x - x_min, y - y_min) for (x, y) in poligono]
        
        return poligono_ajustado 

    
def rotate_point(x, y, angle):
    rad = math.radians(angle)
    x_new = x * math.cos(rad) - y * math.sin(rad)
    y_new = x * math.sin(rad) + y * math.cos(rad)
    return x_new, y_new
def transformar_coordenadas_inversa(x, y, x_min, x_max, y_min, y_max, width, height):
    x_scaled = (x / width) * (x_max - x_min) + x_min
    y_scaled = ((height - y) / height) * (y_max - y_min) + y_min
    return x_scaled, y_scaled
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
def tratar_lista(lista_poligonos, Escala):
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
def transformar_coordenadas(x, y, x_min, x_max, y_min, y_max, width, height):
    x_scaled = (x - x_min) / (x_max - x_min) * width
    y_scaled = (y - y_min) / (y_max - y_min) * height
    y_scaled = height - y_scaled

    return x_scaled, y_scaled
def calcular_area(vertices):
    num_vertices = len(vertices)
    area = 0

    for i in range(num_vertices):
        x1, y1 = vertices[i]
        x2, y2 = vertices[(i + 1) % num_vertices]  # Próximo vértice (considerando o último vértice conectado ao primeiro)

        area += (x1 * y2) - (x2 * y1)

    return abs(area) / 2
def calcular_fecho_e_area(poligonos):
    if poligonos:
        # Concatenar todos os vértices de todos os polígonos
        vertices = np.concatenate(poligonos, axis=0)
        
        # Calcular o fecho convexo
        fecho = ConvexHull(vertices)
        
        # Obter os vértices do fecho convexo
        vertices_fecho = vertices[fecho.vertices]
        
        # Calcular a área do fecho convexo
        area = fecho.volume
        
        return vertices_fecho, area
    else:
        return None, 0
    
def desenhar_fecho(vertices_fecho, p):
    if vertices_fecho is not None:
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
        
#ambiente = CSP(dataset='fu', render=True, plot=True)
#ambiente.click()

