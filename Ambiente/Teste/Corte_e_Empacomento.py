import numpy as np
from PIL import Image, ImageDraw
import matplotlib.pyplot as plt
from geradores import Area
import turtle
import math
from scipy.spatial import ConvexHull
#import pyautogui
from botao import Botao
from shapely.geometry import Polygon, Point
import random
import copy

class CSP():
    def __init__(self,dataset='fu',Base=None,Altura=None,Escala=None,render=False,plot=True, x=-200, y=200):
        
        if render:
            turtle.tracer(0)
        self.x = x
        self.y = y
        
        self.render = render
        self.plot = plot
        
        self.pecas = dataset
        self.instancias(Base,Altura,Escala)
        self.area = self.base * self.altura
        
        self.lista_original = ler_poligonos('/home/fsilvestre/Cutting_Stock_Problem/' + self.pecas +'.dat')
        self.nova_lista_original, self.nova_lista_completa_original = tratar_lista(self.lista_original, self.Escala)

        self.tamanho_lista = len(self.lista_original)
        
        self.lista = ler_poligonos('/home/fsilvestre/Cutting_Stock_Problem/' + self.pecas +'.dat')
        self.nova_lista, self.nova_lista_completa = tratar_lista(self.lista, self.Escala)
        self.pecas_posicionadas = []
        self.lista_removida = []
        
        self.cordenadas_area = ((self.x, self.y),((self.x + self.base), self.y),((self.x + self.base), (self.y - self.altura)),(self.x, (self.y - self.altura)))
        
        self._img = Image.new('L', (self.base, self.altura), 0)
        self._matriz = np.array(self._img)
        self._fig, self._ax = plt.subplots()
        
        if self.render:
            self.botao_remover = Botao(350, 300, 200, 150)
            self.turtle_fecho = turtle.Turtle()
            self.turtle_dados = turtle.Turtle()
            self.lista_turtle = []
            self.renderizar()
            
        self.area_ocupada = 0
        self.fecho = 0
        self.vertices_fecho = []

        self.fr = False

        self.peca_maior = max([calcular_area(peca) for peca in self.nova_lista])

        #print([round(calcular_area(peca)/self.area, 4) for peca in self.nova_lista])

            
    def acao(self,peca,x,y,grau_indice):
        if peca < len(self.lista):     
            pol_verificacao = self.nova_lista_completa[peca]
            pol_posicionar = self.nova_lista[peca]
            graus = [0, 90, 180, 270]
            grau = graus[grau_indice]

            x, y = transformar_coordenadas_inversa(x, y, self.x, self.x + self.base, self.y, self.y - self.altura, self.base, self.altura)
    
            pontos_verificacao = [(int(x) + int(rotate_point(cor[0], cor[1], grau)[0]), int(y) + int(rotate_point(cor[0], cor[1], grau)[1])) for cor in pol_verificacao]
            pontos_posicionar = [(int(x) + int(rotate_point(cor[0], cor[1], grau)[0]), int(y) + int(rotate_point(cor[0], cor[1], grau)[1])) for cor in pol_posicionar]
            
            if self.posicao_disponivel(pontos_verificacao):
                self.area_ocupada += calcular_area(pol_posicionar)
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
                self.sobreposicao = 0
                for pol in self.pecas_posicionadas:
                    self.sobreposicao+=porcentagem_sobreposicao(pontos_posicionar, pol)


                cordenadas_area  = [(self.x, self.y),((self.x + self.base), self.y),((self.x + self.base), (self.y - self.altura)),(self.x, (self.y - self.altura))]
                sobreposicao_fora = 1 - porcentagem_sobreposicao(pontos_posicionar, cordenadas_area)
                self.sobreposicao += sobreposicao_fora



        else:
            print("indice maior que a lista")
     

    def atualizar_dados(self):
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
        if self.render:
            self.lista_removida.append([self.lista[peca], self.nova_lista[peca], self.nova_lista_completa[peca]])
 
        self.lista.pop(peca)
        self.nova_lista.pop(peca)
        self.nova_lista_completa.pop(peca)
  
    def remover_da_area(self):
        peca = -1
        if self.pecas_posicionadas:
            self.pecas_posicionadas.pop(peca)
            
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
        lista_restaurada_1 = []
        for item in self.lista_original:
            if item in self.lista:
                lista_restaurada_1.append(item)
                
        lista_restaurada_2 = []
        for item in self.nova_lista_original:
            if item in self.nova_lista:
                lista_restaurada_2.append(item)
                
        lista_restaurada_3 = []
        for item in self.nova_lista_completa_original:
            if item in self.nova_lista_completa:
                lista_restaurada_3.append(item)
                
        self.lista = lista_restaurada_1
        self.nova_lista = lista_restaurada_2
        self.nova_lista_completa = lista_restaurada_3
          
    def posicao_disponivel(self,pontos_pol):

        pontos = pontos_no_poligono(pontos_pol)

        if any(self.CordenadaProibida(px, py) for px, py in pontos):
            return False
        
        
        if any(self.CordenadaProibida(px, py) for px, py in pontos_pol):
            return False

        if any(any(ponto_dentro_poligono(cor[0], cor[1], pontos_pol, False) for cor in pol) for pol in self.pecas_posicionadas):
            return False
        
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

    def calcular_recompensa(self, tamanho, n, acao):
        if tamanho > len(self.lista):
            if acao == 1:
                area = round(calcular_area(self.pecas_posicionadas[-1])/self.peca_maior, 3)
                self.peca_maior = max([calcular_area(peca) for peca in self.nova_lista])
                if self.maior_peca(self.pecas_posicionadas[-1]):
                    min_dist = 10000
                    for cor in self.pecas_posicionadas[-1]:
                        for ponto in self.cordenadas_area:
                            distancia = normalized_distance(cor, ponto, self.base, self.altura)
                            if distancia < min_dist:
                                min_dist = distancia

                    recompensa = math.log2(-min_dist + 1) + 1

                    self.dist = round(recompensa, 4)
                    if recompensa >= 0.95:
                        recompensa += 1
                        self.fr = True
                else:
                    recompensa = -1 * (1 - area)


            else:
                if self.fr:
                    area = round(calcular_area(self.pecas_posicionadas[-1])/self.peca_maior, 3)
                    if self.nova_lista:
                        self.peca_maior = max([calcular_area(peca) for peca in self.nova_lista])
                    if self.maior_peca(self.pecas_posicionadas[-1]):
                        recompensa = 1 + self.area_ocupada/self.fecho
                        
                    else:
                        recompensa = self.area_ocupada/self.fecho
                else:
                    recompensa = -1

                                
            
        else:
            if self.fr:
                recompensa = -1
            else:
                recompensa = -1

            #self.remover_lista(n)
    
        return round(recompensa,2)

    def maior_peca(self, peca):
        area = calcular_area(peca)
        for pol in self.nova_lista:
            if area < calcular_area(pol):
                return False
        return True
        
    def instancias(self,Base, Altura, escala):
        if self.pecas == 'fu':
            self.base = 114
            self.altura = 102
            self.Escala = 3
            
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
            self.base = 414
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

    def draw_click(self,x1,y1, grau_indice,n):
        graus = [0,90,180,270]
        grau = graus[grau_indice]
        x, y = x1,y1
        x_escala = 0
        y_escala = 0
        if len(self.lista) > 0:
            pol = self.lista[n]
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
        self.acao(self.n,x,y,self.g)
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
            #self.draw_click(x1,y1,self.g,self.n)
            
            if len(self.nova_lista)==0:
                self.loop = False
                turtle.clear()
                break






    
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
        
def polygon_features(lista_poligonos):
    features = []
    
    for coords in lista_poligonos:
        
        
        
        poly = Polygon(coords)
        
      
        num_lados = len(coords)
        
        
        area = poly.area
        
        
        perimetro = poly.length
        
        
        centroide = poly.centroid
        centroide_x, centroide_y = centroide.x, centroide.y
        
       
        minx, miny, maxx, maxy = poly.bounds
        width = maxx - minx
        height = maxy - miny
        
        
        features.append([
            area,
            num_lados, 
            area, 
            perimetro, 
            centroide_x, 
            centroide_y, 
            width, 
            height
        ])
    
    return features
        

def max_coordinates_count(lista_poligonos):
    max_count = 0
    for coords in lista_poligonos:
        max_count = max(max_count, len(coords))
    return max_count

def matriz_cordenadas(lista):
    lista_poligonos = copy.deepcopy(lista)
    max_len = max_coordinates_count(lista_poligonos)  # Encontra o número máximo de coordenadas
    for pol in lista_poligonos:
        if len(pol) < max_len:  # Verifica se o polígono tem menos coordenadas do que o máximo
            pol.extend([(-1, -1)] * (max_len - len(pol)))  # Preenche com (-1, -1) se necessário
    return lista_poligonos

def pontos_no_poligono(vertices, num_pontos=5):
    """
    Gera pontos dentro de um polígono incluindo o centroide.

    :param vertices: Lista de tuplas com as coordenadas dos vértices do polígono [(x1, y1), (x2, y2), ...]
    :param num_pontos: Número total de pontos a serem gerados dentro do polígono
    :return: Lista de tuplas com as coordenadas dos pontos gerados [(x1, y1), (x2, y2), ...]
    """
    
    # Criar o polígono a partir dos vértices
    poligono = Polygon(vertices)
    
    # Calcular o centroide do polígono
    centroide = poligono.centroid
    pontos = [(centroide.x, centroide.y)]  # Incluir o ponto central
    
    # Gerar mais pontos dentro do polígono
    min_x, min_y, max_x, max_y = poligono.bounds  # Limites do polígono

    while len(pontos) < num_pontos:
        # Gerar pontos aleatórios dentro do bounding box do polígono
        x_rand = random.uniform(min_x, max_x)
        y_rand = random.uniform(min_y, max_y)
        ponto_aleatorio = Point(x_rand, y_rand)

        # Verificar se o ponto está dentro do polígono
        if poligono.contains(ponto_aleatorio):
            pontos.append((x_rand, y_rand))

    return pontos

def normalized_distance(point1, point2, rect_width, rect_height):

    distance = np.linalg.norm(np.array(point1) - np.array(point2))
    max_distance = np.linalg.norm(np.array([rect_width, rect_height]))
    normalized_dist = distance / max_distance
    
    return normalized_dist

def porcentagem_sobreposicao(pontos_A, pontos_B):
    pontos_A.append(pontos_A[0])
    pontos_B.append(pontos_B[0])

    poligono_A = Polygon(pontos_A)
    poligono_B = Polygon(pontos_B)
    
    if not poligono_A.is_valid:
        poligono_A = poligono_A.buffer(0)
    if not poligono_B.is_valid:
        poligono_B = poligono_B.buffer(0)
    area_A = poligono_A.area

    interseccao = poligono_A.intersection(poligono_B)
    
    area_interseccao = interseccao.area
    porcentagem = round((area_interseccao / area_A),2)
    
    return porcentagem

ambieente = CSP(render=False, plot=False)