import sys
import os

sys.path.append(os.path.abspath("/home/fsilvestre"))



import numpy as np
from PIL import Image, ImageDraw
import matplotlib.pyplot as plt
import turtle
import math
from Cutting_Stock_Problem.Ambiente.Main.botao import Botao
from Cutting_Stock_Problem.Ambiente.Teste.nfp_teste import combinar_poligonos, triangulate_shapely,NoFitPolygon

from scipy.spatial import ConvexHull
import numpy as np
import cv2
import copy

def NFP(PecaA,grauA,PecaB,grauB):
    graus = [0,90,180,270]
    grauA = graus[grauA]
    grauB = graus[grauB]
    pontos_pol_A = [(int(rotate_point(cor[0], cor[1], grauA)[0]), int(rotate_point(cor[0], cor[1], grauA)[1])) for cor in PecaA]
    pontos_pol_B = [(int(rotate_point(cor[0], cor[1], grauB)[0]), int(rotate_point(cor[0], cor[1], grauB)[1])) for cor in PecaB]
    nfps_CB_CA = []
    convex_partsB = triangulate_shapely(pontos_pol_B)
    
    convex_partsA = triangulate_shapely(pontos_pol_A)

    nfps_convx = []
    for CB in convex_partsB:
        for convex in convex_partsA:
            nfps_convx.append(Polygon(NoFitPolygon(convex, CB)))


    nfp_final = list(combinar_poligonos(nfps_convx).exterior.coords)
    #print(nfp_final)

    return nfp_final

class CSP():
    def __init__(self,dataset='fu',Base=None,Altura=None,Escala=None,render=False,plot=True, x=-200, y=200, suavizar = True, pre_processar = None, margem = 5, ajuste = False):
        
        self.x = x
        self.y = y
        
        self.render = render
        self.plot = plot
        
        self.pecas = dataset
        self.rotacoes = [0, 1, 2, 3]

        if Base == None and Altura == None:
            self.instancias()
        else:
            self.base = Base
            self.altura = Altura
            self.Escala = 0.3
        self.area = self.base * self.altura

        
       
        iteracoes = 10  # Número de iterações para a abertura
        tamanho_kernel = 15  # Tamanho do kernel para suavização
        area_minima = 50  # Área mínima para considerar o polígono válido
        epsilon = 1  # Parâmetro de suavização para simplificação de vértices


        lista = []
        try:
            self.lista_original = ler_poligonos('/home/fsilvestre/Cutting_Stock_Problem/Datasets/' + self.pecas +'.dat')
        except FileNotFoundError:
            self.lista_original = ler_poligonos(self.pecas)
        
        if suavizar:
            lista = [suavizar_poligono(idx, self.lista_original, iteracoes, tamanho_kernel, area_minima, epsilon) for idx in range(len(self.lista_original))]
        else:
            lista = self.lista_original
        

        
        self.lista_original = copy.deepcopy(lista)
        self.nova_lista_original, self.nova_lista_completa_original = tratar_lista(self.lista_original, self.Escala)
        
        if ajuste:
            self.lista_original = [ajustar_poligono(pol) for pol in self.lista_original]
            self.nova_lista_original = [ajustar_poligono(pol) for pol in self.nova_lista_original]
            self.nova_lista_completa_original = [ajustar_poligono(pol) for pol in self.nova_lista_completa_original]

        self.max_pecas = len(self.lista_original)
        

        
        self.lista, self.nova_lista, self.nova_lista_completa = copy.deepcopy(self.lista_original), copy.deepcopy(self.nova_lista_original), copy.deepcopy(self.nova_lista_completa_original)
        
        if pre_processar is None:
            self.tabela_nfps = pre_processar_NFP(self.rotacoes, self.nova_lista, margem)

        else:
            self.tabela_nfps = pre_processar


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
            self.botao_reset = Botao(self.x + 100, self.y + 300, 200, 150)

            self.turtle_fecho = turtle.Turtle()
            self.turtle_dados = turtle.Turtle()
            self.lista_turtle = []
            self.renderizar()
            
        self.area_ocupada = 0
        self.fecho = 0
        self.vertices_fecho = []
            
    def acao(self,peca,x,y,grau_indice,flip = False,verificado=False):
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
                
                if self.posicao_disponivel(peca,grau_indice,x,y,pontos_verificacao):
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
     
    def nfp(self, peca, grau_indice):
        nfps = []
        for x2, y2, grau1, pol in self.indices_pecas_posicionadas:
        # Converta as peças para tuplas ao acessar a tabela
            chave = (tuple(self.nova_lista_original[pol]), grau1, 
                    tuple(self.nova_lista[peca]), grau_indice)
        
        # Acesse o NFP com a chave correta
            nfp1 = [(x1 + x2, y1 + y2) for x1, y1 in self.tabela_nfps[chave]]
            #print(x2,y2)
            # Adicione o polígono ao resultado
            nfps.append(Polygon(copy.deepcopy(nfp1)))

    # Combine os polígonos
        nfp = list(combinar_poligonos(nfps).exterior.coords)

        return nfp

    def resetar(self):
        for i in range(len(self.pecas_posicionadas)):
            self.remover_da_area()

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
      
    def posicao_disponivel(self,peca,grau,x,y, pontos_verificação):
        
        Dentro = True
        for x1,y1 in pontos_verificação:
            if not ponto_dentro_poligono(x1,y1,self.cordenadas_area, True):
                Dentro = False
        
        if len(self.pecas_posicionadas) == 0:
            return Dentro
        else:
            nfp = self.nfp(peca, grau)
            if ponto_dentro_poligono(x,y,nfp,True):
                return False
            else:
                return Dentro
  
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
        self._img = Image.new('L', (int(self.base), int(self.altura)), 0)
        for poligono in self.pecas_posicionadas:
            poligono_transformado = [transformar_coordenadas(x, y, self.x, (self.x + self.base), (self.y - self.altura), self.y, self.base, self.altura) for x, y in poligono]
            ImageDraw.Draw(self._img).polygon(poligono_transformado, outline=1, fill=1)
        self._matriz = np.array(self._img)
        
        if self.plot:
            self._plot()
        
    def instancias(self):
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
            self.base = 500
            self.altura = 500    
            self.Escala = 0.3

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
            if self.botao_reset.clicou(x,y):
                self.n = 0
                self.resetar()
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

            # x1,y1 = pyautogui.position()
            x1 = 963
            y1 = 533
             
            self.draw_click(x1,y1,self.g,self.n, False)
            
            if len(self.nova_lista)==0:
                self.loop = False
                turtle.clear()
                break


def offset_polygon(vertices, offset):
    """
    Cria um novo polígono com offset a partir de um polígono original usando Shapely
    
    Args:
        vertices: Lista de tuplas (x, y) representando os vértices do polígono
        offset: Valor do offset (positivo para expandir, negativo para contrair)
    
    Returns:
        Lista de tuplas (x, y) representando os vértices do novo polígono
    """
    if offset > 0:
        # Cria um polígono Shapely
        poly = Polygon(vertices)
        
        # Verifica se o polígono é válido
        if not poly.is_valid:
            return vertices
        
        # Aplica o buffer (offset)
        # join_style=1 (round) para suavizar cantos
        # mitre_limit controla quanto os cantos podem se estender
        buffered = poly.buffer(offset, join_style=1, mitre_limit=2.0)
        
        # Se o resultado for vazio ou inválido, retorna o original
        if buffered.is_empty or not buffered.is_valid:
            return vertices
        
        # Extrai os vértices do polígono resultante
        if buffered.geom_type == 'Polygon':
            # Pega apenas o exterior do polígono
            new_vertices = list(buffered.exterior.coords)[:-1]  # Remove o último ponto (duplicado)
        else:
            # Se o resultado for um MultiPolygon, pega o maior polígono
            largest = max(buffered.geoms, key=lambda x: x.area)
            new_vertices = list(largest.exterior.coords)[:-1]
        
        return new_vertices
    
    else:
        return vertices

def pre_processar_NFP(rotacoes, lista_pecas,offset):
    tabela_nfps = {}
    
    # Calcula o total de iterações
    total = len(lista_pecas) * len(rotacoes) * len(lista_pecas) * len(rotacoes)
    atual = 0
    
    for pecaA in lista_pecas:
        for grauA in rotacoes:
            for pecaB in lista_pecas:
                for grauB in rotacoes:
                    # Atualiza e mostra o progresso
                    atual += 1
                    porcentagem = (atual / total) * 100
                    print(f"\rPré-processando NFPs: {porcentagem:.1f}% concluído", end="")
                    
                   
                    chave = (tuple(pecaA), grauA, tuple(pecaB), grauB)
                    nfp = NFP(pecaA, grauA, pecaB, grauB)
                 
                    tabela_nfps[chave] = offset_polygon(nfp,offset)


    
   
    return tabela_nfps


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

def suavizar_poligono(index, poligonos, iteracoes=1, tamanho_kernel=5, area_minima=50, epsilon=5.0):
    """
    Suaviza um polígono e garante que a geometria resultante seja válida.
    
    Args:
        index: índice do polígono na lista
        poligonos: lista de polígonos
        iteracoes: número de iterações de suavização
        tamanho_kernel: tamanho do kernel para operações morfológicas
        area_minima: área mínima para filtrar polígonos
        epsilon: parâmetro de aproximação do polígono
    
    Returns:
        Lista de polígonos suavizados e validados
    """
    if index < 0 or index >= len(poligonos):
        print(f"Índice inválido. Escolha um índice entre 0 e {len(poligonos) - 1}.")
        return []

    # Selecionar o polígono
    poligono = poligonos[index]
    
    # Verificar se o polígono tem pontos suficientes
    if len(poligono) < 3:
        print("Polígono inválido: precisa ter pelo menos 3 pontos")
        return []

    # Determinar o tamanho da matriz de visualização com margem
    max_x = int(max(pt[0] for pt in poligono)) + 50
    max_y = int(max(pt[1] for pt in poligono)) + 50
    min_x = int(min(pt[0] for pt in poligono)) - 50
    min_y = int(min(pt[1] for pt in poligono)) - 50
    
    width = max_x - min_x
    height = max_y - min_y

    # Criar uma matriz binária
    canvas = np.zeros((height, width), dtype=np.uint8)
    
    # Ajustar coordenadas para o novo sistema
    adjusted_poly = np.array([(pt[0] - min_x, pt[1] - min_y) for pt in poligono], dtype=np.int32)
    
    # Preencher o polígono na matriz
    cv2.fillPoly(canvas, [adjusted_poly], 255)
    
    # Aplicar abertura morfológica para remover ruídos
    kernel = np.ones((tamanho_kernel, tamanho_kernel), np.uint8)
    for i in range(iteracoes):
        canvas = cv2.morphologyEx(canvas, cv2.MORPH_OPEN, kernel)
        canvas = cv2.morphologyEx(canvas, cv2.MORPH_CLOSE, kernel)

    # Encontrar novos contornos após a abertura
    contornos, _ = cv2.findContours(canvas, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # Filtrar e simplificar os contornos
    novos_poligonos = []
    for contorno in contornos:
        area = cv2.contourArea(contorno)
        if area >= area_minima:
            # Suavizar o contorno com aproximação de polígonos
            epsilon_value = epsilon * (cv2.arcLength(contorno, True) / 100)
            contorno_suavizado = cv2.approxPolyDP(contorno, epsilon_value, True)
            
            # Converter de volta para o sistema de coordenadas original
            pontos_ajustados = [(pt[0] + min_x, pt[1] + min_y) for pt in contorno_suavizado.reshape(-1, 2)]
            
            # Validar o polígono resultante
            if len(pontos_ajustados) >= 3:  # Verificar se tem pelo menos 3 pontos
                # Verificar se o polígono não tem auto-interseções
                poly = Polygon(pontos_ajustados)
                if poly.is_valid:
                    novos_poligonos.append(pontos_ajustados)
                else:
                    # Tentar corrigir o polígono
                    poly_corrigido = poly.buffer(0)
                    if poly_corrigido.is_valid:
                        coords = list(poly_corrigido.exterior.coords)[:-1]  # Remover o último ponto (duplicado)
                        novos_poligonos.append(coords)

   
        
    
    return novos_poligonos

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
from shapely.geometry import Polygon

def tratar_lista(lista_poligonos, Escala):
    def remove_vertices_repetidos(polygon):
        seen = set()
        unique_polygon = []
        for vertex in polygon:
            if vertex not in seen:
                unique_polygon.append(vertex)
                seen.add(vertex)
        return unique_polygon
    
    nova_lista = []
    nova_lista_completa = []
    
    for pol in lista_poligonos:
        novo_pol = []
        
        # Escalando os vértices
        for cor in pol:
            novo_pol.append((int(cor[0] * Escala), int(cor[1] * Escala)))
        
        # Removendo qualquer vértice repetido
        novo_pol = remove_vertices_repetidos(novo_pol)
        
        
        # Gerando pontos intermediários
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





