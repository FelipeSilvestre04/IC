import turtle
from geradores import Area
from botao import Botao
import pyautogui
import math
import numpy as np
from scipy.spatial import ConvexHull
import time
import random


grau = 0
CordenadasProibidas = []
AreaPreenchida = 0
l = 0
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

def clicl2(x,y):
    click_meio()

def click_meio():
    global l
    l+=1
    if l > (len(lista_poligonos)-1):
        l = 0

def click(x,y):
    click_direito()
        
def click_direito():
    global grau
    grau +=30


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
    vertices = [vertices[-1]] + vertices

    for i in range(len(vertices) - 1):
        x1, y1 = vertices[i]
        x2, y2 = vertices[i + 1]
        # Calcula as diferenças
        dx = x2 - x1
        dy = y2 - y1
        # Calcula o MDC das diferenças
        passo_mdc = mdc(abs(dx), abs(dy))
        # Calcula os incrementos
        inc_x = dx // passo_mdc
        inc_y = dy // passo_mdc
        # Adiciona os pontos à lista
        for j in range(0, int(passo_mdc)):  # Começa em 1 para evitar duplicação dos vértices
            ponto = (x1 + j * inc_x, y1 + j * inc_y)
            pontos.append(ponto)

    return pontos




def rotate_point(x, y, angle):
    rad = math.radians(angle)
    x_new = x * math.cos(rad) - y * math.sin(rad)
    y_new = x * math.sin(rad) + y * math.cos(rad)
    return x_new, y_new

def draw_piece(i,j, grau):
    x, y = i,j
    pol = lista_poligonos[0]
    turtle.tracer(0)
    turtle.clear()
    turtle.begin_fill()
    turtle.fillcolor("green")
    turtle.penup()
    x_rot, y_rot = rotate_point(pol[0][0], pol[0][1], grau)
    turtle.goto(x+x_rot*10, y+y_rot*10)
    turtle.pendown()
    for cor in pol:

        x_rot, y_rot = rotate_point(cor[0], cor[1], grau)
        turtle.goto(x+x_rot*10, y+y_rot*10)
    turtle.goto(x+x_rot*10, y+y_rot*10)
    turtle.end_fill()
    turtle.hideturtle()
    turtle.tracer(1)

def calcular_area(vertices):
    """
    Calcula a área de um polígono definido pelos vértices.

    Args:
        vertices (list): Uma lista de tuplas (x, y) representando os vértices.

    Returns:
        float: A área do polígono.
    """
    num_vertices = len(vertices)
    area = 0

    for i in range(num_vertices):
        x1, y1 = vertices[i]
        x2, y2 = vertices[(i + 1) % num_vertices]  # Próximo vértice (considerando o último vértice conectado ao primeiro)

        area += (x1 * y2) - (x2 * y1)

    return abs(area) / 2

fecho_turtle = turtle.Turtle()  
dados_turtle = turtle.Turtle()
turtles = []

def on_click(x,y):
    on_click2(x,y,grau,fecho_turtle,dados_turtle,Bot,l)
    global CordenadasProibidas
    poligonos = []
    x = []
    y = []
    if CordenadasProibidas:
        for pol in CordenadasProibidas:
            for cor in pol:
                x.append(cor[0])
                y.append(cor[1])
            poligonos.append([x,y])
            x = []
            y = []
            
        for pol in poligonos:
            print(max(pol[0]),min(pol[0]))
            print(max(pol[1]),min(pol[1]))
            







def on_click2(x, y, grau,fecho_turtle, dados_turtle, botao, n):
    global CordenadasProibidas
    global AreaPreenchida
    global turtles
    global l
    # Certifique-se de que há um polígono para trabalhar
    if botao.clicou(x,y):
        fecho_turtle.clear()
        dados_turtle.clear()
        turtles[-1][0].clear()
        lista_poligonos.insert(n,turtles[-1][1])
        nova_lista.insert(n,turtles[-1][2])
        nova_lista_completa.insert(n,turtles[-1][3])

        turtles.pop()
        turtle.tracer(0)
        area_poligono = calcular_area(CordenadasProibidas[-1])
                
        AreaPreenchida -= area_poligono
        area_total = 340*380
        CordenadasProibidas.pop()
        vertices_fecho = []
        if len(CordenadasProibidas)>0:
            vertices_fecho, area_fecho = calcular_fecho_e_area(CordenadasProibidas)
        else:
            area_fecho = 0
        dados_turtle.penup()
        dados_turtle.goto(-925,200)
        dados_turtle.pendown()
        função = ((((area_fecho/area_total)**(0))*((AreaPreenchida/area_fecho)**(1))))
        dados_turtle.write(f"Qualidade Preenchimento = {round(função, 2)}", font=("Arial",16, "normal"))
        dados_turtle.penup()
        dados_turtle.goto(-925,350)
        dados_turtle.pendown()
        
        porcentagem_fecho = (area_fecho-AreaPreenchida)/area_total
        dados_turtle.write(f"Preenchimento total: {(round((AreaPreenchida/area_total)*100,2))}%", font=("Arial",16, "normal"))
        dados_turtle.penup()
        dados_turtle.goto(-925,300)
        dados_turtle.pendown()
        dados_turtle.write(f"Peças Disponiveis: {len(lista_poligonos)}", font=("Arial",16, "normal"))
        dados_turtle.hideturtle()
        l = 0
        if len(CordenadasProibidas)>0:       
            desenhar_fecho(vertices_fecho, fecho_turtle)
        turtle.tracer(1)
    else:
        if len(lista_poligonos) > n:
            posição_disponivel = True
            pol = nova_lista_completa[n]
            
            pontos_pol = []
            for cor in pol:
                x_rot, y_rot = rotate_point(cor[0], cor[1], grau)
                pontos_pol.append((int(x) + int(x_rot),int(y) + int(y_rot)))
                if CordenadaProibida(int(x) + int(x_rot),int(y) + int(y_rot)):
                    posição_disponivel = False
            
            if any(poligonos_se_sobrepoem(p,pontos_pol) for p in CordenadasProibidas):
                    print(1)
                    posição_disponivel = False
            if posição_disponivel:
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
 
                CordenadasProibidas.append(pontos_pol)
                area_poligono = calcular_area(pontos_pol)
                
                AreaPreenchida += area_poligono
                vertices_fecho = []
                vertices_fecho, area_fecho = calcular_fecho_e_area(CordenadasProibidas)
                #area_desperdiçada = area_desperdiçada_fecho(area_fecho)
                    #print(f"Area desperdicada do fecho convexo: {round(area_desperdiçada, 2)}")

                
                turtle.tracer(0)
                #SumQ_local = calcular_Qs_Local(Area.GreedysOcupados, Area.CordenadasGreedy, Area.AreaOcupada, matriz._matriz, base, altura, greedys, p)
                #punicao = calcular_punição_fecho(area_desperdiçada, len(Area.GreedysOcupados), p, 40000)


                area_total = 340*380  
                dados_turtle.penup()
                dados_turtle.goto(-925,200)
                dados_turtle.pendown()
                função = ((((area_fecho/area_total)**(0))*((AreaPreenchida/area_fecho)**(1))))
                dados_turtle.write(f"Qualidade Preenchimento = {round(função, 2)}", font=("Arial",16, "normal"))
                dados_turtle.penup()
                dados_turtle.goto(-925,350)
                dados_turtle.pendown()
                
                porcentagem_fecho = (area_fecho-AreaPreenchida)/area_total
                dados_turtle.write(f"Preenchimento total: {(round((AreaPreenchida/area_total)*100,2))}%", font=("Arial",16, "normal"))
                dados_turtle.penup()
                dados_turtle.goto(-925,300)
                dados_turtle.pendown()
                dados_turtle.write(f"Peças Disponiveis: {len(lista_poligonos)-1}", font=("Arial",16, "normal"))
                dados_turtle.hideturtle()
                
                desenhar_fecho(vertices_fecho, fecho_turtle)
                turtle.tracer(1)
                turtles.append([p,lista_poligonos[n],nova_lista[n],nova_lista_completa[n]])
                lista_poligonos.pop(n)
                nova_lista.pop(n)
                nova_lista_completa.pop(n)
                l = 0
                return função
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
    


def calcular_punição_fecho(area_desperdiçada, n,p,area_greedy):
    punicao = (area_desperdiçada/(area_greedy*n))

    p.penup()
    p.goto(-925,350)
    p.pendown()
    p.write(f"P = {abs(round(punicao, 2))}", font=("Arial",16, "normal"))
    return punicao



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

    # Divide o conteúdo em blocos, separados por duas novas linhas
    blocos = conteudo.split('\n\n')

    poligonos = []

    for bloco in blocos:
        linhas = bloco.strip().split('\n')
        # Descarta a primeira linha que contém o número de lados
        pontos = [tuple(map(int, linha.split())) for linha in linhas[1:]]
        poligonos.append(pontos)

    return poligonos


# Exemplo de uso:
lista_poligonos = ler_poligonos('fu.dat')
lista_poligonos.pop(0)
#random.shuffle(lista_poligonos)
def tratar_lista(lista_poligonos):
    nova_lista = []
    nova_lista_completa = []
    for pol in lista_poligonos:
        novo_pol = []
        for cor in pol:
            novo_pol.append((cor[0]*10,cor[1]*10))
        pol_completo = pontos_entre_vertices(novo_pol)
        nova_lista_completa.append(pol_completo)
        nova_lista.append(novo_pol)

    return nova_lista, nova_lista_completa

def transformar_coordenadas_inversa(x, y, x_min, x_max, y_min, y_max, width, height):
    x_scaled = (x / width) * (x_max - x_min) + x_min
    y_scaled = ((height - y) / height) * (y_max - y_min) + y_min

    return x_scaled, y_scaled


t = turtle.Turtle()
turtle.tracer(0)
t.penup()
t.goto(-900,375)
t.pendown()
turtle.tracer(1)
inicio = time.time()
tempo_anterior = 0

fim = time.time()
def tempo(time):
    # Limpe a escrita anterior antes de escrever o novo tempo
    t.clear()
    t.write(f"{int(time/60)} min {round(time%60,0)} s",font=("Arial", 16, "normal"))
    t.hideturtle()
    turtle.update()

def draw_piece2(loop, grau,n):
    x, y = pyautogui.position()
    if len(lista_poligonos) > 0:
        pol = lista_poligonos[n]
        turtle.tracer(0)
        turtle.clear()
        turtle.begin_fill()
        turtle.fillcolor("green")
        turtle.penup()
        x_rot, y_rot = rotate_point(pol[0][0], pol[0][1], grau)
        turtle.goto(x+x_rot*10-963, -y+y_rot*10+528)
        turtle.pendown()
        for cor in pol:

            x_rot, y_rot = rotate_point(cor[0], cor[1], grau)
            turtle.goto(x+x_rot*10-963, -y+y_rot*10+528)
        turtle.goto(x+x_rot*10-963, -y+y_rot*10+528)
        turtle.end_fill()
        turtle.hideturtle()
        turtle.tracer(1)
    else:
        loop = False
        return loop
# Agora 'lista_poligonos' é uma lista de listas de vértices (tuplas)
Bot = Botao(350, 300, 200, 150)

turtle.tracer(0)

area = Area(340,380,6)
area.area(-200, 200)
area.greedy(-200,200)
loop = True
nova_lista, nova_lista_completa = tratar_lista(lista_poligonos)

while loop:
    
    turtle.onscreenclick(click,3)
    turtle.onscreenclick(clicl2,2)
    turtle.onscreenclick(on_click)
    
    draw_piece2(loop,grau,l)
    if len(lista_poligonos)==0:
        loop = False
        turtle.clear()
        break

turtle.done()


            
    
    