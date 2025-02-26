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
def draw_piece(i,j, grau,l):
    x, y = i,j
    pol = lista_poligonos[l]
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

    # Divide o conteúdo em blocos, separados por duas novas linhas
    blocos = conteudo.split('\n\n')

    poligonos = []

    for bloco in blocos:
        linhas = bloco.strip().split('\n')
        # Descarta a primeira linha que contém o número de lados
        pontos = [tuple(map(int, linha.split())) for linha in linhas[1:]]
        poligonos.append(pontos)

    return poligonos
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

def on_click2(x, y, grau,fecho_turtle, dados_turtle, botao, n,colunas,base,altura):
    global CordenadasProibidas
    global lista_colunas
    global fecho
    global AreaPreenchida
    global turtles
    global l
    global lista_area_fechos
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
        area_total = base*altura
        CordenadasProibidas.pop()
        fecho.pop()
        vertices_fecho = []
        if len(fecho)>0:
            vertices_fecho, area_fecho = calcular_fecho_e_area(fecho)
        else:
            area_fecho = 0
        dados_turtle.penup()
        dados_turtle.goto(-925,200)
        dados_turtle.pendown()
        if len(fecho)>0:
            função = (((((AreaPreenchida/(variacao_total_eixo_x(fecho)*altura))*0.5)**(0))*((AreaPreenchida/area_fecho)**(1))))        
            dados_turtle.write(f"Qualidade Preenchimento = {round(função, 2)}", font=("Arial",16, "normal"))
        else:
            dados_turtle.write(f"Qualidade Preenchimento = {0}", font=("Arial",16, "normal"))
        dados_turtle.penup()
        dados_turtle.goto(-925,350)
        dados_turtle.pendown()
        
        porcentagem_fecho = (area_fecho-AreaPreenchida)/area_total
        dados_turtle.write(f"Preenchimento total: {(round(((max(AreaPreenchida,sum(lista_area_fechos)+AreaPreenchida))/area_total)*100,2))}%", font=("Arial",16, "normal"))
        dados_turtle.penup()
        dados_turtle.goto(-925,300)
        dados_turtle.pendown()
        dados_turtle.write(f"Peças Disponiveis: {len(lista_poligonos)}", font=("Arial",16, "normal"))
        dados_turtle.penup()
        dados_turtle.goto(-925,250)
        dados_turtle.pendown()
        dados_turtle.write(f"Preenchimento coluna {len(lista_colunas)}: {AreaPreenchida}", font=("Arial",16, "normal"))
        dados_turtle.hideturtle()
        l = 0
        if len(fecho)>0:       
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
                    print("Error - CordenadaProibida")
                    posição_disponivel = False
                    
            
            for pol in CordenadasProibidas:
                if any(ponto_dentro_poligono(cor[0],cor[1],pontos_pol, False) for cor in pol):
                    #posição_disponivel = False
                    print(1)
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

                dados_turtle.write(f"Qualidade Preenchimento = {round(função, 2)}", font=("Arial",16, "normal"))
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
                dados_turtle.write(f"Preenchimento coluna {len(lista_colunas)}: {AreaPreenchida}", font=("Arial",16, "normal"))
                dados_turtle.hideturtle()
                
                desenhar_fecho(vertices_fecho, fecho_turtle)
                turtle.tracer(1)
                turtles.append([p,lista_poligonos[n],nova_lista[n],nova_lista_completa[n]])
                lista_poligonos.pop(n)
                nova_lista.pop(n)
                nova_lista_completa.pop(n)
                l = 0

                return função
def exact_fit(n, adicional, adicional2, m, draw, porcentagem, base, altura, real):
    print(n)
    nova_base = adicional2 + adicional
    Area.CordenadasProibidas[0][0] += adicional
    Area.CordenadasProibidas[1][0] = X + adicional + adicional2
    Area.CordenadasProibidas[2][0] = X + adicional2 + adicional
    Area.CordenadasProibidas[3][0] += adicional
    
    global lista_colunas, lista_area_fechos, fecho, AreaPreenchida
    tempo_anterior = 0
    nova_base = adicional2 + adicional + 1
    graus = [0, 90, 180, 270]

    def process_position(i, j, l, grau):
        global lista_colunas, lista_area_fechos, fecho, AreaPreenchida
        cordenada_disponivel = True
        poli = nova_lista[l]
        x_transformado, y_transformado = transformar_coordenadas_inversa(i, j, X, X+base, Y, Y-altura, base, altura)

        if draw:
            draw_piece(x_transformado, y_transformado, grau, l)
        
        pontos = [(x_transformado + int(rotate_point(cor[0], cor[1], grau)[0]), 
                   y_transformado + int(rotate_point(cor[0], cor[1], grau)[1])) 
                  for cor in poli]
        if calcular_largura(pontos) > adicional2:
            return -1
        
        
        if any(not ponto_dentro_poligono(cor[0],cor[1],Area.CordenadasProibidas,True)for cor in pontos):
            return -1
        for pol in CordenadasProibidas:
            if any(ponto_dentro_poligono(cor[0],cor[1],pontos, False) for cor in pol):
                return None

        
    
        #Area_Disponivel = (adicional2*altura)-AreaPreenchida
        #A_Poligono = calcular_area(pontos)
        #if A_Poligono > Area_Disponivel:
        #    return None
        
        if any(CordenadaProibida(cor[0], cor[1]) for cor in pontos):
            return None

        x3, y3 = calcular_centroide(pontos)
        if CordenadaProibida(x3, y3):
            return None


        
        pol = nova_lista_completa[l]
        pontos_pol = [(x_transformado + int(rotate_point(cor[0], cor[1], grau)[0]), 
                       y_transformado + int(rotate_point(cor[0], cor[1], grau)[1])) 
                      for cor in pol]
        
            
        if any(CordenadaProibida(cor[0], cor[1]) for cor in pontos_pol):
            return None

        fecho.append(pontos_pol)
        area_poligono = calcular_area(pontos_pol)
        AreaPreenchida += area_poligono
        vertices_fecho, area_fecho = calcular_fecho_e_area(fecho)
        area_desperdiçada =  (AreaPreenchida/area_fecho)
        AreaPreenchida -= area_poligono
        fecho.pop()
        
        return (l, area_desperdiçada, x_transformado, y_transformado, grau)

    while lista_poligonos:
        possiveis_posições = []
        for l in range(len(lista_poligonos)):
            for grau in graus:
                novo_inicio = adicional
                if grau == 90:
                    novo_inicio+=calcular_larguraY(lista_poligonos[l])*10
                if grau == 180:
                    novo_inicio+=calcular_largura(lista_poligonos[l])*10
                encontrou_posicao = False
                for i in range(int(novo_inicio), nova_base, 5):
                    if encontrou_posicao:
                        break
                    for j in range(0, altura, 5):
                        fim = time.time()
                        tempo_execucao = fim - inicio

                        if int(tempo_execucao) != int(tempo_anterior):
                            tempo(tempo_execucao)
                            tempo_anterior = int(tempo_execucao)
                            
                        pos = process_position(i, j, l, grau)
                        if pos == -1:
                            encontrou_posicao = True   
                        elif pos is not None:
                            possiveis_posições.append(pos)
                            encontrou_posicao = True
                            break
                        
                    if encontrou_posicao:
                        break
        if possiveis_posições:
            melhor_posicao = max(possiveis_posições, key=lambda x: x[1])
            on_click2(melhor_posicao[2], melhor_posicao[3], melhor_posicao[4], fecho_turtle, dados_turtle, Bot, melhor_posicao[0], n, base, altura)
            if len(lista_poligonos) == 0:
                return 0,tempo_execucao
        elif real:
            lista_colunas.append(fecho)
            lista_area_fechos.append(AreaPreenchida)
            AreaPreenchida = 0
            fecho = []
            if n == -1:
                return tempo_execucao
            else:
                return int(variacao_total_eixo_x(lista_colunas[n])),tempo_execucao
        else:
            return
def gerar_coluna(n,adicional,draw):
    global lista_poligonos,lista_colunas,nova_lista,nova_lista_completa, lista_area_fechos, fecho, AreaPreenchida
    lista_original = copy.deepcopy(lista_poligonos)
    nova_lista_original = copy.deepcopy(nova_lista)
    nova_lista_completa_original = copy.deepcopy(nova_lista_completa)
    tamanho = len(lista_poligonos)
    graus = [0]
    colunas = []
    for j in range(0, len(lista_poligonos), 1):
        acabou_coluna = False
        for grau in graus:
            for i in range(adicional, base, 5):
                #draw_piece(-200 + i, -180, grau, j)
                on_click2((X) + i, Y-altura, grau, fecho_turtle, dados_turtle, Bot, j, 3, base, altura)
                if tamanho > len(lista_poligonos):
                    adicional2 = variacao_total_eixo_x(fecho)
                    exact_fit(0, int(adicional), int(adicional2), 0, draw, 0.83, base, altura, False)
                    acabou_coluna = True
                    colunas.append([j, (AreaPreenchida/(altura*adicional2))*(len(CordenadasProibidas)), len(CordenadasProibidas)])
                    for indice in range(len(fecho)):
                        on_click2(400, 200, grau, fecho_turtle, dados_turtle, Bot, j, 0, base, altura)

                    lista_poligonos = copy.deepcopy(lista_original)  # Restaura a lista original usando deepcopy
                    nova_lista_completa = copy.deepcopy(nova_lista_completa_original)
                    nova_lista = copy.deepcopy(nova_lista_original)


                    Area.CordenadasProibidas[0][0]= X
                    Area.CordenadasProibidas[1][0]= X+base
                    Area.CordenadasProibidas[2][0]= X+base
                    Area.CordenadasProibidas[3][0]= X
                    break
                if acabou_coluna:
                    if len(colunas)==len(lista_poligonos):
                        acabou = True
                    
                    break
   
    melhor_coluna = max(colunas, key=lambda x: x[1])
    tamanho = len(lista_poligonos)
    for f in range(adicional, base, 5):

        on_click2((X)+f, Y-altura, grau, fecho_turtle, dados_turtle, Bot, melhor_coluna[0], 3, base, altura)
        if tamanho > len(lista_poligonos):
            break
    adicional2 = variacao_total_eixo_x(fecho)
    adicional1,tem = exact_fit(n, int(adicional), int(adicional2), 0, draw, 0.83, base, altura, True)
    Area.CordenadasProibidas[0][0]= X
    Area.CordenadasProibidas[1][0]= X + base
    Area.CordenadasProibidas[2][0]= X + base
    Area.CordenadasProibidas[3][0]= X
    return adicional1

def Greedy_Column_Generation(Base,Altura,x,y,peças,draw):
    global X,Y,tempoF,grau,CordenadasProibidas,lista_colunas,fecho,lista_area_fechos,AreaPreenchida,l,fecho_turtle,dados_turtle,turtles,lista_poligonos,t,inicio,tempo_anterior,fim,Bot,nova_lista,nova_lista_completa,base,altura
   
    X = x
    Y = y  
    base = Base
    altura = Altura
    
    tempoF = 0
    grau = 0
    CordenadasProibidas = []
    lista_colunas = []
    fecho = []
    lista_area_fechos = []
    AreaPreenchida = 0
    l = 0

    fecho_turtle = turtle.Turtle()  
    dados_turtle = turtle.Turtle()
    turtles = []

    instancias = peças + '.dat'
    lista_poligonos = ler_poligonos(instancias)
    lista_poligonos.pop(0)
    random.shuffle(lista_poligonos)

    t = turtle.Turtle()
    turtle.tracer(0)
    t.penup()
    t.goto(-900,375)
    t.pendown()
    turtle.tracer(1)
    inicio = time.time()
    tempo_anterior = 0

    fim = time.time()
    Bot = Botao((x+base+200), 300, 200, 150)

    turtle.tracer(0)
    
    area = Area(base,altura,6)
    area.area(x, y)
    print(Area.CordenadasProibidas)
    area.greedy(x,y)
    nova_lista, nova_lista_completa = tratar_lista(lista_poligonos)




    adicional = 0
    i=0
    acabou = False
    while not acabou: 
        res = gerar_coluna(i,int(adicional),draw)
        i+=1
        if res is not None:
            adicional+=res
            if (base - adicional)<(0.1*base):
                acabou = True
                break
        else:
            acabou = True
            break
        if len(lista_poligonos)==0:
            acabou = True
            break

    turtle.done()


Greedy_Column_Generation(340,380,-200,200,'fu',False)
            
    
    