from Cutting_Stock_Problem.Ambiente.Main.CSP_concavo import CSP
from Cutting_Stock_Problem.Ambiente.Main.CSP_concavo import calcular_area, rotate_point, ponto_dentro_poligono
from NFP_concavo import combinar_poligonos,triangulate_shapely
from NFP_concavo import no_fit_polygon
import numpy as np
from shapely.geometry import Polygon, MultiPoint
import turtle
import numpy as np
import random
import turtle
import copy
from numba import njit



def GRASP(mIter, x,q1,q2,q3,rotacoes,pecas_inisol, instancia,tabela_nfps):
    for i in range(mIter):
        ambiente = CSP(dataset=instancia, render=False, plot=False)
        I1, I2 = IniSOL(ambiente,rotacoes,pecas_inisol, False, tabela_nfps)
        if len(I1) == ambiente.max_pecas:
            return I1
        fecho = ambiente.area_ocupada/ambiente.fecho
        I1_novo = LocSEARCH(I1,I2,x,q1,q2,q3,fecho, 0,instancia,tabela_nfps)
        if len(I1_novo) == ambiente.max_pecas:
            return I1_novo
        
        if i == 0:
            I = copy.deepcopy(I1_novo)

        elif porcentagem_area(I1_novo,ambiente) > porcentagem_area(I,ambiente):
            I = copy.deepcopy(I1_novo)


    return I


def LocSEARCH(I1, I2, x_dado, q1, q2, q3,fecho, pq,instancia,tabela_nfps):
    melhorou_global = True
    
    while melhorou_global:
        ambiente = CSP(dataset=instancia,render=False, plot=False)        
        melhorou_global = False
        area_global = porcentagem_area(I1,ambiente)
        for idk in range(0,x_dado,1):
            tipo = random.choices(['permutar', 'trocar', 'adicionar'], weights=[q1, q2, q3])[0]
            
            I1_original = copy.deepcopy(I1)
            I2_original = copy.deepcopy(I2)
            ambiente = CSP(dataset=instancia,render=False, plot=False)
            if ambiente.render:
                ambiente.turtle_dados.clear()

            
            
            if tipo == 'permutar':
                I1_linha = trocar_elementos_aleatorios(I1,pq)
            elif tipo == 'trocar':
                I1_linha, I2_linha = trocar_elementos_entre_listas(I1, I2,pq)
            elif tipo == 'adicionar':
                I1_linha, I2_linha = adicionar_elemento(I1, I2,pq)
            
            x, y = ambiente.cordenadas_area[3]

            minx = min([x for x,y in I1_linha[0][-1]]) 
            miny = min([y for x,y in I1_linha[0][-1]]) 

            minx*=ambiente.Escala
            miny*=ambiente.Escala

            index = ambiente.lista.index(I1_linha[0][-1])
            ambiente.acao(index, x - minx, y - miny, 0, False)


            I1_novo = []
            I1_novo = [I1_linha[0]]
            
            for peca in I1_linha[1:]:
                if peca[-1] in ambiente.lista:
                    index = ambiente.lista.index(peca[-1])
                    peca_posicao = PackS(ambiente, index, rotacoes, False,tabela_nfps)
                    if peca_posicao is not False:
                        I1_novo.append(peca)

                        ambiente.acao(peca_posicao[0], peca_posicao[1], peca_posicao[2], peca_posicao[3],peca_posicao[4],True)
                        
            area_nova = porcentagem_area(I1_novo,ambiente)
            melhorou = area_nova > area_global
            if melhorou:
                I1 = copy.deepcopy(I1_novo)
                I2 = copy.deepcopy(ambiente.lista)
                area_global = area_nova
                fecho = ambiente.area_ocupada/ambiente.fecho
                melhorou_global = True
                I1_novo.clear()
                
            elif (abs(round(area_nova,3) - round(area_global,3)) < 0.01*round(area_global,3)) and ambiente.area_ocupada/ambiente.fecho > fecho:
                I1 = copy.deepcopy(I1_novo)
                I2 = copy.deepcopy(ambiente.lista)
                area_global = area_nova
                fecho = ambiente.area_ocupada/ambiente.fecho
                melhorou = True
                melhorou_global = True
                I1_novo.clear()                
            else:
                I1 = I1_original
                I2 = I2_original
                I1_novo.clear()
            
            print(f"{tipo[0]}, {'M' if melhorou else 'NM'} {area_global:.3f} - {len(I1)} - {len(I2)}, {area_nova:.3f} - {ambiente.max_pecas - len(ambiente.lista)} - {len(ambiente.lista)}")
            if len(I1) == ambiente.max_pecas:
                if False:
                    turtle.clear()
                    ambiente = CSP(dataset=instancia,render=True, plot=False)
                    x, y = ambiente.cordenadas_area[3]
                    index = ambiente.lista.index(I1[0][-1])
                    ambiente.acao(index, x, y, 0, False)
                    for peca in I1[1:]:
                        if peca[-1] in ambiente.lista:
                            index = ambiente.lista.index(peca[-1])
                            peca_posicao = PackS(ambiente, index, rotacoes, True)
                            if peca_posicao is not False:
                                ambiente.acao(peca_posicao[0], peca_posicao[1], peca_posicao[2], peca_posicao[3],peca_posicao[4],True)
                    turtle.done

                print(f"{len(ambiente.pecas_posicionadas)}, {len(I1)}")
                return I1
                
        
    print(f"Iteração concluída. Tamanho de I1: {len(I1)}, Área: {area_global:.2f}")
    return I1

import random
from scipy.spatial import ConvexHull
def calcular_fecho_area(poligonos):
    if poligonos:
        # Concatenar todos os vértices de todos os polígonos
        vertices = np.concatenate(poligonos, axis=0)
        
        # Calcular o fecho convexo
        fecho = ConvexHull(vertices)
        
        # Obter os vértices do fecho convexo
        vertices_fecho = vertices[fecho.vertices]
        
        # Calcular a área do fecho convexo
        area = fecho.volume
        
        return area
    else:
        return 0
def trocar_elementos_entre_listas(lista1, lista2,pq):
    I1 = lista1.copy()
    I2 = lista2.copy()
    if not I1 or not I2:
        return I1, I2  # Se alguma lista estiver vazia, retorna as mesmas listas
    
    # Seleciona um índice aleatório em I1 e I2
    indice = round(pq*len(lista1),1)
    idx1 = random.randint(int(indice), len(I1) - 1)
    idx2 = random.randint(0, len(I2) - 1)
    
    # Troca os elementos entre as listas
    I1[idx1][-1], I2[idx2] = I2[idx2], I1[idx1][-1]
    
    return I1.copy(), I2.copy()

def trocar_elementos_aleatorios(I1,pq):
    lista = I1.copy()
    indice = round(pq * len(lista),1)
    idx1, idx2 = random.sample(range(int(indice),len(lista)-1), 2)
    
    # Troca os elementos nas posições idx1 e idx2
    lista[idx1], lista[idx2] = lista[idx2], lista[idx1]
    
    return lista.copy()

def porcentagem_area(I1,ambiente):
    area = 0
    for peca in I1:
        area+=calcular_area(peca[-1])
    area = round(area/(ambiente.area/ambiente.Escala**2), 3)
    return area    
    

def porcentagem_fecho(I1):
    lista = []
    for peca in I1:
        lista.append(peca[-1])
    
    return calcular_fecho_area(lista)

def indices_sem_copias(lista):
    vistos = set()  # Conjunto para rastrear objetos já encontrados
    indices = []    # Lista para armazenar os índices únicos

    for i, obj in enumerate(lista):
        # Converte a lista de vértices em uma tupla para torná-la imutável
        obj_tupla = tuple(tuple(v) for v in obj)

        if obj_tupla not in vistos:
            indices.append(i)  # Adiciona o índice à lista
            vistos.add(obj_tupla)  # Marca o objeto como visto

    return indices

def IniSOL(ambiente, rotacoes, pq,draw, tabela_nfps):
    I1 = []  # Peças posicionadas
    I2 = []  # Peças não posicionadas
    # Inicializa a primeira peça na coordenada inicial (x, y) do ambiente
    grau = 0
    x, y = ambiente.cordenadas_area[3 - grau]


    index = random.randint(0,len(ambiente.lista) -1)
    minx = min([x for x,y in ambiente.lista[index]]) 
    miny = min([y for x,y in ambiente.lista[index]]) 

    minx*=ambiente.Escala
    miny*=ambiente.Escala

    peca1 = [index,x - minx,y - miny,grau,0,ambiente.lista[index]]      
    ambiente.acao(index, x - minx, y - miny, grau, False)

    
    for i in range(pq):
        ix = random.randint(0, len(ambiente.lista) - 1)
        I2.append(ambiente.lista[ix])
        ambiente.remover_lista(ix)
    # I contém os índices das peças disponíveis no ambiente
    I = [i for i in range(len(ambiente.lista))]
 
    I1.append(peca1)
    acabou = False
    while not acabou:
        Possiveis_posicoes = []
        
        indexs = indices_sem_copias(ambiente.lista)
        for peca in indexs:
            peca_posicao = PackS(ambiente, peca, rotacoes,draw,tabela_nfps)
            if peca_posicao is not False:
                Possiveis_posicoes.append(peca_posicao)
        
        # Se houver posições possíveis para a peça, escolha a melhor
        if Possiveis_posicoes:
            print(1)
            max_valor = best_posicao(Possiveis_posicoes)[5]
            melhores_posicoes = [p for p in Possiveis_posicoes if p[5] == max_valor]

            melhor_posicao = max(melhores_posicoes, key=lambda p: calcular_area(ambiente.lista[p[0]]))
            I1.append(melhor_posicao)
            if draw:
                for i in range(1000):
                    ambiente.draw_click(melhor_posicao[1], -1*melhor_posicao[2], melhor_posicao[3],melhor_posicao[0], melhor_posicao[4])
            ambiente.acao(melhor_posicao[0], melhor_posicao[1], melhor_posicao[2], melhor_posicao[3], melhor_posicao[4], True)
        else:
            if ambiente.lista:
                for peca in ambiente.lista:
                    I2.append(peca)
            acabou = True

    print("Peças posicionadas:", len(I1))
    print("Peças restantes:", len(I2))
    
    return I1, I2  # Retorna a lista de peças posicionadas e as restantes

        
        
import math
        
def normalizar_poligono(vertices):
    # Encontra o menor valor de x e y para normalizar o polígono
    min_x = min(p[0] for p in vertices)
    min_y = min(p[1] for p in vertices)
    
    # Ajusta as coordenadas para que o menor ponto seja (0,0)
    return [(x - min_x, y - min_y) for x, y in vertices]

def rotacoes_distintas(rotacoes):
    indices_distintas = [0]  # Começamos com a rotação inicial como referência
    forma_original_normalizada = Polygon(normalizar_poligono(rotacoes[0][0]))

    for rotacao in rotacoes[1:]:
        forma_rotacionada_normalizada = Polygon(normalizar_poligono(rotacao[0]))
        
        # Verifica se a forma rotacionada é diferente da forma original
        if not forma_rotacionada_normalizada.equals(forma_original_normalizada):
            indices_distintas.append(rotacao[1])

    return indices_distintas
     
def PackS(ambiente, peca, rotacoes,draw,tabela_nfps):
    nfp = []
    lista = []
    melhores_posicoes = []

    

    pol_posicionar = ambiente.nova_lista[peca]
    graus = [0,90,180,270]
    pecas_rot = []
    for i in rotacoes:
        grau = graus[i] 
        pecas_rot.append([[(round(rotate_point(cor[0], cor[1], grau)[0], 1), round(rotate_point(cor[0], cor[1], grau)[1],1)) for cor in pol_posicionar],i])
    
    rotacoes_validas = rotacoes_distintas(pecas_rot)
    for grau_indice in rotacoes_validas:
        #print(peca, len(ambiente.lista))
        melhor_posicao = ExSEARCH(ambiente, peca, grau_indice,draw,tabela_nfps)
        if melhor_posicao is not None:
            melhores_posicoes.append(melhor_posicao)
    
    if melhores_posicoes:
        melhor = best_posicao(melhores_posicoes) 
        #print(melhor)
        return melhor
    else:
        return False




def pre_processar_NFP(rotacoes, lista_pecas):
    tabela_nfps = {}
    i = 0
    for pecaA in lista_pecas:
        i+=1
        print(f"iteracao {i}")
        for grauA in rotacoes:
            for pecaB in lista_pecas:
                for grauB in rotacoes:
                    
                    # Converta as peças para tuplas (se forem listas)
                    chave = (tuple(pecaA), grauA, tuple(pecaB), grauB)
                    tabela_nfps[chave] = NFP(pecaA, grauA, pecaB, grauB)
    
    return tabela_nfps


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
def ExSEARCH(ambiente, peca, grau_indice, draw,tabela_nfps):
    global fecho
    possiveis_posicoes = []
    nfp = []
    flip = False
    pol_posicionar = ambiente.nova_lista[peca]
    graus = [0,90,180,270]
    grau = graus[grau_indice]
    pontos_pol = [(int(rotate_point(cor[0], cor[1], grau)[0]), int(rotate_point(cor[0], cor[1], grau)[1])) for cor in pol_posicionar]
    lista = []

    nfps = []
    for x2, y2, grau1, pol in ambiente.indices_pecas_posicionadas:
        # Converta as peças para tuplas ao acessar a tabela
        chave = (tuple(ambiente.nova_lista_original[pol]), grau1, 
                tuple(ambiente.nova_lista[peca]), grau_indice)
        
        # Acesse o NFP com a chave correta
        nfp1 = [(x1 + x2, y1 + y2) for x1, y1 in tabela_nfps[chave]]
        #print(x2,y2)
        # Adicione o polígono ao resultado
        nfps.append(Polygon(copy.deepcopy(nfp1)))

    # Combine os polígonos
    nfp = list(combinar_poligonos(nfps).exterior.coords)

    for x,y in nfp:
                pontos = [(int(x + rotate_point(cor[0], cor[1], grau)[0]), int(y + rotate_point(cor[0], cor[1], grau)[1])) for cor in pol_posicionar]       
        #for flip in [True, False]:


                
                posicao_valida = True
                for x1,y1 in pontos:
                    if not ponto_dentro_poligono(x1,y1,ambiente.cordenadas_area, True):
                        posicao_valida = False
                if posicao_valida:
                    if draw:
                        for i in range(10):
                            ambiente.draw_click(x,-y,grau_indice,peca,flip)
                    pecas_posicionadas = copy.deepcopy(ambiente.pecas_posicionadas)
                    pecas_posicionadas.append(pontos)

                    area_ocupada = calcular_area_preenchida(pecas_posicionadas)
                    area_FC = calcular_fecho_area(pecas_posicionadas)
                    area_FR = area_fecho_retangular(pecas_posicionadas)

                    #fecho_convexo = (1 - (area_ocupada/area_FR))
                    fecho_convexo_max = area_ocupada/area_FC
                    fecho_retangular_max = area_ocupada/area_FR
                    #fechoR = (1 - (area_ocupada/area_FR))
                    #nfp_meu = (fecho_convexo_max**10)*(area_FC/(ambiente.area))

                    pecas_posicionadas = []
                    if fecho == "convexo":
                        possiveis_posicoes.append([peca,x,y,grau_indice,flip,round(fecho_convexo_max,2),round(fecho_retangular_max,2),ambiente.lista[peca]])
                    elif fecho == "retangular":
                        possiveis_posicoes.append([peca,x,y,grau_indice,flip,round(fecho_retangular_max,2),round(fecho_convexo_max,2),ambiente.lista[peca]])

    if possiveis_posicoes:
        maior = float('-inf')
        possiveis_pecas = []
                
        for posicao1 in possiveis_posicoes:
            if posicao1[4] > maior:
                possiveis_pecas = [copy.deepcopy(posicao1)]
                maior = posicao1[4]
            elif posicao1[4] == maior:
                possiveis_pecas.append(copy.deepcopy(posicao1))
                
                # Encontrar a peça com a maior área entre as possíveis peças
        area = float('-inf')
        melhor_posicao = None
        for pc in possiveis_pecas:
            area_atual = pc[5]
            if area_atual > area:
                area = area_atual
                melhor_posicao = copy.deepcopy(pc)
        
        return melhor_posicao
    else:
        return None


def best_posicao(possiveis_posicoes):
    return max(possiveis_posicoes, key=lambda p: p[5])

def calcular_area_preenchida(poligonos):
    area = 0
    for pol in poligonos:
        area += calcular_area(pol)
    return area

def NoFitPolygon(poligono_fixo, poligono_orbital):
    global inter
    polyA = Polygon(poligono_fixo)
    polyB = Polygon(poligono_orbital)

    return interpolar_pontos_poligono(no_fit_polygon(polyA, polyB), inter)

import random

def adicionar_elemento(I1, I2,pq):
    lista1 = copy.deepcopy(I1)
    lista2 = copy.deepcopy(I2)
    if not lista2:
        print("A lista2 está vazia. Nenhum item pode ser adicionado.")
        return lista1, lista2
    
    # Escolhe um elemento aleatório de lista2
    elemento = random.choice(lista2)
    lista2.remove(elemento)
    
    # Escolhe uma posição aleatória em lista1 para inserir o elemento
    indice = round(pq * len(lista1), 1)
    posicao = random.randint(int(indice), len(lista1) - 1)  
    
    lista1.insert(posicao, [0,0,0,0,0,0,elemento])
    

    
    return lista1.copy(), lista2.copy()
    
def area_fecho_retangular(poligonos):
    # Inicializar valores mínimos e máximos
    min_x = float('inf')
    max_x = float('-inf')
    min_y = float('inf')
    max_y = float('-inf')
    
    # Iterar sobre todos os vértices de todos os polígonos
    for poligono in poligonos:
        for (x, y) in poligono:
            # Atualizar os valores mínimos e máximos de x e y
            if x < min_x:
                min_x = x
            if x > max_x:
                max_x = x
            if y < min_y:
                min_y = y
            if y > max_y:
                max_y = y
    
    # Calcular a largura e altura do fecho retangular
    largura = max_x - min_x
    altura = max_y - min_y
    
    # Calcular a área do fecho retangular
    area = largura * altura
    
    return area

def escolher_percentual(lista, percentual):
    num_itens = int(len(lista) * percentual / 100)
    return random.sample(lista, num_itens)

def interpolar_pontos_poligono(vertices, pontos_entre_vertices=1):
    """
    Interpola pontos adicionais ao longo das arestas de um polígono.
    
    Args:
        vertices: Lista de vértices do polígono no formato [[x1,y1], [x2,y2], ...]
        pontos_entre_vertices: Número de pontos a serem adicionados entre cada par de vértices
    
    Returns:
        Lista com todos os pontos (originais + interpolados) que formam o polígono
    """
    def interpolar_pontos(p1, p2, num_pontos):
        """
        Interpola num_pontos entre p1 e p2.
        """
        x1, y1 = p1
        x2, y2 = p2
        pontos = []
        
        # Adiciona o ponto inicial
        pontos.append([float(x1), float(y1)])
        
        # Se não precisa adicionar pontos, retorna só o ponto inicial
        if num_pontos <= 0:
            return pontos
            
        # Adiciona os pontos intermediários
        for i in range(1, num_pontos + 1):
            t = i / (num_pontos + 1)
            x = x1 + t * (x2 - x1)
            y = y1 + t * (y2 - y1)
            pontos.append([float(x), float(y)])
            
        return pontos
    
    # Verifica se o polígono tem vértices suficientes
    if not vertices or len(vertices) < 3:
        return vertices

    pontos_finais = []
    num_vertices = len(vertices)
    
    # Para cada aresta do polígono
    for i in range(num_vertices):
        p1 = vertices[i]
        p2 = vertices[(i + 1) % num_vertices]
        
        # Interpola pontos entre p1 e p2
        pontos_aresta = interpolar_pontos(p1, p2, pontos_entre_vertices)
        pontos_finais.extend(pontos_aresta)
    
    return pontos_finais




import time
import turtle

def append_to_txt(data, filename):
    with open(filename, "a") as file:
        file.write(f"Iteracao {data['iteration']}: {data['num_pieces']} pecas, {data['area_percent']*100}% , {data['execution_time']:.2f} segundos\n")

rotacoes = [0, 1, 2]

vizinhanca = 30
q1 = 50
q2 = 20
q3 = 30
pecas_inisol = 0
max_iter = 80
import time
import math
import turtle

def append_to_txt(data, filename):
    with open(filename, "a") as file:
        file.write(f"Iteracao {data['iteration']}: {data['num_pieces']} pecas, {data['area_percent']}% area, {data['execution_time']:.2f} segundos, pecas{data['pieces']}\n")

def append_summary_to_txt(media, desvio_padrao, filename):
    with open(filename, "a") as file:
        file.write(f"\nTempo medio: {media:.2f} segundos\n")
        file.write(f"Desvio padrao: {desvio_padrao:.2f} segundos\n")
        file.write("---------------------------------------------------------------------------------------------\n")

def calcular_tempo_medio_desvio(tempos):
    # Calculando o tempo médio
    tempo_medio = sum(tempos) / len(tempos)
        
    # Calculando o desvio padrão
    soma_dos_quadrados = sum((x - tempo_medio) ** 2 for x in tempos)
    desvio_padrao = math.sqrt(soma_dos_quadrados / len(tempos))
    
    return tempo_medio, desvio_padrao


rotacoes = [0, 1, 2]

import time
import math
import datetime


def append_summary_to_txt(media, desvio_padrao, filename):
    with open(filename, "a") as file:
        file.write(f"\nTempo médio: {media:.2f} segundos\n")
        file.write(f"Desvio padrão: {desvio_padrao:.2f} segundos\n")
        file.write("---------------------------------------------------------------------------------------------\n")

def calcular_tempo_medio_desvio(tempos):
    # Calculando o tempo médio
    tempo_medio = sum(tempos) / len(tempos)
        
    # Calculando o desvio padrão
    soma_dos_quadrados = sum((x - tempo_medio) ** 2 for x in tempos)
    desvio_padrao = math.sqrt(soma_dos_quadrados / len(tempos))
    
    return tempo_medio, desvio_padrao

# Gerar um nome de arquivo com data e hora
current_time = datetime.datetime.now().strftime("%d-%m-%H_%M")
filename = f"/home/fsilvestre/Cutting_Stock_Problem/resultados_GRASP_original_{current_time}.txt"

#rotacoes = [0, 1, 2]
instancias = ['fu','dighe1','dighe2','jackobs1','fu','shapes0','shapes1','shapes2','albano',]
ins = ['dagli','jackobs2','jackobs1','fu','dighe1','dighe2','mao','albano','marques']
vizinhancas = [50,30,10]
pesos = [[30,20,50]]#[70,20,10]]
interpolacoes = [0,1]
from datetime import datetime
global fecho
# Obter a data e hora atuais formatadas
data_hora = datetime.now().strftime("%d/%m/%Y %H:%M")
#inisols = [0,1,2]
with open(filename, "a") as file:
    file.write("---------------------------------------------------------------------------------------------\n")
    file.write(f"Data e Hora: {data_hora}\n")
fechos = ["convexo", "retangular"]
#fecho = "convexo"

global inter
inter = 0
for instancia in ins:
    if instancia == "fu":
        rotacoes = [0, 1, 2]  # 0, 90, 180
    elif instancia == "jackobs1":
        rotacoes = [0, 1, 2]  # 0, 90, 180
    elif instancia == "jackobs2":
        rotacoes = [0, 1, 2]  # 0, 90, 180
    elif instancia == "shapes0":
        rotacoes = [0]  # Apenas 0
    elif instancia == "shapes1":
        rotacoes = [0, 2]  # 0, 180
    elif instancia == "shapes2":
        rotacoes = [0, 2]  # 0, 180
    elif instancia == "dighe1":
        rotacoes = [0]  # Apenas 0
    elif instancia == "dighe2":
        rotacoes = [0]  # Apenas 0
    elif instancia == "albano":
        rotacoes = [0, 2]  # 0, 180
    elif instancia == "dagli":
        rotacoes = [0]  # Apenas 0
    elif instancia == "mao":
        rotacoes = [0, 1, 2]  # 0, 90, 180
    elif instancia == "marques":
        rotacoes = [0, 1, 2]  # 0, 90, 180
    
    #fecho = 'convexo'
    for intp in interpolacoes:
                        inter = intp
        #for vin in vizinhancas:
                        for fecho1 in fechos:
                            fecho = fecho1
                
                
                #for p in inisols:    
                #for inter in interpolacoes:
                    #interpo = inter
                    #for Q1,Q2,Q3 in pesos:
                        #for vin in vizinhancas:   
                            env = CSP(instancia,render=False,plot=False)
                            Stime = time.time()
                            tabela_nfps = pre_processar_NFP(rotacoes, env.nova_lista)
                            Etime = time.time()
                            with open(filename, "a") as file:
                                file.write("---------------------------------------------------------------------------------------------\n")
                                file.write(f"GRASP ORIGINAL {instancia}: pre processar nfp {round(Etime-Stime,1)}, iteracoes maxima {max_iter}, {vizinhanca} vizinhancas, q1,q2,q2 = {q1},{q2},{q3}, interpolacao = {inter}, fecho = {fecho}, pecas inisol = {env.max_pecas - int(round(len(env.lista)*0.0,0))}/{env.max_pecas}, rotacoes{rotacoes}\n")

                            tempos = []

                            for i in range(5):
                                
                                start_time = time.time()  # Inicia a medição do 
                                pecas = GRASP(max_iter, vizinhanca, q1, q2, q3, rotacoes,int(round(len(env.lista)*0.0,0)),instancia, tabela_nfps)
                                end_time = time.time()  # Termina a medição do tempo
                                execution_time = end_time - start_time  # Calcula o tempo de execução
                                #env = CSP(instancia,render=False,plot=False)
                                print(f"Iteracao {i+1}: {len(pecas)} pecas, {porcentagem_area(pecas,env)} , {execution_time:.2f} segundos")

                                resultado = {
                                    "iteration": i + 1,
                                    "num_pieces": len(pecas),
                                    "area_percent": porcentagem_area(pecas,env) * 100,
                                    "execution_time": execution_time,
                                    "pieces": [pol[-1] for pol in pecas]
                                }
                                print(resultado)
                                # Armazena o tempo de execução
                                tempos.append(execution_time)
                                
                                # Salvar os resultados em um arquivo .txt a cada iteração
                                append_to_txt(resultado, filename)

                            # Calcula a média e o desvio padrão dos tempos
                            media, desvio_padrao = calcular_tempo_medio_desvio(tempos)
                            
                            # Salvar a média e o desvio padrão no arquivo
                            append_summary_to_txt(media, desvio_padrao, filename)

turtle.done()


    
        