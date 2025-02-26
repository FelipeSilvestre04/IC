
from CSP_GCG import CSP
from CSP_GCG import calcular_area, rotate_point, ponto_dentro_poligono
from NFP_GCG import combinar_poligonos,triangulate_shapely
from NFP_GCG import no_fit_polygon
import numpy as np
from shapely.geometry import Polygon, MultiPoint
import turtle
import numpy as np
import random
import turtle
import copy
import random
from scipy.spatial import ConvexHull

def GCG(tabela_nfps, peças='fu',draw=False,escala=None,plot=False,render=False, pesos = [2,2,1]):
    ambiente = CSP(peças,render=render,plot=plot)

    Colunas_Possivel = True
    inicio_coluna = 0
    Solucao = []
    Indices = []
    i = 0
    while Colunas_Possivel:
        i+=1
        print(f'Iteracao {i}')
        
        col,fecho,res = gerar_coluna(tabela_nfps,ambiente,inicio_coluna,pesos)

        

        if res is not None:
            inicio_coluna+=res
            for pol in col:
                
                    ambiente.acao(ambiente.lista.index(pol[-1]),pol[0],pol[1],pol[2],False,True)

            for peca in ambiente.pecas_posicionadas:
                Solucao.append(peca)

            for idx in ambiente.indices_pecas_posicionadas:
                Indices.append(idx)
            ambiente.pecas_posicionadas.clear()
            ambiente.indices_pecas_posicionadas.clear()



        else:
            Colunas_Possivel = False
            break
    ambiente.pecas_posicionadas = Solucao
    ambiente.indices_pecas_posicionadas = Indices
    ambiente.atualizar_dados()

    if ambiente.lista:
        melhor_peca = []
        while melhor_peca is not None:
            melhor_peca = varredura(tabela_nfps,ambiente,0,ambiente.base,False,False)
            if melhor_peca is not None:
                ambiente.acao(ambiente.lista.index(melhor_peca[-1]), melhor_peca[0], melhor_peca[1], melhor_peca[2],False, True)
            else:
                break

    return Solucao,ambiente.area_ocupada/ambiente.area

def gerar_coluna(tabela_nfps,ambiente,inicio_culuna,pesos):
    x, y = ambiente.cordenadas_area[3]
    grau = 0
    graus = [0,90,180,270]

    colunas_possiveis = []
    ind = indices_sem_copias(ambiente.lista)
    for i in ind:
        for g in range(2):        
            coluna = []
            peca = ambiente.lista[i]
            peca_rot = [(round(rotate_point(cor[0], cor[1], graus[g])[0],0), round(rotate_point(cor[0], cor[1], graus[g])[1],0)) for cor in peca]
            
            minx = min([x for x,y in peca_rot]) 
            miny = min([y for x,y in peca_rot]) 

            maxx = max([x for x,y in peca_rot])

            minx*=ambiente.Escala
            miny*=ambiente.Escala
            maxx*=ambiente.Escala

            peca1 = [x + inicio_culuna - minx,y - miny, g ,ambiente.lista[i]]
            
            tam = len(ambiente.lista)
            ambiente.acao(i, x + inicio_culuna - minx, y - miny, g, False)
            
            if tam > len(ambiente.lista):
                coluna.append(copy.deepcopy(peca1))
            else:
                continue
    

            largura_coluna = maxx - minx

            melhor_peca = []
            while melhor_peca is not None:
                melhor_peca = varredura(tabela_nfps,ambiente,inicio_culuna,largura_coluna)
                tamanh = len(ambiente.lista)
                

                if melhor_peca is None:
                    break
                else:
                    ambiente.acao(ambiente.lista.index(melhor_peca[-1]), melhor_peca[0], melhor_peca[1], melhor_peca[2],False, True)
                    coluna.append(copy.deepcopy(melhor_peca))

            if coluna:
                # Número de peças na coluna
                pecas_coluna = len(coluna)

                # Calcula a área preenchida das peças na coluna
                area_preenchida = calcular_area_preenchida([ambiente.pecas_posicionadas[-1 * (f + 1)] for f in range(len(coluna))])

                # Calcula a área do fecho convexo das peças na coluna
                area_fecho = calcular_fecho_area([ambiente.pecas_posicionadas[-1 * (f + 1)] for f in range(len(coluna))])

                # Calcula o fecho da coluna
                fecho_coluna = area_preenchida / area_fecho

                preenchimento_coluna = area_preenchida / (largura_coluna * ambiente.altura)

                colunas_possiveis.append((copy.deepcopy(coluna),(fecho_coluna**pesos[0]) * (preenchimento_coluna**pesos[1]) * (pecas_coluna**pesos[2]),largura_coluna))
                for h in range(len(coluna)):
                    ambiente.remover_da_area()


            coluna.clear()

    maior = 0
    if colunas_possiveis:
        for coluna1 in colunas_possiveis:
            if coluna1[1] > maior:
                maior = coluna1[1]
                melhor_coluna = coluna1

        print(f"Coluna Escolhida! custo: {round(melhor_coluna[1], 2)}")
        return melhor_coluna
    else:
        return None,None,None
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



def pre_processar_NFP(rotacoes, lista_pecas):
    tabela_nfps = {}
    i = 0
    for pecaA in lista_pecas:
        i+=1
        print(f"iteracao {i} de {len(lista_pecas)}")
        for grauA in rotacoes:
            for pecaB in lista_pecas:
                for grauB in rotacoes:
                    chave = (tuple(pecaA), grauA, tuple(pecaB), grauB)
                    #print(chave)
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
def varredura(tabela_nfps,ambiente,inicio_coluna,largura_coluna,draw = False, maior_peca = False):
    Cordenadas_coluna = copy.deepcopy(ambiente.cordenadas_area)
    Cordenadas_coluna[0][0] += inicio_coluna
    if ambiente.x + inicio_coluna + largura_coluna <= ambiente.cordenadas_area[1][0]:
        Cordenadas_coluna[1][0] = ambiente.x + inicio_coluna + largura_coluna
        Cordenadas_coluna[2][0] = ambiente.x + largura_coluna + inicio_coluna
    else:
        return None
    Cordenadas_coluna[3][0] += inicio_coluna

    graus = [0,90,180,270]
    possiveis_posicoes = []
    posicao_valida = False
    idxs = indices_sem_copias(ambiente.lista)
    while ambiente.lista:
        for peca in idxs:
            #print(f"peca{peca}")
            pol_posicionar = ambiente.nova_lista[peca]
            pecas_rot = []
            for grau1 in graus:                
                pecas_rot.append([[(round(rotate_point(cor[0], cor[1], grau1)[0], 1), round(rotate_point(cor[0], cor[1], grau1)[1],1)) for cor in pol_posicionar],grau1])
            
            rotacoes_validas = rotacoes_distintas(pecas_rot)
                    
            for grau in rotacoes_validas:
                pontos_pol = [(round(rotate_point(cor[0], cor[1], grau)[0]), round(rotate_point(cor[0], cor[1], grau)[1])) for cor in pol_posicionar]
                
                nfps = []
                for x2, y2, rot, pol in ambiente.indices_pecas_posicionadas:
                    # Converta as peças para tuplas ao acessar a tabela
                    chave = (tuple(ambiente.nova_lista_original[pol]), rot, 
                            tuple(pol_posicionar), graus.index(grau))
                    
                    # Acesse o NFP com a chave correta
                    nfp1 = [(x1 + x2, y1 + y2) for x1, y1 in tabela_nfps[chave]]
                    #print(x2,y2)
                    # Adicione o polígono ao resultado
                    nfps.append(Polygon(copy.deepcopy(nfp1)))

                # Combine os polígonos
                nfp = list(combinar_poligonos(nfps).exterior.coords)
                #print(grau)
                posicao_valida = False
                acabou = False

                novo_inicio = inicio_coluna
                nova_base = largura_coluna + inicio_coluna + 1
                if grau == 90:
                    novo_inicio+=calcular_larguraY(ambiente.nova_lista[peca])
                if grau == 180:
                    novo_inicio+=calcular_largura(ambiente.nova_lista[peca])
                
                if ambiente.Escala < 1:
                    Passo = 1
                else:
                    Passo = ambiente.Escala
                for i in range(int(novo_inicio),int(nova_base),Passo):
                    if posicao_valida:
                        break
                    if acabou:
                        break
                    for j in range(0,ambiente.altura,Passo):
                        x_t,y_t = transformar_coordenadas_inversa(i, j, ambiente.x, ambiente.x+ambiente.base, ambiente.y, ambiente.y-ambiente.altura, ambiente.base, ambiente.altura)
                        if draw:
                            for k in range(10):
                                ambiente.draw_click(x_t,-y_t,graus.index(grau),peca,False)

                        #print(i,j)
                        pontos = [(round(x_t + rotate_point(cor[0], cor[1], grau)[0]), round(y_t + rotate_point(cor[0], cor[1], grau)[1])) for cor in pol_posicionar]
                        posicao_invalida = False
                        acabou = False
                        for x1,y1 in pontos:
                            if y1 > ambiente.y:
                                posicao_invalida = True
                                break

                            if x1 > Cordenadas_coluna[1][0]:
                                acabou = True
                        if posicao_invalida:
                            break

                        if acabou:
                            break
                        posicao_valida = True
                        for x1,y1 in pontos:
                            if not ponto_dentro_poligono(x1,y1,Cordenadas_coluna, True) or not ponto_dentro_poligono(x1,y1,ambiente.cordenadas_area, True):
                                posicao_valida = False
                        if ponto_dentro_poligono(x_t,y_t,nfp, False):
                            posicao_valida = False
                        if posicao_valida:

                            #print(1)
                            pecas_posicionadas = copy.deepcopy(ambiente.pecas_posicionadas)
                            pecas_posicionadas.append(pontos)

                            area_ocupada = calcular_area_preenchida(pecas_posicionadas)
                            area_FC = calcular_fecho_area(pecas_posicionadas)
                            area_FR = area_fecho_retangular(pecas_posicionadas)

                            fecho_convexo = (1 - (area_ocupada/area_FR))
                            fecho_convexo_max = area_ocupada/area_FC
                            fecho_retangular_max = area_ocupada/area_FR
                            fechoR = (1 - (area_ocupada/area_FR))
                            nfp_meu = (fecho_convexo_max**10)*(area_FC/(ambiente.area))

                            pecas_posicionadas = []
                            #if inicio_coluna > 0:
                            #    for i in range(200):
                             #       ambiente.draw_click(x_t,-y_t,graus.index(grau),peca,False)
                              #      print(fecho_convexo_max)
                            if maior_peca:
                                possiveis_posicoes.append((peca,x_t,y_t,graus.index(grau),round(fecho_convexo_max * fecho_retangular_max,3),calcular_area(pol_posicionar),ambiente.lista[peca]))
                            else:
                                possiveis_posicoes.append((peca,x_t,y_t,graus.index(grau),round(fecho_convexo_max,3),ambiente.lista[peca]))
                            break


    

        if possiveis_posicoes:
            if maior_peca:
                maior = float('-inf')
                possiveis_pecas = []
                
                for posicao1 in possiveis_posicoes:
                    if posicao1[5] > maior:
                        possiveis_pecas = [copy.deepcopy(posicao1)]
                        maior = posicao1[5]
                    elif posicao1[5] == maior:
                        possiveis_pecas.append(copy.deepcopy(posicao1))
                
                # Encontrar a peça com a maior área entre as possíveis peças
                area = float('-inf')
                melhor_posicao = None
                for pc in possiveis_pecas:
                    area_atual = pc[4]
                    if area_atual > area:
                        area = area_atual
                        melhor_posicao = copy.deepcopy(pc)


                peca = copy.deepcopy([melhor_posicao[1],melhor_posicao[2],melhor_posicao[3],melhor_posicao[-1]])
                #print(peca)
                possiveis_posicoes.clear()
                return copy.deepcopy(peca)
            else:
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
                    area_atual = calcular_area(pc[-1])
                    if area_atual > area:
                        area = area_atual
                        melhor_posicao = copy.deepcopy(pc)


                peca = copy.deepcopy([melhor_posicao[1],melhor_posicao[2],melhor_posicao[3],melhor_posicao[-1]])
                #print(peca)
                possiveis_posicoes.clear()
                return copy.deepcopy(peca)
        else:
            return None
            

def NoFitPolygon(poligono_fixo, poligono_orbital):
    global inter
    polyA = Polygon(poligono_fixo)
    polyB = Polygon(poligono_orbital)

    return interpolar_pontos_poligono(no_fit_polygon(polyA, polyB), inter)


def best_posicao(possiveis_posicoes):
    return max(possiveis_posicoes, key=lambda p: p[4])

def transformar_coordenadas_inversa(x, y, x_min, x_max, y_min, y_max, width, height):
    x_scaled = (x / width) * (x_max - x_min) + x_min
    y_scaled = ((height - y) / height) * (y_max - y_min) + y_min

    return x_scaled, y_scaled

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
    
def calcular_area_preenchida(poligonos):
    area = 0
    for pol in poligonos:
        area += calcular_area(pol)
    return area

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

def calcular_largura(coordenadas):
    x_valores = [coord[0] for coord in coordenadas]
    largura = max(x_valores) - min(x_valores)
    return int(largura)
def calcular_larguraY(coordenadas):
    x_valores = [coord[1] for coord in coordenadas]
    largura = max(x_valores) - min(x_valores)
    return int(largura)
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

def append_to_txt(data, filename):
    with open(filename, "a") as file:
        file.write(f"Iteracao {data['iteration']}: {data['num_pieces']} pecas, {data['area_percent']}% area, {data['execution_time']:.2f} segundos, pecas:{data['pieces']}\n")

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


#rotacoes = [0, 1, 2]

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
filename = f"/home/fsilvestre/Cutting_Stock_Problem/resultados_GCG.txt"

rotacoes = [0, 1, 2,3]
instancias = ['fu','jackobs1','jackobs2','dighe1','dighe2','albano','mao','marques','dagli','shapes0','shapes1','shapes2']
pesos = [[2,2,1],[2,2,0],[1,1,1],[1,1,0],[0,1,0],[1,0,0]]
global inter
inter = 1
for peso in pesos:
    for instancia in instancias:
        env = CSP(instancia, plot=False,render=False)
        Stime = time.time()
        tabela_nfps = pre_processar_NFP(rotacoes, env.nova_lista)
        Etime = time.time()    
       
        with open(filename, "a") as file:
            file.write("---------------------------------------------------------------------------------------------\n")
            file.write(f"GCG ORIGINAL {instancia}, pre_computar nfp = {round(Etime - Stime,1)} pesos = {peso}, interpolacao - {inter}\n")

        tempos = []

        for i in range(5):
            start_time = time.time()  # Inicia a medição do tempo
            pecas,area = GCG(tabela_nfps, instancia,plot=False,render=False,pesos=peso)
            end_time = time.time()  # Termina a medição do tempo
            execution_time = end_time - start_time  # Calcula o tempo de execução
            print(f"Iteracao {i+1}: {len(pecas)} pecas, {round(area*100,2)}% , {execution_time:.2f} segundos")

            resultado = {
                "iteration": i + 1,
                "num_pieces": len(pecas),
                "area_percent": round(area*100, 2),
                "execution_time": execution_time,
                "pieces":pecas
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
    

    
        