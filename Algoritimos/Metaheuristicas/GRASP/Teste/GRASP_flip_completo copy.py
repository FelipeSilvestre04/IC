from Corte_e_Empacomento_PC_flip import CSP
from Corte_e_Empacomento_PC_flip import calcular_area, rotate_point, ponto_dentro_poligono
from nfp_teste import no_fit_polygon,combinar_poligonos,convex_partition
import numpy as np
from shapely.geometry import Polygon, MultiPoint
import turtle
import numpy as np
import random
import turtle
import copy



def GRASP(mIter, x,q1,q2,q3,rotacoes,pecas_inisol):
    for i in range(mIter):
        solucao = None
        ambiente = CSP(render=False, plot=False)
        I1, I2 = IniSOL(ambiente,rotacoes,pecas_inisol, False)
        fecho = ambiente.area_ocupada/ambiente.fecho
        I1_novo = LocSEARCH(I1,I2,x,q1,q2,q3,fecho, 6)
        
        if solucao is None or porcentagem_area(I1_novo) > porcentagem_area(solucao):
            solucao = copy.deepcopy(I1_novo)
            if len(solucao) == 12:
                return solucao

    return solucao


def LocSEARCH(I1, I2, x_dado, q1, q2, q3,fecho, pq):
    melhorou_global = True
    
    while melhorou_global:        
        melhorou_global = False
        area_global = porcentagem_area(I1)
        for idk in range(0,x_dado,1):
            tipo = random.choices(['permutar', 'trocar', 'adicionar'], weights=[q1, q2, q3])[0]
            
            I1_original = copy.deepcopy(I1)
            I2_original = copy.deepcopy(I2)
            ambiente = CSP(render=False, plot=False)
            if ambiente.render:
                ambiente.turtle_dados.clear()

            
            
            if tipo == 'permutar':
                I1_linha = trocar_elementos_aleatorios(I1)
            elif tipo == 'trocar':
                I1_linha, I2_linha = trocar_elementos_entre_listas(I1, I2)
            elif tipo == 'adicionar':
                I1_linha, I2_linha = adicionar_elemento(I1, I2)
            
            x, y = ambiente.cordenadas_area[3]
            index = ambiente.lista.index(I1_linha[0][-1])
            ambiente.acao(index, x, y, 0, False)
            I1_novo = []
            I1_novo = [I1_linha[0]]
            
            for peca in I1_linha[1:]:
                if peca[-1] in ambiente.lista:
                    index = ambiente.lista.index(peca[-1])
                    peca_posicao = PackS(ambiente, index, rotacoes, False)
                    if peca_posicao is not False:
                        I1_novo.append(peca)

                        ambiente.acao(peca_posicao[0], peca_posicao[1], peca_posicao[2], peca_posicao[3],peca_posicao[4],True)
                        


            
            
            #for k in range(len(ambiente.lista)):
             #   Possiveis_posicoes = []
              #  for peca in range(len(ambiente.lista)):
               #     peca_posicao = PackS(ambiente, peca, rotacoes)
                #    if peca_posicao is not False:
                 #       Possiveis_posicoes.append(peca_posicao)

                #if Possiveis_posicoes:
                 #   melhor_posicao = max(Possiveis_posicoes, key=lambda p: p[5])
                  #  I1_novo.append(melhor_posicao) 
                   # i+=1
                    #ambiente.acao(melhor_posicao[0], melhor_posicao[1], melhor_posicao[2], melhor_posicao[3], melhor_posicao[4])
                    

            area_nova = porcentagem_area(I1_novo)
            melhorou = area_nova > area_global
            if melhorou:
                I1 = copy.deepcopy(I1_novo)
                I2 = copy.deepcopy(ambiente.lista)
                area_global = area_nova
                fecho = ambiente.area_ocupada/ambiente.fecho
                melhorou_global = True
                I1_novo.clear()
                
            elif (abs(round(area_nova,3) - round(area_global,3)) < 0.03*round(area_global,3)) and ambiente.area_ocupada/ambiente.fecho > fecho + (0.01 * fecho):
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
            
            print(f"{tipo[0]}, {'M' if melhorou else 'NM'} {area_global:.3f} - {len(I1)} - {len(I2)}, {area_nova:.3f} - {12 - len(ambiente.lista)} - {len(ambiente.lista)}")
            if len(I1) == 12:
                if False:
                    turtle.clear()
                    ambiente = CSP(render=True, plot=False)
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
                print(f"---------------------------------------\nACABOU {len(ambiente.pecas_posicionadas)}\n--------------------------------------- ")
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
def trocar_elementos_entre_listas(lista1, lista2):
    I1 = lista1.copy()
    I2 = lista2.copy()
    if not I1 or not I2:
        return I1, I2  # Se alguma lista estiver vazia, retorna as mesmas listas
    
    # Seleciona um índice aleatório em I1 e I2
    idx1 = random.randint(0, len(I1) - 1)
    idx2 = random.randint(0, len(I2) - 1)
    
    # Troca os elementos entre as listas
    I1[idx1][-1], I2[idx2] = I2[idx2], I1[idx1][-1]
    
    return I1.copy(), I2.copy()

def trocar_elementos_aleatorios(I1):
    lista = I1.copy()
    idx1, idx2 = random.sample(range(len(lista)), 2)
    
    # Troca os elementos nas posições idx1 e idx2
    lista[idx1], lista[idx2] = lista[idx2], lista[idx1]
    
    return lista.copy()

def porcentagem_area(I1):
    area = 0
    for peca in I1:
        area+=calcular_area(peca[-1])
    area = round(area/(34*38), 3)
    return area    
    

def porcentagem_fecho(I1):
    lista = []
    for peca in I1:
        lista.append(peca[-1])
    
    return calcular_fecho_area(lista)

def IniSOL(ambiente, rotacoes, pq,draw):
    I1 = []  # Peças posicionadas
    I2 = []  # Peças não posicionadas
    # Inicializa a primeira peça na coordenada inicial (x, y) do ambiente
    grau = 0
    x, y = ambiente.cordenadas_area[3 - grau]
    index = random.randint(0,len(ambiente.lista) -1)
    peca1 = [index,x,y,grau,0,ambiente.lista[index]]
    
    ambiente.acao(index, x, y, grau, False)
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
        
        # Tentar posicionar cada peça restante
        for peca in range(len(ambiente.lista)):
            peca_posicao = PackS(ambiente, peca, rotacoes,draw)
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

        
        
def PackS(ambiente, peca, rotacoes,draw):
    nfp = []
    lista = []
    melhores_posicoes = []
    

    for grau_indice in rotacoes:
        #print(peca, len(ambiente.lista))
        melhor_posicao = ExSEARCH(ambiente, peca, grau_indice,draw)
        if melhor_posicao is not None:
            melhores_posicoes.append(melhor_posicao)
    
    if melhores_posicoes:
        melhor = best_posicao(melhores_posicoes) 
        return melhor
    else:
        return False





def ExSEARCH(ambiente, peca, grau_indice, draw):
    possiveis_posicoes = []
    nfp = []
    flip = False
    pol_posicionar = ambiente.nova_lista[peca]
    graus = [0,90,180,270]
    grau = graus[grau_indice]
    pontos_pol = [(int(rotate_point(cor[0], cor[1], grau)[0]), int(rotate_point(cor[0], cor[1], grau)[1])) for cor in pol_posicionar]
    lista = []
    for item in ambiente.pecas_posicionadas:
        convex_parts = convex_partition(item)

        nfps_convx = []
        for convex in convex_parts:
            nfps_convx.append(Polygon(NoFitPolygon(convex, pontos_pol)))


        nfp_final = list(combinar_poligonos(nfps_convx).exterior.coords)
        lista.append(Polygon(nfp_final))

    nfp = list(combinar_poligonos(lista).exterior.coords)


            
    for x,y in nfp:
        #for flip in [True, False]:

                pontos = [(int(x + rotate_point(cor[0], cor[1], grau)[0]), int(y + rotate_point(cor[0], cor[1], grau)[1])) for cor in pol_posicionar]
                posicao_valida = True
                for x1,y1 in pontos:
                    if not ponto_dentro_poligono(x1,y1,ambiente.cordenadas_area, True):
                        posicao_valida = False
                if posicao_valida:
                    if draw:
                        for i in range(100):
                            ambiente.draw_click(x,-y,grau_indice,peca,flip)
                    pecas_posicionadas = copy.deepcopy(ambiente.pecas_posicionadas)
                    pecas_posicionadas.append(pontos)

                    area_ocupada = calcular_area_preenchida(pecas_posicionadas)
                    area_FC = calcular_fecho_area(pecas_posicionadas)
                    area_FR = area_fecho_retangular(pecas_posicionadas)

                    fecho_convexo = (1 - (area_ocupada/area_FR))
                    fecho_convexo_max = area_ocupada/area_FC
                    fecho_retangular_max = area_ocupada/area_FR
                    fechoR = (1 - (area_ocupada/area_FR))
                    nfp_meu = (fecho_convexo_max**10)*(area_FC/(340*380))

                    pecas_posicionadas = []

                    possiveis_posicoes.append([peca,x,y,grau_indice,flip,round(fecho_convexo_max,2),ambiente.lista[peca]])
                    
    if possiveis_posicoes:
        melhor_posicao = best_posicao(possiveis_posicoes)
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
    polyA = Polygon(poligono_fixo)
    polyB = Polygon(poligono_orbital)

    return interpolar_pontos_poligono(no_fit_polygon(polyA, polyB))

import random

def adicionar_elemento(I1, I2):
    lista1 = copy.deepcopy(I1)
    lista2 = copy.deepcopy(I2)
    if not lista2:
        print("A lista2 está vazia. Nenhum item pode ser adicionado.")
        return lista1, lista2
    
    # Escolhe um elemento aleatório de lista2
    elemento = random.choice(lista2)
    lista2.remove(elemento)
    
    # Escolhe uma posição aleatória em lista1 para inserir o elemento
    posicao = random.randint(0, len(lista1) - 1)  
    
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

def append_to_txt(data, filename="resultados_GRASP_original.txt"):
    with open(filename, "a") as file:
        file.write(f"Iteracao {data['iteration']}: {data['num_pieces']} pecas, {data['area_percent']*100}% , {data['execution_time']:.2f} segundos\n")

rotacoes = [0, 1, 2]

vizinhanca = 15
q1 = 50
q2 = 20
q3 = 30
pecas_inisol = 0
max_iter = 30
import time
import math
import turtle

def append_to_txt(data, filename):
    with open(filename, "a") as file:
        file.write(f"Iteracao {data['iteration']}: {data['num_pieces']} pecas, {data['area_percent']}% area, {data['execution_time']:.2f} segundos\n")

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

def append_to_txt(data, filename):
    with open(filename, "a") as file:
        file.write(f"Iteração {data['iteration']}: {data['num_pieces']} peças, {data['area_percent']}% área, {data['execution_time']:.2f} segundos\n")

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

rotacoes = [0, 1, 2]

for j in range(3):
    with open(filename, "a") as file:
        file.write("---------------------------------------------------------------------------------------------\n")
        file.write(f"GRASP ORIGINAL: iteracoes maxima {max_iter}, {vizinhanca} vizinhancas, {12 - j} pecas inisol, q1,q2,q2 = {q1},{q2},{q3}\n")

    tempos = []

    for i in range(1000):
        start_time = time.time()  # Inicia a medição do tempo
        pecas = GRASP(max_iter, vizinhanca, q1, q2, q3, rotacoes, j)
        end_time = time.time()  # Termina a medição do tempo
        execution_time = end_time - start_time  # Calcula o tempo de execução
        print(f"Iteracao {i+1}: {len(pecas)} pecas, {porcentagem_area(pecas)} , {execution_time:.2f} segundos")

        resultado = {
            "iteration": i + 1,
            "num_pieces": len(pecas),
            "area_percent": porcentagem_area(pecas) * 100,
            "execution_time": execution_time
        }
        
        # Armazena o tempo de execução
        tempos.append(execution_time)
        
        # Salvar os resultados em um arquivo .txt a cada iteração
        append_to_txt(resultado, filename)

    # Calcula a média e o desvio padrão dos tempos
    media, desvio_padrao = calcular_tempo_medio_desvio(tempos)
    
    # Salvar a média e o desvio padrão no arquivo
    append_summary_to_txt(media, desvio_padrao, filename)

turtle.done()




    


    
        