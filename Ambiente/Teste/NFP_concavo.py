import numpy as np
from shapely.geometry import Polygon, MultiPoint
from shapely.ops import unary_union, triangulate
import turtle

from numba import njit
import math
from shapely.geometry import MultiPoint, Polygon, MultiPolygon
from shapely.ops import triangulate

import numpy as np
from shapely.geometry import MultiPoint, Polygon
import numpy as np
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon, MultiPoint

def convex_hull_to_concave(points):
    """
    Recebe uma lista de pontos e retorna o menor polígono côncavo que os abrange.
    
    Parâmetros:
    - points: lista de tuplas (x, y) representando os pontos.
    
    Retorna:
    - Um objeto Polygon representando o polígono côncavo que abrange todos os pontos.
    """
    # Calcula o fecho convexo dos pontos
    hull = ConvexHull(points)
    
    # Extrai os vértices do fecho convexo
    vertices = [points[idx] for idx in hull.vertices]
    
    # Cria um polígono a partir dos vértices
    convex_polygon = Polygon(vertices)
    
    # Encontra o polígono côncavo que envolve todos os pontos
    concave_polygon = convex_polygon.envelope
    
    return convex_polygon

def minkowski_sum(polygonA, polygonB, ref_point_B):
    """Calcula a soma de Minkowski entre dois polígonos usando um ponto de referência de B e retorna apenas os pontos da soma."""
    
    sum = []
    sum_conv = []
    for p1 in polygonA.exterior.coords:
        sum_points = []
        for p2 in polygonB.exterior.coords:
            # A soma é ajustada para que p2 seja relativo ao ponto de referência ref_point_B
            sum_points.append((p1[0] + p2[0] - ref_point_B[0], p1[1] + p2[1] - ref_point_B[1]))
            sum_conv.append((p1[0] + p2[0] - ref_point_B[0], p1[1] + p2[1] - ref_point_B[1]))
        sum.append(Polygon(sum_points))
        
    
    # Cria um polígono a partir dos pontos da soma
    #print(len(sum))
    #if polygonB.equals(polygonB.convex_hull) and polygonA.equals(polygonA.convex_hull):
        #print("convexo")
    if True:
        return MultiPoint(sum_conv).convex_hull
        
    else:
        return Polygon(combinar_poligonos(sum))

def no_fit_polygon(polyA, polyB):

    ref_point_B = polyB.exterior.coords[0]  # Usa o primeiro vértice de B como ponto de referência
    
    # Inverte as coordenadas de B para simular a diferença de Minkowski
    polyB_inverted = Polygon([(-x, -y) for x, y in polyB.exterior.coords])
    


    nfp = minkowski_sum(polyA, polyB_inverted, ref_point_B)

    
    # Retorna o contorno do NFP combinado
    return [(x + ref_point_B[0], y + ref_point_B[1]) for x,y in list(nfp.exterior.coords)]
def combinar_poligonos(poligonos):
    """
    Recebe uma lista de polígonos e retorna um único polígono que representa
    a união sem sobreposições dos contornos dos polígonos.
    
    Parâmetros:
    - poligonos: lista de objetos Polygon representando cada polígono.

    Retorna:
    - Um polígono que representa a união dos contornos dos polígonos.
    """
    # Realiza a união geométrica dos polígonos para obter apenas a borda externa
    poligono_combinado = unary_union(poligonos)
   
    
    # Verifica se o resultado é um único polígono ou uma coleção
    if poligono_combinado.geom_type == 'Polygon':
        return poligono_combinado  # Retorna o polígono unido diretamente
    elif poligono_combinado.geom_type == 'MultiPolygon':
        # Se for MultiPolygon, converte para um único polígono unindo as partes
        return unary_union([p for p in poligono_combinado.geoms])
    else:
        print(poligonos)
        print(poligono_combinado)
        raise ValueError("Erro: a geometria resultante não é um polígono válido.")

from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import triangulate

import numpy as np
from shapely.geometry import Polygon, LineString

from shapely.geometry import Polygon, LineString
import numpy as np

from shapely.geometry import Polygon, LineString
@njit
def interpolar_pontos_poligono(vertices, pontos_entre_vertices=0):
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

from Corte_e_Empacomento_PC_flip import tratar_lista, ler_poligonos
def desenhar(cords, x_offset, y_offset,color, fill = False, draw = False):
        if not draw:
            p = turtle.Turtle()
            p.penup()
            if color is not None:
                p.color(color)
            if fill:
                p.begin_fill()
            p.goto(cords[0][0] * 50 + x_offset, cords[0][1] * 50 + y_offset)
            p.pendown()
            for x, y in cords:
                p.goto(x * 50 + x_offset, y * 50 + y_offset)
            p.goto(cords[0][0] * 50 + x_offset, cords[0][1] * 50 + y_offset)
            if fill:
                p.end_fill()
            p.hideturtle()
        elif draw:
            fill = True
            
            p = turtle.Turtle()
            p.penup()
            if color is not None:
                p.color(color)
            if fill:
                p.begin_fill()
            p.goto(cords[0][0] * 50 + x_offset, cords[0][1] * 50 + y_offset)
            p.pendown()
            for x, y in cords:
                p.goto(x * 50 + x_offset, y * 50 + y_offset)
            p.goto(cords[0][0] * 50 + x_offset, cords[0][1] * 50 + y_offset)
            if fill:
                p.end_fill()
            p.hideturtle()
            for i in range(200000):
                pass
            p.clear()
            
import math
def ajustar_poligono(poligono):
    # Encontrar o primeiro vértice
    x_min = poligono[0][0]
    y_min = poligono[0][1]
    
    # Transladar todos os vértices para garantir que o primeiro ponto seja (0,0)
    poligono_ajustado = [(x - x_min, y - y_min) for (x, y) in poligono]
    
    return poligono_ajustado


def filter_triangles_inside_polygon(polygon, triangles):
    """
    Filtra os triângulos que estão dentro do polígono original.
    
    Parâmetros:
    polygon (Polygon): Objeto Polygon representando o polígono original.
    triangles (list of tuples): Lista de triângulos, onde cada triângulo é representado por uma lista de 3 pontos.
    
    Retorna:
    list of tuples: Lista de triângulos que estão dentro do polígono original.
    """
    valid_triangles = []
    for triangle in triangles:
        triangle_polygon = Polygon(triangle)
        if polygon.contains(triangle_polygon):
            valid_triangles.append(list(triangle.exterior.coords))
    return valid_triangles

from shapely.geometry import Polygon
from shapely.ops import triangulate

def triangulate_shapely(points):
    """
    Triangula um polígono usando a função `triangulate()` do Shapely e retorna a lista de vértices únicos.
    
    Parâmetros:
    points (list of tuples): Lista de coordenadas (x, y) dos pontos que formam o polígono.
    
    Retorna:
    list of tuples: Lista de vértices únicos.
    """
    # Cria um objeto Polygon a partir dos pontos
    polygon = Polygon(points)
    
    # Verifica se o polígono é convexo
    if polygon.equals(polygon.convex_hull):
        #print("O polígono é convexo.")
        return [points]  # Retorna o próprio polígono como uma lista única de vértices
    
    else:
        # Triangula o polígono usando a função `triangulate()`
        triangles = list(triangulate(polygon))
        
        # Filtro dos triângulos que estão dentro do polígono original
        return filter_triangles_inside_polygon(polygon, triangles)

if __name__ == "__main__":
    pol2 = [(0, 0), (4,0),(4,5),(0,5)]
    
    poly4 = [(4,2),(10,2),(10,5),(4,5)]
    pol3 = [(10,0),(14,0),(14,5),(10,5)]
    # Exemplo de uso:
    #polyB = Polygon(pol1)  # Polígono A
    polyA = Polygon(pol2)  # Polígono B
    polyC = Polygon(pol3)
    polyD = Polygon(poly4)
    polygon_points = [(0, 0), (3, 0), (3, 6),(2,6),(2,3),(0,3)]
    pecas = 'jackobs1'
    
    Escala = 1
    lista_original = ler_poligonos('C:/Users/felip/OneDrive/Documentos/IC 2/' + pecas +'.dat')
    nova_lista_original, nova_lista_completa_original = tratar_lista(lista_original, Escala)
    print(nova_lista_original[-8])
    convex_parts = triangulate_shapely(ajustar_poligono(nova_lista_original[0]))
    pol1 = [(0, 0), (2, 0),(2,3),(4,3),(4,6),(0,6)]
    pol = [(0, 0), (2, 0),(2,6),(0,6)]
    po = ajustar_poligono(nova_lista_original[-7])
    print(po)
    cv = triangulate_shapely(pol1)

    polyB = Polygon(pol1)
    print(len(convex_parts))
    #print(convex_parts)
    

    nfps = []
    nfps2 = []
    nfps3 = []
    #nfp = no_fit_polygon(polyA, polyB)
    #nfp2 = no_fit_polygon(polyC, polyB)
    #nfp3 = no_fit_polygon(polyD, polyB)
    cores = ["red","green","black","pink","orange","blue"]
    i=0
    
    for nfp in convex_parts:
            print(nfp)
            desenhar(nfp,0,0,cores[i],True)
            i+=1
            if i == 6:
                i = 0
            nfps.append(Polygon(no_fit_polygon(Polygon(nfp), Polygon(ajustar_poligono(po)))))
            #nfps2.append(Polygon(no_fit_polygon(Polygon(nfp), Polygon(ajustar_poligono(pol)))))
            #nfps3.append(Polygon(no_fit_polygon(Polygon(nfp), Polygon(ajustar_poligono(pol)))))
            
            #desenhar(list((nfps[-1]).exterior.coords), 0, 0,"blue")


    nfp_final = interpolar_pontos_poligono(list(combinar_poligonos(nfps).exterior.coords),0)
    #nfp2_final = interpolar_pontos_poligono(list(combinar_poligonos(nfps2).exterior.coords),0)
    #nfp3_final = interpolar_pontos_poligono(list(combinar_poligonos([Polygon(nfp_final),Polygon(nfp2_final)]).exterior.coords),0)

    print("NFP:", nfp_final)

    # Função para desenhar polígonos no Turtle

    # Desenhar a forma original (B) na posição de referência
    #desenhar(polygon_points,0,0)



    #desenhar(pol2, 0, 0)  # Posição original de B
    #desenhar(pol3, 0, 0)
    #desenhar(poly4,0,0)
    # Desenhar o NFP
    desenhar(nfp_final, 0, 0,"blue")  # NFP centrado
    #desenhar(nfp2_final, 0, 0, "red")
    #desenhar(nfp3_final, 0, 0, "green")

    # Desenhar várias cópias do polígono B nas posições calculadas pelo NFP
    for x, y in nfp_final:
    
        desenhar(po, x * 50, y * 50, "green", draw=True)  # Mover B para as posições do NFP
    
    turtle.done()




