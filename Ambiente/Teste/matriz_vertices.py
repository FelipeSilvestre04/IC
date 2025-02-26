import numpy as np
from PIL import Image, ImageDraw
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def transformar_coordenadas(x, y, x_min, x_max, y_min, y_max, width, height):
    x_scaled = (x - x_min) / (x_max - x_min) * width
    y_scaled = (y - y_min) / (y_max - y_min) * height
    y_scaled = height - y_scaled

    return x_scaled, y_scaled
class Matriz:
    def __init__(self, base, altura, x, y):
        self._base = base
        self._altura = altura
        self.x = x
        self.y = y
        self._img = Image.new('L', (self._base, self._altura), 0)
        self._matriz = np.array(self._img)
        self._fig, self._ax = plt.subplots()
        self._poligonos = []

    def _criar_poligono(self, vertices):
        self._poligonos.append(vertices)
        self._atualizar_matriz()

    def _remover_poligono(self):
        if self._poligonos:
            self._poligonos.pop()
            self._atualizar_matriz()

    def _atualizar_matriz(self):
        self._img = Image.new('L', (self._base, self._altura), 0)
        for poligono in self._poligonos:
            poligono_transformado = [transformar_coordenadas(x, y, self.x, self.x + self._base, self.y - self._altura, self.y, self._base, self._altura) for x, y in poligono]
            ImageDraw.Draw(self._img).polygon(poligono_transformado, outline=1, fill=1)
        self._matriz = np.array(self._img)

    def _plot(self):
        self._ax.clear()
        self._ax.imshow(self._matriz, cmap='Purples')
        plt.pause(0.001)

    def fechar_janela(self):
        plt.close(self._fig)

