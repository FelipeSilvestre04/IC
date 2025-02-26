import turtle 
from geradores import Area, Gerador, posicionar, on_click
from matriz import Matriz
from botao import Botao
import pyautogui
import math
import time
import gymnasium as gym
from gymnasium import spaces
import numpy as np
from stable_baselines3 import PPO
from stable_baselines3.common.env_checker import check_env
import matplotlib.pyplot as plt
from stable_baselines3.common.callbacks import BaseCallback

def transformar_coordenadas_inversa(x, y, x_min, x_max, y_min, y_max, width, height):
    x_scaled = (x / width) * (x_max - x_min) + x_min
    y_scaled = ((height - y) / height) * (y_max - y_min) + y_min

    return x_scaled, y_scaled

class env2KP(gym.Env):
    def __init__(self, base, altura, x, y, grids, numero_peças, tipo):
        super(env2KP, self).__init__()
        self.base = base
        self.altura = altura
        self.x = x
        self.y = y
        self.grid = grids
        self.numero = numero_peças
        self.tipo = tipo

        self.done = False
        
        self.pecas = []


        for s in Gerador._geracao:
            numeros_string = s.split(", ")
            lista_inteiros = [int(numero_string) for numero_string in numeros_string]
            self.pecas.append(lista_inteiros)
        
        self.observation_space = spaces.Box(low=0, high=1, shape=(400*600,5), dtype=np.float32)
        self.action_space = spaces.MultiDiscrete([600, 400])
     
        
    def reset(self, seed=None):
        self.done = False
        if seed is not None:
            np.random.seed(seed)
            
        turtle.clearscreen()
        turtle.reset()
        self.matriz = Matriz(self.base,self.altura,self.x,self.y)
        self.tecido = Area(self.base,self.altura,self.grid)
        self.bot = Botao(350, 300, 200, 150)
        self.text_turtle = turtle.Turtle()
        self.fecho_turtle = turtle.Turtle()
        
        Area.AreaTotal = []
        Area.AreaGridN = [0.00,0.00,0.00,0.00,0.00,0.00]
        Area.AreaOcupada = 0
        Area.CordenadasProibidas = []
        Area.CordenadasGreedy = []
        Area.AreasGredys = []
        Area.GreedysOcupados = []
        Area.Q_total = -6

        Gerador.Colocou = False
        Gerador.N = False
        Gerador._count = 0
        Gerador._geracaopop = []
        Gerador._CordenadasProibidas = []
        Gerador._desenhando = False
        Gerador._poligonos_posicionados = []
        Gerador._Desenahdos = []
        Gerador._PolEGrid = []
        
        turtle.tracer(0)
        self.tecido.area(self.x, self.y)
        self.tecido.greedy(self.x,self.y)
        posicionar(self.numero, self.tipo, True)
        turtle.tracer(1)
        matriz_com_canal = self.matriz._matriz[..., np.newaxis]
  

        caracteristicas_peca = []


        numeros_string = Gerador._geracao[-1].split(", ")
            
        
        lista_inteiros = [int(numero_string) for numero_string in numeros_string]
            

        caracteristicas_peca.append(lista_inteiros)
        caracteristicas_peca_1D = np.array(caracteristicas_peca).flatten()

        # Remodele os arrays para terem a forma correta
        matriz_reshaped = self.matriz._matriz.flatten().reshape((240000, 1))
        caracteristicas_peca_reshaped = caracteristicas_peca_1D.reshape((1, 4))

        # Repita as características da peça 240000 vezes ao longo do eixo 0
        caracteristicas_peca_repeated = np.repeat(caracteristicas_peca_reshaped, 240000, axis=0)

# Agora você pode concatenar os dois arrays ao longo do eixo 1
        estado_com_peca = np.concatenate([matriz_reshaped, caracteristicas_peca_repeated], axis=1)


        return estado_com_peca, {}

    
    def step(self, action):
        x_transformado, y_transformado = transformar_coordenadas_inversa(int(action[0]), int(action[1]), self.x, 250, self.y, -100, self.base, self.altura)
        on_click(x_transformado,y_transformado,self.matriz,self.text_turtle, self.fecho_turtle, self.base, self.altura, self.grid, self.bot)
        if not Gerador.Colocou:
            Gerador._geracao.pop()
        if len(Gerador._geracao)>0:
            matriz_com_canal = self.matriz._matriz[..., np.newaxis]
            caracteristicas_peca = []


            numeros_string = Gerador._geracao[-1].split(", ")
                
            
            lista_inteiros = [int(numero_string) for numero_string in numeros_string]
                

            caracteristicas_peca.append(lista_inteiros)
            caracteristicas_peca_1D = np.array(caracteristicas_peca).flatten()

            # Remodele os arrays para terem a forma correta
            matriz_reshaped = self.matriz._matriz.flatten().reshape((240000, 1))
            caracteristicas_peca_reshaped = caracteristicas_peca_1D.reshape((1, 4))

            # Repita as características da peça 240000 vezes ao longo do eixo 0
            caracteristicas_peca_repeated = np.repeat(caracteristicas_peca_reshaped, 240000, axis=0)

    # Agora você pode concatenar os dois arrays ao longo do eixo 1
            estado_com_peca = np.concatenate([matriz_reshaped, caracteristicas_peca_repeated], axis=1)
                
        if len(Gerador._geracao) > 1:
            self.done = False
        else:
            self.done = True
            return estado_com_peca, Area.Last_Reward[-1], self.done,False,{}
 
        
        return estado_com_peca, Area.Last_Reward[-1], self.done,False,{}
        
env = env2KP(600,400,-350,300,6,100,"M")
model = PPO('MlpPolicy', env, seed=1, verbose=1)




import matplotlib.pyplot as plt
from IPython import display

plt.ion()

def plot(scores, mean_scores):
    display.clear_output(wait=True)
    display.display(plt.gcf())
    plt.clf()
    plt.title('Training...')
    plt.xlabel('Epocas')
    plt.ylabel('Preenchimento, em %')
    plt.plot(scores)
    plt.plot(mean_scores)
    plt.ylim(ymin=0)
    plt.show(block=False)
    plt.pause(.1)
    
class CustomCallback(BaseCallback):
    def _on_step(self):
        if 'done' in self.locals and self.locals['done']:  # Se o episódio terminou
            Gerador._preenchimentos.append(round((Area.AreaOcupada/(600*400))*100, 1))
            Area.recompensas.append(Area.Last_Reward[-1])
            print(Area.AreaOcupada, Area.recompensas[-1])
            plot(Gerador._preenchimentos,Area.recompensas)
        return True


# Crie uma instância do callback
callback = CustomCallback()

# Passe o callback para a função de treinamento
model.learn(total_timesteps=10000, callback=callback)
    

