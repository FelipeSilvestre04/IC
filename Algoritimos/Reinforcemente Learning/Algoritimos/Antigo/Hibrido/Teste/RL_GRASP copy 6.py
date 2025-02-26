from Cutting_Stock_Problem.Ambiente.Main.CSP_concavo import CSP
from Corte_e_Empacomento import calcular_area, pontos_entre_vertices, polygon_features, max_coordinates_count, matriz_cordenadas
import numpy as np
import random
import math
import time
import gymnasium as gym
from gymnasium import spaces
import numpy as np
from stable_baselines3 import PPO,DQN
from stable_baselines3.common.env_checker import check_env
import matplotlib.pyplot as plt
from stable_baselines3.common.callbacks import BaseCallback
import turtle
from geradores import Area
from botao import Botao
import math
import numpy as np
from scipy.spatial import ConvexHull
import time
import random
import csv
import copy
import matplotlib.pyplot as plt
from matriz_vertices import Matriz
from shapely.geometry import Polygon
import os
from stable_baselines3.common.callbacks import EvalCallback
from stable_baselines3 import PPO
from stable_baselines3.common.vec_env import SubprocVecEnv
import tensorflow as tf
import torch
import torch.nn as nn
from stable_baselines3.common.torch_layers import BaseFeaturesExtractor
from GRASP_RL import PackS, pre_processar_NFP,IniSOL

class envCSP(gym.Env):
    def __init__(self):

        self.ambiente = CSP(dataset='fu', render=False, plot=False)
        self.rotacoes = self.rot(self.ambiente.pecas)

        self.tabela_nfps = pre_processar_NFP(self.rotacoes, self.ambiente.nova_lista)

        self.observation_space = spaces.Dict({
            'cutting_area': spaces.Box(low=0, high=1, shape=(self.ambiente.altura, self.ambiente.base), dtype=np.int8),
            'pieces': spaces.Box(low=-1, high=42, shape=(len(self.ambiente.lista), 4, 2), dtype=np.int8)
        })
           
        self.action_space = spaces.MultiDiscrete([self.ambiente.max_pecas, self.ambiente.max_pecas])
        
        self.acao = 0
        self.epoca = 0

        self.preenchimentos = []
        self.recompensas = []
        self.pecas_colocadas = []

        self.MediaPreenchimentos = []
        self.MediaRecompensas = []
        self.Mpecas_colocadas = []
 
        self.lista_recompensa = []

     

    def reset(self, seed=None):
        if seed is not None:
            np.random.seed(seed)
            set_seed(seed)

        self.acao = 0
        self.epoca += 1

        self.ambiente = CSP(dataset='fu', render=False, plot=False)

        I1,I2 = IniSOL(self.ambiente, self.rotacoes,0,False,self.tabela_nfps)

        self.pecas_iniciais = len(I1)
        self.pecas_posicionadas = len(I1)
        self.melhorou_global = 0

        self.pinicial = round(self.ambiente.area_ocupada/self.ambiente.area, 2)
        self.porcentagem = round(self.ambiente.area_ocupada/self.ambiente.area, 2)
        self.fecho = round(self.ambiente.area_ocupada/self.ambiente.fecho, 2)

        I1 = [peca[-1] for peca in I1]
        self.lista_recompensa = []
        
        pecas = []

        for piece in I1:
            pecas.append(piece)

        for piec in I2:
            pecas.append(piec)

        self.pecas_locsearch = copy.deepcopy(pecas)
        max_len = 4
        for pol in pecas:
            if len(pol) < max_len:
                for i in range(max_len - len(pol)):
                    pol.append((-1,-1))
                
        self.pecas = np.array(pecas, dtype=np.int8)
        self.matriz = np.array(self.ambiente._matriz, dtype=np.int8)

        self.estado = {
            'cutting_area': self.matriz,
            'pieces': self.pecas
        }

        return self.estado, {}

    def step(self, action):
        self.ambiente = CSP(dataset='fu', render=False, plot=False)
        self.acao += 1
        
        idx1,idx2 = action[0],action[1]

        teste_lista = copy.deepcopy(self.pecas_locsearch)


        peca_mover = copy.deepcopy(teste_lista[idx1])
        teste_lista.pop(idx1)
        teste_lista.insert(idx2,peca_mover)


        x, y = self.ambiente.cordenadas_area[3]

        minx = min([x for x,y in teste_lista[0]]) 
        miny = min([y for x,y in teste_lista[0]]) 

        minx*=self.ambiente.Escala
        miny*=self.ambiente.Escala

        index = self.ambiente.lista.index(teste_lista[0])
        self.ambiente.acao(index, x - minx, y - miny, 0, False)



            
        for peca in teste_lista[1:]:
            #if peca in self.ambiente.lista:
            index = self.ambiente.lista.index(peca)
            peca_posicao = PackS(self.ambiente, index, self.rotacoes, False,self.tabela_nfps)
            if peca_posicao is not False:
                self.ambiente.acao(peca_posicao[0], peca_posicao[1], peca_posicao[2], peca_posicao[3],peca_posicao[4],True)
                        
        area_nova = round(self.ambiente.area_ocupada/self.ambiente.area, 2)
        fecho_novo = round(self.ambiente.area_ocupada/self.ambiente.fecho, 2)
        melhorou = area_nova > self.porcentagem
        
        if melhorou:
            self.pecas_locsearch = copy.deepcopy(teste_lista)
            teste_lista = []
            self.recompensa = (area_nova - self.porcentagem)*100
            self.porcentagem = area_nova
            self.fecho = fecho_novo
            self.pecas_posicionadas = len(self.ambiente.pecas_posicionadas)
            self.melhorou_global = 0
               
        elif (abs(round(area_nova,3) - round(self.porcentagem,3)) < 0.01*round(self.porcentagem,3)) and self.ambiente.area_ocupada/self.ambiente.fecho > self.fecho:
            self.pecas_locsearch = copy.deepcopy(teste_lista)
            teste_lista = []
            self.recompensa = fecho_novo - self.fecho
            self.porcentagem = area_nova
            self.fecho = fecho_novo
            self.pecas_posicionadas = len(self.ambiente.pecas_posicionadas)
            self.melhorou_global = 0               
        else:

            if idx1 == idx2:
                self.recompensa = -1
            else:
                self.recompensa = (area_nova - self.porcentagem)
            
            self.melhorou_global += 1
            teste_lista = []
        
        if len(self.ambiente.pecas_posicionadas) == self.ambiente.max_pecas:
            self.recompensa = 100
            self.done = True
        else:
            self.done = False
        self.lista_recompensa.append(round(self.recompensa,2))           
        #print(self.acao, action, self.recompensa)

        if not self.done:
            if self.melhorou_global < 25:
                self.done = False
                #print(self.acao)
            else:
                #print(self.acao)
                self.done = True
                print(f"F:{round(self.porcentagem * 100, 2)}% , {self.pecas_posicionadas}/12 --- I:{self.pinicial * 100}% , {self.pecas_iniciais}/12 ||| R:{[res for res in self.lista_recompensa if res > 0]}, SR:{round(sum(self.lista_recompensa),2)} | NA {self.acao}")

        
        

        if self.done:
            self.preenchimentos.append((self.ambiente.area_ocupada/self.ambiente.area) * 100)
            self.pecas_colocadas.append(len(self.ambiente.pecas_posicionadas))
            self.recompensas.append(sum(self.lista_recompensa))

            if (self.epoca % 100) == 0:
                self.calcular_medias()

        pecas = copy.deepcopy(self.pecas_locsearch)
        max_len = 4 
        for pol in pecas:
            if len(pol) < max_len:
                for i in range(max_len - len(pol)):
                    pol.append((-1,-1))
            


        self.pecas = np.array(pecas, dtype=np.int8)
        self.matriz = np.array(self.ambiente._matriz, dtype=np.int8)

        self.estado = {
            'cutting_area': self.matriz,
            'pieces': self.pecas
        }

        return self.estado, float(self.recompensa), self.done, False, {}
    
    def calcular_medias(self):

        print(f"Media dos ultimos 100 preenchimentos : {round((sum(self.preenchimentos)) / (len(self.preenchimentos)), 2)}%, {round((sum(self.pecas_colocadas)) / (len(self.pecas_colocadas)), 2)}/12")
        print(f"Media das ultimas 100 recompensas: {round((sum(self.recompensas)) / (len(self.recompensas)), 2)}")
        self.MediaPreenchimentos.append(round((sum(self.preenchimentos)) / (len(self.preenchimentos)), 2))
        self.MediaRecompensas.append(round((sum(self.recompensas)) / (len(self.recompensas)), 2))
        self.Mpecas_colocadas.append(round((sum(self.pecas_colocadas)) / (len(self.pecas_colocadas)), 2))
 
        
        if len(self.MediaPreenchimentos) > 10:
            self.MediaPreenchimentos.pop(0)
            self.MediaRecompensas.pop(0)
            self.Mpecas_colocadas.pop(0)
        
        print(f"ultimas 10 preenchimentos{self.MediaPreenchimentos}")
        print(f"ultimas 10 recompensas{self.MediaRecompensas}")
        self.preenchimentos = []
        self.recompensas = []
        self.pecas_colocadas = []
    
    def seed(self, seed=None):
        self.np_random, seed = gym.utils.seeding.np_random(seed)
        return [seed]
    
    def rot(self,instancia):
        if instancia == "fu":
            return [0, 1, 2]  # 0, 90, 180
        elif instancia == "jackobs1":
            return [0, 1, 2]  # 0, 90, 180
        elif instancia == "jackobs2":
            return [0, 1, 2]  # 0, 90, 180
        elif instancia == "shapes0":
            return [0]  # Apenas 0
        elif instancia == "shapes1":
            return [0, 2]  # 0, 180
        elif instancia == "shapes2":
            return [0, 2]  # 0, 180
        elif instancia == "dighe1":
            return [0]  # Apenas 0
        elif instancia == "dighe2":
            return [0]  # Apenas 0
        elif instancia == "albano":
            return [0, 2]  # 0, 180
        elif instancia == "dagli":
            return [0]  # Apenas 0
        elif instancia == "mao":
            return [0, 1, 2]  # 0, 90, 180
        elif instancia == "marques":
            return [0, 1, 2]  # 0, 90, 180
def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)  # Para múltiplas GPUs
    torch.backends.cudnn.deterministic = True  # Para resultados determinísticos
    torch.backends.cudnn.benchmark = False  # Evitar variações do cuDNN

def make_env(env_class, rank, seed=0):
    def _init():
        env = env_class()
        env.seed(seed + rank)
        return env
    return _init

def linear_schedule(initial_value, final):
    def func(progress):
        if initial_value * (progress) <= final:
            return final
        else:
            return initial_value * (progress)
    return func
from stable_baselines3.common.vec_env import DummyVecEnv
if __name__ == '__main__':
    env = envCSP()
    #num_envs = 16


    #envs = DummyVecEnv([make_env(env, i) for i in range(num_envs)])


    log_dir = "./ppo_csp_tensor/"

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    policy_kwargs = dict(
        net_arch=dict(pi=[128, 128], vf=[128, 128]),
    )

    model = PPO(
        "MultiInputPolicy",  # Política para entradas múltiplas
        env,
        verbose=1,
        tensorboard_log="log_dir",  # Substitua pelo caminho desejado
        n_steps=120,               # Maior número de passos para capturar mais transições
        batch_size=30,             # Tamanho de lote ajustado para estabilidade
        n_epochs=15,                # Redução para evitar overfitting
        gamma=0.99,                 # Considerar recompensas de longo prazo
        learning_rate=linear_schedule(3e-4, 1e-4),  # Taxa de aprendizado adaptativa
        ent_coef=0.01,              # Coeficiente de entropia para aumentar exploração
        clip_range=0.1,             # Range de clipping para estabilidade
        normalize_advantage=True    # Normalização das vantagens
    )


  
    model.learn(total_timesteps=1_000_000)
    import os
    from datetime import datetime

    caminho_base = "/home/fsilvestre/Cutting_Stock_Problem"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    caminho_para_salvar = f"{caminho_base}_{timestamp}.zip"
    os.makedirs(os.path.dirname(caminho_para_salvar), exist_ok=True)
    model.save(caminho_para_salvar)


    

  


