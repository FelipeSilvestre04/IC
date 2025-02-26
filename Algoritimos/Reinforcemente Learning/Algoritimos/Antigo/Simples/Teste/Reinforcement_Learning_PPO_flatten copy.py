from Corte_e_Empacomento import CSP
from Corte_e_Empacomento import calcular_area, pontos_entre_vertices, polygon_features, max_coordinates_count, matriz_cordenadas
import numpy as np
import random
import math
import time
import gymnasium as gym
from gymnasium import spaces
import numpy as np
from stable_baselines3 import PPO
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

class envCSP(gym.Env):
    def __init__(self):
        self.ambiente = CSP(dataset='fu', render=False, plot=False)
        cutting_area_shape = self.ambiente.altura * self.ambiente.base
        self.observation_space = spaces.Dict({
            'cutting_area': spaces.Box(low=0, high=1, shape=(cutting_area_shape,), dtype=np.int8),  
            
            'pieces': spaces.Box(low=-1, high=42, shape=(len(self.ambiente.lista) * max_coordinates_count(self.ambiente.nova_lista) * 2,), dtype=np.int8)
        })


        self.action_space = spaces.MultiDiscrete([self.ambiente.base, self.ambiente.altura, len(self.ambiente.lista), 4])

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

        self.acao = 0
        self.epoca += 1

        self.ambiente = CSP(dataset='fu', render=False, plot=False)

        self.lista_recompensa = []
        
        pecas = copy.deepcopy(self.ambiente.nova_lista)
        max_len = 4
        for pol in pecas:
            if len(pol) < max_len:
                for i in range(max_len - len(pol)):
                    pol.append((-1,-1))
                
        # Certifique-se de preencher a lista até atingir o tamanho esperado
        while len(pecas) < self.ambiente.tamanho_lista:
            pecas.append([(-1, -1),(-1,-1),(-1,-1),(-1,-1)])  # Adiciona novos polígonos com o número correto de coordenadas
        

        self.pecas = np.array(pecas, dtype=np.int8).flatten()
        self.matriz = np.array(self.ambiente._matriz, dtype=np.int8).flatten()

        self.estado = {
            'cutting_area': self.matriz,
            'pieces': self.pecas
        }

        return self.estado, {}

    def step(self, action):
        self.acao += 1

        n = action[2]
        grau_indice = action[3]
        x = action[0]
        y = action[1]

        tamanho = len(self.ambiente.lista)
        if n < len(self.ambiente.lista):
            self.ambiente.acao(n, x, y, grau_indice)
            self.recompensa = self.ambiente.calcular_recompensa(tamanho, n, self.acao)
            
        else:
            self.recompensa = -2
            #self.ambiente.remover_lista(-1)
            

        self.lista_recompensa.append(self.recompensa)
        
        #if len(self.ambiente.lista) > 0:
        if self.acao < 12:
            self.done = False
        else:
            self.done = True
            print(f"{round((self.ambiente.area_ocupada / (self.ambiente.base * self.ambiente.altura)) * 100, 2)}%, {len(self.ambiente.pecas_posicionadas)}/12, R:{self.lista_recompensa}, SR:{round(sum(self.lista_recompensa),2)}")

        
        

        if self.done:
            self.preenchimentos.append((self.ambiente.area_ocupada/self.ambiente.area) * 100)
            self.pecas_colocadas.append(len(self.ambiente.pecas_posicionadas))
            self.recompensas.append(sum(self.lista_recompensa))

            if (self.epoca % 100) == 0:
                self.calcular_medias()

        pecas = copy.deepcopy(self.ambiente.nova_lista)
        max_len = 4 
        for pol in pecas:
            if len(pol) < max_len:
                for i in range(max_len - len(pol)):
                    pol.append((-1,-1))
            

        while len(pecas) < self.ambiente.tamanho_lista:
            pecas.append([(-1, -1),(-1,-1),(-1,-1),(-1,-1)])

        self.pecas = np.array(pecas, dtype=np.int8).flatten()
        self.matriz = np.array(self.ambiente._matriz, dtype=np.int8).flatten()

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

if __name__ == '__main__':
    env = envCSP
    num_envs = 32
    envs = SubprocVecEnv([make_env(env, i) for i in range(num_envs)])

    log_dir = "./ppo_csp_tensor/"

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    policy_kwargs = dict(
        net_arch=dict(pi=[128, 128], vf=[128, 128]),
    )
    model = PPO('MultiInputPolicy', envs, verbose=1, tensorboard_log=log_dir,
                n_steps=120,  # Ajustado para ser múltiplo de batch_size
                batch_size=240,  # Aumentado para melhor estabilidade
                n_epochs=20,  # Reduzido para treinamento mais rápido
                gamma=0.9,  # Aumentado para considerar recompensas de longo prazo
                learning_rate=linear_schedule(4e-4, 2e-4),
                #learning_rate=3e-4,
                #policy_kwargs=policy_kwargs,
                ent_coef=0.02,
                clip_range=0.3,
                seed = 42)



    model.learn(total_timesteps=1_000_000)
    import os
    from datetime import datetime

    caminho_base = "/home/fsilvestre/Cutting_Stock_Problem"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    caminho_para_salvar = f"{caminho_base}_{timestamp}.zip"
    os.makedirs(os.path.dirname(caminho_para_salvar), exist_ok=True)
    model.save(caminho_para_salvar)


    

  


