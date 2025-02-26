# Importação das bibliotecas

import numpy as np
import random
import os
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torch.autograd as autograd
from torch.autograd import Variable
from collections import deque
from gym.spaces import Discrete

class ConvDQN(nn.Module):
    
    def __init__(self, input_dim, output_dim):
        super(ConvDQN, self).__init__()
        self.input_dim = input_dim
        self.output_dim = output_dim
        
        self.conv = nn.Sequential(
            nn.Conv2d(1, 32, kernel_size=8, stride=4),  # Mude para 1 canal de entrada
            nn.ReLU(),
            nn.Conv2d(32, 64, kernel_size=4, stride=2),
            nn.ReLU(),
            nn.Conv2d(64, 64, kernel_size=3, stride=1),
            nn.ReLU()
        )


        self.fc = nn.Linear(self.feature_size(), self.output_dim)  # Uma única camada linear

    def forward(self, state):
        features = self.conv(state)
        qvals = self.fc(features)
        return qvals

    def feature_size(self):
        return self.conv(autograd.Variable(torch.zeros(1, *self.input_dim))).view(1, -1).size(1)

class BasicBuffer:

    def __init__(self, max_size):
        self.max_size = max_size
        self.buffer = deque(maxlen=max_size)

    def push(self, state, action, reward, next_state, done):
        experience = (state, action, np.array([reward]), next_state, done)
        self.buffer.append(experience)

    def sample(self, batch_size):
        state_batch = []
        action_batch = []
        reward_batch = []
        next_state_batch = []
        done_batch = []

        batch = random.sample(deque(self.capacidade), batch_size)

        for experience in batch:
            state, action, reward, next_state, done = experience
            state_batch.append(state)
            action_batch.append(action)
            reward_batch.append(reward)
            next_state_batch.append(next_state)
            done_batch.append(done)

        return (state_batch, action_batch, reward_batch, next_state_batch, done_batch)

    def __len__(self):
        return len(self.buffer)



class DQNAgent:

    def __init__(self, env, learning_rate=3e-4, gamma=0.99, buffer_size=10000):
        self.env = env
        self.learning_rate = learning_rate
        self.gamma = gamma
        self.replay_buffer = BasicBuffer(max_size=buffer_size)
	
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model1 = ConvDQN((1, 200, 300), 2).to(self.device)  # entrada é uma matriz 1x200x300, saída é um par ordenado
        self.model2 = ConvDQN((1, 200, 300), 2).to(self.device)  # entrada é uma matriz 1x200x300, saída é um par ordenado

        self.optimizer1 = torch.optim.Adam(self.model1.parameters())
        self.optimizer2 = torch.optim.Adam(self.model2.parameters())
        
    def get_action(self, state, eps=0.20):
        state = torch.FloatTensor(state).float().unsqueeze(0).to(self.device)
        qvals = self.model1.forward(state)  # Use self.model1 aqui
        action = np.unravel_index(np.argmax(qvals.cpu().detach().numpy()), qvals.shape)

        if(np.random.randn() < eps):
            return (self.env.action_space.sample(), self.env.action_space.sample())

        return action


    def compute_loss(self, batch):     
        states, actions, rewards, next_states, dones = batch
        state = torch.FloatTensor(state.to_numpy()).float().unsqueeze(0).to(self.device)
        actions = torch.LongTensor(actions).to(self.device)
        rewards = torch.FloatTensor(rewards).to(self.device)
        next_states = torch.FloatTensor(next_states).to(self.device)
        dones = torch.FloatTensor(dones)

        # resize tensors
        actions = actions.view(actions.size(0), 1)
        dones = dones.view(dones.size(0), 1)

        # compute loss
        curr_Q1 = self.model1.forward(states).gather(1, actions)
        curr_Q2 = self.model2.forward(states).gather(1, actions)

        next_Q1 = self.model1.forward(next_states)
        next_Q2 = self.model2.forward(next_states)
        next_Q = torch.min(
            torch.max(self.model1.forward(next_states), 1)[0],
            torch.max(self.model2.forward(next_states), 1)[0]
        )
        next_Q = next_Q.view(next_Q.size(0), 1)
        expected_Q = rewards + (1 - dones) * self.gamma * next_Q

        loss1 = F.mse_loss(curr_Q1, expected_Q.detach())
        loss2 = F.mse_loss(curr_Q2, expected_Q.detach())

        return loss1, loss2
        
    def update(self, batch_size):
        batch = self.replay_buffer.sample(batch_size)
        loss1, loss2 = self.compute_loss(batch)

        self.optimizer1.zero_grad()
        loss1.backward()
        self.optimizer1.step()

        self.optimizer2.zero_grad()
        loss2.backward()
        self.optimizer2.step()