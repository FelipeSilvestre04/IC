import turtle 
from geradores import Area, Gerador, posicionar, on_click
from matriz import Matriz
from botao import Botao
import pyautogui
import math
import time

base = 600
altura = 400
x = -350
y = 300
greedys = 6
turtle.Turtle()
turtle.tracer(0)
turtle.speed(200)
tecido = Area(base, altura, greedys)
matriz2 = Matriz(base, altura, x, y)
tecido.area(x, y)
tecido.greedy(x, y)
isdrawing = False
posicionar(54, "M", True)
Bot = Botao(350, 300, 200, 150)
turtle.tracer(1)



text_turtle = turtle.Turtle()
fecho_turtle = turtle.Turtle()
def on_click2(x,y):
    on_click(x,y,matriz2,text_turtle, fecho_turtle, base, altura, greedys, Bot)




def raio_circunscrito(n, s):
    r = s / (2 * math.sin(math.pi / n))
    return r

def draw_piece(loop):
    # Pega a posição atual do mouse
    x, y = pyautogui.position()
    if len(Gerador._geracao) > 0:
        shape_name = Gerador._geracao[-1]
    else:
        loop = False
        return loop
    n = int(shape_name.split(',')[0])  
    s = int(shape_name.split(',')[2])  
    turtle.tracer(0)
    turtle.clear()
    turtle.penup()
    turtle.goto(x-963,-y+528)
    turtle.pendown()
    turtle.right(90)  #
    turtle.right(180/n)  
    turtle.forward(raio_circunscrito(n, s))  
    turtle.left(180/n)  
    turtle.left(90)
    turtle.color('green')
    turtle.begin_fill()
    for _ in range(int(shape_name.split(',')[0])):
            turtle.forward(int(shape_name.split(',')[2]))
            turtle.left(360 / int(shape_name.split(',')[0]))
    turtle.end_fill()
    turtle.hideturtle()
    turtle.tracer(1)

t = turtle.Turtle()
turtle.tracer(0)
t.color("black")
t.penup()
t.goto(-925,250)
t.pendown()

def tempo(time):
    # Limpe a escrita anterior antes de escrever o novo tempo
    t.clear()
    t.write(f"{int(time/60)} min {round(time%60,0)} s",font=("Arial", 16, "normal"))
    t.hideturtle()
    turtle.update()

# Loop principal
loop = True
inicio = time.time()
tempo_anterior = 0
while loop:
    fim = time.time()
    tempo_execucao = fim - inicio 
    # Atualize o tempo a cada segundo
    if int(tempo_execucao) != int(tempo_anterior):
        tempo(tempo_execucao)
        tempo_anterior = int(tempo_execucao)
    draw_piece(loop)
    turtle.onscreenclick(on_click2)









