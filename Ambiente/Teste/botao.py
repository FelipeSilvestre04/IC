import turtle

def ponto_dentro_poligono(x, y, poligono):
    n = len(poligono)
    dentro = False

    p1x, p1y = poligono[0]
    for i in range(n+1):
        p2x, p2y = poligono[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        dentro = not dentro
        p1x, p1y = p2x, p2y

    return dentro

class Botao():
    def __init__(self, x, y, largura, altura):
        self.x = x
        self.y = y
        self.largura = largura
        self.altura = altura
        self.tur = turtle.Turtle()
        turtle.tracer(0)
        self.tur.fillcolor("violet")
        self.tur.begin_fill()
        self.tur.penup()
        self.tur.goto(x,y)
        self.tur.pendown()
        self.tur.forward(largura)
        self.tur.right(90)
        self.tur.forward(altura)
        self.tur.right(90)
        self.tur.forward(largura)
        self.tur.right(90)
        self.tur.forward(altura)
        self.tur.right(90)
        self.tur.end_fill()
        turtle.tracer(1)
        self.tur.hideturtle()

    def clicou(self,x,y):
        pol = ((self.x, self.y), (self.x + self.largura, self.y),(self.x+self.largura, self.y-self.altura),(self.x, self.y-self.altura))
        if ponto_dentro_poligono(x,y, pol):
            return True
        else:
            return False

        