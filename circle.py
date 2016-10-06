import math

class Circle:
      #construct a circle object
      def __init__(self, radius = 1):
          self.radius = radius
      def getPerimeter(self):
          return 2*self.radius*math.pi
      def getArea(self):
          return self.radius * self.radius*math.pi
      def setRadius(self, radius):
          self.radius = raduis

def run_t():
    c = Circle(0.5)
    print c.radius
    print c.getPerimeter()
    print c.getArea()

if __name__ == '__main__':
    run_t()
