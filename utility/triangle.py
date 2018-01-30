from .vec3 import vec3

class triangle:
  def __init__(self, a, b, c):
    self.a = a
    self.b = b
    self.c = c
  def center(self, mesh):
    A = mesh.points[a]
    B = mesh.points[b]
    C = mesh.points[c]
    return (A+B+C)/3
  # Centroid?
  def center2(self, mesh):
    A = mesh.points[a]
    B = mesh.points[b]
    C = mesh.points[c]
    AB = A+((B-A)*.5) # midpoint of AB
    return AB + ((C-AB)*.5) # midpoint of AB midpoint and C
  def normal(self, mesh):
    A = mesh.points[a]
    B = mesh.points[b]
    C = mesh.points[c]
    AB = B-A
    AC = C-A
    # CCW
    return AB.cross(AC).normalize()
  def invert(self):
    self.a,self.c = self.c,self.a
  def isFront(self, mesh, point):
    n = self.normal(mesh)
    c = self.center(mesh)
    if n.dot(point-c) > 0.:
      return True
    else:
      return False