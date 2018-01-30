# Rory Vander Valk
# Python3

import math

# Alternatives that may be faster or slower in comments
# Most of this should be boring, straightforward

"""
vec3 class
Relatively simple class for handling 3D space coordinates
Access x,y,z as members
Get coordinates as a list with vec3.xyz
TODO? No setter is included for .xyz
"""
class vec3:
  def __init__(self, x = [0.,0.,0.]):
    if not hasattr(x, '__getitem__') or len(x) < 3:
      # Not a vector, ignore it
      self.x,self.y,self.z = 0.,0.,0.
    else:
      self.x = float(x[0]); self.y = float(x[1]); self.z = float(x[2])
  def lengthSqr(self):
    # return self.dot(self)
    return self.x*self.x + self.y*self.y + self.z*self.z
  def length(self):
    # return math.sqrt(self.lengthSqr())
    return math.sqrt(self.x*self.x+self.y*self.y+self.z*self.z)
  def vectorTo(self, target):
    return vec3([
    target.x - self.x,
    target.y - self.y,
    target.z - self.z])
  def distanceSqr(self, other):
    # return self.vectorTo(other).lengthSqr()
    x = self.x - other.x; y = self.y - other.y; z = self.z - other.z
    return x*x+y*y+z*z
  def distance(self, other):
    return math.sqrt(self.distanceSqr(other))
  def __add__(self,other):
    return vec3([self.x+other.x,self.y+other.y,self.z+other.z])
  __radd__ = __add__
  def __sub__(self,other):
    return vec3([self.x-other.x,self.y-other.y,self.z-other.z])
  # TODO This may be wrong...
  # A -= B should be A = A-B so __rsub__ = __sub__ ???
  __rsub__ = __sub__;
  #def __rsub__(self, x):
  #  return vec3([-self.x+x.x,-self.y+x.y,-self.z+x.z])
  # Still needs testing
  # TODO Add a hasattr check here?
  def __mul__(self,other):
    return vec3([self.x*float(other),self.y*float(other),self.z*float(other)])
  __rmul__ = __mul__
  def __div__(self,other):
    return vec3([self.x/float(other),self.y/float(other),self.z/float(other)])
  def __eq__(self,other):
    return self.x == other.x and self.y == other.y and self.z == other.z
  def __str__(self):
    return "{0:< 8.3e}, {1:< 8.3e}, {2:< 8.3e}".format(self.x,self.y,self.z)
  def __neg__(self):
    return vec3([-self.x, -self.y, -self.z])
  def getList(self):
    return [self.x, self.y, self.z]
  def vectorMultiply(self, b):
    return vec3([self.x*b.x, self.y*b.y, self.z*b.z])
  def normalize(self):
    l = self.length()
    # TODO Better if we let a program fail silently?
    #      by returning [1e-30, 0., 0.] when length=0.
    #      Or should we return a value error
    #      Or give them a nice divide by zero error?
    if l == 0.:
      return vec3([0.,0.,0.])
    return vec3([self.x/l, self.y/l, self.z/l])
  def dot(self, b):
    return (self.x*b.x + self.y*b.y + self.z*b.z)
  # self x B -> AxB -> right-handed
  def cross(self, b):
    return vec3([
    self.y*b.z-self.z*b.y,
    -(self.x*b.z-self.z*b.x),
    self.x*b.y-self.y*b.x])
  #
  # [ i  j  k  ]
  # [ ax ay az ]
  # [ bx by bz ]
  #
  # { i*(ay*bz-az*by) - j*(ax*bz-az*bx) + k*(ax*by-bx*ay) }
  # i = (ay*bz-az*by); j = -(ax*bz-az*bx); k = (ax*by-bx*ay)
  def angleDegrees(self, b):
    return math.degrees( math.acos( self.dot(b)/(self.length()*b.length()) ) )
  def __truediv__(self,other):
    return vec3([self.x/float(other),self.y/float(other),self.z/float(other)])
  """Assuming this is in cartesian coordinates, convert to fractional with
  coordinate vectors A, B and C"""
  def getFractional(self, A,B,C):
      la = A.length(); lb = B.length(); lc = C.length()
      #volume = A.Dot(B.Cross(C)) # Triple Scalar Product
      # volume = math.sqrt(1 - cosa**2 - cosb**2 - cosg**2 + 2*ca*cb*cg)
      # 
      #
      cosA = B.normalize().dot(C.normalize())
      cosB = A.normalize().dot(C.normalize())
      cosG = A.normalize().dot(B.normalize())
      sinA = math.sqrt(1. - cosA**2)
      sinB = math.sqrt(1. - cosB**2)
      sinG = math.sqrt(1. - cosG**2)
      asg = la*sinG
      volume = math.sqrt(1. - cosA**2 - cosB**2 - cosG**2 + 2*cosA*cosB*cosG)
      
      r1 = vec3([ 1./la, -cosG/(la*sinG), (cosA*cosG-cosB)/(la*sinG*volume) ])
      r2 = vec3([ 0.,    1./(lb* sinG),         (cosB*cosG-cosA)/(lb*sinG*volume) ])
      r3 = vec3([ 0.,    0.,                          sinG/(lc*volume) ])
      # Cartesian to fractional from C++ code reference
      # frac = transpose(conversion_matrix) * xyz
      # conversion_matrix = {
      #   ( 1./Length(A), -cos(gamma)/(Length(A)*sin(gamma)), (cosa*cosg-cosb)/(la*vol*sing) ),
      #   ( 0., 1./(lb*sing), (cosb * cosg - cosa)/(lb*vol*sing) ),
      #   ( 0., 0., sing/(c*vol) )
      # }
      # Volume = math.sqrt(1 - cosa**2 - cosb**2 - cosg**2 + 2*cosa*cosb*cosg)
      #   ??? Triple product?
      #     = | (axb) dot c |
      #     = l(axb) * lc * costheta
      #
      # Matrix * vector = v.x*column1, v.y*column2, v.z*column3 ;  sum rows
      #  = Dot(vector, column)
      # CosA = 
      #
      return vec3([ self.dot(r1), self.dot(r2), self.dot(r3) ])
  xyz = property(getList, None)