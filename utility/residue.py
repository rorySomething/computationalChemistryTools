
# TODO complete implementation


class Residue:
  def __init__(self):
    self.name = ""
    self.charge = 0.
    self.atoms = dict() # Name:name,type,charge
    self.bonds = []
    self.impropers = []
    self.delete = [] # For patches (PRES)
    self.donors = [] # Pairs, H N
    self.acceptors = [] # Pairs, O C
    self.internalCoordinates = [] # I, J, K, L, r_ij, t_ijk, p_ijkl, t_jkl, r_kl
  def CheckCharge(self, verbose = False):
    charge = 0.
    for i in self.atoms.values():
      charge += i.charge
    dQ = abs(self.charge - charge)
    if dQ > 1e-4:
      print("Residue",self.name,"has bad topology data.")
      print("\tCharge is off by", dQ)
