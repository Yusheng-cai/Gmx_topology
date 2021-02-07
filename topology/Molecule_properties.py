class atom:
    """
    Class that describes how an atom is represented in the gromacs topology file
    """
    def __init__(self,line):
        # nr type  resnr residue  atom   cgnr    charge       mass
        # get rid of all the comments in the line that holds the atom information
        if ";" in line:
            line = line.split(";")[0].rstrip()

        line = line.split()
        self.nr = line[0]
        self.type = line[1]
        self.resnr = line[2]
        self.residue = line[3]
        self.atom = line[4]
        self.cgnr = line[5]
        self.charge = line[6]
        self.mass = line[7]

    def __eq__(self,other):
        if (self.type == other.type) & (self.mass == other.mass):
            return True
        else:
            return False

    def __str__(self):
        return "Atom of type {} with residue {} of mass {} and charge {}".format(self.type,\
                self.residue,self.mass,self.charge)

    def __repr__(self):
        return self.__str__()

class bond:
    """
    Class that represents how a bond is reprenseted in the gromacs topology file
    """
    def __init__(self,line,atominfo):
        # ai aj funct c0 c1 c2 c3
        if ";" in line:
            line = line.split(";")[0].rstrip()
        line = line.split()
        N = len(line)
        Nparams = N - 3
        
        self.atom1 = atominfo[int(line[0])]
        self.atom2 = atominfo[int(line[1])]
        self.funct = int(line[2])
        self.params = []
        self.Nparams = Nparams

        for i in range(3,N):
            self.params.append(line[i])

        if self.funct == 1 or self.funct == 2:
            self.str = "{0:>6}{1:>6}{2:>6}{3:>9.3f}{4:>12.3f}\n".format(self.atom1.type,\
                    self.atom2.type,\
                    self.funct,float(self.params[0]),float(self.params[1]))
       
    def __eq__(self,other):
        if (self.atom1==other.atom1) and (self.atom2==other.atom2) and \
                (self.funct == other.funct) and (self.params==other.params):
            return True
        elif (self.atom1==other.atom2) and (self.atom2==other.atom1) and \
                (self.funct == other.funct) and (self.params==other.params): 
            return True
        else:
            return False

    def __str__(self):
        return "Bond between {}-{} with funct {} and parameters".format(self.atom1.type,self.atom2.type,self.funct,self.params)
    
    def __repr__(self):
        return self.__str__()

class angle:
    """
    Class that representes an angle based on the representation of angles in GROMACS
    
    Args:
    ----
    line(str): A str that represents a line in the topology file of GROMACS that represents an angle 
    atominfo(dict): A dictionary that holds all the atoms in the molecule which is passed in as a dictionary 
    """
    def __init__(self,line,atominfo):
        # ai aj  ak funct c0 c1 c2  c3
        # delete comments in a line
        if ";" in line:
            line = line.split(";")[0].rstrip()

        line = line.split()
        N = len(line)
        Nparams = N - 4
        
        self.atom1 = atominfo[int(line[0])]
        self.atom2 = atominfo[int(line[1])]
        self.atom3 = atominfo[int(line[2])]
        self.funct = int(line[3])
        self.params = []
        self.Nparams = Nparams

        for i in range(4,N):
            self.params.append(line[i])

        if self.funct == 1 or self.funct == 2:
            self.str = "{0:>6}{1:>6}{2:>6}{3:>6}{4:>9.3f}{5:>12.3f}\n".format(self.atom1.type,\
                    self.atom2.type,self.atom3.type,\
                    self.funct,float(self.params[0]),float(self.params[1]))

    def __eq__(self,other):
        if (self.atom1==other.atom1) and (self.atom2==other.atom2) and \
                (self.atom3==other.atom3) and (self.funct == other.funct) and (self.params==other.params):
            return True
        elif (self.atom1==other.atom3) and (self.atom2==other.atom2) and \
                (self.atom3==other.atom1) and (self.funct == other.funct) and (self.params==other.params): 
            return True
        else:
            return False

    def __str__(self):
        return "Angle between {}-{}-{} with funct {} and params {}".format(self.atom1.type,self.atom2.type,\
                self.atom3.type,self.funct, self.params)
    
    def __repr__(self):
        return self.__str__()


class dihedral:
    """
    Class that represents the dihedrals in GROMACS topology files
    """
    def __init__(self,line,atominfo):
        # ai aj ak al funct  c0  c1  c2  c3  c4  c5
        if ";" in line:
            line = line.split(";")[0].rstrip()
        line = line.split()
        N = len(line)
        
        self.atom1 = atominfo[int(line[0])]
        self.atom2 = atominfo[int(line[1])]
        self.atom3 = atominfo[int(line[2])]
        self.atom4 = atominfo[int(line[3])]
        self.funct = int(line[4])
        self.params = []
        for i in range(5,N):
            self.params.append(line[i])

        if self.funct == 1 or self.funct==4 or self.funct==9:
            # funct = 1 is section 4.2.13 of gromacs manual 5.1.4 which has parameters phis,kphi,multiplicity
            # funct = 4 is section 4.2.12 of gromacs manual 5.1.4 which has parameters phis,kphi and multiplicity
            self.str = "{0:>6}{1:>6}{2:>6}{3:>6}{4:>6}{5:>9.3f}{6:>9.3f}{7:>9d}\n".format(self.atom1.type,\
                    self.atom2.type,self.atom3.type,self.atom4.type,\
                    self.funct,float(self.params[0]),float(self.params[1]),int(self.params[2]))
        
        if self.funct == 3:
            # funct = 3 is RB dihedral in section 4.2.13 of gromacs manual 5.1.4 which has parameters C0, C1, C2, C3, C4, C5
            self.str = "{0:>6}{1:>6}{2:>6}{3:>6}{4:>6}{5:>9.3}{6:>9.3f}{6:>9.3f}{7:9.3f}{8:9.3f}\n".format(self.atom1.type,\
                    self.atom2.type,self.atom3.type,self.atom4.type,\
                    self.funct,float(self.params[0]),float(self.params[1]),float(self.params[2]),\
                    float(self.params[3]),float(self.params[4]))

    def __eq__(self,other):
        if (self.atom2==other.atom2) and (self.atom3==other.atom3) and \
                (self.atom1==other.atom1) and (self.atom4==self.atom4) and (self.funct == other.funct) and (self.params==other.params):
            return True
        elif (self.atom2==other.atom2) and \
                (self.atom3==other.atom3) and (self.funct == other.funct) and (self.params==other.params): 
            return True
        else:
            return False

    def __str__(self):
       return "Dihedral between {}-{}-{}-{} with funct {} with params {}".format(self.atom1.type,self.atom2.type,\
                self.atom3.type,self.atom4.type,self.funct,self.params)
    
    def __repr__(self):
        return self.__str__()
