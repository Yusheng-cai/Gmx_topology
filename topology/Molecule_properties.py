class atom_type:
    """
    Class that describes the atomtype representation in the GROMACS topology file

    Args:
    ----
    line(str): A string which is the line in the GROMACS topology file that holds information about an atom
    """
    def __init__(self,line):
        # name  at.num  mass charge ptype  sigma epsilon
        if ";" in line:
            line = line.split(";")[0].rstrip()

        line = line.split()
        self.name = line[0]
        self.atnum = int(line[1])
        self.mass = float(line[2])
        self.charge = float(line[3])
        self.ptype = line[4]
        self.sigma = float(line[5])
        self.epsilon = float(line[6])

        self.strff = "{0:>6}{1:>6d}{2:>9.3f}{3:>9.3f}{4:>6}{5:>12.5f}{6:>12.5f}\n".format(self.name,self.atnum,\
                self.mass,self.charge,self.ptype,self.sigma,self.epsilon)

    def __str__(self):
        return "Atom type {} with mass {}.Sigma:{}, Epsilon:{}".format(self.name, self.mass, self.sigma,self.epsilon)
    
    def __repr__(self):
        return self.__str__()

class atom:
    """
    Class that describes how an atom is represented in the gromacs topology file. This only goes into the file for the molecular definition of a molecule

    Args:
    ----
    line(str): The line from the GROMACS topology file that holds information about the atom
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
        self.cgnr = int(line[5])
        self.charge = float(line[6])
        self.mass = float(line[7])

        self.strmol = "{0:>6}{1:>6}{2:>6}{3:>7}{4:>7}{5:>6d}{6:>12.5f}{7:>12.5f}\n".format(self.nr,self.type,self.resnr,\
                self.residue,self.atom,self.cgnr,self.charge,self.mass)

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
    def __init__(self,line,atominfo,mol_params=False):
        # ai aj funct c0 c1 c2 c3
        if ";" in line:
            line = line.split(";")[0].rstrip()
        line = line.split()
        N = len(line)
        Nparams = N - 3
        
        self.atom1 = atominfo[int(line[0])]
        self.atnum1 = int(line[0])
        self.atom2 = atominfo[int(line[1])]
        self.atnum2 = int(line[1])

        self.funct = int(line[2])
        self.params = []
        self.Nparams = Nparams
        self.mol_params = mol_params

        for i in range(3,N):
            self.params.append(line[i])

        if self.funct == 1 or self.funct == 2:
            self.strff = "{0:>6}{1:>6}{2:>6}{3:>9.3f}{4:>12.3f}\n".format(self.atom1.type,\
                    self.atom2.type,\
                    self.funct,float(self.params[0]),float(self.params[1]))
        else:
            raise NotImplementedError("The given function {} is not implemented yet".format(self.funct))

            
        if self.mol_params:
            if self.funct == 1 or self.funct == 2:
                self.strmol = "{0:>6}{1:>6}{2:>6}{3:>9.3f}{4:>12.3f}\t;{5:>4}{6:>4}\n".format(self.atnum1,\
                    self.atnum2,self.funct,\
                    float(self.params[0]),float(self.params[1]),\
                    self.atom1.type,self.atom2.type)
            else:
                raise NotImplementedError("The given function {} is not implemented yet".format(self.funct))
        else:
            self.strmol = "{0:>6}{1:>6}\t;{2:>4}{3:>4}\n".format(self.atnum1,\
                    self.atnum2,self.atom1.type,self.atom2.type)
 
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
        return "Bond between {}-{} with funct {} and parameters {}".format(self.atom1.type,self.atom2.type,self.funct,self.params)
    
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
    def __init__(self,line,atominfo,mol_params=False):
        # ai aj  ak funct c0 c1 c2  c3
        # delete comments in a line
        if ";" in line:
            line = line.split(";")[0].rstrip()

        line = line.split()
        N = len(line)
        Nparams = N - 4
        
        self.atom1 = atominfo[int(line[0])]
        self.atnum1 = int(line[0])

        self.atom2 = atominfo[int(line[1])]
        self.atnum2 = int(line[1])

        self.atom3 = atominfo[int(line[2])]
        self.atnum3 = int(line[2])

        self.funct = int(line[3])
        self.params = []
        self.Nparams = Nparams
        self.mol_params = mol_params

        for i in range(4,N):
            self.params.append(line[i])

        if self.funct == 1 or self.funct == 2:
            self.strff = "{0:>6}{1:>6}{2:>6}{3:>6}{4:>9.3f}{5:>12.3f}\n".format(self.atom1.type,\
                    self.atom2.type,self.atom3.type,\
                    self.funct,float(self.params[0]),float(self.params[1]))
        else:
            raise NotImplementedError("The given function {} is not implemented yet".format(self.funct))

        # See if parameter needs to be written in molecule file
        if self.mol_params:
            if self.funct == 1 or self.funct == 2:
                self.strmol = "{0:>6}{1:>6}{2:>6}{3:>6}{4:>9.3f}{5:>12.3f}\t;{6:>4}{7:>4}{8:>4}\n".format(self.atnum1,\
                    self.atnum2,self.atnum3,\
                    self.funct,float(self.params[0]),float(self.params[1]),self.atom1.type,self.atom2.type,self.atom3.type)
            else:
                raise NotImplementedError("The given function {} is not implemented yet".format(self.funct))
        else:
            self.strmol = "{0:>6}{1:>6}{2:>6}\t;{3:>4}{4:>4}{5:>4}\n".format(self.atnum1,\
                    self.atnum2,self.atnum3,self.atom1.type,self.atom2.type,self.atom3.type)

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
    def __init__(self,line,atominfo,mol_params=False):
        # ai aj ak al funct  c0  c1  c2  c3  c4  c5
        if ";" in line:
            line = line.split(";")[0].rstrip()
        line = line.split()
        N = len(line)
        
        self.atom1 = atominfo[int(line[0])]
        self.atnum1 = int(line[0])

        self.atom2 = atominfo[int(line[1])]
        self.atnum2 = int(line[1])

        self.atom3 = atominfo[int(line[2])]
        self.atnum3 = int(line[2])

        self.atom4 = atominfo[int(line[3])]
        self.atnum4 = int(line[3])

        self.funct = int(line[4])
        self.params = []
        self.mol_params = mol_params

        for i in range(5,N):
            self.params.append(line[i])

        if self.funct == 1 or self.funct==4 or self.funct==9:
            # funct = 1 is section 4.2.13 of gromacs manual 5.1.4 which has parameters phis,kphi,multiplicity
            # funct = 4 is section 4.2.12 of gromacs manual 5.1.4 which has parameters phis,kphi and multiplicity
            self.strff = "{0:>6}{1:>6}{2:>6}{3:>6}{4:>6}{5:>9.3f}{6:>9.3f}{7:>9d}\n".format(self.atom1.type,\
                    self.atom2.type,self.atom3.type,self.atom4.type,\
                    self.funct,float(self.params[0]),float(self.params[1]),int(self.params[2])) 
        elif self.funct == 3:
            # funct = 3 is RB dihedral in section 4.2.13 of gromacs manual 5.1.4 which has parameters C0, C1, C2, C3, C4, C5
            self.strff = "{0:>6}{1:>6}{2:>6}{3:>6}{4:>6}{5:>9.3}{6:>9.3f}{6:>9.3f}{7:9.3f}{8:9.3f}{9:9.3f}\n".format(self.atom1.type,\
                    self.atom2.type,self.atom3.type,self.atom4.type,\
                    self.funct,float(self.params[0]),float(self.params[1]),float(self.params[2]),\
                    float(self.params[3]),float(self.params[4]),float(self.params[5]))
        else:
            raise NotImplementedError("The given function {} is not implemented yet".format(self.funct))


        # If parameters are included in molecule
        if self.mol_params:
            if self.funct == 3:
                self.strmol = "{0:>6}{1:>6}{2:>6}{3:>6}{4:>6}{5:>12.6f}{6:>12.6f}{7:>12.6f}{8:12.6f}{9:12.6f}{10:12.6f}\t;{11:>4}{12:>4}{13:>4}{14:>4}\n".format(self.atnum1,\
                        self.atnum2,self.atnum3,self.atnum4,\
                        self.funct,float(self.params[0]),float(self.params[1]),float(self.params[2]),\
                        float(self.params[3]),float(self.params[4]),float(self.params[5]),self.atom1.type,self.atom2.type,self.atom3.type,self.atom4.type)
            elif self.funct == 1 or self.funct == 4 or self.funct == 9:
                self.strmol = "{0:>6}{1:>6}{2:>6}{3:>6}{4:>6}{5:>9.3f}{6:>9.3f}{7:>9d}\t;{8:>4}{9:>4}{10:>4}{11:>4}\n".format(self.atnum1,\
                    self.atnum2,self.atnum3,self.atnum4,\
                    self.funct,float(self.params[0]),float(self.params[1]),int(self.params[2]),self.atom1.type,self.atom2.type,\
                    self.atom3.type,self.atom4.type)
            else:
                raise NotImplementedError("The given function {} is not implemented yet".format(self.funct))
        else:
            self.strmol = "{0:>6}{1:>6}{2:>6}{3:>6}\t;{4:>4}{5:>4}{6:>4}{7:>4}\n".format(self.atnum1,\
                    self.atnum2,self.atnum3,self.atnum4,self.atom1.type,self.atom2.type,self.atom3.type,\
                    self.atom4.type)

    def general_dihedral(self):
        if self.funct == 1 or self.funct==4 or self.funct==9:
            # funct = 1 is section 4.2.13 of gromacs manual 5.1.4 which has parameters phis,kphi,multiplicity
            # funct = 4 is section 4.2.12 of gromacs manual 5.1.4 which has parameters phis,kphi and multiplicity
            self.strff = "{0:>6}{1:>6}{2:>6}{3:>6}{4:>6}{5:>9.3f}{6:>9.3f}{7:>9d}\n".format("X",\
                    self.atom2.type,self.atom3.type,"X",\
                    self.funct,float(self.params[0]),float(self.params[1]),int(self.params[2]))
     
    def __eq__(self,other):
        # Only return true if all the atom numbers match (So it is the same dihedral)
        if (self.atnum1==other.atnum1) and (self.atnum2==other.atnum2) and \
                (self.atnum3==other.atnum3) and (self.atnum4==other.atnum4):
            return True
        else:
            return False

    def append(self,other):
        """
        Function that appends another dihedral to the current dihedral, most likely used scenario is when there is one dihedral that has multiple multiplicities 

        Args:
        ----
        other(dihedral): A dihedral object

        Return:
        ------
        strmol, strff and params are extended by the other dihedral object
        """
        # check if passed in is a list of dihedral objects
        if (self.strmol != other.strmol) and (self.strff != other.strff) and (self.params != other.params):
            self.strmol += other.strmol
            self.strff += other.strff
            self.params += other.params

    def compare(self,other):
        """
        Function that returns True if all the atom types as well as all the parameters matches with each other
        """
        if (self.atom1.type==other.atom1.type) and (self.atom2.type==other.atom2.type) and (self.atom3.type==other.atom3.type) and (self.atom4.type==self.atom4.type) and (self.funct==other.funct) and (self.params == other.params):
            return True
        elif (self.atom1.type==other.atom4.type) and (self.atom2.type==other.atom3.type) and (self.atom3.type==other.atom2.type) and (self.atom4.type==self.atom1.type) and (self.funct==other.funct) and (self.params == other.params):
            return True
        else:
            return False   

    def __str__(self):
       return "Dihedral between {}-{}-{}-{} with funct {} with params {}".format(self.atom1.type,self.atom2.type,\
                self.atom3.type,self.atom4.type,self.funct,self.params)
    
    def __repr__(self):
        return self.__str__()
