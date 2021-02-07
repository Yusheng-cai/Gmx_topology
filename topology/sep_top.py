import sys
from .Molecule_properties import atom,bond,angle,dihedral


class topology:
    def __init__(self,file_name):
        self.file_name = file_name
        self.info = self.read_dat()

        self.atoms = self.get_atoms(self.info)
        self.bonds_list,self.unique_bonds = self.get_bonds(self.info,self.atoms)
        self.angles_list,self.unique_angles = self.get_angles(self.info,self.atoms)
        self.dihedrals_list,self.unique_dihedrals = self.get_dihedrals(self.info,self.atoms)

    def read_dat(self):
        """
        Function that reads the data from the gromacs topology file

        Return:
        ------
        info(dict): dictionary that holds all the lines corresponding to each directive
        """
        f = open(self.file_name)
        lines = f.readlines()
        f.close()
        
        # get rid of all the empty lines as well as \n in the data
        lines = [l for l in lines if l !="\n"]
        lines = [l.rstrip("\n").lstrip() for l in lines] 

        # Create an empty list that will hold all the directive names as well as their index in the data read from file
        directive_idx = []
        directives = []

        for index,line in enumerate(lines):
            # The line is a directive if it starts with "[" e.g. [ dihedral ]
            if line.startswith("["):
                directives.append(line)
                directive_idx.append(index)

        # Create the dictionary and put the data of each directive in it 
        info = {}
        for i in range(len(directives)):
            d = directives[i]
            idx = directive_idx[i]
            if i<len(directives)-1:
                next_idx = directive_idx[i+1]
                info[d] = lines[idx:next_idx]
            else:
                info[d] = lines[idx:]
        
        return info
    
    def get_atoms(self,info):
        atoms = info["[ atoms ]"]
        atoms = [l for l in atoms[1:] if not l.startswith(";")]
        atoms_dic = {}

        # The first two lines are comments
        ix = 1
        for l in atoms:
            atoms_dic[ix] = atom(l) 
            ix += 1

        return atoms_dic

    def get_bonds(self,info,atominfo):
        bonds = info["[ bonds ]"]
        bonds = [l for l in bonds[1:] if not l.startswith(";")]
        bonds_list = []
        unique_bonds = []

        # The first two lines are comments
        for l in bonds:
            b = bond(l,atominfo)
            if b not in bonds_list:
                unique_bonds.append(b)
            bonds_list.append(b) 

        return bonds_list,unique_bonds

    def get_angles(self,info,atominfo):
        angles = info["[ angles ]"]
        angles = [l for l in angles[1:] if not l.startswith(";")]
        angles_list = []
        unique_angles = []

        # The first two lines are comments
        for l in angles:
            a = angle(l,atominfo)
            if a not in angles_list:
                unique_angles.append(a)
            angles_list.append(a) 

        return angles_list,unique_angles

    def get_dihedrals(self,info,atominfo):
        dihedrals = info["[ dihedrals ]"]
        dihedrals = [l for l in dihedrals[1:] if not l.startswith(";")]
        dihedrals_list = []
        unique_dihedrals = []

        # The first two lines are comments
        for l in dihedrals:
            d = dihedral(l,atominfo)
            if d not in dihedrals_list:
                unique_dihedrals.append(d)
            dihedrals_list.append(d) 

        return dihedrals_list,unique_dihedrals
