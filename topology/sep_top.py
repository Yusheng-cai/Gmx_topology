from .Molecule_properties import atom_type,atom,bond,angle,dihedral


class topology:
    def __init__(self,file_name,mol_params=False):
        self.file_name = file_name
        self.info = self.read_dat()
        self.mol_params = mol_params

        self.atoms = self.get_atoms(self.info)
        self.bonds_list,self.unique_bonds = self.get_bonds(self.info,self.atoms,self.mol_params)
        self.angles_list,self.unique_angles = self.get_angles(self.info,self.atoms,self.mol_params)
        self.dihedrals_list,self.unique_dihedrals = self.get_dihedrals(self.info,self.atoms,self.mol_params)

        
        if "[ atomtypes ]" in self.info:
            self.atom_types_list = self.get_atomtypes(self.info)

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
        """
        Function that obtains all the atoms in the molecule

        Args:
        ----
        info(dict): The key of info is the directives (e.g. [ dihedral ]) and the value are the lines belongs to that directive in the topology file

        Return:
        ------
        atom_dic(dict): An dictionary that contains all the atoms information (type, charge, mass etc.)
        """
        atoms = info["[ atoms ]"]
        atoms = [l for l in atoms[1:] if not l.startswith(";")]
        atoms_dic = {}

        # The first two lines are comments
        ix = 1
        charge_total = 0
        for l in atoms:
            a = atom(l)
            atoms_dic[ix] = a 
            charge_total += a.charge
            ix += 1

        print("Total charge of the molecule is {}".format(charge_total))
        print("Modifying...")
        
        charge_total = 0
        if charge_total != 0:
            for a in atoms_dic:
                atoms_dic[a].charge += charge_total/len(atoms_dic)
                charge_total += atoms_dic[a].charge

        print("Modified charge is {}".format(charge_total))
            


        return atoms_dic

    def get_bonds(self,info,atominfo,mol_params=False):
        """
        Function that obtains all the bonds in the molecule as well as all the unique bonds

        Args:
        ----
        info(dict): The key of info is the directives (e.g. [ dihedral ]) and the value are the lines belongs to that directive in the topology file
        atominfo(dict): The dictionary of all atoms info, the key is the atom number in the molecule and the value is the atom object 

        Return:
        ------
        1. bonds_list(list)= A list of bond objects that contains ALL the bonds in the molecule
        2. unique_bonds(list)=A list of bond objects that contains UNIQUE bonds in the molecule
        """
        bonds = info["[ bonds ]"]
        bonds = [l for l in bonds[1:] if not l.startswith(";")]
        bonds_list = []
        unique_bonds = []

        # The first two lines are comments
        for l in bonds:
            b = bond(l,atominfo,mol_params)
            if b not in bonds_list:
                unique_bonds.append(b)
            bonds_list.append(b) 

        return bonds_list,unique_bonds

    def get_angles(self,info,atominfo,mol_params=False):
        """
        Function that obtains all the angles in the molecule as well as all the unique angles 

        Args:
        ----
        info(dict): The key of info is the directives (e.g. [ dihedral ]) and the value are the lines belongs to that directive in the topology file
        atominfo(dict): The dictionary of all atoms info, the key is the atom number in the molecule and the value is the atom object 

        Return:
        ------
        1. angles_list(list)=A list of angle objects that includes all the angles present in the molecule
        2. unique_angles(list)=A list of angle objects that includes all the unique angles present in the molecule
        """
        angles = info["[ angles ]"]
        angles = [l for l in angles[1:] if not l.startswith(";")]
        angles_list = []
        unique_angles = []

        # The first two lines are comments
        for l in angles:
            a = angle(l,atominfo,mol_params)
            if a not in angles_list:
                unique_angles.append(a)
            angles_list.append(a) 

        return angles_list,unique_angles

    def get_dihedrals(self,info,atominfo,mol_params=False):
        """
        Function that obtains all the dihedrals in the molecule as well as all the unique dihedrals

        Args:
        ----
        info(dict): The key of info is the directives (e.g. [ dihedral ]) and the value are the lines belongs to that directive in the topology file
        atominfo(dict): The dictionary of all atoms info, the key is the atom number in the molecule and the value is the atom object 

        Return:
        ------
        1. dihedrals_list(list) = A list of dihedral objects that are present in the GROMACS topology file.
        2. unique_dihedrals(list) = A list of unique dihedral objects that are present in the GROMACS topology file.
        """
        dihedrals = info["[ dihedrals ]"]
        dihedrals = [l for l in dihedrals[1:] if not l.startswith(";")]
        all_dihedral_fcn = []
        dihedrals_list = []
        unique_dihedrals = []

        # The first two lines are comments
        for l in dihedrals:
            d = dihedral(l,atominfo,mol_params)
            all_dihedral_fcn.append(d)
            # dihedrals list only contains all the unique dihedrals as in the atom numbering
            if d not in dihedrals_list:
                dihedrals_list.append(d)
        
        # Append all the functions to each dihedral
        for d1 in dihedrals_list:
            for d2 in all_dihedral_fcn:
                if d2 == d1:
                    d1.append(d2)

        unique_dihedrals.append(dihedrals_list[0])
        for d1 in dihedrals_list:
            for d2 in unique_dihedrals:
                unique = True
                if d1.compare(d2):
                    unique = False
            if unique:
                unique_dihedrals.append(d1)

        return dihedrals_list,unique_dihedrals


    def get_atomtypes(self,info):
        """
        Function that obtains all the atom types in the GROMACS topology (the Lennard Jones parameters)
        
        Args:
        -----
        info(dict): The key of info is the directives (e.g. [ dihedral ]) and the value are the lines belongs to that directive in the topology file
        atominfo(dict): The dictionary of all atoms info, the key is the atom number in the molecule and the value is the atom object 

        Return:
        ------
        atom_types_list(list): A list of all the atom_type objects 
        """
        # name  at.num  mass charge ptype sigma epsilon
        atom_types = info["[ atomtypes ]"]
        atom_types = [l for l in atom_types[1:] if not l.startswith(";")]
        atom_types_list = []

        for l in atom_types:
            at = atom_type(l)
            atom_types_list.append(at)

        return atom_types_list

    def substitute_dihedrals(self,str_):
        """
        A function that substitutes dihedrals from 
        The dihedral must be passed in the form:

        ai aj ak al funct params

        The number of params depends on the funct
        """
        str_ = str_.rstrip("\n")
        split = str_.split()
        type1 = split[0]
        type2 = split[1]
        type3 = split[2] 
        type4 = split[3]
        funct = int(split[4])
        params = []
        for i in range(5,len(split)):
            params.append(split[i])


        for i in range(len(self.dihedrals_list)):
            d = self.dihedrals_list[i]
            if (d.atom1.type == type1) & (d.atom2.type == type2) &\
                    (d.atom3.type == type3) & (d.atom4.type == type4):
                        if funct == 3:
                            s = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(d.atnum1,d.atnum2,d.atnum3,d.atnum4,\
                                    funct,params[0],params[1],params[2],params[3],params[4],params[5])
                            self.dihedrals_list[i] = dihedral(s,self.atoms,mol_params=True)
            elif (d.atom1.type == type4) & (d.atom2.type == type3) &\
                    (d.atom3.type == type2) & (d.atom4.type == type1):
                        if funct == 3:
                            s = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(d.atnum1,d.atnum2,d.atnum3,d.atnum4,\
                                    funct,params[0],params[1],params[2],params[3],params[4],params[5])
                            self.dihedrals_list[i] = dihedral(s,self.atoms,mol_params=True)

        for i in range(len(self.unique_dihedrals)):
            d = self.unique_dihedrals[i]
            if (d.atom1.type == type1) & (d.atom2.type == type2) &\
                    (d.atom3.type == type3) & (d.atom4.type == type4):
                        if funct == 3:
                            s = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(d.atnum1,d.atnum2,d.atnum3,d.atnum4,\
                                    funct,params[0],params[1],params[2],params[3],params[4],params[5])
                            self.unique_dihedrals[i] = dihedral(s,self.atoms,mol_params=True)
            elif (d.atom1.type == type4) & (d.atom2.type == type3) &\
                    (d.atom3.type == type2) & (d.atom4.type == type1):
                        if funct == 3:
                            s = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(d.atnum1,d.atnum2,d.atnum3,d.atnum4,\
                                    funct,params[0],params[1],params[2],params[3],params[4],params[5])
                            self.unique_dihedrals[i] = dihedral(s,self.atoms,mol_params=True)




    def write_ff(self,o_name):
        """
        Function that writes the force field to a file

        Args:
        ----
        o_name(str): output name of the file

        Return:
        ------
        force field file
        """
        unique_bonds = self.unique_bonds
        unique_angles = self.unique_angles
        unique_dihedrals = self.unique_dihedrals

            
        f = open(o_name,"w")
        
        if "[ atomtypes ]" in self.info:
            atom_types_list = self.atom_types_list
            f.write("[ atomtypes ]\n")
            for at in atom_types_list:
                f.write(at.strff)

            f.write("\n")

        f.write("[ bondtypes ]\n")

        for b in unique_bonds:
            f.write(b.strff)

        f.write("\n")
        f.write("[ angletypes ]\n")
        for a in unique_angles:
            f.write(a.strff)
        
        f.write("\n")
        f.write("[ dihedraltypes ]\n")
        for d in unique_dihedrals:
            f.write(d.strff)

        f.close()

    def write_mol(self,o_name):
        """
        Function that writes the molecule
        
        Args:
        ----
        o_name(str): The name of the output file

        Return:
        ------
        Molecular .itp file
        """
        atoms_dic = self.atoms
        bonds_list = self.bonds_list
        angles_list= self.angles_list
        dihedrals_list = self.dihedrals_list

            
        f = open(o_name,"w")
        for line in self.info["[ moleculetype ]"]:
            f.write(line + "\n")
        
        f.write("\n")
        f.write("[ atoms ]\n")
        for key in atoms_dic:
            f.write(atoms_dic[key].strmol)

        if "[ pairs ]" in self.info:
            f.write("\n")
            for line in self.info["[ pairs ]"]:
                f.write(line + "\n")

        f.write("\n")
        f.write("[ bonds ]\n")
        for b in bonds_list:
            f.write(b.strmol)

        f.write("\n")
        f.write("[ angles ]\n")
        for a in angles_list:
            f.write(a.strmol)
        
        f.write("\n")
        f.write("[ dihedrals ]\n")
        for d in dihedrals_list:
            f.write(d.strmol)

        f.close() 
