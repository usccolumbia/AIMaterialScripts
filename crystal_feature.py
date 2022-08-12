from pymatgen.io.cif import CifParser
from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
import random
import numpy as np
from numpy import arccos

class CrystalFeature:
    def __init__(self, file):
        self._file = file
        self._parser = CifParser(self._file)
        self.s = self._parser.get_structures()[0]

    def get_aav(self): 
        '''S1.1 Average atomic volume, AAV (Å3)'''
        aav = 0

        # the number of atoms in the unit cell
        num_atoms = len(self.s.sites)
        # unit cell volume
        cell_vol = self.s.volume

        aav = num_atoms / cell_vol

        return aav

    def get_sdlc(self):
        '''S1.2 Standard deviation in Li neighbor count, SDLC'''
        sdlc = 0

        total_LiNC = 0 # total neighbor count
        num_Li = 0 # number of Li
        liNC = [] # each individual neighbor count


        # for each of the NLi lithium atoms in the unit cell
        for site in self.s.sites:
            ele = site.specie
            if ele.symbol == 'Li': 

                neighbors = self.s.get_neighbors(site, 4) # count number of atoms within 4A
                total_LiNC += len(neighbors)
                liNC.append(len(neighbors))
                num_Li += 1

        avg_LiNC = total_LiNC / num_Li

        # st dev calculation
        numerator = 0
        for value in liNC:
            numerator += (value - avg_LiNC) ** 2
        sdlc = (numerator / (num_Li - 1)) ** .5

        return sdlc

    def get_sdli(self):
        '''S1.3 Standard deviation in Li bond ionicity, SDLI'''
        sdli = 0 

        total_eldf = 0 # total elneg difference
        num_bonds = 0 # total number of Li bonds
        eldf = [] # elneg difference of each bonds

        for i in range(len(self.s.sites)): # loops through all the sites
            ele = self.s.sites[i].specie 
            if ele.symbol == 'Li': # only doing Li bonds
                
                for j in range(len(self.s.sites)): # loop through all sites again
                    ele2 = self.s.sites[j].specie
                    
                    if i == j or (ele2.symbol == 'Li' and j < i): # skip itself and any Li before (since counted befoe)
                        continue
                    
                    dist = self.s.sites[i].distance(self.s.sites[j])
                    
                    # distance under 4 is a bond
                    if dist < 4:
                        total_eldf += abs(ele.X - ele2.X) 
                        eldf.append(abs(ele.X - ele2.X))
                        num_bonds += 1
                        

        avg_eldf = total_eldf / num_bonds

        #st dev calculation
        numerator = 0
        for value in eldf:
            numerator += (value - avg_eldf) ** 2
        sdli = (numerator / (num_bonds - 1)) ** .5

        return sdli

    def get_lbi(self):
        '''S1.4 Average Li bond ionicity, LBI'''
        lbi = 0 

        total_eldf = 0 # total elneg difference
        num_bonds = 0 # total number of Li bonds
        eldf = [] # elneg difference of each bonds

        for i in range(len(self.s.sites)): # loops through all the sites
            ele = self.s.sites[i].specie
            if ele.symbol == 'Li': # only doing Li bonds
                
                for j in range(len(self.s.sites)): # loop through all sites again
                    ele2 = self.s.sites[j].specie
                    
                    if i == j or (ele2.symbol == 'Li' and j < i): # skip itself and any Li before (since counted befoe)
                        continue
                    
                    dist = self.s.sites[i].distance(self.s.sites[j])
                    
                    # distance under 4 is a bond
                    if dist < 4:
                        total_eldf += abs(ele.X - ele2.X) 
                        eldf.append(abs(ele.X - ele2.X))
                        num_bonds += 1
                        

        lbi = total_eldf / (num_bonds-1)

        return lbi

    def get_lnc(self):
        '''S1.5 Average Li neighbor count, LNC'''
        lnc = 0

        total_LiNC = 0 # total neighbor count
        num_Li = 0 # number of Li
        liNC = [] # each individual neighbor count

        for site in self.s.sites:
            ele = site.specie
            if ele.symbol == 'Li': # only want neighbor of lithium
                neighbors = self.s.get_neighbors(site, 4)
                total_LiNC += len(neighbors)
                liNC.append(len(neighbors))
                num_Li += 1

        lnc = total_LiNC / num_Li

        return lnc

    def get_llb(self):
        '''S1.6 Average Li-Li bonds per Li, LLB'''
        llb = 0

        total_llb = 0
        num_Li = 0

        for i in range(len(self.s.sites)): # loop through the sites
            ele = self.s.sites[i].specie
            
            if ele.symbol == 'Li': # take only Li for the first
                num_Li += 1
                
                for j in range(len(self.s.sites)):
                    ele2 = self.s.sites[j].specie
                    
                    # only Li so ignore all not Li, itself
                    # don't ignore previous Li because you want the number of bonds per Li, it wont care if it repeats
                    if ele2.symbol != 'Li':
                        continue
                    elif i == j: 
                        continue
                    
                    dist = self.s.sites[i].distance(self.s.sites[j])
                    
                    # if the distance is less than four it is a bond
                    if dist < 4:
                        total_llb += 1
                        

        llb = total_llb/num_Li

        return llb

    def get_sbi(self):
        '''S1.7 Average bond ionicity of sublattice, SBI'''
        sbi = 0

        total_eldf = 0 # total elneg diff btw non li bonds
        num_bonds = 0

        for i in range(len(self.s.sites)):
            ele = self.s.sites[i].specie
            if ele.symbol != 'Li': # not Li bond
                
                for j in range(len(self.s.sites)):
                    ele2 = self.s.sites[j].specie
                    
                    if i == j or (ele2.symbol != 'Li' and j < i): # ignore duplicates for this bc its bond electron diff 
                        continue
                    
                    dist = self.s.sites[i].distance(self.s.sites[j])
                    
                    # bond if less than 4, X is electron neg
                    if dist < 4:
                        total_eldf += abs(ele.X - ele2.X) 
                        num_bonds += 1
                        

        sbi = total_eldf / num_bonds

        return sbi

    def get_snc(self):
        '''S1.8 Average sublattice neighbor count, SNC'''
        snc = 0

        total_notLiNC = 0 # total neighbor count
        num_notLi = 0 # number of notLi
        notliNC = [] # each individual neighbor count

        for site in self.s.sites:
            ele = site.specie
            if ele.symbol != 'Li': # only want neighbor of lithium
                # .get neighbor returns array so the len of the array is how many neighbors
                total_notLiNC += len(self.s.get_neighbors(site, 4))
                notliNC.append(len(self.s.get_neighbors(site, 4)))
                num_notLi += 1

        snc = total_notLiNC / num_notLi

        return snc

    def get_lattice_anion(self):
        '''lattice anion'''
        a = self.s.sites[0] # placeholder
        maxEN = 0

        for site in self.s.sites:
            if site.specie.X > maxEN:
                a = site
                maxEN = site.specie.X

        return a   

    def get_afc(self):
        '''S1.9 Anion framework coordination, AFC'''
        afc = 0

        a = self.get_lattice_anion()

        total_fc = 0 # total number that follow 𝑟𝑖𝑗0 ≤ 𝑟𝑖𝑗 ≤ 𝑟𝑖𝑗0 + 1Å per anion
        num_a = 0 # number of anions
        nn_distances = [] # nearest neighbor distances for code below

        for site in self.s.sites:
            if site.specie == a.specie: # the anions are any atom in the thing that is a N in this case
                num_a += 1
                
                nearest_distance = 9223372036854775807
                nearest_anion = a # where is this even used...
                
                # finds nearest distance
                for other in self.s.sites:
                    if other.specie == a.specie: # only want w/ other anions
                        if other.distance(site) < nearest_distance and other != site: #whoops there's a distance formula...
                            nearest_anion = other
                            nearest_distance = other.distance(site)
                            
                nn_distances.append(nearest_distance)
                        
                # does supercell        
                for other in self.s.sites:
                    if other.specie == a.specie: # only want 2/ other anions
                        if other.specie == a.specie:
                            if (other.distance(site) >= nearest_distance) \
                                and (other.distance(site) <= (nearest_distance+1)): # 𝑟𝑖𝑗0 ≤ 𝑟𝑖𝑗 ≤ 𝑟𝑖𝑗0 + 1Å
                                total_fc += 1

        afc = total_fc / num_a
        return afc

    def get_aasd(self):
        '''S1.10 Average shortest anion-anion separation distance, AASD (Å)'''
        aasd = 0

        a = self.get_lattice_anion()
                
        total_dist = 0
        num_a = 0 # number of anions

        nn_distances = [] # nearest neighbor distances for code below

        for site in self.s.sites:
            if site.specie == a.specie: # the anions are any atom in the thing that is a N in this case
                num_a += 1
                
                nearest_distance = 9223372036854775807
                
                # finds nearest distance
                for other in self.s.sites:
                    if other.specie == a.specie: # only want w/ other anions
                        if other.distance(site) < nearest_distance and other != site: #whoops there's a distance formula...
                            nearest_distance = other.distance(site)
                            
                total_dist += nearest_distance

        aasd = total_dist / num_a

        return aasd

    def get_vpa(self):
        '''S1.11 Volume per anion, VPA (Å3)'''
        vpa = 0

        a = self.get_lattice_anion()
                
        num_a = 0 # number of anions
        for site in self.s.sites:
            if site.specie == a.specie: # the anions are any atom in the thing that is a N in this case
                num_a += 1
        cell_vol = self.s.volume

        vpa = cell_vol / num_a

        return vpa

    def get_lasd(self):
        '''S1.12 Average shortest Li-anion separation distance, LASD (Å)'''
        lasd = 0

        a = self.get_lattice_anion()

        total_lasd = 0
        num_Li = 0

        for site in self.s.sites:
            if site.specie.symbol == 'Li':
                num_Li += 1
                
                nearest_distance = 9223372036854775807
                
                # find nearest distance
                for other in self.s.sites:
                    if other.specie == a.specie: # Li to anion
                        if other.distance(site) < nearest_distance and other != site:
                             nearest_distance = other.distance(site)
                                
                total_lasd += nearest_distance

        lasd = total_lasd / num_Li

        return lasd

    def get_llsd(self):
        'S1.13 Average shortest Li-Li separation distance, LLSD (Å)'''
        llsd = 0

        total_llsd = 0
        num_Li = 0

        for site in self.s.sites:
            if site.specie.symbol == 'Li':
                num_Li += 1
                
                nearest_distance = 9223372036854775807
                
                # nearest distance
                for other in self.s.sites:
                    if other.specie.symbol == 'Li' and other != site: # Li to Li
                        if other.distance(site) < nearest_distance:
                             nearest_distance = other.distance(site)
                                
                total_llsd += nearest_distance

        llsd = total_llsd / num_Li

        return llsd

    def get_ens(self):
        '''S1.14 Average electronegativity of sublattice, ENS '''
        ens = 0

        total_en = 0
        num_notLi = 0

        for site in self.s.sites:
            ele = site.specie
            
            if ele.symbol != 'Li': # excludes lithium
                num_notLi += 1
                total_en += ele.X

        ens = total_en / num_notLi

        return ens

    def get_effective_radii(self):
        '''returns effective radius based on bond ionity'''
        ionic_radius_dict = {'H': 0.31, 'He': 0.28, 'Li': 0.9, 'Be': 0.59, 'B': 0.41, 'C': 0.3, 'N': 1.32, 'O': 1.26, 'F': 1.19, 'Ne': 0.58, 'Na': 1.16, 'Mg': 0.86, 'Al': 0.68, 'Si': 0.54, 'P': 0.52, 'S': 1.7, 'Cl': 1.67, 'Ar': 1.06, 'K': 1.52, 'Ca': 1.14, 'Sc': 0.89, 'Ti': 0.85, 'V': 0.78, 'Cr': 0.74, 'Mn': 0.75, 'Fe': 0.77, 'Co': 0.76, 'Ni': 0.72, 'Cu': 0.82, 'Zn': 0.88, 'Ga': 0.76, 'Ge': 0.77, 'As': 0.66, 'Se': 1.84, 'Br': 1.82, 'Kr': 1.16, 'Rb': 1.66, 'Sr': 1.32, 'Y': 1.04, 'Zr': 0.86, 'Nb': 0.82, 'Mo': 0.78, 'Tc': 0.74, 'Ru': 0.76, 'Rh': 0.75, 'Pd': 0.88, 'Ag': 1.09, 'Cd': 1.09, 'In': 0.94, 'Sn': 0.83, 'Sb': 0.82, 'Te': 2.07, 'I': 2.06, 'Xe': 1.40, 'Cs': 1.81, 'Ba': 1.49, 'La': 1.17, 'Ce': 1.08, 'Pr': 1.06, 'Nd': 1.12, 'Pm': 1.11, 'Sm': 1.10, 'Eu': 1.20, 'Gd': 1.08, 'Tb': 0.98, 'Dy': 1.13, 'Ho': 1.04, 'Er': 1.03, 'Tm': 1.10, 'Yb': 1.08, 'Lu': 1, 'Hf': 0.85, 'Ta': 0.82, 'W': 0.77, 'Re': 0.71, 'Os': 0.71, 'Ir': 0.77, 'Pt': 0.81, 'Au': 1.07, 'Hg': 1.24, 'Tl': 1.33, 'Pb': 1.12, 'Bi': 1.04, 'Po': 0.94, 'At': 0.76, 'Rn': 1.5, 'Fr': 1.94, 'Ra': 1.62, 'Ac': 1.26, 'Th': 1.08, 'Pa': 1.04, 'U': 0.99}
        covalent_radius_dict = CovalentRadius().radius

        effective_radii = {}

        # For each atom, first choose an atomic radius based on the bonding environment (I'm guessing electronegativity?).
        for site in self.s.sites:
            bond_ionicity = site.specie.X
            if bond_ionicity > 2: 
                effective_radii[site] = ionic_radius_dict[site.specie.symbol]
            else:
                effective_radii[site] = covalent_radius_dict[site.specie.symbol]
        
        return effective_radii

    def get_pf(self):
        # S1.15 Packing fraction of full crystal, PF 
        pf = 0

        effective_radii = self.get_effective_radii()
                
        prev_pf = 0 
        curr_pf = 0
        convergence = False

        num_points = 0 # total amount of points created
        res = 0 # total number of points that are in a effective radius


        # convergence is reached when the packing fraction changes by less than 1% between successive evaluations
        while not convergence:
            # add a random point inside of the strucure
            fract_point = [[round(random.uniform(0,1),4),round(random.uniform(0,1),4),round(random.uniform(0,1),4)]]
            cart_point = self.s.lattice.get_cartesian_coords(fract_point)
            
            for site in self.s.sites:
                if site.distance_from_point(cart_point) <= effective_radii[site]:
                    res += 1
                    break # include break if only counted once, ignore if counted for every atom it is by
            
            if num_points % 1000 == 1000 - 1: # bc index of first is 0, packing fraction is evaluated after every addition of 1,000 new points
                curr_pf = res / num_points
                if prev_pf != 0: # to avoid division from 0
                    if (abs(prev_pf - curr_pf) / prev_pf) < .01:
                        convergence = True
                prev_pf = curr_pf
                    
                
            num_points += 1

        pf = curr_pf

        return pf

    def get_spf(self):
        '''S1.16 Packing fraction of sublattice, SPF'''
        spf = 0

        effective_radii = self.get_effective_radii()
        
        prev_spf = 0 
        curr_spf = 0
        convergence = False

        num_points = 0 # total amount of points created
        res = 0 # total number of points that are in a effective radius


        # convergence is reached when the packing fraction changes by less than 1% between successive evaluations
        while not convergence:
            # add a random point inside of the strucure
            fract_point = [[round(random.uniform(0,1),4),round(random.uniform(0,1),4),round(random.uniform(0,1),4)]]
            cart_point = self.s.lattice.get_cartesian_coords(fract_point)
            
            for site in self.s.sites:
                if site.specie.symbol == 'Li':
                    continue
                if site.distance_from_point(cart_point) <= effective_radii[site]:
                    res += 1
                    break # include break if only counted once, ignore if counted for every atom it is by
            
            if num_points % 1000 == 1000 - 1: # bc index of first is 0, packing fraction is evaluated after every addition of 1,000 new points
                curr_spf = res / num_points
                if prev_spf != 0: # to avoid division from 0
                    if (abs(prev_spf - curr_spf) / prev_spf) < .01:
                        convergence = True
                prev_spf = curr_spf
                    
                
            num_points += 1

        spf = curr_spf

        return spf

    def get_slpw(self):
        '''S1.17 Average straight-line path width, SLPW (Å)'''
        slpw = 0

        effective_radii = self.get_effective_radii()

        # nearest Li pairs
        li_pairs = []
        num_Li = 0

        for site in self.s.sites:
            if site.specie.symbol == 'Li':
                num_Li += 1
                
                nearest_distance = 9223372036854775807
                nearest = site
                
                for other in self.s.sites:
                    if other.specie.symbol == 'Li' and other != site:
                        if other.distance(site) < nearest_distance:
                            nearest_distance = other.distance(site)
                            nearest = other
                
                li_pairs.append([site, nearest])

        total_slpw = 0
        # create straight line
        # ??????
        for pair in li_pairs:
            min_distance = 9223372036854775807
            closest_atom = pair[0]
            
            for atom in self.s.sites:
                if atom == pair[0] or atom == pair[1]:
                    continue
                
                # linear algebra to get angle -- I now realize numpy had a faster way to do this but whatever
                coord1 = pair[0].coords
                coord2 = pair[1].coords
                coord3 = atom.coords
                
                v1 = [coord1[0]-coord2[0], coord1[1]-coord2[1], coord1[2]-coord2[2]] # vector btw lithiums
                v1mag = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]) ** .5 
                v1norm = [v1[0]/v1mag, v1[1]/v1mag, v1[2]/v1mag] # normalize vector
                
                v2 = [coord3[0]-coord2[0], coord3[1]-coord2[1], coord3[2]-coord2[2]] # vector btw li and atom
                v2mag = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]) ** .5
                v2norm = [v2[0]/v2mag, v2[1]/v2mag, v2[2]/v2mag] #normalize vector
                
                res = v1norm[0] * v2norm[0] + v1norm[1] * v2norm[1] + v1norm[2] * v2norm[2] # dot product
                ang = arccos(res)  # take the angle and change to degrees
                
                if ang * (180/np.pi) > 90: # if angle is greater than 90 it is not btw the pair
                    continue
                    
                # total height soh
                height = np.sin(ang) * pair[1].distance(atom)
                distance = height - effective_radii[atom]
                
                if distance < 0:
                    distance = 0
                
                if distance < min_distance:
                    min_distance = distance
            
            total_slpw += min_distance

        slpw = total_slpw / len(li_pairs)

        return slpw

    def get_slpe(self):
        '''S1.18 Average straight-line path electronegativity, SLPE'''
        slpe = 0

        effective_radii = self.get_effective_radii()

        # nearest Li pairs
        li_pairs = []
        num_Li = 0

        for site in self.s.sites:
            if site.specie.symbol == 'Li':
                num_Li += 1
                
                nearest_distance = 9223372036854775807
                nearest = site
                
                for other in self.s.sites:
                    if other.specie.symbol == 'Li' and other != site:
                        if other.distance(site) < nearest_distance:
                            nearest_distance = other.distance(site)
                            nearest = other
                
                li_pairs.append([site, nearest])

        closest_atoms = []
        # create straight line
        # ??????
        for pair in li_pairs:
            min_distance = 9223372036854775807
            closest_atom = pair[0]
            
            for atom in self.s.sites:
                if atom == pair[0] or atom == pair[1]:
                    continue
                
                # linear algebra to get angle -- I now realize numpy had a faster way to do this but whatever
                coord1 = pair[0].coords
                coord2 = pair[1].coords
                coord3 = atom.coords
                
                v1 = [coord1[0]-coord2[0], coord1[1]-coord2[1], coord1[2]-coord2[2]] # vector btw lithiums
                v1mag = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]) ** .5 
                v1norm = [v1[0]/v1mag, v1[1]/v1mag, v1[2]/v1mag] # normalize vector
                
                v2 = [coord3[0]-coord2[0], coord3[1]-coord2[1], coord3[2]-coord2[2]] # vector btw li and atom
                v2mag = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]) ** .5
                v2norm = [v2[0]/v2mag, v2[1]/v2mag, v2[2]/v2mag] #normalize vector
                
                res = v1norm[0] * v2norm[0] + v1norm[1] * v2norm[1] + v1norm[2] * v2norm[2] # dot product
                ang = arccos(res)  # take the angle and change to degrees
                
                if ang * (180/np.pi) > 90: # if angle is greater than 90 it is not btw the pair
                    continue
                    
                # total height soh
                height = np.sin(ang) * pair[1].distance(atom)
                distance = height - effective_radii[atom]
                
                if distance < 0:
                    distance = 0
                
                if distance < min_distance:
                    min_distance = distance
                    closest_atom = atom
            
            closest_atoms.append(atom)

        total_slpe = 0
        for atom in closest_atoms:
            total_slpe += atom.specie.X
        slpe = total_slpe / len(closest_atoms)

        return slpe

    def get_rbi(self):
        '''S1.19 Ratio of average Li bond ionicity to average sublattice bond ionicity, RBI'''
        rbi = 0

        lbi = self.get_lbi()
        sbi = self.get_sbi()

        rbi = lbi / sbi

        return rbi

    def get_rnc(self):
        '''S1.20 Ratio of average Li neighbor count to average sublattice neighbor count, RNC'''
        rnc = 0

        lnc = self.get_lnc()
        snc = self.get_snc()

        rnc = lnc / snc

        return rnc

    def get_feature_names(self):
        '''returns a list of all features'''
        feature_names = ['aav', 'sdlc', 'sdli', 'lbi', 'lnc', 'llb', 'sbi', 'snc', 'afc', 'aasd', \
            'vpa', 'lasd', 'llsd', 'ens', 'pf', 'spf', 'slpw', 'slpe', 'rbi', 'rnc']
        return feature_names

    def get_all_features(self):
        '''returns a dictionary of all features'''
        features = {}

        features['aav'] = self.get_aav()
        features['sdlc'] = self.get_sdlc()
        features['sdli'] = self.get_sdli()
        features['lbi'] = self.get_lbi()
        features['lnc'] = self.get_lnc()
        features['llb'] = self.get_llb()
        features['sbi'] = self.get_sbi()
        features['snc'] = self.get_snc()
        features['afc'] = self.get_afc()
        features['aasd'] = self.get_aasd()

        features['vpa'] = self.get_vpa()
        features['lasd'] = self.get_lasd()
        features['llsd'] = self.get_llsd()
        features['ens'] = self.get_ens()
        features['pf'] = self.get_pf()
        features['spf'] = self.get_spf()
        features['slpw'] = self.get_slpw()
        features['slpe'] = self.get_slpe()
        features['rbi'] = self.get_rbi()
        features['rnc'] = self.get_rnc()

        return features
