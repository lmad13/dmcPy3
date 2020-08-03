#this script should save your coordinates as a variable, and then find the relationship(plot maybe?) between the position
#of the guest molecule and the position of two oxygen atoms? Or maybe the position of the oxygen atoms and energy?
#Or maybe guest molecules and energy?

#For this code, you will have to know which bond you are scanning or which bond you are interested in


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class find_relationships:
	def __init__(self, fileName):
		self.fileName = fileName
	
		
	def guest_position(self):
		readFile = open(self.fileName,'r')
		lines = readFile.readlines()
		self.guest_position_list = []
		
		#in this specific log file, which is Ar1cage6-311+scan.log, the guest molecule displacement is represented 
		#by B60. In another log file, this target will be different
		
		target = "!       B60"
		
		for line in lines:
			if target in line and "Scan" not in line:
				guest_position = line.split()
				self.guest_position_list.append(guest_position[2])
		readFile.close()
		return self.guest_position_list
	
	def collect_SCF(self):
	
		readFile = open(self.fileName,'r')
		lines = readFile.readlines()
		self.SCF_list = []
		
		for line in lines:
			if "SCF Done" in line:
				energy = line.split()
				self.SCF_list.append(energy[4])
		readFile.close()
		return self.SCF_list
		
		
	#In this case, we might be interested in an atom-to-atom distance that is not connected in the input file. So i think
	# it might be better to use the distance formula
		
	def read_coords(self):
		readFile = open(self.fileName,'r')
		lines = readFile.readlines()
		atom_list = []
		coords_list = []
		self.atom_type = []
		Keep_reading = False
		self.scan_count = 0
		for line in lines:
			if "     1          8           0" in line:
				Keep_reading = True
				
				self.scan_count += 1
			if "-----" in line:
				Keep_reading = False
				
			if Keep_reading:
				coords = line.split()
				atom_list.append(coords[0])
				self.atom_type.append(coords[1])
				for i in coords[3:]:
					atom_coords = float(i)
					coords_list.append(atom_coords)
		
		self.atom_array = np.array(atom_list)
		#print(self.scan_count)
		#print(len(atom_list) / self.scan_count)
		self.atom_array = np.reshape(self.atom_array, (int(len(atom_list) / self.scan_count), self.scan_count))
		
		
		self.coords_array = np.array(coords_list)
		self.atom_type_array = np.array(self.atom_type)
		self.atom_number = len(coords_list) / (self.scan_count * 3)
		#print(self.atom_number)
		#print(self.scan_count)
		#print(len(coords_list) / self.scan_count)
		#print(len(coords_list))
		self.coords_array = self.coords_array.reshape(self.scan_count, int(self.atom_number), 3)	# check on the order of this
		
		self.atom_type_array = self.atom_type_array.reshape(self.scan_count, int(self.atom_number), 1)#print(np.shape(self.coords_array))
		return self.atom_array, self.coords_array, self.scan_count, self.atom_number, self.atom_type_array
	
	def normal_vector_eq(self, atom1, atom2, atom3, atom4, atom5):
		self.coords_array_last_set = self.coords_array[-1]
		normal_vector_list = []
		atom_choices = [atom1, atom2, atom3, atom4, atom5]
		for i in range(5):
		for j in range(i + 1, 5):
			for k in range(j + 1, 5):
				#print(atom_choices[i], atom_choices[j], atom_choices[k])

				#x_value = coords_array_last_set[atom_choices[i] - 1 , 0]
				#y_value = coords_array_last_set[rand_atom_list[0] - 1 , 1]
				#z_value = coords_array_last_set[rand_atom_list[0] - 1 , 2]
				atom1_atom2_vector = self.coords_array_last_set[atom_choices[i] - 1] - self.coords_array_last_set[atom_choices[j] - 1]
				atom1_atom3_vector = self.coords_array_last_set[atom_choices[i] - 1] - self.coords_array_last_set[atom_choices[k] - 1]
				normal_vector = np.cross(atom1_atom2_vector, atom1_atom3_vector)
				normal_vector_list.append(normal_vector)
		avg_normal_vector = np.mean(normal_vector_list, axis = 0)
		return self.coords_array_last_set, avg_normal_vector

	def angle_between_planes(self):
	vector_1 = normal_vector_eq(1, 4, 25, 13, 40)
	vector_2 = normal_vector_eq(46, 49, 34, 37, 55)
	

	dot_product = np.dot(vector_1, vector_2)
	magnitude = np.linalg.norm(vector_1) * np.linalg.norm(vector_2)
	angle_radian = np.arccos(dot_product / magnitude)
	vector_angle_degree = (angle_radian * 180.0) / np.pi
	plane_angle_degree = 180.0 - vector_angle_degree
	print(plane_angle_degree)
	
	def center_of_plane(self, atom1, atom2, atom3, atom4, atom5):
		x_com = (self.coords_array_last_set[atom1 - 1, 0] + self.coords_array_last_set[atom2 - 1, 0] + self.coords_array_last_set[atom3 - 1, 0] + self.coords_array_last_set[atom4 - 1, 0] + self.coords_array_last_set[atom5 - 1, 0]) / 5
		y_com = (self.coords_array_last_set[atom1 - 1, 1] + self.coords_array_last_set[atom2 - 1, 1] + self.coords_array_last_set[atom3 - 1, 1] + self.coords_array_last_set[atom4 - 1, 1] + self.coords_array_last_set[atom5 - 1, 1]) / 5
		z_com = (self.coords_array_last_set[atom1 - 1, 2] + self.coords_array_last_set[atom2 - 1, 2] + self.coords_array_last_set[atom3 - 1, 2] + self.coords_array_last_set[atom4 - 1, 2] + self.coords_array_last_set[atom5 - 1, 2]) / 5
		self.center_of_plane = (x_com, y_com, z_com)
		return self.center_of_plane

	def center_of_cage(self):
		water_x = []
		water_y = []
		water_z = []
		atom_type_array_last_set = self.atom_type_array_last_set[-1]

		for i, atom in enumerate(atom_type_array_last_set[:-1]):
			if atom == "8" and atom_array_last_set[i + 1] == "1":
			#print(i, atom)
				water_x.append(coords_array_last_set[i, 0])
				water_y.append(coords_array_last_set[i, 1])
				water_z.append(coords_array_last_set[i, 2])
		center_of_cage_x = np.mean(water_x)
		center_of_cage_y = np.mean(water_y)
		center_of_cage_z = np.mean(water_z)
		
		center_of_cage = (center_of_cage_x, center_of_cage_y, center_of_cage_z)
		return center_of_cage
	# only enumerate atom_array_last_set to the second to the last atom, so it doesn't break when it gets to the last atom
	# make sure to fix this. It will not work if the pattern is not O-H-H
		
	# use this coords_array object in this calc_Atom_Atom_dist			
	def calc_atom_atom_dist(self):	
		self.atom_atom_distances = []
		self.H_bond = []
		self.Covalent_bond = []
		self.H_gas_covalent = []
		self.H_gas_H = []
		self.CH4_gas_H = []
		self.N2_CO_gas_H = []			
		self.CO2_gas_H = []

		self.shortest_H_bond_identity = []
		self.shortest_covalent_bond_identity = []

		for i in range(int(self.atom_number)):
			for j in range(int(self.atom_number)):
			#print(np.shape(self.coords_array[i]))
			#print(np.shape(self.coords_array[i, atom1]))
				atom_atom_vector = self.coords_array[self.scan_count - 1, i] - self.coords_array[self.scan_count - 1, j]
				distance = np.sqrt(np.sum(atom_atom_vector**2))
				self.atom_atom_distances.append(distance)
		self.atom_atom_distances.sort()
		self.atom_atom_distances = list(dict.fromkeys(self.atom_atom_distances))
		#print(self.atom_atom_distances)

		self.H_bond.append(self.atom_atom_distances[41:61])
		#print(type(self.H_bond))
		#print(len(self.H_bond))
		self.Covalent_bond.append(self.atom_atom_distances[1:41])
		
		self.H_gas_H.append(self.atom_atom_distances[42:62])
		self.H_gas_covalent.append(self.atom_atom_distances[2:42])
		
		self.CH4_gas_H.append(self.atom_atom_distances[45:65])
		self.N2_CO_gas_H.append(self.atom_atom_distances[42:62])
		self.CO2_gas_H.append(self.atom_atom_distances[43:63])
		# change your code and count how many hydrogens there are
		# what O-H distances are they having affect on?
		# plot the delta

		#print(self.H_bond)
		#print(self.Covalent_bond)
		#print(self.H_gas_covalent)
		#print(self.H_gas_H)
		return self.H_bond, self.Covalent_bond, self.H_gas_H, self.H_gas_covalent, self.CH4_gas_H, self.N2_CO_gas_H, self.CO2_gas_H
	
	def calc_nearest_neighbor(self, atom):
		self.nearest_atoms = []
		self.far_atoms = []
		#for i in range(self.scan_count):
		for j in range(int(self.atom_number)):
			atom_atom_vector = self.coords_array[4, atom] - self.coords_array[4, j]
			distance = np.sqrt(np.sum(atom_atom_vector**2))
			if distance <= 1.0: #maybe donut shape? 1.5 lower bound. 5 upper bound
				self.nearest_atoms.append(j)
			else:
				self.far_atoms.append(j)
		print(self.nearest_atoms)
		print(self.far_atoms)
	
	
	def calc_angles(self, atom1, atom2, atom3):
		self.angle_list = []
		for i in range(self.scan_count):
			atom1_atom2_vector = self.coords_array[i, atom1] - self.coords_array[i, atom2]
			atom2_atom3_vector = self.coords_array[i, atom3] - self.coords_array[i, atom2]
			dot_product = np.dot(atom1_atom2_vector, atom2_atom3_vector)
			magnitude_1 = np.linalg.norm(atom1_atom2_vector)
			magnitude_2 = np.linalg.norm(atom2_atom3_vector)
			
			angle_radian = np.arccos(dot_product/(magnitude_1 * magnitude_2))
			angle_degree = (angle_radian * 180.0) / np.pi
			self.angle_list.append(angle_degree)
		print(self.angle_list)
		return self.angle_list
		
	#def calc_dihedrals(self, atom1, atom2, atom3, atom4):
	
		
	def plotting_data(self):
		plt.plot(self.guest_position_list, self.SCF_list, 'o', 'r')
		plt.show()

if __name__ == "__main__":	
	#H2_coords = find_relationships("H2O1cagescan6-311g.log")
	#H2_coords.read_coords()
	#H2_coords.calc_atom_atom_dist()
	#print(H2_coords.H_gas_H[0])
	#print("1")

	CO2_coords = find_relationships("CO26-311opt1cage.log")
	CO2_coords.read_coords()
	CO2_coords.calc_atom_atom_dist()
	print("1")

	CO_coords = find_relationships("CO_Clathrate_Opt.log")
	CO_coords.read_coords()
	CO_coords.calc_atom_atom_dist()

	CH4_coords = find_relationships("CH46-311opt1cage.log")
	CH4_coords.read_coords()
	CH4_coords.calc_atom_atom_dist()
	#print(CH4_coords.Covalent_bond[0, 38:])
	#print(len(CH4_coords.Covalent_bond[0]))
	
	Ar_coords = find_relationships("Ar1cage6-311+.log")
	Ar_coords.read_coords()
	Ar_coords.calc_atom_atom_dist()
	#print(Ar_coords.Covalent_bond[0, 38:])
	
	He_coords = find_relationships("He1cageB3LYP6-311+.log")
	He_coords.read_coords()
	He_coords.calc_atom_atom_dist()

	N2_coords = find_relationships("N2_Calthrate_Opt.log")
	N2_coords.read_coords()
	N2_coords.calc_atom_atom_dist()


	empty_coords = find_relationships("Arempty6-311+.log")
	empty_coords.read_coords()
	empty_coords.calc_atom_atom_dist()
	#print(empty_coords.Covalent_bond[0, 38:])



	plt.scatter(range(40), np.array(CO2_coords.Covalent_bond) - empty_coords.Covalent_bond[0][0], label = "CO2")
	plt.scatter(range(40), np.array(Ar_coords.Covalent_bond) - empty_coords.Covalent_bond[0][0], label = "Ar")
	plt.scatter(range(40), np.array(He_coords.Covalent_bond) - empty_coords.Covalent_bond[0][0], label = "He", s = 50)
	plt.scatter(range(40), np.array(CO_coords.Covalent_bond) - empty_coords.Covalent_bond[0][0], label = "CO")
	plt.scatter(range(40), np.array(N2_coords.Covalent_bond) - empty_coords.Covalent_bond[0][0], label = "N2")
	plt.scatter(range(40), np.array(CH4_coords.Covalent_bond) - empty_coords.Covalent_bond[0][0], label = "CH4")
	plt.scatter(range(40), np.array(empty_coords.Covalent_bond) - empty_coords.Covalent_bond[0][0], label = "empty")
	#plt.scatter(range(40), np.array(H2_coords.H_gas_covalent) - empty_coords.Covalent_bond[0][0], label = "H2")

	plt.legend()
	plt.ylim()
	plt.show()

	plt.scatter(range(20), np.array(CO2_coords.CO2_gas_H) - empty_coords.H_bond[0][0], label = "CO2")
	plt.scatter(range(20), np.array(Ar_coords.H_bond) - empty_coords.H_bond[0][0], label = "Ar")
	plt.scatter(range(20), np.array(He_coords.H_bond) - empty_coords.H_bond[0][0], label = "He", s = 50)
	plt.scatter(range(20), np.array(CO_coords.N2_CO_gas_H) - empty_coords.H_bond[0][0], label = "CO")
	plt.scatter(range(20), np.array(N2_coords.N2_CO_gas_H) - empty_coords.H_bond[0][0], label = "N2")
	plt.scatter(range(20), np.array(CH4_coords.CH4_gas_H) - empty_coords.H_bond[0][0], label = "CH4")
	plt.scatter(range(20), np.array(empty_coords.H_bond) - empty_coords.H_bond[0][0], label = "empty")
	#plt.scatter(range(20), np.array(H2_coords.H_gas_H) - empty_coords.H_bond[0][0], label = "H2")

	plt.legend()
	plt.show()
	#guest_molecules = ["empty", "Ar", "H2", CH4", "CO", "N2"]
	#H_bond_values = [empty_coords.H_bond, Ar_coords.H_bond, CH4_coords.H_bond, CO_coords.H_bond, N2_coords.H_bond]
	#print(type(H_bond_values))
	#print(H_bond_values)
	#Covalent_bond_values = [empty_coords.Covalent_bond, Ar_coords.Covalent_bond, H2_coords.Covalent_bond, CO2_coords.Covalent_bond, CH4_coords.Covalent_bond, CO_coords.Covalent_bond, N2_coords.Covalent_bond]
	#for numbers, values in enumerate()
	#plt.scatter(,H_bond_values)
	#plt.scatter(guest_molecules, Covalent_bond_values)
	#plt.ylim(0.94, 0.97)
	#plt.show()
		
		