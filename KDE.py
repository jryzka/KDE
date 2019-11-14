from __future__ import division
from ROOT import TFile
import numpy as np
import matplotlib.pyplot as plt
import statistics as stat
import scipy
from scipy import integrate
import time
import math
import tkinter
import random
import os

class KDE():
	def __init__(self, wch_file, number_of_slices = 4):
		self.num_of_slices = number_of_slices 	#number of regions Dalitz plot will be divided at
		self.tog = 0				#triangle or gauss, gauss won't be used anymore
		self.reg = number_of_slices
		self.m2pk = []				#data containers
		self.m2kpi = []
		self.ID = []
		self.x_par = [[] for i in range(self.num_of_slices)]
		self.x_apar = [[] for i in range(self.num_of_slices)]
		self.files_root = ["", "_KstarAmp5percentage", "_KstarAmp20percentage"]
		self.files_CP = ["0per", "5per", "20per"]
		self.ToG = ["triangle", "gauss"]
		self.wch_file = wch_file		#flag to choose data file 

	def read_data(self):				#reading data from root file and filling adequate lists
		f = TFile.Open("/home/jryzka/Kernel/toyXic2pKpi_200k" + self.files_root[self.wch_file] +".root")
		for event in f.R0Tree:
			self.m2pk.append(event.m2pk)
			self.m2kpi.append(event.m2kpi)
			self.ID.append(event.id)

	def slices(self):
		#maximum and minimum value of m2kpi - it's plotted at y-axis on Dalitz plot, hence the names related with that values start with y				
		y_max = round(max(self.m2kpi) + 0.01,2)
		y_min = round(min(self.m2kpi) - 0.01,2)
		y_bin = (y_max - y_min)/self.num_of_slices
		y_range = np.arange(y_min,y_max+y_bin,y_bin) 	#ranges (widths) of slices (they are parallel to x-axis)
		#create data containers (arrays) for m2pk values - plotted at x-axis -> names start with x. All other calculations will be performed with that data
		#filling and sorting x_par and x_apar
		for i, y in enumerate(self.m2kpi):
			for j, dy in enumerate(y_range):
				if y > dy and y < dy + y_bin:
					if self.ID[i] == 1:
						self.x_par[j].append(self.m2pk[i])
					else:
						self.x_apar[j].append(self.m2pk[i])
					break

		for j in range(0,self.num_of_slices):
			self.x_par[j].sort()
			self.x_apar[j].sort()
		s = 0
		self.tab_reg = []
		while s<self.num_of_slices:
			self.tab_reg.append(str(s))
			s = s+1
		self.tab_reg.append('All')
	
	def save_to_file_Dalitz(self, file_name):		#saving data for Dalitz_plot.
		if os.path.exists(file_name) == False or os.path.getsize(file_name) == 0:
			plik = open(file_name, "w+")
			for i in range(len(self.m2pk)):
				plik.write(str(self.m2pk[i]) + ' ' + str(self.m2kpi[i]) + ' ' + str(self.ID[i]) + '\n')
			plik.close()
			print("Data for Dalitz has been saved to {}".format(file_name))
	
	def save_to_file(self, input_list, output_file_name):
		if os.path.exists(output_file_name) == False or os.stat(output_file_name).st_size == 0:
			output_file = open(output_file_name, "w+")
			for lista in input_list:
				for elem in lista:
					output_file.write(str(elem) + " ")
				output_file.write("\n")
			output_file.close()
	#There is seperate function for Dalitz, because other data is stored in the list of lists

	def choose_kernel(self, kern):		#changing kernel function (triangle is default)
		self.tog = kern

	def choose_number_of_slices(self, num):  #4 number of slices is default
		self.num_of_slices = num

	def choose_region(self, area):		#one can calucalate one chosen region instead of all of them
		if area < self.num_of_slices:	#condition is needed to avoid situation when chosen number of region is greater
			self.reg = area		#than total number of regions
		else:
			self.reg = num_of_slices
		
	def triangle(self, input_list, xi, h):	#triangle and gauus kernel function
		weight = 0.0
		for x in input_list:
			weight = weight + (1 - abs(x-xi)/h)/h
		return weight

	def gauss(self, input_list, xi, h):
		weight = 0.0
		for x in input_list:
			weight = weight + math.exp(-(x-xi)**2/(2*h**2))/(math.sqrt(2*math.pi)*h)
		return weight

	def kernel(self, input_list): #calculating density function (weight) for every xi point of particle and antiparicle data
		wsk = []
		k = 1.25		#correction parameter
		output_weight = [[] for i in range(self.num_of_slices)]	  #list of lists for output (density function)
		self.h_opt_tab = [[] for i in range(self.num_of_slices)]  #list of lists for output (optimised h parameter)
		for i, x_slice in enumerate(input_list):
			if self.reg >=self.num_of_slices or self.reg == i:  #calulating for all regions or for chosen one
				h = k*pow(len(x_slice), -0.2)*np.std(x_slice) #starting smoothing parameter (bandwidth)
				for xi in x_slice:
					if self.tog == 0:
						wsk = [w for w, x in enumerate(x_slice) if abs(x - xi) < h]
						if len(wsk) == 1:
							x_temp = x_slice[wsk[0]]
							weight = (1 - abs(x_temp-xi)/h)/h
						else:
							x_temp = x_slice[wsk[0]:wsk[-1]]
							weight = self.triangle(x_temp, xi, h)
					elif self.tog == 1:
						x_temp = x_slice
						weight = self.gauss(x_temp, xi, h)
					weight = weight/len(x_slice)
					h_opt = h/np.sqrt(weight)  #optimising bandwidth parameter and calculating density again - can be done 										several times
					self.h_opt_tab[i].append(h_opt)
					if self.tog == 0:
						wsk = [w for w, x in enumerate(x_slice) if abs(x - xi) < h_opt]
						if len(wsk) == 1:
							x_temp = x_slice[wsk[0]]
							weight = (1 - abs(x_temp-xi)/h_opt)/h_opt
						else:	
							x_temp = x_slice[wsk[0]:wsk[-1]]
							weight = self.triangle(x_temp, xi, h_opt)
					elif self.tog == 1:
						x_temp = x_slice
						weight = self.gauss(x_temp, xi, h_opt)
					weight = weight/len(x_slice)
					output_weight[i].append(weight)
		return output_weight

def main():
	kde = KDE(0,4)
	kde.read_data()
	kde.slices()
	kde.choose_kernel(1)
	kde.choose_region(2)

	y_weight_p = kde.kernel(kde.x_par)
	hopt_p = kde.h_opt_tab

	y_weight_a = kde.kernel(kde.x_apar)
	hopt_a = kde.h_opt_tab

	dalitz_file = "toy_200k_" + kde.files_CP[kde.wch_file] + "_Dalitz.txt"
	x_file_par = "toy_200k_" + kde.files_CP[kde.wch_file] + "_par.txt"
	x_file_apar = "toy_200k_" + kde.files_CP[kde.wch_file] + "_apar.txt"
	y_file_weight_p = "density_func_" + kde.files_CP[kde.wch_file] + "_" + kde.ToG[kde.tog] + "_" + kde.tab_reg[kde.reg] + "_of_" + str(kde.num_of_slices) + "reg_par.txt"
	y_file_weight_a = "density_func_" + kde.files_CP[kde.wch_file] + "_" + kde.ToG[kde.tog] + "_" + kde.tab_reg[kde.reg] + "_of_" + str(kde.num_of_slices) + "reg_apar.txt"

	file_hopt_p = "h_opt_" + kde.files_CP[kde.wch_file] + "_" + kde.ToG[kde.tog] + "_" + kde.tab_reg[kde.reg] + "_of_" + str(kde.num_of_slices) + "_par.txt"
	file_hopt_a = "h_opt_" + kde.files_CP[kde.wch_file] + "_" + kde.ToG[kde.tog] + "_" + kde.tab_reg[kde.reg] + "_of_" + str(kde.num_of_slices) + "_apar.txt"

	kde.save_to_file_Dalitz(dalitz_file)
	kde.save_to_file(kde.x_par, x_file_par)
	kde.save_to_file(kde.x_apar, x_file_apar)
	kde.save_to_file(y_weight_p, y_file_weight_p)
	kde.save_to_file(y_weight_a, y_file_weight_a)
	kde.save_to_file(hopt_p, file_hopt_p)
	kde.save_to_file(hopt_a, file_hopt_a)

	print("END OF PROGRAM")
	
if __name__ == "__main__":
	main()

