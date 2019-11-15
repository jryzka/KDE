from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import statistics as stat
import time
import math
from ROOT import TH2F, TCanvas, TLatex, TColor, gStyle, TStyle


pd = open("toy_200k_0per_Dalitz.txt", "r")
pxp = open("toy_200k_0per_par.txt", "r")
pxa = open("toy_200k_0per_apar.txt", "r")


def odczyt(input_file):
	out_list = []
	for i, elem in enumerate(input_file):
		temp = elem.split()
		out_list.append([float(tmp) for tmp in temp])
	input_file.close()
	return out_list

m2pk = []
m2kpi = []
ID = []

for elem in pd:
	temp = elem.split()
	m2pk.append(float(temp[0]))
	m2kpi.append(float(temp[1]))
	ID.append(float(temp[2]))

pd.close()

x_par = odczyt(pxp)
x_apar = odczyt(pxa)

num_of_slices = len(x_par)

y_max = round(max(m2kpi) + 0.01,2)
y_min = round(min(m2kpi) - 0.01,2)

y_bin = (y_max - y_min)/num_of_slices
y_range = np.arange(y_min,y_max+y_bin,y_bin)

which_plot = 0

plt.rcParams['font.size'] = 12

which_plot = int(input("What do you want to plot: \n 1: Dalitz \n 2: KDE Triangle (one region) \n 3: KDE Gauss (one region) \n 4: Parameter h (triangle) (one region) \n 5: KDE Triangle \n 6: KDE Gauss \n 7: Parameter h (triangle) \n 0: None (exit) \n Chosen plot: "))

while(which_plot != 0):
	k = 0
	if which_plot < 5:
		fig = plt.figure(figsize = (20,10))
	else:
		fig, axes = plt.subplots(2, num_of_slices//2, figsize = (20,15))
		plt.subplots_adjust(wspace = 0.4, hspace = 0.4)

	if which_plot == 1:
		plt.hist2d(m2pk, m2kpi, bins = (30,30), cmap = 'YlGnBu')
		plt.colorbar(orientation='vertical')
		plt.xlabel(r'$m^2_{pK} \;\ [GeV^2/c^4]$')
		plt.ylabel(r'$m^2_{K\pi} \;\ [GeV^2/c^4]$')
		plt.axis([2.0, 5.5, 0.25, 2.5])
		plt.show()

	if which_plot == 2:
		k = 2
		ptp_1 = open("density_func_0per_triangle_2_of_4reg_par.txt", "r")
		pta_1 = open("density_func_0per_triangle_2_of_4reg_apar.txt", "r")
		y_weight_p_t_1 = odczyt(ptp_1)
		y_weight_a_t_1 = odczyt(pta_1)
		plt.plot(x_par[k], y_weight_p_t_1[k], 'b-', alpha = 0.5, label = 'particles')
		plt.plot(x_apar[k], y_weight_a_t_1[k], 'g-', alpha = 0.8, label = 'antiparticles')
		plt.xlabel(r'$m^2_{pK} \;\ [GeV^2/c^4]$')
		plt.ylabel(r'$\mathrm{density \ function} \;\ [1/GeV^2/c^4]$')
		plt.show()

	if which_plot == 3:
		k = 2
		pgp_1 = open("density_func_0per_gauss_2_of_4reg_par.txt", "r")
		pga_1 = open("density_func_0per_gauss_2_of_4reg_apar.txt", "r")
		y_weight_p_g_1 = odczyt(pgp_1)
		y_weight_a_g_1 = odczyt(pga_1)
		plt.plot(x_par[k], y_weight_p_g_1[k], 'b-', alpha = 0.5, label = 'particles')
		plt.plot(x_apar[k], y_weight_a_g_1[k], 'g-', alpha = 0.8, label = 'antiparticles')
		plt.xlabel(r'$m^2_{pK} \;\ [GeV^2/c^4]$')
		plt.ylabel(r'$\mathrm{density \ function} \;\ [1/GeV^2/c^4]$')
		plt.show()

	if which_plot == 4:
		k = 2
		php_1 = open("h_opt_0per_triangle_2_of_4_par.txt", "r")
		pha_1 = open("h_opt_0per_triangle_2_of_4_apar.txt", "r")
		hopt_p_g_1 = odczyt(php_1)
		hopt_a_g_1 = odczyt(pha_1)
		plt.plot(x_par[k], hopt_p_g_1[k], 'b-', alpha = 0.5, label = 'particles')
		plt.plot(x_apar[k], hopt_a_g_1[k], 'g-', alpha = 0.8, label = 'antiparticles')
		plt.xlabel(r'$m^2_{pK} \;\ [GeV^2/c^4]$')
		plt.ylabel(r'$h_{i}$')
		plt.show()

	elif which_plot == 5:

		ptp = open("density_func_0per_triangle_All_of_4reg_par_BC_NM.txt", "r")
		pta = open("density_func_0per_triangle_All_of_4reg_apar_BC_NM.txt", "r")
		y_weight_p_t = odczyt(ptp)
		y_weight_a_t = odczyt(pta)


		for i in range(0,2):
			for j in range(0,num_of_slices//2):
				axes[i,j].plot(x_par[k], y_weight_p_t[k], 'b-', label = 'particles')
				axes[i,j].plot(x_apar[k], y_weight_a_t[k], 'g-', label = 'antiparticles')
				axes[i,j].hist(x_par[k], 100, density=True, color='b', alpha = 0.5)
				axes[i,j].hist(x_apar[k], 100, density=True, color='g', alpha = 0.4)
				axes[i,j].set_xlabel(r'$m^2_{pK} \;\ [GeV^2/c^4]$')
				axes[i,j].set_ylabel(r'$ \mathrm{density \ function} \;\ [1/GeV^2/c^4]$')
				axes[i,j].set_title(r'$m^2_{K\pi} \in ($' + str(round(y_range[k],2)) + ',' + str(round(y_range[k+1],2)) + r'$) \;\ GeV^2/c^4$')
				axes[0,0].legend(loc = 'best')
				k = k+1
		plt.show()

	elif which_plot == 6:

		pgp = open("density_func_0per_gauss_All_of_4reg_par.txt", "r")
		pga = open("density_func_0per_gauss_All_of_4reg_apar.txt", "r")
		y_weight_p_g = odczyt(pgp)
		y_weight_a_g = odczyt(pga)


		for i in range(0,2):
			for j in range(0,num_of_slices//2):
				axes[i,j].plot(x_par[k], y_weight_p_g[k], 'b-', alpha = 0.5, label = 'particles')
				axes[i,j].plot(x_apar[k], y_weight_a_g[k], 'g-', alpha = 0.8, label = 'antiparticles')
				axes[i,j].set_xlabel(r'$m^2_{pK} \;\ [GeV^2/c^4]$')
				axes[i,j].set_ylabel(r'$\mathrm{density \ function} \;\ [1/GeV^2/c^4]$')
				axes[i,j].set_title(r'$m^2_{K\pi} \in ($' + str(round(y_range[k],2)) + ',' + str(round(y_range[k+1],2)) + r'$) \;\ GeV^2/c^4$')
				axes[0,0].legend(loc = 'best')
				k = k+1
		plt.show()

	elif which_plot == 7:
		php = open("h_opt_0per_triangle_All_of_4_par.txt", "r")
		pha = open("h_opt_0per_triangle_All_of_4_apar.txt", "r")
		hopt_par = odczyt(php)
		hopt_apar = odczyt(pha)

		for i in range(0,2):
			for j in range(0,num_of_slices//2):
				axes[i,j].plot(x_par[k], hopt_par[k], 'b-', alpha = 0.5, label = 'particles')
				axes[i,j].plot(x_apar[k], hopt_apar[k], 'g-', alpha = 0.8, label = 'antiparticles')
				axes[i,j].set_xlabel(r'$m^2_{pK} \;\ [GeV^2/c^4]$')
				axes[i,j].set_ylabel(r'$h_{i}$')
				axes[i,j].set_title(r'$m^2_{K\pi} \in ($' + str(round(y_range[k],2)) + ',' + str(round(y_range[k+1],2)) + r'$) \;\ GeV^2/c^4$')
				axes[0,0].legend(loc = 'best')
				k = k+1
		plt.show()

	elif which_plot == 0:
		break

	k = 0
	which_plot = int(input("What do you want to plot: \n 1: Dalitz \n 2: KDE Triangle (one region) \n 3: KDE Gauss (one region) \n 4: Parameter h (triangle) (one region) \n 5: KDE Triangle \n 6: KDE Gauss \n 7: Parameter h (triangle) \n 0: None (exit) \n Chosen plot: "))

