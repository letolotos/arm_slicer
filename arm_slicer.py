# Author: Konstantin Grishunin; k.s.grishunin@gmail.com

import glob
import sys
import subprocess
import os
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import curve_fit, newton
from scipy.stats import skewnorm
import copy
import shutil
	

def fit_the_points(reg_file):
	"""Function to fit the points along a sp. arm with a set of log spirals"""
	datafile = open(reg_file, "r")
	datafile_lines = datafile.readlines()
	datafile.close()
	
	res_file_name = "results_" + galaxy_name + "_" + reg_name + ".dat"
	res_file = open(res_file_name, "a")
	
	x = np.array([]) 
	y = np.array([]) 
	
	for line in datafile_lines[3:]:
		elem = line.split("#")[0]
		if elem.startswith("line"):
			line_coords = elem.strip("line() ").split(",")
			x1l = float(line_coords[0])
			y1l = float(line_coords[1])
			x2l = float(line_coords[2])
			y2l = float(line_coords[3])
			R = np.sqrt((x2l - x1l) ** 2 + (y2l - y1l) ** 2)
			print("R =", R)
			continue
		coords = elem.strip("point() ").split(",") 
		x = np.append(x, float(coords[0]) - xcen)
		y = np.append(y, float(coords[1]) - ycen)
	
	if x[0] >= 0 and y[0] >= 0:
		add_angle = 0
		previous_quadrand = "I"
	elif x[0] < 0:
		add_angle = np.pi
		previous_quadrand = "II or III"
	elif x[0] >= 0 and y[0] < 0:
		add_angle = 2 * np.pi
		previous_quadrand = "IV"
	
	r1 = np.array([x[0], y[0]] / np.sqrt(x[0] * x[0] + y[0] * y[0]))
	r2 = np.array([x[1], y[1]] / np.sqrt(x[1] * x[1] + y[1] * y[1]))
	
	theta = np.arccos(r1[0] * r2[0] + r1[1] * r2[1])
	
	x_clwse = r1[0] * np.cos(theta) + r1[1] * np.sin(theta)
	y_clwse = -r1[0] * np.sin(theta) + r1[1] * np.cos(theta)
	x_cntrclwse = r1[0] * np.cos(theta) - r1[1] * np.sin(theta)
	y_cntrclwse = r1[0] * np.sin(theta) + r1[1] * np.cos(theta)
	
	ep = 1e-5 # threshold 	
	
	if abs(x_clwse - r2[0]) < ep and abs(y_clwse - r2[1]) < ep:
		cont_clwse = -1 # clockwise rotation
	elif abs(x_cntrclwse - r2[0]) < ep and abs(y_cntrclwse - r2[1]) < ep:
		cont_clwse = 1 # counter-clockwise rotation
	else:
		print("Smth is wrong with rot. direction.")
	
	y = y * cont_clwse 		# we always need to move counter clockwise  
	fi_step = 0.025 * np.pi	# setting angle step to move along the spiral arm
	
	r_scaled = np.array([])
	widths = np.array([])
	astry = np.array([])
	astry1 = np.array([])
	
	for i in range(len(x) - 1):
		if x[i] == 0:
			fi1 = np.pi / 2 + add_angle
		else:
			fi1 = np.arctan(y[i] / x[i]) + add_angle 
		r1 = np.sqrt(x[i] * x[i] + y[i] * y[i])
		
		if x[i+1] >= 0 and y[i+1] >= 0:
			current_quadrand = "I"
		elif x[i+1] < 0:
			current_quadrand = "II or III"
		elif x[i+1] >= 0 and y[i+1] < 0:
			current_quadrand = "IV"
		
		if current_quadrand != previous_quadrand and \
            previous_quadrand != "IV":
			add_angle += np.pi
				
		if x[i+1] == 0:
			fi2 = np.pi / 2 + add_angle
		else:
			fi2 = np.arctan(y[i+1] / x[i+1]) + add_angle
		r2 = np.sqrt(x[i+1] * x[i+1] + y[i+1] * y[i+1])
						
		k = np.log(r2/r1) / (fi2 - fi1)
		r0 = r1 * np.exp(-k * fi1)
		
		if i == 0:
			fi0 = fi1
		
		while fi0 < fi2:
			try:
				r_fi0_scaled, width_fi0, astry_fi0, astry1_fi0 = \
                   slice_the_arm(fi0, r0, k, R, cont_clwse)
				width_fi0 = width_fi0 * px_to_arcsec * kpc_to_arcsec
				if np.abs(width_fi0) < R_exp and np.abs(astry_fi0) < 150 \
                      and astry1_fi0 < 150:
					r_scaled = np.append(r_scaled, r_fi0_scaled)
					widths = np.append(widths, width_fi0)
					astry = np.append(astry, np.abs(astry_fi0))
					astry1 = np.append(astry1, astry1_fi0)
			except:
				print("Failed with slicer; fi0 =", fi0)
				pass
			
			fi0 += fi_step
				
		previous_quadrand = current_quadrand
	
	I00_hdu = pf.PrimaryHDU(I00)
	I00_hdu.writeto("i00.fits")
	shutil.move("i00.fits", galaxy_name + "/" + band_name + "/" + reg_name)
	
	width_fitted_again = least_squares_double_fit(r_scaled, widths)
				
	clr = colors.pop(-1)
	ax_band.plot(r_scaled, widths, clr + symb, markersize=2, label=band_name
              + "; " + reg_name)
	ax_band.plot(r_scaled, width_fitted_again, clr + "-")
	
	ax_arm.plot(r_scaled, widths, clr + symb, markersize=2, label=band_name)
	ax_arm.plot(r_scaled, width_fitted_again, clr + "-")
	ax_arm.legend()
	
	astry_fitted_again = least_squares_double_fit(r_scaled, astry)
	
	ax_asym.plot(r_scaled, astry, clr + symb, markersize=2, label=band_name)
	ax_asym.plot(r_scaled, astry_fitted_again, clr + "-")
	ax_asym.legend()
	
	astry1_fitted_again = least_squares_double_fit(r_scaled, astry1)
	
	ax_asym1.plot(r_scaled, astry1, clr + symb, markersize=2, label=band_name)
	ax_asym1.plot(r_scaled, astry1_fitted_again, clr + "-")
	ax_asym1.legend()
	
	res_file.write("Results for " + galaxy_name + ", band " + band_name + 
                    "; ds9-regions: " + reg_name + "\n")
	res_file.write("r/R_s \t width \t abs(asym) \t asym (areas) \n")			
	for i in range(len(r_scaled)):
		res_file.write(str(r_scaled[i]) + " " + str(widths[i]) + 
                        " " + str(astry[i]) + str(astry1[i]) + "\n")
	
	res_file.close()
	shutil.move(res_file_name, galaxy_name + "/" + band_name + "/" + reg_name)
	

def least_squares_double_fit(x, y):
	p0, p1 = np.polyfit(x, y, 1) # fitting with the line
	y_fitted = p0 * x + p1
	weights_for_y = 1 / np.abs(y - y_fitted)
	p00, p11 = np.polyfit(x, y, 1, w = weights_for_y)
	return p00 * x + p11
	


def slice_the_arm(fi0, r0, k, R, cont_clwse):
	"""Function to slice the spiral arm of the galaxy.
	
	The arm is sliced through the point (fi0, r0) - polar coordinates - 
	with the line of length R (in pixels). The width of the arm and 
	asymmetry of the width profile (approximated by a skew gaussian) 
	are found. k and cont_clwse are characteristics of a logarithmic
	spiral used to approximate the arm itself.
	"""
	I0 = copy.copy(galaxy[:])
	r_step = R / 1000.0 # setting a step to move along the slice segment
	x0 = r0 * np.exp(k * fi0) * np.cos(fi0)
	y0 = r0 * np.exp(k * fi0) * np.sin(fi0)
	r_fi0 = np.sqrt(x0 * x0 + y0 * y0)
	r_fi0_scaled = r_fi0 / R_exp
	delta_fi = .5 / r_fi0 # setting a step for slices smoothing
	I = []
	
	for m in range(-2,3):
		theta = fi0 + m * delta_fi
		ksi1 = -(k * np.cos(theta) - np.sin(theta)) / (k * np.sin(theta) + 
                  np.cos(theta)) # corrected for orientation
		alfa = np.arctan(ksi1)
		x_theta = r0 * np.exp(k * theta) * np.cos(theta)
		y_theta = r0 * np.exp(k * theta) * np.sin(theta)
		b = y_theta - ksi1 * x_theta
		xmin = x_theta - R / (2 * np.sqrt(1 + ksi1 * ksi1))		
		ymin = (ksi1 * xmin + b) * cont_clwse + ycen
		xmin += xcen					
		I1 = []
		r = 0
		while r <= R:
			x = xmin + r * np.cos(alfa)
			y = ymin + r * np.sin(alfa)
			a = int(x)
			b = x - a
			c = int(y)
			I00[c,a] = 100 # to keep track
			d = y - c
			s1 = (1 - d) * (1 - b)
			s2 = (1 - b) * d
			s3 = d * b
			s4 = b * (1- d)
			I1.append(s1 * I0[c,a] + s2 * I0[c+1,a] + s3 * I0[c+1,a+1] + 
                         s4 * I0[c,a+1])
			r += r_step
		
		if xmin < 0:
			I.append(list(reversed(I1)))
		else:
			I.append(I1)
	
	II = np.array(I)
	m0 = len(II)
	Flux = (II[0] + II[1] + II[2] + II[3] + II[4]) / m0
	m = len(Flux)
	u = np.linspace(0, len(Flux)-1, len(Flux))
	
	flux_sigma = np.std(II, axis=0)
	
	raw_fig_name = str(fi0) + '_raw.eps'
	
	fig1, ax1 = plt.subplots()
	ax1.plot(u, Flux, "go")
	ax1.set(ylabel = "Flux", title = r'$\varphi \approx$'
         + str(round(fi0 / np.pi, 2)) + '$\pi$')
	fig1.savefig(galaxy_name + "/" + band_name + "/" + reg_name + "/" + 
                  raw_fig_name, bbox_inches="tight")
	plt.close(fig1)
	
	im0 = np.argwhere(Flux[int(m/3):int(2*m/3)] == 
                        np.max(Flux[int(m/3):int(2*m/3)]))[0][0] + int(m/3)
	im1 = np.argwhere(Flux[:im0] == np.min(Flux[:im0]))[0][0]
	im2 = np.argwhere(Flux[im0:] == np.min(Flux[im0:]))[0][0] + im0
	
	Flux_cut = Flux[im1:im2+1]
	Flux_cut = (Flux_cut - np.min(Flux_cut))
	
	flux_sigma_cut = flux_sigma[im1:im2+1] / np.max(Flux_cut)
	Flux_cut = Flux_cut / np.max(Flux_cut)
		
	u = u - u[im0]
	u_cut = u[im1:im2+1]
	
	n = len(u_cut)
	mean = np.mean(u_cut)
	sigma = np.sqrt(np.sum((u_cut - mean)**2) / n)
	u_cut = u_cut / sigma
		
	popt, pcov = curve_fit(func, u_cut, Flux_cut, p0 = [1, 1, mean, sigma])
	
	xm = np.linspace(skewnorm.ppf(.001, popt[0], loc=popt[2], scale=popt[3]), 
                  skewnorm.ppf(.999, popt[0], loc=popt[2], scale=popt[3]), 
                               1000) / popt[3]
	ym1 = func(xm, popt[0], popt[1], popt[2], popt[3])	
	
	area1, area1_err = integrate.quad(lambda x: func(x, popt[0], popt[1], 
                                       popt[2], popt[3]), xm[0], 0)
	area2, area2_err = integrate.quad(lambda x: func(x, popt[0], popt[1],
                                       popt[2], popt[3]), 0, xm[-1])
	
	asym1 = (area1 - area2) / (area1 + area2)
	
	y_mid = func(0, popt[0], popt[1], popt[2], popt[3]) / 2
	j_ymax = np.argwhere(ym1 == np.max(ym1))[0][0]
	x0_left = xm[np.argwhere(np.abs(ym1[:j_ymax] - y_mid) == 
                  np.min(np.abs(ym1[:j_ymax] - y_mid)))[0][0]]
	x0_right = xm[np.argwhere(np.abs(ym1[j_ymax:] - y_mid) == 
                   np.min(np.abs(ym1[j_ymax:] - y_mid)))[0][0] + j_ymax]
	
	x_left = newton(func1, x0_left, args=(popt[0], popt[1], popt[2], 
                                           popt[3], y_mid), tol=1e-8)
	x_right = newton(func1, x0_right, args=(popt[0], popt[1], popt[2], 
                                             popt[3], y_mid), tol=1e-8)
	arm_width = (x_right - x_left) * sigma * r_step
	
	fig2, ax2 = plt.subplots()
	ax2.plot(u_cut, Flux_cut, "go", label = "Flux")
	ax2.plot(xm, ym1, "r-", label = "Fit")
	ax2.errorbar(u_cut, Flux_cut, yerr = flux_sigma_cut, linestyle = "None")
	ax2.plot(x_left, func(x_left, popt[0], popt[1], popt[2], popt[3]), "b^")
	ax2.plot(x_right, func(x_right, popt[0], popt[1], popt[2], popt[3]), "b^")
	ax2.set(ylabel = "Flux", title = r'$\varphi \approx$' + 
         str(round(fi0 / np.pi, 2)) + '$\pi$')
	ax2.legend()
	fig2.savefig(galaxy_name + "/" + band_name + "/" + reg_name + "/" + 
                  str(fi0) + ".eps", bbox_inches="tight")		
	plt.close(fig2)
		
	return r_fi0_scaled, arm_width, popt[0], asym1
	


def func(x, a, a1, x0, sigma):
	return a1 * skewnorm.pdf(x, a, loc = x0, scale = sigma)

def func1(x, a, a1, x0, sigma, y_mid):
	return func(x, a, a1, x0, sigma) - y_mid


file_name = sys.argv[1]
set_of_bands = sys.argv[2]

galaxy_name = file_name.split("_")[0]

param_file = open("galaxy_info.dat", "r")
param_file_lines = param_file.readlines()
param_file.close()

kpc_to_arcsec = float(param_file_lines[1].split()[1])
px_to_arcsec = float(param_file_lines[2].split()[1])
dict_for_R_exp = {param_file_lines[3].split()[0] : 
                  float(param_file_lines[3].split()[1])}
dict_for_R_exp[param_file_lines[4].split()[0]] = \
               float(param_file_lines[4].split()[1])
dict_for_R_exp[param_file_lines[5].split()[0]] = \
               float(param_file_lines[5].split()[1])

symbols = ["o", "s", "^", "v", "p", "*"]
colors = ["olive", "sienna", "y", "c", "k", "g", "r", "b"]

figures_for_bands = [plt.figure() for j in range(len(set_of_bands))]			
for band_name in set_of_bands:
	
	fig_band = figures_for_bands[set_of_bands.index(band_name)]
	ax_band = fig_band.add_subplot(111)
	ax_band.set(xlabel = "$r/R_s$", ylabel = "Width, $kpc$", 
             title = galaxy_name + ", band " + band_name)
	R_exp = dict_for_R_exp[band_name]
	symb = symbols.pop(-1)
	print("Processing " + galaxy_name + ", band " + band_name + "...")
	
	hdu_galaxy = pf.open(file_name + "_" + band_name + ".fits")
	galaxy = hdu_galaxy[0].data
	I00 = copy.copy(galaxy[:])
	
	xcen = len(galaxy[0]) / 2.0		# image center coords
	ycen = len(galaxy) / 2.0		#
	
	subprocess.call(["ds9", file_name + "_" + band_name + ".fits"])  
	regfiles = glob.glob("*.reg")
	for reg_file in regfiles:
		if band_name == set_of_bands[0] and reg_file == regfiles[0]:
			figures_for_arms = [plt.figure() for j in range(len(regfiles))]
			figures_for_asym = [plt.figure() for j in range(len(regfiles))]
			figures_for_asym1 = [plt.figure() for j in range(len(regfiles))]
		
		fig_arm = figures_for_arms[regfiles.index(reg_file)]
		ax_arm = fig_arm.add_subplot(111)
		ax_arm.set(xlabel = "$r/R_s$", ylabel = "Width, $kpc$", 
            title = galaxy_name + ", arm " + str(regfiles.index(reg_file) + 1))
		
		fig_asym = figures_for_asym[regfiles.index(reg_file)]
		ax_asym = fig_asym.add_subplot(111)
		ax_asym.set(xlabel = "$r/R_s$", ylabel = r"Asymmetry index $\alpha$",
              title = "Asymmetry index for " + galaxy_name + ", arm " + 
              str(regfiles.index(reg_file) + 1))
		
		fig_asym1 = figures_for_asym1[regfiles.index(reg_file)]
		ax_asym1 = fig_asym1.add_subplot(111)
		ax_asym1.set(xlabel = "$r/R_s$", ylabel = "Asymmetry index", 
               title = "Asymmetry index for " + galaxy_name + ", arm " + 
               str(regfiles.index(reg_file) + 1))
			
		reg_name = reg_file.split(".")[0]
		os.makedirs(galaxy_name + "/" + band_name + "/" + reg_name, 
                      exist_ok=True)
		fit_the_points(reg_file)
		shutil.move(reg_file, galaxy_name + "/" + band_name + "/" + reg_name)
		print("Done for " + band_name + " band; " + reg_name + ".")
		
		hdu_galaxy.close()
	
	ax_band.legend()
	fig_band.savefig("widths_" + band_name + ".eps", bbox_inches="tight")
	plt.close(fig_band)

for j in range(len(regfiles)):
	figures_for_arms[j].savefig("arm" + str(j + 1) + ".eps", 
                 bbox_inches="tight")
	figures_for_asym[j].savefig("arm" + str(j + 1) + "_asym.eps", 
                 bbox_inches="tight")
	figures_for_asym1[j].savefig("arm" + str(j + 1) + "_asym1.eps", 
                  bbox_inches="tight")

shutil.move("galaxy_info.dat", galaxy_name)	
