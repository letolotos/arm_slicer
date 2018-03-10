import matplotlib.pyplot as plt
import numpy as np
import pyfits as pf
import subprocess as sp
import os
import shutil
import sys 
from copy import copy
from scipy.optimize import curve_fit, newton
from scipy.stats import skewnorm
from glob import glob


def skew_mode(a, ksi, omega):
    delta = a / np.sqrt(1 + a**2)
    mu_z = np.sqrt(2 / np.pi) * delta
    sigma_z = np.sqrt(1 - mu_z**2)
    gamma1 = (4 - np.pi) * mu_z**3 / (2 * (1 - mu_z**2)**(3/2))
    m_0 = mu_z - gamma1 * sigma_z / 2 - a / (2 * np.abs(a)) * np.exp(-2*np.pi/np.abs(a))
    return ksi + omega * m_0

    

def foo(x, a, b, loc, scale):
    return skewnorm.pdf(b * x, a, loc, scale)
    

def foo1(x, a, b, loc, scale, d):
    return foo(x, a, b, loc, scale) - d


def fit_the_slice(xmin, ymin, xmax, ymax, R):
    I0 = copy(galaxy[:])
    r_step, l_step = .02, .5
    alfa = np.arctan2((ymax - ymin), (xmax - xmin))
    L = int(R / r_step) + 1
    B = np.zeros((5, L))
    
    for j in range(-2,3):
        xmin1, xmax1 = xmin + j * l_step * np.sin(alfa), xmax + j * l_step * np.sin(alfa)
        ymin1, ymax1 = ymin - j * l_step * np.cos(alfa), ymax + j * l_step * np.cos(alfa)
        I1 = np.zeros(L)
        for k in range(L):
            x = xmin1 + k * r_step * np.cos(alfa) + xcen
            y = ymin1 + k * r_step * np.sin(alfa) + ycen
            l, m = int(x), int(y)
            b, d = x - l, y - m
            s1, s2, s3, s4 = (1 - d) * (1 - b), (1 - b) * d, d * b, b * (1- d)
            I1[k] = s1 * I0[m,l] + s2 * I0[m+1,l] + s3 * I0[m+1,l+1] + s4 * I0[m,l+1]
    
        if np.abs(ymin) > np.abs(ymax):
            I1 = np.flip(I1, axis = 0)
        
        B[j+2] = I1
    
    Flux = np.mean(np.array(B), axis = 0)
    Flux -= Flux.min()
    Flux_sigma = np.std(B, axis=0)
    x_R = np.linspace(-R/2, R/2, np.size(Flux))
    
    popt = curve_fit(foo, x_R, Flux, p0 = [1, 1, 0, 1])[0]
    y_fit = foo(x_R, popt[0], popt[1], popt[2], popt[3])
    mode = skew_mode(popt[0], popt[2], popt[3]) / popt[1]
    y_mode = foo(mode, popt[0], popt[1], popt[2], popt[3])
    
    x_left = newton(foo1, -R/4, args=(popt[0], popt[1], popt[2], popt[3], y_mode / 2), tol=1e-8)
    y_left = foo(x_left, popt[0], popt[1], popt[2], popt[3])
    x_right = newton(foo1, R/4, args=(popt[0], popt[1], popt[2], popt[3], y_mode / 2), tol=1e-8)
    y_right = foo(x_right, popt[0], popt[1], popt[2], popt[3])
    
    fwhm = x_right - x_left # in pixels
    x_peak, y_peak = xmin + (mode + R/2) * np.cos(alfa), ymin + (mode + R/2) * np.sin(alfa)
    
    fig2, ax2 = plt.subplots()
    ax2.plot(x_R, Flux, 'go', markersize = 2)
    ax2.plot(x_R, y_fit, 'r-', x_left, y_left, 'bo', x_right, y_right, 'bo', mode, y_mode, 'cs')
    ax2.fill_between(x_R, Flux - Flux_sigma, Flux + Flux_sigma, 
                        alpha = .1, edgecolor = '.5', facecolor = '.75', linewidth = .5)
    ax2.set(ylabel = 'Flux', title = r'$\varphi \approx$' + str(alfa))
    fig2.savefig(path + '/' + str(alfa + np.pi) + '.eps', bbox_inches="tight")
    plt.close(fig2)
    
    return fwhm, x_peak, y_peak
    

def func(t, a, b):
    return a * np.exp(b*t)  


def slice_the_arm(x, y, Rs):
    os.makedirs(path, exist_ok=True)
    N = np.size(x)
    fis = np.array([np.arctan2(y[i],x[i]) for i in range(N)])
    for i in range(1, N):
        if fis[i] < fis[i-1]:
            fis[i:] += 2 * np.pi
    
    rs = np.array([np.sqrt(x[i]**2 + y[i]**2) for i in range(N)])
    
    step = 10
    M = N - step
    s = 3 # length step
    k = lambda t: (b * np.sin(t) + np.cos(t)) / (b * np.cos(t) - np.sin(t))
    fi = fis[0]
    
    r_to_R_exp = fwhms = x_new = y_new = np.array([])
    for i in range(M):
        n = int(i + step/2 - 1)
        a, b = curve_fit(func, fis[i:i+step], rs[i:i+step])[0]
        while fi < fis[n+1] or i == M - 1:
            r = func(fi, a, b)
            x0, y0 = r * np.cos(fi), r * np.sin(fi)
            
            idx = np.abs(fis - fi).argmin()
            if fi > fis[idx]:
                l = idx
            else:
                l = idx - 1
            
            R = (fi - fis[l]) * (Rs[l+1] - Rs[l]) / (fis[l+1] - fis[l]) + Rs[l]
            if k(fi):
                ksi = -1 / k(fi)
                xmin = x0 - R / (2 * np.sqrt(1 + ksi**2))
                xmax = 2 * x0 - xmin
                ymin, ymax = ksi * (xmin - x0) + y0, ksi * (xmax - x0) + y0
            else:
                print('place0')
                xmin = xmax = x0
                ymin, ymax = y0 - R / 2, y0 + R / 2
            
            try:
                fwhm, x_peak, y_peak = fit_the_slice(xmin, ymin, xmax, ymax, R)
                if fwhm < 100:
                    fwhms, r_to_R_exp = np.append(fwhms, fwhm), np.append(r_to_R_exp, r / R_exp)
                    x_new, y_new = np.append(x_new, x_peak), np.append(y_new, y_peak)
            except:
                print('failed to fit the slice for fi =', fi)
                pass
            
            fi = np.log(np.exp(b * fi) + s * b / (a * np.sqrt(1 + b**2))) / b
            if fi > fis[-1]:
                break
    
    fwhm_kpc = fwhms * px_to_arcsec * arcsec_to_kpc
    
    fig, ax = plt.subplots()
    ax.plot(r_to_R_exp, fwhm_kpc, 'go', markersize = 2)
    ax.set(ylabel = 'Width, $kpc$', xlabel = 'r/$R_s$', title = 'Arm width; ' + galaxy_name
             + ', band ' + band + ', ' + reg_file.split('.')[0])
    fig.savefig(path + '/results.eps', bbox_inches="tight")
    plt.close(fig)
    
    global iter_number
    if iter_number < 4:
        iter_number += 1
        update_path()
        slice_the_arm(x_new, y_new, fwhms * 1.5)


def update_path():
    global path
    path = galaxy_name + '/' + band + '/' + reg_file.split('.')[0] + '/' + str(iter_number)


file_name, bands = sys.argv[1:]
gal_info = open('galaxy_info.dat', 'r').readlines()

arcsec_to_kpc = float(gal_info[1].split()[1])
dict_for_R_exp = {gal_info[i].split()[0] : float(gal_info[i].split()[1]) for i in range(2, 5)}
px_to_arcsec = .396 # for the SDSS camera

for band in bands:
    R_exp = dict_for_R_exp[band]
    fits_name = file_name + '_' + band + '.fits'
    galaxy_name = file_name.split('_')[0]
    galaxy = pf.open(fits_name)[0].data
    xcen, ycen = len(galaxy[0]) / 2, len(galaxy) / 2
    sp.call(["ds9", fits_name])
    for reg_file in glob('*.reg'):
        data_file = open(reg_file, 'r')
        data_lines = data_file.readlines()
        data_file.close()
        
        iter_number = 1
        path = galaxy_name + '/' + band + '/' + reg_file.split('.')[0] + '/' + str(iter_number)
                
        bound_line = data_lines.pop(list(map(lambda x: 'line' in x, data_lines)).index(True))
        R_coord = np.array(bound_line.split('#')[0].strip('line() ').split(',')).astype(float)
        R = np.sqrt((R_coord[0] - R_coord[2])**2 + (R_coord[1] - R_coord[3])**2)

        coordinates = np.array([line.split('#')[0].strip('point() ').split(',') for line in data_lines[3:]]).astype(float) 
        x, y = coordinates[:,0] - xcen, coordinates[:,1] - ycen
        
        Rs = np.linspace(R,R,np.size(x))
        slice_the_arm(x, y, Rs)

