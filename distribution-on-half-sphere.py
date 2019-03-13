##############################################################################################################################################################################
# This script distributes a given number of particles on a half sphere according to a predefined discrete distribution along the arclength on the sphere. We assume the 
# particles to be distributed uniformly within each circular ring section. The script produces a list of x,y,z-coordinates of the relevant points on the surface of the 
# half-sphere. 
# Note definition of theta and phi: theta is the azimuthal angle (0 to 2 pi) and phi is the polar angle which is 0 along the z-axis. 
##############################################################################################################################################################################

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 

#Given two angles, n particles are distributed uniformly on the ring section that is limited by the two angles
def uniform_in_sector(phi_0,phi_1,n):

	x_0=(1-np.cos(phi_0))/2.0
	x_1=(1-np.cos(phi_1))/2.0

	uniform_on_sphere=np.random.rand(2,n)

	uniform_on_sphere[0]=[2*np.pi*(uniform_on_sphere[0][i]) for i in range(n)] #theta
	uniform_on_sphere[1]=[np.arccos(1-2*(x_0+uniform_on_sphere[1][i]*(x_1-x_0))) for i in range(n)] #phi

	return uniform_on_sphere

#The half-sphere corresponding to the AA RBC is devided into 5 sections along the arclength (circular rings)
def generate_coordinates_AA(fraction, N, R, phis):

	xs_total=[]
	ys_total=[]
	zs_total=[]

	points=int(np.round(N*fraction[0]))
	angles=uniform_in_sector(phis[0],phis[1],points)
	xs=[R*np.sin(angles[1][i])*np.cos(angles[0][i]) for i in range(points)]
	ys=[R*np.sin(angles[1][i])*np.sin(angles[0][i]) for i in range(points)]
	zs=[R*np.cos(angles[1][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)

	points=int(np.round(N*fraction[1]))
	angles=uniform_in_sector(phis[1],phis[2],points)
	xs=[R*np.sin(angles[1][i])*np.cos(angles[0][i]) for i in range(points)]
	ys=[R*np.sin(angles[1][i])*np.sin(angles[0][i]) for i in range(points)]
	zs=[R*np.cos(angles[1][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)

	points=int(np.round(N*fraction[2]))
	angles=uniform_in_sector(phis[2],phis[3],points)
	xs=[R*np.sin(angles[1][i])*np.cos(angles[0][i]) for i in range(points)]
	ys=[R*np.sin(angles[1][i])*np.sin(angles[0][i]) for i in range(points)]
	zs=[R*np.cos(angles[1][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)

	points=int(np.round(N*fraction[3]))
	angles=uniform_in_sector(phis[3],phis[4],points)
	xs=[R*np.sin(angles[1][i])*np.cos(angles[0][i]) for i in range(points)]
	ys=[R*np.sin(angles[1][i])*np.sin(angles[0][i]) for i in range(points)]
	zs=[R*np.cos(angles[1][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)

	points=int(np.round(N*fraction[4]))
	angles=uniform_in_sector(phis[4],phis[5],points)
	xs=[R*np.sin(angles[1][i])*np.cos(angles[0][i]) for i in range(points)]
	ys=[R*np.sin(angles[1][i])*np.sin(angles[0][i]) for i in range(points)]
	zs=[R*np.cos(angles[1][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)

	print xs_total, ys_total, zs_total

def generate_coordinates_AS(fraction, N, R, phis):

	xs_total=[]
	ys_total=[]
	zs_total=[]

	points=int(np.round(N*fraction[0]))
	angles=uniform_in_sector(phis[0],phis[1],points)
	xs=[R*np.sin(angles[1][i])*np.cos(angles[0][i]) for i in range(points)]
	ys=[R*np.sin(angles[1][i])*np.sin(angles[0][i]) for i in range(points)]
	zs=[R*np.cos(angles[1][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)

	points=int(np.round(N*fraction[1]))
	angles=uniform_in_sector(phis[1],phis[2],points)
	xs=[R*np.sin(angles[1][i])*np.cos(angles[0][i]) for i in range(points)]
	ys=[R*np.sin(angles[1][i])*np.sin(angles[0][i]) for i in range(points)]
	zs=[R*np.cos(angles[1][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)

	points=int(np.round(N*fraction[2]))
	angles=uniform_in_sector(phis[2],phis[3],points)
	xs=[R*np.sin(angles[1][i])*np.cos(angles[0][i]) for i in range(points)]
	ys=[R*np.sin(angles[1][i])*np.sin(angles[0][i]) for i in range(points)]
	zs=[R*np.cos(angles[1][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)

	points=int(np.round(N*fraction[3]))
	angles=uniform_in_sector(phis[3],phis[4],points)
	xs=[R*np.sin(angles[1][i])*np.cos(angles[0][i]) for i in range(points)]
	ys=[R*np.sin(angles[1][i])*np.sin(angles[0][i]) for i in range(points)]
	zs=[R*np.cos(angles[1][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)

	points=int(np.round(N*fraction[4]))
	angles=uniform_in_sector(phis[4],phis[5],points)
	xs=[R*np.sin(angles[1][i])*np.cos(angles[0][i]) for i in range(points)]
	ys=[R*np.sin(angles[1][i])*np.sin(angles[0][i]) for i in range(points)]
	zs=[R*np.cos(angles[1][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)

	points=int(np.round(N*fraction[5]))
	angles=uniform_in_sector(phis[5],phis[6],points)
	xs=[R*np.sin(angles[1][i])*np.cos(angles[0][i]) for i in range(points)]
	ys=[R*np.sin(angles[1][i])*np.sin(angles[0][i]) for i in range(points)]
	zs=[R*np.cos(angles[1][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)

	points=int(np.round(N*fraction[6]))
	angles=uniform_in_sector(phis[6],phis[7],points)
	xs=[R*np.sin(angles[1][i])*np.cos(angles[0][i]) for i in range(points)]
	ys=[R*np.sin(angles[1][i])*np.sin(angles[0][i]) for i in range(points)]
	zs=[R*np.cos(angles[1][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)

	print xs_total, ys_total, zs_total


percentage_in_each_section_AA=[0.047076342223133109, 0.17041750404891309, 0.14288596678141868, 0.3429168902432384, 0.2967032967032967]
percentage_in_each_section_AS=[0.03638455687749774, 0.061833159380116055, 0.18647013167014589, 0.25457432455680268, 0.21226516051360983, 0.097094298202712412, 0.15137836879911531]

phis_AA=[np.arcsin(a) for a in [0,0.2,0.4,0.6,0.8,1]]
phis_AS=[np.arcsin(a/7.0) for a in range(8)]

R_AA=79.206/2.0
R_AS=108.289/2.0

generate_coordinates_AS(percentage_in_each_section_AS, 150, R_AS/R_AA, phis_AS)
generate_coordinates_AA(percentage_in_each_section_AA, 120, 1, phis_AA)


