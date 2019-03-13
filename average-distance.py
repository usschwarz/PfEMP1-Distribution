##############################################################################################################################################################################
# This scipt distributes disc like molecules on the surface of an idealized half-sphere according to a given distribution and calculates the distribution of nearest neighbor distances (measured as surface-to-surface distance). 
# Note definition of theta and phi: theta is the polar angle which is 0 when pointing along the x-axis and pi/2 along the z-axis. phi is the azimuthal angle (0 to 2 pi). 
##############################################################################################################################################################################

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import random

percentage_in_each_section_AA=[0.047076342223133109, 0.17041750404891309, 0.14288596678141868, 0.3429168902432384, 0.2967032967032967]
percentage_in_each_section_AS=[0.03638455687749774, 0.061833159380116055, 0.18647013167014589, 0.25457432455680268, 0.21226516051360983, 0.097094298202712412, 0.15137836879911531]
r_var=np.sqrt(110/np.pi)

#Calculate the distance between two points on the surface of the sphere
def distance_polar(th1,phi1,th2,phi2,R):
	x1=R*np.cos(th1)*np.cos(phi1)
	y1=R*np.cos(th1)*np.sin(phi1)
	z1=R*np.sin(th1)
	x2=R*np.cos(th2)*np.cos(phi2)
	y2=R*np.cos(th2)*np.sin(phi2)
	z2=R*np.sin(th2)
	return distance(x1,y1,z1,x2,y2,z2,R)
	#return R*np.sqrt(2.0-2.0*np.cos(th1)*np.cos(th2)*np.cos(phi1-phi2)-2.0*np.sin(th1)*np.sin(th2))

#Calculate the distance between two points on the surface of the sphere
def distance(x1,y1,z1,x2,y2,z2,r):
	d=np.sqrt(np.power(x1-x2,2)+np.power(y1-y2,2)+np.power(z1-z2,2))
	return r*np.arcsin(d*np.sqrt(4.0*np.power(r,2)-np.power(d,2))*0.5/np.power(r,2))

#When placing a new point, the distance to all previously placed points is calculated and checked if it lies above a given threshold value
def distance_polar_large_enough(theta,phi,R,already_existing_points,new_points,r_thresh):
	#look at previously generated points
	for i in range(len(already_existing_points[0])):
		if already_existing_points[0][i]>=0:
			d=distance_polar(theta,phi,already_existing_points[0][i],already_existing_points[1][i],R)
		else:	
			d=100
		if d<2*r_thresh:
			return False
	#look at points generated in the same ring section
	for i in range(len(new_points[0])):
		if new_points[0][i]>=0:
			d=distance_polar(theta,phi,new_points[0][i],new_points[1][i],R)
		else:	
			d=100
		if d<2*r_thresh:
			return False
	return True

#Given two angles, n particles are distributed uniformly on the ring section that is limited by the two angles. Additionally, points can only be placed if they are far enough form all other points (no overlap of disc like particles). In the rare case that the density is to large to place the particles without overlap, the overlap is allowed. 
def uniform_in_sector(theta_0_real,theta_1_real,n,R,r_thresh,already_existing_points):

	theta_0=np.pi-theta_0_real #conversion to a different definition of theta
	theta_1=np.pi-theta_1_real

	x_0=(1.0-np.cos(theta_0))/2.0
	x_1=(1.0-np.cos(theta_1))/2.0

	uniform_on_sphere=np.zeros((2,n))
	for i in range(n):
		uniform_on_sphere[0][i]=-np.pi/2.0
	points=0
	count=0
	too_small=False

	#make sure n points are generated for the given ring section
	while points < n:
		phi=2*np.pi*(np.random.rand())
		theta=np.pi-np.arccos(1.0-2.0*(x_0+np.random.rand()*(x_1-x_0)))
		if distance_polar_large_enough(theta,phi,R,already_existing_points,uniform_on_sphere,r_thresh) or too_small:
			uniform_on_sphere[0][points]=theta 
			uniform_on_sphere[1][points]=phi 
			points+=1
			count=0

		else:   #try to find another phi placement that does not produce an overlap
			for i in range(100):
				phi=2*np.pi*(np.random.rand())
				if distance_polar_large_enough(theta,phi,R,already_existing_points,uniform_on_sphere,r_thresh):
					uniform_on_sphere[0][points]=theta 
					uniform_on_sphere[1][points]=phi 
					points+=1
					count=0
					break
		count+=1

		#this is executed if the density is too high to prevent overlap
		if count>100:
			too_small=True

	return uniform_on_sphere

#generate a random placement of "number" particles accoring to the given "distribution"
def generate_sample(number,distribution):
	sector_occupation=[0 for i in distribution]
	while sum(sector_occupation)<number:
		x=random.randint(0, len(distribution)-1)
		if random.random()<distribution[x]:
			sector_occupation[x]+=1
	return sector_occupation 

#given the number of PfEMP1 particles these are placed (without overlap) on the surface of an idealized AS knob (half-sphere) and the nearest neighbor distances are calculated
def generate_distances_AS(N_pf,r_thresh):
	thetas=[np.arccos(a/7.0) for a in range(8)]

	R_AA=79.206/2.0
	R_AS=108.289/2.0
	R=R_AS 

	numbers=generate_sample(N_pf,percentage_in_each_section_AS)

	xs_total=[]
	ys_total=[]
	zs_total=[]
	theta_total=[]
	phi_total=[]

	all_angles=[]

	points=numbers[0]
	angles=uniform_in_sector(thetas[0],thetas[1],points,R,r_thresh,all_angles)
	all_angles=angles
	xs=[R*np.cos(angles[0][i])*np.cos(angles[1][i]) for i in range(points)]
	ys=[R*np.cos(angles[0][i])*np.sin(angles[1][i]) for i in range(points)]
	zs=[R*np.sin(angles[0][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)
	theta_total.extend(angles[0])
	phi_total.extend(angles[1])

	points=numbers[1]
	angles=uniform_in_sector(thetas[1],thetas[2],points,R,r_thresh,all_angles)
	if all_angles==[]:
		all_angles=angles
	else:
		all_angles= np.concatenate((all_angles, angles), axis=1)
	xs=[R*np.cos(angles[0][i])*np.cos(angles[1][i]) for i in range(points)]
	ys=[R*np.cos(angles[0][i])*np.sin(angles[1][i]) for i in range(points)]
	zs=[R*np.sin(angles[0][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)
	theta_total.extend(angles[0])
	phi_total.extend(angles[1])

	points=numbers[2]
	angles=uniform_in_sector(thetas[2],thetas[3],points,R,r_thresh,all_angles)
	if all_angles==[]:
		all_angles=angles
	else:
		all_angles= np.concatenate((all_angles, angles), axis=1)
	xs=[R*np.cos(angles[0][i])*np.cos(angles[1][i]) for i in range(points)]
	ys=[R*np.cos(angles[0][i])*np.sin(angles[1][i]) for i in range(points)]
	zs=[R*np.sin(angles[0][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)
	theta_total.extend(angles[0])
	phi_total.extend(angles[1])

	points=numbers[3]
	angles=uniform_in_sector(thetas[3],thetas[4],points,R,r_thresh,all_angles)
	if all_angles==[]:
		all_angles=angles
	else:
		all_angles= np.concatenate((all_angles, angles), axis=1)
	xs=[R*np.cos(angles[0][i])*np.cos(angles[1][i]) for i in range(points)]
	ys=[R*np.cos(angles[0][i])*np.sin(angles[1][i]) for i in range(points)]
	zs=[R*np.sin(angles[0][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)
	theta_total.extend(angles[0])
	phi_total.extend(angles[1])

	points=numbers[4]
	angles=uniform_in_sector(thetas[4],thetas[5],points,R,r_thresh,all_angles)
	if all_angles==[]:
		all_angles=angles
	else:
		all_angles= np.concatenate((all_angles, angles), axis=1)
	xs=[R*np.cos(angles[0][i])*np.cos(angles[1][i]) for i in range(points)]
	ys=[R*np.cos(angles[0][i])*np.sin(angles[1][i]) for i in range(points)]
	zs=[R*np.sin(angles[0][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)
	theta_total.extend(angles[0])
	phi_total.extend(angles[1])

	points=numbers[5]
	angles=uniform_in_sector(thetas[5],thetas[6],points,R,r_thresh,all_angles)
	if all_angles==[]:
		all_angles=angles
	else:
		all_angles= np.concatenate((all_angles, angles), axis=1)
	xs=[R*np.cos(angles[0][i])*np.cos(angles[1][i]) for i in range(points)]
	ys=[R*np.cos(angles[0][i])*np.sin(angles[1][i]) for i in range(points)]
	zs=[R*np.sin(angles[0][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)
	theta_total.extend(angles[0])
	phi_total.extend(angles[1])

	points=numbers[6]
	angles=uniform_in_sector(thetas[6],thetas[7],points,R,r_thresh,all_angles)
	if all_angles==[]:
		all_angles=angles
	else:
		all_angles= np.concatenate((all_angles, angles), axis=1)
	xs=[R*np.cos(angles[0][i])*np.cos(angles[1][i]) for i in range(points)]
	ys=[R*np.cos(angles[0][i])*np.sin(angles[1][i]) for i in range(points)]
	zs=[R*np.sin(angles[0][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)
	theta_total.extend(angles[0])
	phi_total.extend(angles[1])

	distances=[]
	for point1 in range(N_pf): 
		ds=[]
		for j in range(N_pf):
			if point1!=j:
				point2=j
				ds.append(distance(xs_total[point1],ys_total[point1],zs_total[point1],xs_total[point2],ys_total[point2],zs_total[point2],R))
		#record the nearest neighbor distance. 2.0*r_thresh is subtracted because we are interested in the surface to surface distance
		value=np.amin(ds)-2.0*r_thresh
		distances.append(value)
	return distances

#given the number of PfEMP1 particles these are placed (without overlap) on the surface of an idealized AA knob (half-sphere) and the nearest neighbor distances are calculated
def generate_distances_AA(N_pf,r_thresh):
	thetas=[np.arccos(a) for a in [0,0.2,0.4,0.6,0.8,1]]

	R_AA=79.206/2.0
	R_AS=108.289/2.0
	R=R_AA

	numbers=generate_sample(N_pf,percentage_in_each_section_AA)

	xs_total=[]
	ys_total=[]
	zs_total=[]
	theta_total=[]
	phi_total=[]

	all_angles=np.zeros((2,1))
	all_angles[0][0]=-np.pi/2.0

	points=numbers[0]
	angles=uniform_in_sector(thetas[0],thetas[1],points,R,r_thresh,all_angles)
	all_angles=angles
	xs=[R*np.cos(angles[0][i])*np.cos(angles[1][i]) for i in range(points)]
	ys=[R*np.cos(angles[0][i])*np.sin(angles[1][i]) for i in range(points)]
	zs=[R*np.sin(angles[0][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)
	theta_total.extend(angles[0])
	phi_total.extend(angles[1])

	points=numbers[1]
	angles=uniform_in_sector(thetas[1],thetas[2],points,R,r_thresh,all_angles)
	if all_angles==[]:
		all_angles=angles
	else:
		all_angles= np.concatenate((all_angles, angles), axis=1)
	xs=[R*np.cos(angles[0][i])*np.cos(angles[1][i]) for i in range(points)]
	ys=[R*np.cos(angles[0][i])*np.sin(angles[1][i]) for i in range(points)]
	zs=[R*np.sin(angles[0][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)
	theta_total.extend(angles[0])
	phi_total.extend(angles[1])

	points=numbers[2]
	angles=uniform_in_sector(thetas[2],thetas[3],points,R,r_thresh,all_angles)
	if all_angles==[]:
		all_angles=angles
	else:
		all_angles= np.concatenate((all_angles, angles), axis=1)
	xs=[R*np.cos(angles[0][i])*np.cos(angles[1][i]) for i in range(points)]
	ys=[R*np.cos(angles[0][i])*np.sin(angles[1][i]) for i in range(points)]
	zs=[R*np.sin(angles[0][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)
	theta_total.extend(angles[0])
	phi_total.extend(angles[1])

	points=numbers[3]
	angles=uniform_in_sector(thetas[3],thetas[4],points,R,r_thresh,all_angles)
	if all_angles==[]:
		all_angles=angles
	else:
		all_angles= np.concatenate((all_angles, angles), axis=1)
	xs=[R*np.cos(angles[0][i])*np.cos(angles[1][i]) for i in range(points)]
	ys=[R*np.cos(angles[0][i])*np.sin(angles[1][i]) for i in range(points)]
	zs=[R*np.sin(angles[0][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)
	theta_total.extend(angles[0])
	phi_total.extend(angles[1])

	points=numbers[4]
	angles=uniform_in_sector(thetas[4],thetas[5],points,R,r_thresh,all_angles)
	if all_angles==[]:
		all_angles=angles
	else:
		all_angles= np.concatenate((all_angles, angles), axis=1)
	xs=[R*np.cos(angles[0][i])*np.cos(angles[1][i]) for i in range(points)]
	ys=[R*np.cos(angles[0][i])*np.sin(angles[1][i]) for i in range(points)]
	zs=[R*np.sin(angles[0][i]) for i in range(points)]
	xs_total.extend(xs)
	ys_total.extend(ys)
	zs_total.extend(zs)
	theta_total.extend(angles[0])
	phi_total.extend(angles[1])

	distances=[]
	too_small=0
	for point1 in range(N_pf): 
		ds=[]
		for j in range(N_pf):
			if point1!=j:
				point2=j
				d=distance(xs_total[point1],ys_total[point1],zs_total[point1],xs_total[point2],ys_total[point2],zs_total[point2],R)
				#d2=distance_polar(theta_total[point1],phi_total[point1],theta_total[point2],phi_total[point2],R)
				ds.append(d)
		#record the nearest neighbor distance. 2.0*r_thresh is subtracted because we are interested in the surface to surface distance
		value=np.amin(ds)-2.0*r_thresh
		distances.append(value)
	return distances

'''
#method to get all distances instead of nearest neighbors
	distances=[]
	for point1 in range(N_pf): 
		for j in range(N_pf-1-point1):
			point2=point1+1+j
			d=distance(xs_total[point1],ys_total[point1],zs_total[point1],xs_total[point2],ys_total[point2],zs_total[point2],R)
			if d<2*r_thresh:
				print d
			distances.append(d)
'''

def routine_AA(N,molecules,r_var,file_write):
	all_distances=[]
	for i in range(N):
		all_distances=all_distances + generate_distances_AA(molecules,r_var)
	for i in all_distances:
		file_write.write( str(i) + " " )
	file_write.write("\n")
	return all_distances

# the distributions of the nearest neighbor distances are written to the given file
f = open("distribution_data_AA.txt", "w")

all_distances_2=routine_AA(2*3*5*6*7*8*9,2,r_var,f)
all_distances_3=routine_AA(4*5*6*7*8*9,3,r_var,f)
all_distances_4=routine_AA(3*5*6*7*8*9,4,r_var,f)
all_distances_5=routine_AA(4*3*6*7*8*9,5,r_var,f)
all_distances_6=routine_AA(4*5*3*7*8*9,6,r_var,f)
all_distances_7=routine_AA(4*5*6*3*8*9,7,r_var,f)
all_distances_8=routine_AA(4*5*6*7*3*9,8,r_var,f)
all_distances_9=routine_AA(4*5*6*7*8*3,9,r_var,f)
all_distances_10=routine_AA(3*2*6*7*8*9,10,r_var,f)

f.close()

N=float(3*4*5*6*7*8*9)
distance_av=[np.mean(all_distances_2),np.mean(all_distances_3),np.mean(all_distances_4),np.mean(all_distances_5),np.mean(all_distances_6),np.mean(all_distances_7),np.mean(all_distances_8),np.mean(all_distances_9),np.mean(all_distances_10)]
distance_std=[np.std(all_distances_2),np.std(all_distances_3),np.std(all_distances_4),np.std(all_distances_5),np.std(all_distances_6),np.std(all_distances_7),np.std(all_distances_8),np.std(all_distances_9),np.std(all_distances_10)]
distance_std_mean=[i/np.sqrt(N) for i in distance_std]
''' 
distance_av=[26.206709974891435, 18.230267186745269, 14.0906708662481, 11.532365731186836, 9.718260452951915, 8.3779768532398613, 7.3160934904063595, 6.5050765140141209, 5.8223726019055446]
distance_std=[14.064025181998504, 11.886542948684225, 10.09912037020964, 8.8186686344861851, 7.8313933213736107, 7.0572033204525386, 6.4149936928579443, 5.8968353393152864, 5.4695682096597471]
'''

#Plot distirbution of nearest neighbor distances
plt.hist([all_distances_3,all_distances_5,all_distances_7,all_distances_9], 50,edgecolor='None', rwidth=0.9,label=['3 molecules', '5 molecules', '7 molecules', '9 molecules'])
plt.legend()
plt.xlabel('nearest neighbour distance (nm)')
plt.ylabel('events')
plt.savefig('average_distance_1.png')
plt.close()

#Plot average distance as function of molecules per knob
plt.errorbar([2,3,4,5,6,7,8,9,10],distance_av,yerr=distance_std_mean,fmt='o')
plt.xlim(1,11)
plt.xlabel('molecules per knob')
plt.ylabel('average nearest neighbour distance (nm)')
plt.savefig('average_distance_2.png')
plt.close()

