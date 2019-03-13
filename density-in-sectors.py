##############################################################################################################################################################################
# This scipt converts a distribution of particles (on a half-sphere) recorded in x-slices to a distribution in discrete bins along the arclength on the half-sphere.
# To do so, the area patches formed by x and phi cuts through the half sphere were found previously with mathematica (data given in "Areas"). 
# The array "fraction_phi" determines the percentage of particles that would lie within each arclength-bin according to the calculated distribution.
##############################################################################################################################################################################
import numpy as np
import matplotlib.pyplot as plt

count_AA=[98,91,73,66,36,0,0]
count_AS=[57,54,55,44,27,14,11]

N_AA=sum(count_AA)
N_AS=sum(count_AS)

R_AA=79.206/2.0
R_AS=108.289/2.0
r_AA=27.2
r_AS=32.2

#each x-slice has the same surface area
slice_area_AA=0.314159*np.power(R_AA,2)*2
slice_area_AS=0.224399*np.power(R_AS,2)*2

density_x_AA=[i/slice_area_AA for i in count_AA]
density_x_AS=[i/slice_area_AS for i in count_AS]

def density_spherical_coordinates_AA(density_x):
	Areas=[0.129722, 0.0575678, 0.047787, 0.0473456, 0.0317364, 0.136454,0.0633804, 0.0622697, 0.0520555, 0.154092, 0.0871022, 0.0729648,0.20805, 0.106109, 0.314159]
	Areas=[i*np.power(R_AA,2) for i in Areas] #areas now correspond to real areas in nm^2

	#from inner to outer
	areas_phi=[Areas[4],Areas[3]+Areas[8],Areas[2]+Areas[7]+Areas[11],Areas[1]+Areas[6]+Areas[10]+Areas[13],Areas[0]+Areas[5]+Areas[9]+Areas[12]+Areas[14]]
	
	density_phi=[0,0,0,0,36.0/slice_area_AA]

	zone_area=0.314159*np.power(R_AA,2)

	#conversion from discrete density in x-slices to discrete density in phi sections
	density_phi[3]=(zone_area*density_x[3]-density_phi[4]*Areas[13-1])/Areas[14-1]
	density_phi[2]=(zone_area*density_x[2]-density_phi[4]*Areas[10-1]-density_phi[3]*Areas[11-1])/Areas[12-1]
	density_phi[1]=(zone_area*density_x[1]-density_phi[4]*Areas[6-1]-density_phi[3]*Areas[7-1]-density_phi[2]*Areas[8-1])/Areas[9-1]
	density_phi[0]=(zone_area*density_x[0]-density_phi[4]*Areas[1-1]-density_phi[3]*Areas[2-1]-density_phi[2]*Areas[3-1]-density_phi[1]*Areas[4-1])/Areas[5-1]

	arc_lengths=[R_AA*np.arcsin(a) for a in [0.2,0.4,0.6,0.8,1]]
	d=[i/N_AA*1000 for i in density_phi]

	#Calulate fraction of particles in each phi-section
	fraction_phi=[]
	for i in range(5):
		fraction_phi.append(density_phi[i]*areas_phi[i]*2.0/N_AA)
	print fraction_phi

	#make data easily plotable 
	ph_plot=[0,0]
	dph_plot=[0]
	for i in range(5):
		ph_plot.append(arc_lengths[i])
		ph_plot.append(arc_lengths[i])
		dph_plot.append(d[i])
		dph_plot.append(d[i])
	dph_plot.append(0)
	return [ph_plot,dph_plot]


def density_spherical_coordinates_AS(density_x):
	Areas=[0.077595, 0.0336289, 0.0269923, 0.023959, 0.0225391, 0.0235741,0.0161111, 0.0794426, 0.0348242, 0.0284968, 0.0264048, 0.0294374,0.0257937, 0.0835817, 0.0376706, 0.0325264, 0.0365067, 0.0341141,0.0912602, 0.0437597, 0.0461, 0.0432795, 0.106072, 0.0626943,0.0556329, 0.146733, 0.0776669, 0.224399]
	Areas=[i*np.power(R_AS,2) for i in Areas] #areas now correspond to real areas in nm^2

	#inner to outer
	areas_phi=[Areas[6],Areas[5]+Areas[12],Areas[4]+Areas[11]+Areas[17],Areas[3]+Areas[10]+Areas[16]+Areas[21],Areas[2]+Areas[9]+Areas[15]+Areas[20]+Areas[24],
Areas[1]+Areas[8]+Areas[14]+Areas[19]+Areas[23]+Areas[26],Areas[0]+Areas[7]+Areas[13]+Areas[18]+Areas[22]+Areas[25]+Areas[27]]

	density_phi=[0,0,0,0,0,0,11.0/slice_area_AS]

	zone_area=0.224399*np.power(R_AS,2)

	#conversion from discrete density in x-slices to discrete density in phi sections
	density_phi[5]=(zone_area*density_x[5]-density_phi[6]*Areas[26-1])/Areas[27-1]
	density_phi[4]=(zone_area*density_x[4]-density_phi[6]*Areas[23-1]-density_phi[5]*Areas[24-1])/Areas[25-1]
	density_phi[3]=(zone_area*density_x[3]-density_phi[6]*Areas[19-1]-density_phi[5]*Areas[20-1]-density_phi[4]*Areas[21-1])/Areas[22-1]
	density_phi[2]=(zone_area*density_x[2]-density_phi[6]*Areas[14-1]-density_phi[5]*Areas[15-1]-density_phi[4]*Areas[16-1]-density_phi[3]*Areas[17-1])/Areas[18-1]
	density_phi[1]=(zone_area*density_x[1]-density_phi[6]*Areas[8-1]-density_phi[5]*Areas[9-1]-density_phi[4]*Areas[10-1]-density_phi[3]*Areas[11-1]-density_phi[2]*Areas[12-1])/Areas[13-1]
	density_phi[0]=(zone_area*density_x[0]-density_phi[6]*Areas[1-1]-density_phi[5]*Areas[2-1]-density_phi[4]*Areas[3-1]-density_phi[3]*Areas[4-1]-density_phi[2]*Areas[5-1]-density_phi[1]*Areas[6-1])/Areas[7-1]

	arc_lengths=[R_AS*np.arcsin((a+1)/7.0) for a in range(7)]

	d=[i/N_AS*1000 for i in density_phi]

	#Calulate fraction of particles in each phi-section
	fraction_phi=[]
	for i in range(7):
		fraction_phi.append(density_phi[i]*areas_phi[i]*2.0/N_AS)
	print fraction_phi

	#make data easily plotable 
	ph_plot=[0,0]
	dph_plot=[0]
	for i in range(7):
		ph_plot.append(arc_lengths[i])
		ph_plot.append(arc_lengths[i])
		dph_plot.append(d[i])
		dph_plot.append(d[i])
	dph_plot.append(0)
	return [ph_plot,dph_plot]

[ph_plot_AA,dph_plot_AA]=density_spherical_coordinates_AA(density_x_AA)
[ph_plot_AS,dph_plot_AS]=density_spherical_coordinates_AS(density_x_AS)

#Plot density in discrete phi bins as function of arclength on the sphere
plt.plot(ph_plot_AA,dph_plot_AA,label="AA",color="blue")
plt.plot(ph_plot_AS,dph_plot_AS,label="AS",color="orange")
plt.axvline(R_AA*np.arcsin(r_AA/R_AA), linestyle='--',color="blue")
plt.axvline(R_AS*np.arcsin(r_AS/R_AS), linestyle='--',color="orange")
plt.xlabel("arc length (nm)", fontsize='14')
plt.ylabel("normalized density of PfEMP1s", fontsize='14')
plt.legend()
plt.tight_layout()
plt.show()


#Plot histogram of counts in x-sections
'''
bar_width=0.35
d_x_AA=[i for i in density_x_AA]
plt.bar([i+1 for i in np.arange(len(density_x_AA))],d_x_AA,bar_width,label="AA")
d_x_AS=[i for i in density_x_AS]
plt.bar([i+1+ bar_width for i in np.arange(len(density_x_AS))],d_x_AS,bar_width,label="AS")
plt.xlabel("x-slice", fontsize='14')
plt.xticks( [i+1+bar_width/2.0 for i in np.arange(len(density_x_AS))], [i+1 for i in np.arange(len(density_x_AS))])
plt.ylabel("# of PfEMP1s counted", fontsize='14')
plt.legend()
plt.tight_layout()
plt.show()
'''






