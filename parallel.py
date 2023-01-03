import os
import time
import shutil

def copying_files(i):
	os.mkdir('output'+str(i))
	#shutil.copy("iesh.f90",'output'+str(i))
	#shutil.copy("iesh.inp",'output'+str(i))
	#shutil.copy("iesh.o",'output'+str(i))
	shutil.copy("15_knot_points_w.txt",'output'+str(i))
	shutil.copy("15_knot_points_x.txt",'output'+str(i))
	shutil.copy("knot_points_x.inp",'output'+str(i))
	shutil.copy("knot_points_w.inp",'output'+str(i))
	#shutil.copy("mod_iesh.f90",'output'+str(i))
	#shutil.copy("mod_iesh.mod",'output'+str(i))
	#shutil.copy("mod_iesh.o",'output'+str(i))
	#shutil.copy("iesh_d_classical_broadening.f90",'output'+str(i))
	shutil.copy("iesh_d_classical.f90",'output'+str(i))
	shutil.copy("sub.sh",'output'+str(i))
	shutil.copy("a.out",'output'+str(i))
	os.chdir('output'+str(i))
	file1=open("ifolder.inp","w")
	file1.write(str(i))
	file1.close()
	os.system('qsub sub.sh')
	#os.system('qsub -l nodes=node5 sub.sh')
	os.chdir('..')

#os.system("rm -r output*")

for i in range(1,101):
	copying_files(i)

