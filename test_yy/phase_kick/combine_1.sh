 awk -f get_dipole_1.awk rt_x.out > dipole_x
 awk -f get_dipole_1.awk rt_y.out > dipole_y
 awk -f get_dipole_1.awk rt_z.out > dipole_z
 awk -f get_dipole_1.awk rt_0.out > dipole_0
 python subtract.py dipole_x dipole_0 > dipole_x_0
 python subtract.py dipole_y dipole_0 > dipole_y_0
 python subtract.py dipole_z dipole_0 > dipole_z_0
 cut -d " " -f 1 dipole_x_0 > 1
 cut -d " " -f 2 dipole_x_0 > 2
 cut -d " " -f 3 dipole_y_0 > 3
 cut -d " " -f 4 dipole_z_0 > 4
 paste 1 2 3 4 > dipole_xyz_1
