 awk -f get_dipole.awk rt_x.out > dipole_x
 awk -f get_dipole.awk rt_y.out > dipole_y
 awk -f get_dipole.awk rt_z.out > dipole_z
 cut -d " " -f 1 dipole_x > 1
 cut -d " " -f 2 dipole_x > 2
 cut -d " " -f 3 dipole_y > 3
 cut -d " " -f 4 dipole_z > 4
 paste 1 2 3 4 > dipole_xyz
