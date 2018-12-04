export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/nas/longleaf/home/yiy/softwares-dogwood/qbox/xerces-c-src_2_8_0/lib
QB=/nas/longleaf/home/yiy/softwares-dogwood/qball/qbach_absorbing_potential/qbach_test/src/qball
mpiexec -n 44 $QB rt_x.in > rt_x.out  
mpiexec -n 44 $QB rt_y.in > rt_y.out  
mpiexec -n 44 $QB rt_z.in > rt_z.out  
mpiexec -n 44 $QB rt_0.in > rt_0.out  
