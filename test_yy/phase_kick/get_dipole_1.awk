BEGIN{counter = 0}
{
if ($4 == "TD_dt") dt = $5;
if($1 == "total_dipole:")
  {
   if (counter == 0){ xinit = $2; yinit = $3; zinit = $4}
   if (counter > 0) print (counter) * dt, $2-xinit, $3-yinit, $4-zinit;
   counter += 1
  }
}
