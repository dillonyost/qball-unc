&control
    calculation = 'scf',
    prefix = 'Si54_P_FCC',
    wf_collect=.true.,
    tstress = .true.
    tprnfor = .true.
    pseudo_dir = '.',
    outdir='WF/',
    etot_conv_thr = 1.0D-5
    forc_conv_thr = 1.0D-4
    max_seconds = 41400
    restart_mode = 'from_scratch'
    verbosity = 'high'
 /
 &system
    ibrav =  2,
    celldm(1) = 31.05
    nat =  54,
    ntyp = 2,
    ecutwfc = 65,
    occupations='smearing',
    smearing='fermi-dirac',
    degauss=0.0001
 /
 &electrons
    mixing_beta = 0.7
    conv_thr = 1.0D-10
 /
 &ions
 /
 &cell
 /

ATOMIC_SPECIES
 Si 28.086  Si.pbe-n-van.UPF
 P  30.97   P.pbe-n-van.UPF

ATOMIC_POSITIONS (alat)
Si 0.00000000 0.00000000 0.00000000
P 0.08333333 0.08333333 0.08333333
Si 0.16666667 0.00000000 0.16666667
Si 0.25000000 0.08333333 0.25000000
Si 0.33333333 0.00000000 0.33333333
Si 0.41666667 0.08333333 0.41666667
Si 0.00000000 0.16666667 0.16666667
Si 0.08333333 0.25000000 0.25000000
Si 0.16666667 0.16666667 0.33333333
Si 0.25000000 0.25000000 0.41666667
Si 0.33333333 0.16666667 0.50000000
Si 0.41666667 0.25000000 0.58333333
Si 0.00000000 0.33333333 0.33333333
Si 0.08333333 0.41666667 0.41666667
Si 0.16666667 0.33333333 0.50000000
Si 0.25000000 0.41666667 0.58333333
Si 0.33333333 0.33333333 0.66666667
Si 0.41666667 0.41666667 0.75000000
Si 0.16666667 0.16666667 0.00000000
Si 0.25000000 0.25000000 0.08333333
Si 0.33333333 0.16666667 0.16666667
Si 0.41666667 0.25000000 0.25000000
Si 0.50000000 0.16666667 0.33333333
Si 0.58333333 0.25000000 0.41666667
Si 0.16666667 0.33333333 0.16666667
Si 0.25000000 0.41666667 0.25000000
Si 0.33333333 0.33333333 0.33333333
Si 0.41666667 0.41666667 0.41666667
Si 0.50000000 0.33333333 0.50000000
Si 0.58333333 0.41666667 0.58333333
Si 0.16666667 0.50000000 0.33333333
Si 0.25000000 0.58333333 0.41666667
Si 0.33333333 0.50000000 0.50000000
Si 0.41666667 0.58333333 0.58333333
Si 0.50000000 0.50000000 0.66666667
Si 0.58333333 0.58333333 0.75000000
Si 0.33333333 0.33333333 0.00000000
Si 0.41666667 0.41666667 0.08333333
Si 0.50000000 0.33333333 0.16666667
Si 0.58333333 0.41666667 0.25000000
Si 0.66666667 0.33333333 0.33333333
Si 0.75000000 0.41666667 0.41666667
Si 0.33333333 0.50000000 0.16666667
Si 0.41666667 0.58333333 0.25000000
Si 0.50000000 0.50000000 0.33333333
Si 0.58333333 0.58333333 0.41666667
Si 0.66666667 0.50000000 0.50000000
Si 0.75000000 0.58333333 0.58333333
Si 0.33333333 0.66666667 0.33333333
Si 0.41666667 0.75000000 0.41666667
Si 0.50000000 0.66666667 0.50000000
Si 0.58333333 0.75000000 0.58333333
Si 0.66666667 0.66666667 0.66666667
Si 0.75000000 0.75000000 0.75000000

K_POINTS (automatic)
 8 8 8 0 0 0
