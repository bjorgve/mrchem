set(_h2o_energy_blyp_launcher " ")

if(ENABLE_MPI)
    set(_h2o_energy_blyp_launcher "mpirun -np 2")
endif()

add_integration_test(
  NAME "H2O_energy_BLYP"
  LABELS "mrchem;h2o_energy_blyp;H2O_energy_BLYP;dft;gga;energy"
  COST 100
  LAUNCH_AGENT ${_h2o_energy_blyp_launcher}
  )
