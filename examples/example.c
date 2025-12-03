// Example of MOPAC C API usage

#include <stdio.h>
#include <stdlib.h>

// MOPAC API library header
#include "mopac.h"

int main(void)
{
  // define a water molecule (C structs do not set default values)
  struct mopac_system water;
  water.natom = 3;
  water.natom_move = 0;
  water.charge = 0;
  water.spin = 0;
  water.model = 0; // PM7
  water.epsilon = 1.0;
  water.nlattice = 0;
  water.nlattice_move = 0;
  water.pressure = 0.0;
  water.tolerance = 1.0;
  water.max_time = 3600;
  water.atom = (int*)malloc(sizeof(int)*water.natom);
  water.coord = (double*)malloc(sizeof(double)*water.natom*3);
  water.atom[0] = 1;
  water.coord[0] = -0.02110;
  water.coord[1] = -0.00200;
  water.coord[2] = 0.00000;
  water.atom[1] = 8;
  water.coord[3] = 0.83450;
  water.coord[4] = 0.45190;
  water.coord[5] = 0.00000;
  water.atom[2] = 1;
  water.coord[6] = 1.47690;
  water.coord[7] = -0.27300;
  water.coord[8] = 0.00000;

  // self-consistent field (SCF) calculation
  struct mopac_state state;
  struct mopac_properties output;
  state.mpack = 0; // uninitialized state designation
  mopac_scf(&water, &state, &output);
  if(output.nerror > 0)
  {
    for(int i=0 ; i<output.nerror ; i++)
    {
      printf("SCF ERROR #%d: %s\n", i+1, output.error_msg[i]);
    }
    return 1;
  }
  printf("PM7 heat of water = %lf kcal/mol\n", output.heat);
  destroy_mopac_properties(&output); // clean up memory for output reuse

  // SCF calculation in COSMO solvent
  water.epsilon = 78.4;
  mopac_scf(&water, &state, &output);
  if(output.nerror > 0)
  {
    for(int i=0 ; i<output.nerror ; i++)
    {
      printf("COSMO ERROR #%d: %s\n", i+1, output.error_msg[i]);
    }
    return 1;
  }
  printf("PM7 heat of water in COSMO solvent = %lf kcal/mol\n", output.heat);
  destroy_mopac_properties(&output); // clean up memory for output reuse

  // geometry relaxation
  water.epsilon = 1.0; // reset dielectric constant back to vacuum
  water.natom_move = 3;
  mopac_relax(&water, &state, &output);
  printf("geometry relaxation:\n");
  if(output.nerror > 0)
  {
    for(int i=0 ; i<output.nerror ; i++)
    {
      printf("RELAX ERROR #%d: %s\n", i+1, output.error_msg[i]);
    }
    return 1;
  }
  for(int i=0 ; i<water.natom*3 ; i++)
  {
    printf("%d: %lf -> %lf\n", i+1, water.coord[i], output.coord_update[i]);
    water.coord[i] = output.coord_update[i];
  }
  destroy_mopac_properties(&output); // clean up memory for output reuse

  // vibrational calculation of relaxed geometry
  mopac_vibe(&water, &state, &output);
  if(output.nerror > 0)
  {
    for(int i=0 ; i<output.nerror ; i++)
    {
      printf("VIBE ERROR #%d: %s\n", i+1, output.error_msg[i]);
    }
    return 1;
  }
  printf("vibrational modes:\n");
  for(int i=0 ; i<water.natom_move*3 ; i++)
  {
    printf("%d: %lf 1/cm\n", i+1, output.freq[i]);
  }
  destroy_mopac_properties(&output); // clean up memory for output reuse

  // memory hygiene to not leak memory in a larger program
  free(water.atom);
  free(water.coord);
  destroy_mopac_state(&state); // clean up of remaining memory allocated by MOPAC API

  return 0;
}
