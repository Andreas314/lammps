/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(spin/magnon,ComputeMagnon);
// clang-format on
#else

#ifndef LMP_COMPUTE_MAGNON_H
#define LMP_COMPUTE_MAGNON_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeMagnon : public Compute {
 public:
  ComputeMagnon(class LAMMPS *, int, char **);
  ~ComputeMagnon() override;
  void init() override {}
  void compute_local() override;
  void compute_array() override;

 private:
  void create_path();
  int knum, pointsnum, timestep, tsize;
  double omegamax, omegamin, omegastep;
  double ***spint;
  double **kpoints, **kvecs;
  char **knames;
  char *id;
  double *kdistances;
};

}    // namespace LAMMPS_NS

#endif
#endif
