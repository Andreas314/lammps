/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_magnon.h"

#include "angle.h"
#include "atom.h"
#include "atom_masks.h"
#include "bond.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "improper.h"
#include "kspace.h"
#include "modify.h"
#include "memory.h"
#include "pair.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeMagnon::ComputeMagnon(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  if (narg % 4 != 3 || narg < 10 ) error->all(FLERR, 1, "Wrong number of arguments");
  omegamin = atof(arg[0]);
  omegastep = atof(arg[1]);
  omegamax = atof(arg[2]);
  lmp->memory->create(kpoints, narg / 4, 3, "magnon:kpoints");
  lmp->memory->create(kdistances, narg / 4, "magnon:kdistance");
  lmp->memory->create(knames, narg / 4, 1, "magnon:knames");
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */
ComputeMagnon::~ComputeMagnon()
{
	lmp->memory->destroy(kpoints);
	lmp->memory->destroy(kdistances);
	lmp->memory->destroy(knames);
}

void ComputeMagnon::compute_local()
{

}
void ComputeMagnon::compute_array()
{

}
