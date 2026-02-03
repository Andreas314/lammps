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
#include "comm.h"
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
  if (narg % 4 != 3 || narg < 15 ) error->all(FLERR, 1, "Wrong number of arguments");
  //TODO: errors on omegas, knum;
  id = arg[0];
  int newnarg = narg - 3;
  char **newarg = arg + 3;
  omegamin = atof(newarg[0]);
  omegastep = atof(newarg[1]);
  omegamax = atof(newarg[2]);
  knum = atoi(newarg[3]);
  pointsnum = newnarg / 4 - 1;
  tsize = 200;
  lmp->memory->create(kpoints, pointsnum, 3, "magnon:kpoints");
  lmp->memory->create(kdistances, pointsnum * knum, "magnon:kdistance");
  lmp->memory->create(kvecs, pointsnum * knum, 3, "magnon:kvecs");
  lmp->memory->create(knames, pointsnum, 1, "magnon:knames");
  lmp->memory->create(spint, tsize, atom->nmax, 3, "magnon:spint");
  timestep = 0;
  int argnum = 0;
  int i = 0;
  while (argnum < (newnarg - 4) ){
	  //TODO: Names of points with one char
  	strcpy(knames[i], newarg[4 + argnum]);
	kpoints[i][0] = atof(newarg[4 + argnum + 1]);
	kpoints[i][1] = atof(newarg[4 + argnum + 2]);
	kpoints[i][2] = atof(newarg[4 + argnum + 3]);
	i++;
	argnum += 4;
  }
  create_path();
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */
ComputeMagnon::~ComputeMagnon()
{
  lmp->memory->destroy(kpoints);
  lmp->memory->destroy(kvecs);
  lmp->memory->destroy(kdistances);
  lmp->memory->destroy(knames);
  lmp->memory->destroy(spint);
  lmp->memory->destroy(S);
  lmp->memory->destroy(C);
}
void ComputeMagnon::create_path()
{
  double *point1, *point2;
  double dist[3];
  for (int point = 1; point < pointsnum; point++)
  {
    point1 = kpoints[point - 1];
    point2 = kpoints[point];
    dist[0] = (point2[0] - point1[0]) / (double)(knum);
    dist[1] = (point2[1] - point1[1]) / (double)(knum);
    dist[2] = (point2[2] - point1[2]) / (double)(knum);
    double ds = std::sqrt(dist[0]*dist[0] +
			  dist[1]*dist[1] +
			  dist[2]*dist[2]);
    for (int distance = 0; distance <= knum; distance++)
    {
      int indx = distance + (point - 1) * knum;
      if (point == 1 && distance == 0) kdistances[0] = 0;
      else kdistances[indx] = kdistances[indx - 1] + ds;
      kvecs[indx][0] = point1[0] + dist[0] * distance; 
      kvecs[indx][1] = point1[1] + dist[1] * distance; 
      kvecs[indx][2] = point1[2] + dist[2] * distance;
    }
  }
}
void ComputeMagnon::compute_local()
{
  if (timestep > tsize)
  {
    tsize += 10;
    memory->grow(spint, tsize, atom->nmax, 3, "magnon:spint");
  }
  for (int i = 0; i < atom->nmax; i++)
  {
    spint[timestep][i][0] = atom->sp[i][0];
    spint[timestep][i][1] = atom->sp[i][1];
    spint[timestep][i][2] = atom->sp[i][2];
  }
  timestep++;
 }
void ComputeMagnon::compute_array()
{
  //substarct extra step, gives us size of the array
  timestep-=1;

  double dt = update->dt;
  int nomegas = std::ceil((omegamax - omegamin) / omegastep);
  //magnons
  memory->create(S, nomegas, knum * pointsnum, "magnon:S");
  //correletation function
  memory->create(C,atom->nlocal,  timestep, "magnon:C");
  compute_C();
}
void ComputeMagnon::compute_C()
{


}
