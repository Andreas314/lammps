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
#include "lattice.h"
#include "modify.h"
#include "memory.h"
#include "pair.h"
#include "update.h"

#include <cstring>
#include <complex>
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
  lmp->memory->destroy(C);
  lmp->memory->destroy(C_omega_real);
  lmp->memory->destroy(C_omega_imag);
  lmp->memory->destroy(S_real);
  lmp->memory->destroy(S_imag);
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

  nomegas = std::ceil((omegamax - omegamin) / omegastep);
  //magnons
  memory->create(S_real, 3, nomegas, knum * pointsnum, "magnon:S_real");
  memory->create(S_imag, 3, nomegas, knum * pointsnum, "magnon:S_imag");
  //correletation function
  memory->create(C, 3 * timestep / 4, "magnon:C");
  memory->create(C_omega_real, nomegas, "magnon:C_omega_real");
  memory->create(C_omega_imag, nomegas, "magnon:C_omega_imag");
  //Actual calculation
  calculate_reciprocal();
  for (int omega = 0; omega < nomegas; omega++)
  {
    for (int k = 0; k < knum * pointsnum; k++)
    {
      for (int comp = 0; comp < 3; comp++)
      {
        S_imag[comp][omega][k] = 0.0;
        S_real[comp][omega][k] = 0.0;
        calculate_S_entry(omega, k, comp);
      }
    }
  }
  //MPI reduction over ranks
  MPI_Comm world = lmp->world;
  MPI_Allreduce(MPI_IN_PLACE, S_real, 
		3 * nomegas * knum * pointsnum,
                MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, S_imag, 
		3 * nomegas * knum * pointsnum,
                MPI_DOUBLE, MPI_SUM, world);
  //only rank 0 writes to the file
  if (comm->me == 0) write_result();
}
void ComputeMagnon::write_result()
{

}
void ComputeMagnon::calculate_S_entry(int omega, int k, int comp)
{
  std::complex<double> I{0.0, 1.0};
  for (int i = 0; i < atom->nlocal; i++)
  {
    for (int j = 0; j < atom->nmax; j++)
    {
      compute_C(i, j, comp);
      transform_C();
      double* qnounits = kvecs[k];
      double x, q, sum = 0;
      for (int xx = 0; xx < 3; xx++)
      {
        x = atom->x[i][xx] - atom->x[j][xx];
	q = qnounits[0] * b1[xx];
	q += qnounits[1] * b2[xx];
	q += qnounits[2] * b3[xx];
	sum += x * 1;
      }
      std::complex<double> z = std::exp(I * q);
      S_real[comp][omega][k] += z.real() * C_omega_real[omega] - z.imag() * C_omega_imag[omega];
      S_imag[comp][omega][k] += z.real() * C_omega_imag[omega] + z.imag() * C_omega_real[omega];
    }
  }
  S_real[comp][omega][k] /= (std::sqrt(2 * M_PI) * atom->nlocal);
}

void ComputeMagnon::transform_C()
{
  std::complex<double> I{0.0, 1.0};
  for (int j = 0; j < nomegas; j++)
  {
    double omega = omegamin + omegastep * j;
    C_omega_real[j] = 0.0;
    C_omega_imag[j] = 0.0;
    for (int t = 0; t < 3 * timestep / 4; t++)
    {
      std::complex<double> z = std::exp(I *(update->dt * t) * omega ) * C[t] * update->dt;
      C_omega_real[j] += z.real();
      C_omega_imag[j] += z.imag();
    }
  }
}

void ComputeMagnon::compute_C(int i, int j, int comp)
{
  for (int tau = 0; tau < 3 * timestep / 4; tau++)
  {
    double sisj = 0.0;
    double si = 0.0;
    double sj = 0.0;
    for (int t = tau; t < timestep; t++)
    {
      sisj += spint[t + tau][i][comp] * spint[t][j][comp];
      si += spint[t + tau][i][comp];
      sj += spint[t][j][comp];
    }
    sisj /=  (double)(timestep - tau);
    si /= (double)(timestep - tau);
    sj /= (double)(timestep - tau);
    C[tau] = sisj - si * sj;
  }
}
void ComputeMagnon::calculate_reciprocal()
{
    //TODO:Cooking of ChatGPT, check if correct
    double *a1 = domain->lattice->a1;
    double *a2 = domain->lattice->a2;
    double *a3 = domain->lattice->a3;
    double cross23[3], cross31[3], cross12[3];

    // a2 x a3
    cross23[0] = a2[1]*a3[2] - a2[2]*a3[1];
    cross23[1] = a2[2]*a3[0] - a2[0]*a3[2];
    cross23[2] = a2[0]*a3[1] - a2[1]*a3[0];

    // a3 x a1
    cross31[0] = a3[1]*a1[2] - a3[2]*a1[1];
    cross31[1] = a3[2]*a1[0] - a3[0]*a1[2];
    cross31[2] = a3[0]*a1[1] - a3[1]*a1[0];

    // a1 x a2
    cross12[0] = a1[1]*a2[2] - a1[2]*a2[1];
    cross12[1] = a1[2]*a2[0] - a1[0]*a2[2];
    cross12[2] = a1[0]*a2[1] - a1[1]*a2[0];

    // Volume
    double V =
        a1[0]*cross23[0] +
        a1[1]*cross23[1] +
        a1[2]*cross23[2];

    double factor = 2.0 * M_PI / V;

    for (int i=0;i<3;i++) {
        b1[i] = factor * cross23[i];
        b2[i] = factor * cross31[i];
        b3[i] = factor * cross12[i];
    }
}

