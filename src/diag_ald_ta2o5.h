/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   ALD application of HfO2 was developed by:
   mahdi shirazi: m.shirazi@tue.nl, TU/e department of applied physics,
   Simon D. Elliott: simon.elliott@schrodinger.com, Schrodinger Materials Science.
   This application is a part of SPPARKS and authors retian the above term.
   See the manual-app-ald and examples folders for more information.
------------------------------------------------------------------------- */

#ifdef DIAG_CLASS
DiagStyle(ald/ta2o5,DiagAldta2o5)

#else

#ifndef SPK_DIAG_ALD_Ta2o5_H
#define SPK_DIAG_ALD_Ta2o5_H

#include "diag.h"

namespace SPPARKS_NS {

class DiagAldta2o5 : public Diag {
 public:
  DiagAldta2o5(class SPPARKS *, int, char **);
  ~DiagAldta2o5();
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);

 private:
  class AppAldTa2o5 *appaldta2o5;
  int nlist;
  char **list;
  int *which,*index,*ivector;
  int siteflag;
};

}

#endif
#endif
