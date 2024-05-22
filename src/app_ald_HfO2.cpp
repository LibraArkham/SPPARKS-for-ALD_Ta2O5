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

#include "math.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "app_ald_HfO2.h"
#include "solve.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include <iostream>

using namespace SPPARKS_NS;

enum{VACANCY,O,Hf,
    HfX4O,HfX3O,HfX2O,HfXO,HfO,
    HfX3,HfX2,HfX,
    OHf,XOHf};//13 same as DiagAld

#define DELTAEVENT 100000

AppAldHfO2::AppAldHfO2(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  ninteger = 2;
  ndouble = 0;
  delpropensity = 1;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;
  allow_masking = 0;
  

  create_arrays();

  if (narg != 1) error->all(FLERR,"Illegal app_style command");

  cycle = 0;
  pressureOn = 1;
  hello = 1;
  firsttime = 1;
  esites = NULL;
  echeck = NULL;
  events = NULL;
  maxevent = 0;
  firstevent = NULL;

  // reaction lists

  none = ntwo = nthree = 0;
  srate = drate = vrate = NULL;
  spropensity = dpropensity = vpropensity = NULL;
  sinput = soutput = NULL;
  dinput = doutput = NULL;
  vinput = voutput = NULL;
  comneigh = NULL;
  scount = dcount = vcount = NULL;
  sA = dA = vA = NULL;
  scoord = dcoord = vcoord = NULL;
  sexpon = dexpon = vexpon = NULL;
  spresson = dpresson = vpresson = NULL;
}
/* ---------------------------------------------------------------------- */

AppAldHfO2::~AppAldHfO2()
{
  delete [] esites;
  delete [] echeck;
  memory->sfree(events);
  memory->sfree(firstevent);
  memory->sfree(srate);
  memory->sfree(drate);
  memory->sfree(vrate);
  memory->sfree(frate);
  memory->sfree(spropensity);
  memory->sfree(dpropensity);
  memory->sfree(vpropensity);
  memory->sfree(fpropensity);
  memory->sfree(sinput);
  memory->sfree(soutput);
  memory->sfree(dinput);
  memory->sfree(doutput);
  memory->sfree(vinput);
  memory->sfree(voutput);
  memory->sfree(finput);
  memory->sfree(foutput);
  memory->sfree(comneigh);
  memory->sfree(scount);
  memory->sfree(dcount);
  memory->sfree(vcount);
  memory->sfree(fcount);
  memory->sfree(sA);
  memory->sfree(dA);
  memory->sfree(vA);
  memory->sfree(fA);
  memory->sfree(scoord);
  memory->sfree(dcoord);
  memory->sfree(vcoord);
  memory->sfree(fcoord);
  memory->sfree(sexpon);
  memory->sfree(dexpon);
  memory->sfree(vexpon);
  memory->sfree(fexpon);
  memory->sfree(spresson);
  memory->sfree(dpresson);
  memory->sfree(vpresson);
  memory->sfree(fpresson);
}

void AppAldHfO2::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"event") == 0) {
    if (narg < 1) error->all(FLERR,"Illegal event command E1");
    int rstyle = atoi(arg[0]);
    grow_reactions(rstyle);

    if (rstyle == 1) {
      if (narg != 9) error->all(FLERR,"Illegal event arg command");
//type I
      if (strcmp(arg[1],"O") == 0) sinput[none] = O;
      else if (strcmp(arg[1],"VAC") == 0) sinput[none] = VACANCY; 
      else if (strcmp(arg[1],"Hf") == 0) sinput[none] = Hf; 
      else if (strcmp(arg[1],"HfX4O") == 0) sinput[none] = HfX4O;
      else if (strcmp(arg[1],"HfX3O") == 0) sinput[none] = HfX3O;
      else if (strcmp(arg[1],"HfX2O") == 0) sinput[none] = HfX2O;
      else if (strcmp(arg[1],"HfXO") == 0) sinput[none] = HfXO;
      else if (strcmp(arg[1],"HfO") == 0) sinput[none] = HfO;
      else if (strcmp(arg[1],"HfX3") == 0) sinput[none] = HfX3;
      else if (strcmp(arg[1],"HfX2") == 0) sinput[none] = HfX2;
      else if (strcmp(arg[1],"HfX") == 0) sinput[none] = HfX;
      else if (strcmp(arg[1],"OHf") == 0) sinput[none] = OHf;
      else if (strcmp(arg[1],"XOHf") == 0) sinput[none] = XOHf;
   
	
      else error->all(FLERR,"Illegal event arg1 command");

      if (strcmp(arg[2],"O") == 0) soutput[none] = O;
      else if (strcmp(arg[2],"VAC") == 0) soutput[none] = VACANCY; 
      else if (strcmp(arg[2],"Hf") == 0) soutput[none] = Hf; 
      else if (strcmp(arg[2],"HfX4O") == 0) soutput[none] = HfX4O;
      else if (strcmp(arg[2],"HfX3O") == 0) soutput[none] = HfX3O;
      else if (strcmp(arg[2],"HfX2O") == 0) soutput[none] = HfX2O;
      else if (strcmp(arg[2],"HfXO") == 0) soutput[none] = HfXO;
      else if (strcmp(arg[2],"HfO") == 0) soutput[none] = HfO;
      else if (strcmp(arg[2],"HfX3") == 0) soutput[none] = HfX3;
      else if (strcmp(arg[2],"HfX2") == 0) soutput[none] = HfX2;
      else if (strcmp(arg[2],"HfX") == 0) soutput[none] = HfX;
      else if (strcmp(arg[2],"OHf") == 0) soutput[none] = OHf;
      else if (strcmp(arg[2],"XOHf") == 0) soutput[none] = XOHf;
      
      else error->all(FLERR,"Illegal event command E2");
      
      sA[none] = atof(arg[3]);
      if (sA[none] == 0.0) error->warning(FLERR,"Illegal coef during reading command");
      sexpon[none] = atoi(arg[4]);
      srate[none] = atof(arg[5]);
      scoord[none] = atoi(arg[6]);
      spresson[none] = atoi(arg[7]);
      none++;
      
//type II 
    } else if (rstyle == 2) {
      if (narg != 11) error->all(FLERR,"Illegal event command E3");

      if (strcmp(arg[1],"O") == 0) dinput[ntwo][0] = O;
      else if (strcmp(arg[1],"VAC") == 0) dinput[ntwo][0] = VACANCY; 
      else if (strcmp(arg[1],"Hf") == 0) dinput[ntwo][0] = Hf; 
      else if (strcmp(arg[1],"HfX4O") == 0) dinput[ntwo][0] = HfX4O;
      else if (strcmp(arg[1],"HfX3O") == 0) dinput[ntwo][0] = HfX3O;
      else if (strcmp(arg[1],"HfX2O") == 0) dinput[ntwo][0] = HfX2O;
      else if (strcmp(arg[1],"HfXO") == 0) dinput[ntwo][0] = HfXO;
      else if (strcmp(arg[1],"HfO") == 0) dinput[ntwo][0] = HfO;
      else if (strcmp(arg[1],"HfX3") == 0) dinput[ntwo][0] = HfX3;
      else if (strcmp(arg[1],"HfX2") == 0) dinput[ntwo][0] = HfX2;
      else if (strcmp(arg[1],"HfX") == 0) dinput[ntwo][0] = HfX;
      else if (strcmp(arg[1],"OHf") == 0) dinput[ntwo][0] = OHf;
      else if (strcmp(arg[1],"XOHf") == 0) dinput[ntwo][0] = XOHf;

	
      else error->all(FLERR,"Illegal event command E4");
      
      if (strcmp(arg[2],"O") == 0) doutput[ntwo][0] = O;
      else if (strcmp(arg[2],"VAC") == 0) doutput[ntwo][0] = VACANCY; 
      else if (strcmp(arg[2],"Hf") == 0) doutput[ntwo][0] = Hf; 
      else if (strcmp(arg[2],"HfX4O") == 0) doutput[ntwo][0] = HfX4O;
      else if (strcmp(arg[2],"HfX3O") == 0) doutput[ntwo][0] = HfX3O;
      else if (strcmp(arg[2],"HfX2O") == 0) doutput[ntwo][0] = HfX2O;
      else if (strcmp(arg[2],"HfXO") == 0) doutput[ntwo][0] = HfXO;
      else if (strcmp(arg[2],"HfO") == 0) doutput[ntwo][0] = HfO;
      else if (strcmp(arg[2],"HfX3") == 0) doutput[ntwo][0] = HfX3;
      else if (strcmp(arg[2],"HfX2") == 0) doutput[ntwo][0] = HfX2;
      else if (strcmp(arg[2],"HfX") == 0) doutput[ntwo][0] = HfX;
      else if (strcmp(arg[2],"OHf") == 0) doutput[ntwo][0] = OHf;
      else if (strcmp(arg[2],"XOHf") == 0) doutput[ntwo][0] = XOHf;

	
      else error->all(FLERR,"Illegal event command2");

      if (strcmp(arg[3],"O") == 0) dinput[ntwo][1] = O;
      else if (strcmp(arg[3],"VAC") == 0) dinput[ntwo][1] = VACANCY; 
      else if (strcmp(arg[3],"Hf") == 0) dinput[ntwo][1] = Hf; 
      else if (strcmp(arg[3],"HfX4O") == 0) dinput[ntwo][1] = HfX4O;
      else if (strcmp(arg[3],"HfX3O") == 0) dinput[ntwo][1] = HfX3O;
      else if (strcmp(arg[3],"HfX2O") == 0) dinput[ntwo][1] = HfX2O;
      else if (strcmp(arg[3],"HfXO") == 0) dinput[ntwo][1] = HfXO;
      else if (strcmp(arg[3],"HfO") == 0) dinput[ntwo][1] = HfO;
      else if (strcmp(arg[3],"HfX3") == 0) dinput[ntwo][1] = HfX3;
      else if (strcmp(arg[3],"HfX2") == 0) dinput[ntwo][1] = HfX2;
      else if (strcmp(arg[3],"HfX") == 0) dinput[ntwo][1] = HfX;
      else if (strcmp(arg[3],"OHf") == 0) dinput[ntwo][1] = OHf;
      else if (strcmp(arg[3],"XOHf") == 0) dinput[ntwo][1] = XOHf;

	
      else error->all(FLERR,"Illegal event command2");

      if (strcmp(arg[4],"O") == 0) doutput[ntwo][1] = O;
      else if (strcmp(arg[4],"VAC") == 0) doutput[ntwo][1] = VACANCY; 
      else if (strcmp(arg[4],"Hf") == 0) doutput[ntwo][1] = Hf; 
      else if (strcmp(arg[4],"HfX4O") == 0) doutput[ntwo][1] = HfX4O;
      else if (strcmp(arg[4],"HfX3O") == 0) doutput[ntwo][1] = HfX3O;
      else if (strcmp(arg[4],"HfX2O") == 0) doutput[ntwo][1] = HfX2O;
      else if (strcmp(arg[4],"HfXO") == 0) doutput[ntwo][1] = HfXO;
      else if (strcmp(arg[4],"HfO") == 0) doutput[ntwo][1] = HfO;
      else if (strcmp(arg[4],"HfX3") == 0) doutput[ntwo][1] = HfX3;
      else if (strcmp(arg[4],"HfX2") == 0) doutput[ntwo][1] = HfX2;
      else if (strcmp(arg[4],"HfX") == 0) doutput[ntwo][1] = HfX;
      else if (strcmp(arg[4],"OHf") == 0) doutput[ntwo][1] = OHf;
      else if (strcmp(arg[4],"XOHf") == 0) doutput[ntwo][1] = XOHf;
	
      else error->all(FLERR,"Illegal event command2");

      dA[ntwo] = atof(arg[5]);
      dexpon[ntwo] = atoi(arg[6]);
      if (dexpon[ntwo] != 0.0) error->warning(FLERR,"Illegal expon command2");
      drate[ntwo] = atof(arg[7]);
      dcoord[ntwo] = atoi(arg[8]);
      dpresson[ntwo] = atoi(arg[9]);
      ntwo++;
// type III
    }else if (rstyle == 3) {
      if (narg != 11) error->all(FLERR,"Illegal event command31");
 
      if (strcmp(arg[1],"O") == 0) vinput[nthree][0] = O;
      else if (strcmp(arg[1],"VAC") == 0) vinput[nthree][0] = VACANCY; 
      else if (strcmp(arg[1],"Hf") == 0) vinput[nthree][0] = Hf; 
      else if (strcmp(arg[1],"HfX4O") == 0) vinput[nthree][0] = HfX4O;
      else if (strcmp(arg[1],"HfX3O") == 0) vinput[nthree][0] = HfX3O;
      else if (strcmp(arg[1],"HfX2O") == 0) vinput[nthree][0] = HfX2O;
      else if (strcmp(arg[1],"HfXO") == 0) vinput[nthree][0] = HfXO;
      else if (strcmp(arg[1],"HfO") == 0) vinput[nthree][0] = HfO;
      else if (strcmp(arg[1],"HfX3") == 0) vinput[nthree][0] = HfX3;
      else if (strcmp(arg[1],"HfX2") == 0) vinput[nthree][0] = HfX2;
      else if (strcmp(arg[1],"HfX") == 0) vinput[nthree][0] = HfX;
      else if (strcmp(arg[1],"OHf") == 0) vinput[nthree][0] = OHf;
      else if (strcmp(arg[1],"XOHf") == 0) vinput[nthree][0] = XOHf;

      else error->all(FLERR,"Illegal event command32");

      if (strcmp(arg[2],"O") == 0) voutput[nthree][0] = O;
      else if (strcmp(arg[2],"VAC") == 0) voutput[nthree][0]  = VACANCY; 
      else if (strcmp(arg[2],"Hf") == 0) voutput[nthree][0]  = Hf; 
      else if (strcmp(arg[2],"HfX4O") == 0) voutput[nthree][0] = HfX4O;
      else if (strcmp(arg[2],"HfX3O") == 0) voutput[nthree][0] = HfX3O;
      else if (strcmp(arg[2],"HfX2O") == 0) voutput[nthree][0] = HfX2O;
      else if (strcmp(arg[2],"HfXO") == 0) voutput[nthree][0] = HfXO;
      else if (strcmp(arg[2],"HfO") == 0) voutput[nthree][0] = HfO;
      else if (strcmp(arg[2],"HfX3") == 0) voutput[nthree][0] = HfX3;
      else if (strcmp(arg[2],"HfX2") == 0) voutput[nthree][0] = HfX2;
      else if (strcmp(arg[2],"HfX") == 0) voutput[nthree][0] = HfX;
      else if (strcmp(arg[2],"OHf") == 0) voutput[nthree][0] = OHf;
      else if (strcmp(arg[2],"XOHf") == 0) voutput[nthree][0] = XOHf;

      else error->all(FLERR,"Illegal event command33");
      
      if (strcmp(arg[3],"O") == 0) vinput[nthree][1] = O;
      else if (strcmp(arg[3],"VAC") == 0) vinput[nthree][1] = VACANCY; 
      else if (strcmp(arg[3],"Hf") == 0) vinput[nthree][1] = Hf; 
      else if (strcmp(arg[3],"HfX4O") == 0) vinput[nthree][1] = HfX4O;
      else if (strcmp(arg[3],"HfX3O") == 0) vinput[nthree][1] = HfX3O;
      else if (strcmp(arg[3],"HfX2O") == 0) vinput[nthree][1] = HfX2O;
      else if (strcmp(arg[3],"HfXO") == 0) vinput[nthree][1] = HfXO;
      else if (strcmp(arg[3],"HfO") == 0) vinput[nthree][1] = HfO;
      else if (strcmp(arg[3],"HfX3") == 0) vinput[nthree][1] = HfX3;
      else if (strcmp(arg[3],"HfX2") == 0) vinput[nthree][1] = HfX2;
      else if (strcmp(arg[3],"HfX") == 0) vinput[nthree][1] = HfX;
      else if (strcmp(arg[3],"OHf") == 0) vinput[nthree][1] = OHf;
      else if (strcmp(arg[3],"XOHf") == 0) vinput[nthree][1] = XOHf;

    
      else error->all(FLERR,"Illegal event command34");

      if (strcmp(arg[4],"O") == 0) voutput[nthree][1] = O;
      else if (strcmp(arg[4],"VAC") == 0) voutput[nthree][1] = VACANCY; 
      else if (strcmp(arg[4],"Hf") == 0) voutput[nthree][1] = Hf; 
      else if (strcmp(arg[4],"HfX4O") == 0) voutput[nthree][1] = HfX4O;
      else if (strcmp(arg[4],"HfX3O") == 0) voutput[nthree][1] = HfX3O;
      else if (strcmp(arg[4],"HfX2O") == 0) voutput[nthree][1] = HfX2O;
      else if (strcmp(arg[4],"HfXO") == 0) voutput[nthree][1] = HfXO;
      else if (strcmp(arg[4],"HfO") == 0) voutput[nthree][1] = HfO;
      else if (strcmp(arg[4],"HfX3") == 0) voutput[nthree][1] = HfX3;
      else if (strcmp(arg[4],"HfX2") == 0) voutput[nthree][1] = HfX2;
      else if (strcmp(arg[4],"HfX") == 0) voutput[nthree][1] = HfX;
      else if (strcmp(arg[4],"OHf") == 0) voutput[nthree][1] = OHf;
      else if (strcmp(arg[4],"XOHf") == 0) voutput[nthree][1] = XOHf;

      else error->all(FLERR,"Illegal event command35");

      vA[nthree] = atof(arg[5]);
      vexpon[nthree] = atoi(arg[6]);
      if (vexpon[nthree] != 0.0) error->warning(FLERR,"Illegal vexpon command36");
      vrate[nthree] = atof(arg[7]);
      vcoord[nthree] = atoi(arg[8]);
      vpresson[nthree] = atoi(arg[9]);
      nthree++;
    } 
else if (rstyle == 4) {
      if (narg != 11) error->all(FLERR,"Illegal event command31"); 
      if (strcmp(arg[1],"O") == 0) finput[nfour][0] = O;
      else if (strcmp(arg[1],"VAC") == 0) finput[nfour][0] = VACANCY; 
      else if (strcmp(arg[1],"Hf") == 0) finput[nfour][0] = Hf; 
      else if (strcmp(arg[1],"HfX4O") == 0) finput[nfour][0] = HfX4O;
      else if (strcmp(arg[1],"HfX3O") == 0) finput[nfour][0] = HfX3O;
      else if (strcmp(arg[1],"HfX2O") == 0) finput[nfour][0] = HfX2O;
      else if (strcmp(arg[1],"HfXO") == 0) finput[nfour][0] = HfXO;
      else if (strcmp(arg[1],"HfO") == 0) finput[nfour][0] = HfO;
      else if (strcmp(arg[1],"HfX3") == 0) finput[nfour][0] = HfX3;
      else if (strcmp(arg[1],"HfX2") == 0) finput[nfour][0] = HfX2;
      else if (strcmp(arg[1],"HfX") == 0) finput[nfour][0] = HfX;
      else if (strcmp(arg[1],"OHf") == 0) finput[nfour][0] = OHf;
      else if (strcmp(arg[1],"XOHf") == 0) finput[nfour][0] = XOHf;

      else error->all(FLERR,"Illegal event command32");

      if (strcmp(arg[2],"O") == 0) foutput[nfour][0] = O;
      else if (strcmp(arg[2],"VAC") == 0) foutput[nfour][0]  = VACANCY; 
      else if (strcmp(arg[2],"Hf") == 0) foutput[nfour][0]  = Hf; 
      else if (strcmp(arg[2],"HfX4O") == 0) foutput[nfour][0] = HfX4O;
      else if (strcmp(arg[2],"HfX3O") == 0) foutput[nfour][0] = HfX3O;
      else if (strcmp(arg[2],"HfX2O") == 0) foutput[nfour][0] = HfX2O;
      else if (strcmp(arg[2],"HfXO") == 0) foutput[nfour][0] = HfXO;
      else if (strcmp(arg[2],"HfO") == 0) foutput[nfour][0] = HfO;
      else if (strcmp(arg[2],"HfX3") == 0) foutput[nfour][0] = HfX3;
      else if (strcmp(arg[2],"HfX2") == 0) foutput[nfour][0] = HfX2;
      else if (strcmp(arg[2],"HfX") == 0) foutput[nfour][0] = HfX;
      else if (strcmp(arg[2],"OHf") == 0) foutput[nfour][0] = OHf;
      else if (strcmp(arg[2],"XOHf") == 0) foutput[nfour][0] = XOHf;

      else error->all(FLERR,"Illegal event command33");
      
      if (strcmp(arg[3],"O") == 0) finput[nfour][1] = O;
      else if (strcmp(arg[3],"VAC") == 0) finput[nfour][1] = VACANCY; 
      else if (strcmp(arg[3],"Hf") == 0) finput[nfour][1] = Hf; 
      else if (strcmp(arg[3],"HfX4O") == 0) finput[nfour][1] = HfX4O;
      else if (strcmp(arg[3],"HfX3O") == 0) finput[nfour][1] = HfX3O;
      else if (strcmp(arg[3],"HfX2O") == 0) finput[nfour][1] = HfX2O;
      else if (strcmp(arg[3],"HfXO") == 0) finput[nfour][1] = HfXO;
      else if (strcmp(arg[3],"HfO") == 0) finput[nfour][1] = HfO;
      else if (strcmp(arg[3],"HfX3") == 0) finput[nfour][1] = HfX3;
      else if (strcmp(arg[3],"HfX2") == 0) finput[nfour][1] = HfX2;
      else if (strcmp(arg[3],"HfX") == 0) finput[nfour][1] = HfX;
      else if (strcmp(arg[3],"OHf") == 0) finput[nfour][1] = OHf;
      else if (strcmp(arg[3],"XOHf") == 0) finput[nfour][1] = XOHf;

    
      else error->all(FLERR,"Illegal event command34");

      if (strcmp(arg[4],"O") == 0) foutput[nfour][1] = O;
      else if (strcmp(arg[4],"VAC") == 0) foutput[nfour][1] = VACANCY; 
      else if (strcmp(arg[4],"Hf") == 0) foutput[nfour][1] = Hf; 
      else if (strcmp(arg[4],"HfX4O") == 0) foutput[nfour][1] = HfX4O;
      else if (strcmp(arg[4],"HfX3O") == 0) foutput[nfour][1] = HfX3O;
      else if (strcmp(arg[4],"HfX2O") == 0) foutput[nfour][1] = HfX2O;
      else if (strcmp(arg[4],"HfXO") == 0) foutput[nfour][1] = HfXO;
      else if (strcmp(arg[4],"HfO") == 0) foutput[nfour][1] = HfO;
      else if (strcmp(arg[4],"HfX3") == 0) foutput[nfour][1] = HfX3;
      else if (strcmp(arg[4],"HfX2") == 0) foutput[nfour][1] = HfX2;
      else if (strcmp(arg[4],"HfX") == 0) foutput[nfour][1] = HfX;
      else if (strcmp(arg[4],"OHf") == 0) foutput[nfour][1] = OHf;
      else if (strcmp(arg[4],"XOHf") == 0) foutput[nfour][1] = XOHf;

      else error->all(FLERR,"Illegal event command35");

      fA[nfour] = atof(arg[5]);
      fexpon[nfour] = atoi(arg[6]);
      if (fexpon[nfour] != 0.0) error->warning(FLERR,"Illegal vexpon command36");
      frate[nfour] = atof(arg[7]);
      fcoord[nfour] = atoi(arg[8]);
      fpresson[nfour] = atoi(arg[9]);
      nfour++;
    }
  }
  else if (strcmp(command,"pulse_time") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal pulse time");
      T1 = atof(arg[0]);
      T3 = atof(arg[1]);
      std::cout<<"T1 = "<<T1<<" T3 = "<<T3<<std::endl;
  }
  else if (strcmp(command,"purge_time") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal purge time");
      T2 = atof(arg[0]);
      T4 = atof(arg[1]);
      std::cout<<"T2 = "<<T2<<" T4 = "<<T4<<std::endl;
  }else error->all(FLERR,"Unrecognized command38");
}
/* ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppAldHfO2::grow_app()
{
  element = iarray[0];
  coord = iarray[1];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppAldHfO2::init_app()
{
  if (firsttime) {
    firsttime = 0;

    echeck = new int[nlocal];
    firstevent = (int *) memory->smalloc(nlocal*sizeof(int),"app:firstevent");
    //comneigh was defined to avoid double counting of common neighbor in site_propensity
    comneigh = memory->grow(comneigh,12*maxneigh,2,"app/ald/hfo2:comneigh");
    //comneigh = memory->grow_2d_double_array(comneigh,12*maxneigh,2,"app/ald:comneigh");
    // esites must be large enough for 3 sites and their 1st neighbors
    
    // esites = (int *) memory->smalloc(12*maxneigh*sizeof(int),"app:esites");
    esites = (int *) memory->smalloc(nlocal*sizeof(int),"app:esites");
    //esites = new int[12*maxneigh]; 
  }
  // site validity

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (coord[i] < -1 || coord[i] > 8) flag = 1;
    if (element[i] < VACANCY) flag = 1;
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}
/* ---------------------------------------------------------------------- */

void AppAldHfO2::setup_app()
{
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;
  

  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;


  if (temperature == 0.0)
    error->all(FLERR,"Temperature cannot be 0.0 for app_ald");
  for (int m = 0; m < none; m++) {
    spropensity[m] = sA[m]*pow(temperature,sexpon[m])*exp(-srate[m]/temperature);
    scount[m] = 0;
  if (spropensity[m] == 0.0) error->warning(FLERR," spropensity cannot be 0.0 for app_ald");
  }
  for (int m = 0; m < ntwo; m++) {
    dpropensity[m] = dA[m]*pow(temperature,dexpon[m])*exp(-drate[m]/temperature);
    dcount[m] = 0;
  if (dpropensity[m] == 0.0) error->warning(FLERR,"dpropensity cannot be 0.0 for app_ald");
  }
  for (int m = 0; m < nthree; m++) {
    vpropensity[m] = vA[m]*pow(temperature,vexpon[m])*exp(-vrate[m]/temperature);
    vcount[m] = 0;
  if (vpropensity[m] == 0.0) error->warning(FLERR,"vpropensity cannot be 0.0 for app_ald");
  }
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppAldHfO2::site_energy(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppAldHfO2::site_propensity(int i)
{
   int j,k,m;


  clear_events(i);
  // std::cout<<"clear_events"<<std::endl;
  double proball = 0.0;


  //type I, check species of sites and consider possible events

  // count_coordO was added here to prevent adsorption of metal precursor
  // on the low coordinate oxygen at sublayer

  for (m = 0; m < none; m++) {
    if (element[i] != sinput[m]) continue;
    if ((coord[i] == scoord[m] || scoord[m] == -1) && (spresson[m] == pressureOn || spresson[m] == 0) && (coord[i]<numneigh[i])) {
      // std::cout<<"coord[i]"<<coord[i]<<std::endl;
      // std::cout<<"numneigh[i]"<<numneigh[i]<<std::endl;
	    add_event(i,1,m,spropensity[m],-1,-1,-1);
	    proball += spropensity[m];
    }
  }
  // std::cout<<"test1"<<std::endl;
  // type II, check species of sites and second neighbor 
  // comneigh variable used to avoid double counting,
  // we consider more than one event between sites, therefore we need 2d comneigh array.

  int nextneib = 1;
  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
       for (int kk = 0; kk < numneigh[j]; kk++) {
            k = neighbor[j][kk];
            if (i == k) continue;
            for (m = 0; m < ntwo; m++) {
	    if ((element[i] == dinput[m][0] && element[k] == dinput[m][1]) && (dpresson[m] == pressureOn || dpresson[m] == 0) && (coord[i] == dcoord[m] || dcoord[m] == -1)) {
  	      comevent = 1;
  	      for (int ii = 0; ii < nextneib; ii++) {
  	      if ( comneigh[ii][0] == k && comneigh[ii][1] == dpropensity[m]) comevent = 0;
  	      }
  	      if (comevent){	
                add_event(i,2,m,dpropensity[m],-1,k,-1);
                proball += dpropensity[m];
  	        comneigh[nextneib][0] = k;
                comneigh[nextneib][1] = dpropensity[m];
  	        nextneib++;
  	    } 
          }
        }
      }
    }
    for (m = 0; m < nextneib; m++) comneigh[m][0] = comneigh[m][1] = 0;
    // std::cout<<"test2"<<std::endl;
  //type III, check species of sites and first neighbour 
  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
      for (m =0; m < nthree; m++) {
      // std::cout<<"test3!!!"<<std::endl;
      if (element[i] == vinput[m][0] && element[j] == vinput[m][1] && (coord[i] == vcoord[m] || vcoord[m] == -1 ) && (vpresson[m] == pressureOn || vpresson[m] == 0)) {
        add_event(i,3,m,vpropensity[m],j,-1,-1);
        proball += vpropensity[m];
      }
    }
  }
  // type IIII, check species of sites and third neighbor 
  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
      for (int kk = 0; kk < numneigh[j]; kk++) {
        k = neighbor[j][kk];
          for (int gg = 0; gg < numneigh[k]; gg++) {
            int g = neighbor[k][gg];
            for (m = 0; m < nfour; m++) {
              if (element[i] == finput[m][0] && element[g] == finput[m][1] && (coord[i] == fcoord[m] || fcoord[m] == -1 ) && (fpresson[m] == pressureOn || fpresson[m] == 0)) {
                add_event(i,4,m,fpropensity[m],-1,-1,g);
                proball += fpropensity[m];
              }
            }
        }
    }
  }
  // vacant event:5
  add_event(i,5,-1,0.1,-1,-1,-1);
  proball += 0.1;
  // std::cout<<"test3"<<std::endl;
  return proball;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */
//修改
void AppAldHfO2::site_event(int i, class RandomPark *random)
{
  int j,k,g,m,n,mm,jj;
  int elcoord = element[i];

  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
  }

  int rstyle = events[ievent].style;
  int which = events[ievent].which;
  j = events[ievent].jpartner;
  k = events[ievent].kpartner;
  g = events[ievent].gpartner;


  if (rstyle == 1) {
    element[i] = soutput[which];
    scount[which]++;
    // std::cout<<"event1"<<std::endl;
    } 
  else if (rstyle == 2 && j == -1) {
    element[i] = doutput[which][0];
    element[k] = doutput[which][1];
    dcount[which]++;
    // std::cout<<"event2"<<std::endl;
    }
  else if (rstyle == 3 && k == -1) {
    element[i] = voutput[which][0];
    element[j] = voutput[which][1];
    vcount[which]++;
    // std::cout<<"event3"<<std::endl;
    }
  else if (rstyle == 4) {
    element[i] = foutput[which][0];
    element[g] = foutput[which][1];
    fcount[which]++;
  }
  else if (rstyle == 5) {}
  else { error->all(FLERR,"Illegal execution event"); }

  update_coord(elcoord,i,j,k,which);

  // sequence of ALD, 
  // 1 is metal pulse, 3 purge, 2 oxygen pulse.
  if ((cycle+T1) > time ) {pressureOn = 1;
  // std::cout<<"metal pulse"<<std::endl;
  }
  else if ((cycle+T1)<= time && time < (cycle+T1+T2)) {pressureOn = 3;}
  else if ((cycle+T1+T2) <= time && time < (cycle+T1+T2+T3)) {pressureOn = 2;}
  else if ((cycle+T1+T2+T3) <= time && time < (cycle+T1+T2+T3+T4)) {pressureOn = 3;}
  else {cycle += T1+T2+T3+T4; }
  // std::cout<<"pressureOn:"<<pressureOn<<std::endl;
  // std::cout<<"cycle:"<<cycle<<std::endl;
  int nsites = 0;
  int isite = i2site[i];
  // std::cout<<"site event"<<std::endl;
  // 全局更新
  // propensity[isite] = site_propensity(i);
  // esites[nsites++] = isite;
  // echeck[isite] = 1;
  // for (int n = 0; n < nlocal; n++) {
  //   isite = i2site[n];
  //   propensity[isite] = site_propensity(n);
  //   esites[nsites++] = isite;
  // }
  //更新影响的site
  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;
  for (int n = 0; n < numneigh[i]; n++) {
    m = neighbor[i][n];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
    propensity[isite] = site_propensity(m);
    esites[nsites++] = isite;
    // std::cout<<"nsites:"<<nsites<<std::endl;
    echeck[isite] = 1;
    }
    for (int jj = 0; jj< numneigh[m];jj++) {
    mm = neighbor[m][jj];
    isite = i2site[mm];
        if (isite >= 0 && echeck[isite] == 0) {
            propensity[isite] = site_propensity(mm);
            esites[nsites++] = isite;
            // std::cout<<"nsites:"<<nsites<<std::endl;
            echeck[isite] = 1;
        }
        for (int ss = 0; ss< numneigh[mm]; ss++) {
            int s = neighbor[mm][ss];
            isite = i2site[s];
            if (isite >= 0 && echeck[isite] == 0) {
              propensity[isite] = site_propensity(s);
              esites[nsites++] = isite;
              // std::cout<<"nsites:"<<nsites<<std::endl;
              echeck[isite] = 1;
            }
            for (int gg = 0; gg< numneigh[s]; gg++) {
              int g = neighbor[s][gg];
              isite = i2site[g];
              if (isite >= 0 && echeck[isite] == 0) {
                propensity[isite] = site_propensity(g);
                esites[nsites++] = isite;
                // std::cout<<"nsites:"<<nsites<<std::endl;
                echeck[isite] = 1;
              }
            }
        }
    }
  }
  solve->update(nsites,esites,propensity);

   // clear echeck array

  for (m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppAldHfO2::clear_events(int i)
{
  int next;
  int index = firstevent[i];
  while (index >= 0) {
    next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    nevents--;
    index = next;
  }
  firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
   add an event to list for site I
   event = exchange with site J with probability = propensity
------------------------------------------------------------------------- */

void AppAldHfO2::add_event(int i, int rstyle, int which, double propensity,
			  int jpartner, int kpartner, int gpartner)
{
  if (nevents == maxevent) {
    maxevent += DELTAEVENT;
    events = 
      (Event *) memory->srealloc(events,maxevent*sizeof(Event),"app:events");
    for (int m = nevents; m < maxevent; m++) events[m].next = m+1;
    freeevent = nevents;
  }

  int next = events[freeevent].next;

  events[freeevent].style = rstyle;
  events[freeevent].which = which;
  events[freeevent].jpartner = jpartner;
  events[freeevent].kpartner = kpartner;
  events[freeevent].gpartner = gpartner;
  events[freeevent].propensity = propensity;

  if ( propensity == 0 ) error->all(FLERR,"propensity in add_event wrong app ald");
  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ----------------------------------------------------------------------
   grow list of stored reactions for single and double
------------------------------------------------------------------------- */

void AppAldHfO2::grow_reactions(int rstyle)
{
  if (rstyle == 1) {
    int n = none + 1;
    srate = (double *) 
      memory->srealloc(srate,n*sizeof(double),"app/ald/hfo2:srate");
    spropensity = (double *) 
      memory->srealloc(spropensity,n*sizeof(double),"app/ald/hfo2:spropensity");
    sinput = (int *) 
      memory->srealloc(sinput,n*sizeof(int),"app/ald/hfo2:sinput");
    soutput = (int *) 
      memory->srealloc(soutput,n*sizeof(int),"app/ald/hfo2:soutput");
    scount = (int *) 
      memory->srealloc(scount,n*sizeof(int),"app/ald/hfo2:scount");
    sA = (double *) 
      memory->srealloc(sA,n*sizeof(double),"app/ald/hfo2:sA");
    sexpon = (int *) 
      memory->srealloc(sexpon,n*sizeof(int),"app/ald/hfo2:sexpon");
    scoord = (int *) 
      memory->srealloc(scoord,n*sizeof(int),"app/ald/hfo2:scoord");
    spresson = (int *) 
      memory->srealloc(spresson,n*sizeof(int),"app/ald/hfo2:spresson");

  } else if (rstyle == 2) {
    int n = ntwo + 1;
    drate = (double *) 
      memory->srealloc(drate,n*sizeof(double),"app/ald/hfo2:drate");
    dpropensity = (double *) 
      memory->srealloc(dpropensity,n*sizeof(double),"app/ald/hfo2:dpropensity");
    dinput = memory->grow(dinput,n,2,"app/ald/hfo2:dinput");
    doutput = memory->grow(doutput,n,2,"app/ald/hfo2:doutput");
    dcount = (int *) 
      memory->srealloc(dcount,n*sizeof(int),"app/ald/hfo2:dcount");
    dA = (double *) 
      memory->srealloc(dA,n*sizeof(double),"app/ald/hfo2:dA");
    dexpon = (int *) 
      memory->srealloc(dexpon,n*sizeof(int),"app/ald/hfo2:dexpon");
    dcoord = (int *) 
      memory->srealloc(dcoord,n*sizeof(int),"app/ald/hfo2:dcoord");
    dpresson = (int *) 
      memory->srealloc(dpresson,n*sizeof(int),"app/ald/hfo2:dpresson");

  } else if (rstyle == 3) {
    int n = nthree + 1;
    vrate = (double *)
      memory->srealloc(vrate,n*sizeof(double),"app/ald/hfo2:vrate");
    vpropensity = (double *)
      memory->srealloc(vpropensity,n*sizeof(double),"app/ald/hfo2:vpropensity");
    vinput = memory->grow(vinput,n,2,"app/ald/hfo2:vinput");
    voutput = memory->grow(voutput,n,2,"app/ald/hfo2:voutput");
    vcount = (int *)
      memory->srealloc(vcount,n*sizeof(int),"app/ald/hfo2:vcount");
    vA = (double *)
      memory->srealloc(vA,n*sizeof(double),"app/ald/hfo2:vA");
    vexpon = (int *)
      memory->srealloc(vexpon,n*sizeof(int),"app/ald/hfo2:vexpon");
    vcoord = (int *)
      memory->srealloc(vcoord,n*sizeof(int),"app/ald/hfo2:vcoord");
    vpresson = (int *)
      memory->srealloc(vpresson,n*sizeof(int),"app/ald/hfo2:vpresson");
  }
  else if (rstyle == 4) {
    int n = nfour + 1;
    frate = (double *)
      memory->srealloc(frate,n*sizeof(double),"app/ald/hfo2:frate");
    fpropensity = (double *)
      memory->srealloc(fpropensity,n*sizeof(double),"app/ald/hfo2:fpropensity");
    finput = memory->grow(finput,n,2,"app/ald/hfo2:finput");
    foutput = memory->grow(foutput,n,2,"app/ald/hfo2:foutput");
    fcount = (int *)
      memory->srealloc(fcount,n*sizeof(int),"app/ald/hfo2:fcount");
    fA = (double *)
      memory->srealloc(fA,n*sizeof(double),"app/ald/hfo2:fA");
    fexpon = (int *)
      memory->srealloc(fexpon,n*sizeof(int),"app/ald/hfo2:fexpon");
    fcoord = (int *)
      memory->srealloc(fcoord,n*sizeof(int),"app/ald/hfo2:fcoord");
    fpresson = (int *)
      memory->srealloc(fpresson,n*sizeof(int),"app/ald/hfo2:fpresson");
  }
}
/* ----------------------------------------------------------------------
   update c.n. for Hf and O, put and remove mask for relative sites
------------------------------------------------------------------------- */
void AppAldHfO2::update_coord(int elcoord, int i, int j, int k, int which)
{ 
  //type1
  if (j==-1 && k==-1){
    if (elcoord == O  && element[i] == HfX4O ) {
      coord[i]++;
    }
    else if (elcoord == HfX4O && element[i] == O ){
      coord[i]--;
    }
    else if (elcoord == Hf  && element[i] == OHf ) {
      coord[i]++;
    }
    else if (elcoord == OHf && element[i] == Hf ) {
      coord[i]--;
    }
  }
  //type2
  else if (j!=-1 && k==-1) {
    count_coord(i);
    count_coord(j);
  }
  //type3
  else if (j==-1 && k!=-1) {
    // Hf densification
  count_coord(i);
  count_coord(k);
  }
}

void AppAldHfO2::count_coord(int i) {
  coord[i] = 0;
  for (int s = 0; s < numneigh[i]; s++) {
    int nn = neighbor[i][s];
    if (element[nn] != VACANCY) coord[i]++;
  }
  // if (element[i] == HfX3) 
    // coord[i] += 3;
  // else if (element[i] == HfX2) 
    // coord[i] += 2;
  // else if (element[i] == HfX)
    // coord[i] += 1;
  if (element[i] == OHf)
    coord[i] += 1;
  // else if (element[i] == HfX3O)
    // coord[i] += 1;
  // else if (element[i] == HfX2O)
    // coord[i] += 1;
  // else if (element[i] == HfXO)
    // coord[i] += 1;
  else if (element[i] == HfO)
    coord[i] += 1;
}