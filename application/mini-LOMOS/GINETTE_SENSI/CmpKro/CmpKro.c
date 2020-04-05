#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <strings.h>
//#include <libprint.h>
#include "libprint.h"
#include "time_series.h"


int main(int argc,char **argv) {
  FILE *fp,*fpr;
  char cmd_name[STRING_LENGTH_LP];
  char name[NKRO_TS][STRING_LENGTH_LP];
  char name_mv[STRING_LENGTH_LP];
  int eof=CODE_TS;
  s_ft *pft[NKRO_TS];
  s_ft *pft_amv;
  s_ft *pft_mv,*VCN;
  s_ft *pftmp[NKRO_TS];
  int i;
  double yr;
  double val;
  int nclass=CODE_TS;//Number of class for the distribution
  //int nq=CODE_TS;//Number of quantiles
  //double *quant;
  int nq=7;//Number of quantiles
  double quant[7] = {5,10, 25, 50, 75, 90,95};
  double *valq;
  float lneigh,theta;
  int nval_min,length;
  double t0,tf;
  double *useless;

  fp = fopen(argv[2],"w");

  if (argc != 3){ 
    printf("WRONG COMMAND LINE : ccp Cmd_file output_file\n");
    exit(1); 
  }
  else printf("Entering CmpKro\n");
  if ((fpr = fopen(argv[1],"r")) == 0) 
    LP_error(fp,"Impossible to open CMD_file %s\n",argv[2]);

  for (i=0;i<NKRO_TS;i++){
    eof=fscanf(fpr,"%s\n",name[i]);
    printf("Opening %s\n",name[i]);
    pft[i]=TS_create_ts_from_file(name[i],fp);
  }
  TS_cmp_kro(pft[SIM_TS],pft[OBS_TS],fp);

  
  for (i=0;i<NKRO_TS;i++)
    pft[i]=TS_free_ts(pft[i],fp);
  fclose(fpr);
  fclose(fp);
  
  return 1;
}
;
