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
  int nq=CODE_TS;//Number of quantiles
  double *quant;
  double *valq;
  float lneigh,theta;
  int nval_min,length;
  double t0,tf;
  for(i=0;i<NKRO_TS;i++){
    pft[i]=NULL;
    pftmp[i]=NULL;
  }

 if (argc != 3){ 
   printf("WRONG COMMAND LINE : lts.%4.2f Cmd_file output_file\n",VERSION_TS);
    exit(1); 
  }
  if ((fp = fopen(argv[1],"r")) == 0) 
    LP_error(fp,"Impossible to open CMD_file %s\n",argv[2]);

  for(i=0;i<NKRO_TS;i++){
    eof=fscanf(fp,"%s\n",cmd_name);
    sprintf(name[i],"%s",cmd_name);
  }

  eof=fscanf(fp,"lneigh=%f\n",&lneigh);
  eof=fscanf(fp,"theta=%f\n",&theta);
  eof=fscanf(fp,"nval_min=%d\n",&nval_min);
  eof=fscanf(fp,"%s\n",cmd_name);
  sprintf(name_mv,"%s",cmd_name);

  fclose(fp);

  fp = fopen(argv[2],"w");
  
  LP_printf(fp,"Opening : %s\n",name_mv);
    fpr=fopen(name_mv,"r");
    if(fpr==NULL)
      LP_error(fp,"Impossible to open %s\n",name_mv);
    pft_amv=TS_create_ts_from_file(name_mv,fp);
    fclose(fpr);
 
    if(pft_amv!=NULL)
      {
	length=TS_length_ts(pft_amv);
	LP_printf(fp,"function t created of length : %d\n",length);
	pft_amv=TS_browse_ft(pft_amv,END_TS);
	tf=pft_amv->t;
	pft_amv=TS_browse_ft(pft_amv,BEGINNING_TS);
	t0=pft_amv->t;
	LP_printf(fp,"\t INPUT DATA MV_AVRG\n lneigh = %f\n theta = %f \n t0 = %f\n tf = %f\n",(double)lneigh,(double)theta,t0,tf);
	pft_mv=TS_moving_average(pft_amv,(double)lneigh,(double)theta,t0,tf,nval_min,fp);
	sprintf(name_mv,"%s/Programmes/LIBS/libts/trunk/test/O2_mv%4.2fd_%3.2f.txt",getenv("HOME"),lneigh,theta);
	LP_printf(fp,"Printing moving average in %s\n",name_mv);
	fpr=fopen(name_mv,"w");
	TS_print_ts(pft_mv,fpr); 
	fclose(fpr);
	LP_printf(stderr,"lts%4.2f: Freeing memory for moving average mv\n",VERSION_TS);
	pft_mv=TS_free_ts(pft_mv,fp);
	
	
	VCN=TS_VCN(pft_amv,(double)lneigh,(double)theta,(double)CODE_TS,(double)CODE_TS,nval_min,fp);
	sprintf(name_mv,"%s/Programmes/LIBS/libts/trunk/test/O2_VCN%4.2fd_%3.2f.txt",getenv("HOME"),lneigh,theta);
	LP_printf(fp,"Printing VCN in %s\n",name_mv);
	fpr=fopen(name_mv,"w");
	TS_print_ts(VCN,fpr); 
	fclose(fpr);
	
	LP_printf(stderr,"lts%4.2f: Freeing memory for moving average amv\n",VERSION_TS);
	pft_amv=TS_free_ts(pft_amv,fp);
	VCN=TS_free_ts(VCN,fp);
	
	LP_printf(stderr,"lts%4.2f: closing log file\n",VERSION_TS);
	fclose(fp);
	LP_printf(stderr,"lts%4.2f: End of Test OK\n",VERSION_TS);
      }
    else{
      LP_error(fp,"function t not created \n");
    }
  return 1;
}
;
