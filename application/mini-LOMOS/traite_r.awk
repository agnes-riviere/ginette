BEGIN{name="awk.ginette";}
{
if (nbzone==1) {
if (zone==1) {
  if(NR==4) {
    printf("%s %s %s\n",min,max,nbpas)>name;
  }
  else{
     print $0>name;
  }
}
}
if (nbzone==2) {
if (zone==1) {
  if(NR==7) {
    printf("%s %s %s\n",min,max,nbpas)>name;
  }
  else{
     print $0>name;
  }
}
if (zone==2) {
  if(NR==8) {
    printf("%s %s %s\n",min,max,nbpas)>name;
  }
  else{
     print $0>name;
  }
}
}
}
END{printf("zone %s param rho on %s zones \n",zone,nbzone);}
