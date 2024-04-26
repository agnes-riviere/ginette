BEGIN{name="awk.ginette";}
{
if (nbzone==1) {
  if(NR==1) {
    printf("%s %s %s\n",min,max,nbpas)>name;
  }
  else{
     print $0>name;
  }
}
if (nbzone==2) {
if (zone==1) {
  if(NR==1) {
    printf("%s %s %s\n",min,max,nbpas)>name;
  }
  else{
     print $0>name;
  }
}
if (zone==2) {
  if(NR==2) {
    printf("%s %s %s\n",min,max,nbpas)>name;
  }
  else{
     print $0>name;
  }
}
}
}
END{printf("zone %s param k on %s zones \n",zone,nbzone);}
