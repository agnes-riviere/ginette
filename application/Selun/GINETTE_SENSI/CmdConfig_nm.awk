BEGIN{ name="awk.out";
    printf("Configuring Ginette's number cell with nm=%s\n",nm);
 printf("Configuring Ginette's with az = %s\n",L);
}
{
  if((NR==8)||(NR==14)||(NR==16)||(NR==18))	{
    if(NR==8)
	printf("az=%4.3e       profondeur                     C//////////////////////////  CONSTRUCTION MODEL\n",L)>name;
    if(NR==14)
	printf("repbot=%3.2e       m	altitude bottom  C//////////////////////////  CONSTRUCTION MODEL\n",-L)>name;
    if(NR==16)
	printf("nmi=%05d       		nombre de mailles actives             C//////////////////////////  CONSTRUCTION MODEL\n",nm)>name;
    if(NR==18)
	printf("nri=%05d       		nombre de lignes                        C//////////////////////////  CONSTRUCTION MODEL\n",nm)>name;
	}
    else{
	printf("%s\n",$0)>name;
	}
}
