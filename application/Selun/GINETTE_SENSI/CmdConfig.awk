BEGIN{
    name="awk.out";
    printf("Configuring Ginette's time with nitt=%s\n",nitt);
}
{
    if(NR!=5)
	printf("%s\n",$0)>name;
    else
	printf("nitt=%s         nombre d'iterations en temps                                         C//////////////////////////  CONSTRUCTION MODEL\n",int(nitt))>name;
}
