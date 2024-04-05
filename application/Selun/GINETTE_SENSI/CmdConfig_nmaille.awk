BEGIN{ name="awk.out";
    printf("Configuring Ginette's obs cell with nmaille= %s %s  %s \n",nmaille1,nmaille2,nmaille3);
}
{
  if((NR==68)||(NR==69)||(NR==70)){
    if(NR==68)
	printf("nmaille1=%05d       		 numero de la maille enregistrement pression temps                       C//////////////////////////  C FICHIERS DE SORTIE\n",nmaille1)>name;
    if(NR==69)
	printf("nmaille2=%05d       		 numero de la maille enregistrement pression temps                       C//////////////////////////  C FICHIERS DE SORTIE\n",nmaille2)>name;
    if(NR==70)
	printf("nmaille3=%05d        numero de la maille enregistrement pression temps                       C//////////////////////////  C FICHIERS DE SORTIE\n",nmaille3)>name;
	}
    else{
	printf("%s\n",$0)>name;
	}
}
