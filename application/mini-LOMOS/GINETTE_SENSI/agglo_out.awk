BEGIN{
    name="awk_agglo.out";
    namet="tmp.out";
    #nametmp="tmp";
    char *namef;
    char *tmp;
}
{
    nametmp=sprintf("tmp_%d",NR);
   # cmd=sprintf("rm %s",nametmp);
   # system(cmd);
    printf("%d",NR)>nametmp;
    cmd=sprintf("Treating line %d: %s %s %s %s %s",NR,$1,$2,$3,$4,$5);
    printf("%s\n",cmd);
    for (i=1;i<=nsensors;i++){
	namef=sprintf("Sensi_temperature_%s_%d.dat",NR,i);
	cmd=sprintf("./remove_first_line.sh %s %i",namef,i);
	system(cmd);
    }
    for (i=1;i<=nsensors;i++){
	cmd=sprintf("./agglo_line.sh %s tmp%i",nametmp,i);
	system(cmd);
    }

    cmd=sprintf("./cat_agglo.sh %s %s",name,nametmp);
    system(cmd);

    printf("%d %s\n",NR,$0)>namet;
}
END{ printf("temperature files agglomerated in %s\n",name);}
