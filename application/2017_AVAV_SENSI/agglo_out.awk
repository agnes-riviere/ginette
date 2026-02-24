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
    cmd=sprintf("Treating couple of parameter %d",NR);
    printf("%s\n",cmd);
	namef=sprintf("Sensi_%s_%s.dat",NR,nfile1);
	cmd=sprintf("./remove_first_line.sh %s %s",namef,nfile1);
	system(cmd);
	namef=sprintf("Sensi_%s_%s.dat",NR,nfile2);
	cmd=sprintf("./remove_first_line.sh %s %s",namef,nfile2);
	system(cmd);
	namef=sprintf("Sensi_%s_%s.dat",NR,nfile3);
	cmd=sprintf("./remove_first_line.sh %s %s",namef,nfile3);
	system(cmd);


	cmd=sprintf("./agglo_line.sh %s tmp%s",nametmp,nfile1);
	system(cmd);
	cmd=sprintf("./agglo_line.sh %s tmp%s",nametmp,nfile2);
	system(cmd);
	cmd=sprintf("./agglo_line.sh %s tmp%s",nametmp,nfile3);
	system(cmd);

    cmd=sprintf("./cat_agglo.sh %s %s",name,nametmp);
    system(cmd);

    printf("%d %s\n",NR,$0)>namet;
}
END{ printf("Sensi files agglomerated in %s\n",name);}
