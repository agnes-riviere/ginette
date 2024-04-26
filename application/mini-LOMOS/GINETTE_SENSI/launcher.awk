# appele par HZ1D.sh

BEGIN{
  name="awk.out";
}
{
 
    printf("Treating parameter set %d\n",NR);
    system("./treat_onePset.sh "$1" "$2" "$3" "$4" "$5" "$6" "$7""); 
  
}
END{ printf("ALL PARAMETER SETS TREATED BY GINETTE\n");}
