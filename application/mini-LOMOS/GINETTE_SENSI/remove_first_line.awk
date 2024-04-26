BEGIN{ name="awk.out";}
{
    if(NR>1)
	printf("%s\n",$0)>name;
}
