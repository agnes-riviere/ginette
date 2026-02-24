BEGIN{
    name="awk.out";
}
{

}
END{
#First print the obs' file name
    printf("%s%s.dat\n",obs,id)>name;
#Then print the sim' file name
    printf("%s%s_%d.dat\n",out,id,id_sim)>name;
    print("Ginette fed\n");
}
