BEGIN{
    name="awk.out";
}
{

}
END{
#First print the obs' file name
    printf("%s%d_t.dat\n",obs,id)>name;
#Then print the sim' file name
    printf("%s%d_%s.dat\n",out,id,id_sim)>name;
    print("Ginette fed\n");
}
