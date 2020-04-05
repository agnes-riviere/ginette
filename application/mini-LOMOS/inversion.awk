function digits(x)
{
  char *digit;
  if(x>9)
    digit=sprintf("%d",x);
  else
    digit = sprintf("0%d",x);

  return digit;
}

function date(dd,mm,yy)
{
  char *c
  c=sprintf("%s_%s_%s",digits(dd),digits(mm),digits(yy));
  return c
}

BEGIN{ printf("Starting parameter screening process from 1D distributed data at multiple locations\n");
  char *cmd;
}
{
#  if(NR>1){
    cmd=sprintf("./format_ginette.sh %s %s %s %s %s %s %s %s  %s %s %s %s %s  ",$1,$5,date($4,$3,$2),$6,date($4,$3,$2),$8,$7,$9,$10,$11,$12,$13,$14);
#  };



  printf("Launching %s\n",cmd);
  system(cmd);

}
END{printf("Multiple point screening done\n");
}

