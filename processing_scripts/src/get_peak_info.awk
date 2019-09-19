# Short awk program to extract peak table information from ascii text output
# from the Shimadzu GC-2014

BEGIN{
    N2O_area=0
    CO2_area=0
    CH4_area=0
    FS="\t"
}

# If the tenth field contains "Nitrous"
$11~/Nitrous/{
    N2O_area=$5
}

# If the tenth field is "CO2"
$11=="CO2"{
    CO2_area=$5
}

# If the tenth field is "Methane"
$11=="Methane"{
    CH4_area=$5
}

# 
ENDFILE{
printf"%s\t%s\t%s\t%s\n",FILENAME,CH4_area,CO2_area,N2O_area
}
