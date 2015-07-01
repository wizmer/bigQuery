if (( $# < 1 )); then
    echo Missing input file argument
    exit
fi

cat $1 | awk -F',' '{if(NR > 1){ printf("%s,",$1); for(i = 2;i<NF;i++){if(NR == i){printf("1,")}else{printf("0,")}} if(NF == NR){ printf("1\n")} else{ printf("0\n")}} else{ print $0}}'
