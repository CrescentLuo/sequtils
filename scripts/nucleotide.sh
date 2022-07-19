#!/bin/bash
awk -v FS="" '(/^>/){{for(i in a){print i"\t"a[i]"\t"a[i]/cnt}};printf("%s\n",$0);cnt=0;delete a;next;} {for(i=1;i<=NF;i++){a[$i]++;cnt++;};} END {for(i in a){print i"\t"a[i]"\t"a[i]/cnt}}' $1
