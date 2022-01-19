#!/usr/bin/awk -f
BEGIN {OFS="\t";}
{
    split($10,a,"-"); 
    split($11,b,"-");
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,a[1],b[2],$12;
}
