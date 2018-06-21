awk 'FNR == 1 { nfiles++; ncols = NF }
     { for (i = 1; i < NF; i++) sum[FNR,i] += $i
       if (FNR > maxnr) maxnr = FNR
     }
     END {
         for (line = 1; line <= maxnr; line++)
         {
             for (col = 1; col < ncols; col++)
                  printf "  %f", sum[line,col]/nfiles;
             printf "\n"
         }
     }' pr350k*.txt
