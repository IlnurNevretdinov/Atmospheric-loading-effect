gmt begin ../plots/expansion_errors.pdf
    file="../data/sims/expansion_errors.nc"
    I="I+ne0.5+a200"
    gmt makecpt -Cjet -T-20/20/1

    gmt grdimage $file -$I -R0/190/40/89 -JM22c -Bx30 -By10 -BWsNE+t$"Errors of spherical expansion" 
    
  #  gmt grdcontour $file -C20 -A40+f7p -B
    gmt colorbar -DJBC+o0/0.4c+w6i/0.25c -I --FONT_ANNOT_PRIMARY=12p -Bx5f2 -By+lhPa	


gmt end
