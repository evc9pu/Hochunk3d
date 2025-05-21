#! /bin/csh  


/bin/rm -rf idl_3d
/bin/rm -rf python
/bin/rm -rf raytracing
/bin/rm -rf srcdust_findwall
/bin/rm -rf srcdust_old

cd models
/bin/rm -rf eclump1
/bin/rm -rf eclump2
/bin/rm -rf modcII_fspot*
/bin/rm -rf modcII_gap
/bin/rm -rf modcII_warp
/bin/rm -rf modcI_fract
/bin/rm -rf modcI_fract_outflow
/bin/rm -rf rcrb
#mv modcII tmp
#/bin/rm -rf modcII*
#mv tmp modcII
#/bin/rm -rf modcI_*
#/bin/rm -rf modcI-*
#/bin/rm -rf modc0*

exit

