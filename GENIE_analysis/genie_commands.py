''' Script for GENIE generation of events from the fluxes, and subsequent analysis. '''

'''test'''
gevgen -o test1.ghep.root -n 1 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/mutristan_s12_nue.data --cross-sections $GENIE_LUC_XSEC & gevgen -o test2.ghep.root -n 1 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/mutristan_s12_numubar.data --cross-sections $GENIE_LUC_XSEC

cd $GENLUC

''' muTRISTAN s'''

gevgen -o muTs_MD_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/muTs_MD_nue.data --cross-sections $GENIE_LUC_XSEC & gevgen -o muTs_MD_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/muTs_MD_numubar.data --cross-sections $GENIE_LUC_XSEC & gevgen -o muTs_SB_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/muTs_SB_nue.data --cross-sections $GENIE_LUC_XSEC & gevgen -o muTs_SB_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/muTs_SB_numubar.data --cross-sections $GENIE_LUC_XSEC & gevgen -o muTs_SM_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000130270 -f $GENIE_FLUXES/muTs_SM_nue.data --cross-sections $GENIE_LUC_XSEC & gevgen -o muTs_SM_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000130270 -f $GENIE_FLUXES/muTs_SM_numubar.data --cross-sections $GENIE_LUC_XSEC & gevgen -o muTs_HC_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/muTs_SM_nue.data --cross-sections $GENIE_LUC_XSEC & gevgen -o muTs_HC_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/muTs_SM_numubar.data --cross-sections $GENIE_LUC_XSEC & gevgen -o muTs_EC_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/muTs_EC_nue.data --cross-sections $GENIE_LUC_XSEC & gevgen -o muTs_EC_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/muTs_EC_numubar.data --cross-sections $GENIE_LUC_XSEC

gevgen -o muTs_NO_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000741820 -f $GENIE_FLUXES/muTs_NO_nue.data --cross-sections $GENIE_LUC_XSEC & gevgen -o muTs_NO_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000741820 -f $GENIE_FLUXES/muTs_NO_numubar.data --cross-sections $GENIE_LUC_XSEC & gevgen -o mokhov_MD_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/mokhov_MD_nue.data --cross-sections $GENIE_LUC_XSEC & gevgen -o mokhov_MD_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000260560 -f $GENIE_FLUXES/mokhov_MD_nuebar.data --cross-sections $GENIE_LUC_XSEC & gevgen -o mokhov_MD_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000260560 -f $GENIE_FLUXES/mokhov_MD_numu.data --cross-sections $GENIE_LUC_XSEC & gevgen -o mokhov_MD_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/mokhov_MD_numubar.data --cross-sections $GENIE_LUC_XSEC & gevgen -o mokhov_SB_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/mokhov_SB_nue.data --cross-sections $GENIE_LUC_XSEC & gevgen -o mokhov_SB_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000260560 -f $GENIE_FLUXES/mokhov_SB_nuebar.data --cross-sections $GENIE_LUC_XSEC & gevgen -o mokhov_SB_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000260560 -f $GENIE_FLUXES/mokhov_SB_numu.data --cross-sections $GENIE_LUC_XSEC & gevgen -o mokhov_SB_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/mokhov_SB_numubar.data --cross-sections $GENIE_LUC_XSEC & gevgen -o mokhov_SM_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000130270 -f $GENIE_FLUXES/mokhov_SM_nue.data --cross-sections $GENIE_LUC_XSEC & gevgen -o mokhov_SM_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000130270 -f $GENIE_FLUXES/mokhov_SM_nuebar.data --cross-sections $GENIE_LUC_XSEC

gevgen -o lol.ghep.root -n 1 -e 0.1,9999 -p 14 -t 1000260560 -f $GENIE_FLUXES/mokhov_SB_numu.data --cross-sections $GENIE_LUC_XSEC
# done above

gevgen -o mokhov_SM_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000130270 -f $GENIE_FLUXES/mokhov_SM_numu.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_SM_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000130270 -f $GENIE_FLUXES/mokhov_SM_numubar.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_HC_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mokhov_SM_nue.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_HC_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mokhov_SM_nuebar.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_HC_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mokhov_SM_numu.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_HC_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mokhov_SM_numubar.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_EC_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mokhov_EC_nue.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_EC_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mokhov_EC_nuebar.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_EC_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mokhov_EC_numu.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_EC_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mokhov_EC_numubar.data --cross-sections $GENIE_LUC_XSEC &

#done above

gevgen -o mokhov_NO_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000741820 -f $GENIE_FLUXES/mokhov_NO_nue.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_NO_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000741820 -f $GENIE_FLUXES/mokhov_NO_nuebar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_NO_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000741820 -f $GENIE_FLUXES/mokhov_NO_numu.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_NO_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000741820 -f $GENIE_FLUXES/mokhov_NO_numubar.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mucols2_MD_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/mucols2_MD_nue.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mucols2_MD_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000260560 -f $GENIE_FLUXES/mucols2_MD_nuebar.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mucols2_MD_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000260560 -f $GENIE_FLUXES/mucols2_MD_numu.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mucols2_MD_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/mucols2_MD_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SB_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/mucols2_SB_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SB_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000260560 -f $GENIE_FLUXES/mucols2_SB_nuebar.data --cross-sections $GENIE_LUC_XSEC


gevgen -o mucols2_SB_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000260560 -f $GENIE_FLUXES/mucols2_SB_numu.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SB_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/mucols2_SB_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SM_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000130270 -f $GENIE_FLUXES/mucols2_SM_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SM_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000130270 -f $GENIE_FLUXES/mucols2_SM_nuebar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SM_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000130270 -f $GENIE_FLUXES/mucols2_SM_numu.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SM_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000130270 -f $GENIE_FLUXES/mucols2_SM_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_HC_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mucols2_SM_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_HC_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mucols2_SM_nuebar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_HC_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mucols2_SM_numu.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_HC_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mucols2_SM_numubar.data --cross-sections $GENIE_LUC_XSEC


gevgen -o mucols2_EC_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mucols2_EC_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_EC_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mucols2_EC_nuebar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_EC_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mucols2_EC_numu.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_EC_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mucols2_EC_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_NO_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000741820 -f $GENIE_FLUXES/mucols2_NO_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_NO_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000741820 -f $GENIE_FLUXES/mucols2_NO_nuebar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_NO_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000741820 -f $GENIE_FLUXES/mucols2_NO_numu.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_NO_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000741820 -f $GENIE_FLUXES/mucols2_NO_numubar.data --cross-sections $GENIE_LUC_XSEC 

#wrote above

'''Subsequent analysis'''
cd $GENANA

GetPL4P -e muTs -p 12,-14 -o muTs.txt
GetPL4P -e mokhov -p 12,-12 -o mokhovnue.txt
GetPL4P -e mokhov -p 14,-14 -o mokhovnumu.txt
GetPL4P -e mucols2 -p 12,-12 -o mucols2nue.txt
GetPL4P -e mucols2 -p 14,-14 -o mucols2numu.txt



''' to check if properly done (all) '''
gevdump -f muTs_MD_nue.ghep.root -n 0
gevdump -f muTs_MD_numubar.ghep.root -n 0

gevdump -f muTs_SB_nue.ghep.root -n 0
gevdump -f muTs_SB_numubar.ghep.root -n 0

gevdump -f muTs_SM_nue.ghep.root -n 0
gevdump -f muTs_SM_numubar.ghep.root -n 0

gevdump -f muTs_HC_nue.ghep.root -n 0
gevdump -f muTs_HC_numubar.ghep.root -n 0

gevdump -f muTs_EC_nue.ghep.root -n 0
gevdump -f muTs_EC_numubar.ghep.root -n 0

gevdump -f muTs_NO_nue.ghep.root -n 0
gevdump -f muTs_NO_numubar.ghep.root -n 0


gevdump -f mokhov_MD_nue.ghep.root -n 0
gevdump -f mokhov_MD_nuebar.ghep.root -n 0
gevdump -f mokhov_MD_numu.ghep.root -n 0
gevdump -f mokhov_MD_numubar.ghep.root -n 0

gevdump -f mokhov_SB_nue.ghep.root -n 0
gevdump -f mokhov_SB_nuebar.ghep.root -n 0
gevdump -f mokhov_SB_numu.ghep.root -n 0
gevdump -f mokhov_SB_numubar.ghep.root -n 0

gevdump -f mokhov_SM_nue.ghep.root -n 0
gevdump -f mokhov_SM_nuebar.ghep.root -n 0
gevdump -f mokhov_SM_numu.ghep.root -n 0
gevdump -f mokhov_SM_numubar.ghep.root -n 0

gevdump -f mokhov_HC_nue.ghep.root -n 0
gevdump -f mokhov_HC_nuebar.ghep.root -n 0
gevdump -f mokhov_HC_numu.ghep.root -n 0
gevdump -f mokhov_HC_numubar.ghep.root -n 0

gevdump -f mokhov_EC_nue.ghep.root -n 0
gevdump -f mokhov_EC_nuebar.ghep.root -n 0
gevdump -f mokhov_EC_numu.ghep.root -n 0
gevdump -f mokhov_EC_numubar.ghep.root -n 0

gevdump -f mokhov_NO_nue.ghep.root -n 0
gevdump -f mokhov_NO_nuebar.ghep.root -n 0
gevdump -f mokhov_NO_numu.ghep.root -n 0
gevdump -f mokhov_NO_numubar.ghep.root -n 0


gevdump -f mucols2_MD_nue.ghep.root -n 0
gevdump -f mucols2_MD_nuebar.ghep.root -n 0
gevdump -f mucols2_MD_numu.ghep.root -n 0
gevdump -f mucols2_MD_numubar.ghep.root -n 0

gevdump -f mucols2_SB_nue.ghep.root -n 0
gevdump -f mucols2_SB_nuebar.ghep.root -n 0
gevdump -f mucols2_SB_numu.ghep.root -n 0
gevdump -f mucols2_SB_numubar.ghep.root -n 0

gevdump -f mucols2_SM_nue.ghep.root -n 0
gevdump -f mucols2_SM_nuebar.ghep.root -n 0
gevdump -f mucols2_SM_numu.ghep.root -n 0
gevdump -f mucols2_SM_numubar.ghep.root -n 0

gevdump -f mucols2_HC_nue.ghep.root -n 0
gevdump -f mucols2_HC_nuebar.ghep.root -n 0
gevdump -f mucols2_HC_numu.ghep.root -n 0
gevdump -f mucols2_HC_numubar.ghep.root -n 0

gevdump -f mucols2_EC_nue.ghep.root -n 0
gevdump -f mucols2_EC_nuebar.ghep.root -n 0
gevdump -f mucols2_EC_numu.ghep.root -n 0
gevdump -f mucols2_EC_numubar.ghep.root -n 0

gevdump -f mucols2_NO_nue.ghep.root -n 0
gevdump -f mucols2_NO_nuebar.ghep.root -n 0
gevdump -f mucols2_NO_numu.ghep.root -n 0
gevdump -f mucols2_NO_numubar.ghep.root -n 0





(
gevgen -o muTs_MD_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/muTs_MD_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o muTs_MD_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/muTs_MD_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o muTs_SB_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/muTs_SB_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o muTs_SB_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/muTs_SB_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o muTs_SM_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000130270 -f $GENIE_FLUXES/muTs_SM_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o muTs_SM_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000130270 -f $GENIE_FLUXES/muTs_SM_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o muTs_HC_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/muTs_SM_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o muTs_HC_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/muTs_SM_numubar.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o muTs_EC_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/muTs_EC_nue.data --cross-sections $GENIE_LUC_XSEC & gevgen -o muTs_EC_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/muTs_EC_numubar.data --cross-sections $GENIE_LUC_XSEC
wait
gevgen -o muTs_NO_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000741820 -f $GENIE_FLUXES/muTs_NO_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o muTs_NO_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000741820 -f $GENIE_FLUXES/muTs_NO_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_MD_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/mokhov_MD_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_MD_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000260560 -f $GENIE_FLUXES/mokhov_MD_nuebar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_MD_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000260560 -f $GENIE_FLUXES/mokhov_MD_numu.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_MD_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/mokhov_MD_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_SB_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/mokhov_SB_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_SB_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000260560 -f $GENIE_FLUXES/mokhov_SB_nuebar.data --cross-sections $GENIE_LUC_XSEC 
wait)

(
gevgen -o mokhov_SB_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000260560 -f $GENIE_FLUXES/mokhov_SB_numu.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_SB_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/mokhov_SB_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_SM_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000130270 -f $GENIE_FLUXES/mokhov_SM_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_SM_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000130270 -f $GENIE_FLUXES/mokhov_SM_nuebar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_SM_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000130270 -f $GENIE_FLUXES/mokhov_SM_numu.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_SM_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000130270 -f $GENIE_FLUXES/mokhov_SM_numubar.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_HC_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mokhov_SM_nue.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_HC_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mokhov_SM_nuebar.data --cross-sections $GENIE_LUC_XSEC
wait
gevgen -o mokhov_HC_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mokhov_SM_numu.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_HC_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mokhov_SM_numubar.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_EC_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mokhov_EC_nue.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_EC_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mokhov_EC_nuebar.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_EC_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mokhov_EC_numu.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_EC_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mokhov_EC_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_NO_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000741820 -f $GENIE_FLUXES/mokhov_NO_nue.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_NO_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000741820 -f $GENIE_FLUXES/mokhov_NO_nuebar.data --cross-sections $GENIE_LUC_XSEC
wait
)
(
gevgen -o mokhov_NO_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000741820 -f $GENIE_FLUXES/mokhov_NO_numu.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_NO_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000741820 -f $GENIE_FLUXES/mokhov_NO_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_MD_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/mucols2_MD_nue.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mucols2_MD_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000260560 -f $GENIE_FLUXES/mucols2_MD_nuebar.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mucols2_MD_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000260560 -f $GENIE_FLUXES/mucols2_MD_numu.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mucols2_MD_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/mucols2_MD_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SB_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/mucols2_SB_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SB_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000260560 -f $GENIE_FLUXES/mucols2_SB_nuebar.data --cross-sections $GENIE_LUC_XSEC
wait
gevgen -o mucols2_SB_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000260560 -f $GENIE_FLUXES/mucols2_SB_numu.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SB_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/mucols2_SB_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SM_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000130270 -f $GENIE_FLUXES/mucols2_SM_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SM_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000130270 -f $GENIE_FLUXES/mucols2_SM_nuebar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SM_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000130270 -f $GENIE_FLUXES/mucols2_SM_numu.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SM_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000130270 -f $GENIE_FLUXES/mucols2_SM_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_HC_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mucols2_SM_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_HC_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mucols2_SM_nuebar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_HC_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mucols2_SM_numu.data --cross-sections $GENIE_LUC_XSEC
wait
gevgen -o mucols2_HC_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mucols2_SM_numubar.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mucols2_EC_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mucols2_EC_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_EC_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mucols2_EC_nuebar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_EC_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mucols2_EC_numu.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_EC_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mucols2_EC_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_NO_nue.ghep.root -n 100000 -e 0.1,9999 -p 12 -t 1000741820 -f $GENIE_FLUXES/mucols2_NO_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_NO_nuebar.ghep.root -n 100000 -e 0.1,9999 -p -12 -t 1000741820 -f $GENIE_FLUXES/mucols2_NO_nuebar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_NO_numu.ghep.root -n 100000 -e 0.1,9999 -p 14 -t 1000741820 -f $GENIE_FLUXES/mucols2_NO_numu.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_NO_numubar.ghep.root -n 100000 -e 0.1,9999 -p -14 -t 1000741820 -f $GENIE_FLUXES/mucols2_NO_numubar.data --cross-sections $GENIE_LUC_XSEC 
wait    
gevgen -o lol2.ghep.root -n 1 -e 0.1,9999 -p 14 -t 1000260560 -f $GENIE_FLUXES/mokhov_SB_numu.data --cross-sections $GENIE_LUC_XSEC
)



(
gevgen -o muTs_MD_nue.ghep.root -n 1000 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/muTs_MD_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o muTs_MD_numubar.ghep.root -n 1000 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/muTs_MD_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o muTs_SB_nue.ghep.root -n 1000 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/muTs_SB_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o muTs_SB_numubar.ghep.root -n 1000 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/muTs_SB_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o muTs_SM_nue.ghep.root -n 1000 -e 0.1,9999 -p 12 -t 1000130270 -f $GENIE_FLUXES/muTs_SM_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o muTs_SM_numubar.ghep.root -n 1000 -e 0.1,9999 -p -14 -t 1000130270 -f $GENIE_FLUXES/muTs_SM_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o muTs_HC_nue.ghep.root -n 1000 -e 0.1,9999 -p 12 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/muTs_SM_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o muTs_HC_numubar.ghep.root -n 1000 -e 0.1,9999 -p -14 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/muTs_SM_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o muTs_EC_nue.ghep.root -n 1000 -e 0.1,9999 -p 12 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/muTs_EC_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o muTs_EC_numubar.ghep.root -n 1000 -e 0.1,9999 -p -14 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/muTs_EC_numubar.data --cross-sections $GENIE_LUC_XSEC
wait
gevgen -o muTs_NO_nue.ghep.root -n 1000 -e 0.1,9999 -p 12 -t 1000741820 -f $GENIE_FLUXES/muTs_NO_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o muTs_NO_numubar.ghep.root -n 1000 -e 0.1,9999 -p -14 -t 1000741820 -f $GENIE_FLUXES/muTs_NO_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_MD_nue.ghep.root -n 1000 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/mokhov_MD_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_MD_nuebar.ghep.root -n 1000 -e 0.1,9999 -p -12 -t 1000260560 -f $GENIE_FLUXES/mokhov_MD_nuebar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_MD_numu.ghep.root -n 1000 -e 0.1,9999 -p 14 -t 1000260560 -f $GENIE_FLUXES/mokhov_MD_numu.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_MD_numubar.ghep.root -n 1000 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/mokhov_MD_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_SB_nue.ghep.root -n 1000 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/mokhov_SB_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_SB_nuebar.ghep.root -n 1000 -e 0.1,9999 -p -12 -t 1000260560 -f $GENIE_FLUXES/mokhov_SB_nuebar.data --cross-sections $GENIE_LUC_XSEC 
wait
gevgen -o mokhov_SB_numu.ghep.root -n 1000 -e 0.1,9999 -p 14 -t 1000260560 -f $GENIE_FLUXES/mokhov_SB_numu.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_SB_numubar.ghep.root -n 1000 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/mokhov_SB_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_SM_nue.ghep.root -n 1000 -e 0.1,9999 -p 12 -t 1000130270 -f $GENIE_FLUXES/mokhov_SM_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_SM_nuebar.ghep.root -n 1000 -e 0.1,9999 -p -12 -t 1000130270 -f $GENIE_FLUXES/mokhov_SM_nuebar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_SM_numu.ghep.root -n 1000 -e 0.1,9999 -p 14 -t 1000130270 -f $GENIE_FLUXES/mokhov_SM_numu.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_SM_numubar.ghep.root -n 1000 -e 0.1,9999 -p -14 -t 1000130270 -f $GENIE_FLUXES/mokhov_SM_numubar.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_HC_nue.ghep.root -n 1000 -e 0.1,9999 -p 12 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mokhov_SM_nue.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_HC_nuebar.ghep.root -n 1000 -e 0.1,9999 -p -12 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mokhov_SM_nuebar.data --cross-sections $GENIE_LUC_XSEC
wait
gevgen -o mokhov_HC_numu.ghep.root -n 1000 -e 0.1,9999 -p 14 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mokhov_SM_numu.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_HC_numubar.ghep.root -n 1000 -e 0.1,9999 -p -14 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mokhov_SM_numubar.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_EC_nue.ghep.root -n 1000 -e 0.1,9999 -p 12 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mokhov_EC_nue.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_EC_nuebar.ghep.root -n 1000 -e 0.1,9999 -p -12 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mokhov_EC_nuebar.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_EC_numu.ghep.root -n 1000 -e 0.1,9999 -p 14 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mokhov_EC_numu.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_EC_numubar.ghep.root -n 1000 -e 0.1,9999 -p -14 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mokhov_EC_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mokhov_NO_nue.ghep.root -n 1000 -e 0.1,9999 -p 12 -t 1000741820 -f $GENIE_FLUXES/mokhov_NO_nue.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_NO_nuebar.ghep.root -n 1000 -e 0.1,9999 -p -12 -t 1000741820 -f $GENIE_FLUXES/mokhov_NO_nuebar.data --cross-sections $GENIE_LUC_XSEC
wait
gevgen -o mokhov_NO_numu.ghep.root -n 1000 -e 0.1,9999 -p 14 -t 1000741820 -f $GENIE_FLUXES/mokhov_NO_numu.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mokhov_NO_numubar.ghep.root -n 1000 -e 0.1,9999 -p -14 -t 1000741820 -f $GENIE_FLUXES/mokhov_NO_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_MD_nue.ghep.root -n 1000 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/mucols2_MD_nue.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mucols2_MD_nuebar.ghep.root -n 1000 -e 0.1,9999 -p -12 -t 1000260560 -f $GENIE_FLUXES/mucols2_MD_nuebar.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mucols2_MD_numu.ghep.root -n 1000 -e 0.1,9999 -p 14 -t 1000260560 -f $GENIE_FLUXES/mucols2_MD_numu.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mucols2_MD_numubar.ghep.root -n 1000 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/mucols2_MD_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SB_nue.ghep.root -n 1000 -e 0.1,9999 -p 12 -t 1000260560 -f $GENIE_FLUXES/mucols2_SB_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SB_nuebar.ghep.root -n 1000 -e 0.1,9999 -p -12 -t 1000260560 -f $GENIE_FLUXES/mucols2_SB_nuebar.data --cross-sections $GENIE_LUC_XSEC
wait
gevgen -o mucols2_SB_numu.ghep.root -n 1000 -e 0.1,9999 -p 14 -t 1000260560 -f $GENIE_FLUXES/mucols2_SB_numu.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SB_numubar.ghep.root -n 1000 -e 0.1,9999 -p -14 -t 1000260560 -f $GENIE_FLUXES/mucols2_SB_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SM_nue.ghep.root -n 1000 -e 0.1,9999 -p 12 -t 1000130270 -f $GENIE_FLUXES/mucols2_SM_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SM_nuebar.ghep.root -n 1000 -e 0.1,9999 -p -12 -t 1000130270 -f $GENIE_FLUXES/mucols2_SM_nuebar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SM_numu.ghep.root -n 1000 -e 0.1,9999 -p 14 -t 1000130270 -f $GENIE_FLUXES/mucols2_SM_numu.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_SM_numubar.ghep.root -n 1000 -e 0.1,9999 -p -14 -t 1000130270 -f $GENIE_FLUXES/mucols2_SM_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_HC_nue.ghep.root -n 1000 -e 0.1,9999 -p 12 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mucols2_SM_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_HC_nuebar.ghep.root -n 1000 -e 0.1,9999 -p -12 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mucols2_SM_nuebar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_HC_numu.ghep.root -n 1000 -e 0.1,9999 -p 14 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mucols2_SM_numu.data --cross-sections $GENIE_LUC_XSEC
wait
gevgen -o mucols2_HC_numubar.ghep.root -n 1000 -e 0.1,9999 -p -14 -t 1000260560[0.84],1000130270[0.03],1000290630[0.004],1000060120[0.126] -f $GENIE_FLUXES/mucols2_SM_numubar.data --cross-sections $GENIE_LUC_XSEC & 
gevgen -o mucols2_EC_nue.ghep.root -n 1000 -e 0.1,9999 -p 12 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mucols2_EC_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_EC_nuebar.ghep.root -n 1000 -e 0.1,9999 -p -12 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mucols2_EC_nuebar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_EC_numu.ghep.root -n 1000 -e 0.1,9999 -p 14 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mucols2_EC_numu.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_EC_numubar.ghep.root -n 1000 -e 0.1,9999 -p -14 -t 1000741820[0.404],1000290630[0.489],1000140280[0.107] -f $GENIE_FLUXES/mucols2_EC_numubar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_NO_nue.ghep.root -n 1000 -e 0.1,9999 -p 12 -t 1000741820 -f $GENIE_FLUXES/mucols2_NO_nue.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_NO_nuebar.ghep.root -n 1000 -e 0.1,9999 -p -12 -t 1000741820 -f $GENIE_FLUXES/mucols2_NO_nuebar.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_NO_numu.ghep.root -n 1000 -e 0.1,9999 -p 14 -t 1000741820 -f $GENIE_FLUXES/mucols2_NO_numu.data --cross-sections $GENIE_LUC_XSEC &
gevgen -o mucols2_NO_numubar.ghep.root -n 1000 -e 0.1,9999 -p -14 -t 1000741820 -f $GENIE_FLUXES/mucols2_NO_numubar.data --cross-sections $GENIE_LUC_XSEC 
wait    
gevgen -o lol2.ghep.root -n 1 -e 0.1,9999 -p 14 -t 1000260560 -f $GENIE_FLUXES/mokhov_SB_numu.data --cross-sections $GENIE_LUC_XSEC
)