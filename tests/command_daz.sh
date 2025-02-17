#1. run the below with the name of the frame:
python daz_esd_01_prepare_input.py 027A_04887_262625

#2. run the SET correction:
python daz_esd_02_extract_SET.py 027A_04887_262625.csv  earthtides.csv  027A_04887_262625.tide.csv

#3.run the Iono correction:
python daz_esd_03_extract_iono.py 027A_04887_262625.tide.csv 027A_04887_262625.ion.csv

#4.run the PMM correctio:
python daz_esd_04_extract_PMM.py 027A_04887_262625.ion_out_frame.csv 027A_04887_262625.framewithitrf.csv vel_gps_kreemer.nc --add_eu

