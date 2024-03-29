#!/bin/bash

sizeKDT=100000

#./smallpm -film-name img_mem/ej_foggy_dielc_50 -scene 5 -pm-nb-nearest-neighbor 50 -pm-max-photons-shot 1000000 -pm-photons-volume 15000 -pm-photons-caustic 10000 -pm-photons-global 10000 -volumetric 

#./smallpm -film-name img_mem/ej_foggy_50 -scene 2 -pm-nb-nearest-neighbor 50 -pm-max-photons-shot 1000000 -pm-photons-volume 15000 -pm-photons-caustic 10000 -pm-photons-global 10000 -volumetric

#./smallpm -film-name img_mem/ej_foggy_200 -scene 2 -pm-nb-nearest-neighbor 200 -pm-max-photons-shot 1000000 -pm-photons-volume 15000 -pm-photons-caustic 10000 -pm-photons-global 10000 -volumetric 

#./smallpm -film-name img_mem/ej_foggy_350 -scene 2 -pm-nb-nearest-neighbor 350 -pm-max-photons-shot 1000000 -pm-photons-volume 15000 -pm-photons-caustic 10000 -pm-photons-global 10000 -volumetric 


#./smallpm -film-name img_mem/Q23_n_1000 -scene 1 -pm-nb-nearest-neighbor 10 -pm-max-photons-shot 1000 -pm-photons-caustic $sizeKDT -pm-photons-global $sizeKDT -pm-raytraced-direct
#./smallpm -film-name img_mem/Q23_n_10000 -scene 1 -pm-nb-nearest-neighbor 10 -pm-max-photons-shot 10000 -pm-photons-caustic $sizeKDT -pm-photons-global $sizeKDT -pm-raytraced-direct
#./smallpm -film-name img_mem/Q23_n_100000 -scene 1 -pm-nb-nearest-neighbor 10 -pm-max-photons-shot 100000 -pm-photons-caustic $sizeKDT -pm-photons-global $sizeKDT -pm-raytraced-direct

#./smallpm -film-name img_mem/Q23_n_100000_k_1 -scene 1 -pm-nb-nearest-neighbor 1 -pm-max-photons-shot 100000 -pm-photons-caustic $sizeKDT -pm-photons-global $sizeKDT -pm-raytraced-direct
#./smallpm -film-name img_mem/Q23_n_100000_k_10 -scene 1 -pm-nb-nearest-neighbor 10 -pm-max-photons-shot 100000 -pm-photons-caustic $sizeKDT -pm-photons-global $sizeKDT -pm-raytraced-direct
#./smallpm -film-name img_mem/Q23_n_100000_k_50 -scene 1 -pm-nb-nearest-neighbor 50 -pm-max-photons-shot 100000 -pm-photons-caustic $sizeKDT -pm-photons-global $sizeKDT -pm-raytraced-direct
#./smallpm -film-name img_mem/Q23_n_100000_k_100 -scene 1 -pm-nb-nearest-neighbor 100 -pm-max-photons-shot 100000 -pm-photons-caustic $sizeKDT -pm-photons-global $sizeKDT -pm-raytraced-direct

#./smallpm -film-name img_mem/Q22_n_100000_k_1000 -scene 1 -pm-nb-nearest-neighbor 1000 -pm-max-photons-shot 100000 -pm-photons-caustic $sizeKDT -pm-photons-global $sizeKDT -pm-raytraced-direct

#./smallpm -film-name img_mem/Q24_toroid_rt_ID -scene 4 -pm-nb-nearest-neighbor 1000 -pm-max-photons-shot 100000 -pm-photons-caustic $sizeKDT -pm-photons-global $sizeKDT -pm-raytraced-direct
#./smallpm -film-name img_mem/Q24_toroid_pt_ID -scene 4 -pm-nb-nearest-neighbor 1000 -pm-max-photons-shot 100000 -pm-photons-caustic $sizeKDT -pm-photons-global $sizeKDT 

#./smallpm -film-name img_mem/Q21_rt_ID -scene 2 -pm-nb-nearest-neighbor 1000 -pm-max-photons-shot 100000 -pm-photons-caustic $sizeKDT -pm-photons-global $sizeKDT -pm-raytraced-direct
#./smallpm -film-name img_mem/Q21_pm_ID -scene 2 -pm-nb-nearest-neighbor 1000 -pm-max-photons-shot 100000 -pm-photons-caustic $sizeKDT -pm-photons-global $sizeKDT 

#./smallpm -film-name img_mem/RR_SAMP2 -scene 6 -pm-nb-nearest-neighbor 1000 -pm-max-photons-shot 100000 -pm-photons-caustic $sizeKDT -pm-photons-global $sizeKDT -pm-raytraced-direct

#./smallpm -film-name img_mem/normal -scene 0 -pm-nb-nearest-neighbor 50 -pm-max-photons-shot 100000 -pm-photons-caustic $sizeKDT -pm-photons-global $sizeKDT -pm-raytraced-direct
#./smallpm -film-name img_mem/gauss -scene 0 -pm-nb-nearest-neighbor 50 -pm-max-photons-shot 100000 -pm-photons-caustic $sizeKDT -pm-photons-global $sizeKDT -pm-raytraced-direct
#./smallpm -film-name img_mem/cone -scene 0 -pm-nb-nearest-neighbor 50 -pm-max-photons-shot 100000 -pm-photons-caustic $sizeKDT -pm-photons-global $sizeKDT -pm-raytraced-direct
