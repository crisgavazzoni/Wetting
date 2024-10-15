#!/bin/bash  

gcc -Wall -O3 potts.c -lm -o gota3d.out 

for L in 240
	do

for R in 50
	do

for Omg in 12.0
	do

for w in 5
	do

for a in 5
	do

for h in 10
	do

for fo in 0.00 
	do

for CI in 2
	do


# ==============================================================================
# =                               Rodando simulação                            = 
#===============================================================================

./gota3d.out -L ${L} -R ${R} -Omg ${Omg} -fo ${fo} -CI ${CI} -w ${w} -h ${h} -a ${a} -s  12345567 


# ==============================================================================
# =                            Movendo os aqrquivos CB                         = 
#===============================================================================
#if [ ${CI} -eq 1 ] 
#then

#rm ../dsf/CB_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}.dsf
#rm ../xyz/CB_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}*
#rm ../conf/CB_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}*
#rm ../out/CB_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}.out
#rm ../ener/CB_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}.ener

#mv CB_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}.dsf ../dsf/
#mv CB_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}.xyz ../xyz/
#mv CB_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}_init.xyz ../xyz/
#mv CB_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}_conf.dsf ../conf/
#mv CB_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}_LAST.dsf ../conf/
#mv CB_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}.out ../out/
#mv CB_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}.ener ../ener/

#fi

# ==============================================================================
# =                           Movendo os arquivos WE                           = 
#===============================================================================
#if [ ${CI} -eq 2 ] 
#then

#rm ../dsf/WE_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}.dsf
#rm ../xyz/WE_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}*
#rm ../conf/WE_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}*
#rm ../out/WE_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}.out
#rm ../ener/WE_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}.ener

#mv WE_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}.dsf ../dsf/
#mv WE_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}.xyz ../xyz/
#mv WE_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}_init.xyz ../xyz/
#mv WE_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}_conf.dsf ../conf/
#mv WE_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}_LAST.dsf ../conf/
#mv WE_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}.out ../out/
#mv WE_gota_3d_L_${L}_R_${R}_Thc_${Thc}_Omg_${Omg}_w_${w}_h_${h}_a_${a}.ener ../ener/

#fi 
#===============================================================================

done
done
done
done
done
done
done
done

