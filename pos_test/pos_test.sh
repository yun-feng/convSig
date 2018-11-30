#!/bin/bash

pos=(398 448 548 598 748)
probability=(0.3 0.5 0.7)

for prob in ${probability[*]}
do
	for pos_c in ${pos[*]}
	do
		cd /data/ted/multi
		perl ./prepare_open_rc_spikein_pos.pl $pos_c $prob > out_esopre_pos
		perl ./convert_pos.pl $pos_c $prob > out_esocon_pos
		cd ./eso/deep
		/apps/well/torch/20170221-p100-gcc5.4.0/bin/th main_open_bg_pos.lua $pos_c $prob > out_bg_pos
		cd ../..
		/apps/well/torch/20170221-p100-gcc5.4.0/bin/th feature_pos.lua $pos_c $prob
	done
done
