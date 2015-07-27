#! /bin/bash
STR="plompt_levelt_pierce_0.1.1"
./arst_neio -o $STR
echo arst_neio executed
#writing .wav
csound -d -o "arst_neio_$STR".wav "arst_neio_$STR".orc "arst_neio_$STR".sco
#csound --format=24bit -d -o compo_pl.wav compo_pl.orc compo_pl.sco
echo $STR.wav generated
