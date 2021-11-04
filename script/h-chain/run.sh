make rmdat
python ../script/h-chain/generateInputFileHChain.py
./ghf 
python ../script/dataAnalyze.py 1 
python ../script/h-chain/checkSpin_threeband.py Sz_1_Average.dat
python ../script/h-chain/checkSpin_threeband.py Sx_1_Average.dat
python ../script/h-chain/checkSpin_threeband.py Sy_1_Average.dat
python ../script/h-chain/checkSpin_threeband.py Nup_1_Average.dat
python ../script/h-chain/checkSpin_threeband.py Ndn_1_Average.dat
