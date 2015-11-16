@echo off
cls

echo .
echo .


	echo finding feature pairs
	#python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\_1K.grp"                               -l "E:\Cambridge\JanApr_2015\_1K.grp"                                   -x -s -e
	#python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\_10K.grp"                              -l "E:\Cambridge\JanApr_2015\_10K.grp"                                  -x -s -e
	#python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\_20K.grp"                              -l "E:\Cambridge\JanApr_2015\_20K.grp"                                  -x -s -e
	#python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\_1K_fastSwitch.grp"                    -l "E:\Cambridge\JanApr_2015\_1K_fastSwitch.grp"                        -x -s -e
	#python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\_10K_fastSwitch.grp"                   -l "E:\Cambridge\JanApr_2015\_10K_fastSwitch.grp"                       -x -s -e
	#python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\_20K_fastSwitch.grp"                   -l "E:\Cambridge\JanApr_2015\_20K_fastSwitch.grp"                       -x -s -e
	
	echo bracketing
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\MouseES.grp"                           -l "E:\Cambridge\JanApr_2015\MouseES.grp"                               -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\MouseES-fastswitch.grp"                -l "E:\Cambridge\JanApr_2015\MouseES-fastswitch.grp"                    -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\QE_PR493_PV_051014-AIF-02.grp"         -l "E:\Cambridge\JanApr_2015\QE_PR493_PV_051014-AIF-02.grp"             -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\QE_PR493_PV_051014-AIF-03.grp"         -l "E:\Cambridge\JanApr_2015\QE_PR493_PV_051014-AIF-03.grp"             -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\QE_PR493_PV_051014-AIF-04.grp"         -l "E:\Cambridge\JanApr_2015\QE_PR493_PV_051014-AIF-04.grp"             -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\QE_PR493_PV_051014-AIF-05.grp"         -l "E:\Cambridge\JanApr_2015\QE_PR493_PV_051014-AIF-05.grp"             -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\QE_PR493_PV_051014-AIF-06.grp"         -l "E:\Cambridge\JanApr_2015\QE_PR493_PV_051014-AIF-06.grp"             -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\QE_PR493_PV_051014-AIF-11.grp"         -l "E:\Cambridge\JanApr_2015\QE_PR493_PV_051014-AIF-11.grp"             -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\QE_PR493_PV_051014-AIF-12.grp"         -l "E:\Cambridge\JanApr_2015\QE_PR493_PV_051014-AIF-12.grp"             -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\QE_PR533_PV_201214-01.grp"             -l "E:\Cambridge\JanApr_2015\QE_PR533_PV_201214-01.grp"                 -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\QE_PR533_PV_201214-03.grp"             -l "E:\Cambridge\JanApr_2015\QE_PR533_PV_201214-03.grp"                 -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\QE_PR533_PV_201214-05.grp"             -l "E:\Cambridge\JanApr_2015\QE_PR533_PV_201214-05.grp"                 -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\QE_PR533_PV_201214-06.grp"             -l "E:\Cambridge\JanApr_2015\QE_PR533_PV_201214-06.grp"                 -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\QE_PR536_PV_180115_AIF-02.grp"         -l "E:\Cambridge\JanApr_2015\QE_PR536_PV_180115_AIF-02.grp"             -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\QE_PR536_PV_180115_AIF-06.grp"         -l "E:\Cambridge\JanApr_2015\QE_PR536_PV_180115_AIF-06.grp"             -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\QE_PR536_PV_180115_AIF-09.grp"         -l "E:\Cambridge\JanApr_2015\QE_PR536_PV_180115_AIF-09.grp"             -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\RNA-elegans-large_150328005031.grp"    -l "E:\Cambridge\JanApr_2015\RNA-elegans-large_150328005031.grp"        -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\RNA-elegans-large13C.grp"              -l "E:\Cambridge\JanApr_2015\RNA-elegans-large13C.grp"                  -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\RNA-elegans-small.grp"                 -l "E:\Cambridge\JanApr_2015\RNA-elegans-small.grp"                     -x -s -e
	python MExtract.py -m TracExtract    -g "E:\Cambridge\JanApr_2015\RNA-elegans-small_150328003042.grp"    -l "E:\Cambridge\JanApr_2015\RNA-elegans-small_150328003042.grp"        -x -s -e

echo .
echo .
echo finished
