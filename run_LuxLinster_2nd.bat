@echo off
cls

echo .
echo .


	echo Finding feature pairs
	paython MExtract.py -m AllExtract   -g "E:\Lux\Linster\2nd experiment\29hours_neg.grp"   -l group      -u true   -i false  -o false  -p false      -x -e -s
	paython MExtract.py -m AllExtract   -g "E:\Lux\Linster\2nd experiment\29hours_pos.grp"   -l group      -u true   -i false  -o false  -p false      -x -e -s
	paython MExtract.py -m AllExtract   -g "E:\Lux\Linster\2nd experiment\54hours_neg.grp"   -l group      -u true   -i false  -o false  -p false      -x -e -s
	paython MExtract.py -m AllExtract   -g "E:\Lux\Linster\2nd experiment\54hours_pos.grp"   -l group      -u true   -i false  -o false  -p false      -x -e -s

	echo Grouping feature pairs
	paython MExtract.py -m AllExtract   -g "E:\Lux\Linster\2nd experiment\29hours_neg.grp"   -l group      -u false  -i true   -o true   -p false      -x -e -s
	paython MExtract.py -m AllExtract   -g "E:\Lux\Linster\2nd experiment\29hours_pos.grp"   -l group      -u false  -i true   -o true   -p false      -x -e -s
	paython MExtract.py -m AllExtract   -g "E:\Lux\Linster\2nd experiment\54hours_neg.grp"   -l group      -u false  -i true   -o true   -p false      -x -e -s
	paython MExtract.py -m AllExtract   -g "E:\Lux\Linster\2nd experiment\54hours_pos.grp"   -l group      -u false  -i true   -o true   -p false      -x -e -s

	echo Re-integrating missed feature pairs
	paython MExtract.py -m AllExtract   -g "E:\Lux\Linster\2nd experiment\29hours_neg.grp"   -l group      -u false  -i true   -o false  -p true       -x -e -s
	paython MExtract.py -m AllExtract   -g "E:\Lux\Linster\2nd experiment\29hours_pos.grp"   -l group      -u false  -i true   -o false  -p true       -x -e -s
	paython MExtract.py -m AllExtract   -g "E:\Lux\Linster\2nd experiment\54hours_neg.grp"   -l group      -u false  -i true   -o false  -p true       -x -e -s
	paython MExtract.py -m AllExtract   -g "E:\Lux\Linster\2nd experiment\54hours_pos.grp"   -l group      -u false  -i true   -o false  -p true       -x -e -s

	
	
	
	
	echo Generating statistical reports (currently not working)
	"C:\Program Files\R\R-3.1.0\bin\Rsc#ript.exe" -e "require(knitr); setwd('E:/Lux/Linster/2nd experiment/BSK_LuxLin_exp1_29hours_neg/'); file='BSK_LuxLin_exp1_29hours_neg'; knit(paste(file, '.Rmd', sep='')); pandoc(paste(file, '.md', sep=''), format='html')"
	"C:\Program Files\R\R-3.1.0\bin\Rsc#ript.exe" -e "require(knitr); setwd('E:/Lux/Linster/2nd experiment/BSK_LuxLin_exp1_29hours_pos/'); file='BSK_LuxLin_exp1_29hours_pos'; knit(paste(file, '.Rmd', sep='')); pandoc(paste(file, '.md', sep=''), format='html')"
	"C:\Program Files\R\R-3.1.0\bin\Rsc#ript.exe" -e "require(knitr); setwd('E:/Lux/Linster/2nd experiment/BSK_LuxLin_exp2_54hours_neg/'); file='BSK_LuxLin_exp2_54hours_neg'; knit(paste(file, '.Rmd', sep='')); pandoc(paste(file, '.md', sep=''), format='html')"
	"C:\Program Files\R\R-3.1.0\bin\Rsc#ript.exe" -e "require(knitr); setwd('E:/Lux/Linster/2nd experiment/BSK_LuxLin_exp2_54hours_pos/'); file='BSK_LuxLin_exp2_54hours_pos'; knit(paste(file, '.Rmd', sep='')); pandoc(paste(file, '.md', sep=''), format='html')"
	

echo .
echo .
echo finished