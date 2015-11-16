@echo off
cls

echo .
echo .


	echo finding feature pairs
	python MExtract.py -m AllExtract -g "E:\exp134_part123\groups.grp"  -l "E:\exp134_part123\groups.grp"   -u true -i false -o false -p false -e -x

	echo grouping feature pairs
	python MExtract.py -m AllExtract -g "E:\exp134_part123\groups.grp"  -l "E:\exp134_part123\groups.grp"   -u false -i true -o true -p false -e -x -s


echo .
echo .
echo re-integrating files
echo .
echo .


	python MExtract.py -m AllExtract -g "E:\exp134_part123\1-14.grp"    -l "E:\exp134_part123\groups.grp"   -u false -i true -o false -p true -e -x -s
	python MExtract.py -m AllExtract -g "E:\exp134_part123\15-28.grp"   -l "E:\exp134_part123\groups.grp"   -u false -i true -o false -p true -e -x -s
	python MExtract.py -m AllExtract -g "E:\exp134_part123\29-42.grp"   -l "E:\exp134_part123\groups.grp"   -u false -i true -o false -p true -e -x -s
	python MExtract.py -m AllExtract -g "E:\exp134_part123\43-56.grp"   -l "E:\exp134_part123\groups.grp"   -u false -i true -o false -p true -e -x -s
	python MExtract.py -m AllExtract -g "E:\exp134_part123\57-70.grp"   -l "E:\exp134_part123\groups.grp"   -u false -i true -o false -p true -e -x -s
	python MExtract.py -m AllExtract -g "E:\exp134_part123\71-84.grp"   -l "E:\exp134_part123\groups.grp"   -u false -i true -o false -p true -e -x -s
	python MExtract.py -m AllExtract -g "E:\exp134_part123\85-98.grp"   -l "E:\exp134_part123\groups.grp"   -u false -i true -o false -p true -e -x -s
	python MExtract.py -m AllExtract -g "E:\exp134_part123\99-112.grp"  -l "E:\exp134_part123\groups.grp"   -u false -i true -o false -p true -e -x -s
	python MExtract.py -m AllExtract -g "E:\exp134_part123\113-126.grp" -l "E:\exp134_part123\groups.grp"   -u false -i true -o false -p true -e -x -s
	python MExtract.py -m AllExtract -g "E:\exp134_part123\127-140.grp" -l "E:\exp134_part123\groups.grp"   -u false -i true -o false -p true -e -x -s
	python MExtract.py -m AllExtract -g "E:\exp134_part123\141-154.grp" -l "E:\exp134_part123\groups.grp"   -u false -i true -o false -p true -e -x -s
	python MExtract.py -m AllExtract -g "E:\exp134_part123\155-162.grp" -l "E:\exp134_part123\groups.grp"   -u false -i true -o false -p true -e -x -s


echo .
echo .
echo finished