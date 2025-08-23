for rep in `seq 1 10`;
do
	python dfe.py 8 1.0 1.0
	python dfe.py 20 1.0 1.0
	python dfe.py 8 0.05 0.01
	python dfe.py 20 0.05 0.01
done