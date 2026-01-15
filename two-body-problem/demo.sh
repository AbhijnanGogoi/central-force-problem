gcc two-body-prob.c -o two-body-prob.o -lm
./two-body-prob.o --m1=2 --m2=2 --r=100 --v=1 --vt=0.0005 --k=-1e3 --n=-1 --steps=100000000 --step_size=1e-5 --read_steps=1000 --name=demo
python two-body-prob-plot.py --name=demo
