set terminal png
set output 'Sol.png'
set ticslevel 0
splot 'Result/sol0.dat', 'Result/sol1.dat', 'Result/sol2.dat', 'Result/sol3.dat', 'Result/sol4.dat', 'Result/sol5.dat', sin(x)+cos(y)
