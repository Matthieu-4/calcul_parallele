set terminal png
set output 'Sol.png'
set ticslevel 0
splot 'Result/sol0.dat', 'Result/sol1.dat', 'Result/sol2.dat', x*(1-x)*y*(1-y)
