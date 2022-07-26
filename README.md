# QMC_SSE
stochastic series expansion(QMC) for 2d Heisenberg model

It's a naive julia version of Anders Sandvik's code(http://physics.bu.edu/~sandvik/trieste12/ssebasic.f90)

I am new to julia, so the code may seem unmature. 

Amazingly, my tests show that the julia version code is almost good as the fortran version code(for small lattice size lx=ly=8,16)!

For lx=ly=64,beta=128,nbins=10,msteps=10000,isteps=10000, the fortran code(Anders Sandvik) costs about 11619 second, 
while the julia code costs about 14038 second.

Any comments are welocome, particularly the advice about improving the code's performance .

email 202021140008@mail.bnu.edu.cn


