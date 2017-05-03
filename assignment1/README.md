With make and alias:
1. Call ``make clean && make`` in root directory 
2. Run with ``pda``

Or without make:
1. Compile in root directory with ``mpicc -o assignment src/*.c -I ./include``
2. Run with ``mpiexec -np [number of processes] ./assignment [program parameters]``

