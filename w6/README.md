# **Simulation**

### **Task:**
Programming on Random Number Generation
1. Uniform distribution in [a,b]
2. Gaussian distribution with mean=m, standard deviation=s

Generate 1000 samples and draw a histogram; 100 intervals for each distribution (a=-3, b=2, m=0.5, s=1.5). 

Repeat the same task by varying the number of samples (i.e., 100, 10000, 100000).

<br/>

### **Compilation Method**

All files are compiled using GCC with the following input at command prompt.

File: hw6.c, ran1.c, gasdev.c

    gcc hw6.c -o hw6.o -c
    gcc ran1.c -o ran1.o -c
    gcc -o hw6 hw6.o ran1.o -lm
    gcc gasdev.c -o gasdev.o -c
    gcc -o hw6 hw6.o gasdev.o ran1.o -lm

The compiled codes can be executed with the following format in command prompt.

    ./hw6