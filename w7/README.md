# **Eigenvalue & Eigenvector**

### **Task:**
Generate a 11x11 **symmetric matrix** A by using random number generator (Gaussian distribution with mean=0 and standard deviation=1.0). Then, compute all eigenvalues and eigenvectors of A. Print the eigenvalues and their corresponding eigenvectors in descending order.

Functions:
- jacobi(): Obtain eigenvalues using Jacobi transformation
- eigsrt(): Sort results of jacobi()

<br/>

### **Compilation Method**

All files are compiled using GCC with the following input at command prompt.

File: hw7.c, ran1.c, gasdev.c, jacobi.c, eigsrt.c, nrutil.c 

    gcc hw7.c -o hw7.o -c
    gcc ran1.c -o ran1.o -c
    gcc gasdev.c -o gasdev.o -c
    gcc jacobi.c -o jacobi.o -c
    gcc eigsrt.c -o eigsrt.o -c
    gcc nrutil.c -o nrutil.o -c

The compiled codes can be executed with the following format in command prompt.

    ./hw7