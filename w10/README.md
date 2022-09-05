# **Nonlinear Data Fitting**

### **Task:**
Model: 2D transformation between given images (C++ implementation)

![](tmp.PNG)

1. Establish the feature correspondence between the two images. (Output: x<sub>i</sub>, y<sub>i</sub>, x<sub>i</sub>', y<sub>i</sub>', i=1~N)
2. Find the parameters using the correspondence data.

<br/>

### **Compilation Method**

The following file is compiled using G++ with the following input at command prompt.

File: hw10.cpp
    
    g++ hw10.cpp -o hw10 `pkg-config --cflags --libs opencv4`

The compiled code can be executed with the 

    ./hw10