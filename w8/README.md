# **Interpolation**

### **Task:**
Image Resampling (C++ implementation)
- Read an image file and identify the resolution and resample it to a specified resolution
- Input: Target resolution (M'xN')
- Output: Resampled image
- Method: Bilinear interpolation

<br/>

### **Compilation Method**

The following file is compiled using G++ with the following input at command prompt.

File: hw8.c

    g++ hw8.cpp -o hw8 `pkg-config --cflags --libs opencv4`

The compiled code can be executed with the respective format in command prompt.

    ./hw8