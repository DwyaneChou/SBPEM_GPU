# SBPEM_GPU
This is the GPU version of SBPEM, the model framework is the same as SBPEM, but the calculation will be implemented in GPU.

For a coarse resolution CPU version will fast than GPU version, but when the resolution gets finer, GPU version will be much more faster than CPU version. 

In instance, my PC

CPU:I7-4770K

GPU:NVIDIA GeForce GTX 970

In my test, when I defined the resolution to 2 degrees and time_step=18 , the CPU version needed about 3s to run 1 hour, but the GPU version needed over 12s,

After changing resolution to 0.5 degree and time_step=1, the CPU version needed over 1800s to run 1hour, but the GPU version needed only about 239s.

That is, GPU needs the large matrices to make full use of its strength.
