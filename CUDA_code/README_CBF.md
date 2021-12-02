## Compiling and running standalone coherent beamformer

Before compiling and running, if standalone, uncomment main() function at the bottom of the coherent_beamformer_char_in.cu script.
Change sm architecture depending on GPU e.g. for NVIDIA A4000, sm_86 is appropriate.

# To compile
```
nvcc -o coherent_beamformer_char_in.exe -arch=sm_86 coherent_beamformer_char_in.cu
```

# To run
```
./coherent_beamformer_char_in.exe
```

# Generate shared object library
To generate shared object file for hashpipe plugins, comment main() function in .cu script and compile with Makefile using
```
make
```

And install the library and header with:
```
make install
```

Don't forget to change the paths in the Makefile if necessary, and add current library path to LD_LIBRARY_PATH
