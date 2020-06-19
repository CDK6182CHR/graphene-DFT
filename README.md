# graphene-DFT
This is the source code for project of *Material Physics A*, an undergraduate course in Nanjing University. This project calculates **2D** band structure of graphene with DFT &amp; CA-LDA, expand wavefunction with OPW basic-set, programmed with C++. 

To run the project, you need GNU Scientific Library, known as GSL, of which document is available at https://www.gnu.org/software/gsl/doc/html/index.html

Your C++ compiler must support C++11 or higher standard, for using of C++11 features like lambda expression, right value reference etc.

I have provided `Makefile` to build the project with GNU compiler, however you can also use other tools. To build and run this project with GNU compiler, just type 

```shell
$ make
$ ./graphene.exe
```

You may need to change something in `Makefile`, according to your platform. 

The band structure data is stored in text file `EIGENVAL.txt`, you can run `plot.py`  to see the figure, which requires Python3 interpreter. 

```shell
$ python3 plot.py
```

I have built and run the project successfully with following compilers, on Microsoft Windows 10 (64 bits):

- Microsoft Visual Studio 2019 (Community)
- GNU g++ 8.1.0 (MinGW64)

To change some arguments, like number of cells, number of basic sets, look at `base.h`  and `base.cpp`.



## Note

Till now, the outputs seems not to be very well, which means, the band structure plotted by this project is quite different from real band structure of graphene. I haven't got whether the errors from the 2D model and OPW approximation, or some mistakes I have made in the C++ codes leads to the failure.