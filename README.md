# p-median

# Description

This is a framework in `C++` to solve the **p-median problem**.

# Build

You can use the make file

```console
$ cd project_path 
$ make
```

Or use command line compilation

```console
$ g++ -static -std=c++11 -static-libgcc -static-libstdc++ -O3 *.cpp -o run
```

It is also recomended the use of the flags: `-Wall` `-Werror` `-Wextra` `-pedantic`

# Usage

After compile the project in the `run` file, you can use as follows:

```console
$ run P N F FLAG_LIST
```

where:

* _P_ the parameter _p_ of the p-median problem.
* _N_ the number of customers.
* _F_ the number of facilities.
* _FLAG_LIST_ the rest of the parameters to configurate the instance. The parameter options are:
  * `--input ` the inputs files.
  * `--output ` path to the output files.
  * `--no-output ` if it is necessary to write output files. Default value is TRUE.
  * `--debug ` print the logs. Default value is FALSE.
  * `--seed `_`NUMBER`_ seed to the random number generator.
  * `--algo `_`OPT`_ algorithm to execute, the implemented algorithms are: `VNS`, `GA`.
  * `--pop `_`SIZE`_ population size.
  * `--iter ` maximum number of iterations
  * `--neighborhood `_`OPT`_ neighbourhood strategy.
  * `--k `_`SIZE`_ neighbourhood size.
  * `--gen `_`OPT`_  process to generate the initial population.
  * `--kmax `_`NUMBER`_ parameter of the VNS.
  * `--Kmayus `_`NUMBER`_ parameter of the VNS.
  * `--lambda `_`NUMBER`_ parameter of the GA.
  * `--m `_`NUMBER`_ parameter of the VNS.
  * `--accept `_`OPT NUMBER`_ parameter of the VNS.
  * `--next `_`OPT NUMBER`_ parameter of the VNS.
  * `--change `_`OPT`_ parameter of the VNS.
  * `--shake `_`OPT`_ parameter of the VNS.
  * `--ls1 `_`OPT PARAM`_ local search 1.
  * `--ls2 `_`OPT PARAM`_ local search 2.
  * `--ga-sel `_`OPT`_ selection operator of the GA.
  * `--ga-cross `_`OPT`_ crossover operator of the GA.
  * `--ga-mut `_`OPT PROB`_  mutation operator of the GA.
  * `--ga-repl `_`OPT`_ replacement operator of the GA.

# Publications

Cintrano C., Chicano F., Stützle T., Alba E. (2018) Studying Solutions of the p-Median Problem for the Location of Public Bike Stations. Advances in Artificial Intelligence. CAEPIA 2018. Lecture Notes in Computer Science, vol 11160. Springer, Cham. DOI: [https://doi.org/10.1007/978-3-030-00374-6_19](https://doi.org/10.1007/978-3-030-00374-6_19)

# Thanks

* @mbcrawfo for the [MakeFile example](https://github.com/mbcrawfo/GenericMakefile).
