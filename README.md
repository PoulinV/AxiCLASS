AxiCLASS is based on CLASS by Julien Lesgourgues and Thomas Tram

AxiCLASS authors: Vivian Poulin, Tristan Smith, Tanvi Karwal
the code is available at https://github.com/PoulinV/AxiCLASS



-------------------------------------------------------

AxiCLASS is a code based on CLASS to compute the effect of Axion-like
particles onto the CMB and linear matter power spectra, as well as related
observables.

It models a field with a potential

V(phi)= m^2f^2(1-cos(phi/f))^n

where m [H0] is the field mass, f [mpl] the field decay constant, phi the field value, 
and n an exponent [1,infty] that controls the equation of state of the field (i.e.
w = p_phi/rho_phi = (n-1)/(n+1)) when it dilutes.

The code can be run with the exact (scalar field) dynamics, as well as a fluid
approximation, similar to the "axionCAMB" code by Dan Grin https://github.com/dgrin1/axionCAMB

For information on how to run the code in either "mode" see the files:

example_axiCLASS.ini
example_axiCLASS_fld.ini


we provide a ".param" file and some bestfit / covmat in the folder 
montepython_param_files

for any questions with the code, feel free to write to 

vivian.poulin@umontpellier.fr

If you use the code, (on top of standard CLASS citations), 
please cite the following works (to which we refer for all relevant definitions)


@article{Smith:2019ihp,
    author = "Smith, Tristan L. and Poulin, Vivian and Amin, Mustafa A.",
    title = "{Oscillating scalar fields and the Hubble tension: a resolution with novel signatures}",
    eprint = "1908.06995",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1103/PhysRevD.101.063523",
    journal = "Phys. Rev. D",
    volume = "101",
    number = "6",
    pages = "063523",
    year = "2020"
}


@article{Poulin:2018dzj,
    author = "Poulin, Vivian and Smith, Tristan L. and Grin, Daniel and Karwal, Tanvi and Kamionkowski, Marc",
    title = "{Cosmological implications of ultralight axionlike fields}",
    eprint = "1806.10608",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1103/PhysRevD.98.083525",
    journal = "Phys. Rev. D",
    volume = "98",
    number = "8",
    pages = "083525",
    year = "2018"
}


and consider citing related works:

@article{Poulin:2018cxd,
    author = "Poulin, Vivian and Smith, Tristan L. and Karwal, Tanvi and Kamionkowski, Marc",
    title = "{Early Dark Energy Can Resolve The Hubble Tension}",
    eprint = "1811.04083",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1103/PhysRevLett.122.221301",
    journal = "Phys. Rev. Lett.",
    volume = "122",
    number = "22",
    pages = "221301",
    year = "2019"
}

@article{Murgia:2020ryi,
    author = "Murgia, Riccardo and Abell\'an, Guillermo F. and Poulin, Vivian",
    title = "{Early dark energy resolution to the Hubble tension in light of weak lensing surveys and lensing anomalies}",
    eprint = "2009.10733",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1103/PhysRevD.103.063502",
    journal = "Phys. Rev. D",
    volume = "103",
    number = "6",
    pages = "063502",
    year = "2021"
}

Thanks a lot for using the code!
It is still consider "work in progress" so, our apologies if it is not as user-friendly and bug-free as we'd like.


-------------------updates------------------------

2022-03-15

Found a bug in perturbation when the user request switching from scf to fluid. the code seem to get stuck.
Probably an issue with tables handling, to be investigated. 
For now, we encourage the user to use either "full scalar field" or "full fluid", as presented in the 
example files.


2022-02-28

* Improving readability of the folder. 

-- Moved scripts in python_scripts folder. Most of the scripts are old and deprecated, need to be removed / updated.
-- Created a folder with MP param files / covmat / bestfit / conf file.

* Updated the log10_fraction_axion_ac parameter to fraction_axion_ac to ease the compatibility with montepython, which perform searches on "fraction_axion_ac".

-- updated the .ini file
-- added an ini file based on SPT+ACT+PlanckTT650TEEE run

* Implemented the Gamma function in background.h, removed dependency to GSL.

* to do: create a notebook / scripts.


--------known issues with compilation-------------

2022-02-28
The code used to make use of GSL to compute the Gamma function. This was creating problem sometimes when "cythoning" classy if clang was used.
We now implemented the Gamma function explicitely in background.h. Many thanks to Pierre Zhang for finding out this bug, and helping with debugging.

------- THE REST IS THE STANDARD CLASS README -------


CLASS: Cosmic Linear Anisotropy Solving System  {#mainpage}
==============================================

Authors: Julien Lesgourgues, Thomas Tram, Nils Schoeneberg

with several major inputs from other people, especially Benjamin
Audren, Simon Prunet, Jesus Torrado, Miguel Zumalacarregui, Francesco
Montanari, Deanna Hooper, Samuel Brieden, Daniel Meinert, Matteo Lucca, etc.

For download and information, see http://class-code.net

Compiling CLASS and getting started
-----------------------------------

(the information below can also be found on the webpage, just below
the download button)

Download the code from the webpage and unpack the archive (tar -zxvf
class_vx.y.z.tar.gz), or clone it from
https://github.com/lesgourg/class_public. Go to the class directory
(cd class/ or class_public/ or class_vx.y.z/) and compile (make clean;
make class). You can usually speed up compilation with the option -j:
make -j class. If the first compilation attempt fails, you may need to
open the Makefile and adapt the name of the compiler (default: gcc),
of the optimization flag (default: -O4 -ffast-math) and of the OpenMP
flag (default: -fopenmp; this flag is facultative, you are free to
compile without OpenMP if you don't want parallel execution; note that
you need the version 4.2 or higher of gcc to be able to compile with
-fopenmp). Many more details on the CLASS compilation are given on the
wiki page

https://github.com/lesgourg/class_public/wiki/Installation

(in particular, for compiling on Mac >= 10.9 despite of the clang
incompatibility with OpenMP).

To check that the code runs, type:

    ./class explanatory.ini

The explanatory.ini file is THE reference input file, containing and
explaining the use of all possible input parameters. We recommend to
read it, to keep it unchanged (for future reference), and to create
for your own purposes some shorter input files, containing only the
input lines which are useful for you. Input files must have a *.ini
extension. We provide an example of an input file containing a
selection of the most used parameters, default.ini, that you may use as a
starting point.

If you want to play with the precision/speed of the code, you can use
one of the provided precision files (e.g. cl_permille.pre) or modify
one of them, and run with two input files, for instance:

    ./class test.ini cl_permille.pre

The files *.pre are suppposed to specify the precision parameters for
which you don't want to keep default values. If you find it more
convenient, you can pass these precision parameter values in your *.ini
file instead of an additional *.pre file.

The automatically-generated documentation is located in

    doc/manual/html/index.html
    doc/manual/CLASS_manual.pdf

On top of that, if you wish to modify the code, you will find lots of
comments directly in the files.

Python
------

To use CLASS from python, or ipython notebooks, or from the Monte
Python parameter extraction code, you need to compile not only the
code, but also its python wrapper. This can be done by typing just
'make' instead of 'make class' (or for speeding up: 'make -j'). More
details on the wrapper and its compilation are found on the wiki page

https://github.com/lesgourg/class_public/wiki

Plotting utility
----------------

Since version 2.3, the package includes an improved plotting script
called CPU.py (Class Plotting Utility), written by Benjamin Audren and
Jesus Torrado. It can plot the Cl's, the P(k) or any other CLASS
output, for one or several models, as well as their ratio or percentage
difference. The syntax and list of available options is obtained by
typing 'pyhton CPU.py -h'. There is a similar script for MATLAB,
written by Thomas Tram. To use it, once in MATLAB, type 'help
plot_CLASS_output.m'

Developing the code
--------------------

If you want to develop the code, we suggest that you download it from
the github webpage

https://github.com/lesgourg/class_public

rather than from class-code.net. Then you will enjoy all the feature
of git repositories. You can even develop your own branch and get it
merged to the public distribution. For related instructions, check

https://github.com/lesgourg/class_public/wiki/Public-Contributing

Using the code
--------------

You can use CLASS freely, provided that in your publications, you cite
at least the paper `CLASS II: Approximation schemes <http://arxiv.org/abs/1104.2933>`. Feel free to cite more CLASS papers!

Support
-------

To get support, please open a new issue on the

https://github.com/lesgourg/class_public

webpage!
