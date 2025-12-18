# MMCACovid19Vac

MMCACovid19Vac is a Julia package for simulating epidemic dynamics in a metapopulation, with built-in support for modeling vaccination campaigns. The package is based on two previous implementations:
* MMCAcovid19 Julia package [https://docs.juliahub.com/MMCAcovid19/]([url](https://docs.juliahub.com/MMCAcovid19/))
* MMCA_with_vaccination [https://github.com/PGcastioni/MMCA_with_vaccination/]([url](https://github.com/PGcastioni/MMCA_with_vaccination/))

# Brief history

## MMCAcovid19
This package [MMCAcovid19](https://github.com/jtmatamalas/MMCAcovid19), written in the [Julia](https://julialang.org) language, implements the epidemic model for COVID-19 developed by a group of researchers from [Universitat Rovira i Virgili](https://www.urv.cat) and [Universidad de Zaragoza](http://unizar.es) [[1](#References-1)]. The model makes use of a Microscopic Markov Chain Approach (MMCA) to describe mathematically the dynamics of a so-called metapopulation model of epidemic spreading [[2-4](#References-1)]. 

## MMCA_with_vaccination
The original was later extended by Piergiorgio Castioni to allow modelling of the effect of vaccines, herd immunity and reinfections. Each of the previous compartments has been triplicated to account for the fact that there are vaccinated and unvaccinated individuals. Having received the vaccine changes some of the transition probabilities related to the most negative aspects of the disease, such as transmission, hospitalization and death.

## MMCACovid19Vac

This project is an attempt to extend the functionalities of the previous version. Currently, this version provides a standard configuration format for defining a model and provides a simple set of command line scripts to generate configuration templates and run the model.

1. Alex Arenas, Wesley Cota, Jesús Gómez-Gardeñes, Sergio Gómez, Clara Granell, Joan T. Matamalas, David Soriano-Paños and Benjamin Steinegger: Modeling the spatiotemporal epidemic spreading of COVID-19 and the impact of mobility and social distancing interventions, _Physical Review X_ **10** (2020) 041055 ([doi](https://doi.org/10.1103/PhysRevX.10.041055))

2. Sergio Gómez, Alex Arenas, Javier Borge-Holthoefer, Sandro Meloni and Yamir Moreno: Discrete-time Markov chain approach to contact-based disease spreading in complex networks, _Europhysics Letters_ **89** (2010) 38009 ([doi](https://doi.org/10.1209/0295-5075/89/38009))

3. Jesús Gómez-Gardeñes, David Soriano-Paños and Alex Arenas: Critical regimes driven by recurrent mobility patterns of reaction-diffusion processes in networks, _Nature Physics_ **14** (2018) 391–395 ([doi](https://doi.org/10.1101/2020.03.21.20040022))

4. David Soriano-Paños, L. Lotero, Alex Arenas and Jesús Gómez-Gardeñes: Spreading processes in multiplex metapopulations containing different mobility networks, _Physical Review X_ **8** (2018) 031039 ([doi](https://doi.org/10.1103/PhysRevX.8.031039))
