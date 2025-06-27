---
title: 'MieSimulatorGUI: A Mie Simulation tool to predict light scattering with spherical particles'

tags:
- Mie scattering
- Mie theory
- scattering
- spherical particles
  authors:
- name: Janaka Ranasinghesagara
  orcid: 0000-0002-6069-6527
  equal-contrib: false
  affiliation: "1, 2" # (Multiple affiliations must be quoted)
- name: Vasan Venugopalan
  orcid: 0000-0003-4781-1049
  equal-contrib: false
  affiliation: "1, 2" # (Multiple affiliations must be quoted)
  affiliations:
  - name: Department of Chemical and Biomolecular Engineering, University of California, Irvine, California 92697, United States
    index: 1
  - name: Beckman Laser Institute and Medical Clinic, University of California, Irvine, California 92697, United States
    index: 2
    date: 27 June 2025
    bibliography: paper.bib

---

[comment]: https://joss.theoj.org/]

# Summary

Mie theory is a mathematical framework derived from Maxwell's equations that models electromagnetic scattering by spherical dielectric particles.  Mie theory predictions are used in a wide variety of scientific and engineering fields spanning fields such as atmospheric optics, particle characterization, computer graphics, nanofluids, remote sensing, and biomedical optics. We developed the `MieSimulatorGUI` to provide rapid computation of Mie theory predictions while using a user-friendly interface to perform complex Mie calculations. By utilizing a C/C++ computational engine with the Qt framework, our cross-platform tool calculates key optical properties like scattering coefficients, angular scattering distributions, and scattering anisotropy for both monodisperse and polydisperse particle distributions. With six interactive panels, the GUI enables users to specify input parameters, visualize particle distributions, and analyze scattering phenomena, including the fitting of spectrally-varying reduced scattering coefficients parameters, which is particularly valuable in fields like tissue optics. The `MieSimulatorGUI` is open source, hosted on [GitHub](https://github.com/VirtualPhotonics/MieSimulatorGUI), available from a [download page](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki/Downloads) and distributed under permissive [MIT license](https://opensource.org/license/mit).

# Statement of need

Mie theory, derived from Maxwell's equations, provides a comprehensive framework for modeling electromagnetic scattering by spherical particles, proving invaluable in fields ranging from nanomaterials and biomedical optics to atmospheric science and astronomy `[@Goody1989; @Saidi1995; @Wang2005; @Horvath2009]`. Despite its broad applicability, the theory's reliance on complex mathematical constructs, such as infinite series and special functions`[@VandeHulst1957; @Bohren1983; @Majic2020]`, demands advanced computational implementations, as exemplified by [SCATTPORT.org](https://scattport.org). Given its mathematical complexity, a user-friendly Graphical User Interface (GUI) is essential to simplify its application, enabling researchers and engineers to rapidly compute critical metrics associated with scattering phenomena. By combining a powerful C/C++ computational engine with [Qt](https://www.qt.io/) intuitive graphical interface, the `MieSimulatorGUI` provides user-friendly access to complex Mie theory results for streamlined analysis. 

# Main Features

The `MieSimulatorGUI` employs the `BHMIE` code `[@Bohren1983]` and procedures described by `[@Wiscombe1979]` to simulatate light scattering by spheres. This user-friendly tool calculates spectral dependencies of key optical properties including scattering, extinction and back scattering cross sections, scattering coefficient, reduced scattering coefficient, scattering amplitude matrix entries (S1, S2), phase function, average cosine of the phase function, and forward / backward scattering percentage `[@Bohren1983]`. These calculations can be performed for both monodisperse and polydisperse particle distributions `[@Gelebart1996]`. Moreover, the tool facilitates parameter estimation by fitting reduced scattering coefficient curves, a technique of significant value in tissue optics`[@Jacques2013]`. 

# Design and Functionality

The `MieSimulatorGUI` is a cross-platform GUI application developed using the C/C++ programming language within the [Qt](https://www.qt.io/) framework. The GUI contains six interactive panels as shown in Fig. \autoref{fig:Figure1}. 

![Figure 1: Six interactive panels of `MieSimulatorGUI`: (a) Input selection, (b) Number density, (c) µ<sub>s</sub>, (d) Phase function, (e) µ<sub>s</sub>',  and (f) Anisotropy \label{fig:Figure1}](figure1.png){ width=100% }
Figure 1: Six interactive panels of `MieSimulatorGUI`: (a) Input selection, (b) Number density, (c) Scattering coefficient, (d) Phase function, (e) Reduced Scattering,  and (f) Scattering Asymmetry

## Input Selection Panel

This panel provides the user with the ability to specify inputs necessary to define a Mie simulation for single wavelength or series of wavelengths. The distribution of spheres is described using sphere concentration (spheres/mm³) or volume fraction. The volume fraction, the proportion of space occupied by spheres, is determined by the product of sphere concentration and individual sphere volume. Currently, sphere diameters between 0.1 nm and 300 μm and wavelengths between 50 nm and 3000 nm (3 µm) are supported. 
For absorbing spheres, the complex refractive index (m<sub>sphere</sub>) is defined as m<sub>real</sub> – j m<sub>imag</sub>, where m<sub>real</sub> and m<sub>imag</sub> represent real and imaginary components, respectively `[@VandeHulst1957; @Wiscombe1979]`. 

This tool allows users to choose between mono-disperse particle (uniform-sized spheres) and poly-disperse particle (varied-sized spheres) distribution, depending on the application. Mono-disperse distributions restrict analysis to spheres with identical size and refractive index. In contrast, poly-disperse distributions enable simulations of spheres with diverse attributes and supports three distribution models: Log-normal, Gaussian, and user-defined (custom). In this selection, users have the freedom to specify different refractive indexes for different spheres. 

## Number Density Panel

This panel graphically presents the number of spheres (N<sub>s</sub>) used in the simulation.  The subsequent tab displays the 'Size Parameter' defined as 2πR n<sub>med</sub> / λ<sub>vacuum</sub> `[@Bohren1983]`, where R denotes the particle radius, n<sub>med</sub>  the medium's refractive index, and λ<sub>vacuum</sub>  wavelength in vacuum. 

## Scattering Coefficient Panel

The Mie calculations provide three important efficiency factors: the scattering efficiency (Q<sub>sca</sub>), the extinction efficiency (Q<sub>ext</sub>), and the back-scattering efficiency (Q<sub>back</sub>) ([Mie Scattering Efficiencies](https://miepython.readthedocs.io/en/latest/02_efficiencies.html)). These dimensionless quantities combined with the particle's cross section area (πR<sup>2</sup>) yeild the corresponding scattering cross section (C<sub>sca</sub>), extinction cross section (C<sub>ext</sub>) and back-scattering cross section (C<sub>back</sub>). The calculated cross sections are displayed across three separate tabs. 
For mono-disperse distribution, the scattering coefficient (µ<sub>s</sub>)  is simply the product of the scattering cross-section and N<sub>s</sub>.  In poly disperse distribution, however, the scattering coefficient is calculated using the discrete particle model detailed in Schmitt and Kumar `[@Schmitt1998]`. 

## Phase Function Panel

The phase function, representing the angular distribution of scattered light, is displayed in this panel using polar and linear plots. These plots are derived from the complex amplitude scattering matrix elements, S<sub>1</sub>, and  S<sub>2</sub> , which describe the transformation of incident electromagnetic field to far-field scattered field `[@VandeHulst1957; @Bohren1983]`. The Phase function or S<sub>1</sub> and S<sub>2</sub> data at a specific wavelength can be found by using the wavelength slider.

## Reduced Scattering Panel

This panel shows the reduced scattering coefficient (µs'), calculated as the product of the scattering coefficient (µs) and (1-g). µs' is crucial in various fields, particularly in biomedical optics, because it allows for the non-invasive quantification of tissue properties. Users can use µs' 'Power Law Fitting' tab to compute the fitting parameters given in Steve L. Jacques's review paper `[@Jacques2013]`.

## Scattering Asymmetry Panel

The scattering asymmetry (anisotropy) panel displays the directional properties of the scattering phase function. The first tab presents the average cosine of the phase function (g), which quantifies the average scattering angle and indicates forward or backward scattering dominance. The second tab provides the integrated forward and backward scattering fractions, offering a detailed analysis of the angular scattering distribution.

  

# Acknowledgements

We would like to acknowledge Dr. Carole Hayakawa and Lisa Glover for their support and guidance during the code deployment to GitHub.

# References
