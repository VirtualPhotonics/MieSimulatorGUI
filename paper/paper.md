---
title: >-
  MieSimulatorGUI: A user-friendly tool to compute and visualize light
  scattering by spherical dielectric particles
tags:
  - Mie scattering
  - Lorenz-Mie theory
  - Scattering
  - Spherical particles
  - Scattering efficiency
  - Phase function
authors:
  - name: Janaka C. Ranasinghesagara
    orcid: 0000-0002-6069-6527
    affiliation: '1, 2'
  - name: Carole Hayakawa
    orcid: 0000-0001-5696-160X
    affiliation: '1, 2'
  - name: Lisa Glover
    orcid: 0009-0009-0291-5099
    affiliation: 2
  - name: Vasan Venugopalan
    orcid: 0000-0003-4781-1049
    corresponding: true
    email: vvenugop@uci.edu
    affiliation: '1, 2'
affiliations:
  - name: >-
      Department of Chemical and Biomolecular Engineering, University of
      California, Irvine, California 92697, United States
    index: 1
  - name: >-
      Beckman Laser Institute and Medical Clinic, University of California,
      Irvine, California 92697, United States
    index: 2
date: 8 July 2025
bibliography: paper.bib

---

[comment]: https://joss.theoj.org/]

# Summary

Mie theory is a mathematical framework derived from Maxwell's equations that models electromagnetic scattering by spherical dielectric particles.  Mie theory predictions are used in a wide variety of scientific and engineering fields spanning fields such as atmospheric optics, particle characterization, computer graphics, nanofluids, remote sensing, and biomedical optics. We developed the `MieSimulatorGUI` to provide a user-friendly interface for users to perform and visualize computations of optical scattering by particles using Mie theory. By utilizing a C/C++ computational engine with the Qt framework, our cross-platform tool calculates key optical properties for instance scattering coefficients, angular scattering distributions, and scattering asymmetry for both monodisperse and polydisperse particle distributions. With six interactive panels, the GUI enables users to specify input parameters, visualize particle distributions, and analyze scattering phenomena, including the fitting of spectrally-varying reduced scattering coefficient parameters, which is particularly valuable in fields such as tissue optics. The `MieSimulatorGUI` is open source, hosted on [GitHub](https://github.com/VirtualPhotonics/MieSimulatorGUI), available from a [download page](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki/Downloads) and distributed under permissive [MIT license](https://opensource.org/license/mit).

# Statement of need

Mie theory, derived from Maxwell's equations, provides a comprehensive framework for modeling electromagnetic scattering by spherical particles, proving invaluable in fields ranging from nanomaterials and biomedical optics to atmospheric science and astronomy `[@Goody1989; @Saidi1995; @Wang2005; @Horvath2009]`. Despite its broad applicability, the theory's reliance on complex mathematical constructs, such as infinite series and special functions`[@VandeHulst1957; @Bohren1983; @Majic2020]`, demands advanced computational implementations, as exemplified by [SCATTPORT.org](https://scattport.org). Given its mathematical complexity, a user-friendly Graphical User Interface (GUI) is essential to simplify its application, enabling researchers and engineers to rapidly compute critical metrics associated with scattering phenomena. By combining a powerful C/C++ computational engine with [Qt](https://www.qt.io/) intuitive graphical interface, the `MieSimulatorGUI` provides user-friendly access to complex Mie theory computations for streamlined analysis. 

# Main Features

The `MieSimulatorGUI` employs the `BHMIE` code `[@Bohren1983]` and procedures described by `[@Wiscombe1979]` to simulate light scattering by spherical dielectric particles. This user-friendly tool calculates spectral dependencies of key optical properties including scattering, extinction and back scattering cross sections, scattering coefficient, reduced scattering coefficient, scattering amplitude matrix entries ($\text{S}_1$, $\text{S}_2$), phase function, average cosine of the phase function, and forward / backward scattering percentage `[@VandeHulst1957; @Bohren1983]`. These calculations can be performed for both monodisperse and polydisperse particle distributions `[@Gelebart1996]`. Moreover, the tool facilitates parameter estimation by fitting reduced scattering coefficient curves, a technique of significant value in tissue optics`[@Jacques2013]`. All relevant documentation for this tool is available on the [Mie Simulator GUI Wiki](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki) page.

# Design and Functionality

The `MieSimulatorGUI` is a cross-platform GUI application developed using the C/C++ programming language within the [Qt](https://www.qt.io/) framework. The GUI contains six interactive panels as shown in \autoref{fig:Figure1}. 

![Figure 1: Six interactive panels of `MieSimulatorGUI`: (a) Input selection, (b) Particle size distribution, (c) Scattering coefficient, (d) Phase function, (e) Reduced Scattering,  and (f) Scattering Asymmetry \label{fig:Figure1}](Figure1.png)
Figure 1: Six interactive panels of `MieSimulatorGUI`: (a) Input selection, (b) Particle size distribution, (c) Scattering coefficient, (d) Phase function, (e) Reduced Scattering, and (f) Scattering Asymmetry

## Input Selection Panel

This panel provides the user with the ability to specify inputs necessary to define a Mie simulation for single wavelength or series of wavelengths. The distribution of spheres is described using sphere concentration ($\text{spheres/mm}^3$) or volume fraction. The volume fraction, the proportion of space occupied by spheres, is determined by the product of sphere concentration in 1 $\text{mm}^3$ and individual sphere volume. Sphere diameters are presently restricted to 0.1 $\text{nm}$ and 300 $\text{µm}$ and wavelengths to 50 $\text{nm}$ and 3000 $\text{nm}$ (3 $\text{µm}$ ). For absorbing spheres, the complex refractive index ($m_{sphere}$) is defined as $m_{real}$ – j $m_{imag}$, where $m_{real}$ and $m_{imag}$ represent real and imaginary components, respectively `[@VandeHulst1957; @Wiscombe1979]`. 

This tool allows users to choose between mono-disperse (uniform-sized) or poly-disperse (variable-sized) particle distributions. Mono-disperse distributions restrict analysis to spheres with uniform size and refractive index. In contrast, poly-disperse distributions enable simulations of spheres with diverse attributes and support three size distribution models: 1. Log-normal, 2. Gaussian, and 3. user-defined (custom). The user-defined option allows the user to specify different refractive indexes for different spheres, following the examples in the `CustomDataSamples` folder. 

## Number Density Panel

This panel graphically presents the number density of spheres $N_s$ [\#/ $\text{mm}^3$] used in the simulation.  The subsequent tab displays the 'Size Parameter' defined as $\text{2π}Rn_{med} / λ_{vacuum}$ `[@Bohren1983]`, where $R$ [ $\text{µm}$ ] denotes the particle radius, $n_{med}$  the medium's refractive index, and $λ_{vacuum}$ [ $\text{µm}$ ] is the wavelength in vacuum. 

## Scattering Coefficient Panel

The Mie calculations provide three important efficiency factors: the scattering efficiency ($Q_{sca}$), the extinction efficiency ($Q_{ext}$), and the back-scattering efficiency ($Q_{back}$) ([Mie Scattering Efficiencies](https://miepython.readthedocs.io/en/latest/02_efficiencies.html)). These dimensionless quantities combined with the particle's cross section area ($\text{π}R^2$) yield the corresponding scattering cross section $C_{sca}$ [ $\text{/mm}^2$ ], extinction cross section $C_{ext}$ [ $\text{/mm}^2$ ] and back-scattering cross section $C_{back}$ [ $\text{/mm}^2$ ] `[@VandeHulst1957; @Bohren1983]`. The calculated cross sections are displayed across three separate tabs. 
For mono-disperse distribution, the scattering coefficient ($µ_s$)  is simply the product of the scattering cross-section $C_{sca}$ and the number density $N_s$.  For poly-disperse distributions the scattering coefficient is calculated using the discrete particle model detailed in Schmitt and Kumar `[@Schmitt1998]`. 

## Phase Function Panel

The phase function represents the angular distribution of scattered light.  The calculated results for the phase function are displayed in this panel using polar and linear plots. These plots are derived from the complex amplitude scattering matrix elements, $\text{S}_1$, and $\text{S}_2$, which describe the transformation of incident electromagnetic field to far-field scattered field `[@VandeHulst1957; @Bohren1983]`. The wavelength slider enables the user to visualize the phase function or $\text{S}_1$ and $\text{S}_2$ data at a specific wavelength.

## Scattering Asymmetry Panel

The scattering asymmetry (anisotropy) panel displays the directional properties of the scattering phase function. The first tab presents the average cosine of the phase function ( $g$ ), which quantifies the average scattering angle which indicates the prevalence of forward ($g>0$) vs backward scattering ($g<0$). The second tab provides the integrated forward and backward scattering fractions, offering a detailed analysis of the angular scattering distribution.

## Reduced Scattering Panel

This panel shows the reduced scattering coefficient ( $µ_{s}'$ ), calculated as the product of the scattering coefficient ($µ_{s}$) and $(1-g)$ is crucial in various fields, particularly in biomedical optics, because it allows for the non-invasive quantification of tissue properties. Users can use `µs' Power Law Fitting`  tab to compute the fitting parameters that provide a simplified functional form for the wavelength dependence of $µ_{s}'$ as described in `[@Jacques2013]`.



# Acknowledgments

We acknowledge support from the Laser Microbeam and Medical Program (LAMMP), a NIH Biomedical Technology Resource (P41-EB015890). JCR and VV acknowledge the support from the NIH (R21-GM128135) and the NSF (CBET-1805082).

# References
