---
title: >-
  MieSimulatorGUI: A user-friendly tool to compute and visualize light scattering by spherical dielectric particles
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
  - name: Carole K. Hayakawa
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
date: 11 December 2025
bibliography: paper.bib

---

[comment]: https://joss.theoj.org/]

# Summary

Mie theory is a mathematical framework derived from Maxwell's equations that models electromagnetic scattering by spherical dielectric particles.  Its predictions are essential across various scientific and engineering disciplines, including biomedical optics, atmospheric optics, particle characterization, nanofluids, computer graphics, and remote sensing. We developed `MieSimulatorGUI` to bridge the gap for researchers who require Mie simulations but lack specialized programming expertise. By integrating a high-performance C/C++ computational engine with the Qt framework, this user-friendly cross-platform tool calculates key optical properties such as scattering coefficients, cross-sections, angular scattering distributions, and scattering asymmetry for both monodisperse and polydisperse particle distributions. The graphical user interface (GUI) features six interactive panels that allow users to specify input parameters, visualize distributions, and analyze scattering phenomena. Furthermore, it enables the fitting of spectrally varying reduced scattering coefficient parameters, a feature particularly valuable in fields like tissue optics.
`MieSimulatorGUI` is open source, hosted on [GitHub](https://github.com/VirtualPhotonics/MieSimulatorGUI), and distributed under [the MIT license](https://opensource.org/license/mit). It is available via the [project's download page](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki/Downloads).

# Statement of need

Mie theory is derived from Maxwell's equations and provides a comprehensive framework for modeling electromagnetic scattering by spherical particles [@Mie1908]. Mie theory is utilized across diverse fields, ranging from nanomaterials and biomedical optics to atmospheric science and astronomy [@Goody1989; @Saidi1995; @Wang2005; @Chalut2008; @Horvath2009; @Bhandari11]. Despite its broad applicability, the theory's reliance on complex mathematical constructs, such as infinite series and special functions [@VandeHulst1957; @Bohren1983; @Wiscombe1980; @Majic2020], demands advanced computational implementation. 

While numerous Mie simulation packages are available (many of which are listed on [SCATTPORT.org](https://scattport.org) and [Wikipedia](https://en.wikipedia.org/wiki/Codes_for_electromagnetic_scattering_by_spheres)), they generally fall into two categories: older, established codes focusing on computational efficiency [@Wiscombe1980; @Bohren1983], and newer, object-oriented libraries typically hosted on version-control platforms [@Sumlin2018; @PoinsinetdeSivry-Houle2023; @Prahl_mie; @MieScattering]. Although both categories provide robust computational engines, they usually demand significant programming proficiency. This requirement creates a barrier for experimentalists, clinical scientists, and educators who need these analytical capabilities but may lack the specialized coding expertise to integrate such libraries into their workflows.

`MieSimulatorGUI` bridges this gap by providing an intuitive, cross-platform desktop application that computes and fits scattering parameters for monodisperse and polydisperse distributions without requiring any coding. Unlike standard implementations, it supports heterogeneous polydispersity, allowing users to assign bin-specific complex refractive indices via custom data inputs, a feature often absent in simplified GUI tools. The tool facilitates high-impact use cases such as biomedical optics [@Mourant1997; @Wang2005; @Jacques2013] and atmospheric research [@Seinfeld1998; @Teri2022], where users can define complex configurations and directly fit spectrally varying reduced scattering coefficients. By integrating a powerful C/C++ computational engine with intuitive [Qt](https://www.qt.io/) interface, `MieSimulatorGUI` offers accessible yet powerful Mie theory computations, facilitating both streamlined research analysis and interactive pedagogical demonstrations. 

# Main Features

Built on the BHMIE [@Bohren1983] and Wiscombe [@Wiscombe1979] frameworks, `MieSimulatorGUI` calculates spectral optical properties, including scattering coefficients, cross-sections, scattering amplitude matrix entries ($\text{S}_1$, $\text{S}_2$), phase functions, and scattering asymmetry. The tool supports both monodisperse and polydisperse particle distributions [@Gelebart1996] and facilitates parameter estimation by fitting reduced scattering coefficient curves, a technique of significant value in tissue optics [@Jacques2013]. To ensure computational accuracy and GUI stability, the software includes an automated test suite integrated via GitHub Actions. Its high-performance C++ engine enables near-instantaneous computation and plotting, providing real-time visual updates across the spectral range. Comprehensive documentation, including cross-platform installation guides (Windows, Linux, and macOS), dependency specifications, command-line test execution, dependent scattering warning trigger conditions, and several examples, is available on the [Mie Simulator GUI Wiki](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki) page. 

![Figure 1: Six interactive panels of `MieSimulatorGUI`: (a) Input selection, (b) Particle size distribution, (c) Scattering coefficient, (d) Phase function, (e) Reduced Scattering,  and (f) Scattering Asymmetry\label{fig:Figure1}](Figure1.png){ width=100% }
Figure 1: Six interactive panels of `MieSimulatorGUI`: (a) Input selection, (b) Particle size distribution, (c) Scattering coefficient, (d) Phase function, (e) Reduced Scattering, and (f) Scattering Asymmetry (Anisotropy)


# Design and Functionality

The tool is distributed as portable binaries for Windows, macOS, and Linux. For local compilation, the project utilizes `qmake`, with dependencies (`Qt6` and `QCustomPlot`) managed via an automated build script. The Qt GUI contains six interactive panels (Figure 1). 


## Input Selection Panel

This panel enables the user to define inputs for Mie simulations at either a single wavelength or across a spectral range. The distribution of spheres is described using sphere concentration (`Conc`) ($\text{spheres/mm}^3$) or volume fraction (`Vol Frac`). The volume fraction represents the ratio of total sphere volume to the total container volume. For polydisperse systems, it is calculated by summing the products of each sphere size's volume and its respective concentration per unit volume. `MieSimulatorGUI` utilizes the independent scattering approximation, a framework valid for dilute suspensions where particles are sufficiently separated to ignore coherent interactions [@VandeHulst1957; @Schmitt1998]. Since the validity of this approximation depends on the inter-particle spacing, size parameter, and relative refractive index [@Tien1987; @Ivezic1996; @Galy2020; @Yalcin2022], its accuracy diminishes as the volume fraction increases. Consequently, the tool is best suited for dilute systems and the results obtained for concentrated regimes may deviate from physical reality and should be interpreted with caution. `MieSimulatorGUI` triggers a warning when parameters violate dependent scattering conditions.

To maintain numerical stability in the BHMIE algorithm and ensure UI responsiveness, sphere diameters are restricted to a range of 0.1 $\text{nm}$ to 300 $\text{µm}$ and wavelengths are limited to $50\text{ nm}$ – $3000\text{ nm}$. These ranges cover the primary biomedical and atmospheric spectral windows. For absorbing spheres, the complex refractive index ($m_{sphere}$) is defined as $m_{real}$ – j $m_{imag}$, where $m_{real}$ and $m_{imag}$ represent real and imaginary components, respectively [@VandeHulst1957; @Wiscombe1979]. 

The tool provides options for either monodisperse (uniform-sized) or polydisperse (variable-sized) particle distributions. Monodisperse distributions restrict analysis to spheres with uniform size and refractive index. In contrast, polydisperse distributions enable simulations of spheres with diverse attributes and support three size distribution models: 1. Log-normal, 2. Gaussian, and 3. user-defined (custom). The user-defined option allows for the specification of different refractive indices for different spheres, as demonstrated in the examples provided in the `CustomDataSamples` folder. 

## Number Density Panel

This panel graphically presents the number density of spheres $N_s$ [\#/ $\text{mm}^3$] used in the simulation.  The subsequent tab displays the `Size Parameter` defined as $\text{2π}Rn_{med} / λ_{vacuum}$ [@Bohren1983], where $R$ [ $\text{µm}$ ] denotes the particle radius, $n_{med}$  the medium's refractive index, and $λ_{vacuum}$ [ $\text{µm}$ ] is the wavelength in vacuum. 

## Scattering Coefficient Panel

The Mie calculations provide three important efficiency factors: the scattering efficiency ($Q_{sca}$), the extinction efficiency ($Q_{ext}$), and the backscattering efficiency ($Q_{back}$) (see [Mie Scattering Efficiencies](https://miepython.readthedocs.io/en/latest/02_efficiencies.html)). These dimensionless quantities combined with the particle's cross sectional area ($\text{π}R^2$) yield the corresponding scattering cross section $C_{sca}$ [ $\text{/mm}^2$ ], extinction cross section $C_{ext}$ [ $\text{/mm}^2$ ] and backscattering cross section $C_{back}$ [ $\text{/mm}^2$ ] [@VandeHulst1957; @Bohren1983]. The calculated cross sections are displayed across three separate tabs. 
For monodisperse distribution, the [scattering coefficient](https://omlc.org/classroom/ece532/class3/musdefinition.html) ($\text{µ}_s$) is simply the product of the scattering cross-section $C_{sca}$ and the number density $N_s$. For polydisperse distributions, the scattering coefficient is computed via a discrete summation of the cross-sections across individual particle size bins, as detailed by [@Schmitt1998]. 

## Phase Function Panel

The [phase function](https://omlc.org/classroom/ece532/class3/ptheta.html) represents the angular distribution of scattered light.  The [calculated phase function](https://miepython.readthedocs.io/en/latest/03_angular_scattering.html) results are displayed in this panel using both polar and linear plots. These plots are derived from the [complex amplitude scattering matrix elements](https://omlc.org/classroom/ece532/class3/mie_math.html), $\text{S}_1$, and $\text{S}_2$, which describe the transformation of incident electromagnetic field to far-field scattered field [@VandeHulst1957; @Bohren1983]. The wavelength slider allows the user to visualize the phase function or $\text{S}_1$ and $\text{S}_2$ data at any specific wavelength.

## Scattering Asymmetry Panel

The scattering asymmetry ([Anisotropy](https://omlc.org/classroom/ece532/class3/gdefinition.html)) panel displays the directional properties of the scattering phase function. The first tab presents the average cosine of the single scattering phase function ( $g$ ), which quantifies the prevalence of forward ($g>0$) vs backward scattering ($g<0$). The second tab provides the integrated forward and backward scattering fractions, offering a detailed analysis of the angular scattering distribution.

## Reduced Scattering Panel

This panel shows the [reduced scattering coefficient](https://omlc.org/classroom/ece532/class3/musp.html) ( $\text{µ}_s'$ ), which is computed as the product of the scattering coefficient ($\text{µ}_s$) and $(1-g)$ [@Jacques2013]. This parameter of particular interest in biomedical optics, because it allows for the non-invasive quantification of tissue properties. Users can use `$\text{µ}_s'$ Power Law Fitting` tab to compute the fitting parameters that provide a simplified functional form for the wavelength dependence of $\text{µ}_s'$ as described in [@Jacques2013].

## Example Application: Scattering of Intralipid Phantoms

To demonstrate the tool's scientific utility, consider the characterization of Intralipid 20%, a common tissue phantom in biomedical optics [@DiNinni2011]. Users can define the medium by inputting polydisperse particle distribution parameters using a `Log Normal` distribution. For Intralipid 20% (w/w), which corresponds to a volume fraction (`Vol Frac`) of 0.227 [@Aernouts2013], users may assume a mean particle diameter of 0.22 $\text{µm}$ and a standard deviation of 0.36 $\text{µm}$, following Raju (2017) [@Raju2017]. After setting the `Num. sph. sizes` field to 101, users can specify the refractive index of the soybean oil droplets as 1.47 and the surrounding medium as 1.33. Upon executing the simulation across the 600–1000 $\text{nm}$ spectral range, `MieSimulatorGUI` calculates $\text{µ}_s$, $\text{µ}_s'$ and g. Although these inputs may deviate from the independent scattering regime, these results are comparable to the bulk scattering properties shown in Figure 14 of Aernouts et al. [@Aernouts2013]. Consequently, the generated data can be exported as text files for further analysis.

# Acknowledgments

We acknowledge support from the Laser Microbeam and Medical Program (LAMMP), a NIH Biomedical Technology Resource (P41-EB015890). JCR and VV acknowledge the support from the NIH (R21-GM128135) and the NSF (CBET-1805082). JCR was the primary contributor, handling the core software development, UI design, validation, GitHub upload, and drafting the original manuscript. VV provided overall supervision and critically reviewed and edited the final manuscript. CKH developed the initial Mie scattering code, and LG assisted with the GitHub upload.

# References
