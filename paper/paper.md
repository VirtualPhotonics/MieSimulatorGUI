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
  - name: Lisa M. Glover
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
date: 4 February 2026
bibliography: paper.bib

---

[comment]: https://joss.theoj.org/]

# Summary

Mie theory is a mathematical framework derived from Maxwell's equations that models electromagnetic scattering by spherical dielectric particles.  Its predictions are essential across various scientific and engineering disciplines, including biomedical optics, atmospheric optics, particle characterization, nanofluids, computer graphics, and remote sensing. We developed `MieSimulatorGUI` to bridge the gap for researchers who require Mie simulations but lack specialized programming expertise. By integrating a high-performance C/C++ computational engine with the Qt framework, this user-friendly cross-platform tool calculates key optical properties such as scattering coefficients, cross-sections, angular scattering distributions, and scattering asymmetry for both monodisperse and polydisperse particle distributions. The graphical user interface (GUI) features six interactive panels that allow users to specify optical and particle input parameters, visualize particle distributions, and compute scattering metrics. Furthermore, it enables the fitting of the spectral dependence of the reduced scattering coefficient, a feature particularly valuable in fields like tissue optics.
`MieSimulatorGUI` is open source, hosted on [GitHub](https://github.com/VirtualPhotonics/MieSimulatorGUI), and distributed under [the MIT license](https://opensource.org/license/mit). It is available via the [project's download page](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki/Downloads).

# Statement of need

Mie theory is derived from Maxwell's equations and provides a comprehensive framework for modeling electromagnetic scattering by spherical particles [@Mie1908]. Mie theory is utilized across diverse fields, ranging from nanomaterials and biomedical optics to atmospheric science and astronomy [@Goody1989; @Saidi1995; @Wang2005; @Chalut2008; @Horvath2009; @Bhandari11]. Despite its broad applicability, the theory's reliance on complex mathematical constructs, such as infinite series and special functions [@VandeHulst1957; @Bohren1983; @Wiscombe1980; @Majic2020], demands advanced computational implementation. 

While numerous Mie simulation packages are available (many of which are listed on [SCATTPORT.org](https://scattport.org) and [Wikipedia](https://en.wikipedia.org/wiki/Codes_for_electromagnetic_scattering_by_spheres)), they generally fall into two categories: older, established codes focusing on computational efficiency [@Wiscombe1980; @Bohren1983], and newer, object-oriented libraries typically hosted on version-control platforms [@Sumlin2018; @PoinsinetdeSivry-Houle2023; @Prahl_mie; @MieScattering]. Although both categories provide robust computational engines, they usually demand significant programming proficiency. This requirement creates a barrier for experimentalists, clinical scientists, and educators who need these analytical capabilities but may lack the specialized coding expertise to integrate such libraries into their workflows.

`MieSimulatorGUI` bridges this gap by providing an intuitive, cross-platform desktop application that computes and fits scattering parameters for monodisperse and polydisperse distributions without any coding. Unlike standard implementations, it supports heterogeneous polydispersity, allowing users to assign bin-specific complex refractive indices via custom data inputs, a feature often absent in simplified GUI tools. The tool facilitates high-impact use cases such as biomedical optics [@Mourant1997; @Wang2005; @Jacques2013] and atmospheric research [@Seinfeld1998; @Teri2022], where users can define complex particle configurations and directly fit spectrally-varying reduced scattering coefficients. By integrating a powerful C/C++ computational engine with intuitive [Qt](https://www.qt.io/) interface, `MieSimulatorGUI` offers accessible, yet powerful Mie theory computations, facilitating both streamlined research analysis and interactive pedagogical demonstrations. 

# Main Features

Built on the BHMIE [@Bohren1983] and Wiscombe [@Wiscombe1979] frameworks, `MieSimulatorGUI` calculates spectral optical properties, including scattering coefficients, cross-sections, scattering amplitude matrix entries ($S_1$, $S_2$), phase functions, and scattering asymmetry. The tool supports both monodisperse and polydisperse particle distributions [@Gelebart1996] and facilitates parameter estimation by fitting reduced scattering coefficient curves, a technique of significant value in tissue optics [@Jacques2013]. To ensure computational accuracy and GUI stability, the software includes an automated test suite integrated via GitHub Actions. Its high-performance C++ engine enables near-instantaneous computation and plotting, providing real-time visual updates across the spectral range. Comprehensive documentation, including cross-platform installation guides (Windows, Linux, and macOS), dependency specifications, command-line test execution, dependent scattering warning trigger conditions, and several examples, is available on the [Mie Simulator GUI Wiki](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki) page. 

![Figure 1: Six interactive panels of `MieSimulatorGUI`: (a) Input selection, (b) Particle size distribution, (c) Scattering coefficient, (d) Phase function, (e) Reduced Scattering,  and (f) Scattering Asymmetry\label{fig:Figure1}](Figure1.png){ width=100% }
Figure 1: Six interactive panels of `MieSimulatorGUI`: (a) Input selection, (b) Particle size distribution, (c) Scattering coefficient, (d) Phase function, (e) Reduced Scattering, and (f) Scattering Asymmetry (Anisotropy)


# Design and Functionality

The tool is distributed as portable binaries for Windows, macOS, and Linux. For local compilation, the project utilizes `qmake`, with dependencies (`Qt6` and `QCustomPlot`) managed via an automated build script. The Qt GUI contains six interactive panels (Figure 1). 


## Input Selection Panel

This panel enables the user to define inputs for Mie simulations at either a single wavelength or across a spectral range. The distribution of spheres is described using sphere concentration (`Conc`) ($\text{spheres/mm}^3$) or volume fraction (`Vol Frac`). The volume fraction represents the ratio of the volume occupied by the spherical particles to the total solution volume. For polydisperse systems, `Vol Frac` is calculated by multiplying the volume of each sphere size with its corresponding concentration per unit volume, and then summing across all sphere sizes. `MieSimulatorGUI` utilizes the independent scattering approximation, a framework valid for dilute suspensions where particles are sufficiently separated to ignore coherent interactions [@VandeHulst1957; @Schmitt1998]. The accuracy of this approximation decreases as volume fraction increases, as it is sensitive to both inter-particle spacing and size parameters [@Tien1987; @Galy2020; @Yalcin2022]. Consequently, the tool is best suited for dilute systems and the results obtained for concentrated regimes may deviate from physical reality and should be interpreted with caution. `MieSimulatorGUI v2.0` triggers a warning if the inputs exceed the limits of independent scattering.

To maintain numerical stability in the BHMIE algorithm and ensure UI responsiveness, sphere diameters are restricted to a range of 0.1 $\text{nm}$ to 300 $\mu\text{m}$ and wavelengths are limited to 50 $\text{nm}$ – 3000 $\text{ nm}$. These ranges cover the primary biomedical and atmospheric spectral windows. For absorbing spheres, the complex refractive index ($m_{sphere}$) is defined as $m_{real}$ – j $m_{imag}$, where $m_{real}$ and $m_{imag}$ represent real and imaginary components, respectively [@VandeHulst1957; @Wiscombe1979]. 

The tool provides options for either monodisperse (uniform-sized) or polydisperse (variable-sized) particle distributions. Monodisperse distributions restrict analysis to spheres with uniform size and refractive index. In contrast, polydisperse distributions enable simulations of spheres with diverse attributes and support three size distribution models: 1. Log-normal, 2. Gaussian, and 3. user-defined (custom). The user-defined option allows for the specification of different refractive indices for different spheres, as demonstrated in the examples provided in the `CustomDataSamples` folder. 

## Number Density Panel

This panel graphically presents the number density of spheres $N_s$ [\#/ $\text{mm}^3$] used in the simulation.  The subsequent tab displays the `Size Parameter` defined as 2$\pi Rn_{med} / \lambda_{vacuum}$ [@Bohren1983], where $R$ [ $\mu\text{m}$ ] denotes the particle radius, $n_{med}$  the medium's refractive index, and $\lambda_{vacuum}$ [ $\mu\text{m}$ ] is the wavelength in vacuum. 

## Scattering Coefficient Panel

The Mie calculations provide three important efficiency factors: the scattering efficiency ($Q_{sca}$), the extinction efficiency ($Q_{ext}$), and the backscattering efficiency ($Q_{back}$) (see [Mie Scattering Efficiencies](https://miepython.readthedocs.io/en/latest/02_efficiencies.html)). These dimensionless quantities combined with the particle's cross sectional area ($\pi R^2$) yield the corresponding scattering cross section $C_{sca}$ [ $\text{/mm}^2$ ], extinction cross section $C_{ext}$ [ $\text{/mm}^2$ ] and backscattering cross section $C_{back}$ [ $\text{/mm}^2$ ] [@VandeHulst1957; @Bohren1983]. The calculated cross sections are displayed across three separate tabs. 
For monodisperse distribution, the [scattering coefficient](https://omlc.org/classroom/ece532/class3/musdefinition.html) ($\mu_{s}$) is simply the product of the scattering cross-section $C_{sca}$ and the number density $N_{s}$. For polydisperse distributions, the scattering coefficient is computed via a discrete summation of the cross-sections across individual particle size bins, as detailed by [@Schmitt1998]. 

## Phase Function Panel

The [phase function](https://omlc.org/classroom/ece532/class3/ptheta.html) represents the angular distribution of scattered light.  The [calculated phase function](https://miepython.readthedocs.io/en/latest/03_angular_scattering.html) results are displayed in this panel using both polar and linear plots. These plots are derived from the [complex amplitude scattering matrix elements](https://omlc.org/classroom/ece532/class3/mie_math.html), $S_1$, and $S_2$, which describe the transformation of incident electromagnetic field to far-field scattered field [@VandeHulst1957; @Bohren1983]. The wavelength slider allows the user to visualize the phase function or $S_1$ and $S_2$ data at any specific wavelength.

## Scattering Asymmetry Panel

The scattering asymmetry ([Anisotropy](https://omlc.org/classroom/ece532/class3/gdefinition.html)) panel displays the directional properties of the scattering phase function. The first tab presents the average cosine of the single scattering phase function ( $g$ ), which quantifies the prevalence of forward ($g>$0) vs backward scattering ($g<$0). The second tab provides the integrated forward and backward scattering fractions, offering a detailed analysis of the angular scattering distribution.

## Reduced Scattering Panel

This panel shows the [reduced scattering coefficient](https://omlc.org/classroom/ece532/class3/musp.html) ( $\mu_{s}'$ ), which is computed as the product of the scattering coefficient ($\mu_{s}$) and $(1-g)$ [@Jacques2013]. This parameter of particular interest in biomedical optics, because it allows for the non-invasive quantification of tissue properties. Users can use $\mu_{s}'$ `Power Law Fitting` tab to compute the fitting parameters that provide a simplified functional form for the wavelength dependence of $\mu_{s}'$ as described in [@Jacques2013].

## Example Application: Scattering of Intralipid Phantoms

To demonstrate the tool's scientific utility, we considered the characterization of Intralipid, a standard tissue phantom in biomedical optics [@vanStaveren1991; @DiNinni2011]. Based on Intralipid particle distribution profiles in the literature [@Kodach2011; @Raju2017], we assumed a polydisperse `Log Normal` particle distribution with a mean diameter of 0.22 $\mu\text{m}$ and a standard deviation of 0.37 $\mu\text{m}$. We set the refractive indices to 1.47 for the soybean oil droplets and 1.33 for the surrounding medium, while assigning a value of 101 to the `Num. sph. sizes` field. To analyze different concentrations ranging from 0.2% to 20% [@Aernouts2013; @vanStaveren1991], volume fractions were scaled using a baseline value of 0.227 for a 20% (w/w) Intralipid concentration [@Aernouts2013]. Upon executing the simulation across the 400–2250 $\text{nm}$ spectral range, `MieSimulatorGUI` calculates $\mu_{s}$, $\mu_{s}'$ and $g$. While the selected volume fractions may exceed independent scattering limits established in the literature [@Tien1987; @Galy2020; @Yalcin2022], the results show strong agreement with established bulk scattering properties [@vanStaveren1991; @Aernouts2013]. Figures can be exported as .png files and results as text files for further analysis. Beyond preset distributions, the `Custom` option can be utilized to upload specific sphere profiles [@Raju2017].

# Acknowledgments

We acknowledge support from the Laser Microbeam and Medical Program (LAMMP), a NIH Biomedical Technology Resource (P41-EB015890). JCR and VV acknowledge the support from the NIH (R21-GM128135) and the NSF (CBET-1805082). JCR was the primary contributor, handling the core software development, UI design, validation, GitHub upload, and drafting the original manuscript. VV provided overall supervision and critically reviewed and edited the final manuscript. CKH developed the initial Mie scattering code, and LMG assisted with the GitHub and Zenodo uploads.

# References
