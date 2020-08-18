## Mie Simulator GUI
The [Mie Simulator GUI](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki) is an open source tool for calculating the characteristics of Mie scatterers. This user-friendly tool calculates the spectral dependence of the scattering coefficient, scattering cross section, reduced scattering coefficient, phase function, average cosine of the phase function, forward / backward scattering percentage and scattering matrix entries for a single or series of wavelengths. These calculations can be run for mono-disperse and poly-disperse distributions based on Lorentz-Mie solution. The tool is also capable of estimating fitting parameters from the reduced scattering coefficient curve and is limited to the characteristics of the far-field scattering pattern. Mie scattering algorithm in the simulator is based on the [BHMIE code](https://onlinelibrary.wiley.com/doi/book/10.1002/9783527618156). For absorbing spheres, we considered complex refractive index definition (m<sub>sphere</sub> = m<sub>real</sub> – j m<sub>imag</sub>) in van de Hulst's [Light Scattering by Small Particles (1957) book](https://books.google.com/books?id=PlHfPMVAFRcC&source=gbs_navlinks_s). To learn Mie algorithm and other parameters in Mie Simulator GUI, read Bohren and Huffman's [Absorption and Scattering of Light by Small Particles, (1983) book](https://onlinelibrary.wiley.com/doi/book/10.1002/9783527618156), Maetzler's [research report](http://www.atmo.arizona.edu/students/courselinks/spring09/atmo656b/maetzler_mie_v2.pdf) and Scott Prahl's ['Mie Basics' web pages](https://miepython.readthedocs.io/en/latest/01_basics.html#). The phase function calculation for a mixture of spheres is based on Gélébart et al's [Pure and Applied Optics, 5 (1996) paper](https://iopscience.iop.org/article/10.1088/0963-9659/5/4/005/meta), and the reduced scattering power law fiiting functions are based on Steve L. Jacques's [Phys. Med & Bio. 58 (2013) review paper](https://iopscience.iop.org/article/10.1088/0031-9155/58/11/R37).

![Mie-Application](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki/images/Mie-Simulator-GUI-Screenshot.png)

## Downloads
To download the latest version, click [here](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki/Downloads)

## Releases and Source Codes
To access all releases and source codes, click [here](https://github.com/VirtualPhotonics/MieSimulatorGUI/releases). 

## Documentation
 * Getting started in Windows, click [here](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki/Getting-Started-in-Windows)
 * Getting started in Linux, click [here](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki/Getting-Started-in-Linux)
 * To access breif description of inputs, outputs, and execution of the tool, click [here](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki/Mie-Simulator-GUI-Running-Application)

## Acknowledgement
Use the following acknowledgement in publications or applications that make use of this open source software or underlying technology and research:

"__We acknowledge open-source software resources offered by the Virtual Photonics Technology Initiative, at the Beckman Laser Institute, University of California, Irvine.__"

## Software
_Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:_

_The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software._

_THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE._

## Copyright
Virtual Tissue Simulator (VTS) Project.

Copyright ©2020 [Virtual Photonics Technology Initiative](https://virtualphotonics.org/).
