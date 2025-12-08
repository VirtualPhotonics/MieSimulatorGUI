## Mie Simulator GUI
The [Mie Simulator GUI](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki) is an open-source, user-friendly tool designed to calculate the characteristics of Mie scatterers. It determines various properties, including the spectral dependence of the scattering coefficient, scattering cross-section, reduced scattering coefficient, phase function, average cosine of the phase function, and forward/backward scattering percentage. Additionally, it calculates scattering matrix entries for single or multiple wavelengths. Based on the Lorentz-Mie solution, these calculations can be performed for both mono-disperse and poly-disperse distributions. The tool can also estimate fitting parameters from the reduced scattering coefficient curve, though its capabilities are limited to the characteristics of the far-field scattering pattern. 

This GUI tool is based on the [BHMIE code](https://onlinelibrary.wiley.com/doi/book/10.1002/9783527618156). For absorbing spheres, the tool uses the complex refractive index definition (m<sub>sphere</sub> = m<sub>real</sub> – j m<sub>imag</sub>) as described in van de Hulst's [Light Scattering by Small Particles (1957) book](https://books.google.com/books?id=PlHfPMVAFRcC&source=gbs_navlinks_s). To understand the Mie algorithm and other parameters within the GUI, refer to, read Bohren and Huffman's [Absorption and Scattering of Light by Small Particles, (1983) book](https://onlinelibrary.wiley.com/doi/book/10.1002/9783527618156), Maetzler's [research report](http://www.atmo.arizona.edu/students/courselinks/spring09/atmo656b/maetzler_mie_v2.pdf) and Scott Prahl's ['Mie Basics'](https://miepython.readthedocs.io/en/latest/01_basics.html#) web page. Phase function calculations for mixtures of spheres are based on the work by Gélébart et al. in [Pure and Applied Optics, 5 (1996)](https://iopscience.iop.org/article/10.1088/0963-9659/5/4/005/meta). The reduced scattering power law fitting functions are derived from  Steve L. Jacques's review paper in [Phys. Med & Bio. 58 (2013)](https://iopscience.iop.org/article/10.1088/0031-9155/58/11/R37). All relevant documentation, including complete installation guides for Windows, Linux and macOS is available on the [Mie Simulator GUI Wiki](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki) page.

![Mie-Application](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki/images/Mie-Simulator-GUI-Screenshot.png)

## Downloads
To download the latest version, click [here](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki/Downloads)

## Releases and Source Codes
To access all releases and source codes, click [here](https://github.com/VirtualPhotonics/MieSimulatorGUI/releases). 

## Documentation

 * Getting started with the application and source code	
	
	- For Windows, click [here](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki/Getting-Started-on-Windows).
	- For Linux, click [here](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki/Getting-Started-on-Linux).
	- For MacOS, click [here](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki/Getting-Started-on-Mac).
 * Description of GUI panels and running application, click [here](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki/Mie-Simulator-GUI-Running-Application)
 
  * Guidelines for third parties
    - Contributing to the Software, click [here](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki/Coding-Conventions).
	- Reporting Issues, and Seeking Support, click [here](https://github.com/VirtualPhotonics/MieSimulatorGUI/discussions).
  * Comparison
    - For an example demonstrating the comparison between MieSimulatorGUI and Scott Prahl's [Mie Scattering Calculator](https://omlc.org/calc/mie_calc.html), click [here](https://github.com/VirtualPhotonics/MieSimulatorGUI/wiki/Examples).
 
## Acknowledgement
Use the following acknowledgement in publications or applications that make use of this open source software or underlying technology and research:

"__We acknowledge open-source software resources offered by the Virtual Photonics Technology Initiative, at the Beckman Laser Institute, University of California, Irvine.__"

## Software
_Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:_

_The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software._

_THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE._

## Copyright
Virtual Tissue Simulator (VTS) Project.

Copyright ©2025 [Virtual Photonics Technology Initiative](https://virtualphotonics.org/).
