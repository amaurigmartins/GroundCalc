[![View GroundCalc on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/98739-groundcalc) [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=amaurigmartins/groundcalc ) 

# GroundCalc

Made for educational purposes, GroundCalc enables the user to carry out simplified (but not that much) grounding studies involving arbitrary grids and two-layered soils, using the microsegmentation method proposed by Dawalibi et al. It is possible to perform simulations using the good-ol' equipotential assumption, as well as accounting for the internal impedances of conductors.

The main technical highlights are:
- Intuitive, self-explaining UI, which requires a minimum of data entry.
- Ability to load conductor coordinates from CSV or TXT files, including AutoCAD data extractions and DXF files.
- Calculates grounding impedance, GPR, leakage currents, touch voltage and step voltage distributions.
- Calculates safe limits according to IEEE Std. 80.
- It is possible to specify the point where the fault is applied and the neutral return point. 
- Interactive plot control and grid preview.


### Important information

GroundCalc uses external tools that must be downloaded from the corresponding sources and placed in the appropriate directories, listed below:

- Sebastian (2023). Read DXF File Data (https://www.mathworks.com/matlabcentral/fileexchange/24572-read-dxf-file-data), MATLAB Central File Exchange. Retrieved October 26, 2023.

These codes are properties of the respective authors, with all due credits given. Observe any restrictions and licensing/usage requirements in the websites above.


### Acknowledgements

GroundCalc has been actively maintained and improved by [Amauri Martins-Britto](mailto:amaurigmartins@gmail.com) and [João Pedro Ivo Finamore](mailto:jpe.ivo@gmail.com).


### Disclaimer

This is not a professional power grounding simulation tool. The takeaway here is to put in perspective the design errors that may occur when overly simplified methods are adopted, in particular when conductor impedances and neutral return currents are neglected. The equations and methods implemented are fairly accurate for small to medium sized geometries at the mains frequency and that's all. DO NOT USE THIS FOR A REAL DESIGN. Go use CDEGS instead. Do not use spaghetti-copycat-generic-cheap software either. This is about human lifes, for God's sake.


### Restrictions of use

We appreciate the interest in our work and we invite the interested users to use our codes as necessary, as long as they are not embedded in any commercial software, which is **strictly prohibited**. If you use the GroundCalc as a part of scientific research, we kindly ask you to refer to one of our published papers:

- J. P. C. F. Ivo, C. M. Moraes and A. G. Martins-Britto, "Enhanced Circuit-Theory-Based Model for Accurate Simulation of Electrically Large Grounding Systems," 2023 Workshop on Communication Networks and Power Systems (WCNPS), Brasilia, Brazil, 2023, pp. 1-7, doi: 10.1109/WCNPS60622.2023.10344583.
- C. M. Moraes, G. H. S. Matos, A. G. Martins-Britto, K. L. M. Silva, F. V. Lopes, "Total AC Interferences Between a Power Line Subject to a Single-Phase Fault and a Nearby Pipeline With Multilayered Soil," in IEEE Transactions on Electromagnetic Compatibility, pp. 585-594, Vol. 65, Issue: 2, ISSN 0018-9375, Feb 2023, doi: 10.1109/TEMC.2023.3244095.
