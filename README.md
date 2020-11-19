# Phonon-DOS-based-on-MD-trajectory
Normal mode density from molecular dynamics velocity trajectory
The DOS is obtained by summing up the Fourier transform of the atomic velocity.
Users will have to provide a velocity trajectory file (I use CP2K, so the default velocity unit is bohr/au_time), and a mass file that contains the atomic mass.
The trajectory file is in .xyz format, users can check for that at cp2k.org
The mass file will be something like:
1.008
1.008
1.008
1.008
12.011
This is a CH4 mass file.
Users will need to input four parameters in the dos.py and thermo.py, they are:
atom_num (atom number of the system)
frame_num (number of MD steps)
t_total (simulation time in fs)
T (temperature in Kelvin)

Usage: python dos.py your_trajectory
A DoS file will be generated, you can plot it with gnuplot
Usage: python thermo.py
Vibrational energy will be printed
