"""
201511_approxchirp.py

Support material for the blog post "Approximating the area of a chirp by
fitting a polynomial", on DSPrelated.com.

* Author: Alexandre 'Jaguar' Fioravante de Siqueira
* Contact: http://www.programandociencia.com/sobre/
* Support material:
    http://www.github.com/alexandrejaguar/programandociencia

* In order to cite this material, please use the reference below
(this is a Chicago-like style):
de Siqueira, Alexandre Fioravante. “Approximating the area of a chirp by
fitting a polynomial”. DSPrelated. 2015, November 15. Available at
http://www.dsprelated.com/showarticle/863.php.
Access date: <ACCESS DATE>.

Copyright (C) Alexandre Fioravante de Siqueira

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# generating the x-axis and the chirp.
time_interval = np.linspace(0, 10, num=100, endpoint=True)
example_chirp = sp.signal.chirp(time_interval, 0, 1, 1, method='quadratic')

# plotting the discrete chirp.
plt.stem(time_interval, example_chirp)

# fitting three Chebyshev polynomimals, with degrees equal to 20, 45 and 115.
cheby_chirp20 = np.polynomial.Chebyshev.fit(time_interval,
                                            example_chirp,
                                            deg=20)

cheby_chirp45 = np.polynomial.Chebyshev.fit(time_interval,
                                            example_chirp,
                                            deg=45)

cheby_chirp115 = np.polynomial.Chebyshev.fit(time_interval,
                                             example_chirp,
                                             deg=115)

# plotting the Chebyshev polynomials.
f, ax_array = plt.subplots(3, sharex=True)
ax_array[0].stem(time_interval, example_chirp)
ax_array[0].plot(time_interval, cheby_chirp20(time_interval), 'g', lw=2, label='cheby deg 20')
ax_array[0].legend()

ax_array[1].stem(time_interval, example_chirp)
ax_array[1].plot(time_interval, cheby_chirp45(time_interval), 'g', lw=2, label='cheby deg 45')
ax_array[1].legend()

ax_array[2].stem(time_interval, example_chirp)
ax_array[2].plot(time_interval, cheby_chirp115(time_interval), 'g', lw=2, label='cheby deg 115')
ax_array[2].legend()

# defining the integration interval and calculating the areas.

integ_int = [0, np.pi]
cheby20_area, cheby20_err = sp.integrate.quad(cheby_chirp20, *integ_int)
cheby45_area, cheby45_err = sp.integrate.quad(cheby_chirp45, *integ_int)
cheby115_area, cheby115_err = sp.integrate.quad(cheby_chirp115, *integ_int)

# nicely printing the results.
print('The area of cheby20 in this interval is {0}. The absolute error is {1}'.format(cheby20_area, cheby20_err))
print('The area of cheby45 in this interval is {0}. The absolute error is {1}'.format(cheby45_area, cheby45_err))
print('The area of cheby115 in this interval is {0}. The absolute error is {1}'.format(cheby115_area, cheby115_err))

# checking the integration areas.
integ_full = np.linspace(*integ_int)

f, ax_array = plt.subplots(3, sharex=True)
ax_array[0].stem(time_interval, example_chirp)
ax_array[0].plot(time_interval, cheby_chirp20(time_interval), 'g', lw=2, label='cheby deg 20')
ax_array[0].fill_between(integ_full, cheby_chirp20(integ_full), color='orchid', label=r'$\int_{0}^{\pi}\,$cheby_chirp20(x)')
ax_array[0].legend()

ax_array[1].stem(time_interval, example_chirp)
ax_array[1].plot(time_interval, cheby_chirp45(time_interval), 'g', lw=2, label='cheby deg 45')
ax_array[1].fill_between(integ_full, cheby_chirp45(integ_full), color='orchid', label=r'$\int_{0}^{\pi}\,$cheby_chirp45(x)')
ax_array[1].legend()

ax_array[2].stem(time_interval, example_chirp)
ax_array[2].plot(time_interval, cheby_chirp115(time_interval), 'g', lw=2, label='cheby deg 115')
ax_array[2].fill_between(integ_full, cheby_chirp115(integ_full), color='orchid', label=r'$\int_{0}^{\pi}\,$cheby_chirp115(x)')
ax_array[2].legend()
