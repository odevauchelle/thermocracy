
#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# Olivier Devauchelle
#
# Dislike of general opinion makes for tight elections, O. Devauchelle, P. Szymczak, P. Nowakowski, Physical Review E, 109, 044106, 2024

"""
pyFreeFem
"""

__author__ = "Olivier Devauchelle"
__copyright__ = "Copyright 2020"
__license__ = "GPL"
__version__ = "0.0"

__all__ = ['Ising', 'graphics', 'lattice', 'statistics', 'Hamiltonian']

from .Ising import *
from .graphics import *
from .lattice import *
from .statistics import *
from .Hamiltonian import *
