
import numpy as np

from .regular_grid import RegularGrid

class Jra55Grid(RegularGrid):

    def __init__(self, description='JRA55 regular grid'):
        self.type = 'Arakawa A'
        super(Jra55Grid, self).__init__(360, 180, description=description)
