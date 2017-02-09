
from .regular_grid import RegularGrid

class Core2Grid(RegularGrid):

    def __init__(self, description='CORE2 regular grid'):
        self.type = 'Arakawa A'
        super(Core2Grid, self).__init__(192, 94, description=description)
