
from .regular_grid import RegularGrid

class Core2Grid(RegularGrid):

    def __init__(self, description='CORE2 regular grid'):
        self.type = 'Arakawa A'
        self.full_name = 'CORE2'
        super(Core2Grid, self).__init__(192, 94, description=description)
