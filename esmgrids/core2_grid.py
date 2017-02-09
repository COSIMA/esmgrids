
from .regular_grid import RegularGrid

class Core2Grid(RegularGrid):

    def __init__(self, mask=None,
                 description='COREII regular grid'):
        super(Core2Grid, self).__init__(192, 94, mask=mask,
                                        levels=levels, description=description)
