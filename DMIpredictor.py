import DMIDB

# Have to define the unit of DMI prediction, either a pair of proteins or a file that contains many protein pairs
# The function should work with both options, and the arguments given to the function should deal with this.
# can set up options and defaults to deal with this too
class DMIpredictor(DMIDB.InterfaceHandling):
    def __init__(self):
        super().__init__()

# Set up another script to do ML training and learning and outputs a model
# read in the ML model
# Most ML functions like scoring DMI matches will come in here...
