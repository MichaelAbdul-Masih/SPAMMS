# SPAMMS 1.0.0

Introduction
------------
SPAMMS stands for Spectroscopic PAtch Model for Massive Stars and is designed with geometrically deformed systems in mind.  SPAMMS combines the eclipsing binary modelling code PHOEBE 2 and the NLTE radiative transfer code FASTWIND to produce synthetic spectra for systems at given phases, orientations and geometries.  For details on how the code works, please see the corresponding release paper: [Abdul-Masih et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020arXiv200309008A/abstract)

Installation
------------
*   SPAMMS is written in python 2.7 but in a future release, we will support python 3.7 as well.

*   Clone the git repository to create a local copy.

        $ git clone https://github.com/MichaelAbdul-Masih/SPAMMS.git
        
*   SPAMMS relies heavily on PHOEBE 2 so you will need to download PHOEBE 2 as well.  Instructions for installation can be found at http://phoebe-project.org/

*   Minimum package dependencies include:
        
        astropy 2.0.12
        numpy 1.16.3
        phoebe 2.1.15
        schwimmbad 0.3.0
        scipy 1.2.0
        tqdm 4.31.1
        

Getting Started
---------------
### Input file
SPAMMS works by refering to an input file (input.txt) which contains all of the settings.  There are separate input files for single star systems, detached binary systems and contact binary systems, as they contain different arguments.  These input files are broken up into 4 sections:

*   Object type: This defines the morphology of the system ('contact_binary', 'binary', or 'single').  For example:

        object_type = contact_binary

*   IO information: These are used to specify the paths for the input spectra (if applicable), FASTWIND grid and output paths.  If you do not wish to compare to an input spectrum, it can be set to 'None' and instead a times arguement can be passed with an array of times you wish to compute syntehtic spectra for.  For example:

        path_to_obs_spectra = None
        times = [0, 0.1, 0.2, 0.3]

*   System parameters: SPAMMS is built in such a way that these parameters can be given as single values or as arrays, in which case SPAMMS will compute a grid of models. The array of values can be given explicitly using square brackets or given in the form of a np.linspace arguement using parantheses.  All three possibilities are shown below:

        teff_primary =        44000                                     # single value
        teff_primary =        [42000, 43000, 44000, 45000,  46000]      # explicit array
        teff_primary =        (42000, 46000, 5)                         # using np.linspace arguement (returns same as explicit method above)

*   Selected line list: This specifies which lines you wish to compute.  A full line list for the LMC computed grid can be found in the settings.py script.

### FASTWIND Grid
Before using the code, an input grid will need to be either downloaded or created. A grid has been computed for the LMC and is available on request (michael.abdulmasih@gmail.com).  In a future release, we plan to make the input grid calculation scripts available.

### Running the code
SPAMMS can be run on a single core using a normal python call or on multiple cores using an mpiexec call.  By default, SPAMMS uses input.txt as the input file but this can be changed using the '-i' flag and specifying your chosen input file.  Additionally, to ensure that your radii stay within the boundaries of the computed LMC grid, a '-b' flag exists.  An example call can be found below:

        $ python spamms.py -i input.txt -b
        $ mpiexec -n 4 python spamms.py -i input.txt -b
