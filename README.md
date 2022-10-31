# SPAMMS 1.0.dev

Introduction
------------
SPAMMS stands for Spectroscopic PAtch Model for Massive Stars and is designed with geometrically deformed systems in mind.  SPAMMS combines the eclipsing binary modelling code PHOEBE 2 and the NLTE radiative transfer code FASTWIND to produce synthetic spectra for systems at given phases, orientations and geometries.  For details on how the code works, please see the corresponding release paper: [Abdul-Masih et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020arXiv200309008A/abstract)

Installation
------------
*   Clone the git repository to create a local copy.

        $ git clone https://github.com/MichaelAbdul-Masih/SPAMMS.git

*   SPAMMS is written in Python 3 and has several dependencies that are needed to make it run.  For your convenience, we have provided an environment file (SPAMMS_environment.yml), which you can use to create an environment that already has all of these dependencies.  If you are using Anaconda you can create the environment as follows:

        conda env create -f SPAMMS_environment.yml

*   If you prefer to create the environment yourself, the minimum package requirements can be found below:

        astropy 2.0.12
        numpy 1.16.3
        phoebe 2.1.15
        schwimmbad 0.3.0
        scipy 1.2.0
        tqdm 4.31.1
        phoebe 2.3.0

Getting Started
---------------
### Test the Installation
To make sure the installation was successful, cd into the SPAMMS git directory and run the following command:

        python spamms.py -i input.txt

If the final message you get from the run is 'All models ran successfully!' then SPAMMS has been installed successfully.

### Input file
SPAMMS works by referring to an input file (input.txt) which contains all of the settings.  There are separate input files for single star systems, detached binary systems and contact binary systems, as they contain different arguments.  These input files are broken up into 5 sections:

*   Object type: This defines the morphology of the system ('contact_binary', 'binary', or 'single').  For example:

        object_type = contact_binary

*   IO information: These are used to specify the paths for the input spectra (if applicable), FASTWIND grid and output paths.  If you do not wish to compare to an input spectrum, it can be set to 'None' and instead a times argument can be passed with an array of times you wish to compute synthetic spectra for.  For example:

        path_to_obs_spectra = None
        times = [0, 0.1, 0.2, 0.3]

*   System parameters: SPAMMS is built in such a way that these parameters can be given as single values or as arrays, in which case SPAMMS will compute a grid of models. The array of values can be given explicitly using square brackets or given in the form of a np.linspace argument using parentheses.  All three possibilities are shown below:

        teff_primary =        44000                                     # single value
        teff_primary =        [42000, 43000, 44000, 45000,  46000]      # explicit array
        teff_primary =        (42000, 46000, 5)                         # using np.linspace argument (returns same as explicit method above)

*   Grid information: These values specify the setup of the grid.  The way in which the calculations are carried out allows for multiple abundance combinations to be computed simultaneously, but doing so requires the dimensions of the grid to be known.  The first two arguments are the helium and the CNO abundances, and the third is the number of wavelength bins in each line profile (by default FASTWIND uses 161).  An example is shown below:

        he_abundances =       [0.06, 0.1, 0.15, 0.2]
        cno_abundances =      [6.5, 7.0, 7.5, 8.0, 8.5]
        lp_bins =             161


*   Selected line list: This specifies which lines you wish to compute.  A full line list for the LMC computed grid can be found in the settings.py script.

### FASTWIND Grid
Before using the code, an input grid will need to be either downloaded or created. A grid has been computed for the LMC and is available on request (michael.abdulmasih@gmail.com).  In a future release, we plan to make the input grid calculation scripts available.

### Running the code
SPAMMS can be run on a single core or on multiple cores by passing the '-n' flag followed by the number of cores you wish to use.  SPAMMS is parallelized on a model basis meaning that when running a grid of models with different input parameters, each can be run on a different core.  Currently, it is not possible to parallelize based on the time points for a given model.  By default, SPAMMS uses input.txt as the input file but this can be changed using the '-i' flag and specifying your chosen input file. Before running SPAMMS, it is advised to run a grid check to make sure that your model does not fall out of the bounds of the grid.  This can be done by passing the '-c' flag.  An example call can be found below:

        $ python spamms.py -i input.txt -c
        $ python spamms.py -i input.txt -c -n 4

        $ python spamms.py -i input.txt
        $ python spamms.py -i input.txt -n 4
