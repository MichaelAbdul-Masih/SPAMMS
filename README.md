# SPAMMS 1.1.1

Introduction
------------
SPAMMS stands for Spectroscopic PAtch Model for Massive Stars and is designed with geometrically deformed systems in mind.  SPAMMS combines the eclipsing binary modelling code PHOEBE 2 and the NLTE radiative transfer code FASTWIND to produce synthetic spectra for systems at given phases, orientations and geometries.  For details on how the code works, please see the corresponding release paper: [Abdul-Masih et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020arXiv200309008A/abstract)

Installation
------------
*   Clone the git repository to create a local copy.

        git clone https://github.com/MichaelAbdul-Masih/SPAMMS.git

*   SPAMMS is written in Python 3 and has several dependencies that are needed to make it run.  For your convenience, we have provided an environment file (SPAMMS_environment.yml), which you can use to create an environment that already has all of these dependencies.  If you are using Anaconda you can create the environment as follows:

        cd SPAMMS/
        conda env create -f SPAMMS_environment.yml
        conda activate SPAMMS_env

*   If you prefer to create the environment yourself, the minimum package requirements can be found below:

        astropy 2.0.12
        numpy 1.16.3
        matplotlib 3.6.1
        schwimmbad 0.3.0
        scipy 1.2.0
        tqdm 4.31.1
        phoebe 2.3.0

### Test the Installation
To make sure the installation was successful, cd into the SPAMMS git directory and run the following command:

        python spamms.py -i input.txt

If the final message you get from the run is 'All models ran successfully!' then SPAMMS has been installed successfully.


Getting Started
---------------
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

### Input spectra
If you wish to compare SPAMMS synthetic line profiles with observed spectra, SPAMMS has the ability to do so built in.  By providing the path to a folder containing the input spectra and corresponding times, SPAMMS will generate line profiles corresponding to the times provided and calculate the total chi-square for each combination of parameters specified in the input file.

The folder itself must contain two files, one corresponding to the times of observation and one corresponding to the spectra.  SPAMMS requires that these files be named with the same base name and a specific extension for each file (i.e., BASENAME_hjd.txt and BASENAME_spec.txt).  The times file (BASENAME_hjd.txt) should contain a list of times with the same format that is used for the t0 in the input file.  The spectra file (BASENAME_spec.txt) should contain one wavelength column and then several spectra columns corresponding to the number of times in the BASENAME_hjd.txt file.  It is assumed that the wavelength arrays for each of the spectra are identical.  In instances where this is not the case, several pairs of _spec.txt and _hjd.txt files can be provided.

### FASTWIND Grid
Before using the code, an input grid will need to be either downloaded or created. A large grid has been computed for the LMC and is available on request (michael.abdulmasih@gmail.com).  Alternatively, a smaller (~12GB, ~5.5GB in tar format) demo grid is also available to get started with the code (https://www.dropbox.com/s/cclg8en16n1qkmc/Demo_Grid.tar.gz?dl=0).  In a future release, we plan to make the input grid calculation scripts available.

### Running the code
SPAMMS can be run on a single core or on multiple cores by passing the '-n' flag followed by the number of cores you wish to use.  SPAMMS is parallelized on a model basis meaning that when running a grid of models with different input parameters, each can be run on a different core.  Currently, it is not possible to parallelize based on the time points for a given model.  By default, SPAMMS uses input.txt as the input file but this can be changed using the '-i' flag and specifying your chosen input file. Before running SPAMMS, it is advised to run a grid check to make sure that your model does not fall out of the bounds of the grid.  This can be done by passing the '-c' flag.  An example call can be found below:

        $ python spamms.py -i input.txt -c
        $ python spamms.py -i input.txt -c -n 4

        $ python spamms.py -i input.txt
        $ python spamms.py -i input.txt -n 4


Example input files and Jupyter notebooks
---------------
In the Demo folder, we have provided 4 example input files for different types of systems.  This is a great place to start if you are a new user.  These input files all use the demo grid discussed in the previous section, so make sure that it is downloaded and in the Grids/ folder before running them.  We recommend copying each of the files into the main SPAMMS folder and running them there to preserve the originals for future reference.  In order of increasing complexity, the input files are:

*   input_demo_rapid_rotator.txt: Here we demonstrate how to model a rapidly rotating single star at different rotation rates

*   input_demo_semidetached.txt: Here we demonstrate how to model a semidetached binary at different times.

*   input_demo_Rossiter-McLaughlin.txt: Here we show an example of how SPAMMS can model the Rossiter-McLaughlin effect

*   input_demo_overcontact_fit.txt: Here we demonstrate how you can use SPAMMS to fit observed spectra.  In this case, we simulated an overcontact binary using SPAMMS, added noise and degraded the spectra to XSHOOTER resolution (R = 6700) to form the "observed" spectra.


In addition to the sample input files, we've also included two jupyter notebooks, which demonstrate how to handle and work with the SPAMMS outputs.  The first shows how to access and plot the line profiles and uses the output from the input_demo_rapid_rotator.txt for this example.  The second shows how to create animations using a combination of the line profile outputs from SPAMMS and the geometry information generated using PHOEBE 2.  This notebook uses the output generated from the input_demo_semidetached.txt file to create the animation.
