The code and data associated with this repository are sufficient to replicate the analyses reported in Foley et al, "Basic reversal-learning capacity in flies suggests rudiments of complex cognition."

fly_reversal.csv: a csv file containing all the data collected from the experiment
    index - a unique index for each observation
    mon - the month the experimental replicate was performed
    day - the day the experimental replicate was performed
    exp - a unique identifier for each experimental replicate
    stn - the video station at which the experimental replicate was performed
    camera - the camera (0-3) at the station which recorded the observation
    geno - the genotype of the flies observed (either B,C,F [cosmopolitan]; or X,Y,Z [caribbean])
    treat - one of 2 treatments. Either the pineapple had caffeine (pc) or the grapefruit did (gc)
    juice - the juice observed, either pineapple (p) or grapefruit (g)
    period - either 1 (the initial training period, morning) or 2 (the reversal period, evening)
    min - the number of minutes, from the beginnining of recording, at which the experiment was observed
    n - the number of flies counted in the observation
    flies - empty column


code_for_foley_et_al.r: an R file with all the code necessary to conduct the analysis and generate the figures found in the paper.
    - you will need to set your own working directory, and possibly install the appropriate dependencies.
    - there are comments, hopefully sufficient, throughout.    
    
