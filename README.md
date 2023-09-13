# State_of_Tomales_Bay_birds
Report and data visualizations for monitoring data on waterbirds (shorebirds, herons and egrets, seaducks, loons, grebes, and waterfowl) at Tomales Bay, Marin County, California.

Long term monitoring of these birds began around 1990 and 


Within this broad group of waterbirds, we can split things up a little bit more based on what time of year these birds are here and what habitats they use. During the winter we monitor the shorebirds like sandpipers, plovers, and yellowlegs. These birds mostly use the tidal wetlands and mudflats on the bay. Also in the winter we count birds using the open water of the bay, including seaducks, loons, grebes, cormorants, pelicans, and other waterfowl. During the summer we monitor nesting activity of herons and egrets.  Because of these groupings, we use slightly different data collection and analysis methods, so there is a separate data processing script for each these three groups:

herons_and_egrets.R
shorebirds.R
open_water_birds.R

Each of these scripts must be run whenever new data comes in, and each script saves processed data in the data folder. Once these processed data files are updated, then output_report.RMD can be knitted. This creates a .doc mini-report with plots and brief text, and also saves each plot as a stand alone .png (saved to the figures folder).

For each group of waterbirds, a plot is made showing estimated average abundance from a general, good fitting linear model (slightly different model for each group, but same model for each species within each group), plotted over the raw data, for all species in that group combined plus some individual species. Plots are also made showing the effect strength of predictor variables included in the model.

