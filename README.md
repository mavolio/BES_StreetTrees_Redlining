# BES_StreetTrees_Redlining

3 files are necessary to replicate the analysis found in:

Karin Burghardt, Meghan Avolio, Dexter Locke , Morgan Grove, Nancy F. Sonti, and Christopher M Swan. Current Street Tree Communities Reflect Race-based Housing Policy and Modern Attempts to Remedy Environmental Injustice. Ecology XX:XXX (accepted, in press)


List of files:

      link_street_trees_hex_holc.Rmd
      
      05.09.2022_RST_hex_occupancy_and_div_accum_code.R
      
      beta diversity street tree code.R

Detailed description of workflow

      1). link_street_trees_hex_holc.Rmd: 
      
      This code uses polygon files of HOLC risk grades and joins it with tree inventory data from the Baltimore tree inventory to prune the data to just trees in rated neighborhoods and form the hexagons used in subsequent analysis.

            input data: 
                  
                  A 2017-2018 tree inventory of Baltimore, MD publicly available to download at       https://baltimore.maps.arcgis.com/apps/webappviewer/index.html?id=d2cfbbe9a24b4d988de127852e6c26c8
                                
                  GIS shapefiles for HOLC grade designations publicly available under a creative commons license download from the “Mapping Inequality,” American Panorama Project at https://dsl.richmond.edu/panorama/redlining/#loc=10/39.676/-77.019&text=downloads by searching for Baltimore, MD.
     
            output:street_trees_Baltimore_w_HOLC_grades_HEX_2021-03-15.csv

      2). 05.09.2022_RST_hex_occupancy_and_div_accum_code.R: 
      
      This code a) cleans the tree data down to only street trees trees and fixes species issues b) generates an output file used as an input file of clean data for the beta diversity street tree code.R c) performs an occupancy analysis d) generates species accumulation curves at the HOLC grade scale
      
             input: street_trees_Baltimore_w_HOLC_grades_HEX_2021-03-15.csv
      
             output: st_tree_inHEX2021_03_15.csv

      3). beta diversity street tree code.R :
      
      This code performs the beta diversity analysis including the generation of rank abundance curves, NMDS ordinations, and codyn compositional analysis of turnover and re-ordering
      
            input: st_tree_inHEX2021_03_15.csv
