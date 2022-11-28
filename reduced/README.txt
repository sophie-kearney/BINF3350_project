This directory contains reduced GO data from the top 400 GO entries for each cancancer based on
adjusted p-values. 

The <cancer>_reduced_all.csv files contain the direct output from the table provided by Revigo. The
rows that contain a null value for the Representative column are cluster representative values. Rows
that contain an integer refer to the GO termID that is that row's representative. For example, if a
row has 280 under the Representative attribute, that row's parent value is "GO:0000280".
