# ViewMotif
A python script to visualize DNA sequence motifs within your promoter sequences. Generates a HTML file as an output. Allows you to visualize overlapping DNA motifs.

# Usage

python visualmotif.py  [path to promoter sequences in fasta format] [ path to fimo.txt file] 

example : python visualMotif.py promoters.fasta sample_fimo.txt 


# Dependencies

you will need biopython, pandas installed on your machine.

you can use the sample promoter sequences file and sample fimo file to see the output.

VisualMotif.html is the sample output that will be generated if you run the script on the provided sample files.


# Details 

To understand the output use the legend below.

legend :

solid white top and bottom border : implies that there is a self overlap within the motif and the region with border is the overlap.

black background with different top and bottom border : implies that ther is an ovelap between two different motifs.



# Note

you can generate the fimo.txt using FIMO tool in the meme-suite

meme-suite :  http://meme-suite.org/
fimo	   :  http://meme-suite.org/tools/fimo

if you can download the entire meme-suite from http://meme-suite.org/doc/download.html?man_type=web then you can run FIMO locally on your machine, refer its documentation for usage.


