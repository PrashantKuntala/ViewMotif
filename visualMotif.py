#!/usr/bin/python
#
# author : prashant kumar kuntala
# date   : 19th April, 2017
#
# last modified : 19th April, 2017
#

# script to highlight the motif region within the promoter sequences. Generates a html file as an output showing the visualization.

# usage : python visualMotif.py < path to fimo.txt > < path to the promoter fasta file >

import sys,os,re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


# list of hexcodes for coloring the regions
# hexcodes = ["#C0392B","#e67e22","#f1c40f","#2ecc71","#3498db","#2980b9","#8e44ad","#FF0000","#FFFF00","#008000","#000080","#800080"]

hexcodes = ["#f44336","#9c27b0","#e91e63","#2196f3","#009688","#4caf50","#8bc34a","#ffeb3b","#ffc107","#03a9f4","#ff9800","#ff5722","#795548","#00bcd4","#607d8b","#673ab7","#3f51b5","#cddc39"]

# html head
head = """
<html>
<title>VisualMotif</title>
<body>
<style>



html {
        padding : 6px;
        font-family: Arial, Helvetica, sans-serif;
        margin: 12px;
        background-color: #ede3de;
}

.maindiv { padding-top: 250px;
}

h1 {
  padding-top: 20px;
  font-family: Garamond;
}

.seqdiv { word-wrap: break-word;
      font-size: 0.9em;
      padding: 12px;
      margin-right:  10px;
      text-align:justify;
      background-color: #273746;
      color: white;
}

span { font-size : 1em;
        font-style : italic;
        padding-left: 3px;
        padding-right: 3px;
        color:black;
}

.mydiv {
  display: inline-block;
  background-color: #273746;
  padding: 10px;
  position: fixed;
}

.mydiv div {
  display: inline-block;
  /*background-color: #a0a0aa;*/
  margin-left: 15px;
  margin-right: 15px;
  /*border: 2px solid green;*/
}

.mydiv div p {
            display: inline-block;
            margin-left: 4px;
            color: white;
            font-family: Garamond;
}

.mydiv div div{
            display:inline-block;
            padding: 5px;
            margin:0px;
}



</style>
"""

# html tail
tail = """
</body>
</html>
"""

# dictionary containing  [ seq-> motif -> [ starting value ] ]
startmain = dict()
stopmain = dict()
ranges = dict()

# contains only the start locations and only stop locations.
allstarts = list()
allstops = list()
both = dict()

# to hold all the unique motif names and sequence names
motifs = set()
sequence = set()

# contains html for the legend  [ Motifname : color ].
legend = ""

# dictionary that holds the [ motif: hexvalue ]
color = dict()

# list to store all the gene names
allgenes = []

def createDictionary(mydict,motifs,sequence,allgenes):
    """
    assigns default values to the dictionary passed so that you populate the dictionary with right values

    structure :  each sequence has 'm' motifs in them. Each motif 'm' has a start and stop index.
    """
    for s in sequence:
        seq = s.split("|",1)[0]
        allgenes.append(seq)
        mydict[seq] = {}
        for m in motifs:
            mydict[seq][m] = []

def createSequenceDictionary():
    """
    Function that returns the dictionary with keys as sequence names and values as Sequences.
    """

    # To store the genenames.
    genes = []

    # reading the input fasta file and creating a temp file to create the promoter dictionary
    openfile =  open(sys.argv[1],'r').readlines()
    outfile = open("./temp.fasta","w")
    for line in openfile:
        if(re.search(">",line)):
            ln = line.split("|",1)[0]
            genes.append(ln.split('>',1)[1])
            outfile.write(ln+"\n")
        else:
            outfile.write(line)
    outfile.close()

    try:
        # dictionary that stores the genenames and corresponding promoter sequences.
        seqdict = SeqIO.to_dict(SeqIO.parse("./temp.fasta","fasta"))
    except ValueError as ve:
        print "The fasta file contains duplicate genenames\n" + str(ve)
        exit()

    os.remove("./temp.fasta")   # deleting the file once we have the dictionary

    return seqdict

def processFimo(stm,spm,m,s,ran,allgenes):
    """
    Function to read and process fimo.txt file. To setup the start and stop dictionaries.

    parameters :
    stm = startmain
    spm = stopmain
    ran = ranges
    """

    # reading the fimo.txt file
    f = pd.read_csv(sys.argv[2],delimiter="\t")

    # to hold all the unique motif names and sequence names
    m = set(f["#pattern name"])
    s = set(f["sequence name"])

    # populating the respective dictionaries
    createDictionary(stm,m,s,allgenes)
    createDictionary(spm,m,s,allgenes)
    createDictionary(ran,m,s,allgenes)

    allgenes = set(allgenes)
    allgenes = list(allgenes)
    allgenes.sort()

    return [list(m),s,f,allgenes]

def assignColors(motifs,hexcodes,fimo):
    """
    function to assign a unique color to list of unique motifs.
    """
    color = {}
    count = 0
    legend = ""
    for l in set(fimo["#pattern name"]):
        color[l] = hexcodes[count]
        count = count+1
        legend = legend + '<div><div  style=\"height:10px;width:10px;background-color:'+color[l]+'\"></div>' +'<p>'+ str(l)+'</p></br></div>'
        # print color[l]
    legend = '<div class="mydiv">' + legend + "</div>"
    return [color,legend]

def populateDictionary(startmain,stopmain,fimo,ranges):
    """
    Function to populate the start and stop dictionaries
    """

    for i in range(0,len(fimo)):
        sname = fimo["sequence name"][i].split("|",1)[0]
        motif = fimo["#pattern name"][i]
        start = fimo["start"][i]
        stop = fimo["stop"][i]
        # print "{},{}".format(sname,motif)
        startmain[sname][motif].append(start)
        stopmain[sname][motif].append(stop)
        ranges[sname][motif].append(range(start,stop+1))

def lookUp(position,gene):
    """
    function to look up if the particular position has an overlap or not in that particular gene.
    returns the motif name and no.of overlaps at that position if there is no overlap.
    returns the no.of overlap at the position by default
    """
    count = 0
    m = []
    for motif in ranges[gene]:
        # print k
        mylist = ranges[gene][motif]
        for l in mylist:
            if(position in l):
                count = count + 1
                m.append(motif)
                # print " value : {}\n".format(l)
    # print "count in gene : {} is {}\n\n".format(gene,count)

    return count,m

def createHTML(head,legend,promoter,tail):
    """
    Function to create the FINAL HTML.
    """
    # adding the header to the html file.
    htmlfile = open("./VisualMotif.html",'wb')
    htmlfile.write(head+"\n")
    htmlfile.write(legend+"\n")
    htmlfile.write("<div class =\"maindiv\">")

    # for each promoter sequence
    for genes in allgenes:
        # print genes

        htmlfile.write("<h1>"+str(genes)+"</h1>")

        # retrive the sequence as a list
        mypromoter = list(promoter[genes].seq)
        string = ""

        # for each nucleotide in the promoter.
        for i in range(0,len(mypromoter)):
            c,mot = lookUp(int(i+1),str(genes))
            # print "count :{} \t motif : {}".format(c,mot)     ## DEBUG to check the overlap and motif name

            # adding tags according to the count of overlap
            if(c == 0):
                string = string + mypromoter[i]
            if(c == 1):
                string = string + '<span style=\"background-color:'+color[mot[0]]+';\">' + str(mypromoter[i]) + "</span>"
            if(c > 1 and len(mot) >1):
                if(mot[0] != mot[1]):
                    string = string + '<span style=\"border-top: 2px solid '+color[mot[0]]+';border-bottom:2px solid '+ color[mot[1]]+';font-style: oblique;color:white;background-color:black;\">' + str(mypromoter[i]) + "</span>"
                else:
                    string = string + '<span style=\"border-top: 2px solid white;border-bottom:2px solid white;background-color:'+color[mot[0]]+';font-style: oblique;\">' + str(mypromoter[i]) + "</span>"

        htmlfile.write("<div class=\"seqdiv\">"+string+"</div>"+"\n")

    htmlfile.write("</div></br>")
    htmlfile.write(tail)
        # print " \n\n"
    # print list(promoter['Arhgap25'].seq)


#################################################################################################################
# Function calls.
################################################################################################################

promoter = createSequenceDictionary()

# # # DEBUG
# # printing the promoter sequence using the genenames
# print promoter["Mir3572"].seq[1:20]
# # print promoter["AT2G18470"].seq[1:20]

[motifs,sequence,fimo,allgenes] = processFimo(startmain,stopmain,motifs,sequence,ranges,allgenes)

# DEBUG : To check the structure is properly created.
# for k,v in stopmain.items():
#     print k
#     for i,j in v.items():
#         print "{}:{}".format(i,j)
#
# # DEBUG
# for m in motifs:
#     print m

#
[color,legend] = assignColors(motifs,hexcodes,fimo)

populateDictionary(startmain,stopmain,fimo,ranges)
#
# # # DEBUG : To check whether the start and stop are being retrieved correctly for a particular sequence.
# for k,v in startmain["Mir3572"].items():
#     # print "{}\t{}".format(k,startmain["Cnot3"][k])  ## Ttotal of 8 starts
#     # print "{}\t{}".format(k,stopmain["Cnot3"][k])  ## Ttotal of 8 starts
#     print "{}\n{}\n{}".format(k,ranges["Mir3572"][k],v)
#     if(1385 in range(12,1388)):
#         print "\n\nyes"

createHTML(head,legend,promoter,tail)
