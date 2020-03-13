#import packages

#third party regex with overlapping capability finditer
import regex as re
import cairo
import argparse
import math


#argparse to allow user to pass in files
def get_args():
    parser = argparse.ArgumentParser(description="Marker of motif locations in a gene. Exons appear as boxes, and motifs as hashmarks.")
    parser.add_argument("-f", "--fasta_genes", help="Pass in a fasta file containing the genes to be scanned with -f. The sequence for each gene in the file should be on a single line. The file must have no whitespace above the first entry.", required=True)
    parser.add_argument("-m", "--motif_file", help="Use -m to pass in a file containing the motifs to be searched for, with each motif separated by a new line.", required=True)
    return parser.parse_args()

args = get_args()
genes_file = args.fasta_genes
motifs_file = args.motif_file


#initialize dictionaries and lists
gene_name_dic = {}
motif_dic={}
gene_dic_startpos={}
gene_dic_endpos={}
gene_len_dic={}
motifs_list=[]
pattern_list=[]

#Make list for each motif sequence
with open(motifs_file, "r") as fh:
    for line in fh:
        stripped=line.strip()
        motifs_list.append(stripped)

#fill dictionary with the regex matches for each possible character in the motifs
charachter_dic = {}
charachter_dic["c"]=["c,C"]
charachter_dic["C"]=["c,C"]
charachter_dic["g"]=["g,G"]
charachter_dic["G"]=["g,G"]
charachter_dic["t"]=["t,T"]
charachter_dic["T"]=["t,T"]
charachter_dic["a"]=["a,A"]
charachter_dic["A"]=["a,A"]
charachter_dic["u"]=["t,T"]
charachter_dic["U"]=["t,T"]
charachter_dic["y"]=["c,t,C,T"]
charachter_dic["Y"]=["c,t,C,T"]

#Use both the character dictionary and the lists of motifs to create a single string to be used in a regex search
for item in motifs_list:
    pattern=""
    i=0
    for charachter in item:
        list = charachter_dic[charachter]
        string = ','.join(list)
        pattern += "["+re.escape(string)+"]"
    pattern_list.append(pattern)
    i+=1
#Find motif and exon locations within gene files
with open(genes_file, "r") as fh:
    i=1
    while True:
        L1 = fh.readline().strip()
        if L1 == "":
            break
        name="name"+str(i)
        gene_name_dic[name]=str(L1)
        L2 = fh.readline().strip()
        motif_dic[name]={}
        for index in range(len(pattern_list)):
            search = pattern_list[index]
            for m in re.finditer(rf'{search}', L2, overlapped=True):
                number=index+1
                label="m"+str(number)
                motif_dic[name][m.start()]=label
        #find exon_pos
        for gene in re.finditer(r'[G,C,A,T]+', L2):
            gene_startpos = gene.start()
            gene_dic_startpos[name] = gene_startpos
            gene_endpos = gene.end()
            gene_dic_endpos[name] = gene_endpos
        gene_len_dic[name] = len(L2)
        i+=1



#Create coordinates normalized to gene length
#normal_dic = {}
#for n in range(len(motif_dic)):
#    index = n +1
#    keyy = "name" + str(index)
#    normal_dic[keyy] = {}
#    for key in motif_dic[keyy]:
#        key2 = key*795/gene_len_dic[keyy]
#        normal_dic[keyy][key2]=motif_dic[keyy][key]

#Create exon coordinates normalized to gene length

#gene_dic_startpos_normal = {}
#n=0
#for key in gene_dic_startpos:
#    index = n +1
#    keyy = "name" + str(index)
#    value2 = gene_dic_startpos[key]*795/gene_len_dic[keyy]
#    gene_dic_startpos_normal[key]=value2
#    n+=1
#n=0
#gene_dic_endpos_normal = {}
#for key in gene_dic_endpos:
#    index = n +1
#    keyy = "name" + str(index)
#    value2 = gene_dic_endpos[key]*795/gene_len_dic[keyy]
#    gene_dic_endpos_normal[key]=value2
#    n+=1

#make cairo surface
ylen_surface = len(gene_name_dic)*100+100+len(motifs_list)*20
surface = cairo.SVGSurface("plot.svg", 795, ylen_surface)
ctx = cairo.Context(surface)


#generate colors based on number of motifs
color_dic = {}
for n in range(len(motifs_list)):
    i = n+1
    R = abs(math.cos(math.cos(i*1.7)+math.cos(i*1.7+(1/len(motifs_list)))))
    G = abs(math.cos(math.cos(i*2.1)+math.tan(i*2+(1/len(motifs_list)))))
    B= abs(math.sin(math.sin(i*2.8)+math.sin(i*2.8+(1/len(motifs_list)))))
    key = "m"+str(i)
    color_dic[key]=[R,G,B]

for n in range(len(motif_dic)):

    index=n+1
    keyy = "name" + str(index)

    #Draw Line for Gene
    ctx.set_line_width(2)
    ypos = n*100+50
    ctx.move_to(0,ypos)
    ctx.line_to(gene_len_dic[keyy],ypos)
    ctx.stroke()

    #draw gene box
    ctx.set_line_width(2)
    gene_startpos = gene_dic_startpos[keyy]
    gene_endpos = gene_dic_endpos[keyy]
    ctx.rectangle(gene_startpos,ypos-15,gene_endpos-gene_startpos,30)
    ctx.stroke()

    #Type gene name
    ctx.set_font_size(14)
    ctx.move_to(10,ypos-30)
    gene_name = gene_name_dic[keyy]
    ctx.show_text(gene_name)


#Make motif key boxes
n=0
for m in color_dic:
    ypos=len(motif_dic)*100+10+n*20
    R = color_dic[m][0]
    G = color_dic[m][1]
    B = color_dic[m][2]
    ctx.set_source_rgb(R,G,B)
    ctx.rectangle(15,ypos,15,15)
    ctx.fill()
    ctx.stroke()
    ctx.move_to(40,ypos+10)
    motif_name=motifs_list[n]
    ctx.show_text(motif_name)
    n+=1

#Make Legend
ypos=(len(motif_dic)-1)*100+100
ctx.set_source_rgb(0,0,0)
y_span = len(motifs_list)*20-20+10+15+10
ctx.rectangle(10,ypos,300,y_span)
ctx.stroke()

#Draw motif marks
for n in range(len(motif_dic)):
    index=n+1
    ypos = n*100+50
    keyy = "name" + str(index)
    for key in motif_dic[keyy]:
        ctx.set_line_width(2)
        R = color_dic[motif_dic[keyy][key]][0]
        G = color_dic[motif_dic[keyy][key]][1]
        B = color_dic[motif_dic[keyy][key]][2]
        ctx.set_source_rgb(R,G,B)
        xpos = key
        ctx.move_to(xpos,ypos)
        ctx.line_to(xpos,ypos+10)
        ctx.line_to(xpos,ypos-10)
        ctx.stroke()
#write out file
surface.finish()
