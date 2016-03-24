'''
Created on Mar 24, 2016

@author: s3cha
'''

import sys


ForwardCode = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
               "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
               "TAT":"Y", "TAC":"Y", "TAA":"X", "TAG":"X",
               "TGT":"C", "TGC":"C", "TGA":"X", "TGG":"W",
               "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
               "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
               "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
               "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
               "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
               "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
               "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
               "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
               "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
               "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
               "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
               "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
               }
ReverseCode =   {"AAA":"F", "AAG":"F", "AAT":"L", "AAC":"L",
               "AGA":"S", "AGG":"S", "AGT":"S", "AGC":"S",
               "ATA":"Y", "ATG":"Y", "ATT":"X", "ATC":"X",
               "ACA":"C", "ACG":"C", "ACT":"X", "ACC":"W",
               "GAA":"L", "GAG":"L", "GAT":"L", "GAC":"L",
               "GGA":"P", "GGG":"P", "GGT":"P", "GGC":"P",
               "GTA":"H", "GTG":"H", "GTT":"Q", "GTC":"Q",
               "GCA":"R", "GCG":"R", "GCT":"R", "GCC":"R",
               "TAA":"I", "TAG":"I", "TAT":"I", "TAC":"M",
               "TGA":"T", "TGG":"T", "TGT":"T", "TGC":"T",
               "TTA":"N", "TTG":"N", "TTT":"K", "TTC":"K",
               "TCA":"S", "TCG":"S", "TCT":"R", "TCC":"R",
               "CAA":"V", "CAG":"V", "CAT":"V", "CAC":"V",
               "CGA":"A", "CGG":"A", "CGT":"A", "CGC":"A",
               "CTA":"D", "CTG":"D", "CTT":"E", "CTC":"E",
               "CCA":"G", "CCG":"G", "CCT":"G", "CCC":"G",
               }

def getAA(output_array):
    output_seq = []
    for j in range(3):
        AA_seq = ''
        for i in range(len(output_array[j:])/3):
            if 'N' in output_array[3*i+j:3*(i+1)+j]:
                AA_seq += 'X'
            else:
                AA_seq += ForwardCode[output_array[3*i+j:3*(i+1)+j]]
        output_seq.append(AA_seq)
    return output_seq


#Returns a dictionary by parsing the GFF gff
# Key is chromosome -> a list of genes in that chromosome
# For each of these genes, there is a list of attributes
# The attributes are begin, end, name, and a list of CDS
# Each CDS element represents an exon and will have a beginning, end, and sequence of each exon
def parseGFF(gff_filename):
    f = open(gff_filename,'r')
    ref_seq_gff_dic = {}
    pseudo_gff_dic = {}
    gene_gff_dic = {}
    
    prev_CDS_parent_ID = ""
    prev_CDS_start = -1
    prev_CDS_end = -1
    prev_pseudo_ID = ''
    
    for line in f:
        curr_part = line.strip()
        curr_part = curr_part.split('\t')
        if len(curr_part) != 9:
            continue
        chr_str = curr_part[0]
        curr_part[3] = int(curr_part[3]) - 1  # 0 base inclusive transfer
        curr_part[4] = int(curr_part[4])
    
        if curr_part[2] == "mRNA":
            data = curr_part[-1].split(";")
            for i in range(0, len(data)):
                if data[i].startswith("ID="):
                    gene_name_str = data[i].replace("ID=", "")
                if data[i].startswith("Name="):
                    gene_readable_name_str = data[i].replace("Name=", "")        
                    break
            if ref_seq_gff_dic.has_key(chr_str):
                ref_seq_gff_dic[chr_str].append([curr_part[3], curr_part[4], gene_name_str , curr_part, [], [], [], [], gene_readable_name_str])
                curr_ind = len(ref_seq_gff_dic[chr_str]) - 1
            else:
                ref_seq_gff_dic[chr_str] = []
                ref_seq_gff_dic[chr_str].append([curr_part[3], curr_part[4], gene_name_str , curr_part, [], [], [], [], gene_readable_name_str])
                curr_ind = 0
                # [start, end, gene_name_str, mRNA, CDS, five_prime_UTR, three_prime_UTR, splice junctions, gene_readable_name_str]
        else:
            data = curr_part[-1].split(";")
            parent_name_str = ""
            for i in range(0, len(data)):
                if data[i].startswith("Parent="):
                    parent_name_str = data[i].replace("Parent=", "")
                    break
    
        if curr_part[2] == "CDS":
            if parent_name_str != ref_seq_gff_dic[chr_str][curr_ind][2]:
                print "No Parent", ref_seq_gff_dic[chr_str][curr_ind][2]
                break
            ref_seq_gff_dic[chr_str][curr_ind][4].append(curr_part)
    
            if prev_CDS_parent_ID == parent_name_str:
                splice_info = []
                if prev_CDS_start < int(curr_part[3]):
                    j_start = prev_CDS_end
                    j_end = int(curr_part[3])
                else:
                    j_start = int(curr_part[4])
                    j_end = prev_CDS_start
                splice_info.append(j_start)
                splice_info.append(j_end)
                ref_seq_gff_dic[chr_str][curr_ind][7].append(splice_info)
    
            prev_CDS_parent_ID = parent_name_str
            prev_CDS_start = int(curr_part[3])
            prev_CDS_end = int(curr_part[4])
    
        elif curr_part[2] == "five_prime_UTR":
            if parent_name_str != ref_seq_gff_dic[chr_str][curr_ind][2]:
                print "No Parent", ref_seq_gff_dic[chr_str][curr_ind][2]
                break
            ref_seq_gff_dic[chr_str][curr_ind][5].append(curr_part)
    
        elif curr_part[2] == "three_prime_UTR":
            if parent_name_str != ref_seq_gff_dic[chr_str][curr_ind][2]:
                print "No Parent", ref_seq_gff_dic[chr_str][curr_ind][2]
                break
            ref_seq_gff_dic[chr_str][curr_ind][6].append(curr_part)
        
        elif curr_part[1].find('pseudogene')>-1 and curr_part[2] == 'exon':
            if curr_part[8].find('Parent=') < 0:
                continue
            if not pseudo_gff_dic.has_key(curr_part[0]):
                pseudo_gff_dic[curr_part[0]] = [] #{chr:[[start1,start2],[end1,end2],'parent_name','pseudo name']
            curr_parent_ID = curr_part[8].split('Parent=')[1]
            if curr_parent_ID != prev_pseudo_ID:
                prev_pseudo_ID = curr_parent_ID
                pseudo_gff_dic[curr_part[0]].append([[int(curr_part[3])],[int(curr_part[4])],curr_parent_ID,curr_part[1]])
            else:
                pseudo_gff_dic[curr_part[0]][-1][0].append(int(curr_part[3]))
                pseudo_gff_dic[curr_part[0]][-1][1].append(int(curr_part[4]))
        
        elif curr_part[2] == "gene":
            if not gene_gff_dic.has_key(curr_part[0]):
                gene_gff_dic[curr_part[0]] = [] #{chr:[[start],[end],'parent_name}
            curr_part_ID = curr_part[8].split('Name=')[1]
            gene_gff_dic[curr_part[0]].append([[int(curr_part[3])],[int(curr_part[4])],curr_part_ID])
        else:
            continue
    return ref_seq_gff_dic

# gff = parseGFF('/home/s3cha/data/ProteoSaFe/ProteoSAFe-1.2.6_beta-linux64/server/resources/ref_gff/GRCh37.70/Homo_sapiens.GRCh37.70_formatted.gff')
# for item in gff.get('chrX'):
#     print item
#     sys.exit()

