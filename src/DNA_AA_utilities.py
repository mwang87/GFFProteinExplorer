'''
Created on Mar 24, 2016

@author: s3cha
@updater: Mingxun Wang

This class contains utility functions that allow us to transform
DNA to AA sequences

'''



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

#Given input list of dna base pairs, spits out
#a list of amino acid sequences
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
