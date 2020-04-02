import xml.etree.ElementTree as ET
import xml_files

# Opens an XML output containing aligned sequences, produced by NCBI's BLAST software:
## https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome

aligned_sequences=[]
query_name=''
subject_name=''

# parse an xml file by name - files are stored in the 'xml_files' folder.
tree = ET.parse('xml_files/blast3.xml')
root = tree.getroot()

query_sequence=''
subject_sequence=''
similarity_score=[]
similarity=''

# Iterates through the XML nodes and stores the query sequence and its name,
## and the subject sequence and its name. 
for elem in root.iter():
    if(elem.tag=='BlastOutput_query-def'):
        query_name=elem.text
    if(elem.tag=='Hit_def'):
        subject_name=elem.text
    if(elem.tag=='Hsp_qseq'):
        query_sequence=elem.text
    if(elem.tag=='Hsp_hseq'):
        subject_sequence=elem.text
    if(elem.tag=='Hsp_identity'):
        similarity_score.append(elem.text)
    if(elem.tag=='Hsp_positive'):
        similarity=str(int(elem.text)/int(similarity_score[0]))
        similarity_score.clear()
    if(similarity!='' and query_name!=''and subject_name!='' and query_sequence!='' and subject_sequence!=''):
        aligned_sequences.append([query_name, query_sequence,subject_name, subject_sequence, similarity])
        subject_name=query_sequence=subject_sequence=''

## A method for identifying the name of the amino acid (by DNA codon).
## Taken from: https://www.geeksforgeeks.org/dna-protein-python-3/
def translate(seq): 
       
        table = { 
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
        }
        protein =""
        if len(seq)%3 == 0: 
            for i in range(0, len(seq), 3): 
                codon = seq[i:i + 3]
                try:
                    protein+= table[codon]
                except:
                    protein=seq[i:i + 3]
                    break
        return protein 
