from Bio.Blast import NCBIWWW

## Biopython method to run the BLAST algorithm on a sequence
## remotely from the NCBI database.

## Runs the blast algorithm on a query sequence and returns
## 200 subject nucleotide sequences.


## Enter NCBI accession number as the 3rd argument below (e.g. MT263381).
result_handle = NCBIWWW.qblast("blastn", "nt", "MT263381", hitlist_size=200,  format_type='XML')
f=open("xml_files/blast3.xml","w+")
f.write(result_handle.read())
f.close()
