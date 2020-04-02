# DNA Sequence Analysis

An aligned DNA sequence analyzer that remotely runs the BLAST algorithm on the [NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) database and runs analysis on the retrieved XML data.

## Requirements

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install [Biopython](https://biopython.org). 

```bash
pip install biopython
```

## Usage

Enter the accession number (ex. "MT263381") of the nucleotide you would like to run BLAST on to NCBI-BLAST.py

```python
NCBIWWW.qblast("blastn", "nt", "MT263381", hitlist_size=200,  format_type='XML')
```

Then use the Genome.py methods to run analysis.

```python
g=Genome(aligned_sequences)
g.protein_differences()
g.amino_acid_differences()
g.consensus_sequence()
```

## License
[MIT](https://choosealicense.com/licenses/mit/)