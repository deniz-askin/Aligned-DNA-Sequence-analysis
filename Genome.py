from xml_methods import translate
from xml_methods import aligned_sequences

## Genome class takes aligned nucleotide sequences as argument,
## and runs statistical analysis on them.

class Genome:
    def __init__(self, sequence):
        self.sequence=sequence

    ## Prints names and similarities of the sequences being analyzed.
    def report(self):
        LIM = 6
        for s in self.sequence:
            if(s[1]==s[3]):
                print("Names of Aligned Sequences: " + str(s[0]) + "|||||||" + str(s[2]) + '\n' + "Aligned sequences are identical")
            else:
                same=[]
                different=[]
                for k,z in zip(s[1],s[3]):
                    if(k==z):
                        same.append(k)
                    else:
                        different.append((k,z))
                        pass
                print("Names of Aligned Sequences: " + str(s[0]) + "|||||||" + str(s[2]) + '\n' + "Aligned sequences are not identical" + '\n' + "Similarity: "+str(len(same)/len(s[1]))[:LIM] + '\n' + "Protein Mutation: " + str(different))

    ## Returns SNPs (single-nucleotide polymorphisms) and their number of occurence.
    ## Prints the most and least frequent SNPs and their frequency of occurences.
    def protein_differences(self):
        different=[]
        count=[]
        amino_acids=[]
        for s in self.sequence:
            if(s[1]==s[3]):
                pass
            else:
                for k,z in zip(s[1],s[3]):
                    if(k!=z):
                        different.append((k,z))
                        pass
                    if(k=="-" or z=="-"):
                        amino_acids.append((s[1][s[1].index(k)-1],s[1][s[1].index(k)], s[1][s[1].index(k)+1], s[3][s[3].index(z)-1],s[3][s[3].index(z)], s[3][s[3].index(z)+1]))
                    
        def maximum_occuring():
            maximum=[max(set(different), key = different.count), " Frequency of occurence: ", str(different.count(max(set(different), key = different.count))/len(self.sequence))]    
            print("Most frequent SNPs: ",  maximum)
            
        def minimum_occuring():
            minimum=[min(set(different), key = different.count), str(different.count(min(set(different), key = different.count))/len(self.sequence))]
            print("Least frequent SNPs: ", minimum)

        maximum_occuring()
        minimum_occuring()

        for s in list(set(different)):
            count.append(s)
        count.sort(key=different.count, reverse=True)
        
        for s in range(0,len(count)):
            count[s]=[count[s], different.count(count[s])]
            
        return count

    ## Returns amino acid mutations and their number of occurence.
    ## Prints the most and least frequent amino acid mutations and their frequency of occurences.
    def amino_acid_differences(self):
        different=[]
        count=[]
        for s in self.sequence:
            if(s[1]==s[3]):
                pass
            else:
                for k,z in zip(s[1],s[3]):
                    if(k!=z):
                        different.append((s[1][s[1].index(k)-1],s[1][s[1].index(k)], s[1][s[1].index(k)+1], s[3][s[3].index(z)-1],s[3][s[3].index(z)], s[3][s[3].index(z)+1]))
                        pass
                    
        def maximum_occuring():
            maximum=[max(sorted(set(different), key=different.index), key = different.count), " Frequency of occurence: ", str(different.count(max(set(different), key = different.count))/(len(self.sequence)/3))]
            amino_acid1=''
            amino_acid2=''
            for s in range(0,len(maximum[0])):
                if(s<=2):
                    amino_acid1+=maximum[0][s]
                else:
                    amino_acid2+=maximum[0][s]
            print("Most frequent Amino acid Mutations: ", translate(amino_acid1), " to ", translate(amino_acid2), maximum)
        
        def minimum_occuring():
            minimum=[min(sorted(set(different), key=different.index), key = different.count), " Frequency of occurence: ", str(different.count(min(set(different), key = different.count))/(len(self.sequence)/3))]
            amino_acid1=''
            amino_acid2=''
            for s in range(0,len(minimum[0])):
                if(s<=2):
                    amino_acid1+=minimum[0][s]
                else:
                    amino_acid2+=minimum[0][s]
            print("Least frequent Amino acid Mutations: ", translate(amino_acid1), " to ", translate(amino_acid2), minimum)

        maximum_occuring()
        minimum_occuring()
        
        for s in list(set(different)):
            count.append(s)
        count.sort(key=different.count, reverse=True)
        for s in range(0,len(count)):
            count[s]=[count[s], different.count(count[s])]
        return count

    ## Returns the consensus sequence for all the aligned sequences.
    ## Prints the similarity of the query sequence to the consensus sequence.
    def consensus_sequence(self):
        different_sequences=[]
        most_repeated_protein=[]
        cons=''
        
        for s in self.sequence:
            if(s[1]!=s[3]):
                different_sequences.append(s[1])
                different_sequences.append(s[3])

        for k in range(0, len(self.sequence[0][1])):
            tup=()
            for s in different_sequences:
                try:
                    tup+=(s[k],)
                except:
                    pass
            most_repeated_protein.append(tup)

        for s in most_repeated_protein:
            most=max(s, key = s.count)
            cons+=most

        def similarity():
            count=0
            for k,s in zip(self.sequence[0][1],cons):
                if(k!=s):
                    count+=1
                    
            print("Similarity of the query sequence to the consensus sequence: ", round(100*(len(self.sequence[0][1])-count)/len(self.sequence[0][1]),2))

        similarity()

        return cons

g=Genome(aligned_sequences)
##g.report()
g.protein_differences()
g.amino_acid_differences()
g.consensus_sequence()

