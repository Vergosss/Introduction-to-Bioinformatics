from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, GC_skew
import re
file = open("sequence.fasta")#anigo ena arxeio fasta
records = parse(file, "fasta")#
for record in records:
    print('Sequence name: %s ' % record.name)
    print("Sequence Data: %s" % record.seq)
    #print("Sequence Alphabet: %s" % record.seq.alphabet)

a_seq = Seq("ATGCATGCATGCATGCATGCATG")

def complementarySequence(seq):
    complement = seq.complement()

    with open('complementaryAndTranscriptingSequence.fasta', 'w') as file:
        file.write(">complementarySequence\n")
        file.write(str(complement))
    file.close()
    ####
    transcribe = seq.transcribe()
    with open('complementaryAndTranscriptingSequence.fasta','a') as file:
        file.write("\n>TranscriptingSequence\n")
        file.write(str(transcribe))
    file.close()
#######

# Fills Z array for given string str[]



######
def stats(seq):
    with open('sequenceStats.txt','w') as file:
        file.write("Sequence Length: " + str(len(seq)) + "\n" + "gc_ratio: " + str((seq.count('G') + seq.count('C'))/len(seq)) + '\n')
    file.close()
    return (seq.count('G') + seq.count('C'))/len(seq)


complementarySequence(a_seq)
print(len(record.seq))
print(stats(record.seq))
print(gc_fraction(record.seq))

###
new_seq = str(Seq("ATGCATGCATGCATGCATGCATG"))
sub_seq = "ATG"
#print(new_seq.search(sub_seq))
runx1 = re.findall("[CGT]TGTGGT[CT][AT]", str(record.seq))
tgif1 = re.finditer("[AT]GACAG[CGT]",str(record.seq))
ikzf1 = re.finditer("[CGT]TGGGA[AG][AGT]",str(record.seq))
y = re.search("[AT]GACAG[CGT]",str(record.seq))
print(tgif1)
print(y.start())
print(ikzf1)
print(runx1)
for match in tgif1:
    print(match.start())
for match in ikzf1:
    print(match)
    print(match.start())
for match in runx1:
    print(match)