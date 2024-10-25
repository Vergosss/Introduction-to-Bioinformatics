import regex as re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction  # gia to pososto GC se sxesh me thn akolouthia

seq = []
records = []
with open('sequence.fasta', 'r') as file:
    # seq = file.read()
    records = parse(file, "fasta")
    for record in records:
        seq = record.seq  # einai sequence

# print(seq)
seq = str(seq)  # thn kano string
sequence = Seq(seq)
# input('...')
print(sequence)
#######
print('Base count of the sequence: ' + str(len(sequence)) + '\nGc percentage: ' + str(100 * gc_fraction(sequence)) + '%')

#######
RUNX1 = re.finditer('[CGT][ACT]TGTGGT[CT][AT]', seq, overlapped=True)
# BHTGTGGTYW
TGIF1 = re.finditer('[AT]GACAG[CGT]', seq, overlapped=True)
# WGACAGB
IKZF1 = re.finditer('[CGT]TGGGA[AG][AGT]', seq, overlapped=True)
# BTGGGARD
# print('RUNX1 : ',RUNX1)
with open('destination.txt', 'w') as destination:
    for match in RUNX1:
        print(match.group(), "start: ", match.start())
        destination.write(str(match.group()) + " ,position: " + str(match.start()) + '\n')

    for match in TGIF1:
        print(match.group(), "start: ", match.start())
        destination.write(str(match.group()) + " ,position: " + str(match.start()) + '\n')

    for match in IKZF1:
        print(match.group(), "start: ", match.start())
        destination.write(str(match.group()) + " ,position: " + str(match.start()) + '\n')
    # ousiastika to start() metra posoi xaraktires yparxoun prin.Diladi h arxh tis akolouthias prin einai 0 den yparxei kanenas
    # xaraktiras prin, eno to 505 px simenei edo pou einai o kersoras exo 505 xaraktires prin kai h zhtoumenh ypoakolouthia einai
    # amesos meta

# print('IKZF1 : ',IKZF1)

###


new_sequences = []
new_sequences.append(SeqRecord(sequence.complement(), name='Complementary Sequence'))
new_sequences.append(SeqRecord(sequence.transcribe(), name='Transcribed Sequence'))

# input('...')

# grafo ta statistika(arithmos vaseon tis akolouthias kai to gc pososto) tis akolouthias se ena ksexoristo arxeio
# epishs grafo tis dyo nees akolouthies(th sybliromatikh kai thn metagrafomenh ) se ksexoristo arxeio
with open('statistics.txt', 'w') as statistics:
    statistics.write('Base count of the sequence: ' + str(len(sequence)) + '\nGc percentage: ' + str(100 * gc_fraction(sequence)) + '%')
with open('complementary-transcription-sequence.fasta', 'w') as sequences:
    SeqIO.write(new_sequences, sequences, 'fasta')

