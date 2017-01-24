import sys, re
from subprocess import call, Popen, PIPE
from string import upper
import os
import tempfile
import forgi.graph.bulge_graph as fgb

def translateIntoContexts(dotBracketString):
    bg = fgb.BulgeGraph()
    bg.from_dotbracket(dotBracketString)
    rawContextString = upper(bg.to_element_string())
    contextString1 = rawContextString.replace('F', 'E')
    contextString = contextString1.replace('T', 'E')
    return contextString

def calculate_rna_shapes_from_file(output, fastaFileName):
    contextsFile = open(output, 'w')
    print>>sys.stderr, '-->Calculate RNA shapes from', fastaFileName

    proc = Popen(['RNAshapes', '-r', '-o', '1', '-f', fastaFileName], stdout=PIPE, bufsize=-1)
    dotBracketPattern = re.compile('[\(\)\.]+')
    while True:
        data = proc.stdout.readline()
        if not data:
            break
        fields = data.strip().split()
        firstField = fields[0]
        if firstField.startswith('>'):
            print>>contextsFile, firstField
        else:
            match1 = dotBracketPattern.search(firstField)
            if match1 and match1.end() == len(firstField):
                contexts = translateIntoContexts(firstField)
                prob = float(fields[2][1:-1])
                print>>contextsFile, contexts, prob

    print>>sys.stderr, 'Written RNA shapes', output


def calculate_rna_shapes_from_sequence(nucleotide_sequence):
    """Calculates RNAshapes from a given nucleotide sequence
    :param nucleotide_sequence:
    :return: list of (shape, score) tuples
    """
    proc = Popen(['RNAshapes', '-r', '-o', '1', nucleotide_sequence], stdout=PIPE, bufsize=-1)
    dotBracketPattern = re.compile('[\(\)\.]+')
    shapes = []
    while True:
        data = proc.stdout.readline()
        if not data:
            break
        fields = data.strip().split()
        firstField = fields[0]
        match1 = dotBracketPattern.search(firstField)
        if match1 and match1.end() == len(firstField):
            contexts = translateIntoContexts(firstField)
            prob = float(fields[2][1:-1])
            shapes.append((contexts, prob))
    return shapes

def calculate_rna_structures_from_file(output, fastaFileName):
    contextsFile = open(output, 'w')
    print>>sys.stderr, 'Calculate RNA structures from', fastaFileName

    fastaFile = open(fastaFileName, 'r')
    for line in fastaFile:
        if line.startswith('>'):
            sequenceID = line.strip()
        else:
            if sequenceID:
                sequence = line.strip()

                seq_file = tempfile.NamedTemporaryFile('w', delete=False)
                seq_filename = seq_file.name
                print>>seq_file, sequenceID
                print>>seq_file, sequence.upper()
                seq_file.close()

                ct_file = tempfile.NamedTemporaryFile('w', delete=False)
                ct_filename = ct_file.name
                ct_file.close()
                devnull = open(os.devnull, 'w')
                call(['Fold', seq_filename, ct_filename],
                     stdout=devnull)

                print>>contextsFile, sequenceID
                ct_file = open(ct_filename, 'r')
                ct_content = '\n'.join(ct_file.readlines())
                structure_count = ct_content.count('ENERGY')
                ct_file.close()

                dot_file = tempfile.NamedTemporaryFile('w', delete=False)
                dot_filename = dot_file.name
                dot_file.close()
                for i in xrange(structure_count):
                    call(['ct2dot', ct_filename, str(i+1), dot_filename], stdout=devnull)

                    dot_file = open(dot_filename, 'r')
                    dot_file.readline()
                    dot_file.readline()
                    dot_bracket_pattern = dot_file.readline().strip()
                    dot_file.close()
                    contexts = translateIntoContexts(dot_bracket_pattern)
                    print>>contextsFile, contexts, '1'

                sequenceID = None
                os.remove(seq_filename)
                os.remove(ct_filename)
                os.remove(dot_filename)
    print>>sys.stderr, 'Written RNA structures', output


def calculate_rna_structures_from_sequence(nucleotide_sequence):
    """Calculates RNAstructures from a given nucleotide sequence
    :param nucleotide_sequence:
    :return: list of structure strings
    """
    seq_file = tempfile.NamedTemporaryFile('w', delete=False)
    seq_filename = seq_file.name
    print>>seq_file, ">Name"
    print>>seq_file, nucleotide_sequence.upper()
    seq_file.close()

    ct_file = tempfile.NamedTemporaryFile('w', delete=False)
    ct_filename = ct_file.name
    ct_file.close()
    devnull = open(os.devnull, 'w')
    call(['Fold', seq_filename, ct_filename],
         stdout=devnull)

    ct_file = open(ct_filename, 'r')
    ct_content = '\n'.join(ct_file.readlines())
    structure_count = ct_content.count('ENERGY')
    ct_file.close()

    dot_file = tempfile.NamedTemporaryFile('w', delete=False)
    dot_filename = dot_file.name
    dot_file.close()
    structures = []
    for i in xrange(structure_count):
        call(['ct2dot', ct_filename, str(i+1), dot_filename], stdout=devnull)

        dot_file = open(dot_filename, 'r')
        dot_file.readline()
        dot_file.readline()
        dot_bracket_pattern = dot_file.readline().strip()
        contexts = translateIntoContexts(dot_bracket_pattern)
        structures.append(contexts)
        dot_file.close()

    os.remove(seq_filename)
    os.remove(ct_filename)
    os.remove(dot_filename)

    return structures