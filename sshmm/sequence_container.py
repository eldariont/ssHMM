from collections import defaultdict
import random, re

def readSequencesAndShapes(sequencePath, shapePath, motifLength, logger, onlyBestShape = False):
    """Reads sequences and shapes from given FASTA files and returns them as a SequenceContainer object.
    """
    sequenceFile = open(sequencePath, 'r')
    shapeFile = open(shapePath, 'r')

    sequenceLines = sequenceFile.readlines()
    shapeLines = shapeFile.readlines()
    sequenceContainer = SequenceContainer(motifLength, logger)

    uppercasePattern = re.compile("[A-Z]+")

    for i in xrange(len(sequenceLines)/2):
        id = sequenceLines[i*2].split()[0][1:]
        readSequence = sequenceLines[i*2+1].strip()
        sequenceWindow = uppercasePattern.search(readSequence)
        if sequenceWindow:
            newSequence = sequenceWindow.group(0)
        else:
            continue

        valid = True
        for symbol in newSequence:
            if not(symbol in 'ACGT'):
                valid = False

        if valid and len(newSequence) - 2 >= motifLength:
            sequenceContainer.set_sequence_with_id(newSequence, id, sequenceWindow.start(), sequenceWindow.end())
        else:
            logger.warning("Deleted sequence %s because it is either too short or contains symbols other than A, C, G and T.", id)

    currentID = ""
    for line in shapeLines:
        if line[0] == '>':
            #id line
            currentID = line.strip()[1:]
        else:
            if sequenceContainer.contains_id(currentID):
                if onlyBestShape and ('shapes' in sequenceContainer.container[currentID]) and len(sequenceContainer.get_full_shapes_for_id(currentID)) > 0:
                    continue
                (start, stop) = sequenceContainer.get_boundaries_for_id(currentID)
                fields = line.strip().split()
                sequenceContainer.add_shape_with_id(fields[0][start:stop], float(fields[1]), currentID)

    sequenceFile.close()
    shapeFile.close()
    return sequenceContainer


class InvalidMotifBoundaries(Exception):
    """Base class for exceptions in this module."""
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


class NucleotideSequence:
    def __init__(self, seq):
        self.sequence = seq

    def get_sequence(self):
        return self.sequence


class Shape:
    def __init__(self, shape, weight, motif_start):
        self.shape = shape
        self.motifStart = motif_start
        self.weight = weight
        self.uptodate = False
        self.stateSequence = None
        self.indexSequence = None

    def get_shape(self):
        return self.shape

    def get_motif_start(self):
        return self.motifStart

    def set_motif_start(self, motif_start):
        self.motifStart = motif_start
        self.uptodate = False

    def get_weight(self):
        return self.weight

    def update_motif(self, motif_length, index_function):
        if self.motifStart < 0 or (self.motifStart + motif_length) > len(self.shape):
            raise InvalidMotifBoundaries("Invalid arguments motifstart (={0}) or motifLength (={1}) to method "
                                         "updateStateSequences()".format(self.motifStart, motif_length))
        else:
            state_sequence = ["St"]
            state_sequence += [self.shape[self.motifStart + pos] + str(pos + 1) for pos in xrange(motif_length)]
            state_sequence += ["En"]
            self.stateSequence = state_sequence
            self.indexSequence = map(index_function, state_sequence)

    def get_motif_state_sequence(self, motif_length, index_function):
        if not self.uptodate:
            self.update_motif(motif_length, index_function)
            self.uptodate = True
        return self.stateSequence

    def get_motif_index_sequence(self, motif_length, index_function):
        if not self.uptodate:
            self.update_motif(motif_length, index_function)
            self.uptodate = True
        return self.indexSequence


class Shapes:
    def __init__(self):
        self.container = []
        self.chosenShape = 0

    def get_full_shapes(self):
        return [(shapeObj.get_shape(), shapeObj.get_weight(), shapeObj.get_motif_start()) for shapeObj in self.container]

    def get_shapes(self):
        return [shapeObj.get_shape() for shapeObj in self.container]

    def contains_shape(self, shape):
        return shape in self.get_shapes()

    def add_shape(self, shape, weight, motif_start = 1):
        if not(self.contains_shape(shape)):
            new_shape = Shape(shape, weight, motif_start)
            self.container.append(new_shape)

    def set_random_motif_starts(self, motif_length):
        for shapeObj in self.container:
            sequence_length = len(shapeObj.get_shape())
            random_position = random.randint(1, sequence_length - motif_length - 1)
            #initialize motif start for each shape with random position
            shapeObj.set_motif_start(random_position)

    def get_chosen_shape(self):
        #chosenShape must not be set to a non-existent shape
        if self.chosenShape < len(self.container):
            return self.chosenShape
        else:
            raise Exception

    def set_chosen_shape(self, chosen_shape):
        self.chosenShape = chosen_shape

    def get_chosen_motif_state_sequence(self, motif_length, index_function):
        return self.container[self.get_chosen_shape()].get_motif_state_sequence(motif_length, index_function)

    def get_chosen_motif_index_sequence(self, motif_length, index_function):
        return self.container[self.get_chosen_shape()].get_motif_index_sequence(motif_length, index_function)

    def get_chosen_motif_start(self):
        return self.container[self.get_chosen_shape()].get_motif_start()

    def get_motif_start_for_shape_id(self, shape_id):
        return self.container[shape_id].get_motif_start()

    def set_motif_start_for_shape_id(self, motif_start, shape_id):
        self.container[shape_id].set_motif_start(motif_start)

    def set_motif_start_for_all_shapes(self, motif_start):
        for shapeID in xrange(len(self.container)):
            self.container[shapeID].set_motif_start(motif_start)


class SequenceContainer:

    def __init__(self, motif_length, logger):
        self.container = defaultdict(dict)
        self.motifLength = motif_length
        self.logger = logger

    #ID accessors
    def contains_id(self, sequence_id):
        return sequence_id in self.container

    def get_ids(self):
        #return only ids that have sequence and shapes
        return [sequence_id for sequence_id in self.container.iterkeys() if ('seq' in self.container[sequence_id] and 'shapes' in self.container[sequence_id])]

    #Length accessors
    def get_length(self):
        return len(self.get_ids())

    def get_average_sequence_length(self):
        #return only sequences that also have shapes
        sequence_lengths = [len(self.get_sequence_with_id(sequence_id)) for sequence_id in self.container.iterkeys() if ('seq' in self.container[sequence_id] and 'shapes' in self.container[sequence_id])]
        return sum(sequence_lengths) / self.get_length()

    def get_motif_length(self):
        return self.motifLength

    #Motif start Accessors
    def get_motif_start_for_id_and_shape(self, sequence_id, shape_id):
        return self.container[sequence_id]['shapes'].get_motif_start_for_shape_id(shape_id)

    def set_motif_start_for_shape_and_id(self, motif_start, shape_id, sequence_id):
        self.container[sequence_id]['shapes'].set_motif_start_for_shape_id(motif_start, shape_id)

    def set_motif_start_for_id(self, motif_start, sequence_id):
        self.container[sequence_id]['shapes'].set_motif_start_for_all_shapes(motif_start)

    def set_random_motif_starts_for_id(self, sequence_id, motif_length):
        try:
            self.container[sequence_id]['shapes'].set_random_motif_starts(motif_length)
        except KeyError:
            self.logger.warning("No random motif start could be set for %s", sequence_id)

    #Sequence Accessors
    def set_sequence_with_id(self, sequence, sequence_id, start, end):
        self.container[sequence_id]['seq'] = NucleotideSequence(sequence)
        self.container[sequence_id]['core_start'] = start
        self.container[sequence_id]['core_end'] = end

    def get_sequence_with_id(self, sequence_id):
        try:
            return self.container[sequence_id]['seq'].get_sequence()
        except KeyError:
            self.logger.warning("Sequence could not be retrieved for %s", sequence_id)

    def get_boundaries_for_id(self, sequence_id):
        try:
            return self.container[sequence_id]['core_start'], self.container[sequence_id]['core_end']
        except KeyError:
            self.logger.error("Boundaries could not be retrieved for %s", sequence_id)
            raise

    #Shape Accessors
    def add_shape_with_id(self, shape, weight, sequence_id, motif_start = 1):
        try:
            if 'shapes' not in self.container[sequence_id]:
                self.container[sequence_id]['shapes'] = Shapes()
            self.container[sequence_id]['shapes'].add_shape(shape, weight, motif_start)
        except KeyError:
            self.logger.warning("Shape could not be added to %s", sequence_id)

    def get_full_shapes_for_id(self, sequence_id):
        try:
            return self.container[sequence_id]['shapes'].get_full_shapes()
        except KeyError:
            self.logger.warning("Shapes could not be retrieved for %s", sequence_id)

    def get_chosen_shape_for_id(self, sequence_id):
        return self.container[sequence_id]['shapes'].get_chosen_shape()

    def set_chosen_shape_for_id(self, chosen_shape, sequence_id):
        self.container[sequence_id]['shapes'].set_chosen_shape(chosen_shape)

    #List Accessors
    def get_sequences_from_ids(self, id_list):
        return [self.container[sequence_id]['seq'].get_sequence() for sequence_id in id_list]

    def get_motif_list(self, id_list):
        results = []
        for sequence_id in id_list:
            chosen_motif_start = self.container[sequence_id]['shapes'].get_chosen_motif_start()
            results.append('+' + self.container[sequence_id]['seq'].get_sequence()[chosen_motif_start:chosen_motif_start + self.get_motif_length()] + '#')
        return results

    #not in use
    #def get_chosen_motif_state_sequences(self, idList, indexFunction):
    #    return [self.container[id]['shapes'].get_chosen_motif_state_sequence(self.get_motif_length(), indexFunction) for id in idList]

    def get_chosen_motif_index_sequences(self, id_list, index_function):
        return [self.container[sequence_id]['shapes'].get_chosen_motif_index_sequence(self.get_motif_length(), index_function) for sequence_id in id_list]

    #not in use
    #def getSequenceAndStateSequenceList(self, idList, indexFunction):
    #    return [(self.container[id]['seq'].get_sequence(), self.container[id]['shapes'].get_chosen_motif_state_sequence(self.get_motif_length(), indexFunction)) for id in idList]

    #Random
    def sample_k_sequence_ids(self, k):
        all_sequence_ids = self.get_ids()
        sampled_sequence_ids = random.sample(all_sequence_ids, k)
        rest = [sequence_id for sequence_id in all_sequence_ids if not(sequence_id in sampled_sequence_ids)]
        return sampled_sequence_ids, rest

    """
    #Entire object accessors
    def addObjectWithID(self, object, id):
        if self.contains_id(id):
            raise Exception
        else:
            self.container[id] = object

    def moveObjectWithIDToOtherSequenceContainer(self, id, otherSequenceContainer):
        otherSequenceContainer.addObjectWithID(self.container[id], id)
        self.removeObjectWithID(id)

    def removeObjectWithID(self, id):
        try:
            del self.container[id]
        except KeyError:
            self.logger.warning("ID %s could not be removed", id)"""
