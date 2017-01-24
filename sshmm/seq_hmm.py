from MyHMM import DiscreteDistribution, HMMFromMatrices, Alphabet, SequenceSet, log

import random
import logging
import itertools

class HMMOption:
    def __init__(self):
        pass

    def setArgumentToValue(self, argument, value):
        setattr(self, argument, value)

def get_trans_for_state(state, stateNumber):
    trans = [0] * stateNumber
    trans[state + 1] = 1
    return trans

def generate_sequence_hmm(averageSequenceLength, alphabet, options):
    """
    Creates, initializes and returns a Hidden Markov Model based on the length of the motif that shall be captured,
    the average length of the training sequences and the alphabet of the training sequences. The motif is modelled
    as a simple chain of states without loops. Before the motif states comes a looping pre state, after the motif a
    looping post state and finally an end state.
    The returned model is not yet trained but can be trained with the '.baumWelch(<trainingSequences>)' method.'
    :param motifLength: length of the motif to be captured by the model
    :param averageSequenceLength: average length of training sequences
    :param alphabet: alphabet of the HMM
    :param options: parameters parsed by argparse
    :param options.force_end: shall the model be designed so that it must finish in the end state?
    :param options.randomize: shall the initial emission probabilities be randomized?
    :return: initialized HMM of class ghmm.HMM
    """
    loopProb = 1.001 - (2.0 / (averageSequenceLength - options.motif_length))

    if options.force_end:
        stateNumber = options.motif_length + 3
        #don't count end symbol
        alphabetSize = len(alphabet) - 1

        preStateTrans = [loopProb, 1 - loopProb] + [0] * (options.motif_length + 1)
        postStateTrans = [0] * (options.motif_length + 1) + [loopProb, 1 - loopProb]
        endStateTrans = [0] * (options.motif_length + 3)
        transMat = [preStateTrans] + [get_trans_for_state(i, stateNumber) for i in range(1, options.motif_length + 1)] + [postStateTrans, endStateTrans]

        if options.randomize:
            #generate random emission matrix
            emissionMat = []
            for state in range(stateNumber):
                currentEmissions = [random.random() for i in range(alphabetSize)]
                emissionSum = sum(currentEmissions)
                emissionMat.append([emission / emissionSum for emission in currentEmissions] + [0])
            #end state emits end symbol with probability 1
            emissionMat[-1] = [0] * (alphabetSize) + [1]
        else:
            emissionMat = [([1.0 / alphabetSize] * alphabetSize) + [0] for x in range(stateNumber)]
            #end state emits end symbol with probability 1
            emissionMat[-1] = [0] * (alphabetSize) + [1]
    else:
        stateNumber = options.motif_length + 2
        alphabetSize = len(alphabet)

        preStateTrans = [loopProb, 1 - loopProb] + [0] * options.motif_length
        postStateTrans = [0] * (options.motif_length + 1) + [1]
        transMat = [preStateTrans] + [get_trans_for_state(i, stateNumber) for i in range(1, options.motif_length + 1)] + [postStateTrans]

        if options.randomize:
            #generate random emission matrix
            emissionMat = []
            for state in range(stateNumber):
                currentEmissions = [random.random() for i in range(alphabetSize)]
                emissionSum = sum(currentEmissions)
                emissionMat.append([emission / emissionSum for emission in currentEmissions])
        else:
            emissionMat = [([1.0 / alphabetSize] * alphabetSize) for x in range(stateNumber)]

    pi = [1] + [0] * (stateNumber - 1)
    m = HMMFromMatrices(alphabet, DiscreteDistribution(alphabet), transMat, emissionMat, pi)

    return m


def find_best_baumwelch_sequence_models(motif_length, sequence_container, logger, num_models = 100):
        #Train pure sequence HMM with Baum-Welch algorithm and training sequences
        logger.info('Analyzing sequence motif..')
        #suppress annoying warnings by ghmm
        log.setLevel(logging.ERROR)

        options = HMMOption()
        options.setArgumentToValue('motif_length', motif_length)
        options.setArgumentToValue('force_end', True)
        options.setArgumentToValue('randomize', True)

        extendedAlphabetString = 'ACGT#'
        extendedAlphabet = Alphabet(extendedAlphabetString)

        trainingSequences = sequence_container.get_sequences_from_ids(sequence_container.get_ids())
        trainingSequencesWithEndSymbol = [sequence + '#' for sequence in trainingSequences]
        trainingSequenceSet = SequenceSet(extendedAlphabet,trainingSequencesWithEndSymbol)

        best_loglikelihood = None
        best_viterbi_paths = None
        for index in xrange(num_models):
            cur_sequence_model = generate_sequence_hmm(sequence_container.get_average_sequence_length(), extendedAlphabet, options)
            cur_sequence_model.baumWelch(trainingSequenceSet)
            cur_loglikelihood = cur_sequence_model.loglikelihood(trainingSequenceSet)
            if best_loglikelihood == None or cur_loglikelihood > best_loglikelihood:
                best_viterbi_paths = cur_sequence_model.viterbi(trainingSequenceSet)[0]
                best_loglikelihood = cur_loglikelihood
        logger.info('Found sequence motif.')
        return (best_loglikelihood, best_viterbi_paths)
