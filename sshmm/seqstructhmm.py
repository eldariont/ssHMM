from operator import itemgetter
from collections import Counter
import random
import os
import sys
import math

from MyHMM import Alphabet, MyHMMFromMatrices, DiscreteDistribution, LabelDomain, EmissionSequence, SequenceSet
from pwm import PositionWeightMatrix
from SimpleQueue import SimpleQueue


class SeqStructHMM:

    @classmethod
    def weighted_choice_sub(cls, weights):
        #from: http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/
        rnd = random.random() * sum(weights)
        for i, w in enumerate(weights):
            rnd -= w
            if rnd < 0:
                return i


    @classmethod
    def drawDistribution(cls, probabilities):
        drawnChoice = SeqStructHMM.weighted_choice_sub(probabilities.values())
        return probabilities.keys()[drawnChoice]


    @classmethod
    def drawTopNDistribution(cls, probabilities, n):
        sortedProbs = sorted(probabilities.items(), key=itemgetter(1), reverse=True)
        top10Probs = dict(sortedProbs[:min(n, len(sortedProbs))])
        drawnChoice = SeqStructHMM.weighted_choice_sub(top10Probs.values())
        return top10Probs.keys()[drawnChoice]


    @classmethod
    def drawSimulatedAnnealing(cls, probabilities, temperature):
        normalized_values = [value / float(sum(probabilities.values())) for value in probabilities.values()]
        try:
            annealed_values = map(lambda val: val**(1/float(temperature)), normalized_values)
        except OverflowError:
            print 'Overflow Error in drawSimulatedAnnealing()'
            print 'Probability values:', probabilities.values()
            print 'Temperature:', temperature
            raise
        drawnChoice = SeqStructHMM.weighted_choice_sub(annealed_values)
        return probabilities.keys()[drawnChoice]


    @classmethod
    def getMotifStateTransistions(cls, position, motif_length):
        if position < motif_length:
            trans = [0] + [0] * position * 5 + [1/5.0] * 5 + [0] * (motif_length - position - 1) * 5 + [0]
        elif position == motif_length:
            trans = [0] + [0] * motif_length * 5 + [1]
        else:
            return False
        return trans


    @classmethod
    def generateSequenceStructureHMM(cls, alphabet, motif_length):
        stateNumber = motif_length * 5 + 2
        alphabetSize = len(alphabet) - 2

        #create label list
        labelList = ['St']
        labelList += [item for i in xrange(1, motif_length + 1) for item in ('E'+str(i), 'I'+str(i), 'S'+str(i), 'H'+str(i), 'M'+str(i))]
        labelList += ['En']

        preStateTrans = [0] +  [1.0 / 5] * 5 + [0] * ((motif_length - 1) * 5 + 1)
        postStateTrans = [0] * (motif_length * 5 + 1) + [1]
        transMat = [preStateTrans]
        transMat += [item for i in xrange(1, motif_length + 1) for item in ([SeqStructHMM.getMotifStateTransistions(i, motif_length)] * 5)]
        transMat += [postStateTrans]

        emissionMat = [[0 for i in xrange(alphabetSize)] + [1, 0]]
        emissionMat += [[1.0 / alphabetSize] * alphabetSize + [0, 0] for x in xrange(stateNumber - 2)]
        emissionMat += [[0 for i in xrange(alphabetSize)] + [0, 1]]

        pi = [1] + [0] * (stateNumber - 1)
        m = MyHMMFromMatrices(alphabet, DiscreteDistribution(alphabet), transMat, emissionMat, pi, "stateHMM", labelDomain=LabelDomain(labelList), labelList=labelList)
        return m


    def __init__(self, identifier, jobDirectory, mainLogger, numbersLogger, sequenceContainer, motifLength, flexibility, sample_size, temperature = 1.0, delta = 0.001):
        self.identifier = identifier
        self.job_directory = jobDirectory
        self.model_directory = '{jobDirectory}/model{identifier}/'.format(jobDirectory=self.job_directory, identifier=self.identifier)
        self.main_logger = mainLogger
        self.numbers_logger = numbersLogger
        self.sequence_container = sequenceContainer
        self.motif_length = motifLength
        self.sample_size = sample_size
        self.flexibility = flexibility

        if self.flexibility < 0:
            self.temperature = temperature
            self.delta = delta
        self.iteration = 0

        #Create output directory for this model
        os.mkdir(self.model_directory)
        os.mkdir('{modelDirectory}graphs'.format(modelDirectory=self.model_directory))
        os.mkdir('{modelDirectory}models'.format(modelDirectory=self.model_directory))

        #Prepare alphabets
        extendedAlphabetString = 'ACGT#'
        self.extended_alphabet = Alphabet(extendedAlphabetString)
        extendedAlphabetString2 = 'ACGT+#'
        self.extended_alphabet2 = Alphabet(extendedAlphabetString2)


    def getIteration(self):
        return self.iteration


    def prepare_model_randomly(self):
        self.model = SeqStructHMM.generateSequenceStructureHMM(self.extended_alphabet2, self.motif_length)
        
        #Randomly select motif start positions s1..sm for all m sequences
        for id in self.sequence_container.get_ids():
            self.sequence_container.set_random_motif_starts_for_id(id, self.motif_length)

        self.main_logger.info('Finished setting random motif starts')

    def prepare_model_with_viterbi(self, viterbi_paths):
        self.model = SeqStructHMM.generateSequenceStructureHMM(self.extended_alphabet2, self.motif_length)

        #Viterbi path can tell us where motif starts most likely in each sequence,
        #Use that information to set motif start positions s1..sm for all m sequences
        sequenceIDs = self.sequence_container.get_ids()
        try:
            for i in xrange(len(sequenceIDs)):
                self.sequence_container.set_motif_start_for_id(viterbi_paths[i].index(1), sequenceIDs[i])
        except ValueError:
            print>>sys.stderr, "No motif start could be detected in Viterbi path:", i, viterbi_paths[i]
            print>>sys.stderr, "Training Sequence:", self.trainingSequenceSet[i]
        self.main_logger.info('Finished setting motif starts using viterbi paths')


    def do_iteration(self):
        #Step 2: Randomly take out r sequences
        (sample, rest) = self.sequence_container.sample_k_sequence_ids(self.sample_size)

        #Step 3: Train the HMM M on the remaining m-r sequences
        trainingSequences = self.sequence_container.get_motif_list(rest)
        trainingShapes = self.sequence_container.get_chosen_motif_index_sequences(rest, self.model.internalLabel)
        self.model.estimationTraining(trainingSequences, trainingShapes)

        #Step 4: For each sequence x from the r sequences: select the shape sh and motif start position sr that
        # maximizes P(sh, sr, x | M)
        for sequenceID in sample:
            self.select_shape_and_motif_start(sequenceID)

        if self.flexibility < 0:
            self.temperature -= self.delta
        self.iteration += 1


    def select_shape_and_motif_start(self, sequence_id):
        currentSequence = self.sequence_container.get_sequence_with_id(sequence_id)
        sequenceLength = len(currentSequence)

        probabilities = {}
        multiplicator = None

        for (shapeIndex, (currentShape, shapeWeight, currentMotifStart)) in \
                enumerate(self.sequence_container.get_full_shapes_for_id(sequence_id)):
            #collect probabilities of different shapes and motif start positions
            for motifStart in xrange(0, sequenceLength - self.motif_length + 1):
                #generate emission sequence as list
                emissionSequence = ['+'] + [letter for letter in currentSequence[motifStart:motifStart + self.motif_length]] + ['#']

                #generate state sequence
                stateSequence = [self.model.internalLabel("St")]
                stateSequence += [self.model.internalLabel(currentShape[motifStart + pos] + str(pos + 1))
                                  for pos in xrange(self.motif_length)]
                stateSequence += [self.model.internalLabel("En")]

                joinedProb = self.model.joined(EmissionSequence(self.extended_alphabet2, emissionSequence), stateSequence)
                if (multiplicator == None):
                    multiplicator = -1 * joinedProb
                #log(p) has to be converted to p
                probabilities[(shapeIndex, motifStart)] = math.exp(joinedProb + multiplicator)

        #draw one combination of shape and motif start position according to its probability
        if self.flexibility == 0:
            drawnShape, drawnMotifStart = SeqStructHMM.drawDistribution(probabilities)
        elif self.flexibility > 0:
            drawnShape, drawnMotifStart = SeqStructHMM.drawTopNDistribution(probabilities, self.flexibility)
        else:
            drawnShape, drawnMotifStart = SeqStructHMM.drawSimulatedAnnealing(probabilities, self.temperature)
        """
        if drawnShape != self.sequenceContainer.getChosenShapeForID(sequenceID):
            self.logger.warning("Changed shape from {} to {}".format(self.sequenceContainer.getChosenShapeForID(sequenceID), drawnShape))
        if drawnMotifStart != self.sequenceContainer.getMotifStartForIDAndShape(sequenceID, self.sequenceContainer.getChosenShapeForID(sequenceID)):
            self.logger.warning("Changed motif start from {} to {}".format(self.sequenceContainer.getMotifStartForIDAndShape(sequenceID, self.sequenceContainer.getChosenShapeForID(sequenceID)), drawnMotifStart))
        """
        #store the drawn shape and motif start position
        self.sequence_container.set_chosen_shape_for_id(drawnShape, sequence_id)
        self.sequence_container.set_motif_start_for_shape_and_id(drawnMotifStart, drawnShape, sequence_id)


    def do_training(self, output_frequency, termination_threshold, write_model_state=False):
        """Train the model by repeatedly executing self.do_iteration() until termination.

        Execute self.do_iteration() repeatedly. Every <output_frequency> iterations, print statistics into
        self.main_logger and self.numbers_logger. Every <output_frequency> iterations, also check whether
        the sequence-structure loglikelihood was increased by at least <termination_threshold> compared to any of the 3
        last checks. If it was not, terminate the training and return.

        Keyword arguments:
        output_frequency -- number of iterations between two outputs and termination checks
        termination_threshold -- minimum sequence-structure loglikelihood difference to count as "progress"
        """
        last_models = SimpleQueue()
        while True:
            self.do_iteration()

            if (self.getIteration() % output_frequency == 0) or (self.getIteration() == 1):
                self.output_current_statistics()
                if write_model_state:
                    self.write_current_model()

                current_key_statistics = self.get_key_statistics()
                current_pwm_global = self.get_pwm_global()

                if last_models.length() >= 3:
                    progressMade = False
                    for old_key_statistics, old_pwm_global in last_models.as_list():
                        if (current_key_statistics[2] - old_key_statistics[2]) >= termination_threshold:
                            progressMade = True

                    if not progressMade:
                        self.main_logger.info("Break model %s after %s iterations", 0, self.getIteration())
                        last_model = last_models.as_list()[-1]
                        return last_model

                    if last_models.length() == 3:
                        last_models.pop()
                    elif last_models.length() > 3:
                        raise Exception
                last_models.push((current_key_statistics, current_pwm_global))


    def do_training_for_iterations(self, output_frequency, iterations):
        """Train the model by repeatedly executing self.do_iteration() until termination.

        Execute self.do_iteration() repeatedly. Every <output_frequency> iterations, print statistics into
        self.main_logger and self.numbers_logger. Every <output_frequency> iterations, also check whether
        the sequence-structure loglikelihood was increased by at least <termination_threshold> compared to any of the 3
        last checks. If it was not, terminate the training and return.

        Keyword arguments:
        output_frequency -- number of iterations between two outputs and termination checks
        termination_threshold -- minimum sequence-structure loglikelihood difference to count as "progress"
        """
        while self.getIteration() < iterations:
            self.do_iteration()
            if (self.getIteration() % output_frequency == 0) or (self.getIteration() == 1):
                self.output_current_statistics()

        return (self.get_key_statistics(), self.get_pwm_global())

    def write_current_model(self):

        #Optional: output model parameters
        #self.__printModelParameters()
        #self.__printMotifStarts()
        #self.__printMotifs()

        #output model as graph and xml
        self.model.printAsGraph('{modelDirectory}graphs/graph_{iteration}.png'.format(modelDirectory=self.model_directory, iteration=str(self.iteration).zfill(6)), self.sequence_container.get_length())
        self.model.write('{modelDirectory}models/model_{iteration}.xml'.format(modelDirectory=self.model_directory, iteration=str(self.iteration).zfill(6)))


    def output_current_statistics(self):
        self.main_logger.info('Generate output for iteration %s', self.iteration)

        #output loglikelihoods
        sequence_likelihood = self.calculate_sequence_loglikelihood()
        seq_str_likelihood = self.calculateSequenceStructureLoglikelihood()
        self.main_logger.info('Current sequence loglikelihood: %s', sequence_likelihood)
        self.main_logger.info('Current sequence+structure loglikelihood: %s', seq_str_likelihood)

        #output information contents
        (ic_1000_sequence, ic_1000_structure, ic_1000_combined) = self.calculate_information_contents_1000_sequences()
        ic_global_sequence = self.model.calculateInformationContentSequenceGlobal(self.sequence_container.get_length())
        ic_global_structure = self.model.calculateInformationContentStructure(self.sequence_container.get_length())
        ic_global_combined = self.model.calculateInformationContentCombined(self.sequence_container.get_length())

        self.main_logger.info('Average IC of sequence per position (1000 sequences): %s', sum(ic_1000_sequence) / float(self.motif_length))
        self.main_logger.info('Average IC of structure per position (1000 sequences): %s', sum(ic_1000_structure) / float(self.motif_length))
        self.main_logger.info('Average IC of sequence + structure per position (1000 sequences): %s', sum(ic_1000_combined) / float(self.motif_length))

        self.main_logger.info('Average IC of sequence per position (global weighted average from HMM): %s', sum(ic_global_sequence) / float(self.motif_length))
        self.main_logger.info('Average IC of structure per position (global weighted average from HMM): %s', sum(ic_global_structure) / float(self.motif_length))
        self.main_logger.info('Average IC of sequence + structure per position (global weighted average from HMM): %s', sum(ic_global_combined) / float(self.motif_length))

        #output key statistics into concise logger
        self.numbers_logger.info("{0},{1},{2},{3},{4},{5},{6},{7},{8}".format(self.iteration, sequence_likelihood, seq_str_likelihood,
                                                                  sum(ic_1000_sequence) / float(self.motif_length),
                                                                  sum(ic_1000_structure) / float(self.motif_length),
                                                                  sum(ic_1000_combined) / float(self.motif_length),
                                                                  sum(ic_global_sequence) / float(self.motif_length),
                                                                  sum(ic_global_structure) / float(self.motif_length),
                                                                  sum(ic_global_combined) / float(self.motif_length)))


    def get_key_statistics(self):
        """Return key statistics of the model: iteration nr., loglikelihoods"""
        return [self.iteration, self.calculate_sequence_loglikelihood(), self.calculateSequenceStructureLoglikelihood()]


    def get_pwm_most_likely_path(self):
        mostLikelyPath = self.model.getMostLikelyPath()[1:-1]
        pwm_values = []
        for stateIndex in mostLikelyPath:
            distribution = self.model.getEmission(stateIndex)[0:4]
            pwm_values.append(distribution)
        pwm = PositionWeightMatrix(['A', 'C', 'G', 'T'], pwm_values)
        return pwm


    def get_pwm_global(self):
        probabilityToGetIntoState = self.model._getStateProbabilities()

        pwm_values = []
        for position in xrange((self.model.N - 2 ) / 5):
            probabilityOfBaseAtThisPosition = Counter()
            for stateIndex in xrange(1 + (position * 5), 1 + ((position + 1) * 5)):
                distribution = self.model.getEmission(stateIndex)
                for baseIndex in xrange(4):
                    probabilityOfBaseAtThisPosition[baseIndex] += probabilityToGetIntoState[stateIndex] * distribution[baseIndex]

            pwm_values.append([probabilityOfBaseAtThisPosition[baseIndex] for baseIndex in range(4)])
        pwm = PositionWeightMatrix(['A', 'C', 'G', 'T'], pwm_values)
        return pwm


    def get_pwm_hairpin(self):
        pwm_values = []
        for position in xrange((self.model.N - 2 ) / 5):
            stateIndex = 4 + (position * 5)
            distribution = self.model.getEmission(stateIndex)[0:4]
            pwm_values.append(distribution)
        pwm = PositionWeightMatrix(['A', 'C', 'G', 'T'], pwm_values)
        return pwm


    def get_pwm_best_sequences(self):
        best1000 = self.get_best_n_sequences(1000)

        sequenceCounters = [Counter() for i in range(self.motif_length)]

        #add pseudo-counts
        for position in range(self.motif_length):
            for base in ['A', 'C', 'G', 'T']:
                sequenceCounters[position][base] = 1

        for (sequence, structure) in best1000:
            for position in range(self.motif_length):
                sequenceCounters[position][sequence[position]] += 1

        pwm_values = []
        for position in range(self.motif_length):
            sequenceCounterSum = float(sum(sequenceCounters[position].values()))
            distribution = []

            for symbol in ['A', 'C', 'G', 'T']:
                valueFraction = sequenceCounters[position][symbol] / sequenceCounterSum
                distribution.append(valueFraction)

            pwm_values.append(distribution)

        pwm = PositionWeightMatrix(['A', 'C', 'G', 'T'], pwm_values)
        return pwm


    def __printModelParameters(self):
        file = open(self.model_directory + '_modelParameters.txt', 'w')
        file.write(str(self.model))
        file.close()


    def __printMotifStarts(self):
        file = open(self.model_directory + '_motifStarts.txt', 'w')

        for id in self.sequence_container.get_ids():
            chosenShape = self.sequence_container.get_chosen_shape_for_id(id)
            motifStart = self.sequence_container.get_motif_start_for_id_and_shape(id, chosenShape)
            sequence = self.sequence_container.get_sequence_with_id(id)
            sequenceLength = len(sequence)
            file.write(sequence + "\n")
            file.write("-"*motifStart + "|" * self.sequence_container.get_motif_length() + "-" * (sequenceLength-motifStart-self.sequence_container.get_motif_length()) + "\n")

        file.close()


    def __printMotifs(self):
        file = open(self.model_directory + '_motifs.txt', 'w')

        for id in self.sequence_container.get_ids():
            chosenShape = self.sequence_container.get_chosen_shape_for_id(id)
            motifStart = self.sequence_container.get_motif_start_for_id_and_shape(id, chosenShape)
            sequence = self.sequence_container.get_sequence_with_id(id)
            file.write(sequence[motifStart : motifStart + self.sequence_container.get_motif_length()] + "\n")

        file.close()


    def calculate_sequence_loglikelihood(self):
        """Calculate loglikelihood of all sequences given the model"""
        return self.model.loglikelihood(
            SequenceSet(self.extended_alphabet2, self.sequence_container.get_motif_list(self.sequence_container.get_ids())))


    def calculateSequenceStructureLoglikelihood(self):
        """Calculate joint loglikelihood of all sequences and all (currently chosen) shapes given the model"""
        sequenceIDs = self.sequence_container.get_ids()
        motifs = self.sequence_container.get_motif_list(sequenceIDs)
        indexSequences = self.sequence_container.get_chosen_motif_index_sequences(sequenceIDs, self.model.internalLabel)

        totalJoinedProbs = 0

        for i in xrange(len(sequenceIDs)):
            #generate emission sequence as list
            emissionSequence = [letter for letter in motifs[i]]

            #generate state sequence
            indexSequence = indexSequences[i]

            joinedProb = self.model.joined(EmissionSequence(self.extended_alphabet2, emissionSequence), indexSequence)
            totalJoinedProbs += joinedProb
        return totalJoinedProbs


    def get_best_n_sequences(self, n):
        """Retrieve the n sequences (with corresponding structures) that reach highest likelihood under the model"""
        sequenceIDs = self.sequence_container.get_ids()
        sequenceMotifs = self.sequence_container.get_motif_list(sequenceIDs)
        structureMotifs = self.sequence_container.get_chosen_motif_index_sequences(sequenceIDs, self.model.internalLabel)

        joinedProbabilities = dict()

        for i in xrange(len(sequenceIDs)):
            #generate emission sequence as list
            emissionSequence = [letter for letter in sequenceMotifs[i]]

            #generate state sequence
            indexSequence = structureMotifs[i]

            joinedProb = self.model.joined(EmissionSequence(self.extended_alphabet2, emissionSequence), indexSequence)
            joinedProbabilities[i] = joinedProb

        sortedJoinedProbabilities = sorted(joinedProbabilities.items(), key=itemgetter(1), reverse=True)

        return [(sequenceMotifs[i][1:-1], structureMotifs[i][1:-1]) for (i, joinedProb) in sortedJoinedProbabilities[0:n]]


    def calculate_information_contents_1000_sequences(self):
        """Calculate the information contents based on the 1000 best sequences"""
        best1000 = self.get_best_n_sequences(1000)

        sequenceCounters = [Counter() for i in range(self.motif_length)]
        structureCounters = [Counter() for i in range(self.motif_length)]
        combinedCounters = [Counter() for i in range(self.motif_length)]

        for (sequence, structure) in best1000:
            for position in range(self.motif_length):
                sequenceCounters[position][sequence[position]] += 1
                structureCounters[position][structure[position]] += 1
                combinedCounters[position][sequence[position] + str(structure[position])] += 1

        informationContentSequence = []
        informationContentStructure = []
        informationContentCombined = []
        for position in range(self.motif_length):
            sequenceCounterSum = float(sum(sequenceCounters[position].values()))
            structureCounterSum = float(sum(structureCounters[position].values()))
            combinedCounterSum = float(sum(combinedCounters[position].values()))

            shannonEntropy = 0
            for value in sequenceCounters[position].values():
                valueFraction = value / sequenceCounterSum
                shannonEntropy -= valueFraction * math.log(valueFraction, 2)
            smallSampleCorrection = 3.0 / (math.log1p(2) * 2 * 1000)
            informationContentSequence.append(2 - (shannonEntropy + smallSampleCorrection))

            shannonEntropy = 0
            for value in structureCounters[position].values():
                valueFraction = value / structureCounterSum
                shannonEntropy -= valueFraction * math.log(valueFraction, 2)
            smallSampleCorrection = 4.0 / (math.log1p(2) * 2 * 1000)
            informationContentStructure.append(math.log(5, 2) - (shannonEntropy + smallSampleCorrection))

            shannonEntropy = 0
            for value in combinedCounters[position].values():
                valueFraction = value / combinedCounterSum
                shannonEntropy -= valueFraction * math.log(valueFraction, 2)
            smallSampleCorrection = 19.0 / (math.log1p(2) * 2 * 1000)
            informationContentCombined.append(math.log(20, 2) - (shannonEntropy + smallSampleCorrection))

        return (informationContentSequence, informationContentStructure, informationContentCombined)
