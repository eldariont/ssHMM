import warnings
from collections import Counter
from ghmm import *
import pygraphviz as pgv
import cEstimate
from pkg_resources import resource_filename

class MyHMMError(GHMMError):
    """Base class for exceptions in this module."""
    def __init__(self, message):
        self.message = message
    def __str__(self):
        return repr(self.message)

class InequalLengthError(MyHMMError):
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)

class EstimationError(MyHMMError):
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)

class MyHMM(StateLabelHMM):
    def __init__(self, emissionDomain, distribution, labelDomain, cmodel):
        StateLabelHMM.__init__(self, emissionDomain, distribution, labelDomain, cmodel)

    def estimationTraining(self, emissionSequenceList, stateSequenceList):
        """ Reestimates the model with the sequences in 'trainingSequences'.

        @note that training for models including silent states is not yet
        supported.

        @param emissionAndStateList list of (emission sequence, state sequence) tuples
        @param nrSteps the maximal number of BW-steps
        @param loglikelihoodCutoff the least relative improvement in likelihood
        with respect to the last iteration required to continue.

        """
        if self.hasFlags(kSilentStates):
            raise NotImplementedError("Sorry, training of models containing silent states not yet supported.")

        estimationResults = cEstimate.cEstimate(emissionSequenceList, stateSequenceList, self.N)
        try:
            (emissionCounts, transitionCounts) = estimationResults
        except TypeError:
            raise EstimationError, 'Estimation results: ' + str(estimationResults)

        #estimate new transition probabilities
        for originStateIndex, counts in enumerate(transitionCounts):
            totalTransitionsFromOrigin = sum(counts)
            for destinationStateIndex, count in enumerate(counts):
                self.setTransition(originStateIndex, destinationStateIndex, float(count) / totalTransitionsFromOrigin)

        #estimate new emission probabilites
        for stateIndex, counts in enumerate(emissionCounts):
            totalEmissionsFromState = sum(counts)
            newDistribution = []
            for count in counts:
                newDistribution.append(float(count) / totalEmissionsFromState)
            self.setEmission(stateIndex, newDistribution)

    #fix bug in implementation of StateLabelHMM
    def externalLabel(self, internal):
        """
        @returns label representation of an int or list of ints
        """

        if type(internal) is int:
            return self.labelDomain.external(internal) # return Label
        elif type(internal) is list:
            return self.labelDomain.externalSequence(internal)
        else:
            raise TypeError('int or list needed')

    def __str__(self):
        hmm = self.cmodel
        strout = [str(self.__class__.__name__)]
        if self.cmodel.name:
            strout.append( " " + str(self.cmodel.name))
        strout.append(  "(N= "+ str(hmm.N))
        strout.append(  ", M= "+ str(hmm.M)+')\n')

        f = lambda x: "%.2f" % (x,) # float rounding function

        if self.hasFlags(kHigherOrderEmissions):
            order = ghmmwrapper.int_array2list(self.cmodel.order, self.N)
        else:
            order = [0]*hmm.N
        label = ghmmwrapper.int_array2list(hmm.label, self.N)

        if hmm.N <= 40:
            iter_list = xrange(self.N)
        else:
            iter_list = [0,1,'X',hmm.N-2,hmm.N-1]

        for k in iter_list:
            if k == 'X':
                strout.append('\n  ...\n\n')
                continue

            state = hmm.getState(k)
            strout.append( "  State "+ str(self.labelDomain.external(label[k])) + ' (')
            if order[k] > 0:
                strout.append( 'order= '+ str(order[k])+',')

            strout.append( "initial= " + f(state.pi)+', id= ' + str(k) + ')\n')
            strout.append( "    Emissions: ")
            for outp in xrange(hmm.M**(order[k]+1)):
                strout.append(f(ghmmwrapper.double_array_getitem(state.b,outp)))
                if outp < hmm.M**(order[k]+1)-1:
                    strout.append( ', ')
                else:
                    strout.append('\n')

            strout.append( "    Transitions:")
            #trans = [0.0] * hmm.N
            for i in xrange( state.out_states):
                strout.append( " ->" + str( state.getOutState(i))+' ('+ f(ghmmwrapper.double_array_getitem(state.out_a,i) ) +')' )
                if i < state.out_states-1:
                    strout.append( ',')
                #strout.append(" with probability " + str(ghmmwrapper.double_array_getitem(state.out_a,i)))

            strout.append('\n')

        return join(strout,'')

    #Calculate probabilities of reaching each state
    def _getStateProbabilities(self):
        stateNames = map(self.externalLabel, range(self.N))
        probabilityToGetIntoState = Counter({0: 1})

        for originIndex, originState in enumerate(stateNames):
            for destinationIndex, destinationState in enumerate(stateNames):
                prob = self.getTransition(self.internalLabel(originState), self.internalLabel(destinationState))
                probabilityToGetIntoState[destinationIndex] += probabilityToGetIntoState[originIndex] * prob

        for key, value in probabilityToGetIntoState.iteritems():
            if value > 1.0:
                probabilityToGetIntoState[key] = 1.0

        return probabilityToGetIntoState

    #Calculate information content of sequence over all possible paths
    def calculateInformationContentSequenceGlobal(self, sequenceNumber):
        probabilityToGetIntoState = self._getStateProbabilities()

        informationContent = []
        for position in xrange((self.N - 2 ) / 5):
            probabilityOfBaseAtThisPosition = Counter()
            for stateIndex in xrange(1 + (position * 5), 1 + ((position + 1) * 5)):
                stateName = self.externalLabel(stateIndex)
                distribution = self.getEmission(self.internalLabel(stateName))
                for baseIndex in xrange(4):
                    probabilityOfBaseAtThisPosition[baseIndex] += probabilityToGetIntoState[stateIndex] * distribution[baseIndex]

            shannonEntropy = 0
            for baseIndex in xrange(4):
                shannonEntropy -= probabilityOfBaseAtThisPosition[baseIndex] * math.log(probabilityOfBaseAtThisPosition[baseIndex], 2)
            smallSampleCorrection = 3.0 / (math.log1p(2) * 2 * sequenceNumber)
            informationContent.append(2 - (shannonEntropy + smallSampleCorrection))

        return informationContent

    #Compile pointer dictionary; for each state s, pointer points to last state in most likely path ending in state s
    def _compileBackwardPointers(self):
        bestPathProbabilities = Counter({0: 1})
        backwardPointers = dict()

        #First position
        for stateIndex in xrange(1, 6):
            bestPathProbabilities[stateIndex] = self.getTransition(0, stateIndex)
            backwardPointers[stateIndex] = 0

        for position in xrange(1, (self.N - 2 ) / 5):
                for currentStateIndex in xrange(1 + (position * 5), 1 + ((position + 1) * 5)):
                    bestProbability = 0
                    for lastStateIndex in xrange(1 + ((position - 1) * 5), 1 + (position * 5)):
                        currentProbability = bestPathProbabilities[lastStateIndex] * self.getTransition(lastStateIndex, currentStateIndex)
                        if currentProbability > bestProbability:
                            bestProbability = currentProbability
                            backwardPointers[currentStateIndex] = lastStateIndex
                    bestPathProbabilities[currentStateIndex] = bestProbability

        #Last state
        bestProbability = 0
        indexOfLastState = self.N - 1
        for lastStateIndex in xrange(self.N - 6, indexOfLastState):
            currentProbability = bestPathProbabilities[lastStateIndex] * self.getTransition(lastStateIndex, indexOfLastState)
            if currentProbability > bestProbability:
                bestProbability = currentProbability
                backwardPointers[indexOfLastState] = lastStateIndex
        bestPathProbabilities[indexOfLastState] = bestProbability

        return backwardPointers

    #Calculate the most likely state path independent from emissions
    def getMostLikelyPath(self):
        backwardPointers = self._compileBackwardPointers()
        currentIndex = self.N - 1
        bestPath = []
        while currentIndex > 0:
            bestPath.append(currentIndex)
            currentIndex = backwardPointers[currentIndex]
        bestPath.append(currentIndex)
        bestPath.reverse()

        return bestPath

    def calculateInformationContentOfMostLikelySequence(self, sequenceNumber):
        mostLikelyPath = self.getMostLikelyPath()[1:-1]

        informationContent = []
        for stateIndex in mostLikelyPath:
            distribution = self.getEmission(stateIndex)

            shannonEntropy = 0
            for baseIndex in xrange(4):
                shannonEntropy -= distribution[baseIndex] * math.log(distribution[baseIndex], 2)
            smallSampleCorrection = 3.0 / (math.log1p(2) * 2 * sequenceNumber)
            informationContent.append(2 - (shannonEntropy + smallSampleCorrection))

        return informationContent

    def calculateInformationContentStructure(self, sequenceNumber):
        probabilityToGetIntoState = self._getStateProbabilities()

        informationContent = []
        for position in xrange((self.N - 2 ) / 5):
            probs = []
            for structureIndex in xrange(5):
                probs.append(probabilityToGetIntoState[1 + (5 * position) + structureIndex])

            shannonEntropy = 0
            for structureIndex in xrange(5):
                scaledProb = probs[structureIndex] / float(sum(probs))
                shannonEntropy -= scaledProb * math.log(scaledProb, 2)
            smallSampleCorrection = 4.0 / (math.log1p(2) * 2 * sequenceNumber)
            informationContent.append(math.log(5, 2) - (shannonEntropy + smallSampleCorrection))

        return informationContent

    def calculateInformationContentCombined(self, sequenceNumber):
        probabilityToGetIntoState = self._getStateProbabilities()

        informationContent = []
        for position in xrange((self.N - 2 ) / 5):
            probs = []
            for stateIndex in xrange(1 + (5 * position), 6 + (5 * position)):
                distribution = self.getEmission(stateIndex)
                for baseIndex in xrange(4):
                    probs.append(probabilityToGetIntoState[stateIndex] * distribution[baseIndex])

            shannonEntropy = 0
            for combinedIndex in xrange(20):
                scaledProb = probs[combinedIndex] / float(sum(probs))
                shannonEntropy -= scaledProb * math.log(scaledProb, 2)
            smallSampleCorrection = 19.0 / (math.log1p(2) * 2 * sequenceNumber)
            informationContent.append(math.log(20, 2) - (shannonEntropy + smallSampleCorrection))

        return informationContent

    def printAsGraph(self, fileName, sequenceNumber):
        graph = pgv.AGraph(strict=False,directed=True)
        stateNames = map(self.externalLabel, range(self.N))
        coordinate_multiplier = 200

        #Add states
        for index, state in enumerate(stateNames):
            if state == 'St':
                html = "<<font point-size=\"50\">Start</font>>"
            elif state == 'En':
                html = "<<font point-size=\"50\">End</font>>"
            else:
                html = None
                distribution = self.getEmission(self.internalLabel(state))
                bases = ["A", "C", "G", "U"]

                fontSizeScaleFactor = 100

                basesInOrderOfProbability = sorted(range(4), key=lambda k: distribution[k], reverse=True)
                shannonEntropy = 0
                for baseIndex in basesInOrderOfProbability:
                    shannonEntropy -= distribution[baseIndex] * math.log(distribution[baseIndex], 2)
                smallSampleCorrection = 3.0 / (math.log1p(2) * 2 * sequenceNumber)
                informationContent = 2 - (shannonEntropy + smallSampleCorrection)

                for baseIndex in basesInOrderOfProbability:
                    currentFontsize = int(round(distribution[baseIndex] * informationContent * fontSizeScaleFactor))
                    currentBase = bases[baseIndex]
                    #newCell = "<tr><td width=\"40\" border=\"0\"><font point-size=\"{0}\" color=\"{1}\">{2}</font></td></tr>".format(currentFontsize, baseColors[currentBase], currentBase)
                    if currentFontsize > 0:
                        if html == None:
                            html = "<<table border=\"5\" color=\"transparent\" cellpadding=\"0\" cellspacing=\"0\">"
                        #Dirty hack because PyGraphViz throws mysterious IOError when height is set to 4 (no explanation given and hours of debugging could not find a reason)
                        if currentFontsize == 4:
                            currentFontsize = 5
                        newCell = "<tr><td border=\"0\" width=\"80\" height=\"{0}\" fixedsize=\"true\"><img scale=\"both\" src=\"{1}\"/></td></tr>".format(currentFontsize, resource_filename('sshmm', 'img/' + currentBase + '.png'))
                        html += newCell

                if html == None:
                    html = "<<font point-size=\"30\">              </font>>"
                else:
                    html += "</table>>"

            #calculate layout position for node
            if index == 0:
                row = 2
                col = 0
            elif index == self.N - 1:
                row = 2
                col = ((index - 1) / 5) + 1
            else:
                row = (index - 1) % 5
                col = ((index - 1) / 5) + 1

            coordinates = "%d,%d!" % (coordinate_multiplier*col, coordinate_multiplier*row)
            graph.add_node(state, label=html, pos=coordinates, shape="box", penwidth = 1)

        #Add transitions
        for originIndex, originState in enumerate(stateNames):
            for destinationIndex, destinationState in enumerate(stateNames):
                prob = self.getTransition(self.internalLabel(originState), self.internalLabel(destinationState))
                if prob > 0.05:
                    graph.add_edge(originState, destinationState, penwidth = prob*10)

        #Add labels
        graph.add_node("Multiloop", label="<<font point-size=\"50\">Multiloop</font>>", pos="{0},{1}!".format(-coordinate_multiplier, 4*coordinate_multiplier), shape="plaintext")
        graph.add_node("Hairpin", label="<<font point-size=\"50\">Hairpin</font>>", pos="{0},{1}!".format(-coordinate_multiplier, 3*coordinate_multiplier), shape="plaintext")
        graph.add_node("Stem", label="<<font point-size=\"50\">Stem</font>>", pos="{0},{1}!".format(-coordinate_multiplier, 2*coordinate_multiplier), shape="plaintext")
        graph.add_node("Internal loop", label="<<font point-size=\"50\">Internal loop</font>>", pos="{0},{1}!".format(-coordinate_multiplier, coordinate_multiplier), shape="plaintext")
        graph.add_node("Exterior loop", label="<<font point-size=\"50\">Exterior loop</font>>", pos="{0},{1}!".format(-coordinate_multiplier, 0), shape="plaintext")

        for position in xrange(1, ((self.N - 2) / 5) + 1):
            graph.add_node(str(position), label="<<font point-size=\"50\">{0}</font>>".format(position), pos="{0},{1}!".format(coordinate_multiplier * position, 4.5*coordinate_multiplier), shape="plaintext")

        #do not print annoying warning of pygraphviz that cell size is too small for content
        warnings.filterwarnings('ignore', '.*cell size too small for content.*', RuntimeWarning)
        graph.layout(prog='neato', args='-n2')
        graph.draw(fileName)


class MyHMMFromMatricesFactory(HMMFactory):
    """ @todo Document matrix formats """

    # XXX TODO: this should use the editing context
    def __call__(self, emissionDomain, distribution, A, B, pi, hmmName = None, labelDomain= None, labelList = None, densities = None):
        log.setLevel(logging.ERROR)

        if isinstance(emissionDomain, Alphabet):

            if not emissionDomain == distribution.alphabet:
                raise TypeError("emissionDomain and distribution must be compatible")

            # checking matrix dimensions and argument validation, only some obvious errors are checked
            if not len(A) == len(A[0]):
                raise InvalidModelParameters("A is not quadratic.")
            if not len(pi) == len(A):
                raise InvalidModelParameters("Length of pi does not match length of A.")
            if not len(A) == len(B):
                raise InvalidModelParameters("Different number of entries in A and B.")

            if (labelDomain is None and labelList is not None) or (labelList is None and labelList is not None):
                raise InvalidModelParameters("Specify either both labelDomain and labelInput or neither.")

            if isinstance(distribution,DiscreteDistribution):
                # HMM has discrete emissions over finite alphabet: DiscreteEmissionHMM
                cmodel = ghmmwrapper.ghmm_dmodel(len(A), len(emissionDomain))

                # assign model identifier (if specified)
                if hmmName != None:
                    cmodel.name = hmmName
                else:
                    cmodel.name = ''

                states = ghmmwrapper.dstate_array_alloc(cmodel.N)
                silent_states = []
                tmpOrder = []

                #initialize states
                for i in xrange(cmodel.N):
                    state = ghmmwrapper.dstate_array_getRef(states, i)
                    # compute state order
                    if cmodel.M > 1:
                        order = math.log(len(B[i]), cmodel.M)-1
                    else:
                        order = len(B[i]) - 1

                    log.debug( "order in state " + str(i) + " = " + str(order) )
                    # check or valid number of emission parameters
                    order = int(order)
                    if  cmodel.M**(order+1) == len(B[i]):
                        tmpOrder.append(order)
                    else:
                        raise InvalidModelParameters("The number of " + str(len(B[i])) +
                                                     " emission parameters for state " +
                                                     str(i) + " is invalid. State order can not be determined.")

                    state.b = ghmmwrapper.list2double_array(B[i])
                    state.pi = pi[i]

                    if sum(B[i]) == 0.0:
                        silent_states.append(1)
                    else:
                        silent_states.append(0)

                    #set out probabilities
                    state.out_states, state.out_id, state.out_a = ghmmhelper.extract_out(A[i])

                    #set "in" probabilities
                    A_col_i = map(lambda x: x[i], A)
                    # Numarray use A[,:i]
                    state.in_states, state.in_id, state.in_a = ghmmhelper.extract_out(A_col_i)
                    #fix probabilities in reestimation, else 0
                    state.fix = 0

                cmodel.s = states
                if sum(silent_states) > 0:
                    cmodel.model_type |= kSilentStates
                    cmodel.silent = ghmmwrapper.list2int_array(silent_states)

                cmodel.maxorder = max(tmpOrder)
                if cmodel.maxorder > 0:
                    log.debug( "Set kHigherOrderEmissions.")
                    cmodel.model_type |= kHigherOrderEmissions
                    cmodel.order = ghmmwrapper.list2int_array(tmpOrder)

                # initialize lookup table for powers of the alphabet size,
                # speeds up models with higher order states
                powLookUp = [1] * (cmodel.maxorder+2)
                for i in xrange(1,len(powLookUp)):
                    powLookUp[i] = powLookUp[i-1] * cmodel.M
                cmodel.pow_lookup = ghmmwrapper.list2int_array(powLookUp)

                # check for state labels
                if labelDomain is not None and labelList is not None:
                    if not isinstance(labelDomain,LabelDomain):
                        raise TypeError("LabelDomain object required.")

                    cmodel.model_type |= kLabeledStates
                    m = MyHMM(emissionDomain, distribution, labelDomain, cmodel)
                    m.setLabels(labelList)
                    return m
                else:
                    return MyHMM(emissionDomain, distribution, labelDomain, cmodel)
            else:
                raise GHMMError(type(distribution), "Not a valid distribution for Alphabet")

        elif isinstance(emissionDomain, Float):
            # determining number of transition classes
            cos = ghmmhelper.classNumber(A)
            if cos == 1:
                A = [A]

            cmodel = ghmmwrapper.ghmm_cmodel(len(A[0]), cos)
            log.debug("cmodel.cos = " + str(cmodel.cos))

            self.constructSwitchingTransitions(cmodel, pi, A)

            if isinstance(distribution, GaussianDistribution):
                #initialize emissions
                for i in xrange(cmodel.N):
                    state = ghmmwrapper.cstate_array_getRef(cmodel.s, i)
                    state.M = 1

                    # set up emission(s), density type is normal
                    emissions = ghmmwrapper.c_emission_array_alloc(1)
                    emission = ghmmwrapper.c_emission_array_getRef(emissions, 0)
                    emission.type = ghmmwrapper.normal
                    emission.dimension = 1
                    (mu, sigma) = B[i]
                    emission.mean.val = mu #mu = mue in GHMM C-lib.
                    emission.variance.val = sigma
                    emission.fixed = 0  # fixing of emission deactivated by default
                    emission.setDensity(0)

                    # append emission to state
                    state.e = emissions
                    state.c = ghmmwrapper.list2double_array([1.0])

                return GaussianEmissionHMM(emissionDomain, distribution, cmodel)

            elif isinstance(distribution, GaussianMixtureDistribution):
                # Interpretation of B matrix for the Gaussian mixture case
                # (Example with three states and two components each):
                #  B = [
                #      [ ["mu11","mu12"],["sig11","sig12"],["w11","w12"]   ],
                #      [  ["mu21","mu22"],["sig21","sig22"],["w21","w22"]  ],
                #      [  ["mu31","mu32"],["sig31","sig32"],["w31","w32"]  ],
                #      ]

                log.debug( "*** mixture model")

                cmodel.M = len(B[0][0])

                #initialize states
                for i in xrange(cmodel.N):
                    state = ghmmwrapper.cstate_array_getRef(cmodel.s, i)
                    state.M = len(B[0][0])

                    # allocate arrays of emmission parameters
                    mu_list = B[i][0]
                    sigma_list = B[i][1]
                    weight_list = B[i][2]

                    state.c = ghmmwrapper.list2double_array(weight_list)

                    # set up emission(s), density type is normal
                    emissions = ghmmwrapper.c_emission_array_alloc(state.M)

                    for j in xrange(state.M):
                        emission = ghmmwrapper.c_emission_array_getRef(emissions, j)
                        emission.type = ghmmwrapper.normal
                        emission.dimension = 1
                        mu = mu_list[j]
                        sigma = sigma_list[j]
                        emission.mean.val = mu #mu = mue in GHMM C-lib.
                        emission.variance.val = sigma
                        emission.fixed = 0  # fixing of emission deactivated by default
                        emission.setDensity(0)

                    # append emissions to state
                    state.e = emissions

                return GaussianMixtureHMM(emissionDomain, distribution, cmodel)

            elif isinstance(distribution, ContinuousMixtureDistribution):
                # Interpretation of B matrix for the continuous mixture case
                # (Example with three states and two components each):
                #  B = [
                #      [["mu11","mu12"], ["sig11","sig12"], ["a11","a12"], ["w11","w12"]],
                #      [["mu21","mu22"], ["sig21","sig22"], ["a21","a22"], ["w21","w22"]],
                #      [["mu31","mu32"], ["sig31","sig32"], ["a31","a32"], ["w31","w32"]],
                #      ]
                #
                # ghmmwrapper.uniform: mu = min, sig = max
                # ghmmwrapper.normal_right or ghmmwrapper.normal_left: a = cutoff

                log.debug( "*** general mixture model")

                cmodel.M = len(B[0][0])

                #initialize states
                for i in xrange(cmodel.N):
                    state = ghmmwrapper.cstate_array_getRef(cmodel.s, i)
                    state.M = len(B[i][0])

                    # set up emission(s), density type is normal
                    emissions = ghmmwrapper.c_emission_array_alloc(state.M)
                    weight_list = B[i][3]

                    combined_map = [(first, B[i][0][n], B[i][1][n], B[i][2][n])
                                    for n, first  in enumerate(densities[i])]

                    for j, parameters in enumerate(combined_map):
                        emission = ghmmwrapper.c_emission_array_getRef(emissions, j)
                        emission.type = densities[i][j]
                        emission.dimension = 1
                        if (emission.type == ghmmwrapper.normal
                            or emission.type == ghmmwrapper.normal_approx):
                            emission.mean.val = parameters[1]
                            emission.variance.val = parameters[2]
                        elif emission.type == ghmmwrapper.normal_right:
                            emission.mean.val = parameters[1]
                            emission.variance.val = parameters[2]
                            emission.min = parameters[3]
                        elif emission.type == ghmmwrapper.normal_left:
                            emission.mean.val = parameters[1]
                            emission.variance.val = parameters[2]
                            emission.max = parameters[3]
                        elif emission.type == ghmmwrapper.uniform:
                            emission.max = parameters[1]
                            emission.min = parameters[2]
                        else:
                            raise TypeError("Unknown Distribution type:" + str(emission.type))

                    # append emissions to state
                    state.e = emissions
                    state.c = ghmmwrapper.list2double_array(weight_list)

                return ContinuousMixtureHMM(emissionDomain, distribution, cmodel)

            elif isinstance(distribution, MultivariateGaussianDistribution):
                log.debug( "*** multivariate gaussian distribution model")

                # this is being extended to also support mixtures of multivariate gaussians
                # Interpretation of B matrix for the multivariate gaussian case
                # (Example with three states and two mixture components with two dimensions):
                #  B = [
                #       [["mu111","mu112"],["sig1111","sig1112","sig1121","sig1122"],
                #        ["mu121","mu122"],["sig1211","sig1212","sig1221","sig1222"],
                #        ["w11","w12"] ],
                #       [["mu211","mu212"],["sig2111","sig2112","sig2121","sig2122"],
                #        ["mu221","mu222"],["sig2211","sig2212","sig2221","sig2222"],
                #        ["w21","w22"] ],
                #       [["mu311","mu312"],["sig3111","sig3112","sig3121","sig3122"],
                #        ["mu321","mu322"],["sig3211","sig3212","sig3221","sig3222"],
                #        ["w31","w32"] ],
                #      ]
                #
                # ["mu311","mu312"] is the mean vector of the two dimensional
                # gaussian in state 3, mixture component 1
                # ["sig1211","sig1212","sig1221","sig1222"] is the covariance
                # matrix of the two dimensional gaussian in state 1, mixture component 2
                # ["w21","w22"] are the weights of the mixture components
                # in state 2
                # For states with only one mixture component, a implicit weight
                # of 1.0 is assumed

                cmodel.addModelTypeFlags(ghmmwrapper.kMultivariate)
                cmodel.dim = len(B[0][0]) # all states must have same dimension

                #initialize states
                for i in xrange(cmodel.N):
                    # set up state parameterss
                    state = ghmmwrapper.cstate_array_getRef(cmodel.s, i)
                    state.M = len(B[i])/2
                    if state.M > cmodel.M:
                        cmodel.M = state.M

                    # multiple mixture components
                    if state.M > 1:
                        state.c = ghmmwrapper.list2double_array(B[i][state.M*2]) # Mixture weights.
                    else:
                        state.c = ghmmwrapper.list2double_array([1.0])

                    # set up emission(s), density type is normal
                    emissions = ghmmwrapper.c_emission_array_alloc(state.M) # M emission components in this state

                    for em in xrange(state.M):
                        emission = ghmmwrapper.c_emission_array_getRef(emissions,em)
                        emission.dimension = len(B[0][0]) # dimension must be same in all states and emissions
                        mu = B[i][em*2]
                        sigma = B[i][em*2+1]
                        emission.mean.vec = ghmmwrapper.list2double_array(mu)
                        emission.variance.mat = ghmmwrapper.list2double_array(sigma)
                        emission.sigmacd = ghmmwrapper.list2double_array(sigma) # just for allocating the space
                        emission.sigmainv = ghmmwrapper.list2double_array(sigma) # just for allocating the space
                        emission.fixed = 0  # fixing of emission deactivated by default
                        emission.setDensity(6)
                        # calculate inverse and determinant of covariance matrix
                        determinant = ghmmwrapper.list2double_array([0.0])
                        ghmmwrapper.ighmm_invert_det(emission.sigmainv, determinant,
                                                     emission.dimension, emission.variance.mat)
                        emission.det = ghmmwrapper.double_array_getitem(determinant, 0)

                    # append emissions to state
                    state.e = emissions

                return MultivariateGaussianMixtureHMM(emissionDomain, distribution, cmodel)

            else:
                raise GHMMError(type(distribution),
                                "Cannot construct model for this domain/distribution combination")
        else:
            raise TypeError("Unknown emission doamin" + str(emissionDomain))

    def constructSwitchingTransitions(self, cmodel, pi, A):
        """ @internal function: creates switching transitions """

        #initialize states
        for i in xrange(cmodel.N):

            state = ghmmwrapper.cstate_array_getRef(cmodel.s, i)
            state.pi = pi[i]

            #set out probabilities
            trans = ghmmhelper.extract_out_cos(A, cmodel.cos, i)
            state.out_states = trans[0]
            state.out_id = trans[1]
            state.out_a = trans[2]

            #set "in" probabilities
            trans = ghmmhelper.extract_in_cos(A,cmodel.cos, i)
            state.in_states = trans[0]
            state.in_id = trans[1]
            state.in_a = trans[2]


MyHMMFromMatrices = MyHMMFromMatricesFactory()
