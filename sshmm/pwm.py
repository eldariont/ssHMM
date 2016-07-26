import math
import random
import numpy as np
from corebio.seq import Alphabet
from weblogolib import LogoData, LogoFormat, LogoOptions, classic, png_print_formatter

"""Randomly draws an integer from the interval [0,i) based on a given list of i weights.
"""
def weighted_random(weights):
#from: http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/
    rnd = random.random() * sum(weights)
    for index, weight in enumerate(weights):
        rnd -= weight
        if rnd < 0:
            return index

class PositionWeightMatrix:
    """A sample implementation of a position weight matrix (PWM)"""
    def __init__(self, alphabet = None, values = None):
        self.alphabet = alphabet
        self.values = values


    def is_initialized(self):
        if self.alphabet != None and self.values != None:
            return True
        else:
            return False

    def get_average_information_content(self):
        if self.is_initialized():
            information_content = []
            for position in xrange(len(self.values)):
                shannon_entropy = 0
                for symbolIndex in xrange(len(self.alphabet)):
                    shannon_entropy += self.values[position][symbolIndex] * math.log(self.values[position][symbolIndex], 2)
                information_content.append(math.log(len(self.alphabet), 2) + shannon_entropy)
            average_IC = sum(information_content) / len(self.values)
            return average_IC

    """Draw a nucleotide motif from the PWM.
    """
    def draw_motif(self):
        if self.is_initialized():
            sequence = ""
            for distribution in self.values:
                drawn_index = weighted_random(distribution)
                drawn_base = self.alphabet[drawn_index]
                sequence += drawn_base
            return sequence


    """Print the PWM into a given file.
    """
    def write_to_file(self, filepath):
        if self.is_initialized():
            pwm_file = open(filepath, 'w')
            print>>pwm_file, "\t".join(self.alphabet)
            for distribution in self.values:
                print>>pwm_file, "\t".join(map(str, distribution))
            pwm_file.close()


    """Reads a PWM from a file
    """
    def read_from_file(self, filepath):
        pwm_file = open(filepath, 'r')
        alphabet_line = pwm_file.readline().strip()
        self.alphabet = alphabet_line.split('\t')

        pwm_values = []
        for line in pwm_file:
            fields = line.strip().split('\t')
            if len(fields) == len(self.alphabet):
                pwm_values.append(map(float, fields))
        self.values = pwm_values

        pwm_file.close()


    """Writes a weblogo of the PWM.
    """
    def write_weblogo(self, filepath):
        matrix_tuple = []
        for distribution in self.values:
            matrix_tuple.append(tuple(distribution))

        dataArray = np.array( tuple(matrix_tuple) )

        alph = Alphabet(''.join(self.alphabet))

        weblogoData = LogoData.from_counts(alph, dataArray)
        weblogoOptions = LogoOptions(color_scheme=classic)
        weblogoOptions.title = "PWM"
        weblogoFormat = LogoFormat(weblogoData, weblogoOptions)
        weblogo_file = open(filepath, 'w')
        weblogo_file.write(png_print_formatter(weblogoData, weblogoFormat))
        weblogo_file.close()

    def calculate_KL_divergence_to_recovered_pwm(self, pwm):
        if self.alphabet != pwm.alphabet:
            raise Exception

        divergence = 0
        for position, distribution in enumerate(self.values):
            for symbol_index in xrange(len(self.alphabet)):
                try:
                    divergence += distribution[symbol_index] * math.log(distribution[symbol_index] /
                                                                        float(pwm.values[position][symbol_index]), 2)
                except IndexError:
                    print "IndexError raised while calculating KL-divergence"
                    print "self pwm:", self.values
                    print "other pwm:", pwm.values
                    raise
                except ZeroDivisionError:
                    print "ZeroDivisionError raised while calculating KL-divergence"
                    print "self pwm:", self.values
                    print "other pwm:", pwm.values
                    raise
        return divergence

    def calculate_KL_divergences_to_set(self, pwms):
        return [self.calculate_KL_divergence_to_recovered_pwm(pwm) for pwm in pwms]


