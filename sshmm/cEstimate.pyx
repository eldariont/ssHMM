from cpython.string cimport PyString_AsString
from cpython.int cimport PyInt_AsLong
from libc.stdlib cimport malloc, free
from cpython.mem cimport PyMem_Malloc, PyMem_Free

def cEstimate(emissionList, stateList, stateNumber):
    cdef int error = 0
    if len(emissionList) != len(stateList):
        error = 1
        return error

    cdef int nrSequences = len(emissionList)
    cdef char **emissionlst = <char **> PyMem_Malloc(nrSequences * sizeof(char*))
    cdef int **statelst = <int **> PyMem_Malloc(nrSequences * sizeof(int*))

    if not emissionlst or not statelst:
        raise MemoryError()

    #convert python variables into C variables
    cdef int nrStates = stateNumber
    cdef int i, j, k, sequencelen
    for i in range(nrSequences):
        emissionlst[i] = PyString_AsString(emissionList[i])
        sequencelen = len(stateList[i])
        if sequencelen != len(emissionList[i]):
            error = 2
            return error
        statelst[i] = <int*> PyMem_Malloc((sequencelen + 1) * sizeof(int))
        if not statelst[i]:
            raise MemoryError
        for j in range(sequencelen):
            statelst[i][j] = <int>PyInt_AsLong(stateList[i][j])
        statelst[i][len(stateList[i])] = -1

    cdef int **emissioncts = <int **> PyMem_Malloc(nrStates * sizeof(int*))
    cdef int **transitioncts = <int **> PyMem_Malloc(nrStates * sizeof(int*))

    if not emissioncts or not transitioncts:
        raise MemoryError()

    cdef int state
    cdef char emission
        
    #add pseudo counts
    for i in range(nrStates):
        emissioncts[i] = <int *> PyMem_Malloc(6 * sizeof(int))
        if not emissioncts[i]:
            raise MemoryError

        if i == 0:
            emissioncts[i][0] = 0
            emissioncts[i][1] = 0
            emissioncts[i][2] = 0
            emissioncts[i][3] = 0
            emissioncts[i][4] = 1
            emissioncts[i][5] = 0
        elif i == nrStates - 1:
            emissioncts[i][0] = 0
            emissioncts[i][1] = 0
            emissioncts[i][2] = 0
            emissioncts[i][3] = 0
            emissioncts[i][4] = 0
            emissioncts[i][5] = 1
        else:
            emissioncts[i][0] = 1
            emissioncts[i][1] = 1
            emissioncts[i][2] = 1
            emissioncts[i][3] = 1
            emissioncts[i][4] = 0
            emissioncts[i][5] = 0

        transitioncts[i] = <int *> PyMem_Malloc(nrStates * sizeof(int))
        if not transitioncts[i]:
            raise MemoryError

        for j in range(nrStates):
            if ((i+4) / 5) + 1 == ((j+4) / 5):
                transitioncts[i][j] = 1
            elif i == j == nrStates - 1:
                transitioncts[i][j] = 1
            else:
                transitioncts[i][j] = 0

    for i in range(nrSequences):
        j = 0
        while emissionlst[i][j] != 0:
            emission = emissionlst[i][j]
            state = statelst[i][j]
            #count emission
            if state < nrStates:
                if emission == 'A':
                    emissioncts[state][0] = emissioncts[state][0] + 1
                elif emission == 'C':
                    emissioncts[state][1] = emissioncts[state][1] + 1
                elif emission == 'G':
                    emissioncts[state][2] = emissioncts[state][2] + 1
                elif emission == 'T':
                    emissioncts[state][3] = emissioncts[state][3] + 1
                elif emission == '+':
                    emissioncts[state][4] = emissioncts[state][4] + 1
                elif emission == '#':
                    emissioncts[state][5] = emissioncts[state][5] + 1
                else:
                    error = 3

            #count transition
            if statelst[i][j+1] != -1:
                transitioncts[state][statelst[i][j+1]] = transitioncts[state][statelst[i][j+1]] + 1
            j = j + 1
        PyMem_Free(statelst[i])
    PyMem_Free(statelst)
    PyMem_Free(emissionlst)

    if error != 0:
        return error

    #convert back to python return type
    emissionCounts = []
    transitionCounts = []
    for i in range(nrStates):
        newEmission = []
        for j in range(6):
            newEmission.append(emissioncts[i][j])
        emissionCounts.append(newEmission)

        newTransition = []
        for j in range(nrStates):
            newTransition.append(transitioncts[i][j])
        transitionCounts.append(newTransition)

        free(emissioncts[i])
        free(transitioncts[i])
    free(emissioncts)
    free(transitioncts)
    return (emissionCounts, transitionCounts)