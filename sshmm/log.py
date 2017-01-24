import logging, sys


def prepareLogger(name, loggerFile, verbose = False, stdout = False):
    newLogger = logging.getLogger(name)

    if stdout:
        stdoutHandler = logging.StreamHandler(sys.stdout)
        if verbose:
            fmt = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
            stdoutHandler.setFormatter(fmt)
        newLogger.addHandler(stdoutHandler)

    fileHandler = logging.FileHandler(loggerFile)
    if verbose:
        fmt = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        fileHandler.setFormatter(fmt)
    newLogger.addHandler(fileHandler)
    # Set the minimal severity of a message to be shown. The levels in
    # increasing severity are: DEBUG, INFO, WARNING, ERROR, CRITICAL
    newLogger.setLevel(logging.INFO)

    return newLogger