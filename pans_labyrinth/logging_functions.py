import logging

def create_logger():
    """
    Create the logger for pans_labyrinth
    :return: The root logger for the program
    """

    log = logging.getLogger('pans_labyrinth')
    formatter = logging.Formatter(
        '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    log.setLevel(logging.DEBUG)

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    console.setLevel(logging.INFO)
    log.addHandler(console)

    return log
