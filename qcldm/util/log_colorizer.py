import sys, logging

def add_coloring_to_emit_ansi(fn):
    # add methods we need to the class
    def new(*args):
        levelno = args[1].levelno
        if(levelno>=50):
            color = '\x1b[31m' # red
        elif(levelno>=40):
            color = '\x1b[31m' # red
        elif(levelno>=30):
            color = '\x1b[33m' # yellow
        elif(levelno>=20):
            color = '\x1b[32m' # green 
        elif(levelno>=10):
            color = '\x1b[35m' # pink
        else:
            color = '\x1b[0m' # normal
        args[1].msg = color + args[1].msg +  '\x1b[0m'  # normal
        #print "after"
        return fn(*args)
    return new
def colorize():
	logging.StreamHandler.emit = add_coloring_to_emit_ansi(logging.StreamHandler.emit)

def init_log(argv):
	level = logging.WARN
	if '-vv' in argv:
		level = logging.DEBUG
	elif '-v' in argv:
		level = logging.INFO
	logging.basicConfig(format = '%(asctime)s %(levelname)-8s %(message)s', level = level,
                              datefmt='%Y-%m-%d %H:%M:%S')
	colorize()
