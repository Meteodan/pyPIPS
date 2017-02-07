import sys

logdest = sys.stdout

def log(msg):
    logdest.write("%s\n" % msg)
    return

def warning(msg):
    logdest.write("Warning: %s\n" % msg)
    return

def fatal(msg):
    logdest.write("Error: %s\n" % msg)
    sys.exit()


