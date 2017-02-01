import time
import getpass
import datetime
import socket


def getRuntime():
	runtime = dict()
	runtime['version']=getVersion()
	runtime['datetime']=curTime
	runtime['hostname']=socket.gethostname()
	runtime['username']=getpass.getuser()
	return runtime


def getVersion():
	version = 'unknown'
	try:
		version = subprocess.check_output(["git", "describe", "--always"]).split()[0]
	except:
		pass
	return version


def printArgs(args):
	argdict=vars(args)
	print 'tool =', str(argdict['func']).split()[1]#[4:]
	for arg in sorted(argdict.keys()):
		if arg != 'func':
			print arg,'=',argdict[arg]
