#!/usr/bin/python

#Copyright (c) 2017 Brigham Young University

#See the file license.txt for copying permission.

import re
import atom
from copy import deepcopy
from collections import OrderedDict
import functools
import numpy
from numpy import matlib

class MECPFormatException (Exception):
	pass

_blank = lambda line: (line == '' or line.isspace())

_float_pat = r'[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?'
_geom_regex = re.compile(r'\s*(?:' + '|'.join(atom.elements) + r'|\d+)\s+' +
	_float_pat + r'\s+' + _float_pat + r'\s+' + _float_pat + r'\s*')
_float_regex = re.compile(_float_pat)

_alt_geom_header = re.compile(r'[-+]?\d+\s+\d+/\d+')

#_spin_delimiter = re.compile(r"[/,~\-]|\s+")
_sect_regex = re.compile(r'\[(\w+)\]')
_param_delimiter = re.compile(r'\s*[:=]\s*|\s+')

def _parsegeom(section, molecule = None):
	if molecule is None:
		molecule = atom.Molecule()
	print molecule
	for line in section:
		if _blank(line):
			continue
		elif _geom_regex.match(line):
			molecule.append( atom.Atom(*line.split()) )
		else:
			raise MECPFormatException("Invalid Geometry Format: " + line)
	return molecule

def _parsehessian(section, hessian = None):
	if hessian is not None:
		raise MECPFormatException("You may only include one hessian matrix")

	size = int(section.pop(0))
	hessian = matlib.empty((size, size), dtype=numpy.float64)
	r, c = 0, 0
	for line in section:
		for entry in line.split():
			hessian[r,c] = float(entry)
			if c != r:
				hessian[c,r] = hessian[r,c]
			c += 1
			if c > r:
				r += 1
				c = 0
			if r >= size:
				return hessian

	# r < size
	raise MECPFormatException("Not enough values in the Hessian matrix.")

def _parsegradient(section, gradient = None):
	if gradient is not None:
		raise MECPFormatException("You may only include one effective gradient")

	gradient = []
	for line in section:
		gradient += line.split()

	return numpy.array(gradient, numpy.float64)

_multi_route = re.compile(r'(\w+)\((.+)\)')
_single_route = re.compile(r'(\w+)=(\S+)')
_comma_delimiter = re.compile(r',\s*')

def _parsecmd(line):
	result = OrderedDict()

	for item in line.split():
		key, value = None, None
		multi = _multi_route.match(item)
		single = _single_route.match(item)
		if multi:
			key = multi.group(1)
			value = tuple(_comma_delimiter.split(multi.group(2)))
		elif single:
			key = single.group(1)
			value = single.group(2)
		else:
			key = item
			value = None

		result[key] = value

	return result

def _parseroute(section, route = None):
	common = OrderedDict()
	only_a = OrderedDict()
	only_b = OrderedDict()

	if route is None:
		common['force'] = None
		common['integral'] = 'ultrafinegrid'

	for line in section:
		if _blank(line):
			continue

		group = None
		if line.startswith("A:"):
			only_a.update(_parsecmd(line[2:]))
		elif line.startswith("B:"):
			only_b.update(_parsecmd(line[2:]))
		else:
			common.update(_parsecmd(line))

	a, b = None, None
	if route is None:
		a = OrderedDict()
		b = OrderedDict()
	else:
		a, b = route

	a.update(common)
	a.update(only_a)

	b.update(common)
	b.update(only_b)

	return a, b

def _addtoroute(line, route, target=None):
	a, b = route
	if target == 'A':
		a.update(_parsecmd(line))
	elif target == 'B':
		b.update(_parsecmd(line))
	else:
		common = _parsecmd(line)
		a.update(common)
		b.update(common)

	return a, b

def _parsespin(value):
	spins = value.strip().split('/')
	if len(spins) != 2:
		raise MECPFormatException("Invalid format for spin states: " + value)
	try:
		return tuple(map(int, spins))
	except ValueError:
		raise MECPFormatException("Invalid format for spin states: " + value)

def _parsechrgspin(line, params):
	chrg, spin = line.split()
	params['general']['charge'] = int(chrg)
	params['general']['spinstates'] = _parsespin(spin)

def switch(value):
	value = value.strip()
	val = value.lower()
	if val in ['t', 'true', '1', 'on']:
		return True
	elif val in ['f', 'false', '0', 'off']:
		return False
	else:
		raise MECPFormatException("Invalid switch value: " + value)

def enum(*values):
	def ret(value):
		value = value.strip()
		val = value.lower()
		if val in values:
			return val
		else:
			raise MECPFormatException("Invalid value for enum parameter: " + value)
	return ret

def map_type(**mapping):
	def ret(value):
		if value in mapping:
			return mapping[value]
		else:
			raise MECPFormatException("Invalid value for parameter: " + value)
	return ret

def any_type(*types):
	def ret(value):
		for t in types:
			try:
				result = t(value)
				return result
			except:
				continue
		raise MECPFormatException("Invalid value: " + value)
	return ret

def split_type(delimiter, *types):
	def ret(value):
		s = value.split(delimiter)
		if len(s) != len(types):
			raise MECPFormatException("Invalid number of arguments: " + value)
		return tuple(map(lambda x: x[0](x[1]), zip(types, s)))
	return ret

def strip_comments(line):
	index = line.find('!')
	while index >= 0:
		if index >= 1 and line[index - 1] == '\\':
			line = line[:index - 1] + line[index:]
			index = line.find('!', index)
		else:
			return line[:index]
	return line

def linestrip(section, existing = None):

	while len(section) > 0 and _blank(section[0]):
		del section[0]

	while len(section) > 0 and _blank(section[-1]):
		del section[-1]

	result = [s.rstrip() for s in section]

	# tack it on to the existing sections
	if existing is not None:
		result = existing + [''] + result

	return result

##This may deserve its own file
#'''Parses a script section into a dictionary mapping step numbers to runnable data'''
#def parse_script(section, script = None):
	#return []

# Format specification
# Each entry denotes a section, which can either be a dict or a function
# If it's a function, the function tells how to read each line in the section
#  Said function is given an iterator over the lines. It can quit whenever it wants
# That function takes a list of strings (the lines) as its only argument
# If it's a dict, it lists all the available options. If its value is a:
# * callable- its type
# * string- the option is an alias for another parameter
# * tuple of strings- the option is an alias for multiple parameters (puts the same thing in all of them)
# If there is a True key in a section, that indicates the type of parameters not listed
_schema = {
	"cutoffs": {
		'max_grad':     float,
		'rms_grad':     float,
		'max_chg':      float,
		'rms_chg':      float,
		'energy_diff':  float,
		# Aliases
		'maxgrad':      'max_grad',
		'rmsgrad':      'rms_grad',
		'maxchg':       'max_chg',
		'rmschg':       'rms_chg',
		'energydiff':   'energy_diff'
	},
	"link": {
		'nproc':        int,
		'nprocshared':  int,
		'mem':          str,
		True:           str
	},
	"route":    _parseroute,
	"extra":    linestrip,     # footer- extra information
	"geometry": _parsegeom,
	'general': {
		'basis_a':      str,
		'basis_b':      str,
		'max_steps':    int,
		'max_stepsize': float,
		'current_step': int,
		'charge':       int,
		'spinstates':   _parsespin,
		'pre_opt':      split_type(None, any_type(enum('none', 'a', 'b'), int), str),
		'show_hessian': enum('none', 'full', 'sign'),
		'read_later':   switch,  #setting that removes guess keyword if the first step is a preopt
		# Aliases
		'method':       'basis',
		'basis':       ('basis_a', 'basis_b'),
		'basisa':       'basis_a',
		'basisb':       'basis_b',
		'maxsteps':     'max_steps',
		'maxstepsize':  'max_stepsize',
		'currentstep':  'current_step',
		'spin':         'spinstates',
		'preopt':       'pre_opt',
		'readlater':    'read_later'
	},
	# The script section allows you to modify options between steps
	#'script':   parse_script,
	# Data for resuming a run
	'hessian':  _parsehessian,
	'gradient': _parsegradient,
	# Aliases
	'command':  'route',
	"footer":   'extra',
}

# A specification for alternate section starters
# section: (check function, handler for that line, section to start)
# Alternate section starters can affect any of the parameters from any section
_altstart = [
	(lambda line: _alt_geom_header.match(line), _parsechrgspin, 'geometry')
]

# A specification for section enders (defaults to beginning of next section)
# The section immediately switches when the given condition is true
# section: (check function, section to start)
_altend = {
	'geometry':     (_blank, 'extra')
}

# A specification for inlined section contents
# if any line is prefixed with the given string, it will be handled in a certain way.
# prefix:  section
# prefix: (section, handler function)
_inline = [
	('%',   'link'),
	('#',   ('route', _addtoroute)),
	('a#',  ('route', functools.partial(_addtoroute, target='A'))),
	('A#',  ('route', functools.partial(_addtoroute, target='A'))),
	('b#',  ('route', functools.partial(_addtoroute, target='B'))),
	('B#',  ('route', functools.partial(_addtoroute, target='B')))
]

# Default values
_default = {
	"cutoffs": {
		'max_grad':     0.00070,
		'rms_grad':     0.00050,
		'max_chg':      0.00400,
		'rms_chg':      0.00250,
		'energy_diff':  0.00005
	},
	"link": {
	},
	"extra":  [],
	"route":  _parseroute([]),
	"general": {
		'max_steps':    40,
		'max_stepsize': 0.1,
		'current_step':   0,
		'charge':       0,
		'spinstates':   (1,3),
		'basis_a':      'um06/gen',
		'basis_b':      'um06/gen',
		'pre_opt':      ('none', 'm06/gen'),
		'show_hessian': 'none',
		'read_later':   False
	},
	"geometry": None,
	"hessian": None,
	"gradient": None
}

# --========================================================================--

''' verifies the final mecp format after reading the whole file '''
def verify(params):
	# check read_later processed correctly
	if params['general']['current_step'] is 0:
		_checkread(params, 0)
		_checkread(params, 1)
	
	# add guess=mix to singlet if needed
	_addmix(params, 0)
	_addmix(params, 1)

	# verify geometry exists
	if 'geometry' not in params or params['geometry'] is None:
		raise MECPFormatException("Input file is missing geometry")

	# verify basis is handled correctly
	return _verifyBasis(params)

def _verifyUDFT(basis):
	if basis[0] is not 'u':
		raise MECPFormatException("All spin states must use a 'u' DFT method. You are using '"+basis+"'.")

def _verifyBasis(opts):
	# Make sure the basis is only listed once
	routeA, routeB = opts['route']
	basisA = opts['general']['basis_a']
	basisB = opts['general']['basis_b']

	for keyword in routeA:
		if '/' in keyword:
			basisA = keyword
			routeA.pop(keyword, None)

	for keyword in routeB:
		if '/' in keyword:
			basisB = keyword
			routeB.pop(keyword, None)

	_verifyUDFT(basisA)
	_verifyUDFT(basisB)

	opts['route'] = routeA, routeB
	opts['general']['basis_a'] = basisA
	opts['general']['basis_b'] = basisB

	return opts

def _checkread(params, index):
	target = _gettarget(index)
	if 'guess' in params['route'][index]:
		value = params['route'][index]['guess']
		if 'read' in value:
			raise MECPFormatException("To use guess=read you must turn on read_later")

def _addmix(params, index):
	target = _gettarget(index)
	if params['general']['spinstates'][index] is 1:
		_addtoguess(params, index, target, "mix")

def addreadtoboth(params):
	_addtoguess(params, 0, 'A', "read")
	_addtoguess(params, 1, 'B', "read")
	params['general']['read_later'] = False

def _addtoguess(params, index, target, key):
	if 'guess' in params['route'][index]:
		value = params['route'][index]['guess']
		if key in value:
			return
		del params['route'][index]['guess']
		if "," in value:
			_addtoroute("guess=("+key+","+value[1:-1]+")", params['route'], target)
		else:
			_addtoroute("guess=("+key+","+value+")", params['route'], target)
	else:
		_addtoroute("guess="+key, params['route'], target)    

def _gettarget(index):
	if index is 0:
		return 'A'
	elif index is 1:
		return 'B'
	else:
		raise MECPFormatException("Invalid target while adding guess")

# --========================================================================--

'''Reads the input file. Using the format defined above'''
def read(filename):
	params = deepcopy(_default)

	def default_handler(section, line):
		if _blank(line):
			return

		def process_opt(opt, value):

			#Recursively handle aliases
			if opt in _schema[section] and isinstance(_schema[section][opt], tuple):
				for o in _schema[section][opt]:
					process_opt(o, value)

			elif opt in _schema[section] and isinstance(_schema[section][opt], str):
				process_opt(_schema[section][opt], value)

			else:
				func = str
				if opt in _schema[section]:
					func = _schema[section][opt]
				else:
					if True in _schema[section]:
						func = _schema[section][True]

				params[section][opt] = func(value)

		parts = _param_delimiter.split(line, 1)
		opt = parts[0].lower().strip()
		value = parts[1].strip()
		process_opt(opt, value)

	with open(filename) as f:
		section = 'general'
		sectlines = []

		for line in f:
			# Comment stripping
			line = strip_comments(line)

			# Inline handler
			inlined = False
			for prefix, response in _inline:
				if line.startswith(prefix):
					if isinstance(response, str):
						default_handler(response, line[len(prefix):])
					elif isinstance(response, tuple) and len(response) == 2:
						params[response[0]] = response[1](line[len(prefix):], params[response[0]])
					inlined = True
					break
			if inlined:
				continue

			# Look for alternate section starts
			altstart = False
			for check, handler, newsect in _altstart:
				if check(line):
					handler(line, params)
					section = newsect
					altstart = True
					break
			if altstart:
				continue

			# Look for normal section headers
			sectmatch = _sect_regex.match(line)
			if sectmatch:
				if callable(_schema[section]):
					params[section] = _schema[section](sectlines, params[section])
					sectlines = []

				section = sectmatch.group(1).lower()
				if section not in _schema:
					raise MECPFormatException("Invalid section name")
				# Handle Aliases
				while isinstance(_schema[section], str):
					section = _schema[section]

			elif callable(_schema[section]):
				# Handle alternate section enders
				if section in _altend and _altend[section][0](line):
					params[section] = _schema[section](sectlines, params[section])
					sectlines = []
					section = _altend[section][1]
					continue

				sectlines.append(line)

			else:
				default_handler(section, line)

		if callable(_schema[section]):
			params[section] = _schema[section](sectlines)

	return verify(params)

# --========================================================================--

def _writelist(f, data):
	for line in data:
		f.write(line + '\n')

def _writegeom(f, data):
	f.write(str(data) + '\n')

def _writeroute(f, data):
	a, b = data
	for spin, route in zip(['A', 'B'], data):
		for opt, val in route.iteritems():
			if val is None:
				f.write('%s: %s\n' % (spin, opt) )
			elif not isinstance(val, (list, tuple)):
				f.write('%s: %s=%s\n' % (spin, opt, val))
			else:
				f.write('%s: %s(%s)\n' % (spin, opt, ','.join(val)))

def _writehessian(f, data):
	size = data.shape[0]
	f.write('%d\n' % size)
	for row in xrange(size):
		for col in xrange(row + 1):
			f.write(str(data[row, col]))
			if (col + 1) % 5 == 0 or row == col:
				f.write('\n')
			else:
				f.write('\t')
	f.write('\n')

def _writegradient(f, data):
	i = 0
	for val in data:
		f.write(str(val))
		i += 1
		if i >= 3:
			f.write('\n')
			i = 0
		else:
			f.write('\t')

_writeswitch = lambda x: 'on' if x else 'off'

def joinlist(sep, *formatters):
	return lambda data: sep.join(map(lambda x: x[0](x[1]), zip(formatters, data)))

# Format spec for writing a file.
# Each section is either a function that formats it or a dictionary that
# tells how to format each value
# Values are formatted with str() by default
# You can provide a function for complex formats
# Or else, provide a string that can be %-formatted with the value
_writeformat = {
	"cutoffs": {
	},
	"link": {
	},
	"route":    _writeroute,
	"extra":    _writelist,
	"geometry": _writegeom,
	'general': {
		'spinstates': '%d/%d',
		'pre_opt': joinlist(' ', str, str)
	},
	# Used only for resuming an incomplete run
	'hessian':  _writehessian,
	'gradient': _writegradient
}

def _format_std(format, value):
	return format % value

'''Escape all ! on the line.'''
def escape(line):
	return line.replace('!', '\\!')

def write(file, params):
	with open(file, 'w') as f:
		for section, data in params.iteritems():
			if data is not None:
				f.write( escape("[%s]\n" % section) )
				if callable(_writeformat[section]):
					_writeformat[section](f, data)
				else:
					for key, value in data.iteritems():
						formatter = str
						if key in _writeformat[section]:
							formatter = _writeformat[section][key]
							if isinstance(formatter, str):
								formatter = functools.partial(_format_std, formatter)
						f.write( escape("%s = %s\n" % (key, formatter(value))) )

# --====================================================================--

def is_valid_file(parse, arg):
	import os.path
	if not os.path.exists(arg):
		parser.error("The file %s does not exist!" % arg)
	else:
		return arg #return file name

def _parsearguments():
	import argparse
	parser = argparse.ArgumentParser(prog='mecpinput', description='read and write mecp files', add_help=True)
	parser.add_argument("-f", "--filename", dest="filename", type=lambda x: is_valid_file(parser, x), help='name of .mecp file', required=True)
	return parser.parse_args()

# "Unit test"
if __name__ == '__main__':
	from pprint import pprint
	try:
		args = _parsearguments()
		opts = read(args.filename)
		pprint(opts)
		write('testformat-w.mecp', opts)
		pprint(read('testformat-w.mecp'))
		with open('testformat-w.mecp') as f:
			for line in f:
				print line.rstrip()
	except MECPFormatException as e:
		print("  MECPFormatException: " +e.message)
		exit(1)

