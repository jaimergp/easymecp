#!/apps/python/2.7.5/bin/python

# Copyright (c) 2017 Brigham Young University

# See the file license.txt for copying permission.

# Renamed mecpopt and placed in the user's bin directory
# Handles the mecp calculations at each step

import sys

import atom
import re, math

# Needed to interface with Gaussian
import gaussian
import subprocess
import signal

# Efficiency imports
from itertools import *
from copy import deepcopy

# Matrices and Vectors
# I trust it for efficiency and accuracy and it IS available on FSL
import numpy
from numpy import linalg, matlib

# Reporting tools
# import cProfile
from datetime import *
_g09_time = timedelta()

_bohr = 0.529177

class MECPException(Exception):
	pass

import mecpinput

def read_prog(filename):
	try:
		opts = mecpinput.read(filename)
	except mecpinput.MECPFormatException as e:
		raise e
	except Exception as e:
		raise MECPException("Error when reading input file:\n" + str(e))

	molecule = opts['geometry']
	del opts['geometry']

	# Do some processing to get the data ready for sending to gaussian

	# Resuming options - Assert validity
	if opts['hessian'] is not None:
		if opts['gradient'] is None:
			raise MECPException("Input file includes a hessian but no gradient.")
		elif opts['hessian'].shape[0] != opts['gradient'].size:
			raise MECPException("Input file's hessian and gradient that are inconsistent with each other.")
		elif opts['gradient'].size != len(molecule) * 3:
			raise MECPException("The gradient and hessian are inconsistent with the molecule.")
	elif opts['gradient'] is not None:
		raise MECPException("Input file includes a gradient but no hessian.")

	return molecule, opts

_g09_force_pattern = re.compile(r'\s*Center\s+Atomic\s+Forces')

'''Extracts energies and gradients from a gaussian logfile.
Based on extract_energy and extract_gradient awk scripts'''
def extract_data(filename, natoms):
	energy = None
	gradient = []

	with open(filename) as logfile:
		f = iter(logfile)
		for line in f:
			line = line.strip()
			if 'Convergence failure' in line:
				raise MECPException("Convergence Error")
			elif line.startswith('SCF Done:'):
				energy = float(line.split()[4])
			elif _g09_force_pattern.match(line):
				try:
					f.next()
					f.next()
					for i in xrange(natoms):
						fields = f.next().split()
						gradient.extend( map(float, fields[2:5]) )
				except StopIteration:
					raise MECPException("Failed to extract forces (end of file)")

	if energy is None or len(gradient) != natoms * 3:
		raise MECPException("ERROR: Gaussian run produced an invalid logfile. There is missing or extra data.")

	# Oh, and convert from Hartree/Bohr to Hartree/Ang while we're at it
	return energy, -numpy.array(gradient, numpy.float64) / _bohr

'''Extracts geometry from a gaussian optimization job.'''
def extract_geometry(filename):
	molecule = atom.Molecule()

	with open(filename) as logfile:
		f = iter(logfile)
		for line in f:
			if "Standard orientation:" in line:
				try:
					#skip 4 lines
					[f.next() for x in range(4)]
					for line in f:
						if '-----' in line:
							return molecule
						fields = line.split()
						molecule.append(atom.Atom(fields[1], *fields[3:]))
				except StopIteration:
					break

	raise MECPException("ERROR: Gaussian produced an invalid logfile. " +
		"Could not find atomic coordinates")

dict_union = lambda x, y: dict(x.viewitems() | y.viewitems())

'''Calls gaussian to do the heavy lifting, then extracts energies and gradients from the resulting log file.
Based on sub_script'''
def call_gaussian(job, step, options, molecule):
	prefix = job + str(step)

	genopts = options['general']
	spinstates = genopts['spinstates']
	link = options['link']

	global _g09_time

	status, time =  gaussian.run_g09(prefix + '_A',
					"Step %d of MECP optimization - First Spin State (%d)" % (step, spinstates[0]),
					genopts['charge'], spinstates[0], molecule, 'n', genopts['basis_a'],
					options['route'][0], options['extra'],
					dict_union(link, {'chk': job + '_A'}))

	if status != 0:
		raise MECPException("Gaussian '09 exited with non-zero status.")

	_g09_time += time

	energyA, gradientA = extract_data(prefix + "_A.log", len(molecule))

	status, time =  gaussian.run_g09(prefix + '_B',
					"Step %d of MECP optimization - Second Spin State (%d)" % (step, spinstates[1]),
					genopts['charge'], spinstates[1], molecule, 'n', genopts['basis_b'],
					options['route'][1], options['extra'],
					dict_union(link, {'chk': job + '_B'}))

	if status != 0:
		raise MECPException("Gaussian '09 exited with non-zero status.")

	_g09_time += time

	energyB, gradientB = extract_data(prefix + "_B.log", len(molecule))

	return (energyA, energyB), (gradientA, gradientB)

'''Optimize a molecule's geometry using Gaussian'''
def gaussian_opt(job, options, molecule):
	jobname = job + "-opt"

	genopts = options['general']

	route = {}

	spin = 1
	if genopts['pre_opt'][0] == 'a':
		spin = genopts['spinstates'][0]
		route = deepcopy(options['route'][0])
	elif genopts['pre_opt'][0] == 'b':
		spin = genopts['spinstates'][1]
		route = deepcopy(options['route'][1])
	else: #pre_opt is an int
		spin = genopts['pre_opt'][0]
		if spin == genopts['spinstates'][1]:
			route = deepcopy(options['route'][1])
		else:
			route = deepcopy(options['route'][0])

	basis = genopts['pre_opt'][1]

	# Previously: if 'guess' in route and route['guess'] == 'mix':
	if 'guess' in route: 
		del route['guess']

	del route['force']
	route['opt'] = None

	global _g09_time

	status, time = gaussian.run_g09(jobname, "Pre-Optimization of " + job,
									genopts['charge'], spin,
									molecule, 'n', basis, route,
									options['extra'], options['link'])

	if status != 0:
		raise MECPException("Gaussian '09 exited with non-zero status.")

	_g09_time += time

	return extract_geometry(jobname + '.log')

def rms(iterable):
	return math.sqrt(numpy.mean(numpy.square(iterable)))

'''Calculates the gradients needed for optimization.
Calculations based on the paper in section 3.'''
def effective_gradient(energies, gradients):
	# Get the energy difference:
	# (E1 - E2)
	energy_delta = energies[0] - energies[1]
	g1, g2 = gradients

	# Get the difference gradient- i.e. the difference of the derivatives of energy with respect to position
	# The derivatives are obtained via force calculations in Gaussian
	# This corresponds to part of equation 2:
	# x1 = (dE1/dq - dE2/dq)
	diff = g1 - g2

	# Parallel Gradient
	# Get the parallel gradient via the Graham-Schmidt process
	normdiff = linalg.norm(diff)
	par = g1 - (numpy.inner(g1, diff) / normdiff) * (diff / normdiff)

	# Perpendicular Gradient
	# (E1 - E2) * x1
	perp = energy_delta * diff

	# Effective (Composite) gradient
	# This is a linear combination of the perpendicular and parallel gradients
	effective = 140.0 * perp + 1.0 * par
	# These are the same coefficients as the original program
	# POSSIBILITY: customizing these coefficients? Would that help?

	# Be sure that par and diff are orthogonal
	# dot = numpy.inner(diff, par)
	# print "f and g are " + ("" if abs(dot) < 0.0001 else "not ") + "orthogonal (" + str(dot) + ')'

	return diff, par, effective

'''Translated loosely from UpdateX.f, with help from Wikipedia for BFGS searches'''
def update_coords(step, stpmax, molecules, gradients, hessian_in):

	c1 = numpy.array(molecules[0].components, numpy.float64)
	c2 = numpy.array(molecules[1].components, numpy.float64)
	g1, g2 = gradients

	stpmax_loose = stpmax * len(g2)

	#Short-circuit for 0th step
	if step == 0:
		hessian_out = numpy.matrix(hessian_in)
	else:
		# get the deltas of the gradient and coordinates
		dgrad   = g2 - g1 # y[k]
		dcoords = c2 - c1 # s[k]: will be zero on first step.

		# calculate y[k].T * B[k].I * y[k] and an associative half
		lcollapse = (numpy.asmatrix(dgrad) * hessian_in).A.flatten()   # y[k].T * B[k].I
		# We don't need B[k].I * y[k] because it's the same as y[k].T * B[k].I
		# reasoning: A * B = (B.T * A.T).T, the hessian is symmetric, and the
		# representation is flat (rather than column/row vectors)
		collapse = numpy.inner(dgrad, lcollapse)

		# s[k].T * y[k]
		denom = numpy.inner(dgrad, dcoords)

		if denom == 0.0:
			raise MECPException("Could not update coordinates: no difference between forces")

		# Math based on Wikipedia article on BFGS
		u = (denom + collapse) * numpy.outer(dcoords, dcoords) / (denom * denom)
		v = (numpy.outer(lcollapse, dcoords) + numpy.outer(dcoords, lcollapse)) / denom

		hessian_out = hessian_in + u - v

	vector = -(numpy.matrix(g2) * hessian_out).A.flatten()

	stpl = linalg.norm(vector)

	if stpl > stpmax_loose:
		vector = vector / stpl * stpmax_loose

	max_comp = numpy.fabs(vector).max()

	if max_comp > stpmax:
		vector = vector / max_comp * stpmax

	coords_out = c2 + vector

	molecule_out = deepcopy(molecules[1])
	molecule_out.components = coords_out

	#print vector

	return molecule_out, hessian_out

'''Tells you whether there was a convergence or not
returns a bool'''
def test_convergence(logfile, molecules, energies, eff_grad, cutoffs):

	# Do some calculations...
	dcoords = numpy.array(molecules[0].components) - numpy.array(molecules[1].components)

	max_delta = max(numpy.absolute(dcoords))
	max_grad  = max(numpy.absolute(eff_grad))

	dcoords_rms  = rms(dcoords)
	eff_grad_rms = rms(eff_grad)

	energy_diff = abs(energies[1] - energies[0])

	# Test for convergence
	conv = [
		max_grad         < cutoffs['max_grad'],
		eff_grad_rms     < cutoffs['rms_grad'],
		max_delta        < cutoffs['max_chg'],
		dcoords_rms      < cutoffs['rms_chg'],
		energy_diff      < cutoffs['energy_diff']
	]

	yn = {True: 'Yes', False:' No'}

	logfile.write("Convergence Check:   actual  (threshold) status\n")
	logfile.write("Max Gradient El.:  % .6f (% .6f) %s\n" % (max_grad,         cutoffs['max_grad'],    yn[conv[0]]))
	logfile.write("RMS Gradient El.:  % .6f (% .6f) %s\n" % (eff_grad_rms,     cutoffs['rms_grad'],    yn[conv[1]]))
	logfile.write("Max Displacement:  % .6f (% .6f) %s\n" % (max_delta,        cutoffs['max_chg'],     yn[conv[2]]))
	logfile.write("RMS Displacement:  % .6f (% .6f) %s\n" % (dcoords_rms,      cutoffs['rms_chg'],     yn[conv[3]]))
	logfile.write("Difference in E:   % .6f (% .6f) %s\n" % (energy_diff,      cutoffs['energy_diff'], yn[conv[4]]))
	logfile.write('\n')

	return all(conv)

# Some functions for working with the logfile

'''Pretty-prints the gradient to the logfile'''
def write_gradient(logfile, gradient):
	for i in xrange(0, len(gradient), 3):
		logfile.write("\t% .8f \t % .8f \t % .8f \n" % tuple(gradient[i:i+3]))

'''Pretty-prints the hessian to the logfile'''
def write_hessian(logfile, hessian, block = 6):
	size = hessian.shape[0]
	for diag in xrange(0, size, block):
		logfile.write('    ')
		for c in xrange(diag, min(diag + block, size)):
			logfile.write(('XYZ'[c % 3] + str(c // 3 + 1)).center(12))
		logfile.write('\n')
		for row in xrange(diag, size):
			logfile.write(('XYZ'[row % 3] + str(row // 3 + 1)).ljust(4))
			for col in xrange(diag, min(row + 1, diag + block)):
				logfile.write(('% .8f' % hessian[row, col]).rjust(12))
			logfile.write('\n')
	logfile.write('\n')

'''Pretty-prints the hessian to the logfile'''
def write_hessian_sign(logfile, hessian, block = 18):
	size = hessian.shape[0]
	for diag in xrange(0, size, block):
		logfile.write('    ')
		for c in xrange(diag, min(diag + block, size)):
			logfile.write(('XYZ'[c % 3] + str(c // 3 + 1)).center(4))
		logfile.write('\n')
		for row in xrange(diag, size):
			logfile.write(('XYZ'[row % 3] + str(row // 3 + 1)).ljust(4))
			for col in xrange(diag, min(row + 1, diag + block)):
				val = hessian[row, col]
				if val == 0:
					logfile.write('  0 ')
				elif val < -1:
					logfile.write(' ---')
				elif val > 1:
					logfile.write(' +++')
				elif val < -0.35:
					logfile.write(' -- ')
				elif val > 0.35:
					logfile.write(' ++ ')
				elif val < 0:
					logfile.write('  - ')
				elif val > 0:
					logfile.write('  + ')
			logfile.write('\n')
	logfile.write('\n')

'''Performs a MECP optimization of the molecule specified in the filename.

Based loosely on MainProgram.f and the CROSSING shell scripts.

Rather than being controlled by a shell script, this does all its interfacing directly from this code
in order to cut down on unnecessary file access and shell interpretation (which is slow compared to Python),
providing a... marginal... performance benefit. (At least in comparison to the actual Gaussian running time)
At the very least, this simplifies the process of MECP optimization.'''
def main(jobname):

	global _g09_time

	start_time = datetime.now()
	_g09_time = timedelta()

	def timing_report(logfile):
		calc_time = datetime.now() - start_time
		logfile.write("Total running time: " + str(calc_time)[:-4] )
		logfile.write(datetime.now().strftime("\n on %b %d, %Y at %I:%M:%S %p %Z\n") )
		if sys.version_info >= (2, 7):
			py_time = calc_time - _g09_time
			g09_percent = _g09_time.total_seconds() * 100.0 / calc_time.total_seconds()
			py_percent = py_time.total_seconds() * 100.0 / calc_time.total_seconds()
			logfile.write("  Gaussian running time: %s (%.3f%%)\n" % (str(_g09_time)[:-4], g09_percent) )
			logfile.write("  MECP running time:     %s (%.3f%%)\n" % (str(py_time)[:-4], py_percent) )

	# read the program file (geometry and flags)
	cur_mol, options = read_prog(jobname + '.mecp')

	genopts = options['general']
	cutoffs = options['cutoffs']

	spinstates = genopts['spinstates']
	max_steps = genopts['max_steps']
	max_stepsize = genopts['max_stepsize']
	start_step = genopts['current_step']

	if genopts['pre_opt'][0] not in ['none', 0]:
		cur_mol = gaussian_opt(jobname, options, cur_mol)

	initial_mol = last_mol = cur_mol

	natoms = len(last_mol)
	nx = natoms * 3

	if options['hessian'] is not None and options['gradient'] is not None:
		# Resume from a previous state (the last completed step)
		last_grad = options['gradient']
		last_hess  = options['hessian']
		# Gotta recalculate the geometry
		cur_mol, cur_hess = update_coords(0, max_stepsize,
			(last_mol, last_mol), (last_grad, last_grad), last_hess)
	else:
		# Start from scratch: make an empty vector and an identity matrix
		last_grad = numpy.zeros(nx)
		last_hess  = 0.7 * matlib.identity(nx)

	with open(jobname + ".log", 'a') as logfile:
		logfile.write("Geometry Optimization of an MECP\n")
		logfile.write(datetime.now().strftime("Optimization started on %b %d, %Y at %I:%M:%S %p %Z") )
		logfile.write("\n\n")

		def save_state(name, step):
			options['geometry'] = last_mol
			options['hessian']  = last_hess
			options['gradient'] = last_grad
			options['general']['current_step'] = step
			mecpinput.write(name + '.mecp', options)

		def cleanup(signum, frame):
			gaussian.terminate()
			logfile.write('\n\n')
			logfile.write('Process received a TERM signal.\n')
			logfile.write('Run terminated.\n')
			timing_report(logfile)
			save_state(jobname + '.next')
			sys.exit(-2)

		signal.signal(signal.SIGTERM, cleanup)

		# Do up to the maximum optimization steps
		try:
			step = start_step;
			while step < (start_step + max_steps):

				logfile.write("Geometry at beginning of step %d:\n" % step)
				logfile.write(str(cur_mol))
				logfile.write('\n\n')

				# Call gaussian to do the heavy lifting and extract the results
				energies, gradients = call_gaussian(jobname, step, options, cur_mol)

				logfile.write("Energy of first spin state (%d): %.8f\n" % (spinstates[0], energies[0]) )
				logfile.write("Energy of second spin state (%d): %.8f\n" % (spinstates[1], energies[1]) )
				logfile.write('\n')

				# Calculate gradients needed to determine if there is a convergence
				diff_grad, par_grad, cur_grad = effective_gradient(energies, gradients)

				logfile.write("Overall Effective Gradient:\n")
				write_gradient(logfile, cur_grad)
				logfile.write('\n')

				logfile.write("Difference Gradient: (RMS: %.6f)\n" % rms(diff_grad))
				write_gradient(logfile, diff_grad)
				logfile.write('\n')

				logfile.write("Parallel Gradient: (RMS: %.6f)\n" % rms(par_grad))
				write_gradient(logfile, par_grad)
				logfile.write('\n')

				# Update the geometry of the molecule
				next_mol, cur_hess = update_coords(step, max_stepsize,
					(last_mol, cur_mol), (last_grad, cur_grad), last_hess)

				if genopts['show_hessian'] == 'full':
					logfile.write("Inverse Hessian Matrix at this step:\n")
					write_hessian(logfile, cur_hess)
				elif genopts['show_hessian'] == 'sign':
					logfile.write("Sign of Inverse Hessian Matrix at this step:\n")
					write_hessian_sign(logfile, cur_hess)

				# If there is a convergence, then we can stop.
				if test_convergence(logfile, (cur_mol, next_mol), energies, cur_grad, cutoffs):
					logfile.write("CONVERGENCE FOUND!\n")
					logfile.write("Initial Geometry\n")
					logfile.write(str(initial_mol))
					logfile.write('\n\n')

					logfile.write("Final Geometry\n")
					logfile.write(str(cur_mol))
					logfile.write('\n\n')

					# next_mol.normalize()
					# logfile.write("Normalized Final Geometry\n")
					# logfile.write(str(next_mol))
					# logfile.write("\n\n")

					logfile.write("Haikus are easy\nBut sometimes they don't make sense\nRefrigerator\n\n")

					timing_report(logfile)

					return 0

				# Prepare for the next step
				last_mol = cur_mol
				cur_mol = next_mol
				last_grad = cur_grad
				last_hess = cur_hess

				if options['general']['read_later']:
					mecpinput.addreadtoboth(options);

				step += 1

				print "STEP " + str(step)
				sys.stdout.flush()
				save_state(jobname + '.next', step)

				logfile.write('=%s==\n\n' % ( ("=STEP" + str(step)) * 10 ) )
				logfile.flush()

		except MECPException as e:
			print str(e)
			logfile.write('\n\n')
			logfile.write(str(e))
			logfile.write('\nRun terminated.\n')
			timing_report(logfile)
			return -1

		if step >= (start_step + max_steps):
			logfile.write('\n\n')
			logfile.write('Reached maximum number of steps.\n')
			logfile.write('Run terminated.\n')
			timing_report(logfile)
			return 1

	return 0

if __name__ == "__main__":
	job = sys.argv[1]
	job = job[:-5] if job.endswith(".mecp") else job
	status = main(job)
	sys.exit(status)
