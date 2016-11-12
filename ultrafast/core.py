"""Ultrafast core module

This module contains the core functionality required for performing ultrafast
optics calculations. It currently focusses on providing classes for describing
materials commonly employed in ultrafast optics.

This file is part of ultrafast. ultrafast is free software: you can
redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the
License, or any later version. ultrafast is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details. You should have received a copy of the GNU
General Public License along with ultrafast. If not, see
<http://www.gnu.org/licenses/>.

Copyright © 2016 Marcelo J P Alcocer
"""

# 	Package Structure
# 	==================
#
# 	Aim of package is to provide tools for performing ultrafast optical
# 	calculations. These typically involve several different parameters, e.g.
# 	material properties, geometrical constructions, etc.
#
# 	Whilst eventually wish the package to encompass many of these parameters,
# 	perhaps best to start by only handling dispersive materials. These form the
# 	basis of many calculations (e.g. compression, phase-matching, etc) and are
# 	often the most annoying to work with analytically. As such, numerical tools
# 	could be of real use.
#
# 	Can keep the rest of the current classes as examples to show how to use the
# 	materials in more complex constructions (e.g. optical elements, optical
# 	constructions). Upon polishing, these can be brought into sub-packages.
#
# 	Desired functionality:
#
# 		>>> import ultrafast
# 		>>> a = ultrafast.RIIDMaterial()
# 			<ultrafast.RIIDMaterial class>
#
# 	...and eventually, something like:
#
# 		>>> import ultrafast.elements
# 		>>> a = ultrafast.elements.Prism()
# 			<ultrafast.elements.Prism class>
#
# 	Proposed package structure:
#
# 		ultrafast/							# SCM directory
# 			README.rst
# 			LICENSE.txt
# 			DEPENDENCIES.txt
# 			setup.py						# Setup file (for PyPi)
# 			ultrafast/						# Top level package (actual import)
# 				__init__.py
# 				core.py						# Core module (Material, Error, etc)
# 				...							# ...more modules...
# 				elements.py					# Optical elements (prisms, lenses, etc)
# 			tests/							# Unit tests (repo only, not imported)
# 				tests_1.py
# 				tests_2.py
# 				...
# 			docs/							# Dox (repo only, not imported)
# 				index.rst
# 			examples/						# Examples (not imported)
# 				compressor.py				# Prism compressor example
# 				phasematching.py			# Phase matching exmaple
# 				...							# ... more useful examples...
#

# Imports
from scipy.constants import pi, speed_of_light
from math import atan, sqrt, pow
import yaml
import urllib.request
from urllib.parse import urlparse


class Material:
	"""Material class"""

	name = None
	"""Material name

	Name of the material. This can be used for reference, e.g. in case the
	instance name is not related to the material
	"""

	_range_ = None
	"""Frequency range

	Property attribute. See setter and getter methods for further details.
	"""

	references = None
	"""Reference string

	Reference(s) for material properties, e.g. journal articles
	"""

	comments = None
	"""Comments"""

	_n = None
	"""Dispersion function

	Property attribute. See setter and getter methods for further details.
	"""

	def __init__(
		self,
		n,
		range_,
		name=None,
		references=None,
		comments=None
	):
		"""Material class init

		:param n:	Dispersion function
		:type n:	callable
		:param range_:	Valid range for dispersion function (low,high)
		:type range_: 	tuple
		:param name:	Material name
		:type name:		string
		:param references:	Reference(s) for material properties
		:type references:	string
		:param comments:	Comments
		:type comments:		string

		Base class describing a (dispersive) material utilized in ultrafast optics.

		The dispersion function *n* is a callable which takes one argument, the
		angular frequency in :math:`rad/fs`, and returns the refractive index at
		this angular frequency.

		The frequency range *range_* is a tuple of length 2 describing the lower and
		upper angular frequencies for which *n* is valid.
		"""
		self.n = n
		self.range_ = range_
		self.name = name
		self.references = references
		self.comments = comments

	def _assert_frequency(self, omega):
		"""Frequency assertion

		:param omega:	Angular frequency in :math:`rad / fs`
		:type omega:	numeric

		Asserts the angular frequency *omega* is within the valid range of the
		dispersion function as defined by *range_*

		"""
		if(not (self.range_[0] <= omega <= self.range_[1])):
			raise RangeError(
				omega,
				self.range_,
				"Angular frequency out of material range"
			)

	'''
	def _assert_incidence_angle(self, phi):
		range_ = (radians(-90), radians(90))
		if(not (range_[0] < phi < range_[1])):
			raise RangeError(
				phi,
				range_,
				"Incidence angle out of range"
			)
		pass
	'''

	@property
	def range_(self):
		"""Frequency range

		Range of angular frequencies (:math:`rad/fs`) over which dispersion function
		:attr:`n` is valid. A numeric tuple of the form (low,high).

		"""
		return(self._range_)

	@range_.setter
	def range_(self, value):
		'''Frequency range setter method

		- Asserts range length (tuple of length 2)
		- Asserts low-high ordering

		'''

		# Assert length
		if(len(value) is not 2):
			raise PropertySetError(
				"range_",
				"Frequency range is not a 2-tuple"
			)

		# Assert ordering
		if(value[0] > value[1]):
			value = value[::-1]

		# Set frequency range
		self._range_ = value

	@property
	def n(self):
		"""Dispersion function

		A callable which takes one argument, the angular frequency in :math:`rad/fs`,
		and returns the refractive index at this angular frequency.
		"""
		return(self._n)

	@n.setter
	def n(self, value):
		"""Dispersion function setter method

		- Assert callable
		- Prepends frequency assertion of first argument to function
		"""

		# Assert callable
		if(not callable(value)):
			raise PropertySetError(
				"n",
				"Refractive index function is not callable"
			)

		# Set refractive index function
		def n(*args):
			self._assert_frequency(args[0])
			return(value(*args))

		# Set refractive index function
		self._n = n

	def wavevector(self, omega):
		"""Effective wavevector

		:param omega:	Angular frequency in :math:`rad / fs`
		:type omega:	float

		Returns the effective wavevector (:math:`\\omega n / c`) at the angular
		frequency *omega*
		"""

		# Assert frequency
		self._assert_frequency(omega)

		# Return wavevector
		return(omega * self.n(omega) / c)

	def brewster(self, omega, inc_mat=None):
		"""Brewster angle

		:param omega:	Angular frequency in :math:`rad / fs`
		:type omega:	float
		:param inc_mat:	Incident material
		:type inc_mat:	:class:`ultrafast.core.Material`

		Returns the Brewster angle for light rays incident from the material
		*inc_mat*. If None, *inc_mat* is assumed to be :attr:`air`.
		"""

		# Assert frequency
		self._assert_frequency(omega)

		# Define default external material
		if(inc_mat is None):
			inc_mat = air

		# Return brewster angle
		return(atan(self.n(omega) / inc_mat.n(omega)))


class RIIDMaterial(Material):
	"""RefractiveIndex.info material class"""

	type_ = None
	"""Dispersion data type

	String describing the dispersion data type as defined in the
	RefractiveIndex.info database, e.g.: ``formula 1``, ``tabulated k``, etc
	"""

	def __init__(self, db):
		"""RIIDMaterial class init

		:param db:	Database entry
		:type db:	string

		Class describing a dispersive material catalogued in the
		`RefractiveIndex.info <http://www.refractiveindex.info>`_ database.

		An instance is built from the RefractiveIndex.info database entry *db*.
		Entries are stored as YAML files. As such, *db* may be a path to a local YAML
		file, or a URL to a remote YAML file accessible via HTTP.
		"""

		# YAML keys
		keys = {
			"entry": {
				"data": "DATA",
				"ref": "REFERENCES",
				"com": "COMMENTS"
			},
			"data": {
				"type": "type",
				"range": "range",
				"coeff": "coefficients"
			}
		}

		# Fetch remote YAML database entry
		if(urlparse(db).scheme != ""):
			file = urllib.request.urlopen(db)

		# Fetch local YAML database entry
		else:
			file = open(db, newline="\r\n")

		entry = yaml.safe_load(file)
		file.close()

		# Extract dispersion function and range
		n = None
		range_ = None
		if((keys["entry"]["data"] in entry)):
			for datum in entry[keys["entry"]["data"]]:
				self.type_ = datum[keys["data"]["type"]]

				# Analytical dispersion
				if(self.type_.startswith("formula")):

					# Parse range
					range_ = [
						frequency(float(x)) for x in
						reversed(datum[keys["data"]["range"]].split())
					]

					# Parse coefficients
					coefficients = [
						float(x) for x in
						datum[keys["data"]["coeff"]].split()
					]

					# Define dispersion function frequency wrapper
					def frequency_wrap(fun_of_lambda):
						def fun_of_omega(omega):
							# self._assert_frequency(omega)
							return(fun_of_lambda(wavelength(omega)))
						return(fun_of_omega)

					# Construct frequency wrapped dispersion function
					formula = int(self.type_.split()[1])
					if formula == 1:
						@frequency_wrap
						def n(lambda_):
							"""Formula 1 - Sellmeier (preferred)"""
							n2 = 1 + coefficients[0]
							for x in [
								coefficients[i:i + 2] for i in
								range(1, len(coefficients), 2)
							]:
								n2 += (
									(x[0] * pow(lambda_, 2)) /
									(pow(lambda_, 2) - pow(x[1], 2))
								)
							return(sqrt(n2))
					elif formula == 2:
						@frequency_wrap
						def n(lambda_):
							"""Formula 2 - Sellmeier-2"""
							n2 = 1 + coefficients[0]
							for x in [
								coefficients[i:i + 2] for i in
								range(1, len(coefficients), 2)
							]:
								n2 += (
									(x[0] * pow(lambda_, 2)) /
									(pow(lambda_, 2) - x[1])
								)
							return(sqrt(n2))
					elif formula == 3:
						@frequency_wrap
						def n(lambda_):
							"""Formula 3 - Polynomial"""
							n2 = coefficients[0]
							for x in [
								coefficients[i:i + 2] for i in
								range(1, len(coefficients), 2)
							]:
								n2 += x[0] * pow(lambda_, x[1])
							return(sqrt(n2))
					elif formula == 4:
						@frequency_wrap
						def n(lambda_):
							"""Formula 4 - RefractiveIndex.info"""
							n2 = coefficients[0]
							for x in [
								coefficients[i:i + 4] for i in
								[1, 5]
							]:
								n2 += (
									(x[0] * pow(lambda_, x[1])) /
									(pow(lambda_, 2) - pow(x[2], x[3]))
								)
							for x in [
								coefficients[i:i + 2] for i in
								range(9, len(coefficients), 2)
							]:
								n2 += x[0] * pow(lambda_, x[1])
							return(sqrt(n2))
					elif formula == 5:
						@frequency_wrap
						def n(lambda_):
							"""Formula 5 - Cauchy"""
							n = coefficients[0]
							for x in [
								coefficients[i:i + 2] for i in
								range(1, len(coefficients), 2)
							]:
								n += x[0] * pow(lambda_, x[1])
							return(n)
					elif formula == 6:
						@frequency_wrap
						def n(lambda_):
							"""Formula 6 - Gases"""
							n = 1 + coefficients[0]
							for x in [
								coefficients[i:i + 2] for i in
								range(1, len(coefficients), 2)
							]:
								n += (
									x[0] /
									(x[1] - pow(lambda_, -2))
								)
							return(n)
					elif formula == 7:
						@frequency_wrap
						def n(lambda_):
							"""Formula 7 - Herzberger"""

							n = coefficients[0]
							n += coefficients[1] / (pow(lambda_, 2) - 0.028)
							n += coefficients[2] * pow(pow(lambda_, 2) - 0.028, -2)
							for i, x in enumerate(coefficients[3:], 1):
								n += x * pow(lambda_, 2 * i)
							return(n)
					elif formula == 8:
						@frequency_wrap
						def n(lambda_):
							"""Formula 8 - Retro"""
							alpha = (
								coefficients[0] +
								(
									(coefficients[1] * pow(lambda_, 2)) /
									(pow(lambda_, 2) - coefficients[2])
								) +
								(coefficients[3] * pow(lambda_, 2))
							)
							n2 = -(((2 * alpha) + 1) / (alpha - 1))
							return(sqrt(n2))
					elif formula == 9:
						@frequency_wrap
						def n(lambda_):
							"""Formula 9 - Exotic"""
							n2 = coefficients[0]
							n2 += (
								coefficients[1] /
								(pow(lambda_, 2) - coefficients[2])
							)
							n2 += (
								(coefficients[3] * (lambda_ - coefficients[4])) /
								(pow(lambda_ - coefficients[4], 2) + coefficients[5])
							)
							return(sqrt(n2))
					else:
						raise RangeError(
							formula,
							(1, 9),
							"RIID dispersion formula out of range"
						)

					# Break out of datum loop once dispersion function found
					break

				# Tabulated dispersion
				elif(self.type_.startswith("tabulated n")):
					print("Tabulated dispersion function")

					# Break out of datum loop once dispersion function found
					break

			# No dispersion function found (neither formula nor tabulated)
			if(n is None):
				raise UltrafastError("No dispersion data found in RIID entry")

		# Extract references
		references = None
		if(keys["entry"]["ref"] in entry):
			references = entry[keys["entry"]["ref"]]

		# Extract comments
		comments = None
		if(keys["entry"]["com"] in entry):
			comments = entry[keys["entry"]["com"]]

		# Call Material constructor
		Material.__init__(
			self, n, range_, name=db,
			references=references, comments=comments
		)


class UltrafastError(Exception):
	"""Ultrafast module base exception class

	Ultrafast module base exception class
	"""
	pass


class RangeError(UltrafastError):
	"""Out of range exception

	"""

	def __init__(self, value, valid, message):
		"""Out of range exception

		:param value:	Invalid value
		:type value:	float, int
		:param valid:	Valid range
		:type valid:	tuple
		:param message:	Error message
		:type message:	string

		Raised when the numeric value *value* is found to be out of the valid range
		defined by *valid*. *valid* is a tuple of length 2 in the form of (low,high).
		"""
		self.value = value
		self.valid = valid
		self.message = message

	def __str__(self):
		return(
			"{}. Value: {}. Valid range: {}".format(
				self.message,
				self.value,
				self.valid
			)
		)


class PropertySetError(UltrafastError):
	"""Property set exception"""

	def __init__(self, property_, message):
		"""Property set exception

		:param property_:	Property attribute
		:type property_:	string
		:param message:		Error message
		:type message:		string

		Raised when a attempting to set the property attribute *property_* to an
		invalid value.

		"""
		self.property_ = property_
		self.message = message

	def __str__(self):
		return(
			"{}. Property: {}.".format(
				self.message,
				self.property_
			)
		)


def _converter(value):
	"""Wavelength <-> angular frequency conversion

	:param value:	Wavelength in :math:`\\mu m` or angular frequency in :math:`rad
					/ fs`
	:type value:	float

	General method for conversion between wavelength and angular frequency.
	"""

	return(c * 2 * pi / value)


def frequency(lambda_):
	"""Wavelength to angular frequency conversion

	:param lambda_:	Wavelength in :math:`\\mu m`
	:type lambda_:	float

	Wavelength to angular frequency conversion.

	Returns the corresponding angular frequency in :math:`rad / fs`
	"""
	return(_converter(lambda_))


def wavelength(omega):
	"""Angular frequency to wavelength conversion

	:param omega:	Angular frequency in :math:`rad / fs`
	:type omega:	float

	Angular frequency to wavelength conversion.

	Returns the corresponding wavelength in :math:`\\mu m`
	"""
	return(_converter(omega))

c = speed_of_light * (1e-9)
"""Speed of light

Speed of light in :math:`\\mu m / fs`. Defined for convenience

"""

air = RIIDMaterial(
	"http://refractiveindex.info/database/other/mixed%20gases/air/Ciddor.yml"
)
"""Air

:class:`ultrafast.core.RIIDMaterial` describing air (Ciddor 1996). Defined for
convenience

"""
