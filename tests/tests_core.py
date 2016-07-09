"""Tests for core ultrafast functionality"""

import unittest
import ultrafast
import math
from scipy.constants import speed_of_light


class TestCoreMaterial(unittest.TestCase):

	def setUp(self):
		'''Instantiate test material'''

		def n(omega):
			'''Dummy refractive index function'''
			return(omega)

		# Instantiate test material
		range_ = (1, 10)
		self.mat = ultrafast.Material(
			n,
			range_,
			name="Test material name",
			references="Test material references",
			comments="Test material comments"
		)

	def test_init(self):
		'''Test initialization'''

		# Successful initialization
		self.assertIsInstance(
			self.mat,
			ultrafast.Material
		)

	def test_range_(self):
		'''Test range attribute'''

		# Correct length (2)
		self.assertEqual(
			len(self.mat.range_),
			2
		)

		# Fail set on bad length
		def set_long_range():
			self.mat.range_ = (0, 1, 2)
		self.assertRaises(
			ultrafast.PropertySetError,
			set_long_range
		)

		# Reverse ordering of decscending range
		range_ = self.mat.range_
		self.mat.range_ = range_[::-1]
		self.assertEqual(
			self.mat.range_,
			range_
		)

	def test_assert_frequency(self):
		'''Test _assert_frequency method'''

		# Frequency in range
		omega = (self.mat.range_[0] + self.mat.range_[1]) / 2
		self.assertIsNone(self.mat._assert_frequency(omega))

		# Fail on frequency out of range (low)
		omega_lo = 0.9 * self.mat.range_[0]
		self.assertRaises(
			ultrafast.RangeError,
			self.mat._assert_frequency,
			omega_lo
		)

		# Fail on frequency out of range (high)
		omega_hi = 1.1 * self.mat.range_[1]
		self.assertRaises(
			ultrafast.RangeError,
			self.mat._assert_frequency,
			omega_hi
		)

	def test_n(self):
		'''Test refractive index method'''

		# Correct return value
# 		omega = (self.mat.range_[0] + self.mat.range_[1]) / 2
# 		self.assertEqual(
# 			self.mat.n(omega),
# 			omega
# 		)

		# Greater than 1
		omega = (self.mat.range_[0] + self.mat.range_[1]) / 2
		self.assertGreater(
			self.mat.n(omega),
			1
		)

		# Fail on frequency out of range
		omega_lo = 0.9 * self.mat.range_[0]
		self.assertRaises(
			ultrafast.RangeError,
			self.mat.n,
			omega_lo
		)

		# Fail set on bad function
		def set_non_callable_n():
			self.mat.n = None
		self.assertRaises(
			ultrafast.PropertySetError,
			set_non_callable_n
		)

	def test_wavevector(self):
		'''Test wavevector method'''

		# Correct return value
		omega = (self.mat.range_[0] + self.mat.range_[1]) / 2
		self.assertEqual(
			self.mat.wavevector(omega),
			omega * self.mat.n(omega) / ultrafast.c
		)

	def test_brewster(self):
		'''Test brewster method'''

		# Correct return value - Air
		omega = (self.mat.range_[0] + self.mat.range_[1]) / 2
		self.assertEqual(
			self.mat.brewster(omega),
			math.atan(self.mat.n(omega) / ultrafast.air.n(omega))
		)

		# Correct return value - Self
		self.assertEqual(
			self.mat.brewster(omega, self.mat),
			math.atan(1)
		)


class TestCoreRIIDMaterial(TestCoreMaterial):

	def setUp(self):
		'''Instantiate test RIID material (air)'''

		self.mat = ultrafast.RIIDMaterial(
			"http://refractiveindex.info/database/other/mixed%20gases/air/Ciddor.yml"
		)

	def test_init(self):
		'''Test initialization'''

		# Remote database entry
		self.assertIsInstance(
			ultrafast.RIIDMaterial(
				"http://refractiveindex.info/database/other/mixed%20gases/air/Ciddor.yml"
			),
			ultrafast.RIIDMaterial
		)

		# Local database entry
		self.assertIsInstance(
			ultrafast.RIIDMaterial("../examples/Ciddor.yml"),
			ultrafast.RIIDMaterial
		)

		# Analytical dispersion
		self.assertIsInstance(
			ultrafast.RIIDMaterial("../examples/Ciddor.yml"),
			ultrafast.RIIDMaterial
		)

		# Fail on formula out of range
		self.assertRaises(
			ultrafast.RangeError,
			ultrafast.RIIDMaterial,
			"../examples/BadFormula.yml"
		)

		# Tabulated dispersion
		#
		# 	Fails until tabulated dispersions are handled properly
		#
# 		self.assertIsInstance(
# 			ultrafast.RIIDMaterial(
# 				"http://refractiveindex.info/database/main/H2/Leonard.yml"
# 			),
# 			ultrafast.RIIDMaterial
# 		)

		# Fail on no dispersion found
		self.assertRaises(
			ultrafast.UltrafastError,
			ultrafast.RIIDMaterial,
			"../examples/BadType.yml"
		)

		# Reference extraction
		self.assertIsNotNone(self.mat.references)

		# Comments extraction
		self.assertIsNotNone(self.mat.comments)


class TestErrors(unittest.TestCase):

	def test_UltrafastError(self):
		'''Test base error class'''
		pass

	def test_RangeError(self):
		'''Test range error class'''
		pass

	def test_PropertySetError(self):
		'''Test property set error class'''
		pass


class TestCoreFunctions(unittest.TestCase):

	def test_converter(self):
		'''Test converter function'''

		# Correct return value
		self.assertEqual(
			ultrafast.core._converter(1),
			ultrafast.c * 2 * math.pi
		)

	def test_frequency(self):
		'''Test frequency function'''

		# Correct return value
		omega = 1
		lambda_ = ultrafast.wavelength(omega)
		self.assertEqual(
			ultrafast.frequency(lambda_),
			omega
		)

	def test_wavelength(self):
		'''Test wvelength function'''

		# Correct return value
		lambda_ = 1
		omega = ultrafast.frequency(lambda_)
		self.assertEqual(
			ultrafast.wavelength(omega),
			lambda_
		)


class TestCoreAttributes(unittest.TestCase):

	def test_c(self):
		'''Test speed of light'''

		# Correct value
		self.assertEqual(
			ultrafast.c,
			speed_of_light * (1e-9)
		)

	def test_air(self):
		'''Test air material'''

		# Successful instantialization
		self.assertIsInstance(
			ultrafast.air,
			ultrafast.RIIDMaterial
		)


if __name__ == "__main__":
	unittest.main()
