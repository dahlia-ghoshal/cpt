import numpy as np

class Constants():
	def __init__(self,units='SI'):
		self.units = units
		print("Using {} units".format(units))

	def select_units(self,SIValue,AUValue):
		if self.units == 'SI':
			return SIValue
		elif self.units == 'AU':
			return AUValue
		else:
			raise ValueError

	@property
	def boltzmann(self):
		SIValue = 1.38064852e-23 #J/K
		AUValue = 1
		return self.select_units(SIValue,AUValue)
	
	@property
	def hbar(self): # Reduced Planck Constant
		SIValue = 1.054571596e-34 #[Js]
		AUValue = 1
		return self.select_units(SIValue,AUValue)

	@property
	def elementary_charge(self): #Electron Charge
		SIValue = 1.602176462e-19 #[C]
		AUValue = 1
		return self.select_units(SIValue,AUValue)

	@property
	def electron_mass(self): #electron rest mass
		SIValue = 9.10938188e-31 #[kg]
		AUValue = 1
		return self.select_units(SIValue,AUValue)

	@property
	def permittivity_free_space(self): #(4*pi*esp0 = 1/ke)
		SIValue = 8.854187817e-12 #[F⋅m−1]
		AUValue = 1/(4*np.pi)
		return self.select_units(SIValue,AUValue)
	
	@property
	def bohr_radius(self): # Bohr radius
		return 4*np.pi*self.permittivity_free_space*self.hbar**2/(self.electron_mass*self.elementary_charge**2)

	@property
	def c(self): #Speed of light
		SIValue = 2.99792458e8 #[m/s]
		AUValue = 137.036
		return self.select_units(SIValue,AUValue)

	
	