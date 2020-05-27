from constants import Constants

class ThreeLevelLambdaSystem(Constants):
	def __init__(self,intensity1=1,intensity2=1,units='SI'):
		super().__init__(units)

		self._gamma_coherence 	= 300 								#[Hz] relaxation rate of the ground-state coherence
		self._gamma_width	 	= 1e9 								#[Hz] the homogeneous width of the optical transition
		self._dlightshift_di 	= 36.06								#[Hz/W/m^2]
		self.hz_per_au			= 4.13413732e16						#[Hz/au]  self.electron_mass*self.elementary_charge^4/self.hbar^3  
		self.i_per_au			= 3.50944758e20						#[W/m^2/au]
		self._intensity1 		= intensity1 						#[W/m^2]
		self._intensity2 		= intensity2 						#[W/m^2]

	@property
	def intensity1(self):
		SIValue = self._intensity1
		AUValue = self._intensity1 / self.i_per_au
		return self.select_units(SIValue,AUValue)

	@intensity1.setter
	def intensity1(self,value):
		self._intensity1 = value #[W/m^2]

	@property
	def intensity2(self):
		SIValue = self._intensity2
		AUValue = self._intensity2 / self.i_per_au
		return self.select_units(SIValue,AUValue)

	@intensity2.setter
	def intensity2(self,value):
		self._intensity2 = value #[W/m^2]

	@property
	def total_intensity(self):
		return self.intensity1+self.intensity2
	
	@property
	def rabi1_squared(self):
		n = 1 #index of refraction
		dipole = 2.992*self.elementary_charge*self.bohr_radius #87Rb D1 (52S1/2 −→ 52P1/2) Transition Dipole Matrix Element
		e_field_squared = 2/(self.permittivity_free_space*self.c*n) * self.intensity1
		return (dipole/self.hbar)**2 * e_field_squared
	
	@property
	def rabi2_squared(self):
		n = 1 #index of refraction
		dipole = 2.992*self.elementary_charge*self.bohr_radius #87Rb D1 (52S1/2 −→ 52P1/2) Transition Dipole Matrix Element
		e_field_squared = 2/(self.permittivity_free_space*self.c*n) * self.intensity2
		return (dipole/self.hbar)**2 * e_field_squared

	@property
	def normalized_gamma(self): 
		return self.gamma_coherence/self.gamma_width

	@property
	def normalized_rabi1_squared(self):
		return self.rabi1_squared/self.gamma_width**2

	
	@property
	def normalized_rabi2_squared(self):
		return self.rabi2_squared/self.gamma_width**2

	@property
	def normalized_dlightshift_di(self):
		SIValue=self.dlightshift_di/ self.gamma_width
		AUValue=self.dlightshift_di/(self.gamma_width*self.i_per_au)
		return self.select_units(SIValue,AUValue)

	@property
	def dlightshift_di(self):
		SIValue = self._dlightshift_di
		AUValue = self._dlightshift_di*(self.i_per_au/self.hz_per_au)
		return self.select_units(SIValue,AUValue)
	
	@property
	def gamma_coherence(self):
		SIValue = self._gamma_coherence
		AUValue = self._gamma_coherence / self.hz_per_au
		return self.select_units(SIValue,AUValue)

	@property
	def gamma_width(self):
		SIValue = self._gamma_width
		AUValue = self._gamma_width / self.hz_per_au
		return self.select_units(SIValue,AUValue)

	def excited_state_population(self,normalized_raman_detuning,light_shift_slope=None):
		if not light_shift_slope:
			ls_slope = self.normalized_dlightshift_di
		else:
			ls_slope = light_shift_slope
		num = self.normalized_rabi1_squared*self.normalized_rabi2_squared*(self.normalized_gamma + 2*(self.normalized_rabi1_squared+self.normalized_rabi2_squared))
		denom = (self.normalized_gamma + 2*(self.normalized_rabi1_squared+self.normalized_rabi2_squared))**2 + (normalized_raman_detuning- ls_slope*self.total_intensity)**2 
		return num/denom
