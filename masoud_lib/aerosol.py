# For testing the module functions
if __name__ == '__main__':  # This is to test the module functions
    print("This is a module that contains functions that could be useful for my projects.")

# Importing libraries
import numpy as np
from pyfluids import Fluid, FluidsList, Input
from scipy import constants as const
import pandas as pd

# Global constants
PLANCK_CONSTANT = const.Planck  # J·s, Planck constant
GAS_CONSTANT = const.R  # J/(mol·K), molar gas constant
BOLTZMANN_CONSTANT = const.Boltzmann  # J/K, Boltzmann constant
GRAVITY_ACCELERATION = const.g  # m/s², acceleration due to gravity
AVOGADRO_NUMBER = const.Avogadro  # 1/mol, Avogadro's number

# Default air standard properties
STANDARD_PRESSURE = 101325  # Pa
STANDARD_TEMPERATURE = 273.15 + 25  # K
STANDARD_AIR = Fluid(FluidsList.Air).with_state(
    Input.pressure(STANDARD_PRESSURE),
    Input.temperature(STANDARD_TEMPERATURE - 273.15)
)

# Functions

def calc_cunningham_correction_factor(particle_diameter, mean_free_path=65e-9):
    """
    Calculate the Cunningham correction factor.

    Parameters
    ----------
    particle_diameter : float
        Particle diameter in meters (SI units).
    mean_free_path : float, optional
        Mean free path of the gas in meters (SI units).

    Default Values
    --------------
    mean_free_path : float
        Default is 65e-9 m.

    Returns
    -------
    float
        Cunningham correction factor (dimensionless).
    """
    knudsen_number = 2 * mean_free_path / particle_diameter
    correction_factor = 1 + knudsen_number * (1.257 + 0.4 * np.exp(-1.1 / knudsen_number))
    return correction_factor

def calc_mean_free_path_air(temperature: float, pressure: float):
    """
    Calculate the mean free path of air.

    Parameters
    ----------
    temperature : float
        Temperature in Kelvin (SI units).
    pressure : float
        Pressure in Pascals (SI units).

    Returns
    -------
    float
        Mean free path in meters (SI units).
    """
    R = GAS_CONSTANT  # J/(mol·K), gas constant
    M = 0.0289647  # kg/mol, molar mass of air
    air = Fluid(FluidsList.Air).with_state(
        Input.pressure(pressure),
        Input.temperature(temperature - 273.15)
    )
    viscosity = air.dynamic_viscosity  # Pa·s, dynamic viscosity
    mean_free_path = 2 * viscosity / (pressure * np.sqrt(8 * M / (np.pi * R * temperature)))
    return mean_free_path

def calc_reynolds_particle(
    particle_diameter: float,
    velocity: float,
    fluid_density=STANDARD_AIR.density,
    dynamic_viscosity=STANDARD_AIR.dynamic_viscosity
):
    """
    Calculate the Reynolds number for a particle.

    Parameters
    ----------
    particle_diameter : float
        Particle diameter in meters (SI units).
    velocity : float
        Particle velocity in meters per second (SI units).
    fluid_density : float, optional
        Density of the fluid in kg/m³ (SI units).
    dynamic_viscosity : float, optional
        Dynamic viscosity of the fluid in Pa·s (SI units).

    Default Values
    --------------
    fluid_density : float
        Default is standard air density.
    dynamic_viscosity : float
        Default is standard air dynamic viscosity.

    Returns
    -------
    float
        Reynolds number (dimensionless).
    """
    reynolds_number = (particle_diameter * fluid_density * velocity) / dynamic_viscosity
    return reynolds_number

def calc_settling_velocity_stokes(
    particle_diameter,
    density_particle: float,
    temperature: float,
    pressure: float
):
    """
    Calculate the settling velocity of a particle using Stokes' law.

    Parameters
    ----------
    particle_diameter : float or array_like
        Particle diameter(s) in meters (SI units).
    density_particle : float
        Density of the particle in kg/m³ (SI units).
    temperature : float
        Temperature in Kelvin (SI units).
    pressure : float
        Pressure in Pascals (SI units).

    Returns
    -------
    float or ndarray
        Settling velocity in meters per second (SI units).
    """
    if np.isscalar(particle_diameter):
        particle_diameter_array = np.array([particle_diameter])
    else:
        particle_diameter_array = particle_diameter

    velocities = []
    for pd in particle_diameter_array:
        g = GRAVITY_ACCELERATION  # m/s²
        mean_free_path = calc_mean_free_path_air(temperature, pressure)
        cunningham_factor = calc_cunningham_correction_factor(pd, mean_free_path)
        air = Fluid(FluidsList.Air).with_state(
            Input.pressure(pressure),
            Input.temperature(temperature - 273.15)
        )
        mu_f = air.dynamic_viscosity
        rho_f = air.density
        s_velocity = cunningham_factor * (density_particle * g * pd ** 2) / (18 * mu_f)
        Re = calc_reynolds_particle(pd, s_velocity, fluid_density=rho_f, dynamic_viscosity=mu_f)
        if Re < 1:
            velocities.append(s_velocity)
        else:
            m_p = np.pi * density_particle * pd ** 3 / 6
            for _ in range(100):
                c_d = 24 / Re * (1 + 3 / 16 * 0.43 * Re)
                s_velocity = np.sqrt((m_p * g) / (1 / 8 * np.pi * c_d * rho_f * pd ** 2))
                Re_new = calc_reynolds_particle(pd, s_velocity, fluid_density=rho_f, dynamic_viscosity=mu_f)
                if abs(Re_new - Re) < 0.01:
                    break
                else:
                    Re = Re_new
            velocities.append(s_velocity)

    velocities_array = np.array(velocities)

    if np.isscalar(particle_diameter):
        return velocities_array[0]
    else:
        return velocities_array

def calc_condensation_diameter_growth_rate(
    particle_diameter,
    vapor_concentration_bulk,
    vapor_concentration_saturation,
    diffusion_coefficient,
    density_particle=1e3
):
    """
    Calculate the condensation diameter growth rate of a particle.

    Parameters
    ----------
    particle_diameter : float
        Particle diameter in meters (SI units).
    vapor_concentration_bulk : float
        Concentration of the condensing vapor in the bulk gas (far from the particle) in molecules per cubic meter (SI units).
        Represents the ambient vapor concentration, also known as c_inf.
    vapor_concentration_saturation : float
        Saturation concentration of the condensing vapor at the particle surface in molecules per cubic meter (SI units).
        Represents the vapor concentration at equilibrium, also known as c_saturation.
    diffusion_coefficient : float
        Diffusion coefficient of the condensing vapor in m²/s (SI units).
    density_particle : float, optional
        Density of the particle in kg/m³ (SI units).

    Default Values
    --------------
    density_particle : float
        Default is 1e3 kg/m³.

    Returns
    -------
    float
        Diameter growth rate in meters per second (SI units).
    """
    beta_factor = calc_beta_correction_condensation(particle_diameter)
    diameter_growth_rate = beta_factor * 4 * diffusion_coefficient * \
        (vapor_concentration_bulk - vapor_concentration_saturation) / (density_particle * particle_diameter)
    return diameter_growth_rate

def calc_condensation_mass_growth_rate(
    particle_diameter,
    vapor_concentration_bulk,
    vapor_concentration_saturation,
    diffusion_coefficient_vapor,
    density_particle=1e3
):
    """
    Calculate the condensation mass growth rate of a particle.

    Parameters
    ----------
    particle_diameter : float
        Particle diameter in meters (SI units).
    vapor_concentration_bulk : float
        Concentration of the condensing vapor in the bulk gas (far from the particle) in molecules per cubic meter (SI units).
        Represents the ambient vapor concentration, also known as c_inf.
    vapor_concentration_saturation : float
        Saturation concentration of the condensing vapor at the particle surface in molecules per cubic meter (SI units).
        Represents the vapor concentration at equilibrium, also known as c_saturation.
    diffusion_coefficient_vapor : float
        Diffusion coefficient of the condensing vapor in m²/s (SI units).
    density_particle : float, optional
        Density of the particle in kg/m³ (SI units).

    Default Values
    --------------
    density_particle : float
        Default is 1e3 kg/m³.

    Returns
    -------
    float
        Mass growth rate in kilograms per second (SI units).
    """
    beta_factor = calc_beta_correction_condensation(particle_diameter)
    mass_growth_rate = 2 * np.pi * diffusion_coefficient_vapor * particle_diameter * \
        (vapor_concentration_bulk - vapor_concentration_saturation) * beta_factor
    return mass_growth_rate

def calc_beta_correction_condensation(particle_diameter, mean_free_path=65e-9):
    """
    Calculate the Dahneke correction factor for condensation.

    Parameters
    ----------
    particle_diameter : float
        Particle diameter in meters (SI units).
    mean_free_path : float, optional
        Mean free path of the gas in meters (SI units).

    Default Values
    --------------
    mean_free_path : float
        Default is 65e-9 m.

    Returns
    -------
    float
        Correction factor (dimensionless).
    """
    knudsen_number = 2 * mean_free_path / particle_diameter
    beta_factor = (1 + knudsen_number) / (1 + 2 * knudsen_number * (1 + knudsen_number))
    return beta_factor

def convert_units(conversion_key, value=1):
    """
    Convert units based on the provided conversion key.

    Parameters
    ----------
    conversion_key : str
        Key representing the conversion (e.g., 'cm3_to_m3').
    value : float, optional
        Value to be converted.

    Default Values
    --------------
    value : float
        Default is 1.

    Returns
    -------
    float
        Converted value.
    """
    conversion_dict = {
        # Length
        'cm3_to_m3': value * 1e-6,
        'm3_to_cm3': value * 1e6,
        # Time
        'hr_to_s': value * 3600,
        's_to_hr': value / 3600,
        # Volume
        'm3_to_L': value * 1000,
        'L_to_m3': value / 1000,
        # Mass
        'kg_to_g': value * 1000,
        'g_to_kg': value / 1000,
        'kg_to_ug': value * 1e9,
        'ug_to_kg': value * 1e-9,
        # Pressure
        'Pa_to_kPa': value * 1e-3,
        'kPa_to_Pa': value * 1e3,
        'atm_to_Pa': value * 101325,
        'Pa_to_atm': value / 101325,
        # Temperature
        'C_to_K': value + 273.15,
        'K_to_C': value - 273.15,
    }
    return conversion_dict.get(conversion_key, 'Invalid conversion')

def calc_mass_particle_from_diameter(particle_diameter, density_particle=1e3):
    """
    Calculate the mass of a spherical particle from its diameter.

    Parameters
    ----------
    particle_diameter : float
        Particle diameter in meters (SI units).
    density_particle : float, optional
        Density of the particle in kg/m³ (SI units).

    Default Values
    --------------
    density_particle : float
        Default is 1e3 kg/m³.

    Returns
    -------
    float
        Particle mass in kilograms (SI units).
    """
    particle_mass = np.pi * density_particle * particle_diameter ** 3 / 6
    return particle_mass

def calc_volume_sphere(diameter):
    """
    Calculate the volume of a sphere given its diameter.

    Parameters
    ----------
    diameter : float
        Diameter of the sphere in meters (SI units).

    Returns
    -------
    float
        Volume of the sphere in cubic meters (SI units).
    """
    radius = diameter / 2
    volume = (4 / 3) * np.pi * radius ** 3
    return volume

def calc_surface_area_sphere(diameter):
    """
    Calculate the surface area of a sphere given its diameter.

    Parameters
    ----------
    diameter : float
        Diameter of the sphere in meters (SI units).

    Returns
    -------
    float
        Surface area of the sphere in square meters (SI units).
    """
    radius = diameter / 2
    surface_area = 4 * np.pi * radius ** 2
    return surface_area

def calc_coagulation_coefficient(
    particle_diameter1,
    particle_diameter2,
    temperature=STANDARD_TEMPERATURE,
    density_particle1=1000,
    density_particle2=1000
):
    """
    Calculate the coagulation coefficient between two particles.

    Parameters
    ----------
    particle_diameter1 : float
        Diameter of particle 1 in meters (SI units).
    particle_diameter2 : float
        Diameter of particle 2 in meters (SI units).
    temperature : float, optional
        Temperature in Kelvin (SI units).
    density_particle1 : float, optional
        Density of particle 1 in kg/m³ (SI units).
    density_particle2 : float, optional
        Density of particle 2 in kg/m³ (SI units).

    Default Values
    --------------
    temperature : float
        Default is standard temperature (298.15 K).
    density_particle1 : float
        Default is 1000 kg/m³.
    density_particle2 : float
        Default is 1000 kg/m³.

    Returns
    -------
    float
        Coagulation coefficient in cubic meters per second (SI units).
    """
    # Source: Seinfeld, J. H., & Pandis, S. N. (2006). Atmospheric Chemistry and Physics.
    # 2nd Edition, Table 13.1

    particle_diameter1_array = np.array([particle_diameter1])
    particle_diameter2_array = np.array([particle_diameter2])

    air = Fluid(FluidsList.Air).with_state(
        Input.pressure(STANDARD_PRESSURE),
        Input.temperature(temperature - 273.15)
    )
    viscosity_air = air.dynamic_viscosity  # Pa·s

    cunningham_factors1 = np.array([
        calc_cunningham_correction_factor(diameter) for diameter in particle_diameter1_array
    ])
    cunningham_factors2 = np.array([
        calc_cunningham_correction_factor(diameter) for diameter in particle_diameter2_array
    ])
    diffusion_coefficient1 = BOLTZMANN_CONSTANT * temperature * cunningham_factors1 / (
        3 * np.pi * viscosity_air * particle_diameter1_array
    )
    diffusion_coefficient2 = BOLTZMANN_CONSTANT * temperature * cunningham_factors2 / (
        3 * np.pi * viscosity_air * particle_diameter2_array
    )

    mass1 = density_particle1 * np.pi * particle_diameter1_array ** 3 / 6
    mass2 = density_particle2 * np.pi * particle_diameter2_array ** 3 / 6

    mean_speed1 = np.sqrt(8 * BOLTZMANN_CONSTANT * temperature / (np.pi * mass1))
    mean_speed2 = np.sqrt(8 * BOLTZMANN_CONSTANT * temperature / (np.pi * mass2))

    mean_free_path1 = 8 * diffusion_coefficient1 / (np.pi * mean_speed1)
    mean_free_path2 = 8 * diffusion_coefficient2 / (np.pi * mean_speed2)

    g1 = 1 / (3 * particle_diameter1_array * mean_free_path1) * (
        (particle_diameter1_array + mean_free_path1) ** 3 -
        (particle_diameter1_array ** 2 + mean_free_path1 ** 2) ** (1.5)
    ) - particle_diameter1_array
    g2 = 1 / (3 * particle_diameter2_array * mean_free_path2) * (
        (particle_diameter2_array + mean_free_path2) ** 3 -
        (particle_diameter2_array ** 2 + mean_free_path2 ** 2) ** (1.5)
    ) - particle_diameter2_array

    denominator = (
        ((particle_diameter1_array + particle_diameter2_array) /
         (particle_diameter1_array + particle_diameter2_array + 2 * np.sqrt(g1 ** 2 + g2 ** 2))) +
        8 * (diffusion_coefficient1 + diffusion_coefficient2) /
        np.sqrt(mean_speed1 ** 2 + mean_speed2 ** 2) / (particle_diameter1_array + particle_diameter2_array)
    )

    coagulation_coefficient = 2 * np.pi * (particle_diameter1_array + particle_diameter2_array) * \
        (diffusion_coefficient1 + diffusion_coefficient2) / denominator

    return coagulation_coefficient
