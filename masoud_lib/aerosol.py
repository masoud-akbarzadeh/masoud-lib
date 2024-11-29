def c_cunningham(d_p: float, mean_free_path: float = 65e-9) -> float:
    """
    Calculate the Cunningham correction factor.

    Parameters
    ----------
    d_p : float
        Particle diameter in meters.
    mean_free_path : float, optional
        Mean free path of the gas in meters. Default is 65e-9 m.

    Returns
    -------
    float
        Cunningham correction factor (dimensionless).

    Examples
    --------
    >>> c_cunningham(1e-6)
    1.0000814
    """
    kn = 2 * lamda / d_p
    return 1 + kn * (1.257 + 0.4 * np.exp(-1.1 / kn))


def mean_free_path(temperature: float, pressure: float):
    # mean free path of air calculator
    # T is the temperature in Kelvin
    # P is the pressure in Pascals
    # lamda is the mean free path of the gas in meters
    R = 8.314  # J/(mol K) gas constant
    M = 0.0289647  # kg/mol molar mass of air
    air = Fluid(FluidsList.Air).with_state(Input.pressure(pressure), Input.temperature(temperature - 273.15))
    viscosity = air.dynamic_viscosity  # Pa s dynamic viscosity
    return 2 * viscosity / (pressure * np.sqrt(8 * M / (np.pi * R * temperature)))
