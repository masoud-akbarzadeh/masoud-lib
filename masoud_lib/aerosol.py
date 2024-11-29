def c_cunningham(d_p, lamda=65E-9):
    # defining the function
    # Cunningham correction factor
    # Dp is the particle diameter in meters
    # lamda is the mean free path of the gas in meters
    # C is the Cunningham correction factor
    kn = 2 * lamda / d_p
    return 1 + kn * (1.257 + 0.4 * np.exp(-1.1 / kn))