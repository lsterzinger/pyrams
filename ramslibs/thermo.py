""" Contains functions for derived thermodynamic variables """
cp = 1004
R = 287
p0 = 1000


def temperature(theta, pi, celsius=False):
    """
    Calculates temperature from RAMS output (Exner function)

    Arguments
    ---------
    theta : ndarray
        Potential temperature (Kelvin)

    pi : ndarray
        Exner function * Cp (J/(kg*K))

    celsius : bool
        Default=False, return values in degrees Celsius

    Returns
    -------
    temperature : ndarray
        Temperature (Kelvin or Celsius)

    """

    T = theta * pi / cp
    if celsius:
        T = T - 273.15
    return T


def pressure(pi):
    """
    Calculates pressure from Exner function

    Arguments
    ---------
    pi : ndarray
        Exner function * Cp (J/(kg*K))

    Returns
    -------
    pressure : ndarray
        Pressure (hPa)
    """

    p = p0*(pi/cp)**(cp/R)
    return p


def wsat(theta, pi):
    """
    Calculates saturation water vapor pressure from
    temperature and pressure (RAMS specific calculation)

    Arguments
    ---------
    theta : ndarray
        Potential temperature (Kelvin)

    pi : ndarray
        Exner function * Cp (J/(kg*K))

    Returns
    -------
    wsat : ndarray
        Saturation water vapor pressure (hPa)
    """

    # RAMS specific calcuation

    T = temperature(theta, pi, celsius=True)
    p = pressure(pi)

    c0 = 0.6105851e3
    c1 = 0.4440316e2
    c2 = 0.1430341e1
    c3 = 0.2641412e-1
    c4 = 0.2995057e-3
    c5 = 0.2031998e-5
    c6 = 0.6936113e-8
    c7 = 0.2564861e-11
    c8 = -0.3704404e-13

    T[T < -80] = -80
    x = T
    es = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

    ws = 0.622 * es/(p*100 - es)
    return ws


def rh(rv, theta, pi):
    """
    Calculate relative humidity

    Arguments
    ---------
    rv : ndarray
        Water vapor pressure from RAMS output (hPa)

    theta : ndarray
        Potential temperature (Kelvin)

    pi : ndarray
        Exner function * Cp (J/(kg*K))

    Returns
    -------
    relative humidity : ndarray
        Relative humidity (fraction)
    """

    ws = wsat(theta, pi)
    rh = rv/ws * 100

    return rh


def mslp(temp, press, height):
    """
    Calculate the mean sea level pressure

    Parameters
    ------
    temp : numpy.ndarray
        Temperature in Celsius

    press: numpy.ndarray
        Pressure in hPa

    height: numpy.ndarray
        height of the terrain in meters

    Returns
    -------
    p0 : numpy.ndarray
        The MSLP pressure (hPa)
    """
    # p0 = press*(1-(0.0065*height)/(temp + 0.0065*height + 273.15))**-5.257
    p0 = press * (1-(0.0065*height)/(temp+0.0065*height + 273.15))**-5.257

    return(p0)
