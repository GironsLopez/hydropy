#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
hydropy.evap
================

**A package to perform calculations related to evaporation.**

This package is intended to provide functions and methods to calculate
mass and energy fluxes between the surface of the earth and the atmsophere.

.. author:: Marc Girons Lopez

"""

import math

# Physical constants

CP = 1.013e-3   # Specific heat at constant pressure (MJ kg-1 deg C-1)
E = 0.622       # Ratio molecular weight of water vapour / dry air
LHV = 2.45      # Latent heat of vaporisation (MJ kg-1)
SC = 0.0820     # Solar constant (MJ m-2 min-1)
WD = 1000       # Density of water (kg m-3)
SBC = 4.903e-9  # Stefan-Boltzmann constant (MJ K-4 m-2 day-1)
SHC = 2.1       # Soil heat capacity (MJ m-3 deg C-1)


def atmospheric_pressure(elev):
    """
    Approximation of the atmospheric pressure
    as a function of elevation (Allen et al. 1998).

    NOTE: This equation is a simplification of the ideal gas law,
    assuming 20 deg C and a standard atmosphere.

    Parameters
    ----------
    elev: float
        Elevation (m a.s.l.).

    Returns
    -------
    float
        Atmospheric pressure (kPa).

    """
    return 101.3 * math.pow((293 - 0.0065 * elev) / 293, 5.26)


def psychrometric_constant(elev):
    """
    Approximation of the psychrometric constant as a function
    of elevation (atmospheric pressure) and temperature (Allen et al. 1998).

    NOTE: As an average atmospheric pressure is used for each location,
    the physchrometric constant is kept constant for each location.

    Parameters
    ----------
    elev : float
        Elevation (m a.s.l.).

    Returns
    -------
    float
        Psychrometric constant (kPa deg C-1).

    """
    # Calculate the atmospheric pressure
    p = atmospheric_pressure(elev)

    return (CP * p) / (E * LHV)


def celsius2kelvin(temp_c):
    """
    Convert degrees Celsius to degrees Kelvin (Allen et al. 1998).

    NOTE: In some calculation procedures, temperature is required in
    degrees Kelvin. The Kelvin and Celsius temperature scales have
    the same interval.

    Parameters
    ----------
    temp_c : float
        Temperature (deg C).

    Returns
    -------
    float
        Temperature (deg K).

    """
    return temp_c + 273.16


def mean_daily_air_temperature(t_min, t_max):
    """
    Calculation of the mean daily air temperature (Allen et al. 1998).

    NOTE: "For standardization, t_mean for 24-hour periods is defined as the
    mean of the daily maximum and minimum temperatures rather than as the
    average of hourly temperature measurements."

    Parameters
    ----------
    t_min : float
        Daily minimum air temperature (deg C).
    t_max : float
        Daily maximum air temperature (deg C).

    Returns
    -------
        Daily average air temperature (deg C).

    """
    return (t_max + t_min) / 2


def relative_humidity(e_a, e_s):
    """
    Calculate the relative humidity of air (Allen et al., 1998).

    Relative humidity expresses the degree of saturation of the air as a ratio
    of the actual (e_a) to the saturation (e_s) vapour pressure at the same
    temperature.

    Parameters
    ----------
    e_a : float
        Actual vapour pressure at a given temperature (kPa).
    e_s : float
        Saturation vapour pressure at a given temperature (kPa).

    Returns
    -------
    float
        Relative humidity (%).

    """
    return (e_a / e_s) * 100


def saturation_vapour_pressure(temp):
    """
    Approximation of the mean saturation vapour pressure as a function
    of air temperature (Allen et al. 1998).

    NOTE: If we input dewpoint temperature instead of air temperature,
    this equation will return the actual vapour pressure since dewpoint
    temperature is the temperature to which the air needs to be cooled to
    make the air saturated. The relationship t_dew ~ t_min holds for locations
    where the cover crop of the station is well watered. However, particularly
    for arid regions, the air might not be saturated when its temperature is
    at its minimum. Hence, t_min might be greater than t_dew and a further
    calibration may be required to estimate dewpoint temperatures. In these
    situations, t_dew might be better approximated by subtracting 2-3 deg C
    from t_min. In humid and subhumid climates, t_min and t_dew measured in
    early morning may be less than t_dew measured during the daytime because
    of condensation of dew during the night.

    Parameters
    ----------
    temp : float
        Daily average air temperature (deg C).

    Returns
    -------
    float
        Mean saturation vapour pressure (kPa).

    """
    return 0.6108 * math.exp((17.27 * temp) / (temp + 237.3))


def mean_saturation_vapour_pressure(t_min, t_max):
    """
    Calculation of the mean daily saturation vapour pressure
    (Allen et al. 1998).

    NOTE: "Due to the non-linearity of the saturation vapour pressure equation,
    the mean saturation vapour pressure for a day, week, decade or month should
    be computed as the mean between the saturation vapour pressure at the mean
    daily maximum and minimum air temperatures for that period"

    Parameters
    ----------
    t_min : float
        Daily minimum air temperature for the period of interest (deg C).
    t_max : float
        Daily maximum air temperature for the period of interest (deg C).

    Returns
    -------
    float
        Mean saturation vapour pressure for the period of interest (kPa).

    """
    # Calculate the saturation vapour pressure for t_min
    es_min = saturation_vapour_pressure(t_min)

    # Calculate the saturation vapour pressure for t_max
    es_max = saturation_vapour_pressure(t_max)

    return (es_max + es_min) / 2


def saturation_vapour_pressure_slope(temp):
    """
    Approximation of the slope of the saturation vapour pressure curve
    as a function of temperature (Allen et al. 1998).

    Parameters
    ----------
    temp : float
        Daily average air temperature (deg C).

    Returns
    -------
    float
        Slope of the saturation vapour pressure curve (kPa deg C-1).

    """
    # Calculate the saturation vapour pressure
    es = saturation_vapour_pressure(temp)

    return (4098 * es) / math.pow(temp + 237.3, 2)


def actual_vapour_pressure(
        t_min, rh_max, t_max=None, rh_min=None, rh_mean=None):
    """
    Approximation of the actual vapour pressure as derived from relative
    humidity data (Allen et al. 1998).

    Parameters
    ----------
    t_min : float
        Daily minimum air temperature (deg C).
    rh_max : float
        Daily maximum relative humidity (%).
    t_max : float, optional
        Daily maximum air temperature (deg C), default is None.
    rh_min : float, optional
        Daily minimum relative humidity (%), default is None.
    rh_mean : float, optional
        Daily mean relative humidity (%), default is None.

    Returns
    -------
    float
        Actual vapour pressure (kPa).

    Raises
    ------
    ValueError
        If neither rh_min nor rh_mean are provided.

    """
    # Calculate the saturation vapour pressure for t_min
    es_min = saturation_vapour_pressure(t_min)

    if t_max is None:
        # When using equipment where errors in estimating rh_min can be large,
        # or when rh data integrity are in doubt, then one should use only
        # rh_max.
        return es_min * (rh_max / 100)

    else:
        # Calculate the saturation vapour pressure for t_max
        es_max = saturation_vapour_pressure(t_max)

        if rh_min is not None:
            # For periods of a week, ten days or a month, rh_max and rh_min
            # are obtained by dividing the sum of the daily values by the
            # number of days in the period.
            return (es_min * (rh_max / 100) + es_max * (rh_min / 100)) / 2

        elif rh_mean is not None:
            # In absence of rh_max and rh_min, one can use the following eq
            # where rh_mean is the average between rh_max and rh_min. However,
            # this equation is less desirable than the previous ones.
            return (rh_mean / 100) * ((es_max + es_min) / 2)

        else:
            raise ValueError('Either rh_min or rh_mean must be provided.')


def vapour_pressure_deficit(
        t_max, t_min, t_dew=None, rh_max=None, rh_min=None, rh_mean=None):
    """
    Approximation of the vapour pressure deficit as a function of min and
    max temperature and relative humidity (Allen et al. 1998).

    Parameters
    ----------
    t_max : float
        Daily maximum air temperature (deg C).
    t_min : float
        Daily minimum air temperature (deg C).
    t_dew : float, optional
        Dewpoint temperature (deg C), default is None.
    rh_max : float, optional
        Daily maximum relative humidity (%), default is None.
    rh_min : float, optional
        Daily minimum relative humidity (%), default is None.
    rh_mean : float, optional
        Daily mean relative humidity (%), default is None.

    Returns
    -------
    float
        Vapour pressure deficit (kPa).

    Raises
    ------
    ValueError
        If neither t_dew nor t_max are provided.

    """
    # Calculate the saturation vapour pressure
    es = mean_saturation_vapour_pressure(t_min, t_max)

    # Calculate the actual vapour pressure
    if t_dew is not None:
        ea = saturation_vapour_pressure(t_dew)

    else:
        if rh_max is None:
            raise ValueError('Either t_dew or rh_max must be provided.')

        ea = actual_vapour_pressure(
                t_min, rh_max, t_max=t_max, rh_min=rh_min, rh_mean=rh_mean)

    return es - ea


def deg2rad(deg):
    """
    Conversion between decimal degrees and radians (Allen et al. 1998).

    Parameters
    ----------
    deg : float
        Decimal degrees (decimal deg).

    Returns
    -------
    float
        Radians (Rad).

    """
    return (math.pi / 180) * deg


def relative_shortwave_radiation(rs, rso):
    """
    Calculate the ratio of the solar radiation to the clear-sky solar
    radiation (Allen et al. 1998).

    Parameters
    ----------
    rs : float
        Solar radiation (MJ m-2 day-1).
    rso : float
        Clear-sky solar radiation (MJ m-2 day-1).

    Returns
    -------
    float
        Relative shortwave radiation (-).

    """
    return rs / rso


def relative_sunshine_duration(n_act, n_max):
    """
    Calculate the ratio of the actual duration of sunshine to the maximum
    possible duration of sunshine or daylight hours (Allen et al. 1998).

    This ratio expresses the cloudiness of the atmosphere.

    Parameters
    ----------
    n_act : float
        Actual duration of sunshine (hours).
    n_max : float
        Maximum possible duration of sunshine (hours).

    Returns
    -------
    float
        Relative sunshine duration (-).

    """
    return n_act / n_max


def radiation2evaporation(radiation):
    """
    Convert radiation expressed in MJ m-2 day-1 to equivalent evaporation
    in mm day-1 by using a conversion factor equal to the inverse of the
    latent heat of vaporisation (Allen et al. 1998).

    Parameters
    ----------
    radiation : float
        Radiation (MJ m-2 day-1).

    Returns
    -------
    float
        Evaporation (mm day-1).

    """
    return (radiation / (LHV * WD)) * 1000


def inverse_relative_distance_earth_sun(doy):
    """
    Calculate the inverse relative distance Earth-Sun for a given
    day of the year (Allen et al. 1998).

    Parameters
    ----------
    doy : int
        Julian day of the year (-).

    Returns
    -------
    float
        Inverse relative distance Earth-Sun.

    """
    return 1 + 0.033 * math.cos(((2 * math.pi) / 365) * doy)


def solar_declination(doy):
    """
    Calculate the solar declination for a given day of the year
    (Allen et al. 1998).

    Parameters
    ----------
    doy : int
        Julian day of the year (-).

    Returns
    -------
    float
        Solar declination (rad).

    """
    return 0.409 * math.sin(((2 * math.pi) / 365) * doy - 1.39)


def sunset_hour_angle(lat, dec):
    """
    Calculate the sunset hour angle for a given latitude and solar declination
    (Allen et al. 1998).

    Parameters
    ----------
    lat : float
        Latitude (rad).
    dec : float
        Solar declination (rad).

    Returns
    -------
    float
        Sunset hour angle (rad).

    """
    return math.acos(-math.tan(lat) * math.tan(dec))


def extraterrestrial_radiation(lat, doy):
    """
    Calculate the extraterrestrial radiation for daily periods
    for a given latitude and day of the year (Allen et al. 1998).

    NOTE: "For the winter months in latitudes greater than 55° (N or S),
    the equations for Ra have limited validity. Reference should be made to
    the Smithsonian Tables to assess possible deviations."

    Parameters
    ----------
    lat : float
        Latitude (decimal deg), positive (negative) in the
        northern (southern) hemisphere.
    doy : int
        Julian day of the year (-).

    Returns
    -------
    float
        Extraterrestrial radiation (MJ m-2 day-1).

    """
    # Convert the latitude from decimal deg to rad
    lat = deg2rad(lat)

    # Calculate the inverse relative distance Earth-Sun
    dr = inverse_relative_distance_earth_sun(doy)

    # Calculate the solar declination
    dec = solar_declination(doy)

    # Calculate the sunse hour angle
    omega = sunset_hour_angle(lat, dec)

    return ((24 * 60) / math.pi) * SC * dr * (
            omega * math.sin(lat) * math.sin(dec) +
            math.cos(lat) * math.cos(dec) * math.sin(omega))


def daylight_hours(lat, doy):
    """
    Calculate the maximum possible duration of sunshine or daylight hours
    (Allen et al. 1998).

    Parameters
    ----------
    lat : float
        Latitude (decimal deg), positive (negative) in the
        northern (southern) hemisphere.
    doy : int
        Julian day of the year (-).

    Returns
    -------
    float
        Daylight hours (hours).

    """
    # Convert the latitude from decimal deg to rad
    lat = deg2rad(lat)

    # Calculate the solar declination
    dec = solar_declination(doy)

    # Calculate the sunse hour angle
    omega = sunset_hour_angle(lat, dec)

    return omega * (24 / math.pi)


def solar_radiation(
        lat, doy, n_act=None, a_s=None, b_s=None,
        loc=None, t_max=None, t_min=None):
    """
    Calculate the solar radiation with the Angstrom formula, which relates
    solar radiation to extraterrestrial radiation and relative sumshine
    duration (Allen et al. 1998).

    NOTE: a_s + b_s represents the fraction of extraterrestrial radiation
    reaching the earth on clear days (n_act = n_max). Where no actual data
    are available and no calibration has been carried out for improved a_s
    and b_s parameters, the values a_s = 0.25 and b_s = 0.50 are
    recommended.

    Parameters
    ----------
    lat : float
        Latitude (decimal deg), positive (negative) in the
        northern (southern) hemisphere.
    doy : int
        Julian day of the year (-).
    n_act : int or float
        Actual duration of sunshine (hours).
    a_s : float, optional
        Regression constant, expressing the fraction of extraterrestrial
        radiation reaching the earth on overcast days (n_act = 0),
        default is None.
    b_s : float, optional
        Regression constant, default is None.
    loc : {'interior', 'coastal', 'island'}, optional
        Selection criterium for the adjustment coefficient if the sunlight
        hours parameter is not known, default is None.
    t_max : float, optional
        Maximum air temperature (deg C), default is None.
    t_min : float, optional
        Minimum air temperature (deg C), default is None.

    Returns
    -------
    float
        Solar radiation (MJ m-2 day-1).

    Raises
    ------
    ValueError
        If neither n_act nor loc are provided.
    ValueError
        If the value of loc is not recognised.
    ValueError
        If values for t_max and t_min are not provided.

    """
    # Calculate the extraterrestrial radiation.
    ra = extraterrestrial_radiation(lat, doy)

    if n_act is not None:

        # Calculate the daylight hours.
        n_max = daylight_hours(lat, doy)

        # Calculate the relative sunshine duration.
        n_rel = relative_sunshine_duration(n_act, n_max)

        if a_s is None:
            a_s = 0.25

        if b_s is None:
            b_s = 0.50

        return ra * (a_s + b_s * n_rel)

    else:

        if loc is None:
            raise ValueError('Either n_act or loc must be provided.')

        elif loc == 'island':
            # For island locations, where the land mass has a width
            # perpendicular to the coastline of 20km or less.
            # This relationship is only applicable for low altitudes (0-100m)
            b = 4  # MJ m-2 day-1 (empirical constant).
            return 0.7 * ra - b

        elif loc in ['interior', 'coastal']:

            if loc == 'interior':
                krs = 0.16

            elif loc == 'coastal':
                krs = 0.19

            if t_max is None or t_min is None:
                raise ValueError('Values for t_max and t_min ' +
                                 'must be provided.')

            return krs * math.sqrt(t_max - t_min) * ra

        else:
            raise ValueError('The value of loc is not recognised.')


def clear_sky_solar_radiation(ra, a_s=None, b_s=None, elev=None):
    """
    Calculate the clear-sky solar radiation (Allen et al. 1998).

    NOTE: The calculation of this parameter, when n_act = n_max, is
    required for computing the net longwave radiation.

    NOTE: a_s + b_s represents the fraction of extraterrestrial radiation
    reaching the earth on clear-sky days (n_act = n_max).

    Parameters
    ----------
    ra : float
        Extraterrestrial radiation (MJ m-2 day-1)
    a_s : float, optional
        Regression constant, expressing the fraction of extraterrestrial
        radiation reaching the earth on overcast days (n_act = 0),
        default is None.
    a_s : float, optional
        Regression constant, default is None.
    elev : int or float, optional
        Station elevation above the sea level (m a.s.l.).

    Returns
    -------
    float
        Clear-sky solar radiation (MJ m-2 day-1).

    Raises
    ------
    ValueError
        If neither a_s and b_s nor elev are provided.

    """
    if a_s is not None and b_s is not None:
        # For near sea level or when calibrated values for a_s and b_s are
        # available.
        return ra * (a_s + b_s)

    else:
        if elev is None:
            raise ValueError('Either a_s and b_s or elev need to be provided.')

        return ra * (0.75 + 2e-5 * elev)


def net_solar_radiation(rs, albedo=0.23):
    """
    Calculate the fraction of the solar radiation that is not reflected from
    the surface (Allen et al. 1998).

    Parameters
    ----------
    rs : float
        Solar radiation (MJ m-2 day-1).
    albedo : float, optional
        Albedo (-), default is 0.23 (value for a hypothetical
        grass reference crop.

    Returns
    -------
    float
        Net solar radiation (MJ m-2 day-1).

    """
    return (1 - albedo) * rs


def net_longwave_radiation(
        lat, doy, n_act, t_min, t_max, rh_max, rh_min=None,
        rh_mean=None, rs=None, a_s=None, b_s=None, elev=None):
    """
    Calculate the net longwave radiation (Allen et al. 1998).

    Parameters
    ----------
    lat : float
        Latitude (decimal deg), positive (negative) in the
        northern (southern) hemisphere.
    doy : int
        Julian day of the year (-).
    n_act : int or float
        Actual duration of sunshine (hours).
    t_min : float
        Daily minimum air temperature (deg C).
    t_max : float, optional
        Daily maximum air temperature (deg C), default is None.
    rh_max : float
        Daily maximum relative humidity (%).
    rh_min : float, optional
        Daily minimum relative humidity (%), default is None.
    rh_mean : float, optional
        Daily mean relative humidity (%), default is None.
    rs : float, optional
        Solar radiation (MJ m-2 day-1), default is None.
    a_s : float, optional
        Regression constant, expressing the fraction of extraterrestrial
        radiation reaching the earth on overcast days (n_act = 0),
        default is None.
    b_s : float, optional
        Regression constant, default is None.
    elev : int or float, optional
        Station elevation above the sea level (m a.s.l.).

    Returns
    -------
    float
        Net longwave radiation (MJ m-2 day-1).

    """
    # Transform the temperature values to deg K.
    t_min_k = celsius2kelvin(t_min)
    t_max_k = celsius2kelvin(t_max)

    # Calculate the actual vapour pressure
    ea = actual_vapour_pressure(
            t_min, rh_max, t_max=t_max, rh_min=rh_min, rh_mean=rh_mean)

    if rs is None:
        # Calculate the solar radiation if measurements are not available.
        rs = solar_radiation(lat, doy, n_act, a_s=a_s, b_s=b_s)

    # Calculate the extraterrestrial radiation
    ra = extraterrestrial_radiation(lat, doy)

    # Calculate the clear-sky solar radiation
    rso = clear_sky_solar_radiation(ra, a_s=a_s, b_s=b_s, elev=elev)

    # Calculate the relative shortwave radiation
    rs_rel = relative_shortwave_radiation(rs, rso)

    if rs_rel > 1:
        # The relative shortwave radiation must be limited to <= 1.
        rs_rel = 1

    return (SBC * ((math.pow(t_max_k, 4) + math.pow(t_min_k, 4)) / 2) *
            (0.34 - 0.14 * math.sqrt(ea)) * (1.35 * rs_rel - 0.35))


def net_radiation(rns, rnl):
    """
    Calculate the net radiation (Allen et al. 1998).

    Parameters
    ----------
    rns : float
        Incoming net shortwave radiation (MJ m-2 day-1).
    rnl : float
        Outgoing net longwave radiation (MJ m-2 day-1).

    Returns
    -------
    float
        Net radiation (MJ m-2 day-1).

    """
    return rns - rnl


def soil_heat_flux(
        period=None, temp_i=None, temp_i_1=None, delta_t=None, delta_z=None,
        t_month_i=None, t_month_i_minus=None, t_month_i_plus=None,
        hours=None, rn=None):
    """
    Calculate the soil heat flux (Allen et al. 1998).

    The soil heat flux is small compared to the net radiation, particularly
    when the surface is covered by vegetation and calculation time steps are
    24 hours or longer. For this reason a simple calculation procedure
    based on the idea that the soil temperature follows air temperature can
    be used for long time steps.

    NOTE: The effective soil depth, delta_z, is only 0.10-0.20 m for a time
    interval of one or a few days but might be 2 m or more for monthly periods.
    The soil heat capacity is related to its mineral composition and water
    content.

    Parameters
    ----------
    period : {'day', '10-day', 'month', 'hour'}, optional
        Period considered for the calculation of the soil heat flux,
        default is None.
    temp_i : float, optional
        Air temperature at time i (deg C), default is None.
    temp_i_1 : float, optional
        Air temperature at time i-1 (deg C), default is None.
    delta_t : int or float, optional
        Length of time interval (day), default is None.
    delta_z : float, optional
        Effective soil depth (m), default is None.
    t_month_i : float, optional
        Mean air temperature of the month i (deg C), default is None.
    t_month_i_minus : float, optional
        Mean air temperature of the previous month (dg C), default is None.
    t_month_i_plus : float, optional
        Mean air temperature of the next month (deg C), default is None.
    hours : {'night', 'day'}, optional
        Phase of the diurnal cycle, default is None.
    rn : float, optional
        Net radiation (MJ m-2 day-1), default is None.

    Returns
    -------
    float
        Soil heat flux (MJ m-2 day-1).

    Raises
    ------
    ValueError
        If the needed parameters are not provided or not recognised.

    """
    if period is None:
        # Simple calculation procedure for long time steps, based on the idea
        # that the soil temperature follows air temperature.
        if (temp_i is None or temp_i_1 is None or
                delta_t is None or delta_z is None):
            raise ValueError('If no period is provided, the parameters '
                             '"temp_i", "temp_i_1", "delta_t", and "delta_z" '
                             'need to be provided.')
        return SHC * ((temp_i + temp_i_1) / delta_t) * delta_z

    elif period in ['day', '10-day']:
        return 0

    elif period == 'month':
        # Considering a constant soil heat capacity of 2.1 MJ m-3 deg C-1.
        if t_month_i_plus is not None:
            if t_month_i_minus is None:
                raise ValueError('For monthly periods, either '
                                 '"t_month_i_minus", and "t_month_i" or '
                                 '"t_month_i_plus" need to be provided.')
            return 0.07 * (t_month_i_plus - t_month_i_minus)
        else:
            if t_month_i is None or t_month_i_minus is None:
                raise ValueError('For monthly periods, either '
                                 '"t_month_i_minus", and "t_month_i" or '
                                 '"t_month_i_plus" need to be provided.')
            return 0.14 * (t_month_i - t_month_i_minus)

    elif period == 'hour':
        # For hourly (or shorter) calculations, G beneath a dense cover of
        # grass does not correlate well with air temperature.
        if hours is None or rn is None:
            raise ValueError('For hourly periods, both "hours" and "rn" '
                             'parameters need to be provided.')
        if hours == 'day':
            return 0.1 * rn
        elif hours == 'night':
            return 0.5 * rn

    else:
        raise ValueError('The provided value for period is not recognised.')


def normalised_wind_speed(elev, u):
    """
    Adjust wind speed data obtained from instruments placed at
    elevations other than the standard height of 2 m (Allen et al. 1998).

    Parameters
    ----------
    elev : float
        Station elevation (m a.s.l.)
    u : float
        Wind speed at the predefined elevation (m s-1).

    Returns
    -------
    float
        Wind speed at 2 m (m s-1).

    """
    return u * (4.87 / (math.log(67.8 * elev - 5.42)))


def hargreaves_equation(lat, doy, t_min, t_max):
    """
    Calculate the reference evapotranspiration using the Hargreaves equation
    (Allen et al. 1998).

    This equation should be used when solar radiation data, relative humidity
    data, and/or wind speed data are missing.

    NOTE: This equation has a tendency to underpredict under high wind
    conditions (u2 > 3 m s-1) and to overpredict under conditions of high
    relative humidity.

    Parameters
    ----------
    lat : float
        Latitude (decimal deg), positive (negative) in the
        northern (southern) hemisphere.
    doy : int
        Julian day of the year (-).
    t_min : float
        Daily minimum air temperature (deg C).
    t_max : float, optional
        Daily maximum air temperature (deg C), default is None.

    Returns
    -------
    float
        Reference evapotranspiration (mm day-1)

    """
    # Calculate the mean temperature
    t_mean = mean_daily_air_temperature(t_min, t_max)

    # Calculate the extraterrestrial radiation.
    ra = extraterrestrial_radiation(lat, doy)

    return 0.0023 * (t_mean + 17.8) * math.pow(t_max - t_min, 0.5) * ra


def penman_monteith_equation(lat, elev, doy, albedo=0.23):
    """
    Calculate the crop reference evapotranspiration using the FAO
    Penman-Monteith equation (Allen et al. 1998).

    NOTE: The FAO Penman-Monteith equation determines the evapotranspiration
    from the hypothetical grass reference surface and provides a standard to
    which evapotranspiration in different periods of the year or in other
    regions can be compared and to which the evapotranspiration from other
    crops can be related.

    Parameters
    ----------

    Returns
    -------

    """
    # Calculate the solar radiation.
    rs = solar_radiation(lat, doy, n_act=None, a_s=None, b_s=None, loc=None, t_max=None, t_min=None)

    # Calculate the net solar radiation.
    rns = net_solar_radiation(rs, albedo=albedo)

    # Calculate the net longwave radiation.
    rnl = net_longwave_radiation(lat, doy, n_act, t_min, t_max, rh_max, rh_min=None, rh_mean=None, rs=None, a_s=None, b_s=None, elev=None)

    # calculate the vapour pressure deficit.
    e = vapour_pressure_deficit(t_max, t_min, t_dew=None, rh_max=None, rh_min=None, rh_mean=None)

    # Calculate the slope vapour pressure curve.
    delta = saturation_vapour_pressure_slope(temp)

    # Calculate the pshychrometric constant
    gamma = psychrometric_constant(elev, temp)

    # Calculate the soil heat flux density.
    g = soil_heat_flux(temp_i, temp_i_1, delta_t, delta_z)

    # Calculate the net radiation.
    rn = net_radiation(rns, rnl)

    return (0.408 * delta * (rn - g) + gamma * (900 / (temp + 273)) * u2 * (es - ea)) / (delta + gamma * (1 + 0.34 * u2))
