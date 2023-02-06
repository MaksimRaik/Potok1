import numpy as np
import matplotlib.pyplot as plt

R = 8.31

mu = np.asarray( [ 26.0, 32.0, 28.0, 44.0, 18.0 ] ) * 10 ** -3

w_after = np.asarray( [ 0.0703, 0.2166, 0.7130, 0.0, 0.0 ] )

w_before = np.asarray( [ 0.0, 0.0, 0.7130, 0.2383, 0.0487 ] )

Q = 3.38 * 10 ** 6

p0= 10.0 ** 5

T0 = 300.0

nu0 = 1.0 / p0 * np.sum( w_after / mu ) * R * T0

p01= 10.0 ** 5 + 5000

T01 = 250

nu01 = 1.0 / p01 * np.sum( w_after / mu ) * R * T01


a1 = np.asarray( [ 1.598133250 * 10 ** 5, -3.425563420 * 10 ** 4, 2.210371497 * 10 ** 4, 4.943783640 * 10 ** 4, -3.947960830 * 10 ** 4 ] )

a2 = np.asarray( [ -2.216675050 * 10 ** 3, 4.847000970 * 10 ** 2, -3.818461820 * 10 ** 2, -6.264292080 * 10 ** 2, 5.755731020 * 10 ** 2 ] )

a3 = np.asarray( [ 1.2657256100, 1.119010961, 6.082738360, 5.301813360, 9.31782653 ] )

a4 = np.asarray( [ -7.980170140 * 10 ** -3, 4.293889240 * 10 ** -3, -8.530914410 * 10 ** -3, 2.503600571 * 10 ** -3, 7.222712860 * 10 ** -3] )

a5 = np.asarray( [ 8.055801510 * 10 ** -6, -6.836300520 * 10 ** -7, 1.384646189 * 10 ** -5, -2.124700099 * 10 ** -7, -7.342557370 * 10 ** -6 ] )

a6 = np.asarray( [ -2.433944610 * 10 ** -9, -2.023372700 * 10 ** -9, -9.625793620 * 10 ** -9, -7.691486800 * 10 ** -10, 4.955043490 * 10 ** -9 ] )

a7 = np.asarray( [ -7.509454610 * 10 ** -14, 1.039040018 * 10 ** -12, 2.519705809 * 10 ** -12, 2.849979913 * 10 ** -13, -1.336933246 * 10 ** -12 ] )



b1 = np.asarray( [ 1.713792120 * 10 ** 6, -1.037939022 * 10 ** 6, 5.877124060 * 10 ** 5, 1.176969434 * 10 ** 5, 1.034972096 * 10 ** 6 ] )

b2 = np.asarray( [ -5.928970950 * 10 ** 3, 2.344830282 * 10 ** 3, -2.239249073 * 10 ** 3, -1.788801467 * 10 ** 3, -2.412698562 * 10 ** 3 ] )

b3 = np.asarray( [ 1.236115640 * 10 ** 1, 1.819732036, 6.066949220, 8.291543530, 4.646110780 ] )

b4 = np.asarray( [ 1.314706250 * 10 ** -4, 1.267847582 * 10 ** -3, -6.139685500 * 10 ** -4, -9.224778310 * 10 ** -5, 2.291998307  * 10 ** -3 ] )

b5 = np.asarray( [ -1.362869040 * 10 ** -7, -2.188067988 * 10 ** -7, 1.491806679 * 10 ** -7, 4.869635410 * 10 ** -9, -6.836830480 * 10 ** -7 ] )

b6 = np.asarray( [ 2.712746060 * 10 ** -11, 2.053719572 * 10 ** -11, -1.923105485 * 10 ** -11, -1.892063841 * 10 ** -12, 9.426468930 * 10 ** -11 ] )

b7 = np.asarray( [ -1.302086850 * 10 ** -15, -8.193467050 * 10 ** -16, 1.061954386 * 10 ** -15, 6.330675090 * 10 ** -16, -4.822380530 * 10 ** -15 ] )

def c_v( T ):

    global a1, a2, a3, a4, a5, a6, a7

    global b1, b2, b3, b4, b5, b6, b7

    global R

    if T <= 10.0 ** 3:

        return ( a1 / T ** 2 + a2 / T + a3 + a4 * T + a5 * T ** 2 + a6 * T ** 3 + a7 * T ** 4 - 1.0 ) * R

    elif T > 10.0 ** 3:

        return ( b1 / T ** 2 + b2 / T + b3 + b4 * T + b5 * T ** 2 + b6 * T ** 3 + b7 * T ** 4 - 1.0 ) * R

def gama( T ):

    global c_v, w_before, mu, R

    return 1.0 + R * np.sum( w_before / mu ) / np.sum( w_before * c_v( T ) / mu )

def c_v_( T ):

    global a1, a2, a3, a4, a5, a6, a7

    global b1, b2, b3, b4, b5, b6, b7

    global R

    if T <= 10.0 ** 3:

        return (- 2.0 * a1 / T ** 3 - a2 / T ** 2 + a4 + 2.0 * a5 * T + 3.0 * a6 * T ** 2 + 4.0 * a7 * T ** 3 ) * R

    elif T > 10.0 ** 3:

        return (- 2.0 * b1 / T ** 3 - b2 / T ** 2 + b4 + 2.0 * b5 * T + 3.0 * b6 * T ** 2 + 4.0 * b7 * T ** 3 ) * R

def c_v__( T ):

    global a1, a2, a3, a4, a5, a6, a7

    global b1, b2, b3, b4, b5, b6, b7

    global R

    if T <= 10.0 ** 3:

        return ( 6.0 * a1 / T ** 4 + a2 / T ** 3 + 2.0 * a5 + 6.0 * a6 * T + 12.0 * a7 * T ** 2 ) * R

    elif T > 10.0 ** 3:

        return ( 6.0 * b1 / T ** 4 - b2 / T ** 3 + 2.0 * b5 + 6.0 * b6 * T + 12.0 * b7 * T ** 2 ) * R

def gama__( T ):

    global c_v__, w_before, mu, R, c_v

    return R * np.sum(w_before / mu) / ( np.sum(w_before * c_v( T ) / mu) ) ** 2 * np.sum(w_before * c_v__( T ) / mu)
    - 2.0 * R * R * np.sum(w_before / mu) / ( np.sum(w_before * c_v( T ) / mu) ) ** 3 * ( np.sum(w_before * c_v__( T ) / mu) ) ** 2

def gama_( T ):

    global c_v_, w_before, mu, R, c_v

    return R * np.sum(w_before / mu) / ( np.sum(w_before * c_v( T ) / mu) ) ** 2 * np.sum(w_before * c_v_( T ) / mu)

def T_( T, nu ):

    global gama, p0, w_before, mu, R, nu0, Q, T0

    return np.sum(w_before / mu) * R / (gama(T) - 1.0) - np.sum(w_before / mu) * R * T * gama_(T) / (gama(T) - 1.0) ** 2 - 0.5 * np.sum(w_before / mu) * R * 1 / nu * (nu0 - nu)

def T__( T ):

    global gama, p0, w_before, mu, R, nu0, Q, T0

    return - 2.0 * np.sum(w_before / mu) * R * gama_(T) / (gama(T) - 1.0) ** 2 - np.sum(w_before / mu) * R * T * gama__(T) / (gama(T) - 1.0) ** 2
    + 2.0 * np.sum(w_before / mu) * R * T * gama__(T) ** 2 / (gama(T) - 1.0) ** 3

def iteration_min( p, nu, T ):

    global gama, p0, w_before, mu, R, nu0, Q, T0, T_

    #Tk_1 = p * nu / ( np.sum( w_before / mu ) * R )

    Tk_1 = ( np.sum( w_before / mu ) * R * T / ( gama( T ) - 1.0 ) - np.sum( w_before / mu ) * R * T0 / ( gama( T0 ) - 1.0 )
           - 0.5 * ( np.sum( w_before / mu ) * 1.0 / nu * R * T - np.sum( w_before / mu ) * 1.0 / nu0 * R * T0 ) * ( nu0 - nu ) + Q )\
           / T_( T, nu ) - T

    pk_1 = ( gama( Tk_1 ) - 1.0 ) / nu * ( p0 * nu0 / ( gama( T0 ) - 1.0 ) + Q + 1 / 2.0 * ( p0 + p ) * ( nu0 - nu ) )

    nuk_1 = ( p0 * nu + gama( Tk_1 ) * p * nu0 ) / ( p * ( 1 + gama( Tk_1 ) ) )

    return pk_1, nuk_1, Tk_1

def iteration_max( p, nu, T ):

    global gama, p0, w_before, mu, R, nu0, Q, T0

    pk_mass = [p0, ]
    Tk_mass = [T0, ]
    nuk_mass = [nu0, ]

    pk_1 = 10.0 ** 6

    Tk_1 = 10.0

    nuk_1 = (p0 * nu + gama(Tk_1) * p * nu0) / (p * (1 + gama(Tk_1)))

    k = 1

    # Tk_1 = p * nu / ( np.sum( w_before / mu ) * R )

    while (abs((float(p) - float(pk_1)) / float(pk_1)) >= 0.001 and T__( Tk_1 ) < 0.0 and T_( Tk_1, nuk_1 ) ):

        k += 1

        p = pk_1

        T = Tk_1

        Tk_1 = (np.sum(w_before / mu) * R * T / (gama(T) - 1.0) - np.sum(w_before / mu) * R * T0 / (gama(T0) - 1.0)
            - 0.5 * (np.sum(w_before / mu) * 1.0 / nu * R * T - np.sum(w_before / mu) * 1.0 / nu0 * R * T0) * (nu0 - nu) + Q) \
           / (np.sum(w_before / mu) * R / (gama(T) - 1.0) - np.sum(w_before / mu) * R * gama_(T) / (gama(T) - 1.0) ** 2 \
              + 0.5 * np.sum(w_before / mu) * R * 1 / nu * (nu0 - nu)) - T

        pk_1 = (gama(Tk_1) - 1.0) / nu * (p0 * nu0 / (gama(T0) - 1.0) + Q + 1 / 2.0 * (p0 + p) * (nu0 - nu))

        nuk_1 = (p0 * nu + gama(Tk_1) * p * nu0) / (p * (1 + gama(Tk_1)))

        pk_mass.append(pk_1)
        Tk_mass.append(Tk_1)
        nuk_mass.append(nuk_1)

    return pk_mass, nuk_mass, Tk_mass, k


def D( p, nu ):

    global p0, nu0

    return nu0 * np.sqrt( ( p - p0 ) / ( nu0 - nu ) )

def v( p, nu ):

    global p0, nu0

    return ( nu0 - nu ) * np.sqrt( ( p - p0 ) / ( nu0 - nu ) )

def gugonio( gamma, nu ):

    global nu0, p0

    return p0 * ( ( gamma + 1.0 ) * nu0 - ( gamma - 1.0 ) * nu ) / ( ( gamma + 1.0 ) * nu - ( gamma - 1.0 ) * nu0 )

def draw( name, size ):

    plt.figure( figsize = size )
    plt.grid()
    plt.xlabel( name[ 0 ] )
    plt.ylabel( name[ 1 ] )


