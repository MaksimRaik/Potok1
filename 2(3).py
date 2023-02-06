import numpy as np
import matplotlib.pyplot as plt
import function as f

p, nu, T = 0.0, 0.0, 0.0

p_mass = [f.p0, ]
T_mass = [f.T0, ]
nu_mass = [f.nu0, ]

f.draw( [ 'number of iteration', 'p, atm' ], ( 12, 6 ) )

plt.plot( 0.0, f.p0/ 10 ** 5, 'ro' )

print(r' начальные значения = ', f.p0 / 10 ** 5, 'Па,' , f.nu0, 'м^3 / кг,' , f.T0, 'K' )
print(r' начальные значения = ', f.p01 / 10 ** 5, 'Па,' , f.nu01, 'м^3 / кг,' , f.T01, 'K' )

p, nu, T = f.iteration_min( f.p0, f.nu0, f.T0 )

p_mass.append( p )
T_mass.append( T )
nu_mass.append( nu )

p1 = f.p0

k = 1

while(  abs( ( float( p ) - float( p1 ) )  / float( p1 ) )  >= 0.001 ):

    p1 = p

    p, nu, T = f.iteration_min(p, nu, T)

    p_mass.append( p )
    T_mass.append( T )
    nu_mass.append( nu )

    k += 1

plt.plot(range(k + 1), np.asarray(p_mass) / 10 ** 5, 'bo-', label = r'детонационная волна( невязка )')

Pm, num, Tm, k1 = f.iteration_max(f.p01, f.nu01, f.T01)

plt.plot( range( k1 ), np.asarray( Pm ) / 10 ** 5, 'go-', label = r'дефлаграционная волна( невязка )' )

legend = plt.legend( loc = 'lower right', shadow = True, fontsize = 'x-large' )

plt.show()





f.draw( [ r'$ \eta $', 'p, atm' ], ( 12, 6 ) )

x = ( f.gama(  T_mass[ -1 ]  ) - 1.0 ) * f.nu0 / ( f.gama(  T_mass[ -1 ]  ) + 1.0 ) + 0.01

#for i in range( len(  T_mass ) ):

    #print( f.gugonio( f.gama(  T_mass[i] ), np.linspace( x , 0.9, 150 ) ), i )
plt.plot( np.linspace( x , 0.9, 150 ), f.gugonio( f.gama(  T_mass[ -1 ]  ), np.linspace( x, 0.9, 150 )) /10 ** 5, 'b-', label = r'детонационная волна кривая Гюгона' )
plt.plot( [ nu_mass[ 0 ], nu_mass[ -1 ] ], f.gugonio( f.gama(  T_mass[ -1 ]  ), np.asarray( [nu_mass[ 0 ], nu_mass[ -1 ] ] ) ) / 10 ** 5, 'bo' )

legend = plt.legend( loc = 'upper right', shadow = True, fontsize = 'x-large' )

plt.show()

f.draw( [ r'$ \eta $', 'p, atm' ], ( 12, 6 ) )

plt.plot( np.linspace( x , 0.9, 150 ), f.gugonio( f.gama(  Tm[ -1 ]  ), np.linspace( x, 0.9, 150 ) ) / 10 ** 5 , 'g-', label = r'дефлаграционная волна кривая Гюгона' )
plt.plot( [ num[ 0 ], num[ -1 ] ], f.gugonio( f.gama(  Tm[ -1 ]  ), np.asarray( [num[ 0 ], num[ -1 ] ] ) ) / 10 ** 5, 'go' )

legend = plt.legend( loc = 'upper right', shadow = True, fontsize = 'x-large' )

plt.show()

print( 'D =', np.around( f.D( p, nu ), 2), 'm/c ', 'v =', np.around( f.v( p, nu ), 2), 'm/c',
       'p=', p / 10 ** 5, 'ro=', 1 / nu, 'T=', T, 'gamma=',  f.gama(  T_mass[ -1 ] ) )

print( 'D =', np.around( f.D( Pm[ -1 ], num[ -1 ] ), 2), 'm/c ', 'v =', np.around( f.v( Pm[ -1 ], num[ -1 ] ), 2), 'm/c',
       'p=', Pm[-1] / 10 ** 5, 'ro=', 1 / num[ -1 ], 'T=', Tm[ -1 ], 'gamma=',  f.gama(  Tm[ -1 ] ) )



