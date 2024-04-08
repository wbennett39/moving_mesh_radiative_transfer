from numba import njit
import math
@njit

#if there are issues with accuracy later on, look at turning the decimal numbers rational


def GMatrix(i, j, a, b, ap, bp):
    if i == 0:
        if j == 0:

            return -0.16666666666666666*((a**2 + a*b + b**2)*(ap - bp))/(a - b)

        elif j == 1:

            return ((a + b)*(ap - bp))/(6.*math.sqrt(2))

        elif j == 2:

            return ((a**2 + 3*a*b + b**2)*(ap - bp))/(15.*math.sqrt(2)*(a - b))

        elif j == 3:

            return -0.1*((a + b)*(ap - bp))/math.sqrt(2)

        elif j == 4:

            return ((5*a**2 - 3*a*b + 5*b**2)*(ap - bp))/(105.*math.sqrt(2)*(a - b))

        elif j == 5:

            return ((a + b)*(-ap + bp))/(42.*math.sqrt(2))

        elif j == 6:

            return ((5*a**2 - a*b + 5*b**2)*(ap - bp))/(315.*math.sqrt(2)*(a - b))

        elif j == 7:

            return ((a + b)*(-ap + bp))/(90.*math.sqrt(2))

        elif j == 8:

            return ((29*a**2 - 3*a*b + 29*b**2)*(ap - bp))/(3465.*math.sqrt(2)*(a - b))

        elif j == 9:

            return ((a + b)*(-ap + bp))/(154.*math.sqrt(2))

        elif j == 10:

            return ((47*a**2 - 3*a*b + 47*b**2)*(ap - bp))/(9009.*math.sqrt(2)*(a - b))

    elif i == 1:

        if j == 0:

            return (4*a*b*(ap + bp) + a**2*(7*ap + bp) + b**2*(ap + 7*bp))/(6.*math.sqrt(2)*(a - b))

        elif j == 1:

            return (-(b**2*(ap - 11*bp)) + a**2*(-11*ap + bp) + 3*a*b*(-ap + bp))/(15.*(a - b))

        elif j == 2:

            return -0.03333333333333333*(12*a*b*(ap + bp) + b**2*(7*ap + bp) + a**2*(ap + 7*bp))/(a - b)

        elif j == 3:

            return (9*a*b*(ap - bp) + a**2*(13*ap + bp) - b**2*(ap + 13*bp))/(35.*(a - b))

        elif j == 4:

            return (b**2*(19*ap - 59*bp) + 12*a*b*(ap + bp) + a**2*(-59*ap + 19*bp))/(210.*(a - b))

        elif j == 5:

            return (b**2*(ap - 3*bp) + a**2*(3*ap - bp) + a*b*(-ap + bp))/(21.*(a - b))

        elif j == 6:

            return (b**2*(13*ap - 53*bp) + 4*a*b*(ap + bp) + a**2*(-53*ap + 13*bp))/(630.*(a - b))

        elif j == 7:

            return (b**2*(7*ap - 29*bp) + a**2*(29*ap - 7*bp) + 3*a*b*(-ap + bp))/(495.*(a - b))

        elif j == 8:

            return (b**2*(67*ap - 299*bp) + 12*a*b*(ap + bp) + a**2*(-299*ap + 67*bp))/(6930.*(a - b))

        elif j == 9:

            return (b**2*(37*ap - 167*bp) + a**2*(167*ap - 37*bp) + 9*a*b*(-ap + bp))/(5005.*(a - b))

        elif j == 10:

            return (b**2*(103*ap - 479*bp) + 12*a*b*(ap + bp) + a**2*(-479*ap + 103*bp))/(18018.*(a - b))

    elif i == 2:

        if j == 0:

            return (a*b*(-ap + bp) - a**2*(7*ap + bp) + b**2*(ap + 7*bp))/(3.*math.sqrt(2)*(a - b))

        elif j == 1:

            return (16*a*b*(ap + bp) + a**2*(57*ap + 7*bp) + b**2*(7*ap + 57*bp))/(30.*(a - b))

        elif j == 2:

            return (11*a*b*(-ap + bp) + a**2*(-75*ap + 19*bp) + b**2*(-19*ap + 75*bp))/(105.*(a - b))

        elif j == 3:

            return -0.014285714285714285*(48*a*b*(ap + bp) + a**2*(35*ap + 29*bp) + b**2*(29*ap + 35*bp))/(a - b)

        elif j == 4:

            return (119*a*b*(ap - bp) + a**2*(281*ap + 31*bp) - b**2*(31*ap + 281*bp))/(315.*(a - b))

        elif j == 5:

            return (b**2*(15*ap - 79*bp) + 16*a*b*(ap + bp) + a**2*(-79*ap + 15*bp))/(126.*(a - b))

        elif j == 6:

            return (b**2*(221*ap - 1189*bp) + a**2*(1189*ap - 221*bp) + 299*a*b*(-ap + bp))/(3465.*(a - b))

        elif j == 7:

            return (b**2*(23*ap - 215*bp) + 16*a*b*(ap + bp) + a**2*(-215*ap + 23*bp))/(990.*(a - b))

        elif j == 8:

            return (b**2*(711*ap - 7055*bp) + a**2*(7055*ap - 711*bp) + 551*a*b*(-ap + bp))/(45045.*(a - b))

        elif j == 9:

            return (b**2*(101*ap - 1189*bp) + 48*a*b*(ap + bp) + a**2*(-1189*ap + 101*bp))/(10010.*(a - b))

        elif j == 10:

            return (b**2*(69*ap - 845*bp) + a**2*(845*ap - 69*bp) + 35*a*b*(-ap + bp))/(9009.*(a - b))

    elif i == 3:

        if j == 0:

            return (-(b**2*(ap - 25*bp)) + a**2*(25*ap - bp) - 4*a*b*(ap + bp))/(10.*math.sqrt(2)*(a - b))

        elif j == 1:

            return (a*b*(-ap + bp) - a**2*(13*ap + bp) + b**2*(ap + 13*bp))/(5.*(a - b))

        elif j == 2:

            return (68*a*b*(ap + bp) + a**2*(155*ap + 29*bp) + b**2*(29*ap + 155*bp))/(70.*(a - b))

        elif j == 3:

            return (47*a*b*(-ap + bp) + a**2*(-215*ap + 53*bp) + b**2*(-53*ap + 215*bp))/(315.*(a - b))

        elif j == 4:

            return -0.0015873015873015873*(636*a*b*(ap + bp) + a**2*(617*ap + 367*bp) + b**2*(367*ap + 617*bp))/(a - b)

        elif j == 5:

            return (351*a*b*(ap - bp) + a**2*(1007*ap + 115*bp) - b**2*(115*ap + 1007*bp))/(693.*(a - b))

        elif j == 6:

            return (b**2*(349*ap - 2341*bp) + 452*a*b*(ap + bp) + a**2*(-2341*ap + 349*bp))/(2310.*(a - b))

        elif j == 7:

            return (b**2*(527*ap - 3725*bp) + a**2*(3725*ap - 527*bp) + 807*a*b*(-ap + bp))/(6435.*(a - b))

        elif j == 8:

            return (b**2*(2447*ap - 34295*bp) + 2364*a*b*(ap + bp) + a**2*(-34295*ap + 2447*bp))/(90090.*(a - b))

        elif j == 9:

            return (b**2*(275*ap - 4217*bp) + a**2*(4217*ap - 275*bp) + 283*a*b*(-ap + bp))/(15015.*(a - b))

        elif j == 10:

            return (b**2*(1003*ap - 19555*bp) + 732*a*b*(ap + bp) + a**2*(-19555*ap + 1003*bp))/(90090.*(a - b))

    elif i == 4:

        if j == 0:

            return (29*a*b*(ap - bp) + a**2*(-235*ap + 11*bp) + b**2*(-11*ap + 235*bp))/(105.*math.sqrt(2)*(a - b))

        elif j == 1:

            return (a**2*(531*ap - 19*bp) - 64*a*b*(ap + bp) + b**2*(-19*ap + 531*bp))/(210.*(a - b))

        elif j == 2:

            return (43*a*b*(-ap + bp) - a**2*(321*ap + 31*bp) + b**2*(31*ap + 321*bp))/(105.*(a - b))

        elif j == 3:

            return (832*a*b*(ap + bp) + a**2*(1681*ap + 367*bp) + b**2*(367*ap + 1681*bp))/(630.*(a - b))

        elif j == 4:

            return (547*a*b*(-ap + bp) + a**2*(-2339*ap + 579*bp) + b**2*(-579*ap + 2339*bp))/(3465.*(a - b))

        elif j == 5:

            return -0.0007215007215007215*(1856*a*b*(ap + bp) + a**2*(2033*ap + 1039*bp) + b**2*(1039*ap + 2033*bp))/(a - b)

        elif j == 6:

            return (4103*a*b*(ap - bp) + a**2*(13061*ap + 1499*bp) - b**2*(1499*ap + 13061*bp))/(6435.*(a - b))

        elif j == 7:

            return (b**2*(2367*ap - 18239*bp) + 3392*a*b*(ap + bp) + a**2*(-18239*ap + 2367*bp))/(12870.*(a - b))

        elif j == 8:

            return (b**2*(647*ap - 5351*bp) + a**2*(5351*ap - 647*bp) - 1057*a*b*(ap - bp))/(6435.*(a - b))

        elif j == 9:

            return (b**2*(947*ap - 16819*bp) + 1088*a*b*(ap + bp) + a**2*(-16819*ap + 947*bp))/(30030.*(a - b))

        elif j == 10:

            return (b**2*(16211*ap - 321395*bp) + a**2*(321395*ap - 16211*bp) - 19519*a*b*(ap - bp))/(765765.*(a - b))

    elif i == 5:

        if j == 0:

            return (-(b**2*(ap - 89*bp)) + a**2*(89*ap - bp) - 4*a*b*(ap + bp))/(42.*math.sqrt(2)*(a - b))

        elif j == 1:

            return (13*a*b*(ap - bp) + a**2*(-143*ap + 5*bp) + b**2*(-5*ap + 143*bp))/(63.*(a - b))

        elif j == 2:

            return (a**2*(359*ap - 15*bp) - 44*a*b*(ap + bp) + b**2*(-15*ap + 359*bp))/(126.*(a - b))

        elif j == 3:

            return (35*a*b*(-ap + bp) - a**2*(227*ap + 23*bp) + b**2*(23*ap + 227*bp))/(63.*(a - b))

        elif j == 4:

            return (2300*a*b*(ap + bp) + a**2*(4361*ap + 1039*bp) + b**2*(1039*ap + 4361*bp))/(1386.*(a - b))

        elif j == 5:

            return (-1453*a*b*(ap - bp) + a**2*(-6053*ap + 1503*bp) + b**2*(-1503*ap + 6053*bp))/(9009.*(a - b))

        elif j == 6:

            return -0.0003885003885003885*(4300*a*b*(ap + bp) + a**2*(5041*ap + 2359*bp) + b**2*(2359*ap + 5041*bp))/(a - b)

        elif j == 7:

            return (1651*a*b*(ap - bp) + a**2*(5607*ap + 643*bp) - b**2*(643*ap + 5607*bp))/(2145.*(a - b))

        elif j == 8:

            return (b**2*(931*ap - 7851*bp) + 1420*a*b*(ap + bp) + a**2*(-7851*ap + 931*bp))/(4290.*(a - b))

        elif j == 9:

            return (5*b**2*(6095*ap - 55837*bp) + 5*a**2*(55837*ap - 6095*bp) - 51787*a*b*(ap - bp))/(255255.*(a - b))

        elif j == 10:

            return (b**2*(335*ap - 6951*bp) + 428*a*b*(ap + bp) + a**2*(-6951*ap + 335*bp))/(9282.*(a - b))

    elif i == 6:

        if j == 0:

            return (23*a*b*(ap - bp) + a**2*(-655*ap + 7*bp) + b**2*(-7*ap + 655*bp))/(315.*math.sqrt(2)*(a - b))

        elif j == 1:

            return (a**2*(1357*ap - 13*bp) - 48*a*b*(ap + bp) + b**2*(-13*ap + 1357*bp))/(630.*(a - b))

        elif j == 2:

            return (789*a*b*(ap - bp) + a**2*(-8507*ap + 323*bp) + b**2*(-323*ap + 8507*bp))/(3465.*(a - b))

        elif j == 3:

            return (a**2*(7453*ap - 349*bp) - 944*a*b*(ap + bp) + b**2*(-349*ap + 7453*bp))/(2310.*(a - b))

        elif j == 4:

            return (-2403*a*b*(ap - bp) - a**2*(14461*ap + 1499*bp) + b**2*(1499*ap + 14461*bp))/(3465.*(a - b))

        elif j == 5:

            return (5136*a*b*(ap + bp) + a**2*(9353*ap + 2359*bp) + b**2*(2359*ap + 9353*bp))/(2574.*(a - b))

        elif j == 6:

            return (-1049*a*b*(ap - bp) + a**2*(-4313*ap + 1073*bp) + b**2*(-1073*ap + 4313*bp))/(6435.*(a - b))

        elif j == 7:

            return -0.0006993006993006993*(2864*a*b*(ap + bp) + a**2*(3507*ap + 1549*bp) + b**2*(1549*ap + 3507*bp))/(a - b)

        elif j == 8:

            return (997*a*b*(ap - bp) + a**2*(3539*ap + 405*bp) - b**2*(405*ap + 3539*bp))/(1105.*(a - b))

        elif j == 9:

            return (b**2*(553*ap - 4969*bp) + 880*a*b*(ap + bp) + a**2*(-4969*ap + 553*bp))/(2210.*(a - b))

        elif j == 10:

            return (b**2*(20329*ap - 200145*bp) + a**2*(200145*ap - 20329*bp) - 35471*a*b*(ap - bp))/(146965.*(a - b))

    elif i == 7:

        if j == 0:

            return (-(b**2*(ap - 185*bp)) + a**2*(185*ap - bp) - 4*a*b*(ap + bp))/(90.*math.sqrt(2)*(a - b))

        elif j == 1:

            return (29*a*b*(ap - bp) + a**2*(-1043*ap + 9*bp) + b**2*(-9*ap + 1043*bp))/(495.*(a - b))

        elif j == 2:

            return (a**2*(2255*ap - 23*bp) - 76*a*b*(ap + bp) + b**2*(-23*ap + 2255*bp))/(990.*(a - b))

        elif j == 3:

            return (1673*a*b*(ap - bp) + a**2*(-17275*ap + 713*bp) + b**2*(-713*ap + 17275*bp))/(6435.*(a - b))

        elif j == 4:

            return (a**2*(46663*ap - 2367*bp) - 6076*a*b*(ap + bp) + b**2*(-2367*ap + 46663*bp))/(12870.*(a - b))

        elif j == 5:

            return (-1067*a*b*(ap - bp) - a**2*(6119*ap + 643*bp) + b**2*(643*ap + 6119*bp))/(1287.*(a - b))

        elif j == 6:

            return (3332*a*b*(ap + bp) + a**2*(5899*ap + 1549*bp) + b**2*(1549*ap + 5899*bp))/(1430.*(a - b))

        elif j == 7:

            return (-5981*a*b*(ap - bp) + a**2*(-24405*ap + 6079*bp) + b**2*(-6079*ap + 24405*bp))/(36465.*(a - b))

        elif j == 8:

            return -0.00015082956259426848*(15484*a*b*(ap + bp) + a**2*(19545*ap + 8287*bp) + b**2*(8287*ap + 19545*bp))/(a - b)

        elif j == 9:

            return (5015*a*b*(ap - bp) + a**2*(18383*ap + 2099*bp) - b**2*(2099*ap + 18383*bp))/(4845.*(a - b))

        elif j == 10:

            return (b**2*(2747*ap - 25875*bp) + 4508*a*b*(ap + bp) + a**2*(-25875*ap + 2747*bp))/(9690.*(a - b))

    elif i == 8:

        if j == 0:

            return (25*a*b*(ap - bp) + a**2*(-1415*ap + 7*bp) + b**2*(-7*ap + 1415*bp))/(693.*math.sqrt(2)*(a - b))

        elif j == 1:

            return (a**2*(14403*ap - 67*bp) - 256*a*b*(ap + bp) + b**2*(-67*ap + 14403*bp))/(6930.*(a - b))

        elif j == 2:

            return (2609*a*b*(ap - bp) + a**2*(-99045*ap + 869*bp) + b**2*(-869*ap + 99045*bp))/(45045.*(a - b))

        elif j == 3:

            return -0.0000111000111000111*(b**2*(2447*ap - 219535*bp) + 7424*a*b*(ap + bp) + a**2*(-219535*ap + 2447*bp))/(a - b)

        elif j == 4:

            return (4435*a*b*(ap - bp) + a**2*(-44053*ap + 1941*bp) + b**2*(-1941*ap + 44053*bp))/(15015.*(a - b))

        elif j == 5:

            return (a**2*(17315*ap - 931*bp) - 2304*a*b*(ap + bp) + b**2*(-931*ap + 17315*bp))/(4290.*(a - b))

        elif j == 6:

            return (689*a*b*(-ap + bp) - 3*a**2*(1273*ap + 135*bp) + 3*b**2*(135*ap + 1273*bp))/(715.*(a - b))

        elif j == 7:

            return (17664*a*b*(ap + bp) + a**2*(30625*ap + 8287*bp) + b**2*(8287*ap + 30625*bp))/(6630.*(a - b))

        elif j == 8:

            return (-10371*a*b*(ap - bp) + a**2*(-42115*ap + 10499*bp) + b**2*(-10499*ap + 42115*bp))/(62985.*(a - b))

        elif j == 9:

            return -0.00010319917440660474*(25856*a*b*(ap + bp) + a**2*(33377*ap + 13727*bp) + b**2*(13727*ap + 33377*bp))/(a - b)

        elif j == 10:

            return (7923*a*b*(ap - bp) + a**2*(29761*ap + 3391*bp) - b**2*(3391*ap + 29761*bp))/(6783.*(a - b))

    elif i == 9:

        if j == 0:

            return (-(b**2*(ap - 313*bp)) + a**2*(313*ap - bp) - 4*a*b*(ap + bp))/(154.*math.sqrt(2)*(a - b))

        elif j == 1:

            return (151*a*b*(ap - bp) + a**2*(-10313*ap + 43*bp) + b**2*(-43*ap + 10313*bp))/(5005.*(a - b))

        elif j == 2:

            return (a**2*(21517*ap - 101*bp) - 356*a*b*(ap + bp) + b**2*(-101*ap + 21517*bp))/(10010.*(a - b))

        elif j == 3:

            return (917*a*b*(ap - bp) + a**2*(-34831*ap + 325*bp) + b**2*(-325*ap + 34831*bp))/(15015.*(a - b))

        elif j == 4:

            return (a**2*(78491*ap - 947*bp) - 2700*a*b*(ap + bp) + b**2*(-947*ap + 78491*bp))/(30030.*(a - b))

        elif j == 5:

            return (84741*a*b*(ap - bp) + a**2*(-815335*ap + 37789*bp) + b**2*(-37789*ap + 815335*bp))/(255255.*(a - b))

        elif j == 6:

            return (a**2*(9841*ap - 553*bp) - 1332*a*b*(ap + bp) + b**2*(-553*ap + 9841*bp))/(2210.*(a - b))

        elif j == 7:

            return (-3639*a*b*(ap - bp) - a**2*(19663*ap + 2099*bp) + b**2*(2099*ap + 19663*bp))/(3315.*(a - b))

        elif j == 8:

            return (29052*a*b*(ap + bp) + a**2*(49561*ap + 13727*bp) + b**2*(13727*ap + 49561*bp))/(9690.*(a - b))

        elif j == 9:

            return (-5599*a*b*(ap - bp) + a**2*(-22663*ap + 5653*bp) + b**2*(-5653*ap + 22663*bp))/(33915.*(a - b))

        elif j == 10:

            return -0.00007371369600471768*(40716*a*b*(ap + bp) + a**2*(53473*ap + 21479*bp) + b**2*(21479*ap + 53473*bp))/(a - b)

    elif i == 10:

        if j == 0:

            return (197*a*b*(ap - bp) + a**2*(-18253*ap + 53*bp) + b**2*(-53*ap + 18253*bp))/(9009.*math.sqrt(2)*(a - b))

        elif j == 1:

            return (a**2*(36903*ap - 103*bp) - 400*a*b*(ap + bp) + b**2*(-103*ap + 36903*bp))/(18018.*(a - b))

        elif j == 2:

            return (1297*a*b*(ap - bp) + a**2*(-95391*ap + 391*bp) + b**2*(-391*ap + 95391*bp))/(45045.*(a - b))

        elif j == 3:

            return -0.0000111000111000111*(17*b**2*(59*ap - 11899*bp) + 3280*a*b*(ap + bp) + 17*a**2*(-11899*ap + 59*bp))/(a - b)

        elif j == 4:

            return (16771*a*b*(ap - bp) + a**2*(-627075*ap + 6235*bp) + 5*b**2*(-1247*ap + 125415*bp))/(255255.*(a - b))

        elif j == 5:

            return (a**2*(25999*ap - 335*bp) - 912*a*b*(ap + bp) + b**2*(-335*ap + 25999*bp))/(9282.*(a - b))

        elif j == 6:

            return (7751*a*b*(ap - bp) + a**2*(-72665*ap + 3505*bp) + 5*b**2*(-701*ap + 14533*bp))/(20995.*(a - b))

        elif j == 7:

            return (a**2*(47227*ap - 2747*bp) - 6480*a*b*(ap + bp) + b**2*(-2747*ap + 47227*bp))/(9690.*(a - b))

        elif j == 8:

            return (-5967*a*b*(ap - bp) - a**2*(31609*ap + 3391*bp) + b**2*(3391*ap + 31609*bp))/(4845.*(a - b))

        elif j == 9:

            return (45200*a*b*(ap + bp) + a**2*(76121*ap + 21479*bp) + b**2*(21479*ap + 76121*bp))/(13566.*(a - b))

        elif j == 10:

            return (-25803*a*b*(ap - bp) + a**2*(-104203*ap + 26003*bp) + b**2*(-26003*ap + 104203*bp))/(156009.*(a - b))
        
@njit
def MPRIME(i, j, a, b, ap, bp):
    if i == 0:
        if j == 0:
            return (a*(2*ap + bp) + b*(ap + 2*bp))/3.
        elif j == 1:
            return (math.sqrt(2)*(-(a*ap) + b*bp))/3.
        elif j == 2:
            return -0.06666666666666667*(math.sqrt(2)*(2*a*ap + 3*ap*b + 3*a*bp + 2*b*bp))
        elif j == 3:
            return (math.sqrt(2)*(a*ap - b*bp))/5.
        elif j == 4:
            return (math.sqrt(2)*(-10*a*ap + 3*ap*b + 3*a*bp - 10*b*bp))/105.
        elif j == 5:

            return (math.sqrt(2)*(a*ap - b*bp))/21.

        elif j == 6:

            return (math.sqrt(2)*(b*(ap - 10*bp) + a*(-10*ap + bp)))/315.

        elif j == 7:

            return (math.sqrt(2)*(a*ap - b*bp))/45.

        elif j == 8:

            return (math.sqrt(2)*(-58*a*ap + 3*ap*b + 3*a*bp - 58*b*bp))/3465.

        elif j == 9:

            return (math.sqrt(2)*(a*ap - b*bp))/77.

        elif j == 10:

            return (math.sqrt(2)*(-94*a*ap + 3*ap*b + 3*a*bp - 94*b*bp))/9009.

    elif i == 1:

        if j == 0:

            return (math.sqrt(2)*(-(a*ap) + b*bp))/3.

        elif j == 1:

            return (2*(a*(4*ap + bp) + b*(ap + 4*bp)))/15.

        elif j == 2:

            return (-2*(a*ap - b*bp))/15.

        elif j == 3:

            return (-2*(4*a*ap + 3*ap*b + 3*a*bp + 4*b*bp))/35.

        elif j == 4:

            return (26*(a*ap - b*bp))/105.

        elif j == 5:

            return (-2*(4*a*ap - ap*b - a*bp + 4*b*bp))/63.

        elif j == 6:

            return (22*(a*ap - b*bp))/315.

        elif j == 7:

            return (-2*(12*a*ap - ap*b - a*bp + 12*b*bp))/495.

        elif j == 8:

            return (122*(a*ap - b*bp))/3465.

        elif j == 9:

            return (-2*(68*a*ap - 3*ap*b - 3*a*bp + 68*b*bp))/5005.

        elif j == 10:

            return (194*(a*ap - b*bp))/9009.

    elif i == 2:

        if j == 0:

            return -0.06666666666666667*(math.sqrt(2)*(2*a*ap + 3*ap*b + 3*a*bp + 2*b*bp))

        elif j == 1:

            return (-2*(a*ap - b*bp))/15.

        elif j == 2:

            return (2*(30*a*ap + 19*ap*b + 19*a*bp + 30*b*bp))/105.

        elif j == 3:

            return (2*(-(a*ap) + b*bp))/7.

        elif j == 4:

            return (-2*(26*a*ap + 31*ap*b + 31*a*bp + 26*b*bp))/315.

        elif j == 5:

            return (2*(a*ap - b*bp))/9.

        elif j == 6:

            return (-2*(194*a*ap - 51*ap*b - 51*a*bp + 194*b*bp))/3465.

        elif j == 7:

            return (2*(a*ap - b*bp))/33.

        elif j == 8:

            return (-2*(950*a*ap - 79*ap*b - 79*a*bp + 950*b*bp))/45045.

        elif j == 9:

            return (2*(a*ap - b*bp))/65.

        elif j == 10:

            return (-2*(538*a*ap - 23*ap*b - 23*a*bp + 538*b*bp))/45045.

    elif i == 3:

        if j == 0:

            return (math.sqrt(2)*(a*ap - b*bp))/5.

        elif j == 1:

            return (-2*(4*a*ap + 3*ap*b + 3*a*bp + 4*b*bp))/35.

        elif j == 2:

            return (2*(-(a*ap) + b*bp))/7.

        elif j == 3:

            return (2*(100*a*ap + 53*ap*b + 53*a*bp + 100*b*bp))/315.

        elif j == 4:

            return (-14*(a*ap - b*bp))/45.

        elif j == 5:

            return (-2*(52*a*ap + 69*ap*b + 69*a*bp + 52*b*bp))/693.

        elif j == 6:

            return (82*(a*ap - b*bp))/385.

        elif j == 7:

            return (-2*(340*a*ap - 93*ap*b - 93*a*bp + 340*b*bp))/6435.

        elif j == 8:

            return (46*(a*ap - b*bp))/819.

        elif j == 9:

            return (-584*a*ap + 50*ap*b + 50*a*bp - 584*b*bp)/15015.

        elif j == 10:

            return (14*(a*ap - b*bp))/495.

    elif i == 4:

        if j == 0:

            return (math.sqrt(2)*(-10*a*ap + 3*ap*b + 3*a*bp - 10*b*bp))/105.

        elif j == 1:

            return (26*(a*ap - b*bp))/105.

        elif j == 2:

            return (-2*(26*a*ap + 31*ap*b + 31*a*bp + 26*b*bp))/315.

        elif j == 3:

            return (-14*(a*ap - b*bp))/45.

        elif j == 4:

            return (2*(1126*a*ap + 579*ap*b + 579*a*bp + 1126*b*bp))/3465.

        elif j == 5:

            return (-74*(a*ap - b*bp))/231.

        elif j == 6:

            return (-2*(3238*a*ap + 4497*ap*b + 4497*a*bp + 3238*b*bp))/45045.

        elif j == 7:

            return (122*(a*ap - b*bp))/585.

        elif j == 8:

            return (-2*(2306*a*ap - 647*ap*b - 647*a*bp + 2306*b*bp))/45045.

        elif j == 9:

            return (62*(a*ap - b*bp))/1155.

        elif j == 10:

            return (-2*(14150*a*ap - 1247*ap*b - 1247*a*bp + 14150*b*bp))/765765.

    elif i == 5:

        if j == 0:

            return (math.sqrt(2)*(a*ap - b*bp))/21.

        elif j == 1:

            return (-2*(4*a*ap - ap*b - a*bp + 4*b*bp))/63.

        elif j == 2:

            return (2*(a*ap - b*bp))/9.

        elif j == 3:

            return (-2*(52*a*ap + 69*ap*b + 69*a*bp + 52*b*bp))/693.

        elif j == 4:

            return (-74*(a*ap - b*bp))/231.

        elif j == 5:

            return (2*(2956*a*ap + 1503*ap*b + 1503*a*bp + 2956*b*bp))/9009.

        elif j == 6:

            return (-38*(a*ap - b*bp))/117.

        elif j == 7:

            return (-2*(452*a*ap + 643*ap*b + 643*a*bp + 452*b*bp))/6435.

        elif j == 8:

            return (34*(a*ap - b*bp))/165.

        elif j == 9:

            return (-2*(12820*a*ap - 3657*ap*b - 3657*a*bp + 12820*b*bp))/255255.

        elif j == 10:

            return (242*(a*ap - b*bp))/4641.

    elif i == 6:

        if j == 0:

            return (math.sqrt(2)*(b*(ap - 10*bp) + a*(-10*ap + bp)))/315.

        elif j == 1:

            return (22*(a*ap - b*bp))/315.

        elif j == 2:

            return (-2*(194*a*ap - 51*ap*b - 51*a*bp + 194*b*bp))/3465.

        elif j == 3:

            return (82*(a*ap - b*bp))/385.

        elif j == 4:

            return (-2*(3238*a*ap + 4497*ap*b + 4497*a*bp + 3238*b*bp))/45045.

        elif j == 5:

            return (-38*(a*ap - b*bp))/117.

        elif j == 6:

            return (2*(2122*a*ap + 1073*ap*b + 1073*a*bp + 2122*b*bp))/6435.

        elif j == 7:

            return (-18*(a*ap - b*bp))/55.

        elif j == 8:

            return (-2*(842*a*ap + 1215*ap*b + 1215*a*bp + 842*b*bp))/12155.

        elif j == 9:

            return (226*(a*ap - b*bp))/1105.

        elif j == 10:

            return (-6*(2430*a*ap - 701*ap*b - 701*a*bp + 2430*b*bp))/146965.

    elif i == 7:

        if j == 0:

            return (math.sqrt(2)*(a*ap - b*bp))/45.

        elif j == 1:

            return (-2*(12*a*ap - ap*b - a*bp + 12*b*bp))/495.

        elif j == 2:

            return (2*(a*ap - b*bp))/33.

        elif j == 3:

            return (-2*(340*a*ap - 93*ap*b - 93*a*bp + 340*b*bp))/6435.

        elif j == 4:

            return (122*(a*ap - b*bp))/585.

        elif j == 5:

            return (-2*(452*a*ap + 643*ap*b + 643*a*bp + 452*b*bp))/6435.

        elif j == 6:

            return (-18*(a*ap - b*bp))/55.

        elif j == 7:

            return (2*(12060*a*ap + 6079*ap*b + 6079*a*bp + 12060*b*bp))/36465.

        elif j == 8:

            return (-218*(a*ap - b*bp))/663.

        elif j == 9:

            return (-2*(4324*a*ap + 6297*ap*b + 6297*a*bp + 4324*b*bp))/62985.

        elif j == 10:

            return (58*(a*ap - b*bp))/285.

    elif i == 8:

        if j == 0:

            return (math.sqrt(2)*(-58*a*ap + 3*ap*b + 3*a*bp - 58*b*bp))/3465.

        elif j == 1:

            return (122*(a*ap - b*bp))/3465.

        elif j == 2:

            return (-2*(950*a*ap - 79*ap*b - 79*a*bp + 950*b*bp))/45045.

        elif j == 3:

            return (46*(a*ap - b*bp))/819.

        elif j == 4:

            return (-2*(2306*a*ap - 647*ap*b - 647*a*bp + 2306*b*bp))/45045.

        elif j == 5:

            return (34*(a*ap - b*bp))/165.

        elif j == 6:

            return (-2*(842*a*ap + 1215*ap*b + 1215*a*bp + 842*b*bp))/12155.

        elif j == 7:

            return (-218*(a*ap - b*bp))/663.

        elif j == 8:

            return (2*(20870*a*ap + 10499*ap*b + 10499*a*bp + 20870*b*bp))/62985.

        elif j == 9:

            return (-94*(a*ap - b*bp))/285.

        elif j == 10:

            return (-2*(2314*a*ap + 3391*ap*b + 3391*a*bp + 2314*b*bp))/33915.

    elif i == 9:

        if j == 0:

            return (math.sqrt(2)*(a*ap - b*bp))/77.

        elif j == 1:

            return (-2*(68*a*ap - 3*ap*b - 3*a*bp + 68*b*bp))/5005.

        elif j == 2:

            return (2*(a*ap - b*bp))/65.

        elif j == 3:

            return (-584*a*ap + 50*ap*b + 50*a*bp - 584*b*bp)/15015.

        elif j == 4:

            return (62*(a*ap - b*bp))/1155.

        elif j == 5:

            return (-2*(12820*a*ap - 3657*ap*b - 3657*a*bp + 12820*b*bp))/255255.

        elif j == 6:

            return (226*(a*ap - b*bp))/1105.

        elif j == 7:

            return (-2*(4324*a*ap + 6297*ap*b + 6297*a*bp + 4324*b*bp))/62985.

        elif j == 8:

            return (-94*(a*ap - b*bp))/285.

        elif j == 9:

            return (2*(11252*a*ap + 5653*ap*b + 5653*a*bp + 11252*b*bp))/33915.

        elif j == 10:

            return (59*(-2*a*ap + 2*b*bp))/357.

    elif i == 10:

        if j == 0:

            return (math.sqrt(2)*(-94*a*ap + 3*ap*b + 3*a*bp - 94*b*bp))/9009.

        elif j == 1:

            return (194*(a*ap - b*bp))/9009.

        elif j == 2:

            return (-2*(538*a*ap - 23*ap*b - 23*a*bp + 538*b*bp))/45045.

        elif j == 3:

            return (14*(a*ap - b*bp))/495.

        elif j == 4:

            return (-2*(14150*a*ap - 1247*ap*b - 1247*a*bp + 14150*b*bp))/765765.

        elif j == 5:

            return (242*(a*ap - b*bp))/4641.

        elif j == 6:

            return (-6*(2430*a*ap - 701*ap*b - 701*a*bp + 2430*b*bp))/146965.

        elif j == 7:

            return (58*(a*ap - b*bp))/285.

        elif j == 8:

            return (-2*(2314*a*ap + 3391*ap*b + 3391*a*bp + 2314*b*bp))/33915.

        elif j == 9:

            return (59*(-2*a*ap + 2*b*bp))/357.

        elif j == 10:

            return (2*(51806*a*ap + 26003*ap*b + 26003*a*bp + 51806*b*bp))/156009.