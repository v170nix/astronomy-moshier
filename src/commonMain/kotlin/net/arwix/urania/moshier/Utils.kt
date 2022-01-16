package net.arwix.urania.moshier

import net.arwix.urania.core.*
import net.arwix.urania.core.calendar.JT
import net.arwix.urania.core.calendar.MJD
import net.arwix.urania.core.calendar.times
import net.arwix.urania.core.calendar.toJT
import net.arwix.urania.core.kepler.KeplerElementsObject
import net.arwix.urania.core.kepler.getSimonJ2000KeplerElements
import net.arwix.urania.core.math.ARCSEC_TO_RAD
import net.arwix.urania.core.math.JD_2000
import net.arwix.urania.core.math.JULIAN_DAYS_PER_CENTURY
import net.arwix.urania.core.math.angle.rad
import net.arwix.urania.core.math.mod3600
import kotlin.math.abs
import kotlin.math.cos
import kotlin.math.sin

/**
 * Prepare lookup table of sin and cos ( i*Lj ) for required multiple
 * angles.
 *
 * @param k
 * @param arg
 * @param n
 */
private fun sscc(k: Int, arg: Double, n: Int, ss: Array<DoubleArray>, cc: Array<DoubleArray>) {
    var cv: Double
    var sv: Double
    var s: Double
    val su = sin(arg)
    val cu = cos(arg)
    ss[k][0] = su /* sin(L) */
    cc[k][0] = cu /* cos(L) */
    sv = 2.0 * su * cu
    cv = cu * cu - su * su
    ss[k][1] = sv /* sin(2L) */
    cc[k][1] = cv
    var i = 2
    while (i < n) {
        s = su * cv + cu * sv
        cv = cu * cv - su * sv
        sv = s
        ss[k][i] = sv /* sin( i+1 L ) */
        cc[k][i] = cv
        i++
    }
//
//    for (i in 2 until n) {
//        s = su * cv + cu * sv
//        cv = cu * cv - su * sv
//        sv = s
//        ss[k][i] = sv /* sin( i+1 L ) */
//        cc[k][i] = cv
//        i++
//    }
}

/**
 * Generic program to accumulate sum of trigonometric series in three
 * variables (e.g., longitude, latitude, radius) of the same list of
 * arguments.
 *
 * @param J Julian day.
 * @param arg_tbl
 * @param distance
 * @param lat_tbl
 * @param lon_tbl
 * @param rad_tbl
 * @param max_harmonic
 * @param max_power_of_t
 * @param maxargs
 * @param timescale
 * @param trunclvl
 * @return An array with x, y, z (AU).
 */
internal fun gplan(
    tt: JT,
    arg_tbl: IntArray,
    distance: Double,
    lat_tbl: DoubleArray,
    lon_tbl: DoubleArray,
    rad_tbl: DoubleArray,
    max_harmonic: IntArray,
    max_power_of_t: Int,
    maxargs: Int,
    timescale: Double,
    trunclvl: Double
): DoubleArray {
    var i: Int
    var j: Int
    var k: Int
    var m: Int
    var k1: Int
    var ip: Int
    var np: Int
    var nt: Int
    val p: IntArray
    val pl: DoubleArray
    val pb: DoubleArray
    val pr: DoubleArray
    var su: Double
    var cu: Double
    var sv: Double
    var cv: Double
    val T: Double
    var t: Double
    var sl: Double
    var sb: Double
    var sr: Double
    T = (tt * JULIAN_DAYS_PER_CENTURY) / timescale

    /* From Simon et al (1994) */
    val freqs = doubleArrayOf( /* Arc sec per 10000 Julian years. */
        53810162868.8982, 21066413643.3548, 12959774228.3429, 6890507749.3988, 1092566037.7991, 439960985.5372,
        154248119.3933, 78655032.0744, 52272245.1795
    )
    val phases = doubleArrayOf( /* Arc sec. */
        252.25090552 * 3600.0,
        181.97980085 * 3600.0,
        100.46645683 * 3600.0,
        355.43299958 * 3600.0,
        34.35151874 * 3600.0,
        50.07744430 * 3600.0,
        314.05500511 * 3600.0,
        304.34866548 * 3600.0,
        860492.1546
    )

    /* Calculate sin( i*MM ), etc. for needed multiple angles. */
    val ss = Array(20) { DoubleArray(41) }
    val cc = Array(20) { DoubleArray(41) }
    i = 0
    while (i < 9) {
        if (max_harmonic[i] > 0) {
            sr = ((freqs[i] * T).mod3600() + phases[i]) * ARCSEC_TO_RAD
            sscc(i, sr, max_harmonic[i], ss, cc)
        }
        i++
    }

    /* Point to start of table of arguments. */
    p = arg_tbl

    /* Point to tabulated cosine and sine amplitudes. */
    pl = lon_tbl
    pb = lat_tbl
    pr = rad_tbl
    sl = 0.0
    sb = 0.0
    sr = 0.0
    np = 0
    nt = 0
    cu = 0.0
    var p_index = -1
    var pl_index = -1
    var pb_index = -1
    var pr_index = -1
    while (true) {

        /* argument of sine and cosine */
        /* Number of periodic arguments. */
        p_index++
        np = p[p_index]
        if (np < 0) break
        if (np == 0) { /* It is a polynomial term. */
            p_index++
            nt = p[p_index]
            /* "Longitude" polynomial (phi). */
            pl_index++
            cu = pl[pl_index]
            ip = 0
            while (ip < nt) {
                pl_index++
                cu = cu * T + pl[pl_index]
                ip++
            }
            sl += cu.mod3600()
            /* "Latitude" polynomial (theta). */pb_index++
            cu = pb[pb_index]
            ip = 0
            while (ip < nt) {
                pb_index++
                cu = cu * T + pb[pb_index]
                ip++
            }
            sb += cu
            /* Radius polynomial (psi). */
            pr_index++
            cu = pr[pr_index]
            ip = 0
            while (ip < nt) {
                pr_index++
                cu = cu * T + pr[pr_index]
                ip++
            }
            sr += cu
            continue
        }
        k1 = 0
        cv = 0.0
        sv = 0.0
        ip = 0
        while (ip < np) {

            /* What harmonic. */
            p_index++
            j = p[p_index]
            /* Which planet. */
            p_index++
            m = p[p_index] - 1
            if (j != 0) {
                k = abs(j) - 1
                su = ss[m][k] /* sin(k*angle) */
                if (j < 0) su = -su
                cu = cc[m][k]
                if (k1 == 0) { /* set first angle */
                    sv = su
                    cv = cu
                    k1 = 1
                } else { /* combine angles */
                    t = su * cv + cu * sv
                    cv = cu * cv - su * sv
                    sv = t
                }
            }
            ip++
        }

        /* Highest power of T. */
        p_index++
        nt = p[p_index]
        /* Longitude. */
        pl_index++
        cu = pl[pl_index]
        pl_index++
        su = pl[pl_index]
        ip = 0
        while (ip < nt) {
            pl_index++
            cu = cu * T + pl[pl_index]
            pl_index++
            su = su * T + pl[pl_index]
            ip++
        }
        sl += cu * cv + su * sv
        /* Latitude. */
        pb_index++
        cu = pb[pb_index]
        pb_index++
        su = pb[pb_index]
        ip = 0
        while (ip < nt) {
            pb_index++
            cu = cu * T + pb[pb_index]
            pb_index++
            su = su * T + pb[pb_index]
            ip++
        }
        sb += cu * cv + su * sv
        /* Radius. */
        pr_index++
        cu = pr[pr_index]
        pr_index++
        su = pr[pr_index]
        ip = 0
        while (ip < nt) {
            pr_index++
            cu = cu * T + pr[pr_index]
            pr_index++
            su = su * T + pr[pr_index]
            ip++
        }
        sr += cu * cv + su * sv
    }
    if (distance == 0.0) return doubleArrayOf(
        (ARCSEC_TO_RAD * sl).rad.normalize().value,
        (ARCSEC_TO_RAD * sb).rad.normalize().value,
        (ARCSEC_TO_RAD * sr).rad.normalize().value
    )
    val pobj = DoubleArray(3)
    pobj[0] = ARCSEC_TO_RAD * sl
    pobj[1] = ARCSEC_TO_RAD * sb
    pobj[2] = distance * (1.0 + ARCSEC_TO_RAD * sr)
    val x: Double = pobj[2] * cos(pobj[0]) * cos(pobj[1])
    val y: Double = pobj[2] * sin(pobj[0]) * cos(pobj[1])
    val z: Double = pobj[2] * sin(pobj[1])
    return doubleArrayOf(x, y, z)
}

/**
 * Generic program to accumulate sum of trigonometric series in one
 * variables (e.g., latitude) of the same list of arguments.
 *
 * @param J Julian day.
 * @param arg_tbl
 * @param distance
 * @param lat_tbl
 * @param lon_tbl
 * @param rad_tbl
 * @param max_harmonic
 * @param max_power_of_t
 * @param maxargs
 * @param timescale
 * @param trunclvl
 * @return Latitude (rad).
 */
internal fun g1plan(
    tt: JT, arg_tbl: IntArray, lat_tbl: DoubleArray, max_harmonic: IntArray, maxargs: Int, timescale: Double,
): Double {
    var i: Int
    var j: Int
    var k: Int
    var m: Int
    var k1: Int
    var ip: Int
    var np: Int
    var nt: Int
    var su: Double
    var cu: Double
    var sv: Double
    var cv: Double
    val T: Double
    var t: Double
    var sb: Double
    val args: DoubleArray = meanElements(tt)
    T = (tt * JULIAN_DAYS_PER_CENTURY) / timescale


    /* Calculate sin( i*MM ), etc. for needed multiple angles. */
    val ss = Array(20) { DoubleArray(41) }
    val cc = Array(20) { DoubleArray(41) }
    i = 0
    while (i < maxargs) {
        if (max_harmonic[i] > 0) {
            sscc(i, args[i], max_harmonic[i], ss, cc)
        }
        i++
    }
    sb = 0.0
    np = 0
    nt = 0
    cu = 0.0
    var p_index = -1
    var pb_index = -1
    while (true) {

        /* argument of sine and cosine */
        /* Number of periodic arguments. */
        p_index++
        np = arg_tbl[p_index]
        if (np < 0) break
        if (np == 0) { /* It is a polynomial term. */
            p_index++
            nt = arg_tbl[p_index]
            /* "Latitude" polynomial (theta). */
            pb_index++
            cu = lat_tbl[pb_index]
            ip = 0
            while (ip < nt) {
                pb_index++
                cu = cu * T + lat_tbl[pb_index]
                ip++
            }
            sb += cu
            continue
        }
        k1 = 0
        cv = 0.0
        sv = 0.0
        ip = 0
        while (ip < np) {

            /* What harmonic. */
            p_index++
            j = arg_tbl[p_index]
            /* Which planet. */
            p_index++
            m = arg_tbl[p_index] - 1
            if (j != 0) {
                k = abs(j) - 1
                su = ss[m][k] /* sin(k*angle) */
                if (j < 0) su = -su
                cu = cc[m][k]
                if (k1 == 0) { /* set first angle */
                    sv = su
                    cv = cu
                    k1 = 1
                } else { /* combine angles */
                    t = su * cv + cu * sv
                    cv = cu * cv - su * sv
                    sv = t
                }
            }
            ip++
        }

        /* Highest power of T. */p_index++
        nt = arg_tbl[p_index]
        /* Latitude. */pb_index++
        cu = lat_tbl[pb_index].toDouble()
        pb_index++
        su = lat_tbl[pb_index].toDouble()
        ip = 0
        while (ip < nt) {
            pb_index++
            cu = cu * T + lat_tbl[pb_index]
            pb_index++
            su = su * T + lat_tbl[pb_index]
            ip++
        }
        sb += cu * cv + su * sv
    }
    return ARCSEC_TO_RAD * sb * 0.0001
}

/**
 * Generic program to accumulate sum of trigonometric series in two
 * variables (e.g., longitude, radius) of the same list of arguments.
 *
 * @param J Julian day.
 * @param arg_tbl
 * @param distance
 * @param lat_tbl
 * @param lon_tbl
 * @param rad_tbl
 * @param max_harmonic
 * @param max_power_of_t
 * @param maxargs
 * @param timescale
 * @param trunclvl
 * @return An array with x, y, z (AU).
 */
internal fun g2plan(
    tt: JT,
    arg_tbl: IntArray,
    distance: Double,
    lat_tbl: DoubleArray,
    lon_tbl: DoubleArray,
    rad_tbl: DoubleArray,
    max_harmonic: IntArray,
    max_power_of_t: Int,
    maxargs: Int,
    timescale: Double,
    trunclvl: Double,
    lat: Double
): DoubleArray {
    var i: Int
    var j: Int
    var k: Int
    var m: Int
    var k1: Int
    var ip: Int
    var np: Int
    var nt: Int
    var su: Double
    var cu: Double
    var sv: Double
    var cv: Double
    val T: Double
    var t: Double
    var sl: Double
    var sr: Double

    //TODO https://bitbucket.org/talonsoalbi/jparsec/raw/372f81ebe5e66d530bf89a0ced19ca2fd8570f24/src/main/java/jparsec/ephem/planets/PlanetEphem.java

    T = (tt * JULIAN_DAYS_PER_CENTURY) / timescale
    val args: DoubleArray = meanElements(tt)

    /* Calculate sin( i*MM ), etc. for needed multiple angles. */
    val ss = Array(20) {
        DoubleArray(
            41
        )
    }
    val cc = Array(20) { DoubleArray(41) }
    i = 0
    while (i < maxargs) {
        if (max_harmonic[i] > 0) {
            sscc(i, args[i], max_harmonic[i], ss, cc)
        }
        i++
    }
    sl = 0.0
    sr = 0.0
    np = 0
    nt = 0
    cu = 0.0
    var p_index = -1
    var pl_index = -1
    var pr_index = -1
    while (true) {

        /* argument of sine and cosine */
        /* Number of periodic arguments. */p_index++
        np = arg_tbl[p_index]
        if (np < 0) break
        if (np == 0) { /* It is a polynomial term. */
            p_index++
            nt = arg_tbl[p_index]
            /* "Longitude" polynomial (phi). */pl_index++
            cu = lon_tbl[pl_index].toDouble()
            ip = 0
            while (ip < nt) {
                pl_index++
                cu = cu * T + lon_tbl[pl_index]
                ip++
            }
            sl += cu
            /* Radius polynomial (psi). */pr_index++
            cu = rad_tbl[pr_index].toDouble()
            ip = 0
            while (ip < nt) {
                pr_index++
                cu = cu * T + rad_tbl[pr_index]
                ip++
            }
            sr += cu
            continue
        }
        k1 = 0
        cv = 0.0
        sv = 0.0
        ip = 0
        while (ip < np) {

            /* What harmonic. */p_index++
            j = arg_tbl[p_index]
            /* Which planet. */p_index++
            m = arg_tbl[p_index] - 1
            if (j != 0) {
                k = abs(j) - 1
                su = ss[m][k] /* sin(k*angle) */
                if (j < 0) su = -su
                cu = cc[m][k]
                if (k1 == 0) { /* set first angle */
                    sv = su
                    cv = cu
                    k1 = 1
                } else { /* combine angles */
                    t = su * cv + cu * sv
                    cv = cu * cv - su * sv
                    sv = t
                }
            }
            ip++
        }

        /* Highest power of T. */p_index++
        nt = arg_tbl[p_index]
        /* Longitude. */pl_index++
        cu = lon_tbl[pl_index].toDouble()
        pl_index++
        su = lon_tbl[pl_index].toDouble()
        ip = 0
        while (ip < nt) {
            pl_index++
            cu = cu * T + lon_tbl[pl_index]
            pl_index++
            su = su * T + lon_tbl[pl_index]
            ip++
        }
        sl += cu * cv + su * sv
        /* Radius. */pr_index++
        cu = rad_tbl[pr_index].toDouble()
        pr_index++
        su = rad_tbl[pr_index].toDouble()
        ip = 0
        while (ip < nt) {
            pr_index++
            cu = cu * T + rad_tbl[pr_index]
            pr_index++
            su = su * T + rad_tbl[pr_index]
            ip++
        }
        sr += cu * cv + su * sv
    }
    val pobj = DoubleArray(3)
    sl = sl * 0.0001
    sr = sr * 0.0001
    if (distance == 0.0) return doubleArrayOf(
        ARCSEC_TO_RAD * sl + args[13], //  + PlanetEphem.LP_equinox,
        lat, ARCSEC_TO_RAD * sr
    )
    pobj[0] = ARCSEC_TO_RAD * sl + args[13] //  + PlanetEphem.LP_equinox
    pobj[1] = lat
    pobj[2] = distance * (1.0 + ARCSEC_TO_RAD * sr)
    val x: Double = pobj[2] * cos(pobj[0]) * cos(pobj[1])
    val y: Double = pobj[2] * sin(pobj[0]) * cos(pobj[1])
    val z: Double = pobj[2] * sin(pobj[1])
    return doubleArrayOf(x, y, z)
}

/**
 * Generic program to accumulate sum of trigonometric series in three
 * variables (e.g., longitude, latitude, radius) of the same list of
 * arguments.
 *
 * @param J Julian day.
 * @param arg_tbl
 * @param distance
 * @param lat_tbl
 * @param lon_tbl
 * @param rad_tbl
 * @param max_harmonic
 * @param max_power_of_t
 * @param maxargs
 * @param timescale
 * @param trunclvl
 * @return An array with x, y, z (AU).
 */
internal fun g3plan(
    tt: JT, arg_tbl: IntArray, distance: Double, lat_tbl: DoubleArray, lon_tbl: DoubleArray,
    rad_tbl: DoubleArray, max_harmonic: IntArray, max_power_of_t: Int, maxargs: Int, timescale: Double, trunclvl: Double,
    libration: Boolean
): DoubleArray {
    var i: Int
    var j: Int
    var k: Int
    var m: Int
    var k1: Int
    var ip: Int
    var np: Int
    var nt: Int
    var su: Double
    var cu: Double
    var sv: Double
    var cv: Double
    val T: Double
    var t: Double
    var sl: Double
    var sb: Double
    var sr: Double

    T = (tt * JULIAN_DAYS_PER_CENTURY) / timescale
    val J = T * timescale + JD_2000

    val args: DoubleArray = meanElements(tt)
    if (libration) args[13] -= getPA_precession(tt) // Only librations
//    T = (J - J2000) / timescale

    /* Calculate sin( i*MM ), etc. for needed multiple angles. */
    val ss = Array(20) {
        DoubleArray(
            41
        )
    }
    val cc = Array(20) { DoubleArray(41) }
    i = 0
    while (i < maxargs) {
        if (max_harmonic[i] > 0) {
            sscc(i, args[i], max_harmonic[i], ss, cc)
        }
        i++
    }
    sl = 0.0
    sb = 0.0
    sr = 0.0
    np = 0
    nt = 0
    cu = 0.0
    var p_index = -1
    var pl_index = -1
    var pb_index = -1
    var pr_index = -1
    while (true) {

        /* argument of sine and cosine */
        /* Number of periodic arguments. */p_index++
        np = arg_tbl[p_index]
        if (np < 0) break
        if (np == 0) { /* It is a polynomial term. */
            p_index++
            nt = arg_tbl[p_index]
            /* "Longitude" polynomial (phi). */pl_index++
            cu = lon_tbl[pl_index].toDouble()
            ip = 0
            while (ip < nt) {
                pl_index++
                cu = cu * T + lon_tbl[pl_index]
                ip++
            }
            sl += cu
            /* "Latitude" polynomial (theta). */pb_index++
            cu = lat_tbl[pb_index].toDouble()
            ip = 0
            while (ip < nt) {
                pb_index++
                cu = cu * T + lat_tbl[pb_index]
                ip++
            }
            sb += cu
            /* Radius polynomial (psi). */pr_index++
            cu = rad_tbl[pr_index].toDouble()
            ip = 0
            while (ip < nt) {
                pr_index++
                cu = cu * T + rad_tbl[pr_index]
                ip++
            }
            sr += cu
            continue
        }
        k1 = 0
        cv = 0.0
        sv = 0.0
        ip = 0
        while (ip < np) {

            /* What harmonic. */p_index++
            j = arg_tbl[p_index]
            /* Which planet. */p_index++
            m = arg_tbl[p_index] - 1
            if (j != 0) {
                k = abs(j) - 1
                su = ss[m][k] /* sin(k*angle) */
                if (j < 0) su = -su
                cu = cc[m][k]
                if (k1 == 0) { /* set first angle */
                    sv = su
                    cv = cu
                    k1 = 1
                } else { /* combine angles */
                    t = su * cv + cu * sv
                    cv = cu * cv - su * sv
                    sv = t
                }
            }
            ip++
        }

        /* Highest power of T. */p_index++
        nt = arg_tbl[p_index]
        /* Longitude. */pl_index++
        cu = lon_tbl[pl_index].toDouble()
        pl_index++
        su = lon_tbl[pl_index].toDouble()
        ip = 0
        while (ip < nt) {
            pl_index++
            cu = cu * T + lon_tbl[pl_index]
            pl_index++
            su = su * T + lon_tbl[pl_index]
            ip++
        }
        sl += cu * cv + su * sv
        /* Latitude. */pb_index++
        cu = lat_tbl[pb_index].toDouble()
        pb_index++
        su = lat_tbl[pb_index].toDouble()
        ip = 0
        while (ip < nt) {
            pb_index++
            cu = cu * T + lat_tbl[pb_index]
            pb_index++
            su = su * T + lat_tbl[pb_index]
            ip++
        }
        sb += cu * cv + su * sv
        /* Radius. */pr_index++
        cu = rad_tbl[pr_index].toDouble()
        pr_index++
        su = rad_tbl[pr_index].toDouble()
        ip = 0
        while (ip < nt) {
            pr_index++
            cu = cu * T + rad_tbl[pr_index]
            pr_index++
            su = su * T + rad_tbl[pr_index]
            ip++
        }
        sr += cu * cv + su * sv
    }
    sl = sl * 0.0001
    sb = sb * 0.0001
    sr = sr * 0.0001
    if (distance == 0.0) return doubleArrayOf(
        ARCSEC_TO_RAD * sl + args[2], // + PlanetEphem.Ea_arcsec,
        ARCSEC_TO_RAD * sb,
        ARCSEC_TO_RAD * sr
    )
    val pobj = DoubleArray(3)
    pobj[0] = ARCSEC_TO_RAD * sl + args[2] // + PlanetEphem.Ea_arcsec
    pobj[1] = ARCSEC_TO_RAD * sb
    pobj[2] = distance * (1.0 + ARCSEC_TO_RAD * sr)
    val x: Double = pobj[2] * cos(pobj[0]) * cos(pobj[1])
    val y: Double = pobj[2] * sin(pobj[0]) * cos(pobj[1])
    val z: Double = pobj[2] * sin(pobj[1])
    return doubleArrayOf(x, y, z)
}

/**
 * Obtain mean elements of the planets.
 * @param T Julian centuries
 * @return An array with the mean longitudes.
 */
private fun meanElements(T: JT): DoubleArray {
    val Args = DoubleArray(20)

    /** Mean longitudes of planets (Simon et al, 1994) .047" subtracted from
     * constant term for offset to DE403 origin.
     */

    val delta = (-.047 * ARCSEC_TO_RAD).rad

    Args[0] = (getSimonJ2000KeplerElements(KeplerElementsObject.Mercury).getLongitude(T) + delta).value
    Args[1] = (getSimonJ2000KeplerElements(KeplerElementsObject.Venus).getLongitude(T) + delta).value
    Args[2] = (getSimonJ2000KeplerElements(KeplerElementsObject.Earth).getLongitude(T) + delta).value
    Args[3] = (getSimonJ2000KeplerElements(KeplerElementsObject.Mars).getLongitude(T) + delta).value
    Args[4] = (getSimonJ2000KeplerElements(KeplerElementsObject.Jupiter).getLongitude(T) + delta).value
    Args[5] = (getSimonJ2000KeplerElements(KeplerElementsObject.Saturn).getLongitude(T) + delta).value
    Args[6] = (getSimonJ2000KeplerElements(KeplerElementsObject.Uranus).getLongitude(T) + delta).value
    Args[7] = (getSimonJ2000KeplerElements(KeplerElementsObject.Neptune).getLongitude(T) + delta).value

    /* Copied from cmoon.c, DE404 version. */
    /* Mean elongation of moon = elongation */
    Args[9] = getMeanElongationOfMoon(T)
    Args[10] = getAscendingNode(T)
    Args[11] = getMeanAnomalyOfSun(T)
    Args[12] = getMeanAnomalyOfMoon(T)
    Args[13] = getLongitudeOfMoon(T)
    Args[14] = getLunarFreeLibrations(T)
    Args[15] = getLB(T)
    Args[16] = getLC(T)
    Args[17] = Args[13] - Args[10]
    Args[18] = getNB(Args[17], T)
    return Args
}

private fun getMeanElongationOfMoon(T: JT): Double {
    val T2 = T * T
    /* Mean elongation of moon = D */
    var x = (1.6029616009939659e+09 * T + 1.0722612202445078e+06)
    x += (((((-3.207663637426e-013 * T + 2.555243317839e-011) * T + 2.560078201452e-009) * T - 3.702060118571e-005) * T + 6.9492746836058421e-03) * T /* D, t^3 */
            - 6.7352202374457519e+00) * T2; /* D, t^2 */

    return ARCSEC_TO_RAD * x
}

private fun getAscendingNode(T: JT): Double {
    val T2 = T * T
    /* Mean distance of moon from its ascending node = F */
    var x = (1.7395272628437717e+09 * T + 3.3577951412884740e+05);
    x += (((((4.474984866301e-013 * T + 4.189032191814e-011) * T - 2.790392351314e-009) * T - 2.165750777942e-006) * T - 7.5311878482337989e-04) * T /* F, t^3 */
            - 1.3117809789650071e+01) * T2; /* F, t^2 */
    return ARCSEC_TO_RAD * x
}

private fun getMeanAnomalyOfSun(T: JT): Double {
    /* Mean anomaly of sun = l' (J. Laskar) */
    val T2 = T * T
    var x = (1.2959658102304320e+08 * T + 1.2871027407441526e+06);
    x += ((((((((1.62e-20 * T - 1.0390e-17) * T - 3.83508e-15) * T + 4.237343e-13) * T + 8.8555011e-11) * T - 4.77258489e-8) * T - 1.1297037031e-5) * T + 8.7473717367324703e-05) * T - 5.5281306421783094e-01) * T2;
    return ARCSEC_TO_RAD * x
}

private fun getMeanAnomalyOfMoon(T: JT): Double {
    val T2 = T * T
    var x = (1.7179159228846793e+09 * T + 4.8586817465825332e+05);
    x += (((((-1.755312760154e-012 * T + 3.452144225877e-011) * T - 2.506365935364e-008) * T - 2.536291235258e-004) * T + 5.2099641302735818e-02) * T /* l, t^3 */
            + 3.1501359071894147e+01) * T2; /* l, t^2 */
    return ARCSEC_TO_RAD * x
}

private fun getLongitudeOfMoon(T: JT): Double {
    val T2 = T * T
    /* Mean longitude of moon, re mean ecliptic and equinox of date = L */
    var x = (1.7325643720442266e+09 * T + 7.8593980921052420e+05);
    x += (((((7.200592540556e-014 * T + 2.235210987108e-010) * T - 1.024222633731e-008) * T - 6.073960534117e-005) * T + 6.9017248528380490e-03) * T /* L, t^3 */
            - 5.6550460027471399e+00) * T2; /* L, t^2 */
    return ARCSEC_TO_RAD * x
}

private fun getLunarFreeLibrations(T: JT): Double {
    /* Lunar free librations. */
    /* 74.7 years. Denoted W or LA. */
    val x = (-0.112 * T + 1.73655499e6) * T - 389552.81
    return ARCSEC_TO_RAD * x
}

private fun getLB(T: JT): Double {
    return ARCSEC_TO_RAD * (4.48175409e7 * T + 806045.7);
}

private fun getLC(T: JT): Double {
    return ARCSEC_TO_RAD * (5.36486787e6 * T - 391702.8);
}

private fun getNB(lc: Double, T: JT): Double {
    val x = (((-0.000004 * T + 0.000026) * T + 0.153382) * T - 867.919986) * T + 629543.967373;

    var pA_precession = (((((((((-8.66e-20 * T - 4.759e-17) * T + 2.424e-15) * T + 1.3095e-12) * T + 1.7451e-10) * T - 1.8055e-8) * T - 0.0000235316) * T + 0.000076) * T + 1.105414) * T + 5028.791959) * T;
    /* Moon's longitude re fixed J2000 equinox. */
    pA_precession = ARCSEC_TO_RAD * x;

    return lc + ARCSEC_TO_RAD * (3.24e5 - x) - pA_precession
}

private fun getPA_precession(T: JT): Double {
    val x = (((-0.000004 * T + 0.000026) * T + 0.153382) * T - 867.919986) * T + 629543.967373;

    var pA_precession = (((((((((-8.66e-20 * T - 4.759e-17) * T + 2.424e-15) * T + 1.3095e-12) * T + 1.7451e-10) * T - 1.8055e-8) * T - 0.0000235316) * T + 0.000076) * T + 1.105414) * T + 5028.791959) * T;
    /* Moon's longitude re fixed J2000 equinox. */
    return ARCSEC_TO_RAD * x;
}

