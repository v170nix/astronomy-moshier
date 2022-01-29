package net.arwix.urania.moshier

import net.arwix.urania.core.calendar.JT
import net.arwix.urania.core.calendar.times
import net.arwix.urania.core.kepler.KeplerElementsObject
import net.arwix.urania.core.kepler.getKeplerElementsSimonJ2000
import net.arwix.urania.core.math.ARCSEC_TO_RAD
import net.arwix.urania.core.math.JULIAN_DAYS_PER_CENTURY
import net.arwix.urania.core.math.angle.Radian
import net.arwix.urania.core.math.angle.rad
import net.arwix.urania.core.math.angle.times
import net.arwix.urania.core.math.mod3600
import net.arwix.urania.core.math.vector.SphericalVector
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
private inline fun sscc(k: Int, arg: Double, n: Int, ss: Array<DoubleArray>, cc: Array<DoubleArray>) {
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

    for (i in 2 until n) {
        s = su * cv + cu * sv
        cv = cu * cv - su * sv
        sv = s
        ss[k][i] = sv /* sin( i+1 L ) */
        cc[k][i] = cv
    }
}

/**
 * Generic program to accumulate sum of trigonometric series in three
 * variables (e.g., longitude, latitude, radius) of the same list of
 * arguments.
 *
 * @param tt
 * @param arg_tbl
 * @param distance
 * @param lat_tbl
 * @param lon_tbl
 * @param rad_tbl
 * @param max_harmonic
 * @param timescale
 * @return SphericalVector r (AU).
 */
internal fun gplan(
    tt: JT,
    arg_tbl: IntArray,
    distance: Double,
    lat_tbl: DoubleArray,
    lon_tbl: DoubleArray,
    rad_tbl: DoubleArray,
    max_harmonic: IntArray,
    timescale: Double
): SphericalVector {

    var j: Int
    var k: Int
    var m: Int
    var np: Int
    var su: Double
    var cu: Double
    var buffer: Double

    if (distance == 0.0) {
        throw IllegalArgumentException()
    }

    val jT100: Double = (tt * JULIAN_DAYS_PER_CENTURY) / timescale

    /* Calculate sin( i*MM ), etc. for needed multiple angles. */
    val ss = Array(20) { DoubleArray(41) }
    val cc = Array(20) { DoubleArray(41) }

    for (i in 0 until 9) {
        if (max_harmonic[i] > 0) {
            buffer = ((freqs[i] * jT100).mod3600() + phases[i]) * ARCSEC_TO_RAD
            sscc(i, buffer, max_harmonic[i], ss, cc)
        }
    }

    var pIndex = -1
    var plIndex = -1
    var pbIndex = -1
    var prIndex = -1
    var count: Int
    var accumulator: Double
    var isFirst: Boolean
    val arrayAccumulator = doubleArrayOf(0.0, 0.0)

    val vector = SphericalVector.Zero

    while (true) {

        np = arg_tbl[++pIndex]
        if (np < 0) break
        if (np == 0) {

            /* It is a polynomial term. */
            count = arg_tbl[++pIndex]

            /* "Longitude" polynomial (phi). */
            accumulator = lon_tbl[++plIndex]
            for (i in (0 until count)) {
                accumulator = accumulator * jT100 + lon_tbl[++plIndex]
            }
            vector.phi += accumulator.mod3600().rad


            /* "Latitude" polynomial (theta). */
            accumulator = lat_tbl[++pbIndex]
            for (i in (0 until count)) {
                accumulator = accumulator * jT100 + lat_tbl[++pbIndex]
            }
            vector.theta += accumulator.rad

            /* Radius polynomial (psi). */
            accumulator = rad_tbl[++prIndex]
            for (i in (0 until count)) {
                accumulator = accumulator * jT100 + rad_tbl[++prIndex]
            }
            vector.r += accumulator
            continue
        }


        isFirst = true

        for (i in (0 until np)) {
            /* What harmonic. */
            j = arg_tbl[++pIndex]
            /* Which planet. */
            m = arg_tbl[++pIndex] - 1
            if (j != 0) {
                k = abs(j) - 1
                su = ss[m][k] /* sin(k*angle) */
                if (j < 0) su = -su
                cu = cc[m][k]
                if (isFirst) { /* set first angle */
                    arrayAccumulator[0] = su
                    arrayAccumulator[1] = cu
                    isFirst = false
                } else { /* combine angles */
                    buffer = su * arrayAccumulator[1] + cu * arrayAccumulator[0]
                    arrayAccumulator[1] = cu * arrayAccumulator[1] - su * arrayAccumulator[0]
                    arrayAccumulator[0] = buffer
                }
            }
        }
        val (sv, cv) = arrayAccumulator

        /* Highest power of T. */
        count = arg_tbl[++pIndex]

        /* Longitude. */

        arrayAccumulator[0] = lon_tbl[++plIndex]
        arrayAccumulator[1] = lon_tbl[++plIndex]
        for (i in (0 until count)) {
            arrayAccumulator[0] = arrayAccumulator[0] * jT100 + lon_tbl[++plIndex]
            arrayAccumulator[1] = arrayAccumulator[1] * jT100 + lon_tbl[++plIndex]
        }
        vector.phi += arrayAccumulator.let { it[0] * cv + it[1] * sv }.rad

        /* Latitude. */
        arrayAccumulator[0] = lat_tbl[++pbIndex]
        arrayAccumulator[1] = lat_tbl[++pbIndex]
        for (i in (0 until count)) {
            arrayAccumulator[0] = arrayAccumulator[0] * jT100 + lat_tbl[++pbIndex]
            arrayAccumulator[1] = arrayAccumulator[1] * jT100 + lat_tbl[++pbIndex]
        }
        vector.theta += arrayAccumulator.let { it[0] * cv + it[1] * sv }.rad

        /* Radius. */
        arrayAccumulator[0] = rad_tbl[++prIndex]
        arrayAccumulator[1] = rad_tbl[++prIndex]
        for (i in (0 until count)) {
            arrayAccumulator[0] = arrayAccumulator[0] * jT100 + rad_tbl[++prIndex]
            arrayAccumulator[1] = arrayAccumulator[1] * jT100 + rad_tbl[++prIndex]
        }
        vector.r += arrayAccumulator.let { it[0] * cv + it[1] * sv }
    }

//    if (distance == 0.0) return doubleArrayOf(
//        (ARCSEC_TO_RAD * vector.phi).rad.normalize().value,
//        (ARCSEC_TO_RAD * vector.theta).rad.normalize().value,
//        (ARCSEC_TO_RAD * vector.r).rad.normalize().value
//    )

    vector.phi *= ARCSEC_TO_RAD
    vector.theta *= ARCSEC_TO_RAD
    vector.r = distance * (1.0 + vector.r * ARCSEC_TO_RAD)

    return vector
}

/* From Simon et al (1994) */
private val freqs  by lazy {
    doubleArrayOf( /* Arc sec per 10000 Julian years. */
        53810162868.8982, 21066413643.3548, 12959774228.3429, 6890507749.3988, 1092566037.7991, 439960985.5372,
        154248119.3933, 78655032.0744, 52272245.1795
    )
}
private val phases by lazy {
    doubleArrayOf( /* Arc sec. */
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
}

/**
 * Generic program to accumulate sum of trigonometric series in one variable (e.g., latitude) of the same list of arguments.
 *
 * @param tt
 * @param arg_tbl
 * @param lat_tbl
 * @param max_harmonic
 * @param timescale
 * @return Latitude (rad)
 */
internal fun g1plan(
    tt: JT, arg_tbl: IntArray, lat_tbl: DoubleArray, max_harmonic: IntArray, timescale: Double,
): Radian {
    var j: Int
    var k: Int
    var m: Int
    var buffer: Double

    var np: Int

    var su: Double
    var cu: Double

    val args: DoubleArray = getMeanElements(tt)

    val jT100: Double = (tt * JULIAN_DAYS_PER_CENTURY) / timescale

    /* Calculate sin( i*MM ), etc. for needed multiple angles. */
    val ss = Array(20) { DoubleArray(41) }
    val cc = Array(20) { DoubleArray(41) }

    max_harmonic.forEachIndexed { i, harmonic -> if (harmonic > 0) sscc(i, args[i], max_harmonic[i], ss, cc) }

    var sb = 0.0
    var pIndex = -1
    var pbIndex = -1

    var count: Int

    var isFirst: Boolean

    var accumulator: Double
    val arrayAccumulator = doubleArrayOf(0.0, 0.0)

    while (true) {

        /* argument of sine and cosine */
        /* Number of periodic arguments. */
        np = arg_tbl[++pIndex]

        if (np < 0) break
        if (np == 0) {
            /* It is a polynomial term. */
            count = arg_tbl[++pIndex]
            accumulator = lat_tbl[++pbIndex]
            for (i in (0 until  count)) {
                accumulator = accumulator * jT100 + lat_tbl[++pbIndex]
            }
            sb += accumulator
            continue
        }

        isFirst = true

        for (i in (0 until np)) {
            /* What harmonic. */
            j = arg_tbl[++pIndex]
            /* Which planet. */
            m = arg_tbl[++pIndex] - 1
            if (j != 0) {
                k = abs(j) - 1
                su = ss[m][k] /* sin(k*angle) */
                if (j < 0) su = -su
                cu = cc[m][k]
                if (isFirst) { /* set first angle */
                    arrayAccumulator[0] = su
                    arrayAccumulator[1] = cu
                    isFirst = false
                } else { /* combine angles */
                    buffer = su * arrayAccumulator[1] + cu * arrayAccumulator[0]
                    arrayAccumulator[1] = cu * arrayAccumulator[1] - su * arrayAccumulator[0]
                    arrayAccumulator[0] = buffer
                }
            }
        }
        val (sv, cv) = arrayAccumulator

        count = arg_tbl[++pIndex]

        /* Latitude. */
        arrayAccumulator[0] = lat_tbl[++pbIndex]
        arrayAccumulator[1] = lat_tbl[++pbIndex]
        for (i in (0 until count)) {
            arrayAccumulator[0] = arrayAccumulator[0] * jT100 + lat_tbl[++pbIndex]
            arrayAccumulator[1] = arrayAccumulator[1] * jT100 + lat_tbl[++pbIndex]
        }
        sb += arrayAccumulator.let { it[0] * cv + it[1] * sv }

    }
    return (ARCSEC_TO_RAD * sb * 0.0001).rad
}

/**
 * Generic program to accumulate sum of trigonometric series in two
 * variables (e.g., longitude, radius) of the same list of arguments.
 *
 * @param tt
 * @param arg_tbl
 * @param distance
 * @param lon_tbl
 * @param rad_tbl
 * @param max_harmonic
 * @param timescale
 * @return vector with phi, theta, r (AU).
 */
internal fun g2plan(
    tt: JT,
    arg_tbl: IntArray,
    distance: Double,
    lon_tbl: DoubleArray,
    rad_tbl: DoubleArray,
    max_harmonic: IntArray,
    timescale: Double,
    lat: Radian
): SphericalVector {

    var j: Int
    var k: Int
    var m: Int
    var np: Int
    var su: Double
    var cu: Double
    var buffer: Double

    // TODO moon libration
    if (distance == 0.0) throw IllegalArgumentException()

    val jT100: Double = (tt * JULIAN_DAYS_PER_CENTURY) / timescale
    val args: DoubleArray = getMeanElements(tt)
    val lpEquinox = args[13]

    /* Calculate sin( i*MM ), etc. for needed multiple angles. */
    val ss = Array(20) { DoubleArray(41) }
    val cc = Array(20) { DoubleArray(41) }

    max_harmonic.forEachIndexed { i, harmonic ->
        if (harmonic > 0) sscc(i, args[i], max_harmonic[i], ss, cc)
    }

    val vector = SphericalVector.Zero

    var pIndex = -1
    var plIndex = -1
    var prIndex = -1

    var accumulator: Double
    val arrayAccumulator = doubleArrayOf(0.0, 0.0)
    var count: Int
    var isFirst: Boolean

    while (true) {
        /* argument of sine and cosine */
        /* Number of periodic arguments. */
        np = arg_tbl[++pIndex]
        if (np < 0) break
        if (np == 0) {
            /* It is a polynomial term. */
            count = arg_tbl[++pIndex]

            /* "Longitude" polynomial (phi). */
            accumulator = lon_tbl[++plIndex]
            for (i in (0 until count)) {
                accumulator = accumulator * jT100 + lon_tbl[++plIndex]
            }
            vector.phi += accumulator.rad

            /* Radius polynomial (psi). */
            accumulator = rad_tbl[++prIndex]
            for (i in (0 until count)) {
                accumulator = accumulator * jT100 + rad_tbl[++prIndex]
            }
            vector.r += accumulator
            continue
        }

        isFirst = true

        for (i in (0 until np)) {
            /* What harmonic. */
            j = arg_tbl[++pIndex]
            /* Which planet. */
            m = arg_tbl[++pIndex] - 1
            if (j != 0) {
                k = abs(j) - 1
                su = ss[m][k] /* sin(k*angle) */
                if (j < 0) su = -su
                cu = cc[m][k]
                if (isFirst) { /* set first angle */
                    arrayAccumulator[0] = su
                    arrayAccumulator[1] = cu
                    isFirst = false
                } else { /* combine angles */
                    buffer = su * arrayAccumulator[1] + cu * arrayAccumulator[0]
                    arrayAccumulator[1] = cu * arrayAccumulator[1] - su * arrayAccumulator[0]
                    arrayAccumulator[0] = buffer
                }
            }
        }

        val (sv, cv) = arrayAccumulator

        /* Highest power of T. */
        count = arg_tbl[++pIndex]

        /* Longitude. */
        arrayAccumulator[0] = lon_tbl[++plIndex]
        arrayAccumulator[1] = lon_tbl[++plIndex]
        for (i in (0 until count)) {
            arrayAccumulator[0] = arrayAccumulator[0] * jT100 + lon_tbl[++plIndex]
            arrayAccumulator[1] = arrayAccumulator[1] * jT100 + lon_tbl[++plIndex]
        }
        vector.phi += arrayAccumulator.let { it[0] * cv + it[1] * sv }.rad

        /* Radius. */
        arrayAccumulator[0] = rad_tbl[++prIndex]
        arrayAccumulator[1] = rad_tbl[++prIndex]
        for (i in (0 until count)) {
            arrayAccumulator[0] = arrayAccumulator[0] * jT100 + rad_tbl[++prIndex]
            arrayAccumulator[1] = arrayAccumulator[1] * jT100 + rad_tbl[++prIndex]
        }
        vector.r += arrayAccumulator.let { it[0] * cv + it[1] * sv }
    }
    vector.phi *= 0.0001
    vector.r *= 0.0001
    // TODO libration angles
//    if (distance == 0.0) return doubleArrayOf(
//        ARCSEC_TO_RAD * sl + lpEquinox,
//        lat, ARCSEC_TO_RAD * sr
//    )

    vector.phi = (ARCSEC_TO_RAD * vector.phi + lpEquinox).rad
    vector.theta = lat
    vector.r = distance * (1.0 + ARCSEC_TO_RAD * vector.r)
    return vector
}

/**
 * Generic program to accumulate sum of trigonometric series in three
 * variables (e.g., longitude, latitude, radius) of the same list of
 * arguments.
 *
 * @param tt
 * @param arg_tbl
 * @param distance
 * @param lat_tbl
 * @param lon_tbl
 * @param rad_tbl
 * @param max_harmonic
 * @param timescale
 * @return An SphericalVector r (AU)
 */
internal fun g3plan(
    tt: JT, arg_tbl: IntArray, distance: Double, lat_tbl: DoubleArray, lon_tbl: DoubleArray,
    rad_tbl: DoubleArray, max_harmonic: IntArray, timescale: Double, libration: Boolean
): SphericalVector {
    var j: Int
    var k: Int
    var m: Int
    var isFirst: Boolean
    var np: Int
    var su: Double
    var cu: Double
    var buffer: Double

    if (distance == 0.0 || libration) {
        //TODO libration
        throw IllegalArgumentException()
    }

    val jT = (tt * JULIAN_DAYS_PER_CENTURY) / timescale

    val args: DoubleArray = getMeanElements(tt)

    if (libration) args[13] -= getPrecessionOfEquinox(tt) // Only librations

    /* Calculate sin( i*MM ), etc. for needed multiple angles. */
    val ss = Array(20) { DoubleArray(41) }
    val cc = Array(20) { DoubleArray(41) }

    max_harmonic.forEachIndexed { i, harmonic -> if (harmonic > 0) sscc(i, args[i], max_harmonic[i], ss, cc) }

    var pIndex = -1
    var plIndex = -1
    var pbIndex = -1
    var prIndex = -1

    var accumulator: Double
    val arrayAccumulator = doubleArrayOf(0.0, 0.0)
    var count: Int

    val vector = SphericalVector.Zero
    while (true) {

        /* argument of sine and cosine */
        /* Number of periodic arguments. */
        np = arg_tbl[++pIndex]
        if (np < 0) break
        if (np == 0) {

            /* It is a polynomial term. */
            count = arg_tbl[++pIndex]

            /* "Longitude" polynomial (phi). */
            accumulator = lon_tbl[++plIndex]
            for (i in (0 until count)) {
                accumulator = accumulator * jT + lon_tbl[++plIndex]
            }
            vector.phi += accumulator.rad


            /* "Latitude" polynomial (theta). */
            accumulator = lat_tbl[++pbIndex]
            for (i in (0 until count)) {
                accumulator = accumulator * jT + lat_tbl[++pbIndex]
            }
            vector.theta += accumulator.rad

            /* Radius polynomial (psi). */
            accumulator = rad_tbl[++prIndex]
            for (i in (0 until count)) {
                accumulator = accumulator * jT + rad_tbl[++prIndex]
            }
            vector.r += accumulator
            continue
        }

        isFirst = true

        for (i in (0 until np)) {
            /* What harmonic. */
            j = arg_tbl[++pIndex]
            /* Which planet. */
            m = arg_tbl[++pIndex] - 1
            if (j != 0) {
                k = abs(j) - 1
                su = ss[m][k] /* sin(k*angle) */
                if (j < 0) su = -su
                cu = cc[m][k]
                if (isFirst) { /* set first angle */
                    arrayAccumulator[0] = su
                    arrayAccumulator[1] = cu
                    isFirst = false
                } else { /* combine angles */
                    buffer = su * arrayAccumulator[1] + cu * arrayAccumulator[0]
                    arrayAccumulator[1] = cu * arrayAccumulator[1] - su * arrayAccumulator[0]
                    arrayAccumulator[0] = buffer
                }
            }
        }
        val (sv, cv) = arrayAccumulator

        /* Highest power of T. */
        count = arg_tbl[++pIndex]

        /* Longitude. */

        arrayAccumulator[0] = lon_tbl[++plIndex]
        arrayAccumulator[1] = lon_tbl[++plIndex]
        for (i in (0 until count)) {
            arrayAccumulator[0] = arrayAccumulator[0] * jT + lon_tbl[++plIndex]
            arrayAccumulator[1] = arrayAccumulator[1] * jT + lon_tbl[++plIndex]
        }
        vector.phi += arrayAccumulator.let { it[0] * cv + it[1] * sv }.rad

        /* Latitude. */
        arrayAccumulator[0] = lat_tbl[++pbIndex]
        arrayAccumulator[1] = lat_tbl[++pbIndex]
        for (i in (0 until count)) {
            arrayAccumulator[0] = arrayAccumulator[0] * jT + lat_tbl[++pbIndex]
            arrayAccumulator[1] = arrayAccumulator[1] * jT + lat_tbl[++pbIndex]
        }
        vector.theta += arrayAccumulator.let { it[0] * cv + it[1] * sv }.rad

        /* Radius. */
        arrayAccumulator[0] = rad_tbl[++prIndex]
        arrayAccumulator[1] = rad_tbl[++prIndex]
        for (i in (0 until count)) {
            arrayAccumulator[0] = arrayAccumulator[0] * jT + rad_tbl[++prIndex]
            arrayAccumulator[1] = arrayAccumulator[1] * jT + rad_tbl[++prIndex]
        }
        vector.r += arrayAccumulator.let { it[0] * cv + it[1] * sv }

    }
    vector.phi *= 0.0001
    vector.theta *= 0.0001
    vector.r *= 0.0001

    // TODO Libration angles
//    if (distance == 0.0) return
//    doubleArrayOf(
//        ARCSEC_TO_RAD * sl + args[2], // + Ea_arcsec,
//        ARCSEC_TO_RAD * sb,
//        ARCSEC_TO_RAD * sr
//    )

    vector.phi = vector.phi * ARCSEC_TO_RAD + args[2].rad // Ea_rad
    vector.theta = vector.theta * ARCSEC_TO_RAD
    vector.r = distance * (1.0 + vector.r * ARCSEC_TO_RAD)

    return vector
}

/**
 * Obtain mean elements of the planets.
 * @param jT Julian centuries
 * @return An array with the mean longitudes.
 */
private fun getMeanElements(jT: JT): DoubleArray {
    val args = DoubleArray(20)

    /** Mean longitudes of planets (Simon 1994) .047" subtracted from
     * constant term for offset to DE403 origin.
     */

    val delta = (-.047 * ARCSEC_TO_RAD).rad

    args[0] = (getKeplerElementsSimonJ2000(KeplerElementsObject.Mercury).getLongitude(jT) + delta).value
    args[1] = (getKeplerElementsSimonJ2000(KeplerElementsObject.Venus).getLongitude(jT) + delta).value
    args[2] = (getKeplerElementsSimonJ2000(KeplerElementsObject.Earth).getLongitude(jT) + delta).value
    args[3] = (getKeplerElementsSimonJ2000(KeplerElementsObject.Mars).getLongitude(jT) + delta).value
    args[4] = (getKeplerElementsSimonJ2000(KeplerElementsObject.Jupiter).getLongitude(jT) + delta).value
    args[5] = (getKeplerElementsSimonJ2000(KeplerElementsObject.Saturn).getLongitude(jT) + delta).value
    args[6] = (getKeplerElementsSimonJ2000(KeplerElementsObject.Uranus).getLongitude(jT) + delta).value
    args[7] = (getKeplerElementsSimonJ2000(KeplerElementsObject.Neptune).getLongitude(jT) + delta).value

    /* Copied from cmoon.c, DE404 version. */
    /* Mean elongation of moon = elongation */
    args[9] = getMeanElongationOfMoon(jT)
    args[10] = getAscendingNode(jT)
    args[11] = getMeanAnomalyOfSun(jT)
    args[12] = getMeanAnomalyOfMoon(jT)
    args[13] = getLongitudeOfMoon(jT)
    args[14] = getLunarFreeLibrations(jT)
    args[15] = getLB(jT)
    args[16] = getLC(jT)
    args[17] = args[13] - args[10]
    args[18] = getNB(args[17], jT)
    return args
}

private fun getMeanElongationOfMoon(jT: JT): Double {
    val jT2 = jT * jT
    /* Mean elongation of moon = D */
    var x = (1.6029616009939659e+09 * jT + 1.0722612202445078e+06)
    x += (((((-3.207663637426e-013 * jT + 2.555243317839e-011) * jT + 2.560078201452e-009) * jT - 3.702060118571e-005) * jT + 0.006949274683605842) * jT /* D, t^3 */
            - 6.735220237445752) * jT2 /* D, t^2 */

    return ARCSEC_TO_RAD * x
}

internal fun getAscendingNode(jT: JT): Double {
    val jT2 = jT * jT
    /* Mean distance of moon from its ascending node = F */
    var x = (1.7395272628437717e+09 * jT + 3.3577951412884740e+05)
    x += (((((4.474984866301e-013 * jT + 4.189032191814e-011) * jT - 2.790392351314e-009) * jT - 2.165750777942e-006) * jT - 7.531187848233799E-4) * jT /* F, t^3 */
            - 1.3117809789650071e+01) * jT2 /* F, t^2 */
    return ARCSEC_TO_RAD * x
}

private fun getMeanAnomalyOfSun(jT: JT): Double {
    /* Mean anomaly of sun = l' (J. Laskar) */
    val jT2 = jT * jT
    var x = (1.2959658102304320e+08 * jT + 1.2871027407441526e+06)
    x += ((((((((1.62e-20 * jT - 1.0390e-17) * jT - 3.83508e-15) * jT + 4.237343e-13) * jT + 8.8555011e-11) * jT - 4.77258489e-8) * jT - 1.1297037031e-5) * jT + 8.74737173673247E-5) * jT - 0.5528130642178309) * jT2
    return ARCSEC_TO_RAD * x
}

private fun getMeanAnomalyOfMoon(jT: JT): Double {
    val jT2 = jT * jT
    var x = (1.7179159228846793e+09 * jT + 485868.1746582533)
    x += (((((-1.755312760154e-012 * jT + 3.452144225877e-011) * jT - 2.506365935364e-008) * jT - 2.536291235258e-004) * jT + 0.05209964130273582) * jT /* l, t^3 */
            + 3.1501359071894147e+01) * jT2 /* l, t^2 */
    return ARCSEC_TO_RAD * x
}

internal fun getLongitudeOfMoon(jT: JT): Double {
    val jT2 = jT * jT
    /* Mean longitude of moon, re mean ecliptic and equinox of date = L */
    var x = (1.7325643720442266e+09 * jT + 7.8593980921052420e+05)
    x += (((((7.200592540556e-014 * jT + 2.235210987108e-010) * jT - 1.024222633731e-008) * jT - 6.073960534117e-005) * jT + 6.9017248528380490e-03) * jT /* L, t^3 */
            - 5.65504600274714) * jT2 /* L, t^2 */
    return ARCSEC_TO_RAD * x
}

private fun getLunarFreeLibrations(jT: JT): Double {
    /* Lunar free librations. */
    /* 74.7 years. Denoted W or LA. */
    val x = (-0.112 * jT + 1.73655499e6) * jT - 389552.81
    return ARCSEC_TO_RAD * x
}

private fun getLB(jT: JT): Double {
    return ARCSEC_TO_RAD * (4.48175409e7 * jT + 806045.7)
}

private fun getLC(jT: JT): Double {
    return ARCSEC_TO_RAD * (5.36486787e6 * jT - 391702.8)
}

private fun getNB(lc: Double, jT: JT): Double {
    val x = (((-0.000004 * jT + 0.000026) * jT + 0.153382) * jT - 867.919986) * jT + 629543.967373
    return lc + ARCSEC_TO_RAD * (3.24e5 - x) - getPrecessionOfEquinox(jT)
}

private fun getPrecessionOfEquinox(jT: JT): Double {
    val pAPrecession =
        (((((((((-8.66e-20 * jT - 4.759e-17) * jT + 2.424e-15) * jT + 1.3095e-12) * jT + 1.7451e-10) * jT - 1.8055e-8) * jT - 0.0000235316) * jT + 0.000076) * jT + 1.105414) * jT + 5028.791959) * jT
    /* Moon's longitude re fixed J2000 equinox. */
    return ARCSEC_TO_RAD * pAPrecession
}

