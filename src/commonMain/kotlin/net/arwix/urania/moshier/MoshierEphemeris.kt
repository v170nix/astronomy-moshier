package net.arwix.urania.moshier

import net.arwix.urania.core.calendar.MJD
import net.arwix.urania.core.calendar.mJD
import net.arwix.urania.core.calendar.toJT
import net.arwix.urania.core.ephemeris.Epoch
import net.arwix.urania.core.math.vector.RectangularVector
import net.arwix.urania.core.math.vector.Vector
import net.arwix.urania.core.transformation.precession.Precession
import net.arwix.urania.core.transformation.precession.createElements

enum class MoshierId {
    MARS, Earth
}

fun getHeliocentricEclipticPositionJ2000(mjd: MJD, id: MoshierId): Vector {

    val result = when (id) {
        MoshierId.MARS -> {
            gplan(
                mjd.toJT(),
                MoshierMarsData.args,
                MoshierMarsData.distance,
                MoshierMarsData.tabb,
                MoshierMarsData.tabl,
                MoshierMarsData.tabr,
                MoshierMarsData.max_harmonic,
                MoshierMarsData.max_power_of_t,
                MoshierMarsData.maxargs,
                MoshierMarsData.timescale,
                MoshierMarsData.trunclvl
            )
        }
        MoshierId.Earth -> {
            val p = g3plan(
                mjd.toJT(),
                MoshierEarthMoonBarycenterData.args,
                MoshierEarthMoonBarycenterData.distance,
                MoshierEarthMoonBarycenterData.tabb,
                MoshierEarthMoonBarycenterData.tabl,
                MoshierEarthMoonBarycenterData.tabr,
                MoshierEarthMoonBarycenterData.max_harmonic,
                MoshierEarthMoonBarycenterData.max_power_of_t,
                MoshierEarthMoonBarycenterData.maxargs,
                MoshierEarthMoonBarycenterData.timescale,
                MoshierEarthMoonBarycenterData.trunclvl,
                false
            )

            val moonLat = g1plan(
                mjd.toJT(),
                MoshierMoonLatitudeData.args,
                MoshierMoonLatitudeData.tabl,
                MoshierMoonLatitudeData.max_harmonic,
                MoshierMoonLatitudeData.maxargs,
                MoshierMoonLatitudeData.timescale,
            )

            val pMoon = g2plan(
                mjd.toJT(),
                MoshierMoonLongitudeData.args,
                MoshierMoonLongitudeData.distance,
                MoshierMoonLongitudeData.tabb,
                MoshierMoonLongitudeData.tabl,
                MoshierMoonLongitudeData.tabr,
                MoshierMoonLongitudeData.max_harmonic,
                MoshierMoonLongitudeData.max_power_of_t,
                MoshierMoonLongitudeData.maxargs,
                MoshierMoonLongitudeData.timescale,
                MoshierMoonLongitudeData.trunclvl,
                moonLat
            )

            val pMoonVector = Precession.Williams1994.createElements(mjd.toJT()).changeEpoch(
                RectangularVector(pMoon), Epoch.J2000
            )

            val earthMoonRatio = 2.7068700387534E7 / 332946.050895

            return RectangularVector(p) - pMoonVector * (1.0 / (earthMoonRatio + 1.0))
        }
    }

    return RectangularVector(result)

}

fun getGeocentricEclipticPositionJ2000(mjd: MJD, id: MoshierId, lightTime: Double = 0.0): Vector {
    val helioObject = getHeliocentricEclipticPositionJ2000(mjd - lightTime.mJD, id)
    val helioEarth = getHeliocentricEclipticPositionJ2000(mjd, MoshierId.Earth)

    val timeStep = 0.1
    val helioEarthPlus = getHeliocentricEclipticPositionJ2000(mjd + timeStep.mJD, MoshierId.Earth)
    val helioEarthVelocity = (helioEarthPlus - helioEarth) / timeStep
    val geoPosition = -helioEarth + helioObject
    return geoPosition
}