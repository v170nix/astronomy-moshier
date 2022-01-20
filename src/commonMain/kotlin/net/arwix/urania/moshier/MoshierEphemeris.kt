package net.arwix.urania.moshier

import net.arwix.urania.core.calendar.*
import net.arwix.urania.core.ephemeris.*
import net.arwix.urania.core.math.LIGHT_TIME_DAYS_PER_AU
import net.arwix.urania.core.math.vector.*
import net.arwix.urania.core.spherical
import net.arwix.urania.core.transformation.nutation.Nutation
import net.arwix.urania.core.transformation.nutation.createElements
import net.arwix.urania.core.transformation.obliquity.Obliquity
import net.arwix.urania.core.transformation.obliquity.createElements
import net.arwix.urania.core.transformation.precession.Precession
import net.arwix.urania.core.transformation.precession.createElements

private val defaultMetadata = Metadata(
    orbit = Orbit.Heliocentric,
    plane = Plane.Ecliptic,
    epoch = Epoch.J2000
)

enum class MoshierIdBody {
    MARS, Earth, Moon
}

fun getMoonGeocentricEclipticApparentPosition(mjd: MJD): Vector {
    val moonLat = g1plan(
        mjd.toJT(),
        InnerMoshierMoonLatitudeData.args,
        InnerMoshierMoonLatitudeData.tabl,
        InnerMoshierMoonLatitudeData.max_harmonic,
        InnerMoshierMoonLatitudeData.timescale,
    )
    val pMoon = g2plan(
        mjd.toJT(),
        InnerMoshierMoonLongitudeData.args,
        InnerMoshierMoonLongitudeData.distance,
        InnerMoshierMoonLongitudeData.tabl,
        InnerMoshierMoonLongitudeData.tabr,
        InnerMoshierMoonLongitudeData.max_harmonic,
        InnerMoshierMoonLongitudeData.timescale,
        moonLat
    )
    return pMoon
}

fun getMoonGeocentricEquatorialJ2000Position(mjd: MJD): Vector {
    val geoEclipticMoonResult = getMoonGeocentricEclipticApparentPosition(mjd).spherical
    var delta = geoEclipticMoonResult.r * LIGHT_TIME_DAYS_PER_AU
    val geoEclipticMoonResult1 = getMoonGeocentricEclipticApparentPosition(mjd - delta.mJD - (38.0 / 60.0 / 60.0 / 24.0).mJD )

    // moon -mjd.toJT() true  J2000 ra 2h 11m 27.08s; lat 10deg 50m 37.1s
    // moon -mjd.toJT() false J2000 ra 2h 11m 26.82s; lat 10deg 50m 49.8s
    // moon  mjd.toJT() false J2000 ra 2h 13m 48.81s; lat 11deg 02m 54.8s
    // moon  mjd.toJT() true  J2000 ra 2h 13m 48.55s; lat 11deg 03m 07.5s

    return geoEclipticMoonResult1
        .let {
            EphemerisVector(
                it,
                metadata = Metadata(
                    orbit = Orbit.Geocentric,
                    plane = Plane.Ecliptic,
                    epoch = Epoch.Apparent
                )
            )
        }
        .let {
            Obliquity.Williams1994.createElements(mjd.toJT()).rotatePlane(it.value , Plane.Equatorial)
        }
        .let {
            Precession.Vondrak2011.createElements(mjd.toJT()).changeEpoch(it, Epoch.J2000)
        }

}

fun getMoonGeocentricEquatorialApparentPosition(mjd: MJD): Vector {
    val geoEclipticMoonResult = getMoonGeocentricEclipticApparentPosition(mjd).spherical
    val delta = geoEclipticMoonResult.r * LIGHT_TIME_DAYS_PER_AU
    val geoEclipticMoonResult1 = getMoonGeocentricEclipticApparentPosition(mjd - delta.mJD)

//    val deltaMatrix = Matrix.getRotateX((-0.1 * 0.001 * ARCSEC_TO_RAD).rad) *
//            Matrix.getRotateY((3 * 0.001 * ARCSEC_TO_RAD).rad) *
//            Matrix.getRotateZ((-5.2 * 0.001 * ARCSEC_TO_RAD).rad)

    return geoEclipticMoonResult1

        .let {
            Obliquity.Vondrak2011.createElements(mjd.toJT()).rotatePlane(it, Plane.Equatorial)
        }
        .let {
            Nutation.IAU2006.createElements(mjd.toJT(), Obliquity.Vondrak2011).apply(it, Plane.Equatorial)
        }
}

suspend fun getHeliocentricEclipticPositionJ2000(mjd: MJD, id: MoshierIdBody): Vector {
    return when (id) {
        MoshierIdBody.MARS -> {
            MoshierMarsEphemeris(mjd.toJT())
        }
        MoshierIdBody.Earth -> {
            MoshierEarthEphemeris(mjd.toJT()).invoke(mjd.toJT())
        }
        MoshierIdBody.Moon -> {
            MoshierMoonEphemeris(mjd.toJT())
        }
    }
}

suspend fun getGeocentricEclipticPositionJ2000(mjd: MJD, id: MoshierIdBody, lightTime: Double = 0.0): Vector {
    val helioObject = getHeliocentricEclipticPositionJ2000(mjd - lightTime.mJD, id)
    if (id == MoshierIdBody.Moon) return helioObject

    val helioEarth = getHeliocentricEclipticPositionJ2000(mjd, MoshierIdBody.Earth)

//    val timeStep = 0.1
//    val helioEarthPlus = getHeliocentricEclipticPositionJ2000(mjd + timeStep.mJD, MoshierId.Earth)
//    val helioEarthVelocity = (helioEarthPlus - helioEarth) / timeStep
    val geoPosition = -helioEarth + helioObject
    return geoPosition
}