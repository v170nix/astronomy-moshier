package net.arwix.urania.moshier

import kotlinx.coroutines.test.runTest
import kotlinx.datetime.LocalDate
import kotlinx.datetime.TimeZone
import kotlinx.datetime.atStartOfDayIn
import net.arwix.urania.core.calendar.JT
import net.arwix.urania.core.calendar.toJT
import net.arwix.urania.core.calendar.toMJD
import net.arwix.urania.core.ephemeris.*
import net.arwix.urania.core.math.RAD_TO_ARCSEC
import net.arwix.urania.core.math.angle.*
import net.arwix.urania.core.math.vector.SphericalVector
import net.arwix.urania.core.rectangular
import net.arwix.urania.core.spherical
import net.arwix.urania.core.toDeg
import net.arwix.urania.core.toRad
import net.arwix.urania.core.transformation.nutation.Nutation
import net.arwix.urania.core.transformation.nutation.createElements
import net.arwix.urania.core.transformation.obliquity.Obliquity
import net.arwix.urania.core.transformation.obliquity.createElements
import net.arwix.urania.core.transformation.precession.Precession
import net.arwix.urania.core.transformation.precession.createElements
import kotlin.test.Test
import kotlin.test.assertEquals

class MoshierEphemerisTest {

    @Test
    fun getHeliocentricEclipticPositionJ2000Test() = runTest {

        val mjd = LocalDate(1922, 6, 11).atStartOfDayIn(TimeZone.UTC).toMJD()

        val result = getHeliocentricEclipticPositionJ2000(mjd, MoshierIdBody.MARS).spherical

        println("lon ${result.phi.toDeg()}")
        println("lat ${result.theta.toDeg()}")
        println("distance ${result.r}")

        val marsEphemeris = MoshierEphemerisFactory(
            mjd.toJT(),
            earthEphemeris = MoshierEarthEphemeris(mjd.toJT())
        ).createGeocentricEquatorialEphemeris(MoshierMarsEphemeris)

        val geoEclipticResult = getGeocentricEclipticPositionJ2000(mjd, MoshierIdBody.MARS, 3.83453362  / 60.0 / 24.0).spherical
        var geoResult = geoEclipticResult
            .let {
                Obliquity.Williams1994.createElements(JT.J2000).rotatePlane(it, Plane.Equatorial)
            }
            .spherical

        geoResult = marsEphemeris(mjd.toJT()).spherical

        println("mars ra ${geoResult.phi.toRA()}")
        println("mars lat ${geoResult.theta.toDec()}")
        println("mars distance ${geoResult.r}")

        assertEquals("17h 14m 44.87s", geoResult.phi.toRA().toString())
        assertEquals("-26deg 1m 43.2s", geoResult.theta.toDec().toString())
        assertEquals(0.46106149461094065, geoResult.r, absoluteTolerance = 1e-10)

        val geoMoonResult = getMoonGeocentricEquatorialApparentPosition(mjd).spherical

        println("moon apparent ra ${geoMoonResult.phi.toRA()}")
        println("moon apparent lat ${geoMoonResult.theta.toDec()}")
        println("moon apparent distance ${geoMoonResult.r}")

        assertEquals("18h 24m 10.91s", geoMoonResult.phi.toRA().toString())
        assertEquals("-18deg 17m 40.2s", geoMoonResult.theta.toDec().toString())
        assertEquals(0.0025851170027739352, geoMoonResult.r)

        val geoMoonJ2000Result = getMoonGeocentricEquatorialJ2000Position(mjd).spherical

        println("moon J2000 ra ${geoMoonJ2000Result.phi.toRA()}")
        println("moon J2000 lat ${geoMoonJ2000Result.theta.toDec()}")
        println("moon J2000 distance ${geoMoonJ2000Result.r}")

        assertEquals("18h 28m 41.92s", geoMoonJ2000Result.phi.toRA().toString())
        assertEquals("-18deg 14m 52.1s", geoMoonJ2000Result.theta.toDec().toString())
        assertEquals(0.0025851042460009704, geoMoonJ2000Result.r)

    }
}