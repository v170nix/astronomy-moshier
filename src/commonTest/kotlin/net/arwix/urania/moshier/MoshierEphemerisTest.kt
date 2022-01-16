package net.arwix.urania.moshier

import kotlinx.datetime.LocalDate
import kotlinx.datetime.TimeZone
import kotlinx.datetime.atStartOfDayIn
import net.arwix.urania.core.calendar.JT
import net.arwix.urania.core.calendar.toJT
import net.arwix.urania.core.calendar.toMJD
import net.arwix.urania.core.ephemeris.Epoch
import net.arwix.urania.core.ephemeris.Plane
import net.arwix.urania.core.math.angle.toDec
import net.arwix.urania.core.math.angle.toRA
import net.arwix.urania.core.spherical
import net.arwix.urania.core.toDeg
import net.arwix.urania.core.transformation.obliquity.Obliquity
import net.arwix.urania.core.transformation.obliquity.createElements
import net.arwix.urania.core.transformation.precession.Precession
import net.arwix.urania.core.transformation.precession.createElements
import kotlin.test.Test

class MoshierEphemerisTest {

    @Test
    fun getHeliocentricEclipticPositionJ2000Test() {
        val mjd = LocalDate(2022, 1, 11).atStartOfDayIn(TimeZone.UTC).toMJD()

        val result = getHeliocentricEclipticPositionJ2000(mjd, MoshierId.MARS).spherical

        println("lon ${result.phi.toDeg()}")
        println("lat ${result.theta.toDeg()}")
        println("distance ${result.r}")

        val geoEclipticResult = getGeocentricEclipticPositionJ2000(mjd, MoshierId.MARS, 19.02186655 / 60.0 / 24.0).spherical
        val geoResult = geoEclipticResult
            .let {
                Obliquity.Williams1994.createElements(JT.J2000).rotatePlane(it, Plane.Equatorial)
            }
            .spherical

        println("ra ${geoResult.phi.toRA()}")
        println("lat ${geoResult.theta.toDec()}")
        println("distance ${geoResult.r}")

    }
}