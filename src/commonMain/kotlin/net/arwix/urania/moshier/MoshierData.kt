package net.arwix.urania.moshier

import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.async
import kotlinx.coroutines.coroutineScope
import net.arwix.urania.core.calendar.JT
import net.arwix.urania.core.calendar.jT
import net.arwix.urania.core.calendar.toJT
import net.arwix.urania.core.ephemeris.*
import net.arwix.urania.core.math.JD_2000
import net.arwix.urania.core.math.JULIAN_DAYS_PER_CENTURY
import net.arwix.urania.core.math.LIGHT_TIME_DAYS_PER_AU
import net.arwix.urania.core.math.vector.Vector
import net.arwix.urania.core.spherical
import net.arwix.urania.core.transformation.obliquity.Obliquity
import net.arwix.urania.core.transformation.obliquity.createElements
import net.arwix.urania.core.transformation.precession.Precession
import net.arwix.urania.core.transformation.precession.PrecessionElements
import net.arwix.urania.core.transformation.precession.createElements

typealias MoshierId = Int

const val MOSHIER_ID_SUN: MoshierId = 0
const val MOSHIER_ID_MERCURY: MoshierId = 1
const val MOSHIER_ID_VENUS: MoshierId = 2
const val MOSHIER_ID_EARTH: MoshierId = 3
const val MOSHIER_ID_LIBRATION: MoshierId = 31
const val MOSHIER_ID_BARYCENTER: MoshierId = 32
const val MOSHIER_ID_MOON: MoshierId = 33
const val MOSHIER_ID_MARS: MoshierId = 4
const val MOSHIER_ID_JUPITER: MoshierId = 5
const val MOSHIER_ID_SATURN: MoshierId = 6
const val MOSHIER_ID_URANUS: MoshierId = 7
const val MOSHIER_ID_NEPTUNE: MoshierId = 8
const val MOSHIER_ID_PLUTO: MoshierId = 900


class MoshierEphemerisFactory(
    private val jT0: JT,
    private val earthEphemeris: Ephemeris
) {

    private val obliquity by lazy { Obliquity.IAU2006.createElements(JT.J2000) }

    init {
        if (earthEphemeris.metadata.plane != Plane.Ecliptic ||
            earthEphemeris.metadata.epoch != Epoch.J2000 ||
            earthEphemeris.metadata.orbit != Orbit.Heliocentric
        ) throw IllegalArgumentException()
    }

    fun createGeocentricEquatorialEphemeris(bodyEphemeris: MoshierEphemeris): Ephemeris {

//        return when (bodyEphemeris) {
//            is MainBodyEphemerisImplementation -> {
                return object : Ephemeris {
                    override val metadata: Metadata
                        get() = TODO("Not yet implemented")

                    override suspend fun invoke(jT: JT): Vector = coroutineScope {
                        val body = async(Dispatchers.Default) { bodyEphemeris(jT) }
                        val earth = async(Dispatchers.Default) { earthEphemeris(jT) }
                        val geoBody = body.await() - earth.await()
                        val oneWayDown = geoBody.spherical.r * LIGHT_TIME_DAYS_PER_AU

                        (bodyEphemeris(jT - (oneWayDown / JULIAN_DAYS_PER_CENTURY).jT) - earth.await())
                            .let { obliquity.rotatePlane(it, Plane.Equatorial) }
                    }

                }
//            }
//            else -> TODO()
//        }
    }

}


interface MoshierEphemeris: Ephemeris

object MoshierMarsEphemeris : MoshierEphemeris by MainBodyEphemerisImplementation(InnerMoshierMarsData)

object MoshierEarthMoonBarycenterEphemeris : MoshierEphemeris {
    override val metadata: Metadata = defaultMetadata

    override suspend fun invoke(jT: JT): Vector {
        return g3plan(
            jT,
            InnerMoshierEarthMoonBarycenterData.args,
            InnerMoshierEarthMoonBarycenterData.distance,
            InnerMoshierEarthMoonBarycenterData.tabb,
            InnerMoshierEarthMoonBarycenterData.tabl,
            InnerMoshierEarthMoonBarycenterData.tabr,
            InnerMoshierEarthMoonBarycenterData.max_harmonic,
            InnerMoshierEarthMoonBarycenterData.timescale,
            false
        )
    }
}

class MoshierEarthEphemeris(
    private val precessionElements: PrecessionElements,
    private val earthMoonBarycenterEphemeris: MoshierEarthMoonBarycenterEphemeris = MoshierEarthMoonBarycenterEphemeris,
    private val moonEphemeris: MoshierMoonEphemeris = MoshierMoonEphemeris
) : MoshierEphemeris {

    constructor(precessionJT0: JT,
                earthMoonBarycenterEphemeris: MoshierEarthMoonBarycenterEphemeris = MoshierEarthMoonBarycenterEphemeris,
                moonEphemeris: MoshierMoonEphemeris = MoshierMoonEphemeris
    ) : this(
        Precession.Williams1994.createElements(precessionJT0),
        earthMoonBarycenterEphemeris,
        moonEphemeris
    )

    init {
        if (precessionElements.id != Precession.Williams1994) throw IllegalArgumentException()
    }

    override val metadata: Metadata = defaultMetadata

    override suspend fun invoke(jT: JT): Vector = coroutineScope {
        val earth = async { earthMoonBarycenterEphemeris(jT) }
        val moon = async {  precessionElements.changeEpoch(moonEphemeris(jT), Epoch.J2000) }
        earth.await() - moon.await() * (1.0 / (earthMoonRatio + 1.0))
    }

    private companion object {
        private const val earthMoonRatio = 2.7068700387534E7 / 332946.050895
    }

}

object MoshierMoonEphemeris : MoshierEphemeris {
    override val metadata: Metadata
        get() = Metadata(
            orbit = Orbit.Geocentric,
            plane = Plane.Ecliptic,
            epoch = Epoch.Apparent
        )

    override suspend fun invoke(jT: JT): Vector {
        val moonLat = g1plan(
            jT,
            InnerMoshierMoonLatitudeData.args,
            InnerMoshierMoonLatitudeData.tabl,
            InnerMoshierMoonLatitudeData.max_harmonic,
            InnerMoshierMoonLatitudeData.timescale,
        )
        return g2plan(
            jT,
            InnerMoshierMoonLongitudeData.args,
            InnerMoshierMoonLongitudeData.distance,
            InnerMoshierMoonLongitudeData.tabl,
            InnerMoshierMoonLongitudeData.tabr,
            InnerMoshierMoonLongitudeData.max_harmonic,
            InnerMoshierMoonLongitudeData.timescale,
            moonLat
        )
    }
}


private class MainBodyEphemerisImplementation(private val data: InnerMoshierData) : MoshierEphemeris {
    override val metadata = defaultMetadata

    override suspend fun invoke(jT: JT): Vector {
        return gplan(jT, data.args, data.distance, data.tabb, data.tabl, data.tabr, data.max_harmonic, data.timescale)
    }
}

private val defaultMetadata = Metadata(
    orbit = Orbit.Heliocentric,
    plane = Plane.Ecliptic,
    epoch = Epoch.J2000
)