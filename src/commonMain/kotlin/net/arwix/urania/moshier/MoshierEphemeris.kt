package net.arwix.urania.moshier

import kotlinx.coroutines.async
import kotlinx.coroutines.coroutineScope
import net.arwix.urania.core.calendar.JT
import net.arwix.urania.core.ephemeris.*
import net.arwix.urania.core.math.vector.RectangularVector
import net.arwix.urania.core.math.vector.Vector
import net.arwix.urania.core.transformation.precession.Precession
import net.arwix.urania.core.transformation.precession.PrecessionElements
import net.arwix.urania.core.transformation.precession.createElements

interface MoshierEphemerisJ2000 : MoshierEphemeris {
    suspend fun fromJ2000ToApparent(jT: JT, body: Vector): Vector
    override val metadata get() = moshierEphemerisJ2000Metadata
}

private val moshierEphemerisJ2000Metadata: Metadata by lazy {
    Metadata(
        orbit = Orbit.Geocentric,
        plane = Plane.Equatorial,
        epoch = Epoch.J2000
    )
}

interface MoshierEphemerisApparent : MoshierEphemeris {
    suspend fun fromApparentToJ2000(jT: JT, body: Vector): Vector
    override val metadata get() = moshierEphemerisApparentMetadata
}

private val moshierEphemerisApparentMetadata: Metadata by lazy {
    Metadata(
        orbit = Orbit.Geocentric,
        plane = Plane.Equatorial,
        epoch = Epoch.Apparent
    )
}

interface MoshierEphemeris : Ephemeris {
    val moshierId: MoshierId
}

object MoshierSunEphemeris : MoshierEphemeris {
    override val moshierId: MoshierId = MOSHIER_ID_SUN
    override val metadata: Metadata = defaultMetadata
    override suspend fun invoke(jT: JT) = RectangularVector.Zero
}

object MoshierMercuryEphemeris :
    MoshierEphemeris by MainBodyEphemerisImplementation(MOSHIER_ID_MERCURY, InnerMoshierMercuryData)

object MoshierVenusEphemeris :
    MoshierEphemeris by MainBodyEphemerisImplementation(MOSHIER_ID_VENUS, InnerMoshierVenusData)

object MoshierMarsEphemeris : MoshierEphemeris by MainBodyEphemerisImplementation(MOSHIER_ID_MARS, InnerMoshierMarsData)
object MoshierJupiterEphemeris :
    MoshierEphemeris by MainBodyEphemerisImplementation(MOSHIER_ID_JUPITER, InnerMoshierJupiterData)

object MoshierSaturnEphemeris :
    MoshierEphemeris by MainBodyEphemerisImplementation(MOSHIER_ID_SATURN, InnerMoshierSaturnData)

object MoshierUranusEphemeris :
    MoshierEphemeris by MainBodyEphemerisImplementation(MOSHIER_ID_URANUS, InnerMoshierUranusData)

object MoshierNeptuneEphemeris :
    MoshierEphemeris by MainBodyEphemerisImplementation(MOSHIER_ID_NEPTUNE, InnerMoshierNeptuneData)

object MoshierPlutoEphemeris :
    MoshierEphemeris by MainBodyEphemerisImplementation(MOSHIER_ID_PLUTO, InnerMoshierPlutoData)

object MoshierEarthMoonBarycenterEphemeris : MoshierEphemeris {
    override val moshierId: MoshierId = MOSHIER_ID_EARTH_MOON_BARYCENTER
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

    constructor(
        precessionJT0: JT,
        earthMoonBarycenterEphemeris: MoshierEarthMoonBarycenterEphemeris = MoshierEarthMoonBarycenterEphemeris,
        moonEphemeris: MoshierMoonEphemeris = MoshierMoonEphemeris
    ) : this(
        Precession.Williams1994.createElements(precessionJT0),
        earthMoonBarycenterEphemeris,
        moonEphemeris
    )

    override val moshierId: MoshierId = MOSHIER_ID_EARTH

    init {
        if (precessionElements.id != Precession.Williams1994) throw IllegalArgumentException()
    }

    override val metadata: Metadata = defaultMetadata

    override suspend fun invoke(jT: JT): Vector = coroutineScope {
        val earth = async { earthMoonBarycenterEphemeris(jT) }
        val moon = async { precessionElements.changeEpoch(moonEphemeris(jT), Epoch.J2000) }
        earth.await() - moon.await() * (1.0 / (earthMoonRatio + 1.0))
    }

    private companion object {
        private const val earthMoonRatio = 2.7068700387534E7 / 332946.050895
    }
}

object MoshierMoonEphemeris : MoshierEphemeris {
    override val moshierId: MoshierId = MOSHIER_ID_MOON
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

private class MainBodyEphemerisImplementation(
    override val moshierId: MoshierId,
    private val data: InnerMoshierData
) : MoshierEphemeris {
    override val metadata = defaultMetadata

    override suspend fun invoke(jT: JT): Vector {
        return gplan(
            jT,
            data.args,
            data.distance,
            data.tabb,
            data.tabl,
            data.tabr,
            data.max_harmonic,
            data.timescale
        )
    }
}

private val defaultMetadata by lazy {
    Metadata(
        orbit = Orbit.Heliocentric,
        plane = Plane.Ecliptic,
        epoch = Epoch.J2000
    )
}