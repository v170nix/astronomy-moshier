package net.arwix.urania.moshier

import kotlinx.coroutines.CoroutineDispatcher
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.async
import kotlinx.coroutines.coroutineScope
import net.arwix.urania.core.calendar.JT
import net.arwix.urania.core.calendar.jT
import net.arwix.urania.core.ephemeris.Ephemeris
import net.arwix.urania.core.ephemeris.Epoch
import net.arwix.urania.core.ephemeris.Orbit
import net.arwix.urania.core.ephemeris.Plane
import net.arwix.urania.core.math.LIGHT_TIME_DAYS_PER_AU
import net.arwix.urania.core.math.vector.Vector
import net.arwix.urania.core.spherical
import net.arwix.urania.core.transformation.nutation.Nutation
import net.arwix.urania.core.transformation.nutation.NutationElements
import net.arwix.urania.core.transformation.nutation.createElements
import net.arwix.urania.core.transformation.obliquity.Obliquity
import net.arwix.urania.core.transformation.obliquity.ObliquityElements
import net.arwix.urania.core.transformation.obliquity.createElements
import net.arwix.urania.core.transformation.precession.Precession
import net.arwix.urania.core.transformation.precession.PrecessionElements
import net.arwix.urania.core.transformation.precession.createElements

class MoshierEphemerisFactory(
    precessionJT: JT,
    private val dispatcher: CoroutineDispatcher = Dispatchers.Default,
    private val earthEphemeris: Ephemeris = MoshierEarthEphemeris(precessionJT),
    private val obliquity: ObliquityElements = Obliquity.Williams1994.createElements(precessionJT),
    private val precession: PrecessionElements = Precession.Williams1994.createElements(precessionJT),
    private val nutation: NutationElements = Nutation.IAU1980.createElements(precessionJT, Obliquity.Williams1994)
) {

    private val pnMatrix by lazy {
        nutation.eclipticMatrix * precession.fromJ2000Matrix
    }

    private val pnoMatrix by lazy {
        obliquity.eclipticToEquatorialMatrix * pnMatrix
    }

    private val pnoTransposeMatrix by lazy {
        precession.toJ2000Matrix * nutation.eclipticMatrix.transpose() * obliquity.equatorialToEclipticMatrix
    }

    init {
        if (precessionJT < (-50).jT || precessionJT > 10.jT) {
            throw IllegalArgumentException("Invalid date")
        }

        if (earthEphemeris.metadata.plane != Plane.Ecliptic ||
            earthEphemeris.metadata.epoch != Epoch.J2000 ||
            earthEphemeris.metadata.orbit != Orbit.Heliocentric
        ) throw IllegalArgumentException()
    }

    private suspend fun fromJ2000ToApparent(jT: JT, body: Vector): Vector {
        val earth = earthEphemeris(jT)
        val earthVelocity = earthEphemeris.getVelocity(earth, jT)
        val oneWayDown = body.spherical.r * LIGHT_TIME_DAYS_PER_AU

        return body
            .let { obliquityJ2000.rotatePlane(body, Plane.Ecliptic) }
            .let { it + earthVelocity * oneWayDown }
            .let { pnoMatrix * it }
    }

    private suspend fun fromApparentToJ2000(jT: JT, body: Vector): Vector {
        val earth = earthEphemeris(jT)
        val earthVelocity = earthEphemeris.getVelocity(earth, jT)
        val oneWayDown = body.spherical.r * LIGHT_TIME_DAYS_PER_AU

        return body
            .let { pnoTransposeMatrix * it }
            .let { it - earthVelocity * oneWayDown }
            .let { obliquityJ2000.rotatePlane(it, Plane.Equatorial) }
    }

    fun createGeocentricEquatorialEphemeris(
        bodyEphemeris: MoshierEphemeris,
        epoch: Epoch,
        plane: Plane
    ): MoshierEphemeris = when (epoch) {
        Epoch.J2000 -> when (bodyEphemeris.moshierId) {

            MOSHIER_ID_EARTH -> throw IllegalArgumentException("Observer for observer=target disallowed")

            MOSHIER_ID_MOON -> {
                object : MoshierEphemerisJ2000 {

                    override suspend fun fromJ2000ToApparent(jT: JT, body: Vector): Vector {
                        return this@MoshierEphemerisFactory.fromJ2000ToApparent(jT, body)
                    }

                    override val moshierId: MoshierId = bodyEphemeris.moshierId

                    override suspend fun invoke(jT: JT): Vector = coroutineScope {

                        val body = async(dispatcher) {
                            bodyEphemeris(jT)
                        }

                        val earthVelocity = async(dispatcher) {
                            earthEphemeris.getVelocity(earthEphemeris(jT), jT)
                        }

                        val geoBody = body.await()
                        val oneWayDown = geoBody.spherical.r * LIGHT_TIME_DAYS_PER_AU
                        val bodyVelocity = bodyEphemeris.getVelocity(geoBody, jT)
                        val result = (geoBody - bodyVelocity * oneWayDown)
                            .let { precession.changeEpoch(it, Epoch.J2000) }
                            .let { (it + (-earthVelocity.await()) * oneWayDown) }
                        if (plane == Plane.Ecliptic) return@coroutineScope result
                        obliquityJ2000.rotatePlane(result, Plane.Equatorial)
                    }
                }
            }
            else -> {
                object : MoshierEphemerisJ2000 {

                    override suspend fun fromJ2000ToApparent(jT: JT, body: Vector): Vector {
                        return this@MoshierEphemerisFactory.fromJ2000ToApparent(jT, body)
                    }

                    override val moshierId: MoshierId = bodyEphemeris.moshierId

                    override suspend fun invoke(jT: JT): Vector = coroutineScope {
                        val body = async(dispatcher) { bodyEphemeris(jT) }
                        val earth = async(dispatcher) { earthEphemeris(jT) }
                        val geoBody = body.await() - earth.await()
                        val oneWayDown = geoBody.spherical.r * LIGHT_TIME_DAYS_PER_AU
                        val bodyVelocity = bodyEphemeris.getVelocity(body.await(), jT)
                        val result = (geoBody - bodyVelocity * oneWayDown)
                        if (plane == Plane.Ecliptic) return@coroutineScope result
                        obliquityJ2000.rotatePlane(result, Plane.Equatorial)
                    }

                }
            }
        }
        Epoch.Apparent -> when (bodyEphemeris.moshierId) {

            MOSHIER_ID_EARTH -> throw IllegalArgumentException("Observer for observer=target disallowed")

            MOSHIER_ID_MOON -> {
                object : MoshierEphemerisApparent {

                    private val noMatrix by lazy {
                        obliquity.eclipticToEquatorialMatrix * nutation.eclipticMatrix
                    }

                    override suspend fun fromApparentToJ2000(jT: JT, body: Vector): Vector {
                        return this@MoshierEphemerisFactory.fromApparentToJ2000(jT, body)
                    }

                    override val moshierId: MoshierId = bodyEphemeris.moshierId

                    override suspend fun invoke(jT: JT): Vector = coroutineScope {

                        val moon = bodyEphemeris(jT)
                        val oneWayDown = moon.spherical.r * LIGHT_TIME_DAYS_PER_AU
                        val moonVelocity = bodyEphemeris.getVelocity(moon, jT)
                        val matrix = if (plane == Plane.Equatorial) noMatrix else nutation.eclipticMatrix
                        return@coroutineScope matrix * (moon - moonVelocity * oneWayDown)
                    }
                }
            }
            else -> {
                object : MoshierEphemerisApparent {

                    override suspend fun fromApparentToJ2000(jT: JT, body: Vector): Vector {
                        return this@MoshierEphemerisFactory.fromApparentToJ2000(jT, body)
                    }

                    override val moshierId: MoshierId = bodyEphemeris.moshierId

                    override suspend fun invoke(jT: JT): Vector = coroutineScope {

                        val body = async(dispatcher) { bodyEphemeris(jT) }
                        val earth = async(dispatcher) { earthEphemeris(jT) }
                        val geoBody = body.await() - earth.await()
                        val oneWayDown = geoBody.spherical.r * LIGHT_TIME_DAYS_PER_AU

                        val bodyVelocity = async(dispatcher) { bodyEphemeris.getVelocity(body.await(), jT) }
                        val earthVelocity = async(dispatcher) { earthEphemeris.getVelocity(earth.await(), jT) }

                        val oneWayBody = geoBody + (earthVelocity.await() - bodyVelocity.await()) * oneWayDown
                        val matrix = if (plane == Plane.Equatorial) pnoMatrix else pnMatrix
                        return@coroutineScope matrix * oneWayBody
                    }

                }
            }


        }
    }

    private companion object {
        private val obliquityJ2000 by lazy { Obliquity.Williams1994.createElements(JT.J2000) }
    }

}