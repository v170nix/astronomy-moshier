package net.arwix.urania.moshier

import kotlinx.coroutines.test.runTest
import kotlinx.datetime.LocalDate
import kotlinx.datetime.TimeZone
import kotlinx.datetime.atStartOfDayIn
import net.arwix.urania.core.calendar.toJT
import net.arwix.urania.core.ephemeris.Epoch
import net.arwix.urania.core.ephemeris.Plane
import net.arwix.urania.core.math.angle.toDec
import net.arwix.urania.core.math.angle.toRA
import net.arwix.urania.core.spherical
import kotlin.test.Test
import kotlin.test.assertEquals

class MoshierEphemerisTest {

    @Test
    fun sunTest() = runTest {

        val jT = LocalDate(2022, 5, 11).atStartOfDayIn(TimeZone.UTC).toJT()
        val factory = MoshierEphemerisFactory(jT)

        val j2000Ephemeris =
            factory.createGeocentricEphemeris(MoshierSunEphemeris, Epoch.J2000, Plane.Equatorial)
        val j2000Body = j2000Ephemeris.invoke(jT).spherical

        assertEquals("3h 10m 13.88s", j2000Body.phi.toRA().toString())
        assertEquals("17deg 44m 22.0s", j2000Body.theta.toDec().toString())
        assertEquals(1.00982382016339, j2000Body.r, absoluteTolerance = 1e-7)

        val apparentEphemeris =
            factory.createGeocentricEphemeris(MoshierSunEphemeris, Epoch.Apparent, Plane.Equatorial)
        val apparentBody = apparentEphemeris(jT).spherical

        assertEquals("3h 11m 27.28s", apparentBody.phi.toRA().toString())
        assertEquals("17deg 49m 18.3s", apparentBody.theta.toDec().toString())
        assertEquals(1.00982382016339, apparentBody.r, absoluteTolerance = 1e-5)
    }

    @Test
    fun lunarTest() = runTest {

        val jT = LocalDate(1922, 6, 11).atStartOfDayIn(TimeZone.UTC).toJT()
        val factory = MoshierEphemerisFactory(jT)

        val lunarJ2000Ephemeris =
            factory.createGeocentricEphemeris(MoshierMoonEphemeris, Epoch.J2000, Plane.Equatorial)
        val lunarJ2000Body = lunarJ2000Ephemeris(jT).spherical

        assertEquals("18h 28m 41.91s", lunarJ2000Body.phi.toRA().toString())
        assertEquals("-18deg 14m 51.4s", lunarJ2000Body.theta.toDec().toString())
        assertEquals(0.002585104, lunarJ2000Body.r, absoluteTolerance = 1e-7)

        val lunarApparentEphemeris =
            factory.createGeocentricEphemeris(MoshierMoonEphemeris, Epoch.Apparent, Plane.Equatorial)
        val lunarApparentBody = lunarApparentEphemeris(jT).spherical

        assertEquals("18h 24m 10.91s", lunarApparentBody.phi.toRA().toString())
        assertEquals("-18deg 17m 40.2s", lunarApparentBody.theta.toDec().toString())
        assertEquals(0.002585104, lunarApparentBody.r, absoluteTolerance = 1e-7)
    }

    @Test
    fun mercuryTest() = runTest {

        val jT = LocalDate(2022, 5, 11).atStartOfDayIn(TimeZone.UTC).toJT()
        val factory = MoshierEphemerisFactory(jT)

        val j2000Ephemeris =
            factory.createGeocentricEphemeris(MoshierMercuryEphemeris, Epoch.J2000, Plane.Equatorial)
        val j2000Body = j2000Ephemeris.invoke(jT).spherical

        assertEquals("4h 9m 1.14s", j2000Body.phi.toRA().toString())
        assertEquals("22deg 39m 9.4s", j2000Body.theta.toDec().toString())
        assertEquals(0.6334991, j2000Body.r, absoluteTolerance = 1e-6)

        val apparentEphemeris =
            factory.createGeocentricEphemeris(MoshierMercuryEphemeris, Epoch.Apparent, Plane.Equatorial)
        val apparentBody = apparentEphemeris(jT).spherical

        assertEquals("4h 10m 18.45s", apparentBody.phi.toRA().toString())
        assertEquals("22deg 42m 35.7s", apparentBody.theta.toDec().toString())
        assertEquals(0.63349, apparentBody.r, absoluteTolerance = 1e-4)
    }

    @Test
    fun venusTest() = runTest {

        val jT = LocalDate(2022, 5, 11).atStartOfDayIn(TimeZone.UTC).toJT()
        val factory = MoshierEphemerisFactory(jT)

        val j2000Ephemeris =
            factory.createGeocentricEphemeris(MoshierVenusEphemeris, Epoch.J2000, Plane.Equatorial)
        val j2000Body = j2000Ephemeris.invoke(jT).spherical

        assertEquals("0h 36m 36.86s", j2000Body.phi.toRA().toString())
        assertEquals("2deg 5m 31.2s", j2000Body.theta.toDec().toString())
        assertEquals(1.07247543362964, j2000Body.r, absoluteTolerance = 1e-6)

        val apparentEphemeris =
            factory.createGeocentricEphemeris(MoshierVenusEphemeris, Epoch.Apparent, Plane.Equatorial)
        val apparentBody = apparentEphemeris(jT).spherical

        assertEquals("0h 37m 43.92s", apparentBody.phi.toRA().toString())
        assertEquals("2deg 12m 42.8s", apparentBody.theta.toDec().toString())
        assertEquals(1.07247543362964, apparentBody.r, absoluteTolerance = 1e-4)
    }

    @Test
    fun marsTest() = runTest {

        val jT = LocalDate(1922, 6, 11).atStartOfDayIn(TimeZone.UTC).toJT()
        val factory = MoshierEphemerisFactory(jT)

        val marsJ2000Ephemeris =
            factory.createGeocentricEphemeris(MoshierMarsEphemeris, Epoch.J2000, Plane.Equatorial)
        val marsJ2000Body = marsJ2000Ephemeris.invoke(jT).spherical

        assertEquals("17h 14m 44.87s", marsJ2000Body.phi.toRA().toString())
        assertEquals("-26deg 1m 43.2s", marsJ2000Body.theta.toDec().toString())
        assertEquals(0.46106174610382, marsJ2000Body.r, absoluteTolerance = 1e-6)

        val marsApparentTransformBody =
            (marsJ2000Ephemeris as MoshierEphemerisJ2000).fromJ2000ToApparent(jT, marsJ2000Body).spherical

        assertEquals("17h 9m 58.67s", marsApparentTransformBody.phi.toRA().toString())
        assertEquals("-25deg 56m 14.6s", marsApparentTransformBody.theta.toDec().toString())
        assertEquals(0.46106174610382, marsApparentTransformBody.r, absoluteTolerance = 1e-6)

        val marsApparentEphemeris =
            factory.createGeocentricEphemeris(MoshierMarsEphemeris, Epoch.Apparent, Plane.Equatorial)
        val marsApparentBody = marsApparentEphemeris(jT).spherical

        assertEquals("17h 9m 58.67s", marsApparentBody.phi.toRA().toString())
        assertEquals("-25deg 56m 14.6s", marsApparentBody.theta.toDec().toString())
        assertEquals(0.46106174610382, marsApparentBody.r, absoluteTolerance = 1e-6)

        val marsJ2000TransformBody =
            (marsApparentEphemeris as MoshierEphemerisApparent).fromApparentToJ2000(jT, marsApparentBody).spherical

        assertEquals("17h 14m 44.87s", marsJ2000TransformBody.phi.toRA().toString())
        assertEquals("-26deg 1m 43.2s", marsJ2000TransformBody.theta.toDec().toString())
        assertEquals(0.46106174610382, marsJ2000TransformBody.r, absoluteTolerance = 1e-6)
    }

    @Test
    fun jupiterTest() = runTest {

        val jT = LocalDate(2022, 5, 11).atStartOfDayIn(TimeZone.UTC).toJT()
        val factory = MoshierEphemerisFactory(jT)

        val j2000Ephemeris =
            factory.createGeocentricEphemeris(MoshierJupiterEphemeris, Epoch.J2000, Plane.Equatorial)
        val j2000Body = j2000Ephemeris.invoke(jT).spherical

        assertEquals("0h 0m 38.88s", j2000Body.phi.toRA().toString())
        assertEquals("-1deg 7m 35.1s", j2000Body.theta.toDec().toString())
        assertEquals(5.55341, j2000Body.r, absoluteTolerance = 1e-5)

        val apparentEphemeris =
            factory.createGeocentricEphemeris(MoshierJupiterEphemeris, Epoch.Apparent, Plane.Equatorial)
        val apparentBody = apparentEphemeris(jT).spherical

        assertEquals("0h 1m 45.92s", apparentBody.phi.toRA().toString())
        assertEquals("-1deg 0m 17.8s", apparentBody.theta.toDec().toString())
        assertEquals(5.553, apparentBody.r, absoluteTolerance = 1e-3)
    }

    @Test
    fun saturnTest() = runTest {

        val jT = LocalDate(2022, 5, 11).atStartOfDayIn(TimeZone.UTC).toJT()
        val factory = MoshierEphemerisFactory(jT)

        val j2000Ephemeris =
            factory.createGeocentricEphemeris(MoshierSaturnEphemeris, Epoch.J2000, Plane.Equatorial)
        val j2000Body = j2000Ephemeris.invoke(jT).spherical

        assertEquals("21h 48m 23.94s", j2000Body.phi.toRA().toString())
        assertEquals("-14deg 21m 51.0s", j2000Body.theta.toDec().toString())
        assertEquals(9.91839192988159, j2000Body.r, absoluteTolerance = 1e-6)

        val apparentEphemeris =
            factory.createGeocentricEphemeris(MoshierSaturnEphemeris, Epoch.Apparent, Plane.Equatorial)
        val apparentBody = apparentEphemeris(jT).spherical

        assertEquals("21h 49m 35.81s", apparentBody.phi.toRA().toString())
        assertEquals("-14deg 15m 42.3s", apparentBody.theta.toDec().toString())
        assertEquals(9.91839192988159, apparentBody.r, absoluteTolerance = 1e-3)
    }

    @Test
    fun uranusTest() = runTest {

        val jT = LocalDate(2022, 5, 11).atStartOfDayIn(TimeZone.UTC).toJT()
        val factory = MoshierEphemerisFactory(jT)

        val j2000Ephemeris =
            factory.createGeocentricEphemeris(MoshierUranusEphemeris, Epoch.J2000, Plane.Equatorial)
        val j2000Body = j2000Ephemeris.invoke(jT).spherical

        assertEquals("2h 49m 52.38s", j2000Body.phi.toRA().toString())
        assertEquals("15deg 56m 17.1s", j2000Body.theta.toDec().toString())
        assertEquals(20.7097944761007, j2000Body.r, absoluteTolerance = 1e-5)

        val apparentEphemeris =
            factory.createGeocentricEphemeris(MoshierUranusEphemeris, Epoch.Apparent, Plane.Equatorial)
        val apparentBody = apparentEphemeris(jT).spherical

        assertEquals("2h 51m 4.53s", apparentBody.phi.toRA().toString())
        assertEquals("16deg 1m 40.4s", apparentBody.theta.toDec().toString())
        assertEquals(20.7097944761007, apparentBody.r, absoluteTolerance = 1e-3)
    }

    @Test
    fun neptuneTest() = runTest {

        val jT = LocalDate(2022, 5, 11).atStartOfDayIn(TimeZone.UTC).toJT()
        val factory = MoshierEphemerisFactory(jT)

        val j2000Ephemeris =
            factory.createGeocentricEphemeris(MoshierNeptuneEphemeris, Epoch.J2000, Plane.Equatorial)
        val j2000Body = j2000Ephemeris.invoke(jT).spherical

        assertEquals("23h 41m 41.92s", j2000Body.phi.toRA().toString())
        assertEquals("-3deg 13m 35.4s", j2000Body.theta.toDec().toString())
        assertEquals(30.4782075560213, j2000Body.r, absoluteTolerance = 1e-4)

        val apparentEphemeris =
            factory.createGeocentricEphemeris(MoshierNeptuneEphemeris, Epoch.Apparent, Plane.Equatorial)
        val apparentBody = apparentEphemeris(jT).spherical

        assertEquals("23h 42m 49.19s", apparentBody.phi.toRA().toString())
        assertEquals("-3deg 6m 19.2s", apparentBody.theta.toDec().toString())
        assertEquals(30.4782075560213, apparentBody.r, absoluteTolerance = 1e-2)
    }

    @Test
    fun plutoTest() = runTest {

        val jT = LocalDate(2022, 5, 11).atStartOfDayIn(TimeZone.UTC).toJT()
        val factory = MoshierEphemerisFactory(jT)

        val j2000Ephemeris =
            factory.createGeocentricEphemeris(MoshierPlutoEphemeris, Epoch.J2000, Plane.Equatorial)
        val j2000Body = j2000Ephemeris.invoke(jT).spherical

        assertEquals("20h 3m 11.26s", j2000Body.phi.toRA().toString())
        assertEquals("-22deg 27m 23.6s", j2000Body.theta.toDec().toString())
        assertEquals(34.1337376442542, j2000Body.r, absoluteTolerance = 1e-4)

        val apparentEphemeris =
            factory.createGeocentricEphemeris(MoshierPlutoEphemeris, Epoch.Apparent, Plane.Equatorial)
        val apparentBody = apparentEphemeris(jT).spherical

        assertEquals("20h 4m 30.13s", apparentBody.phi.toRA().toString())
        assertEquals("-22deg 23m 38.7s", apparentBody.theta.toDec().toString())
        assertEquals(34.1337376442542, apparentBody.r, absoluteTolerance = 1e-2)
    }
}