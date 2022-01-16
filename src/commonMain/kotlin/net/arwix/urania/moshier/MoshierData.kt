package net.arwix.urania.moshier

abstract class MoshierData {
    abstract internal val tabl: DoubleArray
    abstract internal val tabb: DoubleArray
    abstract internal val tabr: DoubleArray
    abstract internal val args: IntArray
    abstract internal val maxargs: Int
    abstract internal val max_harmonic: IntArray
    abstract internal val max_power_of_t: Int
    abstract internal val distance: Double
    internal val timescale = 3652500.0
    internal val trunclvl = 1.0
}