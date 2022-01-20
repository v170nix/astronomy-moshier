package net.arwix.urania.moshier

abstract class InnerMoshierData {
    internal abstract val tabl: DoubleArray
    internal abstract val tabb: DoubleArray
    internal abstract val tabr: DoubleArray
    internal abstract val args: IntArray
    internal abstract val maxargs: Int
    internal abstract val max_harmonic: IntArray
    internal abstract val max_power_of_t: Int
    internal abstract val distance: Double
    internal val timescale = 3652500.0
    internal val trunclvl = 1.0
}