using OmicsTables

using Aqua: test_all, test_ambiguities

using Test: @test

test_all(OmicsTables; ambiguities = false, deps_compat = false)

test_ambiguities(OmicsTables)

# ----------------------------------------------------------------------------------------------- #
