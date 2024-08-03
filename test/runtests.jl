using OmicsIO

using Aqua: test_all, test_ambiguities

using Test: @test

test_all(OmicsIO; ambiguities = false, deps_compat = false)

test_ambiguities(OmicsIO)

# ----------------------------------------------------------------------------------------------- #

# ---- #

di = joinpath(pkgdir(OmicsIO, "data"))

# ---- #

da = OmicsIO.read_dataframe(joinpath(di, "titanic.tsv"))

# ---- #

ts = OmicsIO.write(joinpath(tempdir(), "write.tsv"), da)

# ---- #

@test all(map(isequal, eachcol(da), eachcol(OmicsIO.read_dataframe(ts))))

# ---- #

da = OmicsIO.read_dataframe(joinpath(di, "enst_gene.tsv.gz"))

# ---- #

ts = OmicsIO.write(joinpath(tempdir(), "write.tsv"), da)

# ---- #

@test da == OmicsIO.read_dataframe(ts)

# ---- #

da = OmicsIO.read_dataframe(
    joinpath(di, "12859_2019_2886_MOESM2_ESM.xlsx"),
    "HumanSpecific Genes",
)

# ---- #

ts = OmicsIO.write(joinpath(tempdir(), "write.tsv"), da)

# ---- #

# TODO: When reading an excel sheet, infer types, 
OmicsIO.read_dataframe(ts)
