using LeMoIO

using Aqua: test_all, test_ambiguities

using Test: @test

test_all(LeMoIO; ambiguities = false, deps_compat = false)

test_ambiguities(LeMoIO)

# ----------------------------------------------------------------------------------------------- #

# ---- #

di = joinpath(pkgdir(LeMoIO, "data"))

# ---- #

da = LeMoIO.read_dataframe(joinpath(di, "titanic.tsv"))

# ---- #

ts = LeMoIO.write(joinpath(tempdir(), "write.tsv"), da)

# ---- #

@test all(map(isequal, eachcol(da), eachcol(LeMoIO.read_dataframe(ts))))

# ---- #

da = LeMoIO.read_dataframe(joinpath(di, "enst_gene.tsv.gz"))

# ---- #

ts = LeMoIO.write(joinpath(tempdir(), "write.tsv"), da)

# ---- #

@test da == LeMoIO.read_dataframe(ts)

# ---- #

da = LeMoIO.read_dataframe(
    joinpath(di, "12859_2019_2886_MOESM2_ESM.xlsx"),
    "HumanSpecific Genes",
)

# ---- #

ts = LeMoIO.write(joinpath(tempdir(), "write.tsv"), da)

# ---- #

# TODO: When reading an excel sheet, infer types, 
LeMoIO.read_dataframe(ts)

# ---- #
const DA = joinpath(Nucleus._DA, "CLS")
@test Nucleus.Path.read(DA) ==
      ["CCLE_mRNA_20Q2_no_haem_phen.cls", "GSE76137.cls", "LPS_phen.cls"]
for n in (10, 100, 1000)
    st_ = string.(1:n)
    # 172.545 ns (1 allocation: 144 bytes)
    # 1.346 μs (1 allocation: 896 bytes)
    # 13.334 μs (1 allocation: 7.94 KiB)
    #@btime Nucleus.CLS._make_matrix(st -> parse(Float64, st), $st_)
end
for (cl, ta, re) in (
    (
        "CCLE_mRNA_20Q2_no_haem_phen.cls",
        "HER2",
        [1.087973, -1.358492, -1.178614, -0.77898, 0.157222, 1.168224, -0.360195, 0.608629],
    ),
    ("GSE76137.cls", "Proliferating_Arrested", [1, 2, 1, 2, 1, 2]),
    ("LPS_phen.cls", "CNTRL_LPS", [1, 1, 1, 2, 2, 2]),
)
    cl = joinpath(DA, cl)
    nar, ro_, co_, ta_x_sa_x_nu = Nucleus.DataFrame.separate(Nucleus.CLS.read(cl))
    @test nar === "Target"
    @test ro_ == [ta]
    @test all(startswith("Sample "), co_)
    @test eltype(ta_x_sa_x_nu) === eltype(re)
    @test view(ta_x_sa_x_nu, 1, eachindex(re)) == re
    # 387.667 μs (6234 allocations: 530.46 KiB)
    # 9.875 μs (98 allocations: 7.70 KiB)
    # 9.666 μs (98 allocations: 7.64 KiB)
    #@btime Nucleus.CLS.read($cl)
end


const DA = joinpath(Nucleus._DA, "GMT")
@test Nucleus.Path.read(DA) == ["c2.all.v7.1.symbols.gmt", "h.all.v7.1.symbols.gmt"]
for (gm, re) in (("h.all.v7.1.symbols.gmt", 50), ("c2.all.v7.1.symbols.gmt", 5529))
    gm = joinpath(DA, gm)
    @test length(Nucleus.GMT.read(gm)) === re
    # 336.375 μs (7717 allocations: 974.86 KiB)
    # 26.149 ms (515903 allocations: 51.28 MiB)
    #@btime Nucleus.GMT.read($gm)
end



# ---- #
const DA = joinpath(Nucleus._DA, "GCT")
@test Nucleus.Path.read(DA) == ["a.gct"]
const GC = joinpath(DA, "a.gct")
@test size(Nucleus.GCT.read(GC)) === (13321, 190)
# 97.734 ms (71117 allocations: 23.64 MiB)
#@btime Nucleus.GCT.read(GC);















using OrderedCollections: OrderedDict

# ---- #

const DA = joinpath(Nucleus._DA, "Dict")

# ---- #

@test Nucleus.Path.read(DA) == [
    "c2.all.v7.1.symbols.gmt",
    "example.toml",
    "example_1.json",
    "example_2.json",
    "gene_x_statistic_x_number.tsv",
]

# ---- #

const KES_VA = Dict("Existing" => 1)

# ---- #

for (ke, va, re) in (
    (
        "Existing",
        2,
        Dict("Existing" => 1, "Existing.2" => 2, "Existing.3" => 2, "Existing.4" => 2),
    ),
    ("New", 2, Dict("Existing" => 1, "New" => 2, "New.2" => 2, "New.3" => 2)),
)

    ke_va = copy(KES_VA)

    for _ in 1:3

        Nucleus.Dict.set_with_suffix!(ke_va, ke, va)

    end

    @test ke_va == re

    # 138.633 μs (4121 allocations: 275.38 KiB)
    # 138.961 μs (4112 allocations: 270.93 KiB)
    #@btime Nucleus.Dict.set_with_suffix!(ke_va, $ke, $va) setup = (ke_va = copy(KES_VA)) evals =
    #    1000

end

# ---- #

for (ke1_va, ke2_va, re) in (
    (Dict(1 => 'a'), Dict(2 => 'b'), Dict{Int, Char}),
    (Dict(1.0 => 'a'), Dict(2 => "Bb"), Dict{Float64, Any}),
    (Dict(1 => "Aa"), Dict(2.0 => view("Bb", 1:2)), Dict{Float64, AbstractString}),
)

    @test typeof(Nucleus.Dict.merge(ke1_va, ke2_va)) === re

end

# ---- #

const KE1_VA = Dict("1A" => 1, "B" => Dict("C" => 1, "1D" => 1))

# ---- #

const KE2_VA = Dict("2A" => 2, "B" => Dict("C" => 2, "2D" => 2))

# ---- #

for (ke1_va, ke2_va, re) in (
    (
        KE1_VA,
        KE2_VA,
        Dict("1A" => 1, "2A" => 2, "B" => Dict("C" => 2, "1D" => 1, "2D" => 2)),
    ),
    (
        KE2_VA,
        KE1_VA,
        Dict("1A" => 1, "2A" => 2, "B" => Dict("C" => 1, "1D" => 1, "2D" => 2)),
    ),
)

    @test Nucleus.Dict.merge(ke1_va, ke2_va) == re

    # 1.679 μs (32 allocations: 2.86 KiB)
    # 1.654 μs (32 allocations: 2.86 KiB)
    #@btime Nucleus.Dict.merge($ke1_va, $ke2_va)

end

# ---- #

const AN1_ = ['A', 'B']

# ---- #

for (an_id, re) in (
    (Dict(), BitVector()),
    (Dict('A' => 1), [true]),
    (Dict('B' => 1), [true]),
    (Dict('Z' => 1), [false]),
    (Dict('A' => 1, 'B' => 2, 'Z' => 3), [true, true, false]),
    (Dict('A' => 1, 'Z' => 2, 'B' => 3), [true, false, true]),
)

    is_ = Nucleus.Dict.is_in(an_id, AN1_)

    @test typeof(is_) === BitVector

    @test is_ == re

    # 37.387 ns (2 allocations: 96 bytes)
    # 33.157 ns (2 allocations: 96 bytes)
    # 33.199 ns (2 allocations: 96 bytes)
    # 32.235 ns (2 allocations: 96 bytes)
    # 32.277 ns (2 allocations: 96 bytes)
    # 32.109 ns (2 allocations: 96 bytes)
    #@btime Nucleus.Dict.is_in($an_id, AN1_)

end

# ---- #

const FE_ = reverse!(Nucleus.DataFrame.read(
    joinpath(DA, "gene_x_statistic_x_number.tsv");
    select = [1],
)[
    !,
    1,
],)

# ---- #

const FE1_ =
    Nucleus.GMT.read(joinpath(DA, "c2.all.v7.1.symbols.gmt"))["COLLER_MYC_TARGETS_UP"]

# ---- #

# 695.792 μs (6 allocations: 6.91 KiB)
#@btime in(FE1_).(FE_);
# 689.250 μs (3 allocations: 6.84 KiB)
#@btime in($FE1_).(FE_);
#@btime in(FE1_).($FE_);
#@btime in($FE1_).($FE_);

# ---- #

const FE1S = Set(FE1_)

# ---- #

# 441.500 ns (7 allocations: 1.13 KiB)
#@btime Set(FE1_);

# ---- #

# 484.334 μs (6 allocations: 6.91 KiB)
#@btime in(FE1S).(FE_);
# 461.583 μs (3 allocations: 6.84 KiB)
#@btime in($FE1S).(FE_);
#@btime in(FE1S).($FE_);
#@btime in($FE1S).($FE_);

# ---- #

const FE_ID = Dict(fe => id for (id, fe) in enumerate(FE_))

# ---- #

# 510.875 μs (7 allocations: 800.92 KiB)
#@btime Dict(fe => id for (id, fe) in enumerate(FE_));

# ---- #

# 362.577 ns (2 allocations: 2.66 KiB)
#@btime Nucleus.Dict.is_in(FE_ID, FE1_);

# ---- #

const JS1 = joinpath(DA, "example_1.json")

# ---- #

for ty in (Dict{String, Any}, OrderedDict{String, Any}, OrderedDict{String, String})

    @test typeof(Nucleus.Dict.read(JS1, ty)) === ty

end

# ---- #

for (fi, ty, re) in (
    (
        JS1,
        OrderedDict{String, Any},
        Dict("fruit" => "Apple", "color" => "Red", "size" => "Large"),
    ),
    (
        joinpath(DA, "example_2.json"),
        OrderedDict{String, Any},
        Dict("quiz" => Dict(
            "sport" => Dict("q1" => Dict(
                "options" => [
                    "New York Bulls",
                    "Los Angeles Kings",
                    "Golden State Warriros",
                    "Huston Rocket",
                ],
                "question" => "Which one is correct team name in NBA?",
                "answer" => "Huston Rocket",
            ),),
            "maths" => Dict(
                "q1" => Dict(
                    "options" => ["10", "11", "12", "13"],
                    "question" => "5 + 7 = ?",
                    "answer" => "12",
                ),
                "q2" => Dict(
                    "options" => ["1", "2", "3", "4"],
                    "question" => "12 - 8 = ?",
                    "answer" => "4",
                ),
            ),
        ),),
    ),
    (
        joinpath(DA, "example.toml"),
        Dict{String, Any},
        Dict(
            "servers" => Dict(
                "alpha" => Dict("dc" => "eqdc10", "ip" => "10.0.0.1"),
                "beta" => Dict("dc" => "eqdc10", "ip" => "10.0.0.2"),
            ),
            "clients" =>
                Dict("hosts" => ["alpha", "omega"], "data" => [["gamma", "delta"], [1, 2]]),
            "owner" => Dict("name" => "Tom Preston-Werner"),
            "title" => "TOML Example",
            "database" => Dict(
                "server" => "192.168.1.1",
                "connection_max" => 5000,
                "ports" => [8000, 8001, 8002],
                "enabled" => true,
            ),
        ),
    ),
)

    ke_va = Nucleus.Dict.read(fi)

    @test typeof(ke_va) === ty

    @test ke_va == re

end

# ---- #

const JSW = joinpath(Nucleus.TE, "write_read.json")

# ---- #

const KEW_VA = Dict(
    "Luffy" => "Pirate King",
    "Crews" => [
        "Luffy",
        "Zoro",
        "Nami",
        "Usopp",
        "Sanji",
        "Chopper",
        "Robin",
        "Franky",
        "Brook",
        "Jinbe",
    ],
    "episode" => 1030,
)

# ---- #

Nucleus.Dict.write(JSW, KEW_VA)

# ---- #

const KEWR_VA = Nucleus.Dict.read(JSW)

# ---- #

@test typeof(KEW_VA) != typeof(KEWR_VA)

# ---- #

@test KEW_VA == KEWR_VA
