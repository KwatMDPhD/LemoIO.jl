module LeMoIO

using CSV: read as CSV_read, write as CSV_write

using CodecZlib: GzipDecompressor, transcode

using DataFrames: AbstractDataFrame, DataFrame, insertcols!

using Mmap: mmap

using XLSX: readtable

using JSON: parsefile as json_parsefile, print

using OrderedCollections: OrderedDict

using TOML: parsefile as toml_parsefile

function read(fi, dicttype = OrderedDict; ke_ar...)

    if splitext(fi) == ".toml"

        toml_parsefile(fi; ke_ar...)

    else

        json_parsefile(fi; dicttype, ke_ar...)

    end

end

function write(js, ke_va, id = 2)

    open(js, "w") do io

        return print(io, ke_va, id)

    end

end

function read_table(fi; ke_ar...)

    it_ = mmap(fi)

    if fi[(end - 2):end] == ".gz"

        it_ = transcode(GzipDecompressor, it_)

    end

    return CSV_read(it_, DataFrame; ke_ar...)

end

function read_table(xl, sh; ke_ar...)

    @error readtable(xl, sh; ke_ar...) typeof(readtable(xl, sh; ke_ar...))
    return DataFrame(readtable(xl, sh; ke_ar...))

end

function write(ts, da::AbstractDataFrame)

    return CSV_write(ts, da; delim = '\t')

end

function write(ts, nr, ro_, nc, co_, ma)

    return write(ts, insertcols!(DataFrame(ma, co_), 1, nr => nr_))

end

function _make_matrix(fu, st_)

    return [fu(st) for _ in 1:1, st in st_]

end

function read(cl)

    li1, li2, li3 = readlines(cl)

    li22 = view(li2, 2:lastindex(li2))

    li3_ = split(li3)

    n_sa3 = lastindex(li3_)

    nar = "Target"

    co_ = (id -> "Sample $id").(1:n_sa3)

    if li1 == "#numeric"

        Nucleus.DataFrame.make(nar, li22, co_, _make_matrix(st -> parse(Float64, st), li3_))

    else

        li1_ = split(li1)

        n_sa1 = parse(Int, li1_[1])

        if n_sa1 != n_sa3

            error("Numbers of samples differ. $n_sa1 (line 1) != $n_sa3 (line 3).")

        end

        n_gr1 = parse(Int, li1_[2])

        gr_ = split(li22)

        n_gr2 = lastindex(gr_)

        n_gr3 = lastindex(unique(li3_))

        if !(n_gr1 == n_gr2 == n_gr3)

            error("Numbers of groups differ. !($n_gr1 (line 1) == $n_gr2 (line 2) == $n_gr3 (line 3)).",)

        end

        gr_id = Dict(gr => id for (id, gr) in enumerate(gr_))

        Nucleus.DataFrame.make(
            nar,
            join(gr_, '_'),
            co_,
            _make_matrix(st -> gr_id[st], li3_),
        )

    end

end

function read(gm)

    se_ge_ = Dict{String, Vector{String}}()

    for li in eachline(gm)

        sp_ = split(li, '\t')

        se = sp_[1]

        Nucleus.Error.error_has_key(se_ge_, se)

        se_ge_[se] = filter!(!isempty, view(sp_, 3:lastindex(sp_)))

    end

    return se_ge_

end

function read(gc)
    return Nucleus.DataFrame.read(gc; header = 3, drop = ["Description"])
end

end
