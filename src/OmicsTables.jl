module OmicsTables

using CSV: read as CSV_read, write as CSV_write

# TODO: Do not depend on `DataFrames`.
using DataFrames: DataFrame

using CodecZlib: GzipDecompressor, transcode

using Mmap: mmap

using XLSX: readtable

import Base: show

struct OmicsTable{T}

    row::Symbol

    rows::Vector{Symbol}

    column::Symbol

    columns::Vector{Symbol}

    name::Symbol

    matrix::Matrix{T}

end

function OmicsTable(ro, ro_, co_, ma; co = "Column", na = "Name")

    return OmicsTable(ro, ro_, co, co_, na, ma)

end

function OmicsTable(da; ke_ar...)

    na_ = names(da)

    id_ = 2:lastindex(na_)

    return OmicsTable(na_[1], da[:, 1], a2_[id_], Matrix(da[!, id_]); ke_ar...)

end

function _make_dataframe(ta)

    return insertcols!(DataFrame(ta.matrix, ta.columns), 1, ta.row => ta.rows)

end

function show(io::IO, ta::OmicsTable)

    return print(io, "$(ta.row)-x-$(ta.column)-x-$(ta.name)\n$(_make_dataframe(ta))")

end

function write(ts, ta)

    return CSV_write(ts, _make_dataframe(ta); delim = '\t')

end

function read_dataframe(fi; ke_ar...)

    # TODO: Test not making an empty file and then remove.
    Nucleus.Error.error_missing(fi)

    it_ = mmap(fi)

    if fi[(end - 2):end] == ".gz"

        it_ = transcode(GzipDecompressor, it_)

    end

    return CSV_read(it_, DataFrame; ke_ar...)

end

function read_dataframe(xl, sh; ke_ar...)

    return DataFrame(readtable(xl, sh; ke_ar...))

end

end
