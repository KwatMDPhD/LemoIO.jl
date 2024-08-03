module OmicsIO

using CSV: read as CSV_read, write as CSV_write

using CodecZlib: GzipDecompressor, transcode

using DataFrames: AbstractDataFrame, DataFrame, insertcols!

using Mmap: mmap

using XLSX: readtable

function read_dataframe(fi; ke_ar...)

    it_ = mmap(fi)

    if fi[(end - 2):end] == ".gz"

        it_ = transcode(GzipDecompressor, it_)

    end

    return CSV_read(it_, DataFrame; ke_ar...)

end

function read_dataframe(xl, sh; ke_ar...)

    @error readtable(xl, sh; ke_ar...) typeof(readtable(xl, sh; ke_ar...))
    return DataFrame(readtable(xl, sh; ke_ar...))

end

function write(ts, da::AbstractDataFrame)

    return CSV_write(ts, da; delim = '\t')

end

function write(ts, nr, ro_, nc, co_, ma)

    return write(ts, insertcols!(DataFrame(ma, co_), 1, nr => nr_))

end

end
