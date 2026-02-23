import pandas as pd


def parse_bed(log_file, out):
    print("chrom\tleft_start\tleft_end\tright_start\tright_end", file=out)
    for data_primers in pd.read_csv(
        snakemake.input[0],
        sep="\t",
        header=None,
        chunksize=chunksize,
        usecols=[0, 1, 2, 5],
    ):
        invalid_mask = ~data_primers[5].isin(["+", "-"])
        if invalid_mask.any():
            for row_id in data_primers.index[invalid_mask]:
                print("Invalid strand in row {}".format(row_id), file=log_file)

        fwd = data_primers[data_primers[5] == "+"]
        rev = data_primers[data_primers[5] == "-"]

        fwd_out = fwd[[0]].copy()
        fwd_out["left_start"] = fwd[1] + 1
        fwd_out["left_end"] = fwd[2]
        fwd_out["right_start"] = -1
        fwd_out["right_end"] = -1

        rev_out = rev[[0]].copy()
        rev_out["left_start"] = -1
        rev_out["left_end"] = -1
        rev_out["right_start"] = rev[1] + 1
        rev_out["right_end"] = rev[2]

        pd.concat([fwd_out, rev_out]).sort_index().to_csv(
            out, sep="\t", index=False, header=False
        )


def parse_bedpe(log_file, out_path):
    first_chunk = True
    for data_primers in pd.read_csv(
        snakemake.input[0],
        sep="\t",
        header=None,
        chunksize=chunksize,
        usecols=[0, 1, 2, 3, 4, 5],
    ):
        valid_primers = data_primers[0] == data_primers[3]
        valid_data = data_primers[valid_primers].copy()
        valid_data.iloc[:, [1, 4]] += 1
        valid_data.drop(columns=[3], inplace=True)
        valid_data.dropna(how="all", inplace=True)

        valid_data.to_csv(
            out_path,
            sep="\t",
            index=False,
            header=["chrom", "left_start", "left_end", "right_start", "right_end"]
            if first_chunk
            else False,
            mode="w" if first_chunk else "a",
        )
        first_chunk = False

        print(
            data_primers[~valid_primers].to_csv(sep="\t", index=False, header=False),
            file=log_file,
            end="",
        )


chunksize = 10**6
with open(snakemake.log[0], "w") as log_file:
    if snakemake.input[0].endswith("bedpe"):
        parse_bedpe(log_file, snakemake.output[0])
    else:
        with open(snakemake.output[0], "w") as out:
            parse_bed(log_file, out)