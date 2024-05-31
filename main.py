import numpy as np
import pysam
from pathlib import Path
import sys


def read_bed(file_path: Path):
    data = np.loadtxt(
        file_path,
        dtype={"names": ("chr", "start", "stop"), "formats": ("U50", "i4", "i4")},
        delimiter="\t",
    )
    chrs = data["chr"]
    starts = data["start"]
    stops = data["stop"]
    region_lengths = stops - starts
    depths = [np.zeros(length) for length in region_lengths]

    return chrs, starts, stops, depths


def process_bedgraph(file_path: str, chrs, starts, stops, depths):
    tbx = pysam.TabixFile(file_path, index=file_path + ".csi")

    for i in range(len(chrs)):
        chr = chrs[i]
        start = starts[i]
        stop = stops[i]
        region_array = depths[i]

        region_length = stop - start

        rows = tbx.fetch(chr, start, stop, parser=pysam.asBed())
        for r in rows:
            norm_start = max(0, r.start - start)
            norm_end = min(region_length, r.end - start)
            region_array[norm_start:norm_end] += int(r.name)


def process_multiple_bedgraphs(bedgraph_files: list, chrs, starts, stops, depths):
    for bedgraph_file in bedgraph_files:
        process_bedgraph(bedgraph_file, chrs, starts, stops, depths)
    return chrs, starts, stops, depths


# Example usage
chrs, starts, stops, depths = read_bed(sys.argv[1])
bedgraph_files = Path(sys.argv[2]).glob("*.per-base.bed.gz")
chrs, starts, stops, depths = process_multiple_bedgraphs(
    bedgraph_files, chrs, starts, stops, depths
)

# Display the results
for i in range(len(chrs)):
    print(f"Region {i}: {chrs[i]}:{starts[i]}-{stops[i]}, Depth: {depths[i]}")
