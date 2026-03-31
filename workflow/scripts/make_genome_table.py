import os

SUFFIXES = (".faa.gz", ".fasta.gz", ".fa.gz", ".faa", ".fasta", ".fa")


def file_code(path: str) -> str:
    name = os.path.basename(path)
    for suffix in SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def main() -> None:
    output_path = snakemake.output[0]
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    seen = set()
    resolved_paths = [
        full_path
        for input_path in snakemake.input
        if (full_path := os.path.normpath(os.path.abspath(input_path))) not in seen
        and not seen.add(full_path)
    ]

    with open(output_path, 'w', encoding='utf-8', newline='') as out_file:
        out_file.writelines(f"{path}\t{file_code(path)}\n" for path in resolved_paths)


if __name__ == '__main__':
    main()
