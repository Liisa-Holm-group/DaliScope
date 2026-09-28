"""Create a DaliScope data pack using explicitly configured local databases."""

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

NAMESPACES = {"PDB": 1, "AFDB2": 2, "BFVD": 3, "VIRO3D": 4}


def parse_namespace(value):
    text = str(value).upper()
    if text in NAMESPACES:
        return NAMESPACES[text]
    if text in {"1", "2", "3", "4"}:
        return int(text)
    raise argparse.ArgumentTypeError("namespace must be PDB, AFDB2, BFVD, VIRO3D, or 1-4")


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tsvfile", type=Path, help="DALI alignment output (TSV or text)")
    parser.add_argument("ns_id", type=parse_namespace, help="Target database namespace")
    parser.add_argument("testpack_name", help="Output basename, without .tar.gz")
    parser.add_argument("pdbfile", type=Path, help="Query PDB structure, optionally gzipped")
    parser.add_argument("chainid", help="Query chain identifier")
    parser.add_argument("datfile", type=Path, help="Query DALI .dat metadata")
    parser.add_argument("--ledger-db", type=Path, required=True, help="Protein SQLite ledger")
    parser.add_argument("--shard-dir", type=Path, required=True, help="Directory of protein NPZ shards")
    parser.add_argument("--pfam-db", type=Path, required=True, help="Pfam annotation SQLite database")
    parser.add_argument("--pfam-names", type=Path, required=True, help="Full Pfam description TSV")
    parser.add_argument("--ligand-db", type=Path, help="Optional PDB ligand SQLite database")
    parser.add_argument("--output-dir", type=Path, default=Path.cwd())
    parser.add_argument("--overwrite", action="store_true", help="Replace an existing output pack")
    parser.add_argument("--upload-to", help="Optional scp destination, such as user@host:/directory/")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if Path(args.testpack_name).name != args.testpack_name or args.testpack_name in {"", ".", ".."}:
        parser.error("testpack_name must be a filename basename")
    for field in ("tsvfile", "pdbfile", "datfile", "ledger_db", "pfam_db", "pfam_names", "ligand_db"):
        value = getattr(args, field)
        if value is not None:
            value = value.resolve()
            if not value.is_file():
                parser.error(f"{field}: file not found: {value}")
            setattr(args, field, value)
    args.shard_dir = args.shard_dir.resolve()
    if not args.shard_dir.is_dir():
        parser.error(f"shard directory not found: {args.shard_dir}")
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / f"{args.testpack_name}.tar.gz"
    if destination.exists() and not args.overwrite:
        parser.error(f"output already exists: {destination}; choose a new name or use --overwrite")
    scripts = Path(__file__).resolve().parent

    with tempfile.TemporaryDirectory(prefix="daliscope-pack-", dir=output_dir) as temp:
        work = Path(temp)

        def run(script, *arguments):
            print(f"Running {script}", flush=True)
            subprocess.run([sys.executable, str(scripts / script), *map(str, arguments)], cwd=work, check=True)

        run("minimal_parser.py", args.tsvfile, args.ns_id)
        run("colab_package.py", args.ns_id, "--database", args.ledger_db, "--shards", args.shard_dir)
        run("pfam_query.py", args.ns_id, "--database", args.pfam_db)
        run("renumber_pdb.py", args.pdbfile, args.chainid, "Query.pdb")
        run("query_meta.py", args.datfile, "query_meta.tsv")
        run("mini_pfam_names.py", args.pfam_names, "mini_pfam_names.tsv")
        required = ["summary.tsv", "id_mapping.tsv", "segments.tsv", "pfam.tsv",
                    "colab_pack.npz", "colab_pack_index.tsv", "mini_pfam_names.tsv",
                    "Query.pdb", "query_meta.tsv"]
        if args.ligand_db:
            if args.ns_id != 1:
                parser.error("--ligand-db currently supports the PDB namespace only")
            run("extract_mini_ligands.py", "colab_pack_index.tsv", args.ligand_db,
                "mini_ligand.tsv", "mini_ligand_names.tsv")
            required += ["mini_ligand.tsv", "mini_ligand_names.tsv"]
        missing = [name for name in required if not (work / name).is_file()]
        if missing:
            raise FileNotFoundError(f"pack generation did not produce: {', '.join(missing)}")
        archive = work / "result.tar.gz"
        with tarfile.open(archive, "w:gz") as tar:
            for name in required:
                tar.add(work / name, arcname=name)
        os.replace(archive, destination)
    print(f"Created {destination}")
    if args.upload_to:
        subprocess.run(["scp", str(destination), args.upload_to], check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
