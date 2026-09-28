#!/usr/bin/env python3

import csv
import sys


def extract_quoted(line):
    s = line.split(maxsplit=1)[1].strip()
    if s.startswith('"'):
        s = s[1:]
    if s.endswith('"'):
        s = s[:-1]
    return s


def parse_line(line):
    """Parses fixed-width domain architecture lines."""
    if not line.strip() or len(line) < 23:
        return None

    try:
        node_id = int(line[0:4])
        child1 = int(line[7:11])
        child2 = int(line[11:15])
        unit_size = int(line[15:19])
        num_segs = int(line[19:23])

        seg_str = line[23:]
        raw_nums = []

        for i in range(0, num_segs * 2 * 4, 4):
            chunk = seg_str[i : i + 4]
            if chunk.strip():
                raw_nums.append(int(chunk))

        segments = []
        for i in range(0, len(raw_nums) - 1, 2):
            segments.append((raw_nums[i], raw_nums[i + 1]))

        return {
            "node_id": node_id,
            "child1": child1,
            "child2": child2,
            "unit_size": unit_size,
            "num_segs": num_segs,
            "segments": segments,
            "raw_line": line.rstrip(),
        }
    except (ValueError, IndexError):
        return None


def get_partition_leaves(nodes, min_size=80):
    """Recursively finds partition nodes starting from root (node 1).

    Stops splitting a node if either of its children is smaller than min_size.
    """
    leaves = []

    def traverse(node_id):
        node = nodes.get(node_id)
        if not node:
            return

        c1_id = node["child1"]
        c2_id = node["child2"]

        if c1_id == 0 and c2_id == 0:
            leaves.append(node)
            return

        c1_node = nodes.get(c1_id)
        c2_node = nodes.get(c2_id)

        if not c1_node or not c2_node:
            leaves.append(node)
            return

        m = 2 * c1_node["unit_size"] * c2_node["unit_size"] / node["unit_size"]
        if m < min_size:
        #if c1_node["unit_size"] < min_size or c2_node["unit_size"] < min_size:
            leaves.append(node)
        else:
            if c1_node["unit_size"] >= min_size: traverse(c1_id)
            if c2_node["unit_size"] >= min_size: traverse(c2_id)

    traverse(1)
    return leaves


def extract_domain_string(nodes, default_length):
    """Computes the domain string from parsed node data."""
    if not nodes:
        return f"1-{default_length}" if default_length else ""

    partition_leaves = get_partition_leaves(nodes, min_size=80)

    domains = []
    for node in sorted(partition_leaves, key=lambda x: x["node_id"]):
        seg_fmt = "_".join([f"{s}-{e}" for s, e in node["segments"]])
        domains.append(seg_fmt)

    domain_string = ", ".join(domains)
    if len(domain_string) < 1:
        domain_string = f"1-{default_length}" if default_length else ""

    return domain_string


def parse_dat(dat_path):
    """Parses metadata headers and domain structure from a .dat file."""
    sequence = None
    dssp = None
    compnd = None
    length = None
    query_id = None
    nodes = {}

    with open(dat_path, "r") as f:
        for line in f:
            line_str = line.strip()

            if line_str.startswith(">>>>") and query_id is None:
                parts = line_str.split()
                if len(parts) >= 3:
                    query_id = parts[1]
                    try:
                        length = int(parts[2])
                    except ValueError:
                        pass

            elif line_str.startswith("-sequence"):
                sequence = extract_quoted(line_str)

            elif line_str.startswith("-dssp"):
                dssp = extract_quoted(line_str)

            elif line_str.startswith("-compnd"):
                compnd = extract_quoted(line_str)

            elif not line_str.startswith("-") and not line_str.startswith(
                ">>>>"
            ):
                parsed = parse_line(line)
                if parsed:
                    nodes[parsed["node_id"]] = parsed

    domain_string = extract_domain_string(nodes, length)

    return {
        "id": query_id,
        "length": length,
        "sequence": sequence,
        "dssp": dssp,
        "compnd": compnd,
        "domain_string": domain_string,
    }


def write_tsv(meta, out_path):
    """Writes metadata dictionary to a TSV file including domain_string."""
    fieldnames = ["id", "length", "sequence", "dssp", "compnd", "domain_string"]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(meta)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <input_dat_file> <output_tsv_file>")
        sys.exit(1)

    dat_file = sys.argv[1]
    out_file = sys.argv[2]

    meta = parse_dat(dat_file)
    write_tsv(meta, out_file)
