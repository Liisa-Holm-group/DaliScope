import infomap
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Union, Any

from dataclasses import dataclass
from typing import Iterator, Set, Tuple


@dataclass(frozen=True, slots=True)
class Segment:
    """
    Represents a continuous Secondary Structure Element (SSE) segment.
    
    Attributes:
        sse_type (str): SSE code ('H' for helix, 'E' for strand, 'C' for coil/loop, etc.)
        start (int): 0-based or 1-based start residue index (inclusive)
        end (int): Start residue index (inclusive)
    """
    sse_type: str
    start: int
    end: int

    def __post_init__(self):
        """Validates boundaries upon instantiation."""
        if self.end < self.start:
            raise ValueError(
                f"Invalid Segment bounds: end ({self.end}) cannot be less than start ({self.start})"
            )

    # ------------------------------------------------------------------
    # Properties & Derived Views (Zero String-Parsing Required)
    # ------------------------------------------------------------------

    @property
    def id(self) -> str:
        """Canonical string identifier used for registry keys and provenance logs."""
        return f"seg_{self.sse_type}_{self.start}_{self.end}"

    @property
    def length(self) -> int:
        """Number of residues contained in this segment."""
        return self.end - self.start + 1

    @property
    def range(self) -> Tuple[int, int]:
        """Returns bounds as a (start, end) tuple."""
        return self.start, self.end

    @property
    def residue_indices(self) -> range:
        """Returns a python range of the residue indices [start, end]."""
        return range(self.start, self.end + 1)

    @property
    def residue_ids(self) -> Set[str]:
        """Set of formatted residue IDs ('r_12', 'r_13', etc.)."""
        return {f"r_{i}" for i in range(self.start, self.end + 1)}

    # ------------------------------------------------------------------
    # Constructors & Utility Methods
    # ------------------------------------------------------------------

    @classmethod
    def from_tuple(cls, seg_tuple: Tuple[str, int, int]) -> "Segment":
        """Constructs a Segment instance from a ('H', start, end) tuple."""
        return cls(sse_type=str(seg_tuple[0]), start=int(seg_tuple[1]), end=int(seg_tuple[2]))

    @classmethod
    def from_id(cls, seg_id: str) -> "Segment":
        """
        Parses a 'seg_H_720_745' string ID back into a Segment instance.
        Useful when deserializing pipeline results or reading legacy records.
        """
        parts = seg_id.split("_")
        if len(parts) != 4 or parts[0] != "seg":
            raise ValueError(f"Invalid segment ID format: '{seg_id}'. Expected 'seg_TYPE_START_END'.")
        return cls(sse_type=parts[1], start=int(parts[2]), end=int(parts[3]))

    def to_tuple(self) -> Tuple[str, int, int]:
        """Exports back to traditional ('H', start, end) tuple format."""
        return self.sse_type, self.start, self.end

    def contains_residue(self, res_idx: int) -> bool:
        """Fast check if a residue index falls within this segment."""
        return self.start <= res_idx <= self.end

    def overlaps_with(self, other: "Segment") -> bool:
        """Checks if this segment overlaps in sequence space with another segment."""
        return max(self.start, other.start) <= min(self.end, other.end)

    def __iter__(self) -> Iterator[int]:
        """Allows iterating directly over residue indices: `for r in segment:`"""
        return iter(range(self.start, self.end + 1))

    def __repr__(self) -> str:
        return f"Segment({self.sse_type}, {self.start}..{self.end}, len={self.length})"


class StructuralRegistry:
    def __init__(self):
        # Maps string ID -> rich Segment object
        self.segments: dict[str, Segment] = {}

        # Core membership mappings
        self.seg_to_res: dict[str, set[str]] = defaultdict(set)
        self.res_to_seg: dict[str, str] = {}

        self._id_to_level: Dict[str, str] = {}

        # Parent-Child mapping tables
        self.sheet_to_seg: Dict[str, Set[str]] = defaultdict(set)
        self.comm_to_res: Dict[str, Set[str]] = defaultdict(set)

        # Reverse mappings (Child -> Parent)
        self.seg_to_sheet: Dict[str, str] = {}
        self.res_to_comm: Dict[str, str] = {}
        self.comm_to_domain: Dict[str, str] = {}

    def register_segment(self, sse_type: str, start: int, end: int) -> Segment:
        seg = Segment(type=sse_type, start=start, end=end)
        seg_id = seg.id

        self.segments[seg_id] = seg
        self._id_to_level[seg_id] = "segment"

        for r_idx in range(start, end + 1):
            r_id = f"r_{r_idx}"
            self.seg_to_res[seg_id].add(r_id)
            self.res_to_seg[r_id] = seg_id

        return seg

    def get_segment_bounds(self, seg_id: str) -> tuple[int, int]:
        """Zero-parsing property access."""
        seg = self.segments[seg_id]
        return seg.start, seg.end

    def build_from_pipeline(
        self,
        sequence: str,
        segments: List[Tuple[str, int, int]],
        sheets: List[Any],
        hierarchy: Dict[int, Tuple[int, ...]]
    ):
        """Populates the registry from raw pipeline outputs."""

        # 1. Register Residues and Segments (L0 & L1)
        for sse_type, start, end in segments:
            seg_id = f"seg_{sse_type}_{start}_{end}"
            self._id_to_level[seg_id] = "segment"

            for res_idx in range(start, end + 1):
                res_id = f"r_{res_idx}"
                self._id_to_level[res_id] = "residue"

                self.seg_to_res[seg_id].add(res_id)
                self.res_to_seg[res_id] = seg_id

        # 2. Register Sheets (L2)
        for sheet_idx, strand_list in enumerate(sheets):
            sheet_id = f"sheet_{sheet_idx:02d}"
            self._id_to_level[sheet_id] = "sheet"

            # Extract (start, end) ranges from strands in this sheet
            strand_ranges = _extract_sheet_strand_tuples(strand_list)

            for s_start, s_end in strand_ranges:
                for seg_type, seg_start, seg_end in segments:
                    if seg_start == s_start and seg_end == s_end:
                        seg_id = f"seg_{seg_type}_{s_start}_{s_end}"
                        self.sheet_to_seg[sheet_id].add(seg_id)
                        self.seg_to_sheet[seg_id] = sheet_id

        # 3. Register Infomap Hierarchical Communities (L3 & L4)
        for res_idx, path in hierarchy.items():
            res_id = f"r_{res_idx}"

            comm_path_str = "_".join(map(str, path))
            comm_id = f"comm_{comm_path_str}"
            self._id_to_level[comm_id] = "community"

            domain_id = f"dom_{path[0]}"
            self._id_to_level[domain_id] = "domain"

            self.comm_to_res[comm_id].add(res_id)
            self.res_to_comm[res_id] = comm_id
            self.comm_to_domain[comm_id] = domain_id

    def get_residues(self, entity_id: str) -> Set[str]:
        """Resolves ANY entity down to its raw residue members (L0)."""
        level = self._id_to_level.get(entity_id)

        if level == "residue":
            return {entity_id}
        elif level == "segment":
            return self.seg_to_res.get(entity_id, set())
        elif level == "sheet":
            res = set()
            for seg_id in self.sheet_to_seg.get(entity_id, set()):
                res.update(self.seg_to_res[seg_id])
            return res
        elif level == "community":
            return self.comm_to_res.get(entity_id, set())
        elif level == "domain":
            res = set()
            for comm_id, dom_id in self.comm_to_domain.items():
                if dom_id == entity_id:
                    res.update(self.comm_to_res[comm_id])
            return res
        raise ValueError(f"Unknown entity ID: {entity_id}")

    def get_parents(self, entity_id: str) -> Dict[str, Any]:
        """Gets all higher-level containers for a given entity."""
        res_set = self.get_residues(entity_id)
        if not res_set:
            return {"segment": None, "sheet": None, "community": None, "domain": None}

        sample_res = next(iter(res_set))

        seg = self.res_to_seg.get(sample_res)
        sheet = self.seg_to_sheet.get(seg) if seg else None
        comm = self.res_to_comm.get(sample_res)
        dom = self.comm_to_domain.get(comm) if comm else None

        return {
            "segment": seg,
            "sheet": sheet,
            "community": comm,
            "domain": dom
        }


def _extract_sheet_strand_tuples(strand_list: Any) -> List[Tuple[int, int]]:
    """Helper to extract (start, end) tuples from various sheet/strand data formats."""
    ranges = []
    for item in strand_list:
        if isinstance(item, tuple) and len(item) == 2 and isinstance(item[0], int):
            ranges.append(item)
        elif hasattr(item, 'start') and hasattr(item, 'end'):
            ranges.append((item.start, item.end))
        elif isinstance(item, (list, tuple)) and len(item) == 2 and isinstance(item[1], (tuple, list)):
            ranges.append(item[1])
    return ranges
