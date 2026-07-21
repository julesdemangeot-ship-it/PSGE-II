"""PSGE-II v1.2 validation subpackage.

Analytical oracles, intrinsic reference meshes, and a runnable validation
campaign. The campaign is exposed programmatically as ``run_campaign`` and as a
module entry point:

    python -m psge.validation.suite
"""

from psge.validation.oracles_v12 import (
    cone_deficit,
    regular_cone,
    regular_tetrahedron,
    flat_fan,
)
from psge.validation.campaign import run_campaign

__all__ = [
    # analytical oracles
    "cone_deficit",
    # intrinsic reference meshes
    "regular_cone",
    "regular_tetrahedron",
    "flat_fan",
    # campaign runner
    "run_campaign",
]
