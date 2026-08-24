#!/usr/bin/env python3

"""Generate the triangular prism mesh used by the sine-wave study.

This file is intentionally a standalone Gmsh-Python entry point.  The normal
system Python in this project does not import ``gmsh``; ``Allrun`` invokes
this file with the VibeFlow Python environment that provides Gmsh 4.15.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import gmsh


EPS = 1.0e-10


def _physical_group(dim: int, entities: list[int], name: str) -> None:
    if not entities:
        raise RuntimeError(f"No Gmsh entities found for physical group {name}")
    tag = gmsh.model.addPhysicalGroup(dim, entities)
    gmsh.model.setPhysicalName(dim, tag, name)


def _classify_surface(surface: int, thickness: float) -> str:
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, surface)
    if abs(zmax - zmin) < EPS and abs(zmin) < EPS:
        return "zMinSource"
    if abs(zmax - zmin) < EPS and abs(zmin - thickness) < EPS:
        return "zMaxSource"
    if abs(xmax - xmin) < EPS and abs(xmin) < EPS:
        return "xMinSource"
    if abs(xmax - xmin) < EPS and abs(xmin - 1.0) < EPS:
        return "xMaxSource"
    if abs(ymax - ymin) < EPS and abs(ymin) < EPS:
        return "yMinSource"
    if abs(ymax - ymin) < EPS and abs(ymin - 1.0) < EPS:
        return "yMaxSource"
    raise RuntimeError(
        "Could not classify extruded surface "
        f"{surface} with bounds {(xmin, ymin, zmin, xmax, ymax, zmax)}"
    )


def _write_geometry_metadata(
    output: Path,
    resolution: int,
    thickness: float,
    z_min_surface: int,
) -> None:
    """Save the 2-D Gmsh surface geometry for inspection and plotting tools."""
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes(
        2, z_min_surface, includeBoundary=True
    )
    nodes = [
        [float(node_coords[3 * i]), float(node_coords[3 * i + 1])]
        for i in range(len(node_tags))
    ]
    node_index = {int(tag): i for i, tag in enumerate(node_tags)}
    element_types, _, element_nodes = gmsh.model.mesh.getElements(2, z_min_surface)
    triangles: list[list[int]] = []
    for element_type, connectivity in zip(element_types, element_nodes):
        _, dim, _, nodes_per_element, _, _ = gmsh.model.mesh.getElementProperties(
            element_type
        )
        if dim != 2 or nodes_per_element != 3:
            continue
        for offset in range(0, len(connectivity), nodes_per_element):
            triangles.append(
                [
                    node_index[int(connectivity[offset])],
                    node_index[int(connectivity[offset + 1])],
                    node_index[int(connectivity[offset + 2])],
                ]
            )

    output.write_text(
        json.dumps(
            {
                "resolution": resolution,
                "thickness": thickness,
                "nodes": nodes,
                "triangles": triangles,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def generate(case: Path, resolution: int, thickness: float) -> Path:
    """Generate ``case/mesh/mesh.msh`` and its geometry metadata."""
    if resolution <= 0:
        raise ValueError("resolution must be positive")
    if thickness <= 0.0:
        raise ValueError("thickness must be positive")

    output_dir = case / "mesh"
    output_dir.mkdir(parents=True, exist_ok=True)
    msh_path = output_dir / "mesh.msh"
    metadata_path = output_dir / "mesh_geometry.json"

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.option.setNumber("Mesh.RecombineAll", 0)
        gmsh.model.add(f"triangular_advection_N{resolution}")

        geo = gmsh.model.geo
        points = [
            geo.addPoint(0.0, 0.0, 0.0),
            geo.addPoint(1.0, 0.0, 0.0),
            geo.addPoint(1.0, 1.0, 0.0),
            geo.addPoint(0.0, 1.0, 0.0),
        ]
        lines = [
            geo.addLine(points[0], points[1]),
            geo.addLine(points[1], points[2]),
            geo.addLine(points[2], points[3]),
            geo.addLine(points[3], points[0]),
        ]
        loop = geo.addCurveLoop(lines)
        surface = geo.addPlaneSurface([loop])
        geo.synchronize()

        for line in lines:
            gmsh.model.mesh.setTransfiniteCurve(line, resolution + 1)
        gmsh.model.mesh.setTransfiniteSurface(
            surface, "Alternate", points
        )

        extruded = geo.extrude(
            [(2, surface)],
            0.0,
            0.0,
            thickness,
            numElements=[1],
            recombine=True,
        )
        geo.synchronize()

        volume_tags = [tag for dim, tag in extruded if dim == 3]
        if not volume_tags:
            raise RuntimeError("Gmsh extrusion did not create a volume")

        # The transfinite constraints on the original surface do not
        # automatically constrain every entity created by extrusion.  Apply
        # them again to the generated curves and surfaces so that N means
        # exactly N intervals in x and y and one interval in z.
        for _, curve in gmsh.model.getEntities(1):
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(
                1, curve
            )
            curve_nodes = resolution + 1 if abs(zmax - zmin) < EPS else 2
            gmsh.model.mesh.setTransfiniteCurve(curve, curve_nodes)

        for _, generated_surface in gmsh.model.getEntities(2):
            corner_tags = [
                tag
                for dim, tag in gmsh.model.getBoundary(
                    [(2, generated_surface)],
                    oriented=False,
                    recursive=True,
                )
                if dim == 0
            ]
            gmsh.model.mesh.setTransfiniteSurface(
                generated_surface, "Alternate", corner_tags
            )

        for volume in volume_tags:
            corner_tags = [
                tag
                for dim, tag in gmsh.model.getBoundary(
                    [(3, volume)],
                    oriented=False,
                    recursive=True,
                )
                if dim == 0
            ]
            gmsh.model.mesh.setTransfiniteVolume(volume, corner_tags)

        surfaces = [tag for dim, tag in gmsh.model.getEntities(2)]
        surface_names: dict[str, list[int]] = {}
        for tag in surfaces:
            name = _classify_surface(tag, thickness)
            surface_names.setdefault(name, []).append(tag)
        expected_names = {
            "zMinSource",
            "zMaxSource",
            "xMinSource",
            "xMaxSource",
            "yMinSource",
            "yMaxSource",
        }
        if set(surface_names) != expected_names:
            raise RuntimeError(
                f"Unexpected Gmsh boundary surfaces: {sorted(surface_names)}"
            )

        for name, tags in sorted(surface_names.items()):
            _physical_group(2, tags, name)
        _physical_group(3, volume_tags, "fluid")

        gmsh.model.mesh.generate(3)
        gmsh.write(str(msh_path))

        _write_geometry_metadata(
            metadata_path,
            resolution,
            thickness,
            surface_names["zMinSource"][0],
        )
    finally:
        gmsh.finalize()

    print(f"mesh={msh_path}")
    print(f"meshGeometry={metadata_path}")
    return msh_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--resolution", type=int, required=True)
    parser.add_argument("--thickness", type=float, required=True)
    args = parser.parse_args()
    generate(args.case.resolve(), args.resolution, args.thickness)


if __name__ == "__main__":
    main()
