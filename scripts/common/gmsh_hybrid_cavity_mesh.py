#!/usr/bin/env python3

"""Generate a square-cavity hybrid quad/triangle prism mesh with Gmsh."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import gmsh


def _physical(dim: int, entities: list[int], name: str) -> None:
    if not entities:
        raise RuntimeError(f"No entities found for physical group {name}")
    tag = gmsh.model.addPhysicalGroup(dim, entities)
    gmsh.model.setPhysicalName(dim, tag, name)


def _line_count(length: float, target_size: float) -> int:
    return max(1, int(round(length / target_size)))


def generate(
    case: Path,
    resolution: int,
    boundary_layer_thickness: float = 0.12,
    boundary_layer_count: int | None = None,
    thickness: float = 0.1,
) -> Path:
    if resolution < 4:
        raise ValueError("resolution must be at least 4")
    if not 0.0 < boundary_layer_thickness < 0.5:
        raise ValueError("boundary_layer_thickness must be in (0, 0.5)")
    if boundary_layer_count is None:
        boundary_layer_count = max(2, round(resolution * boundary_layer_thickness))

    target_size = 1.0 / resolution
    b = boundary_layer_thickness
    q = 1.0 - 2.0 * b
    outer_count = _line_count(1.0, target_size)
    corner_count = int(boundary_layer_count)
    core_count = outer_count - 2 * corner_count
    if core_count < 2:
        raise RuntimeError("Resolution is too low for the selected boundary layer")

    output_dir = case / "mesh"
    output_dir.mkdir(parents=True, exist_ok=True)
    msh_path = output_dir / "mesh.msh"
    metadata_path = output_dir / "mesh_metadata.json"

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.option.setNumber("Mesh.RecombineAll", 0)
        gmsh.model.add(f"hybrid_lid_cavity_N{resolution}")
        geo = gmsh.model.geo

        xs = [0.0, b, 1.0 - b, 1.0]
        ys = [0.0, b, 1.0 - b, 1.0]
        points = [[geo.addPoint(x, y, 0.0) for x in xs] for y in ys]

        horizontal = [[geo.addLine(points[j][i], points[j][i + 1])
                       for i in range(3)] for j in range(4)]
        vertical = [[geo.addLine(points[j][i], points[j + 1][i])
                     for j in range(3)] for i in range(4)]

        surfaces: list[list[int]] = []
        for j in range(3):
            row: list[int] = []
            for i in range(3):
                loop = geo.addCurveLoop(
                    [horizontal[j][i], vertical[i + 1][j],
                     -horizontal[j + 1][i], -vertical[i][j]]
                )
                row.append(geo.addPlaneSurface([loop]))
            surfaces.append(row)
        geo.synchronize()

        horizontal_counts = [corner_count, core_count, corner_count]
        vertical_counts = [boundary_layer_count, core_count, boundary_layer_count]
        for line in horizontal[0] + horizontal[3]:
            gmsh.model.mesh.setTransfiniteCurve(
                line, horizontal_counts[horizontal[0].index(line)]
                if line in horizontal[0] else horizontal_counts[horizontal[3].index(line)]
            )
        for row_index in (1, 2):
            for i, line in enumerate(horizontal[row_index]):
                gmsh.model.mesh.setTransfiniteCurve(line, horizontal_counts[i])
        for i in range(4):
            for j, line in enumerate(vertical[i]):
                count = vertical_counts[j] if i in (0, 3) else vertical_counts[j]
                gmsh.model.mesh.setTransfiniteCurve(line, count)

        quad_surfaces = [
            surfaces[0][0], surfaces[0][1], surfaces[0][2],
            surfaces[1][0], surfaces[1][2],
            surfaces[2][0], surfaces[2][1], surfaces[2][2],
        ]
        for surface in quad_surfaces:
            gmsh.model.mesh.setTransfiniteSurface(surface)
            gmsh.model.mesh.setRecombine(2, surface)

        # The central surface intentionally remains triangular.
        gmsh.model.mesh.setSize(
            gmsh.model.getBoundary([(2, surfaces[1][1])], recursive=True),
            target_size,
        )

        extruded = geo.extrude(
            [(2, surface) for row in surfaces for surface in row],
            0.0,
            0.0,
            thickness,
            numElements=[1],
            recombine=True,
        )
        geo.synchronize()

        volume_tags = [tag for dim, tag in extruded if dim == 3]
        if not volume_tags:
            raise RuntimeError("Gmsh extrusion did not create volumes")

        # All source 2-D surfaces are meshed explicitly; the center remains
        # triangular while the eight outer surfaces remain quadrilateral.
        for dim, tag in gmsh.model.getEntities(2):
            if tag in quad_surfaces:
                gmsh.model.mesh.setRecombine(dim, tag)

        gmsh.model.mesh.generate(3)

        zmin_faces: list[int] = []
        zmax_faces: list[int] = []
        moving_faces: list[int] = []
        left_faces: list[int] = []
        right_faces: list[int] = []
        bottom_faces: list[int] = []
        for dim, tag in gmsh.model.getEntities(2):
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
            if abs(zmin) < 1e-8 and abs(zmax) < 1e-8:
                zmin_faces.append(tag)
            elif abs(zmin - thickness) < 1e-8 and abs(zmax - thickness) < 1e-8:
                zmax_faces.append(tag)
            elif abs(xmin) < 1e-8 and abs(xmax) < 1e-8:
                left_faces.append(tag)
            elif abs(xmin - 1.0) < 1e-8 and abs(xmax - 1.0) < 1e-8:
                right_faces.append(tag)
            elif abs(ymin) < 1e-8 and abs(ymax) < 1e-8:
                bottom_faces.append(tag)
            elif abs(ymin - 1.0) < 1e-8 and abs(ymax - 1.0) < 1e-8:
                moving_faces.append(tag)

        _physical(2, zmin_faces, "frontSource")
        _physical(2, zmax_faces, "backSource")
        _physical(2, left_faces, "leftWall")
        _physical(2, right_faces, "rightWall")
        _physical(2, bottom_faces, "bottomWall")
        _physical(2, moving_faces, "movingTop")
        _physical(3, volume_tags, "fluid")

        gmsh.write(str(msh_path))

        element_counts: dict[str, int] = {"triangle": 0, "quadrilateral": 0}
        for element_type, _, connectivity in zip(
            *gmsh.model.mesh.getElements(3)
        ):
            name = gmsh.model.mesh.getElementProperties(element_type)[0].lower()
            if "prism" in name:
                element_counts["triangle"] += len(connectivity) // 6
            elif "hexahedron" in name:
                element_counts["quadrilateral"] += len(connectivity) // 8
        metadata_path.write_text(
            json.dumps(
                {
                    "meshType": "hybrid",
                    "resolution": resolution,
                    "targetSize": target_size,
                    "boundaryLayerThickness": b,
                    "boundaryLayerCount": boundary_layer_count,
                    "coreTargetSize": target_size,
                    "thickness": thickness,
                    "elementCounts": element_counts,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    finally:
        gmsh.finalize()

    print(f"mesh={msh_path}")
    print(f"meshMetadata={metadata_path}")
    return msh_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--resolution", type=int, required=True)
    parser.add_argument("--boundary-layer-thickness", type=float, default=0.12)
    parser.add_argument("--boundary-layer-count", type=int, default=None)
    parser.add_argument("--thickness", type=float, default=0.1)
    args = parser.parse_args()
    generate(
        args.case.resolve(),
        args.resolution,
        args.boundary_layer_thickness,
        args.boundary_layer_count,
        args.thickness,
    )


if __name__ == "__main__":
    main()
