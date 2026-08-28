#!/usr/bin/env python3

"""Generate the equilateral-triangle lid-driven-cavity hybrid prism mesh."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import gmsh


EPS = 1.0e-7


def _physical(dim: int, entities: list[int], name: str) -> None:
    if not entities:
        raise RuntimeError(f"No entities found for physical group {name}")
    tag = gmsh.model.addPhysicalGroup(dim, entities)
    gmsh.model.setPhysicalName(dim, tag, name)


def _distance_to_line(point: tuple[float, float], a: tuple[float, float], b: tuple[float, float]) -> float:
    px, py = point
    ax, ay = a
    bx, by = b
    return abs((bx - ax) * (ay - py) - (ax - px) * (by - ay)) / math.hypot(bx - ax, by - ay)


def _surface_center(surface: int) -> tuple[float, float, float]:
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, surface)
    return ((xmin + xmax) / 2.0, (ymin + ymax) / 2.0, (zmin + zmax) / 2.0)


def generate(
    case: Path,
    resolution: int,
    boundary_layer_thickness: float = 0.10,
    boundary_layer_count: int | None = None,
    thickness: float = 0.1,
) -> Path:
    """Generate ``case/mesh/mesh.msh`` for the triangular cavity.

    The two-dimensional domain is the problem-statement equilateral triangle:

        A = (-sqrt(3), 0), B = (sqrt(3), 0), C = (0, -3).

    A homothetic inner triangle creates three structured quadrilateral wall
    strips. The central triangle is left unrecombined, so the final prism mesh
    contains hexahedra in the wall strips and triangular prisms in the core.
    """
    if resolution < 12:
        raise ValueError("resolution must be at least 12 for the hybrid triangular cavity")
    if not 0.0 < boundary_layer_thickness < 0.35:
        raise ValueError("boundary_layer_thickness must be in (0, 0.35)")
    if thickness <= 0.0:
        raise ValueError("thickness must be positive")
    if boundary_layer_count is None:
        boundary_layer_count = max(4, round(resolution * boundary_layer_thickness))

    output_dir = case / "mesh"
    output_dir.mkdir(parents=True, exist_ok=True)
    msh_path = output_dir / "mesh.msh"
    metadata_path = output_dir / "mesh_metadata.json"

    sqrt3 = math.sqrt(3.0)
    a = (-sqrt3, 0.0)
    b = (sqrt3, 0.0)
    c = (0.0, -3.0)
    centroid = (0.0, -1.0)
    scale = 1.0 - boundary_layer_thickness
    ai = (centroid[0] + scale * (a[0] - centroid[0]), centroid[1] + scale * (a[1] - centroid[1]))
    bi = (centroid[0] + scale * (b[0] - centroid[0]), centroid[1] + scale * (b[1] - centroid[1]))
    ci = (centroid[0] + scale * (c[0] - centroid[0]), centroid[1] + scale * (c[1] - centroid[1]))

    side_count = resolution + 1
    layer_count = boundary_layer_count + 1
    core_size = 3.0 / resolution

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.option.setNumber("Mesh.RecombineAll", 0)
        gmsh.model.add(f"triangular_lid_cavity_hybrid_N{resolution}")
        geo = gmsh.model.geo

        pa = geo.addPoint(*a, 0.0)
        pb = geo.addPoint(*b, 0.0)
        pc = geo.addPoint(*c, 0.0)
        pai = geo.addPoint(*ai, 0.0)
        pbi = geo.addPoint(*bi, 0.0)
        pci = geo.addPoint(*ci, 0.0)

        ab = geo.addLine(pa, pb)
        bc = geo.addLine(pb, pc)
        ca = geo.addLine(pc, pa)
        abi = geo.addLine(pai, pbi)
        bci = geo.addLine(pbi, pci)
        cai = geo.addLine(pci, pai)
        a_ai = geo.addLine(pa, pai)
        b_bi = geo.addLine(pb, pbi)
        c_ci = geo.addLine(pc, pci)

        top_loop = geo.addCurveLoop([ab, b_bi, -abi, -a_ai])
        right_loop = geo.addCurveLoop([bc, c_ci, -bci, -b_bi])
        left_loop = geo.addCurveLoop([ca, a_ai, -cai, -c_ci])
        core_loop = geo.addCurveLoop([abi, bci, cai])
        top_strip = geo.addPlaneSurface([top_loop])
        right_strip = geo.addPlaneSurface([right_loop])
        left_strip = geo.addPlaneSurface([left_loop])
        core = geo.addPlaneSurface([core_loop])
        geo.synchronize()

        for line in (ab, bc, ca, abi, bci, cai):
            gmsh.model.mesh.setTransfiniteCurve(line, side_count)
        for line in (a_ai, b_bi, c_ci):
            gmsh.model.mesh.setTransfiniteCurve(line, layer_count)
        for surface in (top_strip, right_strip, left_strip):
            gmsh.model.mesh.setTransfiniteSurface(surface)
            gmsh.model.mesh.setRecombine(2, surface)
        gmsh.model.mesh.setSize(gmsh.model.getBoundary([(2, core)], recursive=True), core_size)

        extruded = geo.extrude(
            [(2, top_strip), (2, right_strip), (2, left_strip), (2, core)],
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

        for dim, tag in gmsh.model.getEntities(2):
            if tag in (top_strip, right_strip, left_strip):
                gmsh.model.mesh.setRecombine(dim, tag)

        gmsh.model.mesh.generate(3)

        front_faces: list[int] = []
        back_faces: list[int] = []
        moving_faces: list[int] = []
        left_faces: list[int] = []
        right_faces: list[int] = []
        for dim, tag in gmsh.model.getEntities(2):
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
            cx, cy, _ = _surface_center(tag)
            if abs(zmin) < EPS and abs(zmax) < EPS:
                front_faces.append(tag)
            elif abs(zmin - thickness) < EPS and abs(zmax - thickness) < EPS:
                back_faces.append(tag)
            elif abs(ymin) < EPS and abs(ymax) < EPS:
                moving_faces.append(tag)
            elif _distance_to_line((cx, cy), a, c) < 1.0e-5:
                left_faces.append(tag)
            elif _distance_to_line((cx, cy), b, c) < 1.0e-5:
                right_faces.append(tag)

        _physical(2, front_faces, "frontSource")
        _physical(2, back_faces, "backSource")
        _physical(2, moving_faces, "movingTop")
        _physical(2, left_faces, "leftWall")
        _physical(2, right_faces, "rightWall")
        _physical(3, volume_tags, "fluid")
        gmsh.write(str(msh_path))

        element_counts: dict[str, int] = {"triangularPrism": 0, "hexahedron": 0}
        for element_type, _, connectivity in zip(*gmsh.model.mesh.getElements(3)):
            name = gmsh.model.mesh.getElementProperties(element_type)[0].lower()
            if "prism" in name:
                element_counts["triangularPrism"] += len(connectivity) // 6
            elif "hexahedron" in name:
                element_counts["hexahedron"] += len(connectivity) // 8

        metadata_path.write_text(
            json.dumps(
                {
                    "meshType": "hybrid",
                    "problem": "tri_lid_driven_cavity",
                    "resolution": resolution,
                    "vertices": {"A": a, "B": b, "C": c},
                    "innerVertices": {"A": ai, "B": bi, "C": ci},
                    "boundaryLayerThickness": boundary_layer_thickness,
                    "boundaryLayerCount": boundary_layer_count,
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
    parser.add_argument("--case", required=True, type=Path)
    parser.add_argument("--resolution", required=True, type=int)
    parser.add_argument("--boundary-layer-thickness", default=0.10, type=float)
    parser.add_argument("--boundary-layer-count", default=None, type=int)
    parser.add_argument("--thickness", default=0.1, type=float)
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
