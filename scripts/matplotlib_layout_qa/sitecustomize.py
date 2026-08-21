"""Conservative draw-time warnings for matplotlib teaching figures.

Put this directory on ``PYTHONPATH`` while exporting notebooks. Python imports
``sitecustomize`` automatically, so the check observes figures without adding
QA code or a runtime dependency to the browser-hosted notebooks.
"""

from __future__ import annotations

import warnings

from matplotlib.figure import Figure


_ORIGINAL_DRAW = Figure.draw


def _visible_artists(figure: Figure):
    artists = []
    for axis_index, axis in enumerate(figure.axes, start=1):
        prefix = "colorbar" if axis.get_label() == "<colorbar>" else f"axis {axis_index}"
        artists.extend(
            [
                (f"{prefix} title", "title", axis.title),
                (f"{prefix} x label", "label", axis.xaxis.label),
                (f"{prefix} y label", "label", axis.yaxis.label),
            ]
        )
        artists.extend(
            (f"{prefix} x tick", "tick", artist)
            for artist in axis.get_xticklabels()
        )
        artists.extend(
            (f"{prefix} y tick", "tick", artist)
            for artist in axis.get_yticklabels()
        )
        artists.extend(
            (f"{prefix} annotation", "annotation", artist)
            for artist in axis.texts
        )
        legend = axis.get_legend()
        if legend is not None:
            artists.append((f"{prefix} legend", "legend", legend))
    artists.extend(
        ("figure annotation", "annotation", artist) for artist in figure.texts
    )
    return [
        (name, kind, artist)
        for name, kind, artist in artists
        if artist.get_visible()
        and (not hasattr(artist, "get_text") or artist.get_text().strip())
    ]


def _layout_warnings(figure: Figure, renderer) -> tuple[str, ...]:
    boxes = [
        (name, kind, artist.get_window_extent(renderer))
        for name, kind, artist in _visible_artists(figure)
    ]
    frame = figure.bbox
    issues = [
        f"{name} may be clipped"
        for name, _kind, box in boxes
        if box.x0 < frame.x0 - 12
        or box.y0 < frame.y0 - 12
        or box.x1 > frame.x1 + 12
        or box.y1 > frame.y1 + 12
    ]
    for index, (name_a, kind_a, box_a) in enumerate(boxes):
        for name_b, kind_b, box_b in boxes[index + 1 :]:
            if {kind_a, kind_b} <= {"tick", "annotation"}:
                continue
            width = max(0.0, min(box_a.x1, box_b.x1) - max(box_a.x0, box_b.x0))
            height = max(0.0, min(box_a.y1, box_b.y1) - max(box_a.y0, box_b.y0))
            smaller_area = min(box_a.width * box_a.height, box_b.width * box_b.height)
            if smaller_area > 0 and width * height / smaller_area > 0.25:
                issues.append(f"{name_a} intersects {name_b}")
    return tuple(dict.fromkeys(issues))


def _draw_with_layout_qa(self: Figure, renderer) -> None:
    _ORIGINAL_DRAW(self, renderer)
    self._ssi_layout_warnings = _layout_warnings(self, renderer)
    if self._ssi_layout_warnings:
        warnings.warn(
            "Possible matplotlib layout issue: "
            + "; ".join(self._ssi_layout_warnings),
            RuntimeWarning,
            stacklevel=2,
        )


if not getattr(Figure, "_ssi_layout_qa_installed", False):
    Figure.draw = _draw_with_layout_qa
    Figure._ssi_layout_qa_installed = True