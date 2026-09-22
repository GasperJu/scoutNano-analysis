"""Lightweight ordered cutflow utility."""

from __future__ import annotations

import json
from collections import OrderedDict


class CutflowCounter:
    def __init__(self):
        self._counts = OrderedDict()

    def register_stages(self, stages):
        for stage in stages:
            self._counts.setdefault(stage, 0)

    def increment(self, stage, value=1):
        if stage not in self._counts:
            self._counts[stage] = 0
        self._counts[stage] += int(value)

    def rows(self):
        stages = list(self._counts.items())
        if not stages:
            return []
        total = stages[0][1]
        prev = total
        out = []
        for idx, (name, count) in enumerate(stages):
            if idx == 0:
                incr = 100.0 if count > 0 else 0.0
            else:
                incr = (100.0 * count / prev) if prev > 0 else 0.0
            cum = (100.0 * count / total) if total > 0 else 0.0
            out.append(
                {
                    "stage": name,
                    "count": int(count),
                    "incremental_eff": incr,
                    "cumulative_eff": cum,
                }
            )
            prev = count
        return out

    def format_report(self):
        rows = self.rows()
        if not rows:
            return ""
        lines = []
        lines.append("cut name | pass | eff | cumulative eff")
        lines.append("-" * 64)
        for row in rows:
            lines.append(
                f"{row['stage']:<20} | {row['count']:>8d} | "
                f"eff-{row['incremental_eff']:6.2f}% | cumulative eff-{row['cumulative_eff']:6.2f}%"
            )
        return "\n".join(lines)

    def to_dict(self):
        return {"stages": self.rows()}

    def write_json(self, path):
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(self.to_dict(), handle, indent=2)
