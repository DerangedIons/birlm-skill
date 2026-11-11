"""The session clock: one deadline for a whole run of the skill, preflight to generated code

It is a deadline, not a ledger. ``bimgrl/clock.json`` holds when the session started and how
long it may take; the time left is ``started_at + limit_s - now``, so the time spent at a gate
or between two commands is charged without anyone measuring it. Each command sizes its work to
the time left, and the ladder below says which optional steps are still allowed. A tail reserve
is held back for ``codegen`` and the hand-back. No clock file means no limit.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

LIMIT_S = 3 * 3600.0
TAIL_RESERVE_S = 15 * 60.0

# (time left in seconds below which the step is cut, key, what is cut), coarsest first
LADDER = (
    (75 * 60.0, "long", "the 100-beat confirmation (sweep --beats above 10)"),
    (60 * 60.0, "timing", "time and report"),
    (45 * 60.0, "tuning", "further tuning rounds (sweep --force)"),
    (30 * 60.0, "extras", "the comparison schemes (only BIRL-M or BIMGRL is swept) and the stability bisection"),
    (TAIL_RESERVE_S, "sweep", "every sweep (codegen and the hand-back only)"),
)

# the schemes still swept once "extras" is cut
MAIN_SCHEMES = ("BIRL-M", "BIMGRL")


def parse_duration(text: str) -> float:
    """``3h``, ``20m``, ``90s`` or a bare number of seconds"""
    text = text.strip().lower()
    unit = {"h": 3600.0, "m": 60.0, "s": 1.0}.get(text[-1:])
    return float(text[:-1]) * unit if unit else float(text)


def hm(seconds: float) -> str:
    minutes = int(max(seconds, 0.0) // 60)
    return f"{minutes // 60}h{minutes % 60:02d}"


class Clock:
    def __init__(self, path: Path, data: dict):
        self.path = path
        self.data = data

    @classmethod
    def load(cls, out: Path) -> Clock | None:
        """The clock of the outputs directory ``<model_dir>/bimgrl``, or None if none was started"""
        path = Path(out) / "clock.json"
        return cls(path, json.loads(path.read_text())) if path.exists() else None

    @classmethod
    def start(cls, out: Path, started_at: float | None = None, limit_s: float = LIMIT_S) -> Clock:
        clock = cls(Path(out) / "clock.json", dict(started_at=time.time() if started_at is None else started_at,
                                                   limit_s=limit_s, commands=[], skipped=[]))
        clock.save()
        return clock

    def save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(self.data, indent=2) + "\n")

    def used(self) -> float:
        return time.time() - self.data["started_at"]

    def remaining(self) -> float:
        return self.data["limit_s"] - self.used()

    def allowance(self) -> float:
        """Seconds a command may still spend: the time left less the tail reserve"""
        return max(self.remaining() - TAIL_RESERVE_S, 0.0)

    def cut(self, key: str) -> str | None:
        """Why the step ``key`` of the ladder is cut at the time left, or None if it is still allowed"""
        left = self.remaining()
        for mark, k, what in LADDER:
            if k == key and left < mark:
                return f"{hm(left)} left, below the {hm(mark)} mark for {what}"
        return None

    def cuts(self) -> list[str]:
        left = self.remaining()
        return [what for mark, _, what in LADDER if left < mark]

    def record(self, command: str, started: float, wall_s: float) -> None:
        self.data.setdefault("commands", []).append(dict(command=command, started=started, wall_s=round(wall_s, 3)))
        self.save()

    def skip(self, command: str, why: str) -> None:
        self.data.setdefault("skipped", []).append(dict(command=command, at=time.time(), why=why))
        self.save()

    def line(self) -> str:
        left = self.remaining()
        over = f" (over by {hm(-left)})" if left < 0 else ""
        cuts = self.cuts()
        return (f"clock: {hm(self.used())} used, {hm(left)} left of {hm(self.data['limit_s'])}{over} | "
                + ("cut: " + "; ".join(cuts) if cuts else "full"))
